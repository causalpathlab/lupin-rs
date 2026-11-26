//! `lupin ask`: annotation without a marker panel, through any AI chat.
//!
//! It writes a first round whose clusters carry no labels yet (with the
//! statistics cache, so later rounds are rescored like any enrichment run)
//! and a prompt listing each cluster's top specific genes. Pasted into an AI
//! chat, the prompt asks for decision lines only; pasting the answer into
//! `lupin relabel -f <round> -d - --next` applies it, with every decision's
//! rationale kept in the round's history. Nothing is sent anywhere by lupin.
//!
//! With `--blind`, on an annotated round, it writes no round: it writes the
//! round's blind prompt beside it ([`crate::manifest::reading::BLIND_PROMPT`])
//! and prints it. The prompt lists each cluster's most specific genes and
//! nothing of the panel, the labels or the methods, and asks for `reading`
//! lines: a reader's labels, recorded as evidence beside the methods'
//! (`lupin relabel`), never as decisions.

use crate::annotate::rounds::ClusterId;
use crate::manifest::annotate::load_enrichment_inputs;
use crate::manifest::recalibrate::write_cache;
use crate::manifest::rounds::{self, write_argmax, write_clusters, write_summary, CLUSTERS};
use crate::manifest::run::{self, annotated_path, rel_to_manifest, StatsCache};
use anyhow::{Context, Result};
use clap::Args;
use legume_numeric::matrix::common_io::mkdir_parent;
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::traits::IoOps;
use log::info;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;

#[derive(Args, Debug)]
pub struct AskArgs {
    #[arg(long, short = 'f', help = "Run manifest or its output prefix")]
    pub from: Box<str>,

    #[arg(
        long,
        short = 'o',
        required_unless_present = "blind",
        help = "Output prefix for the first round and the prompt"
    )]
    pub out: Option<Box<str>>,

    #[arg(
        long,
        conflicts_with_all = ["out", "clusters"],
        help = "On an annotated round: write no round, and print the round's blind prompt \
                (its clusters' most specific genes, nothing of its panel or labels), \
                answered with `reading` lines for `lupin relabel`"
    )]
    pub blind: bool,

    #[arg(
        long,
        help = "What the cells are (tissue, species, condition), passed to the AI as context"
    )]
    pub context: Option<String>,

    #[arg(
        long,
        default_value_t = DEFAULT_TOP,
        help = "Top specific genes listed per cluster"
    )]
    pub top: usize,

    #[arg(
        long,
        help = "Cluster parquet (cells × `cluster`); defaults to the run's clusters"
    )]
    pub clusters: Option<Box<str>>,
}

const PROMPT: &str = ".ask_prompt.md";

/// Specific genes listed per cluster unless asked otherwise.
pub const DEFAULT_TOP: usize = 15;

/// Write the first round and print the prompt.
pub fn run_ask(args: &AskArgs) -> Result<()> {
    let loaded = run::load(&args.from)?;
    if args.blind {
        let prompt = write_blind(&loaded, args.context.as_deref(), args.top)?;
        println!("{prompt}");
        return Ok(());
    }
    let out = args
        .out
        .as_deref()
        .context("--out names the first round to write")?
        .to_string();
    run::may_replace(&annotated_path(&loaded.file, &out))?;
    mkdir_parent(&out)?;
    let mut eargs = crate::annotate_cmd::default_enrichment_args(&out);
    eargs.clusters.clone_from(&args.clusters);
    // `ask` aggregates and asks; it runs no SuSiE stage, and records none.
    eargs.susie = None;
    let inputs = load_enrichment_inputs(&eargs, &loaded, None)?;

    // The round: clusters, no labels, the profile, and the cache.
    let clusters_path = format!("{out}{CLUSTERS}");
    let ids = rounds::cluster_ids(&inputs.cluster_labels);
    write_clusters(&clusters_path, &inputs.cell_names, &ids)?;
    let n = inputs.cell_names.len();
    let argmax_path = write_argmax(&out, &inputs.cell_names, &vec![None; n], &vec![f32::NAN; n])?;
    let names: Vec<Box<str>> = (0..inputs.n_clusters)
        .map(|k| format!("K{k}").into())
        .collect();
    let profile_path = format!("{out}.cluster_expression.parquet");
    inputs.profile_gk.to_parquet_with_names(
        &profile_path,
        (Some(&inputs.gene_names), Some("gene")),
        Some(&names),
    )?;
    let cache = write_cache(&out, &inputs)?;

    let mut next = loaded.copy_to(annotated_path(&loaded.file, &out))?;
    let rel = |p: &str| rel_to_manifest(&next.dir, p);
    next.manifest.cluster.clusters = Some(rel(&clusters_path));
    let a = &mut next.manifest.annotate;
    a.argmax = Some(rel(&argmax_path));
    a.cluster_expression = Some(rel(&profile_path));
    a.expression_clusters = Some(rel(&clusters_path));
    a.stats_cache = Some(StatsCache {
        gene_sum: rel(&cache.gene_sum),
        batch_profile: rel(&cache.batch_profile),
        gene_weight: rel(&cache.gene_weight),
        cell_batch: rel(&cache.cell_batch),
    });
    a.settings = Some(serde_json::json!({ "enrichment": eargs }));
    a.drop_susie_tables();
    for p in [
        &mut a.markers,
        &mut a.log,
        &mut a.history,
        &mut a.celltype_tree,
        &mut a.fine_argmax,
        &mut a.cluster_celltype_q,
        &mut a.cluster_celltype_q_values,
        &mut a.cluster_celltype_p,
        &mut a.cluster_celltype_nes,
        &mut a.marker_score,
        &mut a.units_call,
        &mut a.agent_readings,
    ] {
        *p = None;
    }
    a.called_by = None;
    a.stats = None;
    write_summary(&mut next, &out)?;
    next.manifest.save(&next.file)?;

    let sizes = cluster_sizes(&inputs.cluster_labels, inputs.n_clusters);
    let top = top_specific_genes(&inputs.profile_gk, &inputs.gene_names, args.top);
    let prompt = prompt(
        &next.file.to_string_lossy(),
        args.context.as_deref(),
        &sizes,
        &top,
    );
    let prompt_path = format!("{out}{PROMPT}");
    fs::write(&prompt_path, &prompt)?;
    info!(
        "wrote {prompt_path}; paste the AI's answer into `lupin relabel -f {} -d - --next`",
        next.file.display()
    );
    println!("{prompt}");
    Ok(())
}

/// Write the round `loaded`'s blind prompt beside it and return it, with
/// its id: each cluster's `top` most specific genes from the round's
/// profile (over its clusters as they now are), and nothing of its panel,
/// labels or methods' readings.
pub fn write_blind(loaded: &run::Loaded, context: Option<&str>, top: usize) -> Result<String> {
    use super::reading::{blind_prompt_path, prompt_id, PROMPT_ID};
    let cells = rounds::read_cells(loaded)?;
    let profile = rounds::round_profile(loaded, &cells)?.with_context(|| {
        format!(
            "{} records no cluster expression profile to read genes from",
            loaded.file.display()
        )
    })?;
    let ids: Vec<ClusterId> = (profile.cols.iter())
        .map(|c| {
            crate::annotate::rounds::parse_cluster_id(c)
                .context("a profile column is not a cluster")
        })
        .collect::<Result<_>>()?;
    let mut sizes: BTreeMap<ClusterId, usize> = BTreeMap::new();
    for id in cells.clusters.iter().flatten() {
        *sizes.entry(*id).or_default() += 1;
    }
    let genes = top_specific_genes(&profile.mat, &profile.rows, top);
    // Absolute, so a reading names it from wherever its answer is saved.
    let file = loaded
        .file
        .canonicalize()
        .unwrap_or_else(|_| loaded.file.clone());
    let round = file.to_string_lossy();
    let body = blind_prompt(&round, context, &ids, &sizes, &genes);
    let prompt = format!("{body}{PROMPT_ID}{}\n", prompt_id(&body));
    let path = blind_prompt_path(loaded);
    fs::write(&path, &prompt).with_context(|| format!("writing {path}"))?;
    info!(
        "wrote {path}; give it to a reader that has seen nothing else of this analysis, then \
         `lupin relabel -f {round} -d - --next` with its answer"
    );
    Ok(prompt)
}

/// The blind prompt's text: context, the answer format, then the clusters.
/// Its id goes after it ([`super::reading::PROMPT_ID`]).
fn blind_prompt(
    round: &str,
    context: Option<&str>,
    ids: &[ClusterId],
    sizes: &BTreeMap<ClusterId, usize>,
    top: &[Vec<(String, f32)>],
) -> String {
    let mut p = String::new();
    let _ = writeln!(
        p,
        "You are reading cell clusters from single-cell RNA-seq.{}",
        context.map_or(String::new(), |c| format!(" The cells are: {c}."))
    );
    let _ = writeln!(
        p,
        "For each cluster below, say what cell type it is from its most specific genes alone \
         (in brackets: log fold change of its expression over the other clusters'), with the \
         Cell Ontology term for it. Use only these genes and what you know; do not look at any \
         other output of this analysis. Very small clusters may be noise or doublets: say so \
         rather than guess.\n\n\
         Answer ONLY with JSON lines, one per cluster, no other text:\n\
         {{\"cluster\": 0, \"action\": \"reading\", \"label\": \"<cell type>\", \"cl\": \"CL:<7 digits>\", \"features\": [\"<gene>\", \"<gene>\"], \"rationale\": \"<one sentence naming the genes>\", \"prompt\": \"<the prompt id at the end>\", \"model\": \"<your model>\", \"decided_by\": \"agent\", \"round\": \"{round}\"}}\n\
         If the genes name no cell type, give your best broad label and say why in the rationale.\n\n\
         Clusters:"
    );
    cluster_list(
        &mut p,
        (ids.iter().zip(top))
            .map(|(id, g)| (id, sizes.get(id).copied().unwrap_or(0), g.as_slice())),
    );
    p
}

fn cluster_sizes(labels: &[usize], k: usize) -> Vec<usize> {
    let mut sizes = vec![0; k];
    for &l in labels {
        if l < k {
            sizes[l] += 1;
        }
    }
    sizes
}

/// Per cluster, the genes most specific to it: highest log fold change of
/// its `log1p` counts-per-10k over the mean of the other clusters, among the genes
/// it expresses strongly (at or above the 90th percentile of its expressed
/// genes), so a faint transcript cannot top the list.
fn top_specific_genes(profile: &Mat, genes: &[Box<str>], top: usize) -> Vec<Vec<(String, f32)>> {
    let (g, k) = (profile.nrows(), profile.ncols());
    // The profile holds each cluster's share of its counts; per 10k, as usual.
    let logp = |i: usize, c: usize| (1e4 * profile[(i, c)].max(0.0)).ln_1p();
    let totals: Vec<f32> = (0..g).map(|i| (0..k).map(|c| logp(i, c)).sum()).collect();
    (0..k)
        .map(|c| {
            let mut expressed: Vec<f32> = (0..g).map(|i| logp(i, c)).filter(|v| *v > 0.0).collect();
            expressed.sort_by(f32::total_cmp);
            let floor = expressed
                .get(expressed.len() * 9 / 10)
                .copied()
                .unwrap_or(f32::INFINITY);
            let others = (k.max(2) - 1) as f32;
            let mut fc: Vec<(usize, f32)> = (0..g)
                .filter(|&i| logp(i, c) >= floor)
                .map(|i| (i, logp(i, c) - (totals[i] - logp(i, c)) / others))
                .collect();
            fc.sort_by(|a, b| b.1.total_cmp(&a.1));
            fc.into_iter()
                .take(top)
                .map(|(i, f)| (genes[i].to_string(), f))
                .collect()
        })
        .collect()
}

/// The prompt: context, the answer format, then the clusters.
fn prompt(
    round: &str,
    context: Option<&str>,
    sizes: &[usize],
    top: &[Vec<(String, f32)>],
) -> String {
    let mut p = String::new();
    let _ = writeln!(
        p,
        "You are helping annotate cell clusters from single-cell RNA-seq.{}",
        context.map_or(String::new(), |c| format!(" The cells are: {c}."))
    );
    let _ = writeln!(
        p,
        "For each cluster below, propose a cell-type label from its most specific genes \
         (in brackets: log fold change of its expression over the other clusters'). \
         Very small clusters may be noise or doublets; say so rather than guess.\n\n\
         Answer ONLY with JSON lines, one decision per line, no other text:\n\
         {{\"cluster\": 0, \"action\": \"label\", \"label\": \"<cell type>\", \"rationale\": \"<one sentence naming the genes>\", \"evidence\": [{{\"kind\": \"marker\", \"term\": \"<gene>\"}}], \"decided_by\": \"agent_proposed_user_accepted\"}}\n\
         If the genes do not support a label, answer for that cluster:\n\
         {{\"cluster\": 0, \"action\": \"keep\", \"rationale\": \"<why>\", \"decided_by\": \"agent_proposed_user_accepted\"}}\n\
         For each label you propose, also list its marker genes from the lists below:\n\
         {{\"action\": \"markers_add\", \"label\": \"<cell type>\", \"features\": [\"<gene>\", \"<gene>\"], \"rationale\": \"<why>\", \"decided_by\": \"agent_proposed_user_accepted\"}}\n\n\
         (Round: {round})\n\nClusters:"
    );
    cluster_list(
        &mut p,
        top.iter()
            .enumerate()
            .map(|(c, g)| (c, sizes[c], g.as_slice())),
    );
    p
}

/// A prompt's clusters, one line each: id, cells, and the most specific
/// genes with their log fold change.
fn cluster_list<'a, I: std::fmt::Display>(
    p: &mut String,
    clusters: impl IntoIterator<Item = (I, usize, &'a [(String, f32)])>,
) {
    for (id, n, genes) in clusters {
        let list: Vec<String> = genes.iter().map(|(g, z)| format!("{g} ({z:.1})")).collect();
        let _ = writeln!(p, "- cluster {id} ({n} cells): {}", list.join(", "));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn top_genes_are_the_clusters_own_specific_ones() {
        let mut m = Mat::zeros(4, 2);
        // GENE0 high in cluster 0, GENE1 in cluster 1, GENE2 everywhere.
        for (g, c, v) in [
            (0, 0, 50.0),
            (1, 1, 50.0),
            (2, 0, 20.0),
            (2, 1, 20.0),
            (3, 0, 1.0),
            (3, 1, 1.0),
        ] {
            m[(g, c)] = v;
        }
        let genes: Vec<Box<str>> = (0..4).map(|i| format!("GENE{i}").into()).collect();
        let top = top_specific_genes(&m, &genes, 1);
        assert_eq!(top[0][0].0, "GENE0");
        assert_eq!(top[1][0].0, "GENE1");
    }

    #[test]
    fn the_prompt_asks_for_decision_lines_and_lists_every_cluster() {
        let top = vec![
            vec![("GENE1".to_string(), 3.0)],
            vec![("GENE2".to_string(), 2.5)],
        ];
        let p = prompt(
            "r/run.senna.json",
            Some("placeholder tissue"),
            &[10, 20],
            &top,
        );
        assert!(p.contains("Answer ONLY with JSON lines"));
        assert!(p.contains("placeholder tissue"));
        assert!(p.contains("- cluster 1 (20 cells): GENE2 (2.5)"));
    }
}
