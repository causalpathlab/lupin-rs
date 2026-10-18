//! A trajectory's order in words, on made-up facts and a stub model.

use super::*;
use crate::describe::consensus::tests::Stub;
use crate::trajectory::edges::Verdict;

fn edge(a: &str, b: &str, verdict: Option<Verdict>) -> EdgeRow {
    EdgeRow {
        a: a.into(),
        b: b.into(),
        connectivity: 0.8,
        in_prior: matches!(verdict, Some(Verdict::Supported | Verdict::Unsupported)),
        verdict,
        order_agreement: 0.9,
    }
}

/// HSC_like → MPP_like, which parts toward ERY and toward MONO; DC is a
/// component of its own; NK_cells pairs with MONO as a candidate only.
fn facts() -> OrderFacts {
    OrderFacts {
        lineages: parse_lineages(
            "# lineage\tcomponent\tpath\n\
             L0\t0\tHSC_like > MPP_like > ERY\n\
             L1\t0\tHSC_like > MPP_like > MONO\n\
             L2\t1\tDC\n",
        ),
        edges: vec![
            edge("HSC_like", "MPP_like", Some(Verdict::Supported)),
            edge("MPP_like", "ERY", Some(Verdict::Supported)),
            edge("MPP_like", "MONO", Some(Verdict::Unsupported)),
            edge("MONO", "NK_cells", Some(Verdict::Candidate)),
            edge("ERY", "DC", None),
        ],
        medians: [
            ("HSC_like", 0.01),
            ("MPP_like", 0.2),
            ("ERY", 0.7),
            ("MONO", 0.6),
        ]
        .iter()
        .map(|(t, m)| ((*t).to_string(), *m))
        .collect(),
    }
}

fn types() -> BTreeSet<String> {
    facts().types().into_iter().collect()
}

#[test]
fn the_facts_give_roots_branch_points_and_ends() {
    let f = facts();
    assert_eq!(f.roots(), ["HSC_like", "DC"]);
    assert_eq!(f.branch_points(), ["MPP_like"], "where L0 and L1 part");
    assert_eq!(f.ends(), ["ERY", "MONO", "DC"]);
    assert_eq!(f.components(), 2);
    assert_eq!(
        f.types(),
        ["HSC_like", "MPP_like", "ERY", "MONO", "DC", "NK_cells"],
        "the lineages' types, then those only an edge names"
    );
}

#[test]
fn the_templates_say_the_order_and_each_types_place() {
    let f = facts();
    let run = run_template(&f);
    assert!(
        run.contains("The order starts at HSC_like and DC."),
        "{run}"
    );
    assert!(run.contains("It branches after MPP_like."), "{run}");
    assert!(run.contains("end at ERY, MONO and DC."), "{run}");
    assert!(run.contains("L0: HSC_like > MPP_like > ERY"), "{run}");
    let mpp = type_template(&f, "MPP_like");
    assert!(
        mpp.starts_with("MPP_like is on L0 and L1, after HSC_like and before ERY and MONO."),
        "{mpp}"
    );
    assert!(mpp.contains("MPP_like → MONO unsupported"), "{mpp}");
    let nk = type_template(&f, "NK_cells");
    assert!(nk.starts_with("NK_cells is on no lineage."), "{nk}");
    assert!(nk.contains("MONO → NK_cells candidate"), "{nk}");
}

/// The grounding test for the order: a type's line names only the type,
/// its edges' other ends, the types either side of it and its lineages'
/// roots.
#[test]
fn a_types_line_may_name_only_its_neighbours_and_roots() {
    let f = facts();
    let s = subjects(&f);
    let ery = s
        .iter()
        .find(|s| s.subject.as_deref() == Some("ERY"))
        .unwrap();
    assert_eq!(
        ery.own,
        BTreeSet::from(["MPP_like".to_string(), "HSC_like".to_string()])
    );
    let check = |t: &str| check_lm(t, &types(), &ery.named, &ery.own);
    assert_eq!(
        check("ERY follows MPP like on the way from HSC like."),
        Ok(())
    );
    assert_eq!(
        check("ERY sits beside MONO."),
        Err(super::super::CitationError::Foreign("MONO".into()))
    );
    assert!(check("ERY is the third step.").is_ok());
    assert!(check("ERY is step 3.").is_err(), "no numbers");
    // The run's paragraph names roots, branch points and ends only.
    let run = &s[0];
    assert!(run.subject.is_none());
    assert!(run.named.contains("MPP_like") && run.named.contains("DC"));
    assert!(
        check_lm(
            "It starts at HSC like and ends at NK cells.",
            &types(),
            &run.named,
            &run.own
        )
        .is_err(),
        "NK_cells is no root, branch point or end"
    );
    assert!(
        run.prompt.contains("BRANCH POINTS: MPP_like"),
        "{}",
        run.prompt
    );
    assert!(
        run.prompt.ends_with("<think>\n\n</think>\n\n"),
        "{}",
        run.prompt
    );
}

#[test]
fn a_model_text_that_strays_falls_back_to_the_template() {
    let f = facts();
    let s = &subjects(&f)[0];
    let mut stub = Stub::new(
        &[],
        &["It starts at NK cells.", "It starts at NK cells again."],
    );
    let d = word(&mut stub, "m", "r", s, &types()).unwrap();
    assert_eq!(d.check, Check::Template);
    assert_eq!(d.text, run_template(&f));
    assert_eq!(stub.asked[1].1.words, ["NK_cells"]);
    let mut stub = Stub::new(&[], &["The order runs from HSC like through MPP like."]);
    let d = word(&mut stub, "m", "r", s, &types()).unwrap();
    assert_eq!(d.check, Check::Ok);
    assert_eq!(d.cited, ["HSC_like", "MPP_like"]);
}

#[test]
fn the_texts_are_cached_and_current_falls_back_to_the_templates() {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("t.consensus.jsonl");
    let f = facts();
    let n = subjects(&f).len();
    // Before any model: the templates.
    let now = current(&f, &cache);
    assert_eq!(now.len(), n);
    assert!(now.iter().all(|d| d.check == Check::Template));
    let spec = LmSpec {
        model: super::super::lm::LmModel::Qwen3_0_6b,
        local_dir: None,
        device: crate::gene_text::cli::ComputeDevice::Cpu,
    };
    let lines: Vec<String> = f
        .types()
        .iter()
        .map(|t| format!("{t} sits here."))
        .collect();
    let mut all = vec!["The order runs from HSC like.".to_string()];
    all.extend(lines);
    let load = || -> Result<Box<dyn Lm>> {
        let refs: Vec<&str> = all.iter().map(String::as_str).collect();
        Ok(Box::new(Stub::new(&[], &refs)))
    };
    let made = describe(&f, Some(&spec), &cache, load).unwrap();
    assert!(made.iter().all(|d| d.check == Check::Ok), "{made:?}");
    assert_eq!(current(&f, &cache), made);
    let again = describe(&f, Some(&spec), &cache, || -> Result<Box<dyn Lm>> {
        panic!("all cached")
    })
    .unwrap();
    assert_eq!(again, made);
    let md = markdown("t.senna.json", &made);
    assert!(md.starts_with("# Trajectory description\n\nFrom `t.senna.json`."));
    assert!(md.contains("## MPP_like\n\nMPP_like sits here.\n"), "{md}");
    // No model asked for: the templates, nothing cached.
    let plain = describe(
        &f,
        None,
        &dir.path().join("none.jsonl"),
        || -> Result<Box<dyn Lm>> { panic!("no model") },
    )
    .unwrap();
    assert_eq!(plain[0].text, run_template(&f));
    assert!(!dir.path().join("none.jsonl").exists());
}

/// A trajectory manifest in `dir` with the facts' lineages and edges and a
/// pseudotime per cell.
pub(crate) fn fixture_trajectory(dir: &Path) -> PathBuf {
    use crate::manifest::run::RunManifest;
    use legume_numeric::matrix::parquet::{write_named_table, Column};
    let f = facts();
    let at = |name: &str| dir.join(name).to_string_lossy().into_owned();
    std::fs::write(
        dir.join("t.trajectory_lineages.tsv"),
        "# lineage\tcomponent\tpath\nL0\t0\tHSC_like > MPP_like > ERY\nL1\t0\tHSC_like > MPP_like > MONO\nL2\t1\tDC\n",
    )
    .unwrap();
    edges::write(&at("t.trajectory_edges.parquet"), &f.edges).unwrap();
    let cells: Vec<Box<str>> = (0..6).map(|i| format!("c{i}").into()).collect();
    let types: Vec<Box<str>> = ["HSC_like", "HSC_like", "MPP_like", "ERY", "ERY", "MONO"]
        .iter()
        .map(|s| Box::from(*s))
        .collect();
    let pt = [0.0f32, 0.02, 0.2, 0.6, 0.8, 0.6];
    write_named_table(
        &at("t.cell_pseudotime.parquet"),
        "cell",
        &cells,
        &[
            ("pseudotime".into(), Column::F32(&pt)),
            ("type".into(), Column::Str(&types)),
        ],
    )
    .unwrap();
    let mut m = RunManifest::new(crate::manifest::run::RunKind::Topic, "t");
    m.trajectory.lineages = Some("t.trajectory_lineages.tsv".into());
    m.trajectory.edges = Some("t.trajectory_edges.parquet".into());
    m.trajectory.pseudotime = Some("t.cell_pseudotime.parquet".into());
    let file = dir.join("t.senna.json");
    m.save(&file).unwrap();
    file
}

#[test]
fn the_facts_read_from_a_trajectory_run_with_each_types_median() {
    let dir = tempfile::tempdir().unwrap();
    let file = fixture_trajectory(dir.path());
    let loaded = crate::manifest::run::load(&file.to_string_lossy()).unwrap();
    let f = OrderFacts::read(&loaded).unwrap();
    assert_eq!(f.lineages, facts().lineages);
    assert_eq!(f.edges.len(), 5);
    assert!((f.medians["ERY"] - 0.7).abs() < 1e-6);
    assert!((f.medians["HSC_like"] - 0.01).abs() < 1e-6);
}
