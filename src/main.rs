//! `lupin`: Label, Unfold, Place, Interpret, Narrate.

mod annotate;
mod annotate_cmd;
mod cell_labels;
mod child;
mod describe;
mod docs;
mod figures;
mod gene_text;
mod hub;
mod manifest;
mod marker_embedding;
mod mcp;
mod progress;
#[cfg(test)]
#[path = "tests/session_cli.rs"]
mod session_cli_tests;
mod trajectory;
mod tui;

use crate::annotate_cmd::{run_annotate, AnnotateCliArgs};
use crate::gene_text::cli::{run_knn_graph, run_qc, KnnGraphCmd, QcCmd};
use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use describe::{run_describe, DescribeArgs};
use legume_tui::cli::SessionArgs;
use trajectory::run::{run_trajectory, TrajectoryArgs};

/// `bytes` as lower-case hex.
pub(crate) fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

#[derive(Parser)]
#[command(
    name = "lupin",
    version,
    about = "Label, Unfold, Place, Interpret, Narrate: text graphs, cell-type annotation, review, trajectories and description."
)]
struct Cli {
    #[arg(short, long, global = true, help = "Verbose logging")]
    verbose: bool,
    #[arg(
        long,
        global = true,
        help = "Replace existing manifests (.json) without asking; otherwise lupin asks, or refuses when it cannot ask"
    )]
    overwrite: bool,
    #[command(subcommand)]
    cmd: Commands,
}

#[derive(Subcommand)]
enum Commands {
    #[command(
        name = "text-qc",
        about = "Inspect the word vocabulary and its frequency cuts, without encoding",
        long_about = "The vocabulary step of `word-graph` on its own,\n\
                      so the cuts can be inspected before paying for the model pass.\n\
                      \n\
                      Tokenise every description, drop stopwords and filler,\n\
                      then cut both tails of the document-frequency distribution by quantile.\n\
                      Prints the df histogram, the cuts and the words on each side of them,\n\
                      and writes {out}.vocab.tsv.\n\
                      Tune the stopword list and the quantiles here,\n\
                      then hand the file to `word-graph --vocab-file`.\n\
                      `word-graph` runs this step itself when no file is given."
    )]
    TextQc(QcCmd),
    #[command(
        name = "word-graph",
        alias = "vocab-graph",
        about = "Encode the descriptions and write the text graph: feature–word and feature–feature edges",
        long_about = "Runs the vocabulary step,\n\
                      then a BERT-family encoder from the Hugging Face Hub over every description,\n\
                      and writes the text graph.\n\
                      \n\
                      {out}.feature_word.edges.tsv maps each feature to the words of its text\n\
                      (weight = contextual cosine × TF-IDF).\n\
                      {out}.knn_graph.edges.tsv lists the nearest features by text similarity.\n\
                      Both are typed edge files for `senna fne --edges`.\n\
                      Also writes {out}.text_embedding.parquet (pooled, centred),\n\
                      {out}.vocab.tsv and {out}.feature_text.tsv."
    )]
    WordGraph(KnnGraphCmd),
    #[command(
        name = "annotate",
        about = "Cell-type annotation by enrichment, embedding projection, or auto-dispatch",
        long_about = format!("Unified annotation entry point.\n\
                      \n\
                      `--method enrichment` scores the marker panel on per-cluster expression.\n\
                      `--method projection` calls cells on the run's embeddings: `--projection score`\n\
                      (the embedding marker score on the decoder gene table, per cell with its margin)\n\
                      where the run decodes cells, else `--projection units` (cell types, and GO terms\n\
                      with --go/--gaf/--gmt, trained as gene sets on the run's fixed genes).\n\
                      `--feature-embedding` + `--cell-embedding` run projection on those files directly\n\
                      (the units unless `--projection score`).\n\
                      `--method auto` (default) runs every method the run supports, their readings in\n\
                      one round: the enrichment (with SuSiE on raw counts), the projection score where\n\
                      the run decodes cells, and the units where its genes are placed among the cells\n\
                      (aligning the GO terms the enrichment tests). `--call` picks whose reading makes\n\
                      the round's labels (default: SuSiE's, else the enrichment's); a method that does\n\
                      not run is refused. `--projection` under auto runs that projection alone; a run\n\
                      with nothing for the enrichment to read (no counts, cache or decoder) is projected.\n\
                      Runs without a panel are the enrichment's.\n\
                      Ontology-only follow-up: `--from` + `--obo` + `--label-cl` without markers.\n\
                      \n\
                      Without `-o` it opens the TUI, which asks for what is missing.\n\
                      With `-f` on a round and `--pdf`, it saves the round's cell map, as the\n\
                      TUI opens on it, to a PDF instead (`--width` inches, the points at\n\
                      `--dpi`, the labels as text), without a terminal.\n\
                      \n\
                      {}",
            session_help("a pass or a save goes on", "`lupin annotate` with no --from")
        )
    )]
    Annotate {
        #[command(flatten)]
        args: AnnotateCliArgs,
        #[command(flatten)]
        session: SessionArgs,
        #[command(flatten)]
        pdf: legume_tui::export::PdfArgs,
    },
    #[command(
        name = "describe",
        about = "Short citation-checked sentence from annotate evidence",
        long_about = "Builds structured evidence from `{from}.annot.parquet` (or argmax).\n\
                      Optionally fishes per-cluster keywords from a `word-graph` prefix\n\
                      (`feature_word.edges` over each cluster's markers;\n\
                      `--feature-embedding` adds nearest-neighbour genes first).\n\
                      When `{from}.cluster_term_q.parquet` is present,\n\
                      a second FDR-significant contender is named if one exists.\n\
                      Writes `{out}.describe.json` and `{out}.describe.md`. On a trajectory\n\
                      run, code also writes its order in words from the run's outputs,\n\
                      `{out}.trajectory.describe.md`; its clusters are described as above\n\
                      when it carries a round's evidence, and only its order otherwise.\n\
                      \n\
                      The composer never invents labels: sentences are citation-checked\n\
                      templates; no language model writes them."
    )]
    Describe(DescribeArgs),
    #[command(
        about = "Show an annotation round cluster by cluster: calls, terms, CL placement, history",
        long_about = "Reads an annotated manifest (a round) and prints, per cluster,\n\
                      its size and current label, the candidate labels with q and support,\n\
                      the top GO/GMT terms, the Cell Ontology placement,\n\
                      and every decision made on it so far with its rationale.\n\
                      `--json` prints the same for a program or an agent to read."
    )]
    Review(manifest::rounds::ReviewArgs),
    #[command(
        about = "Apply a decisions file to an annotation round and write the next round",
        long_about = "Reads a round and a decisions file (JSONL; see `lupin review --help`),\n\
                      relabels or merges clusters, and writes a new manifest at `-o`\n\
                      whose `annotate.source` points back to the round it started from.\n\
                      Every decision needs a rationale, kept in `{out}.annotation_history.json`.\n\
                      Counts are not re-read; merged clusters take fresh ids."
    )]
    Relabel(manifest::rounds::RelabelArgs),
    #[command(about = "Print the method write-ups (omit the topic to list them)")]
    Docs(docs::DocsArgs),
    #[command(
        name = "trajectory",
        about = "Order cells along a supervised prior over their cell types: Cell Ontology and precedence statements, checked against the kNN graph, then diffusion pseudotime",
        long_about = format!(
            "Order cells along a supervised prior over their cell types: Cell Ontology\n\
             and precedence statements, checked against the kNN graph, then diffusion\n\
             pseudotime.\n\
             \n\
             Without both `--from` and `--out` it opens the TUI on the order view,\n\
             which asks for what is missing.\n\
             \n\
             {}",
            session_help("a pass or a trajectory run goes on", "`lupin trajectory` with no --from")
        )
    )]
    Trajectory {
        #[command(flatten)]
        args: TrajectoryArgs,
        #[command(flatten)]
        session: SessionArgs,
    },
    #[command(
        about = "Where lupin's data files (Cell Ontology, matching rules, aliases) come from; fetch them for offline use"
    )]
    Data(manifest::data_files::DataArgs),
    #[command(
        about = "Without a marker panel: write an unlabelled first round and a prompt for any AI chat",
        long_about = "Writes a first round whose clusters carry no labels yet,\n\
                      and prints a prompt listing each cluster's most specific genes.\n\
                      Paste it into an AI chat; the prompt asks for decision lines only.\n\
                      Paste the answer into `lupin relabel -f <round> -d - --next`:\n\
                      labels, rationales and suggested markers become the next round,\n\
                      rescored like any enrichment run. Nothing is sent anywhere by lupin."
    )]
    Ask(manifest::ask::AskArgs),
    #[command(
        about = "Return to a detached `lupin annotate` or `lupin trajectory`",
        long_about = "Return to a detached `lupin annotate` or `lupin trajectory` in this\n\
                      terminal.\n\
                      \n\
                      A session keeps running when the terminal goes (a dropped connection,\n\
                      or the detach key: ctrl+alt+d unless LEGUME_TUI_DETACH_KEY sets\n\
                      another; `?` shows it) and is picked up again here, or by starting it\n\
                      again with no --from. The first detach by the key asks for a label\n\
                      (Enter keeps the id); `lupin rename` or `r` in the picker changes it\n\
                      later. A session is found by its label or its id; with one live the\n\
                      name is optional, with several a picker opens. `lupin sessions` lists\n\
                      them. Inside one, ctrl+alt+s (or LEGUME_TUI_SWITCH_KEY) switches to\n\
                      another. --no-detach on either starts it plain."
    )]
    Attach(legume_tui::cli::AttachArgs),
    #[command(
        about = "List the live `lupin annotate` and `lupin trajectory` sessions",
        long_about = "List the live `lupin annotate` and `lupin trajectory` sessions, one a\n\
                      line: label and id, attached or detached, age, pid, command line and\n\
                      the run on screen, with the job running and the unsaved edits."
    )]
    Sessions,
    #[command(
        about = "Stop a `lupin annotate` or `lupin trajectory` session",
        long_about = "Stop a session as q does, without asking: a save under way finishes\n\
                      first, a running pass, trajectory, rescoring or term ranking is\n\
                      stopped, and unsaved edits are dropped. Waits up to a minute for it\n\
                      to end. --force kills it at once, and a pass it started may then\n\
                      keep running."
    )]
    Kill(legume_tui::cli::KillArgs),
    #[command(
        about = "Give a `lupin annotate` or `lupin trajectory` session a label",
        long_about = "Give a session a label, found by its label or its id. The label is\n\
                      unique among the live sessions; `attach`, `kill` and `rename` take it\n\
                      in place of the id, and `sessions` shows it."
    )]
    Rename(legume_tui::cli::RenameArgs),
    #[command(
        about = "Serve annotate, review, relabel and trajectory to an agent over MCP (stdio)",
        long_about = "A Model Context Protocol server on stdin/stdout for an agent\n\
                      (`claude mcp add lupin -- lupin mcp`). Its tools run annotation passes,\n\
                      review rounds, preview and apply decisions, and run and read trajectories.\n\
                      lupin states the facts; the agent's decisions carry a rationale, kept\n\
                      with the round. No tool replaces an existing manifest."
    )]
    Mcp,
}

/// The session paragraph of `lupin annotate`'s and `lupin trajectory`'s
/// long help: what goes on when the connection drops, who offers the live
/// sessions to return to, and where `?` shows the keys.
fn session_help(what_goes_on: &str, resume_hint: &str) -> String {
    legume_tui::cli::session_help(
        &tui::TOOL,
        legume_tui::cli::SessionHelp::new()
            .what_goes_on(what_goes_on)
            .resume_hint(Some(resume_hint))
            .more("`?` in the TUI shows the keys in effect."),
    )
}

impl Commands {
    /// Whether the command opens the TUI.
    fn opens_tui(&self) -> bool {
        match self {
            Self::Annotate { args, pdf, .. } => args.opens_tui() && pdf.pdf.is_none(),
            Self::Trajectory { args, .. } => args.opens_tui(),
            _ => false,
        }
    }
}

/// `lupin annotate`: with `--pdf`, the round's cell map saved without the
/// TUI; else a pass or the TUI.
fn annotate_command(
    args: &AnnotateCliArgs,
    session: &SessionArgs,
    pdf: &legume_tui::export::PdfArgs,
) -> Result<()> {
    if pdf.pdf.is_none() {
        return run_annotate(args, session);
    }
    let from = args
        .from
        .as_deref()
        .context("--pdf draws a round's cell map: give the round with -f")?;
    anyhow::ensure!(
        args.out.is_empty(),
        "--pdf draws a round's cell map and runs no pass: drop -o"
    );
    let wrote = tui::round_pdf(from, pdf)?;
    log::info!("wrote {}", wrote.display());
    Ok(())
}

fn main() -> Result<()> {
    // A session's host runs this binary again: it serves the session and
    // exits here, before any of lupin's own work.
    legume_tui::serve_if_host();
    let cli = Cli::parse();
    if cli.overwrite {
        manifest::run::allow_overwrite();
    }
    let interactive = cli.cmd.opens_tui();
    // On stderr, as RUST_LOG says; the TUI's log pane takes it too.
    let mut stderr = env_logger::Builder::from_env(
        env_logger::Env::default().default_filter_or(if cli.verbose { "debug" } else { "info" }),
    );
    if interactive {
        tui::init_logger(cli.verbose, stderr.build());
    } else {
        stderr.init();
    }
    match cli.cmd {
        Commands::TextQc(c) => run_qc(&c),
        Commands::WordGraph(c) => run_knn_graph(&c),
        Commands::Annotate { args, session, pdf } => annotate_command(&args, &session, &pdf),
        Commands::Describe(c) => run_describe(&c),
        Commands::Review(c) => manifest::rounds::run_review(&c),
        Commands::Relabel(c) => manifest::rounds::run_relabel(&c),
        Commands::Ask(c) => manifest::ask::run_ask(&c),
        Commands::Docs(c) => docs::run_docs(&c),
        Commands::Trajectory { args, session } => run_trajectory(&args, &session),
        Commands::Data(c) => manifest::data_files::run_data(&c),
        Commands::Attach(a) => Ok(legume_tui::cli::attach(&tui::TOOL, &a)?),
        Commands::Sessions => Ok(legume_tui::cli::sessions(&tui::TOOL)?),
        Commands::Kill(k) => Ok(legume_tui::cli::kill(&tui::TOOL, &k)?),
        Commands::Rename(r) => Ok(legume_tui::cli::rename(&tui::TOOL, &r)?),
        Commands::Mcp => mcp::run_mcp(),
    }
}
