//! `lupin`: Label, Unfold, Place, Interpret, Narrate.

mod annotate;
mod annotate_cmd;
mod cell_labels;
mod describe;
mod docs;
mod gene_text;
mod manifest;
mod marker_embedding;
mod plot;
mod trajectory;
mod tui;

use crate::annotate_cmd::{run_annotate, AnnotateCliArgs};
use crate::gene_text::cli::{run_knn_graph, run_qc, KnnGraphCmd, QcCmd};
use anyhow::Result;
use clap::{Parser, Subcommand};
use describe::{run_describe, DescribeArgs};
use plot::scatter::{fit_plot, PlotArgs};
use plot::strand::{fit_plot_strand, PlotStrandArgs};
use plot::topic::{fit_plot_topic, PlotTopicArgs};
use trajectory::run::{run_trajectory, TrajectoryArgs};

#[derive(Parser)]
#[command(
    name = "lupin",
    version,
    about = "Label, Unfold, Place, Interpret, Narrate — text graphs, cell-type annotation, review, plotting and description."
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
        long_about = "Unified annotation entry point.\n\
                      \n\
                      `--method enrichment` scores the marker panel on per-cluster expression (topic/svd runs).\n\
                      `--method projection` projects cells onto the run's co-embedded gene space.\n\
                      `--feature-embedding` + `--cell-embedding` run projection on those files directly.\n\
                      `--method auto` (default) picks projection when a co-embedding exists, else enrichment.\n\
                      Ontology-only follow-up: `--from` + `--obo` + `--label-cl` without markers."
    )]
    Annotate(AnnotateCliArgs),
    #[command(
        name = "describe",
        about = "Short citation-checked sentence from annotate evidence",
        long_about = "Builds structured evidence from `{from}.annot.parquet` (or argmax).\n\
                      Optionally fishes per-cluster keywords from a `word-graph` prefix\n\
                      (`feature_word.edges` over each cluster's markers;\n\
                      `--feature-embedding` adds nearest-neighbour genes first).\n\
                      When `{from}.cluster_term_q.parquet` is present,\n\
                      a second FDR-significant contender is named if one exists.\n\
                      Writes `{out}.describe.json` and `{out}.describe.md`.\n\
                      \n\
                      The composer never invents labels: sentences are citation-checked.\n\
                      Default composer is a citation-checked template (candle decoder TBD)."
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
        about = "Order cells along a supervised prior over their cell types: Cell Ontology and precedence statements, checked against the kNN graph, then diffusion pseudotime"
    )]
    Trajectory(TrajectoryArgs),
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
        name = "plot",
        about = "Publication-quality scatter over a senna layout embedding",
        long_about = "Rasterized scatter with vector labels over a transparent background.\n\
                      Preferred: `lupin plot --from {prefix}.senna.json` after `senna layout`.\n\
                      Explicit flags override manifest defaults."
    )]
    Plot(PlotArgs),
    #[command(
        name = "plot-topic",
        about = "Structure-bar and dictionary plots from a senna topic run"
    )]
    PlotTopic(PlotTopicArgs),
    #[command(
        name = "plot-strand",
        about = "Watson/Crick mirrored genomic-activity ideograms per cell type"
    )]
    PlotStrand(PlotStrandArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    if cli.overwrite {
        manifest::run::allow_overwrite();
    }
    if matches!(&cli.cmd, Commands::Annotate(c) if c.tui) {
        tui::init_logger(cli.verbose);
    } else {
        env_logger::Builder::from_env(
            env_logger::Env::default().default_filter_or(if cli.verbose {
                "debug"
            } else {
                "info"
            }),
        )
        .init();
    }
    match cli.cmd {
        Commands::TextQc(c) => run_qc(&c),
        Commands::WordGraph(c) => run_knn_graph(&c),
        Commands::Annotate(c) => run_annotate(&c),
        Commands::Describe(c) => run_describe(&c),
        Commands::Review(c) => manifest::rounds::run_review(&c),
        Commands::Relabel(c) => manifest::rounds::run_relabel(&c),
        Commands::Ask(c) => manifest::ask::run_ask(&c),
        Commands::Docs(c) => docs::run_docs(&c),
        Commands::Trajectory(c) => run_trajectory(&c),
        Commands::Data(c) => manifest::data_files::run_data(&c),
        Commands::Plot(c) => fit_plot(&c),
        Commands::PlotTopic(c) => fit_plot_topic(&c),
        Commands::PlotStrand(c) => fit_plot_strand(&c),
    }
}
