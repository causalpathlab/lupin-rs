//! `lupin`: Label, Unfold, Place, Interpret, Narrate.

mod annotate;
mod annotate_manifest;
mod describe;
mod lineage_manifest;
mod plot;

use anyhow::Result;
use clap::{Parser, Subcommand};
use gene_text::cli::{run_knn_graph, run_qc, KnnGraphCmd, QcCmd};
use lineage::assoc::run::{run_assoc, AssocArgs};
use lineage::lineage::args::LineageArgs;
use lineage::lineage_plot::{run_lineage_plot, LineagePlotArgs};
use lineage::pseudotime::PseudotimeArgs;
// Resolved paths and manifest mutation stay in senna adapters.
use lineage_manifest::{run_lineage_from_manifest, run_pseudotime_from_manifest};

use annotate::{run_annotate, AnnotateCliArgs};
use describe::{run_describe, DescribeArgs};
use plot::scatter::{fit_plot, PlotArgs};
use plot::strand::{fit_plot_strand, PlotStrandArgs};
use plot::topic::{fit_plot_topic, PlotTopicArgs};

#[derive(Parser)]
#[command(
    name = "lupin",
    version,
    about = "Label, Unfold, Place, Interpret, Narrate — annotation, lineage, pseudotime, association, description."
)]
struct Cli {
    #[arg(short, long, global = true, help = "Verbose logging")]
    verbose: bool,
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
                      `--method enrichment` runs the senna topic/svd enrichment pipeline.\n\
                      `--method projection` runs senna co-embed projection.\n\
                      With `--feature-embedding` / `--cell-embedding` (or pinto parquets),\n\
                      it runs pinto-style ORA instead.\n\
                      `--method auto` (default) picks projection when embeddings resolve, else enrichment.\n\
                      Ontology-only follow-up: `--from` + `--obo` + `--label-cl` without markers."
    )]
    Annotate(AnnotateCliArgs),
    #[command(
        name = "lineage",
        about = "Geometry-first lineage and principal curves over a senna gem run"
    )]
    Lineage(LineageArgs),
    #[command(
        name = "lineage-plot",
        aliases = ["plot-lineage", "trajectory-plot"],
        about = "Publication-style figure of a senna lineage trajectory over its 2D embedding"
    )]
    LineagePlot(LineagePlotArgs),
    #[command(
        name = "dyn-assoc",
        about = "Bayesian between-branch modality contrast along a senna lineage"
    )]
    DynAssoc(AssocArgs),
    #[command(
        name = "pseudotime",
        about = "Monocle-style principal-graph pseudotime from a senna latent embedding"
    )]
    Pseudotime(PseudotimeArgs),
    #[command(
        name = "describe",
        about = "Short citation-checked sentence from annotate / lineage_annot evidence",
        long_about = "Builds structured evidence from `{from}.annot.parquet` (or argmax / lineage_annot).\n\
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
    env_logger::Builder::from_env(
        env_logger::Env::default().default_filter_or(if cli.verbose { "debug" } else { "info" }),
    )
    .init();
    match cli.cmd {
        Commands::TextQc(c) => run_qc(&c),
        Commands::WordGraph(c) => run_knn_graph(&c),
        Commands::Annotate(c) => run_annotate(&c),
        Commands::Lineage(c) => run_lineage_from_manifest(&c),
        Commands::LineagePlot(c) => run_lineage_plot(&c),
        Commands::DynAssoc(c) => run_assoc(&c),
        Commands::Pseudotime(c) => run_pseudotime_from_manifest(&c),
        Commands::Describe(c) => run_describe(&c),
        Commands::Plot(c) => fit_plot(&c),
        Commands::PlotTopic(c) => fit_plot_topic(&c),
        Commands::PlotStrand(c) => fit_plot_strand(&c),
    }
}
