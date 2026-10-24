//! [`super::without_dropped_evidence`]: the GAF rows every term is built from.

use super::*;
use std::io::Read;

fn gaf_line(gene: &str, go: &str, evidence: &str) -> String {
    let mut f = vec![""; 15];
    (f[0], f[1], f[2], f[4], f[6], f[11]) = ("UniProtKB", gene, gene, go, evidence, "protein");
    f.join("\t")
}

#[test]
fn untraced_and_no_data_annotations_are_left_out_once() {
    let dir = tempfile::tempdir().unwrap();
    let gaf = dir.path().join("x.gaf");
    let lines = [
        "!gaf-version: 2.2".to_string(),
        gaf_line("RPL10", "GO:1901740", "NAS"),
        gaf_line("CFLAR", "GO:1901740", "IDA"),
        gaf_line("GENE1", "GO:0008150", "ND"),
        gaf_line("GENE2", "GO:0006915", "IEA"),
    ];
    std::fs::write(&gaf, lines.join("\n") + "\n").unwrap();
    let gaf = gaf.to_string_lossy().into_owned();

    let kept = without_dropped_evidence(&gaf).unwrap();
    let mut text = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(&kept).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    let genes: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with('!'))
        .map(|l| l.split('\t').nth(2).unwrap())
        .collect();
    assert_eq!(genes, ["CFLAR", "GENE2"], "NAS and ND rows are left out");
    assert!(text.starts_with("!gaf-version"), "the header stays");
    assert_eq!(without_dropped_evidence(&gaf).unwrap(), kept, "made once");
    std::fs::remove_file(kept).unwrap();
}
