//! `lupin mcp` over the protocol, in process, for [`super`]. Tools that run
//! a child `lupin` are driven by hand (`lupin mcp`), not here: the test
//! binary is not lupin.

use super::*;
use crate::manifest::rounds::tests::{round_with_profile, round_with_profile_rows};
use std::io::Cursor;
use std::sync::{Arc, Mutex};

#[derive(Clone, Default)]
struct Sink(Arc<Mutex<Vec<u8>>>);

impl std::io::Write for Sink {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.lock().unwrap().write(buf)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Each call's result, in order.
fn calls(calls: &[(&str, Value)]) -> Vec<Value> {
    let input: String = calls
        .iter()
        .enumerate()
        .map(|(i, (name, args))| {
            let m = json!({ "jsonrpc": "2.0", "id": i, "method": "tools/call",
                            "params": { "name": name, "arguments": args } });
            format!("{m}\n")
        })
        .collect();
    let sink = Sink::default();
    legume_mcp::serve_on(
        Cursor::new(input),
        sink.clone(),
        ServerInfo::default(),
        Lupin,
    )
    .unwrap();
    let text = String::from_utf8(sink.0.lock().unwrap().clone()).unwrap();
    text.lines()
        .map(|l| serde_json::from_str::<Value>(l).unwrap()["result"].clone())
        .collect()
}

fn error_text(r: &Value) -> &str {
    assert_eq!(r["isError"], true, "{r}");
    r["content"][0]["text"].as_str().unwrap()
}

#[test]
fn every_tool_names_its_required_arguments() {
    let specs = tools();
    let names: BTreeSet<&str> = specs.iter().map(|t| t.name.as_str()).collect();
    assert_eq!(names.len(), specs.len(), "tool names are unique");
    for t in &specs {
        let props = t.input_schema["properties"].as_object().unwrap();
        for r in t.input_schema["required"].as_array().unwrap() {
            assert!(props.contains_key(r.as_str().unwrap()), "{}: {r}", t.name);
        }
    }
}

#[test]
fn a_session_is_labelled_by_its_outputs_name_cut_short_and_a_hash() {
    let long = format!("/results/{}/run.L1.T2_with_a_long_name", "a".repeat(200));
    let label = session_label(&long);
    let (name, hash) = label.rsplit_once('-').unwrap();
    assert_eq!(name, "run.L1", "{label}");
    assert_eq!(hash.len(), 8, "{label}");
    assert!(hash.chars().all(|c| c.is_ascii_hexdigit()), "{label}");
    let odd = session_label("/r/.x y");
    assert!(
        odd.starts_with("x-y-"),
        "a session's characters, no dot first: {odd}"
    );
    assert_eq!(
        label,
        session_label(&long),
        "the same output, the same label"
    );
    assert_ne!(label, session_label("/results/run.L2"));
}

#[test]
fn review_leaves_the_types_out_unless_asked() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path())
        .to_string_lossy()
        .into_owned();
    let r = calls(&[
        ("annotate_review", json!({ "round": round })),
        (
            "annotate_review",
            json!({ "round": round, "clusters": [1], "types": true }),
        ),
    ]);
    let all = r[0]["structuredContent"].as_object().unwrap();
    assert_eq!(all.len(), 3);
    assert_eq!(all["0"]["digest"]["label"], "CT1");
    assert!(all["0"].get("types").is_none());
    let one = r[1]["structuredContent"].as_object().unwrap();
    assert_eq!(one.keys().collect::<Vec<_>>(), ["1"]);
    assert!(one["1"]["types"].is_array());
}

#[test]
fn genes_show_what_sets_a_cluster_apart_and_its_labels_markers() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path())
        .to_string_lossy()
        .into_owned();
    let r = calls(&[(
        "annotate_genes",
        json!({ "round": round, "clusters": [1], "types": ["CT3"], "top": 2 }),
    )]);
    let c1 = &r[0]["structuredContent"]["1"];
    let specific: Vec<&str> = c1["specific"]
        .as_array()
        .unwrap()
        .iter()
        .map(|g| g[0].as_str().unwrap())
        .collect();
    assert_eq!(specific, ["GENE2", "GENE4"]);
    let ct2 = &c1["markers"]["CT2"];
    assert_eq!(
        (ct2["markers"].as_u64(), ct2["up"].as_u64()),
        (Some(1), Some(1))
    );
    assert_eq!(c1["markers"]["CT2"]["top"][0][0], "GENE2");
    assert_eq!(
        c1["markers"]["CT3"]["up"], 0,
        "a type asked for is shown too"
    );
}

#[test]
fn markers_match_genes_as_a_faba_file_names_them() {
    let dir = tempfile::tempdir().unwrap();
    let rows = [
        "GENE1/count/spliced",
        "ENSG2_gene2/count/spliced",
        "GENE3/count/spliced",
        "GENE4/count/spliced",
    ];
    let round = round_with_profile_rows(dir.path(), rows)
        .to_string_lossy()
        .into_owned();
    let r = calls(&[(
        "annotate_genes",
        json!({ "round": round, "clusters": [1], "types": ["CT2"] }),
    )]);
    let ct2 = &r[0]["structuredContent"]["1"]["markers"]["CT2"];
    assert_eq!(
        (ct2["measured"].as_u64(), ct2["up"].as_u64()),
        (Some(1), Some(1)),
        "{ct2}"
    );
    assert_eq!(
        r[0]["content"].as_array().unwrap().len(),
        1,
        "no warning: {}",
        r[0]
    );
}

#[test]
fn few_markers_measured_is_a_warning_with_what_did_not_match() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile_rows(dir.path(), ["X1", "X2", "X3", "X4"])
        .to_string_lossy()
        .into_owned();
    let r = calls(&[
        ("annotate_genes", json!({ "round": round, "clusters": [0] })),
        (
            "annotate_genes",
            json!({ "round": round, "channel": "intronic" }),
        ),
    ]);
    let content = r[0]["content"].as_array().unwrap();
    assert_eq!(content.len(), 2, "{}", r[0]);
    let warning = content[1]["text"].as_str().unwrap();
    for said in ["0 of the panel's 3 markers", "GENE1", "X1"] {
        assert!(warning.contains(said), "{said}: {warning}");
    }
    assert!(error_text(&r[1]).contains("spliced"), "{}", r[1]);
}

#[test]
fn a_coarse_label_is_measured_on_its_member_types_markers() {
    use crate::annotate::celltype_tree::{TreeSource, TypeGroup, TypeTree};
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path());
    let tree = TypeTree {
        source: TreeSource::CellOntology,
        release: None,
        label_cl: Default::default(),
        groups: vec![
            TypeGroup {
                name: "lineage23".into(),
                cl_id: None,
                members: vec!["CT2".into(), "CT3".into()],
            },
            TypeGroup {
                name: "empty".into(),
                cl_id: None,
                members: vec![],
            },
        ],
    };
    let (mut m, _) = run::RunManifest::load(&round).unwrap();
    let file = run::parent_dir(&round).join("tree.json");
    std::fs::write(&file, serde_json::to_string(&tree).unwrap()).unwrap();
    m.annotate.celltype_tree = Some("tree.json".into());
    m.save(&round).unwrap();
    let round = round.to_string_lossy().into_owned();
    let r = calls(&[(
        "annotate_genes",
        json!({ "round": round, "clusters": [1], "types": ["lineage23", "empty"] }),
    )]);
    let group = &r[0]["structuredContent"]["1"]["markers"]["lineage23"];
    assert_eq!(group["group_of"], json!(["CT2", "CT3"]), "{group}");
    let empty = &r[0]["structuredContent"]["1"]["markers"]["empty"];
    assert_eq!(empty["markers"], 0, "a group of no types: {empty}");
    let ct2 = &r[0]["structuredContent"]["1"]["markers"]["CT2"];
    assert!(
        group["markers"].as_u64() > ct2["markers"].as_u64(),
        "its members' markers together: {group} vs {ct2}"
    );
}

#[test]
fn why_answers_a_contrast_and_gives_what_a_decision_cites() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path())
        .to_string_lossy()
        .into_owned();
    let r = calls(&[
        (
            "annotate_why",
            json!({ "round": round, "cluster": 1, "label": "CT2", "over": "CT1" }),
        ),
        (
            "annotate_why",
            json!({ "round": round, "cluster": 1, "gene": "GENE1", "over": "GENE2" }),
        ),
        (
            "annotate_why",
            json!({ "round": round, "cluster": 1, "gene": "GENE1", "term": "x" }),
        ),
    ]);
    let label = &r[0]["structuredContent"];
    assert_eq!(label["question"]["kind"], "label", "{}", r[0]);
    assert!(label["favours"].is_string(), "{label}");
    assert!(label["verdict"]["fact"].is_u64(), "{label}");
    assert_eq!(label["cite"]["why"], label["question"]);
    assert_eq!(label["children"][0]["line"], "panel");
    let gene = &r[1]["structuredContent"];
    assert_eq!(gene["question"]["list"], "specific", "{gene}");
    assert!(error_text(&r[2]).contains("one of `label`, `gene` and `term`"));
}

#[test]
fn the_instructions_ask_in_contrasts_as_the_other_doors_do() {
    assert!(INSTRUCTIONS.contains(crate::manifest::why::ASKING));
    assert!(include_str!("../README.md").contains(crate::manifest::why::ASKING));
}

#[test]
fn outputs_are_never_replaced_and_overwrite_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path());
    let prefix = derive_out_prefix(&round.to_string_lossy());
    let r = calls(&[
        ("annotate_run", json!({ "from": round, "out": prefix })),
        (
            "annotate_run",
            json!({ "from": round, "options": ["--overwrite"] }),
        ),
        (
            "trajectory_run",
            json!({ "from": round, "statements": [
            { "from": "CT1", "to": "CT2", "relation": "precedes", "note": "a\tb" } ] }),
        ),
        (
            "trajectory_run",
            json!({ "from": round, "prior": "p.tsv", "statements": [
            { "from": "CT1", "to": "CT2", "relation": "precedes", "note": "n" } ] }),
        ),
        ("trajectory_read", json!({ "run": round })),
        (
            "annotate_relabel",
            json!({ "round": round, "decisions": [] }),
        ),
        (
            "annotate_relabel",
            json!({ "round": dir.path().join("nope"), "decisions": [{}] }),
        ),
    ]);
    assert!(error_text(&r[0]).contains("exists; not replacing it"));
    assert!(error_text(&r[1]).contains("--overwrite"));
    assert!(error_text(&r[2]).contains("tabs or newlines"));
    assert!(error_text(&r[3]).contains("not both"));
    assert!(error_text(&r[4]).contains("records no trajectory"));
    assert!(error_text(&r[5]).contains("no decisions"));
    assert!(error_text(&r[6]).contains("no such round"));
    assert!(!dir.path().join("nope.mcp_decisions.1.jsonl").exists());
}

#[test]
fn docs_list_and_print_the_write_ups() {
    let r = calls(&[
        ("docs", json!({})),
        ("docs", json!({ "topic": "grouping" })),
        ("docs", json!({ "topic": "nope" })),
    ]);
    assert!(r[0]["content"][0]["text"]
        .as_str()
        .unwrap()
        .contains("trajectory-plan"));
    assert!(r[1]["content"][0]["text"]
        .as_str()
        .unwrap()
        .starts_with('#'));
    error_text(&r[2]);
}

#[test]
fn the_readme_gives_an_agent_the_servers_own_instructions() {
    // For clients that do not pass the MCP instructions on to their model.
    let readme = include_str!("../README.md");
    assert!(
        readme.contains(INSTRUCTIONS),
        "README's agent block differs from INSTRUCTIONS"
    );
}

#[test]
fn docs_names_its_topics_so_an_agent_need_not_guess() {
    let docs = tools().into_iter().find(|t| t.name == "docs").unwrap();
    let topics = docs.input_schema["properties"]["topic"]["enum"]
        .as_array()
        .unwrap()
        .clone();
    assert_eq!(
        topics,
        crate::docs::slugs()
            .into_iter()
            .map(Value::from)
            .collect::<Vec<_>>()
    );
    assert!(topics.contains(&json!("annotation")));
}

#[test]
fn a_figure_of_a_run_with_no_trajectory_is_refused_before_anything_runs() {
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path())
        .to_string_lossy()
        .into_owned();
    let r = calls(&[
        ("trajectory_figure", json!({ "run": round })),
        (
            "annotate_heatmap",
            json!({ "round": round, "out": "x.txt" }),
        ),
        (
            "annotate_heatmap",
            json!({ "round": round, "channel": "intronic" }),
        ),
    ]);
    assert!(error_text(&r[2]).contains("intronic"), "{}", r[2]);
    assert!(
        error_text(&r[0]).contains("trajectory.pseudotime"),
        "{}",
        r[0]
    );
    assert!(error_text(&r[1]).contains(".pdf, .svg or .png"), "{}", r[1]);
    let left: Vec<_> = std::fs::read_dir(dir.path().join("r0"))
        .unwrap()
        .flatten()
        .map(|e| e.file_name())
        .filter(|n| n.to_string_lossy().contains("figure"))
        .collect();
    assert!(left.is_empty(), "{left:?}");
}
