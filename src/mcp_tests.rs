//! `lupin mcp` over the protocol, in process, for [`super`]. Tools that run
//! a child `lupin` are driven by hand (`lupin mcp`), not here: the test
//! binary is not lupin.

use super::*;
use crate::manifest::rounds::tests::round_with_profile;
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
fn a_coarse_label_is_measured_on_its_member_types_markers() {
    use crate::annotate::celltype_tree::{TreeSource, TypeGroup, TypeTree};
    let dir = tempfile::tempdir().unwrap();
    let round = round_with_profile(dir.path());
    let tree = TypeTree {
        source: TreeSource::CellOntology,
        release: None,
        label_cl: Default::default(),
        groups: vec![TypeGroup {
            name: "lineage23".into(),
            cl_id: None,
            members: vec!["CT2".into(), "CT3".into()],
        }],
    };
    let (mut m, _) = run::RunManifest::load(&round).unwrap();
    let file = run::parent_dir(&round).join("tree.json");
    std::fs::write(&file, serde_json::to_string(&tree).unwrap()).unwrap();
    m.annotate.celltype_tree = Some("tree.json".into());
    m.save(&round).unwrap();
    let round = round.to_string_lossy().into_owned();
    let r = calls(&[(
        "annotate_genes",
        json!({ "round": round, "clusters": [1], "types": ["lineage23"] }),
    )]);
    let group = &r[0]["structuredContent"]["1"]["markers"]["lineage23"];
    assert_eq!(group["group_of"], json!(["CT2", "CT3"]), "{group}");
    let ct2 = &r[0]["structuredContent"]["1"]["markers"]["CT2"];
    assert!(
        group["markers"].as_u64() > ct2["markers"].as_u64(),
        "its members' markers together: {group} vs {ct2}"
    );
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
