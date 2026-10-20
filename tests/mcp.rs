//! `lupin mcp` end to end: the server binary over stdio, with the tools that
//! run a child `lupin` (preview and relabel), the project's records, and a
//! replay of the day's `runall.sh` with this `lupin` on the path.

use graph_embedding_util::type_annotation::write_cell_labels;
use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::traits::IoOps;
use serde_json::{json, Value};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

const LUPIN: &str = env!("CARGO_BIN_EXE_lupin");

/// A first round in `dir`: cells a–d in clusters 0, 0, 1, 2 labelled CT1,
/// CT1, CT2, CT3. Its manifest.
fn round(dir: &Path) -> PathBuf {
    let cells: Vec<Box<str>> = ["a", "b", "c", "d"].map(Box::from).into();
    let mut ids = Mat::zeros(4, 1);
    for (i, k) in [0.0, 0.0, 1.0, 2.0].into_iter().enumerate() {
        ids[(i, 0)] = k;
    }
    let col: Vec<Box<str>> = vec!["cluster".into()];
    let clusters = dir.join("run.clusters.parquet");
    ids.to_parquet_with_names(
        &clusters.to_string_lossy(),
        (Some(&cells), Some("cell")),
        Some(&col),
    )
    .unwrap();
    let labels: Vec<Box<str>> = ["CT1", "CT1", "CT2", "CT3"].map(Box::from).into();
    let prefix = dir.join("run").to_string_lossy().into_owned();
    write_cell_labels(&prefix, &cells, &labels, &[0.9, 0.8, 0.7, 0.6]).unwrap();
    let argmax = std::fs::read_dir(dir)
        .unwrap()
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .find(|n| n.starts_with("run.argmax"))
        .expect("the cell labels");
    let manifest = dir.join("run.senna.json");
    let m = json!({
        "version": 2, "kind": "topic", "prefix": "run",
        "cluster": { "clusters": "run.clusters.parquet" },
        "annotate": { "argmax": argmax },
    });
    std::fs::write(&manifest, m.to_string()).unwrap();
    manifest
}

struct Server {
    child: Child,
    stdin: ChildStdin,
    stdout: BufReader<ChildStdout>,
    id: u64,
}

impl Server {
    fn start(project: &Path) -> Self {
        let mut child = Command::new(LUPIN)
            .arg("mcp")
            .current_dir(project)
            .env_remove("LEGUME_MCP_JOURNAL")
            .env_remove("LEGUME_MCP_PROJECT")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .unwrap();
        let stdin = child.stdin.take().unwrap();
        let stdout = BufReader::new(child.stdout.take().unwrap());
        let mut s = Self {
            child,
            stdin,
            stdout,
            id: 0,
        };
        s.request("initialize", json!({ "protocolVersion": "2025-06-18" }));
        s
    }

    fn request(&mut self, method: &str, params: Value) -> Value {
        self.id += 1;
        let m = json!({ "jsonrpc": "2.0", "id": self.id, "method": method, "params": params });
        writeln!(self.stdin, "{m}").unwrap();
        loop {
            let mut line = String::new();
            assert!(
                self.stdout.read_line(&mut line).unwrap() > 0,
                "the server ended"
            );
            let v: Value = serde_json::from_str(&line).unwrap();
            if v["id"] == json!(self.id) {
                return v;
            }
        }
    }

    /// A tool's result: its structured content, or the error text.
    fn call(&mut self, name: &str, args: Value) -> Result<Value, String> {
        let r = self.request("tools/call", json!({ "name": name, "arguments": args }));
        let r = &r["result"];
        if r["isError"] == true {
            return Err(r["content"][0]["text"].as_str().unwrap_or_default().into());
        }
        Ok(r.get("structuredContent").cloned().unwrap_or(Value::Null))
    }

    fn stop(mut self) {
        drop(self.stdin);
        self.child.wait().unwrap();
    }
}

/// Cluster `id`'s label in `round`, as `annotate_review` shows it.
fn label(s: &mut Server, round: &str, id: u32) -> String {
    let r = s
        .call("annotate_review", json!({ "round": round }))
        .unwrap();
    r[id.to_string()]["digest"]["label"]
        .as_str()
        .unwrap()
        .into()
}

#[test]
fn an_agents_relabel_is_recorded_and_replays() {
    let project = tempfile::tempdir().unwrap();
    let data = project.path().join("data");
    std::fs::create_dir_all(&data).unwrap();
    let r0 = round(&data).to_string_lossy().into_owned();
    let mut s = Server::start(project.path());

    let init = s.request("tools/list", json!({}));
    let names: Vec<&str> = init["result"]["tools"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| t["name"].as_str().unwrap())
        .collect();
    for t in [
        "annotate_run",
        "annotate_relabel",
        "trajectory_run",
        "notebook",
    ] {
        assert!(names.contains(&t), "{t} in {names:?}");
    }
    s.call(
        "notebook",
        json!({ "kind": "plan", "text": "Relabel cluster 2." }),
    )
    .unwrap();

    let decisions =
        json!([{ "cluster": 2, "action": "label", "label": "CT9", "rationale": "test" }]);
    let preview = s
        .call(
            "annotate_preview",
            json!({ "round": r0, "decisions": decisions }),
        )
        .unwrap();
    assert_eq!(preview["clusters"]["2"]["label_after"], "CT9", "{preview}");
    let r1 = s
        .call(
            "annotate_relabel",
            json!({ "round": r0, "decisions": decisions }),
        )
        .unwrap()["round"]
        .as_str()
        .unwrap()
        .to_string();
    assert!(r1.ends_with("run.r1.senna.json"), "{r1}");
    assert_eq!(label(&mut s, &r1, 2), "CT9");
    assert_eq!(
        label(&mut s, &r0, 2),
        "CT3",
        "the old round is left as it was"
    );
    let again = s.call(
        "annotate_relabel",
        json!({ "round": r0, "decisions": decisions }),
    );
    assert!(again.unwrap_err().contains("not the latest round"));
    let failed = s.call(
        "annotate_run",
        json!({ "from": r0, "markers": "no-such-panel.tsv" }),
    );
    assert!(!failed.unwrap_err().is_empty(), "a failed child says why");
    s.stop();

    // The records.
    let results = project.path().join("results");
    let day = std::fs::read_dir(&results)
        .unwrap()
        .flatten()
        .map(|e| e.path())
        .find(|p| p.is_dir())
        .expect("today's folder");
    let script = std::fs::read_to_string(day.join("runall.sh")).unwrap();
    assert!(
        script.contains("if [ ! -e data/run.r1.senna.json ]; then"),
        "{script}"
    );
    assert!(
        script.contains(
            "lupin relabel -f data/run.senna.json -d data/run.mcp_decisions.1.jsonl -o data/run.r1"
        ),
        "{script}"
    );
    assert!(script.contains("# did not succeed ("), "{script}");
    assert!(
        !script.contains(&*project.path().to_string_lossy()),
        "paths are relative: {script}"
    );
    let index = std::fs::read_to_string(results.join("INDEX.tsv")).unwrap();
    assert!(index.contains("\tlupin\tannotate_relabel\tdata/run.r1.senna.json\trelabel run.senna.json with 1 decision(s)"), "{index}");
    let notebook = std::fs::read_to_string(day.join("notebook.md")).unwrap();
    assert!(
        notebook.contains("· lupin · plan\n\nRelabel cluster 2."),
        "{notebook}"
    );

    // Replay: without the new round and its decisions file, the script
    // writes both again, with this lupin first on the path.
    for f in std::fs::read_dir(&data).unwrap().flatten() {
        let n = f.file_name().to_string_lossy().into_owned();
        if n.starts_with("run.r1") || n.contains("mcp_decisions") {
            std::fs::remove_file(f.path()).unwrap();
        }
    }
    let bin = Path::new(LUPIN).parent().unwrap();
    let path = format!(
        "{}:{}",
        bin.display(),
        std::env::var("PATH").unwrap_or_default()
    );
    let ran = Command::new("bash")
        .arg(day.join("runall.sh"))
        .env("PATH", path)
        .status()
        .unwrap();
    assert!(ran.success());
    let mut s = Server::start(project.path());
    assert_eq!(
        label(&mut s, &r1, 2),
        "CT9",
        "the replay made the same round"
    );
    s.stop();
}
