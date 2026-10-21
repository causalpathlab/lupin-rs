//! `lupin mcp`'s background steps run in a detached session of lupin's own
//! screen: a person attaches to it, the agent follows it as a job, it
//! outlives the server, and a pass that succeeds is recorded.

#![cfg(unix)]

use legume_numeric::matrix::dense_mat_io::Mat;
use legume_numeric::matrix::traits::IoOps;
use legume_tui::testing::Private;
use legume_tui::Size;
use serde_json::{json, Value};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Duration;

const LUPIN: &str = env!("CARGO_BIN_EXE_lupin");

fn private() -> Private {
    Private::new(LUPIN).patience(Duration::from_secs(60))
}

/// A planted run in `dir`: six genes (GENE0..GENE2 load on axis 0, GENE3..GENE5
/// on axis 1), twenty cells out on the two axes in clusters 0 and 1, and a
/// panel of CT1 (GENE0..GENE2) and CT2 (GENE3..GENE5). The manifest and the
/// panel.
fn run(dir: &Path) -> (PathBuf, PathBuf) {
    let at = |n: &str| dir.join(n).to_string_lossy().into_owned();
    let genes: Vec<Box<str>> = (0..6).map(|g| format!("GENE{g}").into()).collect();
    let mut rho = Mat::zeros(6, 2);
    for g in 0..6 {
        rho[(g, g / 3)] = 1.0 + 0.1 * g as f32;
        rho[(g, 1 - g / 3)] = 0.05 * g as f32;
    }
    let cells: Vec<Box<str>> = (0..20).map(|i| format!("c{i:02}").into()).collect();
    let mut z = Mat::zeros(20, 2);
    let mut ids = Mat::zeros(20, 1);
    for i in 0..20 {
        z[(i, i / 10)] = 2.0 + 0.1 * (i % 10) as f32;
        z[(i, 1 - i / 10)] = 0.2 * ((i % 3) as f32 - 1.0);
        ids[(i, 0)] = (i / 10) as f32;
    }
    rho.to_parquet_with_names(
        &at("run.feature_embedding.parquet"),
        (Some(&genes), Some("gene")),
        None,
    )
    .unwrap();
    z.to_parquet_with_names(
        &at("run.cell_embedding.parquet"),
        (Some(&cells), Some("cell")),
        None,
    )
    .unwrap();
    let col: Vec<Box<str>> = vec!["cluster".into()];
    ids.to_parquet_with_names(
        &at("run.clusters.parquet"),
        (Some(&cells), Some("cell")),
        Some(&col),
    )
    .unwrap();
    let manifest = json!({
        "version": 2, "kind": "bge", "prefix": "run",
        "outputs": {"feature_embedding": "run.feature_embedding.parquet",
                    "cell_embedding": "run.cell_embedding.parquet"},
        "cluster": {"clusters": "run.clusters.parquet"}
    });
    std::fs::write(at("run.senna.json"), manifest.to_string()).unwrap();
    let panel = dir.join("panel.tsv");
    let rows: String = (0..6)
        .map(|g| format!("GENE{g}\tCT{}\n", 1 + g / 3))
        .collect();
    std::fs::write(&panel, rows).unwrap();
    (dir.join("run.senna.json"), panel)
}

/// A `lupin mcp` in `project`, its sessions in `sessions`, given `calls`
/// (tool, arguments) after the handshake: each call's result.
fn server(sessions: &Private, project: &Path, calls: &[(&str, Value)]) -> Vec<Value> {
    let mut child = Command::new(LUPIN)
        .arg("mcp")
        .current_dir(project)
        .env_remove("LEGUME_MCP_RECORDS")
        .env_remove("LEGUME_MCP_PROJECT")
        .env("LEGUME_MCP_WATCH", project.join(".watch"))
        .env(legume_tui::DIR_VAR, sessions.path())
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    let mut stdin = child.stdin.take().unwrap();
    let hello = json!({"jsonrpc": "2.0", "id": 0, "method": "initialize", "params": {}});
    writeln!(stdin, "{hello}").unwrap();
    for (id, (tool, args)) in calls.iter().enumerate() {
        let msg = json!({"jsonrpc": "2.0", "id": id + 1, "method": "tools/call",
                         "params": {"name": tool, "arguments": args}});
        writeln!(stdin, "{msg}").unwrap();
    }
    drop(stdin);
    let results = BufReader::new(child.stdout.take().unwrap())
        .lines()
        .map(|l| serde_json::from_str::<Value>(&l.unwrap()).unwrap())
        .filter(|m| m.get("id").is_some_and(|id| id != 0))
        .map(|m| m["result"].clone())
        .collect();
    assert!(child.wait().unwrap().success());
    results
}

fn structured(r: &Value) -> &Value {
    &r["structuredContent"]
}

#[test]
fn a_background_pass_runs_on_lupin_s_screen_and_is_recorded() {
    let sessions = private();
    let tmp = tempfile::tempdir().unwrap();
    let project = tmp.path();
    let (manifest, panel) = run(project);
    let out = project.join("pass");
    // Held on a panel that is a pipe no one writes yet, so the session is
    // there to attach to; the server ends meanwhile.
    let pipe = project.join("pipe.tsv");
    assert!(Command::new("mkfifo")
        .arg(&pipe)
        .status()
        .unwrap()
        .success());
    let pass = json!({ "from": manifest, "markers": pipe, "out": out, "background": true });
    let results = server(&sessions, project, &[("annotate_run", pass)]);
    let started = structured(&results[0]);
    let session = started["session"]
        .as_str()
        .unwrap_or_else(|| panic!("{}", results[0]))
        .to_string();
    assert_eq!(session, "lupin-annotate-pass", "{started}");
    sessions.wait_listed(&session);
    let size = Size {
        rows: 40,
        cols: 120,
        width: 1200,
        height: 800,
    };
    // Attached: the session takes the terminal (nothing is drawn while
    // lupin waits on its panel); detached again, it runs on.
    let mut t = sessions.term_sized(&["attach", &session], &[], size);
    sessions.wait_listed(&format!("{session})\tattached"));
    t.type_(b"\x1b\x04");
    assert_eq!(t.exit_code(), 0);
    sessions.wait_listed(&format!("{session})\tdetached"));
    // The panel comes: through the pipe to lupin waiting on it, and as a
    // file in its place for what reads it after. The pass runs on, and
    // ends with no server there.
    let mut waiting = std::fs::OpenOptions::new().write(true).open(&pipe).unwrap();
    std::fs::remove_file(&pipe).unwrap();
    std::fs::copy(&panel, &pipe).unwrap();
    waiting.write_all(&std::fs::read(&panel).unwrap()).unwrap();
    drop(waiting);
    sessions.wait_said(&["sessions"], "no live lupin sessions");

    // The next server takes it up: done, and recorded.
    let results = server(&sessions, project, &[("job_status", json!({}))]);
    let jobs = &structured(&results[0])["jobs"];
    assert_eq!(jobs[0]["session"], session.as_str(), "{jobs}");
    assert_eq!(jobs[0]["state"], "done", "{jobs}");
    let index = std::fs::read_to_string(project.join("results/INDEX.tsv")).unwrap();
    assert!(index.contains("pass"), "{index}");
}
