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
    /// Its background steps' sessions, ended with it.
    _sessions: legume_tui::testing::Private,
    child: Child,
    stdin: ChildStdin,
    stdout: BufReader<ChildStdout>,
    id: u64,
}

impl Server {
    fn start(project: &Path) -> Self {
        let sessions = legume_tui::testing::Private::new(LUPIN);
        let mut child = Command::new(LUPIN)
            .env(legume_tui::DIR_VAR, sessions.path())
            .arg("mcp")
            .current_dir(project)
            .env_remove("LEGUME_MCP_RECORDS")
            .env_remove("LEGUME_MCP_PROJECT")
            .env("LEGUME_MCP_WATCH", project.join(".watch"))
            // The user's own data files stay out of reach.
            .env("LUPIN_CONFIG_DIR", project.join(".config"))
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .unwrap();
        let stdin = child.stdin.take().unwrap();
        let stdout = BufReader::new(child.stdout.take().unwrap());
        let mut s = Self {
            _sessions: sessions,
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
    // A folder of the repository's own, not the system's temporary one: on
    // macOS that is reached through a link (/var to /private/var), so the
    // project's paths would not read as relative to where it is.
    let local = Path::new(env!("CARGO_MANIFEST_DIR")).join("results");
    std::fs::create_dir_all(&local).unwrap();
    let project = tempfile::tempdir_in(&local).unwrap();
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
    // A person following the server sees the round the agent made.
    let live = legume_mcp::watch::live_in(&project.path().join(".watch"));
    let view = live.iter().find_map(|l| legume_mcp::watch::state(&l.dir));
    assert_eq!(
        view.as_ref().and_then(|v| v["round"].as_str()),
        Some(r1.as_str())
    );
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
    // As a background job: it returns at once, and waiting tells how it went.
    let job = s
        .call(
            "annotate_run",
            json!({ "from": r0, "markers": "no-such-panel.tsv", "background": true }),
        )
        .unwrap();
    let id = job["job"].as_str().unwrap().to_string();
    assert!(job["round"].as_str().is_some(), "{job}");
    let ended = s
        .call("job_wait", json!({ "id": id, "seconds": 60 }))
        .unwrap();
    assert_eq!(ended["state"], "failed", "{ended}");
    assert!(
        ended["error"].as_str().unwrap().contains("no-such-panel"),
        "{ended}"
    );
    // While a job runs (held here on a panel that is a pipe no one writes),
    // its name is not chosen again, and no step writes over it: a
    // trajectory at its `out` is refused before its statements are written.
    let pipe = data.join("panel.tsv");
    assert!(Command::new("mkfifo")
        .arg(&pipe)
        .status()
        .unwrap()
        .success());
    let jobs: Vec<Value> = (0..2)
        .map(|_| {
            s.call(
                "annotate_run",
                json!({ "from": r0, "markers": pipe, "background": true }),
            )
            .unwrap()
        })
        .collect();
    assert_ne!(jobs[0]["round"], jobs[1]["round"], "{jobs:?}");
    let out = jobs[0]["round"]
        .as_str()
        .unwrap()
        .trim_end_matches(".senna.json");
    let statements = json!([{ "from": "CT1", "to": "CT2", "relation": "precedes", "note": "x" }]);
    let refused = s.call(
        "trajectory_run",
        json!({ "from": r0, "out": out, "statements": statements }),
    );
    assert!(
        refused.unwrap_err().contains("is being written by"),
        "a held name is refused"
    );
    assert!(!Path::new(&format!("{out}.precedence.tsv")).exists());
    for j in &jobs {
        s.call("job_cancel", json!({ "id": j["job"] })).unwrap();
        let ended = s
            .call("job_wait", json!({ "id": j["job"], "seconds": 60 }))
            .unwrap();
        assert_eq!(ended["state"], "cancelled", "{ended}");
    }
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
            "lupin relabel -f data/run.senna.json -d data/run.mcp_decisions.1.jsonl --next --expect data/run.r1"
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

/// [`round`] with a profile whose rows a faba file names
/// (`GENE1/count/spliced` …) and a panel naming the genes bare.
fn round_with_profile(dir: &Path) -> PathBuf {
    let manifest = round(dir);
    let rows: Vec<Box<str>> = (1..=3)
        .map(|g| format!("GENE{g}/count/spliced").into())
        .collect();
    let cols: Vec<Box<str>> = ["K0", "K1", "K2"].map(Box::from).into();
    let mut m = Mat::from_element(3, 3, 1.0);
    for g in 0..3 {
        m[(g, g)] = 50.0;
    }
    m.to_parquet_with_names(
        &dir.join("run.cluster_expression.parquet").to_string_lossy(),
        (Some(&rows), Some("gene")),
        Some(&cols),
    )
    .unwrap();
    std::fs::write(
        dir.join("markers.tsv"),
        "gene\tcelltype\nGENE1\tCT1\nGENE2\tCT2\nGENE3\tCT3\n",
    )
    .unwrap();
    let mut v: Value = serde_json::from_str(&std::fs::read_to_string(&manifest).unwrap()).unwrap();
    v["annotate"]["cluster_expression"] = json!("run.cluster_expression.parquet");
    v["annotate"]["expression_clusters"] = json!("run.clusters.parquet");
    v["annotate"]["markers"] = json!("markers.tsv");
    std::fs::write(&manifest, v.to_string()).unwrap();
    manifest
}

#[test]
fn a_heatmap_is_drawn_returned_and_recorded() {
    let local = Path::new(env!("CARGO_MANIFEST_DIR")).join("results");
    std::fs::create_dir_all(&local).unwrap();
    let project = tempfile::tempdir_in(&local).unwrap();
    let data = project.path().join("data");
    std::fs::create_dir_all(&data).unwrap();
    let r0 = round_with_profile(&data).to_string_lossy().into_owned();
    let mut s = Server::start(project.path());
    let r = s.request(
        "tools/call",
        json!({ "name": "annotate_heatmap", "arguments": { "round": r0, "value": "log2fc" } }),
    );
    let r = &r["result"];
    assert_ne!(r["isError"], true, "{r}");
    let image = r["content"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["type"] == "image");
    assert!(image.is_some(), "the image comes back: {r}");
    let got = &r["structuredContent"];
    let figure = got["figure"].as_str().unwrap();
    assert!(figure.ends_with("run.heatmap.pdf"), "{figure}");
    assert!(Path::new(figure).is_file());
    assert_eq!(
        got["rows"].as_array().unwrap().len(),
        3,
        "faba rows match: {got}"
    );
    assert_eq!(got["rows"][0]["gene"], "GENE1");
    let refused = s.call("annotate_heatmap", json!({ "round": r0, "out": "h.txt" }));
    assert!(refused.unwrap_err().contains(".pdf, .svg or .png"));
    s.stop();

    let results = project.path().join("results");
    let index = std::fs::read_to_string(results.join("INDEX.tsv")).unwrap();
    assert!(
        index.contains("run.heatmap.pdf\tmarker heatmap (log2fc) of run.senna.json"),
        "{index}"
    );
    let day = Path::new(figure).parent().unwrap();
    let script = std::fs::read_to_string(day.join("runall.sh")).unwrap();
    assert!(script.contains("lupin figure heatmap"), "{script}");
}

/// The Cell Ontology for [`round_with_profile`]: CT1 and CT2 under group1.
const OBO: &str = "format-version: 1.2

[Term]
id: CL:0
name: root cell

[Term]
id: CL:1
name: group1
subset: cellxgene_subset
is_a: CL:0

[Term]
id: CL:10
name: CT1
is_a: CL:1

[Term]
id: CL:12
name: CT2
is_a: CL:1
";

#[test]
fn a_mixed_label_is_named_and_the_export_places_it_recorded() {
    let local = Path::new(env!("CARGO_MANIFEST_DIR")).join("results");
    std::fs::create_dir_all(&local).unwrap();
    let project = tempfile::tempdir_in(&local).unwrap();
    let data = project.path().join("data");
    std::fs::create_dir_all(&data).unwrap();
    let r0 = round_with_profile(&data).to_string_lossy().into_owned();
    let obo = data.join("cl.obo");
    std::fs::write(&obo, OBO).unwrap();
    let obo = obo.to_string_lossy().into_owned();
    let mut s = Server::start(project.path());
    let mix = json!([{ "cluster": 0, "action": "label", "label": "CT1+CT2",
                       "parts": ["CT1", "CT2"], "rationale": "both" }]);
    let short = json!([{ "cluster": 0, "action": "label", "label": "CT1+",
                         "parts": ["CT1"], "rationale": "one" }]);
    let refused = s.call(
        "annotate_relabel",
        json!({ "round": r0, "decisions": short }),
    );
    assert!(
        refused.unwrap_err().contains("two or more"),
        "a mix has parts"
    );
    let r1 = s
        .call("annotate_relabel", json!({ "round": r0, "decisions": mix }))
        .unwrap()["round"]
        .as_str()
        .unwrap()
        .to_string();
    let mixed = std::fs::read_to_string(data.join("lupin/mixed_labels.tsv")).unwrap();
    assert!(mixed.contains("CT1+CT2\tCT1\tCT2"), "{mixed}");
    let got = s
        .call("annotate_export", json!({ "round": r1, "obo": obo }))
        .unwrap();
    let cells = got["cells"].as_str().unwrap();
    assert!(cells.ends_with("run.r1.cell_annotation.parquet"), "{got}");
    assert!(Path::new(got["markers"].as_str().unwrap()).is_file());
    let ct = &got["labels"]["CT1+CT2"];
    assert_eq!(ct["cells"], 2, "{got}");
    assert_eq!(ct["cl_id"], "CL:10|CL:12", "its parts' terms: {got}");
    assert!(
        ct["lineage"]
            .as_str()
            .unwrap()
            .ends_with("group1 > CT1+CT2"),
        "under their common ancestor: {got}"
    );
    let again = s
        .call("annotate_export", json!({ "round": r1, "obo": obo }))
        .unwrap();
    assert_ne!(again["cells"], got["cells"], "never over an export");
    let family = s.call("annotate_family", json!({ "from": r0 })).unwrap();
    assert_eq!(family["members"].as_array().unwrap().len(), 2, "{family}");
    s.stop();
    let index = std::fs::read_to_string(project.path().join("results/INDEX.tsv")).unwrap();
    assert!(index.contains("\tlupin\tannotate_export\t"), "{index}");
    let day = Path::new(cells).parent().unwrap();
    let script = std::fs::read_to_string(day.join("runall.sh")).unwrap();
    assert!(
        script.contains("lupin export --from=data/run.r1.senna.json --out=results/"),
        "{script}"
    );
}

/// A trajectory run in `dir` on 40 cells of CT1 and CT2 over two lineages:
/// pseudotime along DC1, DC3 parting the lineages, and no layout.
fn trajectory(dir: &Path) -> PathBuf {
    use legume_numeric::matrix::parquet::{write_named_table, Column};
    let n = 40;
    let cells: Vec<Box<str>> = (0..n).map(|i| format!("c{i}").into()).collect();
    let types: Vec<Box<str>> = (0..n)
        .map(|i| if i < 20 { "CT1" } else { "CT2" }.into())
        .collect();
    let pt: Vec<f32> = (0..n).map(|i| i as f32 / n as f32).collect();
    let component = vec![0i32; n];
    let l0: Vec<f32> = (0..n).map(|i| f32::from(u8::from(i < 20))).collect();
    let l1: Vec<f32> = l0.iter().map(|w| 1.0 - w).collect();
    let at = |name: &str| dir.join(name).to_string_lossy().into_owned();
    write_named_table(
        &at("t.cell_pseudotime.parquet"),
        "cell",
        &cells,
        &[
            ("pseudotime".into(), Column::F32(&pt)),
            ("type".into(), Column::Str(&types)),
            ("component".into(), Column::I32(&component)),
            ("L0".into(), Column::F32(&l0)),
            ("L1".into(), Column::F32(&l1)),
        ],
    )
    .unwrap();
    let mut d = Mat::zeros(n, 4);
    for i in 0..n {
        d[(i, 0)] = 1.0;
        d[(i, 1)] = i as f32;
        d[(i, 2)] = (i % 7) as f32;
        d[(i, 3)] = f32::from(u8::from(i >= 20)) + 0.01 * (i % 3) as f32;
    }
    let dcs: Vec<Box<str>> = (0..4).map(|j| format!("DC{j}").into()).collect();
    d.to_parquet_with_names(
        &at("t.diffusion.parquet"),
        (Some(&cells), Some("cell")),
        Some(&dcs),
    )
    .unwrap();
    let manifest = dir.join("t.senna.json");
    std::fs::write(
        &manifest,
        r#"{"version":2,"kind":"topic","prefix":"t",
            "trajectory":{"pseudotime":"t.cell_pseudotime.parquet","diffusion":"t.diffusion.parquet"}}"#,
    )
    .unwrap();
    manifest
}

/// The images in a tool call's result.
fn images(r: &Value) -> usize {
    (r["content"].as_array().unwrap().iter())
        .filter(|c| c["type"] == "image")
        .count()
}

#[test]
fn a_rounds_charts_and_a_trajectorys_diffusion_maps_are_drawn_returned_and_recorded() {
    let local = Path::new(env!("CARGO_MANIFEST_DIR")).join("results");
    std::fs::create_dir_all(&local).unwrap();
    let project = tempfile::tempdir_in(&local).unwrap();
    let data = project.path().join("data");
    std::fs::create_dir_all(&data).unwrap();
    let r0 = round_with_profile(&data).to_string_lossy().into_owned();
    let t = trajectory(&data).to_string_lossy().into_owned();
    let mut s = Server::start(project.path());

    let r = s.request(
        "tools/call",
        json!({ "name": "annotate_figure", "arguments": { "round": r0, "charts": ["markers"] } }),
    );
    let r = &r["result"];
    assert_ne!(r["isError"], true, "{r}");
    assert_eq!(images(r), 1, "{r}");
    let got = &r["structuredContent"];
    let files: Vec<&str> = (got["files"].as_array().unwrap().iter())
        .map(|f| f.as_str().unwrap())
        .collect();
    assert!(
        files[0].ends_with("run.figure.heatmap_markers.png"),
        "{files:?}"
    );
    assert!(files.iter().all(|f| Path::new(f).is_file()), "{files:?}");
    assert_eq!(
        got["numbers"]["markers"]["rows"].as_array().unwrap().len(),
        3,
        "{got}"
    );
    let terms = (got["offers"]["charts"].as_array().unwrap().iter())
        .find(|c| c["chart"] == "terms")
        .unwrap();
    assert!(
        terms["missing"]
            .as_str()
            .unwrap()
            .contains("go, gaf or gmt"),
        "{terms}"
    );

    // What the run offers, and a map or colouring it lacks refused before
    // anything is drawn.
    let listed = s
        .call("trajectory_figure", json!({ "run": t, "list": true }))
        .unwrap();
    let maps = &listed["maps"];
    assert_eq!(maps["diffusion_opens_on"], "dc1,dc3", "{maps}");
    assert_eq!(maps["components"].as_array().unwrap().len(), 3, "{maps}");
    assert!(
        maps["hint"].as_str().unwrap().contains("senna layout"),
        "{maps}"
    );
    for (args, says) in [
        (json!({ "run": t, "maps": ["dc0,dc1"] }), "trivial"),
        (json!({ "run": t, "maps": ["phate"] }), "dc1..dc3"),
        (
            json!({ "run": t, "maps": ["diffusion"], "colour": "cluster" }),
            "lineage",
        ),
    ] {
        let e = s.call("trajectory_figure", args).unwrap_err();
        assert!(e.contains(says), "{says} in {e}");
    }
    let drawn: Vec<_> = std::fs::read_dir(&data)
        .unwrap()
        .flatten()
        .filter(|e| e.file_name().to_string_lossy().contains("figure"))
        .collect();
    assert!(drawn.is_empty(), "nothing drawn: {drawn:?}");

    let r = s.request(
        "tools/call",
        json!({ "name": "trajectory_figure",
                "arguments": { "run": t, "maps": ["suggested", "dc2,dc3"], "colour": "lineage" } }),
    );
    let r = &r["result"];
    assert_ne!(r["isError"], true, "{r}");
    assert_eq!(images(r), 2, "one image per map: {r}");
    let got = &r["structuredContent"];
    let files: Vec<&str> = (got["files"].as_array().unwrap().iter())
        .map(|f| f.as_str().unwrap())
        .collect();
    assert!(
        files[0].ends_with("t.figure.diffusion_dc1_dc3.lineage.png")
            && files[2].ends_with("t.figure.diffusion_dc2_dc3.lineage.png"),
        "{files:?}"
    );
    assert_eq!(got["maps"]["diffusion_opens_on"], "dc1,dc3", "{got}");
    // The look: the labels' size and the page, as the child is told it.
    let small = s
        .call(
            "trajectory_figure",
            json!({ "run": t, "maps": ["dc1,dc3"], "labels": "off", "page": "paper", "width": 3.5, "dpi": 100 }),
        )
        .unwrap();
    let png = small["files"][0].as_str().unwrap();
    let (w, _) = image::image_dimensions(png).unwrap();
    assert_eq!(w, 350, "3.5 in at 100 dpi");
    let refused = s.call(
        "trajectory_figure",
        json!({ "run": t, "maps": ["dc1,dc3"], "focus": 99 }),
    );
    assert!(
        refused.unwrap_err().contains("K99"),
        "a cluster the run lacks"
    );
    s.stop();

    let results = project.path().join("results");
    let index = std::fs::read_to_string(results.join("INDEX.tsv")).unwrap();
    assert!(
        index.contains("charts and maps of run.senna.json")
            && index.contains("trajectory figures of t.senna.json on suggested, dc2,dc3"),
        "{index}"
    );
    let day = Path::new(files[0]).parent().unwrap();
    let script = std::fs::read_to_string(day.join("runall.sh")).unwrap();
    assert!(
        script.contains("lupin figure round") && script.contains("--map=dc2,dc3"),
        "{script}"
    );
    assert!(
        script.contains("--labels=off --width=3.5 --dpi=100"),
        "{script}"
    );
    assert!(
        script.contains("--page=slide"),
        "an agent's figures are for the report's slides unless asked: {script}"
    );
}
