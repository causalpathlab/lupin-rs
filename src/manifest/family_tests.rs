use super::*;
use std::fs;

/// A run `X` with a round `X.L0`, a round made from it `X.L0.L1`, a
/// trajectory `X.T1` on `X.L0`, one under another prefix (`other.T9`) on
/// it too, and an unrelated run `Y`; each round loadable as a round.
pub(crate) fn family_dir(dir: &Path) {
    let cells: Vec<Box<str>> = ["c1", "c2", "c3", "c4"]
        .iter()
        .map(|c| (*c).into())
        .collect();
    crate::manifest::rounds::write_clusters(
        &dir.join("X.clusters.parquet").to_string_lossy(),
        &cells,
        &[Some(0), Some(0), Some(1), Some(1)],
    )
    .unwrap();
    fs::write(
        dir.join("X.L0.argmax.tsv"),
        "cell\tcell_type\tprobability\nc1\tCT1\t0.9\nc2\tCT1\t0.9\nc3\tCT2\t0.9\nc4\tCT2\t0.9\n",
    )
    .unwrap();
    let manifest = |name: &str, extra: &str| {
        fs::write(
            dir.join(format!("{name}.senna.json")),
            format!(
                r#"{{"version":2,"kind":"topic","prefix":"{name}",
                    "cluster":{{"clusters":"X.clusters.parquet"}}{extra}}}"#
            ),
        )
        .unwrap();
    };
    let annotated =
        |source: &str| format!(r#","annotate":{{"argmax":"X.L0.argmax.tsv","source":"{source}"}}"#);
    manifest("X", "");
    manifest("X.L0", &annotated("X.senna.json"));
    manifest("X.L0.L1", &annotated("X.L0.senna.json"));
    let trajectory = r#","trajectory":{"prior":"t.tsv"}"#;
    manifest(
        "X.T1",
        &format!("{}{trajectory}", annotated("X.L0.senna.json")),
    );
    manifest(
        "other.T9",
        &format!("{}{trajectory}", annotated("X.L0.senna.json")),
    );
    manifest("Y", "");
}

#[test]
fn the_family_is_the_run_its_rounds_and_its_trajectories() {
    let tmp = tempfile::tempdir().unwrap();
    family_dir(tmp.path());
    // From any member, the same family.
    for opened in ["X", "X.L0.L1", "X.T1"] {
        let got = family(&tmp.path().join(format!("{opened}.senna.json")));
        let rows: Vec<(String, Kind)> = got.iter().map(|m| (m.name(), m.pick.kind)).collect();
        assert_eq!(
            rows,
            vec![
                ("X".into(), Kind::Run),
                ("X.L0".into(), Kind::Round),
                ("X.L0.L1".into(), Kind::Round),
                ("X.T1".into(), Kind::Trajectory),
                ("other.T9".into(), Kind::Trajectory),
            ],
            "opened {opened}"
        );
        assert!(has_round(&tmp.path().join(format!("{opened}.senna.json"))));
    }
    // A run whose family has no round says so.
    assert!(!has_round(&tmp.path().join("Y.senna.json")));
}

#[test]
fn a_pinto_run_lists_with_its_lupin_rounds() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    fs::write(dir.join("P.pinto.json"), "{}").unwrap();
    fs::write(
        dir.join("P.L0.lupin.json"),
        r#"{"version":2,"kind":"topic","prefix":"P.L0",
            "annotate":{"argmax":"P.L0.argmax.tsv","source":"P.pinto.json"}}"#,
    )
    .unwrap();
    let got = family(&dir.join("P.L0.lupin.json"));
    let names: Vec<String> = got.iter().map(Member::name).collect();
    // The pinto run is listed when it loads as a run; its round always is.
    assert!(names.contains(&"P.L0".to_string()), "{names:?}");
    if let Some(run) = got.iter().find(|m| m.name() == "P") {
        assert_eq!(run.pick.kind, Kind::Run);
    }
}

#[test]
fn the_tags_format_and_parse_the_same_names() {
    for tag in [Tag::Round, Tag::Trajectory, Tag::Chain] {
        let name = tag.name("X.L0", 12);
        assert_eq!(tag.parse(&name), Some(("X.L0", 12)), "{name}");
    }
    assert_eq!(Tag::Chain.parse("a.rb"), None);
    assert_eq!(Tag::Chain.parse("a.r"), None);
    let tmp = tempfile::tempdir().unwrap();
    let run = tmp.path().join("X.senna.json");
    fs::write(&run, "{}").unwrap();
    fs::write(tmp.path().join("X.L1.senna.json"), "{}").unwrap();
    let (from, out) = pass_origin(&run);
    assert_eq!(from, run);
    assert!(out.ends_with("X.L2"), "the next free round, not L1: {out}");
}

#[test]
fn the_stem_drops_round_and_trajectory_tags_only() {
    for (name, want) in [
        ("X.senna.json", "X"),
        ("X.L0.L1.senna.json", "X"),
        ("X.T2.senna.json", "X"),
        ("X.r3.senna.json", "X"),
        ("run.v2.senna.json", "run.v2"),
        ("X.Late.senna.json", "X.Late"),
    ] {
        assert_eq!(stem(Path::new(name)), want, "{name}");
    }
}

#[test]
fn the_one_on_screen_is_found_by_file_whatever_the_spelling() {
    let tmp = tempfile::tempdir().unwrap();
    family_dir(tmp.path());
    let members = family(&tmp.path().join("X.senna.json"));
    // Another spelling of the same file: through `..`.
    let sub = tmp.path().join("sub");
    fs::create_dir(&sub).unwrap();
    let other = sub.join("..").join("X.L0.senna.json");
    assert_eq!(current(&members, &other), Some(1));
    assert!(members[1].pick.labelled);
}

#[test]
fn a_pass_still_staged_is_no_member_of_the_family() {
    let tmp = tempfile::tempdir().unwrap();
    family_dir(tmp.path());
    // What a pass killed before it was promoted leaves.
    fs::copy(
        tmp.path().join("X.L0.senna.json"),
        crate::manifest::staging::staging_prefix(&tmp.path().join("X.L0").to_string_lossy())
            + ".senna.json",
    )
    .unwrap();
    let got = family(&tmp.path().join("X.senna.json"));
    assert!(got.iter().all(|m| !m.name().contains("staging")));
}
