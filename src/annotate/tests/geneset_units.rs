//! Gene sets trained as units against a run's frozen gene table: planted
//! genes in two directions, sets drawn from each.

use super::*;

/// Twelve genes in 4 dimensions: GENE0-5 along the first axis, GENE6-11
/// along the second, each with a small planted offset.
fn genes() -> GeneTable {
    let mut rows = Vec::new();
    for g in 0..12 {
        let mut r = [0.0f32; 4];
        r[usize::from(g >= 6)] = 1.0;
        r[2] = 0.05 * (g % 3) as f32;
        r[3] = -0.05 * (g % 2) as f32;
        rows.extend(r);
    }
    GeneTable {
        names: (0..12).map(|g| format!("GENE{g}")).collect(),
        rows,
        dim: 4,
    }
}

fn set(name: &str, kind: SetKind, genes: &[usize]) -> GeneSet {
    GeneSet {
        name: name.into(),
        kind,
        genes: genes.iter().map(|g| format!("GENE{g}")).collect(),
    }
}

fn cos(a: &[f32], b: &[f32]) -> f32 {
    let dot: f32 = a.iter().zip(b).map(|(x, y)| x * y).sum();
    let n = |v: &[f32]| v.iter().map(|x| x * x).sum::<f32>().sqrt();
    dot / (n(a) * n(b))
}

fn quick() -> UnitConfig {
    UnitConfig {
        epochs: 60,
        seed: 7,
        ..UnitConfig::default()
    }
}

#[test]
fn a_set_is_trained_toward_its_genes_while_the_genes_stay_fixed() {
    let g = genes();
    let sets = [
        set("CT1", SetKind::CellType, &[0, 1, 2, 3, 4]),
        set("CT2", SetKind::CellType, &[6, 7, 8, 9, 10]),
        set("TERM1", SetKind::Term, &[1, 2, 3]),
        set("TERM2", SetKind::Term, &[7, 8, 9]),
    ];
    // Training checks the gene table came back as it went in.
    let u = train_units(&g, &sets, &quick()).unwrap();
    let (a, b) = ([1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]);
    for (name, near, far) in [
        ("CT1", a, b),
        ("CT2", b, a),
        ("TERM1", a, b),
        ("TERM2", b, a),
    ] {
        let c = u.centre(name).unwrap();
        assert!(
            cos(c, &near) > cos(c, &far) + 0.5,
            "{name} at {c:?} is not toward its genes"
        );
    }
}

#[test]
fn coverage_counts_the_genes_the_table_has_and_too_few_drop_the_set() {
    let g = genes();
    let mut thin = set("CT3", SetKind::CellType, &[5]);
    thin.genes.push("GENE99".into());
    let mut part = set("CT1", SetKind::CellType, &[0, 1, 2]);
    part.genes.push("GENE98".into());
    let sets = [part, set("CT2", SetKind::CellType, &[6, 7, 8]), thin];
    let u = train_units(&g, &sets, &quick()).unwrap();
    assert_eq!(u.coverage("CT1"), Some((3, 4)));
    assert!(
        u.centre("CT3").is_none(),
        "one gene in the table is too few"
    );
    assert_eq!(u.dropped, [("CT3".to_string(), 1, 2)]);
}

/// Cells around three planted centres in 4 dimensions: CT1 tight on the
/// first axis, CT2 wide on the second, and CT3 a little off CT2.
fn planted_regions() -> (Units, Vec<f32>, usize) {
    let unit = |v: [f32; 4]| {
        let n = v.iter().map(|x| x * x).sum::<f32>().sqrt();
        v.map(|x| x / n)
    };
    let centres = [
        unit([1.0, 0.0, 0.0, 0.0]),
        unit([0.0, 1.0, 0.0, 0.0]),
        unit([0.0, 1.0, 0.6, 0.0]),
    ];
    let units = Units {
        origin: vec![0.0; 4],
        names: vec!["CT1".into(), "CT2".into(), "CT3".into()],
        kinds: vec![SetKind::CellType; 3],
        centres: centres.concat(),
        dim: 4,
        covered: vec![(3, 3); 3],
        bins: vec![3; 3],
        dropped: Vec::new(),
        decoys: Vec::new(),
    };
    // A fixed spread of offsets, scaled per type: no random draws.
    let offsets: Vec<[f32; 4]> = (0..40)
        .map(|i| {
            let t = i as f32 * 0.7;
            [t.sin(), (1.3 * t).cos(), (0.9 * t).sin(), (1.7 * t).cos()]
        })
        .collect();
    let mut cells = Vec::new();
    for (c, spread) in [(0usize, 0.05f32), (1, 0.35), (2, 0.1)] {
        for o in &offsets {
            let mut v = centres[c];
            for d in 0..4 {
                v[d] += spread * o[d];
            }
            cells.extend(v);
        }
    }
    let n = cells.len() / 4;
    (units, cells, n)
}

#[test]
fn each_type_is_fitted_a_radius_from_its_cells() {
    let (units, cells, n) = planted_regions();
    let r = fit_regions(&units, SetKind::CellType, &cells, n, 100).unwrap();
    let radius = |name: &str| r.radius[r.names.iter().position(|x| x == name).unwrap()];
    assert!(radius("CT1") < radius("CT3"), "{:?}", r.radius);
    assert!(radius("CT3") < radius("CT2"), "{:?}", r.radius);
    // The tight type's own cells are its, with near certainty.
    let k = r.names.iter().position(|x| x == "CT1").unwrap();
    assert!((0..40).all(|c| r.membership[c * r.names.len() + k] > 0.95));
}

#[test]
fn types_that_share_cells_are_connected_and_apart_ones_are_not() {
    let (units, cells, n) = planted_regions();
    let r = fit_regions(&units, SetKind::CellType, &cells, n, 100).unwrap();
    let k = r.names.len();
    let link = |a: usize, b: usize| r.links[a * k + b];
    assert!(link(1, 2) > 10.0 * link(0, 1), "{:?}", r.links);
    assert_eq!(link(1, 2), link(2, 1), "a connection has no direction");
}

#[test]
fn a_pass_calls_each_cell_by_its_region_and_writes_the_units_and_links() {
    use crate::annotate::marker_score::toy;
    let dir = tempfile::tempdir().unwrap();
    let (feat, cells, groups) = toy::matrices();
    let panel = toy::panel(dir.path());
    let out = dir.path().join("p").to_string_lossy().into_owned();
    let terms = Terms {
        sets: vec![
            GeneSet {
                name: "TERM1".into(),
                kind: SetKind::Term,
                genes: ["g0", "g1"].map(String::from).to_vec(),
            },
            GeneSet {
                name: "TERM2".into(),
                kind: SetKind::Term,
                genes: ["g3", "g4"].map(String::from).to_vec(),
            },
        ],
        pool: Vec::new(),
    };
    let cfg = UnitConfig {
        decoys: 9,
        ..quick()
    };
    let o = run_pass(&panel, &out, &feat, &cells, &groups, &terms, &cfg).unwrap();

    let labels = std::fs::read_to_string(format!("{out}.argmax.tsv")).unwrap();
    let call = |cell: &str| {
        (labels.lines())
            .find(|l| l.starts_with(cell))
            .and_then(|l| l.split('\t').nth(1))
            .unwrap()
            .to_string()
    };
    assert_eq!(call("c00"), "A");
    assert_eq!(call("c19"), "B");
    assert_eq!(
        o.argmax.as_deref(),
        Some(format!("{out}.argmax.parquet").as_str())
    );

    use legume_numeric::matrix::traits::IoOps;
    let t = Mat::from_parquet_with_row_names(&format!("{out}{UNITS_PARQUET}"), Some(0)).unwrap();
    assert_eq!(t.rows, ["A".into(), "B".into()] as [Box<str>; 2]);
    assert!(t.cols.iter().any(|c| &**c == "radius"), "{:?}", t.cols);
    let links = std::fs::read_to_string(format!("{out}{TYPE_LINKS_TSV}")).unwrap();
    assert!(links.lines().any(|l| l.starts_with("A\tB\t")), "{links}");

    // Every cluster against every term, with its decoys' p and q.
    let (strings, nums) = legume_numeric::matrix::parquet::read_table_columns(
        &format!("{out}{TERM_ALIGNMENT_PARQUET}"),
        &["cluster", "term"],
        &["score", "p", "q"],
    )
    .unwrap();
    assert_eq!(strings[0].len(), 4, "two clusters × two terms");
    assert!(
        nums[1].iter().all(|p| (0.0..=1.0).contains(p)),
        "{:?}",
        nums[1]
    );
    assert_eq!(
        o.term_alignment.as_deref(),
        Some(format!("{out}{TERM_ALIGNMENT_PARQUET}").as_str())
    );
}

#[test]
fn a_pass_without_terms_writes_no_alignment() {
    use crate::annotate::marker_score::toy;
    let dir = tempfile::tempdir().unwrap();
    let (feat, cells, groups) = toy::matrices();
    let panel = toy::panel(dir.path());
    let out = dir.path().join("p").to_string_lossy().into_owned();
    let o = run_pass(
        &panel,
        &out,
        &feat,
        &cells,
        &groups,
        &Terms::default(),
        &quick(),
    )
    .unwrap();
    assert!(o.term_alignment.is_none());
    assert!(!std::path::Path::new(&format!("{out}{TERM_ALIGNMENT_PARQUET}")).exists());
}

/// Sets trained by contrast point away from the genes' common origin, not
/// at their genes: in two dimensions with every gene in one quadrant, the
/// cells are read from that origin too, or each sits nearer the wrong set.
#[test]
fn cells_are_read_from_the_genes_origin_as_the_sets_are() {
    use crate::annotate::marker_score::toy;
    let (feat, cells, _) = toy::matrices();
    let sets = [
        GeneSet {
            name: "A".into(),
            kind: SetKind::CellType,
            genes: ["g0", "g1", "g2"].map(String::from).to_vec(),
        },
        GeneSet {
            name: "B".into(),
            kind: SetKind::CellType,
            genes: ["g3", "g4", "g5"].map(String::from).to_vec(),
        },
    ];
    let u = train_units(&table_of(&feat), &sets, &quick()).unwrap();
    let r = fit_regions(&u, SetKind::CellType, &table_of(&cells).rows, 20, 100).unwrap();
    for c in 0..20 {
        let a = r.membership[c * 2];
        assert_eq!(a > 0.5, c < 10, "cell {c}: membership of A {a}");
    }
}

/// Two hundred genes in 24 dimensions, ten along each of the first twenty
/// axes, with small planted offsets on the other four; and cells in four
/// clusters, one along each of the first four axes. Random sets of these
/// genes average out, as in a real gene table.
fn four_groups() -> (GeneTable, Vec<f32>, Vec<usize>) {
    const DIM: usize = 24;
    let offset = |i: usize, d: usize| 0.1 * (((i * 7 + d * 3) % 5) as f32 - 2.0);
    let mut rows = Vec::new();
    for g in 0..200 {
        let mut r = [0.0f32; DIM];
        r[g / 10] = 1.0;
        for (d, x) in r.iter_mut().enumerate().skip(20) {
            *x = offset(g, d);
        }
        rows.extend(r);
    }
    let mut cells = Vec::new();
    let mut groups = Vec::new();
    for c in 0..80 {
        let k = c / 20;
        let mut v = [0.0f32; DIM];
        v[k] = 2.0;
        for (d, x) in v.iter_mut().enumerate().skip(20) {
            *x = offset(c, d);
        }
        cells.extend(v);
        groups.push(k);
    }
    let table = GeneTable {
        names: (0..200).map(|g| format!("GENE{g}")).collect(),
        rows,
        dim: DIM,
    };
    (table, cells, groups)
}

#[test]
fn a_term_aligned_with_a_cluster_beats_random_sets_of_its_size_and_others_do_not() {
    let (g, cells, groups) = four_groups();
    let sets = [
        set("TERM1", SetKind::Term, &[0, 1, 2, 3, 4]),
        set("TERM2", SetKind::Term, &[10, 11, 12, 13, 14]),
    ];
    let cfg = UnitConfig {
        decoys: 99,
        ..quick()
    };
    let u = train_units(&g, &sets, &cfg).unwrap();
    assert_eq!(
        u.names,
        ["TERM1", "TERM2"],
        "decoys are kept apart from the sets"
    );
    let a = align(&u, SetKind::Term, &cells, &groups).unwrap();
    let p = |cluster: usize, term: &str| {
        let t = a.terms.iter().position(|x| x == term).unwrap();
        a.p[cluster * a.terms.len() + t]
    };
    // Read from the decoys' spread, not their count: a clear alignment
    // goes below one in a hundred, the floor of 99 decoys counted.
    assert!(p(0, "TERM1") < 0.001, "{:?}", a.p);
    assert!(p(1, "TERM2") < 0.001, "{:?}", a.p);
    assert!(p(0, "TERM2") > 0.5, "{:?}", a.p);
    assert!(p(1, "TERM1") > 0.5, "{:?}", a.p);
    assert_eq!(a.q.len(), a.p.len());
}

#[test]
fn only_terms_within_the_size_window_are_trained() {
    let dir = tempfile::tempdir().unwrap();
    let obo = dir.path().join("t.obo");
    std::fs::write(&obo, "format-version: 1.2\n\n[Term]\nid: X:0\nname: root\n").unwrap();
    // Seven terms over sixty genes, of 2, 4, 6, 8, 16, 32 and 60 genes.
    let gene_names: Vec<Box<str>> = (0..60).map(|g| format!("GENE{g}").into()).collect();
    let gmt: String = [2usize, 4, 6, 8, 16, 32, 60]
        .iter()
        .map(|&n| {
            let genes: Vec<String> = (0..n).map(|g| format!("GENE{g}")).collect();
            format!("TERM{n}\tdesc\t{}\n", genes.join("\t"))
        })
        .collect();
    let path = dir.path().join("t.gmt");
    std::fs::write(&path, gmt).unwrap();
    let args = crate::annotate::args::UnitsArgs {
        gmt: Some(path.to_string_lossy().into()),
        go_obo: Some(obo.to_string_lossy().into()),
        go_min: 4,
        go_max: 32,
        ..crate::annotate::args::UnitsArgs::new("m.tsv", "out")
    };
    let terms = args.terms(&gene_names).unwrap();
    assert!(
        !terms.pool.is_empty(),
        "the genes the terms were annotated from"
    );
    let mut names: Vec<String> = terms.sets.into_iter().map(|t| t.name).collect();
    names.sort();
    assert_eq!(names, ["TERM16", "TERM32", "TERM4", "TERM6", "TERM8"]);
}

#[test]
fn decoys_are_drawn_from_the_pool_their_kind_is_annotated_from() {
    let g = genes();
    let sets = [
        set("TERM1", SetKind::Term, &[0, 1, 2]),
        set("TERM2", SetKind::Term, &[6, 7, 8]),
    ];
    // Only the first nine genes carry annotations.
    let pool: Vec<String> = (0..9).map(|i| format!("GENE{i}")).collect();
    let cfg = UnitConfig {
        decoys: 20,
        ..quick()
    };
    let u = train_units_from(&g, &sets, &[(SetKind::Term, pool)], &cfg).unwrap();
    assert_eq!(u.decoys.len(), 20, "one size bin of twenty");
    assert!(
        u.decoys.iter().flat_map(|d| &d.genes).all(|&r| r < 9),
        "a decoy drew a gene outside its pool"
    );
}

/// Three dimensions: the first a nuisance axis every gene and cell spreads
/// widely along, the other two telling types A and B apart by a little.
fn nuisance_axis() -> (MatWithNames<Mat>, MatWithNames<Mat>, Vec<usize>) {
    let wide = |i: usize| 5.0 * (((i * 7) % 11) as f32 / 5.0 - 1.0);
    let mut feat = Mat::zeros(20, 3);
    for g in 0..20 {
        feat[(g, 0)] = 3.0 + wide(g);
        feat[(g, 1 + g / 10)] = 0.4;
    }
    let mut cells = Mat::zeros(60, 3);
    for c in 0..60 {
        cells[(c, 0)] = wide(c + 3);
        cells[(c, 1 + c / 30)] = 0.3;
    }
    let names =
        |p: &str, n: usize| -> Vec<Box<str>> { (0..n).map(|i| format!("{p}{i}").into()).collect() };
    let h: Vec<Box<str>> = names("h", 3);
    (
        MatWithNames {
            mat: feat,
            rows: names("GENE", 20),
            cols: h.clone(),
        },
        MatWithNames {
            mat: cells,
            rows: names("c", 60),
            cols: h,
        },
        (0..60).map(|c| c / 30).collect(),
    )
}

#[test]
fn a_wide_axis_every_type_shares_does_not_decide_the_calls() {
    let dir = tempfile::tempdir().unwrap();
    let (feat, cells, groups) = nuisance_axis();
    let panel = dir.path().join("panel.tsv");
    let rows: String = (0..20)
        .map(|g| format!("GENE{g}\t{}\n", if g < 10 { "A" } else { "B" }))
        .collect();
    std::fs::write(&panel, rows).unwrap();
    let out = dir.path().join("w").to_string_lossy().into_owned();
    run_pass(
        &panel.to_string_lossy(),
        &out,
        &feat,
        &cells,
        &groups,
        &Terms::default(),
        &quick(),
    )
    .unwrap();
    let labels = std::fs::read_to_string(format!("{out}.argmax.tsv")).unwrap();
    let wrong: Vec<&str> = (labels.lines())
        .filter_map(|l| {
            let mut f = l.split('\t');
            let (cell, label) = (f.next()?, f.next()?);
            let i: usize = cell.strip_prefix('c')?.parse().ok()?;
            (label != if i < 30 { "A" } else { "B" }).then_some(cell)
        })
        .collect();
    assert!(wrong.len() <= 3, "{} miscalled: {wrong:?}", wrong.len());
}
