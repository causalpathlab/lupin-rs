//! `{out}.trajectory_edges.parquet`: one row per unordered pair of node types
//! with its connectivity and, for a prior edge or a candidate, the data's
//! verdict. The run writes it; the order view and the figures read it.

use anyhow::{Context, Result};
use legume_numeric::matrix::parquet::{read_table_columns, write_named_table, Column};

/// What the connectivity check says about a pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Verdict {
    /// A prior edge whose connectivity reaches the threshold.
    Supported,
    /// A prior edge below it.
    Unsupported,
    /// A pair the prior does not order whose connectivity reaches it.
    Candidate,
    /// An edge the run added: a node type the prior leaves without an edge
    /// joins it along its strongest connectivity.
    Inferred,
}

impl Verdict {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Supported => "supported",
            Self::Unsupported => "unsupported",
            Self::Candidate => "candidate",
            Self::Inferred => "inferred",
        }
    }

    /// The verdict `s` names; `None` for an unremarkable pair (written as "").
    pub(crate) fn parse(s: &str) -> Option<Self> {
        match s {
            "supported" => Some(Self::Supported),
            "unsupported" => Some(Self::Unsupported),
            "candidate" => Some(Self::Candidate),
            "inferred" => Some(Self::Inferred),
            _ => None,
        }
    }
}

/// One pair. A prior or inferred edge runs `a → b`; any other pair is
/// written `a < b`.
#[derive(Debug, Clone)]
pub(crate) struct EdgeRow {
    pub(crate) a: Box<str>,
    pub(crate) b: Box<str>,
    pub(crate) connectivity: f32,
    pub(crate) in_prior: bool,
    pub(crate) verdict: Option<Verdict>,
    /// Fraction of `b`'s cells beyond `a`'s median pseudotime; NaN unless the
    /// pair is a prior edge and pseudotime was computed.
    pub(crate) order_agreement: f32,
}

impl EdgeRow {
    /// `true` when the row is about `x` and `y`, either way round.
    pub(crate) fn is_pair(&self, x: &str, y: &str) -> bool {
        (self.a.as_ref() == x && self.b.as_ref() == y)
            || (self.a.as_ref() == y && self.b.as_ref() == x)
    }

    /// The verdict and numbers in a few words, for a table cell.
    pub(crate) fn summary(&self) -> String {
        let verdict = self.verdict.map_or("", Verdict::as_str);
        if self.order_agreement.is_finite() {
            format!(
                "{verdict} {:.2} · order {:.2}",
                self.connectivity, self.order_agreement
            )
        } else {
            format!("{verdict} {:.2}", self.connectivity)
        }
    }
}

pub(crate) fn write(path: &str, rows: &[EdgeRow]) -> Result<()> {
    let col = |f: fn(&EdgeRow) -> Box<str>| rows.iter().map(f).collect::<Vec<_>>();
    let a = col(|r| r.a.clone());
    let b = col(|r| r.b.clone());
    let in_prior = col(|r| if r.in_prior { "true" } else { "false" }.into());
    let verdict = col(|r| r.verdict.map_or("", Verdict::as_str).into());
    let conn: Vec<f32> = rows.iter().map(|r| r.connectivity).collect();
    let agree: Vec<f32> = rows.iter().map(|r| r.order_agreement).collect();
    write_named_table(
        path,
        "a",
        &a,
        &[
            ("b".into(), Column::Str(&b)),
            ("connectivity".into(), Column::F32(&conn)),
            ("in_prior".into(), Column::Str(&in_prior)),
            ("verdict".into(), Column::Str(&verdict)),
            ("order_agreement".into(), Column::F32(&agree)),
        ],
    )
}

pub(crate) fn read(path: &str) -> Result<Vec<EdgeRow>> {
    let (s, n) = read_table_columns(
        path,
        &["a", "b", "in_prior", "verdict"],
        &["connectivity", "order_agreement"],
    )
    .with_context(|| format!("reading {path}"))?;
    Ok((0..s[0].len())
        .map(|i| EdgeRow {
            a: s[0][i].clone(),
            b: s[1][i].clone(),
            connectivity: n[0][i] as f32,
            in_prior: s[2][i].as_ref() == "true",
            verdict: Verdict::parse(&s[3][i]),
            order_agreement: n[1][i] as f32,
        })
        .collect())
}
