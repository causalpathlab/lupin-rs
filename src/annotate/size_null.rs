//! **The size null**: a cell type's statistic read against what panels of its
//! size score by chance in the same round.
//!
//! # The bias
//!
//! A panel's genes are not independent: a large panel's markers co-vary (a
//! neuronal programme in Visium spots, say), and a sum over many weakly lifted
//! genes builds up evidence that a few strongly lifted ones cannot. The
//! enrichment's random sets match a panel's size and abundance (GOseq's fix,
//! `enrichment::gene_strata`) but their genes are independent, so the gap
//! grows with size. A large panel then wins clusters it has no business in,
//! and a small one loses those it should win.
//!
//! # The fix
//!
//! GOseq's weighting function, with the panel's size as the covariate and the
//! round's own clusters as the data. Most cluster × type pairs are not the
//! type, so a type's statistic over the clusters shows its chance level. Per
//! type, its spread σ_c from the lower part of that distribution (where
//! clusters that are not the type dominate, so a type present in up to three
//! quarters of the clusters still shows its chance level).
//!
//! Markers that move together widen a set's null as CAMERA's variance
//! inflation does, σ² ∝ 1 + ρ(m − 1): flat over small panels, steep over
//! large ones, so not a line in log m. Across types, one curve
//! `σ²(m) = a + b·m` by Theil–Sen on the effective sizes (Kish, over the IDF
//! weights), kept only where Kendall's τ says the spread rises with size. One curve over every type lets the mid-sized
//! panels explain the largest, where a type's own distribution could not
//! tell size from a type truly in most clusters.
//!
//! Each statistic is then divided by its size's inflation,
//! `s* = s / max(1, √(σ²(m_c) / σ²(m_ref)))`, with `m_ref` the panel's
//! median size: a panel above the reference is deflated, and one below it
//! keeps its own test's statistic rather than being boosted past it. A flat
//! curve is the identity.

use statrs::distribution::{ContinuousCDF, Normal};

/// Types a fit needs.
const MIN_TYPES: usize = 4;
/// Clusters (values of one type's statistic) a type needs to enter the fit.
const MIN_VALUES: usize = 5;
/// A type's spread is floored here, so a type whose chance values tie still
/// has one.
const MIN_SPREAD: f64 = 1e-2;
/// The lower part the null is read from, and its standard-normal quantiles.
const LOW: (f64, f64) = (0.10, 0.40);
const Z_LOW: (f64, f64) = (-1.281_551_57, -0.253_347_10);
/// One-sided level of the Kendall τ that keeps the curve's slope.
const TREND_ALPHA: f64 = 0.05;

/// Whether to read the panel calls against the size null (`--size-null`).
#[derive(
    Copy,
    Clone,
    Debug,
    Default,
    PartialEq,
    Eq,
    clap::ValueEnum,
    serde::Serialize,
    serde::Deserialize,
)]
#[serde(rename_all = "lowercase")]
pub enum SizeNull {
    /// The size trend fitted on the round's clusters, the default.
    Trend,
    /// The raw statistics. A record from before the size null reads as this.
    #[default]
    Off,
}

impl SizeNull {
    /// As the command line spells it.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Trend => "trend",
            Self::Off => "off",
        }
    }
}

/// A line `a + b·x`, with Kendall's τ of the points it was fitted on.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Line {
    pub a: f64,
    pub b: f64,
    /// Kendall's τ between x and y, and its one-sided p.
    pub tau: f64,
    pub tau_p: f64,
}

impl Line {
    fn at(&self, x: f64) -> f64 {
        self.a + self.b * x
    }
}

/// One type's input to the fit: its effective size and its statistic over the
/// clusters that count (non-finite values are skipped).
pub struct TypeValues {
    pub size: f64,
    pub values: Vec<f64>,
}

/// One type's chance level, as read and as the fit puts it.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TypeNull {
    /// The type, once [`SizeTrend::named`].
    #[serde(default)]
    pub name: String,
    pub size: f64,
    /// Its own chance spread, where it had enough values.
    pub spread: Option<f64>,
    /// What [`SizeTrend::calibrate`] divides by: σ(m_c)/σ(m_ref), at least 1.
    pub scale: f64,
}

/// The fitted trend, or why there is none.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SizeTrend {
    /// `fitted`, `no trend` (a flat curve), or why no fit was made.
    pub status: String,
    /// m_ref.
    pub reference: f64,
    /// σ²(m) = a + b·m.
    pub variance: Option<Line>,
    pub types: Vec<TypeNull>,
}

impl SizeTrend {
    /// No fit: every type read as it is.
    #[must_use]
    pub fn none(status: &str, sizes: &[f64]) -> Self {
        Self {
            status: status.into(),
            reference: 0.0,
            variance: None,
            types: sizes
                .iter()
                .map(|&size| TypeNull {
                    name: String::new(),
                    size,
                    spread: None,
                    scale: 1.0,
                })
                .collect(),
        }
    }

    /// The types named `names`, in order.
    #[must_use]
    pub fn named(mut self, names: &[impl AsRef<str>]) -> Self {
        for (t, n) in self.types.iter_mut().zip(names) {
            t.name = n.as_ref().to_string();
        }
        self
    }

    /// Whether any type is moved.
    #[must_use]
    pub fn moves(&self) -> bool {
        self.types.iter().any(|t| t.scale != 1.0)
    }

    /// Say in the log what the fit found for `method`.
    pub fn report(&self, method: &str) {
        // Types the test left out have size 0.
        let largest = (self.types.iter()).max_by(|a, b| a.size.total_cmp(&b.size));
        match (self.variance, largest) {
            (Some(v), Some(t)) if self.moves() => log::info!(
                "{method} size null: {}; chance variance {:.3} + {:.4}·m (τ {:.2}, p {:.2e}), \
                 so the largest panel ({:.0} effective markers) is read at 1/{:.2} of its z",
                self.status,
                v.a,
                v.b,
                v.tau,
                v.tau_p,
                t.size,
                t.scale
            ),
            _ => log::info!("{method} size null: {}", self.status),
        }
    }

    /// Type c's statistic `s` with its size's inflation taken out.
    #[must_use]
    pub fn calibrate(&self, c: usize, s: f64) -> f64 {
        s / self.types[c].scale
    }

    /// The Spearman correlation between the types' log sizes and `y` (one
    /// value per type; non-finite left out): how far a statistic follows
    /// size. `None` with fewer than three types.
    #[must_use]
    pub fn size_correlation(&self, y: &[f64]) -> Option<f64> {
        let (x, y): (Vec<f32>, Vec<f32>) = (self.types.iter().zip(y))
            .filter(|(t, y)| t.size > 0.0 && y.is_finite())
            .map(|(t, &y)| (t.size.ln() as f32, y as f32))
            .unzip();
        let rho = legume_numeric::matrix::agreement::spearman(&x, &y);
        (x.len() >= 3 && rho.is_finite()).then_some(f64::from(rho))
    }
}

/// Fit the size trend over `types`.
#[must_use]
pub fn fit(types: &[TypeValues]) -> SizeTrend {
    let sizes: Vec<f64> = types.iter().map(|t| t.size).collect();
    let spread: Vec<Option<f64>> = types.iter().map(chance_spread).collect();
    let (m, var): (Vec<f64>, Vec<f64>) = (types.iter().zip(&spread))
        .filter_map(|(t, s)| s.filter(|_| t.size > 0.0).map(|s| (t.size, s * s)))
        .unzip();
    if m.len() < MIN_TYPES {
        return SizeTrend::none(
            &format!(
                "too few types to fit: {} with {MIN_VALUES}+ clusters, {MIN_TYPES} needed",
                m.len()
            ),
            &sizes,
        );
    }
    let variance = theil_sen(&m, &var);
    let reference = median(&m);
    let at_ref = variance.at(reference);
    let status = if variance.b == 0.0 {
        "no trend"
    } else {
        "fitted"
    };
    let out_types = (types.iter().zip(&spread))
        .map(|(t, s)| TypeNull {
            name: String::new(),
            size: t.size,
            spread: *s,
            // Never below 1: the null only deflates the panels above the
            // reference, and a small panel keeps its own test's p.
            scale: if variance.b == 0.0 || t.size <= 0.0 {
                1.0
            } else {
                (variance.at(t.size) / at_ref).sqrt().max(1.0)
            },
        })
        .collect();
    SizeTrend {
        status: status.into(),
        reference,
        variance: Some(variance),
        types: out_types,
    }
}

/// A panel's effective size over its genes' weights `w`: (Σw)²/Σw², the
/// number of equally weighted genes that carry as much.
#[must_use]
pub fn effective_size(w: impl IntoIterator<Item = f64>) -> f64 {
    let (s, s2) = w
        .into_iter()
        .filter(|w| *w > 0.0)
        .fold((0.0, 0.0), |(s, s2), w| (s + w, s2 + w * w));
    if s2 > 0.0 {
        s * s / s2
    } else {
        0.0
    }
}

/// The upper-tail p of a standard normal at `z`.
fn upper_p(z: f64) -> f64 {
    Normal::standard().sf(z)
}

/// A type's chance spread σ from the lower part of its values, as if they
/// were normal there: from the 10th to the 40th percentile.
fn chance_spread(t: &TypeValues) -> Option<f64> {
    let mut v: Vec<f64> = t.values.iter().copied().filter(|v| v.is_finite()).collect();
    if v.len() < MIN_VALUES {
        return None;
    }
    v.sort_by(f64::total_cmp);
    let (lo, hi) = (quantile(&v, LOW.0), quantile(&v, LOW.1));
    Some(((hi - lo) / (Z_LOW.1 - Z_LOW.0)).max(MIN_SPREAD))
}

/// The `p` quantile of sorted `v`, by linear interpolation.
fn quantile(v: &[f64], p: f64) -> f64 {
    let h = p * (v.len() - 1) as f64;
    let (lo, hi) = (h.floor() as usize, h.ceil() as usize);
    v[lo] + (h - lo as f64) * (v[hi] - v[lo])
}

fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(f64::total_cmp);
    quantile(&s, 0.5)
}

/// Theil–Sen: the median of the pairwise slopes, taken as 0 unless Kendall's
/// τ is positive at one-sided [`TREND_ALPHA`] (a panel's chance spread does
/// not shrink with its size), and the median residual as the intercept,
/// floored at [`MIN_SPREAD`]².
fn theil_sen(x: &[f64], y: &[f64]) -> Line {
    let n = x.len();
    let mut slopes = Vec::with_capacity(n * (n - 1) / 2);
    let (mut concordant, mut discordant) = (0i64, 0i64);
    for i in 0..n {
        for j in i + 1..n {
            let (dx, dy) = (x[j] - x[i], y[j] - y[i]);
            if dx != 0.0 {
                slopes.push(dy / dx);
            }
            match (dx * dy).partial_cmp(&0.0) {
                Some(std::cmp::Ordering::Greater) => concordant += 1,
                Some(std::cmp::Ordering::Less) => discordant += 1,
                _ => {}
            }
        }
    }
    let tau = (concordant - discordant) as f64 / (n * (n - 1) / 2) as f64;
    // The normal approximation to τ's null, var = 2(2n+5)/(9n(n−1)).
    let nf = n as f64;
    let tau_p = upper_p(tau / (2.0 * (2.0 * nf + 5.0) / (9.0 * nf * (nf - 1.0))).sqrt());
    let b = if slopes.is_empty() || tau_p >= TREND_ALPHA {
        0.0
    } else {
        median(&slopes).max(0.0)
    };
    let a = median(&x.iter().zip(y).map(|(x, y)| y - b * x).collect::<Vec<_>>())
        .max(MIN_SPREAD * MIN_SPREAD);
    Line { a, b, tau, tau_p }
}

#[cfg(test)]
#[path = "tests/size_null.rs"]
mod tests;
