//! Which scatters a trajectory's figures are drawn on, by name: a layout's
//! method (`phate`, `umap`, …), `best`, `all`, `suggested` or `diffusion`
//! (the diffusion map at the pair it opens on), or two diffusion components
//! as `dc1,dc3`; and what each diffusion component shows, so a pair can be
//! chosen from the numbers: its spread, how closely it follows pseudotime
//! (Spearman's ρ) and how well it parts the lineages (η²).
//!
//! The diffusion map opens on the suggested pair: of the components the run
//! kept (spread over at least `dc_min_share` of the cells), the one that
//! follows pseudotime best on x, and on y the one that best parts the
//! lineages, else the next best by ρ.

use super::look::Look;
use super::{CellFigures, Colouring, Panel, Style};
use crate::trajectory::diffusion::participation;
use anyhow::{Context, Result};
use legume_numeric::matrix::agreement::spearman;
use legume_numeric::matrix::dense_mat_io::Mat;
use serde::Serialize;
use std::path::PathBuf;

/// What `lupin figure` draws: its maps, in the style's colouring.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Drawing {
    pub(crate) maps: Maps,
    pub(crate) style: Style,
    /// Its page and the look the style was set from.
    pub(crate) look: Look,
}

/// Which maps a drawing has.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Maps {
    /// A trajectory's figures as they open: the best scatter as
    /// `{base}.map`, the order and the connectivity.
    Opening,
    /// These scatters, each as `{base}.{slug}.{colouring}`.
    These(Vec<Panel>),
}

/// The maps and colourings a trajectory's figures offer, for an agent to
/// choose from.
#[derive(Debug, Serialize)]
pub(crate) struct MapsOffered {
    /// The layouts by name, best first.
    layouts: Vec<String>,
    /// The pair `diffusion` and `suggested` draw, as `dcX,dcY`.
    diffusion_opens_on: Option<String>,
    /// Each diffusion component past the trivial DC0.
    components: Vec<ComponentFacts>,
    colourings: Vec<String>,
    /// What to run for a layout the run lacks.
    hint: Option<String>,
}

/// What a diffusion component shows.
#[derive(Clone, Debug, Serialize)]
pub(crate) struct ComponentFacts {
    dc: usize,
    /// Its standard deviation over the cells.
    sd: f64,
    /// The share of the cells it is spread over (its participation ratio).
    share: f64,
    /// Whether the run kept it in pseudotime: its share is at least the
    /// run's `dc_min_share`.
    kept: bool,
    /// |Spearman's ρ| with pseudotime over the cells it reaches.
    rho: Option<f64>,
    /// The share of its variance between the lineages (η²), over the cells
    /// on one lineage; none without two lineages.
    lineage_eta2: Option<f64>,
}

impl Colouring {
    /// Its name in file names and on the command line.
    pub(crate) fn slug(self) -> String {
        use clap::ValueEnum;
        (self.to_possible_value())
            .expect("no variant is skipped")
            .get_name()
            .to_string()
    }
}

impl CellFigures {
    /// A scatter's name as [`Self::map_panels`] reads it back, so a command
    /// line can name the pair `suggested` stood for.
    pub(crate) fn map_name(&self, panel: Panel) -> String {
        match panel {
            Panel::Layout { k } => (self.layouts.get(k))
                .and_then(|l| l.method.clone())
                .unwrap_or_else(|| "layout".into()),
            Panel::Diffusion { x, y } => format!("dc{x},dc{y}"),
            _ => panel.slug(),
        }
    }

    /// The scatters `spec` names: a layout's method, `best` (the first
    /// scatter), `all` (every one the screen steps through), `diffusion` or
    /// `suggested` (the diffusion map at the pair it opens on), or two
    /// diffusion components (`dc1,dc3`; also `dc1 x dc3`, `dc1_dc3`).
    pub(crate) fn map_panels(&self, spec: &str) -> Result<Vec<Panel>> {
        let s = spec.trim().to_lowercase();
        let layout = (0..self.layouts.len()).find(|&k| self.map_name(Panel::Layout { k }) == s);
        if let Some(k) = layout {
            return Ok(vec![Panel::Layout { k }]);
        }
        match s.as_str() {
            "all" => {
                let all = self.scatters();
                anyhow::ensure!(!all.is_empty(), "{}", self.no_map());
                Ok(all)
            }
            "best" => Ok(vec![self.best_scatter().with_context(|| self.no_map())?]),
            "diffusion" | "suggested" => Ok(vec![self.opening_pair()?]),
            _ if s.starts_with("dc") => Ok(vec![self.diffusion_named(&s)?]),
            _ => anyhow::bail!("no map {spec:?} in this run: {}", self.map_names()),
        }
    }

    /// Why the run has no map at all.
    fn no_map(&self) -> String {
        let mut said = "the run has no layout and no diffusion map".to_string();
        if let Some(h) = &self.layout_hint {
            said += &format!("; {h}");
        }
        said
    }

    /// The maps the run offers, in words.
    fn map_names(&self) -> String {
        let mut names = vec!["best".to_string(), "all".into()];
        names.extend((0..self.layouts.len()).map(|k| self.map_name(Panel::Layout { k })));
        let mut said = names.join(", ");
        if let Some(last) = self.last_component() {
            said += &format!(", diffusion, suggested, or two of dc1..dc{last} as dcX,dcY");
        }
        if let Some(h) = &self.layout_hint {
            said += &format!(" ({h})");
        }
        said
    }

    /// The diffusion map, with two components past the trivial one.
    fn diffusion_map(&self) -> Result<&Mat> {
        let d = self
            .diffusion
            .as_ref()
            .context("the run has no diffusion map")?;
        anyhow::ensure!(
            d.ncols() >= 3,
            "the run's diffusion map has no two components past the trivial DC0"
        );
        Ok(d)
    }

    /// The last diffusion component, when there are two past the trivial one.
    fn last_component(&self) -> Option<usize> {
        self.diffusion_map().ok().map(|d| d.ncols() - 1)
    }

    /// The pair `s` names, as `dcX,dcY`.
    fn diffusion_named(&self, s: &str) -> Result<Panel> {
        let last = self.diffusion_map()?.ncols() - 1;
        let parts: Vec<&str> = s
            .split([',', '_', ' ', 'x', '×'])
            .filter(|p| !p.is_empty())
            .collect();
        let picked: Vec<usize> = (parts.iter())
            .map(|p| {
                (p.strip_prefix("dc"))
                    .and_then(|n| n.parse().ok())
                    .with_context(|| format!("{p:?} is not a diffusion component (dc1..dc{last})"))
            })
            .collect::<Result<_>>()?;
        let [x, y] = picked[..] else {
            anyhow::bail!("a diffusion map is two components, as dc1,dc2: not {s:?}");
        };
        anyhow::ensure!(
            x != 0 && y != 0,
            "DC0 is the trivial component, constant over the cells: pick from dc1..dc{last}"
        );
        anyhow::ensure!(
            x <= last && y <= last,
            "the run keeps dc1..dc{last}: not {s:?}"
        );
        anyhow::ensure!(x != y, "a diffusion map is two different components");
        Ok(Panel::Diffusion { x, y })
    }

    /// The pair the diffusion map opens on ([`suggested_pair`]).
    pub(super) fn opening_pair(&self) -> Result<Panel> {
        self.diffusion_map()?;
        let (x, y) = suggested_pair(self.component_facts()).unwrap_or((1, 2));
        Ok(Panel::Diffusion { x, y })
    }

    /// Each diffusion component past the trivial one, with what it shows;
    /// made once.
    fn component_facts(&self) -> &[ComponentFacts] {
        self.dc_facts.get_or_init(|| {
            let Ok(d) = self.diffusion_map() else {
                return Vec::new();
            };
            (1..d.ncols())
                .map(|j| {
                    let v: Vec<f32> = d.column(j).iter().copied().collect();
                    let wide: Vec<f64> = v.iter().map(|&v| f64::from(v)).collect();
                    let share = participation(wide.iter().copied());
                    ComponentFacts {
                        dc: j,
                        sd: sd(&wide),
                        share,
                        kept: share >= f64::from(self.dc_min_share),
                        rho: self.rho_with_pseudotime(&v),
                        lineage_eta2: eta2(&wide, &self.lineage),
                    }
                })
                .collect()
        })
    }

    /// |Spearman's ρ| of `v` with pseudotime over the cells both reach.
    fn rho_with_pseudotime(&self, v: &[f32]) -> Option<f64> {
        let (a, b): (Vec<f32>, Vec<f32>) = (self.pseudotime.iter().zip(v))
            .filter(|(t, v)| t.is_finite() && v.is_finite())
            .unzip();
        let rho = spearman(&a, &b);
        rho.is_finite().then(|| f64::from(rho).abs())
    }

    /// The colouring `name` names (`type`, or as the screen names it, `cell
    /// type`), when the run has it.
    pub(crate) fn colouring_named(&self, name: &str) -> Result<Colouring> {
        use clap::ValueEnum;
        let have = self.colourings();
        let listed = || have.iter().map(|c| c.slug()).collect::<Vec<_>>().join(", ");
        let c = Colouring::from_str(name.trim(), true)
            .map_err(|_| anyhow::anyhow!("no colouring {name:?}: the run's are {}", listed()))?;
        anyhow::ensure!(
            have.contains(&c),
            "this run has no {} colouring: its colourings are {}",
            c.slug(),
            listed()
        );
        Ok(c)
    }

    /// What a trajectory's figures draw: the maps `maps` names, in `colour`
    /// (the opening colouring when none); no maps, the opening set.
    pub(crate) fn drawing<S: AsRef<str>>(
        &self,
        maps: &[S],
        colour: Option<&str>,
        look: &Look,
    ) -> Result<Drawing> {
        Ok(Drawing {
            maps: if maps.is_empty() {
                Maps::Opening
            } else {
                Maps::These(self.maps_named(maps)?)
            },
            style: self.style_in(colour, look)?,
            look: look.clone(),
        })
    }

    /// The scatters `maps` name, each once.
    pub(crate) fn maps_named<S: AsRef<str>>(&self, maps: &[S]) -> Result<Vec<Panel>> {
        let mut panels = Vec::new();
        for m in maps {
            panels.extend(self.map_panels(m.as_ref())?);
        }
        Ok(once_each(panels))
    }

    /// The opening style in `colour`, when one is named, as `look` asks.
    pub(crate) fn style_in(&self, colour: Option<&str>, look: &Look) -> Result<Style> {
        let mut style = Style::opening(self);
        if let Some(c) = colour {
            style.colouring = self.colouring_named(c)?;
        }
        look.apply(&mut style, self)?;
        Ok(style)
    }

    /// Each panel `d` draws under `base` and its two files: the opening
    /// set's ([`Self::panel_files`]), or each map as
    /// `{base}.{slug}.{colouring}` `.png` and `.pdf`.
    pub(super) fn drawn_files(&self, d: &Drawing, base: &str) -> Vec<(Panel, [PathBuf; 2])> {
        match &d.maps {
            Maps::Opening => (self.panels().into_iter())
                .map(|p| (p, Self::panel_files(p, base)))
                .collect(),
            Maps::These(maps) => (maps.iter())
                .map(|&p| {
                    let b = format!("{base}.{}.{}", self.slug(p), d.style.colouring.slug());
                    (p, ["png", "pdf"].map(|e| PathBuf::from(format!("{b}.{e}"))))
                })
                .collect(),
        }
    }

    /// The maps `d` draws as a command line names them, so the pair
    /// `suggested` stood for is named; none for the opening set.
    pub(crate) fn map_specs(&self, d: &Drawing) -> Vec<String> {
        match &d.maps {
            Maps::Opening => Vec::new(),
            Maps::These(maps) => maps.iter().map(|&p| self.map_name(p)).collect(),
        }
    }

    /// The maps and colourings the run offers.
    pub(crate) fn maps_offered(&self) -> MapsOffered {
        MapsOffered {
            layouts: (0..self.layouts.len())
                .map(|k| self.map_name(Panel::Layout { k }))
                .collect(),
            diffusion_opens_on: (self.opening_pair().ok()).map(|p| self.map_name(p)),
            components: self.component_facts().to_vec(),
            colourings: self.colourings().into_iter().map(Colouring::slug).collect(),
            hint: self.layout_hint.clone(),
        }
    }
}

/// `items` in order, each kept once.
pub(super) fn once_each<T: PartialEq>(items: Vec<T>) -> Vec<T> {
    let mut out: Vec<T> = Vec::with_capacity(items.len());
    for t in items {
        if !out.contains(&t) {
            out.push(t);
        }
    }
    out
}

/// The pair to show of `facts`' kept components: on x the one that follows
/// pseudotime best (largest |ρ|), on y the one that best parts the lineages
/// (largest η²), else the next best by |ρ|. Without pseudotime, the first
/// two kept. None with fewer than two kept.
fn suggested_pair(facts: &[ComponentFacts]) -> Option<(usize, usize)> {
    let kept: Vec<&ComponentFacts> = facts.iter().filter(|f| f.kept).collect();
    if kept.len() < 2 {
        return None;
    }
    let best = |key: &dyn Fn(&ComponentFacts) -> Option<f64>, skip: Option<usize>| {
        (kept.iter())
            .filter(|f| Some(f.dc) != skip)
            .filter_map(|f| Some((f.dc, key(f)?)))
            .filter(|(_, v)| v.is_finite())
            .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)))
            .map(|(dc, _)| dc)
    };
    let Some(x) = best(&|f| f.rho, None) else {
        return Some((kept[0].dc, kept[1].dc));
    };
    let y = best(&|f| f.lineage_eta2, Some(x)).or_else(|| best(&|f| f.rho, Some(x)))?;
    Some((x, y))
}

/// The standard deviation of the finite values of `v`.
fn sd(v: &[f64]) -> f64 {
    let finite: Vec<f64> = v.iter().copied().filter(|v| v.is_finite()).collect();
    let n = finite.len() as f64;
    if n == 0.0 {
        return f64::NAN;
    }
    let mean = finite.iter().sum::<f64>() / n;
    (finite.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n).sqrt()
}

/// The share of `v`'s variance between the groups of `group` (η²), over
/// the cells in a group (code ≥ 0) with a finite value; none without two
/// groups or with no variance.
fn eta2(v: &[f64], group: &[i32]) -> Option<f64> {
    let cells: Vec<(i32, f64)> = (group.iter().zip(v))
        .filter(|(&g, v)| g >= 0 && v.is_finite())
        .map(|(&g, &v)| (g, v))
        .collect();
    let mut sums: std::collections::BTreeMap<i32, (f64, f64)> = Default::default();
    for &(g, v) in &cells {
        let s = sums.entry(g).or_default();
        *s = (s.0 + v, s.1 + 1.0);
    }
    if sums.len() < 2 {
        return None;
    }
    let mean = cells.iter().map(|c| c.1).sum::<f64>() / cells.len() as f64;
    let total: f64 = cells.iter().map(|c| (c.1 - mean).powi(2)).sum();
    let between: f64 = (sums.values())
        .map(|&(s, n)| n * (s / n - mean).powi(2))
        .sum();
    (total > 0.0).then(|| between / total)
}
