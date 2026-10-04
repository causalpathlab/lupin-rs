//! Progress a batch command reports to the TUI that started it. With
//! `LUPIN_PROGRESS=1` in its environment (the TUI sets it for its
//! children) a command writes lines to stderr of the form
//!
//! ```text
//! @progress <done> <total> <stage>
//! ```
//!
//! `done` and `total` are non-negative numbers in any unit (the fraction is
//! what counts) and `stage` is the rest of the line, free text. Without the
//! variable nothing is written, so batch users see no change.

use std::sync::OnceLock;

/// The environment variable that turns progress lines on.
pub const ENV: &str = "LUPIN_PROGRESS";

const PREFIX: &str = "@progress ";

fn enabled() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var(ENV).is_ok_and(|v| v == "1"))
}

/// The line for `done` of `total` at `stage`.
fn line(done: f64, total: f64, stage: &str) -> String {
    format!("{PREFIX}{done:.4} {total:.4} {}", stage.replace('\n', " "))
}

/// `(done, total, stage)` from a progress line; `None` for any other line.
pub fn parse(line: &str) -> Option<(f64, f64, String)> {
    let rest = line.strip_prefix(PREFIX)?;
    let mut it = rest.splitn(3, ' ');
    let done: f64 = it.next()?.parse().ok()?;
    let total: f64 = it.next()?.parse().ok()?;
    let stage = it.next().unwrap_or("").to_string();
    (done.is_finite() && total.is_finite() && done >= 0.0 && total > 0.0).then_some((
        done.min(total),
        total,
        stage,
    ))
}

/// A command's stages, each weighted by its rough share of the run time.
pub struct Stages {
    stages: &'static [(&'static str, f64)],
    total: f64,
}

impl Stages {
    pub const fn new(stages: &'static [(&'static str, f64)]) -> Self {
        let mut total = 0.0;
        let mut i = 0;
        while i < stages.len() {
            total += stages[i].1;
            i += 1;
        }
        Self { stages, total }
    }

    /// The weight of the stages before stage `i`.
    fn before(&self, i: usize) -> f64 {
        self.stages[..i.min(self.stages.len())]
            .iter()
            .map(|s| s.1)
            .sum()
    }

    /// Stage `i` has started.
    pub fn start(&self, i: usize) {
        self.within(i, 0, 1, None);
    }

    /// `k` of `n` steps of stage `i` are done; `note` adds to its name.
    pub fn within(&self, i: usize, k: usize, n: usize, note: Option<&str>) {
        if !enabled() {
            return;
        }
        eprintln!("{}", self.line(i, k, n, note));
    }

    fn line(&self, i: usize, k: usize, n: usize, note: Option<&str>) -> String {
        let (name, w) = self.stages.get(i).copied().unwrap_or(("done", 0.0));
        let part = if n == 0 {
            0.0
        } else {
            k.min(n) as f64 / n as f64
        };
        let stage = match note {
            Some(note) => format!("{name} ({note})"),
            None => name.to_string(),
        };
        line(self.before(i) + w * part, self.total, &stage)
    }

    /// Every stage is done.
    pub fn finish(&self) {
        if enabled() {
            eprintln!("{}", line(self.total, self.total, "done"));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_progress_line_reads_back_and_others_do_not() {
        assert_eq!(
            parse("@progress 1.5000 6.0000 kNN graph (3 of 4)"),
            Some((1.5, 6.0, "kNN graph (3 of 4)".into()))
        );
        assert_eq!(parse("@progress 2 4"), Some((2.0, 4.0, String::new())));
        for bad in [
            "[INFO] wrote x",
            "@progress",
            "@progress a b c",
            "@progress 1 0 zero total",
            "@progress -1 4 negative",
            "@progress NaN 4 nan",
            " @progress 1 4 indented",
        ] {
            assert_eq!(parse(bad), None, "{bad}");
        }
        // More done than the total is capped.
        assert_eq!(parse("@progress 9 4 x").unwrap().0, 4.0);
    }

    #[test]
    fn nothing_is_written_without_the_variable() {
        // The test process is not a TUI child.
        if std::env::var(ENV).is_err() {
            assert!(!enabled());
        }
    }

    #[test]
    fn the_trajectory_stages_run_in_order_and_diffusion_is_where_it_says() {
        use crate::trajectory::run::{STAGES, STAGE_DIFFUSION};
        let read = |i| parse(&STAGES.line(i, 0, 1, None)).unwrap();
        let starts: Vec<f64> = (0..STAGES.stages.len()).map(|i| read(i).0).collect();
        assert!(starts.windows(2).all(|w| w[0] < w[1]), "{starts:?}");
        assert_eq!(read(STAGE_DIFFUSION).2, "diffusion components");
        assert_eq!(
            parse(&STAGES.line(STAGES.stages.len(), 0, 1, None))
                .unwrap()
                .0,
            STAGES.total
        );
    }

    #[test]
    fn stages_add_up_in_order() {
        static S: &[(&str, f64)] = &[("one", 1.0), ("two", 3.0)];
        let s = Stages::new(S);
        let read = |l: String| parse(&l).unwrap();
        assert_eq!(read(s.line(0, 0, 1, None)).0, 0.0);
        assert_eq!(read(s.line(1, 0, 1, None)).0, 1.0);
        let (d, t, stage) = read(s.line(1, 1, 3, Some("1 of 3")));
        assert_eq!((d, t), (2.0, 4.0));
        assert_eq!(stage, "two (1 of 3)");
        assert_eq!(read(s.line(1, 3, 3, None)).0, 4.0);
    }
}
