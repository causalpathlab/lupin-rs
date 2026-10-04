//! Progress a batch command reports to the TUI that started it. With
//! `LUPIN_PROGRESS=1` in its environment (the TUI sets it for its
//! children) a command writes lines to stderr of the form
//!
//! ```text
//! @progress <done> <end> <total> <stage>
//! ```
//!
//! `done`, `end` (where the current stage ends) and `total` are
//! non-negative numbers in any unit (the fractions are what count) and
//! `stage` is the rest of the line, free text. Without the
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
fn line(done: f64, end: f64, total: f64, stage: &str) -> String {
    format!(
        "{PREFIX}{done:.4} {end:.4} {total:.4} {}",
        stage.replace('\n', " ")
    )
}

/// A progress line read back.
#[derive(Debug, Clone, PartialEq)]
pub struct Report {
    pub done: f64,
    /// Where the current stage ends.
    pub end: f64,
    pub total: f64,
    pub stage: String,
}

/// The report on a progress line; `None` for any other line.
pub fn parse(line: &str) -> Option<Report> {
    let rest = line.strip_prefix(PREFIX)?;
    let mut it = rest.splitn(4, ' ');
    let mut num = || -> Option<f64> {
        let v: f64 = it.next()?.parse().ok()?;
        (v.is_finite() && v >= 0.0).then_some(v)
    };
    let (done, end, total) = (num()?, num()?, num()?);
    let stage = it.next().unwrap_or("").to_string();
    (total > 0.0).then(|| {
        let done = done.min(total);
        Report {
            done,
            end: end.clamp(done, total),
            total,
            stage,
        }
    })
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
        line(
            self.before(i) + w * part,
            self.before(i) + w,
            self.total,
            &stage,
        )
    }

    /// Every stage is done.
    pub fn finish(&self) {
        if enabled() {
            eprintln!("{}", line(self.total, self.total, self.total, "done"));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_progress_line_reads_back_and_others_do_not() {
        assert_eq!(
            parse("@progress 1.5000 2.0000 6.0000 kNN graph (3 of 4)"),
            Some(Report {
                done: 1.5,
                end: 2.0,
                total: 6.0,
                stage: "kNN graph (3 of 4)".into()
            })
        );
        assert_eq!(parse("@progress 2 3 4").unwrap().stage, "");
        for bad in [
            "[INFO] wrote x",
            "@progress",
            "@progress 1 2",
            "@progress a b c d",
            "@progress 1 1 0 zero total",
            "@progress -1 1 4 negative",
            "@progress NaN 1 4 nan",
            " @progress 1 2 4 indented",
        ] {
            assert_eq!(parse(bad), None, "{bad}");
        }
        // More done than the total is capped, and the stage's end sits
        // between the two.
        let r = parse("@progress 9 1 4 x").unwrap();
        assert_eq!((r.done, r.end), (4.0, 4.0));
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
        let starts: Vec<f64> = (0..STAGES.stages.len()).map(|i| read(i).done).collect();
        assert!(starts.windows(2).all(|w| w[0] < w[1]), "{starts:?}");
        assert_eq!(read(STAGE_DIFFUSION).stage, "diffusion components");
        // Each stage ends where the next starts.
        for i in 1..STAGES.stages.len() {
            assert_eq!(read(i - 1).end, read(i).done);
        }
        assert_eq!(
            parse(&STAGES.line(STAGES.stages.len(), 0, 1, None))
                .unwrap()
                .done,
            STAGES.total
        );
    }

    #[test]
    fn stages_add_up_in_order() {
        static S: &[(&str, f64)] = &[("one", 1.0), ("two", 3.0)];
        let s = Stages::new(S);
        let read = |l: String| parse(&l).unwrap();
        assert_eq!(read(s.line(0, 0, 1, None)).done, 0.0);
        assert_eq!(read(s.line(1, 0, 1, None)).done, 1.0);
        let r = read(s.line(1, 1, 3, Some("1 of 3")));
        assert_eq!((r.done, r.end, r.total), (2.0, 4.0, 4.0));
        assert_eq!(r.stage, "two (1 of 3)");
        assert_eq!(read(s.line(1, 3, 3, None)).done, 4.0);
    }
}
