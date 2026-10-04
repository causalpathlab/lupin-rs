//! `lupin trajectory`: an explicit, supervised prior over the run's cell
//! types, checked against the data and used to order cells by diffusion
//! pseudotime (`docs/trajectory-plan.md`): the scanpy-matched diffusion map,
//! pseudotime and PAGA connectivity, the prior, and the command over them.
//! The TUI order view and the plot build on these.

pub(crate) mod diffusion;
pub(crate) mod edges;
pub(crate) mod figures;
pub(crate) mod prior;
pub(crate) mod run;
pub(crate) mod type_connectivity;

#[cfg(test)]
mod reference;

/// The distinct labels, sorted, and each cell's index into them.
pub(crate) fn encode_groups<S: AsRef<str>>(labels: &[S]) -> (Vec<Box<str>>, Vec<usize>) {
    let mut names: Vec<&str> = labels.iter().map(AsRef::as_ref).collect();
    names.sort_unstable();
    names.dedup();
    let group = labels
        .iter()
        .map(|l| {
            names
                .binary_search(&l.as_ref())
                .expect("label is among the names")
        })
        .collect();
    (names.into_iter().map(Into::into).collect(), group)
}
