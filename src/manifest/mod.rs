//! The run manifest: its schema and paths ([`run`]), a pinto run read as one
//! ([`pinto`]), and the commands that read a run from it and record their
//! outputs ([`annotate`], [`lineage`]), and review / relabel rounds of an
//! annotation ([`rounds`]).

pub mod annotate;
pub mod lineage;
pub mod pinto;
pub mod rounds;
pub mod run;
