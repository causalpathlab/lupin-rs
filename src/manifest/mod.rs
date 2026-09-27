//! The run manifest: its schema and paths ([`run`]), a pinto run read as one
//! ([`pinto`]), and the commands that read a run from it and record their
//! outputs ([`annotate`], [`lineage`]).

pub mod annotate;
pub mod lineage;
pub mod pinto;
pub mod run;
