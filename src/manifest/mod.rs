//! The run manifest: its schema and paths ([`run`]), and the commands that
//! read a run from it and record their outputs back ([`annotate`], [`lineage`]).

pub mod annotate;
pub mod lineage;
pub mod run;
