//! The run manifest: its schema and paths ([`run`]), a pinto run read as one
//! ([`pinto`]), and the commands that read a run from it and record their
//! outputs ([`annotate`], [`lineage`]), and review / relabel rounds of an
//! annotation ([`rounds`]).

pub mod annotate;
pub mod ask;
pub mod data_files;
pub mod first_round;
pub mod lineage;
pub mod ontology;
pub mod pinto;
pub mod recalibrate;
pub mod rounds;
pub mod run;
