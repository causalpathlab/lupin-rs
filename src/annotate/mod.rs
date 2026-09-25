//! Cell-type annotation (enrichment, projection, ontology) over plain file
//! paths and loaded matrices. [`crate::annotate_manifest`] resolves a
//! `run.senna.json` into the input structs here and records the outputs.

pub mod aggregate;
pub mod args;
pub mod by_enrichment;
pub mod by_projection;
mod go_signature;
pub mod inputs;
pub mod markers;
pub mod ontology;
pub mod outputs;
