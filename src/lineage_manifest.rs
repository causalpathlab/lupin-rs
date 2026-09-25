//! Run-manifest glue for the lineage commands: resolve `-f run.senna.json`
//! into the plain paths and loaded matrices [`crate::lineage`] takes, then
//! write the artifact paths it returns back into the manifest.

use crate::run_manifest::resolve;
use anyhow::Result;
use log::info;

use crate::lineage::args::LineageArgs;
use crate::lineage::pseudotime::{
    run_pseudotime, PseudotimeArgs, PseudotimeInputs, PseudotimeOutputs,
};
use crate::lineage::run::{run_lineage, LineageInputs};
use crate::lineage::{LatentContract, RunTables};

use crate::marker_embedding::load_marker_feature_embedding;
use crate::run_manifest::{derive_out_prefix, load, manifest_file, rel_to_manifest, Loaded};

/// What `{prefix}`'s manifest says about its per-cell tables, in the form the
/// θ resolver takes, plus the loaded manifest. A missing or unreadable
/// manifest is not fatal: `--theta-from auto` then assumes the geometry rather
/// than reading it, which is what [`LatentContract::unknown`] means.
fn load_contract(prefix: &str) -> (LatentContract, Option<Loaded>) {
    match load(prefix) {
        Ok(loaded) => {
            let contract = LatentContract {
                latent_is_log_simplex: loaded.manifest.kind.latent_is_log_simplex(),
                kind: Some(loaded.manifest.kind.to_string().into_boxed_str()),
                source: loaded.file.to_string_lossy().into(),
            };
            (contract, Some(loaded))
        }
        Err(_) => {
            let source = manifest_file(prefix).to_string_lossy().into_owned();
            (LatentContract::unknown(source), None)
        }
    }
}

/// `lupin lineage`: resolve the manifest-derived inputs,
/// then run the fit.
pub fn run_lineage_from_manifest(args: &LineageArgs) -> Result<()> {
    let prefix = args.from.as_ref();
    let (contract, manifest) = load_contract(prefix);
    let (tables, feature_embedding) = match manifest {
        Some(loaded) => {
            let feature_embedding = args
                .markers
                .is_some()
                .then(|| load_marker_feature_embedding(&loaded.manifest, &loaded.dir, prefix))
                .transpose()?;
            (run_tables(&loaded), feature_embedding)
        }
        None => {
            anyhow::ensure!(
                args.markers.is_none(),
                "--markers needs a readable run.senna.json under {prefix} to locate \
                 feature_coembedding"
            );
            (RunTables::at_prefix(&derive_out_prefix(prefix)), None)
        }
    };
    let inputs = LineageInputs {
        tables,
        contract,
        feature_embedding,
    };
    run_lineage(args, &inputs)
}

/// The per-cell tables the manifest records, resolved against its directory;
/// a slot it leaves empty, and velocity (never recorded), fall back to the
/// file-name convention next to the manifest itself.
fn run_tables(loaded: &Loaded) -> RunTables {
    let guess = RunTables::at_prefix(&loaded.run_prefix());
    let out = &loaded.manifest.outputs;
    let pick =
        |rel: Option<&str>, fallback: String| rel.map_or(fallback, |r| resolve(&loaded.dir, r));
    RunTables {
        latent: pick(out.latent.as_deref(), guess.latent),
        cell_embedding: pick(out.cell_embedding.as_deref(), guess.cell_embedding),
        velocity: guess.velocity,
    }
}

/// `lupin pseudotime`: resolve `--latent` / `--from`,
/// run the fit, and record the artifacts it wrote in the manifest.
pub fn run_pseudotime_from_manifest(args: &PseudotimeArgs) -> Result<()> {
    let manifest_ctx = args.from.as_deref().map(load).transpose()?;
    let inputs = PseudotimeInputs {
        latent: resolve_latent_path(args, manifest_ctx.as_ref())?,
        cell_coords: manifest_ctx.as_ref().and_then(resolve_cell_coords),
    };

    let outputs = run_pseudotime(args, &inputs)?;

    if let Some(ctx) = manifest_ctx {
        update_manifest(ctx, &outputs)?;
    }
    Ok(())
}

fn resolve_latent_path(args: &PseudotimeArgs, ctx: Option<&Loaded>) -> Result<String> {
    if let Some(p) = args.latent.as_deref() {
        return Ok(p.to_string());
    }
    let ctx = ctx.ok_or_else(|| anyhow::anyhow!("either --latent or --from is required"))?;
    let rel = ctx.manifest.outputs.geometry_latent().ok_or_else(|| {
        anyhow::anyhow!(
            "manifest {} has no outputs.cell_embedding or outputs.latent",
            ctx.file.display()
        )
    })?;
    Ok(resolve(&ctx.dir, rel))
}

/// The run's 2D cell layout, which is what lets the pseudotime fit also
/// project its centroids into layout space. Absent is ordinary — the layout is
/// a separate command — so this says how to get one rather than failing.
fn resolve_cell_coords(ctx: &Loaded) -> Option<String> {
    let Some(rel) = ctx.manifest.layout.cell_coords.as_deref() else {
        info!(
            "manifest has no layout.cell_coords; skipping 2D centroid projection \
             (run `senna layout phate --from ...` first to enable plot overlay)"
        );
        return None;
    };
    Some(resolve(&ctx.dir, rel))
}

fn update_manifest(mut ctx: Loaded, outputs: &PseudotimeOutputs) -> Result<()> {
    let rel = |p: &str| rel_to_manifest(&ctx.dir, p);
    ctx.manifest.pseudotime.pseudotime = Some(rel(&outputs.pseudotime));
    ctx.manifest.pseudotime.nodes_latent = Some(rel(&outputs.nodes_latent));
    ctx.manifest.pseudotime.nodes_2d = outputs.nodes_2d.as_deref().map(rel);
    ctx.manifest.pseudotime.edges = Some(rel(&outputs.edges));
    ctx.manifest.pseudotime.root_node = Some(outputs.root);
    // Tree layout is now produced by `senna layout tree`; clear any stale
    // paths from a previous run so downstream readers don't pick up a tree
    // that no longer matches this pseudotime fit.
    ctx.manifest.pseudotime.tree_cell_coords = None;
    ctx.manifest.pseudotime.tree_nodes_2d = None;
    ctx.manifest.save(&ctx.file)?;
    Ok(())
}

#[cfg(test)]
#[path = "lineage_manifest_tests.rs"]
mod tests;
