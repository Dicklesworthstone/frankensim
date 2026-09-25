//! A retained coarse endpoint is an input to a separately funded fine study.
//! Reuse the original ledger seals, prolongator and projected mechanics owner.
use super::*;
use fs_topols::refinement::prolongate_level_set;
use fs_topols::refinement::projected::{VolumeRefinementStage, refine_controlled};

#[derive(Debug, Clone, Copy)]
pub(super) enum Stage {
    Read(ContentHash),
    Regions(DesignRegionStage),
    Transfer,
    SourceSetup(ProjectedSetupStage),
    Refine(VolumeRefinementStage),
}

pub(super) fn parse_ref(fields: &[Node]) -> Result<Option<ContentHash>> {
    let Some(pair) = fields.windows(2).find(|pair|
        matches!(&pair[0].kind, NodeKind::Keyword(key) if key == "refine-from"))
    else { return Ok(None) };
    let NodeKind::Str(pointer) = &pair[1].kind
        else { return Err(malformed("refine-from requires a quoted study receipt identifier")); };
    let hash = pointer.strip_prefix("study-").and_then(ContentHash::from_hex)
        .ok_or_else(|| malformed("refine-from requires study- followed by a 64-digit receipt hash"))?;
    Ok(Some(hash))
}

fn policy(spec: &ElasticitySpec) -> Result<&Controls> {
    match &spec.projected {
        Some(ProjectedControls::Volume(policy)) => Ok(policy),
        _ => Err(malformed("refinement currently requires projected-volume on both studies")),
    }
}

struct Source {
    loaded: Loaded,
    spec: ElasticitySpec,
    phi: GridSdf,
    report: OptimizeReport,
    current: Measured,
}

fn source(target: &ElasticitySpec, ledger: &Ledger) -> Result<Source> {
    let target_policy = policy(target)?;
    let hash = target_policy.refine_from.ok_or_else(|| malformed("missing refinement source"))?;
    let loaded = load(ledger, &format!("study-{}", hash.to_hex()))?;
    if loaded.value.str_field("driver") != Some(DRIVER) {
        return Err(malformed("refinement source is not the native 2-D elasticity producer"));
    }
    let source_bytes = linked(ledger, &loaded.value, "source", "study-source")?;
    let text = std::str::from_utf8(&source_bytes).map_err(|_| malformed("source project is not UTF-8"))?;
    let spec = crate::study::elasticity::parse(text)?;
    if loaded.value.str_field("study_id") != Some(spec.id.to_hex().as_str()) {
        return Err(malformed("refinement source identity differs from its retained project"));
    }
    let old = policy(&spec)?;
    if settings(target, target.steps).level != settings(&spec, spec.steps).level + 1 {
        return Err(malformed("refine-from requires exactly one finer mesh level"));
    }
    // Compare the complete canonical numerical declaration. Only metadata,
    // new work budgets/update count, the grid level and mesh-check policy may
    // differ. This also preserves hole definitions and every protected region.
    let mut comparison = spec.clone();
    comparison.base.metadata = target.base.metadata.clone();
    comparison.base.budgets = target.base.budgets.clone();
    comparison.base.physics.as_mut().expect("admitted physics").mesh_level =
        target.base.physics.as_ref().expect("admitted physics").mesh_level;
    comparison.steps = target.steps;
    comparison.wall_s = target.wall_s;
    comparison.memory_bytes = target.memory_bytes;
    comparison.max_iterations = target.max_iterations;
    let Some(ProjectedControls::Volume(comparison_policy)) = &mut comparison.projected
        else { return Err(malformed("refinement source policy changed")); };
    comparison_policy.refine_from = target_policy.refine_from;
    comparison_policy.resolution = target_policy.resolution;
    if crate::study::elasticity::canonical(&comparison) != target.canonical {
        return Err(malformed("refinement must preserve loads, material, seed, geometry declarations, area and search policy"));
    }
    let design = document(&linked(ledger, &loaded.value, "design", "study-design")?)?;
    let rows = document(&linked(ledger, &loaded.value, "iterations", "study-iterations")?)?;
    let (phi, report) = decode(&spec, &loaded.value, &design, &rows)?;
    let history = VolumeEvidence::read(loaded.value.path(&["continuation", "constraints"])
        .ok_or_else(|| malformed("source has no volume-only accepted history"))?, &report, old)?;
    let current = history.current();
    if current.snapshot != snapshot(&phi) {
        return Err(malformed("refinement source geometry differs from its measured endpoint"));
    }
    if let Some(origin) = &history.origin {
        if origin.fine_level != settings(&spec, spec.steps).level {
            return Err(malformed("source refinement level differs from its project"));
        }
    }
    if !matches!(loaded.value.str_field("status"),
        Some("running" | "completed" | "cancelled" | "budget-exhausted" | "no-feasible-descent" | "mesh-unresolved")) {
        return Err(malformed("refinement source has no admitted retained state"));
    }
    Ok(Source { loaded, spec, phi, report, current })
}

/// Reconstruct original prescriptions, never derive new pins from a free
/// endpoint. Each ancestry edge decreases the level, with an independent cap.
/// No PDE is run while reconstructing fixed nodes for normal fine-grid resume.
pub(super) fn fixed<B>(spec: &ElasticitySpec, ledger: &Ledger,
    mut control: impl FnMut(Stage) -> ControlFlow<B>) -> Result<ControlFlow<B, Vec<(usize, f64)>>> {
    fixed_inner(spec, ledger, 0, &mut control)
}
fn fixed_inner<B>(spec: &ElasticitySpec, ledger: &Ledger, depth: usize,
    control: &mut impl FnMut(Stage) -> ControlFlow<B>) -> Result<ControlFlow<B, Vec<(usize, f64)>>> {
    if depth >= 8 { return Err(malformed("refinement ancestry exceeds the bounded grid envelope")); }
    let current = policy(spec)?;
    let Some(hash) = current.refine_from else {
        return match regions::prepare(spec, &current.regions, |stage| control(Stage::Regions(stage)))? {
            ControlFlow::Continue(prepared) => Ok(ControlFlow::Continue(prepared.fixed_nodes)),
            ControlFlow::Break(reason) => Ok(ControlFlow::Break(reason)),
        };
    };
    if let ControlFlow::Break(reason) = control(Stage::Read(hash)) {
        return Ok(ControlFlow::Break(reason));
    }
    let source = source(spec, ledger)?;
    let inherited = match fixed_inner(&source.spec, ledger, depth + 1, control)? {
        ControlFlow::Continue(fixed) => fixed,
        ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
    };
    if let ControlFlow::Break(reason) = control(Stage::Transfer) {
        return Ok(ControlFlow::Break(reason));
    }
    let transfer = prolongate_level_set(&source.phi, &inherited)
        .map_err(|error| malformed(&error.to_string()))?;
    Ok(ControlFlow::Continue(transfer.fixed_nodes))
}

/// The source endpoint is independently re-admitted before refining; no old
/// objective, source schedule, or cross-grid improvement becomes a fine result.
pub(super) fn start<B>(spec: &ElasticitySpec, ledger: &Ledger,
    mut control: impl FnMut(Stage) -> ControlFlow<B>)
    -> Result<ControlFlow<B, (ProjectedOptimizer, Origin)>> {
    let hash = policy(spec)?.refine_from.ok_or_else(|| malformed("missing refinement source"))?;
    if let ControlFlow::Break(reason) = control(Stage::Read(hash)) {
        return Ok(ControlFlow::Break(reason));
    }
    let source = source(spec, ledger)?;
    let inherited = match fixed(&source.spec, ledger, &mut control)? {
        ControlFlow::Continue(fixed) => fixed,
        ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
    };
    let old = policy(&source.spec)?;
    let checkpoint = OptimizeCheckpoint::restore(source.phi.clone(), fixture(&source.spec),
        settings(&source.spec, source.spec.steps), source.report.rows.len(),
        source.report.ell.last().copied().unwrap_or(source.spec.ell0))
        .map_err(|error| malformed(&error.to_string()))?;
    let coarse = match ProjectedOptimizer::from_checkpoint_controlled(&checkpoint, inherited,
        old.area, old.search, |stage| control(Stage::SourceSetup(stage)))
        .map_err(|error| malformed(&error.to_string()))? {
        ControlFlow::Continue(coarse) => coarse,
        ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
    };
    if !Measured::from(coarse.current()).same(source.current) {
        return Err(malformed("source mechanics do not reproduce the retained endpoint; no refined study published"));
    }
    let (fine, report) = match refine_controlled(&coarse, spec.steps,
        |stage| control(Stage::Refine(stage))).map_err(|error| malformed(&error.to_string()))? {
        ControlFlow::Continue(result) => result,
        ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
    };
    let origin = Origin { receipt: source.loaded.hash, coarse_level: report.coarse_level,
        fine_level: report.fine_level, source_updates: report.source_updates,
        coarse: report.coarse_endpoint.into(), transferred_area: report.transferred_area,
        projection_change: report.max_projection_change };
    Ok(ControlFlow::Continue((fine, origin)))
}

#[derive(Debug, Clone)]
pub(super) struct Origin {
    receipt: ContentHash,
    coarse_level: u32,
    pub(super) fine_level: u32,
    source_updates: usize,
    coarse: Measured,
    transferred_area: f64,
    projection_change: f64,
}
impl Origin {
    pub(super) fn json_field(&self) -> String {
        format!(concat!(",\"refinement_origin\":{{\"source_receipt\":\"{}\",",
            "\"coarse_level\":{},\"fine_level\":{},\"source_updates\":{},",
            "\"coarse_endpoint\":{},\"transferred_area_m2\":{:.17e},",
            "\"max_projection_field_change\":{:.17e},\"comparison\":\"new-fine-baseline\"}}"),
            self.receipt.to_hex(), self.coarse_level, self.fine_level, self.source_updates,
            self.coarse.json(), self.transferred_area, self.projection_change)
    }
    pub(super) fn html(&self) -> String {
        format!("<p>Refined from retained study-{} (level {}, {} accepted updates) to level {}. Transferred numerical area {:.8e} m²; largest area-restoration field correction {:.8e}. This study has a new work budget and an independently solved fine-grid baseline. Cross-grid compliance changes and baseline projection are not optimization improvement.</p>",
            self.receipt.to_hex(), self.coarse_level, self.source_updates, self.fine_level,
            self.transferred_area, self.projection_change)
    }
}

pub(super) fn read(value: &JsonValue, policy: &Controls) -> Result<Option<Origin>> {
    let Some(hash) = policy.refine_from else {
        if value.get("refinement_origin").is_some() { return Err(malformed("undeclared refinement origin")); }
        return Ok(None);
    };
    let origin = value.get("refinement_origin").ok_or_else(|| malformed("missing refinement origin"))?;
    let coarse_level = u32::try_from(integer(origin, "coarse_level")?).map_err(|_| malformed("invalid coarse level"))?;
    let fine_level = u32::try_from(integer(origin, "fine_level")?).map_err(|_| malformed("invalid fine level"))?;
    let source_updates = integer(origin, "source_updates")?;
    let transferred_area = number(origin, "transferred_area_m2")?.0;
    let projection_change = number(origin, "max_projection_field_change")?.0;
    if origin.str_field("source_receipt") != Some(hash.to_hex().as_str())
        || origin.str_field("comparison") != Some("new-fine-baseline")
        || !(1..=7).contains(&coarse_level) || fine_level != coarse_level + 1
        || source_updates > 32 || transferred_area <= 0.0 || projection_change < 0.0 {
        return Err(malformed("refinement origin differs from the declared source or numerical envelope"));
    }
    let coarse = Measured::read(origin.get("coarse_endpoint")
        .ok_or_else(|| malformed("missing coarse endpoint"))?, policy)?;
    Ok(Some(Origin { receipt: hash, coarse_level, fine_level, source_updates,
        coarse, transferred_area, projection_change }))
}

/// Make the canonical fine source an actual ledger-derived artifact of the
/// coarse receipt. Ordinary study operations consume these same source bytes,
/// so lineage/GC follows real edges without modifying any sealed study operation.
pub(super) fn retain_source(spec: &ElasticitySpec, ledger: &Ledger, origin: &Origin) -> Result<()> {
    ledger.begin()?;
    let result = (|| -> Result<()> {
        let seed = spec.base.seeds.as_ref().expect("admitted seed").root.to_le_bytes();
        let versions = format!("{{\"driver\":{DRIVER:?},\"crate\":{:?}}}", env!("CARGO_PKG_VERSION"));
        let budget = format!("{{\"wall_s\":{},\"memory_bytes\":{},\"max_iterations\":{}}}",
            spec.wall_s, spec.memory_bytes, spec.max_iterations);
        let ir = format!("{{\"operation\":\"projected-volume-refinement-source\",\"units\":\"SI\",\"study_id\":\"{}\",\"source_receipt\":\"{}\"}}",
            spec.id.to_hex(), origin.receipt.to_hex());
        let op = ledger.begin_op(Some(spec.id.as_bytes()), &ir,
            &FiveExplicits { seed: &seed, versions: &versions, budget: &budget,
                capability: "{\"ops\":[\"optimization.marquee-topopt\",\"geometry.sdf\",\"physics.cutfem\"]}" }, 0)?;
        ledger.link(op, &origin.receipt, EdgeRole::In)?;
        let source = ledger.put_artifact("study-source", spec.canonical.as_bytes(), None)?;
        ledger.link(op, &source.hash, EdgeRole::Out)?;
        if ledger.artifact_output_seal(&source.hash)?.is_none() {
            ledger.seal_artifact_output(&source.hash, op)?;
        }
        ledger.finish_op(op, OpOutcome::Ok, None, 1)?;
        Ok(())
    })();
    match result {
        Ok(()) => match ledger.commit() {
            Ok(()) => Ok(()),
            Err(error) => { ledger.rollback()?; Err(error.into()) }
        },
        Err(error) => { ledger.rollback()?; Err(error) }
    }
}
