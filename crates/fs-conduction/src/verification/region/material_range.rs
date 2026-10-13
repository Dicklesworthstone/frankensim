//! Admit finite constant-material domains using an independent continuum barrier.
//! The extension below is a proof device, never a replacement physical material.
use super::*;
use crate::{ConductivityModel, ElementMaterials, MaterialId, MaterialTable, TemperatureSpan};
use std::collections::BTreeMap;

/// Regional mean enclosure with the continuum range that established material
/// validity. The original materials and their property receipts are unchanged.
#[derive(Debug)]
pub struct SlabValidatedRegionMean {
    /// Original regional-integral dual and equilibrated primal/dual majorants.
    pub region: RegionMeanFieldBound,
    /// Outward lower/upper temperatures for the entire continuum solution, K.
    /// These are NOT extrema of the supplied nodal field.
    pub continuum_temperature_range_k: [f64; 2],
    /// Intersection of the temperature domains of all USED material models.
    pub material_temperature_span: TemperatureSpan,
    /// Axis of the independently checked comparison slab, 0, 1 or 2.
    pub axis: usize,
    /// True when the cooling end is the low-coordinate end of this axis.
    pub reverse: bool,
}

fn constant_extension(model: &ConductivityModel) -> Result<ConductivityModel> {
    if model.is_temperature_dependent() {
        return Err(unsupported("slab validity cannot freeze a temperature-dependent conductivity"));
    }
    let witness = match model.temperature_span() {
        TemperatureSpan::Unbounded => 0.0,
        TemperatureSpan::Sampled { low, high } => {
            if !low.is_finite() || !high.is_finite() || low > high {
                return Err(TetError::Invalid("constant material has an empty or nonfinite domain").into());
            }
            low * 0.5 + high * 0.5
        }
    };
    // is_temperature_dependent checks every retained knot, without a tolerance.
    // Preserve the complete tensor and frame, never average principal values.
    Ok(ConductivityModel::constant_tensor(model.tensor_at(witness)?)?)
}

/// Bound a component mean when constant conductivities have finite temperature
/// domains. The independently verified slab barrier must lie inside EVERY used
/// material's domain before any dual solve can run.
///
/// Construct a mathematical constant extension, establish its whole-continuum
/// temperature range with `slab_temperature_range`, then reuse the regional
/// majorant. Since that range lies within the original laws, the extended
/// solution also solves the original model. No original table is mutated or
/// evaluated outside its span, and no primal solve or extrapolated property
/// receipt is fabricated. A flat bounded table is not called globally valid.
///
/// The original slab producer admits bounding-box-plane exteriors, zero
/// non-cooling Neumann flux, a fully Dirichlet/Robin cooling end, and one common
/// uncoupled axial conductivity. Tangential tensor blocks may differ. Contact,
/// radiation, convection feedback and varying conductivity need separate proofs.
/// As with the underlying verifier, conforming non-overlapping volume geometry
/// is a caller obligation. Supplied boundary rows must belong to this mesh.
///
/// `axis` and `reverse` explicitly select the proposed barrier; unsupported
/// geometry or a barrier wider than material validity yields no bound. It does
/// not prove that the physical solution actually leaves the validity domain.
/// This is a nominal PDE bound, not a material-uncertainty or physical certificate.
/// Existing matrix, flux and cancellation limits remain in force.
#[allow(clippy::too_many_arguments)]
pub fn bound_temperature_region_mean_with_slab_validity(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, temperature: &[f64],
    cells: &[usize], axis: usize, reverse: bool, config: RegionMeanConfig,
) -> Result<SlabValidatedRegionMean> {
    let selected = selection(cx, problem, cells, config)?;
    if axis >= 3 { return Err(TetError::Invalid("comparison axis must be 0, 1 or 2").into()); }
    let ne = problem.mesh.element_count();
    let n = problem.mesh.vertex_count();
    if ne > config.flux.max_cells || n > config.flux.max_cells.saturating_mul(4) {
        return Err(TetError::Budget.into());
    }
    if temperature.len() != n { return Err(TetError::Invalid("candidate temperature length").into()); }
    if let Some(assigned) = problem.element_materials { assigned.validate_for(problem.mesh)?; }
    let mut models = BTreeMap::new();
    let mut labels = Vec::with_capacity(ne);
    let mut span = TemperatureSpan::Unbounded;
    for e in 0..ne {
        poll(cx)?;
        let (id, model) = match problem.element_materials {
            Some(assigned) => (assigned.of_element()[e], assigned.model_for(e)?),
            None => (MaterialId(0), problem.material),
        };
        if let std::collections::btree_map::Entry::Vacant(entry) = models.entry(id) {
            entry.insert(constant_extension(model)?);
            span = span.intersect(model.temperature_span());
        }
        labels.push(id);
    }
    if let TemperatureSpan::Sampled { low, high } = span {
        if low > high { return Err(unsupported("used materials have no common temperature domain")); }
    }
    let extension = ElementMaterials::new(MaterialTable::new(models)?, labels)?;
    let extended = ConductionProblem {
        material: extension.model_for(0)?, element_materials: Some(&extension), ..problem
    };
    // The original binder rejects undeclared contacts; an empty contact list
    // cannot silently turn independent coincident traces into continuity.
    let (admitted, interfaces) = contact::admit_contacts(cx, extended, &[], config.flux)?;
    let range = admitted.problem(problem.mesh.positions()).slab_temperature_range(
        temperature, axis, reverse, config.flux, || cx.checkpoint().is_ok())?;
    if let TemperatureSpan::Sampled { low, high } = span {
        if range.lo < low || range.hi > high {
            return Err(unsupported("slab barrier does not establish the original material temperature validity"));
        }
    }
    poll(cx)?;
    let region = bound_field(cx, extended, &admitted, &interfaces,
        temperature, &selected, config)?;
    poll(cx)?;
    Ok(SlabValidatedRegionMean {
        region, continuum_temperature_range_k: [range.lo, range.hi],
        material_temperature_span: span, axis, reverse,
    })
}
