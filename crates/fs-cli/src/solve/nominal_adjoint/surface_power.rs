//! Pull back declared watts through the actual native chip-footprint load.
//!
//! Surface power lowers to outward flux -P*duty/A_patch. Its derivative is
//! therefore +duty/A_patch times the consistent P1 face load contracted with
//! the COMPLETE thermal dual. No regional volume, nodal source mixing or
//! replacement geometry belongs in this path. Prescribed-node duals are zero.

use super::{bad, finite, poll, Cx, RungSolved, ScalarField, SolveRefusal, ThermalBc};

pub(super) fn pullback(
    cx: &Cx<'_>, solved: &RungSolved, target: &str, duty: f64, lambda: &[f64],
) -> Result<f64, SolveRefusal> {
    poll(cx)?;
    let mesh = &solved.mesh;
    if lambda.len() != mesh.vertex_count() || !duty.is_finite() || !(0.0..=1.0).contains(&duty) {
        return Err(bad("surface-power adjoint needs a full native dual and duty in [0,1]"));
    }
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("missing surface-power boundary operator"))?;
    let slots = solved.surface_slots.get(target)
        .ok_or_else(|| bad("surface-power target has no retained native boundary slots"))?;
    let region = data.boundary.region_names().iter().position(|name| name == target)
        .ok_or_else(|| bad("surface-power target has no lowered boundary owner"))?;
    if !matches!(&data.boundary.conditions()[region],
        ThermalBc::Neumann { outward_flux: ScalarField::Uniform(q) } if q.is_finite() && *q <= 0.0) {
        return Err(bad("surface-power target must retain its uniform inward Neumann load"));
    }
    let mut count = 0;
    let mut area = 0.0;
    let mut load_bar = 0.0;
    for (slot, face) in mesh.boundary().iter().enumerate() {
        if slot % 512 == 0 { poll(cx)?; }
        if data.boundary.region_for(slot) != Some(region) { continue; }
        // The retained producer lists slots in mesh order. Require the exact
        // complete trace, not an unweighted node set, subset or duplicate face.
        if slots.get(count) != Some(&slot) {
            return Err(bad("surface-power slots differ from the retained Neumann trace"));
        }
        count += 1;
        // Share the primal normalization and the assembled P1 load's native
        // physical face measure; no second geometric evaluation may drift.
        area = finite(area + face.area)?;
        for &v in &face.vertices {
            load_bar = finite(load_bar + (face.area/3.0)*finite(lambda[v as usize])?)?;
        }
    }
    if count == 0 || count != slots.len() || area <= 0.0 {
        return Err(bad("surface-power adjoint requires one complete positive-area native trace"));
    }
    poll(cx)?;
    finite((load_bar/area)*duty)
}
