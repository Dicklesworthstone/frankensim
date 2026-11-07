//! Named coating storage lowered onto each rung's original boundary faces.

use fs_conduction::ConductionProblem;
use fs_conduction::transient::backward_euler::{
    ArealHeatCapacity, BackwardEuler, SurfaceHeatCapacity,
};
use fs_exec::Cx;
use fs_project::TransientSurfaceCapacity;
use std::time::Instant;

use super::{SolveRefusal, bad, json_string, lower, number, poll};

/// Bound preparation and the shared owner's temporary mass assembly before
/// either time grid starts. No extra allocation is charged for an empty list.
pub(super) fn workspace_bytes(
    vertices: usize,
    elements: usize,
    boundary_faces: usize,
    patches: usize,
) -> Result<usize, SolveRefusal> {
    if patches == 0 {
        return Ok(0);
    }
    vertices
        .checked_mul(128)
        .and_then(|bytes| {
            elements
                .checked_mul(1024)
                .and_then(|n| bytes.checked_add(n))
        })
        .and_then(|bytes| {
            boundary_faces
                .checked_mul(128)
                .and_then(|n| bytes.checked_add(n))
        })
        .and_then(|bytes| patches.checked_mul(512).and_then(|n| bytes.checked_add(n)))
        .ok_or_else(|| bad("surface heat-storage workspace estimate overflows"))
}

/// The original boundary partition resolves region/surface precedence once.
/// One linear face scan then prepares the actual geometric storage on this
/// mesh, shared by both time grids and every coupled boundary response.
pub(super) fn prepare<'m>(
    cx: &Cx<'_>,
    engine: BackwardEuler<'m>,
    problem: ConductionProblem<'m>,
    rows: &[TransientSurfaceCapacity],
    deadline: Option<(Instant, f64)>,
) -> Result<(BackwardEuler<'m>, String), SolveRefusal> {
    if rows.is_empty() {
        return Ok((engine, String::new()));
    }
    poll(cx, deadline)?;
    let mut by_region = vec![None; problem.boundary.conditions().len()];
    for (row_index, row) in rows.iter().enumerate() {
        poll(cx, deadline)?;
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|name| name == &row.target)
            .ok_or_else(|| {
                bad(format!(
                    "surface capacity target `{}` is absent from the retained exterior boundary",
                    row.target,
                ))
            })?;
        let capacity = ArealHeatCapacity::declared(row.areal_heat_capacity.value).map_err(lower)?;
        if by_region[region].replace((row_index, capacity)).is_some() {
            return Err(bad(
                "surface capacity targets overlap; combine coincident layers explicitly",
            ));
        }
    }
    let mut faces = Vec::new();
    faces
        .try_reserve_exact(problem.mesh.boundary().len())
        .map_err(|_| bad("surface heat-storage face allocation refused"))?;
    let mut geometry = vec![(0_usize, 0.0_f64); rows.len()];
    for (face_slot, face) in problem.mesh.boundary().iter().enumerate() {
        if face_slot % 512 == 0 {
            poll(cx, deadline)?;
        }
        // An explicitly adiabatic remainder has no named owner and therefore
        // receives no coating from a declaration on another boundary.
        let Some(region) = problem.boundary.region_for(face_slot) else {
            continue;
        };
        let Some((row_index, capacity)) = by_region[region] else {
            continue;
        };
        faces.push(SurfaceHeatCapacity::uniform(face_slot, capacity));
        geometry[row_index].0 += 1;
        geometry[row_index].1 += face.area;
    }
    let mut receipts = Vec::with_capacity(rows.len());
    let mut total = 0.0;
    for (row, &(face_count, area)) in rows.iter().zip(&geometry) {
        if face_count == 0 || !area.is_finite() || area <= 0.0 {
            return Err(bad(format!(
                "surface capacity target `{}` has no finite positive retained exterior area",
                row.target,
            )));
        }
        let capacity = area * row.areal_heat_capacity.value;
        total += capacity;
        receipts.push(format!(
            "{{\"target\":{},\"areal_heat_capacity_j_m2_k\":{},\"source\":{},\"boundary_faces\":{},\"area_m2\":{},\"capacity_j_k\":{}}}",
            json_string(&row.target), number(row.areal_heat_capacity.value)?,
            json_string(&row.source), face_count, number(area)?, number(capacity)?,
        ));
    }
    let total = number(total)?;
    // The existing owner adds the exact P1 row moments to its capacity
    // matrix, including prescribed-node storage and its energy reactions.
    // All-zero coatings return that original matrix unchanged.
    let engine = engine
        .with_surface_heat_capacity(cx, &faces)
        .map_err(lower)?;
    poll(cx, deadline)?;
    let fragment = format!(
        ",\"surface_storage\":{{\"model\":\"perfectly-bonded-temperature-independent-coating\",\"discretization\":\"p1-surface-row-lumped\",\"surfaces\":[{}],\"total_capacity_j_k\":{total},\"scope\":{}}}",
        receipts.join(","),
        json_string(
            "Declared areal capacity shares the original boundary temperature and adds storage to the same volume capacity operator. Both time grids and every thermal, contact, convection and patch-radiation response use that operator. Prescribed nodes retain their storage reactions. The coating adds no separate temperature, contact resistance, tangential conduction or material-source authority. It must not also be represented by volume storage. Coating trajectory adjoints and total-enthalpy storage are not admitted."
        ),
    );
    Ok((engine, fragment))
}
