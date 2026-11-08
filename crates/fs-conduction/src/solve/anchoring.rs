//! A steady thermal datum belongs to each component of the actual operator.
//!
//! Mesh connectivity alone is insufficient: finite contact resistance joins
//! independent meshes, while an unrelated body's reservoir cannot determine
//! an insulated body's absolute temperature. This check uses the assembled
//! bulk, sheet and contact graph before eliminating prescribed coordinates.

use fs_exec::Cx;
use fs_sparse::Csr;

use super::{ConductionError, ConductionProblem, DofMap, ThermalBc};

/// Reject an undetermined thermal component before residual-based acceptance.
/// Storage/initial data close isolated transient problems, so the assembler
/// and transient owners deliberately do not call this steady-only gate.
pub(crate) fn require_anchored_components(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    dofs: &DofMap,
    operator: &Csr,
) -> Result<(), ConductionError> {
    let n = operator.nrows();
    let mut parents: Vec<_> = (0..n).collect();
    for row in 0..n {
        poll(cx, row)?;
        let (columns, values) = operator.row(row);
        for (&column, &value) in columns.iter().zip(values) {
            crate::require_finite("steady conduction operator", value)?;
            if column != row && value != 0.0 {
                let a = root(&mut parents, row);
                let b = root(&mut parents, column);
                parents[a.max(b)] = a.min(b);
            }
        }
    }

    let mut anchored = vec![false; n];
    for &vertex in dofs.fixed() {
        poll(cx, vertex)?;
        anchored[root(&mut parents, vertex)] = true;
    }
    for (slot, face) in problem.mesh.boundary().iter().enumerate() {
        poll(cx, slot)?;
        let Some(ThermalBc::Robin { htc, .. }) = problem.boundary.condition_for(slot) else {
            continue;
        };
        let vertices = face.vertices.map(|vertex| vertex as usize);
        // An enum tag or a positive coefficient elsewhere in a nodal field
        // is not a local reservoir. Use the same integrated face block as
        // assembly, including underflow to an unrepresented zero coupling.
        let block = crate::assemble::robin_face_matrix(face.area, vertices.map(|v| htc.at(v)));
        for (local, &vertex) in vertices.iter().enumerate() {
            let diagonal = crate::require_finite("Robin anchor", block[local][local])?;
            if diagonal > 0.0 {
                anchored[root(&mut parents, vertex)] = true;
            }
        }
    }
    for vertex in 0..n {
        poll(cx, vertex)?;
        if !anchored[root(&mut parents, vertex)] {
            return Err(ConductionError::SingularPureNeumann);
        }
    }
    Ok(())
}

fn root(parents: &mut [usize], mut vertex: usize) -> usize {
    while parents[vertex] != vertex {
        parents[vertex] = parents[parents[vertex]];
        vertex = parents[vertex];
    }
    vertex
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    if at % crate::ASSEMBLY_TILE == 0 {
        cx.checkpoint().map_err(|_| ConductionError::Cancelled {
            stage: "steady-thermal-anchoring",
            at,
        })?;
    }
    Ok(())
}
