//! Actual native quadrature work, kept outside the caloric heat account.

use super::*;
use fs_flux::MaterialStokesSolution;
use fs_flux::bdm::{edge_gauss_pub, eval_basis, tri_quad};
use fs_math::dd::Dd;

/// Signed original-interval exchanges in joules for the declared extrusion.
/// This is not a complete mechanical-energy or passivity certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminarMechanicalWork {
    /// Integral of additional force density dot accepted BDM velocity.
    /// Boussinesq force has its own existing receipt and is excluded here.
    pub body_work_j: f64,
    /// Integral of supplied open traction dot accepted fluid velocity.
    /// Excludes velocity constraints and any added backflow resistance.
    pub imposed_traction_work_j: f64,
    /// Original upwind carried kinetic energy, signed positive outward.
    pub outward_kinetic_energy_j: f64,
    /// Original numerical boundary upwind loss; never material heating.
    pub boundary_upwind_loss_j: f64,
    /// None when velocity or Navier-slip conditions are present. The original
    /// public reaction owner exposes rigid virtual work, not generalized work
    /// for arbitrary spatial velocity data; no average-velocity proxy is used.
    pub velocity_boundary_work_j: Option<f64>,
    /// None when the separately selected dissipative backflow law is present.
    /// Its stress work is not part of the supplied-traction channel above.
    pub backflow_resistance_work_j: Option<f64>,
}

pub(super) fn measure<D: LaminarMechanicalDrive, E>(
    cx: &Cx<'_>,
    system: &FluxSystem<'_>,
    flow: &MaterialStokesSolution,
    drive: &D,
    depth: f64,
    duration: f64,
) -> Result<LaminarMechanicalWork, LaminarEnthalpyError<E>> {
    let velocity = flow.velocity_dofs();
    let mut body = Dd::from_f64(0.0);
    for (cell, vertices) in system.mesh.tris.iter().enumerate() {
        checkpoint(cx)?;
        let points = vertices.map(|vertex| system.mesh.verts[vertex]);
        for (point, weight) in tri_quad(points, system.mesh.areas[cell]) {
            let force = drive.force_density(point);
            let u = value_at(system, velocity, cell, point);
            body = body + dot(force, u)? * Dd::from_f64(weight);
        }
    }
    let mut imposed = Dd::from_f64(0.0);
    let mut velocity_boundary_work_j = Some(0.0);
    let mut backflow_resistance_work_j = Some(0.0);
    for boundary in drive.boundaries() {
        checkpoint(cx)?;
        let edge = boundary.edge();
        match *boundary {
            LaminarDrivenBoundary::Velocity { .. } | LaminarDrivenBoundary::NavierSlip { .. } => {
                velocity_boundary_work_j = None;
                continue;
            }
            LaminarDrivenBoundary::Open {
                backflow: OpenBoundaryPolicy::DissipativeBackflow,
                ..
            } => backflow_resistance_work_j = None,
            _ => {}
        }
        let geometry = &system.mesh.edges[edge];
        let a = system.mesh.verts[geometry.verts.0];
        let b = system.mesh.verts[geometry.verts.1];
        for (point, weight) in edge_gauss_pub(a, b) {
            let traction = drive.traction(edge, point, geometry.normal);
            let u = value_at(system, velocity, geometry.tris.0, point);
            imposed = imposed + dot(traction, u)? * Dd::from_f64(weight);
        }
    }
    let edges: Vec<_> = (0..system.mesh.edges.len())
        .filter(|&edge| system.mesh.is_boundary(edge))
        .collect();
    let transported = flow
        .boundary_transport(&edges, [0.0; 2], [0.0; 2])
        .map_err(LaminarEnthalpyError::Momentum)?;
    let volume_time = Dd::from_f64(depth) * Dd::from_f64(duration);
    let values = [
        (body * volume_time).to_f64(),
        (imposed * volume_time).to_f64(),
        (Dd::from_f64(transported.outward_kinetic_energy_w_per_m) * volume_time).to_f64(),
        (Dd::from_f64(transported.upwind_kinetic_dissipation_w_per_m) * volume_time).to_f64(),
    ];
    if values.iter().any(|value| !value.is_finite()) {
        return Err(DepositionError::Arithmetic.into());
    }
    checkpoint(cx)?;
    Ok(LaminarMechanicalWork {
        body_work_j: values[0],
        imposed_traction_work_j: values[1],
        outward_kinetic_energy_j: values[2],
        boundary_upwind_loss_j: values[3],
        velocity_boundary_work_j,
        backflow_resistance_work_j,
    })
}

fn value_at(system: &FluxSystem<'_>, dofs: &[f64], cell: usize, point: [f64; 2]) -> [Dd; 2] {
    let mut result = [Dd::from_f64(0.0); 2];
    for local in 0..6 {
        let (edge, _) = system.mesh.tri_edges[cell][local / 2];
        let basis = eval_basis(&system.bases[cell], local, point);
        let coefficient = Dd::from_f64(dofs[2 * edge + local % 2]);
        for axis in 0..2 {
            result[axis] = result[axis] + coefficient * Dd::from_f64(basis[axis]);
        }
    }
    result
}
fn dot<E>(load: [f64; 2], velocity: [Dd; 2]) -> Result<Dd, LaminarEnthalpyError<E>> {
    if load.iter().any(|x| !x.is_finite()) {
        return Err(LaminarEnthalpyError::Invalid(
            "mechanical work source is nonfinite",
        ));
    }
    Ok(Dd::from_f64(load[0]) * velocity[0] + Dd::from_f64(load[1]) * velocity[1])
}
