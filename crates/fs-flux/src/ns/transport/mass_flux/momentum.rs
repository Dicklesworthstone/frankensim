//! Conservative momentum convection on the original BDM1 velocity space.
//!
//! The P2 volume reconstruction and original split face flux belong to one
//! operator: replacing either by rho times an unrelated cell velocity loses
//! the integration-by-parts identity when density crosses a face.

use super::{DensityTransportError, PreparedDensityMassFlux, Result, invalid, poll};
use crate::FluxSystem;
use crate::bdm::{eval_basis, tri_quad};
use crate::ns::convection_boundary_trace;
use fs_math::dd::Dd;

/// Reusable accumulation for the uneliminated native velocity rows.
/// A warm application allocates no storage. Scratch may change on refusal;
/// the caller's output is published only after the final cancellation poll.
pub struct DensityMomentumWorkspace {
    sums: Vec<Dd>,
}

impl PreparedDensityMassFlux<'_> {
    /// Allocate compensated sums for exactly two velocity moments per edge.
    /// `max_storage_bytes` bounds retained vector elements, excluding caller
    /// input/output and the already prepared mass-flux and velocity owners.
    pub fn momentum_workspace(&self, max_storage_bytes: usize) -> Result<DensityMomentumWorkspace> {
        let count = self
            .mesh
            .edges
            .len()
            .checked_mul(2)
            .ok_or(invalid("momentum workspace size overflow"))?;
        let required = count
            .checked_mul(size_of::<Dd>())
            .ok_or(invalid("momentum workspace size overflow"))?;
        if required > max_storage_bytes {
            return Err(DensityTransportError::Storage {
                required,
                maximum: max_storage_bytes,
            });
        }
        let mut sums = Vec::new();
        sums.try_reserve_exact(count)
            .map_err(|_| DensityTransportError::Allocation)?;
        sums.resize(count, Dd::ZERO);
        Ok(DensityMomentumWorkspace { sums })
    }

    /// Apply conservative momentum convection to original BDM1 velocity DOFs.
    /// For each test velocity v this returns
    /// `-sum_K integral((u tensor m):grad(v)) + sum_F integral(m.n u_up.[v])`,
    /// per unit physical depth. The mass flux m already includes donor density;
    /// do NOT multiply the result by another density. No timestep is included.
    /// Dotting the returned rows with velocity moments gives power [W/m].
    ///
    /// `inflow` follows original edge order. Every incoming boundary stream
    /// requires two finite velocity samples at that edge's ORIGINAL Gauss
    /// sites, representing one affine exterior trace. Interior entries must
    /// be None. Outflow uses the interior velocity; a reversing boundary keeps
    /// both directions. The velocity being transported need not equal the
    /// frozen mass advector (e.g. a Picard momentum iteration).
    ///
    /// All velocity rows remain uneliminated, including support rows. Pressure,
    /// strong boundary conditions, viscosity, density/time updates and solves
    /// remain with their owners. This entry imposes no mass or energy repair.
    /// In exact arithmetic its power equals half `integral(div(m) |u|^2)` plus
    /// outward upwind kinetic-energy flux and nonnegative upwind jump loss.
    /// Together with `rho_dot + div(m)=0`, this is the conservative momentum
    /// identity; it does not on its own establish time-discrete stability or
    /// deposit the mixing loss into heat.
    ///
    /// Work is linear in cells plus edges, with fixed-size local quadrature.
    /// Cancellation is polled per cell/edge, per 1024 input/accumulator entries
    /// and before publication. Failure preserves output, while workspace
    /// remains reusable for a retry. The final output copy is indivisible.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_momentum_into(
        &self,
        system: &FluxSystem<'_>,
        velocity: &[f64],
        inflow: &[Option<[[f64; 2]; 2]>],
        output: &mut [f64],
        workspace: &mut DensityMomentumWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<()> {
        poll(&mut cancelled)?;
        if !core::ptr::eq(self.mesh, system.mesh)
            || system.n_u != self.mesh.edges.len() * 2
            || system.bases.len() != self.cells.len()
            || velocity.len() != system.n_u
            || output.len() != system.n_u
            || workspace.sums.len() != system.n_u
            || inflow.len() != self.faces.len()
        {
            return Err(invalid(
                "original momentum mesh and velocity dimensions required",
            ));
        }
        for chunk in velocity.chunks(1024) {
            poll(&mut cancelled)?;
            if chunk.iter().any(|value| !value.is_finite()) {
                return Err(invalid("finite transported velocity moments required"));
            }
        }
        for (edge_id, edge) in self.mesh.edges.iter().enumerate() {
            poll(&mut cancelled)?;
            let samples = inflow[edge_id];
            if let Some(samples) = samples {
                if edge.tris.1 != usize::MAX
                    || samples.iter().flatten().any(|value| !value.is_finite())
                {
                    return Err(invalid(
                        "finite exterior velocity belongs only to boundary edges",
                    ));
                }
            }
            if edge.tris.1 == usize::MAX
                && self.faces[edge_id].directional[1] > 0.0
                && samples.is_none()
            {
                return Err(invalid(
                    "every incoming mass stream needs an exterior velocity",
                ));
            }
        }
        for chunk in workspace.sums.chunks_mut(1024) {
            poll(&mut cancelled)?;
            chunk.fill(Dd::ZERO);
        }
        for cell in 0..self.cells.len() {
            poll(&mut cancelled)?;
            let basis = &system.bases[cell];
            let p = self.mesh.tris[cell].map(|vertex| self.mesh.verts[vertex]);
            for (point, weight) in tri_quad(p, self.mesh.areas[cell]) {
                let mass = self.cell_value(cell, point)?;
                let u = cell_velocity(system, velocity, cell, point)?;
                for local in 0..6 {
                    let gradient = basis.grad[local];
                    let derivative = [
                        mass[0] * gradient[0][0] + mass[1] * gradient[0][1],
                        mass[0] * gradient[1][0] + mass[1] * gradient[1][1],
                    ];
                    add(
                        &mut workspace.sums[dof(system, cell, local)],
                        -weight * (u[0] * derivative[0] + u[1] * derivative[1]),
                    )?;
                }
            }
        }
        for (edge_id, edge) in self.mesh.edges.iter().enumerate() {
            poll(&mut cancelled)?;
            let face = &self.faces[edge_id];
            for point in &face.points[..face.count] {
                let rate = point.mass_rate_kg_m_s;
                if rate == 0.0 {
                    continue;
                }
                let upwind = if rate > 0.0 {
                    cell_velocity(system, velocity, edge.tris.0, point.point_m)?
                } else if edge.tris.1 != usize::MAX {
                    cell_velocity(system, velocity, edge.tris.1, point.point_m)?
                } else {
                    convection_boundary_trace(
                        inflow[edge_id].expect("incoming exterior velocity admitted"),
                        point.boundary_fraction,
                    )
                };
                for (cell, sign) in [(edge.tris.0, 1.0), (edge.tris.1, -1.0)] {
                    if cell == usize::MAX {
                        continue;
                    }
                    for local in 0..6 {
                        let test = eval_basis(&system.bases[cell], local, point.point_m);
                        add(
                            &mut workspace.sums[dof(system, cell, local)],
                            sign * rate * (upwind[0] * test[0] + upwind[1] * test[1]),
                        )?;
                    }
                }
            }
        }
        for chunk in workspace.sums.chunks(1024) {
            poll(&mut cancelled)?;
            if chunk.iter().any(|sum| !sum.to_f64().is_finite()) {
                return Err(DensityTransportError::Arithmetic(
                    "momentum transport accumulation",
                ));
            }
        }
        poll(&mut cancelled)?;
        for (out, sum) in output.iter_mut().zip(&workspace.sums) {
            *out = sum.to_f64();
        }
        Ok(())
    }
}

fn dof(system: &FluxSystem<'_>, cell: usize, local: usize) -> usize {
    2 * system.mesh.tri_edges[cell][local / 2].0 + local % 2
}

fn cell_velocity(
    system: &FluxSystem<'_>,
    velocity: &[f64],
    cell: usize,
    point: [f64; 2],
) -> Result<[f64; 2]> {
    let mut value = [0.0; 2];
    for local in 0..6 {
        let basis = eval_basis(&system.bases[cell], local, point);
        let coefficient = velocity[dof(system, cell, local)];
        value[0] += coefficient * basis[0];
        value[1] += coefficient * basis[1];
    }
    if value.iter().any(|x| !x.is_finite()) {
        return Err(DensityTransportError::Arithmetic(
            "transported momentum velocity",
        ));
    }
    Ok(value)
}

fn add(sum: &mut Dd, value: f64) -> Result<()> {
    if !value.is_finite() {
        return Err(DensityTransportError::Arithmetic(
            "momentum transport quadrature",
        ));
    }
    *sum = *sum + Dd::from_f64(value);
    Ok(())
}
