//! Joint native reservoir reconstruction and implicit porous forcing.

use super::*;
use crate::d3q19::boundary::partial_saturation::open::OpenTransport;
use fs_la::prepared_lu::LuWorkspace;

pub(super) fn transport_error(error: PartialSaturationError3) -> BoundaryStepError3 {
    match error {
        PartialSaturationError3::Cancelled => BoundaryStepError3::Cancelled,
        PartialSaturationError3::Invalid { reason } => BoundaryStepError3::EndpointForce { reason },
        PartialSaturationError3::Cell { cell, .. } => BoundaryStepError3::Unphysical { cell },
        PartialSaturationError3::Collision { cell, source } => {
            BoundaryStepError3::Collision { cell, source }
        }
        PartialSaturationError3::Motion { .. } => unreachable!("transport has no body motion"),
    }
}

impl BoundaryGrid3 {
    pub(super) fn reconstruct_resistance_open_faces(
        &mut self,
        next: &mut [Vec<Tile>; Q3],
        forces: &mut [[f64; 3]],
        prepared: &[PreparedResistance3],
        exchange: &mut ResistanceMomentumExchange3,
        origin: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), BoundaryStepError3> {
        if self.open_cells.is_none() {
            self.open_cells = Some(self.compile_open_cells());
        }
        let cells = self.open_cells.take().unwrap_or_default();
        let outcome = (|| {
            let mut transport = OpenTransport::default();
            // Only mixed wall/open cells need a factorization. Reuse its nine
            // factors across every iteration/cell; face interiors solve directly.
            let mut workspace = None;
            for (index, &(cell, interior, face)) in cells.iter().enumerate() {
                if index.is_multiple_of(TILE_CELLS) {
                    endpoint_checkpoint(gate)?;
                }
                let (tile, lane) = self.addr(cell[0], cell[1], cell[2]);
                let flat = tile * TILE_CELLS + lane;
                let before = core::array::from_fn(|q| next[q][tile].0[lane]);
                let post = core::array::from_fn(|q| self.post[q][tile].0[lane]);
                let open_mask = self.open_link_masks[tile].0[lane];
                transport.add_links(&post, open_mask, cell, origin);
                let (it, il) = self.addr(interior[0], interior[1], interior[2]);
                let (neighbor_rho, neighbor_velocity, stress) = self
                    .regularized_source_with_force(next, interior, forces[it * TILE_CELLS + il])?;
                let (rho, velocity) = self.open_target(cell, face, neighbor_rho, neighbor_velocity);
                let source = OpenSource {
                    template: regularized_populations(rho, velocity, stress),
                    velocity,
                    before,
                    wall_mask: self.link_masks[tile].0[lane],
                    open_mask,
                    law: prepared[flat].law,
                };
                let accepted = source.solve(&mut workspace, gate).ok_or_else(|| {
                    if gate.is_requested() {
                        BoundaryStepError3::Cancelled
                    } else {
                        BoundaryStepError3::Unphysical { cell }
                    }
                })?;
                let contribution = prepared[flat]
                    .exchange(accepted.velocity, accepted.drag, cell, origin)
                    .ok_or(BoundaryStepError3::Unphysical { cell })?;
                transport.add_reconstruction(before, accepted.populations, cell, origin);
                for (field, value) in next.iter_mut().zip(accepted.populations) {
                    field[tile].0[lane] = value;
                }
                forces[flat] = accepted.force;
                exchange.add(contribution);
            }
            exchange.open = transport.finish().map_err(transport_error)?;
            if !exchange.is_finite() {
                return Err(BoundaryStepError3::EndpointForce {
                    reason: "resistance aggregate impulse or work is not representable",
                });
            }
            endpoint_checkpoint(gate)
        })();
        self.open_cells = Some(cells);
        outcome
    }
}

struct OpenSource {
    template: [f64; Q3],
    velocity: [f64; 3],
    before: [f64; Q3],
    wall_mask: u32,
    open_mask: u32,
    law: DarcyForchheimer3,
}

struct Evaluation {
    populations: [f64; Q3],
    rho: f64,
    velocity: [f64; 3],
    force: [f64; 3],
    drag: [f64; 3],
    gradient: [[f64; 3]; 3],
    residual: [f64; 3],
    scaled_error: f64,
}

impl OpenSource {
    fn evaluate(&self, force: [f64; 3]) -> Option<Evaluation> {
        if !force.iter().all(|f| f.is_finite()) {
            return None;
        }
        let populations = reconstruct_open_populations(
            self.template,
            self.velocity,
            force,
            self.before,
            self.wall_mask,
            self.open_mask,
        );
        let (rho, momentum) = raw_moments(&populations)?;
        let velocity = core::array::from_fn(|a| (momentum[a] + 0.5 * force[a]) / rho);
        let drag = self.law.drag(velocity);
        let constitutive = self.law.force(velocity);
        let residual: [f64; 3] = core::array::from_fn(|a| force[a] - constitutive[a]);
        let gradient = drag_gradient(self.law, velocity);
        if !velocity
            .iter()
            .chain(&drag)
            .chain(&residual)
            .all(|v| v.is_finite())
            || !gradient.iter().flatten().all(|v| v.is_finite())
        {
            return None;
        }
        // Resolve weak drag at its own force scale, while admitting the actual
        // raw-moment roundoff propagated through its derivative. No unit-sized
        // absolute force tolerance can swallow a small constitutive response.
        let density_scale: f64 = populations.iter().map(|f| f.abs()).sum();
        let velocity_scale: [f64; 3] = core::array::from_fn(|a| {
            let moment_scale: f64 = populations
                .iter()
                .zip(E3)
                .map(|(f, c)| f.abs() * f64::from([c.0, c.1, c.2][a]).abs())
                .sum();
            (moment_scale + 0.5 * force[a].abs() + velocity[a].abs() * density_scale) / rho
        });
        let mut scaled_error: f64 = 0.0;
        for a in 0..3 {
            let rounding_scale: f64 = (0..3)
                .map(|b| gradient[a][b].abs() * velocity_scale[b])
                .sum();
            let scale = force[a].abs()
                + self.law.linear.body_force[a].abs()
                + drag[a].abs()
                + rounding_scale;
            if !scale.is_finite() {
                return None;
            }
            let tolerance = 64.0 * f64::EPSILON * scale;
            let error = if residual[a] == 0.0 {
                0.0
            } else {
                residual[a].abs() / tolerance
            };
            scaled_error = scaled_error.max(error);
        }
        Some(Evaluation {
            populations,
            rho,
            velocity,
            force,
            drag,
            gradient,
            residual,
            scaled_error,
        })
    }

    fn solve(&self, workspace: &mut Option<LuWorkspace>, gate: &CancelGate) -> Option<Evaluation> {
        let first = self.evaluate(self.law.force(self.velocity));
        if self.wall_mask == 0 {
            // The complete regularized trace has exactly the prescribed
            // physical velocity, independent of the Guo storage frame.
            return first.filter(|candidate| candidate.scaled_error <= 1.0);
        }
        let mut current = first.or_else(|| self.evaluate([0.0; 3]))?;
        if current.scaled_error <= 1.0 {
            return Some(current);
        }
        if workspace.is_none() {
            *workspace = Some(LuWorkspace::new(3).ok()?);
        }
        let workspace = workspace.as_mut()?;
        // At fixed template velocity, the native Guo transform is linear in
        // force. Differentiate its OWN rows, then retain only open-owned rows.
        // rho(F)=rho0+a.F; physical momentum(F)=j0+C.F. The rest-population
        // mass correction is wall-owned here and is therefore excluded.
        let mut density_derivative = [0.0; 3];
        let mut physical_derivative = [[0.0; 3]; 3];
        for axis in 0..3 {
            let mut unit = [0.0; 3];
            unit[axis] = 1.0;
            let mut derivative = [0.0; Q3];
            crate::thermal::forcing::reframe3(&mut derivative, self.velocity, [0.0; 3], unit);
            physical_derivative[axis][axis] = 0.5;
            for q in 1..Q3 {
                if self.open_mask & (1u32 << q) != 0 {
                    density_derivative[axis] += derivative[q];
                    for row in 0..3 {
                        physical_derivative[row][axis] +=
                            f64::from([E3[q].0, E3[q].1, E3[q].2][row]) * derivative[q];
                    }
                }
            }
        }
        for _ in 0..40 {
            if gate.is_requested() {
                return None;
            }
            // d(F-law(u(F)))/dF = I + K (C-u a^T)/rho.
            let mut matrix = [0.0; 9];
            let mut rhs = current.residual;
            for row in 0..3 {
                for column in 0..3 {
                    matrix[3 * row + column] = f64::from(row == column)
                        + (0..3)
                            .map(|a| {
                                current.gradient[row][a]
                                    * (physical_derivative[a][column]
                                        - current.velocity[a] * density_derivative[column])
                                    / current.rho
                            })
                            .sum::<f64>();
                }
                let scale = matrix[3 * row..3 * row + 3]
                    .iter()
                    .fold(0.0_f64, |m, x| m.max(x.abs()));
                if !scale.is_finite() || scale == 0.0 {
                    return None;
                }
                for coefficient in &mut matrix[3 * row..3 * row + 3] {
                    *coefficient /= scale;
                }
                rhs[row] /= scale;
            }
            let mut correction = [0.0; 3];
            workspace.solve_into(&matrix, &rhs, &mut correction).ok()?;
            let mut alpha = 1.0;
            let mut improved = None;
            for _ in 0..24 {
                let trial = core::array::from_fn(|a| current.force[a] - alpha * correction[a]);
                if let Some(candidate) = self.evaluate(trial) {
                    if candidate.scaled_error <= 1.0 {
                        // Accept only the represented populations and their
                        // force/macros, never an unverified small Newton step.
                        return Some(candidate);
                    }
                    if candidate.scaled_error < current.scaled_error {
                        improved = Some(candidate);
                        break;
                    }
                }
                alpha *= 0.5;
            }
            current = improved?;
        }
        None
    }
}

fn drag_gradient(law: DarcyForchheimer3, velocity: [f64; 3]) -> [[f64; 3]; 3] {
    let linear = law.linear;
    let mut gradient = if linear.principal.iter().all(|&d| d == linear.principal[0]) {
        core::array::from_fn(|a| {
            core::array::from_fn(|b| if a == b { linear.principal[0] } else { 0.0 })
        })
    } else {
        core::array::from_fn(|a| {
            core::array::from_fn(|b| {
                (0..3)
                    .map(|k| linear.principal[k] * linear.axes[k][a] * linear.axes[k][b])
                    .sum()
            })
        })
    };
    if law.quadratic != 0.0 {
        let relative: [f64; 3] =
            core::array::from_fn(|a| velocity[a] - linear.skeleton_velocity[a]);
        let speed = relative[0].hypot(relative[1]).hypot(relative[2]);
        if speed != 0.0 {
            let direction = relative.map(|v| v / speed);
            let coefficient = law.quadratic * speed;
            for a in 0..3 {
                for b in 0..3 {
                    gradient[a][b] +=
                        coefficient * (f64::from(a == b) + direction[a] * direction[b]);
                }
            }
        }
    }
    gradient
}
