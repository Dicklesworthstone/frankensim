//! Two original resolved species fields exchanging work through one piston.

use super::*;
use crate::multicomponent::{Totals, sum};

/// Original extensive reports on either side of one impermeable piston.
///
/// Each side retains its own species order, caloric source and energy datum.
/// The interface is `left.right_boundary_change` and
/// `right.left_boundary_change`; the other two ends are fixed supports.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedPistonGasReport {
    /// Resolved chamber to the left of the piston, with a fixed left end.
    pub left: StepReport,
    /// Resolved chamber to the right of the piston, with a fixed right end.
    pub right: StepReport,
}

impl MixtureEuler1d {
    /// Advance two resolved gas chambers and their shared inertial piston.
    ///
    /// `self` is the left chamber and `right` is the right chamber in one axial
    /// coordinate system. Their interface positions, physical times and constant
    /// cross sections must match exactly. The two sealed chambers may have
    /// different caloric sources, species axes, temperatures and resolutions;
    /// their original source owners are retained and no species crosses the wall.
    /// The piston has negligible axial thickness. Additional force [N] is positive
    /// to the right and excludes BOTH resolved gas pressure forces.
    ///
    /// Every midpoint-speed proposal runs the original extensive moving-volume
    /// RK2 operator on both sides, from their accepted fields at one duration.
    /// Positive speed expands the left volume and compresses the right volume.
    /// Both evolving interface pressures supply the same mechanical equation:
    /// `M*(v1-v0) = -(I_left_gas + I_right_gas) + dt*F_external`.
    /// The same face impulses supply gas work. Actual endpoint gas energies,
    /// piston kinetic change and represented force work close independently.
    /// Outer fixed-wall reactions close combined momentum and carry no work.
    /// `PistonStep::fixed_wall_impulse_n_s` is their sum; the original side
    /// reports retain each separate reaction. Exterior-pressure work is zero.
    ///
    /// Initial acoustic/mesh-contraction limits cover both chambers and the
    /// complete speed bracket. A failed stage or root halves and restarts BOTH
    /// fields. The original safeguarded scalar root and each field's reusable
    /// candidate buffers own every proposal; no full field or source table is
    /// cloned per trial. Both fields, meshes, clocks and piston publish together
    /// after the final cancellation poll. Rejected scratch remains reusable.
    /// `gas_trials` counts complete coupled proposals, each with two field solves.
    ///
    /// This is sealed adiabatic multicomponent Euler transport. Separately bound
    /// chemistry, molecular transport, finite thermal walls, gravity and leakage
    /// are outside this entry. The original pressure-oscillation limitation at
    /// advecting unequal-caloric contacts remains; the shared root establishes
    /// neither a general FSI stability theorem nor a multidimensional flow rung.
    ///
    /// # Errors
    /// Invalid coupling geometry, clocks, controls, source states, failed balances,
    /// exhausted retries or cancellation leave all accepted physical state intact.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_piston(
        &mut self,
        right: &mut MixtureEuler1d,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        external_force_n: f64,
        controls: PistonControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep<SharedPistonGasReport>, PistonError> {
        poll(&mut cancelled)?;
        let load = PistonLoad {
            exterior_pressure_pa: 0.0,
            external_force_n,
        };
        validate(*piston, maximum_duration, gas_controls, load, controls)?;
        let closed = Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        };
        let area = self.mesh.uniform_cross_section().ok_or(PistonError::Invalid(
            "shared piston requires constant-area gas chambers",
        ))?;
        let interface = self.mesh.positions()[self.mesh.cell_count()];
        if self.boundaries != closed
            || right.boundaries != closed
            || right.mesh.uniform_cross_section() != Some(area)
            || right.mesh.positions()[0] != interface
            || self.time != right.time
        {
            return Err(PistonError::Invalid(
                "shared piston requires closed chambers with equal areas, interface coordinates and clocks",
            ));
        }
        let original_piston = *piston;
        let [lo, hi] = controls.midpoint_velocity_bracket_m_s;
        let (left_before, left_rate) = prepare(self, [lo, hi], &mut cancelled)?;
        let (right_before, right_rate) = prepare(right, [lo, hi], &mut cancelled)?;
        let mut dt = maximum_duration
            .value()
            .min(0.9 * gas_controls.cfl / left_rate.max(right_rate));
        if lo < 0.0 {
            let length = interface - self.mesh.positions()[0];
            dt = dt.min(0.9 * gas_controls.cfl * (length / -lo));
        }
        if hi > 0.0 {
            let length = right.mesh.positions()[right.mesh.cell_count()] - interface;
            dt = dt.min(0.9 * gas_controls.cfl * (length / hi));
        }
        let mut last = PistonError::Invalid("no admissible shared-piston increment");
        let mut first_order_retry = false;
        let mut gas_trials = 0;
        for halvings in 0..=gas_controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt)
                || !positive(self.time + dt - self.time)
                || !(self.time + dt).is_finite()
            {
                return Err(PistonError::Invalid(
                    "shared-piston duration cannot advance both clocks",
                ));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                gas_controls.reconstruction
            };
            let result = solve(original_piston.velocity_m_s(), controls, |speed| {
                gas_trials += 1;
                let left = trial(
                    self,
                    left_before,
                    dt,
                    gas_controls.cfl,
                    reconstruction,
                    [0.0, speed],
                    &mut cancelled,
                )?;
                let right_report = trial(
                    right,
                    right_before,
                    dt,
                    gas_controls.cfl,
                    reconstruction,
                    [speed, 0.0],
                    &mut cancelled,
                )?;
                let left_mesh = &self.moving.as_ref().expect("prepared left motion").mesh;
                let right_mesh = &right.moving.as_ref().expect("prepared right motion").mesh;
                let next_interface = left_mesh.positions()[left_mesh.cell_count()];
                if next_interface != right_mesh.positions()[0] {
                    return Err(PistonError::Invalid(
                        "shared piston candidate endpoints do not represent one wall",
                    ));
                }
                let exchange = MechanicalExchange {
                    duration_s: dt,
                    momentum_change: sum([
                        left.after.momentum - left.before.momentum,
                        right_report.after.momentum - right_report.before.momentum,
                    ]),
                    energy_change: sum([
                        left.after.energy - left.before.energy,
                        right_report.after.energy - right_report.before.energy,
                    ]),
                    left_impulse: sum([
                        left.left_boundary_change.momentum,
                        right_report.right_boundary_change.momentum,
                    ]),
                    right_impulse: sum([
                        left.right_boundary_change.momentum,
                        right_report.left_boundary_change.momentum,
                    ]),
                    right_work: sum([
                        left.right_boundary_change.energy,
                        right_report.left_boundary_change.energy,
                    ]),
                };
                mechanical_report(
                    SharedPistonGasReport {
                        left,
                        right: right_report,
                    },
                    exchange,
                    original_piston,
                    load,
                    speed,
                    next_interface - interface,
                    area,
                )
            });
            match result {
                Ok(mut report) => {
                    for side in [&mut report.gas.left, &mut report.gas.right] {
                        side.halvings = halvings;
                        side.first_order_retry = first_order_retry
                            && gas_controls.reconstruction == Reconstruction::Muscl;
                    }
                    report.gas_trials = gas_trials;
                    poll(&mut cancelled)?;
                    for flow in [self, right] {
                        let moving = flow.moving.as_mut().expect("prepared accepted motion");
                        core::mem::swap(&mut flow.states, &mut flow.work.candidate);
                        core::mem::swap(&mut flow.mesh, &mut moving.mesh);
                        flow.time += dt;
                    }
                    piston.accept(&report);
                    return Ok(report);
                }
                Err(PistonError::Cancelled) => return Err(PistonError::Cancelled),
                Err(PistonError::MixtureGas(Error::Allocation)) => {
                    return Err(PistonError::MixtureGas(Error::Allocation));
                }
                Err(error) => {
                    if matches!(&error, PistonError::MixtureGas(cause) if !matches!(cause, Error::CflExceeded { .. })) {
                        first_order_retry = true;
                    }
                    last = error;
                    dt *= 0.5;
                }
            }
        }
        Err(PistonError::StepRetries {
            attempts: gas_controls.maximum_halvings + 1,
            last: Box::new(last),
        })
    }
}

fn prepare(
    flow: &mut MixtureEuler1d,
    bracket: [f64; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Totals, f64), PistonError> {
    if flow.moving.is_none() {
        flow.moving = Some(MovingWorkspace::new(&flow.mesh).map_err(geometry_error)?);
    }
    let n = flow.gas.species_count();
    let area = flow.mesh.face_areas()[0];
    let mut rate = 0.0_f64;
    for (i, q) in flow.states.chunks_exact(n + 2).enumerate() {
        let p = flow.gas.decode(q, Some(flow.work.primitive[i].temperature), cancelled)?;
        flow.work.primitive[i] = p;
        let relative = p.velocity.abs()
            .max((p.velocity - bracket[0]).abs())
            .max((p.velocity - bracket[1]).abs());
        let local = (relative + p.sound_speed) * area / flow.mesh.volumes()[i];
        if !positive(local) {
            return Err(PistonError::Invalid("shared-piston relative-wave CFL overflow"));
        }
        rate = rate.max(local);
    }
    Ok((totals(&flow.states, flow.mesh.volumes(), n, cancelled)?, rate))
}

#[allow(clippy::too_many_arguments)]
fn trial(
    flow: &mut MixtureEuler1d,
    before: Totals,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    wall_velocity: [f64; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StepReport, PistonError> {
    let moving = flow.moving.as_mut().expect("prepared shared-piston motion");
    moving.prepare_velocities(&flow.mesh, wall_velocity, cancelled).map_err(geometry_error)?;
    let (a, b) = attempt(
        &flow.gas,
        &flow.states,
        &flow.mesh,
        flow.boundaries,
        dt,
        cfl,
        reconstruction,
        &mut flow.work,
        moving,
        cancelled,
    )?;
    let n = flow.gas.species_count();
    let after = totals(&flow.work.candidate, moving.mesh.volumes(), n, cancelled)?;
    Ok(reported_step(before, after, [a, b], dt, flow.mesh.face_areas()[0], n)?)
}

#[cfg(test)]
mod tests;
