//! Inertial right-wall motion driven by the original evolving species field.

use super::evolve::totals;
use super::moving::{attempt, geometry_error, reported_step};
use super::{
    Boundaries, Boundary, Error, MixtureEuler1d, Reconstruction, StepControl, StepReport, poll,
    positive,
};
use crate::dynamics::MovingWorkspace;
use crate::dynamics::piston::{
    MechanicalExchange, Piston, PistonControl, PistonError, PistonLoad, PistonStep,
    mechanical_report, solve, validate,
};
use fs_qty::Time;

mod two_chambers;
pub use two_chambers::SharedPistonGasReport;

impl MixtureEuler1d {
    /// Advance evolving species, resolved pressure waves and one inertial piston.
    ///
    /// This uses the SAME [`Piston`], controls, bounded midpoint-speed root and
    /// independent momentum/work acceptance as [`crate::dynamics::Euler1d::step_with_piston`].
    /// The original mixture moving-volume operator supplies every trial. Each
    /// species advances extensively with its own formation reference, while
    /// the actual stage composition and caloric source determine wall pressure.
    /// All species fluxes through either wall are exactly zero; absent species
    /// are never seeded and no pressure or composition correction is applied.
    ///
    /// The result retains the native [`StepReport`], including every canonical
    /// species balance. `PistonError::MixtureGas` retains original mixture
    /// failures, while cancellation from source, geometry or scalar iteration
    /// returns `PistonError::Cancelled` before any retry or publication.
    ///
    /// The left wall is fixed; the right wall and interior affine geometry
    /// move at the solved midpoint velocity. Every trial at a given duration
    /// starts from the accepted field. A retry halves the complete coupled
    /// attempt. Only the final admitted trial publishes all species, geometry,
    /// clock and piston velocity. Source/RK workspace is reused, with geometry
    /// workspace allocated fallibly at first use.
    ///
    /// Requires a constant-area closed interval. This does not advance the
    /// separately bound chemistry, diffusion, finite-wall heat or network
    /// owners. Unequal-caloric advecting material contacts retain the original
    /// conservative pressure-oscillation limitation; no exact interface
    /// equilibrium, entropy inequality or general FSI stability is claimed.
    ///
    /// # Errors
    /// Refuses invalid loads/controls, unsupported geometry/source states,
    /// nonlinear or balance failure, exhausted retries and cancellation without
    /// changing accepted gas, mesh, time or piston state.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_piston(
        &mut self,
        piston: &mut Piston,
        maximum_duration: Time,
        gas_controls: StepControl,
        load: PistonLoad,
        controls: PistonControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep<StepReport>, PistonError> {
        poll(&mut cancelled)?;
        validate(*piston, maximum_duration, gas_controls, load, controls)?;
        if self.boundaries
            != (Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            })
            || self.mesh.uniform_cross_section().is_none()
        {
            return Err(PistonError::Invalid(
                "mixture piston requires a constant-area closed interval",
            ));
        }
        if self.moving.is_none() {
            self.moving = Some(MovingWorkspace::new(&self.mesh).map_err(geometry_error)?);
        }
        let moving = self
            .moving
            .as_mut()
            .expect("moving mixture workspace initialized");
        let area = self.mesh.face_areas()[0];
        let force = load.external_force_n - area * load.exterior_pressure_pa;
        if !force.is_finite() {
            return Err(PistonError::Invalid("external piston force overflow"));
        }
        let [lo, hi] = controls.midpoint_velocity_bracket_m_s;
        let n = self.gas.species_count();
        let m = n + 2;
        let count = self.mesh.cell_count();
        let mut rate = 0.0_f64;
        for (i, q) in self.states.chunks_exact(m).enumerate() {
            let p = self
                .gas
                .decode(q, Some(self.work.primitive[i].temperature), &mut cancelled)?;
            self.work.primitive[i] = p;
            let relative = p
                .velocity
                .abs()
                .max((p.velocity - lo).abs())
                .max((p.velocity - hi).abs());
            let local = (relative + p.sound_speed) * area / self.mesh.volumes()[i];
            if !positive(local) {
                return Err(PistonError::Invalid(
                    "mixture piston relative-wave CFL overflow",
                ));
            }
            rate = rate.max(local);
        }
        let mut dt = maximum_duration.value().min(0.9 * gas_controls.cfl / rate);
        if lo < 0.0 {
            let length = self.mesh.positions()[count] - self.mesh.positions()[0];
            dt = dt.min(0.9 * gas_controls.cfl * (length / -lo));
        }
        let before = totals(&self.states, self.mesh.volumes(), n, &mut cancelled)?;
        let mut last = PistonError::Invalid("no admissible coupled mixture piston increment");
        let mut first_order_retry = false;
        let mut gas_trials = 0;
        for halvings in 0..=gas_controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt)
                || !positive(self.time + dt - self.time)
                || !(self.time + dt).is_finite()
            {
                return Err(PistonError::Invalid(
                    "mixture piston step cannot advance the clock",
                ));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                gas_controls.reconstruction
            };
            let result = solve(piston.velocity_m_s(), controls, |speed| {
                gas_trials += 1;
                moving
                    .prepare_velocities(&self.mesh, [0.0, speed], &mut cancelled)
                    .map_err(geometry_error)?;
                let (a, b) = attempt(
                    &self.gas,
                    &self.states,
                    &self.mesh,
                    self.boundaries,
                    dt,
                    gas_controls.cfl,
                    reconstruction,
                    &mut self.work,
                    moving,
                    &mut cancelled,
                )?;
                let after = totals(
                    &self.work.candidate,
                    moving.mesh.volumes(),
                    n,
                    &mut cancelled,
                )?;
                let gas = reported_step(before, after, [a, b], dt, area, n)?;
                let exchange = MechanicalExchange::from(&gas);
                let displacement = moving.mesh.positions()[count] - self.mesh.positions()[count];
                mechanical_report(gas, exchange, *piston, load, speed, displacement, area)
            });
            match result {
                Ok(mut report) => {
                    report.gas.halvings = halvings;
                    report.gas.first_order_retry =
                        first_order_retry && gas_controls.reconstruction == Reconstruction::Muscl;
                    report.gas_trials = gas_trials;
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.states, &mut self.work.candidate);
                    core::mem::swap(&mut self.mesh, &mut moving.mesh);
                    self.time += dt;
                    piston.accept(&report);
                    return Ok(report);
                }
                Err(PistonError::Cancelled) => return Err(PistonError::Cancelled),
                Err(PistonError::MixtureGas(Error::Allocation)) => {
                    return Err(PistonError::MixtureGas(Error::Allocation));
                }
                Err(error) => {
                    if matches!(&error, PistonError::MixtureGas(cause) if !matches!(cause, Error::CflExceeded { .. }))
                    {
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

#[cfg(test)]
mod tests;
