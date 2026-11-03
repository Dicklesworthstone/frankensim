//! Two finite end walls driven by the original resolved species field.
//!
//! Species, caloric sources and energy references stay on the original ALE
//! owner. Only its actual extensive momentum/work exchanges enter the SAME
//! simultaneous wall solve and physical acceptance as the frozen-gas route.

use super::evolve::totals;
use super::moving::{attempt, geometry_error, reported_step};
use super::{
    Boundaries, Boundary, Error, MixtureEuler1d, Reconstruction, StepControl, StepReport, poll,
    positive,
};
use crate::dynamics::MovingWorkspace;
use crate::dynamics::inertial_walls::{
    InertialWallLoad, InertialWallsControl, InertialWallsStep, MechanicalExchange,
    mechanical_report, solve, validate_controls,
};
use crate::dynamics::piston::{Piston, PistonError};
use fs_qty::Time;

impl MixtureEuler1d {
    /// Advance every resolved species and BOTH freely moving finite end walls.
    ///
    /// Wall order is [left, right]; both speeds and additional forces use the
    /// increasing axial coordinate. Exterior pressure pushes the left wall
    /// rightwards and the right wall leftwards. The shared
    /// [`InertialWallsControl`] and [`InertialWallsStep`] retain exactly the
    /// frozen owner's simultaneous midpoint-speed solve, momentum/work checks
    /// and independent combined gas-plus-two-wall energy acceptance.
    ///
    /// Every trial advances the original species partial densities, momentum
    /// and referenced total energy with the source-dependent ALE RK2 operator.
    /// Both walls are impermeable: each species boundary flux is zero. Local
    /// composition and caloric sources determine each stage's pressure; they
    /// are never frozen to a nominal mixture or shared constant gamma. The
    /// complete native report retains every canonical species balance and
    /// independent formation reference. Absent species are never seeded.
    ///
    /// All proposals at one duration start from the same accepted field. The
    /// existing contiguous species buffers and source workspace are reused;
    /// there is no full-field/source clone per trial. All species, geometry,
    /// time and both wall velocities publish together after the final gate.
    /// Both speed intervals bound the initial relative-wave/contraction CFL;
    /// native stage CFL and source refusal retry the entire coupled interval.
    ///
    /// This is sealed adiabatic Euler transport. Chemistry, molecular diffusion,
    /// viscosity, wall heat, leakage and transverse motion remain separate.
    /// The original pressure-oscillation limitation at advecting unequal-caloric
    /// contacts remains; neither exact interface equilibrium nor a general FSI
    /// stability or entropy theorem follows from the shared mechanical solve.
    ///
    /// # Errors
    /// Invalid coupling inputs, unsupported geometry/source states, failed
    /// physical/root admission, exhausted retries or cancellation leave every
    /// accepted owner unchanged. Original mixture failures retain their causes.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_inertial_walls(
        &mut self,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        gas_controls: StepControl,
        loads: [InertialWallLoad; 2],
        controls: InertialWallsControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<InertialWallsStep<StepReport>, PistonError> {
        poll(&mut cancelled)?;
        validate_controls(*walls, maximum_duration, gas_controls, loads, controls)?;
        if self.boundaries
            != (Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            })
            || self.mesh.uniform_cross_section().is_none()
        {
            return Err(PistonError::Invalid(
                "two inertial walls require a constant-area closed mixture interval",
            ));
        }
        let area = self.mesh.face_areas()[0];
        let force = [
            loads[0].external_force_n + area * loads[0].exterior_pressure_pa,
            loads[1].external_force_n - area * loads[1].exterior_pressure_pa,
        ];
        if force.iter().any(|value| !value.is_finite()) {
            return Err(PistonError::Invalid("external wall force overflow"));
        }
        let bounds = controls.midpoint_velocity_bracket_m_s;
        let low = bounds[0][0].min(bounds[1][0]);
        let high = bounds[0][1].max(bounds[1][1]);
        let species = self.gas.species_count();
        let components = species + 2;
        let count = self.mesh.cell_count();
        let mut rate = 0.0_f64;
        for (i, state) in self.states.chunks_exact(components).enumerate() {
            let p = self.gas.decode(
                state,
                Some(self.work.primitive[i].temperature),
                &mut cancelled,
            )?;
            self.work.primitive[i] = p;
            let relative = (p.velocity - low).abs().max((p.velocity - high).abs());
            let local = (relative + p.sound_speed) * (area / self.mesh.volumes()[i]);
            if !positive(local) {
                return Err(PistonError::Invalid(
                    "two-wall mixture relative-wave CFL overflow",
                ));
            }
            rate = rate.max(local);
        }
        let mut dt = maximum_duration.value().min(0.9 * gas_controls.cfl / rate);
        let contraction = bounds[0][1] - bounds[1][0];
        if !contraction.is_finite() {
            return Err(PistonError::Invalid("two-wall contraction bound overflow"));
        }
        if contraction > 0.0 {
            let length = self.mesh.positions()[count] - self.mesh.positions()[0];
            dt = dt.min(0.9 * gas_controls.cfl * (length / contraction));
        }
        if self.moving.is_none() {
            self.moving = Some(MovingWorkspace::new(&self.mesh).map_err(geometry_error)?);
        }
        let moving = self
            .moving
            .as_mut()
            .expect("moving mixture workspace initialized");
        let original = *walls;
        let before = totals(&self.states, self.mesh.volumes(), species, &mut cancelled)?;
        let mut last = PistonError::Invalid("no admissible two-wall mixture increment");
        let mut first_order_retry = false;
        let mut gas_trials = 0;
        for halvings in 0..=gas_controls.maximum_halvings {
            poll(&mut cancelled)?;
            if !positive(dt)
                || !positive(self.time + dt - self.time)
                || !(self.time + dt).is_finite()
            {
                return Err(PistonError::Invalid(
                    "two-wall mixture step cannot advance the clock",
                ));
            }
            let reconstruction = if first_order_retry {
                Reconstruction::FirstOrder
            } else {
                gas_controls.reconstruction
            };
            let result = solve(
                original.map(|wall| wall.velocity_m_s()),
                controls,
                |speed| {
                    gas_trials += 1;
                    moving
                        .prepare_velocities(&self.mesh, speed, &mut cancelled)
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
                        species,
                        &mut cancelled,
                    )?;
                    let gas = reported_step(before, after, [a, b], dt, area, species)?;
                    let exchange = MechanicalExchange::from(&gas);
                    let displacement = [0, count]
                        .map(|end| moving.mesh.positions()[end] - self.mesh.positions()[end]);
                    mechanical_report(gas, exchange, original, loads, speed, displacement, area)
                },
            );
            match result {
                Ok(mut report) => {
                    let next = [
                        Piston::new(original[0].mass_kg(), report.velocity_after_m_s[0])?,
                        Piston::new(original[1].mass_kg(), report.velocity_after_m_s[1])?,
                    ];
                    report.gas.halvings = halvings;
                    report.gas.first_order_retry =
                        first_order_retry && gas_controls.reconstruction == Reconstruction::Muscl;
                    report.gas_trials = gas_trials;
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.states, &mut self.work.candidate);
                    core::mem::swap(&mut self.mesh, &mut moving.mesh);
                    self.time += dt;
                    *walls = next;
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
