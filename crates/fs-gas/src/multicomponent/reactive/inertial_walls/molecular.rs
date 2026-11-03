//! Pressure, molecular stress and finite end heat in one two-speed wall solve.

use super::super::piston::molecular::{MolecularPistonReport, MovingMolecularTrial};
use super::super::{
    FlowError, PistonThermalEnds, ReactingEuler1d, ReactiveError, ReactiveStepControl,
    ReactiveStepReport, ThermalReactiveStepReport, poll, positive,
};
use crate::dynamics::inertial_walls::{
    InertialWallLoad, InertialWallsControl, InertialWallsStep, MechanicalExchange,
    mechanical_report, solve_with, validate_controls,
};
use crate::dynamics::piston::{Piston, PistonError};
use crate::multicomponent::dissipation::{
    DissipationControl, DissipationError, DissipationOperator, TransportLaw,
};
use crate::multicomponent::wall_exchange::{WallControl, WallExchangeError};
use crate::multicomponent::{Boundaries, Boundary};
use core::cell::Cell;
use fs_qty::Time;

impl ReactingEuler1d {
    /// Advance both inertial walls with original molecular transport and chemistry.
    ///
    /// Every two-speed proposal reuses the original D/C/A/C/D trial, with old
    /// volumes before A and moved volumes after A. Both candidate wall speeds
    /// supply longitudinal viscous gradients and the same laboratory-frame
    /// stress work. Pressure and viscous impulses jointly drive each finite wall;
    /// gas momentum remains resolved. Chemical conversion and species enthalpy
    /// transport stay with their original source owners.
    ///
    /// The native root admits each wall's impulse/work and combined physical
    /// energy/momentum. Adiabatic impermeable ends carry no conductive boundary
    /// heat. The complete report retains both molecular half steps. Shortened
    /// acoustic durations restart every stage and the whole two-speed root.
    ///
    /// # Errors
    /// Invalid source/geometry, work-budget, nonlinear, balance or cancellation
    /// refusal preserves accepted gas, both walls and the molecular binding.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_inertial_walls_and_dissipation<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        loads: [InertialWallLoad; 2],
        wall_controls: InertialWallsControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<InertialWallsStep<ReactiveStepReport>, ReactiveError> {
        self.step_with_inertial_walls_molecular(
            operator,
            None,
            walls,
            maximum_duration,
            controls,
            dissipation_controls,
            loads,
            wall_controls,
            cancelled,
        )
        .map(|report| report.map_gas(|gas| gas.flow))
    }

    /// Include both original finite end stores or prescribed thermal reservoirs.
    ///
    /// Each proposal runs H/D/C/A/C/D/H from accepted gas and thermal states.
    /// The same complete moving-gas trial supplies the pressure/stress impulses
    /// and work used by both mechanical equations. Conductive heat is separate:
    /// total energy includes actual chemical gas change, original finite-store
    /// energy differences, both wall kinetic changes and exterior/load work,
    /// minus prescribed-reservoir heat. No thermal mass is added to either
    /// Piston and no second chemical or viscous heat source is introduced.
    ///
    /// Gas, mesh, clock, both velocities, molecular geometry and thermal stores
    /// publish after one final cancellation poll. Original caloric conventions,
    /// source domains and explicit acoustic/parabolic budgets still apply.
    /// The held-stage splitting establishes neither a general second-order
    /// method nor an entropy, FSI-stability or combustion-validity theorem.
    ///
    /// # Errors
    /// Source, geometry, work-budget, nonlinear, balance or cancellation failure
    /// leaves all accepted physical owners unchanged and reusable.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_inertial_walls_dissipation_and_heat<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        thermal: &mut PistonThermalEnds<'_>,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        thermal_controls: WallControl,
        loads: [InertialWallLoad; 2],
        wall_controls: InertialWallsControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<InertialWallsStep<ThermalReactiveStepReport>, ReactiveError> {
        self.step_with_inertial_walls_molecular(
            operator,
            Some((thermal, thermal_controls)),
            walls,
            maximum_duration,
            controls,
            dissipation_controls,
            loads,
            wall_controls,
            cancelled,
        )
        .map(|report| {
            report.map_gas(|gas| ThermalReactiveStepReport {
                reactive: gas.flow,
                before: gas.before,
                after: gas.after,
                heat: gas
                    .heat
                    .expect("thermal two-wall entry supplies both heat halves"),
                gas_energy_change_j: gas.gas_energy_change_j,
            })
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with_inertial_walls_molecular<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        mut thermal: Option<(&mut PistonThermalEnds<'_>, WallControl)>,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        loads: [InertialWallLoad; 2],
        wall_controls: InertialWallsControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<InertialWallsStep<MolecularPistonReport>, ReactiveError> {
        poll(&mut cancelled)?;
        validate_controls(
            *walls,
            maximum_duration,
            controls.transport,
            loads,
            wall_controls,
        )?;
        dissipation_controls.validate()?;
        if controls.maximum_coupling_retries > 32
            || self.accepted.boundaries
                != (Boundaries::Ends {
                    left: Boundary::Reflecting,
                    right: Boundary::Reflecting,
                })
            || self.accepted.mesh.uniform_cross_section().is_none()
        {
            return Err(ReactiveError::Invalid(
                "molecular inertial walls require bounded retries and a constant-area closed interval",
            ));
        }
        if let Some((ends, heat_controls)) = thermal.as_ref() {
            heat_controls.validate()?;
            ends.admit(&self.accepted)?;
        }
        operator.admit_motion(&self.accepted, &mut cancelled)?;
        let original = *walls;
        let before = self.accepted.totals();
        let area = self.accepted.mesh.face_areas()[0];
        let bounds = wall_controls.midpoint_velocity_bracket_m_s;
        let low = bounds[0][0].min(bounds[1][0]);
        let high = bounds[0][1].max(bounds[1][1]);
        let mut duration = maximum_duration.value();
        for i in 0..self.accepted.cell_count() {
            let p = self.accepted.gas().primitive(
                self.accepted.cell(i).expect("admitted cell"),
                &mut cancelled,
            )?;
            let relative = (p.velocity - low).abs().max((p.velocity - high).abs());
            let rate = (relative + p.sound_speed) * area / self.accepted.mesh.volumes()[i];
            if !positive(rate) {
                return Err(ReactiveError::Invalid(
                    "molecular two-wall acoustic bound overflow",
                ));
            }
            duration = duration.min(0.8 * controls.transport.cfl / rate);
        }
        let right = self.accepted.cell_count();
        let contraction = bounds[0][1] - bounds[1][0];
        if !contraction.is_finite() {
            return Err(ReactiveError::Invalid(
                "molecular two-wall contraction bound overflow",
            ));
        }
        if contraction > 0.0 {
            let length = self.accepted.mesh.positions()[right] - self.accepted.mesh.positions()[0];
            duration = duration.min(0.8 * controls.transport.cfl * length / contraction);
        }
        let mut last = ReactiveError::Invalid("no admissible molecular two-wall increment");
        let mut gas_trials = 0;
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if !positive(duration)
                || !positive(0.5 * duration)
                || !positive(self.accepted.time + duration - self.accepted.time)
                || !(self.accepted.time + duration).is_finite()
            {
                return Err(ReactiveError::Invalid(
                    "molecular two-wall step cannot advance the clock",
                ));
            }
            // Both closures borrow only the Cell; no mutable borrow escapes a
            // trial. A shortened field cannot enter line-search acceptance.
            let shorter = Cell::new(None);
            let attempt = solve_with(
                original.map(|wall| wall.velocity_m_s()),
                wall_controls,
                |speed| {
                    poll(&mut cancelled)?;
                    gas_trials += 1;
                    let trial = self.moving_molecular_trial(
                        operator,
                        thermal
                            .as_mut()
                            .map(|(ends, control)| (&mut **ends, *control)),
                        before,
                        duration,
                        controls,
                        dissipation_controls,
                        speed,
                        retry,
                        &mut cancelled,
                    )?;
                    let gas = match trial {
                        MovingMolecularTrial::Shorter(elapsed) => {
                            shorter.set(Some(elapsed.value()));
                            return Err(ReactiveError::Invalid(
                                "moving transport shortened the complete molecular two-wall duration",
                            ));
                        }
                        MovingMolecularTrial::Complete(gas) => gas,
                    };
                    let exchange = MechanicalExchange {
                        duration_s: duration,
                        momentum_change: gas.after.momentum - gas.before.momentum,
                        energy_change: gas.accounted_energy_change_j,
                        boundary_impulse: gas.boundary_impulse_n_s,
                        boundary_work: gas.boundary_work_j,
                    };
                    let displacement = [0, right].map(|end| {
                        self.trial.mesh.positions()[end] - self.accepted.mesh.positions()[end]
                    });
                    mechanical_report(gas, exchange, original, loads, speed, displacement, area)
                        .map_err(ReactiveError::from)
                },
                |error| {
                    shorter.get().is_some()
                        || matches!(error, ReactiveError::Cancelled)
                        || allocation_failure(error)
                },
            );
            match attempt {
                Ok(mut report) => {
                    let next = [
                        Piston::new(original[0].mass_kg(), report.velocity_after_m_s[0])?,
                        Piston::new(original[1].mass_kg(), report.velocity_after_m_s[1])?,
                    ];
                    report.gas_trials = gas_trials;
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.accepted, &mut self.trial);
                    operator.publish_motion();
                    if let Some((ends, _)) = thermal.as_mut() {
                        ends.publish(report.gas.heat.expect("accepted thermal two-wall trial"));
                    }
                    *walls = next;
                    return Ok(report);
                }
                Err(ReactiveError::Cancelled) => return Err(ReactiveError::Cancelled),
                Err(error) if allocation_failure(&error) => return Err(error),
                Err(error) => {
                    last = error;
                    duration = shorter.get().unwrap_or(0.5 * duration);
                }
            }
        }
        Err(ReactiveError::RetryLimit {
            attempts: controls.maximum_coupling_retries + 1,
            last: Box::new(last),
        })
    }
}

fn allocation_failure(error: &ReactiveError) -> bool {
    matches!(
        error,
        ReactiveError::Flow(FlowError::Allocation)
            | ReactiveError::Piston(PistonError::MixtureGas(FlowError::Allocation))
            | ReactiveError::Dissipation(DissipationError::Flow(FlowError::Allocation))
            | ReactiveError::Wall(WallExchangeError::Flow(FlowError::Allocation))
    )
}
