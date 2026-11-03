//! One common-speed root for pressure, longitudinal stress and chemical energy.

use super::super::{
    ChemistryReport, ReactingEuler1d, ReactiveError, ReactiveStepControl, ReactiveStepReport,
    SplitDissipationReport, coupled_report, poll, positive, react_cells, sum,
};
use crate::dynamics::piston::{
    MechanicalExchange, Piston, PistonControl, PistonLoad, PistonStep, mechanical_report, solve,
    validate,
};
use crate::multicomponent::dissipation::{DissipationControl, DissipationOperator, TransportLaw};
use crate::multicomponent::wall_exchange::WallControl;
use crate::multicomponent::{Boundaries, Boundary, Totals};
use fs_qty::Time;

use super::thermal::{PistonThermalEnds, SplitPistonThermalReport};

/// Private common-root result; public entries expose their original report.
pub(in crate::multicomponent::reactive) struct MolecularPistonReport {
    pub flow: ReactiveStepReport,
    pub before: Totals,
    pub after: Totals,
    pub heat: Option<SplitPistonThermalReport>,
    pub gas_energy_change_j: f64,
    /// Signed pressure plus viscous impulses entering the gas, [left, right].
    pub boundary_impulse_n_s: [f64; 2],
    /// Mechanical pressure/stress work entering the gas, excluding end heat.
    pub boundary_work_j: [f64; 2],
    /// Gas and finite-wall storage changes minus prescribed reservoir heat.
    pub accounted_energy_change_j: f64,
}

/// Inline reports reuse native workspaces without allocating per root proposal.
#[allow(clippy::large_enum_variant)]
pub(in crate::multicomponent::reactive) enum MovingMolecularTrial {
    Shorter(Time),
    Complete(MolecularPistonReport),
}

impl ReactingEuler1d {
    /// Reuse the original molecular/source/heat stages for one proposed pair
    /// of end speeds. Every proposal starts from accepted physical endpoints.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::multicomponent::reactive) fn moving_molecular_trial<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        mut thermal: Option<(&mut PistonThermalEnds<'_>, WallControl)>,
        before: Totals,
        duration: f64,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        wall_velocity: [f64; 2],
        retry: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingMolecularTrial, ReactiveError> {
        self.trial
            .copy_physical_with_motion_from(&self.accepted, cancelled)?;
        let first_heat = if let Some((walls, wall_controls)) = thermal.as_mut() {
            walls.begin_trial();
            Some(walls.advance_trial(
                &mut self.trial,
                Time::new(0.5 * duration),
                *wall_controls,
                cancelled,
            )?)
        } else {
            None
        };
        let middle_before = if first_heat.is_some() {
            self.trial.totals()
        } else {
            before
        };
        let first = operator.advance_moving_trial(
            &mut self.trial,
            Time::new(0.5 * duration),
            dissipation_controls,
            wall_velocity,
            cancelled,
        )?;
        let mut chemistry = ChemistryReport::default();
        react_cells(
            &mut self.trial,
            &mut self.reactor,
            &mut self.amounts,
            0.5 * duration,
            &mut chemistry,
            cancelled,
        )?;
        let transport = self.trial.step_with_moving_walls(
            Time::new(duration),
            controls.transport,
            wall_velocity,
            &mut *cancelled,
        )?;
        if transport.elapsed.value() != duration {
            return Ok(MovingMolecularTrial::Shorter(transport.elapsed));
        }
        react_cells(
            &mut self.trial,
            &mut self.reactor,
            &mut self.amounts,
            0.5 * duration,
            &mut chemistry,
            cancelled,
        )?;
        let second = operator.advance_moving_trial(
            &mut self.trial,
            Time::new(0.5 * duration),
            dissipation_controls,
            wall_velocity,
            cancelled,
        )?;
        let middle_after = self.trial.totals();
        let heat = if let Some((walls, wall_controls)) = thermal.as_mut() {
            let second_heat = walls.advance_trial(
                &mut self.trial,
                Time::new(0.5 * duration),
                *wall_controls,
                cancelled,
            )?;
            Some(walls.complete(
                first_heat.expect("both heat halves share one thermal owner"),
                second_heat,
            )?)
        } else {
            None
        };
        let interval = self.reactor.controls();
        for i in 0..self.trial.cell_count() {
            let primitive = self.trial.gas().primitive(
                self.trial.cell(i).expect("admitted cell index"),
                &mut *cancelled,
            )?;
            if !(interval.minimum_temperature_k..=interval.maximum_temperature_k)
                .contains(&primitive.temperature)
            {
                return Err(ReactiveError::Invalid(
                    "final molecular piston state is outside chemistry interval",
                ));
            }
        }
        let after = if heat.is_some() {
            self.trial.totals()
        } else {
            middle_after
        };
        let flow = coupled_report(
            &self.trial.gas,
            middle_before,
            middle_after,
            transport,
            chemistry,
            Some(SplitDissipationReport { first, second }),
            retry,
        )?;
        // Closed ends carry no species. Restore each source reference once,
        // using actual accepted-to-candidate species and energy differences.
        let gas_energy_change_j = sum(core::iter::once(after.energy - before.energy).chain(
            (0..self.trial.gas.species_count()).map(|k| {
                (after.species_masses[k] - before.species_masses[k])
                    * self.trial.gas.energy_references()[k]
            }),
        ));
        let boundary_impulse_n_s = [
            transport.left_boundary_change.momentum
                + first.left_boundary_change.momentum
                + second.left_boundary_change.momentum,
            transport.right_boundary_change.momentum
                + first.right_boundary_change.momentum
                + second.right_boundary_change.momentum,
        ];
        let boundary_work_j = [
            transport.left_boundary_change.energy
                + first.left_boundary_change.energy
                + second.left_boundary_change.energy,
            transport.right_boundary_change.energy
                + first.right_boundary_change.energy
                + second.right_boundary_change.energy,
        ];
        // Finite heat storage is internal. Prescribed reservoir heat is
        // external and never part of the mechanical wall-work arrays.
        let accounted_energy_change_j = match heat {
            Some(heat) => sum([
                gas_energy_change_j,
                heat.wall_energy_change_j,
                -heat.external_heat_to_gas_j,
            ]),
            None => gas_energy_change_j,
        };
        Ok(MovingMolecularTrial::Complete(MolecularPistonReport {
            flow,
            before,
            after,
            heat,
            gas_energy_change_j,
            boundary_impulse_n_s,
            boundary_work_j,
            accounted_energy_change_j,
        }))
    }

    /// Advance chemistry, molecular transport and pressure-driven piston motion.
    ///
    /// Each midpoint-speed proposal starts from the accepted state and runs
    /// D(dt/2), C(dt/2), the original moving-volume Euler A(dt), C(dt/2), D(dt/2).
    /// The same candidate wall speed supplies the normal viscous gradient,
    /// viscous wall work, moving-face pressure work and mechanical impulse root.
    /// Both pressure and longitudinal viscous impulses drive the piston. The
    /// molecular stages freeze the old/new physical mesh respectively; only A
    /// changes volumes. Diffusion carries original species enthalpies, while C
    /// alone owns chemical conversion. No independent heat source is added.
    ///
    /// Requires two impermeable adiabatic ends, a fixed left end and constant
    /// area. Thermal conduction and Maxwell--Stefan/Soret/Dufour transport act
    /// inside the gas according to the supplied original law; there is no wall
    /// reservoir heat or finite wall storage in this entry. Thus both boundary
    /// energy transfers are mechanical work, not conductive reservoir heat.
    /// Symmetric splitting is not by itself a second-order or entropy proof.
    ///
    /// Shorter acoustic durations restart every stage and the entire velocity
    /// root. Existing field/source/molecular workspaces and scalar solve are
    /// reused. The operator lazily allocates one candidate mesh; rejected root
    /// proposals never rebind its accepted mesh. Gas, geometry, clock, piston
    /// velocity and molecular geometry publish together after the final poll.
    ///
    /// The returned native reactive report retains both molecular half reports,
    /// including their separate signed left/right stress impulse and work.
    /// Piston energy acceptance uses actual endpoint chemical gas energy plus
    /// piston kinetic energy, external-force work and exterior-pressure work.
    /// Original source domains, explicit parabolic budgets and unequal-caloric
    /// material-contact limitations remain in force.
    ///
    /// # Errors
    /// Identity/geometry, source, acoustic/parabolic, nonlinear, balance or
    /// cancellation refusal preserves every accepted physical owner.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_piston_and_dissipation<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        piston: &mut Piston,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        load: PistonLoad,
        piston_controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep<ReactiveStepReport>, ReactiveError> {
        self.step_with_piston_molecular(
            operator,
            None,
            piston,
            maximum_duration,
            controls,
            dissipation_controls,
            load,
            piston_controls,
            cancelled,
        )
        .map(|report| report.map_gas(|gas| gas.flow))
    }

    /// One original pressure/stress root, with optional original heat halves.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn step_with_piston_molecular<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        mut thermal: Option<(&mut PistonThermalEnds<'_>, WallControl)>,
        piston: &mut Piston,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        load: PistonLoad,
        piston_controls: PistonControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep<MolecularPistonReport>, ReactiveError> {
        poll(&mut cancelled)?;
        validate(
            *piston,
            maximum_duration,
            controls.transport,
            load,
            piston_controls,
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
                "molecular reacting piston requires bounded retries and a constant-area closed interval",
            ));
        }
        if let Some((walls, wall_controls)) = thermal.as_ref() {
            wall_controls.validate()?;
            walls.admit(&self.accepted)?;
        }
        operator.admit_motion(&self.accepted, &mut cancelled)?;
        let original_piston = *piston;
        let before = self.accepted.totals();
        let area = self.accepted.mesh.face_areas()[0];
        let [lo, hi] = piston_controls.midpoint_velocity_bracket_m_s;
        let mut duration = maximum_duration.value();
        for i in 0..self.accepted.cell_count() {
            let primitive = self.accepted.gas().primitive(
                self.accepted.cell(i).expect("admitted cell index"),
                &mut cancelled,
            )?;
            let relative = primitive
                .velocity
                .abs()
                .max((primitive.velocity - lo).abs())
                .max((primitive.velocity - hi).abs());
            let rate =
                (relative + primitive.sound_speed) * (area / self.accepted.mesh.volumes()[i]);
            if !positive(rate) {
                return Err(ReactiveError::Invalid(
                    "molecular piston relative-wave CFL overflow",
                ));
            }
            duration = duration.min(0.8 * controls.transport.cfl / rate);
        }
        let right = self.accepted.cell_count();
        if lo < 0.0 {
            let length = self.accepted.mesh.positions()[right] - self.accepted.mesh.positions()[0];
            duration = duration.min(0.8 * controls.transport.cfl * (length / -lo));
        }
        let mut last = ReactiveError::Invalid("no admissible molecular reacting piston increment");
        let mut gas_trials = 0;
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if !positive(duration)
                || !positive(self.accepted.time + duration - self.accepted.time)
                || !(self.accepted.time + duration).is_finite()
            {
                return Err(ReactiveError::Invalid(
                    "molecular piston cannot advance the clock",
                ));
            }
            let mut shorter_duration = None;
            let attempt = solve(original_piston.velocity_m_s(), piston_controls, |speed| {
                poll(&mut cancelled)?;
                gas_trials += 1;
                let gas = match self.moving_molecular_trial(
                    operator,
                    thermal
                        .as_mut()
                        .map(|(walls, controls)| (&mut **walls, *controls)),
                    before,
                    duration,
                    controls,
                    dissipation_controls,
                    [0.0, speed],
                    retry,
                    &mut cancelled,
                )? {
                    MovingMolecularTrial::Shorter(elapsed) => {
                        shorter_duration = Some(elapsed.value());
                        return Err(ReactiveError::Invalid(
                            "moving transport shortened the molecular piston duration",
                        ));
                    }
                    MovingMolecularTrial::Complete(report) => report,
                };
                let exchange = MechanicalExchange {
                    duration_s: duration,
                    momentum_change: gas.after.momentum - gas.before.momentum,
                    energy_change: gas.accounted_energy_change_j,
                    left_impulse: gas.boundary_impulse_n_s[0],
                    right_impulse: gas.boundary_impulse_n_s[1],
                    right_work: gas.boundary_work_j[1],
                };
                let displacement =
                    self.trial.mesh.positions()[right] - self.accepted.mesh.positions()[right];
                mechanical_report(
                    gas,
                    exchange,
                    original_piston,
                    load,
                    speed,
                    displacement,
                    area,
                )
                .map_err(ReactiveError::from)
            });
            match attempt {
                Ok(mut report) => {
                    report.gas_trials = gas_trials;
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.accepted, &mut self.trial);
                    operator.publish_motion();
                    if let Some((walls, _)) = thermal.as_mut() {
                        walls.publish(report.gas.heat.expect("accepted thermal piston heat"));
                    }
                    piston.accept(&report);
                    return Ok(report);
                }
                Err(ReactiveError::Cancelled) => return Err(ReactiveError::Cancelled),
                Err(error) => {
                    last = error;
                    duration = shorter_duration.unwrap_or(0.5 * duration);
                }
            }
        }
        Err(ReactiveError::RetryLimit {
            attempts: controls.maximum_coupling_retries + 1,
            last: Box::new(last),
        })
    }
}
