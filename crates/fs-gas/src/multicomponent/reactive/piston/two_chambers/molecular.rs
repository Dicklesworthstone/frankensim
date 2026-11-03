//! Both original molecular/source owners in one shared-piston speed solve.

use super::initial_duration;
use crate::dynamics::piston::{
    MechanicalExchange, Piston, PistonControl, PistonLoad, PistonStep, mechanical_report, solve,
    validate,
};
use crate::multicomponent::dissipation::{DissipationControl, DissipationOperator, TransportLaw};
use crate::multicomponent::reactive::piston::molecular::MovingMolecularTrial;
use crate::multicomponent::reactive::{
    ReactingEuler1d, ReactiveError, ReactiveStepControl, ReactiveStepReport, poll, positive, sum,
};
use fs_qty::Time;

/// Both complete original chamber reports from the same accepted speed trial.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveSharedPistonGasReport {
    /// Left chemistry/transport report, including both molecular half reports.
    pub left: ReactiveStepReport,
    /// Right chemistry/transport report, in its own original source order.
    pub right: ReactiveStepReport,
    /// Physical chemical-energy changes [J], ordered [left, right].
    /// Each chamber uses its own species-reference changes exactly once.
    pub gas_energy_change_j: [f64; 2],
}

impl ReactingEuler1d {
    /// Advance two reacting, diffusing chambers through one finite-mass piston.
    ///
    /// Each candidate midpoint speed runs the original D/C/A/C/D composition
    /// independently in both chambers at one common duration. The left ends
    /// move at [0, speed], the right ends at [speed, 0]. Pressure and longitudinal
    /// molecular stress use those same velocities and their actual integrated
    /// impulses and work in the original scalar mechanical solve. The gas
    /// pressures are both resolved; the supplied force excludes their loads.
    ///
    /// Both original transport laws, species banks, energy references and
    /// chemistry controls remain distinct. D uses the old then moved physical
    /// meshes; C uses the corresponding actual cell volumes. Fourier heat and
    /// Maxwell--Stefan/Soret/Dufour transport retain the supplied original law.
    /// Walls are impermeable and adiabatic. No species, conductive heat or
    /// thermal storage crosses the shared piston; finite piston heat requires
    /// a separate jointly owned thermal inventory.
    ///
    /// Complete endpoint chemical gas energies, piston kinetic storage and
    /// actual applied-force displacement work must balance. Momentum includes
    /// both fixed outer support impulses. No duplicate heat-release or viscous
    /// heating source is added. Symmetric splitting alone proves neither an
    /// accuracy order nor entropy stability.
    ///
    /// A shortened acoustic interval restarts both whole field compositions
    /// and the speed root. Reusable native trial fields, reactors and molecular
    /// workspaces own all proposals. Both fields, clocks, moved meshes, bound
    /// molecular geometries and piston publish together after the final poll.
    /// The two molecular controls are ordered [left, right]; original source
    /// domains and explicit parabolic work budgets remain authoritative.
    ///
    /// # Errors
    /// Invalid source/geometry, clock mismatch, source or transport refusal,
    /// unresolved balances, exhausted retries or cancellation preserves every
    /// accepted physical owner and both molecular geometry bindings.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_piston_and_dissipation<L: TransportLaw, R: TransportLaw>(
        &mut self,
        right: &mut ReactingEuler1d,
        left_operator: &mut DissipationOperator<L>,
        right_operator: &mut DissipationOperator<R>,
        piston: &mut Piston,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: [DissipationControl; 2],
        external_force_n: f64,
        piston_controls: PistonControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep<ReactiveSharedPistonGasReport>, ReactiveError> {
        poll(&mut cancelled)?;
        let load = PistonLoad {
            exterior_pressure_pa: 0.0,
            external_force_n,
        };
        validate(
            *piston,
            maximum_duration,
            controls.transport,
            load,
            piston_controls,
        )?;
        if controls.maximum_coupling_retries > 32 {
            return Err(ReactiveError::Invalid(
                "shared molecular reacting piston retry budget exceeds 32",
            ));
        }
        for molecular_controls in dissipation_controls {
            molecular_controls.validate()?;
        }
        let mut duration = initial_duration(
            &self.accepted,
            &right.accepted,
            maximum_duration,
            controls,
            piston_controls,
            &mut cancelled,
        )?;
        left_operator.admit_motion(&self.accepted, &mut cancelled)?;
        right_operator.admit_motion(&right.accepted, &mut cancelled)?;
        let before = [self.accepted.totals(), right.accepted.totals()];
        let original_piston = *piston;
        let area = self.accepted.mesh.face_areas()[0];
        let interface = self.accepted.cell_count();
        let mut gas_trials = 0;
        let mut last =
            ReactiveError::Invalid("no admissible shared molecular reacting piston increment");
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if !positive(duration)
                || !positive(0.5 * duration)
                || !positive(self.accepted.time + duration - self.accepted.time)
                || !(self.accepted.time + duration).is_finite()
            {
                return Err(ReactiveError::Invalid(
                    "shared molecular reacting piston cannot advance both clocks",
                ));
            }
            let mut shorter_duration = None;
            let attempt = solve(original_piston.velocity_m_s(), piston_controls, |speed| {
                poll(&mut cancelled)?;
                gas_trials += 1;
                let left = match self.moving_molecular_trial(
                    left_operator,
                    None,
                    before[0],
                    duration,
                    controls,
                    dissipation_controls[0],
                    [0.0, speed],
                    retry,
                    &mut cancelled,
                )? {
                    MovingMolecularTrial::Complete(report) => report,
                    MovingMolecularTrial::Shorter(elapsed) => {
                        shorter_duration = Some(elapsed.value());
                        return Err(ReactiveError::Invalid(
                            "left molecular chamber shortened the shared duration",
                        ));
                    }
                };
                let right_report = match right.moving_molecular_trial(
                    right_operator,
                    None,
                    before[1],
                    duration,
                    controls,
                    dissipation_controls[1],
                    [speed, 0.0],
                    retry,
                    &mut cancelled,
                )? {
                    MovingMolecularTrial::Complete(report) => report,
                    MovingMolecularTrial::Shorter(elapsed) => {
                        shorter_duration = Some(elapsed.value());
                        return Err(ReactiveError::Invalid(
                            "right molecular chamber shortened the shared duration",
                        ));
                    }
                };
                if self.trial.mesh.positions()[interface] != right.trial.mesh.positions()[0]
                    || self.trial.time != right.trial.time
                {
                    return Err(ReactiveError::Invalid(
                        "shared molecular chamber endpoints or clocks diverged",
                    ));
                }
                let exchange = MechanicalExchange {
                    duration_s: duration,
                    momentum_change: sum([
                        left.after.momentum - before[0].momentum,
                        right_report.after.momentum - before[1].momentum,
                    ]),
                    energy_change: sum([
                        left.gas_energy_change_j,
                        right_report.gas_energy_change_j,
                    ]),
                    // MechanicalExchange's left slot is the aggregate fixed
                    // support, and its right slot is the one shared piston.
                    left_impulse: sum([
                        left.boundary_impulse_n_s[0],
                        right_report.boundary_impulse_n_s[1],
                    ]),
                    right_impulse: sum([
                        left.boundary_impulse_n_s[1],
                        right_report.boundary_impulse_n_s[0],
                    ]),
                    right_work: sum([left.boundary_work_j[1], right_report.boundary_work_j[0]]),
                };
                let displacement = self.trial.mesh.positions()[interface]
                    - self.accepted.mesh.positions()[interface];
                mechanical_report(
                    ReactiveSharedPistonGasReport {
                        left: left.flow,
                        right: right_report.flow,
                        gas_energy_change_j: [
                            left.gas_energy_change_j,
                            right_report.gas_energy_change_j,
                        ],
                    },
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
                    core::mem::swap(&mut right.accepted, &mut right.trial);
                    left_operator.publish_motion();
                    right_operator.publish_motion();
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
