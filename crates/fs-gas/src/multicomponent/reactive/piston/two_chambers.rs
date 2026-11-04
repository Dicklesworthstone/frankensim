//! Original chemistry on both sides of one resolved, finite-mass piston.

use super::super::{
    ChemistryReport, ReactingEuler1d, ReactiveError, ReactiveStepControl, ReactiveStepReport,
    coupled_report, poll, positive, react_cells, sum,
};
use crate::dynamics::piston::{
    Piston, PistonControl, PistonError, PistonLoad, PistonStep, PistonSuspension,
    SuspendedPistonStep, validate,
};
use crate::multicomponent::{Boundaries, Boundary, MixtureEuler1d, SharedPistonGasReport};
use fs_qty::Time;

mod molecular;
pub use molecular::ReactiveSharedPistonGasReport;

/// One accepted pair of reacting chambers and their common mechanical response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveSharedPistonStep<P = PistonStep<SharedPistonGasReport>> {
    /// Complete left-chamber chemistry and transport, in its original species order.
    pub left: ReactiveStepReport,
    /// Complete right-chamber chemistry and transport, in its original species order.
    pub right: ReactiveStepReport,
    /// Original simultaneous two-field pressure/piston solve between chemistry halves.
    pub piston: P,
    /// Actual physical chemical-energy changes [J], ordered [left, right].
    /// These include each source bank's species-reference changes.
    pub gas_energy_change_j: [f64; 2],
    /// Whole-step gas/piston momentum change plus outer support impulses,
    /// minus the independently applied piston and suspension impulses [N s].
    pub combined_momentum_residual_n_s: f64,
    /// Both physical gas-energy changes plus piston kinetic change, spring
    /// storage and damper loss when present, minus authored drive work [J].
    pub total_energy_residual_j: f64,
}

impl ReactingEuler1d {
    /// Advance two reacting gas chambers through their shared finite piston.
    ///
    /// Each chamber retains its own original mechanism, species bank, caloric
    /// references and chemistry controls. Both first chemistry half steps use
    /// accepted physical cell volumes. The original two-field pressure solve
    /// advances the shared piston and both meshes at one common duration;
    /// both second chemistry halves use those moved volumes. Only chemistry
    /// converts formation energy. No separate heat-release term is introduced.
    ///
    /// A shorter acoustic or mechanical interval restarts both chambers'
    /// complete source/pressure/source composition from their accepted states.
    /// Reusable native trial fields, reactors and amount buffers own the work.
    /// Neither chamber, clock nor piston publishes before whole-system physical
    /// energy and momentum admission and the final cancellation checkpoint.
    ///
    /// The left chamber is `self`; the right chamber is `right`. Closed ends,
    /// matching constant area, coincident interface coordinates and identical
    /// accepted clocks are required. Positive applied force points right and
    /// excludes both resolved gas pressure loads. The two outer fixed supports
    /// retain their explicit momentum transfer and do no work.
    ///
    /// This composes adiabatic chemistry and resolved Euler transport. Molecular
    /// transport, finite-wall heat, leakage and gravity retain their separate
    /// owners. Symmetric splitting alone does not prove second-order accuracy,
    /// entropy stability or validity of a combustion mechanism.
    ///
    /// # Errors
    /// Invalid geometry/controls, original source refusals, exhausted retries,
    /// unresolved balances or cancellation leave all accepted states unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_piston(
        &mut self,
        right: &mut ReactingEuler1d,
        piston: &mut Piston,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        external_force_n: f64,
        piston_controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveSharedPistonStep, ReactiveError> {
        self.step_with_shared_piston_suspension(
            right, piston, maximum_duration, controls, external_force_n,
            None, piston_controls, cancelled,
        ).map(|step| ReactiveSharedPistonStep {
            left: step.left,
            right: step.right,
            piston: step.piston.piston,
            gas_energy_change_j: step.gas_energy_change_j,
            combined_momentum_residual_n_s: step.combined_momentum_residual_n_s,
            total_energy_residual_j: step.total_energy_residual_j,
        })
    }

    /// Advance both reacting chambers with a shared spring-and-damper piston.
    ///
    /// The original C/P/C composition retains each mechanism and chemical
    /// energy datum. Its pressure stage solves both resolved fields together
    /// with the stationary [`PistonSuspension`]. Whole-step admission includes
    /// both physical chemical energies, kinetic and spring storage, damper loss,
    /// authored drive work, and both fixed-support and suspension impulses.
    /// The nested suspended report covers the pressure stage; the outer balances
    /// additionally include both chemistry halves at their actual cell volumes.
    ///
    /// The suspension is held during the complete step. Changing its anchor or
    /// stiffness between calls requires caller accounting of that work. Damper
    /// loss leaves the mechanical system and is not added as gas heat.
    ///
    /// # Errors
    /// The same source, geometry, budget and cancellation refusals as
    /// [`Self::step_with_shared_piston`] preserve both fields and the piston.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_suspended_piston(
        &mut self,
        right: &mut ReactingEuler1d,
        piston: &mut Piston,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        external_force_n: f64,
        suspension: PistonSuspension,
        piston_controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveSharedPistonStep<SuspendedPistonStep<SharedPistonGasReport>>, ReactiveError> {
        self.step_with_shared_piston_suspension(
            right, piston, maximum_duration, controls, external_force_n,
            Some(suspension), piston_controls, cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with_shared_piston_suspension(
        &mut self,
        right: &mut ReactingEuler1d,
        piston: &mut Piston,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        external_force_n: f64,
        suspension: Option<PistonSuspension>,
        piston_controls: PistonControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveSharedPistonStep<SuspendedPistonStep<SharedPistonGasReport>>, ReactiveError> {
        poll(&mut cancelled)?;
        validate(
            *piston,
            maximum_duration,
            controls.transport,
            PistonLoad {
                exterior_pressure_pa: 0.0,
                external_force_n,
            },
            piston_controls,
        )?;
        if controls.maximum_coupling_retries > 32 {
            return Err(ReactiveError::Invalid(
                "shared reacting piston retry budget exceeds 32",
            ));
        }
        let mut duration = initial_duration(
            &self.accepted,
            &right.accepted,
            maximum_duration,
            controls,
            piston_controls,
            &mut cancelled,
        )?;
        let before = [self.accepted.totals(), right.accepted.totals()];
        let original_piston = *piston;
        let mut last = ReactiveError::Invalid("no admissible shared reacting piston increment");
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if !positive(duration)
                || !positive(0.5 * duration)
                || !positive(self.accepted.time + duration - self.accepted.time)
                || !(self.accepted.time + duration).is_finite()
            {
                return Err(ReactiveError::Invalid(
                    "shared reacting piston cannot advance both clocks",
                ));
            }
            self.trial
                .copy_physical_with_motion_from(&self.accepted, &mut cancelled)?;
            right
                .trial
                .copy_physical_with_motion_from(&right.accepted, &mut cancelled)?;
            let mut trial_piston = original_piston;
            let mut chemistry = [ChemistryReport::default(); 2];
            let attempt = (|| {
                react_cells(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut chemistry[0],
                    &mut cancelled,
                )?;
                react_cells(
                    &mut right.trial,
                    &mut right.reactor,
                    &mut right.amounts,
                    0.5 * duration,
                    &mut chemistry[1],
                    &mut cancelled,
                )?;
                let suspended = self.trial.step_with_shared_piston_suspension(
                    &mut right.trial,
                    &mut trial_piston,
                    Time::new(duration),
                    controls.transport,
                    external_force_n,
                    suspension,
                    piston_controls,
                    &mut cancelled,
                )?;
                let mechanical = &suspended.piston;
                let elapsed = mechanical.gas.left.elapsed;
                if elapsed != mechanical.gas.right.elapsed {
                    return Err(ReactiveError::Invalid(
                        "shared pressure solve returned unequal durations",
                    ));
                }
                if elapsed.value() != duration {
                    return Ok((elapsed.value(), None));
                }
                react_cells(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut chemistry[0],
                    &mut cancelled,
                )?;
                react_cells(
                    &mut right.trial,
                    &mut right.reactor,
                    &mut right.amounts,
                    0.5 * duration,
                    &mut chemistry[1],
                    &mut cancelled,
                )?;
                let left = coupled_report(
                    &self.trial.gas,
                    before[0],
                    self.trial.totals(),
                    mechanical.gas.left,
                    chemistry[0],
                    None,
                    retry,
                )?;
                let right_report = coupled_report(
                    &right.trial.gas,
                    before[1],
                    right.trial.totals(),
                    mechanical.gas.right,
                    chemistry[1],
                    None,
                    retry,
                )?;
                // Each original source report reconstructs physical chemical
                // energy using its OWN species references. Never combine gauges
                // or infer chemical heat from a second constitutive formula.
                let gas_energy_change_j = [&left, &right_report].map(|flow| {
                    sum([
                        flow.absolute_energy_residual_j,
                        flow.transport.boundary_change.energy,
                    ])
                });
                let mut total_energy_residual_j = sum([
                    gas_energy_change_j[0],
                    gas_energy_change_j[1],
                    mechanical.kinetic_energy_change_j,
                    -mechanical.external_force_work_j,
                ]);
                let mut combined_momentum_residual_n_s = sum([
                    left.after.momentum - left.before.momentum,
                    right_report.after.momentum - right_report.before.momentum,
                    original_piston.mass_kg()
                        * (trial_piston.velocity_m_s() - original_piston.velocity_m_s()),
                    mechanical.fixed_wall_impulse_n_s,
                    -duration * external_force_n,
                ]);
                if suspension.is_some() {
                    total_energy_residual_j = sum([
                        gas_energy_change_j[0],
                        gas_energy_change_j[1],
                        mechanical.kinetic_energy_change_j,
                        suspended.spring_energy_change_j,
                        suspended.damping_energy_j,
                        -suspended.applied_force_work_j,
                    ]);
                    combined_momentum_residual_n_s -= suspended.suspension_impulse_n_s;
                }
                if !gas_energy_change_j
                    .into_iter()
                    .chain([total_energy_residual_j, combined_momentum_residual_n_s])
                    .all(f64::is_finite)
                {
                    return Err(ReactiveError::Invalid(
                        "non-finite shared chemical/piston balance",
                    ));
                }
                if total_energy_residual_j.abs() > piston_controls.energy_tolerance_j
                    || combined_momentum_residual_n_s.abs() > piston_controls.momentum_tolerance_n_s
                {
                    return Err(ReactiveError::Piston(PistonError::Balance {
                        momentum_residual_n_s: combined_momentum_residual_n_s.abs(),
                        energy_residual_j: total_energy_residual_j.abs(),
                    }));
                }
                Ok::<_, ReactiveError>((
                    duration,
                    Some(ReactiveSharedPistonStep {
                        left,
                        right: right_report,
                        piston: suspended,
                        gas_energy_change_j,
                        combined_momentum_residual_n_s,
                        total_energy_residual_j,
                    }),
                ))
            })();
            match attempt {
                Ok((shorter, None)) => {
                    duration = shorter;
                    last = ReactiveError::Invalid(
                        "pressure shortened both reacting chamber intervals",
                    );
                }
                Ok((_, Some(report))) => {
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.accepted, &mut self.trial);
                    core::mem::swap(&mut right.accepted, &mut right.trial);
                    *piston = trial_piston;
                    return Ok(report);
                }
                Err(ReactiveError::Cancelled) => return Err(ReactiveError::Cancelled),
                Err(error) => {
                    last = error;
                    duration *= 0.5;
                }
            }
        }
        Err(ReactiveError::RetryLimit {
            attempts: controls.maximum_coupling_retries + 1,
            last: Box::new(last),
        })
    }
}

fn initial_duration(
    left: &MixtureEuler1d,
    right: &MixtureEuler1d,
    maximum_duration: Time,
    controls: ReactiveStepControl,
    piston: PistonControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ReactiveError> {
    let closed = Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    };
    let area = left
        .mesh
        .uniform_cross_section()
        .ok_or(ReactiveError::Invalid(
            "shared reacting piston requires constant-area chambers",
        ))?;
    let interface = left.mesh.positions()[left.cell_count()];
    if left.boundaries != closed
        || right.boundaries != closed
        || right.mesh.uniform_cross_section() != Some(area)
        || right.mesh.positions()[0] != interface
        || left.time != right.time
    {
        return Err(ReactiveError::Invalid(
            "shared reacting piston requires closed chambers with equal areas, interface coordinates and clocks",
        ));
    }
    let [lo, hi] = piston.midpoint_velocity_bracket_m_s;
    let mut duration = maximum_duration.value();
    for flow in [left, right] {
        for i in 0..flow.cell_count() {
            let p = flow.gas().primitive(
                flow.cell(i).expect("admitted chamber cell"),
                &mut *cancelled,
            )?;
            let relative = p
                .velocity
                .abs()
                .max((p.velocity - lo).abs())
                .max((p.velocity - hi).abs());
            let rate = (relative + p.sound_speed) * area / flow.mesh.volumes()[i];
            if !positive(rate) {
                return Err(ReactiveError::Invalid(
                    "shared reacting piston acoustic bound overflow",
                ));
            }
            duration = duration.min(0.8 * controls.transport.cfl / rate);
        }
    }
    if lo < 0.0 {
        duration = duration
            .min(0.8 * controls.transport.cfl * ((interface - left.mesh.positions()[0]) / -lo));
    }
    if hi > 0.0 {
        duration = duration.min(
            0.8 * controls.transport.cfl
                * ((right.mesh.positions()[right.cell_count()] - interface) / hi),
        );
    }
    Ok(duration)
}
