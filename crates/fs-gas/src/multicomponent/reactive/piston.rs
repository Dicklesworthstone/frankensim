//! Finite-rate source chemistry driving resolved pressure and inertial wall work.

use super::{
    ChemistryReport, ReactingEuler1d, ReactiveError, ReactiveStepControl, ReactiveStepReport,
    TransportStepReport, coupled_report, poll, positive, react_cells,
};
use crate::dynamics::piston::{Piston, PistonControl, PistonError, PistonLoad, PistonStep};
use crate::multicomponent::{Boundaries, Boundary};
use fs_qty::Time;

/// One atomic chemistry/pressure-wave/piston step on moving physical volumes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactivePistonStep {
    /// Full source/boundary balances, including both chemistry half steps and
    /// the actual endpoint geometry. Its transport lies between chemistry.
    pub flow: ReactiveStepReport,
    /// Original common-duration pressure/piston solve between chemistry steps.
    /// Its energy residual covers that middle stage, not chemical conversion.
    pub piston: PistonStep<TransportStepReport>,
    /// Actual gas physical chemical-energy change [J], including each changed
    /// species reference. It is not the stored sensible/gauged-energy change.
    pub gas_energy_change_j: f64,
    /// Whole-step gas plus piston momentum residual [N s], including the fixed
    /// left-wall reaction and exterior-pressure/additional-force impulse.
    pub combined_momentum_residual_n_s: f64,
    /// Whole-step physical chemical gas energy plus piston kinetic change and
    /// exterior-pressure work, minus additional-force work [J].
    pub total_energy_residual_j: f64,
}

impl ReactingEuler1d {
    /// Advance finite-rate chemistry, resolved pressure and one inertial piston.
    ///
    /// The symmetric composition is C(dt/2), the original spatial mixture
    /// piston step P(dt), then C(dt/2). The first source solve uses old physical
    /// cell volumes; the second uses accepted candidate volumes after motion.
    /// Formation-energy conversion belongs only to the original reactor and
    /// species-reference correction. No additional heat-release term is added.
    ///
    /// A shorter pressure/CFL/root duration restarts BOTH chemistry halves and
    /// piston motion from the original accepted field and velocity. A source or
    /// final balance refusal halves the complete attempt. The original source
    /// controls and finite work budget remain authoritative for every cell.
    /// Gas, moved mesh, clock and piston velocity publish together after actual
    /// chemical-energy/kinetic/external-work and momentum acceptance. The
    /// existing two field workspaces, reactor and amount buffers are reused.
    ///
    /// `piston_controls` admits whole-step energy and momentum defects as well
    /// as the middle pressure solve. Its energy tolerance must account for the
    /// caller's original per-cell chemistry accuracy. The left wall is fixed;
    /// area is constant and both ends are impermeable. The separately bound
    /// fixed-geometry molecular and finite-wall operators are not advanced.
    /// Symmetric composition alone does not establish second-order accuracy,
    /// pressure-equilibrium preservation, entropy stability or combustion validity.
    ///
    /// # Errors
    /// Source, geometry, mechanical, work-budget, balance or cancellation failure
    /// leaves every accepted cell, mesh coordinate, clock and piston unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_piston(
        &mut self,
        piston: &mut Piston,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        load: PistonLoad,
        piston_controls: PistonControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactivePistonStep, ReactiveError> {
        poll(&mut cancelled)?;
        crate::dynamics::piston::validate(
            *piston,
            maximum_duration,
            controls.transport,
            load,
            piston_controls,
        )?;
        if controls.maximum_coupling_retries > 32
            || self.accepted.boundaries
                != (Boundaries::Ends {
                    left: Boundary::Reflecting,
                    right: Boundary::Reflecting,
                })
            || self.accepted.mesh.uniform_cross_section().is_none()
        {
            return Err(ReactiveError::Invalid(
                "reacting piston requires bounded retries and a constant-area closed interval",
            ));
        }
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
                    "reacting piston relative-wave CFL overflow",
                ));
            }
            duration = duration.min(0.8 * controls.transport.cfl / rate);
        }
        if lo < 0.0 {
            let length = self.accepted.mesh.positions()[self.accepted.cell_count()]
                - self.accepted.mesh.positions()[0];
            duration = duration.min(0.8 * controls.transport.cfl * (length / -lo));
        }
        let original_piston = *piston;
        let mut last = ReactiveError::Invalid("no admissible reacting piston increment");
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            if !positive(duration)
                || !positive(self.accepted.time + duration - self.accepted.time)
                || !(self.accepted.time + duration).is_finite()
            {
                return Err(ReactiveError::Invalid(
                    "reacting piston cannot advance the clock",
                ));
            }
            self.trial
                .copy_physical_with_motion_from(&self.accepted, &mut cancelled)?;
            let mut trial_piston = original_piston;
            let mut chemistry = ChemistryReport::default();
            let attempt = (|| {
                react_cells(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut chemistry,
                    &mut cancelled,
                )?;
                let mechanical = self.trial.step_with_piston(
                    &mut trial_piston,
                    Time::new(duration),
                    controls.transport,
                    load,
                    piston_controls,
                    &mut cancelled,
                )?;
                if mechanical.gas.elapsed.value() != duration {
                    return Ok((mechanical.gas.elapsed.value(), None));
                }
                react_cells(
                    &mut self.trial,
                    &mut self.reactor,
                    &mut self.amounts,
                    0.5 * duration,
                    &mut chemistry,
                    &mut cancelled,
                )?;
                let flow = coupled_report(
                    &self.trial.gas,
                    before,
                    self.trial.totals(),
                    mechanical.gas,
                    chemistry,
                    None,
                    retry,
                )?;
                // This original reactive residual is measured from endpoint
                // species and energy differences, with the boundary energy
                // subtracted once. Restoring it gives actual chemical Delta E.
                let gas_energy_change_j =
                    flow.absolute_energy_residual_j + flow.transport.boundary_change.energy;
                let total_energy_residual_j = (gas_energy_change_j
                    + mechanical.kinetic_energy_change_j)
                    + mechanical.exterior_pressure_work_j
                    - mechanical.external_force_work_j;
                let external_impulse =
                    duration * (load.external_force_n - area * load.exterior_pressure_pa);
                let combined_momentum_residual_n_s = (flow.after.momentum - flow.before.momentum)
                    + original_piston.mass_kg()
                        * (trial_piston.velocity_m_s() - original_piston.velocity_m_s())
                    + mechanical.fixed_wall_impulse_n_s
                    - external_impulse;
                if ![
                    gas_energy_change_j,
                    total_energy_residual_j,
                    combined_momentum_residual_n_s,
                ]
                .iter()
                .all(|x| x.is_finite())
                {
                    return Err(ReactiveError::Invalid(
                        "non-finite complete chemical/piston balance",
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
                    Some(ReactivePistonStep {
                        flow,
                        piston: mechanical,
                        gas_energy_change_j,
                        combined_momentum_residual_n_s,
                        total_energy_residual_j,
                    }),
                ))
            })();
            match attempt {
                Ok((shorter, None)) => {
                    duration = shorter;
                    last =
                        ReactiveError::Invalid("piston shortened the complete reacting duration");
                }
                Ok((_, Some(report))) => {
                    poll(&mut cancelled)?;
                    core::mem::swap(&mut self.accepted, &mut self.trial);
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
