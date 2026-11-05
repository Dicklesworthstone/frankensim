//! Chemistry, resolved compression and heat through one finite piston store.

use super::*;
use crate::multicomponent::wall_exchange::WallControl;
use crate::multicomponent::{
    SharedPistonThermal, SharedPistonThermalError, SharedPistonThermalStepReport, Totals,
};

/// Complete physical accounts for two reacting chambers and one thermal piston.
///
/// The nested piston report covers H/2, pressure/motion, H/2 between the two
/// chemistry halves. These outer accounts additionally include the original
/// mechanism in each chamber. Heat, mechanical work and species-reference
/// changes retain separate entries and each enters the balance exactly once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveSharedPistonThermalStep<P = PistonStep<SharedPistonThermalStepReport>> {
    /// Common elapsed duration of pressure, motion and both chemical/heat halves.
    pub elapsed: Time,
    /// Complete reacting compositions restarted from the accepted owners.
    pub coupling_retries: usize,
    /// Independent original-source chemistry in [left, right] chamber order.
    pub chemistry: [ChemistryReport; 2],
    /// Native resolved pressure, wall heat and mechanical response.
    pub piston: P,
    /// Complete chamber inventories before any chemical, thermal or mechanical work.
    pub before: SharedPistonGasReport<Totals>,
    /// Complete chamber inventories after all stages.
    pub after: SharedPistonGasReport<Totals>,
    /// Species, momentum and referenced-energy residuals after pressure transfer,
    /// the original chemistry changes and separately reported conductive heat.
    pub conservation_residual: SharedPistonGasReport<Totals>,
    /// Total mass changes after external transfer in [left, right] order [kg].
    pub mass_residual_kg: [f64; 2],
    /// Actual physical chemical-energy changes [J], including each chamber's
    /// OWN source-species reference changes, in [left, right] order.
    pub gas_energy_change_j: [f64; 2],
    /// Physical chemical-energy change minus pressure work and heat [J].
    /// These differences avoid subtracting large absolute formation energies.
    pub absolute_energy_residual_j: [f64; 2],
    /// Both gas momentum changes, piston momentum and fixed-support reactions,
    /// minus the applied and suspension impulses [N s].
    pub combined_momentum_residual_n_s: f64,
    /// Both chemical gas energies, ONE wall store, piston kinetic and spring
    /// storage, plus external damping loss, minus authored drive work [J].
    pub total_energy_residual_j: f64,
}

impl<P> ReactiveSharedPistonThermalStep<P> {
    fn map_piston<Q>(self, map: impl FnOnce(P) -> Q) -> ReactiveSharedPistonThermalStep<Q> {
        ReactiveSharedPistonThermalStep {
            elapsed: self.elapsed,
            coupling_retries: self.coupling_retries,
            chemistry: self.chemistry,
            piston: map(self.piston),
            before: self.before,
            after: self.after,
            conservation_residual: self.conservation_residual,
            mass_residual_kg: self.mass_residual_kg,
            gas_energy_change_j: self.gas_energy_change_j,
            absolute_energy_residual_j: self.absolute_energy_residual_j,
            combined_momentum_residual_n_s: self.combined_momentum_residual_n_s,
            total_energy_residual_j: self.total_energy_residual_j,
        }
    }
}

impl From<SharedPistonThermalError> for ReactiveError {
    fn from(error: SharedPistonThermalError) -> Self {
        match error {
            SharedPistonThermalError::Piston(error) => error.into(),
            SharedPistonThermalError::Heat(error) => error.into(),
            SharedPistonThermalError::Cancelled => Self::Cancelled,
            SharedPistonThermalError::StepRetries { attempts, last } => Self::RetryLimit {
                attempts,
                last: Box::new((*last).into()),
            },
        }
    }
}

impl ReactingEuler1d {
    /// Couple both original chemical mechanisms to heat through one finite piston.
    ///
    /// The complete composition is C/2, H/2, pressure/motion, H/2, C/2.
    /// The first chemistry and heat stages use original physical cell volumes;
    /// the second stages use actual moved volumes. Both heat contacts share
    /// the SAME finite inventory. Its original calorics and two authored hA
    /// values supply heat; no gas heat capacity, contact area or wall mass is
    /// inferred from piston inertia or chamber compression.
    ///
    /// Chemistry changes the pressure seen by the joint midpoint-speed solve.
    /// Heat changes both actual end-cell pressures before that solve; the second
    /// chemical half sees the moved, thermally updated states. Each chamber may
    /// have a different original species bank, mechanism and formation gauge.
    /// Chemical energy is converted only by its original fixed-U,V reactor.
    ///
    /// A shorter inner thermal/mechanical interval restarts BOTH chemical halves
    /// and BOTH heat halves from accepted storage at the common duration. No
    /// intermediate wall temperature or gas state publishes. Fields, meshes,
    /// clocks, piston, finite wall and heat report publish after complete energy,
    /// momentum and source admission and the final cancellation checkpoint.
    ///
    /// This retains the existing closed constant-area chambers, one uniformly
    /// heated translating wall and first-order implicit heat stages. It adds no
    /// molecular transport, contact correction, leakage, wall gradients, source
    /// certification or general temporal-order/entropy-stability claim.
    ///
    /// # Errors
    /// Original source, heat, geometry and mechanical refusals, invalid controls,
    /// exhausted complete retries or cancellation leave every accepted owner unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_piston_and_heat(
        &mut self,
        right: &mut ReactingEuler1d,
        thermal: &mut SharedPistonThermal<'_>,
        piston: &mut Piston,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        heat_controls: WallControl,
        external_force_n: f64,
        piston_controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveSharedPistonThermalStep, ReactiveError> {
        self.step_with_shared_piston_thermal_suspension(
            right,
            thermal,
            piston,
            maximum_duration,
            controls,
            heat_controls,
            external_force_n,
            None,
            0.0,
            piston_controls,
            cancelled,
        )
        .map(|step| step.map_piston(|mechanical| mechanical.piston))
    }

    /// Add a held spring/damper and an explicit destination for its dissipated heat.
    ///
    /// The original suspension enters the same pressure/motion root. The authored
    /// fraction in [0,1] of its actual accepted loss is deposited into the ONE
    /// wall before H2; the remaining loss leaves this thermal system. H2 and C2
    /// therefore see the returned heat at their actual moved cell volumes.
    /// The whole balance includes spring storage and only external damping loss.
    /// Native mechanical reports still expose the full damper loss separately.
    ///
    /// A failed chemistry, thermal, pressure or balance stage pays no lasting
    /// damper heat. The wall's accepted storage and report remain unchanged until
    /// both chamber mechanisms admit the complete candidate. The suspension's
    /// stiffness and anchor are fixed for this step; changes between calls need
    /// their own caller work account.
    ///
    /// # Errors
    /// The same source, thermal, mechanical and cancellation refusals as the free
    /// thermal piston apply. Nonfinite or out-of-range returned-heat fractions refuse.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_shared_suspended_piston_and_heat(
        &mut self,
        right: &mut ReactingEuler1d,
        thermal: &mut SharedPistonThermal<'_>,
        piston: &mut Piston,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        heat_controls: WallControl,
        external_force_n: f64,
        suspension: PistonSuspension,
        damping_heat_fraction_to_wall: f64,
        piston_controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        ReactiveSharedPistonThermalStep<SuspendedPistonStep<SharedPistonThermalStepReport>>,
        ReactiveError,
    > {
        self.step_with_shared_piston_thermal_suspension(
            right,
            thermal,
            piston,
            maximum_duration,
            controls,
            heat_controls,
            external_force_n,
            Some(suspension),
            damping_heat_fraction_to_wall,
            piston_controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with_shared_piston_thermal_suspension(
        &mut self,
        right: &mut ReactingEuler1d,
        thermal: &mut SharedPistonThermal<'_>,
        piston: &mut Piston,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        heat_controls: WallControl,
        external_force_n: f64,
        suspension: Option<PistonSuspension>,
        damping_heat_fraction_to_wall: f64,
        piston_controls: PistonControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        ReactiveSharedPistonThermalStep<SuspendedPistonStep<SharedPistonThermalStepReport>>,
        ReactiveError,
    > {
        poll(&mut cancelled)?;
        controls.validate()?;
        heat_controls.validate()?;
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
        if !damping_heat_fraction_to_wall.is_finite()
            || !(0.0..=1.0).contains(&damping_heat_fraction_to_wall)
        {
            return Err(ReactiveError::Invalid(
                "shared reacting piston damping heat fraction must be between zero and one",
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
        if let Some(law) = suspension {
            law.bind(self.accepted.mesh.positions()[self.accepted.cell_count()])?;
        }
        let before = SharedPistonGasReport {
            left: self.accepted.totals(),
            right: right.accepted.totals(),
        };
        let original_piston = *piston;
        let mut last =
            ReactiveError::Invalid("no complete shared reacting thermal piston increment");
        for retry in 0..=controls.maximum_coupling_retries {
            poll(&mut cancelled)?;
            let time = self.accepted.time;
            if !positive(duration)
                || !positive(0.5 * duration)
                || !positive(time + duration - time)
                || !(time + duration).is_finite()
            {
                return Err(ReactiveError::Invalid(
                    "shared reacting thermal piston cannot advance both clocks",
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
                let mechanical = thermal.advance_trial(
                    &mut self.trial,
                    &mut right.trial,
                    &mut trial_piston,
                    Time::new(duration),
                    controls.transport,
                    heat_controls,
                    external_force_n,
                    suspension,
                    damping_heat_fraction_to_wall,
                    piston_controls,
                    None,
                    &mut cancelled,
                )?;
                let transport = mechanical.piston.gas.transport;
                let elapsed = transport.left.elapsed;
                if elapsed != transport.right.elapsed {
                    return Err(ReactiveError::Invalid(
                        "shared thermal pressure solve returned unequal durations",
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
                let after = SharedPistonGasReport {
                    left: self.trial.totals(),
                    right: right.trial.totals(),
                };
                let sides = [
                    coupled_report(
                        &self.trial.gas,
                        before.left,
                        after.left,
                        transport.left,
                        chemistry[0],
                        None,
                        retry,
                    )?,
                    coupled_report(
                        &right.trial.gas,
                        before.right,
                        after.right,
                        transport.right,
                        chemistry[1],
                        None,
                        retry,
                    )?,
                ];
                let heat = mechanical.piston.gas.heat;
                // The original chamber balance has no conductive source input.
                // Its physical difference still contains this chamber's heat;
                // remove that once in the complete thermal balance, while keeping
                // its original transport and mechanism reports unmodified.
                let gas_energy_change_j = sides.map(|side| {
                    sum([
                        side.absolute_energy_residual_j,
                        side.transport.boundary_change.energy,
                    ])
                });
                let absolute_energy_residual_j = [0, 1]
                    .map(|j| sum([sides[j].absolute_energy_residual_j, -heat.heat_to_gas_j[j]]));
                let mut residual = sides.map(|side| side.conservation_residual);
                for (j, r) in residual.iter_mut().enumerate() {
                    r.energy = sum([r.energy, -heat.heat_to_gas_j[j]]);
                }
                let mass_residual_kg = sides.map(|side| side.mass_residual_kg);
                let total_energy_residual_j = sum([
                    gas_energy_change_j[0],
                    gas_energy_change_j[1],
                    heat.wall_energy_change_j,
                    mechanical.piston.kinetic_energy_change_j,
                    mechanical.spring_energy_change_j,
                    heat.damping_loss_external_j,
                    -mechanical.applied_force_work_j,
                ]);
                let combined_momentum_residual_n_s = sum([
                    after.left.momentum - before.left.momentum,
                    after.right.momentum - before.right.momentum,
                    original_piston.mass_kg()
                        * (trial_piston.velocity_m_s() - original_piston.velocity_m_s()),
                    mechanical.piston.fixed_wall_impulse_n_s,
                    -duration * external_force_n,
                    -mechanical.suspension_impulse_n_s,
                ]);
                if !gas_energy_change_j
                    .into_iter()
                    .chain(absolute_energy_residual_j)
                    .chain([total_energy_residual_j, combined_momentum_residual_n_s])
                    .all(f64::is_finite)
                {
                    return Err(ReactiveError::Invalid(
                        "non-finite complete shared chemical/thermal/piston balance",
                    ));
                }
                let energy_error = absolute_energy_residual_j[0]
                    .abs()
                    .max(absolute_energy_residual_j[1].abs())
                    .max(total_energy_residual_j.abs());
                if energy_error > piston_controls.energy_tolerance_j
                    || combined_momentum_residual_n_s.abs() > piston_controls.momentum_tolerance_n_s
                {
                    return Err(PistonError::Balance {
                        momentum_residual_n_s: combined_momentum_residual_n_s.abs(),
                        energy_residual_j: energy_error,
                    }
                    .into());
                }
                Ok::<_, ReactiveError>((
                    duration,
                    Some(ReactiveSharedPistonThermalStep {
                        elapsed,
                        coupling_retries: retry,
                        chemistry,
                        piston: mechanical,
                        before,
                        after,
                        conservation_residual: SharedPistonGasReport {
                            left: residual[0],
                            right: residual[1],
                        },
                        mass_residual_kg,
                        gas_energy_change_j,
                        absolute_energy_residual_j,
                        combined_momentum_residual_n_s,
                        total_energy_residual_j,
                    }),
                ))
            })();
            match attempt {
                Ok((shorter, None)) => {
                    duration = shorter;
                    last = ReactiveError::Invalid(
                        "heat/pressure shortened both reacting chamber intervals",
                    );
                }
                Ok((_, Some(report))) => {
                    poll(&mut cancelled)?;
                    // No fallible work follows: publish the one thermal inventory
                    // only when both final chemistry halves have been admitted.
                    thermal.publish_trial_heat(report.piston.piston.gas.heat);
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
