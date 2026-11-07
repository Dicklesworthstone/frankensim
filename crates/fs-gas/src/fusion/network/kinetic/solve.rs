use super::super::solve as burn;
use super::*;
use crate::fusion::emission::{self, NeutronMoments};
use crate::fusion::kinetic::{AlphaCollisionOperator, PreparedAlphaEscape};
use crate::fusion::solve::Sum;
use crate::fusion::volume::{adiabatic_work_on_pools, log_volume_ratio};
use fs_math::det::{exp, pow};
use fs_qty::{Temperature, Volume};
use fs_thermochem::coulomb::test_particle::PreparedTestParticleCoulombV1;
use fs_thermochem::coulomb::{CoulombErrorV1, CoulombScreeningV1};
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateErrorV1;
use fs_thermochem::fusion::plasma::DtFusionPlasmaErrorV1;
use fs_thermochem::plasma_radiation::BremsstrahlungErrorV1;

const KB: f64 = 1.380_649e-23;

#[derive(Clone, Copy)]
pub(super) struct Point {
    core: burn::Point,
    alphas: [f64; MAX_ALPHA_BINS],
    alpha_energy_scale: f64,
    alpha_boundary_work: Sum,
}
struct Sources {
    core: burn::Sources,
    collision: AlphaCollisionOperator,
    escape: Option<PreparedAlphaEscape>,
}
pub(super) struct Problem<'a> {
    initial: &'a KineticFusionState,
    core: burn::Problem<'a>,
    thermal: &'a FusionNetworkTransport,
    alpha: &'a KineticAlphaTransport,
    volume_history: Option<&'a FusionVolumeHistory>,
    controls: FusionNetworkSolve,
}
impl<'a> Problem<'a> {
    pub(super) fn new(
        initial: &'a KineticFusionState,
        channels: FusionNetworkChannels,
        heating: FusionNetworkHeating,
        thermal: &'a FusionNetworkTransport,
        alpha: &'a KineticAlphaTransport,
        volume_history: Option<&'a FusionVolumeHistory>,
        controls: FusionNetworkSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FusionNetworkError> {
        let source = alpha.collisions();
        if source.background().initial_state().thermodynamics() != initial.plasma().thermodynamics()
            || source.projectile_index() != initial.plasma().indices()[4]
            || source.background().screening() != CoulombScreeningV1::ElectronDebye
            || thermal
                .coulomb
                .as_ref()
                .is_some_and(|m| m.screening() != CoulombScreeningV1::ElectronDebye)
            || initial
                .energy_grid()
                .energies()
                .last()
                .is_none_or(|e| *e > source.maximum_energy())
        {
            return Err(FusionNetworkError::Invalid(
                "original He4 collision source, ElectronDebye screening and nonrelativistic alpha grid required",
            ));
        }
        let core = burn::Problem::new(
            &initial.reference,
            initial.plasma(),
            channels,
            heating,
            thermal,
            volume_history,
            controls,
            cancelled,
        )?
        .with_kinetic_alphas(
            initial.alpha_amount().value(),
            initial.thermal_helium4_amount,
            [
                initial.thermal_heavy_energy().value(),
                initial.thermal_electron_energy().value(),
            ],
            cancelled,
        )?;
        let problem = Self {
            initial,
            core,
            thermal,
            alpha,
            volume_history,
            controls,
        };
        problem.sources(problem.start(), cancelled)?;
        Ok(problem)
    }

    fn start(&self) -> Point {
        Point {
            core: self.core.start(),
            alphas: self.initial.alpha_amounts,
            alpha_energy_scale: 1.0,
            alpha_boundary_work: Sum::default(),
        }
    }
    fn moments(&self, point: &Point) -> (f64, f64) {
        state::moments_scaled(
            self.initial.energy_grid(),
            &point.alphas,
            point.alpha_energy_scale,
        )
    }
    fn admit(
        &self,
        point: &Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FusionNetworkError> {
        let temperatures = self.core.admit(point.core, cancelled)?;
        let grid = self.initial.energy_grid();
        grid.admit_scale(point.alpha_energy_scale)?;
        if grid.energies()[grid.len() - 1].value() * point.alpha_energy_scale
            > self.alpha.collisions().maximum_energy().value()
        {
            return Err(FusionNetworkError::OutsideDomain);
        }
        if point.alphas.iter().any(|n| !n.is_finite() || *n < 0.0)
            || point.alphas[grid.len()..].iter().any(|n| *n != 0.0)
            || !point.alpha_boundary_work.value.is_finite()
        {
            return Err(FusionNetworkError::Arithmetic(
                "finite nonnegative canonical alpha distribution",
            ));
        }
        let (number, energy) = self.moments(point);
        let thermal_amounts = self.core.thermal_amounts(point.core);
        let thermal_ions: f64 = self
            .initial
            .plasma()
            .indices()
            .iter()
            .map(|&i| thermal_amounts[i].value())
            .sum();
        let thermal_energy: f64 = self.core.pool_energy(point.core).iter().sum();
        if !number.is_finite()
            || !energy.is_finite()
            || number > self.alpha.maximum_number_ratio() * thermal_ions
            || energy > self.alpha.maximum_energy_ratio() * thermal_energy
            || (number > 0.0
                && point.alphas[grid.len() - 1] > self.alpha.maximum_upper_bin_fraction() * number)
        {
            return Err(FusionNetworkError::OutsideDomain);
        }
        // The radiation source sees all ion charges; explicitly check the fast
        // population against the stationary-ion approximation it otherwise
        // admits using only the Maxwellian heavy temperature.
        if number > 0.0 {
            if let Some(radiation) = &self.thermal.bremsstrahlung {
                let maximum_energy =
                    grid.energies()[grid.len() - 1].value() * point.alpha_energy_scale;
                let ion_variance =
                    (2.0 / 3.0) * maximum_energy / self.alpha.collisions().projectile_mass_kg();
                let electron_mass = self.initial.plasma().thermodynamics().models()
                    [self.initial.plasma().electron_index()]
                .molar_mass()
                .value()
                    / AVOGADRO;
                let electron_variance = KB * temperatures[1] / electron_mass;
                if ion_variance / electron_variance
                    > radiation
                        .controls()
                        .maximum_ion_to_electron_thermal_variance_ratio
                {
                    return Err(FusionNetworkError::Radiation(
                        BremsstrahlungErrorV1::Domain(
                            "kinetic alpha motion exceeds stationary-ion radiation limit",
                        ),
                    ));
                }
            }
        }
        Ok(())
    }

    fn sources(
        &self,
        point: Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Sources, FusionNetworkError> {
        let (core, prepared) = self.physical_sources(point, cancelled)?;
        let collision = AlphaCollisionOperator::new_scaled(
            self.initial.energy_grid(),
            point.alpha_energy_scale,
            &prepared,
            &mut *cancelled,
        )?;
        let escape = self
            .alpha
            .escape()
            .map(|table| {
                table.prepare(
                    self.initial.energy_grid(),
                    point.alpha_energy_scale,
                    cancelled,
                )
            })
            .transpose()?;
        Ok(Sources {
            core,
            collision,
            escape,
        })
    }

    /// Admit complete source coefficients at accepted endpoints and authored
    /// knots without rebuilding their energy-grid collision matrix. Midpoint
    /// samples alone can miss a source-domain excursion at a volume turnaround.
    fn physical_sources(
        &self,
        point: Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(burn::Sources, PreparedTestParticleCoulombV1), FusionNetworkError> {
        self.admit(&point, cancelled)?;
        let [th, te] = self.core.temperatures(point.core)?.map(Temperature::new);
        let prepared = self
            .alpha
            .collisions()
            .prepare_point(
                &self.core.thermal_amounts(point.core),
                point.core.volume,
                th,
                te,
                &mut *cancelled,
            )
            .map_err(coulomb_error)?;
        let core = self.core.sources(point.core, cancelled)?;
        let energies = self.initial.energy_grid().energies();
        if let Some(table) = self.alpha.escape() {
            // The original table covers a contiguous interval. Admit both
            // physical endpoints also at accepted geometry knots, where
            // interior midpoint sampling alone can miss a domain excursion.
            for endpoint in [energies[0], energies[energies.len() - 1]] {
                table.hazard_at(Energy::new(endpoint.value() * point.alpha_energy_scale))?;
            }
        }
        for channel in [0, 3] {
            if let Some(event) = core.events[channel] {
                let birth = event.product_kinetic_energies_j()[0];
                if birth < energies[0].value() * point.alpha_energy_scale
                    || birth > energies[energies.len() - 1].value() * point.alpha_energy_scale
                {
                    // A trial's changing temperature can move a birth beyond
                    // the authored support. Permit a shorter physical trial;
                    // initial support is still admitted before any evolution.
                    return Err(FusionNetworkError::OutsideDomain);
                }
            }
        }
        Ok((core, prepared))
    }

    fn collide(
        &self,
        old: Point,
        dt: f64,
        source: &AlphaCollisionOperator,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, FusionNetworkError> {
        let step = source.advance(
            &old.alphas,
            Time::new(dt),
            self.alpha.maximum_matrix_applications(),
            &mut *cancelled,
        )?;
        let mut next = old;
        next.alphas = step.amounts;
        self.core
            .deposit_alpha_energy(&mut next.core, step.deposited_j);
        self.admit(&next, cancelled)?;
        Ok(next)
    }

    fn frozen_step(
        &self,
        old: Point,
        dt: f64,
        source: &Sources,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, FusionNetworkError> {
        if let Some(escape) = &source.escape {
            let has_alphas_or_births = old.core.kinetic_alpha_amount.value > 0.0
                || source.core.events[0].is_some()
                || source.core.events[3].is_some();
            // Exact separate maps do not make a stiff composition exact:
            // saturated coarse/fine removals can both erase a continuously
            // replenished alpha stock and hide the splitting error. Resolve at
            // most one escape e-fold per coupled frozen interval. The original
            // bounded retry owner shortens the whole physical step; a source
            // too stiff for its declared budget refuses rather than certifying
            // that false equilibrium. An empty alpha-free solve has no split.
            if has_alphas_or_births && escape.maximum_rate_s_inv() * dt > 1.0 {
                return Err(FusionNetworkError::StepBudget);
            }
        }
        let escaped = self.escape(old, 0.5 * dt, source.escape.as_ref(), cancelled)?;
        let mut next = self.collide(escaped, 0.5 * dt, &source.collision, cancelled)?;
        let (core, births) =
            self.core
                .frozen_step_with_births(next.core, dt, source.core, cancelled)?;
        next.core = core;
        for (amount, channel) in births.into_iter().zip([0, 3]) {
            if amount == 0.0 {
                continue;
            }
            let event = source.core.events[channel].ok_or(FusionNetworkError::Arithmetic(
                "alpha birth without selected channel",
            ))?;
            self.initial.energy_grid().inject_scaled(
                &mut next.alphas,
                Amount::new(amount),
                Energy::new(event.product_kinetic_energies_j()[0]),
                next.alpha_energy_scale,
            )?;
        }
        self.admit(&next, cancelled)?;
        let next = self.collide(next, 0.5 * dt, &source.collision, cancelled)?;
        self.escape(next, 0.5 * dt, source.escape.as_ref(), cancelled)
    }

    fn escape(
        &self,
        old: Point,
        dt: f64,
        source: Option<&PreparedAlphaEscape>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, FusionNetworkError> {
        let Some(source) = source else {
            return Ok(old);
        };
        let step = source.advance(
            self.initial.energy_grid(),
            old.alpha_energy_scale,
            &old.alphas,
            dt,
            cancelled,
        )?;
        let mut next = old;
        next.alphas = step.amounts;
        next.core = self.core.escape_alphas(
            old.core,
            step.removed_mol,
            step.retained_mol,
            step.removed_energy_j,
            cancelled,
        )?;
        self.admit(&next, cancelled)?;
        Ok(next)
    }

    fn midpoint(
        &self,
        old: Point,
        end: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, FusionNetworkError> {
        let dt = end - old.core.time;
        let middle_time = old.core.time + 0.5 * dt;
        if middle_time <= old.core.time || middle_time >= end {
            return Err(FusionNetworkError::StepBudget);
        }
        let prepared = if let Some(history) = self.volume_history {
            self.volume_step(old, history.volume_at(Time::new(middle_time))?, cancelled)?
        } else {
            old
        };
        let first = self.sources(prepared, cancelled)?;
        let predicted = self.frozen_step(prepared, 0.5 * dt, &first, cancelled)?;
        let middle = self.sources(predicted, cancelled)?;
        let mut next = self.frozen_step(prepared, dt, &middle, cancelled)?;
        next.core.time = end;
        if let Some(history) = self.volume_history {
            self.volume_step(next, history.volume_at(Time::new(end))?, cancelled)
        } else {
            Ok(next)
        }
    }

    /// Exact isotropic adiabatic motion in a comoving energy coordinate. No
    /// remap, collision, reaction, or source clock is applied by this operation.
    fn volume_step(
        &self,
        old: Point,
        volume: Volume,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, FusionNetworkError> {
        poll(cancelled)?;
        if old.core.volume == volume {
            return Ok(old);
        }
        let energy = self.moments(&old).1;
        let work = adiabatic_work_on_pools([energy], old.core.volume, volume)?[0];
        let mut next = old;
        next.core = self.core.volume_step(old.core, volume, cancelled)?;
        // Anchor every stage to the same initial coordinate. Full and refined
        // trials at an identical volume therefore compare identical bins;
        // returning to the reference volume restores its exact energy scale.
        next.alpha_energy_scale =
            exp((-2.0 / 3.0) * log_volume_ratio(self.initial.volume(), volume)?);
        next.alpha_boundary_work.add(work);
        self.admit(&next, cancelled)?;
        Ok(next)
    }

    fn error(&self, old: Point, full: Point, fine: Point) -> Result<f64, FusionNetworkError> {
        let c = self.controls;
        let mut error = self.core.error(old.core, full.core, fine.core)?;
        for (&a, &b) in full.alphas[..self.initial.energy_grid().len()]
            .iter()
            .zip(&fine.alphas)
        {
            error = error.max(
                (a - b).abs()
                    / (3.0
                        * (c.absolute_amount_tolerance_mol
                            + c.relative_tolerance * a.abs().max(b.abs()))),
            );
        }
        let (a, b) = (self.moments(&full).1, self.moments(&fine).1);
        error = error.max(
            (a - b).abs()
                / (3.0 * (c.energy_tolerance_j + c.relative_tolerance * a.abs().max(b.abs()))),
        );
        let a = full.alpha_boundary_work.value - old.alpha_boundary_work.value;
        let b = fine.alpha_boundary_work.value - old.alpha_boundary_work.value;
        error = error.max(
            (a - b).abs()
                / (3.0 * (c.energy_tolerance_j + c.relative_tolerance * a.abs().max(b.abs()))),
        );
        if !error.is_finite() {
            return Err(FusionNetworkError::Arithmetic(
                "kinetic network local error",
            ));
        }
        Ok(error)
    }

    pub(super) fn integrate(
        &self,
        duration: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Point, FusionNetworkBalance, emission::NeutronHistory), FusionNetworkError> {
        let mut point = self.start();
        let mut clock = Sum::default();
        let mut balance = FusionNetworkBalance::default();
        let mut neutron_history = self.controls.record_neutron_history.then(Vec::new);
        let mut step = self.controls.initial_step_s;
        for _ in 0..self.controls.maximum_steps {
            poll(cancelled)?;
            if point.core.time == duration {
                return Ok((point, balance, neutron_history));
            }
            let next_boundary = if let Some(history) = self.volume_history {
                history
                    .next_knot_after(Time::new(point.core.time))?
                    .map_or(duration, |time| time.value())
            } else {
                duration
            };
            let mut next_clock = clock;
            next_clock.add(step.min(self.controls.maximum_step_s));
            let end = next_clock.value.min(next_boundary);
            if end == next_boundary {
                next_clock = sum(next_boundary);
            }
            let dt = end - point.core.time;
            let middle_time = point.core.time + 0.5 * dt;
            if middle_time <= point.core.time || middle_time >= end {
                return Err(FusionNetworkError::StepBudget);
            }
            let attempt = (|| {
                let mut start = point;
                if let Some(moments) = &mut start.core.neutron_interval {
                    *moments = NeutronMoments::default();
                }
                let full = self.midpoint(start, end, cancelled)?;
                let half = self.midpoint(start, middle_time, cancelled)?;
                let fine = self.midpoint(half, end, cancelled)?;
                self.physical_sources(fine, cancelled)?;
                Ok::<_, FusionNetworkError>((fine, self.error(point, full, fine)?))
            })();
            let error = match attempt {
                Ok((next, error)) if error <= 1.0 => {
                    emission::record_interval(
                        &mut neutron_history,
                        point.core.time,
                        next.core.time,
                        next.core.neutron_interval,
                    )?;
                    point = next;
                    clock = next_clock;
                    balance.accepted_steps += 1;
                    balance.maximum_normalized_local_error =
                        balance.maximum_normalized_local_error.max(error);
                    let factor = if error == 0.0 {
                        2.0
                    } else {
                        (0.9 / pow(error, 1.0 / 3.0)).clamp(0.5, 2.0)
                    };
                    step = (dt * factor).min(self.controls.maximum_step_s);
                    continue;
                }
                Ok((_, error)) => Some(error),
                Err(FusionNetworkError::StepBudget) => {
                    // Frozen collision work scales with interval exposure.
                    // Retry within both the per-call and network work caps.
                    balance.rejected_steps += 1;
                    step = 0.5 * dt;
                    if step < self.controls.minimum_step_s {
                        return Err(FusionNetworkError::StepBudget);
                    }
                    continue;
                }
                Err(
                    FusionNetworkError::OutsideDomain
                    | FusionNetworkError::Plasma(DtFusionPlasmaErrorV1::Domain(_))
                    | FusionNetworkError::Coulomb(CoulombErrorV1::Domain(_))
                    | FusionNetworkError::Radiation(BremsstrahlungErrorV1::Domain(_)),
                ) => None,
                Err(error) => return Err(error),
            };
            balance.rejected_steps += 1;
            step = dt * error.map_or(0.5, |e| (0.8 / pow(e, 1.0 / 3.0)).clamp(0.1, 0.5));
            if step < self.controls.minimum_step_s {
                return Err(if error.is_some() {
                    FusionNetworkError::StepBudget
                } else {
                    FusionNetworkError::OutsideDomain
                });
            }
        }
        if point.core.time == duration {
            Ok((point, balance, neutron_history))
        } else {
            Err(FusionNetworkError::StepBudget)
        }
    }

    pub(super) fn publish(
        &self,
        point: Point,
        balance: &mut FusionNetworkBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(KineticFusionState, KineticAlphaBalance), FusionNetworkError> {
        self.sources(point, cancelled)?;
        let [th, te] = self.core.temperatures(point.core)?.map(Temperature::new);
        let reference = FrozenIonizedGasStateV1::new(
            self.initial.plasma().thermodynamics(),
            &self.core.amounts(point.core),
            point.core.volume,
            th,
            te,
            64.0 * f64::EPSILON,
            &mut *cancelled,
        )
        .map_err(|e| match e {
            FrozenIonizedGasStateErrorV1::Cancelled => FusionNetworkError::Cancelled,
            e => FusionNetworkError::State(e),
        })?;
        let physical_grid = self
            .initial
            .energy_grid()
            .scaled(point.alpha_energy_scale, &mut *cancelled)?;
        let state = KineticFusionState::from_parts(
            &reference,
            self.initial.plasma(),
            &physical_grid,
            point.alphas,
            point.core.thermal_helium4_amount.value,
            point.core.kinetic_alpha_amount.value,
            Time::new(self.initial.elapsed_time().value() + point.core.time),
            self.controls.plasma_limits,
            cancelled,
        )?;
        let energy_change = state.alpha_energy().value() - self.initial.alpha_energy().value();
        self.core.balance_at(
            point.core,
            [
                state.thermal_heavy_energy().value(),
                state.thermal_electron_energy().value(),
            ],
            energy_change,
            point.alpha_boundary_work.value,
            balance,
            cancelled,
        )?;
        let born_amount = balance.reaction_extents_mol[0] + balance.reaction_extents_mol[3];
        let born_energy =
            balance.product_kinetic_energies_j[0][0] + balance.product_kinetic_energies_j[3][0];
        let deposited = point.core.alpha_deposited.map(|v| v.value);
        let (number, _) = self.moments(&point);
        let amount_residual = number - point.core.kinetic_alpha_amount.value;
        let energy_residual =
            energy_change + deposited[0] + deposited[1] + point.core.escaped_alpha_energy.value
                - born_energy
                - point.alpha_boundary_work.value;
        let c = self.controls;
        if !energy_residual.is_finite()
            || !amount_residual.is_finite()
            || energy_residual.abs()
                > c.energy_tolerance_j
                    + c.balance_relative_tolerance
                        * (born_energy
                            + self.initial.alpha_energy().value()
                            + state.alpha_energy().value()
                            + deposited[0].abs()
                            + deposited[1].abs()
                            + point.core.escaped_alpha_energy.value
                            + point.alpha_boundary_work.value.abs())
            || amount_residual.abs()
                > c.absolute_amount_tolerance_mol
                    + c.balance_relative_tolerance
                        * number.max(point.core.kinetic_alpha_amount.value)
        {
            return Err(FusionNetworkError::Balance(
                "kinetic alpha distribution, births and reciprocal heat",
            ));
        }
        let alpha_balance = KineticAlphaBalance {
            born_alpha_amount_mol: born_amount,
            born_alpha_energy_j: born_energy,
            escaped_alpha_amount_mol: point.core.escaped_alpha_amount.value,
            escaped_alpha_energy_j: point.core.escaped_alpha_energy.value,
            escaped_electrons_mol: 2.0 * point.core.escaped_alpha_amount.value,
            escaped_electron_thermal_energy_j: point.core.escaped_alpha_electron_energy.value,
            deposited_heavy_energy_j: deposited[0],
            deposited_electron_energy_j: deposited[1],
            alpha_energy_change_j: energy_change,
            boundary_work_on_alphas_j: point.alpha_boundary_work.value,
            alpha_energy_residual_j: energy_residual,
            alpha_amount_residual_mol: amount_residual,
            upper_bin_fraction: if number == 0.0 {
                0.0
            } else {
                point.alphas[self.initial.energy_grid().len() - 1] / number
            },
        };
        poll(cancelled)?;
        Ok((state, alpha_balance))
    }
}
fn sum(value: f64) -> Sum {
    let mut result = Sum::default();
    result.add(value);
    result
}
fn coulomb_error(e: CoulombErrorV1) -> FusionNetworkError {
    match e {
        CoulombErrorV1::Cancelled => FusionNetworkError::Cancelled,
        e => FusionNetworkError::Coulomb(e),
    }
}
