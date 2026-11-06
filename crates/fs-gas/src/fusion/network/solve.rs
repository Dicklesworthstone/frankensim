use super::super::solve::{Sum, frozen_depletion};
use super::super::thermal::{self, ThermalSources};
use super::super::volume::{adiabatic_work_on_pools, log_volume_ratio};
use super::*;
use fs_math::det::{exp, pow};
use fs_qty::{Temperature, Volume};
use fs_thermochem::coulomb::CoulombErrorV1;
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateErrorV1;
use fs_thermochem::fusion::FusionThermalEventV1;
use fs_thermochem::fusion::plasma::DtFusionPlasmaErrorV1;
use fs_thermochem::plasma_radiation::BremsstrahlungErrorV1;

#[derive(Clone, Copy)]
pub(super) struct Point {
    pub(super) time: f64,
    pub(super) volume: Volume,
    // Moving-volume states carry the positive thermal survivor directly.
    // Reconstructing it as U_initial + W can erase it after deep expansion.
    thermal_energy: [Sum; 2],
    work_on_pools: [Sum; 2],
    pub(super) kinetic_alpha_amount: Sum,
    pub(super) alpha_deposited: [Sum; 2],
    amounts: [Sum; 6],
    extents: [Sum; 4],
    products: [[Sum; 2]; 4],
    deposited: [Sum; 2], // heavy, electron
    ion_removal: Sum,
    neutrons: Sum,
    exchange: Sum,
    photons: Sum,
}

#[derive(Clone, Copy)]
pub(super) struct Sources {
    thermal: ThermalSources,
    pub(super) events: [Option<FusionThermalEventV1>; 4],
}

pub(super) struct Problem<'a> {
    initial: &'a FrozenIonizedGasStateV1,
    plasma: &'a LightIonFusionPlasmaV1,
    channels: FusionNetworkChannels,
    heating: FusionNetworkHeating,
    transport: &'a FusionNetworkTransport,
    volume_history: Option<&'a FusionVolumeHistory>,
    controls: FusionNetworkSolve,
    initial_amounts: [f64; 6],
    initial_energy: [f64; 2],
    q_molar: [f64; 4],
    kinetic_alphas: bool,
    initial_alpha_amount: f64,
    thermal_helium4_amount: f64,
}

impl<'a> Problem<'a> {
    pub(super) fn new(
        initial: &'a FrozenIonizedGasStateV1,
        plasma: &'a LightIonFusionPlasmaV1,
        channels: FusionNetworkChannels,
        heating: FusionNetworkHeating,
        transport: &'a FusionNetworkTransport,
        volume_history: Option<&'a FusionVolumeHistory>,
        controls: FusionNetworkSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FusionNetworkError> {
        plasma
            .admit_state_with_limits(initial, controls.plasma_limits, &mut *cancelled)
            .map_err(plasma_error)?;
        plasma
            .admit_temperature_interval(
                Temperature::new(controls.minimum_heavy_temperature_k),
                Temperature::new(controls.maximum_heavy_temperature_k),
                Temperature::new(controls.minimum_electron_temperature_k),
                Temperature::new(controls.maximum_electron_temperature_k),
                &mut *cancelled,
            )
            .map_err(plasma_error)?;
        for reaction in FusionReactionV1::ALL {
            if channels.contains(reaction) {
                for t in [
                    controls.minimum_heavy_temperature_k,
                    controls.maximum_heavy_temperature_k,
                ] {
                    reaction
                        .reactivity(Temperature::new(t))
                        .map_err(FusionNetworkError::Reaction)?;
                }
            }
        }
        if transport
            .coulomb
            .as_ref()
            .is_some_and(|m| m.initial_state().thermodynamics() != plasma.thermodynamics())
            || transport
                .bremsstrahlung
                .as_ref()
                .is_some_and(|m| m.thermodynamics() != plasma.thermodynamics())
        {
            return Err(FusionNetworkError::Invalid(
                "network thermal operators must retain the original source bank",
            ));
        }
        let result = Self {
            initial,
            plasma,
            channels,
            heating,
            transport,
            volume_history,
            controls,
            initial_amounts: core::array::from_fn(|i| initial.amounts()[i].value()),
            initial_energy: [
                initial.heavy_internal_energy().value(),
                initial.electron_internal_energy().value(),
            ],
            q_molar: FusionReactionV1::ALL.map(|r| r.q_j() * AVOGADRO),
            kinetic_alphas: false,
            initial_alpha_amount: 0.0,
            thermal_helium4_amount: initial.amounts()[plasma.indices()[4]].value(),
        };
        result.sources(result.start(), cancelled)?;
        Ok(result)
    }

    /// Activate the separate alpha inventory after the original bank has been
    /// admitted. The reference state supplies source identity, total inventories
    /// and temperatures; it is not published as the hybrid physical state.
    pub(super) fn with_kinetic_alphas(
        mut self,
        alpha_amount: f64,
        thermal_helium4_amount: f64,
        thermal_energy: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FusionNetworkError> {
        let helium4 = self.initial.amounts()[self.plasma.indices()[4]].value();
        if self.volume_history.is_some()
            || !alpha_amount.is_finite()
            || alpha_amount < 0.0
            || !thermal_helium4_amount.is_finite()
            || thermal_helium4_amount < 0.0
            || (helium4 - alpha_amount - thermal_helium4_amount).abs()
                > 64.0 * f64::EPSILON * helium4.max(alpha_amount)
            || thermal_energy.iter().any(|&u| !positive(u))
        {
            return Err(FusionNetworkError::Invalid(
                "rigid kinetic-alpha component and positive thermal stores",
            ));
        }
        self.kinetic_alphas = true;
        self.initial_alpha_amount = alpha_amount;
        self.thermal_helium4_amount = thermal_helium4_amount;
        self.initial_energy = thermal_energy;
        self.sources(self.start(), cancelled)?;
        Ok(self)
    }

    pub(super) fn start(&self) -> Point {
        Point {
            time: 0.0,
            volume: self.initial.volume(),
            thermal_energy: self.initial_energy.map(sum),
            work_on_pools: [Sum::default(); 2],
            kinetic_alpha_amount: sum(self.initial_alpha_amount),
            alpha_deposited: [Sum::default(); 2],
            amounts: self.initial_amounts.map(sum),
            extents: [Sum::default(); 4],
            products: [[Sum::default(); 2]; 4],
            deposited: [Sum::default(); 2],
            ion_removal: Sum::default(),
            neutrons: Sum::default(),
            exchange: Sum::default(),
            photons: Sum::default(),
        }
    }
    pub(super) fn amounts(&self, point: Point) -> [Amount; 6] {
        point.amounts.map(|n| Amount::new(n.value))
    }
    pub(super) fn thermal_amounts(&self, point: Point) -> [Amount; 6] {
        let mut amounts = self.amounts(point);
        if self.kinetic_alphas {
            // He4 is not a reactant. Retain the original thermal stock directly
            // rather than subtracting two large nearly equal alpha inventories.
            amounts[self.plasma.indices()[4]] = Amount::new(self.thermal_helium4_amount);
        }
        amounts
    }
    fn capacities(&self, point: Point) -> [f64; 2] {
        [
            CV_MOLAR
                * self
                    .plasma
                    .indices()
                    .iter()
                    .map(|&i| {
                        if self.kinetic_alphas && i == self.plasma.indices()[4] {
                            self.thermal_helium4_amount
                        } else {
                            point.amounts[i].value
                        }
                    })
                    .sum::<f64>(),
            CV_MOLAR * point.amounts[self.plasma.electron_index()].value,
        ]
    }
    pub(super) fn pool_energy(&self, point: Point) -> [f64; 2] {
        if self.volume_history.is_some() || self.kinetic_alphas {
            return point.thermal_energy.map(|x| x.value);
        }
        self.accounted_energy(point)
    }
    /// Independent cumulative transfer accounts, also used to verify moving
    /// states. The rigid path retains its original arithmetic and fast path.
    fn accounted_energy(&self, point: Point) -> [f64; 2] {
        [
            self.initial_energy[0]
                + self.heating.heavy_power.value() * point.time
                + (point.deposited[0].value - point.ion_removal.value)
                + point.exchange.value
                + point.work_on_pools[0].value
                + point.alpha_deposited[0].value,
            self.initial_energy[1]
                + self.heating.electron_power.value() * point.time
                + point.deposited[1].value
                - point.exchange.value
                - point.photons.value
                + point.work_on_pools[1].value
                + point.alpha_deposited[1].value,
        ]
    }
    fn add_pool_energy(&self, point: &mut Point, increments: [f64; 2]) {
        if self.volume_history.is_some() || self.kinetic_alphas {
            for (energy, increment) in point.thermal_energy.iter_mut().zip(increments) {
                energy.add(increment);
            }
        }
    }
    pub(super) fn temperatures(&self, point: Point) -> Result<[f64; 2], FusionNetworkError> {
        // Keep exact authored endpoints when no physical transfer has occurred.
        // A thermal substep can act at t=0, so elapsed time alone is insufficient.
        if point.extents.iter().all(|x| x.value == 0.0)
            && point.exchange.value == 0.0
            && point.photons.value == 0.0
            && point.work_on_pools.iter().all(|x| x.value == 0.0)
            && point.alpha_deposited.iter().all(|x| x.value == 0.0)
            && point.volume == self.initial.volume()
            && (point.time == 0.0
                || (self.heating.heavy_power.value() == 0.0
                    && self.heating.electron_power.value() == 0.0))
        {
            return Ok([
                self.initial.heavy_temperature().value(),
                self.initial.electron_temperature().value(),
            ]);
        }
        let capacity = self.capacities(point);
        let energy = self.pool_energy(point);
        let result = [energy[0] / capacity[0], energy[1] / capacity[1]];
        if result.iter().any(|&t| !positive(t)) {
            return Err(FusionNetworkError::OutsideDomain);
        }
        Ok(result)
    }
    pub(super) fn admit(
        &self,
        point: Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], FusionNetworkError> {
        poll(cancelled)?;
        if !point.time.is_finite()
            || point.time < 0.0
            || !point.exchange.value.is_finite()
            || !positive(point.volume.value())
            || point.work_on_pools.iter().any(|x| !x.value.is_finite())
            || point.alpha_deposited.iter().any(|x| !x.value.is_finite())
            || !point.kinetic_alpha_amount.value.is_finite()
            || point.kinetic_alpha_amount.value < 0.0
            || ((self.volume_history.is_some() || self.kinetic_alphas)
                && point.thermal_energy.iter().any(|x| !x.value.is_finite()))
            || point
                .amounts
                .iter()
                .chain(point.extents.iter())
                .chain(point.products.iter().flatten())
                .chain(point.deposited.iter())
                .chain([&point.ion_removal, &point.neutrons, &point.photons])
                .any(|x| !x.value.is_finite() || x.value < 0.0)
        {
            return Err(FusionNetworkError::Arithmetic(
                "finite nonnegative network stocks and accounts",
            ));
        }
        let [th, te] = self.temperatures(point)?;
        let c = self.controls;
        if th < c.minimum_heavy_temperature_k
            || th > c.maximum_heavy_temperature_k
            || te < c.minimum_electron_temperature_k
            || te > c.maximum_electron_temperature_k
        {
            return Err(FusionNetworkError::OutsideDomain);
        }
        self.plasma
            .admit_point(
                Temperature::new(th),
                Temperature::new(te),
                &self.amounts(point),
                point.volume,
                c.plasma_limits,
                &mut *cancelled,
            )
            .map_err(plasma_error)?;
        Ok([th, te])
    }
    pub(super) fn sources(
        &self,
        point: Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Sources, FusionNetworkError> {
        let temperatures = self.admit(point, cancelled)?;
        let mut events = [None; 4];
        for (r, reaction) in FusionReactionV1::ALL.iter().enumerate() {
            poll(cancelled)?;
            if self.channels.enabled[r] {
                events[r] = Some(
                    reaction
                        .thermal_event(Temperature::new(temperatures[0]))
                        .map_err(FusionNetworkError::Reaction)?,
                );
            }
        }
        let thermal = thermal::evaluate_components(
            self.transport,
            &self.thermal_amounts(point),
            &self.amounts(point),
            point.volume,
            temperatures,
            cancelled,
        )?;
        Ok(Sources { thermal, events })
    }

    /// One exact positive binary-reaction map at fixed coefficients. Remaining
    /// stocks are evaluated directly, preserving tiny survivors after deep burn.
    fn reaction_step(
        &self,
        old: Point,
        r: usize,
        dt: f64,
        event: FusionThermalEventV1,
    ) -> Result<(Point, f64), FusionNetworkError> {
        let reaction = FusionReactionV1::ALL[r];
        let indices = self.plasma.indices();
        let [a, b] = reaction.reactants().map(|s| indices[species_slot(s)]);
        let (na, nb) = (old.amounts[a].value, old.amounts[b].value);
        if na == 0.0 || nb == 0.0 || dt == 0.0 {
            return Ok((old, 0.0));
        }
        let rate = event.reactivity().value_m3_s();
        let volume = old.volume.value();
        let mut next = old;
        let extent = if a == b {
            let (remaining, extent) = dd_depletion(na, rate, volume, dt)?;
            next.amounts[a] = consume(old.amounts[a], 2.0 * extent, remaining);
            extent
        } else {
            let delta = (na - nb).abs();
            let (remaining, extent) = frozen_depletion(na.min(nb), delta, rate, volume, dt)?;
            let (small, large) = if na <= nb { (a, b) } else { (b, a) };
            next.amounts[small] = consume(old.amounts[small], extent, remaining);
            next.amounts[large] = consume(old.amounts[large], extent, remaining + delta);
            extent
        };
        next.extents[r].add(extent);
        let ion_removal = (event.reacting_ion_energy_j() * AVOGADRO) * extent;
        next.ion_removal.add(ion_removal);
        let mut thermal_increment = [-ion_removal, 0.0];
        for (p, product) in reaction.products().iter().enumerate() {
            let energy = (event.product_kinetic_energies_j()[p] * AVOGADRO) * extent;
            next.products[r][p].add(energy);
            if *product == NuclearSpeciesV1::Neutron {
                next.neutrons.add(energy);
            } else {
                next.amounts[indices[species_slot(*product)]].add(extent);
                if self.kinetic_alphas && *product == NuclearSpeciesV1::Helium4 {
                    next.kinetic_alpha_amount.add(extent);
                    continue;
                }
                let f = self.heating.product_electron_fractions[r][p];
                next.deposited[0].add((1.0 - f) * energy);
                next.deposited[1].add(f * energy);
                thermal_increment[0] += (1.0 - f) * energy;
                thermal_increment[1] += f * energy;
            }
        }
        self.add_pool_energy(&mut next, thermal_increment);
        Ok((next, extent))
    }

    fn radiate(&self, old: Point, dt: f64, coefficient: f64) -> Result<Point, FusionNetworkError> {
        if coefficient == 0.0 || dt == 0.0 {
            return Ok(old);
        }
        let loss = thermal::radiated_energy(
            self.temperatures(old)?[1],
            self.capacities(old)[1],
            dt,
            coefficient,
        )?;
        let mut next = old;
        next.photons.add(loss);
        self.add_pool_energy(&mut next, [0.0, -loss]);
        Ok(next)
    }
    fn thermal_step(
        &self,
        old: Point,
        dt: f64,
        sources: ThermalSources,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, FusionNetworkError> {
        poll(cancelled)?;
        if sources.conductance_w_k == 0.0 && sources.radiation_w_per_sqrt_k == 0.0 {
            return Ok(old);
        }
        let first = self.radiate(old, 0.5 * dt, sources.radiation_w_per_sqrt_k)?;
        self.admit(first, cancelled)?;
        let transfer = thermal::exchange_energy(
            self.temperatures(first)?,
            self.capacities(first),
            dt,
            sources.conductance_w_k,
        )?;
        let mut exchanged = first;
        exchanged.exchange.add(transfer);
        self.add_pool_energy(&mut exchanged, [transfer, -transfer]);
        self.admit(exchanged, cancelled)?;
        let next = self.radiate(exchanged, 0.5 * dt, sources.radiation_w_per_sqrt_k)?;
        self.admit(next, cancelled)?;
        Ok(next)
    }
    fn frozen_step(
        &self,
        old: Point,
        dt: f64,
        sources: Sources,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, FusionNetworkError> {
        self.frozen_step_with_births(old, dt, sources, cancelled)
            .map(|(point, _)| point)
    }

    /// Shared source step with independently accumulated DT/DHe3 alpha birth
    /// amounts. The kinetic owner injects those particles at the retained
    /// channel energy without differencing large cumulative event accounts.
    pub(super) fn frozen_step_with_births(
        &self,
        old: Point,
        dt: f64,
        sources: Sources,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Point, [f64; 2]), FusionNetworkError> {
        let mut births = [Sum::default(); 2];
        let mut next = self.thermal_step(old, 0.5 * dt, sources.thermal, cancelled)?;
        // Palindromic channel composition: each selected binary map is exact
        // at fixed T. Products can feed later maps without an overdraw or seed.
        for r in 0..4 {
            if let Some(event) = sources.events[r] {
                poll(cancelled)?;
                let (updated, extent) =
                    self.reaction_step(next, r, if r == 3 { dt } else { 0.5 * dt }, event)?;
                next = updated;
                if self.kinetic_alphas && (r == 0 || r == 3) {
                    births[usize::from(r == 3)].add(extent);
                }
                self.admit(next, cancelled)?;
            }
        }
        for r in (0..3).rev() {
            if let Some(event) = sources.events[r] {
                poll(cancelled)?;
                let (updated, extent) = self.reaction_step(next, r, 0.5 * dt, event)?;
                next = updated;
                if self.kinetic_alphas && r == 0 {
                    births[0].add(extent);
                }
                self.admit(next, cancelled)?;
            }
        }
        next.time = old.time + dt; // held external heat enters exactly once
        self.add_pool_energy(
            &mut next,
            [
                self.heating.heavy_power.value() * dt,
                self.heating.electron_power.value() * dt,
            ],
        );
        self.admit(next, cancelled)?;
        let next = self.thermal_step(next, 0.5 * dt, sources.thermal, cancelled)?;
        self.admit(next, cancelled)?;
        Ok((next, births.map(|n| n.value)))
    }

    pub(super) fn deposit_alpha_energy(&self, point: &mut Point, deposited: [f64; 2]) {
        self.add_pool_energy(point, deposited);
        for (account, energy) in point.alpha_deposited.iter_mut().zip(deposited) {
            account.add(energy);
        }
    }
    fn midpoint(
        &self,
        old: Point,
        end_time: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, FusionNetworkError> {
        let dt = end_time - old.time;
        let middle_time = old.time + 0.5 * dt;
        if middle_time <= old.time || middle_time >= end_time {
            return Err(FusionNetworkError::StepBudget);
        }
        let prepared = if let Some(history) = self.volume_history {
            self.volume_step(old, history.volume_at(Time::new(middle_time))?, cancelled)?
        } else {
            old
        };
        let first = self.sources(prepared, cancelled)?;
        let middle = self.frozen_step(prepared, 0.5 * dt, first, cancelled)?;
        let sources = self.sources(middle, cancelled)?;
        let mut burned = self.frozen_step(prepared, dt, sources, cancelled)?;
        burned.time = end_time;
        if let Some(history) = self.volume_history {
            self.volume_step(burned, history.volume_at(Time::new(end_time))?, cancelled)
        } else {
            Ok(burned)
        }
    }

    /// Exact fixed-inventory geometric map; its source clock never advances.
    fn volume_step(
        &self,
        old: Point,
        volume: Volume,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, FusionNetworkError> {
        poll(cancelled)?;
        if old.volume == volume {
            return Ok(old);
        }
        let energies = self.pool_energy(old);
        let work = adiabatic_work_on_pools(energies, old.volume, volume)?;
        let factor = exp((-2.0 / 3.0) * log_volume_ratio(old.volume, volume)?);
        let mut next = old;
        next.volume = volume;
        for i in 0..2 {
            next.work_on_pools[i].add(work[i]);
            if work[i] < -0.5 * energies[i] {
                // A direct positive survivor avoids subtracting almost all U.
                next.thermal_energy[i] = sum(energies[i] * factor);
            } else {
                // Compensate small work against repeated finite-volume motion.
                next.thermal_energy[i].add(work[i]);
            }
        }
        self.admit(next, cancelled)?;
        Ok(next)
    }

    pub(super) fn error(
        &self,
        old: Point,
        full: Point,
        fine: Point,
    ) -> Result<f64, FusionNetworkError> {
        let c = self.controls;
        let mut norm = 0.0_f64;
        let mut add = |a: f64, b: f64, absolute: f64| {
            norm = norm.max(
                (a - b).abs() / (3.0 * (absolute + c.relative_tolerance * a.abs().max(b.abs()))),
            );
        };
        for i in 0..6 {
            add(
                full.amounts[i].value,
                fine.amounts[i].value,
                c.absolute_amount_tolerance_mol,
            );
        }
        for r in 0..4 {
            add(
                full.extents[r].value - old.extents[r].value,
                fine.extents[r].value - old.extents[r].value,
                c.absolute_amount_tolerance_mol,
            );
            for p in 0..2 {
                add(
                    full.products[r][p].value - old.products[r][p].value,
                    fine.products[r][p].value - old.products[r][p].value,
                    c.energy_tolerance_j,
                );
            }
        }
        for (a, b) in self
            .temperatures(full)?
            .iter()
            .zip(self.temperatures(fine)?)
        {
            add(*a, b, c.absolute_temperature_tolerance_k);
        }
        for (a, b, before) in [
            (full.ion_removal, fine.ion_removal, old.ion_removal),
            (full.exchange, fine.exchange, old.exchange),
            (full.photons, fine.photons, old.photons),
            (
                full.alpha_deposited[0],
                fine.alpha_deposited[0],
                old.alpha_deposited[0],
            ),
            (
                full.alpha_deposited[1],
                fine.alpha_deposited[1],
                old.alpha_deposited[1],
            ),
            (
                full.work_on_pools[0],
                fine.work_on_pools[0],
                old.work_on_pools[0],
            ),
            (
                full.work_on_pools[1],
                fine.work_on_pools[1],
                old.work_on_pools[1],
            ),
        ] {
            add(
                a.value - before.value,
                b.value - before.value,
                c.energy_tolerance_j,
            );
        }
        if !norm.is_finite() {
            return Err(FusionNetworkError::Arithmetic(
                "network local error estimate",
            ));
        }
        Ok(norm)
    }

    pub(super) fn integrate(
        &self,
        duration: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Point, FusionNetworkBalance), FusionNetworkError> {
        let mut point = self.start();
        // Compensate accepted elapsed increments. Repeatedly rounding t+h
        // can otherwise leave an unresolvable sliver at an authored endpoint.
        let mut clock = Sum::default();
        let mut balance = FusionNetworkBalance::default();
        let mut step = self.controls.initial_step_s;
        for _ in 0..self.controls.maximum_steps {
            poll(cancelled)?;
            if point.time == duration {
                return Ok((point, balance));
            }
            let next_boundary = if let Some(history) = self.volume_history {
                history
                    .next_knot_after(Time::new(point.time))?
                    .map_or(duration, |time| time.value())
            } else {
                duration
            };
            let mut next_clock = clock;
            next_clock.add(step.min(self.controls.maximum_step_s));
            let end_time = next_clock.value.min(next_boundary);
            if end_time == next_boundary {
                next_clock = sum(next_boundary);
            }
            let dt = end_time - point.time;
            let middle_time = point.time + 0.5 * dt;
            if dt <= 0.0 || middle_time <= point.time || middle_time >= end_time {
                return Err(FusionNetworkError::StepBudget);
            }
            let trial = (|| {
                let full = self.midpoint(point, end_time, cancelled)?;
                let half = self.midpoint(point, middle_time, cancelled)?;
                let fine = self.midpoint(half, end_time, cancelled)?;
                self.admit(fine, cancelled)?;
                let error = self.error(point, full, fine)?;
                Ok::<_, FusionNetworkError>((fine, error))
            })();
            let error = match trial {
                Ok((next, error)) if error <= 1.0 => {
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
        if point.time == duration {
            Ok((point, balance))
        } else {
            Err(FusionNetworkError::StepBudget)
        }
    }

    pub(super) fn publish(
        &self,
        point: Point,
        balance: &mut FusionNetworkBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FrozenIonizedGasStateV1, FusionNetworkError> {
        let sources = self.sources(point, cancelled)?;
        let [th, te] = sources.thermal.temperatures;
        let amounts = self.amounts(point);
        let state = if point.time == 0.0 {
            self.initial.clone()
        } else {
            FrozenIonizedGasStateV1::new(
                self.plasma.thermodynamics(),
                &amounts,
                point.volume,
                Temperature::new(th),
                Temperature::new(te),
                64.0 * f64::EPSILON,
                &mut *cancelled,
            )
            .map_err(|e| match e {
                FrozenIonizedGasStateErrorV1::Cancelled => FusionNetworkError::Cancelled,
                e => FusionNetworkError::State(e),
            })?
        };
        self.balance_at(
            point,
            [
                state.heavy_internal_energy().value(),
                state.electron_internal_energy().value(),
            ],
            0.0,
            balance,
            cancelled,
        )?;
        Ok(state)
    }

    /// Complete nuclear/source accounts with an independently evaluated thermal
    /// state and, for a hybrid owner, its additional kinetic storage change.
    pub(super) fn balance_at(
        &self,
        point: Point,
        thermal_energy: [f64; 2],
        kinetic_energy_change_j: f64,
        balance: &mut FusionNetworkBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FusionNetworkError> {
        let amounts = self.amounts(point);
        let extents = point.extents.map(|x| x.value);
        let mut released = Sum::default();
        for (q, extent) in self.q_molar.iter().zip(extents) {
            released.add(q * extent);
        }
        let neutron_amount = extents[0] + extents[1];
        let expected = self.accounted_energy(point);
        let heavy_residual = thermal_energy[0] - expected[0];
        let electron_residual = thermal_energy[1] - expected[1];
        let external_heavy = self.heating.heavy_power.value() * point.time;
        let external_electron = self.heating.electron_power.value() * point.time;
        let change = (thermal_energy[0] - self.initial_energy[0])
            + (thermal_energy[1] - self.initial_energy[1])
            + kinetic_energy_change_j;
        let total_residual = change + point.neutrons.value + point.photons.value
            - released.value
            - external_heavy
            - external_electron
            - point.work_on_pools[0].value
            - point.work_on_pools[1].value;
        let nuclei = [
            NuclearSpeciesV1::Proton,
            NuclearSpeciesV1::Deuterium,
            NuclearSpeciesV1::Tritium,
            NuclearSpeciesV1::Helium3,
            NuclearSpeciesV1::Helium4,
        ];
        let mut baryons = [0.0, neutron_amount];
        let mut charge = [0.0; 2];
        let mut stoichiometry = 0.0_f64;
        for (species, index) in nuclei.iter().zip(self.plasma.indices()) {
            let before = self.initial_amounts[index];
            let after = amounts[index].value();
            baryons[0] += f64::from(species.baryon_number()) * before;
            baryons[1] += f64::from(species.baryon_number()) * after;
            charge[0] += f64::from(species.charge_number()) * before;
            charge[1] += f64::from(species.charge_number()) * after;
            let mut expected_amount = sum(before);
            let mut activity = before;
            for (reaction, extent) in FusionReactionV1::ALL.iter().zip(extents) {
                let transfer = f64::from(reaction.stoichiometric_coefficient(*species)) * extent;
                expected_amount.add(transfer);
                activity += transfer.abs();
            }
            let residual = (after - expected_amount.value).abs();
            stoichiometry = stoichiometry.max(residual);
            if !residual.is_finite()
                || residual
                    > self.controls.absolute_amount_tolerance_mol
                        + self.controls.balance_relative_tolerance * activity
            {
                return Err(FusionNetworkError::Balance(
                    "network species and channel extents",
                ));
            }
        }
        let baryon_residual = (baryons[1] - baryons[0]) / baryons[0];
        let charge_residual = (charge[1] - charge[0]) / charge[0];
        let tolerance = self.controls.energy_tolerance_j
            + self.controls.balance_relative_tolerance
                * (released.value.abs()
                    + external_heavy.abs()
                    + external_electron.abs()
                    + point.work_on_pools[0].value.abs()
                    + point.work_on_pools[1].value.abs()
                    + kinetic_energy_change_j.abs()
                    + self.initial_energy[0].abs()
                    + self.initial_energy[1].abs());
        if [
            heavy_residual,
            electron_residual,
            total_residual,
            baryon_residual,
            charge_residual,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || [heavy_residual, electron_residual, total_residual]
                .iter()
                .any(|x| x.abs() > tolerance)
            || baryon_residual.abs() > self.controls.balance_relative_tolerance
            || charge_residual.abs() > self.controls.balance_relative_tolerance
        {
            return Err(FusionNetworkError::Balance(
                "network source energy, baryons or charge",
            ));
        }
        *balance = FusionNetworkBalance {
            reaction_extents_mol: extents,
            product_kinetic_energies_j: point.products.map(|row| row.map(|x| x.value)),
            nuclear_energy_released_j: released.value,
            charged_product_electron_energy_j: point.deposited[1].value,
            charged_product_heavy_energy_j: point.deposited[0].value,
            reacting_ion_energy_removed_j: point.ion_removal.value,
            neutron_amount_mol: neutron_amount,
            neutron_energy_escaped_j: point.neutrons.value,
            electron_to_heavy_coulomb_energy_j: point.exchange.value,
            photon_energy_escaped_j: point.photons.value,
            external_electron_energy_j: external_electron,
            external_heavy_energy_j: external_heavy,
            boundary_work_on_heavy_j: point.work_on_pools[0].value,
            boundary_work_on_electrons_j: point.work_on_pools[1].value,
            electron_energy_residual_j: electron_residual,
            heavy_energy_residual_j: heavy_residual,
            total_energy_residual_j: total_residual,
            relative_baryon_residual: baryon_residual,
            relative_charge_residual: charge_residual,
            maximum_stoichiometry_residual_mol: stoichiometry,
            ..*balance
        };
        poll(cancelled)?;
        Ok(())
    }
}

/// For one DD branch, event rate is half nD²<sigma v>, and two D disappear
/// per event. Thus D'=-NA<sigma v>D²/V, extent=(D_before-D_after)/2.
pub(super) fn dd_depletion(
    amount: f64,
    rate: f64,
    volume: f64,
    dt: f64,
) -> Result<(f64, f64), FusionNetworkError> {
    if amount == 0.0 || dt == 0.0 {
        return Ok((amount, 0.0));
    }
    let exposure = ((amount / volume) * AVOGADRO) * rate * dt;
    if exposure.is_nan() || exposure < 0.0 {
        return Err(FusionNetworkError::Arithmetic("DD finite-stock exposure"));
    }
    if exposure.is_infinite() {
        return Ok((0.0, 0.5 * amount));
    }
    let remaining = amount / (1.0 + exposure);
    let extent = 0.5 * amount * (exposure / (1.0 + exposure));
    if !remaining.is_finite()
        || !extent.is_finite()
        || remaining < 0.0
        || extent < 0.0
        || remaining > amount
        || extent > 0.5 * amount
    {
        return Err(FusionNetworkError::Arithmetic("DD finite-stock map"));
    }
    Ok((remaining, extent))
}
fn sum(value: f64) -> Sum {
    let mut result = Sum::default();
    result.add(value);
    result
}

/// Compensate small withdrawals against the same event extents credited to
/// products. Replacing a large stock by a freshly rounded analytic remainder
/// at every tiny substep would accumulate artificial inventory/mass drift.
/// For a deep withdrawal use the direct positive survivor formula instead:
/// subtracting nearly the entire stock would erase a representable trace.
/// The branch keeps a factor-two positivity margin; no stock is clipped.
fn consume(old: Sum, removed: f64, analytic_remaining: f64) -> Sum {
    if removed <= 0.5 * old.value {
        let mut next = old;
        next.add(-removed);
        next
    } else {
        sum(analytic_remaining)
    }
}
fn plasma_error(error: DtFusionPlasmaErrorV1) -> FusionNetworkError {
    match error {
        DtFusionPlasmaErrorV1::Cancelled => FusionNetworkError::Cancelled,
        error => FusionNetworkError::Plasma(error),
    }
}
