use super::super::solve::{Sum, frozen_depletion};
use super::super::thermal::{self, ThermalSources};
use super::*;
use fs_math::det::pow;
use fs_qty::Temperature;
use fs_thermochem::coulomb::CoulombErrorV1;
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateErrorV1;
use fs_thermochem::fusion::FusionThermalEventV1;
use fs_thermochem::fusion::plasma::DtFusionPlasmaErrorV1;
use fs_thermochem::plasma_radiation::BremsstrahlungErrorV1;

#[derive(Clone, Copy)]
pub(super) struct Point {
    time: f64,
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
struct Sources {
    thermal: ThermalSources,
    events: [Option<FusionThermalEventV1>; 4],
}

pub(super) struct Problem<'a> {
    initial: &'a FrozenIonizedGasStateV1,
    plasma: &'a LightIonFusionPlasmaV1,
    channels: FusionNetworkChannels,
    heating: FusionNetworkHeating,
    transport: &'a FusionNetworkTransport,
    controls: FusionNetworkSolve,
    initial_amounts: [f64; 6],
    initial_energy: [f64; 2],
    q_molar: [f64; 4],
}

impl<'a> Problem<'a> {
    pub(super) fn new(
        initial: &'a FrozenIonizedGasStateV1,
        plasma: &'a LightIonFusionPlasmaV1,
        channels: FusionNetworkChannels,
        heating: FusionNetworkHeating,
        transport: &'a FusionNetworkTransport,
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
            controls,
            initial_amounts: core::array::from_fn(|i| initial.amounts()[i].value()),
            initial_energy: [
                initial.heavy_internal_energy().value(),
                initial.electron_internal_energy().value(),
            ],
            q_molar: FusionReactionV1::ALL.map(|r| r.q_j() * AVOGADRO),
        };
        result.sources(result.start(), cancelled)?;
        Ok(result)
    }

    fn start(&self) -> Point {
        Point {
            time: 0.0,
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
    fn amounts(&self, point: Point) -> [Amount; 6] {
        point.amounts.map(|n| Amount::new(n.value))
    }
    fn capacities(&self, point: Point) -> [f64; 2] {
        [
            CV_MOLAR
                * self
                    .plasma
                    .indices()
                    .iter()
                    .map(|&i| point.amounts[i].value)
                    .sum::<f64>(),
            CV_MOLAR * point.amounts[self.plasma.electron_index()].value,
        ]
    }
    fn pool_energy(&self, point: Point) -> [f64; 2] {
        [
            self.initial_energy[0]
                + self.heating.heavy_power.value() * point.time
                + (point.deposited[0].value - point.ion_removal.value)
                + point.exchange.value,
            self.initial_energy[1]
                + self.heating.electron_power.value() * point.time
                + point.deposited[1].value
                - point.exchange.value
                - point.photons.value,
        ]
    }
    fn temperatures(&self, point: Point) -> Result<[f64; 2], FusionNetworkError> {
        // Keep exact authored endpoints when no physical transfer has occurred.
        // A thermal substep can act at t=0, so elapsed time alone is insufficient.
        if point.extents.iter().all(|x| x.value == 0.0)
            && point.exchange.value == 0.0
            && point.photons.value == 0.0
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
    fn admit(
        &self,
        point: Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], FusionNetworkError> {
        poll(cancelled)?;
        if !point.time.is_finite()
            || point.time < 0.0
            || !point.exchange.value.is_finite()
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
                self.initial.volume(),
                c.plasma_limits,
                &mut *cancelled,
            )
            .map_err(plasma_error)?;
        Ok([th, te])
    }
    fn sources(
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
        let thermal = thermal::evaluate(
            self.transport,
            &self.amounts(point),
            self.initial.volume(),
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
    ) -> Result<Point, FusionNetworkError> {
        let reaction = FusionReactionV1::ALL[r];
        let indices = self.plasma.indices();
        let [a, b] = reaction.reactants().map(|s| indices[species_slot(s)]);
        let (na, nb) = (old.amounts[a].value, old.amounts[b].value);
        if na == 0.0 || nb == 0.0 || dt == 0.0 {
            return Ok(old);
        }
        let rate = event.reactivity().value_m3_s();
        let volume = self.initial.volume().value();
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
        next.ion_removal
            .add((event.reacting_ion_energy_j() * AVOGADRO) * extent);
        for (p, product) in reaction.products().iter().enumerate() {
            let energy = (event.product_kinetic_energies_j()[p] * AVOGADRO) * extent;
            next.products[r][p].add(energy);
            if *product == NuclearSpeciesV1::Neutron {
                next.neutrons.add(energy);
            } else {
                next.amounts[indices[species_slot(*product)]].add(extent);
                let f = self.heating.product_electron_fractions[r][p];
                next.deposited[0].add((1.0 - f) * energy);
                next.deposited[1].add(f * energy);
            }
        }
        Ok(next)
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
        let mut next = self.thermal_step(old, 0.5 * dt, sources.thermal, cancelled)?;
        // Palindromic channel composition: each selected binary map is exact
        // at fixed T. Products can feed later maps without an overdraw or seed.
        for r in 0..4 {
            if let Some(event) = sources.events[r] {
                poll(cancelled)?;
                next = self.reaction_step(next, r, if r == 3 { dt } else { 0.5 * dt }, event)?;
                self.admit(next, cancelled)?;
            }
        }
        for r in (0..3).rev() {
            if let Some(event) = sources.events[r] {
                poll(cancelled)?;
                next = self.reaction_step(next, r, 0.5 * dt, event)?;
                self.admit(next, cancelled)?;
            }
        }
        next.time = old.time + dt; // held external work enters exactly once
        self.admit(next, cancelled)?;
        let next = self.thermal_step(next, 0.5 * dt, sources.thermal, cancelled)?;
        self.admit(next, cancelled)?;
        Ok(next)
    }
    fn midpoint(
        &self,
        old: Point,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, FusionNetworkError> {
        let first = self.sources(old, cancelled)?;
        let middle = self.frozen_step(old, 0.5 * dt, first, cancelled)?;
        let sources = self.sources(middle, cancelled)?;
        self.frozen_step(old, dt, sources, cancelled)
    }

    fn error(&self, old: Point, full: Point, fine: Point) -> Result<f64, FusionNetworkError> {
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
        let mut balance = FusionNetworkBalance::default();
        let mut step = self.controls.initial_step_s;
        for _ in 0..self.controls.maximum_steps {
            poll(cancelled)?;
            if point.time == duration {
                return Ok((point, balance));
            }
            let dt = step
                .min(self.controls.maximum_step_s)
                .min(duration - point.time);
            if dt <= 0.0 || point.time + dt == point.time {
                return Err(FusionNetworkError::StepBudget);
            }
            let trial = (|| {
                let full = self.midpoint(point, dt, cancelled)?;
                let half = self.midpoint(point, 0.5 * dt, cancelled)?;
                let mut fine = self.midpoint(half, 0.5 * dt, cancelled)?;
                fine.time = point.time + dt;
                self.admit(fine, cancelled)?;
                let error = self.error(point, full, fine)?;
                Ok::<_, FusionNetworkError>((fine, error))
            })();
            let error = match trial {
                Ok((next, error)) if error <= 1.0 => {
                    point = next;
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
                self.initial.volume(),
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
        let extents = point.extents.map(|x| x.value);
        let mut released = Sum::default();
        for (q, extent) in self.q_molar.iter().zip(extents) {
            released.add(q * extent);
        }
        let neutron_amount = extents[0] + extents[1];
        let expected = self.pool_energy(point);
        let heavy_residual = state.heavy_internal_energy().value() - expected[0];
        let electron_residual = state.electron_internal_energy().value() - expected[1];
        let external_heavy = self.heating.heavy_power.value() * point.time;
        let external_electron = self.heating.electron_power.value() * point.time;
        let change = (state.heavy_internal_energy().value() - self.initial_energy[0])
            + (state.electron_internal_energy().value() - self.initial_energy[1]);
        let total_residual = change + point.neutrons.value + point.photons.value
            - released.value
            - external_heavy
            - external_electron;
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
            electron_energy_residual_j: electron_residual,
            heavy_energy_residual_j: heavy_residual,
            total_energy_residual_j: total_residual,
            relative_baryon_residual: baryon_residual,
            relative_charge_residual: charge_residual,
            maximum_stoichiometry_residual_mol: stoichiometry,
            ..*balance
        };
        poll(cancelled)?;
        Ok(state)
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
