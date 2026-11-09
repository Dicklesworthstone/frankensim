use super::emission::NeutronMoments;
use super::*;
use fs_math::det::{exp, expm1, pow};

mod thermal;

/// Compensated cumulative physical transfers, retained across adaptive steps.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Sum {
    pub(super) value: f64,
    correction: f64,
}
impl Sum {
    pub(super) const fn new(value: f64) -> Self {
        Self {
            value,
            correction: 0.0,
        }
    }

    pub(super) fn add(&mut self, x: f64) {
        let y = x - self.correction;
        let next = self.value + y;
        self.correction = (next - self.value) - y;
        self.value = next;
    }

    /// Retain a small local transfer even after a large accumulated history.
    /// Kahan's low part has the opposite sign of `correction`.
    pub(super) fn increment_since(self, previous: Self) -> f64 {
        (self.value - previous.value) + (previous.correction - self.correction)
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct Point {
    time: f64,
    // Confinement can remove almost all initial energy. Retain its positive
    // survivor directly instead of reconstructing initial minus total loss.
    thermal_energy: [Sum; 2],
    transport_energy_out: [Sum; 2],
    limiting_fuel: f64,
    extent: Sum,
    alpha_energy: Sum,
    ion_removal: Sum,
    neutron_energy: Sum,
    neutron_interval: Option<NeutronMoments>,
    coulomb_transfer: Sum,
    photon_energy: Sum,
    // Only the open-particle path uses these directly evolved canonical stocks.
    particle_amounts: [Sum; 4],
    particles: particles::Accounts,
}

pub(super) struct Problem<'a> {
    initial: &'a FrozenIonizedGasStateV1,
    plasma: &'a DtFusionPlasmaV1,
    heating: DtFusionHeating,
    transport: &'a DtFusionTransport,
    controls: DtFusionSolve,
    initial_amounts: [f64; 4],
    deuterium_is_limiting: bool,
    excess_fuel: f64,
    initial_pool_energy: [f64; 2],
    q_molar: f64,
    particles: Option<particles::Prepared>,
}

impl<'a> Problem<'a> {
    pub fn new(
        initial: &'a FrozenIonizedGasStateV1,
        plasma: &'a DtFusionPlasmaV1,
        heating: DtFusionHeating,
        transport: &'a DtFusionTransport,
        controls: DtFusionSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DtFusionError> {
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
        for t in [
            controls.minimum_heavy_temperature_k,
            controls.maximum_heavy_temperature_k,
        ] {
            FusionReactionV1::DeuteriumTritium
                .reactivity(Temperature::new(t))
                .map_err(DtFusionError::Reaction)?;
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
            return Err(DtFusionError::Invalid(
                "thermal operators must retain the original fusion source bank",
            ));
        }
        let amounts = core::array::from_fn(|i| initial.amounts()[i].value());
        let [d, t, he] = plasma.indices();
        let particles = transport
            .particles
            .map(|law| {
                particles::Prepared::new(
                    law,
                    plasma.thermodynamics(),
                    [None, Some(d), Some(t), None, Some(he)],
                    plasma.electron_index(),
                    controls.plasma_limits,
                    cancelled,
                )
            })
            .transpose()?;
        let result = Self {
            initial,
            plasma,
            heating,
            transport,
            controls,
            initial_amounts: amounts,
            deuterium_is_limiting: amounts[d] <= amounts[t],
            excess_fuel: (amounts[d] - amounts[t]).abs(),
            initial_pool_energy: [
                initial.heavy_internal_energy().value(),
                initial.electron_internal_energy().value(),
            ],
            q_molar: FusionReactionV1::DeuteriumTritium.q_j() * AVOGADRO,
            particles,
        };
        result.sources(result.start(), cancelled)?;
        Ok(result)
    }

    fn start(&self) -> Point {
        let [d, t, _] = self.plasma.indices();
        Point {
            time: 0.0,
            thermal_energy: self.initial_pool_energy.map(Sum::new),
            transport_energy_out: [Sum::default(); 2],
            limiting_fuel: self.initial_amounts[d].min(self.initial_amounts[t]),
            extent: Sum::default(),
            alpha_energy: Sum::default(),
            ion_removal: Sum::default(),
            neutron_energy: Sum::default(),
            neutron_interval: self
                .controls
                .record_neutron_history
                .then(NeutronMoments::default),
            coulomb_transfer: Sum::default(),
            photon_energy: Sum::default(),
            particle_amounts: self.initial_amounts.map(Sum::new),
            particles: particles::Accounts::default(),
        }
    }

    fn amounts(&self, point: Point) -> [Amount; 4] {
        if self.particles.is_some() {
            return point.particle_amounts.map(|n| Amount::new(n.value));
        }
        let mut amounts = self.initial_amounts.map(Amount::new);
        let [d, t, he] = self.plasma.indices();
        let (limiting, other) = if self.deuterium_is_limiting {
            (d, t)
        } else {
            (t, d)
        };
        amounts[limiting] = Amount::new(point.limiting_fuel);
        amounts[other] = Amount::new(point.limiting_fuel + self.excess_fuel);
        amounts[he] = Amount::new(self.initial_amounts[he] + point.extent.value);
        amounts
    }

    fn pool_energy(&self, point: Point) -> [f64; 2] {
        if self.retains_thermal_energy() {
            return point.thermal_energy.map(|energy| energy.value);
        }
        self.accounted_energy(point)
    }

    fn accounted_energy(&self, point: Point) -> [f64; 2] {
        let f = self.heating.alpha_electron_fraction;
        [
            self.initial_pool_energy[0]
                + self.heating.heavy_power.value() * point.time
                + ((1.0 - f) * point.alpha_energy.value - point.ion_removal.value)
                + point.coulomb_transfer.value
                - point.transport_energy_out[0].value
                + point.particles.energy_in[0].value
                - point.particles.energy_out[0].value,
            self.initial_pool_energy[1]
                + self.heating.electron_power.value() * point.time
                + f * point.alpha_energy.value
                - point.coulomb_transfer.value
                - point.photon_energy.value
                - point.transport_energy_out[1].value
                + point.particles.energy_in[1].value
                - point.particles.energy_out[1].value,
        ]
    }

    fn add_pool_energy(&self, point: &mut Point, increments: [f64; 2]) {
        if self.retains_thermal_energy() {
            for (energy, increment) in point.thermal_energy.iter_mut().zip(increments) {
                energy.add(increment);
            }
        }
    }

    fn retains_thermal_energy(&self) -> bool {
        self.transport.energy_confinement.is_some() || self.particles.is_some()
    }

    fn exchange_particles(
        &self,
        old: Point,
        dt: f64,
        other_sources: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, DtFusionError> {
        let Some(source) = &self.particles else {
            return Ok(old);
        };
        let mut next = old;
        source.apply(
            &mut next.particle_amounts,
            None,
            0.0,
            &mut next.thermal_energy,
            &mut next.particles,
            dt,
            other_sources,
            cancelled,
        )?;
        self.admit(next, cancelled)?;
        Ok(next)
    }

    fn temperatures(&self, point: Point) -> Result<[f64; 2], DtFusionError> {
        // Preserve the authored starting temperatures exactly. Reconstructing
        // NASA-9 (h - RT) through a separately rounded Cv*T can move an exact
        // source-domain endpoint one ULP outside before any time has elapsed.
        if point.time == 0.0
            && point.particles.is_zero()
            && point.coulomb_transfer.value == 0.0
            && point.photon_energy.value == 0.0
            && point
                .transport_energy_out
                .iter()
                .all(|energy| energy.value == 0.0)
        {
            return Ok([
                self.initial.heavy_temperature().value(),
                self.initial.electron_temperature().value(),
            ]);
        }
        let amounts = self.amounts(point);
        let heavy = self
            .plasma
            .indices()
            .iter()
            .map(|&i| amounts[i].value())
            .sum::<f64>();
        let electrons = amounts[self.plasma.electron_index()].value();
        let energy = self.pool_energy(point);
        let temperatures = [
            energy[0] / (CV_MOLAR * heavy),
            energy[1] / (CV_MOLAR * electrons),
        ];
        if temperatures.iter().any(|&t| !positive(t)) {
            return Err(DtFusionError::OutsideDomain);
        }
        Ok(temperatures)
    }

    fn admit(
        &self,
        point: Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], DtFusionError> {
        poll(cancelled)?;
        if [
            point.time,
            point.limiting_fuel,
            point.extent.value,
            point.alpha_energy.value,
            point.ion_removal.value,
            point.neutron_energy.value,
            point.photon_energy.value,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
            || !point.coulomb_transfer.value.is_finite()
            || point
                .transport_energy_out
                .iter()
                .any(|energy| !energy.value.is_finite())
            || !point.particles.finite()
            || point
                .particle_amounts
                .iter()
                .any(|n| !n.value.is_finite() || n.value < 0.0)
            || (self.retains_thermal_energy()
                && point
                    .thermal_energy
                    .iter()
                    .any(|energy| !energy.value.is_finite()))
        {
            return Err(DtFusionError::Arithmetic("finite nonnegative burn state"));
        }
        let [th, te] = self.temperatures(point)?;
        let c = self.controls;
        if th < c.minimum_heavy_temperature_k
            || th > c.maximum_heavy_temperature_k
            || te < c.minimum_electron_temperature_k
            || te > c.maximum_electron_temperature_k
        {
            return Err(DtFusionError::OutsideDomain);
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

    /// Exact finite-fuel map with rate and per-event energies frozen at T.
    fn frozen_step(&self, old: Point, dt: f64, temperature: f64) -> Result<Point, DtFusionError> {
        let reaction = FusionReactionV1::DeuteriumTritium;
        let event = reaction
            .thermal_event(Temperature::new(temperature))
            .map_err(DtFusionError::Reaction)?;
        let [d, t, he] = self.plasma.indices();
        let (limiting, excess) = if self.particles.is_some() {
            let d = old.particle_amounts[d].value;
            let t = old.particle_amounts[t].value;
            (d.min(t), (d - t).abs())
        } else {
            (old.limiting_fuel, self.excess_fuel)
        };
        let (remaining, extent) = frozen_depletion(
            limiting,
            excess,
            event.reactivity().value_m3_s(),
            self.initial.volume().value(),
            dt,
        )?;
        let removed = event.reacting_ion_energy_j();
        let [alpha, neutron] = event.product_kinetic_energies_j();
        let alpha_energy = (alpha * AVOGADRO) * extent;
        let removed_energy = (removed * AVOGADRO) * extent;
        let mut next = old;
        next.time = old.time + dt;
        next.limiting_fuel = remaining;
        if self.particles.is_some() {
            let (small, large) = if old.particle_amounts[d].value <= old.particle_amounts[t].value {
                (d, t)
            } else {
                (t, d)
            };
            next.particle_amounts[small] =
                particles::withdraw(old.particle_amounts[small], extent, remaining);
            next.particle_amounts[large] =
                particles::withdraw(old.particle_amounts[large], extent, remaining + excess);
            next.particle_amounts[he].add(extent);
        }
        next.extent.add(extent);
        next.alpha_energy.add(alpha_energy);
        next.ion_removal.add(removed_energy);
        let neutron_energy = (neutron * AVOGADRO) * extent;
        next.neutron_energy.add(neutron_energy);
        if let Some(moments) = &mut next.neutron_interval {
            moments.add(0, extent, neutron_energy);
        }
        let fraction = self.heating.alpha_electron_fraction;
        self.add_pool_energy(
            &mut next,
            [
                (1.0 - fraction) * alpha_energy - removed_energy
                    + self.heating.heavy_power.value() * dt,
                fraction * alpha_energy + self.heating.electron_power.value() * dt,
            ],
        );
        Ok(next)
    }

    fn midpoint(
        &self,
        old: Point,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, DtFusionError> {
        let initial_sources = self.sources(old, cancelled)?;
        let middle = self.frozen_coupled_step(old, 0.5 * dt, initial_sources, cancelled)?;
        let middle_sources = self.sources(middle, cancelled)?;
        self.frozen_coupled_step(old, dt, middle_sources, cancelled)
    }

    fn error(&self, old: Point, full: Point, fine: Point) -> Result<f64, DtFusionError> {
        let c = self.controls;
        let mut norm = 0.0_f64;
        let mut add = |a: f64, b: f64, absolute: f64| {
            // For smooth admitted evolution the accepted two-half midpoint
            // solution has an O(dt³) local defect; the divisor is 2²-1.
            // At a cutoff crossover this remains an indicator, not a bound.
            norm = norm.max(
                (a - b).abs() / (3.0 * (absolute + c.relative_tolerance * a.abs().max(b.abs()))),
            );
        };
        for (a, b) in self.amounts(full).iter().zip(self.amounts(fine)) {
            add(a.value(), b.value(), c.absolute_amount_tolerance_mol);
        }
        // Resolve newly generated trace ash even with a large preexisting pool.
        add(
            full.extent.value - old.extent.value,
            fine.extent.value - old.extent.value,
            c.absolute_amount_tolerance_mol,
        );
        for (a, b) in self
            .temperatures(full)?
            .iter()
            .zip(self.temperatures(fine)?)
        {
            add(*a, b, c.absolute_temperature_tolerance_k);
        }
        add(
            full.neutron_energy.value - old.neutron_energy.value,
            fine.neutron_energy.value - old.neutron_energy.value,
            c.energy_tolerance_j,
        );
        add(
            full.alpha_energy.value - old.alpha_energy.value,
            fine.alpha_energy.value - old.alpha_energy.value,
            c.energy_tolerance_j,
        );
        add(
            full.coulomb_transfer.value - old.coulomb_transfer.value,
            fine.coulomb_transfer.value - old.coulomb_transfer.value,
            c.energy_tolerance_j,
        );
        add(
            full.photon_energy.value - old.photon_energy.value,
            fine.photon_energy.value - old.photon_energy.value,
            c.energy_tolerance_j,
        );
        for i in 0..2 {
            add(
                full.transport_energy_out[i].value - old.transport_energy_out[i].value,
                fine.transport_energy_out[i].value - old.transport_energy_out[i].value,
                c.energy_tolerance_j,
            );
        }
        old.particles.compare(
            full.particles,
            fine.particles,
            c.absolute_amount_tolerance_mol,
            c.energy_tolerance_j,
            &mut add,
        );
        if !norm.is_finite() {
            return Err(DtFusionError::Arithmetic("local error estimate"));
        }
        Ok(norm)
    }

    pub fn integrate(
        &self,
        duration: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Point, DtFusionBalance, emission::NeutronHistory), DtFusionError> {
        let mut point = self.start();
        let mut balance = DtFusionBalance::default();
        let mut neutron_history = self.controls.record_neutron_history.then(Vec::new);
        let mut step = self.controls.initial_step_s;
        for _ in 0..self.controls.maximum_steps {
            poll(cancelled)?;
            if point.time == duration {
                return Ok((point, balance, neutron_history));
            }
            let remaining = duration - point.time;
            let dt = step.min(self.controls.maximum_step_s).min(remaining);
            if dt <= 0.0 || point.time + dt == point.time {
                return Err(DtFusionError::StepBudget);
            }
            let trial = (|| {
                let mut start = point;
                if let Some(moments) = &mut start.neutron_interval {
                    *moments = NeutronMoments::default();
                }
                let full = self.midpoint(start, dt, cancelled)?;
                let half = self.midpoint(start, 0.5 * dt, cancelled)?;
                let mut fine = self.midpoint(half, 0.5 * dt, cancelled)?;
                // Hit the authored interval endpoint despite half-step rounding.
                fine.time = point.time + dt;
                self.admit(fine, cancelled)?;
                let error = self.error(point, full, fine)?;
                Ok::<_, DtFusionError>((fine, error))
            })();
            let error = match trial {
                Ok((next, error)) if error <= 1.0 => {
                    emission::record_interval(
                        &mut neutron_history,
                        point.time,
                        next.time,
                        next.neutron_interval,
                    )?;
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
                Err(DtFusionError::StepBudget) => {
                    // This is a local frozen-source resolution refusal. The
                    // outer attempt/minimum-step budgets remain authoritative.
                    balance.rejected_steps += 1;
                    step = 0.5 * dt;
                    if step < self.controls.minimum_step_s {
                        return Err(DtFusionError::StepBudget);
                    }
                    continue;
                }
                Err(
                    DtFusionError::OutsideDomain
                    | DtFusionError::Plasma(DtFusionPlasmaErrorV1::Domain(_))
                    | DtFusionError::Coulomb(CoulombErrorV1::Domain(_))
                    | DtFusionError::Radiation(BremsstrahlungErrorV1::Domain(_)),
                ) => None,
                Err(error) => return Err(error),
            };
            balance.rejected_steps += 1;
            let factor = error.map_or(0.5, |e| (0.8 / pow(e, 1.0 / 3.0)).clamp(0.1, 0.5));
            step = dt * factor;
            if step < self.controls.minimum_step_s {
                return Err(if error.is_some() {
                    DtFusionError::StepBudget
                } else {
                    DtFusionError::OutsideDomain
                });
            }
        }
        if point.time == duration {
            Ok((point, balance, neutron_history))
        } else {
            Err(DtFusionError::StepBudget)
        }
    }

    pub fn publish(
        &self,
        point: Point,
        balance: &mut DtFusionBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FrozenIonizedGasStateV1, DtFusionError> {
        let sources = self.sources(point, cancelled)?;
        let [th, te] = sources.temperatures;
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
            .map_err(|error| match error {
                FrozenIonizedGasStateErrorV1::Cancelled => DtFusionError::Cancelled,
                error => DtFusionError::State(error),
            })?
        };
        let released = self.q_molar * point.extent.value;
        let external_heavy = self.heating.heavy_power.value() * point.time;
        let external_electron = self.heating.electron_power.value() * point.time;
        let expected = self.accounted_energy(point);
        let heavy_residual = state.heavy_internal_energy().value() - expected[0];
        let electron_residual = state.electron_internal_energy().value() - expected[1];
        // Source differences are computed separately by pool, avoiding subtraction
        // of the huge rest-energy datum (which is deliberately not in the gas U).
        let actual_change = (state.heavy_internal_energy().value() - self.initial_pool_energy[0])
            + (state.electron_internal_energy().value() - self.initial_pool_energy[1]);
        let total_residual = actual_change
            + point.neutron_energy.value
            + point.photon_energy.value
            + point.transport_energy_out[0].value
            + point.transport_energy_out[1].value
            + point.particles.energy_out[0].value
            + point.particles.energy_out[1].value
            - point.particles.energy_in[0].value
            - point.particles.energy_in[1].value
            - released
            - external_heavy
            - external_electron;
        let [d, t, he] = self.plasma.indices();
        let before_baryons = 2.0 * self.initial_amounts[d]
            + 3.0 * self.initial_amounts[t]
            + 4.0 * self.initial_amounts[he];
        let after_baryons = 2.0 * amounts[d].value()
            + 3.0 * amounts[t].value()
            + 4.0 * amounts[he].value()
            + point.extent.value;
        let before_charge =
            self.initial_amounts[d] + self.initial_amounts[t] + 2.0 * self.initial_amounts[he];
        let after_charge = amounts[d].value() + amounts[t].value() + 2.0 * amounts[he].value();
        let particles = point.particles;
        let net = particles.nuclei_net_in();
        let activity = particles.nuclei_activity();
        let baryon_residual =
            (after_baryons - before_baryons - net[0]) / (before_baryons + activity[0]);
        let charge_residual =
            (after_charge - before_charge - net[1]) / (before_charge + activity[1]);
        for (slot, index, coefficient) in [(1, d, -1.0), (2, t, -1.0), (4, he, 1.0)] {
            let mut expected_amount = Sum::new(self.initial_amounts[index]);
            expected_amount.add(coefficient * point.extent.value);
            expected_amount.add(particles.injected[slot].value);
            expected_amount.add(-particles.exhausted[slot].value);
            let activity = self.initial_amounts[index]
                + point.extent.value
                + particles.injected[slot].value
                + particles.exhausted[slot].value;
            if (amounts[index].value() - expected_amount.value).abs()
                > self.controls.absolute_amount_tolerance_mol
                    + self.controls.balance_relative_tolerance * activity
            {
                return Err(DtFusionError::Balance(
                    "D-T reaction and particle-source inventory",
                ));
            }
        }
        let electron = self.plasma.electron_index();
        let electron_change = amounts[electron].value() - self.initial_amounts[electron];
        let electron_particle_residual =
            electron_change - particles.electrons_in.value + particles.electrons_out.value;
        if !electron_particle_residual.is_finite()
            || electron_particle_residual.abs()
                > self.controls.absolute_amount_tolerance_mol
                    + self.controls.balance_relative_tolerance
                        * (self.initial_amounts[electron]
                            + particles.electrons_in.value
                            + particles.electrons_out.value)
        {
            return Err(DtFusionError::Balance(
                "particle-source companion electrons",
            ));
        }
        let tolerance = self.controls.energy_tolerance_j
            + self.controls.balance_relative_tolerance
                * (released.abs()
                    + external_heavy.abs()
                    + external_electron.abs()
                    + point.transport_energy_out[0].value.abs()
                    + point.transport_energy_out[1].value.abs()
                    + particles.energy_in[0].value
                    + particles.energy_in[1].value
                    + particles.energy_out[0].value
                    + particles.energy_out[1].value
                    + self.initial_pool_energy[0].abs()
                    + self.initial_pool_energy[1].abs());
        if [
            heavy_residual,
            electron_residual,
            total_residual,
            baryon_residual,
            charge_residual,
        ]
        .iter()
        .any(|r| !r.is_finite())
            || [heavy_residual, electron_residual, total_residual]
                .iter()
                .any(|r| r.abs() > tolerance)
            || baryon_residual.abs() > self.controls.balance_relative_tolerance
            || charge_residual.abs() > self.controls.balance_relative_tolerance
        {
            return Err(DtFusionError::Balance(
                "source energy, baryon or charge account",
            ));
        }
        *balance = DtFusionBalance {
            particles: particles.balance(),
            reaction_extent_mol: point.extent.value,
            nuclear_energy_released_j: released,
            alpha_energy_deposited_j: point.alpha_energy.value,
            reacting_ion_energy_removed_j: point.ion_removal.value,
            neutron_energy_escaped_j: point.neutron_energy.value,
            electron_to_heavy_coulomb_energy_j: point.coulomb_transfer.value,
            photon_energy_escaped_j: point.photon_energy.value,
            heavy_transport_energy_out_j: point.transport_energy_out[0].value,
            electron_transport_energy_out_j: point.transport_energy_out[1].value,
            external_electron_energy_j: external_electron,
            external_heavy_energy_j: external_heavy,
            electron_energy_residual_j: electron_residual,
            heavy_energy_residual_j: heavy_residual,
            total_energy_residual_j: total_residual,
            relative_baryon_residual: baryon_residual,
            relative_charge_residual: charge_residual,
            ..*balance
        };
        poll(cancelled)?;
        Ok(state)
    }
}

/// Exact `dD/dt=dT/dt=-NA <sigma v> D T/V` at one frozen temperature.
/// `a` is the smaller stock, `delta` the conserved excess of the other fuel.
/// The direct remaining-stock formula preserves a trace after a nearly complete
/// burn; subtracting a clipped extent from a would erase that information.
pub(super) fn frozen_depletion(
    a: f64,
    delta: f64,
    rate: f64,
    volume: f64,
    dt: f64,
) -> Result<(f64, f64), DtFusionError> {
    if a == 0.0 || dt == 0.0 {
        return Ok((a, 0.0));
    }
    let hazard = ((a / volume) * AVOGADRO) * rate * dt;
    let z = ((delta / volume) * AVOGADRO) * rate * dt;
    if hazard.is_nan() || z.is_nan() || hazard < 0.0 || z < 0.0 {
        return Err(DtFusionError::Arithmetic("finite-fuel exposure"));
    }
    if hazard.is_infinite() || z.is_infinite() {
        return Ok((0.0, a));
    }
    let (remaining, extent) = if z < 1e-5 {
        let mean_exponential = if z == 0.0 { 1.0 } else { -expm1(-z) / z };
        let denominator = 1.0 + hazard * mean_exponential;
        (
            a * exp(-z) / denominator,
            a * ((hazard + z) * mean_exponential / denominator),
        )
    } else {
        let b = a + delta;
        let delta_fraction = delta / b;
        let burn_exponential = -expm1(-z);
        let survival = exp(-z);
        // Construct the denominator as numerator plus a nonnegative term;
        // separate a/b and delta/b can sum below one by rounding at large z.
        let denominator = burn_exponential + delta_fraction * survival;
        (
            a * (delta_fraction * survival / denominator),
            a * (burn_exponential / denominator),
        )
    };
    if !remaining.is_finite()
        || !extent.is_finite()
        || remaining < 0.0
        || extent < 0.0
        || remaining > a
        || extent > a
    {
        return Err(DtFusionError::Arithmetic("analytic finite-fuel map"));
    }
    Ok((remaining, extent))
}

fn plasma_error(error: DtFusionPlasmaErrorV1) -> DtFusionError {
    match error {
        DtFusionPlasmaErrorV1::Cancelled => DtFusionError::Cancelled,
        error => DtFusionError::Plasma(error),
    }
}
