use super::thermal::{PreparedThermal, ThermalStage};
use super::*;
use fs_math::det::expm1;
use fs_thermochem::equilibrium::composition::MAX_EQUILIBRIUM_SPECIES_V1;
use fs_thermochem::kinetics::atomic_ionization::photoionization::PhotoionizationEvaluationV1;

const CV: f64 = 1.5 * R;

#[derive(Clone, Copy, Default)]
pub(super) struct Sum {
    pub value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let increment = value - self.correction;
        let next = self.value + increment;
        self.correction = (next - self.value) - increment;
        self.value = next;
    }
}

#[derive(Clone, Copy, Default)]
pub(in crate::chemistry) struct Point {
    neutral: f64,
    minor: f64,
    te: f64,
    th: f64,
    pub(super) extent: Sum,
    photo: Sum,
    collisional: Sum,
    radiative: Sum,
    absorbed: Sum,
    captured: Sum,
    heavy_energy: Sum,
    exchange: Sum,
    neutral_exchange: Sum,
    coulomb_exchange: Sum,
    elapsed: Sum,
}
impl Point {
    pub fn neutral_amount(self) -> f64 {
        self.neutral
    }
    pub fn electron_temperature(self) -> f64 {
        self.te
    }
    pub fn heavy_temperature(self) -> f64 {
        self.th
    }
    pub fn reaction_extent_mol(self) -> f64 {
        self.extent.value
    }
}

#[derive(Clone, Copy)]
pub(in crate::chemistry) struct Rates {
    production: f64,
    destruction: f64,
    radiative_destruction: f64,
    captured_energy_j_mol: f64,
    thermal: Option<ThermalStage>,
}
impl Rates {
    pub fn inactive(self, point: Point) -> bool {
        (self.production == 0.0 || point.neutral == 0.0)
            && (self.destruction == 0.0 || point.minor == 0.0)
            && self.thermal.is_none_or(ThermalStage::inactive)
    }
}

pub(in crate::chemistry) struct Problem<'a> {
    initial: &'a FrozenIonizedGasStateV1,
    source: &'a AtomicRadiativeRecombinationV1,
    controls: CollisionalRadiativeSolve,
    photon: PhotoionizationEvaluationV1,
    indices: [usize; 3],
    excess: [f64; 2],
    total: f64,
    thermal: Option<PreparedThermal<'a>>,
}
impl<'a> Problem<'a> {
    pub fn new(
        initial: &'a FrozenIonizedGasStateV1,
        source: &'a AtomicRadiativeRecombinationV1,
        photons: &[PhotonNumberFluxV1],
        controls: CollisionalRadiativeSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, CollisionalRadiativeError> {
        Self::new_with_thermal(initial, source, photons, None, controls, cancelled)
    }

    pub fn new_with_thermal(
        initial: &'a FrozenIonizedGasStateV1,
        source: &'a AtomicRadiativeRecombinationV1,
        photons: &[PhotonNumberFluxV1],
        thermal: Option<&'a CollisionalRadiativeThermal>,
        controls: CollisionalRadiativeSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, CollisionalRadiativeError> {
        if let Some(inputs) = thermal {
            inputs.validate()?;
            if !(inputs.minimum_heavy_temperature_k..=inputs.maximum_heavy_temperature_k)
                .contains(&initial.heavy_temperature().value())
            {
                return Err(CollisionalRadiativeError::OutsideDomain);
            }
        }
        let channel = source.source().channel();
        if initial.thermodynamics() != channel.thermodynamics()
            || initial.electron_index() != Some(channel.electron_index())
            || initial.relative_charge_residual() > controls.relative_charge_tolerance
        {
            return Err(CollisionalRadiativeError::Invalid(
                "state and original atomic bank",
            ));
        }
        let indices = [
            channel.neutral_index(),
            channel.ion_index(),
            channel.electron_index(),
        ];
        if !(controls.minimum_electron_temperature_k..=controls.maximum_electron_temperature_k)
            .contains(&initial.electron_temperature().value())
        {
            return Err(CollisionalRadiativeError::OutsideDomain);
        }
        for (index, amount) in initial.amounts().iter().enumerate() {
            if !indices.contains(&index) && (thermal.is_none() || amount.value() == 0.0) {
                continue;
            }
            let (lo, hi) = if index == indices[2] {
                (
                    Temperature::new(controls.minimum_electron_temperature_k),
                    Temperature::new(controls.maximum_electron_temperature_k),
                )
            } else if let Some(inputs) = thermal {
                (
                    Temperature::new(inputs.minimum_heavy_temperature_k),
                    Temperature::new(inputs.maximum_heavy_temperature_k),
                )
            } else {
                (initial.heavy_temperature(), initial.heavy_temperature())
            };
            admit_caloric_interval(
                &channel.thermodynamics().models()[index],
                lo,
                hi,
                &mut *cancelled,
            )
            .map_err(|error| match error {
                KineticsErrorV1::Cancelled => CollisionalRadiativeError::Cancelled,
                error => CollisionalRadiativeError::Admission(error),
            })?;
        }
        let photon =
            source
                .source()
                .evaluate(photons, &mut *cancelled)
                .map_err(|error| match error {
                    PhotoionizationErrorV1::Cancelled => CollisionalRadiativeError::Cancelled,
                    error => CollisionalRadiativeError::Photoionization(error),
                })?;
        let ni = initial.amounts()[indices[1]].value();
        let ne = initial.amounts()[indices[2]].value();
        let total = initial.amounts()[indices[0]].value() + ni.min(ne);
        if !nonnegative(total) {
            return Err(CollisionalRadiativeError::Arithmetic("active inventory"));
        }
        let thermal = thermal
            .map(|inputs| PreparedThermal::new(initial, inputs, indices[0], controls, cancelled))
            .transpose()?;
        Ok(Self {
            initial,
            source,
            controls,
            photon,
            indices,
            excess: [(ni - ne).max(0.0), (ne - ni).max(0.0)],
            total,
            thermal,
        })
    }

    pub fn initial_point(&self) -> Point {
        Point {
            neutral: self.initial.amounts()[self.indices[0]].value(),
            minor: self.initial.amounts()[self.indices[1]]
                .value()
                .min(self.initial.amounts()[self.indices[2]].value()),
            te: self.initial.electron_temperature().value(),
            th: self.initial.heavy_temperature().value(),
            ..Point::default()
        }
    }

    pub fn rates(
        &self,
        point: Point,
        balance: &mut CollisionalRadiativeBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Rates, CollisionalRadiativeError> {
        poll(cancelled)?;
        let volume = self.initial.volume().value();
        let [ni, ne] = self.excess.map(|excess| point.minor + excess);
        // At exact electron absence use the photoelectron birth-temperature
        // limit, not an arbitrary stored temperature for an empty pool.
        let te = if ne == 0.0 && self.photon.rate_per_neutral_per_s > 0.0 {
            (NA / CV)
                * (self.photon.electron_heating_per_neutral_w / self.photon.rate_per_neutral_per_s)
        } else {
            point.te
        };
        self.admit_temperature(te)?;
        let channel = self.source.source().channel();
        let atomic = channel
            .evaluate(
                Temperature::new(point.th),
                Temperature::new(te),
                [point.neutral / volume, ni / volume, ne / volume],
                &mut *cancelled,
            )
            .map_err(|error| match error {
                AtomicIonizationErrorV1::Cancelled
                | AtomicIonizationErrorV1::Kinetics(KineticsErrorV1::Cancelled) => {
                    CollisionalRadiativeError::Cancelled
                }
                error => CollisionalRadiativeError::Atomic(error),
            })?;
        let radiative = self
            .source
            .evaluate(
                Temperature::new(te),
                self.controls.quadrature,
                &mut *cancelled,
            )
            .map_err(|error| match error {
                RadiativeRecombinationErrorV1::Cancelled => CollisionalRadiativeError::Cancelled,
                error => CollisionalRadiativeError::Recombination(error),
            })?;
        balance.maximum_relative_quadrature_error = balance
            .maximum_relative_quadrature_error
            .max(radiative.relative_quadrature_error_estimate);
        balance.maximum_relative_tail_bound = balance
            .maximum_relative_tail_bound
            .max(radiative.relative_tail_bound);
        // Evolve the SMALLER ion/electron inventory. Their conserved difference
        // may be nonzero because other original charged species are spectators.
        let three_body = if self.excess[0] == 0.0 {
            -atomic.concentration_derivatives_per_s[1]
        } else if point.minor > 0.0 {
            atomic.reverse_rate_mol_m3_s / (point.minor / volume)
        } else {
            0.0
        };
        let capture = (radiative.coefficient_m3_s * NA) * (ni.max(ne) / volume);
        let thermal = if let Some(thermal) = &self.thermal {
            let len = self.initial.amounts().len();
            let mut amounts = [Amount::new(0.0); MAX_EQUILIBRIUM_SPECIES_V1];
            amounts[..len].copy_from_slice(self.initial.amounts());
            for (index, amount) in self.indices.into_iter().zip([point.neutral, ni, ne]) {
                amounts[index] = Amount::new(amount);
            }
            Some(thermal.stage(&amounts[..len], point.th, te, cancelled)?)
        } else {
            None
        };
        let result = Rates {
            production: self.photon.rate_per_neutral_per_s
                + atomic.concentration_derivatives_per_s[0],
            destruction: three_body + capture,
            radiative_destruction: capture,
            captured_energy_j_mol: radiative.mean_captured_energy_j().unwrap_or(0.0) * NA,
            thermal,
        };
        if [
            result.production,
            result.destruction,
            capture,
            result.captured_energy_j_mol,
        ]
        .iter()
        .any(|v| !nonnegative(*v))
        {
            return Err(CollisionalRadiativeError::Arithmetic(
                "coupled production/capture coefficients",
            ));
        }
        Ok(result)
    }

    pub fn midpoint(
        &self,
        point: Point,
        rates: Rates,
        dt: f64,
        balance: &mut CollisionalRadiativeBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, CollisionalRadiativeError> {
        let predictor = self.frozen(point, rates, 0.5 * dt, cancelled)?;
        let middle = self.rates(predictor, balance, cancelled)?;
        self.frozen(point, middle, dt, cancelled)
    }

    fn frozen(
        &self,
        point: Point,
        rates: Rates,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, CollisionalRadiativeError> {
        if rates.inactive(point) || dt == 0.0 {
            return Ok(point);
        }
        let lambda = rates.production + rates.destruction;
        if !nonnegative(lambda) {
            return Err(CollisionalRadiativeError::Arithmetic("total turnover rate"));
        }
        let (neutral, minor, neutral_time, minor_time, extent) = if lambda == 0.0 {
            (
                point.neutral,
                point.minor,
                point.neutral * dt,
                point.minor * dt,
                0.0,
            )
        } else {
            let dose = lambda * dt;
            let depleted = if dose.is_finite() { -expm1(-dose) } else { 1.0 };
            let surviving = if dose.is_finite() { exp(-dose) } else { 0.0 };
            let relaxation = depleted / lambda;
            // dt - (1-exp(-lambda*dt))/lambda, without small-dose cancellation.
            let remainder = if dose < 1e-3 {
                dt * dose
                    * (0.5
                        + dose
                            * (-1.0 / 6.0
                                + dose * (1.0 / 24.0 + dose * (-1.0 / 120.0 + dose / 720.0))))
            } else {
                dt - relaxation
            };
            let birth = rates.production / lambda;
            let death = rates.destruction / lambda;
            let neutral_time = point.neutral * relaxation + (self.total * death) * remainder;
            let minor_time = point.minor * relaxation + (self.total * birth) * remainder;
            let mut neutral = point.neutral * surviving + (self.total * death) * depleted;
            let mut minor = point.minor * surviving + (self.total * birth) * depleted;
            // Use the exact invariant to reconstruct the larger side, keeping a
            // representable trace neutral or charged tail on the smaller side.
            if minor <= neutral {
                neutral = self.total - minor;
            } else {
                minor = self.total - neutral;
            }
            let extent = (birth * point.neutral - death * point.minor) * depleted;
            (neutral, minor, neutral_time, minor_time, extent)
        };
        let photo = self.photon.rate_per_neutral_per_s * neutral_time;
        let captured = rates.radiative_destruction * minor_time;
        // The coupled analytic transfer fixes NET collisional progress. Gross
        // forward/reverse traffic is never charged as chemical storage twice.
        let collisional = (extent - photo) + captured;
        let photo_heat = (self.photon.electron_heating_per_neutral_w * NA) * neutral_time;
        let absorbed = (self.photon.absorbed_power_per_neutral_w * NA) * neutral_time;
        let kinetic = captured * rates.captured_energy_j_mol;
        let mut sensible = photo_heat
            - self.source.source().channel().ionization_energy_j_mol() * collisional
            - kinetic;
        if [
            neutral,
            minor,
            neutral_time,
            minor_time,
            photo,
            captured,
            absorbed,
            kinetic,
        ]
        .iter()
        .any(|v| !nonnegative(*v))
            || !sensible.is_finite()
            || !collisional.is_finite()
        {
            return Err(CollisionalRadiativeError::Arithmetic(
                "coupled event/energy integral",
            ));
        }
        let electron = minor + self.excess[1];
        let initial_electron = point.minor + self.excess[1];
        let mut thermal_update = None;
        if let (Some(thermal), Some(stage)) = (&self.thermal, rates.thermal) {
            let update = thermal.advance(
                stage,
                [point.th, point.te],
                initial_electron,
                point.heavy_energy.value,
                sensible,
                extent,
                dt,
                cancelled,
            )?;
            sensible = update.electron_energy;
            thermal_update = Some(update);
        }
        if !nonnegative(electron)
            || (electron == 0.0 && (initial_electron != 0.0 || sensible != 0.0))
        {
            return Err(CollisionalRadiativeError::OutsideDomain);
        }
        let te = if electron == 0.0 {
            point.te
        } else if initial_electron == 0.0 {
            sensible / (CV * electron)
        } else {
            let increment = (sensible - CV * point.te * extent) / (CV * electron);
            if increment.is_finite() && increment.abs() < 0.5 * point.te {
                point.te + increment
            } else {
                (CV * initial_electron * point.te + sensible) / (CV * electron)
            }
        };
        self.admit_temperature(te)?;
        let mut next = Point {
            neutral,
            minor,
            te,
            ..point
        };
        next.extent.add(extent);
        next.photo.add(photo);
        next.collisional.add(collisional);
        next.radiative.add(captured);
        next.absorbed.add(absorbed);
        next.captured.add(kinetic);
        if let Some(update) = thermal_update {
            next.th = update.heavy_temperature;
            next.heavy_energy.add(update.heavy_increment);
            next.exchange.add(update.exchange);
            next.neutral_exchange.add(update.neutral_exchange);
            next.coulomb_exchange.add(update.coulomb_exchange);
            next.elapsed.add(dt);
        }
        if [
            next.extent.value,
            next.collisional.value,
            next.heavy_energy.value,
            next.exchange.value,
            next.neutral_exchange.value,
            next.coulomb_exchange.value,
            next.elapsed.value,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || [
                next.photo.value,
                next.radiative.value,
                next.absorbed.value,
                next.captured.value,
            ]
            .iter()
            .any(|v| !nonnegative(*v))
        {
            return Err(CollisionalRadiativeError::Arithmetic(
                "cumulative event/energy integral",
            ));
        }
        Ok(next)
    }

    pub fn error(
        &self,
        start: Point,
        full: Point,
        next: Point,
    ) -> Result<f64, CollisionalRadiativeError> {
        let amount = self.controls.absolute_amount_tolerance_mol;
        let relative = self.controls.relative_tolerance;
        let mut error = ((full.neutral - next.neutral).abs()
            / (amount + relative * start.neutral.max(next.neutral)))
        .max((full.minor - next.minor).abs() / (amount + relative * start.minor.max(next.minor)))
        .max(
            (full.te - next.te).abs()
                / (self.controls.temperature_tolerance_k + relative * full.te.max(next.te)),
        )
        .max(
            (full.th - next.th).abs()
                / (self.controls.temperature_tolerance_k + relative * full.th.max(next.th)),
        );
        for (old, coarse, fine) in [
            (start.photo.value, full.photo.value, next.photo.value),
            (
                start.collisional.value,
                full.collisional.value,
                next.collisional.value,
            ),
            (
                start.radiative.value,
                full.radiative.value,
                next.radiative.value,
            ),
        ] {
            // Event integrals are augmented ODE states. Normalizing by only
            // this step's O(dt) increment would impose a spurious dt² relative
            // restriction instead of controlling their cumulative state error.
            error = error.max(
                (coarse - fine).abs()
                    / (amount + relative * old.abs().max(coarse.abs()).max(fine.abs())),
            );
        }
        for (old, coarse, fine) in [
            (
                start.absorbed.value,
                full.absorbed.value,
                next.absorbed.value,
            ),
            (
                start.captured.value,
                full.captured.value,
                next.captured.value,
            ),
            (
                start.heavy_energy.value,
                full.heavy_energy.value,
                next.heavy_energy.value,
            ),
            (
                start.exchange.value,
                full.exchange.value,
                next.exchange.value,
            ),
            (
                start.coulomb_exchange.value,
                full.coulomb_exchange.value,
                next.coulomb_exchange.value,
            ),
            (
                start.neutral_exchange.value,
                full.neutral_exchange.value,
                next.neutral_exchange.value,
            ),
        ] {
            error = error.max(
                (coarse - fine).abs()
                    / (self.controls.energy_tolerance_j
                        + relative * old.abs().max(coarse.abs()).max(fine.abs())),
            );
        }
        if !nonnegative(error) {
            return Err(CollisionalRadiativeError::Arithmetic(
                "local error estimate",
            ));
        }
        // A second-order full step has four times the leading error of two
        // half steps; their difference is three times the retained fine error.
        Ok(error / 3.0)
    }

    pub(in crate::chemistry) fn admit_temperature(
        &self,
        te: f64,
    ) -> Result<(), CollisionalRadiativeError> {
        if !positive(te)
            || te < self.controls.minimum_electron_temperature_k
            || te > self.controls.maximum_electron_temperature_k
        {
            Err(CollisionalRadiativeError::OutsideDomain)
        } else {
            Ok(())
        }
    }

    /// Apply one admitted photon increment supplied by the spatial radiation
    /// owner. Both energies are TOTAL joules. No time or heavy-pool work is
    /// advanced here. Interval admission may follow the complete spectral
    /// composition; every intermediate electron temperature must be positive.
    pub fn absorb_photons(
        &self,
        point: Point,
        remaining_neutral_mol: f64,
        absorbed_mol: f64,
        excess_energy_j: f64,
        absorbed_photon_energy_j: f64,
    ) -> Result<Point, CollisionalRadiativeError> {
        if [
            remaining_neutral_mol,
            absorbed_mol,
            excess_energy_j,
            absorbed_photon_energy_j,
        ]
        .iter()
        .any(|v| !nonnegative(*v))
            || remaining_neutral_mol > point.neutral
            || absorbed_mol > point.neutral
        {
            return Err(CollisionalRadiativeError::Invalid("photon increment"));
        }
        if absorbed_mol == 0.0 {
            if remaining_neutral_mol != point.neutral
                || excess_energy_j != 0.0
                || absorbed_photon_energy_j != 0.0
            {
                return Err(CollisionalRadiativeError::Invalid(
                    "photon energy without absorption",
                ));
            }
            return Ok(point);
        }
        let amount_allowance =
            self.controls.absolute_amount_tolerance_mol + 64.0 * f64::EPSILON * point.neutral;
        if ((point.neutral - remaining_neutral_mol) - absorbed_mol).abs() > amount_allowance {
            return Err(CollisionalRadiativeError::Invalid(
                "photon/neutral inventory disagreement",
            ));
        }
        let chi = self.source.source().channel().ionization_energy_j_mol();
        let energy_residual = excess_energy_j + chi * absorbed_mol - absorbed_photon_energy_j;
        if !energy_residual.is_finite() || energy_residual.abs() > self.controls.energy_tolerance_j
        {
            return Err(CollisionalRadiativeError::EnergyBalance {
                residual_j: energy_residual,
            });
        }
        let minor = if remaining_neutral_mol <= point.minor + absorbed_mol {
            self.total - remaining_neutral_mol
        } else {
            point.minor + absorbed_mol
        };
        let electron = minor + self.excess[1];
        let increment = (excess_energy_j - CV * point.te * absorbed_mol) / (CV * electron);
        let te = if point.minor + self.excess[1] == 0.0 {
            excess_energy_j / (CV * electron)
        } else if increment.is_finite() && increment.abs() < 0.5 * point.te {
            point.te + increment
        } else {
            (CV * (point.minor + self.excess[1]) * point.te + excess_energy_j) / (CV * electron)
        };
        if !nonnegative(minor) || !positive(electron) || !positive(te) {
            return Err(CollisionalRadiativeError::OutsideDomain);
        }
        let mut next = Point {
            neutral: remaining_neutral_mol,
            minor,
            te,
            ..point
        };
        next.extent.add(absorbed_mol);
        next.photo.add(absorbed_mol);
        next.absorbed.add(absorbed_photon_energy_j);
        if [next.extent.value, next.photo.value, next.absorbed.value]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(CollisionalRadiativeError::Arithmetic(
                "cumulative photon increment",
            ));
        }
        Ok(next)
    }

    pub fn publish(
        &self,
        point: Point,
        balance: &mut CollisionalRadiativeBalance,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FrozenIonizedGasStateV1, CollisionalRadiativeError> {
        let chi = self.source.source().channel().ionization_energy_j_mol();
        let endpoint = [
            point.neutral,
            point.minor + self.excess[0],
            point.minor + self.excess[1],
        ];
        let mut changes = [0.0; 3];
        for (slot, index) in self.indices.into_iter().enumerate() {
            let before = self.initial.amounts()[index].value();
            changes[slot] = endpoint[slot] - before;
            let expected = if slot == 0 {
                -point.extent.value
            } else {
                point.extent.value
            };
            let allowance = self.controls.absolute_amount_tolerance_mol
                + 64.0 * f64::EPSILON * before.max(endpoint[slot]);
            if (changes[slot] - expected).abs() > allowance {
                return Err(CollisionalRadiativeError::Arithmetic(
                    "published population and reaction extent disagree",
                ));
            }
        }
        let initial_ne = self.initial.amounts()[self.indices[2]].value();
        // Use actual endpoint population differences, independently of the
        // integrated reaction counters. Conserved formation-reference shifts
        // cancel on the admitted atomic stoichiometry before this subtraction.
        let sensible = CV
            * (initial_ne * (point.te - self.initial.electron_temperature().value())
                + changes[2] * point.te);
        let chemical = chi * changes[1];
        let emitted = chi * point.radiative.value + point.captured.value;
        // The admitted ground-state channel gives neutral and ion the same
        // heavy molar Cv. Their total is conserved while every spectator keeps
        // its amount, so the ORIGINAL heavy caloric difference remains valid.
        // Reconstruct it from the published temperature, then account for any
        // representable endpoint roundoff in the active heavy total as well.
        let heavy = if let Some(thermal) = &self.thermal {
            thermal.calorics(point.th, cancelled)?.0
        } else {
            0.0
        } + CV * point.th * (changes[0] + changes[1]);
        let [heavy_work, electron_work] = self.thermal.as_ref().map_or([0.0; 2], |thermal| {
            thermal.powers().map(|power| power * point.elapsed.value)
        });
        let mut electron_residual = Sum::default();
        for value in [
            sensible,
            chemical,
            emitted,
            -point.absorbed.value,
            point.exchange.value,
            -electron_work,
        ] {
            electron_residual.add(value);
        }
        let mut heavy_residual = Sum::default();
        for value in [heavy, -point.exchange.value, -heavy_work] {
            heavy_residual.add(value);
        }
        let residual = electron_residual.value + heavy_residual.value;
        if [
            residual,
            electron_residual.value,
            heavy_residual.value,
            heavy,
            heavy_work,
            electron_work,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || !nonnegative(emitted)
            || !nonnegative(point.photo.value * NA)
            || !nonnegative(point.radiative.value * NA)
        {
            return Err(CollisionalRadiativeError::Arithmetic(
                "final photon/gas energy account",
            ));
        }
        for account in [electron_residual.value, heavy_residual.value, residual] {
            if account.abs() > self.controls.energy_tolerance_j {
                return Err(CollisionalRadiativeError::EnergyBalance {
                    residual_j: account,
                });
            }
        }
        let state = if changes == [0.0; 3]
            && point.te == self.initial.electron_temperature().value()
            && point.th == self.initial.heavy_temperature().value()
        {
            self.initial.clone()
        } else {
            let mut amounts = self.initial.amounts().to_vec();
            for (slot, index) in self.indices.into_iter().enumerate() {
                amounts[index] = Amount::new(endpoint[slot]);
            }
            FrozenIonizedGasStateV1::new(
                self.initial.thermodynamics(),
                &amounts,
                self.initial.volume(),
                Temperature::new(point.th),
                Temperature::new(point.te),
                self.controls.relative_charge_tolerance,
                &mut *cancelled,
            )
            .map_err(|error| match error {
                FrozenIonizedGasStateErrorV1::Cancelled => CollisionalRadiativeError::Cancelled,
                error => CollisionalRadiativeError::State(error),
            })?
        };
        balance.photoionization_extent_mol = point.photo.value;
        balance.collisional_extent_mol = point.collisional.value;
        balance.radiative_recombination_extent_mol = point.radiative.value;
        balance.absorbed_photons = NA * point.photo.value;
        balance.emitted_photons = NA * point.radiative.value;
        balance.absorbed_radiation_energy_j = point.absorbed.value;
        balance.emitted_radiation_energy_j = emitted;
        balance.chemical_storage_change_j = chemical;
        balance.electron_sensible_energy_change_j = sensible;
        balance.heavy_sensible_energy_change_j = heavy;
        balance.electron_deposited_energy_j = electron_work;
        balance.heavy_deposited_energy_j = heavy_work;
        balance.electron_to_heavy_energy_j = point.exchange.value;
        balance.electron_to_heavy_neutral_energy_j = point.neutral_exchange.value;
        balance.electron_to_heavy_coulomb_energy_j = point.coulomb_exchange.value;
        balance.electron_sensible_energy_residual_j = electron_residual.value;
        balance.heavy_sensible_energy_residual_j = heavy_residual.value;
        balance.energy_residual_j = residual;
        Ok(state)
    }
}
