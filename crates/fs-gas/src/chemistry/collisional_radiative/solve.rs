use super::*;
use fs_math::det::expm1;
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
pub(super) struct Point {
    neutral: f64,
    minor: f64,
    te: f64,
    pub extent: Sum,
    photo: Sum,
    collisional: Sum,
    radiative: Sum,
    absorbed: Sum,
    captured: Sum,
}

#[derive(Clone, Copy)]
pub(super) struct Rates {
    production: f64,
    destruction: f64,
    radiative_destruction: f64,
    captured_energy_j_mol: f64,
}
impl Rates {
    pub fn inactive(self, point: Point) -> bool {
        (self.production == 0.0 || point.neutral == 0.0)
            && (self.destruction == 0.0 || point.minor == 0.0)
    }
}

pub(super) struct Problem<'a> {
    initial: &'a FrozenIonizedGasStateV1,
    source: &'a AtomicRadiativeRecombinationV1,
    controls: CollisionalRadiativeSolve,
    photon: PhotoionizationEvaluationV1,
    indices: [usize; 3],
    excess: [f64; 2],
    total: f64,
}
impl<'a> Problem<'a> {
    pub fn new(
        initial: &'a FrozenIonizedGasStateV1,
        source: &'a AtomicRadiativeRecombinationV1,
        photons: &[PhotonNumberFluxV1],
        controls: CollisionalRadiativeSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, CollisionalRadiativeError> {
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
        for index in indices {
            let (lo, hi) = if index == indices[2] {
                (
                    Temperature::new(controls.minimum_electron_temperature_k),
                    Temperature::new(controls.maximum_electron_temperature_k),
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
        Ok(Self {
            initial,
            source,
            controls,
            photon,
            indices,
            excess: [(ni - ne).max(0.0), (ne - ni).max(0.0)],
            total,
        })
    }

    pub fn initial_point(&self) -> Point {
        Point {
            neutral: self.initial.amounts()[self.indices[0]].value(),
            minor: self.initial.amounts()[self.indices[1]]
                .value()
                .min(self.initial.amounts()[self.indices[2]].value()),
            te: self.initial.electron_temperature().value(),
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
                self.initial.heavy_temperature(),
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
        let result = Rates {
            production: self.photon.rate_per_neutral_per_s
                + atomic.concentration_derivatives_per_s[0],
            destruction: three_body + capture,
            radiative_destruction: capture,
            captured_energy_j_mol: radiative.mean_captured_energy_j().unwrap_or(0.0) * NA,
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
        let predictor = self.frozen(point, rates, 0.5 * dt)?;
        let middle = self.rates(predictor, balance, cancelled)?;
        self.frozen(point, middle, dt)
    }

    fn frozen(
        &self,
        point: Point,
        rates: Rates,
        dt: f64,
    ) -> Result<Point, CollisionalRadiativeError> {
        if rates.inactive(point) || dt == 0.0 {
            return Ok(point);
        }
        let lambda = rates.production + rates.destruction;
        if !positive(lambda) {
            return Err(CollisionalRadiativeError::Arithmetic("total turnover rate"));
        }
        let dose = lambda * dt;
        let depleted = if dose.is_finite() { -expm1(-dose) } else { 1.0 };
        let surviving = if dose.is_finite() { exp(-dose) } else { 0.0 };
        let relaxation = depleted / lambda;
        // dt - (1-exp(-lambda*dt))/lambda, without small-dose cancellation.
        let remainder = if dose < 1e-3 {
            dt * dose
                * (0.5
                    + dose
                        * (-1.0 / 6.0 + dose * (1.0 / 24.0 + dose * (-1.0 / 120.0 + dose / 720.0))))
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
        let photo = self.photon.rate_per_neutral_per_s * neutral_time;
        let captured = rates.radiative_destruction * minor_time;
        // The coupled analytic transfer fixes NET collisional progress. Gross
        // forward/reverse traffic is never charged as chemical storage twice.
        let collisional = (extent - photo) + captured;
        let photo_heat = (self.photon.electron_heating_per_neutral_w * NA) * neutral_time;
        let absorbed = (self.photon.absorbed_power_per_neutral_w * NA) * neutral_time;
        let kinetic = captured * rates.captured_energy_j_mol;
        let sensible = photo_heat
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
        if !positive(electron) {
            return Err(CollisionalRadiativeError::OutsideDomain);
        }
        let initial_electron = point.minor + self.excess[1];
        let increment = (sensible - CV * point.te * extent) / (CV * electron);
        let te = if initial_electron == 0.0 {
            sensible / (CV * electron)
        } else if increment.is_finite() && increment.abs() < 0.5 * point.te {
            point.te + increment
        } else {
            (CV * initial_electron * point.te + sensible) / (CV * electron)
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
        if [next.extent.value, next.collisional.value]
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

    fn admit_temperature(&self, te: f64) -> Result<(), CollisionalRadiativeError> {
        if !positive(te)
            || te < self.controls.minimum_electron_temperature_k
            || te > self.controls.maximum_electron_temperature_k
        {
            Err(CollisionalRadiativeError::OutsideDomain)
        } else {
            Ok(())
        }
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
        let residual = sensible + chemical + emitted - point.absorbed.value;
        if !residual.is_finite()
            || !nonnegative(emitted)
            || !nonnegative(point.photo.value * NA)
            || !nonnegative(point.radiative.value * NA)
        {
            return Err(CollisionalRadiativeError::Arithmetic(
                "final photon/gas energy account",
            ));
        }
        if residual.abs() > self.controls.energy_tolerance_j {
            return Err(CollisionalRadiativeError::EnergyBalance {
                residual_j: residual,
            });
        }
        let state =
            if changes == [0.0; 3] && point.te == self.initial.electron_temperature().value() {
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
                    self.initial.heavy_temperature(),
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
        balance.energy_residual_j = residual;
        Ok(state)
    }
}
