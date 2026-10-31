use super::*;

#[derive(Clone, Copy)]
pub(super) struct Point {
    pub time: f64,
    pub exchange: f64,
    pub temperatures: [f64; 2],
    pub capacities: [f64; 2],
}

pub(super) struct Problem<'a> {
    initial: &'a FrozenIonizedGasStateV1,
    sources: [Option<&'a ElectronCrossSectionV1>; MAX_EQUILIBRIUM_SPECIES_V1],
    powers: [f64; 2],
    controls: ElectronHeavyEnergySolve,
    energy_limits: [[f64; 2]; 2],
    has_electrons: bool,
}
impl<'a> Problem<'a> {
    pub fn new(
        initial: &'a FrozenIonizedGasStateV1,
        sources: &'a [ElectronCrossSectionV1],
        powers: [f64; 2],
        controls: ElectronHeavyEnergySolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ElectronHeavyEnergyError> {
        let bank = initial.thermodynamics();
        let species = bank.elemental_matrix().species();
        if initial.relative_charge_residual() > controls.relative_charge_tolerance {
            return Err(ElectronHeavyEnergyError::Invalid(
                "initial charge residual exceeds publication tolerance",
            ));
        }
        let has_electrons = initial
            .electron_index()
            .is_some_and(|i| initial.amounts()[i].value() > 0.0);
        if sources.len() > species.len()
            || (!has_electrons && (!sources.is_empty() || powers[1] != 0.0))
        {
            return Err(ElectronHeavyEnergyError::Invalid(
                "absent electrons cannot receive power or require collision sources",
            ));
        }
        let mut result = Self {
            initial,
            sources: [None; MAX_EQUILIBRIUM_SPECIES_V1],
            powers,
            controls,
            energy_limits: [[0.0; 2]; 2],
            has_electrons,
        };
        for source in sources {
            let i = species
                .binary_search(source.target_species())
                .map_err(|_| {
                    ElectronHeavyEnergyError::Invalid("collision target absent from source bank")
                })?;
            if source.kind() != ElectronCollisionKindV1::ElasticMomentumTransfer
                || bank.charge_vector().get(i) != Some(0)
                || initial.amounts()[i].value() == 0.0
                || result.sources[i].is_some()
            {
                return Err(ElectronHeavyEnergyError::Invalid(
                    "each present neutral target needs one elastic momentum source",
                ));
            }
            result.sources[i] = Some(source);
        }
        let mut total = 0.0;
        let mut charged = 0.0;
        for (i, amount) in initial.amounts().iter().enumerate() {
            if amount.value() == 0.0 {
                continue;
            }
            total += amount.value();
            if bank.charge_vector().get(i) != Some(0) {
                charged += amount.value();
            } else if has_electrons && result.sources[i].is_none() {
                return Err(ElectronHeavyEnergyError::MissingSource(species[i].clone()));
            }
            let t = initial
                .species_temperature(i)
                .expect("canonical source index")
                .value();
            if t < controls.minimum_temperature_k || t > controls.maximum_temperature_k {
                return Err(ElectronHeavyEnergyError::OutsideDomain);
            }
            admit_caloric_interval(
                &bank.models()[i],
                Temperature::new(controls.minimum_temperature_k),
                Temperature::new(controls.maximum_temperature_k),
                &mut *cancelled,
            )
            .map_err(|error| match error {
                KineticsErrorV1::Cancelled => ElectronHeavyEnergyError::Cancelled,
                error => ElectronHeavyEnergyError::Caloric(error),
            })?;
        }
        if !positive(total) || charged / total > controls.maximum_charged_mole_fraction {
            return Err(ElectronHeavyEnergyError::Invalid(
                "electron-neutral dilute charged-fraction domain",
            ));
        }
        if has_electrons {
            let i = initial.electron_index().expect("present electron index");
            let mass = bank.models()[i].molar_mass().value() / AVOGADRO_PER_MOL;
            if ((mass - ELECTRON_MASS_KG) / ELECTRON_MASS_KG).abs() > 16.0 * f64::EPSILON {
                return Err(ElectronHeavyEnergyError::Invalid(
                    "source electron mass differs from Maxwellian collision mass",
                ));
            }
        }
        for pool in 0..2 {
            result.energy_limits[pool] = [
                result
                    .calorics(pool, controls.minimum_temperature_k, cancelled)?
                    .0,
                result
                    .calorics(pool, controls.maximum_temperature_k, cancelled)?
                    .0,
            ];
        }
        Ok(result)
    }

    pub fn targets(&self, time: f64, exchange: f64) -> Result<[f64; 2], ElectronHeavyEnergyError> {
        let target = [
            self.powers[0] * time + exchange,
            self.powers[1] * time - exchange,
        ];
        if target.iter().any(|v| !v.is_finite()) {
            return Err(ElectronHeavyEnergyError::Arithmetic("pool energy targets"));
        }
        Ok(target)
    }

    pub fn calorics(
        &self,
        pool: usize,
        temperature: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), ElectronHeavyEnergyError> {
        let t = Temperature::new(temperature);
        let reference = if pool == 1 {
            self.initial.electron_temperature()
        } else {
            self.initial.heavy_temperature()
        };
        let mut energy = Sum::default();
        let mut cv = Sum::default();
        for (i, (amount, model)) in self
            .initial
            .amounts()
            .iter()
            .zip(self.initial.thermodynamics().models())
            .enumerate()
        {
            if amount.value() == 0.0
                || usize::from(Some(i) == self.initial.electron_index()) != pool
            {
                continue;
            }
            poll(cancelled)?;
            let error = |source| ElectronHeavyEnergyError::Species {
                species: model.species().clone(),
                source,
            };
            let increment = model
                .molar_internal_energy_change(reference, t)
                .map_err(error)?
                .value();
            let capacity = model.evaluate_properties(t).map_err(error)?.cp().value() - R;
            energy.add(amount.value() * increment);
            cv.add(amount.value() * capacity);
        }
        if !energy.value.is_finite()
            || !cv.value.is_finite()
            || (cv.value <= 0.0 && (pool == 0 || self.has_electrons))
        {
            return Err(ElectronHeavyEnergyError::Arithmetic(
                "source pool sensible energy or heat capacity",
            ));
        }
        Ok((energy.value, cv.value))
    }

    fn invert(
        &self,
        pool: usize,
        target: f64,
        guess: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), ElectronHeavyEnergyError> {
        if pool == 1 && !self.has_electrons {
            return if target == 0.0 {
                Ok((self.initial.electron_temperature().value(), 0.0))
            } else {
                Err(ElectronHeavyEnergyError::Invalid(
                    "energy supplied to absent electrons",
                ))
            };
        }
        let [minimum, maximum] = self.energy_limits[pool];
        if target < minimum || target > maximum {
            return Err(ElectronHeavyEnergyError::OutsideDomain);
        }
        let mut lower = self.controls.minimum_temperature_k;
        let mut upper = self.controls.maximum_temperature_k;
        let mut t = if (lower..=upper).contains(&guess) {
            guess
        } else {
            lower + 0.5 * (upper - lower)
        };
        for _ in 0..self.controls.maximum_iterations {
            let (energy, cv) = self.calorics(pool, t, cancelled)?;
            let residual = energy - target;
            // A trace electron pool can store far less than the parcel's
            // absolute energy tolerance. Resolve its temperature as well:
            // otherwise an energetically tiny but very large Te change could
            // be accepted as exactly zero by both step-doubling paths.
            let temperature_scale = self.controls.temperature_tolerance_k
                + self.controls.relative_tolerance * t.abs().max(guess.abs());
            let tolerance =
                (self.controls.energy_tolerance_j / 32.0).min(cv * (temperature_scale / 64.0));
            if residual.abs() <= tolerance {
                return Ok((t, cv));
            }
            if residual > 0.0 {
                upper = t;
            } else {
                lower = t;
            }
            let next = t - residual / cv;
            t = if next > lower && next < upper {
                next
            } else {
                lower + 0.5 * (upper - lower)
            };
        }
        Err(ElectronHeavyEnergyError::Nonconvergence)
    }

    pub fn coupling(
        &self,
        point: &Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, ElectronHeavyEnergyError> {
        if !self.has_electrons {
            return Ok(0.0);
        }
        let mut frequency = Sum::default();
        let te = point.temperatures[1];
        let th = point.temperatures[0];
        for (i, source) in self
            .sources
            .iter()
            .enumerate()
            .filter_map(|(i, s)| s.map(|s| (i, s)))
        {
            poll(cancelled)?;
            let mass = self.initial.thermodynamics().models()[i]
                .molar_mass()
                .value()
                / AVOGADRO_PER_MOL;
            if fs_math::det::sqrt(ELECTRON_MASS_KG / mass)
                > self.controls.maximum_target_thermal_speed_ratio
            {
                return Err(ElectronHeavyEnergyError::OutsideDomain);
            }
            let ratio_squared = (th / te) * (ELECTRON_MASS_KG / mass);
            if !ratio_squared.is_finite()
                || ratio_squared > self.controls.maximum_target_thermal_speed_ratio.powi(2)
            {
                return Err(ElectronHeavyEnergyError::OutsideDomain);
            }
            let evaluation = source
                .evaluate_properties(
                    Temperature::new(te),
                    MaxwellianRateWeightingV1::WeakDriftMomentum,
                    self.controls.rate_controls,
                    &mut *cancelled,
                )
                .map_err(|error| match error {
                    ElectronCollisionErrorV1::Cancelled => ElectronHeavyEnergyError::Cancelled,
                    error => ElectronHeavyEnergyError::Collision {
                        species: source.target_species().clone(),
                        source: error,
                    },
                })?;
            let coefficient =
                evaluation
                    .coefficient_m3_s()
                    .ok_or(ElectronHeavyEnergyError::Arithmetic(
                        "unrepresentable elastic momentum coefficient",
                    ))?;
            let density = self
                .initial
                .number_density_per_m3(i)
                .expect("canonical target index");
            frequency.add((density * coefficient) * (ELECTRON_MASS_KG / mass));
        }
        let electron_amount = self.initial.amounts()
            [self.initial.electron_index().expect("present electron")]
        .value();
        // R=kB*NA is exact in the shared SI convention. Keeping the molar
        // amount avoids overflowing an unnecessary total particle count.
        let g = (3.0 * R * electron_amount) * frequency.value;
        if !g.is_finite() || g < 0.0 {
            return Err(ElectronHeavyEnergyError::Arithmetic(
                "elastic electron-heavy coupling",
            ));
        }
        Ok(g)
    }

    pub fn step(
        &self,
        start: &Point,
        end: f64,
        coupling: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, ElectronHeavyEnergyError> {
        let alpha = (end - start.time) * coupling;
        if !alpha.is_finite() {
            return Err(ElectronHeavyEnergyError::Arithmetic(
                "implicit exchange coefficient",
            ));
        }
        let external = self.targets(end, 0.0)?;
        let mut lower =
            (self.energy_limits[0][0] - external[0]).max(external[1] - self.energy_limits[1][1]);
        let mut upper =
            (self.energy_limits[0][1] - external[0]).min(external[1] - self.energy_limits[1][0]);
        if lower > upper {
            return Err(ElectronHeavyEnergyError::OutsideDomain);
        }
        let mut exchange = if (lower..=upper).contains(&start.exchange) {
            start.exchange
        } else {
            lower + 0.5 * (upper - lower)
        };
        for _ in 0..self.controls.maximum_iterations {
            poll(cancelled)?;
            let targets = self.targets(end, exchange)?;
            let (th, ch) = self.invert(0, targets[0], start.temperatures[0], cancelled)?;
            let (te, ce) = self.invert(1, targets[1], start.temperatures[1], cancelled)?;
            let increment = exchange - start.exchange;
            let difference = te - th;
            let residual = increment - alpha * difference;
            if !residual.is_finite() {
                return Err(ElectronHeavyEnergyError::Arithmetic("exchange residual"));
            }
            let temperature_scale = self.controls.temperature_tolerance_k
                + self.controls.relative_tolerance * th.abs().min(te.abs());
            let capacity_scale = if ce == 0.0 { ch } else { ch.min(ce) };
            let tolerance = (self.controls.energy_tolerance_j / 8.0)
                .min(capacity_scale * (temperature_scale / 64.0));
            if residual.abs() <= tolerance {
                if increment * difference < -tolerance * difference.abs() {
                    return Err(ElectronHeavyEnergyError::Balance(
                        "elastic heat transfer against temperature difference",
                    ));
                }
                return Ok(Point {
                    time: end,
                    exchange,
                    temperatures: [th, te],
                    capacities: [ch, ce],
                });
            }
            if residual > 0.0 {
                upper = exchange;
            } else {
                lower = exchange;
            }
            let derivative = 1.0
                + if alpha == 0.0 {
                    0.0
                } else {
                    alpha * (1.0 / ch + 1.0 / ce)
                };
            let next = exchange - residual / derivative;
            exchange = if next > lower && next < upper {
                next
            } else {
                lower + 0.5 * (upper - lower)
            };
        }
        Err(ElectronHeavyEnergyError::Nonconvergence)
    }
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, x: f64) {
        let y = x - self.correction;
        let sum = self.value + y;
        self.correction = (sum - self.value) - y;
        self.value = sum;
    }
}
