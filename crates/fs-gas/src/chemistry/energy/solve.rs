use super::*;

/// Original-population sensible energy change and heat capacity for one pool.
pub(in crate::chemistry) fn source_pool_calorics(
    initial: &FrozenIonizedGasStateV1,
    electrons: bool,
    temperature: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64), ElectronHeavyEnergyError> {
    let t = Temperature::new(temperature);
    let reference = if electrons {
        initial.electron_temperature()
    } else {
        initial.heavy_temperature()
    };
    let mut energy = Sum::default();
    let mut cv = Sum::default();
    for (i, (amount, model)) in initial
        .amounts()
        .iter()
        .zip(initial.thermodynamics().models())
        .enumerate()
    {
        if amount.value() == 0.0 || (Some(i) == initial.electron_index()) != electrons {
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
    let has_electrons = initial
        .electron_index()
        .is_some_and(|i| initial.amounts()[i].value() > 0.0);
    if !energy.value.is_finite()
        || !cv.value.is_finite()
        || (cv.value <= 0.0 && (!electrons || has_electrons))
    {
        return Err(ElectronHeavyEnergyError::Arithmetic(
            "source pool sensible energy or heat capacity",
        ));
    }
    Ok((energy.value, cv.value))
}

#[derive(Clone, Copy)]
pub(super) struct Point {
    pub time: f64,
    pub exchange: f64,
    pub coulomb_exchange: f64,
    pub temperatures: [f64; 2],
    pub capacities: [f64; 2],
}

#[derive(Clone, Copy)]
pub(super) struct Coupling {
    neutral: f64,
    coulomb: f64,
}
impl Coupling {
    fn total(self) -> f64 {
        self.neutral + self.coulomb
    }
}

pub(super) struct Problem<'a> {
    initial: &'a FrozenIonizedGasStateV1,
    sources: [Option<&'a ElectronCrossSectionV1>; MAX_EQUILIBRIUM_SPECIES_V1],
    coulomb: Option<&'a CoulombModelV1>,
    powers: [f64; 2],
    controls: ElectronHeavyEnergySolve,
    energy_limits: [[f64; 2]; 2],
    has_electrons: bool,
}
impl<'a> Problem<'a> {
    pub fn new(
        initial: &'a FrozenIonizedGasStateV1,
        sources: &'a [ElectronCrossSectionV1],
        coulomb: Option<&'a CoulombModelV1>,
        powers: [f64; 2],
        controls: ElectronHeavyEnergySolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ElectronHeavyEnergyError> {
        let bank = initial.thermodynamics();
        let species = bank.elemental_matrix().species();
        if coulomb.is_some_and(|model| !model.is_bound_to_populations(initial)) {
            return Err(ElectronHeavyEnergyError::Invalid(
                "Coulomb law must retain the same source models, amounts and volume",
            ));
        }
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
            coulomb,
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
                "declared charged-fraction domain",
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
        source_pool_calorics(self.initial, pool == 1, temperature, cancelled)
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
    ) -> Result<Coupling, ElectronHeavyEnergyError> {
        let te = point.temperatures[1];
        let th = point.temperatures[0];
        let coulomb = if let Some(model) = self.coulomb {
            let properties = model
                .evaluate(Temperature::new(th), Temperature::new(te), &mut *cancelled)
                .map_err(coulomb_error)?;
            properties.electron_heavy_thermal_conductance_w_m3_k() * self.initial.volume().value()
        } else {
            0.0
        };
        if !coulomb.is_finite() || coulomb < 0.0 {
            return Err(ElectronHeavyEnergyError::Arithmetic(
                "Coulomb thermal conductance",
            ));
        }
        let g = super::exchange::neutral_conductance(
            self.initial,
            &self.sources,
            th,
            te,
            self.controls.maximum_target_thermal_speed_ratio,
            self.controls.rate_controls,
            cancelled,
        )?;
        if !g.is_finite() || g < 0.0 || !(g + coulomb).is_finite() {
            return Err(ElectronHeavyEnergyError::Arithmetic(
                "elastic electron-heavy coupling",
            ));
        }
        Ok(Coupling {
            neutral: g,
            coulomb,
        })
    }

    pub fn step(
        &self,
        start: &Point,
        end: f64,
        coupling: Coupling,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Point, ElectronHeavyEnergyError> {
        let total_coupling = coupling.total();
        let alpha = (end - start.time) * total_coupling;
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
                    coulomb_exchange: start.coulomb_exchange
                        + if total_coupling == 0.0 {
                            0.0
                        } else {
                            increment * (coupling.coulomb / total_coupling)
                        },
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
