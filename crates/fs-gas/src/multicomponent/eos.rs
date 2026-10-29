use super::{
    CaloricDomain, Cell, Conserved, Error, MAX_COMPONENTS, MAX_SPECIES, Primitive, allocated, cell,
    poll, positive, sum,
};
use fs_qty::Temperature;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::kinetics::{KineticsErrorV1, admit_caloric_interval};

/// Immutable complete species bank for a composition-dependent caloric EOS.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureGas {
    source: IdealGasEquilibriumV1,
    domain: CaloricDomain,
    references: Vec<f64>,
    maximum_energies: Vec<f64>,
    inverse_masses: Vec<f64>,
}
impl MixtureGas {
    /// Admit every original species over the complete declared temperature path.
    /// No equilibrium calculation or frozen-composition approximation is made.
    pub fn new(
        source: IdealGasEquilibriumV1,
        domain: CaloricDomain,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if !positive(domain.minimum_temperature)
            || !positive(domain.maximum_temperature)
            || domain.minimum_temperature >= domain.maximum_temperature
            || !positive(domain.absolute_energy_tolerance)
            || !positive(domain.relative_energy_tolerance)
            || domain.relative_energy_tolerance > 1e-3
            || !(1..=128).contains(&domain.maximum_iterations)
            || source.models().is_empty()
            || source.models().len() > MAX_SPECIES
        {
            return Err(Error::Invalid("caloric domain, controls or species count"));
        }
        let n = source.models().len();
        let mut references = allocated(n, 0.0)?;
        let mut maximum_energies = allocated(n, 0.0)?;
        let mut inverse_masses = allocated(n, 0.0)?;
        let lo = Temperature::new(domain.minimum_temperature);
        let hi = Temperature::new(domain.maximum_temperature);
        for (i, model) in source.models().iter().enumerate() {
            admit_caloric_interval(model, lo, hi, &mut cancelled).map_err(|e| match e {
                KineticsErrorV1::Cancelled => Error::Cancelled,
                other => Error::CaloricPath(other),
            })?;
            inverse_masses[i] = 1.0 / model.molar_mass().value();
            references[i] = model.evaluate_properties(lo)?.u().value() * inverse_masses[i];
            maximum_energies[i] =
                model.molar_internal_energy_change(lo, hi)?.value() * inverse_masses[i];
            if !positive(inverse_masses[i])
                || !positive(maximum_energies[i])
                || !references[i].is_finite()
            {
                return Err(Error::Invalid(
                    "unrepresentable species energy or molar mass",
                ));
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            source,
            domain,
            references,
            maximum_energies,
            inverse_masses,
        })
    }
    /// Original, unmodified elemental/charge/NASA-9 source owner.
    pub const fn source(&self) -> &IdealGasEquilibriumV1 {
        &self.source
    }
    /// Explicit source-temperature path and inversion budget.
    pub const fn domain(&self) -> CaloricDomain {
        self.domain
    }
    /// Number of canonical species carried by every cell.
    pub fn species_count(&self) -> usize {
        self.references.len()
    }
    /// Original specific internal energy at Tmin for every species [J/kg].
    pub fn energy_references(&self) -> &[f64] {
        &self.references
    }
    pub(super) fn inverse_molar_masses(&self) -> &[f64] {
        &self.inverse_masses
    }
    /// Construct a conservative state from partial mass densities, velocity and pressure.
    pub fn conserved(
        &self,
        partial_densities: &[f64],
        velocity: f64,
        pressure: f64,
    ) -> Result<Conserved, Error> {
        let n = self.species_count();
        let mut values = [0.0; MAX_COMPONENTS];
        self.state_at(partial_densities, velocity, pressure, &mut values[..n + 2])?;
        let mut densities = allocated(n, 0.0)?;
        densities.copy_from_slice(partial_densities);
        Ok(Conserved {
            partial_densities: densities,
            momentum: values[n],
            energy: values[n + 1],
        })
    }
    /// Decode one conservative cell using the current composition's actual Cv.
    pub fn primitive(
        &self,
        state: Cell<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        self.decode_cell(state, None, &mut cancelled)
    }
    /// Physical total energy, including the retained formation references [J/m³].
    pub fn absolute_energy_density(&self, state: Cell<'_>) -> Result<f64, Error> {
        self.admit_densities(state.partial_densities)?;
        let energy = state.energy
            + sum(state
                .partial_densities
                .iter()
                .zip(&self.references)
                .map(|(a, b)| a * b));
        if !energy.is_finite() {
            return Err(Error::Invalid("unrepresentable absolute energy"));
        }
        Ok(energy)
    }
    pub(super) fn admit_densities(&self, densities: &[f64]) -> Result<(f64, f64), Error> {
        if densities.len() != self.species_count()
            || densities.iter().any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(Error::Invalid(
                "species axis or negative/non-finite partial density",
            ));
        }
        let rho = sum(densities.iter().copied());
        let concentration = sum(densities
            .iter()
            .zip(&self.inverse_masses)
            .map(|(a, b)| a * b));
        if !positive(rho) || !positive(concentration) {
            return Err(Error::Invalid("zero or unrepresentable mixture density"));
        }
        Ok((rho, concentration))
    }
    pub(super) fn state_at(
        &self,
        densities: &[f64],
        velocity: f64,
        pressure: f64,
        out: &mut [f64],
    ) -> Result<Primitive, Error> {
        let (rho, concentration) = self.admit_densities(densities)?;
        if !velocity.is_finite() || !positive(pressure) || out.len() != self.species_count() + 2 {
            return Err(Error::Invalid("primitive velocity, pressure or workspace"));
        }
        let t = (pressure / concentration) / R;
        let (thermal, cv) = self.caloric(densities, t)?;
        let kinetic = 0.5 * rho * velocity * velocity;
        let n = self.species_count();
        out[..n].copy_from_slice(densities);
        out[n] = rho * velocity;
        out[n + 1] = thermal + kinetic;
        let recovered = out[n + 1] - 0.5 * out[n] * (out[n] / rho);
        let tolerance = rho * self.domain.absolute_energy_tolerance
            + self.domain.relative_energy_tolerance * thermal.abs();
        if !self.admissible(out) || (recovered - thermal).abs() > tolerance {
            return Err(Error::Invalid(
                "conservative state loses declared caloric precision",
            ));
        }
        self.primitive_at(rho, concentration, velocity, t, cv)
    }
    // Returns sensible internal-energy density and Cv density. Formation
    // constants cancel within each original species owner before rounding.
    fn caloric(&self, densities: &[f64], t: f64) -> Result<(f64, f64), Error> {
        if !(self.domain.minimum_temperature..=self.domain.maximum_temperature).contains(&t) {
            return Err(Error::EnergyDomain);
        }
        let (mut energy, mut energy_correction) = (0.0, 0.0);
        let (mut cv, mut cv_correction) = (0.0, 0.0);
        let temperature = Temperature::new(t);
        let reference = Temperature::new(self.domain.minimum_temperature);
        for (i, model) in self.source.models().iter().enumerate() {
            if densities[i] == 0.0 {
                continue;
            }
            let c = densities[i] * self.inverse_masses[i];
            let energy_term = c * model
                .molar_internal_energy_change(reference, temperature)?
                .value();
            let cv_term = c * (model.evaluate_properties(temperature)?.cp().value() - R);
            let y = energy_term - energy_correction;
            let next = energy + y;
            energy_correction = (next - energy) - y;
            energy = next;
            let y = cv_term - cv_correction;
            let next = cv + y;
            cv_correction = (next - cv) - y;
            cv = next;
        }
        if !energy.is_finite() || !positive(cv) {
            return Err(Error::Invalid(
                "non-finite caloric energy or nonpositive mixture Cv",
            ));
        }
        Ok((energy, cv))
    }
    fn primitive_at(
        &self,
        rho: f64,
        concentration: f64,
        velocity: f64,
        t: f64,
        cv: f64,
    ) -> Result<Primitive, Error> {
        let pressure = concentration * (R * t);
        let a2 = (1.0 + concentration * R / cv) * (pressure / rho);
        if !positive(pressure) || !positive(a2) {
            return Err(Error::Invalid("unrepresentable pressure or sound speed"));
        }
        Ok(Primitive {
            density: rho,
            velocity,
            pressure,
            temperature: t,
            sound_speed: fs_math::det::sqrt(a2),
        })
    }
    pub(super) fn admissible(&self, values: &[f64]) -> bool {
        let n = self.species_count();
        values.len() == n + 2 && self.admissible_cell(cell(values, n))
    }
    pub(super) fn admissible_cell(&self, state: Cell<'_>) -> bool {
        if !state.momentum.is_finite() || !state.energy.is_finite() {
            return false;
        }
        let Ok((rho, _)) = self.admit_densities(state.partial_densities) else {
            return false;
        };
        let thermal = state.energy - 0.5 * state.momentum * (state.momentum / rho);
        let upper = sum(state
            .partial_densities
            .iter()
            .zip(&self.maximum_energies)
            .map(|(a, b)| a * b));
        thermal.is_finite() && upper.is_finite() && thermal >= 0.0 && thermal <= upper
    }
    pub(super) fn decode(
        &self,
        values: &[f64],
        guess: Option<f64>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        if values.len() != self.species_count() + 2 {
            return Err(Error::Invalid("conservative component count"));
        }
        self.decode_cell(cell(values, self.species_count()), guess, cancelled)
    }
    pub(super) fn decode_cell(
        &self,
        state: Cell<'_>,
        guess: Option<f64>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        poll(cancelled)?;
        let (rho, c) = self.admit_densities(state.partial_densities)?;
        let velocity = state.momentum / rho;
        let thermal = state.energy - 0.5 * state.momentum * velocity;
        let upper = sum(state
            .partial_densities
            .iter()
            .zip(&self.maximum_energies)
            .map(|(a, b)| a * b));
        if !velocity.is_finite()
            || !thermal.is_finite()
            || !upper.is_finite()
            || thermal < 0.0
            || thermal > upper
        {
            return Err(Error::EnergyDomain);
        }
        let (mut lo, mut hi) = (
            self.domain.minimum_temperature,
            self.domain.maximum_temperature,
        );
        let mut t = if thermal == 0.0 {
            lo
        } else if thermal == upper {
            hi
        } else {
            guess
                .filter(|x| *x > lo && *x < hi)
                .unwrap_or(lo + (hi - lo) * (thermal / upper))
        };
        let tolerance = rho * self.domain.absolute_energy_tolerance
            + self.domain.relative_energy_tolerance * thermal.abs();
        for _ in 0..self.domain.maximum_iterations {
            poll(cancelled)?;
            let (energy, cv) = self.caloric(state.partial_densities, t)?;
            let residual = energy - thermal;
            if residual.abs() <= tolerance {
                return self.primitive_at(rho, c, velocity, t, cv);
            }
            if residual < 0.0 {
                lo = t;
            } else {
                hi = t;
            }
            let trial = t - residual / cv;
            t = if trial > lo && trial < hi {
                trial
            } else {
                lo + 0.5 * (hi - lo)
            };
        }
        Err(Error::TemperatureIterations)
    }
}
