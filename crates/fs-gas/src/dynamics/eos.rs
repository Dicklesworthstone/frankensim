use super::{Conserved, Error, Primitive, poll, positive};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::mixture::{FrozenIdealGasMixtureEvaluationV1, FrozenIdealGasMixtureModelV1};

/// Explicit source-temperature path and bounded caloric inversion controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricDomain {
    /// Lower supported temperature [K], also the sensible-energy reference.
    pub minimum_temperature: f64,
    /// Upper supported temperature [K].
    pub maximum_temperature: f64,
    /// Absolute specific-energy inversion tolerance [J/kg].
    pub absolute_energy_tolerance: f64,
    /// Relative specific-energy inversion tolerance.
    pub relative_energy_tolerance: f64,
    /// Maximum Newton/bisection iterations, in 1..=128.
    pub maximum_iterations: usize,
}

impl CaloricDomain {
    /// Bounded, tight default inversion for a caller-declared temperature path.
    pub fn new(minimum_temperature: Temperature, maximum_temperature: Temperature) -> Self {
        Self {
            minimum_temperature: minimum_temperature.value(),
            maximum_temperature: maximum_temperature.value(),
            absolute_energy_tolerance: 1e-10,
            relative_energy_tolerance: 1e-12,
            maximum_iterations: 64,
        }
    }
}

/// Immutable, source-owning frozen-mixture equation of state for compressible flow.
///
/// All cells share this exact composition. The source energy at Tmin is retained
/// separately, while the finite-volume state evolves its reference-invariant
/// sensible energy. No constant-gamma closure replaces the NASA-9 caloric law.
#[derive(Clone, Debug, PartialEq)]
pub struct FrozenGas {
    source: FrozenIdealGasMixtureModelV1,
    domain: CaloricDomain,
    gas_constant: f64,
    energy_reference: f64,
    maximum_energy: f64,
    reference_pressure: f64,
}

impl FrozenGas {
    /// Admit a continuous source path and positive component Cv over that path.
    ///
    /// Crossing source regions requires the source owner's explicit caloric
    /// reconciliation. A bounded interval subdivision of the original Cv/R
    /// polynomial refuses paths whose positive Cv cannot be established.
    pub fn new(
        source: FrozenIdealGasMixtureModelV1,
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
        {
            return Err(Error::Invalid("caloric domain or inversion controls"));
        }
        for component in source.components() {
            poll(&mut cancelled)?;
            let regions = component.regions();
            let selected = |t| {
                regions.iter().rposition(|r| {
                    r.temperature_min().value() <= t && t <= r.temperature_max().value()
                })
            };
            let first = selected(domain.minimum_temperature)
                .ok_or(Error::CaloricPath("lower endpoint unsupported"))?;
            let last = selected(domain.maximum_temperature)
                .ok_or(Error::CaloricPath("upper endpoint unsupported"))?;
            if first != last && component.caloric_reference().is_none() {
                return Err(Error::CaloricPath(
                    "cross-region flow requires explicit source caloric continuity",
                ));
            }
            for i in first..=last {
                if i > first
                    && regions[i - 1].temperature_max().value()
                        != regions[i].temperature_min().value()
                {
                    return Err(Error::CaloricPath("source temperature gap"));
                }
                let lo = regions[i]
                    .temperature_min()
                    .value()
                    .max(domain.minimum_temperature);
                let hi = regions[i]
                    .temperature_max()
                    .value()
                    .min(domain.maximum_temperature);
                positive_cv(*regions[i].coefficients(), lo, hi, &mut cancelled)?;
            }
        }
        let reference_pressure = source.components()[0]
            .convention()
            .reference_pressure()
            .value();
        let low = source.evaluate_properties(
            Temperature::new(domain.minimum_temperature),
            Pressure::new(reference_pressure),
        )?;
        let energy_reference = low.u().value() / low.molar_mass().value();
        let gas_constant = low.mass_specific_gas_constant().value();
        let maximum_energy = source
            .specific_internal_energy_change(
                Temperature::new(domain.minimum_temperature),
                Temperature::new(domain.maximum_temperature),
            )?
            .value();
        if !energy_reference.is_finite() || !positive(gas_constant) || !positive(maximum_energy) {
            return Err(Error::CaloricPath(
                "unrepresentable source energy or gas constant",
            ));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            source,
            domain,
            gas_constant,
            energy_reference,
            maximum_energy,
            reference_pressure,
        })
    }
    /// Exact original source models and frozen composition.
    pub const fn source(&self) -> &FrozenIdealGasMixtureModelV1 {
        &self.source
    }
    /// Declared source path and inversion controls.
    pub const fn domain(&self) -> CaloricDomain {
        self.domain
    }
    /// Original source internal energy at the declared lower temperature [J/kg].
    pub const fn energy_reference(&self) -> f64 {
        self.energy_reference
    }
    /// Mixture-specific ideal-gas constant [J/(kg K)].
    pub const fn gas_constant(&self) -> f64 {
        self.gas_constant
    }

    /// Construct a conservative state from density [kg/m³], velocity [m/s], and pressure [Pa].
    pub fn conserved(
        &self,
        density: f64,
        velocity: f64,
        pressure: f64,
    ) -> Result<Conserved, Error> {
        Ok(self.state_at(density, velocity, pressure)?.0)
    }
    pub(super) fn state_at(
        &self,
        density: f64,
        velocity: f64,
        pressure: f64,
    ) -> Result<(Conserved, Primitive), Error> {
        let primitive = self.from_rho_u_p(density, velocity, pressure)?;
        let e = self.sensible(primitive.temperature)?;
        let q = Conserved {
            density,
            momentum: density * velocity,
            energy: density * (e + 0.5 * velocity * velocity),
        };
        if !q.finite() {
            return Err(Error::Invalid("unrepresentable conservative state"));
        }
        // The cached primitive and the stored conserved state must describe
        // the same caloric state. At high kinetic/sensible-energy ratios a
        // finite addition can erase the heat entirely, so finiteness alone
        // would admit a state that decodes to a different temperature.
        let recovered_velocity = q.momentum / q.density;
        let recovered_energy = q.energy / q.density - 0.5 * recovered_velocity * recovered_velocity;
        let tolerance =
            self.domain.absolute_energy_tolerance + self.domain.relative_energy_tolerance * e.abs();
        if !self.admissible_state(q) || (recovered_energy - e).abs() > tolerance {
            return Err(Error::Invalid(
                "conservative state cannot retain the declared caloric precision",
            ));
        }
        Ok((q, primitive))
    }

    /// Resolve one conservative cell through the original caloric law.
    pub fn primitive(
        &self,
        state: Conserved,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        self.decode(state, None, &mut cancelled)
    }

    /// Source-bearing thermodynamic evaluation of an accepted conservative cell.
    pub fn thermodynamics(
        &self,
        state: Conserved,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FrozenIdealGasMixtureEvaluationV1, Error> {
        let p = self.decode(state, None, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(self
            .source
            .evaluate(Temperature::new(p.temperature), Pressure::new(p.pressure))?)
    }

    /// Recover physical total energy density, including the original signed formation reference.
    pub fn absolute_energy_density(&self, state: Conserved) -> Result<f64, Error> {
        let energy = state.energy + state.density * self.energy_reference;
        if !energy.is_finite() {
            return Err(Error::Invalid("unrepresentable absolute energy"));
        }
        Ok(energy)
    }

    pub(super) fn sensible(&self, t: f64) -> Result<f64, Error> {
        Ok(self
            .source
            .specific_internal_energy_change(
                Temperature::new(self.domain.minimum_temperature),
                Temperature::new(t),
            )?
            .value())
    }
    // Continuous positive Cv makes this a complete thermodynamic interval
    // membership check. Riemann trial states need no temperature iteration.
    pub(super) fn admissible_state(&self, state: Conserved) -> bool {
        if !state.finite() || !positive(state.density) {
            return false;
        }
        let velocity = state.momentum / state.density;
        let energy = state.energy / state.density - 0.5 * velocity * velocity;
        energy.is_finite() && energy >= 0.0 && energy <= self.maximum_energy
    }
    pub(super) fn from_rho_u_p(
        &self,
        density: f64,
        velocity: f64,
        pressure: f64,
    ) -> Result<Primitive, Error> {
        if !positive(density) || !velocity.is_finite() || !positive(pressure) {
            return Err(Error::Invalid("primitive density, velocity or pressure"));
        }
        let t = (pressure / density) / self.gas_constant;
        if !(self.domain.minimum_temperature..=self.domain.maximum_temperature).contains(&t) {
            return Err(Error::CaloricPath(
                "primitive temperature outside declared interval",
            ));
        }
        self.at_temperature(density, velocity, t)
    }
    fn at_temperature(&self, density: f64, velocity: f64, t: f64) -> Result<Primitive, Error> {
        let pressure = density * (self.gas_constant * t);
        if !positive(pressure) {
            return Err(Error::Invalid("unrepresentable EOS pressure"));
        }
        let p = self
            .source
            .evaluate_properties(Temperature::new(t), Pressure::new(pressure))?;
        let gamma = p.mass_specific_cp().value() / p.mass_specific_cv().value();
        let a2 = gamma * self.gas_constant * t;
        if !positive(p.mass_specific_cv().value()) || !positive(a2) {
            return Err(Error::CaloricPath("nonpositive Cv or frozen sound speed"));
        }
        Ok(Primitive {
            density,
            velocity,
            pressure,
            temperature: t,
            sound_speed: fs_math::det::sqrt(a2),
        })
    }
    pub(super) fn decode(
        &self,
        state: Conserved,
        guess: Option<f64>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        poll(cancelled)?;
        if !state.finite() || !positive(state.density) {
            return Err(Error::Invalid("conservative density or non-finite state"));
        }
        let velocity = state.momentum / state.density;
        let energy = state.energy / state.density - 0.5 * velocity * velocity;
        if !energy.is_finite() || energy < 0.0 || energy > self.maximum_energy {
            return Err(Error::EnergyDomain {
                requested: energy,
                maximum: self.maximum_energy,
            });
        }
        if energy == 0.0 {
            return self.at_temperature(state.density, velocity, self.domain.minimum_temperature);
        }
        if energy == self.maximum_energy {
            return self.at_temperature(state.density, velocity, self.domain.maximum_temperature);
        }
        let (mut lo, mut hi) = (
            self.domain.minimum_temperature,
            self.domain.maximum_temperature,
        );
        let mut t = guess
            .filter(|t| *t > lo && *t < hi)
            .unwrap_or(lo + (hi - lo) * (energy / self.maximum_energy));
        let tolerance = self.domain.absolute_energy_tolerance
            + self.domain.relative_energy_tolerance * energy.abs();
        for _ in 0..self.domain.maximum_iterations {
            poll(cancelled)?;
            let residual = self.sensible(t)? - energy;
            if residual.abs() <= tolerance {
                return self.at_temperature(state.density, velocity, t);
            }
            if residual < 0.0 {
                lo = t;
            } else {
                hi = t;
            }
            let p = self
                .source
                .evaluate_properties(Temperature::new(t), Pressure::new(self.reference_pressure))?;
            let cv = p.mass_specific_cv().value();
            if !positive(cv) {
                return Err(Error::CaloricPath("nonpositive trial Cv"));
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

// Each term is monotone on positive T. Directed rounding bounds each original
// NASA-9 Cv/R monomial; subdivision resolves cancellation without sampling a
// polynomial at a few points and mistaking that for positive-path admission.
fn cv_lower_bound(a: [f64; 9], lo: f64, hi: f64) -> f64 {
    let mut low = (a[2] - 1.0).next_down();
    for (coefficient, lower, upper) in [
        (a[0], (1.0 / hi).next_down(), (1.0 / lo).next_up()),
        (a[1], (1.0 / hi).next_down(), (1.0 / lo).next_up()),
        (a[3], lo, hi),
        (a[4], lo, hi),
        (a[5], lo, hi),
        (a[6], lo, hi),
    ]
    .into_iter()
    .enumerate()
    .map(|(i, (c, l, h))| {
        let power = match i {
            0 => 2,
            1 | 2 => 1,
            3 => 2,
            4 => 3,
            _ => 4,
        };
        let (mut lower, mut upper) = (1.0_f64, 1.0_f64);
        for _ in 0..power {
            lower = (lower * l).next_down();
            upper = (upper * h).next_up();
        }
        (c, lower, upper)
    }) {
        if coefficient == 0.0 {
            continue;
        }
        let term = coefficient * if coefficient > 0.0 { lower } else { upper };
        low = (low + term.next_down()).next_down();
    }
    low
}
fn positive_cv(
    a: [f64; 9],
    lo: f64,
    hi: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    // Depth-first fixed-size stack: at most 13 active pieces, 4095 work items.
    let mut stack = [(0.0, 0.0, 0_u8); 14];
    stack[0] = (lo, hi, 0);
    let mut n = 1;
    let mut work = 0;
    while n > 0 {
        poll(cancelled)?;
        n -= 1;
        let (l, h, depth) = stack[n];
        work += 1;
        let bound = cv_lower_bound(a, l, h);
        if bound.is_finite() && bound > 0.0 {
            continue;
        }
        let mid = l + 0.5 * (h - l);
        if depth == 12 || work >= 4095 || mid <= l || mid >= h {
            return Err(Error::CaloricPath(
                "positive component Cv could not be bounded on the source path",
            ));
        }
        stack[n] = (mid, h, depth + 1);
        stack[n + 1] = (l, mid, depth + 1);
        n += 2;
    }
    Ok(())
}

#[cfg(test)]
#[path = "eos_regression_tests.rs"]
mod eos_regression_tests;
