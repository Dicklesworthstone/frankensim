//! Prepared reacting isentropes. This module owns thermodynamics, not flow areas,
//! vehicle forces, or a transport law. Queries reuse bounded numerical storage.

use super::*;

/// Allocation-free scalar result of one converged equilibrium isentrope query.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumIsentropeSampleV1 {
    /// Actual absolute temperature [K].
    pub temperature: Temperature,
    /// Actual absolute pressure [Pa].
    pub pressure: Pressure,
    /// Reacting response at this same composition, T and p.
    pub response: EquilibriumThermodynamicResponseV1,
    /// H-H_reference divided by the original conserved feed mass [J/kg].
    /// Includes composition-dependent formation energy exactly once.
    pub specific_enthalpy_change_j_kg: f64,
    /// Actual S-S_reference [J/K], evaluated before subtracting large totals.
    pub entropy_change_j_k: f64,
    /// Actual product-minus-feed mass divided by original feed mass.
    pub relative_mass_residual: f64,
    /// Composition iterations in this query, including its temperature search.
    pub composition_iterations: usize,
    /// Temperature iterations in this query.
    pub temperature_iterations: usize,
}

struct Reference {
    temperature: Temperature,
    pressure: Pressure,
    amounts: Vec<f64>,
    log_amounts: Vec<f64>,
    enthalpies: Vec<f64>,
    partial_entropies: Vec<f64>,
    total_amount: f64,
}
impl Reference {
    fn differences<const WITH_ENTHALPY: bool>(&self, owner: &IdealGasEquilibriumV1, work: &Workspace,
        temperature: Temperature, cancelled: &mut impl FnMut() -> bool)
        -> Result<(f64, f64), EquilibriumExpansionErrorV1> {
        let mut total = Sum::default();
        for &i in &work.active { total.add(work.n[i]); }
        let log_total_change = crate::entropy::log_ratio(total.value(), self.total_amount);
        let log_pressure_change = crate::entropy::log_ratio(work.pressure.value(), self.pressure.value());
        let mut dh = Sum::default();
        let mut ds = Sum::default();
        for &i in &work.active {
            poll(cancelled)?;
            let n = work.n[i];
            let dn = n - self.amounts[i];
            let model = &owner.models[i];
            if WITH_ENTHALPY {
                let h = model.molar_enthalpy_change(self.temperature, temperature)
                    .map_err(|source| owner.species_error(i, source))?.value();
                dh.add(dn * self.enthalpies[i]);
                dh.add(n * h);
            }
            let s = model.molar_entropy_change(self.temperature, temperature)
                .map_err(|source| owner.species_error(i, source))?.value();
            let log_activity_change = (work.log_n[i] - self.log_amounts[i])
                - log_total_change + log_pressure_change;
            ds.add(dn * self.partial_entropies[i]);
            ds.add(n * (s - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * log_activity_change));
        }
        // Composition changes are real: unlike a frozen-temperature difference,
        // their formation and partial-entropy reference terms cannot be omitted.
        let dh = dh.value() / work.initial_mass;
        let ds = ds.value() * work.scale;
        if !dh.is_finite() || !ds.is_finite() {
            return Err(E::Arithmetic("relative reacting enthalpy or entropy").into());
        }
        Ok((dh, ds))
    }
}

/// One immutable source/reference state and a reusable mutable numerical workspace.
///
/// `sample` changes only warm-start guesses, never the source, element inventory,
/// reference state or a physical parcel. A failed query invalidates publication
/// of the current trial. Repeating the same query sequence is deterministic on
/// one arithmetic target; arbitrary query-order bit identity is not asserted.
/// Workspace allocation is bounded by the original 128-species/16-element caps.
/// The original full-KKT fallback may allocate once when first needed.
pub struct EquilibriumIsentropeV1<'a> {
    owner: &'a IdealGasEquilibriumV1,
    work: Workspace,
    reference: Reference,
    reference_state: EquilibriumThermodynamicStateV1,
    controls: EquilibriumEntropySolveV1,
    current: Option<EquilibriumIsentropeSampleV1>,
}
impl<'a> EquilibriumIsentropeV1<'a> {
    /// Establish an equilibrium reference at explicit T,p and feed elements.
    /// All eligible species must support the complete declared temperature path.
    pub fn new(owner: &'a IdealGasEquilibriumV1, temperature: Temperature,
        pressure: Pressure, feed: &[(SpeciesId, Amount)], controls: EquilibriumEntropySolveV1,
        mut cancelled: impl FnMut() -> bool) -> Result<Self, EquilibriumExpansionErrorV1> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !(controls.minimum_temperature_k..=controls.maximum_temperature_k).contains(&temperature.value()) {
            return Err(EquilibriumExpansionErrorV1::Invalid("isentrope reference is outside the temperature bracket"));
        }
        let mut work = Workspace::new(owner, pressure, feed)?;
        work.admit_temperature_path(owner, controls)?;
        let iterations = work.solve(owner, temperature, controls.equilibrium, &mut cancelled)?;
        let reference_state = work.publish_response(owner, temperature, controls.equilibrium,
            iterations, 0, None, &mut cancelled)?;
        let feed_mass = work.initial_mass * work.scale;
        if !feed_mass.is_finite() || feed_mass <= 0.0 {
            return Err(E::Arithmetic("original feed mass").into());
        }
        let ns = work.n.len();
        let mut reference = Reference {
            temperature, pressure, amounts: allocated(ns, 0.0)?, log_amounts: allocated(ns, 0.0)?,
            enthalpies: allocated(ns, 0.0)?, partial_entropies: allocated(ns, 0.0)?, total_amount: 0.0,
        };
        reference.amounts.copy_from_slice(&work.n);
        reference.log_amounts.copy_from_slice(&work.log_n);
        reference.enthalpies.copy_from_slice(&work.h);
        let mut total = Sum::default();
        for &i in &work.active { total.add(work.n[i]); }
        reference.total_amount = total.value();
        let log_total = ln(reference.total_amount);
        for &i in &work.active {
            poll(&mut cancelled)?;
            let p = owner.models[i].evaluate_properties(temperature)
                .map_err(|source| owner.species_error(i, source))?;
            reference.partial_entropies[i] = p.s().value() - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                * (work.log_n[i] - log_total + work.log_pressure_ratio);
            if !reference.partial_entropies[i].is_finite() {
                return Err(E::Arithmetic("reference partial entropy").into());
            }
        }
        let current = Some(EquilibriumIsentropeSampleV1 {
            temperature, pressure, response: reference_state.response(),
            specific_enthalpy_change_j_kg: 0.0, entropy_change_j_k: 0.0,
            relative_mass_residual: reference_state.equilibrium().relative_mass_residual(),
            composition_iterations: iterations, temperature_iterations: 0,
        });
        poll(&mut cancelled)?;
        Ok(Self { owner, work, reference, reference_state, controls, current })
    }
    /// Original complete source-bearing reference state.
    pub const fn reference_state(&self) -> &EquilibriumThermodynamicStateV1 { &self.reference_state }
    /// Original conserved feed mass [kg], independent of numerical product residuals.
    pub fn feed_mass_kg(&self) -> f64 { self.work.initial_mass * self.work.scale }

    /// Re-equilibrate at one positive pressure while preserving reference entropy.
    /// No receipt is allocated for temperature/root trial states. Relative source
    /// differences retain small sensible changes; composition-reference terms and
    /// their floating-point conservation sensitivity remain explicit.
    pub fn sample(&mut self, pressure: Pressure, mut cancelled: impl FnMut() -> bool)
        -> Result<EquilibriumIsentropeSampleV1, EquilibriumExpansionErrorV1> {
        let previous = self.current.take();
        poll(&mut cancelled)?;
        if !pressure.value().is_finite() || pressure.value() <= 0.0 {
            return Err(EquilibriumExpansionErrorV1::Invalid("isentrope pressure must be finite and positive"));
        }
        self.work.pressure = pressure;
        self.work.log_pressure_ratio = ln(pressure.value()) - ln(self.owner.convention.reference_pressure().value());
        let preferred = previous.map(|s| s.temperature.value() * fs_math::det::exp(
            s.response.isentropic_temperature_exponent
                * crate::entropy::log_ratio(pressure.value(), s.pressure.value())));
        let reference = &self.reference;
        let owner = self.owner;
        let (temperature, iterations, thermal_iterations) = self.work.solve_entropy(
            owner, 0.0, self.controls, preferred, &mut cancelled,
            |work, t, cancel| Ok(reference.differences::<false>(owner, work, t, cancel)?.1),
        )?;
        let response = self.work.response(owner, temperature, self.controls.equilibrium, &mut cancelled)?;
        let mut mass = Sum::default();
        for &i in &self.work.active { mass.add(self.work.n[i] * owner.models[i].molar_mass().value()); }
        let relative_mass_residual = (mass.value() - self.work.initial_mass) / self.work.initial_mass;
        if !relative_mass_residual.is_finite()
            || relative_mass_residual.abs() > self.controls.equilibrium.relative_balance_tolerance {
            return Err(E::MassBalance { relative_residual: relative_mass_residual }.into());
        }
        let (dh, ds) = reference.differences::<true>(owner, &self.work, temperature, &mut cancelled)?;
        if ds.abs() > self.controls.entropy_tolerance_j_k {
            return Err(EquilibriumExpansionErrorV1::EntropyNoConvergence {
                iterations: thermal_iterations, residual_j_k: ds,
            });
        }
        let sample = EquilibriumIsentropeSampleV1 { temperature, pressure, response,
            specific_enthalpy_change_j_kg: dh, entropy_change_j_k: ds, relative_mass_residual,
            composition_iterations: iterations, temperature_iterations: thermal_iterations };
        poll(&mut cancelled)?;
        self.current = Some(sample);
        Ok(sample)
    }

    /// Publish the last successful query with original per-species source receipts.
    /// A failed or cancelled query cannot leave a publishable stale endpoint.
    /// The scalar sample's direct entropy difference is the accurate SP residual;
    /// this endpoint's absolute entropy remains on the original source reference.
    pub fn publish_current(&mut self, mut cancelled: impl FnMut() -> bool)
        -> Result<EquilibriumThermodynamicStateV1, EquilibriumExpansionErrorV1> {
        poll(&mut cancelled)?;
        let sample = self.current.ok_or(EquilibriumExpansionErrorV1::Invalid("no accepted isentrope query to publish"))?;
        let state = self.work.publish_response(self.owner, sample.temperature, self.controls.equilibrium,
            sample.composition_iterations, sample.temperature_iterations,
            Some(self.reference_state.entropy_j_k()), &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(state)
    }
}

/// Density inversion on this same source-bound reacting isentrope.
pub mod density {
use super::*;

/// Positive pressure interval and density/root tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumDensitySolveV1 {
    /// Inclusive lower pressure [Pa].
    pub minimum_pressure_pa: f64,
    /// Inclusive upper pressure [Pa], greater than the lower pressure.
    pub maximum_pressure_pa: f64,
    /// Allowed absolute rho/target-1, in (0,0.01).
    pub relative_density_tolerance: f64,
    /// At most 256 pressure queries; each retains the original SP work bounds.
    pub maximum_pressure_iterations: usize,
}
/// A failed density query cannot leave a publishable stale endpoint.
#[derive(Clone, Debug, PartialEq)]
pub enum EquilibriumDensityErrorV1 {
    /// Original source, composition, entropy or cancellation refusal.
    Thermochemistry(EquilibriumExpansionErrorV1),
    /// Invalid target, pressure interval or work controls.
    Invalid(&'static str),
    /// The pressure interval/work budget did not close density.
    NoConvergence { pressure_queries: usize, relative_density_residual: f64 },
}
impl From<EquilibriumExpansionErrorV1> for EquilibriumDensityErrorV1 {
    fn from(value: EquilibriumExpansionErrorV1) -> Self { Self::Thermochemistry(value) }
}
impl core::fmt::Display for EquilibriumDensityErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "equilibrium isentrope density: {self:?}")
    }
}
impl std::error::Error for EquilibriumDensityErrorV1 {}
/// Accepted isentrope sample and actual accumulated work, not certified error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumDensityQueryV1 {
    /// Same source-relative enthalpy, entropy and reacting response as sample(p).
    pub sample: EquilibriumIsentropeSampleV1,
    /// Actual (rho-target)/target.
    pub relative_density_residual: f64,
    /// Number of pressure queries, zero for an already closed retained endpoint.
    pub pressure_queries: usize,
    /// Sum over all queries, not only the accepted sample.
    pub composition_iterations: usize,
    /// Sum over all queries, not only the accepted sample.
    pub temperature_iterations: usize,
}
impl EquilibriumIsentropeV1<'_> {
    /// Last accepted scalar sample; None follows a failed/cancelled query.
    pub const fn current_sample(&self) -> Option<EquilibriumIsentropeSampleV1> { self.current }

    /// Solve rho(p,s,elements)=target using the original reacting response.
    /// Positive compressibility supplies dln(rho)/dln(p)=1/Gamma_s.
    /// Safeguarded log-pressure Newton uses this derivative and the previous
    /// endpoint as predictor. No new receipt or workspace is allocated. Every
    /// evaluated state retains original source/entropy admission. Interval ends
    /// are not evaluated just to initialize a good interior predictor. An
    /// out-of-interval target still cannot pass the actual density-residual gate.
    pub fn sample_density(&mut self, target: f64, controls: EquilibriumDensitySolveV1,
        mut cancelled: impl FnMut() -> bool) -> Result<EquilibriumDensityQueryV1, EquilibriumDensityErrorV1> {
        let result = self.density_candidate(target, controls, &mut cancelled);
        if result.is_err() { self.current = None; }
        result
    }
    fn density_candidate(&mut self, target: f64, c: EquilibriumDensitySolveV1,
        cancelled: &mut impl FnMut() -> bool) -> Result<EquilibriumDensityQueryV1, EquilibriumDensityErrorV1> {
        poll(cancelled).map_err(EquilibriumExpansionErrorV1::from)?;
        if !target.is_finite() || target <= 0.0
            || !c.minimum_pressure_pa.is_finite() || c.minimum_pressure_pa <= 0.0
            || !c.maximum_pressure_pa.is_finite() || c.maximum_pressure_pa <= c.minimum_pressure_pa
            || !c.relative_density_tolerance.is_finite() || c.relative_density_tolerance <= 0.0
            || c.relative_density_tolerance >= 0.01 || !(1..=256).contains(&c.maximum_pressure_iterations) {
            return Err(EquilibriumDensityErrorV1::Invalid("invalid density, pressure bounds or root controls"));
        }
        let (mut low, mut high) = (ln(c.minimum_pressure_pa), ln(c.maximum_pressure_pa));
        let mut x = low + 0.5 * (high - low);
        if let Some(previous) = self.current {
            let residual = (previous.response.density_kg_m3 - target) / target;
            if (c.minimum_pressure_pa..=c.maximum_pressure_pa).contains(&previous.pressure.value())
                && residual.abs() <= c.relative_density_tolerance {
                poll(cancelled).map_err(EquilibriumExpansionErrorV1::from)?;
                return Ok(EquilibriumDensityQueryV1 { sample: previous, relative_density_residual: residual,
                    pressure_queries: 0, composition_iterations: 0, temperature_iterations: 0 });
            }
            let prediction = ln(previous.pressure.value()) + previous.response.isentropic_pressure_exponent
                * crate::entropy::log_ratio(target, previous.response.density_kg_m3);
            if prediction.is_finite() && prediction >= low && prediction <= high { x = prediction; }
        }
        let mut composition_iterations = 0;
        let mut temperature_iterations = 0;
        let mut residual = f64::INFINITY;
        let mut queries = 0;
        for query in 1..=c.maximum_pressure_iterations {
            queries = query;
            let pressure = if x == low && low == ln(c.minimum_pressure_pa) { c.minimum_pressure_pa }
                else if x == high && high == ln(c.maximum_pressure_pa) { c.maximum_pressure_pa } else { exp(x) };
            if !pressure.is_finite() || pressure < c.minimum_pressure_pa || pressure > c.maximum_pressure_pa {
                return Err(EquilibriumDensityErrorV1::Invalid("unrepresentable pressure iterate"));
            }
            let sample = self.sample(Pressure::new(pressure), &mut *cancelled)?;
            composition_iterations += sample.composition_iterations;
            temperature_iterations += sample.temperature_iterations;
            residual = (sample.response.density_kg_m3 - target) / target;
            if residual.is_finite() && residual.abs() <= c.relative_density_tolerance {
                poll(cancelled).map_err(EquilibriumExpansionErrorV1::from)?;
                return Ok(EquilibriumDensityQueryV1 { sample, relative_density_residual: residual,
                    pressure_queries: query, composition_iterations, temperature_iterations });
            }
            if !residual.is_finite() { break; }
            if residual < 0.0 { low = x; } else { high = x; }
            let next = x - sample.response.isentropic_pressure_exponent
                * crate::entropy::log_ratio(sample.response.density_kg_m3, target);
            let midpoint = low + 0.5 * (high - low);
            let candidate = if next.is_finite() && next > low && next < high { next } else { midpoint };
            if candidate == x || candidate == low || candidate == high { break; }
            x = candidate;
        }
        Err(EquilibriumDensityErrorV1::NoConvergence {
            pressure_queries: queries, relative_density_residual: residual,
        })
    }
}
}
