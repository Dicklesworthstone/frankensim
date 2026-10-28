//! Neutral ideal-gas composition equilibrium from the original NASA-9 models.
//!
//! Minimize G/(RT)=sum n_i[g_i^0/(RT)+ln(y_i p/p0)] subject to A n=b.
//! Damped log-amount Newton iterations eliminate species increments into a
//! bounded element-space system. A full KKT fallback avoids squaring the
//! conditioning when trace species make that reduction numerically singular.
//! HP brackets temperature and uses the differentiated equilibrium equations
//! for safeguarded Newton steps. Formation enthalpy enters the source h once.
//!
//! Reference: Gordon and McBride, NASA RP-1311 Part I (1994), sections 2.3–2.4,
//! <https://ntrs.nasa.gov/citations/19950013764>. Neutral gas only: no condensed
//! phases, ionization, finite-rate chemistry, real-gas EOS or phase-stability
//! verdict. A product list is a caller declaration, not mechanism validation.

use crate::{
    ChargeVector, ElementalMatrix, Nasa9EvaluationReceiptV1, Nasa9StandardStateModelV1, SpeciesId,
    StandardStateConventionV1, ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
};
use core::fmt;
use fs_qty::{Amount, Energy, Mass, Pressure, Temperature};

mod solve;
/// Maximum retained gas species, including reactants and inerts.
pub const MAX_EQUILIBRIUM_SPECIES_V1: usize = 128;
/// Maximum elemental rows; dependent rows remain checked at publication.
pub const MAX_EQUILIBRIUM_ELEMENTS_V1: usize = 16;

/// Explicit nonlinear tolerances and work bounds. No concentration cutoff:
/// finite logarithmic trace abundances survive direct exponential underflow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumSolveV1 {
    /// Relative tolerance on each element separately, total moles, and mass.
    pub relative_balance_tolerance: f64,
    /// Maximum dimensionless chemical-potential stationarity residual.
    pub dimensionless_potential_tolerance: f64,
    /// At most 1024 Newton iterations per temperature.
    pub maximum_iterations: usize,
    /// At most 64 backtracking trials per iteration.
    pub maximum_backtracks: usize,
}
impl Default for EquilibriumSolveV1 {
    fn default() -> Self {
        Self {
            relative_balance_tolerance: 1e-10,
            dimensionless_potential_tolerance: 1e-9,
            maximum_iterations: 128,
            maximum_backtracks: 32,
        }
    }
}

/// Temperature bracket and total-energy tolerance for HP equilibrium.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumEnthalpySolveV1 {
    /// Each warm-started TP solve uses these controls.
    pub equilibrium: EquilibriumSolveV1,
    /// Source-supported lower temperature [K].
    pub minimum_temperature_k: f64,
    /// Source-supported upper temperature [K].
    pub maximum_temperature_k: f64,
    /// Absolute total enthalpy tolerance [J] for this feed amount.
    pub enthalpy_tolerance_j: f64,
    /// At most 256 temperature iterations, excluding bracket endpoints.
    pub maximum_temperature_iterations: usize,
}

/// Physical or numerical refusal before publication.
#[derive(Clone, Debug, PartialEq)]
pub enum CompositionEquilibriumErrorV1 {
    /// Invalid dimensions, identities, controls or input range.
    Invalid(&'static str),
    /// Supplied feed identity has no model.
    UnknownFeedSpecies(SpeciesId),
    /// Feed entries repeat an identity, including a zero entry.
    DuplicateFeedSpecies(SpeciesId),
    /// Models do not share the exact standard-state convention.
    ConventionMismatch(SpeciesId),
    /// Original source model refused the queried state.
    SpeciesEvaluation {
        species: SpeciesId,
        source: ThermochemErrorV1,
    },
    /// Bounded workspace reservation failed.
    AllocationRefused,
    /// Both element reduction and full KKT system are numerically singular.
    SingularSystem,
    /// An arithmetic result cannot be represented.
    Arithmetic(&'static str),
    /// Nonlinear work budget exhausted without closing all equations.
    NoConvergence {
        iterations: usize,
        residuals: EquilibriumResidualsV1,
    },
    /// Source-declared molecular masses disagree across the element inventory.
    MassBalance { relative_residual: f64 },
    /// The temperature bracket does not contain the requested enthalpy.
    EnthalpyNotBracketed {
        lower_enthalpy_j: f64,
        upper_enthalpy_j: f64,
    },
    /// Temperature convergence alone never replaces enthalpy closure.
    EnthalpyNoConvergence { iterations: usize, residual_j: f64 },
    /// Cancellation before publication.
    Cancelled,
}
impl fmt::Display for CompositionEquilibriumErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ideal-gas composition equilibrium: {self:?}")
    }
}
impl std::error::Error for CompositionEquilibriumErrorV1 {}

/// Actual residuals, not interval enclosures or data uncertainty estimates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumResidualsV1 {
    /// Maximum relative element discrepancy, scaled by that element's inventory.
    pub maximum_relative_element_residual: f64,
    /// Independently iterated total moles versus the species sum.
    pub relative_total_moles_residual: f64,
    /// Maximum normalized chemical-potential stationarity discrepancy.
    pub maximum_dimensionless_potential_residual: f64,
}

/// One product in the owning elemental matrix's canonical species order.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumSpeciesAmountV1 {
    /// Original species identity.
    pub species: SpeciesId,
    /// Amount [mol]; direct underflow is distinguished by the logarithm.
    pub amount: Amount,
    /// Mole fraction, with no caller-selected trace cutoff.
    pub mole_fraction: f64,
    /// None means an absent feed element structurally excludes the species.
    /// Some retains the finite logarithm even when the direct fraction is zero.
    pub log_mole_fraction: Option<f64>,
    /// Accepted source receipt; absent only for structurally excluded species.
    pub source: Option<Nasa9EvaluationReceiptV1>,
}

/// Complete equilibrium endpoint. Original source models remain immutable.
#[derive(Clone, Debug, PartialEq)]
pub struct IdealGasEquilibriumSolutionV1 {
    temperature: Temperature,
    pressure: Pressure,
    species: Vec<EquilibriumSpeciesAmountV1>,
    total_amount: Amount,
    mass: Mass,
    enthalpy: Energy,
    gibbs_energy: Energy,
    equilibrium_heat_capacity_p_j_k: Option<f64>,
    residuals: EquilibriumResidualsV1,
    relative_mass_residual: f64,
    composition_iterations: usize,
    temperature_iterations: usize,
    target_enthalpy: Option<Energy>,
}
impl IdealGasEquilibriumSolutionV1 {
    /// Product absolute temperature.
    #[must_use]
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    /// Product total pressure, used in the ideal-gas activities.
    #[must_use]
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    /// Canonical amounts and accepted sources.
    #[must_use]
    pub fn species(&self) -> &[EquilibriumSpeciesAmountV1] {
        &self.species
    }
    /// Total molecular amount, which may change during reaction.
    #[must_use]
    pub const fn total_amount(&self) -> Amount {
        self.total_amount
    }
    /// Mass from source-declared molecular masses.
    #[must_use]
    pub const fn mass(&self) -> Mass {
        self.mass
    }
    /// Total enthalpy, including source formation references once.
    #[must_use]
    pub const fn enthalpy(&self) -> Energy {
        self.enthalpy
    }
    /// Total Gibbs energy at the actual product T, p and composition.
    #[must_use]
    pub const fn gibbs_energy(&self) -> Energy {
        self.gibbs_energy
    }
    /// dH/dT at fixed pressure and element inventory, including composition
    /// response. None means the derivative system was singular/unrepresentable.
    /// This is not frozen cp and implies no reacting cp-cv=R identity.
    #[must_use]
    pub const fn equilibrium_heat_capacity_p_j_k(&self) -> Option<f64> {
        self.equilibrium_heat_capacity_p_j_k
    }
    /// Actual element, molecular-total and potential closure.
    #[must_use]
    pub const fn residuals(&self) -> EquilibriumResidualsV1 {
        self.residuals
    }
    /// Product-minus-feed mass divided by initial mass.
    #[must_use]
    pub const fn relative_mass_residual(&self) -> f64 {
        self.relative_mass_residual
    }
    /// Total Newton iterations, including HP temperature searches.
    #[must_use]
    pub const fn composition_iterations(&self) -> usize {
        self.composition_iterations
    }
    /// Temperature iterations; zero for TP or an accepted HP bracket endpoint.
    #[must_use]
    pub const fn temperature_iterations(&self) -> usize {
        self.temperature_iterations
    }
    /// Explicit HP target; absent for TP.
    #[must_use]
    pub const fn target_enthalpy(&self) -> Option<Energy> {
        self.target_enthalpy
    }
    /// Product enthalpy minus the explicit HP target.
    #[must_use]
    pub fn enthalpy_residual(&self) -> Option<Energy> {
        self.target_enthalpy
            .map(|h| Energy::new(self.enthalpy.value() - h.value()))
    }
}

/// Neutral product species sharing one exact ideal-gas standard-state convention.
#[derive(Clone, Debug, PartialEq)]
pub struct IdealGasEquilibriumV1 {
    elemental: ElementalMatrix,
    charges: ChargeVector,
    models: Vec<Nasa9StandardStateModelV1>,
    coefficients: Vec<f64>,
    convention: StandardStateConventionV1,
}
impl IdealGasEquilibriumV1 {
    /// Bind the original elemental, charge and NASA-9 owners. All charges must
    /// be zero; omitting a charge-conservation constraint is never implicit.
    pub fn new(
        elemental: ElementalMatrix,
        charges: ChargeVector,
        mut models: Vec<Nasa9StandardStateModelV1>,
    ) -> Result<Self, CompositionEquilibriumErrorV1> {
        use CompositionEquilibriumErrorV1 as E;
        let (m, n) = (elemental.row_count(), elemental.column_count());
        if m == 0 || m > MAX_EQUILIBRIUM_ELEMENTS_V1 || n == 0 || n > MAX_EQUILIBRIUM_SPECIES_V1 {
            return Err(E::Invalid("require 1..16 elements and 1..128 species"));
        }
        if charges.species() != elemental.species() || (0..n).any(|i| charges.get(i) != Some(0)) {
            return Err(E::Invalid(
                "charge species axis must match and every charge must be zero",
            ));
        }
        if models.len() != n {
            return Err(E::Invalid(
                "one model is required for every product species",
            ));
        }
        models.sort_by(|a, b| a.species().cmp(b.species()));
        if models
            .iter()
            .zip(elemental.species())
            .any(|(m, id)| m.species() != id)
        {
            return Err(E::Invalid(
                "model identities must exactly match elemental species",
            ));
        }
        let convention = models[0].convention().clone();
        for model in &models {
            if model.convention() != &convention {
                return Err(E::ConventionMismatch(model.species().clone()));
            }
        }
        let mut coefficients = allocated(m * n, 0.0)?;
        for i in 0..n {
            let mut any = false;
            for k in 0..m {
                let value = elemental
                    .get(k, i)
                    .ok_or(E::Invalid("missing elemental entry"))?;
                if !(0..=(1_i128 << 53)).contains(&value) {
                    return Err(E::Invalid(
                        "element counts must be exactly representable nonnegative integers",
                    ));
                }
                coefficients[k * n + i] = value as f64;
                any |= value > 0;
            }
            if !any {
                return Err(E::Invalid("every neutral species must contain an element"));
            }
        }
        Ok(Self {
            elemental,
            charges,
            models,
            coefficients,
            convention,
        })
    }
    /// Original complete elemental source, including dependent rows.
    #[must_use]
    pub const fn elemental_matrix(&self) -> &ElementalMatrix {
        &self.elemental
    }
    /// Original neutral-charge declaration.
    #[must_use]
    pub const fn charge_vector(&self) -> &ChargeVector {
        &self.charges
    }
    /// Canonical unmodified source models.
    #[must_use]
    pub fn models(&self) -> &[Nasa9StandardStateModelV1] {
        &self.models
    }
    /// Exact common phase/EOS/pressure/elemental-reference convention.
    #[must_use]
    pub const fn convention(&self) -> &StandardStateConventionV1 {
        &self.convention
    }

    /// Equilibrate at explicit T and p, conserving the supplied feed elements.
    /// Eligible products may appear without positive feed seeds. Interior
    /// iterates are numerical guesses and are never published physical states.
    pub fn equilibrate_tp(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        feed: &[(SpeciesId, Amount)],
        controls: EquilibriumSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IdealGasEquilibriumSolutionV1, CompositionEquilibriumErrorV1> {
        poll(&mut cancelled)?;
        controls.validate()?;
        let mut work = solve::Workspace::new(self, pressure, feed)?;
        let iterations = work.solve(self, temperature, controls, &mut cancelled)?;
        work.publish(
            self,
            temperature,
            controls,
            iterations,
            0,
            None,
            &mut cancelled,
        )
    }
    /// ORIGINAL feed enthalpy at one inlet temperature, without equilibration.
    /// Sum separate streams' values and call equilibrate_hp_enthalpy to model
    /// distinct oxidizer/fuel inlet temperatures in the same reference system.
    pub fn feed_enthalpy(
        &self,
        temperature: Temperature,
        feed: &[(SpeciesId, Amount)],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Energy, CompositionEquilibriumErrorV1> {
        poll(&mut cancelled)?;
        let initial = self.feed_amounts(feed)?;
        let mut value = Sum::default();
        for (i, amount) in initial.into_iter().enumerate() {
            if amount == 0.0 {
                continue;
            }
            poll(&mut cancelled)?;
            let p = self.models[i]
                .evaluate_properties(temperature)
                .map_err(|source| self.species_error(i, source))?;
            value.add(amount * p.h().value());
        }
        let value = value.value();
        if !value.is_finite() {
            return Err(CompositionEquilibriumErrorV1::Arithmetic("feed enthalpy"));
        }
        poll(&mut cancelled)?;
        Ok(Energy::new(value))
    }
    /// Adiabatic fixed-pressure products for a common reactant inlet temperature.
    /// Formation energy is already in h; no second heat-release term is added.
    pub fn equilibrate_hp(
        &self,
        pressure: Pressure,
        feed: &[(SpeciesId, Amount)],
        inlet_temperature: Temperature,
        controls: EquilibriumEnthalpySolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IdealGasEquilibriumSolutionV1, CompositionEquilibriumErrorV1> {
        let h = self.feed_enthalpy(inlet_temperature, feed, &mut cancelled)?;
        self.equilibrate_hp_enthalpy(pressure, feed, h, controls, cancelled)
    }
    /// Core fixed-pressure solve for an explicit total-enthalpy target [J].
    /// Distinct inlet temperatures or explicit heat/work contributions enter
    /// this one target, under the owner's common formation reference.
    pub fn equilibrate_hp_enthalpy(
        &self,
        pressure: Pressure,
        feed: &[(SpeciesId, Amount)],
        target_enthalpy: Energy,
        controls: EquilibriumEnthalpySolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IdealGasEquilibriumSolutionV1, CompositionEquilibriumErrorV1> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !target_enthalpy.value().is_finite() {
            return Err(CompositionEquilibriumErrorV1::Invalid(
                "target enthalpy must be finite",
            ));
        }
        let mut work = solve::Workspace::new(self, pressure, feed)?;
        work.solve_hp(self, target_enthalpy, controls, &mut cancelled)
    }
    fn feed_amounts(
        &self,
        feed: &[(SpeciesId, Amount)],
    ) -> Result<Vec<f64>, CompositionEquilibriumErrorV1> {
        use CompositionEquilibriumErrorV1 as E;
        let n = self.models.len();
        if feed.is_empty() || feed.len() > n {
            return Err(E::Invalid("feed requires 1..species_count entries"));
        }
        let mut amounts = allocated(n, 0.0)?;
        let mut seen = allocated(n, false)?;
        let mut total = Sum::default();
        for (id, amount) in feed {
            let i = self
                .elemental
                .species()
                .binary_search(id)
                .map_err(|_| E::UnknownFeedSpecies(id.clone()))?;
            if seen[i] {
                return Err(E::DuplicateFeedSpecies(id.clone()));
            }
            seen[i] = true;
            let x = amount.value();
            if !x.is_finite() || x < 0.0 {
                return Err(E::Invalid("feed amounts must be finite and nonnegative"));
            }
            amounts[i] = x;
            total.add(x);
        }
        if !total.value().is_finite() || total.value() <= 0.0 {
            return Err(E::Invalid("total feed amount must be positive and finite"));
        }
        Ok(amounts)
    }
    fn species_error(&self, i: usize, source: ThermochemErrorV1) -> CompositionEquilibriumErrorV1 {
        CompositionEquilibriumErrorV1::SpeciesEvaluation {
            species: self.models[i].species().clone(),
            source,
        }
    }
}
impl EquilibriumSolveV1 {
    fn validate(self) -> Result<(), CompositionEquilibriumErrorV1> {
        if !self.relative_balance_tolerance.is_finite()
            || self.relative_balance_tolerance <= 0.0
            || self.relative_balance_tolerance >= 0.01
            || !self.dimensionless_potential_tolerance.is_finite()
            || self.dimensionless_potential_tolerance <= 0.0
            || self.dimensionless_potential_tolerance >= 0.1
            || !(1..=1024).contains(&self.maximum_iterations)
            || !(1..=64).contains(&self.maximum_backtracks)
        {
            return Err(CompositionEquilibriumErrorV1::Invalid(
                "invalid equilibrium tolerances or work limits",
            ));
        }
        Ok(())
    }
}
impl EquilibriumEnthalpySolveV1 {
    fn validate(self) -> Result<(), CompositionEquilibriumErrorV1> {
        self.equilibrium.validate()?;
        if !self.minimum_temperature_k.is_finite()
            || self.minimum_temperature_k <= 0.0
            || !self.maximum_temperature_k.is_finite()
            || self.maximum_temperature_k <= self.minimum_temperature_k
            || !self.enthalpy_tolerance_j.is_finite()
            || self.enthalpy_tolerance_j <= 0.0
            || !(1..=256).contains(&self.maximum_temperature_iterations)
        {
            return Err(CompositionEquilibriumErrorV1::Invalid(
                "invalid HP bracket, energy tolerance or budget",
            ));
        }
        Ok(())
    }
}
fn allocated<T: Clone>(len: usize, value: T) -> Result<Vec<T>, CompositionEquilibriumErrorV1> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(len)
        .map_err(|_| CompositionEquilibriumErrorV1::AllocationRefused)?;
    result.resize(len, value);
    Ok(result)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CompositionEquilibriumErrorV1> {
    if cancelled() {
        Err(CompositionEquilibriumErrorV1::Cancelled)
    } else {
        Ok(())
    }
}
#[derive(Default, Clone, Copy)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let sum = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - sum) + value
        } else {
            (value - sum) + self.sum
        };
        self.sum = sum;
    }
    fn value(self) -> f64 {
        self.sum + self.correction
    }
}
#[cfg(test)]
mod tests;
