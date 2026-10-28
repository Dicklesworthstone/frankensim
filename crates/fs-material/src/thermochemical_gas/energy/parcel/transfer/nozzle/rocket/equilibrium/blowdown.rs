//! Finite, well-mixed adiabatic reservoirs with instantaneous chemical equilibrium.
//!
//! Reaction conserves elements; well-mixed exhaust removes the same elemental
//! ratios as the remaining gas. At chemical equilibrium the reaction contribution
//! to the Gibbs relation vanishes, and d(m*u)=h*dm implies ds=0 for the remaining
//! gas at fixed volume. Molecular species fractions need NOT stay constant.
//!
//! The original Gibbs/NASA owner supplies this isentrope. Physical tank mass is
//! distinct from its arbitrary feed normalization. One fixed birth reference is
//! retained across all steps, preventing entropy/reference drift on restart.
//! This is not finite-rate chemistry, combustion injection, wall heat or slosh.

use super::*;
use fs_thermochem::equilibrium::composition::isentrope::density::{
    EquilibriumDensityErrorV1, EquilibriumDensityQueryV1, EquilibriumDensitySolveV1,
};
use std::sync::Arc;

/// Source-supported depletion pressure range and physical balance tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumBlowdownSolve {
    /// Positive tank pressure floor [Pa], below its current pressure.
    pub minimum_pressure_pa: f64,
    /// Relative density-inversion tolerance, in (0,0.01).
    pub relative_density_tolerance: f64,
    /// Absolute mass/EOS discrepancy allowed [kg].
    pub mass_tolerance_kg: f64,
    /// Absolute energy arithmetic and mean-exhaust-enthalpy tolerance [J].
    pub energy_tolerance_j: f64,
    /// At most 256 pressure queries per requested mass.
    pub maximum_density_iterations: usize,
}
impl EquilibriumBlowdownSolve {
    fn validate(self, pressure: f64) -> Result<(), EquilibriumBlowdownError> {
        if !positive(self.minimum_pressure_pa) || self.minimum_pressure_pa >= pressure
            || !positive(self.relative_density_tolerance) || self.relative_density_tolerance >= 0.01
            || !positive(self.mass_tolerance_kg) || !positive(self.energy_tolerance_j)
            || !(1..=256).contains(&self.maximum_density_iterations) {
            return Err(EquilibriumBlowdownError::Invalid("invalid equilibrium depletion controls"));
        }
        Ok(())
    }
    fn density(self, maximum_pressure_pa: f64) -> EquilibriumDensitySolveV1 {
        EquilibriumDensitySolveV1 { minimum_pressure_pa: self.minimum_pressure_pa,
            maximum_pressure_pa, relative_density_tolerance: self.relative_density_tolerance,
            maximum_pressure_iterations: self.maximum_density_iterations }
    }
}

/// Complete refusal before a reservoir candidate is published.
#[derive(Clone, Debug, PartialEq)]
pub enum EquilibriumBlowdownError {
    /// Invalid inventory, physical domain or unrepresentable scalar.
    Invalid(&'static str),
    /// Original equilibrium/source/cancellation failure.
    Thermochemistry(EquilibriumExpansionErrorV1),
    /// Original equilibrium density-inversion failure.
    Density(EquilibriumDensityErrorV1),
    /// Actual fixed-volume mass differs from the source EOS beyond tolerance.
    DensityBalance { residual_kg: f64 },
    /// Actual energy/mean-enthalpy accounting failed the declared budget.
    EnergyBalance { residual_j: f64 },
    /// Cancellation at the finite-inventory boundary.
    Cancelled,
}
impl From<EquilibriumExpansionErrorV1> for EquilibriumBlowdownError {
    fn from(value: EquilibriumExpansionErrorV1) -> Self { Self::Thermochemistry(value) }
}
impl From<EquilibriumDensityErrorV1> for EquilibriumBlowdownError {
    fn from(value: EquilibriumDensityErrorV1) -> Self { Self::Density(value) }
}
impl core::fmt::Display for EquilibriumBlowdownError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "equilibrium reservoir blowdown: {self:?}")
    }
}
impl std::error::Error for EquilibriumBlowdownError {}
fn check_cancel(cancelled: &mut impl FnMut() -> bool) -> Result<(), EquilibriumBlowdownError> {
    if cancelled() { Err(EquilibriumBlowdownError::Cancelled) } else { Ok(()) }
}

#[derive(Clone, Debug, PartialEq)]
struct Source {
    chemistry: Arc<IdealGasEquilibriumV1>,
    feed: Vec<(SpeciesId, Amount)>,
    thermochemistry: EquilibriumEntropySolveV1,
    reference: EquilibriumThermodynamicStateV1,
    basis_mass_kg: f64,
    reference_enthalpy_j_kg: f64,
}

/// Immutable physical reservoir. Clones share the original chemistry/feed source.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumGasReservoir {
    source: Arc<Source>,
    volume_m3: f64,
    mass_kg: f64,
    state: EquilibriumThermodynamicStateV1,
    sample: EquilibriumIsentropeSampleV1,
}

/// Physical extensive exchange. Individual species are not conserved in this
/// reacting system; differences of tank species inventories are NOT exhaust fluxes.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumBlowdownStep {
    /// Actual mass_before-mass_after [kg].
    pub expelled_mass_kg: f64,
    /// Integral of outgoing stagnation enthalpy, U_before-U_after [J]. Signed
    /// formation references survive. Nozzle kinetic energy is included, not extra.
    pub exhausted_enthalpy_j: f64,
    /// Exhaust mass-average stagnation enthalpy [J/kg]; absent for no withdrawal.
    pub mean_stagnation_enthalpy_j_kg: Option<f64>,
    /// Reference-relative correction to dm*h_before [J], avoiding a large-datum subtraction.
    pub enthalpy_correction_j: f64,
    /// Arithmetic closure of the extensive energy expression [J].
    pub energy_arithmetic_residual_j: f64,
    /// Actual endpoint rho*V minus assigned remaining mass [kg].
    pub density_mass_residual_kg: f64,
    /// Endpoint specific entropy minus the fixed BIRTH reference [J/(kg K)].
    /// This bound does not accumulate with the number of accepted intervals.
    pub specific_entropy_residual_j_kg_k: f64,
    /// Number of actual pressure queries in the density solve.
    pub pressure_queries: usize,
}

/// The next source-bearing reservoir and exactly its exchanged energy/mass.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumBlowdownCandidate {
    endpoint: EquilibriumGasReservoir,
    step: EquilibriumBlowdownStep,
}
impl EquilibriumBlowdownCandidate {
    /// Tentative next reservoir; the original remains unchanged.
    pub const fn endpoint(&self) -> &EquilibriumGasReservoir { &self.endpoint }
    /// Same withdrawal's integrated exchange.
    pub const fn step(&self) -> &EquilibriumBlowdownStep { &self.step }
    /// Accept the physical state and its accounts together.
    pub fn into_parts(self) -> (EquilibriumGasReservoir, EquilibriumBlowdownStep) {
        (self.endpoint, self.step)
    }
}

impl EquilibriumGasReservoir {
    /// Assemble from one original equilibrium source, explicit volume and T,p.
    /// Feed AMOUNTS set elemental ratios and numerical normalization; actual tank
    /// mass is rho(T,p,elements)*volume, never an independently supplied number.
    /// The thermochemical controls and birth reference stay fixed through history.
    pub fn new(
        chemistry: Arc<IdealGasEquilibriumV1>, feed: Vec<(SpeciesId, Amount)>,
        volume_m3: f64, temperature: Temperature, pressure: Pressure,
        thermochemistry: EquilibriumEntropySolveV1, mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, EquilibriumBlowdownError> {
        check_cancel(&mut cancelled)?;
        if !positive(volume_m3) { return Err(EquilibriumBlowdownError::Invalid("positive finite tank volume required")); }
        let (state, sample, basis_mass_kg) = {
            let path = EquilibriumIsentropeV1::new(&chemistry, temperature, pressure, &feed,
                thermochemistry, &mut cancelled)?;
            (path.reference_state().clone(), path.current_sample().expect("new isentrope has its reference"), path.feed_mass_kg())
        };
        let mass_kg = sample.response.density_kg_m3 * volume_m3;
        let reference_enthalpy_j_kg = state.equilibrium().enthalpy().value() / basis_mass_kg;
        let source = Arc::new(Source { chemistry, feed, thermochemistry, reference: state.clone(),
            basis_mass_kg, reference_enthalpy_j_kg });
        let reservoir = Self { source, volume_m3, mass_kg, state, sample };
        reservoir.validate_inventory()?;
        check_cancel(&mut cancelled)?;
        Ok(reservoir)
    }
    /// Actual physical inventory [kg].
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    /// Fixed physical volume [m3].
    pub const fn volume_m3(&self) -> f64 { self.volume_m3 }
    /// Original source-bearing equilibrium state, on the RETAINED FEED BASIS.
    /// Its extensive amounts/H/S are not the tank's amounts/H/S. Multiply by
    /// physical_amount_scale() for physical amounts; T,p and fractions are intensive.
    pub const fn thermodynamics(&self) -> &EquilibriumThermodynamicStateV1 { &self.state }
    /// Actual tank mass divided by the conserved feed-basis mass.
    pub fn physical_amount_scale(&self) -> f64 { self.mass_kg / self.source.basis_mass_kg }
    /// Actual physical tank internal energy [J] including formation references.
    /// The explicit fixed-volume identity is U=m*h-p*V.
    pub fn internal_energy_j(&self) -> f64 {
        self.mass_kg * self.specific_enthalpy_j_kg() - self.sample.pressure.value() * self.volume_m3
    }
    /// Original species formation reference plus current relative enthalpy [J/kg].
    pub fn specific_enthalpy_j_kg(&self) -> f64 {
        self.source.reference_enthalpy_j_kg + self.sample.specific_enthalpy_change_j_kg
    }
    /// Immutable numerical/source policy used by this whole reservoir history.
    pub fn thermochemical_controls(&self) -> EquilibriumEntropySolveV1 {
        self.source.thermochemistry
    }

    fn validate_inventory(&self) -> Result<(), EquilibriumBlowdownError> {
        if !positive(self.mass_kg) || !positive(self.physical_amount_scale())
            || !self.specific_enthalpy_j_kg().is_finite() || !self.internal_energy_j().is_finite() {
            return Err(EquilibriumBlowdownError::Invalid("unrepresentable equilibrium reservoir inventory or energy"));
        }
        Ok(())
    }
    fn path(&self, cancelled: &mut impl FnMut() -> bool)
        -> Result<EquilibriumIsentropeV1<'_>, EquilibriumBlowdownError> {
        let reference = self.source.reference.equilibrium();
        let path = EquilibriumIsentropeV1::new(&self.source.chemistry,
            reference.temperature(), reference.pressure(), &self.source.feed,
            self.source.thermochemistry, cancelled)?;
        if path.reference_state() != &self.source.reference || path.feed_mass_kg() != self.source.basis_mass_kg {
            return Err(EquilibriumBlowdownError::Invalid("original isentrope reference changed"));
        }
        Ok(path)
    }
    fn density_at_mass(&self, path: &mut EquilibriumIsentropeV1<'_>, mass: f64,
        c: EquilibriumBlowdownSolve, cancelled: &mut impl FnMut() -> bool)
        -> Result<EquilibriumDensityQueryV1, EquilibriumBlowdownError> {
        let query = path.sample_density(mass / self.volume_m3,
            c.density(self.sample.pressure.value()), cancelled)?;
        let residual = query.sample.response.density_kg_m3 * self.volume_m3 - mass;
        if !residual.is_finite() || residual.abs() > c.mass_tolerance_kg {
            return Err(EquilibriumBlowdownError::DensityBalance { residual_kg: residual });
        }
        Ok(query)
    }
    fn candidate_at_mass(&self, path: &mut EquilibriumIsentropeV1<'_>, mass: f64,
        c: EquilibriumBlowdownSolve, cancelled: &mut impl FnMut() -> bool)
        -> Result<EquilibriumBlowdownCandidate, EquilibriumBlowdownError> {
        if !positive(mass) || mass >= self.mass_kg {
            return Err(EquilibriumBlowdownError::Invalid("withdrawal must leave a distinct positive mass"));
        }
        let query = self.density_at_mass(path, mass, c, cancelled)?;
        let sample = query.sample;
        let dm = self.mass_kg - mass;
        let dh = sample.specific_enthalpy_change_j_kg - self.sample.specific_enthalpy_change_j_kg;
        let dp = sample.pressure.value() - self.sample.pressure.value();
        // U=m*h-p*V. Keep all finite changes on one common reference, rather
        // than subtracting two large extensive formation-energy inventories.
        let correction = -mass * dh + self.volume_m3 * dp;
        let baseline = dm * self.specific_enthalpy_j_kg();
        let exhausted = baseline + correction;
        let arithmetic = (exhausted - baseline) - correction;
        let mean = self.specific_enthalpy_j_kg() + correction / dm;
        if !dh.is_finite() || dh >= 0.0 || !dp.is_finite() || dp >= 0.0
            || !correction.is_finite() || !exhausted.is_finite() || !mean.is_finite() {
            return Err(EquilibriumBlowdownError::Invalid("invalid reacting adiabatic depletion or exported energy"));
        }
        let defect = arithmetic.abs().max(correction).max(dm * dh - correction);
        if !defect.is_finite() || defect > c.energy_tolerance_j {
            return Err(EquilibriumBlowdownError::EnergyBalance { residual_j: defect });
        }
        let state = path.publish_current(&mut *cancelled)?;
        let endpoint = Self { source: self.source.clone(), volume_m3: self.volume_m3,
            mass_kg: mass, state, sample };
        endpoint.validate_inventory()?;
        check_cancel(cancelled)?;
        Ok(EquilibriumBlowdownCandidate { endpoint, step: EquilibriumBlowdownStep {
            expelled_mass_kg: dm, exhausted_enthalpy_j: exhausted,
            mean_stagnation_enthalpy_j_kg: Some(mean), enthalpy_correction_j: correction,
            energy_arithmetic_residual_j: arithmetic,
            density_mass_residual_kg: sample.response.density_kg_m3 * self.volume_m3 - mass,
            specific_entropy_residual_j_kg_k: sample.entropy_change_j_k / self.source.basis_mass_kg,
            pressure_queries: query.pressure_queries,
        } })
    }

    /// Prescribed finite withdrawal with actual chemical re-equilibration.
    /// No mass rate is inferred. Export the returned energy once; adding another
    /// reaction-heat or exhaust kinetic-energy debit would count it twice.
    /// Every failure/cancellation leaves the original mass and source state intact.
    pub fn exhaust_mass_adiabatic(&self, expelled_mass_kg: f64,
        controls: EquilibriumBlowdownSolve, mut cancelled: impl FnMut() -> bool)
        -> Result<EquilibriumBlowdownCandidate, EquilibriumBlowdownError> {
        controls.validate(self.sample.pressure.value())?;
        check_cancel(&mut cancelled)?;
        if !expelled_mass_kg.is_finite() || expelled_mass_kg < 0.0 || expelled_mass_kg >= self.mass_kg {
            return Err(EquilibriumBlowdownError::Invalid("withdrawal must be nonnegative and below tank inventory"));
        }
        if expelled_mass_kg == 0.0 {
            return Ok(EquilibriumBlowdownCandidate { endpoint: self.clone(), step: EquilibriumBlowdownStep {
                expelled_mass_kg: 0.0, exhausted_enthalpy_j: 0.0, mean_stagnation_enthalpy_j_kg: None,
                enthalpy_correction_j: 0.0, energy_arithmetic_residual_j: 0.0,
                density_mass_residual_kg: self.sample.response.density_kg_m3 * self.volume_m3 - self.mass_kg,
                specific_entropy_residual_j_kg_k: self.sample.entropy_change_j_k / self.source.basis_mass_kg,
                pressure_queries: 0,
            } });
        }
        let mut path = self.path(&mut cancelled)?;
        self.candidate_at_mass(&mut path, self.mass_kg - expelled_mass_kg, controls, &mut cancelled)
    }
}
