//! Electron cross sections and bounded Maxwellian collision moments.
//!
//! Incident energy is per electron in joules; cross section is in m². The
//! stationary massive-target approximation gives the event coefficient
//! `k = sqrt(8/(pi me)) (kB Te)^(-3/2) integral(E sigma exp(-E/kB Te) dE)`.
//! Multiplication by target number density gives a per-electron event rate.
//! Molar mass-action coefficients instead use `NA*k`, not `k`.
//!
//! A small displaced Maxwellian gives the distinct momentum-drag coefficient
//! `km = (2/3) <sigma_mt v E/(kB Te)>`. This is a fluid momentum closure,
//! not the exact energy-dependent Lorentz-Boltzmann electrical mobility.
//! It assumes `Te/me >> Th/m_target`; consumers must admit that approximation
//! against their actual heavy-particle state. Ionization cross sections cannot
//! be used as elastic momentum-transfer cross sections.
//!
//! The analytic moments and their temperature derivatives follow directly
//! by differentiating the Maxwellian under the energy integral. For the
//! momentum weighting see Cao et al., Earth Planet. Phys. 10, 82–91 (2026),
//! DOI 10.26464/epp2026017, equation (1), and Banks, Planet. Space Sci. 14,
//! 1085–1103 (1966), DOI 10.1016/0032-0633(66)90024-9.
//!
//! Tabulation interpolation is an explicit piecewise-linear source model.
//! A caller-supplied upper envelope beyond its final energy bounds the omitted
//! Maxwellian tail; no missing values are silently replaced with zero.
//! Quadrature differences are numerical estimates, not interval certificates
//! or experimental cross-section uncertainties. Successful trial evaluation
//! allocates nothing; the publication evaluator retains the original source.

use crate::SpeciesId;
use fs_matdb::Provenance;
use fs_math::det::{exp, ln};
use fs_qty::{Area, Energy, Temperature};

mod integrate;
/// Published NIST SRD 107 generalized BEB ionization models.
pub mod nist;

/// CODATA 2022 electron mass [kg], used in every velocity average here.
pub const ELECTRON_MASS_KG: f64 = 9.109_383_713_9e-31;
/// Exact SI Boltzmann constant [J/K].
pub const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;
/// Exact SI Avogadro constant [1/mol].
pub const AVOGADRO_PER_MOL: f64 = 6.022_140_76e23;
/// Exact conversion of electronvolts to joules.
pub const ELECTRONVOLT_J: f64 = 1.602_176_634e-19;
/// Maximum retained tabulation size; construction never grows it implicitly.
pub const MAX_CROSS_SECTION_POINTS_V1: usize = 4096;

/// Physical quantity represented by the original cross section.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ElectronCollisionKindV1 {
    /// Total single-ionization event cross section, not momentum relaxation.
    Ionization,
    /// Elastic momentum-transfer cross section, including its angular weight.
    ElasticMomentumTransfer,
}

/// Maxwellian integral selected by its consuming physical equation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaxwellianRateWeightingV1 {
    /// Ordinary event average `<sigma v>`.
    Event,
    /// Linear momentum drag of a displaced Maxwellian, `(2/3)<sigma v E/kT>`.
    WeakDriftMomentum,
}

/// Explicit numerical work and error-estimate controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaxwellianRateControlsV1 {
    /// Upper bound for each relative quadrature estimate plus omitted-tail
    /// bound, for both the rate moment and its temperature-derivative moment.
    /// This does not bound uncertainty of the original physical model.
    pub relative_tolerance: f64,
    /// Maximum Gauss integration panels in one evaluation, in `1..=16384`.
    pub maximum_subdivisions: usize,
}
impl MaxwellianRateControlsV1 {
    pub(crate) fn validate(self) -> Result<(), ElectronCollisionErrorV1> {
        if !self.relative_tolerance.is_finite()
            || !(64.0 * f64::EPSILON..0.1).contains(&self.relative_tolerance)
            || !(1..=16_384).contains(&self.maximum_subdivisions)
        {
            return Err(ElectronCollisionErrorV1::Invalid("Maxwellian integration controls"));
        }
        Ok(())
    }
}

/// Refused source, unsupported domain, unresolved tail, work budget or range.
#[derive(Clone, Debug, PartialEq)]
pub enum ElectronCollisionErrorV1 {
    /// Invalid physical quantities, provenance, weighting or controls.
    Invalid(&'static str),
    /// Direct cross-section query exceeds the original source interval.
    OutsideEnergyDomain,
    /// Available source coverage/envelope cannot resolve the requested rate.
    UnresolvedTail {
        /// Conservative relative upper bound for the omitted rate/moment.
        relative_bound: f64,
    },
    /// Integration did not resolve the source within the explicit panel bound.
    QuadratureBudget,
    /// A logarithm or needed moment cannot be represented in binary64.
    Arithmetic(&'static str),
    /// Cancellation before a complete result was available.
    Cancelled,
}
impl core::fmt::Display for ElectronCollisionErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "electron collision: {self:?}")
    }
}
impl std::error::Error for ElectronCollisionErrorV1 {}

/// One source sample, with incident energy per electron and SI cross section.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectronCrossSectionPointV1 {
    /// Incident electron energy [J].
    pub energy: Energy,
    /// Original cross section [m²].
    pub cross_section: Area,
}

#[derive(Clone, Debug, PartialEq)]
enum CrossSectionLaw {
    Linear(std::sync::Arc<[ElectronCrossSectionPointV1]>),
    Beb {
        binding_energy: f64,
        orbital_kinetic_energy: f64,
        occupation: f64,
        dipole_constant: f64,
    },
}

/// Immutable source-bound collision law for one explicitly named target.
#[derive(Clone, Debug, PartialEq)]
pub struct ElectronCrossSectionV1 {
    target: SpeciesId,
    kind: ElectronCollisionKindV1,
    provenance: Provenance,
    threshold: Energy,
    maximum_energy: Energy,
    high_energy_upper_bound: Area,
    scale: f64,
    law: CrossSectionLaw,
}
impl ElectronCrossSectionV1 {
    /// Retain a complete piecewise-linear table on `[threshold, maximum]`.
    /// Ionization requires a positive threshold and zero cross section there.
    /// Elastic momentum-transfer tables begin at zero incident energy.
    /// `high_energy_upper_bound` is an explicit source/caller assertion that
    /// `0 <= sigma(E) <= bound` above the last sample, not an extrapolation.
    /// Its provenance must support the intended use. Zero declares compact
    /// support and requires a zero final sample; it is never chosen implicitly.
    pub fn tabulated(
        target: SpeciesId,
        kind: ElectronCollisionKindV1,
        provenance: Provenance,
        points: Vec<ElectronCrossSectionPointV1>,
        high_energy_upper_bound: Area,
    ) -> Result<Self, ElectronCollisionErrorV1> {
        if provenance.source.trim().is_empty() || provenance.license.trim().is_empty() {
            return Err(ElectronCollisionErrorV1::Invalid("collision provenance and usage terms"));
        }
        if !(2..=MAX_CROSS_SECTION_POINTS_V1).contains(&points.len())
            || !nonnegative(high_energy_upper_bound.value())
            || points.iter().any(|p| !nonnegative(p.energy.value()) || !nonnegative(p.cross_section.value()))
            || points.windows(2).any(|p| p[0].energy.value() >= p[1].energy.value())
        {
            return Err(ElectronCollisionErrorV1::Invalid("ordered SI cross-section table"));
        }
        let first = points[0];
        let last = points[points.len() - 1];
        match kind {
            ElectronCollisionKindV1::Ionization
                if !positive(first.energy.value()) || first.cross_section.value() != 0.0 =>
                return Err(ElectronCollisionErrorV1::Invalid("ionization threshold must have zero cross section")),
            ElectronCollisionKindV1::ElasticMomentumTransfer if first.energy.value() != 0.0 =>
                return Err(ElectronCollisionErrorV1::Invalid("momentum-transfer table must cover zero energy")),
            _ => {}
        }
        if high_energy_upper_bound.value() < last.cross_section.value() {
            return Err(ElectronCollisionErrorV1::Invalid("tail envelope below final cross section"));
        }
        let scale = points.iter().map(|p| p.cross_section.value()).fold(high_energy_upper_bound.value(), f64::max);
        Ok(Self {
            target, kind, provenance, threshold: first.energy,
            maximum_energy: last.energy, high_energy_upper_bound, scale,
            law: CrossSectionLaw::Linear(points.into()),
        })
    }

    /// Original target identity; no element or charge is inferred from its name.
    pub const fn target_species(&self) -> &SpeciesId { &self.target }
    /// Original physical cross-section meaning.
    pub const fn kind(&self) -> ElectronCollisionKindV1 { self.kind }
    /// Original source citation, usage terms and any acquired artifact identity.
    pub const fn provenance(&self) -> &Provenance { &self.provenance }
    /// Zero cross section below this threshold [J per incident electron].
    pub const fn threshold(&self) -> Energy { self.threshold }
    /// Upper energy of original direct source evaluation [J].
    pub const fn maximum_energy(&self) -> Energy { self.maximum_energy }
    /// Declared upper envelope of the unknown/unintegrated high energy tail.
    pub const fn high_energy_upper_bound(&self) -> Area { self.high_energy_upper_bound }
    /// Original interpolation samples, absent for the analytic NIST models.
    pub fn points(&self) -> Option<&[ElectronCrossSectionPointV1]> {
        match &self.law { CrossSectionLaw::Linear(p) => Some(p), _ => None }
    }

    /// Evaluate only within the original covered source interval. Below an
    /// ionization threshold the result is exactly zero, not extrapolated.
    pub fn cross_section(&self, energy: Energy) -> Result<Area, ElectronCollisionErrorV1> {
        if !nonnegative(energy.value()) {
            return Err(ElectronCollisionErrorV1::Invalid("incident energy"));
        }
        if energy.value() > self.maximum_energy.value() {
            return Err(ElectronCollisionErrorV1::OutsideEnergyDomain);
        }
        Ok(Area::new(self.value(energy.value())))
    }

    fn value(&self, energy: f64) -> f64 {
        if energy < self.threshold.value() { return 0.0; }
        match &self.law {
            CrossSectionLaw::Linear(points) => {
                let upper = points.partition_point(|p| p.energy.value() <= energy).min(points.len() - 1).max(1);
                let a = points[upper - 1];
                let b = points[upper];
                let f = (energy - a.energy.value()) / (b.energy.value() - a.energy.value());
                (1.0 - f) * a.cross_section.value() + f * b.cross_section.value()
            }
            CrossSectionLaw::Beb { binding_energy, orbital_kinetic_energy, occupation, dipole_constant } => {
                nist::beb_value(energy, *binding_energy, *orbital_kinetic_energy, *occupation, *dipole_constant)
            }
        }
    }

    /// Allocation-free Maxwellian trial evaluation. Event and drag meanings
    /// cannot be interchanged. Finite logarithmic rates survive direct
    /// underflow; an unrepresentable direct value is `None`, never a false zero.
    pub fn evaluate_properties(
        &self,
        electron_temperature: Temperature,
        weighting: MaxwellianRateWeightingV1,
        controls: MaxwellianRateControlsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaxwellianRatePropertiesV1, ElectronCollisionErrorV1> {
        controls.validate()?;
        if !positive(electron_temperature.value()) {
            return Err(ElectronCollisionErrorV1::Invalid("electron temperature"));
        }
        if matches!((self.kind, weighting),
            (ElectronCollisionKindV1::Ionization, MaxwellianRateWeightingV1::WeakDriftMomentum)
            | (ElectronCollisionKindV1::ElasticMomentumTransfer, MaxwellianRateWeightingV1::Event)) {
            return Err(ElectronCollisionErrorV1::Invalid("cross-section kind does not match rate weighting"));
        }
        if cancelled() { return Err(ElectronCollisionErrorV1::Cancelled); }
        if self.scale == 0.0 {
            return Ok(MaxwellianRatePropertiesV1 { log_coefficient: f64::NEG_INFINITY,
                logarithmic_temperature_derivative: 0.0, coefficient: Some(0.0),
                derivative: Some(0.0), relative_quadrature_error: 0.0, relative_tail_bound: 0.0 });
        }
        let t = electron_temperature.value();
        let thermal_energy = BOLTZMANN_J_PER_K * t;
        let a = self.threshold.value() / thermal_energy;
        let normalization = a + 1.0;
        if !positive(thermal_energy) || !a.is_finite() || !positive(normalization) {
            return Err(ElectronCollisionErrorV1::Arithmetic("threshold-to-thermal-energy ratio"));
        }
        let power = if weighting == MaxwellianRateWeightingV1::Event { 1 } else { 2 };
        let integral = integrate::moments(self, thermal_energy, a, power, controls, &mut cancelled)?;
        let log_coefficient = 0.5 * (ln(8.0 / core::f64::consts::PI) + ln(thermal_energy) - ln(ELECTRON_MASS_KG))
            - a + ln(self.scale) + power as f64 * ln(normalization) + ln(integral.value[power - 1])
            + if power == 2 { ln(2.0 / 3.0) } else { 0.0 };
        let logarithmic_temperature_derivative =
            (normalization * (integral.value[power] / integral.value[power - 1]) - (power as f64 + 0.5)) / t;
        if !log_coefficient.is_finite() || !logarithmic_temperature_derivative.is_finite() {
            return Err(ElectronCollisionErrorV1::Arithmetic("Maxwellian rate or derivative logarithm"));
        }
        let coefficient = direct(log_coefficient);
        let derivative = if logarithmic_temperature_derivative == 0.0 { Some(0.0) } else {
            direct(log_coefficient + ln(logarithmic_temperature_derivative.abs()))
                .map(|d| d.copysign(logarithmic_temperature_derivative))
        };
        Ok(MaxwellianRatePropertiesV1 { log_coefficient, logarithmic_temperature_derivative,
            coefficient, derivative, relative_quadrature_error: integral.relative_error,
            relative_tail_bound: integral.relative_tail })
    }

    /// Publish one evaluated rate while retaining its exact original source,
    /// electron temperature, physical weighting and integration controls.
    pub fn evaluate(
        &self,
        electron_temperature: Temperature,
        weighting: MaxwellianRateWeightingV1,
        controls: MaxwellianRateControlsV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ElectronCollisionEvaluationV1, ElectronCollisionErrorV1> {
        let properties = self.evaluate_properties(electron_temperature, weighting, controls, cancelled)?;
        Ok(ElectronCollisionEvaluationV1 { source: self.clone(), electron_temperature, weighting, controls, properties })
    }
}

/// Numerical coefficient and derivative; no duplicated source metadata.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaxwellianRatePropertiesV1 {
    log_coefficient: f64,
    logarithmic_temperature_derivative: f64,
    coefficient: Option<f64>,
    derivative: Option<f64>,
    relative_quadrature_error: f64,
    relative_tail_bound: f64,
}
impl MaxwellianRatePropertiesV1 {
    /// Direct rate coefficient [m³/s], or `None` outside binary64 range.
    pub const fn coefficient_m3_s(&self) -> Option<f64> { self.coefficient }
    /// Temperature derivative [m³/(s K)], or `None` outside binary64 range.
    pub const fn derivative_m3_s_k(&self) -> Option<f64> { self.derivative }
    /// Natural logarithm of the coherent-SI rate coefficient.
    pub const fn log_coefficient_m3_s(&self) -> f64 { self.log_coefficient }
    /// `d ln(k)/d Te` [1/K], evaluated with the same energy integral.
    pub const fn logarithmic_temperature_derivative_per_k(&self) -> f64 { self.logarithmic_temperature_derivative }
    /// Numerical relative quadrature estimate for rate and derivative moments.
    pub const fn estimated_relative_quadrature_error(&self) -> f64 { self.relative_quadrature_error }
    /// Conservative omitted-tail relative bound, conditional on the declared
    /// source envelope; applies to both rate and derivative integral moments.
    pub const fn relative_tail_bound(&self) -> f64 { self.relative_tail_bound }
}

/// Immutable publication result bound to original source and electron state.
#[derive(Clone, Debug, PartialEq)]
pub struct ElectronCollisionEvaluationV1 {
    source: ElectronCrossSectionV1,
    electron_temperature: Temperature,
    weighting: MaxwellianRateWeightingV1,
    controls: MaxwellianRateControlsV1,
    properties: MaxwellianRatePropertiesV1,
}
impl ElectronCollisionEvaluationV1 {
    /// Original source, interpolation/model parameters and provenance.
    pub const fn source(&self) -> &ElectronCrossSectionV1 { &self.source }
    /// Original target identity.
    pub const fn target_species(&self) -> &SpeciesId { self.source.target_species() }
    /// Actual electron temperature used by the Maxwellian average.
    pub const fn electron_temperature(&self) -> Temperature { self.electron_temperature }
    /// Physical mass used by the speed average, for source-bank admission.
    pub const fn electron_mass_kg(&self) -> f64 { ELECTRON_MASS_KG }
    /// Event versus displaced-Maxwellian momentum weighting.
    pub const fn weighting(&self) -> MaxwellianRateWeightingV1 { self.weighting }
    /// Original numerical admission controls.
    pub const fn controls(&self) -> MaxwellianRateControlsV1 { self.controls }
    /// Complete numerical rate and integration diagnostics.
    pub const fn properties(&self) -> MaxwellianRatePropertiesV1 { self.properties }
    /// Direct rate coefficient [m³/s], when representable.
    pub const fn coefficient_m3_s(&self) -> Option<f64> { self.properties.coefficient_m3_s() }
    /// Direct derivative [m³/(s K)], when representable.
    pub const fn derivative_m3_s_k(&self) -> Option<f64> { self.properties.derivative_m3_s_k() }
    /// Natural logarithm of the coherent-SI coefficient.
    pub const fn log_coefficient_m3_s(&self) -> f64 { self.properties.log_coefficient_m3_s() }
    /// `d ln(k)/d Te` [1/K].
    pub const fn logarithmic_temperature_derivative_per_k(&self) -> f64 { self.properties.logarithmic_temperature_derivative_per_k() }
}

fn direct(log: f64) -> Option<f64> {
    let value = exp(log);
    positive(value).then_some(value)
}
fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn nonnegative(value: f64) -> bool { value.is_finite() && value >= 0.0 }

#[cfg(test)]
mod tests;
