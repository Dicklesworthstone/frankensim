//! Source-bound transfer from the original reactive solid to a gas compartment.
//!
//! [`VolatileMixture`] implements the solid owner's [`VolatileEnthalpy`] using
//! the original frozen-mixture/NASA evaluator, including formation enthalpy.
//! Bind it as `PyrolysisMaterial::volatile` when constructing the solid, and
//! supply virgin/char enthalpies on that same physical reference. An effective
//! heat-of-pyrolysis chart with a different datum cannot be silently substituted.
//!
//! Each original [`VolatileRelease`] contains a time integral of outgoing mass
//! and enthalpy, possibly over many temperatures. The adapter retains both
//! integrals. It maps mass onto the bound canonical species axis and subtracts
//! those species' original internal-energy references for [`LayerRates`].
//! No final-solid-temperature proxy, extra heat of pyrolysis, heat of combustion,
//! soot gas species, or instantaneous reaction is introduced.
//!
//! In a `ZoneNetwork::step_with_sources` prepare callback:
//! 1. Check that accepted solid and room clocks agree.
//! 2. Advance the solid *candidate* to `solid.time_s() + duration.value()`.
//! 3. Convert each cell release with its bound `VolatileMixture::release_rates`,
//!    and add the result to the receiving layer's original `LayerRates`.
//! 4. Return the solid candidate and audit as the provider payload.
//! 5. After the network accepts, install that returned candidate by an infallible
//!    move, without another cancellation point. Failed attempts publish neither
//!    owner; a shorter attempt recomputes its solid release from accepted state.
//!
//! Prepared rates already enter the network's pressure equation. Subsequent
//! gas chemistry can use the existing source-inclusive pressure-stage API.
//! Each material/cell release must use the model actually passed to that solid;
//! the plain thermal release record cannot authenticate that caller obligation.

use fs_gas::multicomponent::compartment::{Layer, LayerRates, TwoZone, ZoneError};
use fs_gas::multicomponent::{MAX_SPECIES, MixtureGas};
use fs_qty::{Energy, Temperature, Time};
use fs_thermal::enthalpy::DepositionError;
pub use fs_thermal::enthalpy::pyrolysis::{VolatileEnthalpy, VolatileRelease};
use fs_thermochem::mixture::{FrozenIdealGasMixtureErrorV1, FrozenIdealGasMixtureModelV1};
use fs_thermochem::{Composition, CompositionBasis, MassSpecificEnergyQuantityV1};

/// Numerical transfer gates, separate from source/model uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VolatileTransferControl {
    /// Absolute discrepancy between emitted and rate-integrated mass [kg].
    pub mass_tolerance_kg: f64,
    /// Absolute discrepancy between carried and reconstructed physical enthalpy [J].
    pub enthalpy_tolerance_j: f64,
}
impl Default for VolatileTransferControl {
    fn default() -> Self {
        Self {
            mass_tolerance_kg: 1e-12,
            enthalpy_tolerance_j: 1e-7,
        }
    }
}

/// Source, chart, conservation or cancellation refusal; no owner is mutated.
#[derive(Clone, Debug, PartialEq)]
pub enum VolatileGasError {
    /// Invalid mass composition, release, duration or acceptance budget.
    Invalid(&'static str),
    /// Original full gas bank, caloric domain or energy references differ.
    SourceMismatch,
    /// Original frozen-mixture or NASA owner refused construction/evaluation.
    Mixture(FrozenIdealGasMixtureErrorV1),
    /// Original solid owner refused a candidate.
    Thermal(DepositionError),
    /// Original layer-rate addition refused.
    Zone(ZoneError),
    /// Actual rate arithmetic exceeds the supplied conservation gate.
    Conservation {
        /// Actual injected minus released mass [kg].
        mass_kg: f64,
        /// Reconstructed physical input minus carried enthalpy [J].
        enthalpy_j: f64,
    },
    /// Cancellation before the complete result was returned.
    Cancelled,
}
impl core::fmt::Display for VolatileGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pyrolysis volatile gas transfer: {self:?}")
    }
}
impl std::error::Error for VolatileGasError {}
impl From<FrozenIdealGasMixtureErrorV1> for VolatileGasError {
    fn from(value: FrozenIdealGasMixtureErrorV1) -> Self {
        Self::Mixture(value)
    }
}
impl From<DepositionError> for VolatileGasError {
    fn from(value: DepositionError) -> Self {
        match value {
            DepositionError::Cancelled => Self::Cancelled,
            other => Self::Thermal(other),
        }
    }
}
impl From<ZoneError> for VolatileGasError {
    fn from(value: ZoneError) -> Self {
        match value {
            ZoneError::Cancelled => Self::Cancelled,
            other => Self::Zone(other),
        }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), VolatileGasError> {
    if cancelled() {
        Err(VolatileGasError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, VolatileGasError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(VolatileGasError::Invalid(
            "non-finite volatile transfer arithmetic",
        ))
    }
}
fn sum(values: impl Iterator<Item = f64>) -> Result<f64, VolatileGasError> {
    let (mut total, mut correction) = (0.0, 0.0);
    for value in values {
        let adjusted = finite(value - correction)?;
        let next = finite(total + adjusted)?;
        correction = (next - total) - adjusted;
        total = next;
    }
    Ok(total)
}

/// Immutable volatile composition and calorics bound to a complete original gas.
/// The original frozen-mixture owner validates composition and shared species
/// conventions; no species is inferred from its name or silently normalized.
#[derive(Clone, Debug, PartialEq)]
pub struct VolatileMixture {
    gas: MixtureGas,
    calorics: FrozenIdealGasMixtureModelV1,
    fractions: [f64; MAX_SPECIES],
    enthalpy_range_j_kg: [f64; 2],
}

/// One time-averaged source with the actual integrals it represents.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VolatileInjection {
    /// Positive canonical species rates and species-referenced enthalpy [W].
    pub rates: LayerRates,
    /// Original integrated release, including retained-char accounting.
    pub release: VolatileRelease,
    /// Interval over which the returned constant rates must be applied once.
    pub duration: Time,
    /// Actual rate-integrated energy in the original gas's numerical gauge [J].
    pub referenced_energy: Energy,
    /// Actual injected minus emitted mass [kg], retained without repair.
    pub mass_residual_kg: f64,
    /// Actual reconstructed physical input minus carried enthalpy [J].
    pub enthalpy_residual_j: f64,
}
impl VolatileInjection {
    /// Atomically add to one original layer's source account. Apply each release
    /// once; the caller owns the material/cell routing and acceptance lifecycle.
    pub fn add_to(
        &self,
        layer: Layer,
        rates: &mut [LayerRates; 2],
    ) -> Result<(), VolatileGasError> {
        rates[layer.index()].add(&self.rates)?;
        Ok(())
    }
}

impl VolatileMixture {
    /// Bind explicit mass fractions in the complete gas's canonical species order.
    /// Zero fractions are absent from the frozen evaluator, while every rate
    /// retains the full gas axis. Original exact composition admission applies.
    pub fn new(
        gas: &MixtureGas,
        mass_fractions: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, VolatileGasError> {
        poll(&mut cancelled)?;
        if mass_fractions.len() != gas.species_count() {
            return Err(VolatileGasError::Invalid(
                "volatile mass-fraction axis differs from gas",
            ));
        }
        let composition = Composition::new(CompositionBasis::MassFraction, mass_fractions.to_vec())
            .map_err(FrozenIdealGasMixtureErrorV1::from)?;
        let mut components = Vec::new();
        let mut active = Vec::new();
        for (model, &fraction) in gas.source().models().iter().zip(composition.fractions()) {
            poll(&mut cancelled)?;
            if fraction > 0.0 {
                components.push(model.clone());
                active.push(fraction);
            }
        }
        let active = Composition::new(CompositionBasis::MassFraction, active)
            .map_err(FrozenIdealGasMixtureErrorV1::from)?;
        let calorics = FrozenIdealGasMixtureModelV1::new(components, active)?;
        let mut fractions = [0.0; MAX_SPECIES];
        fractions[..mass_fractions.len()].copy_from_slice(mass_fractions);
        let mut result = Self {
            gas: gas.clone(),
            calorics,
            fractions,
            enthalpy_range_j_kg: [0.0; 2],
        };
        let domain = gas.domain();
        result.enthalpy_range_j_kg = [
            result
                .specific_enthalpy(Temperature::new(domain.minimum_temperature), &mut cancelled)?
                .value(),
            result
                .specific_enthalpy(Temperature::new(domain.maximum_temperature), &mut cancelled)?
                .value(),
        ];
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Complete original gas bank, species axis, domain and reference energies.
    pub const fn gas(&self) -> &MixtureGas {
        &self.gas
    }
    /// Unmodified fractions in the full original gas axis.
    pub fn mass_fractions(&self) -> &[f64] {
        &self.fractions[..self.gas.species_count()]
    }
    /// Original admitted frozen-mixture evaluator, retaining all active sources.
    pub const fn calorics(&self) -> &FrozenIdealGasMixtureModelV1 {
        &self.calorics
    }

    /// Physical enthalpy from original NASA species, including formation terms.
    /// Ideal-gas enthalpy is independent of pressure; the evaluator receives its
    /// own source standard pressure, not a newly chosen physical reservoir.
    pub fn specific_enthalpy(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MassSpecificEnergyQuantityV1, VolatileGasError> {
        poll(&mut cancelled)?;
        let domain = self.gas.domain();
        if !(domain.minimum_temperature..=domain.maximum_temperature).contains(&temperature.value())
        {
            return Err(VolatileGasError::Invalid(
                "volatile enthalpy outside gas caloric domain",
            ));
        }
        let standard_pressure = self.calorics.components()[0]
            .convention()
            .reference_pressure();
        let properties = self
            .calorics
            .evaluate_properties(temperature, standard_pressure)?;
        let h = finite(properties.h().value() / properties.molar_mass().value())?;
        poll(&mut cancelled)?;
        Ok(MassSpecificEnergyQuantityV1::new(h))
    }

    /// Convert the original carried integrals to this exact receiver's rates.
    /// No temperature is inferred from an integrated release. The mean carried
    /// enthalpy must lie in the bound source's admitted caloric range (roundoff
    /// allowed without clipping), and the actual rate arithmetic must satisfy
    /// the supplied mass/energy gates. Empty releases are exactly zero.
    pub fn release_rates(
        &self,
        room: &TwoZone,
        release: VolatileRelease,
        duration: Time,
        controls: VolatileTransferControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VolatileInjection, VolatileGasError> {
        poll(&mut cancelled)?;
        if room.gas() != &self.gas {
            return Err(VolatileGasError::SourceMismatch);
        }
        let dt = duration.value();
        if !dt.is_finite()
            || dt <= 0.0
            || !release.mass_kg.is_finite()
            || release.mass_kg < 0.0
            || !release.enthalpy_j.is_finite()
            || !release.consumed_virgin_mass_kg.is_finite()
            || release.consumed_virgin_mass_kg < release.mass_kg
            || !controls.mass_tolerance_kg.is_finite()
            || controls.mass_tolerance_kg < 0.0
            || !controls.enthalpy_tolerance_j.is_finite()
            || controls.enthalpy_tolerance_j < 0.0
        {
            return Err(VolatileGasError::Invalid(
                "invalid volatile release, duration or conservation gate",
            ));
        }
        if release.mass_kg == 0.0 {
            if release.enthalpy_j != 0.0 || release.consumed_virgin_mass_kg != 0.0 {
                return Err(VolatileGasError::Invalid(
                    "empty volatile release carries material or energy",
                ));
            }
        } else {
            let h = finite(release.enthalpy_j / release.mass_kg)?;
            let [lo, hi] = self.enthalpy_range_j_kg;
            let roundoff = 64.0 * f64::EPSILON * h.abs().max(lo.abs()).max(hi.abs()).max(1.0);
            if h < lo - roundoff || h > hi + roundoff {
                return Err(VolatileGasError::Invalid(
                    "carried volatile enthalpy outside bound source domain",
                ));
            }
        }
        let mut rates = LayerRates::default();
        let n = self.gas.species_count();
        let mut integrated = [0.0; MAX_SPECIES];
        for i in 0..n {
            poll(&mut cancelled)?;
            let mass = finite(release.mass_kg * self.fractions[i])?;
            rates.species_kg_s[i] = finite(mass / dt)?;
            integrated[i] = finite(dt * rates.species_kg_s[i])?;
            if release.mass_kg > 0.0 && self.fractions[i] > 0.0 && integrated[i] == 0.0 {
                return Err(VolatileGasError::Invalid(
                    "positive volatile species mass underflow",
                ));
            }
        }
        let reference_energy = sum(integrated[..n]
            .iter()
            .zip(self.gas.energy_references())
            .map(|(m, r)| m * r))?;
        rates.referenced_energy_w = finite(finite(release.enthalpy_j - reference_energy)? / dt)?;
        let actual_referenced = finite(dt * rates.referenced_energy_w)?;
        let mass_residual_kg = finite(sum(integrated[..n].iter().copied())? - release.mass_kg)?;
        let enthalpy_residual_j =
            finite(finite(actual_referenced + reference_energy)? - release.enthalpy_j)?;
        if mass_residual_kg.abs() > controls.mass_tolerance_kg
            || enthalpy_residual_j.abs() > controls.enthalpy_tolerance_j
        {
            return Err(VolatileGasError::Conservation {
                mass_kg: mass_residual_kg,
                enthalpy_j: enthalpy_residual_j,
            });
        }
        poll(&mut cancelled)?;
        Ok(VolatileInjection {
            rates,
            release,
            duration,
            referenced_energy: Energy::new(actual_referenced),
            mass_residual_kg,
            enthalpy_residual_j,
        })
    }
}

impl VolatileEnthalpy for VolatileMixture {
    /// Fixed-work original caloric query; the solid owner polls cancellation
    /// around each local temperature evaluation (at most 128 original species).
    fn specific_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, DepositionError> {
        self.specific_enthalpy(Temperature::new(temperature_k), || false)
            .map(|h| h.value())
            .map_err(|error| match error {
                VolatileGasError::Cancelled => DepositionError::Cancelled,
                _ => DepositionError::MaterialDomain,
            })
    }
}

#[cfg(test)]
mod tests;
