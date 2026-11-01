//! Source-bound multicomponent vapor/liquid calorics from the original PR flash.
//!
//! Species inventories and their ideal formation/entropy references use the
//! exact feed composition once. Phase redistribution changes ideal mixing
//! entropy and the original EOS departures. This avoids amplifying a rounded
//! phase material-balance residual by an arbitrary species formation datum.
//! No gas transport, mixed Cp or acoustic relaxation closure is inferred.

use super::ThermochemicalPengRobinsonModel;
use crate::MaterialError;
use crate::thermochemical_gas::failure;
use fs_math::det;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::eos::peng_robinson::{
    PengRobinsonMixtureV1, PengRobinsonPhaseRootV1,
    flash::{PengRobinsonFlashErrorV1, PengRobinsonFlashOptionsV1, PengRobinsonFlashV1},
};
use fs_thermochem::mixture::{
    FrozenIdealGasMixtureEvaluationV1, FrozenIdealGasMixturePropertiesV1,
};

/// Fixed-pressure enthalpy and source-increment phase-equilibrium solves.
pub mod enthalpy;
/// Phase-screened homogeneous and separated states, including boundary-crossing enthalpy.
pub mod equilibrium;
/// Conserved mixture inventories with continuing rigid-vessel heat inputs.
pub mod parcel;

/// A source, phase or bounded-work refusal; no partial flash state is returned.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermochemicalFlashError {
    /// Invalid physical input, controls or unrepresentable closure.
    Invalid(&'static str),
    /// Original species-source refusal.
    Material(MaterialError),
    /// Original bounded phase-equilibrium refusal, including absence of a split.
    Phase(PengRobinsonFlashErrorV1),
    /// Cancellation before a bounded evaluation or result publication.
    Cancelled,
    /// Whole-call PT flash budget exhausted, including the reference and endpoints.
    WorkLimit {
        /// Number of attempted original PT flashes.
        flash_evaluations: usize,
    },
    /// The admitted endpoints do not bracket the requested enthalpy.
    OutsideEnthalpyBracket,
    /// The admitted endpoints do not bracket the requested entropy.
    OutsideEntropyBracket,
    /// The remaining temperature interval cannot be resolved at the requested accuracy.
    UnresolvedTemperature,
    /// The pressure endpoints do not bracket the requested conserved volume.
    OutsideVolumeBracket,
    /// The pressure interval cannot resolve both conserved coordinates.
    UnresolvedPressure,
    /// The caller's bounded number of pressure/caloric inversions was exhausted.
    PressureWorkLimit {
        /// Number of attempted pressure/caloric inversions.
        pressure_evaluations: usize,
    },
}
impl core::fmt::Display for ThermochemicalFlashError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermochemical flash refusal: {self:?}")
    }
}
impl std::error::Error for ThermochemicalFlashError {}
impl From<MaterialError> for ThermochemicalFlashError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}
type E = ThermochemicalFlashError;

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), E> {
    if cancelled() {
        Err(E::Cancelled)
    } else {
        Ok(())
    }
}
fn x_log_x(fractions: &[f64]) -> f64 {
    fractions
        .iter()
        .filter(|x| **x > 0.0)
        .map(|x| x * det::ln(*x))
        .sum()
}
fn residual_entropy(root: &PengRobinsonPhaseRootV1, fractions: &[f64], t: f64) -> f64 {
    root.residual_enthalpy_j_mol() / t
        - R * fractions
            .iter()
            .zip(root.log_fugacity_coefficients())
            .map(|(fraction, coefficient)| fraction * coefficient)
            .sum::<f64>()
}

#[derive(Clone, Debug, PartialEq)]
struct CaloricSample<P = PengRobinsonFlashV1> {
    phase: P,
    temperature: Temperature,
    pressure: Pressure,
    eos_evaluation_upper_bound: usize,
    molar_volume: f64,
    liquid_mass_fraction: f64,
    vapor_mass_fraction: f64,
    vapor_volume_fraction: f64,
    residual_h: f64,
    residual_u: f64,
    entropy_correction: f64,
    specific_h: f64,
    specific_u: f64,
    specific_s: f64,
}
impl<P> CaloricSample<P> {
    fn with_phase<Q>(self, phase: Q, eos_evaluation_upper_bound: usize) -> CaloricSample<Q> {
        CaloricSample {
            phase,
            temperature: self.temperature,
            pressure: self.pressure,
            eos_evaluation_upper_bound,
            molar_volume: self.molar_volume,
            liquid_mass_fraction: self.liquid_mass_fraction,
            vapor_mass_fraction: self.vapor_mass_fraction,
            vapor_volume_fraction: self.vapor_volume_fraction,
            residual_h: self.residual_h,
            residual_u: self.residual_u,
            entropy_correction: self.entropy_correction,
            specific_h: self.specific_h,
            specific_u: self.specific_u,
            specific_s: self.specific_s,
        }
    }
}

/// One admitted two-phase state with the original critical/mixing law and
/// complete NASA feed receipt. Fractions retain the original PT owner's
/// numerical material residual; the exact feed is the caloric reference basis.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermochemicalFlashState {
    law: PengRobinsonMixtureV1,
    reference: FrozenIdealGasMixtureEvaluationV1,
    sample: CaloricSample,
}
impl ResolvedThermochemicalFlashState {
    /// Original critical properties, interactions and canonical feed composition.
    #[must_use]
    pub const fn constitutive_law(&self) -> &PengRobinsonMixtureV1 {
        &self.law
    }
    /// Exact feed thermodynamics, conventions and all original species receipts.
    #[must_use]
    pub const fn reference(&self) -> &FrozenIdealGasMixtureEvaluationV1 {
        &self.reference
    }
    /// Original phase roots, mole compositions, vapor mole amount and residuals.
    #[must_use]
    pub const fn flash(&self) -> &PengRobinsonFlashV1 {
        &self.sample.phase
    }
    /// Common absolute temperature [K].
    #[must_use]
    pub fn temperature_k(&self) -> f64 {
        self.sample.phase.temperature().value()
    }
    /// Common pressure [Pa].
    #[must_use]
    pub fn pressure_pa(&self) -> f64 {
        self.sample.phase.pressure().value()
    }
    /// Additive phase volume per mole of the exact original feed [m³/mol].
    #[must_use]
    pub const fn molar_volume_m3_mol(&self) -> f64 {
        self.sample.molar_volume
    }
    /// Conserved feed mass divided by additive phase volume [kg/m³].
    #[must_use]
    pub fn density_kg_m3(&self) -> f64 {
        self.law.molar_mass().value() / self.sample.molar_volume
    }
    /// Liquid mass per unit feed mass; distinct from the liquid mole fraction.
    #[must_use]
    pub const fn liquid_mass_fraction(&self) -> f64 {
        self.sample.liquid_mass_fraction
    }
    /// Vapor mass per unit feed mass; distinct from the vapor mole fraction.
    #[must_use]
    pub const fn vapor_mass_fraction(&self) -> f64 {
        self.sample.vapor_mass_fraction
    }
    /// Vapor volume divided by additive total phase volume.
    #[must_use]
    pub const fn vapor_volume_fraction(&self) -> f64 {
        self.sample.vapor_volume_fraction
    }
    /// Full original-reference specific enthalpy [J/kg].
    #[must_use]
    pub const fn specific_enthalpy_j_kg(&self) -> f64 {
        self.sample.specific_h
    }
    /// Full original-reference specific internal energy [J/kg].
    #[must_use]
    pub const fn specific_internal_energy_j_kg(&self) -> f64 {
        self.sample.specific_u
    }
    /// Full entropy including phase-resolved ideal mixing and PR departures [J/(kg K)].
    #[must_use]
    pub const fn specific_entropy_j_kg_k(&self) -> f64 {
        self.sample.specific_s
    }
}

impl ThermochemicalPengRobinsonModel {
    /// Resolve source-bound PT two-phase calorics without requiring gas transport.
    ///
    /// This retains the original bounded flash's fugacity, material and Gibbs
    /// checks; it does not certify a global minimum. An absent/unresolved split
    /// refuses, including pure-fluid/azeotropic cases with indeterminate amount.
    pub fn evaluate_two_phase_flash(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        options: PengRobinsonFlashOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedThermochemicalFlashState, E> {
        let sample = self.caloric_flash_sample(temperature, pressure, options, &mut cancelled)?;
        self.publish_caloric_flash(sample, &mut cancelled)
    }

    fn caloric_flash_sample(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        options: PengRobinsonFlashOptionsV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CaloricSample, E> {
        poll(cancelled)?;
        let reference = self
            .source()
            .thermodynamics()
            .evaluate_properties(temperature, pressure)
            .map_err(|error| failure("flash feed thermodynamics", error))?;
        let phase = self
            .mechanics()
            .two_phase_flash_pt_cancellable(temperature, pressure, options, &mut *cancelled)
            .map_err(|error| match error {
                PengRobinsonFlashErrorV1::Cancelled => E::Cancelled,
                other => E::Phase(other),
            })?;
        let values =
            self.two_phase_caloric_sample(temperature, pressure, &phase, reference, cancelled)?;
        Ok(values.with_phase(phase, 6 * options.maximum_iterations + 2))
    }

    fn two_phase_caloric_sample(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        phase: &PengRobinsonFlashV1,
        reference: FrozenIdealGasMixturePropertiesV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CaloricSample<()>, E> {
        poll(cancelled)?;
        let beta = phase.vapor_fraction_mol();
        let x = phase.liquid_mole_fractions();
        let y = phase.vapor_mole_fractions();
        let mass = self.mechanics().molar_mass().value();
        let mut liquid_mass = 0.0;
        let mut vapor_mass = 0.0;
        let mut liquid_ideal_cv = 0.0;
        let mut vapor_ideal_cv = 0.0;
        for (index, species) in self
            .source()
            .thermodynamics()
            .components()
            .iter()
            .enumerate()
        {
            poll(cancelled)?;
            let cp = species
                .evaluate_properties(temperature)
                .map_err(|error| failure("flash phase caloric stability", error))?
                .cp()
                .value();
            liquid_mass += x[index] * species.molar_mass().value();
            vapor_mass += y[index] * species.molar_mass().value();
            liquid_ideal_cv += x[index] * (cp - R);
            vapor_ideal_cv += y[index] * (cp - R);
        }
        // Phase mechanics and heat capacities use the SAME phase composition;
        // the feed Cv is not a legitimate substitute for either separated phase.
        if phase
            .liquid()
            .caloric_mechanics(liquid_ideal_cv, liquid_mass)
            .is_none()
            || phase
                .vapor()
                .caloric_mechanics(vapor_ideal_cv, vapor_mass)
                .is_none()
        {
            return Err(E::Invalid(
                "both original flash phases require stable caloric mechanics",
            ));
        }
        let liquid_mass_fraction = (1.0 - beta) * (liquid_mass / mass);
        let vapor_mass_fraction = beta * (vapor_mass / mass);
        if (liquid_mass_fraction + vapor_mass_fraction - 1.0).abs()
            > 4096.0 * f64::EPSILON * x.len() as f64
        {
            return Err(E::Invalid(
                "flash phase masses do not close the retained feed mass",
            ));
        }
        let molar_volume = (1.0 - beta) * phase.liquid().molar_volume_m3_per_mol()
            + beta * phase.vapor().molar_volume_m3_per_mol();
        let vapor_volume_fraction = beta * (phase.vapor().molar_volume_m3_per_mol() / molar_volume);
        let residual_h = (1.0 - beta) * phase.liquid().residual_enthalpy_j_mol()
            + beta * phase.vapor().residual_enthalpy_j_mol();
        let residual_u = (1.0 - beta) * phase.liquid().residual_internal_energy_j_mol()
            + beta * phase.vapor().residual_internal_energy_j_mol();
        let split_mixing = R
            * (x_log_x(self.mechanics().mole_fractions())
                - (1.0 - beta) * x_log_x(x)
                - beta * x_log_x(y));
        let entropy_correction = split_mixing
            + (1.0 - beta) * residual_entropy(phase.liquid(), x, temperature.value())
            + beta * residual_entropy(phase.vapor(), y, temperature.value());
        // The exact feed inventory carries every species reference once. Never
        // form this by subtracting/recombining large absolute phase enthalpies.
        let specific_h = (reference.h().value() + residual_h) / mass;
        let specific_u = (reference.u().value() + residual_u) / mass;
        let specific_s = (reference.s().value() + entropy_correction) / mass;
        if [
            residual_h,
            residual_u,
            entropy_correction,
            specific_h,
            specific_u,
            specific_s,
        ]
        .into_iter()
        .any(|value| !value.is_finite())
            || !positive(molar_volume)
            || !positive(mass / molar_volume)
            || [
                liquid_mass_fraction,
                vapor_mass_fraction,
                vapor_volume_fraction,
            ]
            .into_iter()
            .any(|fraction| !positive(fraction) || fraction >= 1.0)
        {
            return Err(E::Invalid(
                "source-bound two-phase calorics are unrepresentable",
            ));
        }
        poll(cancelled)?;
        Ok(CaloricSample {
            phase: (),
            temperature,
            pressure,
            eos_evaluation_upper_bound: 0,
            molar_volume,
            liquid_mass_fraction,
            vapor_mass_fraction,
            vapor_volume_fraction,
            residual_h,
            residual_u,
            entropy_correction,
            specific_h,
            specific_u,
            specific_s,
        })
    }

    fn publish_caloric_flash(
        &self,
        sample: CaloricSample,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ResolvedThermochemicalFlashState, E> {
        poll(cancelled)?;
        let reference = self
            .source()
            .thermodynamics()
            .evaluate(sample.phase.temperature(), sample.phase.pressure())
            .map_err(|error| failure("accepted flash source receipt", error))?;
        let state = ResolvedThermochemicalFlashState {
            law: self.mechanics().clone(),
            reference,
            sample,
        };
        poll(cancelled)?;
        Ok(state)
    }
}
