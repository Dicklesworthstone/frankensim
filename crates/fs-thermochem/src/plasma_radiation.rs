//! Electron–ion thermal free-free emissive power from original plasma species.
//!
//! Itoh, Kusano, Sakamoto & Nozawa, arXiv:astro-ph/0011044, equation (11),
//! gives the frequency-integrated nonrelativistic emissivity
//! `1.426e-27 * g * n_e * sum(n_i Z_i²) * sqrt(Te)` in erg/(s cm³),
//! with number densities in cm⁻³ and Te in kelvin. The SI coefficient is
//! therefore `1.426e-40`, not the eV-temperature NRL coefficient with its
//! implicit Gaunt approximation. The caller selects a constant integrated
//! Gaunt factor or the paper's bounded nonrelativistic fit, equations (10),
//! (13)–(14), table 2. The latter evaluates each populated positive ion's
//! `g(Z²*1.579e5 K/Te)` separately, including its analytic temperature slope.
//! For mixed ions the sum is `sum(n_i Z_i² g_i)`. The source's stated fit
//! residual is not a physical uncertainty or derivative-error certificate.
//!
//! This is a nondegenerate Maxwellian electron population scattering from
//! effectively stationary Coulomb ions. Actual electron degeneracy, electron
//! relativity and ion/electron thermal-variance ratios must pass caller-selected
//! limits. Weak coupling, negligible collective spectral effects and validity
//! of a bare-Coulomb-ion Gaunt approximation remain the plasma owner's scope.
//! Source mass, charge and caloric support are retained; ionization state is
//! never inferred from a name, and helium ash contributes Z²=4.
//!
//! The result is emitted power. Equating it with escaped power requires an
//! explicit optically thin closure from the consumer; no opacity, absorption,
//! radiation field, electron–electron, recombination, line or synchrotron
//! emission model is supplied here. This L1 module evolves no state.

use crate::electron_collision::{AVOGADRO_PER_MOL, BOLTZMANN_J_PER_K, ELECTRON_MASS_KG};
use crate::equilibrium::composition::{FrozenIonizedGasStateV1, IdealGasEquilibriumV1};
use crate::{SpeciesId, ThermochemErrorV1};
use fs_math::det::{exp, ln, sqrt};
use fs_qty::{Amount, Temperature, Volume};

#[cfg(test)]
mod tests;

/// Primary normalization and definition of the integrated Maxwellian Gaunt factor.
pub const BREMSSTRAHLUNG_SOURCE_URL_V1: &str = "https://arxiv.org/abs/astro-ph/0011044";
/// Published coefficient converted to SI: W m³ K⁻¹ᐟ², before the Gaunt factor.
/// The four-digit source normalization is retained, not an exact SI constant.
pub const ELECTRON_ION_BREMSSTRAHLUNG_COEFFICIENT_SI_V1: f64 = 1.426e-40;

const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;
const PLANCK_J_S: f64 = 6.626_070_15e-34;

/// Explicit frequency-integrated Maxwellian Gaunt law; no automatic fallback.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BremsstrahlungGauntFactorV1 {
    /// Caller-authored positive finite factor, shared by all positive ions.
    Constant(f64),
    /// Itoh et al., equations (10), (13)–(14), table 2: degree 10 in
    /// `(log10(gamma²)+0.5)/2.5`, with `gamma²=Z²*1.579e5 K/Te`.
    /// Every populated ion must satisfy `1e-3 <= gamma² <= 100`.
    /// Thus D/T require 1,579–157,900,000 K; He++ requires 6,316–631,600,000 K,
    /// intersecting the original source and the other physical limits.
    ItohNonrelativistic,
}

// Literal eleven coefficients from table 2 of arXiv:astro-ph/0011044v2.
// Retain the paper's rounded Rydberg-temperature scale with its source fit.
const ITOH_RYDBERG_TEMPERATURE_K: f64 = 1.579e5;
const ITOH_NR_GAUNT_COEFFICIENTS: [f64; 11] = [
    1.42752, 0.170677, -0.813616, -0.488059, 1.33141, 0.727788, -1.50917, -0.599317, 0.986146,
    0.199586, -0.272392,
];

/// Return (g, dg/dlnTe); charge comes from the retained positive-ion axis.
fn itoh_nonrelativistic_gaunt(
    charge: f64,
    temperature_k: f64,
) -> Result<(f64, f64), BremsstrahlungErrorV1> {
    let scale = (charge * charge) * ITOH_RYDBERG_TEMPERATURE_K;
    // Compare in the authored Kelvin basis. Converting an exact source
    // endpoint to gamma² can otherwise round it an ULP outside the interval.
    if temperature_k < scale / 100.0 || temperature_k > scale / 1.0e-3 {
        return Err(BremsstrahlungErrorV1::Domain(
            "outside the populated ion's Itoh nonrelativistic Gaunt fit",
        ));
    }
    let x = (ln(scale / temperature_k) / core::f64::consts::LN_10 + 0.5) / 2.5;
    let mut value = ITOH_NR_GAUNT_COEFFICIENTS[10];
    let mut derivative_x = 0.0_f64;
    for &coefficient in ITOH_NR_GAUNT_COEFFICIENTS[..10].iter().rev() {
        derivative_x = derivative_x.mul_add(x, value);
        value = value.mul_add(x, coefficient);
    }
    let derivative_log_temperature = -derivative_x / (2.5 * core::f64::consts::LN_10);
    if !positive(value) || !derivative_log_temperature.is_finite() {
        return Err(BremsstrahlungErrorV1::Arithmetic(
            "Itoh Gaunt fit or analytic slope",
        ));
    }
    Ok((value, derivative_log_temperature))
}

/// Explicit Gaunt approximation and physical validity screens, with no defaults.
/// These limits are admission policies, not accuracy or uncertainty certificates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BremsstrahlungControlsV1 {
    /// Selected frequency-integrated Gaunt law, retained as part of model identity.
    pub frequency_integrated_gaunt_factor: BremsstrahlungGauntFactorV1,
    /// Inclusive caller-selected lower temperature bound [K], intersecting any
    /// per-ion published fit domain and the original material source support.
    pub minimum_electron_temperature_k: f64,
    /// Inclusive caller-selected upper temperature bound [K].
    pub maximum_electron_temperature_k: f64,
    /// Largest admitted kB*Te/(me*c²), strictly inside (0,1).
    pub maximum_electron_thermal_rest_energy_ratio: f64,
    /// Largest admitted (Th/Te)*(me/mi) among populated positive ions, strictly
    /// inside (0,1). This bounds the neglected ion contribution to relative
    /// thermal velocity variance; it is not a finite-ion-mass correction.
    pub maximum_ion_to_electron_thermal_variance_ratio: f64,
    /// Largest admitted ne*[h/sqrt(2*pi*me*kB*Te)]³, strictly inside (0,1).
    /// No spin degeneracy factor is inferred.
    pub maximum_electron_thermal_phase_space_density: f64,
}

impl BremsstrahlungControlsV1 {
    /// Reject invalid Gaunt data, temperature intervals and physical limits.
    pub fn validate(self) -> Result<(), BremsstrahlungErrorV1> {
        if matches!(self.frequency_integrated_gaunt_factor, BremsstrahlungGauntFactorV1::Constant(g) if !positive(g))
            || !positive(self.minimum_electron_temperature_k)
            || !positive(self.maximum_electron_temperature_k)
            || self.minimum_electron_temperature_k > self.maximum_electron_temperature_k
            || [
                self.maximum_electron_thermal_rest_energy_ratio,
                self.maximum_ion_to_electron_thermal_variance_ratio,
                self.maximum_electron_thermal_phase_space_density,
            ]
            .iter()
            .any(|&x| !positive(x) || x >= 1.0)
        {
            return Err(BremsstrahlungErrorV1::Invalid(
                "positive Gaunt factor, ordered positive temperature interval and physical maxima in (0,1)",
            ));
        }
        Ok(())
    }
}

/// Source, physical-domain, arithmetic or cancellation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum BremsstrahlungErrorV1 {
    /// Invalid bank, source binding, population or authored controls.
    Invalid(&'static str),
    /// The temperature approximation or an explicit physical limit is exceeded.
    Domain(&'static str),
    /// The original populated source refuses its own actual temperature.
    Species {
        /// Unmodified original species identity.
        species: SpeciesId,
        /// Original source error.
        source: ThermochemErrorV1,
    },
    /// A required density or positive emission result is not representable.
    Arithmetic(&'static str),
    /// Cancellation before publication of a complete result.
    Cancelled,
}
impl core::fmt::Display for BremsstrahlungErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "electron-ion bremsstrahlung: {self:?}")
    }
}
impl std::error::Error for BremsstrahlungErrorV1 {}

/// Immutable original species bank and explicitly selected Gaunt radiation law.
/// The bank's existing bounded canonical charge axis is the only species map.
/// Successful evaluations allocate nothing and scan at most128 source species.
#[derive(Clone, Debug, PartialEq)]
pub struct ElectronIonBremsstrahlungV1 {
    thermodynamics: IdealGasEquilibriumV1,
    electron_index: usize,
    controls: BremsstrahlungControlsV1,
}

impl ElectronIonBremsstrahlungV1 {
    /// Retain the original bank. One nucleus-free charge −1 electron with the
    /// physical electron mass is required; negative heavy ions are unsupported.
    /// Neutral rows may coexist, but do not contribute to this free-free term.
    pub fn new(
        thermodynamics: IdealGasEquilibriumV1,
        controls: BremsstrahlungControlsV1,
    ) -> Result<Self, BremsstrahlungErrorV1> {
        use BremsstrahlungErrorV1 as E;
        controls.validate()?;
        let elements = thermodynamics.elemental_matrix();
        let electron_index = (0..thermodynamics.models().len())
            .find(|&i| {
                thermodynamics.charge_vector().get(i) == Some(-1)
                    && (0..elements.row_count()).all(|row| elements.get(row, i) == Some(0))
            })
            .ok_or(E::Invalid("original bank requires a structural electron"))?;
        let electron_mass =
            thermodynamics.models()[electron_index].molar_mass().value() / AVOGADRO_PER_MOL;
        if ((electron_mass - ELECTRON_MASS_KG) / ELECTRON_MASS_KG).abs() > 64.0 * f64::EPSILON {
            return Err(E::Invalid(
                "original structural-electron mass is not physical",
            ));
        }
        if (0..thermodynamics.models().len())
            .any(|i| i != electron_index && thermodynamics.charge_vector().get(i).unwrap_or(-1) < 0)
        {
            return Err(E::Invalid(
                "negative heavy ions are outside this electron-positive-ion law",
            ));
        }
        Ok(Self {
            thermodynamics,
            electron_index,
            controls,
        })
    }

    /// Complete unmodified material source bank and canonical species axes.
    pub const fn thermodynamics(&self) -> &IdealGasEquilibriumV1 {
        &self.thermodynamics
    }
    /// Canonical original structural-electron index.
    pub const fn electron_index(&self) -> usize {
        self.electron_index
    }
    /// Retained explicit Gaunt approximation and physical admission controls.
    pub const fn controls(&self) -> BremsstrahlungControlsV1 {
        self.controls
    }

    /// Evaluate an original frozen state; changing its source bank is refused.
    pub fn evaluate(
        &self,
        state: &FrozenIonizedGasStateV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<BremsstrahlungEvaluationV1, BremsstrahlungErrorV1> {
        if state.thermodynamics() != &self.thermodynamics
            || state.electron_index() != Some(self.electron_index)
        {
            return Err(BremsstrahlungErrorV1::Invalid(
                "state belongs to another source bank",
            ));
        }
        self.evaluate_point(
            state.amounts(),
            state.volume(),
            state.heavy_temperature(),
            state.electron_temperature(),
            cancelled,
        )
    }

    /// Allocation-free evaluation at evolving original-bank populations.
    /// Both temperatures are needed to admit effectively stationary ions.
    /// Original populated caloric sources must admit their respective Th/Te.
    /// Charge neutrality and weak coupling remain the consuming state owner's
    /// responsibility; this local emissivity does not solve electrostatic fields.
    pub fn evaluate_point(
        &self,
        amounts: &[Amount],
        volume: Volume,
        heavy_temperature: Temperature,
        electron_temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BremsstrahlungEvaluationV1, BremsstrahlungErrorV1> {
        use BremsstrahlungErrorV1 as E;
        poll(&mut cancelled)?;
        let (th, te, v) = (
            heavy_temperature.value(),
            electron_temperature.value(),
            volume.value(),
        );
        if amounts.len() != self.thermodynamics.models().len()
            || !positive(v)
            || !positive(th)
            || !positive(te)
            || amounts
                .iter()
                .any(|n| !n.value().is_finite() || n.value() < 0.0)
        {
            return Err(E::Invalid(
                "canonical nonnegative amounts and positive finite V, Th, Te required",
            ));
        }
        if te < self.controls.minimum_electron_temperature_k
            || te > self.controls.maximum_electron_temperature_k
        {
            return Err(E::Domain("outside declared radiation temperature interval"));
        }
        let ne = number_density(amounts[self.electron_index].value(), v)?;
        let mut ion_weight = Sum::default();
        let mut fitted_ion_weight = Sum::default();
        let mut fitted_ion_log_derivative = Sum::default();
        let mut maximum_ion_variance = 0.0_f64;
        for (index, (model, amount)) in self.thermodynamics.models().iter().zip(amounts).enumerate()
        {
            poll(&mut cancelled)?;
            if amount.value() == 0.0 {
                continue;
            }
            let temperature = if index == self.electron_index {
                electron_temperature
            } else {
                heavy_temperature
            };
            model
                .evaluate_properties(temperature)
                .map_err(|source| E::Species {
                    species: model.species().clone(),
                    source,
                })?;
            let z = self
                .thermodynamics
                .charge_vector()
                .get(index)
                .expect("retained charge axis");
            if z > 0 {
                let mass = model.molar_mass().value() / AVOGADRO_PER_MOL;
                let variance = (th / te) * (ELECTRON_MASS_KG / mass);
                if !variance.is_finite()
                    || variance > self.controls.maximum_ion_to_electron_thermal_variance_ratio
                {
                    return Err(E::Domain(
                        "ion thermal motion exceeds the stationary-ion limit",
                    ));
                }
                maximum_ion_variance = maximum_ion_variance.max(variance);
                let z = z as f64;
                let weight = (amount.value() * z) * z;
                if !positive(weight) {
                    return Err(E::Arithmetic("charge-squared ion amount"));
                }
                ion_weight.add(weight);
                if self.controls.frequency_integrated_gaunt_factor
                    == BremsstrahlungGauntFactorV1::ItohNonrelativistic
                {
                    let (gaunt, derivative_log_temperature) = itoh_nonrelativistic_gaunt(z, te)?;
                    fitted_ion_weight.add(weight * gaunt);
                    fitted_ion_log_derivative.add(weight * derivative_log_temperature);
                }
            }
        }
        if !ion_weight.value.is_finite() {
            return Err(E::Arithmetic("summed charge-squared ion amount"));
        }
        let ions = number_density(ion_weight.value, v)?;
        let rest_ratio =
            BOLTZMANN_J_PER_K * te / (ELECTRON_MASS_KG * SPEED_OF_LIGHT_M_S * SPEED_OF_LIGHT_M_S);
        let wavelength = PLANCK_J_S
            / sqrt((2.0 * core::f64::consts::PI * ELECTRON_MASS_KG) * (BOLTZMANN_J_PER_K * te));
        let phase_space = ((ne * wavelength) * wavelength) * wavelength;
        if ne > 0.0
            && (!rest_ratio.is_finite()
                || rest_ratio > self.controls.maximum_electron_thermal_rest_energy_ratio)
        {
            return Err(E::Domain(
                "electron temperature exceeds the nonrelativistic limit",
            ));
        }
        if ne > 0.0
            && (!phase_space.is_finite()
                || phase_space > self.controls.maximum_electron_thermal_phase_space_density)
        {
            return Err(E::Domain(
                "electrons exceed the nondegenerate phase-space limit",
            ));
        }
        let mut result = BremsstrahlungEvaluationV1 {
            total_power: 0.0,
            power_density: 0.0,
            temperature_coefficient: 0.0,
            derivative: 0.0,
            electron_density: ne,
            charge_squared_ion_density: ions,
            maximum_ion_variance,
            electron_phase_space: if ne == 0.0 { 0.0 } else { phase_space },
        };
        if ne > 0.0 && ions > 0.0 {
            let (density_coefficient, logarithmic_gaunt_slope) =
                match self.controls.frequency_integrated_gaunt_factor {
                    BremsstrahlungGauntFactorV1::Constant(gaunt) => (
                        product(&[
                            ELECTRON_ION_BREMSSTRAHLUNG_COEFFICIENT_SI_V1,
                            gaunt,
                            ne,
                            ions,
                        ])?,
                        0.0,
                    ),
                    BremsstrahlungGauntFactorV1::ItohNonrelativistic => {
                        if !positive(fitted_ion_weight.value)
                            || !fitted_ion_log_derivative.value.is_finite()
                        {
                            return Err(E::Arithmetic("summed fitted Gaunt contributions"));
                        }
                        (
                            product(&[
                                ELECTRON_ION_BREMSSTRAHLUNG_COEFFICIENT_SI_V1,
                                ne,
                                number_density(fitted_ion_weight.value, v)?,
                            ])?,
                            fitted_ion_log_derivative.value / fitted_ion_weight.value,
                        )
                    }
                };
            result.temperature_coefficient = product(&[density_coefficient, v])?;
            result.power_density = product(&[density_coefficient, sqrt(te)])?;
            result.total_power = product(&[result.temperature_coefficient, sqrt(te)])?;
            let logarithmic_power_slope = 0.5 + logarithmic_gaunt_slope;
            result.derivative = if logarithmic_power_slope == 0.0 {
                0.0
            } else {
                logarithmic_power_slope.signum()
                    * product(&[logarithmic_power_slope.abs(), result.total_power, 1.0 / te])?
            };
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
}

/// Emission and actual admitted population scales for one trial point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BremsstrahlungEvaluationV1 {
    total_power: f64,
    power_density: f64,
    temperature_coefficient: f64,
    derivative: f64,
    electron_density: f64,
    charge_squared_ion_density: f64,
    maximum_ion_variance: f64,
    electron_phase_space: f64,
}
impl BremsstrahlungEvaluationV1 {
    /// Frequency-integrated electron–ion emissive power [W], before transport.
    pub const fn total_power_w(self) -> f64 {
        self.total_power
    }
    /// Frequency-integrated local emissive power density [W/m³].
    pub const fn power_density_w_m3(self) -> f64 {
        self.power_density
    }
    /// Instantaneous P/sqrt(Te) [W/K¹ᐟ²]. It is temperature-independent only
    /// for the Constant Gaunt law at frozen populations. A fitted-law consumer
    /// may explicitly freeze this coefficient at a predictor midpoint; that
    /// thermal substep is a numerical approximation to the changing Gaunt law.
    /// Subsequent temperatures must still pass every source/domain check.
    pub const fn total_power_over_sqrt_temperature_w_k_sqrt(self) -> f64 {
        self.temperature_coefficient
    }
    /// Partial dP/dTe at fixed populations and volume [W/K]. This includes the
    /// selected fitted Gaunt law's analytic derivative when present.
    pub const fn derivative_w_per_k(self) -> f64 {
        self.derivative
    }
    /// Original structural-electron number density [m⁻³].
    pub const fn electron_density_per_m3(self) -> f64 {
        self.electron_density
    }
    /// Sum of positive-ion number densities weighted by charge squared [m⁻³].
    pub const fn charge_squared_ion_density_per_m3(self) -> f64 {
        self.charge_squared_ion_density
    }
    /// Largest populated-ion (Th/Te)*(me/mi), against the caller's limit.
    pub const fn maximum_ion_to_electron_thermal_variance_ratio(self) -> f64 {
        self.maximum_ion_variance
    }
    /// Actual ne*lambda_e³, without an inferred spin factor.
    pub const fn electron_thermal_phase_space_density(self) -> f64 {
        self.electron_phase_space
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
        let next = self.value + y;
        self.correction = (next - self.value) - y;
        self.value = next;
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), BremsstrahlungErrorV1> {
    if cancelled() {
        Err(BremsstrahlungErrorV1::Cancelled)
    } else {
        Ok(())
    }
}
fn number_density(amount: f64, volume: f64) -> Result<f64, BremsstrahlungErrorV1> {
    if amount == 0.0 {
        return Ok(0.0);
    }
    let direct = (amount / volume) * AVOGADRO_PER_MOL;
    if positive(direct) {
        return Ok(direct);
    }
    let density = exp(ln(amount) + ln(AVOGADRO_PER_MOL) - ln(volume));
    if positive(density) {
        Ok(density)
    } else {
        Err(BremsstrahlungErrorV1::Arithmetic("particle number density"))
    }
}
fn product(factors: &[f64]) -> Result<f64, BremsstrahlungErrorV1> {
    let direct = factors.iter().product::<f64>();
    if positive(direct) {
        return Ok(direct);
    }
    // Normal states take the inexpensive arithmetic path. A log fallback
    // avoids rejecting a representable result because an intermediate product
    // overflowed/underflowed; an unrepresentable positive result is never zero.
    let value = exp(factors.iter().map(|&x| ln(x)).sum());
    if positive(value) {
        Ok(value)
    } else {
        Err(BremsstrahlungErrorV1::Arithmetic(
            "positive emission product",
        ))
    }
}
