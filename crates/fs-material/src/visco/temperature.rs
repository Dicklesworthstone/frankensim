//! Temperature-dependent viscous kinetics for thermorheologically simple solids.
//!
//! The convention is `a_T = tau(T)/tau(T_ref)`: a positive activation energy
//! makes relaxation faster on heating. At a fixed temperature,
//! `C_T*(omega) = C_ref*(a_T omega)` and `C_T(t) = C_ref(t/a_T)`.
//! Consequently BOTH Maxwell times and instantaneous Kelvin viscosity scale
//! by `a_T`. Elastic tensors and their admitted factors stay unchanged.
//!
//! References: Williams, Landel & Ferry (1955), doi:10.1021/ja01619a008;
//! [IUPAC WLF definition](https://doi.org/10.1351/goldbook.12821);
//! [Arrhenius/reduced-time convention](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEMATRefMap/simamat-c-timevisco.htm).
//! Coefficients, reference temperature and calibrated validity are caller
//! declarations: no material-name constants, glass transition, source accuracy,
//! moisture shift, vertical modulus shift or thermal expansion is inferred.
//! This supplies fixed-temperature laws, not integration along a changing
//! thermal history or a heat-capacity/entropy model.

use super::{TensorProny, ViscoError};
use fs_math::det;

mod scalar;
pub use scalar::*;

use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as GAS_CONSTANT_J_MOL_K;

/// Explicit empirical temperature dependence of the relaxation times.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TemperatureShiftModel {
    /// `log10(a_T) = -C1 (T-T_ref)/(C2+T-T_ref)`.
    WilliamsLandelFerry {
        /// Positive dimensionless coefficient of the base-ten law.
        c1: f64,
        /// Positive temperature interval [K]; the pole is `T_ref-C2`.
        c2_k: f64,
    },
    /// `ln(a_T) = (E/R) (1/T-1/T_ref)` with absolute kelvin temperatures.
    Arrhenius {
        /// Positive molar activation energy [J/mol], never a per-particle energy.
        activation_energy_j_mol: f64,
    },
}

/// One declared calibration and its complete inclusive absolute-temperature band.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemperatureShift {
    model: TemperatureShiftModel,
    reference_temperature_k: f64,
    validity_k: [f64; 2],
}

/// Checked shift and analytic logarithmic derivative at one temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemperatureShiftState {
    declaration: TemperatureShift,
    temperature_k: f64,
    factor: f64,
    log_factor: f64,
    log_derivative_per_k: f64,
}

/// The existing tensor constitutive law with the exact temperature declaration.
/// Its stiffness, branch order, anisotropy and strain bound remain the original
/// law's. The source frequency band is transported with the same shift.
pub type ShiftedTensorProny = ShiftedViscoelastic<TensorProny>;

fn invalid(what: &'static str) -> ViscoError {
    ViscoError::Parameters { what }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

impl TemperatureShift {
    /// Admit the caller's calibration, with reference temperature inside its
    /// positive finite validity interval. All WLF support must lie strictly
    /// above its pole; both models must remain representable at both ends.
    /// No universal WLF constants or automatic extrapolation are supplied.
    ///
    /// # Errors
    /// Refuses invalid parameters, a WLF pole inside support, or a shift or
    /// derivative that overflows or vanishes by underflow over that support.
    pub fn new(
        model: TemperatureShiftModel,
        reference_temperature_k: f64,
        validity_k: [f64; 2],
    ) -> Result<Self, ViscoError> {
        if !positive(reference_temperature_k)
            || !validity_k.into_iter().all(positive)
            || validity_k[0] > validity_k[1]
            || !(validity_k[0]..=validity_k[1]).contains(&reference_temperature_k)
        {
            return Err(invalid(
                "temperature shift requires a positive finite band containing its reference temperature",
            ));
        }
        match model {
            TemperatureShiftModel::WilliamsLandelFerry { c1, c2_k } => {
                if !positive(c1)
                    || !positive(c2_k)
                    || !positive(c2_k + (validity_k[0] - reference_temperature_k))
                {
                    return Err(invalid(
                        "WLF requires positive C1/C2 and a validity band strictly above its pole",
                    ));
                }
            }
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol,
            } => {
                if !positive(activation_energy_j_mol) {
                    return Err(invalid(
                        "Arrhenius relaxation requires positive finite molar activation energy",
                    ));
                }
            }
        }
        let shift = Self {
            model,
            reference_temperature_k,
            validity_k,
        };
        // Both laws and their negative logarithmic slopes are monotone on this
        // admitted domain; extrema of their magnitude occur at the endpoints.
        shift.evaluate(validity_k[0])?;
        shift.evaluate(validity_k[1])?;
        Ok(shift)
    }

    /// Original empirical law, with the caller's exact coefficients.
    pub const fn model(&self) -> TemperatureShiftModel {
        self.model
    }
    /// Temperature [K] at which the original relaxation spectrum was declared.
    pub const fn reference_temperature_k(&self) -> f64 {
        self.reference_temperature_k
    }
    /// Original inclusive calibrated support [K].
    pub const fn validity_k(&self) -> [f64; 2] {
        self.validity_k
    }

    /// Evaluate the SAME shift convention and its analytic `d ln(a_T)/dT`.
    /// The reference is an exact fixed point, `a_T(T_ref) = 1`.
    ///
    /// # Errors
    /// Refuses temperatures outside support and unrepresentable responses.
    pub fn evaluate(&self, temperature_k: f64) -> Result<TemperatureShiftState, ViscoError> {
        if !temperature_k.is_finite()
            || !(self.validity_k[0]..=self.validity_k[1]).contains(&temperature_k)
        {
            return Err(invalid(
                "temperature leaves the declared shift calibration band",
            ));
        }
        let delta = temperature_k - self.reference_temperature_k;
        let (log_factor, log_derivative_per_k) = match self.model {
            TemperatureShiftModel::WilliamsLandelFerry { c1, c2_k } => {
                let denominator = c2_k + delta;
                if !positive(denominator) {
                    return Err(invalid(
                        "WLF denominator must remain finite and strictly positive",
                    ));
                }
                (
                    -core::f64::consts::LN_10 * (c1 * (delta / denominator)),
                    -core::f64::consts::LN_10 * (c1 / denominator) * (c2_k / denominator),
                )
            }
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol,
            } => {
                let activation_k = activation_energy_j_mol / GAS_CONSTANT_J_MOL_K;
                (
                    -activation_k * ((delta / temperature_k) / self.reference_temperature_k),
                    -(activation_k / temperature_k) / temperature_k,
                )
            }
        };
        let (log_factor, factor) = if temperature_k == self.reference_temperature_k {
            (0.0, 1.0)
        } else {
            (log_factor, det::exp(log_factor))
        };
        if !log_factor.is_finite()
            || !positive(factor)
            || !log_derivative_per_k.is_finite()
            || log_derivative_per_k >= 0.0
        {
            return Err(invalid(
                "temperature shift or logarithmic slope is unrepresentable",
            ));
        }
        Ok(TemperatureShiftState {
            declaration: *self,
            temperature_k,
            factor,
            log_factor,
            log_derivative_per_k,
        })
    }

    /// Shift the original reference-temperature tensor law without refitting,
    /// diagonalizing, dropping branches or changing stored elastic energy.
    /// Use the original reference law for every call, not a previous result.
    ///
    /// This transports the frequency support to `[omega_lo/a_T, omega_hi/a_T]`,
    /// intersected with finite representable frequencies. It is not authority
    /// to apply time-temperature superposition to an arbitrary material: the
    /// caller must qualify the calibration against that same original source.
    /// Work is O(branch count); no tensor eigensolve or factorization is repeated.
    ///
    /// # Errors
    /// Refuses invalid temperature, unrepresentable positive relaxation times,
    /// lost nonzero viscosity entries or an unrepresentable frequency band.
    /// The original law is never modified on success or refusal.
    pub fn shift_tensor(
        &self,
        reference: &TensorProny,
        temperature_k: f64,
    ) -> Result<ShiftedTensorProny, ViscoError> {
        let state = self.evaluate(temperature_k)?;
        let mut law = reference.clone();
        if state.factor == 1.0 {
            return Ok(ShiftedTensorProny { law, state });
        }
        for branch in &mut law.branches {
            branch.tau_s = scaled(branch.tau_s, state.factor)?;
        }
        if let Some(viscosity) = &mut law.kv_pa_s {
            for value in viscosity.iter_mut().flatten() {
                *value = scaled(*value, state.factor)?;
            }
        }
        if let Some(factor) = &mut law.kv_factor {
            let root = det::sqrt(state.factor);
            for value in factor.iter_mut().flatten() {
                *value = scaled(*value, root)?;
            }
        }
        let low = reference.band_rad_s.0 / state.factor;
        let high = (reference.band_rad_s.1 / state.factor).min(f64::MAX);
        if !low.is_finite()
            || low < 0.0
            || !high.is_finite()
            || low > high
            || reference.band_rad_s.0 > 0.0 && low == 0.0
            || reference.band_rad_s.1 > 0.0 && high == 0.0
        {
            return Err(invalid(
                "shifted source frequency support is unrepresentable",
            ));
        }
        law.band_rad_s = (low, high);
        Ok(ShiftedTensorProny { law, state })
    }
}

// Preserve exact zeros and signs, including legitimate negative off-diagonal
// tensor entries. A nonzero coupling cannot silently vanish through underflow.
fn scaled(value: f64, factor: f64) -> Result<f64, ViscoError> {
    let next = value * factor;
    if !next.is_finite() || value != 0.0 && next == 0.0 {
        return Err(invalid(
            "shifted relaxation time or viscosity is unrepresentable",
        ));
    }
    Ok(next)
}

impl TemperatureShiftState {
    /// Exact calibration declaration used for this evaluation.
    pub const fn declaration(&self) -> &TemperatureShift {
        &self.declaration
    }
    /// Actual absolute evaluation temperature [K].
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    /// Positive relaxation-time ratio `tau(T)/tau(T_ref)`.
    pub const fn factor(&self) -> f64 {
        self.factor
    }
    /// Natural logarithm of that ratio; WLF base-ten data is converted once.
    pub const fn log_factor(&self) -> f64 {
        self.log_factor
    }
    /// Analytic `d ln(a_T)/dT` [1/K], strictly negative for the admitted models.
    pub const fn log_derivative_per_k(&self) -> f64 {
        self.log_derivative_per_k
    }
}
