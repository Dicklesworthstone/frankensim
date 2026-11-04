//! Local thermodynamics from a fluid's original `rho(T, p)` derivatives.
//!
//! For a homogeneous, simple compressible phase at fixed composition,
//! `beta = -rho_T/rho`, `kappa_T = rho_p/rho` and the Gibbs identities give
//! `Cp - Cv = T beta² / rho_p`, `kappa_s = kappa_T Cv/Cp` and
//! `c² = 1/(rho kappa_s)`. In particular, an isothermal density slope alone
//! does not give the adiabatic sound speed when thermal expansion is nonzero.
//!
//! Derivation: `ds = (Cp/T) dT - (beta/rho) dp`; substituting its isentropic
//! temperature response into `d rho = rho_T dT + rho_p dp` yields
//! `(d rho/dp)_s = rho_p - T beta²/Cp`. These are local equilibrium identities,
//! not an ideal-gas approximation. The original owner supplies consistent Cp,
//! density and both derivatives at the SAME state, on its selected phase.
//!
//! This fixed-size, allocation-free algebra neither selects sources nor fits
//! a new EOS. It supplies no phase-stability/critical-point extrapolation,
//! derivative uncertainty, or proof that independently supplied Cp and density
//! functions integrate to a common thermodynamic potential over a finite path.

use core::fmt;
use fs_qty::{Density, Qty, Temperature};

use crate::MassSpecificThermalQuantityV1;

/// Isobaric mass-density derivative [kg/(m³ K)].
pub type DensityTemperatureDerivativeV1 = Qty<-3, 1, 0, -1, 0, 0>;
/// Isothermal mass-density derivative [kg/(m³ Pa)] = [s²/m²].
pub type DensityPressureDerivativeV1 = Qty<-2, 0, 2, 0, 0, 0>;

/// Invalid local inputs, unstable response or unrepresentable arithmetic.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DensityResponseErrorV1 {
    Invalid(&'static str),
    Unstable(&'static str),
    Arithmetic(&'static str),
}

impl fmt::Display for DensityResponseErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "density thermodynamic response: {self:?}")
    }
}
impl std::error::Error for DensityResponseErrorV1 {}

/// One immutable density, thermal-expansion, caloric and acoustic response.
/// Numerical inputs confer no measured material or source authority.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensityThermodynamicResponseV1 {
    temperature_k: f64,
    density_kg_m3: f64,
    specific_heat_cp_j_kg_k: f64,
    density_temperature_derivative: f64,
    density_pressure_derivative: f64,
    volumetric_expansion_per_k: f64,
    isothermal_compressibility_per_pa: f64,
    adiabatic_compressibility_per_pa: f64,
    specific_heat_cv_j_kg_k: f64,
    heat_capacity_difference_j_kg_k: f64,
    heat_capacity_ratio: f64,
    adiabatic_bulk_modulus_pa: f64,
    sound_speed_m_s: f64,
}

impl DensityThermodynamicResponseV1 {
    /// Derive the local response from explicitly isobaric Cp, `rho_T|p` and
    /// `rho_p|T`. Negative or zero thermal expansion is retained. Mechanical
    /// stability requires `rho_p > 0`; caloric stability requires `Cv > 0`.
    /// The exact incompressible and critical limits have no finite acoustic
    /// response here and refuse instead of acquiring an arbitrary sound speed.
    /// A residual Cv no larger than 16 rounding units of Cp is unresolved at
    /// the stability boundary and also refuses; it cannot establish a finite
    /// positive acoustic response from these rounded inputs.
    pub fn try_new(
        temperature: Temperature,
        density: Density,
        specific_heat_cp: MassSpecificThermalQuantityV1,
        density_temperature_derivative: DensityTemperatureDerivativeV1,
        density_pressure_derivative: DensityPressureDerivativeV1,
    ) -> Result<Self, DensityResponseErrorV1> {
        let (t, rho, cp) = (
            temperature.value(),
            density.value(),
            specific_heat_cp.value(),
        );
        let rho_t = density_temperature_derivative.value();
        let rho_p = density_pressure_derivative.value();
        if ![t, rho, cp].into_iter().all(positive) || !rho_t.is_finite() || !rho_p.is_finite() {
            return Err(DensityResponseErrorV1::Invalid(
                "temperature, density and Cp must be positive finite; density derivatives finite",
            ));
        }
        if rho_p <= 0.0 {
            return Err(DensityResponseErrorV1::Unstable(
                "positive finite isothermal density response is required",
            ));
        }
        let beta = -rho_t / rho;
        let kappa_t = rho_p / rho;
        if !beta.is_finite() || rho_t != 0.0 && beta == 0.0 || !positive(kappa_t) {
            return Err(DensityResponseErrorV1::Arithmetic(
                "density expansion or isothermal compressibility",
            ));
        }
        // Square only after division by sqrt(rho_p): this avoids losing a
        // small, physically relevant beta² before a small pressure response
        // rescales it. Keep the nonnegative correction separately even when
        // subtraction from Cp cannot resolve it in the returned Cv.
        let thermal_amplitude = (t.sqrt() * beta) / rho_p.sqrt();
        let delta_cp = thermal_amplitude * thermal_amplitude;
        if !delta_cp.is_finite() || beta != 0.0 && delta_cp == 0.0 {
            return Err(DensityResponseErrorV1::Arithmetic("Cp minus Cv"));
        }
        if delta_cp >= cp {
            return Err(DensityResponseErrorV1::Unstable(
                "density derivatives and Cp imply nonpositive Cv or adiabatic compressibility",
            ));
        }
        let cv = cp - delta_cp;
        if cv / cp <= 16.0 * f64::EPSILON {
            return Err(DensityResponseErrorV1::Arithmetic(
                "adiabatic response is unresolved at the caloric stability boundary",
            ));
        }
        let ratio = cp / cv;
        let kappa_s = kappa_t / ratio;
        let bulk = kappa_s.recip();
        // Root each positive factor independently. bulk/rho can overflow even
        // when its square root (the physical sound speed) is representable.
        let sound = bulk.sqrt() / rho.sqrt();
        if ![cv, ratio, kappa_s, bulk, sound].into_iter().all(positive) {
            return Err(DensityResponseErrorV1::Arithmetic(
                "adiabatic density, caloric or sound response",
            ));
        }
        Ok(Self {
            temperature_k: t,
            density_kg_m3: rho,
            specific_heat_cp_j_kg_k: cp,
            density_temperature_derivative: rho_t,
            density_pressure_derivative: rho_p,
            volumetric_expansion_per_k: beta,
            isothermal_compressibility_per_pa: kappa_t,
            adiabatic_compressibility_per_pa: kappa_s,
            specific_heat_cv_j_kg_k: cv,
            heat_capacity_difference_j_kg_k: delta_cp,
            heat_capacity_ratio: ratio,
            adiabatic_bulk_modulus_pa: bulk,
            sound_speed_m_s: sound,
        })
    }

    /// Original absolute temperature [K].
    pub const fn temperature_k(self) -> f64 {
        self.temperature_k
    }
    /// Original mass density [kg/m³].
    pub const fn density_kg_m3(self) -> f64 {
        self.density_kg_m3
    }
    /// Original isobaric specific heat [J/(kg K)].
    pub const fn specific_heat_cp_j_kg_k(self) -> f64 {
        self.specific_heat_cp_j_kg_k
    }
    /// Original signed isobaric density slope [kg/(m³ K)].
    pub const fn density_temperature_derivative(self) -> DensityTemperatureDerivativeV1 {
        DensityTemperatureDerivativeV1::new(self.density_temperature_derivative)
    }
    /// Original positive isothermal density slope [kg/(m³ Pa)].
    pub const fn density_pressure_derivative(self) -> DensityPressureDerivativeV1 {
        DensityPressureDerivativeV1::new(self.density_pressure_derivative)
    }
    /// Signed isobaric volumetric expansion [1/K], including a density maximum.
    pub const fn volumetric_expansion_per_k(self) -> f64 {
        self.volumetric_expansion_per_k
    }
    /// Isothermal compressibility [1/Pa].
    pub const fn isothermal_compressibility_per_pa(self) -> f64 {
        self.isothermal_compressibility_per_pa
    }
    /// Isentropic compressibility [1/Pa], including the thermal correction.
    pub const fn adiabatic_compressibility_per_pa(self) -> f64 {
        self.adiabatic_compressibility_per_pa
    }
    /// Derived isochoric specific heat [J/(kg K)].
    pub const fn specific_heat_cv_j_kg_k(self) -> f64 {
        self.specific_heat_cv_j_kg_k
    }
    /// Directly evaluated Cp-Cv [J/(kg K)], without subtracting large heat capacities.
    pub const fn heat_capacity_difference_j_kg_k(self) -> f64 {
        self.heat_capacity_difference_j_kg_k
    }
    /// Cp/Cv, at least one, from the same density response.
    pub const fn heat_capacity_ratio(self) -> f64 {
        self.heat_capacity_ratio
    }
    /// Adiabatic bulk modulus [Pa].
    pub const fn adiabatic_bulk_modulus_pa(self) -> f64 {
        self.adiabatic_bulk_modulus_pa
    }
    /// Adiabatic small-signal sound speed [m/s].
    pub const fn sound_speed_m_s(self) -> f64 {
        self.sound_speed_m_s
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

#[cfg(test)]
mod tests;
