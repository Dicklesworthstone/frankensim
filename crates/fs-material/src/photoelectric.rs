//! Cold, single-photon photoemission from an explicitly characterized surface.
//!
//! Photon energy is `h*c/lambda` for **vacuum** wavelength. Emission requires
//! `h*c/lambda > work_function`; the excess is the maximum, not the mean,
//! electron kinetic energy. The caller supplies both the yield per **absorbed**
//! photon and the mean fraction of the excess carried away by an electron.
//! Neither quantity is inferred from a bulk material name. Reflection and
//! transmission belong to the optical owner before this law is evaluated;
//! external quantum efficiency must first be converted to absorbed-photon yield.
//!
//! Removal work, emitted kinetic power and heat are separate accounts. Removal
//! work is not automatically thermalized. Charge has equal and opposite surface
//! and emitted-electron currents. These are mean emission rates, not stochastic
//! integer counts or collector currents. The stopping potential is the ideal
//! vacuum retarding barrier for the fastest electron, not a measured electrode
//! voltage with contact potentials.
//!
//! This bounded constitutive evaluation allocates nothing. It does not model
//! electron transport, finite-temperature tails, field/Schottky emission,
//! multiphoton excitation, space charge, secondary gain or semiconductor bands.
//! Physical accuracy requires independent calibration of the supplied surface
//! parameters, including wavelength/temperature support, by the caller.
//!
//! References: [NIST SI defining constants](https://www.nist.gov/pml/special-publication-330/sp-330-section-2),
//! [UCLA photoelectric-effect experiment](https://demoweb.physics.ucla.edu/node/118),
//! and Dowell & Schmerge, *Quantum efficiency and thermal emittance of metal
//! photocathodes*, [revised SLAC-PUB-13535](https://www.slac.stanford.edu/pubs/slacpubs/13500/slac-pub-13535.pdf).
//! The latter distinguishes emission thresholds from electron energy statistics;
//! its particular distribution model is not silently assumed here.

/// Exact SI Planck constant, in joule seconds.
pub const PLANCK_J_S: f64 = 6.626_070_15e-34;
/// Exact SI speed of light in vacuum, in metres per second.
pub const LIGHT_M_S: f64 = 299_792_458.0;
/// Exact SI magnitude of the electron charge, in coulombs.
pub const ELEMENTARY_CHARGE_C: f64 = 1.602_176_634e-19;

/// Invalid physical parameters or an unrepresentable response.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhotoelectricError {
    Invalid(&'static str),
}

impl core::fmt::Display for PhotoelectricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "invalid photoelectric response: {reason}"),
        }
    }
}

impl core::error::Error for PhotoelectricError {}

/// Immutable numeric law at one explicitly supplied surface condition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotoelectricLaw {
    work_function_j: f64,
    absorbed_quantum_yield: f64,
    mean_excess_fraction: f64,
}

impl PhotoelectricLaw {
    /// Construct a law from positive removal work and fractions in `[0, 1]`.
    ///
    /// `absorbed_quantum_yield` is emitted electrons per absorbed photon above
    /// threshold. `mean_excess_fraction` is mean kinetic energy divided by the
    /// maximum excess energy. No material-universal default is supplied.
    pub fn try_new(
        work_function_j: f64,
        absorbed_quantum_yield: f64,
        mean_excess_fraction: f64,
    ) -> Result<Self, PhotoelectricError> {
        if !work_function_j.is_finite() || work_function_j <= 0.0 {
            return Err(PhotoelectricError::Invalid(
                "work function must be finite and positive in J",
            ));
        }
        if !absorbed_quantum_yield.is_finite() || !(0.0..=1.0).contains(&absorbed_quantum_yield) {
            return Err(PhotoelectricError::Invalid(
                "absorbed-photon yield must be finite and in [0, 1]",
            ));
        }
        if !mean_excess_fraction.is_finite() || !(0.0..=1.0).contains(&mean_excess_fraction) {
            return Err(PhotoelectricError::Invalid(
                "mean excess fraction must be finite and in [0, 1]",
            ));
        }
        Ok(Self {
            work_function_j,
            absorbed_quantum_yield,
            mean_excess_fraction,
        })
    }

    pub const fn work_function_j(&self) -> f64 {
        self.work_function_j
    }

    pub const fn absorbed_quantum_yield(&self) -> f64 {
        self.absorbed_quantum_yield
    }

    pub const fn mean_excess_fraction(&self) -> f64 {
        self.mean_excess_fraction
    }

    /// Evaluate a monochromatic absorbed-power bin, in vacuum metres and watts.
    ///
    /// Evaluate distinct spectral bins independently; averaging photon energies
    /// before the threshold test changes the physics. Zero power is valid.
    /// Nonfinite arithmetic and a positive output rounded to zero are refused.
    /// Energy partition closes to floating-point roundoff, without subtracting
    /// nearly equal powers to recover a small heat contribution.
    pub fn evaluate(
        &self,
        wavelength_m: f64,
        absorbed_power_w: f64,
    ) -> Result<PhotoelectricResponse, PhotoelectricError> {
        if !wavelength_m.is_finite() || wavelength_m <= 0.0 {
            return Err(PhotoelectricError::Invalid(
                "vacuum wavelength must be finite and positive in m",
            ));
        }
        if !absorbed_power_w.is_finite() || absorbed_power_w < 0.0 {
            return Err(PhotoelectricError::Invalid(
                "absorbed power must be finite and nonnegative in W",
            ));
        }
        let photon_energy_j = (PLANCK_J_S * LIGHT_M_S) / wavelength_m;
        if !photon_energy_j.is_finite() || photon_energy_j <= 0.0 {
            return Err(PhotoelectricError::Invalid(
                "photon energy is outside the positive finite J range",
            ));
        }
        let absorbed_photon_rate_s = absorbed_power_w / photon_energy_j;
        if !absorbed_photon_rate_s.is_finite()
            || (absorbed_power_w > 0.0 && absorbed_photon_rate_s == 0.0)
        {
            return Err(PhotoelectricError::Invalid(
                "absorbed photon rate is outside the finite per-second range",
            ));
        }
        let mut response = PhotoelectricResponse {
            wavelength_m,
            photon_energy_j,
            maximum_kinetic_energy_j: 0.0,
            mean_kinetic_energy_j: 0.0,
            stopping_potential_v: 0.0,
            absorbed_power_w,
            absorbed_photon_rate_s,
            emitted_electron_rate_s: 0.0,
            surface_charging_current_a: 0.0,
            emitted_electron_current_a: 0.0,
            kinetic_power_w: 0.0,
            removal_power_w: 0.0,
            heat_power_w: absorbed_power_w,
        };
        if photon_energy_j <= self.work_function_j {
            return Ok(response);
        }

        let excess_j = photon_energy_j - self.work_function_j;
        let stopping_potential_v = excess_j / ELEMENTARY_CHARGE_C;
        if !stopping_potential_v.is_finite() {
            return Err(PhotoelectricError::Invalid(
                "stopping potential is outside the finite V range",
            ));
        }
        let mean_kinetic_energy_j = product(
            excess_j,
            self.mean_excess_fraction,
            "mean electron kinetic energy is outside the finite J range",
        )?;
        let emitted_electron_rate_s = product(
            absorbed_photon_rate_s,
            self.absorbed_quantum_yield,
            "emitted electron rate is outside the finite per-second range",
        )?;
        let surface_charging_current_a = product(
            emitted_electron_rate_s,
            ELEMENTARY_CHARGE_C,
            "emitted charge rate is outside the finite A range",
        )?;
        let removal_power_w = product(
            emitted_electron_rate_s,
            self.work_function_j,
            "removal power is outside the finite W range",
        )?;
        let excess_power_w = product(
            emitted_electron_rate_s,
            excess_j,
            "emitted excess power is outside the finite W range",
        )?;
        let kinetic_power_w = product(
            excess_power_w,
            self.mean_excess_fraction,
            "kinetic power is outside the finite W range",
        )?;
        // Both heat contributions are nonnegative: absorption without escape,
        // and the excess energy not carried by escaping electrons. Do not form
        // P_abs - P_remove - P_kin, which loses small heat near unit efficiency.
        let nonemitting_heat_w = product(
            absorbed_power_w,
            1.0 - self.absorbed_quantum_yield,
            "nonemitting absorption heat is outside the finite W range",
        )?;
        let relaxation_heat_w = product(
            excess_power_w,
            1.0 - self.mean_excess_fraction,
            "electron relaxation heat is outside the finite W range",
        )?;
        let heat_power_w = nonemitting_heat_w + relaxation_heat_w;
        if !heat_power_w.is_finite() {
            return Err(PhotoelectricError::Invalid(
                "heat power is outside the finite W range",
            ));
        }

        response.maximum_kinetic_energy_j = excess_j;
        response.mean_kinetic_energy_j = mean_kinetic_energy_j;
        response.stopping_potential_v = stopping_potential_v;
        response.emitted_electron_rate_s = emitted_electron_rate_s;
        response.surface_charging_current_a = surface_charging_current_a;
        response.emitted_electron_current_a = -surface_charging_current_a;
        response.kinetic_power_w = kinetic_power_w;
        response.removal_power_w = removal_power_w;
        response.heat_power_w = heat_power_w;
        Ok(response)
    }
}

/// A mean emission response for one absorbed monochromatic power bin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotoelectricResponse {
    /// Vacuum wavelength, in metres.
    pub wavelength_m: f64,
    pub photon_energy_j: f64,
    /// Kinematic upper bound, independent of intensity; zero below threshold.
    pub maximum_kinetic_energy_j: f64,
    /// Authored mean fraction of that bound, not a material-independent result.
    pub mean_kinetic_energy_j: f64,
    /// Ideal vacuum retarding voltage for the fastest electron, in volts.
    pub stopping_potential_v: f64,
    pub absorbed_power_w: f64,
    pub absorbed_photon_rate_s: f64,
    pub emitted_electron_rate_s: f64,
    /// Positive charging of the surface as negative electrons escape, in A.
    pub surface_charging_current_a: f64,
    /// Equal and opposite signed charge flux carried by emitted electrons, in A.
    pub emitted_electron_current_a: f64,
    pub kinetic_power_w: f64,
    /// Electron removal work, in W; must not also be booked as surface heat.
    pub removal_power_w: f64,
    pub heat_power_w: f64,
}

fn product(a: f64, b: f64, reason: &'static str) -> Result<f64, PhotoelectricError> {
    let value = a * b;
    if !value.is_finite() || (a > 0.0 && b > 0.0 && value == 0.0) {
        Err(PhotoelectricError::Invalid(reason))
    } else {
        Ok(value)
    }
}

#[cfg(test)]
#[path = "photoelectric/tests.rs"]
mod tests;
