//! Solid-vapor constitutive laws for fs_thermal::sublimation.
//!
//! These scalar laws allocate no memory and return pressure, its temperature
//! derivative, and a consistent latent heat together. Density is always supplied
//! by the caller; the pressure fit does not determine bulk density or porosity.
//! Named fits refuse extrapolation. The generic law has an explicit caller-owned
//! validity interval and does not infer a material or a phase from its name.
//!
//! For the named fits, L = (R/M) T² d(ln p_sat)/dT uses the ideal-vapor,
//! negligible-solid-volume approximation to Clausius-Clapeyron. This is not an
//! exact real-fluid enthalpy difference. No accommodation, surface curvature,
//! transport resistance, nucleation, or melting law is inferred here. Those
//! interface/environment choices belong to the receiving physical owner.
//!
//! Determinism: fixed-order scalar arithmetic on the same ISA and toolchain.
//! Transcendentals use the shared fs_math deterministic kernels. G0/G1 checks
//! cover published pressure values, derivatives and domains; no cross-platform
//! execution claim is made by these tests.

use fs_thermal::sublimation::{SublimationError, SublimationLaw, SublimationProperties};

const GAS_CONSTANT_J_MOL_K: f64 = 8.314_462_618_153_24;
const WATER_ICE_RANGE_K: [f64; 2] = [50.0, 273.16];

/// IAPWS ice-Ih sublimation pressure on 50..=273.16 K.
///
/// Implements equation (4), Table 3 of Wagner et al. (2011),
/// <https://www.nist.gov/system/files/documents/srd/jpcrd402011New.pdf>.
/// The reference point is 273.16 K, 611.657 Pa. Molar mass is 0.01801528 kg/mol,
/// matching the water-vapor species convention used elsewhere in this crate.
/// The supplied density stays constant over this law's interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WaterIceSublimation {
    solid_density_kg_m3: f64,
}

impl WaterIceSublimation {
    /// Construct with an explicitly selected finite positive solid density.
    pub fn try_new(solid_density_kg_m3: f64) -> Result<Self, SublimationError> {
        positive(
            solid_density_kg_m3,
            "solid density must be finite and positive",
        )?;
        Ok(Self {
            solid_density_kg_m3,
        })
    }

    /// Pure-phase saturation pressure in Pa, independent of bulk solid density.
    ///
    /// Humidity and vapor consumers can share this curve without constructing
    /// an arbitrary density or a phase-transfer owner.
    pub fn saturation_pressure_pa(temperature_k: f64) -> Result<f64, SublimationError> {
        Ok(water_ice_curve(temperature_k)?.0)
    }
}

impl SublimationLaw for WaterIceSublimation {
    fn temperature_range_k(&self) -> [f64; 2] {
        WATER_ICE_RANGE_K
    }

    fn properties(&self, temperature_k: f64) -> Result<SublimationProperties, SublimationError> {
        let (pressure, first_moment, latent_derivative_moment) = water_ice_curve(temperature_k)?;
        let specific_gas_constant = GAS_CONSTANT_J_MOL_K / 0.018_015_28;
        checked_properties(SublimationProperties {
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: pressure * first_moment / temperature_k,
            latent_heat_j_kg: specific_gas_constant * temperature_k * first_moment,
            latent_heat_derivative_j_kg_k: specific_gas_constant * latent_derivative_moment,
            solid_density_kg_m3: self.solid_density_kg_m3,
            molar_mass_kg_mol: 0.018_015_28,
        })
    }
}

// Returns p_sat, T d(ln p_sat)/dT, and d[T² d(ln p_sat)/dT]/dT.
// Reuse the three evaluated powers for pressure and both thermal derivatives.
fn water_ice_curve(temperature_k: f64) -> Result<(f64, f64, f64), SublimationError> {
    admit_temperature(temperature_k, WATER_ICE_RANGE_K)?;
    let theta = temperature_k / 273.16;
    let coefficients = [-21.214_400_6, 27.320_381_9, -6.105_981_30];
    let exponents = [0.003_333_333_33, 1.206_666_67, 1.703_333_33];
    let mut log_pressure_ratio = 0.0;
    let mut first_moment = 0.0;
    let mut latent_derivative_moment = 0.0;
    for (a, b) in coefficients.into_iter().zip(exponents) {
        let term = a * fs_math::det::pow(theta, b - 1.0);
        log_pressure_ratio += term;
        first_moment += term * (b - 1.0);
        latent_derivative_moment += term * (b - 1.0) * b;
    }
    Ok((
        611.657 * fs_math::det::exp(log_pressure_ratio),
        first_moment,
        latent_derivative_moment,
    ))
}

/// Solid CO₂ pressure from the NIST Antoine fit on 154.26..=195.89 K.
///
/// NIST Chemistry WebBook, Giauque and Egan (1937),
/// <https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=384>:
/// log10(p / bar) = 6.81228 - 1301.679 / (T/K - 3.494).
/// Conversion to Pa is explicit; no extension to the CO₂ triple point is made.
/// The latent heat follows the ideal-vapor approximation documented above,
/// rather than being represented as measured calorimetry. Caller density is
/// held constant; molar mass is the WebBook's 0.0440095 kg/mol.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CarbonDioxideSublimation {
    solid_density_kg_m3: f64,
}

impl CarbonDioxideSublimation {
    /// Construct with an explicitly selected finite positive solid density.
    pub fn try_new(solid_density_kg_m3: f64) -> Result<Self, SublimationError> {
        positive(
            solid_density_kg_m3,
            "solid density must be finite and positive",
        )?;
        Ok(Self {
            solid_density_kg_m3,
        })
    }
}

impl SublimationLaw for CarbonDioxideSublimation {
    fn temperature_range_k(&self) -> [f64; 2] {
        [154.26, 195.89]
    }

    fn properties(&self, temperature_k: f64) -> Result<SublimationProperties, SublimationError> {
        admit_temperature(temperature_k, self.temperature_range_k())?;
        let denominator = temperature_k - 3.494;
        let pressure = 100_000.0 * fs_math::det::pow(10.0, 6.812_28 - 1301.679 / denominator);
        let slope_scale = std::f64::consts::LN_10 * 1301.679;
        let ratio = temperature_k / denominator;
        let latent_heat = (GAS_CONSTANT_J_MOL_K / 0.044_009_5) * slope_scale * ratio * ratio;
        checked_properties(SublimationProperties {
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: pressure * slope_scale / denominator / denominator,
            latent_heat_j_kg: latent_heat,
            latent_heat_derivative_j_kg_k: -2.0 * 3.494 * latent_heat / temperature_k / denominator,
            solid_density_kg_m3: self.solid_density_kg_m3,
            molar_mass_kg_mol: 0.044_009_5,
        })
    }
}

/// Constant-latent-heat Clausius-Clapeyron law on a declared temperature range.
///
/// ln(p/p_ref) = (L M/R) (1/T_ref - 1/T). This assumes an ideal vapor,
/// negligible solid specific volume, constant latent heat, and constant solid
/// density. The caller owns calibration and phase validity; a bounded interval
/// is not an accuracy certificate. All inputs are in SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClausiusClapeyronSublimation {
    reference_temperature_k: f64,
    reference_pressure_pa: f64,
    log_reference_pressure: f64,
    latent_heat_j_kg: f64,
    characteristic_temperature_k: f64,
    molar_mass_kg_mol: f64,
    solid_density_kg_m3: f64,
    temperature_range_k: [f64; 2],
}

impl ClausiusClapeyronSublimation {
    /// Admit an explicit reference state, latent heat, species mass and density.
    ///
    /// The finite positive reference temperature must lie in the positive,
    /// strictly ordered validity interval. Nonrepresentable endpoint pressures
    /// and derivatives refuse construction rather than silently producing zero
    /// pressure or infinity. Every subsequent evaluation is checked as well.
    pub fn try_new(
        reference_temperature_k: f64,
        reference_pressure_pa: f64,
        latent_heat_j_kg: f64,
        molar_mass_kg_mol: f64,
        solid_density_kg_m3: f64,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, SublimationError> {
        positive(
            reference_pressure_pa,
            "reference pressure must be finite and positive",
        )?;
        positive(latent_heat_j_kg, "latent heat must be finite and positive")?;
        positive(molar_mass_kg_mol, "molar mass must be finite and positive")?;
        positive(
            solid_density_kg_m3,
            "solid density must be finite and positive",
        )?;
        let [minimum_k, maximum_k] = temperature_range_k;
        positive(minimum_k, "minimum temperature must be finite and positive")?;
        positive(maximum_k, "maximum temperature must be finite and positive")?;
        if minimum_k >= maximum_k {
            return Err(SublimationError::Invalid(
                "temperature range must be strictly increasing",
            ));
        }
        admit_temperature(reference_temperature_k, temperature_range_k)?;
        let characteristic_temperature_k =
            latent_heat_j_kg / (GAS_CONSTANT_J_MOL_K / molar_mass_kg_mol);
        positive(
            characteristic_temperature_k,
            "Clausius-Clapeyron scale is not representable",
        )?;
        let law = Self {
            reference_temperature_k,
            reference_pressure_pa,
            log_reference_pressure: fs_math::det::ln(reference_pressure_pa),
            latent_heat_j_kg,
            characteristic_temperature_k,
            molar_mass_kg_mol,
            solid_density_kg_m3,
            temperature_range_k,
        };
        law.properties(minimum_k)?;
        law.properties(maximum_k)?;
        Ok(law)
    }
}

impl SublimationLaw for ClausiusClapeyronSublimation {
    fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }

    fn properties(&self, temperature_k: f64) -> Result<SublimationProperties, SublimationError> {
        admit_temperature(temperature_k, self.temperature_range_k)?;
        // The difference form retains small changes near the reference state.
        // Log pressure avoids an intermediate exponential overflow when p_ref
        // is small, and returning p_ref directly preserves the exact anchor.
        let pressure = if temperature_k == self.reference_temperature_k {
            self.reference_pressure_pa
        } else {
            let exponent = self.characteristic_temperature_k
                * ((temperature_k - self.reference_temperature_k) / temperature_k)
                / self.reference_temperature_k;
            fs_math::det::exp(self.log_reference_pressure + exponent)
        };
        checked_properties(SublimationProperties {
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: pressure
                * (self.characteristic_temperature_k / temperature_k)
                / temperature_k,
            latent_heat_j_kg: self.latent_heat_j_kg,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: self.solid_density_kg_m3,
            molar_mass_kg_mol: self.molar_mass_kg_mol,
        })
    }
}

fn positive(value: f64, reason: &'static str) -> Result<(), SublimationError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(SublimationError::Invalid(reason))
    }
}

fn admit_temperature(temperature_k: f64, range: [f64; 2]) -> Result<(), SublimationError> {
    positive(temperature_k, "temperature must be finite and positive")?;
    if temperature_k < range[0] || temperature_k > range[1] {
        return Err(SublimationError::TemperatureOutOfRange {
            temperature_k,
            minimum_k: range[0],
            maximum_k: range[1],
        });
    }
    Ok(())
}

fn checked_properties(
    properties: SublimationProperties,
) -> Result<SublimationProperties, SublimationError> {
    positive(
        properties.saturation_pressure_pa,
        "saturation pressure is not representable",
    )?;
    positive(
        properties.saturation_pressure_derivative_pa_k,
        "pressure derivative is not representable",
    )?;
    positive(
        properties.latent_heat_j_kg,
        "latent heat is not representable",
    )?;
    if !properties.latent_heat_derivative_j_kg_k.is_finite() {
        return Err(SublimationError::Invalid(
            "latent heat derivative is not representable",
        ));
    }
    Ok(properties)
}

#[cfg(test)]
mod tests;
