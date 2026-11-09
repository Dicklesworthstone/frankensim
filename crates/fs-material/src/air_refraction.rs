//! Ciddor atmospheric-air refraction with explicit pressure, temperature and composition.
//!
//! Implements the NIST Engineering Metrology Toolbox's Ciddor calculation,
//! Appendix A-III, equations A21--A41:
//! <https://emtoolbox.nist.gov/Wavelength/Documentation.asp>.
//! This names the historical formulation, including its 1998 gas constant;
//! it does not claim the refinements of subsequent refractometry experiments.
//!
//! Vacuum wavelength is 300--1700 nm, temperature 233.15--373.15 K, absolute
//! pressure 60--120 kPa, dry-air CO2 0--2000 micromol/mol, and water-vapor mole
//! fraction 0--0.2. These numerical bounds do not prove a homogeneous gas phase:
//! the state owner must exclude condensation, droplets and other constituents.
//! Water fraction is a mole fraction, NOT relative humidity. Composition stays
//! fixed during pressure/temperature differentiation. The empirical density
//! closure is part of this named optical model; it never replaces a fluid EOS.
//! No absorption, scattering, moving-medium correction or uncertainty reduction
//! is inferred. Preparation and queries use fixed storage without allocation.

use fs_blake3::{ContentHash, DomainHasher};

pub const SOURCE_URL: &str = "https://emtoolbox.nist.gov/Wavelength/Documentation.asp";
pub const WAVELENGTH_RANGE_NM: [f64; 2] = [300.0, 1700.0];
pub const TEMPERATURE_RANGE_K: [f64; 2] = [233.15, 373.15];
pub const PRESSURE_RANGE_PA: [f64; 2] = [60_000.0, 120_000.0];

// NIST A21--A30. Retain the source's conventions rather than mixing releases.
const DRY: [f64; 4] = [238.0185, 5_792_105.0, 57.362, 167_917.0];
const VAPOR: [f64; 4] = [295.235, 2.6422, -0.03238, 0.004028];
const A: [f64; 3] = [1.58123e-6, -2.9331e-8, 1.1043e-10];
const B: [f64; 2] = [5.707e-6, -2.051e-8];
const C: [f64; 2] = [1.9898e-4, -2.376e-6];
const D: f64 = 1.83e-11;
const E: f64 = -0.765e-8;
const DRY_SCALE: f64 = 0.9995922115 * 288.15 / 101_325.0;
const VAPOR_SCALE: f64 = 0.018015 / (0.00985938 * 8.314472);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AirRefractionError {
    CompositionOutsideRange,
    WavelengthOutsideRange,
    TemperatureOutsideRange,
    PressureOutsideRange,
    InvalidRelativeIndex,
    NumericalRange,
}
impl core::fmt::Display for AirRefractionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "NIST Ciddor atmospheric-air refraction: {self:?}")
    }
}
impl core::error::Error for AirRefractionError {}

/// One explicitly declared homogeneous atmospheric composition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CiddorAir {
    water_mole_fraction: f64,
    dry_air_co2_micromol_per_mol: f64,
}
impl CiddorAir {
    pub fn try_new(
        water_mole_fraction: f64,
        dry_air_co2_micromol_per_mol: f64,
    ) -> Result<Self, AirRefractionError> {
        if !in_range(water_mole_fraction, [0.0, 0.2])
            || !in_range(dry_air_co2_micromol_per_mol, [0.0, 2000.0])
        {
            return Err(AirRefractionError::CompositionOutsideRange);
        }
        Ok(Self {
            water_mole_fraction: water_mole_fraction + 0.0,
            dry_air_co2_micromol_per_mol: dry_air_co2_micromol_per_mol + 0.0,
        })
    }
    pub const fn water_mole_fraction(self) -> f64 {
        self.water_mole_fraction
    }
    pub const fn dry_air_co2_micromol_per_mol(self) -> f64 {
        self.dry_air_co2_micromol_per_mol
    }
    /// Original equations, coefficient bits, admitted domains and composition.
    pub fn identity(self) -> ContentHash {
        let mut hash = DomainHasher::new("fs-material:nist-ciddor-air:v1");
        hash.update(SOURCE_URL.as_bytes());
        for value in DRY
            .into_iter()
            .chain(VAPOR)
            .chain(A)
            .chain(B)
            .chain(C)
            .chain([
                D,
                E,
                DRY_SCALE,
                VAPOR_SCALE,
                1.022e-8,
                5.34e-7,
                450.0,
                273.15,
            ])
            .chain(WAVELENGTH_RANGE_NM)
            .chain(TEMPERATURE_RANGE_K)
            .chain(PRESSURE_RANGE_PA)
            .chain([
                0.2,
                2000.0,
                self.water_mole_fraction,
                self.dry_air_co2_micromol_per_mol,
            ])
        {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.finalize()
    }
    /// Prepare exact phase, group and curvature coefficients once per vacuum wavelength.
    pub fn prepare(self, wavelength_nm: f64) -> Result<PreparedCiddorAir, AirRefractionError> {
        if !in_range(wavelength_nm, WAVELENGTH_RANGE_NM) {
            return Err(AirRefractionError::WavelengthOutsideRange);
        }
        let s = (1000.0 / wavelength_nm).powi(2);
        let mut dry = [0.0; 3];
        for (pole, strength) in [(DRY[0], DRY[1]), (DRY[2], DRY[3])] {
            let q = pole - s;
            dry[0] += strength / q;
            dry[1] += -2.0 * s * strength / (q * q);
            dry[2] += strength / (q * q) * (6.0 * s + 8.0 * s * s / q);
        }
        let vapor = [
            VAPOR[0] + s * (VAPOR[1] + s * (VAPOR[2] + s * VAPOR[3])),
            -2.0 * s * (VAPOR[1] + s * (2.0 * VAPOR[2] + 3.0 * s * VAPOR[3])),
            6.0 * s * VAPOR[1] + 20.0 * s * s * VAPOR[2] + 42.0 * s * s * s * VAPOR[3],
        ];
        let x = self.water_mole_fraction;
        let dry_weight = DRY_SCALE
            * (1.0 - x)
            * (1.0 + 5.34e-7 * (self.dry_air_co2_micromol_per_mol - 450.0))
            * 1e-8;
        let vapor_weight = VAPOR_SCALE * x * 1.022e-8;
        let spectral: [f64; 3] =
            core::array::from_fn(|i| dry_weight.mul_add(dry[i], vapor_weight * vapor[i]));
        Ok(PreparedCiddorAir {
            composition: self,
            wavelength_nm,
            spectral,
        })
    }
}

/// Fixed spectral preparation, reusable across original pressure/temperature states.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedCiddorAir {
    composition: CiddorAir,
    wavelength_nm: f64,
    // K, lambda K', lambda² K'' in n-1 = K p/(T Z).
    spectral: [f64; 3],
}

/// Absolute indices and derivatives at fixed original molecular composition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AirRefractiveSample {
    pub phase_index: f64,
    pub phase_index_minus_one: f64,
    pub group_index: f64,
    pub wavelength_derivative_per_nm: f64,
    pub wavelength_second_derivative_per_nm2: f64,
    pub temperature_derivative_per_k: f64,
    pub pressure_derivative_per_pa: f64,
    pub temperature_second_derivative_per_k2: f64,
    pub pressure_second_derivative_per_pa2: f64,
    pub temperature_pressure_derivative_per_k_pa: f64,
}

impl PreparedCiddorAir {
    pub const fn wavelength_nm(self) -> f64 {
        self.wavelength_nm
    }
    pub const fn composition(self) -> CiddorAir {
        self.composition
    }
    /// Analytic source evaluation: pressure is absolute Pa and temperature is K.
    /// No spatial differences, source fitting or altered fluid state is used.
    pub fn sample(
        self,
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<AirRefractiveSample, AirRefractionError> {
        if !in_range(temperature_k, TEMPERATURE_RANGE_K) {
            return Err(AirRefractionError::TemperatureOutsideRange);
        }
        if !in_range(pressure_pa, PRESSURE_RANGE_PA) {
            return Err(AirRefractionError::PressureOutsideRange);
        }
        let t = temperature_k - 273.15;
        let p = pressure_pa;
        let inv_t = temperature_k.recip();
        let x = self.composition.water_mole_fraction;
        let c = A[0] + t * (A[1] + t * A[2]) + x * (B[0] + B[1] * t + x * (C[0] + C[1] * t));
        let c_t = A[1] + 2.0 * A[2] * t + x * (B[1] + x * C[1]);
        let b = D + E * x * x;
        // d = T Z. Differentiating this rational form retains compressibility
        // and all mixed pressure/temperature terms without an ideal-gas shortcut.
        let d = temperature_k - p * c + p * p * b * inv_t;
        let dp = -c + 2.0 * p * b * inv_t;
        let dt = 1.0 - p * c_t - p * p * b * inv_t * inv_t;
        let dpp = 2.0 * b * inv_t;
        let dtt = -2.0 * p * A[2] + 2.0 * p * p * b * inv_t.powi(3);
        let dpt = -c_t - 2.0 * p * b * inv_t * inv_t;
        let ap = dp / d;
        let at = dt / d;
        let density_factor = p / d;
        let r = self.spectral[0] * density_factor;
        let sample = AirRefractiveSample {
            phase_index: 1.0 + r,
            phase_index_minus_one: r,
            group_index: 1.0 + (self.spectral[0] - self.spectral[1]) * density_factor,
            wavelength_derivative_per_nm: self.spectral[1] * density_factor / self.wavelength_nm,
            wavelength_second_derivative_per_nm2: self.spectral[2] * density_factor
                / self.wavelength_nm.powi(2),
            temperature_derivative_per_k: -r * at,
            pressure_derivative_per_pa: r * (p.recip() - ap),
            temperature_second_derivative_per_k2: r * (2.0 * at * at - dtt / d),
            pressure_second_derivative_per_pa2: r * (2.0 * ap * ap - dpp / d - 2.0 * ap / p),
            temperature_pressure_derivative_per_k_pa: r * (2.0 * ap * at - dpt / d - at / p),
        };
        if d <= 0.0
            || ![
                sample.phase_index,
                sample.phase_index_minus_one,
                sample.group_index,
                sample.wavelength_derivative_per_nm,
                sample.wavelength_second_derivative_per_nm2,
                sample.temperature_derivative_per_k,
                sample.pressure_derivative_per_pa,
                sample.temperature_second_derivative_per_k2,
                sample.pressure_second_derivative_per_pa2,
                sample.temperature_pressure_derivative_per_k_pa,
            ]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(AirRefractionError::NumericalRange);
        }
        Ok(sample)
    }
}

impl AirRefractiveSample {
    /// Convert n_glass/n_air and its group quantity to absolute vacuum indices.
    /// Both relative inputs must be functions of the SAME vacuum wavelength as
    /// this sample; a fit parameterized by wavelength in air needs a separate
    /// coordinate conversion. Missing relative dispersion remains unavailable.
    /// The result can enter the original SCHOTT thermo-optic `prepare` directly.
    pub fn absolute_from_relative(
        self,
        relative_phase: f64,
        relative_group: Option<f64>,
    ) -> Result<(f64, Option<f64>), AirRefractionError> {
        if !relative_phase.is_finite()
            || relative_phase <= 0.0
            || relative_group.is_some_and(|v| !v.is_finite() || v <= 0.0)
        {
            return Err(AirRefractionError::InvalidRelativeIndex);
        }
        let phase = self.phase_index * relative_phase;
        let group = relative_group.map(|g| {
            self.phase_index
                .mul_add(g, relative_phase * (self.group_index - self.phase_index))
        });
        if !phase.is_finite() || group.is_some_and(|v| !v.is_finite() || v <= 0.0) {
            return Err(AirRefractionError::NumericalRange);
        }
        Ok((phase, group))
    }
}

fn in_range(value: f64, bounds: [f64; 2]) -> bool {
    value.is_finite() && (bounds[0]..=bounds[1]).contains(&value)
}

#[cfg(test)]
mod tests;
