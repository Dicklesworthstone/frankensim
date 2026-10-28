//! Lower U.S. Standard Atmosphere 1976, with explicit geometric height.
//!
//! NOAA/NASA/USAF, NASA-TM-X-74335 (1976), table 4 and equations 17--19,
//! 23, 33a/b: <https://ntrs.nasa.gov/citations/19770009539>.
//! The seven molecular-scale temperature layers determine hydrostatic pressure;
//! density, sound speed and transport use the existing USSA dry-air owner.
//! Layer base pressures are derived once, not rounded independently.
//!
//! This well-mixed lower-atmosphere model admits geometric heights 0..80 km.
//! Above 80 km the source introduces a varying molecular mass (table 8), so
//! extrapolation is refused. This is a standard reference, not local weather,
//! wind, humidity, terrain, oblateness or an inertial Earth gravity model.

use super::{GasSpec, GasState, R_USSA_1976};
use crate::MaterialError;
use fs_math::det;

/// Effective radius used by the standard's geopotential mapping [m].
pub const USSA_EFFECTIVE_RADIUS_M: f64 = 6_356_766.0;
/// Standard sea-level effective gravity [m/s²], including its reference convention.
pub const USSA_SURFACE_GRAVITY_M_S2: f64 = 9.80665;
/// Top of the constant-molecular-mass branch [geometric m].
pub const USSA_MAX_GEOMETRIC_ALTITUDE_M: f64 = 80_000.0;

const HEIGHTS_M: [f64; 7] = [
    0.0, 11_000.0, 20_000.0, 32_000.0, 47_000.0, 51_000.0, 71_000.0,
];
const LAPSE_K_M: [f64; 7] = [-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002];

/// One hydrostatic atmosphere and the common complete numerical gas state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StandardAtmosphereState {
    /// Geometric height above the standard's sea-level surface [m].
    pub geometric_altitude_m: f64,
    /// Geopotential height used by the temperature-layer equations [m].
    pub geopotential_altitude_m: f64,
    /// Effective standard gravity, consistent with the hydrostatic profile [m/s²].
    /// A flight model must separately declare its inertial/rotating frame.
    pub reference_gravity_m_s2: f64,
    /// Existing coherent gas-property snapshot at the actual profile T/p.
    pub gas: GasState,
}

/// Immutable seven-layer reference atmosphere with precomputed base states.
/// Construction and each query allocate no heap storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StandardAtmosphere1976 {
    temperatures_k: [f64; 7],
    pressures_pa: [f64; 7],
}

impl Default for StandardAtmosphere1976 {
    fn default() -> Self {
        Self::new()
    }
}

fn pressure_ratio(temperature_k: f64, lapse_k_m: f64, height_m: f64) -> f64 {
    let r = R_USSA_1976 / GasSpec::dry_air_ussa1976().molar_mass;
    let exponent = if lapse_k_m == 0.0 {
        -USSA_SURFACE_GRAVITY_M_S2 * height_m / (r * temperature_k)
    } else {
        // ln1p retains small height changes close to an exact layer boundary.
        -USSA_SURFACE_GRAVITY_M_S2 / (r * lapse_k_m)
            * (lapse_k_m * height_m / temperature_k).ln_1p()
    };
    det::exp(exponent)
}

impl StandardAtmosphere1976 {
    /// Derive all layer bases from the standard's 288.15 K and 101325 Pa surface.
    #[must_use]
    pub fn new() -> Self {
        let mut temperatures_k = [288.15; 7];
        let mut pressures_pa = [101_325.0; 7];
        for i in 1..7 {
            let dz = HEIGHTS_M[i] - HEIGHTS_M[i - 1];
            temperatures_k[i] = temperatures_k[i - 1] + LAPSE_K_M[i - 1] * dz;
            pressures_pa[i] = pressures_pa[i - 1]
                * pressure_ratio(temperatures_k[i - 1], LAPSE_K_M[i - 1], dz);
        }
        Self {
            temperatures_k,
            pressures_pa,
        }
    }

    /// Query geometric altitude, never silently interpreting it as geopotential.
    pub fn at_geometric_altitude(
        &self,
        altitude_m: f64,
    ) -> Result<StandardAtmosphereState, MaterialError> {
        if !altitude_m.is_finite() || !(0.0..=USSA_MAX_GEOMETRIC_ALTITUDE_M).contains(&altitude_m) {
            return Err(MaterialError::Parameters {
                what: "USSA well-mixed profile requires geometric altitude in [0,80000] m".into(),
            });
        }
        let ratio = USSA_EFFECTIVE_RADIUS_M / (USSA_EFFECTIVE_RADIUS_M + altitude_m);
        self.sample(altitude_m, altitude_m * ratio)
    }

    /// Query the source's geopotential height with the same physical top bound.
    pub fn at_geopotential_altitude(
        &self,
        altitude_m: f64,
    ) -> Result<StandardAtmosphereState, MaterialError> {
        let top = USSA_MAX_GEOMETRIC_ALTITUDE_M * USSA_EFFECTIVE_RADIUS_M
            / (USSA_EFFECTIVE_RADIUS_M + USSA_MAX_GEOMETRIC_ALTITUDE_M);
        if !altitude_m.is_finite() || !(0.0..=top).contains(&altitude_m) {
            return Err(MaterialError::Parameters {
                what: "geopotential altitude lies outside the USSA 0..80 km geometric branch".into(),
            });
        }
        let geometric = altitude_m / (1.0 - altitude_m / USSA_EFFECTIVE_RADIUS_M);
        self.sample(geometric, altitude_m)
    }

    fn sample(
        &self,
        geometric: f64,
        geopotential: f64,
    ) -> Result<StandardAtmosphereState, MaterialError> {
        let layer = HEIGHTS_M.partition_point(|&h| h <= geopotential) - 1;
        let dh = geopotential - HEIGHTS_M[layer];
        let temperature = self.temperatures_k[layer] + LAPSE_K_M[layer] * dh;
        let pressure = self.pressures_pa[layer]
            * pressure_ratio(self.temperatures_k[layer], LAPSE_K_M[layer], dh);
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure)?;
        let radius_ratio = USSA_EFFECTIVE_RADIUS_M / (USSA_EFFECTIVE_RADIUS_M + geometric);
        Ok(StandardAtmosphereState {
            geometric_altitude_m: geometric,
            geopotential_altitude_m: geopotential,
            reference_gravity_m_s2: USSA_SURFACE_GRAVITY_M_S2 * radius_ratio * radius_ratio,
            gas,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, relative: f64) {
        assert!(
            (a - b).abs() <= relative * b.abs().max(1e-15),
            "{a} != {b}"
        );
    }

    #[test]
    fn published_layer_values_and_original_surface_gas_agree() {
        let model = StandardAtmosphere1976::new();
        let zero = model.at_geometric_altitude(0.0).unwrap();
        assert_eq!(
            zero.gas,
            GasState::try_new(&GasSpec::dry_air_ussa1976(), 288.15, 101325.0).unwrap()
        );
        // Source table 4 temperatures; pressures are rounded standard layer
        // values, whose last displayed digit bounds the comparison.
        for (h, t, p) in [
            (11000.0, 216.65, 22632.1),
            (20000.0, 216.65, 5474.89),
            (32000.0, 228.65, 868.019),
            (47000.0, 270.65, 110.906),
            (51000.0, 270.65, 66.9389),
            (71000.0, 214.65, 3.95642),
        ] {
            let state = model.at_geopotential_altitude(h).unwrap();
            close(state.gas.temperature, t, 1e-14);
            close(state.gas.pressure, p, 5e-6);
            state.gas.validate_coherence().unwrap();
        }
    }

    #[test]
    fn geometric_mapping_and_hydrostatic_derivative_use_the_same_gravity() {
        let model = StandardAtmosphere1976::new();
        for z in [1000.0, 15000.0, 30000.0, 50000.0, 79999.0] {
            let state = model.at_geometric_altitude(z).unwrap();
            let inverse = model
                .at_geopotential_altitude(state.geopotential_altitude_m)
                .unwrap();
            close(inverse.geometric_altitude_m, z, 3e-15);
            close(inverse.gas.pressure, state.gas.pressure, 3e-14);
            let dp_dz = (model.at_geometric_altitude(z + 0.1).unwrap().gas.pressure
                - model.at_geometric_altitude(z - 0.1).unwrap().gas.pressure)
                / 0.2;
            close(-dp_dz, state.gas.density * state.reference_gravity_m_s2, 2e-9);
        }
    }

    #[test]
    fn layers_are_continuous_and_upper_composition_extrapolation_refuses() {
        let model = StandardAtmosphere1976::new();
        for h in HEIGHTS_M.into_iter().skip(1) {
            let left = model.at_geopotential_altitude(h - 1e-6).unwrap();
            let right = model.at_geopotential_altitude(h + 1e-6).unwrap();
            close(left.gas.temperature, right.gas.temperature, 1e-10);
            close(left.gas.pressure, right.gas.pressure, 1e-9);
        }
        assert!(model.at_geometric_altitude(80000.0).is_ok());
        for h in [-1e-12, 80000.00001, f64::NAN, f64::INFINITY] {
            assert!(model.at_geometric_altitude(h).is_err());
        }
        assert!(model.at_geopotential_altitude(80000.0).is_err());
    }
}
