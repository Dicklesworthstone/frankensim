//! First-principles ideal-gas state: the universal ambient-medium
//! primitive (musical-acoustics program doctrine, bead
//! frankensim-fsim-duct-acoustics-zdmm1 slice 1).
//!
//! Every acoustic and viscothermal facility needs the same tuple of
//! medium properties — density, sound speed, dynamic viscosity,
//! thermal conductivity, heat-capacity ratio, Prandtl number,
//! characteristic impedance — and every one of them is a FUNCTION of
//! ambient temperature, pressure, and gas identity, not a constant.
//! This module derives the complete tuple bottoms-up:
//!
//! - `rho = p M / (R T)` (ideal gas),
//! - `c = sqrt(gamma R T / M)`,
//! - `mu = beta T^{3/2} / (T + S)` (Sutherland),
//! - `kappa` by the gas's declared model: the US Standard Atmosphere
//!   1976 empirical air fit, or the kinetic-theory Eucken relation
//!   `kappa = mu (cp + 5/4 R/M)` for a generic gas,
//! - `cp = gamma R / ((gamma - 1) M)` (calorically perfect),
//! - `Pr = mu cp / kappa` — DERIVED, never an input: for USSA-1976 dry
//!   air the classic 0.71 emerges from the derivation, which is the
//!   module's built-in kinetic-theory falsifier.
//!
//! Provenance: the dry-air constants and both empirical fits are from
//! the *U.S. Standard Atmosphere, 1976* (NOAA-S/T 76-1562, NASA/NOAA/
//! USAF — U.S. Government work, public domain): `R* = 8.31432`
//! J/(mol K) (the document's own constant, retained so its printed
//! sea-level values reproduce exactly; CODATA 2018 differs in the 5th
//! decimal), `M_air = 28.9644e-3` kg/mol, Sutherland
//! `beta = 1.458e-6`, `S = 110.4 K`, and the conductivity fit
//! `kappa = 2.64638e-3 T^{3/2} / (T + 245.4 * 10^{-12/T})`.
//! Pinned printed values at sea level (288.15 K, 101325 Pa):
//! `rho = 1.2250`, `c = 340.29`, `mu = 1.7894e-5`.
//!
//! Validity (refused outside, documented inside): `T` in [50, 2000] K
//! and `p` in (0, 1e7] Pa admit evaluation; the ideal-gas law is good
//! to well under 1% below ~10 bar away from condensation, and the
//! calorically-perfect `gamma` for diatomic air degrades above ~600 K
//! as vibrational modes excite (order 1-2% at 600 K, growing with T) —
//! stated as an accuracy boundary, not silently hidden, because the
//! extreme-regime doctrine (hot enclosures, melting strings) wants
//! evaluation there with honest error language rather than refusal.
//!
//! Determinism: pure `f64` arithmetic with `fs_math::det` for the
//! non-IEEE transcendentals; repeat evaluations are bitwise identical.

use fs_math::det;

use crate::MaterialError;

/// USSA-1976 universal gas constant [J/(mol K)] — the document's own
/// value, kept so its printed tables reproduce bit-honestly.
pub const R_USSA_1976: f64 = 8.31432;

/// How a [`GasSpec`] derives thermal conductivity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConductivityModel {
    /// The USSA-1976 empirical dry-air fit
    /// `kappa = 2.64638e-3 T^{3/2} / (T + 245.4 * 10^{-12/T})`.
    Ussa1976AirFit,
    /// Kinetic-theory Eucken relation `kappa = mu (cp + 5/4 R/M)` —
    /// the generic route for any gas with a Sutherland viscosity.
    Eucken,
}

/// A calorically-perfect gas: molar mass, heat-capacity ratio,
/// Sutherland viscosity constants, and the conductivity model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasSpec {
    /// Molar mass M [kg/mol].
    pub molar_mass: f64,
    /// Heat-capacity ratio gamma (calorically perfect).
    pub gamma: f64,
    /// Sutherland coefficient beta [kg/(s m K^0.5)].
    pub sutherland_beta: f64,
    /// Sutherland temperature S [K].
    pub sutherland_s: f64,
    /// Thermal-conductivity derivation.
    pub conductivity: ConductivityModel,
}

impl GasSpec {
    /// Dry air per the U.S. Standard Atmosphere, 1976 (public domain):
    /// `M = 28.9644e-3`, `gamma = 1.4`, Sutherland `beta = 1.458e-6`,
    /// `S = 110.4`, USSA conductivity fit.
    #[must_use]
    pub const fn dry_air_ussa1976() -> GasSpec {
        GasSpec {
            molar_mass: 28.9644e-3,
            gamma: 1.4,
            sutherland_beta: 1.458e-6,
            sutherland_s: 110.4,
            conductivity: ConductivityModel::Ussa1976AirFit,
        }
    }
}

/// The complete derived ambient state of a gas at (T, p): everything a
/// linear-acoustic, viscothermal, or convective facility consumes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasState {
    /// Temperature [K].
    pub temperature: f64,
    /// Pressure [Pa].
    pub pressure: f64,
    /// Density rho [kg/m^3].
    pub density: f64,
    /// Sound speed c [m/s].
    pub sound_speed: f64,
    /// Dynamic viscosity mu [Pa s].
    pub dynamic_viscosity: f64,
    /// Thermal conductivity kappa [W/(m K)].
    pub thermal_conductivity: f64,
    /// Heat-capacity ratio gamma.
    pub gamma: f64,
    /// Specific gas constant R/M [J/(kg K)].
    pub specific_gas_constant: f64,
    /// Isobaric specific heat cp [J/(kg K)].
    pub specific_heat_cp: f64,
    /// Prandtl number mu cp / kappa (derived).
    pub prandtl: f64,
    /// Characteristic acoustic impedance rho c [Pa s/m].
    pub characteristic_impedance: f64,
}

impl GasState {
    /// Derive the full state from first principles.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] for a non-physical spec or a
    /// (T, p) outside the documented validity window (T in [50, 2000]
    /// K, p in (0, 1e7] Pa).
    pub fn try_new(spec: &GasSpec, temperature: f64, pressure: f64) -> Result<Self, MaterialError> {
        for (value, what) in [
            (spec.molar_mass, "molar mass must be positive and finite"),
            (
                spec.sutherland_beta,
                "Sutherland beta must be positive and finite",
            ),
            (
                spec.sutherland_s,
                "Sutherland S must be positive and finite",
            ),
        ] {
            if !(value > 0.0 && value.is_finite()) {
                return Err(MaterialError::Parameters { what: what.into() });
            }
        }
        if !(spec.gamma > 1.0 && spec.gamma < 5.0 / 3.0 + 1e-9) {
            return Err(MaterialError::Parameters {
                what: format!(
                    "gamma = {} outside (1, 5/3] (monatomic ceiling)",
                    spec.gamma
                ),
            });
        }
        if !(50.0..=2000.0).contains(&temperature) {
            return Err(MaterialError::Parameters {
                what: format!(
                    "temperature {temperature} K outside the [50, 2000] K validity window"
                ),
            });
        }
        if !(pressure > 0.0 && pressure <= 1.0e7) {
            return Err(MaterialError::Parameters {
                what: format!("pressure {pressure} Pa outside the (0, 1e7] Pa validity window"),
            });
        }
        let r_specific = R_USSA_1976 / spec.molar_mass;
        let density = pressure / (r_specific * temperature);
        let sound_speed = (spec.gamma * r_specific * temperature).sqrt();
        let t32 = temperature * temperature.sqrt();
        let dynamic_viscosity = spec.sutherland_beta * t32 / (temperature + spec.sutherland_s);
        let specific_heat_cp = spec.gamma * r_specific / (spec.gamma - 1.0);
        let thermal_conductivity = match spec.conductivity {
            ConductivityModel::Ussa1976AirFit => {
                // kappa = 2.64638e-3 T^{3/2} / (T + 245.4 * 10^{-12/T}).
                let exponent = -12.0 / temperature;
                let ten_pow = det::exp(exponent * core::f64::consts::LN_10);
                2.64638e-3 * t32 / (temperature + 245.4 * ten_pow)
            }
            ConductivityModel::Eucken => dynamic_viscosity * (specific_heat_cp + 1.25 * r_specific),
        };
        let prandtl = dynamic_viscosity * specific_heat_cp / thermal_conductivity;
        Ok(GasState {
            temperature,
            pressure,
            density,
            sound_speed,
            dynamic_viscosity,
            thermal_conductivity,
            gamma: spec.gamma,
            specific_gas_constant: r_specific,
            specific_heat_cp,
            prandtl,
            characteristic_impedance: density * sound_speed,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ussa_1976_sea_level_pins() {
        // The document's printed sea-level values (288.15 K, 101325
        // Pa): rho = 1.2250 kg/m^3, c = 340.29 m/s, mu = 1.7894e-5
        // Pa s; kappa ~ 2.5326e-2 W/(m K). Pinned to the printed
        // precision.
        let air = GasSpec::dry_air_ussa1976();
        let state = GasState::try_new(&air, 288.15, 101_325.0).expect("sea level");
        assert!(
            (state.density - 1.2250).abs() < 5e-5,
            "rho {}",
            state.density
        );
        assert!(
            (state.sound_speed - 340.294).abs() < 5e-3,
            "c {}",
            state.sound_speed
        );
        assert!(
            (state.dynamic_viscosity - 1.7894e-5).abs() < 5e-9,
            "mu {}",
            state.dynamic_viscosity
        );
        assert!(
            (state.thermal_conductivity - 2.5326e-2).abs() < 5e-5,
            "kappa {}",
            state.thermal_conductivity
        );
        // Ideal-gas round trip is exact by construction.
        let p_back = state.density * state.specific_gas_constant * state.temperature;
        assert!((p_back - 101_325.0).abs() < 1e-6);
        println!(
            "{{\"suite\":\"fs-material-gas\",\"case\":\"ussa-sea-level\",\"rho\":{:.5},\"c\":{:.3},\"mu\":{:.6e},\"kappa\":{:.6e},\"verdict\":\"pass\"}}",
            state.density, state.sound_speed, state.dynamic_viscosity, state.thermal_conductivity
        );
    }

    #[test]
    fn prandtl_emerges_from_kinetic_theory() {
        // Pr = mu cp / kappa is DERIVED. For USSA dry air at ordinary
        // temperatures the classic ~0.71 must EMERGE — if any of the
        // three ingredient models is wrong, this cross-relation breaks
        // (the built-in falsifier).
        let air = GasSpec::dry_air_ussa1976();
        // Measured emergence (2026-08-08): Pr = 0.7378 / 0.7214 /
        // 0.7099 / 0.6934 at 200 / 250 / 288.15 / 350 K — the derived
        // curve tracks real air's published mild decrease with T, and
        // the classic 0.71 appears exactly at the standard atmosphere's
        // reference temperature.
        for &(t, lo, hi) in &[
            (200.0, 0.72, 0.75),
            (250.0, 0.71, 0.73),
            (288.15, 0.70, 0.72),
            (350.0, 0.68, 0.70),
        ] {
            let state = GasState::try_new(&air, t, 101_325.0).expect("state");
            assert!(
                (lo..hi).contains(&state.prandtl),
                "Pr({t}) = {} outside the emergent band [{lo}, {hi})",
                state.prandtl
            );
        }
        println!(
            "{{\"suite\":\"fs-material-gas\",\"case\":\"prandtl-emerges\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn twenty_celsius_matches_the_hardcoded_acoustic_air_constants() {
        // fs-bem's Medium::air() and fs-couple's AcousticMedium::air()
        // carry (1.204, 343.0) — those constants ARE this primitive at
        // (293.15 K, 101325 Pa), which is the doctrine's point: derive
        // the medium, don't hardcode it.
        let state =
            GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("20 C");
        assert!(
            (state.density - 1.204).abs() < 1.5e-3,
            "rho(20C) {}",
            state.density
        );
        assert!(
            (state.sound_speed - 343.0).abs() < 0.5,
            "c(20C) {}",
            state.sound_speed
        );
        println!(
            "{{\"suite\":\"fs-material-gas\",\"case\":\"acoustic-air-consistency\",\"rho\":{:.4},\"c\":{:.2},\"verdict\":\"pass\"}}",
            state.density, state.sound_speed
        );
    }

    #[test]
    fn eucken_cross_checks_the_ussa_air_fit() {
        // Two INDEPENDENT conductivity routes — the USSA empirical fit
        // and the kinetic-theory Eucken relation — must agree within a
        // few percent for air across the instrument-relevant band.
        // This is a real cross-model check, not a golden.
        let ussa = GasSpec::dry_air_ussa1976();
        let eucken = GasSpec {
            conductivity: ConductivityModel::Eucken,
            ..ussa
        };
        // Measured divergence (2026-08-08): 0.1% / 2% / 4% / 7% / 12%
        // at 200 / 250 / 293 / 400 / 600 K — Eucken increasingly
        // UNDERPREDICTS a diatomic's conductivity as internal modes
        // carry more of the heat, exactly the textbook boundary of the
        // relation. The authored envelope grows with T and the
        // monotone-worsening trend itself is asserted.
        let mut previous_rel = 0.0f64;
        for &(t, cap) in &[
            (200.0, 0.01),
            (250.0, 0.03),
            (293.15, 0.05),
            (400.0, 0.09),
            (600.0, 0.14),
        ] {
            let a = GasState::try_new(&ussa, t, 101_325.0).expect("ussa");
            let b = GasState::try_new(&eucken, t, 101_325.0).expect("eucken");
            let rel =
                (a.thermal_conductivity - b.thermal_conductivity).abs() / a.thermal_conductivity;
            assert!(
                rel < cap,
                "kappa routes at {t} K: USSA {} vs Eucken {} (rel {rel:.3} vs cap {cap})",
                a.thermal_conductivity,
                b.thermal_conductivity
            );
            assert!(
                rel >= previous_rel,
                "the Eucken defect must worsen monotonically with T"
            );
            previous_rel = rel;
        }
        println!(
            "{{\"suite\":\"fs-material-gas\",\"case\":\"eucken-vs-ussa\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn extreme_regime_trends_are_physical() {
        // The doctrine's hot-enclosure regime (lead melts at 600.6 K):
        // evaluation must stay finite with the right monotone physics —
        // heating thins the gas, speeds up sound, and thickens both
        // transport coefficients.
        let air = GasSpec::dry_air_ussa1976();
        let cold = GasState::try_new(&air, 293.15, 101_325.0).expect("cold");
        let hot = GasState::try_new(&air, 700.0, 101_325.0).expect("hot");
        assert!(hot.density < cold.density);
        assert!(hot.sound_speed > cold.sound_speed);
        assert!(hot.dynamic_viscosity > cold.dynamic_viscosity);
        assert!(hot.thermal_conductivity > cold.thermal_conductivity);
        assert!(hot.characteristic_impedance < cold.characteristic_impedance);
        // Pressure at fixed T only moves density/impedance, never the
        // sound speed or transport coefficients (ideal gas).
        let altitude = GasState::try_new(&air, 293.15, 70_000.0).expect("altitude");
        assert_eq!(altitude.sound_speed.to_bits(), cold.sound_speed.to_bits());
        assert_eq!(
            altitude.dynamic_viscosity.to_bits(),
            cold.dynamic_viscosity.to_bits()
        );
        assert!(altitude.density < cold.density);
        println!(
            "{{\"suite\":\"fs-material-gas\",\"case\":\"extreme-regime-trends\",\"c_700k\":{:.1},\"verdict\":\"pass\"}}",
            hot.sound_speed
        );
    }

    #[test]
    fn refusals_fire_and_repeats_are_bitwise() {
        let air = GasSpec::dry_air_ussa1976();
        assert!(GasState::try_new(&air, 20.0, 101_325.0).is_err());
        assert!(GasState::try_new(&air, 293.15, -5.0).is_err());
        assert!(GasState::try_new(&air, 293.15, 2.0e7).is_err());
        assert!(GasState::try_new(&air, f64::NAN, 101_325.0).is_err());
        let bad_gamma = GasSpec { gamma: 1.8, ..air };
        assert!(GasState::try_new(&bad_gamma, 293.15, 101_325.0).is_err());
        let a = GasState::try_new(&air, 311.7, 98_000.0).expect("a");
        let b = GasState::try_new(&air, 311.7, 98_000.0).expect("b");
        assert_eq!(a.prandtl.to_bits(), b.prandtl.to_bits());
        assert_eq!(
            a.thermal_conductivity.to_bits(),
            b.thermal_conductivity.to_bits()
        );
    }
}
