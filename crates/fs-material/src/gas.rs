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
//! and `p` in (0, 1e7] Pa admit evaluation. Accuracy inside the
//! window, stated honestly: the ideal-gas law is good to well under 1%
//! below ~10 bar; toward the 100 bar ceiling air's real-gas
//! corrections reach the 1-3% class in density at 300 K (the
//! first-order virial estimate is Z ~ 0.97, higher virials pull it
//! back toward ~0.99; worse when cold either way). The calorically-perfect
//! `gamma` for diatomic air degrades above ~600 K as vibrational modes
//! excite (order 1-2% in gamma at 600 K, growing with T). The USSA
//! source defines its transport fits for the atmosphere (~187-288 K),
//! and Sutherland-for-air is conventionally quoted over roughly
//! 170-1900 K — outside those bands this module EXTRAPOLATES the fit
//! forms, estimated tier.
//!
//! NO-CLAIM — phase validity is NOT checked: an ideal-gas model cannot
//! know about condensation, so a state in the admitted window can
//! describe a phase that does not exist (air at 60 K and 1 atm is not
//! a gas). Callers own the phase check until vapor-pressure/melting
//! rows from the material-property taxonomy land; the returned state
//! is the gas-phase EXTRAPOLATION, never evidence the gas phase
//! exists there.
//!
//! MOIST AIR (music bead 3ez8g.3.5): [`GasState::try_new_moist_air`]
//! derives the ideal dry-air + water-vapor mixture from relative
//! humidity — a player's warm humid breath vs a cold dry hall is a
//! real ~0.3%-class sound-speed effect, physics rather than a detune
//! knob. Mixture rules, both EXACT for an ideal mixture: molar mass
//! linear in the vapor mole fraction, and the isochoric-heat identity
//! `1/(γ_mix−1) = Σ x_i/(γ_i−1)`. The vapor fraction comes from the
//! Buck 1996 saturation fit over liquid water
//! (`e_s = 611.21 exp((18.678 − t/234.5) t/(257.14 + t))` Pa, t in
//! °C; quoted for −20..+50 °C — refused outside when RH > 0). Water
//! vapor spec provenance: `M = 18.01528e-3` kg/mol (CODATA/IUPAC),
//! `γ = 1.3291` from NIST-JANAF `cp(H2O g, 298.15 K) = 33.58`
//! J/(mol K) via `γ = cp/(cp − R)`, Sutherland `β = 2.418e-6`,
//! `S = 1064 K` (White, *Viscous Fluid Flow*, steam constants;
//! reproduces `μ(373 K) ≈ 1.21e-5` Pa s).
//!
//! MOIST-AIR NO-CLAIMS, stated honestly: (a) transport coefficients
//! (μ, κ) REMAIN THE DRY-AIR FITS — humidity's transport effect is
//! sub-1% at musical vapor fractions and ~2%-class at the admitted
//! ceiling `x_w ≤ 0.15` (refused above: the disclosed approximation
//! band, estimated tier). (b) The saturation ENHANCEMENT FACTOR of
//! moist air over the pure-phase `e_s` (~0.5% at 1 atm) is neglected.
//! (c) RH > 1 is REFUSED — a supersaturated input describes a
//! condensing state the ideal-gas model cannot represent, so the
//! phase no-claim is structural at this constructor's input; the
//! plain [`GasState::try_new`] window no-claim above still stands.
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

    /// Water vapor: `M = 18.01528e-3` kg/mol (CODATA/IUPAC),
    /// `gamma = 1.3291` from NIST-JANAF `cp(H2O g, 298.15 K) = 33.58`
    /// J/(mol K) via `gamma = cp/(cp - R)`, Sutherland steam constants
    /// `beta = 2.418e-6`, `S = 1064 K` (White, *Viscous Fluid Flow*;
    /// gives `mu(373 K) ~ 1.21e-5` Pa s vs the tabulated ~1.2e-5),
    /// Eucken conductivity.
    #[must_use]
    pub const fn water_vapor_nist() -> GasSpec {
        GasSpec {
            molar_mass: 18.01528e-3,
            gamma: 1.3291,
            sutherland_beta: 2.418e-6,
            sutherland_s: 1064.0,
            conductivity: ConductivityModel::Eucken,
        }
    }
}

/// Saturation vapor pressure of water over the LIQUID phase [Pa] —
/// the Buck 1996 fit `e_s = 611.21 exp((18.678 - t/234.5) t /
/// (257.14 + t))` with `t` in Celsius, quoted for −20..+50 °C
/// (253.15..323.15 K). Named provenance: A. L. Buck, *New equations
/// for computing vapor pressure and enhancement factor* (J. Appl.
/// Meteorol. 20, 1981; revised constants 1996). The pure-phase value:
/// the moist-air enhancement factor (~0.5% at 1 atm) is NOT applied.
///
/// # Errors
/// [`MaterialError::Parameters`] outside the fit's quoted window.
pub fn saturation_pressure_water_pa(temperature: f64) -> Result<f64, MaterialError> {
    if !((253.15..=323.15).contains(&temperature)) {
        return Err(MaterialError::Parameters {
            what: format!(
                "temperature {temperature} K outside the Buck-1996 liquid-water \
                 window [253.15, 323.15] K (−20..+50 °C)"
            ),
        });
    }
    let t = temperature - 273.15;
    Ok(611.21 * det::exp((18.678 - t / 234.5) * t / (257.14 + t)))
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
    /// Water-vapor mole fraction `x_w` (0 for every dry
    /// construction; set by [`GasState::try_new_moist_air`]).
    pub water_mole_fraction: f64,
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
            water_mole_fraction: 0.0,
        })
    }

    /// Moist air from `(T, p, relative humidity)`: the ideal dry-air +
    /// water-vapor mixture (see the module doc's MOIST AIR section for
    /// provenance and no-claims). Both mixture rules are exact for an
    /// ideal mixture: `M_mix = M_a + x_w (M_w − M_a)` and the
    /// isochoric-heat identity `1/(γ_mix−1) = (1−x_w)/(γ_a−1) +
    /// x_w/(γ_v−1)`. Transport coefficients remain the DRY-AIR fits
    /// (disclosed, estimated tier, `x_w ≤ 0.15` refused above).
    ///
    /// `relative_humidity == 0` takes the plain dry-air path — the dry
    /// limit is the SAME CODE, bitwise.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] for RH outside `[0, 1]`
    /// (supersaturation is a condensing state the ideal-gas model
    /// cannot represent — the structural phase refusal), a temperature
    /// outside the Buck fit's window when RH > 0, a vapor fraction
    /// above the disclosed 0.15 approximation ceiling, or the plain
    /// [`GasState::try_new`] window refusals.
    pub fn try_new_moist_air(
        temperature: f64,
        pressure: f64,
        relative_humidity: f64,
    ) -> Result<Self, MaterialError> {
        if !(relative_humidity.is_finite() && (0.0..=1.0).contains(&relative_humidity)) {
            return Err(MaterialError::Parameters {
                what: format!(
                    "relative humidity {relative_humidity} outside [0, 1] (RH > 1 is a \
                     condensing state the ideal-gas mixture cannot represent)"
                ),
            });
        }
        let dry = GasSpec::dry_air_ussa1976();
        if relative_humidity == 0.0 {
            return Self::try_new(&dry, temperature, pressure);
        }
        let vapor = GasSpec::water_vapor_nist();
        let e_sat = saturation_pressure_water_pa(temperature)?;
        if !(pressure > 0.0 && pressure.is_finite()) {
            return Err(MaterialError::Parameters {
                what: format!("pressure {pressure} Pa must be positive and finite"),
            });
        }
        let x_w = relative_humidity * e_sat / pressure;
        if x_w > 0.15 {
            return Err(MaterialError::Parameters {
                what: format!(
                    "water mole fraction {x_w:.4} above the 0.15 ceiling of the \
                     dry-air-transport approximation (disclosed band)"
                ),
            });
        }
        let molar_mass = dry.molar_mass + x_w * (vapor.molar_mass - dry.molar_mass);
        let inv_gm1 = (1.0 - x_w) / (dry.gamma - 1.0) + x_w / (vapor.gamma - 1.0);
        let gamma = 1.0 + 1.0 / inv_gm1;
        let mixture = GasSpec {
            molar_mass,
            gamma,
            ..dry
        };
        let mut state = Self::try_new(&mixture, temperature, pressure)?;
        state.water_mole_fraction = x_w;
        Ok(state)
    }

    /// Classical Stokes–Kirchhoff absorption coefficient [1/m].
    ///
    /// `α = ω² / (2 ρ c³) · [4μ/3 + (γ−1) κ / c_p]`. This is the
    /// viscothermal bulk law for a calorically perfect gas, not ISO
    /// 9613 humidity relaxation and not a measured outdoor curve.
    #[must_use]
    pub fn stokes_kirchhoff_absorption(&self, omega: f64) -> f64 {
        if !(omega.is_finite() && omega > 0.0) {
            return 0.0;
        }
        let shear = 4.0 / 3.0 * self.dynamic_viscosity;
        let thermal = (self.gamma - 1.0) * self.thermal_conductivity / self.specific_heat_cp;
        let num = omega * omega * (shear + thermal);
        let den = 2.0 * self.density * self.sound_speed * self.sound_speed * self.sound_speed;
        if den > 0.0 { num / den } else { 0.0 }
    }

    /// ISO 9613-1 absorption [Np/m] at this `(T, p)` and an explicit
    /// relative humidity in `[0, 1]`.
    ///
    /// Transport coefficients on this state are **not** used — ISO's
    /// classical term is the air fit. Do not add
    /// [`Self::stokes_kirchhoff_absorption`] on top.
    ///
    /// # Errors
    /// Forwards the ISO meteorological-window refusals.
    pub fn iso9613_absorption(
        &self,
        relative_humidity: f64,
        omega: f64,
    ) -> Result<f64, crate::MaterialError> {
        crate::iso9613::iso9613_absorption(self, relative_humidity, omega)
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
        let a1 = state.stokes_kirchhoff_absorption(1.0e3);
        let a2 = state.stokes_kirchhoff_absorption(2.0e3);
        assert!(a1 > 0.0 && (a2 / a1 - 4.0).abs() < 1.0e-9);
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
    fn monatomic_eucken_prandtl_is_exactly_two_thirds() {
        // Absolute pin on the Eucken coefficient (review round 3: the
        // relative USSA cross-check tolerances would let a ~1% slip in
        // the 5/4 constant survive). Under Eucken, Pr = cp/(cp + 5R/4)
        // = 4 gamma / (9 gamma - 5): mu cancels, so for a monatomic
        // gas (gamma = 5/3) Pr = 2/3 EXACTLY, independent of the
        // Sutherland constants, T, and p. Any perturbation of the 1.25
        // coefficient breaks this to first order.
        let monatomic = GasSpec {
            molar_mass: 39.948e-3, // argon
            gamma: 5.0 / 3.0,
            sutherland_beta: 1.93e-6,
            sutherland_s: 142.0,
            conductivity: ConductivityModel::Eucken,
        };
        for &t in &[100.0, 293.15, 1200.0] {
            let state = GasState::try_new(&monatomic, t, 101_325.0).expect("argon");
            assert!(
                (state.prandtl - 2.0 / 3.0).abs() < 1e-14,
                "monatomic Eucken Pr({t}) = {} must be exactly 2/3",
                state.prandtl
            );
        }
        println!(
            "{{\"suite\":\"fs-material-gas\",\"case\":\"monatomic-eucken-exact\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn moist_saturation_fit_two_sources_and_published_pins() {
        // TWO-SOURCE RULE on the RH -> vapor step: the implemented Buck
        // 1996 fit against the INDEPENDENT Magnus fit with the
        // Alduchov–Eskridge 1996 constants (e_s = 610.94 exp(17.625 t /
        // (243.04 + t))) across the full quoted window, plus published
        // absolute pins: IAPWS gives e_s(0.01 C) = 611.657 Pa and
        // e_s(20 C) = 2339.2 Pa.
        // (-20 C exactly is 253.14999... in floats — the window edges
        // are asserted in Kelvin below instead.)
        for t_c in -19..=50 {
            let t = f64::from(t_c);
            let buck = saturation_pressure_water_pa(t + 273.15).expect("in window");
            let magnus = 610.94 * det::exp(17.625 * t / (243.04 + t));
            let rel = ((buck - magnus) / magnus).abs();
            assert!(
                rel < 5.0e-3,
                "Buck vs Magnus at {t} C: {buck:.2} vs {magnus:.2} (rel {rel:.2e})"
            );
        }
        let e0 = saturation_pressure_water_pa(273.16).expect("triple point");
        assert!((608.0..615.0).contains(&e0), "e_s(0.01 C) = {e0:.2} Pa");
        let e20 = saturation_pressure_water_pa(293.15).expect("20 C");
        assert!((2329.0..2349.0).contains(&e20), "e_s(20 C) = {e20:.1} Pa");
        assert!(saturation_pressure_water_pa(253.15).is_ok());
        assert!(saturation_pressure_water_pa(323.15).is_ok());
        assert!(saturation_pressure_water_pa(250.0).is_err());
        assert!(saturation_pressure_water_pa(330.0).is_err());
        println!(
            "{{\"suite\":\"fs-material-gas\",\"case\":\"moist-sat-two-sources\",\"e0\":{e0:.2},\"e20\":{e20:.1},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn moist_density_hits_the_oiml_conventional_anchor() {
        // The metrology community's CONVENTIONAL air density (OIML
        // R 111 / CIPM): rho = 1.2 kg/m^3 is DEFINED at exactly
        // (20 C, 101325 Pa, 50% RH). An independent published anchor
        // for the whole mixture chain — our ideal mixture (no CO2, no
        // real-gas Z, no enhancement factor) must land within 0.2%.
        // And the classic counterintuitive: humid air is LIGHTER.
        let moist = GasState::try_new_moist_air(293.15, 101_325.0, 0.5).expect("moist");
        let dry = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("dry");
        assert!(
            (moist.density - 1.2).abs() < 2.4e-3,
            "rho(20C, 50% RH) = {} vs the OIML conventional 1.2",
            moist.density
        );
        assert!(
            moist.density < dry.density,
            "humid air must be LIGHTER ({} vs {})",
            moist.density,
            dry.density
        );
        println!(
            "{{\"suite\":\"fs-material-gas\",\"case\":\"moist-oiml-anchor\",\"rho\":{:.5},\"x_w\":{:.5},\"verdict\":\"pass\"}}",
            moist.density, moist.water_mole_fraction
        );
    }

    #[test]
    fn moist_dry_limit_is_bitwise() {
        // DONE-WHEN clause: RH = 0 degenerates to the existing dry
        // path EXACTLY (same code path, so every field is bit-equal
        // and no golden anywhere can move).
        let via_moist = GasState::try_new_moist_air(288.15, 101_325.0, 0.0).expect("rh0");
        let dry = GasState::try_new(&GasSpec::dry_air_ussa1976(), 288.15, 101_325.0).expect("dry");
        for (a, b, what) in [
            (via_moist.density, dry.density, "rho"),
            (via_moist.sound_speed, dry.sound_speed, "c"),
            (via_moist.dynamic_viscosity, dry.dynamic_viscosity, "mu"),
            (
                via_moist.thermal_conductivity,
                dry.thermal_conductivity,
                "kappa",
            ),
            (via_moist.gamma, dry.gamma, "gamma"),
            (via_moist.prandtl, dry.prandtl, "Pr"),
        ] {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "dry limit must be bitwise in {what}"
            );
        }
        println!(
            "{{\"suite\":\"fs-material-gas\",\"case\":\"moist-dry-limit-bitwise\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn moist_mixture_rules_match_independent_expressions() {
        // TAUTOLOGY GUARD (the bead's recorded lesson): never compare
        // two quantities derived from the same GasState. Here the
        // molar mass and gamma are recomputed IN THE TEST from raw
        // literals and the public saturation fn, then matched against
        // what the constructor derived — the mu-cancelling identity
        // route for gamma (cv mixing) vs the constructor's own path.
        let (t, p, rh) = (303.15, 98_000.0, 0.73);
        let state = GasState::try_new_moist_air(t, p, rh).expect("moist");
        let x = rh * saturation_pressure_water_pa(t).expect("es") / p;
        assert!(
            (state.water_mole_fraction - x).abs() < 1e-15,
            "x_w must be recorded on the state"
        );
        let m_indep = 28.9644e-3 * (1.0 - x) + 18.01528e-3 * x;
        let m_state = R_USSA_1976 / state.specific_gas_constant;
        assert!(
            ((m_state - m_indep) / m_indep).abs() < 1e-12,
            "mixture molar mass: state {m_state:.9e} vs independent {m_indep:.9e}"
        );
        let gamma_indep = 1.0 + 1.0 / ((1.0 - x) / 0.4 + x / 0.3291);
        assert!(
            ((state.gamma - gamma_indep) / gamma_indep).abs() < 1e-12,
            "mixture gamma: state {} vs independent {gamma_indep}",
            state.gamma
        );
        // The humidity direction on both audible quantities, at fixed
        // (T, p): faster sound, lighter air — monotonically in RH.
        let mut last_c = 0.0;
        let mut last_rho = f64::INFINITY;
        for rh_step in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let s = GasState::try_new_moist_air(293.15, 101_325.0, rh_step).expect("sweep");
            assert!(s.sound_speed > last_c, "c must rise with RH");
            assert!(s.density < last_rho, "rho must fall with RH");
            last_c = s.sound_speed;
            last_rho = s.density;
        }
        // The published magnitude class (Wong 1986 / Cramer 1993:
        // saturation at 20 C raises c by ~0.35%): authored Estimate
        // band 0.8..1.8 m/s.
        let dry = GasState::try_new_moist_air(293.15, 101_325.0, 0.0).expect("dry");
        let wet = GasState::try_new_moist_air(293.15, 101_325.0, 1.0).expect("wet");
        let uplift = wet.sound_speed - dry.sound_speed;
        assert!(
            (0.8..1.8).contains(&uplift),
            "saturated 20 C uplift {uplift:.3} m/s outside the published class"
        );
        println!(
            "{{\"suite\":\"fs-material-gas\",\"case\":\"moist-mixture-rules\",\"gamma\":{:.6},\"uplift_ms\":{uplift:.3},\"verdict\":\"pass\"}}",
            state.gamma
        );
    }

    #[test]
    fn moist_refusals_fire_and_repeats_are_bitwise() {
        // RH outside [0,1] (the structural phase refusal), the Buck
        // window when RH > 0 (while RH = 0 keeps the full dry window),
        // and the vapor-fraction ceiling of the dry-transport
        // approximation.
        assert!(GasState::try_new_moist_air(293.15, 101_325.0, -0.1).is_err());
        assert!(GasState::try_new_moist_air(293.15, 101_325.0, 1.1).is_err());
        assert!(GasState::try_new_moist_air(293.15, 101_325.0, f64::NAN).is_err());
        assert!(
            GasState::try_new_moist_air(240.0, 101_325.0, 0.3).is_err(),
            "240 K with RH > 0 is outside the Buck window"
        );
        assert!(
            GasState::try_new_moist_air(240.0, 101_325.0, 0.0).is_ok(),
            "RH = 0 keeps the full dry validity window"
        );
        // 50 C saturated at 0.6 atm: x_w = e_s/p ~ 0.20 > 0.15.
        assert!(
            GasState::try_new_moist_air(323.15, 60_000.0, 1.0).is_err(),
            "the vapor-fraction ceiling must refuse"
        );
        let a = GasState::try_new_moist_air(310.15, 99_000.0, 0.9).expect("a");
        let b = GasState::try_new_moist_air(310.15, 99_000.0, 0.9).expect("b");
        assert_eq!(a.sound_speed.to_bits(), b.sound_speed.to_bits());
        assert_eq!(a.density.to_bits(), b.density.to_bits());
        assert_eq!(a.gamma.to_bits(), b.gamma.to_bits());
        println!(
            "{{\"suite\":\"fs-material-gas\",\"case\":\"moist-refusals\",\"verdict\":\"pass\"}}"
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
