//! ISO 9613-1 atmospheric absorption for *air*.
//!
//! Molecular relaxation of O₂ and N₂ plus the standard's classical
//! term. This is **not** a generic-gas law and it is **not**
//! [`crate::gas::GasState::stokes_kirchhoff_absorption`]: ISO already
//! includes a classical viscothermal piece. Adding the two double-counts.
//!
//! Humidity is an explicit argument in `[0, 1]`, never a hidden
//! scenario field. The standard is for terrestrial air; callers that
//! have a [`GasState`] may pass its `(T, p)` but the relaxation
//! frequencies stay the ISO-9613 air fit.
//!
//! Formula (ISO 9613-1:1993 / Bass–Sutherland–Zuckerwar):
//!
//! `α [Np/m] = f² { 1.84e-11 (T/T₀)^{1/2} (p₀/p)
//!   + (T/T₀)^{-5/2} [ 0.01275 e^{-2239.1/T} / (f_{rO} + f²/f_{rO})
//!                   + 0.1068 e^{-3352/T} / (f_{rN} + f²/f_{rN}) ] }`
//!
//! with molar humidity `h` (%) from relative humidity and the ISO
//! saturation-vapour fit, and relaxation frequencies `f_{rO}`, `f_{rN}`
//! as written in the standard. `T₀ = 293.15 K`, `p₀ = 101325 Pa`.
//!
//! Determinism: `fs_math::det` for every transcendental.

use fs_math::det;

use crate::MaterialError;
use crate::gas::GasState;

/// ISO reference temperature [K].
pub const T0_K: f64 = 293.15;
/// ISO reference pressure [Pa].
pub const P0_PA: f64 = 101_325.0;
/// Water triple point used by the ISO saturation fit [K].
const T01_K: f64 = 273.16;
/// Neper to decibel: `20 / ln(10) ≈ 8.685889`.
const NP_TO_DB: f64 = 8.685889638065037;

/// Geometric spreading of a free-field observer path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Spreading {
    /// No geometric decay (plane wave / guided).
    Plane,
    /// Cylindrical `1/√r` (2-D line source).
    Cylindrical,
    /// Spherical `1/r` (compact 3-D monopole).
    Spherical,
}

/// ISO 9613-1 absorption [Np/m] for air at `(T, p, h_r)`.
///
/// `relative_humidity` is a fraction in `[0, 1]`. `omega` is angular
/// frequency [rad/s].
///
/// # Errors
/// [`MaterialError::Parameters`] when `T`, `p`, humidity, or `ω` are
/// outside the ISO evaluation window (`T` in [200, 400] K — the
/// standard's meteorological band — `p` in (0, 2e5] Pa, humidity in
/// `[0, 1]`, `ω > 0`).
pub fn iso9613_absorption_neper_per_m(
    temperature_k: f64,
    pressure_pa: f64,
    relative_humidity: f64,
    omega: f64,
) -> Result<f64, MaterialError> {
    if !(200.0..=400.0).contains(&temperature_k) || !temperature_k.is_finite() {
        return Err(MaterialError::Parameters {
            what: format!("ISO 9613 temperature {temperature_k} K outside the [200, 400] K band"),
        });
    }
    if !(pressure_pa > 0.0 && pressure_pa <= 2.0e5 && pressure_pa.is_finite()) {
        return Err(MaterialError::Parameters {
            what: format!("ISO 9613 pressure {pressure_pa} Pa outside (0, 2e5] Pa"),
        });
    }
    if !(relative_humidity >= 0.0 && relative_humidity <= 1.0 && relative_humidity.is_finite()) {
        return Err(MaterialError::Parameters {
            what: format!(
                "relative humidity {relative_humidity} must be an explicit fraction in [0, 1]"
            ),
        });
    }
    if !(omega > 0.0 && omega.is_finite()) {
        return Err(MaterialError::Parameters {
            what: format!("ISO 9613 omega {omega} must be positive and finite"),
        });
    }

    let t = temperature_k;
    let t_over_t0 = t / T0_K;
    let p_over_p0 = pressure_pa / P0_PA;
    let freq_hz = omega / (2.0 * core::f64::consts::PI);

    // ISO saturation vapour pressure: psat/p0 = 10^C with
    // C = -6.8346 (T01/T)^1.261 + 4.6151.
    let t01_over_t = T01_K / t;
    let c_sat = -6.8346 * det::exp(1.261 * det::ln(t01_over_t)) + 4.6151;
    let psat_over_p0 = det::exp(c_sat * core::f64::consts::LN_10);
    // Molar concentration of water vapour, *percent*.
    let h = 100.0 * relative_humidity * psat_over_p0 * P0_PA / pressure_pa;

    let fr_o = p_over_p0 * (24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h));
    let t_ratio_cbrt = det::exp(det::ln(t_over_t0) / 3.0);
    let fr_n = p_over_p0
        * det::exp(-0.5 * det::ln(t_over_t0))
        * (9.0 + 280.0 * h * det::exp(-4.170 * (1.0 / t_ratio_cbrt - 1.0)));

    let f2 = freq_hz * freq_hz;
    let x_o = 0.01275 * det::exp(-2239.1 / t) / (fr_o + f2 / fr_o);
    let x_n = 0.1068 * det::exp(-3352.0 / t) / (fr_n + f2 / fr_n);
    let classical = 1.84e-11 * det::exp(0.5 * det::ln(t_over_t0)) / p_over_p0;
    let molecular = det::exp(-2.5 * det::ln(t_over_t0)) * (x_o + x_n);
    Ok(f2 * (classical + molecular))
}

/// Same as [`iso9613_absorption_neper_per_m`] but taking `(T, p)` from
/// a [`GasState`]. The state's transport coefficients are **not**
/// used — ISO's classical term is the air fit, not Stokes–Kirchhoff.
///
/// # Errors
/// Forwards the ISO window refusals.
pub fn iso9613_absorption(
    state: &GasState,
    relative_humidity: f64,
    omega: f64,
) -> Result<f64, MaterialError> {
    iso9613_absorption_neper_per_m(state.temperature, state.pressure, relative_humidity, omega)
}

/// ISO 9613-1 absorption [dB/km] — the unit of the published tables.
///
/// # Errors
/// Forwards [`iso9613_absorption_neper_per_m`].
pub fn iso9613_absorption_db_per_km(
    temperature_k: f64,
    pressure_pa: f64,
    relative_humidity: f64,
    omega: f64,
) -> Result<f64, MaterialError> {
    let np = iso9613_absorption_neper_per_m(temperature_k, pressure_pa, relative_humidity, omega)?;
    Ok(np * NP_TO_DB * 1000.0)
}

/// Free-field amplitude factor `e^{-α r} / r^n` for a compact source.
///
/// Plane: `n = 0`. Cylindrical: `n = 1/2`. Spherical: `n = 1`.
/// `range_m` must be positive. `alpha` is [Np/m].
#[must_use]
pub fn range_factor(range_m: f64, alpha_np_per_m: f64, spreading: Spreading) -> f64 {
    if !(range_m > 0.0 && range_m.is_finite())
        || !(alpha_np_per_m >= 0.0 && alpha_np_per_m.is_finite())
    {
        return 0.0;
    }
    let absorb = det::exp(-alpha_np_per_m * range_m);
    match spreading {
        Spreading::Plane => absorb,
        Spreading::Cylindrical => absorb / det::sqrt(range_m),
        Spreading::Spherical => absorb / range_m,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gas::{GasSpec, GasState};

    #[test]
    fn humid_air_absorbs_more_than_dry_at_one_kilohertz() {
        let omega = 2.0 * core::f64::consts::PI * 1_000.0;
        let dry = iso9613_absorption_neper_per_m(293.15, 101_325.0, 0.0, omega).expect("dry");
        let wet = iso9613_absorption_neper_per_m(293.15, 101_325.0, 0.50, omega).expect("wet");
        assert!(dry > 0.0 && wet > dry, "dry={dry} wet={wet}");
        // Published meteorological tables sit near 5 dB/km at 20 °C,
        // 50 % RH, 1 kHz, 1 atm. Band, not a golden digit.
        let db_km = iso9613_absorption_db_per_km(293.15, 101_325.0, 0.50, omega).expect("db");
        assert!(
            (4.0..7.0).contains(&db_km),
            "20C/50%/1kHz absorption {db_km} dB/km outside the ISO table band"
        );
    }

    #[test]
    fn molecular_term_breaks_pure_omega_squared() {
        let o1 = 2.0 * core::f64::consts::PI * 1_000.0;
        let o4 = 2.0 * core::f64::consts::PI * 4_000.0;
        let a1 = iso9613_absorption_neper_per_m(293.15, 101_325.0, 0.50, o1).expect("1");
        let a4 = iso9613_absorption_neper_per_m(293.15, 101_325.0, 0.50, o4).expect("4");
        // Classical-only would be exactly 16 (ω², 4× frequency).
        // Relaxation knocks the ratio off 16 — that is the molecular
        // signature, not a fitting knob.
        let ratio = a4 / a1;
        assert!(
            (ratio - 16.0).abs() > 1.0,
            "humid ISO ratio {ratio} collapsed onto Stokes ω²"
        );
        assert!(ratio > 2.0 && ratio < 30.0, "ratio {ratio}");
    }

    #[test]
    fn iso_includes_classical_do_not_add_stokes() {
        let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air");
        let omega = 2.0 * core::f64::consts::PI * 1_000.0;
        let iso = iso9613_absorption(&air, 0.0, omega).expect("iso");
        let stokes = air.stokes_kirchhoff_absorption(omega);
        // ISO's classical floor plus O₂/N₂ relaxation sits above
        // Stokes–Kirchhoff. They are alternate packages, not addends.
        assert!(stokes > 0.0);
        assert!(
            iso > stokes,
            "ISO {iso} must include a molecular surplus over Stokes {stokes}"
        );
        assert!(
            iso / stokes < 40.0,
            "ISO/Stokes {} is a surplus, not a different unit",
            iso / stokes
        );
    }

    #[test]
    fn farther_spherical_path_kills_high_frequency_more() {
        let lo = iso9613_absorption_neper_per_m(293.15, 101_325.0, 0.50, 2.0e3).expect("lo");
        let hi = iso9613_absorption_neper_per_m(293.15, 101_325.0, 0.50, 2.0e4).expect("hi");
        let near_hi = range_factor(1.0, hi, Spreading::Spherical);
        let far_hi = range_factor(200.0, hi, Spreading::Spherical);
        let far_lo = range_factor(200.0, lo, Spreading::Spherical);
        let near_lo = range_factor(1.0, lo, Spreading::Spherical);
        let hi_drop = far_hi / near_hi;
        let lo_drop = far_lo / near_lo;
        assert!(
            hi_drop < lo_drop,
            "2 kHz drop {hi_drop} should beat 200 Hz-class drop {lo_drop}"
        );
    }

    #[test]
    fn refusals_are_typed_and_repeats_are_bitwise() {
        let omega = 1.0e3;
        assert!(iso9613_absorption_neper_per_m(100.0, 101_325.0, 0.5, omega).is_err());
        assert!(iso9613_absorption_neper_per_m(293.15, 101_325.0, 1.5, omega).is_err());
        assert!(iso9613_absorption_neper_per_m(293.15, 101_325.0, 0.5, 0.0).is_err());
        let a = iso9613_absorption_neper_per_m(293.15, 101_325.0, 0.37, 4.2e3).expect("a");
        let b = iso9613_absorption_neper_per_m(293.15, 101_325.0, 0.37, 4.2e3).expect("b");
        assert_eq!(a.to_bits(), b.to_bits());
    }
}
