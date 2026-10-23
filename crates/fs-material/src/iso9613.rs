//! ISO 9613-1 atmospheric absorption for *air*.
//!
//! Molecular relaxation of O₂ and N₂ plus the standard's classical
//! term. This is **not** a generic-gas law and it is **not**
//! [`crate::gas::GasState::stokes_kirchhoff_absorption`]: ISO already
//! includes a classical viscothermal piece. Adding the two double-counts.
//!
//! Humidity is an explicit argument in `[0, 1]`, never a hidden
//! scenario field. The standard is for terrestrial air. A [`GasState`] is
//! evaluated only with a [`TerrestrialAir`] admission of its explicit
//! composition, so heliox or a CO2-rich mixture refuses (naming the species)
//! instead of receiving air's absorption.
//!
//! Formula (ISO 9613-1:1993 / Bass–Sutherland–Zuckerwar):
//!
//! ```text
//! α [Np/m] = f² { 1.84e-11 (T/T₀)^{1/2} (p₀/p)
//!   + (T/T₀)^{-5/2} [ 0.01275 e^{-2239.1/T} / (f_{rO} + f²/f_{rO})
//!                   + 0.1068 e^{-3352/T} / (f_{rN} + f²/f_{rN}) ] }
//! ```
//!
//! with molar humidity `h` (%) from relative humidity and the ISO
//! saturation-vapour fit, and relaxation frequencies `f_{rO}`, `f_{rN}`
//! as written in the standard. `T₀ = 293.15 K`, `p₀ = 101325 Pa`.
//!
//! Determinism: `fs_math::det` for every transcendental.

use fs_math::det;

use fs_thermochem::SpeciesId;

use crate::MaterialError;
use crate::gas::GasState;
use crate::relaxation::RelaxingMedium;

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
/// [`MaterialError::State`] when those inputs exceed the calculation's finite
/// representable range. No silent absorption or unit-conversion fallback is used.
pub fn iso9613_absorption_neper_per_m(
    temperature_k: f64,
    pressure_pa: f64,
    relative_humidity: f64,
    omega: f64,
) -> Result<f64, MaterialError> {
    admit_temperature_pressure(temperature_k, pressure_pa)?;
    if !((0.0..=1.0).contains(&relative_humidity) && relative_humidity.is_finite()) {
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
    let h = water_percent(t, pressure_pa, relative_humidity);

    let (fr_o, fr_n) = relaxation_frequencies(t_over_t0, p_over_p0, h);

    let f2 = freq_hz * freq_hz;
    let x_o = 0.01275 * det::exp(-2239.1 / t) / (fr_o + f2 / fr_o);
    let x_n = 0.1068 * det::exp(-3352.0 / t) / (fr_n + f2 / fr_n);
    let classical = 1.84e-11 * det::exp(0.5 * det::ln(t_over_t0)) / p_over_p0;
    let molecular = det::exp(-2.5 * det::ln(t_over_t0)) * (x_o + x_n);
    finite_absorption(f2 * (classical + molecular), "Np/m")
}

// ISO's molar water concentration in percent, from its saturation fit
// psat/p0 = 10^C with C = -6.8346 (T01/T)^1.261 + 4.6151.
fn water_percent(temperature_k: f64, pressure_pa: f64, relative_humidity: f64) -> f64 {
    let t01_over_t = T01_K / temperature_k;
    let c_sat = -6.8346 * det::exp(1.261 * det::ln(t01_over_t)) + 4.6151;
    let psat_over_p0 = det::exp(c_sat * core::f64::consts::LN_10);
    100.0 * relative_humidity * psat_over_p0 * P0_PA / pressure_pa
}

/// ISO 9613-1's absorption at `(T, p, h_r)` as its own terms:
/// `alpha(f) = classical f^2 + sum_j strength_j f^2 / (f_r,j + f^2 / f_r,j)`,
/// O2 first, then N2.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Iso9613Terms {
    /// `1.84e-11 (T/T0)^{1/2} (p0/p)` [Np/(m Hz^2)].
    pub classical_np_per_m_hz2: f64,
    /// Per relaxing species: `(T/T0)^{-5/2} a_j e^{-theta_j/T}` [Np/(m Hz)]
    /// with ISO's `a_O = 0.01275`, `a_N = 0.1068`, and `f_r,j` [Hz].
    pub relaxation: [(f64, f64); 2],
}

/// [`Iso9613Terms`] at `(T, p, h_r)`, with the same window as
/// [`iso9613_absorption_neper_per_m`].
///
/// # Errors
/// As [`iso9613_absorption_neper_per_m`], without the frequency.
pub fn iso9613_terms(
    temperature_k: f64,
    pressure_pa: f64,
    relative_humidity: f64,
) -> Result<Iso9613Terms, MaterialError> {
    iso9613_absorption_neper_per_m(temperature_k, pressure_pa, relative_humidity, 1.0)?;
    let t = temperature_k;
    let t_over_t0 = t / T0_K;
    let (fr_o, fr_n) = relaxation_frequencies(
        t_over_t0,
        pressure_pa / P0_PA,
        water_percent(t, pressure_pa, relative_humidity),
    );
    let scale = det::exp(-2.5 * det::ln(t_over_t0));
    Ok(Iso9613Terms {
        classical_np_per_m_hz2: 1.84e-11 * det::exp(0.5 * det::ln(t_over_t0))
            / (pressure_pa / P0_PA),
        relaxation: [
            (scale * 0.01275 * det::exp(-2239.1 / t), fr_o),
            (scale * 0.1068 * det::exp(-3352.0 / t), fr_n),
        ],
    })
}

/// ISO's relative-humidity coordinate of a gas whose water mole fraction is
/// `water_mole_fraction`, by ISO's own saturation fit: the inverse of the
/// humidity-to-water conversion inside [`iso9613_absorption_neper_per_m`].
///
/// # Errors
/// [`MaterialError::Parameters`] outside ISO's window, or for a water fraction
/// that is negative or above ISO's saturation.
pub fn iso9613_relative_humidity(
    temperature_k: f64,
    pressure_pa: f64,
    water_mole_fraction: f64,
) -> Result<f64, MaterialError> {
    admit_temperature_pressure(temperature_k, pressure_pa)?;
    let relative_humidity =
        100.0 * water_mole_fraction / water_percent(temperature_k, pressure_pa, 1.0);
    if !(relative_humidity.is_finite() && (0.0..=1.0).contains(&relative_humidity)) {
        return Err(MaterialError::Parameters {
            what: format!(
                "water mole fraction {water_mole_fraction} is not between dry and ISO 9613's saturation"
            ),
        });
    }
    Ok(relative_humidity)
}

fn admit_temperature_pressure(temperature_k: f64, pressure_pa: f64) -> Result<(), MaterialError> {
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
    Ok(())
}

// The standard's empirical relaxation-frequency laws, shared bit-for-bit by
// the ISO absorption formula and the generic relaxation law's ISO rate variant.
fn relaxation_frequencies(t_over_t0: f64, p_over_p0: f64, h: f64) -> (f64, f64) {
    let fr_o = p_over_p0 * (24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h));
    let t_ratio_cbrt = det::exp(det::ln(t_over_t0) / 3.0);
    let fr_n = p_over_p0
        * det::exp(-0.5 * det::ln(t_over_t0))
        * (9.0 + 280.0 * h * det::exp(-4.170 * (1.0 / t_ratio_cbrt - 1.0)));
    (fr_o, fr_n)
}

/// ISO 9613-1 relaxation frequencies `(f_rO, f_rN)` [Hz] of oxygen and
/// nitrogen at `(T, p)` and an explicit molar water concentration `h` in
/// percent (`h = 100 x_H2O`, not relative humidity):
///
/// ```text
/// f_rO = (p/p0) (24 + 4.04e4 h (0.02 + h) / (0.391 + h))
/// f_rN = (p/p0) (T/T0)^{-1/2} (9 + 280 h exp(-4.170 ((T/T0)^{-1/3} - 1)))
/// ```
///
/// These are the standard's empirical fits for terrestrial air, in which dry
/// air's collision partners (N2, O2, Ar, trace CO2) are implicit. They are not
/// a collision-partner law and say nothing about any other mixture.
///
/// # Errors
/// [`MaterialError::Parameters`] for `T` outside [200, 400] K, `p` outside
/// (0, 2e5] Pa, or a negative/non-finite water concentration.
pub fn iso9613_relaxation_frequencies_hz(
    temperature_k: f64,
    pressure_pa: f64,
    water_molar_percent: f64,
) -> Result<(f64, f64), MaterialError> {
    admit_temperature_pressure(temperature_k, pressure_pa)?;
    if !(water_molar_percent.is_finite() && (0.0..=100.0).contains(&water_molar_percent)) {
        return Err(MaterialError::Parameters {
            what: format!(
                "ISO 9613 molar water concentration {water_molar_percent} % must be in [0, 100]"
            ),
        });
    }
    let (fr_o, fr_n) = relaxation_frequencies(
        temperature_k / T0_K,
        pressure_pa / P0_PA,
        water_molar_percent,
    );
    if !(fr_o.is_finite() && fr_o > 0.0 && fr_n.is_finite() && fr_n > 0.0) {
        return Err(MaterialError::State {
            what: "ISO 9613 relaxation frequencies are not finite and positive".into(),
        });
    }
    Ok((fr_o, fr_n))
}

fn finite_absorption(value: f64, unit: &str) -> Result<f64, MaterialError> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(MaterialError::State {
            what: format!("ISO 9613 absorption is not finite and nonnegative in {unit}"),
        })
    }
}

/// Admission policy for terrestrial air: dry-basis O2 and N2 mole fractions
/// within this relative distance of ISO's 0.209 and 0.781. It identifies
/// terrestrial air; it is not an accuracy statement.
pub const ISO9613_DRY_AIR_RELATIVE_TOLERANCE: f64 = 0.01;
/// Admission policy for terrestrial air: dry-basis CO2 ceiling.
pub const ISO9613_MAX_DRY_CO2: f64 = 1.0e-3;
/// ISO 9613-1's dry-air oxygen mole fraction.
pub const ISO9613_DRY_O2: f64 = 0.209;
/// ISO 9613-1's dry-air nitrogen mole fraction.
pub const ISO9613_DRY_N2: f64 = 0.781;

// U.S. Standard Atmosphere 1976, Table 3: the dry composition behind
// `GasSpec::dry_air_ussa1976`, and so behind every legacy (T, p, RH) record.
const USSA1976_DRY_O2: f64 = 0.209_476;
const USSA1976_DRY_N2: f64 = 0.780_84;
const USSA1976_DRY_CO2: f64 = 0.000_314;

/// Why ISO 9613-1 refused a gas. ISO 9613-1 is a law for terrestrial air.
/// It is applied only to a gas whose explicit composition was admitted.
#[derive(Debug, Clone, PartialEq)]
pub enum Iso9613Refusal {
    /// A species the standard does not model is present.
    ForeignSpecies {
        /// The species.
        species: SpeciesId,
        /// Its mole fraction.
        mole_fraction: f64,
    },
    /// O2 or N2 is off ISO's dry-basis fraction by more than the tolerance.
    OffStandard {
        /// `"O2"` or `"N2"`.
        species: &'static str,
        /// The gas's dry-basis mole fraction.
        dry_mole_fraction: f64,
        /// ISO's dry-basis mole fraction.
        iso_dry_mole_fraction: f64,
        /// [`ISO9613_DRY_AIR_RELATIVE_TOLERANCE`].
        relative_tolerance: f64,
    },
    /// CO2 above the admitted dry-basis ceiling.
    AboveCeiling {
        /// `"CO2"`.
        species: &'static str,
        /// The gas's dry-basis mole fraction.
        dry_mole_fraction: f64,
        /// [`ISO9613_MAX_DRY_CO2`].
        ceiling: f64,
    },
    /// The gas is all water vapour.
    NoDryAir,
    /// The admission was made for a different gas state.
    DifferentState,
}

impl Iso9613Refusal {
    /// The species the refusal names, if any.
    #[must_use]
    pub fn species(&self) -> Option<&str> {
        match self {
            Self::ForeignSpecies { species, .. } => Some(species.as_str()),
            Self::OffStandard { species, .. } | Self::AboveCeiling { species, .. } => {
                Some(*species)
            }
            Self::NoDryAir | Self::DifferentState => None,
        }
    }
}

impl core::fmt::Display for Iso9613Refusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ForeignSpecies {
                species,
                mole_fraction,
            } => write!(
                f,
                "ISO 9613-1 is for terrestrial air only; {species} is present (mole fraction {mole_fraction})"
            ),
            Self::OffStandard {
                species,
                dry_mole_fraction,
                iso_dry_mole_fraction,
                relative_tolerance,
            } => write!(
                f,
                "ISO 9613-1 is for terrestrial air only; dry-basis {species} is {dry_mole_fraction}, not ISO's {iso_dry_mole_fraction} within {relative_tolerance} relative"
            ),
            Self::AboveCeiling {
                species,
                dry_mole_fraction,
                ceiling,
            } => write!(
                f,
                "ISO 9613-1 is for terrestrial air only; dry-basis {species} {dry_mole_fraction} exceeds the {ceiling} ceiling"
            ),
            Self::NoDryAir => {
                f.write_str("ISO 9613-1 is for terrestrial air only; the gas has no dry air")
            }
            Self::DifferentState => f.write_str(
                "ISO 9613-1 terrestrial-air admission was made for a different gas state",
            ),
        }
    }
}

impl std::error::Error for Iso9613Refusal {}

/// The dry-basis air check shared by [`TerrestrialAir`] and the relaxation
/// law's ISO rate variant. Returns the dry-basis (O2, N2, CO2) fractions.
pub(crate) fn admit_dry_air(
    o2: f64,
    n2: f64,
    co2: f64,
    h2o: f64,
) -> Result<[f64; 3], Iso9613Refusal> {
    let dry = 1.0 - h2o;
    if !(dry > 0.0) {
        return Err(Iso9613Refusal::NoDryAir);
    }
    let fractions = [o2 / dry, n2 / dry, co2 / dry];
    for (species, value, iso) in [
        ("O2", fractions[0], ISO9613_DRY_O2),
        ("N2", fractions[1], ISO9613_DRY_N2),
    ] {
        if !((value / iso - 1.0).abs() <= ISO9613_DRY_AIR_RELATIVE_TOLERANCE) {
            return Err(Iso9613Refusal::OffStandard {
                species,
                dry_mole_fraction: value,
                iso_dry_mole_fraction: iso,
                relative_tolerance: ISO9613_DRY_AIR_RELATIVE_TOLERANCE,
            });
        }
    }
    if !(fractions[2] <= ISO9613_MAX_DRY_CO2) {
        return Err(Iso9613Refusal::AboveCeiling {
            species: "CO2",
            dry_mole_fraction: fractions[2],
            ceiling: ISO9613_MAX_DRY_CO2,
        });
    }
    Ok(fractions)
}

/// Admit an explicit composition as terrestrial air. Species are identified
/// by their canonical ids N2, O2, Ar, CO2 and H2O; zero fractions are ignored.
/// Returns the dry-basis (O2, N2, CO2) fractions and the water mole fraction.
pub(crate) fn admit_composition(
    components: &[(SpeciesId, f64)],
) -> Result<([f64; 3], f64), Iso9613Refusal> {
    let (mut n2, mut o2, mut co2, mut h2o) = (0.0, 0.0, 0.0, 0.0);
    for (species, fraction) in components {
        if *fraction == 0.0 {
            continue;
        }
        match species.as_str() {
            "N2" => n2 = *fraction,
            "O2" => o2 = *fraction,
            "CO2" => co2 = *fraction,
            "H2O" => h2o = *fraction,
            "Ar" => {}
            _ => {
                return Err(Iso9613Refusal::ForeignSpecies {
                    species: species.clone(),
                    mole_fraction: *fraction,
                });
            }
        }
    }
    Ok((admit_dry_air(o2, n2, co2, h2o)?, h2o))
}

/// A gas admitted as terrestrial air, the only gas ISO 9613-1 describes.
///
/// The witness can be made in two ways:
/// - from the gas's own explicit composition ([`Self::admit`]);
/// - together with the legacy moist-air state it describes ([`Self::moist_air`]).
///
/// It is bound to that exact gas state, so it does not transfer to another
/// gas. It admits the dry composition only. Humidity stays ISO's own
/// relative-humidity coordinate, declared by the caller.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TerrestrialAir {
    state: GasState,
    dry: [f64; 3],
}

impl TerrestrialAir {
    /// Admit a gas by its explicit composition.
    ///
    /// # Errors
    /// The [`Iso9613Refusal`] naming the first species that is not
    /// terrestrial air: a foreign species (in species-id order), then O2, N2
    /// and CO2.
    pub fn admit(medium: &RelaxingMedium) -> Result<Self, Iso9613Refusal> {
        let (dry, _) = admit_composition(medium.components())?;
        Ok(Self {
            state: *medium.state(),
            dry,
        })
    }

    /// The legacy moist-air state at `(T, p, RH)`
    /// ([`GasState::try_new_moist_air`]), with its witness. Its dry gas is
    /// U.S. Standard Atmosphere 1976 air.
    ///
    /// # Errors
    /// The state's own refusals.
    pub fn moist_air(
        temperature_k: f64,
        pressure_pa: f64,
        relative_humidity: f64,
    ) -> Result<(GasState, Self), MaterialError> {
        let state = GasState::try_new_moist_air(temperature_k, pressure_pa, relative_humidity)?;
        let dry = admit_dry_air(USSA1976_DRY_O2, USSA1976_DRY_N2, USSA1976_DRY_CO2, 0.0).map_err(
            |refusal| MaterialError::State {
                what: refusal.to_string(),
            },
        )?;
        Ok((state, Self { state, dry }))
    }

    /// Whether this admission was made for exactly `state` (every coefficient).
    #[must_use]
    pub fn admits(&self, state: &GasState) -> bool {
        self.state == *state
    }

    /// [`Self::admits`] as a refusal.
    ///
    /// # Errors
    /// [`Iso9613Refusal::DifferentState`].
    pub fn check(&self, state: &GasState) -> Result<(), Iso9613Refusal> {
        if self.admits(state) {
            Ok(())
        } else {
            Err(Iso9613Refusal::DifferentState)
        }
    }

    /// Dry-basis (O2, N2, CO2) mole fractions of the admitted gas.
    #[must_use]
    pub const fn dry_mole_fractions(&self) -> [f64; 3] {
        self.dry
    }
}

/// ISO 9613-1 absorption [Np/m] of an admitted terrestrial-air `state`, at an
/// explicit relative humidity. The state's transport coefficients are **not**
/// used: ISO's classical term is the air fit, not Stokes–Kirchhoff.
///
/// # Errors
/// [`MaterialError::Parameters`] when `air` was not admitted for `state`, and
/// the ISO window refusals.
pub fn iso9613_absorption(
    air: &TerrestrialAir,
    state: &GasState,
    relative_humidity: f64,
    omega: f64,
) -> Result<f64, MaterialError> {
    air.check(state)
        .map_err(|refusal| MaterialError::Parameters {
            what: refusal.to_string(),
        })?;
    iso9613_absorption_neper_per_m(state.temperature, state.pressure, relative_humidity, omega)
}

/// ISO 9613-1 absorption [dB/km] — the unit of the published tables.
///
/// # Errors
/// Forwards [`iso9613_absorption_neper_per_m`] and refuses an unrepresentable
/// dB/km conversion with [`MaterialError::State`].
pub fn iso9613_absorption_db_per_km(
    temperature_k: f64,
    pressure_pa: f64,
    relative_humidity: f64,
    omega: f64,
) -> Result<f64, MaterialError> {
    let np = iso9613_absorption_neper_per_m(temperature_k, pressure_pa, relative_humidity, omega)?;
    finite_absorption(np * NP_TO_DB * 1000.0, "dB/km")
}

/// Free-field amplitude factor `e^{-α r} / r^n` for a compact source.
///
/// Plane: `n = 0`. Cylindrical: `n = 1/2`. Spherical: `n = 1`.
/// `range_m` must be positive. `alpha` is [Np/m].
#[must_use]
pub fn range_factor(range_m: f64, alpha_np_per_m: f64, spreading: Spreading) -> f64 {
    if !(range_m > 0.0
        && range_m.is_finite()
        && alpha_np_per_m >= 0.0
        && alpha_np_per_m.is_finite())
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
        let (air, witness) = TerrestrialAir::moist_air(293.15, 101_325.0, 0.0).expect("air");
        assert_eq!(
            air,
            GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("dry air")
        );
        let omega = 2.0 * core::f64::consts::PI * 1_000.0;
        let iso = iso9613_absorption(&witness, &air, 0.0, omega).expect("iso");
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

    #[test]
    fn finite_admitted_inputs_cannot_return_nonfinite_absorption() {
        for (pressure, humidity, omega) in [
            (101_325.0, 0.5, f64::MAX),
            (f64::MIN_POSITIVE, 0.5, 1_000.0),
        ] {
            let error = iso9613_absorption_neper_per_m(293.15, pressure, humidity, omega)
                .expect_err("unrepresentable absorption must refuse");
            assert!(matches!(error, MaterialError::State { .. }));
        }
        // Positive subnormal frequency may underflow to representable zero;
        // this guard must not reject finite nonnegative attenuation.
        let tiny = iso9613_absorption_neper_per_m(293.15, 101_325.0, 0.5, f64::from_bits(1))
            .expect("representable underflow");
        assert_eq!(tiny.to_bits(), 0.0_f64.to_bits());
    }

    fn medium(state: GasState, parts: &[(&str, f64)]) -> RelaxingMedium {
        RelaxingMedium::new(
            state,
            parts
                .iter()
                .map(|(name, x)| (SpeciesId::new(*name).expect("id"), *x))
                .collect(),
        )
        .expect("medium")
    }

    /// Terrestrial air is admitted; heliox, oxygen-poor air and CO2-rich air
    /// refuse, each naming its species. A witness does not transfer to a
    /// different gas state, and the legacy moist-air witness is bitwise the
    /// formula it always used.
    #[test]
    fn iso_refuses_non_standard_composition_naming_the_species() {
        let (humid, legacy) = TerrestrialAir::moist_air(293.15, 101_325.0, 0.5).expect("air");
        let omega = 2.0 * core::f64::consts::PI * 4_000.0;
        assert_eq!(
            iso9613_absorption(&legacy, &humid, 0.5, omega)
                .expect("legacy")
                .to_bits(),
            iso9613_absorption_neper_per_m(293.15, 101_325.0, 0.5, omega)
                .expect("formula")
                .to_bits()
        );
        let water = humid.water_mole_fraction;
        let dry = 1.0 - water;
        let air = medium(
            humid,
            &[
                ("N2", 0.7808 * dry),
                ("O2", 0.2095 * dry),
                ("Ar", 0.0093 * dry),
                ("CO2", 0.0004 * dry),
                ("H2O", water),
            ],
        );
        let admitted = TerrestrialAir::admit(&air).expect("terrestrial air");
        assert!(admitted.admits(&humid));
        let cases: [(&[(&str, f64)], &str); 4] = [
            (&[("He", 0.8), ("O2", 0.2)], "He"),
            (&[("N2", 0.85), ("O2", 0.15)], "O2"),
            (&[("N2", 0.7908), ("O2", 0.2092)], "N2"),
            (&[("N2", 0.786), ("O2", 0.209), ("CO2", 0.005)], "CO2"),
        ];
        for (parts, named) in cases {
            let refusal =
                TerrestrialAir::admit(&medium(humid, parts)).expect_err("not terrestrial air");
            assert_eq!(refusal.species(), Some(named), "{refusal}");
            assert!(refusal.to_string().contains(named), "{refusal}");
        }
        let (dry_air, _) = TerrestrialAir::moist_air(293.15, 101_325.0, 0.0).expect("dry");
        assert_eq!(
            admitted.check(&dry_air),
            Err(Iso9613Refusal::DifferentState)
        );
        assert!(matches!(
            iso9613_absorption(&admitted, &dry_air, 0.5, omega),
            Err(MaterialError::Parameters { .. })
        ));
    }

    /// ISO's own terms reassemble its absorption to rounding across its domain.
    #[test]
    fn iso_terms_reassemble_the_absorption() {
        for t in [253.15, 293.15, 323.15] {
            for rh in [0.0, 0.2, 0.8, 1.0] {
                let terms = iso9613_terms(t, 101_325.0, rh).expect("terms");
                for f in [50.0, 1_000.0, 8_000.0, 20_000.0] {
                    let f2 = f * f;
                    let relaxation: f64 = terms
                        .relaxation
                        .iter()
                        .map(|(strength, fr)| strength * f2 / (fr + f2 / fr))
                        .sum();
                    let rebuilt = terms.classical_np_per_m_hz2 * f2 + relaxation;
                    let direct = iso9613_absorption_neper_per_m(
                        t,
                        101_325.0,
                        rh,
                        2.0 * core::f64::consts::PI * f,
                    )
                    .expect("iso");
                    assert!(
                        (rebuilt / direct - 1.0).abs() <= 1e-13,
                        "T={t} RH={rh} f={f}"
                    );
                }
            }
        }
        assert!(iso9613_terms(100.0, 101_325.0, 0.5).is_err());
    }

    #[test]
    fn db_conversion_cannot_promote_finite_absorption_to_infinity() {
        let np = iso9613_absorption_neper_per_m(293.15, 1e-300, 0.0, 1e7)
            .expect("finite base absorption");
        assert!(np.is_finite() && np > 0.0);
        assert!(matches!(
            iso9613_absorption_db_per_km(293.15, 1e-300, 0.0, 1e7),
            Err(MaterialError::State { .. })
        ));
    }
}
