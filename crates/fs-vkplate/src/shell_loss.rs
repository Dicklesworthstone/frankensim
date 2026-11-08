//! A cymbal shell's loss from print (opt-in, MB 12552): B20 bronze's thermoelastic damping plus a
//! borrowed anelastic loss factor, in place of the two T60 points a preset declares.
//!
//! The kit's cymbal shells (crash, ride and the hi-hat's two) each carry two declared T60 points
//! (88cccbb78's estimates; the card ride's 10 kHz point from ruling 7486); the plate turns them
//! into its two-term loss (`s0 + s1 beta^2`, proportional to `a + b omega` on a tension-free plate).
//! [`printed_shell_t60`] replaces that pair, at the same two frequencies, with:
//! - thermoelastic damping from fs-material's insulated-slab law (every odd thermal mode, the
//!   fixed-strain heat capacity of Norris and Photiadis, `SlabThermoelasticBending::new_coupled`),
//!   as the cylindrical-bending loss factor `Im(Delta D_xx)/D_xx` of the shell's free edge, the
//!   thinnest region and the thickness the plate's two-point fit already uses;
//! - an anelastic loss factor measured in vacuum on cast bells (van Houten, Schoofs and van
//!   Campen 1997), borrowed.
//!
//! Above its relaxation peak the thermoelastic loss is nearly a constant decay rate and a
//! constant loss factor is a rate proportional to frequency, so the plate's two-term law carries
//! their sum between the two points (within 2 % from 500 Hz to 10 kHz on the crash and ride);
//! below the lower point it decays faster than the printed loss.
//!
//! What the pair does not hold, named in every receipt: radiation loss (none is carried for a
//! cymbal); the bells' own thermoelastic loss inside the borrowed figure (possibly counted twice);
//! the session temperature (declared). The thermoelastic part is an upper bound for its
//! conductivity: thicker regions lose less (as 1/h^2 above the peak), trace-free (twisting)
//! bending loses nothing thermoelastically, and the membrane strain of a shallow shell's lowest
//! modes loses far less than bending.

use crate::kit::{CymbalSpec, edge_thickness};
use fs_material::visco::slab_ted::SlabThermoelasticBending;
use std::f64::consts::PI;

/// Sommerfeld's Lorenz number, `pi^2/3 (k_B/e)^2` [W Ohm/K^2].
pub const SOMMERFELD_LORENZ: f64 = 2.443e-8;
/// 100 % IACS (International Annealed Copper Standard), `1/(1.7241 micro-ohm cm)` [S/m].
pub const IACS_S_PER_M: f64 = 5.8001e7;
/// C91300 (B20 bronze, 79-82 Cu, 18-20 Sn), CDA data reprinted by Concast: electrical
/// conductivity 7 % IACS at 68 F (293.15 K).
pub const C91300_IACS: (f64, f64) = (0.07, 293.15);
/// International Critical Tables Vol. V (1929) p. 120, Heat Capacity - Alloys, Cu-Sn: "Cu, 80%,
/// bell bronze; c_m(15-98 deg) = 0.36" J/g (Regnault 1843); Smithsonian Physical Tables
/// (Fowle) Table 262 agrees, 0.0862 cal/(g C) for 80 Cu + 20 Sn [J/(kg K)].
pub const BELL_BRONZE_SPECIFIC_HEAT: f64 = 360.0;
/// Hidnert, NBS J. Res. 12, RP665 (1934), Table 4: cast 84.84Cu-14.95Sn, 18.0e-6/C over
/// 20-100 C. (tin wt %, expansion 1/K)
pub const HIDNERT_1934_EXPANSION: (f64, f64) = (14.95, 18.0e-6);
/// Smithsonian Physical Tables (Fowle) Table 251: bronze 3Cu+1Sn, 25 % Sn (Daniell),
/// 0.1844e-4/C over 16.6-100 C. (tin wt %, expansion 1/K)
pub const DANIELL_EXPANSION: (f64, f64) = (25.0, 18.44e-6);
/// B20's nominal tin [wt %].
pub const B20_TIN: f64 = 20.0;
/// van Houten, Schoofs and van Campen, Damping of bells using experimental and numerical
/// methods, ICSV5 (1997): material loss factors of four Eijsbouts bells A-D in vacuum, each the
/// mean over its partials from the hum to 3 kHz; bronze composition and temperature not stated.
pub const VAN_HOUTEN_1997_BELLS: [f64; 4] = [1.4e-4, 1.3e-4, 2.6e-4, 1.6e-4];
/// TPRC Vol. 10, Thermal Diffusivity (Touloukian et al.), Table 89 p. 237: 87.56Cu-12.25Sn
/// (GZ-SnBz 12, Boehm and Wachtel 1969), 0.209 cm^2/s at 273 K and 0.234 at 373 K. A
/// different alloy: less tin conducts better, so for B20 an upper bound. ([K, m^2/s] pairs)
pub const TPRC_12_25_SN_DIFFUSIVITY: [(f64, f64); 2] = [(273.0, 0.209e-4), (373.0, 0.234e-4)];
/// The recording's room temperature, declared: no session temperature is printed [K].
pub const DECLARED_TEMPERATURE_K: f64 = 293.15;
/// The bells' own thermoelastic loss as a share of the smallest printed bell loss, at its hum
/// (220 Hz), with the upper-bound diffusivity, for any wall of 10 mm or more (the paper prints no
/// wall thickness): the most of the borrowed figure that can be counted twice.
pub const VAN_HOUTEN_DOUBLE_COUNT_SHARE: f64 = 0.094;
/// `ln 1000`: a T60 is the time the decay rate takes to lose 60 dB.
const LN_1000: f64 = 6.907_755_278_982_137;

/// How a row states heat conduction through the shell's thickness.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Conduction {
    /// Through-thickness conductivity [W/(m K)].
    Conductivity(f64),
    /// Thermal diffusivity [m^2/s], turned into conductivity with the shell's own `rho c_p`.
    Diffusivity(f64),
}

/// One printed-loss row: the thermal data of the thermoelastic part and the anelastic factor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrintedShellLoss {
    /// The `--shell-loss` value that selects it.
    pub label: &'static str,
    /// Heat conduction through the thickness.
    pub conduction: Conduction,
    /// Linear thermal expansion [1/K].
    pub expansion_per_k: f64,
    /// Specific heat at constant stress [J/(kg K)].
    pub specific_heat_j_kg_k: f64,
    /// Absolute temperature [K].
    pub temperature_k: f64,
    /// Anelastic (non-thermoelastic) loss factor, constant in frequency.
    pub anelastic_loss_factor: f64,
    /// The sources and what the row leaves out; joins the cymbal's notes.
    pub receipt: &'static str,
}

/// Wiedemann-Franz with Sommerfeld's Lorenz number on a conductivity printed in % IACS: the
/// electronic part of the thermal conductivity [W/(m K)]. The lattice adds to it, so for a
/// bronze this is a lower estimate.
#[must_use]
pub fn wiedemann_franz_conductivity(iacs_fraction: f64, temperature_k: f64) -> f64 {
    SOMMERFELD_LORENZ * iacs_fraction * IACS_S_PER_M * temperature_k
}

fn linear(a: (f64, f64), b: (f64, f64), x: f64) -> f64 {
    a.1 + (b.1 - a.1) * (x - a.0) / (b.0 - a.0)
}

/// B20's expansion, linear in tin between Hidnert's 14.95 % and Daniell's 25 % Sn [1/K].
#[must_use]
pub fn b20_expansion() -> f64 {
    linear(HIDNERT_1934_EXPANSION, DANIELL_EXPANSION, B20_TIN)
}

/// TPRC's 12.25 % Sn diffusivity, linear between its printed 273 K and 373 K points [m^2/s].
#[must_use]
pub fn tprc_upper_diffusivity(temperature_k: f64) -> f64 {
    linear(TPRC_12_25_SN_DIFFUSIVITY[0], TPRC_12_25_SN_DIFFUSIVITY[1], temperature_k)
}

/// The mean of van Houten's four bells.
#[must_use]
pub fn van_houten_mean() -> f64 {
    VAN_HOUTEN_1997_BELLS.iter().sum::<f64>() / VAN_HOUTEN_1997_BELLS.len() as f64
}

/// The receipt of the row that may post (`--shell-loss b20`).
pub const B20_RECEIPT: &str = "cymbal shell loss from print, b20 (opt-in, MB 12552): the declared T60 pair \
replaced at its two frequencies by B20 bronze's thermoelastic damping (fs-material's insulated-slab law, \
cylindrical bending at the free edge; c_p 0.36 J/g, ICT Vol. V 1929 p. 120, 80Cu-20Sn bell bronze, Regnault; \
expansion 18.22e-6/K, linear in tin between Hidnert 1934 at 14.95 % Sn and Daniell's 25 % Sn in Smithsonian \
Table 251; conductivity 29.1 W/(m K) by Wiedemann-Franz with Sommerfeld's Lorenz number on C91300's printed 7 % \
IACS, CDA via Concast, the electronic part only, so a lower estimate) plus the anelastic loss factor 1.725e-4, the \
mean of four cast bells measured in vacuum (van Houten, Schoofs and van Campen 1997, band 1.3e-4 to 2.6e-4), \
borrowed from cast bells of unstated composition and carried past its printed 3 kHz; it holds those bells' own \
thermoelastic loss, possibly counted twice up to 9.4 % of it. Temperature 20 C, declared. Radiation loss missing. \
The thermoelastic part is an upper bound for its conductivity: thinnest region, cylindrical bending for every \
shape (twisting loses nothing thermoelastically), membrane strain treated as bending";
/// The receipt of the low band end (`--shell-loss b20-low`, never posts).
pub const B20_LOW_RECEIPT: &str = "cymbal shell loss from print, b20-low band end (never posts): as b20 \
with expansion 18.0e-6/K (Hidnert 1934, 14.95 % Sn) and the smallest bell's anelastic loss factor 1.3e-4 \
(van Houten 1997, bell B); conductivity 29.1 W/(m K), Wiedemann-Franz on C91300's 7 % IACS. Radiation loss missing";
/// The receipt of the high band end (`--shell-loss b20-high`, never posts).
pub const B20_HIGH_RECEIPT: &str = "cymbal shell loss from print, b20-high band end (never posts): as b20 \
with the thermal diffusivity of a DIFFERENT alloy, 87.56Cu-12.25Sn (TPRC Vol. 10 Table 89, Boehm and Wachtel \
1969, linear between 273 K and 373 K), borrowed as an upper bound on B20's (less tin conducts better), expansion \
18.44e-6/K (Daniell's 25 % Sn) and the largest bell's anelastic loss factor 2.6e-4 (van Houten 1997, bell C). \
Radiation loss missing";

/// The row selected by a `--shell-loss` value: `b20` (posts by the count rule), `b20-low` or
/// `b20-high` (band ends, never post).
#[must_use]
pub fn b20(label: &str) -> Option<PrintedShellLoss> {
    let (iacs, at) = C91300_IACS;
    let wf = Conduction::Conductivity(wiedemann_franz_conductivity(iacs, at));
    let (lo, hi) = VAN_HOUTEN_1997_BELLS.iter().fold((f64::INFINITY, 0.0f64), |(l, h), &x| (l.min(x), h.max(x)));
    let row = |label, conduction, expansion_per_k, anelastic_loss_factor, receipt| PrintedShellLoss {
        label,
        conduction,
        expansion_per_k,
        specific_heat_j_kg_k: BELL_BRONZE_SPECIFIC_HEAT,
        temperature_k: DECLARED_TEMPERATURE_K,
        anelastic_loss_factor,
        receipt,
    };
    match label {
        "b20" => Some(row("b20", wf, b20_expansion(), van_houten_mean(), B20_RECEIPT)),
        "b20-low" => Some(row("b20-low", wf, HIDNERT_1934_EXPANSION.1, lo, B20_LOW_RECEIPT)),
        "b20-high" => Some(row("b20-high", Conduction::Diffusivity(tprc_upper_diffusivity(DECLARED_TEMPERATURE_K)),
            DANIELL_EXPANSION.1, hi, B20_HIGH_RECEIPT)),
        _ => None,
    }
}

fn check(loss: &PrintedShellLoss, elastic: (f64, f64, f64)) -> Result<(), String> {
    let (e, nu, rho) = elastic;
    let conduction = match loss.conduction {
        Conduction::Conductivity(x) | Conduction::Diffusivity(x) => x,
    };
    let positive = |x: f64| x.is_finite() && x > 0.0;
    if !(positive(e) && positive(rho) && nu.is_finite() && nu > -1.0 && nu < 0.5) {
        return Err(format!("shell loss {}: elastic values (E {e} Pa, nu {nu}, rho {rho} kg/m^3) are not an isotropic solid", loss.label));
    }
    if !(positive(conduction) && positive(loss.expansion_per_k) && positive(loss.specific_heat_j_kg_k)
        && positive(loss.temperature_k) && positive(loss.anelastic_loss_factor) && loss.anelastic_loss_factor < 1.0)
    {
        return Err(format!("shell loss {}: conduction, expansion, specific heat and temperature must be positive and finite, and the anelastic loss factor in (0, 1)", loss.label));
    }
    Ok(())
}

fn slab_loss_factor(loss: &PrintedShellLoss, elastic: (f64, f64, f64), thickness_m: f64, frequency_hz: f64, coupled: bool) -> Result<f64, String> {
    check(loss, elastic)?;
    if !(thickness_m.is_finite() && thickness_m > 0.0 && frequency_hz.is_finite() && frequency_hz > 0.0) {
        return Err(format!("shell loss {}: thickness {thickness_m} m and frequency {frequency_hz} Hz must be positive and finite", loss.label));
    }
    let (e, nu, rho) = elastic;
    let c = e / (1.0 - nu * nu);
    let q = [c, c * nu, 0.0, c * nu, c, 0.0, 0.0, 0.0, c * (1.0 - nu) / 2.0];
    let alpha = loss.expansion_per_k;
    let capacity = rho * loss.specific_heat_j_kg_k;
    let k = match loss.conduction {
        Conduction::Conductivity(k) => k,
        Conduction::Diffusivity(d) => d * capacity,
    };
    let fail = |err: fs_material::visco::ViscoError| format!("shell loss {}: {err}", loss.label);
    let slab = if coupled {
        SlabThermoelasticBending::new_coupled(q, [alpha, alpha, 0.0], loss.temperature_k, capacity, k, thickness_m)
    } else {
        SlabThermoelasticBending::new(q, [alpha, alpha, 0.0], loss.temperature_k, capacity, k, thickness_m)
    }
    .map_err(fail)?;
    let imaginary = slab.loss_matrix(2.0 * PI * frequency_hz).map_err(fail)?[0];
    Ok(imaginary / (thickness_m.powi(3) / 12.0 * c))
}

/// Thermoelastic loss factor of cylindrical bending across a slab of `thickness_m` at
/// `frequency_hz`: `Im(Delta D_xx)/D_xx` from fs-material's coupled insulated-slab law, for the
/// isotropic plane-stress solid `elastic` = (E Pa, nu, rho kg/m^3) with the row's thermal data.
///
/// # Errors
/// Invalid elastic, thermal or geometric data, or a slab the law refuses.
pub fn thermoelastic_loss_factor(loss: &PrintedShellLoss, elastic: (f64, f64, f64), thickness_m: f64, frequency_hz: f64) -> Result<f64, String> {
    slab_loss_factor(loss, elastic, thickness_m, frequency_hz, true)
}

/// The T60 pair the row gives the shell `c`, at the frequencies of its declared pair: the
/// thermoelastic loss at the free edge plus the anelastic factor, `T60 = ln 1000/(eta pi f)`.
///
/// # Errors
/// As [`thermoelastic_loss_factor`].
pub fn printed_shell_t60(c: &CymbalSpec, elastic: (f64, f64, f64), loss: &PrintedShellLoss) -> Result<[(f64, f64); 2], String> {
    let h = edge_thickness(c);
    let mut pair = [(0.0, 0.0); 2];
    for (slot, &(f, _)) in pair.iter_mut().zip(c.t60.iter()) {
        let eta = thermoelastic_loss_factor(loss, elastic, h, f)? + loss.anelastic_loss_factor;
        *slot = (f, LN_1000 / (eta * PI * f));
    }
    Ok(pair)
}

/// The shell `c` with its T60 pair from print and the row's receipt as its provenance, and one
/// receipt line with the numbers (what each point holds and the pair it replaces).
///
/// # Errors
/// As [`printed_shell_t60`].
pub fn with_printed_shell_loss(c: &CymbalSpec, elastic: (f64, f64, f64), loss: &PrintedShellLoss) -> Result<(CymbalSpec, String), String> {
    let t60 = printed_shell_t60(c, elastic, loss)?;
    let h = edge_thickness(c);
    let mut parts = Vec::with_capacity(2);
    for &(f, t) in &t60 {
        let ted = thermoelastic_loss_factor(loss, elastic, h, f)?;
        parts.push(format!("{t:.3} s at {f:.0} Hz (thermoelastic {ted:.3e} + anelastic {:.3e})", loss.anelastic_loss_factor));
    }
    let line = format!(
        "{}: shell loss {} at the free edge, {:.3} mm: T60 {}; replaces the declared {} s at {:.0} Hz and {} s at {:.0} Hz; shell dimensions: {}",
        c.name, loss.label, h * 1e3, parts.join(", "), c.t60[0].1, c.t60[0].0, c.t60[1].1, c.t60[1].0,
        c.provenance.unwrap_or("the preset's"),
    );
    Ok((CymbalSpec { t60, provenance: Some(loss.receipt), ..*c }, line))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kit::{BRONZE, card_ride, crash, hihat_spec};
    use fs_material::visco::slab_ted::slab_relaxation_shape;

    fn row(label: &str) -> PrintedShellLoss {
        b20(label).expect("a b20 row")
    }

    #[test]
    fn printed_inputs_reproduce_their_sources() {
        let k = wiedemann_franz_conductivity(C91300_IACS.0, C91300_IACS.1);
        assert!((k - 29.08).abs() < 0.01, "Wiedemann-Franz conductivity {k}");
        assert!((b20_expansion() - 18.221e-6).abs() < 1e-9, "expansion {}", b20_expansion());
        assert!((van_houten_mean() - 1.725e-4).abs() < 1e-12);
        let d = tprc_upper_diffusivity(DECLARED_TEMPERATURE_K);
        assert!((d - 0.21404e-4).abs() < 1e-8, "diffusivity {d}");
        assert!(b20("b21").is_none());
    }

    #[test]
    fn weak_coupling_is_the_closed_form_plate_strength_times_the_slab_shape() {
        // Delta_p = E alpha^2 T (1 + nu)/(rho c_p (1 - nu)) for cylindrical bending, tau_0 = h^2 rho c_p/(pi^2 k)
        let loss = row("b20");
        let (e, nu, rho) = BRONZE;
        let Conduction::Conductivity(k) = loss.conduction else { panic!("b20 states a conductivity") };
        let capacity = rho * loss.specific_heat_j_kg_k;
        let strength = e * loss.expansion_per_k.powi(2) * loss.temperature_k * (1.0 + nu) / (capacity * (1.0 - nu));
        for (h, f) in [(0.95e-3, 500.0), (0.95e-3, 10_000.0), (1.6e-3, 500.0), (1.6e-3, 10_000.0), (0.95e-3, 20.0)] {
            let tau = h * h * capacity / (PI * PI * k);
            let closed = strength * slab_relaxation_shape(2.0 * PI * f, tau);
            let weak = slab_loss_factor(&loss, BRONZE, h, f, false).unwrap();
            assert!((weak / closed - 1.0).abs() < 1e-12, "h {h} f {f}: {weak} vs {closed}");
            let coupled = thermoelastic_loss_factor(&loss, BRONZE, h, f).unwrap();
            assert!((coupled / weak - 1.0).abs() < 0.03, "h {h} f {f}: coupled {coupled} vs weak {weak}");
        }
    }

    #[test]
    fn above_the_peak_the_loss_falls_as_thickness_squared() {
        let loss = row("b20");
        let thin = thermoelastic_loss_factor(&loss, BRONZE, 0.8e-3, 10_000.0).unwrap();
        let thick = thermoelastic_loss_factor(&loss, BRONZE, 1.6e-3, 10_000.0).unwrap();
        assert!((thin / thick / 4.0 - 1.0).abs() < 0.05, "ratio {}", thin / thick);
    }

    #[test]
    fn the_pairs_are_the_plan_lines_and_the_band_holds_the_row() {
        // plan lines 12547 and 12922 (weak coupling): crash (10.68 s, 1.19 s), ride (16.92 s, 1.24 s), hi-hat top
        // (9.39 s, 1.17 s) and bottom (12.45 s, 1.21 s), each at its own free edge; the coupled law moves them a little
        let (hat_top, hat_bottom) = hihat_spec();
        let shells = [(crash(), [10.68, 1.19]), (card_ride(), [16.92, 1.24]), (hat_top, [9.39, 1.17]), (hat_bottom, [12.45, 1.21])];
        for (spec, plan) in shells {
            let pair = printed_shell_t60(&spec, BRONZE, &row("b20")).unwrap();
            let low = printed_shell_t60(&spec, BRONZE, &row("b20-low")).unwrap();
            let high = printed_shell_t60(&spec, BRONZE, &row("b20-high")).unwrap();
            for i in 0..2 {
                assert_eq!(pair[i].0, spec.t60[i].0);
                assert!((pair[i].1 / plan[i] - 1.0).abs() < 0.03, "{}: {:?} vs plan {plan:?}", spec.name, pair);
                assert!(high[i].1 < pair[i].1 && pair[i].1 < low[i].1, "{}: band {:?} {:?} {:?}", spec.name, high, pair, low);
            }
        }
    }

    #[test]
    fn the_spec_keeps_its_shape_and_names_its_receipt() {
        let c = crash();
        let (printed, line) = with_printed_shell_loss(&c, BRONZE, &row("b20")).unwrap();
        assert_eq!(printed.bow_thickness, c.bow_thickness);
        assert_eq!(printed.radius, c.radius);
        assert_eq!(printed.grid, c.grid);
        assert_eq!(printed.provenance, Some(B20_RECEIPT));
        assert_ne!(printed.t60, c.t60);
        assert!(line.contains("replaces the declared 10 s at 500 Hz and 2 s at 10000 Hz"), "{line}");
        assert!(B20_RECEIPT.contains("Radiation loss missing") && B20_RECEIPT.contains("possibly counted twice"));
    }

    #[test]
    fn invalid_data_refuse() {
        let mut loss = row("b20");
        loss.anelastic_loss_factor = -1e-4;
        assert!(printed_shell_t60(&crash(), BRONZE, &loss).is_err());
        let loss = row("b20");
        assert!(thermoelastic_loss_factor(&loss, BRONZE, 0.0, 500.0).is_err());
        assert!(thermoelastic_loss_factor(&loss, BRONZE, 1e-3, f64::NAN).is_err());
        assert!(thermoelastic_loss_factor(&loss, (1.05e11, 0.6, 8600.0), 1e-3, 500.0).is_err());
        let mut loss = row("b20");
        loss.conduction = Conduction::Diffusivity(0.0);
        assert!(printed_shell_t60(&crash(), BRONZE, &loss).is_err());
    }
}
