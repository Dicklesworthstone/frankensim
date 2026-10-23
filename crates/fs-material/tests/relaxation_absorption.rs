//! Generic molecular-relaxation law (IRP M0.1-W2, bead pexfp.18): G0 algebra,
//! causality and passivity, refusals, and the ISO 9613-1 rate variant.
//!
//! Mathematical fixtures use synthetic species ids (`A`..`D`) and declared
//! times: they test the law, not gas data. The ISO comparison uses the
//! standard's own characteristic temperatures 2239.1 K and 3352 K (Bass et al.
//! 1995, JASA 97:680, Eq. 3; a U.S. Government work) and ISO's dry-air O2/N2
//! fractions 0.209/0.781.

use fs_material::MaterialError;
use fs_material::gas::{ConductivityModel, GasSpec, GasState};
use fs_material::iso9613::{iso9613_absorption_neper_per_m, iso9613_relaxation_frequencies_hz};
use fs_material::relaxation::{
    InternalReservoir, IsoRelaxingSpecies, MAX_RELAXATION_PROCESSES, PartnerRate,
    RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium, SpeciesRelaxation,
    VibrationalMode,
};
use fs_thermochem::SpeciesId;

const TAU: f64 = core::f64::consts::TAU;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).expect("species id")
}

/// A calorically perfect test gas whose c_v = R_s / (gamma - 1) has room for
/// the declared reservoirs.
fn test_gas(gamma: f64, temperature: f64, pressure: f64) -> GasState {
    let spec = GasSpec {
        molar_mass: 28.9644e-3,
        gamma,
        sutherland_beta: 1.458e-6,
        sutherland_s: 110.4,
        conductivity: ConductivityModel::Eucken,
    };
    GasState::try_new(&spec, temperature, pressure).expect("test gas")
}

fn timed(tau_s: f64, pressure: f64) -> RelaxationRateLaw {
    RelaxationRateLaw::PartnerIndependent {
        p_tau_ref_pa_s: tau_s * pressure,
        t_ref_k: 293.15,
        exponent: 0.0,
        t_min_k: 100.0,
        t_max_k: 2000.0,
    }
}

fn vib(theta_k: f64, degeneracy: u8, rate: RelaxationRateLaw) -> RelaxationProcess {
    RelaxationProcess {
        reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
            theta_k,
            degeneracy,
        }]),
        rate,
    }
}

fn rot(degrees_of_freedom: u8, rate: RelaxationRateLaw) -> RelaxationProcess {
    RelaxationProcess {
        reservoir: InternalReservoir::Rotational { degrees_of_freedom },
        rate,
    }
}

fn species(name: &str, processes: Vec<RelaxationProcess>) -> SpeciesRelaxation {
    SpeciesRelaxation {
        species: id(name),
        processes,
    }
}

fn medium(gas: GasState, parts: &[(&str, f64)]) -> RelaxingMedium {
    RelaxingMedium::new(gas, parts.iter().map(|(name, x)| (id(name), *x)).collect())
        .expect("medium")
}

fn refusal_text(error: &MaterialError) -> String {
    error.to_string()
}

#[test]
fn g0_classical_limit_is_stokes_kirchhoff_bitwise() {
    let gas = test_gas(1.4, 293.15, 101_325.0);
    let law = RelaxationAbsorption::new(vec![species("A", vec![])]).expect("law");
    let state = law.evaluate(&medium(gas, &[("A", 1.0)])).expect("state");
    assert!(state.processes().is_empty());
    for f in [20.0, 440.0, 1_000.0, 8_000.0, 20_000.0, 1.0e6] {
        let omega = TAU * f;
        let alpha = state.absorption_np_per_m(omega).expect("alpha");
        let stokes = gas.stokes_kirchhoff_absorption(omega);
        assert_eq!(
            alpha.to_bits(),
            stokes.to_bits(),
            "f={f}: {alpha} vs {stokes}"
        );
    }
    println!(
        "{{\"bead\":\"M0.1-W2\",\"test\":\"g0_classical_limit\",\"verdict\":\"pass\",\"points\":6}}"
    );
}

/// Exact single relaxation in an otherwise frozen gas. With `e^{-iwt}`,
/// `c_v(w) = (c_v0 - i w tau c_vinf)/(1 - i w tau)` and likewise for `c_p`, so
/// `k^2 = (w/c0)^2 (1 - i w tau_v)/(1 - i w tau_p)` with
/// `tau_v = tau c_vinf/c_v0` and `tau_p = tau c_pinf/c_p0`. The phase of `k^2`,
/// `phi = atan(w tau_p) - atan(w tau_v)`, peaks at `w* = 1/sqrt(tau_v tau_p)`,
/// where `tan(phi*) = (r - 1/r)/2`, `r = sqrt(tau_p/tau_v)`. Since
/// `alpha lambda = 2 pi tan(phi/2)`, the peak absorption per wavelength is
/// `2 pi tan(phi*/2)`. Small-strength check (Pierce, Acoustics, Sec. 10-7):
/// `(alpha lambda)_max -> (pi/2) R_s C / (c_v c_p)`.
#[test]
fn g0_single_process_closed_form() {
    let p = 101_325.0;
    let t = 293.15;
    let gas = test_gas(1.4, t, p);
    let tau = 1.0e-5;
    let law = RelaxationAbsorption::new(vec![species("A", vec![vib(960.0, 2, timed(tau, p))])])
        .expect("law");
    let state = law.evaluate(&medium(gas, &[("A", 1.0)])).expect("state");
    let r_s = gas.specific_gas_constant;
    let cv0 = gas.specific_heat_cp - r_s;
    let c = state.processes()[0].heat_capacity_j_kg_k;
    let cvi = state.frozen_cv_j_kg_k();
    assert!(
        ((cv0 - cvi) - c).abs() <= 1e-12 * cv0,
        "frozen c_v bookkeeping"
    );
    let (cp0, cpi) = (cv0 + r_s, cvi + r_s);
    let tau_v = tau * cvi / cv0;
    let tau_p = tau * cpi / cp0;
    let omega_peak = 1.0 / (tau_v * tau_p).sqrt();
    let r = (tau_p / tau_v).sqrt();
    let phi = ((r - 1.0 / r) / 2.0).atan();
    let closed = TAU * (phi / 2.0).tan();
    let alpha_lambda = |omega: f64| {
        let k = state.relaxing_wavenumber(omega).expect("k");
        TAU * k.im / k.re
    };
    let at_peak = alpha_lambda(omega_peak);
    let rel = (at_peak - closed).abs() / closed;
    assert!(
        rel <= 1e-12,
        "peak alpha*lambda {at_peak} vs closed form {closed} (rel {rel:e})"
    );
    for side in [1.0 - 1e-3, 1.0 + 1e-3] {
        assert!(
            alpha_lambda(omega_peak * side) < at_peak,
            "peak is not a maximum at side {side}"
        );
    }
    let small_strength = core::f64::consts::FRAC_PI_2 * r_s * c / (cv0 * cp0);
    println!(
        "{{\"bead\":\"M0.1-W2\",\"test\":\"g0_single_process_closed_form\",\"inputs\":{{\"theta_K\":960,\"g\":2,\"tau_s\":{tau},\"T_K\":{t}}},\"observed\":{{\"alpha_lambda_max\":{at_peak:e},\"closed_form\":{closed:e},\"small_strength_pierce\":{small_strength:e},\"f_peak_Hz\":{}}},\"relative_error\":{rel:e},\"tolerance\":1e-12,\"verdict\":\"pass\"}}",
        omega_peak / TAU
    );
    // Equilibrium and frozen limits of the phase speed. The dispersion departs
    // from either limit as (w tau)^2 or (w tau)^-2 times (c_inf - c0)/c0 (about
    // 0.06 here), so the probes sit six and seven decades from the relaxation.
    let c0 = (r_s * t * cp0 / cv0).sqrt();
    let low = state.phase_speed_m_s(1e-6 / tau).expect("low");
    let high = state.phase_speed_m_s(1e7 / tau).expect("high");
    assert!(
        (low - c0).abs() / c0 < 1e-9,
        "low-frequency speed {low} vs equilibrium {c0}"
    );
    let ci = state.frozen_sound_speed_m_s();
    assert!(
        (high - ci).abs() / ci < 1e-9,
        "high-frequency speed {high} vs frozen {ci}"
    );
    // Sign: decay grows with path length (Im k > 0 under e^{i(kx - wt)}).
    let k = state.complex_wavenumber(omega_peak).expect("k");
    assert!(k.im > 0.0 && (-k.im * 2.0).exp() < (-k.im * 1.0).exp());
}

#[test]
fn g0_passive_over_random_admitted_states() {
    let mut rng = fs_rand::StreamKey {
        seed: 0x0018_5eed,
        kernel: 18,
        tile: 0,
    }
    .stream();
    let names = ["A", "B", "C", "D"];
    let mut checked = 0_u64;
    for _ in 0..10_000 {
        let t = 200.0 + 200.0 * rng.next_f64();
        let p = 1.0e3 * (200.0_f64).powf(rng.next_f64());
        let gas = test_gas(1.2, t, p);
        let n = 1 + rng.next_below(4) as usize;
        let mut data = Vec::new();
        let mut weights = Vec::new();
        for name in &names[..n] {
            let mut processes = Vec::new();
            if rng.next_below(2) == 1 {
                let dof = 2 + rng.next_below(2) as u8;
                let p_tau = 10.0_f64.powf(-8.0 + 10.0 * rng.next_f64());
                processes.push(rot(
                    dof,
                    RelaxationRateLaw::PartnerIndependent {
                        p_tau_ref_pa_s: p_tau,
                        t_ref_k: 300.0,
                        exponent: -2.0 + 4.0 * rng.next_f64(),
                        t_min_k: 100.0,
                        t_max_k: 500.0,
                    },
                ));
            }
            if rng.next_below(2) == 1 {
                let modes = (0..1 + rng.next_below(2))
                    .map(|_| VibrationalMode {
                        theta_k: 200.0 + 5_800.0 * rng.next_f64(),
                        degeneracy: 1,
                    })
                    .collect();
                let p_tau = 10.0_f64.powf(-8.0 + 10.0 * rng.next_f64());
                processes.push(RelaxationProcess {
                    reservoir: InternalReservoir::Vibrational(modes),
                    rate: RelaxationRateLaw::PartnerIndependent {
                        p_tau_ref_pa_s: p_tau,
                        t_ref_k: 300.0,
                        exponent: -2.0 + 4.0 * rng.next_f64(),
                        t_min_k: 100.0,
                        t_max_k: 500.0,
                    },
                });
            }
            data.push(species(name, processes));
            weights.push(0.05 + rng.next_f64());
        }
        let total: f64 = weights.iter().sum();
        let parts: Vec<(SpeciesId, f64)> = names[..n]
            .iter()
            .zip(&weights)
            .map(|(name, w)| (id(name), w / total))
            .collect();
        let law = RelaxationAbsorption::new(data).expect("law");
        let state = law
            .evaluate(&RelaxingMedium::new(gas, parts).expect("medium"))
            .expect("state");
        let r_s = gas.specific_gas_constant;
        let cv0 = gas.specific_heat_cp - r_s;
        let c0 = (r_s * t * (cv0 + r_s) / cv0).sqrt();
        let ci = state.frozen_sound_speed_m_s();
        for f in [10.0, 1.0e3, 1.0e5, 1.0e7] {
            let omega = TAU * f;
            let k = state.relaxing_wavenumber(omega).expect("k");
            let alpha = state.absorption_np_per_m(omega).expect("alpha");
            let speed = omega / k.re;
            assert!(
                k.im >= 0.0 && alpha >= 0.0,
                "active medium: Im k {} alpha {alpha}",
                k.im
            );
            assert!(
                speed >= c0 * (1.0 - 1e-12) && speed <= ci * (1.0 + 1e-12),
                "phase speed {speed} outside [{c0}, {ci}]"
            );
            checked += 1;
        }
    }
    println!(
        "{{\"bead\":\"M0.1-W2\",\"test\":\"g0_passive_over_random_admitted_states\",\"states\":10000,\"points\":{checked},\"verdict\":\"pass\"}}"
    );
}

// Eight-point Gauss-Legendre nodes and weights on [-1, 1].
const GL_X: [f64; 4] = [
    0.183_434_642_495_649_8,
    0.525_532_409_916_329,
    0.796_666_477_413_626_7,
    0.960_289_856_497_536_3,
];
const GL_W: [f64; 4] = [
    0.362_683_783_378_362,
    0.313_706_645_877_887_3,
    0.222_381_034_453_374_5,
    0.101_228_536_290_376_3,
];

fn gauss_legendre(a: f64, b: f64, panels: usize, f: &mut impl FnMut(f64) -> f64) -> f64 {
    let h = (b - a) / panels as f64;
    let mut sum = 0.0;
    for i in 0..panels {
        let mid = a + (i as f64 + 0.5) * h;
        for (x, w) in GL_X.iter().zip(GL_W) {
            sum += w * (f(mid - 0.5 * h * x) + f(mid + 0.5 * h * x));
        }
    }
    0.5 * h * sum
}

/// `1/c(w) - 1/c_inf = (2/pi) P int_0^inf alpha(w')/(w'^2 - w^2) dw'` for the
/// relaxing reservoirs (the Stokes-Kirchhoff term is a first-order attenuation
/// whose dispersion partner is outside this model and is excluded). The
/// singular part is subtracted (`P int_0^inf dw'/(w'^2 - w^2) = 0`) and the
/// regular integrand is integrated in `ln w'` with Gauss-Legendre panels on
/// [w_min, w_max], split at `w`. The tails are added analytically for the
/// `alpha(w)` part; the residual tail terms are bounded and reported.
#[test]
fn g0_kramers_kronig_consistency() {
    let p = 101_325.0;
    let gas = test_gas(1.25, 293.15, p);
    let law = RelaxationAbsorption::new(vec![
        species(
            "A",
            vec![rot(2, timed(2.0e-10, p)), vib(3352.0, 1, timed(5.0e-4, p))],
        ),
        species(
            "B",
            vec![
                vib(2239.1, 1, timed(5.0e-6, p)),
                vib(960.0, 2, timed(2.0e-7, p)),
            ],
        ),
    ])
    .expect("law");
    let state = law
        .evaluate(&medium(gas, &[("A", 0.7), ("B", 0.3)]))
        .expect("state");
    let alpha = |omega: f64| state.relaxing_wavenumber(omega).expect("k").im;
    let (w_min, w_max): (f64, f64) = (1.0e-4 / 5.0e-4, 1.0e6 / 2.0e-10);
    let inv_ci = 1.0 / state.frozen_sound_speed_m_s();
    let mut worst = 0.0_f64;
    for i in 0..13 {
        let f = 20.0 * 1000.0_f64.powf(i as f64 / 12.0);
        let w = TAU * f;
        let a_w = alpha(w);
        let mut integrand = |u: f64| {
            let wp = u.exp();
            wp * (alpha(wp) - a_w) / (wp * wp - w * w)
        };
        let (u_min, u0, u_max) = (w_min.ln(), w.ln(), w_max.ln());
        let lower = gauss_legendre(
            u_min,
            u0,
            ((u0 - u_min) / 0.05).ceil() as usize,
            &mut integrand,
        );
        let upper = gauss_legendre(
            u0,
            u_max,
            ((u_max - u0) / 0.05).ceil() as usize,
            &mut integrand,
        );
        let log_ratio = |big: f64| ((big + w) / (big - w)).ln() / (2.0 * w);
        // Analytic -alpha(w) parts of both tails, plus the alpha(w_max) plateau estimate above.
        let tail_low = a_w * ((w + w_min) / (w - w_min)).ln() / (2.0 * w);
        let tail_high = (alpha(w_max) - a_w) * log_ratio(w_max);
        let bound_low = alpha(w_min) * w_min / (w * w - w_min * w_min);
        let bound_high = (alpha(w_max) - alpha(w_max / 10.0)).abs() * log_ratio(w_max);
        let rhs = 2.0 / core::f64::consts::PI * (lower + upper + tail_low + tail_high);
        let lhs = state.relaxing_wavenumber(w).expect("k").re / w - inv_ci;
        let rel = (rhs - lhs).abs() / lhs.abs();
        let bound_rel = 2.0 / core::f64::consts::PI * (bound_low + bound_high) / lhs.abs();
        worst = worst.max(rel);
        println!(
            "{{\"bead\":\"M0.1-W2\",\"test\":\"g0_kramers_kronig_consistency\",\"f_Hz\":{f:.3},\"lhs\":{lhs:e},\"rhs\":{rhs:e},\"relative_error\":{rel:e},\"tail_bound_rel\":{bound_rel:e},\"tolerance\":1e-6}}"
        );
        assert!(bound_rel < 1e-8, "tail bound {bound_rel:e} at {f} Hz");
        assert!(
            rel <= 1e-6,
            "KK mismatch {rel:e} at {f} Hz: lhs {lhs:e} rhs {rhs:e}"
        );
    }
    println!(
        "{{\"bead\":\"M0.1-W2\",\"test\":\"g0_kramers_kronig_consistency\",\"worst\":{worst:e},\"verdict\":\"pass\"}}"
    );
}

fn iso_air_law() -> RelaxationAbsorption {
    RelaxationAbsorption::new(vec![
        species(
            "O2",
            vec![vib(
                2239.1,
                1,
                RelaxationRateLaw::Iso9613TerrestrialAir(IsoRelaxingSpecies::Oxygen),
            )],
        ),
        species(
            "N2",
            vec![vib(
                3352.0,
                1,
                RelaxationRateLaw::Iso9613TerrestrialAir(IsoRelaxingSpecies::Nitrogen),
            )],
        ),
        species("Ar", vec![]),
        // Declared non-relaxing for this vibrational comparison only: ISO's
        // formula carries no water or CO2 reservoir of its own.
        species("H2O", vec![]),
        species("CO2", vec![]),
    ])
    .expect("ISO air law")
}

#[test]
fn iso_rate_variant_is_the_standards_relaxation_frequency() {
    let law = iso_air_law();
    for (t, p, x_w) in [
        (293.15, 101_325.0, 0.0),
        (283.15, 90_000.0, 0.008),
        (313.15, 101_325.0, 0.03),
    ] {
        let gas = test_gas(1.4, t, p);
        let dry = 1.0 - x_w;
        let state = law
            .evaluate(&medium(
                gas,
                &[
                    ("N2", 0.781 * dry),
                    ("O2", 0.209 * dry),
                    ("Ar", 0.010 * dry),
                    ("H2O", x_w),
                ],
            ))
            .expect("state");
        let (fr_o, fr_n) = iso9613_relaxation_frequencies_hz(t, p, 100.0 * x_w).expect("iso fr");
        let tau_of = |name: &str| {
            state
                .processes()
                .iter()
                .find(|e| e.species.as_str() == name)
                .expect("process")
                .relaxation_time_s
        };
        assert_eq!(tau_of("O2").to_bits(), (1.0 / (TAU * fr_o)).to_bits());
        assert_eq!(tau_of("N2").to_bits(), (1.0 / (TAU * fr_n)).to_bits());
    }
}

/// The generic law with ISO's rate variant and ISO's characteristic
/// temperatures against ISO 9613-1's own vibrational terms (classical term
/// excluded on both sides). Logged per point. The remaining differences are
/// model-form terms of the standard: its strengths are the small-strength
/// linearization with `(1 - e^{-theta/T})^2 -> 1`, `c = 343.23 (T/T0)^{1/2}`,
/// fixed undiluted O2/N2 fractions, and `tau = 1/(2 pi f_r)`. The asserted band
/// is the standard's own stated accuracy (+-10 % in its domain); it is not the
/// bead's 1e-3 G1 gate, which stays open and is reported from these numbers.
#[test]
fn iso_rate_variant_vibrational_terms_track_iso_formula() {
    let law = iso_air_law();
    let mut worst_dry = 0.0_f64;
    let mut worst_humid = 0.0_f64;
    for t in [253.15, 273.15, 293.15, 313.15, 323.15] {
        for x_w in [0.0, 0.005, 0.01, 0.02] {
            let p = 101_325.0;
            let gas = test_gas(1.4, t, p);
            let dry = 1.0 - x_w;
            let state = law
                .evaluate(&medium(
                    gas,
                    &[
                        ("N2", 0.781 * dry),
                        ("O2", 0.209 * dry),
                        ("Ar", 0.010 * dry),
                        ("H2O", x_w),
                    ],
                ))
                .expect("state");
            let (fr_o, fr_n) = iso9613_relaxation_frequencies_hz(t, p, 100.0 * x_w).expect("fr");
            for i in 0..13 {
                let f = 50.0 * 200.0_f64.powf(i as f64 / 12.0);
                let f2 = f * f;
                let iso = f2
                    * (t / 293.15).powf(-2.5)
                    * (0.01275 * (-2239.1 / t).exp() / (fr_o + f2 / fr_o)
                        + 0.1068 * (-3352.0 / t).exp() / (fr_n + f2 / fr_n));
                let generic = state.relaxing_wavenumber(TAU * f).expect("k").im;
                let rel = (generic - iso) / iso;
                if x_w == 0.0 {
                    worst_dry = worst_dry.max(rel.abs());
                } else {
                    worst_humid = worst_humid.max(rel.abs());
                }
                println!(
                    "{{\"bead\":\"M0.1-W2\",\"test\":\"iso_vibrational_terms\",\"T_K\":{t},\"x_H2O\":{x_w},\"f_Hz\":{f:.2},\"generic_Np_m\":{generic:e},\"iso_Np_m\":{iso:e},\"relative\":{rel:e}}}"
                );
                assert!(
                    rel.abs() <= 0.10,
                    "outside ISO's stated accuracy at T={t} x_w={x_w} f={f}: {rel:e}"
                );
            }
        }
    }
    // Context row: the full ISO formula (classical included) at 20 C, 50 % RH, 1 kHz.
    let full = iso9613_absorption_neper_per_m(293.15, 101_325.0, 0.5, TAU * 1_000.0).expect("iso");
    println!(
        "{{\"bead\":\"M0.1-W2\",\"test\":\"iso_vibrational_terms\",\"worst_dry\":{worst_dry:e},\"worst_humid\":{worst_humid:e},\"iso_full_20C_50RH_1kHz_Np_m\":{full:e},\"band\":0.10}}"
    );
}

#[test]
fn refuses_species_without_relaxation_data() {
    let gas = test_gas(1.3, 293.15, 101_325.0);
    let law = RelaxationAbsorption::new(vec![species("N2", vec![])]).expect("law");
    let error = law
        .evaluate(&medium(gas, &[("N2", 0.96), ("CO2", 0.04)]))
        .expect_err("CO2 without data");
    assert!(refusal_text(&error).contains("CO2"), "{error}");
    // Zero fraction: absent, so no data is required.
    law.evaluate(&medium(gas, &[("N2", 1.0), ("CO2", 0.0)]))
        .expect("zero fraction is absent");
}

#[test]
fn refuses_collision_partner_without_rate() {
    let p = 101_325.0;
    let gas = test_gas(1.3, 293.15, p);
    let partner = |name: &str| PartnerRate {
        partner: id(name),
        ln_p_tau_pa_s: -4.0,
        landau_teller_slope_k13: 10.0,
        t_min_k: 200.0,
        t_max_k: 300.0,
    };
    let law = RelaxationAbsorption::new(vec![
        species(
            "O2",
            vec![vib(
                2239.1,
                1,
                RelaxationRateLaw::CollisionPartners(vec![partner("O2"), partner("N2")]),
            )],
        ),
        species("N2", vec![]),
        species("He", vec![]),
    ])
    .expect("law");
    let state = law
        .evaluate(&medium(gas, &[("O2", 0.2), ("N2", 0.8)]))
        .expect("all partners rated");
    // Mixture rule: 1/(p tau) = sum x_k/(p tau_k); both partners share one law here.
    let p_tau = (-4.0_f64 + 10.0 * (293.15_f64).powf(-1.0 / 3.0)).exp();
    let expected = p_tau / p;
    let tau = state.processes()[0].relaxation_time_s;
    assert!(
        (tau - expected).abs() <= 1e-12 * expected,
        "{tau} vs {expected}"
    );
    let error = law
        .evaluate(&medium(gas, &[("O2", 0.2), ("He", 0.8)]))
        .expect_err("He partner");
    assert!(refusal_text(&error).contains("He"), "{error}");
}

#[test]
fn refuses_outside_rate_validity() {
    let p = 101_325.0;
    let law = RelaxationAbsorption::new(vec![species(
        "A",
        vec![vib(
            960.0,
            2,
            RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: 1.0,
                t_ref_k: 300.0,
                exponent: -1.0,
                t_min_k: 250.0,
                t_max_k: 300.0,
            },
        )],
    )])
    .expect("law");
    law.evaluate(&medium(test_gas(1.3, 300.0, p), &[("A", 1.0)]))
        .expect("bound passes");
    law.evaluate(&medium(test_gas(1.3, 250.0, p), &[("A", 1.0)]))
        .expect("lower bound passes");
    let error = law
        .evaluate(&medium(test_gas(1.3, 300.0 + 1e-9, p), &[("A", 1.0)]))
        .expect_err("bound + 1e-9 refuses");
    assert!(refusal_text(&error).contains("valid for T"), "{error}");
}

#[test]
fn process_cap_and_cap_plus_one() {
    let p = 101_325.0;
    let build = |count: usize| {
        let names = ["A", "B", "C", "D", "E", "F", "G", "H", "I"];
        let data: Vec<SpeciesRelaxation> = (0..count.div_ceil(4))
            .map(|s| {
                let n = (count - 4 * s).min(4);
                species(
                    names[s],
                    (0..n)
                        .map(|j| vib(1000.0 + 100.0 * j as f64, 1, timed(1e-5, p)))
                        .collect(),
                )
            })
            .collect();
        RelaxationAbsorption::new(data)
    };
    let law = build(MAX_RELAXATION_PROCESSES).expect("cap passes");
    let total: usize = law.species().iter().map(|s| s.processes.len()).sum();
    assert_eq!(total, MAX_RELAXATION_PROCESSES);
    let error = build(MAX_RELAXATION_PROCESSES + 1).expect_err("cap + 1 refuses");
    assert!(
        refusal_text(&error).contains("MAX_RELAXATION_PROCESSES"),
        "{error}"
    );
}

#[test]
fn iso_variant_refuses_non_air() {
    let law = RelaxationAbsorption::new(vec![
        species(
            "O2",
            vec![vib(
                2239.1,
                1,
                RelaxationRateLaw::Iso9613TerrestrialAir(IsoRelaxingSpecies::Oxygen),
            )],
        ),
        species(
            "N2",
            vec![vib(
                3352.0,
                1,
                RelaxationRateLaw::Iso9613TerrestrialAir(IsoRelaxingSpecies::Nitrogen),
            )],
        ),
        species("He", vec![]),
        species("Ar", vec![]),
        species("H2O", vec![]),
        species("CO2", vec![]),
    ])
    .expect("law");
    let gas = test_gas(1.4, 293.15, 101_325.0);
    let heliox = law
        .evaluate(&medium(gas, &[("He", 0.8), ("O2", 0.2)]))
        .expect_err("heliox");
    assert!(refusal_text(&heliox).contains("He"), "{heliox}");
    let breath = law
        .evaluate(&medium(
            gas,
            &[
                ("N2", 0.70),
                ("O2", 0.16),
                ("Ar", 0.0083),
                ("CO2", 0.04),
                ("H2O", 0.0917),
            ],
        ))
        .expect_err("breath");
    assert!(
        refusal_text(&breath).contains("terrestrial air"),
        "{breath}"
    );
    let enriched = law
        .evaluate(&medium(gas, &[("N2", 0.75), ("O2", 0.24), ("Ar", 0.01)]))
        .expect_err("O2-enriched air");
    assert!(
        refusal_text(&enriched).contains("terrestrial air"),
        "{enriched}"
    );
    // The ISO law belongs to the named species' vibration only.
    let misplaced = RelaxationAbsorption::new(vec![species(
        "N2",
        vec![vib(
            3352.0,
            1,
            RelaxationRateLaw::Iso9613TerrestrialAir(IsoRelaxingSpecies::Oxygen),
        )],
    )])
    .expect_err("ISO O2 law on N2");
    assert!(refusal_text(&misplaced).contains("O2"), "{misplaced}");
}

#[test]
fn refuses_overdeclared_heat_capacity_and_bad_composition() {
    let p = 101_325.0;
    // gamma 1.4: c_v = 2.5 R_s; rotation (1 R_s) plus three fully excited
    // modes cannot fit above translation (1.5 R_s).
    let law = RelaxationAbsorption::new(vec![species(
        "A",
        vec![
            rot(2, timed(1e-9, p)),
            RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 10.0,
                    degeneracy: 3,
                }]),
                rate: timed(1e-5, p),
            },
        ],
    )])
    .expect("law");
    let error = law
        .evaluate(&medium(test_gas(1.4, 293.15, p), &[("A", 1.0)]))
        .expect_err("over-declared");
    assert!(refusal_text(&error).contains("translation"), "{error}");
    let gas = test_gas(1.4, 293.15, p);
    assert!(
        RelaxingMedium::new(gas, vec![(id("A"), 0.5), (id("B"), 0.4)]).is_err(),
        "sum 0.9"
    );
    assert!(
        RelaxingMedium::new(gas, vec![(id("A"), 0.5), (id("A"), 0.5)]).is_err(),
        "duplicate"
    );
    assert!(
        RelaxingMedium::new(gas, vec![(id("A"), 1.5), (id("B"), -0.5)]).is_err(),
        "outside [0,1]"
    );
}
