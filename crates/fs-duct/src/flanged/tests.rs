//! Battery for the multimodal flanged-pipe load (bead
//! `frankensim-irp-v5-instrument-realism-pexfp.355`). Every case prints one
//! stderr JSON receipt.

use super::*;
use fs_material::gas::GasSpec;

const PI: f64 = core::f64::consts::PI;
const RADIUS: f64 = 0.01;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

fn omega_of(gas: &GasState, ka: f64) -> f64 {
    ka * gas.sound_speed / RADIUS
}

/// Dimensionless load `Z_L/z0` at `ka`.
fn load(gas: &GasState, ka: f64, modes: usize) -> C64 {
    let z0 = gas.density * gas.sound_speed / (PI * RADIUS * RADIUS);
    flanged_pipe_impedance(gas, RADIUS, omega_of(gas, ka), modes)
        .expect("flanged load")
        .scale(1.0 / z0)
}

fn log_grid(lo: f64, hi: f64, count: usize) -> impl Iterator<Item = f64> {
    let span = det::ln(hi / lo);
    (0..count).map(move |i| {
        if i + 1 == count {
            hi
        } else {
            lo * det::exp(span * i as f64 / (count - 1) as f64)
        }
    })
}

/// G2: the quadrature's own `Z_00` reproduces the closed-form piston. That
/// certifies the seam every `Z_mn` uses: the real part, the near panels and
/// both halves of the rotated tail. With one mode the load is the piston.
#[test]
fn flanged_n1_is_the_closed_form_piston() {
    let start = std::time::Instant::now();
    let gas = air();
    let z0 = gas.density * gas.sound_speed / (PI * RADIUS * RADIUS);
    let mut worst = 0.0_f64;
    let mut worst_ka = 0.0_f64;
    for ka in log_grid(1.0e-3, FLANGED_MAX_KA, 81) {
        let (r, x) =
            fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, RADIUS, omega_of(&gas, ka))
                .expect("piston");
        let closed = C64::new(r / z0, x / z0);
        let quad = radiation_matrix(ka, 1).expect("matrix")[0];
        let rel = (quad - closed).abs() / closed.abs();
        if rel > worst {
            worst = rel;
            worst_ka = ka;
        }
        let one = flanged_pipe_impedance(&gas, RADIUS, omega_of(&gas, ka), 1).expect("N=1");
        assert_eq!(one.re.to_bits(), r.to_bits(), "N=1 re at ka={ka}");
        assert_eq!(one.im.to_bits(), x.to_bits(), "N=1 im at ka={ka}");
    }
    let verdict = worst < 1.0e-12;
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"flanged_n1_is_the_closed_form_piston\",\"axis\":\"ka\",\
         \"inputs\":\"81-point log grid 1e-3..{FLANGED_MAX_KA}\",\"observed\":{worst:e},\"at_ka\":{worst_ka:e},\
         \"oracle\":\"fs_phs::baffled_piston_impedance\",\"tol\":1e-12,\"verdict\":\"{}\",\"elapsed_ms\":{}}}",
        if verdict { "pass" } else { "fail" },
        start.elapsed().as_millis()
    );
    assert!(verdict, "quadrature Z_00 vs closed form: {worst:e} at ka={worst_ka}");
}

/// G2: the low-`ka` end correction `δ/a = −Im(Z/ρc)/ka` on a mode ladder.
/// One mode is the piston's `8/(3π)`; the ladder moves it to the
/// flanged-pipe `0.8216` (Norris & Sheng 1989; Zorumski 1973).
#[test]
fn flanged_mode_ladder_converges_to_0p8216() {
    let start = std::time::Instant::now();
    let gas = air();
    let ka = 1.0e-3;
    let rungs = [1_usize, 2, 4, 8, 16, 24];
    let delta: Vec<f64> = rungs.iter().map(|&n| -load(&gas, ka, n).im / ka).collect();
    let finest = delta[rungs.len() - 1];
    for (i, (&n, &d)) in rungs.iter().zip(&delta).enumerate() {
        let order = if i + 2 < rungs.len() {
            let e0 = (delta[i] - finest).abs();
            let e1 = (delta[i + 1] - finest).abs();
            det::ln(e0 / e1) / det::ln(rungs[i + 1] as f64 / n as f64)
        } else {
            f64::NAN
        };
        eprintln!(
            "{{\"suite\":\"fs-duct\",\"case\":\"flanged_mode_ladder\",\"axis\":\"modes\",\"rung\":{n},\
             \"ka\":{ka:e},\"observed\":{d:.9},\"err_vs_0p8216\":{:e},\"err_vs_finest\":{:e},\
             \"observed_order\":{order:.3}}}",
            d - 0.8216,
            d - finest
        );
    }
    assert!((delta[0] - 8.0 / (3.0 * PI)).abs() < 1.0e-6, "N=1 is the piston: {}", delta[0]);
    assert!(
        delta.windows(2).all(|w| w[1] < w[0]),
        "every added mode must shorten the end correction: {delta:?}"
    );
    let production = -load(&gas, ka, FLANGED_MODES).im / ka;
    let verdict = (finest - 0.8216).abs() < 1.0e-3 && (production - 0.8216).abs() < 1.0e-3;
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"flanged_mode_ladder_converges_to_0p8216\",\
         \"observed\":{{\"finest\":{finest:.9},\"production\":{production:.9}}},\"oracle\":0.8216,\
         \"tol\":1e-3,\"verdict\":\"{}\",\"elapsed_ms\":{}}}",
        if verdict { "pass" } else { "fail" },
        start.elapsed().as_millis()
    );
    assert!(verdict, "finest {finest} / production {production} vs 0.8216");
}

/// G1: one model, no switch. Second differences on a fine `ka` grid around
/// the old `ka = 1` fit/piston boundary stay `O(h²)`, and the two sides of
/// `ka = 1` agree. The old switch jumped by about 0.19 there.
#[test]
fn flanged_load_is_continuous_through_the_old_ka_switch() {
    let start = std::time::Instant::now();
    let gas = air();
    let h = 1.0e-3;
    let zs: Vec<C64> = (0..=1000)
        .map(|i| load(&gas, 0.5 + h * f64::from(i), FLANGED_MODES))
        .collect();
    let scale = zs.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
    let curvature = zs
        .windows(3)
        .map(|w| (w[2] - w[1].scale(2.0) + w[0]).abs())
        .fold(0.0_f64, f64::max)
        / scale;
    let eps = 1.0e-9;
    let jump = (load(&gas, 1.0 + eps, FLANGED_MODES) - load(&gas, 1.0 - eps, FLANGED_MODES)).abs()
        / scale;
    let verdict = curvature < 1.0e-5 && jump < 1.0e-8;
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"flanged_load_is_continuous\",\"axis\":\"ka\",\
         \"observed\":{{\"second_difference\":{curvature:e},\"jump_at_1\":{jump:e}}},\
         \"tol\":{{\"second_difference\":1e-5,\"jump_at_1\":1e-8}},\"verdict\":\"{}\",\"elapsed_ms\":{}}}",
        if verdict { "pass" } else { "fail" },
        start.elapsed().as_millis()
    );
    assert!(verdict, "second difference {curvature:e}, jump {jump:e}");
}

/// G2: passive over the whole admitted band. `Re Z ≥ 0`, and the matched-tube
/// reflectance has `|R| ≤ 1`.
#[test]
fn flanged_load_is_passive_over_the_admitted_band() {
    let start = std::time::Instant::now();
    let gas = air();
    let mut worst_r = 0.0_f64;
    for ka in log_grid(1.0e-4, FLANGED_MAX_KA, 401) {
        let z = load(&gas, ka, FLANGED_MODES);
        assert!(z.re >= 0.0, "Re Z = {} < 0 at ka={ka}", z.re);
        let r = ((z - C64::ONE) / (z + C64::ONE)).abs();
        worst_r = worst_r.max(r);
        assert!(r <= 1.0 + 1.0e-12, "|R| = {r} > 1 at ka={ka}");
    }
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"flanged_load_is_passive\",\"axis\":\"ka\",\
         \"inputs\":\"401-point log grid 1e-4..{FLANGED_MAX_KA}\",\"observed\":{{\"max_abs_r\":{worst_r:.12}}},\
         \"tol\":\"Re Z >= 0, |R| <= 1+1e-12\",\"verdict\":\"pass\",\"elapsed_ms\":{}}}",
        start.elapsed().as_millis()
    );
}

/// G0: refusals by name, at and just past each cap.
#[test]
fn flanged_refusals_fire_by_name() {
    let gas = air();
    let at_cap = omega_of(&gas, FLANGED_MAX_KA);
    assert!(flanged_pipe_impedance(&gas, RADIUS, at_cap, FLANGED_MODES).is_ok());
    assert!(matches!(
        flanged_pipe_impedance(&gas, RADIUS, at_cap * (1.0 + 1.0e-12), FLANGED_MODES),
        Err(DuctError::RadiationKaTooLarge { .. })
    ));
    assert!(flanged_pipe_impedance(&gas, RADIUS, omega_of(&gas, 0.5), FLANGED_MAX_MODES).is_ok());
    for modes in [0, FLANGED_MAX_MODES + 1] {
        assert!(matches!(
            flanged_pipe_impedance(&gas, RADIUS, omega_of(&gas, 0.5), modes),
            Err(DuctError::BadParameter { .. })
        ));
    }
    for (radius, omega) in [(0.0, 1.0), (RADIUS, 0.0), (f64::NAN, 1.0), (RADIUS, f64::INFINITY)] {
        assert!(matches!(
            flanged_pipe_impedance(&gas, radius, omega, FLANGED_MODES),
            Err(DuctError::BadParameter { .. })
        ));
    }
}

/// G5: bitwise repeat, a memo hit equals a fresh computation, and the
/// uncached per-call cost is on record.
#[test]
fn flanged_load_bitwise_repeat() {
    let gas = air();
    let omega = omega_of(&gas, 2.3);
    // The same ka the public entry derives from omega.
    let ka = omega / gas.sound_speed * RADIUS;
    let a = flanged_pipe_impedance(&gas, RADIUS, omega, FLANGED_MODES).expect("first");
    let hit = flanged_pipe_impedance(&gas, RADIUS, omega, FLANGED_MODES).expect("memo hit");
    let start = std::time::Instant::now();
    let calls = 50_u32;
    let mut b = C64::ZERO;
    for _ in 0..calls {
        b = flanged_load_uncached(&gas, RADIUS, omega, ka, FLANGED_MODES).expect("fresh");
    }
    let per_call_us = start.elapsed().as_secs_f64() * 1.0e6 / f64::from(calls);
    assert_eq!((a.re.to_bits(), a.im.to_bits()), (b.re.to_bits(), b.im.to_bits()));
    assert_eq!((a.re.to_bits(), a.im.to_bits()), (hit.re.to_bits(), hit.im.to_bits()));
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"flanged_load_bitwise_repeat\",\"observed\":{{\"re\":{:e},\"im\":{:e}}},\
         \"modes\":{FLANGED_MODES},\"per_call_us\":{per_call_us:.1},\"verdict\":\"pass\"}}",
        a.re, a.im
    );
}

/// Measurement behind the cap (bead `frankensim-irp-v5-instrument-realism-pexfp.324`). The cap is a
/// quadrature-design limit, not an edge of the model, so what matters is where the design stops
/// resolving. This calls the cap-free seams above `FLANGED_MAX_KA`: the quadrature `Z_00` against the
/// closed-form piston, the 16-mode load against the 24-mode load, and passivity. It records each
/// `ka` (an error is recorded, not unwrapped) and asserts only that the cap itself is inside what
/// it measures. Receipts: `IRPDUCT_FLANGED_ABOVE_CAP`.
#[test]
fn flanged_design_is_measured_above_the_cap() {
    let gas = air();
    let z0 = gas.density * gas.sound_speed / (PI * RADIUS * RADIUS);
    let mut at_cap = None;
    for ka in [FLANGED_MAX_KA, 20.5, 21.0, 22.0, 24.0, 26.0, 28.0, 30.0, 34.0, 38.0] {
        let omega = omega_of(&gas, ka);
        let (r, x) = fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, RADIUS, omega)
            .expect("piston");
        let closed = C64::new(r / z0, x / z0);
        let quad = radiation_matrix(ka, 1).expect("matrix")[0];
        let piston_rel = (quad - closed).abs() / closed.abs();
        let l16 = flanged_load_uncached(&gas, RADIUS, omega, ka, FLANGED_MODES);
        let l24 = flanged_load_uncached(&gas, RADIUS, omega, ka, FLANGED_MAX_MODES);
        let line = match (&l16, &l24) {
            (Ok(a), Ok(b)) => {
                let ladder = (*a - *b).abs() / b.abs();
                let refl = |z: &C64| {
                    let num = C64::new(z.re - z0, z.im);
                    let den = C64::new(z.re + z0, z.im);
                    (num / den).abs()
                };
                format!(
                    "modes16_vs_24={ladder:e} abs_R_16={:.12} abs_R_24={:.12} re_16={:e}",
                    refl(a),
                    refl(b),
                    a.re / z0
                )
            }
            (a, b) => format!("REFUSED 16: {:?} 24: {:?}", a.as_ref().err(), b.as_ref().err()),
        };
        eprintln!("IRPDUCT_FLANGED_ABOVE_CAP ka={ka} piston_rel={piston_rel:e} {line}");
        if ka == FLANGED_MAX_KA {
            at_cap = Some(piston_rel);
        }
        assert!(piston_rel.is_finite(), "quadrature Z_00 not finite at ka={ka}");
    }
    assert!(at_cap.expect("the cap is measured") < 1.0e-12, "the cap itself must hold the 1e-12 gate");
}
