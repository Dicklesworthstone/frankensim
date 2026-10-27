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
    let mut omega = ka * gas.sound_speed / RADIUS;
    // The load recomputes `ka = omega / c * a`; a request at the cap must not round-trip an ulp above it.
    while omega / gas.sound_speed * RADIUS > ka {
        omega = f64::from_bits(omega.to_bits() - 1);
    }
    omega
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
    // One whole unit of ka past the cap, at the production and the largest mode count.
    for modes in [FLANGED_MODES, FLANGED_MAX_MODES] {
        let past =
            flanged_pipe_impedance(&gas, RADIUS, omega_of(&gas, FLANGED_MAX_KA + 1.0), modes);
        assert!(
            matches!(past, Err(DuctError::RadiationKaTooLarge { ka }) if (ka - (FLANGED_MAX_KA + 1.0)).abs() < 1.0e-9),
            "cap + 1 must refuse and name its ka: {past:?}"
        );
    }
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
        b = flanged_load_uncached(&gas, &gas, RADIUS, omega, ka, FLANGED_MODES).expect("fresh");
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

#[test]
fn two_medium_flange_matches_pressure_velocity_continuity() {
    let internal = air();
    let external = GasState::try_new(
        &fs_material::gas_species::GasSpecies::Helium.spec(), 293.15, internal.pressure,
    ).unwrap();
    let alpha = RULES.as_ref().unwrap().roots[0].0;
    let area = PI * RADIUS * RADIUS;
    for internal_ka in [0.01, 0.7, alpha * 0.999, alpha * 1.001, 7.0] {
        let omega = omega_of(&internal, internal_ka);
        let external_ka = omega / external.sound_speed * RADIUS;
        let (r, x) = fs_phs::baffled_piston_impedance(
            external.density, external.sound_speed, RADIUS, omega,
        ).unwrap();
        let piston = C64::new(r, x);
        assert_eq!(flanged_pipe_impedance_with_exterior(
            &internal, &external, RADIUS, omega, 1,
        ).unwrap(), piston);

        // Solve the two pressure/velocity equations directly in dimensional
        // modal coordinates. Only the interior axial wavenumber changes at
        // duct-mode cutoff; the exterior radiation matrix remains smooth.
        let z = radiation_matrix(external_ka, 2).unwrap();
        let exterior_rho_c = external.density * external.sound_speed;
        let q_squared = (omega / internal.sound_speed).powi(2) - (alpha / RADIUS).powi(2);
        let internal_mode = if q_squared > 0.0 {
            C64::from_re(internal.density * omega / det::sqrt(q_squared))
        } else {
            C64::new(0.0, -internal.density * omega / det::sqrt(-q_squared))
        };
        let z01 = z[1].scale(exterior_rho_c);
        let z11 = z[2].scale(exterior_rho_c);
        let expected = piston - (z01 * z01 / (z11 + internal_mode)).scale(1.0 / area);
        let actual = flanged_pipe_impedance_with_exterior(
            &internal, &external, RADIUS, omega, 2,
        ).unwrap();
        assert!((actual - expected).abs() < 2.0e-12 * expected.abs());
        assert_ne!(actual, flanged_pipe_impedance(&external, RADIUS, omega, 2).unwrap(),
            "the interior higher modes cannot be replaced by exterior gas");

        let same = flanged_pipe_impedance_with_exterior(
            &internal, &internal, RADIUS, omega, FLANGED_MODES,
        ).unwrap();
        assert_eq!(same, flanged_pipe_impedance(&internal, RADIUS, omega, FLANGED_MODES).unwrap());
        let production = flanged_pipe_impedance_with_exterior(
            &internal, &external, RADIUS, omega, FLANGED_MODES,
        ).unwrap();
        assert!(production.re >= 0.0);
        let zc = C64::from_re(internal.density * internal.sound_speed / area);
        assert!(((production - zc) / (production + zc)).abs() <= 1.0 + 1.0e-12);
    }
    let wrong_pressure = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 50000.0).unwrap();
    assert!(flanged_pipe_impedance_with_exterior(
        &internal, &wrong_pressure, RADIUS, 1000.0, FLANGED_MODES,
    ).is_err());
}

/// Measurement behind the cap (bead `frankensim-irp-v5-instrument-realism-pexfp.324`). The cap is a
/// quadrature-design limit, not an edge of the model, so what matters is where the design stops
/// resolving. This calls the cap-free seams across and above `FLANGED_MAX_KA`: the quadrature
/// `Z_00` against the closed-form piston, the 16-mode load against the 24-mode load, and
/// passivity. It records each `ka` (an error is recorded, not unwrapped) and asserts that every
/// measured point at or below the cap holds the 1e-12 gate, that the cap itself is measured, and
/// that the first point above the cap still holds the gate, so the cap is never the edge itself.
/// Measured natively (qtest f67a93ae9): 3.5e-14 at 20, 6.5e-14 at 24, 6.0e-14 at 26, 2.9e-13 at
/// 28, 2.0e-11 at 30. Moving the cap past 28 fails the last assertion. Receipts:
/// `IRPDUCT_FLANGED_ABOVE_CAP`.
#[test]
fn flanged_design_is_measured_across_and_beyond_the_cap() {
    let gas = air();
    let z0 = gas.density * gas.sound_speed / (PI * RADIUS * RADIUS);
    let mut at_cap = None;
    let mut first_above = None;
    for ka in [20.0, 20.5, 21.0, 22.0, 24.0, 26.0, 28.0, 30.0, 34.0, 38.0] {
        let omega = omega_of(&gas, ka);
        let (r, x) = fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, RADIUS, omega)
            .expect("piston");
        let closed = C64::new(r / z0, x / z0);
        let quad = radiation_matrix(ka, 1).expect("matrix")[0];
        let piston_rel = (quad - closed).abs() / closed.abs();
        let l16 = flanged_load_uncached(&gas, &gas, RADIUS, omega, ka, FLANGED_MODES);
        let l24 = flanged_load_uncached(&gas, &gas, RADIUS, omega, ka, FLANGED_MAX_MODES);
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
            (a, b) => format!(
                "REFUSED 16: {:?} 24: {:?}",
                a.as_ref().err(),
                b.as_ref().err()
            ),
        };
        eprintln!("IRPDUCT_FLANGED_ABOVE_CAP ka={ka} piston_rel={piston_rel:e} {line}");
        if ka == FLANGED_MAX_KA {
            at_cap = Some(piston_rel);
        }
        if ka > FLANGED_MAX_KA && first_above.is_none() {
            first_above = Some((ka, piston_rel));
        }
        assert!(
            piston_rel.is_finite(),
            "quadrature Z_00 not finite at ka={ka}"
        );
        if ka <= FLANGED_MAX_KA {
            assert!(
                piston_rel < 1.0e-12,
                "admitted ka={ka} loses the 1e-12 gate: {piston_rel:e}"
            );
        }
    }
    assert!(at_cap.expect("the cap is one of the measured points") < 1.0e-12);
    let (ka_above, rel_above) = first_above.expect("a measured point lies above the cap");
    assert!(
        rel_above < 1.0e-12,
        "the cap must sit inside the plateau: ka={ka_above} just above it already loses the gate ({rel_above:e})"
    );
}

/// The band the lift admitted, `ka` from 20 to the cap, held to the module's own gates on a dense
/// linear grid that straddles the one cutoff inside it (`α_7 = 22.76`): the quadrature `Z_00`
/// against the closed-form piston at 1e-12, and a finite, passive load at the production mode
/// count. A mode exactly at cutoff drops out of the Schur complement, so the grid also lands the
/// cutoff to within a few ulp.
#[test]
fn flanged_lifted_band_holds_the_gates_across_a_cutoff() {
    let gas = air();
    let z0 = gas.density * gas.sound_speed / (PI * RADIUS * RADIUS);
    let alpha7 = RULES.as_ref().expect("rules").roots[6].0;
    assert!(
        (20.0..FLANGED_MAX_KA).contains(&alpha7),
        "the cutoff this test straddles must lie inside the lifted band: {alpha7}"
    );
    let steps = ((FLANGED_MAX_KA - 20.0) / 0.05).round() as usize;
    let mut grid: Vec<f64> = (0..=steps).map(|i| 20.0 + 0.05 * i as f64).collect();
    grid.extend([
        alpha7 * (1.0 - 1.0e-6),
        alpha7 * (1.0 - 1.0e-9),
        alpha7,
        alpha7 * (1.0 + 1.0e-9),
        alpha7 * (1.0 + 1.0e-6),
    ]);
    grid.retain(|&ka| ka < FLANGED_MAX_KA);
    let (mut worst_piston, mut worst_piston_ka, mut worst_r, mut worst_r_ka) =
        (0.0_f64, 0.0, 0.0_f64, 0.0);
    for &ka in &grid {
        let omega = omega_of(&gas, ka);
        let (r, x) = fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, RADIUS, omega)
            .expect("piston");
        let closed = C64::new(r / z0, x / z0);
        let piston_rel =
            (radiation_matrix(ka, 1).expect("matrix")[0] - closed).abs() / closed.abs();
        let z = flanged_pipe_impedance(&gas, RADIUS, omega, FLANGED_MODES).expect("admitted load");
        let refl = (C64::new(z.re - z0, z.im) / C64::new(z.re + z0, z.im)).abs();
        assert!(
            z.re.is_finite() && z.im.is_finite(),
            "load not finite at ka={ka}"
        );
        assert!(refl <= 1.0, "load not passive at ka={ka}: |R|={refl}");
        if piston_rel > worst_piston {
            (worst_piston, worst_piston_ka) = (piston_rel, ka);
        }
        if refl > worst_r {
            (worst_r, worst_r_ka) = (refl, ka);
        }
    }
    assert!(
        worst_piston < 1.0e-12,
        "piston gate lost at ka={worst_piston_ka}: {worst_piston:e}"
    );
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"flanged_lifted_band_holds_the_gates_across_a_cutoff\",\
         \"inputs\":\"{} points, ka 20..{FLANGED_MAX_KA} step 0.05 plus alpha_7={alpha7:.6} +-1e-6,1e-9,0\",\
         \"observed\":{{\"worst_piston_rel\":{worst_piston:e},\"at_ka\":{worst_piston_ka},\
         \"max_abs_r\":{worst_r:.12},\"at_ka_r\":{worst_r_ka}}},\"verdict\":\"pass\"}}",
        grid.len()
    );
}
