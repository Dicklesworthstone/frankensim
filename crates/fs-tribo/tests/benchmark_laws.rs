//! M4.2-W4 benchmark tangential laws (IRP pexfp.150): classic oracles for
//! rate-and-state (Dieterich 1979; Ruina 1983), LuGre (Canudas de Wit et
//! al. 1995) and the Smith & Woodhouse (2000) plastic thermal model of
//! rosin. Gauntlet tiers are named per test; each prints JSON lines.
//! The rate-and-state and LuGre parameters are synthetic fixture values,
//! not material data. The Smith & Woodhouse fixtures are that paper's
//! published inputs (Table 1, Eq. 2, sections 5.1 and 5.3), cited where
//! they are defined.

use fs_math::det;
use fs_rand::StreamKey;
use fs_tribo::benchmark::{
    Evolution, LuGre, LuGreState, MAX_HISTORY_STEPS, RateState, RateStateState, RosinContact,
    SmithWoodhouse, SmithWoodhouseParams, SteadySlidingFit, ThermalSolid, half_space_unit_weight,
};
use fs_tribo::{
    FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput, TriboError,
};
use std::fs;
use std::path::{Path, PathBuf};

const NORMAL_N: f64 = 100.0;

fn solve(max_iterations: usize) -> TangentialSolveConfig {
    TangentialSolveConfig {
        max_iterations,
        maximum_force_n: 1.0e4,
        force_absolute_tolerance_n: 1e-12,
        force_relative_tolerance: 1e-14,
    }
}
fn step(slip_m: f64, dt_s: f64, compliance: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: NORMAL_N,
        slip_increment_m: [slip_m, 0.0],
        dt_s,
        contact_temperature_k: None,
        nominal_pressure_pa: None,
        compliance_m_per_n: compliance,
        solve: solve(128),
    }
}
fn rate_state(evolution: Evolution) -> RateState {
    RateState {
        mu0: 0.6,
        a: 0.010,
        b: 0.015,
        v0_mps: 1e-6,
        dc_m: 1e-5,
        evolution,
        speed_range_mps: [1e-15, 1.0],
    }
}
fn lugre() -> LuGre {
    LuGre {
        sigma0_n_per_m: 1e5,
        sigma1_ns_per_m: 316.2,
        sigma2_ns_per_m: 0.4,
        coulomb_n: 1.0,
        static_n: 1.5,
        stribeck_mps: 1e-3,
    }
}

#[test]
fn g1_rate_state_steady_state_log_law() {
    // Prescribed speed for 40 Dc of slip: the exact evolution leaves
    // e^{-40} of the initial state offset, so mu reaches mu_ss.
    let mut worst = 0.0_f64;
    for evolution in [Evolution::Aging, Evolution::Slip] {
        let law = rate_state(evolution);
        for v in [1e-8, 1e-6, 1e-4, 1e-2] {
            let dt = law.dc_m / v;
            let mut s = law.initial_state();
            for _ in 0..40 {
                s = law.advance(&s, &step(v * dt, dt, 0.0)).unwrap().state;
            }
            let mu = law.advance(&s, &step(v * dt, dt, 0.0)).unwrap().traction_n[0] / NORMAL_N;
            let oracle = law.mu0 + (law.a - law.b) * det::ln(v / law.v0_mps);
            let e = (mu - oracle).abs() / oracle.abs();
            worst = worst.max(e);
            println!(
                "{{\"law\":\"{}\",\"test\":\"rate_state_steady\",\"v_mps\":{v:e},\"mu\":{mu},\"oracle\":{oracle},\"rel_err\":{e:e},\"tol\":1e-10}}",
                law.model_id()
            );
        }
    }
    assert!(worst < 1e-10, "steady log law {worst:e}");
}

#[test]
fn g1_rate_state_velocity_step_transient_length() {
    // A step from V1 to V2 after steady sliding at V1. The direct effect is
    // the immediate jump a ln(V2/V1) (its sign is the sign of a); the state
    // then relaxes over slip distance Dc. The e-folding slip is measured on a
    // three-rung step ladder by log-linear interpolation between samples.
    let (v1, v2) = (1e-6, 1e-5);
    for evolution in [Evolution::Aging, Evolution::Slip] {
        let law = rate_state(evolution);
        let mut measured = Vec::new();
        for per_dc in [12.5_f64, 25.0, 50.0] {
            let dt = law.dc_m / (v2 * per_dc);
            let mut s = RateStateState::new(law.dc_m / v1).unwrap();
            let mu_before = law.coefficient(v1, s.theta_s()).unwrap();
            let theta_ss = law.dc_m / v2;
            let first = law.advance(&s, &step(v2 * dt, dt, 0.0)).unwrap();
            let direct = first.traction_n[0] / NORMAL_N - mu_before;
            assert!(
                (direct - law.a * det::ln(v2 / v1)).abs() < 1e-12,
                "direct effect a ln(V2/V1): {direct}"
            );
            // The state deviation decays exactly as e^{-slip/Dc} in ln-space
            // for the slip law and linearly for the aging law.
            let deviation = |theta: f64| match evolution {
                Evolution::Aging => (theta - theta_ss).abs(),
                Evolution::Slip => det::ln(v2 * theta / law.dc_m).abs(),
            };
            let d0 = deviation(s.theta_s());
            let target = d0 * det::exp(-1.0);
            let mut slip = 0.0;
            let mut found = None;
            for _ in 0..10_000 {
                let next = law.advance(&s, &step(v2 * dt, dt, 0.0)).unwrap();
                let (da, db) = (deviation(s.theta_s()), deviation(next.state.theta_s()));
                if db <= target {
                    let frac = (det::ln(da) - det::ln(target)) / (det::ln(da) - det::ln(db));
                    found = Some(slip + frac * v2 * dt);
                    break;
                }
                slip += v2 * dt;
                s = next.state;
            }
            let length = found.expect("transient decays below 1/e");
            measured.push(length / law.dc_m);
        }
        println!(
            "{{\"law\":\"{}\",\"test\":\"velocity_step\",\"rungs_per_dc\":[12.5,25,50],\"efold_slip_over_dc\":{measured:?},\"tol\":0.01}}",
            law.model_id()
        );
        for m in &measured {
            assert!((m - 1.0).abs() < 0.01, "e-folding slip {m} Dc");
        }
    }
}

/// Amplitude of `ln(V/V_lp)` in an early and a late slip window of a
/// quasi-static spring-slider run through the compliant `advance` path, or
/// `None` when the run leaves the admitted speed range (a slip event).
fn spring_slider(law: &RateState, k: f64, v_lp: f64) -> Option<(f64, f64)> {
    let dt = 2.5e-3 * law.dc_m / v_lp;
    let mut s = RateStateState::new(1.001 * law.dc_m / v_lp).unwrap();
    let mut force = NORMAL_N * law.steady_coefficient(v_lp).unwrap();
    let (mut early, mut late) = (0.0_f64, 0.0_f64);
    let total = 300.0 * law.dc_m;
    let mut slip = 0.0;
    while slip < total {
        // F_new = F_old + k (dx_lp - du): released motion dx_lp + F_old/k.
        // The residual's slope here is 1 + N (a/V)(c/dt) ~ 800, so one ulp of
        // a 60 N force moves it by ~6e-12 N: the root budget must be looser
        // than that (1e-9 N, ~1e-11 relative), unlike the stiff-free fixtures.
        let mut x = step(v_lp * dt + force / k, dt, 1.0 / k);
        x.solve.force_absolute_tolerance_n = 1e-9;
        let next = match law.advance(&s, &x) {
            Ok(next) => next,
            Err(TriboError::OutsideApplicability { .. }) => return None,
            Err(e) => panic!("spring-slider step: {e:?}"),
        };
        force = next.traction_n[0];
        let du = next.slip_increment_m[0];
        slip += du;
        let amp = det::ln(du / dt / v_lp).abs();
        if slip > 20.0 * law.dc_m && slip < 60.0 * law.dc_m {
            early = early.max(amp);
        }
        if slip > total - 40.0 * law.dc_m {
            late = late.max(amp);
        }
        s = next.state;
    }
    Some((early, late))
}

#[test]
fn g1_spring_slider_critical_stiffness() {
    // Ruina (1983): quasi-static steady sliding is unstable for
    // k < k_c = N (b - a)/Dc. A stiffness sweep must bracket k_c within 2 %.
    let law = rate_state(Evolution::Aging);
    let k_c = law.critical_stiffness_n_per_m(NORMAL_N);
    let mut rows = Vec::new();
    for ratio in [0.90, 0.94, 0.98, 1.02, 1.06, 1.10] {
        let outcome = spring_slider(&law, ratio * k_c, law.v0_mps);
        let (unstable, amplitudes) = match outcome {
            None => (true, "null".to_string()),
            Some((early, late)) => (late > early, format!("[{early},{late}]")),
        };
        rows.push(format!(
            "{{\"k_over_kc\":{ratio},\"early_late_amplitude\":{amplitudes},\"unstable\":{unstable}}}"
        ));
        assert_eq!(unstable, ratio < 1.0, "k/k_c = {ratio}: {outcome:?}");
    }
    println!(
        "{{\"law\":\"{}\",\"test\":\"spring_slider\",\"k_c_n_per_m\":{k_c},\"sweep\":[{}],\"bracket\":[0.98,1.02]}}",
        law.model_id(),
        rows.join(",")
    );
}

#[test]
fn g1_lugre_steady_state_curve() {
    let law = lugre();
    let mut worst = 0.0_f64;
    for v in [-0.1_f64, -1e-2, -1e-3, -1e-4, 1e-4, 1e-3, 1e-2, 0.1] {
        // Each step relaxes the bristle by e^{-40} or more, so z is steady.
        let dt = 40.0 * law.static_n / (law.sigma0_n_per_m * v.abs());
        let mut s = law.initial_state();
        let mut force = 0.0;
        for _ in 0..3 {
            let next = law.advance(&s, &step(v * dt, dt, 0.0)).unwrap();
            force = next.traction_n[0];
            s = next.state;
        }
        let oracle = law.steady_force_n(v);
        let e = (force - oracle).abs() / oracle.abs();
        worst = worst.max(e);
        println!(
            "{{\"law\":\"benchmark-lugre\",\"test\":\"steady\",\"v_mps\":{v},\"force_n\":{force},\"oracle_n\":{oracle},\"rel_err\":{e:e},\"tol\":1e-10}}"
        );
    }
    assert!(worst < 1e-10, "LuGre steady curve {worst:e}");
}

#[test]
fn g1_lugre_presliding_stiffness() {
    // A 1e-14 m quasi-static displacement from rest: F/x -> sigma0.
    let law = lugre();
    let (delta, dt) = (1e-14, 1e8);
    let next = law
        .advance(&law.initial_state(), &step(delta, dt, 0.0))
        .unwrap();
    let stiffness = next.traction_n[0] / delta;
    let e = (stiffness - law.sigma0_n_per_m).abs() / law.sigma0_n_per_m;
    println!(
        "{{\"law\":\"benchmark-lugre\",\"test\":\"presliding\",\"stiffness\":{stiffness},\"rel_err\":{e:e},\"tol\":1e-8}}"
    );
    assert!(e < 1e-8, "presliding stiffness {e:e}");
}

#[test]
fn g0_lugre_passivity_seeded_histories() {
    // Canudas de Wit et al. (1995): the map v -> z is passive,
    // int z v dt >= (z_T^2 - z_0^2)/2, for every history; checked per
    // interval with the exact interval integral.
    let law = lugre();
    let mut worst = f64::INFINITY;
    let mut intervals = 0;
    for history in 0..64_u32 {
        let mut stream = StreamKey {
            seed: 0x7a1b_0150,
            kernel: 150,
            tile: history,
        }
        .stream();
        let mut z = 0.0;
        for _ in 0..200 {
            let v = 0.02 * (stream.next_f64() - 0.5);
            let dt = 1e-4 + 1e-3 * stream.next_f64();
            let z1 = law.bristle_after(z, v, dt);
            let work = law.bristle_work(z, v, dt);
            let storage = 0.5 * (z1 * z1 - z * z);
            let scale = work.abs().max(storage.abs()).max(1e-30);
            worst = worst.min((work - storage) / scale);
            z = z1;
            intervals += 1;
        }
    }
    println!(
        "{{\"law\":\"benchmark-lugre\",\"test\":\"z_map_passivity\",\"intervals\":{intervals},\"min_normalized_margin\":{worst:e}}}"
    );
    assert!(worst >= -1e-12, "z-map passivity margin {worst:e}");
}

#[test]
fn g1_lugre_drive_spring_stick_slip_exists() {
    // Mass m on the LuGre contact, pulled through a spring k by a point at
    // speed v_d (semi-implicit Euler: compliance dt^2/m). The Stribeck drop
    // must produce a stick-slip limit cycle: sticks, and periodic onsets.
    let law = lugre();
    let (m, k, v_d, dt) = (1.0, 2.0, 0.1, 1e-4);
    let (mut x, mut v, mut t) = (0.0_f64, 0.0_f64, 0.0_f64);
    let mut s = law.initial_state();
    let mut onsets = Vec::new();
    let mut stuck_since_onset = false;
    while t < 40.0 {
        let spring = k * (v_d * t - x);
        let input = TangentialStepInput {
            slip_increment_m: [dt * v + dt * dt * spring / m, 0.0],
            compliance_m_per_n: dt * dt / m,
            dt_s: dt,
            ..step(0.0, dt, 0.0)
        };
        let next = law.advance(&s, &input).unwrap();
        let du = next.slip_increment_m[0];
        let v_new = du / dt;
        if v_new.abs() < 1e-3 * v_d {
            stuck_since_onset = true;
        }
        if stuck_since_onset && v < v_d && v_new >= v_d {
            onsets.push(t);
            stuck_since_onset = false;
        }
        x += du;
        v = v_new;
        t += dt;
        s = next.state;
    }
    let periods: Vec<f64> = onsets.windows(2).map(|w| w[1] - w[0]).collect();
    println!(
        "{{\"law\":\"benchmark-lugre\",\"test\":\"drive_spring_stick_slip\",\"onsets\":{},\"periods_s\":{periods:?}}}",
        onsets.len()
    );
    assert!(periods.len() >= 3, "stick-slip cycles: {periods:?}");
    let (p1, p2) = (periods[periods.len() - 2], periods[periods.len() - 1]);
    assert!(
        (p1 - p2).abs() < 0.01 * p2,
        "limit cycle period settles: {p1} vs {p2}"
    );
}

#[test]
fn g0_parameter_cap_and_cap_plus_one() {
    let base = rate_state(Evolution::Aging);
    assert!(RateState { a: 0.0, ..base }.validate().is_ok());
    assert!(
        RateState {
            a: 0.0_f64.next_down(),
            ..base
        }
        .validate()
        .is_err()
    );
    assert!(RateState { b: 0.0, ..base }.validate().is_ok());
    assert!(
        RateState {
            b: 0.0_f64.next_down(),
            ..base
        }
        .validate()
        .is_err()
    );
    assert!(
        RateState {
            dc_m: f64::MIN_POSITIVE,
            ..base
        }
        .validate()
        .is_ok()
    );
    assert!(RateState { dc_m: 0.0, ..base }.validate().is_err());
    assert!(
        RateState {
            speed_range_mps: [1e-9, 1e-9_f64.next_up()],
            ..base
        }
        .validate()
        .is_ok()
    );
    assert!(
        RateState {
            speed_range_mps: [1e-9, 1e-9],
            ..base
        }
        .validate()
        .is_err()
    );
    let theta = base.dc_m / base.v0_mps;
    assert!(base.coefficient(base.speed_range_mps[1], theta).is_ok());
    assert!(
        base.coefficient(base.speed_range_mps[1].next_up(), theta)
            .is_err()
    );
    assert!(
        base.coefficient(0.0, theta).is_err(),
        "zero speed is outside the log law"
    );
    let l = lugre();
    assert!(
        LuGre {
            static_n: l.coulomb_n,
            ..l
        }
        .validate()
        .is_ok()
    );
    assert!(
        LuGre {
            static_n: l.coulomb_n.next_down(),
            ..l
        }
        .validate()
        .is_err()
    );
    assert!(
        LuGre {
            sigma1_ns_per_m: 0.0,
            ..l
        }
        .validate()
        .is_ok()
    );
    assert!(
        LuGre {
            sigma1_ns_per_m: 0.0_f64.next_down(),
            ..l
        }
        .validate()
        .is_err()
    );
    let mut x = step(1e-9, 1e-3, 0.0);
    x.solve = solve(128);
    assert!(l.advance(&LuGreState::new(0.0).unwrap(), &x).is_ok());
    x.solve = solve(129);
    assert!(l.advance(&LuGreState::new(0.0).unwrap(), &x).is_err());
    let mut planar = step(1e-9, 1e-3, 0.0);
    planar.slip_increment_m[1] = 1e-9;
    assert!(
        l.advance(&l.initial_state(), &planar).is_err(),
        "LuGre is one-dimensional"
    );
}

// Smith & Woodhouse (2000), J. Mech. Phys. Solids 48:1633-1681.
/// Rosin density [kg/m^3], Table 1.
const SW_ROSIN_DENSITY: f64 = 1010.0;
/// Rosin specific heat [J/kg/K], Table 1 (from Cobbold & Jackson 1992; the
/// paper notes it may be low by a factor of 3 or more, section 5.3).
const SW_ROSIN_SPECIFIC_HEAT: f64 = 250.0;
/// Perspex, Table 1: the rod and the wedge of the steady-sliding and
/// sinusoidal rig (sections 3.1, 5.1).
const SW_PERSPEX: ThermalSolid = ThermalSolid {
    conductivity_w_per_m_k: 0.2,
    density_kg_per_m3: 1200.0,
    specific_heat_j_per_kg_k: 1460.0,
};
/// Eq. (2), fitted to the steady-sliding data of Fig. 4.
const SW_FIT: SteadySlidingFit = SteadySlidingFit {
    amplitudes: [0.4, 0.45],
    decay_speeds_mps: [0.01, 0.1],
    asymptote: 0.35,
};
/// Fig. 4 spans 0 to 0.3 m/s.
const SW_FIT_MAX_SPEED: f64 = 0.3;
/// Normal force of the steady-sliding and sinusoidal tests, section 5.1.
const SW_RIG_NORMAL_N: f64 = 3.0;

/// Eq. (2) written out independently of the law.
fn sw_eq2(v: f64) -> f64 {
    0.4 * det::exp(-v / 0.01) + 0.45 * det::exp(-v / 0.1) + 0.35
}

fn sw_law(radius_m: f64, layer_m: f64, step_s: f64, history_capacity: usize) -> SmithWoodhouse {
    let contact = RosinContact {
        radius_m,
        reference_normal_n: SW_RIG_NORMAL_N,
        layer_thickness_m: layer_m,
        stationary: SW_PERSPEX,
        moving: SW_PERSPEX,
    };
    SmithWoodhouse::new(SmithWoodhouseParams {
        rosin_density_kg_per_m3: SW_ROSIN_DENSITY,
        rosin_specific_heat_j_per_kg_k: SW_ROSIN_SPECIFIC_HEAT,
        steady_fit: SW_FIT,
        calibration: contact,
        calibration_speed_max_mps: SW_FIT_MAX_SPEED,
        contact,
        step_s,
        history_capacity,
    })
    .unwrap()
}

fn sw_step(slip_m: f64, dt_s: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: SW_RIG_NORMAL_N,
        slip_increment_m: [slip_m, 0.0],
        dt_s,
        contact_temperature_k: None,
        nominal_pressure_pa: None,
        compliance_m_per_n: 0.0,
        solve: solve(128),
    }
}

/// Prescribed velocity history sampled at step midpoints: per-step signed
/// friction coefficient and temperature rise.
fn sw_drive(
    law: &SmithWoodhouse,
    steps: usize,
    velocity: impl Fn(f64) -> f64,
) -> (Vec<f64>, Vec<f64>) {
    let h = law.params().step_s;
    let (mut mu, mut temperature) = (Vec::with_capacity(steps), Vec::with_capacity(steps));
    let mut s = law.initial_state();
    for k in 0..steps {
        #[allow(clippy::cast_precision_loss)]
        let v = velocity((k as f64 + 0.5) * h);
        let next = law.advance(&s, &sw_step(v * h, h)).unwrap();
        mu.push(next.traction_n[0] / SW_RIG_NORMAL_N);
        temperature.push(next.state.temperature_rise_k());
        s = next.state;
    }
    (mu, temperature)
}

#[test]
fn g1_smith_woodhouse_discrete_green_function_matches_eq13() {
    // Eq. (13) divided by 2 lambda h^{-1/2}: w_0 = 1, w_n = sqrt(n+1) -
    // 2 sqrt(n) + sqrt(n-1). Three independent checks: the literal form
    // where its cancellation is harmless, the telescoped partial sums
    // sum_{n<=N} w_n = sqrt(N+1) - sqrt(N), and Eq. (15).
    assert_eq!(half_space_unit_weight(0).to_bits(), 1.0_f64.to_bits());
    let mut worst_literal = 0.0_f64;
    for n in 1..=20_u32 {
        let x = f64::from(n);
        let literal = det::sqrt(x + 1.0) - 2.0 * det::sqrt(x) + det::sqrt(x - 1.0);
        let e = (half_space_unit_weight(n as usize) - literal).abs() / literal.abs();
        worst_literal = worst_literal.max(e);
    }
    let mut worst_sum = 0.0_f64;
    for big_n in [1_usize, 10, 100, 1_000, 10_000, 100_000] {
        // Smallest terms first; the rounding bound is eps * sum |w| / S_N.
        let tail: f64 = (1..=big_n).rev().map(half_space_unit_weight).sum();
        #[allow(clippy::cast_precision_loss)]
        let x = big_n as f64;
        let exact = 1.0 / (det::sqrt(x + 1.0) + det::sqrt(x));
        let e = ((1.0 + tail) - exact).abs() / exact;
        worst_sum = worst_sum.max(e);
        println!(
            "{{\"law\":\"benchmark-smith-woodhouse-plastic-thermal\",\"test\":\"eq13_partial_sum\",\"N\":{big_n},\"sum\":{},\"oracle\":{exact},\"rel_err\":{e:e},\"tol\":1e-12}}",
            1.0 + tail
        );
    }
    let mut worst_asymptote = 0.0_f64;
    for n in [10_usize, 100, 1_000, 100_000] {
        #[allow(clippy::cast_precision_loss)]
        let x = n as f64;
        let ratio = half_space_unit_weight(n) / (-0.25 / (x * det::sqrt(x)));
        // w_n = -n^{-3/2}/4 (1 + 5/(16 n^2) + ...).
        worst_asymptote = worst_asymptote.max((ratio - 1.0).abs() * x * x);
    }
    println!(
        "{{\"law\":\"benchmark-smith-woodhouse-plastic-thermal\",\"test\":\"eq13_weights\",\"literal_rel_err\":{worst_literal:e},\"partial_sum_rel_err\":{worst_sum:e},\"eq15_scaled_err\":{worst_asymptote},\"tol\":[1e-12,1e-12,1.0]}}"
    );
    assert!(worst_literal < 1e-12, "Eq. (13) literal {worst_literal:e}");
    assert!(worst_sum < 1e-12, "telescoped sums {worst_sum:e}");
    assert!(
        worst_asymptote < 1.0,
        "Eq. (15) n^2-scaled {worst_asymptote}"
    );
}

#[test]
fn g1_smith_woodhouse_calibration_reproduces_eq2_through_eq20() {
    // Section 5.3: Eq. (20) inverted along the fit of Eq. (2) defines
    // k_y(T). Along the calibration curve the law must return T from an
    // independent evaluation of Eq. (20), and A k_y(T)/N must return Eq. (2).
    let mut worst_t = 0.0_f64;
    let mut worst_mu = 0.0_f64;
    for (radius, layer) in [
        (0.5e-3, 10e-6),
        (0.25e-3, 10e-6),
        (0.5e-3, 1e-6),
        (0.2e-3, 30e-6),
    ] {
        let law = sw_law(radius, layer, 1e-4, 1);
        let area = core::f64::consts::PI * radius * radius;
        let d = SW_PERSPEX.conductivity_w_per_m_k
            / (SW_PERSPEX.density_kg_per_m3 * SW_PERSPEX.specific_heat_j_per_kg_k);
        for i in 0..=40 {
            let v = if i == 40 {
                SW_FIT_MAX_SPEED
            } else {
                1e-4 * det::exp(f64::from(i) / 40.0 * det::ln(SW_FIT_MAX_SPEED / 1e-4))
            };
            // Eq. (20): mu N v = v a delta rho_r c_r T + T A K_2 sqrt(3v/(8 a D_2)).
            let oracle_t = sw_eq2(v) * SW_RIG_NORMAL_N * v
                / (v * radius * layer * SW_ROSIN_DENSITY * SW_ROSIN_SPECIFIC_HEAT
                    + area
                        * SW_PERSPEX.conductivity_w_per_m_k
                        * det::sqrt(3.0 * v / (8.0 * radius * d)));
            let t = law.calibration_steady_temperature_k(v).unwrap();
            let mu = area * law.yield_stress_pa(t).unwrap() / SW_RIG_NORMAL_N;
            worst_t = worst_t.max((t - oracle_t).abs() / oracle_t);
            worst_mu = worst_mu.max((mu - sw_eq2(v)).abs() / sw_eq2(v));
        }
        let static_mu = area * law.yield_stress_pa(0.0).unwrap() / SW_RIG_NORMAL_N;
        assert!((static_mu - 1.2).abs() < 1e-15, "k_y(0) A/N {static_mu}");
        println!(
            "{{\"law\":\"benchmark-smith-woodhouse-plastic-thermal\",\"test\":\"calibration\",\"radius_m\":{radius},\"layer_m\":{layer},\"ky0_pa\":{},\"temperature_limit_k\":{},\"ky_at_limit_pa\":{}}}",
            law.yield_stress_pa(0.0).unwrap(),
            law.temperature_limit_k(),
            law.yield_stress_pa(law.temperature_limit_k()).unwrap()
        );
    }
    println!(
        "{{\"law\":\"benchmark-smith-woodhouse-plastic-thermal\",\"test\":\"calibration\",\"oracle\":\"eq20,eq2\",\"t_rel_err\":{worst_t:e},\"mu_rel_err\":{worst_mu:e},\"tol\":1e-12}}"
    );
    assert!(worst_t < 1e-12, "Eq. (20) temperature {worst_t:e}");
    assert!(worst_mu < 1e-12, "Eq. (2) through k_y(T) {worst_mu:e}");
}

#[test]
fn g1_smith_woodhouse_steady_sliding_approaches_eq2() {
    // Steady sliding on the section 5.1 rig (a = 0.5 mm, delta = 10 um,
    // N = 3 N, section 5.3 values). The moving-surface flux becomes steady
    // within a transit time, but the stationary half-space keeps drawing
    // K T/sqrt(pi D t), so mu(t) = mu_eq2(v) + c t^{-1/2} + O(t^{-1}).
    // Error model: the next term is O(sqrt(t_c/t)) relative to the first,
    // with t_c = 8a/(3 pi v). Oracles: the deviation halves over each 4x in
    // time to within that size, and one Richardson step in t^{-1/2}
    // removes at least 90 % of the deviation left at 4 s.
    let radius = 0.5e-3;
    for v in [0.05, 0.1, 0.2] {
        let transit = 8.0 * radius / (3.0 * core::f64::consts::PI * v);
        let h = transit / 16.0;
        let marks = [0.25, 1.0, 4.0];
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let ends: Vec<usize> = marks.iter().map(|t| (t / h).round() as usize).collect();
        let law = sw_law(radius, 10e-6, h, ends[2]);
        let (mu, _) = sw_drive(&law, ends[2], |_| v);
        let oracle = sw_eq2(v);
        #[allow(clippy::cast_precision_loss)]
        let times: Vec<f64> = ends.iter().map(|&n| n as f64 * h).collect();
        let dev: Vec<f64> = ends.iter().map(|&n| mu[n - 1] - oracle).collect();
        let (r1, r2) = (det::sqrt(times[1]), det::sqrt(times[2]));
        let extrapolated = (r2 * mu[ends[2] - 1] - r1 * mu[ends[1] - 1]) / (r2 - r1);
        let ratio = dev[2] / dev[1];
        let next_order = det::sqrt(transit / times[1]);
        let residual = (extrapolated - oracle).abs();
        println!(
            "{{\"law\":\"benchmark-smith-woodhouse-plastic-thermal\",\"test\":\"steady_sliding\",\"v_mps\":{v},\"h_s\":{h:e},\"t_s\":{times:?},\"mu\":[{},{},{}],\"oracle_eq2\":{oracle},\"deviation\":{dev:?},\"ratio_4s_1s\":{ratio},\"ratio_band\":{next_order},\"richardson_mu\":{extrapolated},\"richardson_residual\":{residual:e},\"residual_bound\":{}}}",
            mu[ends[0] - 1],
            mu[ends[1] - 1],
            mu[ends[2] - 1],
            0.1 * dev[2].abs()
        );
        assert!(
            dev[0] > dev[1] && dev[1] > dev[2] && dev[2] > 0.0,
            "approach from above: {dev:?}"
        );
        assert!((ratio - 0.5).abs() <= next_order, "t^-1/2 law: {ratio}");
        assert!(
            residual <= 0.1 * dev[2].abs(),
            "Richardson limit {extrapolated} vs Eq. (2) {oracle}"
        );
    }
}

#[test]
fn g1_smith_woodhouse_fig18_sinusoidal_hysteresis() {
    // Section 5.3, Fig. 18(b): forced sinusoidal sliding at 4 Hz, peak
    // 0.08 m/s, N = 3 N, a = 0.5 mm, delta = 10 um, both surfaces perspex.
    // Claims tested (sections 3.2, 5.3): the loops run clockwise, with
    // higher friction just after sticking than when the speed falls back
    // towards zero; they are confined to low speeds; and the two loops
    // close at the vertical line, because the friction magnitude is the same
    // just before and just after a reversal. Measured on the eighth
    // half-cycle (0.875-1 s) on a three-rung step ladder. The scheme is
    // first order a priori.
    let (peak, omega) = (0.08, 2.0 * core::f64::consts::PI * 4.0);
    let levels = [0.1, 0.25, 0.5, 0.75, 0.9];
    let rungs = [2e-4_f64, 1e-4, 5e-5];
    let mut gaps = Vec::new();
    let mut jumps = Vec::new();
    let mut finest_temperature = Vec::new();
    let mut middle_mu = Vec::new();
    for (rung, h) in rungs.into_iter().enumerate() {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let steps = (1.0 / h).round() as usize;
        let law = sw_law(0.5e-3, 10e-6, h, steps);
        let (mu, temperature) = sw_drive(&law, steps, |t| peak * det::sin(omega * t));
        let magnitude = |t: f64| -> f64 {
            // Linear interpolation of |mu| between step midpoints.
            let x = t / h - 0.5;
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let k = x.floor() as usize;
            #[allow(clippy::cast_precision_loss)]
            let f = x - k as f64;
            mu[k].abs() * (1.0 - f) + mu[k + 1].abs() * f
        };
        let gap: Vec<f64> = levels
            .iter()
            .map(|&s| {
                let rise = (det::asin(s) + 7.0 * core::f64::consts::PI) / omega;
                let fall = (8.0 * core::f64::consts::PI - det::asin(s)) / omega;
                magnitude(rise) - magnitude(fall)
            })
            .collect();
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let reversal = (7.0 * core::f64::consts::PI / omega / h).round() as usize;
        let jump = mu[reversal].abs() - mu[reversal - 1].abs();
        println!(
            "{{\"law\":\"benchmark-smith-woodhouse-plastic-thermal\",\"test\":\"fig18_loop\",\"h_s\":{h:e},\"speed_fractions\":{levels:?},\"rise_minus_fall\":{gap:?},\"reversal_jump\":{jump:e}}}"
        );
        gaps.push(gap);
        jumps.push(jump);
        if rung == 1 {
            middle_mu.clone_from(&mu);
        }
        finest_temperature = temperature;
    }
    let order = |e: [f64; 3]| det::ln((e[0] - e[1]).abs() / (e[1] - e[2]).abs()) / det::ln(2.0);
    let gap_order = order([gaps[0][1], gaps[1][1], gaps[2][1]]);
    let gap_error = (gaps[1][1] - gaps[2][1]).abs();
    let jump_order = det::ln(jumps[1].abs() / jumps[2].abs()) / det::ln(2.0);
    let finest = &gaps[2];
    let smallest = finest.iter().copied().fold(f64::INFINITY, f64::min);
    println!(
        "{{\"law\":\"benchmark-smith-woodhouse-plastic-thermal\",\"test\":\"fig18_loop\",\"gap_order\":{gap_order},\"gap_error_estimate\":{gap_error:e},\"smallest_gap\":{smallest},\"jump_order\":{jump_order},\"order_band\":[0.8,1.2]}}"
    );
    assert!(
        finest.iter().all(|&g| g > 0.0),
        "clockwise loops: rise minus fall {finest:?}"
    );
    assert!(
        smallest > 10.0 * gap_error,
        "loop sense resolved: {smallest} vs {gap_error}"
    );
    assert!(finest[4] < finest[1], "confined to low speed: {finest:?}");
    assert!(
        (0.8..=1.2).contains(&gap_order),
        "loop ladder order {gap_order}"
    );
    assert!(
        jumps[0].abs() > jumps[1].abs() && jumps[1].abs() > jumps[2].abs(),
        "loops close at the vertical line: {jumps:?}"
    );
    assert!(
        (0.8..=1.2).contains(&jump_order),
        "closure order {jump_order}"
    );

    sw_fig18_temperature_disclosure(peak, omega, rungs[2], &finest_temperature, &middle_mu);
}

/// Disclosure, not a gate: Fig. 18 temperature levels against the figure
/// and text, and the a = 0.25 mm run of Fig. 18(d).
fn sw_fig18_temperature_disclosure(
    peak: f64,
    omega: f64,
    h: f64,
    finest_temperature: &[f64],
    middle_mu: &[f64],
) {
    let half = |j: usize| {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let (a, b) = (
            (j as f64 * 0.125 / h).round() as usize,
            ((j + 1) as f64 * 0.125 / h).round() as usize,
        );
        let slice = &finest_temperature[a..b];
        (
            slice.iter().copied().fold(f64::INFINITY, f64::min),
            slice.iter().copied().fold(0.0, f64::max),
        )
    };
    let ranges: Vec<(f64, f64)> = (4..8).map(half).collect();
    println!(
        "{{\"law\":\"benchmark-smith-woodhouse-plastic-thermal\",\"test\":\"fig18_temperature\",\"radius_m\":0.5e-3,\"half_cycles_4_to_7_min_max_k\":{ranges:?},\"paper\":\"Fig. 18(a) reads about 27.5-40.5 K; text: below 40 K\",\"gated\":false}}"
    );
    let wide = sw_law(0.25e-3, 10e-6, 1e-4, 20_000);
    let (wide_mu, wide_temperature) = sw_drive(&wide, 20_000, |t| peak * det::sin(omega * t));
    let tail = &wide_temperature[10_000..];
    // Same step as the middle rung, so the two runs share sample times.
    let similar = middle_mu
        .iter()
        .zip(&wide_mu)
        .skip(5_000)
        .map(|(a, b)| (a.abs() - b.abs()).abs())
        .fold(0.0, f64::max);
    println!(
        "{{\"law\":\"benchmark-smith-woodhouse-plastic-thermal\",\"test\":\"fig18_temperature\",\"radius_m\":0.25e-3,\"t_1_to_2_s_min_max_k\":[{},{}],\"paper\":\"text: between 75 and 115 K\",\"max_abs_mu_difference_vs_0.5mm_0.5_to_1_s\":{similar},\"gated\":false}}",
        tail.iter().copied().fold(f64::INFINITY, f64::min),
        tail.iter().copied().fold(0.0, f64::max)
    );
}

#[test]
fn g0_smith_woodhouse_cap_and_cap_plus_one() {
    let base = *sw_law(0.5e-3, 10e-6, 1e-4, 4).params();
    let ok = |p: SmithWoodhouseParams| SmithWoodhouse::new(p).is_ok();
    // Fit amplitudes are non-negative; the asymptote is positive.
    let mut p = base;
    p.steady_fit.amplitudes[0] = 0.0;
    assert!(ok(p));
    p.steady_fit.amplitudes[0] = 0.0_f64.next_down();
    assert!(!ok(p));
    p = base;
    p.steady_fit.asymptote = 0.0;
    assert!(!ok(p));
    p = base;
    p.contact.layer_thickness_m = 0.0;
    assert!(!ok(p));
    p = base;
    p.calibration.moving.conductivity_w_per_m_k = 0.0;
    assert!(!ok(p));
    p = base;
    p.step_s = 0.0;
    assert!(!ok(p));
    // History capacity: 1..=MAX_HISTORY_STEPS.
    p = base;
    p.history_capacity = 1;
    assert!(ok(p));
    p.history_capacity = 0;
    assert!(!ok(p));
    p.history_capacity = MAX_HISTORY_STEPS;
    assert!(ok(p));
    p.history_capacity = MAX_HISTORY_STEPS + 1;
    assert!(!ok(p));
    // A calibration whose T_ss(v) is not certified increasing refuses.
    p = base;
    p.steady_fit = SteadySlidingFit {
        amplitudes: [10.0, 0.0],
        decay_speeds_mps: [0.01, 0.1],
        asymptote: 0.01,
    };
    assert!(matches!(
        SmithWoodhouse::new(p),
        Err(TriboError::InvalidInput { .. })
    ));

    let law = SmithWoodhouse::new(base).unwrap();
    // Calibration speed and temperature ranges.
    assert!(
        law.calibration_steady_temperature_k(SW_FIT_MAX_SPEED)
            .is_ok()
    );
    assert!(
        law.calibration_steady_temperature_k(SW_FIT_MAX_SPEED.next_up())
            .is_err()
    );
    let limit = law.temperature_limit_k();
    assert!(law.yield_stress_pa(limit).is_ok());
    assert!(matches!(
        law.yield_stress_pa(limit.next_up()),
        Err(TriboError::OutsideApplicability { .. })
    ));
    assert!(law.yield_stress_pa(0.0).is_ok());
    assert!(law.yield_stress_pa(0.0_f64.next_down()).is_err());

    // Step admission: the declared step only, a loaded contact, and the
    // history capacity at cap and cap + 1.
    let h = base.step_s;
    let s0 = law.initial_state();
    assert!(law.advance(&s0, &sw_step(1e-6, h)).is_ok());
    assert!(law.advance(&s0, &sw_step(1e-6, h.next_up())).is_err());
    let mut unloaded = sw_step(1e-6, h);
    unloaded.normal_n = 0.0;
    assert!(law.advance(&s0, &unloaded).is_err());
    let mut s = s0.clone();
    for _ in 0..base.history_capacity {
        s = law.advance(&s, &sw_step(1e-6, h)).unwrap().state;
    }
    assert!(matches!(
        law.advance(&s, &sw_step(1e-6, h)),
        Err(TriboError::OutsideApplicability { .. })
    ));

    // Force budget at the sliding limit A k_y(0) and one ulp below it.
    let sliding = law.advance(&s0, &sw_step(1e-6, h)).unwrap();
    let force = sliding.traction_n[0];
    assert!(
        (force / SW_RIG_NORMAL_N - 1.2).abs() < 1e-14,
        "mu_s {force}"
    );
    let mut x = sw_step(1e-6, h);
    x.solve.maximum_force_n = force;
    assert!(law.advance(&s0, &x).is_ok());
    x.solve.maximum_force_n = force.next_down();
    assert!(law.advance(&s0, &x).is_err());

    // Set-valued sgn(v): a compliant step sticks at or below the limit,
    // with no heat, and slides at the limit above it.
    let compliance = 1e-6;
    let mut stick = sw_step(0.5 * force * compliance, h);
    stick.compliance_m_per_n = compliance;
    let held = law.advance(&s0, &stick).unwrap();
    assert_eq!(held.regime, FrictionRegime::Sticking);
    assert_eq!(held.slip_increment_m.map(f64::to_bits), [0; 2]);
    assert_eq!(held.state.temperature_rise_k().to_bits(), 0.0_f64.to_bits());
    let mut slide = sw_step(2.0 * force * compliance, h);
    slide.compliance_m_per_n = compliance;
    let slid = law.advance(&s0, &slide).unwrap();
    assert_eq!(slid.regime, FrictionRegime::Sliding);
    assert_eq!(slid.traction_n[0].to_bits(), force.to_bits());
    assert!(slid.state.temperature_rise_k() > 0.0);
    assert!(
        (slid.dissipated_j - force * slid.slip_increment_m[0]).abs() <= 1e-15 * slid.dissipated_j
    );
}

fn files_under(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            files_under(&path, out);
        } else {
            out.push(path);
        }
    }
}

#[test]
fn g0_benchmark_laws_not_reachable_from_render_paths() {
    // No crate may enable the feature in [dependencies], and no library
    // source outside fs-tribo may name the module.
    let crates = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_owned();
    let mut checked = 0;
    for entry in fs::read_dir(&crates).unwrap().flatten() {
        let manifest = entry.path().join("Cargo.toml");
        let Ok(text) = fs::read_to_string(&manifest) else {
            continue;
        };
        let mut section = String::new();
        for line in text.lines() {
            let line = line.trim();
            if line.starts_with('[') {
                section = line.to_owned();
            } else if line.contains("benchmark-laws") && !line.starts_with('#') {
                let own = entry.file_name() == "fs-tribo";
                assert!(
                    section == "[dev-dependencies]"
                        || (own && (section == "[features]" || section == "[[test]]")),
                    "{} enables benchmark-laws under {section}",
                    manifest.display()
                );
            }
        }
        if entry.file_name() != "fs-tribo" {
            let mut sources = Vec::new();
            files_under(&entry.path().join("src"), &mut sources);
            for source in sources
                .iter()
                .filter(|p| p.extension().is_some_and(|e| e == "rs"))
            {
                let text = fs::read_to_string(source).unwrap_or_default();
                assert!(
                    !text.contains("fs_tribo::benchmark"),
                    "{} names fs_tribo::benchmark",
                    source.display()
                );
            }
        }
        checked += 1;
    }
    println!("{{\"test\":\"benchmark_laws_not_reachable\",\"crates_checked\":{checked}}}");
    assert!(checked > 100, "the scan must see the workspace ({checked})");
}
