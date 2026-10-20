//! M4.2-W4 benchmark tangential laws (IRP pexfp.150): classic oracles for
//! rate-and-state (Dieterich 1979; Ruina 1983) and LuGre (Canudas de Wit et
//! al. 1995). Gauntlet tiers are named per test; each prints JSON lines.
//! All parameters are synthetic fixture values, not material data.

use fs_math::det;
use fs_rand::StreamKey;
use fs_tribo::benchmark::{Evolution, LuGre, LuGreState, RateState, RateStateState};
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStepInput, TriboError};
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
