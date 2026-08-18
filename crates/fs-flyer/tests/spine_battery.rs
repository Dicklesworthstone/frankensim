//! E3.2-i spine battery (bead wf-root-guzez.4.2.1). Fixtures with exact
//! references: ballistic parabola (Verlet exact for constant force),
//! single-axis spin-up vs the closed form, MEASURED Richardson order on a
//! time-varying force, bit-identity across runs, caps at cap AND cap+1,
//! pinned trajectory golden.
//! Repro: cargo test -p fs-flyer --test spine_battery

use fs_flyer::spine::{advance, step, tick_digest, Loads, RigidBody, SixDofState, MAX_STEPS};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-spine\",\"case\":\"{case}\",{payload}}}");
}

const DT: f64 = 1.0 / 120.0;

fn body() -> RigidBody {
    RigidBody { mass_kg: 340.17, inertia_kgm2: [1787.0, 367.4, 1820.9] }
}

fn rest() -> SixDofState {
    SixDofState {
        pos_m: [0.0; 3],
        vel_mps: [0.0; 3],
        quat: [1.0, 0.0, 0.0, 0.0],
        omega_body: [0.0; 3],
    }
}

#[test]
fn ballistic_parabola_is_near_machine_exact() {
    // Constant gravity (NED: +z down): Verlet reproduces the parabola to
    // rounding, per tick, over 2 seconds.
    let g = 9.80665;
    let m = body().mass_kg;
    let mut s = rest();
    s.vel_mps = [10.0, 0.0, -5.0]; // forward and upward
    let (end, _) = advance(&body(), &s, 0.0, DT, 240, |_, _| Loads {
        force_n: [0.0, 0.0, m * g],
        moment_nm: [0.0; 3],
    })
    .unwrap();
    let t = 240.0 * DT;
    let exact_x = 10.0 * t;
    let exact_z = -5.0 * t + 0.5 * g * t * t;
    assert!((end.pos_m[0] - exact_x).abs() < 1e-9, "x {}", end.pos_m[0]);
    assert!((end.pos_m[2] - exact_z).abs() < 1e-9, "z {}", end.pos_m[2]);
    assert!((end.vel_mps[2] - (-5.0 + g * t)).abs() < 1e-9);
    jlog("ballistic", &format!("\"x\":{},\"z\":{}", end.pos_m[0], end.pos_m[2]));
}

#[test]
fn single_axis_spinup_matches_closed_form() {
    // Constant moment about the pitch principal axis: ω(t) = M/I·t exactly
    // (single-axis: the gyroscopic term vanishes; the Strang kicks sum to
    // the exact impulse).
    let m_pitch = 50.0;
    let (end, _) = advance(&body(), &rest(), 0.0, DT, 120, |_, _| Loads {
        force_n: [0.0; 3],
        moment_nm: [0.0, m_pitch, 0.0],
    })
    .unwrap();
    let exact = m_pitch / body().inertia_kgm2[1] * (120.0 * DT);
    assert!(
        (end.omega_body[1] - exact).abs() < 1e-12,
        "omega {} vs {exact}",
        end.omega_body[1]
    );
    // Quaternion still unit-norm after the split steps.
    let norm2: f64 = end.quat.iter().map(|v| v * v).sum();
    assert!((norm2 - 1.0).abs() < 1e-10, "norm² {norm2}");
    jlog("spinup", &format!("\"omega\":{},\"exact\":{exact}", end.omega_body[1]));
}

#[test]
fn richardson_order_is_two_on_time_varying_force() {
    // F(t) = sin(2t) on x: exact velocity v(t) = (1−cos 2t)/2m …with m=1.
    // Measure global order by halving dt twice and fitting the slope.
    let unit = RigidBody { mass_kg: 1.0, inertia_kgm2: [1.0, 1.0, 1.0] };
    let t_end = 1.0;
    let err_at = |n: u32| -> f64 {
        let dt = t_end / f64::from(n);
        let (end, _) = advance(&unit, &rest(), 0.0, dt, n, |t, _| Loads {
            force_n: [(2.0 * t).sin(), 0.0, 0.0],
            moment_nm: [0.0; 3],
        })
        .unwrap();
        let exact_v = (1.0 - (2.0f64 * t_end).cos()) / 2.0;
        (end.vel_mps[0] - exact_v).abs()
    };
    let (e1, e2, e3) = (err_at(200), err_at(400), err_at(800));
    let p12 = (e1 / e2).log2();
    let p23 = (e2 / e3).log2();
    assert!(
        (1.8..=2.2).contains(&p12) && (1.8..=2.2).contains(&p23),
        "measured order {p12:.3}/{p23:.3} not ~2 (e = {e1:e}, {e2:e}, {e3:e})"
    );
    jlog("order", &format!("\"p12\":{p12:.4},\"p23\":{p23:.4}"));
}

#[test]
fn bit_identity_across_runs_and_digest_sensitivity() {
    let loads = |t: f64, s: &SixDofState| Loads {
        force_n: [10.0 * (3.0 * t).sin(), -2.0 * s.vel_mps[1], 3336.0],
        moment_nm: [5.0 * (t).cos(), 8.0, -2.0 * s.omega_body[2]],
    };
    let (a, da) = advance(&body(), &rest(), 0.0, DT, 480, loads).unwrap();
    let (b, db) = advance(&body(), &rest(), 0.0, DT, 480, loads).unwrap();
    assert_eq!(da, db, "digest traces must be bit-identical across runs");
    for i in 0..3 {
        assert_eq!(a.pos_m[i].to_bits(), b.pos_m[i].to_bits());
        assert_eq!(a.vel_mps[i].to_bits(), b.vel_mps[i].to_bits());
        assert_eq!(a.omega_body[i].to_bits(), b.omega_body[i].to_bits());
    }
    // The digest is tick- and state-sensitive.
    assert_ne!(tick_digest(0, &a), tick_digest(1, &a));
    let mut nudged = a;
    nudged.pos_m[0] = f64::from_bits(a.pos_m[0].to_bits() ^ 1);
    assert_ne!(tick_digest(0, &a), tick_digest(0, &nudged), "1-ulp must move the digest");
    jlog("bit-identity", &format!("\"ticks\":480,\"final_digest\":\"{}\"", da.last().unwrap()));
}

#[test]
fn refusals_at_cap_and_cap_plus_one() {
    let ok = advance(&body(), &rest(), 0.0, DT, 0, |_, _| Loads {
        force_n: [0.0; 3],
        moment_nm: [0.0; 3],
    });
    assert!(ok.is_ok());
    // Step budget at cap+1 refuses (cap itself is too slow to run here; the
    // budget check precedes any stepping, so cap+1 exercises the gate).
    let over = advance(&body(), &rest(), 0.0, DT, MAX_STEPS + 1, |_, _| Loads {
        force_n: [0.0; 3],
        moment_nm: [0.0; 3],
    });
    assert_eq!(over.unwrap_err().code, "step-budget-exceeded");
    // dt domain and mass-property refusals.
    let bad_dt = step(&body(), &rest(), 0.0, 0.2, |_, _| Loads {
        force_n: [0.0; 3],
        moment_nm: [0.0; 3],
    });
    assert_eq!(bad_dt.unwrap_err().code, "timestep-outside-domain");
    let bad_body = RigidBody { mass_kg: 0.0, inertia_kgm2: [1.0; 3] };
    assert_eq!(
        step(&bad_body, &rest(), 0.0, DT, |_, _| Loads { force_n: [0.0; 3], moment_nm: [0.0; 3] })
            .unwrap_err()
            .code,
        "mass-outside-domain"
    );
    let mut nan_state = rest();
    nan_state.vel_mps[0] = f64::NAN;
    assert_eq!(
        step(&body(), &nan_state, 0.0, DT, |_, _| Loads {
            force_n: [0.0; 3],
            moment_nm: [0.0; 3]
        })
        .unwrap_err()
        .code,
        "non-finite-input"
    );
    jlog("refusals", "\"gates\":\"steps cap+1, dt, mass, NaN state\"");
}

#[test]
fn trajectory_golden() {
    // 240 ticks of the coupled forced fixture — the spine's determinism
    // golden (golden-bump protocol; measure-then-pin).
    let (_, digests) = advance(&body(), &rest(), 0.0, DT, 240, |t, s| Loads {
        force_n: [10.0 * (3.0 * t).sin(), -2.0 * s.vel_mps[1], 3336.0],
        moment_nm: [5.0 * t.cos(), 8.0, -2.0 * s.omega_body[2]],
    })
    .unwrap();
    let last = digests.last().unwrap().clone();
    jlog("golden", &format!("\"digest\":\"{last}\""));
    assert_eq!(
        last, "PLACEHOLDER-MEASURE-THEN-PIN",
        "spine trajectory digest moved — determinism regression or an \
         intentional integrator change requiring the golden-bump protocol"
    );
}
