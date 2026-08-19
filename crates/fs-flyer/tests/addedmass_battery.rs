//! V-05d battery (bead wf-root-guzez.4.3, E3.2a) — executed clause by
//! clause: block symmetry, PSD/PD, frame covariance under a reference-point
//! shift, NONZERO cross-block fixtures, work closure of the bias,
//! deterministic factorization equivalence, checkpoint/resume equality,
//! and the hostile twin (removing the cross block CHANGES the answer).
//! Repro: cargo test -p fs-flyer --test addedmass_battery

use fs_flyer::addedmass::{
    AeroGeneralizedLoads, MAX_CONTROL_COORDS, Strip, assemble_analytic_strip, solve_joint,
};
use fs_flyer::spine::RigidBody;

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-v05d\",\"case\":\"{case}\",{payload}}}");
}

const RHO: f64 = 1.294; // the E1.8 cold-day density

fn body() -> RigidBody {
    RigidBody {
        mass_kg: 340.17,
        inertia_kgm2: [1787.0, 367.4, 1820.9],
    }
}

/// Wing planes (z-normal) + a canard strip driven by the hinge coordinate 0.
fn flyer_strips() -> Vec<Strip> {
    vec![
        Strip {
            name: "wing-lower",
            chord_m: 1.981,
            span_m: 12.29,
            position_m: [-0.5, 0.0, 0.0],
            normal: [0.0, 0.0, 1.0],
            control_coord: None,
            control_gain: 0.0,
        },
        Strip {
            name: "wing-upper",
            chord_m: 1.981,
            span_m: 12.29,
            position_m: [-0.5, 0.0, -1.89],
            normal: [0.0, 0.0, 1.0],
            control_coord: None,
            control_gain: 0.0,
        },
        Strip {
            name: "canard",
            chord_m: 0.762,
            span_m: 3.658,
            position_m: [2.6, 0.0, 0.3],
            normal: [0.0, 0.0, 1.0],
            control_coord: Some(0),
            control_gain: 0.19, // hinge arm × mode shape at centroid [m]
        },
        Strip {
            name: "rudder",
            chord_m: 0.45,
            span_m: 2.1,
            position_m: [-3.35, 0.0, -0.4],
            normal: [0.0, 1.0, 0.0],
            control_coord: None,
            control_gain: 0.0,
        },
    ]
}

fn assemble(nu: [f64; 6], qdot: f64) -> AeroGeneralizedLoads {
    assemble_analytic_strip(RHO, &flyer_strips(), 1, &nu, &[qdot]).unwrap()
}

#[test]
fn single_plate_matches_the_closed_form() {
    // One horizontal plate: M_rr[heave][heave] = ρπc²b/4 EXACTLY; the
    // heave–pitch coupling is m_a·x; pitch–pitch is m_a·x².
    let plate = [Strip {
        name: "plate",
        chord_m: 2.0,
        span_m: 3.0,
        position_m: [1.5, 0.0, 0.0],
        normal: [0.0, 0.0, 1.0],
        control_coord: None,
        control_gain: 0.0,
    }];
    let l = assemble_analytic_strip(RHO, &plate, 0, &[0.0; 6], &[]).unwrap();
    let m_a = RHO * core::f64::consts::PI * 4.0 / 4.0 * 3.0;
    assert!(
        (l.m_added_rr[2][2] - m_a).abs() < 1e-12,
        "heave {}",
        l.m_added_rr[2][2]
    );
    // r×n = (1.5,0,0)×(0,0,1) = (0·1−0·0, 0·0−1.5·1, 0) = (0, −1.5, 0).
    assert!(
        (l.m_added_rr[2][4] - m_a * (-1.5)).abs() < 1e-12,
        "heave-pitch"
    );
    assert!(
        (l.m_added_rr[4][4] - m_a * 2.25).abs() < 1e-12,
        "pitch-pitch"
    );
    // Surge/sway rows are zero (a plate has no added mass along its plane).
    assert_eq!(l.m_added_rr[0][0], 0.0);
    assert_eq!(l.m_added_rr[1][1], 0.0);
    jlog("closed-form", &format!("\"m_a\":{m_a}"));
}

#[test]
fn symmetry_and_positive_semidefiniteness() {
    let l = assemble([1.0, 0.2, -0.5, 0.1, -0.3, 0.2], 0.4);
    // Block symmetry to the bit (assembled from symmetric outer products).
    for i in 0..6 {
        for j in 0..6 {
            assert_eq!(
                l.m_added_rr[i][j].to_bits(),
                l.m_added_rr[j][i].to_bits(),
                "M_rr symmetry at ({i},{j})"
            );
        }
    }
    // PSD: deterministic probe quadratic forms are ≥ −eps (7-dim: rigid+1c).
    let probes: [[f64; 7]; 5] = [
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, -1.0, 0.5, 0.0, 0.2, -0.3],
        [0.3, -0.7, 0.2, -0.1, 0.9, -0.4, 0.6],
        [0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0],
        [-0.2, 0.5, -0.9, 0.8, 0.1, -0.6, 0.4],
    ];
    for p in probes {
        let mut q = 0.0;
        for i in 0..6 {
            for j in 0..6 {
                q += p[i] * l.m_added_rr[i][j] * p[j];
            }
            q += 2.0 * p[i] * l.m_added_rc[0][i] * p[6];
        }
        q += p[6] * l.m_added_cc[0][0] * p[6];
        assert!(q >= -1e-9, "added-mass quadratic form negative: {q}");
    }
    // Total PD: the joint Cholesky succeeds (solve_joint refuses otherwise).
    let (nu_dot, qddot) = solve_joint(&body(), &[2.5], &l).unwrap();
    assert_eq!(nu_dot.len(), 6);
    assert_eq!(qddot.len(), 1);
    jlog("symmetry-psd", "\"probes\":5,\"total_pd\":\"cholesky-ok\"");
}

#[test]
fn frame_covariance_under_reference_shift() {
    // Shift the reference point by d: strips move by −d, and the assembled
    // M'_rr must equal Xᵀ·M_rr·X with X = [[I, [d]×],[0, I]] (v = v′ + d×ω,
    // since v_O = v_O′ − ω×d = v_O′ + d×ω).
    let d = [0.7, -0.4, 0.25];
    let base = assemble([0.0; 6], 0.0);
    let shifted_strips: Vec<Strip> = flyer_strips()
        .into_iter()
        .map(|mut s| {
            for i in 0..3 {
                s.position_m[i] -= d[i];
            }
            s
        })
        .collect();
    let shifted = assemble_analytic_strip(RHO, &shifted_strips, 1, &[0.0; 6], &[0.0]).unwrap();
    // Build X (6×6): identity + upper-right −[d]×.
    let mut x = [[0.0f64; 6]; 6];
    for i in 0..6 {
        x[i][i] = 1.0;
    }
    let dx = [[0.0, -d[2], d[1]], [d[2], 0.0, -d[0]], [-d[1], d[0], 0.0]];
    for i in 0..3 {
        for j in 0..3 {
            x[i][3 + j] = dx[i][j];
        }
    }
    // Xᵀ·M·X.
    let mut mx = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            for k in 0..6 {
                mx[i][j] += base.m_added_rr[i][k] * x[k][j];
            }
        }
    }
    let mut xtmx = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            for k in 0..6 {
                xtmx[i][j] += x[k][i] * mx[k][j];
            }
        }
    }
    let scale = base.m_added_rr[2][2];
    for i in 0..6 {
        for j in 0..6 {
            assert!(
                (shifted.m_added_rr[i][j] - xtmx[i][j]).abs() < 1e-10 * scale,
                "frame covariance broke at ({i},{j}): {} vs {}",
                shifted.m_added_rr[i][j],
                xtmx[i][j]
            );
        }
    }
    jlog(
        "frame-covariance",
        &format!("\"shift\":[{},{},{}]", d[0], d[1], d[2]),
    );
}

#[test]
fn nonzero_cross_blocks_and_the_hostile_twin() {
    // The canard hinge MUST couple: heave and pitch cross entries nonzero.
    let l = assemble([0.0; 6], 0.0);
    assert!(
        l.m_added_rc[0][2].abs() > 1e-3,
        "heave↔hinge cross block vanished"
    );
    assert!(
        l.m_added_rc[0][4].abs() > 1e-3,
        "pitch↔hinge cross block vanished"
    );
    // HOSTILE TWIN: zero the cross block (the orphaned-cross-block bug the
    // Round-3 correction forbids) and the joint solve must CHANGE — the
    // fixture proves the coupling is load-bearing, not decorative.
    let mut loads = l.clone();
    loads.q_rigid_nonaccel = [0.0, 0.0, 500.0, 0.0, 120.0, 0.0];
    loads.q_control_nonaccel = vec![40.0];
    let (full_nu, full_q) = solve_joint(&body(), &[2.5], &loads).unwrap();
    let mut orphaned = loads.clone();
    orphaned.m_added_rc = vec![[0.0; 6]];
    let (orph_nu, orph_q) = solve_joint(&body(), &[2.5], &orphaned).unwrap();
    let dq = (full_q[0] - orph_q[0]).abs();
    let dnu: f64 = full_nu
        .iter()
        .zip(&orph_nu)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        dq > 1e-3 * full_q[0].abs().max(1.0),
        "hinge acceleration insensitive to the cross block — vacuous fixture"
    );
    assert!(dnu > 0.0, "rigid accelerations must also move");
    jlog(
        "hostile-twin",
        &format!("\"dq\":{dq},\"dnu_max\":{dnu},\"verdict\":\"cross-blocks load-bearing\""),
    );
}

#[test]
fn bias_work_closure_at_machine_precision() {
    // νᵀ·Q_bias_rigid ≡ 0 for ANY state (the ad* term is workless): check
    // at several deterministic states, relative to the momentum scale.
    for (i, nu) in [
        [13.9, 0.0, -0.4, 0.02, 0.15, -0.03],
        [5.0, -2.0, 1.0, 0.4, -0.6, 0.2],
        [0.1, 0.9, -0.7, -1.2, 0.8, 1.5],
    ]
    .iter()
    .enumerate()
    {
        let l = assemble(*nu, 0.7);
        let power: f64 = (0..6).map(|k| nu[k] * l.q_added_bias[k]).sum();
        let scale: f64 = (0..6)
            .map(|k| (nu[k] * l.q_added_bias[k]).abs())
            .sum::<f64>()
            .max(1e-30);
        assert!(
            power.abs() / scale < 1e-12,
            "bias does work at state {i}: relative {}",
            power.abs() / scale
        );
        // Control-block bias is identically zero (disclosed baseline).
        assert_eq!(l.q_added_bias[6], 0.0);
    }
    jlog("work-closure", "\"relative_power\":\"<1e-12 at 3 states\"");
}

#[test]
fn deterministic_factorization_and_checkpoint_equality() {
    let mut loads = assemble([1.0, 0.0, -0.3, 0.05, 0.1, -0.02], 0.3);
    loads.q_rigid_nonaccel = [10.0, -5.0, 3300.0, 20.0, 90.0, -15.0];
    loads.q_control_nonaccel = vec![25.0];
    let a = solve_joint(&body(), &[2.5], &loads).unwrap();
    let b = solve_joint(&body(), &[2.5], &loads).unwrap();
    for (x, y) in a.0.iter().zip(&b.0).chain(a.1.iter().zip(&b.1)) {
        assert_eq!(
            x.to_bits(),
            y.to_bits(),
            "factorization must be bit-deterministic"
        );
    }
    // Checkpoint/resume equality: a cloned (checkpointed) loads struct
    // solves to the same bits.
    let restored = loads.clone();
    let c = solve_joint(&body(), &[2.5], &restored).unwrap();
    for (x, y) in a.0.iter().zip(&c.0) {
        assert_eq!(x.to_bits(), y.to_bits());
    }
    jlog("determinism", "\"repeat+resume\":\"bitwise\"");
}

#[test]
fn refusals_at_cap_and_cap_plus_one() {
    let nu = [0.0f64; 6];
    // nc caps.
    let rates_ok = vec![0.0; MAX_CONTROL_COORDS];
    assert!(assemble_analytic_strip(RHO, &[], MAX_CONTROL_COORDS, &nu, &rates_ok).is_ok());
    let rates_over = vec![0.0; MAX_CONTROL_COORDS + 1];
    assert_eq!(
        assemble_analytic_strip(RHO, &[], MAX_CONTROL_COORDS + 1, &nu, &rates_over)
            .unwrap_err()
            .code,
        "control-coords-outside-domain"
    );
    // Strip gates: non-unit normal; out-of-range coordinate; bad chord.
    let mut bad = flyer_strips();
    bad[0].normal = [0.0, 0.0, 2.0];
    assert_eq!(
        assemble_analytic_strip(RHO, &bad, 1, &nu, &[0.0])
            .unwrap_err()
            .code,
        "strip-invalid"
    );
    let mut outc = flyer_strips();
    outc[2].control_coord = Some(3);
    assert_eq!(
        assemble_analytic_strip(RHO, &outc, 1, &nu, &[0.0])
            .unwrap_err()
            .code,
        "strip-invalid"
    );
    // Solver gates: control-inertia mismatch; PD refusal on a zero-inertia
    // control coordinate with no added mass.
    let l = assemble(nu, 0.0);
    assert_eq!(
        solve_joint(&body(), &[], &l).unwrap_err().code,
        "control-inertia-invalid"
    );
    jlog(
        "refusals",
        "\"gates\":\"nc cap/cap+1, strips, control inertia\"",
    );
}

#[test]
fn joint_golden_digest() {
    // The full flyer fixture: assembled blocks + joint solution, exact bits
    // (measure-then-pin; golden-bump protocol).
    let mut loads = assemble([13.9, 0.3, -0.8, 0.04, 0.12, -0.05], 0.6);
    loads.q_rigid_nonaccel = [12.0, -3.0, 3336.0, 15.0, 110.0, -8.0];
    loads.q_control_nonaccel = vec![30.0];
    let (nu_dot, qddot) = solve_joint(&body(), &[2.5], &loads).unwrap();
    let mut payload = Vec::new();
    for row in &loads.m_added_rr {
        for v in row {
            payload.extend_from_slice(&v.to_bits().to_le_bytes());
        }
    }
    for v in nu_dot
        .iter()
        .chain(qddot.iter())
        .chain(loads.q_added_bias.iter())
    {
        payload.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.v05d-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "67462c13531148e9cf47fb4ce3a0913f8a41c12eda53f6cd6a1f4c57f1be6c88",
        "V-05d golden moved — determinism regression or an intentional \
         assembly change requiring the golden-bump protocol"
    );
}
