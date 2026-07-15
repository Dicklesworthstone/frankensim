//! SE(3) lane battery (bead frankensim-ext-time-se3-lanes-3ol0).
//!
//! Detailed-logging discipline: every assertion window prints the
//! measured quantity it gates so a failure names the physics, not just
//! the boolean.

use fs_ga::{Motor, Point as GaPoint};
use fs_time::se3::{
    DepSolveParams, RenormPolicy, Se3ClaimClass, Twist, canonicalize_motor, dep_free_step,
    dep_momentum_adjoint, run_dep_free, se3_exp_step, se3_exp_step_renorm, se3_rigid_body_step,
};
use fs_time::{quat_exp_step, quat_rotate};

const INERTIA: [f64; 3] = [1.0, 2.0, 3.0];

#[test]
fn se3_001_pure_rotation_agrees_with_so3_quaternion_lane() {
    let omega = [0.3, -0.7, 0.5];
    let h = 1e-2;
    let steps = 200;
    let mut motor = Motor::identity();
    let mut q = [1.0, 0.0, 0.0, 0.0];
    for _ in 0..steps {
        motor = se3_exp_step(
            &motor,
            &Twist {
                omega,
                vel: [0.0; 3],
            },
            h,
        )
        .expect("exp step");
        q = quat_exp_step(q, omega, h);
    }
    let probe = [0.8, -0.4, 1.2];
    let via_motor = motor
        .transform_point(GaPoint {
            x: probe[0],
            y: probe[1],
            z: probe[2],
        })
        .expect("finite point");
    let via_quat = quat_rotate(q, probe);
    let err = ((via_motor.x - via_quat[0]).powi(2)
        + (via_motor.y - via_quat[1]).powi(2)
        + (via_motor.z - via_quat[2]).powi(2))
    .sqrt();
    println!("se3-001: SO(3)-lane agreement error {err:e}");
    assert!(err < 1e-11, "SE(3) and SO(3) lanes disagree: {err:e}");
}

#[test]
fn se3_002_constant_twist_steps_compose_to_one_exponential() {
    // One-parameter subgroup: N steps of h equal one step of N·h.
    let twist = Twist {
        omega: [0.4, 0.1, -0.2],
        vel: [0.5, -0.3, 0.8],
    };
    let h = 0.05;
    let steps = 40;
    let mut walked = Motor::identity();
    for _ in 0..steps {
        walked = se3_exp_step(&walked, &twist, h).expect("exp step");
    }
    let jumped = se3_exp_step(&Motor::identity(), &twist, h * f64::from(steps))
        .expect("single exponential");
    let probe = GaPoint {
        x: -0.6,
        y: 0.9,
        z: 0.3,
    };
    let a = walked.transform_point(probe).expect("finite");
    let b = jumped.transform_point(probe).expect("finite");
    let err = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt();
    println!("se3-002: screw one-parameter composition error {err:e}");
    assert!(err < 1e-10, "constant twist is not exact: {err:e}");
    let defect = {
        let mut m = walked;
        m.renormalize()
    };
    println!("se3-002: walked-motor renormalization drift {defect:e}");
}

#[test]
fn se3_003_double_cover_canonicalization_is_deterministic() {
    let m = Motor::rotor([0.0, 0.0, 1.0], 1.3).compose(&Motor::translator(0.2, -0.5, 0.7));
    let flipped = Motor(m.0.scale(-1.0));
    let (ca, fa) = canonicalize_motor(&m).expect("canonical");
    let (cb, fb) = canonicalize_motor(&flipped).expect("canonical");
    assert_ne!(fa, fb, "exactly one representative must flip");
    for (x, y) in ca.0.0.iter().zip(cb.0.0.iter()) {
        assert_eq!(x.to_bits(), y.to_bits(), "canonical representatives differ");
    }
    // Bitwise replay across two identical long runs that cross the
    // scalar-zero surface (rotation through pi).
    let twist = Twist {
        omega: [0.0, 0.0, 2.5],
        vel: [0.1, 0.0, 0.0],
    };
    let run = || -> Vec<u64> {
        let mut m = Motor::identity();
        let mut bits = Vec::new();
        for _ in 0..2000 {
            m = se3_exp_step(&m, &twist, 1e-3).expect("step");
            bits.push(m.0.0[0].to_bits());
        }
        bits
    };
    assert_eq!(run(), run(), "se3-003: replay is not bitwise deterministic");
}

#[test]
fn se3_004_dep_free_body_conserves_momentum_and_bounds_energy() {
    let omega0 = [0.7, 1.1, -0.4];
    let h = 1e-3;
    let steps = 10_000;
    let (_, _, receipt) = run_dep_free(
        [1.0, 0.0, 0.0, 0.0],
        omega0,
        INERTIA,
        h,
        steps,
        0.0,
        &DepSolveParams::default(),
    )
    .expect("conservative run completes");
    println!(
        "se3-004: energy drift {:e} (E0 {:e}), momentum drift {:e}, worst iters {}",
        receipt.energy_max_abs_drift,
        receipt.energy_start,
        receipt.momentum_max_abs_drift,
        receipt.max_solver_iters
    );
    assert_eq!(receipt.claim, Se3ClaimClass::ConservativeVariationalTheorem);
    assert!(receipt.all_solves_converged);
    // Spatial angular momentum is conserved by construction: drift is
    // pure roundoff accumulation.
    assert!(
        receipt.momentum_max_abs_drift < 1e-10,
        "momentum drift {:e}",
        receipt.momentum_max_abs_drift
    );
    // Energy oscillates bounded (no secular growth at this horizon).
    let rel = receipt.energy_max_abs_drift / receipt.energy_start;
    assert!(rel < 1e-4, "relative energy drift {rel:e}");
}

#[test]
fn se3_005_dep_adjoint_matches_finite_differences() {
    let omega0 = [0.9, -0.6, 0.3];
    let h = 5e-3;
    let steps = 25;
    let params = DepSolveParams::default();
    let bar_n = [0.25, -1.0, 0.5];
    let bar_0 = dep_momentum_adjoint(omega0, INERTIA, h, steps, &params, bar_n)
        .expect("adjoint runs");
    // Directional FD of <bar_n, omega_N(omega_0)> along d.
    let forward = |w0: [f64; 3]| -> [f64; 3] {
        let mut q = [1.0, 0.0, 0.0, 0.0];
        let mut w = w0;
        for _ in 0..steps {
            let (q1, w1, _) = dep_free_step(q, w, INERTIA, h, &params).expect("step");
            q = q1;
            w = w1;
        }
        w
    };
    let dirs = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.6, -0.8, 0.2]];
    for d in dirs {
        let eps = 1e-6;
        let mut plus = omega0;
        let mut minus = omega0;
        for a in 0..3 {
            plus[a] += eps * d[a];
            minus[a] -= eps * d[a];
        }
        let wp = forward(plus);
        let wm = forward(minus);
        let fd: f64 = (0..3)
            .map(|a| bar_n[a] * (wp[a] - wm[a]) / (2.0 * eps))
            .sum();
        let adj: f64 = (0..3).map(|a| bar_0[a] * d[a]).sum();
        let denom = fd.abs().max(adj.abs()).max(1e-12);
        let rel = (fd - adj).abs() / denom;
        println!("se3-005: dir {d:?} fd {fd:e} adjoint {adj:e} rel {rel:e}");
        assert!(rel < 1e-6, "adjoint-vs-FD gate failed: rel {rel:e}");
    }
}

#[test]
fn se3_006_damped_run_is_demoted_to_measured_only() {
    let (_, _, receipt) = run_dep_free(
        [1.0, 0.0, 0.0, 0.0],
        [0.7, 1.1, -0.4],
        INERTIA,
        1e-3,
        2_000,
        1e-4,
        &DepSolveParams::default(),
    )
    .expect("damped run completes");
    println!(
        "se3-006: damped claim {:?}, measured energy drift {:e}",
        receipt.claim, receipt.energy_max_abs_drift
    );
    // The honesty gate: dissipation must NOT inherit the theorem even
    // though every solve converged and the drift is smooth.
    assert_eq!(receipt.claim, Se3ClaimClass::MeasuredOnly);
    assert!(receipt.all_solves_converged);
    assert!(
        receipt.energy_max_abs_drift > 0.0,
        "damping must show up in the measured receipt"
    );
}

#[test]
fn se3_007_renormalization_receipts_bound_long_run_drift() {
    let twist = Twist {
        omega: [1.3, -0.8, 0.6],
        vel: [0.4, 0.9, -0.2],
    };
    let policy = RenormPolicy::default();
    let mut m = Motor::identity();
    let mut worst_defect = 0.0f64;
    let mut renorm_count = 0usize;
    let mut worst_drift = 0.0f64;
    for _ in 0..100_000 {
        let (next, receipt) = se3_exp_step_renorm(&m, &twist, 1e-3, &policy).expect("step");
        m = next;
        worst_defect = worst_defect.max(receipt.defect_before);
        if receipt.renormalized {
            renorm_count += 1;
            worst_drift = worst_drift.max(receipt.drift);
        }
    }
    println!(
        "se3-007: worst pre-decision defect {worst_defect:e}, renormalizations {renorm_count}, \
         worst drift {worst_drift:e}, final defect {:e}",
        m.unit_defect()
    );
    assert!(
        m.unit_defect() <= 1e-11,
        "final unit defect {:e} exceeds the receipt-controlled bound",
        m.unit_defect()
    );
}

#[test]
fn se3_008_rigid_body_step_tracks_so3_lane_and_conserves_free_velocity() {
    // Pure-rotation SE(3) rigid body vs the SO(3) reference lane.
    let mut motor = Motor::identity();
    let mut twist = Twist {
        omega: [0.7, 1.1, -0.4],
        vel: [0.0; 3],
    };
    let mut q = [1.0, 0.0, 0.0, 0.0];
    let mut w = twist.omega;
    let h = 1e-3;
    for _ in 0..2_000 {
        let (m1, t1) = se3_rigid_body_step(&motor, &twist, INERTIA, h).expect("step");
        motor = m1;
        twist = t1;
        let (q1, w1) = fs_time::rigid_body_step(q, w, INERTIA, h);
        q = q1;
        w = w1;
    }
    let werr = ((twist.omega[0] - w[0]).powi(2)
        + (twist.omega[1] - w[1]).powi(2)
        + (twist.omega[2] - w[2]).powi(2))
    .sqrt();
    println!("se3-008: angular-velocity agreement error {werr:e}");
    assert!(werr < 1e-12, "omega lanes disagree: {werr:e}");
    // Free translation: the SPATIAL velocity R·v_b must stay constant.
    let mut motor2 = Motor::identity();
    let mut twist2 = Twist {
        omega: [0.7, 1.1, -0.4],
        vel: [0.5, -0.2, 0.3],
    };
    let spatial0 = quat_rotate([1.0, 0.0, 0.0, 0.0], twist2.vel);
    for _ in 0..2_000 {
        let (m1, t1) = se3_rigid_body_step(&motor2, &twist2, INERTIA, h).expect("step");
        motor2 = m1;
        twist2 = t1;
    }
    // Reconstruct R from the motor by transporting basis directions.
    let origin = motor2
        .transform_point(GaPoint { x: 0.0, y: 0.0, z: 0.0 })
        .expect("finite");
    let rot_v = {
        let tip = motor2
            .transform_point(GaPoint {
                x: twist2.vel[0],
                y: twist2.vel[1],
                z: twist2.vel[2],
            })
            .expect("finite");
        [tip.x - origin.x, tip.y - origin.y, tip.z - origin.z]
    };
    let drift = ((rot_v[0] - spatial0[0]).powi(2)
        + (rot_v[1] - spatial0[1]).powi(2)
        + (rot_v[2] - spatial0[2]).powi(2))
    .sqrt();
    println!("se3-008: spatial free-velocity drift {drift:e} over 2000 midpoint steps");
    // Midpoint (order-2) integration: the drift is discretization
    // error, not conservation-by-construction; gate it at the
    // measured-order level rather than roundoff.
    assert!(drift < 5e-6, "spatial velocity drift {drift:e}");
}
