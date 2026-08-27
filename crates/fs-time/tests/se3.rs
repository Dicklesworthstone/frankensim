//! SE(3) lane battery (bead frankensim-ext-time-se3-lanes-3ol0).
//!
//! Detailed-logging discipline: every assertion window prints the
//! measured quantity it gates so a failure names the physics, not just
//! the boolean.

use fs_ga::{Motor, Se3, So3, Twist, Vec3};
use fs_time::se3::{
    DepSolveParams, RenormPolicy, Se3ClaimClass, Se3Error, dep_free_step, dep_momentum_adjoint,
    run_dep_free, se3_exp_step, se3_exp_step_renorm, se3_rigid_body_step, se3_space_exp_step,
};
use fs_time::{rigid_body_step, so3_body_exp_step, so3_space_exp_step};

const INERTIA: Vec3 = Vec3 {
    x: 1.0,
    y: 2.0,
    z: 3.0,
};

#[test]
fn se3_001_pure_rotation_agrees_with_so3_quaternion_lane() {
    let omega = Vec3::new(0.3, -0.7, 0.5);
    let h = 1e-2;
    let steps = 200;
    let mut pose = Se3::identity();
    let mut rotation = So3::identity();
    for _ in 0..steps {
        pose =
            se3_exp_step(pose, Twist::new(omega, Vec3::new(0.0, 0.0, 0.0)), h).expect("exp step");
        rotation = so3_body_exp_step(rotation, omega, h).expect("SO(3) step");
    }
    let probe = Vec3::new(0.8, -0.4, 1.2);
    let via_pose = pose.transform_point(probe).expect("finite point");
    let via_rotation = rotation.rotate(probe).expect("finite vector");
    let err = (via_pose - via_rotation).norm();
    println!("se3-001: SO(3)-lane agreement error {err:e}");
    assert!(err < 1e-11, "SE(3) and SO(3) lanes disagree: {err:e}");
}

#[test]
fn se3_002_constant_twist_steps_compose_to_one_exponential() {
    // One-parameter subgroup: N steps of h equal one step of N·h.
    let twist = Twist::new(Vec3::new(0.4, 0.1, -0.2), Vec3::new(0.5, -0.3, 0.8));
    let h = 0.05;
    let steps = 40;
    let mut walked = Se3::identity();
    for _ in 0..steps {
        walked = se3_exp_step(walked, twist, h).expect("exp step");
    }
    let jumped =
        se3_exp_step(Se3::identity(), twist, h * f64::from(steps)).expect("single exponential");
    let probe = Vec3::new(-0.6, 0.9, 0.3);
    let a = walked.transform_point(probe).expect("finite");
    let b = jumped.transform_point(probe).expect("finite");
    let err = (a - b).norm();
    println!("se3-002: screw one-parameter composition error {err:e}");
    assert!(err < 1e-10, "constant twist is not exact: {err:e}");
    let defect = {
        let mut m = *walked.as_motor();
        m.renormalize()
    };
    println!("se3-002: walked-motor renormalization drift {defect:e}");
}

#[test]
fn se3_003_double_cover_canonicalization_is_deterministic() {
    let m = Motor::rotor([0.0, 0.0, 1.0], 1.3).compose(&Motor::translator(0.2, -0.5, 0.7));
    let flipped = Motor(m.0.scale(-1.0));
    let ca = Se3::try_from_motor(m).expect("canonical");
    let cb = Se3::try_from_motor(flipped).expect("canonical");
    for (x, y) in ca.as_motor().0.0.iter().zip(cb.as_motor().0.0.iter()) {
        assert_eq!(x.to_bits(), y.to_bits(), "canonical representatives differ");
    }
    // Bitwise replay across two identical long runs that cross the
    // scalar-zero surface (rotation through pi).
    let twist = Twist::new(Vec3::new(0.0, 0.0, 2.5), Vec3::new(0.1, 0.0, 0.0));
    let run = || -> Vec<u64> {
        let mut pose = Se3::identity();
        let mut bits = Vec::new();
        for _ in 0..2000 {
            pose = se3_exp_step(pose, twist, 1e-3).expect("step");
            bits.push(pose.as_motor().0.0[0].to_bits());
        }
        bits
    };
    assert_eq!(run(), run(), "se3-003: replay is not bitwise deterministic");
}

#[test]
fn se3_004_dep_free_body_conserves_momentum_and_bounds_energy() {
    let omega0 = Vec3::new(0.7, 1.1, -0.4);
    let h = 1e-3;
    let steps = 10_000;
    let (_, _, receipt) = run_dep_free(
        So3::identity(),
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
    let omega0 = Vec3::new(0.9, -0.6, 0.3);
    let h = 5e-3;
    let steps = 25;
    let params = DepSolveParams::default();
    let bar_n = Vec3::new(0.25, -1.0, 0.5);
    let bar_0 =
        dep_momentum_adjoint(omega0, INERTIA, h, steps, &params, bar_n).expect("adjoint runs");
    // Directional FD of <bar_n, omega_N(omega_0)> along d.
    let forward = |w0: Vec3| -> Vec3 {
        let mut rotation = So3::identity();
        let mut w = w0;
        for _ in 0..steps {
            let (next_rotation, w1, _) =
                dep_free_step(rotation, w, INERTIA, h, &params).expect("step");
            rotation = next_rotation;
            w = w1;
        }
        w
    };
    let dirs = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.6, -0.8, 0.2],
    ];
    for direction in dirs {
        let direction = Vec3::new(direction[0], direction[1], direction[2]);
        let eps = 1e-6;
        let plus = omega0 + direction.scale(eps);
        let minus = omega0 - direction.scale(eps);
        let wp = forward(plus);
        let wm = forward(minus);
        let fd = bar_n.dot((wp - wm).scale(1.0 / (2.0 * eps)));
        let adj = bar_0.dot(direction);
        let denom = fd.abs().max(adj.abs()).max(1e-12);
        let rel = (fd - adj).abs() / denom;
        println!("se3-005: dir {direction:?} fd {fd:e} adjoint {adj:e} rel {rel:e}");
        assert!(rel < 1e-6, "adjoint-vs-FD gate failed: rel {rel:e}");
    }
}

#[test]
fn se3_006_damped_run_is_demoted_to_measured_only() {
    let (_, _, receipt) = run_dep_free(
        So3::identity(),
        Vec3::new(0.7, 1.1, -0.4),
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
    let twist = Twist::new(Vec3::new(1.3, -0.8, 0.6), Vec3::new(0.4, 0.9, -0.2));
    let policy = RenormPolicy::default();
    let mut pose = Se3::identity();
    let mut worst_defect = 0.0f64;
    let mut renorm_count = 0usize;
    let mut worst_drift = 0.0f64;
    for _ in 0..100_000 {
        let (next, receipt) = se3_exp_step_renorm(pose, twist, 1e-3, &policy).expect("step");
        pose = next;
        worst_defect = worst_defect.max(receipt.defect_before);
        if receipt.renormalized {
            renorm_count += 1;
            worst_drift = worst_drift.max(receipt.drift);
        }
    }
    println!(
        "se3-007: worst pre-decision defect {worst_defect:e}, renormalizations {renorm_count}, \
         worst drift {worst_drift:e}, final defect {:e}",
        pose.as_motor().unit_defect()
    );
    assert!(
        pose.as_motor().unit_defect() <= 1e-11,
        "final unit defect {:e} exceeds the receipt-controlled bound",
        pose.as_motor().unit_defect()
    );
}

#[test]
fn se3_008_rigid_body_step_tracks_so3_lane_and_conserves_free_velocity() {
    // Pure-rotation SE(3) rigid body vs the SO(3) reference lane.
    let mut pose = Se3::identity();
    let mut twist = Twist::new(Vec3::new(0.7, 1.1, -0.4), Vec3::new(0.0, 0.0, 0.0));
    let mut rotation = So3::identity();
    let mut angular = twist.angular;
    let h = 1e-3;
    for _ in 0..2_000 {
        let (next_pose, next_twist) = se3_rigid_body_step(pose, twist, INERTIA, h).expect("step");
        pose = next_pose;
        twist = next_twist;
        let (next_rotation, next_angular) =
            rigid_body_step(rotation, angular, INERTIA, h).expect("SO(3) step");
        rotation = next_rotation;
        angular = next_angular;
    }
    let werr = (twist.angular - angular).norm();
    println!("se3-008: angular-velocity agreement error {werr:e}");
    assert!(werr < 1e-12, "omega lanes disagree: {werr:e}");
    // Free translation: the SPATIAL velocity R·v_b must stay constant.
    let mut pose2 = Se3::identity();
    let mut twist2 = Twist::new(Vec3::new(0.7, 1.1, -0.4), Vec3::new(0.5, -0.2, 0.3));
    let spatial0 = twist2.linear;
    for _ in 0..2_000 {
        let (next_pose, next_twist) = se3_rigid_body_step(pose2, twist2, INERTIA, h).expect("step");
        pose2 = next_pose;
        twist2 = next_twist;
    }
    let spatial_velocity = pose2
        .rotation()
        .rotate(twist2.linear)
        .expect("finite spatial velocity");
    let drift = (spatial_velocity - spatial0).norm();
    println!("se3-008: spatial free-velocity drift {drift:e} over 2000 midpoint steps");
    // Midpoint (order-2) integration: the drift is discretization
    // error, not conservation-by-construction; gate it at the
    // measured-order level rather than roundoff.
    assert!(drift < 5e-6, "spatial velocity drift {drift:e}");
}

#[test]
fn se3_009_body_right_and_space_left_steps_are_explicit_and_equivalent() {
    let pose = Se3::exp(Twist::new(
        Vec3::new(0.2, -0.3, 0.4),
        Vec3::new(0.8, -0.5, 0.1),
    ))
    .expect("fixture pose");
    let body_twist = Twist::new(Vec3::new(-0.4, 0.2, 0.1), Vec3::new(0.3, 0.7, -0.2));
    let space_twist = body_twist.transform_by(&pose);
    let h = 0.037;
    let via_body = se3_exp_step(pose, body_twist, h).expect("body/right step");
    let via_space = se3_space_exp_step(pose, space_twist, h).expect("space/left step");
    let residual = via_body.body_minus(via_space).expect("group difference");
    assert!(
        residual.angular.norm().max(residual.linear.norm()) < 3.0e-12,
        "SE(3) body/space convention mismatch: {residual:?}"
    );

    let rotation = pose.rotation();
    let omega_body = body_twist.angular;
    let omega_space = rotation
        .rotate(omega_body)
        .expect("finite angular velocity");
    let via_body = so3_body_exp_step(rotation, omega_body, h).expect("body/right SO(3) step");
    let via_space = so3_space_exp_step(rotation, omega_space, h).expect("space/left SO(3) step");
    let residual = via_body
        .body_minus(via_space)
        .expect("SO(3) group difference");
    assert!(
        residual.angular.norm() < 2.0e-12,
        "SO(3) body/space convention mismatch: {residual:?}"
    );

    let run = || {
        let mut body_pose = pose;
        let mut space_pose = pose;
        let mut body_rotation = rotation;
        let mut space_rotation = rotation;
        for _ in 0..256 {
            body_pose = se3_exp_step(body_pose, body_twist, h).expect("body replay step");
            space_pose = se3_space_exp_step(space_pose, space_twist, h).expect("space replay step");
            body_rotation =
                so3_body_exp_step(body_rotation, omega_body, h).expect("SO(3) body replay step");
            space_rotation = so3_space_exp_step(space_rotation, omega_space, h)
                .expect("SO(3) space replay step");
        }
        (
            body_pose.as_motor().0.0.map(f64::to_bits),
            space_pose.as_motor().0.0.map(f64::to_bits),
            {
                let q = body_rotation.as_quat();
                [q.w.to_bits(), q.x.to_bits(), q.y.to_bits(), q.z.to_bits()]
            },
            {
                let q = space_rotation.as_quat();
                [q.w.to_bits(), q.x.to_bits(), q.y.to_bits(), q.z.to_bits()]
            },
        )
    };
    assert_eq!(run(), run(), "body/space stepping did not replay bitwise");
}

#[test]
fn se3_010_invalid_variational_controls_refuse_even_zero_step_runs() {
    let invalid = DepSolveParams {
        tol: f64::NAN,
        max_iters: 0,
    };
    let result = run_dep_free(
        So3::identity(),
        Vec3::new(0.7, 1.1, -0.4),
        INERTIA,
        1e-3,
        0,
        0.0,
        &invalid,
    );
    assert!(
        matches!(
            result,
            Err(Se3Error::InvalidParameter {
                context: "variational solve controls"
            })
        ),
        "invalid controls must refuse before an empty trajectory can bypass validation: {result:?}"
    );

    let excessive_damping = run_dep_free(
        So3::identity(),
        Vec3::new(0.7, 1.1, -0.4),
        INERTIA,
        1e-3,
        1,
        1.01,
        &DepSolveParams::default(),
    );
    assert!(
        matches!(
            excessive_damping,
            Err(Se3Error::InvalidParameter {
                context: "per-step damping fraction"
            })
        ),
        "a damping fraction outside [0, 1] must refuse: {excessive_damping:?}"
    );
}
