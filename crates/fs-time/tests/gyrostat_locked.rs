//! Retained shaft locks: original kinetic form, carrier recoil and native drift.
use fs_ga::{Mat3, Quat, So3, So3Tangent, Vec3};
use fs_time::lie::gyrostat::rotors::{AxialRotor, RotorBank, RotorState};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    assert!((a - b).norm() <= tolerance, "{a:?} != {b:?}");
}
fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn attitude(a: So3, b: So3, tolerance: f64) {
    assert!(a.body_minus(b).unwrap().angular.norm() <= tolerance);
}

// Independent adjugate inverse, used only as a 3x3 analytic reference.
fn inverse(a: Mat3) -> Mat3 {
    let m = a.m;
    let adj = [
        m[4] * m[8] - m[5] * m[7],
        m[2] * m[7] - m[1] * m[8],
        m[1] * m[5] - m[2] * m[4],
        m[5] * m[6] - m[3] * m[8],
        m[0] * m[8] - m[2] * m[6],
        m[2] * m[3] - m[0] * m[5],
        m[3] * m[7] - m[4] * m[6],
        m[1] * m[6] - m[0] * m[7],
        m[0] * m[4] - m[1] * m[3],
    ];
    let determinant = m[0] * adj[0] + m[1] * adj[3] + m[2] * adj[6];
    Mat3 {
        m: adj.map(|x| x / determinant),
    }
}

fn mode_tensor<const N: usize>(bank: &RotorBank<N>, locked: [bool; N]) -> Mat3 {
    let mut matrix = bank.reduced_inertia().matrix().m;
    for (i, rotor) in bank.rotors().iter().enumerate() {
        if locked[i] {
            let axis = rotor.axis_body();
            let a = [axis.x, axis.y, axis.z];
            for r in 0..3 {
                for c in 0..3 {
                    matrix[3 * r + c] += rotor.axial_inertia_kg_m2() * a[r] * a[c];
                }
            }
        }
    }
    Mat3 { m: matrix }
}

fn fixture() -> (RotorBank<3>, RotorState<3>) {
    let reduced = InertiaTensor3::from_principal(
        Vec3::new(2.0, 3.0, 4.0),
        rotation(Vec3::new(0.4, -0.3, 0.2)),
    )
    .unwrap();
    let rotors = [
        AxialRotor::new(Vec3::new(1.0, 0.3, 0.2), 0.6).unwrap(),
        AxialRotor::new(Vec3::new(0.2, 1.0, -0.4), 0.8).unwrap(),
        AxialRotor::new(Vec3::new(0.4, -0.2, 1.0), 0.5).unwrap(),
    ];
    let bank = RotorBank::from_reduced(reduced, rotors).unwrap();
    let state = bank
        .state_from_velocities(
            rotation(Vec3::new(-0.2, 0.4, 0.3)),
            Vec3::new(0.7, -0.5, 0.8),
            [1.4, -0.9, 0.6],
            [0.2, -0.3, 0.4],
        )
        .unwrap();
    (bank, state)
}

#[test]
fn coaxial_engagement_has_recoil_loss_fixed_phase_and_zero_work_holding() {
    let bank = RotorBank::from_reduced(
        InertiaTensor3::new(Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0],
        })
        .unwrap(),
        [AxialRotor::new(Vec3::new(0.0, 0.0, 1.0), 1.0).unwrap()],
    )
    .unwrap();
    let start = bank
        .state_from_velocities(
            rotation(Vec3::new(0.2, -0.3, 0.4)),
            Vec3::new(0.0, 0.0, 0.25),
            [6.0],
            [123.0],
        )
        .unwrap();
    let mode = bank.lock_mode([true]).unwrap();
    let engagement = mode.engage(&start).unwrap();
    let stopped = engagement.impulse.state_after;
    vector(
        bank.angular_velocity(&stopped).unwrap(),
        Vec3::new(0.0, 0.0, 1.75),
        1e-15,
    );
    close(stopped.axial_momentum_nms[0], 1.75, 1e-15);
    close(engagement.impulse.axial_momentum_change_nms[0], -4.5, 1e-15);
    close(engagement.dissipated_energy_j, 13.5, 1e-14);
    close(
        bank.kinetic_energy(&start).unwrap() - bank.kinetic_energy(&stopped).unwrap(),
        13.5,
        1e-14,
    );
    close(engagement.work_energy_residual_j, 0.0, 1e-14);
    assert_eq!(stopped.total_momentum_body, start.total_momentum_body);
    assert_eq!(stopped.relative_angle_rad, start.relative_angle_rad);
    assert_eq!(stopped.rotation, start.rotation);
    let repeated = mode.engage(&stopped).unwrap();
    assert_eq!(repeated.impulse.state_after, stopped);
    assert_eq!(repeated.dissipated_energy_j, 0.0);

    let step = mode.drift(&stopped, 0.75).unwrap();
    let expected = stopped
        .rotation
        .compose(rotation(Vec3::new(0.0, 0.0, 1.75 * 0.75)))
        .unwrap();
    attitude(step.state_after.rotation, expected, 1e-15);
    assert_eq!(step.state_after.relative_angle_rad, [123.0]);
    assert_eq!(
        step.state_after.axial_momentum_nms,
        stopped.axial_momentum_nms
    );
    close(step.drift_energy_defect_j, 0.0, 1e-15);
    let kick = mode
        .apply_impulse(&step.state_after, [0.0], Vec3::new(0.0, 0.0, 0.8))
        .unwrap();
    let after = kick.impulse.state_after;
    vector(
        bank.angular_velocity(&after).unwrap(),
        Vec3::new(0.0, 0.0, 1.95),
        2e-15,
    );
    close(after.axial_momentum_nms[0], 1.95, 2e-15);
    close(kick.impulse.external_work_j, 1.48, 2e-15);
    close(kick.reaction_work_j, 0.0, 2e-15);
    close(kick.work_energy_residual_j, 0.0, 3e-15);
    assert_eq!(after.relative_angle_rad, [123.0]);
    assert_eq!(start.axial_momentum_nms, [6.25]);
}

#[test]
fn simultaneous_nonorthogonal_locks_and_impulse_ports_use_the_complete_tensor() {
    let (bank, start) = fixture();
    let locks = [true, false, true];
    let mode = bank.lock_mode(locks).unwrap();
    let stopped = mode.engage(&start).unwrap();
    let zero = mode
        .apply_impulse(&stopped.impulse.state_after, [0.0; 3], Vec3::default())
        .unwrap();
    assert_eq!(zero.impulse.state_after, stopped.impulse.state_after);
    let expected_inverse = inverse(mode_tensor(&bank, locks));
    let stored = bank.rotors()[1]
        .axis_body()
        .scale(start.axial_momentum_nms[1]);
    let expected = expected_inverse.apply(start.total_momentum_body - stored);
    vector(
        bank.angular_velocity(&stopped.impulse.state_after).unwrap(),
        expected,
        2e-14,
    );
    let rates = bank.relative_rates(&stopped.impulse.state_after).unwrap();
    close(rates[0], 0.0, 3e-14);
    close(rates[2], 0.0, 3e-14);
    assert_eq!(
        stopped.impulse.state_after.axial_momentum_nms[1],
        start.axial_momentum_nms[1]
    );
    assert!(stopped.dissipated_energy_j > 0.0);
    close(
        bank.kinetic_energy(&start).unwrap()
            - bank.kinetic_energy(&stopped.impulse.state_after).unwrap(),
        stopped.dissipated_energy_j,
        5e-14,
    );

    let permutation = [2, 0, 1];
    let other_bank = RotorBank::from_reduced(
        bank.reduced_inertia(),
        permutation.map(|i| bank.rotors()[i]),
    )
    .unwrap();
    let other_start = RotorState {
        axial_momentum_nms: permutation.map(|i| start.axial_momentum_nms[i]),
        relative_angle_rad: permutation.map(|i| start.relative_angle_rad[i]),
        ..start
    };
    let other = other_bank
        .lock_mode(permutation.map(|i| locks[i]))
        .unwrap()
        .engage(&other_start)
        .unwrap();
    vector(
        other_bank
            .angular_velocity(&other.impulse.state_after)
            .unwrap(),
        expected,
        3e-14,
    );
    close(
        other.dissipated_energy_j,
        stopped.dissipated_energy_j,
        3e-14,
    );
    for (index, original) in permutation.into_iter().enumerate() {
        close(
            other.impulse.state_after.axial_momentum_nms[index],
            stopped.impulse.state_after.axial_momentum_nms[original],
            3e-14,
        );
    }

    let external = Vec3::new(0.03, -0.02, 0.04);
    let kick = mode
        .apply_impulse(&stopped.impulse.state_after, [0.0, 0.15, 0.0], external)
        .unwrap();
    let after = kick.impulse.state_after;
    let stored = bank.rotors()[1]
        .axis_body()
        .scale(after.axial_momentum_nms[1]);
    vector(
        bank.angular_velocity(&after).unwrap(),
        expected_inverse.apply(after.total_momentum_body - stored),
        3e-14,
    );
    close(
        after.axial_momentum_nms[1] - start.axial_momentum_nms[1],
        0.15,
        1e-15,
    );
    vector(
        after.total_momentum_body - start.total_momentum_body,
        external,
        1e-15,
    );
    close(kick.reaction_work_j, 0.0, 2e-15);
    close(kick.work_energy_residual_j, 0.0, 5e-15);
    close(kick.drive_work_j, kick.impulse.shaft_work_j[1], 0.0);
}

// Independent original-frame Euler equations, quaternion attitude kinematics
// and unlocked relative phases. No production gyrostat stepping is used.
fn reference(
    bank: &RotorBank<3>,
    state: &RotorState<3>,
    locks: [bool; 3],
    duration: f64,
) -> [f64; 10] {
    let inv = inverse(mode_tensor(bank, locks));
    let mut stored = Vec3::default();
    for i in 0..3 {
        if !locks[i] {
            stored = stored
                + bank.rotors()[i]
                    .axis_body()
                    .scale(state.axial_momentum_nms[i]);
        }
    }
    let q = state.rotation.as_quat();
    let l = state.total_momentum_body;
    let mut y = [
        q.w,
        q.x,
        q.y,
        q.z,
        l.x,
        l.y,
        l.z,
        state.relative_angle_rad[0],
        state.relative_angle_rad[1],
        state.relative_angle_rad[2],
    ];
    let rhs = |s: [f64; 10]| {
        let momentum = Vec3::new(s[4], s[5], s[6]);
        let w = inv.apply(momentum - stored);
        let dl = momentum.cross(w);
        let phase = core::array::from_fn::<_, 3, _>(|i| {
            if locks[i] {
                0.0
            } else {
                state.axial_momentum_nms[i] / bank.rotors()[i].axial_inertia_kg_m2()
                    - bank.rotors()[i].axis_body().dot(w)
            }
        });
        [
            -0.5 * (s[1] * w.x + s[2] * w.y + s[3] * w.z),
            0.5 * (s[0] * w.x + s[2] * w.z - s[3] * w.y),
            0.5 * (s[0] * w.y + s[3] * w.x - s[1] * w.z),
            0.5 * (s[0] * w.z + s[1] * w.y - s[2] * w.x),
            dl.x,
            dl.y,
            dl.z,
            phase[0],
            phase[1],
            phase[2],
        ]
    };
    let h = duration / 8192.0;
    for _ in 0..8192 {
        let k1 = rhs(y);
        let k2 = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]));
        let k3 = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]));
        let k4 = rhs(core::array::from_fn(|i| y[i] + h * k3[i]));
        for i in 0..10 {
            y[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
    y
}

#[test]
fn retained_full_tensor_precession_refines_with_zero_locked_phase_and_no_artificial_heat() {
    let (bank, moving) = fixture();
    let locks = [true, false, true];
    let mode = bank.lock_mode(locks).unwrap();
    let start = mode.engage(&moving).unwrap().impulse.state_after;
    let duration = 0.4;
    let oracle = reference(&bank, &start, locks, duration);
    let expected_rotation = So3::try_from_quat(Quat {
        w: oracle[0],
        x: oracle[1],
        y: oracle[2],
        z: oracle[3],
    })
    .unwrap();
    let expected_momentum = Vec3::new(oracle[4], oracle[5], oracle[6]);
    let spatial = start.rotation.rotate(start.total_momentum_body).unwrap();
    let mut errors = Vec::new();
    for count in [20, 40, 80] {
        let mut state = start;
        let mut measured_change = 0.0;
        let mut measured_native = 0.0;
        for _ in 0..count {
            let step = mode.drift(&state, duration / count as f64).unwrap();
            measured_change += step.kinetic_energy_change_j;
            measured_native += step.drift_energy_defect_j;
            assert!(step.angular_momentum_residual_world.norm() < 2e-14);
            assert!(step.work_energy_residual_j.abs() < 2e-14);
            state = step.state_after;
            assert_eq!(state.relative_angle_rad[0], start.relative_angle_rad[0]);
            assert_eq!(state.relative_angle_rad[2], start.relative_angle_rad[2]);
            assert_eq!(state.axial_momentum_nms[1], start.axial_momentum_nms[1]);
            let rates = bank.relative_rates(&state).unwrap();
            close(rates[0], 0.0, 3e-14);
            close(rates[2], 0.0, 3e-14);
        }
        vector(
            state.rotation.rotate(state.total_momentum_body).unwrap(),
            spatial,
            1e-12,
        );
        close(
            bank.kinetic_energy(&state).unwrap() - bank.kinetic_energy(&start).unwrap(),
            measured_change,
            4e-14,
        );
        close(measured_change, measured_native, 4e-13);
        let error = state
            .rotation
            .body_minus(expected_rotation)
            .unwrap()
            .angular
            .norm()
            + (state.total_momentum_body - expected_momentum).norm()
            + (state.relative_angle_rad[1] - oracle[8]).abs();
        errors.push(error);
        println!(
            "locked Euler refinement n={count}: error={error:.16e}, native energy defect={measured_native:.16e}"
        );
    }
    assert!(errors[1] < 0.27 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.27 * errors[1], "{errors:?}");
    let forward = mode.drift(&start, 0.02).unwrap().state_after;
    let backward = mode.drift(&forward, -0.02).unwrap().state_after;
    vector(
        backward.total_momentum_body,
        start.total_momentum_body,
        4e-14,
    );
    attitude(backward.rotation, start.rotation, 4e-14);
    for i in 0..3 {
        close(
            backward.axial_momentum_nms[i],
            start.axial_momentum_nms[i],
            4e-14,
        );
        close(
            backward.relative_angle_rad[i],
            start.relative_angle_rad[i],
            4e-14,
        );
    }
}

#[test]
fn release_passes_the_accepted_canonical_state_to_free_dynamics_without_a_reset() {
    let (bank, moving) = fixture();
    let mode = bank.lock_mode([true, false, true]).unwrap();
    let start = mode.engage(&moving).unwrap().impulse.state_after;
    let stopped = mode.drift(&start, 0.1).unwrap().state_after;
    let free = bank
        .step(&stopped, [0.0; 3], HeldTorque::Body(Vec3::default()), 0.1)
        .unwrap();
    let retained = mode.drift(&stopped, 0.1).unwrap();
    assert_eq!(
        free.state_after.axial_momentum_nms,
        stopped.axial_momentum_nms
    );
    assert_eq!(
        retained.state_after.relative_angle_rad[0],
        stopped.relative_angle_rad[0]
    );
    assert!((free.state_after.relative_angle_rad[0] - stopped.relative_angle_rad[0]).abs() > 1e-6);
    assert!(
        (free.state_after.axial_momentum_nms[0] - retained.state_after.axial_momentum_nms[0]).abs()
            > 1e-5
    );
    vector(
        free.state_after
            .rotation
            .rotate(free.state_after.total_momentum_body)
            .unwrap(),
        stopped
            .rotation
            .rotate(stopped.total_momentum_body)
            .unwrap(),
        3e-14,
    );
}

#[test]
fn empty_mode_reduces_to_native_free_motion_and_unresolved_engagement_is_atomic() {
    let (bank, start) = fixture();
    let mode = bank.lock_mode([false; 3]).unwrap();
    let unchanged = mode.engage(&start).unwrap();
    assert_eq!(unchanged.impulse.state_after, start);
    assert_eq!(unchanged.dissipated_energy_j, 0.0);
    let ordinary = bank
        .step(&start, [0.0; 3], HeldTorque::Body(Vec3::default()), 0.1)
        .unwrap();
    let empty = mode.drift(&start, 0.1).unwrap();
    assert_eq!(empty.state_after, ordinary.state_after);
    let held = bank.lock_mode([true, false, true]).unwrap();
    assert!(held.drift(&start, 0.1).is_err());
    let stopped = held.engage(&start).unwrap().impulse.state_after;
    assert!(held.drift(&stopped, 0.0).is_err());
    assert!(held.drift(&stopped, f64::NAN).is_err());
    assert!(
        held.apply_impulse(&stopped, [0.1, 0.0, 0.0], Vec3::default())
            .is_err()
    );
    assert!(
        held.apply_impulse(&stopped, [0.0, f64::INFINITY, 0.0], Vec3::default())
            .is_err()
    );
    let rotor = AxialRotor::new(Vec3::new(1.0, 0.0, 0.0), 1.0).unwrap();
    let bank = RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(2.0, 3.0, 4.0), So3::identity()).unwrap(),
        [rotor, rotor],
    )
    .unwrap();
    // Locked complete I=diag(4,3,4) is geometrically plausible. A large free
    // rotor erases the smaller lock transfer in aggregate h. Refuse, rather
    // than claiming zero carrier recoil or inventing the engagement loss.
    let unresolved = bank
        .state_from_velocities(So3::identity(), Vec3::default(), [1e16, 1.0], [0.0; 2])
        .unwrap();
    let saved = unresolved;
    assert!(
        bank.lock_mode([false, true])
            .unwrap()
            .engage(&unresolved)
            .is_err()
    );
    assert_eq!(unresolved, saved);
    assert_eq!(start, fixture().1);
}
