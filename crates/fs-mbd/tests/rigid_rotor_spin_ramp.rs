//! G1/G3/G4: physical spin acceleration, transverse angular momentum and continuation.
use fs_ga::{Mat3, So3, Vec3};
use fs_mbd::articulated::SpatialInertia;
use fs_mbd::rotordynamics::RotorError;
use fs_mbd::rotordynamics::rigid::{
    RadialBearing, RigidRotor, RotorSpinDriveAxis, RotorSpinRamp, RotorSpinState, RotorUnbalance,
};
use fs_mbd::rotordynamics::transient::RotorTransientState;
use std::cell::Cell;

fn rotor(bearings: &[RadialBearing]) -> RigidRotor {
    let inertia = SpatialInertia::new(
        4.0,
        Vec3::default(),
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0],
        },
    )
    .unwrap();
    RigidRotor::new(inertia, So3::identity(), bearings, 1e-12, || false).unwrap()
}

fn ramp(alpha: f64, h: f64, drive_axis: RotorSpinDriveAxis) -> RotorSpinRamp {
    RotorSpinRamp {
        angular_acceleration_rad_s2: alpha,
        duration_seconds: h,
        drive_axis,
        relative_tolerance: 1e-11,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn zero_state() -> RotorTransientState {
    RotorTransientState {
        displacement: vec![0.0; 4],
        velocity: vec![0.0; 4],
    }
}

fn bearings() -> [RadialBearing; 2] {
    [
        RadialBearing {
            station_m: -0.5,
            stiffness_n_per_m: [[17.0, -3.0], [2.0, 23.0]],
            damping_ns_per_m: [[0.4, -0.3], [0.2, -0.1]],
        },
        RadialBearing {
            station_m: 0.8,
            stiffness_n_per_m: [[11.0, 1.5], [-2.5, 19.0]],
            damping_ns_per_m: [[0.2, 0.1], [-0.4, 0.3]],
        },
    ]
}

#[test]
fn g3_constant_spin_is_the_original_native_step_with_all_bearing_entries() {
    let rotor = rotor(&bearings());
    let source = rotor
        .prepare_unbalance(
            &[RotorUnbalance {
                station_m: 0.3,
                mass_eccentricity_kg_m: 0.02,
                phase_rad: -0.4,
            }],
            || false,
        )
        .unwrap();
    let state = RotorTransientState {
        displacement: vec![0.03, -0.02, 0.04, 0.01],
        velocity: vec![-0.07, 0.05, -0.03, 0.02],
    };
    let additional = [0.2, -0.3, 0.1, -0.05];
    for speed in [-3.0, 0.0, 4.0] {
        let spin = RotorSpinState {
            time_seconds: 7.0,
            phase_rad: 0.7,
            spin_rad_s: speed,
        };
        for axis in [RotorSpinDriveAxis::Shaft, RotorSpinDriveAxis::Stationary] {
            let result = rotor
                .spin_ramp_step(
                    &state,
                    spin,
                    &source,
                    additional,
                    ramp(0.0, 0.02, axis),
                    || false,
                )
                .unwrap();
            let unbalance = source
                .force_at(result.spin_midpoint.phase_rad, speed, 0.0, || false)
                .unwrap();
            let force = std::array::from_fn::<_, 4, _>(|i| additional[i] + unbalance[i]);
            let expected = rotor
                .model()
                .at_speed(speed, || false)
                .unwrap()
                .midpoint_stepper(0.02, 1e-11, || false)
                .unwrap()
                .step(&state, &force, || false)
                .unwrap();
            assert_eq!(result.motion, expected);
            assert_eq!(result.unbalance_force, unbalance);
            assert_eq!(result.acceleration_coupling_work_j, 0.0);
            assert_eq!(
                result.bearing_stiffness_work_j,
                expected.stiffness_force_work_j
            );
            near(result.spin_after.phase_rad, 0.7 + 0.02 * speed, 1e-16);
            assert_eq!(result.spin_after.spin_rad_s, speed);
        }
    }
}

fn free_tilt(steps: usize, axis: RotorSpinDriveAxis) -> RotorTransientState {
    let rotor = rotor(&[]);
    let source = rotor.prepare_unbalance(&[], || false).unwrap();
    let mut state = zero_state();
    state.displacement[2] = 0.05;
    state.displacement[3] = -0.03;
    let mut spin = RotorSpinState {
        time_seconds: 0.0,
        phase_rad: 1.3,
        spin_rad_s: 0.0,
    };
    for _ in 0..steps {
        let result = rotor
            .spin_ramp_step(
                &state,
                spin,
                &source,
                [0.0; 4],
                ramp(4.0, 0.8 / steps as f64, axis),
                || false,
            )
            .unwrap();
        near(result.motion.work_energy_residual_j, 0.0, 1e-17);
        near(result.bearing_stiffness_work_j, 0.0, 1e-17);
        near(
            result.motion.stiffness_force_work_j,
            result.acceleration_coupling_work_j,
            1e-17,
        );
        state = result.motion.state_after;
        spin = result.spin_after;
        if axis == RotorSpinDriveAxis::Stationary {
            // No transverse torque: actual inertial angular momentum stays
            // zero even though spin and transverse angular velocity grow.
            near(
                2.0 * state.velocity[2] + 3.0 * spin.spin_rad_s * state.displacement[3],
                0.0,
                2e-15,
            );
            near(
                2.0 * state.velocity[3] - 3.0 * spin.spin_rad_s * state.displacement[2],
                0.0,
                2e-15,
            );
        }
    }
    near(spin.time_seconds, 0.8, 2e-15);
    near(spin.phase_rad, 1.3 + 0.5 * 4.0 * 0.8 * 0.8, 4e-15);
    near(spin.spin_rad_s, 3.2, 4e-15);
    state
}

#[test]
fn g1_drive_direction_controls_euler_acceleration_and_matches_exact_free_precession() {
    // H_perp=0 gives theta_dot=-(Jp/Jd)*Omega*J*theta, hence an
    // independently known planar rotation by (Jp/Jd)*integral(Omega dt).
    let angle: f64 = 1.5 * 0.5 * 4.0 * 0.8 * 0.8;
    let (s, c) = angle.sin_cos();
    let theta = [c * 0.05 + s * 0.03, s * 0.05 - c * 0.03];
    let exact = [
        theta[0],
        theta[1],
        -1.5 * 3.2 * theta[1],
        1.5 * 3.2 * theta[0],
    ];
    let error = |steps| {
        let state = free_tilt(steps, RotorSpinDriveAxis::Stationary);
        [
            state.displacement[2],
            state.displacement[3],
            state.velocity[2],
            state.velocity[3],
        ]
        .into_iter()
        .zip(exact)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
    };
    let coarse = error(20);
    let medium = error(40);
    let fine = error(80);
    assert!(coarse / medium > 3.97, "{coarse:e} / {medium:e}");
    assert!(medium / fine > 3.99, "{medium:e} / {fine:e}");
    assert!(fine < 1.1e-4, "{fine:e}");
    let shaft_driven = free_tilt(20, RotorSpinDriveAxis::Shaft);
    assert_eq!(shaft_driven.displacement, vec![0.0, 0.0, 0.05, -0.03]);
    assert_eq!(shaft_driven.velocity, vec![0.0; 4]);
}

fn translating_unbalance(steps: usize, speed0: f64, alpha: f64) -> RotorTransientState {
    let rotor = rotor(&[]);
    let source = rotor
        .prepare_unbalance(
            &[RotorUnbalance {
                station_m: 0.0,
                mass_eccentricity_kg_m: 0.02,
                phase_rad: 0.3,
            }],
            || false,
        )
        .unwrap();
    let mut state = zero_state();
    let mut spin = RotorSpinState {
        time_seconds: 11.0,
        phase_rad: -0.2,
        spin_rad_s: speed0,
    };
    for _ in 0..steps {
        let h = 0.8 / steps as f64;
        let result = rotor
            .spin_ramp_step(
                &state,
                spin,
                &source,
                [0.0; 4],
                ramp(alpha, h, RotorSpinDriveAxis::Shaft),
                || false,
            )
            .unwrap();
        // Quadratic phase must not be replaced by the mean of endpoints.
        near(
            result.spin_midpoint.phase_rad,
            spin.phase_rad + spin.spin_rad_s * 0.5 * h + 0.5 * alpha * (0.5 * h).powi(2),
            3e-16,
        );
        state = result.motion.state_after;
        spin = result.spin_after;
    }
    near(
        spin.phase_rad,
        -0.2 + speed0 * 0.8 + 0.5 * alpha * 0.8 * 0.8,
        6e-15,
    );
    near(spin.spin_rad_s, speed0 + alpha * 0.8, 6e-15);
    assert_eq!(&state.displacement[2..], &[0.0; 2]);
    assert_eq!(&state.velocity[2..], &[0.0; 2]);
    state
}

#[test]
fn g1_unbalance_startup_braking_and_reversal_match_exact_inertial_translation() {
    // Each case's bound sits just above the midpoint rule's own error for its
    // force history at 80 steps: v += h a(t_mid), q += h (v0 + v1) / 2 run
    // independently against the closed form below gives 1.0934e-6 and
    // 1.3934e-6. One bound of 1.2e-6 lay below the second case's.
    for (speed0, alpha, fine_bound) in [(0.0, 4.0, 1.1e-6), (2.0, -5.0, 1.4e-6)] {
        let duration = 0.8;
        let phase0 = 0.1_f64;
        let phase = phase0 + speed0 * duration + 0.5 * alpha * duration * duration;
        let speed = speed0 + alpha * duration;
        let (s0, c0) = phase0.sin_cos();
        let (s, c) = phase.sin_cos();
        // Integrating M q'' = -U d²(e_r)/dt² from q=q'=0 also
        // retains the initial eccentric velocity when starting at nonzero spin.
        let ratio = 0.02 / 4.0;
        let exact = [
            -ratio * (c - c0 + duration * speed0 * s0),
            -ratio * (s - s0 - duration * speed0 * c0),
            ratio * (speed * s - speed0 * s0),
            -ratio * (speed * c - speed0 * c0),
        ];
        let error = |steps| {
            let state = translating_unbalance(steps, speed0, alpha);
            [
                state.displacement[0],
                state.displacement[1],
                state.velocity[0],
                state.velocity[1],
            ]
            .into_iter()
            .zip(exact)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
        };
        let coarse = error(20);
        let medium = error(40);
        let fine = error(80);
        assert!(coarse / medium > 3.95, "{coarse:e} / {medium:e}");
        assert!(medium / fine > 3.98, "{medium:e} / {fine:e}");
        assert!(fine < fine_bound, "{fine:e}");
    }
}

#[test]
fn g3_ramp_work_keeps_full_bearing_forces_separate_from_prescribed_spin_exchange() {
    let bearings = bearings();
    let rotor = rotor(&bearings);
    let source = rotor.prepare_unbalance(&[], || false).unwrap();
    let before = RotorTransientState {
        displacement: vec![0.03, -0.02, 0.04, 0.01],
        velocity: vec![-0.07, 0.05, -0.03, 0.02],
    };
    let result = rotor
        .spin_ramp_step(
            &before,
            RotorSpinState {
                time_seconds: 0.0,
                phase_rad: 0.0,
                spin_rad_s: 1.0,
            },
            &source,
            [0.2, -0.3, 0.1, -0.05],
            ramp(2.0, 0.02, RotorSpinDriveAxis::Stationary),
            || false,
        )
        .unwrap();
    let after = &result.motion.state_after;
    let q = std::array::from_fn::<_, 4, _>(|i| {
        0.5 * before.displacement[i] + 0.5 * after.displacement[i]
    });
    let v = std::array::from_fn::<_, 4, _>(|i| 0.5 * before.velocity[i] + 0.5 * after.velocity[i]);
    // Independent work at each actual bearing point, with nonsymmetric K
    // and active/cross-coupled C retained in the source coordinates.
    let mut stiffness_work = 0.0;
    let mut damping_work = 0.0;
    for bearing in bearings {
        let r = [
            q[0] + bearing.station_m * q[3],
            q[1] - bearing.station_m * q[2],
        ];
        let rate = [
            v[0] + bearing.station_m * v[3],
            v[1] - bearing.station_m * v[2],
        ];
        for i in 0..2 {
            for j in 0..2 {
                stiffness_work += 0.02 * rate[i] * bearing.stiffness_n_per_m[i][j] * r[j];
                damping_work += 0.02 * rate[i] * bearing.damping_ns_per_m[i][j] * rate[j];
            }
        }
    }
    let acceleration_work = 0.02 * 2.0 * 3.0 * (v[2] * q[3] - v[3] * q[2]);
    near(result.bearing_stiffness_work_j, stiffness_work, 3e-17);
    near(result.motion.velocity_force_work_j, damping_work, 3e-17);
    near(result.acceleration_coupling_work_j, acceleration_work, 3e-17);
    assert!(result.acceleration_coupling_work_j.abs() > 1e-7);
    near(
        result.motion.stiffness_force_work_j,
        stiffness_work + acceleration_work,
        3e-17,
    );
    near(result.motion.gyroscopic_work_j, 0.0, 3e-17);
    near(result.motion.work_energy_residual_j, 0.0, 3e-17);
}

#[test]
fn g4_invalid_ramps_and_early_or_late_cancellation_publish_no_partial_continuation() {
    let rotor = rotor(&bearings());
    let source = rotor.prepare_unbalance(&[], || false).unwrap();
    let state = zero_state();
    let state_before = state.clone();
    let spin = RotorSpinState {
        time_seconds: 0.0,
        phase_rad: 0.2,
        spin_rad_s: 1.0,
    };
    let controls = ramp(2.0, 0.01, RotorSpinDriveAxis::Stationary);
    let force = [1.0, -0.3, 0.2, -0.1];
    let calls = Cell::new(0);
    let expected = rotor
        .spin_ramp_step(&state, spin, &source, force, controls, || {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    for stop in [1, calls.get()] {
        let polls = Cell::new(0);
        assert_eq!(
            rotor.spin_ramp_step(&state, spin, &source, force, controls, || {
                polls.set(polls.get() + 1);
                polls.get() == stop
            }),
            Err(RotorError::Cancelled)
        );
        assert_eq!(state, state_before);
        assert_eq!(
            rotor
                .spin_ramp_step(&state, spin, &source, force, controls, || false)
                .unwrap(),
            expected
        );
    }
    for bad in [f64::NAN, f64::INFINITY] {
        let mut invalid = controls;
        invalid.angular_acceleration_rad_s2 = bad;
        assert!(
            rotor
                .spin_ramp_step(&state, spin, &source, force, invalid, || false)
                .is_err()
        );
        for field in 0..3 {
            let mut invalid = spin;
            match field {
                0 => invalid.time_seconds = bad,
                1 => invalid.phase_rad = bad,
                _ => invalid.spin_rad_s = bad,
            }
            assert!(
                rotor
                    .spin_ramp_step(&state, invalid, &source, force, controls, || false)
                    .is_err()
            );
        }
    }
    for duration in [0.0, -0.1, f64::from_bits(1), f64::INFINITY] {
        let mut invalid = controls;
        invalid.duration_seconds = duration;
        assert!(
            rotor
                .spin_ramp_step(&state, spin, &source, force, invalid, || false)
                .is_err()
        );
    }
    for tolerance in [0.0, 1.0, f64::NAN] {
        let mut invalid = controls;
        invalid.relative_tolerance = tolerance;
        assert!(
            rotor
                .spin_ramp_step(&state, spin, &source, force, invalid, || false)
                .is_err()
        );
    }
    let mut invalid_spin = spin;
    invalid_spin.time_seconds = 1e100;
    assert!(
        rotor
            .spin_ramp_step(&state, invalid_spin, &source, force, controls, || false)
            .is_err()
    );
    let mut invalid_state = state.clone();
    invalid_state.velocity.pop();
    assert!(
        rotor
            .spin_ramp_step(&invalid_state, spin, &source, force, controls, || false)
            .is_err()
    );
    let mut invalid = controls;
    invalid.angular_acceleration_rad_s2 = f64::MAX;
    invalid.duration_seconds = 2.0;
    assert!(
        rotor
            .spin_ramp_step(&state, spin, &source, force, invalid, || false)
            .is_err()
    );
    assert_eq!(state, state_before);
    assert_eq!(
        rotor
            .spin_ramp_step(&state, spin, &source, force, controls, || false)
            .unwrap(),
        expected
    );
}
