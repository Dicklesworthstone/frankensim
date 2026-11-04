//! G1/G3/G4: acceleration-driven physical unbalance on the original rotor.
use fs_ga::{Mat3, So3, Vec3};
use fs_mbd::articulated::SpatialInertia;
use fs_mbd::rotordynamics::RotorError;
use fs_mbd::rotordynamics::rigid::{MAX_ROTOR_STATIONS, RigidRotor, RotorUnbalance};
use fs_mbd::rotordynamics::transient::RotorTransientState;
use std::cell::Cell;

fn rotor() -> RigidRotor {
    let inertia = SpatialInertia::new(
        4.0,
        Vec3::default(),
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0],
        },
    )
    .unwrap();
    RigidRotor::new(inertia, So3::identity(), &[], 1e-12, || false).unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_startup_and_braking_have_tangential_force_and_the_original_bending_moment() {
    let rotor = rotor();
    let source = rotor
        .prepare_unbalance(
            &[RotorUnbalance {
                station_m: 0.25,
                mass_eccentricity_kg_m: 0.03,
                phase_rad: 0.0,
            }],
            || false,
        )
        .unwrap();
    // At rest the eccentric point accelerates +y; its reaction is -y.
    let initial = source.force_at(0.0, 0.0, 4.0, || false).unwrap();
    for (actual, expected) in initial.into_iter().zip([0.0, -0.12, 0.03, 0.0]) {
        near(actual, expected, 2e-17);
    }
    for speed in [-3.0, 3.0] {
        let actual = source.force_at(0.0, speed, -2.0, || false).unwrap();
        for (actual, expected) in actual.into_iter().zip([0.27, 0.06, -0.015, 0.0675]) {
            near(actual, expected, 2e-16);
        }
    }
    let quarter = source
        .force_at(core::f64::consts::FRAC_PI_2, 0.0, 4.0, || false)
        .unwrap();
    for (actual, expected) in quarter.into_iter().zip([0.12, 0.0, 0.0, 0.03]) {
        near(actual, expected, 2e-16);
    }
}

fn stations() -> [RotorUnbalance; 3] {
    [
        RotorUnbalance {
            station_m: -0.7,
            mass_eccentricity_kg_m: 0.03,
            phase_rad: 0.2,
        },
        RotorUnbalance {
            station_m: 0.1,
            mass_eccentricity_kg_m: 0.01,
            phase_rad: -0.4,
        },
        RotorUnbalance {
            station_m: 0.9,
            mass_eccentricity_kg_m: 0.02,
            phase_rad: 0.8,
        },
    ]
}

#[test]
fn g3_prepared_stations_retain_steady_phasor_and_phase_origin_equivalence() {
    let rotor = rotor();
    let points = stations();
    let source = rotor.prepare_unbalance(&points, || false).unwrap();
    for speed in [-7.0, 0.0, 3.0] {
        let harmonic = rotor.unbalance_force(speed, &points, || false).unwrap();
        for phase in [-1.3_f64, 0.0, 2.1] {
            let load = source.force_at(phase, speed, 0.0, || false).unwrap();
            let (s, c) = phase.sin_cos();
            for (actual, phasor) in load.into_iter().zip(harmonic) {
                near(actual, phasor.re * c - phasor.im * s, 3e-15);
            }
        }
    }
    let shift = 0.37;
    let shifted = points.map(|mut p| {
        p.phase_rad += shift;
        p
    });
    let shifted = rotor.prepare_unbalance(&shifted, || false).unwrap();
    let original = source.force_at(1.2, -2.0, 3.0, || false).unwrap();
    let moved = shifted.force_at(1.2 - shift, -2.0, 3.0, || false).unwrap();
    for (a, b) in original.into_iter().zip(moved) {
        near(a, b, 3e-15);
    }
}

#[test]
fn g3_phase_speed_and_acceleration_pullbacks_match_physical_load_perturbations() {
    let source = rotor().prepare_unbalance(&stations(), || false).unwrap();
    let seed = [0.4, -0.7, 1.3, 0.2];
    for input in [[0.2, -3.0, 1.1], [-0.5, 0.0, 2.0], [0.7, 4.0, 0.0]] {
        let pullback = source
            .force_vjp(input[0], input[1], input[2], seed, || false)
            .unwrap();
        let exact = [pullback.phase, pullback.spin, pullback.angular_acceleration];
        let evaluate = |x: [f64; 3]| {
            source
                .force_at(x[0], x[1], x[2], || false)
                .unwrap()
                .into_iter()
                .zip(seed)
                .map(|(f, s)| f * s)
                .sum::<f64>()
        };
        for i in 0..3 {
            let epsilon = 1e-5;
            let mut plus = input;
            let mut minus = input;
            plus[i] += epsilon;
            minus[i] -= epsilon;
            near(
                exact[i],
                (evaluate(plus) - evaluate(minus)) / (2.0 * epsilon),
                2e-10,
            );
        }
        if input[1] == 0.0 {
            assert_eq!(pullback.spin, 0.0);
        }
    }
}

fn startup(steps: usize) -> RotorTransientState {
    let rotor = rotor();
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
    let duration = 0.8;
    let acceleration = 4.0;
    let h = duration / steps as f64;
    let mut state = RotorTransientState {
        displacement: vec![0.0; 4],
        velocity: vec![0.0; 4],
    };
    for k in 0..steps {
        let midpoint_time = (k as f64 + 0.5) * h;
        let phase = 0.5 * acceleration * midpoint_time * midpoint_time;
        let speed = acceleration * midpoint_time;
        let load = source
            .force_at(phase, speed, acceleration, || false)
            .unwrap();
        // Original four-DOF model; the centered unbalance excites translation,
        // so the tilt remains zero throughout this variable-speed trajectory.
        let integrator = rotor
            .model()
            .at_speed(speed, || false)
            .unwrap()
            .midpoint_stepper(h, 1e-11, || false)
            .unwrap();
        let step = integrator.step(&state, &load, || false).unwrap();
        near(step.work_energy_residual_j, 0.0, 2e-16);
        state = step.state_after;
        assert_eq!(&state.displacement[2..], &[0.0; 2]);
        assert_eq!(&state.velocity[2..], &[0.0; 2]);
    }
    state
}

#[test]
fn g1_accelerating_unbalance_drives_native_motion_with_second_order_refinement() {
    let (duration, acceleration, phase0) = (0.8, 4.0, 0.3_f64);
    let phase = phase0 + 0.5 * acceleration * duration * duration;
    let speed = acceleration * duration;
    // Integrate M*q''=-U*d²(e_r)/dt² analytically from physical rest.
    let (s0, c0) = phase0.sin_cos();
    let (s, c) = phase.sin_cos();
    let ratio = 0.02 / 4.0;
    let exact = [
        -ratio * (c - c0),
        -ratio * (s - s0),
        ratio * speed * s,
        -ratio * speed * c,
    ];
    let error = |steps| {
        let state = startup(steps);
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
    assert!(fine < 1.2e-6, "{fine:e}");
}

#[test]
fn g4_sources_are_reusable_after_refusal_and_keep_representable_high_speed_loads() {
    let rotor = rotor();
    let points = stations();
    let source = rotor.prepare_unbalance(&points, || false).unwrap();
    let before = source.force_at(0.2, 3.0, -1.0, || false).unwrap();
    for bad in [f64::NAN, f64::INFINITY] {
        assert!(source.force_at(bad, 1.0, 1.0, || false).is_err());
        assert!(source.force_at(0.0, bad, 1.0, || false).is_err());
        assert!(source.force_at(0.0, 1.0, bad, || false).is_err());
    }
    assert_eq!(
        source.force_at(0.2, 3.0, -1.0, || true),
        Err(RotorError::Cancelled)
    );
    let calls = Cell::new(0);
    assert_eq!(
        rotor.prepare_unbalance(&points, || {
            calls.set(calls.get() + 1);
            calls.get() == 3
        }),
        Err(RotorError::Cancelled)
    );
    assert!(
        rotor
            .prepare_unbalance(&[points[0]; MAX_ROTOR_STATIONS + 1], || false)
            .is_err()
    );
    let mut bad = points;
    bad[1].mass_eccentricity_kg_m = -1.0;
    assert!(rotor.prepare_unbalance(&bad, || false).is_err());
    assert_eq!(source.force_at(0.2, 3.0, -1.0, || false).unwrap(), before);
    assert_eq!(
        rotor
            .prepare_unbalance(&[], || false)
            .unwrap()
            .force_at(1.0, 4.0, 3.0, || false)
            .unwrap(),
        [0.0; 4]
    );

    let tiny = rotor
        .prepare_unbalance(
            &[RotorUnbalance {
                station_m: 0.0,
                mass_eccentricity_kg_m: f64::from_bits(1),
                phase_rad: 0.0,
            }],
            || false,
        )
        .unwrap();
    let force = tiny.force_at(0.0, 1e200, 0.0, || false).unwrap();
    near(force[0] / 4.940656458412466e76, 1.0, 3e-16);
    assert_eq!(&force[1..], &[0.0; 3]);
    assert!(tiny.force_at(0.0, 0.5, 0.0, || false).is_err());
}
