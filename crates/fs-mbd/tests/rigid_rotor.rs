//! Physical inertia, bearing moment, gyroscopic and synchronous-unbalance oracles.
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_math::c64::C64;
use fs_mbd::articulated::SpatialInertia;
use fs_mbd::rotordynamics::RotorError;
use fs_mbd::rotordynamics::rigid::{MAX_ROTOR_STATIONS, RadialBearing, RigidRotor, RotorUnbalance};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}
fn complex(a: C64, b: C64, tolerance: f64) {
    near(a.re, b.re, tolerance);
    near(a.im, b.im, tolerance);
}
fn inertia() -> SpatialInertia {
    SpatialInertia::new(
        4.0,
        Vec3::new(0.2, -0.3, 0.1),
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0],
        },
    )
    .unwrap()
}
fn bearings(damping: f64) -> [RadialBearing; 2] {
    [-0.5, 0.5].map(|station_m| RadialBearing {
        station_m,
        stiffness_n_per_m: [[18.0, 0.0], [0.0, 18.0]],
        damping_ns_per_m: [[damping, 0.0], [0.0, damping]],
    })
}
fn rotor(damping: f64) -> RigidRotor {
    RigidRotor::new(
        inertia(),
        So3::identity(),
        &bearings(damping),
        1e-12,
        || false,
    )
    .unwrap()
}

#[test]
fn physical_inertia_and_bearing_stations_give_analytical_gyroscopic_poles() {
    let rotor = rotor(0.0);
    let model = rotor.model();
    for i in 0..4 {
        for j in 0..4 {
            near(
                model.mass()[4 * i + j],
                if i != j {
                    0.0
                } else if i < 2 {
                    4.0
                } else {
                    2.0
                },
                0.0,
            );
            near(
                model.stiffness()[4 * i + j],
                if i != j {
                    0.0
                } else if i < 2 {
                    36.0
                } else {
                    9.0
                },
                0.0,
            );
        }
    }
    near(model.gyro_per_rad_s()[11], 3.0, 0.0);
    near(model.gyro_per_rad_s()[14], -3.0, 0.0);
    for spin in [-1.7_f64, 0.0, 1.7] {
        let report = model
            .at_speed(spin, || false)
            .unwrap()
            .poles(|| false)
            .unwrap();
        let mut frequencies: Vec<_> = report
            .poles_per_s
            .iter()
            .filter_map(|p| {
                near(p.re, 0.0, 3e-10);
                (p.im > 0.0).then_some(p.im)
            })
            .collect();
        frequencies.sort_by(f64::total_cmp);
        let split = (3.0 * spin).abs();
        let root = (split * split + 4.0 * 2.0 * 9.0).sqrt();
        let mut reference = [3.0, 3.0, (root - split) / 4.0, (root + split) / 4.0];
        reference.sort_by(f64::total_cmp);
        assert_eq!(frequencies.len(), 4);
        for (actual, expected) in frequencies.into_iter().zip(reference) {
            near(actual, expected, 3e-10);
        }
    }
}

#[test]
fn circular_applied_torque_has_the_physical_forward_and_backward_gyro_sign() {
    let rotor = rotor(0.4);
    let (spin, frequency) = (1.7, 1.1);
    let point = rotor.model().at_speed(spin, || false).unwrap();
    for sense in [-1.0, 1.0] {
        // theta_y = -i theta_x is forward whirl in the stationary tilt plane.
        let force = [C64::ZERO, C64::ZERO, C64::ONE, C64::new(0.0, -sense)];
        let response = point
            .harmonic_response(frequency, &force, 1e-12, || false)
            .unwrap();
        let amplitude = C64::ONE
            / C64::new(
                9.0 - 2.0 * frequency * frequency + sense * 3.0 * spin * frequency,
                0.2 * frequency,
            );
        complex(response.displacement[2], amplitude, 2e-12);
        complex(
            response.displacement[3],
            amplitude * C64::new(0.0, -sense),
            2e-12,
        );
        complex(response.displacement[0], C64::ZERO, 2e-13);
        complex(response.displacement[1], C64::ZERO, 2e-13);
        near(response.mean_gyroscopic_power_w, 0.0, 2e-12);
    }
}

#[test]
fn synchronous_unbalance_and_recovered_bearing_work_match_jeffcott_response() {
    let rotor = rotor(0.4);
    let imbalance = RotorUnbalance {
        station_m: 0.0,
        mass_eccentricity_kg_m: 0.02,
        phase_rad: 0.37,
    };
    for spin in [-3.0_f64, -1.4, 0.0, 1.4, 3.0] {
        let force = rotor.unbalance_force(spin, &[imbalance], || false).unwrap();
        let response = rotor
            .model()
            .at_speed(spin, || false)
            .unwrap()
            .harmonic_response(spin, &force, 1e-12, || false)
            .unwrap();
        let amplitude = C64::new(imbalance.phase_rad.cos(), imbalance.phase_rad.sin())
            .scale(0.02 * spin * spin)
            / C64::new(36.0 - 4.0 * spin * spin, 0.8 * spin);
        complex(response.displacement[0], amplitude, 2e-12);
        complex(
            response.displacement[1],
            amplitude * C64::new(0.0, -1.0),
            2e-12,
        );
        complex(response.displacement[2], C64::ZERO, 2e-13);
        complex(response.displacement[3], C64::ZERO, 2e-13);
        let reactions = rotor
            .bearing_reactions(&response.displacement, &response.velocity, || false)
            .unwrap();
        near(
            reactions
                .iter()
                .map(|b| b.mean_power_on_rotor_w)
                .sum::<f64>(),
            -response.mean_input_power_w,
            2e-12,
        );
        for i in 0..2 {
            complex(
                reactions[0].moment_nm[i] + reactions[1].moment_nm[i],
                C64::ZERO,
                2e-13,
            );
        }
    }
    let offset = RotorUnbalance {
        station_m: 0.7,
        phase_rad: 0.0,
        ..imbalance
    };
    let load = rotor.unbalance_force(2.0, &[offset], || false).unwrap();
    complex(load[0], C64::from_re(0.08), 2e-14);
    complex(load[1], C64::new(0.0, -0.08), 2e-14);
    complex(load[2], C64::new(0.0, 0.056), 2e-14);
    complex(load[3], C64::from_re(0.056), 2e-14);
}

#[test]
fn nonsymmetric_bearings_preserve_station_moments_and_virtual_work() {
    let bearings = [
        RadialBearing {
            station_m: -0.4,
            stiffness_n_per_m: [[13.0, 1.2], [-0.7, 16.0]],
            damping_ns_per_m: [[0.3, 0.2], [-0.1, 0.4]],
        },
        RadialBearing {
            station_m: 0.8,
            stiffness_n_per_m: [[17.0, -0.9], [1.3, 20.0]],
            damping_ns_per_m: [[-0.1, 0.3], [0.05, 0.2]],
        },
    ];
    let rotor = RigidRotor::new(inertia(), So3::identity(), &bearings, 1e-12, || false).unwrap();
    let q = [
        C64::new(0.2, -0.1),
        C64::new(0.1, 0.3),
        C64::new(-0.07, 0.04),
        C64::new(0.08, -0.05),
    ];
    let v = [
        C64::new(-0.3, 0.2),
        C64::new(0.5, 0.1),
        C64::new(0.06, -0.09),
        C64::new(-0.02, 0.07),
    ];
    let reactions = rotor.bearing_reactions(&q, &v, || false).unwrap();
    let mut total = [C64::ZERO; 4];
    for (b, actual) in bearings.iter().zip(&reactions) {
        let u = [
            q[0] + q[3].scale(b.station_m),
            q[1] - q[2].scale(b.station_m),
        ];
        let speed = [
            v[0] + v[3].scale(b.station_m),
            v[1] - v[2].scale(b.station_m),
        ];
        let expected: [C64; 2] = std::array::from_fn(|i| {
            -(u[0].scale(b.stiffness_n_per_m[i][0])
                + u[1].scale(b.stiffness_n_per_m[i][1])
                + speed[0].scale(b.damping_ns_per_m[i][0])
                + speed[1].scale(b.damping_ns_per_m[i][1]))
        });
        for i in 0..2 {
            complex(actual.displacement_m[i], u[i], 2e-13);
            complex(actual.force_n[i], expected[i], 2e-13);
            total[i] = total[i] + actual.force_n[i];
            total[i + 2] = total[i + 2] + actual.moment_nm[i];
        }
        complex(actual.moment_nm[0], expected[1].scale(-b.station_m), 2e-13);
        complex(actual.moment_nm[1], expected[0].scale(b.station_m), 2e-13);
    }
    for i in 0..4 {
        let expected = (0..4).fold(C64::ZERO, |sum, j| {
            sum - q[j].scale(rotor.model().stiffness()[4 * i + j])
                - v[j].scale(rotor.model().damping()[4 * i + j])
        });
        complex(total[i], expected, 2e-13);
    }
    let power = 0.5 * (0..4).map(|i| (total[i] * v[i].conj()).re).sum::<f64>();
    near(
        reactions
            .iter()
            .map(|b| b.mean_power_on_rotor_w)
            .sum::<f64>(),
        power,
        2e-13,
    );
    assert_ne!(rotor.model().stiffness()[1], rotor.model().stiffness()[4]);
}

#[test]
fn original_full_body_tensor_is_lowered_in_the_declared_proper_frame() {
    let frame = So3::exp(So3Tangent::new(Vec3::new(0.5, -0.8, 0.3))).unwrap();
    let r = frame.matrix();
    let original = inertia();
    let body_tensor = r.compose(original.inertia_com()).compose(r.transpose());
    assert!(body_tensor.m[1].abs() > 0.01);
    let source =
        SpatialInertia::new(original.mass(), original.center_of_mass(), body_tensor).unwrap();
    let lowered = RigidRotor::new(source, frame, &bearings(0.4), 1e-12, || false).unwrap();
    let reference = rotor(0.4);
    assert_eq!(lowered.inertia(), source);
    assert_eq!(lowered.body_from_rotor(), frame);
    assert!(lowered.axisymmetry_relative_residual() < 1e-14);
    for (a, b) in lowered.model().mass().iter().zip(reference.model().mass()) {
        near(*a, *b, 2e-14);
    }
    for (a, b) in lowered
        .model()
        .gyro_per_rad_s()
        .iter()
        .zip(reference.model().gyro_per_rad_s())
    {
        near(*a, *b, 2e-14);
    }
    let force = [
        C64::new(0.7, -0.3),
        C64::new(0.2, 0.4),
        C64::new(-0.1, 0.2),
        C64::new(0.4, -0.3),
    ];
    let a = lowered
        .model()
        .at_speed(1.7, || false)
        .unwrap()
        .harmonic_response(1.1, &force, 1e-12, || false)
        .unwrap();
    let b = reference
        .model()
        .at_speed(1.7, || false)
        .unwrap()
        .harmonic_response(1.1, &force, 1e-12, || false)
        .unwrap();
    for (x, y) in a.displacement.iter().zip(b.displacement) {
        complex(*x, y, 3e-12);
    }
    assert!(matches!(
        RigidRotor::new(source, So3::identity(), &bearings(0.4), 1e-12, || false),
        Err(RotorError::Invalid(_))
    ));
}

#[test]
fn invalid_rotor_reductions_and_cancelled_queries_return_no_partial_result() {
    let asymmetric = SpatialInertia::new(
        4.0,
        Vec3::new(0.0, 0.0, 0.0),
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 2.2, 0.0, 0.0, 0.0, 3.0],
        },
    )
    .unwrap();
    assert!(matches!(
        RigidRotor::new(asymmetric, So3::identity(), &[], 1e-12, || false),
        Err(RotorError::Invalid(_))
    ));
    assert!(matches!(
        RigidRotor::new(inertia(), So3::identity(), &[], 0.01, || false),
        Err(RotorError::Invalid(_))
    ));
    assert!(matches!(
        RigidRotor::new(
            inertia(),
            So3::identity(),
            &vec![bearings(0.0)[0]; MAX_ROTOR_STATIONS + 1],
            1e-12,
            || false
        ),
        Err(RotorError::Invalid(_))
    ));
    let rotor = rotor(0.4);
    assert!(matches!(
        rotor.unbalance_force(2.0, &[], || true),
        Err(RotorError::Cancelled)
    ));
    assert!(matches!(
        rotor.bearing_reactions(&[C64::ZERO; 4], &[C64::ZERO; 4], || true),
        Err(RotorError::Cancelled)
    ));
    assert!(matches!(
        rotor.unbalance_force(f64::NAN, &[], || false),
        Err(RotorError::Invalid(_))
    ));
    let free = RigidRotor::new(inertia(), So3::identity(), &[], 1e-12, || false).unwrap();
    assert!(matches!(
        free.model()
            .at_speed(0.0, || false)
            .unwrap()
            .harmonic_response(0.0, &[C64::ONE; 4], 1e-12, || false),
        Err(RotorError::Linear(_))
    ));
}
