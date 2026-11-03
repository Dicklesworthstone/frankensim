//! Physical acceleration must remain resolved below absolute-displacement ULPs.

use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};

fn axial(formulation: FrameMassFormulation) -> (Frame2d, f64) {
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, formulation);
    let base = frame.add_node(0.0, 0.0);
    let tip = frame.add_node(1.0, 0.0);
    frame.fix(base, [true; 3]);
    frame.fix(tip, [false, true, true]);
    frame
        .add_element(base, tip, 3, 6.0, &|| SectionModel::Elastic {
            ea: 4000.0,
            ei: 20.0,
        })
        .unwrap();
    frame.add_mass(tip, 2.0);
    frame
        .static_load(&[0.0, 0.0, 0.0, 4.0, 0.0, 0.0], 1)
        .unwrap();
    let mass = match formulation {
        FrameMassFormulation::Lumped => 5.0,
        FrameMassFormulation::ConsistentReference => 4.0,
    };
    (frame, mass)
}

#[test]
fn g1_tiny_steps_retain_elastic_acceleration_and_force_work() {
    for formulation in [
        FrameMassFormulation::Lumped,
        FrameMassFormulation::ConsistentReference,
    ] {
        let (mut frame, mass) = axial(formulation);
        let initial = frame.displacements()[3];
        let acceleration = -4000.0 * initial / mass;
        let h = 1e-12;
        let steps = 16;
        let history = frame
            .newmark(
                &vec![0.0; steps + 1],
                h,
                (0.0, 0.0),
                Rayleigh { a0: 0.0, a1: 0.0 },
            )
            .expect("resolved acceleration cannot depend on absolute-position precision");
        // Every displacement update is below an ULP, but velocity, momentum,
        // kinetic energy and the balancing internal work remain representable.
        assert!(history.u.iter().all(|u| u[3] == initial));
        let velocity = steps as f64 * h * acceleration;
        let kinetic = 0.5 * mass * velocity * velocity;
        assert!((frame.accelerations()[3] / acceleration - 1.0).abs() < 2e-14);
        assert!((frame.velocities()[3] / velocity - 1.0).abs() < 2e-14);
        assert!((history.kinetic_energy[steps] / kinetic - 1.0).abs() < 2e-14);
        assert!((history.internal_work[steps] / kinetic + 1.0).abs() < 2e-14);
    }
}

#[test]
fn g1_changed_base_load_resolves_acceleration_with_physical_damping() {
    for formulation in [
        FrameMassFormulation::Lumped,
        FrameMassFormulation::ConsistentReference,
    ] {
        let (mut frame, mass) = axial(formulation);
        let initial = frame.displacements()[3];
        let restoring = 4000.0 * initial;
        let a0 = -restoring / mass;
        let h = 1e-12;
        let ground = 1.6;
        let damping = Rayleigh { a0: 0.4, a1: 1e-4 };
        let c = damping.a0 * mass + damping.a1 * 4000.0;
        // The full consistent row sums to 5 kg, including the support-side
        // coupling; its free diagonal is 4 kg. Base forcing uses the full row.
        let expected_a = (-5.0 * ground - restoring - c * 0.5 * h * a0) / (mass + 0.5 * h * c);
        let expected_v = 0.5 * h * (a0 + expected_a);
        frame
            .newmark(&[0.0, ground], h, (1.0, 0.0), damping)
            .unwrap();
        assert_eq!(frame.displacements()[3], initial);
        assert!((frame.accelerations()[3] / expected_a - 1.0).abs() < 2e-14);
        assert!((frame.velocities()[3] / expected_v - 1.0).abs() < 2e-14);
    }
}

#[test]
fn g1_implicit_endpoint_survives_an_inadmissible_explicit_predictor() {
    let mut frame = Frame2d::new(Geometry::Corotational);
    let base = frame.add_node(0.0, 0.0);
    let tip = frame.add_node(1.0, 0.0);
    frame.fix(base, [true; 3]);
    frame.fix(tip, [false, true, true]);
    frame.add_mass(tip, 1.0);
    frame
        .add_element(base, tip, 3, 0.0, &|| SectionModel::Elastic {
            ea: 1000.0,
            ei: 1.0,
        })
        .unwrap();
    // Initial a=-2 m/s² predicts u=-1 m and a collapsed chord after one
    // second. The implicit average-acceleration endpoint is only -1/251 m.
    frame
        .newmark_with_load(
            &[0.0; 2],
            1.0,
            (0.0, 0.0),
            Rayleigh { a0: 0.0, a1: 0.0 },
            &[0.0, 0.0, 0.0, -2.0, 0.0, 0.0],
        )
        .unwrap();
    assert!((frame.displacements()[3] + 1.0 / 251.0).abs() < 1e-12);
    assert!((frame.accelerations()[3] + 2.0 + 1000.0 * frame.displacements()[3]).abs() < 1e-8);
}

#[test]
fn g0_new_support_freezes_moving_coordinate_and_clears_static_support_derivatives() {
    for formulation in [
        FrameMassFormulation::Lumped,
        FrameMassFormulation::ConsistentReference,
    ] {
        let (mut frame, _) = axial(formulation);
        frame
            .newmark(&[0.0; 3], 0.01, (0.0, 0.0), Rayleigh { a0: 0.0, a1: 0.0 })
            .unwrap();
        assert!(frame.velocities()[3].abs() > 1e-6);
        assert!(frame.accelerations()[3].abs() > 1e-6);
        let frozen = frame.displacements().to_vec();

        frame.fix(1, [true, false, false]);
        assert_eq!(frame.displacements(), frozen);
        assert_eq!(frame.velocities()[3], 0.0);
        assert_eq!(frame.accelerations()[3], 0.0);
        // False flags must not release an already-established support.
        frame.fix(1, [false; 3]);
        let held = frame
            .newmark(
                &[0.2, -0.4, 0.7],
                0.01,
                (1.0, 0.0),
                Rayleigh { a0: 0.4, a1: 1e-4 },
            )
            .unwrap();
        assert!(held.u.iter().all(|u| u == &frozen));
        assert!(frame.velocities().iter().all(|v| *v == 0.0));
        assert!(frame.accelerations().iter().all(|a| *a == 0.0));
        assert!(held.kinetic_energy.iter().all(|e| *e == 0.0));
        assert!(held.damping_energy.iter().all(|e| *e == 0.0));
        assert!(held.input_energy.iter().all(|e| *e == 0.0));
    }
}
