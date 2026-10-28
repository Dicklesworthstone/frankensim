//! Reference centerline inertia, modal convergence and physical base forcing.

use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};

const E: f64 = 200e9;
const A: f64 = 0.01;
const I: f64 = 8e-5;
const MU: f64 = 78.5;
const L: f64 = 3.0;

fn section() -> SectionModel {
    SectionModel::Elastic {
        ea: E * A,
        ei: E * I,
    }
}

fn cantilever(count: usize, formulation: FrameMassFormulation) -> (Frame2d, usize) {
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, formulation);
    let base = frame.add_node(0.0, 0.0);
    frame.fix(base, [true; 3]);
    let mut previous = base;
    for i in 1..=count {
        let next = frame.add_node(L * i as f64 / count as f64, 0.0);
        frame.add_element(previous, next, 5, MU, &section).unwrap();
        previous = next;
    }
    (frame, previous)
}

fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}

#[test]
fn consistent_mass_integrates_rigid_translation_and_rotation() {
    for angle in [0.0_f64, 0.61, 1.3] {
        let mut frame = Frame2d::with_mass_formulation(
            Geometry::Linear,
            FrameMassFormulation::ConsistentReference,
        );
        let base = frame.add_node(0.0, 0.0);
        let (x, y) = (L * angle.cos(), L * angle.sin());
        let tip = frame.add_node(x, y);
        frame.add_element(base, tip, 5, MU, &section).unwrap();
        frame.add_mass(tip, 17.0);
        let mass = MU * L + 17.0;
        let acceleration = [0.7, -1.2, 0.0, 0.7, -1.2, 0.0];
        let force = frame.inertial_forces(&acceleration).unwrap();
        assert!((force[0] + force[3] - mass * 0.7).abs() < 1e-12 * mass);
        assert!((force[1] + force[4] + mass * 1.2).abs() < 1e-12 * mass);
        let moment = force[2] + force[5] + x * force[4] - y * force[3];
        let exact = (MU * L / 2.0 + 17.0) * (-1.2 * x - 0.7 * y);
        assert!((moment - exact).abs() < 1e-12 * exact.abs().max(1.0));

        let angular = [0.0, 0.0, 1.0, -y, x, 1.0];
        let rotational_force = frame.inertial_forces(&angular).unwrap();
        let exact_inertia = MU * L.powi(3) / 3.0 + 17.0 * L * L;
        assert!((dot(&angular, &rotational_force) - exact_inertia).abs() < 1e-12 * exact_inertia);
    }
}

#[test]
fn consistent_mass_resolves_distributed_cantilever_modes() {
    // Roots of cos(beta)*cosh(beta) = -1 for a uniform Euler–Bernoulli cantilever.
    let roots = [1.875_104_068_711_961, 4.694_091_132_974_175];
    let scale = (E * I / (MU * L.powi(4))).sqrt();
    let exact = roots.map(|beta| beta * beta * scale);
    let lumped = cantilever(1, FrameMassFormulation::Lumped)
        .0
        .natural_frequencies()
        .unwrap();
    assert!(lumped[0] < 0.70 * exact[0]);
    let mut previous = [f64::INFINITY; 2];
    for elements in [1, 2, 4, 8] {
        let frequencies = cantilever(elements, FrameMassFormulation::ConsistentReference)
            .0
            .natural_frequencies()
            .unwrap();
        let error = [0, 1].map(|mode| (frequencies[mode] / exact[mode] - 1.0).abs());
        for mode in 0..2 {
            assert!(
                error[mode] < previous[mode],
                "{elements} elements, mode {mode}: {error:?}"
            );
        }
        if elements == 1 {
            assert!(error[0] < 0.005);
        }
        previous = error;
    }
    assert!(previous[0] < 4e-6 && previous[1] < 3e-4, "{previous:?}");
}

#[test]
fn full_mass_base_excitation_preserves_rigid_body_motion_and_work() {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    let a = frame.add_node(0.0, 0.0);
    let b = frame.add_node(L, 0.0);
    frame.add_element(a, b, 5, MU, &section).unwrap();
    let (dt, ag, direction) = (1e-4, 1.5, (0.3, -0.7));
    let history = frame
        .newmark(&[ag; 101], dt, direction, Rayleigh { a0: 0.0, a1: 0.0 })
        .unwrap();
    for (step, state) in history.u.iter().enumerate() {
        let t = step as f64 * dt;
        for node in [a, b] {
            assert!((state[3 * node] + 0.5 * ag * direction.0 * t * t).abs() < 1e-12);
            assert!((state[3 * node + 1] + 0.5 * ag * direction.1 * t * t).abs() < 1e-12);
            assert!(state[3 * node + 2].abs() < 1e-12);
        }
        let kinetic = 0.5
            * MU
            * L
            * ag
            * ag
            * t
            * t
            * (direction.0 * direction.0 + direction.1 * direction.1);
        assert!((history.kinetic_energy[step] - kinetic).abs() < 1e-9 * kinetic.max(1e-6));
    }
    assert!(history.energy_balance_error() < 1e-8);
}

#[test]
fn fixed_base_reaction_includes_all_distributed_inertia() {
    let (mut frame, tip) = cantilever(1, FrameMassFormulation::ConsistentReference);
    frame.fix(tip, [true; 3]);
    frame.add_mass(tip, 17.0);
    let ground = [0.7, -0.3, 0.9];
    let history = frame
        .newmark(&ground, 0.01, (1.0, 0.0), Rayleigh { a0: 0.0, a1: 0.0 })
        .unwrap();
    for (i, ag) in ground.iter().enumerate() {
        let exact = -(MU * L + 17.0) * ag;
        assert!((history.base_shear[i] - exact).abs() < 1e-12 * exact.abs());
    }
}

// Independent closed-form two-coordinate Hermite beam solution. This uses
// the analytic 2x2 generalized characteristic polynomial, not the frame eigensolver.
fn loaded_step_response(t: f64, ag: f64) -> f64 {
    let m = [
        MU * L * 156.0 / 420.0,
        -MU * L * 22.0 * L / 420.0,
        MU * L * 4.0 * L * L / 420.0,
    ];
    let k = [
        12.0 * E * I / L.powi(3),
        -6.0 * E * I / (L * L),
        4.0 * E * I / L,
    ];
    let determinant_m = m[0] * m[2] - m[1] * m[1];
    let determinant_k = k[0] * k[2] - k[1] * k[1];
    let cross = k[0] * m[2] + k[2] * m[0] - 2.0 * k[1] * m[1];
    let discriminant = (cross * cross - 4.0 * determinant_m * determinant_k).sqrt();
    let lambda = [
        2.0 * determinant_k / (cross + discriminant),
        (cross + discriminant) / (2.0 * determinant_m),
    ];
    let force = [-ag * MU * L / 2.0, ag * MU * L * L / 12.0];
    lambda
        .iter()
        .map(|lambda| {
            let rotation = (k[0] - lambda * m[0]) / (lambda * m[1] - k[1]);
            let modal_mass = m[0] + 2.0 * rotation * m[1] + rotation * rotation * m[2];
            let modal_force = force[0] + rotation * force[1];
            modal_force / (modal_mass * lambda) * (1.0 - (lambda.sqrt() * t).cos())
        })
        .sum()
}

#[test]
fn loaded_consistent_mass_history_matches_independent_modal_response() {
    let (mut base, tip) = cantilever(1, FrameMassFormulation::ConsistentReference);
    let mut gravity_acceleration = vec![0.0; base.ndof()];
    for node in 0..2 {
        gravity_acceleration[3 * node + 1] = -9.806_65;
    }
    let gravity = base.inertial_forces(&gravity_acceleration).unwrap();
    base.static_load(&gravity, 1).unwrap();
    let equilibrium = base.displacements()[3 * tip + 1];
    let exact_static = -MU * 9.806_65 * L.powi(4) / (8.0 * E * I);
    assert!((equilibrium - exact_static).abs() < 1e-12 * exact_static.abs());
    let omega = base.natural_frequencies().unwrap()[0];
    let duration = 0.72 * 2.0 * std::f64::consts::PI / omega;
    let exact = equilibrium + loaded_step_response(duration, 1.0);
    let mut previous = f64::INFINITY;
    for steps in [200, 400, 800] {
        let mut frame = base.clone();
        let history = frame
            .newmark_with_load(
                &vec![1.0; steps + 1],
                duration / steps as f64,
                (0.0, 1.0),
                Rayleigh { a0: 0.0, a1: 0.0 },
                &gravity,
            )
            .unwrap();
        let error = (history.u[steps][3 * tip + 1] - exact).abs();
        assert!(
            error < 0.30 * previous,
            "{steps} steps, error {error}, previous {previous}"
        );
        assert!(history.energy_balance_error() < 1e-8);
        previous = error;
    }
    assert!(previous < 1e-5 * exact.abs());
}

#[test]
fn consistent_selection_preserves_the_nodal_only_path() {
    let build = |formulation| {
        let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, formulation);
        let a = frame.add_node(0.0, 0.0);
        let b = frame.add_node(L, 0.0);
        frame.fix(a, [true; 3]);
        frame.add_element(a, b, 5, 0.0, &section).unwrap();
        frame.add_mass(b, 2000.0);
        frame
    };
    let mut lumped = build(FrameMassFormulation::Lumped);
    let mut consistent = build(FrameMassFormulation::ConsistentReference);
    assert_eq!(
        lumped.natural_frequencies().unwrap(),
        consistent.natural_frequencies().unwrap()
    );
    let ground: Vec<f64> = (0..101).map(|i| (i as f64 * 0.1).sin()).collect();
    let damping = Rayleigh { a0: 0.3, a1: 0.001 };
    assert_eq!(
        lumped.newmark(&ground, 0.001, (0.0, 1.0), damping).unwrap(),
        consistent
            .newmark(&ground, 0.001, (0.0, 1.0), damping)
            .unwrap()
    );
}
