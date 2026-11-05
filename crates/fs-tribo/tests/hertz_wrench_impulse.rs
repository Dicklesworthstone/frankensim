//! Mechanical impacts: a finite pressure patch couples shear and torsion.

use core::f64::consts::PI;
use fs_tribo::sliding_patch::contact::HertzWrenchMode;
use fs_tribo::sliding_patch::contact::impulse::{
    HertzImpulseControls, HertzImpulseError, HertzImpulseInput, HertzImpulsePatch,
};

type Vector = [f64; 3];
type Matrix = [[f64; 3]; 3];

fn controls() -> HertzImpulseControls {
    HertzImpulseControls {
        maximum_static_iterations: 96,
        maximum_iterations: 256,
        maximum_evaluations: 2048,
        impulse_tolerance_ns: 2e-9,
        spin_impulse_tolerance_nms: 5e-10,
        velocity_tolerance_m_s: 2e-13,
        energy_tolerance_j: 2e-10,
        relative_matrix_tolerance: 1e-13,
        maximum_impulse_ns: 20.0,
        maximum_spin_impulse_nms: 5.0,
    }
}

fn mul(matrix: Matrix, x: Vector) -> Vector {
    matrix.map(|row| row.into_iter().zip(x).map(|(a, b)| a * b).sum())
}

fn input(velocity: Vector, mobility: Matrix, radius: f64, normal: f64) -> HertzImpulseInput {
    HertzImpulseInput {
        normal_impulse_ns: normal,
        velocity_before_m_s: [velocity[0], velocity[1]],
        spin_before_rad_s: velocity[2] / radius,
        scaled_mobility_per_kg: mobility,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:e}"
    );
}

#[test]
fn translation_stops_without_rebound_and_loses_the_full_initial_energy() {
    let patch = HertzImpulsePatch::new(0.2, 0.5).unwrap();
    let mobility = [[0.5, 0.0, 0.0], [0.0, 0.25, 0.0], [0.0, 0.0, 0.8]];
    let stopped = patch
        .solve(
            &input([0.2, 0.0, 0.0], mobility, 0.2, 2.0),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(stopped.mode, HertzWrenchMode::Sticking);
    close(stopped.impulse_ns[0], -0.4, 1e-15);
    close(stopped.velocity_after_m_s[0], 0.0, 1e-15);
    close(stopped.dissipated_j, 0.04, 1e-15);
    close(stopped.kinetic_energy_change_j, -0.04, 1e-15);
    assert_eq!(stopped.duality_gap_j, 0.0);

    let sliding = patch
        .solve(
            &input([2.0, 0.0, 0.0], mobility, 0.2, 2.0),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(sliding.mode, HertzWrenchMode::Sliding);
    close(sliding.impulse_ns[0], -1.0, 1e-15);
    close(sliding.velocity_after_m_s[0], 1.5, 1e-15);
    close(sliding.dissipated_j, 1.75, 1e-15);
    close(sliding.energy_residual_j, 0.0, 1e-15);
}

#[test]
fn pure_spin_uses_hertz_first_radial_moment_and_physical_angular_work() {
    let radius = 0.25;
    let patch = HertzImpulsePatch::new(radius, 0.4).unwrap();
    let mobility = [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.5]];
    let step = patch
        .solve(
            &input([0.0, 0.0, 2.0], mobility, radius, 3.0),
            controls(),
            || false,
        )
        .unwrap();
    let angular_impulse = -(3.0 * PI / 16.0) * 1.2 * radius;
    close(step.spin_impulse_nms, angular_impulse, 1e-15);
    assert_eq!(step.impulse_ns, [0.0; 2]);
    let after = (2.0 + 0.5 * angular_impulse / radius) / radius;
    close(step.spin_after_rad_s, after, 1e-15);
    close(
        step.dissipated_j,
        -angular_impulse * (8.0 + after) / 2.0,
        1e-15,
    );
    close(step.energy_residual_j, 0.0, 1e-15);
}

#[test]
fn full_mobility_arrests_a_coupled_shear_and_spin_reaction() {
    let radius = 0.2;
    let patch = HertzImpulsePatch::new(radius, 0.5).unwrap();
    let mobility = [[0.5, 0.08, -0.03], [0.08, 0.3, 0.04], [-0.03, 0.04, 0.7]];
    let reaction = [0.2, -0.1, 0.08];
    let before = mul(mobility, reaction);
    let step = patch
        .solve(&input(before, mobility, radius, 2.0), controls(), || false)
        .unwrap();
    assert_eq!(step.mode, HertzWrenchMode::Sticking);
    close(step.impulse_ns[0], -reaction[0], 2e-15);
    close(step.impulse_ns[1], -reaction[1], 2e-15);
    close(step.spin_impulse_nms, -radius * reaction[2], 2e-15);
    for velocity in step.velocity_after_m_s {
        close(velocity, 0.0, 2e-15);
    }
    close(step.spin_after_rad_s, 0.0, 2e-15);
    let initial_energy: f64 = before
        .into_iter()
        .zip(reaction)
        .map(|(v, p)| 0.5 * v * p)
        .sum();
    close(step.dissipated_j, initial_energy, 2e-15);
    close(step.energy_residual_j, 0.0, 2e-15);
}

// Independent local pressure-impulse integration. The radius substitution
// r=a*sqrt(1-u²) removes the Hertz square-root edge; dP=3N*u² du dphi/(2pi).
fn integrated_impulse(radius: f64, capacity: f64, velocity: Vector) -> Vector {
    const NODES: [f64; 4] = [
        -0.8611363115940526,
        -0.3399810435848563,
        0.3399810435848563,
        0.8611363115940526,
    ];
    const WEIGHTS: [f64; 4] = [
        0.3478548451374538,
        0.6521451548625461,
        0.6521451548625461,
        0.3478548451374538,
    ];
    let mut impulse = [0.0; 3];
    for radial in 0..16 {
        for angular in 0..16 {
            for i in 0..4 {
                let u = (radial as f64 + 0.5 + 0.5 * NODES[i]) / 16.0;
                let r = radius * (1.0 - u * u).sqrt();
                for j in 0..4 {
                    let phi = 2.0 * PI * (angular as f64 + 0.5 + 0.5 * NODES[j]) / 16.0;
                    let x = r * phi.cos();
                    let y = r * phi.sin();
                    let local = [
                        velocity[0] - velocity[2] * y / radius,
                        velocity[1] + velocity[2] * x / radius,
                    ];
                    let speed = local[0].hypot(local[1]);
                    let load =
                        3.0 * capacity * u * u * WEIGHTS[i] * WEIGHTS[j] / (4.0 * 16.0 * 16.0);
                    let px = -load * local[0] / speed;
                    let py = -load * local[1] / speed;
                    impulse[0] += px;
                    impulse[1] += py;
                    impulse[2] += (x * py - y * px) / radius;
                }
            }
        }
    }
    impulse
}

#[test]
fn mixed_impact_matches_local_pressure_integral_with_full_cross_coupling() {
    let radius = 0.2;
    let patch = HertzImpulsePatch::new(radius, 0.5).unwrap();
    let mobility = [[0.5, 0.08, -0.03], [0.08, 0.3, 0.04], [-0.03, 0.04, 0.7]];
    let endpoint = [1.2, -0.6, 0.8];
    let expected = integrated_impulse(radius, 1.0, endpoint);
    let response = mul(mobility, expected);
    let before = core::array::from_fn(|i| endpoint[i] - response[i]);
    let step = patch
        .solve(&input(before, mobility, radius, 2.0), controls(), || false)
        .unwrap();
    assert_eq!(step.mode, HertzWrenchMode::Sliding);
    for i in 0..2 {
        close(step.impulse_ns[i], expected[i], 3e-9);
        close(step.velocity_after_m_s[i], endpoint[i], 3e-9);
    }
    close(step.spin_impulse_nms, radius * expected[2], 5e-10);
    close(radius * step.spin_after_rad_s, endpoint[2], 3e-9);
    assert!(step.dissipated_j > 0.0);
    close(step.energy_residual_j, 0.0, 2e-14);
    // Neither pure-motion capacity can be independently spent while both
    // translation and spin remain nonzero.
    assert!(step.impulse_ns[0].hypot(step.impulse_ns[1]) < 1.0);
    assert!(step.spin_impulse_nms.abs() < 3.0 * PI / 16.0 * radius);
}

#[test]
fn zero_capacity_leaves_velocity_and_energy_unchanged() {
    let mobility = [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]];
    for (mu, normal) in [(0.0, 2.0), (0.5, 0.0)] {
        let patch = HertzImpulsePatch::new(0.2, mu).unwrap();
        let step = patch
            .solve(
                &input([0.3, -0.4, 0.2], mobility, 0.2, normal),
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(step.mode, HertzWrenchMode::Frictionless);
        assert_eq!(step.impulse_ns, [0.0; 2]);
        assert_eq!(step.spin_impulse_nms, 0.0);
        assert_eq!(step.velocity_after_m_s, [0.3, -0.4]);
        assert_eq!(step.spin_after_rad_s, 1.0);
        assert_eq!(step.dissipated_j, 0.0);
        assert_eq!(step.kinetic_energy_change_j, 0.0);
    }
}

#[test]
fn invalid_mobility_cancellation_and_impulse_budgets_are_explicit_refusals() {
    let patch = HertzImpulsePatch::new(0.2, 0.5).unwrap();
    let mobility = [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]];
    let valid = input([2.0, 0.0, 0.0], mobility, 0.2, 2.0);
    assert_eq!(
        patch.solve(&valid, controls(), || true),
        Err(HertzImpulseError::Cancelled)
    );
    let limited = HertzImpulseControls {
        maximum_impulse_ns: 0.1,
        ..controls()
    };
    assert!(matches!(
        patch.solve(&valid, limited, || false),
        Err(HertzImpulseError::ImpulseBudget { .. })
    ));
    for matrix in [
        [[0.0; 3]; 3],
        [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, -0.1]],
        [[0.5, 0.2, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]],
    ] {
        assert!(matches!(
            patch.solve(
                &HertzImpulseInput {
                    scaled_mobility_per_kg: matrix,
                    ..valid
                },
                controls(),
                || false
            ),
            Err(HertzImpulseError::Invalid(_))
        ));
    }
    let mut polls = 0;
    assert_eq!(
        patch.solve(&valid, controls(), || {
            polls += 1;
            polls >= 4
        }),
        Err(HertzImpulseError::Cancelled)
    );
}
