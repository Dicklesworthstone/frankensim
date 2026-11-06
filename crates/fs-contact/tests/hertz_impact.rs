//! Coupled normal restitution against independent finite-patch mechanics.

use core::f64::consts::PI;
use fs_contact::tangential::impact::hertz::{
    HertzImpactControls, HertzImpactError, HertzImpactInput, HertzImpactReport, solve_hertz_impact,
};
use fs_solver::orthant::active_set::SpdOrthantConfig;
use fs_tribo::sliding_patch::contact::HertzWrenchMode;
use fs_tribo::sliding_patch::contact::impulse::{HertzImpulseControls, HertzImpulsePatch};

type Vector = [f64; 4];
type Matrix = [[f64; 4]; 4];

fn controls() -> HertzImpactControls {
    HertzImpactControls {
        patch: HertzImpulseControls {
            maximum_static_iterations: 96,
            maximum_iterations: 256,
            maximum_evaluations: 4096,
            impulse_tolerance_ns: 2e-9,
            spin_impulse_tolerance_nms: 5e-10,
            velocity_tolerance_m_s: 2e-12,
            energy_tolerance_j: 3e-10,
            relative_matrix_tolerance: 1e-13,
            maximum_impulse_ns: 100.0,
            maximum_spin_impulse_nms: 100.0,
        },
        matrix: SpdOrthantConfig {
            max_variables: 4,
            symmetry_tolerance: 1e-13,
            relative_pivot_tolerance: 1e-13,
            ..SpdOrthantConfig::default()
        },
        maximum_iterations: 96,
        maximum_evaluations: 100_000,
        normal_impulse_tolerance_ns: 2e-9,
        normal_velocity_tolerance_m_s: 2e-11,
        energy_tolerance_j: 5e-10,
        maximum_normal_impulse_ns: 100.0,
        minimum_determinacy_margin: 1e-8,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}"
    );
}
fn vector(actual: Vector, expected: Vector, tolerance: f64) {
    for i in 0..4 {
        near(actual[i], expected[i], tolerance);
    }
}
fn multiply(w: Matrix, x: Vector) -> Vector {
    w.map(|row| row.into_iter().zip(x).map(|(a, b)| a * b).sum())
}
fn input(velocity: Vector, mobility: Matrix, restitution: f64) -> HertzImpactInput {
    HertzImpactInput {
        velocity_before_m_s: velocity,
        scaled_mobility_per_kg: mobility,
        restitution_coefficient: restitution,
    }
}
fn diagonal() -> Matrix {
    [
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 0.25, 0.0, 0.0],
        [0.0, 0.0, 0.4, 0.0],
        [0.0, 0.0, 0.0, 0.8],
    ]
}
fn coupled() -> Matrix {
    [
        [0.8, 0.12, -0.05, 0.08],
        [0.12, 0.6, 0.04, -0.02],
        [-0.05, 0.04, 0.7, 0.03],
        [0.08, -0.02, 0.03, 0.5],
    ]
}

// Independent partial-pivot elimination for the original kinetic metric.
// Production uses retained equilibrated Cholesky and polarized momentum work.
fn kinetic(w: Matrix, velocity: Vector) -> f64 {
    let mut a = [[0.0; 5]; 4];
    for i in 0..4 {
        a[i][..4].copy_from_slice(&w[i]);
        a[i][4] = velocity[i];
    }
    for column in 0..4 {
        let pivot = (column..4)
            .max_by(|&i, &j| a[i][column].abs().total_cmp(&a[j][column].abs()))
            .unwrap();
        a.swap(column, pivot);
        let scale = a[column][column];
        assert!(scale.abs() > 1e-8);
        for j in column..5 {
            a[column][j] /= scale;
        }
        for i in 0..4 {
            if i != column {
                let scale = a[i][column];
                for j in column..5 {
                    a[i][j] -= scale * a[column][j];
                }
            }
        }
    }
    (0..4).map(|i| 0.5 * velocity[i] * a[i][4]).sum()
}
fn audit(input: HertzImpactInput, out: HertzImpactReport) {
    let response = multiply(input.scaled_mobility_per_kg, out.scaled_impulse_ns);
    let expected = core::array::from_fn(|i| input.velocity_before_m_s[i] + response[i]);
    vector(out.velocity_after_m_s, expected, 1e-12);
    let difference = kinetic(input.scaled_mobility_per_kg, out.velocity_after_m_s)
        - kinetic(input.scaled_mobility_per_kg, input.velocity_before_m_s);
    near(difference, out.kinetic_energy_change_j, 2e-12);
    near(difference, -out.dissipated_j, 2e-12);
    near(difference, out.normal_work_j + out.friction_work_j, 2e-12);
    near(out.energy_residual_j, 0.0, 2e-12);
    assert!(out.dissipated_j >= -controls().energy_tolerance_j);
    assert!(out.normal_velocity_residual_m_s.abs() <= controls().normal_velocity_tolerance_m_s);
    assert!(out.normal_complementarity_j.abs() <= controls().energy_tolerance_j);
    assert!(out.determinacy.normal_growth_lower_bound_per_kg > 0.0);
    assert!(out.determinacy.friction_feedback_upper_bound < 1.0);
}

#[test]
fn g1_uncoupled_translation_and_spin_match_normal_and_hertz_closed_forms() {
    let patch = HertzImpulsePatch::new(0.2, 0.4).unwrap();
    let data = input([-2.0, 4.0, 0.0, 0.0], diagonal(), 0.5);
    let out = solve_hertz_impact(&patch, &data, controls(), || false).unwrap();
    vector(out.scaled_impulse_ns, [6.0, -2.4, 0.0, 0.0], 2e-15);
    vector(out.velocity_after_m_s, [1.0, 3.4, 0.0, 0.0], 2e-15);
    near(out.normal_work_j, -3.0, 2e-15);
    near(out.friction_work_j, -8.88, 3e-15);
    near(out.dissipated_j, 11.88, 3e-15);
    assert_eq!(out.patch_solves, 1);
    assert_eq!(out.patch.mode, HertzWrenchMode::Sliding);
    audit(data, out);

    let data = input([-1.0, 0.0, 0.0, 3.0], diagonal(), 0.25);
    let out = solve_hertz_impact(&patch, &data, controls(), || false).unwrap();
    let spin = -3.0 * PI / 16.0;
    vector(out.scaled_impulse_ns, [2.5, 0.0, 0.0, spin], 2e-15);
    near(out.patch.spin_impulse_nms, 0.2 * spin, 2e-15);
    near(out.velocity_after_m_s[3], 3.0 + 0.8 * spin, 2e-15);
    near(out.normal_work_j, -0.9375, 2e-15);
    audit(data, out);
}

#[test]
fn g1_frictionless_restitution_opening_and_tiny_closing_motion_are_exact_limits() {
    let patch = HertzImpulsePatch::new(0.2, 0.0).unwrap();
    let w = coupled();
    for e in [0.0, 0.5, 1.0] {
        let data = input([-2.0, 0.7, -0.2, 1.0], w, e);
        let out = solve_hertz_impact(&patch, &data, controls(), || false).unwrap();
        let normal = 2.0 * (1.0 + e) / w[0][0];
        vector(out.scaled_impulse_ns, [normal, 0.0, 0.0, 0.0], 2e-15);
        near(out.velocity_after_m_s[0], 2.0 * e, 2e-15);
        near(out.dissipated_j, 2.0 * (1.0 - e * e) / w[0][0], 3e-15);
        audit(data, out);
    }
    let patch = HertzImpulsePatch::new(0.2, 0.4).unwrap();
    let data = input([0.3, -2.0, 0.4, 1.0], w, 0.7);
    let out = solve_hertz_impact(
        &patch,
        &data,
        HertzImpactControls {
            maximum_normal_impulse_ns: 0.0,
            ..controls()
        },
        || false,
    )
    .unwrap();
    assert_eq!(out.scaled_impulse_ns, [0.0; 4]);
    assert_eq!(out.velocity_after_m_s, data.velocity_before_m_s);
    assert_eq!(out.dissipated_j, 0.0);
    assert_eq!(out.normal_velocity_residual_m_s, 0.0);

    let patch = HertzImpulsePatch::new(0.2, 0.0).unwrap();
    let data = input([-1e-20, 0.0, 0.0, 0.0], diagonal(), 0.0);
    let out = solve_hertz_impact(&patch, &data, controls(), || false).unwrap();
    assert_eq!(out.scaled_impulse_ns[0], 2e-20);
    assert_eq!(out.velocity_after_m_s[0], 0.0);
    assert!(out.dissipated_j > 0.0);
}

// Independent local Hertz pressure integration. r/a=sqrt(1-u²) gives
// dJn=3*Jn*u² du dphi/(2*pi), regular at the edge. The result is radius-scaled.
fn pressure_integral(capacity: f64, velocity: [f64; 3]) -> [f64; 3] {
    const X: [f64; 4] = [
        -0.8611363115940526,
        -0.3399810435848563,
        0.3399810435848563,
        0.8611363115940526,
    ];
    const Q: [f64; 4] = [
        0.3478548451374538,
        0.6521451548625461,
        0.6521451548625461,
        0.3478548451374538,
    ];
    let mut out = [0.0; 3];
    for radial in 0..16 {
        for angular in 0..16 {
            for i in 0..4 {
                let u = (radial as f64 + 0.5 + 0.5 * X[i]) / 16.0;
                let radius = (1.0 - u * u).sqrt();
                for j in 0..4 {
                    let phi = 2.0 * PI * (angular as f64 + 0.5 + 0.5 * X[j]) / 16.0;
                    let x = radius * phi.cos();
                    let y = radius * phi.sin();
                    let local = [velocity[0] - velocity[2] * y, velocity[1] + velocity[2] * x];
                    let speed = local[0].hypot(local[1]);
                    let load = 3.0 * capacity * u * u * Q[i] * Q[j] / (4.0 * 16.0 * 16.0);
                    let px = -load * local[0] / speed;
                    let py = -load * local[1] / speed;
                    out[0] += px;
                    out[1] += py;
                    out[2] += x * py - y * px;
                }
            }
        }
    }
    out
}

fn manufactured() -> (HertzImpulsePatch, HertzImpactInput, Vector, Vector) {
    let w = coupled();
    let tangent = [1.2, -0.6, 0.8];
    let q = pressure_integral(1.0, tangent);
    let impulse = [2.0, q[0], q[1], q[2]];
    let response = multiply(w, impulse);
    let e = 0.35;
    let initial_normal = -response[0] / (1.0 + e);
    let endpoint = [-e * initial_normal, tangent[0], tangent[1], tangent[2]];
    let before = core::array::from_fn(|i| endpoint[i] - response[i]);
    (
        HertzImpulsePatch::new(0.2, 0.5).unwrap(),
        input(before, w, e),
        impulse,
        endpoint,
    )
}

#[test]
fn g1_full_normal_shear_spin_coupling_matches_independent_hertz_pressure() {
    let (patch, data, impulse, endpoint) = manufactured();
    let out = solve_hertz_impact(&patch, &data, controls(), || false).unwrap();
    vector(out.scaled_impulse_ns, impulse, 2e-8);
    vector(out.velocity_after_m_s, endpoint, 2e-8);
    assert_eq!(out.patch.mode, HertzWrenchMode::Sliding);
    assert!(out.patch_solves > 1);
    assert!(out.determinacy.friction_feedback_upper_bound > 0.0);
    assert!(out.determinacy.schur.is_some());
    assert!(out.determinacy.inverse_schur_residual < 1e-12);
    // A normal-first solve would miss the returned normal correction.
    let independent_normal = -(1.0 + data.restitution_coefficient) * data.velocity_before_m_s[0]
        / data.scaled_mobility_per_kg[0][0];
    assert!((independent_normal - out.scaled_impulse_ns[0]).abs() > 0.1);
    // Conditional patch loss is not the simultaneous total event loss.
    assert!((out.dissipated_j - out.patch.dissipated_j).abs() > 0.1);
    audit(data, out);
}

#[test]
fn g1_full_arrest_uses_the_whole_four_coordinate_metric() {
    let w = coupled();
    let impulse = [1.5, 0.07, -0.03, 0.02];
    let before = multiply(w, impulse).map(|v| -v);
    let data = input(before, w, 0.0);
    let patch = HertzImpulsePatch::new(0.2, 0.6).unwrap();
    let out = solve_hertz_impact(&patch, &data, controls(), || false).unwrap();
    vector(out.scaled_impulse_ns, impulse, 4e-11);
    vector(out.velocity_after_m_s, [0.0; 4], 4e-11);
    assert_eq!(out.patch.mode, HertzWrenchMode::Sticking);
    near(out.dissipated_j, kinetic(w, before), 2e-12);
    audit(data, out);
}

fn transform_vector(q: Matrix, v: Vector) -> Vector {
    multiply(q, v)
}
fn transform_matrix(q: Matrix, w: Matrix) -> Matrix {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (0..4)
                .flat_map(|k| (0..4).map(move |l| q[i][k] * w[k][l] * q[j][l]))
                .sum()
        })
    })
}

#[test]
fn g3_tangent_rotation_and_reflection_transport_spin_with_the_pressure_patch() {
    let (patch, data, _, _) = manufactured();
    let reference = solve_hertz_impact(&patch, &data, controls(), || false).unwrap();
    let (s, c) = 0.39_f64.sin_cos();
    for transform in [
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, c, -s, 0.0],
            [0.0, s, c, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, -1.0],
        ],
    ] {
        let changed = input(
            transform_vector(transform, data.velocity_before_m_s),
            transform_matrix(transform, data.scaled_mobility_per_kg),
            data.restitution_coefficient,
        );
        let out = solve_hertz_impact(&patch, &changed, controls(), || false).unwrap();
        vector(
            out.scaled_impulse_ns,
            transform_vector(transform, reference.scaled_impulse_ns),
            2e-8,
        );
        vector(
            out.velocity_after_m_s,
            transform_vector(transform, reference.velocity_after_m_s),
            2e-8,
        );
        near(out.dissipated_j, reference.dissipated_j, 2e-9);
        audit(changed, out);
    }
}

#[test]
fn g1_newton_restitution_that_creates_global_energy_is_refused_without_clipping() {
    let w = [
        [1.0, 0.8, 0.0, 0.0],
        [0.8, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    let data = input([-1.0, -1.0, 0.0, 0.0], w, 1.0);
    let patch = HertzImpulsePatch::new(0.2, 0.3).unwrap();
    // Sliding root: Jn=2/(1-.8*.3), Jt=-.3*Jn. Restitution does zero
    // normal midpoint work, but Jt*(vt_before+vt_after)/2 is positive.
    let normal = 2.0 / (1.0 - 0.8 * 0.3);
    let tangent = -0.3 * normal;
    let after_tangent = -1.0 + 0.8 * normal + tangent;
    let created = tangent * (-1.0 + after_tangent) / 2.0;
    assert!(created > 0.2);
    match solve_hertz_impact(&patch, &data, controls(), || false) {
        Err(HertzImpactError::EnergyBalance {
            dissipated_j,
            kinetic_energy_change_j,
            residual_j,
        }) => {
            near(dissipated_j, -created, 2e-10);
            near(kinetic_energy_change_j, created, 2e-10);
            near(residual_j, 0.0, 1e-12);
        }
        other => panic!("expected original global-energy refusal, got {other:?}"),
    }
}

#[test]
fn g1_singular_or_unproven_feedback_regions_do_not_receive_a_selected_root() {
    let w = [
        [1.0, 0.8, 0.0, 0.0],
        [0.8, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    let data = input([-1.0, 2.0, 0.0, 0.0], w, 0.2);
    for mu in [0.9, 2.0] {
        assert!(matches!(
            solve_hertz_impact(
                &HertzImpulsePatch::new(0.2, mu).unwrap(),
                &data,
                controls(),
                || false
            ),
            Err(HertzImpactError::Determinacy { .. })
        ));
    }
    let patch = HertzImpulsePatch::new(0.2, 0.3).unwrap();
    let mut invalid = data;
    invalid.scaled_mobility_per_kg[0][1] += 0.1;
    assert!(matches!(
        solve_hertz_impact(&patch, &invalid, controls(), || false),
        Err(HertzImpactError::Mobility(_))
    ));
    invalid = data;
    invalid.scaled_mobility_per_kg[3][3] = 0.0;
    assert!(matches!(
        solve_hertz_impact(&patch, &invalid, controls(), || false),
        Err(HertzImpactError::Mobility(_))
    ));
    invalid = data;
    invalid.restitution_coefficient = 1.01;
    assert!(matches!(
        solve_hertz_impact(&patch, &invalid, controls(), || false),
        Err(HertzImpactError::Invalid(_))
    ));
}

#[test]
fn g4_normal_patch_and_global_work_budgets_and_cancellation_refuse_atomically() {
    let (patch, data, _, _) = manufactured();
    let unchanged = data;
    assert!(matches!(
        solve_hertz_impact(
            &patch,
            &data,
            HertzImpactControls {
                maximum_normal_impulse_ns: 0.1,
                ..controls()
            },
            || false
        ),
        Err(HertzImpactError::NormalImpulseBudget { .. })
    ));
    let mut limited = controls();
    limited.patch.maximum_impulse_ns = 0.01;
    assert!(matches!(
        solve_hertz_impact(&patch, &data, limited, || false),
        Err(HertzImpactError::Patch { .. })
    ));
    limited = controls();
    limited.maximum_evaluations = 1;
    assert!(matches!(
        solve_hertz_impact(&patch, &data, limited, || false),
        Err(HertzImpactError::Patch { .. } | HertzImpactError::EvaluationBudget { .. })
    ));
    let mut polls = 0;
    let reference = solve_hertz_impact(&patch, &data, controls(), || {
        polls += 1;
        false
    })
    .unwrap();
    assert!(polls > 4);
    for stop in [0, 1, 3, polls / 2, polls - 1] {
        let mut seen = 0;
        assert_eq!(
            solve_hertz_impact(&patch, &data, controls(), || {
                let cancel = seen == stop;
                seen += 1;
                cancel
            }),
            Err(HertzImpactError::Cancelled)
        );
        assert_eq!(data, unchanged);
    }
    assert_eq!(
        solve_hertz_impact(&patch, &data, controls(), || false).unwrap(),
        reference
    );
}
