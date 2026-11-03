//! Reciprocal finite-body impacts, contact coupling, original support work,
//! independent world momentum and atomic multi-owner refusal/continuation.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactError, NormalVelocityLaw,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::contact_network::{
    SolidContactBody, SolidContactBodyImpulse, SolidContactNetworkControls,
    SolidContactNetworkError, SolidContactOther, SolidContactPair, SolidContactPoint,
    apply_frictional_impacts,
};
use fs_solid::dynamics::frictional_impact::{
    SolidFrictionalImpactControls, SolidFrictionalImpactError,
};
use fs_solid::dynamics::spectral::{CompleteSpectralSolid, CompleteSpectralSolidControls};
use fs_solid::dynamics::supported::{SupportedSolidControls, SupportedSpectralSolid};
use fs_solid::dynamics::{ExactLinearSolid, ExactLinearSolidControls, LinearSolidImpulseError};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::Coo;

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 62,
                kernel_id: 45,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn dense_controls() -> ExactLinearSolidControls {
    ExactLinearSolidControls {
        maximum_free_dofs: 24,
        maximum_preparation_bytes: 32 * 1024 * 1024,
        maximum_exponential_multiply_terms: 1_000_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_work_error_j: 1e-8,
        maximum_steps: 100,
    }
}

fn spectral_controls() -> CompleteSpectralSolidControls {
    CompleteSpectralSolidControls {
        maximum_free_dofs: 24,
        maximum_preparation_bytes: 2 * 1024 * 1024,
        maximum_preparation_multiply_terms: 2_000_000,
        maximum_step_multiply_terms: 2_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_absolute_eigen_residual: 0.0,
        maximum_relative_state_error: 1e-12,
        maximum_work_error_j: 1e-8,
        maximum_steps: 100,
    }
}

fn controls() -> SolidContactNetworkControls {
    SolidContactNetworkControls {
        impact: SolidFrictionalImpactControls {
            contact: FrictionalContactConfig::default(),
            maximum_impulse_error_kg_m_s: 1e-10,
            maximum_contact_speed_error_mps: 1e-10,
            maximum_work_error_j: 1e-9,
            maximum_workspace_bytes: 1024 * 1024,
            maximum_response_multiply_terms: 1_000_000,
        },
        maximum_bodies: 8,
        maximum_total_dofs: 64,
        maximum_workspace_bytes: 1024 * 1024,
    }
}

fn assembly(m: &[f64], n: usize) -> TetElasticAssembly {
    let mut mass = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            mass.push(i, j, m[i * n + j]);
        }
    }
    TetElasticAssembly {
        stiffness: Coo::new(n, n).assemble(),
        mass: mass.assemble(),
        free_dofs: (0..n).collect(),
        total_mass_kg: 1.0,
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0],
    }
}

fn diagonal(value: f64) -> [f64; 9] {
    [value, 0.0, 0.0, 0.0, value, 0.0, 0.0, 0.0, value]
}

fn exact(mass: f64, velocity: [f64; 3], cx: &Cx<'_>) -> ExactLinearSolid {
    ExactLinearSolid::prepare(
        &assembly(&diagonal(mass), 3),
        0.01,
        &[0.0; 3],
        &velocity,
        dense_controls(),
        cx,
    )
    .unwrap()
}

fn spectral(mass: f64, velocity: [f64; 3], cx: &Cx<'_>) -> CompleteSpectralSolid {
    CompleteSpectralSolid::prepare_with_basis(
        &assembly(&diagonal(mass), 3),
        &[0.0; 3],
        &diagonal(1.0 / mass.sqrt()),
        0.01,
        &[0.0; 3],
        &velocity,
        spectral_controls(),
        cx,
    )
    .unwrap()
}

fn point(body: usize) -> SolidContactPoint<'static> {
    SolidContactPoint {
        body,
        velocity_weights: [&IDENTITY[0], &IDENTITY[1], &IDENTITY[2]],
    }
}

fn pair(
    id: u64,
    first: usize,
    second: usize,
    restitution: f64,
    mu: f64,
) -> SolidContactPair<'static> {
    SolidContactPair {
        id,
        first: point(first),
        second: SolidContactOther::Body(point(second)),
        gap_m: 0.0,
        activation_tolerance_m: 0.0,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution {
            coefficient: restitution,
        },
    }
}

fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "actual={actual:e}, expected={expected:e}, allowance={allowance:e}"
    );
}

fn near_slice(actual: &[f64], expected: &[f64], allowance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &e) in actual.iter().zip(expected) {
        near(a, e, allowance);
    }
}

fn exact_bits(body: &ExactLinearSolid) -> Vec<u64> {
    body.displacement_m()
        .iter()
        .chain(body.velocity_mps())
        .copied()
        .chain([
            body.time_s(),
            body.stored_energy_j(),
            body.kinetic_energy_j(),
            body.cumulative_work_j(),
            body.cumulative_dissipated_energy_j(),
            body.accepted_steps() as f64,
        ])
        .map(f64::to_bits)
        .collect()
}

fn spectral_bits(body: &CompleteSpectralSolid) -> Vec<u64> {
    body.displacement_m()
        .iter()
        .chain(body.velocity_mps())
        .copied()
        .chain([
            body.time_s(),
            body.stored_energy_j(),
            body.kinetic_energy_j(),
            body.cumulative_work_j(),
            body.steps() as f64,
        ])
        .map(f64::to_bits)
        .collect()
}

#[test]
fn g1_two_finite_owners_exchange_momentum_work_and_continue_original_motion() {
    with_cx(false, |cx| {
        let mut a = exact(2.0, [-1.0, 2.0, 0.4], cx);
        let mut b = spectral(3.0, [0.2, -0.3, 0.4], cx);
        let report = apply_frictional_impacts(
            &mut [
                SolidContactBody::Exact(&mut a),
                SolidContactBody::Spectral(&mut b),
            ],
            &[pair(2, 0, 1, 0.5, 0.25)],
            controls(),
            cx,
        )
        .unwrap();
        // W=(1/2+1/3)I. Newton normal p=2.16; sliding p_t=-.54.
        near_slice(&report.contact.law.impulses[0], &[2.16, -0.54, 0.0], 2e-9);
        near_slice(a.velocity_mps(), &[0.08, 1.73, 0.4], 2e-9);
        near_slice(b.velocity_mps(), &[-0.52, -0.12, 0.4], 2e-9);
        for (axis, momentum) in [-1.4, 3.1, 2.0].into_iter().enumerate() {
            near(
                2.0 * a.velocity_mps()[axis] + 3.0 * b.velocity_mps()[axis],
                momentum,
                2e-12,
            );
        }
        for (receipt, expected_work) in report.bodies.iter().zip([-2.0007, 0.2322]) {
            let SolidContactBodyImpulse::Free(impulse) = receipt else {
                panic!("free body expected")
            };
            near(impulse.energy.supplied_work_j, expected_work, 2e-9);
            near(impulse.kinetic_change_j, expected_work, 2e-9);
            assert_eq!(impulse.energy.dissipated_j, 0.0);
        }
        near(report.applied_work_j, -1.7685, 2e-9);
        near(report.contact.dissipated_energy_j, 1.7685, 2e-9);
        assert_eq!(report.support_work_j, 0.0);
        assert_eq!(report.contact.prescribed_surface_work_j, 0.0);
        assert_eq!(a.displacement_m(), &[0.0; 3]);
        assert_eq!(b.displacement_m(), &[0.0; 3]);
        assert_eq!((a.accepted_steps(), b.steps()), (0, 0));
        a.advance(&[0.0; 3], cx).unwrap();
        b.advance(&[0.0; 3], cx).unwrap();
        near_slice(a.displacement_m(), &[0.0008, 0.0173, 0.004], 2e-11);
        near_slice(b.displacement_m(), &[-0.0052, -0.0012, 0.004], 2e-11);

        // A subsequent moving external surface supplies separate source work.
        // A disconnected second body retains its exact accepted state.
        let before_b = spectral_bits(&b);
        let mut wall = pair(8, 0, 1, 0.0, 0.0);
        wall.second = SolidContactOther::PrescribedVelocity([0.5, 0.0, 0.0]);
        let source = apply_frictional_impacts(
            &mut [
                SolidContactBody::Exact(&mut a),
                SolidContactBody::Spectral(&mut b),
            ],
            &[wall],
            controls(),
            cx,
        )
        .unwrap();
        near(source.time_s, 0.01, 1e-15);
        near(source.contact.prescribed_surface_work_j, 0.42, 2e-9);
        near(source.contact.relative_work_j, -0.1764, 2e-9);
        near(source.applied_work_j, 0.2436, 2e-9);
        near(a.velocity_mps()[0], 0.5, 2e-9);
        assert_eq!(spectral_bits(&b), before_b);
    });
}

#[test]
fn g2_one_solve_couples_a_three_body_chain_and_closes_a_separating_contact() {
    with_cx(false, |cx| {
        let run = |reversed: bool| {
            let mut a = exact(1.0, [-1.0, 0.8, 0.0], cx);
            let mut b = spectral(1.0, [0.0, -0.4, 0.0], cx);
            let mut c = exact(1.0, [-0.2, 0.1, 0.0], cx);
            let ab = pair(40, 0, 1, 0.0, 2.0);
            let bc = pair(7, 1, 2, 0.0, 2.0);
            let contacts = if reversed { [bc, ab] } else { [ab, bc] };
            let report = apply_frictional_impacts(
                &mut [
                    SolidContactBody::Exact(&mut a),
                    SolidContactBody::Spectral(&mut b),
                    SolidContactBody::Exact(&mut c),
                ],
                &contacts,
                controls(),
                cx,
            )
            .unwrap();
            // Each component has W=[[2,-1],[-1,2]]. Contact BC starts with
            // +.2 separating speed, but AB drives it into simultaneous contact.
            let expected = [[0.6, -19.0 / 30.0, 0.0], [0.2, -1.0 / 15.0, 0.0]];
            for (index, impulse) in report.contact.law.impulses.iter().enumerate() {
                near_slice(
                    impulse,
                    &expected[if reversed { 1 - index } else { index }],
                    3e-9,
                );
            }
            for velocity in [a.velocity_mps(), b.velocity_mps(), c.velocity_mps()] {
                near_slice(velocity, &[-0.4, 1.0 / 6.0, 0.0], 3e-9);
            }
            near(report.applied_work_j, -193.0 / 300.0, 3e-9);
            near(report.work_balance_error_j, 0.0, 1e-12);
            (exact_bits(&a), spectral_bits(&b), exact_bits(&c))
        };
        assert_eq!(run(false), run(true));
    });
}

fn tetra_mass_and_basis(total_mass: f64) -> (Vec<f64>, Vec<f64>) {
    let mut mass = vec![0.0; 144];
    let mut basis = vec![0.0; 144];
    let modes = [
        [1.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, -1.0, -1.0],
        [1.0, -1.0, 1.0, -1.0],
        [1.0, -1.0, -1.0, 1.0],
    ];
    for i in 0..4 {
        for j in 0..4 {
            for axis in 0..3 {
                mass[(3 * i + axis) * 12 + 3 * j + axis] =
                    total_mass / 20.0 * if i == j { 2.0 } else { 1.0 };
                basis[(3 * i + axis) * 12 + 3 * j + axis] = modes[j][i]
                    * if j == 0 {
                        1.0 / total_mass.sqrt()
                    } else {
                        (5.0 / total_mass).sqrt()
                    };
            }
        }
    }
    (mass, basis)
}

fn world_momentum(
    mass: &[f64],
    positions: &[[f64; 3]; 4],
    velocity: &[f64],
    origin: [f64; 3],
) -> [f64; 6] {
    let p: Vec<f64> = (0..12)
        .map(|i| (0..12).map(|j| mass[i * 12 + j] * velocity[j]).sum())
        .collect();
    let mut result = [0.0; 6];
    for i in 0..4 {
        let r: [f64; 3] = core::array::from_fn(|a| positions[i][a] - origin[a]);
        for axis in 0..3 {
            result[axis] += p[3 * i + axis];
        }
        result[3] += r[1] * p[3 * i + 2] - r[2] * p[3 * i + 1];
        result[4] += r[2] * p[3 * i] - r[0] * p[3 * i + 2];
        result[5] += r[0] * p[3 * i + 1] - r[1] * p[3 * i];
    }
    result
}

#[test]
fn g3_consistent_mass_tetrahedra_preserve_world_linear_and_angular_momentum() {
    with_cx(false, |cx| {
        let (ma, _) = tetra_mass_and_basis(7.0);
        let (mb, phi_b) = tetra_mass_and_basis(10.0);
        let va: [f64; 12] = core::array::from_fn(|i| [0.4, 0.3, -0.2][i % 3]);
        let vb: [f64; 12] = core::array::from_fn(|i| [-0.1, 0.0, 0.0][i % 3]);
        let xa = [
            [0.0; 3],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0],
        ];
        let xb = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mut a = ExactLinearSolid::prepare(
            &assembly(&ma, 12),
            0.01,
            &[0.0; 12],
            &va,
            dense_controls(),
            cx,
        )
        .unwrap();
        let mut b = CompleteSpectralSolid::prepare_with_basis(
            &assembly(&mb, 12),
            &[0.0; 12],
            &phi_b,
            0.01,
            &[0.0; 12],
            &vb,
            spectral_controls(),
            cx,
        )
        .unwrap();
        // Both points are vertex zero at the SAME world location. Both maps
        // use n=-x, t1=y, t2=z, so their finite physical impulses are opposite.
        let rows: [[f64; 12]; 3] = core::array::from_fn(|a| {
            core::array::from_fn(|i| {
                if a == i {
                    if a == 0 { -1.0 } else { 1.0 }
                } else {
                    0.0
                }
            })
        });
        let first = SolidContactPoint {
            body: 0,
            velocity_weights: [&rows[0], &rows[1], &rows[2]],
        };
        let second = SolidContactPoint { body: 1, ..first };
        let contact = SolidContactPair {
            first,
            second: SolidContactOther::Body(second),
            ..pair(1, 0, 1, 0.2, 0.2)
        };
        let report = apply_frictional_impacts(
            &mut [
                SolidContactBody::Exact(&mut a),
                SolidContactBody::Spectral(&mut b),
            ],
            &[contact],
            controls(),
            cx,
        )
        .unwrap();
        let pn = 0.6 / (16.0 / 7.0 + 16.0 / 10.0);
        let slip = 0.13_f64.sqrt();
        near_slice(
            &report.contact.law.impulses[0],
            &[pn, -0.2 * pn * 0.3 / slip, 0.2 * pn * 0.2 / slip],
            2e-9,
        );
        for origin in [[0.0; 3], [0.7, -1.2, 0.9]] {
            let before_a = world_momentum(&ma, &xa, &va, origin);
            let before_b = world_momentum(&mb, &xb, &vb, origin);
            let after_a = world_momentum(&ma, &xa, a.velocity_mps(), origin);
            let after_b = world_momentum(&mb, &xb, b.velocity_mps(), origin);
            for axis in 0..6 {
                near(
                    after_a[axis] + after_b[axis],
                    before_a[axis] + before_b[axis],
                    3e-12,
                );
            }
        }
        // A point impulse excites the original nonuniform nodal modes; it is
        // not distributed as a rigid-body or lumped-mass velocity increment.
        assert!((a.velocity_mps()[0] - a.velocity_mps()[3]).abs() > 0.1);
        assert!((b.velocity_mps()[0] - b.velocity_mps()[3]).abs() > 0.1);
        near(report.applied_work_j, report.contact.relative_work_j, 2e-12);
        near(report.kinetic_change_j, report.applied_work_j, 2e-12);
    });
}

fn supported(velocity: [f64; 4], cx: &Cx<'_>) -> SupportedSpectralSolid {
    let mass = [
        2.0, 0.0, 0.0, 0.35, 0.0, 3.0, 0.0, -0.2, 0.0, 0.0, 4.0, 0.1, 0.35, -0.2, 0.1, 1.5,
    ];
    SupportedSpectralSolid::prepare_with_basis(
        &assembly(&mass, 4),
        &[3],
        &[0.0; 3],
        &[
            1.0 / 2.0_f64.sqrt(),
            0.0,
            0.0,
            0.0,
            1.0 / 3.0_f64.sqrt(),
            0.0,
            0.0,
            0.0,
            0.5,
        ],
        0.01,
        &[0.0; 4],
        &velocity,
        &[0.0],
        SupportedSolidControls {
            spectral: spectral_controls(),
            maximum_total_dofs: 8,
            maximum_preparation_bytes: 4 * 1024 * 1024,
            maximum_preparation_multiply_terms: 10_000_000,
            maximum_step_multiply_terms: 2_000_000,
            maximum_free_force_residual_n: 1e-10,
            maximum_work_error_j: 1e-8,
        },
        cx,
    )
    .unwrap()
}

#[test]
fn g4_finite_pair_keeps_support_cross_mass_reaction_work_and_boost_dependence() {
    with_cx(false, |cx| {
        let rows = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ];
        let run = |boost: f64| {
            let av = [-1.0 + boost, 2.0 + boost, boost, 0.4 + boost];
            let bv = [0.2 + boost, -0.3 + boost, boost];
            let mut a = supported(av, cx);
            let mut b = spectral(3.0, bv, cx);
            let contact = SolidContactPair {
                first: SolidContactPoint {
                    body: 0,
                    velocity_weights: [&rows[0], &rows[1], &rows[2]],
                },
                ..pair(1, 0, 1, 0.5, 0.25)
            };
            let report = apply_frictional_impacts(
                &mut [
                    SolidContactBody::Supported(&mut a),
                    SolidContactBody::Spectral(&mut b),
                ],
                &[contact],
                controls(),
                cx,
            )
            .unwrap();
            near_slice(&report.contact.law.impulses[0], &[2.16, -0.54, 0.0], 2e-9);
            near_slice(
                a.velocity_mps(),
                &[0.08 + boost, 1.82 + boost, boost, 0.4 + boost],
                2e-9,
            );
            near_slice(
                b.velocity_mps(),
                &[-0.52 + boost, -0.12 + boost, boost],
                2e-9,
            );
            let SolidContactBodyImpulse::Supported(receipt) = &report.bodies[0] else {
                panic!("supported body expected")
            };
            near(receipt.support_reaction_kg_m_s[0], 0.414, 2e-10);
            near(receipt.momentum_change_kg_m_s[3], 0.414, 2e-10);
            near(receipt.energy.support_work_j, 0.414 * (0.4 + boost), 2e-10);
            near(report.applied_work_j, -1.7928, 3e-9);
            near(report.support_work_j, 0.414 * (0.4 + boost), 3e-9);
            near(
                report.kinetic_change_j,
                -1.7928 + 0.414 * (0.4 + boost),
                3e-9,
            );
            assert_eq!(report.contact.prescribed_surface_work_j, 0.0);
            near(a.cumulative_support_work_j(), report.support_work_j, 1e-12);
            assert_eq!(a.support_acceleration_mps2(), &[0.0]);
            a.advance_affine_supports(&[0.0], &[0.0; 4], &[0.0; 4], cx)
                .unwrap();
            b.advance(&[0.0; 3], cx).unwrap();
            for (position, velocity) in a.displacement_m().iter().zip(a.velocity_mps()) {
                near(*position, 0.01 * velocity, 2e-10);
            }
            report
        };
        let rest = run(0.0);
        let boosted = run(1.3);
        near(
            rest.contact.relative_work_j,
            boosted.contact.relative_work_j,
            2e-12,
        );
        near(
            boosted.support_work_j - rest.support_work_j,
            1.3 * 0.414,
            2e-12,
        );
    });
}

#[test]
fn g5_late_native_readback_failure_rolls_back_all_owners_and_preserves_retry() {
    with_cx(false, |cx| {
        let make_spectral = || {
            let mut limits = spectral_controls();
            limits.maximum_relative_operator_error = 1e-3;
            limits.maximum_relative_state_error = 1e-3;
            limits.maximum_work_error_j = 1e-3;
            CompleteSpectralSolid::prepare_with_basis(
                &assembly(&diagonal(1.0), 3),
                &[0.0; 3],
                &[1.0001, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
                0.01,
                &[0.0; 3],
                &[-0.2, 0.0, 0.0],
                limits,
                cx,
            )
            .unwrap()
        };
        let mut a = exact(1.0, [-0.5, 0.0, 0.0], cx);
        let mut b = make_spectral();
        let (mut reference_a, mut reference_b) =
            (exact(1.0, [-0.5, 0.0, 0.0], cx), make_spectral());
        let (before_a, before_b) = (exact_bits(&a), spectral_bits(&b));
        let contact = pair(1, 0, 1, 0.0, 0.0);
        let mut limits = controls();
        limits.impact.maximum_impulse_error_kg_m_s = 1e-5;
        limits.impact.maximum_work_error_j = 1e-3;
        assert!(matches!(
            apply_frictional_impacts(
                &mut [
                    SolidContactBody::Exact(&mut a),
                    SolidContactBody::Spectral(&mut b)
                ],
                &[contact],
                limits,
                cx,
            ),
            Err(SolidContactNetworkError::Native {
                body: 1,
                error: LinearSolidImpulseError::Momentum { .. },
            })
        ));
        assert_eq!(exact_bits(&a), before_a);
        assert_eq!(spectral_bits(&b), before_b);
        limits.impact.maximum_impulse_error_kg_m_s = 1e-3;
        // Native trials can admit this rounded basis, while the actual pair
        // speed differs from the predictor after the FIRST body's valid trial.
        assert!(matches!(
            apply_frictional_impacts(
                &mut [
                    SolidContactBody::Exact(&mut a),
                    SolidContactBody::Spectral(&mut b)
                ],
                &[contact],
                limits,
                cx
            ),
            Err(SolidContactNetworkError::Contact(
                SolidFrictionalImpactError::Speed { .. }
            ))
        ));
        assert_eq!(exact_bits(&a), before_a);
        assert_eq!(spectral_bits(&b), before_b);
        limits.impact.maximum_contact_speed_error_mps = 1e-3;
        assert!(matches!(
            apply_frictional_impacts(
                &mut [
                    SolidContactBody::Exact(&mut a),
                    SolidContactBody::Spectral(&mut b)
                ],
                &[contact],
                limits,
                cx
            ),
            Err(SolidContactNetworkError::Contact(
                SolidFrictionalImpactError::Contact(FrictionalContactError::ResidualNotSatisfied(
                    _
                ))
            ))
        ));
        for kind in 0..5 {
            let mut bounded = limits;
            match kind {
                0 => bounded.maximum_bodies = 1,
                1 => bounded.maximum_total_dofs = 5,
                2 => bounded.maximum_workspace_bytes = 0,
                3 => bounded.impact.maximum_workspace_bytes = 0,
                _ => bounded.impact.maximum_response_multiply_terms = 0,
            }
            assert!(
                apply_frictional_impacts(
                    &mut [
                        SolidContactBody::Exact(&mut a),
                        SolidContactBody::Spectral(&mut b)
                    ],
                    &[contact],
                    bounded,
                    cx
                )
                .is_err()
            );
        }
        for first in [
            SolidContactPoint {
                body: 2,
                ..point(0)
            },
            SolidContactPoint {
                body: 0,
                velocity_weights: [&[], &IDENTITY[1], &IDENTITY[2]],
            },
        ] {
            assert!(matches!(
                apply_frictional_impacts(
                    &mut [
                        SolidContactBody::Exact(&mut a),
                        SolidContactBody::Spectral(&mut b)
                    ],
                    &[SolidContactPair { first, ..contact }],
                    limits,
                    cx,
                ),
                Err(SolidContactNetworkError::Invalid(_))
            ));
        }
        with_cx(true, |cancelled| {
            assert!(matches!(
                apply_frictional_impacts(
                    &mut [
                        SolidContactBody::Exact(&mut a),
                        SolidContactBody::Spectral(&mut b)
                    ],
                    &[contact],
                    limits,
                    cancelled
                ),
                Err(SolidContactNetworkError::Cancelled)
            ));
        });
        assert_eq!(exact_bits(&a), before_a);
        assert_eq!(spectral_bits(&b), before_b);
        limits.impact.contact.velocity_tolerance = 1e-3;
        limits.impact.contact.complementarity_tolerance = 1e-4;
        for (first, second) in [(&mut a, &mut b), (&mut reference_a, &mut reference_b)] {
            apply_frictional_impacts(
                &mut [
                    SolidContactBody::Exact(first),
                    SolidContactBody::Spectral(second),
                ],
                &[contact],
                limits,
                cx,
            )
            .unwrap();
        }
        assert_eq!(exact_bits(&a), exact_bits(&reference_a));
        assert_eq!(spectral_bits(&b), spectral_bits(&reference_b));
        a.advance(&[0.0; 3], cx).unwrap();
        reference_a.advance(&[0.0; 3], cx).unwrap();
        let (a_before_clock_error, b_before_clock_error) = (exact_bits(&a), spectral_bits(&b));
        assert!(matches!(
            apply_frictional_impacts(
                &mut [
                    SolidContactBody::Exact(&mut a),
                    SolidContactBody::Spectral(&mut b)
                ],
                &[contact],
                limits,
                cx
            ),
            Err(SolidContactNetworkError::UnsynchronizedTime { .. })
        ));
        assert_eq!(exact_bits(&a), a_before_clock_error);
        assert_eq!(spectral_bits(&b), b_before_clock_error);
        b.advance(&[0.0; 3], cx).unwrap();
        reference_b.advance(&[0.0; 3], cx).unwrap();
        assert_eq!(exact_bits(&a), exact_bits(&reference_a));
        assert_eq!(spectral_bits(&b), spectral_bits(&reference_b));
    });
}
