//! Original-mass, measured-velocity and work checks for native solid impacts.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::impact::{
    SolidNormalContact, SolidNormalImpactControls, SolidNormalImpactError,
};
use fs_solid::dynamics::spectral::{CompleteSpectralSolid, CompleteSpectralSolidControls};
use fs_solid::dynamics::{ExactLinearSolid, ExactLinearSolidControls};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::Coo;

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
                kernel_id: 39,
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
        maximum_work_error_j: 1e-9,
        maximum_steps: 100,
    }
}
fn spectral_controls() -> CompleteSpectralSolidControls {
    CompleteSpectralSolidControls {
        maximum_free_dofs: 24,
        maximum_preparation_bytes: 1024 * 1024,
        maximum_preparation_multiply_terms: 1_000_000,
        maximum_step_multiply_terms: 1_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_absolute_eigen_residual: 0.0,
        maximum_relative_state_error: 1e-12,
        maximum_work_error_j: 1e-9,
        maximum_steps: 100,
    }
}
fn impact_controls() -> SolidNormalImpactControls {
    SolidNormalImpactControls {
        contact: MultiImpactConfig::default(),
        maximum_impulse_error_kg_m_s: 1e-11,
        maximum_contact_speed_error_mps: 1e-11,
        maximum_work_error_j: 1e-10,
        maximum_workspace_bytes: 1024 * 1024,
        maximum_response_multiply_terms: 1_000_000,
    }
}
fn contact(row: &[f64], speed: f64, restitution: f64) -> SolidNormalContact<'_> {
    SolidNormalContact {
        normal_velocity_weights: row,
        gap_m: 0.0,
        activation_tolerance_m: 0.0,
        prescribed_surface_speed_mps: speed,
        restitution,
    }
}
fn assembly(k: &[f64], m: &[f64], n: usize) -> TetElasticAssembly {
    let (mut stiffness, mut mass) = (Coo::new(n, n), Coo::new(n, n));
    for i in 0..n {
        for j in 0..n {
            stiffness.push(i, j, k[i * n + j]);
            mass.push(i, j, m[i * n + j]);
        }
    }
    TetElasticAssembly {
        stiffness: stiffness.assemble(),
        mass: mass.assemble(),
        free_dofs: (0..n).collect(),
        total_mass_kg: if n == 12 { 20.0 } else { 1.0 },
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0 / 6.0],
    }
}
fn pencil() -> TetElasticAssembly {
    // Phi=[[1,1/2],[0,1]], M=Phi^-T Phi^-1, lambda=[0,9].
    assembly(&[0.0, 0.0, 0.0, 9.0], &[1.0, -0.5, -0.5, 1.25], 2)
}
fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "actual={actual:e}, expected={expected:e}, allowance={allowance:e}"
    );
}
fn bits(solid: &CompleteSpectralSolid) -> Vec<u64> {
    solid
        .displacement_m()
        .iter()
        .chain(solid.velocity_mps())
        .copied()
        .chain([
            solid.time_s(),
            solid.stored_energy_j(),
            solid.kinetic_energy_j(),
            solid.cumulative_work_j(),
            solid.steps() as f64,
        ])
        .map(f64::to_bits)
        .collect()
}

#[test]
fn g1_simultaneous_cross_mass_impacts_retain_moving_surface_work_and_elastic_continuation() {
    with_cx(false, |cx| {
        let pencil = pencil();
        let mut dense = ExactLinearSolid::prepare(
            &pencil,
            0.03,
            &[0.1, -0.2],
            &[-1.0, -0.1],
            dense_controls(),
            cx,
        )
        .unwrap();
        let mut spectral = CompleteSpectralSolid::prepare_with_basis(
            &pencil,
            &[0.0, 9.0],
            &[1.0, 0.5, 0.0, 1.0],
            0.03,
            &[0.1, -0.2],
            &[-1.0, -0.1],
            spectral_controls(),
            cx,
        )
        .unwrap();
        // W=[[1.25,-.5],[-.5,1]]. Row 2 initially separates at .3 m/s,
        // but row 1's impulse would make it close unless both are solved together.
        let contacts = [
            contact(&[1.0, 0.0], 0.25, 0.0),
            contact(&[0.0, -1.0], -0.2, 0.0),
        ];
        with_cx(true, |cancelled| {
            assert!(matches!(
                dense.apply_normal_impacts(&contacts, impact_controls(), cancelled),
                Err(SolidNormalImpactError::Cancelled)
            ));
        });
        assert_eq!(dense.velocity_mps(), &[-1.0, -0.1]);
        let a = dense
            .apply_normal_impacts(&contacts, impact_controls(), cx)
            .unwrap();
        let b = spectral
            .apply_normal_impacts(&contacts, impact_controls(), cx)
            .unwrap();
        for (velocity, displacement, report) in [
            (dense.velocity_mps(), dense.displacement_m(), &a),
            (spectral.velocity_mps(), spectral.displacement_m(), &b),
        ] {
            near(velocity[0], 0.25, 2e-12);
            near(velocity[1], 0.2, 2e-12);
            assert_eq!(displacement, &[0.1, -0.2]);
            near(report.contact.law.normal_impulses[0], 1.1, 2e-12);
            near(report.contact.law.normal_impulses[1], 0.25, 2e-12);
            near(report.contact.prescribed_surface_work_j, 0.225, 2e-12);
            near(report.contact.relative_work_j, -0.65, 2e-12);
            near(report.contact.dissipated_energy_j, 0.65, 2e-12);
            near(report.impulse.energy.supplied_work_j, -0.425, 2e-12);
            near(report.impulse.kinetic_change_j, -0.425, 2e-12);
            near(report.contact.work_balance_error_j, 0.0, 2e-12);
            assert_eq!(report.impulse.energy.dissipated_j, 0.0);
            assert_eq!(report.impulse.time_s, 0.0);
        }
        assert_eq!(dense.accepted_steps(), 0);
        assert_eq!(spectral.steps(), 0);
        dense.advance(&[0.0; 2], cx).unwrap();
        spectral.advance(&[0.0; 2], cx).unwrap();
        // Independent rigid+elastic analytic motion after the impulse.
        let (sin, cos) = 0.09_f64.sin_cos();
        let q0 = 0.2 + 0.15 * 0.03;
        let q1 = -0.2 * cos + (0.2 / 3.0) * sin;
        let w1 = 0.6 * sin + 0.2 * cos;
        for (u, v) in [
            (dense.displacement_m(), dense.velocity_mps()),
            (spectral.displacement_m(), spectral.velocity_mps()),
        ] {
            near(u[0], q0 + 0.5 * q1, 2e-12);
            near(u[1], q1, 2e-12);
            near(v[0], 0.15 + 0.5 * w1, 2e-12);
            near(v[1], w1, 2e-12);
        }
    });
}

#[test]
fn g3_redundant_normal_rows_preserve_unique_motion_and_relative_work_under_a_boost() {
    with_cx(false, |cx| {
        let pencil = pencil();
        for boost in [0.0, 0.7] {
            let make = || {
                CompleteSpectralSolid::prepare_with_basis(
                    &pencil,
                    &[0.0, 9.0],
                    &[1.0, 0.5, 0.0, 1.0],
                    0.03,
                    &[0.0; 2],
                    &[-0.2 + boost, 0.1],
                    spectral_controls(),
                    cx,
                )
                .unwrap()
            };
            let (mut single, mut duplicate) = (make(), make());
            let row = contact(&[1.0, 0.0], 0.15 + boost, 0.4);
            let a = single
                .apply_normal_impacts(&[row], impact_controls(), cx)
                .unwrap();
            let b = duplicate
                .apply_normal_impacts(&[row, row], impact_controls(), cx)
                .unwrap();
            assert!(matches!(
                b.contact.law.determinacy,
                MultiImpactDeterminacy::UniqueVelocity(_)
            ));
            near(b.contact.law.normal_impulses.iter().sum(), 0.392, 2e-11);
            for axis in 0..2 {
                near(
                    single.velocity_mps()[axis],
                    duplicate.velocity_mps()[axis],
                    2e-11,
                );
            }
            near(duplicate.velocity_mps()[0] - boost, 0.29, 2e-11);
            near(duplicate.velocity_mps()[1], 0.296, 2e-11);
            near(b.contact.relative_work_j, a.contact.relative_work_j, 2e-11);
            near(b.contact.relative_work_j, -0.04116, 2e-11);
            near(
                b.contact.prescribed_surface_work_j,
                0.392 * (0.15 + boost),
                2e-11,
            );
        }
    });
}

#[test]
fn g1_p1_point_contact_transfers_world_linear_and_angular_impulse_through_consistent_mass() {
    with_cx(false, |cx| {
        // The exact reference P1 mass of a 20 kg tetrahedron is
        // M_ab=(1+delta_ab) I_3. Its inverse is (delta_ab-1/5) I_3.
        let mut mass = [0.0; 144];
        for a in 0..4 {
            for b in 0..4 {
                for d in 0..3 {
                    mass[(3 * a + d) * 12 + 3 * b + d] = if a == b { 2.0 } else { 1.0 };
                }
            }
        }
        let pencil = assembly(&[0.0; 144], &mass, 12);
        let mut initial = [0.0; 12];
        for node in 0..4 {
            initial[3 * node + 2] = -0.1;
        }
        let mut solid =
            ExactLinearSolid::prepare(&pencil, 0.02, &[0.0; 12], &initial, dense_controls(), cx)
                .unwrap();
        let positions = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let inventory = |solid: &ExactLinearSolid, origin: [f64; 3]| {
            let (mut p, mut l) = ([0.0; 3], [0.0; 3]);
            for node in 0..4 {
                let nodal_p: [f64; 3] = core::array::from_fn(|d| {
                    (0..12)
                        .map(|j| mass[(3 * node + d) * 12 + j] * solid.velocity_mps()[j])
                        .sum()
                });
                let x: [f64; 3] = core::array::from_fn(|d| {
                    positions[node][d] + solid.displacement_m()[3 * node + d] - origin[d]
                });
                for d in 0..3 {
                    p[d] += nodal_p[d];
                    l[d] += x[(d + 1) % 3] * nodal_p[(d + 2) % 3]
                        - x[(d + 2) % 3] * nodal_p[(d + 1) % 3];
                }
            }
            (p, l)
        };
        let origins = [[0.0; 3], [-3.0, 2.0, 1.0]];
        let before = origins.map(|origin| inventory(&solid, origin));
        let weights = [0.1, 0.2, 0.3, 0.4];
        let mut row = [0.0; 12];
        for node in 0..4 {
            row[3 * node + 2] = weights[node];
        }
        // W=sum(N_a^2)-1/5=.1/kg, so an elastic .1 m/s approach needs p=2.
        let report = solid
            .apply_normal_impacts(&[contact(&row, 0.0, 1.0)], impact_controls(), cx)
            .unwrap();
        near(report.contact.law.normal_impulses[0], 2.0, 2e-12);
        near(report.contact.dissipated_energy_j, 0.0, 2e-12);
        for advance in [false, true] {
            if advance {
                solid.advance(&[0.0; 12], cx).unwrap();
            }
            for (origin, old) in origins.into_iter().zip(before) {
                let actual = inventory(&solid, origin);
                for d in 0..3 {
                    near(actual.0[d] - old.0[d], [0.0, 0.0, 2.0][d], 2e-11);
                    near(
                        actual.1[d] - old.1[d],
                        [2.0 * (0.3 - origin[1]), -2.0 * (0.2 - origin[0]), 0.0][d],
                        2e-11,
                    );
                }
            }
        }
    });
}

#[test]
fn g4_late_contact_readback_refusal_preserves_physical_modal_and_cumulative_state() {
    with_cx(false, |cx| {
        let pencil = assembly(&[0.0], &[1.0], 1);
        let phi = 1.0001;
        let make = || {
            let mut controls = spectral_controls();
            controls.maximum_relative_operator_error = 1e-3;
            controls.maximum_relative_state_error = 1e-3;
            controls.maximum_work_error_j = 1e-3;
            CompleteSpectralSolid::prepare_with_basis(
                &pencil,
                &[0.0],
                &[phi],
                0.05,
                &[0.0],
                &[-0.2],
                controls,
                cx,
            )
            .unwrap()
        };
        let (mut refused, mut reference, mut native) = (make(), make(), make());
        let row = contact(&[1.0], 0.0, 0.0);
        let mut controls = impact_controls();
        controls.maximum_impulse_error_kg_m_s = 1e-5;
        controls.maximum_work_error_j = 1e-3;
        // Native impulse admits this rounded basis at the looser momentum/work
        // budgets, but its ACTUAL contact velocity misses restitution by ~4e-5.
        native
            .apply_momentum_impulse(&[0.2 / (phi * phi)], 1e-5, cx)
            .unwrap();
        assert!(native.velocity_mps()[0] < -1e-5);
        let before = bits(&refused);
        assert!(matches!(
            refused.apply_normal_impacts(&[row], controls, cx),
            Err(SolidNormalImpactError::Speed { .. })
        ));
        assert_eq!(bits(&refused), before);
        let mut budget = controls;
        budget.maximum_workspace_bytes = 0;
        assert!(matches!(
            refused.apply_normal_impacts(&[row], budget, cx),
            Err(SolidNormalImpactError::Budget { .. })
        ));
        let mut invalid = row;
        invalid.gap_m = -1e-20;
        assert!(
            refused
                .apply_normal_impacts(&[invalid], controls, cx)
                .is_err()
        );
        with_cx(true, |cancelled| {
            assert!(matches!(
                refused.apply_normal_impacts(&[row], controls, cancelled),
                Err(SolidNormalImpactError::Cancelled)
            ));
        });
        assert_eq!(bits(&refused), before);
        controls.maximum_contact_speed_error_mps = 1e-3;
        refused.apply_normal_impacts(&[row], controls, cx).unwrap();
        reference
            .apply_normal_impacts(&[row], controls, cx)
            .unwrap();
        assert_eq!(bits(&refused), bits(&reference));
        refused.advance(&[0.0], cx).unwrap();
        reference.advance(&[0.0], cx).unwrap();
        assert_eq!(bits(&refused), bits(&reference));
    });
}
