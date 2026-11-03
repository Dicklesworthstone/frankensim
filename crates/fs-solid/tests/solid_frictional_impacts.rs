//! Physical complete-coordinate Coulomb impulses, continued original dynamics,
//! independent world momentum, and transactional actual-endpoint refusals.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactDeterminacy, FrictionalContactError,
    NormalVelocityLaw,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::frictional_impact::{
    SolidFrictionalContact, SolidFrictionalImpactControls, SolidFrictionalImpactError,
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
                kernel_id: 44,
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
        maximum_preparation_bytes: 1024 * 1024,
        maximum_preparation_multiply_terms: 1_000_000,
        maximum_step_multiply_terms: 1_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_absolute_eigen_residual: 0.0,
        maximum_relative_state_error: 1e-12,
        maximum_work_error_j: 1e-8,
        maximum_steps: 100,
    }
}
fn impact_controls() -> SolidFrictionalImpactControls {
    SolidFrictionalImpactControls {
        contact: FrictionalContactConfig::default(),
        maximum_impulse_error_kg_m_s: 1e-10,
        maximum_contact_speed_error_mps: 1e-10,
        maximum_work_error_j: 1e-9,
        maximum_workspace_bytes: 1024 * 1024,
        maximum_response_multiply_terms: 1_000_000,
    }
}
fn contact<'a>(
    id: u64,
    rows: [&'a [f64]; 3],
    surface: [f64; 3],
    mu: f64,
) -> SolidFrictionalContact<'a> {
    SolidFrictionalContact {
        id,
        velocity_weights: rows,
        gap_m: 0.0,
        activation_tolerance_m: 0.0,
        prescribed_surface_velocity_mps: surface,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
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
fn action(matrix: &[f64], vector: &[f64]) -> Vec<f64> {
    let n = vector.len();
    (0..n)
        .map(|i| (0..n).map(|j| matrix[i * n + j] * vector[j]).sum())
        .collect()
}

#[test]
fn g1_simultaneous_stick_and_slip_keep_cross_mass_source_work_and_all_modes() {
    with_cx(false, |cx| {
        // A complete analytically invertible basis. It produces W_01=.2,
        // W_03=.1 and W_34=-.15: both normal/shear and intercontact coupling.
        let phi = [
            1.0, 0.2, 0.0, 0.1, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 1.0, -0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0,
        ];
        let inverse = [
            1.0, -0.2, 0.0, -0.1, -0.015, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0,
        ];
        let lambda = [0.0, 4.0, 9.0, 16.0, 25.0, 0.0];
        let (mut m, mut k) = ([0.0; 36], [0.0; 36]);
        for i in 0..6 {
            for j in 0..6 {
                for a in 0..6 {
                    m[6 * i + j] += inverse[6 * a + i] * inverse[6 * a + j];
                    k[6 * i + j] += inverse[6 * a + i] * lambda[a] * inverse[6 * a + j];
                }
            }
        }
        let pencil = assembly(&k, &m, 6);
        let u = [0.03, -0.04, 0.05, 0.06, -0.02, 0.01];
        let before = [-1.02, -0.35, 0.2, -1.0824, 0.566, 0.428];
        let expected = [0.1, -0.2, 0.3, -0.15, 0.35, 0.3];
        let mut dense =
            ExactLinearSolid::prepare(&pencil, 0.03, &u, &before, dense_controls(), cx).unwrap();
        let make_spectral = || {
            CompleteSpectralSolid::prepare_with_basis(
                &pencil,
                &lambda,
                &phi,
                0.03,
                &u,
                &before,
                spectral_controls(),
                cx,
            )
            .unwrap()
        };
        let (mut spectral, mut permuted) = (make_spectral(), make_spectral());
        let rows: [[f64; 6]; 6] =
            core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
        let contacts = [
            contact(40, [&rows[0], &rows[1], &rows[2]], [0.1, -0.2, 0.3], 0.4),
            contact(7, [&rows[3], &rows[4], &rows[5]], [-0.15, 0.05, -0.1], 0.2),
        ];
        let a = dense
            .apply_frictional_impacts(&contacts, impact_controls(), cx)
            .unwrap();
        let b = spectral
            .apply_frictional_impacts(&contacts, impact_controls(), cx)
            .unwrap();
        permuted
            .apply_frictional_impacts(&[contacts[1], contacts[0]], impact_controls(), cx)
            .unwrap();
        assert_eq!(bits(&spectral), bits(&permuted));
        for (velocity, displacement, report) in [
            (dense.velocity_mps(), dense.displacement_m(), &a),
            (spectral.velocity_mps(), spectral.displacement_m(), &b),
        ] {
            for (actual, expected) in velocity.iter().zip(expected) {
                near(*actual, expected, 2e-9);
            }
            assert_eq!(displacement, &u);
            for (actual, expected) in report
                .contact
                .law
                .impulses
                .iter()
                .zip([[1.0, -0.05, 0.1], [0.8, -0.096, -0.128]])
            {
                for axis in 0..3 {
                    near(actual[axis], expected[axis], 2e-9);
                }
            }
            near(report.contact.relative_work_j, -1.03277, 2e-9);
            near(report.contact.prescribed_surface_work_j, 0.028, 2e-9);
            near(report.impulse.energy.supplied_work_j, -1.00477, 2e-9);
            near(report.impulse.kinetic_change_j, -1.00477, 2e-9);
            near(report.contact.work_balance_error_j, 0.0, 1e-12);
            assert!(report.contact.measured_residuals.tangential_duality_gap < 1e-9);
            assert_eq!(report.impulse.energy.dissipated_j, 0.0);
            assert_eq!(report.impulse.time_s, 0.0);
        }
        assert_eq!(dense.accepted_steps(), 0);
        assert_eq!(spectral.steps(), 0);
        dense.advance(&[0.0; 6], cx).unwrap();
        spectral.advance(&[0.0; 6], cx).unwrap();
        let q0 = action(&inverse, &u);
        let dq0 = action(&inverse, &expected);
        let (mut q, mut dq) = ([0.0; 6], [0.0; 6]);
        for i in 0..6 {
            let omega: f64 = lambda[i].sqrt();
            if omega == 0.0 {
                q[i] = q0[i] + 0.03 * dq0[i];
                dq[i] = dq0[i];
            } else {
                let (s, c) = (0.03 * omega).sin_cos();
                q[i] = q0[i] * c + dq0[i] * s / omega;
                dq[i] = -omega * q0[i] * s + dq0[i] * c;
            }
        }
        let (expected_u, expected_v) = (action(&phi, &q), action(&phi, &dq));
        for (u, v) in [
            (dense.displacement_m(), dense.velocity_mps()),
            (spectral.displacement_m(), spectral.velocity_mps()),
        ] {
            for i in 0..6 {
                near(u[i], expected_u[i], 2e-9);
                near(v[i], expected_v[i], 2e-9);
            }
        }
    });
}

#[test]
fn g1_p1_barycentric_friction_impulse_closes_world_linear_and_angular_momentum() {
    with_cx(false, |cx| {
        // Exact consistent mass of a 20 kg P1 tetrahedron. M_ab=(1+delta_ab)I;
        // M^-1_ab=(delta_ab-1/5)I. Barycentric point mobility is .1 I/kg.
        let mut mass = [0.0; 144];
        for a in 0..4 {
            for b in 0..4 {
                for axis in 0..3 {
                    mass[(3 * a + axis) * 12 + 3 * b + axis] = if a == b { 2.0 } else { 1.0 };
                }
            }
        }
        let pencil = assembly(&[0.0; 144], &mass, 12);
        let initial: [f64; 12] = core::array::from_fn(|i| [0.3, 0.4, -0.1][i % 3]);
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
        let barycentric = [0.1, 0.2, 0.3, 0.4];
        let mut rows = [[0.0; 12]; 3];
        for node in 0..4 {
            for a in 0..3 {
                rows[a][3 * node + [2, 0, 1][a]] = barycentric[node];
            }
        }
        let mut c = contact(7, [&rows[0], &rows[1], &rows[2]], [0.0, 0.05, -0.1], 0.2);
        c.normal_law = NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 };
        let report = solid
            .apply_frictional_impacts(&[c], impact_controls(), cx)
            .unwrap();
        let slip = 0.25_f64.hypot(0.5);
        let world = [-0.3 * 0.25 / slip, -0.3 * 0.5 / slip, 1.5];
        for (a, d) in [2, 0, 1].into_iter().enumerate() {
            near(report.contact.law.impulses[0][a], world[d], 2e-9);
        }
        near(
            report.contact.prescribed_surface_work_j,
            0.05 * world[0] - 0.1 * world[1],
            2e-9,
        );
        near(
            report.contact.relative_work_j,
            -0.0375 - 0.3 * slip + 0.0045,
            2e-9,
        );
        for advance in [false, true] {
            if advance {
                solid.advance(&[0.0; 12], cx).unwrap();
            }
            for (origin, old) in origins.into_iter().zip(before) {
                let actual = inventory(&solid, origin);
                let arm = [0.2 - origin[0], 0.3 - origin[1], 0.4 - origin[2]];
                for d in 0..3 {
                    near(actual.0[d] - old.0[d], world[d], 5e-9);
                    near(
                        actual.1[d] - old.1[d],
                        arm[(d + 1) % 3] * world[(d + 2) % 3]
                            - arm[(d + 2) % 3] * world[(d + 1) % 3],
                        2e-8,
                    );
                }
            }
        }
    });
}

#[test]
fn g3_redundant_contact_frames_keep_unique_motion_work_and_tangent_gauge_covariance() {
    with_cx(false, |cx| {
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let pencil = assembly(&[0.0; 9], &identity, 3);
        for boost in [[0.0; 3], [1.7, -0.8, 0.6]] {
            for angle in [0.0_f64, 0.71] {
                let (s, c) = angle.sin_cos();
                let rows = [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]];
                let surface_world = [0.3 + boost[0], -0.1 + boost[1], 0.2 + boost[2]];
                let surface =
                    core::array::from_fn(|i| (0..3).map(|a| rows[i][a] * surface_world[a]).sum());
                for mu in [0.0, 0.25] {
                    let initial = [-0.7 + boost[0], 0.6 + boost[1], 0.8 + boost[2]];
                    let make = || {
                        CompleteSpectralSolid::prepare_with_basis(
                            &pencil,
                            &[0.0; 3],
                            &identity,
                            0.02,
                            &[0.0; 3],
                            &initial,
                            spectral_controls(),
                            cx,
                        )
                        .unwrap()
                    };
                    let (mut single, mut redundant) = (make(), make());
                    let c = contact(17, [&rows[0], &rows[1], &rows[2]], surface, mu);
                    let a = single
                        .apply_frictional_impacts(&[c], impact_controls(), cx)
                        .unwrap();
                    let mut other = c;
                    other.id = 5;
                    let b = redundant
                        .apply_frictional_impacts(&[c, other], impact_controls(), cx)
                        .unwrap();
                    assert!(matches!(
                        b.contact.law.determinacy,
                        FrictionalContactDeterminacy::UniqueVelocity(_)
                    ));
                    assert!(b.contact.law.impulse_work_witness.is_some());
                    let speed = 0.7_f64.hypot(0.6);
                    let impulse = [1.0, -mu * 0.7 / speed, -mu * 0.6 / speed];
                    for axis in 0..3 {
                        near(
                            single.velocity_mps()[axis],
                            initial[axis] + impulse[axis],
                            3e-9,
                        );
                        near(
                            redundant.velocity_mps()[axis],
                            single.velocity_mps()[axis],
                            3e-9,
                        );
                    }
                    near(
                        b.contact.relative_work_j,
                        -0.5 - mu * speed + 0.5 * mu * mu,
                        3e-9,
                    );
                    near(b.contact.relative_work_j, a.contact.relative_work_j, 3e-9);
                    near(
                        b.contact.prescribed_surface_work_j,
                        (0..3).map(|i| impulse[i] * surface_world[i]).sum(),
                        3e-9,
                    );
                    near(b.contact.work_balance_error_j, 0.0, 1e-12);
                }
            }
        }
    });
}

#[test]
fn g4_native_readback_refusal_keeps_all_state_and_cannot_relax_the_coulomb_law() {
    with_cx(false, |cx| {
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let pencil = assembly(&[0.0; 9], &identity, 3);
        let phi = 1.0001;
        let make = || {
            let mut controls = spectral_controls();
            controls.maximum_relative_operator_error = 1e-3;
            controls.maximum_relative_state_error = 1e-3;
            controls.maximum_work_error_j = 1e-3;
            CompleteSpectralSolid::prepare_with_basis(
                &pencil,
                &[0.0; 3],
                &[phi, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
                0.05,
                &[0.0; 3],
                &[-0.2, 0.3, 0.4],
                controls,
                cx,
            )
            .unwrap()
        };
        let (mut refused, mut reference, mut native) = (make(), make(), make());
        let c = contact(
            1,
            [&identity[0..3], &identity[3..6], &identity[6..9]],
            [0.0; 3],
            0.2,
        );
        let mut controls = impact_controls();
        controls.maximum_impulse_error_kg_m_s = 1e-5;
        controls.maximum_work_error_j = 1e-3;
        let pn = 0.2 / (phi * phi);
        native
            .apply_momentum_impulse(&[pn, -0.2 * pn * 0.6, -0.2 * pn * 0.8], 1e-5, cx)
            .unwrap();
        assert!(native.velocity_mps()[0] < -1e-5);
        let before = bits(&refused);
        assert!(matches!(
            refused.apply_frictional_impacts(&[c], controls, cx),
            Err(SolidFrictionalImpactError::Speed { .. })
        ));
        assert_eq!(bits(&refused), before);
        let mut loose_readback = controls;
        loose_readback.maximum_contact_speed_error_mps = 1e-3;
        assert!(matches!(
            refused.apply_frictional_impacts(&[c], loose_readback, cx),
            Err(SolidFrictionalImpactError::Contact(
                FrictionalContactError::ResidualNotSatisfied(_)
            ))
        ));
        let mut budget = controls;
        budget.maximum_workspace_bytes = 0;
        assert!(matches!(
            refused.apply_frictional_impacts(&[c], budget, cx),
            Err(SolidFrictionalImpactError::Budget { .. })
        ));
        let mut invalid = c;
        invalid.gap_m = -1e-20;
        assert!(
            refused
                .apply_frictional_impacts(&[invalid], controls, cx)
                .is_err()
        );
        assert!(matches!(
            refused.apply_frictional_impacts(&[c, c], controls, cx),
            Err(SolidFrictionalImpactError::Invalid(_))
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                refused.apply_frictional_impacts(&[c], controls, cancelled),
                Err(SolidFrictionalImpactError::Cancelled)
            ));
        });
        assert_eq!(bits(&refused), before);
        // Only an explicit change of the physical contact allowances admits
        // this known rounded basis; the speed-readback allowance alone did not.
        controls = loose_readback;
        controls.contact.velocity_tolerance = 1e-3;
        controls.contact.complementarity_tolerance = 1e-4;
        refused
            .apply_frictional_impacts(&[c], controls, cx)
            .unwrap();
        reference
            .apply_frictional_impacts(&[c], controls, cx)
            .unwrap();
        assert_eq!(bits(&refused), bits(&reference));
        refused.advance(&[0.0; 3], cx).unwrap();
        reference.advance(&[0.0; 3], cx).unwrap();
        assert_eq!(bits(&refused), bits(&reference));
    });
}
