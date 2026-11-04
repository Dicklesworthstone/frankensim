//! G0/G3 original quadratic force, pressure projection and accepted-history checks.

use super::*;
use fs_solid::hyper3::mixed::MixedUpError;
use fs_solid::hyper3::mixed::visco::dynamics::{
    FiniteViscoContactDriftControls, SolidImpulseControls,
};
use fs_solver::CheckedCgConfig;

fn kick_controls() -> SolidImpulseControls {
    SolidImpulseControls {
        solve: CheckedCgConfig {
            tolerance: 1e-13,
            max_iterations: 64,
            max_corrections: 2,
        },
        maximum_workspace_components: 4096,
        linear_tolerance_kg_m_s: 2e-10,
        angular_tolerance_kg_m2_s: 2e-10,
        work_tolerance_j: 2e-10,
    }
}

fn drift_controls() -> FiniteViscoContactDriftControls {
    FiniteViscoContactDriftControls {
        material_balance_tolerance_j: 2e-9,
        maximum_integration_defect_j: 20.0,
    }
}

fn near_vector(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    for axis in 0..3 {
        near(actual[axis], expected[axis], tolerance);
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn shared_material(before: &MixedViscoTetMotion<'_>, after: &MixedViscoTetMotion<'_>) {
    assert!(core::ptr::eq(
        before.solid().displacement_m().as_ptr(),
        after.solid().displacement_m().as_ptr()
    ));
    assert!(core::ptr::eq(
        before.solid().pressure_pa().as_ptr(),
        after.solid().pressure_pa().as_ptr()
    ));
    assert!(core::ptr::eq(
        before.solid().quadrature_states().as_ptr(),
        after.solid().quadrature_states().as_ptr()
    ));
    assert!(core::ptr::eq(
        before.solid().first_piola_pa().as_ptr(),
        after.solid().first_piola_pa().as_ptr()
    ));
    assert_eq!(before.solid().time_s(), after.solid().time_s());
}

#[test]
fn curved_heterogeneous_supports_keep_signed_loads_actual_torque_and_original_histories() {
    with_cx(|cx| {
        let (mut nodes, elements) = mesh(cx);
        let curvature = 0.2;
        let shear = 0.08;
        // Original reference x = X + c X², current x = X + c X² + s YZ.
        // det(reference J) = 1 + 2cX; the current shear has det F = 1.
        for node in &mut nodes {
            node[0] += curvature * node[0] * node[0];
        }
        let laws: Vec<_> = (0..64)
            .map(|q| {
                ViscoHyperelastic::new(
                    Hyperelastic::new(
                        HyperelasticModel::MooneyRivlin {
                            c10: 6.0 + 0.1 * f64::from(q % 3),
                            c01: 0.0,
                            kappa: 120.0 + f64::from(q),
                        },
                        4.0,
                    )
                    .unwrap(),
                    vec![FiniteStrainBranch {
                        energy: BranchEnergy::NeoHookean {
                            mu_pa: 36.0 + 0.2 * f64::from(q),
                        },
                        tau_s: 0.04 + 0.001 * f64::from(q),
                    }],
                )
                .unwrap()
            })
            .collect();
        let mut state = MixedViscoTetSolid::try_new(
            &nodes,
            &elements,
            MixedViscoLawField::PerQuadrature(&laws),
            VolumetricTreatment::OriginalMaterial,
            TetAssemblyBudget::standard(),
            4,
            64,
            cx,
        )
        .unwrap();
        let prescribed: Vec<_> = nodes
            .iter()
            .flat_map(|x| [shear * x[1] * x[2], 0.0, 0.0])
            .enumerate()
            .collect();
        for dt in [0.0, 0.007] {
            state = state
                .advance(
                    dt,
                    &prescribed,
                    &[],
                    [0.0; 3],
                    None,
                    settings(),
                    controls().material,
                    cx,
                )
                .unwrap()
                .solid;
        }
        assert!(state.quadrature_states().iter().any(|q| q.dissipated > 0.0));
        let motion = MixedViscoTetMotion::try_new(
            state,
            &[DENSITY],
            &vec![[0.0; 3]; nodes.len()],
            MASS_BYTES,
            cx,
        )
        .unwrap();
        let original = motion.solid().quadrature_states().to_vec();
        let held: Vec<_> = (0..3 * nodes.len()).collect();
        let dt = 0.01;
        let origin = [-0.25, 0.3, -0.1];
        let body = [3.0, -4.0, 2.0];
        let loaded_dof = 3 * elements[0][4];
        // Compensation must retain this small net load at the original midside.
        let nodal = [(loaded_dof, 1e16), (loaded_dof, 1.0), (loaded_dof, -1e16)];
        let (_, internal_only) = motion
            .force_kick_candidate(dt, &held, &[], [0.0; 3], origin, kick_controls(), cx)
            .unwrap();
        let (kicked, kick) = motion
            .force_kick_candidate(dt, &held, &nodal, body, origin, kick_controls(), cx)
            .unwrap();
        shared_material(&motion, &kicked);
        assert_eq!(kicked.velocity_mps(), motion.velocity_mps());
        let volume = VOLUME * (1.0 + curvature / 2.0);
        let centroid = [
            (0.25
                + 0.3 * curvature
                + 0.1 * curvature * curvature
                + shear * (1.0 / 20.0 + curvature / 60.0))
                / (1.0 + curvature / 2.0),
            (0.25 + 0.1 * curvature) / (1.0 + curvature / 2.0),
            (0.25 + 0.1 * curvature) / (1.0 + curvature / 2.0),
        ];
        let source_linear = core::array::from_fn(|axis| {
            dt * (volume * body[axis] + if axis == 0 { 1.0 } else { 0.0 })
        });
        let body_moment = cross(core::array::from_fn(|a| centroid[a] - origin[a]), body);
        let loaded = nodes[elements[0][4]];
        let nodal_moment = cross(
            core::array::from_fn(|a| loaded[a] - origin[a]),
            [1.0, 0.0, 0.0],
        );
        let source_angular =
            core::array::from_fn(|a| dt * (volume * body_moment[a] + nodal_moment[a]));
        near_vector(kick.source_impulse.linear_kg_m_s, source_linear, 3e-13);
        near_vector(kick.source_impulse.angular_kg_m2_s, source_angular, 3e-13);
        near_vector(kick.internal_impulse.linear_kg_m_s, [0.0; 3], 3e-12);
        near_vector(kick.internal_impulse.angular_kg_m2_s, [0.0; 3], 3e-12);
        near_vector(
            kick.impulse.support.linear_kg_m_s,
            source_linear.map(|v| -v),
            3e-12,
        );
        near_vector(
            kick.impulse.support.angular_kg_m2_s,
            source_angular.map(|v| -v),
            3e-12,
        );
        near(kick.impulse.support_work_j, 0.0, 0.0);
        // Independently integrated quadratic shapes times (1+2cX), including
        // the NEGATIVE corner resultants. Edges are 01,12,02,03,13,23.
        let edges = [[0, 1], [1, 2], [0, 2], [0, 3], [1, 3], [2, 3]];
        for (local, &node) in elements[0].iter().enumerate() {
            let weight = VOLUME
                * if local < 4 {
                    -1.0 / 20.0 - if local == 1 { 0.0 } else { curvature / 30.0 }
                } else {
                    1.0 / 5.0
                        + if edges[local - 4].contains(&1) {
                            2.0 * curvature / 15.0
                        } else {
                            curvature / 15.0
                        }
                };
            for axis in 0..3 {
                let dof = 3 * node + axis;
                let delta = kick.impulse.support_impulses_kg_m_s[dof].1
                    - internal_only.impulse.support_impulses_kg_m_s[dof].1;
                let load = weight * body[axis] + if dof == loaded_dof { 1.0 } else { 0.0 };
                near(delta, -dt * load, 3e-13);
            }
        }
        let (drifted, drift) = kicked
            .contact_drift_candidate(&kick, settings(), drift_controls(), cx)
            .unwrap();
        assert_eq!(drifted.velocity_mps(), motion.velocity_mps());
        assert_eq!(
            drifted.solid().displacement_m(),
            motion.solid().displacement_m()
        );
        for (index, (before, after)) in original
            .iter()
            .zip(drifted.solid().quadrature_states())
            .enumerate()
        {
            let expected = laws[index]
                .step_remainder_response(before, &before.deformation_gradient, dt)
                .unwrap();
            assert_eq!(*after, expected.state);
        }
        assert!(drift.dissipated_j > 0.0);
        near(drift.source_displacement_work_j, 0.0, 0.0);
        near(drift.kinetic_change_j, 0.0, 0.0);
        near(
            drift.stored_after_j - drift.stored_before_j + drift.dissipated_j,
            0.0,
            2e-10,
        );
        near_vector(drift.linear_residual_kg_m_s, [0.0; 3], 0.0);
        near_vector(drift.angular_residual_kg_m2_s, [0.0; 3], 0.0);
        assert_eq!(motion.solid().quadrature_states(), original);
        assert_eq!(motion.solid().time_s(), 0.007);
        near(drifted.solid().time_s(), 0.017, 3e-18);
    });
}

#[test]
fn pressure_only_drift_retains_original_weak_volume_and_bulk_storage() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let bulk = 120.0;
        let material = law(bulk, true);
        let alpha = 0.08;
        let dt = 0.02;
        let velocities: Vec<_> = nodes
            .iter()
            .map(|x| x.map(|v| alpha / dt * v * v))
            .collect();
        let motion = MixedViscoTetMotion::try_new(
            solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::OriginalMaterial,
                cx,
            ),
            &[DENSITY],
            &velocities,
            MASS_BYTES,
            cx,
        )
        .unwrap();
        let (kicked, kick) = motion
            .force_kick_candidate(dt, &[], &[], [0.0; 3], [0.0; 3], kick_controls(), cx)
            .unwrap();
        shared_material(&motion, &kicked);
        let (after, drift) = kicked
            .contact_drift_candidate(&kick, settings(), drift_controls(), cx)
            .unwrap();
        // x_i' = X_i + alpha X_i², so J = product(1+2 alpha X_i).
        // Exact simplex moments give its P1 L2 projection J(p)-1 = A+B(X+Y+Z).
        // This has negative p at corner 0 even though J(F)>=1 everywhere.
        let a = -2.0 * alpha * alpha / 5.0 - 8.0 * alpha.powi(3) / 105.0;
        let b = 2.0 * alpha + 4.0 * alpha * alpha / 3.0 + 4.0 * alpha.powi(3) / 21.0;
        for (&node, &pressure) in after
            .solid()
            .pressure_nodes()
            .iter()
            .zip(after.solid().pressure_pa())
        {
            let expected = bulk * (a + b * nodes[node].iter().sum::<f64>());
            near(pressure, expected, 2e-8);
        }
        assert!(after.solid().pressure_pa()[0] < -0.3);
        let projected_bulk = 0.5 * bulk * VOLUME * (a * a + 1.5 * a * b + 0.6 * b * b);
        near(after.solid().bulk_energy_j(), projected_bulk, 2e-10);
        // The pointwise original bulk at actual J is a DIFFERENT model here.
        let naive_bulk = 0.5
            * bulk
            * VOLUME
            * (12.0 * alpha.powi(2) / 5.0
                + 2.0 * alpha.powi(3)
                + 24.0 * alpha.powi(4) / 35.0
                + 4.0 * alpha.powi(5) / 35.0
                + 8.0 * alpha.powi(6) / 945.0);
        assert!(naive_bulk - projected_bulk > 1e-6);
        assert_eq!(after.velocity_mps(), kicked.velocity_mps());
        near(drift.kinetic_change_j, 0.0, 0.0);
        let mass = DENSITY * VOLUME;
        near(
            drift.kinetic_after_j,
            3.0 * mass * (alpha / dt).powi(2) / 70.0,
            2e-11,
        );
        near_vector(
            drift.after.linear_kg_m_s,
            [mass * alpha / (10.0 * dt); 3],
            2e-11,
        );
        near_vector(drift.angular_residual_kg_m2_s, [0.0; 3], 2e-11);
        assert!(kick.impulse.support_impulses_kg_m_s.is_empty());
        near(kick.impulse.support_work_j, 0.0, 0.0);
        near(drift.source_displacement_work_j, 0.0, 0.0);
        near(
            kick.impulse.kinetic_change_j + drift.stored_after_j - drift.stored_before_j
                + drift.dissipated_j
                - drift.integration_defect_j,
            0.0,
            2e-9,
        );
        near(
            drift.native_loading_work_j
                - (drift.stored_after_j - drift.stored_before_j)
                - drift.dissipated_j,
            drift.material_balance_error_j,
            2e-10,
        );
        assert!(drift.absolute_integration_defect_j > 0.0);
        assert_eq!(motion.solid().time_s(), 0.0);
        assert!(motion.solid().pressure_pa().iter().all(|p| *p == 0.0));
    });
}

#[test]
fn native_phase_binding_survives_contact_and_ignores_mutable_diagnostics() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(120.0, true);
        let motion = MixedViscoTetMotion::try_new(
            solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::OriginalMaterial,
                cx,
            ),
            &[DENSITY],
            &vec![[0.0; 3]; nodes.len()],
            MASS_BYTES,
            cx,
        )
        .unwrap();
        let (kicked, mut kick) = motion
            .force_kick_candidate(
                0.01,
                &[],
                &[],
                [6.0, 0.0, 0.0],
                [0.0; 3],
                kick_controls(),
                cx,
            )
            .unwrap();
        assert!(matches!(
            motion.contact_drift_candidate(&kick, settings(), drift_controls(), cx),
            Err(MixedViscoTetError::Native(
                MixedUpError::InvalidInput { .. }
            ))
        ));
        let contacted = kicked
            .impulse_candidate(
                &[],
                &[(3 * elements[0][4] + 1, 1e-6)],
                [0.0; 3],
                kick_controls(),
                cx,
            )
            .unwrap()
            .motion;
        shared_material(&motion, &contacted);
        let (expected, expected_drift) = contacted
            .contact_drift_candidate(&kick, settings(), drift_controls(), cx)
            .unwrap();
        kick.duration_s = 123.0;
        kick.source_impulse_work_j = f64::NAN;
        kick.internal_impulse_work_j = f64::INFINITY;
        kick.stored_before_j = -1e99;
        kick.impulse.work_residual_j = 1e99;
        let (actual, actual_drift) = contacted
            .contact_drift_candidate(&kick, settings(), drift_controls(), cx)
            .unwrap();
        assert_eq!(
            actual.solid().displacement_m(),
            expected.solid().displacement_m()
        );
        assert_eq!(actual.solid().pressure_pa(), expected.solid().pressure_pa());
        assert_eq!(
            actual.solid().quadrature_states(),
            expected.solid().quadrature_states()
        );
        assert_eq!(actual_drift.duration_s, 0.01);
        assert_eq!(
            actual_drift.integration_defect_j,
            expected_drift.integration_defect_j
        );
        let (superseding, _) = contacted
            .force_kick_candidate(0.01, &[], &[], [0.0; 3], [0.0; 3], kick_controls(), cx)
            .unwrap();
        assert!(
            superseding
                .contact_drift_candidate(&kick, settings(), drift_controls(), cx)
                .is_err()
        );
        assert!(
            actual
                .contact_drift_candidate(&kick, settings(), drift_controls(), cx)
                .is_err()
        );
        let strict = FiniteViscoContactDriftControls {
            maximum_integration_defect_j: 0.0,
            ..drift_controls()
        };
        assert!(matches!(
            contacted.contact_drift_candidate(&kick, settings(), strict, cx),
            Err(MixedViscoTetError::EnergyGate { .. })
        ));
        let limited = SolidImpulseControls {
            maximum_workspace_components: 1,
            ..kick_controls()
        };
        assert!(
            motion
                .force_kick_candidate(0.01, &[], &[], [0.0; 3], [0.0; 3], limited, cx)
                .is_err()
        );
        let exact = MixedViscoTetMotion::try_new(
            solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::ExactIncompressible,
                cx,
            ),
            &[DENSITY],
            &vec![[0.0; 3]; nodes.len()],
            MASS_BYTES,
            cx,
        )
        .unwrap();
        assert!(matches!(
            exact.force_kick_candidate(0.01, &[], &[], [0.0; 3], [0.0; 3], kick_controls(), cx),
            Err(MixedViscoTetError::ExactIncompressibleSplit)
        ));
        assert_eq!(motion.solid().time_s(), 0.0);
        assert!(motion.velocity_mps().iter().all(|v| *v == [0.0; 3]));
    });
}
