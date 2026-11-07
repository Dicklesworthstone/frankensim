//! Independent physical endpoint fields on the original dynamic fatigue owner.

use super::*;
use fs_material::fatigue::critical_plane::fatemi_socie::{
    FatemiSocieCalibration, FatemiSocieLimits, FatemiSocieNormalStress, FatemiSocieParameters,
};
use fs_material::fatigue::critical_plane::miner::CriticalPlaneMinerLimits;
use fs_solid::plastic3::dynamics::fatigue::J2DynamicFatigueLoadingTarget;
use fs_solid::plastic3::fatigue::{J2FatemiSocieGauge, J2FatigueLoading, J2VariableAmplitudeGauge};

fn translation_values(fixed: &[(usize, f64)], translation: f64) -> Vec<f64> {
    fixed
        .iter()
        .map(|&(dof, _)| if dof % 3 == 0 { translation } else { 0.0 })
        .collect()
}

#[test]
fn independent_dynamic_fields_match_consistent_mass_support_coupling_and_native_energy() {
    let law = material(false, Some(Density::DIMS));
    let supports = fixed(true);
    let p = problem(&law, &supports, &[1.0]);
    let values = [0.001, -0.002, 0.0005, -0.0015]
        .map(|translation| translation_values(&supports, translation));
    let forces = [
        vec![(3, 0.003), (3, 0.002)],
        vec![(3, -0.008)],
        vec![],
        vec![(3, 0.009), (3, -0.001), (0, 0.001)],
    ];
    let times = [0.1, 0.25, 0.35, 0.55];
    let factors = [2.0, 0.5, -1.5, 0.75];
    let targets: Vec<_> = (0..times.len())
        .map(|i| J2DynamicFatigueLoadingTarget {
            time_s: times[i],
            load_factor: factors[i],
            loading: J2FatigueLoading {
                prescribed_values_m: &values[i],
                nodal_forces_n: &forces[i],
            },
        })
        .collect();
    with_cx(&CancelGate::new(), |cx| {
        let curve = curve(&law, 0.5, 1e12);
        let mut run = J2DynamicFatigueRun::new(
            &p,
            FRAME,
            vec![gauge(
                &law,
                curve.clone(),
                StrainLifeMeanStress::Uncorrected,
            )],
            config(targets.len()),
            0.0,
            vec![[0.0; 3]; 4],
            cx,
        )
        .unwrap();
        let mut batch = run.clone();
        let mut native = problem(&law, &supports, &[0.0]).solve(cx).unwrap();
        let mut velocity = vec![[0.0; 3]; 4];
        let mut counted = StrainLifeMiner::new(
            curve,
            StrainLifeMeanStress::Uncorrected,
            config(targets.len()).fatigue.history,
        )
        .unwrap();
        counted.push(0.0).unwrap();
        let (mut time, mut free_u, mut free_v, mut support_u, mut support_v) =
            (0.0, 0.0, 0.0, 0.0, 0.0);
        for target in &targets {
            let dt = target.time_s - time;
            let support = target.load_factor * target.loading.prescribed_values_m[0];
            let force: f64 = target
                .loading
                .nodal_forces_n
                .iter()
                .filter(|(dof, _)| *dof == 3)
                .map(|(_, value)| target.load_factor * value)
                .sum();
            // The unit tet has rho V=1, M_11=1/10 and three moving-support
            // couplings M_1j=1/20. With nu=0 the sole free x DOF has K=E/6=10.
            // This recurrence is independent of the native assembly/BE owner.
            let expected_u = (force
                + 10.0 * support
                + 0.1 * (free_u + dt * free_v) / (dt * dt)
                + 0.15 * (support_u + dt * support_v - support) / (dt * dt))
                / (10.0 + 0.1 / (dt * dt));
            let prescribed: Vec<_> = supports
                .iter()
                .zip(target.loading.prescribed_values_m)
                .map(|(&(dof, _), value)| (dof, target.load_factor * value))
                .collect();
            let endpoint_forces: Vec<_> = target
                .loading
                .nodal_forces_n
                .iter()
                .map(|&(dof, value)| (dof, target.load_factor * value))
                .collect();
            let endpoint = SmallStrainTetProblem {
                nodal_forces_n: &endpoint_forces,
                ..problem(&law, &prescribed, &[1.0])
            };
            let step = endpoint
                .advance_backward_euler(
                    native.displacement_m.as_flattened(),
                    &velocity,
                    &native.element_states,
                    &[6.0],
                    dt,
                    cx,
                )
                .unwrap();
            let report = run
                .advance_to_with_loading(target.time_s, target.load_factor, target.loading, cx)
                .unwrap();
            assert_eq!(run.solution(), &step.solution);
            assert_eq!(run.velocity_mps(), &step.velocity_mps);
            assert_eq!(report.energy, step.energy);
            assert_eq!(report.load_factor, target.load_factor);
            close(run.solution().displacement_m[1][0], expected_u, 1e-12);
            close(
                run.samples()[0].mechanical_strain,
                expected_u - support,
                1e-12,
            );
            close(
                run.samples()[0].normal_stress_pa,
                60.0 * (expected_u - support),
                1e-10,
            );
            counted.push(expected_u - support).unwrap();
            free_v = (expected_u - free_u) / dt;
            support_v = (support - support_u) / dt;
            free_u = expected_u;
            support_u = support;
            time = target.time_s;
            native = step.solution;
            velocity = step.velocity_mps;
        }
        batch.advance_with_loading(&targets, cx).unwrap();
        same_state(&run, &batch);
        let result = run.finish(cx).unwrap();
        assert_eq!(result, batch.finish(cx).unwrap());
        close(
            result[0].miner_damage,
            counted.finish().unwrap().miner_damage,
            1e-12,
        );
        assert!(result[0].miner_damage > 0.0);
    });
}

fn affine_values(xx: f64, xy: f64) -> Vec<f64> {
    NODES
        .iter()
        .flat_map(|x| [xx * x[0] + xy * x[1], xy * x[0], 0.0])
        .collect()
}

fn shear_gauge(law: &ResolvedJ2Plasticity, minimum_cycles: f64) -> J2FatemiSocieGauge {
    let mut scope = curve(law, 0.5, 1e12).scope().clone();
    scope.source_reference =
        "synthetic independent torsional fatigue reference; not measured".into();
    scope.load_basis = "total engineering shear strain amplitude; full cycles".into();
    let shear_modulus = law.law().elastic.lame().1;
    J2FatemiSocieGauge {
        element: 0,
        material_state_identity: law.material().identity(),
        calibration: Arc::new(
            FatemiSocieCalibration::new(
                scope,
                FatemiSocieParameters {
                    shear_modulus_pa: shear_modulus,
                    shear_fatigue_strength_coefficient_pa: 0.1 * shear_modulus,
                    shear_fatigue_strength_exponent: -0.5,
                    shear_fatigue_ductility_coefficient: 0.1,
                    shear_fatigue_ductility_exponent: -0.5,
                    normal_stress_sensitivity: 0.5,
                    initial_yield_stress_pa: law.law().yield_stress,
                    normal_stress: FatemiSocieNormalStress::TensionOnly,
                },
                minimum_cycles,
                1e9,
            )
            .unwrap(),
        ),
        normals_in_mesh_frame: vec![[1.0, 0.0, 0.0]],
        limits: FatemiSocieLimits {
            maximum_samples: 32,
            maximum_planes: 1,
            maximum_projection_work: 32,
            maximum_circle_work: 100_000,
        },
        closure: CycleClosureTolerance {
            absolute_strain: 1e-12,
            absolute_stress_pa: 1e-10,
            relative: 1e-10,
        },
    }
}

fn same_banks(a: &J2DynamicFatigueRun<'_, '_>, b: &J2DynamicFatigueRun<'_, '_>, cx: &Cx<'_>) {
    same_state(a, b);
    assert_eq!(a.fatemi_socie_results(), b.fatemi_socie_results());
    assert_eq!(
        a.variable_amplitude_results(),
        b.variable_amplitude_results()
    );
    assert_eq!(
        a.preview_fatigue(cx).unwrap(),
        b.preview_fatigue(cx).unwrap()
    );
    assert_eq!(
        a.preview_variable_amplitude_fatigue(cx).unwrap(),
        b.preview_variable_amplitude_fatigue(cx).unwrap()
    );
}

#[test]
fn explicit_block_refusals_restore_native_motion_and_all_four_fatigue_banks_before_retry() {
    let law = material(false, Some(Density::DIMS));
    let supports = fixed(false);
    let p = problem(&law, &supports, &[1.0]);
    let mut settings = config(8);
    settings.fatigue.maximum_gauges = 5;
    let initial_values = affine_values(0.0001, 0.0001);
    let values = [
        affine_values(0.002, 0.001),
        affine_values(-0.001, -0.001),
        affine_values(0.0, 0.0),
    ];
    let targets: Vec<_> = values
        .iter()
        .enumerate()
        .map(|(i, values)| J2DynamicFatigueLoadingTarget {
            time_s: 0.1 * (i + 1) as f64,
            load_factor: 1.0,
            loading: J2FatigueLoading {
                prescribed_values_m: values,
                nodal_forces_n: &[],
            },
        })
        .collect();
    with_cx(&CancelGate::new(), |cx| {
        let mut run = J2DynamicFatigueRun::new_with_variable_amplitude(
            &p,
            FRAME,
            vec![gauge(
                &law,
                curve(&law, 0.5, 1e12),
                StrainLifeMeanStress::Uncorrected,
            )],
            vec![critical_gauge(&law, 0, vec![[1.0, 0.0, 0.0]])],
            vec![shear_gauge(&law, 0.5), shear_gauge(&law, 10_000.0)],
            vec![J2VariableAmplitudeGauge {
                element: 0,
                material_state_identity: law.material().identity(),
                curve: curve(&law, 0.5, 1e12),
                mean_stress: StrainLifeMeanStress::Uncorrected,
                normals_in_mesh_frame: vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
                limits: CriticalPlaneMinerLimits {
                    rainflow: settings.fatigue.history,
                    maximum_planes: 2,
                    maximum_projection_work: 64,
                },
            }],
            settings,
            0.0,
            vec![[0.0; 3]; 4],
            cx,
        )
        .unwrap();
        run.advance_to_with_loading(
            0.05,
            1.0,
            J2FatigueLoading {
                prescribed_values_m: &initial_values,
                nodal_forces_n: &[],
            },
            cx,
        )
        .unwrap();
        let mut expected = run.clone();

        // The invalid second field must refuse the complete proposed block.
        let outside = [(12, 1.0)];
        let mut invalid = [targets[0], targets[1]];
        invalid[1].loading.nodal_forces_n = &outside;
        assert!(matches!(
            run.advance_with_loading(&invalid, cx),
            Err(J2FatigueError::Invalid(_))
        ));
        same_banks(&run, &expected, cx);
        let overflowing = vec![f64::MAX; supports.len()];
        invalid[1].loading = J2FatigueLoading {
            prescribed_values_m: &overflowing,
            nodal_forces_n: &[],
        };
        invalid[1].load_factor = 2.0;
        assert!(matches!(
            run.advance_with_loading(&invalid, cx),
            Err(J2FatigueError::Invalid(_))
        ));
        same_banks(&run, &expected, cx);

        // All samples reach native mechanics before closure rejects this prefix.
        assert!(matches!(
            run.advance_cycle_with_loading(&targets[..2], cx),
            Err(J2FatigueError::CriticalPlane {
                error: CriticalPlaneError::NotClosed { .. },
                ..
            })
        ));
        same_banks(&run, &expected, cx);

        // The complete trajectory closes SWT and shear gauge zero, then shear
        // gauge one rejects the life range. No completed result or sample leaks.
        assert!(matches!(
            run.advance_cycle_with_loading(&targets, cx),
            Err(J2FatigueError::FatemiSocie { gauge: 1, .. })
        ));
        same_banks(&run, &expected, cx);

        let smaller_values: Vec<Vec<_>> = values
            .iter()
            .map(|values| values.iter().map(|value| 0.25 * value).collect())
            .collect();
        let retry: Vec<_> = smaller_values
            .iter()
            .zip(&targets)
            .map(|(values, target)| J2DynamicFatigueLoadingTarget {
                loading: J2FatigueLoading {
                    prescribed_values_m: values,
                    nodal_forces_n: &[],
                },
                ..*target
            })
            .collect();
        run.advance_cycle_with_loading(&retry, cx).unwrap();
        expected.advance_with_loading(&retry[..2], cx).unwrap();
        expected
            .advance_cycle_with_loading(&retry[2..], cx)
            .unwrap();
        same_banks(&run, &expected, cx);
        assert_eq!(run.accepted_steps(), 4);
        assert_eq!(run.summaries()[0].samples, 5);
        assert!(
            run.fatemi_socie_results()
                .iter()
                .all(|result| result.is_some_and(|result| result.samples == 5))
        );
        assert_eq!(run.finish(cx).unwrap(), expected.finish(cx).unwrap());
    });
}
