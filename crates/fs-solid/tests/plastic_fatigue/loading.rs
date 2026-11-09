//! Independent physical boundary fields enter every original fatigue bank.
//! Synthetic constitutive and life data are implementation references only.

use super::*;
use fs_material::fatigue::critical_plane::fatemi_socie::{
    FatemiSocieCalibration, FatemiSocieCycle, FatemiSocieError, FatemiSocieLimits,
    FatemiSocieNormalStress, FatemiSocieParameters,
};
use fs_material::fatigue::critical_plane::miner::{CriticalPlaneMiner, CriticalPlaneMinerLimits};
use fs_material::fatigue::critical_plane::{
    CriticalPlaneError, CriticalPlaneLimits, CycleClosureTolerance, TensorFatigueSample,
};
use fs_material::fatigue::cyclic::SourceBoundFatigueLaw;
use fs_material::plastic::J2State;
use fs_solid::plastic3::SmallStrainTetSolution;
use fs_solid::plastic3::fatigue::{
    J2CriticalPlaneGauge, J2FatemiSocieGauge, J2FatigueLoading, J2VariableAmplitudeGauge,
};

const FRAME: ContentHash = ContentHash([73; 32]);
const R: f64 = 2e-4;

fn closure() -> CycleClosureTolerance {
    CycleClosureTolerance {
        absolute_strain: 1e-12,
        absolute_stress_pa: 1e-10,
        relative: 1e-10,
    }
}

fn shear_gauge(law: &ResolvedJ2Plasticity) -> J2FatemiSocieGauge {
    let mut scope = curve(300.0, 1e12).scope().clone();
    scope.source_reference = "synthetic torsional equal-exponent FS reference; not measured".into();
    scope.surface_and_geometry = "synthetic smooth torsion coupon".into();
    scope.load_basis =
        "engineering shear strain amplitude; independently declared torsional fit".into();
    let g = law.law().elastic.lame().1;
    J2FatemiSocieGauge {
        element: 0,
        material_state_identity: law.material().identity(),
        calibration: Arc::new(
            FatemiSocieCalibration::new(
                scope,
                FatemiSocieParameters {
                    shear_modulus_pa: g,
                    shear_fatigue_strength_coefficient_pa: 0.1 * g,
                    shear_fatigue_strength_exponent: -0.5,
                    shear_fatigue_ductility_coefficient: 0.1,
                    shear_fatigue_ductility_exponent: -0.5,
                    normal_stress_sensitivity: 0.2,
                    initial_yield_stress_pa: law.law().yield_stress,
                    normal_stress: FatemiSocieNormalStress::TensionOnly,
                },
                1.0,
                1e8,
            )
            .unwrap(),
        ),
        normals_in_mesh_frame: vec![[0.0, 0.0, 1.0]],
        limits: FatemiSocieLimits {
            maximum_samples: 32,
            maximum_planes: 4,
            maximum_projection_work: 128,
            maximum_circle_work: 4096,
        },
        closure: closure(),
    }
}

fn variable_gauge(law: &ResolvedJ2Plasticity) -> J2VariableAmplitudeGauge {
    J2VariableAmplitudeGauge {
        element: 0,
        material_state_identity: law.material().identity(),
        curve: curve(300.0, 1e12),
        mean_stress: StrainLifeMeanStress::Uncorrected,
        normals_in_mesh_frame: vec![[1.0, 0.0, 1.0], [0.0, 1.0, 1.0]],
        limits: CriticalPlaneMinerLimits {
            rainflow: config().history,
            maximum_planes: 4,
            maximum_projection_work: 128,
        },
    }
}

fn banked_path<'p, 'm>(
    problem: &'p SmallStrainTetProblem<'m, ResolvedJ2Plasticity>,
    law: &ResolvedJ2Plasticity,
    shear: Vec<J2FatemiSocieGauge>,
    cx: &Cx<'_>,
) -> J2FatiguePath<'p, 'm> {
    let mut fixed = gauge(law, curve(300.0, 1e12), 0);
    fixed.normal_in_mesh_frame = [1.0, 0.0, 1.0];
    let critical = J2CriticalPlaneGauge {
        element: 0,
        material_state_identity: law.material().identity(),
        curve: curve(300.0, 1e12),
        normals_in_mesh_frame: vec![[1.0, 0.0, 1.0]],
        limits: CriticalPlaneLimits {
            maximum_samples: 32,
            maximum_planes: 4,
            maximum_projection_work: 128,
        },
        closure: closure(),
    };
    J2FatiguePath::new_with_variable_amplitude(
        problem,
        FRAME,
        vec![fixed],
        vec![critical],
        shear,
        vec![variable_gauge(law)],
        config(),
        cx,
    )
    .unwrap()
}

/// Symmetric displacement gradient with eps_xz=u, eps_yz=v and zero rigid spin.
fn shear_values(u: f64, v: f64) -> Vec<f64> {
    NODES
        .iter()
        .flat_map(|x| [u * x[2], v * x[2], u * x[0] + v * x[1]])
        .collect()
}

fn loading(values: &[f64]) -> J2FatigueLoading<'_> {
    J2FatigueLoading {
        prescribed_values_m: values,
        nodal_forces_n: &[],
    }
}

fn sample(solution: &SmallStrainTetSolution<J2State>) -> TensorFatigueSample {
    TensorFatigueSample {
        stress_pa: solution.element_stress[0],
        total_mechanical_strain: solution.element_strain[0],
    }
}

fn unchanged<L: SourceBoundFatigueLaw>(
    a: &J2FatiguePath<'_, '_, L>,
    b: &J2FatiguePath<'_, '_, L>,
    cx: &Cx<'_>,
) {
    assert_eq!(a.solution(), b.solution());
    assert_eq!(a.accepted_increments(), b.accepted_increments());
    assert_eq!(a.remaining_increments(), b.remaining_increments());
    assert_eq!(a.samples(), b.samples());
    assert_eq!(a.summaries(), b.summaries());
    assert_eq!(a.critical_plane_results(), b.critical_plane_results());
    assert_eq!(a.fatemi_socie_results(), b.fatemi_socie_results());
    assert_eq!(
        a.fatemi_socie_miner_summaries(),
        b.fatemi_socie_miner_summaries()
    );
    for gauge in 0..a.fatemi_socie_results().len() {
        assert_eq!(
            a.fatemi_socie_plane_damage(gauge),
            b.fatemi_socie_plane_damage(gauge)
        );
    }
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
fn independently_phased_native_shear_triangle_retains_vector_geometry_and_every_bank() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let law = material();
        let fixed: Vec<_> = (0..3 * NODES.len()).map(|d| (d, 0.0)).collect();
        let factors = [1.0; 4];
        let p = problem(&law, &fixed, &factors);
        let shear = shear_gauge(&law);
        let mut run = banked_path(&p, &law, vec![shear.clone()], cx);
        let endpoints = [
            (R, 0.0),
            (-0.5 * R, 0.5 * 3.0_f64.sqrt() * R),
            (-0.5 * R, -0.5 * 3.0_f64.sqrt() * R),
            (R, 0.0),
        ];
        let values: Vec<_> = endpoints.iter().map(|&(u, v)| shear_values(u, v)).collect();
        let targets: Vec<_> = values.iter().map(|values| loading(values)).collect();
        let mut native = problem(&law, &fixed, &[0.0]).solve(cx).unwrap();
        let var = variable_gauge(&law);
        let mut direct_variable = CriticalPlaneMiner::new(
            var.curve,
            var.mean_stress,
            var.normals_in_mesh_frame,
            var.limits,
        )
        .unwrap();
        direct_variable.push(sample(&native)).unwrap();
        let mut direct_shear = FatemiSocieCycle::new(
            shear.calibration,
            shear.normals_in_mesh_frame,
            shear.limits,
            shear.closure,
        )
        .unwrap();
        for (index, (&(u, v), values)) in endpoints.iter().zip(&values).enumerate() {
            let prescribed: Vec<_> = fixed
                .iter()
                .zip(values)
                .map(|(&(d, _), &value)| (d, value))
                .collect();
            native = problem(&law, &prescribed, &[1.0])
                .solve_from(
                    native.displacement_m.as_flattened(),
                    &native.element_states,
                    cx,
                )
                .unwrap();
            let expected = [0.0, 0.0, 0.0, 0.0, v, u];
            for (strain, stress) in native.element_strain.iter().zip(&native.element_stress) {
                for component in 0..6 {
                    close(strain[component], expected[component], 1e-14);
                    close(
                        stress[component],
                        2.0 * law.law().elastic.lame().1 * expected[component],
                        1e-10,
                    );
                }
            }
            assert!(native.element_states.iter().all(|state| state.alpha == 0.0));
            direct_variable.push(sample(&native)).unwrap();
            direct_shear.push(sample(&native)).unwrap();
            if index == 0 {
                run.advance_with_loading(&targets[..1], cx).unwrap();
                assert_eq!(run.solution(), &native);
                run.restart_critical_plane_cycle(cx).unwrap();
            }
        }
        run.advance_cycle_with_loading(&targets[1..], cx).unwrap();
        assert_eq!(run.solution(), &native);
        let expected = direct_shear.finish().unwrap();
        assert_eq!(run.fatemi_socie_results(), vec![Some(expected)]);
        // On plane z the engineering shear vector is 2*(u,v). The three
        // vertices form an equilateral triangle: its circle radius is 2R,
        // exceeding the largest half-chord sqrt(3)R. Scalar ranges cannot
        // recover this amplitude; tau_f/G=gamma_f=.1 gives N=125000.
        close(expected.engineering_shear_strain_amplitude, 0.0004, 1e-14);
        assert!(expected.engineering_shear_strain_amplitude > 3.0_f64.sqrt() * R);
        close(expected.maximum_normal_stress_pa, 0.0, 1e-10);
        close(expected.cycles_to_failure.unwrap(), 125000.0, 1e-5);
        assert_eq!(expected.samples, 4);
        assert!(
            run.critical_plane_results()[0]
                .unwrap()
                .cycles_to_failure
                .is_some()
        );
        assert_eq!(run.summaries()[0].samples, 5);
        assert!(run.finish(cx).unwrap()[0].miner_damage > 0.0);
        assert_eq!(
            run.variable_amplitude_results(),
            vec![Some(direct_variable.finish().unwrap())]
        );
    });
}

#[test]
fn independent_supports_and_duplicate_free_forces_match_native_plastic_continuation() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let law = material();
        let fixed = prescribed();
        let factors = [0.5, 1.0, -0.75, 1.25, 0.001];
        let axial = [0.006, -0.004, 0.008, 0.003];
        let forces: Vec<_> = [0.005, -0.012, 0.008, -0.009]
            .into_iter()
            .map(|force| vec![(7, force), (7, 0.25 * force), (19, -0.5 * force)])
            .collect();
        let values: Vec<Vec<_>> = axial
            .iter()
            .map(|&strain| {
                fixed
                    .iter()
                    .map(|&(d, _)| {
                        if d % 3 == 0 {
                            strain * NODES[d / 3][0]
                        } else {
                            0.0
                        }
                    })
                    .collect()
            })
            .collect();
        let targets: Vec<_> = values
            .iter()
            .zip(&forces)
            .map(|(values, forces)| J2FatigueLoading {
                prescribed_values_m: values,
                nodal_forces_n: forces,
            })
            .collect();
        let p = problem(&law, &fixed, &factors);
        let mut whole = J2FatiguePath::new(
            &p,
            FRAME,
            vec![gauge(&law, curve(300.0, 1e12), 0)],
            config(),
            cx,
        )
        .unwrap();
        let mut split = whole.clone();
        whole.advance_with_loading(&targets, cx).unwrap();
        let mut native = problem(&law, &fixed, &[0.0]).solve(cx).unwrap();
        for (index, &factor) in factors[..targets.len()].iter().enumerate() {
            let prescribed: Vec<_> = fixed
                .iter()
                .zip(&values[index])
                .map(|(&(d, _), &value)| (d, value))
                .collect();
            let one_factor = [factor];
            let mut endpoint = problem(&law, &prescribed, &one_factor);
            endpoint.nodal_forces_n = &forces[index];
            native = endpoint
                .solve_from(
                    native.displacement_m.as_flattened(),
                    &native.element_states,
                    cx,
                )
                .unwrap();
            split
                .advance_with_loading(&targets[index..index + 1], cx)
                .unwrap();
            assert_eq!(split.solution(), &native);
            close(
                split.samples()[0].mechanical_strain,
                factor * axial[index],
                1e-12,
            );
            close(
                split.solution().displacement_m[1][0],
                2.0 * factor * axial[index],
                1e-12,
            );
        }
        assert!(native.element_states.iter().all(|state| state.alpha > 0.0));
        unchanged(&whole, &split, cx);
        // The original path remains the reference owner. Ordinary advance
        // restores its support pattern and zero forces, not the last custom fields.
        native = problem(&law, &fixed, &factors[targets.len()..])
            .solve_from(
                native.displacement_m.as_flattened(),
                &native.element_states,
                cx,
            )
            .unwrap();
        whole.advance(1, cx).unwrap();
        split.advance(1, cx).unwrap();
        assert_eq!(split.solution(), &native);
        close(split.samples()[0].mechanical_strain, 0.001, 1e-12);
        unchanged(&whole, &split, cx);
        assert_eq!(whole.finish(cx).unwrap(), split.finish(cx).unwrap());
    });
}

#[test]
fn late_loading_native_and_cycle_failures_restore_open_histories_before_retry() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let law = material();
        let fixed: Vec<_> = (0..3 * NODES.len()).map(|d| (d, 0.0)).collect();
        let p = problem(&law, &fixed, &[1.0; 4]);
        let wide = shear_gauge(&law);
        let mut narrow = wide.clone();
        narrow.calibration = Arc::new(
            FatemiSocieCalibration::new(
                wide.calibration.scope().clone(),
                wide.calibration.parameters(),
                200000.0,
                1e8,
            )
            .unwrap(),
        );
        let mut run = banked_path(&p, &law, vec![wide, narrow], cx);
        let start = shear_values(R, 0.0);
        let second = shear_values(-0.5 * R, 0.5 * 3.0_f64.sqrt() * R);
        let third = shear_values(-0.5 * R, -0.5 * 3.0_f64.sqrt() * R);
        run.advance_with_loading(&[loading(&start)], cx).unwrap();
        run.restart_critical_plane_cycle(cx).unwrap();
        run.enable_fatemi_socie_miner(1, cx).unwrap();
        let prefix = run.clone();
        assert!(matches!(
            run.advance_with_loading(&[loading(&second), loading(&third[..third.len() - 1])], cx),
            Err(J2FatigueError::Invalid(_))
        ));
        unchanged(&run, &prefix, cx);
        let overflowing_forces = [(0, f64::MAX), (0, f64::MAX)];
        let late_native_failure = J2FatigueLoading {
            prescribed_values_m: &third,
            nodal_forces_n: &overflowing_forces,
        };
        // Each force is finite at factor one. Only the native additive force
        // assembly fails, after the first candidate and its observations exist.
        assert!(matches!(
            run.advance_with_loading(&[loading(&second), late_native_failure], cx),
            Err(J2FatigueError::Mechanics(_))
        ));
        unchanged(&run, &prefix, cx);
        assert!(matches!(
            run.advance_cycle_with_loading(&[loading(&second), loading(&third)], cx),
            Err(J2FatigueError::CriticalPlane {
                error: CriticalPlaneError::NotClosed { .. },
                ..
            }) | Err(J2FatigueError::FatemiSocie {
                error: FatemiSocieError::NotClosed { .. },
                ..
            })
        ));
        unchanged(&run, &prefix, cx);
        // SWT and the first shear gauge close successfully before the second
        // shear calibration refuses N=125000 below its declared life interval.
        assert!(matches!(
            run.advance_cycle_with_loading(
                &[loading(&second), loading(&third), loading(&start)],
                cx
            ),
            Err(J2FatigueError::FatemiSocie { gauge: 1, .. })
        ));
        unchanged(&run, &prefix, cx);
        assert_eq!(run.critical_plane_results(), vec![None]);
        assert_eq!(run.fatemi_socie_results(), vec![None, None]);
        let smaller_second = shear_values(0.25 * R, 0.25 * 3.0_f64.sqrt() * R);
        let smaller_third = shear_values(0.25 * R, -0.25 * 3.0_f64.sqrt() * R);
        let retry = [
            loading(&smaller_second),
            loading(&smaller_third),
            loading(&start),
        ];
        let mut clean = prefix;
        clean.advance_with_loading(&retry[..2], cx).unwrap();
        clean.advance_cycle_with_loading(&retry[2..], cx).unwrap();
        run.advance_cycle_with_loading(&retry, cx).unwrap();
        unchanged(&run, &clean, cx);
        for result in run.fatemi_socie_results().into_iter().map(Option::unwrap) {
            assert_eq!(result.samples, 4);
            close(result.engineering_shear_strain_amplitude, R, 1e-14);
            close(result.cycles_to_failure.unwrap(), 500000.0, 1e-4);
        }
        assert_eq!(run.finish(cx).unwrap(), clean.finish(cx).unwrap());
        for summary in run
            .fatemi_socie_miner_summaries()
            .into_iter()
            .map(Option::unwrap)
        {
            assert_eq!(summary.completed_cycles, 1);
            close(summary.miner_damage, 1.0 / 500000.0, 1e-14);
        }
    });
}

#[test]
fn native_shear_miner_keeps_material_planes_across_rotating_cycles() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let law = material();
        let fixed: Vec<_> = (0..3 * NODES.len()).map(|d| (d, 0.0)).collect();
        let p = problem(&law, &fixed, &[1.0; 9]);
        let mut gauge = shear_gauge(&law);
        gauge.normals_in_mesh_frame = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mut run = J2FatiguePath::new_with_fatemi_socie(
            &p,
            FRAME,
            vec![],
            vec![],
            vec![gauge],
            config(),
            cx,
        )
        .unwrap();
        assert_eq!(run.fatemi_socie_miner_summaries(), vec![None]);
        assert!(run.enable_fatemi_socie_miner(0, cx).is_err());
        assert_eq!(run.fatemi_socie_miner_summaries(), vec![None]);
        run.enable_fatemi_socie_miner(3, cx).unwrap();
        let initial = run.fatemi_socie_miner_summaries();
        assert!(run.enable_fatemi_socie_miner(3, cx).is_err());
        assert_eq!(run.fatemi_socie_miner_summaries(), initial);

        // First xz, then weaker yz, then larger yz. Normal stress on both
        // coordinate planes is zero. The independent fit is N=.5*(.2/gamma_a)^2.
        // The two physical plane histories must never be combined into one crack.
        let mut expected = [0.0_f64; 2];
        for (cycle, (u, v, plane)) in [(R, 0.0, 0), (0.0, 0.5 * R, 1), (0.0, R, 1)]
            .into_iter()
            .enumerate()
        {
            if cycle > 0 {
                let previous = run.fatemi_socie_miner_summaries();
                run.restart_critical_plane_cycle(cx).unwrap();
                assert_eq!(run.fatemi_socie_miner_summaries(), previous);
                assert_eq!(run.fatemi_socie_results(), vec![None]);
            }
            let values = [
                shear_values(u, v),
                shear_values(-u, -v),
                shear_values(0.0, 0.0),
            ];
            let targets = values.each_ref().map(|value| loading(value));
            run.advance_cycle_with_loading(&targets, cx).unwrap();
            let gamma = 2.0 * u.abs().max(v.abs());
            expected[plane] += 2.0 * (gamma / 0.2).powi(2);
            let summary = run.fatemi_socie_miner_summaries()[0].unwrap();
            assert_eq!(summary.completed_cycles, cycle as u64 + 1);
            assert_eq!(
                summary.selected_plane,
                if expected[0] >= expected[1] { 0 } else { 1 }
            );
            close(summary.miner_damage, expected[0].max(expected[1]), 1e-14);
            assert_eq!(run.fatemi_socie_plane_damage(0).unwrap().len(), 2);
            assert!(
                run.solution()
                    .element_states
                    .iter()
                    .all(|state| state.alpha == 0.0)
            );
            if cycle == 1 {
                assert!(summary.miner_damage < expected.iter().sum::<f64>());
            }
        }
        let completed = run.fatemi_socie_miner_summaries();
        let mut replay = run.clone();
        run.finish(cx).unwrap();
        run.finish(cx).unwrap();
        replay.finish(cx).unwrap();
        assert_eq!(run.fatemi_socie_miner_summaries(), completed);
        assert_eq!(
            run.fatemi_socie_miner_summaries(),
            replay.fatemi_socie_miner_summaries()
        );
    });
}
