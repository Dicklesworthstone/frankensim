//! Full raw curvature against the independently retained first-order pipeline.
use super::*;
use fs_topopt::RobustPipeline;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::simp_second_order::{SimpComplianceGoal, SimpThermalComplianceGoal};

fn vector_close(actual: &[f64], expected: &[f64], relative: f64) {
    let error = actual
        .iter()
        .zip(expected)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    let scale = expected.iter().map(|v| v.abs()).fold(1e-8, f64::max);
    assert!(
        error <= relative * scale,
        "error={error:e}, scale={scale:e}"
    );
}

#[test]
fn ordinary_and_robust_thermal_curvature_include_both_load_and_design_map_derivatives() {
    for robust in [false, true] {
        let (mesh, positions, elasticity, thermal, a, b) = cube();
        let p = SimpParams {
            e_min: 0.03,
            penal: 2.7,
            beta: 2.3,
            eta: 0.5,
        };
        let ordinary = DesignPipeline {
            filter: DensityFilter::new(&mesh, &positions, 0.15),
            params: p,
        };
        let three_field = RobustPipeline {
            filter: DensityFilter::new(&mesh, &positions, 0.15),
            params: p,
            eta_offset: 0.08,
        };
        let reference = DesignPipeline {
            filter: DensityFilter::new(&mesh, &positions, 0.15),
            params: SimpParams {
                eta: if robust { 0.58 } else { 0.5 },
                ..p
            },
        };
        let loads = [
            ThermalLoadCase {
                force: &a,
                weight: 0.7,
                strains: TetThermalStrainField::PerElement(&thermal),
            },
            // Zero aggregate weight must not erase a separately constrained case.
            ThermalLoadCase {
                force: &b,
                weight: 0.0,
                strains: TetThermalStrainField::PerElement(&thermal),
            },
        ];
        let rho: Vec<_> = (0..elasticity.cells())
            .map(|i| 0.42 + 0.03 * i as f64)
            .collect();
        let direction: Vec<_> = (0..rho.len())
            .map(|i| 0.03 * ((i % 3) as f64 - 1.0))
            .collect();
        let mut callback = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let goal = if robust {
            SimpThermalComplianceGoal::new_robust(
                &three_field,
                &elasticity,
                &rho,
                &loads,
                options(),
                &mut control,
            )
        } else {
            SimpThermalComplianceGoal::new(
                &ordinary,
                &elasticity,
                &rho,
                &loads,
                options(),
                &mut control,
            )
        }
        .unwrap();
        assert_eq!(goal.load_case_count(), 2);
        match goal.physical().thermal_field(0).unwrap() {
            TetThermalStrainField::PerElement(retained) => assert_eq!(retained, thermal),
            _ => panic!("per-element thermal sources must remain attached"),
        }
        let hv = goal
            .hessian_vector_product(&direction, &mut control)
            .unwrap();
        let case_hv = goal
            .load_case_hessian_vector_product(1, &direction, &mut control)
            .unwrap();
        let first = reference
            .try_thermal_compliance_and_gradient(&elasticity, &rho, &loads, options(), &mut control)
            .unwrap();
        close(goal.value(), first.goal.value(), 2e-11);
        vector_close(goal.gradient(), &first.gradient, 2e-10);
        let eps = 1e-4;
        let mut differences = Vec::new();
        for sign in [1.0, -1.0] {
            let perturbed: Vec<_> = rho
                .iter()
                .zip(&direction)
                .map(|(r, v)| r + sign * eps * v)
                .collect();
            let weighted = reference
                .try_thermal_compliance_and_gradient(
                    &elasticity,
                    &perturbed,
                    &loads,
                    options(),
                    &mut control,
                )
                .unwrap();
            let selected = [ThermalLoadCase {
                weight: 1.0,
                ..loads[1]
            }];
            let single = reference
                .try_thermal_compliance_and_gradient(
                    &elasticity,
                    &perturbed,
                    &selected,
                    options(),
                    &mut control,
                )
                .unwrap();
            differences.push((weighted.gradient, single.gradient));
        }
        let fd: Vec<_> = differences[0]
            .0
            .iter()
            .zip(&differences[1].0)
            .map(|(a, b)| (a - b) / (2.0 * eps))
            .collect();
        let case_fd: Vec<_> = differences[0]
            .1
            .iter()
            .zip(&differences[1].1)
            .map(|(a, b)| (a - b) / (2.0 * eps))
            .collect();
        vector_close(&hv, &fd, 3e-5);
        vector_close(&case_hv, &case_fd, 3e-5);
        let volumes: Vec<_> = fs_feec::element_geometry(&mesh, &positions)
            .vol_signed
            .iter()
            .map(|v| v.abs())
            .collect();
        let inventory = DesignPipeline {
            filter: DensityFilter::new(&mesh, &positions, 0.15),
            params: SimpParams {
                eta: if robust { 0.42 } else { 0.5 },
                ..p
            },
        };
        let (_, expected, _) = inventory.try_forward(&rho, &mut control).unwrap();
        vector_close(goal.volume_projected_densities(), &expected, 2e-12);
        let (volume, _) = goal.volume_value_gradient(&volumes, &mut control).unwrap();
        close(
            volume,
            expected
                .iter()
                .zip(&volumes)
                .map(|(r, v)| r * v)
                .sum::<f64>()
                / volumes.iter().sum::<f64>(),
            2e-12,
        );
        let (mass, _) = goal.mass_value_gradient(&mut control).unwrap();
        close(
            mass,
            expected
                .iter()
                .zip(elasticity.reference_cell_masses_kg().unwrap())
                .map(|(r, m)| r * m)
                .sum(),
            2e-12,
        );
        let before = control.work();
        assert_eq!(
            goal.hessian_vector_product(&vec![0.0; rho.len()], &mut control)
                .unwrap(),
            vec![0.0; rho.len()]
        );
        assert_eq!(control.work(), before);
    }
}

#[test]
fn zero_thermal_field_reduces_to_original_dead_load_calculus_and_solve_counts() {
    let (mesh, positions, elasticity, thermal, a, b) = cube();
    let zero: Vec<_> = thermal
        .iter()
        .map(|s| {
            TetThermalStrainState::try_new(
                300.0,
                300.0,
                [0.0; 6],
                s.elastic_material_state_identity(),
                s.thermal_law_identity(),
                0.2,
            )
            .unwrap()
        })
        .collect();
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&mesh, &positions, 0.15),
        params: SimpParams {
            e_min: 0.03,
            penal: 2.7,
            beta: 2.3,
            eta: 0.5,
        },
    };
    let dead = [
        LoadCase {
            force: &a,
            weight: 0.7,
        },
        LoadCase {
            force: &b,
            weight: 1.3,
        },
    ];
    let loads: Vec<_> = dead
        .iter()
        .map(|l| ThermalLoadCase {
            force: l.force,
            weight: l.weight,
            strains: TetThermalStrainField::PerElement(&zero),
        })
        .collect();
    let rho: Vec<_> = (0..elasticity.cells())
        .map(|i| 0.44 + 0.02 * i as f64)
        .collect();
    let direction = vec![0.2; rho.len()];
    let mut callback_a = |_| ControlFlow::Continue(());
    let mut callback_b = |_| ControlFlow::Continue(());
    let mut ca = SolveControl::new(SolveBudget::default(), &mut callback_a);
    let mut cb = SolveControl::new(SolveBudget::default(), &mut callback_b);
    let original =
        SimpComplianceGoal::new(&pipeline, &elasticity, &rho, &dead, options(), &mut ca).unwrap();
    let heated =
        SimpThermalComplianceGoal::new(&pipeline, &elasticity, &rho, &loads, options(), &mut cb)
            .unwrap();
    assert_eq!(ca.work(), cb.work());
    assert_eq!(original.value(), heated.value());
    vector_close(original.gradient(), heated.gradient(), 1e-12);
    let old = original
        .hessian_vector_product(&direction, &mut ca)
        .unwrap();
    let new = heated.hessian_vector_product(&direction, &mut cb).unwrap();
    vector_close(&new, &old, 2e-10);
    assert_eq!(ca.work().linear_solves, cb.work().linear_solves);
}
