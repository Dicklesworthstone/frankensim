//! Thermal design must differentiate its equivalent load as well as stiffness.

#[path = "thermal_compliance/simp.rs"]
mod simp;

#[path = "thermal_compliance/stress.rs"]
mod stress;

#[path = "thermal_compliance/modal.rs"]
mod modal;

use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_rep_mesh::TetComplex;
use fs_solid::linear3::{
    TetElasticMaterial, TetMaterialField, TetThermalStrainField, TetThermalStrainState,
};
use fs_topopt::elasticity::second_order::thermal::{ThermalComplianceGoal, ThermalLoadCase};
use fs_topopt::elasticity::second_order::{ComplianceError, ComplianceSolveOptions};
use fs_topopt::{
    DensityElasticity, DensityFilter, DesignPipeline, EvaluationStop, SimpParams, SolveBudget,
    SolveControl,
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance * expected.abs().max(1e-6),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn state(material: &TetElasticMaterial, strain: [f64; 6], bound: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        320.0,
        300.0,
        strain,
        material.material_state_identity,
        ContentHash([21; 32]),
        bound,
    )
    .unwrap()
}

fn options() -> ComplianceSolveOptions {
    ComplianceSolveOptions {
        tolerance: 1e-11,
        max_iterations: 5000,
    }
}

fn solve<'a>(
    elasticity: &'a DensityElasticity,
    weights: &[f64],
    loads: &[ThermalLoadCase<'_>],
) -> ThermalComplianceGoal<'a> {
    let mut callback = |_| ControlFlow::Continue(());
    ThermalComplianceGoal::new(
        elasticity,
        weights,
        loads,
        options(),
        &mut SolveControl::new(SolveBudget::default(), &mut callback),
    )
    .unwrap()
}

fn one_tet(clamped: bool) -> (DensityElasticity, TetElasticMaterial) {
    let complex = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let material = TetElasticMaterial::try_new(6.0, 120.0, 0.0, [7; 32]).unwrap();
    let elasticity = DensityElasticity::with_material_field(
        &complex,
        &positions,
        TetMaterialField::Uniform(&material),
        &|point| clamped || point[2] == 0.0,
    )
    .unwrap();
    (elasticity, material)
}

#[test]
fn heated_axial_tet_has_positive_thermal_sensitivity_and_exact_tangent_solve() {
    let (elasticity, material) = one_tet(false);
    let alpha = 0.01;
    let thermal = state(&material, [0.0, 0.0, alpha, 0.0, 0.0, 0.0], 0.1);
    let mut force = vec![0.0; elasticity.n()];
    let stiffness = 120.0 / 6.0;
    let m = 0.4;
    for applied in [0.0, 0.02] {
        force[11] = applied;
        let loads = [ThermalLoadCase {
            force: &force,
            weight: 1.0,
            strains: TetThermalStrainField::Uniform(&thermal),
        }];
        let goal = solve(&elasticity, &[m], &loads);
        let displacement = alpha + applied / (m * stiffness);
        close(goal.displacements()[0][11], displacement, 2e-12);
        close(
            goal.value(),
            (applied + m * stiffness * alpha) * displacement,
            2e-12,
        );
        close(
            goal.gradient()[0],
            stiffness * alpha * alpha - applied * applied / (m * m * stiffness),
            2e-12,
        );
        assert!(
            goal.gradient()[0] > 0.0,
            "dead-load-only sensitivity has the wrong sign"
        );
        let mut callback = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let direction = 0.3;
        let hv = goal
            .hessian_vector_product(&[direction], &mut control)
            .unwrap();
        close(
            hv[0],
            2.0 * direction * applied * applied / (m * m * m * stiffness),
            2e-8,
        );
        close(
            goal.load_case_gradient(0, &mut control).unwrap()[0],
            goal.gradient()[0],
            2e-12,
        );
        let before = control.work();
        assert_eq!(
            goal.hessian_vector_product(&[0.0], &mut control).unwrap(),
            [0.0]
        );
        assert_eq!(control.work(), before);
        let response = goal.response(0).unwrap();
        close(
            response.elastic_energy_j,
            0.5 * applied * applied / (m * stiffness),
            1e-9,
        );
        close(
            response.support_reactions_n.iter().skip(2).step_by(3).sum(),
            -applied,
            1e-8,
        );
        close(
            response.elements[0].effective_stress_mandel_pa[2],
            6.0 * applied,
            1e-8,
        );
        assert!(goal.primal_relative_residuals()[0] < options().tolerance);
        match goal.thermal_field(0).unwrap() {
            TetThermalStrainField::Uniform(retained) => assert_eq!(retained, &thermal),
            _ => panic!("uniform thermal field must remain uniform"),
        }
    }
}

#[test]
fn coupled_tensor_free_expansion_preserves_mandel_shear_work() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mut c = [[0.0; 6]; 6];
    for (i, value) in [20.0, 30.0, 50.0, 16.0, 24.0, 32.0].into_iter().enumerate() {
        c[i][i] = value;
    }
    for (i, j, value) in [(0, 4, 2.0), (2, 5, -3.0), (4, 5, 4.0)] {
        c[i][j] = value;
        c[j][i] = value;
    }
    let material = TetElasticMaterial::try_new_mandel(4.0, c, [13; 32]).unwrap();
    let elasticity = DensityElasticity::with_material_field(
        &mesh,
        &positions,
        TetMaterialField::Uniform(&material),
        &|p| p[2] == 0.0,
    )
    .unwrap();
    let free = [0.0, 0.0, 0.003, 0.0, -0.002, 0.004];
    let thermal = state(&material, free, 0.1);
    let force = vec![0.0; elasticity.n()];
    let goal = solve(
        &elasticity,
        &[0.7],
        &[ThermalLoadCase {
            force: &force,
            weight: 1.0,
            strains: TetThermalStrainField::Uniform(&thermal),
        }],
    );
    for (actual, expected) in goal.displacements()[0][9..].iter().zip([
        core::f64::consts::SQRT_2 * free[5],
        core::f64::consts::SQRT_2 * free[4],
        free[2],
    ]) {
        close(*actual, expected, 1e-10);
    }
    let response = goal.response(0).unwrap();
    close(response.elastic_energy_j, 0.0, 1e-9);
    for stress in response.elements[0].effective_stress_mandel_pa {
        close(stress, 0.0, 1e-8);
    }
    for reaction in &response.support_reactions_n {
        close(*reaction, 0.0, 1e-8);
    }
    let free_energy: f64 = (0..6)
        .map(|i| free[i] * c[i].iter().zip(free).map(|(c, e)| c * e).sum::<f64>())
        .sum::<f64>()
        / 6.0;
    close(goal.value(), 0.7 * free_energy, 1e-10);
    close(goal.gradient()[0], free_energy, 1e-10);
}

#[test]
fn full_restraint_has_zero_compliance_but_nonzero_thermal_energy_stress_and_reactions() {
    let (elasticity, material) = one_tet(true);
    let thermal = state(&material, [0.01, -0.005, 0.002, 0.003, -0.004, 0.001], 0.1);
    let mut force = vec![0.0; elasticity.n()];
    force[0] = 0.7;
    let loads = [ThermalLoadCase {
        force: &force,
        weight: 2.0,
        strains: TetThermalStrainField::Uniform(&thermal),
    }];
    let m = 0.6;
    let goal = solve(&elasticity, &[m], &loads);
    assert_eq!(goal.value(), 0.0);
    assert_eq!(goal.gradient(), [0.0]);
    assert!(goal.displacements()[0].iter().all(|u| *u == 0.0));
    let response = goal.response(0).unwrap();
    let strain = thermal.free_strain_mandel();
    let stress: [f64; 6] = core::array::from_fn(|i| {
        -m * material.stiffness_mandel_pa()[i]
            .iter()
            .zip(strain)
            .map(|(c, e)| c * e)
            .sum::<f64>()
    });
    for (actual, expected) in response.elements[0]
        .effective_stress_mandel_pa
        .iter()
        .zip(stress)
    {
        close(*actual, expected, 1e-12);
    }
    close(
        response.elastic_energy_j,
        -strain.iter().zip(stress).map(|(e, s)| e * s).sum::<f64>() / 12.0,
        1e-12,
    );
    for component in 0..3 {
        close(
            response
                .support_reactions_n
                .iter()
                .skip(component)
                .step_by(3)
                .sum(),
            if component == 0 { -0.7 } else { 0.0 },
            1e-5,
        );
    }
}

fn cube() -> (
    TetComplex,
    Vec<[f64; 3]>,
    DensityElasticity,
    Vec<TetThermalStrainState>,
    Vec<f64>,
    Vec<f64>,
) {
    let (mesh, positions) = fs_feec::kuhn_cube(1);
    let mut materials = Vec::new();
    let mut thermal = Vec::new();
    for cell in 0..mesh.tets.len() {
        let mut c = [[0.0; 6]; 6];
        for (i, row) in c.iter_mut().enumerate() {
            row[i] = 30.0 + 3.0 * i as f64 + cell as f64;
        }
        for (i, j, value) in [(0, 3, 3.0), (1, 5, -2.0), (2, 4, 1.5), (3, 4, 2.0)] {
            c[i][j] = value;
            c[j][i] = value;
        }
        let material =
            TetElasticMaterial::try_new_mandel(2.0 + cell as f64, c, [cell as u8 + 1; 32]).unwrap();
        thermal.push(state(
            &material,
            [
                0.002,
                -0.001,
                0.0015,
                0.0003 * (cell + 1) as f64,
                -0.0004,
                0.0005,
            ],
            0.2,
        ));
        materials.push(material);
    }
    let elasticity = DensityElasticity::with_material_field(
        &mesh,
        &positions,
        TetMaterialField::PerElement(&materials),
        &|p| p[0] == 0.0,
    )
    .unwrap();
    let mut a = vec![0.0; elasticity.n()];
    let mut b = a.clone();
    for (node, p) in positions.iter().enumerate() {
        if p[0] == 1.0 {
            a[3 * node + 2] = -0.006;
            b[3 * node + 1] = 0.004 * (1.0 + p[2]);
        }
    }
    (mesh, positions, elasticity, thermal, a, b)
}

#[test]
fn heterogeneous_tensor_thermal_gradients_and_hessian_match_design_perturbations() {
    let (_, _, elasticity, thermal, a, b) = cube();
    let cooling: Vec<_> = thermal
        .iter()
        .map(|s| {
            TetThermalStrainState::try_new(
                290.0,
                300.0,
                s.free_strain_mandel().map(|e| -0.7 * e),
                s.elastic_material_state_identity(),
                s.thermal_law_identity(),
                0.2,
            )
            .unwrap()
        })
        .collect();
    let loads = [
        ThermalLoadCase {
            force: &a,
            weight: 0.7,
            strains: TetThermalStrainField::PerElement(&thermal),
        },
        ThermalLoadCase {
            force: &b,
            weight: 1.3,
            strains: TetThermalStrainField::PerElement(&cooling),
        },
    ];
    let weights: Vec<_> = (0..elasticity.cells())
        .map(|i| 0.5 + 0.05 * i as f64)
        .collect();
    let direction: Vec<_> = (0..elasticity.cells())
        .map(|i| 0.1 * ((i % 3) as f64 - 1.0))
        .collect();
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let goal =
        ThermalComplianceGoal::new(&elasticity, &weights, &loads, options(), &mut control).unwrap();
    assert_eq!(control.work().linear_solves, 2);
    let hv = goal
        .hessian_vector_product(&direction, &mut control)
        .unwrap();
    assert_eq!(control.work().linear_solves, 4);
    let step = 1e-4;
    let plus: Vec<_> = weights
        .iter()
        .zip(&direction)
        .map(|(m, v)| m + step * v)
        .collect();
    let minus: Vec<_> = weights
        .iter()
        .zip(&direction)
        .map(|(m, v)| m - step * v)
        .collect();
    let gp = solve(&elasticity, &plus, &loads);
    let gm = solve(&elasticity, &minus, &loads);
    for i in 0..weights.len() {
        close(
            hv[i],
            (gp.gradient()[i] - gm.gradient()[i]) / (2.0 * step),
            3e-5,
        );
        let mut plus = weights.clone();
        let mut minus = weights.clone();
        plus[i] += step;
        minus[i] -= step;
        close(
            goal.gradient()[i],
            (solve(&elasticity, &plus, &loads).value()
                - solve(&elasticity, &minus, &loads).value())
                / (2.0 * step),
            2e-6,
        );
    }
    let first = goal.load_case_gradient(0, &mut control).unwrap();
    let second = goal.load_case_gradient(1, &mut control).unwrap();
    for i in 0..weights.len() {
        close(goal.gradient()[i], 0.7 * first[i] + 1.3 * second[i], 2e-12);
    }
    let second_direction: Vec<_> = direction
        .iter()
        .enumerate()
        .map(|(i, v)| v + 0.03 * (i + 1) as f64)
        .collect();
    let hw = goal
        .hessian_vector_product(&second_direction, &mut control)
        .unwrap();
    close(
        direction.iter().zip(hw).map(|(v, h)| v * h).sum(),
        second_direction.iter().zip(hv).map(|(w, h)| w * h).sum(),
        1e-8,
    );
}

#[test]
fn complete_filter_projection_simp_pullback_includes_thermal_load_variation() {
    let (mesh, positions, elasticity, thermal, force, _) = cube();
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&mesh, &positions, 0.12),
        params: SimpParams {
            e_min: 0.01,
            penal: 3.0,
            beta: 2.0,
            eta: 0.5,
        },
    };
    let loads = [ThermalLoadCase {
        force: &force,
        weight: 1.0,
        strains: TetThermalStrainField::PerElement(&thermal),
    }];
    let rho: Vec<_> = (0..elasticity.cells())
        .map(|i| 0.6 + 0.02 * i as f64)
        .collect();
    let evaluate = |rho: &[f64]| {
        let mut callback = |_| ControlFlow::Continue(());
        pipeline
            .try_thermal_compliance_and_gradient(
                &elasticity,
                rho,
                &loads,
                options(),
                &mut SolveControl::new(SolveBudget::default(), &mut callback),
            )
            .unwrap()
    };
    let result = evaluate(&rho);
    let step = 1e-5;
    for i in 0..rho.len() {
        let mut plus = rho.clone();
        let mut minus = rho.clone();
        plus[i] += step;
        minus[i] -= step;
        close(
            result.gradient[i],
            (evaluate(&plus).goal.value() - evaluate(&minus).goal.value()) / (2.0 * step),
            5e-6,
        );
    }
    assert_eq!(elasticity.moduli, vec![1.0; elasticity.cells()]);
}

#[test]
fn stale_sources_domain_escape_and_stopped_work_never_publish_a_thermal_goal() {
    let (elasticity, material) = one_tet(false);
    let thermal = state(&material, [0.0, 0.0, 0.01, 0.0, 0.0, 0.0], 0.02);
    let mut force = vec![0.0; elasticity.n()];
    let other = TetElasticMaterial::try_new(6.0, 120.0, 0.0, [8; 32]).unwrap();
    let stale = state(&other, thermal.free_strain_mandel(), 0.02);
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    for strains in [
        TetThermalStrainField::Uniform(&stale),
        TetThermalStrainField::PerElement(&[]),
    ] {
        assert!(matches!(
            ThermalComplianceGoal::new(
                &elasticity,
                &[1.0],
                &[ThermalLoadCase {
                    force: &force,
                    weight: 1.0,
                    strains
                }],
                options(),
                &mut control
            ),
            Err(ComplianceError::InvalidInput(_))
        ));
    }
    assert_eq!(control.work().linear_solves, 0);
    force[11] = 1.0;
    let loads = [ThermalLoadCase {
        force: &force,
        weight: 1.0,
        strains: TetThermalStrainField::Uniform(&thermal),
    }];
    assert!(matches!(
        ThermalComplianceGoal::new(&elasticity, &[1.0], &loads, options(), &mut control),
        Err(ComplianceError::InvalidInput(_))
    ));
    let mut callback = |_| ControlFlow::Continue(());
    let mut zero_budget = SolveControl::new(
        SolveBudget {
            per_solve_iterations: 0,
            total_iterations: 0,
        },
        &mut callback,
    );
    assert!(matches!(
        ThermalComplianceGoal::new(&elasticity, &[1.0], &loads, options(), &mut zero_budget),
        Err(ComplianceError::Stopped(EvaluationStop::TotalBudget { .. }))
    ));
    let mut cancel = |_| ControlFlow::Break(());
    let mut stopped = SolveControl::new(SolveBudget::default(), &mut cancel);
    assert!(matches!(
        ThermalComplianceGoal::new(&elasticity, &[1.0], &loads, options(), &mut stopped),
        Err(ComplianceError::Stopped(EvaluationStop::Cancelled))
    ));
    assert_eq!(stopped.work().linear_solves, 0);
}
