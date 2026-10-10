//! G1 physical endpoint balance and G0 immutable-history radiation refusal.

use super::*;
use fs_conduction::{
    ConductivityTable, LinearConfig,
    transient::{
        VolumetricHeatCapacity,
        backward_euler::{
            BackwardEuler, NonlinearStepConfig, RadiationStepSolution, StepConfig,
        },
    },
};

const CAPACITY_PER_NODE: f64 = 15.0 / 24.0;
const DT: f64 = 0.25;

fn tet() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
    )
    .unwrap()
}

fn driven_boundary(mesh: &ConductionMesh, fixed: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "driven",
            |face| (face.centroid.iter().sum::<f64>() - 1.0).abs() < 1e-12,
            ThermalBc::dirichlet(fixed).unwrap(),
        )
        .unwrap()
        .remainder("cooler", ThermalBc::robin(5.0, 300.0).unwrap())
        .unwrap()
        .finish()
        .unwrap()
}

fn step_controls() -> StepConfig {
    StepConfig {
        linear: LinearConfig {
            tolerance: 1e-12,
            ..LinearConfig::default()
        },
        energy_tolerance_j: 1e-7,
    }
}

/// Independent one-free-node residual. This unit right tet has V=1/6,
/// lumped C_00=rho_cp/24, K_00=k/2, and three coordinate faces of area 1/2.
/// The Robin row is h*(t+fixed-2*reference)/4, rather than a uniform flux.
/// Radiation uses fourth powers here, independently of the production secant.
fn residual(t: f64, old: f64, fixed: f64, ambient: f64, nonlinear: bool) -> f64 {
    let surface_mean = (t + 2.0 * fixed) / 3.0;
    let element_mean = (t + 3.0 * fixed) / 4.0;
    let k = if nonlinear { 0.1 * element_mean - 10.0 } else { 20.0 };
    let h_rad = if surface_mean == ambient {
        0.8 * STEFAN_BOLTZMANN_W_M2_K4 * 4.0 * ambient.powi(3)
    } else {
        0.8 * STEFAN_BOLTZMANN_W_M2_K4
            * (surface_mean.powi(4) - ambient.powi(4)) / (surface_mean - ambient)
    };
    CAPACITY_PER_NODE * (t - old)
        + DT * (
            (k / 2.0) * (t - fixed)
                + (5.0 / 4.0) * (t + fixed - 600.0)
                + (h_rad / 4.0) * (t + fixed - 2.0 * ambient)
        )
}

fn oracle(old: f64, fixed: f64, ambient: f64, nonlinear: bool) -> f64 {
    let mut low = 250.0;
    let mut high = 500.0;
    assert!(residual(low, old, fixed, ambient, nonlinear) < 0.0);
    assert!(residual(high, old, fixed, ambient, nonlinear) > 0.0);
    for _ in 0..90 {
        let mid = (low + high) * 0.5;
        if residual(mid, old, fixed, ambient, nonlinear) > 0.0 {
            high = mid;
        } else {
            low = mid;
        }
    }
    (low + high) * 0.5
}

fn check_endpoint(
    result: &RadiationStepSolution,
    old: &[f64; 4],
    fixed: f64,
    ambient: f64,
    nonlinear: bool,
) {
    let t = oracle(old[0], fixed, ambient, nonlinear);
    assert!((result.conduction.temperature[0] - t).abs() < 2e-8);
    assert_eq!(&result.conduction.temperature[1..], &[fixed; 3]);
    let actual = result.conduction.temperature[0];
    assert!(residual(actual, old[0], fixed, ambient, nonlinear).abs() < 1e-8);
    assert!(result.physical_residual_norm_j <= result.physical_residual_tolerance_j);
    assert!(result.physical_energy_residual_j.abs() <= step_controls().energy_tolerance_j);
    let mean = (actual + 2.0 * fixed) / 3.0;
    let radiation = 1.5 * 0.8 * STEFAN_BOLTZMANN_W_M2_K4
        * (mean.powi(4) - ambient.powi(4));
    let convection = 1.5 * 5.0 * (mean - 300.0);
    let stored = CAPACITY_PER_NODE
        * result.conduction.temperature.iter().zip(old).map(|(t, old)| t - old).sum::<f64>();
    assert!((result.radiation.nonlinear_radiation_out_w - radiation).abs() < 1e-8);
    assert!((result.convective_out_w - convection).abs() < 1e-8);
    assert!((result.conduction.stored_energy_change_j - stored).abs() < 1e-9);
    // This includes storage at ALL prescribed vertices during the first jump.
    assert!((result.physical_dirichlet_in_w - stored / DT - convection - radiation).abs() < 1e-6);
    assert_eq!(result.convective_robin_fluxes[0].mean_htc_w_per_m2_k, 5.0);
    assert!(result.radiation.iterations > 1);
    assert!(result.radiation.krylov_iterations >= result.conduction.krylov_iterations);
    assert!(result.radiation.krylov_iterations
        <= result.radiation.iterations * step_controls().linear.max_iterations);
}

#[test]
fn radiative_prescribed_step_matches_scalar_balance_and_keeps_history() {
    with_cx(|cx| {
        let mesh = tet();
        let boundary = driven_boundary(&mesh, 330.0);
        let material = ConductivityModel::isotropic_declared(20.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let old = [300.0; 4];
        let engine = BackwardEuler::uniform(
            cx, &mesh, VolumetricHeatCapacity::declared(15.0).unwrap(),
        ).unwrap();
        let problem = ConductionProblem {
            mesh: &mesh, boundary: &boundary, material: &material,
            source: &source, element_materials: None,
        };
        let solve = |history: &[f64], ambient| engine.advance_prescribed_with_ambient_radiation(
            cx, problem, None, history, DT, step_controls(), None,
            &[patch(0.8, ambient)], AmbientRadiationConfig::default(),
        ).unwrap();
        let cold = solve(&old, 300.0);
        let hot = solve(&old, 450.0);
        check_endpoint(&cold, &old, 330.0, 300.0, false);
        check_endpoint(&hot, &old, 330.0, 450.0, false);
        assert!(cold.radiation.nonlinear_radiation_out_w > 0.0);
        assert!(hot.radiation.nonlinear_radiation_out_w < 0.0);
        assert!(hot.conduction.temperature[0] > cold.conduction.temperature[0]);
        let bare = engine.advance_prescribed(
            cx, problem, None, &old, DT, step_controls(),
        ).unwrap();
        assert!(cold.conduction.temperature[0] < bare.temperature[0]);
        let repeated = solve(&old, 300.0);
        assert_eq!(cold.conduction.temperature, repeated.conduction.temperature);
        assert_eq!(cold.radiation, repeated.radiation);
        assert_eq!(old, [300.0; 4]);
        let accepted: [f64; 4] = cold.conduction.temperature.clone().try_into().unwrap();
        let next = solve(&accepted, 300.0);
        check_endpoint(&next, &accepted, 330.0, 300.0, false);
        let cumulative = cold.conduction.stored_energy_change_j + next.conduction.stored_energy_change_j;
        let expected = CAPACITY_PER_NODE
            * next.conduction.temperature.iter().map(|t| t - 300.0).sum::<f64>();
        assert!((cumulative - expected).abs() < 1e-9);
    });
}

#[test]
fn radiation_keeps_endpoint_nonlinearity_assignment_and_constant_law_twin() {
    with_cx(|cx| {
        let mesh = tet();
        let boundary = driven_boundary(&mesh, 330.0);
        let fallback = ConductivityModel::isotropic_declared(999.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let curve = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(250.0, 15.0), (500.0, 40.0)]).unwrap(),
        );
        let assigned = ElementMaterials::new(
            MaterialTable::new([(MaterialId(7), curve)]).unwrap(), vec![MaterialId(7)],
        ).unwrap();
        let engine = BackwardEuler::per_element(
            cx, &mesh, &[VolumetricHeatCapacity::declared(15.0).unwrap()],
        ).unwrap();
        let problem = ConductionProblem {
            mesh: &mesh, boundary: &boundary, material: &fallback,
            source: &source, element_materials: Some(&assigned),
        };
        let old = [300.0; 4];
        let solve = |problem: ConductionProblem<'_>, nonlinear: Option<NonlinearStepConfig>| engine.advance_prescribed_with_ambient_radiation(
            cx, problem, None, &old, DT, step_controls(), nonlinear,
            &[patch(0.8, 300.0)], AmbientRadiationConfig::default(),
        );
        assert!(solve(problem, None).is_err());
        let nonlinear = solve(problem, Some(NonlinearStepConfig::default())).unwrap();
        check_endpoint(&nonlinear, &old, 330.0, 300.0, true);
        let accepted = nonlinear.nonlinear.as_ref().unwrap();
        assert_eq!(accepted.step.temperature, nonlinear.conduction.temperature);
        assert!(accepted.residual_j <= accepted.threshold_j);
        assert!(nonlinear.radiation.solid_iterations > accepted.nonlinear_iterations);
        assert!(nonlinear.nonlinear_backtracks >= accepted.backtracks);
        let constant = ConductivityModel::isotropic_declared(20.0).unwrap();
        let constant_problem = ConductionProblem {
            material: &constant, element_materials: None, ..problem
        };
        let linear = solve(constant_problem, None).unwrap();
        let twin = solve(constant_problem, Some(NonlinearStepConfig::default())).unwrap();
        for (a, b) in linear.conduction.temperature.iter().zip(&twin.conduction.temperature) {
            assert!((a - b).abs() < 2e-8);
        }
        assert_eq!(old, [300.0; 4]);
    });
}

#[test]
fn loose_outer_gates_cannot_accept_an_unbalanced_frozen_endpoint() {
    with_cx(|cx| {
        let mesh = tet();
        let boundary = driven_boundary(&mesh, 330.0);
        let material = ConductivityModel::isotropic_declared(20.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let old = [300.0; 4];
        let engine = BackwardEuler::uniform(
            cx, &mesh, VolumetricHeatCapacity::declared(15.0).unwrap(),
        ).unwrap();
        let problem = ConductionProblem {
            mesh: &mesh, boundary: &boundary, material: &material,
            source: &source, element_materials: None,
        };
        let loose = AmbientRadiationConfig {
            max_iterations: 1,
            temperature_tolerance_k: 1e6,
            balance_tolerance_w: 1e9,
            balance_relative_tolerance: 0.0,
            relaxation: 1.0,
        };
        let result = engine.advance_prescribed_with_ambient_radiation(
            cx, problem, None, &old, DT, step_controls(), None,
            &[patch(0.8, 300.0)], loose,
        );
        assert!(matches!(result,
            Err(ConductionError::AmbientRadiationNotConverged { iterations: 1, .. })
        ));
        assert_eq!(old, [300.0; 4]);
        let retry = engine.advance_prescribed_with_ambient_radiation(
            cx, problem, None, &old, DT, step_controls(), None,
            &[patch(0.8, 300.0)], AmbientRadiationConfig::default(),
        ).unwrap();
        check_endpoint(&retry, &old, 330.0, 300.0, false);
    });
}

#[test]
fn radiative_transient_refuses_overflow_bad_binding_validity_and_cancellation() {
    support::with_gate(|gate, cx| {
        let mesh = tet();
        let boundary = driven_boundary(&mesh, 330.0);
        let material = ConductivityModel::isotropic_declared(20.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let old = [300.0; 4];
        let engine = BackwardEuler::uniform(
            cx, &mesh, VolumetricHeatCapacity::declared(15.0).unwrap(),
        ).unwrap();
        let problem = ConductionProblem {
            mesh: &mesh, boundary: &boundary, material: &material,
            source: &source, element_materials: None,
        };
        for patches in [
            Vec::new(),
            vec![patch(0.8, 300.0), patch(0.8, 300.0)],
            vec![AmbientRadiationPatch::new(
                "driven", patch(0.8, 300.0).emissivity().clone(), 300.0,
            ).unwrap()],
        ] {
            assert!(engine.advance_prescribed_with_ambient_radiation(
                cx, problem, None, &old, DT, step_controls(), None,
                &patches, AmbientRadiationConfig::default(),
            ).is_err());
        }
        let overflowing = AmbientRadiationConfig {
            max_iterations: usize::MAX,
            ..AmbientRadiationConfig::default()
        };
        assert!(engine.advance_prescribed_with_ambient_radiation(
            cx, problem, None, &old, DT, step_controls(), None,
            &[patch(0.8, 300.0)], overflowing,
        ).is_err());
        let hot_boundary = driven_boundary(&mesh, 600.0);
        let hot = ConductionProblem { boundary: &hot_boundary, ..problem };
        assert!(engine.advance_prescribed_with_ambient_radiation(
            cx, hot, None, &[600.0; 4], DT, step_controls(), None,
            &[patch(0.8, 300.0)], AmbientRadiationConfig::default(),
        ).is_err());
        gate.request();
        assert!(matches!(engine.advance_prescribed_with_ambient_radiation(
            cx, problem, None, &old, DT, step_controls(), None,
            &[patch(0.8, 300.0)], AmbientRadiationConfig::default(),
        ), Err(ConductionError::Cancelled { .. })));
        assert_eq!(old, [300.0; 4]);
    });
}
