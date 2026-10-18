//! Heterogeneous thermal inertia through the original history and endpoint owners.
mod support;

use fs_conduction::bc::history::{ThermalBoundaryHistory, ThermalRegionHistory};
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryCapacity,
    ThermalHistoryMethod, ThermalHistoryProblem, ThermalHistoryRun, ThermalHistoryStep};
use fs_conduction::transient::{VolumetricHeatCapacity, backward_euler::{BackwardEuler,
    NonlinearStepConfig, StepConfig}};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    ConductivityTable, ElementMaterials, LinearConfig, MaterialId, MaterialTable,
    ScalarField, ThermalBc, ThermalBoundaryBuilder, TEMPERATURE_DIMS};
use fs_exec::Cx;
use fs_qty::QtyAny;
use fs_rep_mesh::TetComplex;
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Physics};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::TimeSignal;
use support::{with_cancelled_cx, with_cx};

fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
            [2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.0, 0.0, 1.0]]).unwrap()
}
fn cap(value: f64) -> VolumetricHeatCapacity { VolumetricHeatCapacity::declared(value).unwrap() }
fn controls() -> StepConfig {
    StepConfig { linear: LinearConfig { tolerance: 1e-11, max_iterations: 4096, restart: 32 },
        energy_tolerance_j: 1e-7 }
}
fn materials(nonlinear: bool) -> ElementMaterials {
    let (a, b) = if nonlinear {
        (ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(250.0, 10.0), (500.0, 50.0)]).unwrap()),
         ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(250.0, 20.0), (500.0, 25.0)]).unwrap()))
    } else {
        (ConductivityModel::isotropic_declared(20.0).unwrap(),
         ConductivityModel::isotropic_declared(40.0).unwrap())
    };
    ElementMaterials::new(MaterialTable::new([(MaterialId(7), a), (MaterialId(2), b)]).unwrap(),
        vec![MaterialId(7), MaterialId(2)]).unwrap()
}
fn history<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh).adiabatic_remainder(),
        &Environment::earth_lab(), Vec::new(), 0.0).unwrap()
}
fn problem<'h, 'm>(history: &'h ThermalBoundaryHistory<'m>, fallback: &'h ConductivityModel,
    assigned: Option<&'h ElementMaterials>) -> ThermalHistoryProblem<'h, 'm>
{
    ThermalHistoryProblem { boundary: history, material: fallback, element_materials: assigned,
        interfaces: None, source: HeatSourceHistory::Fixed(ScalarField::Uniform(30.0)) }
}

#[test]
fn g1_distinct_bodies_heat_at_their_own_capacity_not_an_averaged_rate() {
    let mesh = mesh();
    let fallback = ConductivityModel::isotropic_declared(999.0).unwrap();
    let assigned = materials(false);
    let rows = [(MaterialId(2), cap(150.0)), (MaterialId(7), cap(15.0))];
    with_cx(|cx| {
        let history = history(cx, &mesh);
        for method in [ThermalHistoryMethod::Linear,
            ThermalHistoryMethod::Nonlinear(NonlinearStepConfig::default())] {
            let mut run = ThermalHistoryRun::new_with_capacity(cx,
                problem(&history, &fallback, Some(&assigned)), ThermalHistoryCapacity::ByMaterial(&rows),
                &[300.0; 8], 3, method).unwrap();
            let mut previous = 0.0;
            for time in [0.25, 0.5, 1.0] {
                let solved = run.advance_to(cx, time, controls()).unwrap();
                for (element, capacity) in [15.0, 150.0].into_iter().enumerate() {
                    let expected = 300.0 + 30.0 * time / capacity;
                    assert!(solved.temperature[4*element..4*element+4].iter()
                        .all(|&t| (t-expected).abs() < 1e-8));
                }
                // Two unit right tets, each V=1/6: generation is 10 W total.
                assert!((solved.stored_energy_change_j - 10.0*(time-previous)).abs() < 1e-7);
                assert!(solved.energy_residual_j.abs() <= controls().energy_tolerance_j);
                previous = time;
            }
            assert!(run.temperature()[0] - run.temperature()[4] > 1.0);
        }
    });
}

#[test]
fn g3_capacity_table_order_is_irrelevant_and_uniform_distribution_keeps_the_original_result() {
    let mesh = mesh();
    let fallback = ConductivityModel::isotropic_declared(999.0).unwrap();
    let assigned = materials(false);
    let values = [cap(15.0), cap(150.0)];
    let forward = [(MaterialId(7), values[0]), (MaterialId(2), values[1])];
    let reverse = [forward[1], forward[0]];
    with_cx(|cx| {
        let history = history(cx, &mesh);
        let mut outputs = Vec::new();
        for capacity in [ThermalHistoryCapacity::PerElement(&values),
            ThermalHistoryCapacity::ByMaterial(&forward), ThermalHistoryCapacity::ByMaterial(&reverse)] {
            let mut run = ThermalHistoryRun::new_with_capacity(cx,
                problem(&history, &fallback, Some(&assigned)), capacity, &[300.0; 8], 1,
                ThermalHistoryMethod::Linear).unwrap();
            outputs.push(run.advance_to(cx, 0.25, controls()).unwrap().temperature);
        }
        assert_eq!(outputs[0], outputs[1]);
        assert_eq!(outputs[0], outputs[2]);
        let mut original = ThermalHistoryRun::new(cx, problem(&history, &fallback, Some(&assigned)),
            cap(15.0), &[300.0; 8], 1).unwrap();
        let repeated = [cap(15.0); 2];
        let mut elementwise = ThermalHistoryRun::new_with_capacity(cx,
            problem(&history, &fallback, Some(&assigned)), ThermalHistoryCapacity::PerElement(&repeated),
            &[300.0; 8], 1, ThermalHistoryMethod::Linear).unwrap();
        assert_eq!(original.advance_to(cx, 0.25, controls()).unwrap().temperature,
            elementwise.advance_to(cx, 0.25, controls()).unwrap().temperature);
    });
}

#[test]
fn g1_nonlinear_material_bound_capacity_matches_the_direct_heterogeneous_stepper() {
    let mesh = mesh();
    let fallback = ConductivityModel::isotropic_declared(999.0).unwrap();
    let assigned = materials(true);
    let values = [cap(15.0), cap(150.0)];
    let rows = [(MaterialId(2), values[1]), (MaterialId(7), values[0])];
    with_cx(|cx| {
        let history = history(cx, &mesh);
        let mut old = vec![300.0, 315.0, 310.0, 330.0, 350.0, 360.0, 345.0, 340.0];
        let policy = NonlinearStepConfig::default();
        let mut run = ThermalHistoryRun::new_with_capacity(cx,
            problem(&history, &fallback, Some(&assigned)), ThermalHistoryCapacity::ByMaterial(&rows),
            &old, 3, ThermalHistoryMethod::Nonlinear(policy)).unwrap();
        let manual = BackwardEuler::per_element(cx, &mesh, &values).unwrap();
        let boundary = history.at(cx, 0.0).unwrap();
        let source = ScalarField::Uniform(30.0);
        let mut previous = 0.0;
        for time in [0.01, 0.025, 0.05] {
            let expected = manual.advance_nonlinear_prescribed(cx, ConductionProblem {
                mesh: &mesh, boundary: &boundary, material: &fallback,
                source: &source, element_materials: Some(&assigned),
            }, None, &old, time-previous, controls(), policy).unwrap();
            let ThermalHistoryStep::Nonlinear(actual) = run.advance_to_report(cx, time, controls()).unwrap()
                else { panic!("nonlinear evidence was lost"); };
            assert_eq!(actual.step.temperature, expected.step.temperature);
            assert_eq!(actual.step.stored_energy_change_j, expected.step.stored_energy_change_j);
            assert_eq!(actual.step.energy_residual_j, expected.step.energy_residual_j);
            assert_eq!(actual.nonlinear_iterations, expected.nonlinear_iterations);
            assert!(actual.residual_j <= actual.threshold_j);
            old = expected.step.temperature;
            previous = time;
        }
    });
}

#[test]
fn g1_changing_prescribed_temperatures_account_for_each_materials_stored_heat() {
    let mesh = mesh();
    let fallback = ConductivityModel::isotropic_declared(999.0).unwrap();
    let assigned = materials(false);
    let rows = [(MaterialId(2), cap(150.0)), (MaterialId(7), cap(15.0))];
    with_cx(|cx| {
        let row = BoundaryCondition { region: "surface".into(), physics: Physics::Thermal,
            kind: BcKind::Dirichlet, compatibility: None, frame: 0,
            value: Some(BcValue::Signal(TimeSignal::Ramp { t_start: 0.0, t_end: 1.0,
                from: QtyAny::new(300.0, TEMPERATURE_DIMS), to: QtyAny::new(400.0, TEMPERATURE_DIMS) })) };
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .remainder("surface", ThermalBc::dirichlet(300.0).unwrap()).unwrap(),
            &Environment::earth_lab(), vec![ThermalRegionHistory { row, reference_temperature: None }], 0.0).unwrap();
        let mut p = problem(&history, &fallback, Some(&assigned));
        p.source = HeatSourceHistory::Fixed(ScalarField::Uniform(0.0));
        let mut run = ThermalHistoryRun::new_with_capacity(cx, p, ThermalHistoryCapacity::ByMaterial(&rows),
            &[300.0; 8], 1, ThermalHistoryMethod::Linear).unwrap();
        let step = run.advance_to(cx, 0.5, controls()).unwrap();
        assert_eq!(step.temperature, vec![350.0; 8]);
        let expected_j = (15.0 + 150.0) * 50.0 / 6.0;
        assert!((step.stored_energy_change_j-expected_j).abs() < 1e-7);
        assert!((step.dirichlet_in_w*0.5-expected_j).abs() < 1e-7);
        assert!(step.energy_residual_j.abs() <= controls().energy_tolerance_j);
    });
}

#[test]
fn g0_missing_duplicate_unused_and_mis_sized_capacities_refuse_without_fallback() {
    let mesh = mesh();
    let fallback = ConductivityModel::isotropic_declared(20.0).unwrap();
    let assigned = materials(false);
    let initial = [300.0; 8];
    let history = with_cx(|cx| history(cx, &mesh));
    with_cx(|cx| {
        for values in [Vec::new(), vec![cap(15.0)], vec![cap(15.0); 3]] {
            assert!(ThermalHistoryRun::new_with_capacity(cx, problem(&history, &fallback, Some(&assigned)),
                ThermalHistoryCapacity::PerElement(&values), &initial, 1, ThermalHistoryMethod::Linear).is_err());
        }
        for rows in [vec![], vec![(MaterialId(7), cap(15.0))],
            vec![(MaterialId(7), cap(15.0)), (MaterialId(7), cap(150.0))],
            vec![(MaterialId(7), cap(15.0)), (MaterialId(2), cap(150.0)), (MaterialId(9), cap(10.0))]] {
            assert!(ThermalHistoryRun::new_with_capacity(cx, problem(&history, &fallback, Some(&assigned)),
                ThermalHistoryCapacity::ByMaterial(&rows), &initial, 1, ThermalHistoryMethod::Linear).is_err());
        }
        let valid = [(MaterialId(7), cap(15.0)), (MaterialId(2), cap(150.0))];
        assert!(ThermalHistoryRun::new_with_capacity(cx, problem(&history, &fallback, None),
            ThermalHistoryCapacity::ByMaterial(&valid), &initial, 1, ThermalHistoryMethod::Linear).is_err());
    });
    with_cancelled_cx(|cx| assert!(matches!(ThermalHistoryRun::new_with_capacity(cx,
        problem(&history, &fallback, Some(&assigned)), ThermalHistoryCapacity::PerElement(&[cap(15.0); 2]),
        &initial, 1, ThermalHistoryMethod::Linear), Err(ConductionError::Cancelled { .. }))));
    assert_eq!(initial, [300.0; 8]);
}
