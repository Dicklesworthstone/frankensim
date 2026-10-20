use super::*;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::loss::{
    ScanObservationWithLoss, ScanResidualLoss};

#[test]
fn g3_value_only_matches_physical_scores_and_preserves_shared_work_accounting() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .remainder("surface", ThermalBc::adiabatic()).unwrap(),
            &Environment::earth_lab(), vec![], 0.).unwrap();
        for loss in [ScanResidualLoss::Quadratic, ScanResidualLoss::Huber { threshold: 0.5 }] {
            let rows: Vec<_> = observations(cx, &m).into_iter()
                .map(|o| ScanObservationWithLoss::new(o, loss).unwrap()).collect();
            let mut study = ScanPowerStudy::new(cx, sources([2., 3., 0.]), controls(),
                objective(rows, false), factory(&boundary, &k, &c), 2).unwrap();
            let point = study.initial_coordinates().to_vec();
            let valued = study.evaluate_value(cx, &point).unwrap();
            let differentiated = study.evaluate(cx, &point).unwrap();
            assert_eq!(valued.physical.value, differentiated.physical.value);
            assert_eq!(valued.physical.predictions, differentiated.physical.predictions);
            assert_eq!(valued.physical.standardized_residuals, differentiated.physical.standardized_residuals);
            assert_eq!(valued.physical.final_step.solution().specific_enthalpy_j_kg,
                differentiated.physical.final_step.solution().specific_enthalpy_j_kg);
            assert_eq!(valued.powers_w, differentiated.powers_w);
            assert_eq!(valued.physical.endpoint_count, 2);
            assert_eq!(study.evaluations(), 2);
            assert!(study.evaluate_value(cx, &point).is_err());
            assert_eq!(study.evaluations(), 2);
            study.extend_evaluation_budget(3).unwrap();
            assert_eq!(study.evaluate_value(cx, &point).unwrap().physical.value, valued.physical.value);
        }
    });
}

#[test]
fn g1_value_only_accepts_a_physical_phase_corner_without_fabricating_a_gradient() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|_, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .remainder("surface", ThermalBc::adiabatic()).unwrap(),
            &Environment::earth_lab(), vec![], 0.).unwrap();
        let corner_factory = |cx: &Cx<'_>, sources| {
            let run = EnthalpyHistoryRun::new(cx, ThermalHistoryProblem {
                boundary: &boundary, material: &k, element_materials: None, interfaces: None,
                source: HeatSourceHistory::Fixed(ScalarField::Uniform(0.)) },
                EnthalpyHistoryMaterials { table: &[ReferenceEnthalpyMaterial {
                    curve: &c, reference_density_kg_m3: 1. }], of_element: &[0; 6], phase_conductivity: None },
                &[100.; 8], EnthalpyHistoryBudget {
                    spatial: EnthalpyBudget { max_vertices: 8, max_elements: 6 }, max_steps: 2 })?;
            Ok(ScanningEnthalpyRun::new(cx, run, sources, 1, None)?)
        };
        let rows = vec![ScanObservation { endpoint: 1, vertex: 0,
            observable: ScanObservable::LiquidMassFraction, target: 0.5, scale: 0.1 }];
        let mut study = ScanPowerStudy::new_joint_with_power_coordinates(cx, sources([0.; 3]),
            controls(), vec![], ScanPowerCoordinates::NonnegativeLinear,
            objective(rows, false), corner_factory, 3).unwrap();
        let point = study.initial_coordinates().to_vec();
        let valued = study.evaluate_value(cx, &point).unwrap();
        assert_eq!(valued.physical.predictions, vec![0.]);
        assert_eq!(valued.physical.value, 12.5);
        assert!(valued.physical.final_step.solution().temperature.iter().all(|&t| t == 310.));
        assert!(study.evaluate(cx, &point).is_err());
        assert_eq!(study.evaluations(), 2);
        assert_eq!(study.evaluate_value(cx, &point).unwrap().physical.value, 12.5);
    });
}

#[test]
fn g0_value_only_needs_no_tangent_budget_but_preserves_input_and_cancel_refusals() {
    let m = mesh(); let k = ConductivityModel::isotropic_declared(2.).unwrap(); let c = curve();
    with_cx(|gate, cx| {
        let boundary = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&m)
            .remainder("surface", ThermalBc::adiabatic()).unwrap(),
            &Environment::earth_lab(), vec![], 0.).unwrap();
        let mut spec = objective(observations(cx, &m), false);
        spec.budget.max_retained_vertices = 0;
        spec.budget.max_retained_source_vertices = 0;
        let mut study = ScanPowerStudy::new(cx, sources([2., 3., 0.]), controls(), spec,
            factory(&boundary, &k, &c), 3).unwrap();
        let point = study.initial_coordinates().to_vec();
        assert!(study.evaluate_value(cx, &point).unwrap().physical.value.is_finite());
        assert!(study.evaluate(cx, &point).is_err());
        assert!(study.evaluate_value(cx, &[f64::NAN, 0.]).is_err());
        assert_eq!(study.evaluations(), 2);
        gate.request();
        assert!(study.evaluate_value(cx, &point).is_err());
        assert_eq!(study.evaluations(), 2);
    });
}
