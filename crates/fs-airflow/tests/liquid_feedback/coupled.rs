use super::*;
use fs_airflow::AirflowError;
use fs_airflow::conjugate::{ConjugateConfig, Relaxation, SolidRegionState};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::LiquidFeedbackError;
use fs_couple::iqn_ils::IqnIlsConfig;

fn config() -> ConjugateConfig {
    ConjugateConfig { temperature_tolerance_k: 1e-7, balance_tolerance_w: 1e-6,
        balance_relative_tolerance: 1e-9, max_iterations: 80,
        relaxation: Relaxation::Fixed { omega: 0.7 } }
}

#[test]
fn signed_power_closes_with_updated_materials_and_retained_original_rows() {
    let card = card(true, true, 1.0);
    let law = law(&card);
    let source = caloric(&card);
    with_cx(|cx| {
        for power in [20.0, -15.0] {
            let path = LiquidHeatPath::try_new(cx, &law, &source, MassFlowRate::new(0.001), Temperature::new(330.0),
                &[section("wall")], LiquidPropertyTemperature::BulkMean).unwrap();
            let mut calls = 0;
            let result = path.solve_coupled_iqn(cx, &[330.0], &config(), budget(), IqnIlsConfig::default(), |_, rows| {
                calls += 1;
                let h = rows[0].coefficient().value.value();
                let reference = rows[0].reference_temperature().value();
                let wall = reference + power / (h*0.008);
                Ok(LiquidSolidResponse { regions: vec![SolidRegionState {
                    region: "wall".into(), area_m2: 0.008, mean_wall_temperature_k: wall,
                    heat_rate_w: power, mean_reference_temperature_k: Some(reference),
                }], state: (calls, h, reference) })
            }).unwrap();
            assert_eq!(calls, result.iterations);
            assert_eq!(result.solid.state.0, result.iterations);
            assert_eq!(result.solid.state.1, result.boundary_state.regions[0].boundary.coefficient().value.value());
            assert_eq!(result.solid.state.2, result.boundary_state.regions[0].boundary.reference_temperature().value());
            assert!((result.fluid.total_heat_rate_w - power).abs() <= config().balance_tolerance_w);
            let outlet_h = path.inlet().specific_enthalpy_j_kg() + power / path.mass_flow_kg_per_s();
            let expected = source.curve().state_at_specific_enthalpy(outlet_h).unwrap();
            close(result.fluid.regions[0].outlet.temperature_k(), expected.temperature_k());
            assert!(result.maximum_coefficient_residual_relative <= budget().coefficient_relative_tolerance);
            assert!(result.maximum_heat_imbalance_w <= config().balance_tolerance_w);
        }
    });
}

#[test]
fn production_fem_consumes_temperature_updated_robin_rows_and_returns_its_field() {
    use fs_conduction::fixtures::{box_grid, on_box_face};
    use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, InitialGuess,
        LinearConfig, Nonlinearity, ScalarField, SolveConfig, StopRule, ThermalBoundaryBuilder, solve};
    use fs_rep_mesh::TetComplex;

    let card = card(true, true, 1.0);
    let law = law(&card);
    let source = caloric(&card);
    let (complex, positions) = box_grid([2, 2, 2], [0.1, 0.08, 0.005]);
    let mesh = ConductionMesh::new(TetComplex::from_tets(positions.len(), complex.tets), positions).unwrap();
    let solid_material = ConductivityModel::isotropic_declared(8.0).unwrap();
    let power = ScalarField::nodal("20W", mesh.vertex_count(), vec![20.0 / (0.1*0.08*0.005); mesh.vertex_count()]).unwrap();
    with_cx(|cx| {
        let path = LiquidHeatPath::try_new(cx, &law, &source, MassFlowRate::new(0.001), Temperature::new(300.0),
            &[section("wall")], LiquidPropertyTemperature::BulkMean).unwrap();
        let result = path.solve_coupled_iqn(cx, &[320.0], &config(), budget(), IqnIlsConfig::default(), |cx, rows| {
            let boundary = ThermalBoundaryBuilder::new(&mesh).region("wall", |face|
                face.outward_normal[2] > 0.99 && on_box_face(face.centroid[2], 0.005),
                rows[0].boundary_condition().clone()).unwrap().adiabatic_remainder().finish().unwrap();
            let solution = solve(cx, ConductionProblem { mesh: &mesh, boundary: &boundary, material: &solid_material,
                source: &power, element_materials: None }, SolveConfig {
                    nonlinearity: Nonlinearity::FixedPoint { relaxation: 1.0, max_backtracks: 8 },
                    stop: StopRule { residual_rtol: 1e-10, residual_atol: 1e-20, step_atol: 0.0, max_iterations: 12 },
                    linear: LinearConfig { tolerance: 1e-10, max_iterations: 6000, restart: 80 },
                    initial: InitialGuess::Uniform(320.0),
                }).unwrap();
            let regions = solution.report.robin_fluxes.iter().map(SolidRegionState::from_robin_flux).collect();
            Ok(LiquidSolidResponse { regions, state: solution })
        }).unwrap();
        assert!((result.solid.state.report.energy.robin_out_w - 20.0).abs() < 1e-5);
        assert!((result.fluid.total_heat_rate_w - 20.0).abs() < 1e-5);
        let material = &result.fluid.regions[0].material;
        assert_ne!(material.thermal_conductivity_w_m_k(), law.reference_state().thermal_conductivity_w_m_k());
        close(result.fluid.regions[0].boundary.coefficient().value.value(), 3.66*material.thermal_conductivity_w_m_k()/0.01);
        assert!(result.maximum_wall_residual_k <= config().temperature_tolerance_k);
    });
}

#[test]
fn equal_and_opposite_local_heat_errors_do_not_cancel() {
    let card = card(true, true, 1.0);
    let law = law(&card);
    let source = caloric(&card);
    with_cx(|cx| {
        let walls = [350.0, 320.0];
        let names = ["first", "last"];
        let path = LiquidHeatPath::try_new(cx, &law, &source, MassFlowRate::new(0.001), Temperature::new(300.0),
            &[section(names[0]), section(names[1])], LiquidPropertyTemperature::BulkMean).unwrap();
        let error = path.solve_coupled_iqn(cx, &walls, &config(), budget(), IqnIlsConfig::default(), |_, rows| {
            let regions = rows.iter().enumerate().map(|(i,row)| {
                let reference = row.reference_temperature().value();
                SolidRegionState { region: names[i].into(), area_m2: 0.008, mean_wall_temperature_k: walls[i],
                    heat_rate_w: row.coefficient().value.value()*0.008*(walls[i]-reference) + if i==0 {0.1} else {-0.1},
                    mean_reference_temperature_k: Some(reference) }
            }).collect();
            Ok(LiquidSolidResponse { regions, state: () })
        }).unwrap_err();
        assert!(matches!(error, LiquidFeedbackError::Enthalpy(
            fs_airflow::graph::thermal::medium::series::enthalpy::LiquidEnthalpyError::Flow(
                fs_airflow::graph::thermal::medium::MediumFlowError::Airflow(AirflowError::ConjugateBalanceUnclosed { .. })
            ))));
    });
}

#[test]
fn bad_controls_precede_solid_work_and_bad_responses_are_not_published() {
    let card = card(true, true, 1.0);
    let law = law(&card);
    let source = caloric(&card);
    with_cx(|cx| {
        let path = LiquidHeatPath::try_new(cx, &law, &source, MassFlowRate::new(0.001), Temperature::new(300.0),
            &[section("wall")], LiquidPropertyTemperature::BulkMean).unwrap();
        let mut calls = 0;
        let no_work = path.solve_coupled_iqn(cx, &[320.0], &ConjugateConfig { max_iterations: 0, ..config() },
            budget(), IqnIlsConfig::default(), |_, _| { calls += 1; Ok(LiquidSolidResponse { regions: vec![], state: () }) });
        assert!(no_work.is_err());
        assert_eq!(calls, 0);
        let wrong = path.solve_coupled_iqn(cx, &[320.0], &config(), budget(), IqnIlsConfig::default(), |_, rows| {
            Ok(LiquidSolidResponse { regions: vec![SolidRegionState {
                region: "wall".into(), area_m2: 0.008, mean_wall_temperature_k: 320.0, heat_rate_w: 1.0,
                mean_reference_temperature_k: Some(rows[0].reference_temperature().value()+1.0),
            }], state: () })
        });
        assert!(wrong.is_err());
    });
}
