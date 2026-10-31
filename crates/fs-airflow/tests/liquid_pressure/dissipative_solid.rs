//! Hydraulic heat, live wall feedback and native FEM share one operating point.
use super::*;
use fs_airflow::conjugate::SolidRegionState;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::solid::PressureSolidControls;
use fs_convection::CorrelationRobinBoundary;

fn controls() -> PressureSolidControls {
    let mut controls = solid_controls();
    controls.pressure.mass_flow_bracket_kg_per_s = (0.0008, 0.0013);
    controls.pressure.maximum_neglected_hydraulic_power_w = 0.0;
    controls.pressure.pressure_relative_tolerance = 1e-8;
    controls
}

#[test]
fn pressure_heating_and_fem_wall_heat_close_on_the_same_material_trajectory() {
    use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel,
        InitialGuess, LinearConfig, Nonlinearity, ScalarField, SolveConfig,
        StopRule, ThermalBoundaryBuilder, solve};
    use fs_conduction::fixtures::{box_grid, on_box_face};
    use fs_rep_mesh::TetComplex;

    for variable in [false, true] {
        with_path(variable, 1.0, LiquidPropertyTemperature::FilmMean, |cx, path| {
            let (topology, positions) = box_grid([2, 1, 1], [0.06, 0.1, 0.01]);
            let mesh = ConductionMesh::new(
                TetComplex::from_tets(positions.len(), topology.tets), positions,
            ).unwrap();
            let material = ConductivityModel::isotropic_declared(5.0).unwrap();
            let source = ScalarField::nodal("uniform 20W", mesh.vertex_count(),
                vec![20.0 / (0.06 * 0.1 * 0.01); mesh.vertex_count()]).unwrap();
            let controls = controls();
            let mut calls = 0;
            let answer = path.solve_pressure_coupled_iqn_with_dissipation(
                cx, &[325.0, 330.0], &heated_resistances(5e15),
                Pressure::new(2e7), &controls, |cx, rows| {
                    calls += 1;
                    let mut builder = ThermalBoundaryBuilder::new(&mesh);
                    for (index, name) in ["first", "last"].iter().enumerate() {
                        builder = builder.region(name, |face| {
                            face.outward_normal[2] > 0.99
                                && on_box_face(face.centroid[2], 0.01)
                                && if index == 0 {
                                    face.centroid[0] < 0.03
                                } else { face.centroid[0] >= 0.03 }
                        }, rows[index].boundary_condition().clone()).unwrap();
                    }
                    let boundary = builder.adiabatic_remainder().finish().unwrap();
                    // Each speculative callback solves from the same original
                    // body/source, with no accepted state or source debit changed.
                    let solution = solve(cx, ConductionProblem {
                        mesh: &mesh, boundary: &boundary, material: &material,
                        source: &source, element_materials: None,
                    }, SolveConfig {
                        nonlinearity: Nonlinearity::FixedPoint {
                            relaxation: 1.0, max_backtracks: 8,
                        },
                        stop: StopRule {
                            residual_rtol: 1e-10, residual_atol: 1e-18,
                            step_atol: 0.0, max_iterations: 12,
                        },
                        linear: LinearConfig {
                            tolerance: 1e-10, max_iterations: 10000, restart: 64,
                        },
                        initial: InitialGuess::Uniform(330.0),
                    }).unwrap();
                    let regions = solution.report.robin_fluxes.iter()
                        .map(SolidRegionState::from_robin_flux).collect();
                    Ok(LiquidSolidResponse { regions, state: solution })
                },
            ).unwrap();
            assert_eq!(calls, answer.total_solid_solves);
            assert!(calls <= controls.max_total_solid_solves);
            let operating = answer.operating_point;
            assert!(operating.pressure_residual_pa.abs() <= 0.21);
            let mass = operating.mass_flow_kg_per_s;
            let hydraulic = operating.hydraulics.total_dissipated_power_w;
            let result = operating.thermal;
            let outlet = result.fluid.regions[1].outlet;
            let gain = mass * (outlet.specific_enthalpy_j_kg()
                - result.fluid.regions[0].inlet.specific_enthalpy_j_kg());
            // Independent whole-body balance: 20W enters the real FEM source,
            // and dp*Q enters only the liquid. Neither heat is counted twice.
            assert!((result.solid.state.report.energy.robin_out_w - 20.0).abs() < 2e-5);
            assert!((result.fluid.total_heat_rate_w - 20.0).abs() < 2e-5);
            close(hydraulic, operating.hydraulics.total_pressure_drop_pa * mass / 1000.0);
            assert!((gain - 20.0 - hydraulic).abs() < 3e-5);
            assert!(hydraulic > 15.0);
            assert!(result.maximum_heat_imbalance_w <= controls.coupling.balance_tolerance_w);
            if variable {
                assert!(mass > 0.001);
            } else {
                // Constant source viscosity gives dp = 2 * K * mu * mdot/rho.
                assert!((mass / 0.001 - 1.0).abs() < 2e-8);
                assert!((outlet.temperature_k() - 320.0).abs() < 2e-5);
            }
            for (((sent, fluid), solid), hydraulic) in result.boundary_state.regions.iter()
                .zip(&result.fluid.regions).zip(&result.solid.regions)
                .zip(&operating.hydraulics.regions)
            {
                assert_eq!(sent.boundary.reference_temperature().value(),
                    solid.mean_reference_temperature_k.unwrap());
                assert!((solid.heat_rate_w - fluid.heat_rate_w).abs() <= 1e-5);
                close(fluid.heat_rate_w + hydraulic.loss.dissipated_power_w, mass *
                    (fluid.outlet.specific_enthalpy_j_kg() - fluid.inlet.specific_enthalpy_j_kg()));
                let bulk = 0.5 * (fluid.inlet.temperature_k() + fluid.outlet.temperature_k());
                close(hydraulic.material.dynamic_viscosity_pa_s(),
                    if variable { 0.005 - 1e-5 * bulk } else { 0.002 });
                assert!(hydraulic.material.material().property("dynamic_viscosity")
                    .unwrap().answer().receipt.query_point.contains(&("pressure".into(), 100000.0)));
                if variable {
                    assert_ne!(hydraulic.material.dynamic_viscosity_pa_s(),
                        fluid.material.dynamic_viscosity_pa_s());
                }
            }
        });
    }
}

// A held-wall response is sufficient for work/refusal injections. The positive
// coupled-physics regression above uses the original conduction FEM owner.
fn held_wall(rows: &[CorrelationRobinBoundary], errors: [f64; 2]) -> LiquidSolidResponse<()> {
    LiquidSolidResponse { state: (), regions: rows.iter().enumerate().map(|(i, row)| {
        let reference = row.reference_temperature().value();
        SolidRegionState {
            region: ["first", "last"][i].into(), area_m2: 0.003,
            mean_wall_temperature_k: 300.0,
            heat_rate_w: row.coefficient().value.value() * 0.003 * (300.0 - reference) + errors[i],
            mean_reference_temperature_k: Some(reference),
        }
    }).collect() }
}

#[test]
fn dissipative_solid_rejects_nonzero_omission_before_any_callback() {
    with_path(false, 1.0, LiquidPropertyTemperature::BulkMean, |cx, path| {
        let mut controls = controls();
        controls.pressure.maximum_neglected_hydraulic_power_w = 1.0;
        let result = path.solve_pressure_coupled_iqn_with_dissipation::<(), _>(
            cx, &[300.0, 300.0], &heated_resistances(5e15), Pressure::new(2e7),
            &controls, |_, _| panic!("invalid model must refuse before the solid callback"),
        );
        assert!(matches!(result, Err(PressureDrivenError::Invalid(_))));
    });
}

#[test]
fn dissipative_solid_keeps_one_work_budget_across_pressure_trials() {
    with_path(false, 1.0, LiquidPropertyTemperature::BulkMean, |cx, path| {
        let mut controls = controls();
        controls.max_total_solid_solves = 1;
        let mut calls = 0;
        let result = path.solve_pressure_coupled_iqn_with_dissipation(
            cx, &[300.0, 300.0], &heated_resistances(5e15), Pressure::new(2e7),
            &controls, |_, rows| {
                calls += 1;
                Ok(held_wall(rows, [0.0, 0.0]))
            },
        );
        // The first complete heated-flow/solid candidate converges. The next
        // outer pressure trial must not reset the global callback counter.
        assert_eq!(calls, 1);
        assert!(matches!(result,
            Err(PressureDrivenError::SolidSolveBudget { performed: 1, maximum: 1 })));
    });
}

#[test]
fn opposing_regional_heat_errors_cannot_cancel_after_hydraulic_deposition() {
    with_path(false, 1.0, LiquidPropertyTemperature::BulkMean, |cx, path| {
        let mut calls = 0;
        let result = path.solve_pressure_coupled_iqn_with_dissipation(
            cx, &[300.0, 300.0], &heated_resistances(5e15), Pressure::new(2e7),
            &controls(), |_, rows| {
                calls += 1;
                Ok(held_wall(rows, [1.0, -1.0]))
            },
        );
        assert_eq!(calls, 1);
        assert!(matches!(result, Err(PressureDrivenError::Thermal(_))));
    });
}

#[test]
fn variable_density_cannot_enter_dissipative_solid_coupling() {
    let original = card(false, 1.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "density" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, 1025.0), (450.0, 925.0)],
                dims: Dims([-3, 1, 0, 0, 0, 0]),
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    with_card_path(&card, LiquidPropertyTemperature::BulkMean, |cx, path| {
        let result = path.solve_pressure_coupled_iqn_with_dissipation::<(), _>(
            cx, &[300.0, 300.0], &heated_resistances(5e15), Pressure::new(2e7),
            &controls(), |_, _| panic!("unsupported density must refuse before a solid callback"),
        );
        assert!(result.is_err());
    });
}
