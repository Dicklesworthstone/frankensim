//! The public sourced heat path agrees with an independent backward-Euler
//! solution on a thin, resolved liquid layer above its explicit thermal CFL.
use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::FreeSurfaceImplicitSolve;

#[test]
fn g1_implicit_height_graph_prescribed_heat_matches_independent_stiff_layer() {
    const FILL: f64 = 0.002;
    let card = source_card(1e-6);
    let source = source(&card);
    let law = LiquidTemperatureLaw::try_from_material(
        &card,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let caloric = caloric(&card);
    let gas = gas_source::model(0.0, 250000.0, 1000.0);
    let point = gas_source::point(450.0, 100000.0);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut cells = graph_cells(false);
        for cell in &mut cells {
            if cell.coordinates[2] == 3 {
                cell.fill_fraction = FILL;
            }
        }
        let mut thermal = plan_free_surface_from_material(&source, DX, DT, 0.1)
            .unwrap()
            .free_surface3_from_cells_with_height_graph_thermocapillarity(
                domain(),
                &cells,
                &law,
                &caloric,
                FreeSurfaceEnthalpyControls::default(),
                graph(),
                None,
                None,
                |_| 200000.0,
                cx,
                &gate,
            )
            .unwrap();
        // Distances from the uniform occupied layer's centroid to its plane,
        // and from that plane to the next empty gas cell's centre.
        let distances = [0.5 * FILL * DX, (1.5 - FILL) * DX];
        let capacity = 1000.0 * FILL * DX.powi(3) * 2000.0;
        let gl = 0.5 * DX * DX / distances[0];
        let kg = 0.026 * 1000.0 * (0.5 + 99000.0 / 249000.0);
        let gg = kg * DX * DX / distances[1];
        let conductance = gl * gg / (gl + gg);
        assert!(DT > 2.0 * capacity / conductance);
        let controls = FreeSurfaceConductionControls {
            maximum_faces: 64,
            maximum_substeps: 1,
            maximum_dt_s: DT,
            energy_tolerance_j: 1e-9,
        };
        let mut expected_t = 400.0;
        for interval in 0..2 {
            let intervals = interval + 1;
            let interval_controls = FreeSurfaceConductionControls {
                maximum_substeps: intervals,
                maximum_dt_s: DT / intervals as f64,
                ..controls
            };
            let contacts: Vec<_> = thermal
                .height_graph_interface_patches(&gate)
                .unwrap()
                .into_iter()
                .map(|patch| FreeSurfaceInterfaceContact {
                    surface: patch.with_exchange(distances, 0.0, 400.0),
                    model: &gas,
                    context: &point,
                })
                .collect();
            let before = thermal.total_enthalpy_j().unwrap();
            let mass = thermal.flow().ledger_mass_kg().unwrap();
            if interval == 0 {
                // The unchanged explicit mode must honour its physical CFL
                // and refuse this one-substep budget without publishing heat.
                let populations = thermal.flow().lattice().grid().state_bits();
                assert!(
                    thermal
                        .step_with_surface_exchange(
                            &[],
                            &[],
                            controls,
                            &mut [],
                            &[],
                            &[],
                            &contacts,
                            FreeSurfaceInterfaceSolve::default(),
                            false,
                            &gate,
                        )
                        .is_err()
                );
                assert_eq!(
                    thermal.total_enthalpy_j().unwrap().to_bits(),
                    before.to_bits()
                );
                assert_eq!(populations, thermal.flow().lattice().grid().state_bits());
            }
            let solve = FreeSurfaceInterfaceSolve {
                implicit: Some(FreeSurfaceImplicitSolve::default()),
                ..FreeSurfaceInterfaceSolve::default()
            };
            let report = thermal
                .step_with_surface_exchange(
                    &[],
                    &[],
                    interval_controls,
                    &mut [],
                    &[],
                    &[],
                    &contacts,
                    solve,
                    false,
                    &gate,
                )
                .unwrap()
                .conduction;
            let old_t = expected_t;
            for _ in 0..intervals {
                let dt = DT / intervals as f64;
                expected_t = (capacity * expected_t + dt * conductance * 450.0)
                    / (capacity + dt * conductance);
            }
            let expected_heat = 64.0 * capacity * (expected_t - old_t);
            let surface = report.surface.as_ref().unwrap();
            let work = surface.implicit.as_ref().unwrap();
            assert_eq!(work.intervals, intervals);
            assert_eq!(report.thermal.substeps, intervals);
            assert!(work.residual_evaluations < 24 * intervals);
            assert!(work.maximum_cell_residual_j < 1e-10);
            near(surface.heat_from_gas_j, expected_heat, 1e-8);
            near(report.flow.external_heat_j, expected_heat, 1e-8);
            assert!((thermal.total_enthalpy_j().unwrap() - before - expected_heat).abs() < 1e-9);
            near(thermal.flow().ledger_mass_kg().unwrap(), mass, 3e-12);
            assert_eq!(surface.heat_from_radiation_j, 0.0);
            let expected_ts = (gl * expected_t + gg * 450.0) / (gl + gg);
            let sigma = 0.001 + 1e-6 * (expected_ts - 400.0);
            for sample in &surface.thermal_endpoint {
                near(sample.liquid_temperature_k, expected_t, 2e-10);
                near(sample.surface_temperature_k, expected_ts, 2e-10);
                near(
                    sample.surface_tension.as_ref().unwrap().value_si(),
                    sigma,
                    2e-10,
                );
            }
            for sample in &surface.native_endpoint {
                near(sample.surface_temperature_k, expected_ts, 2e-10);
                let [x, y, z] = sample.surface.face.liquid;
                near(
                    thermal
                        .flow()
                        .lattice()
                        .surface_tension_at(x, y, z)
                        .unwrap(),
                    sigma * DT * DT / (1000.0 * DX.powi(3)),
                    2e-10,
                );
                assert!(
                    thermal
                        .flow()
                        .lattice()
                        .force_density_at(x, y, z)
                        .unwrap()
                        .iter()
                        .all(|f| f.abs() < 1e-15)
                );
            }
        }
    });
}
