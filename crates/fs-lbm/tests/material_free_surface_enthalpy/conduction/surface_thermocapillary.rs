//! Physical surface-temperature gradients drive the original native stress.
//! Synthetic coefficients check implementation, not experimental calibration.
use super::*;

fn patches() -> Vec<FreeSurfaceInterfaceFace> {
    (0..8)
        .flat_map(|y| {
            (0..8).map(move |x| FreeSurfaceInterfaceFace {
                face: FreeSurfaceGasFace {
                    liquid: [x, y, 3],
                    gas: [x, y, 4],
                    area_m2: DX * DX,
                    distances_m: [0.25 * DX, 0.5 * DX],
                },
                emissivity: 0.0,
                surroundings_k: 400.0,
            })
        })
        .collect()
}

fn surface_controls() -> FreeSurfaceConductionControls {
    FreeSurfaceConductionControls {
        maximum_faces: 64,
        ..controls()
    }
}

fn gas_k(scale: f64) -> f64 {
    0.026 * scale * (0.5 + 99000.0 / 249000.0)
}

fn contact_field<'a>(
    gas: &'a fs_material::thermochemical_gas::ThermochemicalGasModel,
    points: &'a [QueryPoint],
    surfaces: &[FreeSurfaceInterfaceFace],
) -> Vec<FreeSurfaceInterfaceContact<'a>> {
    surfaces
        .iter()
        .zip(points)
        .map(|(&surface, context)| FreeSurfaceInterfaceContact {
            surface,
            model: gas,
            context,
        })
        .collect()
}

fn series_surface(tl: f64, tg: f64, kg: f64) -> f64 {
    // Original authored liquid source k(T)=T-200, in an occupied half-cell.
    let gl = (tl - 200.0) * DX * DX / (0.25 * DX);
    let gg = kg * DX * DX / (0.5 * DX);
    (gl * tl + gg * tg) / (gl + gg)
}

#[test]
fn g1_nonuniform_gas_drives_surface_sigma_native_traction_and_pressure_once() {
    let slope = -0.00001;
    let card = thermal_card(slope, 500.0);
    let gas = gas_source::model(0.0, 250000.0, 1000.0);
    let surfaces = patches();
    let points: Vec<_> = surfaces
        .iter()
        .map(|p| gas_source::point(350.0 + 10.0 * p.face.liquid[0] as f64, 100000.0))
        .collect();
    let contacts = contact_field(&gas, &points, &surfaces);
    with_planar(
        &card,
        [0.0; 3],
        plane(),
        FreeSurfaceEnthalpyControls::default(),
        |thermal, gate| {
            let before = thermal.total_enthalpy_j().unwrap();
            let report = thermal
                .step_with_surface_exchange(
                    &[],
                    &[],
                    surface_controls(),
                    &mut [],
                    &[],
                    &[],
                    &contacts,
                    FreeSurfaceInterfaceSolve::default(),
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            let exchange = report.surface.as_ref().unwrap();
            assert_eq!(report.thermal.substeps, 4);
            assert_eq!(exchange.native_endpoint.len(), 64);
            let mut expected_heat = 0.0;
            for x in 0..8 {
                let tg = 350.0 + 10.0 * x as f64;
                let mut tl = 400.0;
                for _ in 0..4 {
                    let ts = series_surface(tl, tg, gas_k(1000.0));
                    let gl = (tl - 200.0) * DX * DX / (0.25 * DX);
                    let heat = 0.00025 * gl * (ts - tl);
                    tl += heat / 0.001; // rho * occupied-volume * Cp = 0.001 J/K
                    expected_heat += 8.0 * heat;
                }
                let sample = exchange
                    .thermal_endpoint
                    .iter()
                    .find(|s| s.surface.face.liquid == [x, 3, 3])
                    .unwrap();
                near(sample.liquid_temperature_k, tl, 2e-11);
                near(
                    sample.surface_temperature_k,
                    series_surface(tl, tg, gas_k(1000.0)),
                    2e-11,
                );
            }
            near(exchange.heat_from_gas_j, expected_heat, 2e-8);
            assert_eq!(exchange.heat_from_radiation_j, 0.0);
            near(report.flow.external_heat_j, expected_heat, 2e-8);
            near(
                thermal.total_enthalpy_j().unwrap() - before,
                expected_heat,
                2e-8,
            );

            let endpoint = |x| {
                exchange
                    .native_endpoint
                    .iter()
                    .find(|s| s.surface.face.liquid == [x, 3, 3])
                    .unwrap()
            };
            let sigma = |x| {
                let s = endpoint(x);
                let ts = series_surface(
                    s.liquid_temperature_k,
                    350.0 + 10.0 * x as f64,
                    gas_k(1000.0),
                );
                near(s.surface_temperature_k, ts, 2e-11);
                let value = 0.01 + slope * (ts - 400.0);
                near(s.surface_tension.as_ref().unwrap().value_si(), value, 2e-12);
                value
            };
            let native = thermal.flow().lattice();
            let fill = native.sample_at(3, 3, 3).unwrap().fill_fraction;
            let expected_force =
                (sigma(4) - sigma(2)) / (2.0 * DX) / (fill * DX) / 1000.0 * DT * DT / DX;
            let actual = native.force_density_at(3, 3, 3).unwrap();
            near(actual[0], expected_force, 2e-9);
            assert!(actual[0] < 0.0);
            assert_eq!(actual[2], 0.0);
            near(
                native.surface_tension_at(3, 3, 3).unwrap(),
                sigma(3) * DT * DT / (1000.0 * DX * DX * DX),
                2e-12,
            );
            assert!(
                thermal
                    .flow()
                    .sample_at_si(3, 3, 3)
                    .unwrap()
                    .unwrap()
                    .velocity_m_s[0]
                    < 0.0
            );

            let state = thermal.material_at(SURFACE).unwrap();
            let bulk_t = 300.0 + thermal.specific_enthalpy_at(SURFACE).unwrap() / 2000.0;
            near(endpoint(3).liquid_temperature_k, bulk_t, 2e-12);
            near(
                state
                    .material()
                    .property(SURFACE_TENSION_PROPERTY)
                    .unwrap()
                    .value_si(),
                0.01 + slope * (bulk_t - 400.0),
                2e-12,
            );
            assert!((bulk_t - endpoint(3).surface_temperature_k).abs() > 0.1);
            assert!(
                exchange
                    .native_endpoint
                    .iter()
                    .zip(&exchange.thermal_endpoint)
                    .any(|(a, b)| a.liquid_temperature_k != b.liquid_temperature_k),
                "advected bulk state must refresh the accepted surface field"
            );
        },
    );
}

#[test]
fn g1_uniform_solved_surface_has_no_marangoni_despite_bulk_temperature_gradient() {
    let card = thermal_card(-0.00001, 500.0);
    let gas = gas_source::model(0.0, 250000.0, 100000.0);
    let surfaces = patches();
    // Choose independent bath temperatures so each interface is exactly 400 K.
    // Authored heat offsets its gas loss and holds the nonuniform bulk field.
    let gg = gas_k(100000.0) * DX * DX / (0.5 * DX);
    let points: Vec<_> = surfaces
        .iter()
        .map(|p| {
            let tl = 400.0 + 3.0 * p.face.liquid[0] as f64;
            let gl = (tl - 200.0) * DX * DX / (0.25 * DX);
            gas_source::point(400.0 - gl * (tl - 400.0) / gg, 100000.0)
        })
        .collect();
    let heat: Vec<_> = surfaces
        .iter()
        .map(|p| {
            let tl = 400.0 + 3.0 * p.face.liquid[0] as f64;
            let gl = (tl - 200.0) * DX * DX / (0.25 * DX);
            FreeSurfaceCellHeat {
                coordinates: p.face.liquid,
                energy_j: DT * gl * (tl - 400.0),
            }
        })
        .collect();
    let contacts = contact_field(&gas, &points, &surfaces);
    with_planar(
        &card,
        [3.0, 0.0, 0.0],
        plane(),
        FreeSurfaceEnthalpyControls::default(),
        |thermal, gate| {
            assert!(thermal.flow().lattice().force_density_at(3, 3, 3).unwrap()[0].abs() > 1e-6);
            let before = thermal.total_enthalpy_j().unwrap();
            let report = thermal
                .step_with_surface_exchange(
                    &[],
                    &heat,
                    surface_controls(),
                    &mut [],
                    &[],
                    &[],
                    &contacts,
                    FreeSurfaceInterfaceSolve::default(),
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            let exchange = report.surface.as_ref().unwrap();
            for sample in &exchange.native_endpoint {
                near(sample.surface_temperature_k, 400.0, 2e-11);
                near(
                    sample.surface_tension.as_ref().unwrap().value_si(),
                    0.01,
                    2e-12,
                );
            }
            let force = thermal.flow().lattice().force_density_at(3, 3, 3).unwrap();
            assert!(force.iter().all(|f| f.abs() < 1e-11), "{force:?}");
            let speed = thermal
                .flow()
                .sample_at_si(3, 3, 3)
                .unwrap()
                .unwrap()
                .velocity_m_s;
            assert!(speed.iter().all(|u| u.abs() < 1e-9), "{speed:?}");
            near(
                thermal
                    .material_at(SURFACE)
                    .unwrap()
                    .material()
                    .property(SURFACE_TENSION_PROPERTY)
                    .unwrap()
                    .value_si(),
                0.01 - 0.00001 * 9.0,
                2e-12,
            );
            assert!((thermal.total_enthalpy_j().unwrap() - before).abs() < 1e-9);
            assert!(report.flow.external_heat_j.abs() < 1e-9);
        },
    );
}

#[test]
fn g4_surface_sigma_domain_complete_field_and_mode_publish_atomically_with_gas() {
    let card = thermal_card(-0.00001, 401.0);
    with_planar(
        &card,
        [0.0; 3],
        plane(),
        FreeSurfaceEnthalpyControls::default(),
        |thermal, gate| {
            let mut gas = IsochoricGasParcel::new(
                gas_source::model(100000.0, 250000.0, 100.0),
                gas_source::point(400.0, 100000.0),
                1e-6,
            )
            .unwrap();
            let gas_controls = GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 690.0,
                energy_tolerance_j_kg: 1e-7,
                maximum_iterations: 128,
            };
            let surfaces = patches();
            let before_bits = thermal.flow().lattice().grid().state_bits();
            let before_energy = thermal.total_enthalpy_j().unwrap();
            let before_gas = gas.state().clone();
            let before_source = thermal.material_at(SURFACE).unwrap().clone();
            for missing in [true, false] {
                let mut offered = surfaces.clone();
                if missing {
                    offered.pop();
                } else {
                    for patch in &mut offered {
                        patch.emissivity = 0.8;
                        patch.surroundings_k = 2000.0;
                    }
                }
                // A later thermal endpoint puts the interface above the sigma
                // claim's 401 K limit while the bulk remains inside that limit.
                assert!(
                    thermal
                        .step_with_finite_gas_surface_exchange(
                            &[],
                            &[],
                            surface_controls(),
                            &mut [],
                            &[],
                            &[],
                            &mut gas,
                            &offered,
                            gas_controls,
                            FreeSurfaceInterfaceSolve::default(),
                            false,
                            gate,
                        )
                        .is_err()
                );
                assert_eq!(thermal.flow().lattice().grid().state_bits(), before_bits);
                assert_eq!(thermal.total_enthalpy_j().unwrap(), before_energy);
                assert_eq!(thermal.material_at(SURFACE).unwrap(), &before_source);
                assert_eq!(gas.state(), &before_gas);
            }
            // Failed mode entry did not latch its persistent boundary requirement.
            thermal.step(gate).unwrap();
            thermal
                .step_with_finite_gas_surface_exchange(
                    &[],
                    &[],
                    surface_controls(),
                    &mut [],
                    &[],
                    &[],
                    &mut gas,
                    &surfaces,
                    gas_controls,
                    FreeSurfaceInterfaceSolve::default(),
                    false,
                    gate,
                )
                .unwrap();
            let accepted = thermal.flow().lattice().grid().state_bits();
            let accepted_energy = thermal.total_enthalpy_j().unwrap();
            let accepted_gas = gas.state().clone();
            assert!(matches!(
                thermal.step(gate),
                Err(FreeSurfaceEnthalpyError::Invalid(
                    "surface-temperature stress requires its complete solved heat boundary"
                ))
            ));
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(
                thermal
                    .step_with_finite_gas_surface_exchange(
                        &[],
                        &[],
                        surface_controls(),
                        &mut [],
                        &[],
                        &[],
                        &mut gas,
                        &surfaces,
                        gas_controls,
                        FreeSurfaceInterfaceSolve::default(),
                        false,
                        &cancelled,
                    )
                    .is_err()
            );
            assert_eq!(thermal.flow().lattice().grid().state_bits(), accepted);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), accepted_energy);
            assert_eq!(gas.state(), &accepted_gas);
            let next = thermal
                .step_with_finite_gas_surface_exchange(
                    &[],
                    &[],
                    surface_controls(),
                    &mut [],
                    &[],
                    &[],
                    &mut gas,
                    &surfaces,
                    gas_controls,
                    FreeSurfaceInterfaceSolve::default(),
                    false,
                    gate,
                )
                .unwrap();
            assert_eq!(next.conduction.surface.unwrap().native_endpoint.len(), 64);
        },
    );
}
