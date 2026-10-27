//! Original electrical/contact/radiative sources reach same-endpoint chemistry.

use super::coupled_core::{binding, implicit_reaction};
use super::reactive_tests::piston;
use super::*;
use fs_couple::PortedModule;
use fs_couple::cavity_phs::piston::thermal::ConductanceSample;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::exchange::gas_law::GasSurface;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::reactive::CoupledChamberStep;

fn surfaces(patches: &[MaterialContact]) -> Vec<GasSurface> {
    patches
        .iter()
        .map(|patch| GasSurface {
            contact: patch.contact,
            area_m2: patch.geometry.area_m2,
        })
        .collect()
}

#[test]
fn material_radiation_changes_same_step_chemistry_circuit_and_continuing_heat() {
    let run = |emissivity| {
        fixture(
            [emissivity; 2],
            150.0,
            |assembly, mut thermal, mut circuit| {
                with_cx(|cx, _| {
                    let prepared = assembly.prepare(cx, &thermal).unwrap();
                    let law = implicit_reaction(10000.0);
                    let mut gases = [piston(&law), gas_fixture::pistons()[1].clone()];
                    let original = gases.clone();
                    let initial_energy = thermal.energy_j() + 50.0;
                    let mut patches = gas_fixture::patches();
                    for patch in &mut patches {
                        // Explicit resolved-layer geometry makes the same-interval
                        // thermal feedback observable above the physical tolerances.
                        patch.geometry.area_m2 = 0.3;
                        patch.geometry.gas_distance_m = 1e-5;
                        patch.geometry.cell_distance_m = 1e-4;
                    }
                    let mut exterior_work = 0.0;
                    let mut extent = 0.0;
                    for _ in 0..2 {
                        let target = thermal.time_s() + 0.00005;
                        let dt = target - thermal.time_s();
                        let mut chambers = [
                            gas_fixture::chamber(20, &gases[1], dt),
                            gas_fixture::chamber(10, &gases[0], dt),
                        ];
                        for chamber in &mut chambers {
                            chamber.exchange.heat_bracket_j = [-10.0, 10.0];
                        }
                        let candidate = prepared
                            .coupled_reactive_step_candidate(
                                cx,
                                &thermal,
                                target,
                                &circuit,
                                &[],
                                &chambers,
                                &[binding(&law)],
                                &patches,
                                controls(),
                                &[0.0; 2],
                            )
                            .unwrap();
                        let (next_thermal, next_circuit, report, account) = candidate.into_parts();
                        let source = &report.source.material.material;
                        assert_eq!(source.contacts.len(), 1);
                        assert!(
                            source.contacts[0].conductance_w_k
                                < source.contacts[0].bulk_conductance_w_k
                        );
                        let radiation = source.radiation.as_ref().unwrap();
                        assert!(radiation.net_heat_w().abs() < 1e-8);
                        for surface in radiation.surfaces() {
                            let cell = source
                                .bulk
                                .cells
                                .iter()
                                .find(|c| c.cell == surface.cell)
                                .unwrap();
                            assert_eq!(surface.phase, cell.phase);
                        }
                        assert_eq!(report.electrical.contacts.len(), 1);
                        assert!(report.electrical.contacts[0].heat_j > 0.0);
                        for heat in &report.electrical.resistors {
                            let (cell_id, base) = if heat.material.branch == 1 {
                                (CellId(30), 1.0)
                            } else {
                                assert_eq!(heat.material.branch, 2);
                                (CellId(10), 2.0)
                            };
                            let cell = source
                                .bulk
                                .cells
                                .iter()
                                .find(|c| c.cell == cell_id)
                                .unwrap();
                            let material = cell.conductivity.as_ref().unwrap();
                            assert_eq!(heat.material.conductor.material(), material);
                            assert_eq!(
                                material.card_identity(),
                                cell.phase.material_card_identity()
                            );
                            assert_eq!(
                                cell.phase.temperature_k(),
                                next_thermal.temperature_k(cell_id).unwrap()
                            );
                            // R=rho(T)*L/A, with the fixture's original 1 m / 1 m^2.
                            near(
                                heat.material.conductor.resistance_ohm(),
                                base * (1.0 + (cell.phase.temperature_k() - 300.0) / 200.0),
                            );
                            near(
                                heat.heat_j,
                                dt * heat.current_a.powi(2)
                                    * heat.material.conductor.resistance_ohm(),
                            );
                        }
                        let mut actual_transfer = 0.0;
                        for chamber in &report.source.chambers {
                            let old_index = if chamber.id == 10 { 0 } else { 1 };
                            let endpoint = chamber.boundary.piston();
                            assert_eq!(endpoint.last_window_duration_s(), Some(dt));
                            assert_eq!(endpoint.gas().mass_kg(), gases[old_index].gas().mass_kg());
                            let evaluated_gas = if let Some(step) = chamber.boundary.reactive_step()
                            {
                                let chemistry = &step.piston.reaction.reaction;
                                extent += chemistry.extent_mol;
                                let residual = chemistry.extent_mol
                                    - dt * endpoint.gas().volume_m3()
                                        * chemistry.rate.net_progress_rate().unwrap().value();
                                assert!(
                                    residual.abs()
                                        <= law.implicit_mechanism().unwrap().extent_tolerance_mol
                                );
                                endpoint.gas().state()
                            } else {
                                &chamber.boundary.nonreactive_step().unwrap().piston.midpoint
                            };
                            for surface in chamber.boundary.source() {
                                let cell = source
                                    .bulk
                                    .cells
                                    .iter()
                                    .find(|c| c.cell == surface.source.0.cell)
                                    .unwrap();
                                assert_eq!(cell.conductivity.as_ref(), Some(&surface.source.1));
                                let k = surface
                                    .source
                                    .1
                                    .property("thermal_conductivity")
                                    .unwrap()
                                    .value_si();
                                let expected_g = 0.3
                                    / (1e-5 / evaluated_gas.state().thermal_conductivity
                                        + 1e-4 / k);
                                near(surface.conductance_w_k, expected_g);
                                near(
                                    surface.heat_j,
                                    dt * expected_g
                                        * (cell.phase.temperature_k()
                                            - evaluated_gas.state().temperature),
                                );
                                actual_transfer += surface.heat_j;
                            }
                        }
                        near(account.transferred_heat_j, actual_transfer);
                        near(account.external_network_heat_j, 0.0);
                        assert!(
                            account.total_energy_residual_j.abs()
                                <= controls().total_energy_tolerance_j
                        );
                        for (index, chamber) in report.source.chambers.into_iter().enumerate() {
                            assert_eq!(chamber.id, if index == 0 { 10 } else { 20 });
                            let (gas, _, step) = chamber.boundary.into_parts();
                            exterior_work += match step {
                                CoupledChamberStep::Reactive(step) => {
                                    step.piston.exterior_pressure_work_j
                                }
                                CoupledChamberStep::Nonreactive(step) => {
                                    step.piston.exterior_pressure_work_j
                                }
                            };
                            gases[index] = gas;
                        }
                        let gas_change: f64 = gases
                            .iter()
                            .zip(&original)
                            .map(|(gas, old)| {
                                gas.gas().mass_kg()
                                    * (gas.gas().state().specific_internal_energy_j_kg()
                                        - old.gas().state().specific_internal_energy_j_kg())
                                    + 0.5 * gas.mass_kg() * gas.velocity_m_s().powi(2)
                            })
                            .sum();
                        assert!(
                            (next_thermal.energy_j()
                                + report.electrical.stored_after_j
                                + gas_change
                                + exterior_work
                                - initial_energy)
                                .abs()
                                < 3e-6
                        );
                        assert_ne!(
                            next_circuit.resistors()[0].conductor.resistance_ohm(),
                            circuit.resistors()[0].conductor.resistance_ohm()
                        );
                        thermal = next_thermal;
                        circuit = next_circuit;
                    }
                    (
                        thermal.temperature_k(CellId(30)).unwrap(),
                        thermal.temperature_k(CellId(10)).unwrap(),
                        circuit.charges_c()[0],
                        extent,
                        gases[0].gas().state().state().temperature,
                    )
                })
            },
        )
    };
    let low = run(0.2);
    let high = run(0.8);
    assert!((low.0 - high.0).abs() > 1e-4);
    assert!((low.1 - high.1).abs() > 1e-4);
    assert!((low.2 - high.2).abs() > 1e-10);
    assert!(low.3 > high.3);
    assert!((low.4 - high.4).abs() > 1e-6);
}

#[test]
fn custom_heat_law_retains_exact_material_and_reacted_endpoint_sources() {
    fixture([0.2, 0.8], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let law = implicit_reaction(1000.0);
            let gases = [piston(&law), gas_fixture::pistons()[1].clone()];
            let chambers = [
                gas_fixture::chamber(10, &gases[0], 0.00005),
                gas_fixture::chamber(20, &gases[1], 0.00005),
            ];
            let patches = gas_fixture::patches();
            let geometry = surfaces(&patches);
            let expected = prepared
                .coupled_reactive_step_candidate(
                    cx,
                    &thermal,
                    0.00005,
                    &circuit,
                    &[],
                    &chambers,
                    &[binding(&law)],
                    &patches,
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let result = prepared
                .coupled_reactive_with_gas_law_candidate(
                    cx,
                    &thermal,
                    0.00005,
                    &circuit,
                    &[],
                    &chambers,
                    &[binding(&law)],
                    &geometry,
                    controls(),
                    &[0.0; 2],
                    |contact, state, source, gas| {
                        let patch = *patches.iter().find(|p| p.contact == contact).unwrap();
                        let sample =
                            prepared.resolve_conductive_gas_contact(patch, state, source, gas)?;
                        let phase = source
                            .bulk
                            .cells
                            .iter()
                            .find(|c| c.cell == contact.cell)
                            .unwrap()
                            .phase;
                        Ok::<_, MaterialNetworkError>(ConductanceSample {
                            conductance_w_k: sample.conductance_w_k,
                            source: (sample.source, gas.clone(), phase),
                        })
                    },
                )
                .unwrap();
            assert_eq!(result.thermal().cells(), expected.thermal().cells());
            assert_eq!(result.circuit(), expected.circuit());
            assert_eq!(result.step(), expected.step());
            for (actual, expected) in result.chambers().iter().zip(expected.chambers()) {
                assert_eq!(actual.boundary.piston(), expected.boundary.piston());
                assert_eq!(
                    actual.boundary.reactive_step(),
                    expected.boundary.reactive_step()
                );
                assert_eq!(
                    actual.boundary.nonreactive_step(),
                    expected.boundary.nonreactive_step()
                );
                let evaluated = match actual.boundary.reactive_step() {
                    Some(_) => actual.boundary.piston().gas().state(),
                    None => &actual.boundary.nonreactive_step().unwrap().piston.midpoint,
                };
                for surface in actual.boundary.source() {
                    let (contact, (material, gas, phase)) = &surface.source;
                    assert_eq!(gas, evaluated);
                    let cell = result
                        .material()
                        .material
                        .bulk
                        .cells
                        .iter()
                        .find(|c| c.cell == contact.cell)
                        .unwrap();
                    assert_eq!(*phase, cell.phase);
                    assert_eq!(Some(material), cell.conductivity.as_ref());
                }
            }
        })
    });
}

#[test]
fn invalid_mappings_zero_area_and_late_source_cancellation_publish_nothing() {
    fixture([0.2; 2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, gate| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let law = implicit_reaction(1000.0);
            let gas = piston(&law);
            let original_gas = gas.clone();
            let original_cells = thermal.cells().clone();
            let original_circuit = circuit.clone();
            let chambers = [gas_fixture::chamber(10, &gas, 0.00005)];
            let reaction = binding(&law);
            for reactions in [vec![reaction, reaction], vec![
            fs_couple::cavity_phs::piston::thermal::boundary::ports::network::reactive::CoupledChamberReaction {
                chamber: 99, ..reaction
            }]] {
            assert!(prepared.coupled_reactive_step_candidate(cx, &thermal, 0.00005,
                &circuit, &[], &chambers, &reactions, &[], controls(), &[0.0; 2]).is_err());
        }
            let mut geometry = surfaces(&gas_fixture::patches()[..1]);
            geometry[0].area_m2 = 0.0;
            let impossible = prepared.coupled_reactive_with_gas_law_candidate(
                cx,
                &thermal,
                0.00005,
                &circuit,
                &[],
                &chambers,
                &[reaction],
                &geometry,
                controls(),
                &[0.0; 2],
                |_, _, _, _| {
                    Ok::<_, &'static str>(ConductanceSample {
                        conductance_w_k: 0.2,
                        source: (),
                    })
                },
            );
            assert!(impossible.is_err());
            geometry[0].area_m2 = 0.1;
            let cancelled = prepared.coupled_reactive_with_gas_law_candidate(
                cx,
                &thermal,
                0.00005,
                &circuit,
                &[],
                &chambers,
                &[reaction],
                &geometry,
                controls(),
                &[0.0; 2],
                |_, _, _, gas| {
                    assert_eq!(gas.components().len(), 2);
                    gate.request();
                    Ok::<_, &'static str>(ConductanceSample {
                        conductance_w_k: 0.0,
                        source: (),
                    })
                },
            );
            assert!(matches!(cancelled, Err(CircuitEnthalpyError::Cancelled)));
            assert_eq!(gas, original_gas);
            assert_eq!(
                gas.last_window_duration_s(),
                original_gas.last_window_duration_s()
            );
            assert_eq!(thermal.cells(), &original_cells);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(circuit, original_circuit);
        })
    });
}
