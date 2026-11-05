//! Complete native path sets retain adaptive selected heat laws and their sources.
use super::*;
use fs_couple::cavity_phs::piston::thermal::ThermalPistonError;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::{Contact, NetworkError};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

fn surfaces(patches: &[MaterialContact]) -> Vec<GasSurface> {
    patches
        .iter()
        .map(|patch| GasSurface {
            contact: patch.contact,
            area_m2: patch.geometry.area_m2,
        })
        .collect()
}

fn same_owners(actual: &MaterialExchangeState<'_>, expected: &MaterialExchangeState<'_>) {
    assert_eq!(actual.time_s(), expected.time_s());
    assert_eq!(actual.chambers(), expected.chambers());
    assert_eq!(actual.circuit(), expected.circuit());
    assert_eq!(actual.thermal().cells(), expected.thermal().cells());
    for (&id, chamber) in actual.chambers() {
        assert_eq!(
            chamber.last_window_transfers(),
            expected.chambers()[&id].last_window_transfers()
        );
        assert_eq!(
            chamber.last_window_duration_s(),
            expected.chambers()[&id].last_window_duration_s()
        );
    }
}

#[test]
fn selected_conduction_reproduces_the_mixed_laminate_adaptive_trajectory() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let card = pair_card(2e-5);
            let point = QueryPoint::new();
            let pairs = binding(&card, &point);
            let medium = laminate(1.0);
            let layers = [laminate_link(&medium, &pairs, &point)];
            let nozzles = [restriction()];
            let paths = GasTransportPaths {
                nozzles: &nozzles,
                layered_dusty_gas: &layers,
                ..Default::default()
            };
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(thermal, circuit, chambers());
            let mut loads = drives(state.chambers());
            loads[0].heat_j = 0.01;
            loads[1].heat_j = 0.02;
            let patches = gas_fixture::patches();
            let geometry = surfaces(&patches);
            let time = CoupledTransportControls {
                maximum_evaluations: 3,
                ..loose_time()
            };
            let (expected, original) = state
                .coupled_transport_paths_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &patches,
                    paths,
                    admission(),
                    coupled(),
                    &[0.0; 2],
                    time,
                )
                .unwrap();
            let resolve = |contact: Contact,
                           thermal: &DepositionSimulation<'_>,
                           material: &MaterialExchangeSample,
                           gas: &ResolvedThermochemicalGasState| {
                let patch = *patches
                    .iter()
                    .find(|patch| patch.contact == contact)
                    .unwrap();
                prepared.resolve_conductive_gas_contact(patch, thermal, material, gas)
            };
            let (preview, selected) = state
                .coupled_transport_paths_with_gas_law_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &geometry,
                    paths,
                    admission(),
                    coupled(),
                    &[0.0; 2],
                    time,
                    &resolve,
                )
                .unwrap();
            same_owners(&preview, &expected);
            assert_eq!(selected.balance, original.balance);
            assert_eq!(selected.estimated_error, original.estimated_error);
            assert_eq!(selected.evaluations, 3);
            assert_eq!(selected.accepted_intervals, 1);
            assert_eq!(selected.rejected_intervals, 0);
            assert_eq!(selected.steps.len(), 2);
            for (selected, original) in selected.steps.iter().zip(&original.steps) {
                assert_eq!(selected.start_s, original.start_s);
                assert_eq!(selected.end_s, original.end_s);
                assert_eq!(selected.advance.flow, original.advance.flow);
                assert!(!selected.advance.flow.layered_dusty_gas.is_empty());
                assert!(!selected.advance.flow.nozzles.is_empty());
                assert_eq!(
                    selected.advance.exchange.step,
                    original.advance.exchange.step
                );
                assert_eq!(
                    selected.advance.total_mass_residual_kg,
                    original.advance.total_mass_residual_kg
                );
                assert_eq!(
                    selected.advance.total_energy_residual_j,
                    original.advance.total_energy_residual_j
                );
                for (selected, original) in selected
                    .advance
                    .exchange
                    .report
                    .source
                    .chambers
                    .iter()
                    .zip(&original.advance.exchange.report.source.chambers)
                {
                    assert_eq!(selected.id, original.id);
                    assert_eq!(selected.boundary, original.boundary);
                    assert_eq!(selected.surfaces, original.surfaces);
                }
            }
            assert_eq!(state.time_s(), 0.0);
            let accepted = state
                .advance_with_coupled_transport_paths_and_gas_law(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &geometry,
                    paths,
                    admission(),
                    coupled(),
                    &[0.0; 2],
                    time,
                    resolve,
                )
                .unwrap();
            same_owners(&state, &preview);
            assert_eq!(accepted.balance, selected.balance);
            for chamber in state.chambers().values() {
                assert_eq!(chamber.last_window_duration_s(), Some(DT));
                assert!(matches!(
                    chamber.audit_last_window().unwrap(),
                    ModuleWindowAudit::Unavailable { .. }
                ));
            }
        });
    });
}

// Neither retained sources nor source refusals need Clone or Debug.
struct HeatReceipt {
    gas: ResolvedThermochemicalGasState,
    time_s: f64,
}
struct HeatRefusal {
    contact: u64,
}

fn mole_a(gas: &ResolvedThermochemicalGasState) -> f64 {
    let receipt = gas.thermodynamics().receipt();
    receipt
        .species()
        .iter()
        .position(|id| id == &SpeciesId::new("A").unwrap())
        .map_or(0.0, |index| {
            f64::from_bits(receipt.mole_fraction_bits()[index])
        })
}

#[test]
fn changed_mixtures_drive_owned_heat_sources_and_late_refusal_preserves_all_owners() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let card = pair_card(2e-5);
            let point = QueryPoint::new();
            let pairs = binding(&card, &point);
            let medium = laminate(1.0);
            let layers = [laminate_link(&medium, &pairs, &point)];
            let nozzles = [restriction()];
            let paths = GasTransportPaths {
                nozzles: &nozzles,
                layered_dusty_gas: &layers,
                ..Default::default()
            };
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(thermal, circuit, chambers());
            let loads = drives(state.chambers());
            let geometry = surfaces(&gas_fixture::patches());
            let resolve = |_: Contact,
                           thermal: &DepositionSimulation<'_>,
                           _: &MaterialExchangeSample,
                           gas: &ResolvedThermochemicalGasState| {
                Ok::<_, HeatRefusal>(ConductanceSample {
                    conductance_w_k: (1.0 + mole_a(gas)) * gas.state().temperature / 400.0,
                    source: HeatReceipt {
                        gas: gas.clone(),
                        time_s: thermal.time_s(),
                    },
                })
            };
            let report = state
                .advance_with_coupled_transport_paths_and_gas_law(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &geometry,
                    paths,
                    admission(),
                    coupled(),
                    &[0.0; 2],
                    CoupledTransportControls {
                        maximum_evaluations: 3,
                        ..loose_time()
                    },
                    &resolve,
                )
                .unwrap_or_else(|_| panic!("composition-dependent selected heat law must advance"));
            assert_eq!(report.steps.len(), 2);
            let mut previous = BTreeMap::new();
            for step in &report.steps {
                for chamber in &step.advance.exchange.report.source.chambers {
                    for surface in &chamber.surfaces {
                        let receipt = &surface.source.1;
                        assert_eq!(receipt.time_s, step.end_s);
                        assert_eq!(receipt.gas, chamber.boundary.piston.midpoint);
                        let fraction = mole_a(&receipt.gas);
                        // Both initial contacted chambers were pure: all selected
                        // coefficients must consume their actual transported mixture.
                        assert!(fraction > 0.0 && fraction < 1.0);
                        if let Some(old) = previous.insert(chamber.id, fraction) {
                            assert_ne!(fraction, old);
                        }
                        let conductance =
                            (1.0 + fraction) * receipt.gas.state().temperature / 400.0;
                        assert_eq!(surface.conductance_w_k, conductance);
                        let expected_heat = (step.end_s - step.start_s)
                            * conductance
                            * (surface.port.temperature_k - receipt.gas.state().temperature);
                        assert!((surface.heat_j - expected_heat).abs() < 1e-12);
                    }
                }
            }
            assert_eq!(previous.len(), 2);
            assert!(
                report.balance.total_energy_residual_j.abs() <= coupled().total_energy_tolerance_j
            );
            assert!(
                report.balance.maximum_species_mass_residual_kg
                    <= admission().exchange.mass_tolerance_kg
            );
            let before = MaterialExchangeState::from_parts(
                state.thermal().clone(),
                state.circuit().clone(),
                state.chambers().clone(),
            );
            let allowed = std::cell::Cell::new(false);
            let refused = std::cell::Cell::new(0usize);
            let failure = state.advance_with_coupled_transport_paths_and_gas_law(
                &prepared,
                cx,
                3.0 * DT,
                &[],
                &loads,
                &[],
                &geometry,
                paths,
                admission(),
                coupled(),
                &[0.0; 2],
                loose_time(),
                |contact, thermal, material, gas| {
                    if thermal.time_s() > 2.0 * DT {
                        refused.set(refused.get() + 1);
                        return Err(HeatRefusal {
                            contact: contact.id,
                        });
                    }
                    allowed.set(true);
                    resolve(contact, thermal, material, gas)
                },
            );
            match failure {
                Err(CoupledTransportGasLawError::Step {
                    start_s,
                    end_s,
                    error:
                        MaterialGasLawFlowError::Exchange(ReactiveExchangeError::Exchange(
                            CircuitEnthalpyError::Source(NetworkError::Boundary {
                                error: ThermalPistonError::Interface(HeatRefusal { contact: 101 }),
                                ..
                            }),
                        )),
                }) => {
                    // The first comparison accepted its two fine PRIVATE steps;
                    // the later original typed source refusal must remain fatal.
                    assert_eq!(start_s, 2.0 * DT);
                    assert_eq!(end_s, 3.0 * DT);
                }
                _ => panic!("expected the original late typed heat-law refusal"),
            }
            assert!(allowed.get());
            assert_eq!(refused.get(), 1);
            same_owners(&state, &before);
        });
    });
}
