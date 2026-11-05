//! G1/G3/G4: selected heat sources follow every accepted molecular fine step.
use super::*;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::{Contact, NetworkError};
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::exchange::gas_law::GasSurface;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state::flow::molecular::coupled::gas_law::*;
use fs_couple::cavity_phs::piston::thermal::{ConductanceSample, ThermalPistonError};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkError, transport::GasTransportError,
};

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
fn g3_selected_conduction_reproduces_the_complete_reacting_adaptive_trajectory() {
    fixture([0.2, 0.8], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut original) = molecular_fixture(MASS, 1.0, -50000.0);
            let (_, mut molecular) = molecular_fixture(MASS, 1.0, -50000.0);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.7, 390.0, 1e5), (0.3, 410.0, 1e5)], false),
            );
            let mut loads = drives(state.chambers());
            loads[0].heat_j = 0.015625;
            loads[1].heat_j = 0.03125;
            loads[0].external_force_n = 2.0;
            let law = reaction(50.0);
            let reactions = [ChamberReaction {
                chamber: 10,
                reaction: &law,
                controls: chemical(),
            }];
            let patches = gas_fixture::patches();
            let geometry = surfaces(&patches);
            let links = [molecular_link()];
            let time = CoupledMolecularControls {
                maximum_evaluations: 3,
                ..loose_time()
            };
            let (expected, original) = state
                .coupled_molecular_transport_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &patches,
                    &mut original,
                    &links,
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
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
                .coupled_molecular_transport_with_gas_law_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &geometry,
                    &mut molecular,
                    &links,
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
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
                let selected = &selected.advance;
                let original = &original.advance;
                assert_eq!(selected.flow, original.flow);
                assert_eq!(selected.exchange.step, original.exchange.step);
                assert_eq!(
                    selected.exchange.report.source.reactions,
                    original.exchange.report.source.reactions
                );
                assert!(
                    selected.exchange.report.source.reactions[0]
                        .reaction
                        .extent_mol
                        > 0.0
                );
                for (selected, original) in selected
                    .exchange
                    .report
                    .source
                    .chambers
                    .iter()
                    .zip(&original.exchange.report.source.chambers)
                {
                    assert_eq!(selected.boundary, original.boundary);
                    assert_eq!(selected.surfaces, original.surfaces);
                }
            }
            assert_eq!(state.time_s(), 0.0);
            let accepted = state
                .advance_with_coupled_molecular_transport_and_gas_law(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &geometry,
                    &mut molecular,
                    &links,
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    time,
                    resolve,
                )
                .unwrap();
            same_owners(&state, &preview);
            assert_eq!(accepted.balance, selected.balance);
            for chamber in state.chambers().values() {
                assert!(chamber.velocity_m_s().abs() > 1e-8);
                assert_eq!(chamber.last_window_duration_s(), Some(DT));
                assert!(matches!(
                    chamber.audit_last_window().unwrap(),
                    ModuleWindowAudit::Unavailable { .. }
                ));
            }
        });
    });
}

// Deliberately not Clone: every accepted receipt is retained by ownership.
#[derive(Debug)]
struct HeatReceipt {
    gas: ResolvedThermochemicalGasState,
    material: ResolvedMaterialStatePoint,
    base_conductance_w_k: f64,
    cell_temperature_k: f64,
    time_s: f64,
}

fn mole_a(gas: &ResolvedThermochemicalGasState) -> f64 {
    let receipt = gas.thermodynamics().receipt();
    let index = receipt
        .species()
        .iter()
        .position(|id| id == &species_id("A"))
        .unwrap();
    f64::from_bits(receipt.mole_fraction_bits()[index])
}

#[test]
fn g1_native_withdrawal_refines_and_each_fine_step_requeries_actual_composition_and_heat() {
    fixture([0.8, 0.8], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.8, 390.0, 1e5), (0.2, 410.0, 1e5)], true),
            );
            let loads = drives(state.chambers());
            let patches = gas_fixture::patches();
            let geometry = surfaces(&patches);
            let links = [molecular_link()];
            let mut physical = molecular_controls();
            physical.physical.exchange.maximum_species_fraction = 5e-4;
            let calls = std::cell::Cell::new(0usize);
            let full_step = state.molecular_transport_with_gas_law_candidate(
                &prepared,
                cx,
                DT,
                &[],
                &loads,
                &[],
                &geometry,
                &mut molecular,
                &links,
                GasTransportPaths::default(),
                physical,
                controls(),
                &[0.0; 2],
                |_, _, _, _| {
                    calls.set(calls.get() + 1);
                    Ok::<_, HeatRefusal>(ConductanceSample {
                        conductance_w_k: 2.0,
                        source: (),
                    })
                },
            );
            assert!(matches!(
                full_step,
                Err(MaterialGasLawFlowError::Molecular(
                    ParcelMolecularError::Transport(GasTransportError::Exchange(
                        SpeciesNetworkError::Withdrawal { .. }
                    ))
                ))
            ));
            assert_eq!(calls.get(), 0);
            let (next, report) = state
                .coupled_molecular_transport_with_gas_law_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &geometry,
                    &mut molecular,
                    &links,
                    GasTransportPaths::default(),
                    physical,
                    controls(),
                    &[0.0; 2],
                    loose_time(),
                    |contact, thermal, material, gas| {
                        let patch = *patches
                            .iter()
                            .find(|patch| patch.contact == contact)
                            .unwrap();
                        let base = prepared
                            .resolve_conductive_gas_contact(patch, thermal, material, gas)?;
                        let conductance_w_k =
                            base.conductance_w_k * (1.0 + mole_a(gas)) * gas.state().temperature
                                / 400.0;
                        Ok::<_, MaterialNetworkError>(ConductanceSample {
                            conductance_w_k,
                            source: HeatReceipt {
                                gas: gas.clone(),
                                material: base.source,
                                base_conductance_w_k: base.conductance_w_k,
                                cell_temperature_k: thermal.temperature_k(contact.cell).unwrap(),
                                time_s: thermal.time_s(),
                            },
                        })
                    },
                )
                .unwrap();
            assert!(report.rejected_intervals > 0);
            assert!(report.accepted_intervals > 1);
            assert_eq!(report.steps.len(), 2 * report.accepted_intervals);
            assert!(report.estimated_error.maximum() <= 1.0);
            let mut a_mass: BTreeMap<_, _> = state
                .chambers()
                .iter()
                .map(|(&id, chamber)| (id, species_mass(chamber.gas(), "A")))
                .collect();
            let mut energy: BTreeMap<_, f64> =
                state.chambers().keys().map(|&id| (id, 0.0)).collect();
            let mut previous_composition = BTreeMap::<u64, f64>::new();
            for step in &report.steps {
                let exchange = &step.advance.flow.transport.exchange;
                for flux in &exchange.fluxes {
                    let transfer = &flux.transfer;
                    if transfer.species == species_id("A") {
                        *a_mass.get_mut(&transfer.donor).unwrap() -= transfer.mass_kg;
                        *a_mass.get_mut(&transfer.receiver).unwrap() += transfer.mass_kg;
                    }
                }
                for parcel in &exchange.parcels {
                    *energy.get_mut(&parcel.id).unwrap() +=
                        parcel.net_enthalpy_j + parcel.net_additional_heat_j;
                }
                let source = &step.advance.exchange.report.source;
                for chamber in &source.chambers {
                    let original = state.chambers()[&chamber.id].gas();
                    let expected_fraction = a_mass[&chamber.id] / original.mass_kg();
                    let actual_fraction = mole_a(&chamber.boundary.piston.midpoint);
                    close(actual_fraction, expected_fraction, 2e-12);
                    if let Some(previous) = previous_composition.insert(chamber.id, actual_fraction)
                    {
                        assert!((actual_fraction - previous).abs() > 1e-6);
                    }
                    for surface in &chamber.surfaces {
                        let (contact, receipt) = &surface.source;
                        assert_eq!(receipt.time_s, step.end_s);
                        assert_eq!(receipt.gas, chamber.boundary.piston.midpoint);
                        close(mole_a(&receipt.gas), expected_fraction, 2e-12);
                        close(
                            receipt.cell_temperature_k,
                            surface.port.temperature_k,
                            1e-12,
                        );
                        let conductance = receipt.base_conductance_w_k
                            * (1.0 + expected_fraction)
                            * receipt.gas.state().temperature
                            / 400.0;
                        close(surface.conductance_w_k, conductance, 1e-12);
                        close(
                            surface.heat_j,
                            (step.end_s - step.start_s)
                                * conductance
                                * (surface.port.temperature_k - receipt.gas.state().temperature),
                            1e-12,
                        );
                        let material = source
                            .material
                            .material
                            .bulk
                            .cells
                            .iter()
                            .find(|row| row.cell == contact.cell)
                            .unwrap();
                        assert_eq!(Some(&receipt.material), material.conductivity.as_ref());
                        *energy.get_mut(&chamber.id).unwrap() += surface.heat_j;
                    }
                }
            }
            for (&id, original) in state.chambers() {
                let after = next.chambers()[&id].gas();
                close(after.mass_kg(), original.gas().mass_kg(), 1e-15);
                close(species_mass(after, "A"), a_mass[&id], 1e-15);
                // Identical frozen Cv, zero formation energies and fixed support
                // give an independent endpoint storage check for the whole window.
                let capacity =
                    original.gas().mass_kg() * 2.5 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MASS;
                close(
                    capacity
                        * (after.state().state().temperature
                            - original.gas().state().state().temperature),
                    energy[&id],
                    2e-7,
                );
            }
            assert!(
                report.balance.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j
            );
            assert!(
                report.balance.maximum_species_mass_residual_kg
                    <= physical.physical.exchange.mass_tolerance_kg
            );
            assert!(
                report.balance.maximum_cell_residual_bound_j
                    <= controls().network.cell_energy_tolerance_j
            );
            assert_eq!(next.time_s(), DT);
            assert_eq!(state.time_s(), 0.0);
        });
    });
}

// Deliberately not Clone: a typed constitutive refusal survives the time wrapper.
#[derive(Debug, PartialEq, Eq)]
enum HeatRefusal {
    OutOfDomain(u64),
}

#[test]
fn g4_late_selected_source_refusal_and_cancellation_publish_no_accepted_prefix() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, gate| {
            let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.7, 390.0, 1e5), (0.3, 410.0, 1e5)], false),
            );
            let loads = drives(state.chambers());
            let geometry = surfaces(&gas_fixture::patches());
            let links = [molecular_link()];
            state
                .advance_with_coupled_molecular_transport_and_gas_law(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &geometry,
                    &mut molecular,
                    &links,
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    loose_time(),
                    |_, _, _, _| {
                        Ok::<_, HeatRefusal>(ConductanceSample {
                            conductance_w_k: 2.0,
                            source: (),
                        })
                    },
                )
                .unwrap();
            let before = MaterialExchangeState::from_parts(
                state.thermal().clone(),
                state.circuit().clone(),
                state.chambers().clone(),
            );
            let allowed = std::cell::Cell::new(false);
            let refused = std::cell::Cell::new(0usize);
            let failure = state.advance_with_coupled_molecular_transport_and_gas_law(
                &prepared,
                cx,
                3.0 * DT,
                &[],
                &loads,
                &[],
                &geometry,
                &mut molecular,
                &links,
                GasTransportPaths::default(),
                molecular_controls(),
                controls(),
                &[0.0; 2],
                loose_time(),
                |contact, thermal, _, _| {
                    if thermal.time_s() > 2.0 * DT {
                        refused.set(refused.get() + 1);
                        return Err(HeatRefusal::OutOfDomain(contact.id));
                    }
                    allowed.set(true);
                    Ok(ConductanceSample {
                        conductance_w_k: 2.0,
                        source: (),
                    })
                },
            );
            match failure {
                Err(CoupledMolecularGasLawError::Step {
                    start_s,
                    end_s,
                    error:
                        MaterialGasLawFlowError::Exchange(ReactiveExchangeError::Exchange(
                            CircuitEnthalpyError::Source(NetworkError::Boundary {
                                error: ThermalPistonError::Interface(HeatRefusal::OutOfDomain(101)),
                                ..
                            }),
                        )),
                }) => {
                    // The first comparison has already accepted its two fine
                    // private steps. A source refusal is fatal, not refinable.
                    assert_eq!(start_s, 2.0 * DT);
                    assert_eq!(end_s, 3.0 * DT);
                }
                other => panic!(
                    "expected the original late typed source refusal, got {}",
                    match other {
                        Ok(_) => "an accepted window".into(),
                        Err(error) => format!("{error:?}"),
                    }
                ),
            }
            assert!(allowed.get());
            assert_eq!(refused.get(), 1);
            same_owners(&state, &before);
            let reached_private_end = std::cell::Cell::new(false);
            let cancelled = state.advance_with_coupled_molecular_transport_and_gas_law(
                &prepared,
                cx,
                3.0 * DT,
                &[],
                &loads,
                &[],
                &geometry,
                &mut molecular,
                &links,
                GasTransportPaths::default(),
                molecular_controls(),
                controls(),
                &[0.0; 2],
                loose_time(),
                |_, thermal, _, _| {
                    if thermal.time_s() == 2.0 * DT {
                        reached_private_end.set(true);
                    }
                    if thermal.time_s() > 2.0 * DT {
                        gate.request();
                    }
                    Ok::<_, HeatRefusal>(ConductanceSample {
                        conductance_w_k: 2.0,
                        source: (),
                    })
                },
            );
            assert!(matches!(
                cancelled,
                Err(CoupledMolecularGasLawError::Cancelled)
            ));
            assert!(reached_private_end.get());
            same_owners(&state, &before);
        });
    });
}
