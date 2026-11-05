//! Selected gas-surface physics consumes the same molecularly changed gas.
use super::*;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::NetworkError;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::exchange::gas_law::GasSurface;
use fs_couple::cavity_phs::piston::thermal::{ConductanceSample, ThermalPistonError};
use fs_material::state_point::ResolvedMaterialStatePoint;
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
fn g3_selected_conduction_matches_molecular_reaction_and_complete_feedback_exactly() {
    fixture([0.2, 0.8], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut original_molecular) = molecular_fixture(MASS, 1.0, -50000.0);
            let (_, mut molecular) = molecular_fixture(MASS, 1.0, -50000.0);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.7, 390.0, 1e5), (0.3, 410.0, 1e5)], false),
            );
            let loads = drives(state.chambers());
            let law = reaction(50.0);
            let reactions = [ChamberReaction {
                chamber: 10,
                reaction: &law,
                controls: chemical(),
            }];
            let patches = gas_fixture::patches();
            let geometry = surfaces(&patches);
            let (expected, original) = state
                .molecular_transport_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &patches,
                    &mut original_molecular,
                    &[molecular_link()],
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let (preview, report) = state
                .molecular_transport_with_gas_law_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &geometry,
                    &mut molecular,
                    &[molecular_link()],
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    |contact, thermal, material, gas| {
                        let patch = *patches
                            .iter()
                            .find(|patch| patch.contact == contact)
                            .unwrap();
                        prepared.resolve_conductive_gas_contact(patch, thermal, material, gas)
                    },
                )
                .unwrap();
            same_owners(&preview, &expected);
            assert_eq!(report.flow, original.flow);
            assert_eq!(report.exchange.step, original.exchange.step);
            assert_eq!(
                report.exchange.report.source.reactions,
                original.exchange.report.source.reactions
            );
            assert!(
                report.exchange.report.source.reactions[0]
                    .reaction
                    .extent_mol
                    > 0.0
            );
            assert_eq!(
                report.total_mass_residual_kg,
                original.total_mass_residual_kg
            );
            assert_eq!(
                report.total_energy_residual_j,
                original.total_energy_residual_j
            );
            for (selected, original) in report
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
            assert_eq!(state.time_s(), 0.0);
            let accepted = state
                .advance_with_molecular_transport_and_gas_law(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &geometry,
                    &mut molecular,
                    &[molecular_link()],
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    |contact, thermal, material, gas| {
                        let patch = *patches
                            .iter()
                            .find(|patch| patch.contact == contact)
                            .unwrap();
                        prepared.resolve_conductive_gas_contact(patch, thermal, material, gas)
                    },
                )
                .unwrap();
            same_owners(&state, &preview);
            assert_eq!(accepted.flow, report.flow);
            assert_eq!(accepted.exchange.step, report.exchange.step);
            for chamber in state.chambers().values() {
                assert!(matches!(
                    chamber.audit_last_window().unwrap(),
                    ModuleWindowAudit::Unavailable { .. }
                ));
            }
        })
    });
}

#[derive(Clone, Debug, PartialEq)]
struct HeatReceipt {
    material: ResolvedMaterialStatePoint,
    gas: ResolvedThermochemicalGasState,
    base_conductance_w_k: f64,
    cell_temperature_k: f64,
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
fn g1_changed_composition_and_actual_midpoint_set_conductance_and_heat_once() {
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
            let parcels = state
                .chambers()
                .iter()
                .map(|(&id, chamber)| (id, chamber.gas().clone()))
                .collect();
            let endpoints = loads
                .iter()
                .map(|drive| (drive.chamber, drive.piston_solve.gas))
                .collect();
            let native = molecular
                .step_candidate(
                    &parcels,
                    &[molecular_link()],
                    GasTransportPaths::default(),
                    DT,
                    &endpoints,
                    molecular_controls(),
                    || false,
                )
                .unwrap();
            let calls = std::cell::Cell::new(0usize);
            let (next, report) = state
                .molecular_transport_with_gas_law_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &geometry,
                    &mut molecular,
                    &[molecular_link()],
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    |contact, thermal, material, gas| {
                        calls.set(calls.get() + 1);
                        let patch = *patches
                            .iter()
                            .find(|patch| patch.contact == contact)
                            .unwrap();
                        let base = prepared
                            .resolve_conductive_gas_contact(patch, thermal, material, gas)?;
                        // Explicit synthetic response: composition and actual gas
                        // temperature multiply the original material conductance.
                        let conductance =
                            base.conductance_w_k * (1.0 + mole_a(gas)) * gas.state().temperature
                                / 400.0;
                        Ok::<_, MaterialNetworkError>(ConductanceSample {
                            conductance_w_k: conductance,
                            source: HeatReceipt {
                                material: base.source,
                                gas: gas.clone(),
                                base_conductance_w_k: base.conductance_w_k,
                                cell_temperature_k: thermal.temperature_k(contact.cell).unwrap(),
                            },
                        })
                    },
                )
                .unwrap();
            assert!(calls.get() >= geometry.len());
            assert_eq!(&report.flow, native.step());
            assert!(report.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
            for chamber in &report.exchange.report.source.chambers {
                let original = state.chambers()[&chamber.id].gas();
                let staged = &native.parcels()[&chamber.id];
                let final_gas = next.chambers()[&chamber.id].gas();
                let mut heat = 0.0;
                for surface in &chamber.surfaces {
                    let (contact, receipt) = &surface.source;
                    assert_eq!(&receipt.gas, &chamber.boundary.piston.midpoint);
                    assert_eq!(mole_a(&receipt.gas), mole_a(staged.state()));
                    assert!((mole_a(&receipt.gas) - mole_a(original.state())).abs() > 1e-6);
                    assert!(
                        (receipt.gas.state().temperature - original.state().state().temperature)
                            .abs()
                            > 1e-4
                    );
                    let conductance = receipt.base_conductance_w_k
                        * (1.0 + mole_a(&receipt.gas))
                        * receipt.gas.state().temperature
                        / 400.0;
                    close(surface.conductance_w_k, conductance, 1e-14);
                    let frozen = receipt.base_conductance_w_k
                        * (1.0 + mole_a(original.state()))
                        * original.state().state().temperature
                        / 400.0;
                    assert!((surface.conductance_w_k - frozen).abs() > 1e-8);
                    close(
                        receipt.cell_temperature_k,
                        surface.port.temperature_k,
                        1e-12,
                    );
                    close(
                        surface.heat_j,
                        DT * conductance
                            * (surface.port.temperature_k - receipt.gas.state().temperature),
                        1e-12,
                    );
                    let material = report
                        .exchange
                        .report
                        .source
                        .material
                        .material
                        .bulk
                        .cells
                        .iter()
                        .find(|cell| cell.cell == contact.cell)
                        .unwrap();
                    assert_eq!(Some(&receipt.material), material.conductivity.as_ref());
                    heat += surface.heat_j;
                }
                // Fixed support and constant original Cv give an independent
                // endpoint check: molecular heat is already in the staged gas;
                // only the selected surface heat is added by this feedback leg.
                let capacity = staged.mass_kg() * 2.5 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MASS;
                close(
                    final_gas.state().state().temperature,
                    staged.state().state().temperature + heat / capacity,
                    2e-8,
                );
                close(
                    species_mass(final_gas, "A"),
                    species_mass(staged, "A"),
                    1e-16,
                );
            }
            close(
                molecular.face().law().last_resolved().temperature().value(),
                next.chambers()[&20].gas().state().state().temperature,
                2e-8,
            );
            assert_eq!(state.time_s(), 0.0);
        })
    });
}

#[derive(Debug, PartialEq, Eq)]
enum HeatRefusal {
    OutOfDomain(u64),
}

#[test]
fn g4_typed_molecular_heat_source_and_late_cancellation_preserve_existing_histories() {
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
            let patches = gas_fixture::patches();
            let geometry = surfaces(&patches);
            state
                .advance_with_molecular_transport_and_gas_law(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &geometry,
                    &mut molecular,
                    &[molecular_link()],
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    |contact, thermal, material, gas| {
                        let patch = *patches
                            .iter()
                            .find(|patch| patch.contact == contact)
                            .unwrap();
                        prepared.resolve_conductive_gas_contact(patch, thermal, material, gas)
                    },
                )
                .unwrap();
            let before = MaterialExchangeState::from_parts(
                state.thermal().clone(),
                state.circuit().clone(),
                state.chambers().clone(),
            );
            let calls = std::cell::Cell::new(0usize);
            let (_, mut wrong) = molecular_fixture(2.0 * MASS, 1.0, 0.0);
            let rejected = state.advance_with_molecular_transport_and_gas_law(
                &prepared,
                cx,
                2.0 * DT,
                &[],
                &loads,
                &[],
                &geometry,
                &mut wrong,
                &[molecular_link()],
                GasTransportPaths::default(),
                molecular_controls(),
                controls(),
                &[0.0; 2],
                |contact, _, _, _| {
                    calls.set(calls.get() + 1);
                    Err::<ConductanceSample<()>, _>(HeatRefusal::OutOfDomain(contact.id))
                },
            );
            assert!(matches!(
                rejected,
                Err(MaterialGasLawFlowError::Molecular(
                    ParcelMolecularError::SourceMismatch(_)
                ))
            ));
            assert_eq!(calls.get(), 0);
            same_owners(&state, &before);
            let failure = state.advance_with_molecular_transport_and_gas_law(
                &prepared,
                cx,
                2.0 * DT,
                &[],
                &loads,
                &[],
                &geometry,
                &mut molecular,
                &[molecular_link()],
                GasTransportPaths::default(),
                molecular_controls(),
                controls(),
                &[0.0; 2],
                |contact, _, _, _| {
                    calls.set(calls.get() + 1);
                    Err::<ConductanceSample<()>, _>(HeatRefusal::OutOfDomain(contact.id))
                },
            );
            assert!(matches!(
                failure,
                Err(MaterialGasLawFlowError::Exchange(
                    ReactiveExchangeError::Exchange(CircuitEnthalpyError::Source(
                        NetworkError::Boundary {
                            error: ThermalPistonError::Interface(HeatRefusal::OutOfDomain(101)),
                            ..
                        }
                    ))
                ))
            ));
            assert!(calls.get() > 0);
            same_owners(&state, &before);
            let cancelled = state.advance_with_molecular_transport_and_gas_law(
                &prepared,
                cx,
                2.0 * DT,
                &[],
                &loads,
                &[],
                &geometry,
                &mut molecular,
                &[molecular_link()],
                GasTransportPaths::default(),
                molecular_controls(),
                controls(),
                &[0.0; 2],
                |_, _, _, _| {
                    gate.request();
                    Ok::<_, HeatRefusal>(ConductanceSample {
                        conductance_w_k: 0.0,
                        source: (),
                    })
                },
            );
            assert!(matches!(cancelled, Err(MaterialGasLawFlowError::Cancelled)));
            same_owners(&state, &before);
        })
    });
}
