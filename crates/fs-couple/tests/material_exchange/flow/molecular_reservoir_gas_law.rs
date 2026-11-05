//! Selected heat laws see the actual gas after simultaneous bath/internal transport.
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

#[test]
fn selected_conduction_matches_held_bath_reaction_and_complete_feedback_exactly() {
    fixture([0.2, 0.8], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut original_molecular) = molecular_fixture(MASS, 1.0, -50000.0);
            let (_, mut molecular) = molecular_fixture(MASS, 1.0, -50000.0);
            let held = bath(&bank, 0.9, 430.0, 1.2e5);
            let retained = held.clone();
            let links = [bath_link(&held)];
            let internal = [molecular_link()];
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.7, 390.0, 1.1e5), (0.3, 410.0, 0.95e5)], false),
            );
            let mut loads = drives(state.chambers());
            loads[0].heat_j = 0.01;
            loads[1].external_force_n = 2.0;
            let law = reaction(50.0);
            let reactions = [ChamberReaction {
                chamber: 10,
                reaction: &law,
                controls: chemical(),
            }];
            let patches = gas_fixture::patches();
            let geometry = surfaces(&patches);
            let (expected, original) = state
                .molecular_reservoir_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &patches,
                    &mut original_molecular,
                    &internal,
                    &links,
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let (preview, report) = state
                .molecular_reservoir_with_gas_law_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &geometry,
                    &mut molecular,
                    &internal,
                    &links,
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
            assert_eq!(report.balance, original.balance);
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
                assert!(preview.chambers()[&selected.id].velocity_m_s().abs() > 1e-8);
            }
            assert_eq!(state.time_s(), 0.0);
            let accepted = state
                .advance_with_molecular_reservoir_and_gas_law(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &geometry,
                    &mut molecular,
                    &internal,
                    &links,
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
            assert_eq!(accepted.balance, report.balance);
            assert_eq!(accepted.exchange.step, report.exchange.step);
            assert_eq!(held, retained);
        })
    });
}

// Deliberately not Clone: the balance carries the original typed receipt without
// imposing an unrelated ownership or erasure requirement on a user's heat law.
#[derive(Debug, PartialEq)]
struct HeatReceipt {
    material: ResolvedMaterialStatePoint,
    gas: ResolvedThermochemicalGasState,
    base_conductance_w_k: f64,
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
fn bath_changed_composition_and_actual_midpoint_set_surface_heat_once() {
    fixture([0.8, 0.8], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let held = bath(&bank, 0.1, 370.0, 1e5);
            let links = [bath_link(&held)];
            let internal = [molecular_link()];
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
                .network_candidate(
                    &parcels,
                    &internal,
                    &links,
                    DT,
                    &endpoints,
                    molecular_controls(),
                    || false,
                )
                .unwrap();
            let (next, report) = state
                .molecular_reservoir_with_gas_law_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &geometry,
                    &mut molecular,
                    &internal,
                    &links,
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    |contact, thermal, material, gas| {
                        let patch = *patches
                            .iter()
                            .find(|patch| patch.contact == contact)
                            .unwrap();
                        let base = prepared
                            .resolve_conductive_gas_contact(patch, thermal, material, gas)?;
                        // Synthetic constitutive response to the actual trial gas.
                        Ok::<_, MaterialNetworkError>(ConductanceSample {
                            conductance_w_k: base.conductance_w_k
                                * (1.0 + mole_a(gas))
                                * gas.state().temperature
                                / 400.0,
                            source: HeatReceipt {
                                material: base.source,
                                gas: gas.clone(),
                                base_conductance_w_k: base.conductance_w_k,
                            },
                        })
                    },
                )
                .unwrap();
            assert_eq!(&report.flow, native.step());
            assert!(
                report.balance.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j
            );
            assert!(report.flow.links[0].mass_into_gas_kg[0].abs() > 1e-10);
            for chamber in &report.exchange.report.source.chambers {
                let original = state.chambers()[&chamber.id].gas();
                let staged = &native.parcels()[&chamber.id];
                let final_gas = next.chambers()[&chamber.id].gas();
                let mut heat = 0.0;
                assert!(!chamber.surfaces.is_empty());
                for surface in &chamber.surfaces {
                    let (contact, receipt) = &surface.source;
                    assert_eq!(&receipt.gas, &chamber.boundary.piston.midpoint);
                    assert_eq!(mole_a(&receipt.gas), mole_a(staged.state()));
                    assert!((mole_a(&receipt.gas) - mole_a(original.state())).abs() > 1e-6);
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
                // Native bath enthalpy and Dufour/Fourier heat are already in
                // the staged parcel. Only selected surface heat enters this
                // fixed-support, constant-Cv feedback increment.
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
                close(
                    species_mass(final_gas, "B"),
                    species_mass(staged, "B"),
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
fn typed_bath_heat_refusal_and_late_cancellation_preserve_accepted_histories() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, gate| {
            let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let held = bath(&bank, 0.9, 430.0, 1e5);
            let retained = held.clone();
            let links = [bath_link(&held)];
            let internal = [molecular_link()];
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
                .advance_with_molecular_reservoir(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &patches,
                    &mut molecular,
                    &internal,
                    &links,
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let before = MaterialExchangeState::from_parts(
                state.thermal().clone(),
                state.circuit().clone(),
                state.chambers().clone(),
            );
            let calls = std::cell::Cell::new(0usize);
            let (_, mut wrong) = molecular_fixture(2.0 * MASS, 1.0, 0.0);
            let rejected = state.advance_with_molecular_reservoir_and_gas_law(
                &prepared,
                cx,
                2.0 * DT,
                &[],
                &loads,
                &[],
                &geometry,
                &mut wrong,
                &internal,
                &links,
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
            let failure = state.advance_with_molecular_reservoir_and_gas_law(
                &prepared,
                cx,
                2.0 * DT,
                &[],
                &loads,
                &[],
                &geometry,
                &mut molecular,
                &internal,
                &links,
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
            let cancelled = state.advance_with_molecular_reservoir_and_gas_law(
                &prepared,
                cx,
                2.0 * DT,
                &[],
                &loads,
                &[],
                &geometry,
                &mut molecular,
                &internal,
                &links,
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
            assert_eq!(held, retained);
        })
    });
}
