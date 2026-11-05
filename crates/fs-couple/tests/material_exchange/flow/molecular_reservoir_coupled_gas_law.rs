//! G3/G4: selected wall laws retain native adaptive bath trajectories and refusals.
use super::*;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::{Contact, NetworkError};
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::exchange::gas_law::GasSurface;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state::flow::molecular::reservoir::coupled::gas_law::CoupledMolecularReservoirGasLawError;
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

#[test]
fn selected_conduction_preserves_complete_adaptive_bath_and_reaction_trajectory() {
    fixture([0.2, 0.8], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut original_molecular) = molecular_fixture(MASS, 1.0, -50000.0);
            let (_, mut molecular) = molecular_fixture(MASS, 1.0, -50000.0);
            let held = [
                bath(&bank, 0.9, 430.0, 1.2e5),
                bath(&bank, 0.1, 370.0, 0.9e5),
            ];
            let retained = held.clone();
            let reservoirs = [
                ReservoirMolecularLink {
                    area_m2: 0.08,
                    ..bath_link(&held[0])
                },
                ReservoirMolecularLink {
                    id: 202,
                    parcel: 20,
                    area_m2: 0.06,
                    ..bath_link(&held[1])
                },
            ];
            let internal = [molecular_link()];
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.7, 390.0, 1.1e5), (0.3, 410.0, 0.95e5)], false),
            );
            let mut loads = drives(state.chambers());
            loads[0].heat_j = 0.015625;
            loads[1].heat_j = 0.03125;
            loads[0].external_force_n = 2.0;
            loads[1].external_force_n = -1.5;
            let law = reaction(50.0);
            let reactions = [ChamberReaction {
                chamber: 10,
                reaction: &law,
                controls: chemical(),
            }];
            let patches = gas_fixture::patches();
            let geometry = surfaces(&patches);
            let mut time = loose_time(DT / 2.0);
            time.state.maximum_evaluations = 6;
            let (expected, original) = state
                .coupled_molecular_reservoir_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &patches,
                    &mut original_molecular,
                    &internal,
                    &reservoirs,
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    time,
                )
                .unwrap();
            let report = state
                .advance_with_coupled_molecular_reservoir_and_gas_law(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &reactions,
                    &geometry,
                    &mut molecular,
                    &internal,
                    &reservoirs,
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    time,
                    |contact, thermal, material, gas| {
                        let patch = *patches
                            .iter()
                            .find(|patch| patch.contact == contact)
                            .unwrap();
                        prepared.resolve_conductive_gas_contact(patch, thermal, material, gas)
                    },
                )
                .unwrap();
            same_owners(&state, &expected);
            assert_eq!(report.balance, original.balance);
            assert_eq!(report.estimated_error, original.estimated_error);
            assert_eq!(report.evaluations, 6);
            assert_eq!(report.accepted_intervals, 2);
            assert_eq!(report.rejected_intervals, 0);
            assert_eq!(report.steps.len(), 4);
            for (selected, default) in report.steps.iter().zip(&original.steps) {
                assert_eq!(
                    (selected.start_s, selected.end_s),
                    (default.start_s, default.end_s)
                );
                assert_eq!(selected.advance.flow, default.advance.flow);
                assert_eq!(selected.advance.balance, default.advance.balance);
                assert_eq!(
                    selected.advance.exchange.step,
                    default.advance.exchange.step
                );
                let a = &selected.advance.exchange.report.source;
                let b = &default.advance.exchange.report.source;
                assert_eq!(a.reactions, b.reactions);
                assert!(a.reactions[0].reaction.extent_mol > 0.0);
                for (selected, default) in a.chambers.iter().zip(&b.chambers) {
                    assert_eq!(selected.boundary, default.boundary);
                    assert_eq!(selected.surfaces, default.surfaces);
                }
            }
            for chamber in state.chambers().values() {
                assert!(chamber.velocity_m_s().abs() > 1e-8);
                assert_eq!(chamber.last_window_duration_s(), Some(DT));
                assert!(matches!(
                    chamber.audit_last_window().unwrap(),
                    ModuleWindowAudit::Unavailable { .. }
                ));
            }
            fixed_budgets(&report.balance);
            assert_eq!(held, retained);
        });
    });
}

// No Clone implementation: accepted source records must move through the
// adaptive controller and window balance without erasure or extra bounds.
#[derive(Debug, PartialEq)]
struct HeatReceipt {
    material: ResolvedMaterialStatePoint,
    gas: ResolvedThermochemicalGasState,
    base_conductance_w_k: f64,
    endpoint_s: f64,
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
fn every_accepted_fine_step_retains_current_composition_and_temperature_heat_sources() {
    fixture([0.8, 0.8], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let (_, mut replay_molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let held = bath(&bank, 0.1, 370.0, 1e5);
            let reservoirs = [bath_link(&held)];
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
            let resolve = |contact: Contact,
                           thermal: &DepositionSimulation<'_>,
                           material: &MaterialExchangeSample,
                           gas: &ResolvedThermochemicalGasState| {
                let patch = *patches
                    .iter()
                    .find(|patch| patch.contact == contact)
                    .unwrap();
                let base =
                    prepared.resolve_conductive_gas_contact(patch, thermal, material, gas)?;
                Ok::<_, MaterialNetworkError>(ConductanceSample {
                    conductance_w_k: base.conductance_w_k
                        * (1.0 + mole_a(gas))
                        * gas.state().temperature
                        / 400.0,
                    source: HeatReceipt {
                        material: base.source,
                        gas: gas.clone(),
                        base_conductance_w_k: base.conductance_w_k,
                        endpoint_s: thermal.time_s(),
                    },
                })
            };
            let mut time = loose_time(DT / 2.0);
            time.state.maximum_evaluations = 6;
            let (next, report) = state
                .coupled_molecular_reservoir_with_gas_law_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &geometry,
                    &mut molecular,
                    &internal,
                    &reservoirs,
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                    time,
                    &resolve,
                )
                .unwrap();
            assert_eq!(report.steps.len(), 4);
            let mut cursor = MaterialExchangeState::from_parts(
                state.thermal().clone(),
                state.circuit().clone(),
                state.chambers().clone(),
            );
            let mut previous_composition = BTreeMap::new();
            for accepted in &report.steps {
                let dt = accepted.end_s - accepted.start_s;
                let parcels = cursor
                    .chambers()
                    .iter()
                    .map(|(&id, chamber)| (id, chamber.gas().clone()))
                    .collect();
                let endpoints = loads
                    .iter()
                    .map(|drive| (drive.chamber, drive.piston_solve.gas))
                    .collect();
                let native = replay_molecular
                    .network_candidate(
                        &parcels,
                        &internal,
                        &reservoirs,
                        dt,
                        &endpoints,
                        molecular_controls(),
                        || false,
                    )
                    .unwrap();
                let (oracle, replay) = cursor
                    .molecular_reservoir_with_gas_law_candidate(
                        &prepared,
                        cx,
                        accepted.end_s,
                        &[],
                        &loads,
                        &[],
                        &geometry,
                        &mut replay_molecular,
                        &internal,
                        &reservoirs,
                        molecular_controls(),
                        controls(),
                        &[0.0; 2],
                        &resolve,
                    )
                    .unwrap();
                assert_eq!(&accepted.advance.flow, native.step());
                assert_eq!(accepted.advance.flow, replay.flow);
                assert_eq!(accepted.advance.exchange.step, replay.exchange.step);
                for (chamber, replay) in accepted
                    .advance
                    .exchange
                    .report
                    .source
                    .chambers
                    .iter()
                    .zip(&replay.exchange.report.source.chambers)
                {
                    assert_eq!(chamber.boundary, replay.boundary);
                    assert_eq!(chamber.surfaces, replay.surfaces);
                    let staged = &native.parcels()[&chamber.id];
                    let original = cursor.chambers()[&chamber.id].gas();
                    let final_gas = oracle.chambers()[&chamber.id].gas();
                    let mut heat = 0.0;
                    assert!(!chamber.surfaces.is_empty());
                    for surface in &chamber.surfaces {
                        let (contact, receipt) = &surface.source;
                        assert_eq!(receipt.endpoint_s, accepted.end_s);
                        assert_eq!(&receipt.gas, &chamber.boundary.piston.midpoint);
                        assert_eq!(mole_a(&receipt.gas), mole_a(staged.state()));
                        assert!((mole_a(&receipt.gas) - mole_a(original.state())).abs() > 1e-6);
                        if let Some(previous) =
                            previous_composition.insert(chamber.id, mole_a(&receipt.gas))
                        {
                            assert!((mole_a(&receipt.gas) - previous).abs() > 1e-6);
                        }
                        let expected = receipt.base_conductance_w_k
                            * (1.0 + mole_a(&receipt.gas))
                            * receipt.gas.state().temperature
                            / 400.0;
                        close(surface.conductance_w_k, expected, 1e-14);
                        close(
                            surface.heat_j,
                            dt * expected
                                * (surface.port.temperature_k - receipt.gas.state().temperature),
                            1e-12,
                        );
                        let cell = accepted
                            .advance
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
                        assert_eq!(Some(&receipt.material), cell.conductivity.as_ref());
                        heat += surface.heat_j;
                    }
                    // Original bath enthalpy and molecular heat already entered
                    // native.parcels(); only selected wall heat changes this
                    // fixed-volume, nonreacting, constant-Cv feedback energy.
                    let capacity =
                        staged.mass_kg() * 2.5 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MASS;
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
                cursor = oracle;
            }
            assert_eq!(next.time_s(), cursor.time_s());
            assert_eq!(next.chambers(), cursor.chambers());
            assert_eq!(next.circuit(), cursor.circuit());
            assert_eq!(next.thermal().cells(), cursor.thermal().cells());
            for chamber in next.chambers().values() {
                assert_eq!(chamber.last_window_duration_s(), Some(DT));
            }
            fixed_budgets(&report.balance);
            assert_eq!(state.time_s(), 0.0);
        });
    });
}

#[derive(Debug)]
enum HeatRefusal {
    OutOfDomain { contact: u64 },
    Native(MaterialNetworkError),
}
impl core::fmt::Display for HeatRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::OutOfDomain { contact } => {
                write!(f, "held-bath heat law outside domain on {contact}")
            }
            Self::Native(error) => write!(f, "original gas heat law: {error}"),
        }
    }
}

#[test]
fn late_typed_heat_refusal_is_fatal_and_cancellation_preserves_published_histories() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, gate| {
            let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let held = bath(&bank, 0.9, 430.0, 1e5);
            let retained = held.clone();
            let reservoirs = [bath_link(&held)];
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
                    &reservoirs,
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
            let cutoff = 2.0 * DT;
            let target = 3.0 * DT;
            let refused = std::cell::Cell::new(0usize);
            let failure = state.advance_with_coupled_molecular_reservoir_and_gas_law(
                &prepared,
                cx,
                target,
                &[],
                &loads,
                &[],
                &geometry,
                &mut molecular,
                &internal,
                &reservoirs,
                molecular_controls(),
                controls(),
                &[0.0; 2],
                loose_time(DT),
                |contact, thermal, material, gas| {
                    if thermal.time_s() > cutoff {
                        refused.set(refused.get() + 1);
                        return Err(HeatRefusal::OutOfDomain {
                            contact: contact.id,
                        });
                    }
                    let patch = *patches
                        .iter()
                        .find(|patch| patch.contact == contact)
                        .unwrap();
                    prepared
                        .resolve_conductive_gas_contact(patch, thermal, material, gas)
                        .map_err(HeatRefusal::Native)
                },
            );
            match failure {
                Err(CoupledMolecularReservoirGasLawError::Step {
                    start_s,
                    end_s,
                    error:
                        MaterialGasLawFlowError::Exchange(ReactiveExchangeError::Exchange(
                            CircuitEnthalpyError::Source(NetworkError::Boundary {
                                error:
                                    ThermalPistonError::Interface(HeatRefusal::OutOfDomain {
                                        contact: 101,
                                    }),
                                ..
                            }),
                        )),
                }) => {
                    // The first coarse/two-fine comparison has already been
                    // accepted privately. A source-domain failure must retain
                    // its original type and must not retry a shorter interval.
                    assert_eq!(start_s, cutoff);
                    assert_eq!(end_s, target);
                }
                other => panic!(
                    "expected typed late heat refusal, got {}",
                    match other {
                        Ok(_) => "an accepted window".into(),
                        Err(error) => format!("{error:?}"),
                    }
                ),
            }
            assert_eq!(refused.get(), 1);
            same_owners(&state, &before);
            let observed_fine_prefix = std::cell::Cell::new(false);
            let cancelled_at = std::cell::Cell::new(0.0);
            let cancelled = state.advance_with_coupled_molecular_reservoir_and_gas_law(
                &prepared,
                cx,
                target,
                &[],
                &loads,
                &[],
                &geometry,
                &mut molecular,
                &internal,
                &reservoirs,
                molecular_controls(),
                controls(),
                &[0.0; 2],
                loose_time(DT),
                |contact, thermal, material, gas| {
                    let patch = *patches
                        .iter()
                        .find(|patch| patch.contact == contact)
                        .unwrap();
                    let sample = prepared
                        .resolve_conductive_gas_contact(patch, thermal, material, gas)
                        .map_err(HeatRefusal::Native)?;
                    if thermal.time_s() > DT && thermal.time_s() < cutoff {
                        observed_fine_prefix.set(true);
                    }
                    if thermal.time_s() > cutoff {
                        cancelled_at.set(thermal.time_s());
                        gate.request();
                    }
                    Ok::<_, HeatRefusal>(sample)
                },
            );
            assert!(matches!(
                cancelled,
                Err(CoupledMolecularReservoirGasLawError::Cancelled)
            ));
            assert!(observed_fine_prefix.get());
            assert_eq!(cancelled_at.get(), target);
            same_owners(&state, &before);
            assert_eq!(held, retained);
        });
    });
}

#[test]
fn native_species_withdrawal_refines_through_the_selected_heat_law_error_type() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut molecular) = molecular_fixture(MASS, 0.0, 0.0);
            let held = bath(&bank, 0.1, 400.0, 1e5);
            let reservoirs = [bath_link(&held)];
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.8, 400.0, 1e5); 2], true),
            );
            let loads = drives(state.chambers());
            let patches = gas_fixture::patches();
            let geometry = surfaces(&patches);
            let resolve = |contact: Contact,
                           thermal: &DepositionSimulation<'_>,
                           material: &MaterialExchangeSample,
                           gas: &ResolvedThermochemicalGasState| {
                let patch = *patches
                    .iter()
                    .find(|patch| patch.contact == contact)
                    .unwrap();
                prepared
                    .resolve_conductive_gas_contact(patch, thermal, material, gas)
                    .map_err(HeatRefusal::Native)
            };
            let mut physical = molecular_controls();
            // Initially dm_A/m_A = (D A dt / V L) (0.8-0.1)/0.8
            // = 0.000875. A half comparison interval fits 0.0006.
            physical.physical.exchange.maximum_species_fraction = 0.0006;
            let coarse = state.molecular_reservoir_with_gas_law_candidate(
                &prepared,
                cx,
                DT,
                &[],
                &loads,
                &[],
                &geometry,
                &mut molecular,
                &[],
                &reservoirs,
                physical,
                controls(),
                &[0.0; 2],
                &resolve,
            );
            assert!(matches!(
                coarse,
                Err(MaterialGasLawFlowError::Molecular(
                    ParcelMolecularError::Transport(GasTransportError::Exchange(
                        SpeciesNetworkError::Withdrawal { .. }
                    ))
                ))
            ));
            let (next, report) = state
                .coupled_molecular_reservoir_with_gas_law_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &loads,
                    &[],
                    &geometry,
                    &mut molecular,
                    &[],
                    &reservoirs,
                    physical,
                    controls(),
                    &[0.0; 2],
                    loose_time(DT),
                    &resolve,
                )
                .unwrap();
            assert!(report.rejected_intervals >= 1);
            assert!(report.accepted_intervals >= 2);
            assert!(report.evaluations > 3 * report.accepted_intervals);
            assert!(
                species_mass(next.chambers()[&10].gas(), "A")
                    < species_mass(state.chambers()[&10].gas(), "A")
            );
            assert!(report.steps.iter().all(|step| {
                step.advance
                    .exchange
                    .report
                    .source
                    .chambers
                    .iter()
                    .all(|chamber| !chamber.surfaces.is_empty())
            }));
            fixed_budgets(&report.balance);
            assert_eq!(next.time_s(), DT);
            assert_eq!(state.time_s(), 0.0);
        });
    });
}
