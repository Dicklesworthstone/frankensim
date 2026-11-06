use super::*;
use std::sync::Arc;

use crate::burnup::transport::schedule::NuclideTransportChange;
use crate::decay::production::{DecayProductionChange, DecayProductionSchedule};
use crate::decay::stochastic::transport::production::{
    StochasticDecayProductionTransportBudget as FedBudget,
    StochasticDecayProductionTransportState as FedState,
};
use crate::decay::stochastic::transport::schedule::StochasticDecayTransportSchedule as MovementSchedule;
use crate::decay::stochastic::transport::{
    NuclideOutflow, NuclideTransfer, NuclideTransport, StochasticDecayTransportBudget as Budget,
    StochasticDecayTransportLimits, StochasticDecayTransportState as State,
};
use crate::decay::stochastic::{StochasticDecayBudget, StochasticDecayInventory, StreamKey};
use crate::decay::{DecayBranch, DecayModel, DecayNuclide};
use crate::radiation::{
    RadiationDeposit, RadiationResponse, RadiationResponseLimits, RadiationSourceResponse,
};
use crate::thermal::RadiationReceiver;
use fs_la::decay::DecayLimits;
use fs_qty::{Time, Volume, VolumetricFlowRate};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, DepositionError, DepositionSimulation,
    Event, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermochem::decay::DecayLawV1;

type Heating = StochasticDecayTransportHeating;
type Source = StochasticDecayTransportCellSource;
type FedHeating = StochasticDecayProductionTransportHeating;
type FedSource = StochasticDecayProductionTransportCellSource;

// Manufactured event energies and the original 2 kg, c=4 J/(kg K) calorimetry.
// The explicit flow volumes match the two thermal cells; no carrier heat is added.
struct Caloric(f64);
impl ThermalLaw for Caloric {
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        let temperature = energy / 4.0;
        if temperature.is_finite() && temperature > 0.0 && temperature <= self.0 {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.25
    }
}
static LAW: Caloric = Caloric(1e6);
static LIMITED: Caloric = Caloric(300.5);

fn volumes() -> [Volume; 2] {
    [Volume::new(0.01), Volume::new(0.02)]
}
fn controls(max_dt_s: f64) -> StepConfig {
    StepConfig {
        max_dt_s,
        energy_tolerance_j: 1e-8,
    }
}
fn budget() -> Budget {
    Budget {
        decay: StochasticDecayBudget {
            maximum_events: 1024,
            energy_tolerance_j: 1e-8,
        },
        maximum_movements: 2048,
        maximum_clock_checks: 1_000_000,
    }
}
fn fed_budget() -> FedBudget {
    FedBudget {
        transport: budget(),
        maximum_births: 128,
        maximum_schedule_changes: 8,
    }
}
fn cells() -> Vec<(CellId, MaterialId)> {
    vec![(CellId(1), MaterialId(1)), (CellId(2), MaterialId(1))]
}
fn thermal(law: &'static dyn ThermalLaw) -> DepositionSimulation<'static> {
    let events = (0..2)
        .map(|index| {
            let id = index as u64 + 1;
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(id),
                    segment_id: id,
                    body_id: 1,
                    material_id: MaterialId(1),
                    mass_kg: 2.0,
                    specific_h_j_kg: 1200.0,
                    geometry: CellGeometry {
                        start_m: [id as f64, 0.0, 0.0],
                        end_m: [id as f64, 1.0, 0.0],
                        volume_m3: volumes()[index].value(),
                        surface_area_m2: 1.0,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            }
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 2,
            max_links: 2,
            max_events: 4,
            max_substeps: 10_000,
        },
    )
    .unwrap()
    .advance_to(0.0, controls(10.0), || Ok(()))
    .unwrap()
    .0
}
fn frame(transfer_m3_s: f64, extraction_m3_s: f64) -> NuclideTransport {
    NuclideTransport {
        transfers: vec![
            NuclideTransfer::from_volumetric_flow(
                0,
                1,
                &volumes(),
                VolumetricFlowRate::new(transfer_m3_s),
                3,
            )
            .unwrap(),
        ],
        outflows: vec![
            NuclideOutflow::from_volumetric_flow(
                1,
                &volumes(),
                VolumetricFlowRate::new(extraction_m3_s),
                3,
            )
            .unwrap(),
        ],
    }
}
fn state(populations: [[u64; 3]; 2]) -> State {
    let model = Arc::new(
        DecayModel::new(
            "manufactured-valve-and-half-life-heat-v1".into(),
            vec![
                DecayNuclide {
                    id: "parent".into(),
                    law: DecayLawV1::from_half_life(Time::new(1.0)).unwrap(),
                },
                DecayNuclide {
                    id: "daughter".into(),
                    law: DecayLawV1::from_half_life(Time::new(0.5)).unwrap(),
                },
                DecayNuclide {
                    id: "stable".into(),
                    law: DecayLawV1::stable(),
                },
            ],
            vec![
                DecayBranch {
                    parent: 0,
                    daughter: 1,
                    probability: 1.0,
                    energy_per_decay: EnergyChannels {
                        deposited_j: 8.0,
                        radiation_j: 2.0,
                        escaped_j: 1.0,
                    },
                },
                DecayBranch {
                    parent: 1,
                    daughter: 2,
                    probability: 1.0,
                    energy_per_decay: EnergyChannels {
                        deposited_j: 4.0,
                        radiation_j: 1.0,
                        escaped_j: 0.5,
                    },
                },
            ],
            DecayLimits::default(),
        )
        .unwrap(),
    );
    let inventories = populations
        .into_iter()
        .enumerate()
        .map(|(index, population)| {
            StochasticDecayInventory::new(
                Arc::clone(&model),
                population.to_vec(),
                StreamKey {
                    seed: 7,
                    kernel: 0x484c46,
                    tile: index as u32 + 1,
                },
            )
            .unwrap()
        })
        .collect();
    State::new(
        inventories,
        frame(0.0, 0.0),
        StreamKey {
            seed: 99,
            kernel: 0x484c46,
            tile: 0,
        },
        StochasticDecayTransportLimits::default(),
    )
    .unwrap()
}
fn valve_history() -> (Arc<State>, f64, f64, f64) {
    let original = state([[1, 0, 0], [0; 3]]);
    // The valve is authored before this fixture's already sampled first decay.
    let opening = original.inventories()[0].next_event_time_s().unwrap() * 0.25;
    let schedule = Arc::new(
        MovementSchedule::new(
            &original,
            frame(0.0, 0.0),
            vec![NuclideTransportChange {
                time_s: opening,
                transport: frame(1e4, 0.0),
            }],
            1,
        )
        .unwrap(),
    );
    let original = Arc::new(original.with_transport_schedule(schedule).unwrap());
    let mut probe = original.as_ref().clone();
    let changed = probe.advance_to(opening, budget(), || Ok(())).unwrap();
    assert_eq!(changed.transport_changes, 1);
    assert_eq!(changed.decay_events, 0);
    assert_eq!(changed.movement_events, 0);
    let movement = probe.next_event_time_s().unwrap();
    let moved = probe.advance_to(movement, budget(), || Ok(())).unwrap();
    assert_eq!(moved.movement_events, 1);
    assert_eq!(moved.decay_events, 0);
    let decay = probe.next_event_time_s().unwrap();
    let decayed = probe.advance_to(decay, budget(), || Ok(())).unwrap();
    assert_eq!(decayed.decay[0].events, 0);
    assert_eq!(decayed.decay[1].branch_decays, [1, 0]);
    (original, opening, movement, decay)
}
fn fed_history() -> FedState {
    let transport = state([[0; 3]; 2]);
    let movement = Arc::new(
        MovementSchedule::new(
            &transport,
            frame(0.0, 0.0),
            vec![
                NuclideTransportChange {
                    time_s: 0.25,
                    transport: frame(0.03, 0.0),
                },
                NuclideTransportChange {
                    time_s: 1.0,
                    transport: frame(0.0, 0.02),
                },
                NuclideTransportChange {
                    time_s: 3.0,
                    transport: frame(0.01, 0.02),
                },
            ],
            3,
        )
        .unwrap(),
    );
    let supply = Arc::new(
        DecayProductionSchedule::new(
            Arc::clone(transport.inventories()[0].model()),
            vec![8.0, 0.0, 0.0],
            vec![DecayProductionChange {
                time_s: 1.0,
                rates_nuclei_s: vec![0.0; 3],
            }],
            1,
        )
        .unwrap(),
    );
    FedState::poisson(transport, vec![Some(supply), None])
        .unwrap()
        .with_transport_schedule(movement)
        .unwrap()
}
fn produced(state: &FedState) -> u64 {
    state.cumulative_produced_nuclei(0).unwrap().iter().sum()
}
fn assert_balance(state: &FedState) {
    let retained: u64 = state
        .inventories()
        .iter()
        .flat_map(|cell| cell.populations())
        .sum();
    let outflow: u64 = state
        .transport_state()
        .cumulative_outflow_counts()
        .iter()
        .flatten()
        .sum();
    assert_eq!(produced(state), retained + outflow);
}

#[test]
fn g1_valve_boundary_adds_no_heat_or_nuclear_draws_then_the_recipient_gets_decay_heat() {
    let (original, opening, movement, decay) = valve_history();
    let initial_clock = original.inventories()[0].next_event_time_s();
    let initial_rng = original.inventories()[0].rng_checkpoint();
    let mut expected = original.as_ref().clone();
    let mut material = thermal(&LAW);
    let mut heating = Heating::new(
        &material,
        Source::bind(original, cells()).unwrap(),
        budget(),
    )
    .unwrap();
    for end in [opening, movement, decay] {
        let step = expected.advance_to(end, budget(), || Ok(())).unwrap();
        let (next_material, next, audit) = heating
            .advance_to(&material, end, controls(100.0), || Ok(()))
            .unwrap();
        assert_eq!(next.source(CellId(1)).unwrap().state().as_ref(), &expected);
        if end == opening {
            assert_eq!(
                next.last_step(CellId(1))
                    .unwrap()
                    .transport()
                    .transport_changes,
                1
            );
            assert_eq!(step.transport_changes, 1);
            assert_eq!(expected.inventories()[0].next_event_time_s(), initial_clock);
            assert_eq!(expected.inventories()[0].rng_checkpoint(), initial_rng);
            assert_eq!(expected.transport().transfers[0].rates_per_s, [1e6; 3]);
        }
        if end < decay {
            assert_eq!(audit.external_heat_j, 0.0);
            assert_eq!(next_material.temperature_k(CellId(1)).unwrap(), 300.0);
            assert_eq!(next_material.temperature_k(CellId(2)).unwrap(), 300.0);
        } else {
            assert_eq!(next.audit().decays, 1.0);
            assert!((audit.external_heat_j - 8.0).abs() < 1e-10);
            assert_eq!(next_material.temperature_k(CellId(1)).unwrap(), 300.0);
            assert!((next_material.temperature_k(CellId(2)).unwrap() - 301.0).abs() < 1e-10);
        }
        material = next_material;
        heating = next;
    }
    for cell in [CellId(1), CellId(2)] {
        assert_eq!(material.cells()[&cell].mass_kg, 2.0);
    }
}

#[test]
fn g1_joint_feed_shutdown_and_valves_keep_local_and_returned_radiation_energy_separate() {
    let initial = thermal(&LAW);
    let original = Arc::new(fed_history());
    let mut expected = original.as_ref().clone();
    let full = expected.advance_to(8.0, fed_budget(), || Ok(())).unwrap();
    assert_eq!(full.transport.transport_changes, 3);
    assert!(full.transport.outflow_counts.iter().flatten().sum::<u64>() > 0);
    let response = Arc::new(
        RadiationResponse::new(
            "manufactured-scheduled-half-radiation-return-v1".into(),
            (1..=2)
                .map(|source_id| RadiationSourceResponse {
                    source_id,
                    deposits: vec![RadiationDeposit {
                        receiver_id: 3 - source_id,
                        fraction: 0.5,
                    }],
                    escape_fraction: 0.5,
                })
                .collect(),
            RadiationResponseLimits {
                max_sources: 2,
                max_receivers: 2,
                max_transfers: 2,
            },
        )
        .unwrap(),
    );
    let heating = FedHeating::new(
        &initial,
        FedSource::bind(original, cells()).unwrap(),
        fed_budget(),
    )
    .unwrap()
    .with_radiation_response(
        &initial,
        response,
        cells()
            .into_iter()
            .map(|(cell_id, material_id)| RadiationReceiver {
                cell_id,
                material_id,
            })
            .collect(),
    )
    .unwrap();
    let (middle_material, middle, first_audit) = heating
        .advance_to(&initial, 1.0, controls(100.0), || Ok(()))
        .unwrap();
    let shutdown = middle.source(CellId(1)).unwrap().state();
    assert_eq!(
        shutdown.schedule(0).unwrap().rates_at(1.0).unwrap(),
        [0.0; 3]
    );
    assert_eq!(
        shutdown.transport_state().transport().transfers[0].rates_per_s,
        [0.0; 3]
    );
    assert_eq!(
        shutdown.transport_state().transport().outflows[0].rates_per_s,
        [1.0; 3]
    );
    assert_eq!(
        shutdown
            .transport_state()
            .transport_schedule()
            .unwrap()
            .next_after(1.0),
        Some(3.0)
    );
    let births = produced(shutdown);
    assert!(births > 0);
    let (after, continued, late_audit) = middle
        .advance_to(&middle_material, 8.0, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(
        continued.source(CellId(1)).unwrap().state().as_ref(),
        &expected
    );
    assert_eq!(continued.audit().produced_nuclei, 0.0);
    assert!(continued.audit().decays > 0.0);
    assert_eq!(produced(&expected), births);
    assert_balance(&expected);
    let energy: Vec<_> = expected
        .inventories()
        .iter()
        .map(|cell| cell.cumulative_released())
        .collect();
    let heat = energy
        .iter()
        .map(|e| e.deposited_j + 0.5 * e.radiation_j)
        .sum::<f64>();
    assert!((first_audit.external_heat_j + late_audit.external_heat_j - heat).abs() < 1e-8);
    for index in 0..2 {
        let expected_k =
            300.0 + (energy[index].deposited_j + 0.5 * energy[1 - index].radiation_j) / 8.0;
        assert!(
            (after.temperature_k(CellId(index as u64 + 1)).unwrap() - expected_k).abs() < 1e-10
        );
    }
    let radiated = energy.iter().map(|e| e.radiation_j).sum::<f64>();
    assert!(
        (middle.audit().radiation_deposited_j + continued.audit().radiation_deposited_j
            - 0.5 * radiated)
            .abs()
            < 1e-10
    );
    assert!(
        (middle.audit().radiation_escaped_j + continued.audit().radiation_escaped_j
            - 0.5 * radiated)
            .abs()
            < 1e-10
    );
}

#[test]
fn g3_valve_and_supply_clocks_survive_partition_restart_and_reversed_cell_ids() {
    let initial = thermal(&LAW);
    let original = Arc::new(fed_history());
    let mut expected = original.as_ref().clone();
    expected.advance_to(8.0, fed_budget(), || Ok(())).unwrap();
    let heating = FedHeating::new(
        &initial,
        FedSource::bind(Arc::clone(&original), cells()).unwrap(),
        fed_budget(),
    )
    .unwrap();
    let (whole_material, whole, _) = heating
        .advance_to(&initial, 8.0, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(whole.audit().produced_nuclei, produced(&expected) as f64);
    let (mut split_material, mut split) = (initial.clone(), heating);
    let mut saved = None;
    for end in [0.1, 0.25, 0.6, 1.0, 2.3, 3.0, 8.0] {
        let (next_material, next, _) = split
            .advance_to(&split_material, end, controls(0.37), || Ok(()))
            .unwrap();
        split_material = next_material;
        split = next;
        if end == 0.6 {
            saved = Some((split_material.clone(), split.clone()));
        }
    }
    let (saved_material, saved) = saved.unwrap();
    let (restarted_material, restarted, _) = saved
        .advance_to(&saved_material, 8.0, controls(10.0), || Ok(()))
        .unwrap();
    for observed in [&whole, &split, &restarted] {
        assert_eq!(
            observed.source(CellId(1)).unwrap().state().as_ref(),
            &expected
        );
        assert!(Arc::ptr_eq(
            observed.source(CellId(1)).unwrap().state(),
            observed.source(CellId(2)).unwrap().state()
        ));
    }
    assert_balance(&expected);
    for (index, inventory) in expected.inventories().iter().enumerate() {
        let expected_k = 300.0 + inventory.cumulative_released().deposited_j / 8.0;
        for material in [&whole_material, &split_material, &restarted_material] {
            assert!(
                (material.temperature_k(CellId(index as u64 + 1)).unwrap() - expected_k).abs()
                    < 1e-10
            );
        }
    }
    let mut reverse = cells();
    reverse.reverse();
    let mut sources = FedSource::bind(original, reverse).unwrap();
    sources.reverse();
    let swapped = FedHeating::new(&initial, sources, fed_budget()).unwrap();
    let (swapped_material, swapped, _) = swapped
        .advance_to(&initial, 8.0, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(swapped.source(CellId(1)).unwrap().compartment_index(), 1);
    assert_eq!(
        swapped.source(CellId(1)).unwrap().state().as_ref(),
        &expected
    );
    for id in [1, 2] {
        assert!(
            (swapped_material.temperature_k(CellId(id)).unwrap()
                - whole_material.temperature_k(CellId(3 - id)).unwrap())
            .abs()
                < 1e-10
        );
    }
}

#[test]
fn g4_refused_flow_work_heat_and_final_checkpoint_restore_rates_clocks_and_materials() {
    let (original, opening, _, decay) = valve_history();
    let (_, opening_step) = original.trial_to(opening, budget(), || Ok(())).unwrap();
    assert_eq!(opening_step.transport_changes, 1);
    assert!(opening_step.clock_checks > 1);
    let initial = thermal(&LAW);
    let strict = Budget {
        maximum_clock_checks: opening_step.clock_checks - 1,
        ..budget()
    };
    let refused = Heating::new(
        &initial,
        Source::bind(Arc::clone(&original), cells()).unwrap(),
        strict,
    )
    .unwrap();
    assert!(matches!(
        refused.advance_to(&initial, opening, controls(100.0), || Ok(())),
        Err(FissionThermalError::Decay(DecayError::Budget(_)))
    ));
    assert_eq!(
        refused.source(CellId(1)).unwrap().state().as_ref(),
        original.as_ref()
    );
    let limited_material = thermal(&LIMITED);
    let heating = Heating::new(
        &limited_material,
        Source::bind(Arc::clone(&original), cells()).unwrap(),
        budget(),
    )
    .unwrap();
    assert!(
        heating
            .advance_to(&limited_material, decay, controls(100.0), || Ok(()))
            .is_err()
    );
    for cell in [CellId(1), CellId(2)] {
        assert_eq!(limited_material.temperature_k(cell).unwrap(), 300.0);
        assert_eq!(
            heating.source(cell).unwrap().state().as_ref(),
            original.as_ref()
        );
    }
    assert_eq!(limited_material.time_s(), 0.0);
    let heating = Heating::new(
        &initial,
        Source::bind(Arc::clone(&original), cells()).unwrap(),
        budget(),
    )
    .unwrap();
    let mut checks = 0;
    let (_, expected, _) = heating
        .advance_to(&initial, decay, controls(100.0), || {
            checks += 1;
            Ok(())
        })
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        heating.advance_to(&initial, decay, controls(100.0), || {
            seen += 1;
            if seen == checks {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(FissionThermalError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(seen, checks);
    assert_eq!(
        heating.source(CellId(1)).unwrap().state().as_ref(),
        original.as_ref()
    );
    assert_eq!(
        heating
            .source(CellId(1))
            .unwrap()
            .state()
            .transport()
            .transfers[0]
            .rates_per_s,
        [0.0; 3]
    );
    assert_eq!(initial.time_s(), 0.0);
    let (_, retry, _) = heating
        .advance_to(&initial, decay, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(
        retry.source(CellId(1)).unwrap().state(),
        expected.source(CellId(1)).unwrap().state()
    );
}
