use super::*;
use std::sync::Arc;

use crate::burnup::transport::{NuclideTransfer, NuclideTransport};
use crate::decay::stochastic::transport::{
    StochasticDecayTransportBudget, StochasticDecayTransportLimits, StochasticDecayTransportState,
};
use crate::decay::stochastic::{StochasticDecayBudget, StochasticDecayInventory, StreamKey};
use crate::decay::{DecayBranch, DecayModel, DecayNuclide};
use fs_la::decay::DecayLimits;
use fs_qty::Time;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, DepositionError, DepositionSimulation,
    Event, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermochem::decay::DecayLawV1;

// Original stochastic calorimetry fixture: 2 kg and c=4 J/(kg K) in each cell.
// Source coefficients, movement rates and joules/event are manufactured.
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

fn controls(max_dt_s: f64) -> StepConfig {
    StepConfig {
        max_dt_s,
        energy_tolerance_j: 1e-8,
    }
}

fn budget() -> StochasticDecayTransportBudget {
    StochasticDecayTransportBudget {
        decay: StochasticDecayBudget {
            maximum_events: 1024,
            energy_tolerance_j: 1e-8,
        },
        maximum_movements: 1024,
        maximum_clock_checks: 1_000_000,
    }
}

fn thermal(law: &'static dyn ThermalLaw) -> DepositionSimulation<'static> {
    let events = (1..=2)
        .map(|id| TimedEvent {
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
                    volume_m3: 0.01,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
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

fn state(populations: [[u64; 2]; 2], rates: [[f64; 2]; 2]) -> StochasticDecayTransportState {
    let model = Arc::new(
        DecayModel::new(
            "manufactured-stochastic-transport-heat-v1".into(),
            vec![
                DecayNuclide {
                    id: "parent".into(),
                    law: DecayLawV1::from_half_life(Time::new(1.0)).unwrap(),
                },
                DecayNuclide {
                    id: "stable".into(),
                    law: DecayLawV1::stable(),
                },
            ],
            vec![DecayBranch {
                parent: 0,
                daughter: 1,
                probability: 1.0,
                energy_per_decay: EnergyChannels {
                    deposited_j: 8.0,
                    radiation_j: 2.0,
                    escaped_j: 1.0,
                },
            }],
            DecayLimits::default(),
        )
        .unwrap(),
    );
    let inventories = populations
        .into_iter()
        .enumerate()
        .map(|(cell, population)| {
            StochasticDecayInventory::new(
                Arc::clone(&model),
                population.to_vec(),
                StreamKey {
                    seed: 7,
                    kernel: 0x484c46,
                    tile: cell as u32 + 1,
                },
            )
            .unwrap()
        })
        .collect();
    let transfers = rates
        .into_iter()
        .enumerate()
        .map(|(cell, rates)| NuclideTransfer {
            from_compartment: cell,
            to_compartment: 1 - cell,
            rates_per_s: rates.to_vec(),
        })
        .collect();
    StochasticDecayTransportState::new(
        inventories,
        NuclideTransport {
            transfers,
            outflows: Vec::new(),
        },
        StreamKey {
            seed: 99,
            kernel: 0x484c46,
            tile: 0,
        },
        StochasticDecayTransportLimits::default(),
    )
    .unwrap()
}

fn cells() -> Vec<(CellId, MaterialId)> {
    vec![(CellId(1), MaterialId(1)), (CellId(2), MaterialId(1))]
}

fn bind(state: Arc<StochasticDecayTransportState>) -> Vec<StochasticDecayTransportCellSource> {
    StochasticDecayTransportCellSource::bind(state, cells()).unwrap()
}

#[test]
fn g1_movement_adds_no_heat_and_the_recipient_receives_the_later_decay_energy() {
    let initial = thermal(&LAW);
    let original = Arc::new(state([[1, 0], [0, 0]], [[1e6, 0.0], [0.0, 0.0]]));
    let movement_time = original.next_event_time_s().unwrap();
    let mut expected = original.as_ref().clone();
    let moved = expected
        .advance_to(movement_time, budget(), || Ok(()))
        .unwrap();
    // Check the sampled role instead of assuming the large transfer rate wins.
    assert_eq!(moved.movement_events, 1);
    assert_eq!(moved.decay_events, 0);
    assert_eq!(moved.transfer_counts, vec![vec![1, 0], vec![0, 0]]);
    assert_eq!(expected.inventories()[0].populations(), [0, 0]);
    assert_eq!(expected.inventories()[1].populations(), [1, 0]);
    let decay_time = expected.next_event_time_s().unwrap();
    assert!(decay_time > movement_time);

    let heating = StochasticDecayTransportHeating::new(&initial, bind(original), budget()).unwrap();
    let (moved_thermal, moved_heating, audit) = heating
        .advance_to(&initial, movement_time, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(moved_heating.audit().decays, 0.0);
    assert_eq!(
        moved_heating.source(CellId(1)).unwrap().state().as_ref(),
        &expected
    );
    assert!(Arc::ptr_eq(
        moved_heating.source(CellId(1)).unwrap().state(),
        moved_heating.source(CellId(2)).unwrap().state(),
    ));
    assert_eq!(
        moved_heating
            .last_step(CellId(1))
            .unwrap()
            .transport()
            .movement_events,
        1
    );
    for cell in [CellId(1), CellId(2)] {
        assert_eq!(moved_thermal.temperature_k(cell).unwrap(), 300.0);
        assert_eq!(moved_thermal.cells()[&cell].mass_kg, 2.0);
    }

    let decayed = expected
        .advance_to(decay_time, budget(), || Ok(()))
        .unwrap();
    assert_eq!(decayed.movement_events, 0);
    assert_eq!(decayed.decay[0].events, 0);
    assert_eq!(decayed.decay[1].events, 1);
    let (heated, heating, audit) = moved_heating
        .advance_to(&moved_thermal, decay_time, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(
        heating.source(CellId(1)).unwrap().state().as_ref(),
        &expected
    );
    assert_eq!(heating.last_step(CellId(1)).unwrap().decay().events, 0);
    assert_eq!(heating.last_step(CellId(2)).unwrap().decay().events, 1);
    assert_eq!(heating.audit().decays, 1.0);
    assert_eq!(heated.temperature_k(CellId(1)).unwrap(), 300.0);
    assert_eq!(heated.temperature_k(CellId(2)).unwrap(), 301.0);
    assert_eq!(audit.external_heat_j, 8.0);
    assert_eq!(
        heating.audit().released,
        EnergyChannels {
            deposited_j: 8.0,
            radiation_j: 2.0,
            escaped_j: 1.0
        }
    );
}

#[test]
fn g1_stable_nuclei_can_cycle_between_cells_without_nuclear_heating() {
    let initial = thermal(&LAW);
    let original = Arc::new(state([[0, 1], [0, 0]], [[0.0, 2.0], [0.0, 2.0]]));
    let mut expected = original.as_ref().clone();
    let report = expected.advance_to(10.0, budget(), || Ok(())).unwrap();
    assert!(report.transfer_counts[0][1] > 0);
    assert!(report.transfer_counts[1][1] > 0);
    assert_eq!(report.decay_events, 0);
    assert_eq!(
        report.transfer_counts[0][1],
        report.transfer_counts[1][1] + expected.inventories()[1].populations()[1],
    );
    let heating = StochasticDecayTransportHeating::new(&initial, bind(original), budget()).unwrap();
    let (after, heating, audit) = heating
        .advance_to(&initial, 10.0, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(
        heating.source(CellId(1)).unwrap().state().as_ref(),
        &expected
    );
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(heating.audit().decays, 0.0);
    assert_eq!(
        heating.cumulative_released().unwrap(),
        EnergyChannels::default()
    );
    for cell in [CellId(1), CellId(2)] {
        assert_eq!(after.temperature_k(cell).unwrap(), 300.0);
        assert_eq!(after.cells()[&cell].mass_kg, 2.0);
        assert_eq!(after.cells()[&cell].boundary.source_w, 0.0);
    }
}

#[test]
fn g3_partitions_restart_and_cell_order_retain_the_complete_cohort_and_local_heat() {
    let initial = thermal(&LAW);
    let original = Arc::new(state([[4, 1], [2, 0]], [[1.0, 2.0], [0.5, 1.5]]));
    let mut expected = original.as_ref().clone();
    expected.advance_to(8.0, budget(), || Ok(())).unwrap();
    let heating =
        StochasticDecayTransportHeating::new(&initial, bind(Arc::clone(&original)), budget())
            .unwrap();
    let (whole_thermal, whole, _) = heating
        .advance_to(&initial, 8.0, controls(100.0), || Ok(()))
        .unwrap();
    let (mut split_thermal, mut split) = (initial.clone(), heating);
    let mut saved = None;
    for end in [0.137, 0.82, 1.93, 8.0] {
        let (next_thermal, next, _) = split
            .advance_to(&split_thermal, end, controls(0.37), || Ok(()))
            .unwrap();
        split_thermal = next_thermal;
        split = next;
        if end == 0.82 {
            saved = Some((split_thermal.clone(), split.clone()));
        }
    }
    let (saved_thermal, saved) = saved.unwrap();
    let (restarted_thermal, restarted, _) = saved
        .advance_to(&saved_thermal, 8.0, controls(10.0), || Ok(()))
        .unwrap();
    for history in [&whole, &split, &restarted] {
        assert_eq!(
            history.source(CellId(1)).unwrap().state().as_ref(),
            &expected
        );
        assert!(Arc::ptr_eq(
            history.source(CellId(1)).unwrap().state(),
            history.source(CellId(2)).unwrap().state(),
        ));
    }
    for (index, inventory) in expected.inventories().iter().enumerate() {
        let cell = CellId(index as u64 + 1);
        let expected_k = 300.0 + inventory.cumulative_released().deposited_j / 8.0;
        for thermal in [&whole_thermal, &split_thermal, &restarted_thermal] {
            assert!((thermal.temperature_k(cell).unwrap() - expected_k).abs() < 1e-10);
        }
    }

    // Thermal CellId ordering differs from the physical compartment ordering.
    // Nuclear histories stay unchanged; each cell receives its bound local heat.
    let mut swapped_cells = cells();
    swapped_cells.reverse();
    let mut sources = StochasticDecayTransportCellSource::bind(original, swapped_cells).unwrap();
    sources.reverse();
    let swapped = StochasticDecayTransportHeating::new(&initial, sources, budget()).unwrap();
    let (swapped_thermal, swapped, _) = swapped
        .advance_to(&initial, 8.0, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(swapped.source(CellId(1)).unwrap().compartment_index(), 1);
    assert_eq!(swapped.source(CellId(2)).unwrap().compartment_index(), 0);
    assert_eq!(
        swapped.source(CellId(1)).unwrap().state().as_ref(),
        &expected
    );
    assert!(
        (swapped_thermal.temperature_k(CellId(1)).unwrap()
            - whole_thermal.temperature_k(CellId(2)).unwrap())
        .abs()
            < 1e-10
    );
    assert!(
        (swapped_thermal.temperature_k(CellId(2)).unwrap()
            - whole_thermal.temperature_k(CellId(1)).unwrap())
        .abs()
            < 1e-10
    );
}

#[test]
fn g4_incomplete_or_mixed_bindings_and_refused_heat_cannot_publish_partial_transport() {
    let initial = thermal(&LAW);
    let original = Arc::new(state([[1, 0], [0, 0]], [[1e6, 0.0], [0.0, 0.0]]));
    assert!(
        StochasticDecayTransportCellSource::bind(
            Arc::clone(&original),
            vec![(CellId(1), MaterialId(1))],
        )
        .is_err()
    );
    assert!(
        StochasticDecayTransportCellSource::bind(
            Arc::clone(&original),
            vec![(CellId(1), MaterialId(1)); 2],
        )
        .is_err()
    );
    let complete = bind(Arc::clone(&original));
    assert!(
        StochasticDecayTransportHeating::new(&initial, vec![complete[0].clone()], budget())
            .is_err()
    );
    let other_cohort = bind(Arc::new(original.as_ref().clone()));
    assert!(
        StochasticDecayTransportHeating::new(
            &initial,
            vec![complete[0].clone(), other_cohort[1].clone()],
            budget(),
        )
        .is_err()
    );

    let movement_time = original.next_event_time_s().unwrap();
    let mut probe = original.as_ref().clone();
    let movement = probe
        .advance_to(movement_time, budget(), || Ok(()))
        .unwrap();
    assert_eq!(movement.movement_events, 1);
    assert_eq!(movement.decay_events, 0);
    let decay_time = probe.next_event_time_s().unwrap();
    let no_movement = StochasticDecayTransportHeating::new(
        &initial,
        bind(Arc::clone(&original)),
        StochasticDecayTransportBudget {
            maximum_movements: 0,
            ..budget()
        },
    )
    .unwrap();
    assert!(matches!(
        no_movement.advance_to(&initial, movement_time, controls(100.0), || Ok(())),
        Err(FissionThermalError::Decay(DecayError::Budget(_)))
    ));
    assert_eq!(
        no_movement.source(CellId(1)).unwrap().state().as_ref(),
        original.as_ref()
    );

    let initial = thermal(&LIMITED);
    let heating =
        StochasticDecayTransportHeating::new(&initial, bind(Arc::clone(&original)), budget())
            .unwrap();
    // A movement succeeds inside the candidate before the receiving material
    // refuses the later 8 J. Neither cell nor the shared history may be published.
    assert!(
        heating
            .advance_to(&initial, decay_time, controls(100.0), || Ok(()))
            .is_err()
    );
    assert_eq!(initial.time_s(), 0.0);
    for cell in [CellId(1), CellId(2)] {
        assert_eq!(initial.temperature_k(cell).unwrap(), 300.0);
        assert_eq!(
            heating.source(cell).unwrap().state().as_ref(),
            original.as_ref()
        );
        assert_eq!(
            heating
                .source(cell)
                .unwrap()
                .state()
                .cumulative_transfer_counts(),
            &[vec![0, 0], vec![0, 0]]
        );
    }
    assert_eq!(
        heating.cumulative_released().unwrap(),
        EnergyChannels::default()
    );
}
