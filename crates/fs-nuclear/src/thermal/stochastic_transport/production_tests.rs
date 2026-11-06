use super::*;
use std::sync::Arc;

use crate::decay::production::{DecayProductionChange, DecayProductionSchedule};
use crate::decay::stochastic::transport::production::{
    StochasticDecayProductionTransportBudget as Budget,
    StochasticDecayProductionTransportState as State,
};
use crate::decay::stochastic::transport::{
    NuclideOutflow, NuclideTransfer, NuclideTransport, StochasticDecayTransportBudget,
    StochasticDecayTransportLimits, StochasticDecayTransportState,
};
use crate::decay::stochastic::{StochasticDecayBudget, StochasticDecayInventory, StreamKey};
use crate::decay::{DecayBranch, DecayModel, DecayNuclide};
use crate::radiation::{
    RadiationDeposit, RadiationResponse, RadiationResponseLimits, RadiationSourceResponse,
};
use crate::thermal::RadiationReceiver;
use fs_la::decay::DecayLimits;
use fs_qty::Time;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, DepositionError, DepositionSimulation,
    Event, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermochem::decay::DecayLawV1;

type Heating = StochasticDecayProductionTransportHeating;
type CellSource = StochasticDecayProductionTransportCellSource;

// Same manufactured calorimetry as the original stochastic coupling: each
// cell has mass 2 kg and c=4 J/(kg K). No isotope or carrier data is inferred.
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

fn budget() -> Budget {
    Budget {
        transport: StochasticDecayTransportBudget {
            decay: StochasticDecayBudget {
                maximum_events: 1024,
                energy_tolerance_j: 1e-8,
            },
            maximum_movements: 2048,
            maximum_clock_checks: 1_000_000,
        },
        maximum_births: 128,
        maximum_schedule_changes: 8,
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

fn history(rates: [[f64; 3]; 2], outflow_rates: [f64; 3]) -> State {
    let model = Arc::new(
        DecayModel::new(
            "manufactured-produced-transported-decay-heat-v1".into(),
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
    let inventories = (0..2)
        .map(|cell| {
            StochasticDecayInventory::new(
                Arc::clone(&model),
                vec![0; 3],
                StreamKey {
                    seed: 7,
                    kernel: 0x484c46,
                    tile: cell + 1,
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
    let transport = StochasticDecayTransportState::new(
        inventories,
        NuclideTransport {
            transfers,
            outflows: vec![NuclideOutflow {
                compartment: 1,
                rates_per_s: outflow_rates.to_vec(),
            }],
        },
        StreamKey {
            seed: 99,
            kernel: 0x484c46,
            tile: 0,
        },
        StochasticDecayTransportLimits::default(),
    )
    .unwrap();
    let schedule = Arc::new(
        DecayProductionSchedule::new(
            model,
            vec![8.0, 0.0, 0.0],
            vec![DecayProductionChange {
                time_s: 1.0,
                rates_nuclei_s: vec![0.0; 3],
            }],
            1,
        )
        .unwrap(),
    );
    State::poisson(transport, vec![Some(schedule), None]).unwrap()
}

fn cells() -> Vec<(CellId, MaterialId)> {
    vec![(CellId(1), MaterialId(1)), (CellId(2), MaterialId(1))]
}

fn bind(state: Arc<State>) -> Vec<CellSource> {
    CellSource::bind(state, cells()).unwrap()
}

fn produced(state: &State) -> u64 {
    (0..state.inventories().len())
        .map(|cell| {
            state
                .cumulative_produced_nuclei(cell)
                .unwrap_or(&[])
                .iter()
                .sum::<u64>()
        })
        .sum()
}

fn assert_nuclei_balance(state: &State) {
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

fn first_events(state: &State) -> (f64, f64, f64) {
    let mut probe = state.clone();
    let (mut birth, mut movement) = (None, None);
    for _ in 0..512 {
        let end = probe
            .next_transition_time_s()
            .expect("a pending source event");
        let step = probe.advance_to(end, budget(), || Ok(())).unwrap();
        if step.births > 0 {
            birth.get_or_insert(end);
        }
        if step.transport.movement_events > 0 {
            movement.get_or_insert(end);
        }
        if step.transport.decay_events > 0 {
            return (
                birth.expect("a nucleus was supplied before decaying"),
                movement.expect("the sampled fixture includes movement before decay"),
                end,
            );
        }
    }
    panic!("the seeded fixture must reach a decay");
}

#[test]
fn g1_supplied_nuclei_move_without_heat_then_release_energy_in_the_recipient() {
    let original = Arc::new(history([[1e6; 3], [0.0; 3]], [0.0; 3]));
    let (birth, movement, decay) = first_events(&original);
    assert!(birth < movement && movement < decay.next_down());
    let mut expected = original.as_ref().clone();
    let mut material = thermal(&LAW);
    let mut heating = Heating::new(&material, bind(original), budget()).unwrap();

    for end in [birth, movement, decay.next_down(), decay] {
        let step = expected.advance_to(end, budget(), || Ok(())).unwrap();
        let (next_material, next, audit) = heating
            .advance_to(&material, end, controls(100.0), || Ok(()))
            .unwrap();
        assert_eq!(next.source(CellId(1)).unwrap().state().as_ref(), &expected);
        assert!(Arc::ptr_eq(
            next.source(CellId(1)).unwrap().state(),
            next.source(CellId(2)).unwrap().state()
        ));
        assert_eq!(next.audit().produced_nuclei, step.births as f64);
        assert_nuclei_balance(&expected);
        if end == birth {
            assert_eq!(step.births, 1);
            assert_eq!(
                next.last_step(CellId(1)).unwrap().report().produced_nuclei,
                vec![vec![1, 0, 0], vec![0; 3]]
            );
        }
        if end == movement {
            assert!(step.transport.movement_events > 0);
        }
        if end < decay {
            assert_eq!(step.transport.decay_events, 0);
            assert_eq!(audit.external_heat_j, 0.0);
            assert_eq!(next_material.temperature_k(CellId(1)).unwrap(), 300.0);
            assert_eq!(next_material.temperature_k(CellId(2)).unwrap(), 300.0);
        } else {
            assert_eq!(step.transport.decay[0].events, 0);
            assert_eq!(step.transport.decay[1].branch_decays, [1, 0]);
            assert_eq!(next.audit().decays, 1.0);
            assert_eq!(
                next.audit().released,
                EnergyChannels {
                    deposited_j: 8.0,
                    radiation_j: 2.0,
                    escaped_j: 1.0
                }
            );
            assert_eq!(audit.external_heat_j, 8.0);
            assert_eq!(next_material.temperature_k(CellId(1)).unwrap(), 300.0);
            assert_eq!(next_material.temperature_k(CellId(2)).unwrap(), 301.0);
        }
        material = next_material;
        heating = next;
    }
    for cell in [CellId(1), CellId(2)] {
        assert_eq!(material.cells()[&cell].mass_kg, 2.0);
        assert_eq!(material.cells()[&cell].boundary.source_w, 0.0);
    }
}

#[test]
fn g1_shutdown_keeps_decay_and_outflow_on_the_same_conservative_source_history() {
    let initial = thermal(&LAW);
    let original = Arc::new(history([[3.0; 3], [0.0; 3]], [0.25; 3]));
    let heating = Heating::new(&initial, bind(original), budget()).unwrap();
    let (at_shutdown, shutdown, _) = heating
        .advance_to(&initial, 1.0, controls(100.0), || Ok(()))
        .unwrap();
    let state = shutdown.source(CellId(1)).unwrap().state();
    let total_produced = produced(state);
    assert!(total_produced > 0);
    assert_eq!(state.schedule(0).unwrap().rates_at(1.0).unwrap(), [0.0; 3]);
    assert!(state.schedule(1).is_none());
    assert_nuclei_balance(state);
    let (expected, step) = state.trial_to(20.0, budget(), || Ok(())).unwrap();
    assert_eq!(step.births, 0);
    assert_eq!(step.produced_nuclei, vec![vec![0; 3]; 2]);
    assert!(step.transport.decay_events > 0);
    assert!(step.transport.outflow_counts.iter().flatten().sum::<u64>() > 0);
    let deposited: f64 = step
        .transport
        .decay
        .iter()
        .map(|cell| cell.released.deposited_j)
        .sum();

    let (late, continued, audit) = shutdown
        .advance_to(&at_shutdown, 20.0, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(
        continued.source(CellId(1)).unwrap().state().as_ref(),
        &expected
    );
    assert_eq!(continued.audit().produced_nuclei, 0.0);
    assert_eq!(continued.audit().decays, step.transport.decay_events as f64);
    assert!((audit.external_heat_j - deposited).abs() < 1e-8);
    assert_eq!(produced(&expected), total_produced);
    assert_nuclei_balance(&expected);
    for (index, inventory) in expected.inventories().iter().enumerate() {
        let expected_k = 300.0 + inventory.cumulative_released().deposited_j / 8.0;
        assert!((late.temperature_k(CellId(index as u64 + 1)).unwrap() - expected_k).abs() < 1e-10);
    }
}

#[test]
fn g1_transported_source_radiation_credits_its_receiver_once_and_keeps_escape_separate() {
    let initial = thermal(&LAW);
    let original = Arc::new(history([[1e6; 3], [0.0; 3]], [0.0; 3]));
    let (_, _, decay) = first_events(&original);
    let mut expected = original.as_ref().clone();
    let step = expected.advance_to(decay, budget(), || Ok(())).unwrap();
    assert_eq!(step.transport.decay[0].events, 0);
    assert_eq!(step.transport.decay[1].events, 1);
    let response = Arc::new(
        RadiationResponse::new(
            "manufactured-half-energy-return-v1".into(),
            vec![
                RadiationSourceResponse {
                    source_id: 1,
                    deposits: Vec::new(),
                    escape_fraction: 1.0,
                },
                RadiationSourceResponse {
                    source_id: 2,
                    deposits: vec![RadiationDeposit {
                        receiver_id: 1,
                        fraction: 0.5,
                    }],
                    escape_fraction: 0.5,
                },
            ],
            RadiationResponseLimits {
                max_sources: 2,
                max_receivers: 1,
                max_transfers: 1,
            },
        )
        .unwrap(),
    );
    let heating = Heating::new(&initial, bind(original), budget())
        .unwrap()
        .with_radiation_response(
            &initial,
            response,
            vec![RadiationReceiver {
                cell_id: CellId(1),
                material_id: MaterialId(1),
            }],
        )
        .unwrap();
    let (after, heating, audit) = heating
        .advance_to(&initial, decay, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(
        heating.source(CellId(1)).unwrap().state().as_ref(),
        &expected
    );
    let nuclear = heating.audit();
    assert_eq!(nuclear.produced_nuclei, produced(&expected) as f64);
    assert_eq!(
        nuclear.released,
        EnergyChannels {
            deposited_j: 8.0,
            radiation_j: 2.0,
            escaped_j: 1.0
        }
    );
    assert_eq!(nuclear.radiation_deposited_j, 1.0);
    assert_eq!(nuclear.radiation_escaped_j, 1.0);
    assert!((nuclear.thermal_credit_j - 9.0).abs() < 1e-10);
    assert!((audit.external_heat_j - 9.0).abs() < 1e-10);
    assert!((after.temperature_k(CellId(1)).unwrap() - 300.125).abs() < 1e-10);
    assert!((after.temperature_k(CellId(2)).unwrap() - 301.0).abs() < 1e-10);
    assert!(
        (nuclear.thermal_credit_j + nuclear.radiation_escaped_j + nuclear.released.escaped_j
            - nuclear.released.total_j())
        .abs()
            < 1e-10
    );
}

#[test]
fn g3_partitions_restart_and_reversed_cell_ids_keep_birth_movement_and_heat_histories() {
    let initial = thermal(&LAW);
    let original = Arc::new(history([[1.2, 2.0, 1.0], [0.3, 0.4, 0.2]], [0.1; 3]));
    let mut expected = original.as_ref().clone();
    expected.advance_to(8.0, budget(), || Ok(())).unwrap();
    let heating = Heating::new(&initial, bind(Arc::clone(&original)), budget()).unwrap();
    let (whole_material, whole, _) = heating
        .advance_to(&initial, 8.0, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(whole.audit().produced_nuclei, produced(&expected) as f64);
    let (mut split_material, mut split) = (initial.clone(), heating);
    let mut saved = None;
    for end in [0.137, 0.82, 1.0, 1.93, 8.0] {
        let (next_material, next, _) = split
            .advance_to(&split_material, end, controls(0.37), || Ok(()))
            .unwrap();
        split_material = next_material;
        split = next;
        if end == 0.82 {
            saved = Some((split_material.clone(), split.clone()));
        }
    }
    let (saved_material, saved) = saved.unwrap();
    let (restarted_material, restarted, _) = saved
        .advance_to(&saved_material, 8.0, controls(10.0), || Ok(()))
        .unwrap();
    for history in [&whole, &split, &restarted] {
        assert_eq!(
            history.source(CellId(1)).unwrap().state().as_ref(),
            &expected
        );
        assert!(Arc::ptr_eq(
            history.source(CellId(1)).unwrap().state(),
            history.source(CellId(2)).unwrap().state()
        ));
    }
    assert_nuclei_balance(&expected);
    for (index, inventory) in expected.inventories().iter().enumerate() {
        let cell = CellId(index as u64 + 1);
        let expected_k = 300.0 + inventory.cumulative_released().deposited_j / 8.0;
        for material in [&whole_material, &split_material, &restarted_material] {
            assert!((material.temperature_k(cell).unwrap() - expected_k).abs() < 1e-10);
        }
    }

    let mut reverse = cells();
    reverse.reverse();
    let mut sources = CellSource::bind(original, reverse).unwrap();
    sources.reverse();
    let swapped = Heating::new(&initial, sources, budget()).unwrap();
    let (swapped_material, swapped, _) = swapped
        .advance_to(&initial, 8.0, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(swapped.source(CellId(1)).unwrap().compartment_index(), 1);
    assert_eq!(swapped.source(CellId(2)).unwrap().compartment_index(), 0);
    assert_eq!(
        swapped.source(CellId(1)).unwrap().state().as_ref(),
        &expected
    );
    assert_eq!(
        swapped.audit().produced_nuclei,
        whole.audit().produced_nuclei
    );
    assert!(
        (swapped_material.temperature_k(CellId(1)).unwrap()
            - whole_material.temperature_k(CellId(2)).unwrap())
        .abs()
            < 1e-10
    );
    assert!(
        (swapped_material.temperature_k(CellId(2)).unwrap()
            - whole_material.temperature_k(CellId(1)).unwrap())
        .abs()
            < 1e-10
    );
}

#[test]
fn g4_cohort_binding_material_work_and_late_cancellation_refusals_preserve_all_owners() {
    let initial = thermal(&LAW);
    let original = Arc::new(history([[1e6; 3], [0.0; 3]], [0.0; 3]));
    let (birth, movement, decay) = first_events(&original);
    assert!(CellSource::bind(Arc::clone(&original), vec![(CellId(1), MaterialId(1))]).is_err());
    assert!(CellSource::bind(Arc::clone(&original), vec![(CellId(1), MaterialId(1)); 2]).is_err());
    let complete = bind(Arc::clone(&original));
    assert!(Heating::new(&initial, vec![complete[0].clone()], budget()).is_err());
    let other = bind(Arc::new(original.as_ref().clone()));
    assert!(
        Heating::new(
            &initial,
            vec![complete[0].clone(), other[1].clone()],
            budget()
        )
        .is_err()
    );

    let limited = [
        (
            Budget {
                maximum_births: 0,
                ..budget()
            },
            birth,
        ),
        (
            Budget {
                transport: StochasticDecayTransportBudget {
                    maximum_movements: 0,
                    ..budget().transport
                },
                ..budget()
            },
            movement,
        ),
        (
            Budget {
                transport: StochasticDecayTransportBudget {
                    decay: StochasticDecayBudget {
                        maximum_events: 0,
                        ..budget().transport.decay
                    },
                    ..budget().transport
                },
                ..budget()
            },
            decay,
        ),
        (
            Budget {
                maximum_schedule_changes: 0,
                ..budget()
            },
            1.0,
        ),
    ];
    for (budget, end) in limited {
        let heating = Heating::new(&initial, bind(Arc::clone(&original)), budget).unwrap();
        assert!(matches!(
            heating.advance_to(&initial, end, controls(100.0), || Ok(())),
            Err(FissionThermalError::Decay(DecayError::Budget(_)))
        ));
        assert_eq!(
            heating.source(CellId(1)).unwrap().state().as_ref(),
            original.as_ref()
        );
        assert_eq!(produced(heating.source(CellId(1)).unwrap().state()), 0);
        assert_eq!(initial.time_s(), 0.0);
    }

    let limited_material = thermal(&LIMITED);
    let heating = Heating::new(&limited_material, bind(Arc::clone(&original)), budget()).unwrap();
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

    let heating = Heating::new(&initial, bind(Arc::clone(&original)), budget()).unwrap();
    let mut checks = 0;
    let (reference_material, reference, _) = heating
        .advance_to(&initial, decay, controls(100.0), || {
            checks += 1;
            Ok(())
        })
        .unwrap();
    assert!(checks > 1);
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
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(
        heating.source(CellId(1)).unwrap().state().as_ref(),
        original.as_ref()
    );
    let (retried_material, retried, _) = heating
        .advance_to(&initial, decay, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(
        retried.source(CellId(1)).unwrap().state(),
        reference.source(CellId(1)).unwrap().state()
    );
    for cell in [CellId(1), CellId(2)] {
        assert_eq!(initial.temperature_k(cell).unwrap(), 300.0);
        assert_eq!(
            retried_material.temperature_k(cell).unwrap(),
            reference_material.temperature_k(cell).unwrap()
        );
    }
}
