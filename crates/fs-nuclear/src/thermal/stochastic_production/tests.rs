use super::*;
use std::sync::Arc;

use crate::decay::production::{DecayProductionChange, DecayProductionSchedule};
use crate::decay::stochastic::production::{
    StochasticDecayProduction, StochasticDecayProductionBudget,
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

// The original stochastic thermal fixture: a 2 kg cell with c=4 J/(kg K).
// These visible joules/event and source rates are manufactured, not isotope data.
struct Caloric {
    maximum_k: f64,
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        let temperature = energy / 4.0;
        if temperature.is_finite() && temperature > 0.0 && temperature <= self.maximum_k {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.25
    }
}
static LAW: Caloric = Caloric { maximum_k: 1e6 };
static LIMITED: Caloric = Caloric { maximum_k: 300.5 };

fn controls(max_dt_s: f64) -> StepConfig {
    StepConfig {
        max_dt_s,
        energy_tolerance_j: 1e-8,
    }
}

fn budget() -> StochasticDecayProductionBudget {
    StochasticDecayProductionBudget {
        decay: StochasticDecayBudget {
            maximum_events: 1024,
            energy_tolerance_j: 1e-8,
        },
        maximum_births: 256,
        maximum_schedule_changes: 16,
    }
}

fn thermal(law: &'static dyn ThermalLaw, cells: u64, time_s: f64) -> DepositionSimulation<'static> {
    let events = (1..=cells)
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
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 10_000,
        },
    )
    .unwrap()
    .advance_to(time_s, controls(10.0), || Ok(()))
    .unwrap()
    .0
}

fn source(
    rate_nuclei_s: f64,
    shutdown_s: Option<f64>,
    chain: bool,
    time_s: f64,
) -> StochasticDecayProductionCellSource {
    let law = DecayLawV1::from_half_life(Time::new(1.0)).unwrap();
    let mut nodes = vec![
        DecayNuclide {
            id: "parent".into(),
            law,
        },
        DecayNuclide {
            id: "daughter".into(),
            law: if chain { law } else { DecayLawV1::stable() },
        },
    ];
    let mut branches = vec![DecayBranch {
        parent: 0,
        daughter: 1,
        probability: 1.0,
        energy_per_decay: EnergyChannels {
            deposited_j: 8.0,
            radiation_j: 2.0,
            escaped_j: 1.0,
        },
    }];
    if chain {
        nodes.push(DecayNuclide {
            id: "stable".into(),
            law: DecayLawV1::stable(),
        });
        branches.push(DecayBranch {
            parent: 1,
            daughter: 2,
            probability: 1.0,
            energy_per_decay: EnergyChannels {
                deposited_j: 4.0,
                radiation_j: 1.0,
                escaped_j: 0.5,
            },
        });
    }
    let count = nodes.len();
    let model = Arc::new(
        DecayModel::new(
            "manufactured-stochastic-production-heat-v1".into(),
            nodes,
            branches,
            DecayLimits::default(),
        )
        .unwrap(),
    );
    let inventory = StochasticDecayInventory::new_at(
        Arc::clone(&model),
        vec![0; count],
        StreamKey {
            seed: 7,
            kernel: 0x484c46,
            tile: 1,
        },
        time_s,
    )
    .unwrap();
    let mut rates = vec![0.0; count];
    rates[0] = rate_nuclei_s;
    let changes = shutdown_s
        .map(|time_s| DecayProductionChange {
            time_s,
            rates_nuclei_s: vec![0.0; count],
        })
        .into_iter()
        .collect();
    let schedule = Arc::new(DecayProductionSchedule::new(model, rates, changes, 1).unwrap());
    StochasticDecayProductionCellSource {
        cell_id: CellId(1),
        material_id: MaterialId(1),
        production: StochasticDecayProduction::poisson(inventory, schedule).unwrap(),
    }
}

fn first_decay_time(production: &StochasticDecayProduction) -> f64 {
    let mut probe = production.clone();
    for _ in 0..256 {
        let end = probe
            .next_transition_time_s()
            .expect("a pending source event");
        let step = probe.advance_to(end, budget(), || Ok(())).unwrap();
        if step.decay.events > 0 {
            return end;
        }
    }
    panic!("the seeded source must produce a subsequent decay");
}

#[test]
fn g1_births_add_no_heat_and_the_first_sampled_decay_credits_its_original_channels() {
    let initial = thermal(&LAW, 1, 0.0);
    let source = source(8.0, Some(1.0), false, 0.0);
    let birth = source.production.next_transition_time_s().unwrap();
    assert!(birth > 0.0 && birth < 1.0);
    let decay = first_decay_time(&source.production);
    assert!(decay > birth);
    let heating = StochasticDecayProductionHeating::new(&initial, vec![source], budget()).unwrap();

    let (before, heating, audit) = heating
        .advance_to(&initial, birth * 0.5, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(heating.audit().produced_nuclei, 0.0);
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(before.temperature_k(CellId(1)).unwrap(), 300.0);

    let (born, heating, audit) = heating
        .advance_to(&before, birth, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(heating.audit().produced_nuclei, 1.0);
    assert_eq!(heating.audit().decays, 0.0);
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(born.temperature_k(CellId(1)).unwrap(), 300.0);
    assert_eq!(
        heating
            .source(CellId(1))
            .unwrap()
            .production
            .inventory()
            .populations(),
        [1, 0]
    );

    let (waiting, heating, audit) = heating
        .advance_to(&born, decay.next_down(), controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(heating.audit().decays, 0.0);
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(waiting.temperature_k(CellId(1)).unwrap(), 300.0);

    let (heated, heating, audit) = heating
        .advance_to(&waiting, decay, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(heating.audit().decays, 1.0);
    assert_eq!(
        heating.audit().released,
        EnergyChannels {
            deposited_j: 8.0,
            radiation_j: 2.0,
            escaped_j: 1.0,
        }
    );
    assert_eq!(audit.external_heat_j, 8.0);
    assert_eq!(heated.temperature_k(CellId(1)).unwrap(), 301.0);
    assert_eq!(heated.cells()[&CellId(1)].mass_kg, 2.0);
    assert_eq!(heated.cells()[&CellId(1)].boundary.source_w, 0.0);
    let (same, heating, audit) = heating
        .advance_to(&heated, decay, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(heating.cumulative_released().unwrap().deposited_j, 8.0);
    assert_eq!(same.temperature_k(CellId(1)).unwrap(), 301.0);
}

#[test]
fn g1_shutdown_stops_births_while_existing_nuclei_keep_heating_the_material() {
    let initial = thermal(&LAW, 1, 0.0);
    let heating = StochasticDecayProductionHeating::new(
        &initial,
        vec![source(8.0, Some(1.0), false, 0.0)],
        budget(),
    )
    .unwrap();
    // This maximum step crosses shutdown; the source owns its exact boundary.
    let (at_shutdown, shutdown, _) = heating
        .advance_to(&initial, 1.0, controls(100.0), || Ok(()))
        .unwrap();
    let production = &shutdown.source(CellId(1)).unwrap().production;
    let produced = production.cumulative_produced_nuclei().to_vec();
    assert!(produced[0] > 0);
    assert!(production.inventory().populations()[0] > 0);
    assert_eq!(production.schedule().rates_at(1.0).unwrap(), [0.0, 0.0]);
    let next_decay = production.inventory().next_event_time_s().unwrap();
    assert_eq!(production.next_transition_time_s(), Some(next_decay));
    let before_j = shutdown.cumulative_released().unwrap().deposited_j;

    let (later, continued, audit) = shutdown
        .advance_to(&at_shutdown, next_decay, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(continued.audit().produced_nuclei, 0.0);
    assert_eq!(continued.audit().decays, 1.0);
    assert_eq!(audit.external_heat_j, 8.0);
    assert_eq!(
        continued.cumulative_released().unwrap().deposited_j,
        before_j + 8.0
    );
    assert!(
        (later.temperature_k(CellId(1)).unwrap()
            - at_shutdown.temperature_k(CellId(1)).unwrap()
            - 1.0)
            .abs()
            < 1e-12
    );
    assert_eq!(
        continued
            .source(CellId(1))
            .unwrap()
            .production
            .cumulative_produced_nuclei(),
        produced
    );

    let (finished, exhausted, _) = continued
        .advance_to(&later, 100.0, controls(100.0), || Ok(()))
        .unwrap();
    let production = &exhausted.source(CellId(1)).unwrap().production;
    assert_eq!(production.cumulative_produced_nuclei(), produced);
    assert_eq!(production.inventory().populations(), [0, produced[0]]);
    assert_eq!(
        exhausted.cumulative_released().unwrap().deposited_j,
        8.0 * produced[0] as f64
    );
    assert!(
        (finished.temperature_k(CellId(1)).unwrap() - 300.0 - produced[0] as f64).abs() < 1e-10
    );
}

#[test]
fn g3_observation_partitions_and_a_cloned_restart_retain_births_daughters_and_energy() {
    let initial = thermal(&LAW, 1, 0.0);
    let heating = StochasticDecayProductionHeating::new(
        &initial,
        vec![source(8.0, Some(1.0), true, 0.0)],
        budget(),
    )
    .unwrap();
    let (whole_thermal, whole, _) = heating
        .advance_to(&initial, 100.0, controls(100.0), || Ok(()))
        .unwrap();
    let (mut split_thermal, mut split) = (initial, heating);
    let mut saved = None;
    for end in [0.137, 0.82, 1.0, 1.93, 100.0] {
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
        .advance_to(&saved_thermal, 100.0, controls(10.0), || Ok(()))
        .unwrap();
    let expected = &whole.source(CellId(1)).unwrap().production;
    assert_eq!(expected, &split.source(CellId(1)).unwrap().production);
    assert_eq!(expected, &restarted.source(CellId(1)).unwrap().production);
    let born = expected.cumulative_produced_nuclei()[0];
    assert!(born > 0);
    assert_eq!(expected.inventory().populations(), [0, 0, born]);
    assert_eq!(
        whole.cumulative_released().unwrap().deposited_j,
        12.0 * born as f64
    );
    for thermal in [&whole_thermal, &split_thermal, &restarted_thermal] {
        assert!(
            (thermal.temperature_k(CellId(1)).unwrap() - 300.0 - 1.5 * born as f64).abs() < 1e-10
        );
    }
}

#[test]
fn g4_material_work_and_late_cancellation_refusals_preserve_the_complete_sampled_history() {
    let initial = thermal(&LIMITED, 1, 0.0);
    let source = source(8.0, Some(1.0), false, 0.0);
    let original = source.production.clone();
    let birth = original.next_transition_time_s().unwrap();
    let decay = first_decay_time(&original);
    let heating =
        StochasticDecayProductionHeating::new(&initial, vec![source.clone()], budget()).unwrap();
    assert!(
        heating
            .advance_to(&initial, decay, controls(10.0), || Ok(()))
            .is_err()
    );
    assert_eq!(heating.source(CellId(1)).unwrap().production, original);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(1)).unwrap(), 300.0);

    let initial = thermal(&LAW, 1, 0.0);
    let limited = [
        (
            StochasticDecayProductionBudget {
                maximum_births: 0,
                ..budget()
            },
            birth,
        ),
        (
            StochasticDecayProductionBudget {
                decay: StochasticDecayBudget {
                    maximum_events: 0,
                    ..budget().decay
                },
                ..budget()
            },
            decay,
        ),
        (
            StochasticDecayProductionBudget {
                maximum_schedule_changes: 0,
                ..budget()
            },
            1.0,
        ),
    ];
    for (budget, end) in limited {
        let heating =
            StochasticDecayProductionHeating::new(&initial, vec![source.clone()], budget).unwrap();
        assert!(matches!(
            heating.advance_to(&initial, end, controls(10.0), || Ok(())),
            Err(FissionThermalError::Decay(DecayError::Budget(_)))
        ));
        assert_eq!(heating.source(CellId(1)).unwrap().production, original);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.temperature_k(CellId(1)).unwrap(), 300.0);
    }

    let heating = StochasticDecayProductionHeating::new(&initial, vec![source], budget()).unwrap();
    let mut checks = 0;
    let (reference_thermal, reference, _) = heating
        .advance_to(&initial, decay, controls(10.0), || {
            checks += 1;
            Ok(())
        })
        .unwrap();
    assert!(checks > 1);
    let mut seen = 0;
    assert!(matches!(
        heating.advance_to(&initial, decay, controls(10.0), || {
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
    assert_eq!(heating.source(CellId(1)).unwrap().production, original);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(1)).unwrap(), 300.0);
    let (retried_thermal, retried, _) = heating
        .advance_to(&initial, decay, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(
        reference.source(CellId(1)).unwrap().production,
        retried.source(CellId(1)).unwrap().production
    );
    assert_eq!(
        reference_thermal.temperature_k(CellId(1)).unwrap(),
        retried_thermal.temperature_k(CellId(1)).unwrap()
    );
}

#[test]
fn g0_bindings_require_independent_streams_and_preserve_the_original_absolute_endpoint() {
    let initial = thermal(&LAW, 2, 0.0);
    assert!(matches!(
        StochasticDecayProductionHeating::new(
            &initial,
            vec![source(8.0, None, false, 1.0)],
            budget()
        ),
        Err(FissionThermalError::ClockMismatch { .. })
    ));
    let first = source(8.0, Some(1.0), false, 0.0);
    let mut second = first.clone();
    second.cell_id = CellId(2);
    assert!(matches!(
        StochasticDecayProductionHeating::new(&initial, vec![first, second], budget()),
        Err(FissionThermalError::Invalid(_))
    ));

    let start = 0.258_207_182_586_193_1;
    let end = 0.889_967_574_761_249_1;
    assert_ne!(start + (end - start), end);
    let initial = thermal(&LAW, 1, start);
    let heating = StochasticDecayProductionHeating::new(
        &initial,
        vec![source(0.0, None, false, start)],
        budget(),
    )
    .unwrap();
    let (after, heating, _) = heating
        .advance_to(&initial, end, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(after.time_s(), end);
    assert_eq!(
        heating
            .source(CellId(1))
            .unwrap()
            .production
            .inventory()
            .time_s(),
        end
    );
    assert_eq!(heating.audit().produced_nuclei, 0.0);
    assert_eq!(
        heating.cumulative_released().unwrap(),
        EnergyChannels::default()
    );
}
