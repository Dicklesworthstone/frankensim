use super::*;
use std::sync::Arc;

use crate::decay::stochastic::StreamKey;
use crate::decay::{DecayBranch, DecayModel, DecayNuclide};
use fs_la::decay::DecayLimits;
use fs_qty::Time;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermochem::decay::DecayLawV1;

// Manufactured per-event energies and calorics make the coupling visible.
// They are not isotope data or a calibrated particle deposition model.
struct Caloric {
    maximum_k: f64,
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        let temperature = energy / 4.0;
        if temperature > 0.0 && temperature <= self.maximum_k {
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
fn budget() -> StochasticDecayBudget {
    StochasticDecayBudget {
        maximum_events: 1024,
        energy_tolerance_j: 1e-8,
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
fn source(nuclei: u64, chain: bool, time_s: f64) -> StochasticDecayCellSource {
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
    let mut populations = vec![nuclei, 0];
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
        populations.push(0);
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
    let model = Arc::new(
        DecayModel::new(
            "synthetic-event-calorimetry".into(),
            nodes,
            branches,
            DecayLimits::default(),
        )
        .unwrap(),
    );
    StochasticDecayCellSource {
        cell_id: CellId(1),
        material_id: MaterialId(1),
        inventory: StochasticDecayInventory::new_at(
            model,
            populations,
            StreamKey {
                seed: 7,
                kernel: 0x484c46,
                tile: 1,
            },
            time_s,
        )
        .unwrap(),
    }
}

#[test]
fn g1_no_heat_precedes_the_sampled_event_and_channels_are_credited_once() {
    let thermal = thermal(&LAW, 1, 0.0);
    let source = source(1, false, 0.0);
    let event = source.inventory.next_event_time_s().unwrap();
    let heating = StochasticDecayHeating::new(&thermal, vec![source], budget()).unwrap();
    let (before, heating, audit) = heating
        .advance_to(&thermal, 0.5 * event, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(before.temperature_k(CellId(1)).unwrap(), 300.0);
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(heating.audit().decays, 0.0);
    let (after, heating, audit) = heating
        .advance_to(&before, event, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(after.temperature_k(CellId(1)).unwrap(), 301.0);
    assert_eq!(
        heating.source(CellId(1)).unwrap().inventory.populations(),
        &[0, 1]
    );
    assert_eq!(heating.audit().decays, 1.0);
    assert_eq!(
        heating.audit().released,
        EnergyChannels {
            deposited_j: 8.0,
            radiation_j: 2.0,
            escaped_j: 1.0
        }
    );
    assert_eq!(audit.external_heat_j, 8.0);
    assert_eq!(after.cells()[&CellId(1)].boundary.source_w, 0.0);
    let (late, heating, _) = heating
        .advance_to(&after, event + 1.0, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(late.temperature_k(CellId(1)).unwrap(), 301.0);
    assert_eq!(heating.audit().decays, 0.0);
    assert_eq!(heating.cumulative_released().unwrap().deposited_j, 8.0);
}

#[test]
fn g3_subdivisions_keep_daughter_events_rng_and_event_energy() {
    let thermal = thermal(&LAW, 1, 0.0);
    let source = source(8, true, 0.0);
    let heating = StochasticDecayHeating::new(&thermal, vec![source], budget()).unwrap();
    let (whole_thermal, whole, _) = heating
        .advance_to(&thermal, 100.0, controls(100.0), || Ok(()))
        .unwrap();
    let (mut split_thermal, mut split) = (thermal, heating);
    for endpoint in [0.137, 0.82, 1.93, 100.0] {
        let (next_thermal, next, _) = split
            .advance_to(&split_thermal, endpoint, controls(0.37), || Ok(()))
            .unwrap();
        split_thermal = next_thermal;
        split = next;
    }
    assert_eq!(
        whole.source(CellId(1)).unwrap().inventory,
        split.source(CellId(1)).unwrap().inventory
    );
    assert_eq!(
        whole.source(CellId(1)).unwrap().inventory.populations(),
        &[0, 0, 8]
    );
    assert_eq!(whole.cumulative_released().unwrap().deposited_j, 96.0);
    assert!((whole_thermal.temperature_k(CellId(1)).unwrap() - 312.0).abs() < 1e-12);
    assert!((split_thermal.temperature_k(CellId(1)).unwrap() - 312.0).abs() < 1e-12);
}

#[test]
fn g0_source_receives_the_original_endpoint_without_reconstructing_its_clock() {
    let start = 0.258_207_182_586_193_1;
    let end = 0.889_967_574_761_249_1;
    assert_ne!(start + (end - start), end);
    let thermal = thermal(&LAW, 1, start);
    let source = source(1, false, start);
    let heating = StochasticDecayHeating::new(&thermal, vec![source], budget()).unwrap();
    let (thermal, heating, _) = heating
        .advance_to(&thermal, end, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(thermal.time_s(), end);
    assert_eq!(heating.source(CellId(1)).unwrap().inventory.time_s(), end);
    heating
        .advance_to(&thermal, end + 1.0, controls(10.0), || Ok(()))
        .unwrap();
}

#[test]
fn g4_refusals_preserve_material_inventory_future_clocks_and_rng() {
    let thermal = thermal(&LIMITED, 1, 0.0);
    let source = source(1, false, 0.0);
    let event = source.inventory.next_event_time_s().unwrap();
    let original = source.inventory.clone();
    let heating = StochasticDecayHeating::new(&thermal, vec![source.clone()], budget()).unwrap();
    assert!(
        heating
            .advance_to(&thermal, event, controls(10.0), || Ok(()))
            .is_err()
    );
    assert_eq!(heating.source(CellId(1)).unwrap().inventory, original);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.temperature_k(CellId(1)).unwrap(), 300.0);

    let zero = StochasticDecayHeating::new(
        &thermal,
        vec![source],
        StochasticDecayBudget {
            maximum_events: 0,
            ..budget()
        },
    )
    .unwrap();
    assert!(matches!(
        zero.advance_to(&thermal, event, controls(10.0), || Ok(())),
        Err(FissionThermalError::Decay(DecayError::Budget(
            "stochastic decay events"
        )))
    ));
    assert_eq!(zero.source(CellId(1)).unwrap().inventory, original);
    assert!(matches!(
        zero.advance_to(&thermal, event, controls(10.0), || Err(
            DepositionError::Cancelled
        )),
        Err(FissionThermalError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(zero.source(CellId(1)).unwrap().inventory, original);
}

#[test]
fn g0_cells_require_matching_clocks_and_independent_random_streams() {
    let thermal = thermal(&LAW, 2, 0.0);
    assert!(matches!(
        StochasticDecayHeating::new(&thermal, vec![source(1, false, 1.0)], budget()),
        Err(FissionThermalError::ClockMismatch { .. })
    ));
    let first = source(1, false, 0.0);
    let mut second = first.clone();
    second.cell_id = CellId(2);
    assert!(matches!(
        StochasticDecayHeating::new(&thermal, vec![first, second], budget()),
        Err(FissionThermalError::Invalid(
            "independent stochastic decay cells require distinct stream keys"
        ))
    ));
}
