//! G1/G3/G4: evolving spectra, original calorics, supplied nuclei and joint rollback.
//! Coefficients are manufactured; the ENDF parser and physical owners are real.

use super::*;
use crate::decay::production::DecayProductionSchedule;
use crate::decay::thermal::DecayCellSource;
use crate::decay::{DecayBudget, DecayInventory};
use crate::endf::decay::model::emissions::{DecayEmissionLimits, PreparedDecayEmission};
use crate::endf::decay::spectra::DecayRadiation;
use crate::radiation::RadiationDeposit;
use crate::radiation::decay::{DecayRadiationClosure, PreparedDecayRadiation};
use crate::radiation::grouped::{
    GroupedRadiationResponse, GroupedRadiationResponseLimits, GroupedRadiationSourceResponse,
    RadiationGroupResponse,
};
use fs_la::decay::DecayLimits;

const EV_J: f64 = 1.602_176_634e-19;
const Q: f64 = 0.5e6 * EV_J;
const INITIAL: f64 = 1e14;

fn partition() -> Arc<PreparedDecayRadiation> {
    let assembly = crate::radiation::decay::tests::gamma_chain(1e6, 2e6, 0.5, 0.25);
    Arc::new(
        PreparedDecayRadiation::new(
            Arc::new(
                PreparedDecayEmission::new(
                    assembly,
                    DecayRadiation::Gamma,
                    vec![0.0, 1.5e6 * EV_J, 3e6 * EV_J],
                    DecayEmissionLimits::default(),
                    || Ok(()),
                )
                .unwrap(),
            ),
            DecayRadiationClosure {
                absolute_j_per_decay: 0.0,
                relative: 1e-12,
            },
            DecayEmissionLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn spectral_budget() -> DecayBudget {
    DecayBudget {
        propagation: DecayLimits::default(),
        nuclei_tolerance: 1e5,
        energy_tolerance_j: 1e-8,
    }
}

fn grouped_response(partition: &PreparedDecayRadiation) -> Arc<GroupedRadiationResponse> {
    let group = |local, remote| RadiationGroupResponse {
        deposits: vec![
            RadiationDeposit {
                receiver_id: 1,
                fraction: local,
            },
            RadiationDeposit {
                receiver_id: 2,
                fraction: remote,
            },
        ],
        escape_fraction: 0.25,
    };
    Arc::new(
        GroupedRadiationResponse::new(
            "manufactured energy-dependent deposited fractions".into(),
            DecayRadiation::Gamma,
            partition.group_edges_j().to_vec(),
            vec![GroupedRadiationSourceResponse {
                source_id: 1,
                groups: vec![group(0.25, 0.5), group(0.5, 0.25)],
            }],
            GroupedRadiationResponseLimits {
                max_sources: 1,
                max_groups: 2,
                max_receivers: 2,
                max_source_groups: 2,
                max_transfers: 4,
            },
        )
        .unwrap(),
    )
}

fn receivers() -> Vec<RadiationReceiver> {
    vec![
        RadiationReceiver {
            cell_id: CellId(1),
            material_id: MaterialId(1),
        },
        RadiationReceiver {
            cell_id: CellId(2),
            material_id: MaterialId(2),
        },
    ]
}

fn initial_thermal() -> DepositionSimulation<'static> {
    thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0), cell(2, 2, 3.0, 600.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    )
}

fn heating(
    thermal: &DepositionSimulation<'_>,
    partition: Arc<PreparedDecayRadiation>,
    initial: f64,
) -> NuclearHeating<DecayCellSource> {
    NuclearHeating::new(
        thermal,
        vec![DecayCellSource {
            cell_id: CellId(1),
            material_id: MaterialId(1),
            inventory: DecayInventory::new(Arc::clone(partition.model()), vec![initial, 0.0, 0.0])
                .unwrap(),
        }],
        spectral_budget(),
    )
    .unwrap()
    .with_grouped_decay_radiation(
        thermal,
        grouped_response(&partition),
        vec![DecayRadiationEmitter {
            cell_id: CellId(1),
            partition,
        }],
        receivers(),
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn g1_daughter_spectral_change_heats_original_distinct_materials_once() {
    let initial = initial_thermal();
    let source = heating(&initial, partition(), INITIAL);
    let (heated, evolved, audit) = source
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let parent = 0.5 * INITIAL;
    let daughter = INITIAL * (1.0 - 0.5 * (1.0 + std::f64::consts::LN_2));
    let local = Q * (parent + 3.0 * daughter);
    let returned = Q * (0.25 * parent + 0.5 * daughter);
    let remote = Q * (0.5 * parent + 0.25 * daughter);
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        300.0 + (local + returned) / 8.0,
        2e-12,
    );
    close(
        heated.temperature_k(CellId(2)).unwrap(),
        300.0 + remote / 6.0,
        2e-12,
    );
    let account = evolved.audit();
    close(account.released.radiation_j, Q * (parent + daughter), 1e-12);
    close(account.radiation_deposited_j, returned + remote, 1e-12);
    close(
        account.radiation_escaped_j,
        0.25 * Q * (parent + daughter),
        1e-12,
    );
    close(account.thermal_credit_j, local + returned + remote, 1e-12);
    close(audit.external_heat_j, account.thermal_credit_j, 1e-12);
    close(account.radiation_ungrouped_j, 0.0, 1e-13);
    assert_eq!(account.radiation_reconciled_spectral_j, 0.0);
    assert!(evolved.radiation_response().is_none());
    assert!(evolved.grouped_radiation_response().is_some());
    // The hard daughter line has a different remote fraction. Reusing the
    // parent's 50% whole-spectrum response would measurably overheat receiver 2.
    assert!((remote - 0.5 * account.released.radiation_j).abs() > 0.1);
    assert_eq!(heated.cells()[&CellId(1)].mass_kg, 2.0);
    assert_eq!(heated.cells()[&CellId(2)].material_id, MaterialId(2));
}

#[test]
fn g3_partitioned_histories_preserve_spectral_heat_and_source_populations() {
    let initial = initial_thermal();
    let source = heating(&initial, partition(), INITIAL);
    let (whole_t, whole_s, _) = source
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let (mid_t, mid_s, _) = source
        .advance_to(&initial, 0.4, controls(1.0), || Ok(()))
        .unwrap();
    let (split_t, split_s, _) = mid_s
        .advance_to(&mid_t, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    for id in [CellId(1), CellId(2)] {
        close(
            whole_t.temperature_k(id).unwrap(),
            split_t.temperature_k(id).unwrap(),
            2e-12,
        );
    }
    for (&a, &b) in whole_s
        .source(CellId(1))
        .unwrap()
        .inventory
        .populations()
        .iter()
        .zip(split_s.source(CellId(1)).unwrap().inventory.populations())
    {
        close(a, b, 32.0);
    }
    close(
        whole_s.audit().radiation_deposited_j,
        mid_s.audit().radiation_deposited_j + split_s.audit().radiation_deposited_j,
        1e-12,
    );
}

#[test]
fn g1_continuously_born_daughters_emit_before_the_interval_endpoint() {
    let initial = initial_thermal();
    let partition = partition();
    let model = Arc::clone(partition.model());
    let schedule = Arc::new(
        DecayProductionSchedule::new(Arc::clone(&model), vec![INITIAL, 0.0, 0.0], vec![], 0)
            .unwrap(),
    );
    let source = DecayProductionCellSource::new(
        CellId(1),
        MaterialId(1),
        DecayInventory::new(model, vec![0.0; 3]).unwrap(),
        schedule,
    )
    .unwrap();
    let source = DecayProductionHeating::new(&initial, vec![source], spectral_budget())
        .unwrap()
        .with_grouped_decay_radiation(
            &initial,
            grouped_response(&partition),
            vec![DecayRadiationEmitter {
                cell_id: CellId(1),
                partition,
            }],
            receivers(),
            || Ok(()),
        )
        .unwrap();
    let (heated, evolved, _) = source
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let lambda = std::f64::consts::LN_2;
    let parent = INITIAL * (1.0 - 0.5 / lambda);
    let daughter = INITIAL * (1.5 - 1.0 / lambda);
    let remote = Q * (0.5 * parent + 0.25 * daughter);
    close(
        heated.temperature_k(CellId(2)).unwrap(),
        300.0 + remote / 6.0,
        2e-12,
    );
    close(
        evolved.audit().released.radiation_j,
        Q * (parent + daughter),
        1e-12,
    );
    assert!(evolved.last_step(CellId(1)).unwrap().decay.decays[1] > 0.0);
    close(evolved.audit().produced_nuclei, INITIAL, 1.0);
}

#[test]
fn g4_late_cancellation_and_receiver_domain_failure_publish_no_history() {
    let initial = initial_thermal();
    let source = heating(&initial, partition(), INITIAL);
    let original_cells = initial.cells().clone();
    let original_nuclei = source
        .source(CellId(1))
        .unwrap()
        .inventory
        .populations()
        .to_vec();
    let mut visits = 0;
    source
        .advance_to(&initial, 1.0, controls(1.0), || {
            visits += 1;
            Ok(())
        })
        .unwrap();
    assert!(visits > 20);
    let mut late = 0;
    let failure = source
        .advance_to(&initial, 1.0, controls(1.0), || {
            late += 1;
            if late == visits - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .err()
        .expect("late cancellation must refuse the joint step");
    assert_eq!(
        failure,
        FissionThermalError::Thermal(DepositionError::Cancelled)
    );
    assert_eq!(initial.cells(), &original_cells);
    assert_eq!(
        source.source(CellId(1)).unwrap().inventory.populations(),
        original_nuclei
    );
    assert_eq!(source.time_s(), 0.0);
    assert!(source.last_step(CellId(1)).is_none());

    let limited = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0), cell(2, 2, 1.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &LIMITED)],
        vec![],
    );
    let source = heating(&limited, partition(), 100.0 * INITIAL);
    assert!(
        source
            .advance_to(&limited, 1.0, controls(1.0), || Ok(()))
            .is_err()
    );
    assert_eq!(
        source.source(CellId(1)).unwrap().inventory.populations()[0],
        100.0 * INITIAL
    );
    assert_eq!(limited.temperature_k(CellId(2)).unwrap(), 300.0);
}

#[test]
fn g0_binding_rejects_relabelled_models_and_mismatched_response_grids() {
    let initial = initial_thermal();
    let first = partition();
    let other = partition();
    let plain = NuclearHeating::new(
        &initial,
        vec![DecayCellSource {
            cell_id: CellId(1),
            material_id: MaterialId(1),
            inventory: DecayInventory::new(Arc::clone(first.model()), vec![INITIAL, 0.0, 0.0])
                .unwrap(),
        }],
        spectral_budget(),
    )
    .unwrap();
    assert!(
        plain
            .clone()
            .with_grouped_decay_radiation(
                &initial,
                grouped_response(&first),
                vec![DecayRadiationEmitter {
                    cell_id: CellId(1),
                    partition: other
                }],
                receivers(),
                || Ok(())
            )
            .is_err()
    );
    let different = GroupedRadiationResponse::new(
        "mismatched group meaning".into(),
        DecayRadiation::Gamma,
        vec![0.0, 3e6 * EV_J],
        vec![GroupedRadiationSourceResponse {
            source_id: 1,
            groups: vec![RadiationGroupResponse {
                deposits: vec![],
                escape_fraction: 1.0,
            }],
        }],
        GroupedRadiationResponseLimits {
            max_sources: 1,
            max_groups: 1,
            max_receivers: 0,
            max_source_groups: 1,
            max_transfers: 0,
        },
    )
    .unwrap();
    assert!(
        plain
            .with_grouped_decay_radiation(
                &initial,
                Arc::new(different),
                vec![DecayRadiationEmitter {
                    cell_id: CellId(1),
                    partition: first
                }],
                vec![],
                || Ok(())
            )
            .is_err()
    );
}

#[test]
fn g3_sampled_decay_spectra_preserve_event_clocks_rng_and_integer_histories() {
    use crate::decay::stochastic::{StochasticDecayBudget, StochasticDecayInventory, StreamKey};

    let initial = thermal(
        vec![cell(1, 1, 1e-9, 1200.0, 0.0), cell(2, 2, 1e-9, 600.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    let partition = partition();
    let state = StochasticDecayInventory::new(
        Arc::clone(partition.model()),
        vec![40, 0, 0],
        StreamKey {
            seed: 17,
            kernel: 0x53504543,
            tile: 1,
        },
    )
    .unwrap();
    let budget = StochasticDecayBudget {
        maximum_events: 128,
        energy_tolerance_j: 1e-20,
    };
    let (reference, _) = state.trial_to(3.0, budget, || Ok(())).unwrap();
    let source = StochasticDecayHeating::new(
        &initial,
        vec![StochasticDecayCellSource {
            cell_id: CellId(1),
            material_id: MaterialId(1),
            inventory: state,
        }],
        budget,
    )
    .unwrap()
    .with_grouped_decay_radiation(
        &initial,
        grouped_response(&partition),
        vec![DecayRadiationEmitter {
            cell_id: CellId(1),
            partition,
        }],
        receivers(),
        || Ok(()),
    )
    .unwrap();
    let (heated, evolved, _) = source
        .advance_to(&initial, 3.0, controls(0.3), || Ok(()))
        .unwrap();
    assert_eq!(evolved.source(CellId(1)).unwrap().inventory, reference);
    let events = reference.cumulative_decays();
    let remote = Q * (0.5 * events[0] as f64 + 0.25 * events[1] as f64);
    close(
        heated.temperature_k(CellId(2)).unwrap(),
        300.0 + remote / 2e-9,
        3e-11,
    );
    close(
        evolved.audit().released.radiation_j,
        reference.cumulative_released().radiation_j,
        1e-24,
    );
    close(
        evolved.audit().radiation_deposited_j,
        0.75 * reference.cumulative_released().radiation_j,
        1e-24,
    );

    let (middle_t, middle_s, _) = source
        .advance_to(&initial, 0.7, controls(0.3), || Ok(()))
        .unwrap();
    let (split_t, split_s, _) = middle_s
        .advance_to(&middle_t, 3.0, controls(0.3), || Ok(()))
        .unwrap();
    assert_eq!(split_s.source(CellId(1)).unwrap().inventory, reference);
    close(
        split_t.temperature_k(CellId(2)).unwrap(),
        heated.temperature_k(CellId(2)).unwrap(),
        3e-11,
    );
}
