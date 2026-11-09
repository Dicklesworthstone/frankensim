//! Manufactured coefficients isolate the coupled equations, not isotope data.

use super::*;
use crate::decay::{DecayBranch, DecayBudget, DecayInventory, DecayModel, DecayNuclide};
use crate::diffusion::decay::{
    DecayNeutronBudget, DecayNeutronEmission, DecayNeutronSource, DecayNeutronState,
};
use crate::diffusion::{
    DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionLink, DiffusionModel,
    DiffusionReaction, DiffusionState,
};
use crate::fusion_neutrons::NeutronKineticModel;
use crate::irradiation::NeutronGroups;
use fs_la::decay::DecayLimits;
use fs_qty::Frequency;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, DepositionError, DepositionSimulation,
    Event, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermochem::decay::DecayLawV1;

fn budget() -> DecayNeutronEnergyBudget {
    DecayNeutronEnergyBudget {
        spatial: DecayNeutronBudget {
            diffusion: DiffusionBudget {
                limits: DiffusionLimits::default(),
                population_tolerance: 1e-8,
            },
            decay: DecayBudget {
                propagation: Default::default(),
                nuclei_tolerance: 1e-8,
                energy_tolerance_j: 1e-7,
            },
        },
        energy_tolerance_j: 1e-7,
    }
}
fn state(cells: usize) -> Arc<DecayNeutronEnergyState> {
    let groups = Arc::new(
        NeutronGroups::new("manufactured-decay-heat-grid".into(), vec![0.0, 8.0]).unwrap(),
    );
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-decay-heat-transport".into(),
                groups: Arc::clone(&groups),
                speeds_m_s: vec![1.0],
                cell_volumes_m3: vec![0.01; cells],
                delayed_groups: vec![],
                reactions: vec![DiffusionReaction {
                    data_id: "manufactured-capture".into(),
                    cell: cells - 1,
                    incident_group: 0,
                    macroscopic_cross_section_m_inv: 0.5,
                    kind: DiffusionReactionKind::Capture,
                }],
                links: if cells == 1 {
                    vec![]
                } else {
                    vec![DiffusionLink {
                        cell_a: 0,
                        cell_b: 1,
                        conductance_m2: vec![0.01],
                    }]
                },
                leakage: vec![],
            },
            DiffusionLimits::default(),
        )
        .unwrap(),
    );
    let kinetic = Arc::new(
        NeutronKineticModel::new(
            "manufactured-capture-partition".into(),
            Arc::clone(&diffusion),
            vec![4.0],
            vec![0.75],
        )
        .unwrap(),
    );
    let decay = Arc::new(
        DecayModel::new(
            "manufactured-decay-heat-source".into(),
            vec![
                DecayNuclide {
                    id: "parent".into(),
                    law: DecayLawV1::from_decay_constant(Frequency::new(0.5)).unwrap(),
                },
                DecayNuclide {
                    id: "daughter".into(),
                    law: DecayLawV1::stable(),
                },
            ],
            vec![DecayBranch {
                parent: 0,
                daughter: 1,
                probability: 1.0,
                energy_per_decay: EnergyChannels {
                    deposited_j: 2.0,
                    radiation_j: 4.0,
                    escaped_j: 1.0,
                },
            }],
            DecayLimits::default(),
        )
        .unwrap(),
    );
    let emission = Arc::new(
        DecayNeutronEmission::new(
            "manufactured-one-neutron".into(),
            Arc::clone(&decay),
            groups,
            vec![vec![1.0], vec![0.0]],
        )
        .unwrap(),
    );
    let spatial = DecayNeutronState::new(
        DiffusionState::new(diffusion, vec![0.0; cells], vec![]).unwrap(),
        vec![DecayNeutronSource {
            cell: 0,
            inventory: DecayInventory::new(decay, vec![100.0, 0.0]).unwrap(),
            emission,
        }],
    )
    .unwrap();
    Arc::new(DecayNeutronEnergyState::new(kinetic, spatial).unwrap())
}

struct Caloric {
    maximum_k: f64,
}
impl ThermalLaw for Caloric {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([7; 32])
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = h / 4.0;
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
static LIMITED: Caloric = Caloric { maximum_k: 310.0 };
fn controls(max_dt_s: f64) -> StepConfig {
    StepConfig {
        max_dt_s,
        energy_tolerance_j: 1e-7,
    }
}
fn thermal(
    cells: usize,
    volume: f64,
    law: &'static dyn ThermalLaw,
) -> DepositionSimulation<'static> {
    let events = [20, 10]
        .into_iter()
        .take(cells)
        .map(|id| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id + 100,
                body_id: id + 200,
                material_id: MaterialId(1),
                mass_kg: 2.0,
                specific_h_j_kg: 1200.0,
                geometry: CellGeometry {
                    start_m: [id as f64, 0.0, 0.0],
                    end_m: [id as f64, 1.0, 0.0],
                    volume_m3: volume,
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
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, controls(1.0), || Ok(()))
    .unwrap()
    .0
}
fn bound(
    thermal: &DepositionSimulation<'_>,
    state: Arc<DecayNeutronEnergyState>,
) -> DecayNeutronHeating {
    let cells = [20, 10]
        .into_iter()
        .take(state.model().diffusion().cell_count())
        .map(|id| (CellId(id), MaterialId(1)))
        .collect();
    DecayNeutronHeating::new(
        thermal,
        DecayNeutronCellSource::bind(state, Time::new(thermal.time_s()), cells).unwrap(),
        budget(),
    )
    .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 1e-7, "{actual} != {expected}");
}

#[test]
fn g1_decay_neutron_age_and_capture_heat_match_equal_rate_solution() {
    let initial = thermal(1, 0.01, &LAW);
    let heating = bound(&initial, state(1));
    let (warm, accepted, _) = heating
        .advance_to(&initial, 1.0, controls(0.25), || Ok(()))
        .unwrap();
    let decays = 100.0 * -(-0.5_f64).exp_m1();
    let neutrons = 50.0 * (-0.5_f64).exp();
    let captures = decays - neutrons;
    let heat = 2.0 * decays + 3.0 * captures;
    let after = accepted.source(CellId(20)).unwrap().state();
    close(after.spatial().neutrons().neutron_counts()[0], neutrons);
    close(accepted.audit().decays, decays);
    close(accepted.audit().captures, captures);
    assert_eq!(accepted.audit().fissions, 0.0);
    close(warm.energy_j() - initial.energy_j(), heat);
    close(warm.temperature_k(CellId(20)).unwrap(), 300.0 + heat / 8.0);
    let released = accepted.cumulative_released().unwrap();
    close(released.deposited_j, heat);
    close(released.radiation_j, captures);
    close(released.escaped_j, decays);
    close(released.total_j() + 4.0 * neutrons, 7.0 * decays);
}

#[test]
fn g3_reverse_cell_ids_share_one_transport_and_deposit_in_original_cells() {
    let initial = thermal(2, 0.01, &LAW);
    let original = state(2);
    let (_, expected) = original.trial(1.0, budget(), || Ok(())).unwrap();
    let heating = bound(&initial, original);
    let (warm, accepted, _) = heating
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let source = accepted.source(CellId(20)).unwrap();
    let receiver = accepted.source(CellId(10)).unwrap();
    assert!(Arc::ptr_eq(source.state(), receiver.state()));
    assert!(Arc::ptr_eq(
        &accepted.last_step(CellId(20)).unwrap().spatial,
        &accepted.last_step(CellId(10)).unwrap().spatial
    ));
    assert!(expected.cell_released[1].deposited_j > 0.0);
    for (index, id) in [20, 10].into_iter().enumerate() {
        close(
            warm.temperature_k(CellId(id)).unwrap(),
            300.0 + expected.cell_released[index].deposited_j / 8.0,
        );
    }
    close(
        accepted.cumulative_released().unwrap().total_j(),
        expected.released.total_j(),
    );
}

#[test]
fn g3_restart_keeps_decaying_source_and_resident_neutron_history() {
    let initial = thermal(1, 0.01, &LAW);
    let heating = bound(&initial, state(1));
    let (whole, whole_source, _) = heating
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let (middle, middle_source, _) = heating
        .advance_to(&initial, 0.4, controls(0.4), || Ok(()))
        .unwrap();
    let (split, split_source, _) = middle_source
        .advance_to(&middle, 1.0, controls(0.6), || Ok(()))
        .unwrap();
    close(split.energy_j(), whole.energy_j());
    // Thermal audits describe one requested advance; physical source accounts
    // retain the entire history across the restart.
    close(
        middle_source.audit().decays + split_source.audit().decays,
        whole_source.audit().decays,
    );
    close(
        middle_source.audit().captures + split_source.audit().captures,
        whole_source.audit().captures,
    );
    close(
        split_source.cumulative_released().unwrap().total_j(),
        whole_source.cumulative_released().unwrap().total_j(),
    );
    close(
        split_source
            .source(CellId(20))
            .unwrap()
            .state()
            .spatial()
            .neutrons()
            .neutron_counts()[0],
        whole_source
            .source(CellId(20))
            .unwrap()
            .state()
            .spatial()
            .neutrons()
            .neutron_counts()[0],
    );
}

#[test]
fn g0_incomplete_mixed_volume_and_clock_bindings_are_refused() {
    let original = state(2);
    assert!(
        DecayNeutronCellSource::bind(
            Arc::clone(&original),
            Time::new(0.0),
            vec![(CellId(20), MaterialId(1))]
        )
        .is_err()
    );
    assert!(
        DecayNeutronCellSource::bind(
            Arc::clone(&original),
            Time::new(0.0),
            vec![(CellId(20), MaterialId(1)); 2]
        )
        .is_err()
    );
    let map = vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(1))];
    let sources =
        DecayNeutronCellSource::bind(Arc::clone(&original), Time::new(0.0), map.clone()).unwrap();
    assert!(DecayNeutronHeating::new(&thermal(2, 0.02, &LAW), sources.clone(), budget()).is_err());
    let mut mixed = sources;
    mixed[1].state = state(2);
    assert!(DecayNeutronHeating::new(&thermal(2, 0.01, &LAW), mixed, budget()).is_err());
    let late = DecayNeutronCellSource::bind(original, Time::new(1.0), map).unwrap();
    assert!(DecayNeutronHeating::new(&thermal(2, 0.01, &LAW), late, budget()).is_err());
}

#[test]
fn g4_late_cancellation_and_caloric_refusal_leave_sources_and_material_unchanged() {
    let initial = thermal(1, 0.01, &LAW);
    let original = state(1);
    let heating = bound(&initial, Arc::clone(&original));
    let mut polls = 0;
    let (expected, _, _) = heating
        .advance_to(&initial, 1.0, controls(0.25), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut remaining = polls / 2;
    let cancelled = heating.advance_to(&initial, 1.0, controls(0.25), || {
        if remaining == 0 {
            Err(DepositionError::Cancelled)
        } else {
            remaining -= 1;
            Ok(())
        }
    });
    assert!(matches!(
        cancelled,
        Err(FissionThermalError::Thermal(DepositionError::Cancelled))
    ));
    assert!(Arc::ptr_eq(
        heating.source(CellId(20)).unwrap().state(),
        &original
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(
        heating.cumulative_released().unwrap(),
        EnergyChannels::default()
    );
    let (retry, _, _) = heating
        .advance_to(&initial, 1.0, controls(0.25), || Ok(()))
        .unwrap();
    close(retry.energy_j(), expected.energy_j());
    let limited = thermal(1, 0.01, &LIMITED);
    let limited_source = bound(&limited, Arc::clone(&original));
    assert!(
        limited_source
            .advance_to(&limited, 1.0, controls(0.25), || Ok(()))
            .is_err()
    );
    assert_eq!(limited.time_s(), 0.0);
    assert!(Arc::ptr_eq(
        limited_source.source(CellId(20)).unwrap().state(),
        &original
    ));
}

#[test]
fn g1_g3_delayed_background_fission_and_radioactivity_heat_the_same_original_material() {
    use crate::diffusion::decay::energy::material_induced_decay_fixture;

    let original = Arc::new(material_induced_decay_fixture(0.6));
    let initial = thermal(1, 0.01, &LAW);
    let heating = bound(&initial, Arc::clone(&original));
    let time = 1.25_f64;
    let (warm, accepted, _) = heating
        .advance_to(&initial, time, controls(0.125), || Ok(()))
        .unwrap();
    let removal = 2.6;
    let response = |rate: f64| ((-rate * time).exp() - (-0.25 * time).exp()) / (0.25 - rate);
    let fissions = 40.0 * -(-time).exp_m1();
    let decays = 60.0 * -(-0.5 * time).exp_m1();
    let residence = 10.0 * -(-removal * time).exp_m1() / removal
        + 20.0 / (removal - 1.0) * (-(-time).exp_m1() + (-removal * time).exp_m1() / removal);
    let neutrons = 40.0 * response(1.0)
        + 40.0 / (removal - 1.0) * (response(1.0) - response(removal))
        + 20.0 * response(removal)
        + 15.0 * response(0.5);
    let captures = fissions + 2.0 * residence + 0.5 * decays - neutrons;
    let heat = 2.0 * decays + 2.0 * fissions + 0.5 * captures;
    close(accepted.audit().fissions, fissions);
    close(accepted.audit().decays, decays);
    close(accepted.audit().captures, captures);
    close(warm.energy_j() - initial.energy_j(), heat);
    close(warm.temperature_k(CellId(20)).unwrap(), 300.0 + heat / 8.0);
    assert_eq!(warm.cells()[&CellId(20)].mass_kg, 2.0);
    let shared = accepted.source(CellId(20)).unwrap().state();
    close(
        shared
            .kinetic()
            .cumulative_induced_precursor_energy_j()
            .unwrap(),
        fissions,
    );
    close(
        shared
            .kinetic()
            .cumulative_exported_precursor_energy_j()
            .unwrap(),
        1.2 * residence,
    );
    close(
        shared.kinetic().stored_energy_j().unwrap()
            + shared.cell_released(0).unwrap().total_j()
            + shared
                .kinetic()
                .cumulative_exported_precursor_energy_j()
                .unwrap(),
        180.0 + 6.0 * decays + 7.0 * fissions,
    );
    let (middle, checkpoint, _) = heating
        .advance_to(&initial, 0.5, controls(0.125), || Ok(()))
        .unwrap();
    let rebound = bound(
        &middle,
        Arc::clone(checkpoint.source(CellId(20)).unwrap().state()),
    );
    let (resumed, continued, _) = rebound
        .advance_to(&middle, time, controls(0.125), || Ok(()))
        .unwrap();
    close(resumed.energy_j(), warm.energy_j());
    close(
        checkpoint.audit().fissions + continued.audit().fissions,
        fissions,
    );
    assert_eq!(original.spatial().neutrons().precursor_counts(), &[10.0]);
    assert_eq!(initial.time_s(), 0.0);
}
