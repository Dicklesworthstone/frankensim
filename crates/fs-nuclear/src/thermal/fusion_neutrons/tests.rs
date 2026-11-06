//! Manufactured neutron kinetics on the original material thermal test seam.

use super::*;
use crate::diffusion::{
    DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionLink, DiffusionModel,
    DiffusionReaction, DiffusionReactionKind, LeakageBoundary,
};
use crate::fusion_neutrons::{
    FusionNeutronSource, NeutronKineticBudget, NeutronKineticModel, NeutronKineticState,
};
use fs_qty::{Amount, Energy, Time};
use fs_thermochem::fusion::neutrons::FusionNeutronEmissionV1;

fn neutron_budget() -> NeutronKineticBudget {
    NeutronKineticBudget {
        diffusion: DiffusionBudget {
            limits: DiffusionLimits::default(),
            population_tolerance: 1e-8,
        },
        energy_tolerance_j: 1e-8,
    }
}

fn capture_state() -> Arc<NeutronKineticState> {
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-one-cell-kinetic-capture".into(),
                groups: Arc::new(
                    NeutronGroups::new("manufactured-energy-grid".into(), vec![0.0, 8.0]).unwrap(),
                ),
                speeds_m_s: vec![1.0],
                cell_volumes_m3: vec![0.01],
                delayed_groups: vec![],
                reactions: vec![DiffusionReaction {
                    data_id: "manufactured-kinetic-capture".into(),
                    cell: 0,
                    incident_group: 0,
                    macroscopic_cross_section_m_inv: 0.7,
                    kind: DiffusionReactionKind::Capture,
                }],
                links: vec![],
                leakage: vec![],
            },
            neutron_budget().diffusion.limits,
        )
        .unwrap(),
    );
    let model = Arc::new(
        NeutronKineticModel::new(
            "manufactured-kinetic-energy-partition".into(),
            diffusion,
            vec![4.0],
            vec![0.75],
        )
        .unwrap(),
    );
    Arc::new(NeutronKineticState::new(model, vec![0.0]).unwrap())
}

fn spatial_state() -> Arc<NeutronKineticState> {
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-scatter-transport-capture".into(),
                groups: Arc::new(
                    NeutronGroups::new("manufactured-two-energy-grid".into(), vec![0.0, 2.0, 8.0])
                        .unwrap(),
                ),
                speeds_m_s: vec![1.0, 1.0],
                cell_volumes_m3: vec![0.01, 0.01],
                delayed_groups: vec![],
                reactions: vec![
                    DiffusionReaction {
                        data_id: "manufactured-downscatter".into(),
                        cell: 0,
                        incident_group: 1,
                        macroscopic_cross_section_m_inv: 0.5,
                        kind: DiffusionReactionKind::Scatter { outgoing_group: 0 },
                    },
                    DiffusionReaction {
                        data_id: "manufactured-slow-capture".into(),
                        cell: 1,
                        incident_group: 0,
                        macroscopic_cross_section_m_inv: 1.0,
                        kind: DiffusionReactionKind::Capture,
                    },
                    DiffusionReaction {
                        data_id: "manufactured-fast-capture".into(),
                        cell: 1,
                        incident_group: 1,
                        macroscopic_cross_section_m_inv: 1.0,
                        kind: DiffusionReactionKind::Capture,
                    },
                ],
                links: vec![DiffusionLink {
                    cell_a: 0,
                    cell_b: 1,
                    conductance_m2: vec![0.005, 0.005],
                }],
                leakage: vec![LeakageBoundary {
                    cell: 0,
                    conductance_m2: vec![0.001, 0.002],
                }],
            },
            neutron_budget().diffusion.limits,
        )
        .unwrap(),
    );
    let model = Arc::new(
        NeutronKineticModel::new(
            "manufactured-scattering-and-capture-partition".into(),
            diffusion,
            vec![1.0, 4.0],
            vec![1.0, 0.25, 0.25],
        )
        .unwrap(),
    );
    Arc::new(NeutronKineticState::new(model, vec![0.0; 4]).unwrap())
}

fn pulse(
    state: &NeutronKineticState,
    start_s: f64,
    duration_s: f64,
    count: f64,
) -> Arc<FusionNeutronPulse> {
    let amount = Amount::new(count / 6.022_140_76e23);
    let emission = FusionNeutronEmissionV1::new(
        Time::new(duration_s),
        [amount, Amount::new(0.0)],
        [
            Energy::new(4.0 * (amount.value() * 6.022_140_76e23)),
            Energy::new(0.0),
        ],
    )
    .unwrap();
    let source =
        FusionNeutronSource::project_interval_mean(Arc::clone(state.model()), emission, 0).unwrap();
    Arc::new(FusionNeutronPulse::new(Time::new(start_s), source).unwrap())
}

fn one_cell(source_w: f64) -> DepositionSimulation<'static> {
    thermal(
        vec![cell(20, 1, 2.0, 1200.0, source_w)],
        &[(MaterialId(1), &HOT)],
        vec![],
    )
}
fn two_cells() -> DepositionSimulation<'static> {
    thermal(
        vec![cell(20, 1, 2.0, 1200.0, 0.0), cell(10, 2, 3.0, 600.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    )
}
fn bound(
    thermal: &DepositionSimulation<'_>,
    state: Arc<NeutronKineticState>,
    pulse: Arc<FusionNeutronPulse>,
) -> NeutronKineticHeating {
    let cells = if state.model().diffusion().cell_count() == 1 {
        vec![(CellId(20), MaterialId(1))]
    } else {
        vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))]
    };
    NeutronKineticHeating::new(
        thermal,
        NeutronKineticCellSource::bind(state, Time::new(thermal.time_s()), cells, pulse).unwrap(),
        neutron_budget(),
    )
    .unwrap()
}

fn accepted_history(
    state: &NeutronKineticState,
) -> Arc<crate::fusion_neutrons::history::FusionNeutronHistory> {
    use crate::fusion_neutrons::history::{
        FusionNeutronHistory, FusionNeutronHistoryEntry, FusionNeutronHistoryLimits,
    };
    use fs_thermochem::fusion::neutrons::FusionNeutronIntervalV1;
    let intervals = [
        (0, 0.0, 0.25, 4.0, 4.0),
        (0, 0.25, 0.75, 8.0, 2.0),
        (1, 0.125, 0.625, 2.0, 4.0),
    ];
    let entries = intervals
        .into_iter()
        .map(|(cell, start, end, count, mean)| {
            let amount = Amount::new(count / 6.022_140_76e23);
            let emission = FusionNeutronEmissionV1::new(
                Time::new(end - start),
                [amount, Amount::new(0.0)],
                [
                    Energy::new((amount.value() * 6.022_140_76e23) * mean),
                    Energy::new(0.0),
                ],
            )
            .unwrap();
            FusionNeutronHistoryEntry {
                cell,
                interval: FusionNeutronIntervalV1::new(Time::new(start), Time::new(end), emission)
                    .unwrap(),
            }
        })
        .collect();
    Arc::new(
        FusionNeutronHistory::project(
            Arc::clone(state.model()),
            entries,
            FusionNeutronHistoryLimits::default(),
        )
        .unwrap(),
    )
}

#[test]
fn accepted_multicell_history_heats_original_materials_from_the_same_neutron_steps() {
    use crate::fusion_neutrons::history::NeutronKineticHistoryBudget;
    let initial = two_cells();
    let state = spatial_state();
    let history = accepted_history(&state);
    let sources = NeutronKineticCellSource::bind_with_history(
        Arc::clone(&state),
        Time::new(0.0),
        vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))],
        Arc::clone(&history),
    )
    .unwrap();
    let heating = NeutronKineticHeating::new(&initial, sources, neutron_budget()).unwrap();
    let (hot, accepted, _) = heating
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let (direct, report) = state
        .trial_with_history(
            Time::new(0.0),
            Time::new(1.0),
            &history,
            NeutronKineticHistoryBudget {
                transport: neutron_budget(),
                maximum_step_s: 1.0,
                maximum_substeps: 100,
            },
            || Ok(()),
        )
        .unwrap();
    let left = accepted.source(CellId(20)).unwrap();
    let right = accepted.source(CellId(10)).unwrap();
    assert!(Arc::ptr_eq(left.state(), right.state()));
    assert!(Arc::ptr_eq(left.history().unwrap(), &history));
    assert!(left.pulse().is_none());
    for (a, b) in left
        .state()
        .neutrons()
        .neutron_counts()
        .iter()
        .zip(direct.neutrons().neutron_counts())
    {
        close(*a, *b, 1e-9);
    }
    close(
        hot.energy_j() - initial.energy_j(),
        report.released.deposited_j,
        1e-8,
    );
    close(
        accepted.audit().thermal_credit_j,
        report.released.deposited_j,
        1e-8,
    );
    close(
        accepted.cumulative_released().unwrap().radiation_j,
        report.released.radiation_j,
        1e-8,
    );
    assert_eq!(accepted.audit().substeps, report.substeps);
    assert_eq!(
        accepted
            .last_step(CellId(20))
            .unwrap()
            .spatial()
            .diffusion
            .injected_neutrons,
        0.0
    );
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn accepted_history_restart_preserves_material_heat_and_rejects_mixed_cohorts() {
    let initial = two_cells();
    let state = spatial_state();
    let history = accepted_history(&state);
    let bindings = vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))];
    let sources = NeutronKineticCellSource::bind_with_history(
        Arc::clone(&state),
        Time::new(0.0),
        bindings.clone(),
        Arc::clone(&history),
    )
    .unwrap();
    let heating = NeutronKineticHeating::new(&initial, sources.clone(), neutron_budget()).unwrap();
    let (whole, end, _) = heating
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let (middle, running, _) = heating
        .advance_to(&initial, 0.4, controls(1.0), || Ok(()))
        .unwrap();
    let (split, continued, _) = running
        .advance_to(&middle, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    close(whole.energy_j(), split.energy_j(), 1e-8);
    close(
        end.cumulative_released().unwrap().deposited_j,
        continued.cumulative_released().unwrap().deposited_j,
        1e-8,
    );
    let another = accepted_history(&state);
    let other =
        NeutronKineticCellSource::bind_with_history(state, Time::new(0.0), bindings, another)
            .unwrap();
    assert!(
        NeutronKineticHeating::new(
            &initial,
            vec![sources[0].clone(), other[1].clone()],
            neutron_budget()
        )
        .is_err()
    );
}

#[test]
fn g1_finite_fusion_pulse_heats_original_material_only_after_capture() {
    let initial = one_cell(3.0);
    let state = capture_state();
    let pulse = pulse(&state, 0.2, 0.4, 8.0);
    let heating = bound(&initial, Arc::clone(&state), Arc::clone(&pulse));
    let (hot, accepted, _) = heating
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();

    let rate = pulse.source().neutron_rates_s()[0];
    let active = pulse.end().value() - pulse.start().value();
    let remaining =
        rate / 0.7 * -(-0.7 * active).exp_m1() * (-0.7 * (1.0 - pulse.end().value())).exp();
    let captures = rate * active - remaining;
    let released = accepted.cumulative_released().unwrap();
    close(
        accepted
            .source(CellId(20))
            .unwrap()
            .state()
            .neutrons()
            .neutron_counts()[0],
        remaining,
        1e-9,
    );
    close(released.deposited_j, 3.0 * captures, 1e-8);
    close(released.radiation_j, captures, 1e-8);
    assert_eq!(released.escaped_j, 0.0);
    close(
        hot.temperature_k(CellId(20)).unwrap(),
        300.0 + (3.0 * captures + 3.0) / 8.0,
        1e-8,
    );
    close(
        accepted.audit().thermal_credit_j,
        released.deposited_j,
        1e-8,
    );
    close(accepted.audit().captures, captures, 1e-9);
    assert_eq!(accepted.audit().fissions, 0.0);
    assert_eq!(accepted.audit().substeps, 3);
    assert_eq!(accepted.time_s(), hot.time_s());
    let last = accepted.last_step(CellId(20)).unwrap();
    assert_eq!(last.start(), pulse.end());
    assert_eq!(last.spatial().diffusion.injected_neutrons, 0.0);
    assert_eq!(last.source_number_residual(), 0.0);
    assert_eq!(last.source_energy_residual_j(), 0.0);
    assert!(Arc::ptr_eq(last.pulse().unwrap(), &pulse));
    close(
        pulse.source().emission().total_kinetic_energy().value(),
        32.0,
        1e-13,
    );
    assert_eq!(state.neutrons().neutron_counts(), &[0.0]);
    assert_eq!(initial.temperature_k(CellId(20)).unwrap(), 300.0);
}

#[test]
fn g1_spatial_slowing_capture_radiation_leakage_and_heat_share_one_advance() {
    let initial = two_cells();
    let state = spatial_state();
    let pulse = pulse(&state, 0.0, 0.4, 10.0);
    let heating = bound(&initial, Arc::clone(&state), Arc::clone(&pulse));
    let (hot, accepted, _) = heating
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let left = accepted.source(CellId(20)).unwrap();
    let right = accepted.source(CellId(10)).unwrap();
    assert!(Arc::ptr_eq(left.state(), right.state()));
    assert!(Arc::ptr_eq(left.pulse().unwrap(), right.pulse().unwrap()));
    assert!(core::ptr::eq(
        accepted.last_step(CellId(20)).unwrap().spatial(),
        accepted.last_step(CellId(10)).unwrap().spatial(),
    ));
    let released = accepted.cumulative_released().unwrap();
    assert!(released.deposited_j > 0.0 && released.radiation_j > 0.0 && released.escaped_j > 0.0);
    let stored = left.state().kinetic_energy_j().unwrap();
    assert!(stored > 0.0 && stored < 40.0);
    close(stored + released.total_j(), 40.0, 1e-8);
    let left_heat = 8.0 * (hot.temperature_k(CellId(20)).unwrap() - 300.0);
    let right_heat = 6.0 * (hot.temperature_k(CellId(10)).unwrap() - 300.0);
    assert!(left_heat > 0.0 && right_heat > 0.0);
    close(
        left_heat,
        left.state().cell_released(0).unwrap().deposited_j,
        1e-8,
    );
    close(
        right_heat,
        right.state().cell_released(1).unwrap().deposited_j,
        1e-8,
    );
    close(
        left_heat + right_heat,
        accepted.audit().thermal_credit_j,
        1e-8,
    );
    assert_eq!(accepted.audit().substeps, 2);

    // An independent use of the original neutron owner advances once during
    // the source interval and once after shutoff. The thermal cohort must not
    // advance the same field once for every bound material cell.
    let (after_pulse, _) = state
        .trial_with_fusion(0.4, &[pulse.source().clone()], neutron_budget(), || Ok(()))
        .unwrap();
    let (reference, _) = after_pulse
        .trial_with_rates(0.6, &[0.0; 4], neutron_budget(), || Ok(()))
        .unwrap();
    for (&actual, &expected) in left
        .state()
        .neutrons()
        .neutron_counts()
        .iter()
        .zip(reference.neutrons().neutron_counts())
    {
        close(actual, expected, 1e-9);
    }
}

#[test]
fn g3_restart_and_native_subdivision_preserve_source_shutoff_and_energy() {
    let initial = one_cell(0.0);
    let state = capture_state();
    let pulse = pulse(&state, 0.2, 0.4, 8.0);
    let heating = bound(&initial, state, pulse);
    let (whole, whole_source, _) = heating
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let (middle, middle_source, _) = heating
        .advance_to(&initial, 0.45, controls(0.07), || Ok(()))
        .unwrap();
    let (restarted, restarted_source, _) = middle_source
        .advance_to(&middle, 1.0, controls(0.07), || Ok(()))
        .unwrap();
    close(
        restarted.temperature_k(CellId(20)).unwrap(),
        whole.temperature_k(CellId(20)).unwrap(),
        1e-8,
    );
    close(
        restarted_source
            .source(CellId(20))
            .unwrap()
            .state()
            .neutrons()
            .neutron_counts()[0],
        whole_source
            .source(CellId(20))
            .unwrap()
            .state()
            .neutrons()
            .neutron_counts()[0],
        1e-9,
    );
    close(
        restarted_source.cumulative_released().unwrap().total_j(),
        whole_source.cumulative_released().unwrap().total_j(),
        1e-8,
    );
    assert_eq!(restarted_source.time_s(), 1.0);
    assert_eq!(middle_source.time_s(), 0.45);
    assert_eq!(heating.time_s(), 0.0);
}

#[test]
fn g0_cell_model_clock_geometry_and_complete_cohort_binding_are_required() {
    let initial = two_cells();
    let state = spatial_state();
    let pulse = pulse(&state, 0.2, 0.4, 8.0);
    let cells = vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))];
    for invalid in [cells[..1].to_vec(), vec![cells[0], cells[0]]] {
        assert!(
            NeutronKineticCellSource::bind(
                Arc::clone(&state),
                Time::new(0.0),
                invalid,
                Arc::clone(&pulse)
            )
            .is_err()
        );
    }
    assert!(
        NeutronKineticCellSource::bind(
            capture_state(),
            Time::new(0.0),
            vec![cells[0]],
            Arc::clone(&pulse)
        )
        .is_err()
    );
    let stale = NeutronKineticCellSource::bind(
        Arc::clone(&state),
        Time::new(0.1),
        cells.clone(),
        Arc::clone(&pulse),
    )
    .unwrap();
    assert!(NeutronKineticHeating::new(&initial, stale, neutron_budget()).is_err());
    let sources = NeutronKineticCellSource::bind(
        Arc::clone(&state),
        Time::new(0.0),
        cells.clone(),
        Arc::clone(&pulse),
    )
    .unwrap();
    let separate_state = NeutronKineticCellSource::bind(
        Arc::new((*state).clone()),
        Time::new(0.0),
        cells.clone(),
        Arc::clone(&pulse),
    )
    .unwrap();
    assert!(
        NeutronKineticHeating::new(
            &initial,
            vec![sources[0].clone(), separate_state[1].clone()],
            neutron_budget()
        )
        .is_err()
    );
    let separate_pulse = NeutronKineticCellSource::bind(
        Arc::clone(&state),
        Time::new(0.0),
        cells,
        Arc::new((*pulse).clone()),
    )
    .unwrap();
    assert!(
        NeutronKineticHeating::new(
            &initial,
            vec![sources[0].clone(), separate_pulse[1].clone()],
            neutron_budget()
        )
        .is_err()
    );
    let mut wrong_volume = cell(20, 1, 2.0, 1200.0, 0.0);
    wrong_volume.geometry.volume_m3 = 0.02;
    let wrong = thermal(
        vec![wrong_volume, cell(10, 2, 3.0, 600.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    assert!(matches!(
        NeutronKineticHeating::new(&wrong, sources.clone(), neutron_budget()),
        Err(FissionThermalError::BindingChanged(CellId(20)))
    ));
    let refs: Vec<_> = sources.iter().collect();
    assert!(
        NeutronKineticCellSource::trial_sources(
            &refs,
            0.0,
            0.4,
            &[300.0; 2],
            neutron_budget(),
            &mut || Ok(())
        )
        .is_err()
    );
    assert!(
        NeutronKineticCellSource::trial_sources(
            &refs,
            0.1,
            0.0,
            &[300.0; 2],
            neutron_budget(),
            &mut || Ok(())
        )
        .is_err()
    );
    for start in [-1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(FusionNeutronPulse::new(Time::new(start), pulse.source().clone()).is_err());
    }
}

#[test]
fn g4_cancellation_keeps_material_and_complete_neutron_checkpoint_unchanged() {
    let initial = two_cells();
    let state = spatial_state();
    let pulse = pulse(&state, 0.2, 0.4, 8.0);
    let heating = bound(&initial, Arc::clone(&state), pulse);
    let mut polls = 0;
    let result = heating.advance_to(&initial, 1.0, controls(0.1), || {
        polls += 1;
        if polls == 80 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        result,
        Err(FissionThermalError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(heating.time_s(), 0.0);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(20)).unwrap(), 300.0);
    assert_eq!(initial.temperature_k(CellId(10)).unwrap(), 300.0);
    assert_eq!(state.neutrons().neutron_counts(), &[0.0; 4]);
    assert_eq!(heating.cumulative_released().unwrap().total_j(), 0.0);
    assert!(Arc::ptr_eq(
        heating.source(CellId(20)).unwrap().state(),
        &state
    ));
    assert!(
        heating
            .advance_to(&initial, 1.0, controls(0.1), || Ok(()))
            .is_ok()
    );
}

#[test]
fn g0_absolute_pulse_clock_cannot_change_the_original_integrated_source() {
    let state = capture_state();
    let original = pulse(&state, 0.0, 0.75, 3.0);
    // At this epoch, adding 0.75 seconds rounds to a full second. Keeping
    // the original rate would emit 4/3 of the declared neutron inventory.
    let large_epoch = 4_503_599_627_370_496.0;
    assert_eq!((large_epoch + 0.75) - large_epoch, 1.0);
    assert!(FusionNeutronPulse::new(Time::new(large_epoch), original.source().clone(),).is_err());
    // Ordinary binary64 endpoint subtraction remains admitted without any
    // rescaling of the original source rate or interval-integrated accounts.
    let ordinary = pulse(&state, 0.0, 0.4, 3.0);
    let shifted = FusionNeutronPulse::new(Time::new(0.2), ordinary.source().clone()).unwrap();
    assert_ne!(shifted.end().value() - shifted.start().value(), 0.4);
    assert_eq!(
        shifted.source().neutron_rates_s(),
        ordinary.source().neutron_rates_s()
    );
    assert_eq!(shifted.source().emission(), ordinary.source().emission());
}
