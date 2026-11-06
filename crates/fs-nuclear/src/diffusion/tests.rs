//! G1 manufactured transport/reaction laws; G3 partition, reciprocity and
//! conservation checks. No evaluated isotope data or reactor benchmark claim.

use super::*;

fn data(speeds: &[f64], volumes: &[f64]) -> DiffusionData {
    DiffusionData {
        data_id: "manufactured-diffusion".into(),
        groups: Arc::new(
            NeutronGroups::new(
                "manufactured-energy-grid".into(),
                (0..=speeds.len()).map(|n| n as f64).collect(),
            )
            .unwrap(),
        ),
        speeds_m_s: speeds.to_vec(),
        cell_volumes_m3: volumes.to_vec(),
        delayed_groups: vec![],
        reactions: vec![],
        links: vec![],
        leakage: vec![],
    }
}

fn reaction(
    cell: usize,
    group: usize,
    sigma: f64,
    kind: DiffusionReactionKind,
) -> DiffusionReaction {
    DiffusionReaction {
        data_id: format!("manufactured-channel-{cell}-{group}"),
        cell,
        incident_group: group,
        macroscopic_cross_section_m_inv: sigma,
        kind,
    }
}

fn family(rate: f64, spectrum: &[f64]) -> DelayedGroup {
    DelayedGroup {
        data_id: "manufactured-precursor-family".into(),
        decay_constant_s_inv: rate,
        spectrum: spectrum.to_vec(),
    }
}

fn state(
    data: DiffusionData,
    neutron_counts: Vec<f64>,
    precursor_counts: Vec<f64>,
) -> DiffusionState {
    DiffusionState::new(
        Arc::new(DiffusionModel::new(data, DiffusionLimits::default()).unwrap()),
        neutron_counts,
        precursor_counts,
    )
    .unwrap()
}

fn budget() -> DiffusionBudget {
    DiffusionBudget {
        limits: DiffusionLimits::default(),
        population_tolerance: 1.0e-7,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2.0e-10 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn accounts_close(step: &DiffusionStep) {
    assert!(step.population_roundoff <= budget().population_tolerance);
    close(step.neutron_balance_residual, 0.0);
    close(step.precursor_balance_residual, 0.0);
    for &residual in &step.population_balance_residuals {
        close(residual, 0.0);
    }
}

#[test]
fn capture_and_vacuum_leakage_share_the_same_finite_inventory() {
    let mut supplied = data(&[4.0], &[3.0]);
    supplied
        .reactions
        .push(reaction(0, 0, 0.1, DiffusionReactionKind::Capture));
    supplied.leakage.push(LeakageBoundary {
        cell: 0,
        conductance_m2: vec![0.3],
    });
    let initial = state(supplied, vec![120.0], vec![]);
    let (next, step) = initial.trial(1.25, budget(), || Ok(())).unwrap();
    let lost_fraction = -(-1.0_f64).exp_m1();
    close(next.neutron_counts()[0], 120.0 * (-1.0_f64).exp());
    close(step.captures, 60.0 * lost_fraction);
    close(step.leaked_neutrons, 60.0 * lost_fraction);
    close(step.reaction_counts[0], step.captures);
    close(step.leakage_counts[0], step.leaked_neutrons);
    close(step.integrated_fluence_m2[0], 200.0 * lost_fraction);
    close(
        step.mean_fluxes[0].values_m2_s()[0],
        step.integrated_fluence_m2[0] / 1.25,
    );
    close(initial.group_flux(0).unwrap().values_m2_s()[0], 160.0);
    assert_eq!(initial.neutron_counts(), &[120.0]);
    assert_eq!(step.fissions, 0.0);
    accounts_close(&step);
}

#[test]
fn prompt_multiplication_uses_explicit_yields_and_competing_capture() {
    let mut supplied = data(&[2.0], &[5.0]);
    supplied.reactions = vec![
        reaction(
            0,
            0,
            0.3,
            DiffusionReactionKind::Fission {
                prompt_yields: vec![2.5],
                delayed_yields: vec![],
            },
        ),
        reaction(0, 0, 0.1, DiffusionReactionKind::Capture),
    ];
    let initial = state(supplied, vec![25.0], vec![]);
    let (next, step) = initial.trial(0.8, budget(), || Ok(())).unwrap();
    let growth = 0.6 * (2.5 - 1.0) - 0.2;
    let integral = 25.0 * (growth * 0.8_f64).exp_m1() / growth;
    close(next.neutron_counts()[0], 25.0 * (growth * 0.8_f64).exp());
    close(step.fissions, 0.6 * integral);
    close(step.captures, 0.2 * integral);
    close(step.prompt_neutrons, 2.5 * step.fissions);
    close(step.integrated_fluence_m2[0], 0.4 * integral);
    assert!(step.mean_fluxes[0].values_m2_s()[0] > initial.group_flux(0).unwrap().values_m2_s()[0]);
    assert!(step.mean_fluxes[0].values_m2_s()[0] < next.group_flux(0).unwrap().values_m2_s()[0]);
    accounts_close(&step);
}

#[test]
fn downscatter_integrates_neutrons_born_in_a_zero_loss_group() {
    let mut supplied = data(&[1.0, 4.0], &[2.0]);
    supplied.reactions.push(reaction(
        0,
        1,
        0.25,
        DiffusionReactionKind::Scatter { outgoing_group: 0 },
    ));
    let initial = state(supplied, vec![0.0, 10.0], vec![]);
    let (next, step) = initial.trial(2.0, budget(), || Ok(())).unwrap();
    let lost = -(-2.0_f64).exp_m1();
    close(next.neutron_counts()[0], 10.0 * lost);
    close(next.neutron_counts()[1], 10.0 * (-2.0_f64).exp());
    close(step.scatters, 10.0 * lost);
    close(step.integrated_fluence_m2[0], 5.0 * (2.0 - lost));
    close(step.integrated_fluence_m2[1], 20.0 * lost);
    assert!(step.integrated_fluence_m2[0] > 0.0);
    let (half, first) = initial.trial(1.0, budget(), || Ok(())).unwrap();
    let (split, second) = half.trial(1.0, budget(), || Ok(())).unwrap();
    for group in 0..2 {
        close(split.neutron_counts()[group], next.neutron_counts()[group]);
        close(
            first.integrated_fluence_m2[group] + second.integrated_fluence_m2[group],
            step.integrated_fluence_m2[group],
        );
    }
    close(first.scatters + second.scatters, step.scatters);
    accounts_close(&step);
}

fn exchange_data() -> DiffusionData {
    let mut supplied = data(&[3.0], &[2.0, 5.0]);
    supplied.links.push(DiffusionLink {
        cell_a: 0,
        cell_b: 1,
        conductance_m2: vec![0.4],
    });
    supplied
}

#[test]
fn reciprocal_exchange_conserves_counts_and_equilibrates_flux_not_cell_counts() {
    let supplied = exchange_data();
    let initial = state(supplied.clone(), vec![14.0, 0.0], vec![]);
    let duration = 2.0;
    let (next, step) = initial.trial(duration, budget(), || Ok(())).unwrap();
    let tail = (-0.84_f64 * duration).exp();
    let tail_integral = -(-0.84_f64 * duration).exp_m1() / 0.84;
    let integrals = [
        4.0 * duration + 10.0 * tail_integral,
        10.0 * duration - 10.0 * tail_integral,
    ];
    close(next.neutron_counts()[0], 4.0 + 10.0 * tail);
    close(next.neutron_counts()[1], 10.0 - 10.0 * tail);
    close(next.neutron_counts().iter().sum(), 14.0);
    close(step.link_a_to_b_counts[0], 0.6 * integrals[0]);
    close(step.link_b_to_a_counts[0], 0.24 * integrals[1]);
    close(
        step.link_a_to_b_counts[0] - step.link_b_to_a_counts[0],
        next.neutron_counts()[1],
    );
    close(step.integrated_fluence_m2[0], 1.5 * integrals[0]);
    close(step.integrated_fluence_m2[1], 0.6 * integrals[1]);
    let mut reversed = supplied.clone();
    reversed.links[0].cell_a = 1;
    reversed.links[0].cell_b = 0;
    let (reversed, reversed_step) = state(reversed, vec![14.0, 0.0], vec![])
        .trial(duration, budget(), || Ok(()))
        .unwrap();
    for cell in 0..2 {
        close(reversed.neutron_counts()[cell], next.neutron_counts()[cell]);
    }
    close(
        reversed_step.link_b_to_a_counts[0],
        step.link_a_to_b_counts[0],
    );
    let equilibrium = state(supplied, vec![4.0, 10.0], vec![]);
    let (unchanged, circulating) = equilibrium.trial(duration, budget(), || Ok(())).unwrap();
    close(unchanged.neutron_counts()[0], 4.0);
    close(unchanged.neutron_counts()[1], 10.0);
    close(circulating.link_a_to_b_counts[0], 2.4 * duration);
    close(circulating.link_b_to_a_counts[0], 2.4 * duration);
    close(circulating.mean_fluxes[0].values_m2_s()[0], 6.0);
    close(circulating.mean_fluxes[1].values_m2_s()[0], 6.0);
    accounts_close(&step);
    accounts_close(&circulating);
}

#[test]
fn existing_precursors_emit_into_the_declared_groups_and_zero_loss_receivers() {
    let mut supplied = data(&[1.0, 4.0], &[2.0]);
    supplied.delayed_groups.push(family(0.5, &[0.25, 0.75]));
    let initial = state(supplied, vec![0.0, 0.0], vec![8.0]);
    let (next, step) = initial.trial(2.0, budget(), || Ok(())).unwrap();
    let lost = -(-1.0_f64).exp_m1();
    close(next.precursor_counts()[0], 8.0 * (-1.0_f64).exp());
    close(next.neutron_counts()[0], 2.0 * lost);
    close(next.neutron_counts()[1], 6.0 * lost);
    close(step.precursor_decays, 8.0 * lost);
    close(step.delayed_neutrons, step.precursor_decays);
    close(step.delayed_decays[0], step.precursor_decays);
    close(step.integrated_fluence_m2[0], 2.0 - 2.0 * lost);
    close(step.integrated_fluence_m2[1], 24.0 - 24.0 * lost);
    assert_eq!(step.delayed_precursors_produced, 0.0);
    assert_eq!(step.fissions, 0.0);
    accounts_close(&step);
}

#[test]
fn fission_precursors_born_during_the_interval_can_decay_and_react_again() {
    let mut supplied = data(&[1.0], &[1.0]);
    supplied.delayed_groups.push(family(0.5, &[1.0]));
    supplied.reactions.push(reaction(
        0,
        0,
        0.5,
        DiffusionReactionKind::Fission {
            prompt_yields: vec![0.0],
            delayed_yields: vec![1.0],
        },
    ));
    let initial = state(supplied, vec![8.0], vec![0.0]);
    let (next, step) = initial.trial(2.0, budget(), || Ok(())).unwrap();
    let tail = (-2.0_f64).exp();
    let lost = -(-2.0_f64).exp_m1();
    close(next.neutron_counts()[0], 4.0 + 4.0 * tail);
    close(next.precursor_counts()[0], 4.0 - 4.0 * tail);
    close(step.fissions, 4.0 + 2.0 * lost);
    close(step.precursor_decays, 4.0 - 2.0 * lost);
    close(step.delayed_precursors_produced, step.fissions);
    close(step.delayed_neutrons, step.precursor_decays);
    close(step.integrated_fluence_m2[0], 8.0 + 4.0 * lost);
    assert!(step.precursor_decays > 0.0);
    close(next.neutron_counts()[0] + next.precursor_counts()[0], 8.0);
    accounts_close(&step);
}

#[test]
fn same_group_scattering_and_zero_duration_keep_inventory_without_losing_events() {
    let mut supplied = data(&[2.0], &[4.0]);
    supplied.reactions.push(reaction(
        0,
        0,
        3.0,
        DiffusionReactionKind::Scatter { outgoing_group: 0 },
    ));
    let initial = state(supplied, vec![7.0], vec![]);
    let (next, step) = initial.trial(0.25, budget(), || Ok(())).unwrap();
    close(next.neutron_counts()[0], 7.0);
    close(step.scatters, 6.0 * 7.0 * 0.25);
    close(step.integrated_fluence_m2[0], 0.5 * 7.0 * 0.25);
    let (identity, zero) = next.trial(0.0, budget(), || Ok(())).unwrap();
    assert_eq!(identity.neutron_counts(), next.neutron_counts());
    assert_eq!(zero.integrated_fluence_m2, vec![0.0]);
    assert_eq!(zero.reaction_counts, vec![0.0]);
    close(zero.mean_fluxes[0].values_m2_s()[0], 3.5);
    accounts_close(&step);
    accounts_close(&zero);
}

#[test]
fn cached_operator_reuse_and_coefficient_rebinding_preserve_the_correct_inventory() {
    let mut supplied = data(&[2.0], &[4.0]);
    supplied
        .reactions
        .push(reaction(0, 0, 0.25, DiffusionReactionKind::Capture));
    let initial = state(supplied.clone(), vec![10.0], vec![]);
    let (first, _) = initial.trial(0.5, budget(), || Ok(())).unwrap();
    let (second, _) = first.trial(0.5, budget(), || Ok(())).unwrap();
    assert!(Arc::ptr_eq(
        first.prepared.as_ref().unwrap(),
        second.prepared.as_ref().unwrap()
    ));
    close(second.neutron_counts()[0], 10.0 * (-0.5_f64).exp());
    let (different_time, _) = first.trial(0.25, budget(), || Ok(())).unwrap();
    assert!(!Arc::ptr_eq(
        first.prepared.as_ref().unwrap(),
        different_time.prepared.as_ref().unwrap()
    ));
    supplied.reactions[0].macroscopic_cross_section_m_inv = 0.5;
    supplied.data_id = "manufactured-updated-coefficient".into();
    let rebound = first
        .rebind_model(Arc::new(
            DiffusionModel::new(supplied.clone(), DiffusionLimits::default()).unwrap(),
        ))
        .unwrap();
    assert_eq!(rebound.neutron_counts(), first.neutron_counts());
    assert!(rebound.prepared.is_none());
    let (changed, _) = rebound.trial(0.5, budget(), || Ok(())).unwrap();
    close(
        changed.neutron_counts()[0],
        first.neutron_counts()[0] * (-0.5_f64).exp(),
    );
    assert!(changed.neutron_counts()[0] < second.neutron_counts()[0]);
    let unchanged = first.rebind_model(Arc::clone(first.model())).unwrap();
    assert!(Arc::ptr_eq(
        first.prepared.as_ref().unwrap(),
        unchanged.prepared.as_ref().unwrap()
    ));
    for change in 0..3 {
        let mut incompatible = supplied.clone();
        match change {
            0 => incompatible.cell_volumes_m3[0] = 5.0,
            1 => incompatible.speeds_m_s[0] = 3.0,
            _ => {
                incompatible.groups = Arc::new(
                    NeutronGroups::new("different-grid-identity".into(), vec![0.0, 1.0]).unwrap(),
                )
            }
        }
        let incompatible =
            Arc::new(DiffusionModel::new(incompatible, DiffusionLimits::default()).unwrap());
        assert!(matches!(
            first.rebind_model(incompatible),
            Err(DiffusionError::Invalid(_))
        ));
    }
}

#[test]
fn cancellation_at_start_during_preparation_and_before_publication_is_atomic() {
    let initial = state(exchange_data(), vec![14.0, 0.0], vec![]);
    let mut calls = 0;
    initial
        .trial(2.0, budget(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut candidate = initial.clone();
        let mut visited = 0;
        let result = candidate.advance(2.0, budget(), || {
            visited += 1;
            if visited == stop {
                Err(FissionError::Budget("manufactured checkpoint refusal"))
            } else {
                Ok(())
            }
        });
        assert_eq!(
            result.unwrap_err(),
            DiffusionError::Fission(FissionError::Budget("manufactured checkpoint refusal"))
        );
        assert_eq!(candidate.neutron_counts(), initial.neutron_counts());
        assert_eq!(candidate.precursor_counts(), initial.precursor_counts());
        assert!(candidate.prepared.is_none());
    }
    let mut calls = 0;
    let result =
        DiffusionModel::new_with_checkpoint(exchange_data(), DiffusionLimits::default(), || {
            calls += 1;
            if calls == 2 {
                Err(FissionError::Cancelled)
            } else {
                Ok(())
            }
        });
    assert_eq!(
        result.unwrap_err(),
        DiffusionError::Fission(FissionError::Cancelled)
    );
}

#[test]
fn budgets_invalid_coefficients_and_incompatible_precursors_refuse_before_acceptance() {
    let supplied = exchange_data();
    for limited in [
        DiffusionLimits {
            max_cells: 1,
            ..DiffusionLimits::default()
        },
        DiffusionLimits {
            max_neutron_groups: 0,
            ..DiffusionLimits::default()
        },
        DiffusionLimits {
            max_links: 0,
            ..DiffusionLimits::default()
        },
        DiffusionLimits {
            max_coefficients: 1,
            ..DiffusionLimits::default()
        },
        DiffusionLimits {
            propagation: DecayLimits {
                maximum_nodes: 1,
                ..DecayLimits::default()
            },
            ..DiffusionLimits::default()
        },
        DiffusionLimits {
            propagation: DecayLimits {
                maximum_edges: 1,
                ..DecayLimits::default()
            },
            ..DiffusionLimits::default()
        },
    ] {
        assert!(matches!(
            DiffusionModel::new(supplied.clone(), limited),
            Err(DiffusionError::Budget(_))
        ));
    }
    let mut malformed = supplied.clone();
    malformed.links[0].cell_b = 0;
    assert!(matches!(
        DiffusionModel::new(malformed, DiffusionLimits::default()),
        Err(DiffusionError::Invalid(_))
    ));
    let mut malformed = supplied.clone();
    malformed.cell_volumes_m3[0] = 0.0;
    assert!(matches!(
        DiffusionModel::new(malformed, DiffusionLimits::default()),
        Err(DiffusionError::Invalid(_))
    ));
    let mut malformed = supplied.clone();
    malformed
        .reactions
        .push(reaction(0, 0, -1.0, DiffusionReactionKind::Capture));
    assert!(matches!(
        DiffusionModel::new(malformed, DiffusionLimits::default()),
        Err(DiffusionError::Invalid(_))
    ));
    let mut malformed = supplied;
    malformed.delayed_groups.push(family(0.5, &[0.9]));
    assert!(matches!(
        DiffusionModel::new(malformed, DiffusionLimits::default()),
        Err(DiffusionError::Invalid(_))
    ));
    let mut initial = state(exchange_data(), vec![14.0, 0.0], vec![]);
    let mut limited = budget();
    limited.limits.propagation.maximum_work = 1;
    assert!(initial.advance(2.0, limited, || Ok(())).is_err());
    assert_eq!(initial.neutron_counts(), &[14.0, 0.0]);
    assert!(initial.group_flux(2).is_err());
    assert!(initial.trial(-1.0, budget(), || Ok(())).is_err());
    assert!(DiffusionState::new(Arc::clone(initial.model()), vec![14.0], vec![]).is_err());
    assert!(DiffusionState::new(Arc::clone(initial.model()), vec![14.0, -1.0], vec![]).is_err());
    let mut delayed = data(&[1.0], &[1.0]);
    delayed.delayed_groups.push(family(0.5, &[1.0]));
    let original = state(delayed.clone(), vec![1.0], vec![2.0]);
    delayed.delayed_groups[0].data_id = "different-family".into();
    assert!(
        original
            .rebind_model(Arc::new(
                DiffusionModel::new(delayed, DiffusionLimits::default()).unwrap()
            ))
            .is_err()
    );
}

#[test]
fn unrepresentable_positive_rates_and_competing_fractions_are_not_dropped() {
    let mut tiny_rate = data(&[0.5], &[1.0]);
    tiny_rate.reactions.push(reaction(
        0,
        0,
        f64::from_bits(1),
        DiffusionReactionKind::Capture,
    ));
    assert!(matches!(
        DiffusionModel::new(tiny_rate, DiffusionLimits::default()),
        Err(DiffusionError::Arithmetic(_))
    ));
    let tiny_flux = data(&[f64::MIN_POSITIVE], &[f64::MAX]);
    assert!(matches!(
        DiffusionModel::new(tiny_flux, DiffusionLimits::default()),
        Err(DiffusionError::Arithmetic(_))
    ));
    let mut tiny_fraction = data(&[1.0], &[1.0]);
    tiny_fraction.reactions = vec![
        reaction(0, 0, 1.0e300, DiffusionReactionKind::Capture),
        reaction(
            0,
            0,
            1.0e-300,
            DiffusionReactionKind::Fission {
                prompt_yields: vec![1.0e300],
                delayed_yields: vec![],
            },
        ),
    ];
    assert!(matches!(
        DiffusionModel::new(tiny_fraction, DiffusionLimits::default()),
        Err(DiffusionError::Arithmetic(_))
    ));
    let mut overflow = data(&[2.0], &[1.0]);
    overflow
        .reactions
        .push(reaction(0, 0, f64::MAX, DiffusionReactionKind::Capture));
    assert!(matches!(
        DiffusionModel::new(overflow, DiffusionLimits::default()),
        Err(DiffusionError::Arithmetic(_))
    ));
}

#[test]
fn exported_defect_bounds_round_up_before_comparison_and_refuse_overflow() {
    let allowance = 2.0_f64.powi(60);
    let mut defect = Accumulator::new(allowance);
    defect.add(1.0).unwrap();
    assert_eq!(defect.value(), allowance);
    let bound = upward_account_bound(defect).unwrap();
    assert_eq!(bound, allowance.next_up());
    assert!(bound > allowance);
    assert_eq!(
        upward_account_bound(Accumulator::new(allowance)).unwrap(),
        allowance
    );
    assert_eq!(
        upward_account_bound(Accumulator {
            hi: allowance,
            lo: -1.0
        })
        .unwrap(),
        allowance
    );
    let mut overflow = Accumulator::new(f64::MAX);
    overflow.add(1.0).unwrap();
    assert!(matches!(
        upward_account_bound(overflow),
        Err(DiffusionError::Arithmetic(_))
    ));
}
