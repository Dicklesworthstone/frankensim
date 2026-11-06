//! G1 prescribed-flow/decay laws; G3 count conservation and interval partition.

use std::sync::Arc;

use crate::FissionError;
use crate::diffusion::*;
use crate::irradiation::NeutronGroups;

fn data(volumes: &[f64], decay_rate: f64) -> DiffusionData {
    DiffusionData {
        data_id: "manufactured-precursor-flow".into(),
        groups: Arc::new(
            NeutronGroups::new("manufactured-flow-grid".into(), vec![0.0, 1.0]).unwrap(),
        ),
        speeds_m_s: vec![1.0],
        cell_volumes_m3: volumes.to_vec(),
        delayed_groups: vec![DelayedGroup {
            data_id: "manufactured-flow-family".into(),
            decay_constant_s_inv: decay_rate,
            spectrum: vec![1.0],
        }],
        reactions: vec![],
        links: vec![],
        leakage: vec![],
    }
}

fn flow(from_cell: usize, to_cell: usize, rate: f64) -> PrecursorTransfer {
    PrecursorTransfer {
        from_cell,
        to_cell,
        volumetric_flow_m3_s: rate,
    }
}

fn state(
    data: DiffusionData,
    transport: PrecursorTransport,
    precursors: Vec<f64>,
) -> DiffusionState {
    let cells = data.cell_volumes_m3.len();
    DiffusionState::new(
        Arc::new(
            DiffusionModel::new_with_transport(data, transport, DiffusionLimits::default())
                .unwrap(),
        ),
        vec![0.0; cells],
        precursors,
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

fn survival_integral(rate: f64, duration: f64) -> f64 {
    -(-rate * duration).exp_m1() / rate
}

fn balances(step: &DiffusionStep) {
    close(step.neutron_balance_residual, 0.0);
    close(step.precursor_balance_residual, 0.0);
    for &residual in &step.population_balance_residuals {
        close(residual, 0.0);
    }
    assert!(step.population_roundoff <= budget().population_tolerance);
}

#[test]
fn transported_precursors_decay_in_their_current_cell_with_the_correct_birth_age_fluence() {
    let transport = PrecursorTransport {
        transfers: vec![flow(0, 1, 1.0)],
        outflows: vec![],
    };
    let initial = state(data(&[2.0, 5.0], 0.5), transport, vec![10.0, 0.0]);
    let duration = 2.0;
    let (next, step) = initial.trial(duration, budget(), || Ok(())).unwrap();
    let fast = survival_integral(1.0, duration);
    let slow = survival_integral(0.5, duration);
    close(next.precursor_counts()[0], 10.0 * (-duration).exp());
    close(
        next.precursor_counts()[1],
        10.0 * ((-0.5 * duration).exp() - (-duration).exp()),
    );
    close(next.neutron_counts()[0], 5.0 * fast);
    close(next.neutron_counts()[1], 5.0 * (slow - fast));
    close(step.precursor_transfer_counts[0], 5.0 * fast);
    close(step.delayed_decays[0], next.neutron_counts()[0]);
    close(step.delayed_decays[1], next.neutron_counts()[1]);
    close(step.integrated_fluence_m2[0], 0.5 * 5.0 * (duration - fast));
    close(
        step.integrated_fluence_m2[1],
        0.2 * (10.0 * (duration - slow) - 5.0 * (duration - fast)),
    );
    close(
        next.neutron_counts()
            .iter()
            .chain(next.precursor_counts())
            .sum(),
        10.0,
    );
    assert!(next.neutron_counts()[1] > 0.0);
    assert_eq!(step.escaped_precursors, 0.0);
    balances(&step);
    let (half, first) = initial.trial(1.0, budget(), || Ok(())).unwrap();
    let (split, second) = half.trial(1.0, budget(), || Ok(())).unwrap();
    for cell in 0..2 {
        close(split.neutron_counts()[cell], next.neutron_counts()[cell]);
        close(
            split.precursor_counts()[cell],
            next.precursor_counts()[cell],
        );
        close(
            first.integrated_fluence_m2[cell] + second.integrated_fluence_m2[cell],
            step.integrated_fluence_m2[cell],
        );
    }
    close(
        first.precursor_transfer_counts[0] + second.precursor_transfer_counts[0],
        step.precursor_transfer_counts[0],
    );
}

#[test]
fn closed_flow_loop_conserves_precursors_plus_emitted_neutrons_in_unequal_volumes() {
    let transport = PrecursorTransport {
        transfers: vec![flow(0, 1, 2.0), flow(1, 0, 2.0)],
        outflows: vec![],
    };
    let initial = state(data(&[2.0, 4.0], 0.25), transport, vec![12.0, 0.0]);
    let duration = 1.25;
    let (next, step) = initial.trial(duration, budget(), || Ok(())).unwrap();
    let surviving = (-0.25_f64 * duration).exp();
    let mixing = (-1.5_f64 * duration).exp();
    let residence_a =
        4.0 * survival_integral(0.25, duration) + 8.0 * survival_integral(1.75, duration);
    let residence_b = 8.0 * (survival_integral(0.25, duration) - survival_integral(1.75, duration));
    close(next.precursor_counts()[0], surviving * (4.0 + 8.0 * mixing));
    close(next.precursor_counts()[1], surviving * (8.0 - 8.0 * mixing));
    close(next.neutron_counts()[0], 0.25 * residence_a);
    close(next.neutron_counts()[1], 0.25 * residence_b);
    close(step.precursor_transfer_counts[0], residence_a);
    close(step.precursor_transfer_counts[1], 0.5 * residence_b);
    close(
        next.neutron_counts()
            .iter()
            .chain(next.precursor_counts())
            .sum(),
        12.0,
    );
    close(step.delayed_neutrons, 12.0 * (1.0 - surviving));
    balances(&step);
}

#[test]
fn continuous_precursor_feed_competes_with_outlet_escape_and_local_decay() {
    let transport = PrecursorTransport {
        transfers: vec![],
        outflows: vec![PrecursorOutflow {
            cell: 0,
            volumetric_flow_m3_s: 1.0,
        }],
    };
    let initial = state(data(&[2.0], 0.5), transport, vec![9.0]);
    let duration = 2.0;
    let (next, step) = initial
        .trial_with_sources(duration, &[0.0], &[3.0], budget(), || Ok(()))
        .unwrap();
    let residence_factor = survival_integral(1.0, duration);
    let precursor_integral = 9.0 * residence_factor + 3.0 * (duration - residence_factor);
    close(
        next.precursor_counts()[0],
        9.0 * (-duration).exp() + 3.0 * residence_factor,
    );
    close(step.precursor_decays, 0.5 * precursor_integral);
    close(step.escaped_precursors, 0.5 * precursor_integral);
    close(step.precursor_outflow_counts[0], step.escaped_precursors);
    close(next.neutron_counts()[0], step.precursor_decays);
    close(
        step.integrated_fluence_m2[0],
        0.25 * (6.0 * (duration - residence_factor) + 1.5 * duration * duration),
    );
    close(
        next.precursor_counts()[0] + next.neutron_counts()[0] + step.escaped_precursors,
        15.0,
    );
    assert_eq!(step.injected_precursors, 6.0);
    assert_eq!(step.leaked_neutrons, 0.0);
    balances(&step);
}

#[test]
fn rebuilding_material_coefficients_preserves_flow_and_rebinding_guards_endpoint_identity() {
    let supplied = data(&[2.0, 5.0], 0.5);
    let transport = PrecursorTransport {
        transfers: vec![flow(0, 1, 1.0)],
        outflows: vec![PrecursorOutflow {
            cell: 1,
            volumetric_flow_m3_s: 0.5,
        }],
    };
    let initial = state(supplied.clone(), transport.clone(), vec![5.0, 2.0]);
    let (accepted, _) = initial.trial(0.5, budget(), || Ok(())).unwrap();
    let mut changed = supplied.clone();
    changed.reactions.push(DiffusionReaction {
        data_id: "manufactured-rebound-capture".into(),
        cell: 0,
        incident_group: 0,
        macroscopic_cross_section_m_inv: 0.25,
        kind: DiffusionReactionKind::Capture,
    });
    let rebuilt = Arc::new(
        accepted
            .model()
            .rebuild_with_data(changed, DiffusionLimits::default(), || Ok(()))
            .unwrap(),
    );
    assert_eq!(rebuilt.precursor_transport(), &transport);
    let rebound = accepted.rebind_model(rebuilt).unwrap();
    assert_eq!(rebound.precursor_counts(), accepted.precursor_counts());
    assert!(rebound.prepared.is_none());
    let mut faster = transport.clone();
    faster.transfers[0].volumetric_flow_m3_s = 2.0;
    let faster = Arc::new(
        DiffusionModel::new_with_transport(supplied.clone(), faster, DiffusionLimits::default())
            .unwrap(),
    );
    assert!(accepted.rebind_model(faster).is_ok());
    let dropped =
        Arc::new(DiffusionModel::new(supplied.clone(), DiffusionLimits::default()).unwrap());
    assert!(accepted.rebind_model(dropped).is_err());
    let mut reversed = transport.clone();
    reversed.transfers[0] = flow(1, 0, 1.0);
    assert!(
        accepted
            .rebind_model(Arc::new(
                DiffusionModel::new_with_transport(
                    supplied.clone(),
                    reversed,
                    DiffusionLimits::default()
                )
                .unwrap()
            ))
            .is_err()
    );
    let mut outlet_changed = transport;
    outlet_changed.outflows[0].cell = 0;
    assert!(
        accepted
            .rebind_model(Arc::new(
                DiffusionModel::new_with_transport(
                    supplied,
                    outlet_changed,
                    DiffusionLimits::default()
                )
                .unwrap()
            ))
            .is_err()
    );
}

#[test]
fn shared_flow_budgets_invalid_paths_and_unrepresentable_rates_refuse_without_mutation() {
    let mut supplied = data(&[2.0, 5.0], 0.5);
    supplied.links.push(DiffusionLink {
        cell_a: 0,
        cell_b: 1,
        conductance_m2: vec![0.1],
    });
    supplied.leakage.push(LeakageBoundary {
        cell: 1,
        conductance_m2: vec![0.1],
    });
    let transport = PrecursorTransport {
        transfers: vec![flow(0, 1, 1.0)],
        outflows: vec![PrecursorOutflow {
            cell: 1,
            volumetric_flow_m3_s: 0.5,
        }],
    };
    for limits in [
        DiffusionLimits {
            max_links: 1,
            ..DiffusionLimits::default()
        },
        DiffusionLimits {
            max_leakage_boundaries: 1,
            ..DiffusionLimits::default()
        },
        DiffusionLimits {
            max_coefficients: 7,
            ..DiffusionLimits::default()
        },
    ] {
        assert!(matches!(
            DiffusionModel::new_with_transport(supplied.clone(), transport.clone(), limits),
            Err(DiffusionError::Budget(_))
        ));
    }
    for transfer in [
        flow(0, 0, 1.0),
        flow(0, 2, 1.0),
        flow(0, 1, -1.0),
        flow(0, 1, f64::INFINITY),
    ] {
        let invalid = PrecursorTransport {
            transfers: vec![transfer],
            outflows: vec![],
        };
        assert!(matches!(
            DiffusionModel::new_with_transport(
                supplied.clone(),
                invalid,
                DiffusionLimits::default()
            ),
            Err(DiffusionError::Invalid(_))
        ));
    }
    let underflow = PrecursorTransport {
        transfers: vec![flow(0, 1, f64::from_bits(1))],
        outflows: vec![],
    };
    assert!(matches!(
        DiffusionModel::new_with_transport(supplied.clone(), underflow, DiffusionLimits::default()),
        Err(DiffusionError::Arithmetic(_))
    ));
    let overflow = PrecursorTransport {
        transfers: vec![flow(0, 1, f64::MAX)],
        outflows: vec![],
    };
    assert!(matches!(
        DiffusionModel::new_with_transport(
            data(&[0.5, 5.0], 0.5),
            overflow,
            DiffusionLimits::default()
        ),
        Err(DiffusionError::Arithmetic(_))
    ));
    let mut callbacks = 0;
    DiffusionModel::new_with_transport_and_checkpoint(
        supplied.clone(),
        transport.clone(),
        DiffusionLimits::default(),
        || {
            callbacks += 1;
            Ok(())
        },
    )
    .unwrap();
    for stop in [1, callbacks / 2, callbacks] {
        let mut visited = 0;
        let result = DiffusionModel::new_with_transport_and_checkpoint(
            supplied.clone(),
            transport.clone(),
            DiffusionLimits::default(),
            || {
                visited += 1;
                if visited == stop {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        );
        assert_eq!(
            result.unwrap_err(),
            DiffusionError::Fission(FissionError::Cancelled)
        );
    }
    let mut initial = state(supplied, transport, vec![5.0, 2.0]);
    let before = initial.precursor_counts().to_vec();
    let result = initial.advance(1.0, budget(), || Err(FissionError::Cancelled));
    assert_eq!(
        result.unwrap_err(),
        DiffusionError::Fission(FissionError::Cancelled)
    );
    assert_eq!(initial.precursor_counts(), before);
}
