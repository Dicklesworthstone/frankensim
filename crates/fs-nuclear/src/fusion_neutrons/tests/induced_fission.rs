//! G1 exact fast-neutron -> fission -> precursor/slow-neutron -> capture.
//! Joule-scale energies, rates and source yields are manufactured fixtures.

use super::*;
use crate::diffusion::DelayedGroup;
use crate::diffusion::transport::{PrecursorOutflow, PrecursorTransport};
use crate::fusion_neutrons::history::{
    FusionNeutronHistory, FusionNeutronHistoryLimits, NeutronKineticHistoryBudget,
};

fn supplied() -> DiffusionData {
    DiffusionData {
        data_id: "manufactured-funded-induced-fission".into(),
        groups: Arc::new(
            NeutronGroups::new(
                "manufactured-induced-fission-grid".into(),
                vec![0.0, 3.0, 6.0],
            )
            .unwrap(),
        ),
        speeds_m_s: vec![1.0; 2],
        cell_volumes_m3: vec![0.01],
        delayed_groups: vec![DelayedGroup {
            data_id: "manufactured-induced-fission-family".into(),
            decay_constant_s_inv: 2.0,
            spectrum: vec![1.0, 0.0],
        }],
        reactions: vec![
            reaction(
                0,
                1,
                1.0,
                DiffusionReactionKind::Fission {
                    prompt_yields: vec![1.0, 0.0],
                    delayed_yields: vec![0.5],
                },
            ),
            reaction(0, 0, 0.25, DiffusionReactionKind::Capture),
        ],
        links: vec![],
        leakage: vec![],
    }
}

fn declarations(cap: f64) -> Vec<Option<InducedFissionPrecursorFunding>> {
    vec![
        Some(
            InducedFissionPrecursorFunding::AdditionalDelayedConversion {
                maximum_energy_j_per_fission: cap,
            },
        ),
        None,
    ]
}

fn energies() -> Vec<NeutronReactionEnergy> {
    vec![
        NeutronReactionEnergy {
            nuclear_q_j: 6.0,
            threshold_j: 0.0,
        },
        NeutronReactionEnergy::default(),
    ]
}

fn model(outflow_rate_s: f64, cap: f64) -> Arc<NeutronKineticModel> {
    let transport = PrecursorTransport {
        transfers: vec![],
        outflows: if outflow_rate_s == 0.0 {
            vec![]
        } else {
            vec![PrecursorOutflow {
                cell: 0,
                volumetric_flow_m3_s: 0.01 * outflow_rate_s,
            }]
        },
    };
    Arc::new(
        NeutronKineticModel::new_with_fission_precursor_funding(
            "manufactured-independent-prompt-and-delayed-conversion".into(),
            Arc::new(
                DiffusionModel::new_with_transport(
                    supplied(),
                    transport,
                    DiffusionLimits::default(),
                )
                .unwrap(),
            ),
            vec![2.0, 4.0],
            vec![0.25; 2],
            energies(),
            declarations(cap),
        )
        .unwrap(),
    )
}

/// Shared original native source for the existing material owner.
pub(crate) fn material_induced_fission_fixture(
    initial_fast: f64,
    precursor_outflow_rate_s: f64,
) -> NeutronKineticState {
    NeutronKineticState::new(
        model(precursor_outflow_rate_s, 2.0),
        vec![0.0, initial_fast],
    )
    .unwrap()
}

fn budget() -> NeutronKineticBudget {
    NeutronKineticBudget {
        diffusion: DiffusionBudget {
            limits: DiffusionLimits::default(),
            population_tolerance: 1e-8,
        },
        energy_tolerance_j: 1e-8,
    }
}

fn near(actual: f64, expected: f64) {
    let tolerance = 1e-9 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:.3e}",
    );
}

// Closed-form Bateman histories are independent of the native propagator.
fn expected(stock: f64, time: f64, outflow: f64) -> [f64; 6] {
    let removal = 2.0 + outflow;
    let e1 = (-time).exp();
    let er = (-removal * time).exp();
    let response = |rate: f64| ((-rate * time).exp() - (-0.25 * time).exp()) / (0.25 - rate);
    let fissions = stock * -(-time).exp_m1();
    let precursors = 0.5 * stock * (e1 - er) / (removal - 1.0);
    let residence =
        0.5 * stock / (removal - 1.0) * (-(-time).exp_m1() + (-removal * time).exp_m1() / removal);
    let slow =
        stock * response(1.0) + stock / (removal - 1.0) * (response(1.0) - response(removal));
    let captures = fissions + 2.0 * residence - slow;
    [
        stock * e1,
        precursors,
        slow,
        fissions,
        captures,
        outflow * residence,
    ]
}

#[test]
fn g1_induced_family_energy_is_separate_from_prompt_q_and_actual_capture_heat() {
    let initial = material_induced_fission_fixture(100.0, 0.0);
    let duration = 1.25;
    let [fast, families, slow, fissions, captures, _] = expected(100.0, duration, 0.0);
    let (after, step) = initial
        .trial_with_rates(duration, &[0.0; 2], budget(), || Ok(()))
        .unwrap();
    near(after.neutrons().neutron_counts()[1], fast);
    near(after.neutrons().neutron_counts()[0], slow);
    near(after.neutrons().precursor_counts()[0], families);
    near(step.diffusion.fissions, fissions);
    near(step.diffusion.captures, captures);
    near(step.diffusion.delayed_precursors_produced, 0.5 * fissions);
    assert_eq!(step.diffusion.injected_precursors, 0.0);
    assert_eq!(step.injected_precursor_energy_j, 0.0);
    assert_eq!(step.injected_kinetic_energy_j, 0.0);
    near(step.nuclear_reaction_energy_j, 6.0 * fissions);
    near(step.induced_precursor_energy_j, fissions);
    near(step.cell_induced_precursor_energy_j[0], fissions);
    near(
        after.cumulative_induced_precursor_energy_j().unwrap(),
        fissions,
    );
    near(after.cell_induced_precursor_energy_j(0).unwrap(), fissions);
    assert_eq!(after.cell_induced_precursor_energy_j(1), None);
    near(step.final_precursor_energy_j, 2.0 * families);
    near(
        step.precursor_conversion_energy_j,
        fissions - 2.0 * families,
    );
    near(step.released.deposited_j, 2.0 * fissions + 0.5 * captures);
    near(step.released.radiation_j, 6.0 * fissions + 1.5 * captures);
    assert_eq!(step.released.escaped_j, 0.0);
    near(
        after.stored_energy_j().unwrap() + step.released.total_j(),
        400.0 + 7.0 * fissions,
    );
    assert!(step.energy_roundoff_j <= budget().energy_tolerance_j);
    assert_eq!(
        initial.cumulative_induced_precursor_energy_j().unwrap(),
        0.0
    );
    assert_eq!(initial.neutrons().precursor_counts(), &[0.0]);

    // A larger cap changes admission only; none of its unused energy is a source.
    let larger = NeutronKineticState::new(model(0.0, 200.0), vec![0.0, 100.0]).unwrap();
    let (same, same_step) = larger
        .trial_with_rates(duration, &[0.0; 2], budget(), || Ok(()))
        .unwrap();
    assert_eq!(
        same.neutrons().neutron_counts(),
        after.neutrons().neutron_counts()
    );
    assert_eq!(same_step.released, step.released);
    assert_eq!(
        same_step.induced_precursor_energy_j,
        step.induced_precursor_energy_j
    );
    near(initial.model().fission_precursor_energy_j()[0], 1.0);
    assert_eq!(initial.model().fission_precursor_energy_j()[1], 0.0);
}

#[test]
fn g0_funding_admits_only_original_fission_layout_and_positive_reserves() {
    let native = Arc::new(DiffusionModel::new(supplied(), DiffusionLimits::default()).unwrap());
    let make = |declarations| {
        NeutronKineticModel::new_with_fission_precursor_funding(
            "manufactured-funding-admission".into(),
            Arc::clone(&native),
            vec![2.0, 4.0],
            vec![0.25; 2],
            energies(),
            declarations,
        )
    };
    for cap in [0.0, 0.99, -1.0, f64::NAN, f64::INFINITY] {
        assert!(make(declarations(cap)).is_err());
    }
    assert!(make(vec![None; 2]).is_err());
    assert!(make(vec![]).is_err());
    let mut wrong_order = declarations(2.0);
    wrong_order.swap(0, 1);
    assert!(make(wrong_order).is_err());
    assert!(
        NeutronKineticModel::new_with_precursor_storage(
            "legacy-explicit-storage".into(),
            Arc::clone(&native),
            vec![2.0, 4.0],
            vec![0.25; 2],
            energies(),
        )
        .is_err()
    );
    assert!(
        NeutronKineticModel::new_with_reaction_energies(
            "legacy-prompt-conversion".into(),
            native,
            vec![2.0, 4.0],
            vec![0.25; 2],
            energies(),
        )
        .is_err()
    );

    let mut inactive = supplied();
    inactive.reactions[0].macroscopic_cross_section_m_inv = 0.0;
    assert!(
        NeutronKineticModel::new_with_fission_precursor_funding(
            "inactive-still-requires-original-funding".into(),
            Arc::new(DiffusionModel::new(inactive, DiffusionLimits::default()).unwrap()),
            vec![2.0, 4.0],
            vec![0.25; 2],
            energies(),
            vec![None; 2],
        )
        .is_err()
    );

    // The smallest positive yield cannot disappear from preparation merely
    // because its per-event energy product rounds to zero.
    let mut tiny = supplied();
    tiny.groups = Arc::new(
        NeutronGroups::new("tiny-representative-grid".into(), vec![0.0, 0.25, 1.0]).unwrap(),
    );
    tiny.reactions[0].kind = DiffusionReactionKind::Fission {
        prompt_yields: vec![1.0, 0.0],
        delayed_yields: vec![f64::from_bits(1)],
    };
    tiny.reactions[0].macroscopic_cross_section_m_inv = 0.0;
    assert!(
        NeutronKineticModel::new_with_fission_precursor_funding(
            "positive-reserve-underflow-refuses".into(),
            Arc::new(DiffusionModel::new(tiny, DiffusionLimits::default()).unwrap()),
            vec![0.125, 0.5],
            vec![0.25; 2],
            energies(),
            declarations(2.0),
        )
        .is_err()
    );
}

#[test]
fn g0_delayed_conversion_cannot_fund_prompt_energy_or_lower_thresholds() {
    let make = |data, energy| {
        NeutronKineticModel::new_with_fission_precursor_funding(
            "separate-prompt-and-future-energy".into(),
            Arc::new(DiffusionModel::new(data, DiffusionLimits::default()).unwrap()),
            vec![2.0, 4.0],
            vec![0.25; 2],
            vec![energy, NeutronReactionEnergy::default()],
            declarations(1e6),
        )
    };
    // Incident4 + promptQ(-3) < outgoing2; future funding cannot fill the gap.
    assert!(
        make(
            supplied(),
            NeutronReactionEnergy {
                nuclear_q_j: -3.0,
                threshold_j: 3.0
            },
        )
        .is_err()
    );
    assert!(
        make(
            supplied(),
            NeutronReactionEnergy {
                nuclear_q_j: 6.0,
                threshold_j: 5.0
            },
        )
        .is_err()
    );
    let mut inactive = supplied();
    inactive.reactions[0].macroscopic_cross_section_m_inv = 0.0;
    let model = make(
        inactive,
        NeutronReactionEnergy {
            nuclear_q_j: 6.0,
            threshold_j: 5.0,
        },
    )
    .unwrap();
    let active = DiffusionModel::new(supplied(), DiffusionLimits::default()).unwrap();
    assert!(model.admit_active_reactions(&active).is_err());
}

#[test]
fn g1_g3_induced_history_retains_storage_exports_and_conversion_across_restart() {
    let initial = material_induced_fission_fixture(100.0, 0.6);
    let history = FusionNeutronHistory::project(
        Arc::clone(initial.model()),
        vec![],
        FusionNeutronHistoryLimits::default(),
    )
    .unwrap();
    let history_budget = NeutronKineticHistoryBudget {
        transport: budget(),
        maximum_step_s: 0.125,
        maximum_substeps: 32,
    };
    let (after, whole) = initial
        .trial_with_history(
            Time::new(0.0),
            Time::new(2.0),
            &history,
            history_budget,
            || Ok(()),
        )
        .unwrap();
    let [fast, families, slow, fissions, captures, escaped] = expected(100.0, 2.0, 0.6);
    near(after.neutrons().neutron_counts()[1], fast);
    near(after.neutrons().neutron_counts()[0], slow);
    near(after.neutrons().precursor_counts()[0], families);
    near(whole.induced_precursor_energy_j, fissions);
    near(whole.cell_induced_precursor_energy_j[0], fissions);
    near(whole.nuclear_reaction_energy_j, 6.0 * fissions);
    near(whole.exported_precursor_energy_j, 2.0 * escaped);
    near(whole.cell_exported_precursor_energy_j[0], 2.0 * escaped);
    near(whole.released.deposited_j, 2.0 * fissions + 0.5 * captures);
    assert_eq!(whole.released.escaped_j, 0.0);
    near(
        after.stored_energy_j().unwrap()
            + whole.released.total_j()
            + whole.exported_precursor_energy_j,
        400.0 + 7.0 * fissions,
    );
    let (middle, first) = initial
        .trial_with_history(
            Time::new(0.0),
            Time::new(0.75),
            &history,
            history_budget,
            || Ok(()),
        )
        .unwrap();
    let (split, second) = middle
        .trial_with_history(
            Time::new(0.75),
            Time::new(2.0),
            &history,
            history_budget,
            || Ok(()),
        )
        .unwrap();
    near(
        split.stored_energy_j().unwrap(),
        after.stored_energy_j().unwrap(),
    );
    near(
        first.induced_precursor_energy_j + second.induced_precursor_energy_j,
        whole.induced_precursor_energy_j,
    );
    near(
        split.cumulative_induced_precursor_energy_j().unwrap(),
        fissions,
    );
    near(
        split.cell_nuclear_reaction_energy_j(0).unwrap(),
        6.0 * fissions,
    );
    near(
        split.cumulative_exported_precursor_energy_j().unwrap(),
        2.0 * escaped,
    );
    near(
        first.released.total_j() + second.released.total_j(),
        whole.released.total_j(),
    );
    assert_eq!(whole.injected_neutrons, 0.0);
    assert_eq!(whole.injected_kinetic_energy_j, 0.0);
}

#[test]
fn g4_induced_energy_and_late_cancellation_leave_the_accepted_history_unchanged() {
    let initial = material_induced_fission_fixture(100.0, 0.6);
    let (accepted, _) = initial
        .trial_with_rates(0.5, &[0.0; 2], budget(), || Ok(()))
        .unwrap();
    let before = accepted.cumulative_induced_precursor_energy_j().unwrap();
    let exported = accepted.cumulative_exported_precursor_energy_j().unwrap();
    let populations = accepted.neutrons().neutron_counts().to_vec();
    let precursors = accepted.neutrons().precursor_counts().to_vec();
    let mut polls = 0;
    let (expected, report) = accepted
        .trial_with_rates(0.75, &[0.0; 2], budget(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    assert!(polls > 1);
    let mut seen = 0;
    assert!(
        accepted
            .trial_with_rates(0.75, &[0.0; 2], budget(), || {
                seen += 1;
                if seen == polls {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    let mut tight = budget();
    tight.energy_tolerance_j = 0.0;
    assert!(
        accepted
            .trial_with_rates(0.75, &[0.0; 2], tight, || Ok(()))
            .is_err()
    );
    assert_eq!(
        accepted.cumulative_induced_precursor_energy_j().unwrap(),
        before
    );
    assert_eq!(
        accepted.cumulative_exported_precursor_energy_j().unwrap(),
        exported
    );
    assert_eq!(accepted.neutrons().neutron_counts(), populations);
    assert_eq!(accepted.neutrons().precursor_counts(), precursors);
    let (retry, retry_step) = accepted
        .trial_with_rates(0.75, &[0.0; 2], budget(), || Ok(()))
        .unwrap();
    assert_eq!(
        retry.neutrons().neutron_counts(),
        expected.neutrons().neutron_counts()
    );
    assert_eq!(
        retry_step.induced_precursor_energy_j,
        report.induced_precursor_energy_j
    );
    near(
        retry.cumulative_induced_precursor_energy_j().unwrap(),
        before + report.induced_precursor_energy_j,
    );
}
