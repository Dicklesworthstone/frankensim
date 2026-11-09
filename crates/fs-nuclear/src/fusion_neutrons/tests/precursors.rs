//! Manufactured G1 delayed-source/capture and directed-flow solutions.
//! These test the shared particle history and disjoint energy accounts.

use super::*;
use crate::diffusion::DelayedGroup;
use crate::diffusion::transport::{PrecursorOutflow, PrecursorTransfer, PrecursorTransport};
use crate::fusion_neutrons::history::{
    FusionNeutronHistory, FusionNeutronHistoryLimits, NeutronKineticHistoryBudget,
};

fn supplied(cells: usize) -> DiffusionData {
    let mut supplied = data(&[1.0, 4.0], cells);
    supplied.delayed_groups.push(DelayedGroup {
        data_id: "manufactured-delayed-family".into(),
        decay_constant_s_inv: 0.4,
        spectrum: vec![0.25, 0.75],
    });
    supplied
}

fn owner(
    supplied: DiffusionData,
    transport: PrecursorTransport,
    fractions: Vec<f64>,
) -> Arc<NeutronKineticModel> {
    let reactions = vec![NeutronReactionEnergy::default(); supplied.reactions.len()];
    let diffusion = Arc::new(
        DiffusionModel::new_with_transport(supplied, transport, DiffusionLimits::default()).unwrap(),
    );
    Arc::new(
        NeutronKineticModel::new_with_precursor_storage(
            "manufactured-reserved-neutron-energy".into(),
            diffusion,
            vec![MEV, 4.0 * MEV],
            fractions,
            reactions,
        )
        .unwrap(),
    )
}

fn capturing() -> Arc<NeutronKineticModel> {
    let mut data = supplied(1);
    for (group, rate) in [1.3, 0.8].into_iter().enumerate() {
        data.reactions
            .push(reaction(0, group, rate, DiffusionReactionKind::Capture));
    }
    owner(data, PrecursorTransport::default(), vec![1.0, 0.5])
}

#[test]
fn initial_precursor_energy_becomes_delayed_neutrons_and_actual_capture_heat() {
    let model = capturing();
    close(model.precursor_energies_j()[0], 3.25 * MEV);
    let stock = 37.0;
    let duration: f64 = 1.7;
    let initial =
        NeutronKineticState::new_with_precursors(model, vec![0.0; 2], vec![stock]).unwrap();
    assert_eq!(initial.kinetic_energy_j().unwrap(), 0.0);
    close(initial.stored_energy_j().unwrap(), stock * 3.25 * MEV);
    let (next, step) = initial
        .trial_with_rates(duration, &[0.0; 2], budget(), || Ok(()))
        .unwrap();
    let remaining = stock * (-0.4 * duration).exp();
    let decays = stock - remaining;
    let mut heat = 0.0;
    let mut radiation = 0.0;
    for (group, (capture, weight, fraction)) in
        [(1.3, 0.25, 1.0), (0.8, 0.75, 0.5)].into_iter().enumerate()
    {
        let neutrons = stock * weight * 0.4
            * ((-0.4 * duration).exp() - (-capture * duration).exp())
            / (capture - 0.4);
        let captures = weight * decays - neutrons;
        let energy = [MEV, 4.0 * MEV][group];
        close(next.neutrons().neutron_counts()[group], neutrons);
        close(step.diffusion.reaction_counts[group], captures);
        heat += fraction * captures * energy;
        radiation += (1.0 - fraction) * captures * energy;
    }
    close(next.neutrons().precursor_counts()[0], remaining);
    close(step.initial_precursor_energy_j, stock * 3.25 * MEV);
    close(step.final_precursor_energy_j, remaining * 3.25 * MEV);
    close(step.precursor_conversion_energy_j, decays * 3.25 * MEV);
    close(step.released.deposited_j, heat);
    close(step.released.radiation_j, radiation);
    assert_eq!(step.released.escaped_j, 0.0);
    assert_eq!(step.nuclear_reaction_energy_j, 0.0);
    assert_eq!(step.injected_precursor_energy_j, 0.0);
    assert_eq!(step.injected_kinetic_energy_j, 0.0);
    assert_eq!(step.exported_precursor_energy_j, 0.0);
    close(
        next.stored_energy_j().unwrap() + step.released.total_j(),
        initial.stored_energy_j().unwrap(),
    );
    assert!(step.energy_roundoff_j <= budget().energy_tolerance_j);
    assert_eq!(initial.neutrons().precursor_counts(), &[stock]);
}

#[test]
fn external_precursor_feed_reserves_energy_and_resolves_birth_residence_time() {
    let initial = NeutronKineticState::new(capturing(), vec![0.0; 2]).unwrap();
    assert_eq!(initial.neutrons().precursor_counts(), &[0.0]);
    let rate = 5.0;
    let duration: f64 = 2.3;
    let (next, step) = initial
        .trial_with_sources(duration, &[0.0; 2], &[rate], budget(), || Ok(()))
        .unwrap();
    let remaining = rate * -(-0.4 * duration).exp_m1() / 0.4;
    let decays = rate * duration - remaining;
    close(next.neutrons().precursor_counts()[0], remaining);
    close(step.diffusion.injected_precursors, rate * duration);
    close(step.injected_precursor_energy_j, rate * duration * 3.25 * MEV);
    close(step.precursor_conversion_energy_j, decays * 3.25 * MEV);
    assert_eq!(step.injected_kinetic_energy_j, 0.0);
    for (group, (capture, weight)) in [(1.3, 0.25), (0.8, 0.75)].into_iter().enumerate() {
        let neutrons = rate * weight
            * (-(-capture * duration).exp_m1() / capture
                - ((-0.4 * duration).exp() - (-capture * duration).exp()) / (capture - 0.4));
        close(next.neutrons().neutron_counts()[group], neutrons);
        close(step.diffusion.reaction_counts[group], weight * decays - neutrons);
    }
    close(
        next.stored_energy_j().unwrap() + step.released.total_j(),
        step.injected_precursor_energy_j,
    );
    let (tail, shutdown) = next
        .trial_with_rates(1.0, &[0.0; 2], budget(), || Ok(()))
        .unwrap();
    assert_eq!(shutdown.injected_precursor_energy_j, 0.0);
    close(tail.neutrons().precursor_counts()[0], remaining * (-0.4_f64).exp());
    close(
        tail.stored_energy_j().unwrap()
            + tail.cumulative_cell_released()[0].total_j(),
        step.injected_precursor_energy_j,
    );
}

fn flowing() -> Arc<NeutronKineticModel> {
    let mut data = supplied(2);
    data.delayed_groups.push(DelayedGroup {
        data_id: "manufactured-second-family".into(),
        decay_constant_s_inv: 0.7,
        spectrum: vec![1.0, 0.0],
    });
    owner(
        data,
        PrecursorTransport {
            transfers: vec![PrecursorTransfer {
                from_cell: 0,
                to_cell: 1,
                volumetric_flow_m3_s: 0.6,
            }],
            outflows: vec![PrecursorOutflow {
                cell: 1,
                volumetric_flow_m3_s: 0.2,
            }],
        },
        vec![],
    )
}

#[test]
fn transported_families_keep_their_energy_and_export_only_the_undecayed_reserve() {
    let model = flowing();
    let initial = NeutronKineticState::new_with_precursors(
        model,
        vec![0.0; 4],
        vec![40.0, 13.0, 0.0, 0.0],
    )
    .unwrap();
    let duration: f64 = 1.9;
    let (next, step) = initial
        .trial_with_rates(duration, &[0.0; 4], budget(), || Ok(()))
        .unwrap();
    let mut expected_export = 0.0;
    let mut expected_neutrons = [0.0; 4];
    for (family, (stock, lambda, energy, spectrum)) in [
        (40.0, 0.4, 3.25 * MEV, [0.25, 0.75]),
        (13.0, 0.7, MEV, [1.0, 0.0]),
    ]
    .into_iter()
    .enumerate()
    {
        let a = lambda + 0.6;
        let b = lambda + 0.2;
        let donor = stock * (-a * duration).exp();
        let receiver = 0.6 * stock * ((-a * duration).exp() - (-b * duration).exp()) / (b - a);
        let donor_integral = stock * -(-a * duration).exp_m1() / a;
        let receiver_integral = 0.6 * stock
            * (-(-a * duration).exp_m1() / a + (-b * duration).exp_m1() / b)
            / (b - a);
        close(next.neutrons().precursor_counts()[family], donor);
        close(next.neutrons().precursor_counts()[2 + family], receiver);
        close(step.diffusion.precursor_outflow_counts[family], 0.2 * receiver_integral);
        expected_export += 0.2 * receiver_integral * energy;
        for (group, weight) in spectrum.into_iter().enumerate() {
            expected_neutrons[group] += lambda * donor_integral * weight;
            expected_neutrons[2 + group] += lambda * receiver_integral * weight;
        }
    }
    for (&actual, expected) in next.neutrons().neutron_counts().iter().zip(expected_neutrons) {
        close(actual, expected);
    }
    close(step.exported_precursor_energy_j, expected_export);
    assert_eq!(step.cell_exported_precursor_energy_j[0], 0.0);
    close(step.cell_exported_precursor_energy_j[1], expected_export);
    close(next.cell_exported_precursor_energy_j(1).unwrap(), expected_export);
    close(next.cumulative_exported_precursor_energy_j().unwrap(), expected_export);
    assert_eq!(step.released, EnergyChannels::default());
    close(
        next.stored_energy_j().unwrap() + expected_export,
        initial.stored_energy_j().unwrap(),
    );
    let (tail, later) = next
        .trial_with_rates(0.8, &[0.0; 4], budget(), || Ok(()))
        .unwrap();
    close(
        tail.cumulative_exported_precursor_energy_j().unwrap(),
        expected_export + later.exported_precursor_energy_j,
    );
    close(
        tail.stored_energy_j().unwrap()
            + tail.cumulative_exported_precursor_energy_j().unwrap(),
        initial.stored_energy_j().unwrap(),
    );
}

#[test]
fn existing_history_accounts_precursor_storage_and_outflow_across_restart() {
    let model = flowing();
    let history = FusionNeutronHistory::project(
        model.clone(),
        vec![],
        FusionNeutronHistoryLimits::default(),
    )
    .unwrap();
    let initial = NeutronKineticState::new_with_precursors(
        model,
        vec![0.0; 4],
        vec![40.0, 13.0, 0.0, 0.0],
    )
    .unwrap();
    let limits = NeutronKineticHistoryBudget {
        transport: budget(),
        maximum_step_s: 0.3,
        maximum_substeps: 20,
    };
    let (whole, report) = initial
        .trial_with_history(Time::new(0.0), Time::new(2.0), &history, limits, || Ok(()))
        .unwrap();
    let (half, first) = initial
        .trial_with_history(Time::new(0.0), Time::new(1.0), &history, limits, || Ok(()))
        .unwrap();
    let (split, second) = half
        .trial_with_history(Time::new(1.0), Time::new(2.0), &history, limits, || Ok(()))
        .unwrap();
    close(report.initial_precursor_energy_j, initial.precursor_energy_j().unwrap());
    close(report.final_precursor_energy_j, whole.precursor_energy_j().unwrap());
    close(report.exported_precursor_energy_j, whole.cumulative_exported_precursor_energy_j().unwrap());
    close(report.exported_precursor_energy_j, first.exported_precursor_energy_j + second.exported_precursor_energy_j);
    close(report.precursor_conversion_energy_j, first.precursor_conversion_energy_j + second.precursor_conversion_energy_j);
    close(whole.stored_energy_j().unwrap(), split.stored_energy_j().unwrap());
    close(
        whole.stored_energy_j().unwrap() + report.exported_precursor_energy_j,
        initial.stored_energy_j().unwrap(),
    );
    assert_eq!(report.released, EnergyChannels::default());
    assert!(report.energy_roundoff_j <= limits.transport.energy_tolerance_j);
}

#[test]
fn explicit_admission_rejects_unfunded_induced_families_and_preserves_legacy_refusals() {
    let native = Arc::new(DiffusionModel::new(supplied(1), DiffusionLimits::default()).unwrap());
    assert!(NeutronKineticModel::new(
        "legacy".into(), native.clone(), vec![MEV, 4.0 * MEV], vec![],
    ).is_err());
    assert!(NeutronKineticModel::new_with_reaction_energies(
        "legacy-explicit-q".into(), native, vec![MEV, 4.0 * MEV], vec![], vec![],
    ).is_err());
    for (yield_per_fission, admitted) in [(0.0, true), (0.03, false)] {
        let mut data = supplied(1);
        data.reactions.push(reaction(
            0, 1, 0.2,
            DiffusionReactionKind::Fission {
                prompt_yields: vec![2.0, 0.0],
                delayed_yields: vec![yield_per_fission],
            },
        ));
        let native = Arc::new(DiffusionModel::new(data, DiffusionLimits::default()).unwrap());
        let result = NeutronKineticModel::new_with_precursor_storage(
            "explicit-precursor-storage".into(),
            native,
            vec![MEV, 4.0 * MEV],
            vec![1.0],
            vec![NeutronReactionEnergy { nuclear_q_j: MEV, threshold_j: 0.0 }],
        );
        assert_eq!(result.is_ok(), admitted);
    }
}

#[test]
fn precursor_sources_and_late_cancellation_refuse_without_changing_storage_or_exports() {
    let initial = NeutronKineticState::new_with_precursors(
        flowing(), vec![0.0; 4], vec![40.0, 13.0, 0.0, 0.0],
    ).unwrap();
    assert!(NeutronKineticState::new_with_precursors(
        initial.model().clone(), vec![0.0; 4], vec![1.0],
    ).is_err());
    assert!(initial.trial_with_sources(
        1.0, &[0.0; 4], &[1.0], budget(), || Ok(()),
    ).is_err());
    assert!(initial.trial_with_sources(
        1.0, &[0.0; 4], &[0.0, -1.0, 0.0, 0.0], budget(), || Ok(()),
    ).is_err());
    let mut calls = 0;
    initial.trial_with_rates(1.0, &[0.0; 4], budget(), || {
        calls += 1;
        Ok(())
    }).unwrap();
    let mut visited = 0;
    let result = initial.trial_with_rates(1.0, &[0.0; 4], budget(), || {
        visited += 1;
        if visited == calls { Err(FissionError::Cancelled) } else { Ok(()) }
    });
    assert!(matches!(result, Err(NeutronKineticError::Source(FissionError::Cancelled))));
    let mut too_tight = budget();
    too_tight.energy_tolerance_j = 0.0;
    assert!(matches!(
        initial.trial_with_rates(1.0, &[0.0; 4], too_tight, || Ok(())),
        Err(NeutronKineticError::Balance(_)),
    ));
    assert_eq!(initial.neutrons().precursor_counts(), &[40.0, 13.0, 0.0, 0.0]);
    assert_eq!(initial.neutrons().neutron_counts(), &[0.0; 4]);
    assert_eq!(initial.cumulative_exported_precursor_energy_j().unwrap(), 0.0);
    assert_eq!(initial.cumulative_cell_released(), vec![EnergyChannels::default(); 2]);
}
