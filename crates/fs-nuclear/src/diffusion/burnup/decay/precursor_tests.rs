//! Manufactured parent -> precursor -> neutron histories with independent
//! analytic balances. Family parameters and yields are not isotope data.

use super::*;
use crate::burnup::transport::{
    BurnupTransportLimits, BurnupTransportState, NuclideOutflow, NuclideTransport,
};
use crate::diffusion::DelayedGroup;
use crate::diffusion::transport::{PrecursorOutflow, PrecursorTransport};

fn family() -> DelayedGroup {
    DelayedGroup {
        data_id: "manufactured-radioactive-family".into(),
        decay_constant_s_inv: 2.0,
        spectrum: vec![1.0],
    }
}

fn emission(model: &Arc<BurnupModel>) -> Arc<BurnupDecayNeutronEmission> {
    let mut neutrons = vec![vec![0.0]; model.nuclides().len()];
    let mut precursors = neutrons.clone();
    neutrons[0][0] = 0.25;
    precursors[0][0] = 0.5;
    Arc::new(
        BurnupDecayNeutronEmission::new_with_precursors(
            "manufactured-immediate-and-delayed-radioactive-yields".into(),
            Arc::clone(model),
            neutrons,
            vec![family()],
            precursors,
        )
        .unwrap(),
    )
}

fn fixture(capture: Option<(usize, f64)>, loss: f64, resident: f64) -> DiffusionBurnupState {
    let original = unbound(&[1.0, 0.0, 0.0], &[100.0, 0.0, 0.0], capture, loss);
    let mut data = original.neutrons().model().data().clone();
    data.delayed_groups = vec![family()];
    let model = Arc::new(DiffusionModel::new(data, budget().diffusion.limits).unwrap());
    let source = emission(original.inventories()[0].model());
    DiffusionBurnupState::new(
        DiffusionState::new(model, vec![resident], vec![0.0]).unwrap(),
        original.inventories().to_vec(),
        original.bindings().to_vec(),
    )
    .unwrap()
    .with_decay_neutrons(vec![source])
    .unwrap()
}

#[test]
fn g1_native_family_births_keep_both_ages_and_distinguish_external_feed() {
    let initial = fixture(None, 3.0, 0.0);
    let duration = 1.2_f64;
    let e1 = (-duration).exp();
    let e2 = (-2.0 * duration).exp();
    let e3 = (-3.0 * duration).exp();
    let decays = 100.0 * -(-duration).exp_m1();
    let internal_precursors = 50.0 * (e1 - e2);
    let internal_neutrons = 12.5 * (e1 - e3) + 50.0 * (e1 - 2.0 * e2 + e3);
    for external in [false, true] {
        let neutron_rate = if external { 3.0 } else { 0.0 };
        let precursor_rate = if external { 7.0 } else { 0.0 };
        let sources = DiffusionBurnupSources {
            neutron_rates_s: vec![neutron_rate],
            precursor_rates_s: vec![precursor_rate],
            nuclide_rates_s: vec![vec![0.0; 3]],
        };
        let (next, report) = initial
            .trial_with_sources(duration, &[300.0], &sources, budget(), || Ok(()))
            .unwrap();
        let precursors = internal_precursors + precursor_rate * (1.0 - e2) / 2.0;
        let neutrons = internal_neutrons + neutron_rate * (1.0 - e3) / 3.0
            + precursor_rate * ((1.0 - e3) / 3.0 - (e2 - e3));
        let precursor_births = 0.5 * decays + precursor_rate * duration;
        let neutron_births = 0.25 * decays + neutron_rate * duration;
        let family_decays = precursor_births - precursors;
        close(next.inventories()[0].populations()[0], 100.0 * e1, 1e-9);
        close(next.neutrons().precursor_counts()[0], precursors, 1e-9);
        close(next.neutrons().neutron_counts()[0], neutrons, 1e-9);
        close(report.decay_neutron_counts.as_ref().unwrap()[0][0], 0.25 * decays, 1e-9);
        close(report.decay_precursor_counts.as_ref().unwrap()[0][0], 0.5 * decays, 1e-9);
        close(report.neutrons.injected_neutron_counts[0], neutron_births, 1e-9);
        close(report.neutrons.injected_precursor_counts[0], precursor_births, 1e-9);
        close(report.neutrons.injected_precursors, precursor_births, 1e-9);
        close(report.neutrons.precursor_decays, family_decays, 1e-9);
        close(report.neutrons.delayed_neutrons, family_decays, 1e-9);
        close(report.neutrons.leaked_neutrons, neutron_births + family_decays - neutrons, 1e-9);
        assert_eq!(report.neutrons.delayed_precursors_produced, 0.0);
        assert!(report.neutrons.population_roundoff <= budget().diffusion.population_tolerance);

        let (middle, first) = initial
            .trial_with_sources(0.4, &[300.0], &sources, budget(), || Ok(()))
            .unwrap();
        let (restarted, second) = middle
            .trial_with_sources(0.8, &[300.0], &sources, budget(), || Ok(()))
            .unwrap();
        close(restarted.neutrons().precursor_counts()[0], precursors, 1e-9);
        close(restarted.neutrons().neutron_counts()[0], neutrons, 1e-9);
        close(first.neutrons.delayed_neutrons + second.neutrons.delayed_neutrons, family_decays, 1e-9);
    }
    let (unchanged, zero) = initial.trial(0.0, &[300.0], budget(), || Ok(())).unwrap();
    assert_eq!(unchanged.inventories(), initial.inventories());
    assert_eq!(zero.decay_precursor_counts, Some(vec![vec![0.0]]));
    assert_eq!(zero.decay_neutron_counts, Some(vec![vec![0.0]]));
    assert_eq!(unchanged.neutrons().precursor_counts(), &[0.0]);
}

// Independent nonlinear reference: parent, precursor, neutron, decay daughter,
// capture product and neutron leakage. No production numerical owner is used.
fn reference(steps: usize) -> [f64; 6] {
    fn rhs(x: [f64; 6]) -> [f64; 6] {
        let decay = x[0];
        let capture = 0.02 * x[0] * x[2];
        let delayed = 2.0 * x[1];
        let leak = 0.3 * x[2];
        [-decay - capture, 0.5 * decay - delayed,
            0.25 * decay + delayed - capture - leak, decay, capture, leak]
    }
    let h = 1.0 / steps as f64;
    let mut x = [100.0, 0.0, 10.0, 0.0, 0.0, 0.0];
    for _ in 0..steps {
        let a = rhs(x);
        let b = rhs(std::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = rhs(std::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = rhs(std::array::from_fn(|i| x[i] + h * c[i]));
        x = std::array::from_fn(|i| x[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    x
}

#[test]
fn g3_capture_consumes_the_same_parent_without_creating_precursors() {
    let expected = reference(2048);
    for (coarse, fine) in reference(1024).into_iter().zip(expected) {
        close(coarse, fine, 1e-9);
    }
    let mut errors = Vec::new();
    for steps in [32, 64] {
        let mut state = fixture(Some((0, 0.02)), 0.3, 10.0);
        let mut leakage = 0.0;
        let mut history_allowance = 0.0;
        for _ in 0..steps {
            let (next, step) = state.trial(1.0 / steps as f64, &[300.0], budget(), || Ok(())).unwrap();
            let decays = step.burnup[0].reaction_counts[0];
            close(step.decay_precursor_counts.as_ref().unwrap()[0][0], 0.5 * decays, 1e-9);
            close(step.decay_neutron_counts.as_ref().unwrap()[0][0], 0.25 * decays, 1e-9);
            close(step.neutrons.captures, step.burnup[0].reaction_counts[1], 1e-8);
            leakage += step.neutrons.leaked_neutrons;
            history_allowance += step.reaction_defect + step.neutrons.population_roundoff
                + step.burnup[0].population_roundoff_nuclei;
            state = next;
        }
        let p = state.inventories()[0].populations();
        assert!(p[2] > 1.0);
        let got = [p[0], state.neutrons().precursor_counts()[0], state.neutrons().neutron_counts()[0], p[1], p[2], leakage];
        close(p.iter().sum(), 100.0, 1e-8);
        // The complete history composes actual accepted defects from every
        // interval; the small final margin covers this assertion's arithmetic.
        close(got[1] + got[2] + p[2] + leakage, 10.0 + 0.75 * p[1], history_allowance + 1e-11);
        errors.push(got.into_iter().zip(expected).map(|(a, b)| (a - b).abs()).sum::<f64>());
    }
    assert!(errors[1] < 0.4 * errors[0], "physical-refinement errors {errors:?}");
    assert!(errors[1] < 2e-2, "physical-refinement errors {errors:?}");
}

#[test]
fn g1_parent_and_precursor_outflow_remove_their_future_local_emission() {
    let original = unbound(&[1.0, 0.0], &[100.0, 0.0], None, 0.0);
    let source = emission(original.inventories()[0].model());
    let mut data = original.neutrons().model().data().clone();
    data.delayed_groups = vec![family()];
    let neutron_model = Arc::new(DiffusionModel::new_with_transport(
        data,
        PrecursorTransport { transfers: vec![], outflows: vec![PrecursorOutflow { cell: 0, volumetric_flow_m3_s: 1.0 }] },
        budget().diffusion.limits,
    ).unwrap());
    let limits = BurnupTransportLimits { burnup: budget().burnup.limits, ..BurnupTransportLimits::default() };
    let material = BurnupTransportState::new(
        original.inventories().to_vec(), vec![1.0],
        NuclideTransport { transfers: vec![], outflows: vec![NuclideOutflow { compartment: 0, rates_per_s: vec![1.0, 0.0] }] },
        limits, || Ok(()),
    ).unwrap();
    let initial = DiffusionBurnupState::new_with_transport(
        DiffusionState::new(neutron_model, vec![0.0], vec![0.0]).unwrap(),
        material, vec![], limits,
    ).unwrap().with_decay_neutrons(vec![source]).unwrap();
    let (next, report) = initial.trial(1.0, &[300.0], budget(), || Ok(())).unwrap();
    let e2 = (-2.0_f64).exp();
    let e3 = (-3.0_f64).exp();
    let decays = 50.0 * (1.0 - e2);
    let precursors = 50.0 * (e2 - e3);
    let precursor_residence = 50.0 * ((1.0 - e2) / 2.0 - (1.0 - e3) / 3.0);
    close(next.inventories()[0].populations()[0], 100.0 * e2, 1e-9);
    close(report.transport.as_ref().unwrap().outflow_counts[0][0], decays, 1e-9);
    close(report.decay_precursor_counts.as_ref().unwrap()[0][0], 0.5 * decays, 1e-9);
    close(next.neutrons().precursor_counts()[0], precursors, 1e-9);
    close(report.neutrons.precursor_outflow_counts[0], precursor_residence, 1e-9);
    close(report.neutrons.precursor_decays, 2.0 * precursor_residence, 1e-9);
    close(next.neutrons().neutron_counts()[0], 0.25 * decays + 2.0 * precursor_residence, 1e-9);
    close(precursors + report.neutrons.precursor_decays + report.neutrons.escaped_precursors, 0.5 * decays, 1e-9);
    assert_eq!(report.neutrons.leaked_neutrons, 0.0);
}

#[test]
fn g0_g4_family_identity_stable_rows_limits_and_cancellation_are_atomic() {
    let initial = fixture(None, 3.0, 0.0);
    let model = initial.inventories()[0].model();
    let make = |families, precursors| BurnupDecayNeutronEmission::new_with_precursors(
        "manufactured-source-admission".into(), Arc::clone(model),
        vec![vec![0.25], vec![0.0], vec![0.0]], families, precursors,
    );
    assert!(make(vec![family()], vec![vec![0.5], vec![1.0], vec![0.0]]).is_err());
    assert!(make(vec![family()], vec![vec![0.5], vec![], vec![0.0]]).is_err());
    assert!(make(vec![family()], vec![vec![f64::NAN], vec![0.0], vec![0.0]]).is_err());
    let mut different = family();
    different.decay_constant_s_inv = 1.0;
    let wrong = Arc::new(make(vec![different], vec![vec![0.5], vec![0.0], vec![0.0]]).unwrap());
    assert!(initial.clone().with_decay_neutrons(vec![wrong]).is_err());
    let mut invalid = family();
    invalid.spectrum = vec![0.5];
    assert!(make(vec![invalid], vec![vec![0.5], vec![0.0], vec![0.0]]).is_err());
    assert!(make(
        vec![family(); fs_la::decay::MAX_DECAY_NODES + 1],
        vec![vec![0.0; fs_la::decay::MAX_DECAY_NODES + 1]; 3],
    ).is_err());
    // Immediate-only declarations remain admissible on a native family model.
    let immediate = Arc::new(BurnupDecayNeutronEmission::new(
        "manufactured-immediate-only".into(), Arc::clone(model), vec![vec![0.25], vec![0.0], vec![0.0]],
    ).unwrap());
    let (_, immediate_step) = initial.clone().with_decay_neutrons(vec![immediate]).unwrap()
        .trial(0.2, &[300.0], budget(), || Ok(())).unwrap();
    assert_eq!(immediate_step.decay_precursor_counts, Some(vec![vec![0.0]]));
    assert_eq!(immediate_step.neutrons.delayed_neutrons, 0.0);

    let (accepted, _) = initial.trial(0.2, &[300.0], budget(), || Ok(())).unwrap();
    let before = accepted.clone();
    let mut tight = budget();
    tight.diffusion.limits.propagation.maximum_edges = 2;
    tight.burnup.limits.propagation.maximum_edges = 2;
    assert!(accepted.trial(0.3, &[300.0], tight, || Ok(())).is_err());
    let mut polls = 0;
    accepted.trial(0.3, &[300.0], budget(), || { polls += 1; Ok(()) }).unwrap();
    let mut seen = 0;
    assert!(accepted.trial(0.3, &[300.0], budget(), || {
        seen += 1;
        if seen == polls - 1 { Err(FissionError::Cancelled) } else { Ok(()) }
    }).is_err());
    assert_eq!(accepted.inventories(), before.inventories());
    assert_eq!(accepted.neutrons().neutron_counts(), before.neutrons().neutron_counts());
    assert_eq!(accepted.neutrons().precursor_counts(), before.neutrons().precursor_counts());
}
