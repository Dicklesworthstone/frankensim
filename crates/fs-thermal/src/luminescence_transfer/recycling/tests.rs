use super::*;

fn model(
    energies: [f64; 2],
    radiative: [f64; 2],
    nonradiative: [f64; 2],
    crossing: [f64; 2],
) -> SingletTripletDecay {
    SingletTripletDecay::try_new(energies, radiative, nonradiative, crossing).unwrap()
}
fn photon(cell: usize, state: usize, energy: f64, escape: f64, unresolved: f64) -> EmissionChannel {
    EmissionChannel {
        source_cell: cell,
        source_state: state,
        photon_energy_j: energy,
        photon_fraction: 1.0,
        escape_fraction: escape,
        unresolved_fraction: unresolved,
    }
}
fn absorption(
    channel: usize,
    receiver: Option<usize>,
    fraction: f64,
    activation: [f64; 2],
) -> PhotonAbsorption {
    PhotonAbsorption {
        channel,
        receiver_cell: receiver,
        absorbed_fraction: fraction,
        excitation_fractions: activation,
    }
}
fn prepare(
    models: &[SingletTripletDecay],
    channels: &[EmissionChannel],
    links: &[PhotonAbsorption],
    duration_s: f64,
) -> SingletTripletRecyclingPlan {
    SingletTripletRecyclingPlan::prepare(
        models,
        channels,
        links,
        duration_s,
        32.0 * f64::EPSILON,
        RecyclingLimits::default(),
        || false,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-11 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}
fn balanced(step: &SingletTripletRecyclingStep, scale: f64) {
    assert!(step.residual_j.abs() <= 4e-11 * scale, "{step:?}");
}

#[test]
fn g0_unequal_energy_recycling_obeys_independent_bateman_counts_and_stokes_heat() {
    let a = 2.0;
    let b = 0.5;
    let t: f64 = 1.3;
    let models = [
        model([5e-19, 4e-19], [a, 0.0], [0.0; 2], [0.0; 2]),
        model([4e-19, 3e-19], [0.0, b], [0.0; 2], [0.0; 2]),
    ];
    let channels = [
        photon(0, 0, 4.5e-19, 0.0, 0.0),
        photon(1, 1, 2.5e-19, 1.0, 0.0),
    ];
    let links = [absorption(0, Some(1), 1.0, [0.0, 1.0])];
    let plan = prepare(&models, &channels, &links, t);
    let step = plan
        .advance(&[[1.0, 0.0], [0.0; 2]], &[[0.0; 2]; 2])
        .unwrap();
    let source = (-a * t).exp();
    let receiver = a * ((-b * t).exp() - source) / (a - b);
    let emitted_source = 1.0 - source;
    let emitted_receiver = emitted_source - receiver;
    close(step.populations[0][0], source);
    close(step.populations[1][1], receiver);
    close(step.radiative_events[0][0], emitted_source);
    close(step.radiative_events[1][1], emitted_receiver);
    close(step.recycled_excitations[1][1], emitted_source);
    close(step.recycled_energy_j[1][1], emitted_source * 3e-19);
    close(step.emitted_energy_j[0], emitted_source * 4.5e-19);
    close(step.absorbed_energy_j[0], emitted_source * 4.5e-19);
    close(step.relaxation_heat_j[0][0], emitted_source * 0.5e-19);
    close(step.absorption_heat_j[0], emitted_source * 1.5e-19);
    close(step.relaxation_heat_j[1][1], emitted_receiver * 0.5e-19);
    close(step.escaped_energy_j[1], emitted_receiver * 2.5e-19);
    balanced(&step, 5e-19);

    // Independent constant-injection solution, including the actual ages of
    // photons absorbed during this interval instead of an endpoint reinjection.
    let pump = 3.0;
    let charged = plan
        .advance(&[[0.0; 2]; 2], &[[pump, 0.0], [0.0; 2]])
        .unwrap();
    let ns = pump * (1.0 - source) / a;
    let nt = pump * ((1.0 - (-b * t).exp()) / b - ((-b * t).exp() - source) / (a - b));
    close(charged.populations[0][0], ns);
    close(charged.populations[1][1], nt);
    close(charged.radiative_events[0][0], pump * t - ns);
    close(charged.radiative_events[1][1], pump * t - ns - nt);
    balanced(&charged, pump * t * 5e-19);
}

#[test]
fn g0_closed_resonant_photon_cycles_retain_population_and_all_reemissions() {
    let t: f64 = 4.0;
    let m = model([4e-19, 3e-19], [2.0, 0.0], [0.0; 2], [0.0; 2]);
    let plan = prepare(
        &[m],
        &[photon(0, 0, 4e-19, 0.0, 0.0)],
        &[absorption(0, Some(0), 1.0, [1.0, 0.0])],
        t,
    );
    let step = plan.advance(&[[7.0, 0.0]], &[[2.0, 0.0]]).unwrap();
    let events = 2.0 * (7.0 * t + 0.5 * 2.0 * t * t);
    close(step.populations[0][0], 7.0 + 2.0 * t);
    close(step.radiative_events[0][0], events);
    close(step.recycled_excitations[0][0], events);
    close(step.emitted_energy_j[0], events * 4e-19);
    close(step.absorbed_energy_j[0], events * 4e-19);
    assert!(step.emitted_energy_j[0] > step.stored_energy_j[0][0]);
    assert_eq!(step.absorption_heat_j, [0.0]);
    assert_eq!(step.relaxation_heat_j, [[0.0; 2]]);
    assert_eq!(step.escaped_energy_j, [0.0]);
    assert_eq!(step.unresolved_energy_j, [0.0]);
    balanced(&step, 15.0 * 4e-19);
}

#[test]
fn g0_intersystem_crossing_and_triplet_photon_return_share_one_history() {
    let t: f64 = 1.2;
    let m = model([5e-19, 4e-19], [0.0, 3.0], [0.0; 2], [2.0, 0.0]);
    let plan = prepare(
        &[m],
        &[photon(0, 1, 4e-19, 0.0, 0.0)],
        &[absorption(0, Some(0), 1.0, [0.0, 1.0])],
        t,
    );
    let step = plan.advance(&[[1.0, 0.0]], &[[0.0; 2]]).unwrap();
    let singlet = (-2.0 * t).exp();
    let triplet = 1.0 - singlet;
    let emitted = 3.0 * (t - triplet / 2.0);
    close(step.populations[0][0], singlet);
    close(step.populations[0][1], triplet);
    close(step.crossing_events[0][0], triplet);
    close(step.radiative_events[0][1], emitted);
    close(step.recycled_excitations[0][1], emitted);
    close(step.crossing_heat_j[0], triplet * 1e-19);
    assert_eq!(step.reverse_crossing_heat_absorbed_j, [0.0]);
    assert_eq!(step.escaped_energy_j, [0.0]);
    balanced(&step, 5e-19);
}

#[test]
fn g0_reverse_crossing_pays_for_recycled_delayed_singlet_photons() {
    let t: f64 = 0.9;
    let m = model([4e-19, 3e-19], [1.0, 0.0], [0.0; 2], [0.0, 2.0]);
    let plan = prepare(
        &[m],
        &[photon(0, 0, 4e-19, 0.5, 0.0)],
        &[absorption(0, Some(0), 0.5, [1.0, 0.0])],
        t,
    );
    let step = plan.advance(&[[0.0, 1.0]], &[[0.0; 2]]).unwrap();
    let triplet = (-2.0 * t).exp();
    let singlet = (2.0 / 1.5) * ((-0.5 * t).exp() - triplet);
    let emitted = (2.0 / 1.5) * ((1.0 - (-0.5 * t).exp()) / 0.5 - (1.0 - triplet) / 2.0);
    close(step.populations[0][0], singlet);
    close(step.populations[0][1], triplet);
    close(step.radiative_events[0][0], emitted);
    close(step.recycled_excitations[0][0], 0.5 * emitted);
    close(step.escaped_energy_j[0], 0.5 * emitted * 4e-19);
    close(
        step.reverse_crossing_heat_absorbed_j[0],
        (1.0 - triplet) * 1e-19,
    );
    assert_eq!(step.crossing_heat_j, [0.0]);
    balanced(&step, 4e-19);
}

#[test]
fn g0_activation_branching_passive_heat_escape_and_unresolved_remain_distinct() {
    let models = [
        model([5e-19, 4e-19], [1.0, 0.0], [0.0; 2], [0.0; 2]),
        model([3e-19, 2e-19], [0.0; 2], [0.0; 2], [0.0; 2]),
    ];
    let plan = prepare(
        &models,
        &[photon(0, 0, 4e-19, 0.2, 0.1)],
        &[
            absorption(0, Some(1), 0.5, [0.2, 0.3]),
            absorption(0, None, 0.2, [0.0; 2]),
        ],
        1.0,
    );
    let step = plan
        .advance(&[[1.0, 0.0], [0.0; 2]], &[[0.0; 2]; 2])
        .unwrap();
    let released = 1.0 - (-1.0_f64).exp();
    close(step.populations[1][0], 0.1 * released);
    close(step.populations[1][1], 0.15 * released);
    close(
        step.recycled_energy_j[1][0] + step.recycled_energy_j[1][1],
        0.6e-19 * released,
    );
    close(step.absorbed_energy_j[0], 2e-19 * released);
    close(step.absorption_heat_j[0], 1.4e-19 * released);
    close(step.absorption_heat_j[1], 0.8e-19 * released);
    close(step.relaxation_heat_j[0][0], 1e-19 * released);
    close(step.escaped_energy_j[0], 0.8e-19 * released);
    close(step.unresolved_energy_j[0], 0.4e-19 * released);
    assert_eq!(step.nonradiative_heat_j, [[0.0; 2]; 2]);
    balanced(&step, 5e-19);
}

#[test]
fn g1_without_return_the_existing_native_singlet_triplet_solution_is_preserved() {
    let m = model([4e-19, 3e-19], [2.0, 0.1], [1.0, 0.2], [3.0, 0.7]);
    let t = 0.8;
    let old = [4.0, 2.0];
    let pump = [1.0, 0.3];
    let independent = m
        .prepare(t, RecyclingLimits::default(), || false)
        .unwrap()
        .advance(old, pump)
        .unwrap();
    let plan = prepare(
        &[m],
        &[photon(0, 0, 4e-19, 1.0, 0.0), photon(0, 1, 3e-19, 1.0, 0.0)],
        &[],
        t,
    );
    let step = plan.advance(&[old], &[pump]).unwrap();
    for state in 0..2 {
        close(step.populations[0][state], independent.populations[state]);
        close(
            step.radiative_events[0][state],
            independent.radiative_events[state],
        );
        close(
            step.nonradiative_events[0][state],
            independent.nonradiative_events[state],
        );
        close(
            step.crossing_events[0][state],
            independent.crossing_events[state],
        );
        close(
            step.stored_energy_j[0][state],
            independent.stored_energy_j[state],
        );
        close(
            step.escaped_energy_j[state],
            independent.radiative_decay_energy_j[state],
        );
    }
    close(step.crossing_heat_j[0], independent.crossing_heat_j);
    close(
        step.reverse_crossing_heat_absorbed_j[0],
        independent.reverse_crossing_heat_absorbed_j,
    );
    balanced(&step, 3e-18);
}

#[test]
fn g3_reprepared_split_intervals_preserve_stores_and_every_integrated_channel() {
    let m = model([4e-19, 3e-19], [1.0, 0.5], [0.3, 0.2], [2.0, 0.7]);
    let channels = [
        photon(0, 0, 3.5e-19, 0.2, 0.1),
        photon(0, 1, 2.5e-19, 0.3, 0.2),
    ];
    let links = [
        absorption(0, Some(0), 0.5, [0.0, 0.8]),
        absorption(0, None, 0.2, [0.0; 2]),
        absorption(1, None, 0.5, [0.0; 2]),
    ];
    let whole = prepare(&[m], &channels, &links, 1.25);
    let first = whole.with_duration(0.4, || false).unwrap();
    let second = whole.with_duration(0.85, || false).unwrap();
    let all = whole.advance(&[[2.0, 1.0]], &[[0.8, 0.3]]).unwrap();
    let a = first.advance(&[[2.0, 1.0]], &[[0.8, 0.3]]).unwrap();
    let b = second.advance(&a.populations, &[[0.8, 0.3]]).unwrap();
    for state in 0..2 {
        close(all.populations[0][state], b.populations[0][state]);
        close(all.stored_energy_j[0][state], b.stored_energy_j[0][state]);
    }
    let all_pairs = [
        &all.injected_excitations,
        &all.radiative_events,
        &all.nonradiative_events,
        &all.crossing_events,
        &all.pumped_energy_j,
        &all.radiative_decay_energy_j,
        &all.nonradiative_heat_j,
        &all.relaxation_heat_j,
        &all.recycled_excitations,
        &all.recycled_energy_j,
    ];
    let first_pairs = [
        &a.injected_excitations,
        &a.radiative_events,
        &a.nonradiative_events,
        &a.crossing_events,
        &a.pumped_energy_j,
        &a.radiative_decay_energy_j,
        &a.nonradiative_heat_j,
        &a.relaxation_heat_j,
        &a.recycled_excitations,
        &a.recycled_energy_j,
    ];
    let second_pairs = [
        &b.injected_excitations,
        &b.radiative_events,
        &b.nonradiative_events,
        &b.crossing_events,
        &b.pumped_energy_j,
        &b.radiative_decay_energy_j,
        &b.nonradiative_heat_j,
        &b.relaxation_heat_j,
        &b.recycled_excitations,
        &b.recycled_energy_j,
    ];
    for ((all, a), b) in all_pairs.into_iter().zip(first_pairs).zip(second_pairs) {
        for ((all, a), b) in all
            .iter()
            .flatten()
            .zip(a.iter().flatten())
            .zip(b.iter().flatten())
        {
            close(*all, *a + *b);
        }
    }
    let all_fields = [
        &all.crossing_heat_j,
        &all.reverse_crossing_heat_absorbed_j,
        &all.emitted_energy_j,
        &all.escaped_energy_j,
        &all.unresolved_energy_j,
        &all.absorbed_energy_j,
        &all.absorption_heat_j,
    ];
    let first_fields = [
        &a.crossing_heat_j,
        &a.reverse_crossing_heat_absorbed_j,
        &a.emitted_energy_j,
        &a.escaped_energy_j,
        &a.unresolved_energy_j,
        &a.absorbed_energy_j,
        &a.absorption_heat_j,
    ];
    let second_fields = [
        &b.crossing_heat_j,
        &b.reverse_crossing_heat_absorbed_j,
        &b.emitted_energy_j,
        &b.escaped_energy_j,
        &b.unresolved_energy_j,
        &b.absorbed_energy_j,
        &b.absorption_heat_j,
    ];
    for ((all, a), b) in all_fields.into_iter().zip(first_fields).zip(second_fields) {
        for ((all, a), b) in all.iter().zip(a).zip(b) {
            close(*all, *a + *b);
        }
    }
    balanced(&all, 2e-18);
    balanced(&a, 2e-18);
    balanced(&b, 2e-18);
}

#[test]
fn g0_energy_threshold_multiplicity_and_missing_transport_are_never_inferred() {
    let m = model([4e-19, 3e-19], [1.0, 0.0], [0.0; 2], [0.0; 2]);
    let admission = |channels: &[EmissionChannel], links: &[PhotonAbsorption]| {
        SingletTripletRecyclingPlan::prepare(
            &[m],
            channels,
            links,
            1.0,
            32.0 * f64::EPSILON,
            RecyclingLimits::default(),
            || false,
        )
    };
    assert!(admission(&[photon(0, 0, 5e-19, 1.0, 0.0)], &[]).is_err());
    assert!(
        admission(
            &[photon(0, 0, 3.5e-19, 0.0, 0.0)],
            &[absorption(0, Some(0), 1.0, [1.0, 0.0])]
        )
        .is_err()
    );
    assert!(
        admission(
            &[photon(0, 0, 4e-19, 0.0, 0.0)],
            &[absorption(0, Some(0), 1.0, [1.0, 1.0])]
        )
        .is_err()
    );
    assert!(
        admission(
            &[photon(0, 0, 4e-19, 0.0, 0.0)],
            &[absorption(0, None, 1.0, [0.0, 1.0])]
        )
        .is_err()
    );
    assert!(admission(&[photon(0, 0, 4e-19, 0.8, 0.0)], &[]).is_err());
    assert!(admission(&[], &[]).is_err());

    // A zero-weight blue bin in a shared spectrum is not an energetic photon.
    let mut inactive = photon(0, 0, 5e-19, 1.0, 0.0);
    inactive.photon_fraction = 0.0;
    let valid = admission(&[photon(0, 0, 4e-19, 1.0, 0.0), inactive], &[]).unwrap();
    let step = valid.advance(&[[1.0, 0.0]], &[[0.0; 2]]).unwrap();
    assert_eq!(step.emitted_energy_j[1], 0.0);
    assert_eq!(step.escaped_energy_j[1], 0.0);

    // An explicitly stable state needs neither a fake lifetime nor a spectrum.
    let stable = model([4e-19, 3e-19], [0.0; 2], [0.0; 2], [0.0; 2]);
    let unchanged = prepare(&[stable], &[], &[], 2.0)
        .advance(&[[1.0, 2.0]], &[[0.5, 1.0]])
        .unwrap();
    assert_eq!(unchanged.populations, [[2.0, 4.0]]);
}

#[test]
fn g0_positive_external_loss_cannot_hide_behind_closed_photon_traffic() {
    let m = model([4e-19, 3e-19], [1.0, 0.0], [0.0; 2], [0.0; 2]);
    let result = SingletTripletRecyclingPlan::prepare(
        &[m],
        &[photon(0, 0, 4e-19, 1e-20, 0.0)],
        &[absorption(0, Some(0), 1.0, [1.0, 0.0])],
        1.0,
        1e-14,
        RecyclingLimits::default(),
        || false,
    );
    assert!(matches!(
        result,
        Err(SingletTripletError::Unrepresentable(
            "positive excitation loss erased by photon recycling"
        ))
    ));

    let models = [m];
    let plan = prepare(
        &models,
        &[photon(0, 0, 4e-19, 0.0, 0.0)],
        &[absorption(0, Some(0), 1.0, [1.0, 0.0])],
        1.0,
    );
    let mut step = plan.advance(&[[1.0, 0.0]], &[[0.0; 2]]).unwrap();
    step.emitted_energy_j[0] = 1e20;
    step.absorbed_energy_j[0] = 1e20;
    step.recycled_energy_j[0][0] = 1e20;
    step.escaped_energy_j[0] = 1e-19;
    assert!(matches!(
        admit_energy_balance(&[[1.0, 0.0]], &models, &step, 2e-12),
        Err(SingletTripletError::EnergyBalance { .. })
    ));
}

#[test]
fn g0_zero_time_invalid_shapes_budgets_and_cancellation_are_transactional() {
    let m = model([4e-19, 3e-19], [1.0, 0.0], [0.0; 2], [0.0; 2]);
    let channels = [photon(0, 0, 4e-19, 1.0, 0.0)];
    let zero = prepare(&[m], &channels, &[], 0.0)
        .advance(&[[2.0, 3.0]], &[[f64::MAX; 2]])
        .unwrap();
    assert_eq!(zero.populations, [[2.0, 3.0]]);
    assert_eq!(zero.injected_excitations, [[0.0; 2]]);
    assert_eq!(zero.emitted_energy_j, [0.0]);
    let plan = prepare(&[m], &channels, &[], 1.0);
    assert!(plan.advance(&[], &[[0.0; 2]]).is_err());
    assert!(plan.advance(&[[1.0, -1.0]], &[[0.0; 2]]).is_err());
    assert!(
        plan.advance(&[[1.0, 0.0]], &[[f64::INFINITY, 0.0]])
            .is_err()
    );
    assert!(matches!(
        plan.advance_with_cancel(&[[1.0, 0.0]], &[[0.0; 2]], || true),
        Err(SingletTripletError::Evolution(DecayError::Cancelled))
    ));
    assert!(matches!(
        SingletTripletRecyclingPlan::prepare(
            &[m],
            &channels,
            &[],
            1.0,
            1e-14,
            RecyclingLimits::default(),
            || true
        ),
        Err(SingletTripletError::Evolution(DecayError::Cancelled))
    ));
    let too_many = vec![m; MAX_DECAY_NODES / 2 + 1];
    assert!(matches!(
        SingletTripletRecyclingPlan::prepare(
            &too_many,
            &[],
            &[],
            1.0,
            1e-14,
            RecyclingLimits::default(),
            || false
        ),
        Err(SingletTripletError::Evolution(DecayError::Budget { .. }))
    ));
    let limits = RecyclingLimits {
        maximum_work: 1,
        ..RecyclingLimits::default()
    };
    assert!(matches!(
        SingletTripletRecyclingPlan::prepare(&[m], &channels, &[], 1.0, 1e-14, limits, || false),
        Err(SingletTripletError::Evolution(DecayError::Budget { .. }))
    ));
}

#[test]
fn g0_closed_bath_photon_cycles_cannot_inflate_the_external_energy_allowance() {
    let models = [model([4e-19, 3e-19], [1.0, 0.0], [0.0; 2], [0.0, 2.0])];
    let plan = prepare(
        &models,
        &[photon(0, 0, 4e-19, 0.0, 0.0)],
        &[absorption(0, Some(0), 1.0, [0.0, 1.0])],
        1.0,
    );
    // At this kinetic equilibrium every bath-funded T->S event is followed by
    // an S photon whose reabsorption returns the same energy gap as heat.
    let initial = [[2.0 / 3.0, 1.0 / 3.0]];
    let mut step = plan.advance(&initial, &[[0.0; 2]]).unwrap();
    close(step.populations[0][0], initial[0][0]);
    close(step.populations[0][1], initial[0][1]);
    close(
        step.absorption_heat_j[0],
        step.reverse_crossing_heat_absorbed_j[0],
    );
    balanced(&step, 4e-19);

    // Enormous internally balanced heat traffic is not permission to invent an
    // external photon. The external scale remains the finite physical stores.
    step.absorption_heat_j[0] = 1e20;
    step.reverse_crossing_heat_absorbed_j[0] = 1e20;
    step.escaped_energy_j[0] = 1e-19;
    assert!(matches!(
        admit_energy_balance(&initial, &models, &step, 2e-12),
        Err(SingletTripletError::EnergyBalance { .. })
    ));
}
