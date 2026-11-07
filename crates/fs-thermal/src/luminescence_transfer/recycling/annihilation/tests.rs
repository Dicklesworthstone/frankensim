use super::super::{EmissionChannel, PhotonAbsorption, RecyclingLimits, SingletTripletDecay};
use super::*;
use fs_la::decay::DecayError;
use fs_math::det;

fn plan(
    model: SingletTripletDecay,
    channels: &[EmissionChannel],
    absorption: &[PhotonAbsorption],
    duration: f64,
) -> SingletTripletRecyclingPlan {
    SingletTripletRecyclingPlan::prepare(
        &[model],
        channels,
        absorption,
        duration,
        1e-12,
        RecyclingLimits::default(),
        || false,
    )
    .unwrap()
}
fn stable(duration: f64) -> SingletTripletRecyclingPlan {
    plan(
        SingletTripletDecay::try_new([3.0, 2.0], [0.0; 2], [0.0; 2], [0.0; 2]).unwrap(),
        &[],
        &[],
        duration,
    )
}
fn channel(state: usize, energy: f64, escape: f64, unresolved: f64) -> EmissionChannel {
    EmissionChannel {
        source_cell: 0,
        source_state: state,
        photon_energy_j: energy,
        photon_fraction: 1.0,
        escape_fraction: escape,
        unresolved_fraction: unresolved,
    }
}
fn controls(relative: f64, maximum: usize) -> AnnihilationLimits {
    AnnihilationLimits {
        absolute_energy_tolerance_j: 1e-12,
        relative_tolerance: relative,
        maximum_substeps: maximum,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn isolated_annihilation_uses_the_exact_reciprocal_law_and_separate_heat() {
    let p = stable(1.0);
    let law = TripletAnnihilation::try_new(0.25, 2.0, 0.5).unwrap();
    let result = p
        .advance_with_annihilation(
            &[[0.0, 8.0]],
            &[[0.0; 2]],
            &[Some(law)],
            controls(1e-12, 1),
            || false,
        )
        .unwrap();
    close(result.linear.populations[0][1], 8.0 / 3.0, 2e-14);
    close(result.linear.populations[0][0], 4.0 / 3.0, 2e-14);
    close(result.annihilation_events[0], 8.0 / 3.0, 2e-14);
    close(result.annihilation_singlets[0], 4.0 / 3.0, 2e-14);
    close(result.annihilation_heat_j[0], 20.0 / 3.0, 2e-14);
    assert_eq!(result.linear.nonradiative_events, [[0.0; 2]]);
    assert_eq!(result.linear.nonradiative_heat_j, [[0.0; 2]]);
    close(
        result.linear.residual_j,
        result.annihilation_heat_j[0],
        2e-14,
    );
    assert!(result.residual_j.abs() < 1e-13);
    assert_eq!(result.substeps, 1);
    assert_eq!(result.max_energy_discrepancy_j, 0.0);
}

#[test]
fn inactive_laws_preserve_the_exact_full_linear_result_and_validate_limits() {
    let model =
        SingletTripletDecay::try_new([2.0, 1.0], [0.7, 0.3], [0.1, 0.2], [0.4, 0.05]).unwrap();
    let p = plan(
        model,
        &[channel(0, 2.0, 1.0, 0.0), channel(1, 1.0, 1.0, 0.0)],
        &[],
        0.5,
    );
    let populations = [[0.2, 3.0]];
    let pumping = [[0.4, 0.7]];
    let exact = p.advance(&populations, &pumping).unwrap();
    for law in [
        None,
        Some(TripletAnnihilation::try_new(0.0, 1.0, 0.5).unwrap()),
    ] {
        let result = p
            .advance_with_annihilation(&populations, &pumping, &[law], controls(1e-8, 1), || false)
            .unwrap();
        assert_eq!(result.linear, exact);
        assert_eq!(result.annihilation_events, [0.0]);
        assert_eq!(result.substeps, 0);
    }
    let invalid = AnnihilationLimits {
        maximum_substeps: 0,
        ..controls(1e-8, 1)
    };
    assert!(matches!(
        p.advance_with_annihilation(&populations, &pumping, &[None], invalid, || false,),
        Err(SingletTripletError::Invalid(_))
    ));
    assert!(matches!(
        p.advance_with_annihilation(&populations, &pumping, &[], controls(1e-8, 1), || false,),
        Err(SingletTripletError::Invalid(_))
    ));
}

#[test]
fn nonlinear_and_linear_loss_converge_to_the_independent_riccati_solution() {
    // Nt' = -gamma*Nt - beta*Nt^2. Integrating the Riccati solution gives
    // integral Nt dt = ln(1 + beta*N0/gamma*(1-exp(-gamma*t)))/beta.
    let gamma = 1.0;
    let beta = 0.4;
    let n0 = 3.0;
    let duration = 0.8;
    let model = SingletTripletDecay::try_new([2.0, 1.0], [0.0, 0.6], [0.0, 0.4], [0.0; 2]).unwrap();
    let p = plan(model, &[channel(1, 1.0, 1.0, 0.0)], &[], duration);
    let law = TripletAnnihilation::try_new(0.4, 2.0, 0.0).unwrap();
    let decay = det::exp(-gamma * duration);
    let denominator = 1.0 + beta * n0 / gamma * (1.0 - decay);
    let endpoint = n0 * decay / denominator;
    let visits = det::ln(denominator) / beta;
    let pairs = 0.5 * (n0 - endpoint - gamma * visits);
    let run = |tolerance| {
        p.advance_with_annihilation(
            &[[0.0, n0]],
            &[[0.0; 2]],
            &[Some(law)],
            controls(tolerance, 1024),
            || false,
        )
        .unwrap()
    };
    let loose = run(1e-3);
    let tight = run(1e-5);
    assert!(
        (tight.linear.populations[0][1] - endpoint).abs()
            < (loose.linear.populations[0][1] - endpoint).abs()
    );
    close(tight.linear.populations[0][1], endpoint, 5e-5);
    close(tight.linear.radiative_events[0][1], 0.6 * visits, 5e-5);
    close(tight.linear.nonradiative_events[0][1], 0.4 * visits, 5e-5);
    close(tight.annihilation_events[0], pairs, 5e-5);
    close(tight.annihilation_heat_j[0], 2.0 * pairs, 5e-5);
    assert!(tight.substeps >= loose.substeps);
    assert_eq!(tight.linear.crossing_events[0][1], 0.0);
    assert_eq!(tight.linear.reverse_crossing_heat_absorbed_j[0], 0.0);
}

#[test]
fn pumped_closed_photon_traffic_is_resolved_as_an_integrated_energy_channel() {
    // With unit triplet density, pump=1 and beta=1, the exact solution is
    // Nt(t)=1. Self-return cancels radiative removal but still produces 100
    // photons/s. Endpoint-only absolute control cannot establish that history.
    let model = SingletTripletDecay::try_new([2.0, 1.0], [0.0, 100.0], [0.0; 2], [0.0; 2]).unwrap();
    let p = plan(
        model,
        &[channel(1, 1.0, 0.0, 0.0)],
        &[PhotonAbsorption {
            channel: 0,
            receiver_cell: Some(0),
            absorbed_fraction: 1.0,
            excitation_fractions: [0.0, 1.0],
        }],
        0.5,
    );
    let law = TripletAnnihilation::try_new(0.5, 1.0, 0.0).unwrap();
    let result = p
        .advance_with_annihilation(
            &[[0.0, 1.0]],
            &[[0.0, 1.0]],
            &[Some(law)],
            AnnihilationLimits {
                absolute_energy_tolerance_j: 1e-5,
                relative_tolerance: 0.0,
                maximum_substeps: 512,
            },
            || false,
        )
        .unwrap();
    close(result.linear.populations[0][1], 1.0, 2e-5);
    assert!((result.linear.emitted_energy_j[0] - 50.0).abs() < 2e-5);
    assert!((result.linear.absorbed_energy_j[0] - 50.0).abs() < 2e-5);
    close(result.annihilation_events[0], 0.25, 2e-5);
    close(result.annihilation_heat_j[0], 0.5, 2e-5);
    assert!(result.max_energy_discrepancy_j <= 1e-5);
    assert_eq!(result.linear.nonradiative_heat_j, [[0.0; 2]]);
    assert!(result.residual_j.abs() < 2e-10);
}

#[test]
fn nanosecond_singlet_return_recovers_the_correct_slow_annihilation_coefficient() {
    // tauS=1 ns, fISC=.75, phi=.8. In the fast-singlet limit:
    // Nt' = -(2-phi*fISC)*k*Nt^2/V. Singlet fluorescence counts the remaining
    // (1-fISC)*phi fraction of pair events. No ns-sized integration grid is used.
    let fast_rate = 1e9;
    let return_fraction = 0.75;
    let yield_ = 0.8;
    let k = 0.1;
    let n0 = 2.0;
    let duration = 0.5;
    let model = SingletTripletDecay::try_new(
        [1.5, 1.0],
        [(1.0 - return_fraction) * fast_rate, 0.0],
        [0.0; 2],
        [return_fraction * fast_rate, 0.0],
    )
    .unwrap();
    let p = plan(model, &[channel(0, 1.5, 1.0, 0.0)], &[], duration);
    let loss_per_pair = 2.0 - yield_ * return_fraction;
    let triplets = n0 / (1.0 + loss_per_pair * k * n0 * duration);
    let pairs = (n0 - triplets) / loss_per_pair;
    let result = p
        .advance_with_annihilation(
            &[[0.0, n0]],
            &[[0.0; 2]],
            &[Some(TripletAnnihilation::try_new(k, 1.0, yield_).unwrap())],
            AnnihilationLimits {
                absolute_energy_tolerance_j: 1e-9,
                relative_tolerance: 1e-5,
                maximum_substeps: 512,
            },
            || false,
        )
        .unwrap();
    close(result.linear.populations[0][1], triplets, 5e-5);
    close(result.annihilation_events[0], pairs, 5e-5);
    close(
        result.linear.radiative_events[0][0],
        (1.0 - return_fraction) * yield_ * pairs,
        5e-5,
    );
    let singlets_at_endpoint = yield_ * k * triplets * triplets / fast_rate;
    close(result.linear.populations[0][0], singlets_at_endpoint, 8e-5);
    let photon_power = result.linear.stored_energy_j[0][0] * (1.0 - return_fraction) * fast_rate;
    close(
        photon_power,
        1.5 * (1.0 - return_fraction) * yield_ * k * triplets * triplets,
        8e-5,
    );
    assert_eq!(result.linear.nonradiative_events, [[0.0; 2]]);
    assert_eq!(result.linear.nonradiative_heat_j, [[0.0; 2]]);
    assert_eq!(result.linear.crossing_events[0][1], 0.0);
    assert_eq!(result.linear.reverse_crossing_heat_absorbed_j[0], 0.0);
    assert!(result.substeps <= 512);
}

#[test]
fn pumping_crossing_and_photon_return_close_the_independent_external_energy_account() {
    let model =
        SingletTripletDecay::try_new([1.5, 1.0], [0.8, 0.3], [0.1, 0.2], [0.5, 0.1]).unwrap();
    let p = plan(
        model,
        &[channel(0, 1.4, 0.2, 0.1), channel(1, 1.0, 0.3, 0.1)],
        &[
            PhotonAbsorption {
                channel: 0,
                receiver_cell: Some(0),
                absorbed_fraction: 0.7,
                excitation_fractions: [0.0, 0.7],
            },
            PhotonAbsorption {
                channel: 1,
                receiver_cell: Some(0),
                absorbed_fraction: 0.6,
                excitation_fractions: [0.0, 0.6],
            },
        ],
        0.4,
    );
    let result = p
        .advance_with_annihilation(
            &[[0.4, 2.0]],
            &[[0.3, 0.5]],
            &[Some(TripletAnnihilation::try_new(0.2, 1.0, 0.5).unwrap())],
            controls(1e-5, 512),
            || false,
        )
        .unwrap();
    let linear = &result.linear;
    let input = 0.4 * 1.5 + 2.0 + 0.4 * (0.3 * 1.5 + 0.5);
    let stored: f64 = linear.stored_energy_j.iter().flatten().sum();
    let outgoing: f64 = linear
        .escaped_energy_j
        .iter()
        .chain(&linear.unresolved_energy_j)
        .sum();
    let heat: f64 = linear
        .nonradiative_heat_j
        .iter()
        .flatten()
        .chain(linear.relaxation_heat_j.iter().flatten())
        .chain(&linear.crossing_heat_j)
        .chain(&linear.absorption_heat_j)
        .chain(&result.annihilation_heat_j)
        .sum::<f64>()
        - linear.reverse_crossing_heat_absorbed_j.iter().sum::<f64>();
    assert!((input - stored - outgoing - heat).abs() < 3e-10);
    assert!(result.residual_j.abs() < 3e-10);
    close(linear.residual_j, result.annihilation_heat_j[0], 3e-9);
    assert!(linear.crossing_events[0].iter().all(|&count| count > 0.0));
    assert!(linear.relaxation_heat_j[0][0] > 0.0);
    assert!(linear.recycled_excitations[0][1] > 0.0);
}

#[test]
fn budget_and_cancellation_refuse_unresolved_coupling() {
    let model = SingletTripletDecay::try_new([3.0, 2.0], [0.0; 2], [0.0, 1.0], [0.0; 2]).unwrap();
    let p = plan(model, &[], &[], 1.0);
    let law = TripletAnnihilation::try_new(0.25, 1.0, 0.5).unwrap();
    assert!(matches!(
        p.advance_with_annihilation(
            &[[0.0, 2.0]],
            &[[0.0; 2]],
            &[Some(law)],
            controls(1e-8, 1),
            || false,
        ),
        Err(SingletTripletError::Evolution(DecayError::Budget {
            required: 2,
            maximum: 1,
            ..
        }))
    ));
    let mut polls = 0;
    assert!(matches!(
        p.advance_with_annihilation(
            &[[0.0, 2.0]],
            &[[0.0; 2]],
            &[Some(law)],
            controls(1e-5, 512),
            || {
                polls += 1;
                polls >= 8
            },
        ),
        Err(SingletTripletError::Evolution(DecayError::Cancelled))
    ));
    // An over-damped frozen predictor must not create false coarse/fine
    // agreement with the no-annihilation trajectory.
    let strong = TripletAnnihilation::try_new(1000.0, 1.0, 0.5).unwrap();
    assert!(matches!(
        p.advance_with_annihilation(
            &[[0.0, 2.0]],
            &[[0.0; 2]],
            &[Some(strong)],
            controls(1e-5, 8),
            || false,
        ),
        Err(SingletTripletError::Evolution(DecayError::Budget { .. }))
    ));
}

#[test]
fn unrepresentable_split_time_frozen_rate_and_positive_channels_are_refused() {
    let least = f64::from_bits(1);
    let law = TripletAnnihilation::try_new(1.0, 1.0, 0.5).unwrap();
    assert!(matches!(
        stable(least).advance_with_annihilation(
            &[[0.0, 1.0]],
            &[[0.0, 1.0]],
            &[Some(law)],
            controls(1e-5, 16),
            || false,
        ),
        Err(SingletTripletError::Unrepresentable(_))
    ));
    let huge = TripletAnnihilation::try_new(f64::MAX, 1.0, 0.5).unwrap();
    assert!(matches!(
        stable(1.0).advance_with_annihilation(
            &[[0.0, 2.0]],
            &[[0.0, 1.0]],
            &[Some(huge)],
            controls(1e-5, 16),
            || false,
        ),
        Err(SingletTripletError::Unrepresentable(_))
    ));
    let tiny = TripletAnnihilation::try_new(least, 1.0, 0.5).unwrap();
    assert!(matches!(
        stable(1.0).advance_with_annihilation(
            &[[0.0, least]],
            &[[0.0, 1.0]],
            &[Some(tiny)],
            controls(1e-5, 16),
            || false,
        ),
        Err(SingletTripletError::Unrepresentable(_))
    ));
}

#[test]
fn loose_nonlinear_accuracy_cannot_relax_the_external_energy_audit() {
    let p = stable(1.0);
    let initial = [[0.0, 8.0]];
    let law = TripletAnnihilation::try_new(0.25, 2.0, 0.5).unwrap();
    let mut result = p
        .advance_with_annihilation(
            &initial,
            &[[0.0; 2]],
            &[Some(law)],
            AnnihilationLimits {
                absolute_energy_tolerance_j: 1e20,
                relative_tolerance: 0.0,
                maximum_substeps: 1,
            },
            || false,
        )
        .unwrap();
    // A future reconstruction bug that invents one joule must fail the
    // independent physical audit even when the integration request is loose.
    result.annihilation_heat_j[0] += 1.0;
    assert!(matches!(
        balance(&p, &initial, &mut result),
        Err(SingletTripletError::EnergyBalance { .. })
    ));
}
