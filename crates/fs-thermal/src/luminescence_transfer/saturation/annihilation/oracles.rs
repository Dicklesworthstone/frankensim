//! Independent finite-inventory annihilation and delayed-fluorescence oracles.

use crate::luminescence::transfer::{
    SingletTripletDecay,
    saturation::{
        SaturableSingletTripletPlan,
        annihilation::{AnnihilationLimits, TripletAnnihilation},
    },
};
use fs_la::decay::DecayLimits;

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-290),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn controls(relative_tolerance: f64, maximum_substeps: usize) -> AnnihilationLimits {
    AnnihilationLimits {
        absolute_energy_tolerance_j: 0.0,
        relative_tolerance,
        maximum_substeps,
    }
}

#[test]
fn g1_isolated_pairs_return_the_original_ground_inventory_and_physical_pair_heat() {
    let energies = [4e-19, 3e-19];
    let initial = [2e12, 1e12, 8e12];
    let model = SingletTripletDecay::try_new(energies, [0.0; 2], [0.0; 2], [0.0; 2]).unwrap();
    let plan = SaturableSingletTripletPlan::prepare(
        model,
        [0.0; 2],
        [0.0; 2],
        1.0,
        DecayLimits::default(),
        || false,
    )
    .unwrap();

    // a=k/V=1.25e-13 per emitter per second and T0=8e12 give
    // T(1)=T0/(1+2*a*T0)=T0/3, A=(T0-T)/2=8e12/3.
    // Each productive event returns one emitter to ground; an unsuccessful
    // event returns both. The original singlet store is unaffected.
    let pairs = 8e12 / 3.0;
    for yield_ in [0.0, 0.5, 1.0] {
        let law = TripletAnnihilation::try_new(0.25e-18, 2e-6, yield_).unwrap();
        let result = plan
            .advance_with_annihilation(initial, Some(law), controls(1e-10, 1), || false)
            .unwrap();
        let step = &result.kinetics;
        close(
            step.populations[0],
            initial[0] + (2.0 - yield_) * pairs,
            3e-14,
        );
        close(step.populations[1], initial[1] + yield_ * pairs, 3e-14);
        close(step.populations[2], 8e12 / 3.0, 3e-14);
        close(step.populations.iter().sum(), 11e12, 3e-14);
        close(step.annihilation_events, pairs, 3e-14);
        close(step.annihilation_singlets, yield_ * pairs, 3e-14);
        close(
            step.annihilation_heat_j,
            pairs * (2.0 * energies[1] - yield_ * energies[0]),
            3e-14,
        );
        assert_eq!(step.pump_events, [0.0; 2]);
        assert_eq!(step.radiative_events, [0.0; 2]);
        assert_eq!(step.nonradiative_events, [0.0; 2]);
        assert_eq!(step.nonradiative_heat_j, [0.0; 2]);
        assert_eq!(step.crossing_events, [0.0; 2]);
        assert_eq!(step.reverse_crossing_heat_absorbed_j, 0.0);
        assert!(step.residual_j.abs() < 2e-18);
        assert_eq!(result.substeps, 1);
    }
    assert_eq!(initial, [2e12, 1e12, 8e12]);
}

#[test]
fn g3_nanosecond_delayed_fluorescence_matches_the_rational_source_convolution() {
    let count = 1e12;
    let singlet_rate = 1e9;
    let model =
        SingletTripletDecay::try_new([4e-19, 3e-19], [singlet_rate, 0.0], [0.0; 2], [0.0; 2])
            .unwrap();
    let plan = SaturableSingletTripletPlan::prepare(
        model,
        [0.0; 2],
        [0.0; 2],
        1.0,
        DecayLimits::default(),
        || false,
    )
    .unwrap();
    let law = TripletAnnihilation::try_new(1e-18, 1e-6, 0.8).unwrap();
    let result = plan
        .advance_with_annihilation([0.0, 0.0, count], Some(law), controls(1e-5, 1024), || false)
        .unwrap();
    let step = &result.kinetics;

    // T=T0/q and A=(T0-T)/2, q=1+2*a*T0*t. The original S lifetime
    // convolves the production rate: S=phi*integral(a*T(s)^2*exp(-ks*(t-s)),ds).
    // With c=ks/(2*a*T0), the exact result is phi*T0/2 times
    // [c*exp(-c*q)*(Ei(c*q)-Ei(c)) - 1/q + exp(-ks*t)].
    // These constants were independently evaluated at 85 decimal digits.
    let triplets = 333_333_333_333.333_3;
    let singlets = 88.888_889_007_407_41;
    let pairs = 333_333_333_333.333_3;
    let fluorescence = 266_666_666_577.777_77;
    close(step.populations[2], triplets, 3e-5);
    close(step.populations[1], singlets, 3e-5);
    close(step.populations[0], count - triplets - singlets, 3e-5);
    close(step.populations.iter().sum(), count, 2e-12);
    close(step.annihilation_events, pairs, 3e-5);
    close(step.annihilation_singlets, 0.8 * pairs, 3e-5);
    close(step.radiative_events[0], fluorescence, 3e-5);
    close(step.radiative_decay_energy_j[0], fluorescence * 4e-19, 3e-5);
    close(step.annihilation_heat_j, pairs * 2.8e-19, 3e-5);
    // Singlets created by a pair must first survive or decay on their own
    // clock; the pair event itself does not emit the delayed photon.
    close(
        step.annihilation_singlets,
        step.radiative_events[0] + step.populations[1],
        2e-12,
    );
    assert!(step.radiative_events[0] < step.annihilation_singlets);
    assert_eq!(step.pump_events, [0.0; 2]);
    assert_eq!(step.radiative_events[1], 0.0);
    assert_eq!(step.nonradiative_events, [0.0; 2]);
    assert_eq!(step.nonradiative_heat_j, [0.0; 2]);
    assert_eq!(step.crossing_events, [0.0; 2]);
    assert_eq!(step.reverse_crossing_heat_absorbed_j, 0.0);
    assert!(step.residual_j.abs() < 1e-18);
    assert!(result.substeps <= 1024);
    eprintln!(
        "finite dark TTA: singlets={:.9e}, endpoint_photons_s={:.9e}, substeps={}, kinetic_evaluations={}, residual_j={:.9e}",
        step.populations[1],
        step.populations[1] * singlet_rate,
        result.substeps,
        result.kinetic_evaluations,
        step.residual_j,
    );
}

#[test]
fn g1_triplet_pumped_nonlinear_equilibrium_preserves_populations_and_pair_accounts() {
    let model =
        SingletTripletDecay::try_new([4e-19, 3e-19], [1e9, 0.0], [0.0; 2], [0.0; 2]).unwrap();
    let plan = SaturableSingletTripletPlan::prepare(
        model,
        [0.0, 2.0],
        [0.0, 4e-19],
        1.0,
        DecayLimits::default(),
        || false,
    )
    .unwrap();
    let law = TripletAnnihilation::try_new(1e-18, 1e-6, 0.8).unwrap();

    // Pumping G->T at W balances two T removed per pair. With no other T
    // channel, G=2*a*T^2/W and S=phi*a*T^2/ks. Writing
    // c=a*(2/W+phi/ks), inventory conservation gives
    // T=2*N/(1+sqrt(1+4*c*N)). Constants below use N=1e12, a=1e-12,
    // phi=.8, W=2/s, ks=1e9/s, independently evaluated at 85 digits.
    let initial = [
        381_966_011_081.188_66,
        305.572_808_864_950_9,
        618_033_988_613.238_5,
    ];
    let pairs = 381_966_011_081.188_66;
    let result = plan
        .advance_with_annihilation(initial, Some(law), controls(1e-10, 128), || false)
        .unwrap();
    let step = &result.kinetics;
    for (actual, expected) in step.populations.into_iter().zip(initial) {
        close(actual, expected, 5e-12);
    }
    close(step.populations.iter().sum(), 1e12, 2e-12);
    close(step.annihilation_events, pairs, 5e-12);
    close(step.annihilation_singlets, 0.8 * pairs, 5e-12);
    close(step.annihilation_heat_j, 2.8e-19 * pairs, 5e-12);
    close(step.pump_events[1], 2.0 * pairs, 5e-12);
    close(step.absorbed_pump_energy_j[1], 8e-19 * pairs, 5e-12);
    close(step.pump_relaxation_heat_j[1], 2e-19 * pairs, 5e-12);
    close(step.radiative_events[0], 0.8 * pairs, 5e-12);
    close(step.radiative_decay_energy_j[0], 3.2e-19 * pairs, 5e-12);
    assert_eq!(step.pump_events[0], 0.0);
    assert_eq!(step.radiative_events[1], 0.0);
    assert_eq!(step.nonradiative_events, [0.0; 2]);
    assert_eq!(step.nonradiative_heat_j, [0.0; 2]);
    assert_eq!(step.crossing_events, [0.0; 2]);
    assert_eq!(step.reverse_crossing_heat_absorbed_j, 0.0);
    assert!(step.residual_j.abs() < 1e-18);
    assert!(result.substeps <= 128);
    eprintln!(
        "finite pumped TTA equilibrium: substeps={}, kinetic_evaluations={}, pair_events={:.9e}, residual_j={:.9e}",
        result.substeps, result.kinetic_evaluations, step.annihilation_events, step.residual_j,
    );
}
