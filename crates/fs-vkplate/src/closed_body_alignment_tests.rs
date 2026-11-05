//! Exact streaming superposition and geometry oracles; no BEM/table-fit claim.
use super::{CoherentObserverMix, same_translated_plane};

fn pulse(t: usize) -> f64 {
    match t {
        0 => 1.0,
        5 => -0.5,
        19 => 0.25,
        _ => 0.0,
    }
}

fn delayed(t: usize, lag: usize, signal: impl FnOnce(usize) -> f64) -> f64 {
    t.checked_sub(lag).map_or(0.0, signal)
}

#[test]
fn g1_common_latency_preserves_physical_delays_and_counts_the_baseline_once() {
    let latencies = [3, 7, 1];
    let paths = [4, 9, 2];
    let gains = [2.0, -0.5, 1.25];
    let mut mixer = CoherentObserverMix::new(&latencies);
    let mut calls = [0; 3];
    for t in 0..80 {
        let actual = mixer.step(pulse(t), |i| {
            calls[i] += 1;
            gains[i] * delayed(t, latencies[i] + paths[i], pulse)
        });
        // Analytically, the one baseline and every physical response share
        // max(D_i), while each ray's physical path delay remains untouched.
        let expected = delayed(t, 7, pulse)
            + (0..3)
                .map(|i| gains[i] * delayed(t, 7 + paths[i], pulse))
                .sum::<f64>();
        assert_eq!(actual, expected, "physical observation sample {t}");
    }
    assert_eq!(
        calls, [80; 3],
        "every observer advances exactly once per sample"
    );
}

#[test]
fn g1_unequal_filter_latencies_do_not_spoil_destructive_interference() {
    let latencies = [2, 11, 5];
    let mut mixer = CoherentObserverMix::new(&latencies);
    let epsilon = 2.0_f64.powi(-60);
    let gains = [1.0, epsilon, -1.0];
    let signal = |t: usize| (0.37 * t as f64).sin();
    for t in 0..160 {
        let actual = mixer.step(0.0, |i| gains[i] * delayed(t, latencies[i], signal));
        // The weak source is deliberately summed BETWEEN opposite strong
        // fields. Ordinary f64 accumulation erases it before cancellation.
        let expected = epsilon * delayed(t, 11, signal);
        assert_eq!(actual, expected, "false phase offset at sample {t}");
    }
}

#[test]
fn g3_cloned_alignment_preserves_pending_wavefronts_through_ring_wraps() {
    let latencies = [1, 17];
    let mut original = CoherentObserverMix::new(&latencies);
    let step = |mixer: &mut CoherentObserverMix, t: usize| {
        mixer.step(pulse(t), |i| {
            (i + 1) as f64 * delayed(t, latencies[i] + 3, pulse)
        })
    };
    for t in 0..13 {
        step(&mut original, t);
    }
    let mut copied = original.clone();
    for t in 13..96 {
        assert_eq!(
            step(&mut original, t).to_bits(),
            step(&mut copied, t).to_bits()
        );
    }
}

#[test]
fn g0_translated_floor_matches_the_plane_not_an_arbitrary_tangent_point() {
    let baked = ([0.0, 0.0, -0.25], [0.0, 0.0, 1.0]);
    let origin = [0.5, -0.75, 0.125];
    assert!(same_translated_plane(
        baked,
        ([0.0, 0.0, -0.125], baked.1),
        origin
    ));
    assert!(same_translated_plane(
        baked,
        ([8.0, -3.0, -0.125], baked.1),
        origin
    ));
    assert!(!same_translated_plane(
        baked,
        ([0.0, 0.0, -0.25], baked.1),
        origin
    ));
    assert!(!same_translated_plane(
        baked,
        ([0.0, 0.0, -0.125], [0.0, 0.0, -1.0]),
        origin
    ));
    assert!(!same_translated_plane(
        baked,
        ([f64::NAN, 0.0, -0.125], baked.1),
        origin
    ));
    let side_floor = ([0.0, -0.25, 0.0], [0.0, 1.0, 0.0]);
    assert!(same_translated_plane(
        side_floor,
        ([1.0, -1.0, 2.0], side_floor.1),
        origin
    ));
}
