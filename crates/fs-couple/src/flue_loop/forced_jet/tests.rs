//! G0/G1/G3: the original response on an accepted convective clock.
//! These are numerical/source-memory checks, not measured instrument qualification.
use super::*;

fn kernel(rate: u32) -> ForcedJetKernel {
    ForcedJetKernel::new(
        ForcedJetSpec { spacing: 0.5, window: 8.0 },
        JetExit::UNIFORM, 0.001, 0.00225, 1.5e-5, 12.0, 1.0, 4096,
        1.0 / f64::from(rate),
    ).unwrap()
}

fn input(sample: usize) -> f64 {
    ((sample * 17 % 43) as f64 - 21.0) * 0.001
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 1e-10 * expected.abs().max(1e-6),
        "convective response {actual:e} differs from integral oracle {expected:e}");
}

// Independent integration by parts through the response owner's accumulated
// integrals. Production instead integrates local derivative knots, which avoids
// cancellation in old tails. Every interval carries its own physical dt/dtau.
fn integral_oracle(kernel: &ForcedJetKernel, intervals: &[JetInterval]) -> f64 {
    let (mut age, mut result) = (0.0, 0.0);
    for interval in intervals.iter().rev() {
        if age >= kernel.step.window() { break; }
        let end = age + interval.advance;
        let mean = (kernel.step.integral_at(end) - kernel.step.integral_at(age))
            / interval.advance;
        result += kernel.dt / interval.advance * (
            (mean - kernel.step.value_at(age)) * interval.end_velocity_m_s
            + (kernel.step.value_at(end) - mean) * interval.start_velocity_m_s
        );
        age = end;
    }
    result
}

#[test]
fn g0_uniform_speed_keeps_original_hat_arithmetic_at_three_clocks() {
    for rate in [48_000, 96_000, 192_000] {
        let mut actual = kernel(rate);
        let weights = actual.step.hat_weights(12.0, actual.half_width_m, actual.dt).unwrap();
        let mut ring = vec![0.0; actual.ring.len()];
        let mut head = 0;
        for sample in 0..(ring.len() + 17) {
            let mut known = 0.0;
            for (k, g) in weights.iter().enumerate().skip(1) {
                known += g * ring[(head + ring.len() - k) % ring.len()];
            }
            let prepared = actual.prepare(12.0).unwrap().unwrap();
            assert_eq!(prepared.0.to_bits(), known.to_bits());
            assert_eq!(prepared.1.to_bits(), weights[0].to_bits());
            let velocity = input(sample);
            actual.push(velocity);
            ring[head] = velocity;
            head = (head + 1) % ring.len();
        }
    }
}

#[test]
fn g1_convective_intervals_reduce_to_the_uniform_response() {
    let mut uniform = kernel(48_000);
    let mut interval = uniform.clone();
    interval.nonuniform_history = true;
    for sample in 0..(uniform.ring.len() + 17) {
        let a = uniform.prepare(12.0).unwrap().unwrap();
        let b = interval.prepare(12.0).unwrap().unwrap();
        close(b.0, a.0);
        close(b.1, a.1);
        let velocity = input(sample);
        uniform.push(velocity);
        interval.push(velocity);
    }
    assert!(!interval.nonuniform_history, "a completely uniform history reuses cached hats");
}

#[test]
fn g1_changing_speed_integrates_accepted_intervals_without_retiming_them() {
    let mut actual = kernel(48_000);
    let mut intervals = Vec::new();
    let mut velocities = Vec::new();
    let mut previous = 0.0;
    let (mut largest, mut old_model_difference): (f64, f64) = (0.0, 0.0);
    // Cross the allocated ring boundary as well as the response's held tail.
    for sample in 0..(actual.ring.len() + 37) {
        let bulk = [12.0, 3.0, 24.0, 6.0][(sample / 17) % 4];
        let velocity = input(sample);
        let advance = actual.centre_over_bulk * bulk * actual.dt / actual.half_width_m;
        intervals.push(JetInterval {
            start_velocity_m_s: previous, end_velocity_m_s: velocity, advance,
        });
        let (known, current) = actual.prepare(bulk).unwrap().unwrap();
        let displacement = known + current * velocity;
        close(displacement, integral_oracle(&actual, &intervals));

        // This is the previous implementation: all past samples receive THIS
        // sample's speed. A nonzero difference is the falsifier for that shortcut.
        velocities.push(velocity);
        let weights = actual.step.hat_weights(
            actual.centre_over_bulk * bulk, actual.half_width_m, actual.dt).unwrap();
        let retimed: f64 = weights.iter().zip(velocities.iter().rev())
            .map(|(g, v)| g * v).sum();
        largest = largest.max(displacement.abs());
        old_model_difference = old_model_difference.max((retimed - displacement).abs());
        actual.push(velocity);
        previous = velocity;
    }
    assert!(largest > 0.0);
    assert!(old_model_difference > 0.01 * largest,
        "fixture must distinguish accepted convection from instantaneous retiming");
}

#[test]
fn g0_trials_refusals_and_rest_do_not_advance_accepted_convection() {
    let mut original = kernel(48_000);
    for sample in 0..23 {
        original.prepare(12.0).unwrap();
        original.push(input(sample));
    }
    let mut trial = original.clone();
    trial.prepare(24.0).unwrap();
    trial.prepare(0.0).unwrap();
    // Below the declared minimum speed the retained intervals still reach the rest state
    // (ruling 8132), so the trial is admitted.
    assert!(trial.prepare(0.5).unwrap().is_some());
    for invalid in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(trial.prepare(invalid).is_err());
    }
    assert_eq!(trial.ring, original.ring);
    assert_eq!(trial.head, original.head);
    assert_eq!(trial.used, original.used);
    assert_eq!(trial.prepare(12.0).unwrap(), original.prepare(12.0).unwrap());

    let mut one_rest = original.clone();
    let mut many_rest = original.clone();
    one_rest.prepare(0.0).unwrap();
    one_rest.push(17.0); // A stopped jet admits no transverse forcing.
    for _ in 0..(2 * many_rest.ring.len()) {
        assert_eq!(many_rest.prepare(0.0).unwrap(), Some((0.0, 0.0)));
        many_rest.push(-19.0);
    }
    assert_eq!(many_rest.ring, original.ring);
    assert_eq!(many_rest.head, original.head);
    assert_eq!(many_rest.used, original.used);
    assert_eq!(many_rest.last_velocity_m_s, 0.0);
    assert_eq!(many_rest.prepare(6.0).unwrap(), one_rest.prepare(6.0).unwrap());
}

#[test]
fn g1_a_jet_below_the_minimum_reads_rest_until_its_ring_ends_then_freezes() {
    // Ruling 8132. Held at half the declared 1 m/s minimum from rest, the response needs more
    // intervals than the ring keeps. Until the ring fills, older input is the rest state, so a
    // kernel whose ring covers the whole response (minimum 0.25 m/s) agrees bit for bit.
    let mut short = kernel(48_000);
    let mut long = ForcedJetKernel::new(
        ForcedJetSpec { spacing: 0.5, window: 8.0 },
        JetExit::UNIFORM, 0.001, 0.00225, 1.5e-5, 12.0, 0.25, 4096, 1.0 / 48_000.0,
    ).unwrap();
    let n = short.ring.len();
    assert!(long.ring.len() > 2 * n);
    for sample in 0..n {
        let a = short.prepare(0.5).unwrap().expect("the ring still reaches the rest state");
        let b = long.prepare(0.5).unwrap().unwrap();
        assert_eq!((a.0.to_bits(), a.1.to_bits()), (b.0.to_bits(), b.1.to_bits()));
        short.push(input(sample));
        long.push(input(sample));
    }
    // The ring is full and its intervals fall short of the window: nothing is prepared and the
    // clock freezes as at rest, so no later push enters the history.
    let (ring, head) = (short.ring.clone(), short.head);
    for sample in n..(n + 40) {
        assert_eq!(short.prepare(0.5).unwrap(), None);
        assert_eq!(short.current(), (0.0, 0.0));
        short.push(input(sample));
    }
    assert_eq!((short.ring == ring, short.head, short.used), (true, head, n));
}

#[test]
fn measured_exit_reynolds_neighbourhood_is_closed_at_both_ends() {
    use crate::flue_loop::forced_jet::{measured_exit_reynolds_admitted as admitted,
        MEASURED_EXIT_RE_H, MEASURED_EXIT_RE_H_FACTOR};
    assert_eq!(MEASURED_EXIT_RE_H, 200.0);
    assert_eq!(MEASURED_EXIT_RE_H_FACTOR, 2.0);
    for ok in [100.0, 100.000001, 200.0, 399.999999, 400.0] { assert!(admitted(ok), "{ok}"); }
    // Planted negatives: just outside either end, zero, negative, non-finite, and the cards' own values.
    for no in [99.999999, 400.000001, 0.0, -200.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY,
        1224.0, 1755.0] { assert!(!admitted(no), "{no}"); }
}

// The jet-length track (MB 23643) on the plumbing row: the card's station is 3 b, the span 2.8 to
// 3.47 b. A track held at the card's length changes no bit; a parallel span read at any length is
// the kernel built at that length, bit for bit; a length step reads the same retained intervals at
// the new station; a spread span is exact at its far end; refusals change nothing.
const CARD_LENGTH_M: f64 = 0.00225;
const NEAR_LENGTH_M: f64 = 0.0021;
const FAR_LENGTH_M: f64 = 0.0026;

fn kernel_at(cut_up_m: f64) -> ForcedJetKernel {
    ForcedJetKernel::new(ForcedJetSpec { spacing: 0.5, window: 8.0 }, JetExit::UNIFORM, 0.001,
        cut_up_m, 1.5e-5, 12.0, 1.0, 4096, 1.0 / 48_000.0).unwrap()
}

fn spread_at(cut_up_m: f64) -> ForcedJetKernel {
    ForcedJetKernel::new_with_spreading(ForcedJetSpec { spacing: 0.5, window: 8.0 }, JetExitLaw::Kutta,
        JetSpreadingLaw::Bickley, JetExit::UNIFORM, 0.001, cut_up_m, 1.5e-5, 12.0, 1.0, 4096,
        1.0 / 48_000.0).unwrap()
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

// Samples `from..to`: one speed (the cached-weight path) before 250, a speed that changes every
// seven samples after it, so the intervals carry their own convective ages.
fn drive(kernel: &mut ForcedJetKernel, from: usize, to: usize) -> Vec<(u64, u64)> {
    (from..to).map(|n| {
        let bulk = if n < 250 { 12.0 } else { 12.0 + 0.5 * ((n / 7) % 3) as f64 };
        let (known, gain) = kernel.prepare(bulk).unwrap().unwrap();
        kernel.push(input(n));
        (known.to_bits(), gain.to_bits())
    }).collect()
}

#[test]
fn g0_a_track_held_at_the_card_length_changes_no_bit() {
    let mut plain = kernel_at(CARD_LENGTH_M);
    let mut held = kernel_at(CARD_LENGTH_M);
    held.span_jet_lengths(CARD_LENGTH_M, CARD_LENGTH_M).unwrap();
    assert_eq!(held.jet_length_span_m(), Some((CARD_LENGTH_M, CARD_LENGTH_M)));
    assert_eq!(held.jet_length_m().to_bits(), CARD_LENGTH_M.to_bits());
    assert_eq!(bits(held.step.values()), bits(plain.step.values()));
    held.set_jet_length(CARD_LENGTH_M).unwrap();
    assert_eq!(drive(&mut held, 0, 400), drive(&mut plain, 0, 400));
}

#[test]
fn g1_a_parallel_span_read_at_any_length_is_the_kernel_built_there_bit_for_bit() {
    for length in [NEAR_LENGTH_M, 0.0023, FAR_LENGTH_M] {
        let mut spanned = kernel_at(CARD_LENGTH_M);
        spanned.span_jet_lengths(NEAR_LENGTH_M, FAR_LENGTH_M).unwrap();
        assert_eq!(spanned.jet_length_m().to_bits(), CARD_LENGTH_M.to_bits(), "read at the card's length inside");
        spanned.set_jet_length(length).unwrap();
        let mut built = kernel_at(length);
        assert_eq!(bits(spanned.step.values()), bits(built.step.values()), "{length}");
        assert_eq!(spanned.step.station().to_bits(), built.step.station().to_bits());
        assert_eq!(spanned.ring.len(), built.ring.len());
        assert_eq!(drive(&mut spanned, 0, 400), drive(&mut built, 0, 400), "{length}");
    }
}

#[test]
fn g1_a_length_step_reads_the_same_retained_history_at_the_new_station() {
    let mut moved = kernel_at(CARD_LENGTH_M);
    moved.span_jet_lengths(NEAR_LENGTH_M, FAR_LENGTH_M).unwrap();
    moved.set_jet_length(NEAR_LENGTH_M).unwrap();
    let mut short = kernel_at(NEAR_LENGTH_M);
    let mut long = kernel_at(FAR_LENGTH_M);
    // Before the step the moved jet is the short one; the step publishes nothing by itself.
    assert_eq!(drive(&mut moved, 0, 200), drive(&mut short, 0, 200));
    let _ = drive(&mut long, 0, 200);
    let used = moved.used;
    moved.set_jet_length(FAR_LENGTH_M).unwrap();
    assert_eq!(moved.used, used);
    // After it, the same history reads at the long station: on the cached speed, then a varied one.
    assert_eq!(drive(&mut moved, 200, 500), drive(&mut long, 200, 500));
}

#[test]
fn g1_a_spread_span_is_the_kernel_built_at_its_far_end_bit_for_bit() {
    let mut spanned = spread_at(CARD_LENGTH_M);
    spanned.span_jet_lengths(NEAR_LENGTH_M, FAR_LENGTH_M).unwrap();
    spanned.set_jet_length(FAR_LENGTH_M).unwrap();
    let mut built = spread_at(FAR_LENGTH_M);
    assert_eq!(bits(spanned.step.values()), bits(built.step.values()));
    assert_eq!(spanned.station_spread(), built.station_spread());
    assert_eq!(drive(&mut spanned, 0, 400), drive(&mut built, 0, 400));
    // Short of it the base spreads out to the far end: a stated, small difference.
    spanned.set_jet_length(NEAR_LENGTH_M).unwrap();
    let near = spread_at(NEAR_LENGTH_M);
    assert_ne!(bits(spanned.step.values()), bits(near.step.values()));
    let rel = (spanned.step.held_gain() - near.step.held_gain()).abs() / near.step.held_gain().abs();
    assert!(rel < 1e-2, "{rel:e}");
}

#[test]
fn g0_lengths_off_a_prepared_span_are_refused_and_change_nothing() {
    let mut cold = kernel_at(CARD_LENGTH_M);
    assert!(cold.set_jet_length(NEAR_LENGTH_M).is_err(), "no span");
    assert!(cold.set_jet_length(CARD_LENGTH_M).is_ok(), "its own length is no change");
    for (near, far) in [(FAR_LENGTH_M, NEAR_LENGTH_M), (0.0, FAR_LENGTH_M), (-NEAR_LENGTH_M, FAR_LENGTH_M),
        (f64::NAN, FAR_LENGTH_M), (NEAR_LENGTH_M, f64::INFINITY)] {
        assert!(cold.span_jet_lengths(near, far).is_err(), "({near}, {far})");
        assert_eq!(cold.jet_length_span_m(), None);
    }
    let mut warm = kernel_at(CARD_LENGTH_M);
    let _ = drive(&mut warm, 0, 1);
    assert!(warm.span_jet_lengths(NEAR_LENGTH_M, FAR_LENGTH_M).is_err(), "a jet that has accepted an interval");
    cold.span_jet_lengths(NEAR_LENGTH_M, FAR_LENGTH_M).unwrap();
    assert!(cold.span_jet_lengths(NEAR_LENGTH_M, FAR_LENGTH_M).is_err(), "one preparation");
    let values = bits(cold.step.values());
    for length in [NEAR_LENGTH_M * 0.999, FAR_LENGTH_M * 1.001, f64::NAN, -CARD_LENGTH_M] {
        assert!(cold.set_jet_length(length).is_err(), "{length}");
        assert_eq!(cold.jet_length_m().to_bits(), CARD_LENGTH_M.to_bits());
        assert_eq!(bits(cold.step.values()), values);
    }
}
