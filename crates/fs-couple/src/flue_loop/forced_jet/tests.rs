//! G0/G1/G3: the original response on an accepted convective clock.
//! These are numerical/source-memory checks, not measured instrument qualification.
use super::*;

fn kernel(rate: u32) -> ForcedJetKernel {
    ForcedJetKernel::new(
        ForcedJetSpec { spacing: 0.5, window: 8.0 },
        false, 0.001, 0.00225, 1.5e-5, 12.0, 1.0, 4096,
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
            let prepared = actual.prepare(12.0).unwrap();
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
        let a = uniform.prepare(12.0).unwrap();
        let b = interval.prepare(12.0).unwrap();
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
        let (known, current) = actual.prepare(bulk).unwrap();
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
    for invalid in [f64::NAN, f64::INFINITY, -1.0, 0.5] {
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
        assert_eq!(many_rest.prepare(0.0).unwrap(), (0.0, 0.0));
        many_rest.push(-19.0);
    }
    assert_eq!(many_rest.ring, original.ring);
    assert_eq!(many_rest.head, original.head);
    assert_eq!(many_rest.used, original.used);
    assert_eq!(many_rest.last_velocity_m_s, 0.0);
    assert_eq!(many_rest.prepare(6.0).unwrap(), one_rest.prepare(6.0).unwrap());
}
