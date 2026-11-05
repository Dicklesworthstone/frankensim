use super::*;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-15 * expected.abs().max(1e-300),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_half_life_and_mean_lifetime_have_distinct_physical_meanings() {
    let law = DecayLawV1::from_half_life(Time::new(10.0)).unwrap();
    assert_eq!(law.half_life(), Some(Time::new(10.0)));
    close(law.decay_constant().value(), LN_2 / 10.0);
    close(law.mean_lifetime().unwrap().unwrap().value(), 10.0 / LN_2);
    let half = law.fractions(Time::new(10.0)).unwrap();
    assert_eq!(
        half,
        DecayFractionsV1 {
            remaining: 0.5,
            decayed: 0.5
        }
    );
    close(law.fractions(Time::new(30.0)).unwrap().remaining, 0.125);
    close(
        law.fractions(law.mean_lifetime().unwrap().unwrap())
            .unwrap()
            .remaining,
        (-1.0_f64).exp(),
    );
}

#[test]
fn g3_partition_and_time_rescaling_preserve_the_same_decay() {
    let law = DecayLawV1::from_half_life(Time::new(7.0)).unwrap();
    let a = law.fractions(Time::new(3.0)).unwrap();
    let b = law.fractions(Time::new(17.0)).unwrap();
    let total = law.fractions(Time::new(20.0)).unwrap();
    close(a.remaining * b.remaining, total.remaining);
    close(a.decayed + a.remaining * b.decayed, total.decayed);
    for scale in [1e-200, 1e-50, 1.0, 1e50, 1e200] {
        let scaled = DecayLawV1::from_half_life(Time::new(7.0 * scale)).unwrap();
        let result = scaled.fractions(Time::new(20.0 * scale)).unwrap();
        close(result.remaining, total.remaining);
        close(result.decayed, total.decayed);
    }
}

#[test]
fn g0_losses_and_survivors_retain_separate_small_tails() {
    let law = DecayLawV1::from_decay_constant(Frequency::new(1.0)).unwrap();
    let tiny = law.fractions(Time::new(1e-25)).unwrap();
    assert_eq!(tiny.remaining, 1.0);
    close(tiny.decayed, 1e-25);
    let late = law.fractions(Time::new(100.0)).unwrap();
    assert_eq!(late.decayed, 1.0);
    close(late.remaining, (-100.0_f64).exp());
    let saturated = law.fractions(Time::new(f64::MAX)).unwrap();
    assert_eq!(
        saturated,
        DecayFractionsV1 {
            remaining: 0.0,
            decayed: 1.0
        }
    );
    let rapid = DecayLawV1::from_decay_constant(Frequency::new(f64::MAX)).unwrap();
    assert_eq!(rapid.fractions(Time::new(2.0)).unwrap(), saturated);
}

#[test]
fn g0_stability_and_invalid_or_unrepresentable_inputs_are_explicit() {
    let stable = DecayLawV1::stable();
    assert_eq!(stable.half_life(), None);
    assert_eq!(stable.mean_lifetime().unwrap(), None);
    assert_eq!(stable.decay_constant().value(), 0.0);
    assert_eq!(
        stable,
        DecayLawV1::from_decay_constant(Frequency::new(-0.0)).unwrap()
    );
    for t in [0.0, 1.0, f64::MAX] {
        assert_eq!(
            stable.fractions(Time::new(t)).unwrap(),
            DecayFractionsV1 {
                remaining: 1.0,
                decayed: 0.0
            }
        );
    }
    for bad in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(DecayLawV1::from_half_life(Time::new(bad)).is_err());
        assert!(DecayLawV1::from_decay_constant(Frequency::new(bad)).is_err());
        assert!(stable.fractions(Time::new(bad)).is_err());
    }
    assert!(DecayLawV1::from_half_life(Time::new(0.0)).is_err());
    assert!(DecayLawV1::from_half_life(Time::new(f64::from_bits(1))).is_err());
    assert!(DecayLawV1::from_decay_constant(Frequency::new(f64::from_bits(1))).is_err());
    let slow = DecayLawV1::from_half_life(Time::new(f64::MAX)).unwrap();
    assert_eq!(slow.half_life(), Some(Time::new(f64::MAX)));
    assert!(slow.mean_lifetime().is_err());
    assert!(slow.fractions(Time::new(f64::MIN_POSITIVE)).is_err());
    close(slow.fractions(Time::new(f64::MAX)).unwrap().remaining, 0.5);
}

#[test]
fn g3_independent_removal_rates_add_and_survival_probabilities_multiply() {
    let a = DecayLawV1::from_half_life(Time::new(2.0)).unwrap();
    let b = DecayLawV1::from_half_life(Time::new(3.0)).unwrap();
    let combined = a.competing_with(b).unwrap();
    close(combined.half_life().unwrap().value(), 1.2);
    close(
        combined.fractions(Time::new(5.0)).unwrap().remaining,
        a.fractions(Time::new(5.0)).unwrap().remaining
            * b.fractions(Time::new(5.0)).unwrap().remaining,
    );
    assert_eq!(a.competing_with(DecayLawV1::stable()).unwrap(), a);
    assert_eq!(DecayLawV1::stable().competing_with(a).unwrap(), a);
    let rapid = DecayLawV1::from_decay_constant(Frequency::new(f64::MAX)).unwrap();
    assert!(rapid.competing_with(rapid).is_err());
}
