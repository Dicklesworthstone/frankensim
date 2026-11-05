//! G0/G3/G4 prescribed-geometry interpolation, domain and cancellation checks.

use super::*;

fn knot(time: f64, volume: f64) -> FusionVolumeKnot {
    FusionVolumeKnot {
        time: Time::new(time),
        volume: Volume::new(volume),
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn geometric_interpolation_retains_authored_knots() {
    let knots = vec![knot(0.0, 64.0), knot(2.0, 4.0)];
    let history = FusionVolumeHistory::new(knots.clone(), || false).unwrap();
    assert_eq!(history.knots(), knots);
    assert_eq!(history.duration(), Time::new(2.0));
    assert_eq!(history.initial_volume(), Volume::new(64.0));
    assert_eq!(history.final_volume(), Volume::new(4.0));
    for (time, expected) in [(0.5, 32.0), (1.0, 16.0), (1.5, 8.0)] {
        close(
            history.volume_at(Time::new(time)).unwrap().value(),
            expected,
            5e-15,
        );
    }
    let shared = history.clone();
    assert_eq!(shared, history);
    assert_eq!(shared.knots().as_ptr(), history.knots().as_ptr());
}

#[test]
fn nonuniform_holds_turnarounds_and_strict_next_knot() {
    let history = FusionVolumeHistory::new(
        vec![
            knot(0.0, 0.025),
            knot(0.04, 0.0125),
            knot(0.06, 0.0125),
            knot(0.1, 0.02),
        ],
        || false,
    )
    .unwrap();
    for knot in history.knots() {
        assert_eq!(
            history.volume_at(knot.time).unwrap().value().to_bits(),
            knot.volume.value().to_bits()
        );
    }
    assert_eq!(
        history
            .volume_at(Time::new(0.05))
            .unwrap()
            .value()
            .to_bits(),
        0.0125_f64.to_bits()
    );
    close(
        history.volume_at(Time::new(0.02)).unwrap().value(),
        (0.025_f64 * 0.0125).sqrt(),
        5e-15,
    );
    close(
        history.volume_at(Time::new(0.08)).unwrap().value(),
        (0.0125_f64 * 0.02).sqrt(),
        5e-15,
    );
    for (time, next) in [
        (0.0, Some(0.04)),
        (0.03, Some(0.04)),
        (0.04, Some(0.06)),
        (0.06, Some(0.1)),
        (0.1, None),
    ] {
        assert_eq!(
            history.next_knot_after(Time::new(time)).unwrap(),
            next.map(Time::new)
        );
    }
    let instant = FusionVolumeHistory::new(vec![knot(0.0, 1.25)], || false).unwrap();
    assert_eq!(instant.duration(), Time::new(0.0));
    assert_eq!(
        instant.volume_at(Time::new(0.0)).unwrap(),
        Volume::new(1.25)
    );
    assert_eq!(instant.next_knot_after(Time::new(0.0)).unwrap(), None);
    assert_eq!(
        instant.volume_at(Time::new(f64::from_bits(1))),
        Err(DtFusionError::OutsideDomain)
    );
}

#[test]
fn invalid_domains_bounded_work_and_cancellation_refuse() {
    for knots in [
        vec![],
        vec![knot(1.0, 1.0)],
        vec![knot(-1.0, 1.0)],
        vec![knot(0.0, 0.0)],
        vec![knot(0.0, -1.0)],
        vec![knot(0.0, f64::INFINITY)],
        vec![knot(0.0, f64::NAN)],
        vec![knot(0.0, 1.0), knot(0.0, 2.0)],
        vec![knot(0.0, 1.0), knot(2.0, 1.0), knot(1.0, 2.0)],
        vec![knot(0.0, 1.0), knot(f64::INFINITY, 2.0)],
        vec![knot(0.0, 1.0), knot(f64::NAN, 2.0)],
    ] {
        assert!(matches!(
            FusionVolumeHistory::new(knots, || false),
            Err(DtFusionError::Invalid(_))
        ));
    }
    let maximum: Vec<_> = (0..FusionVolumeHistory::MAX_KNOTS)
        .map(|i| knot(i as f64, 1.0))
        .collect();
    let history = FusionVolumeHistory::new(maximum.clone(), || false).unwrap();
    let mut over_budget = maximum.clone();
    over_budget.push(knot(FusionVolumeHistory::MAX_KNOTS as f64, 1.0));
    assert!(matches!(
        FusionVolumeHistory::new(over_budget, || false),
        Err(DtFusionError::Invalid(_))
    ));
    for time in [-1.0, 256.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
        assert!(history.volume_at(Time::new(time)).is_err());
        assert!(history.next_knot_after(Time::new(time)).is_err());
    }
    assert_eq!(
        FusionVolumeHistory::new(maximum.clone(), || true),
        Err(DtFusionError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        FusionVolumeHistory::new(maximum, || {
            polls += 1;
            polls >= 4
        }),
        Err(DtFusionError::Cancelled)
    );
}

#[test]
fn tiny_changes_and_extreme_ratios_remain_resolved() {
    // Taking ln(V1) - ln(V0) at this scale would lose the small variation.
    let initial = 1e200_f64;
    let final_volume = f64::from_bits(initial.to_bits() + 8);
    let tiny =
        FusionVolumeHistory::new(vec![knot(0.0, initial), knot(1.0, final_volume)], || false)
            .unwrap();
    let midpoint = tiny.volume_at(Time::new(0.5)).unwrap().value();
    assert!(midpoint > initial && midpoint < final_volume);
    assert!((midpoint.to_bits() as i128 - (initial.to_bits() + 4) as i128).abs() <= 1);

    // The endpoint ratio overflows; at half time even its square root exceeds
    // f64::MAX. The true interpolated volume is nevertheless finite and positive.
    let low = f64::from_bits(1);
    let high = f64::MAX;
    let expected = (0.5 * low.ln() + 0.5 * high.ln()).exp();
    for (left, right) in [(low, high), (high, low)] {
        let history =
            FusionVolumeHistory::new(vec![knot(0.0, left), knot(1.0, right)], || false).unwrap();
        close(
            history.volume_at(Time::new(0.5)).unwrap().value(),
            expected,
            3e-13,
        );
        for time in [0.01, 0.25, 0.49, 0.51, 0.75, 0.99] {
            let volume = history.volume_at(Time::new(time)).unwrap().value();
            assert!(volume >= low && volume <= high);
        }
    }
}
