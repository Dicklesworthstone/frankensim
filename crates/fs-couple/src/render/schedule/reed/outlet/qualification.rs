//! G0/G3/G5 observation-policy regressions; physical source parameters stay fixed.
use super::*;

const SOURCE: &str = include_str!("../../../../../examples/reed-baffled-mixture.performance");

#[test]
fn legacy_listener_refused_in_acceptance() {
    let legacy = SOURCE.replace("baffled_outlet 0 0 0.5 4 16 1000", "listener_m 0.5");
    let error = ReedPerformance::from_bytes_with_policy(
        legacy.as_bytes(),
        37,
        ObservationPolicy::Acceptance,
    )
    .err()
    .unwrap();
    assert!(matches!(
        error,
        ReedPerformanceError::Observation(ObservationError::LegacyListenerInAcceptance)
    ));
    let diagnostic = ReedPerformance::from_bytes_with_policy(
        legacy.as_bytes(),
        37,
        ObservationPolicy::Diagnostic,
    )
    .unwrap();
    assert!(diagnostic.outlet_receiver().is_none());
    assert_eq!(
        diagnostic.observation_policy(),
        ObservationPolicy::Diagnostic
    );
}

#[test]
fn baffled_outlet_requires_flanged() {
    for termination in ["closed", "ideal-open", "unflanged"] {
        let text = SOURCE.replace("termination flanged", &format!("termination {termination}"));
        for policy in [ObservationPolicy::Acceptance, ObservationPolicy::Diagnostic] {
            let error = ReedPerformance::from_bytes_with_policy(text.as_bytes(), 37, policy)
                .err()
                .unwrap();
            assert!(
                error
                    .to_string()
                    .contains("same flanged terminal feedback load")
            );
        }
    }
    let source = ReedPerformance::from_bytes_with_policy(
        SOURCE.as_bytes(),
        37,
        ObservationPolicy::Acceptance,
    )
    .unwrap();
    assert!(source.outlet_receiver().is_some());
}

#[test]
fn observer_bitwise_repeat() {
    let build = || {
        ReedPerformance::from_bytes_with_policy(
            SOURCE.as_bytes(),
            113,
            ObservationPolicy::Acceptance,
        )
        .unwrap()
    };
    let (mut a, mut b) = (build(), build());
    let mut x = [0.0; 701];
    let mut y = [0.0; 701];
    for part in x.chunks_mut(113) {
        a.pressure_block(part).unwrap();
    }
    for part in y.chunks_mut(37) {
        b.pressure_block(part).unwrap();
    }
    assert!(x.iter().any(|p| *p != 0.0));
    assert!(x.iter().zip(y).all(|(p, q)| p.to_bits() == q.to_bits()));
    assert_eq!(a.observation_policy(), ObservationPolicy::Acceptance);
}
