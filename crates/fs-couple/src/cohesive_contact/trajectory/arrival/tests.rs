use super::*;
use std::convert::Infallible;

fn controls() -> ArrivalControls {
    ArrivalControls {
        minimum_step_s: 1e-12,
        time_tolerance_s: 1e-10,
        gap_tolerance_m: 1e-10,
        maximum_trials: 128,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn solve(
    initial_gap: f64,
    upper: f64,
    controls: ArrivalControls,
    response: impl Fn(f64) -> (f64, f64),
) -> (ArrivalReport, Vec<(f64, f64, f64)>) {
    let upper_values = response(upper);
    let mut calls = Vec::new();
    let report = locate(
        initial_gap,
        upper,
        upper_values.0,
        upper_values.1,
        controls,
        &mut |h, _| {
            assert!(h > 0.0, "no native zero-duration solve");
            let (gap, velocity) = response(h);
            calls.push((h, gap, velocity));
            Ok::<_, Infallible>((h, gap, velocity))
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(
        calls.last(),
        Some(&(report.elapsed_s, report.gap_m, report.normal_velocity_m_s))
    );
    assert_eq!(report.response_evaluations, calls.len());
    assert!(report.gap_m >= 0.0 && report.gap_m <= controls.gap_tolerance_m);
    assert!(report.bracket_s[1] - report.bracket_s[0] <= controls.time_tolerance_s);
    assert_eq!(report.elapsed_s, report.bracket_s[0]);
    assert!(
        calls
            .iter()
            .all(|s| s.0 >= controls.minimum_step_s && s.0 <= upper)
    );
    (report, calls)
}

#[test]
fn original_ballistic_endpoint_reaches_an_exact_contact_without_projection() {
    // m=2,F=-8,v0=-1/2,g0=3/32: q=g0-h/2-2h²; arrival h=1/8.
    let mut c = controls();
    c.gap_tolerance_m = 0.0;
    let (arrival, _) = solve(3.0 / 32.0, 0.25, c, |h| {
        (3.0 / 32.0 - 0.5 * h - 2.0 * h * h, -0.5 - 4.0 * h)
    });
    assert_eq!(arrival.elapsed_s, 0.125);
    assert_eq!(arrival.gap_m, 0.0);
    assert_eq!(arrival.normal_velocity_m_s, -1.0);
    assert_eq!(arrival.bracket_s, [0.125; 2]);
    let gravity_work = 8.0 * (3.0 / 32.0);
    close(
        0.25 + gravity_work,
        arrival.normal_velocity_m_s.powi(2),
        1e-15,
    );
    // The original plastic impact then loses 1 J; supported rest needs 8 N.
    close(-2.0 * arrival.normal_velocity_m_s, 2.0, 1e-15);
}

#[test]
fn moving_surface_arrival_keeps_actual_velocity_and_impulse_source_work() {
    // Same falling body, but a wall approaches at +1/2 m/s and g0=5/32.
    let mut c = controls();
    c.gap_tolerance_m = 0.0;
    let (arrival, _) = solve(5.0 / 32.0, 0.25, c, |h| {
        (5.0 / 32.0 - h - 2.0 * h * h, -1.0 - 4.0 * h)
    });
    assert_eq!(arrival.elapsed_s, 0.125);
    assert_eq!(arrival.normal_velocity_m_s, -1.5);
    let before_world = arrival.normal_velocity_m_s + 0.5;
    let impulse = -2.0 * arrival.normal_velocity_m_s;
    let after_world = before_world + impulse / 2.0;
    close(after_world, 0.5, 1e-15);
    let kinetic_change = after_world.powi(2) - before_world.powi(2);
    let loss = arrival.normal_velocity_m_s.powi(2);
    close(0.5 * impulse, kinetic_change + loss, 1e-15);
}

#[test]
fn nonlinear_original_spring_midpoint_response_retains_energy() {
    // m=1,k=4,q0=1/8,v0=-1/2. This is the exact original midpoint response
    // as a function of duration, not a linear extrapolation of its upper gap.
    let response = |h: f64| {
        let q = ((1.0 - h * h) * 0.125 - 0.5 * h) / (1.0 + h * h);
        let v = 2.0 * (q - 0.125) / h + 0.5;
        (q, v)
    };
    let (arrival, _) = solve(0.125, 0.5, controls(), response);
    close(arrival.elapsed_s, 5.0_f64.sqrt() - 2.0, 1e-10);
    close(
        0.5 * arrival.normal_velocity_m_s.powi(2) + 2.0 * arrival.gap_m.powi(2),
        5.0 / 32.0,
        1e-12,
    );
}

#[test]
fn approximate_arrival_is_the_actual_positive_side_not_a_zeroed_gap() {
    // Fused evaluation retains the nonzero difference for this irrational
    // root; neither neighboring floating-point duration is an exact zero.
    let response = |h: f64| ((-h).mul_add(h, 0.02), -2.0 * h);
    let (arrival, _) = solve(0.02, 0.25, controls(), response);
    assert!(arrival.gap_m > 0.0);
    assert!(arrival.bracket_s[0] < arrival.bracket_s[1]);
    close(arrival.elapsed_s, 0.02_f64.sqrt(), 1e-10);
    let mut exact = controls();
    exact.gap_tolerance_m = 0.0;
    let upper = response(0.25);
    let error = locate(
        0.02,
        0.25,
        upper.0,
        upper.1,
        exact,
        &mut |h, _| {
            let (g, v) = response(h);
            Ok::<_, Infallible>((h, g, v))
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        ArrivalError::Unrepresentable { .. } | ArrivalError::Budget { .. }
    ));
}

#[test]
fn another_contact_is_retargeted_from_the_same_original_predecessor() {
    // The outer owner's all-contact check finds B before the initially
    // selected A. Both searches use original gaps, not a rejected A successor.
    let (a, _) = solve(0.25, 0.5, controls(), |h| (0.25 - h, -1.0));
    assert!(0.125 - a.elapsed_s < 0.0);
    let (b, _) = solve(0.125, a.elapsed_s, controls(), |h| (0.125 - h, -1.0));
    assert_eq!(b.elapsed_s, 0.125);
    assert!(0.25 - b.elapsed_s > 0.0);
}

#[test]
fn grazing_endpoint_retains_zero_normal_speed_and_departure_is_not_arrival() {
    let (grazing, calls) = solve(1.0, 1.0, controls(), |h| {
        ((1.0 - h).powi(2), -2.0 * (1.0 - h))
    });
    assert_eq!(grazing.normal_velocity_m_s, 0.0);
    assert_eq!(grazing.bracket_s, [1.0; 2]);
    assert_eq!(
        calls.len(),
        1,
        "the supplied exact upper endpoint only needs replay"
    );
    let mut count = 0;
    let departure = locate(
        0.0,
        1.0,
        -0.1,
        -1.0,
        controls(),
        &mut |h, _| {
            count += 1;
            Ok::<_, Infallible>((h, 0.0, 1.0))
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(departure, ArrivalError::Invalid(_)));
    assert_eq!(count, 0);
    let unbracketed = locate(
        0.1,
        1.0,
        1e-12,
        -1.0,
        controls(),
        &mut |h, _| Ok::<_, Infallible>((h, 0.0, -1.0)),
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(unbracketed, ArrivalError::NotBracketed { .. }));
}

#[test]
fn earlier_physical_endpoint_marker_is_terminal_and_keeps_its_native_cache() {
    #[derive(Debug, PartialEq)]
    enum Probe {
        EarlierEndpoint,
    }
    let mut cached = None;
    let mut calls = 0;
    let error = locate(
        0.125,
        0.25,
        -0.125,
        -1.0,
        controls(),
        &mut |_, _| {
            calls += 1;
            cached = Some((0.0625, 0.0625, -1.0));
            Err::<(f64, f64, f64), _>(Probe::EarlierEndpoint)
        },
        &mut || false,
    )
    .unwrap_err();
    assert_eq!(calls, 1);
    assert_eq!(cached, Some((0.0625, 0.0625, -1.0)));
    assert!(matches!(
        error,
        ArrivalError::Response(Probe::EarlierEndpoint)
    ));
    let mapped = ArrivalError::Response("original refusal").map_response(str::len);
    assert!(matches!(mapped, ArrivalError::Response(16)));
    let preserved = ArrivalError::<Infallible>::Budget {
        evaluations: 3,
        maximum: 3,
    }
    .map_response(|never| match never {});
    assert!(matches!(
        preserved,
        ArrivalError::Budget {
            evaluations: 3,
            maximum: 3
        }
    ));
}

#[test]
fn replay_rejects_different_gap_inside_tolerance_or_different_velocity() {
    let response = |h: f64| ((-h).mul_add(h, 0.02), -2.0 * h);
    let upper = response(0.25);
    let mut seen = Vec::new();
    let changed_gap = locate(
        0.02,
        0.25,
        upper.0,
        upper.1,
        controls(),
        &mut |h, _| {
            let (g, v) = response(h);
            let repeated = seen.contains(&h);
            seen.push(h);
            Ok::<_, Infallible>((h, if repeated { 0.5 * g } else { g }, v))
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(
        matches!(changed_gap, ArrivalError::Replay { gap_m, .. } if gap_m > 0.0 && gap_m < controls().gap_tolerance_m)
    );
    let changed_velocity = locate(
        0.125,
        0.125,
        0.0,
        -1.0,
        controls(),
        &mut |h, _| Ok::<_, Infallible>((h, 0.0, -2.0)),
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(changed_velocity, ArrivalError::Replay { .. }));
    let changed_time = locate(
        0.125,
        0.125,
        0.0,
        -1.0,
        controls(),
        &mut |h, _| Ok::<_, Infallible>((h + 0.001, 0.0, -1.0)),
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(changed_time, ArrivalError::Replay { .. }));
    let negative_gap = locate(
        0.125,
        0.125,
        0.0,
        -1.0,
        controls(),
        &mut |h, _| Ok::<_, Infallible>((h, -1e-30, -1.0)),
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(negative_gap, ArrivalError::Replay { .. }));
}

#[test]
fn represented_clock_and_minimum_duration_are_physical_limits() {
    let quantized = locate(
        0.2,
        0.5,
        -0.3,
        -1.0,
        controls(),
        &mut |h, _| {
            let actual = (8.0 * h).round() / 8.0;
            Ok::<_, Infallible>((actual, 0.2 - actual, -1.0))
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(quantized, ArrivalError::Unrepresentable { .. }));
    let mut c = controls();
    c.minimum_step_s = 0.2;
    let minimum = locate(
        0.1,
        0.5,
        -0.4,
        -1.0,
        c,
        &mut |h, _| {
            assert!(h >= 0.2);
            Ok::<_, Infallible>((h, 0.1 - h, -1.0))
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(minimum, ArrivalError::MinimumStep { .. }));
}

#[test]
fn replay_is_charged_and_cancellation_and_source_failures_do_not_publish() {
    let mut c = controls();
    c.maximum_trials = 1;
    let mut calls = 0;
    let budget = locate(
        0.125,
        0.25,
        -0.125,
        -1.0,
        c,
        &mut |h, _| {
            calls += 1;
            Ok::<_, Infallible>((h, 0.125 - h, -1.0))
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(
        budget,
        ArrivalError::Budget {
            evaluations: 1,
            maximum: 1
        }
    ));
    assert_eq!(calls, 1);
    let mut polls = 0;
    let cancellation = locate(
        0.125,
        0.25,
        -0.125,
        -1.0,
        controls(),
        &mut |h, cancelled| {
            for _ in 0..8 {
                if cancelled() {
                    break;
                }
            }
            Ok::<_, Infallible>((h, 0.0, -1.0))
        },
        &mut || {
            polls += 1;
            polls > 3
        },
    )
    .unwrap_err();
    assert!(matches!(cancellation, ArrivalError::Cancelled));
    let refused = locate(
        0.125,
        0.25,
        -0.125,
        -1.0,
        controls(),
        &mut |_, _| Err::<(f64, f64, f64), _>("original material refusal"),
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(
        refused,
        ArrivalError::Response("original material refusal")
    ));
    let panic = locate::<Infallible>(
        0.125,
        0.25,
        -0.125,
        -1.0,
        controls(),
        &mut |_, _| panic!("original response panic"),
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(panic, ArrivalError::ResponsePanicked));
}
