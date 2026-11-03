use super::*;
use std::convert::Infallible;

fn controls() -> StoppingControls {
    StoppingControls {
        minimum_step_s: 1e-12,
        time_tolerance_s: 1e-10,
        velocity_tolerance_m_s: 1e-9,
        maximum_trials: 128,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn solve(
    initial: [f64; 2],
    upper: f64,
    velocity: impl Fn(f64) -> [f64; 2],
) -> (StoppingReport, Vec<f64>) {
    let mut calls = Vec::new();
    let result = locate(
        initial,
        upper,
        velocity(upper),
        controls(),
        &mut |h, _| {
            calls.push(h);
            Ok::<_, Infallible>((h, velocity(h)))
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(calls.last().copied(), Some(result.elapsed_s));
    assert_eq!(result.response_evaluations, calls.len());
    assert_eq!(result.velocity_m_s, velocity(result.elapsed_s));
    assert!(
        result.velocity_m_s[0].hypot(result.velocity_m_s[1]) <= controls().velocity_tolerance_m_s
    );
    assert!(result.bracket_s[1] - result.bracket_s[0] <= controls().time_tolerance_s);
    assert!(
        calls
            .iter()
            .all(|h| *h >= controls().minimum_step_s && *h <= upper)
    );
    (result, calls)
}

#[test]
fn coulomb_stop_retains_actual_endpoint_and_exact_work() {
    // m=2, N=20, mu=1/4, v0=3: a=-5/2, stop=6/5, travel=9/5.
    let (stop, calls) = solve([3.0, 0.0], 1.5, |h| [3.0 - 2.5 * h, 0.0]);
    close(stop.elapsed_s, 1.2, 1e-10);
    let travel = 0.5 * stop.elapsed_s * (3.0 + stop.velocity_m_s[0]);
    close(travel, 1.8, 1e-9);
    close(5.0 * travel, 9.0, 1e-8);
    assert!(
        calls.len() < 12,
        "secant should resolve this linear response directly"
    );
}

#[test]
fn long_midpoint_reflection_does_not_become_a_false_zero_slip_stop() {
    // The integrated-slip Coulomb law changes to its interior disk after
    // h=2*t_stop, giving v1=-v0 with zero net slip. It is a private probe only.
    let (stop, _) = solve([3.0, 0.0], 4.0, |h| {
        [if h <= 2.4 { 3.0 - 2.5 * h } else { -3.0 }, 0.0]
    });
    close(stop.elapsed_s, 1.2, 1e-10);
    assert!(stop.elapsed_s < 2.4);
    close(stop.velocity_m_s[0], 0.0, 1e-9);
}

#[test]
fn moving_surface_stop_is_relative_and_internal_pair_preserves_momentum() {
    // A 2 m/s wall brings a stationary m=2 body to its speed under +5 N.
    let (wall, _) = solve([-2.0, 0.0], 1.0, |h| [-2.0 + 2.5 * h, 0.0]);
    close(wall.elapsed_s, 0.8, 1e-10);
    let body_speed = 2.5 * wall.elapsed_s;
    let body_travel = 0.5 * body_speed * wall.elapsed_s;
    let relative_travel = body_travel - 2.0 * wall.elapsed_s;
    close(body_speed, 2.0, 1e-9);
    close(
        5.0 * 2.0 * wall.elapsed_s,
        body_speed * body_speed - 5.0 * relative_travel,
        1e-8,
    );
    // mA=1,mB=2; vA=0,vB=3; equal/opposite 5 N friction preserves p=6.
    let (pair, _) = solve([3.0, 0.0], 0.6, |h| [3.0 - 7.5 * h, 0.0]);
    close(pair.elapsed_s, 0.4, 1e-10);
    let va = 5.0 * pair.elapsed_s;
    let vb = 3.0 - 2.5 * pair.elapsed_s;
    close(va, 2.0, 1e-9);
    close(vb, 2.0, 1e-9);
    close(va + 2.0 * vb, 6.0, 1e-12);
    close(9.0 - (0.5 * va * va + vb * vb), 3.0, 1e-8);
}

#[test]
fn full_inertia_and_viscous_midpoint_responses_locate_their_own_stop() {
    // M=[[2,1],[1,2]], fixed normal, Fn=-20 and mu=1/4 imply
    // N=160/9, Ft=-40/9, a_t=-20/9. A diagonal-mass guess is wrong.
    let (coupled, _) = solve([3.0, 0.0], 1.7, |h| [3.0 - (20.0 / 9.0) * h, 0.0]);
    close(coupled.elapsed_s, 27.0 / 20.0, 1e-10);
    let travel = 0.5 * coupled.elapsed_s * (3.0 + coupled.velocity_m_s[0]);
    close(travel, 81.0 / 40.0, 1e-9);
    close((40.0 / 9.0) * travel, 9.0, 1e-8);
    // Original linear viscous midpoint response: m=2,c=1,Fc=5,v0=3.
    // (4+h)*v1=12-13h, so h_stop=12/13. Friction+viscosity=initial K.
    let (viscous, _) = solve([3.0, 0.0], 2.0, |h| [(12.0 - 13.0 * h) / (4.0 + h), 0.0]);
    close(viscous.elapsed_s, 12.0 / 13.0, 1e-10);
    let mean_speed = 0.5 * (3.0 + viscous.velocity_m_s[0]);
    close(
        5.0 * viscous.elapsed_s * mean_speed + viscous.elapsed_s * mean_speed * mean_speed,
        9.0,
        1e-8,
    );
}

#[test]
fn approximate_zero_requires_a_real_narrow_bracket_and_keeps_residual_velocity() {
    let epsilon = 1e-10;
    let (stop, _) = solve([1.0, 0.0], 1.0 + epsilon, |h| [1.0 - h, epsilon * h]);
    assert!(stop.bracket_s[0] < stop.bracket_s[1]);
    assert!(
        stop.velocity_m_s[1] > 0.0,
        "the original transverse residue must not be zeroed"
    );
    close(stop.elapsed_s, 1.0, 1e-10);
    let mut calls = 0;
    let error = locate(
        [1.0, 0.0],
        1.0,
        [epsilon, 0.0],
        controls(),
        &mut |h, _| {
            calls += 1;
            Ok::<_, Infallible>((h, [0.0; 2]))
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(error, StoppingError::NotBracketed { .. }));
    assert_eq!(calls, 0);
}

#[test]
fn transverse_direction_change_is_not_a_stop() {
    let mut calls = 0;
    let error = locate(
        [1.0, 0.0],
        2.0,
        [-1.0, 1.0],
        controls(),
        &mut |h, _| {
            calls += 1;
            Ok::<_, Infallible>((h, [1.0 - h, 0.5 * h]))
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        StoppingError::TransverseVelocity {
            velocity_m_s: [0.0, 0.5],
            ..
        }
    ));
    assert_eq!(calls, 1);
}

#[test]
fn native_clock_quantization_and_nonrepeatable_replay_refuse() {
    let error = locate(
        [1.0, 0.0],
        1.0,
        [-1.5, 0.0],
        controls(),
        &mut |h, _| {
            let actual = (4.0 * h).round() / 4.0;
            Ok::<_, Infallible>((actual, [1.0 - actual / 0.4, 0.0]))
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(error, StoppingError::Unrepresentable { .. }));
    let mut calls = 0;
    let replay = locate(
        [1.0, 0.0],
        1.0,
        [-1.0, 0.0],
        controls(),
        &mut |h, _| {
            calls += 1;
            Ok::<_, Infallible>((if calls == 2 { h + 0.001 } else { h }, [0.0; 2]))
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(replay, StoppingError::Replay { .. }));
    assert_eq!(calls, 2);
}

#[test]
fn budget_charges_final_replay_and_minimum_duration_is_not_crossed() {
    let mut c = controls();
    c.maximum_trials = 1;
    let mut calls = 0;
    let exhausted = locate(
        [1.0, 0.0],
        1.0,
        [-1.0, 0.0],
        c,
        &mut |h, _| {
            calls += 1;
            Ok::<_, Infallible>((h, [1.0 - 2.0 * h, 0.0]))
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(
        exhausted,
        StoppingError::Budget {
            evaluations: 1,
            maximum: 1
        }
    ));
    assert_eq!(calls, 1);
    let mut c = controls();
    c.minimum_step_s = 0.75;
    let below = locate(
        [1.0, 0.0],
        1.0,
        [-1.0, 0.0],
        c,
        &mut |h, _| {
            assert!(h >= 0.75);
            Ok::<_, Infallible>((h, [1.0 - 2.0 * h, 0.0]))
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(below, StoppingError::MinimumStep { .. }));
}

#[test]
fn exact_upper_stop_is_replayed_and_cancellation_cannot_publish_it() {
    let mut calls = 0;
    let upper = locate(
        [1.0, 0.0],
        0.5,
        [0.0; 2],
        controls(),
        &mut |h, _| {
            calls += 1;
            Ok::<_, Infallible>((h, [0.0; 2]))
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(upper.bracket_s, [0.5; 2]);
    assert_eq!(calls, 1);
    let mut polls = 0;
    let stopped = locate(
        [1.0, 0.0],
        1.0,
        [-1.0, 0.0],
        controls(),
        &mut |h, cancelled| {
            for _ in 0..8 {
                if cancelled() {
                    break;
                }
            }
            Ok::<_, Infallible>((h, [0.0; 2]))
        },
        &mut || {
            polls += 1;
            polls > 3
        },
    )
    .unwrap_err();
    assert!(matches!(stopped, StoppingError::Cancelled));
}

#[test]
fn first_original_error_or_unwind_ends_the_search() {
    let mut calls = 0;
    let refused = locate(
        [1.0, 0.0],
        1.0,
        [-1.0, 0.0],
        controls(),
        &mut |_, _| {
            calls += 1;
            Err::<(f64, [f64; 2]), _>("native material refusal")
        },
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(
        refused,
        StoppingError::Response("native material refusal")
    ));
    assert_eq!(calls, 1);
    let panicked = locate::<Infallible>(
        [1.0, 0.0],
        1.0,
        [-1.0, 0.0],
        controls(),
        &mut |_, _| panic!("original response panic"),
        &mut || false,
    )
    .unwrap_err();
    assert!(matches!(panicked, StoppingError::ResponsePanicked));
}
