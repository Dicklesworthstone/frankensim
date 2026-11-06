//! G0/G1/G3/G4/G5 guard checks against exact scalar trajectories and roots.

use super::{
    CrossingDirection, GuardBounds, GuardError, GuardOptions, GuardScan, GuardWindowKind,
    scan_guard,
};
use crate::validated::trajectory::{TrajectoryOptions, TrajectoryTube, validated_trajectory};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey};
use fs_ivl::Interval;

fn with_cx<R>(budget: Budget, f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn linear_tube(cx: &Cx<'_>, initial: Interval, rate: f64, step: f64) -> TrajectoryTube {
    validated_trajectory(
        cx,
        0.0,
        1.0,
        &[initial],
        &mut |_, _, out, _| {
            out[0] = Interval::point(rate);
            Ok(())
        },
        &TrajectoryOptions {
            initial_step: step,
            ..TrajectoryOptions::default()
        },
    )
    .unwrap()
}

fn options() -> GuardOptions {
    GuardOptions {
        min_width: 1.0 / 128.0,
        ..GuardOptions::default()
    }
}

fn covers(scan: &GuardScan, time: f64) -> bool {
    scan.windows().iter().any(|window| {
        let [start, end] = window.interval();
        start <= time && time <= end
    })
}

#[test]
fn g1_strict_monotone_brackets_localize_rising_and_falling_crossings() {
    with_cx(Budget::INFINITE, |cx, _| {
        for (rate, direction) in [
            (1.0, CrossingDirection::Rising),
            (-1.0, CrossingDirection::Falling),
        ] {
            let tube = linear_tube(cx, Interval::point(0.0), rate, 1.0);
            let threshold = rate / 3.0;
            let result = scan_guard(
                cx,
                &tube,
                &mut |_, state, _| {
                    Ok(GuardBounds::with_derivative(
                        state[0] - Interval::point(threshold),
                        Interval::point(rate),
                    ))
                },
                &options(),
            )
            .unwrap();
            assert_eq!(result.interval(), [0.0, 1.0]);
            assert!(!result.no_events_proved());
            let crossing = result
                .windows()
                .iter()
                .find(|window| {
                    let [start, end] = window.interval();
                    start <= 1.0 / 3.0 && 1.0 / 3.0 <= end
                })
                .expect("the exact crossing must remain covered");
            assert_eq!(
                crossing.kind(),
                GuardWindowKind::UniqueCrossing { direction }
            );
            assert!(crossing.interval()[1] - crossing.interval()[0] <= options().min_width);
            assert!(crossing.value().contains(0.0));
            assert_eq!(crossing.step_index(), Some(0));
        }
    });
}

#[test]
fn g1_sign_separated_whole_tube_proves_no_event_with_one_guard_call() {
    with_cx(Budget::INFINITE, |cx, _| {
        let tube = linear_tube(cx, Interval::new(-0.25, 0.25), 1.0, 1.0);
        let mut calls = 0;
        let result = scan_guard(
            cx,
            &tube,
            &mut |_, state, _| {
                calls += 1;
                Ok(GuardBounds::value(state[0] + Interval::point(2.0)))
            },
            &GuardOptions {
                max_cells: 1,
                max_guard_evaluations: 1,
                ..options()
            },
        )
        .unwrap();
        assert!(result.no_events_proved());
        assert!(result.windows().is_empty());
        assert_eq!(result.interval(), [0.0, 1.0]);
        assert_eq!(result.cells(), 1);
        assert_eq!(result.guard_evaluations(), 1);
        assert_eq!(calls, 1);
    });
}

#[test]
fn g1_monotone_endpoints_exclude_a_loose_whole_cell_guard_bound() {
    with_cx(Budget::INFINITE, |cx, _| {
        let tube = linear_tube(cx, Interval::point(0.0), 0.0, 1.0);
        let mut calls = 0;
        let result = scan_guard(
            cx,
            &tube,
            &mut |time, _, _| {
                calls += 1;
                // The guard is g(t)=t+1. Its deliberately loose interval
                // extension straddles zero, but monotonicity and exact
                // positive endpoint values exclude a crossing.
                let value = if time.lo() == time.hi() {
                    Interval::point(time.lo() + 1.0)
                } else {
                    Interval::new(-2.0, 3.0)
                };
                Ok(GuardBounds::with_derivative(value, Interval::point(1.0)))
            },
            &GuardOptions {
                min_width: 1.0,
                ..options()
            },
        )
        .unwrap();
        assert!(result.no_events_proved());
        assert!(result.windows().is_empty());
        assert_eq!(result.cells(), 1);
        assert_eq!(result.guard_evaluations(), 3);
        assert_eq!(calls, 3);
    });
}

#[test]
fn g1_roots_missed_by_endpoint_and_midpoint_samples_remain_covered() {
    with_cx(Budget::INFINITE, |cx, _| {
        let tube = linear_tube(cx, Interval::point(0.0), 0.0, 1.0);
        for (left_root, right_root) in [(0.375, 0.375), (0.1875, 0.3125)] {
            // A tangency and a pair of crossings all have positive values at
            // the endpoints and midpoint of this integration step.
            for sample in [0.0, 0.5, 1.0] {
                assert!((sample - left_root) * (sample - right_root) > 0.0);
            }
            let result = scan_guard(
                cx,
                &tube,
                &mut |time, _, _| {
                    Ok(GuardBounds::with_derivative(
                        (time - Interval::point(left_root)) * (time - Interval::point(right_root)),
                        Interval::point(2.0) * time - Interval::point(left_root + right_root),
                    ))
                },
                &options(),
            )
            .unwrap();
            assert!(!result.no_events_proved());
            assert!(covers(&result, left_root));
            assert!(covers(&result, right_root));
            // These roots lie on exact binary subdivision boundaries; the
            // strict-sign bracket theorem must not invent an interior crossing.
            for window in result.windows() {
                assert_eq!(window.kind(), GuardWindowKind::Unresolved);
            }
        }
    });
}

#[test]
fn g1_initial_uncertainty_retains_the_full_family_of_possible_root_times() {
    with_cx(Budget::INFINITE, |cx, _| {
        let tube = linear_tube(cx, Interval::new(-0.125, 0.125), 1.0, 1.0);
        let result = scan_guard(
            cx,
            &tube,
            &mut |_, state, _| {
                Ok(GuardBounds::with_derivative(
                    state[0] - Interval::point(0.5),
                    Interval::point(1.0),
                ))
            },
            &options(),
        )
        .unwrap();
        assert!(!result.no_events_proved());
        // x(t)=x0+t gives t_event=0.5-x0 for every admitted initial value.
        for root in [0.375, 0.4375, 0.5, 0.5625, 0.625] {
            assert!(covers(&result, root), "lost possible root at {root}");
        }
        assert!(
            result
                .windows()
                .iter()
                .all(|window| window.kind() == GuardWindowKind::Unresolved)
        );
    });
}

#[test]
fn g1_event_at_a_shared_step_endpoint_is_not_discarded() {
    with_cx(Budget::INFINITE, |cx, _| {
        let tube = linear_tube(cx, Interval::point(0.0), 1.0, 0.5);
        assert_eq!(tube.steps().len(), 2);
        assert_eq!(tube.steps()[0].interval()[1], 0.5);
        let result = scan_guard(
            cx,
            &tube,
            &mut |_, state, _| {
                Ok(GuardBounds::with_derivative(
                    state[0] - Interval::point(0.5),
                    Interval::point(1.0),
                ))
            },
            &options(),
        )
        .unwrap();
        assert!(!result.no_events_proved());
        assert!(covers(&result, 0.5));
        assert!(
            result
                .windows()
                .iter()
                .any(|window| { window.step_index() == Some(0) && window.interval()[1] == 0.5 })
        );
        assert!(
            result
                .windows()
                .iter()
                .any(|window| { window.step_index() == Some(1) && window.interval()[0] == 0.5 })
        );
    });
}

#[test]
fn g0_zero_duration_checks_the_initial_guard_instead_of_proving_a_vacuous_miss() {
    with_cx(Budget::INFINITE, |cx, _| {
        let tube = validated_trajectory(
            cx,
            2.0,
            2.0,
            &[Interval::point(0.0)],
            &mut |_, _, _, _| panic!("a zero-duration trajectory must not evaluate its RHS"),
            &TrajectoryOptions::default(),
        )
        .unwrap();
        for offset in [0.0, 1.0] {
            let mut calls = 0;
            let result = scan_guard(
                cx,
                &tube,
                &mut |time, state, _| {
                    calls += 1;
                    assert_eq!(time, Interval::point(2.0));
                    Ok(GuardBounds::value(state[0] + Interval::point(offset)))
                },
                &options(),
            )
            .unwrap();
            assert_eq!(calls, 1);
            assert_eq!(result.interval(), [2.0, 2.0]);
            assert_eq!(result.no_events_proved(), offset != 0.0);
            if offset == 0.0 {
                assert!(covers(&result, 2.0));
                assert_eq!(result.windows()[0].kind(), GuardWindowKind::Unresolved);
                assert_eq!(result.windows()[0].step_index(), None);
            }
        }
    });
}

#[test]
fn g0_adjacent_representable_times_terminate_with_an_unresolved_window() {
    with_cx(Budget::INFINITE, |cx, _| {
        let start = 1.0;
        let end = fs_math::next_up(start);
        let tube = validated_trajectory(
            cx,
            start,
            end,
            &[Interval::point(0.0)],
            &mut |_, _, out, _| {
                out[0] = Interval::point(0.0);
                Ok(())
            },
            &TrajectoryOptions::default(),
        )
        .unwrap();
        let result = scan_guard(
            cx,
            &tube,
            &mut |time, _, _| Ok(GuardBounds::value(time - Interval::point(start))),
            &GuardOptions {
                min_width: f64::from_bits(1),
                ..options()
            },
        )
        .unwrap();
        assert!(!result.no_events_proved());
        assert_eq!(result.windows().len(), 1);
        assert_eq!(result.windows()[0].interval(), [start, end]);
        assert_eq!(result.windows()[0].kind(), GuardWindowKind::Unresolved);
    });
}

#[test]
fn g0_exact_flat_guard_retains_one_unresolved_full_window() {
    with_cx(Budget::INFINITE, |cx, _| {
        let tube = linear_tube(cx, Interval::point(0.0), 0.0, 1.0);
        let result = scan_guard(
            cx,
            &tube,
            &mut |_, _, _| Ok(GuardBounds::value(Interval::point(0.0))),
            &GuardOptions {
                max_cells: 1,
                max_guard_evaluations: 1,
                ..options()
            },
        )
        .unwrap();
        assert!(!result.no_events_proved());
        assert_eq!(result.cells(), 1);
        assert_eq!(result.guard_evaluations(), 1);
        assert_eq!(result.windows().len(), 1);
        assert_eq!(result.windows()[0].interval(), [0.0, 1.0]);
        assert_eq!(result.windows()[0].kind(), GuardWindowKind::Unresolved);
    });
}

#[test]
fn g0_search_and_output_limits_never_publish_partial_no_event_evidence() {
    with_cx(Budget::INFINITE, |cx, _| {
        let tube = linear_tube(cx, Interval::point(0.0), 0.0, 1.0);
        for bounded in [
            GuardOptions {
                max_cells: 1,
                ..options()
            },
            GuardOptions {
                max_guard_evaluations: 1,
                ..options()
            },
            GuardOptions {
                max_pending_cells: 1,
                ..options()
            },
            GuardOptions {
                max_windows: 1,
                ..options()
            },
            GuardOptions {
                max_work: 1,
                ..options()
            },
        ] {
            let mut calls = 0;
            let result = scan_guard(
                cx,
                &tube,
                &mut |_, _, _| {
                    calls += 1;
                    Ok(GuardBounds::value(Interval::new(-1.0, 1.0)))
                },
                &bounded,
            );
            assert!(result.is_err(), "an incomplete search published {result:?}");
            assert!(calls <= bounded.max_guard_evaluations);
        }
    });
}

#[test]
fn g0_exact_work_budget_is_shared_across_all_trajectory_steps() {
    with_cx(Budget::INFINITE, |cx, _| {
        let tube = linear_tube(cx, Interval::point(0.0), 1.0, 0.25);
        assert_eq!(tube.steps().len(), 4);
        let mut guard = |time: Interval, _: &[Interval], _: &Cx<'_>| {
            Ok(GuardBounds::value(time + Interval::point(1.0)))
        };
        let baseline = scan_guard(cx, &tube, &mut guard, &options()).unwrap();
        assert!(baseline.no_events_proved());
        let exact = GuardOptions {
            max_work: baseline.work(),
            ..options()
        };
        assert_eq!(scan_guard(cx, &tube, &mut guard, &exact).unwrap(), baseline);
        let short = GuardOptions {
            max_work: baseline.work() - 1,
            ..options()
        };
        assert_eq!(
            scan_guard(cx, &tube, &mut guard, &short).unwrap_err(),
            GuardError::WorkLimit {
                limit: short.max_work
            }
        );
    });
}

#[test]
fn g0_invalid_or_refused_guard_outputs_cannot_be_mistaken_for_a_miss() {
    with_cx(Budget::INFINITE, |cx, _| {
        let tube = linear_tube(cx, Interval::point(0.0), 0.0, 1.0);
        for (bounds, expected) in [
            (
                GuardBounds::value(Interval::WHOLE),
                GuardError::NonFiniteGuard,
            ),
            (
                GuardBounds::with_derivative(Interval::point(1.0), Interval::WHOLE),
                GuardError::NonFiniteDerivative,
            ),
        ] {
            let result = scan_guard(cx, &tube, &mut |_, _, _| Ok(bounds), &options());
            assert_eq!(result.unwrap_err(), expected);
        }
        let mut calls = 0;
        let result = scan_guard(
            cx,
            &tube,
            &mut |_, _, _| {
                calls += 1;
                Err(GuardError::GuardRefused("outside guard model domain"))
            },
            &options(),
        );
        assert_eq!(
            result.unwrap_err(),
            GuardError::GuardRefused("outside guard model domain")
        );
        assert_eq!(calls, 1);
    });
}

#[test]
fn g4_cancellation_before_or_inside_a_guard_never_returns_a_certificate() {
    let tube = with_cx(Budget::INFINITE, |cx, _| {
        linear_tube(cx, Interval::point(0.0), 0.0, 1.0)
    });
    for before in [true, false] {
        with_cx(Budget::INFINITE, |cx, gate| {
            if before {
                gate.request();
            }
            let mut calls = 0;
            let result = scan_guard(
                cx,
                &tube,
                &mut |_, _, _| {
                    calls += 1;
                    gate.request();
                    Ok(GuardBounds::value(Interval::point(1.0)))
                },
                &options(),
            );
            assert_eq!(result.unwrap_err(), GuardError::Cancelled);
            assert_eq!(calls, usize::from(!before));
        });
    }
}

#[test]
fn g4_ambient_cost_poll_and_clock_requirements_are_enforced_before_guard_work() {
    let tube = with_cx(Budget::INFINITE, |cx, _| {
        linear_tube(cx, Interval::point(0.0), 0.0, 1.0)
    });
    for budget in [
        Budget {
            cost_quota: Some(1),
            ..Budget::INFINITE
        },
        Budget {
            poll_quota: 0,
            ..Budget::INFINITE
        },
        Budget::ZERO,
    ] {
        with_cx(budget, |cx, _| {
            let result = scan_guard(
                cx,
                &tube,
                &mut |_, _, _| panic!("refused execution budget evaluated the guard"),
                &options(),
            );
            assert!(matches!(
                result,
                Err(GuardError::Budget(
                    BudgetRefusal::CostPlanExceedsQuota { .. }
                        | BudgetRefusal::PollsExhausted { .. }
                        | BudgetRefusal::DeadlineWithoutClock { .. }
                ))
            ));
        });
    }
}

#[test]
fn g3_refinement_narrows_a_regular_root_and_g5_repeats_match() {
    with_cx(Budget::INFINITE, |cx, _| {
        let tube = linear_tube(cx, Interval::point(0.0), 1.0, 1.0);
        let mut guard = |_: Interval, state: &[Interval], _: &Cx<'_>| {
            Ok(GuardBounds::with_derivative(
                state[0] - Interval::point(1.0 / 3.0),
                Interval::point(1.0),
            ))
        };
        let coarse = scan_guard(
            cx,
            &tube,
            &mut guard,
            &GuardOptions {
                min_width: 1.0 / 32.0,
                ..options()
            },
        )
        .unwrap();
        let fine = scan_guard(cx, &tube, &mut guard, &options()).unwrap();
        let repeat = scan_guard(cx, &tube, &mut guard, &options()).unwrap();
        assert_eq!(coarse.windows().len(), 1);
        assert_eq!(fine.windows().len(), 1);
        let coarse_interval = coarse.windows()[0].interval();
        let fine_interval = fine.windows()[0].interval();
        assert!(coarse_interval[0] <= fine_interval[0]);
        assert!(fine_interval[1] <= coarse_interval[1]);
        assert!(fine_interval[1] - fine_interval[0] < coarse_interval[1] - coarse_interval[0]);
        assert!(covers(&coarse, 1.0 / 3.0));
        assert!(covers(&fine, 1.0 / 3.0));
        assert_eq!(fine.interval(), repeat.interval());
        assert_eq!(fine.cells(), repeat.cells());
        assert_eq!(fine.guard_evaluations(), repeat.guard_evaluations());
        assert_eq!(fine.work(), repeat.work());
        assert_eq!(fine.windows().len(), repeat.windows().len());
        for (left, right) in fine.windows().iter().zip(repeat.windows()) {
            assert_eq!(left.interval(), right.interval());
            assert_eq!(left.step_index(), right.step_index());
            assert_eq!(left.value(), right.value());
            assert_eq!(left.kind(), right.kind());
        }
    });
}
