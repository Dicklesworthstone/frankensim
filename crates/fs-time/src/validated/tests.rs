//! G0/G1/G3/G4/G5 tests against analytic flows and exact resource counts.

use super::{ValidatedError, ValidatedOptions, validate_state, validated_step};
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
                seed: 37,
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

#[test]
fn g1_constant_flow_preserves_initial_uncertainty_and_subwindow_ranges() {
    with_cx(Budget::INFINITE, |cx, _| {
        let initial = [Interval::new(-2.0, -1.0), Interval::new(3.0, 4.0)];
        let step = validated_step(
            cx,
            2.0,
            2.25,
            &initial,
            &mut |_, _, out, _| {
                out[0] = Interval::point(2.0);
                out[1] = Interval::point(-3.0);
                Ok(())
            },
            &ValidatedOptions::default(),
        )
        .unwrap();
        assert_eq!(step.initial(), &initial);
        assert_eq!(step.interval(), [2.0, 2.25]);
        assert!(step.endpoint()[0].encloses(Interval::new(-1.5, -0.5)));
        assert!(step.endpoint()[1].encloses(Interval::new(2.25, 3.25)));
        assert_eq!(step.range(Interval::point(2.0), cx).unwrap(), initial);
        assert_eq!(
            step.range(Interval::point(2.25), cx).unwrap(),
            step.endpoint()
        );
        let middle = step.range(Interval::point(2.125), cx).unwrap();
        assert!(middle[0].encloses(Interval::new(-1.75, -0.75)));
        assert!(middle[1].encloses(Interval::new(2.625, 3.625)));
        let window = step.range(Interval::new(2.0625, 2.1875), cx).unwrap();
        assert!(window[0].encloses(Interval::new(-1.875, -0.625)));
        assert!(window[1].encloses(Interval::new(2.4375, 3.8125)));
        for (whole, end) in step.enclosure().iter().zip(step.endpoint()) {
            assert!(whole.encloses(*end));
        }
        assert!(step.work() > 0);
        assert_eq!(step.iterations(), step.rhs_evaluations());
    });
}

#[test]
fn g1_exponential_flows_enclosed_and_g3_step_refinement_reduces_wrapping() {
    with_cx(Budget::INFINITE, |cx, _| {
        for rate in [-1.0, 1.0] {
            let mut widths = Vec::new();
            for h in [0.2, 0.1, 0.05] {
                let step = validated_step(
                    cx,
                    0.0,
                    h,
                    &[Interval::point(1.0)],
                    &mut |_, x, out, _| {
                        out[0] = Interval::point(rate) * x[0];
                        Ok(())
                    },
                    &ValidatedOptions::default(),
                )
                .unwrap();
                assert!(step.endpoint()[0].contains((rate * h).exp()));
                for fraction in [0.0, 0.25, 0.5, 0.75, 1.0] {
                    let t = fraction * h;
                    let range = step.range(Interval::point(t), cx).unwrap();
                    assert!(range[0].contains((rate * t).exp()));
                    assert!(step.enclosure()[0].encloses(range[0]));
                }
                widths.push(step.endpoint()[0].width());
            }
            assert!(widths[1] < widths[0] * 0.3, "{widths:?}");
            assert!(widths[2] < widths[1] * 0.3, "{widths:?}");
        }
    });
}

#[test]
fn g1_coupled_oscillator_and_nonlinear_flow_are_enclosed() {
    with_cx(Budget::INFINITE, |cx, _| {
        let oscillator = validated_step(
            cx,
            0.0,
            0.125,
            &[Interval::point(1.0), Interval::point(0.0)],
            &mut |_, x, out, _| {
                out[0] = x[1];
                out[1] = -x[0];
                Ok(())
            },
            &ValidatedOptions::default(),
        )
        .unwrap();
        assert!(oscillator.endpoint()[0].contains(0.125_f64.cos()));
        assert!(oscillator.endpoint()[1].contains(-0.125_f64.sin()));
        let nonlinear = validated_step(
            cx,
            0.0,
            0.125,
            &[Interval::new(0.5, 1.0)],
            &mut |_, x, out, _| {
                out[0] = x[0] * x[0];
                Ok(())
            },
            &ValidatedOptions::default(),
        )
        .unwrap();
        // x(t) = x0 / (1 - t*x0), monotone over this admitted initial set.
        assert!(nonlinear.endpoint()[0].contains(0.5 / (1.0 - 0.125 * 0.5)));
        assert!(nonlinear.endpoint()[0].contains(1.0 / (1.0 - 0.125)));
    });
}

#[test]
fn g1_nonautonomous_rhs_receives_the_full_absolute_time_interval() {
    with_cx(Budget::INFINITE, |cx, _| {
        let step = validated_step(
            cx,
            3.0,
            3.25,
            &[Interval::point(1.0)],
            &mut |time, _, out, _| {
                assert_eq!(time, Interval::new(3.0, 3.25));
                out[0] = time;
                Ok(())
            },
            &ValidatedOptions::default(),
        )
        .unwrap();
        assert!(step.endpoint()[0].contains(1.0 + (3.25 * 3.25 - 9.0) / 2.0));
    });
}

#[test]
fn g0_zero_flow_handles_subnormal_time_and_preserves_signed_initial_zero() {
    with_cx(Budget::INFINITE, |cx, _| {
        let step = validated_step(
            cx,
            0.0,
            f64::from_bits(1),
            &[Interval::point(-0.0)],
            &mut |_, _, out, _| {
                out[0] = Interval::point(0.0);
                Ok(())
            },
            &ValidatedOptions::default(),
        )
        .unwrap();
        assert_eq!(step.initial()[0].lo().to_bits(), (-0.0_f64).to_bits());
        assert!(step.endpoint()[0].contains(0.0));
        assert!(step.endpoint()[0].width() < 1.0e-300);
    });
}

#[test]
fn g0_vector_larger_than_one_tile_keeps_every_flow_component() {
    with_cx(Budget::INFINITE, |cx, _| {
        let initial = [Interval::point(0.0); 65];
        let step = validated_step(
            cx,
            0.0,
            0.125,
            &initial,
            &mut |_, _, out, _| {
                for (i, derivative) in out.iter_mut().enumerate() {
                    *derivative = Interval::point(i as f64 - 32.0);
                }
                Ok(())
            },
            &ValidatedOptions::default(),
        )
        .unwrap();
        for (i, value) in step.endpoint().iter().enumerate() {
            assert!(value.contains((i as f64 - 32.0) * 0.125));
        }
        let middle = step.range(Interval::point(0.0625), cx).unwrap();
        for (i, value) in middle.iter().enumerate() {
            assert!(value.contains((i as f64 - 32.0) * 0.0625));
        }
    });
}

#[test]
fn g0_unproved_inclusion_and_rhs_limits_count_every_callback() {
    with_cx(Budget::INFINITE, |cx, _| {
        let mut calls = 0;
        let mut rhs = |_: Interval, x: &[Interval], out: &mut [Interval], _: &Cx<'_>| {
            calls += 1;
            out[0] = x[0];
            Ok(())
        };
        let options = ValidatedOptions {
            max_iterations: 3,
            ..ValidatedOptions::default()
        };
        let result = validated_step(cx, 0.0, 2.0, &[Interval::point(1.0)], &mut rhs, &options);
        assert_eq!(
            result,
            Err(ValidatedError::NoInclusion {
                iterations: 3,
                rhs_evaluations: 3
            })
        );
        assert_eq!(calls, 3);
        let mut calls = 0;
        let options = ValidatedOptions {
            max_rhs_evaluations: 2,
            ..options
        };
        let result = validated_step(
            cx,
            0.0,
            2.0,
            &[Interval::point(1.0)],
            &mut |_, x, out, _| {
                calls += 1;
                out[0] = x[0];
                Ok(())
            },
            &options,
        );
        assert_eq!(result, Err(ValidatedError::RhsEvaluationLimit { limit: 2 }));
        assert_eq!(calls, 2);
    });
}

#[test]
fn g0_width_and_component_work_limits_refuse_before_publication() {
    with_cx(Budget::INFINITE, |cx, _| {
        let options = ValidatedOptions {
            max_width: 0.1,
            ..ValidatedOptions::default()
        };
        let result = validated_step(
            cx,
            0.0,
            1.0,
            &[Interval::point(0.0)],
            &mut |_, _, out, _| {
                out[0] = Interval::point(1.0);
                Ok(())
            },
            &options,
        );
        assert!(matches!(
            result,
            Err(ValidatedError::WidthLimit { component: 0, .. })
        ));
        let mut calls = 0;
        let options = ValidatedOptions {
            max_work: 6,
            ..ValidatedOptions::default()
        };
        let result = validated_step(
            cx,
            0.0,
            1.0,
            &[Interval::point(0.0)],
            &mut |_, _, out, _| {
                calls += 1;
                out[0] = Interval::point(1.0);
                Ok(())
            },
            &options,
        );
        assert_eq!(result, Err(ValidatedError::WorkLimit { limit: 6 }));
        assert_eq!(calls, 1);
    });
}

#[test]
fn g0_unwritten_derivative_cannot_reuse_a_previous_iteration() {
    with_cx(Budget::INFINITE, |cx, _| {
        let mut calls = 0;
        let result = validated_step(
            cx,
            0.0,
            1.0,
            &[Interval::point(0.0); 2],
            &mut |_, _, out, _| {
                calls += 1;
                out[0] = Interval::point(1.0);
                if calls == 1 {
                    out[1] = Interval::point(1.0);
                }
                Ok(())
            },
            &ValidatedOptions::default(),
        );
        assert_eq!(calls, 2);
        assert_eq!(
            result,
            Err(ValidatedError::NonFiniteDerivative { component: 1 })
        );
    });
}

#[test]
fn g0_invalid_inputs_and_nonfinite_arithmetic_refuse() {
    with_cx(Budget::INFINITE, |cx, _| {
        let options = ValidatedOptions::default();
        assert_eq!(
            validate_state(cx, &[], &options),
            Err(ValidatedError::EmptyState)
        );
        assert_eq!(
            validate_state(cx, &[Interval::WHOLE], &options),
            Err(ValidatedError::NonFiniteInitial { component: 0 })
        );
        let tiny = ValidatedOptions {
            max_dimension: 1,
            ..options.clone()
        };
        assert_eq!(
            validate_state(cx, &[Interval::point(0.0); 2], &tiny),
            Err(ValidatedError::DimensionLimit {
                dimension: 2,
                limit: 1
            })
        );
        let bad = ValidatedOptions {
            inflation: f64::NAN,
            ..options.clone()
        };
        assert!(matches!(
            validate_state(cx, &[Interval::point(0.0)], &bad),
            Err(ValidatedError::InvalidOptions(_))
        ));
        let mut calls = 0;
        for (start, end) in [
            (1.0, 1.0),
            (1.0, 0.0),
            (f64::NAN, 1.0),
            (0.0, f64::INFINITY),
            (-f64::MAX, f64::MAX),
        ] {
            let result = validated_step(
                cx,
                start,
                end,
                &[Interval::point(0.0)],
                &mut |_, _, _, _| {
                    calls += 1;
                    Ok(())
                },
                &options,
            );
            assert_eq!(result, Err(ValidatedError::InvalidTime));
        }
        assert_eq!(calls, 0);
        let result = validated_step(
            cx,
            0.0,
            1.0,
            &[Interval::point(f64::MAX)],
            &mut |_, _, out, _| {
                out[0] = Interval::point(0.0);
                Ok(())
            },
            &options,
        );
        assert!(matches!(result, Err(ValidatedError::Arithmetic { .. })));
        let result = validated_step(
            cx,
            0.0,
            1.0,
            &[Interval::point(0.0)],
            &mut |_, _, _, _| Err(ValidatedError::RhsRefused("outside material domain")),
            &options,
        );
        assert_eq!(
            result,
            Err(ValidatedError::RhsRefused("outside material domain"))
        );
    });
}

#[test]
fn g4_cancellation_before_and_after_rhs_returns_no_step() {
    with_cx(Budget::INFINITE, |cx, gate| {
        gate.request();
        let result = validated_step(
            cx,
            0.0,
            1.0,
            &[Interval::point(0.0)],
            &mut |_, _, _, _| panic!("cancelled admission called the RHS"),
            &ValidatedOptions::default(),
        );
        assert_eq!(result, Err(ValidatedError::Cancelled));
    });
    with_cx(Budget::INFINITE, |cx, gate| {
        let mut calls = 0;
        let result = validated_step(
            cx,
            0.0,
            1.0,
            &[Interval::point(0.0)],
            &mut |_, _, out, _| {
                calls += 1;
                out[0] = Interval::point(1.0);
                gate.request();
                Ok(())
            },
            &ValidatedOptions::default(),
        );
        assert_eq!(calls, 1);
        assert_eq!(result, Err(ValidatedError::Cancelled));
    });
}

#[test]
fn g4_ambient_cost_poll_and_unenforceable_deadline_are_not_ignored() {
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
            let result = validated_step(
                cx,
                0.0,
                1.0,
                &[Interval::point(0.0)],
                &mut |_, _, _, _| panic!("refused budget called the RHS"),
                &ValidatedOptions::default(),
            );
            assert!(matches!(
                result,
                Err(ValidatedError::Budget(
                    BudgetRefusal::CostPlanExceedsQuota { .. }
                        | BudgetRefusal::PollsExhausted { .. }
                        | BudgetRefusal::DeadlineWithoutClock { .. }
                ))
            ));
        });
    }
}

#[test]
fn g5_repeated_steps_match_and_range_refuses_extrapolation_or_cancellation() {
    with_cx(Budget::INFINITE, |cx, gate| {
        let mut rhs = |_: Interval, x: &[Interval], out: &mut [Interval], _: &Cx<'_>| {
            out[0] = -x[0];
            Ok(())
        };
        let initial = [Interval::new(1.0, 1.25)];
        let options = ValidatedOptions::default();
        let first = validated_step(cx, 0.0, 0.125, &initial, &mut rhs, &options).unwrap();
        let second = validated_step(cx, 0.0, 0.125, &initial, &mut rhs, &options).unwrap();
        assert_eq!(first, second);
        for invalid in [
            Interval::WHOLE,
            Interval::new(-0.1, 0.1),
            Interval::point(0.126),
        ] {
            assert_eq!(first.range(invalid, cx), Err(ValidatedError::InvalidRange));
        }
        gate.request();
        assert_eq!(
            first.range(Interval::point(0.0625), cx),
            Err(ValidatedError::Cancelled)
        );
    });
}
