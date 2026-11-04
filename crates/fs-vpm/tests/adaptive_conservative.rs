//! Independent full-ODE and exact-map checks for local adaptive vortex time.
//! Local acceptance is not used as an assertion of global/continuum error.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, Time, TimeSource};
use fs_vpm::conservative::{
    CompositionOrder,
    adaptive::{self, AdaptiveBudget, AdaptiveError, AdaptiveOptions, AmbientAdaptiveBudget},
    ambient::{self, AffineFlow2},
};
use fs_vpm::{VortexParticle as Particle, VpmBudget, VpmError};
use std::mem::size_of;
use std::sync::atomic::{AtomicU64, Ordering};

const CORE: f64 = 0.03;
const CLOUD: [Particle; 5] = [
    Particle {
        pos: [-0.09, 0.015],
        circulation: 1.0,
    },
    Particle {
        pos: [0.085, -0.02],
        circulation: -0.5,
    },
    Particle {
        pos: [-0.03, 0.22],
        circulation: 0.25,
    },
    Particle {
        pos: [-0.26, -0.1],
        circulation: -0.75,
    },
    Particle {
        pos: [0.2, -0.15],
        circulation: 0.0,
    },
];
const ORIGIN: [f64; 2] = [0.12, -0.07];
const VELOCITY: [f64; 2] = [0.2, -0.1];
const GRADIENT: [f64; 3] = [0.7, 0.4, -0.15];

fn source() -> AffineFlow2 {
    AffineFlow2::new(ORIGIN, VELOCITY, GRADIENT).unwrap()
}

fn with_cx<T>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = ArenaPool::new(ArenaConfig::default());
    let result = pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x4144_4150_5456,
                kernel_id: 0x5650_4d,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        ))
    });
    assert!(pool.stats().quiescent());
    result
}

fn options(tolerance: f64) -> AdaptiveOptions {
    AdaptiveOptions {
        initial_dt_s: 0.2,
        min_dt_s: 1e-8,
        max_dt_s: 0.2,
        absolute_position_tolerance_m: tolerance,
        relative_tolerance: 0.0,
    }
}

fn caps(n: usize) -> AmbientAdaptiveBudget {
    AmbientAdaptiveBudget {
        run: AdaptiveBudget {
            max_particles: n,
            max_attempts: 10_000,
            max_pair_evaluations: 10_000_000,
            max_work_units: 20_000_000,
            max_live_bytes: 3 * n * size_of::<Particle>(),
        },
        max_source_particle_updates: 1_000_000,
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance={tolerance}"
    );
}

fn error(a: &[Particle], b: &[Particle]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(a, b)| (a.pos[0] - b.pos[0]).hypot(a.pos[1] - b.pos[1]))
        .fold(0.0, f64::max)
}

fn bits(a: &[Particle], b: &[Particle]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_eq!(a.pos.map(f64::to_bits), b.pos.map(f64::to_bits));
        assert_eq!(a.circulation.to_bits(), b.circulation.to_bits());
    }
}

fn moments(particles: &[Particle]) -> [f64; 4] {
    let mut result = [0.0; 4];
    for p in particles {
        for (m, x) in
            result
                .iter_mut()
                .zip([1.0, p.pos[0], p.pos[1], p.pos[0].powi(2) + p.pos[1].powi(2)])
        {
            *m += p.circulation * x;
        }
    }
    result
}

// Independent direct Biot-Savart ODE with an authored affine source. No
// production pair solver, splitting coefficient or source map is called.
fn rhs(p: &[Particle; 5], forced: bool) -> [[f64; 2]; 5] {
    let mut v = [[0.0; 2]; 5];
    for i in 0..5 {
        if forced {
            let [x, y] = [p[i].pos[0] - ORIGIN[0], p[i].pos[1] - ORIGIN[1]];
            let [a, b, c] = GRADIENT;
            v[i] = [VELOCITY[0] + a * x + b * y, VELOCITY[1] + c * x - a * y];
        }
        for j in 0..5 {
            if i != j {
                let [x, y] = [p[i].pos[0] - p[j].pos[0], p[i].pos[1] - p[j].pos[1]];
                let factor =
                    p[j].circulation / (std::f64::consts::TAU * (x * x + y * y + CORE * CORE));
                v[i][0] -= factor * y;
                v[i][1] += factor * x;
            }
        }
    }
    v
}

fn reference(time: f64, steps: usize, forced: bool) -> [Particle; 5] {
    let dt = time / steps as f64;
    let mut p = CLOUD;
    for _ in 0..steps {
        let shifted = |k: [[f64; 2]; 5], fraction: f64| {
            std::array::from_fn(|i| {
                Particle::new(
                    std::array::from_fn(|a| p[i].pos[a] + fraction * dt * k[i][a]),
                    p[i].circulation,
                )
            })
        };
        let k1 = rhs(&p, forced);
        let k2 = rhs(&shifted(k1, 0.5), forced);
        let k3 = rhs(&shifted(k2, 0.5), forced);
        let k4 = rhs(&shifted(k3, 1.0), forced);
        for i in 0..5 {
            for a in 0..2 {
                p[i].pos[a] += dt * (k1[i][a] + 2.0 * k2[i][a] + 2.0 * k3[i][a] + k4[i][a]) / 6.0;
            }
        }
    }
    p
}

#[test]
fn g1_tolerance_refinement_resolves_close_vortices_and_noncommuting_strain() {
    let time = 0.3;
    let initial = moments(&CLOUD);
    assert_eq!(initial[0], 0.0);
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for forced in [false, true] {
            let truth = reference(time, 8192, forced);
            assert!(error(&truth, &reference(time, 4096, forced)) < 2e-12);
            let flow = if forced {
                source()
            } else {
                AffineFlow2::default()
            };
            let mut tight_work = [0; 2];
            let mut tight_error = [0.0; 2];
            for (i, order) in [CompositionOrder::Second, CompositionOrder::Fourth]
                .into_iter()
                .enumerate()
            {
                let run = |tolerance| {
                    adaptive::simulate_with_affine_flow_cx(
                        cx,
                        &CLOUD,
                        time,
                        CORE,
                        order,
                        flow,
                        options(tolerance),
                        caps(5),
                    )
                    .unwrap()
                };
                let loose = run(1e-4);
                let tight = run(1e-8);
                let loose_error = error(&loose.run.particles, &truth);
                let fine_error = error(&tight.run.particles, &truth);
                assert!(
                    fine_error < 2e-5,
                    "forced={forced}, order={order:?}, error={fine_error}"
                );
                assert!(
                    loose_error > 20.0 * fine_error,
                    "forced={forced}, order={order:?}, errors={loose_error} -> {fine_error}"
                );
                assert!(tight.run.accepted_steps > loose.run.accepted_steps);
                assert!(tight.run.rejected_steps > 0);
                for result in [&loose, &tight] {
                    let run = &result.run;
                    assert_eq!(run.elapsed_s.to_bits(), time.to_bits());
                    assert_eq!(run.attempts, run.accepted_steps + run.rejected_steps);
                    assert!((0.0..=1.0).contains(&run.max_accepted_error_ratio));
                    let q = if order == CompositionOrder::Second {
                        1
                    } else {
                        3
                    };
                    let sources = if !forced {
                        0
                    } else if q == 1 {
                        25
                    } else {
                        55
                    };
                    assert_eq!(run.pair_evaluations, 60 * q * run.attempts as u64);
                    assert_eq!(run.fraction_evaluations, 0);
                    assert_eq!(
                        result.source_particle_updates,
                        sources * run.attempts as u64
                    );
                    assert_eq!(
                        run.work_units,
                        10 + (60 * q + sources + 19) * run.attempts as u64
                    );
                    assert_eq!(run.work_units, run.budget.cost_charged);
                    assert_eq!(run.peak_live_bytes, 15 * size_of::<Particle>());
                    for (p, original) in run.particles.iter().zip(CLOUD) {
                        assert_eq!(p.circulation.to_bits(), original.circulation.to_bits());
                    }
                    let final_moments = moments(&run.particles);
                    let delta = result.ambient_moment_change;
                    near(
                        delta.first_moment_m3_s[0],
                        final_moments[1] - initial[1],
                        2e-13,
                    );
                    near(
                        delta.first_moment_m3_s[1],
                        final_moments[2] - initial[2],
                        2e-13,
                    );
                    near(
                        delta.quadratic_moment_m4_s,
                        final_moments[3] - initial[3],
                        2e-13,
                    );
                    if forced {
                        // Net circulation is zero, so M'=A M independently of
                        // vortex shape, core, source origin and uniform drift.
                        let [a, b, c] = GRADIENT;
                        let lambda = (a * a + b * c).sqrt();
                        let ch = (lambda * time).cosh();
                        let sh = (lambda * time).sinh() / lambda;
                        near(
                            final_moments[1],
                            ch * initial[1] + sh * (a * initial[1] + b * initial[2]),
                            2e-13,
                        );
                        near(
                            final_moments[2],
                            ch * initial[2] + sh * (c * initial[1] - a * initial[2]),
                            2e-13,
                        );
                    } else {
                        assert_eq!(delta, Default::default());
                        for (a, b) in final_moments.into_iter().zip(initial) {
                            near(a, b, 2e-13);
                        }
                    }
                }
                tight_work[i] = tight.run.work_units;
                tight_error[i] = fine_error;
                assert_eq!(tight, run(1e-8));
            }
            // On this independently resolved encounter, the higher-order map
            // reaches the same local target with less actual admitted work.
            assert!(
                tight_work[1] < tight_work[0],
                "forced={forced}, work={tight_work:?}"
            );
            assert!(tight_error[1] < 2e-5);
        }
    });
}

#[test]
fn g0_accepted_fine_state_is_original_map_without_extrapolation() {
    let time = 0.125;
    let settings = AdaptiveOptions {
        initial_dt_s: time,
        min_dt_s: time,
        max_dt_s: time,
        ..options(1.0)
    };
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for order in [CompositionOrder::Second, CompositionOrder::Fourth] {
            let unforced =
                adaptive::simulate_with_cx(cx, &CLOUD, time, CORE, order, settings, caps(5).run)
                    .unwrap();
            for flow in [AffineFlow2::default(), source()] {
                let adaptive = adaptive::simulate_with_affine_flow_cx(
                    cx,
                    &CLOUD,
                    time,
                    CORE,
                    order,
                    flow,
                    settings,
                    caps(5),
                )
                .unwrap();
                let original = ambient::simulate_with_order_cx(
                    cx,
                    &CLOUD,
                    time / 2.0,
                    2,
                    CORE,
                    flow,
                    order,
                    VpmBudget::new(5, 2, 1_000, 5 * size_of::<Particle>()),
                )
                .unwrap();
                assert_eq!(adaptive.run.accepted_steps, 1);
                assert_eq!(adaptive.run.rejected_steps, 0);
                if flow == AffineFlow2::default() {
                    assert_eq!(adaptive.run, unforced);
                    bits(&adaptive.run.particles, &original.run.particles);
                } else {
                    // Only neighboring exact source maps are combined; their
                    // floating rounding need not match two stored passes.
                    assert!(error(&adaptive.run.particles, &original.run.particles) < 2e-15);
                    near(
                        adaptive.ambient_moment_change.quadratic_moment_m4_s,
                        original.ambient_moment_change.quadratic_moment_m4_s,
                        2e-15,
                    );
                }
            }
        }
    });
}

#[test]
fn g3_adaptive_error_scale_is_objective_under_frame_and_unit_changes() {
    let (s, c) = 0.63_f64.sin_cos();
    let rotation = [[c, -s], [s, c]];
    let rotate = |x: [f64; 2]| [c * x[0] - s * x[1], s * x[0] + c * x[1]];
    let length = 4.0;
    let time_scale = 2.0;
    let place = |p| {
        let p = rotate(p);
        [length * p[0] - 100.0, length * p[1] + 70.0]
    };
    let a = [[GRADIENT[0], GRADIENT[1]], [GRADIENT[2], -GRADIENT[0]]];
    let transformed_a: [[f64; 2]; 2] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let mut sum = 0.0;
            for k in 0..2 {
                for l in 0..2 {
                    sum += rotation[i][k] * a[k][l] * rotation[j][l];
                }
            }
            sum / time_scale
        })
    });
    let transformed_flow = AffineFlow2::new(
        place(ORIGIN),
        rotate(VELOCITY).map(|v| length * v / time_scale),
        [
            transformed_a[0][0],
            transformed_a[0][1],
            transformed_a[1][0],
        ],
    )
    .unwrap();
    let input =
        CLOUD.map(|p| Particle::new(place(p.pos), p.circulation * length * length / time_scale));
    let settings = AdaptiveOptions {
        relative_tolerance: 1e-6,
        ..options(1e-9)
    };
    let transformed_settings = AdaptiveOptions {
        initial_dt_s: settings.initial_dt_s * time_scale,
        min_dt_s: settings.min_dt_s * time_scale,
        max_dt_s: settings.max_dt_s * time_scale,
        absolute_position_tolerance_m: settings.absolute_position_tolerance_m * length,
        ..settings
    };
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for order in [CompositionOrder::Second, CompositionOrder::Fourth] {
            let original = adaptive::simulate_with_affine_flow_cx(
                cx,
                &CLOUD,
                0.3,
                CORE,
                order,
                source(),
                settings,
                caps(5),
            )
            .unwrap();
            let transformed = adaptive::simulate_with_affine_flow_cx(
                cx,
                &input,
                0.3 * time_scale,
                CORE * length,
                order,
                transformed_flow,
                transformed_settings,
                caps(5),
            )
            .unwrap();
            assert_eq!(original.run.attempts, transformed.run.attempts);
            assert_eq!(original.run.work_units, transformed.run.work_units);
            for (p, q) in transformed.run.particles.iter().zip(original.run.particles) {
                for (x, y) in p.pos.into_iter().zip(place(q.pos)) {
                    near(x, y, 2e-11);
                }
            }
        }
    });
}

#[test]
fn g0_inactive_and_subnormal_windows_and_short_endpoint_are_physical() {
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for (n, duration) in [(0, 0.3), (1, 0.3), (5, 0.0)] {
            let mut initial = CLOUD[..n].to_vec();
            if n > 0 {
                initial[0].circulation = -0.0;
            }
            let limits = AdaptiveBudget {
                max_attempts: 0,
                max_pair_evaluations: 0,
                max_work_units: n as u64,
                max_live_bytes: n * size_of::<Particle>(),
                ..caps(n).run
            };
            let run = adaptive::simulate_with_cx(
                cx,
                &initial,
                duration,
                CORE,
                CompositionOrder::Fourth,
                options(1e-8),
                limits,
            )
            .unwrap();
            bits(&run.particles, &initial);
            assert_eq!(run.attempts, 0);
            assert_eq!(run.work_units, n as u64);
            assert_eq!(run.elapsed_s.to_bits(), duration.to_bits());
        }
        let dt = f64::from_bits(1);
        let speed = 1e308;
        let tracer = [Particle::new([0.0; 2], 0.0)];
        let flow = AffineFlow2::new([0.0; 2], [speed, 0.0], [0.0; 3]).unwrap();
        let pair = [
            Particle::new([-0.5, 0.0], speed),
            Particle::new([0.5, 0.0], -speed),
        ];
        let settings = AdaptiveOptions {
            initial_dt_s: dt,
            min_dt_s: dt,
            max_dt_s: dt,
            ..options(1e-27)
        };
        for order in [CompositionOrder::Second, CompositionOrder::Fourth] {
            let moved = adaptive::simulate_with_affine_flow_cx(
                cx,
                &tracer,
                dt,
                CORE,
                order,
                flow,
                settings,
                caps(1),
            )
            .unwrap();
            near(moved.run.particles[0].pos[0], dt * speed, 1e-30);
            assert_eq!(moved.run.accepted_steps, 1);
            let moved_pair =
                adaptive::simulate_with_cx(cx, &pair, dt, CORE, order, settings, caps(2).run)
                    .unwrap();
            let expected = dt * speed / (std::f64::consts::TAU * (1.0 + CORE * CORE));
            for p in moved_pair.particles {
                near(p.pos[1], expected, 1e-30);
            }
        }
        let flow = AffineFlow2::new([0.0; 2], [0.7, -0.2], [0.0; 3]).unwrap();
        let settings = AdaptiveOptions {
            initial_dt_s: 0.3,
            min_dt_s: 0.2,
            max_dt_s: 0.3,
            ..options(1e-14)
        };
        let run = adaptive::simulate_with_affine_flow_cx(
            cx,
            &tracer,
            1.0,
            CORE,
            CompositionOrder::Fourth,
            flow,
            settings,
            caps(1),
        )
        .unwrap();
        assert_eq!(run.run.accepted_steps, 4);
        assert_eq!(run.run.elapsed_s.to_bits(), 1.0_f64.to_bits());
        near(run.run.particles[0].pos[0], 0.7, 1e-15);
        near(run.run.particles[0].pos[1], -0.2, 1e-15);
    });
}

#[test]
fn g4_all_rejected_work_source_visits_and_live_banks_are_capped() {
    let input = CLOUD;
    let saved = input;
    let settings = options(1e-7);
    let gate = CancelGate::new_clock_free();
    let run = |cx: &Cx<'_>, limits| {
        adaptive::simulate_with_affine_flow_cx(
            cx,
            &input,
            0.3,
            CORE,
            CompositionOrder::Fourth,
            source(),
            settings,
            limits,
        )
    };
    let good = with_cx(&gate, Budget::INFINITE, |cx| run(cx, caps(5))).unwrap();
    assert!(good.run.rejected_steps > 0);
    let exact = AmbientAdaptiveBudget {
        run: AdaptiveBudget {
            max_particles: 5,
            max_attempts: good.run.attempts,
            max_pair_evaluations: good.run.pair_evaluations,
            max_work_units: good.run.work_units,
            max_live_bytes: good.run.peak_live_bytes,
        },
        max_source_particle_updates: good.source_particle_updates,
    };
    let replay = with_cx(
        &gate,
        Budget::INFINITE.with_cost_quota(good.run.work_units),
        |cx| run(cx, exact),
    )
    .unwrap();
    // Different declared ceilings change admission metadata, not the path.
    let mut expected = good.clone();
    expected.run.budget = replay.run.budget;
    assert_eq!(expected, replay);
    assert_eq!(replay.run.budget.planned_cost, good.run.work_units);
    assert_eq!(replay.run.budget.cost_charged, good.run.work_units);
    for limited in [
        AmbientAdaptiveBudget {
            run: AdaptiveBudget {
                max_particles: 4,
                ..exact.run
            },
            ..exact
        },
        AmbientAdaptiveBudget {
            run: AdaptiveBudget {
                max_attempts: exact.run.max_attempts - 1,
                ..exact.run
            },
            ..exact
        },
        AmbientAdaptiveBudget {
            run: AdaptiveBudget {
                max_pair_evaluations: exact.run.max_pair_evaluations - 1,
                ..exact.run
            },
            ..exact
        },
        AmbientAdaptiveBudget {
            run: AdaptiveBudget {
                max_work_units: exact.run.max_work_units - 1,
                ..exact.run
            },
            ..exact
        },
        AmbientAdaptiveBudget {
            run: AdaptiveBudget {
                max_live_bytes: exact.run.max_live_bytes - 1,
                ..exact.run
            },
            ..exact
        },
        AmbientAdaptiveBudget {
            max_source_particle_updates: exact.max_source_particle_updates - 1,
            ..exact
        },
    ] {
        assert!(matches!(
            with_cx(&gate, Budget::INFINITE, |cx| run(cx, limited)),
            Err(AdaptiveError::Engine(VpmError::CapExceeded { .. }))
        ));
    }
    assert!(matches!(
        with_cx(
            &gate,
            Budget::INFINITE.with_cost_quota(good.run.work_units - 1),
            |cx| run(cx, exact)
        ),
        Err(AdaptiveError::Engine(VpmError::ExecutionBudget(_)))
    ));
    let retry = with_cx(&gate, Budget::INFINITE, |cx| run(cx, exact)).unwrap();
    expected.run.budget = retry.run.budget;
    assert_eq!(expected, retry);
    bits(&input, &saved);
}

#[test]
fn g4_tolerance_floor_invalid_data_and_late_numeric_failure_publish_nothing() {
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let run = |p: &[Particle], duration, settings| {
            adaptive::simulate_with_cx(
                cx,
                p,
                duration,
                CORE,
                CompositionOrder::Fourth,
                settings,
                caps(p.len()).run,
            )
        };
        let floor = AdaptiveOptions {
            min_dt_s: 0.2,
            ..options(1e-15)
        };
        assert!(matches!(
            run(&CLOUD, 0.3, floor),
            Err(AdaptiveError::ToleranceUnreachable { .. })
        ));
        for settings in [
            AdaptiveOptions {
                initial_dt_s: 0.0,
                ..options(1e-8)
            },
            AdaptiveOptions {
                min_dt_s: 0.3,
                ..options(1e-8)
            },
            AdaptiveOptions {
                max_dt_s: 0.1,
                ..options(1e-8)
            },
            AdaptiveOptions {
                absolute_position_tolerance_m: 0.0,
                ..options(1e-8)
            },
            AdaptiveOptions {
                relative_tolerance: f64::NAN,
                ..options(1e-8)
            },
        ] {
            assert!(matches!(
                run(&CLOUD, 0.0, settings),
                Err(AdaptiveError::Engine(VpmError::InvalidScalar { .. }))
            ));
        }
        for duration in [-1.0, f64::INFINITY, f64::NAN] {
            assert!(matches!(
                run(&CLOUD, duration, options(1e-8)),
                Err(AdaptiveError::Engine(VpmError::InvalidScalar { .. }))
            ));
        }
        let mut invalid = CLOUD;
        invalid[2].circulation = f64::NAN;
        assert!(matches!(
            run(&invalid, 0.0, options(1e-8)),
            Err(AdaptiveError::Engine(VpmError::InvalidScalar { .. }))
        ));
        let huge = [Particle::new([1e307, 0.0], 0.0)];
        let saved = huge;
        let flow = AffineFlow2::new([0.0; 2], [1e308, 0.0], [0.0; 3]).unwrap();
        let settings = AdaptiveOptions {
            initial_dt_s: 1.0,
            min_dt_s: 1.0,
            max_dt_s: 1.0,
            ..options(1e308)
        };
        // Even an equilibrium tracer cannot admit an unrepresentable source
        // exponential. This refuses during preparation, before copy work.
        let extreme = AffineFlow2::new([0.0; 2], [0.0; 2], [10_000.0, 0.0, 0.0]).unwrap();
        assert!(matches!(
            adaptive::simulate_with_affine_flow_cx(
                cx,
                &[Particle::new([0.0; 2], 0.0)],
                1.0,
                CORE,
                CompositionOrder::Fourth,
                extreme,
                settings,
                caps(1),
            ),
            Err(AdaptiveError::Engine(VpmError::NonFiniteResult {
                phase: "fs-vpm.ambient-map",
                ..
            }))
        ));
        assert!(matches!(
            adaptive::simulate_with_affine_flow_cx(
                cx,
                &huge,
                4.0,
                CORE,
                CompositionOrder::Second,
                flow,
                settings,
                caps(1),
            ),
            Err(AdaptiveError::Engine(VpmError::NonFiniteResult { .. }))
        ));
        bits(&huge, &saved);
    });
}

struct TripClock<'a> {
    gate: &'a CancelGate,
    calls: AtomicU64,
    cancel_on: u64,
}
impl TimeSource for TripClock<'_> {
    fn now(&self) -> Time {
        if self.calls.fetch_add(1, Ordering::Relaxed) + 1 == self.cancel_on {
            self.gate.request();
        }
        Time::ZERO
    }
}

#[test]
fn g4_cancellation_inside_trials_sources_and_final_publication_is_atomic() {
    let input: Vec<_> = (0..39)
        .map(|i| {
            let (s, c) = (i as f64 * std::f64::consts::TAU / 39.0).sin_cos();
            Particle::new([2.0 * c, 2.0 * s], if i % 2 == 0 { 0.01 } else { -0.01 })
        })
        .collect();
    let saved = input.clone();
    let budget = Budget::INFINITE.with_deadline(Time::from_nanos(1_000_000));
    let settings = AdaptiveOptions {
        initial_dt_s: 0.01,
        min_dt_s: 0.01,
        max_dt_s: 0.01,
        ..options(1e-3)
    };
    let run = |gate: &CancelGate, clock: &TripClock<'_>| {
        with_cx(gate, budget, |cx| {
            let timed = Cx::new(gate, cx.arena(), cx.stream_key(), budget, cx.mode())
                .with_time_source(clock);
            adaptive::simulate_with_affine_flow_cx(
                &timed,
                &input,
                0.01,
                CORE,
                CompositionOrder::Fourth,
                source(),
                settings,
                caps(39),
            )
        })
    };
    let gate = CancelGate::new_clock_free();
    let counter = TripClock {
        gate: &gate,
        calls: AtomicU64::new(0),
        cancel_on: u64::MAX,
    };
    let good = run(&gate, &counter).unwrap();
    let polls = counter.calls.load(Ordering::Relaxed);
    let (mut pair_seen, mut source_seen, mut publication_seen) = (false, false, false);
    for cancel_on in 1..polls {
        let gate = CancelGate::new_clock_free();
        let clock = TripClock {
            gate: &gate,
            calls: AtomicU64::new(0),
            cancel_on,
        };
        match run(&gate, &clock) {
            Err(AdaptiveError::Engine(VpmError::ExecutionBudget(BudgetRefusal::Cancelled {
                phase,
            }))) => {
                pair_seen |= phase == "fs-vpm.conservative-reverse";
                source_seen |= phase == "fs-vpm.ambient-source";
                publication_seen |= phase == "fs-vpm.adaptive-publication";
            }
            other => panic!("cancel_on={cancel_on}, result={other:?}"),
        }
        bits(&input, &saved);
    }
    assert!(pair_seen && source_seen && publication_seen);
    assert_eq!(good, run(&gate, &counter).unwrap());
}
