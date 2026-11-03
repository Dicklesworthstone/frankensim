//! G0/G1/G3/G4: accuracy and impulse of the optional reversible composition.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, Time, TimeSource};
use fs_vpm::conservative::{
    self, CompositionOrder,
    ambient::{self, AffineFlow2, AmbientVpmRun},
};
use fs_vpm::{VortexParticle as Particle, VpmBudget, VpmError};
use std::sync::atomic::{AtomicU64, Ordering};

const CORE: f64 = 0.12;
const CLOUD: [Particle; 5] = [
    Particle {
        pos: [-0.7, -0.1],
        circulation: 1.0,
    },
    Particle {
        pos: [0.5, -0.3],
        circulation: -0.5,
    },
    Particle {
        pos: [0.2, 0.8],
        circulation: 0.25,
    },
    Particle {
        pos: [-0.4, 0.55],
        circulation: -0.75,
    },
    Particle {
        pos: [0.1, -0.65],
        circulation: 0.0,
    },
];
const ORIGIN: [f64; 2] = [0.2, -0.15];
const VELOCITY: [f64; 2] = [0.12, -0.08];
const GRADIENT: [f64; 3] = [0.25, 0.8, -0.35];

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
                seed: 47,
                kernel_id: 11,
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

fn caps(n: usize, steps: usize) -> VpmBudget {
    VpmBudget::new(n, steps, u64::MAX, n * std::mem::size_of::<Particle>())
}

fn evolve(
    p: &[Particle],
    dt: f64,
    steps: usize,
    flow: AffineFlow2,
    order: CompositionOrder,
) -> AmbientVpmRun {
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        ambient::simulate_with_order_cx(cx, p, dt, steps, CORE, flow, order, caps(p.len(), steps))
    })
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn error(a: &[Particle], b: &[Particle]) -> f64 {
    a.iter()
        .zip(b)
        .flat_map(|(a, b)| (0..2).map(move |k| (a.pos[k] - b.pos[k]).abs()))
        .fold(0.0, f64::max)
}

fn moments(p: &[Particle]) -> [f64; 4] {
    let mut out = [0.0; 4];
    for p in p {
        for (m, x) in out.iter_mut().zip([
            1.0,
            p.pos[0],
            p.pos[1],
            p.pos[0] * p.pos[0] + p.pos[1] * p.pos[1],
        ]) {
            *m += p.circulation * x;
        }
    }
    out
}

// Independent full ODE: direct Biot-Savart plus the authored affine velocity.
// No production pair map, splitting schedule or source exponential is called.
fn rhs(p: &[Particle; 5], forced: bool) -> [[f64; 2]; 5] {
    let mut out = [[0.0; 2]; 5];
    for i in 0..5 {
        if forced {
            let [x, y] = [p[i].pos[0] - ORIGIN[0], p[i].pos[1] - ORIGIN[1]];
            let [a, b, c] = GRADIENT;
            out[i] = [VELOCITY[0] + a * x + b * y, VELOCITY[1] + c * x - a * y];
        }
        for j in 0..5 {
            if i == j {
                continue;
            }
            let x = p[i].pos[0] - p[j].pos[0];
            let y = p[i].pos[1] - p[j].pos[1];
            let factor = p[j].circulation / (std::f64::consts::TAU * (x * x + y * y + CORE * CORE));
            out[i][0] -= factor * y;
            out[i][1] += factor * x;
        }
    }
    out
}

fn reference(time: f64, steps: usize, forced: bool) -> [Particle; 5] {
    let dt = time / steps as f64;
    let mut p = CLOUD;
    for _ in 0..steps {
        let shift = |k: [[f64; 2]; 5], fraction: f64| {
            std::array::from_fn(|i| {
                Particle::new(
                    std::array::from_fn(|axis| p[i].pos[axis] + fraction * dt * k[i][axis]),
                    p[i].circulation,
                )
            })
        };
        let k1 = rhs(&p, forced);
        let k2 = rhs(&shift(k1, 0.5), forced);
        let k3 = rhs(&shift(k2, 0.5), forced);
        let k4 = rhs(&shift(k3, 1.0), forced);
        for i in 0..5 {
            for a in 0..2 {
                p[i].pos[a] += dt * (k1[i][a] + 2.0 * k2[i][a] + 2.0 * k3[i][a] + k4[i][a]) / 6.0;
            }
        }
    }
    p
}

#[test]
fn g1_many_vortex_and_noncommuting_affine_evolution_converge_at_fourth_order() {
    let time = 1.2;
    for forced in [false, true] {
        let flow = if forced {
            source()
        } else {
            AffineFlow2::default()
        };
        let truth = reference(time, 8192, forced);
        assert!(error(&truth, &reference(time, 4096, forced)) < 5e-14);
        let mut errors = [0.0; 3];
        for (i, steps) in [8, 16, 32].into_iter().enumerate() {
            let result = evolve(
                &CLOUD,
                time / steps as f64,
                steps,
                flow,
                CompositionOrder::Fourth,
            );
            errors[i] = error(&result.run.particles, &truth);
        }
        for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
            assert!(
                (12.0..20.0).contains(&ratio),
                "forced={forced}, errors={errors:?}, ratio={ratio}"
            );
        }
        // Equal pair work: 32 fourth-order steps versus 96 original sweeps.
        let second = evolve(&CLOUD, time / 96.0, 96, flow, CompositionOrder::Second);
        let fourth = evolve(&CLOUD, time / 32.0, 32, flow, CompositionOrder::Fourth);
        assert_eq!(second.run.pair_evaluations, fourth.run.pair_evaluations);
        assert!(
            errors[2] < error(&second.run.particles, &truth) / 8.0,
            "equal-work second error={}, fourth error={}",
            error(&second.run.particles, &truth),
            errors[2]
        );
        assert!(fourth.run.work_units <= second.run.work_units);
    }
}

#[test]
fn g0_signed_zero_net_cloud_retains_impulse_and_reversibility_without_projection() {
    let initial = moments(&CLOUD);
    assert_eq!(initial[0], 0.0);
    let result = evolve(
        &CLOUD,
        0.001,
        3000,
        AffineFlow2::default(),
        CompositionOrder::Fourth,
    );
    for (actual, expected) in moments(&result.run.particles).into_iter().zip(initial) {
        near(actual, expected, 3e-12);
    }
    for (actual, expected) in result.run.particles.iter().zip(CLOUD) {
        assert_eq!(actual.circulation.to_bits(), expected.circulation.to_bits());
    }
    for flow in [AffineFlow2::default(), source()] {
        let forward = evolve(&CLOUD, 0.01, 100, flow, CompositionOrder::Fourth);
        let back = evolve(
            &forward.run.particles,
            -0.01,
            100,
            flow,
            CompositionOrder::Fourth,
        );
        assert!(error(&back.run.particles, &CLOUD) < 2e-12);
        assert_eq!(
            forward,
            evolve(&CLOUD, 0.01, 100, flow, CompositionOrder::Fourth)
        );
    }
}

#[test]
fn g0_source_account_closes_signed_stage_moments_and_analytic_first_moment() {
    let rate = 0.3_f64;
    let velocity = [0.12, -0.07];
    let flow = AffineFlow2::new(ORIGIN, velocity, [rate, 0.0, 0.0]).unwrap();
    for extra in [0.0, 0.25] {
        let mut input = CLOUD;
        input[3].circulation += extra;
        let initial = moments(&input);
        let time = 0.9;
        let result = evolve(&input, 0.03, 30, flow, CompositionOrder::Fourth);
        let final_moments = moments(&result.run.particles);
        let source = result.ambient_moment_change;
        for a in 0..2 {
            near(
                source.first_moment_m3_s[a],
                final_moments[1 + a] - initial[1 + a],
                3e-14,
            );
            let k = if a == 0 { rate } else { -rate };
            let c = initial[0] * (velocity[a] - k * ORIGIN[a]);
            let analytic = initial[1 + a] * (k * time).exp() + c * (k * time).exp_m1() / k;
            near(final_moments[1 + a], analytic, 3e-14);
        }
        near(
            source.quadratic_moment_m4_s,
            final_moments[3] - initial[3],
            3e-14,
        );
        assert!(source.quadratic_moment_m4_s.abs() > 1e-3);
        assert_eq!(result.source_particle_updates, 4 * input.len() as u64 * 30);
    }
}

#[test]
fn g3_frames_and_physical_units_transform_the_actual_fourth_order_motion() {
    let angle = 0.63_f64;
    let (s, c) = angle.sin_cos();
    let rotation = [[c, -s], [s, c]];
    let rotate = |x: [f64; 2]| [c * x[0] - s * x[1], s * x[0] + c * x[1]];
    let origin_shift = [-3.0, 7.0];
    let length = 4.0;
    let duration = 2.3;
    let place = |p: [f64; 2]| {
        let q = rotate(p);
        [
            length * q[0] + origin_shift[0],
            length * q[1] + origin_shift[1],
        ]
    };
    let matrix = [[GRADIENT[0], GRADIENT[1]], [GRADIENT[2], -GRADIENT[0]]];
    let rotated_matrix: [[f64; 2]; 2] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let mut sum = 0.0;
            for k in 0..2 {
                for l in 0..2 {
                    sum += rotation[i][k] * matrix[k][l] * rotation[j][l];
                }
            }
            sum / duration
        })
    });
    let transformed_flow = AffineFlow2::new(
        place(ORIGIN),
        rotate(VELOCITY).map(|v| length * v / duration),
        [
            rotated_matrix[0][0],
            rotated_matrix[0][1],
            rotated_matrix[1][0],
        ],
    )
    .unwrap();
    let transformed =
        CLOUD.map(|p| Particle::new(place(p.pos), p.circulation * length * length / duration));
    let original = evolve(&CLOUD, 0.02, 40, source(), CompositionOrder::Fourth);
    let actual = with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        ambient::simulate_with_order_cx(
            cx,
            &transformed,
            0.02 * duration,
            40,
            CORE * length,
            transformed_flow,
            CompositionOrder::Fourth,
            caps(5, 40),
        )
    })
    .unwrap();
    for (got, expected) in actual.run.particles.iter().zip(&original.run.particles) {
        for (x, y) in got.pos.into_iter().zip(place(expected.pos)) {
            near(x, y, 5e-13);
        }
    }
    assert_eq!(actual.run.pair_evaluations, original.run.pair_evaluations);
    assert_eq!(actual.run.work_units, original.run.work_units);
}

#[test]
fn g0_selected_order_and_inactive_paths_preserve_original_receipts() {
    let gate = CancelGate::new_clock_free();
    for (p, dt, steps) in [
        (&CLOUD[..], 0.0, 7),
        (&CLOUD[..], 0.1, 0),
        (&CLOUD[..0], 0.1, 7),
        (&CLOUD[..1], 0.1, 7),
    ] {
        let original = with_cx(&gate, Budget::INFINITE, |cx| {
            conservative::simulate_with_cx(cx, p, dt, steps, CORE, caps(p.len(), steps))
        })
        .unwrap();
        for order in [CompositionOrder::Second, CompositionOrder::Fourth] {
            let run = evolve(p, dt, steps, AffineFlow2::default(), order);
            assert_eq!(run.run, original);
            assert_eq!(run.source_particle_updates, 0);
            assert_eq!(run.ambient_moment_change, Default::default());
        }
    }
    for flow in [AffineFlow2::default(), source()] {
        let original = with_cx(&gate, Budget::INFINITE, |cx| {
            ambient::simulate_with_cx(cx, &CLOUD, 0.03, 5, CORE, flow, caps(5, 5))
        })
        .unwrap();
        assert_eq!(
            original,
            evolve(&CLOUD, 0.03, 5, flow, CompositionOrder::Second)
        );
    }
    let fourth = with_cx(&gate, Budget::INFINITE, |cx| {
        conservative::simulate_with_order_cx(
            cx,
            &CLOUD,
            0.03,
            5,
            CORE,
            CompositionOrder::Fourth,
            caps(5, 5),
        )
    })
    .unwrap();
    assert_eq!(
        fourth,
        evolve(
            &CLOUD,
            0.03,
            5,
            AffineFlow2::default(),
            CompositionOrder::Fourth
        )
        .run
    );
}

#[test]
fn g0_subnormal_time_keeps_resolvable_pair_and_source_motion() {
    let dt = f64::from_bits(1);
    let speed = 1e308;
    let tracer = [Particle::new([0.0; 2], 0.0)];
    let flow = AffineFlow2::new([0.0; 2], [speed, 0.0], [0.0; 3]).unwrap();
    let moved = evolve(&tracer, dt, 1, flow, CompositionOrder::Fourth);
    near(moved.run.particles[0].pos[0], dt * speed, 4e-31);
    let pair = [
        Particle::new([-0.5, 0.0], speed),
        Particle::new([0.5, 0.0], -speed),
    ];
    let result = evolve(
        &pair,
        dt,
        1,
        AffineFlow2::default(),
        CompositionOrder::Fourth,
    );
    let expected = dt * speed / (std::f64::consts::TAU * (1.0 + CORE * CORE));
    for p in &result.run.particles {
        near(p.pos[1], expected, 2e-31);
    }
    // The small signed bridge coefficient cannot be represented here. Its
    // action on a large coordinate would matter, so refuse before evolution.
    let large = [Particle::new([1e308, 0.0], 0.0)];
    let shear = AffineFlow2::new([0.0; 2], [0.0; 2], [0.0, 0.0, 2.0]).unwrap();
    assert!(matches!(
        with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
            ambient::simulate_with_order_cx(
                cx,
                &large,
                dt,
                1,
                CORE,
                shear,
                CompositionOrder::Fourth,
                caps(1, 1),
            )
        }),
        Err(VpmError::NonFiniteResult {
            phase: "fs-vpm.ambient-map",
            ..
        })
    ));
}

#[test]
fn g4_whole_signed_plan_memory_and_map_failures_are_atomic() {
    let gate = CancelGate::new_clock_free();
    let n = CLOUD.len();
    let steps = 7;
    let pairs = 3 * steps as u64 * n as u64 * (n - 1) as u64;
    let sources = 4 * steps as u64 * n as u64;
    let work = pairs + sources + n as u64 + steps as u64;
    let mut local = caps(n, steps);
    local.max_pair_evaluations = pairs;
    let run = with_cx(&gate, Budget::INFINITE.with_cost_quota(work), |cx| {
        ambient::simulate_with_order_cx(
            cx,
            &CLOUD,
            0.01,
            steps,
            CORE,
            source(),
            CompositionOrder::Fourth,
            local,
        )
    })
    .unwrap();
    assert_eq!(run.run.pair_evaluations, pairs);
    assert_eq!(run.source_particle_updates, sources);
    assert_eq!(run.run.work_units, work);
    assert_eq!(run.run.budget.cost_charged, work);
    assert_eq!(run.run.peak_live_bytes, n * std::mem::size_of::<Particle>());
    for (pair_cap, bytes, quota) in [
        (pairs - 1, local.max_live_bytes, work),
        (pairs, local.max_live_bytes - 1, work),
        (pairs, local.max_live_bytes, work - 1),
    ] {
        let limited = VpmBudget::new(n, steps, pair_cap, bytes);
        assert!(
            with_cx(&gate, Budget::INFINITE.with_cost_quota(quota), |cx| {
                ambient::simulate_with_order_cx(
                    cx,
                    &CLOUD,
                    0.01,
                    steps,
                    CORE,
                    source(),
                    CompositionOrder::Fourth,
                    limited,
                )
            })
            .is_err()
        );
    }
    // Finite equilibrium and representable old half-map do not excuse a
    // nonrepresentable fourth-order internal source map.
    let equilibrium = [Particle::new([0.0; 2], 0.0)];
    let extreme = AffineFlow2::new([0.0; 2], [0.0; 2], [1200.0, 0.0, 0.0]).unwrap();
    assert_eq!(
        evolve(&equilibrium, 1.0, 1, extreme, CompositionOrder::Second)
            .run
            .particles,
        equilibrium
    );
    assert!(matches!(
        with_cx(&gate, Budget::INFINITE, |cx| {
            ambient::simulate_with_order_cx(
                cx,
                &equilibrium,
                1.0,
                1,
                CORE,
                extreme,
                CompositionOrder::Fourth,
                caps(1, 1),
            )
        }),
        Err(VpmError::NonFiniteResult { .. })
    ));
    assert_eq!(equilibrium[0].pos, [0.0; 2]);
    let retry = with_cx(&gate, Budget::INFINITE.with_cost_quota(work), |cx| {
        ambient::simulate_with_order_cx(
            cx,
            &CLOUD,
            0.01,
            steps,
            CORE,
            source(),
            CompositionOrder::Fourth,
            local,
        )
    })
    .unwrap();
    assert_eq!(run, retry);
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
fn g4_cancellation_inside_signed_pairs_and_coalesced_sources_publishes_nothing() {
    let input: Vec<_> = (0..39)
        .map(|i| {
            let angle = i as f64 * std::f64::consts::TAU / 39.0;
            let (s, c) = angle.sin_cos();
            Particle::new([2.0 * c, 2.0 * s], if i % 2 == 0 { 0.01 } else { -0.01 })
        })
        .collect();
    let saved = input.clone();
    let budget = Budget::INFINITE.with_deadline(Time::from_nanos(1_000_000));
    let run = |gate: &CancelGate, clock: &TripClock<'_>| {
        with_cx(gate, budget, |cx| {
            let timed = Cx::new(gate, cx.arena(), cx.stream_key(), budget, cx.mode())
                .with_time_source(clock);
            ambient::simulate_with_order_cx(
                &timed,
                &input,
                0.01,
                1,
                CORE,
                source(),
                CompositionOrder::Fourth,
                caps(input.len(), 1),
            )
        })
    };
    let active = CancelGate::new_clock_free();
    let counter = TripClock {
        gate: &active,
        calls: AtomicU64::new(0),
        cancel_on: u64::MAX,
    };
    let reference = run(&active, &counter).unwrap();
    let polls = counter.calls.load(Ordering::Relaxed);
    let (mut source_seen, mut pair_seen) = (false, false);
    for cancel_on in 1..polls {
        let gate = CancelGate::new_clock_free();
        let clock = TripClock {
            gate: &gate,
            calls: AtomicU64::new(0),
            cancel_on,
        };
        match run(&gate, &clock) {
            Err(VpmError::ExecutionBudget(BudgetRefusal::Cancelled { phase })) => {
                source_seen |= phase == "fs-vpm.ambient-source";
                pair_seen |= phase == "fs-vpm.conservative-reverse";
            }
            other => panic!("cancel_on={cancel_on}, result={other:?}"),
        }
        assert_eq!(input, saved);
    }
    assert!(source_seen && pair_seen);
    assert_eq!(run(&active, &counter).unwrap(), reference);
}
