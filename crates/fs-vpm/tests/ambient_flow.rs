//! G0/G1/G3/G4: prescribed affine flow on original conservative/PSE owners.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey};
use fs_vpm::conservative::{
    self,
    ambient::{self, AffineFlow2, AmbientVpmRun},
};
use fs_vpm::diffusion::PseDiffusion;
use fs_vpm::{VortexParticle as Particle, VpmBudget, VpmError};

fn with_cx<T>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = ArenaPool::new(ArenaConfig::default());
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 41,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        );
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn caps(n: usize, steps: usize) -> VpmBudget {
    VpmBudget::new(n, steps, u64::MAX, n * std::mem::size_of::<Particle>())
}

fn flow(origin: [f64; 2], velocity: [f64; 2], gradient: [f64; 3]) -> AffineFlow2 {
    AffineFlow2::new(origin, velocity, gradient).unwrap()
}

fn evolve(p: &[Particle], dt: f64, steps: usize, source: AffineFlow2) -> AmbientVpmRun {
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        ambient::simulate_with_cx(cx, p, dt, steps, 0.1, source, caps(p.len(), steps))
    })
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn moments(p: &[Particle]) -> [f64; 4] {
    p.iter().fold([0.0; 4], |mut m, p| {
        m[0] += p.circulation;
        m[1] += p.circulation * p.pos[0];
        m[2] += p.circulation * p.pos[1];
        m[3] += p.circulation * (p.pos[0] * p.pos[0] + p.pos[1] * p.pos[1]);
        m
    })
}

const CLOUD: [Particle; 3] = [
    Particle {
        pos: [-0.7, -0.1],
        circulation: 0.8,
    },
    Particle {
        pos: [0.5, -0.3],
        circulation: -0.3,
    },
    Particle {
        pos: [0.2, 0.8],
        circulation: 0.2,
    },
];

#[test]
fn g1_single_vortex_matches_independent_affine_solutions() {
    let p = Particle::new([1.2, -0.7], 2.5);
    let origin = [0.4, -0.3];
    let velocity = [0.2, -0.6];
    let time = 1.1_f64;
    let relative = [p.pos[0] - origin[0], p.pos[1] - origin[1]];
    let extension = 0.7_f64;
    let shear = 1.7_f64;
    let omega = 0.9_f64;
    let centre = [
        origin[0] - velocity[1] / omega,
        origin[1] + velocity[0] / omega,
    ];
    let (s, c) = (omega * time).sin_cos();
    let cases = [
        (
            [0.0, 0.0, 0.0],
            [p.pos[0] + time * velocity[0], p.pos[1] + time * velocity[1]],
        ),
        (
            [0.0, shear, 0.0],
            [
                p.pos[0]
                    + time * (velocity[0] + shear * relative[1])
                    + 0.5 * shear * velocity[1] * time * time,
                p.pos[1] + time * velocity[1],
            ],
        ),
        (
            [0.0, -omega, omega],
            [
                centre[0] + c * (p.pos[0] - centre[0]) - s * (p.pos[1] - centre[1]),
                centre[1] + s * (p.pos[0] - centre[0]) + c * (p.pos[1] - centre[1]),
            ],
        ),
        (
            [extension, 0.0, 0.0],
            [
                origin[0]
                    + (extension * time).exp() * relative[0]
                    + velocity[0] * (extension * time).exp_m1() / extension,
                origin[1] + (-extension * time).exp() * relative[1]
                    - velocity[1] * (-extension * time).exp_m1() / extension,
            ],
        ),
        (
            [1.0, 1.0, -1.0],
            [
                p.pos[0]
                    + time * (relative[0] + relative[1] + velocity[0])
                    + 0.5 * time * time * (velocity[0] + velocity[1]),
                p.pos[1] + time * (-relative[0] - relative[1] + velocity[1])
                    - 0.5 * time * time * (velocity[0] + velocity[1]),
            ],
        ),
    ];
    for (gradient, expected) in cases {
        let run = evolve(&[p], time, 1, flow(origin, velocity, gradient));
        for axis in 0..2 {
            close(run.run.particles[0].pos[axis], expected[axis], 3e-15);
        }
        assert_eq!(run.run.particles[0].circulation, p.circulation);
        assert_eq!(run.run.pair_evaluations, 0);
        assert_eq!(run.source_particle_updates, 2);
        assert_eq!(run.run.work_units, 4);
        assert_eq!(run.run.peak_live_bytes, std::mem::size_of::<Particle>());
    }
    // A contracting eigenvector must survive even when cosh-sinh would cancel.
    let contract = evolve(
        &[Particle::new([1.0, 1.0], 0.0)],
        1.0,
        1,
        flow([0.0; 2], [0.0; 2], [80.0, 0.0, 0.0]),
    );
    close(
        contract.run.particles[0].pos[0] / 80.0_f64.exp(),
        1.0,
        8e-16,
    );
    close(
        contract.run.particles[0].pos[1] / (-80.0_f64).exp(),
        1.0,
        8e-16,
    );
    let rotated_contract = evolve(
        &[Particle::new([1.0, -1.0], 1.0)],
        1.0,
        1,
        flow([0.0; 2], [0.0; 2], [0.0, 100.0, 100.0]),
    );
    close(
        rotated_contract.run.particles[0].pos[0] / (-100.0_f64).exp(),
        1.0,
        1e-15,
    );
    close(
        rotated_contract.run.particles[0].pos[1] / (-100.0_f64).exp(),
        -1.0,
        1e-15,
    );
    // A forced expanding mode at equilibrium must remain so while the other
    // mode contracts; large exponential source terms must not be subtracted.
    let forced_contract = evolve(
        &[Particle::new([1.0, 1.0], 1.0)],
        1.0,
        1,
        flow([0.0; 2], [-100.0, 0.0], [100.0, 0.0, 0.0]),
    );
    assert_eq!(forced_contract.run.particles[0].pos[0], 1.0);
    close(
        forced_contract.run.particles[0].pos[1] / (-100.0_f64).exp(),
        1.0,
        1e-15,
    );
    let tiny = evolve(
        &[Particle::new([1e12, 0.0], 0.0)],
        1.0,
        1,
        flow([0.0; 2], [0.0; 2], [0.0, -1e-12, 1e-12]),
    );
    close(tiny.run.particles[0].pos[1], 1.0, 2e-16);
    // Exact affine stagnation is an identity even with nonzero source entries.
    let fixed = Particle::new([1.0, -2.0], 1.0);
    assert_eq!(
        evolve(
            &[fixed],
            0.7,
            1,
            flow([0.0; 2], [-5.0, -8.0], [1.0, -2.0, 6.0])
        )
        .run
        .particles,
        [fixed]
    );
}

#[test]
fn g0_zero_source_and_zero_viscosity_preserve_original_bits_and_resources() {
    let law = PseDiffusion::new(0.2, 1.0).unwrap();
    let areas = [0.5, 0.8, 1.1];
    let zero = AffineFlow2::default();
    let source = flow([0.2, -0.1], [0.3, -0.2], [0.2, 0.3, -0.1]);
    with_cx(
        &CancelGate::new_clock_free(),
        Budget::INFINITE.with_poll_quota(100),
        |cx| {
            for (dt, steps) in [(0.03, 4), (0.0, 3), (0.4, 0)] {
                let native =
                    conservative::simulate_with_cx(cx, &CLOUD, dt, steps, 0.1, caps(3, steps))
                        .unwrap();
                let ambient =
                    ambient::simulate_with_cx(cx, &CLOUD, dt, steps, 0.1, zero, caps(3, steps))
                        .unwrap();
                assert_eq!(native, ambient.run);
                assert_eq!(ambient.ambient_moment_change, Default::default());
                assert_eq!(ambient.source_particle_updates, 0);
                let native = conservative::simulate_viscous_with_cx(
                    cx,
                    &CLOUD,
                    &areas,
                    dt,
                    steps,
                    0.1,
                    law,
                    caps(3, steps),
                )
                .unwrap();
                assert_eq!(
                    native,
                    ambient::simulate_viscous_with_cx(
                        cx,
                        &CLOUD,
                        &areas,
                        dt,
                        steps,
                        0.1,
                        zero,
                        law,
                        caps(3, steps)
                    )
                    .unwrap()
                    .run
                );
                let inviscid =
                    ambient::simulate_with_cx(cx, &CLOUD, dt, steps, 0.1, source, caps(3, steps))
                        .unwrap();
                let viscous = ambient::simulate_viscous_with_cx(
                    cx,
                    &CLOUD,
                    &areas,
                    dt,
                    steps,
                    0.1,
                    source,
                    PseDiffusion::new(0.0, 1.0).unwrap(),
                    caps(3, steps),
                )
                .unwrap();
                assert_eq!(inviscid, viscous);
            }
            let empty =
                ambient::simulate_with_cx(cx, &[], -0.2, 2, 0.1, source, caps(0, 2)).unwrap();
            assert!(empty.run.particles.is_empty());
            assert_eq!(empty.run.work_units, 2);
            assert_eq!(empty.source_particle_updates, 0);
        },
    );
}

#[test]
fn g0_signed_impulse_changes_are_owned_by_the_external_source() {
    let (origin, velocity, rate, time) = ([0.2, -0.3], [0.4, -0.1], 0.13_f64, 1.2_f64);
    let source = flow(origin, velocity, [rate, 0.0, 0.0]);
    for last_gamma in [-0.5, 0.2] {
        let mut input = CLOUD;
        input[2].circulation = last_gamma;
        let before = moments(&input);
        let run = evolve(&input, time / 40.0, 40, source);
        let after = moments(&run.run.particles);
        let expected = [
            (rate * time).exp() * before[1]
                + before[0] * (velocity[0] - rate * origin[0]) * (rate * time).exp_m1() / rate,
            (-rate * time).exp() * before[2]
                - before[0] * (velocity[1] + rate * origin[1]) * (-rate * time).exp_m1() / rate,
        ];
        close(after[0], before[0], 1e-15);
        for axis in 0..2 {
            close(after[axis + 1], expected[axis], 8e-14);
            close(
                after[axis + 1] - before[axis + 1],
                run.ambient_moment_change.first_moment_m3_s[axis],
                8e-14,
            );
        }
        close(
            after[3] - before[3],
            run.ambient_moment_change.quadratic_moment_m4_s,
            8e-14,
        );
        assert!(run.ambient_moment_change.quadratic_moment_m4_s.abs() > 1e-3);
        let reverse = evolve(&run.run.particles, -time / 40.0, 40, source);
        for (actual, expected) in reverse.run.particles.iter().zip(input) {
            for axis in 0..2 {
                close(actual.pos[axis], expected.pos[axis], 2e-13);
            }
        }
    }
}

// Independent coupled ODE reference: no pair maps or production diffusion
// kernels are called. State is [x,y,Gamma] for each of three vortices.
fn derivative(y: [f64; 9], nu: f64) -> [f64; 9] {
    let areas = [0.5, 0.8, 1.1];
    let mut rhs = [0.0; 9];
    for i in 0..3 {
        let (x, z) = (y[3 * i] - 0.2, y[3 * i + 1] + 0.1);
        rhs[3 * i] = 0.3 + 0.21 * x + 0.7 * z;
        rhs[3 * i + 1] = -0.2 - 0.13 * x - 0.21 * z;
        for j in 0..3 {
            if i == j {
                continue;
            }
            let dx = y[3 * i] - y[3 * j];
            let dy = y[3 * i + 1] - y[3 * j + 1];
            let r2 = dx * dx + dy * dy;
            let speed = y[3 * j + 2] / (std::f64::consts::TAU * (r2 + 0.01));
            rhs[3 * i] -= speed * dy;
            rhs[3 * i + 1] += speed * dx;
            rhs[3 * i + 2] += 4.0 * nu / std::f64::consts::PI
                * (-r2).exp()
                * (areas[i] * y[3 * j + 2] - areas[j] * y[3 * i + 2]);
        }
    }
    rhs
}

fn reference(nu: f64) -> [f64; 9] {
    let mut y = core::array::from_fn(|i| {
        if i % 3 == 2 {
            CLOUD[i / 3].circulation
        } else {
            CLOUD[i / 3].pos[i % 3]
        }
    });
    let h = 0.8 / 2048.0;
    for _ in 0..2048 {
        let a = derivative(y, nu);
        let b = derivative(core::array::from_fn(|i| y[i] + 0.5 * h * a[i]), nu);
        let c = derivative(core::array::from_fn(|i| y[i] + 0.5 * h * b[i]), nu);
        let d = derivative(core::array::from_fn(|i| y[i] + h * c[i]), nu);
        y = core::array::from_fn(|i| y[i] + h / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]));
    }
    y
}

#[test]
fn g1_full_inviscid_and_viscous_coupling_converges_at_second_order() {
    let source = flow([0.2, -0.1], [0.3, -0.2], [0.21, 0.7, -0.13]);
    let areas = [0.5, 0.8, 1.1];
    for nu in [0.0, 0.2] {
        let exact = reference(nu);
        let law = PseDiffusion::new(nu, 1.0).unwrap();
        let errors = [8, 16, 32].map(|steps| {
            let run = with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
                ambient::simulate_viscous_with_cx(
                    cx,
                    &CLOUD,
                    &areas,
                    0.8 / steps as f64,
                    steps,
                    0.1,
                    source,
                    law,
                    caps(3, steps),
                )
            })
            .unwrap();
            close(moments(&run.run.particles)[0], 0.7, 3e-15);
            let pairs = if nu == 0.0 { 6 } else { 18 } * steps as u64;
            assert_eq!(run.run.pair_evaluations, pairs);
            assert_eq!(run.source_particle_updates, 6 * steps as u64);
            assert_eq!(run.run.work_units, 3 + 7 * steps as u64 + pairs);
            assert_eq!(run.run.budget.cost_charged, run.run.work_units);
            assert_eq!(run.run.budget.planned_cost, run.run.work_units);
            assert_eq!(run.run.peak_live_bytes, 3 * std::mem::size_of::<Particle>());
            (0..9)
                .map(|i| {
                    let p = run.run.particles[i / 3];
                    ((if i % 3 == 2 {
                        p.circulation
                    } else {
                        p.pos[i % 3]
                    }) - exact[i])
                        .abs()
                })
                .fold(0.0_f64, f64::max)
        });
        assert!(errors[2] < 1e-5, "nu={nu}, errors={errors:?}");
        for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
            assert!((3.8..4.2).contains(&ratio), "nu={nu}, errors={errors:?}");
        }
    }
}

#[test]
fn g3_world_frames_units_and_uniform_boosts_transform_the_actual_motion() {
    let source = flow([0.2, -0.1], [0.3, -0.2], [0.2, 0.3, -0.1]);
    let baseline = evolve(&CLOUD, 0.02, 10, source);
    let shift = [3.0, -7.0];
    let transform = |p: [f64; 2]| [-p[1] + shift[0], p[0] + shift[1]];
    let transformed = CLOUD.map(|p| Particle::new(transform(p.pos), p.circulation));
    let rotated = evolve(
        &transformed,
        0.02,
        10,
        flow(transform([0.2, -0.1]), [0.2, 0.3], [-0.2, 0.1, -0.3]),
    );
    for (actual, expected) in rotated.run.particles.iter().zip(&baseline.run.particles) {
        for axis in 0..2 {
            close(actual.pos[axis], transform(expected.pos)[axis], 3e-14);
        }
    }
    let m = baseline.ambient_moment_change;
    close(
        rotated.ambient_moment_change.first_moment_m3_s[0],
        -m.first_moment_m3_s[1],
        2e-14,
    );
    close(
        rotated.ambient_moment_change.first_moment_m3_s[1],
        m.first_moment_m3_s[0],
        2e-14,
    );
    close(
        rotated.ambient_moment_change.quadratic_moment_m4_s,
        m.quadratic_moment_m4_s
            + 2.0 * (-shift[0] * m.first_moment_m3_s[1] + shift[1] * m.first_moment_m3_s[0]),
        3e-13,
    );
    for (length, time) in [(1e-4, 3.0), (1e4, 0.3)] {
        let strength = length * length / time;
        let scaled =
            CLOUD.map(|p| Particle::new(p.pos.map(|x| length * x), strength * p.circulation));
        let result = with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
            ambient::simulate_with_cx(
                cx,
                &scaled,
                time * 0.02,
                10,
                length * 0.1,
                flow(
                    [length * 0.2, -length * 0.1],
                    [length / time * 0.3, -length / time * 0.2],
                    [0.2 / time, 0.3 / time, -0.1 / time],
                ),
                caps(3, 10),
            )
        })
        .unwrap();
        for (actual, expected) in result.run.particles.iter().zip(&baseline.run.particles) {
            for axis in 0..2 {
                close(actual.pos[axis] / length, expected.pos[axis], 3e-14);
            }
            close(actual.circulation / strength, expected.circulation, 2e-16);
        }
        for axis in 0..2 {
            close(
                result.ambient_moment_change.first_moment_m3_s[axis] / (strength * length),
                m.first_moment_m3_s[axis],
                2e-14,
            );
        }
        close(
            result.ambient_moment_change.quadratic_moment_m4_s / (strength * length * length),
            m.quadratic_moment_m4_s,
            3e-14,
        );
    }
    let boost = [0.8, -0.4];
    let isolated = evolve(&CLOUD, 0.02, 10, Default::default());
    let translated = evolve(&CLOUD, 0.02, 10, flow([0.0; 2], boost, [0.0; 3]));
    for (a, b) in translated.run.particles.iter().zip(isolated.run.particles) {
        for axis in 0..2 {
            close(a.pos[axis], b.pos[axis] + 0.2 * boost[axis], 5e-15);
        }
    }
}

#[test]
fn g4_resources_invalid_sources_and_late_numerical_refusals_are_atomic() {
    let source = flow([0.0; 2], [0.3, -0.1], [0.2, 0.3, -0.1]);
    let input = CLOUD;
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(AffineFlow2::new([bad, 0.0], [0.0; 2], [0.0; 3]).is_err());
        assert!(AffineFlow2::new([0.0; 2], [bad, 0.0], [0.0; 3]).is_err());
        assert!(AffineFlow2::new([0.0; 2], [0.0; 2], [0.0, bad, 0.0]).is_err());
    }
    let active = CancelGate::new_clock_free();
    with_cx(&active, Budget::INFINITE, |cx| {
        for limited in [
            VpmBudget::new(3, 2, 11, 72),
            VpmBudget::new(3, 2, 12, 71),
            VpmBudget::new(2, 2, 12, 72),
            VpmBudget::new(3, 1, 12, 72),
        ] {
            assert!(matches!(
                ambient::simulate_with_cx(cx, &input, 0.1, 2, 0.1, source, limited),
                Err(VpmError::CapExceeded { .. })
            ));
        }
        for bad in [f64::NAN, f64::INFINITY] {
            assert!(matches!(
                ambient::simulate_with_cx(cx, &input, bad, 2, 0.1, source, caps(3, 2)),
                Err(VpmError::InvalidScalar { .. })
            ));
        }
        assert!(matches!(
            ambient::simulate_with_cx(
                cx,
                &input,
                0.1,
                usize::MAX,
                0.1,
                source,
                VpmBudget::new(3, usize::MAX, u64::MAX, 72)
            ),
            Err(VpmError::PlanOverflow { .. })
        ));
        let law = PseDiffusion::new(0.0, 1.0).unwrap();
        assert!(
            ambient::simulate_viscous_with_cx(
                cx,
                &input,
                &[1.0, 0.0, 1.0],
                0.0,
                0,
                0.1,
                source,
                law,
                caps(3, 0)
            )
            .is_err()
        );
        assert!(
            ambient::simulate_viscous_with_cx(
                cx,
                &input,
                &[1.0; 3],
                -0.1,
                2,
                0.1,
                source,
                law,
                caps(3, 2)
            )
            .is_err()
        );
    });
    let planned = 29;
    assert!(matches!(
        with_cx(
            &active,
            Budget::INFINITE.with_cost_quota(planned - 1),
            |cx| ambient::simulate_with_cx(cx, &input, 0.1, 2, 0.1, source, caps(3, 2))
        ),
        Err(VpmError::ExecutionBudget(
            BudgetRefusal::CostPlanExceedsQuota {
                planned: 29,
                quota: 28
            }
        ))
    ));
    let admitted = with_cx(
        &active,
        Budget::INFINITE
            .with_cost_quota(planned)
            .with_poll_quota(100),
        |cx| ambient::simulate_with_cx(cx, &input, 0.1, 2, 0.1, source, caps(3, 2)),
    )
    .unwrap();
    for quota in 0..admitted.run.budget.polls_used {
        assert!(matches!(
            with_cx(&active, Budget::INFINITE.with_poll_quota(quota), |cx| {
                ambient::simulate_with_cx(cx, &input, 0.1, 2, 0.1, source, caps(3, 2))
            }),
            Err(VpmError::ExecutionBudget(
                BudgetRefusal::PollsExhausted { .. }
            ))
        ));
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        with_cx(
            &cancelled,
            Budget::INFINITE,
            |cx| ambient::simulate_with_cx(cx, &input, 0.1, 2, 0.1, source, caps(3, 2))
        ),
        Err(VpmError::ExecutionBudget(BudgetRefusal::Cancelled { .. }))
    ));
    let extreme = [
        Particle::new([1.0, 1.0], 1.0),
        Particle::new([1e308, 0.0], 0.0),
    ];
    let extensional = flow([0.0; 2], [0.0; 2], [1.0, 0.0, 0.0]);
    assert!(matches!(
        with_cx(&active, Budget::INFINITE, |cx| ambient::simulate_with_cx(
            cx,
            &extreme,
            1.0,
            1,
            0.1,
            extensional,
            caps(2, 1)
        )),
        Err(VpmError::NonFiniteResult {
            phase: "fs-vpm.ambient-source",
            particle: 1
        })
    ));
    assert_eq!(input, CLOUD);
    assert_eq!(extreme[1].pos, [1e308, 0.0]);
    assert_eq!(
        evolve(&input, 0.1, 2, source),
        evolve(&input, 0.1, 2, source)
    );
}

#[test]
fn g0_subnormal_durations_and_unresolved_source_motion_are_not_silenced() {
    let smallest = f64::from_bits(1);
    let tracer = Particle::new([0.0; 2], 0.0);
    let run = evolve(
        &[tracer],
        smallest,
        1,
        flow([0.0; 2], [1e308, 0.0], [0.0; 3]),
    );
    close(run.run.particles[0].pos[0], smallest * 1e308, 1e-31);
    let translated = evolve(
        &[Particle::new([1e308, 0.0], 0.0)],
        1.0,
        1,
        flow([-1e308, 1e308], [-1e308, 0.0], [0.0; 3]),
    );
    assert_eq!(translated.run.particles[0].pos, [0.0; 2]);
    let unresolved = [Particle::new([1e16, 0.0], 0.0)];
    assert!(matches!(
        with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
            ambient::simulate_with_cx(
                cx,
                &unresolved,
                1.0,
                1,
                0.1,
                flow([0.0; 2], [2.0, 0.0], [0.0; 3]),
                caps(1, 1),
            )
        }),
        Err(VpmError::NonFiniteResult {
            phase: "fs-vpm.ambient-position-resolution",
            ..
        })
    ));
    assert!(matches!(
        with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
            ambient::simulate_with_cx(
                cx,
                &[tracer],
                smallest,
                1,
                0.1,
                flow([0.0; 2], [0.0; 2], [1.0, 0.0, 0.0]),
                caps(1, 1),
            )
        }),
        Err(VpmError::NonFiniteResult {
            phase: "fs-vpm.ambient-map",
            ..
        })
    ));
}
