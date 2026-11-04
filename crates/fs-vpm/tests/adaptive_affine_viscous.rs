//! Independent finite-cloud ODE and analytic stiff-dipole checks. Temporal
//! refinement does not certify spatial PSE quadrature or continuum accuracy.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, Time, TimeSource};
use fs_vpm::conservative::adaptive::{
    self, AmbientAdaptiveBudget, AmbientViscousAdaptiveRun, ViscousAdaptiveOptions,
};
use fs_vpm::conservative::ambient::{self, AffineFlow2};
use fs_vpm::diffusion::PseDiffusion;
use fs_vpm::diffusion::adaptive::{AdaptiveBudget, AdaptiveError, AdaptiveOptions};
use fs_vpm::{VortexParticle as Particle, VpmBudget, VpmError};
use std::f64::consts::{PI, TAU};
use std::mem::size_of;
use std::sync::atomic::{AtomicU64, Ordering};

const CLOUD: [Particle; 4] = [
    Particle {
        pos: [-0.7, -0.2],
        circulation: 0.9,
    },
    Particle {
        pos: [0.5, 0.1],
        circulation: -0.5,
    },
    Particle {
        pos: [0.1, 0.8],
        circulation: 0.3,
    },
    Particle {
        pos: [-0.1, -0.6],
        circulation: -0.7,
    },
];
const AREAS: [f64; 4] = [0.8, 1.1, 0.7, 1.2];
const ORIGIN: [f64; 2] = [0.1, 0.1];
const VELOCITY: [f64; 2] = [0.2, -0.1];
const GRADIENT: [f64; 3] = [0.4, 0.55, -0.15];
const CORE: f64 = 0.12;
const NU: f64 = 0.1;
const WIDTH: f64 = 1.1;

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
                seed: 0x5649_5343_4f55,
                kernel_id: 0x4146_4649,
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

fn options(tolerance: f64) -> ViscousAdaptiveOptions {
    ViscousAdaptiveOptions {
        time: AdaptiveOptions {
            initial_dt_s: 0.6,
            min_dt_s: 1e-9,
            max_dt_s: 0.6,
            absolute_vorticity_tolerance_s_inv: tolerance,
            relative_tolerance: 0.0,
        },
        absolute_position_tolerance_m: tolerance,
    }
}

fn caps(n: usize) -> AmbientAdaptiveBudget {
    AmbientAdaptiveBudget {
        run: AdaptiveBudget {
            max_particles: n,
            max_attempts: 10_000,
            max_pair_evaluations: 20_000_000,
            max_work_units: 30_000_000,
            max_live_bytes: 3 * n * size_of::<Particle>(),
        },
        max_source_particle_updates: 1_000_000,
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance={tolerance}"
    );
}

fn bits(a: &[Particle], b: &[Particle]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_eq!(a.pos.map(f64::to_bits), b.pos.map(f64::to_bits));
        assert_eq!(a.circulation.to_bits(), b.circulation.to_bits());
    }
}

fn moments(particles: &[Particle]) -> [f64; 4] {
    let mut m = [0.0; 4];
    for p in particles {
        for (total, factor) in m.iter_mut().zip([
            1.0,
            p.pos[0],
            p.pos[1],
            p.pos[0] * p.pos[0] + p.pos[1] * p.pos[1],
        ]) {
            *total += p.circulation * factor;
        }
    }
    m
}

fn account_values(run: &AmbientViscousAdaptiveRun) -> [f64; 6] {
    let source = run.ambient_moment_change;
    let diffusion = run.diffusion_moment_change;
    [
        source.first_moment_m3_s[0],
        source.first_moment_m3_s[1],
        source.quadratic_moment_m4_s,
        diffusion.first_moment_m3_s[0],
        diffusion.first_moment_m3_s[1],
        diffusion.quadratic_moment_m4_s,
    ]
}

fn check_closure(initial: &[Particle], actual: &AmbientViscousAdaptiveRun, tolerance: f64) {
    let before = moments(initial);
    let after = moments(&actual.run.particles);
    let accounts = account_values(actual);
    for i in 0..3 {
        near(
            after[i + 1] - before[i + 1],
            accounts[i] + accounts[i + 3],
            tolerance,
        );
    }
    near(
        after[0] - before[0],
        actual.diffusion_moment_change.circulation_residual_m2_s,
        tolerance,
    );
}

// Independent simultaneous Biot-Savart + affine + PSE differential equations.
// The final six entries integrate the separate physical moment derivatives.
// No production pair, kernel-fraction, source-map or splitting owner is called.
fn rhs(state: &[f64; 18], compact: bool) -> [f64; 18] {
    let mut out = [0.0; 18];
    for i in 0..4 {
        let [x, y, gamma] = [state[3 * i], state[3 * i + 1], state[3 * i + 2]];
        let [a, b, c] = GRADIENT;
        let u = [
            VELOCITY[0] + a * (x - ORIGIN[0]) + b * (y - ORIGIN[1]),
            VELOCITY[1] + c * (x - ORIGIN[0]) - a * (y - ORIGIN[1]),
        ];
        out[3 * i] = u[0];
        out[3 * i + 1] = u[1];
        out[12] += gamma * u[0];
        out[13] += gamma * u[1];
        out[14] += 2.0 * gamma * (x * u[0] + y * u[1]);
        for j in 0..4 {
            if i == j {
                continue;
            }
            let dx = x - state[3 * j];
            let dy = y - state[3 * j + 1];
            let r2 = dx * dx + dy * dy;
            let velocity_factor = state[3 * j + 2] / (TAU * (r2 + CORE * CORE));
            out[3 * i] -= velocity_factor * dy;
            out[3 * i + 1] += velocity_factor * dx;
            let kappa = if compact {
                let r = r2.sqrt() / WIDTH;
                if r < 1.0 {
                    1008.0 * NU / (5.0 * PI * WIDTH.powi(4)) * (1.0 - r).powi(4) * (1.0 + 4.0 * r)
                } else {
                    0.0
                }
            } else {
                4.0 * NU / (PI * WIDTH.powi(4)) * (-r2 / (WIDTH * WIDTH)).exp()
            };
            out[3 * i + 2] += kappa * (AREAS[i] * state[3 * j + 2] - AREAS[j] * gamma);
        }
        out[15] += x * out[3 * i + 2];
        out[16] += y * out[3 * i + 2];
        out[17] += (x * x + y * y) * out[3 * i + 2];
    }
    out
}

fn reference(time: f64, steps: usize, compact: bool) -> [f64; 18] {
    let mut state = [0.0; 18];
    for i in 0..4 {
        state[3 * i] = CLOUD[i].pos[0];
        state[3 * i + 1] = CLOUD[i].pos[1];
        state[3 * i + 2] = CLOUD[i].circulation;
    }
    let dt = time / steps as f64;
    for _ in 0..steps {
        let shifted = |derivative: [f64; 18], fraction: f64| {
            std::array::from_fn(|i| state[i] + fraction * dt * derivative[i])
        };
        let k1 = rhs(&state, compact);
        let k2 = rhs(&shifted(k1, 0.5), compact);
        let k3 = rhs(&shifted(k2, 0.5), compact);
        let k4 = rhs(&shifted(k3, 1.0), compact);
        for i in 0..18 {
            state[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    state
}

#[test]
fn g1_joint_refinement_resolves_nonlinear_motion_and_separate_physical_moments() {
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for compact in [false, true] {
            let exact = reference(0.6, 4096, compact);
            let finer = reference(0.6, 8192, compact);
            for (a, b) in exact.iter().zip(finer) {
                near(*a, b, 3e-13);
            }
            let law = if compact {
                PseDiffusion::compact_wendland(NU, WIDTH).unwrap()
            } else {
                PseDiffusion::new(NU, WIDTH).unwrap()
            };
            let results = [1e-4, 1e-8].map(|tolerance| {
                adaptive::simulate_viscous_with_affine_flow_cx(
                    cx,
                    &CLOUD,
                    &AREAS,
                    0.6,
                    CORE,
                    law,
                    source(),
                    options(tolerance),
                    caps(4),
                )
                .unwrap_or_else(|error| {
                    panic!("compact={compact}, tolerance={tolerance}: {error:?}")
                })
            });
            let [loose, tight] = &results;
            let state_error = |run: &AmbientViscousAdaptiveRun| {
                run.run
                    .particles
                    .iter()
                    .enumerate()
                    .flat_map(|(i, p)| {
                        [
                            (p.pos[0] - exact[3 * i]).abs(),
                            (p.pos[1] - exact[3 * i + 1]).abs(),
                            (p.circulation - exact[3 * i + 2]).abs(),
                        ]
                    })
                    .fold(0.0, f64::max)
            };
            let account_error = |run: &AmbientViscousAdaptiveRun| {
                account_values(run)
                    .iter()
                    .zip(&exact[12..])
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max)
            };
            assert!(
                state_error(tight) < 3e-7,
                "state error={}",
                state_error(tight)
            );
            assert!(
                account_error(tight) < 3e-7,
                "moment error={}",
                account_error(tight)
            );
            assert!(state_error(loose) > 30.0 * state_error(tight));
            assert!(account_error(loose) > 30.0 * account_error(tight));
            assert!(tight.run.accepted_steps > loose.run.accepted_steps);
            assert!(tight.run.rejected_steps > 0);
            for actual in &results {
                check_closure(&CLOUD, actual, 3e-13);
                let r = &actual.run;
                assert_eq!(r.elapsed_s.to_bits(), 0.6_f64.to_bits());
                assert_eq!(r.attempts, r.accepted_steps + r.rejected_steps);
                assert_eq!(r.pair_evaluations, 144 * r.attempts as u64);
                assert_eq!(r.fraction_evaluations, 84 * r.attempts as u64);
                assert_eq!(actual.source_particle_updates, 28 * r.attempts as u64);
                assert_eq!(r.work_units, 4 + 193 * r.attempts as u64);
                assert_eq!(r.budget.cost_charged, r.work_units);
                assert_eq!(r.peak_live_bytes, 12 * size_of::<Particle>());
                assert!(actual.diffusion_moment_change.first_moment_m3_s[1].abs() > 1e-4);
                assert!(
                    actual
                        .diffusion_moment_change
                        .first_moment_residual_m3_s
                        .iter()
                        .all(|v| v.abs() < 1e-13)
                );
                assert!(
                    actual
                        .diffusion_moment_change
                        .quadratic_moment_residual_m4_s
                        .abs()
                        < 1e-13
                );
            }
        }
    });
}

#[test]
fn g1_stiff_rotating_dipole_retains_induced_transport_hidden_by_step_doubling() {
    let initial = [
        Particle::new([0.0, 0.5], 1.0),
        Particle::new([0.0, -0.5], -1.0),
    ];
    let nu = 1000.0;
    let lambda = 8.0 * nu / (PI * std::f64::consts::E);
    let omega: f64 = 0.7;
    let (sin, cos) = omega.sin_cos();
    let flow = AffineFlow2::new([0.0; 2], [0.0; 2], [0.0, -omega, omega]).unwrap();
    // In a frame rotating with the imposed source, self-induced translation is
    // e^(-lambda*t)/(2*pi*(d^2+core^2)) in x. Rotate its integral back at t=1.
    let displacement = -(-lambda).exp_m1() / (lambda * TAU * 1.01);
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let results = [1e-5, 1e-8].map(|tolerance| {
            let settings = ViscousAdaptiveOptions {
                time: AdaptiveOptions {
                    initial_dt_s: 1.0,
                    max_dt_s: 1.0,
                    ..options(tolerance).time
                },
                ..options(tolerance)
            };
            let actual = adaptive::simulate_viscous_with_affine_flow_cx(
                cx,
                &initial,
                &[1.0; 2],
                1.0,
                0.1,
                PseDiffusion::new(nu, 1.0).unwrap(),
                flow,
                settings,
                caps(2),
            )
            .unwrap();
            assert!(actual.run.rejected_steps > 0);
            let mut error: f64 = 0.0;
            for (p, original) in actual.run.particles.iter().zip(initial) {
                let exact = [
                    cos * displacement - sin * original.pos[1],
                    sin * displacement + cos * original.pos[1],
                ];
                error = error.max((p.pos[0] - exact[0]).hypot(p.pos[1] - exact[1]));
                assert!(p.circulation.abs() < 1e-13);
            }
            check_closure(&initial, &actual, 3e-13);
            (error, actual.run.accepted_steps)
        });
        assert!(
            results[1].0 < 3e-8,
            "stiff rotating displacement error={}",
            results[1].0
        );
        assert!(results[0].0 > 50.0 * results[1].0);
        assert!(results[1].1 > results[0].1);
    });
}

#[test]
fn g0_zero_source_preserves_native_adaptation_and_fine_steps_use_original_owners() {
    let law = PseDiffusion::new(NU, WIDTH).unwrap();
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let native = adaptive::simulate_viscous_with_cx(
            cx,
            &CLOUD,
            &AREAS,
            0.6,
            CORE,
            law,
            options(1e-7),
            caps(4).run,
        )
        .unwrap();
        let observed = adaptive::simulate_viscous_with_affine_flow_cx(
            cx,
            &CLOUD,
            &AREAS,
            0.6,
            CORE,
            law,
            AffineFlow2::default(),
            options(1e-7),
            AmbientAdaptiveBudget {
                max_source_particle_updates: 0,
                ..caps(4)
            },
        )
        .unwrap();
        assert_eq!(native, observed.run);
        assert_eq!(observed.source_particle_updates, 0);
        assert_eq!(observed.ambient_moment_change, Default::default());
        check_closure(&CLOUD, &observed, 3e-13);
        let h = 0.02;
        let settings = ViscousAdaptiveOptions {
            time: AdaptiveOptions {
                initial_dt_s: h,
                min_dt_s: h,
                max_dt_s: h,
                ..options(1.0).time
            },
            ..options(1.0)
        };
        for nu in [0.0, NU] {
            let law = PseDiffusion::new(nu, WIDTH).unwrap();
            let actual = adaptive::simulate_viscous_with_affine_flow_cx(
                cx,
                &CLOUD,
                &AREAS,
                h,
                CORE,
                law,
                source(),
                settings,
                caps(4),
            )
            .unwrap();
            let original = ambient::simulate_viscous_with_cx(
                cx,
                &CLOUD,
                &AREAS,
                h / 2.0,
                2,
                CORE,
                source(),
                law,
                VpmBudget::new(4, 2, u64::MAX, 4 * size_of::<Particle>()),
            )
            .unwrap();
            assert_eq!(actual.run.attempts, 1);
            assert_eq!(actual.run.accepted_steps, 1);
            for (p, q) in actual.run.particles.iter().zip(&original.run.particles) {
                near(p.pos[0], q.pos[0], 5e-16);
                near(p.pos[1], q.pos[1], 5e-16);
                near(p.circulation, q.circulation, 5e-16);
            }
            for i in 0..2 {
                near(
                    actual.ambient_moment_change.first_moment_m3_s[i],
                    original.ambient_moment_change.first_moment_m3_s[i],
                    5e-16,
                );
            }
            near(
                actual.ambient_moment_change.quadratic_moment_m4_s,
                original.ambient_moment_change.quadratic_moment_m4_s,
                5e-16,
            );
            check_closure(&CLOUD, &actual, 3e-15);
            if nu == 0.0 {
                assert_eq!(actual.diffusion_moment_change, Default::default());
                assert_eq!(actual.source_particle_updates, 20);
                assert_eq!(actual.run.pair_evaluations, 36);
                assert_eq!(actual.run.fraction_evaluations, 0);
                assert_eq!(actual.run.work_units, 72);
            }
        }
    });
}

#[test]
fn g3_frame_and_length_units_transform_the_separate_accounts_covariantly() {
    let angle: f64 = 0.43;
    let (s, c) = angle.sin_cos();
    let rotate = |x: [f64; 2]| [c * x[0] - s * x[1], s * x[0] + c * x[1]];
    let shift = [3.0, -2.0];
    let [a, b, d] = GRADIENT;
    let gradient = [
        a * (c * c - s * s) - (b + d) * s * c,
        2.0 * a * s * c + b * c * c - d * s * s,
        2.0 * a * s * c + d * c * c - b * s * s,
    ];
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let original = adaptive::simulate_viscous_with_affine_flow_cx(
            cx,
            &CLOUD,
            &AREAS,
            0.6,
            CORE,
            PseDiffusion::new(NU, WIDTH).unwrap(),
            source(),
            options(1e-7),
            caps(4),
        )
        .unwrap();
        for length in [1.0, 1024.0] {
            let transform = |x: [f64; 2]| {
                let r = rotate(x);
                [length * (r[0] + shift[0]), length * (r[1] + shift[1])]
            };
            let cloud =
                CLOUD.map(|p| Particle::new(transform(p.pos), length * length * p.circulation));
            let v = rotate(VELOCITY).map(|v| length * v);
            let actual = adaptive::simulate_viscous_with_affine_flow_cx(
                cx,
                &cloud,
                &AREAS.map(|a| length * length * a),
                0.6,
                length * CORE,
                PseDiffusion::new(length * length * NU, length * WIDTH).unwrap(),
                AffineFlow2::new(transform(ORIGIN), v, gradient).unwrap(),
                ViscousAdaptiveOptions {
                    absolute_position_tolerance_m: length * 1e-7,
                    ..options(1e-7)
                },
                caps(4),
            )
            .unwrap();
            assert_eq!(actual.run.attempts, original.run.attempts);
            assert_eq!(actual.run.accepted_steps, original.run.accepted_steps);
            for (p, q) in actual.run.particles.iter().zip(&original.run.particles) {
                let expected = transform(q.pos);
                near(p.pos[0] / length, expected[0] / length, 3e-13);
                near(p.pos[1] / length, expected[1] / length, 3e-13);
                near(p.circulation / (length * length), q.circulation, 3e-13);
            }
            let old = account_values(&original);
            let new = account_values(&actual);
            for offset in [0, 3] {
                let first = rotate([old[offset], old[offset + 1]]);
                let circulation_change = if offset == 3 {
                    original.diffusion_moment_change.circulation_residual_m2_s
                } else {
                    0.0
                };
                near(
                    new[offset] / length.powi(3),
                    first[0] + shift[0] * circulation_change,
                    3e-13,
                );
                near(
                    new[offset + 1] / length.powi(3),
                    first[1] + shift[1] * circulation_change,
                    3e-13,
                );
                near(
                    new[offset + 2] / length.powi(4),
                    old[offset + 2]
                        + 2.0 * (shift[0] * first[0] + shift[1] * first[1])
                        + (shift[0] * shift[0] + shift[1] * shift[1]) * circulation_change,
                    3e-12,
                );
            }
        }
    });
}

#[test]
fn g0_singleton_subnormal_source_time_and_pair_roundoff_have_honest_accounts() {
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let single = [Particle::new([0.2, -0.1], 0.7)];
        let flow = AffineFlow2::new([0.0; 2], [0.4, -0.2], [0.0; 3]).unwrap();
        let settings = ViscousAdaptiveOptions {
            time: AdaptiveOptions {
                initial_dt_s: 0.9,
                min_dt_s: 0.2,
                max_dt_s: 0.9,
                ..options(1e-8).time
            },
            ..options(1e-8)
        };
        let actual = adaptive::simulate_viscous_with_affine_flow_cx(
            cx,
            &single,
            &[1.0],
            1.0,
            CORE,
            PseDiffusion::new(1.0, 1.0).unwrap(),
            flow,
            settings,
            caps(1),
        )
        .unwrap();
        near(actual.run.particles[0].pos[0], 0.6, 3e-16);
        near(actual.run.particles[0].pos[1], -0.3, 3e-16);
        assert_eq!(actual.run.accepted_steps, 2);
        assert_eq!(actual.run.elapsed_s, 1.0);
        assert_eq!(actual.source_particle_updates, 10);
        assert_eq!(actual.run.work_units, 23);
        assert_eq!(actual.run.pair_evaluations, 0);
        assert_eq!(actual.run.fraction_evaluations, 0);
        assert_eq!(actual.diffusion_moment_change, Default::default());
        check_closure(&single, &actual, 3e-16);
        let h = f64::from_bits(1);
        let tiny = ViscousAdaptiveOptions {
            time: AdaptiveOptions {
                initial_dt_s: h,
                min_dt_s: h,
                max_dt_s: h,
                ..options(1e-8).time
            },
            ..options(1e-8)
        };
        let initial = [Particle::new([0.0; 2], 1.0)];
        let actual = adaptive::simulate_viscous_with_affine_flow_cx(
            cx,
            &initial,
            &[1.0],
            h,
            CORE,
            PseDiffusion::new(0.0, 1.0).unwrap(),
            AffineFlow2::new([0.0; 2], [1e300, 0.0], [0.0; 3]).unwrap(),
            tiny,
            caps(1),
        )
        .unwrap();
        near(actual.run.particles[0].pos[0], h * 1e300, 1e-38);
        assert!(actual.run.particles[0].pos[0] > 0.0);
        assert_eq!(actual.run.elapsed_s.to_bits(), h.to_bits());
        // The native exchange is smaller than either circulation's ULP. The
        // actual PSE moment change is zero; its requested transfer loss must
        // be visible as a representation defect, never invented source work.
        let large = [
            Particle::new([0.0, 0.5], 1e16),
            Particle::new([0.0, -0.5], 1e16),
        ];
        let h = 1e-12;
        let settings = ViscousAdaptiveOptions {
            time: AdaptiveOptions {
                initial_dt_s: h,
                min_dt_s: h,
                max_dt_s: h,
                ..options(1.0).time
            },
            ..options(1.0)
        };
        let actual = adaptive::simulate_viscous_with_affine_flow_cx(
            cx,
            &large,
            &[1.0, 2.0],
            h,
            1e8,
            PseDiffusion::new(1e-4, 1.0).unwrap(),
            AffineFlow2::default(),
            settings,
            caps(2),
        )
        .unwrap();
        for p in &actual.run.particles {
            assert_eq!(p.circulation, 1e16);
        }
        assert_eq!(actual.diffusion_moment_change.first_moment_m3_s, [0.0; 2]);
        assert_eq!(actual.diffusion_moment_change.quadratic_moment_m4_s, 0.0);
        assert_eq!(
            actual.diffusion_moment_change.circulation_residual_m2_s,
            0.0
        );
        assert!(actual.diffusion_moment_change.first_moment_residual_m3_s[1].abs() > 0.1);
        assert_eq!(actual.ambient_moment_change, Default::default());
    });
}

#[test]
fn g4_whole_window_caps_include_rejected_source_and_witness_work() {
    let saved = CLOUD;
    let law = PseDiffusion::new(NU, WIDTH).unwrap();
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let run = |limits| {
            adaptive::simulate_viscous_with_affine_flow_cx(
                cx,
                &CLOUD,
                &AREAS,
                0.6,
                CORE,
                law,
                source(),
                options(1e-7),
                limits,
            )
        };
        let good = run(caps(4)).unwrap();
        assert!(good.run.rejected_steps > 0);
        let exact = AmbientAdaptiveBudget {
            run: AdaptiveBudget {
                max_particles: 4,
                max_attempts: good.run.attempts,
                max_pair_evaluations: good.run.pair_evaluations,
                max_work_units: good.run.work_units,
                max_live_bytes: good.run.peak_live_bytes,
            },
            max_source_particle_updates: good.source_particle_updates,
        };
        let mut expected = good.clone();
        expected.run.budget.planned_cost = exact.run.max_work_units;
        assert_eq!(expected, run(exact).unwrap());
        for limited in [
            AmbientAdaptiveBudget {
                max_source_particle_updates: exact.max_source_particle_updates - 1,
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
        ] {
            assert!(matches!(
                run(limited),
                Err(AdaptiveError::Engine(VpmError::CapExceeded { .. }))
            ));
            bits(&CLOUD, &saved);
        }
        assert_eq!(expected, run(exact).unwrap());
    });
}

#[test]
fn g1_compact_affine_updates_source_geometry_and_matches_dense_physical_accounts() {
    let compressed = [
        Particle::new([0.0, 0.75], 1.0),
        Particle::new([0.0, -0.75], -1.0),
    ];
    let strain = AffineFlow2::new([0.0; 2], [0.0; 2], [0.8, 0.0, 0.0]).unwrap();
    let fixed = ViscousAdaptiveOptions {
        time: AdaptiveOptions {
            initial_dt_s: 1.0,
            min_dt_s: 1.0,
            max_dt_s: 1.0,
            ..options(1.0).time
        },
        ..options(1.0)
    };
    // In the single-step fixture the pair is outside support after the first
    // fine native step, then the middle affine map brings it inside support.
    // Reusing that adjacent diffusion graph loses a real strength exchange.
    assert!(1.5 * (-0.2_f64).exp() > 1.0);
    assert!(1.5 * (-0.6_f64).exp() < 1.0);
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for (initial, areas, time, law, flow, settings) in [
            (
                compressed.as_slice(),
                [1.0; 2].as_slice(),
                1.0,
                PseDiffusion::compact_wendland(0.2, 1.0).unwrap(),
                strain,
                fixed,
            ),
            (
                CLOUD.as_slice(),
                AREAS.as_slice(),
                0.6,
                PseDiffusion::compact_wendland(NU, WIDTH).unwrap(),
                source(),
                options(1e-7),
            ),
        ] {
            let n = initial.len();
            let mut budget = caps(n);
            budget.run.max_live_bytes = 1024 * 1024;
            let neighbors = adaptive::CompactNeighborBudget {
                max_edges: n * (n - 1) / 2,
                max_candidate_visits: 1_000_000,
            };
            let dense = adaptive::simulate_viscous_with_affine_flow_cx(
                cx, initial, areas, time, CORE, law, flow, settings, budget,
            )
            .unwrap();
            let compact = adaptive::simulate_viscous_compact_with_affine_flow_cx(
                cx, initial, areas, time, CORE, law, flow, settings, budget, neighbors,
            )
            .unwrap();
            bits(&dense.run.particles, &compact.run.particles);
            assert_eq!(compact.run.accepted_steps, dense.run.accepted_steps);
            assert_eq!(compact.run.rejected_steps, dense.run.rejected_steps);
            assert_eq!(
                compact.run.max_accepted_error_ratio.to_bits(),
                dense.run.max_accepted_error_ratio.to_bits()
            );
            assert_eq!(compact.ambient_moment_change, dense.ambient_moment_change);
            assert_eq!(
                compact.diffusion_moment_change,
                dense.diffusion_moment_change
            );
            assert_eq!(
                compact.source_particle_updates,
                dense.source_particle_updates
            );
            assert_eq!(compact.graph_rebuilds, 7 * compact.run.attempts);
            assert!(compact.candidate_visits > 0);
            assert_eq!(
                compact.run.pair_evaluations,
                5 * (n * (n - 1)) as u64 * compact.run.attempts as u64
                    + 2 * compact.run.fraction_evaluations
            );
            assert_eq!(compact.run.budget.cost_charged, compact.run.work_units);
            assert!(compact.run.peak_live_bytes > dense.run.peak_live_bytes);
            check_closure(initial, &compact, 3e-13);
            if n == 2 {
                assert_eq!(compact.run.attempts, 1);
                assert!(compact.run.particles[0].circulation < 0.999);
                assert!(matches!(
                    adaptive::simulate_viscous_compact_with_affine_flow_cx(
                        cx,
                        initial,
                        areas,
                        time,
                        CORE,
                        law,
                        flow,
                        settings,
                        budget,
                        adaptive::CompactNeighborBudget {
                            max_edges: 0,
                            ..neighbors
                        },
                    ),
                    Err(AdaptiveError::Engine(VpmError::CapExceeded {
                        resource: "moving compact edges",
                        ..
                    }))
                ));
            } else {
                assert!(compact.run.rejected_steps > 0);
            }
            // Inactive source uses the exact current compact native driver,
            // including its six graph rebuilds per trial and execution receipt.
            let original = adaptive::simulate_viscous_compact_with_cx(
                cx, initial, areas, time, CORE, law, settings, budget.run, neighbors,
            )
            .unwrap();
            let zero = adaptive::simulate_viscous_compact_with_affine_flow_cx(
                cx,
                initial,
                areas,
                time,
                CORE,
                law,
                AffineFlow2::default(),
                settings,
                AmbientAdaptiveBudget {
                    max_source_particle_updates: 0,
                    ..budget
                },
                neighbors,
            )
            .unwrap();
            assert_eq!(original.run, zero.run);
            assert_eq!(original.graph_rebuilds, zero.graph_rebuilds);
            assert_eq!(original.candidate_visits, zero.candidate_visits);
            assert_eq!(zero.ambient_moment_change, Default::default());
        }
    });
}

#[test]
fn g4_compact_source_caps_cover_all_graphs_sources_and_rejected_trials() {
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let law = PseDiffusion::compact_wendland(NU, WIDTH).unwrap();
        let mut budget = caps(4);
        budget.run.max_live_bytes = 1024 * 1024;
        let graph = adaptive::CompactNeighborBudget {
            max_edges: 6,
            max_candidate_visits: 1_000_000,
        };
        let run = |budget, graph| {
            adaptive::simulate_viscous_compact_with_affine_flow_cx(
                cx,
                &CLOUD,
                &AREAS,
                0.6,
                CORE,
                law,
                source(),
                options(1e-7),
                budget,
                graph,
            )
        };
        let good = run(budget, graph).unwrap();
        assert!(good.run.rejected_steps > 0);
        let exact = AmbientAdaptiveBudget {
            run: AdaptiveBudget {
                max_particles: 4,
                max_attempts: good.run.attempts,
                max_pair_evaluations: good.run.pair_evaluations,
                max_work_units: good.run.work_units,
                max_live_bytes: good.run.peak_live_bytes,
            },
            max_source_particle_updates: good.source_particle_updates,
        };
        let graph = adaptive::CompactNeighborBudget {
            max_candidate_visits: good.candidate_visits,
            ..graph
        };
        let mut expected = good.clone();
        expected.run.budget.planned_cost = exact.run.max_work_units;
        assert_eq!(expected, run(exact, graph).unwrap());
        for budget in [
            AmbientAdaptiveBudget {
                max_source_particle_updates: exact.max_source_particle_updates - 1,
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
                    max_pair_evaluations: exact.run.max_pair_evaluations - 1,
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
        ] {
            assert!(matches!(
                run(budget, graph),
                Err(AdaptiveError::Engine(VpmError::CapExceeded { .. }))
            ));
        }
        assert!(matches!(
            run(
                exact,
                adaptive::CompactNeighborBudget {
                    max_candidate_visits: graph.max_candidate_visits - 1,
                    ..graph
                }
            ),
            Err(AdaptiveError::Engine(VpmError::CapExceeded { .. }))
        ));
        assert_eq!(expected, run(exact, graph).unwrap());
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
fn g4_source_diffusion_and_publication_cancellation_leave_inputs_unchanged() {
    let input: Vec<_> = (0..12)
        .map(|i| {
            let (s, c) = (TAU * i as f64 / 12.0).sin_cos();
            Particle::new([2.0 * c, 2.0 * s], if i % 2 == 0 { 0.01 } else { -0.01 })
        })
        .collect();
    let saved = input.clone();
    let areas = [1.0; 12];
    let execution = Budget::INFINITE.with_deadline(Time::from_nanos(1_000_000));
    let h = 0.01;
    let settings = ViscousAdaptiveOptions {
        time: AdaptiveOptions {
            initial_dt_s: h,
            min_dt_s: h,
            max_dt_s: h,
            ..options(1e-2).time
        },
        ..options(1e-2)
    };
    let run = |gate: &CancelGate, clock: &TripClock<'_>| {
        with_cx(gate, execution, |cx| {
            let timed = Cx::new(gate, cx.arena(), cx.stream_key(), execution, cx.mode())
                .with_time_source(clock);
            adaptive::simulate_viscous_with_affine_flow_cx(
                &timed,
                &input,
                &areas,
                h,
                CORE,
                PseDiffusion::new(NU, WIDTH).unwrap(),
                source(),
                settings,
                caps(12),
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
    let (mut diffusion_seen, mut source_seen, mut publication_seen) = (false, false, false);
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
                diffusion_seen |= phase == "fs-vpm.diffusion-pair";
                source_seen |= phase == "fs-vpm.ambient-source";
                publication_seen |= phase == "fs-vpm.viscous-adaptive-publication";
            }
            other => panic!("cancel_on={cancel_on}, result={other:?}"),
        }
        bits(&input, &saved);
    }
    assert!(diffusion_seen && source_seen && publication_seen);
    assert_eq!(good, run(&gate, &counter).unwrap());
}
