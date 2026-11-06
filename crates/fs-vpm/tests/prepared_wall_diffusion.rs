//! G0/G1/G3/G4: compact wall acceleration retains the native wall operator,
//! singleton image losses, physical decay and transactional prepared ownership.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_vpm::diffusion::{PreparationBudget, PreparedPse, PseDiffusion};
use fs_vpm::planar_wall::{self, PlanarWall2};
use fs_vpm::{VortexParticle, VpmBudget, VpmError};
use std::f64::consts::PI;
use std::mem::size_of;

fn with_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x5741_4c4c,
                kernel_id: 0x5053_45,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn budget(n: usize, steps: usize) -> VpmBudget {
    VpmBudget::new(n, steps, u64::MAX, usize::MAX)
}

fn wall() -> PlanarWall2 {
    PlanarWall2::new([0.0; 2], [0.0, 1.0]).unwrap()
}

fn prepare(cx: &Cx<'_>, p: &[VortexParticle], areas: &[f64], nu: f64) -> PreparedPse {
    PreparedPse::prepare_with_cx(
        cx,
        p,
        areas,
        PseDiffusion::compact_wendland(nu, 1.0).unwrap(),
        PreparationBudget::new(p.len(), usize::MAX, 1_000_000, 1_000_000, usize::MAX),
    )
    .unwrap()
}

fn same_bits(a: &[VortexParticle], b: &[VortexParticle]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_eq!(a.pos.map(f64::to_bits), b.pos.map(f64::to_bits));
        assert_eq!(a.circulation.to_bits(), b.circulation.to_bits());
    }
}

#[test]
fn g0_g3_irregular_rotated_cloud_is_bitwise_dense_after_wall_and_viscosity_changes() {
    let origin = [0.25, -0.5];
    let normal = [0.6, 0.8];
    let bounded = PlanarWall2::new(origin, normal).unwrap();
    let local = [
        [-0.8, 0.02],
        [0.0, 0.3],
        [0.7, 0.18],
        [-0.2, 0.75],
        [1.6, 0.08],
        [0.0, 1.5],
        [4.0, 0.1],
    ];
    let gamma = [-0.4, 0.8, 0.0, -0.0, 1e-11, -0.3, 0.7];
    let p: Vec<_> = local
        .iter()
        .zip(gamma)
        .map(|(q, gamma)| {
            VortexParticle::new(
                [
                    origin[0] - 0.8 * q[0] + 0.6 * q[1],
                    origin[1] + 0.6 * q[0] + 0.8 * q[1],
                ],
                gamma,
            )
        })
        .collect();
    let areas = [0.03, 0.7, 0.2, 1.3, 0.05, 0.4, 0.9];
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let mut graph = prepare(cx, &p, &areas, 0.17);
        let positions = graph.positions().as_ptr();
        let edges = graph.edges().as_ptr();
        let prepared = *graph.preparation();
        for nu in [0.17, 0.0, 0.31] {
            graph.rebind_viscosity(nu).unwrap();
            for boundary in [bounded, PlanarWall2::new([-0.35, -1.3], normal).unwrap()] {
                for dt in [0.07, 40.0] {
                    let dense = planar_wall::diffuse_with_cx(
                        cx,
                        &p,
                        &areas,
                        dt,
                        9,
                        graph.law(),
                        boundary,
                        budget(p.len(), 9),
                    )
                    .unwrap();
                    let sparse = graph
                        .diffuse_at_wall_with_cx(cx, &p, dt, 9, boundary, budget(p.len(), 9))
                        .unwrap();
                    same_bits(&sparse.particles, &dense.particles);
                }
            }
        }
        assert_eq!(positions, graph.positions().as_ptr());
        assert_eq!(edges, graph.edges().as_ptr());
        assert_eq!(prepared, *graph.preparation());
    });
}

#[test]
fn g1_edgeless_graph_retains_analytic_self_image_decay_and_inactive_images() {
    let p = [
        VortexParticle::new([0.0, 0.1], -0.7),
        VortexParticle::new([3.0, 0.5], -0.0),
        VortexParticle::new([6.0, 1.0], 0.8),
    ];
    let areas = [0.6, 0.3, 0.9];
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let mut graph = prepare(cx, &p, &areas, 0.12);
        assert!(graph.edges().is_empty());
        let run = graph
            .diffuse_at_wall_with_cx(cx, &p, 0.05, 7, wall(), budget(3, 7))
            .unwrap();
        let q = 0.2_f64;
        let rate =
            2.0 * areas[0] * 1008.0 / (5.0 * PI) * 0.12 * (1.0 - q).powi(4) * (1.0 + 4.0 * q);
        let exact = p[0].circulation * (-rate * 0.35).exp();
        assert!((run.particles[0].circulation - exact).abs() < 3e-14);
        same_bits(&run.particles[1..], &p[1..]);
        assert_eq!(run.pair_evaluations, 42);
        assert_eq!(run.work_units, 55);
        assert_eq!(
            run.peak_live_bytes,
            3 * (size_of::<VortexParticle>() + 2 * size_of::<f64>())
        );
        graph.rebind_viscosity(f64::MAX).unwrap();
        let stiff = graph
            .diffuse_at_wall_with_cx(cx, &p, 1.0, 1, wall(), budget(3, 1))
            .unwrap();
        assert_eq!(
            stiff.particles[0].circulation.to_bits(),
            (-0.0_f64).to_bits()
        );
        same_bits(&stiff.particles[1..], &p[1..]);
    });
}

#[test]
fn g0_sparse_wall_work_scales_with_supported_edges_and_caps_are_exact() {
    let p: Vec<_> = (0..256)
        .map(|i| VortexParticle::new([0.75 * f64::from(i), 0.1], 0.1 + f64::from(i % 7) * 0.03))
        .collect();
    let areas = vec![0.2; p.len()];
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let graph = prepare(cx, &p, &areas, 0.08);
        assert_eq!(graph.edges().len(), 255);
        let run = graph
            .diffuse_at_wall_with_cx(cx, &p, 0.03, 2, wall(), budget(256, 2))
            .unwrap();
        let dense = planar_wall::diffuse_with_cx(
            cx,
            &p,
            &areas,
            0.03,
            2,
            graph.law(),
            wall(),
            budget(256, 2),
        )
        .unwrap();
        same_bits(&run.particles, &dense.particles);
        assert_eq!(run.pair_evaluations, 2 * (2 * 256 + 4 * 255));
        assert_eq!(run.work_units, 256 + 256 + 255 + run.pair_evaluations + 2);
        assert_eq!(run.budget.cost_charged, run.work_units);
        assert_eq!(run.peak_live_bytes, 40 * (256 + 255));
        assert!(dense.pair_evaluations > 85 * run.pair_evaluations);
        let pairs = VpmBudget::new(256, 2, run.pair_evaluations - 1, usize::MAX);
        assert_eq!(
            graph.diffuse_at_wall_with_cx(cx, &p, 0.03, 2, wall(), pairs),
            Err(VpmError::CapExceeded {
                resource: "pair evaluations",
                required: run.pair_evaluations.into(),
                maximum: (run.pair_evaluations - 1).into()
            })
        );
        let bytes = VpmBudget::new(256, 2, u64::MAX, run.peak_live_bytes - 1);
        assert_eq!(
            graph.diffuse_at_wall_with_cx(cx, &p, 0.03, 2, wall(), bytes),
            Err(VpmError::CapExceeded {
                resource: "logical live bytes",
                required: run.peak_live_bytes as u128,
                maximum: (run.peak_live_bytes - 1) as u128
            })
        );
    });
}

#[test]
fn g1_cached_loss_preserves_representable_transfer_after_exponential_underflow() {
    let p = [
        VortexParticle::new([0.0, 0.01], 1e308),
        VortexParticle::new([0.0, 0.7], 0.0),
    ];
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let graph = prepare(cx, &p, &[1.0; 2], 1.0);
        let run = graph
            .diffuse_at_wall_with_cx(cx, &p, 14.0, 1, wall(), budget(2, 1))
            .unwrap();
        let dense = planar_wall::diffuse_with_cx(
            cx,
            &p,
            &[1.0; 2],
            14.0,
            1,
            graph.law(),
            wall(),
            budget(2, 1),
        )
        .unwrap();
        same_bits(&run.particles, &dense.particles);
        // The first self-loss has exp(-rate*dt/2)==0, but multiplying the
        // original large circulation in logarithmic form leaves a finite
        // amount that transfers to the farther site before the reverse loss.
        assert!(run.particles[1].circulation > 0.0);
        assert!(run.particles[1].circulation < 1e-60);
    });
}

#[test]
fn g4_identity_domains_and_cancelled_runs_keep_the_prepared_graph_reusable() {
    let p: Vec<_> = (0..17)
        .map(|i| VortexParticle::new([0.05 * f64::from(i), 0.1], f64::from(i) - 8.0))
        .collect();
    let saved = p.clone();
    let graph = with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        prepare(cx, &p, &[0.1; 17], 0.2)
    });
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        for (dt, steps) in [(0.0, 3), (0.1, 0)] {
            let run = graph
                .diffuse_at_wall_with_cx(cx, &p, dt, steps, wall(), budget(17, steps))
                .unwrap();
            same_bits(&run.particles, &p);
            assert_eq!(run.pair_evaluations, 0);
            assert_eq!(run.work_units, (p.len() + steps) as u64);
            assert_eq!(run.peak_live_bytes, p.len() * size_of::<VortexParticle>());
        }
        let outside = PlanarWall2::new([0.0, 0.1], [0.0, 1.0]).unwrap();
        assert!(matches!(
            graph.diffuse_at_wall_with_cx(cx, &p, 0.0, 0, outside, budget(17, 0)),
            Err(VpmError::InvalidScalar {
                field: "wall clearance",
                ..
            })
        ));
        let mut changed = p.clone();
        changed[0].pos[0] = -0.0;
        assert!(matches!(
            graph.diffuse_at_wall_with_cx(cx, &changed, 0.0, 0, wall(), budget(17, 0)),
            Err(VpmError::InvalidScalar {
                field: "prepared position.x",
                ..
            })
        ));
    });
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let result = with_cx(&cancelled, Budget::INFINITE, |cx| {
        graph.diffuse_at_wall_with_cx(cx, &p, 0.1, 1, wall(), budget(17, 1))
    });
    assert!(matches!(
        result,
        Err(VpmError::ExecutionBudget(BudgetRefusal::Cancelled { .. }))
    ));
    let exhausted = with_cx(
        &CancelGate::new_clock_free(),
        Budget::INFINITE.with_poll_quota(6),
        |cx| graph.diffuse_at_wall_with_cx(cx, &p, 0.1, 1, wall(), budget(17, 1)),
    );
    assert!(matches!(
        exhausted,
        Err(VpmError::ExecutionBudget(
            BudgetRefusal::PollsExhausted { .. }
        ))
    ));
    with_cx(&CancelGate::new_clock_free(), Budget::INFINITE, |cx| {
        let retry = graph
            .diffuse_at_wall_with_cx(cx, &p, 0.1, 1, wall(), budget(17, 1))
            .unwrap();
        let dense = planar_wall::diffuse_with_cx(
            cx,
            &p,
            &[0.1; 17],
            0.1,
            1,
            graph.law(),
            wall(),
            budget(17, 1),
        )
        .unwrap();
        same_bits(&retry.particles, &dense.particles);
    });
    same_bits(&p, &saved);
}
