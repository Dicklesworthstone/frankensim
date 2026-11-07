//! Periodic heat-equation and exact pair oracles for compact vorticity diffusion.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_vpm::diffusion::adaptive::{AdaptiveBudget, AdaptiveOptions};
use fs_vpm::diffusion::{PeriodicBox2, PreparationBudget, PreparedPse, PseDiffusion};
use fs_vpm::planar_wall::PlanarWall2;
use fs_vpm::{VortexParticle, VpmBudget, VpmError, total_circulation};
use std::f64::consts::PI;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    let gate = CancelGate::new_clock_free();
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x544f_5255_53,
                kernel_id: 0x5053_45,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn prep_budget(n: usize) -> PreparationBudget {
    PreparationBudget::new(n, usize::MAX, 50_000_000, 50_000_000, usize::MAX)
}

fn run_budget(n: usize, steps: usize) -> VpmBudget {
    VpmBudget::new(n, steps, u64::MAX, usize::MAX)
}

fn prepare(cx: &Cx<'_>, p: &[VortexParticle], areas: &[f64], nu: f64, radius: f64) -> PreparedPse {
    PreparedPse::prepare_periodic_with_cx(
        cx,
        p,
        areas,
        PseDiffusion::compact_wendland(nu, radius).unwrap(),
        PeriodicBox2::new([1.0; 2]).unwrap(),
        prep_budget(p.len()),
    )
    .unwrap()
}

fn exact_pair(
    initial: [f64; 2],
    areas: [f64; 2],
    nu: f64,
    radius: f64,
    distance: f64,
    time: f64,
) -> [f64; 2] {
    let q = distance / radius;
    let rate = 1008.0 * nu / (5.0 * PI * radius.powi(4)) * (1.0 - q).powi(4) * (1.0 + 4.0 * q);
    let total = initial[0] + initial[1];
    let equilibrium = areas[0] / (areas[0] + areas[1]) * total;
    let first =
        equilibrium + (initial[0] - equilibrium) * (-rate * (areas[0] + areas[1]) * time).exp();
    [first, total - first]
}

#[test]
fn g1_corner_seam_relaxes_at_physical_rate_for_unequal_particle_areas() {
    let initial = [
        VortexParticle::new([0.02, 0.03], 0.6),
        VortexParticle::new([0.98, 0.97], -0.1),
    ];
    let areas = [0.2, 0.7];
    with_cx(|cx| {
        let graph = prepare(cx, &initial, &areas, 0.08, 0.2);
        assert_eq!(graph.edges(), &[[0, 1]]);
        let finite =
            PreparedPse::prepare_with_cx(cx, &initial, &areas, graph.law(), prep_budget(2))
                .unwrap();
        assert!(finite.edges().is_empty());
        assert!(finite.periodic_domain().is_none());
        let expected = exact_pair([0.6, -0.1], areas, 0.08, 0.2, 0.04_f64.hypot(0.06), 0.001);
        for steps in [1, 7] {
            let run = graph
                .diffuse_with_cx(
                    cx,
                    &initial,
                    0.001 / steps as f64,
                    steps,
                    run_budget(2, steps),
                )
                .unwrap();
            for (particle, exact) in run.particles.iter().zip(expected) {
                assert!((particle.circulation - exact).abs() < 2e-14);
            }
            assert!((total_circulation(&run.particles) - 0.5).abs() < 2e-15);
            assert_eq!(run.pair_evaluations, 2 * steps as u64);
            assert_eq!(run.work_units, 3 + 3 * steps as u64);
            assert_eq!(run.budget.cost_charged, run.work_units);
        }
        let adaptive = graph
            .advance_adaptive_with_cx(
                cx,
                &initial,
                0.001,
                AdaptiveOptions {
                    initial_dt_s: 0.001,
                    min_dt_s: 1e-9,
                    max_dt_s: 0.001,
                    absolute_vorticity_tolerance_s_inv: 1e-10,
                    relative_tolerance: 1e-10,
                },
                AdaptiveBudget {
                    max_particles: 2,
                    max_attempts: 100,
                    max_pair_evaluations: 1000,
                    max_work_units: 10_000,
                    max_live_bytes: 1024 * 1024,
                },
            )
            .unwrap();
        assert_eq!(adaptive.accepted_steps, 1);
        for (particle, exact) in adaptive.particles.iter().zip(expected) {
            assert!((particle.circulation - exact).abs() < 2e-14);
        }
    });
}

#[test]
fn g0_periodic_graph_matches_exhaustive_images_without_duplicate_small_grid_edges() {
    // Rectangular, incommensurate periods and 2/3/4/5/10 cells per x-period.
    // Repeated halo residues must not duplicate exchanges on short periods.
    let lengths = [1.0, 1.37];
    let mut particles: Vec<_> = (0..37)
        .map(|i| {
            VortexParticle::new(
                [
                    ((i * 19) % 37) as f64 / 37.0,
                    1.37 * ((i * 13) % 37) as f64 / 37.0,
                ],
                0.1 * i as f64 - 1.0,
            )
        })
        .collect();
    particles.push(VortexParticle::new(
        [1.0_f64.next_down(), 1.37_f64.next_down()],
        0.5,
    ));
    particles.push(VortexParticle::new([0.0, 0.0], -0.5));
    let areas = vec![0.02; particles.len()];
    with_cx(|cx| {
        for radius in [0.5, 0.31, 0.24, 0.19, 0.097] {
            let graph = PreparedPse::prepare_periodic_with_cx(
                cx,
                &particles,
                &areas,
                PseDiffusion::compact_wendland(0.1, radius).unwrap(),
                PeriodicBox2::new(lengths).unwrap(),
                prep_budget(particles.len()),
            )
            .unwrap();
            let mut expected = Vec::new();
            for i in 0..particles.len() {
                for j in i + 1..particles.len() {
                    let mut inside = false;
                    for ix in -1..=1 {
                        for iy in -1..=1 {
                            let dx = particles[i].pos[0] - particles[j].pos[0]
                                + f64::from(ix) * lengths[0];
                            let dy = particles[i].pos[1] - particles[j].pos[1]
                                + f64::from(iy) * lengths[1];
                            inside |= dx.hypot(dy) < radius;
                        }
                    }
                    if inside {
                        expected.push([i, j]);
                    }
                }
            }
            assert_eq!(graph.edges(), expected, "radius={radius}");
            assert_eq!(
                graph.preparation().work_units,
                graph.preparation().budget.cost_charged
            );
            assert!(
                graph.preparation().candidate_visits
                    <= (particles.len() * (particles.len() - 1)) as u64
            );
        }
    });
}

#[test]
fn g1_stiff_periodic_diffusion_conserves_signed_circulation_and_dissipates_enstrophy() {
    let initial: Vec<_> = (0..16)
        .map(|i| {
            VortexParticle::new(
                [(i % 4) as f64 / 4.0, (i / 4) as f64 / 4.0],
                0.2 * (i as f64 - 4.0),
            )
        })
        .collect();
    let areas: Vec<_> = (0..16).map(|i| 0.01 * (i + 1) as f64).collect();
    let enstrophy = |p: &[VortexParticle]| -> f64 {
        p.iter()
            .zip(&areas)
            .map(|(p, area)| p.circulation.powi(2) / area)
            .sum()
    };
    let old_enstrophy = enstrophy(&initial);
    let minimum = initial
        .iter()
        .zip(&areas)
        .map(|(p, a)| p.circulation / a)
        .fold(f64::INFINITY, f64::min);
    let maximum = initial
        .iter()
        .zip(&areas)
        .map(|(p, a)| p.circulation / a)
        .fold(f64::NEG_INFINITY, f64::max);
    with_cx(|cx| {
        let mut graph = prepare(cx, &initial, &areas, 1e5, 0.45);
        let run = graph
            .diffuse_with_cx(cx, &initial, 1.0, 20, run_budget(16, 20))
            .unwrap();
        assert!((total_circulation(&run.particles) - total_circulation(&initial)).abs() < 2e-13);
        assert!(enstrophy(&run.particles) < old_enstrophy);
        for (p, area) in run.particles.iter().zip(&areas) {
            assert!((minimum..=maximum).contains(&(p.circulation / area)));
        }
        let uniform: Vec<_> = initial
            .iter()
            .zip(&areas)
            .map(|(p, a)| VortexParticle::new(p.pos, 3.0 * a))
            .collect();
        graph.rebind_viscosity(0.3).unwrap();
        assert_eq!(graph.periodic_domain().unwrap().lengths_m(), [1.0; 2]);
        let unchanged = graph
            .diffuse_with_cx(cx, &uniform, 0.01, 10, run_budget(16, 10))
            .unwrap();
        for (p, a) in unchanged.particles.iter().zip(&areas) {
            assert!((p.circulation / a - 3.0).abs() < 3e-15);
        }
    });
}

fn fourier_errors(n: usize, radius: f64) -> (f64, f64) {
    let nu = 0.1;
    let time = 0.1;
    let area = 1.0 / (n * n) as f64;
    let particles: Vec<_> = (0..n * n)
        .map(|i| {
            let pos = [(i % n) as f64 / n as f64, (i / n) as f64 / n as f64];
            VortexParticle::new(pos, area * (2.0 * PI * pos[0]).cos())
        })
        .collect();
    let areas = vec![area; particles.len()];
    // Independent Fourier symbol of the periodized nonlocal discrete operator.
    // This separates temporal composition error from compact-kernel consistency.
    let mut eigenvalue = 0.0;
    for p in &particles[1..] {
        let distance = p.pos[0]
            .min(1.0 - p.pos[0])
            .hypot(p.pos[1].min(1.0 - p.pos[1]));
        let q = distance / radius;
        if q < 1.0 {
            let rate =
                1008.0 * nu / (5.0 * PI * radius.powi(4)) * (1.0 - q).powi(4) * (1.0 + 4.0 * q);
            eigenvalue += area * rate * ((2.0 * PI * p.pos[0]).cos() - 1.0);
        }
    }
    let heat_amplitude = (-4.0 * PI * PI * nu * time).exp();
    let discrete_amplitude = (eigenvalue * time).exp();
    with_cx(|cx| {
        let graph = prepare(cx, &particles, &areas, nu, radius);
        let run = graph
            .diffuse_with_cx(
                cx,
                &particles,
                time / 100.0,
                100,
                run_budget(particles.len(), 100),
            )
            .unwrap();
        let error = |amplitude: f64| {
            (run.particles
                .iter()
                .map(|p| (p.circulation / area - amplitude * (2.0 * PI * p.pos[0]).cos()).powi(2))
                .sum::<f64>()
                * area)
                .sqrt()
        };
        assert!(total_circulation(&run.particles).abs() < 2e-14);
        (error(heat_amplitude), error(discrete_amplitude))
    })
}

#[test]
fn g1_periodic_fourier_mode_converges_to_continuum_viscous_decay() {
    let coarse = fourier_errors(8, 0.45);
    let fine = fourier_errors(16, 0.225);
    assert!(fine.0 < 0.007, "continuum errors {coarse:?} -> {fine:?}");
    assert!(
        coarse.0 > 3.0 * fine.0,
        "continuum errors {coarse:?} -> {fine:?}"
    );
    assert!(
        coarse.1 < 1e-5 && fine.1 < 1e-5,
        "temporal errors {coarse:?}, {fine:?}"
    );
}

#[test]
fn g4_periodic_geometry_refuses_ambiguous_images_and_incompatible_walls() {
    for lengths in [
        [0.0, 1.0],
        [-1.0, 1.0],
        [f64::INFINITY, 1.0],
        [1.0, f64::NAN],
    ] {
        assert!(PeriodicBox2::new(lengths).is_err());
    }
    let initial = [
        VortexParticle::new([0.05, 0.05], 1.0),
        VortexParticle::new([0.95, 0.95], 0.0),
    ];
    let areas = [0.1; 2];
    let domain = PeriodicBox2::new([1.0; 2]).unwrap();
    with_cx(|cx| {
        for radius in [0.51, 1e-16] {
            assert!(
                PreparedPse::prepare_periodic_with_cx(
                    cx,
                    &initial,
                    &areas,
                    PseDiffusion::compact_wendland(0.1, radius).unwrap(),
                    domain,
                    prep_budget(2),
                )
                .is_err()
            );
        }
        for outside in [[1.0, 0.0], [-f64::MIN_POSITIVE, 0.1]] {
            let p = [VortexParticle::new(outside, 1.0)];
            assert!(
                PreparedPse::prepare_periodic_with_cx(
                    cx,
                    &p,
                    &[0.1],
                    PseDiffusion::compact_wendland(0.1, 0.2).unwrap(),
                    domain,
                    prep_budget(1),
                )
                .is_err()
            );
        }
        let graph = prepare(cx, &initial, &areas, 0.1, 0.2);
        let too_few_edges = PreparationBudget {
            max_edges: 0,
            ..prep_budget(2)
        };
        assert!(matches!(
            PreparedPse::prepare_periodic_with_cx(
                cx,
                &initial,
                &areas,
                graph.law(),
                domain,
                too_few_edges,
            ),
            Err(VpmError::CapExceeded {
                resource: "prepared edges",
                ..
            })
        ));
        assert!(matches!(
            graph.diffuse_at_wall_with_cx(
                cx,
                &initial,
                0.01,
                1,
                PlanarWall2::new([0.0; 2], [0.0, 1.0]).unwrap(),
                run_budget(2, 1),
            ),
            Err(VpmError::InvalidScalar {
                field: "planar wall on periodic diffusion graph",
                ..
            })
        ));
        // Failed use does not poison the immutable graph or mutate the input.
        let run = graph
            .diffuse_with_cx(cx, &initial, 0.01, 1, run_budget(2, 1))
            .unwrap();
        assert!(run.particles[1].circulation > 0.0);
        assert_eq!(initial[1].circulation.to_bits(), 0.0_f64.to_bits());
    });
}
