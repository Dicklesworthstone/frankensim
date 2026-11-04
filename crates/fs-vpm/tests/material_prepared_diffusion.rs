//! G0/G3: material-state changes reuse the original prepared PSE geometry.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_material::gas::{GasSpec, GasState};
use fs_vpm::diffusion::adaptive::{AdaptiveBudget, AdaptiveOptions, AdaptiveRun};
use fs_vpm::diffusion::{MaterialPseDiffusion, PreparationBudget, PreparedPse, PseDiffusion};
use fs_vpm::{VortexParticle, VpmBudget, VpmError, VpmRun};

const RADIUS: f64 = 0.003;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    let gate = CancelGate::new_clock_free();
    let clock = VirtualClock::new();
    let stream = StreamKey {
        seed: 0x5053_455f_5245_4249,
        kernel_id: 1,
        tile: 0,
        iteration: 0,
    };
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            stream,
            Budget::INFINITE,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn cloud() -> ([VortexParticle; 3], [f64; 3]) {
    (
        [
            VortexParticle::new([0.0, 0.0], 2e-6),
            VortexParticle::new([0.001, 0.0], -0.3e-6),
            VortexParticle::new([0.002, 0.0005], 0.8e-6),
        ],
        [1e-6, 0.8e-6, 1.2e-6],
    )
}

fn prepare(particles: &[VortexParticle], areas: &[f64], law: PseDiffusion) -> PreparedPse {
    with_cx(|cx| {
        PreparedPse::prepare_with_cx(
            cx,
            particles,
            areas,
            law,
            PreparationBudget::new(particles.len(), 16, 1000, 100_000, 16_384),
        )
    })
    .unwrap()
}

fn buffers(graph: &PreparedPse) -> (*const [f64; 2], *const f64, *const [usize; 2]) {
    (
        graph.positions().as_ptr(),
        graph.areas_m2().as_ptr(),
        graph.edges().as_ptr(),
    )
}

fn fixed(graph: &PreparedPse, particles: &[VortexParticle]) -> VpmRun {
    with_cx(|cx| {
        graph.diffuse_with_cx(
            cx,
            particles,
            0.02,
            5,
            VpmBudget::new(particles.len(), 5, 100_000, 16_384),
        )
    })
    .unwrap()
}

fn adaptive(graph: &PreparedPse, particles: &[VortexParticle]) -> AdaptiveRun {
    with_cx(|cx| {
        graph.advance_adaptive_with_cx(
            cx,
            particles,
            0.1,
            AdaptiveOptions {
                initial_dt_s: 0.02,
                min_dt_s: 1e-7,
                max_dt_s: 0.02,
                absolute_vorticity_tolerance_s_inv: 1e-6,
                relative_tolerance: 1e-5,
            },
            AdaptiveBudget {
                max_particles: particles.len(),
                max_attempts: 4096,
                max_pair_evaluations: 1_000_000,
                max_work_units: 2_000_000,
                max_live_bytes: 16_384,
            },
        )
    })
    .unwrap()
}

#[test]
fn g3_material_snapshots_reuse_graph_and_match_fresh_fixed_and_adaptive_runs() {
    let spec = GasSpec::dry_air_ussa1976();
    let states = [
        GasState::try_new(&spec, 293.15, 101_325.0).unwrap(),
        GasState::try_new(&spec, 373.15, 101_325.0).unwrap(),
        GasState::try_new(&spec, 293.15, 202_650.0).unwrap(),
    ];
    let original = MaterialPseDiffusion::compact_wendland(&states[0], RADIUS).unwrap();
    let (particles, areas) = cloud();
    let mut graph = prepare(&particles, &areas, original.law());
    assert_eq!(graph.edges().len(), 3);
    let addresses = buffers(&graph);
    let edges = graph.edges().to_vec();
    let report = *graph.preparation();
    let cold = fixed(&graph, &particles);

    for state in &states[1..] {
        let binding = MaterialPseDiffusion::compact_wendland(state, RADIUS).unwrap();
        graph.rebind_law(binding.law()).unwrap();
        assert_eq!(graph.law(), binding.law());
        assert_eq!(buffers(&graph), addresses);
        assert_eq!(graph.edges(), edges);
        assert_eq!(*graph.preparation(), report);
        assert_eq!(graph.positions(), particles.map(|p| p.pos));
        assert_eq!(graph.areas_m2(), areas);
        let fresh = prepare(&particles, &areas, binding.law());
        let observed = fixed(&graph, &particles);
        assert_eq!(observed, fixed(&fresh, &particles));
        assert_ne!(observed.particles, cold.particles);
        assert_eq!(adaptive(&graph, &particles), adaptive(&fresh, &particles));
    }

    // The scalar update uses the same checked numerical owner and restores the
    // first material snapshot without searching or retaining hot coefficients.
    graph
        .rebind_viscosity(original.transport().kinematic_viscosity_m2_s())
        .unwrap();
    assert_eq!(fixed(&graph, &particles), cold);
    assert_eq!(buffers(&graph), addresses);
    assert_eq!(*graph.preparation(), report);
}

#[test]
fn g0_zero_viscosity_keeps_edges_for_later_diffusion_and_returns_to_identity() {
    let (particles, areas) = cloud();
    let zero = PseDiffusion::compact_wendland(0.0, RADIUS).unwrap();
    let mut graph = prepare(&particles, &areas, zero);
    let addresses = buffers(&graph);
    let edges = graph.edges().to_vec();
    let report = *graph.preparation();
    assert_eq!(edges.len(), 3);
    let zero_run = fixed(&graph, &particles);
    assert_eq!(zero_run.particles, particles);
    assert_eq!(zero_run.pair_evaluations, 0);
    let zero_adaptive = adaptive(&graph, &particles);
    assert_eq!(zero_adaptive.particles, particles);
    assert_eq!(zero_adaptive.attempts, 0);
    assert_eq!(zero_adaptive.fraction_evaluations, 0);

    graph.rebind_viscosity(1.8e-5).unwrap();
    let fresh = prepare(&particles, &areas, graph.law());
    let active = fixed(&graph, &particles);
    assert_ne!(active.particles, particles);
    assert_eq!(active, fixed(&fresh, &particles));
    graph.rebind_viscosity(0.0).unwrap();
    assert_eq!(fixed(&graph, &particles), zero_run);
    assert_eq!(adaptive(&graph, &particles), zero_adaptive);
    assert_eq!(buffers(&graph), addresses);
    assert_eq!(graph.edges(), edges);
    assert_eq!(*graph.preparation(), report);
}

#[test]
fn g0_rebind_refuses_invalid_viscosity_or_changed_kernel_support_atomically() {
    let (particles, areas) = cloud();
    let original = PseDiffusion::compact_wendland(1.5e-5, RADIUS).unwrap();
    let mut graph = prepare(&particles, &areas, original);
    let addresses = buffers(&graph);
    let edges = graph.edges().to_vec();
    let report = *graph.preparation();
    let baseline = fixed(&graph, &particles);
    let unchanged = |graph: &PreparedPse| {
        assert_eq!(graph.law(), original);
        assert_eq!(buffers(graph), addresses);
        assert_eq!(graph.edges(), edges);
        assert_eq!(*graph.preparation(), report);
        assert_eq!(fixed(graph, &particles), baseline);
    };

    for viscosity in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let expected = PseDiffusion::compact_wendland(viscosity, RADIUS).unwrap_err();
        assert_eq!(graph.rebind_viscosity(viscosity), Err(expected));
        unchanged(&graph);
    }

    let gaussian = PseDiffusion::new(2e-5, RADIUS).unwrap();
    assert!(matches!(
        graph.rebind_law(gaussian),
        Err(VpmError::InvalidScalar {
            field: "prepared diffusion kernel",
            ..
        })
    ));
    unchanged(&graph);

    // Even a one-ULP radius change must not reuse edges discovered for another
    // exact support, irrespective of this particular cloud's nearest distances.
    for radius in [RADIUS.next_down(), RADIUS.next_up()] {
        let other = PseDiffusion::compact_wendland(2e-5, radius).unwrap();
        assert!(matches!(
            graph.rebind_law(other),
            Err(VpmError::InvalidScalar {
                field: "prepared diffusion support radius",
                ..
            })
        ));
        unchanged(&graph);
    }
}
