use fs_vpm::remesh::{CartesianGrid2, RemeshBudget, RemeshError, remesh};
use fs_vpm::{VortexParticle, induced_velocity};

fn budget() -> RemeshBudget {
    RemeshBudget {
        max_input_particles: 10000,
        max_grid_particles: 10000,
        max_work_units: 200000,
        max_live_bytes: 400000,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn node(grid: CartesianGrid2, x: f64, y: f64) -> [f64; 2] {
    let origin = grid.origin_m();
    let spacing = grid.spacing_m();
    [origin[0] + spacing[0] * x, origin[1] + spacing[1] * y]
}
fn moments(particles: &[VortexParticle], grid: CartesianGrid2) -> [f64; 6] {
    let mut sum = [0.0; 6];
    for p in particles {
        let x = (p.pos[0] - grid.origin_m()[0]) / grid.spacing_m()[0];
        let y = (p.pos[1] - grid.origin_m()[1]) / grid.spacing_m()[1];
        for (slot, monomial) in sum.iter_mut().zip([1.0, x, y, x * x, x * y, y * y]) {
            *slot += p.circulation * monomial;
        }
    }
    sum
}

#[test]
fn signed_m4_stencil_reproduces_every_quadratic_moment() {
    let grid = CartesianGrid2::new([-2.3, 4.1], [0.37, 0.29], [12, 11]).unwrap();
    let input = [VortexParticle::new(node(grid, 4.27, 5.62), -1.7)];
    let result = remesh(&input, grid, budget(), || false).unwrap();
    for (actual, expected) in moments(&result.particles, grid)
        .into_iter()
        .zip(moments(&input, grid))
    {
        near(actual, expected, 3e-13);
    }
    assert_eq!(
        result
            .particles
            .iter()
            .filter(|p| p.circulation != 0.0)
            .count(),
        16
    );
    assert!(result.particles.iter().any(|p| p.circulation > 0.0));
    assert!(result.particles.iter().any(|p| p.circulation < 0.0));
    assert_eq!(result.particles.len(), grid.node_count());
    assert_eq!(result.areas_m2, vec![0.37 * 0.29; grid.node_count()]);
    assert_eq!(result.stencil_deposits, 16);
    assert_eq!(result.work_units, 17 + 2 * grid.node_count() as u64);
    assert_eq!(
        result.peak_live_bytes,
        grid.node_count() * (std::mem::size_of::<VortexParticle>() + 8)
    );
}

#[test]
fn irregular_zero_net_cloud_retains_impulse_and_repeated_grid_identity() {
    let grid = CartesianGrid2::new([-1.0, -0.8], [0.13, 0.11], [20, 18]).unwrap();
    let mut input = Vec::new();
    for i in 0..17 {
        let k = f64::from(i);
        let gamma = 0.1 + 0.03 * k;
        input.push(VortexParticle::new(
            node(grid, 4.2 + (0.6 * k).sin(), 6.3 + (0.4 * k).cos()),
            gamma,
        ));
        input.push(VortexParticle::new(
            node(grid, 10.7 + (0.3 * k).cos(), 9.4 + (0.7 * k).sin()),
            -gamma,
        ));
    }
    let before = moments(&input, grid);
    assert_eq!(before[0], 0.0);
    let result = remesh(&input, grid, budget(), || false).unwrap();
    for (actual, expected) in moments(&result.particles, grid).into_iter().zip(before) {
        near(actual, expected, 4e-12);
    }
    let mut repeated = result.clone();
    for _ in 0..20 {
        repeated = remesh(&repeated.particles, grid, budget(), || false).unwrap();
        assert_eq!(repeated.particles, result.particles);
        assert_eq!(repeated.areas_m2, result.areas_m2);
    }
    let replay = remesh(&input, grid, budget(), || false).unwrap();
    for (a, b) in result.particles.iter().zip(replay.particles) {
        assert_eq!(a.circulation.to_bits(), b.circulation.to_bits());
        assert_eq!(a.pos.map(f64::to_bits), b.pos.map(f64::to_bits));
    }
}

#[test]
fn exact_decimal_grid_nodes_include_nonzero_boundary_particles() {
    let grid = CartesianGrid2::new([1.2, -3.7], [0.1, 0.03], [9, 8]).unwrap();
    let mut input = Vec::new();
    for y in 0..8 {
        for x in 0..9 {
            input.push(VortexParticle::new(
                node(grid, f64::from(x), f64::from(y)),
                f64::from(2 * x - y) / 7.0,
            ));
        }
    }
    let output = remesh(&input, grid, budget(), || false).unwrap();
    assert_eq!(output.particles, input);
}

#[test]
fn signed_deposition_retains_small_circulation_between_large_canceling_sources() {
    let grid = CartesianGrid2::new([0.0; 2], [1.0; 2], [6; 2]).unwrap();
    let p = node(grid, 2.0, 3.0);
    let input = [
        VortexParticle::new(p, 1e16),
        VortexParticle::new(p, 1.0),
        VortexParticle::new(p, -1e16),
    ];
    let output = remesh(&input, grid, budget(), || false).unwrap();
    assert_eq!(output.particles[3 * 6 + 2].circulation, 1.0);
    assert_eq!(
        output.particles.iter().map(|p| p.circulation).sum::<f64>(),
        1.0
    );
}

#[test]
fn translation_and_length_units_preserve_the_same_interpolation() {
    let grid = CartesianGrid2::new([-0.2, 0.3], [0.11, 0.17], [10, 9]).unwrap();
    let input = [
        VortexParticle::new(node(grid, 3.27, 4.16), 0.8),
        VortexParticle::new(node(grid, 5.43, 2.78), -1.2),
    ];
    let reference = remesh(&input, grid, budget(), || false).unwrap();
    for scale in [1e-80, 2.5, 1e80] {
        let shift = [1.2 * scale, -0.7 * scale];
        let target = CartesianGrid2::new(
            std::array::from_fn(|a| shift[a] + scale * grid.origin_m()[a]),
            grid.spacing_m().map(|h| h * scale),
            grid.shape(),
        )
        .unwrap();
        let transformed = input.map(|p| {
            VortexParticle::new(
                std::array::from_fn(|a| shift[a] + scale * p.pos[a]),
                p.circulation * scale * scale,
            )
        });
        let output = remesh(&transformed, target, budget(), || false).unwrap();
        for (a, b) in output.particles.iter().zip(&reference.particles) {
            near(a.circulation / scale / scale, b.circulation, 5e-14);
        }
        near(
            output.areas_m2[0] / scale / scale,
            reference.areas_m2[0],
            2e-16,
        );
    }
}

#[test]
fn induced_velocity_converges_to_the_original_cloud_when_grid_is_refined() {
    let source = VortexParticle::new([0.2, -0.1], 0.7);
    let probes = [[0.6, 0.3], [0.8, -0.25], [-0.1, 0.4]];
    let reference = probes.map(|p| induced_velocity(&[source], p, 0.3));
    let errors = [0.08, 0.04, 0.02].map(|h| {
        let grid = CartesianGrid2::new(
            [source.pos[0] - 4.27 * h, source.pos[1] - 4.16 * h],
            [h; 2],
            [10; 2],
        )
        .unwrap();
        let output = remesh(&[source], grid, budget(), || false).unwrap();
        probes
            .iter()
            .zip(reference)
            .map(|(&p, expected)| {
                let actual = induced_velocity(&output.particles, p, 0.3);
                (actual[0] - expected[0]).powi(2) + (actual[1] - expected[1]).powi(2)
            })
            .sum::<f64>()
            .sqrt()
    });
    assert!(errors[2] < 2e-6, "errors={errors:?}");
    for pair in errors.windows(2) {
        assert!(pair[0] > 6.0 * pair[1], "errors={errors:?}");
    }
}

#[test]
fn refuses_nonzero_escaping_support_and_keeps_zero_strength_input_harmless() {
    let grid = CartesianGrid2::new([0.0; 2], [1.0; 2], [6; 2]).unwrap();
    for x in [0.2, 4.7, -0.1] {
        let input = [VortexParticle::new([x, 2.7], 1.0)];
        assert!(matches!(
            remesh(&input, grid, budget(), || false),
            Err(RemeshError::SupportOutsideGrid {
                particle: 0,
                axis: 0
            })
        ));
        assert_eq!(input[0].pos, [x, 2.7]);
    }
    let output = remesh(
        &[VortexParticle::new([1e9; 2], 0.0)],
        grid,
        budget(),
        || false,
    )
    .unwrap();
    assert!(output.particles.iter().all(|p| p.circulation == 0.0));
    assert!(
        remesh(
            &[VortexParticle::new([f64::INFINITY, 0.0], 0.0)],
            grid,
            budget(),
            || false
        )
        .is_err()
    );
    assert!(CartesianGrid2::new([0.0; 2], [1e-200; 2], [6; 2]).is_err());
    assert!(CartesianGrid2::new([1e16, 0.0], [1.0; 2], [6; 2]).is_err());
    // This source has a finite nonzero first moment, but normalizing its x
    // coordinate underflows. It must not be snapped to the boundary node.
    let anisotropic = CartesianGrid2::new([0.0; 2], [1e150, 1e-150], [6; 2]).unwrap();
    assert!(matches!(
        remesh(
            &[VortexParticle::new([1e-200, 2.7e-150], 1e200)],
            anisotropic,
            budget(),
            || false
        ),
        Err(RemeshError::Invalid(_))
    ));
}

#[test]
fn explicit_resource_limits_and_every_cancellation_boundary_are_atomic() {
    let grid = CartesianGrid2::new([0.0; 2], [0.1; 2], [10; 2]).unwrap();
    let input = [VortexParticle::new(node(grid, 4.27, 3.16), 0.8); 20];
    let original = input;
    let mut polls = 0;
    let expected = remesh(&input, grid, budget(), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in 1..=polls {
        let mut visits = 0;
        assert!(matches!(
            remesh(&input, grid, budget(), || {
                visits += 1;
                visits == stop
            }),
            Err(RemeshError::Cancelled)
        ));
        assert_eq!(input, original);
    }
    let mut caps = [budget(); 4];
    caps[0].max_input_particles = input.len() - 1;
    caps[1].max_grid_particles = grid.node_count() - 1;
    caps[2].max_work_units = expected.work_units - 1;
    caps[3].max_live_bytes = expected.peak_live_bytes - 1;
    for cap in caps {
        assert!(matches!(
            remesh(&input, grid, cap, || false),
            Err(RemeshError::Limit { .. })
        ));
    }
    assert_eq!(remesh(&input, grid, budget(), || false).unwrap(), expected);
}

#[test]
fn remeshed_areas_feed_the_original_viscous_pse_owner_without_reconstruction() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_vpm::{
        VpmBudget,
        diffusion::{PseDiffusion, simulate_viscous_with_cx},
    };
    let grid = CartesianGrid2::new([-1.2; 2], [0.2; 2], [13; 2]).unwrap();
    let initial = [
        VortexParticle::new([-0.21, 0.07], 0.2),
        VortexParticle::new([0.31, -0.13], -0.2),
    ];
    let cloud = remesh(&initial, grid, budget(), || false).unwrap();
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 4,
                kernel_id: 9,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let run = simulate_viscous_with_cx(
            &cx,
            &cloud.particles,
            &cloud.areas_m2,
            0.002,
            2,
            0.2,
            PseDiffusion::new(0.01, 0.4).unwrap(),
            VpmBudget::new(cloud.particles.len(), 2, 1000000, 1000000),
        )
        .unwrap();
        near(
            run.particles.iter().map(|p| p.circulation).sum(),
            0.0,
            2e-15,
        );
        assert!(
            run.particles
                .iter()
                .all(|p| p.circulation.is_finite() && p.pos.iter().all(|x| x.is_finite()))
        );
        assert!(
            run.particles
                .iter()
                .zip(&cloud.particles)
                .any(|(a, b)| a.circulation != b.circulation)
        );
        // Subsequent finite-cloud PSE has its own spatial error; remeshing
        // does not turn its circulation invariant into exact impulse closure.
        // The moving finite cloud gets explicit padding for its next stencil;
        // even small diffused tails are retained rather than truncated.
        let padded = CartesianGrid2::new([-1.6; 2], [0.2; 2], [17; 2]).unwrap();
        let again = remesh(&run.particles, padded, budget(), || false).unwrap();
        for (a, b) in moments(&again.particles, grid)
            .into_iter()
            .zip(moments(&run.particles, grid))
        {
            near(a, b, 2e-12);
        }
    });
}
