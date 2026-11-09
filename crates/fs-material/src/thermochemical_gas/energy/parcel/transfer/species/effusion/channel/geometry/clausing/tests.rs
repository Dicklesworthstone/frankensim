use super::*;

fn geometry(h: f64) -> CylindricalMolecularGeometry {
    CylindricalMolecularGeometry::try_new(
        1.0,
        h,
        100.0 * h.max(1.0),
        "independently qualified molecular flow",
    )
    .unwrap()
}

fn controls() -> ClausingControls {
    ClausingControls {
        initial_rings: 32,
        maximum_rings: 1024,
        maximum_iterations: 256,
        absolute_tolerance: 3e-5,
        relative_tolerance: 0.0,
    }
}

#[test]
fn resolved_transmission_matches_independent_integral_equation_results() {
    // Matuda, DOI 10.3131/jvsj2.58.311, Table 1: diameter-normalized L/D,
    // independently computed by Gomez-Goni/Lobo and Mohan et al.
    for (h, reference) in [
        (1.0, 0.514_230_535),
        (2.0, 0.356_572_732),
        (5.0, 0.190_942_432),
        (10.0, 0.109_320_714),
    ] {
        let original = geometry(h);
        let resolved = original
            .clone()
            .resolve_transmission(controls(), || false)
            .unwrap();
        let report = resolved.resolution().unwrap();
        assert!((resolved.transmission_probability() - reference).abs() < 1e-5);
        assert!(report.refinements >= 2);
        assert!(report.estimated_absolute_error <= controls().absolute_tolerance);
        assert_eq!(
            report.estimated_absolute_error,
            report.discretization_estimate
                + report.linear_solve_estimate
                + report.previous_linear_solve_estimate,
        );
        assert_eq!(
            report.linear_solve_estimate,
            4.0 * (h / report.rings as f64) * report.residual_l1,
        );
        assert_eq!(resolved.direct_probability(), original.direct_probability());
        assert_eq!(
            resolved.wall_hit_probability(),
            original.wall_hit_probability()
        );
        assert_eq!(
            resolved.minimum_mean_free_path_m(),
            original.minimum_mean_free_path_m()
        );
        assert_eq!(resolved.provenance(), original.provenance());
        assert!(
            (resolved.direct_probability()
                + resolved.wall_transmitted_probability()
                + resolved.reflected_probability()
                - 1.0)
                .abs()
                < 8.0 * f64::EPSILON,
        );
    }
}

#[test]
fn exact_ring_averages_conserve_probability_across_scales() {
    for delta in [1e-18, 0.002, 1.0, 1000.0] {
        let grid = Grid::new(64.0 * delta, 64, &mut || false).unwrap();
        for i in 0..64_usize {
            let total = compensated((0..64).map(|j| grid.wall[i.abs_diff(j)]))
                + grid.end[i]
                + grid.end[63 - i];
            assert!((total - 1.0).abs() < 8e-15, "delta={delta}, row={i}");
            assert!(grid.end[i] >= 0.0);
            assert!(grid.wall[i] >= 0.0);
        }
    }
    let tiny = Grid::new(512e-18, 512, &mut || false).unwrap();
    assert!((tiny.wall[500] / 1e-18 - 1.0).abs() < 1e-14);
    let distant = Grid::new(51.2, 512, &mut || false).unwrap();
    assert!((distant.wall[500] / 5.996_042_043_266_29e-9 - 1.0).abs() < 1e-13);
}

#[test]
fn nested_ring_approximations_converge_from_above() {
    let original = geometry(2.0);
    let values: Vec<_> = [32, 64, 128]
        .into_iter()
        .map(|rings| {
            let grid = Grid::new(2.0, rings, &mut || false).unwrap();
            let result = grid
                .solve(128, 1e-12, original.wall_hit_probability(), &mut || false)
                .unwrap();
            original.direct_probability() + result.wall_transmitted
        })
        .collect();
    assert!(values[0] > values[1] && values[1] > values[2]);
    assert!(((values[0] - values[1]) / (values[1] - values[2]) - 4.0).abs() < 0.01);
}

#[test]
fn zero_and_tiny_lengths_preserve_aperture_and_wall_heat_limits() {
    let zero = geometry(0.0)
        .resolve_transmission(controls(), || false)
        .unwrap();
    assert_eq!(zero.transmission_probability(), 1.0);
    assert_eq!(zero.reflected_probability(), 0.0);
    assert_eq!(zero.wall_transmitted_probability(), 0.0);
    assert_eq!(zero.resolution().unwrap().rings, 0);
    assert_eq!(zero.resolution().unwrap().estimated_absolute_error, 0.0);

    let tiny = geometry(1e-18)
        .resolve_transmission(controls(), || false)
        .unwrap();
    assert_eq!(tiny.transmission_probability(), 1.0);
    assert!((tiny.wall_transmitted_probability() / 1e-18 - 1.0).abs() < 1e-14);
    assert!((tiny.reflected_probability() / 1e-18 - 1.0).abs() < 1e-14);
    assert!(tiny.resolution().unwrap().refinements >= 2);
}

#[test]
fn resolution_is_deterministic_and_invariant_under_uniform_length_scaling() {
    let original = geometry(2.0);
    let first = original
        .clone()
        .resolve_transmission(controls(), || false)
        .unwrap();
    let repeated = original.resolve_transmission(controls(), || false).unwrap();
    assert_eq!(first, repeated);
    let scaled = CylindricalMolecularGeometry::try_new(0.125, 0.25, 25.0, "scaled qualification")
        .unwrap()
        .resolve_transmission(controls(), || false)
        .unwrap();
    assert_eq!(
        first.transmission_probability().to_bits(),
        scaled.transmission_probability().to_bits()
    );
    assert_eq!(
        first.reflected_probability().to_bits(),
        scaled.reflected_probability().to_bits()
    );
    assert_eq!(
        first.wall_transmitted_probability().to_bits(),
        scaled.wall_transmitted_probability().to_bits()
    );
    assert_eq!(first.resolution(), scaled.resolution());
}

#[test]
fn invalid_controls_and_exhausted_budgets_refuse_explicitly() {
    let default = controls();
    for bad in [
        ClausingControls {
            initial_rings: 0,
            ..default
        },
        ClausingControls {
            initial_rings: 5,
            ..default
        },
        ClausingControls {
            maximum_rings: 64,
            ..default
        },
        ClausingControls {
            maximum_rings: 2049,
            ..default
        },
        ClausingControls {
            maximum_iterations: 0,
            ..default
        },
        ClausingControls {
            maximum_iterations: 4097,
            ..default
        },
        ClausingControls {
            absolute_tolerance: -1.0,
            ..default
        },
        ClausingControls {
            absolute_tolerance: f64::NAN,
            ..default
        },
        ClausingControls {
            relative_tolerance: f64::INFINITY,
            ..default
        },
        ClausingControls {
            absolute_tolerance: 0.0,
            relative_tolerance: 0.0,
            ..default
        },
    ] {
        assert_eq!(
            geometry(1.0).resolve_transmission(bad, || false),
            Err(invalid("invalid Clausing accuracy or work budgets")),
        );
    }
    let spatial = ClausingControls {
        initial_rings: 4,
        maximum_rings: 16,
        absolute_tolerance: 1e-12,
        ..default
    };
    assert_eq!(
        geometry(10.0).resolve_transmission(spatial, || false),
        Err(invalid("Clausing spatial refinement budget exhausted")),
    );
    let linear = ClausingControls {
        maximum_iterations: 1,
        absolute_tolerance: 1e-12,
        ..default
    };
    assert_eq!(
        geometry(10.0).resolve_transmission(linear, || false),
        Err(invalid("Clausing linear solve budget exhausted")),
    );
}

#[test]
fn cancellation_is_observed_before_work_and_during_the_linear_solve() {
    let original = geometry(10.0);
    assert_eq!(
        original.clone().resolve_transmission(controls(), || true),
        Err(GasEnergyError::Cancelled),
    );
    let mut polls = 0;
    let result = original.clone().resolve_transmission(controls(), || {
        polls += 1;
        polls >= 12
    });
    assert_eq!(result, Err(GasEnergyError::Cancelled));
    assert_eq!(polls, 12);
    assert!(original.resolution().is_none());
}
