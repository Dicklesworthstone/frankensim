use super::*;

mod batch;

fn geometry(h: f64) -> CylindricalMolecularGeometry {
    CylindricalMolecularGeometry::try_new(1.0, h, 100.0 * h.max(1.0), "test collision support")
        .unwrap()
}

fn controls() -> ClausingControls {
    ClausingControls {
        initial_rings: 8,
        maximum_rings: 512,
        maximum_iterations: 256,
        absolute_tolerance: 1e-6,
        relative_tolerance: 1e-4,
    }
}

#[test]
fn exact_accommodation_and_aperture_limits_preserve_mass_geometry() {
    let original = geometry(2.0);
    for alpha in [0.0, 1.0] {
        let (resolved, response) = original
            .clone()
            .resolve_thermal_accommodation(alpha, controls(), || false)
            .unwrap();
        assert_eq!(resolved, original);
        assert_eq!(response.coefficient(), alpha);
        assert_eq!(
            response.transmitted_heat_probability(),
            alpha * original.wall_transmitted_probability()
        );
        assert_eq!(
            response.reflected_heat_probability(),
            alpha * original.reflected_probability()
        );
        assert!(response.resolution().is_none());
    }
    for alpha in [0.0, 0.4, 1.0] {
        let (resolved, response) = geometry(0.0)
            .resolve_thermal_accommodation(alpha, controls(), || false)
            .unwrap();
        assert_eq!(resolved.transmission_probability(), 1.0);
        assert_eq!(response.transmitted_heat_probability(), 0.0);
        assert_eq!(response.reflected_heat_probability(), 0.0);
        assert!(response.resolution().is_none());
    }
}

#[test]
fn collision_count_enumeration_independently_checks_repeated_accommodation() {
    let original = geometry(2.0);
    let grid = Grid::new(2.0, 64, &mut || false).unwrap();
    let mass = grid
        .solve(256, 1e-13, original.wall_hit_probability(), &mut || false)
        .unwrap();
    for alpha in [0.2, 0.7] {
        let resolved = solve_heat(&grid, &mass, alpha, 256, 1e-13, &mut || false).unwrap();
        // Independent collision-count transport: move the original incident
        // population one diffuse wall flight at a time and accumulate 1-q^N.
        let mut population: Vec<_> = grid.end.iter().map(|e| 4.0 * grid.delta * e).collect();
        let mut thermalized = 0.0;
        let (mut transmitted, mut reflected) = (0.0, 0.0);
        for _ in 0..2000 {
            thermalized += (1.0 - thermalized) * alpha;
            transmitted +=
                thermalized * compensated((0..64).map(|i| population[i] * grid.end[63 - i]));
            reflected += thermalized * compensated((0..64).map(|i| population[i] * grid.end[i]));
            population = (0..64_usize)
                .map(|j| {
                    compensated((0..64_usize).map(|i| population[i] * grid.wall[i.abs_diff(j)]))
                })
                .collect();
            if compensated(population.iter().copied()) < 1e-13 {
                break;
            }
        }
        assert!(compensated(population.iter().copied()) < 1e-13);
        assert!((resolved.transmitted - transmitted).abs() < 2e-11);
        assert!((resolved.reflected - reflected).abs() < 2e-11);
        assert!(resolved.transmitted > alpha * mass.wall_transmitted * 1.1);
        assert!(resolved.reflected > alpha * mass.reflected * 1.1);
        assert!(resolved.transmitted < mass.wall_transmitted);
        assert!(resolved.reflected < mass.reflected);
    }
}

#[test]
fn tiny_coefficients_retain_the_collision_count_first_moment() {
    let original = geometry(1.0);
    let grid = Grid::new(1.0, 64, &mut || false).unwrap();
    let mass = grid
        .solve(256, 1e-13, original.wall_hit_probability(), &mut || false)
        .unwrap();
    let alpha = 1e-18;
    let heat = solve_heat(&grid, &mass, alpha, 256, 1e-13, &mut || false).unwrap();
    // dH/dalpha at zero is E[N 1(exit)]. Reciprocity gives the following
    // independent first-moment integrals from the ordinary escape field.
    let transmitted_moment =
        8.0 * grid.delta * compensated(mass.half_escape.iter().map(|u| 0.25 - u * u));
    let reflected_moment =
        4.0 * grid.delta * compensated(mass.half_escape.iter().map(|u| 0.5 + 2.0 * u * u));
    assert!((heat.transmitted / alpha / transmitted_moment - 1.0).abs() < 1e-10);
    assert!((heat.reflected / alpha / reflected_moment - 1.0).abs() < 1e-10);
    // At equilibrium, wall-collision flux / combined entrance flux gives
    // E[N] = wall_area/(2*mouth_area) = 2 L/D, independent of the exit.
    assert!(((heat.transmitted + heat.reflected) / alpha - 2.0).abs() < 1e-10);
    let (_, thin) = geometry(1e-18)
        .resolve_thermal_accommodation(0.25, controls(), || false)
        .unwrap();
    assert!((thin.transmitted_heat_probability() / 2.5e-19 - 1.0).abs() < 1e-13);
    assert!((thin.reflected_heat_probability() / 2.5e-19 - 1.0).abs() < 1e-13);
}

#[test]
fn accommodation_is_monotone_and_both_numerical_branches_agree() {
    let original = geometry(1.0);
    let grid = Grid::new(1.0, 64, &mut || false).unwrap();
    let mass = grid
        .solve(256, 1e-13, original.wall_hit_probability(), &mut || false)
        .unwrap();
    let mut previous = [0.0; 2];
    for alpha in [
        0.001,
        0.1,
        0.499_999_99,
        0.500_000_01,
        0.9,
        1.0 - f64::EPSILON,
    ] {
        let current = solve_heat(&grid, &mass, alpha, 256, 1e-13, &mut || false).unwrap();
        assert!(current.transmitted >= previous[0]);
        assert!(current.reflected >= previous[1]);
        if alpha == 0.500_000_01 {
            assert!(current.transmitted - previous[0] < 1e-7);
            assert!(current.reflected - previous[1] < 1e-7);
        }
        previous = [current.transmitted, current.reflected];
    }
    assert!((previous[0] / mass.wall_transmitted - 1.0).abs() < 1e-13);
    assert!((previous[1] / mass.reflected - 1.0).abs() < 1e-13);
}

#[test]
fn public_resolution_retains_three_grid_evidence_and_is_scale_invariant() {
    let original = geometry(1.0);
    let (resolved, response) = original
        .clone()
        .resolve_thermal_accommodation(0.3, controls(), || false)
        .unwrap();
    let repeated = original
        .resolve_thermal_accommodation(0.3, controls(), || false)
        .unwrap();
    assert_eq!((&resolved, &response), (&repeated.0, &repeated.1));
    let scaled = CylindricalMolecularGeometry::try_new(0.125, 0.125, 12.5, "scaled test support")
        .unwrap()
        .resolve_thermal_accommodation(0.3, controls(), || false)
        .unwrap();
    assert_eq!(response, scaled.1);
    assert_eq!(
        resolved.transmission_probability(),
        scaled.0.transmission_probability()
    );
    let report = response.resolution().unwrap();
    assert!(report.refinements >= 2);
    assert_eq!(report.rings, resolved.resolution().unwrap().rings);
    assert_eq!(
        report.transmitted_error_estimate,
        report.transmitted_grid_change
            + report.linear_solve_estimate
            + report.previous_linear_solve_estimate
    );
    assert_eq!(
        report.reflected_error_estimate,
        report.reflected_grid_change
            + report.linear_solve_estimate
            + report.previous_linear_solve_estimate
    );
    assert!(
        report.transmitted_error_estimate
            <= controls().absolute_tolerance
                + controls().relative_tolerance * response.transmitted_heat_probability()
    );
    assert!(
        report.reflected_error_estimate
            <= controls().absolute_tolerance
                + controls().relative_tolerance * response.reflected_heat_probability()
    );
}

#[test]
fn invalid_coefficients_exhausted_budgets_and_cancellation_refuse() {
    for coefficient in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert_eq!(
            geometry(1.0).resolve_thermal_accommodation(coefficient, controls(), || false),
            Err(invalid(
                "thermal accommodation coefficient must be finite and in [0,1]"
            ))
        );
    }
    let exhausted = ClausingControls {
        initial_rings: 4,
        maximum_rings: 16,
        absolute_tolerance: 1e-12,
        relative_tolerance: 0.0,
        ..controls()
    };
    assert_eq!(
        geometry(2.0).resolve_thermal_accommodation(0.3, exhausted, || false),
        Err(invalid("accommodation spatial refinement budget exhausted"))
    );
    assert_eq!(
        geometry(1.0).resolve_thermal_accommodation(0.3, controls(), || true),
        Err(GasEnergyError::Cancelled)
    );
    let grid = Grid::new(1.0, 64, &mut || false).unwrap();
    let mass = grid
        .solve(
            256,
            1e-12,
            geometry(1.0).wall_hit_probability(),
            &mut || false,
        )
        .unwrap();
    let mut polls = 0;
    let cancelled = solve_heat(&grid, &mass, 0.3, 256, 1e-12, &mut || {
        polls += 1;
        polls >= 3
    });
    assert!(matches!(cancelled, Err(GasEnergyError::Cancelled)));
    assert_eq!(polls, 3);
}
