//! G0/G1/G3/G4: joint angle/history, independent collision enumeration and rays.
use super::*;

fn split(response: &MolecularAxialAccommodation) -> [f64; 4] {
    [
        response.retained_transmitted_axial_cosine(),
        response.reset_transmitted_axial_cosine(),
        response.retained_reflected_axial_cosine(),
        response.reset_reflected_axial_cosine(),
    ]
}

#[test]
fn exact_limits_share_mass_heat_and_total_angular_geometry() {
    for h in [0.0, 0.01, 2.0] {
        let coefficients = [1.0, 0.0, 0.25, 0.75, 0.25];
        let (g, total, rows) = geometry(h)
            .resolve_axial_accommodation_batch(&coefficients, controls(), || false)
            .unwrap();
        assert_eq!(rows[2], rows[4]);
        assert_eq!(
            rows[0].accommodation().transmitted_heat_probability(),
            g.wall_transmitted_probability()
        );
        assert_eq!(
            rows[0].accommodation().reflected_heat_probability(),
            g.reflected_probability()
        );
        assert_eq!(
            split(&rows[0]),
            [
                total.direct_axial_cosine(),
                total.wall_transmitted_axial_cosine(),
                0.0,
                total.reflected_axial_cosine()
            ]
        );
        assert_eq!(
            split(&rows[1]),
            [
                total.transmitted_axial_cosine(),
                0.0,
                total.reflected_axial_cosine(),
                0.0
            ]
        );
        for (&coefficient, row) in coefficients.iter().zip(&rows) {
            assert_eq!(row.coefficient(), coefficient);
            let [retained_t, reset_t, retained_r, reset_r] = split(row);
            near(
                retained_t + reset_t,
                total.transmitted_axial_cosine(),
                2e-16,
            );
            near(retained_r + reset_r, total.reflected_axial_cosine(), 2e-16);
            assert!(reset_t <= row.accommodation().transmitted_heat_probability());
            assert!(reset_r <= row.accommodation().reflected_heat_probability());
            if let Some(report) = row.resolution() {
                assert_eq!(report.rings, g.resolution().unwrap().rings);
                assert_eq!(
                    report.rings,
                    row.accommodation().resolution().unwrap().rings
                );
                assert_eq!(report.angular_panels, total.resolution().angular_panels);
                assert_eq!(report.batch_size, coefficients.len());
                assert_eq!(report.thermal_coefficients, 2);
                assert_eq!(
                    report.iterations - report.angular_escape_iterations,
                    row.accommodation().resolution().unwrap().iterations
                );
                for (value, error) in [retained_t, retained_r]
                    .into_iter()
                    .zip(report.retained_error_estimate)
                {
                    assert!(
                        error
                            <= controls().absolute_tolerance
                                + controls().relative_tolerance * value
                    );
                }
                for (value, error) in [reset_t, reset_r]
                    .into_iter()
                    .zip(report.reset_error_estimate)
                {
                    assert!(
                        error
                            <= controls().absolute_tolerance
                                + controls().relative_tolerance * value
                    );
                }
            } else {
                assert!(h == 0.0 || coefficient == 0.0 || coefficient == 1.0);
            }
        }
    }
}

// Independent Neumann-series flight propagation, with no inverse, CG solve or
// resolvent identity. Count N on each wall-emission generation and accumulate
// q^N and 1-q^N against the outgoing angle and ordinary exit kernels separately.
fn collision_series(
    grid: &Grid,
    angular: &AngularIntegrals,
    wall_hit: f64,
    alpha: f64,
) -> [f64; 6] {
    let n = grid.end.len();
    let scale = wall_hit / compensated(grid.end.iter().copied());
    let mut population: Vec<_> = grid.end.iter().map(|e| scale * e).collect();
    let mut result = [angular.values[0], 0.0, 0.0, 0.0, 0.0, 0.0];
    let (mut retained, mut reset) = (1.0, 0.0);
    let mut finished = false;
    for _ in 0..3000 {
        retained *= 1.0 - alpha;
        reset += alpha * (1.0 - reset);
        for side in 0..2 {
            let moment = compensated(
                (0..n).map(|i| population[i] * angular.end[if side == 0 { n - 1 - i } else { i }]),
            );
            result[2 * side] += retained * moment;
            result[2 * side + 1] += reset * moment;
            result[4 + side] += reset
                * compensated(
                    (0..n).map(|i| population[i] * grid.end[if side == 0 { n - 1 - i } else { i }]),
                );
        }
        population = (0..n)
            .map(|j| compensated((0..n).map(|i| population[i] * grid.wall[i.abs_diff(j)])))
            .collect();
        if compensated(population.iter().copied()) < 1e-16 {
            finished = true;
            break;
        }
    }
    assert!(
        finished,
        "independent collision series exceeded its explicit work bound"
    );
    result
}

#[test]
fn joint_moments_and_heat_match_independent_collision_count_series() {
    let coefficients = [1e-18, 0.2, 0.7, 1.0 - f64::EPSILON];
    let (g, total, rows) = geometry(2.0)
        .resolve_axial_accommodation_batch(&coefficients, controls(), || false)
        .unwrap();
    let grid = Grid::new(2.0, g.resolution().unwrap().rings, &mut || false).unwrap();
    let mass = grid
        .solve(512, 1e-13, g.wall_hit_probability(), &mut || false)
        .unwrap();
    let angular = angular_moments(
        &grid,
        &mass,
        2.0,
        total.resolution().angular_panels,
        &mut || false,
    )
    .unwrap();
    for (&alpha, row) in coefficients.iter().zip(&rows) {
        let reference = collision_series(&grid, &angular, g.wall_hit_probability(), alpha);
        let actual = [
            split(row)[0],
            split(row)[1],
            split(row)[2],
            split(row)[3],
            row.accommodation().transmitted_heat_probability(),
            row.accommodation().reflected_heat_probability(),
        ];
        for (actual, expected) in actual.into_iter().zip(reference) {
            if expected < 1e-12 {
                assert!(
                    (actual / expected - 1.0).abs() < 2e-7,
                    "alpha={alpha}, actual={actual:e}, expected={expected:e}"
                );
            } else {
                near(
                    actual,
                    expected,
                    2.0 * g.resolution().unwrap().linear_solve_estimate + 2e-10,
                );
            }
        }
    }
    // The new angular reset fraction must not be inferred from its energy
    // fraction: wall collision count and final exit direction are correlated.
    let row = &rows[1];
    let angular_reset_fraction =
        row.reset_transmitted_axial_cosine() / total.wall_transmitted_axial_cosine();
    let energy_reset_fraction =
        row.accommodation().transmitted_heat_probability() / g.wall_transmitted_probability();
    assert!((angular_reset_fraction - energy_reset_fraction).abs() > 1e-4);
}

#[test]
fn actual_three_dimensional_rays_check_joint_angle_and_collision_history() {
    let coefficients = [0.2, 0.7];
    let (_, _, rows) = geometry(2.0)
        .resolve_axial_accommodation_batch(&coefficients, controls(), || false)
        .unwrap();
    let histories = 160_000;
    let mut reference = [[0.0; 6]; 2];
    trace_histories(2.0, histories, |exit, cosine, collisions| {
        let side = usize::from(exit == 2);
        for (index, alpha) in coefficients.into_iter().enumerate() {
            let log_retained = collisions as f64 * (-alpha).ln_1p();
            let retained = log_retained.exp();
            let reset = -log_retained.exp_m1();
            reference[index][2 * side] += retained * cosine;
            reference[index][2 * side + 1] += reset * cosine;
            reference[index][4 + side] += reset;
        }
    });
    for (row, reference) in rows.iter().zip(reference) {
        let actual = [
            split(row)[0],
            split(row)[1],
            split(row)[2],
            split(row)[3],
            row.accommodation().transmitted_heat_probability(),
            row.accommodation().reflected_heat_probability(),
        ];
        for (actual, ray) in actual.into_iter().zip(reference) {
            near(actual, ray / histories as f64, 0.0035);
        }
    }
}

#[test]
fn coefficient_order_and_duplicates_share_geometric_and_discounted_work() {
    let (g, total, rows) = geometry(1.0)
        .resolve_axial_accommodation_batch(&[0.2, 0.8], controls(), || false)
        .unwrap();
    let (other, other_total, repeated) = geometry(1.0)
        .resolve_axial_accommodation_batch(&[0.8, 0.2, 0.2, 1.0, 0.0], controls(), || false)
        .unwrap();
    assert_eq!(g, other);
    assert_eq!(total, other_total);
    assert_eq!(repeated[1], repeated[2]);
    for (original, repeated) in rows.iter().zip([&repeated[1], &repeated[0]]) {
        assert_eq!(split(original), split(repeated));
        assert_eq!(
            original.accommodation().transmitted_heat_probability(),
            repeated.accommodation().transmitted_heat_probability()
        );
        assert_eq!(
            original.accommodation().reflected_heat_probability(),
            repeated.accommodation().reflected_heat_probability()
        );
        assert_eq!(
            original.resolution().unwrap().iterations,
            repeated.resolution().unwrap().iterations
        );
        assert_eq!(
            original.resolution().unwrap().angular_escape_iterations,
            repeated.resolution().unwrap().angular_escape_iterations
        );
    }
    let singleton = geometry(1.0)
        .resolve_axial_accommodation_batch(&[0.2], controls(), || false)
        .unwrap();
    let maximum = geometry(1.0)
        .resolve_axial_accommodation_batch(&[0.2; 32], controls(), || false)
        .unwrap();
    assert_eq!(singleton.0, maximum.0);
    assert_eq!(singleton.1, maximum.1);
    assert_eq!(
        singleton.2[0].resolution().unwrap().iterations,
        maximum.2[31].resolution().unwrap().iterations
    );
}

#[test]
fn tiny_reset_probabilities_and_thin_bores_remain_positive_and_scale_invariant() {
    let coefficients = [1e-18, 0.4, 1.0 - f64::EPSILON];
    let (_, _, thin) = geometry(1e-12)
        .resolve_axial_accommodation_batch(&coefficients, controls(), || false)
        .unwrap();
    for (&alpha, row) in coefficients.iter().zip(&thin) {
        let expected = alpha * 1e-12 * 4.0 / (3.0 * core::f64::consts::PI);
        assert!((row.reset_transmitted_axial_cosine() / expected - 1.0).abs() < 2e-10);
        assert!((row.reset_reflected_axial_cosine() / expected - 1.0).abs() < 2e-10);
        assert!(row.retained_reflected_axial_cosine() > 0.0);
    }
    let relative = ClausingControls {
        absolute_tolerance: 0.0,
        ..controls()
    };
    let (_, _, small) = geometry(1.0)
        .resolve_axial_accommodation_batch(&[1e-18], relative, || false)
        .unwrap();
    for (value, error) in [
        small[0].reset_transmitted_axial_cosine(),
        small[0].reset_reflected_axial_cosine(),
    ]
    .into_iter()
    .zip(small[0].resolution().unwrap().reset_error_estimate)
    {
        assert!(value > 0.0 && error <= relative.relative_tolerance * value);
    }
    let reference = geometry(1.0)
        .resolve_axial_accommodation_batch(&coefficients, controls(), || false)
        .unwrap();
    for scale in [1e-90, 1e90] {
        let result =
            CylindricalMolecularGeometry::try_new(scale, scale, 100.0 * scale, "scaled joint law")
                .unwrap()
                .resolve_axial_accommodation_batch(&coefficients, controls(), || false)
                .unwrap();
        assert_eq!(result.1, reference.1);
        assert_eq!(result.2, reference.2);
    }
}

#[test]
fn invalid_batches_work_exhaustion_and_late_cancellation_publish_nothing() {
    for coefficients in [vec![], vec![0.2; 33], vec![-0.1], vec![1.1], vec![f64::NAN]] {
        assert!(
            geometry(1.0)
                .resolve_axial_accommodation_batch(&coefficients, controls(), || false)
                .is_err()
        );
    }
    let exhausted = ClausingControls {
        initial_rings: 4,
        maximum_rings: 16,
        absolute_tolerance: 1e-15,
        relative_tolerance: 0.0,
        ..controls()
    };
    assert!(
        geometry(2.0)
            .resolve_axial_accommodation_batch(&[0.2], exhausted, || false)
            .is_err()
    );
    let original = geometry(1.0);
    let mut calls = 0;
    original
        .clone()
        .resolve_axial_accommodation_batch(&[0.2, 0.8], controls(), || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, 71, calls] {
        let mut observed = 0;
        let result =
            original
                .clone()
                .resolve_axial_accommodation_batch(&[0.2, 0.8], controls(), || {
                    observed += 1;
                    observed == stop
                });
        assert!(matches!(result, Err(GasEnergyError::Cancelled)));
        assert_eq!(observed, stop);
    }
    assert_eq!(original, geometry(1.0));
}
