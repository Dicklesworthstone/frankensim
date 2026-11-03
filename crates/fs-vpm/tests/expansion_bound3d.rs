//! Absolute continuous-line Taylor remainder, including heterogeneous cores.

use fs_vpm::farfield3d::{
    FarCell, FarField, WakeCoreEvolutionMode, cell_from_segs, cell_from_segs_with_cores,
};

fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
fn difference(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm(std::array::from_fn(|i| a[i] - b[i]))
}
fn mixed_cell() -> FarCell {
    cell_from_segs_with_cores(
        vec![
            ([-0.8, -0.2, 0.3], [0.6, 0.3, -0.4], 1.25),
            ([0.3, -0.7, -0.1], [-0.2, 0.5, 0.4], -0.75),
            ([-0.5, 0.4, -0.3], [0.7, -0.1, 0.2], 0.5),
        ],
        vec![0.015625, 0.25, 1.0],
    )
}
fn check_bound(cell: &FarCell, p: [f64; 3]) -> f64 {
    let bound = cell.expansion_error_bound(p);
    let error = difference(cell.eval(p), cell.exact(p));
    assert!(bound.is_finite() && bound > 0.0, "bound {bound} at {p:?}");
    assert!(error <= bound, "error {error} exceeds {bound} at {p:?}");
    bound
}

#[test]
fn heterogeneous_core_and_continuous_line_remainders_cover_all_directions_and_diffusion() {
    for common_core in [None, Some(0.0), Some(0.125), Some(16.0)] {
        let original = mixed_cell();
        let cell = common_core.map_or(original.clone(), |core| cell_from_segs(original.segs, core));
        let mut far = FarField {
            mode: WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.125 },
            shed_dt_s: 0.01,
            cells: vec![cell],
        };
        for elapsed in [0.0, 0.5, 2.0] {
            far.advance_core(elapsed).unwrap();
            let cell = &far.cells[0];
            for direction in [
                [1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.3, -0.6, 0.74],
            ] {
                for distance in [4.0, 8.0, 16.0, 32.0] {
                    let p = std::array::from_fn(|i| {
                        cell.centroid[i] + distance * direction[i] / norm(direction)
                    });
                    check_bound(cell, p);
                }
            }
        }
    }
}

#[test]
fn intrinsic_segment_quadrupole_is_required_by_an_independent_long_line_oracle() {
    let gamma = 1.25;
    let core2 = 0.125;
    let cell = cell_from_segs(vec![([-2.0, 0.0, 0.0], [2.0, 0.0, 0.0], gamma)], core2);
    // The line integral of x² from -2 to 2 is 16/3, not its midpoint value.
    assert!((cell.quad[0][0][0] - gamma * 16.0 / 3.0).abs() < 2e-15);
    let p = [0.0, 100.0, 0.0];
    let q = p[1] * p[1] + core2;
    let exact = [
        0.0,
        0.0,
        gamma * p[1] / (std::f64::consts::PI * q * (4.0 + q).sqrt()),
    ];
    let bound = cell.expansion_error_bound(p);
    assert!(difference(cell.exact(p), exact) < 2e-19);
    assert!(difference(cell.eval(p), exact) <= bound);
    let midpoint_only = [
        0.0,
        0.0,
        gamma * p[1] / (std::f64::consts::PI * q * q.sqrt()),
    ];
    assert!(
        difference(midpoint_only, exact) > bound,
        "the bound must detect missing within-segment quadrupole"
    );
}

#[test]
fn core_spatial_correlation_needs_more_than_a_pure_spatial_remainder() {
    let cell = cell_from_segs_with_cores(
        vec![
            ([-0.01, -0.02, 0.0], [0.01, -0.02, 0.0], 1.0),
            ([-0.01, 0.02, 0.0], [0.01, 0.02, 0.0], 1.0),
        ],
        vec![0.01, 1.0],
    );
    let p = [1.5, 2.5, 4.0];
    let rho2 = p.iter().map(|v| v * v).sum::<f64>() + cell.core2_m2;
    let theta = cell.radius_m / rho2.sqrt();
    let spatial_only =
        cell.w_abs_sum / (4.0 * std::f64::consts::PI * rho2) * theta.powi(3) * (4.0 - 3.0 * theta)
            / (1.0 - theta).powi(2);
    let error = difference(cell.eval(p), cell.exact(p));
    assert!(
        error > 10.0 * spatial_only,
        "mixed/core remainder must be observable"
    );
    check_bound(&cell, p);
}

#[test]
fn closed_signed_filaments_keep_a_positive_unsigned_error_budget() {
    let vertices = [
        [-1.0, -0.5, 0.0],
        [1.0, -0.5, 0.0],
        [1.0, 0.5, 0.0],
        [-1.0, 0.5, 0.0],
    ];
    let cell = cell_from_segs(
        (0..4)
            .map(|i| (vertices[i], vertices[(i + 1) % 4], 2.0))
            .collect(),
        0.0625,
    );
    assert_eq!(cell.w_sum, [0.0; 3]);
    assert!(cell.w_abs_sum > 0.0);
    let p = [3.0, 5.0, 8.0];
    assert!(difference(cell.eval(p), cell.exact(p)) > 1e-10);
    check_bound(&cell, p);
}

#[test]
fn bounds_shrink_at_cubic_spatial_order_and_preserve_length_units() {
    let cell = mixed_cell();
    let mut bounds = Vec::new();
    for distance in [32.0, 64.0, 128.0] {
        let p = [0.3 * distance, 0.6 * distance, 0.74 * distance];
        bounds.push(check_bound(&cell, p));
        assert!(cell.expansion_error_bound(p) < 0.01 * cell.contribution_bound(p));
        for scale in [0.001, 1000.0] {
            let scaled = cell_from_segs_with_cores(
                cell.segs
                    .iter()
                    .map(|(a, b, g)| (a.map(|v| v * scale), b.map(|v| v * scale), *g))
                    .collect(),
                cell.source_core2_m2
                    .iter()
                    .map(|c| c * scale * scale)
                    .collect(),
            );
            let transformed = check_bound(&scaled, p.map(|v| v * scale)) * scale;
            assert!((transformed - bounds.last().unwrap()).abs() < 5e-13 * bounds.last().unwrap());
        }
    }
    for pair in bounds.windows(2) {
        assert!((27.0..40.0).contains(&(pair[0] / pair[1])), "{bounds:?}");
    }
}

#[test]
fn near_core_spread_invalid_and_unrepresentable_domains_require_exact_sources() {
    let cell = cell_from_segs_with_cores(
        vec![
            ([0.0, -0.01, 0.0], [0.0, 0.01, 0.0], 1.0),
            ([0.02, -0.01, 0.0], [0.02, 0.01, 0.0], 1.0),
        ],
        vec![0.01, 4.0],
    );
    for p in [cell.centroid, [0.1, 0.0, 0.0]] {
        assert_eq!(cell.expansion_error_bound(p), f64::INFINITY);
        assert_eq!(cell.eval(p), cell.exact(p));
    }
    for p in [
        [f64::NAN, 0.0, 0.0],
        [f64::INFINITY, 0.0, 0.0],
        [1e100, 0.0, 0.0],
    ] {
        assert_eq!(cell.expansion_error_bound(p), f64::INFINITY);
    }
    check_bound(&cell, [32.0, 16.0, 8.0]);
    let mut invalid = cell;
    invalid.quad[0][0][0] = f64::NAN;
    assert_eq!(
        invalid.expansion_error_bound([32.0, 16.0, 8.0]),
        f64::INFINITY
    );
}
