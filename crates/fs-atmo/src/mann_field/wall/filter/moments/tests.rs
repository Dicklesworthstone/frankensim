use super::{MannWallBoxMoments, sample};
use crate::mann::MannParams;
use crate::mann_field::fixed_step::FixedStepMannField;
use crate::mann_field::wall::MannPlaneWall;
use crate::mann_field::{HomogeneousMannField, MannQuadrature, MannSpectralModel, SpectralMode};

fn model(seed: u64, modes: usize) -> MannSpectralModel {
    MannSpectralModel::build(
        seed,
        MannParams {
            alpha_eps23: 0.4,
            length_m: 1.0,
            gamma: 3.9,
        },
        MannQuadrature {
            n_modes: modes,
            min_kl: 0.03,
            max_kl: 4.0,
        },
        [0.7, -0.2, 0.1],
    )
    .unwrap()
}

fn field() -> HomogeneousMannField {
    HomogeneousMannField::at_anchor(model(927, 16), 17, vec![0.7; 16]).unwrap()
}

fn prepare(field: &HomogeneousMannField, edges: [[f64; 3]; 3]) -> MannWallBoxMoments {
    MannWallBoxMoments::prepare(
        field.model(),
        MannPlaneWall::new(0.0).unwrap(),
        edges,
        1_000_000,
        10_000,
        || false,
    )
    .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

// Independent Gauss-Legendre integration of the original point observation;
// normalized weights represent a uniform coordinate on [-1,1].
const GAUSS: [(f64, f64); 8] = [
    (-0.9602898564975363, 0.1012285362903763 / 2.0),
    (-0.7966664774136267, 0.2223810344533745 / 2.0),
    (-0.5255324099163290, 0.3137066458778873 / 2.0),
    (-0.1834346424956498, 0.3626837833783620 / 2.0),
    (0.1834346424956498, 0.3626837833783620 / 2.0),
    (0.5255324099163290, 0.3137066458778873 / 2.0),
    (0.7966664774136267, 0.2223810344533745 / 2.0),
    (0.9602898564975363, 0.1012285362903763 / 2.0),
];

fn integrate(
    field: &HomogeneousMannField,
    point: [f64; 3],
    edges: [[f64; 3]; 3],
    mean: [f64; 3],
) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let wall = MannPlaneWall::new(0.0).unwrap();
    let mut stress = [[0.0; 3]; 3];
    let mut flux = [[0.0; 3]; 3];
    for (x, wx) in GAUSS {
        for (y, wy) in GAUSS {
            for (z, wz) in GAUSS {
                let at = core::array::from_fn(|i| {
                    point[i] + x * edges[0][i] + y * edges[1][i] + z * edges[2][i]
                });
                let u = field.sample_wall_blocked(at, wall).unwrap().u;
                let weight = wx * wy * wz;
                for i in 0..3 {
                    for j in 0..3 {
                        stress[i][j] += weight * (u[i] - mean[i]) * (u[j] - mean[j]);
                        flux[i][j] += weight * u[i] * u[j];
                    }
                }
            }
        }
    }
    (stress, flux)
}

#[test]
fn wall_moments_match_point_quadrature_and_are_positive_semidefinite() {
    let mut field = field();
    field.advance_to(-3811).unwrap();
    let cases = [
        ([[0.18, 0.04, 0.0], [-0.02, 0.11, 0.0], [0.0; 3]], 0.0),
        ([[0.17, -0.06, 0.04], [0.0; 3], [0.0; 3]], 0.0),
        ([[0.12, 0.02, 0.03], [-0.03, 0.1, -0.02], [0.0; 3]], 0.0),
        (
            [[0.12, 0.02, 0.03], [-0.03, 0.1, -0.02], [0.02, 0.01, 0.07]],
            0.0,
        ),
        (
            [[0.12, 0.02, 0.03], [-0.03, 0.1, -0.02], [0.02, 0.01, 0.07]],
            0.9,
        ),
    ];
    for (edges, clearance) in cases {
        let observation = prepare(&field, edges);
        let point = [
            0.37,
            -0.29,
            observation.mean_filter().vertical_radius_m() + clearance,
        ];
        let actual = field
            .sample_wall_blocked_filtered_moments(point, &observation, || false)
            .unwrap();
        assert_eq!(
            actual.mean,
            field
                .sample_wall_blocked_filtered(point, observation.mean_filter())
                .unwrap()
        );
        let (expected, flux) = integrate(&field, point, edges, actual.mean.u);
        let actual_flux = actual.momentum_flux_per_density();
        let trace = (0..3).map(|i| expected[i][i]).sum::<f64>();
        assert!(trace > 0.0);
        for i in 0..3 {
            for j in 0..3 {
                close(
                    actual.subfilter_stress_m2_s2[i][j],
                    expected[i][j],
                    2e-11 * trace,
                );
                close(
                    actual_flux[i][j],
                    flux[i][j],
                    3e-12 * (1.0 + flux[i][j].abs()),
                );
                assert_eq!(
                    actual.subfilter_stress_m2_s2[i][j],
                    actual.subfilter_stress_m2_s2[j][i]
                );
            }
        }
        for direction in [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, -2.0, 3.0],
            [3.0, 1.0, -2.0],
        ] {
            let mut quadratic = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    quadratic += direction[i] * actual.subfilter_stress_m2_s2[i][j] * direction[j];
                }
            }
            assert!(
                quadratic >= -1e-12 * trace,
                "negative energy: {quadratic:e}"
            );
        }
        if point[2] == 0.0 {
            assert_eq!(actual.mean.u[2], 0.0);
            assert_eq!(actual.subfilter_stress_m2_s2[2], [0.0; 3]);
        }
    }
}

#[test]
fn small_wall_footprints_retain_the_physical_gradient_limit() {
    let field = field();
    let wall = MannPlaneWall::new(0.0).unwrap();
    for (base, clearance) in [
        ([[0.2, 0.03, 0.0], [-0.04, 0.15, 0.0], [0.0; 3]], 0.0),
        (
            [[0.2, 0.03, 0.04], [-0.04, 0.15, -0.02], [0.01, -0.02, 0.1]],
            0.8,
        ),
        (
            [[0.2, 0.03, 0.04], [-0.04, 0.15, -0.02], [0.01, -0.02, 0.1]],
            0.0,
        ),
    ] {
        for scale in [1e-3, 1e-6, 1e-10] {
            let edges = base.map(|edge| edge.map(|x| x * scale));
            let observation = prepare(&field, edges);
            let point = [
                0.37,
                -0.29,
                observation.mean_filter().vertical_radius_m() + clearance,
            ];
            let gradient = field.sample_wall_blocked(point, wall).unwrap().grad;
            let mut expected = [[0.0; 3]; 3];
            for edge in edges {
                let increment: [f64; 3] =
                    core::array::from_fn(|i| (0..3).map(|j| gradient[i][j] * edge[j]).sum());
                for i in 0..3 {
                    for j in 0..3 {
                        expected[i][j] += increment[i] * increment[j] / 3.0;
                    }
                }
            }
            let trace = (0..3).map(|i| expected[i][i]).sum::<f64>();
            assert!(trace > 0.0);
            let actual = field
                .sample_wall_blocked_filtered_moments(point, &observation, || false)
                .unwrap();
            for i in 0..3 {
                for j in 0..3 {
                    close(
                        actual.subfilter_stress_m2_s2[i][j],
                        expected[i][j],
                        trace * (2e-11 + 10.0 * scale * scale),
                    );
                }
            }
            assert!(actual.subfilter_kinetic_energy() > 0.0);
        }
    }
}

#[test]
fn footprint_subdivision_conserves_momentum_flux() {
    let field = field();
    let edges = [[0.17, 0.03, 0.04], [-0.04, 0.13, -0.02], [0.01, 0.02, 0.11]];
    let point = [0.37, -0.29, 0.6];
    let whole = field
        .sample_wall_blocked_filtered_moments(point, &prepare(&field, edges), || false)
        .unwrap();
    for split in 0..3 {
        let mut half = edges;
        half[split] = half[split].map(|x| x * 0.5);
        let observation = prepare(&field, half);
        let children = [-1.0, 1.0].map(|sign| {
            let center = core::array::from_fn(|i| point[i] + sign * half[split][i]);
            field
                .sample_wall_blocked_filtered_moments(center, &observation, || false)
                .unwrap()
        });
        let flux = children
            .each_ref()
            .map(|child| child.momentum_flux_per_density());
        for i in 0..3 {
            close(
                whole.mean.u[i],
                0.5 * (children[0].mean.u[i] + children[1].mean.u[i]),
                2e-13,
            );
            for j in 0..3 {
                close(
                    whole.momentum_flux_per_density()[i][j],
                    0.5 * (flux[0][i][j] + flux[1][i][j]),
                    2e-12,
                );
            }
        }
    }
}

#[test]
fn wall_touching_deep_footprint_retains_the_thin_wall_layer() {
    // One solenoidal mode u=(-sin(x)e^-z,0,cos(x)(1-e^-z)). The center
    // wall amplitude underflows, but its finite vertical average does not.
    let mut model = model(927, 1);
    model.modes[0] = SpectralMode {
        k: [1.0, 0.0, 0.0],
        weight: 1.0,
        polarization: [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
    };
    let wall = MannPlaneWall::new(0.0).unwrap();
    let height = 1000.0;
    let observation = MannWallBoxMoments::prepare(
        &model,
        wall,
        [[0.0, 0.0, height], [0.0; 3], [0.0; 3]],
        1000,
        3,
        || false,
    )
    .unwrap();
    let actual = sample(
        &model,
        0.0,
        [&[1.0], &[0.0], &[0.0], &[0.0]],
        [0.0, 0.0, height],
        &observation,
        || false,
    )
    .unwrap();
    let exponential_mean = 1.0 / (2.0 * height);
    let exponential_square_mean = 1.0 / (4.0 * height);
    close(actual.mean.u[2], 1.0 - exponential_mean, 1e-15);
    close(
        actual.subfilter_stress_m2_s2[2][2],
        exponential_square_mean - exponential_mean.powi(2),
        1e-17,
    );
    assert_eq!(actual.subfilter_stress_m2_s2[0], [0.0; 3]);
    assert_eq!(actual.subfilter_stress_m2_s2[1], [0.0; 3]);
}

#[test]
fn footprint_orientation_and_edge_order_do_not_change_the_observation() {
    let field = field();
    let edges = [[0.17, 0.03, 0.04], [-0.04, 0.13, -0.02], [0.01, 0.02, 0.11]];
    let point = [0.37, -0.29, 0.6];
    let original = field
        .sample_wall_blocked_filtered_moments(point, &prepare(&field, edges), || false)
        .unwrap();
    for alternate in [
        [edges[2], edges[0], edges[1]],
        [edges[0].map(|x| -x), edges[1], edges[2]],
        edges.map(|edge| edge.map(|x| -x)),
    ] {
        let actual = field
            .sample_wall_blocked_filtered_moments(point, &prepare(&field, alternate), || false)
            .unwrap();
        for i in 0..3 {
            close(actual.mean.u[i], original.mean.u[i], 2e-13);
            for j in 0..3 {
                close(
                    actual.subfilter_stress_m2_s2[i][j],
                    original.subfilter_stress_m2_s2[i][j],
                    2e-13,
                );
            }
        }
    }
}

#[test]
fn tiny_tangential_footprints_keep_quartic_normal_stress_near_the_wall() {
    let mut model = model(927, 1);
    model.modes[0] = SpectralMode {
        k: [1.0, 0.0, 0.0],
        weight: 1.0,
        polarization: [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
    };
    let wall = MannPlaneWall::new(0.0).unwrap();
    let height: f64 = 1e-12;
    let half_width: f64 = 1e-8;
    let observation = MannWallBoxMoments::prepare(
        &model,
        wall,
        [[half_width, 0.0, 0.0], [0.0; 3], [0.0; 3]],
        1000,
        1,
        || false,
    )
    .unwrap();
    let actual = sample(
        &model,
        0.0,
        [&[1.0], &[0.0], &[0.0], &[0.0]],
        [0.0, 0.0, height],
        &observation,
        || false,
    )
    .unwrap();
    // w=(1-exp(-h))*cos(x). Its leading central cosine moment is d^4/45,
    // even though both rounded uncentered cosine moments are exactly one.
    let normal = (-height).exp_m1().powi(2) * half_width.powi(4) / 45.0;
    let tangential = (-2.0 * height).exp() * half_width.powi(2) / 3.0;
    close(actual.subfilter_stress_m2_s2[2][2], normal, 3e-14 * normal);
    close(
        actual.subfilter_stress_m2_s2[0][0],
        tangential,
        3e-14 * tangential,
    );
    assert_eq!(actual.subfilter_stress_m2_s2[0][2], 0.0);
}

#[test]
fn point_budgets_cancellation_and_source_binding_are_explicit() {
    let field = field();
    let wall = MannPlaneWall::new(0.0).unwrap();
    let point =
        MannWallBoxMoments::prepare(field.model(), wall, [[0.0; 3]; 3], 0, 0, || false).unwrap();
    let actual = field
        .sample_wall_blocked_filtered_moments([0.2, 0.3, 0.0], &point, || false)
        .unwrap();
    assert_eq!(
        actual.mean,
        field.sample_wall_blocked([0.2, 0.3, 0.0], wall).unwrap()
    );
    assert_eq!(actual.subfilter_stress_m2_s2, [[0.0; 3]; 3]);
    assert_eq!(point.workspace_bytes(), 0);
    assert_eq!(point.function_pairs(), 0);
    for edges in [
        [[0.2, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0; 3]],
        [[0.2, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.05]],
    ] {
        let observation = prepare(&field, edges);
        let bytes = observation.workspace_bytes();
        let pairs = observation.function_pairs();
        assert_eq!(pairs, if edges[2][2] == 0.0 { 136 } else { 528 });
        assert!(
            MannWallBoxMoments::prepare(field.model(), wall, edges, bytes - 1, pairs, || false)
                .is_err()
        );
        assert!(
            MannWallBoxMoments::prepare(field.model(), wall, edges, bytes, pairs - 1, || false)
                .is_err()
        );
        assert_eq!(
            MannWallBoxMoments::prepare(field.model(), wall, edges, bytes, pairs, || false)
                .unwrap(),
            observation
        );
        let mut polls = 0;
        assert!(
            MannWallBoxMoments::prepare(field.model(), wall, edges, bytes, pairs, || {
                polls += 1;
                polls == 50
            })
            .is_err()
        );
        let before = field.checkpoint();
        let coefficient_address = observation.coefficients.as_ptr();
        let mut polls = 0;
        assert!(
            field
                .sample_wall_blocked_filtered_moments([0.2, 0.3, 0.4], &observation, || {
                    polls += 1;
                    polls == 40
                })
                .is_err()
        );
        assert_eq!(field.checkpoint(), before);
        assert_eq!(observation.coefficients.as_ptr(), coefficient_address);
        let different = HomogeneousMannField::at_anchor(model(928, 16), 17, vec![0.7; 16]).unwrap();
        assert!(
            different
                .sample_wall_blocked_filtered_moments([0.2, 0.3, 0.4], &observation, || false)
                .is_err()
        );
        if edges[2][2] > 0.0 {
            assert!(
                field
                    .sample_wall_blocked_filtered_moments([0.2, 0.3, 0.049], &observation, || false)
                    .is_err()
            );
        }
    }
}

#[test]
fn both_temporal_owners_observe_the_same_snapshot_and_resume_exactly() {
    let mut source = field();
    source.advance_to(-3811).unwrap();
    let edges = [[0.17, 0.03, 0.04], [-0.04, 0.13, -0.02], [0.01, 0.02, 0.11]];
    let observation = prepare(&source, edges);
    let point = [0.37, -0.29, 0.6];
    let source_before = source.checkpoint();
    let source_sample = source
        .sample_wall_blocked_filtered_moments(point, &observation, || false)
        .unwrap();
    let mut fixed = FixedStepMannField::from_snapshot(source.clone(), 91, 0.003, 100_000).unwrap();
    let fixed_before = fixed.checkpoint();
    assert_eq!(
        fixed
            .sample_wall_blocked_filtered_moments(point, &observation, || false)
            .unwrap(),
        source_sample
    );
    assert_eq!(source.checkpoint(), source_before);
    assert_eq!(fixed.checkpoint(), fixed_before);
    assert_eq!(
        HomogeneousMannField::resume(source_before)
            .unwrap()
            .sample_wall_blocked_filtered_moments(point, &observation, || false)
            .unwrap(),
        source_sample
    );
    fixed.advance_to(7, 10_000, || false).unwrap();
    let evolved = fixed
        .sample_wall_blocked_filtered_moments(point, &observation, || false)
        .unwrap();
    assert_ne!(evolved, source_sample);
    let resumed = FixedStepMannField::resume(fixed.checkpoint(), 100_000).unwrap();
    assert_eq!(
        resumed
            .sample_wall_blocked_filtered_moments(point, &observation, || false)
            .unwrap(),
        evolved
    );
}
