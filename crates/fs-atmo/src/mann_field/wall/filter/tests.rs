use super::MannWallBoxFilter;
use crate::FlowSample;
use crate::mann::MannParams;
use crate::mann_field::filter::MannBoxFilter;
use crate::mann_field::fixed_step::FixedStepMannField;
use crate::mann_field::wall::MannPlaneWall;
use crate::mann_field::{HomogeneousMannField, MannQuadrature, MannSpectralModel};
use crate::ou::STATIONARY_ANCHOR_TICK;

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
    HomogeneousMannField::at_anchor(model(927, 32), 17, vec![0.7; 32]).unwrap()
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-12 * (1.0 + b.abs()), "{a} != {b}");
}

fn close_sample(a: FlowSample, b: FlowSample) {
    for (a, b) in a.u.into_iter().zip(b.u) {
        close(a, b);
    }
    for (a, b) in a
        .grad
        .into_iter()
        .flatten()
        .zip(b.grad.into_iter().flatten())
    {
        close(a, b);
    }
}

fn same_bits(a: FlowSample, b: FlowSample) {
    assert_eq!(a.u.map(f64::to_bits), b.u.map(f64::to_bits));
    assert_eq!(
        a.grad.map(|row| row.map(f64::to_bits)),
        b.grad.map(|row| row.map(f64::to_bits))
    );
}

// Independent Gauss-Legendre integration of the existing point observation.
// Weights are normalized for a uniform average on [-1,1].
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

fn quadrature(
    point: [f64; 3],
    edges: [[f64; 3]; 3],
    sample: impl Fn([f64; 3]) -> FlowSample,
) -> FlowSample {
    let mut average = FlowSample {
        u: [0.0; 3],
        grad: [[0.0; 3]; 3],
    };
    for (x, wx) in GAUSS {
        for (y, wy) in GAUSS {
            for (z, wz) in GAUSS {
                let at = core::array::from_fn(|i| {
                    point[i] + x * edges[0][i] + y * edges[1][i] + z * edges[2][i]
                });
                let value = sample(at);
                let weight = wx * wy * wz;
                for i in 0..3 {
                    average.u[i] += weight * value.u[i];
                    for j in 0..3 {
                        average.grad[i][j] += weight * value.grad[i][j];
                    }
                }
            }
        }
    }
    average
}

#[test]
fn g1_wall_face_tilted_face_and_volume_match_independent_point_quadrature() {
    let field = field();
    let wall = MannPlaneWall::new(-0.4).unwrap();
    for (point, edges) in [
        (
            [0.13, -0.3, -0.4],
            [[0.12, 0.09, 0.0], [-0.09, 0.12, 0.0], [0.0; 3]],
        ),
        (
            [0.13, -0.3, -0.18],
            [[0.12, 0.09, 0.04], [-0.09, 0.12, 0.07], [0.0; 3]],
        ),
        (
            [0.13, -0.3, -0.18],
            [[0.12, 0.09, 0.0], [0.04, 0.12, 0.03], [0.0, 0.0, 0.16]],
        ),
    ] {
        let filter =
            MannWallBoxFilter::prepare(field.model(), wall, edges, 32 * 32, || false).unwrap();
        let actual = field.sample_wall_blocked_filtered(point, &filter).unwrap();
        let expected = quadrature(point, edges, |p| {
            field.sample_wall_blocked(p, wall).unwrap()
        });
        close_sample(actual, expected);
        let reversed = MannWallBoxFilter::prepare(
            field.model(),
            wall,
            edges.map(|edge| edge.map(|v| -v)),
            32 * 32,
            || false,
        )
        .unwrap();
        close_sample(
            actual,
            field
                .sample_wall_blocked_filtered(point, &reversed)
                .unwrap(),
        );
    }
}

#[test]
fn g1_averaged_gradients_are_exact_solenoidal_and_preserve_wall_parity() {
    let field = field();
    let wall = MannPlaneWall::new(0.0).unwrap();
    let edges = [[0.2, 0.07, 0.04], [-0.05, 0.1, 0.02], [0.03, 0.0, 0.08]];
    let filter = MannWallBoxFilter::prepare(field.model(), wall, edges, 32 * 32, || false).unwrap();
    let point = [0.13, -0.3, 0.4];
    let value = field.sample_wall_blocked_filtered(point, &filter).unwrap();
    let trace: f64 = (0..3).map(|i| value.grad[i][i]).sum();
    assert!(trace.abs() < 1e-12);
    let step = 1e-5;
    for j in 0..3 {
        let mut plus = point;
        let mut minus = point;
        plus[j] += step;
        minus[j] -= step;
        let plus = field.sample_wall_blocked_filtered(plus, &filter).unwrap();
        let minus = field.sample_wall_blocked_filtered(minus, &filter).unwrap();
        for i in 0..3 {
            let derivative = (plus.u[i] - minus.u[i]) / (2.0 * step);
            assert!((derivative - value.grad[i][j]).abs() < 3e-8);
        }
    }
    let surface = MannWallBoxFilter::prepare(
        field.model(),
        wall,
        [[0.2, 0.07, 0.0], [-0.05, 0.1, 0.0], [0.0; 3]],
        32 * 32,
        || false,
    )
    .unwrap();
    let value = field
        .sample_wall_blocked_filtered([point[0], point[1], 0.0], &surface)
        .unwrap();
    assert_eq!(value.u[2].to_bits(), 0.0f64.to_bits());
    assert_eq!(value.grad[2][0].to_bits(), 0.0f64.to_bits());
    assert_eq!(value.grad[2][1].to_bits(), 0.0f64.to_bits());
    assert!((value.grad[0][0] + value.grad[1][1] + value.grad[2][2]).abs() < 1e-12);
}

#[test]
fn g1_complete_covariance_matches_two_footprint_integrals_at_distinct_heights() {
    let field = field();
    let wall = MannPlaneWall::new(-0.4).unwrap();
    let edge = [0.34, -0.23, 0.18];
    let filter = MannWallBoxFilter::prepare(
        field.model(),
        wall,
        [edge, [0.0; 3], [0.0; 3]],
        32 * 32,
        || false,
    )
    .unwrap();
    let first = [0.13, -0.3, -0.1];
    let second = [0.42, 0.1, 0.8];
    for lag in [0.0, 0.31, -0.21] {
        let actual = field
            .covariance_wall_blocked_filtered(first, second, lag, &filter)
            .unwrap();
        let mut expected = [[0.0; 3]; 3];
        for (x, wx) in GAUSS {
            for (y, wy) in GAUSS {
                let a = core::array::from_fn(|i| first[i] + x * edge[i]);
                let b = core::array::from_fn(|i| second[i] + y * edge[i]);
                let value = field.covariance_wall_blocked(a, b, lag, wall).unwrap();
                for i in 0..3 {
                    for j in 0..3 {
                        expected[i][j] += wx * wy * value[i][j];
                    }
                }
            }
        }
        let reverse = field
            .covariance_wall_blocked_filtered(second, first, -lag, &filter)
            .unwrap();
        for i in 0..3 {
            for j in 0..3 {
                close(actual[i][j], expected[i][j]);
                close(actual[i][j], reverse[j][i]);
            }
        }
        if lag == 0.0 {
            assert_eq!(
                actual,
                field
                    .model()
                    .covariance_wall_blocked_filtered(first, second, &filter)
                    .unwrap()
            );
        }
    }
}

#[test]
fn g0_g5_zero_footprints_and_both_clock_owners_keep_original_bits_and_state() {
    let mut field = field();
    field.advance_to(STATIONARY_ANCHOR_TICK + 11).unwrap();
    let wall = MannPlaneWall::new(0.0).unwrap();
    let point_filter =
        MannWallBoxFilter::prepare(field.model(), wall, [[0.0; 3]; 3], 0, || false).unwrap();
    assert_eq!(point_filter.workspace_bytes(), 0);
    let edges = [[0.12, 0.09, 0.0], [0.04, 0.12, 0.03], [0.0, 0.0, 0.16]];
    let filter = MannWallBoxFilter::prepare(field.model(), wall, edges, 32 * 32, || false).unwrap();
    let mut fixed = FixedStepMannField::from_snapshot(field.clone(), 123, 0.0013, 88 * 32).unwrap();
    let point = [0.1, -0.7, 0.9];
    same_bits(
        fixed.sample_wall_blocked_filtered(point, &filter).unwrap(),
        field.sample_wall_blocked_filtered(point, &filter).unwrap(),
    );
    let before = field.checkpoint();
    fixed.advance_to(3, 4 * 32 * 3, || false).unwrap();
    let fixed_before = fixed.checkpoint();
    for point in [[0.0; 3], [0.1, -0.7, 0.9], [1.2, 2.3, 0.4]] {
        same_bits(
            field
                .sample_wall_blocked_filtered(point, &point_filter)
                .unwrap(),
            field.sample_wall_blocked(point, wall).unwrap(),
        );
        same_bits(
            fixed
                .sample_wall_blocked_filtered(point, &point_filter)
                .unwrap(),
            fixed.sample_wall_blocked(point, wall).unwrap(),
        );
        let covariance = field
            .covariance_wall_blocked_filtered(point, [0.2, 0.1, 1.0], 0.23, &point_filter)
            .unwrap();
        assert_eq!(
            covariance,
            field
                .covariance_wall_blocked(point, [0.2, 0.1, 1.0], 0.23, wall)
                .unwrap()
        );
        assert_eq!(
            covariance,
            fixed
                .covariance_wall_blocked_filtered(point, [0.2, 0.1, 1.0], 0.23, &point_filter)
                .unwrap()
        );
    }
    let expected = quadrature(point, edges, |p| {
        fixed.sample_wall_blocked(p, wall).unwrap()
    });
    close_sample(
        fixed.sample_wall_blocked_filtered(point, &filter).unwrap(),
        expected,
    );
    assert_eq!(before, field.checkpoint());
    assert_eq!(fixed_before, fixed.checkpoint());
    let resumed = FixedStepMannField::resume(fixed_before, 88 * 32).unwrap();
    same_bits(
        fixed.sample_wall_blocked_filtered(point, &filter).unwrap(),
        resumed
            .sample_wall_blocked_filtered(point, &filter)
            .unwrap(),
    );
}

#[test]
fn g1_far_field_recovers_homogeneous_filter_and_deep_boxes_retain_wall_layer_average() {
    let field = field();
    let wall = MannPlaneWall::new(0.0).unwrap();
    let edges = [[0.2, 0.07, 0.04], [-0.05, 0.1, 0.02], [0.03, 0.0, 0.08]];
    let filter = MannWallBoxFilter::prepare(field.model(), wall, edges, 32 * 32, || false).unwrap();
    let homogeneous = MannBoxFilter::prepare(field.model(), edges, 8 * 32, || false).unwrap();
    let point = [0.13, -0.3, 1e8];
    close_sample(
        field.sample_wall_blocked_filtered(point, &filter).unwrap(),
        field.sample_filtered(point, &homogeneous).unwrap(),
    );
    // Independent one-mode integral: w = A cos(theta) + B sin(theta),
    // blocked by (1-exp(-z)); u_x = (-A sin(theta)+B cos(theta))*exp(-z).
    // Across z in [0,2H], the wall gain is (1-exp(-2H))/(2H), even when
    // the center's exp(-H) has already underflowed to exact zero.
    let mut model = model(19, 1);
    model.advection_mps = [0.0; 3];
    model.modes[0].k = [1.0, 0.0, 0.0];
    model.modes[0].polarization = [[0.0, 0.0, 1.0], [0.0; 3]];
    let field = HomogeneousMannField::at_anchor(model, 41, vec![1.0]).unwrap();
    let phase: f64 = 0.3;
    let a = field.paths[0].amplitudes()[0];
    let b = field.paths[1].amplitudes()[0];
    let normal = a * phase.cos() + b * phase.sin();
    let tangent = -a * phase.sin() + b * phase.cos();
    for height in [1e-7, 1e-20, 1000.0] {
        let filter = MannWallBoxFilter::prepare(
            field.model(),
            wall,
            [[0.0, 0.0, height], [0.0; 3], [0.0; 3]],
            32,
            || false,
        )
        .unwrap();
        let gain = -(-2.0f64 * height).exp_m1() / (2.0 * height);
        // This independent integral's series avoids subtracting 1-gain when
        // the wall distance is far below floating-point relative precision.
        let normal_gain = if height < 1e-4 {
            height * (1.0 - height * (2.0 / 3.0 - height / 3.0))
        } else {
            1.0 - gain
        };
        let actual = field
            .sample_wall_blocked_filtered([phase, 0.0, height], &filter)
            .unwrap();
        close(actual.u[0], tangent * gain);
        let expected_normal = normal * normal_gain;
        assert!(
            ((actual.u[2] - expected_normal) / expected_normal).abs() < 3e-12,
            "height {height}: {} != {expected_normal}",
            actual.u[2]
        );
        close(actual.grad[0][0], -normal * gain);
        close(actual.grad[2][2], normal * gain);
        if height > 1.0 {
            assert!(actual.u[0].abs() > 1e-6, "wall-layer average was erased");
        }
    }
}

#[test]
fn g0_source_budget_cancellation_and_complete_footprint_domain_are_enforced() {
    let field = field();
    let wall = MannPlaneWall::new(0.0).unwrap();
    let edges = [[0.2, 0.07, 0.04], [-0.05, 0.1, 0.02], [0.03, 0.0, 0.08]];
    let filter = MannWallBoxFilter::prepare(field.model(), wall, edges, 32 * 32, || false).unwrap();
    assert_eq!(filter.half_edges_m(), edges);
    assert_eq!(filter.wall(), wall);
    close(filter.vertical_radius_m(), 0.14);
    let bytes = filter.workspace_bytes();
    assert!(bytes > 0 && bytes <= 32 * 32);
    assert!(MannWallBoxFilter::prepare(field.model(), wall, edges, bytes - 1, || false).is_err());
    assert!(MannWallBoxFilter::prepare(field.model(), wall, edges, bytes, || true).is_err());
    let mut polls = 0;
    assert!(
        MannWallBoxFilter::prepare(field.model(), wall, edges, bytes, || {
            polls += 1;
            polls == 3
        })
        .is_err()
    );
    assert_eq!(polls, 3);
    let mut invalid_edges = edges;
    invalid_edges[1][0] = f64::NAN;
    assert!(
        MannWallBoxFilter::prepare(field.model(), wall, invalid_edges, bytes, || false).is_err()
    );
    let before = field.checkpoint();
    let radius = filter.vertical_radius_m();
    assert!(
        field
            .sample_wall_blocked_filtered([0.0, 0.0, radius], &filter)
            .is_ok()
    );
    for height in [0.5 * radius, -0.01, f64::NAN] {
        assert!(
            field
                .sample_wall_blocked_filtered([0.0, 0.0, height], &filter)
                .is_err()
        );
        assert!(
            field
                .covariance_wall_blocked_filtered([0.0, 0.0, 1.0], [0.0, 0.0, height], 0.0, &filter)
                .is_err()
        );
    }
    let other = HomogeneousMannField::at_anchor(model(928, 32), 17, vec![0.7; 32]).unwrap();
    assert!(
        other
            .sample_wall_blocked_filtered([0.0, 0.0, 1.0], &filter)
            .is_err()
    );
    assert!(
        other
            .covariance_wall_blocked_filtered([0.0, 0.0, 1.0], [0.1, 0.2, 0.8], 0.0, &filter)
            .is_err()
    );
    assert_eq!(field.checkpoint(), before);
}
