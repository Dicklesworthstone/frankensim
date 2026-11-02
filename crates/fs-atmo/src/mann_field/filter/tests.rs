use super::*;
use crate::mann_field::HomogeneousMannField;

fn field() -> HomogeneousMannField {
    let model = MannSpectralModel::build(
        927,
        MannParams {
            alpha_eps23: 0.4,
            length_m: 1.0,
            gamma: 3.9,
        },
        MannQuadrature {
            n_modes: 32,
            min_kl: 0.03,
            max_kl: 2.5,
        },
        [0.7, -0.2, 0.1],
    )
    .unwrap();
    HomogeneousMannField::at_anchor(model, 17, vec![2.0; 32]).unwrap()
}

// Independent eight-point Gauss-Legendre nodes/weights on [-1,1]. We divide
// weights by two because this is a normalized spatial observation.
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

#[test]
fn g1_rotated_face_and_volume_averages_match_independent_spatial_quadrature() {
    let field = field();
    let point = [0.13, -0.3, 0.72];
    for edges in [
        [[0.12, 0.09, 0.0], [-0.09, 0.12, 0.0], [0.0; 3]],
        [[0.12, 0.09, 0.0], [-0.09, 0.12, 0.0], [0.0, 0.0, 0.16]],
        // The same formula also handles a nonsingular skew parallelepiped.
        [[0.12, 0.09, 0.0], [0.04, 0.12, 0.03], [0.0, 0.0, 0.16]],
    ] {
        let filter = MannBoxFilter::prepare(field.model(), edges, 32 * 8, || false).unwrap();
        let actual = field.sample_filtered(point, &filter).unwrap();
        let mut expected = crate::FlowSample {
            u: [0.0; 3],
            grad: [[0.0; 3]; 3],
        };
        for (x, wx) in GAUSS {
            for (y, wy) in GAUSS {
                for (z, wz) in GAUSS {
                    let at = core::array::from_fn(|i| {
                        point[i] + x * edges[0][i] + y * edges[1][i] + z * edges[2][i]
                    });
                    let value = field.sample(at).unwrap();
                    for i in 0..3 {
                        expected.u[i] += wx * wy * wz * value.u[i];
                        for j in 0..3 {
                            expected.grad[i][j] += wx * wy * wz * value.grad[i][j];
                        }
                    }
                }
            }
        }
        for (a, b) in actual
            .u
            .iter()
            .chain(actual.grad.iter().flatten())
            .zip(expected.u.iter().chain(expected.grad.iter().flatten()))
        {
            assert!((a - b).abs() < 2e-12 * (1.0 + b.abs()), "{a} != {b}");
        }
        let divergence: f64 = (0..3).map(|i| actual.grad[i][i]).sum();
        assert!(divergence.abs() < 1e-13);
    }
}

#[test]
fn g1_filtered_covariance_matches_two_observation_integrals_without_variance_restoration() {
    let field = field();
    let edge = [0.6, -0.3, 0.2];
    let filter =
        MannBoxFilter::prepare(field.model(), [edge, [0.0; 3], [0.0; 3]], 32 * 8, || false)
            .unwrap();
    let separation = [0.17, -0.2, 0.4];
    let lag = 0.31;
    let actual = field.covariance_filtered(separation, lag, &filter).unwrap();
    let mut expected = [[0.0; 3]; 3];
    for (x, wx) in GAUSS {
        for (y, wy) in GAUSS {
            let at = core::array::from_fn(|i| separation[i] + (x - y) * edge[i]);
            let covariance = field.covariance(at, lag).unwrap();
            for i in 0..3 {
                for j in 0..3 {
                    expected[i][j] += wx * wy * covariance[i][j];
                }
            }
        }
    }
    for (a, b) in actual.iter().flatten().zip(expected.iter().flatten()) {
        assert!((a - b).abs() < 2e-12 * (1.0 + b.abs()), "{a} != {b}");
    }
    let point = field.covariance([0.0; 3], 0.0).unwrap();
    let averaged = field.covariance_filtered([0.0; 3], 0.0, &filter).unwrap();
    for direction in [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, -2.0, 3.0],
    ] {
        let variance = |c: [[f64; 3]; 3]| -> f64 {
            (0..3)
                .map(|i| {
                    (0..3)
                        .map(|j| direction[i] * c[i][j] * direction[j])
                        .sum::<f64>()
                })
                .sum()
        };
        assert!(variance(averaged) > 0.0);
        assert!(variance(averaged) < variance(point));
    }
}

#[test]
fn g3_point_filter_is_bit_identical_and_filtered_sampling_preserves_ou_state() {
    let mut field = field();
    let point_filter = MannBoxFilter::prepare(field.model(), [[0.0; 3]; 3], 0, || false).unwrap();
    assert_eq!(point_filter.workspace_bytes(), 0);
    let box_filter = MannBoxFilter::prepare(
        field.model(),
        [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]],
        32 * 8,
        || false,
    )
    .unwrap();
    field.advance_to(field.tick() + 7).unwrap();
    let checkpoint = field.checkpoint();
    let point = [0.4, -0.3, 0.2];
    assert_eq!(
        field.sample(point).unwrap(),
        field.sample_filtered(point, &point_filter).unwrap()
    );
    assert_eq!(
        field.covariance(point, 0.1).unwrap(),
        field
            .covariance_filtered(point, 0.1, &point_filter)
            .unwrap()
    );
    let sample = field.sample_filtered(point, &box_filter).unwrap();
    assert_eq!(field.checkpoint(), checkpoint);
    let resumed = HomogeneousMannField::resume(checkpoint).unwrap();
    assert_eq!(resumed.sample_filtered(point, &box_filter).unwrap(), sample);
}

#[test]
fn g0_filter_admission_binds_every_source_input_and_enforces_budget_cancellation() {
    let field = field();
    let edges = [[0.2, 0.0, 0.0], [0.0; 3], [0.0; 3]];
    let filter = MannBoxFilter::prepare(field.model(), edges, 32 * 8, || false).unwrap();
    assert_eq!(filter.workspace_bytes(), 32 * 8);
    assert_eq!(
        MannBoxFilter::prepare(field.model(), edges, 32 * 8 - 1, || false)
            .unwrap_err()
            .code,
        "mann-filter-budget"
    );
    let mut polls = 0;
    assert_eq!(
        MannBoxFilter::prepare(field.model(), edges, 32 * 8, || {
            polls += 1;
            polls == 4
        })
        .unwrap_err()
        .code,
        "mann-filter-cancelled"
    );
    assert_eq!(polls, 4);
    assert_eq!(
        MannBoxFilter::prepare(field.model(), [[f64::NAN; 3]; 3], 32 * 8, || false)
            .unwrap_err()
            .code,
        "mann-filter-invalid"
    );
    let m = field.model();
    for change in 0..8 {
        let mut p = m.parameters;
        let mut q = m.quadrature;
        let mut seed = m.geometry_seed;
        let mut advection = m.advection_mps;
        match change {
            0 => seed += 1,
            1 => p.alpha_eps23 *= 2.0,
            2 => p.length_m *= 2.0,
            3 => p.gamma += 1.0,
            4 => q.max_kl *= 2.0,
            5 => q.min_kl *= 2.0,
            6 => q.n_modes += 1,
            _ => advection[0] += 1.0,
        }
        let model = MannSpectralModel::build(seed, p, q, advection).unwrap();
        assert_eq!(
            filter.gains(&model).unwrap_err().code,
            "mann-filter-source-mismatch"
        );
    }
}
