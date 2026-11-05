use super::*;
use crate::mann::MannParams;
use crate::mann_field::fixed_step::FixedStepMannField;
use crate::mann_field::{HomogeneousMannField, MannQuadrature};

fn field(seed: u64) -> HomogeneousMannField {
    let model = MannSpectralModel::build(
        927,
        MannParams {
            alpha_eps23: 0.4,
            length_m: 1.0,
            gamma: 3.9,
        },
        MannQuadrature {
            n_modes: 16,
            min_kl: 0.03,
            max_kl: 2.5,
        },
        [0.7, -0.2, 0.1],
    )
    .unwrap();
    HomogeneousMannField::at_anchor(model, seed, vec![2.0; 16]).unwrap()
}

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
fn g1_quadratic_momentum_and_stress_match_independent_spatial_integration() {
    let field = field(17);
    let point = [0.13, -0.3, 0.72];
    for edges in [
        [[0.3, 0.15, -0.1], [0.0; 3], [0.0; 3]],
        [[0.12, 0.09, 0.0], [-0.09, 0.12, 0.0], [0.0; 3]],
        [[0.12, 0.09, 0.0], [0.04, 0.12, 0.03], [0.0, 0.0, 0.16]],
    ] {
        let filter = MannBoxFilter::prepare(field.model(), edges, 128, || false).unwrap();
        let observation =
            MannBoxMoments::prepare(field.model(), edges, 2304, 136, || false).unwrap();
        let actual = field
            .sample_filtered_moments(point, &observation, || false)
            .unwrap();
        assert_eq!(actual.mean, field.sample_filtered(point, &filter).unwrap());
        let mut expected = [[0.0; 3]; 3];
        let mut stress = [[0.0; 3]; 3];
        for (x, wx) in GAUSS {
            for (y, wy) in GAUSS {
                for (z, wz) in GAUSS {
                    let at = core::array::from_fn(|i| {
                        point[i] + x * edges[0][i] + y * edges[1][i] + z * edges[2][i]
                    });
                    let sample = field.sample(at).unwrap();
                    for i in 0..3 {
                        for j in 0..3 {
                            expected[i][j] += wx * wy * wz * sample.u[i] * sample.u[j];
                            stress[i][j] += wx
                                * wy
                                * wz
                                * (sample.u[i] - actual.mean.u[i])
                                * (sample.u[j] - actual.mean.u[j]);
                        }
                    }
                }
            }
        }
        for (a, b) in actual
            .momentum_flux_per_density()
            .iter()
            .flatten()
            .zip(expected.iter().flatten())
        {
            assert!((a - b).abs() < 3e-12 * (1.0 + b.abs()), "{a} != {b}");
        }
        for (a, b) in actual
            .subfilter_stress_m2_s2
            .iter()
            .flatten()
            .zip(stress.iter().flatten())
        {
            assert!((a - b).abs() < 3e-12 * (1.0 + b.abs()), "{a} != {b}");
        }
        assert!(actual.subfilter_kinetic_energy() > 1e-7);
        for direction in [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, -2.0, 3.0],
        ] {
            let energy = (0..3)
                .map(|i| {
                    (0..3)
                        .map(|j| direction[i] * actual.subfilter_stress_m2_s2[i][j] * direction[j])
                        .sum::<f64>()
                })
                .sum::<f64>();
            assert!(energy >= -3e-13, "{energy}");
        }
    }
}

#[test]
fn g0_point_limit_budgets_and_cancellation_are_explicit() {
    let field = field(17);
    let point = [1.0, 2.0, 3.0];
    let observation =
        MannBoxMoments::prepare(field.model(), [[0.0; 3]; 3], 0, 0, || false).unwrap();
    let actual = field
        .sample_filtered_moments(point, &observation, || false)
        .unwrap();
    assert_eq!(actual.mean, field.sample(point).unwrap());
    assert_eq!(actual.subfilter_stress_m2_s2, [[0.0; 3]; 3]);
    let edges = [[0.2, 0.0, 0.0], [0.0; 3], [0.0; 3]];
    for (bytes, pairs) in [(2303, 136), (2304, 135)] {
        assert_eq!(
            MannBoxMoments::prepare(field.model(), edges, bytes, pairs, || false)
                .unwrap_err()
                .code,
            "mann-moments-budget"
        );
    }
    let observation = MannBoxMoments::prepare(field.model(), edges, 2304, 136, || false).unwrap();
    assert_eq!(observation.workspace_bytes(), 2304);
    assert_eq!(observation.mode_pairs(), 136);
    let before = field.checkpoint();
    let mut calls = 0;
    let failure = field
        .sample_filtered_moments(point, &observation, || {
            calls += 1;
            calls == 40
        })
        .unwrap_err();
    assert_eq!(failure.code, "mann-filter-cancelled");
    assert_eq!(field.checkpoint(), before);
    let mut calls = 0;
    assert_eq!(
        MannBoxMoments::prepare(field.model(), edges, 2304, 136, || {
            calls += 1;
            calls == 40
        })
        .unwrap_err()
        .code,
        "mann-filter-cancelled"
    );
}

#[test]
fn g5_both_temporal_owners_share_the_observation_and_replay() {
    let mut field = field(17);
    field.advance_to(-3811).unwrap();
    let edges = [[0.3, 0.0, 0.0], [0.1, 0.2, 0.0], [0.0, 0.0, 0.15]];
    let observation = MannBoxMoments::prepare(field.model(), edges, 2304, 136, || false).unwrap();
    let point = [0.3, -0.2, 1.0];
    let expected = field
        .sample_filtered_moments(point, &observation, || false)
        .unwrap();
    let native = FixedStepMannField::from_snapshot(field.clone(), 73, 0.01, 4096).unwrap();
    assert_eq!(
        expected,
        native
            .sample_filtered_moments(point, &observation, || false)
            .unwrap()
    );
    let restored = HomogeneousMannField::resume(field.checkpoint()).unwrap();
    assert_eq!(
        expected,
        restored
            .sample_filtered_moments(point, &observation, || false)
            .unwrap()
    );
    let mut other = field.model().clone();
    other.advection_mps[0] += 0.1;
    let other = HomogeneousMannField::at_anchor(other, 17, vec![2.0; 16]).unwrap();
    assert_eq!(
        other
            .sample_filtered_moments(point, &observation, || false)
            .unwrap_err()
            .code,
        "mann-filter-source-mismatch"
    );
}
