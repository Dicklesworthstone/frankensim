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
fn g1_refined_footprints_recover_the_full_velocity_gradient_stress() {
    let field = field(17);
    let point = [0.13, -0.3, 0.72];
    let gradient = field.sample(point).unwrap().grad;
    for shape in [
        [[1.0, 0.4, -0.2], [0.0; 3], [0.0; 3]],
        [[0.8, 0.6, 0.0], [-0.6, 0.8, 0.0], [0.0; 3]],
        [[1.0, 0.4, -0.2], [-0.3, 0.8, 0.1], [0.2, 0.1, 0.9]],
    ] {
        for scale in [1e-3, 1e-6, 1e-10] {
            let edges = shape.map(|edge| edge.map(|x| x * scale));
            // Taylor expansion of the physical velocity, with independent
            // uniform coordinates s_j: E[s_j s_k] = delta_jk/3. This oracle
            // uses the analytic velocity gradient, not modal filter gains.
            let slopes: [[f64; 3]; 3] =
                edges.map(|edge| core::array::from_fn(|i| super::super::dot(gradient[i], edge)));
            let expected: [[f64; 3]; 3] = core::array::from_fn(|i| {
                core::array::from_fn(|j| slopes.iter().map(|s| s[i] * s[j] / 3.0).sum())
            });
            let observation =
                MannBoxMoments::prepare(field.model(), edges, 2304, 136, || false).unwrap();
            let actual = field
                .sample_filtered_moments(point, &observation, || false)
                .unwrap();
            let norm = (0..3).map(|i| expected[i][i]).sum::<f64>();
            assert!(norm > 0.0);
            for i in 0..3 {
                for j in 0..3 {
                    let error = (actual.subfilter_stress_m2_s2[i][j] - expected[i][j]).abs();
                    assert!(
                        error < (2e-12 + 20.0 * scale * scale) * norm,
                        "scale={scale}, component={i},{j}: {} != {}",
                        actual.subfilter_stress_m2_s2[i][j],
                        expected[i][j],
                    );
                }
            }
            assert!(actual.subfilter_kinetic_energy() > 0.0);
            for direction in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, -2.0, 3.0]] {
                let energy = (0..3)
                    .map(|i| {
                        (0..3)
                            .map(|j| {
                                direction[i] * actual.subfilter_stress_m2_s2[i][j] * direction[j]
                            })
                            .sum::<f64>()
                    })
                    .sum::<f64>();
                assert!(energy >= -2e-12 * norm, "scale={scale}: {energy}");
            }
        }
    }
}

#[test]
fn g1_small_cosine_and_sine_moments_retain_their_distinct_physical_orders() {
    for scale in [1e-3, 1e-6, 1e-9, 1e-40] {
        let [cosine, sine, _] = edge_moments(scale, -2.0 * scale).unwrap();
        // cos(a s) = 1-a²s²/2+..., sin(a s) = a s+...;
        // Var(s²)=1/5-1/9=4/45 and E[s²]=1/3.
        let expected_cosine = scale * scale * scale * scale * (4.0 / 45.0);
        let expected_sine = scale * scale * (-2.0 / 3.0);
        let tolerance = 2e-14 + 2.0 * scale * scale;
        assert!((cosine / expected_cosine - 1.0).abs() < tolerance);
        assert!((sine / expected_sine - 1.0).abs() < tolerance);
    }
    // Both sides of the series switch agree with an independent integral.
    for (a, b) in [(0.5, -0.37), (0.500_000_000_000_000_1, -0.37)] {
        let actual = edge_moments(a, b).unwrap();
        let mean_a = GAUSS.iter().map(|&(s, w)| w * det::cos(a * s)).sum::<f64>();
        let mean_b = GAUSS.iter().map(|&(s, w)| w * det::cos(b * s)).sum::<f64>();
        let cosine = GAUSS
            .iter()
            .map(|&(s, w)| w * (det::cos(a * s) - mean_a) * (det::cos(b * s) - mean_b))
            .sum::<f64>();
        let sine = GAUSS
            .iter()
            .map(|&(s, w)| w * det::sin(a * s) * det::sin(b * s))
            .sum::<f64>();
        assert!((actual[0] - cosine).abs() < 4e-16);
        assert!((actual[1] - sine).abs() < 4e-16);
    }
}

#[test]
fn g1_edge_order_and_orientation_do_not_change_a_refined_footprint() {
    let field = field(17);
    let point = [0.13, -0.3, 0.72];
    let edges = [[1e-8, 4e-9, -2e-9], [-3e-9, 8e-9, 1e-9], [2e-9, 1e-9, 9e-9]];
    let reordered = [edges[2].map(|x| -x), edges[0], edges[1].map(|x| -x)];
    let observe = |footprint| {
        let observation =
            MannBoxMoments::prepare(field.model(), footprint, 2304, 136, || false).unwrap();
        field
            .sample_filtered_moments(point, &observation, || false)
            .unwrap()
    };
    let a = observe(edges);
    let b = observe(reordered);
    let norm = 2.0 * a.subfilter_kinetic_energy();
    assert!(norm > 0.0);
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (a.subfilter_stress_m2_s2[i][j] - b.subfilter_stress_m2_s2[i][j]).abs()
                    < 2e-13 * norm
            );
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
