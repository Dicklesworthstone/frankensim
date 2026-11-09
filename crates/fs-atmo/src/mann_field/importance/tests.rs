use super::*;
use crate::mann::{MANN_TARGET_V1, stress_integrals};
use crate::mann_field::fixed_step::{FixedStepMannField, MannOuBridge, MannOuBridgeFrame};

fn quadrature(n_modes: usize) -> MannQuadrature {
    MannQuadrature {
        n_modes,
        min_kl: 0.01,
        max_kl: 1000.0,
    }
}

fn target(gamma: f64) -> MannParams {
    MannParams {
        gamma,
        ..MANN_TARGET_V1
    }
}

fn matrix_error(got: [[f64; 3]; 3], expected: [[f64; 3]; 3], scale: [[f64; 3]; 3]) -> f64 {
    let mut sum = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            let e = (got[i][j] - expected[i][j]) / det::sqrt(scale[i][i] * scale[j][j]);
            sum += e * e / 9.0;
        }
    }
    sum
}

#[test]
fn g0_positive_volume_weights_match_an_independent_coordinate_jacobian() {
    for gamma in [0.0, 3.9, 10.0] {
        let prepared = Prepared::new(target(gamma), quadrature(64)).unwrap();
        for unit in [[0.23, 0.31, 0.43], [0.61, 0.73, 0.27], [0.91, 0.47, 0.82]] {
            let (_, weight) = prepared.draw(unit).unwrap();
            let mut derivative = [[0.0; 3]; 3];
            for j in 0..3 {
                let mut plus = unit;
                let mut minus = unit;
                plus[j] += 1e-7;
                minus[j] -= 1e-7;
                let a = prepared.draw(plus).unwrap().0;
                let b = prepared.draw(minus).unwrap().0;
                for i in 0..3 {
                    derivative[i][j] = (a[i] - b[i]) / 2e-7;
                }
            }
            let [a, b, c] = derivative;
            let determinant = a[0] * (b[1] * c[2] - b[2] * c[1])
                - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]);
            // One partner represents its conjugate too; the reflected partner
            // is a separate mode with a separate weight and stochastic state.
            let expected = 2.0 * determinant.abs() / 32.0;
            assert!(weight > 0.0 && (weight / expected - 1.0).abs() < 2e-6);
        }
    }
}

#[test]
fn g0_every_mode_retains_the_full_authoritative_tensor_and_independent_geometry_pairs() {
    for gamma in [0.0, 3.9, 10.0] {
        let p = target(gamma);
        let model = MannSpectralModel::build_with_rule(
            43,
            p,
            quadrature(64),
            [3.0, 0.2, 0.0],
            MannQuadratureRule::BalancedEnergyV1,
        )
        .unwrap();
        assert_eq!(model.mode_count(), 64);
        assert_eq!(
            model.quadrature_rule(),
            MannQuadratureRule::BalancedEnergyV1
        );
        for pair in model.modes.chunks_exact(2) {
            assert_eq!(pair[0].weight, pair[1].weight);
            assert_eq!(pair[1].k, [pair[0].k[0], -pair[0].k[1], pair[0].k[2]]);
            for mode in pair {
                let phi = mann_tensor(&p, mode.k).unwrap();
                let scale = mode.weight * (phi[0][0] + phi[1][1] + phi[2][2]);
                assert!(scale.is_finite() && scale > 0.0);
                for i in 0..3 {
                    for j in 0..3 {
                        let factored = mode.polarization[0][i] * mode.polarization[0][j]
                            + mode.polarization[1][i] * mode.polarization[1][j];
                        assert!((factored - mode.weight * phi[i][j]).abs() < 3e-12 * scale);
                    }
                }
                for vector in mode.polarization {
                    let norm = det::sqrt(dot(mode.k, mode.k) * dot(vector, vector));
                    assert!(dot(mode.k, vector).abs() < 3e-14 * norm);
                }
            }
        }
        let zero = model.covariance([0.0; 3]).unwrap();
        assert_eq!(zero[0][1], 0.0);
        assert_eq!(zero[1][2], 0.0);
        // Reflection forbids these ONLY at zero separation/equal time.
        let oblique = model.covariance([0.4, 0.7, -0.2]).unwrap();
        assert!(oblique[0][1].abs() + oblique[1][2].abs() > 1e-4);
    }
}

#[test]
fn g1_held_out_geometry_seeds_reduce_full_stress_error_at_both_modal_budgets() {
    for gamma in [0.0, 1.0, 3.9, 7.0, 10.0] {
        let p = target(gamma);
        // Independent, much denser Gauss-Legendre continuum integration in
        // mann.rs; neither proposal masses nor synthesized factors enter it.
        let expected = stress_integrals(&p).unwrap();
        for count in [64, 256] {
            let q = quadrature(count);
            let prepared = Prepared::new(p, q).unwrap();
            let mut old_error = 0.0;
            let mut new_error = 0.0;
            let mut diagonal_error = [0.0; 3];
            for seed in 32..64 {
                let old = MannSpectralModel::build(seed, p, q, [0.0; 3])
                    .unwrap()
                    .covariance([0.0; 3])
                    .unwrap();
                let new = prepared
                    .build(seed, [0.0; 3])
                    .unwrap()
                    .covariance([0.0; 3])
                    .unwrap();
                old_error += matrix_error(old, expected, expected) / 32.0;
                new_error += matrix_error(new, expected, expected) / 32.0;
                for i in 0..3 {
                    diagonal_error[i] += (new[i][i] / expected[i][i] - 1.0).powi(2) / 32.0;
                }
                assert_eq!(new[0][1], 0.0);
                assert_eq!(new[1][2], 0.0);
            }
            assert!(
                new_error < 0.5 * old_error,
                "gamma={gamma}, count={count}: {new_error} vs {old_error}"
            );
            let bound = if count == 64 { 0.055 } else { 0.022 };
            assert!(new_error.sqrt() < bound);
            let diagonal_bound = if count == 64 { 0.105 } else { 0.045 };
            assert!(
                diagonal_error
                    .into_iter()
                    .all(|e| e.sqrt() < diagonal_bound)
            );
        }
    }
}

fn symmetric([xx, yy, zz, xy, xz, yz]: [f64; 6]) -> [[f64; 3]; 3] {
    [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]
}

#[test]
fn g1_full_tensor_covariance_error_improves_against_refined_continuum_integrals() {
    // Independent integration of Phi(k) cos(k.r) over 0.05<=kL<=20,
    // using Gauss-Legendre log-radius/polar angle/quarter-azimuth, with
    // explicit reflected and conjugate contributions. Refinement from
    // 64x64x32 to 96x96x48 changed any component by at most 2.43e-8.
    // These are continuum tensor integrals, not snapshots of this sampler.
    let cases = [
        (
            0.0,
            [
                [
                    0.25887063086847717,
                    0.25887063086847706,
                    0.25887063086844359,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    0.09658352049396240,
                    0.06155454020191522,
                    0.06544664912324429,
                    0.01089790497974783,
                    -0.01634685746962351,
                    -0.00467053070560524,
                ],
                [
                    0.02599785151410620,
                    0.00698216192236643,
                    0.01030990760091960,
                    -0.00998323703566281,
                    0.01331098271421741,
                    -0.00570470687752066,
                ],
            ],
        ),
        (
            3.9,
            [
                [
                    0.93461519890175937,
                    0.45297619707776243,
                    0.21795672100091520,
                    0.0,
                    -0.23960691147563956,
                    0.0,
                ],
                [
                    0.68189902989688245,
                    0.14671154289519439,
                    0.06873587958109742,
                    0.03927208114141804,
                    -0.17572616564361507,
                    -0.02023002757081430,
                ],
                [
                    0.52823532465858347,
                    0.21140680415469668,
                    0.03860077008061556,
                    -0.00138733893319954,
                    -0.10305926278320408,
                    0.01074867591806379,
                ],
            ],
        ),
        (
            10.0,
            [
                [
                    3.25034043693285524,
                    0.73521478572476806,
                    0.12459744000671548,
                    0.0,
                    -0.38641120949852359,
                    0.0,
                ],
                [
                    2.58695786898845270,
                    0.25213192055440997,
                    0.04131610122791900,
                    0.10107661857647726,
                    -0.24800769772941347,
                    -0.02005112456324518,
                ],
                [
                    2.05986274791965052,
                    0.31009830917353742,
                    0.02176984038908319,
                    0.13598575395237969,
                    -0.15621438453051473,
                    -0.00942416362969714,
                ],
            ],
        ),
    ];
    for (gamma, values) in cases {
        let expected = values.map(symmetric);
        let p = target(gamma);
        let q = MannQuadrature {
            n_modes: 256,
            min_kl: 0.05,
            max_kl: 20.0,
        };
        let prepared = Prepared::new(p, q).unwrap();
        let separations =
            [[0.0; 3], [0.7, 0.2, -0.3], [1.4, -0.6, 0.8]].map(|r| r.map(|x| x * p.length_m));
        let mut old_error = [0.0; 3];
        let mut new_error = [0.0; 3];
        for seed in 32..64 {
            let old = MannSpectralModel::build(seed, p, q, [0.0; 3]).unwrap();
            let new = prepared.build(seed, [0.0; 3]).unwrap();
            for r in 0..3 {
                old_error[r] += matrix_error(
                    old.covariance(separations[r]).unwrap(),
                    expected[r],
                    expected[0],
                ) / 32.0;
                new_error[r] += matrix_error(
                    new.covariance(separations[r]).unwrap(),
                    expected[r],
                    expected[0],
                ) / 32.0;
            }
        }
        for r in 0..3 {
            assert!(new_error[r] < old_error[r], "gamma={gamma}, separation={r}");
            assert!(new_error[r].sqrt() < 0.03);
        }
    }
}

#[test]
fn g1_reflected_modes_keep_independent_ou_lanes_and_spatially_homogeneous_statistics() {
    let model = MannSpectralModel::build_with_rule(
        7,
        target(3.9),
        quadrature(16),
        [3.0, -0.2, 0.0],
        MannQuadratureRule::BalancedEnergyV1,
    )
    .unwrap();
    let points = [[0.3, -0.6, 0.4], [0.7, 0.2, 0.1]];
    let separation = core::array::from_fn(|i| points[1][i] - points[0][i]);
    let zero = model.covariance([0.0; 3]).unwrap();
    let expected = model.covariance(separation).unwrap();
    let mut measured = [[[0.0; 3]; 3]; 3];
    let mut lane_covariance = [0.0; 4];
    for seed in 0..4096 {
        let field = HomogeneousMannField::at_anchor(model.clone(), seed, vec![0.8; 16]).unwrap();
        for (path, covariance) in lane_covariance.iter_mut().enumerate() {
            let amplitudes = field.paths[path].amplitudes();
            *covariance += amplitudes[0] * amplitudes[1] / 4096.0;
        }
        let a = field.sample(points[0]).unwrap().u;
        let b = field.sample(points[1]).unwrap().u;
        for i in 0..3 {
            for j in 0..3 {
                measured[0][i][j] += a[i] * a[j] / 4096.0;
                measured[1][i][j] += b[i] * b[j] / 4096.0;
                measured[2][i][j] += a[i] * b[j] / 4096.0;
            }
        }
    }
    assert!(lane_covariance.into_iter().all(|x| x.abs() < 0.055));
    for (got, want) in measured.into_iter().zip([zero, zero, expected]) {
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (got[i][j] - want[i][j]).abs() < 0.075 * det::sqrt(zero[i][i] * zero[j][j])
                );
            }
        }
    }
    // Unequal authored partner lifetimes need not retain reflection symmetry
    // at nonzero time lag; the original temporal model remains authoritative.
    let times = (0..16)
        .map(|i| if i % 2 == 0 { 0.2 } else { 2.0 })
        .collect();
    let field = HomogeneousMannField::at_anchor(model, 43, times).unwrap();
    let lagged = field.covariance([1.2, -0.08, 0.0], 0.4).unwrap();
    assert!(lagged[0][1].abs() + lagged[1][2].abs() > 1e-4);
}

#[test]
fn g3_both_clocks_bridges_and_filters_bind_and_replay_the_explicit_rule() {
    let p = target(3.9);
    let q = quadrature(16);
    let legacy = MannSpectralModel::build(7, p, q, [3.0, -0.2, 0.0]).unwrap();
    assert_eq!(
        legacy,
        MannSpectralModel::build_with_rule(
            7,
            p,
            q,
            [3.0, -0.2, 0.0],
            MannQuadratureRule::LegacyLogRadiusV1
        )
        .unwrap()
    );
    let model = MannSpectralModel::build_with_rule(
        7,
        p,
        q,
        [3.0, -0.2, 0.0],
        MannQuadratureRule::BalancedEnergyV1,
    )
    .unwrap();
    assert_eq!(model.modes.capacity(), legacy.modes.capacity());
    let mut field = HomogeneousMannField::at_anchor(model, 43, vec![0.8; 16]).unwrap();
    field.advance_to(field.tick() + 3).unwrap();
    let checkpoint = field.checkpoint();
    assert_eq!(
        checkpoint.quadrature_rule,
        MannQuadratureRule::BalancedEnergyV1
    );
    let mut resumed = HomogeneousMannField::resume(checkpoint.clone()).unwrap();
    assert_eq!(resumed, field);
    let point = [0.3, -0.6, 0.4];
    for edges in [[[0.0; 3]; 3], [[0.2, 0.1, 0.0], [0.0; 3], [0.0; 3]]] {
        let wrong = MannBoxFilter::prepare(&legacy, edges, 128, || false).unwrap();
        assert_eq!(
            field.sample_filtered(point, &wrong).unwrap_err().code,
            "mann-filter-source-mismatch"
        );
        let right = MannBoxFilter::prepare(field.model(), edges, 128, || false).unwrap();
        assert_eq!(
            field.sample_filtered(point, &right),
            resumed.sample_filtered(point, &right)
        );
    }
    let wrong_wall = MannWallBoxFilter::prepare(
        &legacy,
        MannPlaneWall::new(0.0).unwrap(),
        [[0.1, 0.0, 0.0], [0.0; 3], [0.0; 3]],
        1024,
        || false,
    )
    .unwrap();
    assert!(
        field
            .sample_wall_blocked_filtered(point, &wrong_wall)
            .is_err()
    );
    assert_eq!(field.checkpoint(), checkpoint);
    field.advance_to(field.tick() + 4).unwrap();
    resumed.advance_to(resumed.tick() + 4).unwrap();
    assert_eq!(field, resumed);
    let mut fixed = FixedStepMannField::from_snapshot(field.clone(), 191, 0.017, 88 * 16).unwrap();
    assert_eq!(fixed.workspace_bytes(), 88 * 16);
    assert_eq!(fixed.sample(point), field.sample(point));
    fixed.advance_to(3, 4 * 16 * 3, || false).unwrap();
    let checkpoint = fixed.checkpoint();
    let mut restored = FixedStepMannField::resume(checkpoint.clone(), 88 * 16).unwrap();
    assert_eq!(restored, fixed);
    assert_eq!(
        restored.model().quadrature_rule(),
        MannQuadratureRule::BalancedEnergyV1
    );
    let bridge = MannOuBridge::prepare(
        &fixed,
        277,
        128,
        MannOuBridge::required_workspace_bytes(16).unwrap(),
        64,
        || false,
    )
    .unwrap();
    let mut frame = bridge
        .frame(
            MannOuBridgeFrame::required_workspace_bytes(16).unwrap(),
            || false,
        )
        .unwrap();
    bridge
        .sample_frame(bridge.start_time_s(), &mut frame, 0, || false)
        .unwrap();
    assert_eq!(frame.sample(point), fixed.sample(point));
    let middle = 0.5 * (bridge.start_time_s() + bridge.end_time_s());
    let work = bridge.required_normal_draws(middle).unwrap();
    bridge
        .sample_frame(middle, &mut frame, work, || false)
        .unwrap();
    let observed = frame.sample(point).unwrap();
    assert!((observed.grad[0][0] + observed.grad[1][1] + observed.grad[2][2]).abs() < 1e-12);
    bridge
        .sample_frame(bridge.end_time_s(), &mut frame, 0, || false)
        .unwrap();
    bridge
        .sample_frame(middle, &mut frame, work, || false)
        .unwrap();
    assert_eq!(observed, frame.sample(point).unwrap());
    assert_eq!(fixed.checkpoint(), checkpoint);
    fixed.advance_to(5, 4 * 16 * 2, || false).unwrap();
    restored.advance_to(5, 4 * 16 * 2, || false).unwrap();
    assert_eq!(fixed, restored);
}

#[test]
fn g0_physical_scaling_band_admission_and_even_mode_limits_are_explicit() {
    let p = target(3.9);
    let prepared = Prepared::new(p, quadrature(64)).unwrap();
    let a = prepared.build(51, [0.0; 3]).unwrap();
    let factor: f64 = 6.0;
    let scaled = MannParams {
        alpha_eps23: p.alpha_eps23 * 9.0,
        length_m: p.length_m * factor,
        gamma: p.gamma,
    };
    let b = MannSpectralModel::build_with_rule(
        51,
        scaled,
        quadrature(64),
        [0.0; 3],
        MannQuadratureRule::BalancedEnergyV1,
    )
    .unwrap();
    for r in [[0.0; 3], [0.7, -0.4, 0.2]] {
        let left = a.covariance(r).unwrap();
        let right = b.covariance(r.map(|x| x * factor)).unwrap();
        let scale = 9.0 * factor.powf(2.0 / 3.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!((right[i][j] - scale * left[i][j]).abs() < 1e-11);
            }
        }
    }
    for n in [1, 3, 255, 257] {
        assert!(Prepared::new(p, quadrature(n)).is_err());
    }
    assert_eq!(
        Prepared::new(p, quadrature(2))
            .unwrap()
            .build(51, [0.0; 3])
            .unwrap()
            .mode_count(),
        2
    );
    assert_eq!(
        Prepared::new(p, quadrature(256))
            .unwrap()
            .build(51, [0.0; 3])
            .unwrap()
            .mode_count(),
        256
    );
    assert!(MannSpectralModel::build(51, p, quadrature(3), [0.0; 3]).is_ok());
    for (low, high) in [
        (0.0, 1.0),
        (1.0, 1.0),
        (f64::NAN, 1.0),
        (1.0, f64::INFINITY),
        (1.0, f64::MAX),
    ] {
        let q = MannQuadrature {
            n_modes: 16,
            min_kl: low,
            max_kl: high,
        };
        assert!(
            MannSpectralModel::build_with_rule(
                51,
                p,
                q,
                [0.0; 3],
                MannQuadratureRule::BalancedEnergyV1
            )
            .is_err()
        );
    }
    assert!(prepared.build(51, [f64::NAN, 0.0, 0.0]).is_err());
    for u in [[-0.1, 0.5, 0.5], [0.5, f64::NAN, 0.5], [0.5, 0.5, 1.1]] {
        assert!(prepared.draw(u).is_err());
    }
}
