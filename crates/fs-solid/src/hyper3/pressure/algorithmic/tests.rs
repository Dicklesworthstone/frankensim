use super::*;

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= scale,
        "{actual:e} != {expected:e} (allowance {scale:e})"
    );
}

// Independent truncated-power divided difference of the positive affine field.
// This oracle is used only with separated nodal pressures, away from its poles.
fn oracle(p: [f64; 4]) -> f64 {
    (0..4)
        .map(|i| {
            let denominator = (0..4)
                .filter(|&j| i != j)
                .map(|j| p[i] - p[j])
                .product::<f64>();
            p[i].max(0.0).powi(4) / (4.0 * denominator)
        })
        .sum()
}

#[test]
fn coordinate_mean_closes_work_is_reversible_and_has_a_consistent_tangent() {
    for (before, after) in [
        ([0.11, -0.31, 0.83, -0.23], [0.27, 0.15, -0.42, -0.18]),
        ([-0.17, -0.32, -0.53, -0.89], [0.29, 0.47, 0.71, 1.13]),
        ([0.11, -0.31, 0.83, -0.23], [0.11, -0.31, 0.83, -0.23]),
        ([0.11, -0.31, 0.83, -0.23], [0.11, 0.17, 0.83, -0.48]),
    ] {
        let response = evaluate(before, after, true, true).unwrap();
        let reverse = evaluate(after, before, true, false).unwrap();
        let work = (0..4)
            .map(|i| response.weights[i] * (after[i] - before[i]))
            .sum::<f64>();
        near(work, oracle(after) - oracle(before), 3e-14);
        for i in 0..4 {
            assert!(response.weights[i] >= -1e-15 && response.weights[i] <= 0.25 + 1e-15);
            near(response.weights[i], reverse.weights[i], 3e-15);
        }
        for j in 0..4 {
            let h = 1e-6;
            let mut plus = after;
            let mut minus = after;
            plus[j] += h;
            minus[j] -= h;
            let fp = evaluate(before, plus, true, false).unwrap();
            let fm = evaluate(before, minus, true, false).unwrap();
            near(
                response.after_derivative[j],
                (oracle(plus) - oracle(minus)) / (2.0 * h),
                2e-9,
            );
            for i in 0..4 {
                near(
                    response.tangent_after[i][j],
                    (fp.weights[i] - fm.weights[i]) / (2.0 * h),
                    3e-10,
                );
            }
        }
        for permutation in [[1, 0, 2, 3], [3, 2, 1, 0], [2, 0, 3, 1]] {
            let permuted = evaluate(
                permutation.map(|i| before[i]),
                permutation.map(|i| after[i]),
                true,
                true,
            )
            .unwrap();
            for i in 0..4 {
                near(permuted.weights[i], response.weights[permutation[i]], 3e-15);
                for j in 0..4 {
                    near(
                        permuted.tangent_after[i][j],
                        response.tangent_after[permutation[i]][permutation[j]],
                        4e-14,
                    );
                }
            }
        }
    }
}

#[test]
fn thin_caps_and_coincident_pressures_keep_resolved_secants_without_subtraction() {
    let p = 1e-40;
    let ratio = 1.0 + 1e-8;
    let before = [p, -1.0, -2.0, -3.0];
    let after = [p * ratio, -1.0, -2.0, -3.0];
    let response = evaluate(before, after, true, true).unwrap();
    // To far below double precision, denominators are 1*2*3. This exact
    // fourth-power quotient is independent of the implementation's moment split.
    let expected = p.powi(3) * (1.0 + ratio + ratio * ratio + ratio.powi(3)) / 24.0;
    near(response.weights[0], expected, expected.abs() * 3e-15);
    let same = evaluate(before, before, true, true).unwrap();
    near(same.weights[0], p.powi(3) / 6.0, p.powi(3) * 2e-15);
    let crossing = evaluate([-p, -1.0, -2.0, -3.0], before, true, true).unwrap();
    near(crossing.weights[0], p.powi(3) / 48.0, p.powi(3) * 1e-15);
    let zero = evaluate([0.0; 4], [0.0; 4], true, true).unwrap();
    assert_eq!(zero.weights, [0.0; 4]);
    assert_eq!(zero.tangent_after, [[0.0; 4]; 4]);
    let uniform_crossing = evaluate([-p; 4], [p; 4], true, true).unwrap();
    for weight in uniform_crossing.weights {
        near(weight, 0.125, 2e-15);
    }
}

#[test]
fn a_wet_edge_landing_on_another_coplanar_vertex_has_a_finite_complete_limit() {
    for p in [1.0, 1e-40] {
        let before = [p, 0.0, -1.0, -1.0];
        let after = [0.0, 0.0, -1.0, -1.0];
        let response = evaluate(before, after, true, true).unwrap();
        let expected = p * p / (4.0 * (1.0 + p).powi(2));
        near(response.weights[0], expected, expected * 3e-15);
        near(
            response.tangent_after[0][0],
            expected / p,
            expected / p * 3e-15,
        );
        // Differentiating the fixed zero coordinate before collapse gives
        // d(m_before)/dp_1=p^2/[4(1+p)^2], whereas after collapse it is zero.
        near(
            response.tangent_after[0][1],
            0.5 * expected / p,
            expected / p * 4e-15,
        );
    }
    let before = [1e-20, 0.0, 1.0, -1.0];
    let after = [0.0, 0.0, 1.0, -1.0];
    let response = evaluate(before, after, true, true).unwrap();
    let coincident = evaluate(after, after, true, true).unwrap();
    for i in 0..4 {
        near(response.weights[i], coincident.weights[i], 3e-15);
        for j in 0..4 {
            near(
                response.tangent_after[i][j],
                coincident.tangent_after[i][j],
                3e-15,
            );
        }
    }
}
