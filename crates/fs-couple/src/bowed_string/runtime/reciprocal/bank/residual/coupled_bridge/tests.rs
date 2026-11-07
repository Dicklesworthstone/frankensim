use super::*;

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-12 * a.abs().max(b.abs()).max(1e-15),
        "{a:e} != {b:e}"
    );
}

#[test]
fn full_static_and_first_order_viscous_operators_match_independent_inverse() {
    for g in [
        [0.002, 0.0012, 0.001],
        [0.002, 0.002, 0.002],
        [0.002, -0.002, 0.002],
    ] {
        for (k, c) in [
            ([3000.0, 1200.0], [0.1, 0.04]),
            ([0.0, 1200.0], [0.1, 0.04]),
            ([3000.0, 0.0], [0.1, 0.04]),
            ([0.0, 0.0], [0.1, 0.04]),
            ([3000.0, 1200.0], [0.0, 0.04]),
        ] {
            let pair = Pair::new(g, k, c).unwrap();
            let determinant =
                1.0 + k[0] * g[0] + k[1] * g[2] + k[0] * k[1] * (g[0] * g[2] - g[1] * g[1]);
            let stiffness = [
                k[0] * (1.0 + k[1] * g[2]) / determinant,
                -k[0] * k[1] * g[1] / determinant,
                k[1] * (1.0 + k[0] * g[0]) / determinant,
            ];
            let t = [
                (1.0 + k[1] * g[2]) / determinant,
                -k[1] * g[1] / determinant,
                -k[0] * g[1] / determinant,
                (1.0 + k[0] * g[0]) / determinant,
            ];
            for x in [[0.3, -0.7], [0.6, 0.4], [1.0, 0.0], [0.0, 1.0]] {
                let second = x[1] - pair.mix * x[0];
                let spring = pair.first.stiffness_n_m * x[0] * x[0]
                    + pair.second.stiffness_n_m * second * second;
                let loss = pair.first.damping_n_s_m * (x[0] - pair.back * second).powi(2)
                    + pair.second.damping_n_s_m * second * second;
                let expected = stiffness[0] * x[0] * x[0]
                    + 2.0 * stiffness[1] * x[0] * x[1]
                    + stiffness[2] * x[1] * x[1];
                close(spring, expected);
                close(
                    loss,
                    c[0] * (t[0] * x[0] + t[1] * x[1]).powi(2)
                        + c[1] * (t[2] * x[0] + t[3] * x[1]).powi(2),
                );
                assert!(spring >= 0.0 && loss >= 0.0);
                // Arbitrary independent virtual displacement: the force uses
                // the transpose of the SAME factored displacement map.
                let y = [-0.2, 0.8];
                let virtual_work = pair.first.stiffness_n_m * x[0] * y[0]
                    + pair.second.stiffness_n_m * second * (y[1] - pair.mix * y[0]);
                close(
                    virtual_work,
                    y[0] * (stiffness[0] * x[0] + stiffness[1] * x[1])
                        + y[1] * (stiffness[1] * x[0] + stiffness[2] * x[1]),
                );
            }
        }
    }
}

#[test]
fn rank_one_compliance_keeps_the_weak_conditional_direction_without_a_floor() {
    let g = 0.002;
    let k = 1e30;
    let pair = Pair::new([g, g, g], [k, k], [0.1, 0.04]).unwrap();
    close(pair.second.residual_m_n, g / (1.0 + k * g));
    close(pair.second.stiffness_n_m, k / (1.0 + k * g / (1.0 + k * g)));
    assert!(pair.second.residual_m_n > 0.0);
    let cross = 1.0 - 4.0 * f64::EPSILON;
    let weak = Pair::new([1.0, cross, 1.0], [1e20, 1e20], [0.1, 0.04]).unwrap();
    close(
        weak.second.residual_m_n,
        (1.0 - cross) * (1.0 + cross) + cross * cross / (1.0 + 1e20),
    );
    assert!(weak.second.residual_m_n > 1e-16);
    let independent = Pair::new([g, 0.0, 2.0 * g], [3000.0, 1200.0], [0.1, 0.04]).unwrap();
    assert_eq!(independent.first, condense(g, 3000.0, 0.1).unwrap());
    assert_eq!(independent.second, condense(2.0 * g, 1200.0, 0.04).unwrap());
    assert!(!independent.coupled);
}

#[test]
fn inadmissible_compliance_and_support_laws_refuse() {
    for g in [
        [-1.0, 0.0, 1.0],
        [0.0, 0.1, 1.0],
        [1.0, 1.01, 1.0],
        [1.0, f64::NAN, 1.0],
        [f64::INFINITY, 0.0, 1.0],
    ] {
        assert!(Pair::new(g, [3.0, 2.0], [0.1, 0.0]).is_err());
    }
    for k in [
        [-1.0, 1.0],
        [f64::INFINITY, 1.0],
        [f64::NAN, 1.0],
        [0.0, 1.0],
    ] {
        assert!(Pair::new([0.1, 0.04, 0.2], k, [0.0, 0.1]).is_err());
    }
}
