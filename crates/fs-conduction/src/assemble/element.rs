//! Local P1 Fourier transport in the three independent temperature differences.
//!
//! The FEEC owner supplies barycentric gradients. Since their sum is zero,
//! subtracting the first nodal value before differentiation is exact algebra:
//! an absolute temperature cannot act as a heat source. This also avoids
//! cancellation of large Kelvin offsets when the physical gradient is small.

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

fn tensor_action(k: &[[f64; 3]; 3], value: [f64; 3]) -> [f64; 3] {
    k.map(|row| dot(row, value))
}

pub(super) fn gradient(shape: &[[f64; 3]; 4], nodal: [f64; 4]) -> [f64; 3] {
    let delta = [
        nodal[1] - nodal[0],
        nodal[2] - nodal[0],
        nodal[3] - nodal[0],
    ];
    std::array::from_fn(|axis| {
        delta[0].mul_add(
            shape[1][axis],
            delta[1].mul_add(shape[2][axis], delta[2] * shape[3][axis]),
        )
    })
}

/// Symmetric conductivity gives six independent entries in the edge basis.
/// Completing the first row/column enforces the same constant null mode as
/// the local action, without separately rounding sixteen tensor contractions.
pub(super) fn stiffness(volume: f64, shape: &[[f64; 3]; 4], k: &[[f64; 3]; 3]) -> [[f64; 4]; 4] {
    let mut block = [[0.0; 4]; 4];
    for a in 1..4 {
        let action = tensor_action(k, shape[a]);
        for b in a..4 {
            let entry = volume * dot(shape[b], action);
            block[a][b] = entry;
            block[b][a] = entry;
        }
    }
    for a in 1..4 {
        let entry = -(block[a][1] + block[a][2] + block[a][3]);
        block[a][0] = entry;
        block[0][a] = entry;
    }
    block[0][0] = -(block[0][1] + block[0][2] + block[0][3]);
    block
}

/// Integrated internal nodal loads; their sum vanishes by construction.
/// This applies the physical operator without assembling a local matrix.
pub(super) fn action(
    volume: f64,
    shape: &[[f64; 3]; 4],
    k: &[[f64; 3]; 3],
    nodal: [f64; 4],
) -> [f64; 4] {
    let flux = tensor_action(k, gradient(shape, nodal));
    weak_load(volume, shape, flux)
}

pub(super) fn weak_load(volume: f64, shape: &[[f64; 3]; 4], flux: [f64; 3]) -> [f64; 4] {
    let mut result = [0.0; 4];
    for a in 1..4 {
        result[a] = volume * dot(shape[a], flux);
    }
    result[0] = -(result[1] + result[2] + result[3]);
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    // Independent inverse of the edge matrix with columns
    // (3,0,0), (1,7,0), (2,4,11). Non-binary fractions expose cancellation.
    fn shape() -> [[f64; 3]; 4] {
        let g1 = [1.0 / 3.0, -1.0 / 21.0, -10.0 / 231.0];
        let g2 = [0.0, 1.0 / 7.0, -4.0 / 77.0];
        let g3 = [0.0, 0.0, 1.0 / 11.0];
        [std::array::from_fn(|i| -g1[i] - g2[i] - g3[i]), g1, g2, g3]
    }

    const K: [[f64; 3]; 3] = [[4.0, 1.0, -0.5], [1.0, 3.0, 0.25], [-0.5, 0.25, 2.0]];
    const VOLUME: f64 = 231.0 / 6.0;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() <= 8e-14 * b.abs().max(1.0), "{a:e} != {b:e}");
    }

    #[test]
    fn isothermal_state_has_no_gradient_or_internal_heat_at_any_offset() {
        for temperature in [300.1, 873.15, 1e20, 1e200] {
            assert_eq!(gradient(&shape(), [temperature; 4]), [0.0; 3]);
            assert_eq!(action(VOLUME, &shape(), &K, [temperature; 4]), [0.0; 4]);
        }
    }

    #[test]
    fn affine_gradient_and_fourier_dissipation_match_independent_formula() {
        // T=300 + 2*x - 3*y + 4*z on the declared skew tetrahedron.
        let nodal = [300.0, 306.0, 281.0, 336.0];
        let expected_gradient = [2.0, -3.0, 4.0];
        let grad = gradient(&shape(), nodal);
        for i in 0..3 {
            close(grad[i], expected_gradient[i]);
        }
        let loads = action(VOLUME, &shape(), &K, nodal);
        close(loads.iter().sum(), 0.0);
        let dissipation: f64 = loads.iter().zip(nodal).map(|(r, t)| r * (t - 300.0)).sum();
        close(dissipation, VOLUME * 49.0);
    }

    #[test]
    fn representable_temperature_translation_does_not_change_conduction() {
        let perturbation = [
            0.0,
            2.0_f64.powi(-12),
            -2.0_f64.powi(-11),
            2.0_f64.powi(-10),
        ];
        let baseline = gradient(&shape(), perturbation);
        let loads = action(VOLUME, &shape(), &K, perturbation);
        for reference in [300.0, 1024.0, 1_048_576.0] {
            let shifted = perturbation.map(|value| value + reference);
            assert_eq!(gradient(&shape(), shifted), baseline);
            assert_eq!(action(VOLUME, &shape(), &K, shifted), loads);
        }
    }

    #[test]
    fn stiffness_is_symmetric_and_matches_the_independent_weak_form() {
        let shape = shape();
        let block = stiffness(VOLUME, &shape, &K);
        for a in 0..4 {
            for b in 0..4 {
                assert_eq!(block[a][b], block[b][a]);
                let mut expected = 0.0;
                for i in 0..3 {
                    for j in 0..3 {
                        expected += VOLUME * shape[a][i] * K[i][j] * shape[b][j];
                    }
                }
                close(block[a][b], expected);
            }
            close(block[a].iter().sum(), 0.0);
        }
        let nodal = [0.0, 6.0, -19.0, 36.0];
        let loads = action(VOLUME, &shape, &K, nodal);
        for a in 0..4 {
            close(
                block[a].iter().zip(nodal).map(|(k, t)| k * t).sum(),
                loads[a],
            );
        }
    }
}
