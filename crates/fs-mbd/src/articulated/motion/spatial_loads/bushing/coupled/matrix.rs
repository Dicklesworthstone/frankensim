//! Fixed-size admission of explicitly supplied, work-conjugate passive forms.

/// Physical coordinate magnitudes used for numerical congruence scaling.
/// These do not change the supplied SI coefficients or the physical law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BushingCoordinateScales {
    pub translation_m: f64,
    pub rotation_rad: f64,
}

/// Invalid, active, unresolved or unrepresentable supplied matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BushingMatrixError(pub &'static str);

impl core::fmt::Display for BushingMatrixError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "passive bushing matrix: {}", self.0)
    }
}
impl std::error::Error for BushingMatrixError {}

/// Prepared symmetric positive-semidefinite form in native `[q, r]` order.
///
/// q is translation [m], r the principal rotation vector [rad]. For stiffness,
/// translation/translation entries are N/m, translation/rotation entries N/rad,
/// reciprocal entries N m/m, and rotation/rotation entries N m/rad. Damping has
/// the corresponding extra second. Radians use their SI dimensionless value;
/// generalized force is `[force N, rotation-vector force N m]`.
///
/// The supplied matrix is retained exactly. Congruence by explicit physical
/// coordinate scales precedes a diagonally equilibrated, strict semidefinite
/// Cholesky admission. Negative pivots are never clipped, asymmetry is never
/// averaged, and coupled zero diagonals refuse. Exact zero pivots require exact
/// zero residual columns. Numerically unresolved passive matrices may refuse;
/// no condition-independent rank certification is claimed.
///
/// Non-diagonal actions use the admitted Gram factors for nonnegative energy
/// and loss without cancellation near a free mode. The exposed reconstruction
/// residual bounds the factor's entrywise error in equilibrated coordinates.
/// Diagonal forms retain direct coefficient arithmetic, including exact zeros.
#[derive(Debug, Clone, PartialEq)]
pub struct PassiveBushingMatrix {
    native: [[f64; 6]; 6],
    scales: [f64; 6],
    gram: [[f64; 6]; 6],
    diagonal: bool,
    reconstruction_residual: f64,
}

impl PassiveBushingMatrix {
    /// Prepare an explicitly supplied complete matrix, including zero entries.
    /// Work is bounded by the fixed six-coordinate size; no heap allocation.
    ///
    /// # Errors
    /// Refuses invalid scales, asymmetric/nonfinite/active matrices, unresolved
    /// semidefinite pivots, or unrepresentable scaling/factor reconstruction.
    pub fn try_new(
        native: [[f64; 6]; 6],
        scales: BushingCoordinateScales,
    ) -> Result<Self, BushingMatrixError> {
        if [scales.translation_m, scales.rotation_rad]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0.0)
        {
            return Err(BushingMatrixError(
                "finite positive physical coordinate scales required",
            ));
        }
        let scales = [
            scales.translation_m,
            scales.translation_m,
            scales.translation_m,
            scales.rotation_rad,
            scales.rotation_rad,
            scales.rotation_rad,
        ];
        let mut scaled = [[0.0; 6]; 6];
        let mut diagonal = true;
        for i in 0..6 {
            for j in 0..=i {
                if !native[i][j].is_finite() || native[i][j] != native[j][i] {
                    return Err(BushingMatrixError(
                        "finite exactly symmetric native coefficients required",
                    ));
                }
                let value = native[i][j] * scales[i] * scales[j];
                if !value.is_finite() || (native[i][j] != 0.0 && value == 0.0) {
                    return Err(BushingMatrixError(
                        "physical coordinate congruence is unrepresentable",
                    ));
                }
                scaled[i][j] = value;
                scaled[j][i] = value;
                if i != j && value != 0.0 {
                    diagonal = false;
                }
            }
        }
        let mut root = [0.0; 6];
        for i in 0..6 {
            if scaled[i][i] < 0.0 || (scaled[i][i] == 0.0 && scaled[i].iter().any(|x| *x != 0.0)) {
                return Err(BushingMatrixError(
                    "negative diagonal or coupled zero diagonal",
                ));
            }
            root[i] = fs_math::det::sqrt(scaled[i][i]);
        }
        let mut equilibrated = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..=i {
                // Retain the actual scaled diagonal as well. Forcing it to one
                // can erase a representable negative pivot near a dependent mode.
                let value = if root[i] == 0.0 || root[j] == 0.0 {
                    0.0
                } else {
                    (scaled[i][j] / root[i]) / root[j]
                };
                if !value.is_finite() || (scaled[i][j] != 0.0 && value == 0.0) {
                    return Err(BushingMatrixError(
                        "matrix equilibration is unrepresentable",
                    ));
                }
                equilibrated[i][j] = value;
                equilibrated[j][i] = value;
            }
        }
        let mut factor = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..=i {
                let mut value = equilibrated[i][j];
                for k in 0..j {
                    value -= factor[i][k] * factor[j][k];
                }
                if !value.is_finite() {
                    return Err(BushingMatrixError("semidefinite factor is unrepresentable"));
                }
                if i == j {
                    if value < 0.0 {
                        return Err(BushingMatrixError(
                            "negative or unresolved semidefinite pivot",
                        ));
                    }
                    factor[i][j] = fs_math::det::sqrt(value);
                } else if factor[j][j] == 0.0 {
                    if value != 0.0 {
                        return Err(BushingMatrixError("coupled zero semidefinite pivot"));
                    }
                } else {
                    factor[i][j] = value / factor[j][j];
                    if !factor[i][j].is_finite() {
                        return Err(BushingMatrixError("semidefinite factor is unrepresentable"));
                    }
                }
            }
        }
        let mut reconstruction_residual = 0.0_f64;
        let mut gram = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                let reconstructed: f64 = (0..6).map(|k| factor[i][k] * factor[j][k]).sum();
                let residual = (reconstructed - equilibrated[i][j]).abs();
                if !residual.is_finite() || residual > 128.0 * f64::EPSILON {
                    return Err(BushingMatrixError(
                        "factor does not reconstruct the supplied form",
                    ));
                }
                reconstruction_residual = reconstruction_residual.max(residual);
                gram[i][j] = root[i] * factor[i][j];
                if !gram[i][j].is_finite()
                    || (root[i] != 0.0 && factor[i][j] != 0.0 && gram[i][j] == 0.0)
                {
                    return Err(BushingMatrixError("scaled factor is unrepresentable"));
                }
            }
        }
        Ok(Self {
            native,
            scales,
            gram,
            diagonal,
            reconstruction_residual,
        })
    }

    /// Original supplied SI coefficient values, without repair or replacement.
    pub const fn native_coefficients(&self) -> &[[f64; 6]; 6] {
        &self.native
    }

    pub const fn coordinate_scales(&self) -> BushingCoordinateScales {
        BushingCoordinateScales {
            translation_m: self.scales[0],
            rotation_rad: self.scales[3],
        }
    }

    /// Maximum absolute entry defect after diagonal equilibration, dimensionless.
    pub const fn reconstruction_residual(&self) -> f64 {
        self.reconstruction_residual
    }

    pub(super) fn rotation_active(&self) -> bool {
        (3..6).any(|i| self.native[i][i] != 0.0)
    }

    /// Original-coordinate conjugate effort and nonnegative quadratic form.
    pub(super) fn evaluate(&self, value: [f64; 6]) -> Result<([f64; 6], f64), BushingMatrixError> {
        if value.iter().any(|x| !x.is_finite()) {
            return Err(BushingMatrixError("finite native strain/rate required"));
        }
        let (effort, quadratic) = if self.diagonal {
            let force = core::array::from_fn::<_, 6, _>(|i| self.native[i][i] * value[i]);
            let quadratic = (value[0] * force[0] + value[1] * force[1] + value[2] * force[2])
                + (value[3] * force[3] + value[4] * force[4] + value[5] * force[5]);
            (force, quadratic)
        } else {
            let scaled = core::array::from_fn::<_, 6, _>(|i| value[i] / self.scales[i]);
            if scaled
                .iter()
                .zip(value)
                .any(|(x, original)| !x.is_finite() || (original != 0.0 && *x == 0.0))
            {
                return Err(BushingMatrixError(
                    "scaled native strain/rate is unrepresentable",
                ));
            }
            let modes = core::array::from_fn::<_, 6, _>(|j| {
                (0..6).map(|i| self.gram[i][j] * scaled[i]).sum::<f64>()
            });
            let force = core::array::from_fn::<_, 6, _>(|i| {
                (0..6).map(|j| self.gram[i][j] * modes[j]).sum::<f64>() / self.scales[i]
            });
            (force, modes.iter().map(|x| x * x).sum())
        };
        if effort.iter().any(|x| !x.is_finite()) || !quadratic.is_finite() {
            return Err(BushingMatrixError(
                "native effort or quadratic work is unrepresentable",
            ));
        }
        Ok((effort, quadratic))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 3e-12 * a.abs().max(b.abs()).max(1.0),
            "{a} != {b}"
        );
    }

    #[test]
    fn full_native_matrix_preserves_cross_work_under_physical_rescaling() {
        let mode = [1., -2., 0.5, 3., -1., 2.];
        let matrix = core::array::from_fn(|i| {
            core::array::from_fn(|j| mode[i] * mode[j] + if i == j { 2.0 + i as f64 } else { 0.0 })
        });
        let x = [0.02, -0.04, 0.03, 0.3, -0.2, 0.1];
        let expected: [f64; 6] =
            core::array::from_fn(|i| (0..6).map(|j| matrix[i][j] * x[j]).sum());
        let expected_work: f64 = (0..6).map(|i| x[i] * expected[i]).sum();
        for scales in [
            BushingCoordinateScales {
                translation_m: 1.,
                rotation_rad: 1.,
            },
            BushingCoordinateScales {
                translation_m: 0.002,
                rotation_rad: 0.3,
            },
        ] {
            let law = PassiveBushingMatrix::try_new(matrix, scales).unwrap();
            assert_eq!(law.native_coefficients(), &matrix);
            assert_eq!(law.coordinate_scales(), scales);
            assert!(law.reconstruction_residual() <= 128. * f64::EPSILON);
            let (force, work) = law.evaluate(x).unwrap();
            for i in 0..6 {
                close(force[i], expected[i]);
            }
            close(work, expected_work);
            close((0..6).map(|i| force[i] * x[i]).sum(), work);
            let direction = [-0.4, 0.5, 0.1, 0.2, -0.3, 0.6];
            let epsilon = 1e-5;
            let plus = law
                .evaluate(core::array::from_fn(|i| x[i] + epsilon * direction[i]))
                .unwrap()
                .1;
            let minus = law
                .evaluate(core::array::from_fn(|i| x[i] - epsilon * direction[i]))
                .unwrap()
                .1;
            let gradient = 0.5 * (plus - minus) / (2. * epsilon);
            assert!((gradient - (0..6).map(|i| force[i] * direction[i]).sum::<f64>()).abs() < 1e-9);
        }
    }

    #[test]
    fn dependent_mixed_mode_and_free_axes_have_no_artificial_loss_floor() {
        let mode = [1., 2., 0., -4., 0., 0.];
        let matrix = core::array::from_fn(|i| core::array::from_fn(|j| mode[i] * mode[j]));
        let law = PassiveBushingMatrix::try_new(
            matrix,
            BushingCoordinateScales {
                translation_m: 1.,
                rotation_rad: 1.,
            },
        )
        .unwrap();
        assert_eq!(
            law.evaluate([2., -1., 7., 0., 9., -4.]).unwrap(),
            ([0.; 6], 0.)
        );
        let (force, loss) = law.evaluate([1., 0., 0., 0., 0., 0.]).unwrap();
        assert_eq!(force, mode);
        assert_eq!(loss, 1.);
        let (_, loss) = law.evaluate([2., -1. + 1e-10, 0., 0., 0., 0.]).unwrap();
        assert!(loss > 0. && loss < 1e-18);
    }

    #[test]
    fn diagonal_action_retains_original_products_and_zero_modes() {
        let d = [3., 0., 5., 0., 7., 11.];
        let matrix =
            core::array::from_fn(|i| core::array::from_fn(|j| if i == j { d[i] } else { 0. }));
        let law = PassiveBushingMatrix::try_new(
            matrix,
            BushingCoordinateScales {
                translation_m: 0.1,
                rotation_rad: 2.,
            },
        )
        .unwrap();
        let x = [0.1, 0.2, -0.3, 0.4, -0.5, 0.6];
        let expected = core::array::from_fn(|i| d[i] * x[i]);
        let energy = (x[0] * expected[0] + x[1] * expected[1] + x[2] * expected[2])
            + (x[3] * expected[3] + x[4] * expected[4] + x[5] * expected[5]);
        assert_eq!(law.evaluate(x).unwrap(), (expected, energy));
        assert_eq!(
            law.evaluate([0., 8., 0., 9., 0., 0.]).unwrap(),
            ([0.; 6], 0.)
        );
    }

    #[test]
    fn asymmetric_active_coupled_zero_and_invalid_scales_refuse() {
        let scales = BushingCoordinateScales {
            translation_m: 1.,
            rotation_rad: 1.,
        };
        let identity =
            core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1. } else { 0. }));
        let mut asymmetric = identity;
        asymmetric[0][4] = 0.2;
        assert!(PassiveBushingMatrix::try_new(asymmetric, scales).is_err());
        let mut active = identity;
        active[0][4] = 1.0000001;
        active[4][0] = 1.0000001;
        assert!(PassiveBushingMatrix::try_new(active, scales).is_err());
        let mut near = [[0.; 6]; 6];
        near[0][0] = 2.;
        near[4][4] = 2.;
        near[0][4] = 2. + 2. * f64::EPSILON;
        near[4][0] = near[0][4];
        assert!(PassiveBushingMatrix::try_new(near, scales).is_err());
        let mut zero = [[0.; 6]; 6];
        zero[0][4] = 0.1;
        zero[4][0] = 0.1;
        assert!(PassiveBushingMatrix::try_new(zero, scales).is_err());
        let mut nonfinite = identity;
        nonfinite[3][3] = f64::NAN;
        assert!(PassiveBushingMatrix::try_new(nonfinite, scales).is_err());
        assert!(
            PassiveBushingMatrix::try_new(
                identity,
                BushingCoordinateScales {
                    translation_m: 0.,
                    rotation_rad: 1.
                }
            )
            .is_err()
        );
    }
}
