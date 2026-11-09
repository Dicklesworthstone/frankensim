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
    // Explicit source factors stay in physical coordinates without a round trip.
    physical_gram: bool,
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
            physical_gram: false,
            diagonal,
            reconstruction_residual,
        })
    }

    /// Admit an original physical Gram factor without refactorization.
    ///
    /// Rows map the six work-conjugate SI coordinates to energy/loss amplitudes.
    /// All rows, cross terms and exact free directions survive. The native form
    /// is their compensated product L^T L; no pivot, eigenspace, clipping or
    /// diagonal approximation is introduced. Explicit physical scales remain
    /// admission/reporting data and never change these original factor entries.
    pub fn from_gram_rows(
        rows: [[f64; 6]; 6],
        scales: BushingCoordinateScales,
    ) -> Result<Self, BushingMatrixError> {
        use fs_math::dd::Dd;
        if [scales.translation_m, scales.rotation_rad]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0.0)
            || rows.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(BushingMatrixError(
                "finite original factors and positive physical scales required",
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
        let mut native = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..=i {
                let value = rows
                    .iter()
                    .fold(Dd::ZERO, |sum, row| {
                        sum + Dd::from_f64(row[i]) * Dd::from_f64(row[j])
                    })
                    .to_f64();
                let scaled = value * scales[i] * scales[j];
                if !value.is_finite() || !scaled.is_finite() || (value != 0.0 && scaled == 0.0) {
                    return Err(BushingMatrixError(
                        "original factor form or congruence is unrepresentable",
                    ));
                }
                native[i][j] = value;
                native[j][i] = value;
            }
            if native[i][i] == 0.0 && rows.iter().any(|row| row[i] != 0.0) {
                return Err(BushingMatrixError(
                    "nonzero original factor has unrepresentable stiffness",
                ));
            }
            for row in &rows {
                let scaled = row[i] * scales[i];
                if !scaled.is_finite() || (row[i] != 0.0 && scaled == 0.0) {
                    return Err(BushingMatrixError(
                        "original factor congruence is unrepresentable",
                    ));
                }
            }
        }
        Ok(Self {
            native,
            scales,
            gram: core::array::from_fn(|i| core::array::from_fn(|j| rows[j][i])),
            physical_gram: true,
            diagonal: false,
            reconstruction_residual: 0.0,
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

    /// Multiply this admitted form by a finite nonnegative physical coefficient.
    ///
    /// The original Gram directions are retained and scaled by the square root
    /// of the coefficient. No factorization or rank decision is repeated, so a
    /// changing viscosity cannot create a spurious negative pivot in an exactly
    /// dependent passive mode. Original coordinate scales and diagonal action
    /// are preserved. A unit multiplier returns an identical clone; an explicit
    /// zero multiplier supplies the exact zero form.
    ///
    /// Native coefficients, scaled factors and their equilibrated reconstruction
    /// are checked before returning. Nonzero coefficients or factor entries may
    /// not silently disappear through underflow. Work and storage are bounded
    /// by six coordinates; no heap allocation is required.
    pub fn scaled(&self, coefficient: f64) -> Result<Self, BushingMatrixError> {
        if !coefficient.is_finite() || coefficient < 0.0 {
            return Err(BushingMatrixError(
                "finite nonnegative matrix multiplier required",
            ));
        }
        if coefficient == 1.0 {
            return Ok(self.clone());
        }
        if coefficient == 0.0 {
            return Ok(Self {
                native: [[0.0; 6]; 6],
                scales: self.scales,
                gram: [[0.0; 6]; 6],
                physical_gram: self.physical_gram,
                diagonal: true,
                reconstruction_residual: 0.0,
            });
        }
        let product = |original: f64, multiplier: f64| {
            let result = original * multiplier;
            if !result.is_finite() || (original != 0.0 && result == 0.0) {
                Err(BushingMatrixError("scaled passive form is unrepresentable"))
            } else {
                Ok(result)
            }
        };
        let root = fs_math::det::sqrt(coefficient);
        let mut next = self.clone();
        let mut scaled = [[0.0; 6]; 6];
        let mut diagonal_roots = [0.0; 6];
        for i in 0..6 {
            for j in 0..6 {
                next.native[i][j] = product(self.native[i][j], coefficient)?;
                next.gram[i][j] = product(self.gram[i][j], root)?;
                scaled[i][j] =
                    product(product(next.native[i][j], self.scales[i])?, self.scales[j])?;
            }
            diagonal_roots[i] = fs_math::det::sqrt(scaled[i][i]);
        }
        let mut normalized = [[0.0; 6]; 6];
        for i in 0..6 {
            if diagonal_roots[i] == 0.0 {
                continue;
            }
            for j in 0..6 {
                let value = if next.physical_gram {
                    product(next.gram[i][j], self.scales[i])? / diagonal_roots[i]
                } else {
                    next.gram[i][j] / diagonal_roots[i]
                };
                if !value.is_finite() || (next.gram[i][j] != 0.0 && value == 0.0) {
                    return Err(BushingMatrixError("scaled passive factor is unresolved"));
                }
                normalized[i][j] = value;
            }
        }
        next.reconstruction_residual = 0.0;
        for i in 0..6 {
            for j in 0..=i {
                if diagonal_roots[i] == 0.0 || diagonal_roots[j] == 0.0 {
                    continue;
                }
                let expected = (scaled[i][j] / diagonal_roots[i]) / diagonal_roots[j];
                let actual = (0..6).fold(fs_math::dd::Dd::ZERO, |sum, k| {
                    sum + fs_math::dd::Dd::from_f64(normalized[i][k])
                        * fs_math::dd::Dd::from_f64(normalized[j][k])
                });
                let residual = (actual - fs_math::dd::Dd::from_f64(expected))
                    .to_f64()
                    .abs();
                if !expected.is_finite()
                    || (scaled[i][j] != 0.0 && expected == 0.0)
                    || !residual.is_finite()
                    || residual > 128.0 * f64::EPSILON
                {
                    return Err(BushingMatrixError(
                        "scaled factor does not reconstruct the supplied form",
                    ));
                }
                next.reconstruction_residual = next.reconstruction_residual.max(residual);
            }
        }
        Ok(next)
    }

    /// Maximum absolute entry defect after diagonal equilibration, dimensionless.
    pub const fn reconstruction_residual(&self) -> f64 {
        self.reconstruction_residual
    }

    pub(super) fn rotation_active(&self) -> bool {
        (3..6).any(|i| self.native[i][i] != 0.0)
    }

    // Original physical factor; exact zero columns retain undamped modes.
    pub(crate) fn gram_factor(&self) -> Result<[[f64; 6]; 6], BushingMatrixError> {
        let mut result = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                let value = if self.physical_gram {
                    self.gram[i][j]
                } else if self.diagonal {
                    if i == j {
                        fs_math::det::sqrt(self.native[i][i])
                    } else {
                        0.0
                    }
                } else {
                    self.gram[i][j] / self.scales[i]
                };
                if !value.is_finite() || (!self.diagonal && self.gram[i][j] != 0.0 && value == 0.0)
                {
                    return Err(BushingMatrixError(
                        "physical passive factor is unrepresentable",
                    ));
                }
                result[i][j] = value;
            }
        }
        Ok(result)
    }

    /// Original-coordinate conjugate effort and nonnegative quadratic form.
    pub(super) fn evaluate(&self, value: [f64; 6]) -> Result<([f64; 6], f64), BushingMatrixError> {
        if value.iter().any(|x| !x.is_finite()) {
            return Err(BushingMatrixError("finite native strain/rate required"));
        }
        let (effort, quadratic) = if self.physical_gram {
            use fs_math::dd::Dd;
            let modes = core::array::from_fn::<_, 6, _>(|j| {
                (0..6).fold(Dd::ZERO, |sum, i| {
                    sum + Dd::from_f64(self.gram[i][j]) * Dd::from_f64(value[i])
                })
            });
            let force = core::array::from_fn(|i| {
                (0..6)
                    .fold(Dd::ZERO, |sum, j| {
                        sum + Dd::from_f64(self.gram[i][j]) * modes[j]
                    })
                    .to_f64()
            });
            let quadratic = modes.iter().fold(Dd::ZERO, |sum, mode| sum + *mode * *mode);
            (force, quadratic.to_f64())
        } else if self.diagonal {
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
    fn physical_scaling_retains_dependent_directions_and_positive_loss() {
        let mode = [1.0, 2.0, 0.0, -4.0, 0.0, 0.0];
        let law = PassiveBushingMatrix::try_new(
            core::array::from_fn(|i| core::array::from_fn(|j| mode[i] * mode[j])),
            BushingCoordinateScales {
                translation_m: 1.0,
                rotation_rad: 1.0,
            },
        )
        .unwrap();
        assert_eq!(law.scaled(1.0).unwrap(), law);
        for coefficient in [1e-120, 0.1, 0.3, 7.3, 1e120] {
            let shifted = law.scaled(coefficient).unwrap();
            assert_eq!(shifted.coordinate_scales(), law.coordinate_scales());
            for i in 0..6 {
                for j in 0..6 {
                    assert_eq!(
                        shifted.native_coefficients()[i][j],
                        law.native_coefficients()[i][j] * coefficient
                    );
                }
            }
            // This exact free twist must stay free as the viscosity changes.
            assert_eq!(
                shifted.evaluate([2.0, -1.0, 7.0, 0.0, 9.0, -4.0]).unwrap(),
                ([0.0; 6], 0.0)
            );
            let (force, heat_rate) = shifted.evaluate([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]).unwrap();
            assert!((heat_rate - coefficient).abs() < 1e-13 * coefficient);
            for i in 0..6 {
                assert!(
                    (force[i] - coefficient * mode[i]).abs() <= 1e-13 * coefficient * mode[i].abs()
                );
            }
            assert!(shifted.reconstruction_residual() <= 128.0 * f64::EPSILON);
        }
        assert_eq!(
            law.scaled(0.0).unwrap().evaluate(mode).unwrap(),
            ([0.0; 6], 0.0)
        );
    }

    #[test]
    fn physical_scaling_matches_original_cross_work_and_preserves_weak_modes() {
        let direction = [1.0, -2.0, 0.5, 3.0, -1.0, 2.0];
        let native = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                direction[i] * direction[j] + if i == j { 2.0 + i as f64 } else { 0.0 }
            })
        });
        let law = PassiveBushingMatrix::try_new(
            native,
            BushingCoordinateScales {
                translation_m: 0.002,
                rotation_rad: 0.3,
            },
        )
        .unwrap();
        let x = [0.02, -0.04, 0.03, 0.3, -0.2, 0.1];
        let base_force: [f64; 6] =
            core::array::from_fn(|i| (0..6).map(|j| native[i][j] * x[j]).sum());
        let base_work: f64 = (0..6).map(|i| x[i] * base_force[i]).sum();
        for coefficient in [1e-60, 0.37, 6.5, 1e60] {
            let (force, work) = law.scaled(coefficient).unwrap().evaluate(x).unwrap();
            for i in 0..6 {
                assert!(
                    (force[i] - coefficient * base_force[i]).abs()
                        <= 2e-12 * (coefficient * base_force[i]).abs()
                );
            }
            assert!((work - coefficient * base_work).abs() <= 2e-12 * coefficient * base_work);
        }
        let tiny = PassiveBushingMatrix::try_new(
            core::array::from_fn(|i| {
                core::array::from_fn(|j| if i == j { f64::MIN_POSITIVE } else { 0.0 })
            }),
            BushingCoordinateScales {
                translation_m: 1.0,
                rotation_rad: 1.0,
            },
        )
        .unwrap()
        .scaled(0.5)
        .unwrap();
        let (force, work) = tiny.evaluate([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]).unwrap();
        assert_eq!(force[0], 0.5 * f64::MIN_POSITIVE);
        assert_eq!(work, 0.5 * f64::MIN_POSITIVE);
    }

    #[test]
    fn physical_scaling_refuses_lost_coefficients_without_changing_reference() {
        let diagonal = |entry| {
            PassiveBushingMatrix::try_new(
                core::array::from_fn(|i| {
                    core::array::from_fn(|j| if i == j { entry } else { 0.0 })
                }),
                BushingCoordinateScales {
                    translation_m: 1.0,
                    rotation_rad: 1.0,
                },
            )
            .unwrap()
        };
        let original = diagonal(1e-200);
        let saved = original.clone();
        for coefficient in [f64::NAN, f64::INFINITY, -1.0, 1e-200] {
            assert!(original.scaled(coefficient).is_err());
            assert_eq!(original, saved);
        }
        assert!(diagonal(1e200).scaled(1e200).is_err());
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
