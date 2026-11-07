//! Columnwise tensor Kirchhoff extension of the native Cartesian HMM leaf.
//!
//! Let M(K) be the existing linear local face matrix and B its row sums.
//! Since M is linear in K, the cell-outward heat is
//! q_f = -sum_g M_fg(mean_K(Tc,Tg)) (Tg - Tc).
//! The fundamental theorem gives dq_f/dTc = B_f(K(Tc)) and
//! dq_f/dTg = -M_fg(K(Tg)); no K'(T) term is inserted again.
//!
//! This is the exact scalar Kirchhoff transform for K(T)=k(T)K0 and the
//! native integrated half-cell law for independent grid-diagonal K(T).
//! For arbitrary full K(T), the column means generally differ. The resulting
//! Jacobian is not symmetric and the linear HMM coercivity theorem does not
//! supply a discrete entropy/maximum principle or unconditional nonlinear
//! convergence. Smooth-law local flux consistency is tested under refinement.
//! The solve owner must retain bounded iteration and original heat residuals.
//!
//! This numerical leaf owns no material, source, interface or time history.

#[derive(Debug, Clone, Copy)]
pub(super) struct KirchhoffCell {
    widths: [f64; 3],
    areas: [f64; 3],
}

#[derive(Debug, Clone, Copy)]
pub(super) struct KirchhoffHeat {
    /// Physical cell-outward heat in W.
    pub outward_w: [f64; 6],
    /// Arithmetic scale used by the native physical heat-roundoff gate.
    pub arithmetic_w: [f64; 6],
}

#[derive(Debug, Clone, Copy)]
pub(super) struct KirchhoffTangent {
    /// dq_f/dTc, positive, in W/K.
    pub cell: [f64; 6],
    /// dq_f/dTg in W/K, stored as row f, column g.
    pub face: [[f64; 6]; 6],
}

/// Exact cell recovery retained after the local matrix has been assembled.
/// Only one denominator, six face links and one original residual are stored.
#[derive(Debug, Clone, Copy)]
pub(super) struct KirchhoffRecovery {
    diagonal: f64,
    face_link: [f64; 6],
    cell_residual_w: f64,
}

impl KirchhoffRecovery {
    /// Recover the free cell's Newton increment from its six trace increments.
    pub fn cell_increment(&self, face_increment: [f64; 6]) -> Option<f64> {
        if face_increment.iter().any(|v| !v.is_finite()) {
            return None;
        }
        // Scale into Kelvin before the reduction, avoiding an unnecessary
        // large intermediate heat when the final increment is representable.
        let mut increment = -self.cell_residual_w / self.diagonal;
        for (link, delta) in self.face_link.into_iter().zip(face_increment) {
            increment = (link / self.diagonal).mul_add(delta, increment);
        }
        increment.is_finite().then_some(increment)
    }
}

/// Transient local Schur data; the spatial owner assembles and drops the
/// dense matrix immediately, retaining only its small recovery record.
#[derive(Debug, Clone, Copy)]
pub(super) struct KirchhoffCondensation {
    pub matrix: [[f64; 6]; 6],
    pub load: [f64; 6],
    recovery: KirchhoffRecovery,
}

impl KirchhoffCondensation {
    pub fn recovery(&self) -> KirchhoffRecovery {
        self.recovery
    }

    pub fn cell_increment(&self, face_increment: [f64; 6]) -> Option<f64> {
        self.recovery.cell_increment(face_increment)
    }
}

impl KirchhoffTangent {
    /// Eliminate one free cell from the exact original Newton equations.
    ///
    /// With positive endpoint face links b_g=-dq_g/dTg and
    /// D=sum_f dq_f/dTc + added_diagonal, the sign-reversed face Schur
    /// matrix is -dq_f/dTg-(dq_f/dTc)b_g/D. Its residual correction is
    /// -(dq_f/dTc)Rc/D. The added diagonal is the original nonnegative
    /// sink/storage tangent, not a damping or conductivity approximation.
    pub fn condense(
        self,
        added_diagonal: f64,
        cell_residual_w: f64,
    ) -> Option<KirchhoffCondensation> {
        if !added_diagonal.is_finite()
            || added_diagonal < 0.0
            || !cell_residual_w.is_finite()
            || self.cell.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || self.face.iter().flatten().any(|v| !v.is_finite())
        {
            return None;
        }
        let diagonal = self.cell.iter().sum::<f64>() + added_diagonal;
        let face_link: [f64; 6] = std::array::from_fn(|g| -self.face[g][g]);
        if !diagonal.is_finite()
            || diagonal <= 0.0
            || face_link.iter().any(|v| !v.is_finite() || *v <= 0.0)
        {
            return None;
        }
        let mut matrix = [[0.0; 6]; 6];
        let mut load = [0.0; 6];
        for f in 0..6 {
            let weight = self.cell[f] / diagonal;
            load[f] = -weight * cell_residual_w;
            for (g, link) in face_link.into_iter().enumerate() {
                matrix[f][g] = (-weight).mul_add(link, -self.face[f][g]);
            }
        }
        if matrix.iter().flatten().chain(&load).any(|v| !v.is_finite()) {
            return None;
        }
        Some(KirchhoffCondensation {
            matrix,
            load,
            recovery: KirchhoffRecovery {
                diagonal,
                face_link,
                cell_residual_w,
            },
        })
    }
}

impl KirchhoffCell {
    pub fn new(widths: [f64; 3]) -> Option<Self> {
        if widths.iter().any(|h| !h.is_finite() || *h <= 0.0) {
            return None;
        }
        let volume = widths[0] * widths[1] * widths[2];
        let areas = std::array::from_fn(|a| widths[(a + 1) % 3] * widths[(a + 2) % 3]);
        if !volume.is_finite() || volume <= 0.0 || areas.iter().any(|a| !a.is_finite() || *a <= 0.0)
        {
            return None;
        }
        Some(Self { widths, areas })
    }

    /// One HMM column; only the needed entries are formed. The material owner
    /// has already admitted each full symmetric positive-definite tensor.
    fn column(&self, k: [[f64; 3]; 3], g: usize) -> Option<[f64; 6]> {
        let b = g / 2;
        let half = self.areas[b] * (k[b][b] / (0.5 * self.widths[b]));
        if !half.is_finite() || half <= 0.0 {
            return None;
        }
        let mut column = [0.0; 6];
        column[g] = half;
        for (f, value) in column.iter_mut().enumerate() {
            let a = f / 2;
            if a == b {
                continue;
            }
            let sign = if (f & 1) == (g & 1) { 1.0 } else { -1.0 };
            *value = sign * self.widths[3 - a - b] * k[a][b];
            if !value.is_finite() {
                return None;
            }
        }
        Some(column)
    }

    pub fn heat(
        &self,
        cell: f64,
        face: [f64; 6],
        mean: [[[f64; 3]; 3]; 6],
    ) -> Option<KirchhoffHeat> {
        if !cell.is_finite() || face.iter().any(|t| !t.is_finite()) {
            return None;
        }
        let mut result = KirchhoffHeat {
            outward_w: [0.0; 6],
            arithmetic_w: [0.0; 6],
        };
        for g in 0..6 {
            let column = self.column(mean[g], g)?;
            let difference = face[g] - cell;
            let kelvin_scale = cell.abs() + face[g].abs();
            if !difference.is_finite() || !kelvin_scale.is_finite() {
                return None;
            }
            for (f, coefficient) in column.into_iter().enumerate() {
                result.outward_w[f] = (-coefficient).mul_add(difference, result.outward_w[f]);
                result.arithmetic_w[f] += coefficient.abs() * kelvin_scale;
            }
        }
        if result
            .outward_w
            .iter()
            .chain(&result.arithmetic_w)
            .any(|v| !v.is_finite())
        {
            return None;
        }
        Some(result)
    }

    pub fn tangent(
        &self,
        cell: [[f64; 3]; 3],
        face: [[[f64; 3]; 3]; 6],
    ) -> Option<KirchhoffTangent> {
        let mut result = KirchhoffTangent {
            cell: [0.0; 6],
            face: [[0.0; 6]; 6],
        };
        for (g, tensor) in face.into_iter().enumerate() {
            let column = self.column(tensor, g)?;
            let a = g / 2;
            let half = self.areas[a] * (cell[a][a] / (0.5 * self.widths[a]));
            if !half.is_finite() || half <= 0.0 {
                return None;
            }
            result.cell[g] = half;
            for (f, coefficient) in column.into_iter().enumerate() {
                result.face[f][g] = -coefficient;
            }
        }
        if !result.cell.iter().sum::<f64>().is_finite() {
            return None;
        }
        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::local::LocalOperator;
    use super::*;

    const BASE: [[f64; 3]; 3] = [[5.0, 2.0, -1.0], [2.0, 4.0, 0.5], [-1.0, 0.5, 3.0]];

    fn rotate(principal: [f64; 3]) -> [[f64; 3]; 3] {
        let r = [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]];
        std::array::from_fn(|i| {
            std::array::from_fn(|j| (0..3).map(|a| r[i][a] * principal[a] * r[j][a]).sum())
        })
    }

    // An original SPD polynomial law with independently varying eigenvalues.
    fn law(t: f64) -> [[f64; 3]; 3] {
        rotate([1.0 + 0.01 * t, 2.0 + 0.00002 * t * t, 3.0 + 0.003 * t])
    }

    fn mean(from: f64, to: f64) -> [[f64; 3]; 3] {
        rotate([
            1.0 + 0.005 * (from + to),
            2.0 + 0.00002 / 3.0 * (from * from + from * to + to * to),
            3.0 + 0.0015 * (from + to),
        ])
    }

    #[test]
    fn constant_full_tensor_reduces_to_existing_native_hmm() {
        let widths = [0.2, 0.7, 1.3];
        let old = LocalOperator::new(BASE, widths).unwrap();
        let new = KirchhoffCell::new(widths).unwrap();
        let cell = 350.0;
        let faces = [311.0, 402.0, 355.0, 307.0, 381.0, 369.0];
        let heat = new.heat(cell, faces, [BASE; 6]).unwrap();
        let old_heat = old.outward_heat(cell, faces);
        let old_arithmetic = old.heat_arithmetic(cell, faces);
        for f in 0..6 {
            assert!((heat.outward_w[f] - old_heat[f]).abs() < 3e-12);
            assert_eq!(heat.arithmetic_w[f], old_arithmetic[f]);
        }
        let tangent = new.tangent(BASE, [BASE; 6]).unwrap();
        assert_eq!(tangent.cell, old.cell_link);
        for f in 0..6 {
            for g in 0..6 {
                assert_eq!(tangent.face[f][g], -old.matrix[f][g]);
            }
        }
    }

    #[test]
    fn exact_endpoint_tangent_matches_original_nonlinear_residual() {
        let op = KirchhoffCell::new([0.3, 0.7, 1.1]).unwrap();
        let temperature = [342.0, 310.0, 415.0, 389.0, 291.0, 448.0, 374.0];
        let evaluate = |t: [f64; 7]| {
            let face = std::array::from_fn(|g| t[g + 1]);
            op.heat(t[0], face, face.map(|v| mean(t[0], v)))
                .unwrap()
                .outward_w
        };
        let faces = std::array::from_fn(|g| temperature[g + 1]);
        let tangent = op.tangent(law(temperature[0]), faces.map(law)).unwrap();
        let epsilon = 0.001;
        for g in 0..7 {
            let mut plus = temperature;
            let mut minus = temperature;
            plus[g] += epsilon;
            minus[g] -= epsilon;
            let qp = evaluate(plus);
            let qm = evaluate(minus);
            for f in 0..6 {
                let finite_difference = (qp[f] - qm[f]) / (2.0 * epsilon);
                let exact = if g == 0 {
                    tangent.cell[f]
                } else {
                    tangent.face[f][g - 1]
                };
                assert!(
                    (finite_difference - exact).abs() < 2e-8 * exact.abs().max(1.0),
                    "row {f}, unknown {g}: {finite_difference} vs {exact}",
                );
            }
        }
        // Differing endpoint laws produce a nonsymmetric Newton matrix.
        assert!((tangent.face[0][2] - tangent.face[2][0]).abs() > 0.01);
    }

    #[test]
    fn nonlinear_common_factor_transformed_affine_patch_is_exact() {
        // K(T)=(T/300) BASE, with scalar potential U(T)=T^2/600.
        let widths = [0.3, 0.7, 1.1];
        let gradient_u = [0.0, 70.0, -20.0];
        let cell: f64 = 350.0;
        let potential = cell * cell / 600.0;
        let faces: [f64; 6] = std::array::from_fn(|f| {
            let a = f / 2;
            let sign = if f & 1 == 0 { -1.0 } else { 1.0 };
            (600.0 * (potential + 0.5 * sign * widths[a] * gradient_u[a])).sqrt()
        });
        let means = faces.map(|t| BASE.map(|row| row.map(|k| k * (cell + t) / 600.0)));
        let heat = KirchhoffCell::new(widths)
            .unwrap()
            .heat(cell, faces, means)
            .unwrap();
        for f in 0..6 {
            let a = f / 2;
            let sign = if f & 1 == 0 { -1.0 } else { 1.0 };
            let expected = -sign
                * widths[(a + 1) % 3]
                * widths[(a + 2) % 3]
                * BASE[a]
                    .iter()
                    .zip(gradient_u)
                    .map(|(k, g)| k * g)
                    .sum::<f64>();
            assert!((heat.outward_w[f] - expected).abs() < 5e-12);
        }
        assert!(heat.outward_w[0].abs() > 1.0);
        assert!((faces[0] - cell).abs() < 1e-12);
        assert!(heat.outward_w.iter().sum::<f64>().abs() < 8e-12);
    }

    #[test]
    fn nonlinear_rotated_flux_matches_manufactured_field_under_refinement() {
        // T=340+30x-20y+10z+8xy+6z^2, sampled about the origin.
        // This test compares to the continuum Fourier flux, not to another
        // discrete operator. It establishes smooth local flux consistency,
        // not convergence of a whole nonlinear boundary-value solve.
        let field = |x: [f64; 3]| {
            340.0 + 30.0 * x[0] - 20.0 * x[1] + 10.0 * x[2] + 8.0 * x[0] * x[1] + 6.0 * x[2] * x[2]
        };
        let gradient = |x: [f64; 3]| [30.0 + 8.0 * x[1], -20.0 + 8.0 * x[0], 10.0 + 12.0 * x[2]];
        let mut previous = f64::INFINITY;
        let mut errors = Vec::new();
        for h in [0.4, 0.2, 0.1, 0.05] {
            let widths = [h, 0.8 * h, 1.2 * h];
            let coordinates: [[f64; 3]; 6] = std::array::from_fn(|f| {
                let mut x = [0.0; 3];
                x[f / 2] = if f & 1 == 0 { -0.5 } else { 0.5 } * widths[f / 2];
                x
            });
            let cell = field([0.0; 3]);
            let faces = coordinates.map(field);
            let heat = KirchhoffCell::new(widths)
                .unwrap()
                .heat(cell, faces, faces.map(|t| mean(cell, t)))
                .unwrap();
            let mut error: f64 = 0.0;
            for f in 0..6 {
                let a = f / 2;
                let sign = if f & 1 == 0 { -1.0 } else { 1.0 };
                let grad = gradient(coordinates[f]);
                let exact = -sign
                    * law(faces[f])[a]
                        .iter()
                        .zip(grad)
                        .map(|(k, g)| k * g)
                        .sum::<f64>();
                let numerical = heat.outward_w[f] / (widths[(a + 1) % 3] * widths[(a + 2) % 3]);
                error = error.max((numerical - exact).abs());
            }
            assert!(
                error < 0.55 * previous,
                "refinement errors {errors:?}, next {error}"
            );
            previous = error;
            errors.push(error);
        }
        assert!(errors[3] < 1.0, "{errors:?}");
    }

    #[test]
    fn diagonal_nonlinear_links_integrate_each_direction_without_cross_terms() {
        let widths = [0.3, 0.7, 1.1];
        let cell = 325.0;
        let faces = [250.0, 450.0, 270.0, 410.0, 299.0, 501.0];
        let diagonal_mean = |t| {
            let average = [2.0 + 0.005 * (cell + t), 3.0 + 0.01 * (cell + t), 4.0];
            [
                [average[0], 0.0, 0.0],
                [0.0, average[1], 0.0],
                [0.0, 0.0, average[2]],
            ]
        };
        let heat = KirchhoffCell::new(widths)
            .unwrap()
            .heat(cell, faces, faces.map(diagonal_mean))
            .unwrap();
        for f in 0..6 {
            let a = f / 2;
            let (base, slope) = [(2.0, 0.01), (3.0, 0.02), (4.0, 0.0)][a];
            let integral =
                base * (faces[f] - cell) + 0.5 * slope * (faces[f] * faces[f] - cell * cell);
            let expected =
                -widths[(a + 1) % 3] * widths[(a + 2) % 3] / (0.5 * widths[a]) * integral;
            assert!((heat.outward_w[f] - expected).abs() < 4e-12);
        }
    }

    #[test]
    fn refuses_unrepresentable_geometry_heat_and_tangents() {
        for widths in [[0.0, 1.0, 1.0], [f64::NAN, 1.0, 1.0], [1e200; 3]] {
            assert!(KirchhoffCell::new(widths).is_none());
        }
        let op = KirchhoffCell::new([1.0; 3]).unwrap();
        assert!(op.heat(f64::NAN, [300.0; 6], [BASE; 6]).is_none());
        assert!(op.heat(f64::MAX, [f64::MAX; 6], [BASE; 6]).is_none());
        let mut bad = BASE;
        bad[0][0] = f64::INFINITY;
        assert!(op.tangent(bad, [BASE; 6]).is_none());
        assert!(op.tangent(BASE, [bad; 6]).is_none());
    }

    // Independent partial-pivoting solve, used only to compare the six-trace
    // Schur equation to the complete seven-variable Newton equation.
    fn dense_solve<const N: usize>(mut a: [[f64; N]; N], mut b: [f64; N]) -> [f64; N] {
        for column in 0..N {
            let pivot = (column..N)
                .max_by(|&i, &j| a[i][column].abs().total_cmp(&a[j][column].abs()))
                .unwrap();
            assert!(a[pivot][column].abs() > 1e-14);
            a.swap(column, pivot);
            b.swap(column, pivot);
            for row in column + 1..N {
                let factor = a[row][column] / a[column][column];
                for j in column + 1..N {
                    a[row][j] -= factor * a[column][j];
                }
                b[row] -= factor * b[column];
            }
        }
        let mut x = [0.0; N];
        for row in (0..N).rev() {
            let mut rhs = b[row];
            for j in row + 1..N {
                rhs -= a[row][j] * x[j];
            }
            x[row] = rhs / a[row][row];
        }
        x
    }

    #[test]
    fn exact_nonlinear_schur_and_recovery_match_the_full_seven_equations() {
        let cell = KirchhoffCell::new([0.3, 0.7, 1.1]).unwrap();
        let face_t = [310.0, 415.0, 389.0, 291.0, 448.0, 374.0];
        let tangent = cell.tangent(law(342.0), face_t.map(law)).unwrap();
        let residual = 17.25;
        let face_rhs = [11.0, -8.0, 6.0, 3.0, -4.0, 13.0];
        // Boundary conductances independently anchor the face equation.
        let boundary = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
        for added in [0.0, 3.5, 1200.0] {
            let reduced = tangent.condense(added, residual).unwrap();
            let mut full = [[0.0; 7]; 7];
            let mut rhs = [0.0; 7];
            full[0][0] = tangent.cell.iter().sum::<f64>() + added;
            rhs[0] = -residual;
            for f in 0..6 {
                full[0][f + 1] = tangent.face[f][f];
                full[f + 1][0] = -tangent.cell[f];
                rhs[f + 1] = face_rhs[f];
                for g in 0..6 {
                    full[f + 1][g + 1] = -tangent.face[f][g];
                }
                full[f + 1][f + 1] += boundary[f];
            }
            let expected = dense_solve(full, rhs);
            let mut matrix = reduced.matrix;
            let mut load = reduced.load;
            for f in 0..6 {
                matrix[f][f] += boundary[f];
                load[f] += face_rhs[f];
            }
            let face_increment = dense_solve(matrix, load);
            let recovery = reduced.recovery();
            let cell_increment = recovery.cell_increment(face_increment).unwrap();
            assert_eq!(Some(cell_increment), reduced.cell_increment(face_increment));
            assert!((cell_increment - expected[0]).abs() < 2e-12);
            for f in 0..6 {
                assert!((face_increment[f] - expected[f + 1]).abs() < 2e-12);
            }
            let actual: [f64; 7] = std::array::from_fn(|i| {
                if i == 0 {
                    cell_increment
                } else {
                    face_increment[i - 1]
                }
            });
            for i in 0..7 {
                let defect = full[i].iter().zip(actual).map(|(a, x)| a * x).sum::<f64>() - rhs[i];
                assert!(defect.abs() < 3e-11, "row {i}, added {added}: {defect}");
            }
            // Independent nonlinear endpoint tensors retain the nonsymmetric
            // Jacobian; condensation must not mirror or symmetrize it.
            assert!((reduced.matrix[0][2] - reduced.matrix[2][0]).abs() > 0.01);
        }
    }

    #[test]
    fn constant_schur_reduces_to_existing_storage_condensation() {
        let widths = [0.3, 0.7, 1.1];
        let native = LocalOperator::new(BASE, widths).unwrap();
        let tangent = KirchhoffCell::new(widths)
            .unwrap()
            .tangent(BASE, [BASE; 6])
            .unwrap();
        for added in [0.0, 4.0, 900.0] {
            let reduced = tangent.condense(added, 0.0).unwrap();
            let expected = native.schur(added);
            for f in 0..6 {
                for g in 0..6 {
                    assert!((reduced.matrix[f][g] - expected[f][g]).abs() < 2e-14);
                }
                assert_eq!(reduced.load[f], 0.0);
            }
            let recovered = reduced.recovery().cell_increment([2.0; 6]).unwrap();
            let expected = 2.0 * native.cell_diagonal / (native.cell_diagonal + added);
            assert!((recovered - expected).abs() < 4e-15);
        }
    }

    #[test]
    fn condensation_and_recovery_refuse_invalid_or_overflowing_arithmetic() {
        let tangent = KirchhoffCell::new([1.0; 3])
            .unwrap()
            .tangent(BASE, [BASE; 6])
            .unwrap();
        for added in [-1.0, f64::INFINITY, f64::NAN] {
            assert!(tangent.condense(added, 1.0).is_none());
        }
        assert!(tangent.condense(0.0, f64::INFINITY).is_none());
        let mut bad = tangent;
        bad.cell[0] = f64::NAN;
        assert!(bad.condense(0.0, 0.0).is_none());
        bad = tangent;
        bad.cell = [f64::MAX; 6];
        assert!(bad.condense(0.0, 0.0).is_none());
        bad = tangent;
        bad.face[0][0] = 0.0;
        assert!(bad.condense(0.0, 0.0).is_none());
        bad = tangent;
        bad.face[0][2] = f64::NEG_INFINITY;
        assert!(bad.condense(0.0, 0.0).is_none());
        let recovery = tangent.condense(0.0, 1.0).unwrap().recovery();
        assert!(recovery.cell_increment([f64::NAN; 6]).is_none());
        assert!(recovery.cell_increment([f64::INFINITY; 6]).is_none());
        let extreme = KirchhoffRecovery {
            diagonal: 1.0,
            face_link: [f64::MAX; 6],
            cell_residual_w: 0.0,
        };
        assert!(extreme.cell_increment([2.0; 6]).is_none());
    }
}
