//! Cartesian specialization of the hybrid finite-volume energy in
//! Droniou et al., M3AS 20 (2010), equations 2.17--2.21,
//! doi:10.1142/S0218202510004222. Stabilization alpha_f = n_f^T K n_f.
//!
//! This numerical leaf has no material, boundary, time or mesh ownership.

/// Refuse non-finite, nonsymmetric or numerically non-positive tensors.
/// Scaling and LDL pivots avoid a determinant test and arbitrary eigenvalue
/// floors. No off-diagonal component is projected away or symmetrized.
pub(super) fn is_spd(k: [[f64; 3]; 3]) -> bool {
    if k.iter().flatten().any(|v| !v.is_finite())
        || (0..3).any(|a| (0..3).any(|b| k[a][b] != k[b][a]))
    {
        return false;
    }
    let scale = k[0][0].max(k[1][1]).max(k[2][2]);
    if scale <= 0.0 { return false; }
    let a = k.map(|row| row.map(|v| v / scale));
    let d0 = a[0][0];
    if d0 <= 0.0 { return false; }
    let l10 = a[1][0] / d0;
    let l20 = a[2][0] / d0;
    let d1 = (-l10).mul_add(a[1][0], a[1][1]);
    let off = (-l20).mul_add(a[1][0], a[2][1]);
    if !d1.is_finite() || d1 <= 0.0 { return false; }
    let d2 = (-off / d1).mul_add(off, (-l20).mul_add(a[2][0], a[2][2]));
    d2.is_finite() && d2 > 0.0
}

/// Half-cell matrix on the six face-minus-cell temperatures. Its energy is
/// V grad(T)^T K grad(T) + sum_a 4 V K_aa/h_a^2 * r_a^2, where
/// grad_a=(T_a+ - T_a-)/h_a and r_a=(T_a+ + T_a-)/2-T_cell.
/// Constant offsets disappear; affine fields have zero stabilization error.
#[derive(Debug, Clone, Copy)]
pub(super) struct LocalOperator {
    pub matrix: [[f64; 6]; 6],
    pub cell_link: [f64; 6],
    pub cell_diagonal: f64,
}

impl LocalOperator {
    pub fn new(k: [[f64; 3]; 3], widths: [f64; 3]) -> Option<Self> {
        if widths.iter().any(|h| !h.is_finite() || *h <= 0.0) { return None; }
        let volume = widths[0] * widths[1] * widths[2];
        if !volume.is_finite() || volume <= 0.0 { return None; }
        let mut matrix = [[0.0; 6]; 6];
        let mut cell_link = [0.0; 6];
        for f in 0..6 {
            let a = f / 2;
            let half = widths[(a + 1) % 3] * widths[(a + 2) % 3]
                * (k[a][a] / (0.5 * widths[a]));
            if !half.is_finite() || half <= 0.0 { return None; }
            cell_link[f] = half;
            matrix[f][f] = half;
            for g in 0..f {
                let b = g / 2;
                if a == b { continue; }
                // V/(h_a h_b) is the remaining width. This form avoids
                // unnecessary products/quotients on graded box grids.
                let sign = if (f & 1) == (g & 1) { 1.0 } else { -1.0 };
                let value = sign * widths[3 - a - b] * k[a][b];
                if !value.is_finite() { return None; }
                matrix[f][g] = value;
                matrix[g][f] = value;
            }
        }
        let cell_diagonal: f64 = cell_link.iter().sum();
        if !cell_diagonal.is_finite() { return None; }
        Some(Self { matrix, cell_link, cell_diagonal })
    }

    /// Remove the cell unknown with an optional positive sink/storage
    /// diagonal. The caller owns the corresponding source/reservoir load.
    pub fn schur(&self, added_diagonal: f64) -> [[f64; 6]; 6] {
        let diagonal = self.cell_diagonal + added_diagonal;
        let mut reduced = self.matrix;
        for f in 0..6 {
            for g in 0..=f {
                let value = self.matrix[f][g]
                    - self.cell_link[f] / diagonal * self.cell_link[g];
                reduced[f][g] = value;
                reduced[g][f] = value;
            }
        }
        reduced
    }

    /// Physical cell-outward heat in W, evaluated at the returned endpoint.
    /// Centering before multiplication avoids cancelling a Kelvin offset.
    pub fn outward_heat(&self, cell: f64, face: [f64; 6]) -> [f64; 6] {
        let difference = face.map(|t| t - cell);
        self.matrix.map(|row| {
            -row.iter().zip(difference).fold(0.0, |sum, (a, d)| a.mul_add(d, sum))
        })
    }

    /// Explicit arithmetic scale for the native local heat-roundoff gate.
    pub fn heat_arithmetic(&self, cell: f64, face: [f64; 6]) -> [f64; 6] {
        self.matrix.map(|row| {
            row.iter().zip(face).map(|(a, t)| a.abs() * (cell.abs() + t.abs())).sum()
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rotated_affine_patch_keeps_transverse_heat_at_equal_cell_temperatures() {
        let k = [[5.0, 2.0, -1.0], [2.0, 4.0, 0.5], [-1.0, 0.5, 3.0]];
        assert!(is_spd(k));
        let widths = [0.2, 0.7, 1.3];
        let gradient = [0.0, 7.0, -2.0];
        let op = LocalOperator::new(k, widths).unwrap();
        let faces = core::array::from_fn(|f| {
            350.0 + if f & 1 == 0 { -0.5 } else { 0.5 }
                * widths[f / 2] * gradient[f / 2]
        });
        let heat = op.outward_heat(350.0, faces);
        for f in 0..6 {
            let a = f / 2;
            let expected = -if f & 1 == 0 { -1.0 } else { 1.0 }
                * widths[(a + 1) % 3] * widths[(a + 2) % 3]
                * k[a].iter().zip(gradient).map(|(v, g)| v * g).sum::<f64>();
            assert!((heat[f] - expected).abs() < 2e-12, "{heat:?}");
        }
        assert!(heat[0].abs() > 1.0); // grad_x=0 does not imply q_x=0.
        assert!(heat.iter().sum::<f64>().abs() < 2e-12);
    }

    #[test]
    fn diagonal_reduces_exactly_to_half_cell_resistors_and_contact_jump() {
        let a = LocalOperator::new([[2.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,4.0]],
            [0.3,0.4,0.5]).unwrap();
        let b = LocalOperator::new([[7.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,4.0]],
            [0.8,0.4,0.5]).unwrap();
        let area = 0.4 * 0.5;
        let resistance = 0.17;
        let heat = area * (400.0 - 300.0) / (0.15/2.0 + resistance + 0.4/7.0);
        let ta = 400.0 - heat/a.cell_link[1];
        let tb = 300.0 + heat/b.cell_link[0];
        let mut fa = [400.0;6]; fa[1] = ta;
        let mut fb = [300.0;6]; fb[0] = tb;
        assert!((a.outward_heat(400.0, fa)[1] - heat).abs() < 1e-12);
        assert!((b.outward_heat(300.0, fb)[0] + heat).abs() < 1e-12);
        assert!((area / resistance * (ta-tb) - heat).abs() < 1e-12);
        for f in 0..6 {
            for g in 0..6 {
                if f != g { assert_eq!(a.matrix[f][g], 0.0); }
            }
        }
    }

    #[test]
    fn condensation_preserves_rotated_energy_and_cell_source_equation() {
        let op = LocalOperator::new([[4.0,1.0,0.5],[1.0,3.0,-0.2],[0.5,-0.2,2.0]],
            [0.2,0.5,0.9]).unwrap();
        let face = [-3.0, 8.0, 2.0, -5.0, 1.0, 7.0];
        let added = 11.0;
        let load = 9.0;
        let diagonal = op.cell_diagonal + added;
        let cell = (load + op.cell_link.iter().zip(face).map(|(b,t)| b*t).sum::<f64>())
            / diagonal;
        let heat = op.outward_heat(cell, face);
        assert!((heat.iter().sum::<f64>() + added*cell-load).abs() < 1e-12);
        let reduced = op.schur(added);
        for f in 0..6 {
            let residual = reduced[f].iter().zip(face).map(|(a,t)| a*t).sum::<f64>()
                - op.cell_link[f]*load/diagonal;
            assert!((residual + heat[f]).abs() < 1e-12);
        }
        let energy = face.iter().enumerate().map(|(f,t)| t *
            reduced[f].iter().zip(face).map(|(a,u)| a*u).sum::<f64>()).sum::<f64>();
        assert!(energy > 0.0);
        let null = op.schur(0.0);
        for row in null { assert!(row.iter().sum::<f64>().abs() < 1e-14); }
    }

    #[test]
    fn tensor_admission_refuses_indefinite_asymmetric_and_unrepresentable_geometry() {
        assert!(!is_spd([[1.0,2.0,0.0],[2.0,1.0,0.0],[0.0,0.0,1.0]]));
        assert!(!is_spd([[1.0,0.1,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]));
        assert!(!is_spd([[1.0,0.0,0.0],[0.0,f64::NAN,0.0],[0.0,0.0,1.0]]));
        let k = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
        assert!(LocalOperator::new(k, [1e-300;3]).is_none());
        assert!(LocalOperator::new(k, [1e300;3]).is_none());
    }
}
