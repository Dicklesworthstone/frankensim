//! Orthonormal acoustic ports in the complete mass-normalized board basis.
//!
//! Rows are Q^T. The acoustic velocity is Q^T v and its reaction is Q f, so
//! mechanical power is unchanged. Orthonormality also makes acoustic-port
//! kinetic energy the exact selected part of the full board kinetic energy.
//! This is a coordinate map; only comparison with the complete physical load
//! can establish whether a truncated acoustic space is adequate.
use fs_math::c64::C64;

#[derive(Clone, Debug)]
pub struct PortBasis {
    board_ports: usize,
    vectors: Vec<Vec<f64>>,
    pub(super) orthogonality_error: f64,
}

impl PortBasis {
    /// Admit a supplied basis without normalizing, deleting or reordering rows.
    pub fn new(board_ports: usize, vectors: Vec<Vec<f64>>) -> Result<Self, String> {
        let ports = vectors.len();
        if !(1..=crate::linear::MAX_BOARD_MODES).contains(&board_ports)
            || !(1..=32).contains(&ports) || ports > board_ports
            || vectors.iter().any(|row| row.len() != board_ports
                || row.iter().any(|v| !v.is_finite())) {
            return Err("radiation basis requires 1..32 complete rows over 1..128 board coordinates".into());
        }
        let tolerance = 256. * f64::EPSILON * (board_ports + ports + 1) as f64;
        let mut orthogonality_error = 0.0_f64;
        for (i, a) in vectors.iter().enumerate() {
            for (j, b) in vectors.iter().enumerate().take(i + 1) {
                let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
                let error = (dot - if i == j { 1. } else { 0. }).abs();
                if !error.is_finite() || error > tolerance {
                    return Err("radiation basis is not orthonormal in the loaded board energy coordinates".into());
                }
                orthogonality_error = orthogonality_error.max(error);
            }
        }
        Ok(Self { board_ports, vectors, orthogonality_error })
    }

    pub fn board_ports(&self) -> usize { self.board_ports }
    pub fn ports(&self) -> usize { self.vectors.len() }
    pub fn vectors(&self) -> &[Vec<f64>] { &self.vectors }

    /// Q^T Z Q, retaining signed complex cross-port terms.
    pub fn project_impedance(&self, impedance: &[C64]) -> Result<Vec<C64>, String> {
        let n = self.board_ports;
        let r = self.ports();
        self.check_matrix(impedance, n)?;
        let mut right = vec![C64::ZERO; n * r];
        for i in 0..n {
            for (a, row) in self.vectors.iter().enumerate() {
                for (j, &q) in row.iter().enumerate() {
                    right[i * r + a] = right[i * r + a] + impedance[i * n + j].scale(q);
                }
            }
        }
        let mut result = vec![C64::ZERO; r * r];
        for (a, row) in self.vectors.iter().enumerate() {
            for b in 0..r {
                for (i, &q) in row.iter().enumerate() {
                    result[a * r + b] = result[a * r + b] + right[i * r + b].scale(q);
                }
            }
        }
        self.check_matrix(&result, r)?;
        Ok(result)
    }

    /// Q Z Q^T in the original complete board coordinates. This lift preserves
    /// reciprocity and positive-realness, but cannot restore an omitted field.
    pub fn lift_impedance(&self, impedance: &[C64]) -> Result<Vec<C64>, String> {
        let n = self.board_ports;
        let r = self.ports();
        self.check_matrix(impedance, r)?;
        let mut left = vec![C64::ZERO; n * r];
        for i in 0..n {
            for b in 0..r {
                for (a, row) in self.vectors.iter().enumerate() {
                    left[i * r + b] = left[i * r + b] + impedance[a * r + b].scale(row[i]);
                }
            }
        }
        let mut result = vec![C64::ZERO; n * n];
        for i in 0..n {
            for j in 0..n {
                for (b, row) in self.vectors.iter().enumerate() {
                    result[i * n + j] = result[i * n + j] + left[i * r + b].scale(row[j]);
                }
            }
        }
        self.check_matrix(&result, n)?;
        Ok(result)
    }

    fn check_matrix(&self, matrix: &[C64], n: usize) -> Result<(), String> {
        if matrix.len() != n * n
            || matrix.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
            return Err("radiation projection requires a complete finite impedance matrix".into());
        }
        Ok(())
    }
}
