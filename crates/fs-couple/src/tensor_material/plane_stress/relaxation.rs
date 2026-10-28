//! Exact coordinates for the complete condensed Maxwell energy, not a fit.
use super::*;

/// Exact independent relaxation coordinates of one original plane-stress point.
/// The source still owns normal-strain elimination, its full 3D law and strain
/// limit. Every original history coordinate survives this invertible transform.
/// Its energy is `q^T K_eq q/2 + sum (c_j^T q-z_j)^2/2`, with
/// `z_j' = rate_j (c_j^T q-z_j)`. Rates can differ from the source branch rates:
/// enforcing zero normal traction couples the ORIGINAL complete histories.
#[derive(Clone, Debug)]
pub struct PlaneStressRelaxation {
    point: PlaneStressPoint,
    equilibrium: Vec<f64>,
    rows: Vec<Vec<f64>>,
    rates: Vec<f64>,
    original_from_memory: Vec<f64>,
}

impl PlaneStressPoint {
    /// Express this already-condensed energy in independent passive channels.
    /// This is an exact change of history coordinates, using the shared symmetric
    /// eigensolver. No branches are condensed separately, fitted, clipped or
    /// truncated. The complete energy reconstruction must pass at 1e-9 relative
    /// diagonal scale. Cancellation is polled before/after the bounded dense
    /// eigensolve; interruption within that low-level call is not promised.
    pub fn relaxation_image(
        self,
        max_scalar_terms: usize,
        gate: &CancelGate,
    ) -> Result<PlaneStressRelaxation> {
        poll(gate)?;
        let n = self.coordinates;
        let r = self.memory_coordinates();
        let m = n + r;
        let d = n + 3;
        let work = 512 * (r as u128).pow(3)
            + 8 * (m as u128).pow(3)
            + (n as u128).pow(2) * (d as u128).pow(2);
        if work > max_scalar_terms as u128 {
            return Err(bad(
                "complete plane-stress relaxation transform exceeds scalar-work budget",
            ));
        }
        // The equilibrium normal map comes from C_eq alone. It is used only
        // for K_eq; all transient normal strains still use the complete owner.
        let mut t = vec![vec![0.0; d]; n];
        for i in 0..n {
            poll(gate)?;
            t[i][i] = 1.0;
            let rhs =
                core::array::from_fn(|k| -self.source.equilibrium_stiffness()[(n + k) * d + i]);
            let normal = normal_solve(self.source.equilibrium_stiffness(), d, n, &rhs)?;
            t[i][n..].copy_from_slice(&normal);
        }
        let mut equilibrium = vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                let mut value = 0.0;
                for a in 0..d {
                    for b in 0..d {
                        value += t[i][a] * self.source.equilibrium_stiffness()[a * d + b] * t[j][b];
                    }
                }
                equilibrium[i * n + j] = value;
                equilibrium[j * n + i] = value;
            }
        }
        if !finite(&equilibrium) {
            return Err(bad("plane-stress equilibrium projection overflow"));
        }
        if r == 0 {
            return Ok(PlaneStressRelaxation {
                point: self,
                equilibrium,
                rows: Vec::new(),
                rates: Vec::new(),
                original_from_memory: Vec::new(),
            });
        }
        // H=[q,a]^T[A B; B^T D][q,a]/2, a'=-R(B^T q+D a).
        // S=sqrt(R) D sqrt(R)=U Lambda U^T;
        // z=sqrt(Lambda) U^T R^-1/2 a,
        // c=-B sqrt(R) U Lambda^-1/2.
        // Therefore H=q^T K_eq q/2+|c^T q-z|^2/2 and z'=Lambda(c^T q-z).
        let scale = self.rates.iter().copied().fold(0.0_f64, f64::max);
        let roots: Vec<_> = self.rates.iter().map(|v| det::sqrt(v / scale)).collect();
        if roots.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(bad(
                "plane-stress memory rates are unresolved under common scaling",
            ));
        }
        let matrix: Vec<_> = (0..r * r)
            .map(|k| roots[k / r] * self.energy[(n + k / r) * m + n + k % r] * roots[k % r])
            .collect();
        if !finite(&matrix) {
            return Err(bad("plane-stress rate matrix overflow"));
        }
        poll(gate)?;
        let (values, vectors) = fs_la::eigen::jacobi_eigh(&matrix, r);
        poll(gate)?;
        if values.iter().any(|v| !v.is_finite() || *v <= 0.0) || !finite(&vectors) {
            return Err(bad(
                "plane-stress memory has an unresolved positive relaxation rate",
            ));
        }
        for i in 0..r {
            for j in 0..r {
                let orthogonal: f64 = (0..r)
                    .map(|k| vectors[k * r + i] * vectors[k * r + j])
                    .sum();
                let rebuilt: f64 = (0..r)
                    .map(|k| vectors[i * r + k] * values[k] * vectors[j * r + k])
                    .sum();
                let diagonal = det::sqrt(matrix[i * r + i]) * det::sqrt(matrix[j * r + j]);
                if !diagonal.is_finite()
                    || diagonal <= 0.0
                    || (orthogonal - if i == j { 1.0 } else { 0.0 }).abs() > 1e-10
                    || (rebuilt - matrix[i * r + j]).abs() > 1e-10 * diagonal
                {
                    return Err(bad(
                        "plane-stress relaxation eigensolve failed reconstruction",
                    ));
                }
            }
        }
        let rates: Vec<_> = values.iter().map(|v| scale * v).collect();
        let mut original_from_memory = vec![0.0; r * r];
        let mut memory_from_original = vec![0.0; r * r];
        for i in 0..r {
            for j in 0..r {
                original_from_memory[i * r + j] =
                    roots[i] * vectors[i * r + j] / det::sqrt(values[j]);
                memory_from_original[j * r + i] =
                    det::sqrt(values[j]) * vectors[i * r + j] / roots[i];
            }
        }
        let rows: Vec<Vec<f64>> = (0..r)
            .map(|channel| {
                (0..n)
                    .map(|i| {
                        -(0..r)
                            .map(|j| {
                                self.energy[i * m + n + j] * original_from_memory[j * r + channel]
                            })
                            .sum::<f64>()
                    })
                    .collect()
            })
            .collect();
        if rates
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0 || !v.recip().is_finite())
            || !finite(&original_from_memory)
            || !finite(&memory_from_original)
            || rows.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(bad("plane-stress relaxation coordinates overflow"));
        }
        // Reconstruct the COMPLETE original Q, including normal-history cross
        // terms; the frequency/time comparison is not circular with this gate.
        for i in 0..m {
            poll(gate)?;
            for j in 0..m {
                let mut value = if i < n && j < n {
                    equilibrium[i * n + j]
                } else {
                    0.0
                };
                for channel in 0..r {
                    let a = if i < n {
                        rows[channel][i]
                    } else {
                        -memory_from_original[channel * r + i - n]
                    };
                    let b = if j < n {
                        rows[channel][j]
                    } else {
                        -memory_from_original[channel * r + j - n]
                    };
                    value += a * b;
                }
                let diagonal =
                    det::sqrt(self.energy[i * m + i]) * det::sqrt(self.energy[j * m + j]);
                if !value.is_finite()
                    || !diagonal.is_finite()
                    || (value - self.energy[i * m + j]).abs()
                        > 1e-9 * diagonal.max(f64::MIN_POSITIVE)
                {
                    return Err(bad(
                        "plane-stress channels do not reconstruct the complete original energy",
                    ));
                }
            }
        }
        poll(gate)?;
        Ok(PlaneStressRelaxation {
            point: self,
            equilibrium,
            rows,
            rates,
            original_from_memory,
        })
    }
}

impl PlaneStressRelaxation {
    /// The original full-law point and its normal-traction constraint.
    pub fn point(&self) -> &PlaneStressPoint {
        &self.point
    }
    /// Equilibrium stiffness, with no material memory counted twice.
    pub fn equilibrium_stiffness(&self) -> &[f64] {
        &self.equilibrium
    }
    /// Exact strain rows c_j; zero/dependent channels are retained too.
    pub fn coupling_rows(&self) -> &[Vec<f64>] {
        &self.rows
    }
    /// Exact constrained relaxation rates [1/s], one per retained coordinate.
    pub fn relaxation_rates_per_s(&self) -> &[f64] {
        &self.rates
    }
    /// Recover all source viscous factor coordinates for observation/admission.
    pub fn original_memory(&self, memory: &[f64]) -> Result<Vec<f64>> {
        let r = self.rates.len();
        if memory.len() != r || !finite(memory) {
            return Err(bad("complete finite plane-stress memory required"));
        }
        let original: Vec<_> = self
            .original_from_memory
            .chunks(r.max(1))
            .map(|row| row.iter().zip(memory).map(|(a, z)| a * z).sum())
            .collect();
        if !finite(&original) {
            return Err(bad("original plane-stress memory reconstruction overflow"));
        }
        Ok(original)
    }
    /// Reconstruct the full original strain and enforce its material domain.
    pub fn strain(&self, displacement: &[f64], memory: &[f64]) -> Result<Voigt> {
        self.point
            .strain(displacement, &self.original_memory(memory)?)
    }
    /// Independently recompute normal traction from the original 3D law.
    pub fn normal_stress_pa(&self, displacement: &[f64], memory: &[f64]) -> Result<[f64; 3]> {
        self.point
            .normal_stress_pa(displacement, &self.original_memory(memory)?)
    }
}
