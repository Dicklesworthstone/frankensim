//! Finite reflected waves without an intermediate input-impedance inverse.

use super::*;

/// Staircase input reflection in the inlet station's native modal waves.
///
/// Pressure and volume flow satisfy `P = (I + R) a` and
/// `U = diag(1/Zc) (I - R) a`. The response is retained before the impedance
/// conversion, so an inlet flow node (`R` with eigenvalue +1) need not fail a
/// reflection query. Native evanescent waves are not power-normalized; use
/// [`Self::at_references`] for a real-reference scattering port.
#[derive(Debug, Clone)]
pub struct ModalReflectionResponse {
    /// Angular frequency [rad/s].
    pub omega: f64,
    /// Retained radial mode count, including the plane mode.
    pub n_modes: usize,
    /// Native pressure-wave reflection matrix, row-major N x N (`b = R a`).
    pub reflection_matrix: Vec<C64>,
    /// Native per-mode characteristic impedances [Pa s/m^3].
    pub characteristic_impedances: Vec<C64>,
    /// Smallest plane-mode shear number along the chain.
    pub min_shear_number: f64,
    /// Mouth `ka`, using the exterior gas for a full flanged radiation load.
    pub mouth_ka: f64,
    /// Cutoff frequencies of modes 1..N-1 at the narrowest station [Hz].
    pub cutoffs_at_throat_hz: Vec<f64>,
    /// Number of staircase stations used in the recurrence.
    pub steps_used: usize,
}

impl ModalReflectionResponse {
    /// Reflection against positive real per-mode pressure-wave references.
    ///
    /// For `D = diag(references)`, form `A = (P + D U)/2` and
    /// `B = (P - D U)/2` from the same native waves and solve `R_new A = B`.
    /// This right solve preserves modal coupling and never inverts `I - R`.
    /// The returned waves carry signed power
    /// `sum((|a_i|^2 - |b_i|^2) / references[i]) / 2`; unequal references
    /// require that metric when testing passivity or reciprocity.
    ///
    /// # Errors
    /// Invalid shapes/references or a singular wave conversion. Modal cutoff,
    /// load and staircase validity limits are unchanged.
    pub fn at_references(&self, references: &[f64]) -> Result<Vec<C64>, DuctError> {
        let (p, q) = self.pressure_and_scaled_flow(references)?;
        let n = self.n_modes;
        // The factor 1/2 cancels between B and A. Solving A^T x = B_row^T
        // avoids an explicit inverse and the incorrect noncommuting A^-1 B.
        let mut a: Vec<_> = p.iter().zip(&q).map(|(p, q)| *p + *q).collect();
        let mut reflected: Vec<_> = p.iter().zip(&q).map(|(p, q)| *p - *q).collect();
        finite_waves(&a)?;
        finite_waves(&reflected)?;
        // Re-scale the native incident basis in both maps. This leaves B A^-1
        // unchanged and keeps quasi-static evanescent ports from setting the
        // entire factorization's scale through their very small native Zc.
        for j in 0..n {
            let scale = (0..n).map(|i| a[i * n + j].abs()).fold(0.0, f64::max);
            if scale == 0.0 {
                return Err(DuctError::Singular);
            }
            for i in 0..n {
                a[i * n + j] = a[i * n + j] / C64::from_re(scale);
                reflected[i * n + j] = reflected[i * n + j] / C64::from_re(scale);
            }
        }
        let lu = lu_complex(&transpose(&a, n), n).map_err(|_| DuctError::Singular)?;
        for row in reflected.chunks_exact_mut(n) {
            lu.solve(row);
        }
        finite_waves(&reflected)?;
        Ok(reflected)
    }

    /// Plane-source reflection with zero higher-mode inlet volume flow.
    ///
    /// This is the closure associated with the impedance matrix's `Z00`:
    /// `a0 = (P0 + reference U0)/2`, `b0 = (P0 - reference U0)/2`, and
    /// `U_i = 0` for i > 0. It differs from the (0,0) entry of a full
    /// scattering matrix, which instead matches the other incident ports.
    /// Solving those boundary equations directly retains finite plane
    /// reflection at a pressure antinode without forming `Z00` or the full Z.
    ///
    /// # Errors
    /// Invalid response/reference or a singular source-boundary solve.
    pub fn plane_reflection(&self, reference: f64) -> Result<C64, DuctError> {
        let n = self.n_modes;
        let references = [reference; MAX_MODES];
        let references = references.get(..n).ok_or(DuctError::BadParameter {
            what: "modal reflection mode count exceeds MAX_MODES",
        })?;
        let (p, q) = self.pressure_and_scaled_flow(references)?;
        let mut boundary: Vec<_> = self.reflection_matrix.iter().map(|r| -*r).collect();
        for i in 0..n {
            boundary[i * n + i] = boundary[i * n + i] + C64::ONE;
        }
        for j in 0..n {
            boundary[j] = p[j] + q[j];
        }
        let lu = lu_complex(&boundary, n).map_err(|_| DuctError::Singular)?;
        let mut outgoing = vec![C64::ZERO; n];
        outgoing[0] = C64::from_re(2.0);
        lu.solve(&mut outgoing);
        let reflected = (0..n)
            .fold(C64::ZERO, |sum, j| sum + (p[j] - q[j]) * outgoing[j])
            .scale(0.5);
        finite_waves(&[reflected])?;
        Ok(reflected)
    }

    fn pressure_and_scaled_flow(&self, references: &[f64]) -> Result<(Vec<C64>, Vec<C64>), DuctError> {
        let n = self.n_modes;
        if !(1..=MAX_MODES).contains(&n)
            || self.reflection_matrix.len() != n * n
            || self.characteristic_impedances.len() != n
            || references.len() != n
            || references.iter().any(|z| !z.is_finite() || *z <= 0.0)
            || self.characteristic_impedances.iter().any(|z| {
                !z.re.is_finite() || !z.im.is_finite() || z.abs() == 0.0
            })
        {
            return Err(DuctError::BadParameter {
                what: "modal reflection needs an N x N matrix, nonzero native impedances and positive real references",
            });
        }
        finite_waves(&self.reflection_matrix)?;
        let mut p = self.reflection_matrix.clone();
        let mut q = vec![C64::ZERO; n * n];
        for i in 0..n {
            let ratio = C64::from_re(references[i]) / self.characteristic_impedances[i];
            for j in 0..n {
                let delta = if i == j { C64::ONE } else { C64::ZERO };
                p[i * n + j] = delta + p[i * n + j];
                q[i * n + j] = ratio * (delta - self.reflection_matrix[i * n + j]);
            }
        }
        finite_waves(&p)?;
        finite_waves(&q)?;
        Ok((p, q))
    }
}

fn finite_waves(waves: &[C64]) -> Result<(), DuctError> {
    if waves.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
        return Err(DuctError::BadParameter { what: "modal reflection waves must be finite" });
    }
    Ok(())
}

#[cfg(test)]
mod tests;
