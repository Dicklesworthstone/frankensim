//! Signed tangent modes with physical shapes and uniform-base participation.

use super::{Frame2d, SolidError};
use fs_la::eigen::{admit_jacobi_eigh, jacobi_eigh};
use fs_la::factor::cholesky;

/// Numerical interpretation of a signed tangent eigenvalue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameModeStatus {
    /// Positive eigenvalue resolved above the reported numerical scale.
    Oscillatory,
    /// Negative eigenvalue resolved below the reported numerical scale.
    Unstable,
    /// Rigid motion or an eigenvalue too close to zero to resolve its sign.
    /// This does not certify neutral stability.
    Unresolved,
}

/// One mode of the current tangent and the original reference mass.
#[derive(Debug, Clone, PartialEq)]
pub struct FrameMode {
    /// Signed lambda in K*phi = lambda*M*phi [s^-2], never clipped.
    pub eigenvalue_per_s2: f64,
    /// Full node-major [ux, uy, theta] shape, phi^T M phi = 1.
    /// Fixed entries are zero; massless entries are statically restored.
    pub displacement: Vec<f64>,
    /// Gamma = phi^T M r for unit horizontal/vertical base translation.
    pub base_participation: [f64; 2],
    /// Gamma^2 for the mass-normalized mode [kg].
    pub effective_mass_kg: [f64; 2],
    /// Free-pencil backward residual ||Kphi-lambda*Mphi|| /
    /// ((||K||_F + |lambda|*||M||_F)*||phi||).
    pub relative_residual: f64,
    /// Recomputed residual in mass-whitened coordinates [s^-2].
    pub spectral_residual_per_s2: f64,
    /// Residual plus 64*eps*n*||L^-1 K L^-T||_inf [s^-2].
    /// This is a floating-point resolution diagnostic, not a certified bound.
    pub eigenvalue_resolution_per_s2: f64,
}

impl FrameMode {
    /// Classify without modifying the signed eigenvalue.
    #[must_use]
    pub fn status(&self) -> FrameModeStatus {
        if self.eigenvalue_per_s2 > self.eigenvalue_resolution_per_s2 {
            FrameModeStatus::Oscillatory
        } else if self.eigenvalue_per_s2 < -self.eigenvalue_resolution_per_s2 {
            FrameModeStatus::Unstable
        } else {
            FrameModeStatus::Unresolved
        }
    }

    /// Circular frequency of a resolved oscillatory mode [rad/s].
    #[must_use]
    pub fn angular_frequency(&self) -> Option<f64> {
        (self.status() == FrameModeStatus::Oscillatory).then(|| self.eigenvalue_per_s2.sqrt())
    }

    /// Exponential growth rate of a resolved unstable linearized mode [s^-1].
    #[must_use]
    pub fn growth_rate(&self) -> Option<f64> {
        (self.status() == FrameModeStatus::Unstable).then(|| (-self.eigenvalue_per_s2).sqrt())
    }
}

/// Ascending signed tangent modes and their full reference influence inventory.
#[derive(Debug, Clone, PartialEq)]
pub struct FrameModalAnalysis {
    /// Mass-normalized modes ordered by signed eigenvalue.
    pub modes: Vec<FrameMode>,
    /// r^T M r for horizontal/vertical rigid translation [kg], including
    /// constrained mass. Modal effective masses need not exhaust this inventory
    /// when supported inertia also contributes directly to the base reaction.
    pub base_influence_mass_kg: [f64; 2],
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}
fn norm(x: &[f64]) -> f64 {
    x.iter().fold(0.0_f64, |n, x| n.hypot(*x))
}

impl Frame2d {
    /// Mass-normalized modes of the current tangent, including signed unstable
    /// eigenvalues and full base participation. Material state is not changed.
    ///
    /// Massless free coordinates must have positive-definite stiffness before
    /// static condensation; an unstable or unconstrained massless subsystem is
    /// refused. Physical shapes restore those coordinates using the same
    /// condensation map. Repeated eigenvalues define a basis, not unique modes.
    ///
    /// # Errors
    /// Invalid inertia, inadmissible dense work, unresolved eigenpairs or an
    /// unstable/underconstrained massless block. Instability of massive modes is
    /// returned explicitly in the signed spectrum, not converted to zero.
    #[allow(clippy::too_many_lines)]
    pub fn modal_analysis(&self) -> Result<FrameModalAnalysis, SolidError> {
        let n = self.ndof();
        let mass = self.reference_mass()?;
        let (t, r) = self.mass_partition(&mass);
        let (nt, nr) = (t.len(), r.len());
        admit_jacobi_eigh(nt)
            .map_err(|_| invalid("frame modal dimension exceeds the eigensolver work allowance"))?;
        let free: Vec<usize> = (0..n).filter(|&i| !self.fixed[i]).collect();
        let influence: [Vec<f64>; 2] = [0, 1].map(|axis| {
            (0..n)
                .map(|i| if i % 3 == axis { 1.0 } else { 0.0 })
                .collect()
        });
        let mass_influence = influence
            .each_ref()
            .map(|x| (0..n).map(|i| mass.row_action(i, x)).collect::<Vec<_>>());
        let inventory = [0, 1].map(|axis| {
            influence[axis]
                .iter()
                .zip(&mass_influence[axis])
                .map(|(r, mr)| r * mr)
                .sum::<f64>()
        });
        if inventory.iter().any(|m| !m.is_finite()) {
            return Err(invalid("nonfinite frame reference influence mass"));
        }
        if nt == 0 {
            if nr > 0 {
                return Err(invalid("free frame coordinates have no reference inertia"));
            }
            return Ok(FrameModalAnalysis {
                modes: Vec::new(),
                base_influence_mass_kg: inventory,
            });
        }
        let (_, k, _) = self.assemble(&self.u)?;
        if k.iter().any(|v| !v.is_finite()) {
            return Err(invalid("nonfinite frame tangent"));
        }
        let mut kcc = vec![0.0; nt * nt];
        let mut transfer = vec![0.0; nr * nt];
        for i in 0..nt {
            for j in 0..nt {
                kcc[i * nt + j] = k[t[i] * n + t[j]];
            }
        }
        if nr > 0 {
            let mut krr = vec![0.0; nr * nr];
            for i in 0..nr {
                for j in 0..nr {
                    krr[i * nr + j] = k[r[i] * n + r[j]];
                }
            }
            let factor = cholesky(&krr, nr).map_err(|_| invalid(
                "massless frame stiffness is not positive definite; unstable or unconstrained coordinates cannot be condensed"))?;
            for j in 0..nt {
                let mut column: Vec<f64> = r.iter().map(|&i| k[i * n + t[j]]).collect();
                factor.solve(&mut column);
                for i in 0..nr {
                    transfer[i * nt + j] = -column[i];
                }
                for i in 0..nt {
                    for p in 0..nr {
                        kcc[i * nt + j] -= k[t[i] * n + r[p]] * column[p];
                    }
                }
            }
        }
        let (a, transform) = mass.modal_matrix(&kcc, &t)?;
        let a: Vec<f64> = a.into_iter().flatten().collect();
        let scale = a
            .chunks(nt)
            .map(|row| row.iter().map(|x| x.abs()).sum::<f64>())
            .fold(0.0_f64, f64::max);
        if !scale.is_finite() {
            return Err(invalid("nonfinite frame modal scale"));
        }
        let (values, vectors) = jacobi_eigh(&a, nt);
        let mut k_norm = 0.0_f64;
        let mut m_norm = 0.0_f64;
        for &i in &free {
            for &j in &free {
                k_norm = k_norm.hypot(k[i * n + j]);
                m_norm = m_norm.hypot(mass.entry(i, j));
            }
        }
        let mut modes = Vec::with_capacity(nt);
        for (mode, &value) in values.iter().enumerate() {
            let y: Vec<f64> = (0..nt).map(|i| vectors[i * nt + mode]).collect();
            let reduced_residual: Vec<f64> = (0..nt)
                .map(|i| (0..nt).map(|j| a[i * nt + j] * y[j]).sum::<f64>() - value * y[i])
                .collect();
            let spectral_residual = norm(&reduced_residual) / norm(&y);
            let active = transform.restore(&y);
            let mut shape = vec![0.0; n];
            for (i, &d) in t.iter().enumerate() {
                shape[d] = active[i];
            }
            for (i, &d) in r.iter().enumerate() {
                shape[d] = (0..nt).map(|j| transfer[i * nt + j] * active[j]).sum();
            }
            let modal_mass: f64 = (0..n).map(|i| shape[i] * mass.row_action(i, &shape)).sum();
            if !modal_mass.is_finite() || modal_mass <= 0.0 {
                return Err(invalid("frame eigenvector has invalid modal mass"));
            }
            let mut pivot = 0;
            for i in 1..n {
                if shape[i].abs() > shape[pivot].abs() {
                    pivot = i;
                }
            }
            let normalization = modal_mass.sqrt().copysign(shape[pivot]);
            for x in &mut shape {
                *x /= normalization;
            }
            let residual: Vec<f64> = free
                .iter()
                .map(|&i| {
                    (0..n).map(|j| k[i * n + j] * shape[j]).sum::<f64>()
                        - value * mass.row_action(i, &shape)
                })
                .collect();
            let denominator = (k_norm + value.abs() * m_norm) * norm(&shape);
            let relative_residual = if denominator > 0.0 {
                norm(&residual) / denominator
            } else {
                norm(&residual)
            };
            let resolution = spectral_residual + 64.0 * f64::EPSILON * nt as f64 * scale;
            let participation = [0, 1].map(|axis| {
                shape
                    .iter()
                    .zip(&mass_influence[axis])
                    .map(|(phi, mr)| phi * mr)
                    .sum::<f64>()
            });
            let effective = participation.map(|gamma| gamma * gamma);
            if !value.is_finite()
                || !relative_residual.is_finite()
                || !denominator.is_finite()
                || !resolution.is_finite()
                || relative_residual > 1e-8
                || spectral_residual > 1e-8 * scale.max(f64::MIN_POSITIVE)
                || shape
                    .iter()
                    .chain(&participation)
                    .chain(&effective)
                    .any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "frame modal eigenpair did not resolve the full tangent/mass pencil",
                ));
            }
            modes.push(FrameMode {
                eigenvalue_per_s2: value,
                displacement: shape,
                base_participation: participation,
                effective_mass_kg: effective,
                relative_residual,
                spectral_residual_per_s2: spectral_residual,
                eigenvalue_resolution_per_s2: resolution,
            });
        }
        Ok(FrameModalAnalysis {
            modes,
            base_influence_mass_kg: inventory,
        })
    }

    /// Circular natural frequencies of a resolved oscillatory tangent.
    /// Rigid, unresolved and unstable modes are inspected with modal_analysis;
    /// this convenience method never clips a negative eigenvalue to zero.
    ///
    /// # Errors
    /// Modal extraction failures or any non-oscillatory mode.
    pub fn natural_frequencies(&self) -> Result<Vec<f64>, SolidError> {
        self.modal_analysis()?
            .modes
            .iter()
            .map(|mode| {
                mode.angular_frequency().ok_or_else(|| {
                    invalid(&format!(
                        "frame tangent has {:?} mode: signed eigenvalue {}, resolution {}",
                        mode.status(),
                        mode.eigenvalue_per_s2,
                        mode.eigenvalue_resolution_per_s2
                    ))
                })
            })
            .collect()
    }
}
