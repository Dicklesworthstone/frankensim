//! Original signed radial mouth coordinates and one physical matrix sample.

use super::{ModeCount, bessel_j0_certified, j1_roots_unchecked};
use crate::DuctError;
use fs_material::gas::GasState;
use fs_math::c64::C64;

fn bad(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}

/// Area-orthonormal axisymmetric rigid-wall modes on a circular mouth.
///
/// `psi_n(eta) = J0(gamma_n eta)/|J0(gamma_n)|`, with `gamma_0 = 0` and
/// the remaining eigenvalues the positive zeros of J1. Thus
/// `(1/S) integral psi_m psi_n dS = delta_mn`. The pressure coordinate is
/// `P_n = (1/S) integral p psi_n dS`; the outward volume-flow coordinate is
/// `U_n = integral v psi_n dS`. Only U0 is the net volume flow. Higher modes
/// are signed, zero-net-flow traces, never independent uniform pistons.
#[derive(Debug, Clone)]
pub struct AxisymmetricMouthBasis {
    eigenvalues: Vec<f64>,
    normalization: Vec<f64>,
}

impl AxisymmetricMouthBasis {
    /// Prepare the same certified radial basis as continuous duct propagation.
    ///
    /// # Errors
    /// A mode count outside the existing certified Bessel range.
    pub fn new(modes: ModeCount) -> Result<Self, DuctError> {
        modes.validate()?;
        let mut eigenvalues = vec![0.0];
        eigenvalues.extend(j1_roots_unchecked(modes.0 - 1));
        let normalization = eigenvalues
            .iter()
            .map(|&g| bessel_j0_certified(g).abs())
            .collect();
        Ok(Self {
            eigenvalues,
            normalization,
        })
    }

    /// Retained radial modes, including the plane mode.
    #[must_use]
    pub fn modes(&self) -> ModeCount {
        ModeCount(self.eigenvalues.len())
    }

    /// Radial eigenvalues, including the plane mode's zero eigenvalue.
    #[must_use]
    pub fn radial_eigenvalues(&self) -> &[f64] {
        &self.eigenvalues
    }

    /// Original signed mode at `radius_fraction = r/R` in the closed unit disk.
    /// A unit modal-volume-flow boundary trace is this value divided by S.
    ///
    /// # Errors
    /// An absent mode or nonfinite/out-of-disk radial coordinate.
    pub fn shape(&self, mode: usize, radius_fraction: f64) -> Result<f64, DuctError> {
        if mode >= self.eigenvalues.len()
            || !radius_fraction.is_finite()
            || !(0.0..=1.0).contains(&radius_fraction)
        {
            return Err(bad(
                "modal mouth shape needs an admitted mode and radius fraction in [0,1]",
            ));
        }
        Ok(
            bessel_j0_certified(self.eigenvalues[mode] * radius_fraction)
                / self.normalization[mode],
        )
    }
}

/// Full mouth relation `P = Z U_outward` [Pa s/m³], at one positive frequency.
///
/// The row-major matrix uses [`AxisymmetricMouthBasis`] and `exp(-i omega t)`.
/// Its Hermitian part must be positive semidefinite, so average outward power
/// is `Re(U^H Z U)/2 >= 0` for every signed modal flow, including evanescent
/// modes. Admission tests the complete matrix with floating-point tolerance;
/// no clipping, symmetrization, diagonalization or per-mode net-flow division
/// changes the supplied entries. Reciprocity is not imposed on a general load.
///
/// This is a frequency-bound physical sample, not a fit, interpolator or
/// causality certificate. The caller owns its exterior geometry and source's
/// validity band; a new sample is required at every evaluated frequency.
#[derive(Debug, Clone)]
pub struct ModalMouthImpedance {
    omega_rad_s: f64,
    radius_m: f64,
    exterior: GasState,
    basis: AxisymmetricMouthBasis,
    impedance: Vec<C64>,
}

impl ModalMouthImpedance {
    /// Admit and retain the original matrix and its physical coordinate system.
    ///
    /// # Errors
    /// Invalid gas, nonpositive/nonfinite frequency or radius, wrong shape,
    /// nonfinite entries, or a negative eigenvalue of the Hermitian part beyond
    /// `512 eps N` times its largest component. The native eigensolve's residual
    /// is checked at that same normalized tolerance. The tolerance is relative
    /// to resistance, so large lossless reactance cannot hide active loss.
    pub fn new(
        omega_rad_s: f64,
        radius_m: f64,
        exterior: GasState,
        modes: ModeCount,
        impedance: Vec<C64>,
    ) -> Result<Self, DuctError> {
        let basis = AxisymmetricMouthBasis::new(modes)?;
        crate::validate_gas(&exterior)?;
        if !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || !radius_m.is_finite()
            || radius_m <= 0.0
        {
            return Err(bad(
                "modal mouth needs finite positive frequency and radius",
            ));
        }
        if impedance.len() != modes.0 * modes.0
            || impedance
                .iter()
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(bad(
                "modal mouth needs one finite row-major N by N impedance matrix",
            ));
        }
        admit_passive(&impedance, modes.0)?;
        Ok(Self {
            omega_rad_s,
            radius_m,
            exterior,
            basis,
            impedance,
        })
    }

    /// Original angular frequency [rad/s], matched exactly by every evaluation.
    #[must_use]
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }
    /// Circular mouth radius [m] where this matrix is physically defined.
    #[must_use]
    pub const fn radius_m(&self) -> f64 {
        self.radius_m
    }
    /// Original coherent exterior gas, distinct from the interior gas.
    #[must_use]
    pub const fn exterior(&self) -> &GasState {
        &self.exterior
    }
    /// Full retained radial mode count.
    #[must_use]
    pub fn modes(&self) -> ModeCount {
        self.basis.modes()
    }
    /// Original row-major acoustic impedance [Pa s/m³], unchanged by admission.
    #[must_use]
    pub fn impedance(&self) -> &[C64] {
        &self.impedance
    }
    /// Original signed, area-orthonormal pressure/volume-flow coordinates.
    #[must_use]
    pub const fn basis(&self) -> &AxisymmetricMouthBasis {
        &self.basis
    }

    pub(super) fn admit(
        &self,
        interior: &GasState,
        omega: f64,
        radius: f64,
        modes: usize,
    ) -> Result<(), DuctError> {
        if modes != self.modes().0 {
            return Err(bad("modal mouth matrix and retained mode count must agree"));
        }
        self.admit_prefix(interior, omega, radius, modes)
    }

    // Only the continuous N/2 comparison uses a principal submatrix. Public
    // entry points require the complete declared rung before propagation.
    pub(super) fn admit_prefix(
        &self,
        interior: &GasState,
        omega: f64,
        radius: f64,
        modes: usize,
    ) -> Result<(), DuctError> {
        crate::exterior::admit(interior, &self.exterior)?;
        if omega.to_bits() != self.omega_rad_s.to_bits() {
            return Err(bad("modal mouth matrix belongs to a different frequency"));
        }
        if !radius.is_finite()
            || radius <= 0.0
            || (radius - self.radius_m).abs() > 128.0 * f64::EPSILON * radius.max(self.radius_m)
        {
            return Err(bad(
                "modal mouth matrix belongs to a different mouth radius",
            ));
        }
        if modes == 0 || modes > self.modes().0 {
            return Err(bad(
                "modal mouth matrix does not contain the requested radial basis",
            ));
        }
        Ok(())
    }
}

fn admit_passive(z: &[C64], n: usize) -> Result<(), DuctError> {
    // Real symmetric representation of H=(Z+Z^H)/2. Its eigenvalues are those
    // of H, each twice. Midpoint avoids overflow and preserves subnormal loss.
    let size = 2 * n;
    let mut h = vec![0.0; size * size];
    let mut scale: f64 = 0.0;
    for i in 0..n {
        for j in 0..n {
            let real = z[i * n + j].re.midpoint(z[j * n + i].re);
            let imag = z[i * n + j].im.midpoint(-z[j * n + i].im);
            h[i * size + j] = real;
            h[(i + n) * size + j + n] = real;
            h[i * size + j + n] = -imag;
            h[(i + n) * size + j] = imag;
            scale = scale.max(real.abs()).max(imag.abs());
        }
    }
    if scale == 0.0 {
        return Ok(());
    }
    for value in &mut h {
        *value /= scale;
    }
    let (values, vectors) = fs_la::eigen::jacobi_eigh(&h, size);
    if vectors.iter().any(|value| !value.is_finite()) {
        return Err(bad(
            "modal mouth Hermitian passivity eigensolve returned nonfinite vectors",
        ));
    }
    let tolerance = 256.0 * f64::EPSILON * size as f64;
    for (column, &value) in values.iter().enumerate() {
        if !value.is_finite() || value < -tolerance {
            return Err(bad(
                "modal mouth Hermitian resistance must be positive semidefinite",
            ));
        }
        let mut residual: f64 = 0.0;
        for i in 0..size {
            let applied: f64 = (0..size)
                .map(|j| h[i * size + j] * vectors[j * size + column])
                .sum();
            residual = residual.max((applied - value * vectors[i * size + column]).abs());
        }
        if !residual.is_finite() || residual > tolerance {
            return Err(bad(
                "modal mouth Hermitian passivity eigensolve did not resolve its residual",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
