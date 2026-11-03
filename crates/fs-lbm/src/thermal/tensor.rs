//! Full symmetric thermal diffusion in the existing D2Q5/D3Q7 slabs.
//!
//! The flux block is S = (I/2 + alpha/c_T^2)^(-1), with c_T^2 = 1/3
//! in two dimensions and 1/4 in three. Off-diagonal entries retain the
//! supplied physical axes; they are not replaced by directional averages.
//! See Yoshida & Nagaoka, JCP 229 (2010), 7774-7795,
//! <https://doi.org/10.1016/j.jcp.2010.06.037>.
//!
//! These slabs evolve temperature with uniform volumetric heat capacity.
//! For a physical conductivity tensor K and a common rho*Cp, supply
//! alpha_lattice = K/(rho*Cp) * dt/dx^2 on a square/cubic grid. This does
//! not infer material coefficients, variable storage or interface resistance.

use super::ThermalLbm;
use super::sources::{
    ThermalSourceBalance, ThermalSourceError, ThermalSourceLimits, source_step_tensor_mapped,
};
use crate::d3q19::ThermalLbm3;

/// A prepared SPD thermal diffusivity tensor in lattice-cell^2/lattice-step.
///
/// Construction is restricted to D=2 or D=3 and requires exact symmetry.
/// The represented relaxation and 2I-S must both be positive definite; tiny
/// diffusivities lost in the half-step offset and unresolved inverses refuse.
/// No floor, cap, eigenvalue clipping or coefficient symmetrization is applied.
/// Preparation is independent of the cells and may be reused across steps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffusivityTensor<const D: usize> {
    diffusivity: [[f64; D]; D],
    pub(crate) relaxation: [[f64; D]; D],
}

impl<const D: usize> DiffusivityTensor<D> {
    /// Admit a finite, symmetric, strictly positive-definite tensor expressed
    /// in the fixed Cartesian axes of the square/cubic lattice.
    pub fn from_lattice(diffusivity: [[f64; D]; D]) -> Result<Self, ThermalSourceError> {
        let cs2 = match D {
            2 => 1.0 / 3.0,
            3 => 0.25,
            _ => return Err(invalid("thermal tensors require two or three dimensions")),
        };
        for row in 0..D {
            for column in 0..D {
                if !diffusivity[row][column].is_finite()
                    || diffusivity[row][column] != diffusivity[column][row]
                {
                    return Err(invalid(
                        "thermal diffusivity must be finite and exactly symmetric",
                    ));
                }
            }
        }
        let scale = (0..D)
            .map(|axis| diffusivity[axis][axis])
            .fold(0.0_f64, f64::max);
        if scale <= 0.0 {
            return Err(invalid("thermal diffusivity must be positive definite"));
        }
        let normalized: [[f64; D]; D] = core::array::from_fn(|row| {
            core::array::from_fn(|column| diffusivity[row][column] / scale)
        });
        fs_la::factor::cholesky(normalized.as_flattened(), D)
            .map_err(|_| invalid("thermal diffusivity must be positive definite"))?;
        let tau: [[f64; D]; D] = core::array::from_fn(|row| {
            core::array::from_fn(|column| {
                diffusivity[row][column] / cs2 + if row == column { 0.5 } else { 0.0 }
            })
        });
        if !tau.as_flattened().iter().all(|entry| entry.is_finite()) {
            return Err(invalid("thermal tensor relaxation is not representable"));
        }
        let factor = fs_la::factor::cholesky(tau.as_flattened(), D)
            .map_err(|_| invalid("thermal tensor relaxation is not positive definite"))?;
        let mut relaxation = [[0.0; D]; D];
        for column in 0..D {
            let mut basis = [0.0; D];
            basis[column] = 1.0;
            factor.solve(&mut basis);
            // Store the symmetric inverse from one triangle of the solve.
            for row in column..D {
                relaxation[row][column] = basis[row];
                relaxation[column][row] = basis[row];
            }
        }
        let complement: [[f64; D]; D] = core::array::from_fn(|row| {
            core::array::from_fn(|column| {
                (if row == column { 2.0 } else { 0.0 }) - relaxation[row][column]
            })
        });
        fs_la::factor::cholesky(relaxation.as_flattened(), D)
            .and_then(|_| fs_la::factor::cholesky(complement.as_flattened(), D))
            .map_err(|_| {
                invalid("thermal tensor flux relaxation must have eigenvalues in (0, 2)")
            })?;
        for row in 0..D {
            let mut residual = 0.0;
            for column in 0..D {
                let mut entry = if row == column { -1.0 } else { 0.0 };
                for axis in 0..D {
                    entry = tau[row][axis].mul_add(relaxation[axis][column], entry);
                }
                residual += entry.abs();
            }
            if !residual.is_finite() || residual > 256.0 * f64::EPSILON {
                return Err(invalid("thermal tensor inverse is numerically unresolved"));
            }
        }
        Ok(Self {
            diffusivity,
            relaxation,
        })
    }

    /// Original Cartesian tensor, in lattice-cell^2/lattice-step.
    #[must_use]
    pub fn lattice_diffusivity(&self) -> [[f64; D]; D] {
        self.diffusivity
    }
}

/// A borrowed diffusion field, fixed over one thermal interval.
///
/// Uniform fields reuse one factorization for the entire lattice. Per-cell
/// fields include both wall planes; wall tensors are not used in collision.
/// The background tau_g still relaxes even kinetic moments and remains the
/// scalar diagnostic value returned by diffusivity/rayleigh/nusselt.
#[derive(Clone, Copy, Debug)]
pub enum DiffusivityField<'a, const D: usize> {
    /// One prepared law shared by all fluid cells.
    Uniform(&'a DiffusivityTensor<D>),
    /// One prepared law for each cell, including unused wall entries.
    PerCell(&'a [DiffusivityTensor<D>]),
}

impl<'a, const D: usize> DiffusivityField<'a, D> {
    pub(crate) fn admit(self, cells: usize) -> Result<(), ThermalSourceError> {
        if let Self::PerCell(field) = self {
            if field.len() != cells {
                return Err(invalid(
                    "one thermal tensor per cell, including walls, is required",
                ));
            }
        }
        Ok(())
    }

    pub(crate) fn at(self, cell: usize) -> &'a DiffusivityTensor<D> {
        match self {
            Self::Uniform(tensor) => tensor,
            Self::PerCell(field) => &field[cell],
        }
    }
}

impl ThermalLbm {
    /// Advance the original D2Q9/D2Q5 slab with full tensor heat diffusion.
    ///
    /// Uses the existing source kicks, physical velocity/force frame, retained
    /// advective-flux correction, wall heat measurement and atomic acceptance.
    /// Momentum and even thermal moments retain their original relaxation.
    /// Tensor coefficients are borrowed for this step; no material or tensor
    /// remains installed in a later scalar step. Spatial/temporal refinement
    /// remains required; SPD does not establish population positivity or the
    /// stability of the coupled buoyancy problem.
    pub fn step_with_sources_and_tensor_diffusion(
        &mut self,
        heat_per_step: &[f64],
        forces: &[[f64; 2]],
        diffusion: DiffusivityField<'_, 2>,
        limits: ThermalSourceLimits,
    ) -> Result<ThermalSourceBalance, ThermalSourceError> {
        source_step_tensor_mapped::<Self, 2, 5, _>(
            self,
            heat_per_step,
            forces,
            limits,
            diffusion,
            |_, balance| Ok(balance),
        )
    }
}

impl ThermalLbm3 {
    /// Three-dimensional counterpart on the original D3Q19/D3Q7 slab.
    /// The full Cartesian flux block, source balance and rollback are shared.
    pub fn step_with_sources_and_tensor_diffusion(
        &mut self,
        heat_per_step: &[f64],
        forces: &[[f64; 3]],
        diffusion: DiffusivityField<'_, 3>,
        limits: ThermalSourceLimits,
    ) -> Result<ThermalSourceBalance, ThermalSourceError> {
        source_step_tensor_mapped::<Self, 3, 7, _>(
            self,
            heat_per_step,
            forces,
            limits,
            diffusion,
            |_, balance| Ok(balance),
        )
    }
}

fn invalid(reason: &'static str) -> ThermalSourceError {
    ThermalSourceError::Invalid(reason)
}
