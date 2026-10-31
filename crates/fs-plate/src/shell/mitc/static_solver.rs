//! Reusable factorization of one immutable native MITC elastic operator.
//!
//! Prescribed loads and thermal eigenstrains may change between solves. The
//! original model, support map and stiffness must stay paired with their factor.

use super::{MitcShellModel, PlateError, dd_add, two_prod};
use fs_math::det;
use fs_sparse::{Coo, DirectOrdering, LdltFactor, LdltOptions, NdOptions, SymbolicLdlt};
use std::sync::Arc;

/// An immutable shell and its positive-definite static factorization.
///
/// Repeated solves reuse the same diagonally equilibrated sparse LDLT factor.
/// Each right-hand side still receives compensated iterative refinement and
/// the original unscaled true-residual admission check. The model and factor
/// cannot be replaced independently. Sharing the model does not copy its mesh
/// or sparse pencil, so thermal candidates can retain the same elastic owner.
#[derive(Debug)]
pub struct PreparedMitcStatic {
    model: Arc<MitcShellModel>,
    factor: StaticFactor,
}

impl PreparedMitcStatic {
    /// Prepare an owned or already shared shell for repeated world-load solves.
    ///
    /// Free-body, singular and indefinite operators retain the ordinary static
    /// solver's refusals. Cancellation is checked between sparse-kernel calls
    /// and at bounded row intervals; it cannot interrupt an individual sparse
    /// ordering, factorization or triangular-solve kernel.
    pub fn try_new(
        model: impl Into<Arc<MitcShellModel>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PlateError> {
        checkpoint(&mut cancelled)?;
        let model = model.into();
        let factor = StaticFactor::try_new(&model, &mut cancelled)?;
        Ok(Self { model, factor })
    }

    /// The exact immutable shell whose stiffness was factored.
    #[must_use]
    pub fn model(&self) -> &MitcShellModel {
        &self.model
    }

    /// Share the original shell without copying its sparse pencil or geometry.
    #[must_use]
    pub fn shared_model(&self) -> Arc<MitcShellModel> {
        Arc::clone(&self.model)
    }

    /// Solve one finite physical force/moment vector and return world motion.
    ///
    /// Thermal loads must be prepared from [`Self::model`] and added to any
    /// external load before calling. Recovery uses that candidate's thermal
    /// owner, including its affine internal bubble response. Failure or
    /// cancellation leaves the shared model and cached factor unchanged.
    pub fn solve(
        &self,
        full_world_load: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, PlateError> {
        checkpoint(&mut cancelled)?;
        let rhs = self.model.project_load(full_world_load)?;
        self.factor.solve_reduced(&self.model, &rhs, &mut cancelled)
    }
}

/// Borrowed compatibility path: no model or sparse-pencil clone is needed.
pub(super) fn solve_static(
    model: &MitcShellModel,
    full_load_6n: &[f64],
) -> Result<Vec<f64>, PlateError> {
    // Preserve the existing load-admission order before expensive factor work.
    let rhs = model.project_load(full_load_6n)?;
    let mut cancelled = || false;
    let factor = StaticFactor::try_new(model, &mut cancelled)?;
    factor.solve_reduced(model, &rhs, &mut cancelled)
}

#[derive(Debug)]
struct StaticFactor {
    scale: Vec<f64>,
    factor: LdltFactor,
}

impl StaticFactor {
    fn try_new(
        model: &MitcShellModel,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PlateError> {
        checkpoint(cancelled)?;
        let n = model.pencil.free;
        let mut scale = vec![0.0; n];
        for (row, s) in scale.iter_mut().enumerate() {
            row_checkpoint(row, cancelled)?;
            let (cols, vals) = model.pencil.k.row(row);
            let mut diag = 0.0;
            for (&col, &val) in cols.iter().zip(vals) {
                if col == row {
                    diag = val;
                    break;
                }
            }
            if !diag.is_finite() || diag <= 0.0 {
                return Err(PlateError::StaticSolve {
                    what: "static shell stiffness is not positive definite",
                });
            }
            *s = 1.0 / det::sqrt(diag);
        }
        let mut k_eq = Coo::new(n, n);
        for row in 0..n {
            row_checkpoint(row, cancelled)?;
            let sr = scale[row];
            let (cols, vals) = model.pencil.k.row(row);
            for (&col, &val) in cols.iter().zip(vals) {
                k_eq.push(row, col, sr * val * scale[col]);
            }
        }
        checkpoint(cancelled)?;
        let k_eq = k_eq.assemble();
        let ordering = if n >= 4_096 {
            DirectOrdering::NestedDissection(NdOptions {
                leaf_size: 256,
                refinement_passes: 2,
            })
        } else {
            DirectOrdering::Amd
        };
        checkpoint(cancelled)?;
        let symbolic = SymbolicLdlt::analyze(&k_eq, ordering).map_err(PlateError::StaticFactor)?;
        checkpoint(cancelled)?;
        let factor = symbolic
            .factor(&k_eq, &LdltOptions::default())
            .map_err(PlateError::StaticFactor)?;
        checkpoint(cancelled)?;
        if factor.inertia().negative != 0 || factor.inertia().positive != n {
            return Err(PlateError::StaticSolve {
                what: "static shell stiffness is not positive definite",
            });
        }
        Ok(Self { scale, factor })
    }

    fn solve_reduced(
        &self,
        model: &MitcShellModel,
        rhs: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, PlateError> {
        checkpoint(cancelled)?;
        let n = model.pencil.free;
        let scale = &self.scale;
        let mut rhs_eq = vec![0.0; n];
        for row in 0..n {
            row_checkpoint(row, cancelled)?;
            rhs_eq[row] = scale[row] * rhs[row];
        }
        checkpoint(cancelled)?;
        let y0 = self.factor.solve(&rhs_eq);
        checkpoint(cancelled)?;
        let mut x_hi = vec![0.0; n];
        let mut x_lo = vec![0.0; n];
        for i in 0..n {
            row_checkpoint(i, cancelled)?;
            let (hi, lo) = two_prod(scale[i], y0[i]);
            x_hi[i] = hi;
            x_lo[i] = lo;
        }
        let mut r = vec![0.0; n];
        let mut r_eq = vec![0.0; n];
        let load = rhs.iter().fold(0.0_f64, |acc, v| acc.hypot(*v));
        let mut residual = f64::INFINITY;
        for step in 0..5 {
            checkpoint(cancelled)?;
            for (row, &b_i) in rhs.iter().enumerate() {
                row_checkpoint(row, cancelled)?;
                let (cols, vals) = model.pencil.k.row(row);
                let mut sum = (b_i, 0.0_f64);
                for (&col, &val) in cols.iter().zip(vals) {
                    let (p_hi, p_lo) = two_prod(-val, x_hi[col]);
                    sum = dd_add(sum, (p_hi, p_lo - val * x_lo[col]));
                }
                r[row] = sum.0 + sum.1;
                r_eq[row] = scale[row] * r[row];
            }
            residual = r.iter().fold(0.0_f64, |acc, v| acc.hypot(*v));
            if !residual.is_finite() || residual <= 1e-11 * load.max(f64::MIN_POSITIVE) || step == 4
            {
                break;
            }
            checkpoint(cancelled)?;
            let delta_eq = self.factor.solve(&r_eq);
            checkpoint(cancelled)?;
            for i in 0..n {
                row_checkpoint(i, cancelled)?;
                let (d_hi, d_lo) = two_prod(scale[i], delta_eq[i]);
                let (hi, lo) = dd_add((x_hi[i], x_lo[i]), (d_hi, d_lo));
                x_hi[i] = hi;
                x_lo[i] = lo;
            }
        }
        checkpoint(cancelled)?;
        if !residual.is_finite() || residual > 1e-8 * load.max(f64::MIN_POSITIVE) {
            return Err(PlateError::StaticSolve {
                what: "static shell true residual exceeds 1e-8 relative",
            });
        }
        let displacement = model.expand(&x_hi)?;
        checkpoint(cancelled)?;
        Ok(displacement)
    }
}

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), PlateError> {
    if cancelled() {
        Err(PlateError::Cancelled)
    } else {
        Ok(())
    }
}

fn row_checkpoint(row: usize, cancelled: &mut impl FnMut() -> bool) -> Result<(), PlateError> {
    if row % 256 == 0 {
        checkpoint(cancelled)?;
    }
    Ok(())
}
