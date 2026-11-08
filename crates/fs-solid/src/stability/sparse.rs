//! Buckling with bounded sparse factorization.
//!
//! K is admitted by the complete native sparse LDLᵀ factor, with no
//! regularization and 2×2 pivots disabled. Positive scalar pivots supply the
//! whitening S, Sᵀ K S = I. The generic fs-spectral service applies the
//! symmetric operator Sᵀ(−G)S and selects its largest algebraic roots, so
//! tensile directions do not displace the first compressive buckling modes.
//!
//! Sparse fill is explicitly bounded; no global dense K, G or inverse is
//! assembled. This is a sparse-factor route, not a factor-free algorithm.
//! Residual-qualified Ritz estimates do not certify spectrum completeness,
//! mesh convergence, nonlinear collapse or stability of an indefinite base.

use super::{BucklingResult, check_cancelled, scaled_symmetric};
use crate::SolidError;
use fs_sparse::direct::{
    DirectLimits,
    positive::{PositiveFactorError, PositiveWhitening},
};
use fs_sparse::{Csr, DirectOrdering, LdltError, LdltOptions, SparseApplyError, SymbolicLdlt};
use fs_spectral::service::{
    CertifiedEigenvalue, EigenBackend, EigenQuery, EigenService, ServiceError, SymmetricOp,
};
use std::cell::{Cell, RefCell};

/// Explicit input, sparse-factor and operator allowances. Outer eigeniteration
/// is separately limited by the `steps` argument; nothing silently expands.
#[derive(Debug, Clone, Copy)]
pub struct SparseBucklingLimits {
    /// Maximum number of unconstrained coordinates.
    pub max_dofs: usize,
    /// Combined stored entries of the two input matrices.
    pub max_matrix_entries: usize,
    /// Native sparse symbolic/front/work limits.
    pub factor: DirectLimits,
    /// Numeric factor plus scratch payload bytes (not process RSS).
    pub max_factor_bytes: usize,
    /// Total whitened operator applications across all outer iterations.
    pub max_operator_applications: usize,
    /// Conservative scalar-operation cap per whitened operator application.
    pub max_operator_work: usize,
}

impl Default for SparseBucklingLimits {
    fn default() -> Self {
        Self {
            max_dofs: 100_000,
            max_matrix_entries: 4_000_000,
            factor: DirectLimits {
                max_symbolic_entries: 4_000_000,
                max_analysis_work: 400_000_000,
                max_numeric_work: 4_000_000_000,
                max_front_dim: 512,
            },
            max_factor_bytes: 256 * 1024 * 1024,
            max_operator_applications: 100_000,
            max_operator_work: 32_000_000,
        }
    }
}

/// Solve the lowest positive load factors using sparse positive-factor
/// whitening and the existing generic LOBPCG service. This explicit entry is
/// also available for small pencils, for comparison with the legacy dense
/// fixture route. Requires room for a `3*min(nev,n)` native search space.
///
/// # Errors
/// Invalid input or resource exhaustion returns `InvalidInput`. A singular or
/// indefinite base, non-finite factor, or unresolved eigenpair returns
/// `SolveFailed`. Cancellation returns `Cancelled`, without partial modes.
pub fn buckling_loads(
    k: &Csr,
    kg: &Csr,
    dof_map: &[Option<usize>],
    nev: usize,
    steps: usize,
    limits: &SparseBucklingLimits,
) -> Result<BucklingResult, SolidError> {
    buckling_loads_cancellable(k, kg, dof_map, nev, steps, limits, &|| false)
}

/// Cancellable counterpart with polling during symbolic/numeric factorization,
/// triangular actions, sparse matrix rows and each generic spectral tick.
#[allow(clippy::too_many_lines)]
pub fn buckling_loads_cancellable(
    k: &Csr,
    kg: &Csr,
    dof_map: &[Option<usize>],
    nev: usize,
    steps: usize,
    limits: &SparseBucklingLimits,
    cancelled: &impl Fn() -> bool,
) -> Result<BucklingResult, SolidError> {
    check_cancelled(cancelled)?;
    let n = k.nrows();
    if n == 0
        || n > limits.max_dofs
        || k.ncols() != n
        || kg.nrows() != n
        || kg.ncols() != n
        || k.nnz()
            .checked_add(kg.nnz())
            .is_none_or(|v| v > limits.max_matrix_entries)
        || dof_map.len() > limits.max_dofs.saturating_mul(8)
    {
        return Err(invalid(
            "sparse buckling pencil dimensions or input storage exceed admission",
        ));
    }
    if nev == 0
        || steps == 0
        || limits.max_factor_bytes == 0
        || limits.max_operator_applications == 0
        || limits.max_operator_work == 0
        || limits.factor.max_front_dim == 0
    {
        return Err(invalid(
            "sparse buckling needs positive mode, iteration and resource allowances",
        ));
    }
    let count = nev.min(n);
    // This native constructor admits all n*b, n*(3b), reduced-square and
    // unpolled scalar-work products before allocating any search block.
    let mut service = EigenService::new(
        EigenBackend::Lobpcg,
        n,
        EigenQuery {
            k: count,
            largest: true,
            tol: 1e-10,
            steps_per_tick: 1,
        },
    )
    .map_err(|error| spectral_error(error, 0))?;
    let mut seen = vec![false; n];
    for (slot, entry) in dof_map.iter().enumerate() {
        if slot % 256 == 0 {
            check_cancelled(cancelled)?;
        }
        if let Some(free) = *entry {
            if free >= n || seen[free] {
                return Err(invalid(
                    "buckling DOF map must cover each free DOF exactly once",
                ));
            }
            seen[free] = true;
        }
    }
    if seen.contains(&false) {
        return Err(invalid("buckling DOF map is missing a free DOF"));
    }
    let (k_scale, scaled_k) = scaled_symmetric(k, cancelled)?;
    let (g_scale, scaled_g) = scaled_symmetric(kg, cancelled)?;
    if k_scale == 0.0 {
        return Err(failed(0, f64::INFINITY));
    }
    // A front strictly smaller than n prevents a global n-by-n production
    // allocation, even for a small input which happens to be fully dense.
    let mut direct_limits = limits.factor;
    direct_limits.max_front_dim = direct_limits.max_front_dim.min(n.saturating_sub(1));
    let symbolic =
        SymbolicLdlt::analyze_bounded(&scaled_k, DirectOrdering::Amd, &direct_limits, cancelled)
            .map_err(direct_error)?;
    let options = LdltOptions {
        allow_2x2: false,
        max_factor_bytes: Some(limits.max_factor_bytes),
        ..LdltOptions::default()
    };
    let factor = symbolic
        .factor_bounded(&scaled_k, &options, &direct_limits, cancelled)
        .map_err(direct_error)?;
    let whitening = PositiveWhitening::new(&factor, cancelled).map_err(positive_error)?;
    if g_scale == 0.0 {
        check_cancelled(cancelled)?;
        return Ok(BucklingResult {
            loads: Vec::new(),
            modes: Vec::new(),
            relative_residuals: Vec::new(),
            dof_map: dof_map.to_vec(),
            iters: 0,
        });
    }
    // Two triangular actions, one sparse product, copies and finite checks.
    let operator_work = factor
        .stats()
        .nnz_l
        .saturating_mul(8)
        .saturating_add(scaled_g.nnz().saturating_mul(2))
        .saturating_add(n.saturating_mul(32));
    if operator_work > limits.max_operator_work {
        return Err(invalid(
            "sparse buckling operator exceeds its per-application work allowance",
        ));
    }
    let op = BucklingOp {
        whitening: &whitening,
        geometric: &scaled_g,
        cancelled,
        scratch: RefCell::new((vec![0.0; n], vec![0.0; n], vec![0.0; n])),
        applications: Cell::new(0),
        max_applications: limits.max_operator_applications,
        failure: RefCell::new(None),
    };
    let mut last_error = failed(0, f64::INFINITY);
    for iteration in 0..steps {
        check_cancelled(cancelled)?;
        let progress = service.tick_cancellable(&op, cancelled);
        // Preserve the domain's typed budget/cancellation error. An aborted
        // infallible SymmetricOp action never publishes its placeholder image.
        if let Some(error) = op.failure.borrow_mut().take() {
            return Err(error);
        }
        let progress = progress.map_err(|error| spectral_error(error, iteration))?;
        let resolved = progress
            .pairs
            .iter()
            .all(|pair| relative_transformed(pair) <= 1e-10);
        if resolved || iteration + 1 == steps {
            match physical_result(
                &progress.pairs,
                &whitening,
                &scaled_k,
                &scaled_g,
                k_scale,
                g_scale,
                dof_map,
                progress.ticks,
                cancelled,
            ) {
                Ok(result) => return Ok(result),
                Err(error @ SolidError::SolveFailed { .. }) => last_error = error,
                Err(error) => return Err(error),
            }
        }
    }
    Err(last_error)
}

struct BucklingOp<'a, 'f, F: Fn() -> bool> {
    whitening: &'a PositiveWhitening<'f>,
    geometric: &'a Csr,
    cancelled: &'a F,
    scratch: RefCell<(Vec<f64>, Vec<f64>, Vec<f64>)>,
    applications: Cell<usize>,
    max_applications: usize,
    failure: RefCell<Option<SolidError>>,
}

impl<F: Fn() -> bool> SymmetricOp for BucklingOp<'_, '_, F> {
    fn dim(&self) -> usize {
        self.whitening.dim()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        if self.failure.borrow().is_some() {
            y.fill(0.0);
            return;
        }
        if self.applications.get() >= self.max_applications {
            *self.failure.borrow_mut() = Some(invalid(
                "sparse buckling exhausted its total operator-application allowance",
            ));
            y.fill(0.0);
            return;
        }
        self.applications.set(self.applications.get() + 1);
        let result = (|| {
            let mut scratch = self.scratch.borrow_mut();
            let (physical, force, triangular) = &mut *scratch;
            self.whitening
                .apply(x, physical, triangular, self.cancelled)
                .map_err(positive_error)?;
            self.geometric
                .spmv_with_checkpoint(physical, force, || check_cancelled(self.cancelled))
                .map_err(sparse_apply_error)?;
            for value in force.iter_mut() {
                *value = -*value;
            }
            self.whitening
                .apply_transpose(force, y, triangular, self.cancelled)
                .map_err(positive_error)
        })();
        if let Err(error) = result {
            *self.failure.borrow_mut() = Some(error);
            y.fill(0.0);
        }
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn physical_result(
    pairs: &[CertifiedEigenvalue],
    whitening: &PositiveWhitening<'_>,
    k: &Csr,
    g: &Csr,
    k_scale: f64,
    g_scale: f64,
    dof_map: &[Option<usize>],
    iters: usize,
    cancelled: &impl Fn() -> bool,
) -> Result<BucklingResult, SolidError> {
    let n = k.nrows();
    let mut loads = Vec::new();
    let mut modes = Vec::new();
    let mut relative_residuals = Vec::new();
    let mut k_phi = vec![0.0; n];
    let mut g_phi = vec![0.0; n];
    let mut triangular = vec![0.0; n];
    // Generic service presents ascending theta; physical loads are the
    // reciprocals of its descending positive roots.
    for pair in pairs.iter().rev() {
        check_cancelled(cancelled)?;
        let theta = pair.value();
        let transformed = relative_transformed(pair);
        if !transformed.is_finite() || transformed > 1e-8 {
            return Err(failed(iters, transformed));
        }
        if theta <= 0.0 {
            continue;
        }
        // The residual must describe the rounded physical load we return,
        // including subnormal restoration and extreme scale ratios.
        let (load, (ks, gs)) = super::restored_load(k_scale, g_scale, theta)?;
        let mut phi = vec![0.0; n];
        whitening
            .apply(pair.vector(), &mut phi, &mut triangular, cancelled)
            .map_err(positive_error)?;
        let amplitude = phi.iter().fold(0.0_f64, |scale, v| scale.max(v.abs()));
        if !amplitude.is_finite() || amplitude <= 0.0 {
            return Err(invalid("buckling mode is not representable"));
        }
        for value in &mut phi {
            *value /= amplitude;
        }
        k.spmv_with_checkpoint(&phi, &mut k_phi, || check_cancelled(cancelled))
            .map_err(sparse_apply_error)?;
        g.spmv_with_checkpoint(&phi, &mut g_phi, || check_cancelled(cancelled))
            .map_err(sparse_apply_error)?;
        let (mut residual, mut knorm, mut gnorm) = (0.0_f64, 0.0_f64, 0.0_f64);
        for i in 0..n {
            residual = residual.hypot(ks.mul_add(k_phi[i], gs * g_phi[i]));
            knorm = knorm.hypot(k_phi[i]);
            gnorm = gnorm.hypot(g_phi[i]);
        }
        let relative = residual / (ks * knorm + gs * gnorm);
        if !relative.is_finite() || relative > 1e-8 {
            return Err(failed(iters, relative));
        }
        let work = -phi.iter().zip(&g_phi).map(|(a, b)| a * b).sum::<f64>();
        let normalization = work.sqrt() * g_scale.sqrt();
        if !normalization.is_finite() || normalization <= 0.0 {
            return Err(invalid(
                "buckling mode has no representable positive geometric work",
            ));
        }
        for value in &mut phi {
            *value /= normalization;
        }
        if phi.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "work-normalized buckling mode is not representable",
            ));
        }
        loads.push(load);
        modes.push(phi);
        relative_residuals.push(relative);
    }
    check_cancelled(cancelled)?;
    Ok(BucklingResult {
        loads,
        modes,
        relative_residuals,
        dof_map: dof_map.to_vec(),
        iters,
    })
}

fn relative_transformed(pair: &CertifiedEigenvalue) -> f64 {
    if pair.value() == 0.0 && pair.residual() == 0.0 {
        0.0
    } else {
        pair.residual() / pair.value().abs()
    }
}
fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}
fn failed(iters: usize, rel_residual: f64) -> SolidError {
    SolidError::SolveFailed {
        iters,
        rel_residual,
    }
}
fn direct_error(error: LdltError) -> SolidError {
    match error {
        LdltError::Cancelled => SolidError::Cancelled,
        LdltError::PivotBreakdown { .. } => failed(0, f64::INFINITY),
        other => invalid(&format!("sparse buckling factor admission: {other}")),
    }
}
fn positive_error(error: PositiveFactorError) -> SolidError {
    match error {
        PositiveFactorError::Cancelled => SolidError::Cancelled,
        PositiveFactorError::NotPositive | PositiveFactorError::NonFinite => {
            failed(0, f64::INFINITY)
        }
        other => invalid(&format!("sparse buckling whitening: {other}")),
    }
}
fn spectral_error(error: ServiceError, iters: usize) -> SolidError {
    match error {
        ServiceError::Cancelled => SolidError::Cancelled,
        ServiceError::NonFiniteOperator | ServiceError::OperatorFault => {
            failed(iters, f64::INFINITY)
        }
        ServiceError::Unconverged {
            ticks,
            worst_residual,
        } => failed(ticks, worst_residual),
        other => invalid(&format!("sparse buckling spectral admission: {other}")),
    }
}

fn sparse_apply_error(error: SparseApplyError<SolidError>) -> SolidError {
    match error {
        SparseApplyError::Interrupted(error) => error,
        other => invalid(&format!("sparse buckling product admission: {other}")),
    }
}
