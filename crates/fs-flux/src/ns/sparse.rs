//! Sparse velocity/pressure solves for the assembled BDM saddle system.
//!
//! Approximate block LDU: A_u^{-1} is ILU(0), and the pressure factor uses
//! S_hat = C - B diag(A_u)^{-1} G. These blocks come from the actual operator,
//! including variable viscosity, transient mass, convection and constraints.
//! The approximation changes the preconditioner only, never the flow equations.
//! See Benzi, Golub & Liesen (2005), Acta Numerica, sections 3 and 10:
//! https://doi.org/10.1017/S0962492904000212.

use super::{FluxParams, FluxSolution, FluxSystem};
use fs_math::dd::Dd;
use fs_solver::{FgmresState, op::CsrOp};
use fs_sparse::precond::{Ilu0, Precond, ilu0};
use fs_sparse::{Coo, Csr};
use std::cell::RefCell;

fn row_norm(matrix: &Csr, row: usize) -> f64 {
    matrix
        .row(row)
        .1
        .iter()
        .fold(0.0_f64, |norm, value| norm.max(value.abs()))
        .max(f64::MIN_POSITIVE)
}

/// Both diagnostics use the ORIGINAL matrix, including after equilibration.
/// Row normalization keeps large viscous loads from concealing errors in the
/// continuity, pressure-pin and strong boundary equations.
pub(super) fn residuals(matrix: &Csr, rhs: &[f64], x: &[f64], work: &mut [f64]) -> (f64, f64) {
    matrix.spmv(x, work);
    let (mut defect, mut norm, mut scaled_defect, mut scaled_norm) =
        (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    for row in 0..rhs.len() {
        let scale = row_norm(matrix, row);
        let residual = work[row] - rhs[row];
        defect = defect.hypot(residual);
        norm = norm.hypot(rhs[row]);
        scaled_defect = scaled_defect.hypot(residual / scale);
        scaled_norm = scaled_norm.hypot(rhs[row] / scale);
    }
    (
        defect / norm.max(f64::MIN_POSITIVE),
        scaled_defect / scaled_norm.max(f64::MIN_POSITIVE),
    )
}

fn equilibrate(matrix: &Csr, rhs: &[f64]) -> (Csr, Vec<f64>) {
    let mut row_ptr = Vec::with_capacity(matrix.nrows() + 1);
    let mut columns = Vec::with_capacity(matrix.nnz());
    let mut values = Vec::with_capacity(matrix.nnz());
    let mut scaled_rhs = Vec::with_capacity(rhs.len());
    row_ptr.push(0);
    for (row, &right) in rhs.iter().enumerate() {
        let scale = row_norm(matrix, row);
        let (indices, coefficients) = matrix.row(row);
        columns.extend_from_slice(indices);
        values.extend(coefficients.iter().map(|value| value / scale));
        row_ptr.push(columns.len());
        scaled_rhs.push(right / scale);
    }
    (
        Csr::from_parts(matrix.nrows(), matrix.ncols(), row_ptr, columns, values),
        scaled_rhs,
    )
}

/// Translate the pressure unknowns, without changing the assembled operator
/// or the physical traction datum. A large common pressure is an affine load,
/// not a suitable convergence scale for the much smaller flow-driving load.
/// Each connected component gets its own origin; a prescribed pressure pin
/// takes precedence over any momentum row. The largest nonzero constant-
/// pressure response gives a well-scaled origin without a normal equation.
pub(super) fn pressure_origin(
    system: &FluxSystem<'_>,
    matrix: &Csr,
    rhs: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let count = system.pressure_roots.len();
    let mut offsets = vec![0.0; count];
    let mut pinned = vec![false; count];
    let mut largest_response = vec![0.0; count];
    for row in system.n_u..system.n {
        let (columns, values) = matrix.row(row);
        if columns == [row] && values == [1.0] {
            let component = system.pressure_component[row - system.n_u];
            if !pinned[component] {
                offsets[component] = rhs[row];
                pinned[component] = true;
            }
        }
    }
    for (row, &right) in rhs.iter().enumerate() {
        let (columns, values) = matrix.row(row);
        let mut component = None;
        let mut response = Dd::from_f64(0.0);
        let mut mixed = false;
        for (&column, &value) in columns.iter().zip(values) {
            if column < system.n_u {
                continue;
            }
            let current = system.pressure_component[column - system.n_u];
            if component.is_some_and(|previous| previous != current) {
                mixed = true;
                break;
            }
            component = Some(current);
            response = response + Dd::from_f64(value);
        }
        let Some(component) = component else { continue };
        let response = response.to_f64();
        if !mixed && !pinned[component] && response.abs() > largest_response[component] {
            let offset = right / response;
            if offset.is_finite() {
                offsets[component] = offset;
                largest_response[component] = response.abs();
            }
        }
    }
    if offsets.iter().all(|offset| *offset == 0.0) {
        // Preserve the original load bits and arithmetic when no translation
        // is needed, including the usual zero-pinned closed components.
        return (rhs.to_vec(), offsets);
    }
    let shifted_rhs: Vec<_> = rhs
        .iter()
        .enumerate()
        .map(|(row, &right)| {
            let (columns, values) = matrix.row(row);
            let mut shifted = Dd::from_f64(right);
            for (&column, &value) in columns.iter().zip(values) {
                if column >= system.n_u {
                    shifted = shifted
                        - Dd::from_f64(value)
                            * Dd::from_f64(offsets[system.pressure_component[column - system.n_u]]);
                }
            }
            shifted.to_f64()
        })
        .collect();
    if shifted_rhs.iter().any(|value| !value.is_finite()) {
        // An optional origin must never make an otherwise representable
        // original system overflow. Retain the original coordinates then.
        return (rhs.to_vec(), vec![0.0; count]);
    }
    (shifted_rhs, offsets)
}

/// Recheck the pressure values that actually survive datum restoration.
/// Both linear-solver routes reuse an existing work buffer for these returned
/// pressure differences; the original operator and translated load are fixed.
pub(super) fn returned_pressure_residuals(
    system: &FluxSystem<'_>,
    matrix: &Csr,
    shifted_rhs: &[f64],
    pressure_offsets: &[f64],
    x: &[f64],
    returned_shifted: &mut [f64],
    work: &mut [f64],
) -> (f64, f64) {
    returned_shifted.copy_from_slice(x);
    for (cell, &component) in system.pressure_component.iter().enumerate() {
        returned_shifted[system.n_u + cell] -= pressure_offsets[component];
    }
    residuals(matrix, shifted_rhs, returned_shifted, work)
}

struct BlockWork {
    velocity: Vec<f64>,
    pressure_rhs: Vec<f64>,
    gradient: Vec<f64>,
    correction: Vec<f64>,
}

struct FlowBlocks {
    velocity: Ilu0,
    pressure: Ilu0,
    divergence: Csr,
    gradient: Csr,
    work: RefCell<BlockWork>,
}

impl FlowBlocks {
    fn new(matrix: &Csr, n_u: usize) -> Option<Self> {
        let n_p = matrix.nrows() - n_u;
        let mut velocity = Coo::new(n_u, n_u);
        let mut divergence = Coo::new(n_p, n_u);
        let mut gradient = Coo::new(n_u, n_p);
        let mut schur = Coo::new(n_p, n_p);
        let mut diagonal = vec![0.0; n_u];
        for row in 0..matrix.nrows() {
            let (columns, values) = matrix.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                if !value.is_finite() {
                    return None;
                }
                match (row < n_u, column < n_u) {
                    (true, true) => {
                        velocity.push(row, column, value);
                        if row == column {
                            diagonal[row] = value;
                        }
                    }
                    (false, true) => divergence.push(row - n_u, column, value),
                    (true, false) => gradient.push(row, column - n_u, value),
                    (false, false) => schur.push(row - n_u, column - n_u, value),
                }
            }
        }
        if diagonal.iter().any(|d| *d == 0.0) {
            return None;
        }
        let divergence = divergence.assemble();
        let gradient = gradient.assemble();
        for row in 0..n_p {
            let (columns, values) = divergence.row(row);
            for (&velocity_column, &b) in columns.iter().zip(values) {
                let (pressure_columns, gradient_values) = gradient.row(velocity_column);
                for (&column, &g) in pressure_columns.iter().zip(gradient_values) {
                    let value = -b * (g / diagonal[velocity_column]);
                    if !value.is_finite() {
                        return None;
                    }
                    schur.push(row, column, value);
                }
            }
        }
        Some(Self {
            velocity: ilu0(&velocity.assemble()).ok()?,
            pressure: ilu0(&schur.assemble()).ok()?,
            divergence,
            gradient,
            work: RefCell::new(BlockWork {
                velocity: vec![0.0; n_u],
                pressure_rhs: vec![0.0; n_p],
                gradient: vec![0.0; n_u],
                correction: vec![0.0; n_u],
            }),
        })
    }
}

struct FlowPreconditioner(Option<FlowBlocks>);

impl Precond for FlowPreconditioner {
    fn apply(&self, residual: &[f64], output: &mut [f64]) {
        let Some(blocks) = &self.0 else {
            // Breakdown of an approximate factor does not alter the physical
            // operator. The original residual still decides convergence.
            output.copy_from_slice(residual);
            return;
        };
        let n_u = blocks.gradient.nrows();
        let mut work = blocks.work.borrow_mut();
        let BlockWork {
            velocity,
            pressure_rhs,
            gradient,
            correction,
        } = &mut *work;
        blocks.velocity.apply(&residual[..n_u], velocity);
        blocks.divergence.spmv(velocity, pressure_rhs);
        for (value, rhs) in pressure_rhs.iter_mut().zip(&residual[n_u..]) {
            *value = rhs - *value;
        }
        blocks.pressure.apply(pressure_rhs, &mut output[n_u..]);
        blocks.gradient.spmv(&output[n_u..], gradient);
        blocks.velocity.apply(gradient, correction);
        for i in 0..n_u {
            output[i] = velocity[i] - correction[i];
        }
    }
}

impl FluxSystem<'_> {
    /// Sparse block-preconditioned solve, including below the automatic dense
    /// cutoff. Velocity ILU(0) and a sparse approximate pressure Schur factor
    /// are built from this exact assembled operator, so the preconditioner
    /// follows local viscosity, transient mass and convection automatically.
    /// Per-application block work arrays are allocated once and reused.
    ///
    /// Rows and loads are equilibrated by each row's largest absolute entry.
    /// A constant pressure origin is removed independently in each connected
    /// component before iteration and restored in the returned field. The
    /// translated load measures flow accuracy independently of a large ambient
    /// traction pressure. This is an affine change of variables, not a new
    /// pressure boundary condition or a change to the physical pressure datum.
    /// The returned pressure is translated back again for a final defect check:
    /// rounding away a necessary pressure difference during restoration cannot
    /// produce a converged result.
    /// Both the original and row-scaled residuals are independently recomputed
    /// against the original matrix after every restart. Large viscous loads
    /// therefore cannot hide the continuity and strong boundary equations.
    /// Use `FluxSolution::converged` to require both criteria on return.
    /// `max_iters` bounds actual Krylov iterations, including a short final
    /// restart. Zero budget cannot reduce a nonzero translated load. An
    /// incomplete-factorization breakdown falls back to identity
    /// preconditioning; reported residuals never claim unverified convergence.
    /// No dense matrix or dense pressure Schur complement is constructed.
    ///
    /// # Panics
    /// If dimensions disagree or tolerance is not finite and positive.
    #[must_use]
    pub fn solve_sparse(&self, matrix: &Csr, rhs: &[f64], params: FluxParams) -> FluxSolution {
        assert_eq!(matrix.nrows(), self.n, "flow matrix row count mismatch");
        assert_eq!(matrix.ncols(), self.n, "flow matrix column count mismatch");
        assert_eq!(rhs.len(), self.n, "flow load count mismatch");
        assert!(
            params.tol.is_finite() && params.tol > 0.0,
            "positive finite tolerance required"
        );
        let (shifted_rhs, pressure_offsets) = pressure_origin(self, matrix, rhs);
        let (scaled_matrix, scaled_rhs) = equilibrate(matrix, &shifted_rhs);
        let preconditioner = FlowPreconditioner(FlowBlocks::new(&scaled_matrix, self.n_u));
        let operator = CsrOp::general(scaled_matrix);
        let mut state = FgmresState::new(&scaled_rhs, 80);
        let mut work = vec![0.0; self.n];
        let (mut original_residual, mut scaled_residual) =
            residuals(matrix, &shifted_rhs, &state.x, &mut work);
        let mut inner_tolerance = params.tol;
        while state.iters < params.max_iters
            && (original_residual > params.tol || scaled_residual > params.tol)
        {
            state.restart = 80.min(params.max_iters - state.iters);
            let before = state.iters;
            let report = state.run(&operator, &preconditioner, &scaled_rhs, inner_tolerance, 1);
            (original_residual, scaled_residual) =
                residuals(matrix, &shifted_rhs, &state.x, &mut work);
            if state.iters == before
                || report.diagnosis == Some(fs_solver::StallDiagnosis::Breakdown)
            {
                break;
            }
            // Equilibrated convergence can precede convergence in the original
            // norm (or its independently rounded residual). Continue refining
            // rather than publishing either criterion as the other one.
            if report.converged_euclidean() && original_residual.max(scaled_residual) > params.tol {
                inner_tolerance = inner_tolerance.min(
                    0.5 * state.rel_residual() * params.tol
                        / original_residual.max(scaled_residual),
                );
            }
        }
        let shifted_residual = original_residual.max(scaled_residual);
        for (cell, &component) in self.pressure_component.iter().enumerate() {
            state.x[self.n_u + cell] += pressure_offsets[component];
        }
        // Recheck the ACTUALLY RETURNED pressure, not only the more accurately
        // resolved internal pressure differences. Restoring a large datum can
        // round a small gradient away. The original-load norm would conceal
        // that defect, so evaluate the same translated equations again. Reuse
        // the now-unused equilibrated load instead of allocating another field.
        let mut returned_shifted = scaled_rhs;
        let (returned_residual, returned_scaled_residual) = returned_pressure_residuals(
            self,
            matrix,
            &shifted_rhs,
            &pressure_offsets,
            &state.x,
            &mut returned_shifted,
            &mut work,
        );
        let (original_residual, scaled_residual) = residuals(matrix, rhs, &state.x, &mut work);
        FluxSolution {
            x: state.x,
            n_u: self.n_u,
            iters: state.iters,
            rel_residual: original_residual,
            scaled_rel_residual: scaled_residual
                .max(shifted_residual)
                .max(returned_residual)
                .max(returned_scaled_residual),
        }
    }
}
