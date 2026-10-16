//! The plastic-design layout LP with a first-order primal-dual solver.
//!
//! Variables: per member, split tension/compression forces
//! `q⁺, q⁻ ≥ 0`. Objective: material volume `Σ lᵢ(qᵢ⁺ + qᵢ⁻)/σ_y`.
//! Constraints: nodal equilibrium `B·(q⁺ − q⁻) = f` on FREE degrees of
//! freedom (supports drop their rows). Standard form:
//! `min cᵀx  s.t.  A x = b, x ≥ 0`.
//!
//! PDHG (Chambolle–Pock): `x ← Π₊(x − τ(c + Aᵀy))`,
//! `y ← y + σ(A(2x − x_prev) − b)`, with `τσ‖A‖² < 1` from a
//! power-iteration norm estimate. Sparse-matvec dominated (fs-sparse
//! CSR), bitwise deterministic, warm-startable across load cases. Relative
//! primal/dual objective separation and equilibrium residual are tracked at
//! every check interval. These are convergence diagnostics, not a certified
//! optimum interval: the returned primal is only approximately equilibrated,
//! and the floating dual scaling is not outward-verified.

use crate::ground::{GroundStructure, TrussConstructionError};
use fs_blake3::{Blake3, hash_domain};
use fs_exec::Cx;
use fs_sparse::Csr;
use std::convert::Infallible;
use std::fmt::Write as _;
use std::mem::size_of;

/// Hard ceiling for free equilibrium degrees of freedom in one layout LP.
pub const HARD_MAX_LAYOUT_FREE_DOFS: usize = 1_048_576;
/// Hard ceiling for split tension/compression variables in one layout LP.
pub const HARD_MAX_LAYOUT_VARIABLES: usize = 2_000_000;
/// Hard ceiling for COO triplets staged while assembling one layout LP.
pub const HARD_MAX_LAYOUT_STAGED_TRIPLETS: usize = 8_000_000;
/// Hard ceiling for conservatively estimated retained LP storage.
pub const HARD_MAX_LAYOUT_RETAINED_BYTES: usize = 512 * 1024 * 1024;

const LAYOUT_POLL_STRIDE: usize = 256;
const MIN_LAYOUT_LOAD_NORM_SQUARED: f64 = 1e-60;
const PDHG_REPORT_SCALAR_SNAPSHOT_DOMAIN: &str =
    "frankensim.fs-truss.pdhg-report-scalar-snapshot.v1";
const PDHG_REPORT_TRACE_SNAPSHOT_DOMAIN: &str = "frankensim.fs-truss.pdhg-report-trace-snapshot.v1";

fn stable_usize_bytes(value: usize) -> [u8; 8] {
    // Solver-minted values are far below u64::MAX on every admitted target.
    // Saturation keeps hostile future 128-bit mutations deterministic and
    // unable to reproduce any valid solver snapshot.
    u64::try_from(value).unwrap_or(u64::MAX).to_le_bytes()
}

/// Immutable support mask and nodal loads admitted for one ground structure.
#[derive(Debug, Clone, PartialEq)]
pub struct LayoutCase {
    supported: Vec<[bool; 2]>,
    loads: Vec<[f64; 2]>,
}

impl LayoutCase {
    /// Admit one support/load case for an expected number of nodes.
    ///
    /// # Errors
    /// Returns a structured refusal when either vector has the wrong length,
    /// the expected node count is zero, or any load component is non-finite.
    pub fn try_new(
        supported: Vec<[bool; 2]>,
        loads: Vec<[f64; 2]>,
        expected_nodes: usize,
    ) -> Result<Self, TrussConstructionError> {
        if expected_nodes == 0 {
            return Err(TrussConstructionError::InvalidInput {
                field: "expected_nodes",
                requirement: "must be at least one",
            });
        }
        if supported.len() != expected_nodes {
            return Err(TrussConstructionError::VectorLength {
                field: "supported",
                expected: expected_nodes,
                actual: supported.len(),
            });
        }
        if loads.len() != expected_nodes {
            return Err(TrussConstructionError::VectorLength {
                field: "loads",
                expected: expected_nodes,
                actual: loads.len(),
            });
        }
        if loads
            .iter()
            .flatten()
            .any(|component| !component.is_finite())
        {
            return Err(TrussConstructionError::InvalidInput {
                field: "loads",
                requirement: "must contain only finite components",
            });
        }
        Ok(Self { supported, loads })
    }

    /// Per-node support flags, indexed as `[x, y]`.
    #[must_use]
    pub fn supported(&self) -> &[[bool; 2]] {
        &self.supported
    }

    /// Per-node load vectors, indexed as `[x, y]`.
    #[must_use]
    pub fn loads(&self) -> &[[f64; 2]] {
        &self.loads
    }

    /// Number of nodes represented by this admitted case.
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.loads.len()
    }
}

/// Per-call resource limits for deterministic layout construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LayoutLimits {
    free_dofs: usize,
    variables: usize,
    staged_triplets: usize,
    retained_bytes: usize,
}

impl LayoutLimits {
    /// Admit caller-selected limits below the crate's hard ceilings.
    ///
    /// # Errors
    /// Returns a structured refusal when any limit is zero or exceeds its hard
    /// ceiling.
    pub fn try_new(
        max_free_dofs: usize,
        max_variables: usize,
        max_staged_triplets: usize,
        max_retained_bytes: usize,
    ) -> Result<Self, TrussConstructionError> {
        validate_limit("max_free_dofs", max_free_dofs, HARD_MAX_LAYOUT_FREE_DOFS)?;
        validate_limit("max_variables", max_variables, HARD_MAX_LAYOUT_VARIABLES)?;
        validate_limit(
            "max_staged_triplets",
            max_staged_triplets,
            HARD_MAX_LAYOUT_STAGED_TRIPLETS,
        )?;
        validate_limit(
            "max_retained_bytes",
            max_retained_bytes,
            HARD_MAX_LAYOUT_RETAINED_BYTES,
        )?;
        Ok(Self {
            free_dofs: max_free_dofs,
            variables: max_variables,
            staged_triplets: max_staged_triplets,
            retained_bytes: max_retained_bytes,
        })
    }

    /// Maximum admitted free equilibrium degrees of freedom.
    #[must_use]
    pub const fn max_free_dofs(self) -> usize {
        self.free_dofs
    }

    /// Maximum admitted split tension/compression variables.
    #[must_use]
    pub const fn max_variables(self) -> usize {
        self.variables
    }

    /// Maximum admitted staged COO triplets.
    #[must_use]
    pub const fn max_staged_triplets(self) -> usize {
        self.staged_triplets
    }

    /// Maximum conservative retained-storage estimate in bytes.
    #[must_use]
    pub const fn max_retained_bytes(self) -> usize {
        self.retained_bytes
    }
}

impl Default for LayoutLimits {
    fn default() -> Self {
        Self {
            free_dofs: HARD_MAX_LAYOUT_FREE_DOFS,
            variables: HARD_MAX_LAYOUT_VARIABLES,
            staged_triplets: HARD_MAX_LAYOUT_STAGED_TRIPLETS,
            retained_bytes: HARD_MAX_LAYOUT_RETAINED_BYTES,
        }
    }
}

fn validate_limit(
    field: &'static str,
    requested: usize,
    hard_ceiling: usize,
) -> Result<(), TrussConstructionError> {
    if requested == 0 || requested > hard_ceiling {
        return Err(TrussConstructionError::InvalidInput {
            field,
            requirement: "must be positive and no greater than its hard ceiling",
        });
    }
    Ok(())
}

fn poll(cx: &Cx<'_>, stage: &'static str) -> Result<(), TrussConstructionError> {
    cx.checkpoint()
        .map_err(|_| TrussConstructionError::Cancelled { stage })
}

fn empty_with_capacity<T>(
    requested: usize,
    resource: &'static str,
) -> Result<Vec<T>, TrussConstructionError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(requested)
        .map_err(|_| TrussConstructionError::AllocationFailed {
            resource,
            requested,
        })?;
    Ok(values)
}

fn zeroed_f64(
    requested: usize,
    resource: &'static str,
    cx: &Cx<'_>,
    stage: &'static str,
) -> Result<Vec<f64>, TrussConstructionError> {
    let mut values = empty_with_capacity(requested, resource)?;
    for index in 0..requested {
        if index % LAYOUT_POLL_STRIDE == 0 {
            poll(cx, stage)?;
        }
        values.push(0.0);
    }
    Ok(values)
}

fn zeroed_usize(
    requested: usize,
    resource: &'static str,
    cx: &Cx<'_>,
    stage: &'static str,
) -> Result<Vec<usize>, TrussConstructionError> {
    let mut values = empty_with_capacity(requested, resource)?;
    for index in 0..requested {
        if index % LAYOUT_POLL_STRIDE == 0 {
            poll(cx, stage)?;
        }
        values.push(0);
    }
    Ok(values)
}

fn copied_usize(
    source: &[usize],
    resource: &'static str,
    cx: &Cx<'_>,
    stage: &'static str,
) -> Result<Vec<usize>, TrussConstructionError> {
    let mut values = empty_with_capacity(source.len(), resource)?;
    for chunk in source.chunks(LAYOUT_POLL_STRIDE) {
        poll(cx, stage)?;
        values.extend_from_slice(chunk);
    }
    Ok(values)
}

#[allow(clippy::too_many_arguments)] // Owned CSR parts plus publication diagnostics and context.
fn publish_canonical_csr(
    rows: usize,
    columns: usize,
    row_ptr: Vec<usize>,
    col_idx: Vec<usize>,
    values: Vec<f64>,
    field: &'static str,
    requirement: &'static str,
    stage: &'static str,
    cx: &Cx<'_>,
) -> Result<Csr, TrussConstructionError> {
    let mut validation_steps = 0usize;
    Csr::try_from_parts_with_checkpoint(rows, columns, row_ptr, col_idx, values, || {
        if validation_steps.is_multiple_of(LAYOUT_POLL_STRIDE) {
            poll(cx, stage)?;
        }
        validation_steps += 1;
        Ok(())
    })?
    .ok_or(TrussConstructionError::InvalidInput { field, requirement })
}

fn checked_product(a: usize, b: usize) -> Option<usize> {
    a.checked_mul(b)
}

fn checked_sum(parts: &[usize]) -> Option<usize> {
    parts
        .iter()
        .try_fold(0usize, |total, part| total.checked_add(*part))
}

fn retained_layout_bytes(
    node_dofs: usize,
    free_dofs: usize,
    variables: usize,
    staged_triplets: usize,
) -> Option<usize> {
    let dof_map = checked_product(node_dofs, size_of::<Option<usize>>())?;
    let rhs = checked_product(free_dofs, size_of::<f64>())?;
    let costs = checked_product(variables, size_of::<f64>())?;
    let sparse_entries = checked_product(
        staged_triplets,
        size_of::<usize>().checked_add(size_of::<f64>())?,
    )?;
    let a_rows = checked_product(free_dofs.checked_add(1)?, size_of::<usize>())?;
    let at_rows = checked_product(variables.checked_add(1)?, size_of::<usize>())?;
    checked_sum(&[
        dof_map,
        rhs,
        costs,
        sparse_entries,
        sparse_entries,
        a_rows,
        at_rows,
        size_of::<f64>(),
    ])
}

fn csr_is_finite_and_nonzero(
    matrix: &Csr,
    cx: &Cx<'_>,
    stage: &'static str,
) -> Result<bool, TrussConstructionError> {
    let mut nonzero = false;
    let mut visited = 0usize;
    for row in 0..matrix.nrows() {
        if row % LAYOUT_POLL_STRIDE == 0 {
            poll(cx, stage)?;
        }
        for &value in matrix.row(row).1 {
            if visited.is_multiple_of(LAYOUT_POLL_STRIDE) {
                poll(cx, stage)?;
            }
            visited += 1;
            if !value.is_finite() {
                return Ok(false);
            }
            nonzero |= value != 0.0;
        }
    }
    Ok(nonzero)
}

fn first_nonzero_column(matrix: &Csr, cx: &Cx<'_>) -> Result<usize, TrussConstructionError> {
    let mut visited = 0usize;
    for row in 0..matrix.nrows() {
        if row % LAYOUT_POLL_STRIDE == 0 {
            poll(cx, "layout fallback-column search")?;
        }
        let (columns, values) = matrix.row(row);
        for (&column, &value) in columns.iter().zip(values) {
            if visited.is_multiple_of(LAYOUT_POLL_STRIDE) {
                poll(cx, "layout fallback-column search")?;
            }
            visited += 1;
            if value != 0.0 {
                return Ok(column);
            }
        }
    }
    Err(TrussConstructionError::InvalidInput {
        field: "equilibrium matrix",
        requirement: "must expose a nonzero column for norm estimation",
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // One checked canonical CSR construction.
fn assemble_equilibrium(
    nodes: &[[f64; 2]],
    members: &[(usize, usize)],
    lengths: &[f64],
    dof_map: &[Option<usize>],
    free_dofs: usize,
    variables: usize,
    staged_triplets: usize,
    cx: &Cx<'_>,
) -> Result<Csr, TrussConstructionError> {
    let member_count = members.len();
    let mut row_counts = zeroed_usize(
        free_dofs,
        "layout row counts",
        cx,
        "layout row-count allocation",
    )?;
    let mut staged_entries_visited = 0usize;
    for (member, &(a, b)) in members.iter().enumerate() {
        if member % LAYOUT_POLL_STRIDE == 0 {
            poll(cx, "layout sparse row counting")?;
        }
        for dof in [2 * a, 2 * a + 1, 2 * b, 2 * b + 1] {
            if let Some(row) = dof_map[dof] {
                for _ in 0..2 {
                    if staged_entries_visited.is_multiple_of(LAYOUT_POLL_STRIDE) {
                        poll(cx, "layout sparse row counting")?;
                    }
                    staged_entries_visited += 1;
                }
                row_counts[row] =
                    row_counts[row]
                        .checked_add(2)
                        .ok_or(TrussConstructionError::WorkBudget {
                            resource: "layout staged triplets",
                            limit: staged_triplets,
                            observed: usize::MAX,
                        })?;
            }
        }
    }

    let mut row_ptr = zeroed_usize(
        free_dofs
            .checked_add(1)
            .ok_or(TrussConstructionError::WorkBudget {
                resource: "layout sparse row pointers",
                limit: free_dofs,
                observed: usize::MAX,
            })?,
        "layout sparse row pointers",
        cx,
        "layout row-pointer allocation",
    )?;
    let mut row_offset = 0usize;
    for (row, count) in row_counts.iter().copied().enumerate() {
        if row % LAYOUT_POLL_STRIDE == 0 {
            poll(cx, "layout sparse row offsets")?;
        }
        row_offset = row_offset
            .checked_add(count)
            .ok_or(TrussConstructionError::WorkBudget {
                resource: "layout staged triplets",
                limit: staged_triplets,
                observed: usize::MAX,
            })?;
        row_ptr[row + 1] = row_offset;
    }
    if row_ptr[free_dofs] != staged_triplets {
        return Err(TrussConstructionError::InvalidInput {
            field: "equilibrium matrix",
            requirement: "must preserve the exact admitted triplet count",
        });
    }

    let mut col_idx = zeroed_usize(
        staged_triplets,
        "layout sparse columns",
        cx,
        "layout sparse-column allocation",
    )?;
    let mut values = zeroed_f64(
        staged_triplets,
        "layout sparse values",
        cx,
        "layout sparse-value allocation",
    )?;
    let mut plus_cursor = copied_usize(
        &row_ptr[..free_dofs],
        "layout positive cursors",
        cx,
        "layout positive-cursor allocation",
    )?;
    let mut minus_cursor = empty_with_capacity(free_dofs, "layout negative cursors")?;
    for (row, &count) in row_counts.iter().enumerate() {
        if row % LAYOUT_POLL_STRIDE == 0 {
            poll(cx, "layout negative-cursor allocation")?;
        }
        minus_cursor.push(row_ptr[row] + count / 2);
    }

    staged_entries_visited = 0;
    for (member, &(a, b)) in members.iter().enumerate() {
        if member % LAYOUT_POLL_STRIDE == 0 {
            poll(cx, "layout sparse fill")?;
        }
        let dx = (nodes[b][0] - nodes[a][0]) / lengths[member];
        let dy = (nodes[b][1] - nodes[a][1]) / lengths[member];
        for (dof, value) in [(2 * a, dx), (2 * a + 1, dy), (2 * b, -dx), (2 * b + 1, -dy)] {
            if let Some(row) = dof_map[dof] {
                for _ in 0..2 {
                    if staged_entries_visited.is_multiple_of(LAYOUT_POLL_STRIDE) {
                        poll(cx, "layout sparse fill")?;
                    }
                    staged_entries_visited += 1;
                }
                let positive = plus_cursor[row];
                col_idx[positive] = member;
                values[positive] = value;
                plus_cursor[row] += 1;

                let negative = minus_cursor[row];
                col_idx[negative] = member_count + member;
                values[negative] = -value;
                minus_cursor[row] += 1;
            }
        }
    }
    for (row, ((&positive, &negative), &count)) in plus_cursor
        .iter()
        .zip(&minus_cursor)
        .zip(&row_counts)
        .enumerate()
    {
        if row % LAYOUT_POLL_STRIDE == 0 {
            poll(cx, "layout sparse fill validation")?;
        }
        if positive != row_ptr[row] + count / 2 || negative != row_ptr[row + 1] {
            return Err(TrussConstructionError::InvalidInput {
                field: "equilibrium matrix",
                requirement: "must fill every admitted sparse slot exactly once",
            });
        }
    }
    poll(cx, "layout sparse publication")?;
    publish_canonical_csr(
        free_dofs,
        variables,
        row_ptr,
        col_idx,
        values,
        "equilibrium matrix",
        "must preserve checked canonical CSR shape",
        "layout sparse canonical validation",
        cx,
    )
}

fn transpose_checked(matrix: &Csr, cx: &Cx<'_>) -> Result<Csr, TrussConstructionError> {
    let rows = matrix.nrows();
    let columns = matrix.ncols();
    let mut counts = zeroed_usize(
        columns
            .checked_add(1)
            .ok_or(TrussConstructionError::WorkBudget {
                resource: "layout transpose row pointers",
                limit: columns,
                observed: usize::MAX,
            })?,
        "layout transpose row pointers",
        cx,
        "layout transpose-count allocation",
    )?;
    let mut visited = 0usize;
    for row in 0..rows {
        if row % LAYOUT_POLL_STRIDE == 0 {
            poll(cx, "layout transpose counting")?;
        }
        for &column in matrix.row(row).0 {
            if visited.is_multiple_of(LAYOUT_POLL_STRIDE) {
                poll(cx, "layout transpose counting")?;
            }
            visited += 1;
            counts[column + 1] =
                counts[column + 1]
                    .checked_add(1)
                    .ok_or(TrussConstructionError::WorkBudget {
                        resource: "layout transpose entries",
                        limit: matrix.nnz(),
                        observed: usize::MAX,
                    })?;
        }
    }
    let mut transpose_offset = 0usize;
    for (column, count) in counts.iter_mut().enumerate().skip(1) {
        if column % LAYOUT_POLL_STRIDE == 0 {
            poll(cx, "layout transpose row offsets")?;
        }
        transpose_offset =
            transpose_offset
                .checked_add(*count)
                .ok_or(TrussConstructionError::WorkBudget {
                    resource: "layout transpose entries",
                    limit: matrix.nnz(),
                    observed: usize::MAX,
                })?;
        *count = transpose_offset;
    }
    let row_ptr = copied_usize(
        &counts,
        "layout transpose published row pointers",
        cx,
        "layout transpose row-pointer copy",
    )?;
    let mut cursor = counts;
    let mut col_idx = zeroed_usize(
        matrix.nnz(),
        "layout transpose columns",
        cx,
        "layout transpose-column allocation",
    )?;
    let mut values = zeroed_f64(
        matrix.nnz(),
        "layout transpose values",
        cx,
        "layout transpose-value allocation",
    )?;
    visited = 0;
    for row in 0..rows {
        if row % LAYOUT_POLL_STRIDE == 0 {
            poll(cx, "layout transpose fill")?;
        }
        let (source_columns, source_values) = matrix.row(row);
        for (&column, &value) in source_columns.iter().zip(source_values) {
            if visited.is_multiple_of(LAYOUT_POLL_STRIDE) {
                poll(cx, "layout transpose fill")?;
            }
            visited += 1;
            let destination = cursor[column];
            col_idx[destination] = row;
            values[destination] = value;
            cursor[column] += 1;
        }
    }
    poll(cx, "layout transpose publication")?;
    publish_canonical_csr(
        columns,
        rows,
        row_ptr,
        col_idx,
        values,
        "equilibrium transpose",
        "must preserve checked canonical CSR shape",
        "layout transpose canonical validation",
        cx,
    )
}

fn spmv_checked(
    matrix: &Csr,
    input: &[f64],
    output: &mut [f64],
    cx: &Cx<'_>,
    stage: &'static str,
) -> Result<(), TrussConstructionError> {
    if input.len() != matrix.ncols() || output.len() != matrix.nrows() {
        return Err(TrussConstructionError::InvalidInput {
            field: "layout norm multiply",
            requirement: "must have exact matrix/vector dimensions",
        });
    }
    let mut visited = 0usize;
    for (row, destination) in output.iter_mut().enumerate() {
        if row % LAYOUT_POLL_STRIDE == 0 {
            poll(cx, stage)?;
        }
        let (columns, values) = matrix.row(row);
        let mut accumulator = 0.0f64;
        for (&column, &value) in columns.iter().zip(values) {
            if visited.is_multiple_of(LAYOUT_POLL_STRIDE) {
                poll(cx, stage)?;
            }
            visited += 1;
            accumulator = value.mul_add(input[column], accumulator);
        }
        *destination = accumulator;
    }
    Ok(())
}

/// Maximum iterations admitted to one direct PDHG solve.
pub const MAX_PDHG_ITERS: usize = 1_000_000;

/// PDHG controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PdhgSettings {
    /// Iteration cap.
    pub max_iters: usize,
    /// Relative primal/dual objective-separation target.
    pub gap_tol: f64,
    /// Check/ledger interval.
    pub check_every: usize,
}

impl Default for PdhgSettings {
    fn default() -> Self {
        PdhgSettings {
            max_iters: 200_000,
            gap_tol: 1e-6,
            check_every: 500,
        }
    }
}

/// Structured refusal for invalid PDHG controls or warm-start state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PdhgError {
    /// A solver setting is outside its admitted domain.
    InvalidSetting {
        /// Stable field name.
        field: &'static str,
        /// Stable requirement.
        requirement: &'static str,
    },
    /// A solver-state vector has the wrong shape.
    VectorLength {
        /// `x` or `y`.
        vector: &'static str,
        /// Required length.
        expected: usize,
        /// Supplied length.
        actual: usize,
    },
    /// A solver-state entry is outside its numerical domain.
    InvalidVector {
        /// `x` or `y`.
        vector: &'static str,
        /// Offending entry.
        index: usize,
        /// Stable requirement.
        requirement: &'static str,
    },
}

impl core::fmt::Display for PdhgError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidSetting { field, requirement } => {
                write!(formatter, "PDHG setting {field} {requirement}")
            }
            Self::VectorLength {
                vector,
                expected,
                actual,
            } => write!(
                formatter,
                "PDHG state {vector} length {actual}; expected {expected}"
            ),
            Self::InvalidVector {
                vector,
                index,
                requirement,
            } => write!(formatter, "PDHG state {vector}[{index}] {requirement}"),
        }
    }
}

impl std::error::Error for PdhgError {}

/// Solve evidence.
#[derive(Debug, Clone, Default)]
pub struct PdhgReport {
    /// Iterations run.
    pub iters: usize,
    /// Final primal objective (volume).
    pub volume: f64,
    /// Final relative primal/dual objective separation diagnostic.
    pub gap: f64,
    /// Final equilibrium residual ‖Ax − b‖/‖b‖.
    pub eq_residual: f64,
    /// Gap trace (iteration, gap) at check intervals. Kept private so a
    /// solver-minted trace cannot be mutated after it is bound.
    trace: Vec<(usize, f64)>,
    /// Outward lower endpoint from an independently verified dual witness.
    verified_dual_lower_bound: Option<f64>,
    /// Outward upper endpoint from an exactly feasible repaired primal.
    verified_primal_upper_bound: Option<f64>,
    /// Uniform scale used by the retained verified dual witness.
    verified_dual_scale: Option<f64>,
    /// Content identity of the complete retained certificate.
    certificate_identity: Option<[u8; 32]>,
    /// Private binding of this report to one LP/settings/output state.
    solver_state_identity: Option<[u8; 32]>,
    /// Private snapshot preventing public diagnostic substitution.
    solver_scalar_snapshot: Option<[u8; 32]>,
    /// Private snapshot binding the complete solver-minted trace.
    solver_trace_snapshot: Option<[u8; 32]>,
    /// Solver-minted trace length, checked before any snapshot traversal.
    solver_trace_len: Option<usize>,
}

impl PdhgReport {
    /// Read-only solver convergence trace.
    #[must_use]
    pub fn trace(&self) -> &[(usize, f64)] {
        &self.trace
    }

    /// Verified dual lower bound, if a complete certificate was retained.
    #[must_use]
    pub fn verified_dual_lower_bound(&self) -> Option<f64> {
        self.certification_is_current()
            .then_some(self.verified_dual_lower_bound)
            .flatten()
    }

    /// Verified primal upper bound, if a complete certificate was retained.
    #[must_use]
    pub fn verified_primal_upper_bound(&self) -> Option<f64> {
        self.certification_is_current()
            .then_some(self.verified_primal_upper_bound)
            .flatten()
    }

    /// Uniform scale used by the retained verified dual witness.
    #[must_use]
    pub fn verified_dual_scale(&self) -> Option<f64> {
        self.certification_is_current()
            .then_some(self.verified_dual_scale)
            .flatten()
    }

    /// Content identity of the retained certificate.
    #[must_use]
    pub fn certificate_identity(&self) -> Option<[u8; 32]> {
        self.certification_is_current()
            .then_some(self.certificate_identity)
            .flatten()
    }

    fn scalar_snapshot(&self) -> [u8; 32] {
        let mut hasher = Blake3::new();
        hasher.update(&stable_usize_bytes(self.iters));
        hasher.update(&self.volume.to_bits().to_le_bytes());
        hasher.update(&self.gap.to_bits().to_le_bytes());
        hasher.update(&self.eq_residual.to_bits().to_le_bytes());
        let stream_hash = hasher.finalize();
        *hash_domain(PDHG_REPORT_SCALAR_SNAPSHOT_DOMAIN, stream_hash.as_bytes()).as_bytes()
    }

    fn trace_snapshot_checked<E>(
        &self,
        mut checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<[u8; 32], E> {
        let mut hasher = Blake3::new();
        hasher.update(&stable_usize_bytes(self.trace.len()));
        checkpoint()?;
        for (index, &(iteration, gap)) in self.trace.iter().enumerate() {
            hasher.update(&stable_usize_bytes(iteration));
            hasher.update(&gap.to_bits().to_le_bytes());
            if (index + 1) % LAYOUT_POLL_STRIDE == 0 {
                checkpoint()?;
            }
        }
        checkpoint()?;
        let stream_hash = hasher.finalize();
        Ok(*hash_domain(PDHG_REPORT_TRACE_SNAPSHOT_DOMAIN, stream_hash.as_bytes()).as_bytes())
    }

    fn trace_snapshot(&self) -> [u8; 32] {
        match self.trace_snapshot_checked(|| Ok::<(), Infallible>(())) {
            Ok(snapshot) => snapshot,
            Err(never) => match never {},
        }
    }

    fn certification_is_current(&self) -> bool {
        self.verified_dual_lower_bound.is_some()
            && self.verified_primal_upper_bound.is_some()
            && self.verified_dual_scale.is_some()
            && self.certificate_identity.is_some()
            && self.solver_state_identity.is_some()
            && self.solver_scalar_snapshot == Some(self.scalar_snapshot())
            && self.solver_trace_snapshot.is_some()
            && self.solver_trace_len == Some(self.trace.len())
            && self
                .trace
                .last()
                .is_some_and(|&(iteration, _)| iteration == self.iters)
    }

    pub(crate) fn bind_solver_state(&mut self, identity: [u8; 32]) {
        self.solver_state_identity = Some(identity);
        self.solver_scalar_snapshot = Some(self.scalar_snapshot());
        self.solver_trace_len = Some(self.trace.len());
        self.solver_trace_snapshot = Some(self.trace_snapshot());
    }

    pub(crate) fn matches_solver_state_checked<E>(
        &self,
        identity: [u8; 32],
        checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<bool, E> {
        if self.solver_state_identity != Some(identity)
            || self.solver_scalar_snapshot != Some(self.scalar_snapshot())
            || self.solver_trace_len != Some(self.trace.len())
        {
            return Ok(false);
        }
        Ok(
            self.solver_trace_snapshot == Some(self.trace_snapshot_checked(checkpoint)?)
                && self
                    .trace
                    .last()
                    .is_some_and(|&(iteration, _)| iteration == self.iters),
        )
    }

    pub(crate) fn clear_certified_bounds(&mut self) {
        self.verified_dual_lower_bound = None;
        self.verified_primal_upper_bound = None;
        self.verified_dual_scale = None;
        self.certificate_identity = None;
    }

    pub(crate) fn retain_certified_bounds(
        &mut self,
        lower: f64,
        upper: f64,
        dual_scale: f64,
        identity: [u8; 32],
    ) {
        self.verified_dual_lower_bound = Some(lower);
        self.verified_primal_upper_bound = Some(upper);
        self.verified_dual_scale = Some(dual_scale);
        self.certificate_identity = Some(identity);
    }

    /// Ledger row.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::new();
        let _ = write!(
            s,
            "{{\"iters\":{},\"volume\":{:.8e},\"gap\":{:.3e},\"eq_residual\":{:.3e}",
            self.iters, self.volume, self.gap, self.eq_residual
        );
        if self.certification_is_current()
            && let (Some(lower), Some(upper), Some(scale)) = (
                self.verified_dual_lower_bound,
                self.verified_primal_upper_bound,
                self.verified_dual_scale,
            )
        {
            let _ = write!(
                s,
                ",\"verified_dual_lower_bound\":{lower},\"verified_primal_upper_bound\":{upper},\"verified_dual_scale\":{scale}"
            );
            if let Some(identity) = self.certificate_identity {
                s.push_str(",\"certificate_identity\":\"");
                for byte in identity {
                    let _ = write!(s, "{byte:02x}");
                }
                s.push('"');
            }
        }
        s.push('}');
        s
    }
}

/// The assembled layout LP for one ground structure.
pub struct LayoutLp {
    /// Equilibrium matrix on free DOFs over split variables (n_free ×
    /// 2·members): columns `[q⁺ | q⁻]` with `B` and `−B` blocks.
    a: Csr,
    /// Aᵀ (materialized once; PDHG applies both directions).
    at: Csr,
    /// Cost per split variable (length/σ_y).
    c: Vec<f64>,
    /// Free-DOF load vector.
    b: Vec<f64>,
    /// Free-DOF index per (node, component); None = supported.
    dof_map: Vec<Option<usize>>,
    /// Estimated operator norm ‖A‖.
    norm_est: f64,
}

impl LayoutLp {
    /// Equilibrium matrix on the free degrees of freedom.
    #[must_use]
    pub fn a(&self) -> &Csr {
        &self.a
    }

    /// Materialized transpose of the equilibrium matrix.
    #[must_use]
    pub fn at(&self) -> &Csr {
        &self.at
    }

    /// Split-variable objective costs in `[q⁺ | q⁻]` order.
    #[must_use]
    pub fn c(&self) -> &[f64] {
        &self.c
    }

    /// Load vector after supported degrees of freedom are eliminated.
    #[must_use]
    pub fn b(&self) -> &[f64] {
        &self.b
    }

    /// Free-row identity for each node/component pair.
    #[must_use]
    pub fn dof_map(&self) -> &[Option<usize>] {
        &self.dof_map
    }

    /// Deterministic power-iteration estimate of the operator norm.
    #[must_use]
    pub fn norm_est(&self) -> f64 {
        self.norm_est
    }

    /// Assemble a bounded layout LP from admitted ground and load values.
    ///
    /// The constructor validates all dimensions and numerical domains, checks
    /// free-DOF, variable, triplet, and retained-byte budgets before allocating
    /// authoritative output, and polls cancellation at deterministic strides.
    ///
    /// # Errors
    /// Returns a structured refusal for malformed dimensions, non-finite
    /// numerical state, a load eliminated entirely by supports, a degenerate
    /// sparse operator, budget or allocation failure, or cancellation.
    #[allow(clippy::cast_precision_loss, clippy::too_many_lines)]
    pub fn try_assemble(
        gs: &GroundStructure,
        case: &LayoutCase,
        sigma_y: f64,
        limits: LayoutLimits,
        cx: &Cx<'_>,
    ) -> Result<LayoutLp, TrussConstructionError> {
        let nodes = gs.nodes();
        let members = gs.members();
        let lengths = gs.lengths();
        let n = nodes.len();
        let m = members.len();
        if n == 0 {
            return Err(TrussConstructionError::InvalidInput {
                field: "ground nodes",
                requirement: "must not be empty",
            });
        }
        if members.len() != lengths.len() {
            return Err(TrussConstructionError::VectorLength {
                field: "ground lengths",
                expected: members.len(),
                actual: lengths.len(),
            });
        }
        if case.node_count() != n {
            return Err(TrussConstructionError::VectorLength {
                field: "layout case nodes",
                expected: n,
                actual: case.node_count(),
            });
        }
        if !sigma_y.is_finite() || sigma_y <= 0.0 {
            return Err(TrussConstructionError::InvalidInput {
                field: "sigma_y",
                requirement: "must be finite and positive",
            });
        }
        let node_dofs = n.checked_mul(2).ok_or(TrussConstructionError::WorkBudget {
            resource: "layout node degrees of freedom",
            limit: limits.max_free_dofs(),
            observed: usize::MAX,
        })?;
        let variables = m.checked_mul(2).ok_or(TrussConstructionError::WorkBudget {
            resource: "layout split variables",
            limit: limits.max_variables(),
            observed: usize::MAX,
        })?;
        if variables == 0 {
            return Err(TrussConstructionError::InvalidInput {
                field: "ground members",
                requirement: "must not be empty",
            });
        }
        if variables > limits.max_variables() {
            return Err(TrussConstructionError::WorkBudget {
                resource: "layout split variables",
                limit: limits.max_variables(),
                observed: variables,
            });
        }

        let mut free_dofs = 0usize;
        let mut has_surviving_load = false;
        for (node, coordinates) in nodes.iter().enumerate() {
            if node % LAYOUT_POLL_STRIDE == 0 {
                poll(cx, "layout support/load admission")?;
            }
            if coordinates.iter().any(|coordinate| !coordinate.is_finite()) {
                return Err(TrussConstructionError::InvalidInput {
                    field: "ground nodes",
                    requirement: "must contain only finite coordinates",
                });
            }
            for comp in 0..2 {
                let load = case.loads[node][comp];
                if !load.is_finite() {
                    return Err(TrussConstructionError::InvalidInput {
                        field: "loads",
                        requirement: "must contain only finite components",
                    });
                }
                if !case.supported[node][comp] {
                    free_dofs =
                        free_dofs
                            .checked_add(1)
                            .ok_or(TrussConstructionError::WorkBudget {
                                resource: "layout free degrees of freedom",
                                limit: limits.max_free_dofs(),
                                observed: usize::MAX,
                            })?;
                    has_surviving_load |= load != 0.0;
                }
            }
        }
        if free_dofs == 0 {
            return Err(TrussConstructionError::InvalidInput {
                field: "supports",
                requirement: "must leave at least one free degree of freedom",
            });
        }
        if free_dofs > limits.max_free_dofs() {
            return Err(TrussConstructionError::WorkBudget {
                resource: "layout free degrees of freedom",
                limit: limits.max_free_dofs(),
                observed: free_dofs,
            });
        }
        if !has_surviving_load {
            return Err(TrussConstructionError::InvalidInput {
                field: "loads",
                requirement: "must contain a nonzero component on a free degree of freedom",
            });
        }

        let mut staged_triplets = 0usize;
        let mut previous_member = None;
        for (member_index, (&(a, b), &length)) in members.iter().zip(lengths).enumerate() {
            if member_index % LAYOUT_POLL_STRIDE == 0 {
                poll(cx, "layout member admission")?;
            }
            if a >= n || b >= n || a >= b {
                return Err(TrussConstructionError::InvalidInput {
                    field: "ground members",
                    requirement: "must contain distinct in-range canonical endpoints a < b",
                });
            }
            if previous_member.is_some_and(|previous| previous >= (a, b)) {
                return Err(TrussConstructionError::InvalidInput {
                    field: "ground members",
                    requirement: "must be unique and strictly lexicographically ordered",
                });
            }
            previous_member = Some((a, b));
            if !length.is_finite() || length <= 0.0 {
                return Err(TrussConstructionError::InvalidInput {
                    field: "ground lengths",
                    requirement: "must contain only finite positive lengths",
                });
            }
            let dx = (nodes[b][0] - nodes[a][0]) / length;
            let dy = (nodes[b][1] - nodes[a][1]) / length;
            if !dx.is_finite() || !dy.is_finite() {
                return Err(TrussConstructionError::InvalidInput {
                    field: "ground directions",
                    requirement: "must remain finite after length normalization",
                });
            }
            for (node, comp) in [(a, 0usize), (a, 1), (b, 0), (b, 1)] {
                if !case.supported[node][comp] {
                    staged_triplets = staged_triplets.checked_add(2).ok_or(
                        TrussConstructionError::WorkBudget {
                            resource: "layout staged triplets",
                            limit: limits.max_staged_triplets(),
                            observed: usize::MAX,
                        },
                    )?;
                }
            }
        }
        if staged_triplets == 0 {
            return Err(TrussConstructionError::InvalidInput {
                field: "equilibrium matrix",
                requirement: "must contain at least one staged contribution",
            });
        }
        if staged_triplets > limits.max_staged_triplets() {
            return Err(TrussConstructionError::WorkBudget {
                resource: "layout staged triplets",
                limit: limits.max_staged_triplets(),
                observed: staged_triplets,
            });
        }
        let retained_bytes =
            retained_layout_bytes(node_dofs, free_dofs, variables, staged_triplets).ok_or(
                TrussConstructionError::WorkBudget {
                    resource: "layout retained bytes",
                    limit: limits.max_retained_bytes(),
                    observed: usize::MAX,
                },
            )?;
        if retained_bytes > limits.max_retained_bytes() {
            return Err(TrussConstructionError::WorkBudget {
                resource: "layout retained bytes",
                limit: limits.max_retained_bytes(),
                observed: retained_bytes,
            });
        }

        let mut dof_map = empty_with_capacity(node_dofs, "layout DOF map")?;
        let mut next_free = 0usize;
        for node in 0..n {
            if node % LAYOUT_POLL_STRIDE == 0 {
                poll(cx, "layout DOF map construction")?;
            }
            for comp in 0..2 {
                if case.supported[node][comp] {
                    dof_map.push(None);
                } else {
                    dof_map.push(Some(next_free));
                    next_free += 1;
                }
            }
        }
        debug_assert_eq!(next_free, free_dofs);

        let a_mat = assemble_equilibrium(
            nodes,
            members,
            lengths,
            &dof_map,
            free_dofs,
            variables,
            staged_triplets,
            cx,
        )?;
        if a_mat.nrows() != free_dofs
            || a_mat.ncols() != variables
            || a_mat.nnz() != staged_triplets
            || !csr_is_finite_and_nonzero(&a_mat, cx, "layout sparse validation")?
        {
            return Err(TrussConstructionError::InvalidInput {
                field: "equilibrium matrix",
                requirement: "must have exact dimensions/nnz and finite nonzero state",
            });
        }
        let at = transpose_checked(&a_mat, cx)?;
        if at.nrows() != variables
            || at.ncols() != free_dofs
            || at.nnz() != a_mat.nnz()
            || !csr_is_finite_and_nonzero(&at, cx, "layout transpose validation")?
        {
            return Err(TrussConstructionError::InvalidInput {
                field: "equilibrium transpose",
                requirement: "must preserve finite canonical matrix state",
            });
        }

        let mut b_vec = zeroed_f64(
            free_dofs,
            "layout load vector",
            cx,
            "layout load-vector allocation",
        )?;
        for node in 0..n {
            if node % LAYOUT_POLL_STRIDE == 0 {
                poll(cx, "layout load assembly")?;
            }
            for comp in 0..2 {
                if let Some(row) = dof_map[2 * node + comp] {
                    b_vec[row] = case.loads[node][comp];
                }
            }
        }
        let mut squared_load_norm = 0.0f64;
        let mut load_entries_visited = 0usize;
        for (row, &load) in b_vec.iter().enumerate() {
            if row % LAYOUT_POLL_STRIDE == 0 {
                poll(cx, "layout load validation")?;
            }
            let squared = load * load;
            squared_load_norm += squared;
            if !load.is_finite() || !squared.is_finite() || !squared_load_norm.is_finite() {
                return Err(TrussConstructionError::InvalidInput {
                    field: "free load vector",
                    requirement: "must have a finite squared Euclidean norm",
                });
            }
            if load != 0.0 {
                let mut connected = false;
                for &coefficient in a_mat.row(row).1 {
                    if load_entries_visited.is_multiple_of(LAYOUT_POLL_STRIDE) {
                        poll(cx, "layout loaded-row validation")?;
                    }
                    load_entries_visited += 1;
                    connected |= coefficient != 0.0;
                }
                if !connected {
                    return Err(TrussConstructionError::InvalidInput {
                        field: "free load vector",
                        requirement: "must act only on free degrees of freedom connected to a member",
                    });
                }
            }
        }
        if squared_load_norm < MIN_LAYOUT_LOAD_NORM_SQUARED {
            return Err(TrussConstructionError::InvalidInput {
                field: "free load vector",
                requirement: "must have squared Euclidean norm at least 1e-60",
            });
        }

        let mut c = empty_with_capacity(variables, "layout objective costs")?;
        for (index, &length) in lengths.iter().enumerate() {
            if index % LAYOUT_POLL_STRIDE == 0 {
                poll(cx, "layout cost assembly")?;
            }
            let cost = length / sigma_y;
            if !cost.is_finite() || cost <= 0.0 {
                return Err(TrussConstructionError::InvalidInput {
                    field: "layout objective costs",
                    requirement: "must remain finite and positive",
                });
            }
            c.push(cost);
        }
        for (index, &length) in lengths.iter().enumerate() {
            if index % LAYOUT_POLL_STRIDE == 0 {
                poll(cx, "layout cost assembly")?;
            }
            c.push(length / sigma_y);
        }
        // Power iteration for ‖A‖ (deterministic start).
        let mut v = empty_with_capacity(variables, "layout norm input")?;
        for i in 0..variables {
            if i % LAYOUT_POLL_STRIDE == 0 {
                poll(cx, "layout norm-input allocation")?;
            }
            v.push(1.0 + (i % 7) as f64 * 0.1);
        }
        let mut norm_est = 1.0;
        let mut av = zeroed_f64(
            free_dofs,
            "layout norm row workspace",
            cx,
            "layout norm-row allocation",
        )?;
        let mut atv = zeroed_f64(
            variables,
            "layout norm column workspace",
            cx,
            "layout norm-column allocation",
        )?;
        let mut completed_iterations = 0usize;
        let mut fallback_used = false;
        while completed_iterations < 30 {
            poll(cx, "layout norm estimation")?;
            spmv_checked(&a_mat, &v, &mut av, cx, "layout norm forward multiply")?;
            spmv_checked(&at, &av, &mut atv, cx, "layout norm transpose multiply")?;
            for (index, value) in av.iter().enumerate() {
                if index % LAYOUT_POLL_STRIDE == 0 {
                    poll(cx, "layout norm-row validation")?;
                }
                if !value.is_finite() {
                    return Err(TrussConstructionError::InvalidInput {
                        field: "layout norm estimate",
                        requirement: "must retain finite power-iteration state",
                    });
                }
            }
            let mut squared_norm = 0.0f64;
            for (index, value) in atv.iter().enumerate() {
                if index % LAYOUT_POLL_STRIDE == 0 {
                    poll(cx, "layout norm-column validation")?;
                }
                if !value.is_finite() {
                    return Err(TrussConstructionError::InvalidInput {
                        field: "layout norm estimate",
                        requirement: "must retain finite power-iteration state",
                    });
                }
                squared_norm += value * value;
            }
            if squared_norm == 0.0 && !fallback_used {
                let fallback_column = first_nonzero_column(&a_mat, cx)?;
                for (index, value) in v.iter_mut().enumerate() {
                    if index % LAYOUT_POLL_STRIDE == 0 {
                        poll(cx, "layout norm fallback seed")?;
                    }
                    *value = 0.0;
                }
                v[fallback_column] = 1.0;
                fallback_used = true;
                continue;
            }
            if !squared_norm.is_finite() || squared_norm <= 0.0 {
                return Err(TrussConstructionError::InvalidInput {
                    field: "layout norm estimate",
                    requirement: "must have finite positive iteration norm",
                });
            }
            let nrm = squared_norm.sqrt();
            norm_est = nrm.sqrt();
            if !norm_est.is_finite() || norm_est <= 0.0 {
                return Err(TrussConstructionError::InvalidInput {
                    field: "layout norm estimate",
                    requirement: "must be finite and positive",
                });
            }
            for (index, (vi, ai)) in v.iter_mut().zip(&atv).enumerate() {
                if index % LAYOUT_POLL_STRIDE == 0 {
                    poll(cx, "layout norm-vector update")?;
                }
                *vi = ai / nrm;
            }
            completed_iterations += 1;
        }
        poll(cx, "layout publication")?;
        Ok(LayoutLp {
            a: a_mat,
            at,
            c,
            b: b_vec,
            dof_map,
            norm_est,
        })
    }

    /// Run PDHG from a warm start (zeros for cold); returns the
    /// primal solution (split forces) and the report.
    ///
    /// # Errors
    /// Refuses zero iteration/check intervals, non-finite or out-of-range
    /// tolerances, malformed warm-start lengths, non-finite state, and negative
    /// primal warm starts before entering the iteration loop.
    #[allow(clippy::too_many_lines)] // validation plus one diagnostic iteration loop
    pub fn solve(
        &self,
        warm_x: Option<Vec<f64>>,
        warm_y: Option<Vec<f64>>,
        settings: PdhgSettings,
    ) -> Result<(Vec<f64>, Vec<f64>, PdhgReport), PdhgError> {
        let nvar = self.c.len();
        let nrow = self.b.len();
        if settings.max_iters == 0 {
            return Err(PdhgError::InvalidSetting {
                field: "max_iters",
                requirement: "must be at least one",
            });
        }
        if settings.max_iters > MAX_PDHG_ITERS {
            return Err(PdhgError::InvalidSetting {
                field: "max_iters",
                requirement: "exceeds the one-million-iteration direct-solve limit",
            });
        }
        if settings.check_every == 0 {
            return Err(PdhgError::InvalidSetting {
                field: "check_every",
                requirement: "must be at least one",
            });
        }
        if !settings.gap_tol.is_finite() || !(0.0..=1.0).contains(&settings.gap_tol) {
            return Err(PdhgError::InvalidSetting {
                field: "gap_tol",
                requirement: "must be finite and in 0..=1",
            });
        }
        if let Some(values) = &warm_x {
            if values.len() != nvar {
                return Err(PdhgError::VectorLength {
                    vector: "x",
                    expected: nvar,
                    actual: values.len(),
                });
            }
            if let Some(index) = values
                .iter()
                .position(|value| !value.is_finite() || *value < 0.0)
            {
                return Err(PdhgError::InvalidVector {
                    vector: "x",
                    index,
                    requirement: "must be finite and non-negative",
                });
            }
        }
        if let Some(values) = &warm_y {
            if values.len() != nrow {
                return Err(PdhgError::VectorLength {
                    vector: "y",
                    expected: nrow,
                    actual: values.len(),
                });
            }
            if let Some(index) = values.iter().position(|value| !value.is_finite()) {
                return Err(PdhgError::InvalidVector {
                    vector: "y",
                    index,
                    requirement: "must be finite",
                });
            }
        }
        let mut x = warm_x.unwrap_or_else(|| vec![0.0; nvar]);
        let mut y = warm_y.unwrap_or_else(|| vec![0.0; nrow]);
        let step = 0.95 / self.norm_est.max(1e-30);
        let (tau, sigma) = (step, step);
        let bnorm = self.b.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-30);
        let mut report = PdhgReport::default();
        let mut aty = vec![0.0f64; nvar];
        let mut ax = vec![0.0f64; nrow];
        let mut x_prev = x.clone();
        let mut xbar = vec![0.0f64; nvar];
        for it in 0..settings.max_iters {
            // x ← Π₊(x − τ(c + Aᵀy))
            self.at.spmv(&y, &mut aty);
            x_prev.copy_from_slice(&x);
            for i in 0..nvar {
                x[i] = (x[i] - tau * (self.c[i] + aty[i])).max(0.0);
            }
            // y ← y + σ(A(2x − x_prev) − b)
            for ((extrapolated, xi), previous) in xbar.iter_mut().zip(&x).zip(&x_prev) {
                *extrapolated = 2.0 * xi - previous;
            }
            self.a.spmv(&xbar, &mut ax);
            for r in 0..nrow {
                y[r] += sigma * (ax[r] - self.b[r]);
            }
            if (it + 1) % settings.check_every == 0 || it + 1 == settings.max_iters {
                let (gap, eq_res, primal) = self.diagnostics(&x, &y, bnorm)?;
                report.trace.push((it + 1, gap));
                report.iters = it + 1;
                report.volume = primal;
                report.gap = gap;
                report.eq_residual = eq_res;
                if gap < settings.gap_tol && eq_res < settings.gap_tol {
                    break;
                }
            }
        }
        report.bind_solver_state(crate::certificate::solver_state_identity(
            self, &x, &y, settings,
        ));
        Ok((x, y, report))
    }

    /// Return `(relative objective separation, equilibrium residual, primal
    /// objective)`. With the saddle `cᵀx + yᵀ(Ax − b)`, the nominal dual
    /// objective is `−bᵀy` under `c + Aᵀy ≥ 0`; scaling `y` repairs observed
    /// floating violations. This routine does not outward-verify dual
    /// feasibility or repair the primal to exact equilibrium, so its tuple is
    /// diagnostic rather than a finite optimum certificate.
    ///
    /// # Errors
    /// Refuses dimension mismatch, non-finite state, negative primal entries,
    /// or a non-finite/non-positive load norm before sparse operations.
    pub fn diagnostics(
        &self,
        x: &[f64],
        y: &[f64],
        bnorm: f64,
    ) -> Result<(f64, f64, f64), PdhgError> {
        if x.len() != self.c.len() {
            return Err(PdhgError::VectorLength {
                vector: "x",
                expected: self.c.len(),
                actual: x.len(),
            });
        }
        if y.len() != self.b.len() {
            return Err(PdhgError::VectorLength {
                vector: "y",
                expected: self.b.len(),
                actual: y.len(),
            });
        }
        if let Some(index) = x
            .iter()
            .position(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(PdhgError::InvalidVector {
                vector: "x",
                index,
                requirement: "must be finite and non-negative",
            });
        }
        if let Some(index) = y.iter().position(|value| !value.is_finite()) {
            return Err(PdhgError::InvalidVector {
                vector: "y",
                index,
                requirement: "must be finite",
            });
        }
        if !bnorm.is_finite() || bnorm <= 0.0 {
            return Err(PdhgError::InvalidSetting {
                field: "bnorm",
                requirement: "must be finite and positive",
            });
        }
        let primal: f64 = self.c.iter().zip(x).map(|(c, x)| c * x).sum();
        let mut aty = vec![0.0f64; self.c.len()];
        self.at.spmv(y, &mut aty);
        let mut scale = 1.0f64;
        for (a, c) in aty.iter().zip(&self.c) {
            // Violation where c + Aᵀy < 0, i.e. aty < −c.
            if *a < -c && *a < 0.0 {
                scale = scale.min(-c / a);
            }
        }
        let dual: f64 = -(y.iter().zip(&self.b).map(|(y, b)| y * b).sum::<f64>()) * scale.max(0.0);
        let mut ax = vec![0.0f64; self.b.len()];
        self.a.spmv(x, &mut ax);
        let eq_res = ax
            .iter()
            .zip(&self.b)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt()
            / bnorm;
        let gap = (primal - dual).abs() / primal.abs().max(1e-30);
        Ok((gap, eq_res, primal))
    }
}

#[cfg(test)]
mod tests {
    use super::{LAYOUT_POLL_STRIDE, PdhgReport};

    #[test]
    fn checked_report_snapshot_covers_non_tail_trace_entries() {
        let mut report = PdhgReport {
            iters: 2,
            volume: 1.0,
            gap: 0.25,
            eq_residual: 0.125,
            trace: vec![(1, 0.5), (2, 0.25)],
            ..PdhgReport::default()
        };
        report.bind_solver_state([7; 32]);
        report.retain_certified_bounds(0.75, 1.25, 0.5, [9; 32]);
        assert_eq!(report.verified_dual_lower_bound(), Some(0.75));
        assert_eq!(report.trace(), &[(1, 0.5), (2, 0.25)]);
        assert!(report.to_json().contains("\"certificate_identity\""));
        assert_eq!(
            report.matches_solver_state_checked([7; 32], || Ok::<(), ()>(())),
            Ok(true)
        );

        report.trace[0].1 = f64::from_bits(report.trace[0].1.to_bits().wrapping_add(1));
        assert_eq!(
            report.matches_solver_state_checked([7; 32], || Ok::<(), ()>(())),
            Ok(false)
        );
    }

    #[test]
    fn checked_report_snapshot_polls_at_bounded_strides() {
        let trace_len = LAYOUT_POLL_STRIDE * 2 + 1;
        let mut report = PdhgReport {
            iters: trace_len,
            volume: 1.0,
            gap: 0.25,
            eq_residual: 0.125,
            trace: (1..=trace_len).map(|iteration| (iteration, 0.25)).collect(),
            ..PdhgReport::default()
        };
        report.bind_solver_state([7; 32]);

        let mut interrupted_polls = 0usize;
        let interrupted = report.matches_solver_state_checked([7; 32], || {
            interrupted_polls += 1;
            if interrupted_polls == 2 {
                Err(())
            } else {
                Ok(())
            }
        });
        assert_eq!(interrupted, Err(()));
        assert_eq!(interrupted_polls, 2);

        let mut completed_polls = 0usize;
        assert_eq!(
            report.matches_solver_state_checked([7; 32], || {
                completed_polls += 1;
                Ok::<(), ()>(())
            }),
            Ok(true)
        );
        assert_eq!(completed_polls, 4);
    }

    #[test]
    fn checked_report_snapshot_rejects_length_mismatch_before_polling() {
        let mut report = PdhgReport {
            iters: 2,
            volume: 1.0,
            gap: 0.25,
            eq_residual: 0.125,
            trace: vec![(1, 0.5), (2, 0.25)],
            ..PdhgReport::default()
        };
        report.bind_solver_state([7; 32]);
        report.trace.insert(0, (0, 0.0));

        let mut polls = 0usize;
        assert_eq!(
            report.matches_solver_state_checked([7; 32], || {
                polls += 1;
                Ok::<(), ()>(())
            }),
            Ok(false)
        );
        assert_eq!(polls, 0);
    }
}
