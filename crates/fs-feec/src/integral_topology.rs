//! Bounded exact-integer topology algebra for I13.2b.
//!
//! This first tranche is an independently replayable Smith-normal-form
//! *witness verifier*. It accepts a general exact integer matrix only when
//! explicit integer left/right transformations and their explicit inverses
//! prove
//!
//! `U * A * V = D`, `U^-1 * U = U * U^-1 = I`, and
//! `V^-1 * V = V * V^-1 = I`.
//!
//! The diagonal must be canonical: nonnegative invariant factors precede
//! zeros and each nonzero factor divides its successor. All arithmetic is
//! checked `i128`; overflow, allocation pressure, work exhaustion, and
//! cancellation refuse without publishing a partially verified value.
//!
//! A verified value is tagged [`TopologyApplicability::AbstractAlgebraOnly`].
//! It is not yet a terminal-relative homology receipt or physical R3 winding
//! authority. The constructive normal-form solver and the binding to an
//! admitted [`crate::terminal_relative::TerminalRelativePair`] follow in
//! later I13.2b tranches.

use core::fmt;

/// Default maximum row or column extent admitted by the exact checker.
pub const DEFAULT_MAX_MATRIX_EXTENT: usize = 256;
/// Default maximum entries in any one admitted matrix.
pub const DEFAULT_MAX_MATRIX_ENTRIES: usize = DEFAULT_MAX_MATRIX_EXTENT * DEFAULT_MAX_MATRIX_EXTENT;
/// Default maximum entries retained across the source and five witness
/// matrices.
pub const DEFAULT_MAX_RETAINED_ENTRIES: usize = 6 * DEFAULT_MAX_MATRIX_ENTRIES;
/// Default maximum scratch entries used by exact witness multiplication.
pub const DEFAULT_MAX_WORKSPACE_ENTRIES: usize = DEFAULT_MAX_MATRIX_ENTRIES;
/// Exact dot-product terms (one checked multiply/add pair each) admitted by
/// the default checker.
pub const DEFAULT_MAX_SCALAR_OPERATIONS: u128 = 101_000_000;

/// Explicit resource envelope for exact integer witness admission.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExactAlgebraBudget {
    max_rows: usize,
    max_cols: usize,
    max_matrix_entries: usize,
    max_retained_entries: usize,
    max_workspace_entries: usize,
    max_scalar_operations: u128,
}

impl ExactAlgebraBudget {
    /// Construct an exact algebra envelope. Zero limits are valid and admit
    /// only the corresponding empty structures.
    #[must_use]
    pub const fn new(
        max_rows: usize,
        max_cols: usize,
        max_matrix_entries: usize,
        max_retained_entries: usize,
        max_workspace_entries: usize,
        max_scalar_operations: u128,
    ) -> Self {
        Self {
            max_rows,
            max_cols,
            max_matrix_entries,
            max_retained_entries,
            max_workspace_entries,
            max_scalar_operations,
        }
    }

    /// Maximum admitted row count.
    #[must_use]
    pub const fn max_rows(self) -> usize {
        self.max_rows
    }

    /// Maximum admitted column count.
    #[must_use]
    pub const fn max_cols(self) -> usize {
        self.max_cols
    }

    /// Maximum entries in one matrix.
    #[must_use]
    pub const fn max_matrix_entries(self) -> usize {
        self.max_matrix_entries
    }

    /// Maximum entries retained by the source plus complete witness.
    #[must_use]
    pub const fn max_retained_entries(self) -> usize {
        self.max_retained_entries
    }

    /// Maximum internal scratch entries.
    #[must_use]
    pub const fn max_workspace_entries(self) -> usize {
        self.max_workspace_entries
    }

    /// Maximum exact dot-product terms (one checked multiply/add pair each).
    #[must_use]
    pub const fn max_scalar_operations(self) -> u128 {
        self.max_scalar_operations
    }
}

impl Default for ExactAlgebraBudget {
    fn default() -> Self {
        Self::new(
            DEFAULT_MAX_MATRIX_EXTENT,
            DEFAULT_MAX_MATRIX_EXTENT,
            DEFAULT_MAX_MATRIX_ENTRIES,
            DEFAULT_MAX_RETAINED_ENTRIES,
            DEFAULT_MAX_WORKSPACE_ENTRIES,
            DEFAULT_MAX_SCALAR_OPERATIONS,
        )
    }
}

/// Dense row-major exact integer matrix with admitted extents.
///
/// Dense storage is intentional in this first witness checker: every retained
/// entry is hard-capped before any checker allocation, and exact product
/// verification visits a deterministic rectangular domain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExactIntegerMatrix {
    rows: usize,
    cols: usize,
    entries: Vec<i128>,
}

impl ExactIntegerMatrix {
    /// Admit a row-major matrix without normalizing or narrowing any integer.
    pub fn try_new(
        rows: usize,
        cols: usize,
        entries: Vec<i128>,
        budget: ExactAlgebraBudget,
    ) -> Result<Self, IntegralTopologyError> {
        if rows > budget.max_rows || cols > budget.max_cols {
            return Err(IntegralTopologyError::MatrixExtentExceeded {
                rows,
                cols,
                max_rows: budget.max_rows,
                max_cols: budget.max_cols,
            });
        }
        let expected = rows
            .checked_mul(cols)
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "matrix entry count",
            })?;
        if expected > budget.max_matrix_entries {
            return Err(IntegralTopologyError::MatrixEntryBudgetExceeded {
                requested: expected,
                max: budget.max_matrix_entries,
            });
        }
        if entries.len() != expected {
            return Err(IntegralTopologyError::MatrixEntryCount {
                rows,
                cols,
                expected,
                actual: entries.len(),
            });
        }
        Ok(Self {
            rows,
            cols,
            entries,
        })
    }

    /// Row count.
    #[must_use]
    pub const fn rows(&self) -> usize {
        self.rows
    }

    /// Column count.
    #[must_use]
    pub const fn cols(&self) -> usize {
        self.cols
    }

    /// Canonical row-major entries.
    #[must_use]
    pub fn entries(&self) -> &[i128] {
        &self.entries
    }

    /// Exact entry at `(row, col)`, or `None` outside the admitted rectangle.
    #[must_use]
    pub fn get(&self, row: usize, col: usize) -> Option<i128> {
        if row >= self.rows || col >= self.cols {
            return None;
        }
        Some(self.entries[row * self.cols + col])
    }

    fn entry(&self, row: usize, col: usize) -> i128 {
        self.entries[row * self.cols + col]
    }

    fn ensure_within(
        &self,
        role: MatrixRole,
        budget: ExactAlgebraBudget,
    ) -> Result<(), IntegralTopologyError> {
        if self.rows > budget.max_rows
            || self.cols > budget.max_cols
            || self.entries.len() > budget.max_matrix_entries
        {
            return Err(IntegralTopologyError::RetainedMatrixExceedsBudget {
                role,
                rows: self.rows,
                cols: self.cols,
                entries: self.entries.len(),
            });
        }
        Ok(())
    }
}

/// Role of one retained matrix in a Smith witness.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatrixRole {
    /// Original exact matrix.
    Source,
    /// Claimed canonical diagonal matrix.
    Diagonal,
    /// Left transformation `U`.
    LeftTransform,
    /// Explicit inverse `U^-1`.
    LeftInverse,
    /// Right transformation `V`.
    RightTransform,
    /// Explicit inverse `V^-1`.
    RightInverse,
}

/// Exact product or inverse identity being checked.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SmithWitnessStage {
    /// `U * U^-1 = I`.
    LeftTimesInverse,
    /// `U^-1 * U = I`.
    LeftInverseTimesTransform,
    /// `V * V^-1 = I`.
    RightTimesInverse,
    /// `V^-1 * V = I`.
    RightInverseTimesTransform,
    /// Intermediate `U * A`.
    LeftTimesSource,
    /// Final `U * A * V = D`.
    DiagonalTransform,
}

/// Untrusted complete Smith-normal-form witness.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SmithNormalFormWitness {
    diagonal: ExactIntegerMatrix,
    left: ExactIntegerMatrix,
    left_inverse: ExactIntegerMatrix,
    right: ExactIntegerMatrix,
    right_inverse: ExactIntegerMatrix,
}

impl SmithNormalFormWitness {
    /// Assemble an untrusted witness. Mathematical admission occurs only in
    /// [`verify_smith_normal_form`].
    #[must_use]
    pub const fn new(
        diagonal: ExactIntegerMatrix,
        left: ExactIntegerMatrix,
        left_inverse: ExactIntegerMatrix,
        right: ExactIntegerMatrix,
        right_inverse: ExactIntegerMatrix,
    ) -> Self {
        Self {
            diagonal,
            left,
            left_inverse,
            right,
            right_inverse,
        }
    }

    /// Claimed canonical diagonal.
    #[must_use]
    pub const fn diagonal(&self) -> &ExactIntegerMatrix {
        &self.diagonal
    }

    /// Claimed left transformation.
    #[must_use]
    pub const fn left(&self) -> &ExactIntegerMatrix {
        &self.left
    }

    /// Claimed left inverse.
    #[must_use]
    pub const fn left_inverse(&self) -> &ExactIntegerMatrix {
        &self.left_inverse
    }

    /// Claimed right transformation.
    #[must_use]
    pub const fn right(&self) -> &ExactIntegerMatrix {
        &self.right
    }

    /// Claimed right inverse.
    #[must_use]
    pub const fn right_inverse(&self) -> &ExactIntegerMatrix {
        &self.right_inverse
    }
}

/// Scope intentionally carried by the first exact-algebra tranche.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TopologyApplicability {
    /// Algebraic verification only. No physical R3 embedding or winding
    /// conclusion may consume this value as authority.
    AbstractAlgebraOnly,
}

/// Authority classification for an unsuccessful exact verification.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntegralTopologyFailureClass {
    /// Exact supplied structure or witness bytes contradict the claimed Smith
    /// decomposition.
    Refuted,
    /// Resource, cancellation, allocation, or arithmetic limits prevented a
    /// mathematical decision.
    Unknown,
}

/// Opaque successfully verified Smith normal form.
///
/// There is no public constructor. The exact source and every witness matrix
/// remain attached so the authority cannot be replayed against another input.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VerifiedSmithNormalForm {
    source: ExactIntegerMatrix,
    witness: SmithNormalFormWitness,
    invariant_factors: Vec<i128>,
    rank: usize,
    scalar_operations: u128,
}

impl VerifiedSmithNormalForm {
    /// Exact source matrix bound to this verification.
    #[must_use]
    pub const fn source(&self) -> &ExactIntegerMatrix {
        &self.source
    }

    /// Canonical diagonal matrix.
    #[must_use]
    pub const fn diagonal(&self) -> &ExactIntegerMatrix {
        &self.witness.diagonal
    }

    /// Verified left transformation.
    #[must_use]
    pub const fn left_transform(&self) -> &ExactIntegerMatrix {
        &self.witness.left
    }

    /// Verified left inverse.
    #[must_use]
    pub const fn left_inverse(&self) -> &ExactIntegerMatrix {
        &self.witness.left_inverse
    }

    /// Verified right transformation.
    #[must_use]
    pub const fn right_transform(&self) -> &ExactIntegerMatrix {
        &self.witness.right
    }

    /// Verified right inverse.
    #[must_use]
    pub const fn right_inverse(&self) -> &ExactIntegerMatrix {
        &self.witness.right_inverse
    }

    /// Positive canonical invariant factors.
    #[must_use]
    pub fn invariant_factors(&self) -> &[i128] {
        &self.invariant_factors
    }

    /// Rank over the integers/rationals, equal to the number of positive
    /// invariant factors.
    #[must_use]
    pub const fn rank(&self) -> usize {
        self.rank
    }

    /// Exact dot-product terms completed by verification. Each term performs
    /// one checked multiplication followed by one checked addition.
    #[must_use]
    pub const fn scalar_operations(&self) -> u128 {
        self.scalar_operations
    }

    /// This first tranche is deliberately not physical topology authority.
    #[must_use]
    pub const fn applicability(&self) -> TopologyApplicability {
        TopologyApplicability::AbstractAlgebraOnly
    }
}

/// Verify a complete Smith witness without injected cancellation.
pub fn verify_smith_normal_form(
    source: ExactIntegerMatrix,
    witness: SmithNormalFormWitness,
    budget: ExactAlgebraBudget,
) -> Result<VerifiedSmithNormalForm, IntegralTopologyError> {
    verify_smith_normal_form_with_checkpoint(source, witness, budget, &mut |_| true)
}

/// Verify a complete Smith witness with bounded cancellation polling.
///
/// After bounded structural preflight, the callback runs before exact work,
/// before each output scalar, and once more before final publication. A
/// callback returning `false` publishes only [`IntegralTopologyError::Cancelled`].
/// Between polls at most `max(rows, cols)` checked dot-product terms execute.
pub fn verify_smith_normal_form_with_checkpoint(
    source: ExactIntegerMatrix,
    witness: SmithNormalFormWitness,
    budget: ExactAlgebraBudget,
    checkpoint: &mut impl FnMut(&'static str) -> bool,
) -> Result<VerifiedSmithNormalForm, IntegralTopologyError> {
    preflight_shapes(&source, &witness, budget)?;
    let rows = source.rows;
    let cols = source.cols;
    let planned = planned_scalar_operations(rows, cols)?;
    if planned > budget.max_scalar_operations {
        return Err(IntegralTopologyError::ScalarWorkBudgetExceeded {
            requested: planned,
            max: budget.max_scalar_operations,
        });
    }
    poll(checkpoint, "smith witness preflight", 0, planned)?;

    let mut completed = 0_u128;
    verify_identity_product(
        &witness.left,
        &witness.left_inverse,
        SmithWitnessStage::LeftTimesInverse,
        checkpoint,
        &mut completed,
        planned,
    )?;
    verify_identity_product(
        &witness.left_inverse,
        &witness.left,
        SmithWitnessStage::LeftInverseTimesTransform,
        checkpoint,
        &mut completed,
        planned,
    )?;
    verify_identity_product(
        &witness.right,
        &witness.right_inverse,
        SmithWitnessStage::RightTimesInverse,
        checkpoint,
        &mut completed,
        planned,
    )?;
    verify_identity_product(
        &witness.right_inverse,
        &witness.right,
        SmithWitnessStage::RightInverseTimesTransform,
        checkpoint,
        &mut completed,
        planned,
    )?;

    let invariant_factors =
        verify_canonical_diagonal(&witness.diagonal, checkpoint, completed, planned)?;
    verify_transform_product(&source, &witness, checkpoint, &mut completed, planned)?;
    poll(checkpoint, "smith witness finalize", completed, planned)?;
    debug_assert_eq!(completed, planned);

    Ok(VerifiedSmithNormalForm {
        source,
        witness,
        rank: invariant_factors.len(),
        invariant_factors,
        scalar_operations: completed,
    })
}

fn verify_transform_product(
    source: &ExactIntegerMatrix,
    witness: &SmithNormalFormWitness,
    checkpoint: &mut impl FnMut(&'static str) -> bool,
    completed: &mut u128,
    planned: u128,
) -> Result<(), IntegralTopologyError> {
    let rows = source.rows;
    let cols = source.cols;
    let workspace = rows
        .checked_mul(cols)
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "smith workspace entries",
        })?;
    let mut left_times_source = allocate_zeroed(workspace, "left-times-source workspace")?;
    for row in 0..rows {
        for col in 0..cols {
            poll(checkpoint, "smith left-times-source", *completed, planned)?;
            left_times_source[row * cols + col] = checked_dot(
                witness.left(),
                row,
                source,
                col,
                SmithWitnessStage::LeftTimesSource,
                completed,
            )?;
        }
    }
    let left_times_source = ExactIntegerMatrix {
        rows,
        cols,
        entries: left_times_source,
    };
    for row in 0..rows {
        for col in 0..cols {
            poll(checkpoint, "smith diagonal transform", *completed, planned)?;
            let actual = checked_dot(
                &left_times_source,
                row,
                witness.right(),
                col,
                SmithWitnessStage::DiagonalTransform,
                completed,
            )?;
            let expected = witness.diagonal().entry(row, col);
            if actual != expected {
                return Err(IntegralTopologyError::WitnessProductMismatch {
                    stage: SmithWitnessStage::DiagonalTransform,
                    row,
                    col,
                    expected,
                    actual,
                });
            }
        }
    }
    Ok(())
}

fn preflight_shapes(
    source: &ExactIntegerMatrix,
    witness: &SmithNormalFormWitness,
    budget: ExactAlgebraBudget,
) -> Result<(), IntegralTopologyError> {
    for (role, matrix) in [
        (MatrixRole::Source, source),
        (MatrixRole::Diagonal, &witness.diagonal),
        (MatrixRole::LeftTransform, &witness.left),
        (MatrixRole::LeftInverse, &witness.left_inverse),
        (MatrixRole::RightTransform, &witness.right),
        (MatrixRole::RightInverse, &witness.right_inverse),
    ] {
        matrix.ensure_within(role, budget)?;
    }

    let rows = source.rows;
    let cols = source.cols;
    require_shape(&witness.diagonal, MatrixRole::Diagonal, rows, cols)?;
    require_shape(&witness.left, MatrixRole::LeftTransform, rows, rows)?;
    require_shape(&witness.left_inverse, MatrixRole::LeftInverse, rows, rows)?;
    require_shape(&witness.right, MatrixRole::RightTransform, cols, cols)?;
    require_shape(&witness.right_inverse, MatrixRole::RightInverse, cols, cols)?;

    let retained = source
        .entries
        .len()
        .checked_add(witness.diagonal.entries.len())
        .and_then(|value| value.checked_add(witness.left.entries.len()))
        .and_then(|value| value.checked_add(witness.left_inverse.entries.len()))
        .and_then(|value| value.checked_add(witness.right.entries.len()))
        .and_then(|value| value.checked_add(witness.right_inverse.entries.len()))
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "retained smith witness entries",
        })?;
    if retained > budget.max_retained_entries {
        return Err(IntegralTopologyError::RetainedEntryBudgetExceeded {
            requested: retained,
            max: budget.max_retained_entries,
        });
    }
    let workspace = rows
        .checked_mul(cols)
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "smith workspace entries",
        })?;
    if workspace > budget.max_workspace_entries {
        return Err(IntegralTopologyError::WorkspaceEntryBudgetExceeded {
            requested: workspace,
            max: budget.max_workspace_entries,
        });
    }
    Ok(())
}

fn require_shape(
    matrix: &ExactIntegerMatrix,
    role: MatrixRole,
    expected_rows: usize,
    expected_cols: usize,
) -> Result<(), IntegralTopologyError> {
    if matrix.rows != expected_rows || matrix.cols != expected_cols {
        return Err(IntegralTopologyError::WitnessShape {
            role,
            expected_rows,
            expected_cols,
            actual_rows: matrix.rows,
            actual_cols: matrix.cols,
        });
    }
    Ok(())
}

fn planned_scalar_operations(rows: usize, cols: usize) -> Result<u128, IntegralTopologyError> {
    let rows = u128::try_from(rows).map_err(|_| IntegralTopologyError::WorkPlanOverflow {
        phase: "row work units",
    })?;
    let cols = u128::try_from(cols).map_err(|_| IntegralTopologyError::WorkPlanOverflow {
        phase: "column work units",
    })?;
    let left_inverse = rows
        .checked_mul(rows)
        .and_then(|value| value.checked_mul(rows))
        .and_then(|value| value.checked_mul(2))
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "left inverse verification work",
        })?;
    let right_inverse = cols
        .checked_mul(cols)
        .and_then(|value| value.checked_mul(cols))
        .and_then(|value| value.checked_mul(2))
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "right inverse verification work",
        })?;
    let left_times_source = rows
        .checked_mul(rows)
        .and_then(|value| value.checked_mul(cols))
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "left transform work",
        })?;
    let diagonal_transform = rows
        .checked_mul(cols)
        .and_then(|value| value.checked_mul(cols))
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "right transform work",
        })?;
    left_inverse
        .checked_add(right_inverse)
        .and_then(|value| value.checked_add(left_times_source))
        .and_then(|value| value.checked_add(diagonal_transform))
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "total smith verification work",
        })
}

fn verify_identity_product(
    left: &ExactIntegerMatrix,
    right: &ExactIntegerMatrix,
    stage: SmithWitnessStage,
    checkpoint: &mut impl FnMut(&'static str) -> bool,
    completed: &mut u128,
    planned: u128,
) -> Result<(), IntegralTopologyError> {
    for row in 0..left.rows {
        for col in 0..right.cols {
            poll(checkpoint, "smith inverse product", *completed, planned)?;
            let actual = checked_dot(left, row, right, col, stage, completed)?;
            let expected = i128::from(row == col);
            if actual != expected {
                return Err(IntegralTopologyError::WitnessProductMismatch {
                    stage,
                    row,
                    col,
                    expected,
                    actual,
                });
            }
        }
    }
    Ok(())
}

fn checked_dot(
    left: &ExactIntegerMatrix,
    row: usize,
    right: &ExactIntegerMatrix,
    col: usize,
    stage: SmithWitnessStage,
    completed: &mut u128,
) -> Result<i128, IntegralTopologyError> {
    debug_assert_eq!(left.cols, right.rows);
    let mut sum = 0_i128;
    for term in 0..left.cols {
        let product = left
            .entry(row, term)
            .checked_mul(right.entry(term, col))
            .ok_or(IntegralTopologyError::ArithmeticOverflow {
                stage,
                row,
                col,
                term,
            })?;
        sum = sum
            .checked_add(product)
            .ok_or(IntegralTopologyError::ArithmeticOverflow {
                stage,
                row,
                col,
                term,
            })?;
        *completed = completed
            .checked_add(1)
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "completed scalar operations",
            })?;
    }
    Ok(sum)
}

fn verify_canonical_diagonal(
    diagonal: &ExactIntegerMatrix,
    checkpoint: &mut impl FnMut(&'static str) -> bool,
    completed: u128,
    planned: u128,
) -> Result<Vec<i128>, IntegralTopologyError> {
    for row in 0..diagonal.rows {
        poll(checkpoint, "smith canonical diagonal", completed, planned)?;
        for col in 0..diagonal.cols {
            if row != col {
                let value = diagonal.entry(row, col);
                if value != 0 {
                    return Err(IntegralTopologyError::OffDiagonalEntry { row, col, value });
                }
            }
        }
    }

    let mut invariant_factors = Vec::new();
    invariant_factors
        .try_reserve_exact(diagonal.rows.min(diagonal.cols))
        .map_err(|_| IntegralTopologyError::AllocationRefused {
            phase: "smith invariant factors",
            requested_entries: diagonal.rows.min(diagonal.cols),
        })?;
    let mut zero_seen = false;
    poll(checkpoint, "smith invariant factors", completed, planned)?;
    for index in 0..diagonal.rows.min(diagonal.cols) {
        let value = diagonal.entry(index, index);
        if value < 0 {
            return Err(IntegralTopologyError::NegativeInvariantFactor { index, value });
        }
        if value == 0 {
            zero_seen = true;
            continue;
        }
        if zero_seen {
            return Err(IntegralTopologyError::NonzeroAfterZero { index, value });
        }
        if let Some(previous) = invariant_factors.last().copied()
            && !value.is_multiple_of(previous)
        {
            return Err(IntegralTopologyError::InvariantFactorDivisibility {
                index,
                previous,
                value,
            });
        }
        invariant_factors.push(value);
    }
    Ok(invariant_factors)
}

fn allocate_zeroed(
    entries: usize,
    phase: &'static str,
) -> Result<Vec<i128>, IntegralTopologyError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(entries)
        .map_err(|_| IntegralTopologyError::AllocationRefused {
            phase,
            requested_entries: entries,
        })?;
    values.resize(entries, 0);
    Ok(values)
}

fn poll(
    checkpoint: &mut impl FnMut(&'static str) -> bool,
    phase: &'static str,
    completed: u128,
    planned: u128,
) -> Result<(), IntegralTopologyError> {
    if checkpoint(phase) {
        Ok(())
    } else {
        Err(IntegralTopologyError::Cancelled {
            phase,
            completed_scalar_operations: completed,
            planned_scalar_operations: planned,
        })
    }
}

/// Structured fail-closed exact-algebra refusal.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IntegralTopologyError {
    /// Matrix extent exceeded its explicit envelope.
    MatrixExtentExceeded {
        /// Supplied rows.
        rows: usize,
        /// Supplied columns.
        cols: usize,
        /// Maximum rows.
        max_rows: usize,
        /// Maximum columns.
        max_cols: usize,
    },
    /// `rows * cols` exceeded the per-matrix entry envelope.
    MatrixEntryBudgetExceeded {
        /// Requested entries.
        requested: usize,
        /// Maximum entries.
        max: usize,
    },
    /// Row-major entry count did not match the declared rectangle.
    MatrixEntryCount {
        /// Declared rows.
        rows: usize,
        /// Declared columns.
        cols: usize,
        /// Required entries.
        expected: usize,
        /// Supplied entries.
        actual: usize,
    },
    /// A previously admitted retained matrix exceeds the verification budget.
    RetainedMatrixExceedsBudget {
        /// Matrix role.
        role: MatrixRole,
        /// Retained rows.
        rows: usize,
        /// Retained columns.
        cols: usize,
        /// Retained entries.
        entries: usize,
    },
    /// Complete retained source/witness storage exceeded the envelope.
    RetainedEntryBudgetExceeded {
        /// Requested retained entries.
        requested: usize,
        /// Maximum retained entries.
        max: usize,
    },
    /// Verification scratch storage exceeded the envelope.
    WorkspaceEntryBudgetExceeded {
        /// Requested scratch entries.
        requested: usize,
        /// Maximum scratch entries.
        max: usize,
    },
    /// One witness matrix had the wrong exact shape.
    WitnessShape {
        /// Matrix role.
        role: MatrixRole,
        /// Required rows.
        expected_rows: usize,
        /// Required columns.
        expected_cols: usize,
        /// Supplied rows.
        actual_rows: usize,
        /// Supplied columns.
        actual_cols: usize,
    },
    /// Exact scalar work exceeded its admitted cap.
    ScalarWorkBudgetExceeded {
        /// Planned scalar operations.
        requested: u128,
        /// Maximum scalar operations.
        max: u128,
    },
    /// Work accounting overflowed before execution.
    WorkPlanOverflow {
        /// Refusing phase.
        phase: &'static str,
    },
    /// Internal exact workspace allocation refused.
    AllocationRefused {
        /// Refusing allocation phase.
        phase: &'static str,
        /// Requested entries.
        requested_entries: usize,
    },
    /// Cancellation was observed before any verified value was published.
    Cancelled {
        /// Observation phase.
        phase: &'static str,
        /// Exact completed scalar operations.
        completed_scalar_operations: u128,
        /// Exact planned scalar operations.
        planned_scalar_operations: u128,
    },
    /// Checked `i128` multiplication or addition overflowed.
    ArithmeticOverflow {
        /// Product being checked.
        stage: SmithWitnessStage,
        /// Output row.
        row: usize,
        /// Output column.
        col: usize,
        /// Inner-product term.
        term: usize,
    },
    /// An inverse or transformed-product witness disagreed exactly.
    WitnessProductMismatch {
        /// Product being checked.
        stage: SmithWitnessStage,
        /// Output row.
        row: usize,
        /// Output column.
        col: usize,
        /// Exact expected integer.
        expected: i128,
        /// Exact observed integer.
        actual: i128,
    },
    /// Claimed Smith matrix contained a nonzero off-diagonal entry.
    OffDiagonalEntry {
        /// Matrix row.
        row: usize,
        /// Matrix column.
        col: usize,
        /// Rejected value.
        value: i128,
    },
    /// Canonical invariant factors cannot be negative.
    NegativeInvariantFactor {
        /// Diagonal index.
        index: usize,
        /// Rejected value.
        value: i128,
    },
    /// A positive invariant appeared after the diagonal's zero suffix began.
    NonzeroAfterZero {
        /// Diagonal index.
        index: usize,
        /// Rejected value.
        value: i128,
    },
    /// Consecutive positive invariant factors violated divisibility.
    InvariantFactorDivisibility {
        /// Later diagonal index.
        index: usize,
        /// Previous factor.
        previous: i128,
        /// Rejected later factor.
        value: i128,
    },
}

impl IntegralTopologyError {
    /// Distinguish a mathematical/structural counterexample from an
    /// inconclusive resource or arithmetic refusal.
    #[must_use]
    pub const fn failure_class(&self) -> IntegralTopologyFailureClass {
        match self {
            Self::MatrixEntryCount { .. }
            | Self::WitnessShape { .. }
            | Self::WitnessProductMismatch { .. }
            | Self::OffDiagonalEntry { .. }
            | Self::NegativeInvariantFactor { .. }
            | Self::NonzeroAfterZero { .. }
            | Self::InvariantFactorDivisibility { .. } => IntegralTopologyFailureClass::Refuted,
            Self::MatrixExtentExceeded { .. }
            | Self::MatrixEntryBudgetExceeded { .. }
            | Self::RetainedMatrixExceedsBudget { .. }
            | Self::RetainedEntryBudgetExceeded { .. }
            | Self::WorkspaceEntryBudgetExceeded { .. }
            | Self::ScalarWorkBudgetExceeded { .. }
            | Self::WorkPlanOverflow { .. }
            | Self::AllocationRefused { .. }
            | Self::Cancelled { .. }
            | Self::ArithmeticOverflow { .. } => IntegralTopologyFailureClass::Unknown,
        }
    }
}

impl fmt::Display for IntegralTopologyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "integral topology admission refused: {self:?}")
    }
}

impl core::error::Error for IntegralTopologyError {}
