//! Bounded exact-integer topology algebra for I13.2b.
//!
//! The first tranche is an independently replayable Smith-normal-form
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
//! The second tranche extracts the exact cellular boundary for one admitted
//! [`crate::terminal_relative::TerminalRelativePair`], phase, component, and
//! degree while retaining both canonical quotient bases. It includes the
//! unaugmented zero maps at the two chain-complex edges.
//!
//! The third tranche verifies adjacent pair-bound boundaries and computes
//! `V^-1 * A_(k+1)` for one exact outgoing Smith witness. It refuses a nonzero
//! nonkernel prefix and retains only the witness-relative kernel-coordinate
//! image.
//!
//! The fourth tranche deterministically constructs a complete Smith witness
//! using checked elementary unimodular operations, explicit coefficient and
//! data-dependent work caps, and row-major minimum-magnitude pivots. It never
//! publishes construction output until the independent witness verifier above
//! accepts every inverse, product, and canonical-divisibility obligation.
//!
//! The fifth tranche binds a verified Smith form byte-for-byte to the retained
//! lower kernel image and publishes the exact phase-local cellular quotient
//! decomposition into free rank and nontrivial torsion invariant factors.
//! It retains both complete authorities but still makes no generator, period,
//! naturality, embedding, or physical-R3 winding claim.
//!
//! The sixth tranche lifts every nontrivial torsion and free presentation
//! column back into the original pair-bound `C_k` basis. It also retains the
//! torsion filling chains and replays the exact cycle and finite-order
//! boundary equations before publication. These generators remain relative
//! to the admitted Smith witnesses; no canonical-basis or physical claim is
//! added.
//!
//! The seventh tranche constructs the integral cocycles dual to the retained
//! free generators and checks their complete integer period-pairing matrix.
//! This is the witness-relative free `Hom(H_k, Z)` sector only: torsion
//! cohomology, linking, harmonic representatives, continuum periods, and
//! physical winding remain outside its authority.

use core::fmt;

use crate::terminal_relative::{
    CellRef, ConductorComponentId, IncidenceSign, MAX_TERMINAL_RELATIVE_CELLS,
    MAX_TERMINAL_RELATIVE_INCIDENCES, PhaseId, TerminalRelativePair, TerminalRelativePairId,
};

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
/// Default maximum component-cell visits across both canonical basis scans.
pub const DEFAULT_MAX_BOUNDARY_COMPONENT_VISITS: usize = 2 * MAX_TERMINAL_RELATIVE_CELLS;
/// Default maximum admitted incidences visited while binding one pair boundary.
pub const DEFAULT_MAX_BOUNDARY_INCIDENCE_VISITS: usize = MAX_TERMINAL_RELATIVE_INCIDENCES;
/// Default retained entries across a pair boundary matrix and its two bases.
pub const DEFAULT_MAX_BOUNDARY_RETAINED_ENTRIES: usize =
    DEFAULT_MAX_MATRIX_ENTRIES + 2 * DEFAULT_MAX_MATRIX_EXTENT;
/// Default retained integer/cell entries for adjacent-boundary transport.
pub const DEFAULT_MAX_KERNEL_RETAINED_ENTRIES: usize =
    10 * DEFAULT_MAX_MATRIX_ENTRIES + 4 * DEFAULT_MAX_MATRIX_EXTENT;
/// Default exact binding comparisons for adjacent-boundary transport.
pub const DEFAULT_MAX_KERNEL_BINDING_ITEMS: usize =
    DEFAULT_MAX_MATRIX_ENTRIES + DEFAULT_MAX_MATRIX_EXTENT;
/// Default maximum simultaneously live source, diagonal, and transform entries
/// during constructive Smith reduction.
pub const DEFAULT_MAX_CONSTRUCTION_LIVE_ENTRIES: usize = DEFAULT_MAX_RETAINED_ENTRIES;
/// Default maximum elementary unimodular row/column operations.
pub const DEFAULT_MAX_CONSTRUCTION_OPERATIONS: u128 = 1_000_000;
/// Default maximum admitted source/diagonal inspections and destination-slot
/// updates during reduction.
pub const DEFAULT_MAX_CONSTRUCTION_ENTRY_STEPS: u128 = DEFAULT_MAX_SCALAR_OPERATIONS;
/// Default exact source/invariant comparisons for quotient-homology binding.
pub const DEFAULT_MAX_HOMOLOGY_BINDING_ITEMS: usize =
    DEFAULT_MAX_MATRIX_ENTRIES + DEFAULT_MAX_MATRIX_EXTENT;
/// Default retained integer/cell entries across kernel transport and the
/// lower-image Smith authority.
pub const DEFAULT_MAX_HOMOLOGY_RETAINED_ENTRIES: usize =
    20 * DEFAULT_MAX_MATRIX_ENTRIES + 8 * DEFAULT_MAX_MATRIX_EXTENT;
/// Default maximum generator and torsion-filling coefficients retained at
/// once.
pub const DEFAULT_MAX_GENERATOR_OUTPUT_ENTRIES: usize = DEFAULT_MAX_MATRIX_ENTRIES;
/// Default retained entries across homology authority and generator lift.
pub const DEFAULT_MAX_GENERATOR_RETAINED_ENTRIES: usize =
    DEFAULT_MAX_HOMOLOGY_RETAINED_ENTRIES + DEFAULT_MAX_GENERATOR_OUTPUT_ENTRIES;
/// Default maximum free-cocycle and period-pairing coefficients retained at
/// once.
pub const DEFAULT_MAX_COCYCLE_OUTPUT_ENTRIES: usize = 2 * DEFAULT_MAX_MATRIX_ENTRIES;
/// Default retained entries across generator authority, free cocycles, and
/// their period-pairing matrix.
pub const DEFAULT_MAX_COCYCLE_RETAINED_ENTRIES: usize =
    DEFAULT_MAX_GENERATOR_RETAINED_ENTRIES + DEFAULT_MAX_COCYCLE_OUTPUT_ENTRIES;

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

/// Explicit resource envelope for one terminal-relative boundary extraction.
///
/// Component and incidence visits are capped independently because a sparse
/// admitted complex can be much larger than the requested phase-local matrix.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TerminalRelativeBoundaryBudget {
    max_rows: usize,
    max_cols: usize,
    max_matrix_entries: usize,
    max_retained_entries: usize,
    max_component_visits: usize,
    max_incidence_visits: usize,
}

impl TerminalRelativeBoundaryBudget {
    /// Construct an exact extraction envelope.
    #[must_use]
    pub const fn new(
        max_rows: usize,
        max_cols: usize,
        max_matrix_entries: usize,
        max_retained_entries: usize,
        max_component_visits: usize,
        max_incidence_visits: usize,
    ) -> Self {
        Self {
            max_rows,
            max_cols,
            max_matrix_entries,
            max_retained_entries,
            max_component_visits,
            max_incidence_visits,
        }
    }

    /// Maximum target-basis rows.
    #[must_use]
    pub const fn max_rows(self) -> usize {
        self.max_rows
    }

    /// Maximum source-basis columns.
    #[must_use]
    pub const fn max_cols(self) -> usize {
        self.max_cols
    }

    /// Maximum dense matrix entries.
    #[must_use]
    pub const fn max_matrix_entries(self) -> usize {
        self.max_matrix_entries
    }

    /// Maximum entries retained across matrix and bases.
    #[must_use]
    pub const fn max_retained_entries(self) -> usize {
        self.max_retained_entries
    }

    /// Maximum component-support visits across both canonical basis scans.
    #[must_use]
    pub const fn max_component_visits(self) -> usize {
        self.max_component_visits
    }

    /// Maximum ambient incidence rows visited.
    #[must_use]
    pub const fn max_incidence_visits(self) -> usize {
        self.max_incidence_visits
    }
}

impl Default for TerminalRelativeBoundaryBudget {
    fn default() -> Self {
        Self::new(
            DEFAULT_MAX_MATRIX_EXTENT,
            DEFAULT_MAX_MATRIX_EXTENT,
            DEFAULT_MAX_MATRIX_ENTRIES,
            DEFAULT_MAX_BOUNDARY_RETAINED_ENTRIES,
            DEFAULT_MAX_BOUNDARY_COMPONENT_VISITS,
            DEFAULT_MAX_BOUNDARY_INCIDENCE_VISITS,
        )
    }
}

/// Resource envelope for transporting an incoming boundary into one verified
/// Smith kernel coordinate system.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KernelCoordinateBudget {
    max_extent: usize,
    max_output_entries: usize,
    max_retained_entries: usize,
    max_binding_items: usize,
    max_scalar_operations: u128,
}

impl KernelCoordinateBudget {
    /// Construct an adjacent-boundary transport envelope.
    #[must_use]
    pub const fn new(
        max_extent: usize,
        max_output_entries: usize,
        max_retained_entries: usize,
        max_binding_items: usize,
        max_scalar_operations: u128,
    ) -> Self {
        Self {
            max_extent,
            max_output_entries,
            max_retained_entries,
            max_binding_items,
            max_scalar_operations,
        }
    }

    /// Maximum admitted outgoing row, shared-chain, or incoming-column extent.
    #[must_use]
    pub const fn max_extent(self) -> usize {
        self.max_extent
    }

    /// Maximum entries retained in the lower kernel-coordinate image.
    #[must_use]
    pub const fn max_output_entries(self) -> usize {
        self.max_output_entries
    }

    /// Maximum retained integer/cell entries across complete input authority
    /// and the lower image. Stable-identity bytes are bounded by pair admission
    /// but are not counted as entries here.
    #[must_use]
    pub const fn max_retained_entries(self) -> usize {
        self.max_retained_entries
    }

    /// Maximum exact source/basis comparisons before multiplication.
    #[must_use]
    pub const fn max_binding_items(self) -> usize {
        self.max_binding_items
    }

    /// Maximum checked dot-product terms in `V^-1 * B`.
    #[must_use]
    pub const fn max_scalar_operations(self) -> u128 {
        self.max_scalar_operations
    }
}

impl Default for KernelCoordinateBudget {
    fn default() -> Self {
        Self::new(
            DEFAULT_MAX_MATRIX_EXTENT,
            DEFAULT_MAX_MATRIX_ENTRIES,
            DEFAULT_MAX_KERNEL_RETAINED_ENTRIES,
            DEFAULT_MAX_KERNEL_BINDING_ITEMS,
            DEFAULT_MAX_SCALAR_OPERATIONS,
        )
    }
}

/// Data-dependent resource envelope for deterministic constructive Smith
/// reduction. Final witness verification remains independently bounded by an
/// [`ExactAlgebraBudget`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SmithConstructionBudget {
    max_live_entries: usize,
    max_elementary_operations: u128,
    max_entry_steps: u128,
    max_abs_coefficient: u128,
}

impl SmithConstructionBudget {
    /// Construct a reduction envelope. Excluding bounded workspace
    /// initialization, entry steps count admitted-source and active-diagonal
    /// inspections plus each destination slot updated or swapped. Source reads
    /// inside one already-reserved destination update are not counted again.
    #[must_use]
    pub const fn new(
        max_live_entries: usize,
        max_elementary_operations: u128,
        max_entry_steps: u128,
        max_abs_coefficient: u128,
    ) -> Self {
        Self {
            max_live_entries,
            max_elementary_operations,
            max_entry_steps,
            max_abs_coefficient,
        }
    }

    /// Maximum simultaneously live entries across the retained source and
    /// five construction matrices.
    #[must_use]
    pub const fn max_live_entries(self) -> usize {
        self.max_live_entries
    }

    /// Maximum row swaps, column swaps, Euclidean reductions, repairs, and
    /// sign normalizations.
    #[must_use]
    pub const fn max_elementary_operations(self) -> u128 {
        self.max_elementary_operations
    }

    /// Maximum admitted source/diagonal inspections and destination updates.
    #[must_use]
    pub const fn max_entry_steps(self) -> u128 {
        self.max_entry_steps
    }

    /// Maximum unsigned magnitude of any source or generated coefficient.
    /// A cap no greater than `i128::MAX` excludes the unrepresentable positive
    /// magnitude of `i128::MIN`.
    #[must_use]
    pub const fn max_abs_coefficient(self) -> u128 {
        self.max_abs_coefficient
    }
}

impl Default for SmithConstructionBudget {
    fn default() -> Self {
        Self::new(
            DEFAULT_MAX_CONSTRUCTION_LIVE_ENTRIES,
            DEFAULT_MAX_CONSTRUCTION_OPERATIONS,
            DEFAULT_MAX_CONSTRUCTION_ENTRY_STEPS,
            i128::MAX.unsigned_abs(),
        )
    }
}

/// Resource envelope for binding a lower-image Smith authority to one verified
/// terminal-relative kernel transport.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HomologyDecompositionBudget {
    max_binding_items: usize,
    max_retained_entries: usize,
}

impl HomologyDecompositionBudget {
    /// Construct a quotient-decomposition envelope.
    #[must_use]
    pub const fn new(max_binding_items: usize, max_retained_entries: usize) -> Self {
        Self {
            max_binding_items,
            max_retained_entries,
        }
    }

    /// Maximum exact lower-source and invariant-factor inspections.
    #[must_use]
    pub const fn max_binding_items(self) -> usize {
        self.max_binding_items
    }

    /// Maximum retained integer/cell entries across both complete authorities.
    #[must_use]
    pub const fn max_retained_entries(self) -> usize {
        self.max_retained_entries
    }
}

impl Default for HomologyDecompositionBudget {
    fn default() -> Self {
        Self::new(
            DEFAULT_MAX_HOMOLOGY_BINDING_ITEMS,
            DEFAULT_MAX_HOMOLOGY_RETAINED_ENTRIES,
        )
    }
}

/// Resource envelope for lifting quotient-presentation generators into the
/// original pair-bound chain basis and verifying their obligations.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HomologyGeneratorBudget {
    max_output_entries: usize,
    max_retained_entries: usize,
    max_scalar_operations: u128,
}

impl HomologyGeneratorBudget {
    /// Construct a generator-lift envelope.
    #[must_use]
    pub const fn new(
        max_output_entries: usize,
        max_retained_entries: usize,
        max_scalar_operations: u128,
    ) -> Self {
        Self {
            max_output_entries,
            max_retained_entries,
            max_scalar_operations,
        }
    }

    /// Maximum entries across the retained original-chain generator and
    /// torsion-filling matrices.
    #[must_use]
    pub const fn max_output_entries(self) -> usize {
        self.max_output_entries
    }

    /// Maximum retained entries across the complete homology authority and
    /// both generator-witness matrices.
    #[must_use]
    pub const fn max_retained_entries(self) -> usize {
        self.max_retained_entries
    }

    /// Maximum checked lift, cycle, and quotient-order scalar terms.
    #[must_use]
    pub const fn max_scalar_operations(self) -> u128 {
        self.max_scalar_operations
    }
}

impl Default for HomologyGeneratorBudget {
    fn default() -> Self {
        Self::new(
            DEFAULT_MAX_GENERATOR_OUTPUT_ENTRIES,
            DEFAULT_MAX_GENERATOR_RETAINED_ENTRIES,
            DEFAULT_MAX_SCALAR_OPERATIONS,
        )
    }
}

/// Resource envelope for exact free integral cocycles and their period duals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HomologyCocycleBudget {
    max_output_entries: usize,
    max_retained_entries: usize,
    max_scalar_operations: u128,
}

impl HomologyCocycleBudget {
    /// Construct a free-cocycle/period envelope.
    #[must_use]
    pub const fn new(
        max_output_entries: usize,
        max_retained_entries: usize,
        max_scalar_operations: u128,
    ) -> Self {
        Self {
            max_output_entries,
            max_retained_entries,
            max_scalar_operations,
        }
    }

    /// Maximum entries across cocycle and period-pairing matrices.
    #[must_use]
    pub const fn max_output_entries(self) -> usize {
        self.max_output_entries
    }

    /// Maximum retained entries including the complete generator authority.
    #[must_use]
    pub const fn max_retained_entries(self) -> usize {
        self.max_retained_entries
    }

    /// Maximum checked cocycle-lift, coboundary, and pairing scalar terms.
    #[must_use]
    pub const fn max_scalar_operations(self) -> u128 {
        self.max_scalar_operations
    }
}

impl Default for HomologyCocycleBudget {
    fn default() -> Self {
        Self::new(
            DEFAULT_MAX_COCYCLE_OUTPUT_ENTRIES,
            DEFAULT_MAX_COCYCLE_RETAINED_ENTRIES,
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

/// Canonically bound terminal-relative boundary matrix.
///
/// Rows are the phase-local quotient basis in degree `k - 1` and columns are
/// the basis in degree `k`. The unaugmented edge maps are represented too:
/// degree zero has no target basis (`0 x dim(C_0)`), while degree one above
/// the complex has no source basis (`dim(C_top) x 0`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TerminalRelativeBoundaryMatrix {
    pair: TerminalRelativePairId,
    phase: PhaseId,
    component: ConductorComponentId,
    source_degree: u8,
    target_degree: Option<u8>,
    source_basis: Vec<CellRef>,
    target_basis: Vec<CellRef>,
    matrix: ExactIntegerMatrix,
    work_items: u128,
}

impl TerminalRelativeBoundaryMatrix {
    /// Strong identity of the complete admitted pair.
    #[must_use]
    pub const fn pair_id(&self) -> TerminalRelativePairId {
        self.pair
    }

    /// Phase whose component-local quotient basis was extracted.
    #[must_use]
    pub const fn phase(&self) -> &PhaseId {
        &self.phase
    }

    /// Component bound to the retained phase by pair admission.
    #[must_use]
    pub const fn component(&self) -> &ConductorComponentId {
        &self.component
    }

    /// Source chain degree `k`.
    #[must_use]
    pub const fn source_degree(&self) -> u8 {
        self.source_degree
    }

    /// Target chain degree, or `None` for the unaugmented degree-zero edge.
    #[must_use]
    pub const fn target_degree(&self) -> Option<u8> {
        self.target_degree
    }

    /// Canonically ordered source cells corresponding to matrix columns.
    #[must_use]
    pub fn source_basis(&self) -> &[CellRef] {
        &self.source_basis
    }

    /// Canonically ordered target cells corresponding to matrix rows.
    #[must_use]
    pub fn target_basis(&self) -> &[CellRef] {
        &self.target_basis
    }

    /// Exact row-major boundary matrix.
    #[must_use]
    pub const fn matrix(&self) -> &ExactIntegerMatrix {
        &self.matrix
    }

    /// Deterministic component-cell and incidence visits completed.
    #[must_use]
    pub const fn work_items(&self) -> u128 {
        self.work_items
    }

    /// This binds exact admitted incidence, not a homology or winding result.
    #[must_use]
    pub const fn applicability(&self) -> TopologyApplicability {
        TopologyApplicability::TerminalRelativeIncidenceOnly
    }
}

/// Extract a canonical terminal-relative boundary without injected cancellation.
pub fn extract_terminal_relative_boundary_matrix(
    pair: &TerminalRelativePair,
    phase: &PhaseId,
    source_degree: u8,
    budget: TerminalRelativeBoundaryBudget,
) -> Result<TerminalRelativeBoundaryMatrix, IntegralTopologyError> {
    extract_terminal_relative_boundary_matrix_with_checkpoint(
        pair,
        phase,
        source_degree,
        budget,
        &mut |_| true,
    )
}

/// Extract a canonical pair/phase/degree boundary under bounded polling.
///
/// This is the intrinsic cellular incidence. Terminal-current orientation,
/// port trivialization, and phase-current relabel signs are deliberately not
/// applied. The callback runs before both deterministic component scans,
/// before every ambient incidence visit, and before final publication.
#[allow(clippy::too_many_lines)]
pub fn extract_terminal_relative_boundary_matrix_with_checkpoint(
    pair: &TerminalRelativePair,
    phase: &PhaseId,
    source_degree: u8,
    budget: TerminalRelativeBoundaryBudget,
    checkpoint: &mut impl FnMut(&'static str) -> bool,
) -> Result<TerminalRelativeBoundaryMatrix, IntegralTopologyError> {
    let Some(component_id) = pair.phase_component(phase) else {
        return Err(IntegralTopologyError::UnknownTerminalRelativePhase {
            phase: phase.as_str().to_owned(),
        });
    };
    let max_source_degree = pair.complex().dimension().checked_add(1).ok_or(
        IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative edge degree",
        },
    )?;
    if source_degree > max_source_degree {
        return Err(IntegralTopologyError::BoundaryDegreeOutOfRange {
            degree: source_degree,
            max: max_source_degree,
        });
    }
    let Some(component) = pair
        .components()
        .iter()
        .find(|component| component.id() == component_id)
    else {
        return Err(IntegralTopologyError::PhaseComponentBindingLost {
            phase: phase.as_str().to_owned(),
            component: component_id.as_str().to_owned(),
        });
    };
    let component_cells = component.support().cells().len();
    let component_visits =
        component_cells
            .checked_mul(2)
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "terminal-relative component visits",
            })?;
    if component_visits > budget.max_component_visits {
        return Err(IntegralTopologyError::ComponentVisitBudgetExceeded {
            requested: component_visits,
            max: budget.max_component_visits,
        });
    }
    let incidence_visits = pair.complex().incidences().len();
    if incidence_visits > budget.max_incidence_visits {
        return Err(IntegralTopologyError::IncidenceVisitBudgetExceeded {
            requested: incidence_visits,
            max: budget.max_incidence_visits,
        });
    }
    let planned = u128::try_from(component_visits)
        .ok()
        .and_then(|visits| {
            u128::try_from(incidence_visits)
                .ok()
                .and_then(|incidences| visits.checked_add(incidences))
        })
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative boundary extraction",
        })?;
    let retained_phase = phase.clone();
    let retained_component = component_id.clone();
    poll_pair_boundary(
        checkpoint,
        "terminal-relative boundary preflight",
        0,
        planned,
    )?;

    let target_degree = source_degree.checked_sub(1);
    let mut completed = 0_u128;
    let mut source_count = 0_usize;
    let mut target_count = 0_usize;
    for cell in component.support().cells() {
        poll_pair_boundary(
            checkpoint,
            "terminal-relative basis count",
            completed,
            planned,
        )?;
        if !pair.relative().cells().contains(cell) {
            if cell.degree() == source_degree {
                source_count =
                    source_count
                        .checked_add(1)
                        .ok_or(IntegralTopologyError::WorkPlanOverflow {
                            phase: "terminal-relative source basis count",
                        })?;
            }
            if target_degree.is_some_and(|degree| cell.degree() == degree) {
                target_count =
                    target_count
                        .checked_add(1)
                        .ok_or(IntegralTopologyError::WorkPlanOverflow {
                            phase: "terminal-relative target basis count",
                        })?;
            }
        }
        complete_pair_boundary_item(&mut completed)?;
    }
    preflight_boundary_output(target_count, source_count, budget)?;

    let mut source_basis = allocate_cell_refs(source_count, "terminal-relative source basis")?;
    let mut target_basis = allocate_cell_refs(target_count, "terminal-relative target basis")?;
    for cell in component.support().cells() {
        poll_pair_boundary(
            checkpoint,
            "terminal-relative basis materialization",
            completed,
            planned,
        )?;
        if !pair.relative().cells().contains(cell) {
            if cell.degree() == source_degree {
                source_basis.push(*cell);
            }
            if target_degree.is_some_and(|degree| cell.degree() == degree) {
                target_basis.push(*cell);
            }
        }
        complete_pair_boundary_item(&mut completed)?;
    }

    let matrix_entries =
        target_count
            .checked_mul(source_count)
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "terminal-relative boundary entries",
            })?;
    let mut entries = allocate_zeroed(matrix_entries, "terminal-relative boundary matrix")?;
    for incidence in pair.complex().incidences() {
        poll_pair_boundary(
            checkpoint,
            "terminal-relative incidence projection",
            completed,
            planned,
        )?;
        if let (Ok(row), Ok(col)) = (
            target_basis.binary_search(&incidence.lower()),
            source_basis.binary_search(&incidence.upper()),
        ) {
            entries[row * source_count + col] = match incidence.sign() {
                IncidenceSign::Negative => -1,
                IncidenceSign::Positive => 1,
            };
        }
        complete_pair_boundary_item(&mut completed)?;
    }
    poll_pair_boundary(
        checkpoint,
        "terminal-relative boundary finalize",
        completed,
        planned,
    )?;
    debug_assert_eq!(completed, planned);

    Ok(TerminalRelativeBoundaryMatrix {
        pair: pair.identity(),
        phase: retained_phase,
        component: retained_component,
        source_degree,
        target_degree,
        source_basis,
        target_basis,
        matrix: ExactIntegerMatrix {
            rows: target_count,
            cols: source_count,
            entries,
        },
        work_items: completed,
    })
}

fn preflight_boundary_output(
    rows: usize,
    cols: usize,
    budget: TerminalRelativeBoundaryBudget,
) -> Result<(), IntegralTopologyError> {
    if rows > budget.max_rows || cols > budget.max_cols {
        return Err(IntegralTopologyError::MatrixExtentExceeded {
            rows,
            cols,
            max_rows: budget.max_rows,
            max_cols: budget.max_cols,
        });
    }
    let matrix_entries = rows
        .checked_mul(cols)
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative boundary entry count",
        })?;
    if matrix_entries > budget.max_matrix_entries {
        return Err(IntegralTopologyError::MatrixEntryBudgetExceeded {
            requested: matrix_entries,
            max: budget.max_matrix_entries,
        });
    }
    let retained = matrix_entries
        .checked_add(rows)
        .and_then(|entries| entries.checked_add(cols))
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative boundary retained entries",
        })?;
    if retained > budget.max_retained_entries {
        return Err(IntegralTopologyError::RetainedEntryBudgetExceeded {
            requested: retained,
            max: budget.max_retained_entries,
        });
    }
    Ok(())
}

fn allocate_cell_refs(
    entries: usize,
    phase: &'static str,
) -> Result<Vec<CellRef>, IntegralTopologyError> {
    let mut cells = Vec::new();
    cells
        .try_reserve_exact(entries)
        .map_err(|_| IntegralTopologyError::AllocationRefused {
            phase,
            requested_entries: entries,
        })?;
    Ok(cells)
}

fn complete_pair_boundary_item(completed: &mut u128) -> Result<(), IntegralTopologyError> {
    *completed = completed
        .checked_add(1)
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "completed terminal-relative boundary work",
        })?;
    Ok(())
}

fn poll_pair_boundary(
    checkpoint: &mut impl FnMut(&'static str) -> bool,
    phase: &'static str,
    completed: u128,
    planned: u128,
) -> Result<(), IntegralTopologyError> {
    if checkpoint(phase) {
        Ok(())
    } else {
        Err(IntegralTopologyError::PairBoundaryCancelled {
            phase,
            completed_work_items: completed,
            planned_work_items: planned,
        })
    }
}

/// Exact incoming-boundary image in one verified outgoing Smith kernel basis.
///
/// The retained rows are coordinates in the basis given by the columns
/// `V[:, rank..]` of the exact Smith witness. They are not `CellRef` rows and
/// are not canonical across different valid Smith witnesses.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VerifiedTerminalRelativeKernelTransport {
    outgoing: TerminalRelativeBoundaryMatrix,
    incoming: TerminalRelativeBoundaryMatrix,
    outgoing_smith: VerifiedSmithNormalForm,
    kernel_image: ExactIntegerMatrix,
    binding_items: usize,
    scalar_operations: u128,
}

impl VerifiedTerminalRelativeKernelTransport {
    /// Pair identity shared by both adjacent boundaries.
    #[must_use]
    pub const fn pair_id(&self) -> TerminalRelativePairId {
        self.outgoing.pair
    }

    /// Shared phase identity.
    #[must_use]
    pub const fn phase(&self) -> &PhaseId {
        &self.outgoing.phase
    }

    /// Shared phase-owned component.
    #[must_use]
    pub const fn component(&self) -> &ConductorComponentId {
        &self.outgoing.component
    }

    /// Homological chain degree `k` whose outgoing map is `A_k`.
    #[must_use]
    pub const fn degree(&self) -> u8 {
        self.outgoing.source_degree
    }

    /// Complete pair-bound outgoing boundary `A_k`.
    #[must_use]
    pub const fn outgoing_boundary(&self) -> &TerminalRelativeBoundaryMatrix {
        &self.outgoing
    }

    /// Complete pair-bound incoming boundary `A_(k+1)`.
    #[must_use]
    pub const fn incoming_boundary(&self) -> &TerminalRelativeBoundaryMatrix {
        &self.incoming
    }

    /// Exact Smith authority that defines the retained kernel coordinates.
    #[must_use]
    pub const fn outgoing_smith(&self) -> &VerifiedSmithNormalForm {
        &self.outgoing_smith
    }

    /// Shared chain-group extent `dim(C_k)`.
    #[must_use]
    pub const fn chain_extent(&self) -> usize {
        self.outgoing.matrix.cols
    }

    /// Verified rank of `A_k`.
    #[must_use]
    pub const fn outgoing_rank(&self) -> usize {
        self.outgoing_smith.rank
    }

    /// Rank of the verified outgoing kernel basis.
    #[must_use]
    pub const fn kernel_dimension(&self) -> usize {
        self.kernel_image.rows
    }

    /// Incoming image in the retained Smith kernel coordinates.
    #[must_use]
    pub const fn kernel_image(&self) -> &ExactIntegerMatrix {
        &self.kernel_image
    }

    /// Exact source/basis comparisons completed before multiplication.
    #[must_use]
    pub const fn binding_items(&self) -> usize {
        self.binding_items
    }

    /// Checked dot-product terms completed for `V^-1 * A_(k+1)`.
    #[must_use]
    pub const fn scalar_operations(&self) -> u128 {
        self.scalar_operations
    }

    /// These are witness-relative kernel coordinates, not homology authority.
    #[must_use]
    pub const fn applicability(&self) -> TopologyApplicability {
        TopologyApplicability::TerminalRelativeKernelCoordinatesOnly
    }
}

/// Verify adjacent pair boundaries and transport the incoming image into the
/// outgoing Smith kernel coordinates without injected cancellation.
#[allow(clippy::large_types_passed_by_value)]
pub fn verify_terminal_relative_kernel_transport(
    outgoing: TerminalRelativeBoundaryMatrix,
    incoming: TerminalRelativeBoundaryMatrix,
    outgoing_smith: VerifiedSmithNormalForm,
    budget: KernelCoordinateBudget,
) -> Result<VerifiedTerminalRelativeKernelTransport, IntegralTopologyError> {
    verify_terminal_relative_kernel_transport_with_checkpoint(
        outgoing,
        incoming,
        outgoing_smith,
        budget,
        &mut |_| true,
    )
}

/// Verify `im(A_(k+1))` lies in `ker(A_k)` and retain its exact coordinates.
///
/// With `U * A_k * V = D`, old and Smith source coordinates satisfy
/// `x = V * y`, so the incoming boundary must be transformed as
/// `V^-1 * A_(k+1)`. Its first `rank(A_k)` rows must vanish exactly. Only the
/// lower rows are retained; the full transformed matrix is never allocated.
#[allow(clippy::too_many_lines)]
#[allow(clippy::large_types_passed_by_value)]
pub fn verify_terminal_relative_kernel_transport_with_checkpoint(
    outgoing: TerminalRelativeBoundaryMatrix,
    incoming: TerminalRelativeBoundaryMatrix,
    outgoing_smith: VerifiedSmithNormalForm,
    budget: KernelCoordinateBudget,
    checkpoint: &mut impl FnMut(&'static str) -> bool,
) -> Result<VerifiedTerminalRelativeKernelTransport, IntegralTopologyError> {
    preflight_adjacent_boundary_bindings(&outgoing, &incoming)?;
    preflight_boundary_internal_shape(&outgoing, "outgoing boundary")?;
    preflight_boundary_internal_shape(&incoming, "incoming boundary")?;

    let outgoing_matrix = outgoing.matrix();
    let incoming_matrix = incoming.matrix();
    let rows = outgoing_matrix.rows;
    let chain_extent = outgoing_matrix.cols;
    let incoming_cols = incoming_matrix.cols;
    if incoming_matrix.rows != chain_extent {
        return Err(IntegralTopologyError::KernelCoordinateInvariantLost {
            field: "adjacent matrix inner extent",
        });
    }
    if outgoing_smith.source.rows != rows || outgoing_smith.source.cols != chain_extent {
        return Err(IntegralTopologyError::OutgoingSmithSourceShapeMismatch {
            expected_rows: rows,
            expected_cols: chain_extent,
            actual_rows: outgoing_smith.source.rows,
            actual_cols: outgoing_smith.source.cols,
        });
    }
    if outgoing_smith.rank > chain_extent {
        return Err(IntegralTopologyError::KernelCoordinateInvariantLost {
            field: "outgoing Smith rank exceeds chain extent",
        });
    }
    if rows > budget.max_extent
        || chain_extent > budget.max_extent
        || incoming_cols > budget.max_extent
    {
        return Err(IntegralTopologyError::KernelCoordinateExtentExceeded {
            outgoing_rows: rows,
            chain_extent,
            incoming_cols,
            max: budget.max_extent,
        });
    }

    let kernel_rows = chain_extent - outgoing_smith.rank;
    let output_entries =
        kernel_rows
            .checked_mul(incoming_cols)
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "kernel-coordinate output entries",
            })?;
    if output_entries > budget.max_output_entries {
        return Err(IntegralTopologyError::MatrixEntryBudgetExceeded {
            requested: output_entries,
            max: budget.max_output_entries,
        });
    }
    let binding_items = outgoing
        .source_basis
        .len()
        .checked_add(outgoing_matrix.entries.len())
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "kernel-coordinate binding items",
        })?;
    if binding_items > budget.max_binding_items {
        return Err(IntegralTopologyError::KernelBindingBudgetExceeded {
            requested: binding_items,
            max: budget.max_binding_items,
        });
    }
    let scalar_operations = planned_kernel_scalar_operations(chain_extent, incoming_cols)?;
    if scalar_operations > budget.max_scalar_operations {
        return Err(IntegralTopologyError::ScalarWorkBudgetExceeded {
            requested: scalar_operations,
            max: budget.max_scalar_operations,
        });
    }
    let retained_entries =
        retained_kernel_transport_entries(&outgoing, &incoming, &outgoing_smith, output_entries)?;
    if retained_entries > budget.max_retained_entries {
        return Err(IntegralTopologyError::RetainedEntryBudgetExceeded {
            requested: retained_entries,
            max: budget.max_retained_entries,
        });
    }

    let mut completed_binding = 0_usize;
    let mut completed_scalar = 0_u128;
    poll_kernel_transport(
        checkpoint,
        "kernel-coordinate preflight",
        completed_binding,
        binding_items,
        completed_scalar,
        scalar_operations,
    )?;
    for (index, (outgoing_cell, incoming_cell)) in outgoing
        .source_basis
        .iter()
        .zip(&incoming.target_basis)
        .enumerate()
    {
        poll_kernel_transport(
            checkpoint,
            "kernel-coordinate shared basis",
            completed_binding,
            binding_items,
            completed_scalar,
            scalar_operations,
        )?;
        if outgoing_cell != incoming_cell {
            return Err(IntegralTopologyError::AdjacentBoundaryBasisMismatch {
                index,
                outgoing: Some(*outgoing_cell),
                incoming: Some(*incoming_cell),
            });
        }
        completed_binding += 1;
    }
    for (index, (expected, actual)) in outgoing_matrix
        .entries
        .iter()
        .zip(&outgoing_smith.source.entries)
        .enumerate()
    {
        poll_kernel_transport(
            checkpoint,
            "kernel-coordinate Smith source binding",
            completed_binding,
            binding_items,
            completed_scalar,
            scalar_operations,
        )?;
        if expected != actual {
            return Err(IntegralTopologyError::OutgoingSmithSourceEntryMismatch {
                row: index / chain_extent.max(1),
                col: index % chain_extent.max(1),
                expected: *expected,
                actual: *actual,
            });
        }
        completed_binding += 1;
    }
    debug_assert_eq!(completed_binding, binding_items);

    poll_kernel_transport(
        checkpoint,
        "kernel-coordinate output allocation",
        completed_binding,
        binding_items,
        completed_scalar,
        scalar_operations,
    )?;
    let mut lower = allocate_zeroed(output_entries, "kernel-coordinate incoming image")?;
    for row in 0..chain_extent {
        for col in 0..incoming_cols {
            poll_kernel_transport(
                checkpoint,
                "kernel-coordinate incoming transform",
                completed_binding,
                binding_items,
                completed_scalar,
                scalar_operations,
            )?;
            let coordinate = checked_dot(
                outgoing_smith.right_inverse(),
                row,
                incoming_matrix,
                col,
                SmithWitnessStage::KernelCoordinateIncoming,
                &mut completed_scalar,
            )?;
            if row < outgoing_smith.rank {
                if coordinate != 0 {
                    return Err(IntegralTopologyError::IncomingImageOutsideKernel {
                        row,
                        col,
                        coordinate,
                        invariant_factor: outgoing_smith.invariant_factors[row],
                    });
                }
            } else {
                lower[(row - outgoing_smith.rank) * incoming_cols + col] = coordinate;
            }
        }
    }
    poll_kernel_transport(
        checkpoint,
        "kernel-coordinate finalize",
        completed_binding,
        binding_items,
        completed_scalar,
        scalar_operations,
    )?;
    debug_assert_eq!(completed_scalar, scalar_operations);

    Ok(VerifiedTerminalRelativeKernelTransport {
        outgoing,
        incoming,
        outgoing_smith,
        kernel_image: ExactIntegerMatrix {
            rows: kernel_rows,
            cols: incoming_cols,
            entries: lower,
        },
        binding_items: completed_binding,
        scalar_operations: completed_scalar,
    })
}

fn preflight_adjacent_boundary_bindings(
    outgoing: &TerminalRelativeBoundaryMatrix,
    incoming: &TerminalRelativeBoundaryMatrix,
) -> Result<(), IntegralTopologyError> {
    for (field, matches) in [
        ("pair identity", outgoing.pair == incoming.pair),
        ("phase identity", outgoing.phase == incoming.phase),
        (
            "component identity",
            outgoing.component == incoming.component,
        ),
    ] {
        if !matches {
            return Err(IntegralTopologyError::AdjacentBoundaryBindingMismatch { field });
        }
    }
    let expected_incoming =
        outgoing
            .source_degree
            .checked_add(1)
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "adjacent boundary degree",
            })?;
    if incoming.source_degree != expected_incoming
        || incoming.target_degree != Some(outgoing.source_degree)
    {
        return Err(IntegralTopologyError::AdjacentBoundaryDegreeMismatch {
            outgoing: outgoing.source_degree,
            incoming_source: incoming.source_degree,
            incoming_target: incoming.target_degree,
        });
    }
    if outgoing.target_degree != outgoing.source_degree.checked_sub(1) {
        return Err(IntegralTopologyError::KernelCoordinateInvariantLost {
            field: "outgoing target degree",
        });
    }
    if outgoing.source_basis.len() != incoming.target_basis.len() {
        let index = outgoing.source_basis.len().min(incoming.target_basis.len());
        return Err(IntegralTopologyError::AdjacentBoundaryBasisMismatch {
            index,
            outgoing: outgoing.source_basis.get(index).copied(),
            incoming: incoming.target_basis.get(index).copied(),
        });
    }
    Ok(())
}

fn preflight_boundary_internal_shape(
    boundary: &TerminalRelativeBoundaryMatrix,
    field: &'static str,
) -> Result<(), IntegralTopologyError> {
    if boundary.matrix.rows != boundary.target_basis.len()
        || boundary.matrix.cols != boundary.source_basis.len()
    {
        return Err(IntegralTopologyError::KernelCoordinateInvariantLost { field });
    }
    Ok(())
}

fn planned_kernel_scalar_operations(
    chain_extent: usize,
    incoming_cols: usize,
) -> Result<u128, IntegralTopologyError> {
    u128::try_from(chain_extent)
        .ok()
        .and_then(|extent| extent.checked_mul(extent))
        .and_then(|terms| {
            u128::try_from(incoming_cols)
                .ok()
                .and_then(|cols| terms.checked_mul(cols))
        })
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "kernel-coordinate scalar operations",
        })
}

fn retained_kernel_transport_entries(
    outgoing: &TerminalRelativeBoundaryMatrix,
    incoming: &TerminalRelativeBoundaryMatrix,
    smith: &VerifiedSmithNormalForm,
    output_entries: usize,
) -> Result<usize, IntegralTopologyError> {
    [
        outgoing.matrix.entries.len(),
        outgoing.source_basis.len(),
        outgoing.target_basis.len(),
        incoming.matrix.entries.len(),
        incoming.source_basis.len(),
        incoming.target_basis.len(),
        smith.source.entries.len(),
        smith.witness.diagonal.entries.len(),
        smith.witness.left.entries.len(),
        smith.witness.left_inverse.entries.len(),
        smith.witness.right.entries.len(),
        smith.witness.right_inverse.entries.len(),
        smith.invariant_factors.len(),
        output_entries,
    ]
    .into_iter()
    .try_fold(0_usize, |total, entries| total.checked_add(entries))
    .ok_or(IntegralTopologyError::WorkPlanOverflow {
        phase: "kernel-coordinate retained entries",
    })
}

#[allow(clippy::too_many_arguments)]
fn poll_kernel_transport(
    checkpoint: &mut impl FnMut(&'static str) -> bool,
    phase: &'static str,
    completed_binding_items: usize,
    planned_binding_items: usize,
    completed_scalar_operations: u128,
    planned_scalar_operations: u128,
) -> Result<(), IntegralTopologyError> {
    if checkpoint(phase) {
        Ok(())
    } else {
        Err(IntegralTopologyError::KernelCoordinateCancelled {
            phase,
            completed_binding_items,
            planned_binding_items,
            completed_scalar_operations,
            planned_scalar_operations,
        })
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
    /// Incoming boundary transformed into outgoing Smith coordinates,
    /// `V^-1 * B`.
    KernelCoordinateIncoming,
}

/// Elementary phase whose checked constructive Smith arithmetic refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SmithConstructionStage {
    /// Euclidean elimination below the active pivot.
    RowReduction,
    /// Euclidean elimination right of the active pivot.
    ColumnReduction,
    /// Trailing-block mixing needed to repair invariant-factor divisibility.
    DivisibilityRepair,
    /// Conversion of a negative final pivot to its positive canonical form.
    SignNormalization,
}

/// Exact generator-lift calculation or verification phase.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HomologyGeneratorStage {
    /// `V[:, rank..] * P^-1` lift into the original `C_k` basis.
    OriginalChainLift,
    /// Exact outgoing-boundary cycle check `A_k * G = 0`.
    CycleVerification,
    /// Exact finite-presentation relation `d_i G_i = A_(k+1) Q_i`.
    BoundaryOrderVerification,
}

/// Algebraic role of one column in the retained quotient-presentation basis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HomologyGeneratorKind {
    /// Nontrivial finite cyclic summand with exact positive order.
    Torsion {
        /// Exact cyclic order.
        order: i128,
    },
    /// Infinite free summand.
    Free,
}

/// Exact free-cocycle construction or verification phase.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HomologyCocycleStage {
    /// `V_A^-T[:, rank..] * U_L^T[:, boundary_rank..]` lift into the
    /// pair-bound dual `C^k` basis.
    OriginalCochainLift,
    /// Exact cocycle check `A_(k+1)^T * Omega = 0`.
    CocycleVerification,
    /// Exact free-period check `G^T * Omega = [0 | I]^T`.
    PeriodDualVerification,
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
    /// Exact incidence bound to an admitted pair, phase, component, degree,
    /// and canonical quotient bases. No homology or physical R3 conclusion.
    TerminalRelativeIncidenceOnly,
    /// Exact incoming image in one retained Smith witness's kernel coordinates.
    /// The coordinates are witness-dependent and are not canonical homology.
    TerminalRelativeKernelCoordinatesOnly,
    /// Exact invariant-factor decomposition of one admitted phase-local
    /// quotient chain homology group. No generator, period, naturality,
    /// embedding, or physical-R3 conclusion follows.
    TerminalRelativeChainHomologyOnly,
    /// Original-chain presentation generators with exact cycle and finite-
    /// order boundary witnesses. No period, naturality, embedding, or
    /// physical-R3 conclusion follows.
    TerminalRelativeHomologyGeneratorsOnly,
    /// Exact integer cocycles dual to the witness-relative free homology
    /// generators. No torsion-cohomology, metric, continuum, or physical
    /// period conclusion follows.
    TerminalRelativeFreeCocyclePeriodsOnly,
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

/// Deterministically constructed Smith witness after independent exact
/// verification.
///
/// Construction counters remain attached to this receipt so consumers that
/// retain constructive provenance can distinguish reduction work from the
/// independently replayed verification work.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConstructedSmithNormalForm {
    verified: VerifiedSmithNormalForm,
    elementary_operations: u128,
    entry_steps: u128,
}

impl ConstructedSmithNormalForm {
    /// Independently verified source, canonical diagonal, and complete
    /// unimodular witness.
    #[must_use]
    pub const fn verified(&self) -> &VerifiedSmithNormalForm {
        &self.verified
    }

    /// Exact elementary operations completed by deterministic construction.
    #[must_use]
    pub const fn elementary_operations(&self) -> u128 {
        self.elementary_operations
    }

    /// Exact admitted source/diagonal inspections and destination updates
    /// completed by construction.
    #[must_use]
    pub const fn entry_steps(&self) -> u128 {
        self.entry_steps
    }

    /// Consume the construction receipt while preserving its independently
    /// verified Smith authority for a downstream algebraic check.
    #[must_use]
    pub fn into_verified(self) -> VerifiedSmithNormalForm {
        self.verified
    }

    /// Constructive Smith form is still abstract algebra, not physical
    /// terminal-relative topology authority.
    #[must_use]
    pub const fn applicability(&self) -> TopologyApplicability {
        TopologyApplicability::AbstractAlgebraOnly
    }
}

/// Exact invariant-factor decomposition of one admitted terminal-relative
/// cellular chain homology group.
///
/// The result owns the complete adjacent-boundary/kernel authority and the
/// complete independently verified Smith authority for its lower image. It
/// publishes only the abstract quotient invariants, not generators or periods.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VerifiedTerminalRelativeHomology {
    transport: VerifiedTerminalRelativeKernelTransport,
    image_smith: VerifiedSmithNormalForm,
    torsion_start: usize,
    binding_items: usize,
    retained_entries: usize,
}

/// Exact witness-relative generators for one admitted terminal-relative
/// cellular homology group.
///
/// Generator columns are ordered as nontrivial torsion summands followed by
/// free summands. Factors equal to one are omitted. Rows of
/// [`Self::cycle_representatives`] use [`Self::original_chain_basis`], while
/// rows of [`Self::torsion_bounding_chains`] use
/// [`Self::bounding_chain_basis`]. The complete homology and Smith authorities
/// remain attached, so the exact order of each torsion column does not rest on
/// the retained boundary equation alone.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VerifiedTerminalRelativeHomologyGenerators {
    homology: VerifiedTerminalRelativeHomology,
    cycle_representatives: ExactIntegerMatrix,
    torsion_bounding_chains: ExactIntegerMatrix,
    torsion_count: usize,
    work_items: u128,
    scalar_operations: u128,
    retained_entries: usize,
}

impl VerifiedTerminalRelativeHomologyGenerators {
    /// Admitted terminal-relative pair identity.
    #[must_use]
    pub const fn pair_id(&self) -> TerminalRelativePairId {
        self.homology.pair_id()
    }

    /// Admitted phase identity.
    #[must_use]
    pub const fn phase(&self) -> &PhaseId {
        self.homology.phase()
    }

    /// Admitted phase-owned conductor component.
    #[must_use]
    pub const fn component(&self) -> &ConductorComponentId {
        self.homology.component()
    }

    /// Homological degree `k` represented by every generator column.
    #[must_use]
    pub const fn degree(&self) -> u8 {
        self.homology.degree()
    }

    /// Complete quotient-homology authority from which these columns were
    /// lifted.
    #[must_use]
    pub const fn homology(&self) -> &VerifiedTerminalRelativeHomology {
        &self.homology
    }

    /// Exact cycle representatives in the original pair-bound `C_k` basis.
    #[must_use]
    pub const fn cycle_representatives(&self) -> &ExactIntegerMatrix {
        &self.cycle_representatives
    }

    /// Exact incoming-chain witnesses whose boundaries are the corresponding
    /// torsion orders times the first [`Self::torsion_generator_count`]
    /// representative columns.
    #[must_use]
    pub const fn torsion_bounding_chains(&self) -> &ExactIntegerMatrix {
        &self.torsion_bounding_chains
    }

    /// Canonically ordered pair-bound `C_k` cells corresponding to generator
    /// rows.
    #[must_use]
    pub fn original_chain_basis(&self) -> &[CellRef] {
        self.homology.transport().outgoing_boundary().source_basis()
    }

    /// Canonically ordered pair-bound `C_(k+1)` cells corresponding to
    /// torsion-filling rows.
    #[must_use]
    pub fn bounding_chain_basis(&self) -> &[CellRef] {
        self.homology.transport().incoming_boundary().source_basis()
    }

    /// Number of retained nontrivial torsion and free generators.
    #[must_use]
    pub const fn generator_count(&self) -> usize {
        self.cycle_representatives.cols
    }

    /// Number of retained nontrivial torsion generators.
    #[must_use]
    pub const fn torsion_generator_count(&self) -> usize {
        self.torsion_count
    }

    /// Number of retained free generators.
    #[must_use]
    pub const fn free_generator_count(&self) -> usize {
        self.generator_count() - self.torsion_count
    }

    /// Algebraic role of one retained generator column.
    #[must_use]
    pub fn generator_kind(&self, column: usize) -> Option<HomologyGeneratorKind> {
        if column >= self.generator_count() {
            return None;
        }
        if column < self.torsion_count {
            return Some(HomologyGeneratorKind::Torsion {
                order: self.homology.image_smith.invariant_factors
                    [self.homology.torsion_start + column],
            });
        }
        Some(HomologyGeneratorKind::Free)
    }

    /// Deterministic output and equation-verification work units completed.
    #[must_use]
    pub const fn work_items(&self) -> u128 {
        self.work_items
    }

    /// Checked lift, cycle, and torsion-order scalar terms completed.
    #[must_use]
    pub const fn scalar_operations(&self) -> u128 {
        self.scalar_operations
    }

    /// Retained integer/cell entries across homology authority and both
    /// generator-witness matrices.
    #[must_use]
    pub const fn retained_entries(&self) -> usize {
        self.retained_entries
    }

    /// Consume the generator receipt while preserving its complete quotient-
    /// homology authority.
    #[must_use]
    pub fn into_homology(self) -> VerifiedTerminalRelativeHomology {
        self.homology
    }

    /// Exact witness-relative cellular generators only.
    #[must_use]
    pub const fn applicability(&self) -> TopologyApplicability {
        TopologyApplicability::TerminalRelativeHomologyGeneratorsOnly
    }
}

/// Exact free integral cocycles and their periods on the retained homology
/// generators.
///
/// Each column of [`Self::cocycle_representatives`] is expressed in the dual
/// of [`Self::cochain_basis`]; evaluating a chain column `x` means
/// `Omega^T x`. Rows of [`Self::period_pairing`] follow the generator order
/// (torsion, then free), and columns follow the free cocycle order. The checked
/// matrix is exactly `[0 | I]^T`: integer-valued cocycles kill finite torsion
/// and pair primitively with the selected free generators.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VerifiedTerminalRelativeFreeCocyclePeriods {
    generators: VerifiedTerminalRelativeHomologyGenerators,
    cocycle_representatives: ExactIntegerMatrix,
    period_pairing: ExactIntegerMatrix,
    work_items: u128,
    scalar_operations: u128,
    retained_entries: usize,
}

impl VerifiedTerminalRelativeFreeCocyclePeriods {
    /// Complete generator authority to which these cocycles are dual.
    #[must_use]
    pub const fn generators(&self) -> &VerifiedTerminalRelativeHomologyGenerators {
        &self.generators
    }

    /// Admitted terminal-relative pair identity.
    #[must_use]
    pub const fn pair_id(&self) -> TerminalRelativePairId {
        self.generators.pair_id()
    }

    /// Admitted phase identity.
    #[must_use]
    pub const fn phase(&self) -> &PhaseId {
        self.generators.phase()
    }

    /// Admitted phase-owned conductor component.
    #[must_use]
    pub const fn component(&self) -> &ConductorComponentId {
        self.generators.component()
    }

    /// Homological/cohomological degree `k` represented by the pairing.
    #[must_use]
    pub const fn degree(&self) -> u8 {
        self.generators.degree()
    }

    /// Pair-bound ordered `C_k` cells whose induced dual basis supplies the
    /// cocycle rows.
    #[must_use]
    pub fn cochain_basis(&self) -> &[CellRef] {
        self.generators.original_chain_basis()
    }

    /// Integer cocycle columns in the pair-bound dual `C^k` basis.
    #[must_use]
    pub const fn cocycle_representatives(&self) -> &ExactIntegerMatrix {
        &self.cocycle_representatives
    }

    /// Exact generator-by-cocycle integer period matrix `[0 | I]^T`.
    #[must_use]
    pub const fn period_pairing(&self) -> &ExactIntegerMatrix {
        &self.period_pairing
    }

    /// Rank of the verified free period lattice.
    #[must_use]
    pub const fn free_cocycle_count(&self) -> usize {
        self.cocycle_representatives.cols
    }

    /// Generator column paired to one free cocycle, or `None` outside the
    /// retained free-cocycle range.
    #[must_use]
    pub fn paired_generator_column(&self, cocycle: usize) -> Option<usize> {
        if cocycle >= self.free_cocycle_count() {
            return None;
        }
        Some(self.generators.torsion_count + cocycle)
    }

    /// Deterministic output and equation-verification work units completed.
    #[must_use]
    pub const fn work_items(&self) -> u128 {
        self.work_items
    }

    /// Checked cocycle-lift, coboundary, and period-pairing terms completed.
    #[must_use]
    pub const fn scalar_operations(&self) -> u128 {
        self.scalar_operations
    }

    /// Retained integer/cell entries across generators, cocycles, and periods.
    #[must_use]
    pub const fn retained_entries(&self) -> usize {
        self.retained_entries
    }

    /// Consume the period receipt while preserving its complete generator
    /// authority.
    #[must_use]
    pub fn into_generators(self) -> VerifiedTerminalRelativeHomologyGenerators {
        self.generators
    }

    /// Exact witness-relative free integral period duals only.
    #[must_use]
    pub const fn applicability(&self) -> TopologyApplicability {
        TopologyApplicability::TerminalRelativeFreeCocyclePeriodsOnly
    }
}

impl VerifiedTerminalRelativeHomology {
    /// Admitted terminal-relative pair identity.
    #[must_use]
    pub const fn pair_id(&self) -> TerminalRelativePairId {
        self.transport.pair_id()
    }

    /// Admitted phase identity.
    #[must_use]
    pub const fn phase(&self) -> &PhaseId {
        self.transport.phase()
    }

    /// Admitted phase-owned conductor component.
    #[must_use]
    pub const fn component(&self) -> &ConductorComponentId {
        self.transport.component()
    }

    /// Homological degree `k` of this quotient.
    #[must_use]
    pub const fn degree(&self) -> u8 {
        self.transport.degree()
    }

    /// Complete proof that the incoming image lies in the outgoing kernel.
    #[must_use]
    pub const fn transport(&self) -> &VerifiedTerminalRelativeKernelTransport {
        &self.transport
    }

    /// Complete Smith authority bound byte-for-byte to the lower image.
    #[must_use]
    pub const fn image_smith(&self) -> &VerifiedSmithNormalForm {
        &self.image_smith
    }

    /// Rank of `ker(A_k)` before quotienting by incoming boundaries.
    #[must_use]
    pub const fn cycle_rank(&self) -> usize {
        self.transport.kernel_dimension()
    }

    /// Rank of `im(A_(k+1))` inside the verified kernel coordinates.
    #[must_use]
    pub const fn boundary_rank(&self) -> usize {
        self.image_smith.rank()
    }

    /// Rank of the free summand `Z^(cycle_rank - boundary_rank)`.
    #[must_use]
    pub const fn free_rank(&self) -> usize {
        self.cycle_rank() - self.boundary_rank()
    }

    /// All positive presentation factors, including factors equal to one.
    #[must_use]
    pub fn presentation_invariant_factors(&self) -> &[i128] {
        self.image_smith.invariant_factors()
    }

    /// Nontrivial finite cyclic summands. Factors equal to one are omitted.
    #[must_use]
    pub fn torsion_invariant_factors(&self) -> &[i128] {
        &self.image_smith.invariant_factors()[self.torsion_start..]
    }

    /// Exact source/factor inspections completed before publication.
    #[must_use]
    pub const fn binding_items(&self) -> usize {
        self.binding_items
    }

    /// Retained integer/cell entries across both complete authorities.
    #[must_use]
    pub const fn retained_entries(&self) -> usize {
        self.retained_entries
    }

    /// Exact phase-local cellular quotient homology only.
    #[must_use]
    pub const fn applicability(&self) -> TopologyApplicability {
        TopologyApplicability::TerminalRelativeChainHomologyOnly
    }
}

/// Bind a verified Smith form to one exact lower kernel image without injected
/// cancellation and publish its free/torsion invariant decomposition.
#[allow(clippy::large_types_passed_by_value)]
pub fn verify_terminal_relative_homology(
    transport: VerifiedTerminalRelativeKernelTransport,
    image_smith: VerifiedSmithNormalForm,
    budget: HomologyDecompositionBudget,
) -> Result<VerifiedTerminalRelativeHomology, IntegralTopologyError> {
    verify_terminal_relative_homology_with_checkpoint(transport, image_smith, budget, &mut |_| true)
}

/// Bind a lower-image Smith authority with bounded cancellation polling.
///
/// For lower image `L` with `q` rows and verified invariant factors
/// `d_1 | ... | d_s`, this proves only
/// `H_k = Z^(q-s) + direct_sum_i Z/d_i`, omitting factors `d_i = 1` from the
/// reported torsion slice. The result retains the authorities needed to replay
/// every byte of this statement.
#[allow(clippy::large_types_passed_by_value)]
#[allow(clippy::too_many_lines)]
pub fn verify_terminal_relative_homology_with_checkpoint(
    transport: VerifiedTerminalRelativeKernelTransport,
    image_smith: VerifiedSmithNormalForm,
    budget: HomologyDecompositionBudget,
    checkpoint: &mut impl FnMut(&'static str) -> bool,
) -> Result<VerifiedTerminalRelativeHomology, IntegralTopologyError> {
    let lower = transport.kernel_image();
    let source = image_smith.source();
    if source.rows != lower.rows || source.cols != lower.cols {
        return Err(IntegralTopologyError::HomologySmithSourceShapeMismatch {
            expected_rows: lower.rows,
            expected_cols: lower.cols,
            actual_rows: source.rows,
            actual_cols: source.cols,
        });
    }
    if image_smith.rank > lower.rows || image_smith.rank != image_smith.invariant_factors.len() {
        return Err(IntegralTopologyError::HomologyDecompositionInvariantLost {
            field: "lower-image Smith rank",
        });
    }
    let binding_items = lower
        .entries
        .len()
        .checked_add(image_smith.invariant_factors.len())
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative homology binding items",
        })?;
    if binding_items > budget.max_binding_items {
        return Err(IntegralTopologyError::HomologyBindingBudgetExceeded {
            requested: binding_items,
            max: budget.max_binding_items,
        });
    }
    let retained_entries = retained_terminal_relative_homology_entries(&transport, &image_smith)?;
    if retained_entries > budget.max_retained_entries {
        return Err(IntegralTopologyError::RetainedEntryBudgetExceeded {
            requested: retained_entries,
            max: budget.max_retained_entries,
        });
    }

    let mut completed = 0_usize;
    poll_homology_decomposition(
        checkpoint,
        "terminal-relative homology preflight",
        completed,
        binding_items,
    )?;
    for (index, (expected, actual)) in lower.entries.iter().zip(&source.entries).enumerate() {
        poll_homology_decomposition(
            checkpoint,
            "terminal-relative homology source binding",
            completed,
            binding_items,
        )?;
        if expected != actual {
            return Err(IntegralTopologyError::HomologySmithSourceEntryMismatch {
                row: index / lower.cols.max(1),
                col: index % lower.cols.max(1),
                expected: *expected,
                actual: *actual,
            });
        }
        completed += 1;
    }

    let mut torsion_start = image_smith.invariant_factors.len();
    for (index, factor) in image_smith.invariant_factors.iter().copied().enumerate() {
        poll_homology_decomposition(
            checkpoint,
            "terminal-relative homology invariant factors",
            completed,
            binding_items,
        )?;
        if factor <= 0 {
            return Err(IntegralTopologyError::HomologyDecompositionInvariantLost {
                field: "nonpositive verified invariant factor",
            });
        }
        if factor > 1 && torsion_start == image_smith.invariant_factors.len() {
            torsion_start = index;
        }
        completed += 1;
    }
    poll_homology_decomposition(
        checkpoint,
        "terminal-relative homology finalize",
        completed,
        binding_items,
    )?;
    debug_assert_eq!(completed, binding_items);

    Ok(VerifiedTerminalRelativeHomology {
        transport,
        image_smith,
        torsion_start,
        binding_items: completed,
        retained_entries,
    })
}

fn retained_terminal_relative_homology_entries(
    transport: &VerifiedTerminalRelativeKernelTransport,
    image_smith: &VerifiedSmithNormalForm,
) -> Result<usize, IntegralTopologyError> {
    let transport_entries = retained_kernel_transport_entries(
        &transport.outgoing,
        &transport.incoming,
        &transport.outgoing_smith,
        transport.kernel_image.entries.len(),
    )?;
    [
        image_smith.source.entries.len(),
        image_smith.witness.diagonal.entries.len(),
        image_smith.witness.left.entries.len(),
        image_smith.witness.left_inverse.entries.len(),
        image_smith.witness.right.entries.len(),
        image_smith.witness.right_inverse.entries.len(),
        image_smith.invariant_factors.len(),
    ]
    .into_iter()
    .try_fold(transport_entries, |total, entries| {
        total.checked_add(entries)
    })
    .ok_or(IntegralTopologyError::WorkPlanOverflow {
        phase: "terminal-relative homology retained entries",
    })
}

fn poll_homology_decomposition(
    checkpoint: &mut impl FnMut(&'static str) -> bool,
    phase: &'static str,
    completed_binding_items: usize,
    planned_binding_items: usize,
) -> Result<(), IntegralTopologyError> {
    if checkpoint(phase) {
        Ok(())
    } else {
        Err(IntegralTopologyError::HomologyDecompositionCancelled {
            phase,
            completed_binding_items,
            planned_binding_items,
        })
    }
}

/// Lift every nontrivial torsion and free presentation column into the
/// original pair-bound chain basis without injected cancellation.
#[allow(clippy::large_types_passed_by_value)]
pub fn verify_terminal_relative_homology_generators(
    homology: VerifiedTerminalRelativeHomology,
    budget: HomologyGeneratorBudget,
) -> Result<VerifiedTerminalRelativeHomologyGenerators, IntegralTopologyError> {
    verify_terminal_relative_homology_generators_with_checkpoint(homology, budget, &mut |_| true)
}

/// Lift and independently replay the chain-level generator equations under
/// bounded cancellation polling.
///
/// For `U_A A_k V_A = D_A` of rank `r`, lower kernel image `L`, and
/// `U_L L V_L = D_L`, this constructs the selected columns of
/// `G = V_A[:, r..] U_L^-1`. Unit-factor columns are omitted. It also retains
/// `H = V_L[:, torsion_start..rank(L)]` and verifies exactly that
/// `A_k G = 0` and `A_(k+1) H_j = d_j G_j`. Exact torsion order additionally
/// follows from the retained complete Smith authority; the boundary equation
/// alone is not treated as sufficient.
#[allow(clippy::large_types_passed_by_value)]
#[allow(clippy::too_many_lines)]
pub fn verify_terminal_relative_homology_generators_with_checkpoint(
    homology: VerifiedTerminalRelativeHomology,
    budget: HomologyGeneratorBudget,
    checkpoint: &mut impl FnMut(&'static str) -> bool,
) -> Result<VerifiedTerminalRelativeHomologyGenerators, IntegralTopologyError> {
    let transport = homology.transport();
    let outgoing = transport.outgoing_boundary().matrix();
    let incoming = transport.incoming_boundary().matrix();
    let outgoing_smith = transport.outgoing_smith();
    let image_smith = homology.image_smith();

    let outgoing_rows = outgoing.rows;
    let chain_extent = outgoing.cols;
    let incoming_cols = incoming.cols;
    let outgoing_rank = outgoing_smith.rank;
    let cycle_rank = homology.cycle_rank();
    let boundary_rank = homology.boundary_rank();
    let torsion_start = homology.torsion_start;

    let internal_shapes_hold = incoming.rows == chain_extent
        && outgoing_rank <= chain_extent
        && cycle_rank == chain_extent - outgoing_rank
        && image_smith.source.rows == cycle_rank
        && image_smith.source.cols == incoming_cols
        && image_smith.left_inverse().rows == cycle_rank
        && image_smith.left_inverse().cols == cycle_rank
        && image_smith.right_transform().rows == incoming_cols
        && image_smith.right_transform().cols == incoming_cols
        && outgoing_smith.right_transform().rows == chain_extent
        && outgoing_smith.right_transform().cols == chain_extent
        && boundary_rank == image_smith.invariant_factors.len()
        && torsion_start <= boundary_rank
        && boundary_rank <= cycle_rank;
    if !internal_shapes_hold {
        return Err(IntegralTopologyError::HomologyGeneratorInvariantLost {
            field: "retained generator-lift shapes and ranks",
        });
    }
    if image_smith.invariant_factors[..torsion_start]
        .iter()
        .any(|factor| *factor != 1)
        || image_smith.invariant_factors[torsion_start..]
            .iter()
            .any(|factor| *factor <= 1)
    {
        return Err(IntegralTopologyError::HomologyGeneratorInvariantLost {
            field: "unit and nontrivial presentation-factor split",
        });
    }

    let torsion_count = boundary_rank - torsion_start;
    let generator_count = cycle_rank - torsion_start;
    let cycle_entries = chain_extent.checked_mul(generator_count).ok_or(
        IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative generator cycle output entries",
        },
    )?;
    let bounding_entries = incoming_cols.checked_mul(torsion_count).ok_or(
        IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative generator bounding output entries",
        },
    )?;
    let output_entries = cycle_entries.checked_add(bounding_entries).ok_or(
        IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative generator output entries",
        },
    )?;
    if output_entries > budget.max_output_entries {
        return Err(
            IntegralTopologyError::HomologyGeneratorOutputBudgetExceeded {
                requested: output_entries,
                max: budget.max_output_entries,
            },
        );
    }
    let retained_entries = homology
        .retained_entries
        .checked_add(output_entries)
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative generator retained entries",
        })?;
    if retained_entries > budget.max_retained_entries {
        return Err(IntegralTopologyError::RetainedEntryBudgetExceeded {
            requested: retained_entries,
            max: budget.max_retained_entries,
        });
    }
    let scalar_operations = planned_homology_generator_scalar_operations(
        outgoing_rows,
        chain_extent,
        incoming_cols,
        cycle_rank,
        generator_count,
        torsion_count,
    )?;
    if scalar_operations > budget.max_scalar_operations {
        return Err(IntegralTopologyError::ScalarWorkBudgetExceeded {
            requested: scalar_operations,
            max: budget.max_scalar_operations,
        });
    }
    let work_items = planned_homology_generator_work_items(
        outgoing_rows,
        chain_extent,
        incoming_cols,
        generator_count,
        torsion_count,
    )?;

    let mut completed_work = 0_u128;
    let mut completed_scalar = 0_u128;
    poll_homology_generators(
        checkpoint,
        "terminal-relative generator preflight",
        completed_work,
        work_items,
        completed_scalar,
        scalar_operations,
    )?;
    poll_homology_generators(
        checkpoint,
        "terminal-relative generator output allocation",
        completed_work,
        work_items,
        completed_scalar,
        scalar_operations,
    )?;
    let mut cycle_values = allocate_zeroed(cycle_entries, "homology cycle representatives")?;
    let mut bounding_values =
        allocate_zeroed(bounding_entries, "homology torsion bounding chains")?;

    let outgoing_right = outgoing_smith.right_transform();
    let image_left_inverse = image_smith.left_inverse();
    for row in 0..chain_extent {
        for generator in 0..generator_count {
            poll_homology_generators(
                checkpoint,
                "terminal-relative generator original-chain lift",
                completed_work,
                work_items,
                completed_scalar,
                scalar_operations,
            )?;
            let presentation_column = torsion_start + generator;
            let mut sum = 0_i128;
            for term in 0..cycle_rank {
                sum = checked_homology_generator_accumulate(
                    sum,
                    outgoing_right.entry(row, outgoing_rank + term),
                    image_left_inverse.entry(term, presentation_column),
                    HomologyGeneratorStage::OriginalChainLift,
                    row,
                    generator,
                    term,
                    &mut completed_scalar,
                )?;
            }
            cycle_values[row * generator_count + generator] = sum;
            completed_work =
                completed_work
                    .checked_add(1)
                    .ok_or(IntegralTopologyError::WorkPlanOverflow {
                        phase: "completed terminal-relative generator work",
                    })?;
        }
    }
    let cycle_representatives = ExactIntegerMatrix {
        rows: chain_extent,
        cols: generator_count,
        entries: cycle_values,
    };

    let image_right = image_smith.right_transform();
    for row in 0..incoming_cols {
        for torsion in 0..torsion_count {
            poll_homology_generators(
                checkpoint,
                "terminal-relative generator torsion fillings",
                completed_work,
                work_items,
                completed_scalar,
                scalar_operations,
            )?;
            bounding_values[row * torsion_count + torsion] =
                image_right.entry(row, torsion_start + torsion);
            completed_work =
                completed_work
                    .checked_add(1)
                    .ok_or(IntegralTopologyError::WorkPlanOverflow {
                        phase: "completed terminal-relative generator work",
                    })?;
        }
    }
    let torsion_bounding_chains = ExactIntegerMatrix {
        rows: incoming_cols,
        cols: torsion_count,
        entries: bounding_values,
    };

    for row in 0..outgoing_rows {
        for generator in 0..generator_count {
            poll_homology_generators(
                checkpoint,
                "terminal-relative generator cycle verification",
                completed_work,
                work_items,
                completed_scalar,
                scalar_operations,
            )?;
            let actual = checked_homology_generator_dot(
                outgoing,
                row,
                &cycle_representatives,
                generator,
                HomologyGeneratorStage::CycleVerification,
                &mut completed_scalar,
            )?;
            if actual != 0 {
                return Err(
                    IntegralTopologyError::HomologyGeneratorVerificationMismatch {
                        stage: HomologyGeneratorStage::CycleVerification,
                        row,
                        col: generator,
                        expected: 0,
                        actual,
                    },
                );
            }
            completed_work =
                completed_work
                    .checked_add(1)
                    .ok_or(IntegralTopologyError::WorkPlanOverflow {
                        phase: "completed terminal-relative generator work",
                    })?;
        }
    }

    for row in 0..chain_extent {
        for torsion in 0..torsion_count {
            poll_homology_generators(
                checkpoint,
                "terminal-relative generator boundary-order verification",
                completed_work,
                work_items,
                completed_scalar,
                scalar_operations,
            )?;
            let actual = checked_homology_generator_dot(
                incoming,
                row,
                &torsion_bounding_chains,
                torsion,
                HomologyGeneratorStage::BoundaryOrderVerification,
                &mut completed_scalar,
            )?;
            let order = image_smith.invariant_factors[torsion_start + torsion];
            let expected = order
                .checked_mul(cycle_representatives.entry(row, torsion))
                .ok_or(IntegralTopologyError::HomologyGeneratorArithmeticOverflow {
                    stage: HomologyGeneratorStage::BoundaryOrderVerification,
                    row,
                    col: torsion,
                    term: incoming_cols,
                })?;
            completed_scalar =
                completed_scalar
                    .checked_add(1)
                    .ok_or(IntegralTopologyError::WorkPlanOverflow {
                        phase: "completed terminal-relative generator scalar operations",
                    })?;
            if actual != expected {
                return Err(
                    IntegralTopologyError::HomologyGeneratorVerificationMismatch {
                        stage: HomologyGeneratorStage::BoundaryOrderVerification,
                        row,
                        col: torsion,
                        expected,
                        actual,
                    },
                );
            }
            completed_work =
                completed_work
                    .checked_add(1)
                    .ok_or(IntegralTopologyError::WorkPlanOverflow {
                        phase: "completed terminal-relative generator work",
                    })?;
        }
    }

    poll_homology_generators(
        checkpoint,
        "terminal-relative generator finalize",
        completed_work,
        work_items,
        completed_scalar,
        scalar_operations,
    )?;
    debug_assert_eq!(completed_work, work_items);
    debug_assert_eq!(completed_scalar, scalar_operations);

    Ok(VerifiedTerminalRelativeHomologyGenerators {
        homology,
        cycle_representatives,
        torsion_bounding_chains,
        torsion_count,
        work_items: completed_work,
        scalar_operations: completed_scalar,
        retained_entries,
    })
}

fn planned_homology_generator_scalar_operations(
    outgoing_rows: usize,
    chain_extent: usize,
    incoming_cols: usize,
    cycle_rank: usize,
    generator_count: usize,
    torsion_count: usize,
) -> Result<u128, IntegralTopologyError> {
    let lift = checked_homology_generator_plan_product(
        &[chain_extent, cycle_rank, generator_count],
        "terminal-relative generator lift scalar operations",
    )?;
    let cycle = checked_homology_generator_plan_product(
        &[outgoing_rows, chain_extent, generator_count],
        "terminal-relative generator cycle scalar operations",
    )?;
    let boundary = checked_homology_generator_plan_product(
        &[chain_extent, incoming_cols, torsion_count],
        "terminal-relative generator boundary scalar operations",
    )?;
    let scaling = checked_homology_generator_plan_product(
        &[chain_extent, torsion_count],
        "terminal-relative generator scaling operations",
    )?;
    [lift, cycle, boundary, scaling]
        .into_iter()
        .try_fold(0_u128, u128::checked_add)
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative generator scalar operations",
        })
}

fn planned_homology_generator_work_items(
    outgoing_rows: usize,
    chain_extent: usize,
    incoming_cols: usize,
    generator_count: usize,
    torsion_count: usize,
) -> Result<u128, IntegralTopologyError> {
    let cycle_output = checked_homology_generator_plan_product(
        &[chain_extent, generator_count],
        "terminal-relative generator cycle output work",
    )?;
    let bounding_output = checked_homology_generator_plan_product(
        &[incoming_cols, torsion_count],
        "terminal-relative generator bounding output work",
    )?;
    let cycle_checks = checked_homology_generator_plan_product(
        &[outgoing_rows, generator_count],
        "terminal-relative generator cycle check work",
    )?;
    let boundary_checks = checked_homology_generator_plan_product(
        &[chain_extent, torsion_count],
        "terminal-relative generator boundary check work",
    )?;
    [cycle_output, bounding_output, cycle_checks, boundary_checks]
        .into_iter()
        .try_fold(0_u128, u128::checked_add)
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative generator work items",
        })
}

fn checked_homology_generator_plan_product(
    factors: &[usize],
    phase: &'static str,
) -> Result<u128, IntegralTopologyError> {
    factors
        .iter()
        .copied()
        .try_fold(1_u128, |product, factor| {
            u128::try_from(factor)
                .ok()
                .and_then(|factor| product.checked_mul(factor))
        })
        .ok_or(IntegralTopologyError::WorkPlanOverflow { phase })
}

#[allow(clippy::too_many_arguments)]
fn checked_homology_generator_accumulate(
    sum: i128,
    left: i128,
    right: i128,
    stage: HomologyGeneratorStage,
    row: usize,
    col: usize,
    term: usize,
    completed_scalar: &mut u128,
) -> Result<i128, IntegralTopologyError> {
    let product = left.checked_mul(right).ok_or(
        IntegralTopologyError::HomologyGeneratorArithmeticOverflow {
            stage,
            row,
            col,
            term,
        },
    )?;
    let next = sum.checked_add(product).ok_or(
        IntegralTopologyError::HomologyGeneratorArithmeticOverflow {
            stage,
            row,
            col,
            term,
        },
    )?;
    *completed_scalar =
        completed_scalar
            .checked_add(1)
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "completed terminal-relative generator scalar operations",
            })?;
    Ok(next)
}

fn checked_homology_generator_dot(
    left: &ExactIntegerMatrix,
    row: usize,
    right: &ExactIntegerMatrix,
    col: usize,
    stage: HomologyGeneratorStage,
    completed_scalar: &mut u128,
) -> Result<i128, IntegralTopologyError> {
    if left.cols != right.rows {
        return Err(IntegralTopologyError::HomologyGeneratorInvariantLost {
            field: "generator verification inner matrix extent",
        });
    }
    let mut sum = 0_i128;
    for term in 0..left.cols {
        sum = checked_homology_generator_accumulate(
            sum,
            left.entry(row, term),
            right.entry(term, col),
            stage,
            row,
            col,
            term,
            completed_scalar,
        )?;
    }
    Ok(sum)
}

fn poll_homology_generators(
    checkpoint: &mut impl FnMut(&'static str) -> bool,
    phase: &'static str,
    completed_work_items: u128,
    planned_work_items: u128,
    completed_scalar_operations: u128,
    planned_scalar_operations: u128,
) -> Result<(), IntegralTopologyError> {
    if checkpoint(phase) {
        Ok(())
    } else {
        Err(IntegralTopologyError::HomologyGeneratorCancelled {
            phase,
            completed_work_items,
            planned_work_items,
            completed_scalar_operations,
            planned_scalar_operations,
        })
    }
}

/// Construct the exact integer cocycles dual to the retained free generators
/// without injected cancellation.
#[allow(clippy::large_types_passed_by_value)]
pub fn verify_terminal_relative_free_cocycle_periods(
    generators: VerifiedTerminalRelativeHomologyGenerators,
    budget: HomologyCocycleBudget,
) -> Result<VerifiedTerminalRelativeFreeCocyclePeriods, IntegralTopologyError> {
    verify_terminal_relative_free_cocycle_periods_with_checkpoint(generators, budget, &mut |_| true)
}

/// Construct and independently replay the free integer period dual under
/// bounded cancellation polling.
///
/// For outgoing `U_A A_k V_A = D_A` of rank `r` and lower-image
/// `U_L L V_L = D_L` of rank `s`, the cocycle columns are
/// `Omega = V_A^-T[:, r..] U_L^T[:, s..]`. Publication requires both
/// `A_(k+1)^T Omega = 0` and `G^T Omega = [0 | I]^T` exactly. This proves the
/// free `Hom(H_k, Z)` dual selected by the retained witnesses; it does not
/// compute the universal-coefficient `Ext` sector or torsion linking.
#[allow(clippy::large_types_passed_by_value)]
#[allow(clippy::too_many_lines)]
pub fn verify_terminal_relative_free_cocycle_periods_with_checkpoint(
    generators: VerifiedTerminalRelativeHomologyGenerators,
    budget: HomologyCocycleBudget,
    checkpoint: &mut impl FnMut(&'static str) -> bool,
) -> Result<VerifiedTerminalRelativeFreeCocyclePeriods, IntegralTopologyError> {
    let homology = generators.homology();
    let transport = homology.transport();
    let incoming = transport.incoming_boundary().matrix();
    let outgoing_smith = transport.outgoing_smith();
    let image_smith = homology.image_smith();
    let generator_matrix = generators.cycle_representatives();

    let chain_extent = transport.chain_extent();
    let incoming_cols = incoming.cols;
    let outgoing_rank = transport.outgoing_rank();
    let cycle_rank = homology.cycle_rank();
    let boundary_rank = homology.boundary_rank();
    let torsion_count = generators.torsion_count;
    let generator_count = generators.generator_count();
    let free_count = generators.free_generator_count();

    let internal_shapes_hold = incoming.rows == chain_extent
        && outgoing_rank <= chain_extent
        && cycle_rank == chain_extent - outgoing_rank
        && boundary_rank <= cycle_rank
        && free_count == cycle_rank - boundary_rank
        && generator_count == torsion_count + free_count
        && generator_matrix.rows == chain_extent
        && generator_matrix.cols == generator_count
        && generators.original_chain_basis().len() == chain_extent
        && outgoing_smith.right_inverse().rows == chain_extent
        && outgoing_smith.right_inverse().cols == chain_extent
        && image_smith.left_transform().rows == cycle_rank
        && image_smith.left_transform().cols == cycle_rank;
    if !internal_shapes_hold {
        return Err(IntegralTopologyError::HomologyCocycleInvariantLost {
            field: "retained free-cocycle shapes and ranks",
        });
    }

    let cocycle_entries =
        chain_extent
            .checked_mul(free_count)
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "terminal-relative free-cocycle output entries",
            })?;
    let pairing_entries =
        generator_count
            .checked_mul(free_count)
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "terminal-relative period-pairing output entries",
            })?;
    let output_entries = cocycle_entries.checked_add(pairing_entries).ok_or(
        IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative free-cocycle aggregate output entries",
        },
    )?;
    if output_entries > budget.max_output_entries {
        return Err(IntegralTopologyError::HomologyCocycleOutputBudgetExceeded {
            requested: output_entries,
            max: budget.max_output_entries,
        });
    }
    let retained_entries = generators
        .retained_entries
        .checked_add(output_entries)
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative free-cocycle retained entries",
        })?;
    if retained_entries > budget.max_retained_entries {
        return Err(IntegralTopologyError::RetainedEntryBudgetExceeded {
            requested: retained_entries,
            max: budget.max_retained_entries,
        });
    }
    let scalar_operations = planned_homology_cocycle_scalar_operations(
        chain_extent,
        incoming_cols,
        cycle_rank,
        generator_count,
        free_count,
    )?;
    if scalar_operations > budget.max_scalar_operations {
        return Err(IntegralTopologyError::ScalarWorkBudgetExceeded {
            requested: scalar_operations,
            max: budget.max_scalar_operations,
        });
    }
    let work_items = planned_homology_cocycle_work_items(
        chain_extent,
        incoming_cols,
        generator_count,
        free_count,
    )?;

    let mut completed_work = 0_u128;
    let mut completed_scalar = 0_u128;
    poll_homology_cocycles(
        checkpoint,
        "terminal-relative free-cocycle preflight",
        completed_work,
        work_items,
        completed_scalar,
        scalar_operations,
    )?;
    poll_homology_cocycles(
        checkpoint,
        "terminal-relative free-cocycle output allocation",
        completed_work,
        work_items,
        completed_scalar,
        scalar_operations,
    )?;
    let mut cocycle_values = allocate_zeroed(cocycle_entries, "free cocycle representatives")?;
    let mut pairing_values = allocate_zeroed(pairing_entries, "free period pairing")?;

    let outgoing_right_inverse = outgoing_smith.right_inverse();
    let image_left = image_smith.left_transform();
    for row in 0..chain_extent {
        for cocycle in 0..free_count {
            poll_homology_cocycles(
                checkpoint,
                "terminal-relative original-cocycle lift",
                completed_work,
                work_items,
                completed_scalar,
                scalar_operations,
            )?;
            let presentation_row = boundary_rank + cocycle;
            let mut sum = 0_i128;
            for term in 0..cycle_rank {
                sum = checked_homology_cocycle_accumulate(
                    sum,
                    outgoing_right_inverse.entry(outgoing_rank + term, row),
                    image_left.entry(presentation_row, term),
                    HomologyCocycleStage::OriginalCochainLift,
                    row,
                    cocycle,
                    term,
                    &mut completed_scalar,
                )?;
            }
            cocycle_values[row * free_count + cocycle] = sum;
            completed_work =
                completed_work
                    .checked_add(1)
                    .ok_or(IntegralTopologyError::WorkPlanOverflow {
                        phase: "completed terminal-relative free-cocycle work",
                    })?;
        }
    }
    let cocycle_representatives = ExactIntegerMatrix {
        rows: chain_extent,
        cols: free_count,
        entries: cocycle_values,
    };

    for incoming_column in 0..incoming_cols {
        for cocycle in 0..free_count {
            poll_homology_cocycles(
                checkpoint,
                "terminal-relative cocycle verification",
                completed_work,
                work_items,
                completed_scalar,
                scalar_operations,
            )?;
            let actual = checked_homology_cocycle_transposed_dot(
                &cocycle_representatives,
                cocycle,
                incoming,
                incoming_column,
                HomologyCocycleStage::CocycleVerification,
                &mut completed_scalar,
            )?;
            if actual != 0 {
                return Err(IntegralTopologyError::HomologyCocycleVerificationMismatch {
                    stage: HomologyCocycleStage::CocycleVerification,
                    row: incoming_column,
                    col: cocycle,
                    expected: 0,
                    actual,
                });
            }
            completed_work =
                completed_work
                    .checked_add(1)
                    .ok_or(IntegralTopologyError::WorkPlanOverflow {
                        phase: "completed terminal-relative free-cocycle work",
                    })?;
        }
    }

    for generator in 0..generator_count {
        for cocycle in 0..free_count {
            poll_homology_cocycles(
                checkpoint,
                "terminal-relative period-dual verification",
                completed_work,
                work_items,
                completed_scalar,
                scalar_operations,
            )?;
            let actual = checked_homology_cocycle_transposed_dot(
                &cocycle_representatives,
                cocycle,
                generator_matrix,
                generator,
                HomologyCocycleStage::PeriodDualVerification,
                &mut completed_scalar,
            )?;
            let expected = i128::from(generator == torsion_count + cocycle);
            if actual != expected {
                return Err(IntegralTopologyError::HomologyCocycleVerificationMismatch {
                    stage: HomologyCocycleStage::PeriodDualVerification,
                    row: generator,
                    col: cocycle,
                    expected,
                    actual,
                });
            }
            pairing_values[generator * free_count + cocycle] = actual;
            completed_work =
                completed_work
                    .checked_add(1)
                    .ok_or(IntegralTopologyError::WorkPlanOverflow {
                        phase: "completed terminal-relative free-cocycle work",
                    })?;
        }
    }
    let period_pairing = ExactIntegerMatrix {
        rows: generator_count,
        cols: free_count,
        entries: pairing_values,
    };

    poll_homology_cocycles(
        checkpoint,
        "terminal-relative free-cocycle finalize",
        completed_work,
        work_items,
        completed_scalar,
        scalar_operations,
    )?;
    debug_assert_eq!(completed_work, work_items);
    debug_assert_eq!(completed_scalar, scalar_operations);

    Ok(VerifiedTerminalRelativeFreeCocyclePeriods {
        generators,
        cocycle_representatives,
        period_pairing,
        work_items: completed_work,
        scalar_operations: completed_scalar,
        retained_entries,
    })
}

fn planned_homology_cocycle_scalar_operations(
    chain_extent: usize,
    incoming_cols: usize,
    cycle_rank: usize,
    generator_count: usize,
    free_count: usize,
) -> Result<u128, IntegralTopologyError> {
    let lift = checked_homology_generator_plan_product(
        &[chain_extent, cycle_rank, free_count],
        "terminal-relative free-cocycle lift scalar operations",
    )?;
    let cocycle_checks = checked_homology_generator_plan_product(
        &[incoming_cols, chain_extent, free_count],
        "terminal-relative cocycle-check scalar operations",
    )?;
    let period_checks = checked_homology_generator_plan_product(
        &[generator_count, chain_extent, free_count],
        "terminal-relative period-check scalar operations",
    )?;
    [lift, cocycle_checks, period_checks]
        .into_iter()
        .try_fold(0_u128, u128::checked_add)
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative free-cocycle scalar operations",
        })
}

fn planned_homology_cocycle_work_items(
    chain_extent: usize,
    incoming_cols: usize,
    generator_count: usize,
    free_count: usize,
) -> Result<u128, IntegralTopologyError> {
    let cocycle_output = checked_homology_generator_plan_product(
        &[chain_extent, free_count],
        "terminal-relative free-cocycle output work",
    )?;
    let cocycle_checks = checked_homology_generator_plan_product(
        &[incoming_cols, free_count],
        "terminal-relative cocycle-check work",
    )?;
    let period_checks = checked_homology_generator_plan_product(
        &[generator_count, free_count],
        "terminal-relative period-check work",
    )?;
    [cocycle_output, cocycle_checks, period_checks]
        .into_iter()
        .try_fold(0_u128, u128::checked_add)
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "terminal-relative free-cocycle work items",
        })
}

#[allow(clippy::too_many_arguments)]
fn checked_homology_cocycle_accumulate(
    sum: i128,
    left: i128,
    right: i128,
    stage: HomologyCocycleStage,
    row: usize,
    col: usize,
    term: usize,
    completed_scalar: &mut u128,
) -> Result<i128, IntegralTopologyError> {
    let product = left.checked_mul(right).ok_or(
        IntegralTopologyError::HomologyCocycleArithmeticOverflow {
            stage,
            row,
            col,
            term,
        },
    )?;
    let next = sum.checked_add(product).ok_or(
        IntegralTopologyError::HomologyCocycleArithmeticOverflow {
            stage,
            row,
            col,
            term,
        },
    )?;
    *completed_scalar =
        completed_scalar
            .checked_add(1)
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "completed terminal-relative free-cocycle scalar operations",
            })?;
    Ok(next)
}

fn checked_homology_cocycle_transposed_dot(
    left: &ExactIntegerMatrix,
    left_col: usize,
    right: &ExactIntegerMatrix,
    right_col: usize,
    stage: HomologyCocycleStage,
    completed_scalar: &mut u128,
) -> Result<i128, IntegralTopologyError> {
    if left.rows != right.rows {
        return Err(IntegralTopologyError::HomologyCocycleInvariantLost {
            field: "free-cocycle verification row extent",
        });
    }
    let mut sum = 0_i128;
    for term in 0..left.rows {
        sum = checked_homology_cocycle_accumulate(
            sum,
            left.entry(term, left_col),
            right.entry(term, right_col),
            stage,
            right_col,
            left_col,
            term,
            completed_scalar,
        )?;
    }
    Ok(sum)
}

fn poll_homology_cocycles(
    checkpoint: &mut impl FnMut(&'static str) -> bool,
    phase: &'static str,
    completed_work_items: u128,
    planned_work_items: u128,
    completed_scalar_operations: u128,
    planned_scalar_operations: u128,
) -> Result<(), IntegralTopologyError> {
    if checkpoint(phase) {
        Ok(())
    } else {
        Err(IntegralTopologyError::HomologyCocycleCancelled {
            phase,
            completed_work_items,
            planned_work_items,
            completed_scalar_operations,
            planned_scalar_operations,
        })
    }
}

/// Deterministically construct and then independently verify a complete Smith
/// witness without injected cancellation.
pub fn construct_smith_normal_form(
    source: ExactIntegerMatrix,
    construction_budget: SmithConstructionBudget,
    verification_budget: ExactAlgebraBudget,
) -> Result<ConstructedSmithNormalForm, IntegralTopologyError> {
    construct_smith_normal_form_with_checkpoint(
        source,
        construction_budget,
        verification_budget,
        &mut |_| true,
    )
}

/// Deterministically construct and independently verify a complete Smith
/// witness with bounded cancellation polling.
///
/// Pivot selection is minimum unsigned magnitude with row-major tie-breaking.
/// Every transformation is an exact elementary unimodular operation applied
/// to the working matrix, the relevant transform, and the inverse transform
/// on its mathematically required opposite side. A trailing-block entry that
/// is not divisible by the active pivot is mixed into the pivot row and the
/// Euclidean reduction resumes; merely diagonalizing is never accepted as
/// Smith form. The existing witness verifier is the final publication gate.
#[allow(clippy::too_many_lines)]
pub fn construct_smith_normal_form_with_checkpoint(
    source: ExactIntegerMatrix,
    construction_budget: SmithConstructionBudget,
    verification_budget: ExactAlgebraBudget,
    checkpoint: &mut impl FnMut(&'static str) -> bool,
) -> Result<ConstructedSmithNormalForm, IntegralTopologyError> {
    let initial_entry_steps = preflight_smith_construction(
        &source,
        construction_budget,
        verification_budget,
        checkpoint,
    )?;
    poll_smith_construction(
        checkpoint,
        "smith construction preflight",
        0,
        construction_budget.max_elementary_operations,
        initial_entry_steps,
        construction_budget.max_entry_steps,
    )?;

    let rows = source.rows;
    let cols = source.cols;
    let diagonal = try_clone_matrix(&source, "smith construction source clone")?;
    let left = try_identity_matrix(rows, "smith construction left identity")?;
    let left_inverse = try_identity_matrix(rows, "smith construction left inverse identity")?;
    let right = try_identity_matrix(cols, "smith construction right identity")?;
    let right_inverse = try_identity_matrix(cols, "smith construction right inverse identity")?;
    let mut state = SmithConstructionState {
        diagonal,
        left,
        left_inverse,
        right,
        right_inverse,
        elementary_operations: 0,
        entry_steps: initial_entry_steps,
        budget: construction_budget,
    };

    for pivot in 0..rows.min(cols) {
        let Some((pivot_row, pivot_col)) = state.find_pivot(pivot, checkpoint)? else {
            break;
        };
        if pivot_row != pivot {
            state.swap_rows(pivot, pivot_row, checkpoint)?;
        }
        if pivot_col != pivot {
            state.swap_columns(pivot, pivot_col, checkpoint)?;
        }
        if state.inspect(pivot, pivot, "smith construction pivot sign", checkpoint)? < 0 {
            state.normalize_pivot_sign(pivot, checkpoint)?;
        }

        loop {
            let mut restart = false;
            for row in (pivot + 1)..rows {
                let entry =
                    state.inspect(row, pivot, "smith construction column scan", checkpoint)?;
                if entry == 0 {
                    continue;
                }
                let pivot_value =
                    state.inspect(pivot, pivot, "smith construction row quotient", checkpoint)?;
                if pivot_value <= 0 {
                    return Err(IntegralTopologyError::SmithConstructionInvariantLost {
                        field: "row reduction pivot is not positive",
                    });
                }
                let quotient = entry.checked_div_euclid(pivot_value).ok_or(
                    IntegralTopologyError::SmithConstructionArithmeticOverflow {
                        stage: SmithConstructionStage::RowReduction,
                        role: MatrixRole::Diagonal,
                        row,
                        col: pivot,
                    },
                )?;
                state.subtract_row_multiple(row, pivot, quotient, checkpoint)?;
                let remainder =
                    state.inspect(row, pivot, "smith construction row remainder", checkpoint)?;
                if remainder != 0 {
                    if remainder.unsigned_abs() >= pivot_value.unsigned_abs() {
                        return Err(IntegralTopologyError::SmithConstructionInvariantLost {
                            field: "row Euclidean remainder did not reduce pivot",
                        });
                    }
                    state.swap_rows(pivot, row, checkpoint)?;
                }
                restart = true;
                break;
            }
            if restart {
                continue;
            }

            for col in (pivot + 1)..cols {
                let entry = state.inspect(pivot, col, "smith construction row scan", checkpoint)?;
                if entry == 0 {
                    continue;
                }
                let pivot_value = state.inspect(
                    pivot,
                    pivot,
                    "smith construction column quotient",
                    checkpoint,
                )?;
                if pivot_value <= 0 {
                    return Err(IntegralTopologyError::SmithConstructionInvariantLost {
                        field: "column reduction pivot is not positive",
                    });
                }
                let quotient = entry.checked_div_euclid(pivot_value).ok_or(
                    IntegralTopologyError::SmithConstructionArithmeticOverflow {
                        stage: SmithConstructionStage::ColumnReduction,
                        role: MatrixRole::Diagonal,
                        row: pivot,
                        col,
                    },
                )?;
                state.subtract_column_multiple(
                    col,
                    pivot,
                    quotient,
                    SmithConstructionStage::ColumnReduction,
                    "smith construction column reduction",
                    checkpoint,
                )?;
                let remainder = state.inspect(
                    pivot,
                    col,
                    "smith construction column remainder",
                    checkpoint,
                )?;
                if remainder != 0 {
                    if remainder.unsigned_abs() >= pivot_value.unsigned_abs() {
                        return Err(IntegralTopologyError::SmithConstructionInvariantLost {
                            field: "column Euclidean remainder did not reduce pivot",
                        });
                    }
                    state.swap_columns(pivot, col, checkpoint)?;
                }
                restart = true;
                break;
            }
            if restart {
                continue;
            }

            let pivot_value = state.inspect(
                pivot,
                pivot,
                "smith construction divisibility pivot",
                checkpoint,
            )?;
            if pivot_value <= 0 {
                return Err(IntegralTopologyError::SmithConstructionInvariantLost {
                    field: "divisibility pivot is not positive",
                });
            }
            let Some((repair_row, repair_col)) =
                state.find_nondivisible(pivot, pivot_value, checkpoint)?
            else {
                break;
            };
            state.add_row_for_divisibility_repair(pivot, repair_row, checkpoint)?;
            let repair_entry = state.inspect(
                pivot,
                repair_col,
                "smith construction repair entry",
                checkpoint,
            )?;
            let quotient = repair_entry.checked_div_euclid(pivot_value).ok_or(
                IntegralTopologyError::SmithConstructionArithmeticOverflow {
                    stage: SmithConstructionStage::DivisibilityRepair,
                    role: MatrixRole::Diagonal,
                    row: pivot,
                    col: repair_col,
                },
            )?;
            state.subtract_column_multiple(
                repair_col,
                pivot,
                quotient,
                SmithConstructionStage::DivisibilityRepair,
                "smith construction divisibility repair column reduction",
                checkpoint,
            )?;
            let remainder = state.inspect(
                pivot,
                repair_col,
                "smith construction repair remainder",
                checkpoint,
            )?;
            if remainder == 0 || remainder.unsigned_abs() >= pivot_value.unsigned_abs() {
                return Err(IntegralTopologyError::SmithConstructionInvariantLost {
                    field: "divisibility repair did not strictly reduce pivot",
                });
            }
            state.swap_columns(pivot, repair_col, checkpoint)?;
        }
    }

    poll_smith_construction(
        checkpoint,
        "smith construction verification handoff",
        state.elementary_operations,
        construction_budget.max_elementary_operations,
        state.entry_steps,
        construction_budget.max_entry_steps,
    )?;
    let elementary_operations = state.elementary_operations;
    let entry_steps = state.entry_steps;
    let witness = SmithNormalFormWitness::new(
        state.diagonal,
        state.left,
        state.left_inverse,
        state.right,
        state.right_inverse,
    );
    let verified = match verify_smith_normal_form_with_checkpoint(
        source,
        witness,
        verification_budget,
        checkpoint,
    ) {
        Ok(verified) => verified,
        Err(error) if error.failure_class() == IntegralTopologyFailureClass::Refuted => {
            return Err(IntegralTopologyError::SmithConstructionInvariantLost {
                field: "constructed witness failed independent verification",
            });
        }
        Err(error) => return Err(error),
    };
    Ok(ConstructedSmithNormalForm {
        verified,
        elementary_operations,
        entry_steps,
    })
}

#[allow(clippy::too_many_lines)]
fn preflight_smith_construction(
    source: &ExactIntegerMatrix,
    construction_budget: SmithConstructionBudget,
    verification_budget: ExactAlgebraBudget,
    checkpoint: &mut impl FnMut(&'static str) -> bool,
) -> Result<u128, IntegralTopologyError> {
    source.ensure_within(MatrixRole::Source, verification_budget)?;
    if (source.rows > 0 || source.cols > 0) && construction_budget.max_abs_coefficient < 1 {
        return Err(
            IntegralTopologyError::SmithConstructionCoefficientMagnitudeExceeded {
                role: if source.rows > 0 {
                    MatrixRole::LeftTransform
                } else {
                    MatrixRole::RightTransform
                },
                row: 0,
                col: 0,
                magnitude: 1,
                max: construction_budget.max_abs_coefficient,
            },
        );
    }
    let mut initial_entry_steps = 0_u128;
    let mut index = 0_usize;
    while index < source.entries.len() {
        let requested =
            initial_entry_steps
                .checked_add(1)
                .ok_or(IntegralTopologyError::WorkPlanOverflow {
                    phase: "smith construction source coefficient scan",
                })?;
        if requested > construction_budget.max_entry_steps {
            return Err(
                IntegralTopologyError::SmithConstructionEntryStepBudgetExceeded {
                    requested,
                    max: construction_budget.max_entry_steps,
                },
            );
        }
        poll_smith_construction(
            checkpoint,
            "smith construction source coefficient scan",
            0,
            construction_budget.max_elementary_operations,
            initial_entry_steps,
            construction_budget.max_entry_steps,
        )?;
        initial_entry_steps = requested;
        let value = source.entries[index];
        let magnitude = value.unsigned_abs();
        if magnitude > construction_budget.max_abs_coefficient {
            return Err(
                IntegralTopologyError::SmithConstructionCoefficientMagnitudeExceeded {
                    role: MatrixRole::Source,
                    row: index / source.cols.max(1),
                    col: index % source.cols.max(1),
                    magnitude,
                    max: construction_budget.max_abs_coefficient,
                },
            );
        }
        index += 1;
    }
    let rows_squared =
        source
            .rows
            .checked_mul(source.rows)
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "smith construction left entries",
            })?;
    let cols_squared =
        source
            .cols
            .checked_mul(source.cols)
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "smith construction right entries",
            })?;
    for entries in [source.entries.len(), rows_squared, cols_squared] {
        if entries > verification_budget.max_matrix_entries {
            return Err(IntegralTopologyError::MatrixEntryBudgetExceeded {
                requested: entries,
                max: verification_budget.max_matrix_entries,
            });
        }
    }
    let live_entries = source
        .entries
        .len()
        .checked_mul(2)
        .and_then(|entries| {
            rows_squared
                .checked_mul(2)
                .and_then(|left| entries.checked_add(left))
        })
        .and_then(|entries| {
            cols_squared
                .checked_mul(2)
                .and_then(|right| entries.checked_add(right))
        })
        .ok_or(IntegralTopologyError::WorkPlanOverflow {
            phase: "smith construction live entries",
        })?;
    if live_entries > construction_budget.max_live_entries {
        return Err(
            IntegralTopologyError::SmithConstructionLiveEntryBudgetExceeded {
                requested: live_entries,
                max: construction_budget.max_live_entries,
            },
        );
    }
    if live_entries > verification_budget.max_retained_entries {
        return Err(IntegralTopologyError::RetainedEntryBudgetExceeded {
            requested: live_entries,
            max: verification_budget.max_retained_entries,
        });
    }
    if source.entries.len() > verification_budget.max_workspace_entries {
        return Err(IntegralTopologyError::WorkspaceEntryBudgetExceeded {
            requested: source.entries.len(),
            max: verification_budget.max_workspace_entries,
        });
    }
    Ok(initial_entry_steps)
}

fn try_clone_matrix(
    matrix: &ExactIntegerMatrix,
    phase: &'static str,
) -> Result<ExactIntegerMatrix, IntegralTopologyError> {
    let mut entries = Vec::new();
    entries
        .try_reserve_exact(matrix.entries.len())
        .map_err(|_| IntegralTopologyError::AllocationRefused {
            phase,
            requested_entries: matrix.entries.len(),
        })?;
    entries.extend_from_slice(&matrix.entries);
    Ok(ExactIntegerMatrix {
        rows: matrix.rows,
        cols: matrix.cols,
        entries,
    })
}

fn try_identity_matrix(
    extent: usize,
    phase: &'static str,
) -> Result<ExactIntegerMatrix, IntegralTopologyError> {
    let entries = extent
        .checked_mul(extent)
        .ok_or(IntegralTopologyError::WorkPlanOverflow { phase })?;
    let mut matrix = ExactIntegerMatrix {
        rows: extent,
        cols: extent,
        entries: allocate_zeroed(entries, phase)?,
    };
    for index in 0..extent {
        matrix.entries[index * extent + index] = 1;
    }
    Ok(matrix)
}

struct SmithConstructionState {
    diagonal: ExactIntegerMatrix,
    left: ExactIntegerMatrix,
    left_inverse: ExactIntegerMatrix,
    right: ExactIntegerMatrix,
    right_inverse: ExactIntegerMatrix,
    elementary_operations: u128,
    entry_steps: u128,
    budget: SmithConstructionBudget,
}

impl SmithConstructionState {
    fn inspect(
        &mut self,
        row: usize,
        col: usize,
        phase: &'static str,
        checkpoint: &mut impl FnMut(&'static str) -> bool,
    ) -> Result<i128, IntegralTopologyError> {
        self.reserve_scan_steps(1, phase, checkpoint)?;
        Ok(self.diagonal.entry(row, col))
    }

    fn find_pivot(
        &mut self,
        offset: usize,
        checkpoint: &mut impl FnMut(&'static str) -> bool,
    ) -> Result<Option<(usize, usize)>, IntegralTopologyError> {
        let mut best: Option<(u128, usize, usize)> = None;
        for row in offset..self.diagonal.rows {
            for col in offset..self.diagonal.cols {
                let value = self.inspect(row, col, "smith construction pivot scan", checkpoint)?;
                if value == 0 {
                    continue;
                }
                let candidate = (value.unsigned_abs(), row, col);
                if best.is_none_or(|current| candidate < current) {
                    best = Some(candidate);
                }
            }
        }
        Ok(best.map(|(_, row, col)| (row, col)))
    }

    fn find_nondivisible(
        &mut self,
        pivot: usize,
        divisor: i128,
        checkpoint: &mut impl FnMut(&'static str) -> bool,
    ) -> Result<Option<(usize, usize)>, IntegralTopologyError> {
        debug_assert!(divisor > 0);
        for row in (pivot + 1)..self.diagonal.rows {
            for col in (pivot + 1)..self.diagonal.cols {
                let value =
                    self.inspect(row, col, "smith construction divisibility scan", checkpoint)?;
                let remainder = value.checked_rem_euclid(divisor).ok_or(
                    IntegralTopologyError::SmithConstructionArithmeticOverflow {
                        stage: SmithConstructionStage::DivisibilityRepair,
                        role: MatrixRole::Diagonal,
                        row,
                        col,
                    },
                )?;
                if remainder != 0 {
                    return Ok(Some((row, col)));
                }
            }
        }
        Ok(None)
    }

    fn swap_rows(
        &mut self,
        first: usize,
        second: usize,
        checkpoint: &mut impl FnMut(&'static str) -> bool,
    ) -> Result<(), IntegralTopologyError> {
        let steps = self
            .diagonal
            .cols
            .checked_add(self.left.cols)
            .and_then(|count| count.checked_add(self.left_inverse.rows))
            .and_then(|count| count.checked_mul(2))
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "smith construction row swap steps",
            })?;
        self.begin_operation(steps, "smith construction row swap", checkpoint)?;
        swap_matrix_rows(&mut self.diagonal, first, second);
        swap_matrix_rows(&mut self.left, first, second);
        swap_matrix_columns(&mut self.left_inverse, first, second);
        Ok(())
    }

    fn swap_columns(
        &mut self,
        first: usize,
        second: usize,
        checkpoint: &mut impl FnMut(&'static str) -> bool,
    ) -> Result<(), IntegralTopologyError> {
        let steps = self
            .diagonal
            .rows
            .checked_add(self.right.rows)
            .and_then(|count| count.checked_add(self.right_inverse.cols))
            .and_then(|count| count.checked_mul(2))
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "smith construction column swap steps",
            })?;
        self.begin_operation(steps, "smith construction column swap", checkpoint)?;
        swap_matrix_columns(&mut self.diagonal, first, second);
        swap_matrix_columns(&mut self.right, first, second);
        swap_matrix_rows(&mut self.right_inverse, first, second);
        Ok(())
    }

    fn subtract_row_multiple(
        &mut self,
        target: usize,
        source: usize,
        quotient: i128,
        checkpoint: &mut impl FnMut(&'static str) -> bool,
    ) -> Result<(), IntegralTopologyError> {
        let steps = self
            .diagonal
            .cols
            .checked_add(self.left.cols)
            .and_then(|count| count.checked_add(self.left_inverse.rows))
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "smith construction row reduction steps",
            })?;
        self.begin_operation(steps, "smith construction row reduction", checkpoint)?;
        checked_subtract_matrix_row(
            &mut self.diagonal,
            target,
            source,
            quotient,
            SmithConstructionStage::RowReduction,
            MatrixRole::Diagonal,
            self.budget.max_abs_coefficient,
        )?;
        checked_subtract_matrix_row(
            &mut self.left,
            target,
            source,
            quotient,
            SmithConstructionStage::RowReduction,
            MatrixRole::LeftTransform,
            self.budget.max_abs_coefficient,
        )?;
        checked_add_matrix_column(
            &mut self.left_inverse,
            source,
            target,
            quotient,
            SmithConstructionStage::RowReduction,
            MatrixRole::LeftInverse,
            self.budget.max_abs_coefficient,
        )
    }

    fn subtract_column_multiple(
        &mut self,
        target: usize,
        source: usize,
        quotient: i128,
        stage: SmithConstructionStage,
        phase: &'static str,
        checkpoint: &mut impl FnMut(&'static str) -> bool,
    ) -> Result<(), IntegralTopologyError> {
        let steps = self
            .diagonal
            .rows
            .checked_add(self.right.rows)
            .and_then(|count| count.checked_add(self.right_inverse.cols))
            .ok_or(IntegralTopologyError::WorkPlanOverflow { phase })?;
        self.begin_operation(steps, phase, checkpoint)?;
        checked_subtract_matrix_column(
            &mut self.diagonal,
            target,
            source,
            quotient,
            stage,
            MatrixRole::Diagonal,
            self.budget.max_abs_coefficient,
        )?;
        checked_subtract_matrix_column(
            &mut self.right,
            target,
            source,
            quotient,
            stage,
            MatrixRole::RightTransform,
            self.budget.max_abs_coefficient,
        )?;
        checked_add_matrix_row(
            &mut self.right_inverse,
            source,
            target,
            quotient,
            stage,
            MatrixRole::RightInverse,
            self.budget.max_abs_coefficient,
        )
    }

    fn add_row_for_divisibility_repair(
        &mut self,
        target: usize,
        source: usize,
        checkpoint: &mut impl FnMut(&'static str) -> bool,
    ) -> Result<(), IntegralTopologyError> {
        let steps = self
            .diagonal
            .cols
            .checked_add(self.left.cols)
            .and_then(|count| count.checked_add(self.left_inverse.rows))
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "smith construction divisibility repair steps",
            })?;
        self.begin_operation(steps, "smith construction divisibility repair", checkpoint)?;
        checked_add_matrix_row(
            &mut self.diagonal,
            target,
            source,
            1,
            SmithConstructionStage::DivisibilityRepair,
            MatrixRole::Diagonal,
            self.budget.max_abs_coefficient,
        )?;
        checked_add_matrix_row(
            &mut self.left,
            target,
            source,
            1,
            SmithConstructionStage::DivisibilityRepair,
            MatrixRole::LeftTransform,
            self.budget.max_abs_coefficient,
        )?;
        checked_subtract_matrix_column(
            &mut self.left_inverse,
            source,
            target,
            1,
            SmithConstructionStage::DivisibilityRepair,
            MatrixRole::LeftInverse,
            self.budget.max_abs_coefficient,
        )
    }

    fn normalize_pivot_sign(
        &mut self,
        pivot: usize,
        checkpoint: &mut impl FnMut(&'static str) -> bool,
    ) -> Result<(), IntegralTopologyError> {
        let steps = self
            .diagonal
            .cols
            .checked_add(self.left.cols)
            .and_then(|count| count.checked_add(self.left_inverse.rows))
            .ok_or(IntegralTopologyError::WorkPlanOverflow {
                phase: "smith construction sign steps",
            })?;
        self.begin_operation(steps, "smith construction sign normalization", checkpoint)?;
        checked_negate_matrix_row(
            &mut self.diagonal,
            pivot,
            SmithConstructionStage::SignNormalization,
            MatrixRole::Diagonal,
            self.budget.max_abs_coefficient,
        )?;
        checked_negate_matrix_row(
            &mut self.left,
            pivot,
            SmithConstructionStage::SignNormalization,
            MatrixRole::LeftTransform,
            self.budget.max_abs_coefficient,
        )?;
        checked_negate_matrix_column(
            &mut self.left_inverse,
            pivot,
            SmithConstructionStage::SignNormalization,
            MatrixRole::LeftInverse,
            self.budget.max_abs_coefficient,
        )
    }

    fn reserve_scan_steps(
        &mut self,
        steps: usize,
        phase: &'static str,
        checkpoint: &mut impl FnMut(&'static str) -> bool,
    ) -> Result<(), IntegralTopologyError> {
        let steps = u128::try_from(steps).map_err(|_| IntegralTopologyError::WorkPlanOverflow {
            phase: "smith construction scan step conversion",
        })?;
        let requested =
            self.entry_steps
                .checked_add(steps)
                .ok_or(IntegralTopologyError::WorkPlanOverflow {
                    phase: "smith construction scan steps",
                })?;
        if requested > self.budget.max_entry_steps {
            return Err(
                IntegralTopologyError::SmithConstructionEntryStepBudgetExceeded {
                    requested,
                    max: self.budget.max_entry_steps,
                },
            );
        }
        poll_smith_construction(
            checkpoint,
            phase,
            self.elementary_operations,
            self.budget.max_elementary_operations,
            self.entry_steps,
            self.budget.max_entry_steps,
        )?;
        self.entry_steps = requested;
        Ok(())
    }

    fn begin_operation(
        &mut self,
        steps: usize,
        phase: &'static str,
        checkpoint: &mut impl FnMut(&'static str) -> bool,
    ) -> Result<(), IntegralTopologyError> {
        let requested_operations = self.elementary_operations.checked_add(1).ok_or(
            IntegralTopologyError::WorkPlanOverflow {
                phase: "smith construction operation count",
            },
        )?;
        if requested_operations > self.budget.max_elementary_operations {
            return Err(
                IntegralTopologyError::SmithConstructionOperationBudgetExceeded {
                    requested: requested_operations,
                    max: self.budget.max_elementary_operations,
                },
            );
        }
        let steps = u128::try_from(steps).map_err(|_| IntegralTopologyError::WorkPlanOverflow {
            phase: "smith construction update step conversion",
        })?;
        let requested_steps =
            self.entry_steps
                .checked_add(steps)
                .ok_or(IntegralTopologyError::WorkPlanOverflow {
                    phase: "smith construction update steps",
                })?;
        if requested_steps > self.budget.max_entry_steps {
            return Err(
                IntegralTopologyError::SmithConstructionEntryStepBudgetExceeded {
                    requested: requested_steps,
                    max: self.budget.max_entry_steps,
                },
            );
        }
        poll_smith_construction(
            checkpoint,
            phase,
            self.elementary_operations,
            self.budget.max_elementary_operations,
            self.entry_steps,
            self.budget.max_entry_steps,
        )?;
        self.elementary_operations = requested_operations;
        self.entry_steps = requested_steps;
        Ok(())
    }
}

fn poll_smith_construction(
    checkpoint: &mut impl FnMut(&'static str) -> bool,
    phase: &'static str,
    completed_operations: u128,
    max_operations: u128,
    completed_entry_steps: u128,
    max_entry_steps: u128,
) -> Result<(), IntegralTopologyError> {
    if checkpoint(phase) {
        Ok(())
    } else {
        Err(IntegralTopologyError::SmithConstructionCancelled {
            phase,
            completed_operations,
            max_operations,
            completed_entry_steps,
            max_entry_steps,
        })
    }
}

fn swap_matrix_rows(matrix: &mut ExactIntegerMatrix, first: usize, second: usize) {
    for col in 0..matrix.cols {
        matrix
            .entries
            .swap(first * matrix.cols + col, second * matrix.cols + col);
    }
}

fn swap_matrix_columns(matrix: &mut ExactIntegerMatrix, first: usize, second: usize) {
    for row in 0..matrix.rows {
        matrix
            .entries
            .swap(row * matrix.cols + first, row * matrix.cols + second);
    }
}

fn checked_subtract_matrix_row(
    matrix: &mut ExactIntegerMatrix,
    target: usize,
    source: usize,
    multiplier: i128,
    stage: SmithConstructionStage,
    role: MatrixRole,
    max_abs_coefficient: u128,
) -> Result<(), IntegralTopologyError> {
    for col in 0..matrix.cols {
        let source_value = matrix.entries[source * matrix.cols + col];
        let target_index = target * matrix.cols + col;
        let value = matrix.entries[target_index]
            .checked_sub(checked_construction_product(
                source_value,
                multiplier,
                stage,
                role,
                target,
                col,
            )?)
            .ok_or(IntegralTopologyError::SmithConstructionArithmeticOverflow {
                stage,
                role,
                row: target,
                col,
            })?;
        matrix.entries[target_index] =
            checked_construction_coefficient(value, max_abs_coefficient, role, target, col)?;
    }
    Ok(())
}

fn checked_add_matrix_row(
    matrix: &mut ExactIntegerMatrix,
    target: usize,
    source: usize,
    multiplier: i128,
    stage: SmithConstructionStage,
    role: MatrixRole,
    max_abs_coefficient: u128,
) -> Result<(), IntegralTopologyError> {
    for col in 0..matrix.cols {
        let source_value = matrix.entries[source * matrix.cols + col];
        let target_index = target * matrix.cols + col;
        let value = matrix.entries[target_index]
            .checked_add(checked_construction_product(
                source_value,
                multiplier,
                stage,
                role,
                target,
                col,
            )?)
            .ok_or(IntegralTopologyError::SmithConstructionArithmeticOverflow {
                stage,
                role,
                row: target,
                col,
            })?;
        matrix.entries[target_index] =
            checked_construction_coefficient(value, max_abs_coefficient, role, target, col)?;
    }
    Ok(())
}

fn checked_subtract_matrix_column(
    matrix: &mut ExactIntegerMatrix,
    target: usize,
    source: usize,
    multiplier: i128,
    stage: SmithConstructionStage,
    role: MatrixRole,
    max_abs_coefficient: u128,
) -> Result<(), IntegralTopologyError> {
    for row in 0..matrix.rows {
        let source_value = matrix.entries[row * matrix.cols + source];
        let target_index = row * matrix.cols + target;
        let value = matrix.entries[target_index]
            .checked_sub(checked_construction_product(
                source_value,
                multiplier,
                stage,
                role,
                row,
                target,
            )?)
            .ok_or(IntegralTopologyError::SmithConstructionArithmeticOverflow {
                stage,
                role,
                row,
                col: target,
            })?;
        matrix.entries[target_index] =
            checked_construction_coefficient(value, max_abs_coefficient, role, row, target)?;
    }
    Ok(())
}

fn checked_add_matrix_column(
    matrix: &mut ExactIntegerMatrix,
    target: usize,
    source: usize,
    multiplier: i128,
    stage: SmithConstructionStage,
    role: MatrixRole,
    max_abs_coefficient: u128,
) -> Result<(), IntegralTopologyError> {
    for row in 0..matrix.rows {
        let source_value = matrix.entries[row * matrix.cols + source];
        let target_index = row * matrix.cols + target;
        let value = matrix.entries[target_index]
            .checked_add(checked_construction_product(
                source_value,
                multiplier,
                stage,
                role,
                row,
                target,
            )?)
            .ok_or(IntegralTopologyError::SmithConstructionArithmeticOverflow {
                stage,
                role,
                row,
                col: target,
            })?;
        matrix.entries[target_index] =
            checked_construction_coefficient(value, max_abs_coefficient, role, row, target)?;
    }
    Ok(())
}

fn checked_negate_matrix_row(
    matrix: &mut ExactIntegerMatrix,
    row: usize,
    stage: SmithConstructionStage,
    role: MatrixRole,
    max_abs_coefficient: u128,
) -> Result<(), IntegralTopologyError> {
    for col in 0..matrix.cols {
        let index = row * matrix.cols + col;
        let value = matrix.entries[index].checked_neg().ok_or(
            IntegralTopologyError::SmithConstructionArithmeticOverflow {
                stage,
                role,
                row,
                col,
            },
        )?;
        matrix.entries[index] =
            checked_construction_coefficient(value, max_abs_coefficient, role, row, col)?;
    }
    Ok(())
}

fn checked_negate_matrix_column(
    matrix: &mut ExactIntegerMatrix,
    col: usize,
    stage: SmithConstructionStage,
    role: MatrixRole,
    max_abs_coefficient: u128,
) -> Result<(), IntegralTopologyError> {
    for row in 0..matrix.rows {
        let index = row * matrix.cols + col;
        let value = matrix.entries[index].checked_neg().ok_or(
            IntegralTopologyError::SmithConstructionArithmeticOverflow {
                stage,
                role,
                row,
                col,
            },
        )?;
        matrix.entries[index] =
            checked_construction_coefficient(value, max_abs_coefficient, role, row, col)?;
    }
    Ok(())
}

fn checked_construction_product(
    left: i128,
    right: i128,
    stage: SmithConstructionStage,
    role: MatrixRole,
    row: usize,
    col: usize,
) -> Result<i128, IntegralTopologyError> {
    left.checked_mul(right)
        .ok_or(IntegralTopologyError::SmithConstructionArithmeticOverflow {
            stage,
            role,
            row,
            col,
        })
}

fn checked_construction_coefficient(
    value: i128,
    max_abs_coefficient: u128,
    role: MatrixRole,
    row: usize,
    col: usize,
) -> Result<i128, IntegralTopologyError> {
    let magnitude = value.unsigned_abs();
    if magnitude > max_abs_coefficient {
        Err(
            IntegralTopologyError::SmithConstructionCoefficientMagnitudeExceeded {
                role,
                row,
                col,
                magnitude,
                max: max_abs_coefficient,
            },
        )
    } else {
        Ok(value)
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
        // i128::is_multiple_of is unavailable on the pinned nightly; this
        // is its exact semantics (value != 0 on this path).
        if let Some(previous) = invariant_factors.last().copied()
            && (previous == 0 || value % previous != 0)
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
    /// The requested phase has no admitted component binding in the pair.
    UnknownTerminalRelativePhase {
        /// Rejected phase identity.
        phase: String,
    },
    /// Boundary degree exceeded the unaugmented chain-complex edge.
    BoundaryDegreeOutOfRange {
        /// Requested source degree.
        degree: u8,
        /// Highest admitted source degree (`complex dimension + 1`).
        max: u8,
    },
    /// An opaque admitted pair lost its retained phase/component invariant.
    PhaseComponentBindingLost {
        /// Phase whose component could not be recovered.
        phase: String,
        /// Component identity retained by the phase map.
        component: String,
    },
    /// Component-support traversal exceeded its explicit envelope.
    ComponentVisitBudgetExceeded {
        /// Required support-cell visits.
        requested: usize,
        /// Maximum admitted visits.
        max: usize,
    },
    /// Ambient incidence traversal exceeded its explicit envelope.
    IncidenceVisitBudgetExceeded {
        /// Required incidence visits.
        requested: usize,
        /// Maximum admitted visits.
        max: usize,
    },
    /// Cancellation was observed while binding a pair boundary.
    PairBoundaryCancelled {
        /// Observation phase.
        phase: &'static str,
        /// Deterministic component/incidence work items completed.
        completed_work_items: u128,
        /// Deterministic component/incidence work items planned.
        planned_work_items: u128,
    },
    /// Adjacent boundaries disagreed on a retained semantic binding.
    AdjacentBoundaryBindingMismatch {
        /// First mismatched binding field.
        field: &'static str,
    },
    /// Candidate boundaries were not consecutive maps `A_k`, `A_(k+1)`.
    AdjacentBoundaryDegreeMismatch {
        /// Outgoing source degree `k`.
        outgoing: u8,
        /// Supplied incoming source degree.
        incoming_source: u8,
        /// Supplied incoming target degree.
        incoming_target: Option<u8>,
    },
    /// Shared `C_k` bases differed in cell identity or order.
    AdjacentBoundaryBasisMismatch {
        /// First differing index.
        index: usize,
        /// Outgoing source-basis cell, if present.
        outgoing: Option<CellRef>,
        /// Incoming target-basis cell, if present.
        incoming: Option<CellRef>,
    },
    /// An opaque verified/boundary value violated an internal shape invariant.
    KernelCoordinateInvariantLost {
        /// Broken invariant field.
        field: &'static str,
    },
    /// The supplied Smith authority described a different source shape.
    OutgoingSmithSourceShapeMismatch {
        /// Required outgoing rows.
        expected_rows: usize,
        /// Required shared-chain columns.
        expected_cols: usize,
        /// Smith-bound source rows.
        actual_rows: usize,
        /// Smith-bound source columns.
        actual_cols: usize,
    },
    /// The supplied Smith authority described different outgoing matrix bytes.
    OutgoingSmithSourceEntryMismatch {
        /// First differing row.
        row: usize,
        /// First differing column.
        col: usize,
        /// Pair-bound outgoing entry.
        expected: i128,
        /// Smith-bound source entry.
        actual: i128,
    },
    /// One adjacent-boundary extent exceeded the transport envelope.
    KernelCoordinateExtentExceeded {
        /// Outgoing target extent `m`.
        outgoing_rows: usize,
        /// Shared chain extent `n`.
        chain_extent: usize,
        /// Incoming source extent `p`.
        incoming_cols: usize,
        /// Maximum extent for each axis.
        max: usize,
    },
    /// Exact source/basis binding comparisons exceeded their envelope.
    KernelBindingBudgetExceeded {
        /// Required comparison items.
        requested: usize,
        /// Maximum admitted items.
        max: usize,
    },
    /// `V^-1 * A_(k+1)` had a nonzero nonkernel coordinate.
    IncomingImageOutsideKernel {
        /// Smith-coordinate row below `rank(A_k)`.
        row: usize,
        /// Incoming basis column.
        col: usize,
        /// Exact nonzero coordinate.
        coordinate: i128,
        /// Positive Smith invariant factor proving this is nonkernel.
        invariant_factor: i128,
    },
    /// Cancellation was observed before kernel-coordinate publication.
    KernelCoordinateCancelled {
        /// Observation phase.
        phase: &'static str,
        /// Completed binding comparisons.
        completed_binding_items: usize,
        /// Planned binding comparisons.
        planned_binding_items: usize,
        /// Completed checked dot-product terms.
        completed_scalar_operations: u128,
        /// Planned checked dot-product terms.
        planned_scalar_operations: u128,
    },
    /// Simultaneously live source and construction matrices exceeded their
    /// explicit storage envelope.
    SmithConstructionLiveEntryBudgetExceeded {
        /// Required live entries.
        requested: usize,
        /// Maximum admitted live entries.
        max: usize,
    },
    /// Another elementary operation would exceed the construction envelope.
    SmithConstructionOperationBudgetExceeded {
        /// Operation count including the refused next operation.
        requested: u128,
        /// Maximum admitted elementary operations.
        max: u128,
    },
    /// Exact construction inspections/destination updates exceeded their
    /// envelope.
    SmithConstructionEntryStepBudgetExceeded {
        /// Step count including the refused inspection/update block.
        requested: u128,
        /// Maximum admitted entry steps.
        max: u128,
    },
    /// A source or generated coefficient exceeded its explicit magnitude cap.
    SmithConstructionCoefficientMagnitudeExceeded {
        /// Matrix containing the refused coefficient.
        role: MatrixRole,
        /// Coefficient row.
        row: usize,
        /// Coefficient column.
        col: usize,
        /// Unsigned exact magnitude.
        magnitude: u128,
        /// Maximum admitted magnitude.
        max: u128,
    },
    /// Checked constructive row/column arithmetic overflowed `i128`.
    SmithConstructionArithmeticOverflow {
        /// Elementary operation phase.
        stage: SmithConstructionStage,
        /// Matrix whose update or quotient refused.
        role: MatrixRole,
        /// Refusing row.
        row: usize,
        /// Refusing column.
        col: usize,
    },
    /// A runtime construction invariant or its independent witness check
    /// failed. This is inconclusive for the mathematically valid source.
    SmithConstructionInvariantLost {
        /// Broken invariant.
        field: &'static str,
    },
    /// Cancellation was observed during data-dependent constructive reduction.
    SmithConstructionCancelled {
        /// Observation phase.
        phase: &'static str,
        /// Completed elementary operations.
        completed_operations: u128,
        /// Maximum admitted elementary operations.
        max_operations: u128,
        /// Completed source/diagonal inspections and destination updates.
        completed_entry_steps: u128,
        /// Maximum admitted source/diagonal inspections and destination
        /// updates.
        max_entry_steps: u128,
    },
    /// The supplied lower-image Smith authority described another matrix
    /// shape.
    HomologySmithSourceShapeMismatch {
        /// Required kernel-image rows.
        expected_rows: usize,
        /// Required incoming-chain columns.
        expected_cols: usize,
        /// Smith-bound source rows.
        actual_rows: usize,
        /// Smith-bound source columns.
        actual_cols: usize,
    },
    /// The supplied lower-image Smith authority described different matrix
    /// bytes.
    HomologySmithSourceEntryMismatch {
        /// First differing row.
        row: usize,
        /// First differing column.
        col: usize,
        /// Kernel-transport entry.
        expected: i128,
        /// Smith-bound source entry.
        actual: i128,
    },
    /// Exact source/factor binding work exceeded its envelope.
    HomologyBindingBudgetExceeded {
        /// Required binding items.
        requested: usize,
        /// Maximum admitted binding items.
        max: usize,
    },
    /// Opaque verified inputs violated an internal quotient invariant.
    HomologyDecompositionInvariantLost {
        /// Broken invariant.
        field: &'static str,
    },
    /// Cancellation was observed before quotient-homology publication.
    HomologyDecompositionCancelled {
        /// Observation phase.
        phase: &'static str,
        /// Completed source/factor inspections.
        completed_binding_items: usize,
        /// Planned source/factor inspections.
        planned_binding_items: usize,
    },
    /// Generator and torsion-filling output exceeded its aggregate envelope.
    HomologyGeneratorOutputBudgetExceeded {
        /// Required retained output coefficients.
        requested: usize,
        /// Maximum admitted output coefficients.
        max: usize,
    },
    /// An opaque homology authority violated a generator-lift invariant.
    HomologyGeneratorInvariantLost {
        /// Broken invariant field.
        field: &'static str,
    },
    /// Checked generator construction or equation arithmetic overflowed.
    HomologyGeneratorArithmeticOverflow {
        /// Generator calculation or verification phase.
        stage: HomologyGeneratorStage,
        /// Output or verification row.
        row: usize,
        /// Generator column.
        col: usize,
        /// Inner-product term, or the incoming column count for torsion
        /// scaling.
        term: usize,
    },
    /// A constructed generator violated an exact chain-level equation.
    ///
    /// Because the input authority is already opaque and verified, this is an
    /// internal composition failure rather than a refutation of user data.
    HomologyGeneratorVerificationMismatch {
        /// Equation that disagreed.
        stage: HomologyGeneratorStage,
        /// Equation row.
        row: usize,
        /// Generator column.
        col: usize,
        /// Exact required value.
        expected: i128,
        /// Exact observed value.
        actual: i128,
    },
    /// Cancellation was observed before generator publication.
    HomologyGeneratorCancelled {
        /// Observation phase.
        phase: &'static str,
        /// Completed deterministic output/check work items.
        completed_work_items: u128,
        /// Planned deterministic output/check work items.
        planned_work_items: u128,
        /// Completed checked scalar terms.
        completed_scalar_operations: u128,
        /// Planned checked scalar terms.
        planned_scalar_operations: u128,
    },
    /// Free-cocycle and period-pairing output exceeded its aggregate envelope.
    HomologyCocycleOutputBudgetExceeded {
        /// Required retained output coefficients.
        requested: usize,
        /// Maximum admitted output coefficients.
        max: usize,
    },
    /// An opaque generator authority violated a free-cocycle invariant.
    HomologyCocycleInvariantLost {
        /// Broken invariant field.
        field: &'static str,
    },
    /// Checked cocycle construction or pairing arithmetic overflowed.
    HomologyCocycleArithmeticOverflow {
        /// Cocycle calculation or verification phase.
        stage: HomologyCocycleStage,
        /// Output or equation row.
        row: usize,
        /// Cocycle column.
        col: usize,
        /// Inner-product term.
        term: usize,
    },
    /// A constructed cocycle violated its exact coboundary or period equation.
    ///
    /// The input generator authority is already opaque and verified, so this
    /// indicates an internal composition failure rather than refuted user data.
    HomologyCocycleVerificationMismatch {
        /// Equation that disagreed.
        stage: HomologyCocycleStage,
        /// Equation row.
        row: usize,
        /// Cocycle column.
        col: usize,
        /// Exact required value.
        expected: i128,
        /// Exact observed value.
        actual: i128,
    },
    /// Cancellation was observed before free-cocycle publication.
    HomologyCocycleCancelled {
        /// Observation phase.
        phase: &'static str,
        /// Completed deterministic output/check work items.
        completed_work_items: u128,
        /// Planned deterministic output/check work items.
        planned_work_items: u128,
        /// Completed checked scalar terms.
        completed_scalar_operations: u128,
        /// Planned checked scalar terms.
        planned_scalar_operations: u128,
    },
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
    /// Complete retained matrix/basis or source/witness storage exceeded the
    /// envelope.
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
    /// Internal exact matrix, basis, or workspace allocation refused.
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
            Self::UnknownTerminalRelativePhase { .. }
            | Self::BoundaryDegreeOutOfRange { .. }
            | Self::AdjacentBoundaryBindingMismatch { .. }
            | Self::AdjacentBoundaryDegreeMismatch { .. }
            | Self::AdjacentBoundaryBasisMismatch { .. }
            | Self::OutgoingSmithSourceShapeMismatch { .. }
            | Self::OutgoingSmithSourceEntryMismatch { .. }
            | Self::IncomingImageOutsideKernel { .. }
            | Self::HomologySmithSourceShapeMismatch { .. }
            | Self::HomologySmithSourceEntryMismatch { .. }
            | Self::MatrixEntryCount { .. }
            | Self::WitnessShape { .. }
            | Self::WitnessProductMismatch { .. }
            | Self::OffDiagonalEntry { .. }
            | Self::NegativeInvariantFactor { .. }
            | Self::NonzeroAfterZero { .. }
            | Self::InvariantFactorDivisibility { .. } => IntegralTopologyFailureClass::Refuted,
            Self::PhaseComponentBindingLost { .. }
            | Self::ComponentVisitBudgetExceeded { .. }
            | Self::IncidenceVisitBudgetExceeded { .. }
            | Self::PairBoundaryCancelled { .. }
            | Self::KernelCoordinateInvariantLost { .. }
            | Self::KernelCoordinateExtentExceeded { .. }
            | Self::KernelBindingBudgetExceeded { .. }
            | Self::KernelCoordinateCancelled { .. }
            | Self::SmithConstructionLiveEntryBudgetExceeded { .. }
            | Self::SmithConstructionOperationBudgetExceeded { .. }
            | Self::SmithConstructionEntryStepBudgetExceeded { .. }
            | Self::SmithConstructionCoefficientMagnitudeExceeded { .. }
            | Self::SmithConstructionArithmeticOverflow { .. }
            | Self::SmithConstructionInvariantLost { .. }
            | Self::SmithConstructionCancelled { .. }
            | Self::HomologyBindingBudgetExceeded { .. }
            | Self::HomologyDecompositionInvariantLost { .. }
            | Self::HomologyDecompositionCancelled { .. }
            | Self::HomologyGeneratorOutputBudgetExceeded { .. }
            | Self::HomologyGeneratorInvariantLost { .. }
            | Self::HomologyGeneratorArithmeticOverflow { .. }
            | Self::HomologyGeneratorVerificationMismatch { .. }
            | Self::HomologyGeneratorCancelled { .. }
            | Self::HomologyCocycleOutputBudgetExceeded { .. }
            | Self::HomologyCocycleInvariantLost { .. }
            | Self::HomologyCocycleArithmeticOverflow { .. }
            | Self::HomologyCocycleVerificationMismatch { .. }
            | Self::HomologyCocycleCancelled { .. }
            | Self::MatrixExtentExceeded { .. }
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
