//! Original mixture sources refreshed after conservative field transport.
//!
//! A diffuse-interface cell retains its homogeneous branch, including spinodal
//! compositions. The field already represents phase separation, so no cell is
//! independently replaced by a bulk equilibrium split. The borrowed source
//! binding stays unchanged across every accepted transport interval.
//!
//! Source receipt allocations are separate from the numerical graph workspace;
//! an explicit cell-count budget bounds queries and retained material results.
//! Density, viscosity and thermal properties are exposed through the existing
//! liquid consumers. This adapter does not advance momentum or heat feedback.

use super::advection::{BinaryMiscibilityAdvection, BinaryMiscibilityAdvectionControls};
use super::transport::BinaryMiscibilityTransportStep;
use super::{
    BinaryMiscibilityError, BinaryMiscibilityNetwork, BinaryMiscibilitySolve,
    BinaryMiscibilityState, poll,
};
use fs_material::liquid::solution::BinaryLiquidError;
use fs_material::liquid::solution::gradient::{
    BinaryLiquidSquareGradientLaw, ResolvedBinaryLiquidCell,
};
use std::sync::Arc;

/// Field/source incompatibility, work budget, or an atomic solver/source refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum BinaryMiscibilityMaterialError {
    /// The original gradient model and source binding do not agree.
    Invalid(&'static str),
    /// Maximum number of source-resolved cells permitted by the caller.
    CellBudget { required: usize, maximum: usize },
    /// The original numerical field could not produce an admissible candidate.
    Solver(BinaryMiscibilityError),
    /// A local source rejected the unchanged accepted paired inventory.
    Source {
        cell: usize,
        error: BinaryLiquidError,
    },
}
impl core::fmt::Display for BinaryMiscibilityMaterialError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound miscibility field: {self:?}")
    }
}
impl std::error::Error for BinaryMiscibilityMaterialError {}
impl From<BinaryMiscibilityError> for BinaryMiscibilityMaterialError {
    fn from(error: BinaryMiscibilityError) -> Self {
        Self::Solver(error)
    }
}
type Result<T> = core::result::Result<T, BinaryMiscibilityMaterialError>;

/// One accepted field and its original per-cell material evaluations.
/// Neither the source reference nor the numerical inventories are replaceable.
pub struct ResolvedBinaryMiscibilityField<'source, 'liquid, 'card> {
    source: &'source BinaryLiquidSquareGradientLaw<'liquid, 'card>,
    state: BinaryMiscibilityState,
    cells: Vec<ResolvedBinaryLiquidCell>,
    maximum_resolved_cells: usize,
}

impl core::fmt::Debug for ResolvedBinaryMiscibilityField<'_, '_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ResolvedBinaryMiscibilityField")
            .field("state", &self.state)
            .field("cells", &self.cells)
            .field("maximum_resolved_cells", &self.maximum_resolved_cells)
            .finish_non_exhaustive()
    }
}

impl BinaryMiscibilityNetwork {
    /// Bind an original graph state to the exact liquid/gradient source owner.
    /// The complete model (pair, temperature, concentration and kappa) must
    /// agree. Owner and count-budget checks precede state cloning and queries.
    /// Each source receives the original independently represented A/B amounts.
    /// Refusal or cancellation returns no partially resolved field.
    pub fn resolve_liquid_field<'source, 'liquid, 'card>(
        &self,
        state: &BinaryMiscibilityState,
        source: &'source BinaryLiquidSquareGradientLaw<'liquid, 'card>,
        maximum_resolved_cells: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedBinaryMiscibilityField<'source, 'liquid, 'card>> {
        poll(&mut cancelled)?;
        self.admit(state)?;
        if self.model() != source.gradient_model() {
            return Err(BinaryMiscibilityMaterialError::Invalid(
                "field and liquid source require the same original full gradient model",
            ));
        }
        if state.component_moles().len() > maximum_resolved_cells {
            return Err(BinaryMiscibilityMaterialError::CellBudget {
                required: state.component_moles().len(),
                maximum: maximum_resolved_cells,
            });
        }
        ResolvedBinaryMiscibilityField::resolve_owned(
            state.clone(),
            source,
            maximum_resolved_cells,
            cancelled,
        )
    }
}

impl<'source, 'liquid, 'card> ResolvedBinaryMiscibilityField<'source, 'liquid, 'card> {
    pub const fn state(&self) -> &BinaryMiscibilityState {
        &self.state
    }
    pub fn cells(&self) -> &[ResolvedBinaryLiquidCell] {
        &self.cells
    }
    pub const fn source(&self) -> &'source BinaryLiquidSquareGradientLaw<'liquid, 'card> {
        self.source
    }
    pub const fn maximum_resolved_cells(&self) -> usize {
        self.maximum_resolved_cells
    }

    /// Advance advection and molecular/gradient diffusion, then refresh sources.
    /// The same borrowed liquid binding resolves only the completed numerical
    /// state. Any final source-domain/volume failure or cancellation discards
    /// the whole candidate, preserving this field and all its original receipts.
    pub fn transport_candidate(
        &self,
        flow: &BinaryMiscibilityAdvection,
        duration_s: f64,
        advection_controls: BinaryMiscibilityAdvectionControls,
        diffusion_controls: BinaryMiscibilitySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityMaterialCandidate<'source, 'liquid, 'card>> {
        let network = BinaryMiscibilityNetwork {
            owner: Arc::clone(&self.state.owner),
        };
        let candidate = network.transport_candidate(
            &self.state,
            flow,
            duration_s,
            advection_controls,
            diffusion_controls,
            &mut cancelled,
        )?;
        let (state, step) = candidate.into_parts();
        let field =
            Self::resolve_owned(state, self.source, self.maximum_resolved_cells, cancelled)?;
        Ok(BinaryMiscibilityMaterialCandidate { field, step })
    }

    // Both callers have already admitted this exact owner and cell-count budget.
    // Consume the solved state to avoid cloning its numerical vectors again.
    fn resolve_owned(
        state: BinaryMiscibilityState,
        source: &'source BinaryLiquidSquareGradientLaw<'liquid, 'card>,
        maximum_resolved_cells: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        let mut cells = Vec::with_capacity(state.component_moles().len());
        for (cell, &amounts) in state.component_moles().iter().enumerate() {
            poll(&mut cancelled)?;
            cells.push(
                source
                    .resolve_cell(amounts, &mut cancelled)
                    .map_err(|error| BinaryMiscibilityMaterialError::Source { cell, error })?,
            );
        }
        poll(&mut cancelled)?;
        Ok(Self {
            source,
            state,
            cells,
            maximum_resolved_cells,
        })
    }
}

/// Complete conservative transport and refreshed original material evaluations.
#[derive(Debug)]
pub struct BinaryMiscibilityMaterialCandidate<'source, 'liquid, 'card> {
    field: ResolvedBinaryMiscibilityField<'source, 'liquid, 'card>,
    step: BinaryMiscibilityTransportStep,
}
impl<'source, 'liquid, 'card> BinaryMiscibilityMaterialCandidate<'source, 'liquid, 'card> {
    pub const fn field(&self) -> &ResolvedBinaryMiscibilityField<'source, 'liquid, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &BinaryMiscibilityTransportStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedBinaryMiscibilityField<'source, 'liquid, 'card>,
        BinaryMiscibilityTransportStep,
    ) {
        (self.field, self.step)
    }
}
