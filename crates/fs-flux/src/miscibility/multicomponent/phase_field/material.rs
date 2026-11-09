//! Original homogeneous material sources through multicomponent demixing.
//!
//! The numerical phase-field owner retains all independent component amounts,
//! the original nonconvex bulk law and its explicit interfacial energy. Each
//! completed endpoint resolves those actual local inventories through the same
//! borrowed liquid source. A spinodal cell remains homogeneous at mesh scale;
//! this adapter never replaces it with a local bulk-equilibrium flash.
//!
//! Source composition axes, pressure and selected claims remain authoritative.
//! The source density must still agree with the fixed numerical cell volume.
//! Pair diffusivities and the square-gradient coefficient remain explicitly
//! authored numerical data. Material receipts qualify returned endpoints,
//! not every internal adaptive trial, and are bounded by a retained cell cap.

use super::adaptive::{MulticomponentPhaseFieldAdvance, MulticomponentPhaseFieldTimeControl};
use super::{MulticomponentPhaseField, MulticomponentPhaseFieldNetwork, MulticomponentPhaseFieldStep};
use super::super::MulticomponentMiscibilitySolve;
use super::super::material::{admit_source, resolve_cells};
use fs_material::liquid::solution::multicomponent::{
    MulticomponentLiquidSolutionLaw, ResolvedMulticomponentLiquidCell,
};
use std::sync::Arc;

/// The original material admission, source, conservation and work-limit errors.
pub use super::super::material::MulticomponentMiscibilityMaterialError
    as MulticomponentPhaseFieldMaterialError;

type Result<T> = core::result::Result<T, MulticomponentPhaseFieldMaterialError>;

/// Accepted interfacial field and material receipts from its exact original law.
///
/// Only the phase-field numerical owner is exposed. The private ordinary
/// inventory representation cannot escape to a zero-interface-energy solver.
pub struct ResolvedMulticomponentPhaseField<'source, 'card> {
    numerical: MulticomponentPhaseField,
    source: &'source MulticomponentLiquidSolutionLaw<'card>,
    cells: Arc<[ResolvedMulticomponentLiquidCell]>,
    maximum_resolved_cells: usize,
    maximum_relative_volume_mismatch: f64,
}

impl core::fmt::Debug for ResolvedMulticomponentPhaseField<'_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ResolvedMulticomponentPhaseField")
            .field("numerical", &self.numerical)
            .field("cells", &self.cells)
            .field("maximum_resolved_cells", &self.maximum_resolved_cells)
            .field(
                "maximum_relative_volume_mismatch",
                &self.maximum_relative_volume_mismatch,
            )
            .finish_non_exhaustive()
    }
}

impl MulticomponentPhaseField {
    /// Bind original mixture properties to an already admitted interfacial field.
    ///
    /// The complete original model, canonical species and molar masses must
    /// agree. Source admission and density/volume checks are the same as for
    /// ordinary multicomponent fields; no bulk-convexity gate is added here.
    /// The relative volume allowance is numerical, in [0,1e-6]. Source and
    /// resolved-cell budgets are checked before cloning numerical state.
    pub fn resolve_liquid_field<'source, 'card>(
        &self,
        source: &'source MulticomponentLiquidSolutionLaw<'card>,
        maximum_resolved_cells: usize,
        maximum_relative_volume_mismatch: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedMulticomponentPhaseField<'source, 'card>> {
        admit_source(
            self.bulk.network(),
            source,
            maximum_relative_volume_mismatch,
            &mut cancelled,
        )?;
        let cells = resolve_cells(
            &self.bulk,
            source,
            maximum_resolved_cells,
            maximum_relative_volume_mismatch,
            &mut cancelled,
        )?;
        let numerical = self.clone();
        poll(&mut cancelled)?;
        Ok(ResolvedMulticomponentPhaseField {
            numerical,
            source,
            cells,
            maximum_resolved_cells,
            maximum_relative_volume_mismatch,
        })
    }
}

impl<'source, 'card> ResolvedMulticomponentPhaseField<'source, 'card> {
    pub const fn field(&self) -> &MulticomponentPhaseField {
        &self.numerical
    }
    pub fn network(&self) -> &MulticomponentPhaseFieldNetwork {
        self.numerical.network()
    }
    pub const fn source(&self) -> &'source MulticomponentLiquidSolutionLaw<'card> {
        self.source
    }
    pub fn cells(&self) -> &[ResolvedMulticomponentLiquidCell] {
        &self.cells
    }
    pub const fn maximum_resolved_cells(&self) -> usize {
        self.maximum_resolved_cells
    }
    pub const fn maximum_relative_volume_mismatch(&self) -> f64 {
        self.maximum_relative_volume_mismatch
    }

    /// Complete one conservative bulk-plus-interface step, then resolve its
    /// actual local compositions. Source or volume refusal is atomic.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: MulticomponentMiscibilitySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentPhaseFieldMaterialCandidate<'source, 'card>> {
        let (numerical, step) = self
            .numerical
            .step_candidate(duration_s, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(MulticomponentPhaseFieldMaterialCandidate { field, step })
    }

    /// Qualify the complete adaptive endpoint through the retained original
    /// source. Composition/interface accuracy and all physical transfer and
    /// total-free-energy receipts remain those of the numerical owner.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        time: MulticomponentPhaseFieldTimeControl,
        controls: MulticomponentMiscibilitySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentPhaseFieldMaterialAdvanceCandidate<'source, 'card>> {
        let (numerical, advance) = self
            .numerical
            .advance_candidate(duration_s, time, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(MulticomponentPhaseFieldMaterialAdvanceCandidate { field, advance })
    }

    fn refresh(
        &self,
        numerical: MulticomponentPhaseField,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let cells = resolve_cells(
            &numerical.bulk,
            self.source,
            self.maximum_resolved_cells,
            self.maximum_relative_volume_mismatch,
            cancelled,
        )?;
        Ok(Self {
            numerical,
            source: self.source,
            cells,
            maximum_resolved_cells: self.maximum_resolved_cells,
            maximum_relative_volume_mismatch: self.maximum_relative_volume_mismatch,
        })
    }
}

/// One original-source endpoint and its unchanged interfacial step receipt.
#[derive(Debug)]
pub struct MulticomponentPhaseFieldMaterialCandidate<'source, 'card> {
    field: ResolvedMulticomponentPhaseField<'source, 'card>,
    step: MulticomponentPhaseFieldStep,
}

impl<'source, 'card> MulticomponentPhaseFieldMaterialCandidate<'source, 'card> {
    pub const fn field(&self) -> &ResolvedMulticomponentPhaseField<'source, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentPhaseFieldStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentPhaseField<'source, 'card>,
        MulticomponentPhaseFieldStep,
    ) {
        (self.field, self.step)
    }
}

/// Complete adaptive demixing interval and its original-source final state.
#[derive(Debug)]
pub struct MulticomponentPhaseFieldMaterialAdvanceCandidate<'source, 'card> {
    field: ResolvedMulticomponentPhaseField<'source, 'card>,
    advance: MulticomponentPhaseFieldAdvance,
}

impl<'source, 'card> MulticomponentPhaseFieldMaterialAdvanceCandidate<'source, 'card> {
    pub const fn field(&self) -> &ResolvedMulticomponentPhaseField<'source, 'card> {
        &self.field
    }
    pub const fn advance(&self) -> &MulticomponentPhaseFieldAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentPhaseField<'source, 'card>,
        MulticomponentPhaseFieldAdvance,
    ) {
        (self.field, self.advance)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentPhaseFieldMaterialError::Cancelled)
    } else {
        Ok(())
    }
}
