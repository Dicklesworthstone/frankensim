//! Original homogeneous liquid sources for coupled multicomponent transport.
//!
//! Each completed numerical endpoint resolves every actual species inventory
//! through the same borrowed liquid law at the network's fixed temperature.
//! Original pressure, composition-axis assignments and selected claims remain
//! authoritative. Cells retain the homogeneous branch admitted by the numerical
//! owner's convexity gate; this adapter never replaces them by local flashes.
//!
//! The original source density must agree with each fixed numerical cell volume
//! within an explicit retained allowance. Full Maxwell–Stefan pair matrices
//! remain the numerical owner's authored data: neither a binary-only source
//! nor a transport property is reinterpreted as multicomponent molecular data.
//! Material receipts qualify completed endpoints, not every internal adaptive
//! trial. Their allocations are separate from the numerical workspace and are
//! bounded by an explicit resolved-cell count.

use super::adaptive::{
    MulticomponentMiscibilityAdvance, MulticomponentMiscibilityAdvanceError,
    MulticomponentMiscibilityTimeControl,
};
use super::{
    MulticomponentMiscibilityError, MulticomponentMiscibilityField,
    MulticomponentMiscibilityNetwork, MulticomponentMiscibilitySolve,
    MulticomponentMiscibilityStep,
};
use fs_material::liquid::solution::multicomponent::{
    MulticomponentLiquidError, MulticomponentLiquidSolutionLaw, ResolvedMulticomponentLiquidCell,
};
use std::sync::Arc;

/// Original material receipts for prescribed stirring and coupled transport.
pub mod transport;

/// Original source incompatibility or an atomic numerical/material refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentMiscibilityMaterialError {
    Invalid(&'static str),
    CellBudget {
        required: usize,
        maximum: usize,
    },
    Solver(MulticomponentMiscibilityError),
    Advance(MulticomponentMiscibilityAdvanceError),
    Source {
        cell: usize,
        error: MulticomponentLiquidError,
    },
    VolumeMismatch {
        cell: usize,
        expected_m3: f64,
        actual_m3: f64,
        maximum_relative_mismatch: f64,
    },
    Cancelled,
}
impl core::fmt::Display for MulticomponentMiscibilityMaterialError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound multicomponent miscibility: {self:?}")
    }
}
impl std::error::Error for MulticomponentMiscibilityMaterialError {}
impl From<MulticomponentMiscibilityError> for MulticomponentMiscibilityMaterialError {
    fn from(error: MulticomponentMiscibilityError) -> Self {
        Self::Solver(error)
    }
}
impl From<MulticomponentMiscibilityAdvanceError> for MulticomponentMiscibilityMaterialError {
    fn from(error: MulticomponentMiscibilityAdvanceError) -> Self {
        Self::Advance(error)
    }
}
type Result<T> = core::result::Result<T, MulticomponentMiscibilityMaterialError>;

/// Accepted numerical field with its original homogeneous material receipts.
///
/// Binding compares the full thermodynamic model and canonical molecular
/// masses, allowing independently reconstructed identical numerical models.
/// The exact borrowed material law is then retained across every candidate.
pub struct ResolvedMulticomponentMiscibilityField<'source, 'card> {
    numerical: MulticomponentMiscibilityField,
    source: &'source MulticomponentLiquidSolutionLaw<'card>,
    cells: Arc<[ResolvedMulticomponentLiquidCell]>,
    maximum_resolved_cells: usize,
    maximum_relative_volume_mismatch: f64,
}
impl core::fmt::Debug for ResolvedMulticomponentMiscibilityField<'_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ResolvedMulticomponentMiscibilityField")
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

impl MulticomponentMiscibilityField {
    /// Bind the original mixture source to this already admitted numerical field.
    ///
    /// Canonical species, all thermodynamic coefficients and original molar
    /// masses must agree. The caller qualifies the common pressure and fixed
    /// equal partial molar volumes; actual source density is checked at every
    /// returned cell. The relative volume allowance is numerical, in [0,1e-6],
    /// and cannot be used to replace an incompatible density law. Admission and
    /// cell-budget checks precede source queries and numerical state cloning.
    pub fn resolve_liquid_field<'source, 'card>(
        &self,
        source: &'source MulticomponentLiquidSolutionLaw<'card>,
        maximum_resolved_cells: usize,
        maximum_relative_volume_mismatch: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedMulticomponentMiscibilityField<'source, 'card>> {
        poll(&mut cancelled)?;
        admit_source(
            self.network(),
            source,
            maximum_relative_volume_mismatch,
            &mut cancelled,
        )?;
        let cells = resolve_cells(
            self,
            source,
            maximum_resolved_cells,
            maximum_relative_volume_mismatch,
            &mut cancelled,
        )?;
        let numerical = self.clone();
        poll(&mut cancelled)?;
        Ok(ResolvedMulticomponentMiscibilityField {
            numerical,
            source,
            cells,
            maximum_resolved_cells,
            maximum_relative_volume_mismatch,
        })
    }
}

impl<'source, 'card> ResolvedMulticomponentMiscibilityField<'source, 'card> {
    pub const fn field(&self) -> &MulticomponentMiscibilityField {
        &self.numerical
    }
    pub fn network(&self) -> &MulticomponentMiscibilityNetwork {
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

    /// Complete one coupled numerical step, then resolve all actual endpoint
    /// compositions. A late source or volume refusal publishes no candidate.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: MulticomponentMiscibilitySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentMiscibilityMaterialCandidate<'source, 'card>> {
        let (numerical, step) = self
            .numerical
            .step_candidate(duration_s, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(MulticomponentMiscibilityMaterialCandidate { field, step })
    }

    /// Complete the adaptive interval before qualifying its final material
    /// state. Numerical trial/error and conservation receipts remain unchanged.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        time: MulticomponentMiscibilityTimeControl,
        controls: MulticomponentMiscibilitySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentMiscibilityMaterialAdvanceCandidate<'source, 'card>> {
        let (numerical, advance) = self
            .numerical
            .advance_candidate(duration_s, time, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(MulticomponentMiscibilityMaterialAdvanceCandidate { field, advance })
    }

    fn refresh(
        &self,
        numerical: MulticomponentMiscibilityField,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let cells = resolve_cells(
            &numerical,
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

/// One complete original-source endpoint and the unchanged numerical step.
#[derive(Debug)]
pub struct MulticomponentMiscibilityMaterialCandidate<'source, 'card> {
    field: ResolvedMulticomponentMiscibilityField<'source, 'card>,
    step: MulticomponentMiscibilityStep,
}
impl<'source, 'card> MulticomponentMiscibilityMaterialCandidate<'source, 'card> {
    pub const fn field(&self) -> &ResolvedMulticomponentMiscibilityField<'source, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentMiscibilityStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentMiscibilityField<'source, 'card>,
        MulticomponentMiscibilityStep,
    ) {
        (self.field, self.step)
    }
}

/// Complete requested interval with its original-source final composition.
#[derive(Debug)]
pub struct MulticomponentMiscibilityMaterialAdvanceCandidate<'source, 'card> {
    field: ResolvedMulticomponentMiscibilityField<'source, 'card>,
    advance: MulticomponentMiscibilityAdvance,
}
impl<'source, 'card> MulticomponentMiscibilityMaterialAdvanceCandidate<'source, 'card> {
    pub const fn field(&self) -> &ResolvedMulticomponentMiscibilityField<'source, 'card> {
        &self.field
    }
    pub const fn advance(&self) -> &MulticomponentMiscibilityAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentMiscibilityField<'source, 'card>,
        MulticomponentMiscibilityAdvance,
    ) {
        (self.field, self.advance)
    }
}

/// Shared original-source admission for homogeneous and interfacial owners.
/// The caller retains the appropriate numerical energy owner; this helper
/// checks material identity without imposing a homogeneous convexity claim.
pub(super) fn admit_source(
    network: &MulticomponentMiscibilityNetwork,
    source: &MulticomponentLiquidSolutionLaw<'_>,
    maximum_relative_volume_mismatch: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    poll(cancelled)?;
    if network.model() != source.thermodynamics()
        || network.transport().species() != source.species()
        || !network.transport().molar_masses().iter().copied().eq(source
            .components()
            .iter()
            .map(|component| component.molar_mass))
    {
        return Err(MulticomponentMiscibilityMaterialError::Invalid(
            "field and liquid source require the same full model, ordered species and molar masses",
        ));
    }
    if !maximum_relative_volume_mismatch.is_finite()
        || !(0.0..=1e-6).contains(&maximum_relative_volume_mismatch)
    {
        return Err(MulticomponentMiscibilityMaterialError::Invalid(
            "relative source-volume allowance must be in [0,1e-6]",
        ));
    }
    Ok(())
}

pub(super) fn resolve_cells(
    numerical: &MulticomponentMiscibilityField,
    source: &MulticomponentLiquidSolutionLaw<'_>,
    maximum_resolved_cells: usize,
    maximum_relative_volume_mismatch: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Arc<[ResolvedMulticomponentLiquidCell]>> {
    poll(cancelled)?;
    let count = numerical.network().cells().len();
    if count > maximum_resolved_cells {
        return Err(MulticomponentMiscibilityMaterialError::CellBudget {
            required: count,
            maximum: maximum_resolved_cells,
        });
    }
    let mut cells = Vec::with_capacity(count);
    for (cell, geometry) in numerical.network().cells().iter().enumerate() {
        poll(cancelled)?;
        let amounts = numerical
            .cell_component_moles(cell)
            .expect("admitted field retains every original cell inventory");
        let resolved = source
            .homogeneous(numerical.network().temperature(), amounts, &mut *cancelled)
            .map_err(|error| MulticomponentMiscibilityMaterialError::Source { cell, error })?;
        let actual_m3 = resolved.volume_m3();
        let expected_m3 = geometry.volume_m3;
        let mismatch = (actual_m3 - expected_m3).abs() / actual_m3.max(expected_m3);
        if !mismatch.is_finite() || mismatch > maximum_relative_volume_mismatch {
            return Err(MulticomponentMiscibilityMaterialError::VolumeMismatch {
                cell,
                expected_m3,
                actual_m3,
                maximum_relative_mismatch: maximum_relative_volume_mismatch,
            });
        }
        cells.push(resolved);
    }
    poll(cancelled)?;
    Ok(cells.into())
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentMiscibilityMaterialError::Cancelled)
    } else {
        Ok(())
    }
}
