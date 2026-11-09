//! Original-source local-temperature demixing with energetic interfaces.
//!
//! The numerical owner retains the original homogeneous calorics and the
//! original square-gradient energy separately. Each accepted endpoint resolves
//! its actual local temperature and composition through the same liquid source
//! and density/volume allowance. Cells remain homogeneous branches of the
//! evolving diffuse interface; this adapter never flashes them independently.
//!
//! Material receipts qualify completed endpoints only. They do not certify
//! every internal adaptive trial, replace the numerical energy/entropy audit,
//! or infer face conductivity from the resolved transport properties. The
//! authored thermal faces remain unchanged. Resolved mixture Cp remains a
//! transport property, never a substitute for the original pure-liquid calorics.
//! Material query allocations are separate from the numerical workspace and
//! bounded by the original cell cap.

use super::nonisothermal::{
    BinaryMiscibilityMaterialNonisothermalError as Error, resolve_local_cells,
};
use super::{
    BinaryLiquidSquareGradientLaw, BinaryMiscibilityNetwork, ResolvedBinaryLiquidCell,
    ResolvedBinaryMiscibilityField, SourceBinding, poll,
};
use crate::miscibility::nonisothermal::phase_field::{
    BinaryNonisothermalPhaseField, BinaryNonisothermalPhaseFieldStep,
    adaptive::{BinaryNonisothermalPhaseFieldAdvance, BinaryNonisothermalPhaseFieldTimeControl},
};
use crate::miscibility::nonisothermal::{BinaryNonisothermalFace, BinaryNonisothermalSolve};
use fs_material::liquid::solution::caloric::BinaryLiquidCaloricLaw;
use fs_qty::Temperature;
use std::sync::Arc;

type Result<T> = core::result::Result<T, Error>;

/// Accepted energetic phase field and its original local material evaluations.
///
/// The private numerical owner cannot be used as a zero-gradient convex field.
/// Original liquid/caloric bindings and the resolved-cell budget survive every
/// candidate, including fields bound after a prescribed temperature change.
pub struct ResolvedBinaryMiscibilityPhaseField<'caloric, 'source, 'liquid, 'card> {
    numerical: BinaryNonisothermalPhaseField,
    caloric: &'caloric BinaryLiquidCaloricLaw<'liquid, 'card>,
    source: SourceBinding<'source, 'liquid, 'card>,
    cells: Arc<[ResolvedBinaryLiquidCell]>,
    maximum_resolved_cells: usize,
}

impl core::fmt::Debug for ResolvedBinaryMiscibilityPhaseField<'_, '_, '_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ResolvedBinaryMiscibilityPhaseField")
            .field("numerical", &self.numerical)
            .field("cells", &self.cells)
            .field("maximum_resolved_cells", &self.maximum_resolved_cells)
            .finish_non_exhaustive()
    }
}

impl<'source, 'liquid, 'card> ResolvedBinaryMiscibilityField<'source, 'liquid, 'card> {
    /// Initialize local temperatures on the exact original liquid/caloric pair.
    /// This is thermal initialization, not an energy-conserving heating event.
    ///
    /// The original graph requires positive energetic kappa and common positive
    /// species support. Its pair law may be nonconvex. Every initialized cell
    /// must satisfy the original source domain and fixed-volume allowance before
    /// publication; inventories, elapsed time and material identities persist.
    pub fn bind_nonisothermal_phase_field<'caloric>(
        &self,
        caloric: &'caloric BinaryLiquidCaloricLaw<'liquid, 'card>,
        temperatures: Vec<Temperature>,
        faces: Vec<BinaryNonisothermalFace>,
        temperature_bracket: [Temperature; 2],
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedBinaryMiscibilityPhaseField<'caloric, 'source, 'liquid, 'card>> {
        poll(&mut cancelled)?;
        if !std::ptr::eq(self.source().liquid_law(), caloric.liquid_law()) {
            return Err(Error::Invalid(
                "phase field requires its exact original liquid and caloric binding",
            ));
        }
        let numerical = self.network().bind_nonisothermal_phase_field(
            &self.state,
            caloric.model().clone(),
            temperatures,
            faces,
            temperature_bracket,
            maximum_workspace_bytes,
            &mut cancelled,
        )?;
        let cells = resolve_local_cells(
            numerical.component_moles(),
            numerical.temperatures(),
            self.source(),
            self.maximum_resolved_cells,
            &mut cancelled,
        )?;
        Ok(ResolvedBinaryMiscibilityPhaseField {
            numerical,
            caloric,
            source: self.source.clone(),
            cells,
            maximum_resolved_cells: self.maximum_resolved_cells,
        })
    }
}

impl<'caloric, 'source, 'liquid, 'card>
    ResolvedBinaryMiscibilityPhaseField<'caloric, 'source, 'liquid, 'card>
{
    /// Original energetic-interface owner at its actual local temperatures.
    pub const fn field(&self) -> &BinaryNonisothermalPhaseField {
        &self.numerical
    }
    /// Original graph. Its model temperature remains a reference; actual
    /// local temperatures belong to `field().temperatures()` and `cells()`.
    pub fn network(&self) -> &BinaryMiscibilityNetwork {
        self.numerical.network()
    }
    pub fn cells(&self) -> &[ResolvedBinaryLiquidCell] {
        &self.cells
    }
    /// Original source wrapper, including its inherited density/volume allowance.
    pub fn source(&self) -> &BinaryLiquidSquareGradientLaw<'liquid, 'card> {
        self.source.law()
    }
    pub const fn caloric_law(&self) -> &'caloric BinaryLiquidCaloricLaw<'liquid, 'card> {
        self.caloric
    }
    pub const fn maximum_resolved_cells(&self) -> usize {
        self.maximum_resolved_cells
    }

    /// Complete one demixing/heat sweep, then refresh original endpoint sources.
    /// A late source, volume or cancellation refusal publishes no candidate.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: BinaryNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityMaterialPhaseFieldCandidate<'caloric, 'source, 'liquid, 'card>>
    {
        let (numerical, step) = self
            .numerical
            .step_candidate(duration_s, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(BinaryMiscibilityMaterialPhaseFieldCandidate { field, step })
    }

    /// Control composition, local temperature and individual interface energies
    /// over a complete interval, then qualify only its final source evaluations.
    /// The original numerical transfer and conservation receipts are preserved.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        time: BinaryNonisothermalPhaseFieldTimeControl,
        controls: BinaryNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        BinaryMiscibilityMaterialPhaseFieldAdvanceCandidate<'caloric, 'source, 'liquid, 'card>,
    > {
        let (numerical, advance) = self
            .numerical
            .advance_candidate(duration_s, time, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(BinaryMiscibilityMaterialPhaseFieldAdvanceCandidate { field, advance })
    }

    fn refresh(
        &self,
        numerical: BinaryNonisothermalPhaseField,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let cells = resolve_local_cells(
            numerical.component_moles(),
            numerical.temperatures(),
            self.source(),
            self.maximum_resolved_cells,
            cancelled,
        )?;
        Ok(Self {
            numerical,
            caloric: self.caloric,
            source: self.source.clone(),
            cells,
            maximum_resolved_cells: self.maximum_resolved_cells,
        })
    }
}

/// Complete audited phase-field sweep and its qualified endpoint materials.
#[derive(Debug)]
pub struct BinaryMiscibilityMaterialPhaseFieldCandidate<'caloric, 'source, 'liquid, 'card> {
    field: ResolvedBinaryMiscibilityPhaseField<'caloric, 'source, 'liquid, 'card>,
    step: BinaryNonisothermalPhaseFieldStep,
}
impl<'caloric, 'source, 'liquid, 'card>
    BinaryMiscibilityMaterialPhaseFieldCandidate<'caloric, 'source, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedBinaryMiscibilityPhaseField<'caloric, 'source, 'liquid, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &BinaryNonisothermalPhaseFieldStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedBinaryMiscibilityPhaseField<'caloric, 'source, 'liquid, 'card>,
        BinaryNonisothermalPhaseFieldStep,
    ) {
        (self.field, self.step)
    }
}

/// Complete adaptive phase-field interval and its qualified endpoint materials.
#[derive(Debug)]
pub struct BinaryMiscibilityMaterialPhaseFieldAdvanceCandidate<'caloric, 'source, 'liquid, 'card> {
    field: ResolvedBinaryMiscibilityPhaseField<'caloric, 'source, 'liquid, 'card>,
    advance: BinaryNonisothermalPhaseFieldAdvance,
}
impl<'caloric, 'source, 'liquid, 'card>
    BinaryMiscibilityMaterialPhaseFieldAdvanceCandidate<'caloric, 'source, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedBinaryMiscibilityPhaseField<'caloric, 'source, 'liquid, 'card> {
        &self.field
    }
    pub const fn advance(&self) -> &BinaryNonisothermalPhaseFieldAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedBinaryMiscibilityPhaseField<'caloric, 'source, 'liquid, 'card>,
        BinaryNonisothermalPhaseFieldAdvance,
    ) {
        (self.field, self.advance)
    }
}
