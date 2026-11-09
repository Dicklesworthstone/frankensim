//! Original liquid and caloric sources for local-temperature multicomponent interfaces.
//!
//! Every complete endpoint resolves its actual local temperature and independent
//! species inventory through the same original liquid law. The numerical owner
//! retains the original homogeneous calorics and energetic interface separately;
//! this adapter preserves its complete energy/entropy and transport receipts.
//! Spinodal cells remain homogeneous branches of the evolving diffuse interface.
//!
//! The original pressure, composition axes, selected claims, source-volume
//! allowance and cell budget persist. Material receipts qualify returned
//! endpoints, not every adaptive trial. Molecular pair data, conductivity and
//! heats of transport remain explicitly authored numerical coefficients;
//! resolved mixture Cp never replaces the original pure-liquid calorics.

use super::adaptive::{
    MulticomponentNonisothermalPhaseFieldAdvance, MulticomponentNonisothermalPhaseFieldTimeControl,
};
use super::{MulticomponentNonisothermalPhaseField, MulticomponentNonisothermalPhaseFieldStep};
use crate::miscibility::multicomponent::nonisothermal::{
    MulticomponentNonisothermalError, MulticomponentNonisothermalFace,
    MulticomponentNonisothermalSolve, material::resolve_local_cells,
};
use crate::miscibility::multicomponent::phase_field::{
    MulticomponentPhaseFieldNetwork, material::ResolvedMulticomponentPhaseField,
};
use fs_material::liquid::solution::multicomponent::{
    MulticomponentLiquidSolutionLaw, ResolvedMulticomponentLiquidCell,
    caloric::MulticomponentLiquidCaloricLaw,
};
use fs_qty::Temperature;
use std::sync::Arc;

/// Original numerical, source, volume and cancellation failures.
pub use crate::miscibility::multicomponent::nonisothermal::material::MulticomponentMiscibilityMaterialNonisothermalError
    as MulticomponentPhaseFieldMaterialNonisothermalError;

type Result<T> = core::result::Result<T, MulticomponentPhaseFieldMaterialNonisothermalError>;

/// Complete energetic thermal field with its exact original caloric/liquid pair.
///
/// The caloric law itself borrows the liquid source, keeping their identity
/// together across candidates. Only the energetic thermal owner is exposed;
/// its private ordinary inventory and thermal representations cannot escape.
pub struct ResolvedMulticomponentNonisothermalPhaseField<'caloric, 'liquid, 'card> {
    numerical: MulticomponentNonisothermalPhaseField,
    caloric: &'caloric MulticomponentLiquidCaloricLaw<'liquid, 'card>,
    cells: Arc<[ResolvedMulticomponentLiquidCell]>,
    maximum_resolved_cells: usize,
    maximum_relative_volume_mismatch: f64,
}

impl core::fmt::Debug for ResolvedMulticomponentNonisothermalPhaseField<'_, '_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ResolvedMulticomponentNonisothermalPhaseField")
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

impl<'source, 'card> ResolvedMulticomponentPhaseField<'source, 'card> {
    /// Initialize local temperatures from the exact original liquid's calorics.
    /// This initializes thermal state, not an energy-conserving heating event.
    ///
    /// The existing source field has already admitted the complete bulk law,
    /// canonical species/masses and fixed geometry. A different liquid instance
    /// refuses even when its values match. Every actual initialized endpoint
    /// must satisfy the original source and inherited density/volume allowance.
    pub fn bind_nonisothermal<'caloric, 'liquid>(
        &self,
        caloric: &'caloric MulticomponentLiquidCaloricLaw<'liquid, 'card>,
        temperatures: Vec<Temperature>,
        faces: Vec<MulticomponentNonisothermalFace>,
        temperature_bracket: [Temperature; 2],
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedMulticomponentNonisothermalPhaseField<'caloric, 'liquid, 'card>> {
        poll(&mut cancelled)?;
        if !std::ptr::eq(self.source(), caloric.liquid_law()) {
            return Err(MulticomponentPhaseFieldMaterialNonisothermalError::Invalid(
                "thermal phase field requires its exact original liquid and caloric binding",
            ));
        }
        let numerical = self.network().bind_nonisothermal_phase_field(
            self.field(),
            caloric.model().clone(),
            temperatures,
            faces,
            temperature_bracket,
            maximum_workspace_bytes,
            &mut cancelled,
        )?;
        ResolvedMulticomponentNonisothermalPhaseField::from_numerical(
            numerical,
            caloric,
            self.maximum_resolved_cells(),
            self.maximum_relative_volume_mismatch(),
            &mut cancelled,
        )
    }
}

impl<'caloric, 'liquid, 'card>
    ResolvedMulticomponentNonisothermalPhaseField<'caloric, 'liquid, 'card>
{
    pub const fn field(&self) -> &MulticomponentNonisothermalPhaseField {
        &self.numerical
    }
    /// The graph temperature remains a reference. Actual local temperatures
    /// belong to the thermal field and each original material receipt.
    pub fn network(&self) -> &MulticomponentPhaseFieldNetwork {
        self.numerical.network()
    }
    pub fn source(&self) -> &MulticomponentLiquidSolutionLaw<'card> {
        self.caloric.liquid_law()
    }
    pub const fn caloric_law(&self) -> &'caloric MulticomponentLiquidCaloricLaw<'liquid, 'card> {
        self.caloric
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

    /// Complete one coupled interfacial species/heat step, then qualify its
    /// actual local sources. Late source, volume or cancellation refusal is atomic.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: MulticomponentNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentPhaseFieldMaterialNonisothermalCandidate<'caloric, 'liquid, 'card>>
    {
        let (numerical, step) = self
            .numerical
            .step_candidate(duration_s, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(MulticomponentPhaseFieldMaterialNonisothermalCandidate { field, step })
    }

    /// Control composition, temperature and individual interface energies over
    /// the complete interval, then qualify its final original-source state.
    /// Numerical trial/error and physical energy/entropy receipts are preserved.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        time: MulticomponentNonisothermalPhaseFieldTimeControl,
        controls: MulticomponentNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        MulticomponentPhaseFieldMaterialNonisothermalAdvanceCandidate<'caloric, 'liquid, 'card>,
    > {
        let (numerical, advance) = self
            .numerical
            .advance_candidate(duration_s, time, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(MulticomponentPhaseFieldMaterialNonisothermalAdvanceCandidate { field, advance })
    }

    fn refresh(
        &self,
        numerical: MulticomponentNonisothermalPhaseField,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::from_numerical(
            numerical,
            self.caloric,
            self.maximum_resolved_cells,
            self.maximum_relative_volume_mismatch,
            cancelled,
        )
    }

    fn from_numerical(
        numerical: MulticomponentNonisothermalPhaseField,
        caloric: &'caloric MulticomponentLiquidCaloricLaw<'liquid, 'card>,
        maximum_resolved_cells: usize,
        maximum_relative_volume_mismatch: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let cells = resolve_local_cells(
            numerical.component_moles(),
            numerical.temperatures(),
            numerical.network().cells(),
            caloric.liquid_law(),
            maximum_resolved_cells,
            maximum_relative_volume_mismatch,
            cancelled,
        )?;
        Ok(Self {
            numerical,
            caloric,
            cells,
            maximum_resolved_cells,
            maximum_relative_volume_mismatch,
        })
    }
}

/// Complete source-qualified endpoint and unchanged interfacial thermal receipt.
#[derive(Debug)]
pub struct MulticomponentPhaseFieldMaterialNonisothermalCandidate<'caloric, 'liquid, 'card> {
    field: ResolvedMulticomponentNonisothermalPhaseField<'caloric, 'liquid, 'card>,
    step: MulticomponentNonisothermalPhaseFieldStep,
}

impl<'caloric, 'liquid, 'card>
    MulticomponentPhaseFieldMaterialNonisothermalCandidate<'caloric, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedMulticomponentNonisothermalPhaseField<'caloric, 'liquid, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentNonisothermalPhaseFieldStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentNonisothermalPhaseField<'caloric, 'liquid, 'card>,
        MulticomponentNonisothermalPhaseFieldStep,
    ) {
        (self.field, self.step)
    }
}

/// Complete adaptive demixing/heat interval and its original-source endpoint.
#[derive(Debug)]
pub struct MulticomponentPhaseFieldMaterialNonisothermalAdvanceCandidate<'caloric, 'liquid, 'card>
{
    field: ResolvedMulticomponentNonisothermalPhaseField<'caloric, 'liquid, 'card>,
    advance: MulticomponentNonisothermalPhaseFieldAdvance,
}

impl<'caloric, 'liquid, 'card>
    MulticomponentPhaseFieldMaterialNonisothermalAdvanceCandidate<'caloric, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedMulticomponentNonisothermalPhaseField<'caloric, 'liquid, 'card> {
        &self.field
    }
    pub const fn advance(&self) -> &MulticomponentNonisothermalPhaseFieldAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentNonisothermalPhaseField<'caloric, 'liquid, 'card>,
        MulticomponentNonisothermalPhaseFieldAdvance,
    ) {
        (self.field, self.advance)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentNonisothermalError::Cancelled.into())
    } else {
        Ok(())
    }
}
