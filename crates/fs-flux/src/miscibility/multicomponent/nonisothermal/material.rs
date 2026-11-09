//! Original liquid and pure-caloric sources at every actual local temperature.
//!
//! The same original caloric law supplies all chemical enthalpy/entropy
//! references. Only a complete numerical endpoint refreshes material properties
//! at its actual temperature and independent species amounts. The original
//! pressure, composition axes, selected claims and density/volume allowance
//! persist. No cell is independently flashed or assigned a replacement density.
//!
//! Mixture transport Cp remains a separate property; it never substitutes for
//! the original pure-liquid calorics. Face conductivity and heats of transport
//! remain explicitly authored numerical data. Material receipts qualify final
//! endpoints, not every adaptive trial, and their allocations are bounded by
//! the inherited cell-count cap separately from numerical workspace admission.

use super::adaptive::{
    MulticomponentNonisothermalAdvance, MulticomponentNonisothermalAdvanceError,
    MulticomponentNonisothermalTimeControl,
};
use super::advection::{
    MulticomponentNonisothermalAdvectionControls, MulticomponentNonisothermalAdvectionStep,
};
use super::transport::{
    MulticomponentNonisothermalTransportStep, adaptive::MulticomponentNonisothermalTransportAdvance,
};
use super::{
    MulticomponentNonisothermalError, MulticomponentNonisothermalFace,
    MulticomponentNonisothermalField, MulticomponentNonisothermalSolve,
    MulticomponentNonisothermalStep,
};
use crate::miscibility::multicomponent::MulticomponentMiscibilityNetwork;
use crate::miscibility::multicomponent::advection::MulticomponentMiscibilityAdvection;
use crate::miscibility::multicomponent::material::{
    MulticomponentMiscibilityMaterialError as MaterialError, ResolvedMulticomponentMiscibilityField,
};
use fs_material::liquid::solution::multicomponent::{
    MulticomponentLiquidSolutionLaw, ResolvedMulticomponentLiquidCell,
    caloric::MulticomponentLiquidCaloricLaw,
};
use fs_qty::Temperature;
use std::sync::Arc;

/// Original-source incompatibility or an atomic numerical/material refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentMiscibilityMaterialNonisothermalError {
    Invalid(&'static str),
    Numerical(MulticomponentNonisothermalError),
    Advance(MulticomponentNonisothermalAdvanceError),
    Material(MaterialError),
}
impl core::fmt::Display for MulticomponentMiscibilityMaterialNonisothermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "source-bound nonisothermal multicomponent miscibility: {self:?}"
        )
    }
}
impl std::error::Error for MulticomponentMiscibilityMaterialNonisothermalError {}
impl From<MulticomponentNonisothermalError>
    for MulticomponentMiscibilityMaterialNonisothermalError
{
    fn from(error: MulticomponentNonisothermalError) -> Self {
        Self::Numerical(error)
    }
}
impl From<MulticomponentNonisothermalAdvanceError>
    for MulticomponentMiscibilityMaterialNonisothermalError
{
    fn from(error: MulticomponentNonisothermalAdvanceError) -> Self {
        Self::Advance(error)
    }
}
impl From<MaterialError> for MulticomponentMiscibilityMaterialNonisothermalError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
type Result<T> = core::result::Result<T, MulticomponentMiscibilityMaterialNonisothermalError>;

/// Complete numerical thermal state and its exact original source binding.
/// The retained caloric law itself borrows the original liquid law, so their
/// identity cannot drift independently across numerical candidates.
pub struct ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card> {
    numerical: MulticomponentNonisothermalField,
    caloric: &'caloric MulticomponentLiquidCaloricLaw<'liquid, 'card>,
    cells: Arc<[ResolvedMulticomponentLiquidCell]>,
    maximum_resolved_cells: usize,
    maximum_relative_volume_mismatch: f64,
}
impl core::fmt::Debug for ResolvedMulticomponentMiscibilityNonisothermalField<'_, '_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ResolvedMulticomponentMiscibilityNonisothermalField")
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

impl<'source, 'card> ResolvedMulticomponentMiscibilityField<'source, 'card> {
    /// Initialize local temperatures from the exact original liquid's calorics.
    /// This initializes thermal state, rather than applying a conserved heat
    /// transaction. The existing source field has already admitted the complete
    /// model, canonical species/molar masses and original fixed geometry.
    ///
    /// A different liquid-law instance refuses even if its values are identical.
    /// The numerical owner admits the complete caloric bracket; every actual
    /// initialized endpoint must also satisfy the original material sources and
    /// inherited source-volume allowance before publication.
    pub fn bind_nonisothermal<'caloric, 'liquid>(
        &self,
        caloric: &'caloric MulticomponentLiquidCaloricLaw<'liquid, 'card>,
        temperatures: Vec<Temperature>,
        faces: Vec<MulticomponentNonisothermalFace>,
        temperature_bracket: [Temperature; 2],
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card>> {
        poll(&mut cancelled)?;
        if !std::ptr::eq(self.source(), caloric.liquid_law()) {
            return Err(
                MulticomponentMiscibilityMaterialNonisothermalError::Invalid(
                    "thermal field requires its exact original liquid and caloric source binding",
                ),
            );
        }
        let numerical = self.network().bind_nonisothermal_field(
            self.field(),
            caloric.model().clone(),
            temperatures,
            faces,
            temperature_bracket,
            maximum_workspace_bytes,
            &mut cancelled,
        )?;
        ResolvedMulticomponentMiscibilityNonisothermalField::from_numerical(
            numerical,
            caloric,
            self.maximum_resolved_cells(),
            self.maximum_relative_volume_mismatch(),
            &mut cancelled,
        )
    }
}

impl<'caloric, 'liquid, 'card>
    ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card>
{
    pub const fn field(&self) -> &MulticomponentNonisothermalField {
        &self.numerical
    }
    /// The graph's temperature is a reference. Actual local temperatures belong
    /// to field().temperatures() and the original per-cell material receipts.
    pub fn network(&self) -> &MulticomponentMiscibilityNetwork {
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

    /// Complete the original heat/species step, then resolve each accepted
    /// local source. A late domain, volume or cancellation refusal is atomic.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: MulticomponentNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentMiscibilityMaterialNonisothermalCandidate<'caloric, 'liquid, 'card>>
    {
        let (numerical, step) = self
            .numerical
            .step_candidate(duration_s, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(MulticomponentMiscibilityMaterialNonisothermalCandidate { field, step })
    }

    /// Qualify only the complete adaptive endpoint through its original liquid
    /// sources, preserving all numerical transfer/energy/entropy receipts.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        time: MulticomponentNonisothermalTimeControl,
        controls: MulticomponentNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        MulticomponentMiscibilityMaterialNonisothermalAdvanceCandidate<'caloric, 'liquid, 'card>,
    > {
        let (numerical, advance) = self
            .numerical
            .advance_candidate(duration_s, time, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(MulticomponentMiscibilityMaterialNonisothermalAdvanceCandidate { field, advance })
    }

    /// Complete prescribed stirring of every species and its original donor
    /// enthalpy, then qualify the actual endpoint through the retained sources.
    pub fn advect_candidate(
        &self,
        flow: &MulticomponentMiscibilityAdvection,
        duration_s: f64,
        advection: MulticomponentNonisothermalAdvectionControls,
        controls: MulticomponentNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        MulticomponentMiscibilityMaterialNonisothermalAdvectionCandidate<'caloric, 'liquid, 'card>,
    > {
        let (numerical, step) = self
            .numerical
            .advect_candidate(flow, duration_s, advection, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(MulticomponentMiscibilityMaterialNonisothermalAdvectionCandidate { field, step })
    }

    /// Complete stirring and molecular/heat exchange over the same physical
    /// interval. Only the final endpoint receives original material receipts.
    pub fn transport_candidate(
        &self,
        flow: &MulticomponentMiscibilityAdvection,
        duration_s: f64,
        advection: MulticomponentNonisothermalAdvectionControls,
        controls: MulticomponentNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        MulticomponentMiscibilityMaterialNonisothermalTransportCandidate<'caloric, 'liquid, 'card>,
    > {
        let (numerical, step) = self
            .numerical
            .transport_candidate(flow, duration_s, advection, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(MulticomponentMiscibilityMaterialNonisothermalTransportCandidate { field, step })
    }

    /// Complete the bounded adaptive stirring/mixing interval and preserve its
    /// accepted species/energy receipts before final source qualification.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_transport_candidate(
        &self,
        flow: &MulticomponentMiscibilityAdvection,
        duration_s: f64,
        advection: MulticomponentNonisothermalAdvectionControls,
        controls: MulticomponentNonisothermalSolve,
        time: MulticomponentNonisothermalTimeControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        MulticomponentMiscibilityMaterialNonisothermalTransportAdvanceCandidate<
            'caloric,
            'liquid,
            'card,
        >,
    > {
        let (numerical, advance) = self
            .numerical
            .advance_transport_candidate(
                flow,
                duration_s,
                advection,
                controls,
                time,
                &mut cancelled,
            )?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(
            MulticomponentMiscibilityMaterialNonisothermalTransportAdvanceCandidate {
                field,
                advance,
            },
        )
    }

    fn refresh(
        &self,
        numerical: MulticomponentNonisothermalField,
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
        numerical: MulticomponentNonisothermalField,
        caloric: &'caloric MulticomponentLiquidCaloricLaw<'liquid, 'card>,
        maximum_resolved_cells: usize,
        maximum_relative_volume_mismatch: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let count = numerical.network().cells().len();
        if count > maximum_resolved_cells {
            return Err(MaterialError::CellBudget {
                required: count,
                maximum: maximum_resolved_cells,
            }
            .into());
        }
        if numerical.temperatures().len() != count {
            return Err(
                MulticomponentMiscibilityMaterialNonisothermalError::Invalid(
                    "one actual temperature per original material cell required",
                ),
            );
        }
        let mut cells = Vec::with_capacity(count);
        for (cell, (geometry, &temperature)) in numerical
            .network()
            .cells()
            .iter()
            .zip(numerical.temperatures())
            .enumerate()
        {
            poll(cancelled)?;
            let amounts = numerical
                .cell_component_moles(cell)
                .expect("admitted thermal field retains every original cell inventory");
            let resolved = caloric
                .liquid_law()
                .homogeneous(temperature, amounts, &mut *cancelled)
                .map_err(|error| MaterialError::Source { cell, error })?;
            let expected_m3 = geometry.volume_m3;
            let actual_m3 = resolved.volume_m3();
            let mismatch = (actual_m3 - expected_m3).abs() / actual_m3.max(expected_m3);
            if !mismatch.is_finite() || mismatch > maximum_relative_volume_mismatch {
                return Err(MaterialError::VolumeMismatch {
                    cell,
                    expected_m3,
                    actual_m3,
                    maximum_relative_mismatch: maximum_relative_volume_mismatch,
                }
                .into());
            }
            cells.push(resolved);
        }
        let cells = cells.into();
        poll(cancelled)?;
        Ok(Self {
            numerical,
            caloric,
            cells,
            maximum_resolved_cells,
            maximum_relative_volume_mismatch,
        })
    }
}

/// Original-source endpoint and the unchanged numerical step receipt.
#[derive(Debug)]
pub struct MulticomponentMiscibilityMaterialNonisothermalCandidate<'caloric, 'liquid, 'card> {
    field: ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card>,
    step: MulticomponentNonisothermalStep,
}
impl<'caloric, 'liquid, 'card>
    MulticomponentMiscibilityMaterialNonisothermalCandidate<'caloric, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentNonisothermalStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card>,
        MulticomponentNonisothermalStep,
    ) {
        (self.field, self.step)
    }
}

/// Complete adaptive endpoint and the unchanged numerical advance receipt.
#[derive(Debug)]
pub struct MulticomponentMiscibilityMaterialNonisothermalAdvanceCandidate<'caloric, 'liquid, 'card>
{
    field: ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card>,
    advance: MulticomponentNonisothermalAdvance,
}
impl<'caloric, 'liquid, 'card>
    MulticomponentMiscibilityMaterialNonisothermalAdvanceCandidate<'caloric, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card> {
        &self.field
    }
    pub const fn advance(&self) -> &MulticomponentNonisothermalAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card>,
        MulticomponentNonisothermalAdvance,
    ) {
        (self.field, self.advance)
    }
}

/// Complete original-enthalpy stirring and its source-qualified endpoint.
#[derive(Debug)]
pub struct MulticomponentMiscibilityMaterialNonisothermalAdvectionCandidate<
    'caloric,
    'liquid,
    'card,
> {
    field: ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card>,
    step: MulticomponentNonisothermalAdvectionStep,
}
impl<'caloric, 'liquid, 'card>
    MulticomponentMiscibilityMaterialNonisothermalAdvectionCandidate<'caloric, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentNonisothermalAdvectionStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card>,
        MulticomponentNonisothermalAdvectionStep,
    ) {
        (self.field, self.step)
    }
}

/// Complete stirring/molecular step and its source-qualified local temperatures.
#[derive(Debug)]
pub struct MulticomponentMiscibilityMaterialNonisothermalTransportCandidate<
    'caloric,
    'liquid,
    'card,
> {
    field: ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card>,
    step: MulticomponentNonisothermalTransportStep,
}
impl<'caloric, 'liquid, 'card>
    MulticomponentMiscibilityMaterialNonisothermalTransportCandidate<'caloric, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentNonisothermalTransportStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card>,
        MulticomponentNonisothermalTransportStep,
    ) {
        (self.field, self.step)
    }
}

/// Complete adaptive stirring interval and unchanged accepted physical receipts.
#[derive(Debug)]
pub struct MulticomponentMiscibilityMaterialNonisothermalTransportAdvanceCandidate<
    'caloric,
    'liquid,
    'card,
> {
    field: ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card>,
    advance: MulticomponentNonisothermalTransportAdvance,
}
impl<'caloric, 'liquid, 'card>
    MulticomponentMiscibilityMaterialNonisothermalTransportAdvanceCandidate<
        'caloric,
        'liquid,
        'card,
    >
{
    pub const fn field(
        &self,
    ) -> &ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card> {
        &self.field
    }
    pub const fn advance(&self) -> &MulticomponentNonisothermalTransportAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentMiscibilityNonisothermalField<'caloric, 'liquid, 'card>,
        MulticomponentNonisothermalTransportAdvance,
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
