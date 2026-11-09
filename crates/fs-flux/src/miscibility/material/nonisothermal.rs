//! Original liquid material sources for a locally nonisothermal mixture.
//!
//! The numerical owner advances the same homogeneous species and caloric law;
//! each accepted cell then resolves its actual temperature and composition
//! through the original liquid source. The source's existing fixed-volume
//! allowance still applies independently to every cell. No cell is flashed.
//!
//! The exact caloric and liquid bindings survive every candidate. Mixture Cp
//! remains an independent transport property, never a second caloric source.
//! Face conductivity and Soret heat remain explicitly authored numerical data:
//! resolving a local conductivity does not silently select a face mixing law.
//! Material queries and receipts allocate outside the numerical workspace and
//! are bounded by the original maximum resolved-cell count. A late domain,
//! volume or cancellation refusal publishes neither a field nor a partial report.

use super::{
    BinaryLiquidSquareGradientLaw, BinaryMiscibilityError, BinaryMiscibilityMaterialError,
    BinaryMiscibilityNetwork, ResolvedBinaryLiquidCell, ResolvedBinaryMiscibilityField,
    SourceBinding, poll,
};
use crate::miscibility::advection::{
    BinaryMiscibilityAdvection, BinaryMiscibilityAdvectionControls,
};
use crate::miscibility::nonisothermal::{
    BinaryNonisothermalError, BinaryNonisothermalFace, BinaryNonisothermalField,
    BinaryNonisothermalSolve, BinaryNonisothermalStep,
    adaptive::{
        BinaryNonisothermalAdvance, BinaryNonisothermalAdvanceError, BinaryNonisothermalTimeControl,
    },
    advection::BinaryNonisothermalAdvectionStep,
    transport::{
        BinaryNonisothermalTransportStep,
        adaptive::BinaryNonisothermalTransportAdvance,
    },
};
use fs_material::liquid::solution::caloric::BinaryLiquidCaloricLaw;
use fs_qty::Temperature;
use std::sync::Arc;

/// An original-source mismatch or an atomic numerical/material refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum BinaryMiscibilityMaterialNonisothermalError {
    Invalid(&'static str),
    Numerical(BinaryNonisothermalError),
    Advance(BinaryNonisothermalAdvanceError),
    Material(BinaryMiscibilityMaterialError),
}
impl core::fmt::Display for BinaryMiscibilityMaterialNonisothermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound nonisothermal miscibility: {self:?}")
    }
}
impl std::error::Error for BinaryMiscibilityMaterialNonisothermalError {}
impl From<BinaryNonisothermalError> for BinaryMiscibilityMaterialNonisothermalError {
    fn from(error: BinaryNonisothermalError) -> Self {
        Self::Numerical(error)
    }
}
impl From<BinaryNonisothermalAdvanceError> for BinaryMiscibilityMaterialNonisothermalError {
    fn from(error: BinaryNonisothermalAdvanceError) -> Self {
        Self::Advance(error)
    }
}
impl From<BinaryMiscibilityMaterialError> for BinaryMiscibilityMaterialNonisothermalError {
    fn from(error: BinaryMiscibilityMaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<BinaryMiscibilityError> for BinaryMiscibilityMaterialNonisothermalError {
    fn from(error: BinaryMiscibilityError) -> Self {
        Self::Material(error.into())
    }
}
type Result<T> = core::result::Result<T, BinaryMiscibilityMaterialNonisothermalError>;

/// One accepted numerical field and its original local material evaluations.
/// The field privately retains its exact caloric and liquid-source bindings.
pub struct ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card> {
    numerical: BinaryNonisothermalField,
    caloric: &'caloric BinaryLiquidCaloricLaw<'liquid, 'card>,
    source: SourceBinding<'source, 'liquid, 'card>,
    cells: Arc<[ResolvedBinaryLiquidCell]>,
    maximum_resolved_cells: usize,
}
impl core::fmt::Debug for ResolvedBinaryMiscibilityNonisothermalField<'_, '_, '_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ResolvedBinaryMiscibilityNonisothermalField")
            .field("numerical", &self.numerical)
            .field("cells", &self.cells)
            .field("maximum_resolved_cells", &self.maximum_resolved_cells)
            .finish_non_exhaustive()
    }
}

impl<'source, 'liquid, 'card> ResolvedBinaryMiscibilityField<'source, 'liquid, 'card> {
    /// Assign local temperatures using the exact original liquid's pure calorics.
    /// This initializes thermal state; it is not an energy-conserving heat event.
    /// The original numerical owner enforces kappa=0, convexity and positive
    /// common species support over the complete temperature bracket.
    ///
    /// Each local material source is refreshed before publication, including
    /// its density/volume compatibility. Geometry, paired inventories, elapsed
    /// time and the original maximum resolved-cell count remain unchanged.
    pub fn bind_nonisothermal<'caloric>(
        &self,
        caloric: &'caloric BinaryLiquidCaloricLaw<'liquid, 'card>,
        temperatures: Vec<Temperature>,
        faces: Vec<BinaryNonisothermalFace>,
        temperature_bracket: [Temperature; 2],
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card>> {
        poll(&mut cancelled)?;
        if !std::ptr::eq(self.source().liquid_law(), caloric.liquid_law()) {
            return Err(BinaryMiscibilityMaterialNonisothermalError::Invalid(
                "local thermal field requires its exact original liquid and caloric binding",
            ));
        }
        let numerical = self.network().bind_nonisothermal_field(
            &self.state,
            caloric.model().clone(),
            temperatures,
            faces,
            temperature_bracket,
            maximum_workspace_bytes,
            &mut cancelled,
        )?;
        let cells = resolve_local_cells(
            &numerical,
            self.source(),
            self.maximum_resolved_cells,
            &mut cancelled,
        )?;
        Ok(ResolvedBinaryMiscibilityNonisothermalField {
            numerical,
            caloric,
            source: self.source.clone(),
            cells,
            maximum_resolved_cells: self.maximum_resolved_cells,
        })
    }
}

impl<'caloric, 'source, 'liquid, 'card>
    ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card>
{
    /// Original numerical field, including each cell's actual temperature.
    pub const fn field(&self) -> &BinaryNonisothermalField {
        &self.numerical
    }
    /// The network's temperature is a reference; use field().temperatures()
    /// or cells() for actual local thermal states.
    pub fn network(&self) -> &BinaryMiscibilityNetwork {
        self.numerical.network()
    }
    pub fn cells(&self) -> &[ResolvedBinaryLiquidCell] {
        &self.cells
    }
    /// Original reference wrapper. Actual local sources are resolved at the
    /// corresponding field temperature, retaining this wrapper's volume allowance.
    pub fn source(&self) -> &BinaryLiquidSquareGradientLaw<'liquid, 'card> {
        self.source.law()
    }
    pub const fn caloric_law(&self) -> &'caloric BinaryLiquidCaloricLaw<'liquid, 'card> {
        self.caloric
    }
    pub const fn maximum_resolved_cells(&self) -> usize {
        self.maximum_resolved_cells
    }

    /// Advance one complete heat/species sweep, then qualify all final cells
    /// through their original material claims. Any late refusal is atomic.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: BinaryNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityMaterialNonisothermalCandidate<'caloric, 'source, 'liquid, 'card>>
    {
        let (numerical, step) = self
            .numerical
            .step_candidate(duration_s, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(BinaryMiscibilityMaterialNonisothermalCandidate { field, step })
    }

    /// Advance an accuracy-controlled interval using the original temperature
    /// and composition estimator. Only the complete numerical endpoint is
    /// materialized: this is endpoint source qualification, not proof that the
    /// independently sourced transport properties support every internal trial.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        time: BinaryNonisothermalTimeControl,
        controls: BinaryNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityMaterialNonisothermalAdvanceCandidate<'caloric, 'source, 'liquid, 'card>>
    {
        let (numerical, advance) = self
            .numerical
            .advance_candidate(duration_s, time, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(BinaryMiscibilityMaterialNonisothermalAdvanceCandidate { field, advance })
    }

    /// Advect original component amounts and parcel enthalpy together, then
    /// refresh local material properties at the recovered temperatures.
    pub fn advect_candidate(
        &self,
        flow: &BinaryMiscibilityAdvection,
        duration_s: f64,
        advection_controls: BinaryMiscibilityAdvectionControls,
        controls: BinaryNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityMaterialNonisothermalAdvectionCandidate<'caloric, 'source, 'liquid, 'card>>
    {
        let (numerical, step) = self
            .numerical
            .advect_candidate(
                flow,
                duration_s,
                advection_controls,
                controls,
                &mut cancelled,
            )?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(BinaryMiscibilityMaterialNonisothermalAdvectionCandidate { field, step })
    }

    /// Compose heat/species diffusion and stirring over one physical interval.
    /// The complete numerical candidate precedes the original material refresh;
    /// a failed final source query rejects every stage and its report together.
    pub fn transport_candidate(
        &self,
        flow: &BinaryMiscibilityAdvection,
        duration_s: f64,
        advection_controls: BinaryMiscibilityAdvectionControls,
        controls: BinaryNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityMaterialNonisothermalTransportCandidate<'caloric, 'source, 'liquid, 'card>>
    {
        let (numerical, step) = self
            .numerical
            .transport_candidate(
                flow,
                duration_s,
                advection_controls,
                controls,
                &mut cancelled,
            )?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(BinaryMiscibilityMaterialNonisothermalTransportCandidate { field, step })
    }

    /// Control the complete stirring/heat/species interval using both local
    /// temperature and composition, then resolve every accepted endpoint cell.
    /// Material receipts do not replace the numerical owner or its estimator.
    pub fn advance_transport_candidate(
        &self,
        flow: &BinaryMiscibilityAdvection,
        duration_s: f64,
        advection_controls: BinaryMiscibilityAdvectionControls,
        controls: BinaryNonisothermalSolve,
        time: BinaryNonisothermalTimeControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityMaterialNonisothermalTransportAdvanceCandidate<'caloric, 'source, 'liquid, 'card>>
    {
        let (numerical, advance) = self
            .numerical
            .advance_transport_candidate(
                flow,
                duration_s,
                advection_controls,
                controls,
                time,
                &mut cancelled,
            )?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        Ok(BinaryMiscibilityMaterialNonisothermalTransportAdvanceCandidate { field, advance })
    }

    fn refresh(
        &self,
        numerical: BinaryNonisothermalField,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let cells = resolve_local_cells(
            &numerical,
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

fn resolve_local_cells(
    numerical: &BinaryNonisothermalField,
    source: &BinaryLiquidSquareGradientLaw<'_, '_>,
    maximum_resolved_cells: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Arc<[ResolvedBinaryLiquidCell]>> {
    poll(cancelled)?;
    let count = numerical.component_moles().len();
    if count > maximum_resolved_cells {
        return Err(BinaryMiscibilityMaterialError::CellBudget {
            required: count,
            maximum: maximum_resolved_cells,
        }
        .into());
    }
    let mut cells = Vec::with_capacity(count);
    for (cell, (&amounts, &temperature)) in numerical
        .component_moles()
        .iter()
        .zip(numerical.temperatures())
        .enumerate()
    {
        poll(cancelled)?;
        // Reuse the original volume allowance and homogeneous source resolver;
        // neither an additive-density rule nor a local equilibrium is inferred.
        let local = source
            .at_temperature(temperature)
            .and_then(|local| local.resolve_cell(amounts, &mut *cancelled))
            .map_err(|error| BinaryMiscibilityMaterialError::Source { cell, error })?;
        cells.push(local);
    }
    poll(cancelled)?;
    Ok(cells.into())
}

/// Complete numerical sweep with all actual endpoint materials qualified.
#[derive(Debug)]
pub struct BinaryMiscibilityMaterialNonisothermalCandidate<'caloric, 'source, 'liquid, 'card> {
    field: ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card>,
    step: BinaryNonisothermalStep,
}
impl<'caloric, 'source, 'liquid, 'card>
    BinaryMiscibilityMaterialNonisothermalCandidate<'caloric, 'source, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &BinaryNonisothermalStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card>,
        BinaryNonisothermalStep,
    ) {
        (self.field, self.step)
    }
}

/// Complete adaptive interval with all actual endpoint materials qualified.
#[derive(Debug)]
pub struct BinaryMiscibilityMaterialNonisothermalAdvanceCandidate<'caloric, 'source, 'liquid, 'card> {
    field: ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card>,
    advance: BinaryNonisothermalAdvance,
}
impl<'caloric, 'source, 'liquid, 'card>
    BinaryMiscibilityMaterialNonisothermalAdvanceCandidate<'caloric, 'source, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card> {
        &self.field
    }
    pub const fn advance(&self) -> &BinaryNonisothermalAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card>,
        BinaryNonisothermalAdvance,
    ) {
        (self.field, self.advance)
    }
}

/// Conservative parcel advection and its actual endpoint material evaluations.
#[derive(Debug)]
pub struct BinaryMiscibilityMaterialNonisothermalAdvectionCandidate<'caloric, 'source, 'liquid, 'card> {
    field: ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card>,
    step: BinaryNonisothermalAdvectionStep,
}
impl<'caloric, 'source, 'liquid, 'card>
    BinaryMiscibilityMaterialNonisothermalAdvectionCandidate<'caloric, 'source, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &BinaryNonisothermalAdvectionStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card>,
        BinaryNonisothermalAdvectionStep,
    ) {
        (self.field, self.step)
    }
}

/// Complete coupled transport and its actual endpoint material evaluations.
#[derive(Debug)]
pub struct BinaryMiscibilityMaterialNonisothermalTransportCandidate<
    'caloric,
    'source,
    'liquid,
    'card,
> {
    field: ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card>,
    step: BinaryNonisothermalTransportStep,
}
impl<'caloric, 'source, 'liquid, 'card>
    BinaryMiscibilityMaterialNonisothermalTransportCandidate<'caloric, 'source, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &BinaryNonisothermalTransportStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card>,
        BinaryNonisothermalTransportStep,
    ) {
        (self.field, self.step)
    }
}

/// Complete adaptive transport interval and original final material receipts.
#[derive(Debug)]
pub struct BinaryMiscibilityMaterialNonisothermalTransportAdvanceCandidate<
    'caloric,
    'source,
    'liquid,
    'card,
> {
    field: ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card>,
    advance: BinaryNonisothermalTransportAdvance,
}
impl<'caloric, 'source, 'liquid, 'card>
    BinaryMiscibilityMaterialNonisothermalTransportAdvanceCandidate<'caloric, 'source, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card> {
        &self.field
    }
    pub const fn advance(&self) -> &BinaryNonisothermalTransportAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedBinaryMiscibilityNonisothermalField<'caloric, 'source, 'liquid, 'card>,
        BinaryNonisothermalTransportAdvance,
    ) {
        (self.field, self.advance)
    }
}
