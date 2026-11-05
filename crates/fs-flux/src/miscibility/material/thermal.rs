//! Source-bound adiabatic diffusion with the original pure-liquid calorics.
//!
//! The exact liquid-law reference must match the field's original gradient
//! source. Its qualified pure-component Cp laws, pressure, species order and
//! caloric datum survive every accepted interval. Mixture transport Cp remains
//! independently sourced and is never substituted into the energy balance.
//!
//! The numerical owner admits the complete diffusion/temperature correction
//! against energy and entropy budgets before any final material queries. Cells
//! remain homogeneous, and their actual final densities must still satisfy the
//! fixed equal-partial-molar-volume model. Every failure is atomic.

use super::{
    BinaryLiquidSquareGradientLaw, BinaryMiscibilityError, BinaryMiscibilityMaterialError,
    BinaryMiscibilityNetwork, BinaryMiscibilitySolve, BinaryMiscibilityState,
    ResolvedBinaryLiquidCell, ResolvedBinaryMiscibilityField, SourceBinding, poll, resolve_cells,
};
use crate::miscibility::thermal::{
    BinaryMiscibilityThermalError, BinaryMiscibilityThermalField, BinaryMiscibilityThermalSolve,
    BinaryMiscibilityThermalStep,
};
use fs_material::liquid::solution::caloric::BinaryLiquidCaloricLaw;
use std::sync::Arc;

/// Original source mismatch or an incomplete numerical/material candidate.
#[derive(Clone, Debug, PartialEq)]
pub enum BinaryMiscibilityMaterialThermalError {
    Invalid(&'static str),
    Thermal(BinaryMiscibilityThermalError),
    Material(BinaryMiscibilityMaterialError),
}
impl core::fmt::Display for BinaryMiscibilityMaterialThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound adiabatic miscibility: {self:?}")
    }
}
impl std::error::Error for BinaryMiscibilityMaterialThermalError {}
impl From<BinaryMiscibilityThermalError> for BinaryMiscibilityMaterialThermalError {
    fn from(error: BinaryMiscibilityThermalError) -> Self {
        Self::Thermal(error)
    }
}
impl From<BinaryMiscibilityMaterialError> for BinaryMiscibilityMaterialThermalError {
    fn from(error: BinaryMiscibilityMaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<BinaryMiscibilityError> for BinaryMiscibilityMaterialThermalError {
    fn from(error: BinaryMiscibilityError) -> Self {
        Self::Material(error.into())
    }
}
type Result<T> = core::result::Result<T, BinaryMiscibilityMaterialThermalError>;

/// One numerical field, original caloric binding and final material evaluations.
/// Binding shares immutable cell receipts; no second numerical state is stored.
pub struct ResolvedBinaryMiscibilityThermalField<'caloric, 'source, 'liquid, 'card> {
    numerical: BinaryMiscibilityThermalField,
    caloric: &'caloric BinaryLiquidCaloricLaw<'liquid, 'card>,
    source: SourceBinding<'source, 'liquid, 'card>,
    cells: Arc<[ResolvedBinaryLiquidCell]>,
    maximum_resolved_cells: usize,
}
impl core::fmt::Debug for ResolvedBinaryMiscibilityThermalField<'_, '_, '_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ResolvedBinaryMiscibilityThermalField")
            .field("numerical", &self.numerical)
            .field("cells", &self.cells)
            .field("maximum_resolved_cells", &self.maximum_resolved_cells)
            .finish_non_exhaustive()
    }
}

impl<'source, 'liquid, 'card> ResolvedBinaryMiscibilityField<'source, 'liquid, 'card> {
    /// Bind the exact original liquid's qualified pure-component calorics.
    /// A different liquid-law instance refuses even if its numerical pair and
    /// property values compare equal. Existing cell receipts are shared; this
    /// operation does not query sources or advance composition/temperature.
    pub fn bind_caloric<'caloric>(
        &self,
        caloric: &'caloric BinaryLiquidCaloricLaw<'liquid, 'card>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedBinaryMiscibilityThermalField<'caloric, 'source, 'liquid, 'card>> {
        poll(&mut cancelled)?;
        if !std::ptr::eq(self.source().liquid_law(), caloric.liquid_law()) {
            return Err(BinaryMiscibilityMaterialThermalError::Invalid(
                "adiabatic field requires its exact original liquid and caloric source binding",
            ));
        }
        let numerical =
            self.network()
                .bind_caloric_field(&self.state, caloric.model(), &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(ResolvedBinaryMiscibilityThermalField {
            numerical,
            caloric,
            source: self.source.clone(),
            cells: Arc::clone(&self.cells),
            maximum_resolved_cells: self.maximum_resolved_cells,
        })
    }
}

impl<'caloric, 'source, 'liquid, 'card>
    ResolvedBinaryMiscibilityThermalField<'caloric, 'source, 'liquid, 'card>
{
    pub fn state(&self) -> &BinaryMiscibilityState {
        self.numerical.state()
    }
    pub fn network(&self) -> &BinaryMiscibilityNetwork {
        self.numerical.network()
    }
    pub fn cells(&self) -> &[ResolvedBinaryLiquidCell] {
        &self.cells
    }
    pub fn source(&self) -> &BinaryLiquidSquareGradientLaw<'liquid, 'card> {
        self.source.law()
    }
    /// Exact retained source law, including both original pure-Cp receipts.
    pub const fn caloric_law(&self) -> &'caloric BinaryLiquidCaloricLaw<'liquid, 'card> {
        self.caloric
    }
    pub const fn maximum_resolved_cells(&self) -> usize {
        self.maximum_resolved_cells
    }

    /// Advance closed diffusion and its uniform adiabatic temperature response.
    /// Only the numerically accepted final temperature/compositions refresh
    /// material sources. Final density/domain refusal or cancellation publishes
    /// neither an advanced field nor a partial heat/entropy report.
    pub fn adiabatic_candidate(
        &self,
        duration_s: f64,
        diffusion_controls: BinaryMiscibilitySolve,
        thermal_controls: BinaryMiscibilityThermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityMaterialThermalCandidate<'caloric, 'source, 'liquid, 'card>> {
        let candidate = self.numerical.adiabatic_candidate(
            duration_s,
            diffusion_controls,
            thermal_controls,
            &mut cancelled,
        )?;
        let (numerical, step) = candidate.into_parts();
        let temperature = numerical.state().model().temperature();
        let source = if temperature == self.source().gradient_model().temperature() {
            self.source.clone()
        } else {
            SourceBinding::Owned(Arc::new(
                self.source()
                    .at_temperature(temperature)
                    .map_err(BinaryMiscibilityMaterialError::SourceBinding)?,
            ))
        };
        let cells = resolve_cells(numerical.state(), source.law(), &mut cancelled)?;
        let field = Self {
            numerical,
            caloric: self.caloric,
            source,
            cells,
            maximum_resolved_cells: self.maximum_resolved_cells,
        };
        Ok(BinaryMiscibilityMaterialThermalCandidate { field, step })
    }
}

/// Complete energy/entropy-admitted field and refreshed original material cells.
#[derive(Debug)]
pub struct BinaryMiscibilityMaterialThermalCandidate<'caloric, 'source, 'liquid, 'card> {
    field: ResolvedBinaryMiscibilityThermalField<'caloric, 'source, 'liquid, 'card>,
    step: BinaryMiscibilityThermalStep,
}
impl<'caloric, 'source, 'liquid, 'card>
    BinaryMiscibilityMaterialThermalCandidate<'caloric, 'source, 'liquid, 'card>
{
    pub const fn field(
        &self,
    ) -> &ResolvedBinaryMiscibilityThermalField<'caloric, 'source, 'liquid, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &BinaryMiscibilityThermalStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedBinaryMiscibilityThermalField<'caloric, 'source, 'liquid, 'card>,
        BinaryMiscibilityThermalStep,
    ) {
        (self.field, self.step)
    }
}
