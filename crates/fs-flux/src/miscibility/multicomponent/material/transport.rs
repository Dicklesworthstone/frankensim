//! Original material receipts for prescribed multicomponent stirring and mixing.
//!
//! Numerical advection, molecular diffusion and adaptive trials complete before
//! the final endpoint is resolved through the same borrowed homogeneous liquid
//! source. The original species chart, claims, fixed-temperature context, cell
//! budget and volume allowance persist. A late source, volume or cancellation
//! refusal publishes neither a new field nor a partial physical report.

use super::{ResolvedMulticomponentMiscibilityField, Result, poll};
use crate::miscibility::multicomponent::{
    MulticomponentMiscibilitySolve,
    adaptive::MulticomponentMiscibilityTimeControl,
    advection::{
        MulticomponentMiscibilityAdvection, MulticomponentMiscibilityAdvectionControls,
        MulticomponentMiscibilityAdvectionStep,
    },
    transport::{
        MulticomponentMiscibilityTransportStep,
        adaptive::MulticomponentMiscibilityTransportAdvance,
    },
};

impl<'source, 'card> ResolvedMulticomponentMiscibilityField<'source, 'card> {
    /// Complete closed prescribed-flow transport, then refresh all actual
    /// endpoint compositions and original volume compatibility atomically.
    pub fn advect_candidate(
        &self,
        flow: &MulticomponentMiscibilityAdvection,
        duration_s: f64,
        controls: MulticomponentMiscibilityAdvectionControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentMiscibilityMaterialAdvectionCandidate<'source, 'card>> {
        let (numerical, step) = self
            .numerical
            .advect_candidate(flow, duration_s, controls, &mut cancelled)?
            .into_parts();
        let field = self.refresh(numerical, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(MulticomponentMiscibilityMaterialAdvectionCandidate { field, step })
    }

    /// Resolve the source after both stirring and full Maxwell–Stefan exchange
    /// complete over the same physical interval. Numerical receipts are retained.
    pub fn transport_candidate(
        &self,
        flow: &MulticomponentMiscibilityAdvection,
        duration_s: f64,
        advection_controls: MulticomponentMiscibilityAdvectionControls,
        controls: MulticomponentMiscibilitySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentMiscibilityMaterialTransportCandidate<'source, 'card>> {
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
        poll(&mut cancelled)?;
        Ok(MulticomponentMiscibilityMaterialTransportCandidate { field, step })
    }

    /// Complete the entire adaptive stirring/mixing interval before final
    /// material qualification. Rejected numerical trials do not refresh sources.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_transport_candidate(
        &self,
        flow: &MulticomponentMiscibilityAdvection,
        duration_s: f64,
        advection_controls: MulticomponentMiscibilityAdvectionControls,
        controls: MulticomponentMiscibilitySolve,
        time: MulticomponentMiscibilityTimeControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentMiscibilityMaterialTransportAdvanceCandidate<'source, 'card>> {
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
        poll(&mut cancelled)?;
        Ok(MulticomponentMiscibilityMaterialTransportAdvanceCandidate { field, advance })
    }
}

/// Completed advection with the original source-qualified endpoint.
#[derive(Debug)]
pub struct MulticomponentMiscibilityMaterialAdvectionCandidate<'source, 'card> {
    field: ResolvedMulticomponentMiscibilityField<'source, 'card>,
    step: MulticomponentMiscibilityAdvectionStep,
}
impl<'source, 'card> MulticomponentMiscibilityMaterialAdvectionCandidate<'source, 'card> {
    pub const fn field(&self) -> &ResolvedMulticomponentMiscibilityField<'source, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentMiscibilityAdvectionStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentMiscibilityField<'source, 'card>,
        MulticomponentMiscibilityAdvectionStep,
    ) {
        (self.field, self.step)
    }
}

/// Stirring and molecular exchange with original endpoint material receipts.
#[derive(Debug)]
pub struct MulticomponentMiscibilityMaterialTransportCandidate<'source, 'card> {
    field: ResolvedMulticomponentMiscibilityField<'source, 'card>,
    step: MulticomponentMiscibilityTransportStep,
}
impl<'source, 'card> MulticomponentMiscibilityMaterialTransportCandidate<'source, 'card> {
    pub const fn field(&self) -> &ResolvedMulticomponentMiscibilityField<'source, 'card> {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentMiscibilityTransportStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentMiscibilityField<'source, 'card>,
        MulticomponentMiscibilityTransportStep,
    ) {
        (self.field, self.step)
    }
}

/// Full adaptive physical interval and its original-source final field.
#[derive(Debug)]
pub struct MulticomponentMiscibilityMaterialTransportAdvanceCandidate<'source, 'card> {
    field: ResolvedMulticomponentMiscibilityField<'source, 'card>,
    advance: MulticomponentMiscibilityTransportAdvance,
}
impl<'source, 'card> MulticomponentMiscibilityMaterialTransportAdvanceCandidate<'source, 'card> {
    pub const fn field(&self) -> &ResolvedMulticomponentMiscibilityField<'source, 'card> {
        &self.field
    }
    pub const fn advance(&self) -> &MulticomponentMiscibilityTransportAdvance {
        &self.advance
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedMulticomponentMiscibilityField<'source, 'card>,
        MulticomponentMiscibilityTransportAdvance,
    ) {
        (self.field, self.advance)
    }
}
