//! Original mixture sources refreshed after conservative field transport.
//!
//! A diffuse-interface cell retains its homogeneous branch, including spinodal
//! compositions. The field already represents phase separation, so no cell is
//! independently replaced by a bulk equilibrium split. The original liquid
//! sources stay bound across transport and externally prescribed temperature
//! changes; every accepted field retains its actual active gradient model.
//!
//! Source receipt allocations are separate from the numerical graph workspace;
//! an explicit cell-count budget bounds queries and retained material results.
//! Density, viscosity and thermal properties are exposed through the existing
//! liquid consumers. The `thermal` adapter adds a closed, uniform-temperature
//! energy balance with original pure-liquid calorics. The `nonisothermal`
//! adapter resolves each evolving local temperature/composition after heat,
//! species and prescribed-flow transport; momentum remains external.

use super::advection::{BinaryMiscibilityAdvection, BinaryMiscibilityAdvectionControls};
use super::temperature::BinaryMiscibilityTemperatureChange;
use super::transport::BinaryMiscibilityTransportStep;
use super::{
    BinaryMiscibilityError, BinaryMiscibilityNetwork, BinaryMiscibilitySolve,
    BinaryMiscibilityState, poll,
};
use fs_material::liquid::solution::BinaryLiquidError;
use fs_material::liquid::solution::gradient::{
    BinaryLiquidSquareGradientLaw, ResolvedBinaryLiquidCell,
};
use fs_qty::Temperature;
use std::sync::Arc;

/// Original material sources at each evolving local temperature/composition.
pub mod nonisothermal;

/// Original pure-liquid calorics for uniform-temperature adiabatic diffusion.
pub mod thermal;

/// Field/source incompatibility, work budget, or an atomic solver/source refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum BinaryMiscibilityMaterialError {
    /// The original gradient model and source binding do not agree.
    Invalid(&'static str),
    /// Maximum number of source-resolved cells permitted by the caller.
    CellBudget { required: usize, maximum: usize },
    /// The original numerical field could not produce an admissible candidate.
    Solver(BinaryMiscibilityError),
    /// The original source binding refused a prescribed temperature.
    SourceBinding(BinaryLiquidError),
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

// Temperature updates own only a new wrapper around the exact borrowed liquid
// law. Later transport shares that wrapper, without cloning its thermodynamics.
#[derive(Clone)]
enum SourceBinding<'source, 'liquid, 'card> {
    Borrowed(&'source BinaryLiquidSquareGradientLaw<'liquid, 'card>),
    Owned(Arc<BinaryLiquidSquareGradientLaw<'liquid, 'card>>),
}
impl<'liquid, 'card> SourceBinding<'_, 'liquid, 'card> {
    fn law(&self) -> &BinaryLiquidSquareGradientLaw<'liquid, 'card> {
        match self {
            Self::Borrowed(source) => source,
            Self::Owned(source) => source.as_ref(),
        }
    }
}

/// One accepted field and its original per-cell material evaluations.
/// The original liquid sources and numerical inventories are not replaceable.
pub struct ResolvedBinaryMiscibilityField<'source, 'liquid, 'card> {
    source: SourceBinding<'source, 'liquid, 'card>,
    state: BinaryMiscibilityState,
    cells: Arc<[ResolvedBinaryLiquidCell]>,
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
            SourceBinding::Borrowed(source),
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
    /// Active temperature-specific wrapper around the original liquid sources.
    pub fn source(&self) -> &BinaryLiquidSquareGradientLaw<'liquid, 'card> {
        self.source.law()
    }
    pub const fn maximum_resolved_cells(&self) -> usize {
        self.maximum_resolved_cells
    }

    /// Share the accepted numerical owner, including its current temperature.
    /// A prescribed flow can be explicitly rebound through this owner after a
    /// temperature change without replacing the field's geometry or inventory.
    pub fn network(&self) -> BinaryMiscibilityNetwork {
        BinaryMiscibilityNetwork {
            owner: Arc::clone(&self.state.owner),
        }
    }

    /// Impose one uniform temperature and refresh every original cell source.
    /// Composition, compensated inventories, geometry and clock remain intact.
    /// Source density must still agree with the fixed site concentration. A
    /// source-domain/volume refusal or cancellation returns no partial field.
    /// The reported external free-energy change is not a heat balance.
    pub fn temperature_candidate(
        &self,
        temperature: Temperature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityMaterialTemperatureCandidate<'source, 'liquid, 'card>> {
        let candidate =
            self.network()
                .temperature_candidate(&self.state, temperature, &mut cancelled)?;
        let (_, state, change) = candidate.into_parts();
        let source = self
            .source()
            .at_temperature(temperature)
            .map_err(BinaryMiscibilityMaterialError::SourceBinding)?;
        let field = Self::resolve_owned(
            state,
            SourceBinding::Owned(Arc::new(source)),
            self.maximum_resolved_cells,
            cancelled,
        )?;
        Ok(BinaryMiscibilityMaterialTemperatureCandidate { field, change })
    }

    /// Advance advection and molecular/gradient diffusion, then refresh sources.
    /// The same active source binding resolves only the completed numerical
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
        let network = self.network();
        let candidate = network.transport_candidate(
            &self.state,
            flow,
            duration_s,
            advection_controls,
            diffusion_controls,
            &mut cancelled,
        )?;
        let (state, step) = candidate.into_parts();
        let field = Self::resolve_owned(
            state,
            self.source.clone(),
            self.maximum_resolved_cells,
            cancelled,
        )?;
        Ok(BinaryMiscibilityMaterialCandidate { field, step })
    }

    // Callers have already admitted this exact owner and cell-count budget.
    // Consume the solved state to avoid cloning its numerical vectors again.
    fn resolve_owned(
        state: BinaryMiscibilityState,
        source: SourceBinding<'source, 'liquid, 'card>,
        maximum_resolved_cells: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        let cells = resolve_cells(&state, source.law(), cancelled)?;
        Ok(Self {
            source,
            state,
            cells,
            maximum_resolved_cells,
        })
    }
}

// Every caller has admitted its exact owner and resolved-cell budget. Share
// immutable receipts when binding another numerical consumer of the same field.
fn resolve_cells(
    state: &BinaryMiscibilityState,
    source: &BinaryLiquidSquareGradientLaw<'_, '_>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Arc<[ResolvedBinaryLiquidCell]>> {
    let mut cells = Vec::with_capacity(state.component_moles().len());
    for (cell, &amounts) in state.component_moles().iter().enumerate() {
        poll(&mut cancelled)?;
        cells.push(
            source
                .resolve_cell(amounts, &mut cancelled)
                .map_err(|error| BinaryMiscibilityMaterialError::Source { cell, error })?,
        );
    }
    let cells = cells.into();
    poll(&mut cancelled)?;
    Ok(cells)
}

/// Prescribed uniform temperature and refreshed original material evaluations.
#[derive(Debug)]
pub struct BinaryMiscibilityMaterialTemperatureCandidate<'source, 'liquid, 'card> {
    field: ResolvedBinaryMiscibilityField<'source, 'liquid, 'card>,
    change: BinaryMiscibilityTemperatureChange,
}
impl<'source, 'liquid, 'card>
    BinaryMiscibilityMaterialTemperatureCandidate<'source, 'liquid, 'card>
{
    pub const fn field(&self) -> &ResolvedBinaryMiscibilityField<'source, 'liquid, 'card> {
        &self.field
    }
    pub const fn change(&self) -> &BinaryMiscibilityTemperatureChange {
        &self.change
    }
    pub fn into_parts(
        self,
    ) -> (
        ResolvedBinaryMiscibilityField<'source, 'liquid, 'card>,
        BinaryMiscibilityTemperatureChange,
    ) {
        (self.field, self.change)
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
