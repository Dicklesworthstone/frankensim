//! State-dependent heat transport on the original gas/deposition timeline.
//!
//! The supplied observer runs on every actual thermal substep, before the
//! existing molecular exchange. Its accepted history is returned with BOTH
//! physical endpoints; any subsequent heat, phase, budget or cancellation
//! failure discards all three candidates. The observer's `Clone` must therefore
//! preserve the original thermal observer's deep-copy/side-effect contract.
//!
//! This admits the existing conduction, radiation and source owners without
//! another heat equation or an extra heat debit. It is still a first-order
//! heat/phase split: coefficient and splitting error need timestep refinement.
//! Observers see the accepted mass from preceding phase intervals. Their step
//! callbacks see the heat endpoint, before that interval's molecular exchange.
//! The adaptive entry points use the original complete-split controller and
//! retain only the accepted fine observer history, including source transitions.

use super::adaptive::{
    GasTimelineAdaptiveCandidate, GasTimelineAdaptiveError, GasTimelineAdaptiveSolve,
};
use super::{
    GasChamberSolve, GasTimelineCandidate, GasTimelineEquilibriumSurface, GasTimelineError,
    GasTimelineSurface, ThermalAccommodation, advance_deposition_observed, normalize_surfaces,
};
use crate::sublimation::gas::GasSurfaceError;
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, StepConfig,
};
use fs_thermal::observers::ObserverCouplingError;
use std::collections::BTreeMap;

/// Retain the physical source's original error instead of erasing it into text.
#[derive(Clone, Debug, PartialEq)]
pub enum GasTimelineObservedError<E> {
    Timeline(GasTimelineError),
    Observer(E),
}
impl<E> From<DepositionError> for GasTimelineObservedError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Timeline(error.into())
    }
}
impl<E> From<GasSurfaceError> for GasTimelineObservedError<E> {
    fn from(error: GasSurfaceError) -> Self {
        Self::Timeline(error.into())
    }
}
impl<E> From<GasTimelineError> for GasTimelineObservedError<E> {
    fn from(error: GasTimelineError) -> Self {
        Self::Timeline(error)
    }
}
impl<E: core::fmt::Display> core::fmt::Display for GasTimelineObservedError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Timeline(error) => write!(f, "{error}"),
            Self::Observer(error) => write!(f, "deposition heat observer: {error}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for GasTimelineObservedError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Timeline(error) => Some(error),
            Self::Observer(error) => Some(error),
        }
    }
}

/// Preserve both adaptive numerical refusals and the source owner's typed error.
#[derive(Clone, Debug, PartialEq)]
pub enum GasTimelineAdaptiveObservedError<E> {
    Adaptive(GasTimelineAdaptiveError),
    Observer(E),
}
impl<E> From<GasTimelineAdaptiveError> for GasTimelineAdaptiveObservedError<E> {
    fn from(error: GasTimelineAdaptiveError) -> Self {
        Self::Adaptive(error)
    }
}
impl<E> From<DepositionError> for GasTimelineAdaptiveObservedError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Adaptive(error.into())
    }
}
impl<E> From<GasTimelineObservedError<E>> for GasTimelineAdaptiveObservedError<E> {
    fn from(error: GasTimelineObservedError<E>) -> Self {
        match error {
            GasTimelineObservedError::Timeline(error) => Self::Adaptive(error.into()),
            GasTimelineObservedError::Observer(error) => Self::Observer(error),
        }
    }
}
impl<E: core::fmt::Display> core::fmt::Display for GasTimelineAdaptiveObservedError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Adaptive(error) => write!(f, "{error}"),
            Self::Observer(error) => write!(f, "adaptive deposition heat observer: {error}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for GasTimelineAdaptiveObservedError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Adaptive(error) => Some(error),
            Self::Observer(error) => Some(error),
        }
    }
}

pub(super) fn coupled_error<E>(
    error: ObserverCouplingError<DepositionError, E>,
) -> GasTimelineObservedError<E> {
    match error {
        ObserverCouplingError::Thermal(error) | ObserverCouplingError::First(error) => error.into(),
        ObserverCouplingError::Second(error) => GasTimelineObservedError::Observer(error),
        // The branch guard owns no transport coefficients. Keep a defensive
        // refusal if that invariant changes rather than choosing either owner.
        ObserverCouplingError::ConflictingConductance(_)
        | ObserverCouplingError::ConflictingReservoirConductance(_)
        | ObserverCouplingError::ConflictingAmbientConductance(_) => {
            DepositionError::InvalidInput("phase guard and heat observer share a conductance")
                .into()
        }
    }
}

/// Couple a heat/source observer and original pure-solid surfaces to one gas.
///
/// All observer hooks, including transitions, interval sources, conductance
/// overrides and timestep restrictions, retain their original thermal meaning.
/// External source/bath heat is already included in the thermal energy audit;
/// molecular and reflected gas heat is an internal transfer counted once.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_with_gas<'thermal, 'catalog, O: DepositionObserver>(
    thermal: &DepositionSimulation<'thermal>,
    surfaces: &[GasTimelineSurface<'_>],
    accommodation: Option<&BTreeMap<CellId, ThermalAccommodation>>,
    observer: &O,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasChamberSolve,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<(GasTimelineCandidate<'thermal, 'catalog>, O), GasTimelineObservedError<O::Error>> {
    let normalized = normalize_surfaces(thermal, surfaces, phase_config, &mut checkpoint)?;
    advance_deposition_observed(
        thermal,
        &normalized,
        accommodation,
        observer,
        gas,
        catalog,
        target_s,
        thermal_config,
        phase_config,
        None,
        true,
        checkpoint,
    )
}

/// Couple the same observer to explicitly bound solid/liquid equilibrium cells.
///
/// Original caloric identities, positive solid inventory and phase ownership
/// remain enforced after every heat substep and molecular exchange. A growing
/// planar layer can use `ConductionObserver::with_planar_phase_layers` when its
/// explicitly chosen constant-density geometry applies to that phase regime.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_with_equilibrium<'thermal, 'catalog, O: DepositionObserver>(
    thermal: &DepositionSimulation<'thermal>,
    surfaces: &[GasTimelineEquilibriumSurface<'_>],
    accommodation: Option<&BTreeMap<CellId, ThermalAccommodation>>,
    observer: &O,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasChamberSolve,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<(GasTimelineCandidate<'thermal, 'catalog>, O), GasTimelineObservedError<O::Error>> {
    let normalized = normalize_surfaces(thermal, surfaces, phase_config, &mut checkpoint)?;
    advance_deposition_observed(
        thermal,
        &normalized,
        accommodation,
        observer,
        gas,
        catalog,
        target_s,
        thermal_config,
        phase_config,
        None,
        true,
        checkpoint,
    )
}

/// Adapt heat, mass-dependent conduction and molecular exchange together.
///
/// Each trial compares one complete heat/phase step with two half steps using
/// the same committed observer. Only the conservative fine state and its source
/// history survive. `begin` runs once per requested window; rejected source
/// heat, transitions and histories are discarded. All thermal work, including
/// rejected trials, spends the original single substep ceiling.
///
/// The local error norm controls original cell mass, enthalpy and temperature,
/// gas temperature and species inventories. It does not independently estimate
/// private observer-state error or supply a global/spatial error guarantee.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_adaptive_with_gas<'thermal, 'catalog, O: DepositionObserver>(
    thermal: &DepositionSimulation<'thermal>,
    surfaces: &[GasTimelineSurface<'_>],
    accommodation: Option<&BTreeMap<CellId, ThermalAccommodation>>,
    observer: &O,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasChamberSolve,
    controls: GasTimelineAdaptiveSolve,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<
    (GasTimelineAdaptiveCandidate<'thermal, 'catalog>, O),
    GasTimelineAdaptiveObservedError<O::Error>,
> {
    let normalized = normalize_surfaces(thermal, surfaces, phase_config, &mut checkpoint)?;
    super::adaptive::advance_observed(
        thermal,
        &normalized,
        accommodation,
        observer,
        gas,
        catalog,
        target_s,
        thermal_config,
        phase_config,
        controls,
        checkpoint,
    )
}

/// The same observer-aware controller through sharp solid/liquid equilibrium.
///
/// Bulk enthalpy participates in the error estimate even at a constant melting
/// temperature. Caloric identity, positive solid seeds, original source domains
/// and atomic publication retain the fixed-step equilibrium owner's semantics.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_adaptive_with_equilibrium<'thermal, 'catalog, O: DepositionObserver>(
    thermal: &DepositionSimulation<'thermal>,
    surfaces: &[GasTimelineEquilibriumSurface<'_>],
    accommodation: Option<&BTreeMap<CellId, ThermalAccommodation>>,
    observer: &O,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasChamberSolve,
    controls: GasTimelineAdaptiveSolve,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<
    (GasTimelineAdaptiveCandidate<'thermal, 'catalog>, O),
    GasTimelineAdaptiveObservedError<O::Error>,
> {
    let normalized = normalize_surfaces(thermal, surfaces, phase_config, &mut checkpoint)?;
    super::adaptive::advance_observed(
        thermal,
        &normalized,
        accommodation,
        observer,
        gas,
        catalog,
        target_s,
        thermal_config,
        phase_config,
        controls,
        checkpoint,
    )
}

#[cfg(test)]
mod tests;
