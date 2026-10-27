//! Same-endpoint chemistry on the existing circuit/enthalpy network residual.
//!
//! Listed chambers use the original coupled chemical/pressure/heat owner.
//! Unlisted chambers keep their original midpoint gas integrator. Both send
//! individual signed contact heats into the SAME receiving-cell source map;
//! the circuit, finite thermal cells and all chambers publish together.

use super::*;
use crate::cavity_phs::piston::reactive::coupled::CoupledReactivePistonSolve;
use crate::cavity_phs::piston::reactive::thermal::coupled::{
    CoupledReactiveBoundaryCandidate, CoupledReactiveBoundaryStep, CoupledReactiveThermalError,
};
use crate::cavity_phs::piston::reactive::{GasReaction, GasReactionSolve};

/// Optional chemistry on one existing original chamber. Mechanical controls
/// remain on that chamber; two competing piston or gas-inversion policies refuse.
#[derive(Clone, Copy, Debug)]
pub struct CoupledChamberReaction<'a> {
    pub chamber: u64,
    pub reaction: &'a GasReaction,
    pub chemistry: GasReactionSolve,
    pub pressure_interval_pa: [f64; 2],
    pub pressure_tolerance_pa: f64,
}

/// The actual integration result, without converting endpoint chemistry into a
/// midpoint gas receipt or discarding its original kinetic/thermal equations.
#[derive(Clone, Debug, PartialEq)]
pub enum CoupledChamberBoundary<I> {
    Nonreactive(ThermalPortsCandidate<(Contact, I)>),
    Reactive(CoupledReactiveBoundaryCandidate<Vec<ThermalPortSample<(Contact, I)>>>),
}

/// Native detailed balances after moving the original candidate's owned state.
#[derive(Clone, Debug, PartialEq)]
pub enum CoupledChamberStep {
    Nonreactive(BoundaryHeatStep),
    Reactive(CoupledReactiveBoundaryStep),
}

impl<I> CoupledChamberBoundary<I> {
    pub const fn piston(&self) -> &GasPiston {
        match self {
            Self::Nonreactive(boundary) => boundary.piston(),
            Self::Reactive(boundary) => boundary.piston(),
        }
    }
    pub fn source(&self) -> &[ThermalPortSample<(Contact, I)>] {
        match self {
            Self::Nonreactive(boundary) => boundary.source(),
            Self::Reactive(boundary) => boundary.source(),
        }
    }
    pub const fn reactive_step(&self) -> Option<&CoupledReactiveBoundaryStep> {
        match self {
            Self::Reactive(boundary) => Some(boundary.step()),
            Self::Nonreactive(_) => None,
        }
    }
    pub const fn nonreactive_step(&self) -> Option<&BoundaryHeatStep> {
        match self {
            Self::Nonreactive(boundary) => Some(boundary.step()),
            Self::Reactive(_) => None,
        }
    }
    /// All contact heats are already in the circuit/enthalpy candidate.
    /// Move each returned piston alongside that network, without advancing it
    /// again or applying its retained reaction extents a second time.
    pub fn into_parts(
        self,
    ) -> (
        GasPiston,
        Vec<ThermalPortSample<(Contact, I)>>,
        CoupledChamberStep,
    ) {
        match self {
            Self::Nonreactive(boundary) => {
                let (piston, source, step) = boundary.into_parts();
                (piston, source, CoupledChamberStep::Nonreactive(step))
            }
            Self::Reactive(boundary) => {
                let (piston, source, step) = boundary.into_parts();
                (piston, source, CoupledChamberStep::Reactive(step))
            }
        }
    }
}

/// Original native failure, with chamber identity supplied by DrivenNetworkError.
#[derive(Debug)]
pub enum CoupledBoundaryError<F> {
    Nonreactive(ThermalPistonError<F>),
    Reactive(CoupledReactiveThermalError<F>),
}
impl<F: core::fmt::Debug> core::fmt::Display for CoupledBoundaryError<F> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reacting network boundary: {self:?}")
    }
}
impl<F: std::error::Error + 'static> std::error::Error for CoupledBoundaryError<F> {}

/// Typed mixed-chamber specialization of the original shared network owner.
pub type CoupledNetworkCandidate<'m, S, I> =
    DrivenNetworkCandidate<'m, S, CoupledChamberBoundary<I>>;
pub type CoupledNetworkError<E, F> = DrivenNetworkError<E, CoupledBoundaryError<F>>;

struct ReactiveBoundary<'a> {
    reactions: BTreeMap<u64, CoupledChamberReaction<'a>>,
}
impl<I, F> BoundaryDriver<I, F> for ReactiveBoundary<'_> {
    type Boundary = CoupledChamberBoundary<I>;
    type Error = CoupledBoundaryError<F>;
    fn advance(
        &self,
        chamber: &Chamber<'_>,
        ports: &[ThermalPort],
        maximum_ports: usize,
        resolve: &mut dyn FnMut(
            ThermalPort,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<(Contact, I)>, F>,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Self::Boundary, Self::Error> {
        match self.reactions.get(&chamber.id) {
            None => NativeBoundary
                .advance(chamber, ports, maximum_ports, resolve, cancelled)
                .map(CoupledChamberBoundary::Nonreactive)
                .map_err(CoupledBoundaryError::Nonreactive),
            Some(binding) => chamber
                .piston
                .coupled_reactive_thermal_ports_candidate(
                    binding.reaction,
                    chamber.input,
                    binding.chemistry,
                    CoupledReactivePistonSolve {
                        piston: chamber.piston_solve,
                        pressure_interval_pa: binding.pressure_interval_pa,
                        pressure_tolerance_pa: binding.pressure_tolerance_pa,
                    },
                    chamber.exchange,
                    ports,
                    maximum_ports,
                    resolve,
                    cancelled,
                )
                .map(CoupledChamberBoundary::Reactive)
                .map_err(CoupledBoundaryError::Reactive),
        }
    }
    fn source(boundary: &Self::Boundary) -> &[ThermalPortSample<(Contact, I)>] {
        boundary.source()
    }
    fn account(boundary: &Self::Boundary) -> BoundaryAccount {
        match boundary {
            CoupledChamberBoundary::Nonreactive(boundary) => {
                <NativeBoundary as BoundaryDriver<I, F>>::account(boundary)
            }
            CoupledChamberBoundary::Reactive(boundary) => {
                let step = boundary.step();
                let p = &step.piston;
                BoundaryAccount {
                    gas_energy_change_j: p.reaction.reaction.stored_energy_change_j,
                    kinetic_energy_change_j: p.kinetic_energy_change_j,
                    exterior_pressure_work_j: p.exterior_pressure_work_j,
                    external_force_work_j: p.external_force_work_j,
                    support_work_j: p.support_work_j,
                    exchange_trials: step.exchange_trials,
                    material_trials: step.pressure_trials,
                }
            }
        }
    }
    fn is_cancelled(error: &Self::Error) -> bool {
        matches!(
            error,
            CoupledBoundaryError::Nonreactive(ThermalPistonError::Cancelled)
                | CoupledBoundaryError::Reactive(CoupledReactiveThermalError::Cancelled)
        )
    }
}

/// Solve electrical contacts, finite cell enthalpy, every gas-wall contact,
/// native wall motion and declared reactions through ONE outer residual.
///
/// Chemistry sees the proposed heated/cooled/compressed endpoint on EVERY
/// thermal trial. No accepted gas is pre-reacted or replaced before solving.
/// Unlisted chambers use the unchanged ordinary path; no reaction bindings
/// reproduces the original network equations and source ordering.
///
/// Each original chamber has at most one mechanism. That mechanism may contain
/// multiple original reaction channels and must use the coupled implicit policy.
/// Pure callbacks see the actual gas state and the SAME endpoint cell/material
/// snapshot. The shared owner aggregates per-contact heat at its physical cell,
/// retains electrical loss ownership and admits combined energy independently.
///
/// Independent walls, closed homogeneous gases and the original fixed thermal
/// topology/C/L remain required. This adds no gas mixing, common crank, moving
/// contact geometry or spatial flame. Refine the complete assembly interval;
/// chemistry-only subcycling cannot repair an unresolved coupled clock.
/// No state is published on source, pressure, reaction, heat, outer iteration,
/// global work-budget or cancellation refusal.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_with_contacts_candidate<'m, S, I, E, F>(
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    electrical: CircuitEnthalpyInput<'_>,
    heat_routes: &[ContactHeatCells],
    chambers: &[Chamber<'_>],
    reactions: &[CoupledChamberReaction<'_>],
    contacts: &[Contact],
    controls: NetworkSolve,
    guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<ContactThermalCoefficients<S>, E>,
    resolve_port: impl Fn(
        Contact,
        &DepositionSimulation<'m>,
        &S,
        &ResolvedThermochemicalGasState,
    ) -> Result<ConductanceSample<I>, F>,
    cancelled: impl FnMut() -> bool,
) -> Result<
    CoupledNetworkCandidate<'m, ContactThermalSample<S>, I>,
    CircuitEnthalpyError<CoupledNetworkError<E, F>>,
> {
    let cancel = RefCell::new(cancelled);
    let check = || {
        if (*cancel.borrow_mut())() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    check()?;
    if reactions.len() > chambers.len() || chambers.len() > controls.maximum_chambers {
        return Err(CircuitEnthalpyError::Invalid(
            "bounded optional reaction bindings per original chamber required",
        ));
    }
    let mut mapped = BTreeMap::new();
    for &binding in reactions {
        check()?;
        if !chambers.iter().any(|chamber| chamber.id == binding.chamber)
            || mapped.insert(binding.chamber, binding).is_some()
        {
            return Err(CircuitEnthalpyError::Invalid(
                "one optional reaction mechanism per existing chamber required",
            ));
        }
    }
    let driver = ReactiveBoundary { reactions: mapped };
    step_with_driver(
        thermal,
        target_s,
        electrical,
        heat_routes,
        chambers,
        contacts,
        controls,
        guess_j,
        &driver,
        resolve,
        resolve_port,
        || (*cancel.borrow_mut())(),
    )
}
