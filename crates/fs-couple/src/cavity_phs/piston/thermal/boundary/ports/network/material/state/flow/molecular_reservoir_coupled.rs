//! Complete-assembly local accuracy for internal and held-bath molecular paths.
//!
//! Every trial runs the original simultaneous transport and full feedback split.
//! Accepted fine trajectories retain native stream receipts. Local errors cover
//! both physical storage and separate incoming/outgoing exposed streams.
use super::*;
use crate::cavity_phs::piston::thermal::boundary::ports::network::material::state::flow::diffusion::coupled::{
    CoupledDiffusionControls, CoupledDiffusionEstimate, CoupledPhysicalError,
    PhysicalEstimate, PhysicalStepError, PhysicalSubstep, TransportLedger,
};

/// Complete adaptive reservoir feedback with selected gas-surface heat laws.
#[path = "molecular_reservoir_coupled_gas_law.rs"]
pub mod gas_law;
#[path = "molecular_reservoir_streams.rs"]
mod streams;
#[path = "molecular_reservoir_window.rs"]
mod window;
pub use streams::{MolecularReservoirStreamEstimate, MolecularReservoirStreamScales};
pub use window::{MolecularReservoirWindowBalance, MolecularReservoirWindowChamberBalance};

/// Local accuracy and work for one complete physical reservoir window.
/// Absolute stream scales are separate from physical storage-coordinate scales;
/// both use `state.relative_tolerance` on changes, never absolute formation U.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledMolecularReservoirControls {
    pub state: CoupledDiffusionControls,
    pub streams: MolecularReservoirStreamScales,
}

/// Componentwise maximum over accepted coarse/fine comparisons only.
/// A nearly steady chamber still has independently controlled gross streams.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CoupledMolecularReservoirEstimate {
    pub state: CoupledDiffusionEstimate,
    pub streams: MolecularReservoirStreamEstimate,
}
impl CoupledMolecularReservoirEstimate {
    pub fn maximum(self) -> f64 {
        self.state.maximum().max(self.streams.maximum())
    }
}
impl PhysicalEstimate for CoupledMolecularReservoirEstimate {
    fn from_state(state: CoupledDiffusionEstimate) -> Self {
        Self {
            state,
            streams: MolecularReservoirStreamEstimate::default(),
        }
    }
    fn maximum(self) -> f64 {
        CoupledMolecularReservoirEstimate::maximum(self)
    }
    fn accumulate(&mut self, other: Self) {
        self.state.accumulate(other.state);
        self.streams.accumulate(other.streams);
    }
}

/// Original typed physical refusals with the complete state/stream error payload.
pub type CoupledMolecularReservoirError = CoupledPhysicalError<CoupledMolecularReservoirEstimate>;
type Error = CoupledMolecularReservoirError;

/// Accepted original molecular, bath and complete feedback receipts.
/// The balance is independently checked over the original-to-final window;
/// neither physical tolerance nor retained-ledger bounds grow with subdivision.
pub struct CoupledMolecularReservoirReport<I = ResolvedMaterialStatePoint> {
    pub start_s: f64,
    pub end_s: f64,
    pub evaluations: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub minimum_accepted_substep_s: f64,
    pub maximum_accepted_substep_s: f64,
    pub estimated_error: CoupledMolecularReservoirEstimate,
    pub balance: MolecularReservoirWindowBalance,
    /// Two actual fine physical steps per accepted comparison; no coarse receipts.
    pub steps: Vec<PhysicalSubstep<MolecularReservoirAdvance<I>>>,
}

impl<'m> MaterialExchangeState<'m> {
    /// Adapt simultaneous internal/bath molecular transport and the COMPLETE
    /// chemistry, finite-wall, circuit and moving-piston feedback trajectory.
    ///
    /// Each original fine step samples current composition, T/p and gas volume.
    /// The common complete-assembly controller compares every species inventory,
    /// gas temperature/pressure/volume, piston velocity, finite-cell energy,
    /// capacitor charge and inductor flux. Independently it compares each face's
    /// incoming/outgoing species and additional heat, and each bath's sensible
    /// enthalpy in a common native species chart. Cancelling throughflow cannot
    /// hide behind an unchanged chamber or a net heat stream.
    ///
    /// Only the two original fine steps are accepted: there is no extrapolated
    /// mass, energy or state. Whole-window conservation consumes the original
    /// fixed absolute budgets, including absolute accumulated numerical defects.
    /// Native transport keeps its formation-reference cancellation; feedback
    /// retains the original reaction/caloric owner's accuracy limits. Native
    /// transport entropy production is reported separately from feedback storage.
    ///
    /// Bath composition/T/p, authored geometry, material topology and external
    /// mechanical/electrical inputs remain held. Prescribed integrated gas heat
    /// is apportioned at constant power over the window. This controls local
    /// first-order split error, not a global trajectory error or an implicit PDE.
    /// Failed/rejected probes may update the molecular workspace's last-query
    /// source receipts; they cannot publish any physical owner or port history.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_molecular_reservoir_candidate<L: TransportLaw>(
        &self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        molecular: &mut ParcelMolecularTransport<L>,
        internal_links: &[ParcelMolecularLink],
        reservoirs: &[ReservoirMolecularLink<'_>],
        molecular_controls: ParcelMolecularSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
        time: CoupledMolecularReservoirControls,
    ) -> Result<(Self, CoupledMolecularReservoirReport), Error> {
        self.coupled_molecular_reservoir_window_candidate(
            cx,
            target_s,
            ports,
            drives,
            reactions,
            molecular_controls,
            controls,
            guess_j,
            time,
            |state, end, subdrives, guess| {
                state.molecular_reservoir_step_candidate(
                    prepared,
                    cx,
                    end,
                    ports,
                    subdrives,
                    reactions,
                    gas_patches,
                    molecular,
                    internal_links,
                    reservoirs,
                    molecular_controls,
                    controls,
                    guess,
                )
            },
            <ParcelMolecularStep as TransportLedger>::refinable,
        )
    }

    // The complete open-window ledger consumes physical quantities, never the
    // selected law's receipt. Share admission, stream accuracy and publication
    // preparation without copying or erasing that receipt or its error type.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn coupled_molecular_reservoir_window_candidate<I, S>(
        &self,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        molecular_controls: ParcelMolecularSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
        time: CoupledMolecularReservoirControls,
        evolve: impl FnMut(
            &Self,
            f64,
            &[ChamberDrive],
            &[f64],
        ) -> Result<(Self, MolecularReservoirAdvance<I>), S>,
        transport_refinable: impl Fn(&S) -> bool,
    ) -> Result<
        (Self, CoupledMolecularReservoirReport<I>),
        CoupledPhysicalError<CoupledMolecularReservoirEstimate, S>,
    >
    where
        S: PhysicalStepError + From<FlowBalanceError>,
    {
        use CoupledPhysicalError as Failure;
        cx.checkpoint().map_err(|_| Failure::Cancelled)?;
        time.streams
            .validate_for_window()
            .map_err(|error| match error {
                FlowBalanceError::Invalid(reason) => Failure::Invalid(reason),
                FlowBalanceError::Gas { chamber, error } => Failure::Gas { chamber, error },
                FlowBalanceError::Cancelled => Failure::Cancelled,
            })?;
        let (mut next, window) = self.coupled_physical_window_candidate(
            cx,
            target_s,
            drives,
            guess_j,
            time.state,
            evolve,
            |state, coarse, first, second, factor| {
                Ok(CoupledMolecularReservoirEstimate {
                    state,
                    streams: streams::compare(
                        &coarse.flow,
                        &first.flow,
                        &second.flow,
                        time.streams,
                        time.state.relative_tolerance,
                        factor,
                        cx,
                    )?,
                })
            },
            transport_refinable,
        )?;
        let balance = window::summarize(
            self,
            &next,
            cx,
            ports,
            drives,
            reactions,
            molecular_controls,
            controls,
            &window.steps,
        )?;
        window::macro_ports(self, &mut next, cx, controls, &balance)?;
        cx.checkpoint().map_err(|_| Failure::Cancelled)?;
        Ok((
            next,
            CoupledMolecularReservoirReport {
                start_s: window.start_s,
                end_s: window.end_s,
                evaluations: window.evaluations,
                accepted_intervals: window.accepted_intervals,
                rejected_intervals: window.rejected_intervals,
                minimum_accepted_substep_s: window.minimum_accepted_substep_s,
                maximum_accepted_substep_s: window.maximum_accepted_substep_s,
                estimated_error: window.estimated_error,
                balance,
                steps: window.steps,
            },
        ))
    }

    /// Publish one complete accepted window after accuracy, all independent
    /// physical accounts and the final cancellation checkpoint. A late refusal
    /// leaves the live clock, gas, motion, thermal/circuit stores and histories intact.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_coupled_molecular_reservoir<L: TransportLaw>(
        &mut self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        molecular: &mut ParcelMolecularTransport<L>,
        internal_links: &[ParcelMolecularLink],
        reservoirs: &[ReservoirMolecularLink<'_>],
        molecular_controls: ParcelMolecularSolve,
        controls: NetworkSolve,
        guess_j: &[f64],
        time: CoupledMolecularReservoirControls,
    ) -> Result<CoupledMolecularReservoirReport, Error> {
        let (next, report) = self.coupled_molecular_reservoir_candidate(
            prepared,
            cx,
            target_s,
            ports,
            drives,
            reactions,
            gas_patches,
            molecular,
            internal_links,
            reservoirs,
            molecular_controls,
            controls,
            guess_j,
            time,
        )?;
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        *self = next;
        Ok(report)
    }
}
