//! Adaptive kinetic trapping with original-source admission on every probe.
//!
//! Coarse and fine trajectories use the same source-bound implicit owner.
//! Source refusals are hard errors, including on a rejected coarse trajectory;
//! reducing a time step must never hide unsupported material conditions.
//! Only accepted half steps advance inventories, transfer totals and time.

use super::super::slab::adaptive::{
    MembraneDiffusionAdvanceError, MembraneDiffusionAdvanceWorkspace, Reservoirs,
};
use super::super::slab::trapping::{
    TrappingMembraneDiffusionState, TrappingMembraneSolve, TrappingMembraneStep,
    TrappingMembraneWorkspace,
    adaptive::{
        self as kinetics, KineticNative, TrappingMembraneAdvanceCandidate,
        TrappingMembraneAdvanceControls,
    },
};
use super::{
    ResolvedFiniteTrappingMembraneDiffusion, ResolvedTrappingMembraneDiffusion,
    ResolvedTrappingMembraneError,
};

impl<'a> KineticNative for ResolvedTrappingMembraneDiffusion<'a> {
    type Error = ResolvedTrappingMembraneError;

    fn kinetics(&self) -> &TrappingMembraneDiffusionState {
        self.state()
    }

    fn is_cancelled(error: &Self::Error) -> bool {
        *error == ResolvedTrappingMembraneError::Cancelled
    }

    fn kinetic_step(
        &self,
        reservoirs: [f64; 2],
        duration: f64,
        controls: TrappingMembraneSolve,
        work: &mut TrappingMembraneWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, TrappingMembraneStep), Self::Error> {
        Ok(self
            .step_candidate(reservoirs, duration, controls, work, cancelled)?
            .into_parts())
    }
}

impl<'a> ResolvedTrappingMembraneDiffusion<'a> {
    pub fn adaptive_workspace(
        &self,
    ) -> Result<
        MembraneDiffusionAdvanceWorkspace<TrappingMembraneWorkspace>,
        ResolvedTrappingMembraneError,
    > {
        Ok(self.state().adaptive_workspace()?)
    }

    /// Complete a constant prescribed-pressure window with independent mobile
    /// and trapped inventory error estimates. Every original D/S/interface and
    /// adsorption source is checked on every coarse and fine native probe.
    ///
    /// The returned state retains the final accepted step's source support and
    /// all original bindings. Temperature and coefficients remain fixed; split
    /// calls at prescribed boundary changes. These local time-error estimates
    /// are neither whole-window error bounds nor spatial error certificates.
    /// Any refusal leaves this state, its inventories and its time unchanged.
    pub fn advance_candidate(
        &self,
        reservoirs_pa: [f64; 2],
        duration_s: f64,
        controls: TrappingMembraneAdvanceControls,
        work: &mut MembraneDiffusionAdvanceWorkspace<TrappingMembraneWorkspace>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        TrappingMembraneAdvanceCandidate<Self>,
        MembraneDiffusionAdvanceError<ResolvedTrappingMembraneError>,
    > {
        kinetics::advance(self, reservoirs_pa, duration_s, controls, work, cancelled)
    }
}

impl<'a> KineticNative for ResolvedFiniteTrappingMembraneDiffusion<'a> {
    type Error = ResolvedTrappingMembraneError;

    fn kinetics(&self) -> &TrappingMembraneDiffusionState {
        self.membrane()
    }

    fn reservoirs(&self) -> Option<Reservoirs> {
        Some(Reservoirs {
            capacities: self.capacities().map(|capacity| capacity.value_si()),
            pressures: self.reservoir_pressures_pa(),
        })
    }

    fn is_cancelled(error: &Self::Error) -> bool {
        *error == ResolvedTrappingMembraneError::Cancelled
    }

    fn kinetic_step(
        &self,
        _: [f64; 2],
        duration: f64,
        controls: TrappingMembraneSolve,
        work: &mut TrappingMembraneWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, TrappingMembraneStep), Self::Error> {
        Ok(self
            .step_candidate(duration, controls, work, cancelled)?
            .into_parts())
    }
}

impl<'a> ResolvedFiniteTrappingMembraneDiffusion<'a> {
    pub fn adaptive_workspace(
        &self,
    ) -> Result<
        MembraneDiffusionAdvanceWorkspace<TrappingMembraneWorkspace>,
        ResolvedTrappingMembraneError,
    > {
        Ok(self.state().adaptive_workspace()?)
    }

    /// Advance both finite reservoirs and both membrane stores with one error
    /// controller. Original material laws admit every actual coarse/fine cell
    /// and ordered-face pressure envelope. Reservoir depletion cannot silently
    /// extend a source domain or change a constitutive coefficient.
    ///
    /// The report accounts for all four physical inventories. This fixed-
    /// temperature pressure-capacity model supplies no gas-mixture composition,
    /// caloric storage or thermostat-energy balance.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        controls: TrappingMembraneAdvanceControls,
        work: &mut MembraneDiffusionAdvanceWorkspace<TrappingMembraneWorkspace>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        TrappingMembraneAdvanceCandidate<Self>,
        MembraneDiffusionAdvanceError<ResolvedTrappingMembraneError>,
    > {
        kinetics::advance(self, [0.0; 2], duration_s, controls, work, cancelled)
    }
}
