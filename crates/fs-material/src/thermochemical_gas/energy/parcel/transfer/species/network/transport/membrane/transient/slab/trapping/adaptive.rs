//! Accuracy-controlled kinetic trapping through the original shared time driver.
//!
//! Mobile pressure/amount, trapped amount, finite reservoirs and integrated
//! face amounts enter separate local error tests. Original-to-final balances
//! include BOTH mobile and trapped changes. Only accepted native half steps
//! contribute to boundary transfers; the original coupled diffusion/site solve
//! owns every coarse/fine probe. No kinetic splitting or second controller is
//! introduced. Native convergence/source refusals remain hard errors.

use super::super::adaptive::{
    self as driver, MembraneDiffusionAdvanceControls, MembraneDiffusionAdvanceError,
    MembraneDiffusionAdvanceReport, MembraneDiffusionAdvanceWorkspace, Native, Reservoirs,
};
use super::super::{
    MembraneDiffusionError, MembraneDiffusionSolve, MembraneDiffusionState, MembraneDiffusionStep,
    Sum, finite,
};
use super::{
    FiniteTrappingMembraneDiffusionState, TrappingMembraneDiffusionState, TrappingMembraneError,
    TrappingMembraneSolve, TrappingMembraneStep, TrappingMembraneWorkspace,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TrappingMembraneAdvanceControls {
    /// Includes the original algebraic balance controls for each native solve.
    pub temporal: MembraneDiffusionAdvanceControls,
    pub maximum_iterations: usize,
    pub maximum_backtracks: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TrappingMembraneAdvanceReport {
    /// Shared time/work, integrated faces, and mobile-plus-trapped balances.
    /// `last_substep.stored_*` also includes both physical stores.
    pub transport: MembraneDiffusionAdvanceReport,
    pub mobile_before_mol: f64,
    pub mobile_after_mol: f64,
    pub mobile_change_mol: f64,
    pub trapped_before_mol: f64,
    pub trapped_after_mol: f64,
    pub trapped_change_mol: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TrappingMembraneAdvanceCandidate<S = TrappingMembraneDiffusionState> {
    state: S,
    report: TrappingMembraneAdvanceReport,
}
impl<S> TrappingMembraneAdvanceCandidate<S> {
    pub const fn state(&self) -> &S {
        &self.state
    }
    pub const fn report(&self) -> &TrappingMembraneAdvanceReport {
        &self.report
    }
    pub fn into_parts(self) -> (S, TrappingMembraneAdvanceReport) {
        (self.state, self.report)
    }
}

// The source adapter can implement this same small seam while rechecking its
// original selected laws on every native probe. Coefficients never migrate into
// the temporal controller and no unbound replacement kinetics is synthesized.
pub(crate) trait KineticNative: Clone {
    type Error: From<TrappingMembraneError>;
    fn kinetics(&self) -> &TrappingMembraneDiffusionState;
    fn reservoirs(&self) -> Option<Reservoirs> {
        None
    }
    fn is_cancelled(error: &Self::Error) -> bool;
    fn kinetic_step(
        &self,
        reservoirs: [f64; 2],
        duration: f64,
        controls: TrappingMembraneSolve,
        work: &mut TrappingMembraneWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, TrappingMembraneStep), Self::Error>;
}
impl KineticNative for TrappingMembraneDiffusionState {
    type Error = TrappingMembraneError;
    fn kinetics(&self) -> &Self {
        self
    }
    fn is_cancelled(error: &Self::Error) -> bool {
        *error == TrappingMembraneError::Cancelled
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
impl KineticNative for FiniteTrappingMembraneDiffusionState {
    type Error = TrappingMembraneError;
    fn kinetics(&self) -> &TrappingMembraneDiffusionState {
        self.membrane()
    }
    fn reservoirs(&self) -> Option<Reservoirs> {
        Some(Reservoirs {
            capacities: self.capacities().map(|c| c.value_si()),
            pressures: self.reservoir_pressures_pa(),
        })
    }
    fn is_cancelled(error: &Self::Error) -> bool {
        *error == TrappingMembraneError::Cancelled
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

#[derive(Clone)]
struct Controlled<S> {
    state: S,
    maximum_iterations: usize,
    maximum_backtracks: usize,
}
impl<S: KineticNative> Native for Controlled<S> {
    type Workspace = TrappingMembraneWorkspace;
    type Error = S::Error;
    fn membrane(&self) -> &MembraneDiffusionState {
        self.state.kinetics().mobile()
    }
    fn immobile_amounts(&self) -> Option<&[f64]> {
        Some(self.state.kinetics().trapped_moles())
    }
    fn reservoirs(&self) -> Option<Reservoirs> {
        self.state.reservoirs()
    }
    fn is_cancelled(error: &Self::Error) -> bool {
        S::is_cancelled(error)
    }
    fn accuracy_controls(
        &self,
        requested: MembraneDiffusionSolve,
        limit: f64,
    ) -> MembraneDiffusionSolve {
        if requested.absolute_amount_tolerance_mol > 0.0 {
            // Zero relative contribution and a smaller absolute allowance are
            // stricter than the caller's original native row/global controls.
            MembraneDiffusionSolve {
                absolute_amount_tolerance_mol: requested
                    .absolute_amount_tolerance_mol
                    .min(0.5 * limit),
                relative_amount_tolerance: 0.0,
            }
        } else {
            // Retain relative-only native requests exactly; the shared absolute
            // audit still refuses any unresolved algebraic error afterward.
            requested
        }
    }
    fn native_step(
        &self,
        reservoirs: [f64; 2],
        duration: f64,
        balance: MembraneDiffusionSolve,
        work: &mut Self::Workspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, MembraneDiffusionStep), Self::Error> {
        let (state, report) = self.state.kinetic_step(
            reservoirs,
            duration,
            TrappingMembraneSolve {
                balance,
                maximum_iterations: self.maximum_iterations,
                maximum_backtracks: self.maximum_backtracks,
            },
            work,
            cancelled,
        )?;
        let report = transport_report(report).map_err(S::Error::from)?;
        Ok((
            Self {
                state,
                maximum_iterations: self.maximum_iterations,
                maximum_backtracks: self.maximum_backtracks,
            },
            report,
        ))
    }
}

fn transport_report(
    step: TrappingMembraneStep,
) -> Result<MembraneDiffusionStep, TrappingMembraneError> {
    Ok(MembraneDiffusionStep {
        duration_s: step.duration_s,
        reservoir_partial_pressures_pa: step.transport.reservoir_partial_pressures_pa,
        face_molar_fluxes_mol_m2_s: step.transport.face_molar_fluxes_mol_m2_s,
        face_partial_pressures_pa: step.transport.face_partial_pressures_pa,
        stored_before_mol: finite(step.mobile_before_mol + step.trapped_before_mol)?,
        stored_after_mol: finite(step.mobile_after_mol + step.trapped_after_mol)?,
        stored_change_mol: finite(step.mobile_change_mol + step.trapped_change_mol)?,
        left_reservoir_to_membrane_mol: step.boundary_transfers_mol[0],
        membrane_to_right_reservoir_mol: step.boundary_transfers_mol[1],
        amount_residual_mol: step.membrane_amount_residual_mol,
        maximum_cell_amount_residual_mol: step.maximum_cell_amount_residual_mol,
    })
}
fn native_error<E: From<TrappingMembraneError>>(
    error: MembraneDiffusionError,
) -> MembraneDiffusionAdvanceError<E> {
    MembraneDiffusionAdvanceError::Native(E::from(error.into()))
}

pub(crate) fn advance<S: KineticNative>(
    initial: &S,
    reservoirs: [f64; 2],
    duration: f64,
    controls: TrappingMembraneAdvanceControls,
    work: &mut MembraneDiffusionAdvanceWorkspace<TrappingMembraneWorkspace>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<TrappingMembraneAdvanceCandidate<S>, MembraneDiffusionAdvanceError<S::Error>> {
    let controlled = Controlled {
        state: initial.clone(),
        maximum_iterations: controls.maximum_iterations,
        maximum_backtracks: controls.maximum_backtracks,
    };
    let (next, transport) = driver::advance(
        &controlled,
        reservoirs,
        duration,
        controls.temporal,
        work,
        &mut cancelled,
    )?
    .into_parts();
    let before = initial.kinetics();
    let after = next.state.kinetics();
    let mut mobile = Sum::default();
    let mut trapped = Sum::default();
    let mut trapped_before = Sum::default();
    let mut trapped_after = Sum::default();
    for (i, cell) in before.mobile.model.cells.iter().enumerate() {
        if i % 128 == 0 && cancelled() {
            return Err(MembraneDiffusionAdvanceError::Cancelled);
        }
        mobile
            .add(
                finite(
                    (cell.capacity_mol_m2_pa * before.mobile.model.area_m2)
                        * (after.mobile.partial_pressures_pa[i]
                            - before.mobile.partial_pressures_pa[i]),
                )
                .map_err(native_error::<S::Error>)?,
            )
            .map_err(native_error::<S::Error>)?;
        trapped
            .add(after.trapped_moles[i] - before.trapped_moles[i])
            .map_err(native_error::<S::Error>)?;
        trapped_before
            .add(before.trapped_moles[i])
            .map_err(native_error::<S::Error>)?;
        trapped_after
            .add(after.trapped_moles[i])
            .map_err(native_error::<S::Error>)?;
    }
    let report = TrappingMembraneAdvanceReport {
        transport,
        mobile_before_mol: before
            .mobile
            .stored_amount_mol()
            .map_err(native_error::<S::Error>)?,
        mobile_after_mol: after
            .mobile
            .stored_amount_mol()
            .map_err(native_error::<S::Error>)?,
        mobile_change_mol: mobile.value,
        trapped_before_mol: trapped_before.value,
        trapped_after_mol: trapped_after.value,
        trapped_change_mol: trapped.value,
    };
    if cancelled() {
        return Err(MembraneDiffusionAdvanceError::Cancelled);
    }
    Ok(TrappingMembraneAdvanceCandidate {
        state: next.state,
        report,
    })
}

impl TrappingMembraneDiffusionState {
    pub fn adaptive_workspace(
        &self,
    ) -> Result<MembraneDiffusionAdvanceWorkspace<TrappingMembraneWorkspace>, TrappingMembraneError>
    {
        Ok(MembraneDiffusionAdvanceWorkspace {
            full: self.workspace()?,
            half: self.workspace()?,
        })
    }
    /// Complete a fixed-boundary window with separate mobile/trapped accuracy
    /// and shared face-transfer accounting. Positive native absolute balance
    /// tolerances are tightened for Newton convergence; relative-only controls
    /// retain their original meaning and must still pass the absolute audit.
    /// Any refusal leaves this state, both stores and its clock unchanged.
    pub fn advance_candidate(
        &self,
        reservoirs_pa: [f64; 2],
        duration_s: f64,
        controls: TrappingMembraneAdvanceControls,
        work: &mut MembraneDiffusionAdvanceWorkspace<TrappingMembraneWorkspace>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        TrappingMembraneAdvanceCandidate,
        MembraneDiffusionAdvanceError<TrappingMembraneError>,
    > {
        advance(self, reservoirs_pa, duration_s, controls, work, cancelled)
    }
}
impl FiniteTrappingMembraneDiffusionState {
    pub fn adaptive_workspace(
        &self,
    ) -> Result<MembraneDiffusionAdvanceWorkspace<TrappingMembraneWorkspace>, TrappingMembraneError>
    {
        Ok(MembraneDiffusionAdvanceWorkspace {
            full: self.workspace()?,
            half: self.workspace()?,
        })
    }
    /// Both finite reservoir inventories, mobile Henry loading and trapped
    /// loading participate in one accepted adaptive timeline and closed balance.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        controls: TrappingMembraneAdvanceControls,
        work: &mut MembraneDiffusionAdvanceWorkspace<TrappingMembraneWorkspace>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        TrappingMembraneAdvanceCandidate<Self>,
        MembraneDiffusionAdvanceError<TrappingMembraneError>,
    > {
        advance(self, [0.0; 2], duration_s, controls, work, cancelled)
    }
}

#[cfg(test)]
mod tests;
