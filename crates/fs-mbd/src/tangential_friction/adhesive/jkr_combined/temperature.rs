//! Loaded coupled normal/shear/twist storage at original material temperatures.
//!
//! The normal refresh authenticates the material successor. The complete
//! four-coordinate source then relaxes the shared area and supplies the actual
//! signed storage work. Pure-normal work cannot replace that loaded response.

use super::{
    JkrCombinedMotionError, MassProperties, PreparedRigidJkrCombinedPair, RigidBodyState,
    RigidJkrCombinedContact, RigidJkrCombinedError, RigidJkrCombinedState, SourceJkrMixedMode,
    poll,
};
use fs_contact::normal_patch::adhesive::{
    mixed_mode::combined::motion::temperature::JkrCombinedTemperatureReport,
    temperature::SourceJkrThermalRefresh,
};
use fs_exec::CancelGate;

/// Sealed original predecessor, new source binding and retained native history.
/// Source and state are adopted together, accounting for signed parameter work.
pub struct RigidJkrCombinedTemperatureContinuation<'a> {
    previous: RigidJkrCombinedState,
    pair: PreparedRigidJkrCombinedPair<'a>,
    state: RigidJkrCombinedState,
    report: JkrCombinedTemperatureReport,
}
impl<'a> RigidJkrCombinedTemperatureContinuation<'a> {
    pub const fn previous_state(&self) -> RigidJkrCombinedState {
        self.previous
    }
    pub const fn pair(&self) -> &PreparedRigidJkrCombinedPair<'a> {
        &self.pair
    }
    pub const fn state(&self) -> RigidJkrCombinedState {
        self.state
    }
    /// Original A/B poses and momentum vectors; this update applies no impulse.
    pub const fn bodies(&self) -> &[RigidBodyState; 2] {
        &self.pair.before
    }
    pub const fn properties(&self) -> &[MassProperties; 2] {
        &self.pair.properties
    }
    pub const fn report(&self) -> &JkrCombinedTemperatureReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        PreparedRigidJkrCombinedPair<'a>,
        RigidJkrCombinedState,
        JkrCombinedTemperatureReport,
    ) {
        (self.pair, self.state, self.report)
    }
}

impl PreparedRigidJkrCombinedPair<'_> {
    /// Continue the original material source at held approach, both elastic
    /// shear coordinates and elastic twist, retaining actual native momenta.
    ///
    /// The source owner authenticates the original normal thermal witness and
    /// re-evaluates the complete I/II/III equilibrium. All four stores change
    /// together with the relaxed area; a frozen old area or normal-only work
    /// cannot substitute for this response. The source's original stable-branch,
    /// material, static-capacity, effort and numerical gates remain active.
    ///
    /// Native geometry, body association and actual rates must already match
    /// the accepted token. In particular, temperature cannot repair stale held
    /// poses after normal motion. Complete pending COM and angular travel,
    /// frame, material markers, datum, inertia, full mobility, coordinates,
    /// rates, clock, step count and every budget remain unchanged.
    ///
    /// This endpoint parameter update introduces no body impulse, pose motion,
    /// thermal expansion, elapsed time, heat, loss or thermal-path stability
    /// claim. The caller owns the destination of signed parameter work and
    /// publishes the new source binding and native history together.
    pub fn continue_at_temperature<'new>(
        &self,
        previous: RigidJkrCombinedState,
        refresh: &SourceJkrThermalRefresh,
        successor: &'new SourceJkrMixedMode,
        gate: &CancelGate,
    ) -> Result<RigidJkrCombinedTemperatureContinuation<'new>, RigidJkrCombinedError> {
        self.admit_history(previous, gate)?;
        let continuation =
            self.motion
                .continue_at_temperature(previous.motion, refresh, successor, gate)?;
        let (motion, state, report) = continuation.into_parts();
        if !report.energy_residual_j.is_finite()
            || report.energy_residual_j.abs() > self.controls.energy_tolerance_j
        {
            return Err(JkrCombinedMotionError::EnergyBalance {
                residual_j: report.energy_residual_j,
                tolerance_j: self.controls.energy_tolerance_j,
            }
            .into());
        }
        let pair = PreparedRigidJkrCombinedPair {
            contact: RigidJkrCombinedContact {
                binding: self.contact.binding,
                frame: self.contact.frame,
                reference_approach_m: self.contact.reference_approach_m,
                source: successor,
            },
            before: self.before,
            properties: self.properties,
            marker_arms: self.marker_arms,
            reaction_arms: self.reaction_arms,
            axes: self.axes,
            markers_before: self.markers_before,
            reaction_before: self.reaction_before,
            approach_m: self.approach_m,
            mobility: self.mobility,
            motion,
            controls: self.controls,
        };
        poll(gate)?;
        Ok(RigidJkrCombinedTemperatureContinuation {
            previous,
            pair,
            state: RigidJkrCombinedState {
                motion: state,
                ..previous
            },
            report,
        })
    }
}
