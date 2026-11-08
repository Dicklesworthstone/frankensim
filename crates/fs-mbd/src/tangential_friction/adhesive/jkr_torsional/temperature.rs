//! Original thermal parameter work with the native body pair held unchanged.
//!
//! Temperature changes the radius-relaxed normal, torsional and surface
//! storage together. The source's pure-normal refresh authenticates the
//! successor; its pure-normal work cannot replace the actual twisted work.

use super::{
    PreparedRigidJkrTorsionalPair, RigidBodyState, RigidJkrTorsionalContact,
    RigidJkrTorsionalError, RigidJkrTorsionalState, SourceJkrMixedMode, poll,
};
use fs_contact::normal_patch::adhesive::{
    mixed_mode::torsion::motion::temperature::JkrTorsionalTemperatureReport,
    temperature::SourceJkrThermalRefresh,
};
use fs_exec::CancelGate;

/// Complete immutable source successor, native pair and retained chronology.
/// Adopting the proposal requires accounting for its signed parameter work.
pub struct RigidJkrTorsionalTemperatureContinuation<'a> {
    previous: RigidJkrTorsionalState,
    pair: PreparedRigidJkrTorsionalPair<'a>,
    state: RigidJkrTorsionalState,
    report: JkrTorsionalTemperatureReport,
}
impl<'a> RigidJkrTorsionalTemperatureContinuation<'a> {
    pub const fn previous_state(&self) -> RigidJkrTorsionalState {
        self.previous
    }
    pub const fn pair(&self) -> &PreparedRigidJkrTorsionalPair<'a> {
        &self.pair
    }
    pub const fn state(&self) -> RigidJkrTorsionalState {
        self.state
    }
    /// Original A/B poses and both momentum vectors, without an update.
    pub const fn bodies(&self) -> &[RigidBodyState; 2] {
        &self.pair.before
    }
    pub const fn report(&self) -> &JkrTorsionalTemperatureReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        PreparedRigidJkrTorsionalPair<'a>,
        RigidJkrTorsionalState,
        JkrTorsionalTemperatureReport,
    ) {
        (self.pair, self.state, self.report)
    }
}

impl PreparedRigidJkrTorsionalPair<'_> {
    /// Continue the authenticated original thermal source at held approach and
    /// retained elastic twist, with actual native momenta and rates unchanged.
    ///
    /// The source owner checks the sealed held-indentation refresh, pinned
    /// original claims, unchanged model envelope and the refreshed attached
    /// I/III equilibrium. The complete relaxed normal/torsional/surface storage
    /// change is signed parameter work. No normal-only work, frozen old area or
    /// newly initialized twist substitutes for that source calculation.
    ///
    /// Native association, geometry and actual rates must match the accepted
    /// token first. In particular, stale held poses after normal travel must
    /// pass the compatible geometry owner before this handoff. Any still-pending
    /// compatible travel and strict axial eligibility are retained unchanged.
    ///
    /// Body properties, poses, momenta, full mobility, frame, marker arms,
    /// normal datum, coordinates, rates, clock, steps and every budget remain
    /// fixed. Only the authenticated source and cached contact response change.
    /// This endpoint parameter update supplies no impulse, elapsed time, heat,
    /// thermal expansion, intervening thermal-path stability or reattachment.
    pub fn continue_at_temperature<'new>(
        &self,
        previous: RigidJkrTorsionalState,
        refresh: &SourceJkrThermalRefresh,
        successor: &'new SourceJkrMixedMode,
        gate: &CancelGate,
    ) -> Result<RigidJkrTorsionalTemperatureContinuation<'new>, RigidJkrTorsionalError> {
        self.admit_history(previous, gate)?;
        let transfer =
            self.motion
                .continue_at_temperature(previous.motion, refresh, successor, gate)?;
        let (motion, state, report) = transfer.into_parts();
        if !report.energy_residual_j.is_finite()
            || report.energy_residual_j.abs() > self.controls.energy_tolerance_j
        {
            return Err(RigidJkrTorsionalError::ParameterWorkBalance {
                residual_j: report.energy_residual_j,
                tolerance_j: self.controls.energy_tolerance_j,
            });
        }
        let pair = PreparedRigidJkrTorsionalPair {
            contact: RigidJkrTorsionalContact {
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
        Ok(RigidJkrTorsionalTemperatureContinuation {
            previous,
            pair,
            state: RigidJkrTorsionalState {
                motion: state,
                ..previous
            },
            report,
        })
    }
}
