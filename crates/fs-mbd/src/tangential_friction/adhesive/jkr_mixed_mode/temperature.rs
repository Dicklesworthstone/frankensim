//! Original mixed-mode parameter work with the actual native state held fixed.
use super::{
    PreparedRigidJkrMixedModePair, RigidBodyState, RigidJkrMixedModeContact,
    RigidJkrMixedModeError, RigidJkrMixedModeState, SourceAdhesiveShear, SourceJkrMixedMode, poll,
};
use fs_contact::normal_patch::adhesive::{
    mixed_mode::motion::temperature::JkrMixedModeTemperatureReport,
    temperature::SourceJkrThermalRefresh,
};
use fs_exec::CancelGate;

/// Replacement native owner and retained history with the original source's
/// complete signed normal/shear/surface parameter-work report.
pub struct RigidJkrMixedModeTemperatureContinuation<'a> {
    previous: RigidJkrMixedModeState,
    pair: PreparedRigidJkrMixedModePair<'a>,
    state: RigidJkrMixedModeState,
    report: JkrMixedModeTemperatureReport,
}
impl<'a> RigidJkrMixedModeTemperatureContinuation<'a> {
    pub const fn previous_state(&self) -> RigidJkrMixedModeState {
        self.previous
    }
    pub const fn pair(&self) -> &PreparedRigidJkrMixedModePair<'a> {
        &self.pair
    }
    pub const fn state(&self) -> RigidJkrMixedModeState {
        self.state
    }
    /// Original A/B body momenta, poses and orientations, without any update.
    pub const fn bodies(&self) -> &[RigidBodyState; 2] {
        &self.pair.before
    }
    pub const fn report(&self) -> &JkrMixedModeTemperatureReport {
        &self.report
    }
    /// Publish these owners together only after any surrounding parameter-work
    /// or thermal owner admits the same immutable source successor.
    pub fn into_parts(
        self,
    ) -> (
        PreparedRigidJkrMixedModePair<'a>,
        RigidJkrMixedModeState,
        JkrMixedModeTemperatureReport,
    ) {
        (self.pair, self.state, self.report)
    }
}

impl PreparedRigidJkrMixedModePair<'_> {
    /// Continue the original source at held native geometry and elastic shifts.
    ///
    /// Native geometry, actual rates, material markers and tangent frame must
    /// match the accepted token before the original source owner authenticates
    /// its held-indentation refresh and pinned junction-strength successor.
    /// Full native inertia/mobility, momenta, coordinates, time, steps and every
    /// policy are copied unchanged. Only source-owned parameter identities and
    /// the resulting full-stick force, area and stored potential may change.
    ///
    /// This endpoint parameter update supplies no impulse, heat, elapsed time,
    /// thermal expansion, path between temperatures, slip or reattachment.
    /// Source instability, insufficient static capacity, a stale native state,
    /// failed work balance or cancellation leaves every accepted owner intact.
    pub fn continue_at_temperature<'new>(
        &self,
        previous: RigidJkrMixedModeState,
        refresh: &SourceJkrThermalRefresh,
        successor: &'new SourceJkrMixedMode,
        shear: &'new SourceAdhesiveShear,
        gate: &CancelGate,
    ) -> Result<RigidJkrMixedModeTemperatureContinuation<'new>, RigidJkrMixedModeError> {
        self.admit_history(previous, gate)?;
        let source = self.motion.continue_at_temperature(
            previous.motion,
            refresh,
            successor,
            shear,
            gate,
        )?;
        let (motion, state, report) = source.into_parts();
        if !report.energy_residual_j.is_finite()
            || report.energy_residual_j.abs() > self.controls.energy_tolerance_j
        {
            return Err(RigidJkrMixedModeError::ParameterWorkBalance {
                residual_j: report.energy_residual_j,
                tolerance_j: self.controls.energy_tolerance_j,
            });
        }
        let pair = PreparedRigidJkrMixedModePair {
            contact: RigidJkrMixedModeContact {
                binding: self.contact.binding,
                frame: self.contact.frame,
                reference_approach_m: self.contact.reference_approach_m,
                source: successor,
                shear,
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
        Ok(RigidJkrMixedModeTemperatureContinuation {
            previous,
            pair,
            state: RigidJkrMixedModeState {
                motion: state,
                ..previous
            },
            report,
        })
    }
}
