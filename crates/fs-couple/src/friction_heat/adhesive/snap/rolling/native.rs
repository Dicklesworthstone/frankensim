//! Unforced native rolling, original snap release and finite calorics together.
//!
//! One prescribed normal segment precedes the original held-pose rolling phase.
//! The heat owner independently replays that exact source path and deposits its
//! outgoing release at the CURRENT thermal time. Mechanical duration does not
//! advance the thermal clock. Temperature feedback, normal body motion, external
//! wrenches and transport between events retain their separate physical owners.

use super::{
    JkrRollingSnapHeat, JkrRollingSnapHeatError, JkrRollingSnapHeatReport, JkrRollingSnapHeatStamp,
    JkrRollingSnapSegment, JkrSnapHeatControls,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::{
    MassProperties, RigidBodyState,
    tangential_friction::{
        RigidContactBinding,
        adhesive::rolling::{
            RigidAdhesiveRollingControls, RigidAdhesiveRollingError,
            events::{
                RigidAdhesiveRollingEventContact, RigidAdhesiveRollingEventReport,
                RigidAdhesiveRollingEventState, prepare_rigid_adhesive_rolling_events,
            },
        },
    },
};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use fs_tribo::partial_slip::TangentFrame;

#[derive(Clone, Copy, Debug)]
pub struct RigidJkrRollingSnapHeatControls {
    pub mechanical: RigidAdhesiveRollingControls,
    pub heat: JkrSnapHeatControls,
}

#[derive(Debug)]
pub enum RigidJkrRollingSnapHeatError {
    Mechanics(RigidAdhesiveRollingError),
    Heat(JkrRollingSnapHeatError),
    HistoryMismatch,
    EnergyBalance {
        residual_j: f64,
        absolute_stage_residual_j: f64,
    },
}
impl From<RigidAdhesiveRollingError> for RigidJkrRollingSnapHeatError {
    fn from(value: RigidAdhesiveRollingError) -> Self {
        Self::Mechanics(value)
    }
}
impl From<JkrRollingSnapHeatError> for RigidJkrRollingSnapHeatError {
    fn from(value: JkrRollingSnapHeatError) -> Self {
        Self::Heat(value)
    }
}
impl From<DepositionError> for RigidJkrRollingSnapHeatError {
    fn from(value: DepositionError) -> Self {
        Self::Heat(value.into())
    }
}
impl core::fmt::Display for RigidJkrRollingSnapHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native rolling snap heat: {self:?}")
    }
}
impl std::error::Error for RigidJkrRollingSnapHeatError {}

#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrRollingSnapHeatReport {
    pub mechanics: RigidAdhesiveRollingEventReport,
    pub heat: JkrRollingSnapHeatReport,
    /// Delta K + Delta U_contact + Delta H + plastic rolling loss + outside
    /// snap work - prescribed normal work. Internal snap release cancels once.
    pub combined_energy_residual_j: f64,
    /// Independent native and heat-transfer defects; opposite errors cannot
    /// cancel. Only the original endpoint-observation allowance applies here.
    pub absolute_stage_residual_j: f64,
}

/// All four accepted owners remain private until the complete proposal passes.
pub struct RigidJkrRollingSnapHeatCandidate<'s, 'm> {
    bodies_after: [RigidBodyState; 2],
    native_state: RigidAdhesiveRollingEventState,
    history: JkrRollingSnapHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: RigidJkrRollingSnapHeatReport,
}
impl<'s, 'm> RigidJkrRollingSnapHeatCandidate<'s, 'm> {
    pub const fn bodies_after(&self) -> &[RigidBodyState; 2] {
        &self.bodies_after
    }
    pub const fn native_state(&self) -> &RigidAdhesiveRollingEventState {
        &self.native_state
    }
    pub const fn history(&self) -> &JkrRollingSnapHeat<'s> {
        &self.history
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn report(&self) -> &RigidJkrRollingSnapHeatReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        [RigidBodyState; 2],
        RigidAdhesiveRollingEventState,
        JkrRollingSnapHeat<'s>,
        DepositionSimulation<'m>,
        RigidJkrRollingSnapHeatReport,
    ) {
        (
            self.bodies_after,
            self.native_state,
            self.history,
            self.thermal,
            self.report,
        )
    }
}

impl<'s> JkrRollingSnapHeat<'s> {
    /// Stage actual body recoil, exact source-history replay and original
    /// finite-cell release. The binding/frame must match the native state;
    /// the immutable physical source is taken only from this heat owner.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn propose_native<'m>(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        state: &RigidAdhesiveRollingEventState,
        binding: RigidContactBinding,
        frame: TangentFrame,
        thermal: &DepositionSimulation<'m>,
        expected: JkrRollingSnapHeatStamp,
        indentation_m: f64,
        dt_s: f64,
        controls: RigidJkrRollingSnapHeatControls,
        gate: &CancelGate,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<RigidJkrRollingSnapHeatCandidate<'s, 'm>, RigidJkrRollingSnapHeatError> {
        let mut check = || {
            if gate.is_requested() {
                return Err(DepositionError::Cancelled);
            }
            checkpoint()?;
            if gate.is_requested() {
                return Err(DepositionError::Cancelled);
            }
            Ok(())
        };
        check()?;
        if expected != self.stamp || self.state() != state.source_state() {
            return Err(RigidJkrRollingSnapHeatError::HistoryMismatch);
        }
        let native = prepare_rigid_adhesive_rolling_events(
            bodies,
            properties,
            state,
            RigidAdhesiveRollingEventContact {
                binding,
                frame,
                source: self.source(),
            },
            indentation_m,
            dt_s,
            controls.mechanical,
            gate,
        )?;
        let heat = self.propose(
            thermal,
            expected,
            &[
                JkrRollingSnapSegment::NormalIndentation(indentation_m),
                JkrRollingSnapSegment::RollingAngleIncrement(native.report.angle_increment_rad),
            ],
            controls.heat,
            &mut check,
        )?;
        let h = heat.report();
        if heat.history().state() != native.state_after.source_state()
            || h.segments[0] != native.report.normal
            || h.segments[1].rolling != native.report.mechanical.as_ref().map(|r| r.contact)
        {
            return Err(RigidJkrRollingSnapHeatError::HistoryMismatch);
        }
        let (kinetic, storage, loss) =
            native
                .report
                .mechanical
                .as_ref()
                .map_or((0.0, 0.0, 0.0), |r| {
                    (
                        r.kinetic_energy_change_j,
                        r.contact.rolling().storage_change_j,
                        r.contact.rolling().irreversible_loss_j,
                    )
                });
        // Eliminate the original normal-stage work/storage identity without
        // subtracting its rounded totals. Actual native and caloric increments
        // remain independent, and weak rolling release keeps its own channel.
        let residual = (Dd::from_f64(kinetic)
            + Dd::from_f64(storage)
            + Dd::from_f64(loss)
            + Dd::from_f64(h.thermal.stored_energy_change_j)
            + Dd::from_f64(h.work.other_work_j())
            - Dd::from_f64(native.report.normal.normal_released_energy_j)
            - Dd::from_f64(native.report.normal.rolling_released_energy_j))
        .to_f64();
        let absolute = (Dd::from_f64(native.report.absolute_stage_residual_j)
            + Dd::from_f64(h.transfer_residual_j.abs()))
        .to_f64();
        let tolerance = (Dd::from_f64(controls.mechanical.energy_tolerance_j)
            + Dd::from_f64(controls.heat.transfer_energy_tolerance_j))
        .to_f64();
        let observation_limit = tolerance + native.report.storage_observation_roundoff_j;
        if [residual, absolute, tolerance, observation_limit]
            .iter()
            .any(|v| !v.is_finite())
            || residual.abs() > tolerance
            || absolute > observation_limit
        {
            return Err(RigidJkrRollingSnapHeatError::EnergyBalance {
                residual_j: residual,
                absolute_stage_residual_j: absolute,
            });
        }
        check()?;
        let (history, thermal, heat) = heat.into_parts();
        Ok(RigidJkrRollingSnapHeatCandidate {
            bodies_after: native.bodies_after,
            native_state: native.state_after,
            history,
            thermal,
            report: RigidJkrRollingSnapHeatReport {
                mechanics: native.report,
                heat,
                combined_energy_residual_j: residual,
                absolute_stage_residual_j: absolute,
            },
        })
    }

    /// Publish bodies, native lifecycle, heat lifecycle and caloric cells only
    /// after final cancellation admission. Unrelated bodies are untouched.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_native<'m>(
        &mut self,
        bodies: &mut [RigidBodyState],
        properties: &[MassProperties],
        state: &mut RigidAdhesiveRollingEventState,
        binding: RigidContactBinding,
        frame: TangentFrame,
        thermal: &mut DepositionSimulation<'m>,
        expected: JkrRollingSnapHeatStamp,
        indentation_m: f64,
        dt_s: f64,
        controls: RigidJkrRollingSnapHeatControls,
        gate: &CancelGate,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<RigidJkrRollingSnapHeatReport, RigidJkrRollingSnapHeatError> {
        let next = self.propose_native(
            bodies,
            properties,
            state,
            binding,
            frame,
            thermal,
            expected,
            indentation_m,
            dt_s,
            controls,
            gate,
            &mut checkpoint,
        )?;
        if gate.is_requested() {
            return Err(DepositionError::Cancelled.into());
        }
        checkpoint()?;
        if gate.is_requested() {
            return Err(DepositionError::Cancelled.into());
        }
        let (after, native_state, history, next_thermal, report) = next.into_parts();
        bodies[binding.body_a] = after[0];
        bodies[binding.body_b] = after[1];
        *state = native_state;
        *self = history;
        *thermal = next_thermal;
        Ok(report)
    }
}
