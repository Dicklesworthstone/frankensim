//! Driven native rolling and finite heat transport on one accepted interval.
//!
//! Prescribed normal loading and any junction release occur at the interval
//! start, after native thermal events there. The original finite-body solve
//! supplies the accepted rolling angle and plastic loss. Its source-owned heat
//! is transported at interval-average power while the native thermal timeline
//! advances. Both momenta, contact history, heat history and the clock accept
//! together; there is no endpoint heat correction or extra physical store.
//!
//! This retains held poses/frame, prescribed normal approach and frozen source
//! properties. Average-power reconstruction and the original midpoint mechanics
//! require temporal refinement; temperature feedback and finite rotation remain
//! separate physical operations.

pub use super::super::native::RigidJkrRollingSnapHeatError;
use super::super::{JkrRollingDissipationHeat, JkrRollingSnapHeatStamp};
use super::{JkrTimedRollingHeatControls, JkrTimedRollingHeatReport};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::{
    MassProperties, RigidBodyState, Wrench,
    tangential_friction::{
        RigidContactBinding,
        adhesive::rolling::{
            RigidAdhesiveRollingControls,
            events::{
                RigidAdhesiveRollingEventContact, RigidAdhesiveRollingEventReport,
                RigidAdhesiveRollingEventState,
                prepare_rigid_adhesive_rolling_events_with_wrenches,
            },
        },
    },
};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use fs_tribo::partial_slip::TangentFrame;

/// Original mechanical limits plus bounded native thermal transport and transfer.
#[derive(Clone, Copy, Debug)]
pub struct RigidJkrRollingIntervalHeatControls {
    pub mechanical: RigidAdhesiveRollingControls,
    pub heat: JkrTimedRollingHeatControls,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrRollingIntervalHeatReport {
    pub mechanics: RigidAdhesiveRollingEventReport,
    pub heat: JkrTimedRollingHeatReport,
    /// Native kinetic/contact/thermal change plus outside loss, minus prescribed
    /// normal/area work, body work, ordinary heat and incoming material energy.
    /// Internal snap and plastic transfers cancel once. Uses original relative
    /// audits rather than subtracting rounded absolute formation enthalpies.
    pub combined_energy_residual_j: f64,
    pub absolute_stage_residual_j: f64,
}

/// Fully admitted endpoints, including the advanced original thermal clock.
pub struct RigidJkrRollingIntervalHeatCandidate<'s, 'm> {
    bodies_after: [RigidBodyState; 2],
    native_state: RigidAdhesiveRollingEventState,
    history: JkrRollingDissipationHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: RigidJkrRollingIntervalHeatReport,
}
impl<'s, 'm> RigidJkrRollingIntervalHeatCandidate<'s, 'm> {
    pub const fn bodies_after(&self) -> &[RigidBodyState; 2] {
        &self.bodies_after
    }
    pub const fn native_state(&self) -> &RigidAdhesiveRollingEventState {
        &self.native_state
    }
    pub const fn history(&self) -> &JkrRollingDissipationHeat<'s> {
        &self.history
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn report(&self) -> &RigidJkrRollingIntervalHeatReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        [RigidBodyState; 2],
        RigidAdhesiveRollingEventState,
        JkrRollingDissipationHeat<'s>,
        DepositionSimulation<'m>,
        RigidJkrRollingIntervalHeatReport,
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

impl<'s> JkrRollingDissipationHeat<'s> {
    /// Advance actual held-load mechanics and original thermal transport over
    /// the same positive interval. Use zero wrenches for unforced motion.
    ///
    /// Only the internally prepared native motion supplies rolling travel;
    /// caller-authored work and edited mechanical reports cannot create heat.
    /// Source/receiver identities, the complete retained spring and the exact
    /// accepted thermal clock must continue. Scheduled thermal events keep
    /// their original ordering and original material-domain admission.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn propose_native_interval<'m>(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        state: &RigidAdhesiveRollingEventState,
        binding: RigidContactBinding,
        frame: TangentFrame,
        thermal: &DepositionSimulation<'m>,
        expected: JkrRollingSnapHeatStamp,
        indentation_m: f64,
        wrenches: [Wrench; 2],
        dt_s: f64,
        controls: RigidJkrRollingIntervalHeatControls,
        gate: &CancelGate,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<RigidJkrRollingIntervalHeatCandidate<'s, 'm>, RigidJkrRollingSnapHeatError> {
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
        if expected != self.stamp() || self.state() != state.source_state() {
            return Err(RigidJkrRollingSnapHeatError::HistoryMismatch);
        }
        let native = prepare_rigid_adhesive_rolling_events_with_wrenches(
            bodies,
            properties,
            state,
            RigidAdhesiveRollingEventContact {
                binding,
                frame,
                source: self.source(),
            },
            indentation_m,
            wrenches,
            dt_s,
            controls.mechanical,
            gate,
        )?;
        let heat = self.propose_interval(
            thermal,
            expected,
            indentation_m,
            native.report.angle_increment_rad,
            dt_s,
            controls.heat,
            &mut check,
        )?;
        let h = heat.report();
        if heat.history().state() != native.state_after.source_state()
            || h.normal.segments[0] != native.report.normal
            || h.rolling.rolling != native.report.mechanical.as_ref().map(|r| r.contact)
            || h.rolling_loss.dissipated_work_j() != native.report.irreversible_rolling_loss_j
        {
            return Err(RigidJkrRollingSnapHeatError::HistoryMismatch);
        }
        let rolling_storage = native
            .report
            .mechanical
            .as_ref()
            .map_or(0.0, |r| r.contact.rolling().storage_change_j);
        // The original heat audit subtracts the internal source losses and
        // native transport/birth inputs. Combining that relative audit with
        // this independently checked mechanical balance cancels internal loss
        // once, without subtracting large rounded caloric summary fields.
        let mechanical_balance = Dd::from_f64(native.report.kinetic_energy_change_j)
            + Dd::from_f64(rolling_storage)
            + Dd::from_f64(native.report.irreversible_rolling_loss_j)
            - Dd::from_f64(native.report.external_work_j);
        let residual = (mechanical_balance + Dd::from_f64(h.combined_energy_residual_j)).to_f64();
        let absolute = (Dd::from_f64(native.report.absolute_stage_residual_j)
            + Dd::from_f64(h.absolute_stage_energy_residual_j))
        .to_f64();
        let tolerance = (Dd::from_f64(controls.mechanical.energy_tolerance_j)
            + Dd::from_f64(controls.heat.total_energy_tolerance_j))
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
        Ok(RigidJkrRollingIntervalHeatCandidate {
            bodies_after: native.bodies_after,
            native_state: native.state_after,
            history,
            thermal,
            report: RigidJkrRollingIntervalHeatReport {
                mechanics: native.report,
                heat,
                combined_energy_residual_j: residual,
                absolute_stage_residual_j: absolute,
            },
        })
    }

    /// Publish native momenta, both retained histories and the thermal timeline.
    /// Any late transport/caloric/budget/cancellation refusal leaves all owners
    /// at their previous accepted physical endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_native_interval<'m>(
        &mut self,
        bodies: &mut [RigidBodyState],
        properties: &[MassProperties],
        state: &mut RigidAdhesiveRollingEventState,
        binding: RigidContactBinding,
        frame: TangentFrame,
        thermal: &mut DepositionSimulation<'m>,
        expected: JkrRollingSnapHeatStamp,
        indentation_m: f64,
        wrenches: [Wrench; 2],
        dt_s: f64,
        controls: RigidJkrRollingIntervalHeatControls,
        gate: &CancelGate,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<RigidJkrRollingIntervalHeatReport, RigidJkrRollingSnapHeatError> {
        let next = self.propose_native_interval(
            bodies,
            properties,
            state,
            binding,
            frame,
            thermal,
            expected,
            indentation_m,
            wrenches,
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
