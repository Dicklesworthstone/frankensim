//! Actual native rolling recoil with all adhesive loss channels deposited once.

pub use super::super::native::{RigidJkrRollingSnapHeatControls, RigidJkrRollingSnapHeatError};
use super::{
    JkrRollingDissipationHeat, JkrRollingDissipationHeatReport, JkrRollingSnapHeatStamp,
    JkrRollingSnapSegment,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::{
    MassProperties, RigidBodyState,
    tangential_friction::{
        RigidContactBinding,
        adhesive::rolling::events::{
            RigidAdhesiveRollingEventContact, RigidAdhesiveRollingEventReport,
            RigidAdhesiveRollingEventState, prepare_rigid_adhesive_rolling_events,
        },
    },
};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use fs_tribo::partial_slip::TangentFrame;

#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrRollingDissipationHeatReport {
    pub mechanics: RigidAdhesiveRollingEventReport,
    pub heat: JkrRollingDissipationHeatReport,
    /// Delta K + Delta U_contact + Delta H + named outside work - prescribed
    /// normal work. Plastic loss is already in H/outside, never counted twice.
    pub combined_energy_residual_j: f64,
    pub absolute_stage_residual_j: f64,
}

/// Bodies, native spring/branch, all loss ledgers and finite calorics are staged
/// together. Caloric refusal after the native solve cannot publish any recoil.
pub struct RigidJkrRollingDissipationHeatCandidate<'s, 'm> {
    bodies_after: [RigidBodyState; 2],
    native_state: RigidAdhesiveRollingEventState,
    history: JkrRollingDissipationHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: RigidJkrRollingDissipationHeatReport,
}
impl<'s, 'm> RigidJkrRollingDissipationHeatCandidate<'s, 'm> {
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
    pub const fn report(&self) -> &RigidJkrRollingDissipationHeatReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        [RigidBodyState; 2],
        RigidAdhesiveRollingEventState,
        JkrRollingDissipationHeat<'s>,
        DepositionSimulation<'m>,
        RigidJkrRollingDissipationHeatReport,
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
    /// Prescribe normal approach, advance unforced native rolling, and convert
    /// its actual plastic loss plus outgoing junction releases. Geometry and
    /// material source remain the original sealed owners. The mechanical step
    /// duration does not advance this same-time thermal impulse's clock.
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
    ) -> Result<RigidJkrRollingDissipationHeatCandidate<'s, 'm>, RigidJkrRollingSnapHeatError> {
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
        if expected != self.snap.stamp || self.state() != state.source_state() {
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
            || h.rolling_loss.dissipated_work_j() != native.report.irreversible_rolling_loss_j
        {
            return Err(RigidJkrRollingSnapHeatError::HistoryMismatch);
        }
        let (kinetic, storage) = native.report.mechanical.as_ref().map_or((0.0, 0.0), |r| {
            (
                r.kinetic_energy_change_j,
                r.contact.rolling().storage_change_j,
            )
        });
        // The normal work/storage identity is checked independently by the
        // native owner. Eliminate it without subtracting rounded absolute
        // stores. Plastic loss is now internal to the actual caloric/outside
        // increment, so it must NOT be added as another outgoing energy here.
        let residual = (Dd::from_f64(kinetic)
            + Dd::from_f64(storage)
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
        Ok(RigidJkrRollingDissipationHeatCandidate {
            bodies_after: native.bodies_after,
            native_state: native.state_after,
            history,
            thermal,
            report: RigidJkrRollingDissipationHeatReport {
                mechanics: native.report,
                heat,
                combined_energy_residual_j: residual,
                absolute_stage_residual_j: absolute,
            },
        })
    }

    /// All owners publish together after final cancellation admission.
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
    ) -> Result<RigidJkrRollingDissipationHeatReport, RigidJkrRollingSnapHeatError> {
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
