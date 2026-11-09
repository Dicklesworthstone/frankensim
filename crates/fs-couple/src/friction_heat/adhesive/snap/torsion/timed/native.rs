//! Driven finite-body twist and finite calorics on one physical interval.
//!
//! The original native solver supplies actual twist travel, retained storage
//! and plastic loss. The source lifecycle supplies the normal and twist spring
//! release at the initial normal stage. The thermal owner transports declared
//! plastic heat while its actual clock advances. Both native momenta, contact
//! history, loss accounts, cells and time accept together.
//!
//! Poses and frame remain held; normal approach is prescribed at interval start
//! and the source is frozen. The homogeneous twist and mode-III applicability
//! screen remain unchanged. Midpoint motion and average heat power require
//! temporal refinement; this does not supply mixed-mode stability or feedback.

use super::super::{JkrTorsionHeat, JkrTorsionHeatError, JkrTorsionHeatStamp};
use super::{JkrTimedTorsionHeatControls, JkrTimedTorsionHeatReport};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::{
    MassProperties, RigidBodyState, Wrench,
    tangential_friction::{
        RigidContactBinding,
        adhesive::torsion::{
            RigidAdhesiveTorsionControls, RigidAdhesiveTorsionError,
            events::{
                RigidAdhesiveTorsionEventContact, RigidAdhesiveTorsionEventReport,
                RigidAdhesiveTorsionEventState,
                prepare_rigid_adhesive_torsion_events_with_wrenches,
            },
        },
    },
};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use fs_tribo::partial_slip::TangentFrame;

#[derive(Clone, Copy, Debug)]
pub struct RigidJkrTorsionIntervalHeatControls {
    pub mechanical: RigidAdhesiveTorsionControls,
    pub heat: JkrTimedTorsionHeatControls,
}

#[derive(Debug)]
pub enum RigidJkrTorsionIntervalHeatError {
    Mechanics(RigidAdhesiveTorsionError),
    Heat(JkrTorsionHeatError),
    HistoryMismatch,
    EnergyBalance {
        residual_j: f64,
        absolute_stage_residual_j: f64,
    },
}
impl From<RigidAdhesiveTorsionError> for RigidJkrTorsionIntervalHeatError {
    fn from(error: RigidAdhesiveTorsionError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<JkrTorsionHeatError> for RigidJkrTorsionIntervalHeatError {
    fn from(error: JkrTorsionHeatError) -> Self {
        Self::Heat(error)
    }
}
impl From<DepositionError> for RigidJkrTorsionIntervalHeatError {
    fn from(error: DepositionError) -> Self {
        Self::Heat(error.into())
    }
}
impl core::fmt::Display for RigidJkrTorsionIntervalHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native torsion interval heat: {self:?}")
    }
}
impl std::error::Error for RigidJkrTorsionIntervalHeatError {}

#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrTorsionIntervalHeatReport {
    pub mechanics: RigidAdhesiveTorsionEventReport,
    pub heat: JkrTimedTorsionHeatReport,
    /// Delta K + Delta U_contact + Delta H + outside loss, minus prescribed
    /// normal/area work, body work, ordinary heat and incoming material energy.
    /// Internal contact-to-thermal transfers cancel once.
    pub combined_energy_residual_j: f64,
    pub absolute_stage_residual_j: f64,
}

pub struct RigidJkrTorsionIntervalHeatCandidate<'s, 'm> {
    bodies_after: [RigidBodyState; 2],
    native_state: RigidAdhesiveTorsionEventState,
    history: JkrTorsionHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: RigidJkrTorsionIntervalHeatReport,
}
impl<'s, 'm> RigidJkrTorsionIntervalHeatCandidate<'s, 'm> {
    pub const fn bodies_after(&self) -> &[RigidBodyState; 2] {
        &self.bodies_after
    }
    pub const fn native_state(&self) -> &RigidAdhesiveTorsionEventState {
        &self.native_state
    }
    pub const fn history(&self) -> &JkrTorsionHeat<'s> {
        &self.history
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn report(&self) -> &RigidJkrTorsionIntervalHeatReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        [RigidBodyState; 2],
        RigidAdhesiveTorsionEventState,
        JkrTorsionHeat<'s>,
        DepositionSimulation<'m>,
        RigidJkrTorsionIntervalHeatReport,
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

impl<'s> JkrTorsionHeat<'s> {
    /// Prepare driven native twist and original finite-cell thermal transport
    /// over one positive duration. Zero wrenches select unforced motion.
    ///
    /// The native solve owns angular travel and work. Its exact original
    /// normal stage, angular constitutive step and complete retained successor
    /// must replay in the heat owner. Authored energy or edited reports cannot
    /// create heat. All trial endpoints remain private until accepted together.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn propose_native_interval<'m>(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        state: &RigidAdhesiveTorsionEventState,
        binding: RigidContactBinding,
        frame: TangentFrame,
        thermal: &DepositionSimulation<'m>,
        expected: JkrTorsionHeatStamp,
        indentation_m: f64,
        wrenches: [Wrench; 2],
        dt_s: f64,
        controls: RigidJkrTorsionIntervalHeatControls,
        gate: &CancelGate,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<RigidJkrTorsionIntervalHeatCandidate<'s, 'm>, RigidJkrTorsionIntervalHeatError>
    {
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
            return Err(RigidJkrTorsionIntervalHeatError::HistoryMismatch);
        }
        let native = prepare_rigid_adhesive_torsion_events_with_wrenches(
            bodies,
            properties,
            state,
            RigidAdhesiveTorsionEventContact {
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
            || h.torsion.torsion != native.report.mechanical.as_ref().map(|r| r.contact)
            || h.torsion_loss.dissipated_work_j() != native.report.irreversible_torsion_loss_j
        {
            return Err(RigidJkrTorsionIntervalHeatError::HistoryMismatch);
        }
        let twist_storage = native
            .report
            .mechanical
            .as_ref()
            .map_or(0.0, |r| r.contact.torsion().storage_change_j);
        // The heat residual subtracts original source loss; this independent
        // mechanical balance cancels that internal plastic transfer exactly
        // once. Original incremental audits retain tiny heat beside large
        // formation enthalpy, normal work or background thermal supply.
        let mechanical_balance = Dd::from_f64(native.report.kinetic_energy_change_j)
            + Dd::from_f64(twist_storage)
            + Dd::from_f64(native.report.irreversible_torsion_loss_j)
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
            .any(|value| !value.is_finite())
            || residual.abs() > tolerance
            || absolute > observation_limit
        {
            return Err(RigidJkrTorsionIntervalHeatError::EnergyBalance {
                residual_j: residual,
                absolute_stage_residual_j: absolute,
            });
        }
        check()?;
        let (history, thermal, heat) = heat.into_parts();
        Ok(RigidJkrTorsionIntervalHeatCandidate {
            bodies_after: native.bodies_after,
            native_state: native.state_after,
            history,
            thermal,
            report: RigidJkrTorsionIntervalHeatReport {
                mechanics: native.report,
                heat,
                combined_energy_residual_j: residual,
                absolute_stage_residual_j: absolute,
            },
        })
    }

    /// Publish both momenta, native twist, heat history, calorics and clock.
    /// Every fallible phase and final cancellation check precedes publication.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_native_interval<'m>(
        &mut self,
        bodies: &mut [RigidBodyState],
        properties: &[MassProperties],
        state: &mut RigidAdhesiveTorsionEventState,
        binding: RigidContactBinding,
        frame: TangentFrame,
        thermal: &mut DepositionSimulation<'m>,
        expected: JkrTorsionHeatStamp,
        indentation_m: f64,
        wrenches: [Wrench; 2],
        dt_s: f64,
        controls: RigidJkrTorsionIntervalHeatControls,
        gate: &CancelGate,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<RigidJkrTorsionIntervalHeatReport, RigidJkrTorsionIntervalHeatError> {
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
