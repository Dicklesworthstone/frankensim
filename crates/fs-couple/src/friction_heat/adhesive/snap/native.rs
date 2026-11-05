//! Native rigid JKR motion and event-timed caloric transfer in one candidate.
//!
//! The existing rigid owner supplies actual body inertia, held wrench work and
//! contact event times. The existing thermal owner transports the declared snap
//! shares at those times. Both retain their original laws and admission gates.

use super::{
    JkrSnapHeat, JkrSnapHeatError, JkrSnapHeatStamp,
    timed::{JkrTimedSnapHeatControls, JkrTimedSnapHeatReport},
};
use fs_contact::normal_patch::adhesive::motion::JkrNormalMotionState;
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::{
    RigidBodyState, Wrench,
    tangential_friction::adhesive::jkr_normal_motion::{
        PreparedRigidJkrPair, RigidJkrError, RigidJkrReport,
    },
};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

/// Original thermal controls plus an independent whole-system energy gate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidJkrTimedSnapHeatControls {
    pub heat: JkrTimedSnapHeatControls,
    /// Bounds both the combined residual and the sum of absolute native
    /// mechanical and thermal stage defects [J]. Errors cannot cancel to pass.
    pub total_energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum RigidJkrTimedSnapHeatError {
    Mechanics(RigidJkrError),
    Heat(JkrSnapHeatError),
    EnergyBalance {
        residual_j: f64,
        absolute_stage_residual_j: f64,
    },
}
impl From<RigidJkrError> for RigidJkrTimedSnapHeatError {
    fn from(error: RigidJkrError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<JkrSnapHeatError> for RigidJkrTimedSnapHeatError {
    fn from(error: JkrSnapHeatError) -> Self {
        Self::Heat(error)
    }
}
impl core::fmt::Display for RigidJkrTimedSnapHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native rigid JKR timed heat: {self:?}")
    }
}
impl std::error::Error for RigidJkrTimedSnapHeatError {}

/// Actual native mechanical and caloric accounts with snap transfer eliminated
/// from the total-system balance. No reversible contact work becomes heat.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidJkrTimedSnapHeatReport {
    pub mechanics: RigidJkrReport,
    pub heat: JkrTimedSnapHeatReport,
    /// Native kinetic + contact + caloric storage change + named outside work
    /// - actual body-wrench work - thermal boundary/source heat - incoming
    /// material energy [J]. Internal snap release cancels between the owners.
    pub combined_energy_residual_j: f64,
    pub absolute_stage_energy_residual_j: f64,
}

/// Every body endpoint, source history and caloric state is private until all
/// original mechanical, transport, source and combined energy gates pass.
pub struct RigidJkrTimedSnapHeatCandidate<'s, 'm> {
    bodies_after: [RigidBodyState; 2],
    history: JkrSnapHeat<'s>,
    thermal: DepositionSimulation<'m>,
    report: RigidJkrTimedSnapHeatReport,
}
impl<'s, 'm> RigidJkrTimedSnapHeatCandidate<'s, 'm> {
    /// Original A/B order. Poses remain held; the surrounding pose owner must
    /// consume the reported travel before preparing the next geometry.
    #[must_use]
    pub const fn bodies_after(&self) -> &[RigidBodyState; 2] {
        &self.bodies_after
    }
    #[must_use]
    pub const fn motion_state(&self) -> JkrNormalMotionState {
        self.report.mechanics.motion.state()
    }
    #[must_use]
    pub fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    #[must_use]
    pub const fn report(&self) -> &RigidJkrTimedSnapHeatReport {
        &self.report
    }
    /// Publish these parts together after any surrounding native owners admit
    /// their own candidates. A refusal before this handoff changes no owner.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        [RigidBodyState; 2],
        JkrSnapHeat<'s>,
        DepositionSimulation<'m>,
        RigidJkrTimedSnapHeatReport,
    ) {
        (self.bodies_after, self.history, self.thermal, self.report)
    }
}

impl<'s> JkrSnapHeat<'s> {
    /// Advance loaded native bodies and transport their original snap heat.
    ///
    /// The first proposal initializes the native clock from the actual body
    /// rate and this accepted contact branch. Later proposals automatically
    /// continue the exact accepted source/mass/policy/time/step token. A stale
    /// stamp, altered native mobility or inconsistent refreshed geometry/rate
    /// refuses; the mechanical clock is never restarted to make heat replay.
    ///
    /// Wrenches use original A/B order and remain held over this phase. The
    /// caller explicitly maps the mechanical clock to the thermal timeline and
    /// supplies the native pair corresponding to these ordered heat receivers.
    /// Original scheduled thermal events execute before a same-time snap.
    /// Normal coefficients and the declared partition stay frozen. This is the
    /// original frictionless, fixed-Jacobian normal phase; pose integration,
    /// tangential adhesion and temperature feedback retain separate owners.
    #[allow(clippy::too_many_arguments)]
    pub fn propose_timed_rigid_motion<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        expected: JkrSnapHeatStamp,
        native: &PreparedRigidJkrPair<'_>,
        wrenches: [Wrench; 2],
        dt_s: f64,
        thermal_time_at_motion_zero_s: f64,
        controls: RigidJkrTimedSnapHeatControls,
        gate: &CancelGate,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<RigidJkrTimedSnapHeatCandidate<'s, 'm>, RigidJkrTimedSnapHeatError> {
        let mut checkpoint = || {
            if gate.is_requested() {
                return Err(DepositionError::Cancelled);
            }
            checkpoint()
        };
        checkpoint().map_err(JkrSnapHeatError::from)?;
        if expected != self.stamp || native.contact().normal.identity() != self.source.identity() {
            return Err(JkrSnapHeatError::HistoryMismatch.into());
        }
        if !controls.total_energy_tolerance_j.is_finite() || controls.total_energy_tolerance_j < 0.0
        {
            return Err(JkrSnapHeatError::Invalid(
                "native timed heat needs a finite nonnegative combined energy gate",
            )
            .into());
        }
        let previous = match self.normal_motion_state {
            Some(state) => state,
            None => native.initialize_motion(self.state())?,
        };
        let mechanics = native.continue_motion_with_wrenches(previous, wrenches, dt_s, gate)?;
        let heat = self.propose_timed_normal_motion(
            thermal,
            expected,
            &mechanics.report.motion,
            thermal_time_at_motion_zero_s,
            controls.heat,
            &mut checkpoint,
        )?;
        let native_residual = mechanics.report.audit.total_energy_residual_j;
        let thermal_residual = heat.report().combined_energy_residual_j;
        let residual = (Dd::from_f64(native_residual) + Dd::from_f64(thermal_residual)).to_f64();
        let absolute = (Dd::from_f64(native_residual.abs())
            + Dd::from_f64(heat.report().absolute_stage_energy_residual_j))
        .to_f64();
        if !residual.is_finite()
            || !absolute.is_finite()
            || residual.abs() > controls.total_energy_tolerance_j
            || absolute > controls.total_energy_tolerance_j
        {
            return Err(RigidJkrTimedSnapHeatError::EnergyBalance {
                residual_j: residual,
                absolute_stage_residual_j: absolute,
            });
        }
        checkpoint().map_err(JkrSnapHeatError::from)?;
        let (history, thermal, heat) = heat.into_parts();
        Ok(RigidJkrTimedSnapHeatCandidate {
            bodies_after: mechanics.bodies_after,
            history,
            thermal,
            report: RigidJkrTimedSnapHeatReport {
                mechanics: mechanics.report,
                heat,
                combined_energy_residual_j: residual,
                absolute_stage_energy_residual_j: absolute,
            },
        })
    }
}
