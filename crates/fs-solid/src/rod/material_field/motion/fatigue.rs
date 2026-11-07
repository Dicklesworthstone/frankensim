//! Observe fixed material points at the actual accepted motion pose.

use std::sync::Arc;

use fs_exec::CancelGate;

use super::super::fatigue::{
    RodFatigueConfig, RodFatigueError, RodFatigueGauge, RodFatigueObservation, RodFindleyGauge,
    TransverseModel,
};
use super::{PreparedRodMotion, RodMotionPose};
use crate::rod::derivatives::analytic::Kinematics;

impl RodFatigueObservation {
    /// Bind one immutable native motion source and observe its supplied pose
    /// exactly once. A preloaded or already moving specimen starts at that
    /// actual pose; the source's retained initialization pose is not sampled.
    ///
    /// Future observations must retain this source owner. Cloning the prepared
    /// model is permitted; replacing material, reference geometry, temperature
    /// or body-load sources requires a new observation history. Retaining the
    /// source Arc makes this admission constant-time without raw pointer keys.
    pub fn new_motion<I>(
        model: &PreparedRodMotion<I>,
        pose: &RodMotionPose,
        gauges: Vec<RodFatigueGauge>,
        findley: Vec<RodFindleyGauge>,
        config: RodFatigueConfig,
        gate: &CancelGate,
    ) -> Result<Self, RodFatigueError> {
        let mut bank = Self::prepare(
            model.source(),
            gauges,
            findley,
            config,
            TransverseModel::Omitted,
            gate,
        )?;
        bank.motion_source = Some(Arc::clone(&model.source));
        bank.observe_motion(model, pose, gate)?;
        Ok(bank)
    }

    /// Bind the same immutable trajectory while resolving transverse shear at
    /// each fixed gauge through the original circular Saint-Venant section.
    /// The supported source family and small-strain envelope are identical to
    /// [`Self::new_with_transverse_shear`]. Only accepted poses enter histories.
    pub fn new_motion_with_transverse_shear<I>(
        model: &PreparedRodMotion<I>,
        pose: &RodMotionPose,
        gauges: Vec<RodFatigueGauge>,
        findley: Vec<RodFindleyGauge>,
        config: RodFatigueConfig,
        gate: &CancelGate,
    ) -> Result<Self, RodFatigueError> {
        let mut bank = Self::prepare(
            model.source(),
            gauges,
            findley,
            config,
            TransverseModel::Homogeneous,
            gate,
        )?;
        bank.motion_source = Some(Arc::clone(&model.source));
        bank.observe_motion(model, pose, gate)?;
        Ok(bank)
    }

    /// Observe the same accepted motion with exact bonded common-nu circular
    /// layer transmission; see [`Self::new_with_layered_transverse_shear`].
    pub fn new_motion_with_layered_transverse_shear<I>(
        model: &PreparedRodMotion<I>,
        pose: &RodMotionPose,
        gauges: Vec<RodFatigueGauge>,
        findley: Vec<RodFindleyGauge>,
        config: RodFatigueConfig,
        gate: &CancelGate,
    ) -> Result<Self, RodFatigueError> {
        let mut bank = Self::prepare(
            model.source(),
            gauges,
            findley,
            config,
            TransverseModel::LayeredCommonPoisson,
            gate,
        )?;
        bank.motion_source = Some(Arc::clone(&model.source));
        bank.observe_motion(model, pose, gate)?;
        Ok(bank)
    }

    /// Append one admitted actual pose through the original rod strain map.
    /// Only selected segments are recovered; no section-energy sweep, tangent,
    /// material requery or duplicate stress law is needed. All gauge histories
    /// and latest samples publish together or remain at their accepted prefix.
    ///
    /// This method validates the source, pose and supports, not its chronology.
    /// A native time owner must call it once per accepted endpoint and retain
    /// mechanics only after observation succeeds. This keeps the native time
    /// integrator responsible for the joint acceptance boundary.
    pub fn observe_motion<I>(
        &mut self,
        model: &PreparedRodMotion<I>,
        pose: &RodMotionPose,
        gate: &CancelGate,
    ) -> Result<(), RodFatigueError> {
        if self
            .motion_source
            .as_ref()
            .is_none_or(|source| !Arc::ptr_eq(source, &model.source))
        {
            return Err(RodFatigueError::Invalid(
                "motion fatigue must retain its original immutable native source",
            ));
        }
        model.validate_pose(pose, gate)?;
        self.accept_motion_strains(
            model.source(),
            |segment| {
                let kinematics = Kinematics::from_nodes(
                    &pose.positions_m,
                    &pose.frames,
                    model.source.reference_segment_lengths_m()[segment],
                    segment,
                )?;
                model
                    .source
                    .reference_relative_strains(segment, kinematics.strains())
            },
            gate,
        )
    }
}
