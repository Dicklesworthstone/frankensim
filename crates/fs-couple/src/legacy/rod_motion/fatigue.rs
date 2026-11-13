//! Fatigue on the original finite-rotation elastic rod trajectory.
//!
//! The native motion owner supplies inertia, support reactions, energy checks
//! and the actual accepted pose. The solid observer supplies source admission,
//! fixed-layer stress recovery and bounded scalar/Findley histories. Each step
//! publishes both or neither. Applied force alone is never used as the stress:
//! a freely vibrating rod continues accumulating its resolved fatigue history.
//!
//! This path retains fixed material conditions and uses the native conservative
//! elastic step. Intrinsic viscous stresses, temperature refresh, impulses and
//! contact evolution are not inferred. Only immutable native state is exposed;
//! those other motion operations cannot silently bypass the observation.

use fs_exec::CancelGate;
use fs_material::fatigue::FatigueSummary;
use fs_material::fatigue::critical_plane::findley::FindleyCycleResult;
use fs_solid::rod::material_field::fatigue::{
    RodFatigueConfig, RodFatigueError, RodFatigueGauge, RodFatigueObservation, RodFindleyGauge,
};
use fs_solid::rod::material_field::motion::RodNodalInertia;

use super::{
    PreparedRodLoads, RodMotion, RodMotionControls, RodMotionError, RodMotionInertia,
    RodMotionStep, RodNodalWrench,
};

#[derive(Debug, Clone)]
pub enum RodMotionFatigueError {
    Motion(RodMotionError),
    Fatigue(RodFatigueError),
}
impl core::fmt::Display for RodMotionFatigueError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Motion(error) => write!(f, "rod fatigue motion: {error}"),
            Self::Fatigue(error) => write!(f, "rod fatigue motion: {error}"),
        }
    }
}
impl std::error::Error for RodMotionFatigueError {}
impl From<RodMotionError> for RodMotionFatigueError {
    fn from(error: RodMotionError) -> Self {
        Self::Motion(error)
    }
}
impl From<RodFatigueError> for RodMotionFatigueError {
    fn from(error: RodFatigueError) -> Self {
        Self::Fatigue(error)
    }
}

/// One native motion state and the complete fatigue observation since admission.
/// Centered and eccentric reference inertia use the same implementation. Clones
/// share immutable sources and preserve the entire accepted mechanical/history
/// prefix; only bounded rainflow residues are copied during ordinary advances.
#[derive(Debug, Clone)]
pub struct RodMotionFatiguePath<I: RodMotionInertia = RodNodalInertia> {
    motion: RodMotion<I>,
    observation: RodFatigueObservation,
}

impl<I: RodMotionInertia> RodMotionFatiguePath<I> {
    /// Begin at the supplied native state, including preloading and momentum.
    /// Previous damage is not reconstructed and no fictitious zero is inserted.
    pub fn new(
        motion: RodMotion<I>,
        gauges: Vec<RodFatigueGauge>,
        findley: Vec<RodFindleyGauge>,
        config: RodFatigueConfig,
        gate: &CancelGate,
    ) -> Result<Self, RodMotionFatigueError> {
        let observation = RodFatigueObservation::new_motion(
            motion.model(),
            motion.pose(),
            gauges,
            findley,
            config,
            gate,
        )?;
        Ok(Self {
            motion,
            observation,
        })
    }

    /// Include transverse-force shear in both signed S-N and Findley histories.
    /// Every gauged segment must retain an explicitly isotropic, elastically
    /// homogeneous circular section. Its original Saint-Venant stress owner
    /// supplies the local field; native inertia, authored shear correction and
    /// mechanical work remain unchanged. This is quasistatic section recovery,
    /// not an end-effect or cross-sectional-inertia stress solution.
    pub fn new_with_transverse_shear(
        motion: RodMotion<I>,
        gauges: Vec<RodFatigueGauge>,
        findley: Vec<RodFindleyGauge>,
        config: RodFatigueConfig,
        gate: &CancelGate,
    ) -> Result<Self, RodMotionFatigueError> {
        let observation = RodFatigueObservation::new_motion_with_transverse_shear(
            motion.model(),
            motion.pose(),
            gauges,
            findley,
            config,
            gate,
        )?;
        Ok(Self {
            motion,
            observation,
        })
    }

    pub fn motion(&self) -> &RodMotion<I> {
        &self.motion
    }

    /// Include transverse stress from bonded circular layers with different
    /// original E/G and one common isotropic nu. The same native trajectory,
    /// source pins, whole-layer strain bounds and fatigue acceptance remain.
    pub fn new_with_layered_transverse_shear(
        motion: RodMotion<I>,
        gauges: Vec<RodFatigueGauge>,
        findley: Vec<RodFindleyGauge>,
        config: RodFatigueConfig,
        gate: &CancelGate,
    ) -> Result<Self, RodMotionFatigueError> {
        let observation = RodFatigueObservation::new_motion_with_layered_transverse_shear(
            motion.model(),
            motion.pose(),
            gauges,
            findley,
            config,
            gate,
        )?;
        Ok(Self {
            motion,
            observation,
        })
    }
    pub fn observation(&self) -> &RodFatigueObservation {
        &self.observation
    }
    pub fn sample_count(&self) -> u64 {
        self.observation.sample_count()
    }
    pub fn summaries(&self) -> Vec<FatigueSummary> {
        self.observation.summaries()
    }
    pub fn preview_fatigue(
        &self,
        gate: &CancelGate,
    ) -> Result<Vec<FatigueSummary>, RodMotionFatigueError> {
        Ok(self.observation.preview_fatigue(gate)?)
    }

    /// One native kick/drift/kick step with world-frame force and couple.
    /// The original step receipt is returned without modifying its balances.
    /// A late fatigue/source/budget refusal retains the old pose, momenta,
    /// clock, accepted-step count and every gauge history.
    pub fn advance(
        &mut self,
        dt_s: f64,
        loads: &[RodNodalWrench],
        controls: RodMotionControls,
        gate: &CancelGate,
    ) -> Result<RodMotionStep, RodMotionFatigueError> {
        self.admit_advance()?;
        let (candidate, report) = self.motion.step_candidate(dt_s, loads, controls, gate)?;
        self.accept(candidate, report, gate)
    }

    /// Same native motion under held material-point forces and explicitly
    /// body/world free couples, including original eccentric-inertia transport.
    pub fn advance_spatial(
        &mut self,
        dt_s: f64,
        loads: &PreparedRodLoads,
        controls: RodMotionControls,
        gate: &CancelGate,
    ) -> Result<RodMotionStep, RodMotionFatigueError> {
        self.admit_advance()?;
        let (candidate, report) = self
            .motion
            .step_spatial_candidate(dt_s, loads, controls, gate)?;
        self.accept(candidate, report, gate)
    }

    pub fn close_findley_cycle(
        &mut self,
        gate: &CancelGate,
    ) -> Result<Vec<FindleyCycleResult>, RodMotionFatigueError> {
        Ok(self.observation.close_findley_cycle(gate)?)
    }
    pub fn restart_findley_cycle(
        &mut self,
        gate: &CancelGate,
    ) -> Result<(), RodMotionFatigueError> {
        Ok(self.observation.restart_findley_cycle(gate)?)
    }
    /// Finalize endpoint half cycles and any explicitly closed Findley window.
    pub fn finish(
        &mut self,
        gate: &CancelGate,
    ) -> Result<Vec<FatigueSummary>, RodMotionFatigueError> {
        Ok(self.observation.finish(gate)?)
    }

    fn admit_advance(&self) -> Result<(), RodMotionFatigueError> {
        if self.observation.is_finished() {
            return Err(RodFatigueError::Finished.into());
        }
        Ok(())
    }

    fn accept(
        &mut self,
        candidate: RodMotion<I>,
        report: RodMotionStep,
        gate: &CancelGate,
    ) -> Result<RodMotionStep, RodMotionFatigueError> {
        self.observation
            .observe_motion(candidate.model(), candidate.pose(), gate)?;
        self.motion = candidate;
        Ok(report)
    }
}
