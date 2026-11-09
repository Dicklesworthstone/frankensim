//! Authored ambient velocity on the original heat, motion and expansion clock.
use super::*;

/// Uniform ambient velocity endpoints on one actual heated-motion call window.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatedCreepingAffineFlow {
    /// Ambient velocity at the body's current accepted time [m/s].
    pub ambient_velocity_before_world: NativeVec3,
    /// Ambient velocity at this call's target time [m/s].
    pub ambient_velocity_after_world: NativeVec3,
    /// Shortest physical variation scale of the original flow history [s].
    /// Positive and finite for changing endpoints. Keep it across subdivision
    /// and continuation; endpoint samples alone cannot certify bandwidth or a
    /// sharp acceleration onset. Equal endpoints ignore this inactive field.
    pub minimum_flow_variation_timescale_s: f64,
}

impl<'a> HeatedRigidBody<'a> {
    /// Continue a heated sphere under authored affine force AND ambient flow.
    ///
    /// Both histories are sampled at the actual accepted thermal epochs. The
    /// existing material/native sphere integrates slip, applied work, maintained
    /// flow work and positive viscous heat, including its displaced-fluid
    /// acceleration pressure. Expansion then uses the original current source
    /// geometry and conserved-momentum remap. The endpoint load and both physical
    /// source-timescale gates are rechecked at the updated radius and inertia.
    ///
    /// The explicit flow endpoints override the constructor's stored velocity
    /// for this window; continue from the audit's accepted ambient endpoint.
    /// Pressure work is separately visible within background work and is never
    /// counted again as applied work, heat or configuration work. The original
    /// unforced and force-only entries retain their held bath and zero pressure.
    ///
    /// Exactness concerns each frozen-radius interval. Source expansion remains
    /// a first-order split, without radial flow, added mass or Basset history.
    /// The original thermal boundary/heater model remains authoritative: this
    /// entry does not infer surrounding-fluid conduction or advective heat.
    #[allow(clippy::too_many_arguments)] // Independent force/flow sources and numerical controls.
    pub fn advance_creeping_sphere_exact_affine_flow<F: FluidTransport + ?Sized>(
        &self,
        sphere: &CreepingSphere<'_, F>,
        force: HeatedCreepingAffineForce,
        flow: HeatedCreepingAffineFlow,
        controls: HeatedExactCreepingControls,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, HeatedCreepingAudit), HeatedMotionError> {
        let force = AffineForceWindow::new(self.thermal.time_s(), target_s, force)?;
        let flow = AffineFlowWindow::new(self.thermal.time_s(), target_s, flow)?;
        let observer = self.creeping_observer_with_stepper(
            sphere,
            CreepingStepper::ExactAffineFlow(controls, force, flow),
            &mut checkpoint,
        )?;
        let before = observer.snapshot(&mut checkpoint)?;
        let (thermal, observer, audit) =
            self.thermal
                .advance_observed(&observer, target_s, config, &mut checkpoint)?;
        self.accept_creeping(before, thermal, observer, audit, &mut checkpoint)
    }

    /// Compose another heat/field observer on the same force and flow history.
    /// All heat, source, expansion and continuing motion histories commit
    /// together; a late source/regime/work/cancellation refusal commits none.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_creeping_sphere_exact_affine_flow_with<
        O: DepositionObserver,
        F: FluidTransport + ?Sized,
    >(
        &self,
        observer: &O,
        sphere: &CreepingSphere<'_, F>,
        force: HeatedCreepingAffineForce,
        flow: HeatedCreepingAffineFlow,
        controls: HeatedExactCreepingControls,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, O, HeatedCreepingAudit), ObserverCouplingError<O::Error, HeatedMotionError>>
    {
        let force = AffineForceWindow::new(self.thermal.time_s(), target_s, force)
            .map_err(ObserverCouplingError::Second)?;
        let flow = AffineFlowWindow::new(self.thermal.time_s(), target_s, flow)
            .map_err(ObserverCouplingError::Second)?;
        let motion = self
            .creeping_observer_with_stepper(
                sphere,
                CreepingStepper::ExactAffineFlow(controls, force, flow),
                &mut checkpoint,
            )
            .map_err(ObserverCouplingError::Second)?;
        let before = motion
            .snapshot(&mut checkpoint)
            .map_err(ObserverCouplingError::Second)?;
        let pair = CoupledObservers::new(observer.clone(), motion);
        let (thermal, pair, audit) =
            self.thermal
                .advance_observed(&pair, target_s, config, &mut checkpoint)?;
        let (observer, motion) = pair.into_parts();
        let (next, audit) = self
            .accept_creeping(before, thermal, motion, audit, &mut checkpoint)
            .map_err(ObserverCouplingError::Second)?;
        Ok((next, observer, audit))
    }
}

#[derive(Clone, Copy)]
pub(crate) struct AffineFlowWindow {
    window_s: [f64; 2],
    flow: HeatedCreepingAffineFlow,
}

impl AffineFlowWindow {
    fn new(
        start_s: f64,
        stop_s: f64,
        flow: HeatedCreepingAffineFlow,
    ) -> Result<Self, HeatedMotionError> {
        admit_affine_history(
            [start_s, stop_s],
            [
                flow.ambient_velocity_before_world,
                flow.ambient_velocity_after_world,
            ],
            flow.minimum_flow_variation_timescale_s,
        )?;
        Ok(Self {
            window_s: [start_s, stop_s],
            flow,
        })
    }

    pub(crate) fn minimum_variation_timescale_s(self) -> f64 {
        self.flow.minimum_flow_variation_timescale_s
    }

    pub(crate) fn diffusion_ratio<F: FluidTransport + ?Sized>(
        self,
        sphere: &CreepingSphere<'_, F>,
    ) -> Result<f64, CreepingFlowError> {
        if self.flow.ambient_velocity_before_world == self.flow.ambient_velocity_after_world {
            Ok(0.0)
        } else {
            // Preserve the original source gate even when a short numerical
            // interval's rounded samples happen to coincide.
            sphere.exact_drive_diffusion_ratio(self.minimum_variation_timescale_s())
        }
    }

    pub(crate) fn at(self, time_s: f64) -> Result<NativeVec3, CreepingFlowError> {
        affine_vector_at(
            self.window_s,
            [
                self.flow.ambient_velocity_before_world,
                self.flow.ambient_velocity_after_world,
            ],
            time_s,
        )
    }
}
