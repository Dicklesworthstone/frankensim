//! Exact frozen-radius drag on the existing thermal and expansion timeline.

use super::*;
use fs_math::dd::Dd;
use fs_mbd::Vec3 as NativeVec3;

mod flow;
pub use flow::HeatedCreepingAffineFlow;
pub(super) use flow::AffineFlowWindow;

/// Work acceptance for analytic drag on each accepted thermal interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatedExactCreepingControls {
    /// Positive absolute budget for each native mechanical/work defect [J].
    /// The ordinary heated-motion budget also bounds their absolute sum.
    pub maximum_step_energy_defect_j: f64,
}

impl<'a> HeatedRigidBody<'a> {
    /// Continue the same heated sphere using analytic frozen-radius Stokes flow.
    ///
    /// Each actual thermal interval integrates velocity, displacement and free
    /// spin with the original exact native owner, then uses the same source
    /// expansion and momentum-preserving inertia remap as the midpoint entry.
    /// There is no mechanical relaxation-fraction ceiling. Thermal transport,
    /// the heated-motion maximum step and all source/regime/work gates still
    /// apply. Radius-dependent coefficients, Navier slip and buoyancy update
    /// from the accepted source configuration after every interval.
    ///
    /// Exactness concerns held-radius motion only: thermal/geometry staggering
    /// is still first order and requires refinement. Drag heat belongs to the
    /// exterior fluid, never also to this finite solid's caloric storage.
    /// Expansion parameter work and mechanical numerical defects remain
    /// separate. Selected rotational lift refuses. The original thermal,
    /// mechanical and observer histories publish together only on success.
    pub fn advance_creeping_sphere_exact<F: FluidTransport + ?Sized>(
        &self,
        sphere: &CreepingSphere<'_, F>,
        controls: HeatedExactCreepingControls,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, HeatedCreepingAudit), HeatedMotionError> {
        let observer = self.creeping_observer_with_stepper(
            sphere,
            CreepingStepper::Exact(controls),
            &mut checkpoint,
        )?;
        let before = observer.snapshot(&mut checkpoint)?;
        let (thermal, observer, audit) =
            self.thermal
                .advance_observed(&observer, target_s, config, &mut checkpoint)?;
        self.accept_creeping(before, thermal, observer, audit, &mut checkpoint)
    }

    /// Compose an original heat/field observer with the same exact drag stages.
    /// The added observer owns its heat once; expansion and publication use
    /// the identical common path as [`Self::advance_creeping_sphere_exact`].
    pub fn advance_creeping_sphere_exact_with<O: DepositionObserver, F: FluidTransport + ?Sized>(
        &self,
        observer: &O,
        sphere: &CreepingSphere<'_, F>,
        controls: HeatedExactCreepingControls,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, O, HeatedCreepingAudit), ObserverCouplingError<O::Error, HeatedMotionError>>
    {
        let motion = self
            .creeping_observer_with_stepper(
                sphere,
                CreepingStepper::Exact(controls),
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


/// Authored external COM-force endpoints on one heated-motion call window.
///
/// The first endpoint belongs to the body's current accepted time and the last
/// to that call's target time. The same affine history is evaluated at actual
/// thermal substeps and retries. Force excludes gravity, buoyancy and drag;
/// applied couples remain unsupported by the exact spherical flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatedCreepingAffineForce {
    pub force_before_world_n: NativeVec3,
    pub force_after_world_n: NativeVec3,
    /// Shortest physical variation time of the source history [s], independent
    /// of numerical subdivision. Positive and finite for changing endpoints.
    /// Keep this source scale across continuation calls; endpoint samples alone
    /// cannot certify bandwidth. Equal endpoints ignore this inactive field.
    pub minimum_force_variation_timescale_s: f64,
}

impl<'a> HeatedRigidBody<'a> {
    /// Drive the same heating/expansion timeline with an authored force ramp.
    ///
    /// Every accepted thermal interval receives the force at its actual start
    /// and stop times. The original material sphere integrates those endpoints
    /// with held-radius exact drag, then the existing expansion owner preserves
    /// world momentum while updating radius, inertia and sourced Navier loads.
    /// Complete-path Reynolds admission and the physical force-history diffusion
    /// gate remain with that sphere; the latter is also checked after expansion.
    ///
    /// Signed external-force work appears separately in the common audit and
    /// combined energy balance. Viscous heat remains an exterior-fluid transfer.
    /// Thermal/geometry splitting still requires refinement. This entry does
    /// not add conduction from the surrounding fluid: that existing model's
    /// Peclet admission needs its own complete-path treatment for a ramp.
    pub fn advance_creeping_sphere_exact_affine<F: FluidTransport + ?Sized>(
        &self,
        sphere: &CreepingSphere<'_, F>,
        force: HeatedCreepingAffineForce,
        controls: HeatedExactCreepingControls,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, HeatedCreepingAudit), HeatedMotionError> {
        let force = AffineForceWindow::new(self.thermal.time_s(), target_s, force)?;
        let observer = self.creeping_observer_with_stepper(
            sphere,
            CreepingStepper::ExactAffine(controls, force),
            &mut checkpoint,
        )?;
        let before = observer.snapshot(&mut checkpoint)?;
        let (thermal, observer, audit) =
            self.thermal.advance_observed(&observer, target_s, config, &mut checkpoint)?;
        self.accept_creeping(before, thermal, observer, audit, &mut checkpoint)
    }

    /// Compose an original heat/field observer on those same driven intervals.
    /// Source heat, native trajectory and thermal configuration publish together;
    /// a failed later observer, physical-regime or work check commits none.
    #[allow(clippy::too_many_arguments)] // Original source, force history and independent controls.
    pub fn advance_creeping_sphere_exact_affine_with<
        O: DepositionObserver,
        F: FluidTransport + ?Sized,
    >(
        &self,
        observer: &O,
        sphere: &CreepingSphere<'_, F>,
        force: HeatedCreepingAffineForce,
        controls: HeatedExactCreepingControls,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, O, HeatedCreepingAudit), ObserverCouplingError<O::Error, HeatedMotionError>>
    {
        let force = AffineForceWindow::new(self.thermal.time_s(), target_s, force)
            .map_err(ObserverCouplingError::Second)?;
        let motion = self.creeping_observer_with_stepper(
            sphere,
            CreepingStepper::ExactAffine(controls, force),
            &mut checkpoint,
        ).map_err(ObserverCouplingError::Second)?;
        let before = motion.snapshot(&mut checkpoint).map_err(ObserverCouplingError::Second)?;
        let pair = CoupledObservers::new(observer.clone(), motion);
        let (thermal, pair, audit) =
            self.thermal.advance_observed(&pair, target_s, config, &mut checkpoint)?;
        let (observer, motion) = pair.into_parts();
        let (next, audit) = self.accept_creeping(before, thermal, motion, audit, &mut checkpoint)
            .map_err(ObserverCouplingError::Second)?;
        Ok((next, observer, audit))
    }
}

/// Immutable authored window; numerical subdivision never reanchors the ramp.
#[derive(Clone, Copy)]
pub(super) struct AffineForceWindow {
    start_s: f64,
    stop_s: f64,
    force: HeatedCreepingAffineForce,
}

impl AffineForceWindow {
    fn new(
        start_s: f64,
        stop_s: f64,
        force: HeatedCreepingAffineForce,
    ) -> Result<Self, HeatedMotionError> {
        admit_affine_history(
            [start_s, stop_s],
            [force.force_before_world_n, force.force_after_world_n],
            force.minimum_force_variation_timescale_s,
        )?;
        Ok(Self { start_s, stop_s, force })
    }

    pub(super) fn minimum_variation_timescale_s(self) -> f64 {
        self.force.minimum_force_variation_timescale_s
    }

    pub(super) fn drive_diffusion_ratio<F: FluidTransport + ?Sized>(
        self,
        sphere: &CreepingSphere<'_, F>,
    ) -> Result<f64, CreepingFlowError> {
        if self.force.force_before_world_n == self.force.force_after_world_n {
            Ok(0.0)
        } else {
            // Use the original physical gate even if rounding makes the two
            // samples of a very short numerical interval equal.
            sphere.exact_drive_diffusion_ratio(self.minimum_variation_timescale_s())
        }
    }

    pub(super) fn endpoints(
        self,
        time_window_s: [f64; 2],
    ) -> Result<[NativeVec3; 2], CreepingFlowError> {
        if time_window_s[1] < time_window_s[0] {
            return Err(CreepingFlowError::Invalid("reversed thermal force interval"));
        }
        Ok([self.at(time_window_s[0])?, self.at(time_window_s[1])?])
    }

    fn at(self, time_s: f64) -> Result<NativeVec3, CreepingFlowError> {
        affine_vector_at(
            [self.start_s, self.stop_s],
            [self.force.force_before_world_n, self.force.force_after_world_n],
            time_s,
        )
    }
}

/// Common source-window admission; physical meanings stay with each typed history.
fn admit_affine_history(
    [start_s, stop_s]: [f64; 2],
    [before, after]: [NativeVec3; 2],
    minimum_variation_timescale_s: f64,
) -> Result<(), HeatedMotionError> {
    if !start_s.is_finite()
        || !stop_s.is_finite()
        || stop_s < start_s
        || !(stop_s - start_s).is_finite()
        || !before.is_finite()
        || !after.is_finite()
    {
        return Err(HeatedMotionError::Invalid(
            "finite source endpoints and a forward representable window required",
        ));
    }
    if before != after
        && (stop_s == start_s
            || !minimum_variation_timescale_s.is_finite()
            || minimum_variation_timescale_s <= 0.0)
    {
        return Err(HeatedMotionError::Invalid(
            "a changing source needs a positive window and physical variation timescale",
        ));
    }
    Ok(())
}

/// Endpoint-preserving interpolation shared by authored force and fluid velocity.
fn affine_vector_at(
    [start_s, stop_s]: [f64; 2],
    [before, after]: [NativeVec3; 2],
    time_s: f64,
) -> Result<NativeVec3, CreepingFlowError> {
    if !time_s.is_finite() || time_s < start_s || time_s > stop_s {
        return Err(CreepingFlowError::Invalid(
            "thermal interval is outside the authored source window",
        ));
    }
    if time_s == start_s {
        return Ok(before);
    }
    if time_s == stop_s {
        return Ok(after);
    }
    if before == after {
        return Ok(before);
    }
    let fraction = (Dd::from_f64(time_s) - Dd::from_f64(start_s))
        / (Dd::from_f64(stop_s) - Dd::from_f64(start_s));
    let complement = Dd::from_f64(1.0) - fraction;
    let interpolate = |first: f64, last: f64| {
        (complement * Dd::from_f64(first) + fraction * Dd::from_f64(last)).to_f64()
    };
    let value = NativeVec3::new(
        interpolate(before.x, after.x),
        interpolate(before.y, after.y),
        interpolate(before.z, after.z),
    );
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CreepingFlowError::Invalid("unrepresentable authored source interpolation"))
    }
}
