//! Exact frozen-radius drag on the existing thermal and expansion timeline.

use super::*;

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
