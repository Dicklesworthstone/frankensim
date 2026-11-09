//! Explicit slow uniform-flow histories through the original Stokes sources.
use super::*;

/// Physical source scales and numerical acceptance for affine force and flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExactCreepingFlowStepControls {
    /// Native interval and independent work-defect budget.
    pub step: ExactCreepingStepControls,
    /// Shortest physical variation scale of the external force history [s].
    /// Required only for unequal force endpoints; unchanged by subdivision.
    pub minimum_force_variation_timescale_s: f64,
    /// Shortest physical variation scale of the ambient velocity history [s].
    /// Required only for unequal flow endpoints; unchanged by subdivision.
    /// The two endpoint samples alone cannot certify the source's bandwidth
    /// or resolve a discontinuous onset of fluid acceleration.
    pub minimum_flow_variation_timescale_s: f64,
}

/// Quasi-steady source-bound motion with explicit accelerating-fluid pressure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExactCreepingFlowStep {
    /// Original Stokes load/regime receipts, native state and independent work.
    /// The instantaneous Stokes loads exclude the separately reported pressure.
    pub flow: ExactCreepingStep,
    /// Accepted ambient endpoint, to carry explicitly into the next interval.
    pub ambient_velocity_after_world: Vec3,
    /// Original displaced-fluid pressure force `m_fluid*dU/dt` [N].
    pub ambient_pressure_force_world_n: Vec3,
    /// Signed pressure work against the actual body displacement [J].
    /// Included in `flow.background_work_j`, excluded from applied work and heat.
    pub ambient_pressure_work_j: f64,
    /// Conservative speed bound relative to the actual time-varying flow [m/s].
    pub maximum_relative_speed_bound_m_s: f64,
    /// Diameter Reynolds number at that whole-path speed bound.
    pub maximum_translation_reynolds_bound: f64,
    /// Fluid diffusion / declared external-force variation time; zero if held.
    pub drive_diffusion_ratio: f64,
    /// Fluid diffusion / declared ambient-flow variation time; zero if held.
    pub flow_diffusion_ratio: f64,
}

impl<'a, F: FluidTransport + ?Sized> CreepingSphere<'a, F> {
    /// Integrate affine force and uniform flow under the original slow-fluid gates.
    ///
    /// Both ambient endpoints are explicit. The sphere's stored bath velocity
    /// is not silently reused on continuation; pass the accepted receipt's
    /// endpoint as the next starting value. The original material snapshot,
    /// spherical geometry, Navier interface and resistance remain fixed.
    ///
    /// The existing physical pressure owner supplies `m_fluid*dU/dt`; native
    /// mechanics integrates it with external force, buoyancy, gravity and
    /// slip drag. External applied work excludes pressure and buoyancy. Signed
    /// fluid work includes pressure work and the independently integrated drag
    /// source work against U(t). Heat is the original positive slip/spin integral.
    ///
    /// This is an explicit quasi-steady reduction: the original heavy-body
    /// diffusion gate AND both authored source-variation scales must pass.
    /// Added mass and Basset history are not silently approximated; use the
    /// existing unsteady owner when those terms matter. Fast source histories
    /// refuse with `UnsteadyDrive`; interval refinement does not relax the
    /// physical source scale. The complete relative-velocity path is admitted,
    /// including interior maxima. Selected lift, applied couple and nonspherical
    /// inertia still refuse. No evolving-temperature or finite-fluid model is
    /// inferred. Failure/cancellation publishes no candidate.
    #[allow(clippy::too_many_arguments)] // Both authored force/flow histories are explicit.
    pub fn step_exact_with_affine_force_and_flow(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        wrench_before: Wrench,
        force_after_world: Vec3,
        ambient_velocity_before_world: Vec3,
        ambient_velocity_after_world: Vec3,
        controls: ExactCreepingFlowStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExactCreepingFlowStep, Error> {
        poll(&mut cancelled)?;
        self.admit_body(body)?;
        self.step_exact_frame_with_affine_force_and_flow(
            body.frame(),
            state,
            wrench_before,
            force_after_world,
            ambient_velocity_before_world,
            ambient_velocity_after_world,
            controls,
            cancelled,
        )
    }

    /// Same physical history for an already admitted spherical configuration.
    /// Thermal consumers own source/geometry admission and momentum transfer.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn step_exact_frame_with_affine_force_and_flow(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        wrench_before: Wrench,
        force_after_world: Vec3,
        ambient_velocity_before_world: Vec3,
        ambient_velocity_after_world: Vec3,
        controls: ExactCreepingFlowStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExactCreepingFlowStep, Error> {
        vector(ambient_velocity_before_world)?;
        vector(ambient_velocity_after_world)?;
        if ambient_velocity_before_world == self.ambient_velocity_world
            && ambient_velocity_after_world == self.ambient_velocity_world
        {
            let held = self.step_exact_frame_with_affine_force(
                frame,
                state,
                wrench_before,
                force_after_world,
                ExactCreepingAffineStepControls {
                    step: controls.step,
                    minimum_force_variation_timescale_s: controls
                        .minimum_force_variation_timescale_s,
                },
                &mut cancelled,
            )?;
            return Ok(ExactCreepingFlowStep {
                flow: held.flow,
                ambient_velocity_after_world,
                ambient_pressure_force_world_n: Vec3::ZERO,
                ambient_pressure_work_j: 0.0,
                maximum_relative_speed_bound_m_s: held.maximum_relative_speed_bound_m_s,
                maximum_translation_reynolds_bound: held.maximum_translation_reynolds_bound,
                drive_diffusion_ratio: held.drive_diffusion_ratio,
                flow_diffusion_ratio: 0.0,
            });
        }
        let drive_diffusion_ratio = if wrench_before.force_world == force_after_world {
            0.0
        } else {
            self.exact_drive_diffusion_ratio(controls.minimum_force_variation_timescale_s)?
        };
        let flow_diffusion_ratio = if ambient_velocity_before_world == ambient_velocity_after_world
        {
            0.0
        } else {
            let scale = controls.minimum_flow_variation_timescale_s;
            if !scale.is_finite() || scale <= 0.0 {
                return Err(Error::Invalid(
                    "a positive finite physical flow-variation timescale is required",
                ));
            }
            self.exact_drive_diffusion_ratio(scale)?
        };
        vector(force_after_world)?;
        let prepared = self.prepare_exact_frame_in_flow(
            frame,
            state,
            wrench_before,
            ambient_velocity_before_world,
            controls.step,
            &mut cancelled,
        )?;
        let before = components(ambient_velocity_before_world).map(dd);
        let after = components(ambient_velocity_after_world).map(dd);
        let acceleration =
            |i: usize| finite(((after[i] - before[i]) / dd(controls.step.duration_s)).to_f64());
        let pressure = self.ambient_pressure_force(Vec3::new(
            acceleration(0)?,
            acceleration(1)?,
            acceleration(2)?,
        ))?;
        let buoyancy = prepared.load_before.buoyancy_force_world_n;
        let native_before = Wrench {
            force_world: vector(wrench_before.force_world.add(buoyancy).add(pressure))?,
            torque_body: Vec3::ZERO,
        };
        let native_after = vector(force_after_world.add(buoyancy).add(pressure))?;
        let maximum_relative_speed_bound_m_s = prepared
            .native
            .affine_force_and_flow_relative_speed_bound_m_s(
                state,
                native_before,
                native_after,
                ambient_velocity_after_world,
                controls.step.duration_s,
            )
            .map_err(native_error)?;
        let maximum_translation_reynolds_bound = self
            .transport
            .reynolds(maximum_relative_speed_bound_m_s, 2.0 * self.radius_m)?;
        if maximum_translation_reynolds_bound > self.limits.maximum_reynolds {
            return Err(Error::Reynolds {
                translation: maximum_translation_reynolds_bound,
                rotation: prepared.load_before.rotation_reynolds,
                maximum: self.limits.maximum_reynolds,
            });
        }
        poll(&mut cancelled)?;
        let dynamics = prepared
            .native
            .step_affine_force_and_flow(
                state,
                native_before,
                native_after,
                ambient_velocity_after_world,
                controls.step.duration_s,
            )
            .map_err(native_error)?;
        let displacement = components(dynamics.displacement_world).map(dd);
        let pressure_work = force_work(pressure, displacement);
        let applied_work =
            dd(dynamics.applied_work_j) - force_work(buoyancy, displacement) - pressure_work;
        let background_work = dd(dynamics.background_work_j) + pressure_work;
        let flow = self.finish_exact_frame_in_flow(
            frame,
            state,
            dynamics,
            applied_work,
            background_work,
            ambient_velocity_after_world,
            prepared,
            controls.step,
            &mut cancelled,
        )?;
        Ok(ExactCreepingFlowStep {
            flow,
            ambient_velocity_after_world,
            ambient_pressure_force_world_n: pressure,
            ambient_pressure_work_j: finite(pressure_work.to_f64())?,
            maximum_relative_speed_bound_m_s,
            maximum_translation_reynolds_bound,
            drive_diffusion_ratio,
            flow_diffusion_ratio,
        })
    }
}
