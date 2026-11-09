//! Exact prescribed-force Stokes flow through the original native mechanics owner.
//!
//! This explicit path retains the ordinary sphere's material/interface sources,
//! buoyancy and validity checks. It removes the midpoint relaxation restriction
//! for a frozen uniform bath and zero applied couple. It does not turn that
//! numerical improvement into a broader hydrodynamic applicability claim.

use super::{CreepingLoad, CreepingSphere, Error, SphereLift, finite, norm, poll, vector};
use crate::principal::{PrincipalBody, PrincipalFrame};
use fs_material::fluid_transport::FluidTransport;
use fs_math::dd::Dd;
use fs_mbd::linear_drag::{IsotropicLinearDrag, LinearDragError, LinearDragStep};
use fs_mbd::{RigidBodyState, RigidKineticEnergyChange, Vec3, Wrench};

/// Acceptance controls for the exact frozen-reservoir sphere flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExactCreepingStepControls {
    /// Positive finite duration [s]. No body-relaxation fraction limit applies.
    pub duration_s: f64,
    /// Positive absolute budget for the independent mechanical/work residual [J].
    pub maximum_energy_defect_j: f64,
}

/// Controls for a prescribed affine-force history in a quasi-steady fluid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExactCreepingAffineStepControls {
    /// Native duration and independent mechanical/work acceptance budget.
    pub step: ExactCreepingStepControls,
    /// Shortest physical variation timescale of the authored force history [s].
    ///
    /// Positive and finite for unequal force endpoints. This belongs to the
    /// source history and must remain unchanged when that same history is
    /// subdivided; it is not the numerical substep duration. Two endpoint
    /// samples cannot establish this scale or certify a source's bandwidth.
    /// Equal endpoints use the original held entry and ignore this field.
    pub minimum_force_variation_timescale_s: f64,
}

/// Source-bound exact Stokes motion and separately retained work/heat accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExactCreepingStep {
    /// The original native body's exact drag flow and numerical residuals.
    /// Its applied force includes buoyancy and the caller's external force.
    pub dynamics: LinearDragStep,
    /// Original constitutive load and regime admission at the initial state.
    pub load_before: CreepingLoad,
    /// Original constitutive load and regime admission at the final state.
    pub load_after: CreepingLoad,
    /// Integrated translation and rotation dissipation [J], including any
    /// sourced Navier interface loss; no finite-body heat partition is inferred.
    pub viscous_heat_j: f64,
    /// Signed work maintaining the uniform background flow [J].
    pub background_work_j: f64,
    /// Integrated external force work [J], excluding buoyancy, gravity and drag.
    pub applied_work_j: f64,
    /// Actual delta(kinetic + body gravity + buoyancy potential) [J].
    pub stored_energy_change_j: f64,
    /// Stored change plus heat minus background and external-force work [J].
    /// Includes all represented-state rounding; numerical loss is never heat.
    pub energy_defect_j: f64,
    /// Original fluid diffusion time divided by the fastest body relaxation.
    pub diffusion_ratio: f64,
}

/// Exact affine-force motion with explicit admission of the complete slip path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExactCreepingAffineStep {
    /// Same original source, native motion and independent heat/work accounts.
    pub flow: ExactCreepingStep,
    /// Conservative model bound on COM speed relative to the held fluid [m/s].
    /// This may exceed the actual maximum; it is not an interval certificate.
    pub maximum_relative_speed_bound_m_s: f64,
    /// Diameter-based Reynolds number evaluated at that same complete-path bound.
    pub maximum_translation_reynolds_bound: f64,
    /// Fluid diffusion time / declared physical force-variation time.
    /// Zero for equal endpoints, where no varying-drive admission is required.
    pub drive_diffusion_ratio: f64,
}

impl<'a, F: FluidTransport + ?Sized> CreepingSphere<'a, F> {
    /// Exact unforced relaxation in this sphere's frozen uniform reservoir.
    ///
    /// Reuses source-resolved translation/rotation resistance, including a
    /// selected Navier boundary, and the original native mass/inertia/pose.
    /// Body weight and buoyancy act throughout the integrated settling path.
    /// The original midpoint methods and thermal adapters retain their behavior.
    pub fn step_exact(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        controls: ExactCreepingStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ExactCreepingStep, Error> {
        self.step_exact_with_wrench(body, state, Wrench::ZERO, controls, cancelled)
    }

    /// Exact frozen-reservoir motion with a held world force and free spin.
    ///
    /// Applied force excludes weight, buoyancy and fluid resistance. A nonzero
    /// applied couple or selected Rubinow--Keller lift refuses, since either
    /// changes the coupled path outside this exact reduction. Inertia must be
    /// exactly spherical; no near-spherical projection is made.
    ///
    /// Relative translation follows a straight segment parametrized by a
    /// monotone exponential, so its speed cannot exceed the larger endpoint
    /// speed. Free spin decreases monotonically. Checking both endpoint loads
    /// therefore admits the entire exact interval against the original Re
    /// limits; the fluid/body diffusion gate remains unchanged. This numerical
    /// admission is not an interval certificate or evidence for source validity.
    pub fn step_exact_with_wrench(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        applied: Wrench,
        controls: ExactCreepingStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExactCreepingStep, Error> {
        poll(&mut cancelled)?;
        self.admit_body(body)?;
        self.step_exact_frame_with_wrench(body.frame(), state, applied, controls, cancelled)
    }

    /// Integrate a prescribed affine world-force history in this SAME material fluid.
    ///
    /// The endpoints define the external force, excluding weight, buoyancy and
    /// drag. The initial body couple must be zero. Radius, original source
    /// snapshot, selected Navier interface and bath velocity remain held. Native
    /// mechanics integrates the ramp, gravity and the unchanged buoyancy force
    /// simultaneously with the original translation/rotation resistances.
    ///
    /// An affine drive can produce an interior slip maximum even when both
    /// endpoints pass. Admission therefore uses the native positive-convolution
    /// speed bound over the WHOLE interval, then the original source's diameter
    /// Reynolds calculation. The bound may refuse a coarse interval that source-
    /// consistent subdivision can admit. Free-spin speed decreases monotonically;
    /// the initial rotational Reynolds check and existing fluid/body diffusion
    /// gate still apply. The separately declared shortest physical force-variation
    /// scale must also pass that same diffusion allowance. This source-history
    /// admission is independent of numerical subdivision; endpoint samples alone
    /// cannot establish bandwidth. Fast drive refuses with `UnsteadyDrive` and
    /// belongs with the existing unsteady fluid owner.
    ///
    /// Equal force endpoints retain the held entry's sharper
    /// endpoint admission and identical motion/work receipt. Selected lift,
    /// torque and nonspherical inertia refuse.
    ///
    /// Applied work is the native integrated force work less the held buoyancy
    /// contribution; gravity and buoyancy potential remain in the original
    /// stored-energy account. Independent native viscous heat is never an energy
    /// remainder. Receipt rounding remains visible to the existing work-defect
    /// gate. These are authored source endpoints, not samples of an arbitrary
    /// feedback law. No accelerating bath, added mass, Basset history, finite
    /// fluid recoil, or broader hydrodynamic validity is inferred. A failed
    /// source, path-regime, work or cancellation check publishes no candidate.
    pub fn step_exact_with_affine_force(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        wrench_before: Wrench,
        force_after_world: Vec3,
        controls: ExactCreepingAffineStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExactCreepingAffineStep, Error> {
        poll(&mut cancelled)?;
        self.admit_body(body)?;
        self.step_exact_frame_with_affine_force(
            body.frame(),
            state,
            wrench_before,
            force_after_world,
            controls,
            cancelled,
        )
    }

    /// Same physical flow for an already admitted current spherical frame.
    /// Thermal owners retain geometry/source admission and momentum transfer.
    pub(crate) fn step_exact_frame_with_affine_force(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        wrench_before: Wrench,
        force_after_world: Vec3,
        controls: ExactCreepingAffineStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExactCreepingAffineStep, Error> {
        let ExactCreepingAffineStepControls {
            step: controls,
            minimum_force_variation_timescale_s,
        } = controls;
        if wrench_before.force_world == force_after_world {
            let flow = self.step_exact_frame_with_wrench(
                frame,
                state,
                wrench_before,
                controls,
                &mut cancelled,
            )?;
            // A held drive follows one straight velocity segment. Preserve the
            // existing sharper endpoint admission and its exact native receipt.
            let speed = |state: RigidBodyState| {
                norm(
                    state
                        .center_of_mass_velocity_world(frame.properties())?
                        .sub(self.ambient_velocity_world),
                )
            };
            let maximum_relative_speed_bound_m_s =
                speed(state)?.max(speed(flow.dynamics.state_after)?);
            let maximum_translation_reynolds_bound = flow
                .load_before
                .translation_reynolds
                .max(flow.load_after.translation_reynolds);
            poll(&mut cancelled)?;
            return Ok(ExactCreepingAffineStep {
                flow,
                maximum_relative_speed_bound_m_s,
                maximum_translation_reynolds_bound,
                drive_diffusion_ratio: 0.0,
            });
        }
        let drive_diffusion_ratio =
            self.exact_drive_diffusion_ratio(minimum_force_variation_timescale_s)?;
        vector(force_after_world)?;
        let prepared = self.prepare_exact_frame(
            frame,
            state,
            wrench_before,
            controls,
            &mut cancelled,
        )?;
        let buoyancy = prepared.load_before.buoyancy_force_world_n;
        let native_before = Wrench {
            force_world: vector(wrench_before.force_world.add(buoyancy))?,
            torque_body: Vec3::ZERO,
        };
        let native_after = vector(force_after_world.add(buoyancy))?;
        let maximum_relative_speed_bound_m_s = prepared
            .native
            .affine_force_relative_speed_bound_m_s(
                state,
                native_before,
                native_after,
                controls.duration_s,
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
            .step_affine_force(state, native_before, native_after, controls.duration_s)
            .map_err(native_error)?;
        let displacement = components(dynamics.displacement_world).map(dd);
        let applied_work = dd(dynamics.applied_work_j) - force_work(buoyancy, displacement);
        let flow = self.finish_exact_frame(
            frame,
            state,
            dynamics,
            applied_work,
            prepared,
            controls,
            &mut cancelled,
        )?;
        Ok(ExactCreepingAffineStep {
            flow,
            maximum_relative_speed_bound_m_s,
            maximum_translation_reynolds_bound,
            drive_diffusion_ratio,
        })
    }

    /// Admit a changing source history at the CURRENT radius, independently
    /// of numerical subdivision or rounded equality of local force samples.
    pub(crate) fn exact_drive_diffusion_ratio(
        &self,
        minimum_force_variation_timescale_s: f64,
    ) -> Result<f64, Error> {
        if !minimum_force_variation_timescale_s.is_finite()
            || minimum_force_variation_timescale_s <= 0.0
        {
            return Err(Error::Invalid(
                "a positive finite physical force-variation timescale is required",
            ));
        }
        let drive_diffusion_ratio =
            finite(self.diffusion_time_s / minimum_force_variation_timescale_s)?;
        if drive_diffusion_ratio > self.limits.maximum_diffusion_ratio {
            return Err(Error::UnsteadyDrive {
                diffusion_ratio: drive_diffusion_ratio,
                maximum: self.limits.maximum_diffusion_ratio,
            });
        }
        Ok(drive_diffusion_ratio)
    }

    /// Same native flow for a current, already admitted spherical configuration.
    /// Thermal adapters retain source/geometry admission and the momentum map.
    pub(crate) fn step_exact_frame(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        controls: ExactCreepingStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ExactCreepingStep, Error> {
        self.step_exact_frame_with_wrench(frame, state, Wrench::ZERO, controls, cancelled)
    }

    fn step_exact_frame_with_wrench(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        applied: Wrench,
        controls: ExactCreepingStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExactCreepingStep, Error> {
        let prepared =
            self.prepare_exact_frame(frame, state, applied, controls, &mut cancelled)?;
        let dynamics = prepared
            .native
            .step(
                state,
                Wrench {
                    force_world: vector(
                        applied.force_world.add(prepared.load_before.buoyancy_force_world_n),
                    )?,
                    torque_body: Vec3::ZERO,
                },
                controls.duration_s,
            )
            .map_err(native_error)?;
        let applied_work = force_work(
            applied.force_world,
            components(dynamics.displacement_world).map(dd),
        );
        self.finish_exact_frame(
            frame,
            state,
            dynamics,
            applied_work,
            prepared,
            controls,
            &mut cancelled,
        )
    }

    fn prepare_exact_frame(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        applied: Wrench,
        controls: ExactCreepingStepControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedExactStep, Error> {
        poll(cancelled)?;
        if self.selected_lift != SphereLift::None || applied.torque_body != Vec3::ZERO {
            return Err(Error::Invalid(
                "exact Stokes relaxation requires no selected lift and zero applied torque",
            ));
        }
        if !controls.duration_s.is_finite()
            || controls.duration_s <= 0.0
            || !controls.maximum_energy_defect_j.is_finite()
            || controls.maximum_energy_defect_j <= 0.0
        {
            return Err(Error::Invalid(
                "positive finite exact-step duration and work budget required",
            ));
        }
        let properties = frame.properties();
        let (_, diffusion_ratio) = self.timescales_frame(frame)?;
        let load_before = self.load_frame(frame, state, cancelled)?;
        let native = IsotropicLinearDrag::new(
            properties,
            self.gravity,
            self.ambient_velocity_world,
            self.resistance_kg_s,
            self.rotary_resistance_kg_m2_s,
        )
        .map_err(native_error)?;
        Ok(PreparedExactStep { native, load_before, diffusion_ratio })
    }

    #[allow(clippy::too_many_arguments)] // Original state, source receipt and independent work.
    fn finish_exact_frame(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        dynamics: LinearDragStep,
        applied_work: Dd,
        prepared: PreparedExactStep,
        controls: ExactCreepingStepControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ExactCreepingStep, Error> {
        let PreparedExactStep { load_before, diffusion_ratio, .. } = prepared;
        let properties = frame.properties();
        let load_after = self.load_frame(frame, dynamics.state_after, &mut *cancelled)?;
        let mut kinetic = RigidKineticEnergyChange::new();
        kinetic.add_body(properties, state, dynamics.state_after)?;
        let old_position = components(state.pose().position_world());
        let new_position = components(dynamics.state_after.pose().position_world());
        let gravity = components(self.gravity.acceleration_world());
        let effective_mass = dd(properties.mass()) - dd(self.fluid_mass_kg);
        let mut stored_change = kinetic.finish_retained()?;
        for i in 0..3 {
            stored_change = stored_change
                - effective_mass * dd(gravity[i]) * (dd(new_position[i]) - dd(old_position[i]));
        }
        let defect = finite(
            (stored_change + dd(dynamics.dissipated_energy_j)
                - dd(dynamics.background_work_j)
                - applied_work)
                .to_f64(),
        )?;
        if defect.abs() > controls.maximum_energy_defect_j {
            return Err(Error::EnergyDefect {
                defect_j: defect,
                maximum_j: controls.maximum_energy_defect_j,
            });
        }
        poll(cancelled)?;
        Ok(ExactCreepingStep {
            dynamics,
            load_before,
            load_after,
            viscous_heat_j: dynamics.dissipated_energy_j,
            background_work_j: dynamics.background_work_j,
            applied_work_j: finite(applied_work.to_f64())?,
            stored_energy_change_j: finite(stored_change.to_f64())?,
            energy_defect_j: defect,
            diffusion_ratio,
        })
    }
}

struct PreparedExactStep {
    native: IsotropicLinearDrag,
    load_before: CreepingLoad,
    diffusion_ratio: f64,
}

fn force_work(force: Vec3, displacement: [Dd; 3]) -> Dd {
    let force = components(force).map(dd);
    let mut work = Dd::ZERO;
    for i in 0..3 {
        work = work + force[i] * displacement[i];
    }
    work
}

fn native_error(error: LinearDragError) -> Error {
    match error {
        LinearDragError::Invalid(message) => Error::Invalid(message),
        LinearDragError::Dynamics(error) => Error::Dynamics(error),
    }
}
fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}
fn components(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
