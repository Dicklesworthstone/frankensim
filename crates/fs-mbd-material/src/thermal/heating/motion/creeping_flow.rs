//! Source thermal expansion updates the original fluid-loaded sphere.
//!
//! The sphere's reference material inventory, source strain and native momentum
//! remap determine its current radius and inertia. The original creeping-flow
//! owner supplies all drag, selected rotational lift, rotary resistance,
//! buoyancy and regime admission. The ordinary entry uses native midpoint
//! dynamics; the explicit exact entry uses native analytic frozen-radius drag.
//! Both retain the same thermal timeline and momentum-preserving expansion.

/// Source-bound heat exchange with the same surrounding fluid.
pub mod heat_exchange;

mod exact;
pub use exact::{HeatedCreepingAffineFlow, HeatedCreepingAffineForce, HeatedExactCreepingControls};

use super::*;
use crate::creeping_flow::exact::{
    ExactCreepingAffineStepControls, ExactCreepingFlowStepControls, ExactCreepingStepControls,
};
use crate::creeping_flow::{CreepingFlowError, CreepingLoad, CreepingSphere, CreepingStepControls};
use fs_material::fluid_transport::FluidTransport;

/// Controls of the original fluid-loaded step on each actual thermal interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatedCreepingControls {
    /// In (0,1]; dt/min(m/c, I/c_rot) cannot exceed this fraction. The existing
    /// heated-motion maximum timestep still applies; reduce it on refusal.
    pub maximum_relaxation_fraction: f64,
    /// Positive absolute allowance [J] for each original midpoint work defect.
    /// The heated-motion budget also bounds their sum of absolute values.
    pub maximum_step_energy_defect_j: f64,
}

/// One accepted heat, trajectory and external-fluid accounting window.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatedCreepingAudit {
    pub thermal: EnergyAudit,
    pub mechanical_before: DynamicsDiagnostics,
    pub mechanical_after: DynamicsDiagnostics,
    pub radius_before_m: f64,
    pub radius_after_m: f64,
    pub load_before: CreepingLoad,
    pub load_after: CreepingLoad,
    /// Actual authored ambient endpoints [m/s]. Carry the accepted endpoint
    /// into a continuing flow history; ordinary entries retain the frozen bath.
    pub ambient_velocity_before_world: fs_mbd::Vec3,
    pub ambient_velocity_after_world: fs_mbd::Vec3,
    /// Native bulk-viscous plus any Navier interface-friction loss in the
    /// external reservoir [J]; no solid/fluid heat partition is inferred.
    /// This is never also credited to the finite solid's caloric storage.
    pub viscous_heat_j: f64,
    /// Signed work supplied by maintaining the original background flow [J].
    pub background_work_j: f64,
    /// Signed displaced-fluid acceleration-pressure work [J], already included
    /// once in background_work_j. Zero for the original held-flow entries.
    pub ambient_pressure_work_j: f64,
    /// Integrated external COM-force work [J], excluding gravity, buoyancy,
    /// fluid resistance and configuration work. Zero for unforced entries.
    pub applied_work_j: f64,
    /// Buoyancy potential rho_fluid*V*g dot x_COM [J], with the same world-origin
    /// datum as native gravity. This unbounded reservoir has no pressure datum.
    pub buoyancy_potential_before_j: f64,
    pub buoyancy_potential_after_j: f64,
    /// Signed inertia-remap plus buoyancy-potential configuration work [J].
    pub parameter_work_j: f64,
    /// Absolute contributions summed before opposite signs can cancel [J].
    pub absolute_parameter_work_j: f64,
    /// The buoyancy-potential component of parameter_work_j [J]. This is a
    /// configuration account, not a model of pressure-volume expansion heat.
    pub displacement_parameter_work_j: f64,
    /// Signed fixed-configuration mechanical/fluid numerical defect [J].
    pub numerical_energy_change_j: f64,
    pub absolute_numerical_energy_change_j: f64,
    /// Thermal + mechanical + buoyancy storage + fluid loss, less external
    /// heat, maintained-flow and external-force work, configuration work and numerical defect [J].
    pub combined_residual_j: f64,
    /// Largest admitted fluid/body diffusion ratio, including updated endpoints.
    pub maximum_diffusion_ratio: f64,
    /// Largest fluid-diffusion / declared force-variation timescale ratio,
    /// including updated-radius endpoints. Zero for an unforced or held drive.
    pub maximum_drive_diffusion_ratio: f64,
    /// Largest fluid-diffusion / physical flow-variation timescale ratio,
    /// including updated-radius endpoints. Zero for a held ambient history.
    pub maximum_flow_diffusion_ratio: f64,
    /// Largest admitted translational Reynolds bound over mechanical intervals
    /// and updated-radius endpoints. Affine intervals use their complete-path
    /// conservative bound, which may exceed the actual peak.
    pub maximum_translation_reynolds_bound: f64,
}

impl<'a> HeatedRigidBody<'a> {
    /// Continue one uniformly heated spherical region in the same material
    /// fluid. Supply the sphere's ORIGINAL reference radius on every call.
    ///
    /// Every actual thermal interval runs the original frozen-radius midpoint
    /// step, then applies source-integrated expansion and the original
    /// momentum-preserving inertia update. Radius is rebuilt from the reference,
    /// not compounded. The same fluid snapshot supplies density and viscosity;
    /// no source is re-queried during this advance or between radius rebuilds.
    /// A sourced Navier boundary retains physical slip length, recomputes b/a
    /// at each radius, and checks its constant source's surface-temperature
    /// support at every updated endpoint. Slip length is not thermally scaled.
    /// Endpoint Reynolds and fluid-diffusion checks use the UPDATED radius.
    /// Body mass is unchanged, so heating changes drag, torque, displacement,
    /// angular speed and subsequent acceleration without replacing matter.
    /// A sphere selected with `with_rotational_lift` retains that model through
    /// every radius rebuild. The force uses the current radius cubed and the
    /// actual native world spin; its work belongs to the maintained-flow account,
    /// and it adds no viscous heat. This also applies to composed heat observers
    /// and the surrounding-fluid conduction adapter.
    ///
    /// The reservoir's gravity must match the original mechanical integrator.
    /// Continue with the same reservoir: replacing it is an external model
    /// change whose exchange work is not inferred. This first-order split
    /// requires time refinement. The original unbounded, fully immersed
    /// quasi-steady limits apply; expanding-surface flow, its dissipation and
    /// kinetic energy, added mass/history and reciprocal thermoelastic heat are
    /// not supplied. Thermal geometry remains the original heating owner's
    /// reference-configuration approximation. Failed source, regime, work or
    /// cancellation checks publish neither proposed history.
    pub fn advance_creeping_sphere<F: FluidTransport + ?Sized>(
        &self,
        sphere: &CreepingSphere<'_, F>,
        controls: HeatedCreepingControls,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, HeatedCreepingAudit), HeatedMotionError> {
        let observer = self.creeping_observer(sphere, controls, &mut checkpoint)?;
        let before = observer.snapshot(&mut checkpoint)?;
        let (thermal, observer, audit) =
            self.thermal
                .advance_observed(&observer, target_s, config, &mut checkpoint)?;
        self.accept_creeping(before, thermal, observer, audit, &mut checkpoint)
    }

    /// Compose an existing electrical, frictional or other heat/transport
    /// observer on these SAME thermal and fluid-loaded mechanical substeps.
    /// The original observer owns generated heat; Stokes loss remains an
    /// external-fluid transfer. All three histories commit together.
    pub fn advance_creeping_sphere_with<O: DepositionObserver, F: FluidTransport + ?Sized>(
        &self,
        observer: &O,
        sphere: &CreepingSphere<'_, F>,
        controls: HeatedCreepingControls,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, O, HeatedCreepingAudit), ObserverCouplingError<O::Error, HeatedMotionError>>
    {
        let motion = self
            .creeping_observer(sphere, controls, &mut checkpoint)
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

    fn creeping_observer<'m, 'f, F: FluidTransport + ?Sized>(
        &self,
        sphere: &'m CreepingSphere<'f, F>,
        controls: HeatedCreepingControls,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<CreepingExpansion<'a, 'm, 'f, F>, HeatedMotionError> {
        self.creeping_observer_with_stepper(sphere, CreepingStepper::Midpoint(controls), checkpoint)
    }

    fn creeping_observer_with_stepper<'m, 'f, F: FluidTransport + ?Sized>(
        &self,
        sphere: &'m CreepingSphere<'f, F>,
        stepper: CreepingStepper,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<CreepingExpansion<'a, 'm, 'f, F>, HeatedMotionError> {
        checkpoint()?;
        if sphere.gravity() != self.motion.integrator.gravity() {
            return Err(HeatedMotionError::Invalid(
                "fluid and mechanical gravity must agree",
            ));
        }
        stepper.validate()?;
        if matches!(
            stepper,
            CreepingStepper::Exact(_)
                | CreepingStepper::ExactAffine(..)
                | CreepingStepper::ExactAffineFlow(..)
        )
            && sphere.rotational_lift_enabled()
        {
            return Err(HeatedMotionError::Invalid(
                "exact creeping motion does not admit selected rotational lift",
            ));
        }
        if self.reference().regions().len() != 1 {
            return Err(HeatedMotionError::Invalid(
                "heated creeping motion requires one uniform spherical material region",
            ));
        }
        sphere.admit_geometry(self.reference().regions()[0].geometry())?;
        Ok(CreepingExpansion {
            motion: self.motion.clone(),
            sphere,
            stepper,
            viscous_heat_j: 0.0,
            background_work_j: 0.0,
            ambient_pressure_work_j: 0.0,
            applied_work_j: 0.0,
            displacement_work_j: 0.0,
            maximum_diffusion_ratio: 0.0,
            maximum_drive_diffusion_ratio: 0.0,
            maximum_flow_diffusion_ratio: 0.0,
            maximum_translation_reynolds_bound: 0.0,
        })
    }

    fn accept_creeping<F: FluidTransport + ?Sized>(
        &self,
        before: CreepingSnapshot,
        thermal: DepositionSimulation<'a>,
        observer: CreepingExpansion<'a, '_, '_, F>,
        thermal_audit: EnergyAudit,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, HeatedCreepingAudit), HeatedMotionError> {
        let after = observer.snapshot(checkpoint)?;
        let work = observer.motion.work;
        let viscous_heat_j = observer.viscous_heat_j;
        let background_work_j = observer.background_work_j;
        let ambient_pressure_work_j = observer.ambient_pressure_work_j;
        let applied_work_j = observer.applied_work_j;
        let displacement_parameter_work_j = observer.displacement_work_j;
        let maximum_diffusion_ratio = observer
            .maximum_diffusion_ratio
            .max(before.diffusion_ratio)
            .max(after.diffusion_ratio);
        let maximum_drive_diffusion_ratio = observer
            .maximum_drive_diffusion_ratio
            .max(before.drive_diffusion_ratio)
            .max(after.drive_diffusion_ratio);
        let maximum_flow_diffusion_ratio = observer
            .maximum_flow_diffusion_ratio
            .max(before.flow_diffusion_ratio)
            .max(after.flow_diffusion_ratio);
        let maximum_translation_reynolds_bound = observer
            .maximum_translation_reynolds_bound
            .max(before.load.translation_reynolds)
            .max(after.load.translation_reynolds);
        let next = Self {
            thermal,
            motion: observer.motion,
        };
        let mechanical_before = self.diagnostics()?;
        let mechanical_after = next.diagnostics()?;
        let combined_residual_j = finite(
            thermal_audit.residual_j
                + mechanical_change(mechanical_before, mechanical_after)?
                + (after.buoyancy_potential_j - before.buoyancy_potential_j)
                + viscous_heat_j
                - background_work_j
                - applied_work_j
                - work.parameter_j
                - work.numerical_j,
        )?;
        checkpoint()?;
        Ok((
            next,
            HeatedCreepingAudit {
                thermal: thermal_audit,
                mechanical_before,
                mechanical_after,
                radius_before_m: before.radius_m,
                radius_after_m: after.radius_m,
                load_before: before.load,
                load_after: after.load,
                ambient_velocity_before_world: before.ambient_velocity_world,
                ambient_velocity_after_world: after.ambient_velocity_world,
                viscous_heat_j,
                background_work_j,
                ambient_pressure_work_j,
                applied_work_j,
                buoyancy_potential_before_j: before.buoyancy_potential_j,
                buoyancy_potential_after_j: after.buoyancy_potential_j,
                parameter_work_j: work.parameter_j,
                absolute_parameter_work_j: work.absolute_parameter_j,
                displacement_parameter_work_j,
                numerical_energy_change_j: work.numerical_j,
                absolute_numerical_energy_change_j: work.absolute_numerical_j,
                combined_residual_j,
                maximum_diffusion_ratio,
                maximum_drive_diffusion_ratio,
                maximum_flow_diffusion_ratio,
                maximum_translation_reynolds_bound,
            },
        ))
    }
}

struct CreepingSnapshot {
    radius_m: f64,
    load: CreepingLoad,
    buoyancy_potential_j: f64,
    diffusion_ratio: f64,
    drive_diffusion_ratio: f64,
    flow_diffusion_ratio: f64,
    ambient_velocity_world: fs_mbd::Vec3,
}

struct CreepingExpansion<'a, 'm, 'f, F: FluidTransport + ?Sized> {
    motion: RigidExpansion<'a>,
    sphere: &'m CreepingSphere<'f, F>,
    stepper: CreepingStepper,
    viscous_heat_j: f64,
    background_work_j: f64,
    ambient_pressure_work_j: f64,
    applied_work_j: f64,
    displacement_work_j: f64,
    maximum_diffusion_ratio: f64,
    maximum_drive_diffusion_ratio: f64,
    maximum_flow_diffusion_ratio: f64,
    maximum_translation_reynolds_bound: f64,
}
impl<F: FluidTransport + ?Sized> Clone for CreepingExpansion<'_, '_, '_, F> {
    fn clone(&self) -> Self {
        Self {
            motion: self.motion.clone(),
            sphere: self.sphere,
            stepper: self.stepper,
            viscous_heat_j: self.viscous_heat_j,
            background_work_j: self.background_work_j,
            ambient_pressure_work_j: self.ambient_pressure_work_j,
            applied_work_j: self.applied_work_j,
            displacement_work_j: self.displacement_work_j,
            maximum_diffusion_ratio: self.maximum_diffusion_ratio,
            maximum_drive_diffusion_ratio: self.maximum_drive_diffusion_ratio,
            maximum_flow_diffusion_ratio: self.maximum_flow_diffusion_ratio,
            maximum_translation_reynolds_bound: self.maximum_translation_reynolds_bound,
        }
    }
}
impl<'f, F: FluidTransport + ?Sized> CreepingExpansion<'_, '_, 'f, F> {
    fn current_sphere(&self) -> Result<CreepingSphere<'f, F>, HeatedMotionError> {
        let region = &self.motion.expansion.configuration.regions()[0];
        let radius = self.sphere.radius_m() * (1.0 + region.expansion().free_linear_strain());
        let sphere = self.sphere.at_radius(radius)?;
        sphere.admit_geometry(region.geometry())?;
        sphere.admit_surface_material(region.expansion().current().resolved())?;
        Ok(sphere)
    }
    fn snapshot(
        &self,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<CreepingSnapshot, HeatedMotionError> {
        checkpoint()?;
        let sphere = self.current_sphere()?;
        let ambient_velocity_world = match self.stepper {
            CreepingStepper::ExactAffineFlow(_, _, flow) => flow.at(self.motion.time_s)?,
            _ => sphere.ambient_velocity_world(),
        };
        let load = checked(checkpoint, |cancelled| {
            match self.stepper {
                CreepingStepper::ExactAffineFlow(..) => sphere.load_frame_in_flow(
                    &self.motion.principal,
                    self.motion.state,
                    ambient_velocity_world,
                    cancelled,
                ),
                _ => sphere.load_frame(&self.motion.principal, self.motion.state, cancelled),
            }
        })?;
        let (_, diffusion_ratio) = sphere.timescales_frame(&self.motion.principal)?;
        let drive_diffusion_ratio = self.stepper.drive_diffusion_ratio(&sphere)?;
        let flow_diffusion_ratio = self.stepper.flow_diffusion_ratio(&sphere)?;
        Ok(CreepingSnapshot {
            radius_m: sphere.radius_m(),
            load,
            buoyancy_potential_j: sphere.buoyancy_potential_j(self.motion.state)?,
            diffusion_ratio,
            drive_diffusion_ratio,
            flow_diffusion_ratio,
            ambient_velocity_world,
        })
    }
}
impl<F: FluidTransport + ?Sized> DepositionObserver for CreepingExpansion<'_, '_, '_, F> {
    type Error = HeatedMotionError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.viscous_heat_j = 0.0;
        self.background_work_j = 0.0;
        self.ambient_pressure_work_j = 0.0;
        self.applied_work_j = 0.0;
        self.displacement_work_j = 0.0;
        self.maximum_diffusion_ratio = 0.0;
        self.maximum_drive_diffusion_ratio = 0.0;
        self.maximum_flow_diffusion_ratio = 0.0;
        self.maximum_translation_reynolds_bound = 0.0;
        self.motion.begin()
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.motion.event(event, thermal)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.motion
            .limit_timestep_s(thermal, proposed_s, checkpoint)
    }
    fn step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        if start_s != self.motion.time_s {
            return Err(HeatedMotionError::Invalid(
                "fluid and thermal clocks disagree",
            ));
        }
        let sphere = self.current_sphere()?;
        let step = checked(checkpoint, |cancelled| {
            self.stepper.step(
                &sphere,
                &self.motion.principal,
                self.motion.state,
                dt_s,
                [start_s, thermal.time_s()],
                cancelled,
            )
        })?;
        let potential_before_remap = sphere.buoyancy_potential_j(step.state_after)?;
        self.motion.finish_step(
            start_s,
            dt_s,
            old,
            thermal,
            step.state_after,
            step.energy_defect_j,
            checkpoint,
        )?;
        let after = self.snapshot(checkpoint)?;
        let displacement_work_j = finite(after.buoyancy_potential_j - potential_before_remap)?;
        self.displacement_work_j = finite(self.displacement_work_j + displacement_work_j)?;
        self.motion.work.parameter_j = finite(self.motion.work.parameter_j + displacement_work_j)?;
        self.motion.work.absolute_parameter_j =
            finite(self.motion.work.absolute_parameter_j + displacement_work_j.abs())?;
        if self.motion.work.absolute_parameter_j
            > self.motion.controls.transfer.maximum_parameter_work_j
        {
            return Err(HeatedMotionError::Budget(
                "total absolute configuration work including fluid displacement",
            ));
        }
        self.viscous_heat_j = finite(self.viscous_heat_j + step.viscous_heat_j)?;
        self.background_work_j = finite(self.background_work_j + step.background_work_j)?;
        self.ambient_pressure_work_j =
            finite(self.ambient_pressure_work_j + step.ambient_pressure_work_j)?;
        self.applied_work_j = finite(self.applied_work_j + step.applied_work_j)?;
        self.maximum_diffusion_ratio = self
            .maximum_diffusion_ratio
            .max(step.diffusion_ratio)
            .max(after.diffusion_ratio);
        self.maximum_drive_diffusion_ratio = self
            .maximum_drive_diffusion_ratio
            .max(after.drive_diffusion_ratio);
        self.maximum_flow_diffusion_ratio = self
            .maximum_flow_diffusion_ratio
            .max(after.flow_diffusion_ratio);
        self.maximum_translation_reynolds_bound = self
            .maximum_translation_reynolds_bound
            .max(step.maximum_translation_reynolds_bound)
            .max(after.load.translation_reynolds);
        checkpoint()?;
        Ok(())
    }
}

// Both methods pass their actual native state and independent fluid accounts
// through the same expansion, momentum transfer and thermal publication path.
#[derive(Clone, Copy)]
enum CreepingStepper {
    Midpoint(HeatedCreepingControls),
    Exact(HeatedExactCreepingControls),
    ExactAffine(HeatedExactCreepingControls, exact::AffineForceWindow),
    ExactAffineFlow(
        HeatedExactCreepingControls,
        exact::AffineForceWindow,
        exact::AffineFlowWindow,
    ),
}

struct CreepingInterval {
    state_after: RigidBodyState,
    energy_defect_j: f64,
    viscous_heat_j: f64,
    background_work_j: f64,
    ambient_pressure_work_j: f64,
    applied_work_j: f64,
    diffusion_ratio: f64,
    maximum_translation_reynolds_bound: f64,
}

impl CreepingStepper {
    fn validate(self) -> Result<(), HeatedMotionError> {
        let budget = match self {
            Self::Midpoint(controls) => {
                if !controls.maximum_relaxation_fraction.is_finite()
                    || controls.maximum_relaxation_fraction <= 0.0
                    || controls.maximum_relaxation_fraction > 1.0
                {
                    return Err(HeatedMotionError::Invalid(
                        "creeping motion needs a relaxation fraction in (0,1] and a positive step work tolerance",
                    ));
                }
                controls.maximum_step_energy_defect_j
            }
            Self::Exact(controls)
            | Self::ExactAffine(controls, _)
            | Self::ExactAffineFlow(controls, _, _) =>
                controls.maximum_step_energy_defect_j,
        };
        if !budget.is_finite() || budget <= 0.0 {
            return Err(HeatedMotionError::Invalid(
                "creeping motion needs a positive finite step work tolerance",
            ));
        }
        Ok(())
    }

    fn drive_diffusion_ratio<F: FluidTransport + ?Sized>(
        self,
        sphere: &CreepingSphere<'_, F>,
    ) -> Result<f64, CreepingFlowError> {
        match self {
            Self::ExactAffine(_, force) | Self::ExactAffineFlow(_, force, _) =>
                force.drive_diffusion_ratio(sphere),
            Self::Midpoint(_) | Self::Exact(_) => Ok(0.0),
        }
    }

    fn flow_diffusion_ratio<F: FluidTransport + ?Sized>(
        self,
        sphere: &CreepingSphere<'_, F>,
    ) -> Result<f64, CreepingFlowError> {
        match self {
            Self::ExactAffineFlow(_, _, flow) => flow.diffusion_ratio(sphere),
            _ => Ok(0.0),
        }
    }

    fn step<F: FluidTransport + ?Sized>(
        self,
        sphere: &CreepingSphere<'_, F>,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        duration_s: f64,
        time_window_s: [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingInterval, CreepingFlowError> {
        match self {
            Self::Midpoint(controls) => {
                let step = sphere.step_frame(
                    frame,
                    state,
                    CreepingStepControls {
                        duration_s,
                        maximum_relaxation_fraction: controls.maximum_relaxation_fraction,
                        maximum_energy_defect_j: controls.maximum_step_energy_defect_j,
                    },
                    cancelled,
                )?;
                Ok(CreepingInterval {
                    state_after: step.dynamics.state_after,
                    energy_defect_j: step.energy_defect_j,
                    viscous_heat_j: step.viscous_heat_j,
                    background_work_j: step.background_work_j,
                    ambient_pressure_work_j: 0.0,
                    applied_work_j: step.applied_work_j,
                    diffusion_ratio: step.diffusion_ratio,
                    maximum_translation_reynolds_bound: step.load_before.translation_reynolds
                        .max(step.load_after.translation_reynolds),
                })
            }
            Self::Exact(controls) => {
                let step = sphere.step_exact_frame(
                    frame,
                    state,
                    ExactCreepingStepControls {
                        duration_s,
                        maximum_energy_defect_j: controls.maximum_step_energy_defect_j,
                    },
                    cancelled,
                )?;
                Ok(CreepingInterval {
                    state_after: step.dynamics.state_after,
                    energy_defect_j: step.energy_defect_j,
                    viscous_heat_j: step.viscous_heat_j,
                    background_work_j: step.background_work_j,
                    ambient_pressure_work_j: 0.0,
                    applied_work_j: step.applied_work_j,
                    diffusion_ratio: step.diffusion_ratio,
                    maximum_translation_reynolds_bound: step.load_before.translation_reynolds
                        .max(step.load_after.translation_reynolds),
                })
            }
            Self::ExactAffine(controls, force) => {
                let [before, after] = force.endpoints(time_window_s)?;
                let affine = sphere.step_exact_frame_with_affine_force(
                    frame,
                    state,
                    Wrench {
                        force_world: before,
                        torque_body: fs_mbd::Vec3::ZERO,
                    },
                    after,
                    ExactCreepingAffineStepControls {
                        step: ExactCreepingStepControls {
                            duration_s,
                            maximum_energy_defect_j: controls.maximum_step_energy_defect_j,
                        },
                        minimum_force_variation_timescale_s:
                            force.minimum_variation_timescale_s(),
                    },
                    cancelled,
                )?;
                let step = affine.flow;
                Ok(CreepingInterval {
                    state_after: step.dynamics.state_after,
                    energy_defect_j: step.energy_defect_j,
                    viscous_heat_j: step.viscous_heat_j,
                    background_work_j: step.background_work_j,
                    ambient_pressure_work_j: 0.0,
                    applied_work_j: step.applied_work_j,
                    diffusion_ratio: step.diffusion_ratio,
                    maximum_translation_reynolds_bound:
                        affine.maximum_translation_reynolds_bound,
                })
            }
            Self::ExactAffineFlow(controls, force, flow) => {
                let [before, after] = force.endpoints(time_window_s)?;
                let ambient_before = flow.at(time_window_s[0])?;
                let ambient_after = flow.at(time_window_s[1])?;
                let affine = sphere.step_exact_frame_with_affine_force_and_flow(
                    frame,
                    state,
                    Wrench {
                        force_world: before,
                        torque_body: fs_mbd::Vec3::ZERO,
                    },
                    after,
                    ambient_before,
                    ambient_after,
                    ExactCreepingFlowStepControls {
                        step: ExactCreepingStepControls {
                            duration_s,
                            maximum_energy_defect_j: controls.maximum_step_energy_defect_j,
                        },
                        minimum_force_variation_timescale_s:
                            force.minimum_variation_timescale_s(),
                        minimum_flow_variation_timescale_s:
                            flow.minimum_variation_timescale_s(),
                    },
                    cancelled,
                )?;
                let step = affine.flow;
                Ok(CreepingInterval {
                    state_after: step.dynamics.state_after,
                    energy_defect_j: step.energy_defect_j,
                    viscous_heat_j: step.viscous_heat_j,
                    background_work_j: step.background_work_j,
                    ambient_pressure_work_j: affine.ambient_pressure_work_j,
                    applied_work_j: step.applied_work_j,
                    diffusion_ratio: step.diffusion_ratio,
                    maximum_translation_reynolds_bound:
                        affine.maximum_translation_reynolds_bound,
                })
            }
        }
    }
}
