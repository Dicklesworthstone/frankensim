//! Source thermal expansion updates the original fluid-loaded sphere.
//!
//! The sphere's reference material inventory, source strain and native momentum
//! remap determine its current radius and inertia. The original creeping-flow
//! owner supplies all drag, rotary resistance, buoyancy, regime admission and
//! midpoint dynamics. There is no parallel heat or motion integrator here.

use super::*;
use crate::creeping_flow::{CreepingLoad, CreepingSphere, CreepingStepControls};
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
    /// Native Stokes loss deposited in the external fluid reservoir [J].
    /// This is never also credited to the finite solid's caloric storage.
    pub viscous_heat_j: f64,
    /// Signed work supplied by maintaining the original background flow [J].
    pub background_work_j: f64,
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
    /// heat, maintained-flow work, configuration work and numerical defect [J].
    pub combined_residual_j: f64,
    /// Largest admitted fluid/body diffusion ratio, including updated endpoints.
    pub maximum_diffusion_ratio: f64,
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
    /// Endpoint Reynolds and fluid-diffusion checks use the UPDATED radius.
    /// Body mass is unchanged, so heating changes drag, torque, displacement,
    /// angular speed and subsequent acceleration without replacing matter.
    ///
    /// The reservoir's gravity must match the original mechanical integrator.
    /// Continue with the same reservoir: replacing it is an external model
    /// change whose exchange work is not inferred. This first-order split
    /// requires time refinement. The original no-slip, unbounded, fully immersed
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
        checkpoint()?;
        if sphere.gravity() != self.motion.integrator.gravity() {
            return Err(HeatedMotionError::Invalid(
                "fluid and mechanical gravity must agree",
            ));
        }
        if !controls.maximum_relaxation_fraction.is_finite()
            || controls.maximum_relaxation_fraction <= 0.0
            || controls.maximum_relaxation_fraction > 1.0
            || !controls.maximum_step_energy_defect_j.is_finite()
            || controls.maximum_step_energy_defect_j <= 0.0
        {
            return Err(HeatedMotionError::Invalid(
                "creeping motion needs a relaxation fraction in (0,1] and a positive step work tolerance",
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
            controls,
            viscous_heat_j: 0.0,
            background_work_j: 0.0,
            displacement_work_j: 0.0,
            maximum_diffusion_ratio: 0.0,
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
        let displacement_parameter_work_j = observer.displacement_work_j;
        let maximum_diffusion_ratio = observer
            .maximum_diffusion_ratio
            .max(before.diffusion_ratio)
            .max(after.diffusion_ratio);
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
                viscous_heat_j,
                background_work_j,
                buoyancy_potential_before_j: before.buoyancy_potential_j,
                buoyancy_potential_after_j: after.buoyancy_potential_j,
                parameter_work_j: work.parameter_j,
                absolute_parameter_work_j: work.absolute_parameter_j,
                displacement_parameter_work_j,
                numerical_energy_change_j: work.numerical_j,
                absolute_numerical_energy_change_j: work.absolute_numerical_j,
                combined_residual_j,
                maximum_diffusion_ratio,
            },
        ))
    }
}

struct CreepingSnapshot {
    radius_m: f64,
    load: CreepingLoad,
    buoyancy_potential_j: f64,
    diffusion_ratio: f64,
}

struct CreepingExpansion<'a, 'm, 'f, F: FluidTransport + ?Sized> {
    motion: RigidExpansion<'a>,
    sphere: &'m CreepingSphere<'f, F>,
    controls: HeatedCreepingControls,
    viscous_heat_j: f64,
    background_work_j: f64,
    displacement_work_j: f64,
    maximum_diffusion_ratio: f64,
}
impl<F: FluidTransport + ?Sized> Clone for CreepingExpansion<'_, '_, '_, F> {
    fn clone(&self) -> Self {
        Self {
            motion: self.motion.clone(),
            sphere: self.sphere,
            controls: self.controls,
            viscous_heat_j: self.viscous_heat_j,
            background_work_j: self.background_work_j,
            displacement_work_j: self.displacement_work_j,
            maximum_diffusion_ratio: self.maximum_diffusion_ratio,
        }
    }
}
impl<'f, F: FluidTransport + ?Sized> CreepingExpansion<'_, '_, 'f, F> {
    fn current_sphere(&self) -> Result<CreepingSphere<'f, F>, HeatedMotionError> {
        let region = &self.motion.expansion.configuration.regions()[0];
        let radius = self.sphere.radius_m() * (1.0 + region.expansion().free_linear_strain());
        let sphere = self.sphere.at_radius(radius)?;
        sphere.admit_geometry(region.geometry())?;
        Ok(sphere)
    }
    fn snapshot(
        &self,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<CreepingSnapshot, HeatedMotionError> {
        checkpoint()?;
        let sphere = self.current_sphere()?;
        let load = checked(checkpoint, |cancelled| {
            sphere.load_frame(&self.motion.principal, self.motion.state, cancelled)
        })?;
        let (_, diffusion_ratio) = sphere.timescales_frame(&self.motion.principal)?;
        Ok(CreepingSnapshot {
            radius_m: sphere.radius_m(),
            load,
            buoyancy_potential_j: sphere.buoyancy_potential_j(self.motion.state)?,
            diffusion_ratio,
        })
    }
}
impl<F: FluidTransport + ?Sized> DepositionObserver for CreepingExpansion<'_, '_, '_, F> {
    type Error = HeatedMotionError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.viscous_heat_j = 0.0;
        self.background_work_j = 0.0;
        self.displacement_work_j = 0.0;
        self.maximum_diffusion_ratio = 0.0;
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
            sphere.step_frame(
                &self.motion.principal,
                self.motion.state,
                CreepingStepControls {
                    duration_s: dt_s,
                    maximum_relaxation_fraction: self.controls.maximum_relaxation_fraction,
                    maximum_energy_defect_j: self.controls.maximum_step_energy_defect_j,
                },
                cancelled,
            )
        })?;
        let potential_before_remap = sphere.buoyancy_potential_j(step.dynamics.state_after)?;
        self.motion.finish_step(
            start_s,
            dt_s,
            old,
            thermal,
            step.dynamics.state_after,
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
        self.maximum_diffusion_ratio = self
            .maximum_diffusion_ratio
            .max(step.diffusion_ratio)
            .max(after.diffusion_ratio);
        checkpoint()?;
        Ok(())
    }
}
