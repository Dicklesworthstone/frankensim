//! Exact held Stokes relaxation through the original native mechanics owner.
//!
//! This explicit path retains the ordinary sphere's material/interface sources,
//! buoyancy and validity checks. It removes the midpoint relaxation restriction
//! for a frozen uniform bath and zero applied couple. It does not turn that
//! numerical improvement into a broader hydrodynamic applicability claim.

use super::{CreepingLoad, CreepingSphere, Error, SphereLift, finite, poll, vector};
use crate::principal::PrincipalBody;
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

/// Source-bound exact Stokes motion and separately retained work/heat accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExactCreepingStep {
    /// The original native body's exact held drag flow and numerical residuals.
    /// Its held force includes buoyancy and the caller's external force.
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
    /// Held external force work [J], excluding buoyancy, gravity and drag.
    pub applied_work_j: f64,
    /// Actual delta(kinetic + body gravity + buoyancy potential) [J].
    pub stored_energy_change_j: f64,
    /// Stored change plus heat minus background and external-force work [J].
    /// Includes all represented-state rounding; numerical loss is never heat.
    pub energy_defect_j: f64,
    /// Original fluid diffusion time divided by the fastest body relaxation.
    pub diffusion_ratio: f64,
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
        let frame = body.frame();
        let properties = frame.properties();
        let (_, diffusion_ratio) = self.timescales_frame(frame)?;
        let load_before = self.load_frame(frame, state, &mut cancelled)?;
        let native = IsotropicLinearDrag::new(
            properties,
            self.gravity,
            self.ambient_velocity_world,
            self.resistance_kg_s,
            self.rotary_resistance_kg_m2_s,
        )
        .map_err(native_error)?;
        let dynamics = native
            .step(
                state,
                Wrench {
                    force_world: vector(
                        applied.force_world.add(load_before.buoyancy_force_world_n),
                    )?,
                    torque_body: Vec3::ZERO,
                },
                controls.duration_s,
            )
            .map_err(native_error)?;
        let load_after = self.load_frame(frame, dynamics.state_after, &mut cancelled)?;
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
        let force = components(applied.force_world);
        let displacement = components(dynamics.displacement_world);
        let mut applied_work = Dd::ZERO;
        for i in 0..3 {
            applied_work = applied_work + dd(force[i]) * dd(displacement[i]);
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
        poll(&mut cancelled)?;
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
