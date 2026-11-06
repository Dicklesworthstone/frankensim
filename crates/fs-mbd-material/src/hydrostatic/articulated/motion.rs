//! Configuration-dependent buoyancy through the native articulated midpoint step.
use super::{
    ArticulatedFiniteTankLoad, ArticulatedHydrostaticError, ArticulatedHydrostaticLoad,
    PreparedArticulatedHydrostatics, PrescribedHydrostaticReservoir,
};
use crate::hydrostatic::{
    FiniteHydrostaticTank, FiniteLayeredHydrostaticTank, HydrostaticError, HydrostaticStepControls,
    finite, ga, poll,
};
use fs_ga::{Se3, Wrench};
use fs_material::fluid_transport::FluidTransport;
use fs_math::dd::Dd;
use fs_mbd::Gravity;
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
    FreeFloatingMotionStep, step_free_floating_with_forcing,
};
use std::cell::RefCell;

/// A refused step publishes no candidate mechanism, surface or liquid history.
#[derive(Clone, Debug, PartialEq)]
pub enum ArticulatedHydrostaticMotionError {
    /// Original geometry, material, kinematics or hydrostatic cancellation error.
    Load(ArticulatedHydrostaticError),
    /// Original articulated motion, state, joint limit or cancellation error.
    Motion(FreeFloatingMotionError),
    /// Kinetic + body-gravity + fluid-potential change minus actuator work [J].
    /// This is numerical error, not drag, heat or physical fluid dissipation.
    EnergyDefect { defect_j: f64, maximum_j: f64 },
}

impl core::fmt::Display for ArticulatedHydrostaticMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "articulated hydrostatic motion: {self:?}")
    }
}
impl std::error::Error for ArticulatedHydrostaticMotionError {}
impl From<ArticulatedHydrostaticError> for ArticulatedHydrostaticMotionError {
    fn from(error: ArticulatedHydrostaticError) -> Self {
        Self::Load(error)
    }
}
impl From<HydrostaticError> for ArticulatedHydrostaticMotionError {
    fn from(error: HydrostaticError) -> Self {
        Self::Load(error.into())
    }
}
impl From<FreeFloatingMotionError> for ArticulatedHydrostaticMotionError {
    fn from(error: FreeFloatingMotionError) -> Self {
        Self::Motion(error)
    }
}

/// One native trajectory step and its independently evaluated fluid storage.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedHydrostaticMotionStep<L> {
    /// Original Lie-midpoint/ABA motion, momentum, actuator and work accounts.
    /// Its external_work already contains both buoyancy and body-gravity work.
    pub dynamics: FreeFloatingMotionStep,
    /// Actual start-stage load used by the native predictor.
    pub load_before: L,
    /// Actual midpoint-stage load used for native acceleration and impulse.
    pub load_midpoint: L,
    /// Hydrostatics at the returned configuration, for domain and energy admission.
    /// This is not another force kick or another fluid-energy account.
    pub load_after: L,
    /// Endpoint fluid-potential change [J], counted once for a finite tank.
    pub fluid_potential_change_j: f64,
    /// Minus native uniform-gravity work along its relative COM displacement [J].
    /// This avoids subtracting large absolute world-position potentials.
    pub body_gravitational_potential_change_j: f64,
    /// Actual kinetic + both potential changes minus native actuator work [J].
    /// Passing a local absolute gate is not a global stability/accuracy proof.
    pub energy_defect_j: f64,
}

impl PreparedArticulatedHydrostatics<'_> {
    /// Advance a free-floating mechanism in a fixed prescribed reservoir.
    ///
    /// The existing native Lie-midpoint integrator evaluates buoyancy at its
    /// ACTUAL initial and midpoint base/joint configurations. Its ABA and
    /// whole-tree momentum reconstruction retain joint coupling and base recoil.
    /// The same stage loads supply its momentum and work; no frozen-force step,
    /// additional gravity, independent link integration or new solver is used.
    ///
    /// Held generalized forces are internal actuator efforts in original joint
    /// units/order. Their native midpoint work is an external energy supply to
    /// the mechanism-fluid account. Use zeros for passive motion. Body gravity
    /// is taken from the reservoir exactly once. The native fixed-root,
    /// joint/state/effort and motion-budget admissions remain authoritative.
    ///
    /// The endpoint is evaluated again for immersion, material-domain and total
    /// energy admission. Keep the same mechanism, envelopes and fluid sources
    /// between steps; model replacement is not a conservative motion event.
    /// No contacts, drag, added mass, waves or damping are introduced. The native
    /// method is second order for resolved smooth motion, not energy-exact or
    /// stiffly stable. Refine across changing immersion and near interfaces.
    ///
    /// Cancellation is shared with every native and hydrostatic stage. All
    /// inputs are borrowed; a late load, motion or energy refusal returns no step.
    pub fn step_free_floating<R: PrescribedHydrostaticReservoir + ?Sized>(
        &self,
        reservoir: &R,
        state: &FreeFloatingMotionState,
        generalized_force: &[f64],
        controls: HydrostaticStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        ArticulatedHydrostaticMotionStep<ArticulatedHydrostaticLoad>,
        ArticulatedHydrostaticMotionError,
    > {
        self.step_with_hydrostatic_load(
            reservoir.gravity(),
            state,
            generalized_force,
            controls,
            cancelled,
            |pose, q, cancel| self.load(reservoir, pose, q, cancel),
        )
    }

    /// Advance ALL links against one conserved finite liquid inventory.
    ///
    /// Each native force stage and the endpoint solve one shared waterline.
    /// The complete finite-liquid gravitational potential, including cross-link
    /// displacement terms, enters the absolute energy gate exactly once.
    /// Original stationary-tank bounds, floor/wall clearance, 64-link limit and
    /// overflow refusal apply at every sampled configuration. The caller owns
    /// nonintersection and timestep resolution between those samples.
    pub fn step_free_floating_in_finite_tank<F: FluidTransport + ?Sized>(
        &self,
        tank: &FiniteHydrostaticTank<'_, F>,
        state: &FreeFloatingMotionState,
        generalized_force: &[f64],
        controls: HydrostaticStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        ArticulatedHydrostaticMotionStep<ArticulatedFiniteTankLoad>,
        ArticulatedHydrostaticMotionError,
    > {
        self.step_with_hydrostatic_load(
            tank.gravity(),
            state,
            generalized_force,
            controls,
            cancelled,
            |pose, q, cancel| self.load_finite_tank(tank, pose, q, cancel),
        )
    }

    /// Advance against one shared stack of conserved immiscible liquids.
    ///
    /// Every cumulative interface includes all links at each actual native
    /// stage. Density-weighted forces, torques and the complete stack potential
    /// reuse the existing finite-layered owner, with its source/domain bounds.
    pub fn step_free_floating_in_finite_layered_tank(
        &self,
        tank: &FiniteLayeredHydrostaticTank<'_>,
        state: &FreeFloatingMotionState,
        generalized_force: &[f64],
        controls: HydrostaticStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        ArticulatedHydrostaticMotionStep<ArticulatedFiniteTankLoad>,
        ArticulatedHydrostaticMotionError,
    > {
        self.step_with_hydrostatic_load(
            tank.gravity(),
            state,
            generalized_force,
            controls,
            cancelled,
            |pose, q, cancel| self.load_finite_layered_tank(tank, pose, q, cancel),
        )
    }

    fn step_with_hydrostatic_load<L: StageLoad, C: FnMut() -> bool>(
        &self,
        gravity: Gravity,
        state: &FreeFloatingMotionState,
        generalized_force: &[f64],
        controls: HydrostaticStepControls,
        mut cancelled: C,
        mut evaluate: impl FnMut(Se3, &[f64], &mut C) -> Result<L, ArticulatedHydrostaticError>,
    ) -> Result<ArticulatedHydrostaticMotionStep<L>, ArticulatedHydrostaticMotionError> {
        poll(&mut cancelled)?;
        if !controls.duration_s.is_finite()
            || controls.duration_s <= 0.0
            || !controls.maximum_energy_defect_j.is_finite()
            || controls.maximum_energy_defect_j <= 0.0
        {
            return Err(HydrostaticError::Invalid(
                "positive finite timestep and absolute energy allowance required",
            )
            .into());
        }
        let held = vec![Wrench::default(); self.source.model().link_count()];
        let loads = FreeFloatingMotionLoads {
            generalized_force,
            external_body_wrench: &held,
            gravity_world: ga(gravity.acceleration_world()),
        };
        // The native forcing and cancellation callbacks are separate synchronous
        // closures. Each borrows this ONE cancellation source only for its call.
        let cancellation = RefCell::new(cancelled);
        let mut load_error = None;
        let mut load_before = None;
        let mut load_midpoint = None;
        let result = step_free_floating_with_forcing(
            self.source.model(),
            state,
            loads,
            controls.duration_s,
            |_, stage, wrenches, _| {
                let sample = match evaluate(
                    stage.base.world_from_base,
                    &stage.joint_position,
                    &mut *cancellation.borrow_mut(),
                ) {
                    Ok(sample) => sample,
                    Err(error) => {
                        load_error = Some(error);
                        return Err(FreeFloatingMotionError::InvalidLoad(
                            "hydrostatic stage refused",
                        ));
                    }
                };
                wrenches.copy_from_slice(sample.wrenches());
                if load_before.is_none() {
                    load_before = Some(sample);
                } else {
                    load_midpoint = Some(sample);
                }
                Ok(())
            },
            || (*cancellation.borrow_mut())(),
        );
        if let Some(error) = load_error {
            return Err(error.into());
        }
        let dynamics = result?;
        let load_before = load_before.expect("native start forcing stage");
        let load_midpoint = load_midpoint.expect("native midpoint forcing stage");
        let load_after = evaluate(
            dynamics.state_after.base.world_from_base,
            &dynamics.state_after.joint_position,
            &mut *cancellation.borrow_mut(),
        )?;
        let fluid_change =
            Dd::from_f64(load_after.potential()) - Dd::from_f64(load_before.potential());
        let fluid_potential_change_j = finite(fluid_change.to_f64())?;
        let body_gravitational_potential_change_j = finite(-dynamics.gravity_work)?;
        // Native kinetic polarization and relative COM work retain small
        // increments; never subtract rounded total energies with a common datum.
        let energy_defect_j = finite(
            (Dd::from_f64(dynamics.kinetic_energy_change) - Dd::from_f64(dynamics.gravity_work)
                + fluid_change
                - Dd::from_f64(dynamics.actuator_work))
            .to_f64(),
        )?;
        if energy_defect_j.abs() > controls.maximum_energy_defect_j {
            return Err(ArticulatedHydrostaticMotionError::EnergyDefect {
                defect_j: energy_defect_j,
                maximum_j: controls.maximum_energy_defect_j,
            });
        }
        poll(&mut *cancellation.borrow_mut())?;
        Ok(ArticulatedHydrostaticMotionStep {
            dynamics,
            load_before,
            load_midpoint,
            load_after,
            fluid_potential_change_j,
            body_gravitational_potential_change_j,
            energy_defect_j,
        })
    }
}

// Only the two existing hydrostatic receipts can enter the common step owner.
trait StageLoad {
    fn wrenches(&self) -> &[Wrench];
    fn potential(&self) -> f64;
}
impl StageLoad for ArticulatedHydrostaticLoad {
    fn wrenches(&self) -> &[Wrench] {
        &self.external_body_wrenches
    }
    fn potential(&self) -> f64 {
        self.pressure_potential_j
    }
}
impl StageLoad for ArticulatedFiniteTankLoad {
    fn wrenches(&self) -> &[Wrench] {
        &self.external_body_wrenches
    }
    fn potential(&self) -> f64 {
        self.tank_load.surface.liquid_gravitational_energy_change_j
    }
}
