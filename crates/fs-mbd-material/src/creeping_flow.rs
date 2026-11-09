//! A material-defined rigid sphere in an unbounded Newtonian fluid.
//!
//! No-slip Stokes translation gives F=-6*pi*mu*a*(v-U); rotation gives
//! torque=-8*pi*mu*a^3*omega. An explicit source-bound [`navier_slip`]
//! interface changes both resistances. Fully immersed buoyancy is -rho*V*g.
//! The shared fluid snapshot supplies rho and mu; the original material body
//! supplies mass and inertia. Native rigid-body mechanics owns all momentum,
//! position and attitude updates. The ordinary entry closes its midpoint load
//! relation; [`exact`] selects native exact motion for authored force/flow histories.
//!
//! References: R. Fitzpatrick, *Axisymmetric Stokes Flow Around a Solid Sphere*,
//! equations (10.94), (10.125)-(10.128),
//! <https://farside.ph.utexas.edu/teaching/336L/Fluid/node142.html>;
//! D. Tong, *Fluid Dynamics I and II: Questions on Viscous Flows*, rotating
//! sphere problem, <https://www.damtp.cam.ac.uk/user/tong/books/examples/f2.pdf>.
//!
//! This is an incompressible, quasi-steady reduction. It refuses large
//! Reynolds numbers and body relaxation too fast for the caller's declared
//! fluid-diffusion allowance. No Basset history, added mass, rarefaction,
//! free surface, walls, other particles or nonlinear drag is supplied. The
//! material's continuum/Newtonian validity remains an external model obligation.
//! Do not superpose another exterior drag law on this same sphere.
//! The explicit [`finite_reynolds`] and [`unsteady`] extensions retain their
//! own applicability domains; these Stokes-only limits describe [`CreepingSphere`].
//! The explicit rotational-lift methods add the leading Rubinow--Keller force
//! for a spinning sphere in uniform flow, within the same small-Reynolds limits.

use crate::principal::{PrincipalBody, PrincipalFrame};
use crate::{MaterialInertiaError, VolumeMoments};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_mbd::{
    DynamicsError, Gravity, RigidBodyIntegrator, RigidBodyState, StepReceipt, Vec3, Wrench,
};

/// Exact prescribed force/flow Stokes motion and free-spin relaxation.
pub mod exact;
pub use exact::{
    ExactCreepingAffineStep, ExactCreepingAffineStepControls, ExactCreepingFlowStep,
    ExactCreepingFlowStepControls, ExactCreepingStep, ExactCreepingStepControls,
};
/// Nonrotating, quasi-steady sphere motion beyond the Stokes regime.
pub mod finite_reynolds;
/// Uniform Navier slip from the original ordered solid/fluid interface source.
pub mod navier_slip;
/// Added fluid inertia and causal viscous memory on the same native body.
pub mod unsteady;
use navier_slip::{NavierSlipBoundary, ResistanceFactors};

/// Source, geometry, regime, native dynamics or acceptance failure.
#[derive(Clone, Debug, PartialEq)]
pub enum CreepingFlowError {
    /// Invalid controls, incompatible sphere geometry or unrepresentable values.
    Invalid(&'static str),
    /// Original fluid source refused its shared transport snapshot.
    Material(MaterialError),
    /// Original material-geometry owner refused the sphere.
    Geometry(MaterialInertiaError),
    /// Original rigid-body owner refused a state or step.
    Dynamics(DynamicsError),
    /// Translation or rotation left the declared creeping-flow envelope.
    Reynolds {
        /// Reynolds number based on the diameter and relative translation.
        translation: f64,
        /// Reynolds number based on radius squared and angular speed.
        rotation: f64,
        /// Caller-declared maximum for both dimensionless numbers.
        maximum: f64,
    },
    /// Fluid diffusion is too slow relative to the body's fastest relaxation.
    Unsteady {
        /// Fluid diffusion time divided by the fastest body relaxation time.
        diffusion_ratio: f64,
        /// Caller-declared maximum dimensionless ratio.
        maximum: f64,
    },
    /// Prescribed force or ambient flow varies too fast for quasi-steady response.
    UnsteadyDrive {
        /// Fluid diffusion time / caller-declared shortest source-variation time.
        diffusion_ratio: f64,
        /// Existing caller-declared maximum diffusion ratio.
        maximum: f64,
    },
    /// Independent kinetic/potential/viscous-work balance failed.
    EnergyDefect {
        /// Signed mechanical energy plus heat minus external flow work [J].
        defect_j: f64,
        /// Caller-declared absolute tolerance [J].
        maximum_j: f64,
    },
    /// Native body plus retained spectral-fluid impulse balance failed.
    MomentumDefect { residual_ns: Vec3, maximum_ns: f64 },
    /// Native world angular impulse plus retained rotary-fluid impulse failed.
    AngularMomentumDefect {
        residual_nms: Vec3,
        maximum_nms: f64,
    },
    /// Fluid history and the actual native midpoint angular velocity did not close.
    RotationalCoupling {
        relative_residual: f64,
        tolerance: f64,
        iterations: usize,
    },
    /// Cancellation before publication; accepted state remains untouched.
    Cancelled,
}
impl core::fmt::Display for CreepingFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material creeping sphere: {self:?}")
    }
}
impl std::error::Error for CreepingFlowError {}
impl From<MaterialError> for CreepingFlowError {
    fn from(e: MaterialError) -> Self {
        Self::Material(e)
    }
}
impl From<MaterialInertiaError> for CreepingFlowError {
    fn from(e: MaterialInertiaError) -> Self {
        Self::Geometry(e)
    }
}
impl From<DynamicsError> for CreepingFlowError {
    fn from(e: DynamicsError) -> Self {
        Self::Dynamics(e)
    }
}
type Error = CreepingFlowError;

/// Explicit approximation limits, not certified error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CreepingFlowLimits {
    /// In (0,0.1]. Uses 2*a*|v-U|/nu and a^2*|omega|/nu.
    pub maximum_reynolds: f64,
    /// In (0,0.1]. The ratio (a^2/nu)/min(m/c, I/c_rot) must fit.
    /// This conservative gate favors heavy bodies in a light fluid. It does
    /// not certify neglect of the unsteady initial fluid boundary layer.
    pub maximum_diffusion_ratio: f64,
}

/// Numerical controls for the original rigid-body midpoint step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CreepingStepControls {
    /// Positive physical interval [s].
    pub duration_s: f64,
    /// In (0,1]; dt/min(m/c,I/c_rot) cannot exceed this value. Refinement is
    /// still needed; the energy balance alone does not establish time accuracy.
    pub maximum_relaxation_fraction: f64,
    /// Positive absolute allowance for the independent work residual [J].
    pub maximum_energy_defect_j: f64,
}

/// Immutable isothermal reservoir and explicitly declared spherical boundary.
/// Source receipts remain reachable through `source`; no density/viscosity
/// defaults or numerical damping coefficients are introduced.
#[derive(Debug)]
pub struct CreepingSphere<'a, F: FluidTransport + ?Sized> {
    source: &'a F,
    transport: FluidTransportState,
    sphere: VolumeMoments,
    radius_m: f64,
    navier_slip: Option<NavierSlipBoundary<'a>>,
    resistance_kg_s: f64,
    rotary_resistance_kg_m2_s: f64,
    fluid_mass_kg: f64,
    diffusion_time_s: f64,
    gravity: Gravity,
    ambient_velocity_world: Vec3,
    limits: CreepingFlowLimits,
    selected_lift: SphereLift,
}

/// Instantaneous physical load at the actual supplied velocity, with body
/// gravity excluded. Bulk viscous and any Navier interface-friction loss go
/// to the prescribed external reservoir;
/// it is not automatically added to any finite material heat store.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CreepingLoad {
    /// Fluid drag, selected lift, buoyancy and rotary resistance; excludes body weight.
    pub wrench: Wrench,
    /// Translational resistance in the world frame [N].
    pub drag_force_world_n: Vec3,
    /// Transverse Rubinow--Keller force [N]; zero unless rotational lift is selected.
    /// It does no work against the relative translational velocity.
    pub lift_force_world_n: Vec3,
    /// Rotary resistance in the principal body frame [N m].
    pub rotary_torque_body_nm: Vec3,
    /// Fully immersed Archimedes force in the world frame [N].
    pub buoyancy_force_world_n: Vec3,
    /// Nonnegative total bulk-viscous plus interface-friction loss [W].
    /// The external-reservoir account does not infer a solid/fluid heat partition.
    pub dissipated_power_w: f64,
    /// Work supplied by maintaining the admitted background flow. The uniform
    /// sphere uses (F_drag + F_lift) dot U; an ambient-rotation owner also includes
    /// its rotary drag moment dotted with the local fluid angular velocity.
    pub background_power_w: f64,
    /// Diameter-based Reynolds number at the supplied relative velocity.
    pub translation_reynolds: f64,
    /// Radius-squared Reynolds number at the supplied angular velocity.
    pub rotation_reynolds: f64,
}

/// One atomic candidate, retaining the exact native motion receipt.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CreepingStep {
    /// Candidate pose, momenta and diagnostics from the existing integrator.
    pub dynamics: StepReceipt,
    /// Fluid load evaluated at the accepted starting state.
    pub load_before: CreepingLoad,
    /// Fluid load at the algebraically closed native midpoint velocity.
    pub load_midpoint: CreepingLoad,
    /// Fluid load evaluated at the actual native candidate endpoint.
    pub load_after: CreepingLoad,
    /// Total bulk-viscous plus interface-friction energy transferred to the
    /// external reservoir [J], without a finite-body heat-partition model.
    pub viscous_heat_j: f64,
    /// Energy supplied by maintaining the uniform background flow [J].
    pub background_work_j: f64,
    /// Work of the separately supplied, held external wrench [J]. Gravity
    /// and the fluid wrench are excluded so their work is never counted twice.
    pub applied_work_j: f64,
    /// Delta(kinetic + body gravity + buoyancy potential) [J].
    pub stored_energy_change_j: f64,
    /// Stored-energy change plus heat minus background and applied work [J].
    pub energy_defect_j: f64,
    /// Fluid diffusion time divided by the fastest body relaxation time.
    pub diffusion_ratio: f64,
}

impl<'a, F: FluidTransport + ?Sized> CreepingSphere<'a, F> {
    /// Declare a fully immersed no-slip sphere in a spatially uniform, frozen fluid.
    /// Its material body must have one uniform spherical region with these
    /// same volume moments. A geometry declaration is not inferred from a
    /// material name, and matching moments alone do not certify actual shape.
    pub fn new(
        source: &'a F,
        radius_m: f64,
        gravity: Gravity,
        ambient_velocity_world: Vec3,
        limits: CreepingFlowLimits,
    ) -> Result<Self, Error> {
        Self::with_transport(
            source,
            source.transport_state()?,
            radius_m,
            gravity,
            ambient_velocity_world,
            limits,
            None,
        )
    }

    fn with_transport(
        source: &'a F,
        transport: FluidTransportState,
        radius_m: f64,
        gravity: Gravity,
        ambient_velocity_world: Vec3,
        limits: CreepingFlowLimits,
        navier_slip: Option<NavierSlipBoundary<'a>>,
    ) -> Result<Self, Error> {
        let sphere = VolumeMoments::sphere(radius_m)?;
        if !ambient_velocity_world.is_finite()
            || !small_positive(limits.maximum_reynolds, 0.1)
            || !small_positive(limits.maximum_diffusion_ratio, 0.1)
        {
            return Err(Error::Invalid(
                "finite ambient velocity and limits in (0,0.1] required",
            ));
        }
        let resistance_kg_s =
            positive(6.0 * core::f64::consts::PI * transport.dynamic_viscosity_pa_s() * radius_m)?;
        let rotary_resistance_kg_m2_s =
            positive((4.0 / 3.0) * resistance_kg_s * radius_m * radius_m)?;
        // Derive both coefficients from their ORIGINAL no-slip expressions.
        // Scaling the already corrected translation coefficient would apply
        // the wrong Gamma_2,3 factor to rotary resistance. Keep b=0 arithmetic
        // exactly unchanged, including an explicitly sourced no-slip limit.
        let (resistance_kg_s, rotary_resistance_kg_m2_s) =
            if let Some(boundary) = navier_slip.filter(|boundary| boundary.slip_length_m() > 0.0) {
                let factors = ResistanceFactors::new(radius_m, boundary.slip_length_m()).ok_or(
                    Error::Invalid("Navier spherical resistance is unrepresentable"),
                )?;
                (
                    positive(resistance_kg_s * factors.translation)?,
                    positive(rotary_resistance_kg_m2_s * factors.rotation)?,
                )
            } else {
                (resistance_kg_s, rotary_resistance_kg_m2_s)
            };
        let fluid_mass_kg = positive(transport.density_kg_m3() * sphere.volume_m3())?;
        let diffusion_time_s =
            positive(radius_m * (radius_m / transport.kinematic_viscosity_m2_s()))?;
        Ok(Self {
            source,
            transport,
            sphere,
            radius_m,
            navier_slip,
            resistance_kg_s,
            rotary_resistance_kg_m2_s,
            fluid_mass_kg,
            diffusion_time_s,
            gravity,
            ambient_velocity_world,
            limits,
            selected_lift: SphereLift::None,
        })
    }

    /// Select a homogeneous impermeable Navier boundary from its actual
    /// solid/fluid interface source. No material or transport snapshot is
    /// queried again. Physical b stays fixed when thermal expansion changes a;
    /// translation uses (1+2b/a)/(1+3b/a), rotation uses 1/(1+3b/a).
    /// The original Reynolds, diffusion, step and work gates still apply.
    ///
    /// Positive slip is incompatible with the existing no-slip Basset,
    /// rotational-memory and Rubinow--Keller extensions. A selected rotational
    /// lift law is retained for b=0 and causes positive b to refuse here; it
    /// is never silently reset. This adds no rarefaction, permeability or bubble law.
    pub fn with_navier_slip(self, boundary: NavierSlipBoundary<'a>) -> Result<Self, Error> {
        let mut updated = Self::with_transport(
            self.source,
            self.transport,
            self.radius_m,
            self.gravity,
            self.ambient_velocity_world,
            self.limits,
            Some(boundary),
        )?;
        updated.selected_lift = self.selected_lift;
        if updated.rotational_lift_enabled() {
            updated.require_no_slip()?;
        }
        Ok(updated)
    }

    /// Select the leading quasi-steady rotational lift for every ordinary
    /// load/step and the heated-sphere adapters that consume this reservoir.
    ///
    /// Construction is Stokes-only by default. This explicit selection uses
    /// the same law and native midpoint solve as [`Self::step_with_rotational_lift`].
    /// Thermal radius rebuilds retain it together with the frozen fluid snapshot;
    /// the current radius and native spin determine the actual lift at each stage.
    /// This does not model expansion-induced fluid motion or make the law valid
    /// with unsteady added-mass/history dynamics; those entry points refuse it.
    /// Selecting this on a positive-slip Navier sphere makes ordinary loads,
    /// steps and thermal adapters refuse the incompatible no-slip lift law.
    #[must_use]
    pub fn with_rotational_lift(mut self) -> Self {
        self.selected_lift = SphereLift::RubinowKeller;
        self
    }

    /// Whether ordinary loads, steps and thermal adapters include rotational lift.
    pub const fn rotational_lift_enabled(&self) -> bool {
        matches!(self.selected_lift, SphereLift::RubinowKeller)
    }

    /// Rebuild only the spherical boundary, retaining the SAME frozen fluid
    /// snapshot. Thermal configuration adapters must not re-query a mutable
    /// external source halfway through an accepted reservoir interval.
    pub(crate) fn at_radius(&self, radius_m: f64) -> Result<Self, Error> {
        let mut resized = Self::with_transport(
            self.source,
            self.transport,
            radius_m,
            self.gravity,
            self.ambient_velocity_world,
            self.limits,
            self.navier_slip,
        )?;
        resized.selected_lift = self.selected_lift;
        Ok(resized)
    }

    /// Retain the original fluid owner and its provenance or mixture receipts.
    pub const fn source(&self) -> &'a F {
        self.source
    }
    /// Shared fluid snapshot resolved once when this reservoir was declared.
    pub const fn transport(&self) -> FluidTransportState {
        self.transport
    }
    /// Declared spherical boundary radius [m], independent of slip length.
    pub const fn radius_m(&self) -> f64 {
        self.radius_m
    }
    /// Original optional interface source. Absence is the explicit no-slip model.
    pub const fn navier_slip(&self) -> Option<NavierSlipBoundary<'a>> {
        self.navier_slip
    }
    /// The same gravity used for body weight and displaced-fluid buoyancy.
    pub const fn gravity(&self) -> Gravity {
        self.gravity
    }

    /// Prepared uniform reservoir velocity [m/s]. Explicit history consumers
    /// carry their accepted ambient endpoints independently of this value.
    pub const fn ambient_velocity_world(&self) -> Vec3 {
        self.ambient_velocity_world
    }

    fn admit_body(&self, body: &PrincipalBody<'_>) -> Result<(), Error> {
        let regions = body.source().regions();
        if regions.len() != 1 {
            return Err(Error::Invalid(
                "creeping sphere requires one uniform spherical material region",
            ));
        }
        self.admit_surface_material(regions[0].material())?;
        self.admit_geometry(regions[0].geometry())
    }

    /// Thermal adapters pass the current source-backed solid state, not just
    /// updated moments, so source identity and surface-temperature support hold.
    pub(crate) fn admit_surface_material(
        &self,
        solid: &fs_material::state_point::ResolvedMaterialStatePoint,
    ) -> Result<(), Error> {
        if let Some(boundary) = self.navier_slip {
            boundary.require_solid(solid)?;
        }
        Ok(())
    }

    fn require_no_slip(&self) -> Result<(), Error> {
        if self
            .navier_slip
            .is_some_and(|boundary| boundary.slip_length_m() > 0.0)
        {
            return Err(Error::Invalid(
                "this sphere extension derives from a no-slip boundary and does not admit Navier slip",
            ));
        }
        Ok(())
    }

    /// Configuration adapters establish the single uniform region and its
    /// material inventory; this compares its current volume and full moments.
    pub(crate) fn admit_geometry(&self, geometry: VolumeMoments) -> Result<(), Error> {
        let expected = self.sphere.inertia_m5().m[0];
        if (geometry.volume_m3() / self.sphere.volume_m3() - 1.0).abs() > 1e-10
            || geometry
                .inertia_m5()
                .m
                .iter()
                .zip(self.sphere.inertia_m5().m)
                .any(|(actual, target)| (actual - target).abs() > 1e-10 * expected)
        {
            return Err(Error::Invalid(
                "declared sphere and source body's volume/inertia disagree",
            ));
        }
        Ok(())
    }

    /// Evaluate drag, selected lift, rotary resistance and buoyancy without advancing
    /// physical state. Original body weight is supplied ONCE by native gravity.
    pub fn load(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingLoad, Error> {
        poll(&mut cancelled)?;
        self.admit_body(body)?;
        self.load_frame(body.frame(), state, &mut cancelled)
    }

    /// Add the leading small-Reynolds rotational lift to Stokes drag and buoyancy.
    ///
    /// Rubinow & Keller (1961), <https://doi.org/10.1017/S0022112061000640>:
    /// F_lift = pi*rho*a^3 * omega_world cross (v_body - U). The sphere's
    /// spin and translational slip are expressed in the SAME world frame.
    /// Fluid density comes from the original transport snapshot, not body mass.
    ///
    /// This is an explicit quasi-steady, unbounded, uniform-flow asymptote,
    /// not a high-Reynolds spinning-ball law or a shear/Saffman correction.
    /// Requires the no-slip boundary. Both existing Reynolds ceilings still apply. Viscous dissipation is
    /// unchanged; lift contributes to maintained-flow work when U is nonzero.
    pub fn load_with_rotational_lift(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingLoad, Error> {
        poll(&mut cancelled)?;
        self.admit_body(body)?;
        self.load_frame_with_lift(body.frame(), state, SphereLift::RubinowKeller, cancelled)
    }

    /// Same force owner for an already admitted spherical configuration.
    pub(crate) fn load_frame(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingLoad, Error> {
        self.load_frame_with_lift(frame, state, self.selected_lift, cancelled)
    }

    /// Undisturbed uniform-flow pressure force, shared by exact quasi-steady
    /// and unsteady motion. Added inertia and history are separate owners.
    fn ambient_pressure_force(&self, acceleration_world_m_s2: Vec3) -> Result<Vec3, Error> {
        vector(acceleration_world_m_s2)?;
        vector(acceleration_world_m_s2.scale(self.fluid_mass_kg))
    }

    /// Shared instantaneous Stokes admission for a time-evolving uniform flow.
    /// Acceleration pressure and any history are accounted separately.
    pub(crate) fn load_frame_in_flow(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        ambient_velocity_world: Vec3,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingLoad, Error> {
        self.load_frame_in_flow_with_lift(
            frame,
            state,
            ambient_velocity_world,
            SphereLift::None,
            cancelled,
        )
    }

    fn load_frame_with_lift(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        lift_model: SphereLift,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingLoad, Error> {
        self.load_frame_in_flow_with_lift(
            frame,
            state,
            self.ambient_velocity_world,
            lift_model,
            cancelled,
        )
    }

    fn load_frame_in_flow_with_lift(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        ambient_velocity_world: Vec3,
        lift_model: SphereLift,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingLoad, Error> {
        poll(&mut cancelled)?;
        vector(ambient_velocity_world)?;
        if lift_model == SphereLift::RubinowKeller {
            self.require_no_slip()?;
        }
        let properties = frame.properties();
        let relative = state
            .center_of_mass_velocity_world(properties)?
            .sub(ambient_velocity_world);
        let omega = properties.angular_velocity_body_checked(state.angular_momentum_body())?;
        let translation_reynolds = self
            .transport
            .reynolds(norm(relative)?, 2.0 * self.radius_m)?;
        let rotation_reynolds = finite(norm(omega)? * self.diffusion_time_s)?;
        if translation_reynolds > self.limits.maximum_reynolds
            || rotation_reynolds > self.limits.maximum_reynolds
        {
            return Err(Error::Reynolds {
                translation: translation_reynolds,
                rotation: rotation_reynolds,
                maximum: self.limits.maximum_reynolds,
            });
        }
        let drag = vector(relative.scale(-self.resistance_kg_s))?;
        let lift = match lift_model {
            SphereLift::None => Vec3::ZERO,
            SphereLift::RubinowKeller => vector(
                state
                    .pose()
                    .orientation()
                    .rotate_body_to_world(omega)
                    .cross(relative)
                    .scale(0.75 * self.fluid_mass_kg),
            )?,
        };
        let rotary = vector(omega.scale(-self.rotary_resistance_kg_m2_s))?;
        let buoyancy = vector(self.gravity.acceleration_world().scale(-self.fluid_mass_kg))?;
        let dissipated_power_w = finite(-drag.dot(relative) - rotary.dot(omega))?;
        let fluid_force = match lift_model {
            SphereLift::None => drag,
            SphereLift::RubinowKeller => vector(drag.add(lift))?,
        };
        let background_power_w = finite(fluid_force.dot(ambient_velocity_world))?;
        if dissipated_power_w < 0.0
            || (relative != Vec3::ZERO || omega != Vec3::ZERO) && dissipated_power_w == 0.0
        {
            return Err(Error::Invalid(
                "nonzero relative motion needs representable positive viscous dissipation",
            ));
        }
        let wrench = Wrench {
            force_world: vector(fluid_force.add(buoyancy))?,
            torque_body: rotary,
        };
        poll(&mut cancelled)?;
        Ok(CreepingLoad {
            wrench,
            drag_force_world_n: drag,
            lift_force_world_n: lift,
            rotary_torque_body_nm: rotary,
            buoyancy_force_world_n: buoyancy,
            dissipated_power_w,
            background_power_w,
            translation_reynolds,
            rotation_reynolds,
        })
    }

    /// Close the native midpoint load against its actual velocity. Algebraic
    /// elimination of linear drag supplies a wrench; the existing integrator
    /// advances the only physical state. Endpoint and midpoint Reynolds checks,
    /// a fluid/body timescale gate and an independent work gate precede return.
    /// A failed late check or cancellation leaves the input state unchanged.
    ///
    /// Uniform body/shape/fluid data remain fixed during this step. Fluid heat
    /// and maintained-flow work are external reservoir accounts. The finite
    /// body's material density and the displaced fluid density are distinct.
    pub fn step(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        controls: CreepingStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingStep, Error> {
        poll(&mut cancelled)?;
        self.admit_body(body)?;
        self.step_frame(body.frame(), state, controls, &mut cancelled)
    }

    /// Compose drag with an external COM force and body-frame free couple,
    /// held over this interval. This permits driven motion without applying
    /// drag and the drive in separate, inconsistent velocity updates.
    ///
    /// The force is expressed in world axes and the couple in the body's
    /// principal axes, as in the native [`Wrench`]. Gravity and fluid forces
    /// are already supplied by this owner and must not be included again.
    /// A force at an offset point or a world-held torque generally changes
    /// its body-frame moment as the sphere turns; that load needs its own
    /// stage evaluation and is not inferred from a held wrench.
    /// Work is evaluated at the same midpoint as drag and reported separately
    /// from fluid heat and maintained-flow work. All existing regime,
    /// timestep, energy and cancellation checks still precede publication.
    pub fn step_with_wrench(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        applied: Wrench,
        controls: CreepingStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingStep, Error> {
        poll(&mut cancelled)?;
        self.admit_body(body)?;
        self.step_frame_with_wrench(body.frame(), state, applied, controls, &mut cancelled)
    }

    /// Advance with [`Self::load_with_rotational_lift`] at the original native
    /// midpoint, optionally driven by a held external COM force/body couple.
    /// Supply [`Wrench::ZERO`] for unforced motion; gravity is already included.
    ///
    /// The transverse force is solved simultaneously with drag, rather than
    /// applied as an explicit kick that would spuriously increase kinetic energy.
    /// A closed-form three-dimensional skew solve adds no iterations or heap
    /// allocations. The original rigid-body integrator remains the sole motion
    /// owner, with the same Reynolds, diffusion, work and cancellation checks.
    pub fn step_with_rotational_lift(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        applied: Wrench,
        controls: CreepingStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingStep, Error> {
        poll(&mut cancelled)?;
        self.admit_body(body)?;
        self.step_frame_with_lift(
            body.frame(),
            state,
            applied,
            controls,
            SphereLift::RubinowKeller,
            cancelled,
        )
    }

    /// Native midpoint closure shared by frozen and thermally updated bodies.
    /// The caller has already admitted this frame's one uniform sphere.
    pub(crate) fn step_frame(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        controls: CreepingStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingStep, Error> {
        self.step_frame_with_wrench(frame, state, Wrench::ZERO, controls, cancelled)
    }

    fn step_frame_with_wrench(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        applied: Wrench,
        controls: CreepingStepControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingStep, Error> {
        self.step_frame_with_lift(
            frame,
            state,
            applied,
            controls,
            self.selected_lift,
            cancelled,
        )
    }

    fn step_frame_with_lift(
        &self,
        frame: &PrincipalFrame,
        state: RigidBodyState,
        applied: Wrench,
        controls: CreepingStepControls,
        lift_model: SphereLift,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingStep, Error> {
        poll(&mut cancelled)?;
        vector(applied.force_world)?;
        vector(applied.torque_body)?;
        if !small_positive(controls.maximum_relaxation_fraction, 1.0)
            || !controls.duration_s.is_finite()
            || controls.duration_s <= 0.0
            || !controls.maximum_energy_defect_j.is_finite()
            || controls.maximum_energy_defect_j <= 0.0
        {
            return Err(Error::Invalid(
                "positive finite step/work controls and relaxation fraction in (0,1] required",
            ));
        }
        let load_before = self.load_frame_with_lift(frame, state, lift_model, &mut cancelled)?;
        let p = frame.properties();
        let inertia = p.principal_inertia_body();
        let translation_time = positive(p.mass() / self.resistance_kg_s)?;
        let (shortest_time, diffusion_ratio) = self.timescales_frame(frame)?;
        let dt = controls.duration_s;
        if finite(dt / shortest_time)? > controls.maximum_relaxation_fraction {
            return Err(Error::Invalid(
                "step exceeds the declared body relaxation fraction; reduce dt",
            ));
        }
        let half_dt = positive(0.5 * dt)?;
        let velocity = state.center_of_mass_velocity_world(p)?;
        let omega = p.angular_velocity_body_checked(state.angular_momentum_body())?;
        let acceleration = vector(
            self.gravity
                .acceleration_world()
                .scale(1.0 - self.fluid_mass_kg / p.mass())
                .add(applied.force_world.scale(1.0 / p.mass())),
        )?;
        let mut relative_mid = vector(
            velocity
                .sub(self.ambient_velocity_world)
                .add(acceleration.scale(half_dt))
                .scale(1.0 / (1.0 + half_dt / translation_time)),
        )?;
        let angular_mid = vector(Vec3::new(
            (omega.x + half_dt * applied.torque_body.x / inertia.x)
                / (1.0 + half_dt * self.rotary_resistance_kg_m2_s / inertia.x),
            (omega.y + half_dt * applied.torque_body.y / inertia.y)
                / (1.0 + half_dt * self.rotary_resistance_kg_m2_s / inertia.y),
            (omega.z + half_dt * applied.torque_body.z / inertia.z)
                / (1.0 + half_dt * self.rotary_resistance_kg_m2_s / inertia.z),
        ))?;
        if lift_model == SphereLift::RubinowKeller && angular_mid != Vec3::ZERO {
            // (A I - [b]x) r_mid = rhs, A = 1+h*c/m,
            // b = h*pi*rho*a^3*omega_world/m, h = dt/2.
            // With beta=b/A and u=rhs/A, the inverse is
            // r_mid=(u + beta cross u + beta*(beta dot u))/(1+|beta|^2).
            // pi*rho*a^3/c = (a^2/nu)/6. The admitted Reynolds and
            // relaxation limits keep beta small; no speed division is needed.
            // The native attitude increment rotates about angular_mid, so
            // its world direction is identical at the start and midpoint.
            let omega_world = state.pose().orientation().rotate_body_to_world(angular_mid);
            let fraction = half_dt / translation_time;
            let beta = vector(
                omega_world.scale(fraction / (1.0 + fraction) * (self.diffusion_time_s / 6.0)),
            )?;
            relative_mid = vector(
                relative_mid
                    .add(beta.cross(relative_mid))
                    .add(beta.scale(beta.dot(relative_mid)))
                    .scale(1.0 / positive(1.0 + beta.dot(beta))?),
            )?;
        }
        let velocity_mid = vector(relative_mid.add(self.ambient_velocity_world))?;
        let midpoint = RigidBodyState::new(
            state.pose(),
            velocity_mid.scale(p.mass()),
            Vec3::new(
                angular_mid.x * inertia.x,
                angular_mid.y * inertia.y,
                angular_mid.z * inertia.z,
            ),
        )?;
        let load_midpoint =
            self.load_frame_with_lift(frame, midpoint, lift_model, &mut cancelled)?;
        let total = Wrench {
            force_world: vector(load_midpoint.wrench.force_world.add(applied.force_world))?,
            torque_body: vector(load_midpoint.wrench.torque_body.add(applied.torque_body))?,
        };
        let dynamics = RigidBodyIntegrator::new(self.gravity).step(state, p, total, dt)?;
        let load_after =
            self.load_frame_with_lift(frame, dynamics.state_after, lift_model, &mut cancelled)?;
        let viscous_heat_j = finite(dt * load_midpoint.dissipated_power_w)?;
        let background_work_j = finite(dt * load_midpoint.background_power_w)?;
        let applied_power_w =
            finite(applied.force_world.dot(velocity_mid) + applied.torque_body.dot(angular_mid))?;
        let applied_work_j = finite(dt * applied_power_w)?;
        if load_midpoint.dissipated_power_w > 0.0 && viscous_heat_j == 0.0
            || load_midpoint.background_power_w != 0.0 && background_work_j == 0.0
            || applied_power_w != 0.0 && applied_work_j == 0.0
        {
            return Err(Error::Invalid(
                "fluid work or heat vanished in conversion to energy",
            ));
        }
        let dx = dynamics
            .state_after
            .pose()
            .position_world()
            .sub(state.pose().position_world());
        let before = dynamics.diagnostics_before;
        let after = dynamics.diagnostics_after;
        // Evaluate a displacement, not a subtraction of large potential datums.
        let stored_energy_change_j = finite(
            (after.translational_kinetic_energy - before.translational_kinetic_energy)
                + (after.rotational_kinetic_energy - before.rotational_kinetic_energy)
                - (p.mass() - self.fluid_mass_kg) * self.gravity.acceleration_world().dot(dx),
        )?;
        let energy_defect_j =
            finite(stored_energy_change_j + viscous_heat_j - background_work_j - applied_work_j)?;
        if energy_defect_j.abs() > controls.maximum_energy_defect_j {
            return Err(Error::EnergyDefect {
                defect_j: energy_defect_j,
                maximum_j: controls.maximum_energy_defect_j,
            });
        }
        poll(&mut cancelled)?;
        Ok(CreepingStep {
            dynamics,
            load_before,
            load_midpoint,
            load_after,
            viscous_heat_j,
            background_work_j,
            applied_work_j,
            stored_energy_change_j,
            energy_defect_j,
            diffusion_ratio,
        })
    }

    /// Fastest body relaxation and the admitted fluid/body diffusion ratio.
    /// Re-evaluate after configuration changes even at a thermal endpoint that
    /// has no subsequent mechanical step in the current advance.
    pub(crate) fn timescales_frame(&self, frame: &PrincipalFrame) -> Result<(f64, f64), Error> {
        let p = frame.properties();
        let inertia = p.principal_inertia_body();
        let translation_time = positive(p.mass() / self.resistance_kg_s)?;
        let rotation_time =
            positive(inertia.x.min(inertia.y).min(inertia.z) / self.rotary_resistance_kg_m2_s)?;
        let shortest_time = translation_time.min(rotation_time);
        let diffusion_ratio = positive(self.diffusion_time_s / shortest_time)?;
        if diffusion_ratio > self.limits.maximum_diffusion_ratio {
            return Err(Error::Unsteady {
                diffusion_ratio,
                maximum: self.limits.maximum_diffusion_ratio,
            });
        }
        Ok((shortest_time, diffusion_ratio))
    }

    /// Frozen-volume buoyancy potential, using the native gravity's world
    /// origin as datum. With changing volume its parameter change is an
    /// explicit reservoir account, NOT an inferred thermodynamic p dV law.
    pub(crate) fn buoyancy_potential_j(&self, state: RigidBodyState) -> Result<f64, Error> {
        finite(
            self.fluid_mass_kg
                * self
                    .gravity
                    .acceleration_world()
                    .dot(state.pose().position_world()),
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SphereLift {
    None,
    RubinowKeller,
}

fn small_positive(x: f64, maximum: f64) -> bool {
    x.is_finite() && x > 0.0 && x <= maximum
}
fn finite(x: f64) -> Result<f64, Error> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(Error::Invalid("unrepresentable creeping-flow arithmetic"))
    }
}
fn positive(x: f64) -> Result<f64, Error> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err(Error::Invalid(
            "positive creeping-flow scale is unrepresentable",
        ))
    }
}
fn vector(v: Vec3) -> Result<Vec3, Error> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(Error::Invalid("unrepresentable creeping-flow vector"))
    }
}
fn norm(v: Vec3) -> Result<f64, Error> {
    finite(v.x.hypot(v.y).hypot(v.z))
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
