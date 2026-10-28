//! A material-defined rigid sphere in an unbounded Newtonian fluid.
//!
//! Stokes translation gives F=-6*pi*mu*a*(v-U); rotation gives
//! torque=-8*pi*mu*a^3*omega. Fully immersed buoyancy is -rho*V*g.
//! The shared fluid snapshot supplies rho and mu; the original material body
//! supplies mass and inertia. The existing rigid-body integrator owns ALL
//! momentum, position and attitude updates. Only its linear midpoint load
//! relation is closed here, not a second time integrator.
//!
//! References: R. Fitzpatrick, *Axisymmetric Stokes Flow Around a Solid Sphere*,
//! equations (10.94), (10.125)-(10.128),
//! <https://farside.ph.utexas.edu/teaching/336L/Fluid/node142.html>;
//! D. Tong, *Fluid Dynamics I and II: Questions on Viscous Flows*, rotating
//! sphere problem, <https://www.damtp.cam.ac.uk/user/tong/books/examples/f2.pdf>.
//!
//! This is a no-slip, incompressible, quasi-steady reduction. It refuses large
//! Reynolds numbers and body relaxation too fast for the caller's declared
//! fluid-diffusion allowance. No Basset history, added mass, slip/rarefaction,
//! free surface, walls, other particles or nonlinear drag is supplied. The
//! material's continuum/Newtonian validity remains an external model obligation.
//! Do not superpose another exterior drag law on this same sphere.

use crate::principal::PrincipalBody;
use crate::{MaterialInertiaError, VolumeMoments};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_mbd::{
    DynamicsError, Gravity, RigidBodyIntegrator, RigidBodyState, StepReceipt, Vec3, Wrench,
};

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
    /// Independent kinetic/potential/viscous-work balance failed.
    EnergyDefect {
        /// Signed mechanical energy plus heat minus external flow work [J].
        defect_j: f64,
        /// Caller-declared absolute tolerance [J].
        maximum_j: f64,
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
    resistance_kg_s: f64,
    rotary_resistance_kg_m2_s: f64,
    fluid_mass_kg: f64,
    diffusion_time_s: f64,
    gravity: Gravity,
    ambient_velocity_world: Vec3,
    limits: CreepingFlowLimits,
}

/// Instantaneous physical load at the actual supplied velocity, with body
/// gravity excluded. Dissipation goes to the prescribed fluid reservoir;
/// it is not automatically added to any finite material heat store.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CreepingLoad {
    /// Fluid drag, buoyancy and rotary resistance; excludes body weight.
    pub wrench: Wrench,
    /// Translational resistance in the world frame [N].
    pub drag_force_world_n: Vec3,
    /// Rotary resistance in the principal body frame [N m].
    pub rotary_torque_body_nm: Vec3,
    /// Fully immersed Archimedes force in the world frame [N].
    pub buoyancy_force_world_n: Vec3,
    /// Nonnegative viscous loss deposited in the external fluid reservoir [W].
    pub dissipated_power_w: f64,
    /// Work supplied by maintaining the uniform background flow, F_drag dot U.
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
    /// Viscous energy deposited in the external fluid reservoir [J].
    pub viscous_heat_j: f64,
    /// Energy supplied by maintaining the uniform background flow [J].
    pub background_work_j: f64,
    /// Delta(kinetic + body gravity + buoyancy potential) [J].
    pub stored_energy_change_j: f64,
    /// Signed stored-energy change plus heat minus background-flow work [J].
    pub energy_defect_j: f64,
    /// Fluid diffusion time divided by the fastest body relaxation time.
    pub diffusion_ratio: f64,
}

impl<'a, F: FluidTransport + ?Sized> CreepingSphere<'a, F> {
    /// Declare a fully immersed sphere in a spatially uniform, frozen fluid.
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
        let sphere = VolumeMoments::sphere(radius_m)?;
        if !ambient_velocity_world.is_finite()
            || !small_positive(limits.maximum_reynolds, 0.1)
            || !small_positive(limits.maximum_diffusion_ratio, 0.1)
        {
            return Err(Error::Invalid(
                "finite ambient velocity and limits in (0,0.1] required",
            ));
        }
        let transport = source.transport_state()?;
        let resistance_kg_s =
            positive(6.0 * core::f64::consts::PI * transport.dynamic_viscosity_pa_s() * radius_m)?;
        let rotary_resistance_kg_m2_s =
            positive((4.0 / 3.0) * resistance_kg_s * radius_m * radius_m)?;
        let fluid_mass_kg = positive(transport.density_kg_m3() * sphere.volume_m3())?;
        let diffusion_time_s =
            positive(radius_m * (radius_m / transport.kinematic_viscosity_m2_s()))?;
        Ok(Self {
            source,
            transport,
            sphere,
            radius_m,
            resistance_kg_s,
            rotary_resistance_kg_m2_s,
            fluid_mass_kg,
            diffusion_time_s,
            gravity,
            ambient_velocity_world,
            limits,
        })
    }

    /// Retain the original fluid owner and its provenance or mixture receipts.
    pub const fn source(&self) -> &'a F {
        self.source
    }
    /// Shared fluid snapshot resolved once when this reservoir was declared.
    pub const fn transport(&self) -> FluidTransportState {
        self.transport
    }
    /// Declared no-slip spherical boundary radius [m].
    pub const fn radius_m(&self) -> f64 {
        self.radius_m
    }
    /// The same gravity used for body weight and displaced-fluid buoyancy.
    pub const fn gravity(&self) -> Gravity {
        self.gravity
    }

    fn admit_body(&self, body: &PrincipalBody<'_>) -> Result<(), Error> {
        let regions = body.source().regions();
        if regions.len() != 1 {
            return Err(Error::Invalid(
                "creeping sphere requires one uniform spherical material region",
            ));
        }
        let geometry = regions[0].geometry();
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

    /// Evaluate Stokes drag, rotary resistance and buoyancy without advancing
    /// physical state. Original body weight is supplied ONCE by native gravity.
    pub fn load(
        &self,
        body: &PrincipalBody<'_>,
        state: RigidBodyState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CreepingLoad, Error> {
        poll(&mut cancelled)?;
        self.admit_body(body)?;
        let properties = body.properties();
        let relative = state
            .center_of_mass_velocity_world(properties)?
            .sub(self.ambient_velocity_world);
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
        let rotary = vector(omega.scale(-self.rotary_resistance_kg_m2_s))?;
        let buoyancy = vector(self.gravity.acceleration_world().scale(-self.fluid_mass_kg))?;
        let dissipated_power_w = finite(-drag.dot(relative) - rotary.dot(omega))?;
        let background_power_w = finite(drag.dot(self.ambient_velocity_world))?;
        if dissipated_power_w < 0.0
            || (relative != Vec3::ZERO || omega != Vec3::ZERO) && dissipated_power_w == 0.0
        {
            return Err(Error::Invalid(
                "nonzero relative motion needs representable positive viscous dissipation",
            ));
        }
        let wrench = Wrench {
            force_world: vector(drag.add(buoyancy))?,
            torque_body: rotary,
        };
        poll(&mut cancelled)?;
        Ok(CreepingLoad {
            wrench,
            drag_force_world_n: drag,
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
        let load_before = self.load(body, state, &mut cancelled)?;
        let p = body.properties();
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
                .scale(1.0 - self.fluid_mass_kg / p.mass()),
        )?;
        let relative_mid = vector(
            velocity
                .sub(self.ambient_velocity_world)
                .add(acceleration.scale(half_dt))
                .scale(1.0 / (1.0 + half_dt / translation_time)),
        )?;
        let velocity_mid = vector(relative_mid.add(self.ambient_velocity_world))?;
        let angular_mid = Vec3::new(
            omega.x / (1.0 + half_dt * self.rotary_resistance_kg_m2_s / inertia.x),
            omega.y / (1.0 + half_dt * self.rotary_resistance_kg_m2_s / inertia.y),
            omega.z / (1.0 + half_dt * self.rotary_resistance_kg_m2_s / inertia.z),
        );
        let midpoint = RigidBodyState::new(
            state.pose(),
            velocity_mid.scale(p.mass()),
            Vec3::new(
                angular_mid.x * inertia.x,
                angular_mid.y * inertia.y,
                angular_mid.z * inertia.z,
            ),
        )?;
        let load_midpoint = self.load(body, midpoint, &mut cancelled)?;
        let dynamics =
            RigidBodyIntegrator::new(self.gravity).step(state, p, load_midpoint.wrench, dt)?;
        let load_after = self.load(body, dynamics.state_after, &mut cancelled)?;
        let viscous_heat_j = finite(dt * load_midpoint.dissipated_power_w)?;
        let background_work_j = finite(dt * load_midpoint.background_power_w)?;
        if load_midpoint.dissipated_power_w > 0.0 && viscous_heat_j == 0.0
            || load_midpoint.background_power_w != 0.0 && background_work_j == 0.0
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
        let energy_defect_j = finite(stored_energy_change_j + viscous_heat_j - background_work_j)?;
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
            stored_energy_change_j,
            energy_defect_j,
            diffusion_ratio,
        })
    }
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
