//! Finite guided translation driven by the original nonlinear ALE fluid.
//!
//! A native flow candidate solves the moving-domain BDM1/P0 equations for each
//! trial endpoint speed. Its actual pressure/viscous boundary reaction closes
//! the body's backward-Euler momentum equation. A bounded bracketed scalar
//! solve replaces prescribed mesh motion; no drag or added-mass law is fitted.
//! Body momentum, world pose and the original fluid history publish together.
//!
//! The entire reference domain translates with the body, without changing its
//! metric or topology. Exterior callbacks therefore describe conditions on that
//! moving outer boundary, at endpoint WORLD coordinates. They do not represent
//! fixed obstacles left behind by the moving mesh. Rotation, deformation,
//! variable-density flow and general moving-mesh remapping remain separate.
//! The guide supplies transverse forces and a couple explicitly. The reported
//! body work identity does not establish whole-fluid passivity or viscous heat.

/// Finite translation with both momentum components free and orientation restrained.
pub mod planar;

use fs_flux::ale::translating::{
    TranslatingAleError, TranslatingAleFlow, TranslatingAleState, TranslatingAleStep,
};
use fs_flux::{
    FluxSystem, MaterialBoundaryLoad, MaterialFluxNumerics, MaterialNavierStokesBoundary,
    NavierStokesControls,
};
use fs_material::fluid_transport::FluidTransport;
use fs_math::dd::Dd;
use fs_mbd::{DynamicsError, MassProperties, MomentumImpulseReceipt, Pose, RigidBodyState, Vec3};
use std::cell::Cell;

/// Physical whole-body coefficients, including the original 2-D mesh depth.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GuidedTranslationConfig {
    /// Unit translation direction in the mesh's physical x/y frame.
    pub direction: [f64; 2],
    /// Extrusion depth [m], converting native forces per metre into N.
    pub depth_m: f64,
    /// Axial spring [N/m], at rest at the construction pose.
    pub stiffness_n_m: f64,
    /// Mechanical dashpot [N s/m], separate from resolved fluid viscosity.
    pub damping_n_s_m: f64,
    /// Allowed mismatch of body, mesh and represented drift velocities [m/s].
    /// Also bounds unrelaxed native affine normal-velocity trace corrections;
    /// the fluid's independent nonlinear equation residual gates still apply.
    pub velocity_tolerance_m_s: f64,
}

/// Explicit endpoint search and balance controls. The caller chooses a speed
/// interval containing the desired root; non-bracketing data refuse. This is
/// not an assertion that a nonlinear fluid/body problem has a unique root.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GuidedTranslationControls {
    /// Ascending finite endpoint speeds [m/s], including reverse motion.
    pub speed_bracket_m_s: [f64; 2],
    /// Maximum native nonlinear flow solves, including discarded candidates.
    pub maximum_flow_solves: usize,
    pub force_absolute_n: f64,
    pub force_relative: f64,
    pub work_absolute_j: f64,
    pub work_relative: f64,
}

#[derive(Debug)]
pub enum GuidedTranslationError {
    Invalid(&'static str),
    Flow(TranslatingAleError),
    Mechanics(DynamicsError),
    /// The actual residuals at both supplied interval endpoints have one sign.
    Unbracketed {
        residual_n: [f64; 2],
    },
    /// No endpoint passed the momentum tolerance within the declared work.
    SolveLimit {
        flow_solves: usize,
        residual_n: f64,
    },
    /// A momentum-converged endpoint failed its independent mechanical work gate.
    WorkBalance {
        residual_j: f64,
    },
    Cancelled,
}

impl core::fmt::Display for GuidedTranslationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "guided translating fluid/body coupling: {self:?}")
    }
}
impl std::error::Error for GuidedTranslationError {}
impl From<TranslatingAleError> for GuidedTranslationError {
    fn from(value: TranslatingAleError) -> Self {
        Self::Flow(value)
    }
}
impl From<DynamicsError> for GuidedTranslationError {
    fn from(value: DynamicsError) -> Self {
        Self::Mechanics(value)
    }
}

/// One accepted physical endpoint. The native flow report retains its actual
/// world momentum, kinetic energy, nonlinear residuals and moving-surface flux.
pub struct GuidedTranslationStep {
    pub flow: TranslatingAleStep,
    pub boundary_load: MaterialBoundaryLoad,
    /// Original native momentum jump, before the endpoint pose drift.
    pub momentum: MomentumImpulseReceipt,
    pub body_after: RigidBodyState,
    pub speed_m_s: f64,
    pub displacement_m: f64,
    pub flow_solves: usize,
    /// Full physical fluid force projected onto the guide [N].
    pub fluid_force_n: f64,
    /// Transverse guide force on the body [N], excluding axial spring/dashpot.
    pub guide_force_n: [f64; 2],
    /// Guide couple [N m] balancing the fluid moment about the moving COM.
    pub guide_torque_n_m: f64,
    pub interface_work_j: f64,
    pub external_work_j: f64,
    pub body_kinetic_change_j: f64,
    pub spring_energy_change_j: f64,
    pub damping_loss_j: f64,
    /// Backward-Euler body/spring loss; numerical, never assigned as heat.
    pub mechanical_time_discretization_loss_j: f64,
    pub force_residual_n: f64,
    pub mechanical_work_residual_j: f64,
}

/// Native mechanics on the SAME finite translating nonlinear fluid owner.
pub struct GuidedTranslationMotion<'s, 'm, 'a, F: FluidTransport + ?Sized> {
    flow: TranslatingAleFlow<'s, 'm, 'a, F>,
    properties: MassProperties,
    reference_pose: Pose,
    config: GuidedTranslationConfig,
    edges: Vec<usize>,
    body: RigidBodyState,
}

impl<'s, 'm, 'a, F: FluidTransport + ?Sized> GuidedTranslationMotion<'s, 'm, 'a, F> {
    /// Bind the original material, physical rigid inertia and mesh. Body edges
    /// are impermeable no-slip boundaries; all other edges are supplied at each
    /// step. Initial normal BDM moments must follow the original body velocity.
    /// The finite guide admits no initial transverse or angular momentum.
    /// Native `MassProperties` admits only a COM reference, so the original
    /// body's `Pose::position_world` also owns the moving load-moment origin.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        system: &'s FluxSystem<'m>,
        material: &'a F,
        numerics: MaterialFluxNumerics,
        properties: MassProperties,
        body: RigidBodyState,
        initial_world_velocity: &[f64],
        body_edges: &[usize],
        config: GuidedTranslationConfig,
    ) -> Result<Self, GuidedTranslationError> {
        let norm = config.direction[0].hypot(config.direction[1]);
        if !norm.is_finite()
            || (norm - 1.0).abs() > 64.0 * f64::EPSILON
            || !positive(config.depth_m)
            || !nonnegative(config.stiffness_n_m)
            || !nonnegative(config.damping_n_s_m)
            || !positive(config.velocity_tolerance_m_s)
            || initial_world_velocity.len() != system.n_u
            || !initial_world_velocity.iter().all(|x| x.is_finite())
            || body_edges.is_empty()
        {
            return Err(invalid(
                "finite physical coefficients, unit guide and complete initial fluid history required",
            ));
        }
        let mut edges = body_edges.to_vec();
        edges.sort_unstable();
        if edges.windows(2).any(|pair| pair[0] == pair[1])
            || edges
                .iter()
                .any(|&e| e >= system.mesh.edges.len() || !system.mesh.is_boundary(e))
        {
            return Err(invalid(
                "body interface requires unique original boundary edges",
            ));
        }
        let axis = axis(config.direction);
        let momentum = body.linear_momentum_world();
        let transverse = momentum.sub(axis.scale(momentum.dot(axis)));
        if body.angular_momentum_body() != Vec3::ZERO
            || !length(transverse).is_finite()
            || length(transverse) > 128.0 * f64::EPSILON * length(momentum).max(f64::MIN_POSITIVE)
        {
            return Err(invalid(
                "initial momentum must follow the translation guide without rotation",
            ));
        }
        let initial_velocity = body.center_of_mass_velocity_world(properties)?;
        let speed = initial_velocity.dot(axis);
        let off_guide_velocity = initial_velocity.sub(axis.scale(speed));
        if !length(off_guide_velocity).is_finite()
            || length(off_guide_velocity) > config.velocity_tolerance_m_s
        {
            return Err(invalid(
                "initial transverse body velocity exceeds the guide tolerance",
            ));
        }
        let velocity = config.direction.map(|d| d * speed);
        for &edge in &edges {
            let normal = system.mesh.edges[edge].normal;
            let expected = velocity[0] * normal[0] + velocity[1] * normal[1];
            if (initial_world_velocity[2 * edge] - expected).abs() > config.velocity_tolerance_m_s
                || initial_world_velocity[2 * edge + 1].abs() > config.velocity_tolerance_m_s
            {
                return Err(invalid("initial body and fluid normal velocities disagree"));
            }
        }
        let flow =
            TranslatingAleFlow::new(system, material, numerics, velocity, initial_world_velocity)?;
        Ok(Self {
            flow,
            properties,
            reference_pose: body.pose(),
            config,
            edges,
            body,
        })
    }

    #[must_use]
    pub const fn body(&self) -> RigidBodyState {
        self.body
    }
    /// Original accepted relative-flow history, actual world translation and clock.
    #[must_use]
    pub const fn flow_state(&self) -> &TranslatingAleState {
        self.flow.state()
    }
    #[must_use]
    pub const fn material(&self) -> &F {
        self.flow.material()
    }

    /// Close one finite-translation nonlinear fluid/body endpoint. External
    /// force acts along the guide through the COM; include any pressure on the
    /// body's unmeshed side explicitly. `exterior` owns every other boundary.
    /// Source force [N/m³], world velocity [m/s] and traction [Pa] callbacks
    /// must remain pure across all trials; each sees the trial endpoint WORLD
    /// position. Pressure offsets retain the native component-local meaning.
    ///
    /// Every trial uses the SAME accepted fluid and body histories. Safeguarded
    /// secants with periodic bisection close the actual represented momentum,
    /// and the native fluid retains its own independent nonlinear acceptance.
    /// Failure, cancellation or an exhausted solve budget publishes neither
    /// body nor fluid. The explicit speed bracket selects the root interval.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &mut self,
        timestep_s: f64,
        force_world: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialNavierStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: f64,
        fluid_controls: NavierStokesControls,
        controls: GuidedTranslationControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GuidedTranslationStep, GuidedTranslationError> {
        poll(&mut cancelled)?;
        let dt = timestep_s;
        let [lo, hi] = controls.speed_bracket_m_s;
        if !positive(dt)
            || !positive(dt.recip())
            || !external_force_n.is_finite()
            || !lo.is_finite()
            || !hi.is_finite()
            || !positive(hi - lo)
            || controls.maximum_flow_solves == 0
            || ![
                controls.force_absolute_n,
                controls.force_relative,
                controls.work_absolute_j,
                controls.work_relative,
            ]
            .into_iter()
            .all(positive)
        {
            return Err(invalid(
                "finite timestep, ascending speed bracket, nonzero work and positive tolerances required",
            ));
        }
        let config = self.config;
        let direction = config.direction;
        let axis = axis(direction);
        let body0 = self.body;
        let v0 = body0
            .center_of_mass_velocity_world(self.properties)?
            .dot(axis);
        let old_translation = self.flow.state().translation_m();
        let q0 = dot(old_translation, direction);
        let reference = self.reference_pose;
        let mass = self.properties.mass();
        // One boundary vector per step, reused across all nonlinear trials.
        // This callback is constant throughout every individual native solve.
        let wall_velocity = Cell::new([0.0; 2]);
        let wall = |_: [f64; 2]| wall_velocity.get();
        let mut boundaries = Vec::with_capacity(exterior.len() + self.edges.len());
        boundaries.extend_from_slice(exterior);
        boundaries.extend(
            self.edges
                .iter()
                .map(|&edge| MaterialNavierStokesBoundary::Velocity { edge, value: &wall }),
        );
        let mut bracket = Bracket::new(lo, hi);
        let mut last_residual = f64::NAN;
        for index in 0..controls.maximum_flow_solves {
            poll(&mut cancelled)?;
            let requested_speed = bracket.next(index)?;
            let momentum = body0.apply_momentum_impulse(
                self.properties,
                axis.scale(mass * (requested_speed - v0)),
                Vec3::ZERO,
            )?;
            let endpoint_velocity = momentum
                .state_after
                .center_of_mass_velocity_world(self.properties)?;
            let speed = endpoint_velocity.dot(axis);
            let off_guide_velocity = endpoint_velocity.sub(axis.scale(speed));
            if !length(off_guide_velocity).is_finite()
                || length(off_guide_velocity) > config.velocity_tolerance_m_s
            {
                return Err(invalid(
                    "native endpoint transverse velocity exceeds the guide tolerance",
                ));
            }
            let velocity = direction.map(|d| d * speed);
            wall_velocity.set(velocity);
            let candidate = self.flow.candidate(
                dt,
                velocity,
                config.velocity_tolerance_m_s,
                force_world,
                &boundaries,
                fluid_controls,
                config.velocity_tolerance_m_s,
                &mut cancelled,
            )?;
            let flow = candidate.step();
            let translation = flow.translation_m;
            let q1 = dot(translation, direction);
            let dx = q1 - q0;
            let body = RigidBodyState::new(
                Pose::new(
                    reference
                        .position_world()
                        .add(Vec3::new(translation[0], translation[1], 0.0)),
                    reference.orientation(),
                )?,
                momentum.state_after.linear_momentum_world(),
                Vec3::ZERO,
            )?;
            let actual_drift = body
                .pose()
                .position_world()
                .sub(body0.pose().position_world());
            if !q1.is_finite()
                || !dx.is_finite()
                || length(actual_drift.sub(axis.scale(dt * speed)))
                    > dt * config.velocity_tolerance_m_s
                || (dx / dt - speed).abs() > config.velocity_tolerance_m_s
                || (0..2).any(|i| {
                    (flow.mesh_velocity_m_s[i] - velocity[i]).abs() > config.velocity_tolerance_m_s
                })
            {
                return Err(invalid(
                    "body and mesh motion are not jointly representable at the velocity tolerance",
                ));
            }
            let center = body.pose().position_world();
            let load =
                flow.boundary_loads_world(&self.edges, [center.x, center.y], pressure_offsets_pa)?;
            // stress is the ORIGINAL discrete reaction, including strong normal
            // constraints and the full weak tangential SIP flux. Outward
            // advective momentum is never added to this impermeable wall load.
            let full_force = load.stress.force_n_per_m.map(|f| config.depth_m * f);
            let force = dot(full_force, direction);
            let inertia = momentum.linear_increment_world.dot(axis) / dt;
            let spring = config.stiffness_n_m * q1;
            let damping = config.damping_n_s_m * speed;
            let residual = sum([inertia, -force, -external_force_n, spring, damping]);
            let force_scale = maximum([inertia, force, external_force_n, spring, damping]);
            if ![force, inertia, spring, damping, residual, force_scale]
                .into_iter()
                .all(f64::is_finite)
            {
                return Err(invalid("unrepresentable physical momentum equation"));
            }
            last_residual = residual;
            if residual.abs() <= controls.force_absolute_n + controls.force_relative * force_scale {
                let interface_work = force * dx;
                let external_work = external_force_n * dx;
                let kinetic_change = momentum.work.kinetic_energy_change_j;
                let spring_change = 0.5 * config.stiffness_n_m * (q1 + q0) * dx;
                let damping_loss = damping * dt * speed;
                let impulse = momentum.linear_increment_world;
                let numerical_loss =
                    impulse_time_loss(impulse, mass) + 0.5 * config.stiffness_n_m * dx * dx;
                let work = sum([
                    kinetic_change,
                    spring_change,
                    damping_loss,
                    numerical_loss,
                    -interface_work,
                    -external_work,
                ]);
                let work_scale = maximum([
                    kinetic_change,
                    spring_change,
                    damping_loss,
                    numerical_loss,
                    interface_work,
                    external_work,
                ]);
                if ![
                    interface_work,
                    external_work,
                    kinetic_change,
                    spring_change,
                    damping_loss,
                    numerical_loss,
                    work,
                    work_scale,
                ]
                .into_iter()
                .all(f64::is_finite)
                {
                    return Err(invalid("unrepresentable mechanical work"));
                }
                if work.abs() > controls.work_absolute_j + controls.work_relative * work_scale {
                    return Err(GuidedTranslationError::WorkBalance { residual_j: work });
                }
                let guide_force = core::array::from_fn(|i| force * direction[i] - full_force[i]);
                let guide_torque = -config.depth_m * load.stress.moment_n;
                if !guide_force.into_iter().all(f64::is_finite) || !guide_torque.is_finite() {
                    return Err(invalid("unrepresentable guide reaction"));
                }
                poll(&mut cancelled)?;
                let flow = candidate.commit();
                self.body = body;
                return Ok(GuidedTranslationStep {
                    flow,
                    boundary_load: load,
                    momentum,
                    body_after: body,
                    speed_m_s: speed,
                    displacement_m: q1,
                    flow_solves: index + 1,
                    fluid_force_n: force,
                    guide_force_n: guide_force,
                    guide_torque_n_m: guide_torque,
                    interface_work_j: interface_work,
                    external_work_j: external_work,
                    body_kinetic_change_j: kinetic_change,
                    spring_energy_change_j: spring_change,
                    damping_loss_j: damping_loss,
                    mechanical_time_discretization_loss_j: numerical_loss,
                    force_residual_n: residual,
                    mechanical_work_residual_j: work,
                });
            }
            drop(candidate);
            bracket.observe(index, speed, residual)?;
        }
        Err(GuidedTranslationError::SolveLimit {
            flow_solves: controls.maximum_flow_solves,
            residual_n: last_residual,
        })
    }
}

// Two endpoint solves establish a true sign bracket. Secants use normalized
// residual weights, avoiding products of two possibly large forces. Every
// fourth trial bisects, so false-position endpoint stagnation cannot persist.
struct Bracket {
    speeds: [f64; 2],
    residuals: [f64; 2],
}
impl Bracket {
    fn new(lo: f64, hi: f64) -> Self {
        Self {
            speeds: [lo, hi],
            residuals: [0.0; 2],
        }
    }
    fn next(&self, index: usize) -> Result<f64, GuidedTranslationError> {
        if index < 2 {
            return Ok(self.speeds[index]);
        }
        let [lo, hi] = self.speeds;
        let scale = self.residuals[0].abs().max(self.residuals[1].abs());
        let a = self.residuals[0].abs() / scale;
        let b = self.residuals[1].abs() / scale;
        let weight = if index % 4 == 0 { 0.5 } else { a / (a + b) };
        let mut next = lo + (hi - lo) * weight;
        if next <= lo || next >= hi {
            next = lo + 0.5 * (hi - lo);
        }
        if !next.is_finite() || next <= lo || next >= hi {
            return Err(invalid(
                "momentum root interval has no representable interior speed",
            ));
        }
        Ok(next)
    }
    fn observe(
        &mut self,
        index: usize,
        speed: f64,
        residual: f64,
    ) -> Result<(), GuidedTranslationError> {
        if index < 2 {
            self.speeds[index] = speed;
            self.residuals[index] = residual;
            if index == 1 {
                if self.speeds[0] >= self.speeds[1] {
                    return Err(invalid(
                        "distinct bracket endpoints collapse in native body momentum",
                    ));
                }
                if self.residuals[0].is_sign_negative() == residual.is_sign_negative() {
                    return Err(GuidedTranslationError::Unbracketed {
                        residual_n: self.residuals,
                    });
                }
            }
        } else {
            if speed <= self.speeds[0] || speed >= self.speeds[1] {
                return Err(invalid(
                    "native body speed no longer lies inside the momentum bracket",
                ));
            }
            let side =
                usize::from(self.residuals[0].is_sign_negative() != residual.is_sign_negative());
            self.speeds[side] = speed;
            self.residuals[side] = residual;
        }
        Ok(())
    }
}

fn axis(d: [f64; 2]) -> Vec3 {
    Vec3::new(d[0], d[1], 0.0)
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn length(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
// ||Delta p||²/(2m), without squaring dimensional momentum before dividing by
// mass. Scaling each component also admits a finite energy when the unscaled
// three-component norm would overflow. This is always a nonnegative loss.
fn impulse_time_loss(impulse: Vec3, mass: f64) -> f64 {
    let scale = mass.sqrt() * core::f64::consts::SQRT_2;
    let norm = (impulse.x / scale)
        .hypot(impulse.y / scale)
        .hypot(impulse.z / scale);
    norm * norm
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn maximum<const N: usize>(values: [f64; N]) -> f64 {
    values.into_iter().map(f64::abs).fold(0.0, f64::max)
}
fn sum<const N: usize>(values: [f64; N]) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |s, x| s + Dd::from_f64(x))
        .to_f64()
}
fn invalid(message: &'static str) -> GuidedTranslationError {
    GuidedTranslationError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GuidedTranslationError> {
    if cancelled() {
        Err(GuidedTranslationError::Cancelled)
    } else {
        Ok(())
    }
}
