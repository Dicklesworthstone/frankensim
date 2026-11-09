//! Affine prescribed air in the existing three-mode strip kinetic model.
//!
//! For a local mode λ and ambient normal-flow coefficient β(q,t), let
//! y=λ·V−β and T=m*y²/2. A native body-column η receives
//! m[y*(λ·[V,η]−Dβ[η])−(λ·η)*(λ·c−∂tβ−Dβ[V])].
//! Together with the original −M_added*udot, this is the variational reaction.
//! Both advection and configuration-gradient forces are required; replacing
//! only uniform-flow acceleration would violate this model's work identity.

use super::{
    AeroGeneralizedLoads, FlowDiagnostics, GeneralizedKinematics, LinkedStrip, MAX_CONTROL_COORDS,
    Refusal, UniformFlow, UniformFlowAssembly, UniformFlowReaction, assemble,
    assemble_articulated_strips_in_uniform_flow, base_axes_world, finite_vector, native_column,
    refuse,
};
use fs_ga::{Mat3, Twist, Vec3};

/// Instantaneous coefficients of `U(x,t)=U0(t)+A(t)*(x−reference_position)`.
/// The reference point and all vectors/matrices use fixed inertial world axes.
/// Derivatives are partial time derivatives at fixed world position, not a
/// material derivative or a finite difference along the aircraft trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AffineFlow {
    /// Fixed world reference point [m].
    pub reference_position_world_m: Vec3,
    /// Velocity U0 at that reference point [m/s].
    pub velocity_world_mps: Vec3,
    /// Partial time derivative of U0 [m/s²].
    pub partial_time_velocity_world_mps2: Vec3,
    /// A_ij=∂U_i/∂x_j [1/s], including shear, rotation and strain.
    pub gradient_world_per_s: Mat3,
    /// Partial time derivative of A [1/s²].
    pub partial_time_gradient_world_per_s2: Mat3,
}

impl Default for AffineFlow {
    fn default() -> Self {
        Self {
            reference_position_world_m: Vec3::default(),
            velocity_world_mps: Vec3::default(),
            partial_time_velocity_world_mps2: Vec3::default(),
            gradient_world_per_s: Mat3::zero(),
            partial_time_gradient_world_per_s2: Mat3::zero(),
        }
    }
}

fn invalid(message: &str) -> Refusal {
    refuse(
        "affine-flow-invalid",
        message.into(),
        "supply finite coefficients and native motion at one physical instant",
    )
}

impl AffineFlow {
    fn admit(&self) -> Result<(), Refusal> {
        if !finite_vector(self.reference_position_world_m)
            || !finite_vector(self.velocity_world_mps)
            || !finite_vector(self.partial_time_velocity_world_mps2)
            || !self
                .gradient_world_per_s
                .m
                .iter()
                .chain(&self.partial_time_gradient_world_per_s2.m)
                .all(|x| x.is_finite())
        {
            return Err(invalid("nonfinite affine velocity or derivative"));
        }
        Ok(())
    }

    fn at(&self, point: Vec3) -> (Vec3, Vec3) {
        let offset = point - self.reference_position_world_m;
        (
            self.velocity_world_mps + self.gradient_world_per_s.apply(offset),
            self.partial_time_velocity_world_mps2
                + self.partial_time_gradient_world_per_s2.apply(offset),
        )
    }

    /// Fluid material acceleration `∂tU+(U·∇)U` at a fixed world point [m/s²].
    /// The strip law instead differentiates ambient flow along each moving
    /// strip, using its actual point velocity; those derivatives are distinct.
    pub fn material_acceleration_world_at(&self, point: Vec3) -> Result<Vec3, Refusal> {
        self.admit()?;
        if !finite_vector(point) {
            return Err(invalid("nonfinite flow sample point"));
        }
        let (velocity, partial_time) = self.at(point);
        let acceleration = partial_time + self.gradient_world_per_s.apply(velocity);
        if !finite_vector(velocity) || !finite_vector(acceleration) {
            return Err(invalid("unrepresentable fluid material acceleration"));
        }
        Ok(acceleration)
    }

    pub(super) fn prepare(
        &self,
        strip: &LinkedStrip,
        kinematics: &GeneralizedKinematics,
    ) -> Result<[ModeFlow; 3], Refusal> {
        let pose = kinematics.world_from_link(strip.link);
        let rotation = pose.rotation();
        let rotate = |v| {
            rotation
                .rotate(v)
                .map_err(|_| invalid("unrepresentable flow frame"))
        };
        let center = pose
            .transform_point(strip.center)
            .map_err(|_| invalid("unrepresentable strip center"))?;
        let normal = rotate(strip.normal)?;
        let span = rotate(strip.span)?;
        let chord = rotate(strip.span.cross(strip.normal))?;
        let (velocity, partial_time) = self.at(center);
        let a = self.gradient_world_per_s;
        let da = self.partial_time_gradient_world_per_s2;
        let ac = a.apply(chord);
        let aspan = a.apply(span);
        // Normal velocity over the finite plate is affine in chord/span.
        // Its three coefficients share the original translation, pitch and
        // spanwise rotary modes. Thus gradients across the strip are retained.
        let beta = [normal.dot(velocity), -normal.dot(ac), normal.dot(aspan)];
        let time = [
            normal.dot(partial_time),
            -normal.dot(da.apply(chord)),
            normal.dot(da.apply(span)),
        ];
        let derivative = |eta: Twist| -> Result<[f64; 3], Refusal> {
            let omega = rotate(eta.angular)?;
            let point_velocity = rotate(eta.linear + eta.angular.cross(strip.center))?;
            let dn = omega.cross(normal);
            let dc = omega.cross(chord);
            let ds = omega.cross(span);
            Ok([
                dn.dot(velocity) + normal.dot(a.apply(point_velocity)),
                -(dn.dot(ac) + normal.dot(a.apply(dc))),
                dn.dot(aspan) + normal.dot(a.apply(ds)),
            ])
        };
        let along_body = derivative(kinematics.body_twist(strip.link))?;
        let mut modes = std::array::from_fn(|i| ModeFlow {
            speed: beta[i],
            time_rate: time[i],
            body_rate: along_body[i],
            variation: [0.0; 6 + MAX_CONTROL_COORDS],
        });
        for column in 0..kinematics.coordinate_count() {
            let value = derivative(kinematics.body_jacobian(strip.link)[native_column(column)])?;
            for i in 0..3 {
                modes[i].variation[column] = value[i];
            }
        }
        if modes.iter().any(|mode| {
            ![mode.speed, mode.time_rate, mode.body_rate]
                .iter()
                .chain(&mode.variation)
                .all(|x| x.is_finite())
        }) {
            return Err(invalid("unrepresentable strip flow or gradient transport"));
        }
        Ok(modes)
    }
}

#[derive(Clone, Copy)]
pub(super) struct ModeFlow {
    pub(super) speed: f64,
    pub(super) time_rate: f64,
    pub(super) body_rate: f64,
    pub(super) variation: [f64; 6 + MAX_CONTROL_COORDS],
}

#[derive(Clone, Debug, Default, PartialEq)]
pub(super) struct AffinePower {
    bias: f64,
    acceleration: [f64; 6 + MAX_CONTROL_COORDS],
    scale: f64,
}

impl AffinePower {
    pub(super) fn add(
        &mut self,
        ambient: &ModeFlow,
        mass: f64,
        speed: f64,
        bias_speed: f64,
        row: &[f64],
    ) {
        let storage_exchange = -mass * (ambient.speed * bias_speed);
        let gradient_exchange = -mass * (speed * ambient.body_rate);
        self.bias += storage_exchange + gradient_exchange;
        self.scale += storage_exchange.abs() + gradient_exchange.abs();
        for (coefficient, &j) in self.acceleration.iter_mut().zip(row) {
            *coefficient -= mass * (ambient.speed * j);
        }
    }

    pub(super) fn evaluate(&self, acceleration: &[f64], scale: &mut f64) -> f64 {
        let mut power = self.bias;
        *scale += self.scale;
        for (&coefficient, &a) in self.acceleration.iter().zip(acceleration) {
            let term = coefficient * a;
            power += term;
            *scale += term.abs();
        }
        power
    }
}

/// The original coupled inertia blocks and affine-flow work at one native state.
/// This extends the declared finite-strip kinetic model; arbitrary shear or
/// vorticity is NOT a solved nonuniform potential-flow boundary-value problem.
/// Pressure-gradient/displaced-volume forces remain separate external wrenches.
#[derive(Clone, Debug, PartialEq)]
pub struct AffineFlowAssembly {
    snapshot: UniformFlowAssembly,
    power: Option<AffinePower>,
}

impl AffineFlowAssembly {
    /// Native simultaneous body/joint mass blocks and complete fluid bias.
    #[must_use]
    pub fn loads(&self) -> &AeroGeneralizedLoads {
        self.snapshot.loads()
    }

    /// Supply this stage's original blocks to the articulated midpoint solver.
    #[must_use]
    pub fn into_loads(self) -> AeroGeneralizedLoads {
        self.snapshot.into_loads()
    }

    /// Disturbance kinetic energy of the original three relative strip modes [J].
    #[must_use]
    pub fn disturbance_energy_j(&self) -> f64 {
        self.snapshot.disturbance_energy_j()
    }

    /// Affine generalized fluid momentum in Flyer order.
    #[must_use]
    pub fn generalized_fluid_momentum(&self) -> &[f64] {
        self.snapshot.generalized_fluid_momentum()
    }

    /// Native reaction and smooth work balance at the supplied acceleration.
    /// The shared receipt's ambient power includes the spatial-gradient port,
    /// `−Σ m[β*dy/dt + y*Dβ[V]]`; it is not only `U0 dot net_force`.
    pub fn reaction(&self, acceleration: &[f64]) -> Result<UniformFlowReaction, Refusal> {
        self.snapshot
            .reaction_with_ambient(acceleration, self.power.as_ref())
    }
}

/// Assemble linked strips in an affine world flow, including material-point
/// advection, moving strip bases, and conjugate configuration-gradient forces.
/// The same positive mode masses and native joint map supply all mass blocks;
/// no numerical derivative or second inertia solve is introduced. Zero spatial
/// gradient AND zero gradient time derivative delegate exactly to uniform flow.
///
/// Density, strip geometry and the affine coefficients describe one instant.
/// A smooth callback must supply them at its actual evaluation time. This model
/// excludes density derivatives, cross-surface fluid interaction, viscosity,
/// displaced-volume pressure forces, and a full solution for nonuniform flow.
pub fn assemble_articulated_strips_in_affine_flow(
    rho_kg_m3: f64,
    strips: &[LinkedStrip],
    kinematics: &GeneralizedKinematics,
    flow: AffineFlow,
) -> Result<AffineFlowAssembly, Refusal> {
    flow.admit()?;
    if flow.gradient_world_per_s == Mat3::zero()
        && flow.partial_time_gradient_world_per_s2 == Mat3::zero()
    {
        return Ok(AffineFlowAssembly {
            snapshot: assemble_articulated_strips_in_uniform_flow(
                rho_kg_m3,
                strips,
                kinematics,
                UniformFlow {
                    velocity_world_mps: flow.velocity_world_mps,
                    acceleration_world_mps2: flow.partial_time_velocity_world_mps2,
                },
            )?,
            power: None,
        });
    }
    let mut diagnostics = FlowDiagnostics::default();
    let mut power = AffinePower::default();
    let loads = assemble(
        rho_kg_m3,
        strips,
        kinematics,
        None,
        Some(&mut diagnostics),
        Some((&flow, &mut power)),
    )?;
    if ![
        diagnostics.energy,
        diagnostics.storage_bias_rate,
        diagnostics.bias_body_power,
        diagnostics.work_scale,
        power.bias,
        power.scale,
    ]
    .iter()
    .chain(&diagnostics.momentum)
    .chain(&diagnostics.homogeneous_momentum)
    .chain(&power.acceleration)
    .all(|v| v.is_finite())
    {
        return Err(invalid("unrepresentable affine-flow momentum or work"));
    }
    Ok(AffineFlowAssembly {
        snapshot: UniformFlowAssembly {
            loads,
            flow: UniformFlow::default(),
            base_axes_world: base_axes_world(kinematics)?,
            diagnostics,
        },
        power: Some(power),
    })
}
