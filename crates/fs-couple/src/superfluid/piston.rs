//! Reciprocal small-displacement piston on the original two-fluid channel.
//!
//! The native rigid-body owner supplies mass, momentum and position. Both fluid
//! components follow its axial midpoint velocity at one impermeable end; the
//! opposite end is fixed. The endpoint half-cell fluid inertia is included
//! explicitly, because the channel's staggered wave energy excludes it.
//!
//! Forces are pressure PERTURBATIONS about an equilibrium whose reference
//! pressure is balanced by an exterior pressure or preload. This fixed-reference
//! linear model does not provide an unbalanced absolute-pressure force, a moving
//! mesh, rotational piston mechanics, or a nonlinear helium equation of state.

use fs_duct::superfluid::{
    BoundaryVelocities, ChannelEnd, ChannelError, ChannelStep, EndVelocities, HeatStepError,
    MutualFriction, TwoFluidChannel,
};
use fs_math::det;
use fs_mbd::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, StepReceipt, Vec3,
    Wrench,
};

/// Centered translating piston with explicitly balanced reference pressure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BalancedPistonBoundary {
    /// Moving end; the other end remains impermeable and fixed.
    pub end: ChannelEnd,
    /// Channel +x direction in the world frame, normalized by native mechanics.
    pub positive_axis_world: Vec3,
    /// Native center-of-mass position at the balanced reference configuration.
    pub equilibrium_position_world: Vec3,
    /// Maximum displacement from that configuration [m], at most 0.1 length.
    pub max_displacement_m: f64,
}

/// Absolute admission tolerances for the actual coupled increments.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonControls {
    /// Pressure impulse minus actual body plus half-cell momentum change [N s].
    pub momentum_tolerance_n_s: f64,
    /// Combined wave, body, half-cell and exported friction-heat residual [J].
    pub energy_tolerance_j: f64,
    /// Port/native midpoint velocity and displacement-per-duration mismatch [m/s].
    pub velocity_tolerance_m_s: f64,
}

/// Original owner reports and the physical endpoint-inertia correction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonStep {
    /// Original two-fluid step. Its friction heat remains available per cell.
    pub channel: ChannelStep,
    /// Original native body step, including its actual new position and momenta.
    pub mechanical: StepReceipt,
    /// Fluid mass between the adjacent cell center and moving end [kg].
    pub boundary_half_mass_kg: f64,
    /// Native outward end velocity before the step [m/s].
    pub velocity_before_outward_m_s: f64,
    /// Prescribed outward midpoint velocity of both components [m/s].
    pub midpoint_velocity_outward_m_s: f64,
    /// Native outward end velocity after the step [m/s].
    pub velocity_after_outward_m_s: f64,
    /// Adjacent-cell pressure impulse on body plus half-cell, outward [N s].
    pub port_impulse_outward_n_s: f64,
    /// Actual native body momentum increment projected outward [N s].
    pub body_impulse_outward_n_s: f64,
    /// Actual endpoint half-cell fluid momentum increment, outward [N s].
    pub half_cell_impulse_outward_n_s: f64,
    /// Fluid impulse ON the opposite stationary support, in channel +x [N s].
    pub stationary_support_impulse_x_n_s: f64,
    /// Body kinetic change from actual momentum increments and midpoint velocity [J].
    pub body_kinetic_change_j: f64,
    /// Endpoint half-cell kinetic-energy change [J].
    pub half_cell_kinetic_change_j: f64,
    /// Actual body plus half-cell impulse minus port impulse, world frame [N s].
    pub interface_momentum_residual_world_n_s: Vec3,
    /// Wave, body and half-cell change plus exported friction heat [J].
    pub combined_energy_residual_j: f64,
}

/// Any refusal preserves both accepted owners, including late cancellation.
#[derive(Debug)]
pub enum PistonError {
    /// Invalid geometry, unsupported motion, unresolved transfer or arithmetic.
    Invalid(&'static str),
    /// Original two-fluid refusal.
    Channel(ChannelError),
    /// Original native mechanical refusal.
    Mechanical(DynamicsError),
    /// Cancellation before joint publication.
    Cancelled,
}
impl core::fmt::Display for PistonError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "superfluid piston refused: {self:?}")
    }
}
impl std::error::Error for PistonError {}
impl From<ChannelError> for PistonError {
    fn from(error: ChannelError) -> Self {
        match error {
            ChannelError::Cancelled => Self::Cancelled,
            other => Self::Channel(other),
        }
    }
}
impl From<DynamicsError> for PistonError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanical(error)
    }
}
impl From<HeatStepError<PistonError>> for PistonError {
    fn from(error: HeatStepError<PistonError>) -> Self {
        match error {
            HeatStepError::Channel(error) => error.into(),
            HeatStepError::Heat(error) => error,
        }
    }
}

/// Advance the two-fluid waves and a native centered translating piston together.
///
/// The exact scalar endpoint response solves pressure and piston velocity
/// reciprocally. Native body inertia is unchanged; a separate rho*A*dx/2 fluid
/// mass receives its own momentum and kinetic energy. The native integrator
/// advances the body's position as well as its momentum. Both owners publish
/// only after native mechanics, displacement, transfer and cancellation checks.
///
/// `friction` dissipates counterflow only in the original channel. Its exported
/// joules are included once in this balance, but are not deposited in a thermal
/// owner here. They remain in `channel.last_mutual_friction_heat_j()`. The frozen
/// material reference is unchanged. Channel boundary mass describes compression
/// from swept volume on the reference grid, not leakage through this piston.
#[allow(clippy::too_many_arguments)]
pub fn advance_balanced_piston(
    channel: &mut TwoFluidChannel,
    body: &mut RigidBodyState,
    properties: MassProperties,
    boundary: BalancedPistonBoundary,
    dt_s: f64,
    friction: MutualFriction,
    controls: PistonControls,
    cancelled: impl Fn() -> bool,
) -> Result<PistonStep, PistonError> {
    check(&cancelled)?;
    if [
        controls.momentum_tolerance_n_s,
        controls.energy_tolerance_j,
        controls.velocity_tolerance_m_s,
    ]
    .into_iter()
    .any(|value| !value.is_finite() || value < 0.0)
    {
        return Err(PistonError::Invalid(
            "tolerances must be finite and nonnegative",
        ));
    }
    if !boundary.max_displacement_m.is_finite()
        || boundary.max_displacement_m <= 0.0
        || boundary.max_displacement_m > 0.1 * channel.length_m()
        || !boundary.equilibrium_position_world.is_finite()
    {
        return Err(PistonError::Invalid(
            "invalid small-displacement reference or limit",
        ));
    }
    let before = *body;
    let axis = before
        .directional_effective_mass(properties, Vec3::ZERO, boundary.positive_axis_world)?
        .unit_direction_world;
    let velocity_before = before.center_of_mass_velocity_world(properties)?;
    if before.angular_momentum_body() != Vec3::ZERO
        || norm(velocity_before.cross(axis)) > 64.0 * f64::EPSILON * norm(velocity_before)
    {
        return Err(PistonError::Invalid(
            "piston requires centered axial translation without spin",
        ));
    }
    admit_position(before, boundary, axis)?;
    let (end_index, sign) = match boundary.end {
        ChannelEnd::Left => (0, -1.0),
        ChannelEnd::Right => (1, 1.0),
    };
    let moving_face = channel
        .face(end_index * channel.cells())
        .ok_or(PistonError::Invalid("missing moving endpoint"))?;
    let fixed_face = channel
        .face((1 - end_index) * channel.cells())
        .ok_or(PistonError::Invalid("missing fixed endpoint"))?;
    let coflow_scale = moving_face
        .normal_velocity_m_s
        .abs()
        .max(moving_face.superfluid_velocity_m_s.abs());
    if fixed_face.normal_velocity_m_s != 0.0
        || fixed_face.superfluid_velocity_m_s != 0.0
        || moving_face.counterflow_velocity_m_s.abs() > 64.0 * f64::EPSILON * coflow_scale
    {
        return Err(PistonError::Invalid(
            "accepted opposite end must be closed and piston end must already have coflow",
        ));
    }
    let outward = axis.scale(sign);
    let u_before = velocity_before.dot(outward);
    let response =
        channel.impermeable_end_response(dt_s, boundary.end, friction, || cancelled())?;
    let half_mass = response.boundary_half_mass_kg;
    let effective_mass = properties.mass() + half_mass;
    if !effective_mass.is_finite()
        || !half_mass.is_finite()
        || half_mass <= 0.0
        || !response.blocked_pressure_pa.is_finite()
        || !response.pressure_velocity_impedance_pa_s_m.is_finite()
        || response.pressure_velocity_impedance_pa_s_m < 0.0
    {
        return Err(PistonError::Invalid(
            "endpoint mass or pressure response is not representable",
        ));
    }
    let admittance = (0.5 * dt_s) * (channel.area_m2() / effective_mass);
    let denominator = 1.0 + admittance * response.pressure_velocity_impedance_pa_s_m;
    let u_mid = (u_before + admittance * response.blocked_pressure_pa) / denominator;
    if !denominator.is_finite() || !u_mid.is_finite() {
        return Err(PistonError::Invalid(
            "coupled piston velocity is not representable",
        ));
    }
    let moving = EndVelocities {
        normal_m_s: sign * u_mid,
        superfluid_m_s: sign * u_mid,
    };
    let mut ports = BoundaryVelocities::CLOSED;
    match boundary.end {
        ChannelEnd::Left => ports.left = moving,
        ChannelEnd::Right => ports.right = moving,
    }
    let (_, report) = channel.step_with_heat(
        dt_s,
        ports,
        friction,
        || cancelled(),
        |wave, _heat| {
            check(&cancelled)?;
            // The channel reports pressure impulse INTO its interior, in +x.
            let port_impulse = -sign * wave.boundary_pressure_impulse_n_s[end_index];
            let body_impulse = port_impulse * (properties.mass() / effective_mass);
            let mechanical = RigidBodyIntegrator::new(Gravity::ZERO).step(
                before,
                properties,
                Wrench {
                    force_world: outward.scale(body_impulse / dt_s),
                    torque_body: Vec3::ZERO,
                },
                dt_s,
            )?;
            let after = mechanical.state_after;
            admit_position(after, boundary, axis)?;
            let velocity_after = after.center_of_mass_velocity_world(properties)?;
            let u_after = velocity_after.dot(outward);
            let average_velocity = velocity_before.scale(0.5).add(velocity_after.scale(0.5));
            let actual_u_mid = average_velocity.dot(outward);
            let displacement_velocity = after
                .pose()
                .position_world()
                .sub(before.pose().position_world())
                .scale(dt_s.recip());
            if !displacement_velocity.is_finite()
                || !average_velocity.is_finite()
                || (actual_u_mid - u_mid).abs() > controls.velocity_tolerance_m_s
                || norm(displacement_velocity.sub(average_velocity))
                    > controls.velocity_tolerance_m_s
            {
                return Err(PistonError::Invalid(
                    "native midpoint motion does not resolve the fluid port",
                ));
            }
            let delta_momentum = after
                .linear_momentum_world()
                .sub(before.linear_momentum_world());
            let half_impulse = half_mass * (u_after - u_before);
            let momentum_residual = delta_momentum.add(outward.scale(half_impulse - port_impulse));
            let body_energy_change = average_velocity.dot(delta_momentum);
            let half_energy_change = half_impulse * actual_u_mid;
            let combined_energy_residual = wave.wave_balance_defect_j
                + wave.boundary_wave_work_j
                + body_energy_change
                + half_energy_change;
            if ![
                port_impulse,
                half_impulse,
                body_energy_change,
                half_energy_change,
                combined_energy_residual,
                actual_u_mid,
                u_after,
            ]
            .into_iter()
            .all(f64::is_finite)
                || !momentum_residual.is_finite()
                || wave.body_force_work_j != 0.0
            {
                return Err(PistonError::Invalid(
                    "coupled impulse or energy is not representable",
                ));
            }
            if norm(momentum_residual) > controls.momentum_tolerance_n_s {
                return Err(PistonError::Invalid(
                    "body and half-cell do not resolve pressure impulse",
                ));
            }
            if combined_energy_residual.abs() > controls.energy_tolerance_j {
                return Err(PistonError::Invalid(
                    "coupled wave and piston energy does not close",
                ));
            }
            Ok(PistonStep {
                channel: *wave,
                mechanical,
                boundary_half_mass_kg: half_mass,
                velocity_before_outward_m_s: u_before,
                midpoint_velocity_outward_m_s: u_mid,
                velocity_after_outward_m_s: u_after,
                port_impulse_outward_n_s: port_impulse,
                body_impulse_outward_n_s: delta_momentum.dot(outward),
                half_cell_impulse_outward_n_s: half_impulse,
                stationary_support_impulse_x_n_s: -wave.boundary_pressure_impulse_n_s
                    [1 - end_index],
                body_kinetic_change_j: body_energy_change,
                half_cell_kinetic_change_j: half_energy_change,
                interface_momentum_residual_world_n_s: momentum_residual,
                combined_energy_residual_j: combined_energy_residual,
            })
        },
    )?;
    // No fallible operation remains after the native channel publishes.
    *body = report.mechanical.state_after;
    Ok(report)
}

fn admit_position(
    body: RigidBodyState,
    boundary: BalancedPistonBoundary,
    axis: Vec3,
) -> Result<(), PistonError> {
    let displacement = body
        .pose()
        .position_world()
        .sub(boundary.equilibrium_position_world);
    if !displacement.is_finite()
        || norm(displacement) > boundary.max_displacement_m
        || norm(displacement.cross(axis)) > 64.0 * f64::EPSILON * boundary.max_displacement_m
    {
        Err(PistonError::Invalid(
            "piston left its declared axial small-displacement domain",
        ))
    } else {
        Ok(())
    }
}

fn norm(value: Vec3) -> f64 {
    det::hypot(det::hypot(value.x, value.y), value.z)
}

fn check(cancelled: &impl Fn() -> bool) -> Result<(), PistonError> {
    if cancelled() {
        Err(PistonError::Cancelled)
    } else {
        Ok(())
    }
}
