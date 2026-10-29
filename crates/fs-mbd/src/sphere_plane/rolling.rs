//! Material-derived viscoelastic resistance in the original sphere trajectory.
//!
//! For the original first-order sphere law, `M = -A N R omega_t`. With
//! no contact-point slip, force and moment balance give
//! `(m + I/R²) du/dt = F_t - (A N/R) u`, where `u` is centre velocity
//! relative to the prescribed plane. This constant-load equation is solved
//! analytically, including turning and force-driven reversal. There is no
//! finite-time coast-down stop, velocity cutoff, or numerical loss term.
//!
//! The source's actual Hertz patch and its temperature, frequency, wave-speed
//! and relaxation limits are retained. Both endpoint traction norms bound the
//! entire exponential path, so infeasible no-slip motion refuses. A caller
//! must resolve sliding or impact through the corresponding contact owner;
//! this steady no-slip material law supplies neither mechanism. Hertz normal
//! storage is constant. The sphere uses the original rigid contact radius,
//! within the source's small-indentation approximation. Normal-axis spin is
//! retained and has no loss in this law. Heat is generated in the source's
//! deforming sphere; this module performs no thermal transport.

use super::*;
use fs_math::det;
use fs_tribo::InterfaceSystemRef;
use fs_tribo::rolling_loss::viscoelastic::{
    ViscoelasticSphereResponse, ViscoelasticSphereRollingCard,
};
use fs_tribo::rolling_loss::{
    RollingKinematics, RollingLossChannel, RollingLossError, RollingPatchReceipt,
    RollingWorkOwnership,
};

/// Original immutable constitutive and normal-contact sources for one interval.
#[derive(Clone, Copy, Debug)]
pub struct SphereRollingBinding<'a> {
    /// Material-derived Kelvin–Voigt sphere law; no fitted body damping.
    pub card: &'a ViscoelasticSphereRollingCard,
    /// Original Hertz normal patch at the supported load.
    pub patch: &'a RollingPatchReceipt,
    /// Ordered dry material interface and its source authority.
    pub interface: &'a InterfaceSystemRef,
    /// Dedicated rolling-deformation patch/interval/coordinate work owner.
    pub ownership: &'a RollingWorkOwnership,
    /// Held source temperature [K]; generated heat is reported separately.
    pub temperature_kelvin: f64,
    /// Held source excitation frequency [Hz], not inferred from body speed.
    pub excitation_frequency_hz: f64,
}

/// Accepted original-body motion and its separately identified bulk loss.
#[derive(Debug)]
pub struct ViscoelasticSphereAdvance<'a> {
    /// Native momentum state, trajectory work, impulses and balance residuals.
    pub motion: SpherePlaneAdvance,
    /// Exact source identities, work ownership and held material conditions.
    pub binding: SphereRollingBinding<'a>,
    /// Original constitutive response at the initial rolling-speed magnitude.
    pub response_before: ViscoelasticSphereResponse,
    /// Original constitutive response at the final rolling-speed magnitude.
    pub response_after: ViscoelasticSphereResponse,
    /// Actual maximum required static traction on the complete interval [N].
    pub maximum_required_traction_n: f64,
    /// Time integral of the material rolling torque about the centre [N m s].
    pub rolling_angular_impulse_world: Vec3,
    /// Nonnegative continuous torque work generated in the deforming sphere [J].
    /// This is already included in `motion.dissipated_energy`.
    pub deformation_heat_j: f64,
}

/// Atomic refusal; no projected state or partly accepted heat is returned.
#[derive(Debug)]
pub enum SphereRollingError {
    /// Original sphere dynamics or representability check refused.
    Dynamics(SpherePlaneError),
    /// Original material/normal-patch/applicability/work owner refused.
    Constitutive(RollingLossError),
    /// The requested interval needs more static traction than is available.
    TractionExceeded {
        /// Maximum endpoint demand, which bounds the whole exponential path.
        required_n: f64,
        /// Original Coulomb static capacity at this load.
        available_n: f64,
    },
}
impl fmt::Display for SphereRollingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "viscoelastic sphere rolling: {self:?}")
    }
}
impl std::error::Error for SphereRollingError {}
impl From<SpherePlaneError> for SphereRollingError {
    fn from(error: SpherePlaneError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<DynamicsError> for SphereRollingError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error.into())
    }
}
impl From<RollingLossError> for SphereRollingError {
    fn from(error: RollingLossError) -> Self {
        Self::Constitutive(error)
    }
}

/// Advance feasible supported rolling with the original material-derived loss.
///
/// The homogeneous sphere's radius, mass and inertia must agree with its source
/// density. No-slip and normal support must already be resolved. Static
/// traction must remain feasible throughout; an exhausted cone refuses instead
/// of imposing an impossible constraint. Translation, momenta, torque work and
/// force work are analytic. Changing-axis attitude uses the existing bounded
/// Lie-midpoint convention; commuting angular velocities integrate exactly.
///
/// Endpoint source evaluations bound speed applicability over this constant
/// coefficient interval. Constitutive/model-form approximation remains that
/// of the source card; this is not an arbitrary-speed rolling model or a
/// self-heating material update. Every refusal leaves caller state unchanged.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_viscoelastic_sphere_rolling<'a>(
    scene: &SpherePlane,
    state: SpherePlaneState,
    support_friction: &FrictionLaw,
    binding: SphereRollingBinding<'a>,
    duration_s: f64,
    config: SpherePlaneConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ViscoelasticSphereAdvance<'a>, SphereRollingError> {
    poll(&mut cancelled)?;
    if !duration_s.is_finite() || duration_s <= 0.0 {
        return Err(SpherePlaneError::Invalid("positive finite rolling duration required").into());
    }
    progress(state.time_s, duration_s)?;
    let end_time = finite(state.time_s + duration_s)?;
    let mut engine = Engine::new(*scene, config, support_friction)?;
    engine.velocity_scale = state.velocity_scale;
    let body = state.body;
    let x0 = body.pose().position_world();
    let v0 = body.center_of_mass_velocity_world(scene.properties)?;
    let spin0 = world_spin(body);
    let omega0 = spin0.scale(engine.inertia.recip());
    let normal_spin = engine.n.scale(omega0.dot(engine.n));
    let rolling0 = omega0.sub(normal_spin);
    let relative0 = v0.sub(scene.velocity_world);
    let gap = engine.gap(body, state.time_s)?;
    let geometry_allowance = engine.geometry_allowance(body, state.time_s)?;
    if gap < 0.0 && (!state.contact_known || -gap > geometry_allowance) {
        return Err(SpherePlaneError::Penetration { gap }.into());
    }
    if gap != 0.0 && (!state.contact_known || gap.abs() > geometry_allowance) {
        return Err(SpherePlaneError::Invalid("rolling requires resolved plane contact").into());
    }
    let (normal_speed, slip) = engine.kinematics(body)?;
    // Normal-axis spin does not participate in point slip. Its possibly large
    // kinetic energy/rate must not hide a resolved tangential sliding mode.
    let slip_allowance =
        4096.0 * f64::EPSILON * finite(norm(relative0)? + scene.radius * norm(rolling0)?)?;
    let normal_allowance = 4096.0 * f64::EPSILON * finite(norm(v0)? + norm(scene.velocity_world)?)?;
    if normal_speed.abs() > normal_allowance || slip[0].hypot(slip[1]) > slip_allowance {
        return Err(SpherePlaneError::Invalid("rolling requires resolved zero point slip").into());
    }
    let load = finite(-scene.force_world.dot(engine.n))?;
    if load < 0.0 {
        return Err(
            SpherePlaneError::Invalid("rolling cannot retain tensile normal support").into(),
        );
    }
    let parameters = binding.card.parameters();
    same(
        parameters.radius_m,
        scene.radius,
        "rolling source radius differs from body",
    )?;
    let volume = finite((4.0 * core::f64::consts::PI / 3.0) * scene.radius)?;
    let source_mass = finite(((volume * scene.radius) * scene.radius) * parameters.density_kg_m3)?;
    same(
        source_mass,
        engine.mass,
        "rolling source density differs from body mass",
    )?;
    same(
        finite(((0.4 * engine.mass) * scene.radius) * scene.radius)?,
        engine.inertia,
        "rolling source requires homogeneous sphere inertia",
    )?;
    same(
        binding.patch.normal_load_n(),
        load,
        "rolling patch differs from support load",
    )?;
    if binding.ownership.patch_id() != binding.patch.patch_id() {
        return Err(RollingLossError::WorkOwnershipMismatch { field: "patch_id" }.into());
    }
    if binding.ownership.channel() != RollingLossChannel::RollingHysteresis {
        return Err(RollingLossError::WorkOwnershipMismatch { field: "channel" }.into());
    }

    let radius = scene.radius;
    let effective_mass = finite(engine.mass + (engine.inertia / radius) / radius)?;
    let moment_coefficient = finite((binding.card.rolling_time_s() * load) * radius)?;
    let rate = finite(((moment_coefficient / radius) / radius) / effective_mass)?;
    if binding.card.rolling_time_s() > 0.0
        && load > 0.0
        && (moment_coefficient <= 0.0 || rate <= 0.0)
    {
        return Err(SpherePlaneError::Unresolved("unrepresentable rolling relaxation rate").into());
    }
    let z = finite(rate * duration_s)?;
    let moments = exponential_moments(z);
    let [decay, phi1, phi2, phi3, variance] = moments;
    let acceleration_free = engine.force_tangent_world.scale(effective_mass.recip());
    let initial_acceleration = acceleration_free.sub(relative0.scale(rate));
    if decay == 0.0 && initial_acceleration != Vec3::ZERO {
        return Err(
            SpherePlaneError::Unresolved("rolling exponential tail is unrepresentable").into(),
        );
    }
    let relative1 = relative0
        .scale(decay)
        .add(acceleration_free.scale(duration_s * phi1));
    let mean_relative = relative0
        .scale(phi1)
        .add(acceleration_free.scale(duration_s * phi2));
    let relative_displacement = mean_relative.scale(duration_s);
    let displacement = scene
        .velocity_world
        .scale(duration_s)
        .add(relative_displacement);
    let v1 = scene.velocity_world.add(relative1);
    for value in [
        relative1,
        mean_relative,
        displacement,
        v1,
        initial_acceleration,
    ] {
        vector(value)?;
    }
    if relative1 != Vec3::ZERO && v1.sub(scene.velocity_world) == Vec3::ZERO {
        return Err(
            SpherePlaneError::Unresolved("native rolling velocity is unrepresentable").into(),
        );
    }
    let traction0 = initial_acceleration
        .scale(engine.mass)
        .sub(engine.force_tangent_world);
    let traction1 = initial_acceleration
        .scale(engine.mass * decay)
        .sub(engine.force_tangent_world);
    let required_n = norm(traction0)?.max(norm(traction1)?);
    let available_n = finite(engine.support_static_mu * load)?;
    if required_n > available_n {
        return Err(SphereRollingError::TractionExceeded {
            required_n,
            available_n,
        });
    }

    let angular_drive = engine.n.cross(acceleration_free).scale(radius.recip());
    let rolling_at = |time: f64| {
        let [decay, phi1, ..] = exponential_moments(rate * time);
        rolling0.scale(decay).add(angular_drive.scale(time * phi1))
    };
    let rolling1 = rolling_at(duration_s);
    let mean_rolling = rolling0
        .scale(phi1)
        .add(angular_drive.scale(duration_s * phi2));
    let rolling_rotation = mean_rolling.scale(duration_s);
    let response = |rolling: Vec3| -> Result<ViscoelasticSphereResponse, SphereRollingError> {
        let angular_speed = norm(rolling)?;
        Ok(binding.card.response(
            binding.patch,
            binding.interface,
            RollingKinematics::new(
                finite(radius * angular_speed)?,
                angular_speed,
                omega0.dot(engine.n),
                binding.temperature_kelvin,
                binding.excitation_frequency_hz,
                duration_s,
            )?,
        )?)
    };
    let response_before = response(rolling0)?;
    let response_after = response(rolling1)?;
    // The nonnegative variance form integrates the actual squared rolling rate.
    // It stays positive during reversal and avoids subtracting large energy
    // or terminal-velocity terms to infer a small dissipated heat.
    let mean_square = finite(norm(mean_rolling)?.powi(2))?;
    let difference_square = finite(norm(rolling1.sub(rolling0))?.powi(2))?;
    let heat = finite(
        (moment_coefficient * duration_s) * finite(mean_square + variance * difference_square)?,
    )?;
    let rolling_impulse = rolling_rotation.scale(-moment_coefficient);
    vector(rolling_impulse)?;

    let mut orientation = body.pose().orientation();
    let angular_acceleration0 = angular_drive.sub(rolling0.scale(rate));
    if commuting_directions(omega0, angular_acceleration0) {
        let rotation = normal_spin.scale(duration_s).add(rolling_rotation);
        orientation = orientation.right_exp(orientation.rotate_world_to_body(rotation))?;
    } else {
        let count = (duration_s / config.max_attitude_step_s).ceil().max(1.0);
        if !count.is_finite() || count > config.max_attitude_substeps as f64 {
            return Err(SpherePlaneError::AttitudeBudget.into());
        }
        let steps = count as usize;
        let step = duration_s / steps as f64;
        for index in 0..steps {
            poll(&mut cancelled)?;
            let rotation = normal_spin
                .add(rolling_at((index as f64 + 0.5) * step))
                .scale(step);
            orientation = orientation.right_exp(orientation.rotate_world_to_body(rotation))?;
        }
    }
    let spin1 = normal_spin.add(rolling1).scale(engine.inertia);
    let after = RigidBodyState::new(
        Pose::new(x0.add(displacement), orientation)?,
        v1.scale(engine.mass),
        orientation.rotate_world_to_body(spin1),
    )?;
    let tangent_impulse = relative1
        .sub(relative0)
        .scale(engine.mass)
        .sub(engine.force_tangent_world.scale(duration_s));
    let contact_impulse = engine.n.scale(load * duration_s).add(tangent_impulse);
    // Independently integrate x(t) and x(t) cross m*a(t). The exponential
    // displacement parallel to the initial acceleration has zero cross product
    // with that acceleration, leaving these two analytic force-moment terms.
    let position_integral = x0
        .scale(duration_s)
        .add(
            scene
                .velocity_world
                .scale(0.5 * duration_s)
                .scale(duration_s),
        )
        .add(relative0.scale(phi2 * duration_s).scale(duration_s))
        .add(
            acceleration_free
                .scale(phi3 * duration_s)
                .scale(duration_s)
                .scale(duration_s),
        );
    let first_exponential_moment = if z < 0.25 {
        phi1 - phi2
    } else {
        ((1.0 - (1.0 + z) * decay) / z) / z
    };
    let acceleration_angular = x0
        .cross(initial_acceleration)
        .scale(engine.mass * duration_s * phi1)
        .add(
            v0.cross(initial_acceleration)
                .scale(engine.mass * first_exponential_moment)
                .scale(duration_s)
                .scale(duration_s),
        );
    let external_angular = position_integral.cross(scene.force_world);
    let contact_angular = acceleration_angular
        .add(position_integral.cross(engine.n.scale(load).sub(engine.force_tangent_world)))
        .add(engine.arm.cross(tangent_impulse))
        .add(rolling_impulse);
    let force_work = finite(scene.force_world.dot(displacement))?;
    let plane_work = finite(contact_impulse.dot(scene.velocity_world))?;
    let initial_energy = body.kinetic_energy(scene.properties)?;
    let final_energy = after.kinetic_energy(scene.properties)?;
    let energy_residual = finite(final_energy - initial_energy - force_work - plane_work + heat)?;
    let linear_residual = after
        .linear_momentum_world()
        .sub(body.linear_momentum_world())
        .sub(scene.force_world.scale(duration_s))
        .sub(contact_impulse);
    let angular_residual = total_angular(after)
        .sub(total_angular(body))
        .sub(external_angular)
        .sub(contact_angular);
    for value in [
        position_integral,
        contact_impulse,
        external_angular,
        contact_angular,
        linear_residual,
        angular_residual,
    ] {
        vector(value)?;
    }
    let tolerance = 16_384.0 * f64::EPSILON;
    let energy_scale = initial_energy + final_energy + force_work.abs() + plane_work.abs() + heat;
    let linear_scale = norm(body.linear_momentum_world())?
        + norm(after.linear_momentum_world())?
        + norm(scene.force_world)? * duration_s
        + norm(contact_impulse)?;
    let angular_scale = angular_magnitude(body)?
        + angular_magnitude(after)?
        + norm(external_angular)?
        + norm(contact_angular)?;
    if energy_residual.abs() > finite(tolerance * energy_scale)?
        || norm(linear_residual)? > finite(tolerance * linear_scale)?
        || norm(angular_residual)? > finite(tolerance * angular_scale)?
    {
        return Err(
            SpherePlaneError::Unresolved("rolling momentum or physical work balance").into(),
        );
    }
    let (normal_after, slip_after) = engine.kinematics(after)?;
    let slip_allowance =
        8192.0 * f64::EPSILON * finite(norm(relative1)? + radius * norm(rolling1)?)?;
    let normal_allowance = 8192.0 * f64::EPSILON * finite(norm(v1)? + norm(scene.velocity_world)?)?;
    if normal_after.abs() > normal_allowance || slip_after[0].hypot(slip_after[1]) > slip_allowance
    {
        return Err(SpherePlaneError::Unresolved("rolling endpoint point velocity").into());
    }
    let gap_after = engine.gap(after, end_time)?;
    if gap_after.abs() > engine.geometry_allowance(after, end_time)? {
        return Err(SpherePlaneError::Unresolved("rolling supported geometry").into());
    }
    if state.mode != SpherePlaneMode::Rolling {
        engine.event(state.time_s, SpherePlaneEventKind::Support, 0.0, Vec3::ZERO)?;
    }
    poll(&mut cancelled)?;
    Ok(ViscoelasticSphereAdvance {
        motion: SpherePlaneAdvance {
            state_after: SpherePlaneState {
                body: after,
                time_s: end_time,
                mode: SpherePlaneMode::Rolling,
                contact_known: true,
                zero_normal: true,
                velocity_scale: engine.velocity_scale,
            },
            events: engine.events,
            contact_impulse_world: contact_impulse,
            contact_angular_impulse_world: contact_angular,
            external_force_work: force_work,
            plane_work,
            dissipated_energy: heat,
            energy_balance_residual: energy_residual,
            linear_momentum_residual: linear_residual,
            angular_momentum_residual: angular_residual,
            gap_after,
            sliding_substeps: 0,
        },
        binding,
        response_before,
        response_after,
        maximum_required_traction_n: required_n,
        rolling_angular_impulse_world: rolling_impulse,
        deformation_heat_j: heat,
    })
}

fn same(actual: f64, expected: f64, field: &'static str) -> Result<(), SpherePlaneError> {
    finite(actual)?;
    finite(expected)?;
    if (actual - expected).abs() > 512.0 * f64::EPSILON * actual.abs().max(expected.abs()) {
        return Err(SpherePlaneError::Invalid(field));
    }
    Ok(())
}

// exp(-z), phi_1(-z), phi_2(-z), phi_3(-z), and variance of the
// endpoint-normalized exponential on [0,1]. Fixed Taylor work avoids loss of
// the force-driven and nearly elastic limits; no small damping is rounded away.
fn exponential_moments(z: f64) -> [f64; 5] {
    let decay = det::exp(-z);
    if z < 0.25 {
        let mut phi = [1.0, 0.5, 1.0 / 6.0];
        let mut term = phi;
        for order in 1..=16 {
            for index in 0..3 {
                term[index] *= -z / (order + index + 1) as f64;
                phi[index] += term[index];
            }
        }
        let z2 = z * z;
        let variance = 1.0 / 12.0
            + z2 * (-1.0 / 720.0
                + z2 * (1.0 / 30_240.0 + z2 * (-1.0 / 1_209_600.0 + z2 / 47_900_160.0)));
        [decay, phi[0], phi[1], phi[2], variance]
    } else {
        let one_minus_decay = -det::expm1(-z);
        let phi1 = one_minus_decay / z;
        let phi2 = (1.0 - phi1) / z;
        let phi3 = (0.5 - phi2) / z;
        let variance = (0.5 * (1.0 + decay) / one_minus_decay - 1.0 / z) / z;
        [decay, phi1, phi2, phi3, variance]
    }
}
