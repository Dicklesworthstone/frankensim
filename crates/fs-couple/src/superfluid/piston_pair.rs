//! Two moving impermeable ends on the original two-fluid and native body owners.
//!
//! Both pressures and midpoint velocities are solved together using the exact
//! two-end response of the channel's existing wave operators. Each body carries
//! its own native mass and a distinct endpoint half-cell of fluid; the latter
//! is absent from the channel's interior-face kinetic storage. External forces
//! on either component or either body retain their physical impulse and work.
//!
//! This is the same balanced-pressure, axial, small-displacement reference-grid
//! model as [`super::piston`]. With no external forces the combined longitudinal
//! momentum is conserved without a hidden fixed support. It does not introduce
//! moving meshes, nonlinear advection, rotational piston motion or an evolving
//! helium equation of state.

use fs_duct::superfluid::{
    BodyForce, BoundaryVelocities, ChannelEnd, ChannelStep, EndVelocities, HeatStepError,
    MutualFriction, TwoFluidChannel,
};
use fs_math::{dd::Dd, det};
use fs_mbd::{MassProperties, RigidBodyState, StepReceipt, Vec3};

use super::piston::{
    BalancedPistonBoundary, PistonControls, PistonError, admit_boundary, check, norm,
    prepare_motion, validate_controls,
};

/// One original native body's motion and its separate fluid half-cell account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonEndStep {
    /// Original native mechanical receipt; mass and inertia are not modified.
    pub mechanical: StepReceipt,
    /// Native outward velocity before the interval [m/s].
    pub velocity_before_outward_m_s: f64,
    /// Shared impermeable fluid/body midpoint velocity, positive outward [m/s].
    pub midpoint_velocity_outward_m_s: f64,
    /// Native outward velocity after the interval [m/s].
    pub velocity_after_outward_m_s: f64,
    /// Adjacent-cell pressure impulse on body plus half-cell, outward [N s].
    pub port_impulse_outward_n_s: f64,
    /// Actual native body momentum change, including its external load [N s].
    pub body_impulse_outward_n_s: f64,
    /// Fluid contact impulse on the native body, excluding its external load [N s].
    pub fluid_contact_impulse_outward_n_s: f64,
    /// Actual endpoint half-cell momentum increment, outward [N s].
    pub half_cell_impulse_outward_n_s: f64,
    /// Authored external impulse on the native body, outward [N s].
    pub external_body_impulse_outward_n_s: f64,
    /// Authored fluid body-force impulse on the endpoint half-cell, outward [N s].
    pub half_cell_source_impulse_outward_n_s: f64,
    /// Actual native body kinetic-energy increment [J].
    pub body_kinetic_change_j: f64,
    /// Actual half-cell kinetic-energy increment [J].
    pub half_cell_kinetic_change_j: f64,
    /// External work on the native body, at its actual midpoint velocity [J].
    pub external_body_work_j: f64,
    /// External body-force work on the fluid half-cell [J].
    pub half_cell_source_work_j: f64,
    /// Actual body plus half-cell increment minus pressure/source impulses [N s].
    pub interface_momentum_residual_world_n_s: Vec3,
}

/// Joint fluid and two-body result; endpoint order is always left, right.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonPairStep {
    /// Original wave/friction/source receipt, including actual interior momentum.
    pub channel: ChannelStep,
    /// Original native body updates and separate endpoint-fluid storage.
    pub ends: [PistonEndStep; 2],
    /// Fluid mass assigned to EACH endpoint [kg].
    pub boundary_half_mass_kg: f64,
    /// All authored fluid and body source impulses in the world frame [N s].
    pub external_impulse_world_n_s: Vec3,
    /// Interior fluid, endpoint-fluid and native body external work [J].
    pub total_external_work_j: f64,
    /// Actual total momentum increment minus all external impulses [N s].
    pub combined_momentum_residual_world_n_s: Vec3,
    /// Wave/body/half-cell storage change plus friction heat minus external work [J].
    pub combined_energy_residual_j: f64,
}

/// Advance both finite-mass pistons and the driven two-fluid channel atomically.
///
/// `boundaries` and `bodies` are ordered left, right, and share channel +x in
/// world coordinates. Each constant `external_force_outward_n` is an authored
/// force on that native body. `force` separately drives the original fluid
/// components, including the endpoint half-cell source account. Equal and
/// opposite endpoint component forces do no work because both velocities agree.
///
/// The response is anchored at the force-only midpoint predictor, so uniform
/// translation and common acceleration do not first create a fictitious blocked
/// pressure. A constant-size symmetric solve supplies the reciprocal pressure
/// correction. The actual final channel proposal still passes every original
/// physical, residual and cancellation gate. Native body positions and momenta
/// publish only after both endpoint and whole-system accounts have passed.
///
/// Friction heat is exported by the original channel and included once in this
/// account. Channel endpoint face samples are interval prescribed velocities;
/// the instantaneous half-cell momenta are defined by the native body states.
#[allow(clippy::too_many_arguments)]
pub fn advance_balanced_piston_pair(
    channel: &mut TwoFluidChannel,
    bodies: &mut [RigidBodyState; 2],
    properties: [MassProperties; 2],
    boundaries: [BalancedPistonBoundary; 2],
    dt_s: f64,
    force: BodyForce<'_>,
    external_force_outward_n: [f64; 2],
    friction: MutualFriction,
    controls: PistonControls,
    cancelled: impl Fn() -> bool,
) -> Result<PistonPairStep, PistonError> {
    advance_pair_with(
        channel,
        bodies,
        properties,
        boundaries,
        dt_s,
        force,
        external_force_outward_n,
        friction,
        controls,
        &cancelled,
        |_, _| Ok::<(), PistonError>(()),
    )
    .map(|(report, ())| report)
}

// The thermal composition uses this same solve and receives a private candidate
// inside the original channel's publication boundary, never a replayed step.
#[allow(clippy::too_many_arguments)]
fn advance_pair_with<T, E: From<PistonError>>(
    channel: &mut TwoFluidChannel,
    bodies: &mut [RigidBodyState; 2],
    properties: [MassProperties; 2],
    boundaries: [BalancedPistonBoundary; 2],
    dt_s: f64,
    force: BodyForce<'_>,
    external_force_outward_n: [f64; 2],
    friction: MutualFriction,
    controls: PistonControls,
    cancelled: &impl Fn() -> bool,
    prepare: impl FnOnce(&PistonPairStep, &[f64]) -> Result<T, E>,
) -> Result<(PistonPairStep, T), E> {
    check(cancelled)?;
    validate_controls(controls)?;
    if !dt_s.is_finite()
        || dt_s <= 0.0
        || external_force_outward_n.into_iter().any(|v| !v.is_finite())
        || boundaries[0].end != ChannelEnd::Left
        || boundaries[1].end != ChannelEnd::Right
    {
        return Err(PistonError::Invalid(
            "pair needs ordered ends, finite forces and positive duration",
        )
        .into());
    }
    let before = *bodies;
    let admitted = [
        admit_boundary(channel, before[0], properties[0], boundaries[0])?,
        admit_boundary(channel, before[1], properties[1], boundaries[1])?,
    ];
    let axis = admitted[0].0;
    if norm(axis.sub(admitted[1].0)) > 64.0 * f64::EPSILON {
        return Err(PistonError::Invalid(
            "both pistons must share the channel's positive world axis",
        )
        .into());
    }
    for i in 0..2 {
        let face = channel
            .face(i * channel.cells())
            .ok_or(PistonError::Invalid("missing piston endpoint"))?;
        let scale = face
            .normal_velocity_m_s
            .abs()
            .max(face.superfluid_velocity_m_s.abs());
        if face.counterflow_velocity_m_s.abs() > 64.0 * f64::EPSILON * scale {
            return Err(
                PistonError::Invalid("accepted piston endpoints must already have coflow").into(),
            );
        }
    }
    let signs = [-1.0, 1.0];
    let half_mass = channel.boundary_half_mass_kg();
    let effective_mass = properties.map(|p| p.mass() + half_mass);
    if !half_mass.is_finite()
        || half_mass <= 0.0
        || effective_mass
            .into_iter()
            .any(|m| !m.is_finite() || m <= 0.0)
    {
        return Err(PistonError::Invalid("endpoint effective masses are not representable").into());
    }
    let half_force = channel
        .boundary_body_force_outward_n(force, || cancelled())
        .map_err(PistonError::from)?;
    let u_before = [admitted[0].1.dot(axis.scale(-1.0)), admitted[1].1.dot(axis)];
    let external_impulse = external_force_outward_n.map(|f| f * dt_s);
    let half_source_impulse = half_force.map(|f| f * dt_s);
    let mut reference = [0.0; 2];
    for i in 0..2 {
        reference[i] = u_before[i]
            + 0.5 * (sum([external_impulse[i], half_source_impulse[i]]) / effective_mass[i]);
    }
    if !external_impulse
        .into_iter()
        .chain(half_source_impulse)
        .chain(reference)
        .all(f64::is_finite)
    {
        return Err(PistonError::Invalid(
            "endpoint source impulse or predictor is not representable",
        )
        .into());
    }
    let response = channel
        .impermeable_ends_response(dt_s, reference, force, friction, || cancelled())
        .map_err(PistonError::from)?;
    let correction = pressure_correction(
        effective_mass,
        0.5 * dt_s * channel.area_m2(),
        response.pressure_velocity_impedance_pa_s_m,
        response.reference_pressure_pa,
    )?;
    let midpoint = [reference[0] + correction[0], reference[1] + correction[1]];
    let ports = BoundaryVelocities {
        left: EndVelocities {
            normal_m_s: -midpoint[0],
            superfluid_m_s: -midpoint[0],
        },
        right: EndVelocities {
            normal_m_s: midpoint[1],
            superfluid_m_s: midpoint[1],
        },
    };
    let result = channel.step_forced_with_heat(
        dt_s,
        ports,
        force,
        friction,
        || cancelled(),
        |wave, heat| {
            let prepare_end = |i: usize| -> Result<PistonEndStep, PistonError> {
                check(cancelled)?;
                let outward = axis.scale(signs[i]);
                let port = -signs[i] * wave.boundary_pressure_impulse_n_s[i];
                let load = sum([port, external_impulse[i], half_source_impulse[i]]);
                let motion = prepare_motion(
                    before[i],
                    properties[i],
                    boundaries[i],
                    axis,
                    admitted[i].1,
                    midpoint[i],
                    (properties[i].mass() / effective_mass[i]) * load,
                    dt_s,
                    controls,
                )?;
                let actual_midpoint = motion.midpoint_velocity_outward_m_s;
                let body_impulse = motion.delta_momentum_world_n_s.dot(outward);
                let half_impulse = half_mass * (motion.velocity_after_outward_m_s - u_before[i]);
                let residual = motion.delta_momentum_world_n_s.add(outward.scale(sum([
                    half_impulse,
                    -port,
                    -external_impulse[i],
                    -half_source_impulse[i],
                ])));
                let end = PistonEndStep {
                    mechanical: motion.mechanical,
                    velocity_before_outward_m_s: u_before[i],
                    midpoint_velocity_outward_m_s: midpoint[i],
                    velocity_after_outward_m_s: motion.velocity_after_outward_m_s,
                    port_impulse_outward_n_s: port,
                    body_impulse_outward_n_s: body_impulse,
                    fluid_contact_impulse_outward_n_s: body_impulse - external_impulse[i],
                    half_cell_impulse_outward_n_s: half_impulse,
                    external_body_impulse_outward_n_s: external_impulse[i],
                    half_cell_source_impulse_outward_n_s: half_source_impulse[i],
                    body_kinetic_change_j: motion.kinetic_change_j,
                    half_cell_kinetic_change_j: half_impulse * actual_midpoint,
                    external_body_work_j: external_impulse[i] * actual_midpoint,
                    half_cell_source_work_j: half_source_impulse[i] * actual_midpoint,
                    interface_momentum_residual_world_n_s: residual,
                };
                if ![
                    load,
                    body_impulse,
                    half_impulse,
                    end.fluid_contact_impulse_outward_n_s,
                    end.body_kinetic_change_j,
                    end.half_cell_kinetic_change_j,
                    end.external_body_work_j,
                    end.half_cell_source_work_j,
                ]
                .into_iter()
                .all(f64::is_finite)
                    || !residual.is_finite()
                {
                    return Err(PistonError::Invalid(
                        "endpoint impulse or work is not representable",
                    ));
                }
                if norm(residual) > controls.momentum_tolerance_n_s {
                    return Err(PistonError::Invalid(
                        "native body and half-cell do not resolve pressure and source impulses",
                    ));
                }
                Ok(end)
            };
            let ends = [prepare_end(0)?, prepare_end(1)?];
            let sources = sum([
                wave.interior_body_impulse_n_s,
                -external_impulse[0],
                -half_source_impulse[0],
                external_impulse[1],
                half_source_impulse[1],
            ]);
            let mut momentum = axis.scale(sum([
                wave.resolved_momentum_after_kg_m_s,
                -wave.resolved_momentum_before_kg_m_s,
                -wave.interior_body_impulse_n_s,
            ]));
            let mut energy =
                Dd::from_f64(wave.wave_balance_defect_j) + Dd::from_f64(wave.boundary_wave_work_j);
            let mut external_work = Dd::from_f64(wave.body_force_work_j);
            for i in 0..2 {
                // Actual body increments are retained in world coordinates;
                // only the fluid's longitudinal momentum needs an axis lift.
                let dp = ends[i]
                    .mechanical
                    .state_after
                    .linear_momentum_world()
                    .sub(before[i].linear_momentum_world());
                momentum = momentum.add(dp).add(axis.scale(
                    signs[i]
                        * sum([
                            ends[i].half_cell_impulse_outward_n_s,
                            -external_impulse[i],
                            -half_source_impulse[i],
                        ]),
                ));
                energy = energy
                    + Dd::from_f64(ends[i].body_kinetic_change_j)
                    + Dd::from_f64(ends[i].half_cell_kinetic_change_j)
                    - Dd::from_f64(ends[i].external_body_work_j)
                    - Dd::from_f64(ends[i].half_cell_source_work_j);
                external_work = external_work
                    + Dd::from_f64(ends[i].external_body_work_j)
                    + Dd::from_f64(ends[i].half_cell_source_work_j);
            }
            let report = PistonPairStep {
                channel: *wave,
                ends,
                boundary_half_mass_kg: half_mass,
                external_impulse_world_n_s: axis.scale(sources),
                total_external_work_j: external_work.to_f64(),
                combined_momentum_residual_world_n_s: momentum,
                combined_energy_residual_j: energy.to_f64(),
            };
            if !report.external_impulse_world_n_s.is_finite()
                || !momentum.is_finite()
                || !report.total_external_work_j.is_finite()
                || !report.combined_energy_residual_j.is_finite()
            {
                return Err(PistonError::Invalid(
                    "combined momentum or energy is not representable",
                )
                .into());
            }
            if norm(momentum) > controls.momentum_tolerance_n_s {
                return Err(PistonError::Invalid(
                    "fluid and both native bodies do not conserve supplied momentum",
                )
                .into());
            }
            if report.combined_energy_residual_j.abs() > controls.energy_tolerance_j {
                return Err(PistonError::Invalid(
                    "fluid and both native bodies do not conserve supplied energy",
                )
                .into());
            }
            let candidate = prepare(&report, heat)?;
            Ok((report, candidate))
        },
    );
    let (report, candidate) = match result {
        Ok((_, accepted)) => accepted,
        Err(HeatStepError::Channel(error)) => return Err(PistonError::from(error).into()),
        Err(HeatStepError::Heat(error)) => return Err(error),
    };
    // Original channel acceptance has succeeded. No fallible work follows.
    *bodies = report.ends.map(|end| end.mechanical.state_after);
    Ok((report, candidate))
}

// Mass scaling makes the original 2x2 pressure response symmetric. Cholesky
// avoids subtracting nearly equal determinant products or inverting a matrix.
fn pressure_correction(
    masses: [f64; 2],
    pressure_impulse_scale: f64,
    resistance: [[f64; 2]; 2],
    pressure: [f64; 2],
) -> Result<[f64; 2], PistonError> {
    let roots = masses.map(det::sqrt);
    let diagonal = [
        1.0 + pressure_impulse_scale * (resistance[0][0] / masses[0]),
        1.0 + pressure_impulse_scale * (resistance[1][1] / masses[1]),
    ];
    let cross = (pressure_impulse_scale * (resistance[0][1] / roots[0])) / roots[1];
    let rhs = [
        pressure_impulse_scale * (pressure[0] / roots[0]),
        pressure_impulse_scale * (pressure[1] / roots[1]),
    ];
    let l00 = det::sqrt(diagonal[0]);
    let l10 = cross / l00;
    let schur = diagonal[1] - l10 * l10;
    if ![diagonal[0], diagonal[1], cross, rhs[0], rhs[1], schur]
        .into_iter()
        .all(f64::is_finite)
        || diagonal[0] <= 0.0
        || schur <= 0.0
    {
        return Err(PistonError::Invalid(
            "coupled endpoint pressure system is not resolvable",
        ));
    }
    let y0 = rhs[0] / l00;
    let x1 = (rhs[1] - l10 * y0) / schur;
    let x0 = (y0 - l10 * x1) / l00;
    let correction = [x0 / roots[0], x1 / roots[1]];
    if correction.into_iter().any(|v| !v.is_finite()) {
        return Err(PistonError::Invalid(
            "coupled endpoint velocity correction is not representable",
        ));
    }
    Ok(correction)
}

fn sum<const N: usize>(values: [f64; N]) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |total, value| total + Dd::from_f64(value))
        .to_f64()
}
