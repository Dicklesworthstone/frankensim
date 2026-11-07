//! Physical held-wrench candidates; the parent alone applies native motion.

use super::*;
use fs_tribo::sliding_patch::contact::{HertzWrenchInput, HertzWrenchMode, HertzWrenchPhase};
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};

pub(super) struct Step {
    pub dt_s: f64,
    pub force_n: [f64; 2],
    pub spin_moment_nm: f64,
    pub velocity_after: [f64; 3],
    pub resting: bool,
    pub sticking: bool,
}

pub(super) enum Error {
    Mechanics(SpherePlaneError),
    Law(HertzWrenchError),
}
impl From<SpherePlaneError> for Error {
    fn from(error: SpherePlaneError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<HertzWrenchError> for Error {
    fn from(error: HertzWrenchError) -> Self {
        Self::Law(error)
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn resolve(
    engine: &mut Engine,
    law: &HertzWrenchContact,
    velocity: [f64; 3],
    drive: [f64; 3],
    resting: bool,
    time: f64,
    remaining: f64,
    controls: SpherePlaneHertzSupportControls,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Step, Error> {
    let mobility = [
        engine.inverse_tangent_mass,
        engine.inverse_tangent_mass,
        engine.inertia.recip(),
    ];
    if !resting && let Some(mut step) = pure(law, velocity, drive, mobility, remaining, controls)? {
        count(engine, time, step.dt_s, cancel)?;
        if step.resting {
            // A true zero is not automatically static capture: an opposed
            // drive can cross zero outside the complete static wrench set.
            let demand = input(
                [0.0; 3],
                drive,
                mobility,
                law.parameters().contact_radius_m,
                step.dt_s,
            )?;
            let support = law.advance_phase(
                &demand,
                HertzWrenchPhase::StaticAdmissible,
                controls.wrench,
                &mut || cancel(),
            )?;
            step.sticking = support.mode == HertzWrenchMode::Sticking
                || drive == [0.0; 3] && support.mode == HertzWrenchMode::Frictionless;
        }
        return Ok(step);
    }
    let radius = law.parameters().contact_radius_m;
    let phase = if resting {
        HertzWrenchPhase::StaticAdmissible
    } else {
        HertzWrenchPhase::Kinetic
    };
    // From exactly zero, a constant driving wrench selects a constant
    // acceleration ray by homogeneity. Static balance is constant as well.
    let mut h = if resting {
        remaining
    } else {
        remaining.min(controls.motion.max_sliding_step_s)
    };
    loop {
        count(engine, time, h, cancel)?;
        let query = input(velocity, drive, mobility, radius, h)?;
        let mut config = controls.wrench;
        config.displacement_tolerance_m = config
            .displacement_tolerance_m
            .min(h * controls.velocity_tolerance_m_s)
            .min((h * radius) * controls.spin_tolerance_rad_s);
        let response = law.advance_phase(&query, phase, config, &mut || cancel())?;
        let endpoint = [
            finite(response.slip_displacement_m[0] / h)?,
            finite(response.slip_displacement_m[1] / h)?,
            finite(response.spin_angle_rad / h)?,
        ];
        let captured = response.mode == HertzWrenchMode::Sticking || endpoint == [0.0; 3];
        // A moving phase cannot enter the static/kinetic hysteresis band through
        // a coarse projection. Require the original kinetic set to support a
        // strictly larger complete driving demand. The small outward probe is
        // a conservative numerical interior test, never an increased capacity.
        if captured && !resting {
            let strict_drive = drive.map(|a| a * (1.0 + 1024.0 * f64::EPSILON));
            let demand = input([0.0; 3], strict_drive, mobility, radius, h)?;
            let support =
                law.advance_phase(&demand, HertzWrenchPhase::Kinetic, config, &mut || cancel())?;
            if support.mode != HertzWrenchMode::Sticking
                && !(drive == [0.0; 3] && support.mode == HertzWrenchMode::Frictionless)
            {
                h *= 0.5;
                continue;
            }
        }
        // Velocity and spin vary linearly under this same held wrench. A
        // nonpositive generalized power at both endpoints makes the whole
        // represented segment dissipative, including a captured endpoint.
        // The endpoint law was already checked by the original owner.
        if assists_initial_motion(response.force_n, response.spin_moment_nm, velocity, radius)? {
            h *= 0.5;
            continue;
        }
        return Ok(Step {
            dt_s: h,
            force_n: response.force_n,
            spin_moment_nm: response.spin_moment_nm,
            velocity_after: endpoint,
            resting: captured,
            sticking: captured,
        });
    }
}

fn count(
    engine: &mut Engine,
    time: f64,
    h: f64,
    cancel: &mut impl FnMut() -> bool,
) -> Result<(), SpherePlaneError> {
    poll(cancel)?;
    if engine.sliding_steps >= engine.config.max_sliding_substeps {
        return Err(SpherePlaneError::SlidingBudget);
    }
    engine.sliding_steps += 1;
    progress(time, h)
}

fn input(
    velocity: [f64; 3],
    drive: [f64; 3],
    mobility: [f64; 3],
    radius: f64,
    h: f64,
) -> Result<HertzWrenchInput, SpherePlaneError> {
    let free = [0, 1, 2].map(|i| h * drive[i].mul_add(h, velocity[i]));
    let tangent = finite((h * mobility[0]) * h)?;
    let spin = finite(((h * radius) * mobility[2]) * (h * radius))?;
    if tangent <= 0.0 || spin <= 0.0 {
        return Err(SpherePlaneError::Invalid(
            "unrepresentable Hertz support response",
        ));
    }
    Ok(HertzWrenchInput {
        free_displacement_m: [free[0], free[1]],
        free_spin_angle_rad: free[2],
        scaled_compliance_m_per_n: [[tangent, 0.0, 0.0], [0.0, tangent, 0.0], [0.0, 0.0, spin]],
    })
}

/// Constant native pure-motion wrench and its first true zero-speed event.
/// No finite timestep can step across that event and reflect the velocity.
fn pure(
    law: &HertzWrenchContact,
    velocity: [f64; 3],
    drive: [f64; 3],
    mobility: [f64; 3],
    remaining: f64,
    controls: SpherePlaneHertzSupportControls,
) -> Result<Option<Step>, Error> {
    let speed = velocity[0].hypot(velocity[1]);
    let spin_only = speed == 0.0 && drive[0] == 0.0 && drive[1] == 0.0;
    let mut shear_only = velocity[2] == 0.0 && drive[2] == 0.0 && speed > 0.0;
    if shear_only {
        let drive_length = drive[0].hypot(drive[1]);
        if drive_length > 0.0 {
            let cross = (velocity[0] / speed).mul_add(
                drive[1] / drive_length,
                -(velocity[1] / speed) * (drive[0] / drive_length),
            );
            shear_only = cross.abs() <= 256.0 * f64::EPSILON;
        }
    }
    if !spin_only && !shear_only {
        return Ok(None);
    }
    let p = law.parameters();
    let native = HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: p.normal_force_n,
        contact_radius_m: p.contact_radius_m,
        friction_coefficient: p.kinetic_mu,
    })
    .map_err(HertzWrenchError::Law)?;
    let response = native
        .evaluate([velocity[0], velocity[1]], velocity[2])
        .map_err(HertzWrenchError::Law)?;
    let force = response.force_n[0].hypot(response.force_n[1]);
    if force > controls.wrench.maximum_force_n
        || response.spin_moment_nm.abs() > controls.wrench.maximum_spin_moment_nm
    {
        return Err(HertzWrenchError::EffortBudget {
            force_n: force,
            spin_moment_nm: response.spin_moment_nm.abs(),
        }
        .into());
    }
    let acceleration = [
        mobility[0].mul_add(response.force_n[0], drive[0]),
        mobility[1].mul_add(response.force_n[1], drive[1]),
        mobility[2].mul_add(response.spin_moment_nm, drive[2]),
    ];
    let (magnitude, rate) = if spin_only {
        (velocity[2].abs(), velocity[2].signum() * acceleration[2])
    } else {
        (
            speed,
            (velocity[0] / speed).mul_add(acceleration[0], (velocity[1] / speed) * acceleration[1]),
        )
    };
    let stop = if rate < 0.0 {
        magnitude / -rate
    } else {
        f64::INFINITY
    };
    let h = remaining.min(stop);
    let after = [0, 1, 2].map(|i| acceleration[i].mul_add(h, velocity[i]));
    Ok(Some(Step {
        dt_s: h,
        force_n: response.force_n,
        spin_moment_nm: response.spin_moment_nm,
        velocity_after: after,
        resting: stop <= remaining,
        sticking: false,
    }))
}

fn assists_initial_motion(
    force: [f64; 2],
    moment: f64,
    velocity: [f64; 3],
    radius: f64,
) -> Result<bool, SpherePlaneError> {
    let effort = [force[0], force[1], finite(moment / radius)?];
    let motion = [velocity[0], velocity[1], finite(radius * velocity[2])?];
    let e = effort.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    let v = motion.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    if e == 0.0 || v == 0.0 {
        return Ok(false);
    }
    let power = (effort[0] / e).mul_add(
        motion[0] / v,
        (effort[1] / e).mul_add(motion[1] / v, (effort[2] / e) * (motion[2] / v)),
    );
    Ok(power > 256.0 * f64::EPSILON)
}
