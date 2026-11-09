//! Analytic original affine rays, with physical event splitting.
//!
//! Translation has f(v)=N(mu0+b|v|), while radius-scaled pure spin has
//! f(v)=N((3*pi/16)mu0+(2/5)b|v|). These are the original Hertz pressure
//! moments. On a fixed-sign segment r=|v| obeys r'=a-B*r. Its exponential,
//! travel and integrated quadratic power are evaluated without time stepping.
//! Clipped translation splits at its actual coefficient-zero boundary; a
//! partly clipped spinning disk is not an affine law and stays with the
//! original pressure integrator. Static decisions still use the joint owner.

use super::*;

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn advance(
    patch: &HertzFrictionPatch,
    input: &HertzDrivenInput,
    original: V,
    initial: V,
    drive: V,
    demand: V,
    factor: Factor,
    config: HertzWrenchConfig,
    quadrature: HertzFrictionQuadrature,
    controls: HertzRateDrivenControls,
    limit: &HertzWrenchContact,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<HertzRateDrivenReport>, Error> {
    let Some(slope_per_speed) = patch.affine_slope() else {
        return Ok(None);
    };
    let mu_zero = patch.law.zero_speed_kinetic_coefficient()?;
    let direction = if initial[2] == 0.0 && det::hypot(initial[0], initial[1]) > 0.0 {
        let magnitude = det::hypot(initial[0], initial[1]);
        [initial[0] / magnitude, initial[1] / magnitude, 0.0]
    } else if initial[0] == 0.0 && initial[1] == 0.0 && initial[2] != 0.0 {
        [0.0, 0.0, 1.0]
    } else if initial == [0.0; 3] && drive[2] == 0.0 && det::hypot(drive[0], drive[1]) > 0.0 {
        let magnitude = det::hypot(drive[0], drive[1]);
        [drive[0] / magnitude, drive[1] / magnitude, 0.0]
    } else if initial == [0.0; 3] && drive[0] == 0.0 && drive[1] == 0.0 && drive[2] != 0.0 {
        [0.0, 0.0, 1.0]
    } else {
        return Ok(None);
    };
    let response = multiply(input.scaled_mobility_per_kg, direction);
    let mobility = dot(direction, response);
    let acceleration = dot(direction, drive);
    if mobility <= 0.0
        || !mobility.is_finite()
        || norm(add(response, scale(direction, -mobility))) > 256.0 * f64::EPSILON * norm(response)
        || norm(add(drive, scale(direction, -acceleration))) > 256.0 * f64::EPSILON * norm(drive)
    {
        return Ok(None);
    }
    let spin = direction[2] != 0.0;
    let normal = patch.normal_force_n();
    let radius = patch.contact_radius_m();
    let dry = checked(normal * mu_zero * if spin { K } else { 1.0 })?;
    let viscous = checked(normal * slope_per_speed * if spin { 0.4 } else { 1.0 })?;
    let cutoff = if slope_per_speed < 0.0 {
        mu_zero / (-slope_per_speed)
    } else {
        f64::INFINITY
    };
    let mut velocity = dot(direction, initial);
    // Relative collinearity alone is insufficient over a long phase: drive
    // and contact transverse impulses can cancel at its endpoint while the
    // interior path leaves the requested rate/travel accuracy. Bound the
    // complete history before consuming any original support evaluation.
    let initial_error = norm(add(original, scale(direction, -velocity)));
    let free_speed = velocity.abs() + acceleration.abs() * input.dt_s;
    let force_bound = dry + viscous.max(0.0) * free_speed;
    let transverse_acceleration = norm(add(drive, scale(direction, -acceleration)))
        + norm(add(response, scale(direction, -mobility))) * force_bound;
    let rate_bound = initial_error + input.dt_s * transverse_acceleration;
    let travel_bound = input.dt_s * (initial_error + 0.5 * input.dt_s * transverse_acceleration);
    if !rate_bound.is_finite()
        || !travel_bound.is_finite()
        || rate_bound > controls.integration.scaled_rate_tolerance_m_s
        || travel_bound > controls.integration.scaled_travel_tolerance_m
    {
        return Ok(None);
    }
    // This sufficient free-speed bound admits the whole spinning phase only
    // while every local speed remains on its original affine branch.
    if spin && slope_per_speed < 0.0 && velocity.abs() + acceleration.abs() * input.dt_s > cutoff {
        return Ok(None);
    }
    let mut mode = input.mode;
    let mut elapsed = 0.0;
    let mut impulse = [0.0; 3];
    let mut travel = [0.0; 3];
    let mut heat = 0.0;
    let mut first_arrest = None;
    let mut steps = 0;
    let mut attempts = 0;
    let mut evaluations = 0;
    let mut iterations = 0;
    let mut maximum_wrench = 0.0_f64;
    let mut maximum_speed = velocity.abs();
    while elapsed < input.dt_s {
        poll(cancelled)?;
        if attempts >= controls.integration.maximum_attempts {
            return Err(TransientError::AttemptBudget.into());
        }
        attempts += 1;
        let remaining = input.dt_s - elapsed;
        if mode == HertzDrivenMode::Resting {
            if evaluations >= controls.integration.maximum_evaluations {
                return Err(TransientError::EvaluationBudget.into());
            }
            let mut c = config;
            c.maximum_evaluations = c
                .maximum_evaluations
                .min(controls.integration.maximum_evaluations - evaluations);
            let static_step = limit.advance_phase(
                &HertzWrenchInput {
                    free_displacement_m: [drive[0], drive[1]],
                    free_spin_angle_rad: checked(drive[2] / radius)?,
                    scaled_compliance_m_per_n: input.scaled_mobility_per_kg,
                },
                HertzWrenchPhase::StaticAdmissible,
                c,
                &mut *cancelled,
            )?;
            evaluations += static_step.evaluations;
            iterations += static_step.iterations;
            let w = vector([
                static_step.force_n[0],
                static_step.force_n[1],
                static_step.spin_moment_nm / radius,
            ])?;
            if static_step.mode == HertzWrenchMode::Sticking || drive == [0.0; 3] {
                maximum_wrench = maximum_wrench.max(norm(w));
                impulse = vector(add(impulse, scale(w, remaining)))?;
                elapsed = input.dt_s;
                velocity = 0.0;
                steps += 1;
                break;
            }
            let released = [
                static_step.slip_displacement_m[0],
                static_step.slip_displacement_m[1],
                radius * static_step.spin_angle_rad,
            ];
            if acceleration == 0.0
                || norm(add(released, scale(direction, -dot(released, direction))))
                    > controls.integration.scaled_rate_tolerance_m_s
            {
                return Err(event(elapsed, scale(direction, velocity)));
            }
            mode = HertzDrivenMode::Sliding;
            if attempts >= controls.integration.maximum_attempts {
                return Err(TransientError::AttemptBudget.into());
            }
            attempts += 1;
        }
        let sign = if velocity == 0.0 {
            acceleration.signum()
        } else {
            velocity.signum()
        };
        let speed = velocity.abs();
        let drive_along = sign * acceleration;
        let clipped = !spin
            && slope_per_speed < 0.0
            && (speed > cutoff || (speed == cutoff && drive_along >= 0.0));
        let (offset, slope) = if clipped { (0.0, 0.0) } else { (dry, viscous) };
        // At the analytically localized clipping boundary the traction is zero.
        let initial_force = if !spin && slope_per_speed < 0.0 && speed == cutoff {
            0.0
        } else {
            checked(offset + slope * speed)?
        };
        if initial_force < 0.0 {
            return Err(event(elapsed, scale(direction, velocity)));
        }
        let derivative = checked(drive_along - mobility * initial_force)?;
        let a = checked(drive_along - mobility * offset)?;
        let b = checked(mobility * slope)?;
        let stop = if speed > 0.0 && a < 0.0 && derivative < 0.0 {
            if b > 0.0 {
                log_ratio(b, speed, -a)? / b
            } else if b < 0.0 {
                log_ratio(-b, speed, -derivative)? / (-b)
            } else {
                speed / (-a)
            }
        } else {
            f64::INFINITY
        };
        let boundary = if !spin && slope_per_speed < 0.0 {
            if clipped && derivative < 0.0 && speed > cutoff {
                (speed - cutoff) / (-derivative)
            } else if !clipped && derivative > 0.0 && speed < cutoff {
                log_ratio(-b, cutoff - speed, derivative)? / (-b)
            } else {
                f64::INFINITY
            }
        } else {
            f64::INFINITY
        };
        let h = remaining.min(stop).min(boundary);
        if !h.is_finite() || h <= 0.0 || elapsed + h == elapsed {
            return Err(event(elapsed, scale(direction, velocity)));
        }
        let (computed, distance, force_integral, loss) =
            flow(speed, a, b, derivative, offset, slope, initial_force, h)?;
        let stopped = h == stop;
        let crossed = h == boundary;
        let after = if stopped {
            0.0
        } else if crossed {
            cutoff
        } else {
            computed
        };
        if after < 0.0
            || (!stopped && after == 0.0)
            || distance < 0.0
            || force_integral < 0.0
            || loss < 0.0
        {
            return Err(event(elapsed, scale(direction, velocity)));
        }
        for rate in [speed, after] {
            if evaluations >= controls.integration.maximum_evaluations {
                return Err(TransientError::EvaluationBudget.into());
            }
            evaluations += 1;
            let native = if rate == 0.0 && !clipped {
                patch
                    .zero_speed_patch
                    .evaluate([direction[0], direction[1]], direction[2] / radius)?
            } else {
                patch.evaluate(
                    [direction[0] * rate, direction[1] * rate],
                    direction[2] * rate / radius,
                    quadrature,
                )?
            };
            if det::hypot(native.force_n[0], native.force_n[1]) > config.maximum_force_n
                || native.spin_moment_nm.abs() > config.maximum_spin_moment_nm
            {
                return Err(TransientError::EffortBudget.into());
            }
            maximum_wrench = maximum_wrench.max(norm([
                native.force_n[0],
                native.force_n[1],
                native.spin_moment_nm / radius,
            ]));
            poll(cancelled)?;
        }
        impulse = vector(add(impulse, scale(direction, -sign * force_integral)))?;
        travel = vector(add(travel, scale(direction, sign * distance)))?;
        heat = checked(heat + loss)?;
        elapsed = if h == remaining {
            input.dt_s
        } else {
            checked(elapsed + h)?
        };
        steps += 1;
        maximum_speed = maximum_speed.max(after);
        velocity = sign * after;
        if stopped {
            mode = HertzDrivenMode::Resting;
            if first_arrest.is_none() {
                first_arrest = Some(elapsed);
            }
        }
    }
    let endpoint = scale(direction, velocity);
    let actual = vector(add(
        add(original, scale(drive, input.dt_s)),
        multiply(input.scaled_mobility_per_kg, impulse),
    ))?;
    if mode == HertzDrivenMode::Sliding {
        let actual_along = dot(actual, direction);
        if actual_along == 0.0 || actual_along.signum() != velocity.signum() {
            return Err(event(elapsed, endpoint));
        }
    }
    let mut trajectory = finish(
        input,
        original,
        endpoint,
        drive,
        demand,
        factor,
        radius,
        impulse,
        travel,
        heat,
        [0.0; 4],
        mode,
        config,
        controls,
        steps,
        0,
        evaluations,
        iterations,
    )?;
    trajectory.first_arrest_time_s = first_arrest;
    trajectory.arrest_time_is_exact = first_arrest.is_some();
    let speed_bound = checked(maximum_speed * (1.0 + PAD))?;
    let bounds = coefficient_bounds(&patch.law, speed_bound)?;
    let mut result = report(
        trajectory,
        speed_bound,
        bounds,
        checked(maximum_wrench * (1.0 + PAD))?,
        false,
        [0.0; 4],
    );
    result.exact_invariant = true;
    poll(cancelled)?;
    Ok(Some(result))
}

// log(1+x*y/z), with a small-argument deterministic series and an overflow-
// safe large-ratio branch. All callers supply positive scalar magnitudes.
fn log_ratio(x: f64, y: f64, z: f64) -> Result<f64, Error> {
    let ratio = crate::slope_product_ratio([x, y], [z]);
    if !ratio.is_finite() {
        return checked(det::ln(x) + det::ln(y) - det::ln(z));
    }
    if ratio < 0.125 {
        let mut power = ratio;
        let mut result = ratio;
        for n in 2..=20 {
            power *= -ratio;
            result += power / n as f64;
        }
        checked(result)
    } else {
        checked(det::ln(1.0 + ratio))
    }
}

// Analytic integrals of r, r^2 and original affine frictional power. Around
// zero exposure the Taylor forms retain small travel/heat terms. Long positive
// exposures use the shared stable exponential relaxation, preserving viscous
// tails instead of subtracting nearly equal initial/absorbed velocities.
#[allow(clippy::too_many_arguments)]
fn flow(
    r: f64,
    a: f64,
    b: f64,
    derivative: f64,
    dry: f64,
    slope: f64,
    force: f64,
    h: f64,
) -> Result<(f64, f64, f64, f64), Error> {
    if derivative == 0.0 {
        let distance = checked(r * h)?;
        return Ok((r, distance, checked(force * h)?, checked(force * distance)?));
    }
    let z = checked(b * h)?;
    let phi1 = if z == 0.0 {
        1.0
    } else {
        checked(-det::expm1(-z) / z)?
    };
    let (phi2, psi) = if z.abs() <= 0.5 {
        let mut phi = 0.5;
        let mut term = 0.5;
        let mut psi = 1.0 / 3.0;
        let mut power = 1.0;
        let mut factorial = 6.0;
        let mut two = 4.0;
        for n in 1..=18 {
            term *= -z / (n + 2) as f64;
            phi += term;
            power *= -z;
            factorial *= (n + 3) as f64;
            two *= 2.0;
            psi += power * (two - 2.0) / factorial;
        }
        (phi, psi)
    } else {
        let second = checked(-det::expm1(-2.0 * z) / (2.0 * z))?;
        ((1.0 - phi1) / z, (1.0 - 2.0 * phi1 + second) / (z * z))
    };
    let (end, distance, square) = if z.abs() <= 0.5 {
        (
            checked(r + derivative * h * phi1)?,
            checked(r * h + derivative * h * h * phi2)?,
            checked(
                r * r * h
                    + 2.0 * r * derivative * h * h * phi2
                    + derivative * derivative * h * h * h * psi,
            )?,
        )
    } else {
        let equilibrium = checked(a / b)?;
        let delta = checked(r - equilibrium)?;
        let end = if b > 0.0 {
            fs_math::decay::relax_toward(r, equilibrium, z)
        } else {
            r + derivative * h * phi1
        };
        let first = h * phi1;
        let second = checked(-det::expm1(-2.0 * z) / (2.0 * b))?;
        (
            checked(end)?,
            checked(equilibrium * h + delta * first)?,
            checked(
                equilibrium * equilibrium * h
                    + 2.0 * equilibrium * delta * first
                    + delta * delta * second,
            )?,
        )
    };
    let (impulse, heat) = if slope >= 0.0 {
        (
            checked(dry * h + slope * distance)?,
            checked(dry * distance + slope * square)?,
        )
    } else {
        // Near coefficient clipping, dry and viscous contributions almost
        // cancel. Integrating from the actual starting force retains the small
        // positive traction/heat generated after entry into the active region.
        (
            checked(force * h + slope * derivative * h * h * phi2)?,
            checked(
                force * r * h
                    + (force + slope * r) * derivative * h * h * phi2
                    + slope * derivative * derivative * h * h * h * psi,
            )?,
        )
    };
    Ok((end, distance, impulse, heat))
}
