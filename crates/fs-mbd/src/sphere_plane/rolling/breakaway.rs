//! First exit of the material trajectory's static-traction disk.

use super::*;

/// Physical data at the handoff from feasible material rolling.
///
/// The native motion ends with a zero-impulse `Sliding` event and zero point
/// slip. These are instantaneous boundary data, not a law for subsequent slip.
#[derive(Clone, Copy, Debug)]
pub struct SphereRollingBreakaway {
    /// Traction that continued no-slip motion would require [N]. At an
    /// immediate handoff this may already exceed the available static force.
    pub required_static_traction_world_n: Vec3,
    /// Original static capacity at the held normal load [N].
    pub available_static_traction_n: f64,
    /// Actual material torque at the boundary [N m]. The next physical owner
    /// must account for deformation loss; it cannot silently discard this law.
    pub rolling_moment_world_nm: Vec3,
}

/// A complete material-rolling interval or its accepted prefix at breakaway.
#[derive(Debug)]
pub struct SphereRollingEventAdvance<'a> {
    /// Original native state, impulses and work/heat for the accepted prefix.
    /// An immediate handoff preserves the input pose and momenta exactly and
    /// reports zero impulse, work and heat. Only contact-mode metadata changes.
    pub rolling: ViscoelasticSphereAdvance<'a>,
    /// Duration actually advanced, possibly zero at immediate breakaway [s].
    pub accepted_duration_s: f64,
    /// Requested duration still requiring another physical contact owner [s].
    pub remaining_duration_s: f64,
    /// `None` means the complete request stayed feasible. `Some` includes a
    /// native `Sliding` event at `rolling.motion.state_after.time_s()` even
    /// when the event coincides with the requested endpoint.
    pub breakaway: Option<SphereRollingBreakaway>,
}

/// Advance until the first static-traction breakaway or the requested horizon.
///
/// The held force, world torque and material coefficient give an exact
/// exponential traction path. Its first cone exit is located analytically,
/// including a later exit after initially unloading from the static boundary.
/// Floating-point evaluation retains the feasible side of the event; this is
/// not an interval-certified event enclosure. No trajectory timestep controls
/// the event time, translation, momenta or material heat.
///
/// The original source and no-slip requirements still apply to the entire
/// accepted prefix. A source limit reached before breakaway refuses; no
/// extrapolation supplies a sliding or thermal model. All published events
/// count against `config.max_events`, including an immediate handoff.
/// Cancellation or any later refusal publishes neither the prefix nor heat.
///
/// After breakaway, the returned state has `Sliding` mode with zero point slip
/// at onset. Use the explicit remaining duration with a separately admitted
/// deformation/sliding owner. Passing it to an ideal point-contact trajectory
/// would omit the material loss and is not a continuation supplied by this API.
#[allow(clippy::too_many_arguments)]
pub fn advance_viscoelastic_sphere_rolling_until_breakaway<'a>(
    scene: &SpherePlane,
    state: SpherePlaneState,
    external_torque_world: Vec3,
    support_friction: &FrictionLaw,
    binding: SphereRollingBinding<'a>,
    duration_s: f64,
    config: SpherePlaneConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<SphereRollingEventAdvance<'a>, SphereRollingError> {
    advance_rolling(
        scene,
        state,
        external_torque_world,
        support_friction,
        binding,
        duration_s,
        config,
        true,
        cancelled,
    )
}

pub(super) fn handoff(
    traction: Vec3,
    capacity: f64,
    rolling: Vec3,
    moment_coefficient: f64,
) -> Result<SphereRollingBreakaway, SphereRollingError> {
    let moment = rolling.scale(-moment_coefficient);
    vector(traction)?;
    vector(moment)?;
    Ok(SphereRollingBreakaway {
        required_static_traction_world_n: traction,
        available_static_traction_n: capacity,
        rolling_moment_world_nm: moment,
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn first_breakaway(
    acceleration: Vec3,
    mass: f64,
    force: Vec3,
    capacity: f64,
    rate: f64,
    horizon: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<f64>, SphereRollingError> {
    poll(cancelled)?;
    let inertial = acceleration.scale(mass);
    let initial = inertial.sub(force);
    let initial_norm = norm(initial)?;
    if initial_norm > capacity {
        return Ok(Some(0.0));
    }
    // Traction is q(t)=m*a0*exp(-rate*t)-F_t. Its path is a line segment
    // from q(0) toward -F_t, so two feasible endpoints exclude every exit.
    // A terminal point on the disk is approached asymptotically, not crossed.
    if rate == 0.0 || inertial == Vec3::ZERO || norm(force)? <= capacity {
        return Ok(None);
    }
    if capacity == 0.0 {
        return Ok(Some(0.0));
    }
    let magnitude = norm(inertial)?;
    let axis = Vec3::new(
        inertial.x / magnitude,
        inertial.y / magnitude,
        inertial.z / magnitude,
    );
    let start = Vec3::new(
        initial.x / capacity,
        initial.y / capacity,
        initial.z / capacity,
    );
    let axis_square = axis.norm_squared();
    // Use distance along the traction line in units of the capacity. The
    // first disk exit is at most two such units away. This avoids squaring a
    // tiny capacity/drive ratio and losing an early positive event to zero.
    let crossing = first_ballistic_contact(BallisticContactInput {
        gap: ((capacity - initial_norm) / capacity) * (1.0 + initial_norm / capacity),
        relative_normal_speed: 2.0 * start.dot(axis),
        relative_normal_acceleration: -2.0 * axis_square,
        duration_s: 4.0,
    })
    .map_err(SpherePlaneError::Detection)?
    .ok_or(SpherePlaneError::Unresolved(
        "rolling traction exit geometry",
    ))?;
    if crossing.time_s == 0.0 {
        return Ok(Some(0.0));
    }
    let fraction = scaled_distance(crossing.time_s, capacity, magnitude)?;
    let clock = if fraction <= 0.5 {
        negative_log_one_minus(fraction)
    } else {
        // Close to terminal traction, solve for the exponential itself.
        // Subtracting a nearly unit traveled fraction would lose a late exit.
        let terminal_norm = norm(force)?;
        let terminal = Vec3::new(force.x / capacity, force.y / capacity, force.z / capacity);
        let reverse = first_ballistic_contact(BallisticContactInput {
            gap: ((terminal_norm - capacity) / capacity) * (1.0 + terminal_norm / capacity),
            relative_normal_speed: -2.0 * terminal.dot(axis),
            relative_normal_acceleration: 2.0 * axis_square,
            duration_s: 4.0,
        })
        .map_err(SpherePlaneError::Detection)?
        .ok_or(SpherePlaneError::Unresolved(
            "rolling terminal traction exit",
        ))?;
        let decay = scaled_distance(reverse.time_s, capacity, magnitude)?;
        -det::ln(decay)
    };
    if !clock.is_finite() || clock <= 0.0 {
        return Err(SpherePlaneError::Unresolved("rolling breakaway clock").into());
    }
    let mut time = clock / rate;
    let traction = |time: f64| {
        // Keep the same arithmetic as the original full-interval admission.
        let scaled_time = rate * time;
        required_traction(
            acceleration,
            mass,
            force,
            scaled_time,
            det::exp(-scaled_time),
        )
    };
    if time > horizon {
        if norm(traction(horizon))? <= capacity {
            return Ok(None);
        }
        time = horizon;
    }
    if !time.is_finite() || time <= 0.0 {
        return Err(SpherePlaneError::Unresolved("rolling breakaway time").into());
    }
    // The algebraic root and exp/norm each round independently. Preserve the
    // original strict cone admission: do not enlarge friction to admit a root
    // rounded to the infeasible side. Usually a neighboring float suffices.
    for _ in 0..8 {
        poll(cancelled)?;
        if norm(traction(time))? <= capacity {
            return Ok(Some(time));
        }
        time = time.next_down();
    }
    let mut upper = time;
    let mut lower = 0.5 * upper;
    for _ in 0..64 {
        poll(cancelled)?;
        if norm(traction(lower))? <= capacity {
            break;
        }
        lower *= 0.5;
    }
    if norm(traction(lower))? > capacity {
        return Err(SpherePlaneError::Unresolved("rolling breakaway feasible bracket").into());
    }
    for _ in 0..64 {
        poll(cancelled)?;
        let middle = f64::midpoint(lower, upper);
        if middle == lower || middle == upper {
            if lower <= 0.0 {
                return Err(SpherePlaneError::NoProgress.into());
            }
            return Ok(Some(lower));
        }
        if norm(traction(middle))? <= capacity {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    Err(SpherePlaneError::Unresolved("rolling breakaway rounding").into())
}

fn scaled_distance(
    distance: f64,
    capacity: f64,
    magnitude: f64,
) -> Result<f64, SphereRollingError> {
    let value = if capacity <= magnitude {
        distance * (capacity / magnitude)
    } else {
        (distance * capacity) / magnitude
    };
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(SpherePlaneError::Unresolved("rolling breakaway traction scale").into())
    }
}

// -ln(1-rho) for 0 < rho <= 1/2. The atanh series has ratio <= 1/9;
// fixed terms through degree 35 leave less than one rounding unit. Factoring
// rho outside retains a positive subnormal that rho/2 could otherwise erase.
fn negative_log_one_minus(fraction: f64) -> f64 {
    let denominator = 2.0 - fraction;
    let square = (fraction / denominator).powi(2);
    let mut polynomial: f64 = 1.0 / 35.0;
    for order in (0..17).rev() {
        polynomial = polynomial.mul_add(square, 1.0 / (2 * order + 1) as f64);
    }
    fraction * (2.0 / denominator) * polynomial
}

#[cfg(test)]
mod tests {
    use super::*;

    fn locate(inertial: Vec3, force: Vec3, capacity: f64, horizon: f64) -> Option<f64> {
        first_breakaway(inertial, 1.0, force, capacity, 1.0, horizon, &mut || false).unwrap()
    }

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2e-13 * expected.abs(),
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn an_inward_boundary_trajectory_keeps_its_positive_return_time() {
        let inertial = Vec3::new(2.5, 0.0, 0.0);
        let force = Vec3::new(1.5, 0.0, 0.0);
        // q(0)=+1 and q(infinity)=-1.5. The first outward crossing is
        // exp(-t)=0.2, after unloading and reversing static traction.
        assert_eq!(locate(inertial, force, 1.0, 1e-20), None);
        assert_eq!(locate(inertial, force, 1.0, 1.0), None);
        close(locate(inertial, force, 1.0, 2.0).unwrap(), 5.0_f64.ln());
        assert_eq!(locate(Vec3::new(0.5, 0.0, 0.0), force, 1.0, 2.0), Some(0.0));
        // Reaching terminal capacity is asymptotic, never a finite event.
        assert_eq!(
            locate(
                Vec3::new(2.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                1.0,
                1000.0
            ),
            None
        );
    }

    #[test]
    fn early_and_late_breakaway_clocks_preserve_force_scale() {
        for scale in [1e-100, 1.0, 1e100] {
            let force = Vec3::new(scale, 0.0, 0.0);
            for fraction in [1e-12, 0.4, 1.0 - 1e-12] {
                let capacity = fraction * scale;
                let expected = if fraction < 0.5 {
                    -(-capacity / scale).ln_1p()
                } else {
                    -((scale - capacity) / scale).ln()
                };
                let time = locate(force, force, capacity, 100.0).unwrap();
                close(time, expected);
                let traction = required_traction(force, 1.0, force, time, det::exp(-time));
                assert!(norm(traction).unwrap() <= capacity);
                close(norm(traction).unwrap(), capacity);
            }
        }
    }

    #[test]
    fn a_turning_traction_path_uses_the_full_disk_and_elastic_paths_do_not_decay() {
        let inertial = Vec3::new(2.0, 0.75, 0.0);
        let force = Vec3::new(2.0, 0.0, 0.0);
        let decay = (8.0 - 9.25_f64.sqrt()) / 9.125;
        close(locate(inertial, force, 1.0, 2.0).unwrap(), -decay.ln());
        assert_eq!(
            first_breakaway(inertial, 1.0, force, 1.0, 0.0, 1000.0, &mut || false).unwrap(),
            None
        );
    }
}
