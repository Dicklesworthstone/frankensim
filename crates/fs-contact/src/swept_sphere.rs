//! First contact between spheres with constant relative centre velocity.
//!
//! Common acceleration cancels from relative motion. The caller supplies the
//! actual centre difference and relative velocity; rotation does not change
//! spherical geometry. Closest-approach/chord geometry avoids subtracting two
//! large squared distances. This is analytic floating-point detection, not an
//! interval-certified CCD query for arbitrary shapes or curved relative paths.

use core::fmt;
use fs_math::det;

/// Relative straight-line trajectory in coherent length/time units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SweptSphereInput {
    /// Centre A minus centre B, in the world frame.
    pub relative_position_world: [f64; 3],
    /// Centre velocity A minus centre velocity B, in the world frame.
    pub relative_velocity_world: [f64; 3],
    /// Sum of the two positive sphere radii.
    pub radius_sum: f64,
    /// Inclusive finite search horizon, in seconds; zero is allowed.
    pub duration_s: f64,
}

/// First closing or isolated grazing contact along the relative line.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SweptSphereContact {
    /// Time from the initial configuration, in seconds.
    pub time_s: f64,
    /// Unit normal from B toward A at contact.
    pub normal_world: [f64; 3],
    /// Relative normal centre velocity, nonpositive at entry.
    pub relative_normal_speed: f64,
}

/// Invalid geometry, initial overlap, or unresolvable event arithmetic.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SweptSphereError {
    /// Inputs must be finite, radius positive, and horizon nonnegative.
    InvalidInput,
    /// Initial centre distance is smaller than the radius sum.
    InitialPenetration,
    /// The event cannot be resolved within the floating-point range.
    NumericalRange,
}
impl fmt::Display for SweptSphereError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidInput => {
                "swept spheres require finite data, positive radius and nonnegative horizon"
            }
            Self::InitialPenetration => "swept spheres do not repair initial overlap",
            Self::NumericalRange => "swept sphere contact is outside the resolvable numeric range",
        })
    }
}
impl std::error::Error for SweptSphereError {}

/// Find first sphere contact without stepping over the incoming trajectory.
///
/// Initially touching, closing spheres return time zero. A touching pair with
/// zero or positive normal speed has no new entering event. A tangent encounter
/// from positive separation returns a zero-normal-speed grazing event.
///
/// Work and storage are constant. All direction geometry is dimensionless;
/// neither physical distances nor velocities are squared. Positive separation
/// and event time are never clamped to zero. A radius below the floating-point
/// resolution of the initial separation refuses instead of claiming a contact.
pub fn first_sphere_contact(
    input: SweptSphereInput,
) -> Result<Option<SweptSphereContact>, SweptSphereError> {
    let r = input.relative_position_world;
    let v = input.relative_velocity_world;
    let radius = input.radius_sum;
    if !r.into_iter().chain(v).all(f64::is_finite)
        || !radius.is_finite()
        || radius <= 0.0
        || !input.duration_s.is_finite()
        || input.duration_s < 0.0
    {
        return Err(SweptSphereError::InvalidInput);
    }
    let distance = norm(r)?;
    if distance < radius {
        return Err(SweptSphereError::InitialPenetration);
    }
    let speed = norm(v)?;
    if speed == 0.0 {
        return Ok(None);
    }
    let radial = r.map(|x| x / distance);
    let direction = v.map(|x| x / speed);
    let projection = dot(radial, direction);
    if projection >= 0.0 {
        return Ok(None);
    }
    if distance == radius {
        let normal_speed = finite(speed * projection)?;
        if normal_speed == 0.0 {
            return Err(SweptSphereError::NumericalRange);
        }
        return Ok(Some(SweptSphereContact {
            time_s: 0.0,
            normal_world: radial,
            relative_normal_speed: normal_speed,
        }));
    }
    if input.duration_s == 0.0 {
        return Ok(None);
    }
    let radius_ratio = radius / distance;
    if radius_ratio <= 256.0 * f64::EPSILON {
        return Err(SweptSphereError::NumericalRange);
    }
    // Triple product gives the component perpendicular to the trajectory;
    // head-on paths retain an exact zero transverse vector.
    let perpendicular = cross(direction, cross(radial, direction));
    let transverse = norm(perpendicular)?;
    if transverse > radius_ratio {
        return Ok(None);
    }
    let fraction = transverse / radius_ratio;
    let chord_fraction = det::sqrt((1.0 - fraction) * (1.0 + fraction));
    let half_chord = radius_ratio * chord_fraction;
    // Rationalized entry root: (d^2-R^2)/(-r.u + sqrt(R^2-r_perp^2)).
    let gap_ratio = (distance - radius) / distance;
    let travel_ratio = finite(gap_ratio * (1.0 + radius_ratio) / (-projection + half_chord))?;
    let time_s = (travel_ratio * distance) / speed;
    if time_s > input.duration_s {
        return Ok(None);
    }
    if !time_s.is_finite() || time_s <= 0.0 {
        return Err(SweptSphereError::NumericalRange);
    }
    let contact_direction = core::array::from_fn(|i| {
        (-direction[i]).mul_add(chord_fraction, perpendicular[i] / radius_ratio)
    });
    let contact_norm = norm(contact_direction)?;
    if contact_norm == 0.0 {
        return Err(SweptSphereError::NumericalRange);
    }
    let normal_world = contact_direction.map(|x| x / contact_norm);
    let relative_normal_speed = finite(-speed * chord_fraction)?;
    if chord_fraction > 0.0 && relative_normal_speed == 0.0 {
        return Err(SweptSphereError::NumericalRange);
    }
    Ok(Some(SweptSphereContact {
        time_s,
        normal_world,
        relative_normal_speed,
    }))
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -a[2] * b[1]),
        a[2].mul_add(b[0], -a[0] * b[2]),
        a[0].mul_add(b[1], -a[1] * b[0]),
    ]
}
fn norm(v: [f64; 3]) -> Result<f64, SweptSphereError> {
    finite(det::hypot(det::hypot(v[0], v[1]), v[2]))
}
fn finite(value: f64) -> Result<f64, SweptSphereError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SweptSphereError::NumericalRange)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn query(r: [f64; 3], v: [f64; 3], radius: f64, duration: f64) -> Option<SweptSphereContact> {
        first_sphere_contact(SweptSphereInput {
            relative_position_world: r,
            relative_velocity_world: v,
            radius_sum: radius,
            duration_s: duration,
        })
        .unwrap()
    }
    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2e-13 * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn detects_head_on_crossing_between_separated_endpoints() {
        let hit = query([-3.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 6.0).unwrap();
        near(hit.time_s, 1.0);
        assert_eq!(hit.normal_world, [-1.0, 0.0, 0.0]);
        near(hit.relative_normal_speed, -1.0);
        assert!(query([-3.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 0.9).is_none());
        near(
            query([-3.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 1.0)
                .unwrap()
                .time_s,
            1.0,
        );
    }

    #[test]
    fn oblique_entry_normal_and_grazing_are_geometric() {
        let hit = query([-3.0, 0.6, 0.0], [2.0, 0.0, 0.0], 1.0, 2.0).unwrap();
        near(hit.time_s, 1.1);
        near(hit.normal_world[0], -0.8);
        near(hit.normal_world[1], 0.6);
        near(hit.relative_normal_speed, -1.6);
        let grazing = query([-3.0, 2.0, 0.0], [1.0, 0.0, 0.0], 2.0, 4.0).unwrap();
        near(grazing.time_s, 3.0);
        near(grazing.relative_normal_speed, 0.0);
        near(grazing.normal_world[1], 1.0);
        assert!(query([-3.0, 2.1, 0.0], [1.0, 0.0, 0.0], 2.0, 8.0).is_none());
    }

    #[test]
    fn initial_closing_departing_and_stationary_contacts_differ() {
        assert_eq!(
            query([2.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 2.0, 0.0)
                .unwrap()
                .time_s,
            0.0
        );
        assert!(query([2.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 2.0).is_none());
        assert!(query([2.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2.0, 2.0).is_none());
        assert!(query([3.0, 0.0, 0.0], [0.0, 0.0, 0.0], 2.0, 2.0).is_none());
    }

    #[test]
    fn spatial_rotation_and_length_scaling_preserve_event() {
        for scale in [1e-150, 1.0, 1e150] {
            let hit = query(
                [0.6 * scale, 0.0, -3.0 * scale],
                [0.0, 0.0, 2.0 * scale],
                scale,
                2.0,
            )
            .unwrap();
            near(hit.time_s, 1.1);
            near(hit.normal_world[0], 0.6);
            near(hit.normal_world[2], -0.8);
            near(hit.relative_normal_speed / scale, -1.6);
        }
    }

    #[test]
    fn near_touching_entry_avoids_root_cancellation() {
        let gap = 1e-10;
        let hit = query([1.0 + gap, 0.0, 0.0], [-1.0, 0.0, 0.0], 1.0, 1.0).unwrap();
        near(hit.time_s / ((1.0 + gap) - 1.0), 1.0);
    }

    #[test]
    fn rejects_overlap_invalid_and_unresolved_geometry() {
        let input = SweptSphereInput {
            relative_position_world: [1.0, 0.0, 0.0],
            relative_velocity_world: [-1.0, 0.0, 0.0],
            radius_sum: 2.0,
            duration_s: 1.0,
        };
        assert_eq!(
            first_sphere_contact(input),
            Err(SweptSphereError::InitialPenetration)
        );
        assert_eq!(
            first_sphere_contact(SweptSphereInput {
                radius_sum: 0.0,
                ..input
            }),
            Err(SweptSphereError::InvalidInput)
        );
        assert_eq!(
            first_sphere_contact(SweptSphereInput {
                relative_position_world: [1e20, 0.0, 0.0],
                radius_sum: 1.0,
                duration_s: 1e21,
                ..input
            }),
            Err(SweptSphereError::NumericalRange)
        );
    }
}
