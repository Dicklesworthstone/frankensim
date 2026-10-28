//! First contact of a scalar gap under constant relative acceleration.
//!
//! The gap is g(t) = gap + speed*t + acceleration*t*t/2. A sphere against
//! a translating plane supplies these three scalars from its actual geometry
//! and motion. This query selects the first closing zero, including a return
//! after initially departing from contact. It does not repair penetration.
//! The quadratic formula is analytic, with floating-point evaluation; this is
//! not an interval-certified CCD query for general shapes or curved paths.

use core::fmt;

/// Constant-acceleration gap trajectory in coherent length/time units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BallisticContactInput {
    /// Initial separation; strictly negative values are inadmissible.
    pub gap: f64,
    /// Initial relative gap rate; negative means approaching.
    pub relative_normal_speed: f64,
    /// Constant relative gap acceleration.
    pub relative_normal_acceleration: f64,
    /// Inclusive search horizon, in seconds; zero is allowed.
    pub duration_s: f64,
}

/// The earliest closing contact, or an isolated grazing contact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BallisticContact {
    /// Time measured from the supplied initial state, in seconds.
    pub time_s: f64,
    /// Relative normal velocity on the closing root (zero for grazing).
    pub relative_normal_speed: f64,
}

/// Invalid data or a floating-point range that cannot resolve the event.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BallisticContactError {
    /// A scalar is non-finite, or the search duration is negative.
    InvalidInput,
    /// The supplied configuration already penetrates the surface.
    InitialPenetration,
    /// Coefficient scaling or event evaluation loses a nonzero quantity.
    NumericalRange,
}

impl fmt::Display for BallisticContactError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidInput => {
                "ballistic contact requires finite data and a nonnegative horizon"
            }
            Self::InitialPenetration => "ballistic contact does not repair initial penetration",
            Self::NumericalRange => {
                "ballistic contact event is outside the resolvable numeric range"
            }
        })
    }
}
impl std::error::Error for BallisticContactError {}

/// Finds the first closing zero of the gap within the inclusive horizon.
///
/// An initially touching, departing body is allowed to return later. A touching
/// stationary body is reported at time zero only when acceleration points into
/// the surface; sustained support is the response owner's responsibility.
/// A stationary or departing contact with no future closing root returns None.
///
/// Work and memory are constant. Common coefficient scaling and the stable
/// q quadratic formula avoid overflow of the unscaled discriminant and the
/// cancellation of a small root. Nonzero coefficients lost during scaling,
/// and positive event times rounded to zero, produce a typed refusal.
pub fn first_ballistic_contact(
    input: BallisticContactInput,
) -> Result<Option<BallisticContact>, BallisticContactError> {
    let BallisticContactInput {
        gap,
        relative_normal_speed: speed,
        relative_normal_acceleration: acceleration,
        duration_s,
    } = input;
    if !gap.is_finite()
        || !speed.is_finite()
        || !acceleration.is_finite()
        || !duration_s.is_finite()
        || duration_s < 0.0
    {
        return Err(BallisticContactError::InvalidInput);
    }
    if gap < 0.0 {
        return Err(BallisticContactError::InitialPenetration);
    }
    if gap == 0.0 {
        if speed < 0.0 || (speed == 0.0 && acceleration < 0.0) {
            return Ok(Some(BallisticContact {
                time_s: 0.0,
                relative_normal_speed: speed,
            }));
        }
        if speed > 0.0 && acceleration < 0.0 {
            // The zero root is departing; the next root has velocity -speed.
            return admit_time((speed / -acceleration) * 2.0, -speed, duration_s);
        }
        return Ok(None);
    }
    if duration_s == 0.0 || (speed >= 0.0 && acceleration >= 0.0) {
        return Ok(None);
    }
    if acceleration == 0.0 {
        return admit_time(gap / -speed, speed, duration_s);
    }

    // Scaling all coefficients by the same scalar leaves both roots unchanged.
    let scale = gap.max(speed.abs()).max(acceleration.abs());
    let c = gap / scale;
    let b = speed / scale;
    let scaled_acceleration = acceleration / scale;
    let a = scaled_acceleration * 0.5;
    if c == 0.0 || (speed != 0.0 && b == 0.0) || a == 0.0 {
        return Err(BallisticContactError::NumericalRange);
    }
    let product = (4.0 * a) * c;
    if product == 0.0 {
        return Err(BallisticContactError::NumericalRange);
    }
    // Recover the rounding error of 4*a*c before subtracting. This avoids
    // inventing a grazing root through cancellation in b*b-4*a*c.
    let product_error = (4.0 * a).mul_add(c, -product);
    let discriminant = b.mul_add(b, -product) - product_error;
    if discriminant < 0.0 {
        return Ok(None);
    }
    let root = discriminant.sqrt();
    let q = -0.5 * if b >= 0.0 { b + root } else { b - root };
    // For positive gap, select the earliest closing root. Positive
    // acceleration can reach contact only with negative initial speed.
    let time_s = if acceleration < 0.0 && speed >= 0.0 {
        q / a
    } else {
        c / q
    };
    let impact_speed = -(root * scale);
    admit_time(time_s, impact_speed, duration_s)
}

fn admit_time(
    time_s: f64,
    speed: f64,
    horizon: f64,
) -> Result<Option<BallisticContact>, BallisticContactError> {
    // Infinite positive time is beyond every admitted finite horizon.
    if time_s > horizon {
        return Ok(None);
    }
    if !time_s.is_finite() || time_s <= 0.0 || !speed.is_finite() {
        return Err(BallisticContactError::NumericalRange);
    }
    Ok(Some(BallisticContact {
        time_s,
        relative_normal_speed: speed,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn query(gap: f64, speed: f64, acceleration: f64, horizon: f64) -> Option<BallisticContact> {
        first_ballistic_contact(BallisticContactInput {
            gap,
            relative_normal_speed: speed,
            relative_normal_acceleration: acceleration,
            duration_s: horizon,
        })
        .unwrap()
    }

    fn near(actual: f64, expected: f64) {
        assert!((actual - expected).abs() < 2e-14 * expected.abs().max(1.0));
    }

    #[test]
    fn detects_ballistic_collision_before_both_separated_endpoints() {
        // g(t) = (t-1)(t-3): checking only endpoints at 0 and 4 misses contact.
        let contact = query(3.0, -4.0, 2.0, 4.0).unwrap();
        near(contact.time_s, 1.0);
        near(contact.relative_normal_speed, -2.0);
    }

    #[test]
    fn falling_and_departing_return_trajectories() {
        let falling = query(0.4, -1.0, -10.0, 0.3).unwrap();
        near(falling.time_s, 0.2);
        near(falling.relative_normal_speed, -3.0);
        let returning = query(0.0, 2.0, -10.0, 0.5).unwrap();
        near(returning.time_s, 0.4);
        near(returning.relative_normal_speed, -2.0);
        assert!(query(0.0, 2.0, -10.0, 0.3).is_none());
        near(query(0.5, -0.0, -1.0, 2.0).unwrap().time_s, 1.0);
    }

    #[test]
    fn touching_support_grazing_and_departure_are_distinct() {
        assert_eq!(query(0.0, 0.0, -1.0, 0.0).unwrap().time_s, 0.0);
        assert_eq!(query(0.0, -1.0, 1.0, 0.0).unwrap().time_s, 0.0);
        assert!(query(0.0, 0.0, 1.0, 10.0).is_none());
        assert!(query(0.0, 1.0, 0.0, 10.0).is_none());
        let grazing = query(1.0, -2.0, 2.0, 2.0).unwrap();
        near(grazing.time_s, 1.0);
        near(grazing.relative_normal_speed, 0.0);
        assert!(query(1.0, -1.0, 2.0, 10.0).is_none());
    }

    #[test]
    fn linear_root_and_inclusive_horizon() {
        near(query(1.0, -2.0, 0.0, 0.5).unwrap().time_s, 0.5);
        assert!(query(1.0, -2.0, 0.0, 0.49).is_none());
        assert!(query(1.0, 0.0, 0.0, 1.0).is_none());
    }

    #[test]
    fn small_root_avoids_subtractive_cancellation() {
        near(query(1.0, -1e12, 2.0, 1.0).unwrap().time_s / 1e-12, 1.0);
        near(query(1.0, -1e12, -2.0, 1.0).unwrap().time_s / 1e-12, 1.0);
    }

    #[test]
    fn common_length_scaling_preserves_event_time() {
        for scale in [1e-180, 1.0, 1e180] {
            let hit = query(0.4 * scale, -scale, -10.0 * scale, 1.0).unwrap();
            near(hit.time_s, 0.2);
            near(hit.relative_normal_speed / scale, -3.0);
        }
    }

    #[test]
    fn rejects_penetration_and_unresolved_ranges() {
        let input = BallisticContactInput {
            gap: -1e-16,
            relative_normal_speed: 1.0,
            relative_normal_acceleration: 0.0,
            duration_s: 1.0,
        };
        assert_eq!(
            first_ballistic_contact(input),
            Err(BallisticContactError::InitialPenetration)
        );
        assert_eq!(
            first_ballistic_contact(BallisticContactInput {
                gap: 1.0,
                duration_s: f64::INFINITY,
                ..input
            }),
            Err(BallisticContactError::InvalidInput)
        );
        assert_eq!(
            first_ballistic_contact(BallisticContactInput {
                gap: f64::MIN_POSITIVE,
                relative_normal_speed: -f64::MAX,
                ..input
            }),
            Err(BallisticContactError::NumericalRange)
        );
    }
}
