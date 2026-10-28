//! Complete quadratic relative paths of two finite spheres.
//!
//! On the normalized clock s in [0,1], r(s)=r0+d*s+q*s². The derivative
//! of |r(s)|²/2 is cubic. Its quadratic derivative separates at most three
//! monotone intervals, so bounded bisection finds every interior distance
//! extremum without a timestep search. Distance, rather than a difference of
//! large squared lengths, determines entry and penetration. The original
//! straight-line query remains the exact zero-acceleration fast path.
//!
//! This is floating analytic geometry, not an interval CCD certificate.
//! Unresolvable scaling or root brackets refuse; positions are never repaired.

use super::{SweptSphereContact, SweptSphereError, SweptSphereInput, dot, norm};

/// Original relative motion and a constant acceleration difference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcceleratedSphereInput {
    pub motion: SweptSphereInput,
    pub relative_acceleration_world: [f64; 3],
}

/// Complete path observation, including possible penetration after first entry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcceleratedSpherePath {
    /// First entering contact. An isolated tangent touch needs no impulse and
    /// is omitted. A stationary initial contact accelerating inward is included;
    /// its sustained reaction must be supplied by a response owner.
    pub first_entering_contact: Option<SweptSphereContact>,
    /// Minimum separation over both endpoints and every interior extremum.
    /// Negative values describe the unimpeded path AFTER entry, not a repair.
    pub minimum_gap: f64,
    /// Number of interior stationary distances, at most three.
    pub interior_extrema: usize,
}

/// Query a complete curved path with constant memory and bounded work.
///
/// At most three cubic roots and one entry root are bisected, each for at most
/// 128 iterations. Entry returns the nonpenetrating side of an adjacent-float
/// bracket; its time error is bounded by that bracket's floating resolution.
/// Native consumers must still check their actual surface witnesses.
pub fn sphere_path(
    input: AcceleratedSphereInput,
) -> Result<AcceleratedSpherePath, SweptSphereError> {
    let m = input.motion;
    let a = input.relative_acceleration_world;
    if !m
        .relative_position_world
        .into_iter()
        .chain(m.relative_velocity_world)
        .chain(a)
        .all(f64::is_finite)
        || !m.radius_sum.is_finite()
        || m.radius_sum <= 0.0
        || !m.duration_s.is_finite()
        || m.duration_s < 0.0
    {
        return Err(SweptSphereError::InvalidInput);
    }
    let initial_distance = norm(m.relative_position_world)?;
    if initial_distance < m.radius_sum {
        return Err(SweptSphereError::InitialPenetration);
    }
    if m.duration_s == 0.0 {
        let normal = m.relative_position_world.map(|x| x / initial_distance);
        let speed = dot(normal, m.relative_velocity_world);
        return Ok(AcceleratedSpherePath {
            first_entering_contact: (initial_distance == m.radius_sum && speed < 0.0).then_some(
                SweptSphereContact {
                    time_s: 0.0,
                    normal_world: normal,
                    relative_normal_speed: speed,
                },
            ),
            minimum_gap: initial_distance - m.radius_sum,
            interior_extrema: 0,
        });
    }
    let displacement = m.relative_velocity_world.map(|v| v * m.duration_s);
    let curvature = a.map(|v| (v * m.duration_s) * (0.5 * m.duration_s));
    if !displacement
        .into_iter()
        .chain(curvature)
        .all(f64::is_finite)
        || m.relative_velocity_world
            .iter()
            .zip(displacement)
            .any(|(&v, x)| v != 0.0 && x == 0.0)
        || a.iter().zip(curvature).any(|(&v, x)| v != 0.0 && x == 0.0)
    {
        return Err(SweptSphereError::NumericalRange);
    }
    let scale = initial_distance
        .max(norm(displacement)?)
        .max(norm(curvature)?);
    let r = m.relative_position_world.map(|v| v / scale);
    let d = displacement.map(|v| v / scale);
    let q = curvature.map(|v| v / scale);
    let radius = m.radius_sum / scale;
    if radius <= 256.0 * f64::EPSILON
        || !scale.is_finite()
        || m.relative_position_world
            .iter()
            .zip(r)
            .any(|(&v, x)| v != 0.0 && x == 0.0)
    {
        return Err(SweptSphereError::NumericalRange);
    }
    let curve = Curve {
        r: m.relative_position_world,
        d: displacement,
        q: curvature,
        scale,
        radius: m.radius_sum,
    };
    if a == [0.0; 3] {
        let square = dot(d, d);
        let s = if square == 0.0 {
            0.0
        } else {
            (-dot(r, d) / square).clamp(0.0, 1.0)
        };
        return Ok(AcceleratedSpherePath {
            first_entering_contact: super::first_sphere_contact(m)?
                .filter(|c| c.relative_normal_speed < 0.0),
            minimum_gap: curve.gap(s)? * scale,
            interior_extrema: usize::from(s > 0.0 && s < 1.0),
        });
    }
    let coefficients = [
        dot(r, d),
        dot(d, d) + 2.0 * dot(r, q),
        3.0 * dot(d, q),
        2.0 * dot(q, q),
    ];
    let mut split = [0.0; 4];
    let mut split_count = 1;
    for x in quadratic_roots(
        3.0 * coefficients[3],
        2.0 * coefficients[2],
        coefficients[1],
    )? {
        if x > 0.0 && x < 1.0 {
            split[split_count] = x;
            split_count += 1;
        }
    }
    split[split_count] = 1.0;
    split_count += 1;
    let derivative = |s: f64| {
        coefficients[3]
            .mul_add(s, coefficients[2])
            .mul_add(s, coefficients[1])
            .mul_add(s, coefficients[0])
    };
    let mut extrema = [0.0; 3];
    let mut count = 0;
    for pair in split[..split_count].windows(2) {
        let [lo, hi] = [pair[0], pair[1]];
        let flo = derivative(lo);
        let fhi = derivative(hi);
        if lo > 0.0 && flo == 0.0 {
            push_root(&mut extrema, &mut count, lo)?;
        }
        if flo != 0.0 && fhi != 0.0 && flo.is_sign_positive() != fhi.is_sign_positive() {
            let root = bisect(lo, hi, flo, derivative)?;
            push_root(&mut extrema, &mut count, root)?;
        }
        if hi < 1.0 && fhi == 0.0 {
            push_root(&mut extrema, &mut count, hi)?;
        }
    }
    let mut times = [0.0; 5];
    times[1..count + 1].copy_from_slice(&extrema[..count]);
    times[count + 1] = 1.0;
    let mut minimum = initial_distance - m.radius_sum;
    let mut first = None;
    let mut previous_gap = curve.gap(0.0)?;
    for pair in times[..count + 2].windows(2) {
        let lo = pair[0];
        let hi = pair[1];
        let gap = curve.gap(hi)?;
        minimum = minimum.min(gap * scale);
        // Between consecutive stationary distances the norm is monotone.
        // A zero minimum followed by increasing distance is only grazing.
        if first.is_none() && (gap < 0.0 || (gap == 0.0 && derivative(hi) < 0.0)) {
            let entry = if previous_gap == 0.0 {
                lo
            } else {
                if previous_gap < 0.0 {
                    return Err(SweptSphereError::NumericalRange);
                }
                entry_bracket(&curve, lo, hi)?
            };
            first = Some(contact(&curve, input, entry)?);
        }
        previous_gap = gap;
    }
    if first.is_none() && initial_distance == m.radius_sum && dot(r, d) < 0.0 {
        first = Some(contact(&curve, input, 0.0)?);
    }
    Ok(AcceleratedSpherePath {
        first_entering_contact: first,
        minimum_gap: minimum,
        interior_extrema: count,
    })
}

struct Curve {
    r: [f64; 3],
    d: [f64; 3],
    q: [f64; 3],
    scale: f64,
    radius: f64,
}
impl Curve {
    fn position(&self, s: f64) -> [f64; 3] {
        core::array::from_fn(|i| self.q[i].mul_add(s, self.d[i]).mul_add(s, self.r[i]))
    }
    fn gap(&self, s: f64) -> Result<f64, SweptSphereError> {
        Ok((norm(self.position(s))? - self.radius) / self.scale)
    }
}

fn contact(
    curve: &Curve,
    input: AcceleratedSphereInput,
    s: f64,
) -> Result<SweptSphereContact, SweptSphereError> {
    let p = curve.position(s);
    let distance = norm(p)?;
    if distance == 0.0 {
        return Err(SweptSphereError::NumericalRange);
    }
    let normal = p.map(|v| v / distance);
    let time = s * input.motion.duration_s;
    let velocity = core::array::from_fn(|i| {
        input.relative_acceleration_world[i].mul_add(time, input.motion.relative_velocity_world[i])
    });
    let speed = dot(normal, velocity);
    if !speed.is_finite() || (s != 0.0 && time == 0.0) {
        return Err(SweptSphereError::NumericalRange);
    }
    Ok(SweptSphereContact {
        time_s: time,
        normal_world: normal,
        relative_normal_speed: speed,
    })
}

fn quadratic_roots(a: f64, b: f64, c: f64) -> Result<[f64; 2], SweptSphereError> {
    if a == 0.0 {
        return Ok(if b == 0.0 {
            [f64::INFINITY; 2]
        } else {
            [-c / b, f64::INFINITY]
        });
    }
    let scale = a.abs().max(b.abs()).max(c.abs());
    let (a, b, c) = (a / scale, b / scale, c / scale);
    let product = 4.0 * a * c;
    let discriminant = b.mul_add(b, -product) - (4.0 * a).mul_add(c, -product);
    if discriminant < 0.0 {
        return Ok([f64::INFINITY; 2]);
    }
    let q = -0.5 * (b + discriminant.sqrt().copysign(b));
    if q == 0.0 {
        return Ok([0.0, f64::INFINITY]);
    }
    let mut roots = [q / a, c / q];
    if roots[0] > roots[1] {
        roots.swap(0, 1);
    }
    if roots[0] == roots[1] {
        roots[1] = f64::INFINITY;
    }
    Ok(roots)
}

fn push_root(roots: &mut [f64; 3], count: &mut usize, value: f64) -> Result<(), SweptSphereError> {
    if *count > 0 && roots[*count - 1] == value {
        return Ok(());
    }
    if *count == roots.len() {
        return Err(SweptSphereError::NumericalRange);
    }
    roots[*count] = value;
    *count += 1;
    Ok(())
}

fn bisect(
    lo: f64,
    hi: f64,
    sign_at_lo: f64,
    f: impl Fn(f64) -> f64,
) -> Result<f64, SweptSphereError> {
    let (mut lo, mut hi) = (lo, hi);
    for _ in 0..128 {
        let mid = lo + 0.5 * (hi - lo);
        if mid == lo || mid == hi {
            return Ok(mid);
        }
        let value = f(mid);
        if value == 0.0 {
            return Ok(mid);
        }
        if value.is_sign_positive() == sign_at_lo.is_sign_positive() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(SweptSphereError::NumericalRange)
}

fn entry_bracket(curve: &Curve, lo: f64, hi: f64) -> Result<f64, SweptSphereError> {
    let (mut lo, mut hi) = (lo, hi);
    for _ in 0..128 {
        let mid = lo + 0.5 * (hi - lo);
        if mid == lo || mid == hi {
            return Ok(lo);
        }
        if curve.gap(mid)? >= 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(SweptSphereError::NumericalRange)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn query(r: [f64; 3], v: [f64; 3], a: [f64; 3], radius: f64, h: f64) -> AcceleratedSpherePath {
        sphere_path(AcceleratedSphereInput {
            motion: SweptSphereInput {
                relative_position_world: r,
                relative_velocity_world: v,
                radius_sum: radius,
                duration_s: h,
            },
            relative_acceleration_world: a,
        })
        .unwrap()
    }
    fn near(a: f64, b: f64) {
        assert!(
            (a - b).abs() < 2e-13 * a.abs().max(b.abs()).max(1.0),
            "{a} != {b}"
        );
    }
    #[test]
    fn g1_accelerated_entry_hidden_by_separated_endpoints_and_return() {
        // x=3-4t+t² crosses radius1 at 2-sqrt2 and leaves again.
        let p = query([3.0, 0.0, 0.0], [-4.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0, 4.0);
        let c = p.first_entering_contact.unwrap();
        near(c.time_s, 2.0 - 2.0_f64.sqrt());
        near(c.relative_normal_speed, -8.0_f64.sqrt());
        near(p.minimum_gap, -1.0);
        assert_eq!(p.interior_extrema, 3);
        // Initially touching and departing, then returns under acceleration.
        let p = query([1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [-2.0, 0.0, 0.0], 1.0, 2.5);
        near(p.first_entering_contact.unwrap().time_s, 2.0);
    }
    #[test]
    fn g1_curved_three_dimensional_contact_has_independent_root() {
        // r(t)=(2-2t, t-t², 3/5). At t=3/5: r=(4/5,6/25,3/5).
        let radius = (1.0_f64 + 36.0 / 625.0).sqrt();
        let p = query(
            [2.0, 0.0, 0.6],
            [-2.0, 1.0, 0.0],
            [0.0, -2.0, 0.0],
            radius,
            1.0,
        );
        let c = p.first_entering_contact.unwrap();
        near(c.time_s, 0.6);
        near(c.normal_world[0], 0.8 / radius);
        near(c.normal_world[1], 0.24 / radius);
        near(c.relative_normal_speed, (-1.6 - 0.048) / radius);
        assert!(p.minimum_gap < 0.0);
    }
    #[test]
    fn g1_grazing_stationary_and_inward_support_are_distinct() {
        let p = query([2.0, 1.0, 0.0], [-4.0, 0.0, 0.0], [4.0, 0.0, 0.0], 1.0, 2.0);
        assert!(p.first_entering_contact.is_none());
        near(p.minimum_gap, 0.0);
        let p = query([1.0, 0.0, 0.0], [0.0; 3], [-1.0, 0.0, 0.0], 1.0, 0.5);
        let c = p.first_entering_contact.unwrap();
        near(c.time_s, 0.0);
        near(c.relative_normal_speed, 0.0);
        let p = query([1.0, 0.0, 0.0], [0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.5);
        assert!(p.first_entering_contact.is_none());
    }
    #[test]
    fn g3_scaling_and_clock_refinement_preserve_curved_entry() {
        for scale in [1e-100, 1.0, 1e100] {
            let p = query(
                [3.0 * scale, 0.0, 0.0],
                [-4.0 * scale, 0.0, 0.0],
                [2.0 * scale, 0.0, 0.0],
                scale,
                4.0,
            );
            near(
                p.first_entering_contact.unwrap().time_s,
                2.0 - 2.0_f64.sqrt(),
            );
            near(p.minimum_gap / scale, -1.0);
        }
        for h in [0.75, 1.5, 3.0] {
            let p = query([3.0, 0.0, 0.0], [-4.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0, h);
            near(
                p.first_entering_contact.unwrap().time_s,
                2.0 - 2.0_f64.sqrt(),
            );
        }
    }
}
