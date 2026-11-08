//! Dimensionless, monotone radius solves for the original I/III potential.
//! y=a/a_s, d=delta/h, b=8 Gbar/E* (R phi/a_s)^2. Stationarity is
//! (y^2-d)^2+b y^2=16y. No dimensional material or geometry is replaced here.

#[derive(Clone, Copy, Debug)]
pub(super) struct Boundary {
    pub radius: f64,
    pub twist_square: f64,
    pub radial_fold: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Unresolved;

/// `distance` is formed dimensionally as (delta-delta_s)/h, so an input
/// adjacent to the normal fold does not lose its distance by adding 3 to d.
pub(super) fn boundary(d: f64, distance: f64) -> Result<Option<Boundary>, Unresolved> {
    if !d.is_finite() || !distance.is_finite() || distance < 0.0 {
        return Err(Unresolved);
    }
    if d == 0.0 {
        return Ok(None);
    }
    if d > 0.0 {
        let radius = d.sqrt();
        return Ok(Some(Boundary {
            radius,
            twist_square: 16.0 / radius,
            radial_fold: false,
        }));
    }
    // The unique fold solves f^4+8f=d^2. Near d=-3 retain z=1-f
    // and rationalize b_f=8/f-2f^2+2d instead of subtracting equal terms.
    let (radius, offset) = if distance < 1.5 {
        let target = distance * (6.0 - distance);
        let z = root(target, target / 9.0, |z| {
            (
                z * (12.0 + z * (-6.0 + z * (4.0 - z))),
                12.0 + z * (-12.0 + z * (12.0 - 4.0 * z)),
            )
        })?;
        (1.0 - z, z)
    } else {
        let target = d * d;
        if target == 0.0 {
            return Err(Unresolved);
        }
        let f = root(target, target / 8.0, |f| {
            (f * (8.0 + f * f * f), 8.0 + 4.0 * f * f * f)
        })?;
        (f, 1.0 - f)
    };
    let f = radius;
    let square = 32.0 * offset * (1.0 + f + f * f) / (f * (4.0 - f * f * f - f * d));
    if !f.is_normal() || f <= 0.0 || !square.is_finite() || square < 0.0 {
        return Err(Unresolved);
    }
    Ok(Some(Boundary {
        radius: f,
        twist_square: square,
        radial_fold: true,
    }))
}

/// Return both radius and its positive offset from the physical boundary.
/// Distance is b_limit*(1-|phi|/phi_limit)*(1+|phi|/phi_limit).
pub(super) fn attached(
    d: f64,
    b: f64,
    boundary: Boundary,
    distance: f64,
) -> Result<(f64, f64), Unresolved> {
    let f = boundary.radius;
    // Near zero approach b_limit diverges. Away from that remote boundary,
    // evaluate the original quartic using the ORIGINAL b, not b_limit-distance.
    if boundary.twist_square > 64.0 && distance > boundary.twist_square / 64.0 {
        let y = direct_quartic(d, b, f)?;
        return Ok((y, y - f));
    }
    let mut upper = distance.sqrt();
    // These positive terms give tight upper brackets even when f is tiny:
    // fold: D >= (8/f)[t/(f+t)]^2; edge: D >= (16/f)t/(f+t).
    let ratio = if boundary.radial_fold {
        (distance * (f / 8.0)).sqrt()
    } else {
        distance * (f / 16.0)
    };
    if ratio < 1.0 {
        upper = upper.min(f * ratio / (1.0 - ratio));
    }
    let offset = root(distance, upper, |t| {
        let y = f + t;
        let ratio = f / y;
        let common = (1.0 + ratio) * (1.0 + ratio * ratio);
        if boundary.radial_fold {
            (
                t * t * ((1.0 + ratio).powi(2) + 8.0 / (f * y * y)),
                2.0 * t * (common + 8.0 / (y * y * y)),
            )
        } else {
            (
                (t * (1.0 + ratio)).powi(2) + (16.0 / f) * (t / y),
                2.0 * t * common + 16.0 / (y * y),
            )
        }
    })?;
    Ok((f + offset, offset))
}

fn direct_quartic(d: f64, b: f64, lower: f64) -> Result<f64, Unresolved> {
    if !b.is_finite() || b < 0.0 {
        return Err(Unresolved);
    }
    // Both bounds independently have q>=0. q has exactly one root above the
    // physical lower boundary, even if its derivative is initially negative.
    let mut hi = (d.max(0.0).sqrt() + 16.0_f64.cbrt()).min(16.0 / b) * (1.0 + 32.0 * f64::EPSILON);
    let mut lo = lower;
    let mut y = hi;
    for _ in 0..64 {
        let opening = y * y - d;
        let residual = opening * opening + (b * y) * y - 16.0 * y;
        if residual.abs() <= 16.0 * f64::EPSILON * (16.0 * y) {
            return Ok(y);
        }
        if residual > 0.0 {
            hi = y;
        } else {
            lo = y;
        }
        let derivative = 4.0 * y * opening + 2.0 * b * y - 16.0;
        let newton = y - residual / derivative;
        let next = if newton > lo && newton < hi {
            newton
        } else {
            f64::midpoint(lo, hi)
        };
        if next == y {
            break;
        }
        y = next;
    }
    Err(Unresolved)
}

/// At zero approach y^3+b y=16 has one positive root and no finite fold.
pub(super) fn zero_approach(b: f64) -> Result<f64, Unresolved> {
    if !b.is_finite() || b < 0.0 {
        return Err(Unresolved);
    }
    root(16.0, 16.0_f64.cbrt().min(16.0 / b), |y| {
        (y * (y * y + b), 3.0 * y * y + b)
    })
}

fn root(target: f64, upper: f64, f: impl Fn(f64) -> (f64, f64)) -> Result<f64, Unresolved> {
    if target == 0.0 {
        return Ok(0.0);
    }
    let mut hi = upper * (1.0 + 32.0 * f64::EPSILON);
    if !target.is_normal() || target < 0.0 || !hi.is_finite() || hi <= 0.0 || f(hi).0 < target {
        return Err(Unresolved);
    }
    let mut lo = 0.0;
    let mut x = hi;
    for _ in 0..64 {
        let (value, derivative) = f(x);
        let residual = value - target;
        if residual.abs() <= 16.0 * f64::EPSILON * target {
            return Ok(x);
        }
        if residual > 0.0 {
            hi = x;
        } else {
            lo = x;
        }
        let newton = x - residual / derivative;
        let next = if newton > lo && newton < hi {
            newton
        } else {
            f64::midpoint(lo, hi)
        };
        if next == x {
            break;
        }
        x = next;
    }
    Err(Unresolved)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stable_radii_satisfy_original_quartic_and_boundary() {
        for d in [-3.0, -2.99, -2.0, -1.0, -0.01, 0.0, 0.1, 1.0, 10.0] {
            let boundary = boundary(d, d + 3.0).unwrap();
            for fraction in [0.0, 0.1, 0.8, 0.99999, 1.0] {
                let (y, b) = if let Some(boundary) = boundary {
                    let b = boundary.twist_square * fraction * fraction;
                    let distance = boundary.twist_square * (1.0 - fraction) * (1.0 + fraction);
                    let (y, _) = attached(d, b, boundary, distance).unwrap();
                    assert!(y >= boundary.radius);
                    (y, b)
                } else {
                    let b = 100.0 * fraction;
                    (zero_approach(b).unwrap(), b)
                };
                let residual = (y * y - d).powi(2) + b * y * y - 16.0 * y;
                assert!(
                    residual.abs() <= 2e-13 * 16.0 * y,
                    "d={d}, b={b}, y={y}, r={residual}"
                );
                assert!(y * y >= d - 8.0 * f64::EPSILON * d.abs());
            }
        }
    }

    #[test]
    fn distance_to_fold_survives_rounded_dimensionless_approach() {
        let b = boundary(-3.0, 1e-18).unwrap().unwrap();
        assert!(b.twist_square > 0.0);
        assert!((b.twist_square / 8e-18 - 1.0).abs() < 1e-14);
        let (y, offset) = attached(-3.0, b.twist_square, b, b.twist_square * f64::EPSILON).unwrap();
        assert!(offset > 0.0 && y >= b.radius);
        assert!(boundary(-1e-200, 3.0).is_err());
    }

    #[test]
    fn either_side_of_zero_approach_retains_original_quartic_accuracy() {
        for d in [-1e-3, -1e-8, -1e-14, 1e-3, 1e-8, 1e-14] {
            let boundary = boundary(d, d + 3.0).unwrap().unwrap();
            for b in [0.0, 0.5, 3.0, 30.0, 300.0] {
                let (y, _) = attached(d, b, boundary, boundary.twist_square - b).unwrap();
                let residual = (y * y - d).powi(2) + (b * y) * y - 16.0 * y;
                assert!(
                    residual.abs() <= 2e-13 * 16.0 * y,
                    "d={d}, b={b}, r={residual}"
                );
            }
            let fraction = 0.99999;
            let b = boundary.twist_square * fraction * fraction;
            let distance = boundary.twist_square * (1.0 - fraction) * (1.0 + fraction);
            let (y, _) = attached(d, b, boundary, distance).unwrap();
            let residual = (y * y - d).powi(2) + (b * y) * y - 16.0 * y;
            assert!(
                residual.abs() <= 2e-13 * 16.0 * y,
                "near limit: d={d}, r={residual}"
            );
        }
    }
}
