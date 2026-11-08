//! Stable root of the ORIGINAL combined-load quartic, in scaled coordinates.
//!
//! F(y)=(y²-d)²+b y²+t-ell*y. All four coefficients arrive in Dd directly
//! from original dimensional inputs; a rounded radius scale does not imply
//! ell is exactly 16. The unique positive stationary minimum m satisfies
//! 2m³+(b-2d)m=ell/2. The positive stable-root offset obeys
//! F(m+z)-F(m)=z²*(4m²+ell/(2m)+4mz+z²).
//! Neither a quartic closed form nor the vanishing radial tangent is inverted.
//! Coefficient-side opening admission avoids subtracting rounded radii at the
//! closed-crack boundary. Unresolved arithmetic refuses; these are numerical
//! screens, not interval certificates.

use fs_math::dd::Dd;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Error {
    NoEquilibrium,
    OpeningEdge,
    Unresolved,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct Equilibrium {
    pub radius: f64,
    pub derivative: f64,
    pub marginal: bool,
}

pub(super) fn solve(d: Dd, b: Dd, t: Dd, ell: Dd) -> Result<Equilibrium, Error> {
    if [d, b, t, ell].into_iter().any(|x| !finite(x))
        || b.lt(Dd::ZERO)
        || t.lt(Dd::ZERO)
        || !Dd::ZERO.lt(ell)
    {
        return Err(Error::Unresolved);
    }
    admit_opening(d, b, t, ell)?;
    let c = b - dd(2.0) * d;
    let target = ell / dd(2.0);
    let mut m = dd(stationary(c, target)?);
    // F'' is positive at this minimum even at a contact-radius fold.
    for _ in 0..2 {
        m = m - (m * (dd(2.0) * m * m + c) - target) / (dd(6.0) * m * m + c);
    }
    if !finite(m) || !m.hi.is_normal() || m.hi <= 0.0 {
        return Err(Error::Unresolved);
    }
    let square = m * m;
    let opening = square - d;
    let minimum = opening * opening + b * square + t - ell * m;
    let stationary_residual = m * (dd(2.0) * square + c) - target;
    let scale = (opening * opening + b * square + t + ell * m).to_f64();
    if !finite(minimum) || !finite(stationary_residual) || !scale.is_finite() {
        return Err(Error::Unresolved);
    }
    if minimum == Dd::ZERO && stationary_residual == Dd::ZERO {
        // The dimensional owner must independently authenticate Ua=Uaa=0
        // before publishing a marginal response at the rounded output radius.
        return Ok(Equilibrium {
            radius: m.to_f64(),
            derivative: 0.0,
            marginal: true,
        });
    }
    let uncertainty = roundoff(scale);
    let value = minimum.to_f64();
    if value > uncertainty {
        return Err(Error::NoEquilibrium);
    }
    if value >= -uncertainty {
        return Err(Error::Unresolved);
    }
    let gap = -value;
    let alpha = (dd(4.0) * square + target / m).to_f64();
    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(Error::Unresolved);
    }
    let z = offset(gap, m.to_f64(), alpha)?;
    let radius = (m + dd(z)).to_f64();
    let derivative = 2.0 * z * (alpha + 6.0 * m.to_f64() * z + 2.0 * z * z);
    if !radius.is_normal() || radius <= 0.0 || !derivative.is_finite() || derivative <= 0.0 {
        return Err(Error::Unresolved);
    }
    Ok(Equilibrium {
        radius,
        derivative,
        marginal: false,
    })
}

// At y=sqrt(d), sign(F) is sign((b*d+t)^2-ell^2*d). The edge is
// right of the stationary minimum exactly when 4*b^2*d>ell^2.
// Thus a positive edge F on that side locates the stable root below the
// opening edge even when rounded y*y-d appears positive.
fn admit_opening(d: Dd, b: Dd, t: Dd, ell: Dd) -> Result<(), Error> {
    if !Dd::ZERO.lt(d) {
        return Ok(());
    }
    let load = b * d + t;
    let load_square = load * load;
    let attraction = ell * ell * d;
    let edge = load_square - attraction;
    let edge_tolerance = roundoff((load_square + attraction).to_f64());
    let derivative_left = dd(4.0) * b * b * d;
    let derivative_right = ell * ell;
    let derivative = derivative_left - derivative_right;
    let derivative_tolerance = roundoff((derivative_left + derivative_right).to_f64());
    if !finite(edge)
        || !finite(derivative)
        || !edge_tolerance.is_finite()
        || !derivative_tolerance.is_finite()
    {
        return Err(Error::Unresolved);
    }
    if edge.to_f64() < -edge_tolerance || derivative.to_f64() < -derivative_tolerance {
        return Ok(());
    }
    if derivative.to_f64() > derivative_tolerance && edge.to_f64() > edge_tolerance {
        Err(Error::OpeningEdge)
    } else {
        // Do not turn an unresolved sign into an attached opening-edge state.
        Err(Error::Unresolved)
    }
}

fn stationary(c: Dd, target: Dd) -> Result<f64, Error> {
    let coefficient = c.to_f64();
    let rhs = target.to_f64();
    if !coefficient.is_finite() || !rhs.is_finite() || rhs <= 0.0 {
        return Err(Error::Unresolved);
    }
    let cubic_bound = (0.5 * rhs).cbrt();
    let upper = if coefficient < 0.0 {
        (-0.5 * coefficient).sqrt() + cubic_bound
    } else {
        cubic_bound.min(rhs / coefficient)
    };
    let mut hi = upper * (1.0 + 32.0 * f64::EPSILON);
    let mut lo = 0.0;
    let mut x = hi;
    if !hi.is_finite() || hi <= 0.0 || cubic(hi, c, target)?.0 < 0.0 {
        return Err(Error::Unresolved);
    }
    for _ in 0..64 {
        let (residual, derivative) = cubic(x, c, target)?;
        if residual == 0.0 {
            return Ok(x);
        }
        if residual > 0.0 {
            hi = x;
        } else {
            lo = x;
        }
        let proposal = x - residual / derivative;
        let next = if derivative > 0.0 && proposal > lo && proposal < hi {
            proposal
        } else {
            lo.midpoint(hi)
        };
        if next == x || lo.next_up() >= hi {
            let mut best = x;
            let mut error = residual.abs();
            for candidate in [lo, hi] {
                if candidate > 0.0 {
                    let candidate_error = cubic(candidate, c, target)?.0.abs();
                    if candidate_error < error {
                        best = candidate;
                        error = candidate_error;
                    }
                }
            }
            return Ok(best);
        }
        x = next;
    }
    Err(Error::Unresolved)
}

fn cubic(x: f64, c: Dd, target: Dd) -> Result<(f64, f64), Error> {
    let x = dd(x);
    let value = (x * (dd(2.0) * x * x + c) - target).to_f64();
    let derivative = (dd(6.0) * x * x + c).to_f64();
    if value.is_finite() && derivative.is_finite() {
        Ok((value, derivative))
    } else {
        Err(Error::Unresolved)
    }
}

fn offset(target: f64, m: f64, alpha: f64) -> Result<f64, Error> {
    let upper = (target / alpha).sqrt().min(target.sqrt().sqrt());
    let mut hi = upper * (1.0 + 32.0 * f64::EPSILON);
    let mut lo = 0.0;
    let mut z = hi;
    if !target.is_normal() || !hi.is_finite() || hi <= 0.0 || polynomial(hi, m, alpha)?.0 < target {
        return Err(Error::Unresolved);
    }
    for _ in 0..64 {
        let (value, derivative) = polynomial(z, m, alpha)?;
        let residual = value - target;
        if residual.abs() <= 8.0 * f64::EPSILON * target {
            return Ok(z);
        }
        if residual > 0.0 {
            hi = z;
        } else {
            lo = z;
        }
        let proposal = z - residual / derivative;
        let next = if proposal > lo && proposal < hi {
            proposal
        } else {
            lo.midpoint(hi)
        };
        if next == z {
            break;
        }
        z = next;
    }
    Err(Error::Unresolved)
}

fn polynomial(z: f64, m: f64, alpha: f64) -> Result<(f64, f64), Error> {
    let x = dd(z);
    let value = (x * x * (dd(alpha) + dd(4.0) * dd(m) * x + x * x)).to_f64();
    let derivative = 2.0 * z * (alpha + 6.0 * m * z + 2.0 * z * z);
    if value.is_finite() && derivative.is_finite() {
        Ok((value, derivative))
    } else {
        Err(Error::Unresolved)
    }
}

fn roundoff(scale: f64) -> f64 {
    128.0 * f64::EPSILON * f64::EPSILON * scale
}
fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}
fn finite(x: Dd) -> bool {
    x.hi.is_finite() && x.lo.is_finite()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn retained_original_coefficient_bits_control_fold_and_opening_admission() {
        // Each perturbation is erased by conversion to f64 before the solve.
        // Their original signs nevertheless decide whether an attached root exists.
        let small = 1e-20;
        let exact = solve(dd(-1.0), dd(4.0), dd(8.0), dd(16.0)).unwrap();
        assert!(exact.marginal);
        assert_eq!(exact.radius, 1.0);
        for (d, b, t, ell) in [
            (Dd::from_pair(-1.0, -small), dd(4.0), dd(8.0), dd(16.0)),
            (dd(-1.0), Dd::from_pair(4.0, small), dd(8.0), dd(16.0)),
            (dd(-1.0), dd(4.0), Dd::from_pair(8.0, small), dd(16.0)),
            (dd(-1.0), dd(4.0), dd(8.0), Dd::from_pair(16.0, -small)),
        ] {
            assert_eq!(solve(d, b, t, ell).unwrap_err(), Error::NoEquilibrium);
        }
        for (d, b, t, ell) in [
            (Dd::from_pair(-1.0, small), dd(4.0), dd(8.0), dd(16.0)),
            (dd(-1.0), Dd::from_pair(4.0, -small), dd(8.0), dd(16.0)),
            (dd(-1.0), dd(4.0), Dd::from_pair(8.0, -small), dd(16.0)),
            (dd(-1.0), dd(4.0), dd(8.0), Dd::from_pair(16.0, small)),
        ] {
            let root = solve(d, b, t, ell).unwrap();
            assert!(!root.marginal && root.derivative > 0.0);
        }
        // Binary64 input d=.01 is not exact 1/100. A radius rounded to .1
        // gives y*y-d>0, although the original stable root is on the closed side.
        assert!(0.1_f64 * 0.1 - 0.01 > 0.0);
        assert_eq!(
            solve(dd(0.01), dd(120.0), dd(0.4), dd(16.0)).unwrap_err(),
            Error::OpeningEdge,
        );
    }
}
