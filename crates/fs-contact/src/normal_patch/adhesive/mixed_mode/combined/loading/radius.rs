//! First boundary of F=(y²-d)²+lambda²(b*y²+t)-ell*y.
//!
//! H=A'B-AB' is strictly convex for positive radius:
//! H''=40*b*y³+24*t*y+2*ell*b>0. With t>0, H(0)<0, so there
//! is exactly one positive root. With t=0, divide out b*y and solve the
//! strictly increasing 2*y⁴+ell*y-2*d² instead. The first boundary is the
//! opening edge when d>0 and b*d>=t; otherwise it is that radial fold.
//! These bounded Dd sign/root screens are not interval certificates.

use fs_math::dd::Dd;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Kind {
    Fold,
    Opening,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct Boundary {
    pub kind: Kind,
    pub radius: f64,
    pub bracket: [f64; 2],
}

#[derive(Clone, Debug, PartialEq)]
pub(super) enum Error<E> {
    Probe(E),
    Unresolved,
    IterationBudget,
}

/// `upper` is the original attached radius at lambda=0. Every actual
/// polynomial probe polls the caller's cancellation/evaluation owner.
#[allow(clippy::too_many_arguments)]
pub(super) fn locate<E>(
    d: Dd,
    b: Dd,
    t: Dd,
    ell: Dd,
    upper: f64,
    maximum_iterations: usize,
    mut probe: impl FnMut() -> Result<(), E>,
) -> Result<Option<Boundary>, Error<E>> {
    if [d, b, t, ell].into_iter().any(|x| !finite(x))
        || b.lt(Dd::ZERO)
        || t.lt(Dd::ZERO)
        || !Dd::ZERO.lt(b + t)
        || !Dd::ZERO.lt(ell)
        || !upper.is_normal()
        || upper <= 0.0
    {
        return Err(Error::Unresolved);
    }
    if d == Dd::ZERO && t == Dd::ZERO {
        return Ok(None);
    }
    probe().map_err(Error::Probe)?;
    let opening_competition = b * d;
    if !finite(opening_competition) {
        return Err(Error::Unresolved);
    }
    if Dd::ZERO.lt(d) && !opening_competition.lt(t) {
        let edge = d.to_f64().sqrt();
        let radius = edge.next_up();
        let opening = dd(radius) * dd(radius) - d;
        if !radius.is_normal() || !finite(opening) || opening.lt(Dd::ZERO) {
            return Err(Error::Unresolved);
        }
        return Ok(Some(Boundary {
            kind: Kind::Opening,
            radius,
            bracket: [edge.next_down(), radius],
        }));
    }
    let value = |y: f64| {
        let y = dd(y);
        let square = y * y;
        if t == Dd::ZERO {
            dd(2.0) * square * square + ell * y - dd(2.0) * d * d
        } else {
            t * (dd(4.0) * y * (square - d) - ell)
                + b * y * (dd(2.0) * (square * square - d * d) + ell * y)
        }
    };
    let mut lo = 0.0_f64;
    let mut hi = upper * (1.0 + 32.0 * f64::EPSILON);
    probe().map_err(Error::Probe)?;
    let at_upper = value(hi);
    if !hi.is_finite() || !finite(at_upper) || !Dd::ZERO.lt(at_upper) {
        return Err(Error::Unresolved);
    }
    for _ in 0..maximum_iterations {
        let middle = lo.midpoint(hi);
        if lo.next_up() >= hi || middle == lo || middle == hi {
            return Ok(Some(Boundary {
                kind: Kind::Fold,
                radius: hi,
                bracket: [lo, hi],
            }));
        }
        probe().map_err(Error::Probe)?;
        let f = value(middle);
        if !finite(f) {
            return Err(Error::Unresolved);
        }
        if f == Dd::ZERO {
            return Ok(Some(Boundary {
                kind: Kind::Fold,
                radius: middle,
                bracket: [middle; 2],
            }));
        }
        if f.lt(Dd::ZERO) {
            lo = middle;
        } else {
            hi = middle;
        }
    }
    Err(Error::IterationBudget)
}

fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}
fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn boundary(d: f64, b: f64, t: f64) -> Option<Boundary> {
        locate(dd(d), dd(b), dd(t), dd(16.0), 4.0, 128, || Ok::<_, ()>(())).unwrap()
    }

    #[test]
    fn manufactured_mixed_fold_and_original_stationarity_agree() {
        // F(1,lambda=1)=F_y(1,lambda=1)=0 for d=-1,b=4,t=8,ell=16.
        let result = boundary(-1.0, 4.0, 8.0).unwrap();
        assert_eq!(result.kind, Kind::Fold);
        assert!(result.bracket[0] <= 1.0 && result.bracket[1] >= 1.0);
        assert!((result.radius - 1.0).abs() <= 4.0 * f64::EPSILON);
        let y = result.radius;
        let load_square = (16.0 * y - (y * y + 1.0).powi(2)) / (4.0 * y * y + 8.0);
        assert!((load_square - 1.0).abs() <= 4.0 * f64::EPSILON);
        assert!((4.0 * y * (y * y + 1.0) + 8.0 * y * load_square - 16.0).abs() < 1e-13);
    }

    #[test]
    fn opening_competes_with_fold_and_zero_approach_twist_has_no_finite_limit() {
        assert_eq!(boundary(1.0, 8.0, 2.0).unwrap().kind, Kind::Opening);
        assert_eq!(boundary(1.0, 2.0, 8.0).unwrap().kind, Kind::Fold);
        assert_eq!(boundary(1.0, 4.0, 4.0).unwrap().kind, Kind::Opening);
        assert!(boundary(0.0, 4.0, 0.0).is_none());
        assert_eq!(boundary(0.0, 4.0, 1.0).unwrap().kind, Kind::Fold);
        let pure_twist = boundary(-1.0, 4.0, 0.0).unwrap();
        let y = pure_twist.radius;
        assert!((2.0 * y.powi(4) + 16.0 * y - 2.0).abs() < 1e-13);
        let pure_shear = boundary(0.0, 0.0, 4.0).unwrap();
        assert!((pure_shear.radius - 4.0_f64.cbrt()).abs() < 1e-14);
    }

    #[test]
    fn probes_and_iterations_remain_bounded() {
        let stopped = locate(dd(-1.0), dd(4.0), dd(8.0), dd(16.0), 4.0, 128, || Err(7));
        assert_eq!(stopped.unwrap_err(), Error::Probe(7));
        let budget = locate(dd(-1.0), dd(4.0), dd(8.0), dd(16.0), 4.0, 1, || {
            Ok::<_, ()>(())
        });
        assert_eq!(budget.unwrap_err(), Error::IterationBudget);
    }
}
