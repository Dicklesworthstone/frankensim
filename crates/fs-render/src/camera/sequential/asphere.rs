//! Conic/even-polynomial surface geometry in metres, independent of transport.
//!
//! The conic sag and even powers follow the conventional even-asphere model:
//! https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v261/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Even_Asphere.html
//! Here the quadratic departure is excluded so c remains the vertex curvature.
//! Constructor-time bounds are analytic over the complete clear disk, never
//! inferred from radial samples. Numerical range or unproved ray uniqueness
//! refuses; an exhausted solve must not become a fabricated optical miss.

use fs_math::{dd::Dd, det};

type Hit = ([f64; 3], f64);
type Result<T> = core::result::Result<T, ()>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct AsphericSurface {
    curvature: f64,
    conic: f64,
    coefficients: [f64; 7],
    radius: f64,
    curvature_radius: f64,
    conic_shape: f64,
    // B_j=A_(2j)*R^(2j), j=2..8. Horner evaluation then uses u=r²/R².
    departures: [f64; 7],
    bounds: [f64; 2],
    slope_bound: f64,
}

impl AsphericSurface {
    pub(super) fn try_new(
        curvature: f64,
        conic: f64,
        coefficients: [f64; 7],
        radius: f64,
    ) -> Option<Self> {
        if ![curvature, conic, radius].iter().all(|v| v.is_finite())
            || coefficients.iter().any(|v| !v.is_finite())
            || radius <= 0.0
            || !(radius * radius > 0.0)
            || !(radius * radius).is_finite()
        {
            return None;
        }
        let cr = curvature * radius;
        let shape = if conic == -1.0 {
            0.0
        } else {
            ((conic + 1.0) * cr) * cr
        };
        // Pad product rounding before admitting a real, nonvertical conic.
        let shape_upper = (shape + 64.0 * f64::EPSILON * shape.abs()).next_up();
        if !cr.is_finite() || !shape.is_finite() || shape_upper >= 1.0 {
            return None;
        }
        let cr_bound = (cr.abs() * (1.0 + 64.0 * f64::EPSILON)).next_up();
        let root_lower = det::sqrt((1.0 - shape_upper).next_down()).next_down();
        let base_bound = ((cr_bound / (1.0 + root_lower)) * radius).next_up();
        let mut bounds = if curvature < 0.0 {
            [-base_bound, 0.0]
        } else {
            [0.0, base_bound]
        };
        let mut slope_bound = (cr_bound / root_lower.min(1.0)).next_up();
        let mut departures = [0.0; 7];
        for (j, &coefficient) in coefficients.iter().enumerate() {
            if coefficient == 0.0 {
                continue;
            }
            let order = 2 * (j + 2);
            let mut value = coefficient;
            for _ in 0..order {
                value *= radius;
            }
            if !value.is_finite() || value == 0.0 {
                return None;
            }
            departures[j] = value;
            let bound = (value.abs() * (1.0 + 128.0 * f64::EPSILON)).next_up();
            if coefficient < 0.0 {
                bounds[0] = (bounds[0] - bound).next_down();
            } else {
                bounds[1] = (bounds[1] + bound).next_up();
            }
            slope_bound = (slope_bound + ((order as f64 * bound) / radius).next_up()).next_up();
        }
        if bounds.iter().any(|v| !v.is_finite()) || !slope_bound.is_finite() {
            return None;
        }
        Some(Self {
            curvature,
            conic,
            coefficients,
            radius,
            curvature_radius: cr,
            conic_shape: shape,
            departures,
            bounds,
            slope_bound,
        })
    }

    pub(super) const fn conic_constant(&self) -> f64 {
        self.conic
    }
    pub(super) const fn even_coefficients(&self) -> &[f64; 7] {
        &self.coefficients
    }
    pub(super) const fn sag_bounds(&self) -> [f64; 2] {
        self.bounds
    }

    pub(super) fn sag(&self, radius_squared: f64) -> f64 {
        self.evaluate((radius_squared / (self.radius * self.radius)).clamp(0.0, 1.0))
            .0
    }

    // Returns sag and the factor g with (dz/dx,dz/dy)=g*(x/R,y/R).
    fn evaluate(&self, u: f64) -> (f64, f64) {
        let root = det::sqrt((-self.conic_shape).mul_add(u, 1.0));
        let (mut value, mut derivative) = (0.0_f64, 0.0_f64);
        for (j, &coefficient) in self.departures.iter().enumerate().rev() {
            value = value.mul_add(u, coefficient);
            derivative = derivative.mul_add(u, (j + 2) as f64 * coefficient);
        }
        let sag = (self.curvature_radius * u / (1.0 + root)) * self.radius + u * u * value;
        let gradient = self.curvature_radius / root + (2.0 * u * derivative) / self.radius;
        (sag, gradient)
    }

    pub(super) fn normal(&self, point: [f64; 2]) -> Option<[f64; 3]> {
        if point.iter().any(|v| !v.is_finite()) || det::hypot(point[0], point[1]) > self.radius {
            return None;
        }
        let x = point[0] / self.radius;
        let y = point[1] / self.radius;
        let (_, g) = self.evaluate(x.mul_add(x, y * y).min(1.0));
        let normal = [g * x, g * y, -1.0];
        normal.iter().all(|v| v.is_finite()).then_some(normal)
    }

    /// Positive ray parameter; the transport owner supplies unit directions.
    /// The iteration cap bounds work between its existing surface checkpoints.
    pub(super) fn intersect(
        &self,
        vertex: f64,
        origin: [f64; 3],
        direction: [f64; 3],
    ) -> Result<Option<Hit>> {
        if !vertex.is_finite() || origin.iter().chain(&direction).any(|v| !v.is_finite()) {
            return Err(());
        }
        if direction[2] <= 0.0 {
            return Ok(None);
        }
        if self.departures == [0.0; 7] {
            return self.intersect_conic(vertex, origin, direction);
        }
        // Extend sag constantly beyond the clear disk only for bracketing.
        // This preserves its global Lipschitz bound. If F'=dz-grad(sag).dxy
        // is positive everywhere, the extension has exactly one axial root;
        // a root outside the disk then proves an aperture miss.
        let transverse = det::hypot(direction[0], direction[1]).next_up();
        if (direction[2] - self.slope_bound * transverse).next_down() <= 0.0 {
            return Err(());
        }
        let z_lo = (vertex + self.bounds[0]).next_down();
        let z_hi = (vertex + self.bounds[1]).next_up();
        let mut lo = ((z_lo - origin[2]) / direction[2]).next_down().max(0.0);
        let mut hi = ((z_hi - origin[2]) / direction[2]).next_up();
        if !lo.is_finite() || !hi.is_finite() {
            return Err(());
        }
        if hi <= 0.0 {
            return Ok(None);
        }
        if hi < lo {
            return Err(());
        }
        let evaluate = |t: f64| -> Result<(f64, f64)> {
            let p = point_at(origin, direction, t);
            if p.iter().any(|v| !v.is_finite()) {
                return Err(());
            }
            let radius = det::hypot(p[0], p[1]);
            let (sag, slope) = if radius >= self.radius {
                (self.evaluate(1.0).0, 0.0)
            } else {
                let x = p[0] / self.radius;
                let y = p[1] / self.radius;
                let (sag, g) = self.evaluate(x.mul_add(x, y * y).min(1.0));
                (sag, g * x.mul_add(direction[0], y * direction[1]))
            };
            let residual = (Dd::from_f64(origin[2]) - Dd::from_f64(vertex)
                + Dd::from_f64(direction[2]) * Dd::from_f64(t)
                - Dd::from_f64(sag))
            .to_f64();
            let derivative = direction[2] - slope;
            if !residual.is_finite() || !derivative.is_finite() || derivative <= 0.0 {
                return Err(());
            }
            Ok((residual, derivative))
        };
        let (lower, _) = evaluate(lo)?;
        let (upper, _) = evaluate(hi)?;
        if lo == 0.0 && lower > 0.0 {
            return Ok(None);
        }
        if lower > 0.0 || upper < 0.0 {
            return Err(());
        }
        let mut t = lo + 0.5 * (hi - lo);
        for _ in 0..96 {
            let (residual, derivative) = evaluate(t)?;
            if residual == 0.0 || hi - lo <= 4.0 * f64::EPSILON * hi.abs().max(lo.abs()) {
                return self.admit_hit(origin, direction, t);
            }
            if residual < 0.0 {
                lo = t;
            } else {
                hi = t;
            }
            let newton = t - residual / derivative;
            let next = if newton > lo && newton < hi {
                newton
            } else {
                lo + 0.5 * (hi - lo)
            };
            if next == lo || next == hi {
                return self.admit_hit(origin, direction, next);
            }
            t = next;
        }
        Err(())
    }

    fn intersect_conic(
        &self,
        vertex: f64,
        origin: [f64; 3],
        direction: [f64; 3],
    ) -> Result<Option<Hit>> {
        let c = self.curvature;
        let [ox, oy, oz] = [origin[0], origin[1], origin[2] - vertex];
        let [dx, dy, dz] = direction;
        let mut roots = [f64::INFINITY; 2];
        if c == 0.0 {
            roots[0] = -oz / dz;
        } else {
            // c*(x²+y²+(1+k)*z²)-2*z=0 avoids the sphere centre and admits
            // paraboloids (k=-1), including their genuinely linear ray case.
            let k = 1.0 + self.conic;
            let a = c * dx.mul_add(dx, dy.mul_add(dy, k * dz * dz));
            let b = c * ox.mul_add(dx, oy.mul_add(dy, k * oz * dz)) - dz;
            let constant = c * ox.mul_add(ox, oy.mul_add(oy, k * oz * oz)) - 2.0 * oz;
            if ![a, b, constant].iter().all(|v| v.is_finite()) {
                return Err(());
            }
            if a == 0.0 {
                if b == 0.0 {
                    return if constant == 0.0 { Err(()) } else { Ok(None) };
                }
                roots[0] = -0.5 * (constant / b);
            } else {
                let discriminant = b.mul_add(b, -a * constant);
                if !discriminant.is_finite() {
                    return Err(());
                }
                if discriminant < 0.0 {
                    return Ok(None);
                }
                let q = -b - det::sqrt(discriminant).copysign(b);
                if q == 0.0 {
                    roots[0] = -b / a;
                } else {
                    roots = [q / a, constant / q];
                }
            }
        }
        let mut hit = None;
        for t in roots {
            if let Some((point, t)) = self.admit_hit(origin, direction, t)? {
                // Keep the vertex-connected conic sheet, never its far branch.
                if (1.0 + self.conic) * c * (point[2] - vertex) > 1.0 {
                    continue;
                }
                if hit.is_none_or(|(_, old)| t < old) {
                    hit = Some((point, t));
                }
            }
        }
        Ok(hit)
    }

    fn admit_hit(&self, origin: [f64; 3], direction: [f64; 3], t: f64) -> Result<Option<Hit>> {
        if !t.is_finite() || t <= 0.0 {
            return Ok(None);
        }
        let point = point_at(origin, direction, t);
        if point.iter().any(|v| !v.is_finite()) {
            return Err(());
        }
        if det::hypot(point[0], point[1]) > self.radius {
            return Ok(None);
        }
        Ok(Some((point, t)))
    }
}

fn point_at(origin: [f64; 3], direction: [f64; 3], t: f64) -> [f64; 3] {
    core::array::from_fn(|j| direction[j].mul_add(t, origin[j]))
}

#[cfg(test)]
#[path = "asphere/tests.rs"]
mod tests;
