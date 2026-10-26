//! Elliptical Hertz contact between two elastic bodies.
//!
//! The gap between the bodies near first contact is `curv2` x2^2 + `curv3` x3^2. Each body's
//! half-space compliance may vary with direction in the surface: `h2` in the plane through the
//! indentation normal and x2, `h3` in the plane through the normal and x3, and between them as
//! cos 2 phi (point-load deflection F h / r). For isotropic bodies both are the summed
//! (1 - nu^2)/(pi E). The pressure on the contact ellipse is Hertzian (Willis 1966). With the
//! surface Green's function G(k) = 2 pi H(phi) / |k|, the displacement under
//! p0 sqrt(1 - x2^2/a^2 - x3^2/b^2) inside the ellipse is
//!   u = (pi/4) a b p0 Int_0^2pi H(phi) / c (1 - (x2 cos phi + x3 sin phi)^2 / c^2) dphi,
//! where c = sqrt(a^2 cos^2 phi + b^2 sin^2 phi). Its constant part is the approach and its
//! quadratic parts are the two gap curvatures, so the gap alone fixes F = K delta^(3/2). For
//! isotropic bodies this is Johnson's (Contact Mechanics, 3.4) elliptic-integral solution.

use std::f64::consts::PI;

/// Quadrature points over the circle; resolves contact ellipses with b/a in [1e-2, 1e2].
const QUADRATURE: usize = 4096;

/// The angular moments (Int H/c, Int H cos^2/c^3, Int H sin^2/c^3) at a = 1, b = `aspect`.
fn moments(aspect: f64, mean: f64, half: f64) -> (f64, f64, f64) {
    let step = 2.0 * PI / QUADRATURE as f64;
    let (mut i0, mut i2, mut i3) = (0.0, 0.0, 0.0);
    for j in 0..QUADRATURE {
        let (sin, cos) = ((j as f64 + 0.5) * step).sin_cos();
        let h = mean + half * (cos * cos - sin * sin);
        let c = (cos * cos + aspect * aspect * sin * sin).sqrt();
        i0 += h / c;
        i2 += h * cos * cos / (c * c * c);
        i3 += h * sin * sin / (c * c * c);
    }
    (i0 * step, i2 * step, i3 * step)
}

/// An elastic elliptical Hertz contact: semi-axis a along x2, b along x3.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipticHertz {
    k: f64,
    aspect: f64,
    unit_force: f64,
}

impl EllipticHertz {
    /// The contact of the gap `curv2` x2^2 + `curv3` x3^2 [1/m] on the summed compliance of both
    /// bodies, `h2` in the plane through the normal and x2 and `h3` in the plane through the
    /// normal and x3 (point-load convention F h / r).
    pub fn solve(h2: f64, h3: f64, curv2: f64, curv3: f64) -> Result<Self, String> {
        if ![h2, h3, curv2, curv3]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0)
        {
            return Err(format!(
                "elliptical Hertz contact needs positive compliances and curvatures, got {h2}, {h3}, {curv2}, {curv3}"
            ));
        }
        let (mean, half) = (0.5 * (h2 + h3), 0.5 * (h2 - h3));
        let target = (curv3 / curv2).ln();
        let (mut lo, mut hi) = (-20.0_f64, 20.0_f64);
        // the curvature ratio falls as the ellipse lengthens along x3
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            let (_, i2, i3) = moments(mid.exp(), mean, half);
            if (i3 / i2).ln() > target {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let aspect = (0.5 * (lo + hi)).exp();
        if !(1e-2..=1e2).contains(&aspect) {
            return Err(format!(
                "elliptical Hertz contact b/a {aspect:.3e} is outside the quadrature's [1e-2, 1e2]"
            ));
        }
        let (i0, i2, _) = moments(aspect, mean, half);
        // at a = 1: curv2 = (pi/4) b p0 i2, delta = (pi/4) b p0 i0, F = (2 pi/3) b p0
        let p0 = curv2 / (PI / 4.0 * aspect * i2);
        let delta = PI / 4.0 * aspect * p0 * i0;
        let unit_force = 2.0 * PI / 3.0 * aspect * p0;
        Ok(Self {
            k: unit_force / (delta * delta.sqrt()),
            aspect,
            unit_force,
        })
    }

    /// The contact of two isotropic bodies, (E, nu) each, whose gap has principal radii `r2`
    /// along x2 and `r3` along x3 [m] (a flat against a body of those radii).
    pub fn isotropic(a: (f64, f64), b: (f64, f64), r2: f64, r3: f64) -> Result<Self, String> {
        let h = ((1.0 - a.1 * a.1) / a.0 + (1.0 - b.1 * b.1) / b.0) / PI;
        Self::solve(h, h, 0.5 / r2, 0.5 / r3)
    }

    /// F = K delta^(3/2) [N/m^1.5].
    #[must_use]
    pub const fn k(&self) -> f64 {
        self.k
    }

    /// b/a of the contact ellipse.
    #[must_use]
    pub const fn aspect(&self) -> f64 {
        self.aspect
    }

    /// The semi-axes (a along x2, b along x3) at normal force `force_n` [m]; they scale as
    /// F^(1/3).
    #[must_use]
    pub fn semi_axes_m(&self, force_n: f64) -> (f64, f64) {
        let a = (force_n / self.unit_force).cbrt();
        (a, a * self.aspect)
    }

    /// Mean pressure over the contact ellipse at normal force `force_n` [Pa]; the peak is 3/2
    /// of it.
    #[must_use]
    pub fn mean_pressure_pa(&self, force_n: f64) -> f64 {
        let (a, b) = self.semi_axes_m(force_n);
        force_n / (PI * a * b)
    }

    /// Contact area per unit approach [m]: the ellipse's area pi a b is proportional to delta.
    #[must_use]
    pub fn area_per_approach(&self) -> f64 {
        PI * self.aspect * (self.k / self.unit_force).powf(2.0 / 3.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, rel: f64) {
        assert!(
            (actual - expected).abs() <= rel * expected.abs(),
            "{actual} != {expected} (rel {rel})"
        );
    }

    /// Complete elliptic integrals K(k), E(k) by the arithmetic-geometric mean.
    fn elliptic(k: f64) -> (f64, f64) {
        let (mut a, mut b, mut c) = (1.0_f64, (1.0 - k * k).sqrt(), k);
        let (mut sum, mut power) = (0.5 * c * c, 0.5);
        for _ in 0..40 {
            let next = (0.5 * (a + b), (a * b).sqrt(), 0.5 * (a - b));
            (a, b, c) = next;
            power *= 2.0;
            sum += power * c * c;
        }
        let kk = PI / (2.0 * a);
        (kk, kk * (1.0 - sum))
    }

    #[test]
    fn an_isotropic_ellipse_is_classical_hertz() {
        // Johnson, Contact Mechanics, 3.4: semi-axes a > b, e^2 = 1 - (b/a)^2, gap A x^2 + B y^2.
        let e_star = 50e9;
        let h = 1.0 / (PI * e_star);
        for (a, b, p0) in [
            (1e-3_f64, 0.3e-3, 2e8),
            (1e-3, 0.08e-3, 5e8),
            (2e-3, 2e-3, 1e8),
        ] {
            let e = (1.0 - (b / a) * (b / a)).sqrt();
            let (kk, ee) = elliptic(e);
            let (curv_a, curv_b, delta) = if e > 0.0 {
                let pre = p0 * b / (e_star * e * e * a * a);
                (
                    pre * (kk - ee),
                    pre * ((a / b) * (a / b) * ee - kk),
                    p0 * b * kk / e_star,
                )
            } else {
                // a sphere: p0 = 2 a E*/(pi R), curvature 1/(2R) on each axis
                let r = 2.0 * a * e_star / (PI * p0);
                (0.5 / r, 0.5 / r, a * a / r)
            };
            let force = 2.0 * PI / 3.0 * a * b * p0;
            let contact = EllipticHertz::solve(h, h, curv_a, curv_b).unwrap();
            close(contact.aspect(), b / a, 1e-9);
            close(contact.k(), force / delta.powf(1.5), 1e-9);
            close(contact.mean_pressure_pa(force), 2.0 / 3.0 * p0, 1e-9);
            close(contact.area_per_approach() * delta, PI * a * b, 1e-9);
            let (sa, sb) = contact.semi_axes_m(force);
            close(sa, a, 1e-9);
            close(sb, b, 1e-9);
        }
        // isotropic: a sphere of radius R on a flat is 4/3 E* sqrt(R)
        let (steel, film) = ((212e9, 0.3), (4.9e9, 0.38));
        let e_star = 1.0 / ((1.0 - 0.09) / 212e9 + (1.0 - 0.38 * 0.38) / 4.9e9);
        let sphere = EllipticHertz::isotropic(steel, film, 0.25e-3, 0.25e-3).unwrap();
        close(sphere.k(), 4.0 / 3.0 * e_star * 0.25e-3_f64.sqrt(), 1e-9);
    }
}
