//! Certified roots of the absorbing uniaxial extraordinary dispersion.
//!
//! The extraordinary normal wavenumber z solves F(z) = f z^2 + 2 g p z + (h p2 - c) = 0 with, for N = sum a_i^2,
//! o2 = n_o^2, e2 = n_e^2:
//! f = o2 (a_x^2 + a_y^2) + e2 a_z^2, g = (e2 - o2) a_x a_z, h = o2 (a_y^2 + a_z^2) + e2 a_x^2, c = N o2 e2.
//! (The module header's dispersion divided by N.) The double-precision roots of the stable pair have an
//! absolute error of about 1e-16 |z| that the inputs do not carry: with the optic axis in the interface plane and
//! normal incidence f = o2 and c = o2 e2, so the ordinary tensor element is multiplied in and divided out, and a
//! loss of 1e-24 on the extraordinary axis, or none at all, comes out as Im z = -5e-17. A dichroic polarizer at normal
//! incidence was refused as "ambiguous" for that reason.
//!
//! Each root is therefore polished by one Newton step whose residual is evaluated in complex double-double
//! arithmetic from the exact float inputs, and returned with an a-posteriori bound. Where that bound separates Im z
//! from zero the sign decides the outgoing branch; where it does not, the normal flux of the mode decides, with its
//! own bound; where neither does, the roots are refused. The polished root is the one returned.
//!
//! Proof of the certificate. F is quadratic, so Taylor's series at z0 ends at its second term:
//! 0 = F(z*) = F(z0) + F'(z0) (z* - z0) + f (z* - z0)^2. Dividing by F'(z0) gives, exactly,
//! (z0 - F(z0)/F'(z0)) - z* = f (z0 - z*)^2 / F'(z0). With Delta = F(z0)/F'(z0) and kappa = |f| / |F'(z0)|, the
//! Newton-Kantorovich theorem puts one root z* within 2 |Delta| of z0 whenever h = kappa |Delta| <= 1/4, so the
//! Newton iterate is within 4 kappa |Delta|^2 of it.

use super::*;

/// 2^-53.
const U: f64 = f64::EPSILON * 0.5;
/// 2^-98: the bound on one double-double evaluation of F or F' relative to the sum of its absolute terms.
/// Each real double-double addition, subtraction and multiplication is within 2^-104 of its exact result (fs-math
/// dd.rs). Counting the operations on each term of F(z0) (normwise, each real operation contributing at most
/// 2^-104 of the magnitudes it combines):
///  - f = o2 s_xy + e2 a_z^2: o2 and e2 one subtraction each (their squares and 2 n_r n_i are exact products), s_xy
///    one addition, two scalings, one addition: at most 8 on |o2| s_xy + |e2| a_z^2; f z0^2: z0^2 two operations, the
///    complex product four: 14 in all on f's term;
///  - 2 g p z0: g one subtraction and a scaling on (|o2| + |e2|) |a_x a_z|, 9, the scaling by p 1, the complex product
///    4: 14;
///  - h p2: h as f, 8, and one scaling: 9;
///  - c = N o2 e2: o2 and e2 1 each, the complex product 6, N two additions, the scaling 1: 11;
///  - three sums: 3.
/// The total is at most 51 < 64 = 2^6, so 2^6 * 2^-104 = 2^-98. F'(z0) = 2 (f z0 + g p) has fewer operations.
const ETA: f64 = 3.155_443_620_884_047e-30;
/// 2^-104.
const DD_UNIT: f64 = 4.930_380_657_631_324e-32;

/// Complex double-double. Every real operation is within 2^-104 of its exact result.
#[derive(Clone, Copy, Debug)]
pub(super) struct Cdd {
    re: Dd,
    im: Dd,
}

impl Cdd {
    const ZERO: Cdd = Cdd {
        re: Dd::ZERO,
        im: Dd::ZERO,
    };
    fn new(re: Dd, im: Dd) -> Cdd {
        Cdd { re, im }
    }
    fn lift(z: C64) -> Cdd {
        Cdd::new(Dd::from_f64(z.re), Dd::from_f64(z.im))
    }
    fn real(re: Dd) -> Cdd {
        Cdd::new(re, Dd::ZERO)
    }
    fn add(self, o: Cdd) -> Cdd {
        Cdd::new(self.re + o.re, self.im + o.im)
    }
    fn sub(self, o: Cdd) -> Cdd {
        Cdd::new(self.re - o.re, self.im - o.im)
    }
    fn mul(self, o: Cdd) -> Cdd {
        Cdd::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
    fn scale(self, k: Dd) -> Cdd {
        Cdd::new(self.re * k, self.im * k)
    }
    fn conj(self) -> Cdd {
        Cdd::new(self.re, -self.im)
    }
    fn norm_sq(self) -> Dd {
        self.re * self.re + self.im * self.im
    }
    fn recip(self) -> Cdd {
        let norm = self.norm_sq();
        Cdd::new(self.re / norm, -(self.im / norm))
    }
    fn round(self) -> C64 {
        C64::new(self.re.to_f64(), self.im.to_f64())
    }
    fn magnitude(self) -> f64 {
        self.round().abs()
    }
}

fn cross(x: [Cdd; 3], y: [Cdd; 3]) -> [Cdd; 3] {
    [
        x[1].mul(y[2]).sub(x[2].mul(y[1])),
        x[2].mul(y[0]).sub(x[0].mul(y[2])),
        x[0].mul(y[1]).sub(x[1].mul(y[0])),
    ]
}

fn cross_magnitude(x: [f64; 3], y: [f64; 3]) -> [f64; 3] {
    [
        x[1] * y[2] + x[2] * y[1],
        x[2] * y[0] + x[0] * y[2],
        x[0] * y[1] + x[1] * y[0],
    ]
}

/// A polished root and the bound on its distance from the exact root.
#[derive(Clone, Copy, Debug)]
pub(super) struct Root {
    pub(super) z: C64,
    /// Bound on |z1 - z*| before the final rounding of the two components to f64.
    pub(super) error: f64,
}

impl Root {
    /// Bound on |Im z - Im z*|: the double-double error plus the rounding of the imaginary component itself.
    pub(super) fn im_bound(&self) -> f64 {
        self.error + U * self.z.im.abs()
    }
}

/// Which way a certified root carries power.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Direction {
    Outgoing,
    Incoming,
}

pub(super) const NOT_CERTIFIED: &str =
    "absorbing extraordinary roots not separated: no certified branch";
const AMBIGUOUS: &str = "absorbing extraordinary outgoing branch is ambiguous";
const NO_BRANCH: &str = "absorbing extraordinary mode has no resolved outgoing branch";

/// The branch of one root: by the sign of Im z where the bound separates it from zero, else by the sign of the normal
/// flux where its bound separates that from zero, else refused.
pub(super) fn classify(
    im: f64,
    im_bound: f64,
    flux: impl FnOnce() -> (f64, f64),
) -> Result<Direction> {
    if im < -im_bound {
        return Ok(Direction::Outgoing);
    }
    if im > im_bound {
        return Ok(Direction::Incoming);
    }
    let (flux, bound) = flux();
    if flux > bound {
        Ok(Direction::Outgoing)
    } else if flux < -bound {
        Ok(Direction::Incoming)
    } else {
        Err(NOT_CERTIFIED)
    }
}

/// The one outgoing root of two.
pub(super) fn select(directions: [Direction; 2]) -> Result<usize> {
    match directions.map(|d| d == Direction::Outgoing) {
        [true, false] => Ok(0),
        [false, true] => Ok(1),
        [true, true] => Err(AMBIGUOUS),
        [false, false] => Err(NO_BRANCH),
    }
}

/// The N-scaled dispersion function in complex double-double, from the exact float inputs.
pub(super) struct Dispersion {
    o2: Cdd,
    e2: Cdd,
    axis: [f64; 3],
    norm: Dd,
    f: Cdd,
    g: Cdd,
    h: Cdd,
    c: Cdd,
    p: f64,
    p2: Dd,
    // Sums of absolute terms: |o2| s_xy + |e2| a_z^2, (|o2| + |e2|) |a_x a_z|, |o2| s_yz + |e2| a_x^2, N |o2| |e2|.
    sum_f: f64,
    sum_g: f64,
    sum_h: f64,
    sum_c: f64,
}

impl Dispersion {
    pub(super) fn new(indices: [C64; 2], axis: [f64; 3], p: f64, p2: Dd) -> Dispersion {
        let square = |n: C64| {
            let (re, im) = (Dd::from_f64(n.re), Dd::from_f64(n.im));
            Cdd::new(re * re - im * im, re * im * Dd::from_f64(2.0))
        };
        let o2 = square(indices[0]);
        let e2 = square(indices[1]);
        let a2 = axis.map(|v| Dd::from_f64(v) * Dd::from_f64(v));
        let axz = Dd::from_f64(axis[0]) * Dd::from_f64(axis[2]);
        let s_xy = a2[0] + a2[1];
        let s_yz = a2[1] + a2[2];
        let norm = a2[0] + a2[1] + a2[2];
        let f = o2.scale(s_xy).add(e2.scale(a2[2]));
        let g = e2.sub(o2).scale(axz);
        let h = o2.scale(s_yz).add(e2.scale(a2[0]));
        let c = o2.mul(e2).scale(norm);
        let (mo, me) = (o2.magnitude(), e2.magnitude());
        Dispersion {
            o2,
            e2,
            axis,
            norm,
            f,
            g,
            h,
            c,
            p,
            p2,
            sum_f: mo * s_xy.to_f64() + me * a2[2].to_f64(),
            sum_g: (mo + me) * axz.to_f64().abs(),
            sum_h: mo * s_yz.to_f64() + me * a2[0].to_f64(),
            sum_c: norm.to_f64() * mo * me,
        }
    }

    /// F(z) and F'(z), each within `ETA` of the sum of its absolute terms.
    fn evaluate(&self, z: Cdd) -> (Cdd, Cdd) {
        let p = Dd::from_f64(self.p);
        let two = Dd::from_f64(2.0);
        let g_p = self.g.scale(p);
        let residual = self
            .f
            .mul(z.mul(z))
            .add(g_p.mul(z).scale(two))
            .add(self.h.scale(self.p2))
            .sub(self.c);
        let slope = self.f.mul(z).add(g_p).scale(two);
        (residual, slope)
    }

    /// One certified Newton step from `z0`, or `None` where the certificate does not hold (the slope is not
    /// separated from its own error, or 4 kappa |Delta| > 1: the two roots are closer than about 4 |Delta|).
    pub(super) fn refine(&self, z0: C64) -> Option<Root> {
        if !finite(z0) {
            return None;
        }
        let z = Cdd::lift(z0);
        let (residual, slope) = self.evaluate(z);
        let radius = z0.abs();
        let p_abs = self.p.abs();
        let sum_value = self.sum_f * radius * radius
            + 2.0 * self.sum_g * p_abs * radius
            + self.sum_h * self.p2.to_f64()
            + self.sum_c;
        let sum_slope = 2.0 * (self.sum_f * radius + self.sum_g * p_abs);
        let slope_magnitude = slope.magnitude();
        let slope_error = ETA * sum_slope;
        // The exact slope is at least this: its f64 magnitude, the double-double error and the magnitude's own rounding.
        let slope_low = slope_magnitude * (1.0 - 64.0 * U) - slope_error;
        if !(slope_low > 0.0) || !(slope_error <= 0.25 * slope_magnitude) {
            return None;
        }
        let numerator = residual.mul(slope.conj());
        let norm_sq = slope.norm_sq();
        let delta = Cdd::new(numerator.re / norm_sq, numerator.im / norm_sq);
        let delta_magnitude = delta.magnitude();
        // |F/F' - delta|: the residual's error over the slope, the slope's error acting on delta, and the division.
        let delta_error =
            (ETA * sum_value + delta_magnitude * slope_error) / slope_low + ETA * delta_magnitude;
        let step = delta_magnitude + delta_error;
        let kappa = self.f.magnitude() * (1.0 + 64.0 * U) / slope_low;
        if !(4.0 * kappa * step <= 1.0) {
            return None;
        }
        let z1 = z.sub(delta);
        let rounded = z1.round();
        let error =
            (4.0 * kappa * step * step + delta_error + DD_UNIT * rounded.abs()) * (1.0 + 64.0 * U);
        if !finite(rounded) || !error.is_finite() {
            return None;
        }
        Some(Root { z: rounded, error })
    }

    /// The normal flux of the extraordinary mode at `root`, evaluated in double-double from the construction of its
    /// fields in `raw_eigenvectors`, with a bound that separates its sign from zero.
    ///
    /// The fields are q = (p, 0, z), h = (0, 1, 0) for an optic axis in the incidence plane and q x a otherwise,
    /// v = -q x h, e_j = v_j / o2 + ((v . a) / N) (1 / e2 - 1 / o2) a_j, flux = (Re(e_x h_y*) - Re(e_y h_x*)) / 2.
    /// The positive normalisation of the mode does not change its sign. At most 24 successive double-double operations
    /// reach any product (the reciprocals 4 each, their difference and the scaling 2, the projection 5, its product 3,
    /// the sum 1, the cross products 3 and the two products and difference 3 more), each within 2^-104 of the
    /// magnitudes it combines; 24 * 2^-104 < 2^-99. The sum of absolute partial products S bounds each of them, and
    /// the root's own certified error enters through a polynomial of degree four in z: 5 (E / |z| + u) S. So the
    /// separating bound is (2^-99 + 5 (E / |z| + u)) S.
    pub(super) fn flux(&self, root: &Root) -> (f64, f64) {
        let z = root.z;
        let p = Dd::from_f64(self.p);
        let q = [Cdd::real(p), Cdd::ZERO, Cdd::lift(z)];
        let q_magnitude = [self.p.abs(), 0.0, z.abs()];
        let a = self.axis.map(|v| Cdd::real(Dd::from_f64(v)));
        let a_magnitude = self.axis.map(f64::abs);
        let (h, h_magnitude) = if self.axis[1] == 0.0 {
            ([Cdd::ZERO, Cdd::real(Dd::ONE), Cdd::ZERO], [0.0, 1.0, 0.0])
        } else {
            (cross(q, a), cross_magnitude(q_magnitude, a_magnitude))
        };
        let negate = |c: Cdd| Cdd::ZERO.sub(c);
        let v = cross(q, h).map(negate);
        let v_magnitude = cross_magnitude(q_magnitude, h_magnitude);
        let w = self.axis.map(|x| Dd::from_f64(x) / self.norm);
        let projection = (0..3).fold(Cdd::ZERO, |s, j| s.add(v[j].scale(w[j])));
        let projection_magnitude: f64 = (0..3).map(|j| v_magnitude[j] * w[j].to_f64().abs()).sum();
        let r_o = self.o2.recip();
        let r_e = self.e2.recip();
        let difference = r_e.sub(r_o);
        let e: [Cdd; 3] = core::array::from_fn(|j| {
            v[j].mul(r_o)
                .add(projection.mul(difference.scale(Dd::from_f64(self.axis[j]))))
        });
        let e_magnitude: [f64; 3] = core::array::from_fn(|j| {
            v_magnitude[j] * r_o.magnitude()
                + projection_magnitude * (r_e.magnitude() + r_o.magnitude()) * a_magnitude[j]
        });
        let product = |j: usize, k: usize| e[j].re * h[k].re + e[j].im * h[k].im;
        let flux = ((product(0, 1) - product(1, 0)) * Dd::from_f64(0.5)).to_f64();
        let sum = 0.5 * (e_magnitude[0] * h_magnitude[1] + e_magnitude[1] * h_magnitude[0]);
        let relative = 2.0_f64.powi(-99) + 5.0 * (root.error / z.abs().max(f64::MIN_POSITIVE) + U);
        (flux, relative * sum)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Scaled units: the indices of the 1.45 - 0.12i and 2.1 - ik crystal divided by the binary scale 2.
    fn indices(extraordinary_loss: f64) -> [C64; 2] {
        [C64::new(0.725, -0.06), C64::new(1.05, extraordinary_loss)]
    }

    // The two roots of the double-precision quadratic: only a start for the certified step.
    fn starts(d: &Dispersion) -> [C64; 2] {
        let (f, g, h, c) = (d.f.round(), d.g.round(), d.h.round(), d.c.round());
        let gp = g.scale(d.p);
        let discriminant = gp * gp - f * (h.scale(d.p2.to_f64()) - c);
        let root = discriminant.sqrt();
        [(-gp + root) / f, (-gp - root) / f]
    }

    // The inverse-tensor form of the dispersion, F / (N o2 e2), in double-double: the reduction identity.
    fn inverse_form(d: &Dispersion, z: Cdd) -> Cdd {
        let a2 = d.axis.map(|v| Dd::from_f64(v) * Dd::from_f64(v));
        let axz = Dd::from_f64(d.axis[0]) * Dd::from_f64(d.axis[2]);
        let (r_o, r_e) = (d.o2.recip(), d.e2.recip());
        let a = r_e
            .scale((a2[0] + a2[1]) / d.norm)
            .add(r_o.scale(a2[2] / d.norm));
        let b = r_o.sub(r_e).scale(Dd::from_f64(d.p) * axz / d.norm);
        let c = r_e
            .scale((a2[1] + a2[2]) / d.norm)
            .add(r_o.scale(a2[0] / d.norm))
            .scale(d.p2)
            .sub(Cdd::real(Dd::ONE));
        a.mul(z.mul(z))
            .add(b.mul(z).scale(Dd::from_f64(2.0)))
            .add(c)
    }

    #[test]
    fn g0_the_scaled_polynomial_is_the_inverse_tensor_form_times_n_o2_e2() {
        let axes = [
            [1.0, 0.0, 0.0],
            [0.6, 0.0, 0.8],
            [0.3, 0.5, 0.812_403_840_463_596_1],
            [1.0, 0.0, 1.0e-9],
        ];
        for axis in axes {
            for p in [0.0, 0.15, 0.35] {
                for loss in [0.0, -5.0e-25, -5.0e-4] {
                    let pd = Dd::from_f64(p);
                    let d = Dispersion::new(indices(loss), axis, p, pd * pd);
                    for start in starts(&d) {
                        let root = d.refine(start).expect("certified");
                        // Polishing moves the root by the start's rounding only, and the bound is small.
                        assert!(
                            (root.z - start).abs() <= 1.0e-13,
                            "{axis:?} {p} {loss} {start:?} {root:?}"
                        );
                        assert!(root.error < 1.0e-24, "{axis:?} {p} {loss} {root:?}");
                        let z = Cdd::lift(root.z);
                        let (value, _) = d.evaluate(z);
                        let reduced = inverse_form(&d, z).mul(d.o2.mul(d.e2).scale(d.norm));
                        let difference = value.sub(reduced).magnitude();
                        let radius = root.z.abs();
                        let sum = d.sum_f * radius * radius
                            + 2.0 * d.sum_g * p.abs() * radius
                            + d.sum_h * p * p
                            + d.sum_c;
                        assert!(
                            difference <= 2.0_f64.powi(-85) * sum,
                            "{axis:?} {p} {loss} {difference:e} {sum:e}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn g1_a_loss_of_1e_24_survives_where_the_start_root_had_lost_it() {
        let exact = C64::new(2.1 / 2.0, -1.0e-24 / 2.0);
        let d = Dispersion::new(indices(exact.im), [1.0, 0.0, 0.0], 0.0, Dd::ZERO);
        // The roots of the double-precision quadratic carry the Im of about -5e-17 of the 27209 print.
        for start in [C64::new(1.05, -5.17e-17), C64::new(-1.05, -4.93e-17)] {
            let root = d.refine(start).expect("certified");
            let sign = if start.re > 0.0 { 1.0 } else { -1.0 };
            let expected = C64::new(sign * exact.re, sign * exact.im);
            assert!(
                (root.z.im - expected.im).abs() <= root.im_bound(),
                "{root:?} {expected:?}"
            );
            assert!(root.im_bound() < 1.0e-28, "{root:?}");
            assert!((root.z.re - expected.re).abs() <= root.error + U * root.z.re.abs());
        }
    }

    #[test]
    fn g1_an_exactly_transparent_extraordinary_axis_leaves_im_inside_its_bound() {
        let d = Dispersion::new(indices(0.0), [1.0, 0.0, 0.0], 0.0, Dd::ZERO);
        for start in [C64::new(1.05, -5.17e-17), C64::new(-1.05, -4.93e-17)] {
            let root = d.refine(start).expect("certified");
            assert!(root.z.im.abs() <= root.im_bound(), "{root:?}");
            // So the sign of Im decides nothing and the flux separates the two roots.
            let (flux, bound) = d.flux(&root);
            assert!(flux.abs() > bound && bound < 1.0e-14, "{flux} {bound}");
            assert_eq!(flux > 0.0, start.re > 0.0, "{flux}");
        }
    }

    #[test]
    fn g0_a_start_the_step_cannot_certify_is_refused() {
        let d = Dispersion::new(indices(-5.0e-4), [1.0, 0.0, 0.0], 0.0, Dd::ZERO);
        // F'(0) = 0 at normal incidence, and 4 kappa |Delta| = 26 at 0.2.
        assert!(d.refine(C64::ZERO).is_none());
        assert!(d.refine(C64::new(0.2, 0.0)).is_none());
        assert!(d.refine(C64::new(f64::NAN, 0.0)).is_none());
    }

    #[test]
    fn g0_classification_follows_the_bound_then_the_flux_and_refuses_where_neither_separates() {
        let never = || -> (f64, f64) { panic!("the flux is not consulted where Im is certified") };
        assert_eq!(classify(-1.0e-24, 4.0e-30, never), Ok(Direction::Outgoing));
        assert_eq!(classify(1.0e-24, 4.0e-30, never), Ok(Direction::Incoming));
        assert_eq!(
            classify(0.0, 4.0e-30, || (0.4, 1.0e-15)),
            Ok(Direction::Outgoing)
        );
        assert_eq!(
            classify(3.0e-31, 4.0e-30, || (-0.4, 1.0e-15)),
            Ok(Direction::Incoming)
        );
        // Planted negatives: a flux inside its own bound, and an Im inside its bound with no flux to decide.
        assert_eq!(
            classify(0.0, 4.0e-30, || (1.0e-20, 1.0e-15)),
            Err(NOT_CERTIFIED)
        );
        assert_eq!(
            classify(-4.0e-30, 4.0e-30, || (0.0, 0.0)),
            Err(NOT_CERTIFIED)
        );
        // Two roots: one outgoing selects it; equal-sign fluxes, or two certified outgoing roots, refuse.
        let [i, o] = [Direction::Incoming, Direction::Outgoing];
        assert_eq!(select([i, o]), Ok(1));
        assert_eq!(select([o, i]), Ok(0));
        assert_eq!(select([o, o]), Err(AMBIGUOUS));
        assert_eq!(select([i, i]), Err(NO_BRANCH));
        let positive = classify(0.0, 4.0e-30, || (0.4, 1.0e-15)).unwrap();
        assert_eq!(select([positive, positive]), Err(AMBIGUOUS));
        let negative = classify(0.0, 4.0e-30, || (-0.4, 1.0e-15)).unwrap();
        assert_eq!(select([negative, negative]), Err(NO_BRANCH));
    }
}
