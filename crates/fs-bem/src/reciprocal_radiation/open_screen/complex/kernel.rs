//! Entire scalar factors of the declared cosine/antipodal screen continuation.
use core::f64::consts::PI;
use fs_math::{c64::C64, det};

/// (cos(k*r)-1)/(4*pi*r) and its analytic k derivative -sin(k*r)/(4*pi).
/// The coincident value and derivative are both zero. The series retains k²r
/// without forming a tiny (kr)² and then dividing it by a tiny r.
pub(crate) fn cosine_remainder(k: C64, r: f64) -> (C64, C64) {
    if r == 0.0 {
        return (C64::ZERO, C64::ZERO);
    }
    let z = k.scale(r);
    let scale = 1.0 / (4.0 * PI);
    if z.abs() < 0.5 {
        let minus_z2 = -(z * z);
        let mut cosine_term = (k * z).scale(-0.5);
        let mut cosine = cosine_term;
        let mut sine_term = z;
        let mut sine = sine_term;
        for n in 2..=16 {
            let twice = f64::from(2 * n);
            cosine_term = (cosine_term * minus_z2).scale(1.0 / ((twice - 1.0) * twice));
            sine_term = (sine_term * minus_z2).scale(1.0 / ((twice - 2.0) * (twice - 1.0)));
            cosine = cosine + cosine_term;
            sine = sine + sine_term;
        }
        return (cosine.scale(scale), sine.scale(-scale));
    }
    let (sin_x, cos_x) = det::sin_cos(z.re);
    // expm1 keeps a small imaginary component when Re(k*r) is not small.
    let ep = det::expm1(0.5 * z.im);
    let em = det::expm1(-0.5 * z.im);
    let sinh_half = 0.5 * (ep - em);
    let cosh_half = 1.0 + 0.5 * (ep + em);
    let cosh_minus_one = 2.0 * sinh_half * sinh_half;
    let cosh = 1.0 + cosh_minus_one;
    let sinh = 2.0 * sinh_half * cosh_half;
    let sin_half = det::sin(0.5 * z.re);
    let cosine_minus_one = C64::new(
        (-2.0 * sin_half * sin_half) * cosh + cosh_minus_one,
        -sin_x * sinh,
    );
    let derivative = C64::new(-sin_x * cosh, -cos_x * sinh).scale(scale);
    (cosine_minus_one.scale(scale / r), derivative)
}

/// c_+=-ik(n.d)exp(-ik*u), c_-=+ik(n.d)exp(+ik*u), and d/dk of each.
/// Here u=d.(x-origin) for one FIXED geometric origin. No conjugation appears
/// in this analytic extension; both signs are evaluated at the same complex k.
pub(crate) fn plane_factors(k: C64, u: f64, normal_dot: f64) -> [C64; 4] {
    let (sine, cosine) = det::sin_cos(k.re * u);
    let forward_phase = C64::new(cosine, -sine).scale(det::exp(k.im * u));
    let backward_phase = C64::new(cosine, sine).scale(det::exp(-k.im * u));
    let forward_base = C64::new(0.0, -normal_dot) * forward_phase;
    let backward_base = C64::new(0.0, normal_dot) * backward_phase;
    let forward = k * forward_base;
    let backward = k * backward_base;
    [
        forward,
        backward,
        forward_base + forward * C64::new(0.0, -u),
        backward_base + backward * C64::new(0.0, u),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: C64, b: C64, tolerance: f64) {
        assert!(
            (a - b).abs() <= tolerance * (1.0 + a.abs() + b.abs()),
            "{a:?} != {b:?}"
        );
    }

    #[test]
    fn cosine_kernel_has_independent_complex_direction_derivatives_and_zero_limit() {
        for k in [
            C64::ZERO,
            C64::new(0.3, -0.2),
            C64::new(2.1, -1.3),
            C64::new(-1.4, 0.8),
        ] {
            for r in [0.0, 0.01, 0.7, 1.6] {
                let (_, derivative) = cosine_remainder(k, r);
                for step in [C64::from_re(2e-6), C64::new(0.0, 2e-6)] {
                    let estimate = (cosine_remainder(k + step, r).0
                        - cosine_remainder(k - step, r).0)
                        / (step.scale(2.0));
                    close(derivative, estimate, 3e-10);
                }
                let (mirror, d_mirror) = cosine_remainder(-k.conj(), r);
                let (value, d) = cosine_remainder(k, r);
                close(mirror, value.conj(), 2e-15);
                close(d_mirror, -d.conj(), 2e-15);
            }
        }
        assert_eq!(cosine_remainder(C64::ZERO, 2.0), (C64::ZERO, C64::ZERO));
        let k = C64::new(1e-12, -2e-12);
        let (value, derivative) = cosine_remainder(k, 0.3);
        let expected = (k * k).scale(-0.3 / (8.0 * PI));
        assert!((value - expected).abs() < 4e-16 * expected.abs());
        assert!((derivative - k.scale(-0.3 / (4.0 * PI))).abs() < 4e-16 * derivative.abs());
    }

    #[test]
    fn antipodal_factors_are_analytic_and_their_product_is_translation_invariant() {
        for k in [
            C64::ZERO,
            C64::from_re(1.1),
            C64::new(1.1, -0.8),
            C64::new(-0.7, 0.6),
        ] {
            let factors = plane_factors(k, 0.3, -0.7);
            for step in [C64::from_re(2e-6), C64::new(0.0, 2e-6)] {
                let left = plane_factors(k - step, 0.3, -0.7);
                let right = plane_factors(k + step, 0.3, -0.7);
                for side in 0..2 {
                    close(
                        factors[side + 2],
                        (right[side] - left[side]) / step.scale(2.0),
                        3e-10,
                    );
                }
            }
            let a = plane_factors(k, -0.2, 0.8);
            let b = plane_factors(k, 0.4, -0.3);
            let shifted_a = plane_factors(k, 1.1, 0.8);
            let shifted_b = plane_factors(k, 1.7, -0.3);
            close(a[0] * b[1], shifted_a[0] * shifted_b[1], 3e-15);
            close(
                a[2] * b[1] + a[0] * b[3],
                shifted_a[2] * shifted_b[1] + shifted_a[0] * shifted_b[3],
                4e-15,
            );
            if k.im == 0.0 {
                close(factors[1], factors[0].conj(), 2e-15);
                close(factors[3], factors[2].conj(), 2e-15);
            }
        }
    }
}
