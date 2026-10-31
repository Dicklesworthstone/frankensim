//! Normalized circular-pupil imaging with finite quadratic defocus.
//!
//! For dimensionless radius q=2*pi*a*theta/lambda and rim phase phi,
//! A=2 integral_0^1 r J0(q*r) exp(i*phi*r*r) dr. The density is |A|²/(4*pi).
//! Bauer's expansion of the pupil phase in radial Zernike polynomials gives
//! A=2 exp(i*phi/2) sum_n (2n+1)(-i)^n j_n(phi/2) J_(2n+1)(q)/q.
//! This converges with defocus, independently of the unbounded image radius.
//! A fixed sampled pupil would alias its far field and is not used here.
//!
//! With |phi|<=16, 40 terms omit at most 1.849e-23 in amplitude: the
//! coefficient bound is (|phi|/2)^n/(2n-1)!! and |R_(2n)^0|<=1.
//! This bounds series truncation, not floating-point error. Recurrence and
//! Hankel transition tests below check the latter separately.
//!
//! References: Braat et al., Progress in Optics 51 (2008), Eqs. 2.45/2.61,
//! <https://nijboerzernike.nl/_PDF/BraatProgOpt2008_20140413.pdf>;
//! NIST DLMF 10.17, 10.54.1, 10.60.7 and 10.74(iv).

use fs_math::{det, special::bessel_j1};

const PI: f64 = core::f64::consts::PI;
pub(super) const MAX_PHASE_RAD: f64 = 16.0;
const TERMS: usize = 40;
const MILLER_START: usize = 256;
const HANKEL_SWITCH: f64 = 128.0;

pub(super) fn jinc(x: f64) -> f64 {
    if x.abs() < 1.0e-4 {
        1.0 - x * x / 8.0
    } else {
        2.0 * bessel_j1(x) / x
    }
}

/// Stack-only spectral pupil coefficients. The omitted common phase cannot
/// affect incoherent intensity. Both signs of defocus have the same PSF.
pub(super) struct CircularDefocus {
    phase_rad: f64,
    coefficients: [f64; TERMS],
}

impl CircularDefocus {
    pub(super) fn new(phase_rad: f64) -> Self {
        debug_assert!(phase_rad.is_finite() && phase_rad.abs() <= MAX_PHASE_RAD);
        Self {
            phase_rad,
            coefficients: spherical_coefficients(0.5 * phase_rad.abs()),
        }
    }

    /// Unit-integral density per dimensionless two-dimensional image area.
    /// The caller admits finite, nonnegative radius and the bounded rim phase.
    pub(super) fn density(&self, radius: f64) -> f64 {
        if self.phase_rad == 0.0 {
            let amplitude = jinc(radius);
            return amplitude * amplitude / (4.0 * PI);
        }
        if radius == 0.0 {
            return self.coefficients[0] * self.coefficients[0] / (4.0 * PI);
        }
        let radial = odd_bessel_over_radius(radius);
        let (mut re, mut im, mut cr, mut ci) = (0.0, 0.0, 0.0, 0.0);
        for (n, (&coefficient, bessel)) in self.coefficients.iter().zip(radial).enumerate() {
            let term = 2.0 * coefficient * bessel;
            match n % 4 {
                0 => compensated(&mut re, &mut cr, term),
                1 => compensated(&mut im, &mut ci, -term),
                2 => compensated(&mut re, &mut cr, -term),
                _ => compensated(&mut im, &mut ci, term),
            }
        }
        (re * re + im * im) / (4.0 * PI)
    }
}

fn compensated(sum: &mut f64, correction: &mut f64, value: f64) {
    let y = value - *correction;
    let next = *sum + y;
    *correction = (next - *sum) - y;
    *sum = next;
}

// (2n+1) j_n(x), x in [0,8]. Small x uses its rapidly converging series;
// downward recurrence avoids the unstable high-order forward recurrence.
fn spherical_coefficients(x: f64) -> [f64; TERMS] {
    if x == 0.0 {
        let mut values = [0.0; TERMS];
        values[0] = 1.0;
        return values;
    }
    if x < 0.5 {
        let mut leading = 1.0;
        return std::array::from_fn(|n| {
            if n != 0 {
                leading *= x / (2 * n + 1) as f64;
            }
            let (mut term, mut sum) = (1.0, 1.0);
            for k in 1..=12 {
                term *= -x * x / (2 * k * (2 * n + 2 * k + 1)) as f64;
                sum += term;
            }
            (2 * n + 1) as f64 * leading * sum
        });
    }
    let mut values = [0.0; 66];
    values[64] = 1.0;
    for n in (1..=64).rev() {
        values[n - 1] = (2 * n + 1) as f64 / x * values[n] - values[n + 1];
    }
    let (sin_x, cos_x) = det::sin_cos(x);
    let j0 = sin_x / x;
    let j1 = (j0 - cos_x) / x;
    // The two anchors never vanish together. In particular, do not divide
    // by a nearly zero sinc at the axial defocus minima.
    let normalization = if j0.abs() >= j1.abs() {
        j0 / values[0]
    } else {
        j1 / values[1]
    };
    std::array::from_fn(|n| (2 * n + 1) as f64 * values[n] * normalization)
}

fn odd_bessel_over_radius(x: f64) -> [f64; TERMS] {
    if x < 0.5 {
        // Evaluate J_(2n+1)(x)/x directly, including arbitrarily small x.
        // Dividing a high-order underflowed J by x would lose useful digits.
        let mut leading = 0.5;
        return std::array::from_fn(|n| {
            let order = 2 * n + 1;
            if n != 0 {
                leading *= (0.5 * x) / (order - 1) as f64;
                leading *= (0.5 * x) / order as f64;
            }
            let (mut term, mut sum) = (1.0, 1.0);
            for k in 1..=12 {
                term *= -0.25 * x * x / (k * (order + k)) as f64;
                sum += term;
            }
            leading * sum
        });
    }
    if x <= HANKEL_SWITCH {
        let mut values = [0.0; MILLER_START + 2];
        values[MILLER_START] = 1.0;
        for n in (1..=MILLER_START).rev() {
            values[n - 1] = (2 * n) as f64 / x * values[n] - values[n + 1];
            if values[n - 1].abs() > 1.0e200 {
                for value in &mut values {
                    *value *= 1.0e-200;
                }
            }
        }
        // Neumann normalization J0 + 2 J2 + 2 J4 + ... = 1 avoids
        // normalizing at a zero of any one cylindrical Bessel function.
        let (mut normalization, mut correction) = (values[0], 0.0);
        for n in (2..=MILLER_START).step_by(2) {
            compensated(&mut normalization, &mut correction, 2.0 * values[n]);
        }
        return std::array::from_fn(|n| values[2 * n + 1] / normalization / x);
    }
    // All requested orders are <=79<x. Forward recurrence is stable here,
    // with the existing shared J1 owner and an order-zero Hankel start.
    let mut values = [0.0; 2 * TERMS];
    values[0] = j0_hankel(x);
    values[1] = bessel_j1(x);
    for n in 1..2 * TERMS - 1 {
        values[n + 1] = (2 * n) as f64 / x * values[n] - values[n - 1];
    }
    std::array::from_fn(|n| values[2 * n + 1] / x)
}

// Only evaluated for x>128, where the generated order-zero Hankel series
// reaches far below f64 rounding before its terms grow. Angle addition keeps
// huge arguments and neighborhoods of zeros from losing phase in x-pi/4.
fn j0_hankel(x: f64) -> f64 {
    let (mut p, mut q, mut term) = (0.0, 0.0, 1.0_f64);
    for m in 0..32 {
        match m % 4 {
            0 => p += term,
            1 => q += term,
            2 => p -= term,
            _ => q -= term,
        }
        if term.abs() < 1.0e-19 {
            break;
        }
        let m = f64::from(m);
        term *= -(2.0 * m + 1.0) * (2.0 * m + 1.0) / (8.0 * (m + 1.0)) / x;
    }
    let (sin_x, cos_x) = det::sin_cos(x);
    let amplitude = det::sqrt(core::f64::consts::FRAC_2_PI / x) * core::f64::consts::FRAC_1_SQRT_2;
    amplitude * ((p + q) * cos_x + (p - q) * sin_x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn axial_defocus_and_focused_limits() {
        for phase in [
            0.0,
            1.0e-12,
            0.01,
            0.999999,
            1.0,
            2.0,
            2.0 * PI,
            4.0 * PI,
            16.0,
        ] {
            let pupil = CircularDefocus::new(phase);
            let expected = if phase == 0.0 {
                1.0
            } else {
                (0.5 * phase).sin() / (0.5 * phase)
            };
            assert!((4.0 * PI * pupil.density(0.0) - expected * expected).abs() < 3.0e-15);
            for radius in [
                1.0e-150, 0.01, 0.5, 8.0, 79.0, 128.0, 128.001, 1.0e4, 1.0e17,
            ] {
                let value = pupil.density(radius);
                assert!(value.is_finite() && value >= 0.0);
                assert_eq!(
                    value.to_bits(),
                    CircularDefocus::new(-phase).density(radius).to_bits()
                );
            }
        }
        for radius in [0.0, 1.0e-5, 0.5, 3.8317059702075125, 20.0, 128.0, 1.0e5] {
            let expected = jinc(radius).powi(2) / (4.0 * PI);
            assert_eq!(
                CircularDefocus::new(0.0).density(radius).to_bits(),
                expected.to_bits()
            );
        }
    }

    #[test]
    fn defocus_matches_independent_radial_integrals() {
        // Independent adaptive quadrature of the defining complex integral,
        // using SciPy's J0 (QUADPACK epsabs=epsrel=1e-13), not the Zernike
        // expansion or either production recurrence. The anchors straddle
        // every branch and include the geometric blur rim and shadow tail.
        let landmarks = [
            (0.499999, 0.999999, 0.0687094518617747),
            (0.500001, 1.000001, 0.0687093938820057),
            (0.5, 16.0, 0.001143386529250452),
            (3.8317059702075125, 2.0 * PI, 0.003901632748242573),
            (8.0, 4.0 * PI, 0.0003889829028078573),
            (16.0, 16.0, 0.00036218322872995015),
            (31.0, 16.0, 0.00008660626892535932),
            (79.0, 16.0, 1.429443349864231e-7),
            (127.999, 16.0, 1.0991019952003596e-7),
            (128.0, 1.0, 9.664898836290083e-8),
            (128.001, 16.0, 1.0991101981043534e-7),
            (200.0, 2.0 * PI, 2.366170347340552e-8),
            (1000.0, 16.0, 7.331850445888677e-12),
        ];
        for (radius, phase, expected) in landmarks {
            let actual = CircularDefocus::new(phase).density(radius);
            assert!(
                (actual - expected).abs() < 2.0e-11 * expected,
                "q={radius}, phase={phase}: {actual} != {expected}"
            );
        }
    }

    #[test]
    fn defocus_packet_conserves_power_with_complete_tail_support() {
        const SAMPLES: u32 = 65_536;
        for phase in [0.0_f64, 2.0, 2.0 * PI, 10.0] {
            let ratios = [550.0 / 400.0, 1.0, 550.0 / 700.0];
            let pupils = ratios.map(|ratio| CircularDefocus::new(phase * ratio));
            let proposal_radius = 2.0_f64.max(2.0 * phase);
            let mut power = [0.0; 3];
            for sample in 0..SAMPLES {
                let u = (f64::from(sample) + 0.5) / f64::from(SAMPLES);
                let remainder = 1.0 - u;
                let radius = proposal_radius * (u * (2.0 - u)).sqrt() / remainder;
                let proposal = remainder.powi(3) / (2.0 * PI * proposal_radius.powi(2));
                for lane in 0..3 {
                    power[lane] += ratios[lane].powi(2)
                        * pupils[lane].density(ratios[lane] * radius)
                        / proposal
                        / f64::from(SAMPLES);
                }
            }
            for value in power {
                assert!(
                    (value - 1.0).abs() < 0.004,
                    "phase={phase}: power={power:?}"
                );
            }
            // No radial cutoff is allowed, including the last representable
            // open-unit draw. A quadrature-raster alias would grow this ratio.
            for u in [0.0, 0.5, 0.999, 1.0 - 0.5 * f64::EPSILON] {
                let remainder = 1.0 - u;
                let radius = proposal_radius * (u * (2.0 - u)).sqrt() / remainder;
                let proposal = remainder.powi(3) / (2.0 * PI * proposal_radius.powi(2));
                for lane in 0..3 {
                    let weight = ratios[lane].powi(2) * pupils[lane].density(ratios[lane] * radius)
                        / proposal;
                    assert!(weight.is_finite() && (0.0..64.0).contains(&weight));
                }
            }
        }
    }
}
