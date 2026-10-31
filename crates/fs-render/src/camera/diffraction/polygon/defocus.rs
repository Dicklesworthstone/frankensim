//! Fresnel propagation through the actual polygon, including all image tails.
//!
//! F(q,phi)=integral_P exp(i*phi*|r|²-i*q.r) dA. A regularized divergence
//! formula handles the stationary point; an exterior edge formula handles
//! the far field without enormous completed-square phases. Only bounded
//! tangential frequencies use Gauss quadrature, with an explicit derivative
//! remainder. Unbounded frequencies use endpoint integration by parts with
//! a separately bounded remainder. No sampled-pupil alias or cutoff occurs.
//!
//! Boundary identities: Barnett, arXiv:2010.05978, Appendix A, Eqs.19/22;
//! Huang, Christian & McDonald, JOSA A23,2768 (2006), Eqs.14–16.
//! <https://arxiv.org/pdf/2010.05978>
//! <https://www.salfordphysics.com/gsmcdonald/32_JOSAA_v23_p2768_2006.pdf>

use super::{PI, PolygonPupil, cross, dot};
use fs_math::{c64::C64, det};

const TOLERANCE: f64 = 2.0e-13;
const MAX_SEGMENTS: usize = 64;
const ENDPOINT_TERMS: usize = 32;
const ENDPOINT_FREQUENCY: f64 = 96.0;

impl PolygonPupil<'_> {
    /// Unit-integral density per dqx dqy. None refuses unresolved numerical
    /// work instead of substituting a focused or geometrical blur.
    pub(crate) fn defocused_density(&self, q: [f64; 2], phase: f64) -> Option<f64> {
        if !phase.is_finite() || phase.abs() > 16.0 || q.iter().any(|v| !v.is_finite()) {
            return None;
        }
        if phase == 0.0 {
            let value = self.density(q);
            return value.is_finite().then_some(value);
        }
        let radius = det::hypot(q[0], q[1]);
        if !radius.is_finite() {
            return None;
        }
        let normalization = 4.0 * PI * PI * self.area;
        if radius < 1.0e-5 && phase.abs() < 1.0e-5 {
            // For independent uniform pupil points, intensity/area² is the
            // expectation of cos(S(r)-S(r')). The fourth-order remainder is
            // <=[2*(|q|+|phi|)]^4/24 <1.07e-19 in this branch.
            let c = det::cos(2.0 * PI / self.vertices.len() as f64);
            let mean_r2 = (2.0 + c) / 6.0;
            let variance_r2 = (7.0 + 6.0 * c + 2.0 * c * c) / 45.0 - mean_r2 * mean_r2;
            return Some(
                self.area
                    * self.area
                    * (1.0 - self.second_moment * radius * radius - variance_r2 * phase * phase)
                    / normalization,
            );
        }
        let mut sum = ComplexSum::default();
        if radius <= 4.0 * phase.abs() {
            // c=q/(2phi) has |c|<=2. With exprel(z)=(exp(z)-1)/z,
            // F=.5 exp(-i phi|c|²) integral_boundary (r-c).n exprel(i phi|r-c|²) ds.
            // The common unit phase is immaterial to intensity. The removable
            // limit is valid even when c lies exactly on an edge or vertex.
            let center = q.map(|value| 0.5 * value / phase);
            for index in 0..self.vertices.len() {
                let (a, b) = self.edge(index);
                let midpoint = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
                let half = [0.5 * (b[0] - a[0]), 0.5 * (b[1] - a[1])];
                let shift = [midpoint[0] - center[0], midpoint[1] - center[1]];
                let prefactor = 0.5 * cross(shift, half);
                if prefactor == 0.0 {
                    continue;
                }
                let curvature = phase * dot(half, half);
                let linear = 2.0 * phase * dot(shift, half);
                let bound =
                    derivative_bounds(curvature.abs(), linear.abs() + 2.0 * curvature.abs(), 0.0);
                let rule = Quadrature::admit(bound, TOLERANCE)?;
                let integral = rule.integrate(|t| {
                    let r = [shift[0] + t * half[0], shift[1] + t * half[1]];
                    exprel_i(phase * dot(r, r))
                });
                sum.add(integral.scale(prefactor));
            }
        } else {
            // Outside the pupil, div[(q-2phi*r)/|q-2phi*r|²]=0 and winding
            // contributes nothing. F=i integral_boundary w.n exp(iS)/|w|² ds.
            // Normalize w by |q| so huge image samples do not square |q|.
            let unit_q = q.map(|value| value / radius);
            let phase_scale = 2.0 * phase / radius;
            let vertex_radius = self
                .vertices
                .iter()
                .map(|v| det::hypot(v[0] / self.radius_m, v[1] / self.radius_m))
                .fold(1.0_f64, f64::max)
                * (1.0 + 16.0 * f64::EPSILON);
            let separation = radius - 2.0 * phase.abs() * vertex_radius;
            if separation <= 0.0 || !separation.is_finite() {
                return None;
            }
            let far = radius >= 8.0 * phase.abs() * vertex_radius + 32.0;
            for index in 0..self.vertices.len() {
                let (a, b) = self.edge(index);
                let midpoint = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
                let half = [0.5 * (b[0] - a[0]), 0.5 * (b[1] - a[1])];
                let length = det::hypot(half[0], half[1]);
                let w = [
                    unit_q[0] - phase_scale * midpoint[0],
                    unit_q[1] - phase_scale * midpoint[1],
                ];
                let delta = half.map(|v| -phase_scale * v);
                let coefficient = cross(w, half) / radius;
                if coefficient == 0.0 {
                    continue;
                }
                let denominator = [dot(w, w), 2.0 * dot(w, delta), dot(delta, delta)];
                let curvature = phase * dot(half, half);
                // Keep the tiny actual midpoint/edge non-orthogonality. It
                // need not be rounded to zero to use the analytic formulas.
                let frequency = dot(q, half) - 2.0 * phase * dot(midpoint, half);
                let eta = 2.0 * phase.abs() * length / separation;
                let integral = if far && frequency.abs() >= ENDPOINT_FREQUENCY {
                    let bound = derivative_bounds(curvature.abs(), 2.0 * curvature.abs(), eta);
                    let mut remainder = 2.0 * bound[ENDPOINT_TERMS];
                    for n in 1..=ENDPOINT_TERMS {
                        remainder *= n as f64 / frequency.abs();
                    }
                    // The omitted edge amplitude is <=length/separation times
                    // this bound. The threshold makes it smaller than 1e-13.
                    if !remainder.is_finite() || remainder > TOLERANCE {
                        return None;
                    }
                    endpoint_integral(
                        denominator,
                        coefficient,
                        curvature,
                        frequency,
                        a,
                        b,
                        q,
                        phase,
                    )
                } else {
                    let mut bound = derivative_bounds(
                        curvature.abs(),
                        frequency.abs() + 2.0 * curvature.abs(),
                        eta,
                    );
                    let subtract_constant = radius < 1.0;
                    if subtract_constant {
                        let mut power = 1.0;
                        for value in &mut bound {
                            *value += power;
                            power *= eta;
                        }
                    }
                    // The low-q subtraction has O(1) total amplitude, despite
                    // the un-subtracted boundary envelope growing like 1/q.
                    let rule = Quadrature::admit(bound, TOLERANCE * separation.min(1.0))?;
                    let common = exp_i(phase * dot(midpoint, midpoint)) * exp_i(-dot(q, midpoint));
                    rule.integrate(|t| {
                        let den = denominator[0] + t * (denominator[1] + t * denominator[2]);
                        let amplitude = coefficient / den;
                        let wave = if subtract_constant {
                            let r = [midpoint[0] + t * half[0], midpoint[1] + t * half[1]];
                            expm1_i(phase * dot(r, r) - dot(q, r))
                        } else {
                            common * exp_i(curvature * t * t) * exp_i(-frequency * t)
                        };
                        wave.scale(amplitude)
                    })
                };
                // Multiplication by the common i does not affect intensity.
                sum.add(integral);
            }
        }
        let amplitude = sum.value();
        let density = (amplitude.re * amplitude.re + amplitude.im * amplitude.im) / normalization;
        (density.is_finite() && density >= 0.0).then_some(density)
    }
}

/// a_n bounds |exp(i quadratic)^(n)|/n!, then convolution with eta^n
/// bounds the rational edge amplitude as well. All arithmetic is positive.
fn derivative_bounds(curvature: f64, slope: f64, eta: f64) -> [f64; 65] {
    let mut values = [0.0; 65];
    let (mut previous, mut older) = (1.0, 0.0);
    values[0] = 1.0;
    for n in 1..values.len() {
        let coefficient = (slope * previous + 2.0 * curvature * older) / n as f64;
        values[n] = coefficient + eta * values[n - 1];
        older = previous;
        previous = coefficient;
    }
    values
}

#[derive(Clone, Copy)]
struct Quadrature {
    segments: usize,
    order: usize,
}

impl Quadrature {
    fn admit(bounds: [f64; 65], tolerance: f64) -> Option<Self> {
        let mut selected: Option<Self> = None;
        // Gauss N remainder: C_N sup|f^(2N)|/(2N)!,
        // C_N=2^(2N+1)(N!)^4/((2N+1)((2N)!)²). A factor four guards
        // rounding of the positive bound; it does not enclose kernel roundoff.
        for (order, constant) in [(16, 7.204615332388622e-10), (32, 1.690011092571857e-19)] {
            let mut error = 4.0 * constant * bounds[2 * order];
            let mut segments = 1;
            while error > tolerance && segments < MAX_SEGMENTS {
                segments *= 2;
                error *= if order == 16 {
                    1.0 / 4_294_967_296.0
                } else {
                    1.0 / 18_446_744_073_709_551_616.0
                };
            }
            if error.is_finite()
                && error <= tolerance
                && selected.is_none_or(|best| segments * order < best.segments * best.order)
            {
                selected = Some(Self { segments, order });
            }
        }
        selected
    }

    fn integrate(self, f: impl Fn(f64) -> C64) -> C64 {
        let nodes: &[(f64, f64)] = if self.order == 16 { &GL16 } else { &GL32 };
        let half_width = 1.0 / self.segments as f64;
        let mut sum = ComplexSum::default();
        for segment in 0..self.segments {
            let midpoint = -1.0 + (2 * segment + 1) as f64 * half_width;
            for &(node, weight) in nodes {
                let offset = node * half_width;
                sum.add((f(midpoint - offset) + f(midpoint + offset)).scale(weight * half_width));
            }
        }
        sum.value()
    }
}

#[allow(clippy::too_many_arguments)]
fn endpoint_integral(
    d: [f64; 3],
    coefficient: f64,
    curvature: f64,
    frequency: f64,
    a: [f64; 2],
    b: [f64; 2],
    q: [f64; 2],
    phase: f64,
) -> C64 {
    let mut result = ComplexSum::default();
    for (endpoint, vertex, sign) in [(-1.0, a, 1.0), (1.0, b, -1.0)] {
        let d0 = d[0] + endpoint * d[1] + d[2];
        let d1 = d[1] + 2.0 * endpoint * d[2];
        let d2 = d[2];
        let factors = [
            C64::new(-d1, 2.0 * curvature * endpoint * d0),
            C64::new(-2.0 * d2, 2.0 * curvature * (d0 + endpoint * d1)),
            C64::new(0.0, 2.0 * curvature * (d1 + endpoint * d2)),
            C64::new(0.0, 2.0 * curvature * d2),
        ];
        // Taylor coefficients of R(t)*exp(i curvature*(t²-endpoint²))
        // from its polynomial first-order ODE; fixed-size, linear-time work.
        let mut coefficients = [C64::ZERO; ENDPOINT_TERMS];
        coefficients[0] = C64::from_re(coefficient / d0);
        for n in 0..ENDPOINT_TERMS - 1 {
            let mut value = (factors[0] - C64::from_re(d1 * n as f64)) * coefficients[n];
            if n >= 1 {
                value =
                    value + (factors[1] - C64::from_re(d2 * (n - 1) as f64)) * coefficients[n - 1];
            }
            if n >= 2 {
                value = value + factors[2] * coefficients[n - 2];
            }
            if n >= 3 {
                value = value + factors[3] * coefficients[n - 3];
            }
            coefficients[n + 1] = value.scale(1.0 / (d0 * (n + 1) as f64));
        }
        let inverse = C64::new(0.0, -1.0 / frequency);
        let mut factor = inverse;
        let mut boundary = ComplexSum::default();
        for (n, coefficient) in coefficients.into_iter().enumerate() {
            boundary.add(coefficient * factor);
            factor = (factor * inverse).scale((n + 1) as f64);
        }
        // Keep the bounded defocus phase separate from the huge linear phase.
        let wave = exp_i(phase * dot(vertex, vertex)) * exp_i(-dot(q, vertex));
        result.add((boundary.value() * wave).scale(sign));
    }
    result.value()
}

fn exp_i(value: f64) -> C64 {
    let (sin, cos) = det::sin_cos(value);
    C64::new(cos, sin)
}

fn expm1_i(value: f64) -> C64 {
    C64::new(-2.0 * det::sin(0.5 * value).powi(2), det::sin(value))
}

fn exprel_i(value: f64) -> C64 {
    let half = 0.5 * value;
    let sinc = if half.abs() < 1.0e-4 {
        1.0 + half * half * (-1.0 / 6.0 + half * half / 120.0)
    } else {
        det::sin(half) / half
    };
    exp_i(half).scale(sinc)
}

#[derive(Default)]
struct ComplexSum {
    re: f64,
    im: f64,
    cr: f64,
    ci: f64,
}
impl ComplexSum {
    fn add(&mut self, value: C64) {
        super::add_compensated(&mut self.re, &mut self.cr, value.re);
        super::add_compensated(&mut self.im, &mut self.ci, value.im);
    }
    fn value(&self) -> C64 {
        C64::new(self.re, self.im)
    }
}

const GL16: [(f64, f64); 8] = [
    (0.09501250983763744, 0.18945061045506864),
    (0.2816035507792589, 0.18260341504492364),
    (0.45801677765722737, 0.16915651939500265),
    (0.6178762444026438, 0.1495959888165767),
    (0.755404408355003, 0.12462897125553407),
    (0.8656312023878318, 0.0951585116824926),
    (0.9445750230732326, 0.06225352393864746),
    (0.9894009349916499, 0.027152459411754176),
];
const GL32: [(f64, f64); 16] = [
    (0.048307665687738324, 0.09654008851472766),
    (0.1444719615827965, 0.09563872007927471),
    (0.23928736225213706, 0.09384439908080451),
    (0.33186860228212767, 0.09117387869576378),
    (0.42135127613063533, 0.08765209300440378),
    (0.5068999089322294, 0.0833119242269467),
    (0.5877157572407623, 0.07819389578707023),
    (0.6630442669302152, 0.07234579410884834),
    (0.7321821187402897, 0.06582222277636168),
    (0.7944837959679424, 0.058684093478535565),
    (0.84936761373257, 0.05099805926237609),
    (0.8963211557660521, 0.042835898022226836),
    (0.9349060759377397, 0.034273862913021765),
    (0.9647622555875064, 0.025392065309262024),
    (0.9856115115452684, 0.016274394730905743),
    (0.9972638618494816, 0.007018610009470506),
];

#[cfg(test)]
mod tests {
    use super::*;

    fn vertices(count: usize, angle: f64) -> Vec<[f64; 2]> {
        (0..count)
            .map(|n| {
                let (sin, cos) = det::sin_cos(angle + 2.0 * PI * n as f64 / count as f64);
                [cos, sin]
            })
            .collect()
    }

    #[test]
    fn propagated_polygon_matches_independent_square_and_triangle_integrals() {
        let h = core::f64::consts::FRAC_1_SQRT_2;
        let square = [[h, h], [-h, h], [-h, -h], [h, -h]];
        // Independent separable Fresnel C/S integrals on [-1/sqrt(2),1/sqrt(2)]²,
        // evaluated by SciPy special.fresnel; no polygon boundary formula.
        let square_landmarks = [
            ([0.0, 0.0], 8.0, 0.002345119344354164),
            ([0.01, 0.03], 16.0, 0.0008789445415953272),
            ([3.0, 2.0], 6.0, 0.0035869843911316883),
            ([31.999, 0.4], 8.0, 0.00001225799331028901),
            ([32.001, 0.4], 8.0, 0.000012288620969230895),
            ([135.76436056646091, 1.7], 8.0, 1.1357884239691948e-6),
            ([135.76464340917337, 1.7], 8.0, 1.135700817608941e-6),
            ([200.0, 131.0], 16.0, 4.899028581433692e-12),
            ([1000.0, 3.0], 16.0, 4.894724287012833e-10),
            ([1000.0, 777.0], 8.0, 2.48537421358121e-15),
        ];
        let triangle = [
            [1.0, 0.0],
            [-0.5, 3.0_f64.sqrt() / 2.0],
            [-0.5, -3.0_f64.sqrt() / 2.0],
        ];
        // Direct two-dimensional complex pupil-area integration, checked with
        // independent 128x128 and 256x256 Duffy/Gauss rules. Includes stationary
        // points on an actual vertex/edge, and both signs of triangular defocus.
        let triangle_landmarks = [
            ([0.0, 0.0], 8.0, 0.006488674096812713),
            ([3.0, 2.0], 6.0, 0.0053226444319934745),
            ([3.0, 2.0], -6.0, 0.004882815746736282),
            ([16.0, 0.0], 8.0, 0.00003750725080215952),
            ([4.0, 6.928203230275509], 8.0, 0.0007158863380438737),
            ([31.999, 0.1], 8.0, 5.635083747632688e-6),
            ([32.001, 0.1], 8.0, 5.6311044215431065e-6),
            ([96.0, 53.0], 16.0, 1.1269553437948862e-8),
            ([200.0, 113.0], 16.0, 3.2295951320355726e-11),
        ];
        for (points, landmarks) in [
            (&square[..], &square_landmarks[..]),
            (&triangle[..], &triangle_landmarks[..]),
        ] {
            let pupil = PolygonPupil::new(1.0, points);
            for &(q, phase, expected) in landmarks {
                let actual = pupil.defocused_density(q, phase).unwrap();
                let tolerance = if q[0] > 150.0 { 2.0e-8 } else { 2.0e-10 };
                assert!(
                    (actual - expected).abs() < tolerance * expected,
                    "{} blades, q={q:?}, phase={phase}: {actual} != {expected}",
                    points.len()
                );
            }
        }
    }

    #[test]
    fn propagated_polygon_keeps_signed_phase_rotation_and_full_sampling_support() {
        for count in [3, 4, 7, 16, 64] {
            let points = vertices(count, 0.0);
            let rotated = vertices(count, 0.37);
            let pupil = PolygonPupil::new(1.0, &points);
            let other = PolygonPupil::new(1.0, &rotated);
            for phase in [1.0e-8, 0.001, 6.0, 16.0] {
                for q in [[0.0; 2], [1.0e-5, -2.0e-5], [3.0, 2.0], [79.0, -12.0]] {
                    let value = pupil.defocused_density(q, phase).unwrap();
                    let reflected = pupil.defocused_density(q.map(|v| -v), -phase).unwrap();
                    assert!((value - reflected).abs() < 1.0e-12 * value.max(1.0e-8));
                    let transformed = [
                        q[0] * 0.37_f64.cos() - q[1] * 0.37_f64.sin(),
                        q[0] * 0.37_f64.sin() + q[1] * 0.37_f64.cos(),
                    ];
                    let rotated_value = other.defocused_density(transformed, phase).unwrap();
                    assert!((value - rotated_value).abs() < 1.0e-10 * value.max(1.0e-8));
                }
                for sample in [
                    [0.0; 2],
                    [0.5, 0.0],
                    [0.0, 0.5],
                    [1.0 - 0.5 * f64::EPSILON; 2],
                ] {
                    let (q, proposal) = pupil.sample_defocused(sample, phase);
                    let weight = pupil.defocused_density(q, phase).unwrap() / proposal;
                    assert!(weight.is_finite() && weight >= 0.0);
                }
            }
            for q in [[0.0; 2], [1.0e-5, 2.0e-5], [3.0, 2.0], [200.0, 131.0]] {
                assert_eq!(
                    pupil.defocused_density(q, 0.0).unwrap().to_bits(),
                    pupil.density(q).to_bits()
                );
            }
        }
        let points = vertices(3, 0.0);
        let pupil = PolygonPupil::new(1.0, &points);
        let positive = pupil.defocused_density([3.0, 2.0], 6.0).unwrap();
        let negative = pupil.defocused_density([3.0, 2.0], -6.0).unwrap();
        assert!(
            (positive - negative).abs() > 0.05 * positive,
            "odd-blade defocus must not take abs(phase)"
        );
    }

    #[test]
    fn propagated_polygon_packet_conserves_power() {
        const SAMPLES: u32 = 16_384;
        for count in [3, 7] {
            let points = vertices(count, 0.17);
            let pupil = PolygonPupil::new(1.0, &points);
            let phase = 6.0;
            let ratios = [550.0 / 400.0, 1.0, 550.0 / 700.0];
            let mut power = [0.0; 3];
            for i in 0..SAMPLES {
                let sample = [
                    (f64::from(i) + 0.5) / f64::from(SAMPLES),
                    (f64::from(i.reverse_bits()) + 0.5) / 4_294_967_296.0,
                ];
                let (q, proposal) = pupil.sample_defocused(sample, phase);
                for lane in 0..3 {
                    let ratio = ratios[lane];
                    power[lane] += ratio
                        * ratio
                        * pupil
                            .defocused_density(q.map(|v| v * ratio), phase * ratio)
                            .unwrap()
                        / proposal
                        / f64::from(SAMPLES);
                }
            }
            assert!(
                power.iter().all(|value| (value - 1.0).abs() < 0.01),
                "{count}-blade propagated packet power: {power:?}"
            );
        }
    }
}
