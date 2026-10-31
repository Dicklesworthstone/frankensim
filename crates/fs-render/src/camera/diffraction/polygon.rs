//! Focused diffraction of a uniformly illuminated regular polygon pupil.
//!
//! Coordinates are normalized by the physical circumradius. For unit-pupil
//! area A and Fourier transform F(q), Parseval gives the unit-integral density
//! `|F(q)|² / ((2*pi)²*A)`. The divergence theorem evaluates F as one analytic
//! straight-edge integral per blade; no pupil raster or truncated FFT is used.
//! See Gallatin, *J. Math. Phys.* 53, 013509 (2012), equations (5) and (21):
//! <https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=909668>.
//!
//! An edge-aligned product-Cauchy mixture covers the narrow inverse-square
//! blade-normal ridges and the generic inverse-fourth-power intensity tails.
//! An isotropic inverse-cubic proposal alone poorly resolves those ridges.

use fs_math::det;

#[path = "polygon/defocus.rs"]
mod defocus;

const PI: f64 = core::f64::consts::PI;
const PROPOSAL_SCALE: f64 = 2.0;

/// Borrow the already validated physical polygon; construction and evaluation
/// allocate nothing. All loops are bounded by the aperture's 64-blade limit.
pub(super) struct PolygonPupil<'a> {
    radius_m: f64,
    vertices: &'a [[f64; 2]],
    area: f64,
    second_moment: f64,
    frames: usize,
}

impl<'a> PolygonPupil<'a> {
    pub(super) fn new(radius_m: f64, vertices: &'a [[f64; 2]]) -> Self {
        let count = vertices.len();
        let mut twice_area = 0.0;
        let mut correction = 0.0;
        for index in 0..count {
            let a = vertices[index].map(|v| v / radius_m);
            let b = vertices[(index + 1) % count].map(|v| v / radius_m);
            add_compensated(&mut twice_area, &mut correction, cross(a, b));
        }
        Self {
            radius_m,
            vertices,
            area: 0.5 * twice_area,
            // Each coordinate has this second moment under a uniform regular
            // unit-circumradius polygon. It supplies the removable DC limit.
            second_moment: (2.0 + det::cos(2.0 * PI / count as f64)) / 12.0,
            // Product Cauchy is invariant under quarter turns and sign flips.
            // Avoid evaluating geometrically identical mixture components.
            frames: if count % 4 == 0 {
                count / 4
            } else if count % 2 == 0 {
                count / 2
            } else {
                count
            },
        }
    }

    /// Density with respect to dqx dqy, where q=2*pi*a*theta/lambda.
    pub(super) fn density(&self, q: [f64; 2]) -> f64 {
        let radius_squared = dot(q, q);
        let normalization = 4.0 * PI * PI * self.area;
        if radius_squared < 1.0e-8 {
            // Odd moments only affect the imaginary amplitude; their squared
            // contribution starts at sixth order. The omitted fourth-order
            // intensity term here is below f64 rounding at this threshold.
            return self.area * self.area * (1.0 - self.second_moment * radius_squared)
                / normalization;
        }
        let mut real = 0.0;
        let mut imaginary = 0.0;
        let mut real_correction = 0.0;
        let mut imaginary_correction = 0.0;
        for index in 0..self.vertices.len() {
            let (a, b) = self.edge(index);
            let edge = [b[0] - a[0], b[1] - a[1]];
            let midpoint = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
            let half_edge_phase = 0.5 * dot(q, edge);
            let middle_phase = dot(q, midpoint);
            let sinc = sinc(half_edge_phase);
            let (sin, cos) = det::sin_cos(middle_phase);
            let coefficient = cross(q, edge) / radius_squared;
            add_compensated(&mut real, &mut real_correction, coefficient * sinc * sin);
            let imaginary_factor = if radius_squared < 1.0 {
                // Sum of the edge normals is zero. Subtract its exactly zero
                // integral analytically, before the ill-conditioned summation.
                sinc_minus_one(half_edge_phase) * cos - 2.0 * det::sin(0.5 * middle_phase).powi(2)
            } else {
                // In the far tail, avoid subtracting O(1/r) constants from
                // edge amplitudes that may themselves be only O(1/r²).
                sinc * cos
            };
            add_compensated(
                &mut imaginary,
                &mut imaginary_correction,
                coefficient * imaginary_factor,
            );
        }
        (real * real + imaginary * imaginary) / normalization
    }

    /// One exact mixture selection and two conditional Cauchy coordinates from
    /// two independent unit variates; no rejection loop or authored tail cutoff.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    pub(super) fn sample(&self, sample: [f64; 2]) -> ([f64; 2], f64) {
        let component = sample[0] * self.frames as f64;
        let index = component as usize;
        let x = cauchy(component - index as f64);
        let y = cauchy(sample[1]);
        let normal = self.normal(index);
        let q = [normal[0] * x - normal[1] * y, normal[1] * x + normal[0] * y];
        (q, self.proposal_density(q))
    }

    /// Broaden the blade-aligned mixture with the propagated geometrical
    /// pupil footprint, while preserving its full directional tail support.
    pub(super) fn sample_defocused(&self, sample: [f64; 2], phase: f64) -> ([f64; 2], f64) {
        let (q, density) = self.sample(sample);
        let widening = phase.abs().max(1.0);
        (
            q.map(|value| value * widening),
            density / widening / widening,
        )
    }

    fn proposal_density(&self, q: [f64; 2]) -> f64 {
        let mut density = 0.0;
        for index in 0..self.frames {
            let normal = self.normal(index);
            let x = dot(q, normal) / PROPOSAL_SCALE;
            let y = cross(normal, q) / PROPOSAL_SCALE;
            density += (1.0 / (1.0 + x * x)) * (1.0 / (1.0 + y * y));
        }
        density / (self.frames as f64 * PI * PI * PROPOSAL_SCALE * PROPOSAL_SCALE)
    }

    fn edge(&self, index: usize) -> ([f64; 2], [f64; 2]) {
        (
            self.vertices[index].map(|v| v / self.radius_m),
            self.vertices[(index + 1) % self.vertices.len()].map(|v| v / self.radius_m),
        )
    }

    fn normal(&self, index: usize) -> [f64; 2] {
        let (a, b) = self.edge(index);
        let normal = [a[0] + b[0], a[1] + b[1]];
        let norm = det::sqrt(dot(normal, normal));
        normal.map(|v| v / norm)
    }
}

fn cauchy(u: f64) -> f64 {
    // The input RNG has finite precision. Give exact CDF endpoints their
    // nearest open-interval representatives rather than manufacturing infinity.
    let open = u.clamp(0.5 * f64::EPSILON, 1.0 - 0.5 * f64::EPSILON);
    PROPOSAL_SCALE * det::tan(PI * (open - 0.5))
}

fn sinc(x: f64) -> f64 {
    if x.abs() < 0.125 {
        1.0 + sinc_minus_one(x)
    } else {
        det::sin(x) / x
    }
}

fn sinc_minus_one(x: f64) -> f64 {
    if x.abs() < 0.125 {
        let square = x * x;
        square
            * (-1.0 / 6.0
                + square
                    * (1.0 / 120.0
                        + square
                            * (-1.0 / 5040.0 + square * (1.0 / 362_880.0 - square / 39_916_800.0))))
    } else {
        det::sin(x) / x - 1.0
    }
}

fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn add_compensated(sum: &mut f64, correction: &mut f64, value: f64) {
    let corrected = value - *correction;
    let updated = *sum + corrected;
    *correction = (updated - *sum) - corrected;
    *sum = updated;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vertices(count: usize, rotation: f64, radius: f64) -> Vec<[f64; 2]> {
        (0..count)
            .map(|i| {
                let (sin, cos) = det::sin_cos(rotation + 2.0 * PI * i as f64 / count as f64);
                [radius * cos, radius * sin]
            })
            .collect()
    }

    fn rotate(q: [f64; 2], angle: f64) -> [f64; 2] {
        let (sin, cos) = det::sin_cos(angle);
        [q[0] * cos - q[1] * sin, q[0] * sin + q[1] * cos]
    }

    /// Independent rectangular-aperture integral, including DC, edge-parallel
    /// directions and both sides of the small-phase branch.
    #[test]
    fn square_matches_separable_sinc() {
        let points = vertices(4, PI / 4.0, 1.0);
        let pupil = PolygonPupil::new(1.0, &points);
        for x in [
            0.0, 1.0e-10, 0.000_099, 0.000_101, 0.1, 1.0, 3.0, 10.0, 1000.0,
        ] {
            for y in [0.0, 1.0e-5, 0.3, 2.0, 11.0] {
                let oracle_sinc = |t: f64| {
                    if t == 0.0 { 1.0 } else { t.sin() / t }
                };
                let expected = 2.0 / (4.0 * PI * PI)
                    * (oracle_sinc(x / core::f64::consts::SQRT_2)
                        * oracle_sinc(y / core::f64::consts::SQRT_2))
                    .powi(2);
                let actual = pupil.density([x, y]);
                assert!(
                    (actual - expected).abs() < 1.0e-12 * expected.abs().max(1.0e-8),
                    "q=({x}, {y}): {actual} vs {expected}"
                );
            }
        }
        assert!(pupil.density([PI * core::f64::consts::SQRT_2, 0.0]) < 1.0e-28);
    }

    #[test]
    fn polygon_rotation_scale_and_complete_sampling_support() {
        for count in [3, 4, 6, 9, 16, 64] {
            let original = vertices(count, 0.0, 0.007);
            let rotated = vertices(count, 0.37, 0.021);
            let left = PolygonPupil::new(0.007, &original);
            let right = PolygonPupil::new(0.021, &rotated);
            for q in [[0.0; 2], [1.0e-5, 2.0e-5], [1.0, -2.0], [17.0, 4.0]] {
                let a = left.density(q);
                let b = right.density(rotate(q, 0.37));
                assert!((a - b).abs() < 2.0e-13 * a.max(1.0e-8));
            }
            for sample in [[0.0; 2], [0.5, 0.0], [0.0, 0.5], [1.0_f64.next_down(); 2]] {
                let (q, proposal) = left.sample(sample);
                assert!(q.into_iter().all(f64::is_finite));
                assert!(proposal.is_finite() && proposal > 0.0);
                let weight = left.density(q) / proposal;
                assert!(weight.is_finite() && weight >= 0.0);
            }
        }
    }

    /// Fixed Hammersley cubature of the actual sampler and mixture PDF.
    /// Complete tails and all packet wavelengths retain unit integrated power.
    #[test]
    fn polygon_packet_conserves_radiance() {
        const SAMPLES: u32 = 65_536;
        for count in [3, 4, 6, 9, 64] {
            let points = vertices(count, 0.234, 1.0);
            let pupil = PolygonPupil::new(1.0, &points);
            let mut totals = [0.0; 3];
            for i in 0..SAMPLES {
                let u = (f64::from(i) + 0.5) / f64::from(SAMPLES);
                let v = f64::from(i.reverse_bits()) / 4_294_967_296.0 + 0.5 / f64::from(SAMPLES);
                let (q, proposal) = pupil.sample([u, v]);
                for (sum, ratio) in totals.iter_mut().zip([550.0 / 400.0, 1.0, 550.0 / 700.0]) {
                    *sum += ratio * ratio * pupil.density(q.map(|v| ratio * v)) / proposal;
                }
            }
            for total in totals {
                let power = total / f64::from(SAMPLES);
                assert!(
                    (power - 1.0).abs() < 0.004,
                    "{count}-blade integrated power {power}"
                );
            }
        }
    }
}
