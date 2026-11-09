//! Coherent diffraction through a finite array of disjoint rectangular openings.
//!
//! Each component is evaluated in its own circumradius, so a long baseline
//! does not turn a well-conditioned rectangle into a vanishing polygon. For
//! component/outer radius ratio r and normalized half-pitch c, translation
//! of the original Fresnel integral gives
//! `a(q,phi) = r/sqrt(N) sum_j exp(i*phi*|c_j|²-i*q.c_j)
//! a_single(r*(q-2*phi*c_j),phi*r²)`, `c_j=(2*j-N+1)*c`.
//! The normalization uses the sum of transmitting areas; dark gaps transmit
//! nothing. At focus the array factor is `sin(N*t)/sin(t)`, `t=q.c`.
//! A division-free binary composition evaluates that finite sum in O(log N),
//! retaining its signed removable limits at every principal diffraction order.
//! Reference: TU Delft, *Scalar Diffraction Optics*, "Diffraction gratings":
//! <https://interactivetextbooks.tudelft.nl/interactive-optics/content/Chap6_Diffraction/DiffractiveOptics_2022_01Clean.html>.

use super::polygon::PolygonPupil;
use fs_math::{c64::C64, dd::Dd, det};

/// Vertices and half-axes describe one centred opening. Successive physical
/// copies are spaced by twice `half_offset_m`, centred symmetrically as a whole.
pub(super) struct RectangularArrayPupil<'a> {
    component: PolygonPupil<'a>,
    radius_ratio: f64,
    half_offset: [f64; 2],
    openings: u16,
}

impl<'a> RectangularArrayPupil<'a> {
    pub(super) fn new(
        outer_radius_m: f64,
        component_radius_m: f64,
        vertices: &'a [[f64; 2]],
        half_axes_m: [[f64; 2]; 2],
        half_offset_m: [f64; 2],
        openings: u16,
    ) -> Self {
        Self {
            component: PolygonPupil::rectangular(component_radius_m, vertices, half_axes_m),
            radius_ratio: component_radius_m / outer_radius_m,
            half_offset: half_offset_m.map(|value| value / outer_radius_m),
            openings,
        }
    }

    /// Unit-integral density per dqx dqy, including the interference fringes.
    pub(super) fn density(&self, q: [f64; 2]) -> f64 {
        self.normalized_amplitude(q).norm_sq()
    }

    /// Exact focused sinc envelope times the signed finite-array factor.
    pub(super) fn normalized_amplitude(&self, q: [f64; 2]) -> C64 {
        // Preserve the established two-opening arithmetic and its random paths.
        let factor = if self.openings == 2 {
            core::f64::consts::SQRT_2 * det::cos(dot(q, self.half_offset))
        } else {
            array_factor(q, self.half_offset, self.openings) / det::sqrt(f64::from(self.openings))
        };
        self.component
            .normalized_amplitude(q.map(|value| self.radius_ratio * value))
            .scale(self.radius_ratio * factor)
    }

    /// The same global phase admission and numerical refusal as the original
    /// finite-pupil Fresnel propagator; no geometrical-blur fallback is used.
    pub(super) fn defocused_density(&self, q: [f64; 2], phase: f64) -> Option<f64> {
        let value = self.normalized_defocused_amplitude(q, phase)?.norm_sq();
        (value.is_finite() && value >= 0.0).then_some(value)
    }

    pub(super) fn normalized_defocused_amplitude(&self, q: [f64; 2], phase: f64) -> Option<C64> {
        if !phase.is_finite() || phase.abs() > 16.0 || q.iter().any(|value| !value.is_finite()) {
            return None;
        }
        if phase == 0.0 {
            return finite(self.normalized_amplitude(q));
        }
        if !det::hypot(q[0], q[1]).is_finite() {
            return None;
        }
        if self.openings != 2 {
            return self.array_defocused_amplitude(q, phase);
        }
        let local_phase = phase * self.radius_ratio * self.radius_ratio;
        let shift = self.half_offset.map(|value| 2.0 * phase * value);
        let positive = self.component.normalized_defocused_amplitude(
            [
                self.radius_ratio * (q[0] - shift[0]),
                self.radius_ratio * (q[1] - shift[1]),
            ],
            local_phase,
        )?;
        let negative = self.component.normalized_defocused_amplitude(
            [
                self.radius_ratio * (q[0] + shift[0]),
                self.radius_ratio * (q[1] + shift[1]),
            ],
            local_phase,
        )?;
        // Keep the common quadratic phase separate: adding it to an enormous
        // q.c would lose its coherent phase before trigonometric reduction.
        let (sin, cos) = det::sin_cos(dot(q, self.half_offset));
        let wave = positive * C64::new(cos, -sin) + negative * C64::new(cos, sin);
        let (sin, cos) = det::sin_cos(phase * dot(self.half_offset, self.half_offset));
        finite(
            (wave * C64::new(cos, sin)).scale(self.radius_ratio * core::f64::consts::FRAC_1_SQRT_2),
        )
    }

    fn center(&self, index: u16) -> [f64; 2] {
        let offset = 2.0 * f64::from(index) - f64::from(self.openings - 1);
        self.half_offset.map(|v| offset * v)
    }

    fn array_defocused_amplitude(&self, q: [f64; 2], phase: f64) -> Option<C64> {
        let local_phase = phase * self.radius_ratio * self.radius_ratio;
        let (mut real, mut imaginary) = (Dd::ZERO, Dd::ZERO);
        for index in 0..self.openings {
            let center = self.center(index);
            let local = [
                self.radius_ratio * (q[0] - 2.0 * phase * center[0]),
                self.radius_ratio * (q[1] - 2.0 * phase * center[1]),
            ];
            let amplitude = self
                .component
                .normalized_defocused_amplitude(local, local_phase)?;
            // Do not add an O(1) quadratic phase to an unbounded image phase
            // in f64: a finite far-field query would lose the former entirely.
            let argument = -Dd::from_f64(q[0]) * Dd::from_f64(center[0])
                - Dd::from_f64(q[1]) * Dd::from_f64(center[1]);
            let (linear_sin, linear_cos) = det::sin_cos_dd(argument);
            let (quadratic_sin, quadratic_cos) =
                det::sin_cos_dd(Dd::from_f64(phase * dot(center, center)));
            let cos = linear_cos * quadratic_cos - linear_sin * quadratic_sin;
            let sin = linear_sin * quadratic_cos + linear_cos * quadratic_sin;
            let re = Dd::from_f64(amplitude.re);
            let im = Dd::from_f64(amplitude.im);
            real = real + re * cos - im * sin;
            imaginary = imaginary + re * sin + im * cos;
        }
        let scale = self.radius_ratio / det::sqrt(f64::from(self.openings));
        finite(C64::new(real.to_f64(), imaginary.to_f64()).scale(scale))
    }

    /// Reuse the component's complete Cauchy proposal. The focused fringe
    /// multiplies the component PSF by the Fejer kernel, adding at most a
    /// factor N to that proposal's weights. No diffraction order is cropped.
    pub(super) fn sample(&self, sample: [f64; 2]) -> ([f64; 2], f64) {
        let (local, density) = self.component.sample(sample);
        (
            local.map(|value| value / self.radius_ratio),
            density * self.radius_ratio * self.radius_ratio,
        )
    }

    /// Select one translated Fresnel component but return the complete mixture
    /// density. Both the interference region and each geometrical footprint
    /// retain positive proposal density, including their unbounded tails.
    pub(super) fn sample_defocused(&self, sample: [f64; 2], phase: f64) -> ([f64; 2], f64) {
        if phase == 0.0 {
            return self.sample(sample);
        }
        let coordinate = f64::from(self.openings) * sample[0];
        let index = coordinate.floor() as u16;
        let u = coordinate - f64::from(index);
        let center = self.center(index);
        let local_phase = phase * self.radius_ratio * self.radius_ratio;
        let (local, _) = self.component.sample_defocused([u, sample[1]], local_phase);
        let q = [
            local[0] / self.radius_ratio + 2.0 * phase * center[0],
            local[1] / self.radius_ratio + 2.0 * phase * center[1],
        ];
        (q, self.defocused_proposal_density(q, phase))
    }

    /// One source-depth-independent narrow/broad mixture. In component units,
    /// the largest translated-center excursion is 2*max_phase*r*|c|. Add its
    /// half-width to the scale-2 Cauchy envelope's Fresnel broadening; the same
    /// proposal then covers both signs and every admitted intermediate depth.
    pub(super) fn sample_phase_range(
        &self,
        sample: [f64; 2],
        maximum_phase: f64,
    ) -> ([f64; 2], f64) {
        let maximum_phase = maximum_phase.abs();
        let widening = (maximum_phase * self.radius_ratio * self.radius_ratio).max(1.0)
            + maximum_phase
                * self.radius_ratio
                * det::hypot(self.half_offset[0], self.half_offset[1])
                * f64::from(self.openings - 1);
        let (u, selected) = if sample[0] < 0.5 {
            (2.0 * sample[0], 1.0)
        } else {
            (2.0 * sample[0] - 1.0, widening)
        };
        let (local, _) = self.component.sample([u, sample[1]]);
        let q = local.map(|value| value * selected / self.radius_ratio);
        let density = 0.5
            * (self.widened_proposal_density(q, 1.0) + self.widened_proposal_density(q, widening));
        (q, density)
    }

    fn defocused_proposal_density(&self, q: [f64; 2], phase: f64) -> f64 {
        let widening = (phase.abs() * self.radius_ratio * self.radius_ratio).max(1.0);
        let mut density = 0.0;
        for index in 0..self.openings {
            let center = self.center(index);
            density += self.widened_proposal_density(
                [
                    q[0] - 2.0 * phase * center[0],
                    q[1] - 2.0 * phase * center[1],
                ],
                widening,
            );
        }
        density / f64::from(self.openings)
    }

    fn widened_proposal_density(&self, q: [f64; 2], widening: f64) -> f64 {
        let scale = self.radius_ratio / widening;
        self.component
            .proposal_density(q.map(|value| scale * value))
            * scale
            * scale
    }
}

// Compose (cos(n*t), sin(n*t), sin(n*t)/sin(t)) without dividing by sin(t).
// The addition law A_(n+m)=C_n*A_m+C_m*A_n is valid at the removable limits
// too. Double-double phases/sums retain grating minima and integer-order signs;
// the number of bounded compositions grows logarithmically with opening count.
fn array_factor(q: [f64; 2], half_pitch: [f64; 2], mut count: u16) -> f64 {
    let phase = Dd::from_f64(q[0]) * Dd::from_f64(half_pitch[0])
        + Dd::from_f64(q[1]) * Dd::from_f64(half_pitch[1]);
    let (sin, cos) = det::sin_cos_dd(phase);
    let (mut c, mut s, mut a) = (Dd::ONE, Dd::ZERO, Dd::ZERO);
    let (mut pc, mut ps, mut pa) = (cos, sin, Dd::ONE);
    while count != 0 {
        if count & 1 != 0 {
            a = c * pa + pc * a;
            (c, s) = (c * pc - s * ps, s * pc + c * ps);
        }
        count >>= 1;
        if count != 0 {
            pa = Dd::from_f64(2.0) * pc * pa;
            (pc, ps) = (pc * pc - ps * ps, Dd::from_f64(2.0) * pc * ps);
        }
    }
    a.to_f64()
}

fn finite(value: C64) -> Option<C64> {
    (value.re.is_finite() && value.im.is_finite()).then_some(value)
}

fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::{PI, TAU};

    struct Geometry {
        width: f64,
        height: f64,
        separation: f64,
        outer: f64,
        component: f64,
        axes: [[f64; 2]; 2],
        offset: [f64; 2],
        vertices: [[f64; 2]; 4],
    }

    fn rotate(q: [f64; 2], angle: f64) -> [f64; 2] {
        let (sin, cos) = angle.sin_cos();
        [cos * q[0] - sin * q[1], sin * q[0] + cos * q[1]]
    }

    fn corners(axes: [[f64; 2]; 2]) -> [[f64; 2]; 4] {
        [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]].map(|s| {
            [
                s[0] * axes[0][0] + s[1] * axes[1][0],
                s[0] * axes[0][1] + s[1] * axes[1][1],
            ]
        })
    }

    impl Geometry {
        fn new(width: f64, height: f64, separation: f64, angle: f64) -> Self {
            let axes = [
                rotate([0.5 * width, 0.0], angle),
                rotate([0.0, 0.5 * height], angle),
            ];
            Self {
                width,
                height,
                separation,
                outer: (0.5 * (width + separation)).hypot(0.5 * height),
                component: (0.5 * width).hypot(0.5 * height),
                axes,
                offset: rotate([0.5 * separation, 0.0], angle),
                vertices: corners(axes),
            }
        }

        fn pupil(&self) -> RectangularArrayPupil<'_> {
            self.array(2)
        }

        fn array(&self, openings: u16) -> RectangularArrayPupil<'_> {
            RectangularArrayPupil::new(
                self.outer,
                self.component,
                &self.vertices,
                self.axes,
                self.offset,
                openings,
            )
        }

        fn peak(&self) -> f64 {
            (2.0 * self.width * self.height).sqrt() / (TAU * self.outer)
        }
    }

    fn sinc(value: f64) -> f64 {
        if value == 0.0 {
            1.0
        } else {
            value.sin() / value
        }
    }

    fn close(actual: C64, expected: C64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual:?} vs {expected:?}"
        );
    }

    /// G1: the finite-slit envelope and baseline fringes keep their signs,
    /// physical scale, rotation and exact total-area normalization.
    #[test]
    fn focused_pair_matches_signed_sinc_cos_and_dark_fringes() {
        for aspect in [1.0, 8.0, 256.0] {
            let geometry = Geometry::new(2.0, 2.0 / aspect, 5.0, 0.37);
            let pupil = geometry.pupil();
            for q in [
                [0.0; 2],
                [2e-7, -3e-7],
                [0.4, -1.1],
                [10.0, 0.4],
                [71.0, -9.0],
            ] {
                let expected = geometry.peak()
                    * sinc(q[0] * geometry.width / (2.0 * geometry.outer))
                    * sinc(q[1] * geometry.height / (2.0 * geometry.outer))
                    * (q[0] * geometry.separation / (2.0 * geometry.outer)).cos();
                let actual = pupil.normalized_amplitude(rotate(q, 0.37));
                close(actual, C64::from_re(expected), 3e-12 * geometry.peak());
                assert_eq!(actual.im, 0.0);
                assert!(
                    (pupil.density(rotate(q, 0.37)) - expected * expected).abs()
                        < 6e-12 * geometry.peak().powi(2)
                );
            }
            let first_dark = rotate([PI * geometry.outer / geometry.separation, 0.0], 0.37);
            assert!(pupil.normalized_amplitude(first_dark).abs() < 2e-14 * geometry.peak());
            let reversed_lobe = rotate([TAU * geometry.outer / geometry.separation, 0.0], 0.37);
            assert!(pupil.normalized_amplitude(reversed_lobe).re < 0.0);
        }
    }

    // Independent composite Simpson integration in the actual translated
    // coordinate. It does not use the native Fresnel translation identity.
    fn axis_average(center: f64, half: f64, frequency: f64, phase: f64) -> C64 {
        const INTERVALS: usize = 8192;
        let mut sum = C64::ZERO;
        for index in 0..=INTERVALS {
            let coordinate = center + half * (-1.0 + 2.0 * index as f64 / INTERVALS as f64);
            let argument = phase * coordinate * coordinate - frequency * coordinate;
            let (sin, cos) = argument.sin_cos();
            let weight = if index == 0 || index == INTERVALS {
                1.0
            } else if index % 2 == 0 {
                2.0
            } else {
                4.0
            };
            sum = sum + C64::new(cos, sin).scale(weight);
        }
        sum.scale(1.0 / (3.0 * INTERVALS as f64))
    }

    /// G1: shifted Fresnel frequencies and both common/relative phases match
    /// a direct integral over the two physical openings, on either side of focus.
    #[test]
    fn defocused_pair_matches_independent_translated_fresnel_integrals() {
        for aspect in [1.0, 16.0, 256.0] {
            let geometry = Geometry::new(2.0, 2.0 / aspect, 5.0, 0.37);
            let pupil = geometry.pupil();
            let half = [
                geometry.width / (2.0 * geometry.outer),
                geometry.height / (2.0 * geometry.outer),
            ];
            let offset = geometry.separation / (2.0 * geometry.outer);
            for phase in [-16.0, -6.0, 1e-6, 7.0, 16.0] {
                for q in [[0.0; 2], [2e-7, -3e-7], [0.4, -1.1], [20.0, 1.0]] {
                    let expected = (axis_average(offset, half[0], q[0], phase)
                        + axis_average(-offset, half[0], q[0], phase))
                        * axis_average(0.0, half[1], q[1], phase);
                    let expected = expected.scale(0.5 * geometry.peak());
                    let actual = pupil
                        .normalized_defocused_amplitude(rotate(q, 0.37), phase)
                        .unwrap();
                    close(actual, expected, 2e-10 * geometry.peak());
                    let density = pupil.defocused_density(rotate(q, 0.37), phase).unwrap();
                    assert_eq!(density, actual.norm_sq());
                }
            }
            for phase in [f64::NAN, f64::INFINITY, 16.01, -16.01] {
                assert!(
                    pupil
                        .normalized_defocused_amplitude([0.0; 2], phase)
                        .is_none()
                );
            }
            assert!(pupil.defocused_density([f64::MAX; 2], 1.0).is_none());
        }
    }

    /// G3: closing the dark gap yields the same coherent pupil as the one
    /// rectangle containing both openings, including its Fresnel phase.
    #[test]
    fn touching_pair_equals_the_combined_rectangle() {
        for height in [2.0, 0.125] {
            let geometry = Geometry::new(2.0, height, 2.0, 0.37);
            let pupil = geometry.pupil();
            let axes = [geometry.axes[0].map(|v| 2.0 * v), geometry.axes[1]];
            let vertices = corners(axes);
            let combined = PolygonPupil::rectangular(geometry.outer, &vertices, axes);
            for phase in [-16.0, -3.0, 0.0, 1e-6, 8.0, 16.0] {
                for q in [
                    [0.0; 2],
                    [2e-7, -3e-7],
                    [0.4, -1.1],
                    [20.0, 1.0],
                    [100.0, -7.0],
                ] {
                    close(
                        pupil.normalized_defocused_amplitude(q, phase).unwrap(),
                        combined.normalized_defocused_amplitude(q, phase).unwrap(),
                        2e-10 * geometry.peak(),
                    );
                }
            }
        }
    }

    fn sample(index: u32, count: u32) -> [f64; 2] {
        [
            (f64::from(index) + 0.5) / f64::from(count),
            f64::from(index.reverse_bits()) / 4_294_967_296.0 + 0.5 / f64::from(count),
        ]
    }

    fn grating(openings: u16, width: f64, height: f64, pitch: f64, angle: f64) -> Geometry {
        let mut geometry = Geometry::new(width, height, pitch, angle);
        geometry.outer = (0.5 * (width + f64::from(openings - 1) * pitch)).hypot(0.5 * height);
        geometry
    }

    /// G1: independently enumerate the physical opening phases. Integer
    /// diffraction orders have a removable limit with a real, signed amplitude.
    #[test]
    fn grating_factor_matches_finite_phasor_sum_at_orders_and_nulls() {
        for count in [1, 2, 3, 7, 32, 255, 256] {
            for phase in [0.0, 1e-12, -0.23, PI, -2.0 * PI, PI + 1e-9, 27.1] {
                let expected: f64 = (0..count)
                    .map(|j| ((2.0 * f64::from(j) - f64::from(count - 1)) * phase).cos())
                    .sum();
                let actual = array_factor([phase, 0.0], [1.0, 0.0], count);
                assert!(
                    (actual - expected).abs() < 2e-12 * f64::from(count),
                    "N={count}, phase={phase}: {actual} != {expected}"
                );
            }
            for order in [0.0, 1.0, 2.0, -1.0] {
                let first_minimum = (order + 1.0 / f64::from(count)) * PI;
                if count > 1 {
                    assert!(
                        array_factor([first_minimum, 0.0], [1.0, 0.0], count).abs()
                            < 2e-12 * f64::from(count)
                    );
                }
            }
        }
    }

    /// G1: finite slit width modulates every grating order; total OPEN area
    /// fixes the on-axis amplitude, independently of opaque gap area.
    #[test]
    fn grating_matches_signed_fraunhofer_pattern_and_resolving_width() {
        for count in [3, 8, 31, 256] {
            let geometry = grating(count, 1.0, 2.0, 3.0, 0.37);
            let pupil = geometry.array(count);
            let peak = (f64::from(count) * geometry.width * geometry.height).sqrt()
                / (TAU * geometry.outer);
            for period in [0.0, 0.17, 0.73, 1.0, -2.0, 3.25] {
                let q = [TAU * geometry.outer * period / geometry.separation, 0.7];
                let phases: f64 = (0..count)
                    .map(|j| (TAU * period * (f64::from(j) - 0.5 * f64::from(count - 1))).cos())
                    .sum();
                let expected = peak * phases / f64::from(count)
                    * sinc(PI * period * geometry.width / geometry.separation)
                    * sinc(0.7 * geometry.height / (2.0 * geometry.outer));
                close(
                    pupil.normalized_amplitude(rotate(q, 0.37)),
                    C64::from_re(expected),
                    3e-12 * peak,
                );
            }
            let zero = [
                TAU * geometry.outer / (f64::from(count) * geometry.separation),
                0.0,
            ];
            assert!(pupil.normalized_amplitude(rotate(zero, 0.37)).abs() < 2e-12 * peak);
            close(
                pupil.normalized_amplitude([0.0; 2]),
                C64::from_re(peak),
                1e-14 * peak,
            );
        }
    }

    /// G1/G3: integrate the actual finite openings, preserving full absolute
    /// phase and both signs of defocus before taking the squared field.
    #[test]
    fn grating_fresnel_matches_independent_opening_integrals() {
        for count in [3, 7] {
            let geometry = grating(count, 2.0, 0.5, 3.0, 0.37);
            let pupil = geometry.array(count);
            let peak = (f64::from(count) * geometry.width * geometry.height).sqrt()
                / (TAU * geometry.outer);
            for phase in [-16.0, 1e-6, 7.0, 16.0] {
                for q in [[0.0; 2], [1.3, -0.7], [24.0, 8.0]] {
                    let mut x = C64::ZERO;
                    for j in 0..count {
                        let center = (f64::from(j) - 0.5 * f64::from(count - 1))
                            * geometry.separation
                            / geometry.outer;
                        x = x + axis_average(
                            center,
                            geometry.width / (2.0 * geometry.outer),
                            q[0],
                            phase,
                        );
                    }
                    let expected = x * axis_average(
                        0.0,
                        geometry.height / (2.0 * geometry.outer),
                        q[1],
                        phase,
                    );
                    close(
                        pupil
                            .normalized_defocused_amplitude(rotate(q, 0.37), phase)
                            .unwrap(),
                        expected.scale(peak / f64::from(count)),
                        3e-10 * peak,
                    );
                }
            }
        }
    }

    /// G3: with no opaque gaps, the finite grating is exactly one rectangle,
    /// including odd counts, high counts and nonzero Fresnel phase.
    #[test]
    fn touching_grating_matches_single_rectangle() {
        for count in [1, 3, 8, 256] {
            let geometry = grating(count, 2.0, 1.0, 2.0, 0.37);
            let pupil = geometry.array(count);
            let axes = [
                geometry.axes[0].map(|v| f64::from(count) * v),
                geometry.axes[1],
            ];
            let vertices = corners(axes);
            let combined = PolygonPupil::rectangular(geometry.outer, &vertices, axes);
            let peak = (f64::from(count) * geometry.width * geometry.height).sqrt()
                / (TAU * geometry.outer);
            for phase in [-16.0, 0.0, 4.0] {
                for q in [[0.0; 2], [0.3, -0.7], [30.0, 4.0]] {
                    close(
                        pupil.normalized_defocused_amplitude(q, phase).unwrap(),
                        combined.normalized_defocused_amplitude(q, phase).unwrap(),
                        5e-10 * peak,
                    );
                }
            }
        }
    }

    /// G3: complete optical power, including high grating orders, is retained
    /// by the original full-support spectral and source-depth proposals.
    #[test]
    fn grating_proposals_conserve_power_without_order_cutoffs() {
        const SAMPLES: u32 = 32_768;
        let geometry = grating(7, 2.0, 1.0, 3.0, 0.37);
        let pupil = geometry.array(7);
        let mut totals = [0.0; 3];
        for index in 0..SAMPLES {
            let (q, proposal) = pupil.sample(sample(index, SAMPLES));
            for (sum, ratio) in totals.iter_mut().zip([550.0 / 400.0, 1.0, 550.0 / 700.0]) {
                *sum += ratio * ratio * pupil.density(q.map(|v| ratio * v)) / proposal;
            }
        }
        for total in totals {
            let power = total / f64::from(SAMPLES);
            assert!((power - 1.0).abs() < 0.015, "grating power {power}");
        }
        // Endpoints and interior draws must remain in the full proposal even
        // for the largest admitted array and the widest source-depth interval.
        let geometry = grating(256, 1.0, 2.0, 3.0, 0.37);
        let pupil = geometry.array(256);
        for draw in [[0.0; 2], [0.5; 2], [1.0_f64.next_down(); 2]] {
            for phase in [-16.0, 0.0, 16.0] {
                let (q, density) = pupil.sample_defocused(draw, phase);
                assert!(q.into_iter().all(f64::is_finite));
                assert!(density.is_finite() && density > 0.0);
                assert!(pupil.defocused_density(q, phase).unwrap().is_finite());
            }
            let (q, density) = pupil.sample_phase_range(draw, 16.0);
            assert!(q.into_iter().all(f64::is_finite));
            assert!(density.is_finite() && density > 0.0);
        }
    }

    /// G3: complete sinc/fringe tails retain unit power across wavelengths and
    /// slit aspect ratios, without integrating a cropped detector window.
    #[test]
    fn paired_proposals_conserve_spectral_and_defocused_power() {
        const SAMPLES: u32 = 32_768;
        for aspect in [1.0, 128.0] {
            let geometry = Geometry::new(2.0, 2.0 / aspect, 5.0, 0.37);
            let pupil = geometry.pupil();
            let mut totals = [0.0; 3];
            for index in 0..SAMPLES {
                let (q, proposal) = pupil.sample(sample(index, SAMPLES));
                for (sum, ratio) in totals.iter_mut().zip([550.0 / 400.0, 1.0, 550.0 / 700.0]) {
                    *sum += ratio * ratio * pupil.density(q.map(|v| ratio * v)) / proposal;
                }
            }
            for total in totals {
                let power = total / f64::from(SAMPLES);
                assert!((power - 1.0).abs() < 0.01, "aspect {aspect}: power {power}");
            }
        }

        let geometry = Geometry::new(2.0, 0.5, 3.0, 0.37);
        let pupil = geometry.pupil();
        const DEFOCUSED_SAMPLES: u32 = 8192;
        let mut totals = [0.0; 4];
        for index in 0..DEFOCUSED_SAMPLES {
            let draw = sample(index, DEFOCUSED_SAMPLES);
            let (q, proposal) = pupil.sample_defocused(draw, 6.0);
            totals[0] += pupil.defocused_density(q, 6.0).unwrap() / proposal;
            let (q, proposal) = pupil.sample_phase_range(draw, 10.0);
            for (sum, (phase, ratio)) in totals[1..].iter_mut().zip([
                (0.0, 550.0 / 400.0),
                (6.0, 1.0),
                (-6.0, 550.0 / 700.0),
            ]) {
                *sum += ratio
                    * ratio
                    * pupil
                        .defocused_density(q.map(|v| ratio * v), phase * ratio)
                        .unwrap()
                    / proposal;
            }
        }
        for total in totals {
            let power = total / f64::from(DEFOCUSED_SAMPLES);
            assert!((power - 1.0).abs() < 0.015, "propagated pair power {power}");
        }
        for draw in [[0.0; 2], [0.5; 2], [1.0_f64.next_down(); 2]] {
            let (common_q, common_proposal) = pupil.sample_phase_range(draw, 16.0);
            assert!(common_q.into_iter().all(f64::is_finite));
            assert!(common_proposal.is_finite() && common_proposal > 0.0);
            for phase in [-16.0, 0.0, 16.0] {
                let (q, proposal) = pupil.sample_defocused(draw, phase);
                assert!(q.into_iter().all(f64::is_finite));
                assert!(proposal.is_finite() && proposal > 0.0);
                assert!((pupil.defocused_density(q, phase).unwrap() / proposal).is_finite());
                assert!(
                    (pupil.defocused_density(common_q, phase).unwrap() / common_proposal)
                        .is_finite()
                );
            }
        }
    }
}
