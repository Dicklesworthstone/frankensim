//! Coherent diffraction through two identical, disjoint rectangular openings.
//!
//! Each component is evaluated in its own circumradius, so a long baseline
//! does not turn a well-conditioned rectangle into a vanishing polygon. For
//! component/outer radius ratio r and normalized half-offset c, translation
//! of the original Fresnel integral gives
//! `a(q,phi) = r/sqrt(2) sum_(s=+/-1) exp(i*phi*|c|²-i*s*q.c)
//! a_single(r*(q-2*s*phi*c),phi*r²)`.
//! The normalization uses the sum of the two transmitting areas. The dark gap
//! transmits nothing, and the two fields interfere before power is measured.

use super::polygon::PolygonPupil;
use fs_math::{c64::C64, det};

/// Borrow an admitted pair. Vertices and half-axes describe the centered
/// component; the two physical copies are translated by +/-half_offset_m.
pub(super) struct RectangularPairPupil<'a> {
    component: PolygonPupil<'a>,
    radius_ratio: f64,
    half_offset: [f64; 2],
}

impl<'a> RectangularPairPupil<'a> {
    pub(super) fn new(
        outer_radius_m: f64,
        component_radius_m: f64,
        vertices: &'a [[f64; 2]],
        half_axes_m: [[f64; 2]; 2],
        half_offset_m: [f64; 2],
    ) -> Self {
        Self {
            component: PolygonPupil::rectangular(component_radius_m, vertices, half_axes_m),
            radius_ratio: component_radius_m / outer_radius_m,
            half_offset: half_offset_m.map(|value| value / outer_radius_m),
        }
    }

    /// Unit-integral density per dqx dqy, including the interference fringes.
    pub(super) fn density(&self, q: [f64; 2]) -> f64 {
        self.normalized_amplitude(q).norm_sq()
    }

    /// The exact focused sinc envelope times its signed two-opening fringe.
    pub(super) fn normalized_amplitude(&self, q: [f64; 2]) -> C64 {
        self.component
            .normalized_amplitude(q.map(|value| self.radius_ratio * value))
            .scale(
                core::f64::consts::SQRT_2 * self.radius_ratio * det::cos(dot(q, self.half_offset)),
            )
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

    /// Reuse the component's complete Cauchy proposal. The focused fringe
    /// multiplies the component PSF by 2*cos(q.c)^2, so it adds at most a factor
    /// of two to that proposal's original importance weights.
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
        let (u, sign) = if sample[0] < 0.5 {
            (2.0 * sample[0], -1.0)
        } else {
            (2.0 * sample[0] - 1.0, 1.0)
        };
        let local_phase = phase * self.radius_ratio * self.radius_ratio;
        let (local, _) = self.component.sample_defocused([u, sample[1]], local_phase);
        let q = [
            local[0] / self.radius_ratio + sign * 2.0 * phase * self.half_offset[0],
            local[1] / self.radius_ratio + sign * 2.0 * phase * self.half_offset[1],
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
                * det::hypot(self.half_offset[0], self.half_offset[1]);
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
        let shift = self.half_offset.map(|value| 2.0 * phase * value);
        0.5 * (self.widened_proposal_density([q[0] - shift[0], q[1] - shift[1]], widening)
            + self.widened_proposal_density([q[0] + shift[0], q[1] + shift[1]], widening))
    }

    fn widened_proposal_density(&self, q: [f64; 2], widening: f64) -> f64 {
        let scale = self.radius_ratio / widening;
        self.component
            .proposal_density(q.map(|value| scale * value))
            * scale
            * scale
    }
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

        fn pupil(&self) -> RectangularPairPupil<'_> {
            RectangularPairPupil::new(
                self.outer,
                self.component,
                &self.vertices,
                self.axes,
                self.offset,
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
