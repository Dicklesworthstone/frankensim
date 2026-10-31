//! Coherent diffraction of a uniformly transmitting concentric annulus.
//!
//! For inner/outer ratio e, delta=1-e² and the full complex disk amplitude A,
//! B(q,phi)=A(q,phi)-e²*A(e*q,e²*phi). The unit-integral image density is
//! |B|²/(4*pi*delta). Both amplitudes include their different exp(i*phi/2)
//! factors. Subtracting intensities, or omitting the relative phase, is wrong.
//! This follows directly by splitting the circular Fresnel pupil integral.
//!
//! The admitted e<=0.95 keeps delta>=0.0975; a zero-width ring has no finite
//! collection area and its J0² limit is not a normalizable 2D image PSF.
//! Axis and small focused-radius limits avoid removable cancellation.

use super::circular::{CircularDefocus, jinc};
use fs_math::det;

const PI: f64 = core::f64::consts::PI;

pub(super) struct AnnularPupil {
    ratio: f64,
    clear_fraction: f64,
    phase: f64,
    outer: CircularDefocus,
    inner: CircularDefocus,
}

impl AnnularPupil {
    pub(super) fn new(ratio: f64, phase: f64) -> Self {
        debug_assert!(ratio.is_finite() && (0.0..=0.95).contains(&ratio));
        Self {
            ratio,
            clear_fraction: (1.0 - ratio) * (1.0 + ratio),
            phase,
            outer: CircularDefocus::new(phase),
            inner: CircularDefocus::new(ratio * ratio * phase),
        }
    }

    pub(super) fn density(&self, radius: f64) -> f64 {
        if self.ratio == 0.0 {
            return self.outer.density(radius);
        }
        if radius == 0.0 {
            let argument = 0.5 * self.phase * self.clear_fraction;
            let sinc = if argument.abs() < 1.0e-4 {
                1.0 + argument * argument * (-1.0 / 6.0 + argument * argument / 120.0)
            } else {
                det::sin(argument) / argument
            };
            return self.clear_fraction * sinc * sinc / (4.0 * PI);
        }
        let ratio_squared = self.ratio * self.ratio;
        let power = if self.phase == 0.0 {
            let amplitude = if radius < 1.0e-4 {
                let square = radius * radius;
                self.clear_fraction
                    * (1.0 - (1.0 + ratio_squared) * square / 8.0
                        + (1.0 + ratio_squared + ratio_squared * ratio_squared) * square * square
                            / 192.0)
            } else {
                jinc(radius) - ratio_squared * jinc(self.ratio * radius)
            };
            amplitude * amplitude
        } else {
            let amplitude = self.outer.amplitude(radius)
                - self
                    .inner
                    .amplitude(self.ratio * radius)
                    .scale(ratio_squared);
            amplitude.re * amplitude.re + amplitude.im * amplitude.im
        };
        power / (4.0 * PI * self.clear_fraction)
    }
}

/// One radial draw from an exact half/half Cauchy mixture. The narrow scale
/// retains the central lobe of thin rings, whose phase width is phi*(1-e²),
/// while the broad scale covers their ring groups and geometric defocus.
/// The common azimuth is supplied independently by the caller.
pub(super) fn sample_radius(ratio: f64, phase: f64, sample: f64) -> (f64, f64) {
    let broad = 2.0_f64.max(2.0 * phase.abs()).max(1.0 / (1.0 - ratio));
    let (scale, u) = if sample < 0.5 {
        (2.0, 2.0 * sample)
    } else {
        (broad, 2.0 * sample - 1.0)
    };
    let radius = scale * det::sqrt(u * (2.0 - u)) / (1.0 - u);
    let density = 0.5 * (cauchy_density(radius, 2.0) + cauchy_density(radius, broad));
    (radius, density)
}

fn cauchy_density(radius: f64, scale: f64) -> f64 {
    let relative = radius / scale;
    let square = 1.0 + relative * relative;
    1.0 / (2.0 * PI * scale * scale * square * det::sqrt(square))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coherent_annulus_matches_independent_opening_integrals() {
        // Direct adaptive integration over the opening [e,1], with SciPy J0
        // and QUADPACK epsabs=2e-14/epsrel=2e-13. This does not subtract disks
        // or reuse the production Bauer expansion. Cases straddle both disk
        // recurrence transitions, retain negative phase, and reach the tail.
        let references = [
            (0.2, 0.0, 0.0001, 0.07639437248548439),
            (0.6, 0.0, 3.8317059702075125, 0.0035567385140390126),
            (0.95, 0.0, 128.001, 1.1294416980161745e-9),
            (0.2, 1.0, 0.499999, 0.06623095012949863),
            (0.6, 1.0, 0.500001, 0.04516601365568397),
            (0.95, 16.0, 1.0e-8, 0.0063075596032000584),
            (0.6, 16.0, 0.499999, 0.0015031640030211237),
            (0.95, 16.0, 0.500001, 0.0055902767392891984),
            (0.2, -6.0, 4.0, 0.0038466516273179133),
            (0.6, 2.0 * PI, 8.0, 0.0010898403896214158),
            (0.95, 16.0, 16.0, 9.651944292349691e-5),
            (0.6, 16.0, 31.0, 0.0002239702942235605),
            (0.95, 16.0, 79.0, 9.54620543531154e-6),
            (0.2, 16.0, 127.999, 8.771064753478863e-8),
            (0.6, 16.0, 128.001, 3.297754127958416e-7),
            (0.95, -16.0, 128.001, 2.0840100927527273e-6),
            (0.6, 16.0, 1000.0, 1.5206265165147768e-10),
            (0.95, 16.0, 1000.0, 4.887078426521367e-10),
        ];
        for (ratio, phase, radius, expected) in references {
            let actual = AnnularPupil::new(ratio, phase).density(radius);
            assert!(
                (actual - expected).abs() < 5.0e-10 * expected,
                "e={ratio}, phi={phase}, q={radius}: {actual} != {expected}"
            );
        }
    }

    #[test]
    fn annular_axis_phase_and_clear_disk_limits_are_physical() {
        for ratio in [0.0_f64, 0.2, 0.6, 0.95] {
            let fraction = (1.0 - ratio) * (1.0 + ratio);
            for phase in [0.0_f64, 1.0e-8, 2.0 * PI, 16.0] {
                let pupil = AnnularPupil::new(ratio, phase);
                let argument = 0.5 * phase * fraction;
                let sinc = if argument == 0.0 {
                    1.0
                } else {
                    argument.sin() / argument
                };
                let expected = fraction * sinc * sinc / (4.0 * PI);
                assert!((pupil.density(0.0) - expected).abs() < 2.0e-16);
                for radius in [1.0e-150, 0.000099, 0.000101, 0.5, 8.0, 128.0, 1.0e4, 1.0e17] {
                    let density = pupil.density(radius);
                    assert!(density.is_finite() && density >= 0.0);
                    assert_eq!(
                        density.to_bits(),
                        AnnularPupil::new(ratio, -phase).density(radius).to_bits()
                    );
                    if ratio == 0.0 {
                        assert_eq!(
                            density.to_bits(),
                            CircularDefocus::new(phase).density(radius).to_bits()
                        );
                    }
                }
            }
        }
        let ratio = 0.6;
        let axial_null = AnnularPupil::new(ratio, 2.0 * PI / ((1.0 - ratio) * (1.0 + ratio)));
        assert!(axial_null.density(0.0) < 1.0e-28);
        // A disk's first zero is filled by coherent light from the inner edge.
        assert!(AnnularPupil::new(ratio, 0.0).density(3.8317059702075125) > 0.003);
        // The new complex composition helper restores the original phase,
        // including its sign; the old circle density keeps its own bit path.
        for radius in [0.0, 0.5, 8.0, 128.001, 1.0e4] {
            let positive = CircularDefocus::new(6.0).amplitude(radius);
            let negative = CircularDefocus::new(-6.0).amplitude(radius);
            assert_eq!(positive.re.to_bits(), negative.re.to_bits());
            assert_eq!(positive.im.to_bits(), (-negative.im).to_bits());
        }
    }

    #[test]
    fn annular_packet_conserves_power_with_full_mixture_and_complete_tails() {
        const SAMPLES: u32 = 65_536;
        let ratios = [550.0 / 400.0, 1.0, 550.0 / 700.0];
        for obscuration in [0.2, 0.6, 0.95] {
            for phase in [0.0, 6.0, 10.0] {
                let pupils = ratios.map(|ratio| AnnularPupil::new(obscuration, phase * ratio));
                let mut power = [0.0; 3];
                for sample in 0..SAMPLES {
                    let u = (f64::from(sample) + 0.5) / f64::from(SAMPLES);
                    let (radius, proposal) = sample_radius(obscuration, phase, u);
                    for lane in 0..3 {
                        power[lane] += ratios[lane]
                            * ratios[lane]
                            * pupils[lane].density(ratios[lane] * radius)
                            / proposal
                            / f64::from(SAMPLES);
                    }
                }
                for value in power {
                    assert!(
                        (value - 1.0).abs() < 0.006,
                        "e={obscuration}, phi={phase}: {power:?}"
                    );
                }
                for u in [
                    0.0,
                    0.5 - 0.25 * f64::EPSILON,
                    0.5,
                    0.999,
                    1.0 - 0.5 * f64::EPSILON,
                ] {
                    let (radius, proposal) = sample_radius(obscuration, phase, u);
                    assert!(radius.is_finite() && proposal.is_finite() && proposal > 0.0);
                    for lane in 0..3 {
                        let weight = ratios[lane]
                            * ratios[lane]
                            * pupils[lane].density(ratios[lane] * radius)
                            / proposal;
                        assert!(
                            weight.is_finite() && (0.0..256.0).contains(&weight),
                            "unbounded tail weight {weight}"
                        );
                    }
                }
            }
        }
    }
}
