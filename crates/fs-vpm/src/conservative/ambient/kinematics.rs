//! Eulerian samples of the original held source and its current vortex field.

use super::AffineFlow2;
use crate::{ParticleKinematics2, VortexParticle, induced_kinematics};

impl AffineFlow2 {
    /// Exact velocity, Jacobian and counterclockwise-positive curl of this
    /// authored field at a physical point [m]. For `A=[[a,b],[c,-a]]`, the
    /// gradient is A and vorticity is c-b; local fluid spin is HALF this curl.
    ///
    /// A point must be finite. Unrepresentable needed coordinate differences,
    /// field contributions or sums return `None`; cancellation between
    /// unrepresentable terms is not promised. Uniform fields do not depend on
    /// their chart origin. This constant-work query allocates nothing.
    #[must_use]
    pub fn kinematics_at(self, point_m: [f64; 2]) -> Option<ParticleKinematics2> {
        if !point_m.into_iter().all(f64::is_finite) {
            return None;
        }
        let [a, b, c] = self.gradient;
        let gradient = [[a, b], [c, -a]];
        let mut velocity = self.velocity;
        for j in 0..2 {
            // An unused chart coordinate must not make a uniform flow or
            // single-column shear singular at a representable physical point.
            if gradient[0][j] != 0.0 || gradient[1][j] != 0.0 {
                let offset = point_m[j] - self.origin[j];
                if !offset.is_finite() {
                    return None;
                }
                for i in 0..2 {
                    velocity[i] += gradient[i][j] * offset;
                }
            }
        }
        finite(ParticleKinematics2 {
            velocity,
            gradient,
            vorticity: c - b,
        })
    }

    /// Complete prescribed-background plus particle-induced flow at a point.
    ///
    /// The original algebraic-core kernel evaluates the CURRENT particle
    /// positions and circulations in their original source order. This source
    /// is added once after that sum, including its analytic curl; curl is not
    /// reconstructed from rounded gradients or from PSE `Gamma/area`.
    /// A zero background preserves the original induced sample exactly.
    ///
    /// The background must describe an additional maintained field; supplying
    /// the same vorticity both here and in the particles would count it twice.
    /// This is a field query, not a body-force law or a finite-fluid reaction.
    /// Invalid/singular/unrepresentable source fields or total sums return
    /// `None`. Work is O(N), storage O(1), with no heap allocation or mutation.
    /// As for `induced_kinematics`, callers bound bulk work; no `Cx` receipt or
    /// cancellation checkpoint is added to a single query.
    #[must_use]
    pub fn kinematics_with_particles_at(
        self,
        particles: &[VortexParticle],
        point_m: [f64; 2],
        core_m: f64,
    ) -> Option<ParticleKinematics2> {
        let mut sample = induced_kinematics(particles, point_m, core_m)?;
        if self.is_zero() {
            return Some(sample);
        }
        let background = self.kinematics_at(point_m)?;
        for i in 0..2 {
            sample.velocity[i] += background.velocity[i];
            for j in 0..2 {
                sample.gradient[i][j] += background.gradient[i][j];
            }
        }
        sample.vorticity += background.vorticity;
        finite(sample)
    }
}

fn finite(sample: ParticleKinematics2) -> Option<ParticleKinematics2> {
    sample
        .velocity
        .iter()
        .chain(sample.gradient.iter().flatten())
        .chain(core::iter::once(&sample.vorticity))
        .all(|value| value.is_finite())
        .then_some(sample)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    fn flow(origin: [f64; 2], velocity: [f64; 2], gradient: [f64; 3]) -> AffineFlow2 {
        // These finite fixtures are precisely the source data admitted by new.
        AffineFlow2 {
            origin,
            velocity,
            gradient,
        }
    }

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 8e-15 * (1.0 + expected.abs()),
            "{actual:e} != {expected:e}"
        );
    }

    fn same(actual: ParticleKinematics2, expected: ParticleKinematics2) {
        for i in 0..2 {
            near(actual.velocity[i], expected.velocity[i]);
            for j in 0..2 {
                near(actual.gradient[i][j], expected.gradient[i][j]);
            }
        }
        near(actual.vorticity, expected.vorticity);
    }

    #[test]
    fn simple_shear_rotation_and_extension_have_their_analytic_curl() {
        let origin = [2.0, 3.0];
        let point = [5.0, 4.0];
        let velocity = [0.5, -0.2];
        same(
            flow(origin, velocity, [0.0, 0.4, 0.0])
                .kinematics_at(point)
                .unwrap(),
            ParticleKinematics2 {
                velocity: [0.9, -0.2],
                gradient: [[0.0, 0.4], [0.0, 0.0]],
                vorticity: -0.4,
            },
        );
        let rotation = flow(origin, velocity, [0.0, -0.7, 0.7])
            .kinematics_at(point)
            .unwrap();
        same(
            rotation,
            ParticleKinematics2 {
                velocity: [-0.2, 1.9],
                gradient: [[0.0, -0.7], [0.7, 0.0]],
                vorticity: 1.4,
            },
        );
        near(0.5 * rotation.vorticity, 0.7);
        same(
            flow(origin, velocity, [0.3, 0.0, 0.0])
                .kinematics_at(point)
                .unwrap(),
            ParticleKinematics2 {
                velocity: [1.4, -0.5],
                gradient: [[0.3, 0.0], [0.0, -0.3]],
                vorticity: 0.0,
            },
        );
    }

    #[test]
    fn chart_origin_shift_and_world_rotation_preserve_the_physical_field() {
        let sample = flow([0.2, -0.3], [0.6, -0.1], [0.2, 0.7, -0.3])
            .kinematics_at([0.8, 0.4])
            .unwrap();
        // Shift the source chart by (0.4,-0.8), and U by A times that shift.
        same(
            flow([0.6, -1.1], [0.12, -0.06], [0.2, 0.7, -0.3])
                .kinematics_at([0.8, 0.4])
                .unwrap(),
            sample,
        );
        let rotated = flow([0.3, 0.2], [0.1, 0.6], [-0.2, 0.3, -0.7])
            .kinematics_at([-0.4, 0.8])
            .unwrap();
        same(
            rotated,
            ParticleKinematics2 {
                velocity: [-sample.velocity[1], sample.velocity[0]],
                gradient: [
                    [sample.gradient[1][1], -sample.gradient[1][0]],
                    [-sample.gradient[0][1], sample.gradient[0][0]],
                ],
                vorticity: sample.vorticity,
            },
        );
    }

    #[test]
    fn single_vortex_and_shear_match_independent_superposition() {
        let particles = [VortexParticle::new([0.0; 2], 2.0 * PI)];
        let saved = particles;
        let source = flow([0.0; 2], [0.1, -0.2], [0.0, 0.25, 0.0]);
        same(
            source
                .kinematics_with_particles_at(&particles, [1.0, 0.0], 1.0)
                .unwrap(),
            ParticleKinematics2 {
                velocity: [0.1, 0.3],
                gradient: [[0.0, -0.25], [0.0, 0.0]],
                vorticity: 0.25,
            },
        );
        assert_eq!(particles, saved);
    }

    #[test]
    fn zero_background_preserves_original_particle_order_and_bits() {
        let particles = [
            VortexParticle::new([0.1, -0.2], 0.7),
            VortexParticle::new([-0.3, 0.5], -0.2),
            VortexParticle::new([0.4, 0.8], 0.3),
        ];
        let point = [0.2, 0.3];
        let induced = induced_kinematics(&particles, point, 0.17).unwrap();
        let sample = flow([1e308, -1e308], [0.0; 2], [0.0; 3])
            .kinematics_with_particles_at(&particles, point, 0.17)
            .unwrap();
        assert_eq!(sample, induced);
        for i in 0..2 {
            assert_eq!(sample.velocity[i].to_bits(), induced.velocity[i].to_bits());
            for j in 0..2 {
                assert_eq!(
                    sample.gradient[i][j].to_bits(),
                    induced.gradient[i][j].to_bits()
                );
            }
        }
        assert_eq!(sample.vorticity.to_bits(), induced.vorticity.to_bits());
    }

    #[test]
    fn composite_curl_retains_the_resolved_particle_tail() {
        let particles = [VortexParticle::new([0.0; 2], 1e300)];
        let sample = flow([0.0; 2], [0.0; 2], [0.0, 1.0, 1.0])
            .kinematics_with_particles_at(&particles, [1.0, 0.0], 1e-200)
            .unwrap();
        let expected = 1e-100 / PI;
        assert_eq!(sample.gradient[1][0] - sample.gradient[0][1], 0.0);
        assert!((sample.vorticity / expected - 1.0).abs() < 4e-15);
    }

    #[test]
    fn uniform_and_shear_sources_ignore_unused_chart_offsets() {
        let point = [-1e308, 0.5];
        let uniform = flow([1e308, 1e308], [0.2, -0.3], [0.0; 3]);
        assert_eq!(uniform.kinematics_at(point).unwrap().velocity, [0.2, -0.3]);
        let shear = flow([1e308, 0.0], [0.2, -0.3], [0.0, 0.4, 0.0]);
        assert_eq!(shear.kinematics_at(point).unwrap().velocity, [0.4, -0.3]);
    }

    #[test]
    fn singular_invalid_and_unrepresentable_sources_and_sums_refuse() {
        let uniform = flow([0.0; 2], [0.2, -0.3], [0.0; 3]);
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(uniform.kinematics_at([bad, 0.0]).is_none());
            assert!(uniform.kinematics_at([0.0, bad]).is_none());
        }
        let particles = [VortexParticle::new([0.0; 2], 1.0)];
        assert!(
            uniform
                .kinematics_with_particles_at(&particles, [0.0; 2], 0.0)
                .is_none()
        );
        for core in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                uniform
                    .kinematics_with_particles_at(&[], [0.0; 2], core)
                    .is_none()
            );
        }
        assert!(
            flow([1e308, 0.0], [0.0; 2], [1.0, 0.0, 0.0])
                .kinematics_at([-1e308, 0.0])
                .is_none()
        );
        assert!(
            flow([0.0; 2], [0.0; 2], [1e308, 0.0, 0.0])
                .kinematics_at([2.0, 0.0])
                .is_none()
        );
        assert!(
            flow([0.0; 2], [0.0; 2], [0.0, -1e308, 1e308])
                .kinematics_at([0.0; 2])
                .is_none()
        );
        let huge = [VortexParticle::new([0.0; 2], 1e308)];
        assert!(
            flow([0.0; 2], [f64::MAX, 0.0], [0.0; 3])
                .kinematics_with_particles_at(&huge, [0.0, -1.0], 1.0)
                .is_none()
        );
    }
}
