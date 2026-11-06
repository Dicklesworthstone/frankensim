//! Compact physical stress sources on the existing exterior Green owner.
//!
//! Lighthill's source equation gives p = Q_ij d_i d_j G, with
//! Q_ij = integral(T_ij dV) [N m] and G = exp(ikr)/(4 pi r).
//! All amplitudes use exp(-i omega t). Density is already in the supplied
//! stress; it enters the particle velocity only through i rho omega.
//! See Lighthill (1952), <https://doi.org/10.1098/rspa.1952.0060>.
//!
//! Each tensor represents one compact source cell. Its full trace and physical
//! shear entries survive; a trace-free tensor is not assumed. This does not
//! validate source size, resolve turbulence or omitted stresses, or propagate
//! through mean flow. The same scalar Green owner also admits the existing
//! homogeneous molecular medium's complex k, retaining its physical rho*omega.
//!
//! BEM scattering supplies the response of the retained boundary. Do not also
//! add Curle loading from that same passive boundary: doing so would count its
//! acoustic response twice. Independently prescribed physical boundary forces
//! require an explicit source decomposition. No source feedback is inferred.
//! The stationary-boundary equivalence is discussed by Kaltenbacher & Schoder,
//! DAGA 2020, <https://pub.dega-akustik.de/DAGA_2020/data/articles/000031.pdf>.

use super::*;

/// Prescribed symmetric compact Lighthill stress moment at a physical point.
///
/// The layout matches fs-aeroac's QuadrupolePhasor. Transform the physical
/// stress history at its original emission clock; squaring a velocity phasor
/// while retaining its frequency is not a valid quadratic source spectrum.
/// Any number of these sources, up to the existing incident-source cap, may
/// be combined with plane waves, volume flows and force dipoles through
/// `IncidentField::Sources`. All existing source-clearance, boundary-solve,
/// finite-receiver and outgoing far-field routes retain that coherent field.
#[derive(Clone, Copy, Debug)]
pub struct StressQuadrupole {
    /// Fixed source position [m], in the scatterer's original reference frame.
    pub position_m: [f64; 3],
    /// Qxx, Qyy, Qzz, Qxy, Qxz, Qyz [N m], as complex peak amplitudes.
    /// Shear components are physical tensor entries, without a factor of two.
    pub stress_moment_n_m: [C64; 6],
}

fn tensor_vector(moment: [C64; 6], direction: [f64; 3]) -> [C64; 3] {
    [[0, 3, 4], [3, 1, 5], [4, 5, 2]].map(|row| {
        let mut sum = ComplexSum::default();
        for (index, component) in row.into_iter().zip(direction) {
            sum.add(moment[index].scale(component));
        }
        sum.finish()
    })
}

fn complex_projection(vector: [C64; 3], direction: [f64; 3]) -> C64 {
    let mut sum = ComplexSum::default();
    for (value, component) in vector.into_iter().zip(direction) {
        sum.add(value.scale(component));
    }
    sum.finish()
}

impl StressQuadrupole {
    pub(super) fn validate(self) -> Result<(), HelmholtzError> {
        if self
            .position_m
            .iter()
            .any(|v| !v.is_finite() || v.abs() > 1e9)
            || self.stress_moment_n_m.iter().any(|&v| !finite(v))
        {
            return Err(bad("stress quadrupole needs finite position and physical tensor"));
        }
        Ok(())
    }

    pub(super) fn sample_wave<K: WaveNumber>(
        self,
        wave: K,
        omega_rho: f64,
        point: [f64; 3],
    ) -> Result<IncidentSample, HelmholtzError> {
        let delta = sub(point, self.position_m);
        let radius = norm(delta);
        let coordinate_scale = point
            .into_iter()
            .chain(self.position_m)
            .fold(0.0_f64, |m, v| m.max(v.abs()));
        if !radius.is_finite()
            || radius <= 0.0
            || radius <= 128.0 * f64::EPSILON * coordinate_scale
            || !finite(wave.value().scale(radius))
        {
            return Err(bad("stress-quadrupole receiver coincides with its source or is unresolved"));
        }
        let direction = delta.map(|v| v / radius);
        let inverse_radius = 1.0 / radius;
        let inverse_radius_squared = inverse_radius * inverse_radius;
        let (g, gp, gpp) = green_radial_derivatives(wave, radius);
        let gppp = gpp * wave.radial_factor(radius)
            + gp.scale(2.0 * inverse_radius_squared)
            - g.scale(2.0 * inverse_radius_squared * inverse_radius);
        let a = gpp - gp.scale(inverse_radius);
        let b = gp.scale(inverse_radius);
        let c = gppp - gpp.scale(3.0 * inverse_radius)
            + gp.scale(3.0 * inverse_radius_squared);

        // Split Q = Qzz I + (Q - Qzz I). Exactly isotropic input then uses
        // Laplacian(G)=-k²G directly, even when kr is tiny: subtracting three
        // O(1/r³) Hessian terms would erase its physical O(k²/r) pressure.
        // The remainder is not assumed trace-free and keeps all six entries.
        let isotropic = self.stress_moment_n_m[2];
        let mut remainder = self.stress_moment_n_m;
        remainder[0] = remainder[0] - isotropic;
        remainder[1] = remainder[1] - isotropic;
        remainder[2] = C64::ZERO;
        let vector = tensor_vector(remainder, direction);
        let radial = complex_projection(vector, direction);
        let trace = remainder[0] + remainder[1];
        let isotropic_wave = -(wave.value() * wave.value()) * isotropic;

        let mut pressure = ComplexSum::default();
        pressure.add(a * radial);
        pressure.add(b * trace);
        pressure.add(g * isotropic_wave);
        let velocity_m_s = core::array::from_fn(|axis| {
            // d_l(Q:H(G)) = C (n.Q.n)n_l
            //                  + (A/r)(2(Qn)_l + tr(Q)n_l).
            // All direction derivatives survive, including at pressure nulls.
            let mut gradient = ComplexSum::default();
            gradient.add(c * radial.scale(direction[axis]));
            gradient.add(a.scale(inverse_radius) * vector[axis].scale(2.0));
            gradient.add(a.scale(inverse_radius) * trace.scale(direction[axis]));
            gradient.add(gp * isotropic_wave.scale(direction[axis]));
            let gradient = gradient.finish();
            C64::new(gradient.im / omega_rho, -gradient.re / omega_rho)
        });
        let sample = IncidentSample {
            pressure_pa: pressure.finish(),
            velocity_m_s,
        };
        if !finite(sample.pressure_pa) || sample.velocity_m_s.iter().any(|&v| !finite(v)) {
            return Err(bad("stress-quadrupole pressure or velocity overflow"));
        }
        Ok(sample)
    }

    // The source-position phase belongs to the existing outgoing-field owner.
    // Units are Pa m; this is -k² n.Q.n/(4 pi), not a volume-flow amplitude.
    pub(crate) fn far_field_unphased(self, k: f64, direction: [f64; 3]) -> C64 {
        complex_projection(tensor_vector(self.stress_moment_n_m, direction), direction)
            .scale(-k * k / (4.0 * PI))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::helmholtz::laplace::DampedWaveNumber;

    fn source() -> StressQuadrupole {
        StressQuadrupole {
            position_m: [0.0; 3],
            stress_moment_n_m: [
                C64::new(0.7, 0.1),
                C64::new(-0.3, 0.2),
                C64::new(0.9, -0.4),
                C64::new(0.23, -0.15),
                C64::new(-0.11, 0.3),
                C64::new(0.4, 0.2),
            ],
        }
    }

    fn sample(source: StressQuadrupole, k: f64, point: [f64; 3]) -> IncidentSample {
        IncidentField::Sources(&[IncidentSource::StressQuadrupole(source)])
            .pressure_velocity_at(k, Medium::air(), point)
            .unwrap()
    }

    fn close(actual: C64, expected: C64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
            "{actual:?} != {expected:?}"
        );
    }

    // G1: independent central differences of G, including the factor of two
    // on physical shear entries, converge to the analytic tensor contraction.
    #[test]
    fn stress_pressure_matches_spatial_green_derivatives() {
        let point = [1.3, -0.8, 0.6];
        let source = source();
        let k = 2.3;
        let expected = sample(source, k, point).pressure_pa;
        let axes = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
        let errors = [0.002, 0.001, 0.0005].map(|h| {
            let mut sum = C64::ZERO;
            for ((i, j), moment) in axes.into_iter().zip(source.stress_moment_n_m) {
                let derivative = if i == j {
                    let mut plus = point;
                    let mut minus = point;
                    plus[i] += h;
                    minus[i] -= h;
                    (green(k, norm(plus)) - green(k, norm(point)).scale(2.0)
                        + green(k, norm(minus)))
                        .scale(1.0 / (h * h))
                } else {
                    let mut mixed = C64::ZERO;
                    for a in [-1.0, 1.0] {
                        for b in [-1.0, 1.0] {
                            let mut offset = point;
                            offset[i] += a * h;
                            offset[j] += b * h;
                            mixed = mixed + green(k, norm(offset)).scale(a * b);
                        }
                    }
                    mixed.scale(0.5 / (h * h))
                };
                sum = sum + moment * derivative;
            }
            (sum - expected).abs()
        });
        assert!(errors[1] < 0.28 * errors[0], "{errors:?}");
        assert!(errors[2] < 0.29 * errors[1], "{errors:?}");
        assert!(errors[2] < 2e-6 * expected.abs(), "{errors:?}");
    }

    // G1: physical momentum remains rho*omega, independently of complex k.
    #[test]
    fn all_velocity_components_follow_pressure_for_real_and_absorbing_media() {
        let point = [1.3, -0.8, 0.6];
        let omega_rho = 1100.0;
        for k in [C64::new(2.3, 0.0), C64::new(2.3, 0.4)] {
            let wave = DampedWaveNumber { k, phase_origin_m: 0.0 };
            let field = source().sample_wave(wave, omega_rho, point).unwrap();
            for axis in 0..3 {
                let h = 1e-5;
                let mut plus = point;
                let mut minus = point;
                plus[axis] += h;
                minus[axis] -= h;
                let derivative = (source().sample_wave(wave, omega_rho, plus).unwrap().pressure_pa
                    - source().sample_wave(wave, omega_rho, minus).unwrap().pressure_pa)
                    .scale(0.5 / h);
                close(derivative, field.velocity_m_s[axis] * C64::new(0.0, omega_rho), 5e-9);
            }
        }
    }

    // G1: the trace is a physical source, not a traceless multipole convention.
    #[test]
    fn isotropic_stress_retains_its_small_kr_helmholtz_limit() {
        let point = [1.3, -0.8, 0.6];
        let radius = norm(point);
        let direction = point.map(|x| x / radius);
        let moment = C64::new(0.7, -0.3);
        let source = StressQuadrupole {
            position_m: [0.0; 3],
            stress_moment_n_m: [moment, moment, moment, C64::ZERO, C64::ZERO, C64::ZERO],
        };
        for k in [1e-10, 0.7, 2.3] {
            let field = sample(source, k, point);
            let (g, gp, _) = green_radial_derivatives(k, radius);
            close(field.pressure_pa, (g * moment).scale(-k * k), 5e-15);
            let omega_rho = k * Medium::air().sound_speed * Medium::air().density;
            for axis in 0..3 {
                let gradient = (gp * moment).scale(-k * k * direction[axis]);
                close(field.velocity_m_s[axis],
                    C64::new(gradient.im / omega_rho, -gradient.re / omega_rho), 5e-15);
            }
        }
    }

    // G1: shear dipole-pairs can have zero pressure but nonzero transverse flow.
    #[test]
    fn pressure_null_keeps_transverse_shear_velocity() {
        let mut source = source();
        source.stress_moment_n_m = [C64::ZERO; 6];
        source.stress_moment_n_m[3] = C64::ONE;
        let field = sample(source, 2.3, [0.7, 0.0, 0.0]);
        assert_eq!(field.pressure_pa, C64::ZERO);
        assert_eq!(field.velocity_m_s[0], C64::ZERO);
        assert_eq!(field.velocity_m_s[2], C64::ZERO);
        assert!(field.velocity_m_s[1].abs() > 1e-4);
    }

    // G3: rotate the actual stress tensor, not merely its source position.
    #[test]
    fn rotation_preserves_pressure_and_rotates_velocity() {
        let point = [1.3, -0.8, 0.6];
        let original = source();
        let q = original.stress_moment_n_m;
        let rotated = StressQuadrupole {
            position_m: [0.0; 3],
            stress_moment_n_m: [q[2], q[0], q[1], q[4], q[5], q[3]],
        };
        let a = sample(original, 2.3, point);
        let b = sample(rotated, 2.3, [point[2], point[0], point[1]]);
        close(b.pressure_pa, a.pressure_pa, 1e-14);
        for (actual, expected) in b.velocity_m_s.into_iter()
            .zip([a.velocity_m_s[2], a.velocity_m_s[0], a.velocity_m_s[1]])
        {
            close(actual, expected, 2e-14);
        }
    }

    // G1/G3: the SAME boundary solve must commute with a source-position
    // derivative. This compares a quadrupole with two displaced force sources,
    // including its diffracted near field and complete outgoing directivity.
    #[test]
    fn rigid_body_scattering_retains_the_physical_stress_source() {
        let surface = SpherePanels::icosphere(0.25, 1).unwrap();
        let medium = Medium::air();
        let k = 0.8;
        let h = 1e-4;
        let length = 0.15;
        let position = [0.4, 0.3, 0.9];
        let force = C64::new(0.2, -0.05);
        let quadrupoles = [IncidentSource::StressQuadrupole(StressQuadrupole {
            position_m: position,
            stress_moment_n_m: [
                force.scale(length), C64::ZERO, C64::ZERO, C64::ZERO, C64::ZERO, C64::ZERO,
            ],
        })];
        let plus = [ForceDipole {
            position_m: [position[0] + h, position[1], position[2]],
            force_on_fluid_n: [force, C64::ZERO, C64::ZERO],
        }];
        let minus = [ForceDipole {
            position_m: [position[0] - h, position[1], position[2]],
            force_on_fluid_n: [force, C64::ZERO, C64::ZERO],
        }];
        let incidents = [
            IncidentField::Sources(&quadrupoles),
            IncidentField::ForceDipoles(&plus),
            IncidentField::ForceDipoles(&minus),
        ];
        let fields = solve_scattering_batch(
            &surface, k, medium, Admittance::Rigid, &incidents, 0.01, Formulation::PlainCbie,
        ).unwrap();
        let receiver = [[-0.6, 0.25, -0.5]];
        let options = near_field::Options::default();
        let observations: Vec<_> = fields.iter().map(|f|
            f.total_pressure_velocity_at(&receiver, 0.01, options).unwrap()
        ).collect();
        let difference = |a: C64, b: C64| (a - b).scale(length / (2.0 * h));
        close(observations[0].scalar.pressure[0],
            difference(observations[1].scalar.pressure[0], observations[2].scalar.pressure[0]), 3e-6);
        for axis in 0..3 {
            close(observations[0].particle_velocity_m_s[0][axis],
                difference(observations[1].particle_velocity_m_s[0][axis],
                    observations[2].particle_velocity_m_s[0][axis]), 3e-6);
        }
        let direction = [0.6, 0.0, 0.8];
        close(fields[0].total_far_field(direction).unwrap(),
            difference(fields[1].total_far_field(direction).unwrap(),
                fields[2].total_far_field(direction).unwrap()), 3e-6);
    }

    // G0/G4: existing source caps, coherence, singularity and cancellation gates
    // apply to the new source kind without creating partial boundary data.
    #[test]
    fn mixed_sources_keep_interference_and_refuse_invalid_stress() {
        let q = source();
        let negative = StressQuadrupole {
            stress_moment_n_m: q.stress_moment_n_m.map(|v| -v),
            ..q
        };
        let mono = Monopole {
            position_m: [0.1, 0.2, 0.3],
            volume_velocity_m3_s: C64::new(1e-9, -2e-9),
        };
        let sources = [
            IncidentSource::StressQuadrupole(q),
            IncidentSource::Monopole(mono),
            IncidentSource::StressQuadrupole(negative),
        ];
        let point = [1.3, -0.8, 0.6];
        let field = IncidentField::Sources(&sources).pressure_velocity_at(2.3, Medium::air(), point).unwrap();
        let expected = IncidentField::Monopoles(&[mono]).pressure_velocity_at(2.3, Medium::air(), point).unwrap();
        close(field.pressure_pa, expected.pressure_pa, 2e-10);
        for (actual, expected) in field.velocity_m_s.into_iter().zip(expected.velocity_m_s) {
            close(actual, expected, 2e-10);
        }
        assert!(IncidentField::Sources(&sources).pressure_velocity_at(2.3, Medium::air(), [0.0; 3]).is_err());
        let mut invalid = q;
        invalid.stress_moment_n_m[4].re = f64::NAN;
        assert!(IncidentField::Sources(&[IncidentSource::StressQuadrupole(invalid)])
            .pressure_velocity_at(2.3, Medium::air(), point).is_err());
        let many = vec![IncidentSource::StressQuadrupole(q); MAX_INCIDENT_SOURCES + 1];
        assert!(IncidentField::Sources(&many).pressure_velocity_at(2.3, Medium::air(), point).is_err());
        assert!(matches!(
            IncidentField::Sources(&sources).sample_wave(2.3, 1100.0, None, point, &mut || true),
            Err(HelmholtzError::Cancelled)
        ));
        let surface = SpherePanels::icosphere(1.0, 0).unwrap();
        assert!(IncidentTraces::from_sources(&surface, 0.5, Medium::air(),
            &[IncidentSource::StressQuadrupole(q)], 0.01).is_err());
    }
}
