//! Static strain plus one phase-coherent harmonic field on the same facet.
//!
//! The spatial Bernstein hull is unchanged. At each control tensor, maximizing
//! `||mean + real*cos(t) + imaginary*sin(t)||_F` is a two-dimensional quadratic
//! problem on the unit circle. Its trust-region dual gives an upper bound at
//! every admitted multiplier; a bounded secular solve makes that bound tight
//! without phase sampling or adding separately maximized component norms.
//! See Moré and Sorensen, doi:10.1137/0904038, for the quadratic trust-region
//! optimality conditions. Arithmetic has the parent's ordinary floating-point
//! qualification; this is not an interval or continuum-error certificate.

use super::{
    HarmonicElementStrain, PlateError, PreparedElementStrain, RealElementStrain, TensorStrainBound,
    det, empty_bound, hypot, invalid, plane_stress_tensor, poll, tensor_cycle_norm, validate_layer,
};

impl PreparedElementStrain {
    /// Scalar-work allowance for all twelve biased layer control tensors,
    /// including at most 64 secular bisections per tensor. Fields are prepared
    /// separately with `REAL_FIELD_WORK` and `HARMONIC_FIELD_WORK`.
    pub const BIASED_LAYER_BOUND_WORK: usize = 32_768;
}

impl HarmonicElementStrain {
    /// Bound a complete cycle about a supplied static displacement field.
    /// Both fields must come from the SAME original prepared facet, material
    /// section and reference configuration. They add as strains in that linear
    /// trial space; this does not solve a prestressed nonlinear equilibrium or
    /// supply omitted cure/thermal eigenstrain.
    ///
    /// `component_amplitudes` now bounds the absolute TOTAL components, including
    /// the offset. A common scale applied to both input fields scales the result;
    /// scaling only the oscillation generally does not. Faces and Poisson ratio
    /// have exactly the semantics of [`Self::isotropic_layer_bound`].
    pub fn isotropic_layer_bound_with_offset(
        &self,
        offset: &RealElementStrain,
        faces_m: [f64; 2],
        poisson_ratio: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TensorStrainBound, PlateError> {
        poll(cancelled)?;
        validate_layer(faces_m, poisson_ratio)?;
        if offset.controls.iter().flatten().all(|x| *x == 0.0) {
            return self.isotropic_layer_bound(faces_m, poisson_ratio, cancelled);
        }
        if self
            .real
            .iter()
            .chain(&self.imaginary)
            .flatten()
            .all(|x| *x == 0.0)
        {
            return offset.isotropic_layer_bound(faces_m, poisson_ratio, cancelled);
        }
        let mut bound = empty_bound();
        for (face, &z) in faces_m.iter().enumerate() {
            for control in 0..6 {
                poll(cancelled)?;
                let mean = plane_stress_tensor(&offset.controls[control], z, poisson_ratio)?;
                let real = plane_stress_tensor(&self.real[control], z, poisson_ratio)?;
                let imaginary = plane_stress_tensor(&self.imaginary[control], z, poisson_ratio)?;
                let peak = tensor_biased_cycle_norm(&mean, &real, &imaginary)?;
                if peak > bound.frobenius_upper_bound {
                    bound.frobenius_upper_bound = peak;
                    bound.ruling_control_point = control;
                    bound.ruling_z_face = face;
                }
                for i in 0..6 {
                    let peak = mean[i].abs() + hypot(real[i], imaginary[i])?;
                    if !peak.is_finite() {
                        return Err(invalid("unrepresentable biased strain component envelope"));
                    }
                    bound.component_amplitudes[i] = bound.component_amplitudes[i].max(peak);
                }
            }
        }
        Ok(bound)
    }
}

fn tensor_biased_cycle_norm(
    mean: &[f64; 6],
    real: &[f64; 6],
    imaginary: &[f64; 6],
) -> Result<f64, PlateError> {
    if mean
        .iter()
        .chain(real)
        .chain(imaginary)
        .any(|x| !x.is_finite())
    {
        return Err(invalid("finite mean and harmonic strain tensors required"));
    }
    if mean.iter().all(|x| *x == 0.0) {
        return tensor_cycle_norm(real, imaginary);
    }
    let scale = mean
        .iter()
        .chain(real)
        .chain(imaginary)
        .fold(0.0_f64, |maximum, x| maximum.max(x.abs()));
    let (mut mm, mut aa, mut bb, mut ab, mut ma, mut mb) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for i in 0..6 {
        let m = mean[i] / scale;
        let a = real[i] / scale;
        let b = imaginary[i] / scale;
        let weight = if i < 3 { 1.0 } else { 2.0 };
        mm += weight * m * m;
        aa += weight * a * a;
        bb += weight * b * b;
        ab += weight * a * b;
        ma += weight * m * a;
        mb += weight * m * b;
    }
    // max_{|u|=1} mm + u^T A u + 2 g^T u. In A's eigenbasis, any
    // delta > 0 gives the dual upper bound
    // mm + lambda_max + delta + g0²/delta + g1²/(delta + gap).
    // Its minimizer satisfies (g0/delta)²+(g1/(delta+gap))² = 1.
    let diagonal = 0.5 * (aa - bb);
    let radius = hypot(diagonal, ab)?;
    let maximum = 0.5 * (aa + bb) + radius;
    let linear_norm = hypot(ma, mb)?;
    let upper = if linear_norm == 0.0 {
        mm + maximum
    } else if radius == 0.0 {
        mm + maximum + 2.0 * linear_norm
    } else {
        // Stable Jacobi eigenvector: avoid subtracting nearly equal diagonal
        // entries when the off-diagonal term is small. No angle fit is used.
        let c = det::sqrt(0.5 * (1.0 + diagonal.abs() / radius));
        let s = ab / (2.0 * radius * c);
        let vector = if diagonal >= 0.0 { [c, s] } else { [s, c] };
        let g0 = (vector[0] * ma + vector[1] * mb).abs();
        let g1 = (-vector[1] * ma + vector[0] * mb).abs();
        let gap = 2.0 * radius;
        if g0 == 0.0 && g1 <= gap {
            // The singular "hard case": complete the unit vector in the top
            // eigenspace instead of dividing by a zero secular multiplier.
            mm + maximum + g1 * (g1 / gap)
        } else {
            let mut lower = 0.0;
            let mut upper = hypot(g0, g1)?.next_up();
            for _ in 0..64 {
                let middle = 0.5 * lower + 0.5 * upper;
                if middle <= lower || middle >= upper {
                    break;
                }
                let x = g0 / middle;
                let y = g1 / (middle + gap);
                if x > 1.0 || y > 1.0 || x * x + y * y > 1.0 {
                    lower = middle;
                } else {
                    upper = middle;
                }
            }
            mm + maximum + upper + g0 * (g0 / upper) + g1 * (g1 / (upper + gap))
        }
    };
    // Preserve the dual's upper-side convergence with a small roundoff guard.
    // This does not turn ordinary f64 arithmetic into interval arithmetic.
    let peak = scale * det::sqrt(upper * (1.0 + 64.0 * f64::EPSILON));
    if !peak.is_finite() {
        return Err(invalid("unrepresentable biased harmonic tensor norm"));
    }
    Ok(peak)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 3e-12 * expected.abs().max(f64::MIN_POSITIVE),
            "actual={actual:e}, expected={expected:e}"
        );
    }

    fn norm(x: [f64; 6]) -> f64 {
        x.iter()
            .enumerate()
            .map(|(i, x)| if i < 3 { x * x } else { 2.0 * x * x })
            .sum::<f64>()
            .sqrt()
    }

    #[test]
    fn biased_cycle_matches_scalar_circular_and_hard_case_solutions() {
        let zero = [0.0; 6];
        close(
            tensor_biased_cycle_norm(
                &[0.7, 0.0, 0.0, 0.0, 0.0, 0.0],
                &[0.4, 0.0, 0.0, 0.0, 0.0, 0.0],
                &[0.3, 0.0, 0.0, 0.0, 0.0, 0.0],
            )
            .unwrap(),
            1.2,
        );
        close(
            tensor_biased_cycle_norm(
                &[0.0, 0.0, 0.7, 0.0, 0.0, 0.0],
                &[0.4, 0.0, 0.0, 0.0, 0.0, 0.0],
                &[0.0, 0.4, 0.0, 0.0, 0.0, 0.0],
            )
            .unwrap(),
            0.65_f64.sqrt(),
        );
        close(
            tensor_biased_cycle_norm(
                &[0.0, 0.5, 0.0, 0.0, 0.0, 0.0],
                &[2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                &[0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            )
            .unwrap(),
            (13.0_f64 / 3.0).sqrt(),
        );
        close(
            tensor_biased_cycle_norm(&[0.0, 0.0, 0.0, 0.2, 0.0, 0.0], &zero, &zero).unwrap(),
            2.0_f64.sqrt() * 0.2,
        );
        assert_eq!(tensor_biased_cycle_norm(&zero, &zero, &zero).unwrap(), 0.0);
    }

    #[test]
    fn biased_cycle_preserves_phase_and_extreme_representable_scales() {
        let mean = [0.7, -0.2, 0.1, 0.03, 0.04, -0.2];
        let real = [0.2, 0.5, -0.1, -0.04, 0.13, 0.07];
        let imaginary = [-0.4, 0.1, 0.2, 0.1, -0.17, 0.08];
        let expected = tensor_biased_cycle_norm(&mean, &real, &imaginary).unwrap();
        let (sin, cos) = 0.731_f64.sin_cos();
        let rotated_real = core::array::from_fn(|i| cos * real[i] - sin * imaginary[i]);
        let rotated_imaginary = core::array::from_fn(|i| sin * real[i] + cos * imaginary[i]);
        close(
            tensor_biased_cycle_norm(&mean, &rotated_real, &rotated_imaginary).unwrap(),
            expected,
        );
        for scale in [1e-160, 1e160] {
            close(
                tensor_biased_cycle_norm(
                    &mean.map(|x| x * scale),
                    &real.map(|x| x * scale),
                    &imaginary.map(|x| x * scale),
                )
                .unwrap(),
                expected * scale,
            );
        }
        let sampled = (0..20_000)
            .map(|i| {
                let (s, c) = (core::f64::consts::TAU * i as f64 / 20_000.0).sin_cos();
                norm(core::array::from_fn(|j| {
                    mean[j] + c * real[j] + s * imaginary[j]
                }))
            })
            .fold(0.0_f64, f64::max);
        assert!(expected >= sampled);
        assert!(expected - sampled < 1e-7);
    }

    #[test]
    fn static_and_harmonic_strain_share_one_layer_envelope() {
        let mean = RealElementStrain {
            controls: [[0.008, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; 6],
        };
        let harmonic = HarmonicElementStrain {
            real: [[0.003, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; 6],
            imaginary: [[0.004, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; 6],
        };
        let bound = harmonic
            .isotropic_layer_bound_with_offset(&mean, [-0.1, 0.1], 0.0, &mut || false)
            .unwrap();
        close(bound.frobenius_upper_bound, 0.013);
        close(bound.component_amplitudes[0], 0.013);
        assert!(
            mean.isotropic_layer_bound([-0.1, 0.1], 0.0, &mut || false)
                .unwrap()
                .frobenius_upper_bound
                < 0.01
        );
        assert!(
            harmonic
                .isotropic_layer_bound([-0.1, 0.1], 0.0, &mut || false)
                .unwrap()
                .frobenius_upper_bound
                < 0.01
        );
        assert!(bound.frobenius_upper_bound > 0.01);
        let zero = RealElementStrain {
            controls: [[0.0; 8]; 6],
        };
        assert_eq!(
            harmonic
                .isotropic_layer_bound_with_offset(&zero, [-0.1, 0.1], 0.0, &mut || false)
                .unwrap(),
            harmonic
                .isotropic_layer_bound([-0.1, 0.1], 0.0, &mut || false)
                .unwrap()
        );
    }

    #[test]
    fn biased_layer_bound_covers_space_thickness_and_phase() {
        let mean = RealElementStrain {
            controls: core::array::from_fn(|i| {
                core::array::from_fn(|j| ((13 * i + 5 * j) % 19) as f64 * 0.003 - 0.025)
            }),
        };
        let harmonic = HarmonicElementStrain {
            real: core::array::from_fn(|i| {
                core::array::from_fn(|j| ((7 * i + 3 * j) % 17) as f64 * 0.002 - 0.013)
            }),
            imaginary: core::array::from_fn(|i| {
                core::array::from_fn(|j| ((5 * i + 7 * j) % 23) as f64 * 0.001 - 0.011)
            }),
        };
        let faces = [-0.3, 0.4];
        let bound = harmonic
            .isotropic_layer_bound_with_offset(&mean, faces, 0.27, &mut || false)
            .unwrap();
        for i in 0..=8 {
            for j in 0..=8 - i {
                let bary = [i as f64 / 8.0, j as f64 / 8.0, (8 - i - j) as f64 / 8.0];
                let w = [
                    bary[0] * bary[0],
                    bary[1] * bary[1],
                    bary[2] * bary[2],
                    2.0 * bary[0] * bary[1],
                    2.0 * bary[1] * bary[2],
                    2.0 * bary[2] * bary[0],
                ];
                for z_index in 0..=4 {
                    let z = faces[0] + (faces[1] - faces[0]) * z_index as f64 / 4.0;
                    for phase in 0..64 {
                        let (s, c) = (core::f64::consts::TAU * phase as f64 / 64.0).sin_cos();
                        let coefficients = core::array::from_fn(|k| {
                            (0..6)
                                .map(|p| {
                                    w[p] * (mean.controls[p][k]
                                        + c * harmonic.real[p][k]
                                        + s * harmonic.imaginary[p][k])
                                })
                                .sum()
                        });
                        let tensor = plane_stress_tensor(&coefficients, z, 0.27).unwrap();
                        assert!(norm(tensor) <= bound.frobenius_upper_bound * (1.0 + 1e-12));
                        for (value, cap) in tensor.iter().zip(bound.component_amplitudes) {
                            assert!(value.abs() <= cap * (1.0 + 1e-12));
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn biased_layer_refuses_invalid_geometry_and_honors_cancellation() {
        let offset = RealElementStrain {
            controls: [[0.01; 8]; 6],
        };
        let harmonic = HarmonicElementStrain {
            real: [[0.01; 8]; 6],
            imaginary: [[0.02; 8]; 6],
        };
        assert!(
            harmonic
                .isotropic_layer_bound_with_offset(&offset, [1.0, 0.0], 0.2, &mut || false)
                .is_err()
        );
        assert!(
            harmonic
                .isotropic_layer_bound_with_offset(&offset, [0.0, 1.0], 0.5, &mut || false)
                .is_err()
        );
        let mut polls = 0;
        assert!(matches!(
            harmonic.isotropic_layer_bound_with_offset(&offset, [0.0, 1.0], 0.2, &mut || {
                polls += 1;
                polls == 4
            }),
            Err(PlateError::Cancelled)
        ));
        assert!(tensor_biased_cycle_norm(&[f64::NAN; 6], &[0.0; 6], &[0.0; 6]).is_err());
    }
}
