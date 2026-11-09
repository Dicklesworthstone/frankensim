//! Native mixed-quadrature motion: full pressure stress, compatible rates and work.
use super::*;
use fs_material::MaterialError;
use fs_material::hyper::mixed::RemainderResponse;
use fs_material::hyper::reference::{ReferenceDeformation, ReferencedHyperelastic};

fn mesh(cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    promote_p1_to_p2(
        &[
            [0.0, 0.0, 0.0],
            [1.1, 0.1, 0.0],
            [0.2, 1.0, 0.1],
            [0.0, 0.2, 0.9],
        ],
        &[[0, 1, 2, 3]],
        None,
        cx,
    )
    .unwrap()
}
fn multiply(a: [f64; 9], b: [f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}
fn transpose(a: [f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| a[3 * (ij % 3) + ij / 3])
}
fn determinant(a: [f64; 9]) -> f64 {
    a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6])
        + a[2] * (a[3] * a[7] - a[4] * a[6])
}
fn inverse(a: [f64; 9]) -> [f64; 9] {
    let j = determinant(a);
    [
        a[4] * a[8] - a[5] * a[7],
        a[2] * a[7] - a[1] * a[8],
        a[1] * a[5] - a[2] * a[4],
        a[5] * a[6] - a[3] * a[8],
        a[0] * a[8] - a[2] * a[6],
        a[2] * a[3] - a[0] * a[5],
        a[3] * a[7] - a[4] * a[6],
        a[1] * a[6] - a[0] * a[7],
        a[0] * a[4] - a[1] * a[3],
    ]
    .map(|x| x / j)
}
fn field(nodes: &[[f64; 3]], gradient: [f64; 9], shift: [f64; 3]) -> Vec<[f64; 3]> {
    nodes
        .iter()
        .map(|x| {
            core::array::from_fn(|i| {
                shift[i] + (0..3).map(|j| gradient[3 * i + j] * x[j]).sum::<f64>()
            })
        })
        .collect()
}
fn displacement(nodes: &[[f64; 3]], mut f: [f64; 9], shift: [f64; 3]) -> Vec<[f64; 3]> {
    for i in 0..3 {
        f[3 * i + i] -= 1.0;
    }
    field(nodes, f, shift)
}
fn pressure_at(x: [f64; 3]) -> f64 {
    3.0 + 0.7 * x[0] - 0.3 * x[1] + 0.4 * x[2]
}

#[test]
fn affine_mixed_motion_retains_quadrature_materials_full_pressure_and_native_force_power() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let f = [1.15, 0.08, -0.04, -0.03, 0.94, 0.05, 0.02, -0.06, 1.10];
        let fdot = [0.11, -0.07, 0.05, 0.03, -0.04, 0.09, -0.02, 0.06, 0.08];
        let u = displacement(&nodes, f, [11.0, -7.0, 3.0]);
        let v = field(&nodes, fdot, [17.0, -13.0, 5.0]);
        let laws: Vec<_> = (0..64)
            .map(|q| law(2.0 + q as f64 * 0.01, 40.0 + q as f64))
            .collect();
        let j = determinant(f);
        let inverse_transpose = transpose(inverse(f));
        let l = multiply(fdot, inverse(f));
        let c = multiply(transpose(f), f);
        let material_rate = multiply(transpose(f), fdot);
        let jdot = j * (l[0] + l[4] + l[8]);
        for treatment in [
            VolumetricTreatment::OriginalMaterial,
            VolumetricTreatment::ExactIncompressible,
        ] {
            let mut p = problem(&nodes, &elements, &laws[0], &[], &[], treatment);
            p.material = MixedLawField::PerQuadrature(&laws);
            let prepared = p.prepare(cx).unwrap();
            let pressure: Vec<_> = prepared
                .pressure_nodes()
                .iter()
                .map(|&n| pressure_at(nodes[n]))
                .collect();
            let responses = prepared
                .quadrature_motion_responses(&u, &v, &pressure, cx)
                .unwrap();
            let plain = p
                .quadrature_motion_responses(&u, &v, &pressure, cx)
                .unwrap();
            assert_eq!(responses.len(), 64);
            let i1: f64 = f.iter().map(|x| x * x).sum();
            let mut power = 0.0;
            for (q, point) in responses.iter().enumerate() {
                assert_eq!((point.element, point.local_point), (0, q));
                let pressure = pressure_at(point.reference_position_m);
                let mu = 2.0 + q as f64 * 0.01;
                // Independent mixed MR/NH stress: mu J^(-2/3)(F-I1 F^-T/3)+p cof(F).
                let piola: [f64; 9] = core::array::from_fn(|i| {
                    mu * j.powf(-2.0 / 3.0) * (f[i] - i1 / 3.0 * inverse_transpose[i])
                        + pressure * j * inverse_transpose[i]
                });
                let sigma = multiply(piola, transpose(f)).map(|x| x / j);
                for ij in 0..9 {
                    let ji = 3 * (ij % 3) + ij / 3;
                    near(point.deformation_gradient[ij], f[ij], 2e-13);
                    near(point.deformation_gradient_rate_per_s[ij], fdot[ij], 2e-13);
                    near(point.material_state.first_piola_pa[ij], piola[ij], 2e-11);
                    near(point.cauchy_pa[ij], sigma[ij], 2e-11);
                    near(
                        point.green_lagrange_strain[ij],
                        0.5 * (c[ij] - if ij / 3 == ij % 3 { 1.0 } else { 0.0 }),
                        3e-13,
                    );
                    near(point.spatial_velocity_gradient_per_s[ij], l[ij], 3e-13);
                    near(
                        point.spatial_strain_rate_per_s[ij],
                        0.5 * (l[ij] + l[ji]),
                        3e-13,
                    );
                    near(
                        point.green_lagrange_strain_rate_per_s[ij],
                        0.5 * (material_rate[ij] + material_rate[ji]),
                        3e-13,
                    );
                }
                near(point.pressure_pa, pressure, 3e-13);
                near(point.jacobian, j, 3e-13);
                near(point.jacobian_rate_per_s, jdot, 3e-13);
                let bulk_volume = if treatment == VolumetricTreatment::ExactIncompressible {
                    1.0
                } else {
                    1.0 + pressure / (40.0 + q as f64)
                };
                near(
                    point.material_state.pressure_residual,
                    j - bulk_volume,
                    4e-13,
                );
                let work = point.reference_volume_m3
                    * piola.iter().zip(fdot).map(|(a, b)| a * b).sum::<f64>();
                near(point.reference_stress_power_w, work, 2e-12);
                near(point.spatial_stress_power_w, work, 2e-12);
                assert_eq!(point.deformation_gradient, plain[q].deformation_gradient);
                assert_eq!(
                    point.reference_stress_power_w,
                    plain[q].reference_stress_power_w
                );
                power += point.reference_stress_power_w;
            }
            let assembled = prepared.evaluate(&u, &pressure, false, cx).unwrap();
            let nodal_power: f64 = assembled
                .internal_forces_n
                .iter()
                .zip(v.iter().flatten())
                .map(|(a, b)| a * b)
                .sum();
            near(power, nodal_power, 3e-10);
        }
    });
}

#[test]
fn quadratic_velocity_uses_each_original_point_and_preserves_weak_volume_rates() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(2.0, 80.0);
        let p = problem(
            &nodes,
            &elements,
            &material,
            &[],
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        let prepared = p.prepare(cx).unwrap();
        let u: Vec<_> = nodes
            .iter()
            .map(|[x, y, z]| [0.04 * x * x + 0.03 * y * z, -0.02 * x * y, 0.05 * z * z])
            .collect();
        let v: Vec<_> = nodes
            .iter()
            .map(|[x, y, z]| [0.07 * x * x, -0.06 * y * z, 0.04 * x * y])
            .collect();
        let pressure = vec![2.0; prepared.pressure_nodes().len()];
        let responses = prepared
            .quadrature_motion_responses(&u, &v, &pressure, cx)
            .unwrap();
        for point in &responses {
            let [x, y, z] = point.reference_position_m;
            let f = [
                1.0 + 0.08 * x,
                0.03 * z,
                0.03 * y,
                -0.02 * y,
                1.0 - 0.02 * x,
                0.0,
                0.0,
                0.0,
                1.0 + 0.10 * z,
            ];
            let fdot = [
                0.14 * x,
                0.0,
                0.0,
                0.0,
                -0.06 * z,
                -0.06 * y,
                0.04 * y,
                0.04 * x,
                0.0,
            ];
            for ij in 0..9 {
                near(point.deformation_gradient[ij], f[ij], 3e-14);
                near(point.deformation_gradient_rate_per_s[ij], fdot[ij], 3e-14);
            }
            near(
                point.material_state.pressure_residual,
                determinant(f) - 1.0,
                3e-14,
            );
        }
        assert!(responses.iter().any(|r| r.jacobian_rate_per_s.abs() > 1e-3));
        assert!(
            responses
                .iter()
                .any(|r| r.material_state.pressure_residual.abs() > 1e-3)
        );
        let h = 2e-6;
        let endpoint = |sign: f64| {
            let u: Vec<_> = u
                .iter()
                .zip(&v)
                .map(|(u, v)| core::array::from_fn(|i| u[i] + sign * h * v[i]))
                .collect();
            prepared
                .evaluate(&u, &pressure, false, cx)
                .unwrap()
                .pressure_residual_m3
                .iter()
                .sum::<f64>()
        };
        let integrated_rate: f64 = responses.iter().map(|r| r.current_volume_rate_m3_s).sum();
        near(
            (endpoint(1.0) - endpoint(-1.0)) / (2.0 * h),
            integrated_rate,
            2e-10,
        );
    });
}

#[test]
fn mixed_motion_is_objective_under_a_time_dependent_three_dimensional_frame() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(2.0, 80.0);
        let p = problem(
            &nodes,
            &elements,
            &material,
            &[],
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        let prepared = p.prepare(cx).unwrap();
        let pressure = vec![3.0; prepared.pressure_nodes().len()];
        let f = [1.15, 0.08, -0.04, -0.03, 0.94, 0.05, 0.02, -0.06, 1.10];
        let fdot = [0.11, -0.07, 0.05, 0.03, -0.04, 0.09, -0.02, 0.06, 0.08];
        let q = multiply(
            [0.6, -0.8, 0.0, 0.8, 0.6, 0.0, 0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 0.0, 0.8, -0.6, 0.0, 0.6, 0.8],
        );
        let spin = [0.0, -0.7, 0.2, 0.7, 0.0, -0.4, -0.2, 0.4, 0.0];
        let qf = multiply(q, f);
        let qfdot = multiply(q, fdot);
        let spin_f = multiply(spin, qf);
        let moving_rate = core::array::from_fn(|i| qfdot[i] + spin_f[i]);
        let observe = |f, rate, shift, velocity| {
            prepared
                .quadrature_motion_responses(
                    &displacement(&nodes, f, shift),
                    &field(&nodes, rate, velocity),
                    &pressure,
                    cx,
                )
                .unwrap()
        };
        let before = observe(f, fdot, [0.0; 3], [0.0; 3]);
        let after = observe(qf, moving_rate, [4.0, -7.0, 3.0], [13.0, -17.0, 5.0]);
        for (a, b) in before.iter().zip(after) {
            let rotated_d = multiply(multiply(q, a.spatial_strain_rate_per_s), transpose(q));
            let rotated_p = multiply(q, a.material_state.first_piola_pa);
            for ij in 0..9 {
                near(
                    b.green_lagrange_strain_rate_per_s[ij],
                    a.green_lagrange_strain_rate_per_s[ij],
                    4e-12,
                );
                near(b.spatial_strain_rate_per_s[ij], rotated_d[ij], 4e-12);
                near(b.material_state.first_piola_pa[ij], rotated_p[ij], 4e-11);
            }
            near(b.jacobian_rate_per_s, a.jacobian_rate_per_s, 4e-12);
            near(
                b.reference_stress_power_w,
                a.reference_stress_power_w,
                4e-12,
            );
        }
        for rigid in observe(qf, spin_f, [4.0, -7.0, 3.0], [13.0, -17.0, 5.0]) {
            for value in rigid
                .spatial_strain_rate_per_s
                .iter()
                .chain(&rigid.green_lagrange_strain_rate_per_s)
            {
                near(*value, 0.0, 4e-12);
            }
            near(rigid.reference_stress_power_w, 0.0, 4e-12);
        }
    });
}

#[test]
fn mixed_motion_keeps_the_original_natural_reference_and_pressure_convention() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let original = law(2.0, 80.0);
        let g = [1.12, 0.08, -0.03, 0.0, 0.97, 0.04, 0.0, 0.0, 1.08];
        let material =
            ReferencedHyperelastic::new(&original, ReferenceDeformation::try_new(g).unwrap());
        let p = problem(
            &nodes,
            &elements,
            &material,
            &[],
            &[],
            VolumetricTreatment::ExactIncompressible,
        );
        let prepared = p.prepare(cx).unwrap();
        let pressure = vec![1.3; prepared.pressure_nodes().len()];
        let u = displacement(&nodes, g, [0.0; 3]);
        let v = field(
            &nodes,
            [0.04, -0.03, 0.02, 0.01, -0.05, 0.07, -0.02, 0.06, 0.03],
            [0.0; 3],
        );
        for point in prepared
            .quadrature_motion_responses(&u, &v, &pressure, cx)
            .unwrap()
        {
            near(point.jacobian, determinant(g), 3e-14);
            near(
                point.material_state.bulk_volume_ratio,
                determinant(g),
                3e-14,
            );
            near(point.material_state.pressure_residual, 0.0, 3e-14);
            for ij in 0..9 {
                near(
                    point.cauchy_pa[ij],
                    if ij / 3 == ij % 3 { 1.3 } else { 0.0 },
                    3e-13,
                );
            }
            near(
                point.reference_stress_power_w,
                1.3 * point.current_volume_rate_m3_s,
                3e-13,
            );
        }
    });
}

struct OriginalOwner<'a> {
    law: &'a Hyperelastic,
    gate: Option<&'a CancelGate>,
}
impl MixedFiniteStrainPotential for OriginalOwner<'_> {
    fn validate_mixed(&self) -> Result<(), MaterialError> {
        self.law.validate_mixed()
    }
    fn mixed_pressure_potential(
        &self,
        f: &[f64; 9],
        pressure_pa: f64,
        treatment: VolumetricTreatment,
    ) -> Result<fs_material::hyper::mixed::PressurePotential, MaterialError> {
        self.law.mixed_pressure_potential(f, pressure_pa, treatment)
    }
    fn remainder_stress_and_potential(
        &self,
        f: &[f64; 9],
    ) -> Result<([f64; 9], f64), MaterialError> {
        if let Some(gate) = self.gate {
            gate.request();
        }
        let (mut p, w) = self.law.remainder_stress_and_potential(f)?;
        p[1] += 2.0;
        Ok((p, w + 2.0 * f[1]))
    }
    fn remainder_response(&self, _: &[f64; 9]) -> Result<RemainderResponse, MaterialError> {
        panic!("motion observation must not request a constitutive tangent")
    }
}

#[test]
fn mixed_observation_retains_full_source_stress_tiny_rates_and_refuses_partial_results() {
    let gate = CancelGate::new();
    with_gate(&gate, |cx| {
        let (nodes, elements) = mesh(cx);
        let original = law(2.0, 80.0);
        let material = OriginalOwner {
            law: &original,
            gate: None,
        };
        let p = problem(
            &nodes,
            &elements,
            &material,
            &[],
            &[],
            VolumetricTreatment::OriginalMaterial,
        );
        let prepared = p.prepare(cx).unwrap();
        let u = vec![[0.0; 3]; nodes.len()];
        let pressure = vec![0.0; prepared.pressure_nodes().len()];
        let tiny = 1e-30;
        let v = field(
            &nodes,
            [0.0, tiny, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0; 3],
        );
        for point in prepared
            .quadrature_motion_responses(&u, &v, &pressure, cx)
            .unwrap()
        {
            assert!((point.deformation_gradient_rate_per_s[1] / tiny - 1.0).abs() < 2e-13);
            near(point.cauchy_pa[1], 2.0, 2e-14);
            near(point.cauchy_pa[3], 0.0, 2e-14);
            assert!(
                (point.reference_stress_power_w / (point.reference_volume_m3 * tiny) - 2.0).abs()
                    < 4e-13
            );
        }
        assert!(
            prepared
                .quadrature_motion_responses(&u, &v[..v.len() - 1], &pressure, cx)
                .is_err()
        );
        assert!(
            prepared
                .quadrature_motion_responses(&u, &vec![[f64::NAN; 3]; nodes.len()], &pressure, cx)
                .is_err()
        );
        assert!(
            prepared
                .quadrature_motion_responses(&u, &v, &[f64::INFINITY; 4], cx)
                .is_err()
        );
        let inverted = field(
            &nodes,
            [-2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0; 3],
        );
        assert!(matches!(
            prepared.quadrature_motion_responses(&inverted, &v, &pressure, cx),
            Err(MixedUpError::Material { .. })
        ));
        for point in prepared
            .quadrature_motion_responses(&u, &vec![[3.0, -7.0, 5.0]; nodes.len()], &pressure, cx)
            .unwrap()
        {
            assert_eq!(point.deformation_gradient_rate_per_s, [0.0; 9]);
            assert_eq!(point.reference_stress_power_w, 0.0);
        }
        let cancelling = OriginalOwner {
            law: &original,
            gate: Some(&gate),
        };
        let p = problem(
            &nodes,
            &elements,
            &cancelling,
            &[],
            &[],
            VolumetricTreatment::OriginalMaterial,
        );
        let prepared = p.prepare(cx).unwrap();
        assert!(matches!(
            prepared.quadrature_motion_responses(&u, &v, &pressure, cx),
            Err(MixedUpError::Cancelled)
        ));
    });
}
