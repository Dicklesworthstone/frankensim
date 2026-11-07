use fs_material::laminate::{
    LaminaLayer, LaminateSection, LaminateShearRecovery, LayerLaw, ReferenceSurface,
};
use fs_math::c64::C64;
use fs_plate::PlateError;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::shear_recovery::LaminateShearField;
use fs_plate::shell::mitc::strain::integration_rule;
use fs_plate::shell::mitc::thermal::{FacetThermalStrain, ShellThermalStrain};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellSection,
    ShellConstraint, ShellMass, assemble_mitc3p, solve_static,
};
use std::sync::Arc;

fn near(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() < 2e-10 * scale.max(a.abs()).max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn source() -> LaminateShearRecovery {
    let layers = [(0.012, 0.3, 1.0), (0.008, -0.5, 0.2)].map(|(h, angle, scale)| LaminaLayer {
        thickness_m: h,
        orientation_rad: angle,
        density_kg_m3: 1100.0,
        material: None,
        law: LayerLaw::PlaneStress {
            stiffness_pa: [
                2e8 * scale,
                0.1e8 * scale,
                0.0,
                0.1e8 * scale,
                0.8e8 * scale,
                0.0,
                0.0,
                0.0,
                0.3e8 * scale,
            ],
            transverse_shear_pa: [0.4e8 * scale, 0.1e8 * scale],
        },
    });
    LaminateShearRecovery::try_new(Arc::new(
        LaminateSection::try_new_equilibrated_shear(
            layers.to_vec(),
            ReferenceSurface::Midsurface,
            0.12,
        )
        .unwrap(),
    ))
    .unwrap()
}
fn model(source: &LaminateShearRecovery, angle: f64, fixed: bool) -> (ShellMesh, MitcShellModel) {
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.2, 0.6, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.9).unwrap();
    let section = ShearShellSection::from_laminate(source.section(), angle).unwrap();
    let constraints = if fixed {
        vec![
            ShellConstraint {
                node: 0,
                mask: [true; 6],
                frame: ConstraintFrame::Global,
            },
            ShellConstraint {
                node: 1,
                mask: [true; 6],
                frame: ConstraintFrame::Global,
            },
        ]
    } else {
        vec![]
    };
    let model = assemble_mitc3p(
        &mesh,
        &directors,
        &[section],
        &constraints,
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: 3,
            max_triangles: 1,
        },
        1e-3,
        || false,
    )
    .unwrap();
    (mesh, model)
}

// Independently integrate the recovered local material work and traction,
// then compare with the original native assumed-shear field and rotated H.
fn integrate(
    field: &LaminateShearField<'_>,
    source: &LaminateShearRecovery,
    angle: f64,
    area: f64,
) -> ([f64; 2], f64) {
    let (s, c) = angle.sin_cos();
    let h = source.section().shear_pa_m();
    let g = (3.0_f64 / 5.0).sqrt();
    let mut average_q = [0.0; 2];
    let mut energy = 0.0;
    for (lambda, weight) in integration_rule() {
        let gamma = field
            .native_strain()
            .at_barycentric(lambda)
            .unwrap()
            .transverse_shear;
        let e = [c * gamma[0] + s * gamma[1], -s * gamma[0] + c * gamma[1]];
        let q = [h[0] * e[0] + h[1] * e[1], h[2] * e[0] + h[3] * e[1]];
        let expected = [c * q[0] - s * q[1], s * q[0] + c * q[1]];
        let mut local_q = [0.0; 2];
        let mut local_energy = 0.0;
        for (layer, input) in source.section().layers().iter().enumerate() {
            for (fraction, w) in [
                (0.5 * (1.0 - g), 5.0 / 18.0),
                (0.5, 4.0 / 9.0),
                (0.5 * (1.0 + g), 5.0 / 18.0),
            ] {
                let response = field.at_barycentric(lambda, layer, fraction).unwrap();
                for k in 0..2 {
                    local_q[k] += w * input.thickness_m * response.facet_stress_pa[k];
                }
                local_energy += w * input.thickness_m * response.layer.energy_density_j_m3;
                near(
                    response.layer.energy_density_j_m3,
                    0.5 * (response.facet_stress_pa[0] * response.facet_engineering_strain[0]
                        + response.facet_stress_pa[1] * response.facet_engineering_strain[1]),
                    1e-6,
                );
            }
        }
        for k in 0..2 {
            near(local_q[k], expected[k], 1e-4);
            average_q[k] += weight * local_q[k];
        }
        near(local_energy, 0.5 * (e[0] * q[0] + e[1] * q[1]), 1e-8);
        energy += area * weight * local_energy;
    }
    (average_q, energy)
}

fn membrane_bending_work(
    field: &LaminateShearField<'_>,
    source: &LaminateShearRecovery,
    angle: f64,
    area: f64,
) -> f64 {
    let (s, c) = angle.sin_cos();
    let rotate = |v: [f64; 3]| {
        [
            c * c * v[0] + s * s * v[1] + s * c * v[2],
            s * s * v[0] + c * c * v[1] - s * c * v[2],
            -2.0 * s * c * v[0] + 2.0 * s * c * v[1] + (c * c - s * s) * v[2],
        ]
    };
    let product = |a: [f64; 3], matrix: &[f64; 9], b: [f64; 3]| -> f64 {
        (0..3)
            .map(|i| (0..3).map(|j| a[i] * matrix[3 * i + j] * b[j]).sum::<f64>())
            .sum()
    };
    integration_rule()
        .into_iter()
        .map(|(lambda, weight)| {
            let point = field.native_strain().at_barycentric(lambda).unwrap();
            let membrane = rotate(point.membrane);
            let bending = rotate(point.gradient_per_m);
            area * weight
                * (0.5 * product(membrane, source.section().a_pa_m(), membrane)
                    + product(membrane, source.section().b_pa_m2(), bending)
                    + 0.5 * product(bending, source.section().d_pa_m3(), bending))
        })
        .sum()
}

#[test]
fn mixed_shear_ply_work_matches_actual_native_shell_operator() {
    let source = source();
    let angle = 0.43;
    let (mesh, model) = model(&source, angle, false);
    let prepared = model
        .prepare_laminate_shear(0, &source, angle, || false)
        .unwrap();
    for gamma in [[0.002, 0.003], [0.002, -0.003]] {
        let u: Vec<_> = mesh
            .nodes
            .iter()
            .flat_map(|p| [0.0, 0.0, gamma[0] * p[0] + gamma[1] * p[1], 0.0, 0.0, 0.0])
            .collect();
        let field = prepared.evaluate(&u, || false).unwrap();
        let area = mesh.facet(0).unwrap().area_m2;
        let (_, energy) = integrate(&field, &source, angle, area);
        // Condensing the interior directors also creates bending, even though
        // the physical nodal rotations are zero. The recovered ply shear is
        // only H's contribution, not the total A/B/D/H operator's energy.
        let total = energy + membrane_bending_work(&field, &source, angle, area);
        near(total, model.strain_energy(&u).unwrap(), total);
        for lambda in [[0.2, 0.3, 0.5], [1.0 / 3.0; 3]] {
            let lower = field.at_barycentric(lambda, 0, 1.0).unwrap();
            let upper = field.at_barycentric(lambda, 1, 0.0).unwrap();
            for k in 0..2 {
                near(lower.facet_stress_pa[k], upper.facet_stress_pa[k], 1e5);
            }
        }
    }
}

#[test]
fn solved_nonuniform_field_retains_native_shear_and_condensed_directors() {
    let source = source();
    let angle = -0.27;
    let (mesh, model) = model(&source, angle, true);
    let mut load = [0.0; 18];
    load[14] = 0.12;
    load[15] = -0.02;
    let u = solve_static(&model, &load).unwrap();
    let prepared = model
        .prepare_laminate_shear(0, &source, angle, || false)
        .unwrap();
    let field = prepared.evaluate(&u, || false).unwrap();
    let (q, energy) = integrate(&field, &source, angle, mesh.facet(0).unwrap().area_m2);
    assert!(energy > 0.0);
    let zero_thermal = model
        .prepare_thermal(
            &[ShellThermalStrain::try_new([0.0; 3], [0.0; 3])
                .unwrap()
                .into()],
            || false,
        )
        .unwrap();
    let native = zero_thermal.recover(&u, || false).unwrap();
    for k in 0..2 {
        near(q[k], native.facets[0].transverse_shear_n_m[k], 1e-3);
    }
    assert!(energy <= native.strain_energy_j * (1.0 + 1e-10));
}

#[test]
fn thermal_profile_recovers_nonzero_internal_shear_even_at_fixed_nodes() {
    let source = source();
    let angle = 0.0;
    let (mesh, model) = model(&source, angle, false);
    let loading = FacetThermalStrain::AtNodes(core::array::from_fn(|node| {
        ShellThermalStrain::try_new([0.0; 3], [0.02 * mesh.nodes[node][0]; 3]).unwrap()
    }));
    let thermal = model.prepare_thermal(&[loading], || false).unwrap();
    let field = thermal
        .laminate_shear_field(0, &source, angle, &[0.0; 18], || false)
        .unwrap();
    let (_, energy) = integrate(&field, &source, angle, mesh.facet(0).unwrap().area_m2);
    assert!(energy > 1e-14);
    let bare = model
        .prepare_laminate_shear(0, &source, angle, || false)
        .unwrap();
    let cold = bare.evaluate(&[0.0; 18], || false).unwrap();
    assert_eq!(
        cold.at_barycentric([0.2, 0.3, 0.5], 0, 0.5)
            .unwrap()
            .facet_stress_pa,
        [0.0; 2]
    );
    let expected = thermal.strain_field(0, &[0.0; 18], || false).unwrap();
    for (lambda, _) in integration_rule() {
        assert_eq!(
            field.native_strain().at_barycentric(lambda).unwrap(),
            expected.at_barycentric(lambda).unwrap()
        );
    }
}

#[test]
fn source_orientation_shape_and_cancellation_refuse_without_operator_mutation() {
    let source = source();
    let (_, model) = model(&source, 0.17, false);
    assert!(
        model
            .prepare_laminate_shear(0, &source, 0.18, || false)
            .is_err()
    );
    assert!(
        model
            .prepare_laminate_shear(1, &source, 0.17, || false)
            .is_err()
    );
    assert!(matches!(
        model.prepare_laminate_shear(0, &source, 0.17, || true),
        Err(PlateError::Cancelled)
    ));
    let prepared = model
        .prepare_laminate_shear(0, &source, 0.17, || false)
        .unwrap();
    assert!(prepared.evaluate(&[0.0; 17], || false).is_err());
    assert!(prepared.evaluate(&[f64::NAN; 18], || false).is_err());
    assert!(matches!(
        prepared.evaluate(&[0.0; 18], || true),
        Err(PlateError::Cancelled)
    ));
    let field = prepared.evaluate(&[0.0; 18], || false).unwrap();
    assert!(field.at_barycentric([-0.1, 0.6, 0.5], 0, 0.5).is_err());
    assert!(field.at_barycentric([0.2, 0.3, 0.5], 2, 0.5).is_err());
    assert_eq!(model.strain_energy(&[0.0; 18]).unwrap(), 0.0);
}

fn isotropic_source(moduli: &[f64]) -> LaminateShearRecovery {
    let nu = 0.25;
    let layers = moduli
        .iter()
        .map(|&e| {
            let q = e / (1.0 - nu * nu);
            LaminaLayer {
                thickness_m: 0.02 / moduli.len() as f64,
                orientation_rad: 0.0,
                density_kg_m3: 1100.0,
                material: None,
                law: LayerLaw::PlaneStress {
                    stiffness_pa: [
                        q,
                        nu * q,
                        0.0,
                        nu * q,
                        q,
                        0.0,
                        0.0,
                        0.0,
                        0.5 * (1.0 - nu) * q,
                    ],
                    transverse_shear_pa: [e / (2.0 * (1.0 + nu)); 2],
                },
            }
        })
        .collect();
    LaminateShearRecovery::try_new(Arc::new(
        LaminateSection::try_new_equilibrated_shear(layers, ReferenceSurface::Midsurface, 0.12)
            .unwrap(),
    ))
    .unwrap()
}

fn harmonic_motion(mesh: &ShellMesh) -> Vec<C64> {
    mesh.nodes
        .iter()
        .flat_map(|p| {
            [
                C64::new(0.001 * p[0], 0.0004 * p[1]),
                C64::new(-0.0002 * p[1], -0.0006 * p[0]),
                C64::new(0.002 * p[0] + 0.001 * p[1], 0.001 * p[0] - 0.0003 * p[1]),
                C64::new(0.001 + 0.0003 * p[1], -0.0005 * p[1]),
                C64::new(-0.002 + 0.0007 * p[0], 0.0002 + 0.0004 * p[0]),
                C64::new(0.0, 0.0),
            ]
        })
        .collect()
}

fn instantaneous(value: C64, phase: f64) -> f64 {
    value.re * phase.cos() - value.im * phase.sin()
}

#[test]
fn cached_bernstein_maps_reconstruct_original_ply_profiles() {
    // G0: the new coefficient representation reproduces the original
    // quadratic traction/compliance owner, including off-axis cross terms.
    let source = source();
    for layer in 0..source.section().layers().len() {
        let maps = source.strain_bernstein_controls(layer).unwrap();
        for gamma in [[0.002, -0.003], [-0.001, 0.004], [0.0, 0.0]] {
            for fraction in [0.0_f64, 0.17, 0.5, 0.83, 1.0] {
                let w = [
                    (1.0 - fraction).powi(2),
                    2.0 * fraction * (1.0 - fraction),
                    fraction * fraction,
                ];
                let expected = source.from_section_strain(layer, fraction, gamma).unwrap();
                for component in 0..2 {
                    let actual = (0..3)
                        .map(|k| {
                            w[k] * (maps[k][2 * component] * gamma[0]
                                + maps[k][2 * component + 1] * gamma[1])
                        })
                        .sum::<f64>();
                    near(actual, expected.engineering_strain[component], 1e-8);
                }
            }
        }
    }
    assert!(
        source
            .strain_bernstein_controls(source.section().layers().len())
            .is_err()
    );
}

#[test]
fn coherent_harmonic_ply_response_matches_real_motion_and_mean_work() {
    // G3: every common phase reconstructs the original real snapshot; mean
    // energy is checked by independent temporal quadrature of those snapshots.
    let source = source();
    let (mesh, model) = model(&source, 0.43, false);
    let prepared = model
        .prepare_laminate_shear(0, &source, 0.43, || false)
        .unwrap();
    let u = harmonic_motion(&mesh);
    let harmonic = prepared.evaluate_harmonic(&u, || false).unwrap();
    let lambda = [0.2, 0.3, 0.5];
    for layer in 0..source.section().layers().len() {
        let complex = harmonic.at_barycentric(lambda, layer, 0.37).unwrap();
        let mut mean_energy = 0.0;
        for step in 0..32 {
            let phase = core::f64::consts::TAU * step as f64 / 32.0;
            let snapshot: Vec<_> = u.iter().map(|&value| instantaneous(value, phase)).collect();
            let real = prepared.evaluate(&snapshot, || false).unwrap();
            let point = real.at_barycentric(lambda, layer, 0.37).unwrap();
            for component in 0..2 {
                near(
                    instantaneous(complex.facet_stress_pa[component], phase),
                    point.facet_stress_pa[component],
                    1e-5,
                );
                near(
                    instantaneous(complex.facet_engineering_strain[component], phase),
                    point.facet_engineering_strain[component],
                    1e-8,
                );
                near(
                    instantaneous(complex.material_engineering_strain[component], phase),
                    point.layer.material_engineering_strain[component],
                    1e-8,
                );
            }
            mean_energy += point.layer.energy_density_j_m3 / 32.0;
        }
        near(complex.average_energy_density_j_m3, mean_energy, 1e-6);
        let lower = harmonic.at_barycentric(lambda, 0, 1.0).unwrap();
        let upper = harmonic.at_barycentric(lambda, 1, 0.0).unwrap();
        for component in 0..2 {
            near(
                lower.facet_stress_pa[component].re,
                upper.facet_stress_pa[component].re,
                1e-5,
            );
            near(
                lower.facet_stress_pa[component].im,
                upper.facet_stress_pa[component].im,
                1e-5,
            );
        }
    }
}

#[test]
fn product_bounds_cover_triangle_thickness_and_common_harmonic_phase() {
    // G0/G3: exact coefficient conversion plus dense independent point/phase
    // evaluation. The implementation uses convexity, not this sampling grid,
    // to cover the complete represented triangle/thickness/cycle domain.
    let source = isotropic_source(&[2e8, 2e6]);
    let (mesh, model) = model(&source, -0.24, false);
    let prepared = model
        .prepare_laminate_shear(0, &source, -0.24, || false)
        .unwrap();
    let u = harmonic_motion(&mesh);
    let harmonic = prepared.evaluate_harmonic(&u, || false).unwrap();
    let real_u: Vec<_> = u.iter().map(|value| value.re).collect();
    let real = prepared.evaluate(&real_u, || false).unwrap();
    let shift = C64::new(0.61_f64.cos(), 0.61_f64.sin());
    let shifted_u: Vec<_> = u.iter().map(|&value| value * shift).collect();
    let shifted = prepared.evaluate_harmonic(&shifted_u, || false).unwrap();
    for layer in 0..2 {
        let bound = harmonic
            .isotropic_layer_bound(layer, 0.25, &mut || false)
            .unwrap();
        let real_bound = real
            .isotropic_layer_bound(layer, 0.25, &mut || false)
            .unwrap();
        let shifted_bound = shifted
            .isotropic_layer_bound(layer, 0.25, &mut || false)
            .unwrap();
        near(
            bound.frobenius_upper_bound,
            shifted_bound.frobenius_upper_bound,
            1e-8,
        );
        for component in 0..6 {
            near(
                bound.component_amplitudes[component],
                shifted_bound.component_amplitudes[component],
                1e-8,
            );
        }
        for i in 0..=8 {
            for j in 0..=8 - i {
                let lambda = [i as f64 / 8.0, j as f64 / 8.0, (8 - i - j) as f64 / 8.0];
                let native = harmonic.native_strain().at_barycentric(lambda).unwrap();
                for depth in 0..=12 {
                    let response = harmonic
                        .at_barycentric(lambda, layer, depth as f64 / 12.0)
                        .unwrap();
                    let in_plane: [C64; 3] = core::array::from_fn(|k| {
                        C64::new(
                            native.real.membrane[k] + response.z_m * native.real.gradient_per_m[k],
                            native.imaginary.membrane[k]
                                + response.z_m * native.imaginary.gradient_per_m[k],
                        )
                    });
                    for step in 0..24 {
                        let phase = core::f64::consts::TAU * step as f64 / 24.0;
                        let xx = instantaneous(in_plane[0], phase);
                        let yy = instantaneous(in_plane[1], phase);
                        let tensor = [
                            xx,
                            yy,
                            -(xx + yy) / 3.0,
                            0.5 * instantaneous(in_plane[2], phase),
                            0.5 * instantaneous(response.facet_engineering_strain[1], phase),
                            0.5 * instantaneous(response.facet_engineering_strain[0], phase),
                        ];
                        let norm = tensor
                            .iter()
                            .enumerate()
                            .map(|(k, x)| if k < 3 { x * x } else { 2.0 * x * x })
                            .sum::<f64>()
                            .sqrt();
                        assert!(norm <= bound.frobenius_upper_bound * (1.0 + 2e-10));
                        if step == 0 {
                            assert!(norm <= real_bound.frobenius_upper_bound * (1.0 + 2e-10));
                        }
                        for k in 0..6 {
                            assert!(
                                tensor[k].abs()
                                    <= bound.component_amplitudes[k]
                                        + 2e-10 * bound.frobenius_upper_bound
                            );
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn physical_shear_inside_a_ply_is_not_bounded_by_its_free_faces() {
    let source = isotropic_source(&[2e8]);
    let (mesh, model) = model(&source, 0.0, false);
    let prepared = model
        .prepare_laminate_shear(0, &source, 0.0, || false)
        .unwrap();
    let u = harmonic_motion(&mesh);
    let field = prepared.evaluate_harmonic(&u, || false).unwrap();
    let lambda = [0.2, 0.3, 0.5];
    let bottom = field.at_barycentric(lambda, 0, 0.0).unwrap();
    let middle = field.at_barycentric(lambda, 0, 0.5).unwrap();
    let top = field.at_barycentric(lambda, 0, 1.0).unwrap();
    let bound = field.isotropic_layer_bound(0, 0.25, &mut || false).unwrap();
    let mut interior = 0.0_f64;
    for k in 0..2 {
        let amplitude = middle.facet_engineering_strain[k].abs();
        interior = interior.max(amplitude);
        assert!(bottom.facet_engineering_strain[k].abs() < 1e-12 * amplitude.max(1e-12));
        assert!(top.facet_engineering_strain[k].abs() < 1e-12 * amplitude.max(1e-12));
        assert!(0.5 * amplitude <= bound.component_amplitudes[5 - k] * (1.0 + 2e-10));
    }
    assert!(interior > 1e-8);
}

#[test]
fn harmonic_and_whole_ply_refusals_leave_the_observation_reusable() {
    let source = isotropic_source(&[2e8]);
    let (mesh, model) = model(&source, 0.0, false);
    let prepared = model
        .prepare_laminate_shear(0, &source, 0.0, || false)
        .unwrap();
    let u = harmonic_motion(&mesh);
    assert!(prepared.evaluate_harmonic(&u[..17], || false).is_err());
    let mut invalid_u = u.clone();
    invalid_u[2].im = f64::NAN;
    assert!(prepared.evaluate_harmonic(&invalid_u, || false).is_err());
    assert!(matches!(
        prepared.evaluate_harmonic(&u, || true),
        Err(PlateError::Cancelled)
    ));
    let field = prepared.evaluate_harmonic(&u, || false).unwrap();
    let before = field.isotropic_layer_bound(0, 0.25, &mut || false).unwrap();
    assert!(field.at_barycentric([-0.1, 0.5, 0.6], 0, 0.5).is_err());
    assert!(field.at_barycentric([0.2, 0.3, 0.5], 0, 1.01).is_err());
    assert!(field.isotropic_layer_bound(1, 0.25, &mut || false).is_err());
    for nu in [f64::NAN, -1.0, 0.5] {
        assert!(field.isotropic_layer_bound(0, nu, &mut || false).is_err());
    }
    let mut polls = 0;
    assert!(matches!(
        field.isotropic_layer_bound(0, 0.25, &mut || {
            polls += 1;
            polls == 4
        }),
        Err(PlateError::Cancelled)
    ));
    assert_eq!(
        before,
        field.isotropic_layer_bound(0, 0.25, &mut || false).unwrap()
    );
    let zero = prepared
        .evaluate_harmonic(&[C64::ZERO; 18], || false)
        .unwrap();
    assert_eq!(
        zero.isotropic_layer_bound(0, 0.25, &mut || false)
            .unwrap()
            .frobenius_upper_bound,
        0.0
    );
}

#[test]
fn free_thermal_expansion_does_not_receive_a_mechanical_normal_strain_closure() {
    let source = isotropic_source(&[2e8]);
    let (mesh, model) = model(&source, 0.0, false);
    let free_strain = 1e-4;
    let thermal = model
        .prepare_thermal(
            &[
                ShellThermalStrain::try_new([free_strain, free_strain, 0.0], [0.0; 3])
                    .unwrap()
                    .into(),
            ],
            || false,
        )
        .unwrap();
    let u: Vec<_> = mesh
        .nodes
        .iter()
        .flat_map(|p| [free_strain * p[0], free_strain * p[1], 0.0, 0.0, 0.0, 0.0])
        .collect();
    let field = thermal
        .laminate_shear_field(0, &source, 0.0, &u, || false)
        .unwrap();
    let lambda = [0.2, 0.3, 0.5];
    let native = field.native_strain().at_barycentric(lambda).unwrap();
    near(native.membrane[0], free_strain, free_strain);
    near(native.membrane[1], free_strain, free_strain);
    let shear = field.at_barycentric(lambda, 0, 0.5).unwrap();
    for value in shear.facet_engineering_strain {
        near(value, 0.0, free_strain);
    }
    // Free isotropic expansion has positive epsilon_zz = free_strain. Applying
    // the mechanical closure to its total xx/yy would instead give a negative
    // epsilon_zz and understate the full tensor norm for this Poisson ratio.
    // This field does not retain that normal free strain, so it must refuse.
    assert!(matches!(
        field.isotropic_layer_bound(0, 0.25, &mut || false),
        Err(PlateError::BadSection { .. })
    ));
    assert!(matches!(
        field.isotropic_layer_bound(0, 0.25, &mut || true),
        Err(PlateError::Cancelled)
    ));

    // The same supplied displacement may be observed under the explicitly
    // mechanical zero-eigenstrain model, where the Poisson closure is valid.
    let prepared = model
        .prepare_laminate_shear(0, &source, 0.0, || false)
        .unwrap();
    let mechanical = prepared.evaluate(&u, || false).unwrap();
    let real_bound = mechanical
        .isotropic_layer_bound(0, 0.25, &mut || false)
        .unwrap();
    let complex_u: Vec<_> = u.iter().map(|&value| C64::new(value, 0.0)).collect();
    let harmonic = prepared.evaluate_harmonic(&complex_u, || false).unwrap();
    let harmonic_bound = harmonic
        .isotropic_layer_bound(0, 0.25, &mut || false)
        .unwrap();
    near(
        real_bound.frobenius_upper_bound,
        harmonic_bound.frobenius_upper_bound,
        free_strain,
    );
    assert!(real_bound.frobenius_upper_bound > free_strain);
}
