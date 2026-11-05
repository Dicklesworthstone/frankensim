use fs_material::laminate::{
    LaminaLayer, LaminateSection, LaminateShearRecovery, LayerLaw, ReferenceSurface,
};
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
        let (_, energy) = integrate(&field, &source, angle, mesh.facet(0).unwrap().area_m2);
        near(energy, model.strain_energy(&u).unwrap(), energy);
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
