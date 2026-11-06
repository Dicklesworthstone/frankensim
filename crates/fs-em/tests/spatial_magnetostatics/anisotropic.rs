//! G1/G3: tensor flux, interface continuity, energy and physical-frame covariance.
use super::*;
use fs_material::{
    field_storage::spatial::{PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES, SpatialMagneticLaw},
    state_point::ResolvedMaterialStatePoint,
};

const AXES: [[f64; 3]; 3] = [
    [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
    [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0],
    [-2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0],
];

fn principal_material(values: [f64; 3]) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, value) in PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES
        .into_iter()
        .zip(values)
    {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, MAGNETIC_PERMEABILITY_DIMS),
                value: PropertyValue::Scalar {
                    value,
                    dims: MAGNETIC_PERMEABILITY_DIMS,
                },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic orthotropic permeability".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(
                key,
                MAGNETIC_PERMEABILITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic orthotropic magnetic material".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn tensor_flux(values: [f64; 3], axes: [[f64; 3]; 3], h: [f64; 3]) -> [f64; 3] {
    let principal_h = axes.map(|axis| dot(axis, h));
    std::array::from_fn(|j| {
        (0..3)
            .map(|i| axes[i][j] * values[i] * principal_h[i])
            .sum()
    })
}

fn orthotropic(n: usize, values: [f64; 3], axes: [[f64; 3]; 3]) -> SpatialMagnetostatics {
    let (positions, tets) = cube(n);
    let materials = vec![SpatialMagneticMaterial::Linear(principal_material(values)); tets.len()];
    let laws = vec![SpatialMagneticLaw::orthotropic(axes).unwrap(); tets.len()];
    SpatialMagnetostatics::with_magnetic_laws(
        positions,
        tets,
        materials,
        laws,
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn close(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 5e-9 * scale.max(b.abs()), "{a} != {b}");
}

#[test]
fn full_tensor_affine_fields_match_dirichlet_flux_and_source_substitution() {
    let h = [200.0, -150.0, 70.0];
    let values = [2e-6, 5e-6, 13e-6];
    for principal in [values, [3e-6, 7e-6, 11e-6]] {
        let problem = orthotropic(2, principal, AXES);
        assert_eq!(
            problem.magnetic_laws(),
            vec![SpatialMagneticLaw::orthotropic(AXES).unwrap(); problem.mesh().tets.len()]
        );
        assert_eq!(
            problem.materials()[0],
            SpatialMagneticMaterial::Linear(principal_material(principal))
        );
        assert_ne!(
            problem.permeability_coefficients()[0].matrix_h_m()[0][1],
            0.0
        );
        let expected_b = tensor_flux(principal, AXES, h);
        let boundary = SpatialMagneticBoundary {
            potential: vec![],
            flux: problem
                .exterior_faces()
                .into_iter()
                .map(|face| MagneticFluxFace {
                    face,
                    outward_flux_density_t: dot(expected_b, box_normal(&problem, face)),
                })
                .collect(),
        };
        for field in [
            problem
                .solve(
                    &potentials(&problem, |p| -dot(h, p)),
                    controls(),
                    &mut || false,
                )
                .unwrap(),
            problem.solve(&boundary, controls(), &mut || false).unwrap(),
        ] {
            for (actual_h, actual_b) in field
                .field_strength_a_m()
                .iter()
                .zip(field.flux_density_t())
            {
                for j in 0..3 {
                    close(actual_h[j], h[j], 200.0);
                    close(actual_b[j], expected_b[j], 0.003);
                }
            }
            close(field.stored_energy_j(), 0.5 * dot(h, expected_b), 1.0);
            assert!(
                field
                    .interface_flux_jumps_wb()
                    .iter()
                    .all(|(_, jump)| jump.abs() < 1e-11)
            );
        }
    }
}

#[test]
fn layered_tensor_materials_preserve_normal_flux_and_tangential_field() {
    let split = 0.4;
    let principal = [[2e-6, 5e-6, 13e-6], [7e-6, 3e-6, 9e-6]];
    let axes = [AXES, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]];
    let (bn, hy, hz) = (3e-4, 80.0, -35.0);
    let expected_h: [[f64; 3]; 2] = std::array::from_fn(|side| {
        let normal = tensor_flux(principal[side], axes[side], [1.0, 0.0, 0.0])[0];
        let transverse = tensor_flux(principal[side], axes[side], [0.0, hy, hz])[0];
        [(bn - transverse) / normal, hy, hz]
    });
    let expected_b: [[f64; 3]; 2] =
        std::array::from_fn(|side| tensor_flux(principal[side], axes[side], expected_h[side]));
    let (positions, tets) = grid(&[0.0, split, 1.0], &[0.0, 0.5, 1.0], &[0.0, 0.5, 1.0]);
    let side: Vec<_> = tets
        .iter()
        .map(|tet| {
            usize::from(
                tet.iter()
                    .map(|&v| positions[v as usize][0] / 4.0)
                    .sum::<f64>()
                    > split,
            )
        })
        .collect();
    let materials = side
        .iter()
        .map(|&i| SpatialMagneticMaterial::Linear(principal_material(principal[i])))
        .collect();
    let laws = side
        .iter()
        .map(|&i| SpatialMagneticLaw::orthotropic(axes[i]).unwrap())
        .collect();
    let problem = SpatialMagnetostatics::with_magnetic_laws(
        positions,
        tets,
        materials,
        laws,
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = problem
        .solve(
            &potentials(&problem, |p| {
                let i = usize::from(p[0] > split);
                -expected_h[i][0] * (p[0] - split) - hy * p[1] - hz * p[2]
            }),
            controls(),
            &mut || false,
        )
        .unwrap();
    for (cell, &i) in side.iter().enumerate() {
        for j in 0..3 {
            close(field.field_strength_a_m()[cell][j], expected_h[i][j], 200.0);
            close(field.flux_density_t()[cell][j], expected_b[i][j], 0.003);
        }
    }
    close(expected_b[0][0], bn, 0.003);
    close(expected_b[1][0], bn, 0.003);
    let energy = 0.5
        * (split * dot(expected_h[0], expected_b[0])
            + (1.0 - split) * dot(expected_h[1], expected_b[1]));
    close(field.stored_energy_j(), energy, 1.0);
    assert!(
        field
            .interface_flux_jumps_wb()
            .iter()
            .all(|(_, jump)| jump.abs() < 1e-11)
    );
}

#[test]
fn co_rotating_material_and_mesh_preserves_energy_and_rotates_the_actual_field() {
    let values = [2e-6, 5e-6, 13e-6];
    let h = [200.0, -150.0, 70.0];
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let original = orthotropic(2, values, AXES);
    let materials = original.materials().to_vec();
    let positions = original.positions_m().iter().copied().map(rotate).collect();
    let law = SpatialMagneticLaw::orthotropic(AXES.map(rotate)).unwrap();
    let rotated = SpatialMagnetostatics::with_magnetic_laws(
        positions,
        original.mesh().tets.clone(),
        materials,
        vec![law; original.mesh().tets.len()],
        limits(),
        &mut || false,
    )
    .unwrap();
    let first = original
        .solve(
            &potentials(&original, |p| -dot(h, p)),
            controls(),
            &mut || false,
        )
        .unwrap();
    let second = rotated
        .solve(
            &potentials(&rotated, |p| -dot(rotate(h), p)),
            controls(),
            &mut || false,
        )
        .unwrap();
    for (a, b) in first.flux_density_t().iter().zip(second.flux_density_t()) {
        for (expected, actual) in rotate(*a).into_iter().zip(b) {
            close(*actual, expected, 0.003);
        }
    }
    close(first.stored_energy_j(), second.stored_energy_j(), 1.0);
}

#[test]
fn original_isotropic_entrypoint_is_unchanged_and_recoil_never_acquires_an_unrelated_tensor() {
    let problem = linear(2, 2e-6);
    let explicit = SpatialMagnetostatics::with_magnetic_laws(
        problem.positions_m().to_vec(),
        problem.mesh().tets.clone(),
        problem.materials().to_vec(),
        vec![SpatialMagneticLaw::Isotropic; problem.materials().len()],
        limits(),
        &mut || false,
    )
    .unwrap();
    let phi = |p: [f64; 3]| -200.0 * p[0] + 150.0 * p[1] - 70.0 * p[2];
    let first = problem
        .solve(&potentials(&problem, phi), controls(), &mut || false)
        .unwrap();
    let second = explicit
        .solve(&potentials(&explicit, phi), controls(), &mut || false)
        .unwrap();
    assert_eq!(first.flux_density_t(), second.flux_density_t());
    assert_eq!(first.stored_energy_j(), second.stored_energy_j());
    let (p, t) = cube(1);
    let magnet = SpatialRecoilMaterial::new(recoil(2e-6, 0.9, 1e6), [1.0, 0.0, 0.0]).unwrap();
    let laws = vec![SpatialMagneticLaw::orthotropic(AXES).unwrap(); t.len()];
    let materials = vec![SpatialMagneticMaterial::Recoil(magnet); t.len()];
    assert!(
        SpatialMagnetostatics::with_magnetic_laws(p, t, materials, laws, limits(), &mut || false)
            .is_err()
    );
    let (p, t) = cube(1);
    assert!(
        SpatialMagnetostatics::with_magnetic_laws(p, t, vec![], vec![], limits(), &mut || false)
            .is_err()
    );
}

#[test]
fn anisotropic_exterior_refuses_both_maxwell_surface_and_eggshell_stress() {
    use fs_em::spatial::force::SpatialMaxwellForceLimits;
    use fs_ga::Se3;
    let problem = orthotropic(3, [2e-6, 5e-6, 13e-6], AXES);
    let field = problem
        .solve(
            &potentials(&problem, |p| -200.0 * p[0]),
            controls(),
            &mut || false,
        )
        .unwrap();
    let interior = |p: [f64; 3]| p.iter().all(|&v| v > 0.0 && v < 1.0);
    let selected: Vec<_> = problem
        .mesh()
        .tets
        .iter()
        .enumerate()
        .filter_map(|(i, t)| {
            t.iter()
                .all(|&v| interior(problem.positions_m()[v as usize]))
                .then_some(i)
        })
        .collect();
    assert_eq!(selected.len(), 6);
    let force_limits = SpatialMaxwellForceLimits {
        max_cells: 1000,
        max_faces: 3000,
    };
    assert!(matches!(
        field.maxwell_wrench(
            &selected,
            &Se3::identity(),
            [0.0; 3],
            force_limits,
            &mut || false
        ),
        Err(SpatialMagneticError::Invalid(
            "3D Maxwell surface requires isotropic exterior permeability"
        ))
    ));
    let weights: Vec<_> = problem
        .positions_m()
        .iter()
        .copied()
        .map(|p| if interior(p) { 1.0 } else { 0.0 })
        .collect();
    let shell = problem
        .prepare_eggshell(&weights, force_limits, &mut || false)
        .unwrap();
    assert!(matches!(
        field.eggshell_wrench(&shell, &Se3::identity(), [0.0; 3], &mut || false),
        Err(SpatialMagneticError::Invalid(
            "eggshell transition requires isotropic exterior permeability"
        ))
    ));
}
