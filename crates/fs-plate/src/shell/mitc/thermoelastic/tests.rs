use super::*;

fn section() -> ShearShellSection {
    ShearShellSection::from_orthotropic_solid(
        10.0, 8.0, 0.25, 3.0, 2.0, 2.5, 0.1, 2.0, 0.0, 5.0 / 6.0,
    ).unwrap()
}

fn flat() -> (ShellMesh, NodalDirectors) {
    let mesh = ShellMesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    ).unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.9).unwrap();
    (mesh, directors)
}

fn project(
    mesh: &ShellMesh, directors: &NodalDirectors,
    factors: &[Option<ShellThermoelasticFactors>], modes: &[Vec<f64>],
) -> ProjectedThermoelasticRows {
    project_thermoelastic_factors(
        mesh, directors, &vec![section(); mesh.tris.len()], factors, modes,
        MitcBudget { max_nodes: mesh.nodes.len(), max_triangles: mesh.tris.len() },
        1_000_000_000, 1_000_000, || false,
    ).unwrap()
}

fn affine(mesh: &ShellMesh, strain: [f64; 3]) -> Vec<f64> {
    mesh.nodes.iter().flat_map(|p| [
        strain[0] * p[0] + 0.5 * strain[2] * p[1],
        strain[1] * p[1] + 0.5 * strain[2] * p[0],
        0.0, 0.0, 0.0, 0.0,
    ]).collect()
}

fn gram(rows: &[f64], n: usize, i: usize, j: usize) -> f64 {
    rows.chunks_exact(n).map(|r| r[i] * r[j]).sum()
}

fn close(actual: f64, expected: f64, relative: f64) {
    if expected == 0.0 {
        assert_eq!(actual, 0.0);
    } else {
        assert!((actual / expected - 1.0).abs() <= relative,
            "{actual:e} vs {expected:e}");
    }
}

#[test]
fn g1_membrane_thermal_storage_preserves_all_signed_cross_terms_and_weak_modes() {
    let (mesh, directors) = flat();
    let factor = ShellThermoelasticFactors {
        membrane: [2.0, 3.0, 5.0], bending: [0.2, 0.3, 0.5],
    };
    let original = affine(&mesh, [1.0, -0.5, 0.25]);
    let amplitudes = [1.0, -2.0, 1e-120, 0.5];
    let modes: Vec<_> = amplitudes.iter()
        .map(|a| original.iter().map(|v| a * v).collect()).collect();
    let image = project(&mesh, &directors, &[Some(factor)], &modes);
    // Independent constant-strain area integral, not an element-matrix oracle.
    let expected = 0.5_f64 * (2.0_f64 - 1.5 + 1.25).powi(2);
    for i in 0..modes.len() {
        for j in 0..modes.len() {
            close(gram(&image.membrane, modes.len(), i, j),
                expected * amplitudes[i] * amplitudes[j], 4e-14);
        }
    }
    assert!(gram(&image.membrane, modes.len(), 0, 2) > 0.0);
    assert!(image.bending.iter().all(|v| *v == 0.0));
}

#[test]
fn g1_bending_memory_and_uniform_membrane_storage_remain_separate() {
    let (mesh, directors) = flat();
    let membrane = affine(&mesh, [0.5, 0.25, -0.125]);
    let curvature = [0.5, -0.25, 0.125];
    let bending: Vec<f64> = mesh.nodes.iter().flat_map(|p| [
        0.0, 0.0,
        0.5 * curvature[0] * p[0] * p[0] + 0.5 * curvature[1] * p[1] * p[1]
            + 0.5 * curvature[2] * p[0] * p[1],
        curvature[1] * p[1] + 0.5 * curvature[2] * p[0],
        -curvature[0] * p[0] - 0.5 * curvature[2] * p[1],
        0.0,
    ]).collect();
    let mixed = membrane.iter().zip(&bending).map(|(a, b)| a + b).collect();
    let modes = vec![membrane, bending, mixed];
    let factor = ShellThermoelasticFactors {
        membrane: [2.0, 3.0, 5.0], bending: [0.2, 0.3, 0.5],
    };
    let image = project(&mesh, &directors, &[Some(factor)], &modes);
    let expected_a = 0.5_f64 * (2.0_f64 * 0.5 + 3.0 * 0.25 - 5.0 * 0.125).powi(2);
    let expected_d = 0.5_f64 * (0.2_f64 * 0.5 - 0.3 * 0.25 + 0.5 * 0.125).powi(2);
    close(gram(&image.membrane, 3, 0, 0), expected_a, 5e-14);
    close(gram(&image.bending, 3, 1, 1), expected_d, 3e-11);
    close(gram(&image.membrane, 3, 0, 2), expected_a, 5e-14);
    close(gram(&image.bending, 3, 1, 2), expected_d, 3e-11);
    assert!(image.membrane.chunks_exact(3).all(|r| r[1] == 0.0));
    assert!(image.bending.chunks_exact(3).all(|r| r[0] == 0.0));
}

#[test]
fn g0_source_cancellation_keeps_weak_thermal_strain_and_exact_nulls() {
    let (mesh, directors) = flat();
    let weak = 2.0_f64.powi(-55);
    // The original physical strain is exactly [1,1,1]. Forming a dense
    // constitutive Gram first erases this weak factor contribution.
    let image = project(&mesh, &directors, &[Some(ShellThermoelasticFactors {
        membrane: [1.0, weak, -1.0], bending: [0.0; 3],
    })], &[affine(&mesh, [1.0; 3])]);
    close(gram(&image.membrane, 1, 0, 0) / (weak * weak), 0.5, 4e-14);
    let null = project(&mesh, &directors, &[Some(ShellThermoelasticFactors {
        membrane: [1.0, 1.0, 0.0], bending: [1.0, 1.0, 0.0],
    })], &[affine(&mesh, [1.0, -1.0, 0.0]), affine(&mesh, [0.0, 0.0, 1.0])]);
    assert!(null.membrane.iter().chain(&null.bending).all(|v| *v == 0.0));
    let mut original = [0.0; 20];
    original[0] = 1.0;
    original[9] = weak;
    original[19] = -1.0;
    assert_eq!(compensated_sum(original.into_iter()), weak);
}

#[test]
fn g3_curved_directors_preserve_original_bubble_and_rigid_frame_covariance() {
    let mesh = ShellMesh::new(
        vec![[0.0; 3], [0.8, 0.1, 0.04], [0.05, 0.6, -0.02]],
        vec![[0, 1, 2]],
    ).unwrap();
    let geometry = mesh.facet(0).unwrap();
    let vectors: Vec<_> = (0..3).map(|i| {
        let v = [
            geometry.frame[2][0] + 0.02 * (i as f64 - 1.0),
            geometry.frame[2][1] + 0.01 * (i as f64 - 0.5),
            geometry.frame[2][2],
        ];
        let norm = v.iter().map(|v| v * v).sum::<f64>().sqrt();
        v.map(|x| x / norm)
    }).collect();
    let directors = NodalDirectors::new(&mesh, vectors.clone(), 0.9).unwrap();
    let modes: Vec<Vec<f64>> = (0..3).map(|k| (0..18)
        .map(|j| (((j * (7 + k) + k) % 13) as f64 - 6.0) * 0.01).collect()).collect();
    let factor = ShellThermoelasticFactors { membrane: [2.0, 3.0, 5.0], bending: [0.2, -0.3, 0.5] };
    let original = project(&mesh, &directors, &[Some(factor)], &modes);
    let field = element::PreparedStrainField::new(
        &mesh, &mesh.tris[0], geometry, mesh.tris[0].map(|i| vectors[i]),
        &section(), &mut || false,
    ).unwrap();
    let mut expected = [[0.0_f64; 9]; 2];
    for (lambda, weight) in element::stiffness_rule() {
        let channels: Vec<[f64; 2]> = modes.iter().map(|p| {
            let nodal: [f64; 18] = core::array::from_fn(|j| p[j]);
            let (a, d, _) = field.eval(&field.expand(&nodal), lambda[1], lambda[2]);
            [(0..3).map(|k| factor.membrane[k] * a[k]).sum(),
                (0..3).map(|k| factor.bending[k] * d[k]).sum()]
        }).collect();
        for part in 0..2 { for i in 0..3 { for j in 0..3 {
            expected[part][3 * i + j] += geometry.area_m2 * weight * channels[i][part] * channels[j][part];
        }}}
    }
    for (part, actual) in [&original.membrane, &original.bending].into_iter().enumerate() {
        for i in 0..3 { for j in 0..3 {
            close(gram(actual, 3, i, j), expected[part][3 * i + j], 1e-11);
        }}
    }
    // Proper cyclic world-frame rotation leaves every local material axis,
    // strain factor and physical nodal map covariant.
    let rotate = |p: [f64; 3]| [p[2], p[0], p[1]];
    let rotated_mesh = ShellMesh::new(
        mesh.nodes.iter().map(|p| rotate(*p)).collect(), mesh.tris.clone(),
    ).unwrap();
    let rotated_directors = NodalDirectors::new(
        &rotated_mesh, vectors.into_iter().map(rotate).collect(), 0.9,
    ).unwrap();
    let rotated_modes: Vec<Vec<f64>> = modes.iter().map(|p| p.chunks_exact(3)
        .flat_map(|v| rotate([v[0], v[1], v[2]])).collect()).collect();
    let rotated = project(&rotated_mesh, &rotated_directors, &[Some(factor)], &rotated_modes);
    for (a, b) in original.membrane.iter().chain(&original.bending)
        .zip(rotated.membrane.iter().chain(&rotated.bending)) {
        assert!((a - b).abs() < 2e-12, "{a:e} vs {b:e}");
    }

    let rigid: Vec<Vec<f64>> = (0..6).map(|axis| mesh.nodes.iter().flat_map(|p| {
        let mut translation = [0.0; 3];
        let mut rotation = [0.0; 3];
        if axis < 3 { translation[axis] = 1.0; } else {
            rotation[axis - 3] = 1.0;
            translation = crate::shell::cross(rotation, *p);
        }
        [translation[0], translation[1], translation[2], rotation[0], rotation[1], rotation[2]]
    }).collect()).collect();
    let null = project(&mesh, &directors, &[Some(factor)], &rigid);
    assert!(null.membrane.iter().chain(&null.bending).all(|v| v.abs() < 2e-12));
}

#[test]
fn g0_projection_keeps_the_complete_basis_and_refuses_partial_resources() {
    let (mesh, directors) = flat();
    let zero = Some(ShellThermoelasticFactors { membrane: [0.0; 3], bending: [0.0; 3] });
    let modes = vec![vec![0.0; 18]; 65];
    let sections = [section()];
    let budget = MitcBudget { max_nodes: 3, max_triangles: 1 };
    let work = projection_terms(3, 1, 1, 65).unwrap();
    let entries = 2 * 65 * 65;
    let run = |work, entries, cancelled| project_thermoelastic_factors(
        &mesh, &directors, &sections, &[zero], &modes, budget, work, entries, || cancelled,
    );
    let image = run(work, entries, false).unwrap();
    assert_eq!(image.mode_count, 65);
    assert_eq!(image.row_count, 65);
    assert_eq!(&image.row_elements[..7], &[Some(0); 7]);
    assert!(image.row_elements[7..].iter().all(Option::is_none));
    assert!(image.membrane.iter().chain(&image.bending).all(|v| *v == 0.0));
    assert!(run(work - 1, entries, false).is_err());
    assert!(run(work, entries - 1, false).is_err());
    assert!(matches!(run(work, entries, true), Err(PlateError::Cancelled)));
    assert!(project_thermoelastic_factors(
        &mesh, &directors, &sections, &[zero],
        &vec![Vec::new(); MAX_THERMOELASTIC_MODES + 1], budget,
        usize::MAX, usize::MAX, || false,
    ).is_err());
    assert!(projection_terms(usize::MAX, 1, 1, 1).is_none());
    assert!(projection_terms(3, 1, 2, 1).is_none());
    let absent = project(&mesh, &directors, &[None], &[affine(&mesh, [1.0; 3])]);
    assert!(absent.membrane.iter().chain(&absent.bending).all(|v| *v == 0.0));
    assert!(absent.row_elements.iter().all(Option::is_none));
}
