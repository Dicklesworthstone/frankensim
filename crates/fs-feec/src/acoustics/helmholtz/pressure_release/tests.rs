use super::*;
use crate::acoustics::helmholtz::Medium;
use crate::acoustics::helmholtz::heterogeneous::{HeterogeneousError, HeterogeneousHelmholtzFe};
use crate::fixtures::{masked_cube_grid, single_tet, two_tets};
use fs_rep_mesh::TetComplex;

const MEDIUM: Medium = Medium {
    density_kg_m3: 1.3,
    sound_speed_m_s: 2.0,
};

fn brick(cells: [usize; 3], lengths: [f64; 3]) -> (TetComplex, Vec<[f64; 3]>) {
    let (complex, mut positions) = masked_cube_grid(cells[0], cells[1], cells[2], &|_, _, _| true);
    for point in &mut positions {
        for axis in 0..3 {
            point[axis] *= lengths[axis] / cells[axis] as f64;
        }
    }
    (complex, positions)
}

fn plane(fe: &HelmholtzFe<'_>, axis: usize, value: f64) -> Vec<[u32; 3]> {
    fe.boundary_faces()
        .into_iter()
        .filter(|face| {
            face.iter()
                .all(|&v| (fe.positions()[v as usize][axis] - value).abs() < 1e-12)
        })
        .collect()
}

fn uniform<'f, 'c>(fe: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        fe,
        fe.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}

fn normalized(matrix: &Csr, shapes: &[Vec<f64>]) {
    for (i, shape) in shapes.iter().enumerate() {
        let mut applied = vec![0.0; matrix.nrows()];
        matrix.spmv(shape, &mut applied);
        for (j, other) in shapes.iter().enumerate() {
            let pairing: f64 = other.iter().zip(&applied).map(|(a, b)| a * b).sum();
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (pairing - expected).abs() < 2e-9,
                "mass pair {i}/{j}: {pairing}"
            );
        }
    }
}

#[test]
fn pressure_release_duct_modes_converge_to_quarter_and_half_waves() {
    // c=2 and L=1 give (n+1/2) Hz for one release end and (n+1) Hz
    // for two. Transverse rigid modes lie above the chosen windows.
    for both_ends in [false, true] {
        let exact: Vec<_> = (0..3)
            .map(|mode| mode as f64 + if both_ends { 1.0 } else { 0.5 })
            .collect();
        let mut first_errors = Vec::new();
        for cells in [4, 8, 16] {
            let (complex, positions) = brick([cells, 1, 1], [1.0, 0.17, 0.13]);
            let fe = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
            let mut faces = plane(&fe, 0, 1.0);
            if both_ends {
                faces.extend(plane(&fe, 0, 0.0));
            }
            let result = fe
                .modes_with_pressure_release(
                    (0.1, exact[2] + 0.4),
                    &faces,
                    &SliceOptions::default(),
                    &mut || false,
                )
                .unwrap();
            let modes = &result.eigenpairs;
            assert_eq!(modes.frequency_hz.len(), 3);
            assert_eq!(modes.below_low, 0);
            assert_eq!(modes.below_high, 3);
            assert!(result.free_dofs < fe.ndof());
            let errors: Vec<_> = modes
                .frequency_hz
                .iter()
                .zip(&exact)
                .map(|(actual, expected)| (actual / expected - 1.0).abs())
                .collect();
            first_errors.push(errors[0]);
            eprintln!("pressure-release ends={both_ends} nx={cells} frequency_errors={errors:?}");
            if cells == 16 {
                assert!(errors.iter().all(|error| *error < 2e-3));
                normalized(fe.mass(), &modes.shapes);
                let constraint = fe
                    .pressure_release_constraints(&faces, EmptyFreePressure::Refuse, &mut || false)
                    .unwrap();
                for shape in &modes.shapes {
                    for (dof, &fixed) in constraint.constrained.iter().enumerate() {
                        if fixed {
                            assert_eq!(shape[dof].to_bits(), 0.0_f64.to_bits());
                        }
                    }
                }
            }
        }
        for pair in first_errors.windows(2) {
            assert!(
                pair[1] < pair[0] / 6.0,
                "P2 frequency convergence: {first_errors:?}"
            );
        }
    }
}

#[test]
fn mixed_box_pressure_release_keeps_three_dimensional_modes_and_counts() {
    let lengths = [1.0, 0.8, 0.6];
    let (complex, positions) = brick([2, 2, 2], lengths);
    let fe = HelmholtzFe::new(&complex, &positions, 3, MEDIUM).unwrap();
    let mut faces = plane(&fe, 0, 0.0);
    faces.extend(plane(&fe, 0, lengths[0]));
    faces.extend(plane(&fe, 1, 0.0));
    // x has two release ends; y has one release and one rigid end; z is rigid.
    let mut exact = Vec::new();
    for x in 1..=4 {
        for y in 0..4 {
            for z in 0..4 {
                exact.push(
                    ((x as f64 / lengths[0]).powi(2)
                        + ((y as f64 + 0.5) / lengths[1]).powi(2)
                        + (z as f64 / lengths[2]).powi(2))
                    .sqrt(),
                );
            }
        }
    }
    exact.sort_by(f64::total_cmp);
    let high = 0.5 * (exact[3] + exact[4]);
    let result = fe
        .modes_with_pressure_release((0.1, high), &faces, &SliceOptions::default(), &mut || false)
        .unwrap();
    assert_eq!(
        result.eigenpairs.below_high - result.eigenpairs.below_low,
        4
    );
    assert_eq!(result.eigenpairs.frequency_hz.len(), 4);
    for (actual, expected) in result.eigenpairs.frequency_hz.iter().zip(exact) {
        assert!(
            (actual / expected - 1.0).abs() < 0.02,
            "{actual} vs {expected}"
        );
    }
    normalized(fe.mass(), &result.eigenpairs.shapes);
}

#[test]
fn density_jump_pressure_release_modes_use_the_original_physical_pencil() {
    let (complex, positions) = brick([12, 1, 1], [1.0, 0.17, 0.13]);
    let fe = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let faces = plane(&fe, 0, 1.0);
    let options = SliceOptions {
        mass_diagonal_equilibration: true,
        ..SliceOptions::default()
    };
    let homogeneous = fe
        .modes_with_pressure_release((0.1, 0.8), &faces, &options, &mut || false)
        .unwrap();
    let uniform_model = uniform(&fe);
    let uniform_modes = uniform_model
        .modes_with_pressure_release((0.1, 0.8), &faces, &options, &mut || false)
        .unwrap();
    assert_eq!(uniform_modes.eigenpairs.frequency_hz.len(), 1);
    assert!(
        (uniform_modes.eigenpairs.frequency_hz[0] / homogeneous.eigenpairs.frequency_hz[0] - 1.0)
            .abs()
            < 2e-9
    );
    for (left, right) in [(1.0_f64, 4.0_f64), (4.0, 1.0)] {
        let model = HeterogeneousHelmholtzFe::new(
            &fe,
            4,
            &mut |_, position| {
                Ok::<_, AcousticFeError>(Medium {
                    density_kg_m3: if position[0] < 0.5 { left } else { right },
                    sound_speed_m_s: MEDIUM.sound_speed_m_s,
                })
            },
            &mut || false,
        )
        .unwrap();
        // p=A cos(k x) on the rigid half and p=B sin(k(1-x)) on
        // the release half. Pressure and p'/rho continuity at x=1/2
        // give tan(k/2)^2=rho_left/rho_right. The smallest positive root:
        let exact = MEDIUM.sound_speed_m_s * (left / right).sqrt().atan() / core::f64::consts::PI;
        let result = model
            .modes_with_pressure_release((0.1, 0.9), &faces, &options, &mut || false)
            .unwrap();
        assert_eq!(result.eigenpairs.frequency_hz.len(), 1);
        let actual = result.eigenpairs.frequency_hz[0];
        assert!(
            (actual / exact - 1.0).abs() < 2e-4,
            "{left}/{right}: {actual} vs {exact}"
        );
        normalized(model.mass(), &result.eigenpairs.shapes);
        assert_eq!(result.free_dofs, homogeneous.free_dofs);
        assert_eq!(
            result.pressure_release_faces,
            homogeneous.pressure_release_faces
        );
    }
}

#[test]
fn every_hierarchical_trace_dof_is_released_independently_of_face_winding() {
    let (complex, positions) = single_tet();
    for order in 1..=4 {
        let fe = HelmholtzFe::new(&complex, &positions, order, MEDIUM).unwrap();
        let constraints = fe
            .pressure_release_constraints(&[[0, 1, 2]], EmptyFreePressure::Refuse, &mut || false)
            .unwrap();
        assert_eq!(
            constraints
                .constrained
                .iter()
                .filter(|&&value| value)
                .count(),
            (order + 1) * (order + 2) / 2,
        );
        let mut coefficients = vec![0.0; fe.ndof()];
        for &free in &constraints.free_dofs {
            coefficients[free] = (free + 1) as f64;
        }
        for barycentric in [
            [0.13, 0.37, 0.5, 0.0],
            [0.5, 0.5, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ] {
            let trace: f64 = fe
                .space
                .element_dofs(0)
                .iter()
                .map(|(dof, shape)| coefficients[*dof] * shape.eval(barycentric, order))
                .sum();
            assert_eq!(trace, 0.0);
        }
        for face in [[2, 1, 0], [1, 2, 0], [0, 2, 1]] {
            let permuted = fe
                .pressure_release_constraints(&[face], EmptyFreePressure::Refuse, &mut || false)
                .unwrap();
            assert_eq!(permuted.free_dofs, constraints.free_dofs);
            assert_eq!(permuted.constrained, constraints.constrained);
            assert_eq!(permuted.faces, constraints.faces);
        }
    }
}

#[test]
fn release_refusals_are_transactional_and_empty_faces_preserve_rigid_modes() {
    let (complex, positions) = brick([2, 1, 1], [1.0, 0.17, 0.13]);
    let fe = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&fe);
    let window = (-0.1, 1.8);
    let options = SliceOptions::default();
    let rigid = fe.modes(window, &options).unwrap();
    let identity = fe
        .modes_with_pressure_release(window, &[], &options, &mut || false)
        .unwrap();
    assert_eq!(identity.eigenpairs.lambda, rigid.lambda);
    assert_eq!(identity.eigenpairs.shapes, rigid.shapes);
    assert_eq!(identity.eigenpairs.residual, rigid.residual);
    assert_eq!(identity.free_dofs, fe.ndof());
    assert!(identity.pressure_release_faces.is_empty());
    let rigid = model.modes(window, &options, &mut || false).unwrap();
    let identity = model
        .modes_with_pressure_release(window, &[], &options, &mut || false)
        .unwrap();
    assert_eq!(
        identity.eigenpairs.omega_squared_rad2_s2,
        rigid.omega_squared_rad2_s2
    );
    assert_eq!(identity.eigenpairs.shapes, rigid.shapes);
    assert_eq!(identity.eigenpairs.residual, rigid.residual);
    let faces = plane(&fe, 0, 1.0);
    assert!(matches!(
        fe.modes_with_pressure_release(window, &faces, &options, &mut || true),
        Err(AcousticFeError::Cancelled)
    ));
    assert!(matches!(
        model.modes_with_pressure_release(window, &faces, &options, &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    let mut checkpoints = 0;
    assert!(matches!(
        fe.pressure_release_constraints(&faces, EmptyFreePressure::Refuse, &mut || {
            checkpoints += 1;
            checkpoints == 4
        }),
        Err(AcousticFeError::Cancelled)
    ));
    let face = faces[0];
    assert!(matches!(
        fe.pressure_release_constraints(
            &[face, [face[2], face[1], face[0]]],
            EmptyFreePressure::Refuse,
            &mut || false
        ),
        Err(AcousticFeError::InvalidInput { .. })
    ));
    let (complex, positions) = two_tets();
    let fe = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    assert!(matches!(
        fe.pressure_release_constraints(&[[1, 2, 3]], EmptyFreePressure::Refuse, &mut || false),
        Err(AcousticFeError::FaceNotOnBoundary { .. })
    ));
    assert!(matches!(
        fe.pressure_release_constraints(&[[0, 1, 999]], EmptyFreePressure::Refuse, &mut || false),
        Err(AcousticFeError::FaceNotInMesh { .. })
    ));
    let (complex, positions) = single_tet();
    let fe = HelmholtzFe::new(&complex, &positions, 1, MEDIUM).unwrap();
    assert!(matches!(
        fe.pressure_release_constraints(
            &fe.boundary_faces(),
            EmptyFreePressure::Refuse,
            &mut || false
        ),
        Err(AcousticFeError::InvalidInput { .. })
    ));
}

// The empty free pressure set is admitted only for a caller that carries wall states of its own; every other caller, and the
// default of the helper before this variant existed, still refuses it.
#[test]
fn an_empty_free_set_is_admitted_only_when_the_caller_names_wall_states() {
    let (complex, positions) = single_tet();
    let fe = HelmholtzFe::new(&complex, &positions, 1, MEDIUM).unwrap();
    let all = fe.boundary_faces();
    let admitted = fe
        .pressure_release_constraints(&all, EmptyFreePressure::AdmitWithWallStates, &mut || false)
        .unwrap();
    assert!(admitted.free_dofs.is_empty());
    assert_eq!(
        admitted.constrained.iter().filter(|&&c| c).count(),
        fe.ndof()
    );
    assert!(matches!(
        fe.pressure_release_constraints(&all, EmptyFreePressure::Refuse, &mut || false),
        Err(AcousticFeError::InvalidInput { .. })
    ));
}
