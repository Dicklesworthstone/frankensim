//! G1/G3 synthetic tensor weak forms; no experimental dielectric calibration.
use fs_em::electrostatic::tetrahedral::{
    TetDielectricBudget, TetDielectricError, TetDielectricOperator,
};
use fs_feec::fixtures::{kuhn_cube, single_tet};
use fs_rep_mesh::TetComplex;

fn budget() -> TetDielectricBudget {
    TetDielectricBudget {
        maximum_nodes: 1000,
        maximum_tetrahedra: 1000,
        minimum_scaled_jacobian: 1e-12,
    }
}
fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * scale.max(expected.abs()),
        "{actual} != {expected}"
    );
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn epsilon() -> [[f64; 3]; 3] {
    [
        [4e-8, 1e-8, -0.5e-8],
        [1e-8, 3e-8, 0.25e-8],
        [-0.5e-8, 0.25e-8, 2e-8],
    ]
}

#[test]
fn affine_anisotropic_field_matches_native_dielectric_work() {
    let (mesh, positions) = kuhn_cube(2);
    let tensor = epsilon();
    let operator = TetDielectricOperator::assemble(
        &mesh,
        &positions,
        &vec![tensor; mesh.tets.len()],
        budget(),
        &mut || false,
    )
    .unwrap();
    let potential: Vec<_> = positions
        .iter()
        .map(|p| 2.0 * p[0] - 4.0 * p[1] + p[2])
        .collect();
    let fields = operator
        .electric_field_v_m(&potential, &mut || false)
        .unwrap();
    let expected = [-2.0, 4.0, -1.0];
    for field in &fields {
        for j in 0..3 {
            near(field[j], expected[j], 1.0);
        }
    }
    let mut reaction = vec![0.0; potential.len()];
    operator.stiffness().spmv(&potential, &mut reaction);
    let expected_ed = dot(&expected, &tensor.map(|row| dot(&row, &expected)));
    near(dot(&potential, &reaction), expected_ed, 1e-8);
    near(operator.volumes_m3().iter().sum(), 1.0, 1.0);
    near(reaction.iter().sum(), 0.0, 1e-8);
    // Potential offsets disappear before field multiplication, including in
    // anisotropic materials. These fixture values remain exactly represented.
    let shifted: Vec<_> = potential.iter().map(|v| 1e12 + v).collect();
    assert_eq!(
        fields,
        operator
            .electric_field_v_m(&shifted, &mut || false)
            .unwrap()
    );
}

#[test]
fn polarization_weak_form_and_charge_reaction_have_reciprocal_signs() {
    let (mesh, positions) = single_tet();
    let operator =
        TetDielectricOperator::assemble(&mesh, &positions, &[epsilon()], budget(), &mut || false)
            .unwrap();
    let potential = [0.0, 2.0, -4.0, 1.0];
    let polarization = [[0.02, -0.03, 0.04]];
    let load = operator
        .polarization_load_c(&polarization, &mut || false)
        .unwrap();
    let e = operator
        .electric_field_v_m(&potential, &mut || false)
        .unwrap();
    let d = operator
        .displacement_c_m2(&potential, &polarization, &mut || false)
        .unwrap();
    let volume = operator.volumes_m3()[0];
    near(
        dot(&potential, &load),
        -volume * dot(&e[0], &polarization[0]),
        0.1,
    );
    near(load.iter().sum(), 0.0, 0.1);
    let mut reaction = vec![0.0; 4];
    operator.stiffness().spmv(&potential, &mut reaction);
    for a in 0..4 {
        near(
            reaction[a] - load[a],
            -volume * dot(&operator.gradients()[0][a], &d[0]),
            0.1,
        );
    }
    // A uniformly polarized open dielectric has E=-epsilon^-1 P. Choosing
    // P=-epsilon E here makes D=0 and every nodal free-charge reaction vanish.
    let balancing = [epsilon().map(|row| -dot(&row, &e[0]))];
    let balanced = operator
        .polarization_load_c(&balancing, &mut || false)
        .unwrap();
    for a in 0..4 {
        near(reaction[a], balanced[a], 1e-8);
    }
    for component in operator
        .displacement_c_m2(&potential, &balancing, &mut || false)
        .unwrap()[0]
    {
        near(component, 0.0, 1e-8);
    }
}

#[test]
fn reversed_cell_orientation_preserves_fields_and_nodal_operator() {
    let (mesh, positions) = single_tet();
    let reversed = TetComplex::from_tets(4, vec![[0, 2, 1, 3]]);
    let a =
        TetDielectricOperator::assemble(&mesh, &positions, &[epsilon()], budget(), &mut || false)
            .unwrap();
    let b =
        TetDielectricOperator::assemble(&reversed, &positions, &[epsilon()], budget(), &mut || {
            false
        })
        .unwrap();
    let voltage = [0.0, 2.0, -4.0, 1.0];
    assert_eq!(a.volumes_m3(), b.volumes_m3());
    assert_eq!(
        a.electric_field_v_m(&voltage, &mut || false).unwrap(),
        b.electric_field_v_m(&voltage, &mut || false).unwrap()
    );
    let mut av = vec![0.0; 4];
    let mut bv = vec![0.0; 4];
    a.stiffness().spmv(&voltage, &mut av);
    b.stiffness().spmv(&voltage, &mut bv);
    for i in 0..4 {
        near(av[i], bv[i], 1e-8);
    }
}

#[test]
fn invalid_tensors_mesh_bounds_and_cancelled_work_refuse() {
    let (mesh, positions) = single_tet();
    let assemble = |tensor, bounds| {
        TetDielectricOperator::assemble(&mesh, &positions, &[tensor], bounds, &mut || false)
    };
    let mut nonsymmetric = epsilon();
    nonsymmetric[0][1] = 0.0;
    assert!(assemble(nonsymmetric, budget()).is_err());
    let mut indefinite = epsilon();
    indefinite[2][2] = -1.0;
    assert!(assemble(indefinite, budget()).is_err());
    assert!(
        assemble(
            epsilon(),
            TetDielectricBudget {
                maximum_nodes: 3,
                ..budget()
            }
        )
        .is_err()
    );
    assert!(
        assemble(
            epsilon(),
            TetDielectricBudget {
                maximum_tetrahedra: 0,
                ..budget()
            }
        )
        .is_err()
    );
    assert!(
        assemble(
            epsilon(),
            TetDielectricBudget {
                maximum_nodes: 4,
                maximum_tetrahedra: 1,
                ..budget()
            }
        )
        .is_ok()
    );
    assert!(matches!(
        TetDielectricOperator::assemble(&mesh, &positions, &[epsilon()], budget(), &mut || true),
        Err(TetDielectricError::Cancelled)
    ));
    let flat = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]];
    assert!(
        TetDielectricOperator::assemble(&mesh, &flat, &[epsilon()], budget(), &mut || false)
            .is_err()
    );
    let isolated = TetComplex::from_tets(5, mesh.tets.clone());
    let mut extra = positions.clone();
    extra.push([1., 1., 1.]);
    assert!(
        TetDielectricOperator::assemble(&isolated, &extra, &[epsilon()], budget(), &mut || false)
            .is_err()
    );
    let operator = assemble(epsilon(), budget()).unwrap();
    assert!(operator.polarization_load_c(&[], &mut || false).is_err());
    assert!(
        operator
            .electric_field_v_m(&[0.0; 3], &mut || false)
            .is_err()
    );
    assert!(matches!(
        operator.polarization_load_c(&[[0.0; 3]], &mut || true),
        Err(TetDielectricError::Cancelled)
    ));
}
