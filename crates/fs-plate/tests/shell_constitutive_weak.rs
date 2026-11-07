//! Physical weak strain work in the original dimensional constitutive tensor.

use fs_plate::PlateSection;
use fs_plate::shell::viscosity::{ShellViscosity, ShellViscosityBudget, project_viscosity};
use fs_plate::shell::{ShellMesh, ShellModel, ShellSupport, assemble_shell};

fn fixture() -> (ShellMesh, ShellModel) {
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let section = PlateSection::isotropic(1e7, 0.25, 0.01, 800.0).unwrap();
    let model = assemble_shell(&mesh, &section, &[], ShellSupport::Free).unwrap();
    (mesh, model)
}

fn power(
    mesh: &ShellMesh,
    model: &ShellModel,
    law: ShellViscosity,
    membrane: [f64; 3],
    curvature: [f64; 3],
) -> f64 {
    let shape: Vec<_> = mesh
        .nodes
        .iter()
        .flat_map(|&[x, y, _]| {
            [
                membrane[0] * x + 0.5 * membrane[2] * y,
                membrane[1] * y + 0.5 * membrane[2] * x,
                -0.5 * curvature[0] * x * x
                    - 0.5 * curvature[1] * y * y
                    - 0.5 * curvature[2] * x * y,
                -curvature[1] * y - 0.5 * curvature[2] * x,
                curvature[0] * x + 0.5 * curvature[2] * y,
                0.0,
            ]
        })
        .collect();
    project_viscosity(
        mesh,
        &[law],
        model,
        &[shape],
        ShellViscosityBudget {
            max_elements: 1,
            max_full_dofs: 18,
            max_dofs: 18,
            max_modes: 1,
            max_work: 1_000_000,
        },
        || false,
    )
    .unwrap()[0]
}

#[test]
fn nonbinary_weak_membrane_and_offset_coupling_keep_original_source_work() {
    let (mesh, model) = fixture();
    for a in [1.0, 3.7, 3.7e120, 3.7e-120] {
        let b = a + a * 2.0_f64.powi(-48);
        let expected = 0.5 * (b - a); // original area times original source strain power
        let membrane =
            ShellViscosity::new([a, a, 0.0, a, b, 0.0, 0.0, 0.0, 0.0], [0.0; 9]).unwrap();
        let actual = power(&mesh, &model, membrane, [-1.0, 1.0, 0.0], [0.0; 3]);
        assert!(
            (actual / expected - 1.0).abs() < 2e-12,
            "membrane {actual:e} != {expected:e}"
        );
        for sign in [-1.0, 1.0] {
            let coupled = ShellViscosity::new_coupled(
                [a, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [sign * a, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [b, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            )
            .unwrap();
            let actual = power(&mesh, &model, coupled, [-sign, 0.0, 0.0], [1.0, 0.0, 0.0]);
            assert!(
                (actual / expected - 1.0).abs() < 2e-12,
                "coupling {actual:e} != {expected:e}"
            );
        }
    }
}

#[test]
fn exact_nulls_and_independent_small_channels_do_not_follow_strong_scale() {
    let (mesh, model) = fixture();
    for a in [3.7, 1e200, 1e308] {
        let law = ShellViscosity::new([a, a, 0.0, a, a, 0.0, 0.0, 0.0, 1e-100], [0.0; 9]).unwrap();
        assert_eq!(
            power(&mesh, &model, law.clone(), [-1.0, 1.0, 0.0], [0.0; 3]),
            0.0
        );
        let weak = power(&mesh, &model, law, [0.0, 0.0, 1.0], [0.0; 3]);
        assert!(
            (weak / (0.5e-100) - 1.0).abs() < 2e-14,
            "weak source {weak:e}"
        );
    }
    assert!(ShellViscosity::new([1.0, 2.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0], [0.0; 9]).is_err());
    assert!(ShellViscosity::new([0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0], [0.0; 9]).is_err());
}
