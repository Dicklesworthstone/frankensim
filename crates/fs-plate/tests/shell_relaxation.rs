//! Hereditary shell resultants checked against independent physical strain work.

use fs_plate::shell::relaxation::{
    ShellRelaxationBranch, ShellRelaxationBudget, project_relaxation,
};
use fs_plate::shell::viscosity::ShellViscosityBudget;
use fs_plate::shell::{ShellMesh, ShellModel, ShellSupport, assemble_shell};
use fs_plate::{PlateError, PlateSection};

fn fixture() -> (ShellMesh, ShellModel) {
    let mesh = ShellMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
        ],
        vec![[0, 1, 2], [3, 4, 5]],
    )
    .unwrap();
    let section = PlateSection::isotropic(1e7, 0.25, 0.01, 800.0).unwrap();
    let model = assemble_shell(&mesh, &section, &[], ShellSupport::Free).unwrap();
    (mesh, model)
}

fn budget() -> ShellRelaxationBudget {
    ShellRelaxationBudget {
        projection: ShellViscosityBudget {
            max_elements: 2,
            max_full_dofs: 36,
            max_dofs: 36,
            max_modes: 8,
            max_work: 10_000_000,
        },
        max_element_branches: 8,
        max_relaxation_times: 4,
    }
}

fn patch(mesh: &ShellMesh, e: [f64; 3], k: [f64; 3]) -> Vec<f64> {
    mesh.nodes
        .iter()
        .flat_map(|p| {
            let [x, y, _] = *p;
            [
                e[0] * x + 0.5 * e[2] * y,
                e[1] * y + 0.5 * e[2] * x,
                -0.5 * k[0] * x * x - 0.5 * k[1] * y * y - 0.5 * k[2] * x * y,
                -k[1] * y - 0.5 * k[2] * x,
                k[0] * x + 0.5 * k[2] * y,
                0.0,
            ]
        })
        .collect()
}

fn bilinear(a: [f64; 9], x: [f64; 3], y: [f64; 3]) -> f64 {
    (0..3)
        .map(|i| (0..3).map(|j| x[i] * a[3 * i + j] * y[j]).sum::<f64>())
        .sum()
}

#[test]
fn separate_facet_spectra_preserve_signed_coupling_and_exact_equal_times() {
    let (mesh, model) = fixture();
    let a = [8.0, 1.0, 0.0, 1.0, 6.0, 0.0, 0.0, 0.0, 4.0];
    let b = [0.5, 0.1, 0.0, 0.1, 0.2, 0.0, 0.0, 0.0, 0.3];
    let d = [5.0, 0.2, 0.0, 0.2, 3.0, 0.0, 0.0, 0.0, 2.0];
    let strain = [
        ([0.3, -0.2, 0.1], [0.4, 0.2, -0.3]),
        ([-0.1, 0.5, -0.2], [-0.3, 0.1, 0.4]),
    ];
    let modes: Vec<_> = strain.iter().map(|&(e, k)| patch(&mesh, e, k)).collect();
    for sign in [1.0, -1.0] {
        let b = b.map(|x| sign * x);
        let branch = |scale: f64, time| {
            ShellRelaxationBranch::new_coupled(
                a.map(|x| scale * x),
                b.map(|x| scale * x),
                d.map(|x| scale * x),
                time,
            )
            .unwrap()
        };
        let sections = vec![
            vec![branch(1.0, 0.002), branch(0.25, 0.04)],
            vec![
                branch(2.0, 0.002),
                ShellRelaxationBranch::new([0.0; 9], [0.0; 9], 99.0).unwrap(),
            ],
        ];
        let projected =
            project_relaxation(&mesh, &sections, &model, &modes, budget(), || false).unwrap();
        assert_eq!(projected.len(), 2);
        assert_eq!(projected[0].relaxation_time_s(), 0.002);
        assert_eq!(projected[1].relaxation_time_s(), 0.04);
        for (form, area_weight) in projected.iter().zip([1.5, 0.125]) {
            for i in 0..2 {
                for j in 0..2 {
                    let (ei, ki) = strain[i];
                    let (ej, kj) = strain[j];
                    let expected = area_weight
                        * (bilinear(a, ei, ej)
                            + bilinear(b, ei, kj)
                            + bilinear(b, ki, ej)
                            + bilinear(d, ki, kj));
                    let actual = form.matrix()[2 * i + j];
                    assert!((actual - expected).abs() <= 4e-12 * expected.abs().max(1e-12));
                    assert_eq!(actual, form.matrix()[2 * j + i]);
                }
            }
        }
    }
}

#[test]
fn correlated_weak_material_memory_preserves_rigid_and_drilling_nulls() {
    let (mesh, model) = fixture();
    let diagonal = 1.0 + 1e-14;
    let law = ShellRelaxationBranch::new(
        [1.0, 1.0, 0.0, 1.0, diagonal, 0.0, 0.0, 0.0, 0.0],
        [0.0; 9],
        0.03,
    )
    .unwrap();
    let weak = patch(&mesh, [-1.0, 1.0, 0.0], [0.0; 3]);
    let rigid: Vec<_> = mesh
        .nodes
        .iter()
        .flat_map(|_| [0.2, -0.3, 0.4, 0.0, 0.0, 0.0])
        .collect();
    let mut drilling = vec![0.0; model.free];
    for node in 0..mesh.nodes.len() {
        drilling[6 * node + 5] = node as f64 - 2.0;
    }
    let forms = project_relaxation(
        &mesh,
        &[vec![law.clone()], vec![law]],
        &model,
        &[weak, rigid, drilling],
        budget(),
        || false,
    )
    .unwrap();
    let matrix = forms[0].matrix();
    assert!((matrix[0] / (diagonal - 1.0) - 1.0).abs() < 3e-14);
    assert_eq!(matrix[8], 0.0);
    assert!(matrix[4] < 1e-28);
    assert_eq!(matrix[2], 0.0);
    assert_eq!(matrix[5], 0.0);
}

#[test]
fn hereditary_sources_and_complete_projection_budgets_refuse_without_partial_output() {
    let (mesh, model) = fixture();
    let a = [1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.5];
    for time in [0.0, -1.0, f64::INFINITY, f64::NAN, f64::from_bits(1)] {
        assert!(ShellRelaxationBranch::new(a, [0.0; 9], time).is_err());
    }
    assert!(ShellRelaxationBranch::new(a.map(|x| -x), [0.0; 9], 1.0).is_err());
    assert!(ShellRelaxationBranch::new_coupled(a, a.map(|x| 2.0 * x), a, 1.0).is_err());
    let sections = vec![
        vec![ShellRelaxationBranch::new(a, [0.0; 9], 0.01).unwrap()],
        vec![ShellRelaxationBranch::new(a, [0.0; 9], 0.02).unwrap()],
    ];
    let modes = [patch(&mesh, [0.2, 0.1, 0.3], [0.0; 3])];
    let mut small = budget();
    small.max_element_branches = 1;
    assert!(matches!(
        project_relaxation(&mesh, &sections, &model, &modes, small, || false),
        Err(PlateError::Budget { .. })
    ));
    small = budget();
    small.max_relaxation_times = 1;
    assert!(matches!(
        project_relaxation(&mesh, &sections, &model, &modes, small, || false),
        Err(PlateError::Budget { .. })
    ));
    small = budget();
    small.projection.max_work = 1;
    assert!(project_relaxation(&mesh, &sections, &model, &modes, small, || false).is_err());
    let mut polls = 0;
    assert!(matches!(
        project_relaxation(&mesh, &sections, &model, &modes, budget(), || {
            polls += 1;
            polls > 10
        }),
        Err(PlateError::Cancelled)
    ));
    assert!(project_relaxation(&mesh, &sections, &model, &modes, budget(), || false).is_ok());
    assert!(
        project_relaxation(&mesh, &[vec![], vec![]], &model, &modes, budget(), || false)
            .unwrap()
            .is_empty()
    );
}
