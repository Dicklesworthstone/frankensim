use super::*;
use crate::shell::ShellMesh;
use crate::shell::mitc::{
    ConstraintFrame, NodalDirectors, ShearShellSection, ShellConstraint, ShellMass, assemble_mitc3p,
};
use fs_material::sheet::SheetMatrix;

fn section(index: usize) -> ShearShellSection {
    let scale = 1.0 + 0.4 * index as f64;
    ShearShellSection::new_coupled(
        SheetMatrix::orthotropic(1e6 * scale, 1.4e6 * scale, 0.2, 4e5 * scale).unwrap(),
        [100.0, 20.0, 0.0, 20.0, -80.0, 0.0, 0.0, 0.0, 40.0],
        SheetMatrix::orthotropic(400.0 * scale, 500.0 * scale, 0.2, 160.0 * scale).unwrap(),
        [2e5 * scale, 1e4, 1.6e5 * scale],
        0.85,
        2.0,
        0.006,
        0.00016,
        0.23,
    )
    .unwrap()
}

fn fixture(supported: bool) -> MitcShellModel {
    let mesh = ShellMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [0.4, 0.0, 0.02],
            [0.03, 0.3, 0.01],
            [0.35, 0.27, 0.06],
        ],
        vec![[0, 1, 3], [0, 3, 2]],
    )
    .unwrap();
    let normals = [
        [-0.02, -0.04, 1.0],
        [-0.06, -0.08, 1.0],
        [-0.08, -0.10, 1.0],
        [-0.09, -0.13, 1.0],
    ]
    .map(|n| {
        let length = n.iter().map(|x| x * x).sum::<f64>().sqrt();
        n.map(|x| x / length)
    });
    let directors = NodalDirectors::new(&mesh, normals.to_vec(), 0.8).unwrap();
    let mut constraints = Vec::new();
    if supported {
        for node in 0..3 {
            constraints.push(ShellConstraint {
                node,
                mask: [true; 6],
                frame: ConstraintFrame::Global,
            });
        }
        constraints.push(ShellConstraint {
            node: 3,
            mask: [true, false, false, false, false, false],
            frame: ConstraintFrame::Local([[0.8, 0.0, 0.6], [0.0, 1.0, 0.0], [-0.6, 0.0, 0.8]]),
        });
    }
    assemble_mitc3p(
        &mesh,
        &directors,
        &[section(0), section(1)],
        &constraints,
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: 4,
            max_triangles: 2,
        },
        1e-3,
        || false,
    )
    .unwrap()
}

fn budget() -> ShearShellRelaxationBudget {
    ShearShellRelaxationBudget {
        mesh: MitcBudget {
            max_nodes: 4,
            max_triangles: 2,
        },
        max_modes: 24,
        max_element_branches: 8,
        max_relaxation_times: 4,
        max_projection_terms: 100_000_000,
    }
}

fn spectra() -> Vec<Vec<ShearShellRelaxationBranch>> {
    (0..2)
        .map(|facet| {
            let scale = 1.0 + facet as f64;
            vec![
                ShearShellRelaxationBranch::new_coupled(
                    [4e4, 3e3, 1e3, 3e3, 7e4, 2e3, 1e3, 2e3, 2e4].map(|v| scale * v),
                    [8.0, 1.0, 0.0, 1.0, -12.0, 0.0, 0.0, 0.0, 4.0],
                    [6.0, 0.5, 0.2, 0.5, 3.0, 0.1, 0.2, 0.1, 2.0],
                    [2e4, 3e3, 5e4],
                    0.0003,
                )
                .unwrap(),
                // An independent pure-shear arm cannot be represented by one
                // scalar spectrum proportional to the full relaxed section.
                ShearShellRelaxationBranch::new(
                    [0.0; 9],
                    [0.0; 9],
                    [3e4, -2e3, 1e4],
                    if facet == 0 { 0.002 } else { 0.004 },
                )
                .unwrap(),
            ]
        })
        .collect()
}

fn quadratic(a: &[f64], x: &[f64]) -> f64 {
    let n = x.len();
    let mut value = Dd::ZERO;
    for i in 0..n {
        for j in 0..n {
            value = value + Dd::from_f64(x[i]) * Dd::from_f64(a[n * i + j]) * Dd::from_f64(x[j]);
        }
    }
    value.to_f64()
}

// Independent complete physical work: original equilibrium field is evaluated
// once on the actual displacement, then the ORIGINAL source tensors contract
// its strains directly. No material factor or projected modal form is used.
fn direct_work(
    model: &MitcShellModel,
    source: &[Vec<ShearShellRelaxationBranch>],
    q: &[f64],
    time: f64,
) -> f64 {
    let physical = model.expand(q).unwrap();
    let mut integral = Dd::ZERO;
    for (e, tri) in model.mesh.tris.iter().enumerate() {
        let g = model.mesh.facet(e).unwrap();
        let field = element::PreparedStrainField::new(
            &model.mesh,
            tri,
            g,
            tri.map(|i| model.directors.vectors[i]),
            &model.sections[e],
            &mut || false,
        )
        .unwrap();
        let local = field.expand(&core::array::from_fn(|i| physical[6 * tri[i / 6] + i % 6]));
        for (lambda, weight) in element::stiffness_rule() {
            let (m, k, s) = field.eval(&local, lambda[1], lambda[2]);
            for branch in &source[e] {
                if branch.relaxation_time_s() != time {
                    continue;
                }
                let (a, b, d, shear) = branch.coupled_matrices();
                let mut power = Dd::ZERO;
                for i in 0..3 {
                    for j in 0..3 {
                        power = power
                            + Dd::from_f64(m[i]) * Dd::from_f64(a[3 * i + j]) * Dd::from_f64(m[j])
                            + Dd::from_f64(m[i]) * Dd::from_f64(b[3 * i + j]) * Dd::from_f64(k[j])
                            + Dd::from_f64(k[i]) * Dd::from_f64(b[3 * i + j]) * Dd::from_f64(m[j])
                            + Dd::from_f64(k[i]) * Dd::from_f64(d[3 * i + j]) * Dd::from_f64(k[j]);
                    }
                }
                for i in 0..2 {
                    for j in 0..2 {
                        power = power
                            + Dd::from_f64(s[i]) * Dd::from_f64(shear[i][j]) * Dd::from_f64(s[j]);
                    }
                }
                integral = integral + power * Dd::from_f64(g.area_m2 * weight);
            }
        }
    }
    integral.to_f64()
}

#[test]
fn original_oblique_support_and_bubble_space_preserve_every_tensor_time() {
    let model = fixture(true);
    let n = model.pencil.free;
    assert_eq!(n, 5);
    let modes: Vec<_> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let source = spectra();
    let projected = project_relaxation(&model, &source, &modes, budget(), || false).unwrap();
    assert_eq!(projected.len(), 3);
    assert_eq!(
        projected
            .iter()
            .map(|b| b.relaxation_time_s())
            .collect::<Vec<_>>(),
        [0.0003, 0.002, 0.004]
    );
    for q in [[0.2, -0.1, 0.3, -0.4, 0.07], [-0.1, 0.4, 0.2, 0.06, -0.2]] {
        for branch in &projected {
            let expected = direct_work(&model, &source, &q, branch.relaxation_time_s());
            let actual = quadratic(branch.matrix(), &q);
            assert!(expected > 0.0);
            assert!(
                (actual - expected).abs() < 3e-11 * expected,
                "{actual:e} != {expected:e}"
            );
            assert!(
                branch
                    .matrix()
                    .iter()
                    .enumerate()
                    .any(|(i, x)| i / n != i % n && x.abs() > 1e-4)
            );
        }
    }
    // The expansion really differs from inserting reduced numbers into bare
    // physical coordinate slots; this catches accidental loss of support frames.
    let physical = model.expand(&[1.0; 5]).unwrap();
    let bare: Vec<_> = model
        .pencil
        .dof_map
        .iter()
        .map(|dof| dof.map_or(0.0, |_| 1.0))
        .collect();
    assert!(physical.iter().zip(bare).any(|(a, b)| (a - b).abs() > 0.1));
}

#[test]
fn rigid_motion_and_drilling_remain_null_while_shear_carries_memory() {
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.0, 0.3, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 1.0).unwrap();
    let model = assemble_mitc3p(
        &mesh,
        &directors,
        &[section(0)],
        &[],
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: 3,
            max_triangles: 1,
        },
        1e-3,
        || false,
    )
    .unwrap();
    let mut modes = Vec::new();
    for axis in 0..6 {
        let mut v = [0.0; 3];
        v[axis % 3] = 1.0;
        modes.push(
            mesh.nodes
                .iter()
                .flat_map(|p| {
                    let u = if axis < 3 {
                        v
                    } else {
                        [
                            v[1] * p[2] - v[2] * p[1],
                            v[2] * p[0] - v[0] * p[2],
                            v[0] * p[1] - v[1] * p[0],
                        ]
                    };
                    let w = if axis < 3 { [0.0; 3] } else { v };
                    [u[0], u[1], u[2], w[0], w[1], w[2]]
                })
                .collect(),
        );
    }
    for node in 0..3 {
        let mut drill = vec![0.0; 18];
        drill[6 * node + 5] = 1.0;
        modes.push(drill);
    }
    modes.push(
        mesh.nodes
            .iter()
            .flat_map(|p| [0.0, 0.0, p[0] + 0.3 * p[1], 0.0, 0.0, 0.0])
            .collect(),
    );
    let source = vec![spectra().remove(0)];
    let projected = project_relaxation(&model, &source, &modes, budget(), || false).unwrap();
    for branch in projected {
        let n = modes.len();
        for i in 0..9 {
            assert!(
                branch.matrix()[n * i + i] <= 1e-20,
                "physical null {i}: {:e}",
                branch.matrix()[n * i + i]
            );
        }
        assert!(
            branch.matrix()[n * 9 + 9] > 1.0,
            "transverse shear was omitted"
        );
    }
    let mut force = vec![0.0; 18];
    model.pencil.k.spmv(&modes[6], &mut force);
    assert!(
        force[5] > 0.0,
        "fixture must retain the elastic numerical drilling gauge"
    );
}

#[test]
fn original_nonbinary_weak_shear_channel_survives_passive_admission() {
    for scale in [1e-120, 1.0, 1e120] {
        let a = 3.7 * scale;
        let b = a + a * 2.0_f64.powi(-48);
        let branch = ShearShellRelaxationBranch::new([0.0; 9], [0.0; 9], [a, a, b], 0.003).unwrap();
        let channel = branch.factor.channels([0.0; 3], [0.0; 3], [-1.0, 1.0]);
        let actual = channel.iter().map(|x| x * x).sum::<f64>();
        assert!(
            (actual / (b - a) - 1.0).abs() < 2e-12,
            "{actual:e} != {:e}",
            b - a
        );
    }
}

#[test]
fn bounds_zero_spectra_and_invalid_materials_refuse_before_projection() {
    let model = fixture(true);
    let modes = vec![vec![0.2; model.pencil.free]];
    let empty = vec![Vec::new(); 2];
    assert!(
        project_relaxation(&model, &empty, &modes, budget(), || false)
            .unwrap()
            .is_empty()
    );
    let source = spectra();
    for limit in [
        ShearShellRelaxationBudget {
            max_element_branches: 3,
            ..budget()
        },
        ShearShellRelaxationBudget {
            max_relaxation_times: 2,
            ..budget()
        },
        ShearShellRelaxationBudget {
            max_projection_terms: 0,
            ..budget()
        },
    ] {
        assert!(project_relaxation(&model, &source, &modes, limit, || false).is_err());
    }
    assert!(project_relaxation(&model, &source[..1], &modes, budget(), || false).is_err());
    assert!(project_relaxation(&model, &source, &[vec![0.0; 1]], budget(), || false).is_err());
    assert!(
        project_relaxation(
            &model,
            &source,
            &[vec![f64::NAN; model.pencil.free]],
            budget(),
            || false
        )
        .is_err()
    );
    assert!(matches!(
        project_relaxation(&model, &source, &modes, budget(), || true),
        Err(PlateError::Cancelled)
    ));
    for time in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-320] {
        assert!(
            ShearShellRelaxationBranch::new([0.0; 9], [0.0; 9], [1.0, 0.0, 1.0], time).is_err()
        );
    }
    assert!(ShearShellRelaxationBranch::new([0.0; 9], [0.0; 9], [1.0, 2.0, 1.0], 0.1).is_err());
    assert!(projection_terms(usize::MAX, 1, 1, 1, 1).is_none());
}
