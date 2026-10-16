use super::*;

fn section(h: f64) -> ShearShellSection {
    ShearShellSection::from_orthotropic_solid(
        10.0,
        10.0,
        0.3,
        10.0 / 2.6,
        10.0 / 2.6,
        10.0 / 2.6,
        h,
        2.0,
        0.0,
        5.0 / 6.0,
    )
    .unwrap()
}
fn patch(n: usize, curved: bool, distorted: bool) -> (ShellMesh, NodalDirectors) {
    let mut nodes = Vec::new();
    let mut directors = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            let mut x = i as f64 / n as f64;
            let mut y = j as f64 / n as f64;
            if distorted && i > 0 && i < n && j > 0 && j < n {
                x += 0.17 / n as f64 * ((i * 13 + j * 7) % 11) as f64 / 11.0;
                y -= 0.14 / n as f64 * ((i * 5 + j * 17) % 13) as f64 / 13.0;
            }
            let z = if curved { 0.1 * (x * x + y * y) } else { 0.0 };
            nodes.push([x, y, z]);
            let normal = if curved {
                [-0.2 * x, -0.2 * y, 1.0]
            } else {
                [0.0, 0.0, 1.0]
            };
            let length = det::sqrt(dot(normal, normal));
            directors.push(normal.map(|x| x / length));
        }
    }
    let mut tris = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let a = j * (n + 1) + i;
            let b = a + 1;
            let c = a + n + 1;
            let d = c + 1;
            tris.extend([[a, b, d], [a, d, c]]);
        }
    }
    let mesh = ShellMesh::new(nodes, tris).unwrap();
    let directors = NodalDirectors::new(&mesh, directors, 0.5).unwrap();
    (mesh, directors)
}
fn assemble(
    mesh: &ShellMesh,
    directors: &NodalDirectors,
    h: f64,
    constraints: &[ShellConstraint],
    alpha: f64,
) -> MitcShellModel {
    assemble_mitc3p(
        mesh,
        directors,
        &vec![section(h); mesh.tris.len()],
        constraints,
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: mesh.nodes.len(),
            max_triangles: mesh.tris.len(),
        },
        alpha,
        || false,
    )
    .unwrap()
}
fn energy(model: &MitcShellModel, u: &[f64]) -> f64 {
    let mut ku = vec![0.0; u.len()];
    model.pencil.k.spmv(u, &mut ku);
    u.iter().zip(ku).map(|(a, b)| a * b).sum()
}
fn compare(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g0_rigid_body_modes_exact() {
    for (mesh, directors) in [patch(1, false, false), patch(3, true, false)] {
        let model = assemble(&mesh, &directors, 0.2, &[], 1e-3);
        for axis in 0..6 {
            let mut u = vec![0.0; 6 * mesh.nodes.len()];
            for (i, p) in mesh.nodes.iter().enumerate() {
                if axis < 3 {
                    u[6 * i + axis] = 1.0;
                } else {
                    let mut omega = [0.0; 3];
                    omega[axis - 3] = 1.0;
                    u[6 * i..6 * i + 3].copy_from_slice(&cross(omega, *p));
                    u[6 * i + 3..6 * i + 6].copy_from_slice(&omega);
                }
            }
            let mut residual = vec![0.0; u.len()];
            model.pencil.k.spmv(&u, &mut residual);
            let scale = (0..u.len())
                .map(|i| model.pencil.k.row(i).1.iter().map(|v| v.abs()).sum::<f64>())
                .fold(0.0, f64::max)
                * u.iter().map(|x| x.abs()).fold(0.0, f64::max);
            assert!(
                residual.iter().all(|r| r.abs() <= 1e-12 * scale.max(1e-30)),
                "rigid axis {axis}"
            );
        }
        let modes = fs_modal::eigh_gen_dense(
            &model.pencil.k.to_dense(),
            &model.pencil.m.to_dense(),
            model.pencil.free,
        )
        .unwrap();
        let maximum = modes.last().unwrap().lambda;
        assert_eq!(
            modes
                .iter()
                .filter(|m| m.lambda.abs() <= 1e-10 * maximum)
                .count(),
            6
        );
    }
    // The single, unsupported enriched triangle must also have only six.
    let mesh = ShellMesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.2, 0.8, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.5).unwrap();
    let model = assemble(&mesh, &directors, 0.2, &[], 1e-3);
    let modes =
        fs_modal::eigh_gen_dense(&model.pencil.k.to_dense(), &model.pencil.m.to_dense(), 18)
            .unwrap();
    let maximum = modes.last().unwrap().lambda;
    assert_eq!(
        modes
            .iter()
            .filter(|m| m.lambda.abs() <= 1e-10 * maximum)
            .count(),
        6
    );
}

#[test]
fn g0_membrane_patch_constant_strain() {
    let (mesh, directors) = patch(3, false, true);
    let model = assemble(&mesh, &directors, 0.1, &[], 1e-3);
    let strain = [0.07, -0.03, 0.04];
    let mut u = vec![0.0; 6 * mesh.nodes.len()];
    for (i, p) in mesh.nodes.iter().enumerate() {
        u[6 * i] = strain[0] * p[0] + 0.5 * strain[2] * p[1];
        u[6 * i + 1] = 0.5 * strain[2] * p[0] + strain[1] * p[1];
    }
    let a = section(0.1).matrices().0;
    let expected = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| strain[i] * a[3 * i + j] * strain[j])
                .sum::<f64>()
        })
        .sum();
    compare(energy(&model, &u), expected, 1e-12);
}

#[test]
fn g0_bending_patch_constant_curvature() {
    let (mesh, directors) = patch(3, false, true);
    let model = assemble(&mesh, &directors, 0.1, &[], 1e-3);
    let curvature = [0.7, -0.3, 0.2];
    let mut u = vec![0.0; 6 * mesh.nodes.len()];
    for (i, &[x, y, _]) in mesh.nodes.iter().enumerate() {
        u[6 * i + 2] =
            0.5 * curvature[0] * x * x + 0.5 * curvature[1] * y * y + 0.5 * curvature[2] * x * y;
        u[6 * i + 3] = curvature[1] * y + 0.5 * curvature[2] * x;
        u[6 * i + 4] = -curvature[0] * x - 0.5 * curvature[2] * y;
    }
    let d = section(0.1).matrices().1;
    let expected = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| curvature[i] * d[3 * i + j] * curvature[j])
                .sum::<f64>()
        })
        .sum();
    compare(energy(&model, &u), expected, 1e-12);
}

fn supported_square(n: usize, h: f64, alpha: f64) -> (ShellMesh, MitcShellModel, Vec<f64>) {
    let (mesh, directors) = patch(n, false, false);
    let constraints: Vec<_> = mesh
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 0.0 || p[1] == 0.0 || p[0] == 1.0 || p[1] == 1.0)
        .map(|(node, _)| ShellConstraint {
            node,
            mask: [true, true, true, false, false, false],
            frame: ConstraintFrame::Global,
        })
        .collect();
    let model = assemble(&mesh, &directors, h, &constraints, alpha);
    let d = section(h).matrices().1[0];
    let mut load = vec![0.0; 6 * mesh.nodes.len()];
    for (e, t) in mesh.tris.iter().enumerate() {
        let area = mesh.facet(e).unwrap().area_m2;
        for &i in t {
            load[6 * i + 2] += d * area / 3.0;
        }
    }
    (mesh, model, load)
}

#[test]
fn g1_thin_limit_no_shear_locking() {
    let mut values = Vec::new();
    for n in [8, 16, 32] {
        let (mesh, model, load) = supported_square(n, 1e-4, 1e-3);
        let u = solve_static(&model, &load).unwrap();
        let node = mesh
            .nodes
            .iter()
            .position(|p| p[0] == 0.5 && p[1] == 0.5)
            .unwrap();
        values.push(u[6 * node + 2]);
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1a\",\"test\":\"g1_thin_limit_no_shear_locking\",\"tier\":\"G1\",\"n\":{n},\"thickness_ratio\":0.0001,\"observed\":{},\"reference\":0.00406,\"tolerance\":0.01}}",
            u[6 * node + 2]
        );
    }
    assert!((values.last().unwrap() / 0.00406 - 1.0).abs() <= 0.01);
    assert!(
        values
            .windows(2)
            .all(|w| (w[1] - 0.00406).abs() < (w[0] - 0.00406).abs())
    );
}

#[test]
fn g0_drilling_alpha_insensitivity() {
    let mut values = Vec::new();
    for alpha in [1e-4, 1e-3, 1e-2] {
        let (mesh, model, load) = supported_square(4, 0.03, alpha);
        let u = solve_static(&model, &load).unwrap();
        let i = mesh
            .nodes
            .iter()
            .position(|p| p[0] == 0.5 && p[1] == 0.5)
            .unwrap();
        values.push(u[6 * i + 2]);
    }
    let reference = values[1];
    for value in values {
        assert!((value / reference - 1.0).abs() < 1e-4);
    }
}

#[test]
fn g0_consistent_mass_preserves_rigid_mass() {
    let (mesh, directors) = patch(3, true, false);
    let model = assemble(&mesh, &directors, 0.1, &[], 1e-3);
    let expected = (0..mesh.tris.len())
        .map(|e| mesh.facet(e).unwrap().area_m2 * section(0.1).areal_mass_kg_m2())
        .sum::<f64>();
    for c in 0..3 {
        let mut one = vec![0.0; model.pencil.free];
        for n in 0..mesh.nodes.len() {
            one[6 * n + c] = 1.0;
        }
        let mut mone = vec![0.0; one.len()];
        model.pencil.m.spmv(&one, &mut mone);
        compare(
            one.iter().zip(mone).map(|(a, b)| a * b).sum(),
            expected,
            1e-14,
        );
    }
}

#[test]
fn g3_node_renumbering_invariance() {
    let (mesh, directors) = patch(2, true, true);
    let n = mesh.nodes.len();
    let mut nodes = mesh.nodes.clone();
    nodes.reverse();
    let mut normals = directors.vectors.clone();
    normals.reverse();
    let renamed = ShellMesh::new(
        nodes,
        mesh.tris.iter().map(|t| t.map(|i| n - 1 - i)).collect(),
    )
    .unwrap();
    let renamed_dirs = NodalDirectors::new(&renamed, normals, 0.5).unwrap();
    let a = assemble(&mesh, &directors, 0.2, &[], 1e-3);
    let b = assemble(&renamed, &renamed_dirs, 0.2, &[], 1e-3);
    let eigen = |a: &MitcShellModel| {
        fs_modal::eigh_gen_dense(
            &a.pencil.k.to_dense(),
            &a.pencil.m.to_dense(),
            a.pencil.free,
        )
        .unwrap()
    };
    let ea = eigen(&a);
    let eb = eigen(&b);
    let scale = ea.last().unwrap().lambda;
    for (a, b) in ea.iter().zip(eb) {
        assert!((a.lambda - b.lambda).abs() < 1e-12 * scale);
    }
}

#[test]
fn g0_director_refusals() {
    let (mesh, _) = patch(1, false, false);
    let normals = vec![[0.0, 0.0, 1.0]; mesh.nodes.len()];
    for bad in [
        vec![],
        vec![[0.0, 0.0, 2.0]; normals.len()],
        vec![[f64::NAN, 0.0, 1.0]; normals.len()],
        vec![[0.0, 0.0, -1.0]; normals.len()],
    ] {
        assert!(matches!(
            NodalDirectors::new(&mesh, bad, 0.5),
            Err(PlateError::BadDirector { .. })
        ));
    }
    assert!(NodalDirectors::new(&mesh, normals, f64::NAN).is_err());
    let (curved, dirs) = patch(1, true, false);
    assert!(NodalDirectors::new(&curved, dirs.vectors, 1.0).is_err());
}

#[test]
fn g0_constraint_frame_refusals() {
    let (mesh, dirs) = patch(1, false, false);
    let sections = vec![section(0.1); mesh.tris.len()];
    let budget = MitcBudget {
        max_nodes: 4,
        max_triangles: 2,
    };
    let build = |constraints: &[ShellConstraint]| {
        assemble_mitc3p(
            &mesh,
            &dirs,
            &sections,
            constraints,
            ShellMass::Consistent,
            budget,
            1e-3,
            || false,
        )
    };
    let support = ShellConstraint {
        node: 0,
        mask: [true; 6],
        frame: ConstraintFrame::Global,
    };
    for constraints in [
        vec![support, support],
        vec![ShellConstraint { node: 4, ..support }],
        vec![ShellConstraint {
            frame: ConstraintFrame::Local([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]]),
            ..support
        }],
    ] {
        assert!(matches!(
            build(&constraints),
            Err(PlateError::BadConstraint { .. })
        ));
    }
    let rotation = [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]];
    let model = build(&[ShellConstraint {
        frame: ConstraintFrame::Local(rotation),
        mask: [true, false, false, false, false, false],
        ..support
    }])
    .unwrap();
    let reduced: Vec<_> = (0..model.pencil.free)
        .map(|i| (i as f64 + 1.0) * 0.01)
        .collect();
    let full = model.expand(&reduced).unwrap();
    let force: Vec<_> = (0..full.len()).map(|i| (i as f64 - 7.0) * 0.03).collect();
    let pulled = model.project_load(&force).unwrap();
    compare(
        full.iter().zip(&force).map(|(a, b)| a * b).sum(),
        reduced.iter().zip(&pulled).map(|(a, b)| a * b).sum(),
        1e-14,
    );
    assert!(model.expand(&[]).is_err());
    assert!(model.project_load(&[]).is_err());
    assert!(solve_static(&assemble(&mesh, &dirs, 0.1, &[], 1e-3), &vec![0.0; 24]).is_err());
}

#[test]
fn g0_budget_cap_and_cap_plus_one() {
    let (mesh, dirs) = patch(1, false, false);
    let sections = vec![section(0.1); mesh.tris.len()];
    let budget = MitcBudget {
        max_nodes: mesh.nodes.len(),
        max_triangles: mesh.tris.len(),
    };
    let run = |budget| {
        assemble_mitc3p(
            &mesh,
            &dirs,
            &sections,
            &[],
            ShellMass::Consistent,
            budget,
            1e-3,
            || false,
        )
    };
    run(budget).unwrap();
    for budget in [
        MitcBudget {
            max_nodes: budget.max_nodes - 1,
            ..budget
        },
        MitcBudget {
            max_triangles: budget.max_triangles - 1,
            ..budget
        },
    ] {
        assert!(
            matches!(run(budget),Err(PlateError::Budget{requested,limit,..}) if requested==limit+1)
        );
    }
    assert!(matches!(
        assemble_mitc3p(
            &mesh,
            &dirs,
            &sections,
            &[],
            ShellMass::Consistent,
            budget,
            1e-3,
            || true
        ),
        Err(PlateError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        assemble_mitc3p(
            &mesh,
            &dirs,
            &sections,
            &[],
            ShellMass::Consistent,
            budget,
            1e-3,
            || {
                polls += 1;
                polls > 4
            }
        ),
        Err(PlateError::Cancelled)
    ));
}

#[test]
fn g0_section_refusals_and_explicit_mass_choice() {
    let construct = |g13, kappa| {
        ShearShellSection::from_orthotropic_solid(
            10.0, 5.0, 0.2, 2.0, g13, 1.0, 0.1, 2.0, 0.0, kappa,
        )
    };
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(construct(bad, 5.0 / 6.0).is_err());
        assert!(construct(1.0, bad).is_err());
    }
    let (mesh, dirs) = patch(1, false, false);
    let sections = vec![construct(1.0, 5.0 / 6.0).unwrap(); mesh.tris.len()];
    for mass in [ShellMass::Consistent, ShellMass::Lumped] {
        let model = assemble_mitc3p(
            &mesh,
            &dirs,
            &sections,
            &[],
            mass,
            MitcBudget {
                max_nodes: 4,
                max_triangles: 2,
            },
            1e-3,
            || false,
        )
        .unwrap();
        assert_eq!(model.mass_model(), mass);
        let mut one = vec![0.0; 24];
        for node in 0..4 {
            one[6 * node] = 1.0;
        }
        let mut mone = vec![0.0; 24];
        model.pencil.m.spmv(&one, &mut mone);
        compare(one.iter().zip(mone).map(|(a, b)| a * b).sum(), 0.2, 1e-14);
    }
}
