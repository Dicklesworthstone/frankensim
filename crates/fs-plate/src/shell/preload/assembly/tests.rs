use super::*;
use crate::PlateSection;
use crate::shell::junction::{RotationCoupling, TieOffset, TieSide};
use crate::shell::preload::{equilibrate_shell_materials, tests::reference_newton};
use fs_material::sheet::SheetMatrix;

struct Part {
    mesh: ShellMesh,
    sections: Vec<ShellSection>,
    fixed: Vec<usize>,
    loads: Vec<f64>,
    thickness: Vec<f64>,
}

impl Part {
    fn rectangle(first_column: usize, nx: usize, ny: usize, z: f64, supported: bool) -> Self {
        let mut nodes = Vec::new();
        let mut tris = Vec::new();
        let mut fixed = Vec::new();
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push([
                    (first_column + i) as f64 * 0.1,
                    j as f64 * 0.2 / ny as f64,
                    z,
                ]);
                if supported && i == 0 {
                    fixed.push(j * (nx + 1));
                }
            }
        }
        for j in 0..ny {
            for i in 0..nx {
                let a = j * (nx + 1) + i;
                tris.push([a, a + 1, a + nx + 2]);
                tris.push([a, a + nx + 2, a + nx + 1]);
            }
        }
        let mesh = ShellMesh::new(nodes, tris).unwrap();
        let plate = PlateSection::isotropic(9e9, 0.27, 0.004, 430.).unwrap();
        let a = plate.d.map(|d| (12. / plate.thickness.powi(2)) * d);
        let section = ShellSection::new(
            plate,
            SheetMatrix::orthotropic(a[0], a[4], a[1], a[8]).unwrap(),
        )
        .unwrap();
        Self {
            sections: vec![section; mesh.tris.len()],
            loads: vec![0.; 6 * mesh.nodes.len()],
            thickness: vec![0.004; mesh.nodes.len()],
            mesh,
            fixed,
        }
    }

    fn input(&self) -> ShellPreloadComponent<'_> {
        ShellPreloadComponent {
            mesh: &self.mesh,
            sections: &self.sections,
            boundary: &self.fixed,
            support: if self.fixed.is_empty() {
                ShellSupport::Free
            } else {
                ShellSupport::Clamped
            },
            beams: &[],
            loads: &self.loads,
            nodal_thickness_m: &self.thickness,
            reference_state: ReferenceLoadState::Unloaded,
        }
    }
}

fn budget() -> TieBudget {
    TieBudget {
        max_segments: 1000,
        max_slave_dofs: 1000,
        gap_tolerance_m: 1e-10,
        angle_tolerance_rad: 1e-8,
    }
}

fn seam(slave_ny: usize, offset: TieOffset) -> EdgeTie {
    EdgeTie {
        a: TieSide {
            component: 0,
            nodes: vec![2, 5, 8],
        },
        b: TieSide {
            component: 1,
            nodes: (0..=slave_ny).map(|j| 3 * j).collect(),
        },
        offset,
        rotations: RotationCoupling::Full,
        kind: TieKind::Rigid,
    }
}

fn options() -> PreloadOptions {
    PreloadOptions {
        relative_force_tolerance: 1e-10,
        ..PreloadOptions::default()
    }
}

#[test]
fn one_part_reuses_the_single_shell_equilibrium_and_tangent_bitwise() {
    let mut part = Part::rectangle(0, 4, 2, 0., true);
    part.loads[6 * 9 + 2] = -0.5;
    let original = part.mesh.clone();
    let direct = equilibrate_shell_materials(
        &part.mesh,
        &part.sections,
        &part.fixed,
        ShellSupport::Clamped,
        &[],
        &part.loads,
        &[],
        options(),
        &mut reference_newton,
    )
    .unwrap();
    let assembled = equilibrate_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert_eq!(part.mesh, original);
    assert_eq!(assembled.displacements[0], direct.displacement);
    assert_eq!(assembled.residual_force_n, direct.residual_force_n);
    assert_eq!(assembled.stored_energy_j, direct.stored_energy_j);
    assert_eq!(assembled.norm_length_m, direct.norm_length_m);
    assert_eq!(assembled.iterations, direct.iterations);
    for row in 0..direct.model.free {
        assert_eq!(assembled.assembly.model.k.row(row), direct.model.k.row(row));
        assert_eq!(assembled.assembly.model.m.row(row), direct.model.m.row(row));
    }
}

#[test]
fn tied_nonlinear_cantilever_matches_monolithic_body_before_modes() {
    let mut whole = Part::rectangle(0, 4, 2, 0., true);
    let left = Part::rectangle(0, 2, 2, 0., true);
    let mut right = Part::rectangle(2, 2, 2, 0., false);
    // Both a slave-seam force and a force on the unsupported second part must
    // enter the single equilibrium, rather than vanish with slave elimination.
    right.loads[6 * 3 + 2] = -0.25;
    right.loads[6 * 5 + 2] = -0.5;
    whole.loads[6 * 7 + 2] = -0.25;
    whole.loads[6 * 9 + 2] = -0.5;
    let direct = equilibrate_shell_materials(
        &whole.mesh,
        &whole.sections,
        &whole.fixed,
        ShellSupport::Clamped,
        &[],
        &whole.loads,
        &[],
        options(),
        &mut reference_newton,
    )
    .unwrap();
    let tied = equilibrate_shell_assembly(
        &[left.input(), right.input()],
        &[seam(2, TieOffset::Midsurface)],
        &[],
        budget(),
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap();
    assert!(tied.displacements[1][6 * 5 + 2] < 0.);
    assert!(tied.residual_force_n < 1e-8);
    assert!(tied.stored_energy_j > 0.);
    let n = direct.model.free;
    assert_eq!(tied.assembly.model.free, n);
    let mut permutation = vec![usize::MAX; n];
    for component in 0..2 {
        for local in 0..9 {
            let whole_node = 5 * (local / 3) + local % 3 + 2 * component;
            for dof in 0..6 {
                let a = tied.displacements[component][6 * local + dof];
                let b = direct.displacement[6 * whole_node + dof];
                assert!(
                    (a - b).abs() < 2e-10,
                    "component={component} node={local} dof={dof}: {a} vs {b}"
                );
                if let Some(retained) = tied.assembly.map[component][local][dof] {
                    permutation[retained] = direct.model.dof_map[6 * whole_node + dof].unwrap();
                }
            }
        }
    }
    assert!(permutation.iter().all(|&i| i < n));
    for i in 0..n {
        for j in 0..n {
            let a = tied.assembly.model.k.get(i, j);
            let b = direct.model.k.get(permutation[i], permutation[j]);
            assert!(
                (a - b).abs() < 2e-7 + 1e-9 * b.abs(),
                "tangent ({i},{j}): {a} vs {b}"
            );
            let a = tied.assembly.model.m.get(i, j);
            let b = direct.model.m.get(permutation[i], permutation[j]);
            assert!((a - b).abs() < 1e-13, "inertia ({i},{j}): {a} vs {b}");
        }
    }
    assert!((tied.stored_energy_j - direct.stored_energy_j).abs() < 1e-11);
    // Equal mass and physical equilibrium tangents imply the same complete
    // small-vibration pencil, before any frequency-window/modal truncation.
    let unloaded = assemble_stiffened_shell_materials(
        &whole.mesh,
        &whole.sections,
        &whole.fixed,
        ShellSupport::Clamped,
        &[],
    )
    .unwrap();
    let z = direct.model.dof_map[6 * 9 + 2].unwrap();
    assert!((direct.model.k.get(z, z) - unloaded.k.get(z, z)).abs() > 1e-3);
}

#[test]
fn nonmatching_eccentric_tie_preserves_energy_force_and_tangent_derivatives() {
    let left = Part::rectangle(0, 2, 2, 0., true);
    let right = Part::rectangle(2, 2, 3, 0.002, false);
    let mut checked = false;
    let solved = equilibrate_shell_assembly(
        &[left.input(), right.input()],
        &[seam(3, TieOffset::Skin { side: 1 })],
        &[],
        budget(),
        options(),
        &mut |p, x, limit| {
            if !checked {
                checked = true;
                let x: Vec<_> = (0..p.dimension())
                    .map(|i| ((i * 17 % 19) as f64 - 9.) * 1e-6)
                    .collect();
                let direction: Vec<_> = (0..p.dimension())
                    .map(|i| ((i * 11 % 13) as f64 - 6.) / 6.)
                    .collect();
                let eps = 1e-8;
                let plus: Vec<_> = x.iter().zip(&direction).map(|(x, d)| x + eps * d).collect();
                let minus: Vec<_> = x.iter().zip(&direction).map(|(x, d)| x - eps * d).collect();
                let u = p.physical(&x);
                let mut force = vec![0.; p.dimension()];
                p.force_energy(&u, &mut force).unwrap();
                let analytical: f64 = force
                    .iter()
                    .zip(p.physical(&direction))
                    .map(|(f, v)| f * v)
                    .sum();
                let energy = |x: &[f64]| {
                    p.force_energy(&p.physical(x), &mut vec![0.; p.dimension()])
                        .unwrap()
                };
                let difference = (energy(&plus) - energy(&minus)) / (2. * eps);
                assert!(
                    (difference - analytical).abs() < 1e-6 * analytical.abs().max(1.),
                    "T^T force must differentiate component energy: {difference} vs {analytical}"
                );
                let mut a = vec![0.; p.dimension()];
                let mut b = a.clone();
                let mut action = a.clone();
                p.residual(&plus, &mut a);
                p.residual(&minus, &mut b);
                p.jacobian_apply(&x, &direction, &mut action);
                let tangent = p.tangent(&u).unwrap();
                let mut physical_action = vec![0.; p.dimension()];
                tangent.spmv(&p.physical(&direction), &mut physical_action);
                let factored_action = p.base.solve(&physical_action);
                for i in 0..p.dimension() {
                    let (cols, values) = tangent.row(i);
                    for (&j, &value) in cols.iter().zip(values) {
                        assert_eq!(value, tangent.get(j, i), "tied tangent symmetry");
                    }
                    let assembled = factored_action[i] * p.scale[i];
                    assert!(
                        (assembled - action[i]).abs() < 1e-7 * action[i].abs().max(1.),
                        "assembled tangent must agree with matrix-free action"
                    );
                }
                for ((a, b), exact) in a.iter().zip(b).zip(action) {
                    let fd = (a - b) / (2. * eps);
                    assert!(
                        (fd - exact).abs() < 2e-5 * exact.abs().max(1.),
                        "T^T K T must differentiate the original residual: {fd} vs {exact}"
                    );
                }
            }
            reference_newton(p, x, limit)
        },
        || false,
    )
    .unwrap();
    assert!(checked);
    assert!(solved.displacements.iter().flatten().all(|v| *v == 0.));
    assert!(
        solved
            .assembly
            .slave_couplings()
            .iter()
            .any(|s| s.master_terms.len() > 1)
    );
}

#[test]
fn loaded_references_missing_supports_and_false_convergence_return_no_body() {
    let mut part = Part::rectangle(0, 2, 2, 0., true);
    part.loads[6 * 5 + 2] = -0.5;
    let mut loaded = part.input();
    loaded.reference_state = ReferenceLoadState::Loaded {
        description: "survey under string tension".into(),
    };
    let error = equilibrate_shell_assembly(
        &[loaded],
        &[],
        &[],
        budget(),
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap_err();
    assert!(error.contains("unloaded reference"));
    let mut free = part.input();
    free.boundary = &[];
    free.support = ShellSupport::Free;
    let error = equilibrate_shell_assembly(
        &[free],
        &[],
        &[],
        budget(),
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap_err();
    assert!(error.contains("rigid modes"));
    let error = equilibrate_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        options(),
        &mut |_, x, _| Ok((x, 0)),
        || false,
    )
    .unwrap_err();
    assert!(error.contains("physical force balance"));
    let slave = Part::rectangle(2, 2, 2, 0., true);
    let error = equilibrate_shell_assembly(
        &[part.input(), slave.input()],
        &[seam(2, TieOffset::Midsurface)],
        &[],
        budget(),
        options(),
        &mut reference_newton,
        || false,
    )
    .unwrap_err();
    assert!(error.contains("slave supports"));
    let error = equilibrate_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        options(),
        &mut reference_newton,
        || true,
    )
    .unwrap_err();
    assert!(error.contains("cancelled"));
    let mut polls = 0;
    equilibrate_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        options(),
        &mut reference_newton,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 1);
    let mut repeated_polls = 0;
    let error = equilibrate_shell_assembly(
        &[part.input()],
        &[],
        &[],
        budget(),
        options(),
        &mut reference_newton,
        || {
            repeated_polls += 1;
            repeated_polls == polls
        },
    )
    .unwrap_err();
    assert_eq!(repeated_polls, polls);
    assert!(error.contains("cancelled"));
}
