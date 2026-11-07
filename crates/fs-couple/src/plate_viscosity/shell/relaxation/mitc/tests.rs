use super::super::tests::{phs, rk4};
use super::*;
use fs_material::sheet::SheetMatrix;
use fs_modal::ModePair;
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::strain::{PreparedElementStrain, integration_rule};
use fs_plate::shell::mitc::{
    ConstraintFrame, NodalDirectors, ShearShellSection, ShellConstraint, ShellMass, assemble_mitc3p,
};

struct Fixture {
    mesh: ShellMesh,
    directors: NodalDirectors,
    sections: Vec<ShearShellSection>,
    model: MitcShellModel,
    pairs: Vec<ModePair>,
    source: Vec<Vec<ShearShellRelaxationBranch>>,
}

fn fixture(mass: ShellMass) -> Fixture {
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
    .map(|v| {
        let length = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        v.map(|x| x / length)
    });
    let directors = NodalDirectors::new(&mesh, normals.to_vec(), 0.8).unwrap();
    let sections: Vec<_> = (0..2)
        .map(|i| {
            let scale = 1.0 + 0.4 * i as f64;
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
        })
        .collect();
    let mut supports: Vec<_> = (0..3)
        .map(|node| ShellConstraint {
            node,
            mask: [true; 6],
            frame: ConstraintFrame::Global,
        })
        .collect();
    supports.push(ShellConstraint {
        node: 3,
        mask: [true, false, false, false, false, false],
        frame: ConstraintFrame::Local([[0.8, 0.0, 0.6], [0.0, 1.0, 0.0], [-0.6, 0.0, 0.8]]),
    });
    let model = assemble_mitc3p(
        &mesh,
        &directors,
        &sections,
        &supports,
        mass,
        MitcBudget {
            max_nodes: 4,
            max_triangles: 2,
        },
        1e-3,
        || false,
    )
    .unwrap();
    let n = model.pencil().free;
    let k: Vec<_> = (0..n)
        .flat_map(|i| (0..n).map(move |j| (i, j)))
        .map(|(i, j)| model.pencil().k.get(i, j))
        .collect();
    let m: Vec<_> = (0..n)
        .flat_map(|i| (0..n).map(move |j| (i, j)))
        .map(|(i, j)| model.pencil().m.get(i, j))
        .collect();
    let pairs = fs_modal::eigh_gen_dense(&k, &m, n).unwrap();
    let source = (0..2)
        .map(|i| {
            let scale = 1.0 + i as f64;
            vec![
                ShearShellRelaxationBranch::new_coupled(
                    [4e4, 3e3, 1e3, 3e3, 7e4, 2e3, 1e3, 2e3, 2e4].map(|x| scale * x),
                    [8.0, 1.0, 0.0, 1.0, -12.0, 0.0, 0.0, 0.0, 4.0],
                    [6.0, 0.5, 0.2, 0.5, 3.0, 0.1, 0.2, 0.1, 2.0],
                    [2e4, 3e3, 5e4],
                    0.0003,
                )
                .unwrap(),
                ShearShellRelaxationBranch::new(
                    [0.0; 9],
                    [0.0; 9],
                    [3e4, -2e3, 1e4],
                    if i == 0 { 0.002 } else { 0.004 },
                )
                .unwrap(),
            ]
        })
        .collect();
    Fixture {
        mesh,
        directors,
        sections,
        model,
        pairs,
        source,
    }
}

fn limits() -> PlateRelaxationLimits {
    PlateRelaxationLimits {
        max_modes: 12,
        max_projection_terms: 100_000_000,
        max_memory_states: 64,
    }
}
fn modes(f: &Fixture) -> Vec<Vec<f64>> {
    f.pairs.iter().map(|m| m.phi.clone()).collect()
}
fn image(f: &Fixture, gate: &CancelGate) -> ModalShellRelaxation {
    ModalShellRelaxation::from_mitc_shell(&f.model, &f.source, &modes(f), limits(), gate).unwrap()
}
fn near(actual: f64, expected: f64, rel: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= rel * scale.max(1e-25),
        "{actual:e} != {expected:e}, scale={scale:e}"
    );
}

// Observe the complete original displacement through the independent public
// Bernstein strain observer, then contract original source tensors directly.
fn source_storage(f: &Fixture, q: &[f64], time: f64) -> f64 {
    let mut reduced = vec![0.0; f.model.pencil().free];
    for (mode, value) in f.pairs.iter().zip(q) {
        for (out, phi) in reduced.iter_mut().zip(&mode.phi) {
            *out += phi * value;
        }
    }
    let physical = f.model.expand(&reduced).unwrap();
    let mut energy = 0.0;
    for (e, tri) in f.mesh.tris.iter().enumerate() {
        let observer =
            PreparedElementStrain::new(&f.mesh, e, &f.directors, &f.sections[e], &mut || false)
                .unwrap();
        let field = observer
            .real_field(
                &core::array::from_fn(|i| physical[6 * tri[i / 6] + i % 6]),
                &mut || false,
            )
            .unwrap();
        let area = f.mesh.facet(e).unwrap().area_m2;
        for (lambda, weight) in integration_rule() {
            let strain = field.at_barycentric(lambda).unwrap();
            let m = strain.membrane;
            let k = strain.gradient_per_m;
            let s = strain.transverse_shear;
            for branch in &f.source[e] {
                if branch.relaxation_time_s() != time {
                    continue;
                }
                let (a, b, d, shear) = branch.coupled_matrices();
                let mut power = 0.0;
                for i in 0..3 {
                    for j in 0..3 {
                        power += m[i] * a[3 * i + j] * m[j]
                            + m[i] * b[3 * i + j] * k[j]
                            + k[i] * b[3 * i + j] * m[j]
                            + k[i] * d[3 * i + j] * k[j];
                    }
                }
                for i in 0..2 {
                    for j in 0..2 {
                        power += s[i] * shear[i][j] * s[j];
                    }
                }
                energy += 0.5 * area * weight * power;
            }
        }
    }
    energy
}

#[test]
fn mitc_modal_memory_preserves_original_mass_frames_and_full_source_storage() {
    for mass in [ShellMass::Consistent, ShellMass::Lumped] {
        let f = fixture(mass);
        let gate = CancelGate::new_clock_free();
        let image = image(&f, &gate);
        assert_eq!(image.projected_branches().len(), 3);
        assert_eq!(image.mode_count(), 5);
        let q = [0.02, -0.01, 0.007, -0.003, 0.004];
        for branch in image.projected_branches() {
            let expected = source_storage(&f, &q, branch.relaxation_time_s());
            let channels = 0.5
                * image
                    .coupling_rows_per_s()
                    .iter()
                    .zip(image.relaxation_times_s())
                    .filter(|(_, time)| **time == branch.relaxation_time_s())
                    .map(|(row, _)| row.iter().zip(q).map(|(c, x)| c * x).sum::<f64>().powi(2))
                    .sum::<f64>();
            let n = q.len();
            let modal = 0.5
                * (0..n)
                    .map(|i| {
                        (0..n)
                            .map(|j| q[i] * branch.matrix()[n * i + j] * q[j])
                            .sum::<f64>()
                    })
                    .sum::<f64>();
            assert!(expected > 0.0);
            near(channels, expected, 3e-9, expected);
            near(modal, expected, 3e-9, expected);
        }
        let reversed: Vec<_> = modes(&f)
            .into_iter()
            .rev()
            .map(|v| v.into_iter().map(|x| -x).collect())
            .collect();
        let changed =
            ModalShellRelaxation::from_mitc_shell(&f.model, &f.source, &reversed, limits(), &gate)
                .unwrap();
        let n = image.mode_count();
        for (a, b) in image
            .projected_branches()
            .iter()
            .zip(changed.projected_branches())
        {
            let scale = a.matrix().iter().map(|x| x.abs()).fold(0.0, f64::max);
            for i in 0..n {
                for j in 0..n {
                    near(
                        b.matrix()[n * i + j],
                        a.matrix()[n * (n - 1 - i) + n - 1 - j],
                        2e-12,
                        scale,
                    );
                }
            }
        }
    }
}

#[test]
fn mitc_spectra_use_existing_exact_memory_with_independent_heat_and_both_histories() {
    let f = fixture(ShellMass::Consistent);
    let gate = CancelGate::new_clock_free();
    let image = image(&f, &gate);
    let omega2: Vec<_> = f.pairs.iter().map(|mode| mode.lambda).collect();
    let system = phs(&omega2, &image, &gate);
    let n = image.mode_count();
    let states = system.state_dim();
    // Resolve the fastest relaxed mechanical mode, while the exact shared
    // transition handles every independent material rate without substeps.
    let dt = 0.6 / omega2.iter().copied().fold(0.0, f64::max).sqrt();
    let prepared = PreparedExactZoh::new(
        &system,
        dt,
        LinearSplitBudget {
            max_modes: states,
            max_omega_dt: 1.0,
        },
    )
    .unwrap();
    let force = [0.1, -0.03, 0.04, -0.02, 0.01];
    for relaxed in [false, true] {
        let mut initial = vec![0.0; states];
        for i in 0..n {
            initial[2 * i] = (i as f64 + 1.0) * 1e-6;
            initial[2 * i + 1] = if i % 2 == 0 { 0.01 } else { -0.02 };
        }
        if relaxed {
            for (arm, row) in image.coupling_rows_per_s().iter().enumerate() {
                initial[2 * n + arm] = row
                    .iter()
                    .enumerate()
                    .map(|(i, c)| c * initial[2 * i])
                    .sum();
            }
        }
        let mechanical = 0.5
            * (0..n)
                .map(|i| omega2[i] * initial[2 * i].powi(2) + initial[2 * i + 1].powi(2))
                .sum::<f64>();
        let start = system.storage().hamiltonian(&initial);
        if relaxed {
            near(start, mechanical, 3e-13, mechanical);
        } else {
            assert!(start > mechanical);
        }
        let mut exact = initial.clone();
        let (mut heat, mut work) = (0.0, 0.0);
        for _ in 0..48 {
            let step = prepared.step(&exact, &force).unwrap();
            assert!(step.dissipated >= 0.0);
            near(step.balance_residual(), 0.0, 3e-9, start);
            heat += step.dissipated;
            work += step.supplied;
            exact = step.x;
        }
        assert!(heat > 0.0);
        near(
            system.storage().hamiltonian(&exact) + heat - work,
            start,
            3e-9,
            start,
        );
        let mut errors = Vec::new();
        for subdivisions in [16, 32, 64] {
            let mut reference = initial.clone();
            reference.extend([0.0, 0.0]);
            for _ in 0..48 * subdivisions {
                rk4(
                    &mut reference,
                    &omega2,
                    &image,
                    &force,
                    dt / subdivisions as f64,
                );
            }
            let error = exact
                .iter()
                .zip(&reference)
                .enumerate()
                .map(|(i, (a, b))| {
                    let weight = if i < 2 * n && i % 2 == 0 {
                        omega2[i / 2].sqrt()
                    } else {
                        1.0
                    };
                    (weight * (a - b)).powi(2)
                })
                .sum::<f64>()
                .sqrt();
            errors.push(error);
            if subdivisions == 64 {
                near(reference[states], heat, 2e-7, heat);
                near(reference[states + 1], work, 2e-7, work.abs());
                assert!(error < 2e-7 * (2.0 * start).sqrt());
            }
        }
        assert!(errors[1] < 0.15 * errors[0], "RK4 errors {errors:?}");
        assert!(errors[2] < 0.15 * errors[1], "RK4 errors {errors:?}");
    }
}

#[test]
fn mitc_memory_refuses_wrong_basis_missing_facets_and_unbounded_work() {
    let f = fixture(ShellMass::Consistent);
    let gate = CancelGate::new_clock_free();
    let phi = modes(&f);
    let mut wrong = phi.clone();
    wrong[0].iter_mut().for_each(|x| *x *= 1.001);
    assert!(
        ModalShellRelaxation::from_mitc_shell(&f.model, &f.source, &wrong, limits(), &gate)
            .is_err()
    );
    assert!(
        ModalShellRelaxation::from_mitc_shell(&f.model, &f.source[..1], &phi, limits(), &gate)
            .is_err()
    );
    for bound in [
        PlateRelaxationLimits {
            max_modes: 4,
            ..limits()
        },
        PlateRelaxationLimits {
            max_memory_states: 0,
            ..limits()
        },
        PlateRelaxationLimits {
            max_projection_terms: 1,
            ..limits()
        },
    ] {
        assert!(
            ModalShellRelaxation::from_mitc_shell(&f.model, &f.source, &phi, bound, &gate).is_err()
        );
    }
    let empty = ModalShellRelaxation::from_mitc_shell(
        &f.model,
        &[vec![], vec![]],
        &phi,
        PlateRelaxationLimits {
            max_memory_states: 0,
            ..limits()
        },
        &gate,
    )
    .unwrap();
    assert!(empty.coupling_rows_per_s().is_empty());
    assert!(empty.branches(0, phi.len(), &gate).unwrap().is_empty());
    gate.request();
    assert!(matches!(
        ModalShellRelaxation::from_mitc_shell(&f.model, &f.source, &phi, limits(), &gate),
        Err(AcousticRealizeError::Cancelled)
    ));
}
