//! Full original MITC coordinates and causal resultant memory, without modes.
use fs_exec::{CancelGate, Cx};
use fs_material::elastic::OrthotropicElastic;
use fs_material::laminate::{LaminaLayer, LaminateSection, LayerLaw, ReferenceSurface};
use fs_math::det;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::relaxation::{
    ShearShellRelaxationBranch, ShearShellRelaxationBudget, project_relaxation,
};
use fs_plate::shell::mitc::transient::{
    MitcTransientBudget, MitcTransientDomain, MitcTransientError, MitcTransientSource,
};
use fs_plate::shell::mitc::viscosity::{ShearShellViscosity, project_viscosity};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellSection,
    ShellConstraint, ShellMass, assemble_mitc3p,
};

fn with_cx(run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 79,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ));
    });
}
fn mesh_budget() -> MitcBudget {
    MitcBudget {
        max_nodes: 3,
        max_triangles: 1,
    }
}
fn budget() -> MitcTransientBudget {
    MitcTransientBudget {
        maximum_elements: 1,
        maximum_coordinates: 18,
        maximum_memory_channels: 196,
        maximum_payload_bytes: 1 << 28,
        maximum_work: 10_000_000_000,
    }
}
fn domain() -> [MitcTransientDomain; 1] {
    [MitcTransientDomain {
        faces_m: [-0.01, 0.03],
        maximum_reduced_strain: 0.04,
    }]
}
fn make_model(mask: [bool; 6], mass: ShellMass) -> MitcShellModel {
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [0.4, 0.03, 0.05], [0.07, 0.3, 0.04]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let normal = mesh.facet(0).unwrap().frame[2];
    let normals = (0..3)
        .map(|i| {
            let mut n = normal;
            n[0] += 0.01 * (i as f64 - 1.0);
            n[1] += 0.008 * (i as f64 - 0.5);
            let length = det::sqrt(n.iter().map(|v| v * v).sum());
            n.map(|v| v / length)
        })
        .collect();
    let directors = NodalDirectors::new(&mesh, normals, 0.9).unwrap();
    let layers = [0.03, 0.01]
        .into_iter()
        .enumerate()
        .map(|(i, thickness)| LaminaLayer {
            thickness_m: thickness,
            orientation_rad: if i == 0 { 0.37 } else { -0.61 },
            density_kg_m3: if i == 0 { 450.0 } else { 1100.0 },
            material: None,
            law: LayerLaw::Elastic(
                OrthotropicElastic::new(
                    [14000.0, 7000.0, 10000.0],
                    [0.21, 0.27, 0.18],
                    [3000.0, 2000.0, 4000.0],
                    0.04,
                )
                .unwrap(),
            ),
        })
        .collect();
    let section = LaminateSection::try_new(layers, ReferenceSurface::Offset(0.01)).unwrap();
    let (s, c) = (det::sin(0.43), det::cos(0.43));
    let support = [
        ShellConstraint {
            node: 0,
            mask: [true; 6],
            frame: ConstraintFrame::Global,
        },
        ShellConstraint {
            node: 1,
            mask: [true; 6],
            frame: ConstraintFrame::Global,
        },
        ShellConstraint {
            node: 2,
            mask,
            frame: ConstraintFrame::Local([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]]),
        },
    ];
    assemble_mitc3p(
        &mesh,
        &directors,
        &[ShearShellSection::from_laminate(&section, 0.29).unwrap()],
        &support,
        mass,
        mesh_budget(),
        0.08,
        || false,
    )
    .unwrap()
}
fn diagonal([x, y, z]: [f64; 3]) -> [f64; 9] {
    [x, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, z]
}
fn zero() -> ShearShellViscosity {
    ShearShellViscosity::new([0.0; 9], [0.0; 9], [0.0; 3]).unwrap()
}
fn spectra() -> Vec<Vec<ShearShellRelaxationBranch>> {
    vec![vec![
        ShearShellRelaxationBranch::new_coupled(
            diagonal([15.0, 7.0, 3.0]),
            diagonal([0.02, -0.015, 0.005]),
            diagonal([0.0003, 0.0002, 0.0001]),
            [2.0, 0.2, 1.0],
            0.014,
        )
        .unwrap(),
        // Same time, independent material history, and a directional shear arm.
        ShearShellRelaxationBranch::new([0.0; 9], [0.0; 9], [0.3, 0.0, 0.0], 0.014).unwrap(),
    ]]
}
fn basis(n: usize) -> Vec<Vec<f64>> {
    (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-13),
        "{a:e} != {b:e}"
    );
}
fn form(matrix: &[f64], x: &[f64]) -> f64 {
    x.iter()
        .enumerate()
        .map(|(i, a)| {
            x.iter()
                .enumerate()
                .map(|(j, b)| a * matrix[x.len() * i + j] * b)
                .sum::<f64>()
        })
        .sum()
}
fn dense(matrix: &fs_sparse::Csr, n: usize) -> Vec<f64> {
    (0..n)
        .flat_map(|i| (0..n).map(move |j| matrix.get(i, j)))
        .collect()
}

#[test]
fn full_coordinates_preserve_native_resultant_work_inertia_frames_and_independent_memory() {
    with_cx(|cx| {
        for rotations in [false, true] {
            let mask = if rotations {
                [false; 6]
            } else {
                [false, false, false, true, true, true]
            };
            let model = make_model(mask, ShellMass::Consistent);
            let n = model.pencil().free;
            let viscosity = if rotations {
                zero()
            } else {
                ShearShellViscosity::new_coupled(
                    diagonal([0.2, 0.3, 0.1]),
                    diagonal([0.0001, -0.0001, 0.0001]),
                    diagonal([0.0002, 0.0003, 0.0001]),
                    [0.1, 0.01, 0.2],
                )
                .unwrap()
            };
            let sources = spectra();
            let source = MitcTransientSource::prepare(
                &model,
                &[viscosity],
                &sources,
                &domain(),
                budget(),
                cx,
            )
            .unwrap();
            assert_eq!(source.mechanics().coordinates(), n);
            assert_eq!(
                source.mechanics().system().state_dim(),
                2 * n + source.memory_channels().len()
            );
            assert!(
                source
                    .memory_channels()
                    .iter()
                    .any(|r| r.source_branch == 0)
            );
            assert!(
                source
                    .memory_channels()
                    .iter()
                    .any(|r| r.source_branch == 1)
            );
            let basis = basis(n);
            let reference = project_relaxation(
                &model,
                &sources,
                &basis,
                ShearShellRelaxationBudget {
                    mesh: mesh_budget(),
                    max_modes: n,
                    max_element_branches: 2,
                    max_relaxation_times: 1,
                    max_projection_terms: usize::MAX,
                },
                || false,
            )
            .unwrap();
            assert_eq!(reference.len(), 1);
            let mut gram = vec![0.0; n * n];
            for branch in source.mechanics().branches() {
                assert_eq!(branch.stiffness, 1.0);
                for i in 0..n {
                    for j in 0..n {
                        gram[n * i + j] += branch.projection[i] * branch.projection[j];
                    }
                }
            }
            let norm = reference[0]
                .matrix()
                .iter()
                .map(|v| v.abs())
                .fold(0.0, f64::max);
            for (a, b) in gram.iter().zip(reference[0].matrix()) {
                assert!((a - b).abs() < 3e-10 * norm.max(1e-12));
            }
            let q: Vec<_> = (0..n).map(|i| 1e-5 * (i as f64 - 1.7)).collect();
            let v: Vec<_> = (0..n).map(|i| 0.003 * (i as f64 - 2.3)).collect();
            let relaxed = source.relaxed_state(&q, &v, cx).unwrap();
            let loaded = source
                .state(&q, &v, &vec![0.0; source.memory_channels().len()], cx)
                .unwrap();
            let a = source.observe(&relaxed, cx).unwrap();
            let b = source.observe(&loaded, cx).unwrap();
            close(
                a.stored_energy_j,
                0.5 * form(&dense(&model.pencil().k, n), &q)
                    + 0.5 * form(&dense(&model.pencil().m, n), &v),
                3e-10,
            );
            close(
                b.stored_energy_j - a.stored_energy_j,
                0.5 * form(&gram, &q),
                2e-8,
            );
            assert_eq!(a.displacement, model.expand(&q).unwrap());
            for (a, b) in a.velocity.iter().zip(model.expand(&v).unwrap()) {
                close(*a, b, 2e-9);
            }
            let (mesh, directors, sections) = model.original_geometry();
            let physical: Vec<_> = basis.iter().map(|q| model.expand(q).unwrap()).collect();
            let c = project_viscosity(
                mesh,
                directors,
                sections,
                &[viscosity],
                &physical,
                mesh_budget(),
                usize::MAX,
                || false,
            )
            .unwrap();
            close(a.material_heat_w, form(&c, &v), 3e-9);
            close(
                b.material_heat_w - a.material_heat_w,
                form(&gram, &q) / 0.014,
                3e-9,
            );
            let velocity = source.mechanics().system().output(&loaded);
            for (a, b) in velocity.iter().zip(&v) {
                close(*a, *b, 2e-9);
            }
        }
    });
}

#[test]
fn native_time_owner_retains_causal_spring_memory_and_material_heat() {
    with_cx(|cx| {
        let model = make_model([false, true, true, true, true, true], ShellMass::Consistent);
        let branch =
            ShearShellRelaxationBranch::new(diagonal([70.0, 0.0, 0.0]), [0.0; 9], [0.0; 3], 0.04)
                .unwrap();
        let source = MitcTransientSource::prepare(
            &model,
            &[zero()],
            &[vec![branch]],
            &domain(),
            budget(),
            cx,
        )
        .unwrap();
        let mechanics = source.mechanics();
        let system = mechanics.system();
        let k = model.pencil().k.get(0, 0);
        let m = model.pencil().m.get(0, 0);
        let excess: f64 = mechanics
            .branches()
            .iter()
            .map(|b| b.projection[0].powi(2))
            .sum();
        let (mut q, mut v) = (1e-4, 0.0);
        let mut stress = excess * q;
        let mut state = source
            .state(&[q], &[v], &vec![0.0; source.memory_channels().len()], cx)
            .unwrap();
        let initial = system.hamiltonian(&state);
        let mut heat = 0.0;
        let mut supplied = 0.0;
        let dt = 0.001;
        for step in 0..160 {
            let force = if step < 20 { 0.002 } else { 0.0 };
            // Independent scalar midpoint balance after eliminating only the
            // analytic reference's single parallel Maxwell force.
            let a = dt / (2.0 * 0.04);
            let h = dt * dt * (k + excess / (1.0 + a)) / 4.0;
            let next_v =
                ((m - h) * v - dt * k * q - dt * stress / (1.0 + a) + dt * force) / (m + h);
            let next_q = q + 0.5 * dt * (v + next_v);
            stress = ((1.0 - a) * stress + 0.5 * dt * excess * (v + next_v)) / (1.0 + a);
            q = next_q;
            v = next_v;
            let record = fs_phs::step(system, &state, &[force], dt).unwrap();
            source.admit_state(&record.x, cx).unwrap();
            assert!(record.dissipated >= 0.0);
            assert!(record.balance_residual().abs() < 2e-12);
            heat += record.dissipated;
            supplied += record.supplied;
            state = record.x;
            close(state[0], q, 2e-7);
            close(mechanics.velocities(&state).unwrap()[0], v, 2e-7);
            let actual_stress: f64 = mechanics
                .branch_stresses(&state)
                .unwrap()
                .iter()
                .zip(mechanics.branches())
                .map(|(s, b)| s * b.projection[0])
                .sum();
            close(actual_stress, stress, 3e-7);
        }
        assert!(heat > 0.0);
        close(system.hamiltonian(&state) + heat, initial + supplied, 2e-8);
        let snapshot = state.clone();
        let malformed = &snapshot[..snapshot.len() - 1];
        assert!(source.observe(malformed, cx).is_err());
        assert_eq!(state, snapshot);
    });
}

#[test]
fn native_domains_cold_budgets_and_cancellation_refuse_before_publication() {
    with_cx(|cx| {
        let model = make_model([false; 6], ShellMass::Consistent);
        let source =
            MitcTransientSource::prepare(&model, &[zero()], &spectra(), &domain(), budget(), cx)
                .unwrap();
        for b in [
            MitcTransientBudget {
                maximum_coordinates: 5,
                ..budget()
            },
            MitcTransientBudget {
                maximum_memory_channels: 1,
                ..budget()
            },
            MitcTransientBudget {
                maximum_payload_bytes: 0,
                ..budget()
            },
            MitcTransientBudget {
                maximum_work: 0,
                ..budget()
            },
        ] {
            assert!(matches!(
                MitcTransientSource::prepare(&model, &[zero()], &spectra(), &domain(), b, cx),
                Err(MitcTransientError::Budget { .. })
            ));
        }
        assert!(matches!(
            source.relaxed_state(&[2.0; 6], &[0.0; 6], cx),
            Err(MitcTransientError::StrainDomain { .. })
        ));
        assert!(source.state(&[0.0; 6], &[0.0; 6], &[], cx).is_err());
        let lumped = make_model([false; 6], ShellMass::Lumped);
        assert!(
            MitcTransientSource::prepare(&lumped, &[zero()], &spectra(), &domain(), budget(), cx)
                .is_err()
        );
        let lossless =
            MitcTransientSource::prepare(&model, &[zero()], &[vec![]], &domain(), budget(), cx)
                .unwrap();
        assert!(lossless.memory_channels().is_empty());
        let state = lossless.relaxed_state(&[0.0; 6], &[0.001; 6], cx).unwrap();
        assert_eq!(lossless.observe(&state, cx).unwrap().material_heat_w, 0.0);
        let gate = CancelGate::new_clock_free();
        gate.request();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 79,
                    kernel_id: 29,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            assert!(matches!(
                lossless.observe(&state, &cancelled),
                Err(MitcTransientError::Plate(fs_plate::PlateError::Cancelled))
            ));
            assert!(
                MitcTransientSource::prepare(
                    &model,
                    &[zero()],
                    &spectra(),
                    &domain(),
                    budget(),
                    &cancelled
                )
                .is_err()
            );
        });
        assert_eq!(lossless.observe(&state, cx).unwrap().material_heat_w, 0.0);
    });
}
