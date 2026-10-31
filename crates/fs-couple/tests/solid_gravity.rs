//! Native shared gravity on finite-deformation solid quadrature and supports.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::solid_gravity::bind_field;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_mbd::gravity::{CentralGravity, UniformGravity};
use fs_mbd::{Gravity, Vec3};
use fs_solid::{HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField};

const FRAME: u64 = 17;
const NODES: [[f64; 3]; 4] = [
    [4.0, 0.0, 0.0],
    [5.0, 0.0, 0.0],
    [4.0, 1.0, 0.0],
    [4.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn material() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2.0,
            lambda: 3.0,
        },
        5.0,
    )
    .unwrap()
}

fn problem<'a>(law: &'a Hyperelastic, prescribed: &'a [(usize, f64)]) -> HyperTetProblem<'a> {
    HyperTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        material: HyperelasticLawField::Uniform(law),
        prescribed_m: prescribed,
        nodal_forces_n: &[],
        body_force_n_m3: [0.0; 3],
        budget: Default::default(),
        settings: HyperTetSettings::default(),
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn shared_central_field_supplies_current_position_force_potential_and_newton_tangent() {
    let source = CentralGravity::new(FRAME, Vec3::ZERO, 20.0).unwrap();
    let field = bind_field(&source, FRAME);
    let law = material();
    let p = problem(&law, &[]);
    let u = [
        0.1, -0.2, 0.3, 0.25, -0.15, 0.2, 0.2, -0.05, 0.25, 0.0, -0.1, 0.4,
    ];
    with_cx(|cx| {
        let prepared = p.prepare(cx).unwrap();
        let initial = prepared
            .gravity_load(&[6.0], &[0.0; 12], &field, cx)
            .unwrap();
        let moved = prepared.gravity_load(&[6.0], &u, &field, cx).unwrap();
        assert_ne!(initial.force_n, moved.force_n);
        for column in 0..12 {
            let mut plus = u;
            let mut minus = u;
            plus[column] += 1e-5;
            minus[column] -= 1e-5;
            let a = prepared.gravity_load(&[6.0], &plus, &field, cx).unwrap();
            let b = prepared.gravity_load(&[6.0], &minus, &field, cx).unwrap();
            near(
                -(a.potential_energy_j - b.potential_energy_j) / 2e-5,
                moved.force_n[column],
                2e-10,
            );
            for row in 0..12 {
                near(
                    (a.force_n[row] - b.force_n[row]) / 2e-5,
                    moved.force_tangent_n_m.get(row, column),
                    2e-11,
                );
            }
        }
        let wrong_frame = bind_field(&source, FRAME + 1);
        assert!(matches!(
            prepared.gravity_load(&[6.0], &u, &wrong_frame, cx),
            Err(HyperTetError::InvalidInput { .. })
        ));
        let repeated = prepared.gravity_load(&[6.0], &u, &field, cx).unwrap();
        assert_eq!(moved.force_n, repeated.force_n);
        assert_eq!(moved.potential_energy_j, repeated.potential_energy_j);
    });
}

#[test]
fn shared_uniform_field_keeps_reference_weight_and_complete_support_reactions() {
    let acceleration = Vec3::new(0.2, -0.3, -1.7);
    let source = UniformGravity::new(
        FRAME,
        Gravity::new(acceleration).unwrap(),
        Vec3::new(10.0, -3.0, 1.0),
    )
    .unwrap();
    let field = bind_field(&source, FRAME);
    let law = material();
    let fixed: Vec<_> = (0..12).map(|dof| (dof, 0.0)).collect();
    let p = problem(&law, &fixed);
    with_cx(|cx| {
        let solution = p.solve_with_gravity(&[6.0], &field, cx).unwrap();
        // Unit-volume simplex / 6, density 6: exactly one kilogram.
        let gravity = [acceleration.x, acceleration.y, acceleration.z];
        for (axis, g) in gravity.into_iter().enumerate() {
            let reaction: f64 = solution
                .reactions_n
                .iter()
                .filter(|(dof, _)| dof % 3 == axis)
                .map(|(_, force)| force)
                .sum();
            near(reaction, -g, 2e-14);
        }
        // Stretch by two in every direction: eight times current volume,
        // but the same original material inventory and total weight.
        let stretched: Vec<_> = NODES.iter().flatten().copied().collect();
        let load = p.gravity_load(&[6.0], &stretched, &field, cx).unwrap();
        for (axis, g) in gravity.into_iter().enumerate() {
            near(
                load.force_n.chunks_exact(3).map(|f| f[axis]).sum(),
                g,
                2e-14,
            );
        }
        assert!((0..12).all(|row| load.force_tangent_n_m.row(row).1.iter().all(|&v| v == 0.0)));
    });
}
