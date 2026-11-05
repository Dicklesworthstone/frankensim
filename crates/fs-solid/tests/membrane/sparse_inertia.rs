//! Original density-weighted inertia through the existing sparse direct owner.
use super::{close, sheet};
use fs_solid::membrane::{
    MembraneBudget, MembraneDomain, MembraneError, MembranePatch,
    dynamics::{ConsistentMassBudget, MembraneMotion, SparseConsistentMassBudget},
    inertia::MembraneMassImage,
};
use fs_solid::mesh2::Mesh2;
use fs_sparse::{direct::DirectOrdering, ordering_nd::NdOptions};

fn budget() -> SparseConsistentMassBudget {
    SparseConsistentMassBudget {
        max_nodes_per_axis: 2048,
        max_staged_entries_per_axis: 65_536,
        max_symbolic_bytes_per_factor: 32 * 1024 * 1024,
        max_factor_bytes_per_factor: 64 * 1024 * 1024,
        max_factor_nonzeros_per_axis: 250_000,
        max_estimated_factor_flops_per_factor: 2_000_000_000,
    }
}

fn nd() -> DirectOrdering {
    DirectOrdering::NestedDissection(NdOptions {
        leaf_size: 16,
        refinement_passes: 2,
    })
}

fn on_mesh(mesh: Mesh2) -> MembranePatch {
    let elements = mesh.elems.len();
    let nodes = mesh.nodes.len();
    MembranePatch::try_new(
        mesh,
        vec![sheet(8000.0, 2000.0, 500.0, 0.2); elements],
        vec![0.27; elements],
        MembraneDomain {
            minimum_area_ratio: 0.01,
            maximum_abs_strain: 0.2,
        },
        MembraneBudget {
            max_nodes: nodes,
            max_triangles: elements,
        },
        || false,
    )
    .unwrap()
}

fn momentum(motion: &MembraneMotion<'_>) -> [f64; 3] {
    let nodal = motion
        .patch()
        .reference_mass_apply(MembraneMassImage::Consistent, motion.velocity_m_s(), || {
            false
        })
        .unwrap();
    std::array::from_fn(|axis| nodal.iter().map(|p| p[axis]).sum())
}

#[test]
fn all_orderings_preserve_graded_nonlinear_motion_and_supported_mass_impulse() {
    let patch = super::patch(8000.0, 2000.0, 500.0, 0.2, 0.37)
        .with_reference_areal_density(
            vec![
                [2.0, 5.0, 7.0],
                [0.3, 0.9, 0.4],
                [9.0, 4.0, 3.0],
                [0.02, 0.08, 0.1],
            ],
            || false,
        )
        .unwrap();
    let u = super::displaced();
    let v = [
        [0.0, 0.0, 0.03],
        [0.01, -0.02, 0.0],
        [0.02, 0.01, -0.01],
        [-0.01, 0.03, 0.02],
        [0.01, -0.01, 0.04],
    ];
    // x/y share one support pattern; z retains a different principal block.
    let fixed = [(0, u[0][0]), (1, u[0][1]), (5, u[1][2])];
    let load = [
        [0.02, -0.01, 0.03],
        [0.0, 0.01, -0.02],
        [-0.01, 0.0, 0.02],
        [0.01, -0.02, 0.01],
        [0.02, 0.01, -0.01],
    ];
    let dense = MembraneMotion::new_with_consistent_inertia(
        &patch,
        &fixed,
        &u,
        &v,
        0.0,
        ConsistentMassBudget {
            max_nodes_per_axis: 5,
        },
        || false,
    )
    .unwrap();
    assert!(dense.sparse_inertia_statistics().is_none());
    for ordering in [DirectOrdering::Natural, DirectOrdering::Amd, nd()] {
        let mut expected = dense.clone();
        let mut motion = MembraneMotion::new_with_sparse_consistent_inertia(
            &patch,
            &fixed,
            &u,
            &v,
            0.0,
            budget(),
            ordering,
            || false,
        )
        .unwrap();
        assert_eq!(motion.inertia_image(), MembraneMassImage::Consistent);
        let stats = motion.sparse_inertia_statistics().unwrap();
        assert_eq!(stats.distinct_factors, 2);
        assert!(stats.factor_nonzeros_per_inverse > stats.stored_factor_nonzeros);
        let dt = 1e-6;
        for _ in 0..8 {
            let old_p = momentum(&motion);
            let a = motion.step(dt, &load, || false).unwrap();
            let b = expected.step(dt, &load, || false).unwrap();
            assert_eq!(a.time_s, b.time_s);
            for (actual, reference) in [
                (a.strain_energy_j, b.strain_energy_j),
                (a.kinetic_energy_j, b.kinetic_energy_j),
                (a.dead_load_work_j, b.dead_load_work_j),
                (a.energy_defect_j, b.energy_defect_j),
            ] {
                close(actual, reference, 1e-12, 2e-11);
            }
            let new_p = momentum(&motion);
            for axis in 0..3 {
                close(
                    a.support_impulse_n_s[axis],
                    b.support_impulse_n_s[axis],
                    1e-13,
                    2e-11,
                );
                let applied = dt * load.iter().map(|f| f[axis]).sum::<f64>();
                close(
                    new_p[axis] - old_p[axis],
                    applied + a.support_impulse_n_s[axis],
                    1e-12,
                    2e-11,
                );
            }
            for (actual, reference) in motion
                .displacement_m()
                .iter()
                .flatten()
                .zip(expected.displacement_m().iter().flatten())
            {
                close(*actual, *reference, 1e-14, 2e-11);
            }
            for (actual, reference) in motion
                .velocity_m_s()
                .iter()
                .flatten()
                .zip(expected.velocity_m_s().iter().flatten())
            {
                close(*actual, *reference, 1e-11, 2e-11);
            }
            for &(dof, value) in &fixed {
                assert_eq!(motion.displacement_m()[dof / 3][dof % 3], value);
                assert_eq!(motion.velocity_m_s()[dof / 3][dof % 3], 0.0);
            }
            assert_eq!(motion.sparse_inertia_statistics(), Some(stats));
        }
    }
}

#[test]
fn sparse_inertia_exceeds_the_dense_leaf_with_original_uniform_acceleration() {
    let cells = 24;
    let side = cells + 1;
    let mut mesh = Mesh2 {
        nodes: Vec::new(),
        elems: Vec::new(),
        patches: Vec::new(),
    };
    for j in 0..=cells {
        for i in 0..=cells {
            mesh.nodes
                .push([i as f64 / cells as f64, j as f64 / cells as f64]);
        }
    }
    for j in 0..cells {
        for i in 0..cells {
            let n = j * side + i;
            mesh.elems.extend([
                vec![n, n + 1, n + side],
                vec![n + 1, n + side + 1, n + side],
            ]);
        }
    }
    let density: Vec<_> = mesh
        .elems
        .iter()
        .enumerate()
        .map(|(e, nodes)| {
            std::array::from_fn(|a| {
                let [x, y] = mesh.nodes[nodes[a]];
                0.2 + 0.1 * x + 0.3 * y + 0.05 * (e % 3) as f64
            })
        })
        .collect();
    let patch = on_mesh(mesh)
        .with_reference_areal_density(density, || false)
        .unwrap();
    let nodes = patch.mesh().nodes.len();
    assert_eq!(nodes, 625);
    let zero = vec![[0.0; 3]; nodes];
    let initial_v = [0.07, -0.03, 0.05];
    let initial = vec![initial_v; nodes];
    let acceleration = [0.3, -0.2, 0.7];
    let load = patch.body_load_n(acceleration).unwrap();
    assert!(matches!(
        MembraneMotion::new_with_consistent_inertia(
            &patch,
            &[],
            &zero,
            &initial,
            0.0,
            ConsistentMassBudget {
                max_nodes_per_axis: nodes
            },
            || false,
        ),
        Err(MembraneError::Budget(_))
    ));
    let mut motion = MembraneMotion::new_with_sparse_consistent_inertia(
        &patch,
        &[],
        &zero,
        &initial,
        0.0,
        budget(),
        nd(),
        || false,
    )
    .unwrap();
    let stats = motion.sparse_inertia_statistics().unwrap();
    assert_eq!(stats.distinct_factors, 1);
    assert_eq!(
        stats.factor_nonzeros_per_inverse,
        3 * stats.stored_factor_nonzeros
    );
    assert!(stats.scalar_mass_nonzeros < 8 * nodes);
    assert!(stats.stored_factor_nonzeros < nodes * nodes / 4);
    let dt = 1e-5;
    for _ in 0..3 {
        let old_ke = motion.kinetic_energy_j();
        let report = motion.step(dt, &load, || false).unwrap();
        let time = motion.time_s();
        assert_eq!(report.support_impulse_n_s, [0.0; 3]);
        let mut speed_squared = 0.0;
        for axis in 0..3 {
            let velocity = initial_v[axis] + time * acceleration[axis];
            speed_squared += velocity * velocity;
            let displacement = time * initial_v[axis] + 0.5 * time * time * acceleration[axis];
            for node in 0..nodes {
                close(motion.velocity_m_s()[node][axis], velocity, 2e-10, 2e-10);
                close(
                    motion.displacement_m()[node][axis],
                    displacement,
                    1e-14,
                    2e-10,
                );
            }
        }
        close(
            report.kinetic_energy_j,
            0.5 * patch.reference_mass_kg() * speed_squared,
            1e-12,
            2e-10,
        );
        close(
            report.kinetic_energy_j - old_ke,
            report.dead_load_work_j,
            1e-12,
            2e-10,
        );
        close(report.strain_energy_j, 0.0, 1e-18, 0.0);
        close(report.energy_defect_j, 0.0, 1e-12, 0.0);
    }
}

#[test]
fn sparse_caps_and_late_cancellation_keep_checkpoint_replay_atomic() {
    let patch = on_mesh(Mesh2 {
        nodes: vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        elems: vec![vec![0, 1, 2]],
        patches: vec![],
    })
    .with_reference_areal_density(vec![[2.0, 5.0, 7.0]], || false)
    .unwrap();
    let u = [[0.0; 3], [0.001, -0.001, 0.001], [-0.001, 0.002, -0.001]];
    let v = [[0.01, -0.02, 0.03]; 3];
    let zero = [[0.0; 3]; 3];
    let full = budget();
    for limited in [
        SparseConsistentMassBudget {
            max_nodes_per_axis: 2,
            ..full
        },
        SparseConsistentMassBudget {
            max_staged_entries_per_axis: 8,
            ..full
        },
        SparseConsistentMassBudget {
            max_symbolic_bytes_per_factor: 1,
            ..full
        },
        SparseConsistentMassBudget {
            max_factor_bytes_per_factor: 1,
            ..full
        },
        SparseConsistentMassBudget {
            max_factor_nonzeros_per_axis: 1,
            ..full
        },
        SparseConsistentMassBudget {
            max_estimated_factor_flops_per_factor: 1,
            ..full
        },
        SparseConsistentMassBudget {
            max_nodes_per_axis: 0,
            ..full
        },
    ] {
        assert!(matches!(
            MembraneMotion::new_with_sparse_consistent_inertia(
                &patch,
                &[],
                &u,
                &v,
                0.0,
                limited,
                DirectOrdering::Amd,
                || false,
            ),
            Err(MembraneError::Budget(_))
        ));
    }
    let mut polls = 0;
    let mut motion = MembraneMotion::new_with_sparse_consistent_inertia(
        &patch,
        &[],
        &u,
        &v,
        0.0,
        full,
        nd(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(matches!(
            MembraneMotion::new_with_sparse_consistent_inertia(
                &patch,
                &[],
                &u,
                &v,
                0.0,
                full,
                nd(),
                || {
                    seen += 1;
                    seen == stop
                },
            ),
            Err(MembraneError::Cancelled)
        ));
    }
    let before = motion.clone();
    let mut replay = motion.clone();
    let mut step_polls = 0;
    let expected = replay
        .step(1e-5, &zero, || {
            step_polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        motion.step(1e-5, &zero, || {
            seen += 1;
            seen == step_polls
        }),
        Err(MembraneError::Cancelled)
    ));
    assert_eq!(motion.displacement_m(), before.displacement_m());
    assert_eq!(motion.velocity_m_s(), before.velocity_m_s());
    assert_eq!(motion.time_s(), before.time_s());
    assert_eq!(motion.kinetic_energy_j(), before.kinetic_energy_j());
    assert_eq!(
        motion.sparse_inertia_statistics(),
        before.sparse_inertia_statistics()
    );
    assert_eq!(motion.step(1e-5, &zero, || false).unwrap(), expected);
    assert_eq!(motion.displacement_m(), replay.displacement_m());
    assert_eq!(motion.velocity_m_s(), replay.velocity_m_s());
}

#[test]
fn sparse_row_scaling_retains_both_extreme_physical_inventories() {
    let patch = on_mesh(Mesh2 {
        nodes: vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [2.0, 1.0],
        ],
        elems: vec![vec![0, 1, 2], vec![3, 4, 5]],
        patches: vec![],
    })
    .with_reference_areal_density(
        vec![[2e-280, 5e-280, 7e-280], [2e280, 5e280, 7e280]],
        || false,
    )
    .unwrap();
    let velocity = [
        [1e140, -2e140, 3e140],
        [1e140, -2e140, 3e140],
        [1e140, -2e140, 3e140],
        [1e-140, -2e-140, 3e-140],
        [1e-140, -2e-140, 3e-140],
        [1e-140, -2e-140, 3e-140],
    ];
    let zero = [[0.0; 3]; 6];
    let mut motion = MembraneMotion::new_with_sparse_consistent_inertia(
        &patch,
        &[],
        &zero,
        &velocity,
        0.0,
        budget(),
        nd(),
        || false,
    )
    .unwrap();
    // Each disconnected rigid translation carries 49/3 J, even though their
    // reference mass ratio spans 560 orders of magnitude.
    close(motion.kinetic_energy_j(), 98.0 / 3.0, 1e-13, 2e-12);
    for (actual, supplied) in motion
        .velocity_m_s()
        .iter()
        .flatten()
        .zip(velocity.iter().flatten())
    {
        close(*actual, *supplied, 0.0, 2e-12);
    }
    let dt = 1e-145;
    let report = motion.step(dt, &zero, || false).unwrap();
    close(report.kinetic_energy_j, 98.0 / 3.0, 1e-13, 2e-12);
    for (actual, supplied) in motion
        .displacement_m()
        .iter()
        .flatten()
        .zip(velocity.iter().flatten())
    {
        close(*actual, dt * supplied, 0.0, 2e-12);
    }
}
