//! Physical inventory continuation through the ORIGINAL native transport owner.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_flux::ns::transport::density::{
    DensityTransportCell as Cell, DensityTransportError, PreparedDensityTransport,
    remap::{
        DensityRemapControls as Controls, DensityRemapError as Error,
        DensityRemapQuantity as Quantity,
    },
};
use fs_flux::{FluxSystem, TriMesh};
use fs_mesh::{ConservativeRemapError, RemapContribution as Overlap};
use fs_solid::Mesh2;

use super::{counterflow, near};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig {
        chunk_bytes: 4096,
        max_chunk_bytes: 4096,
        limit_bytes: Some(8192),
        free_list_max_bytes: 0,
        ..Default::default()
    });
    pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 28,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn controls<const N: usize>() -> Controls<N> {
    Controls {
        source_coverage_tolerance: 0.0,
        relative_target_volume_tolerance: 0.0,
        volume_tolerance_m3: 0.0,
        mass_tolerance_kg: 1e-12,
        carried_tolerance: [1e-10; N],
        maximum_storage_bytes: 1 << 20,
    }
}

fn mesh(nodes: &[[f64; 2]], elems: &[[usize; 3]]) -> TriMesh {
    TriMesh::from_mesh2(&Mesh2 {
        nodes: nodes.to_vec(),
        elems: elems.iter().map(|tri| tri.to_vec()).collect(),
        patches: Vec::new(),
    })
}

fn square_meshes() -> (TriMesh, TriMesh) {
    let nodes = [[0., 0.], [1., 0.], [1., 1.], [0., 1.], [0.5, 0.5]];
    (
        mesh(&nodes[..4], &[[0, 1, 2], [0, 2, 3]]),
        mesh(&nodes, &[[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]),
    )
}

fn refine_map() -> [Overlap; 4] {
    [
        Overlap::new(0, 0, 0.5),
        Overlap::new(0, 1, 0.5),
        Overlap::new(1, 2, 0.5),
        Overlap::new(1, 3, 0.5),
    ]
}

fn totals(cells: &[Cell<3>], volumes: &[f64]) -> [f64; 4] {
    let mut result = [0.0; 4];
    for (cell, &volume) in cells.iter().zip(volumes) {
        let mass = cell.density_kg_m3 * volume;
        result[0] += mass;
        for component in 0..3 {
            result[1 + component] += mass * cell.specific[component];
        }
    }
    result
}

fn advance<const N: usize>(
    plan: &PreparedDensityTransport<'_>,
    cells: &[Cell<N>],
    edges: usize,
) -> Vec<Cell<N>> {
    let mut output = vec![Cell::default(); cells.len()];
    let mut scratch = plan.workspace::<N>(1 << 20).unwrap();
    plan.advance_into(
        cells,
        &vec![None; edges],
        0.25,
        0.8,
        &mut output,
        &mut scratch,
        || false,
    )
    .unwrap();
    output
}

#[test]
fn advance_refine_advance_coarsen_and_continue_with_enthalpy_and_species() {
    let (coarse_mesh, fine_mesh) = square_meshes();
    let coarse_system = FluxSystem::new(&coarse_mesh);
    let fine_system = FluxSystem::new(&fine_mesh);
    let coarse_edge = coarse_mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let fine_edge = fine_mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX && e.tris.0 < 2 && e.tris.1 >= 2)
        .unwrap();
    let coarse_advector = counterflow(&coarse_system, coarse_edge, 0.1);
    let fine_advector = counterflow(&fine_system, fine_edge, 0.1);
    let coarse = coarse_system
        .prepare_density_transport(&coarse_advector, 2.0, 1e-14, 1 << 20, || false)
        .unwrap();
    let fine = fine_system
        .prepare_density_transport(&fine_advector, 2.0, 1e-14, 1 << 20, || false)
        .unwrap();
    assert_eq!(coarse.cell_volumes_m3(), &[1.0, 1.0]);
    assert_eq!(fine.cell_volumes_m3(), &[0.5; 4]);
    let initial = [
        Cell {
            density_kg_m3: 2.0,
            specific: [10.0, 0.25, 0.75],
        },
        Cell {
            density_kg_m3: 8.0,
            specific: [100.0, 0.75, 0.25],
        },
    ];
    let original_totals = [10.0, 820.0, 6.5, 3.5];
    assert_eq!(totals(&initial, coarse.cell_volumes_m3()), original_totals);
    let accepted = advance(&coarse, &initial, coarse_mesh.edges.len());
    assert_ne!(accepted.as_slice(), &initial);
    let saved = accepted.clone();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let refined = fine
            .remap_from(&coarse, &accepted, &refine_map(), controls(), cx)
            .unwrap();
        let replay = fine
            .remap_from(&coarse, &accepted, &refine_map(), controls(), cx)
            .unwrap();
        assert_eq!(refined, replay);
        assert_eq!(accepted, saved);
        for (child, &parent) in [0, 0, 1, 1].iter().enumerate() {
            near(
                refined.cells[child].density_kg_m3,
                accepted[parent].density_kg_m3,
                1e-13,
            );
            for component in 0..3 {
                near(
                    refined.cells[child].specific[component],
                    accepted[parent].specific[component],
                    1e-12,
                );
            }
        }
        let evolved = advance(&fine, &refined.cells, fine_mesh.edges.len());
        assert_ne!(evolved, refined.cells);
        let coarsened = coarse
            .remap_from(
                &fine,
                &evolved,
                &[
                    Overlap::new(0, 0, 1.0),
                    Overlap::new(1, 0, 1.0),
                    Overlap::new(2, 1, 1.0),
                    Overlap::new(3, 1, 1.0),
                ],
                controls(),
                cx,
            )
            .unwrap();
        let continued = advance(&coarse, &coarsened.cells, coarse_mesh.edges.len());
        for (cells, volumes) in [
            (refined.cells.as_slice(), fine.cell_volumes_m3()),
            (evolved.as_slice(), fine.cell_volumes_m3()),
            (coarsened.cells.as_slice(), coarse.cell_volumes_m3()),
            (continued.as_slice(), coarse.cell_volumes_m3()),
        ] {
            for (actual, expected) in totals(cells, volumes).into_iter().zip(original_totals) {
                near(actual, expected, 2e-10);
            }
            for cell in cells {
                assert!(cell.density_kg_m3 > 0.0);
                assert!(cell.specific[1] >= 0.0 && cell.specific[2] >= 0.0);
                near(cell.specific[1] + cell.specific[2], 1.0, 4e-15);
            }
        }
        near(
            refined.report.mass_after_kg,
            totals(&refined.cells, fine.cell_volumes_m3())[0],
            1e-14,
        );
        for component in 0..3 {
            near(
                refined.report.carried_after[component],
                totals(&refined.cells, fine.cell_volumes_m3())[1 + component],
                1e-12,
            );
        }
    });
}

#[test]
fn uniform_field_refuses_target_volume_mismatch_and_sub_ulp_overcoverage() {
    let (coarse_mesh, _) = square_meshes();
    let system = FluxSystem::new(&coarse_mesh);
    let plan = system
        .prepare_density_transport(&vec![0.0; system.n_u], 2.0, 0.0, 1 << 20, || false)
        .unwrap();
    let cells = [Cell {
        density_kg_m3: 7.0,
        specific: [3.0],
    }; 2];
    with_cx(&CancelGate::new_clock_free(), |cx| {
        assert!(matches!(
            plan.remap_from(
                &plan,
                &cells,
                &[
                    Overlap::new(0, 0, 0.75),
                    Overlap::new(0, 1, 0.25),
                    Overlap::new(1, 1, 1.0),
                ],
                controls(),
                cx
            ),
            Err(Error::TargetVolumeCoverage { target: 0, .. })
        ));
        // The original scalar result rounds 1 + tiny to 1. A compensated
        // local coverage audit must still refuse a zero-tolerance mismatch.
        let tiny = 2.0_f64.powi(-60);
        match plan.remap_from(
            &plan,
            &cells,
            &[
                Overlap::new(0, 0, 1.0),
                Overlap::new(0, 1, tiny),
                Overlap::new(1, 1, 1.0),
            ],
            controls(),
            cx,
        ) {
            Err(Error::TargetVolumeCoverage {
                target, defect_m3, ..
            }) => {
                assert_eq!(target, 1);
                assert_eq!(defect_m3, tiny);
            }
            other => panic!("small target coverage defect was lost: {other:?}"),
        }
        assert!(matches!(
            plan.remap_from(
                &plan,
                &cells,
                &[Overlap::new(0, 0, 0.5), Overlap::new(1, 1, 1.0),],
                controls(),
                cx
            ),
            Err(Error::Remap {
                quantity: Quantity::Volume,
                error: ConservativeRemapError::SourceCoverageExceeded { .. }
            })
        ));
    });
}

fn triangle_meshes() -> (TriMesh, TriMesh) {
    // Two source cells of volume one coarsen into the same physical triangle
    // of volume two at depth one. All overlap fractions are exact ones.
    let nodes = [[0., 0.], [4., 0.], [0., 1.], [2., 0.]];
    (
        mesh(&nodes, &[[0, 3, 2], [3, 1, 2]]),
        mesh(&nodes[..3], &[[0, 1, 2]]),
    )
}

#[test]
fn tiny_mass_with_large_signed_specifics_retains_its_extensive_inventory() {
    let (source_mesh, target_mesh) = triangle_meshes();
    let source_system = FluxSystem::new(&source_mesh);
    let target_system = FluxSystem::new(&target_mesh);
    let source = source_system
        .prepare_density_transport(&vec![0.; source_system.n_u], 1., 0., 1 << 20, || false)
        .unwrap();
    let target = target_system
        .prepare_density_transport(&vec![0.; target_system.n_u], 1., 0., 1 << 20, || false)
        .unwrap();
    let cells = [
        Cell {
            density_kg_m3: 1e-200,
            specific: [1e200, -1e200],
        },
        Cell {
            density_kg_m3: 1e200,
            specific: [0.0, 0.0],
        },
    ];
    let overlaps = [Overlap::new(0, 0, 1.0), Overlap::new(1, 0, 1.0)];
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let mut allowance = controls();
        allowance.mass_tolerance_kg = 1e-199;
        allowance.carried_tolerance = [1e-14; 2];
        let outcome = target
            .remap_from(&source, &cells, &overlaps, allowance, cx)
            .unwrap();
        let state = outcome.cells[0];
        let mass = state.density_kg_m3 * target.cell_volumes_m3()[0];
        for (component, expected) in [1.0, -1.0].into_iter().enumerate() {
            assert_ne!(state.specific[component], 0.0);
            near(mass * state.specific[component], expected, 1e-14);
            near(outcome.report.carried_after[component], expected, 1e-14);
        }
        assert_eq!(outcome.report.mass_before_kg, outcome.report.mass_after_kg);
        assert_eq!(outcome.report.mass_balance_residual_kg, -1e-200);
        // Equal rounded total masses must not conceal a lost small source mass.
        allowance.mass_tolerance_kg = 0.0;
        assert!(
            matches!(target.remap_from(&source, &cells, &overlaps, allowance, cx),
            Err(Error::Conservation { quantity: Quantity::Mass, residual, .. }) if residual == -1e-200)
        );
    });
}

#[test]
fn signed_inventory_audit_retains_losses_hidden_by_rounded_endpoint_totals() {
    let (source_mesh, target_mesh) = triangle_meshes();
    let source_system = FluxSystem::new(&source_mesh);
    let target_system = FluxSystem::new(&target_mesh);
    let source = source_system
        .prepare_density_transport(&vec![0.; source_system.n_u], 1., 0., 1 << 20, || false)
        .unwrap();
    let target = target_system
        .prepare_density_transport(&vec![0.; target_system.n_u], 1., 0., 1 << 20, || false)
        .unwrap();
    let overlaps = [Overlap::new(0, 0, 1.0), Overlap::new(1, 0, 1.0)];
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for sign in [1.0, -1.0] {
            let cells = [
                Cell {
                    density_kg_m3: 1.0,
                    specific: [sign * 1e200],
                },
                Cell {
                    density_kg_m3: 1.0,
                    specific: [sign],
                },
            ];
            let mut allowance = controls();
            allowance.carried_tolerance = [1.0];
            let outcome = target
                .remap_from(&source, &cells, &overlaps, allowance, cx)
                .unwrap();
            assert_eq!(outcome.report.carried_before, outcome.report.carried_after);
            assert_eq!(outcome.report.carried_balance_residual, [-sign]);
            allowance.carried_tolerance = [0.0];
            assert!(
                matches!(target.remap_from(&source, &cells, &overlaps, allowance, cx),
                Err(Error::Conservation { quantity: Quantity::Carried(0), residual, .. }) if residual == -sign)
            );
        }
    });
}

#[test]
fn carried_product_remainders_expose_a_large_common_datum_reconstruction_loss() {
    let mesh = mesh(&[[0., 0.], [3., 0.], [0., 1.]], &[[0, 1, 2]]);
    let system = FluxSystem::new(&mesh);
    let plan = system
        .prepare_density_transport(&vec![0.; system.n_u], 1., 0., 1 << 20, || false)
        .unwrap();
    // Exact mass = 3/2 and q = 2^53 - 2. Rounded m*q is one below the exact
    // product; dividing it by m rounds q down by one. The new rounded product
    // is unchanged, although the represented physical inventory fell by 3/2.
    let cells = [Cell {
        density_kg_m3: 1.0,
        specific: [9_007_199_254_740_990.0],
    }];
    let overlaps = [Overlap::new(0, 0, 1.0)];
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let mut allowance = controls();
        allowance.carried_tolerance = [2.0];
        let outcome = plan
            .remap_from(&plan, &cells, &overlaps, allowance, cx)
            .unwrap();
        assert_eq!(outcome.cells[0].specific[0], cells[0].specific[0] - 1.0);
        assert_eq!(outcome.report.carried_before, outcome.report.carried_after);
        assert_eq!(outcome.report.carried_balance_residual, [-1.5]);
        allowance.carried_tolerance = [1.0];
        assert!(matches!(
            plan.remap_from(&plan, &cells, &overlaps, allowance, cx),
            Err(Error::Conservation {
                quantity: Quantity::Carried(0),
                residual,
                ..
            }) if residual == -1.5
        ));
    });
}

#[test]
fn cancellation_storage_and_invalid_state_publish_no_candidate() {
    let (coarse_mesh, fine_mesh) = square_meshes();
    let coarse_system = FluxSystem::new(&coarse_mesh);
    let fine_system = FluxSystem::new(&fine_mesh);
    let source = coarse_system
        .prepare_density_transport(&vec![0.; coarse_system.n_u], 2., 0., 1 << 20, || false)
        .unwrap();
    let target = fine_system
        .prepare_density_transport(&vec![0.; fine_system.n_u], 2., 0., 1 << 20, || false)
        .unwrap();
    let cells = [Cell {
        density_kg_m3: 2.0,
        specific: [4.0],
    }; 2];
    let saved = cells;
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_cx(&gate, |cx| {
        assert!(matches!(
            target.remap_from(&source, &cells, &refine_map(), controls(), cx),
            Err(Error::Cancelled)
        ));
    });
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let mut allowance = controls();
        let required = 2 * size_of::<f64>() + 4 * (2 * size_of::<f64>() + size_of::<Cell<1>>());
        allowance.maximum_storage_bytes = required - 1;
        assert!(
            matches!(target.remap_from(&source, &cells, &refine_map(), allowance, cx),
            Err(Error::Transport(DensityTransportError::Storage { required: actual, .. })) if actual == required)
        );
        allowance.maximum_storage_bytes = required;
        assert!(
            target
                .remap_from(&source, &cells, &refine_map(), allowance, cx)
                .is_ok()
        );
        let mut invalid = cells;
        invalid[1].density_kg_m3 = 0.0;
        assert!(matches!(
            target.remap_from(&source, &invalid, &refine_map(), controls(), cx),
            Err(Error::Transport(DensityTransportError::Invalid(_)))
        ));
        let mass_only = [Cell::<0> {
            density_kg_m3: 2.0,
            specific: [],
        }; 2];
        let outcome = target
            .remap_from(&source, &mass_only, &refine_map(), controls(), cx)
            .unwrap();
        assert_eq!(outcome.cells.len(), 4);
        assert_eq!(outcome.report.mass_after_kg, 4.0);
    });
    assert_eq!(cells, saved);
}
