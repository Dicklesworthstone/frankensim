//! G0/G1/G4/G5 checks on the native sparse SGS trajectory, not a developed
//! turbulence benchmark or an acceptance-render qualification.

use super::*;
use crate::d3q19::equilibrium3;
use crate::lattice3::{Grid3, GridBudget3};
use crate::sgs3::{SgsModel, SgsUse};
use fs_exec::{Cx, RunError, RunReport, TileKernel, TilePool};
use std::sync::atomic::{AtomicUsize, Ordering};

const DIMS: [usize; 3] = [12, 8, 8];
const ACTIVE: [(u32, u32, u32); 5] = [(1, 0, 1), (0, 0, 0), (1, 0, 0), (1, 1, 0), (2, 0, 0)];

fn settings(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 10.0,
        usage: SgsUse::Exploration,
    }
}

fn coordinates(key: u64, lane: usize) -> [usize; 3] {
    let (x, y, z) = demorton3(key);
    [
        x as usize * TILE + lane % TILE,
        y as usize * TILE + (lane / TILE) % TILE,
        z as usize * TILE + lane / (TILE * TILE),
    ]
}

fn seeded(tiles: &[(u32, u32, u32)]) -> SparseGrid3 {
    let [nx, ny, nz] = DIMS;
    let mut grid = SparseGrid3::new(nx, ny, nz, 0.8, [1e-5, -2e-5, 3e-5]).unwrap();
    grid.activate_tiles(tiles).unwrap();
    for slot in 0..grid.keys.len() {
        for lane in 0..TILE_CELLS {
            let [x, y, z] = coordinates(grid.keys[slot], lane).map(|x| x as f64 + 0.5);
            let rho = 1.0 + 0.001 * (x + 2.0 * y + z);
            let u = [0.004 * (y - 4.0), -0.003 * (z - 4.0), 0.002 * (x - 6.0)];
            grid.set_populations(slot, lane, equilibrium3(rho, u));
        }
    }
    grid
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 713,
                kernel_id: 158,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

#[test]
fn sparse_sgs_matches_dense_fields_and_native_evolution() {
    let gate = CancelGate::new_clock_free();
    for model in [
        SgsModel::Smagorinsky { cs: 0.18 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.07 },
    ] {
        let mut sparse = seeded(&ACTIVE);
        let mut dense = Grid3::<D3q19, 19>::try_new(
            DIMS,
            sparse.model,
            sparse.force,
            [false; 3],
            GridBudget3 {
                maximum_cells: DIMS.iter().product(),
                maximum_mach: 0.2,
            },
            |p| {
                sparse
                    .resolve_source(p[0] as i64, p[1] as i64, p[2] as i64)
                    .is_none()
            },
            |p| {
                sparse
                    .resolve_source(p[0] as i64, p[1] as i64, p[2] as i64)
                    .map_or(equilibrium3(1.0, [0.0; 3]), |(slot, lane)| {
                        sparse.populations(slot, lane)
                    })
            },
        )
        .unwrap();
        let mut sw = sparse.sgs_workspace(&gate).unwrap();
        let mut dw = dense.sgs_workspace().unwrap();
        with_cx(|cx| {
            for _ in 0..3 {
                let sr = sparse.prepare_sgs(&gate, settings(model), &mut sw).unwrap();
                let dr = dense
                    .prepare_sgs(cx, settings(model), None, &mut dw)
                    .unwrap();
                assert_eq!(sr.fluid_cells, ACTIVE.len() * TILE_CELLS);
                assert_eq!(sr.fluid_cells, dr.fluid_cells);
                assert!(sr.maximum_eddy_viscosity_lu > 0.0);
                assert_eq!(sr.maximum_eddy_viscosity_lu, dr.maximum_eddy_viscosity_lu);
                assert!(
                    (sr.sgs_dissipation_lu - dr.sgs_dissipation_lu).abs()
                        < 1e-12 * sr.sgs_dissipation_lu
                );
                for (slot, &key) in sparse.keys.iter().enumerate() {
                    for lane in 0..TILE_CELLS {
                        let [x, y, z] = coordinates(key, lane);
                        let cell = (z * DIMS[1] + y) * DIMS[0] + x;
                        assert_eq!(
                            sw.collisions()[slot * TILE_CELLS + lane],
                            dw.collisions()[cell]
                        );
                    }
                }
                sparse
                    .step_serial_with_sgs(&gate, settings(model), &mut sw)
                    .unwrap();
                dense
                    .step_serial_with_sgs(cx, settings(model), None, &mut dw)
                    .unwrap();
                for (slot, &key) in sparse.keys.iter().enumerate() {
                    for lane in 0..TILE_CELLS {
                        let [x, y, z] = coordinates(key, lane);
                        let cell = (z * DIMS[1] + y) * DIMS[0] + x;
                        assert_eq!(
                            sparse.populations(slot, lane).map(f64::to_bits),
                            dense.populations()[cell].map(f64::to_bits)
                        );
                    }
                }
            }
        });
    }
}

#[test]
fn sparse_halfway_wall_and_tile_halo_recover_quadratic_shear() {
    let gate = CancelGate::new_clock_free();
    let mut grid = SparseGrid3::new(4, 8, 4, 0.8, [0.0; 3]).unwrap();
    grid.activate_tiles(&[(0, 0, 0), (0, 1, 0)]).unwrap();
    for slot in 0..grid.keys.len() {
        for lane in 0..TILE_CELLS {
            let y = coordinates(grid.keys[slot], lane)[1] as f64 + 0.5;
            grid.set_populations(
                slot,
                lane,
                equilibrium3(1.0, [1e-4 * y * (8.0 - y), 0.0, 0.0]),
            );
        }
    }
    let mut workspace = grid.sgs_workspace(&gate).unwrap();
    grid.prepare_sgs(
        &gate,
        settings(SgsModel::Smagorinsky { cs: 0.2 }),
        &mut workspace,
    )
    .unwrap();
    for y in 0..8 {
        let (slot, lane) = grid.resolve_source(1, y, 1).unwrap();
        let derivative = 1e-4 * (8.0 - 2.0 * (y as f64 + 0.5));
        assert!((workspace.gradient(slot, lane)[0][1] - derivative).abs() < 2e-16);
        let eddy = workspace.collisions()[slot * TILE_CELLS + lane].kinematic_viscosity()
            - grid.model.kinematic_viscosity();
        assert!((eddy - 0.04 * derivative.abs()).abs() < 2e-16);
    }
}

#[test]
fn sparse_zero_control_preserves_collision_families_and_workspace() {
    let gate = CancelGate::new_clock_free();
    for law in [
        CollisionModel3::Bgk { tau: 0.8 },
        CollisionModel3::Trt {
            tau: 0.8,
            magic: 0.19,
        },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.3,
            higher_order_rate: 1.1,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.3,
            third_order_rate: 1.2,
            fourth_order_rate: 1.1,
        },
    ] {
        let mut plain = seeded(&ACTIVE);
        plain.set_collision_model(law).unwrap();
        let mut sgs = plain.clone();
        let mut workspace = sgs.sgs_workspace(&gate).unwrap();
        let pointer = workspace.collisions().as_ptr();
        let bytes = workspace.workspace_bytes();
        for model in [
            SgsModel::Smagorinsky { cs: 0.0 },
            SgsModel::Wale { cw: 0.0 },
            SgsModel::Vreman { cv: 0.0 },
        ] {
            plain.step_serial().unwrap();
            let report = sgs
                .step_serial_with_sgs(&gate, settings(model), &mut workspace)
                .unwrap();
            assert_eq!(report.sgs_dissipation_lu, 0.0);
            assert_eq!(sgs.state_bits(), plain.state_bits());
            assert_eq!(workspace.collisions().as_ptr(), pointer);
            assert_eq!(workspace.workspace_bytes(), bytes);
            assert!(workspace.collisions().iter().all(|c| *c == law));
        }
    }
}

#[test]
fn sparse_workspace_tracks_active_topology_and_not_enclosing_volume() {
    let gate = CancelGate::new_clock_free();
    let mut grid = seeded(&ACTIVE);
    let mut workspace = grid.sgs_workspace(&gate).unwrap();
    let mut huge = SparseGrid3::new(4096, 4096, 4096, 0.8, [0.0; 3]).unwrap();
    huge.activate_tiles(&ACTIVE).unwrap();
    assert_eq!(
        huge.sgs_workspace(&gate).unwrap().workspace_bytes(),
        workspace.workspace_bytes()
    );
    assert!(matches!(
        huge.prepare_sgs(&gate, settings(SgsModel::Wale { cw: 0.5 }), &mut workspace),
        Err(SparseSgsError3::TopologyChanged)
    ));
    grid.deactivate_tiles(&[(1, 0, 1)]).unwrap();
    grid.activate_tiles(&[(2, 1, 1)]).unwrap();
    assert_eq!(grid.keys.len(), workspace.keys.len());
    let before = grid.state_bits();
    assert!(matches!(
        grid.step_serial_with_sgs(&gate, settings(SgsModel::Wale { cw: 0.5 }), &mut workspace),
        Err(SparseSgsError3::TopologyChanged)
    ));
    assert_eq!(grid.state_bits(), before);
    let mut rebound = grid.sgs_workspace(&gate).unwrap();
    grid.step_serial_with_sgs(&gate, settings(SgsModel::Wale { cw: 0.5 }), &mut rebound)
        .unwrap();
    let empty = SparseGrid3::new(4, 4, 4, 0.8, [0.0; 3]).unwrap();
    assert!(matches!(
        empty.sgs_workspace(&gate),
        Err(SparseSgsError3::EmptyGrid)
    ));
}

struct CancelAfterStream {
    pool: TilePool,
    passes: AtomicUsize,
}

impl KernelRunner for CancelAfterStream {
    fn workers(&self) -> usize {
        self.pool.workers()
    }
    fn run_with_gate<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> (Result<K::Out, RunError>, RunReport) {
        let output = self.pool.run_with_gate(kernel, gate);
        if self.passes.fetch_add(1, Ordering::Relaxed) == 1 {
            gate.request();
        }
        output
    }
}

#[test]
fn sparse_active_sgs_pool_is_deterministic_and_late_refusal_is_atomic() {
    let tiles: Vec<_> = (0..3)
        .flat_map(|x| (0..2).flat_map(move |y| (0..2).map(move |z| (x, y, z))))
        .filter(|p| *p != (2, 1, 1))
        .collect();
    let settings = settings(SgsModel::Wale { cw: 0.5 });
    let initial = seeded(&tiles);
    let gate = CancelGate::new_clock_free();
    let mut serial = initial.clone();
    let mut workspace = serial.sgs_workspace(&gate).unwrap();
    let mass = serial.total_mass();
    for _ in 0..3 {
        serial
            .step_serial_with_sgs(&gate, settings, &mut workspace)
            .unwrap();
    }
    // The sparse owner sums all populations directly; retain its existing
    // closed-box roundoff envelope rather than requiring a reduction redesign.
    assert!((serial.total_mass() / mass - 1.0).abs() < 1e-12);
    for workers in [1, 2, 4] {
        let pool = TilePool::for_host(workers, 713);
        let mut grid = initial.clone();
        let mut workspace = grid.sgs_workspace(&gate).unwrap();
        for _ in 0..3 {
            grid.step_pooled_with_sgs(&pool, &gate, settings, &mut workspace)
                .unwrap();
        }
        assert_eq!(grid.state_bits(), serial.state_bits());
    }
    let pool = CancelAfterStream {
        pool: TilePool::for_host(2, 713),
        passes: AtomicUsize::new(0),
    };
    let mut grid = initial.clone();
    let mut workspace = grid.sgs_workspace(&gate).unwrap();
    let mut refused = settings;
    refused.maximum_eddy_to_molecular = 0.0;
    assert!(matches!(
        grid.step_pooled_with_sgs(&pool, &gate, refused, &mut workspace),
        Err(SparseSgsError3::Closure(
            SgsError::ViscosityRatioExceeded { .. }
        ))
    ));
    assert_eq!(grid.state_bits(), initial.state_bits());
    assert!(matches!(
        grid.step_pooled_with_sgs(&pool, &gate, settings, &mut workspace),
        Err(SparseSgsError3::Flow(SparseError3::Cancelled))
    ));
    assert_eq!(grid.state_bits(), initial.state_bits());
    assert_eq!(grid.steps(), initial.steps());
    let retry = CancelGate::new_clock_free();
    for _ in 0..3 {
        grid.step_pooled_with_sgs(&pool, &retry, settings, &mut workspace)
            .unwrap();
    }
    assert_eq!(grid.state_bits(), serial.state_bits());
}
