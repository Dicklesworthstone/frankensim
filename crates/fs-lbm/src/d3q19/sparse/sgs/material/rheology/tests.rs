//! Synthetic source/flow invariants, not calibrated rheology or turbulence
//! validation. The test source owns its constitutive response, as callers do.

use super::*;
use crate::d3q19::equilibrium3;
use crate::d3q19::sparse::sgs::tests::{ACTIVE, seeded, settings};
use crate::material::plan_scaling_from_material;
use crate::sgs3::{SgsError, SgsModel};
use fs_exec::{RunError, RunReport, TileKernel, TilePool};
use fs_material::fluid_transport::FluidTransportState;
use std::cell::Cell;
use std::sync::atomic::{AtomicUsize, Ordering};

struct Source {
    state: Result<FluidTransportState, MaterialError>,
    reads: Cell<usize>,
}
impl FluidTransport for Source {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.reads.set(self.reads.get() + 1);
        self.state.clone()
    }
}

fn fluid(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 1000.0, 1.0, mu).unwrap()
}

fn epoch<M>(
    w: &SparseRheologySgsWorkspace3<M>,
) -> Vec<(MaterialImmersedPoint3, FluidTransportState, f64)> {
    w.applied()
        .iter()
        .map(|s| (s.point, s.transport, s.shear_relaxation_time))
        .collect()
}

#[test]
fn live_constant_sources_match_frozen_molecular_fields_including_zero_sgs() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, 0.01, 0.01, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    let count = ACTIVE.len() * TILE_CELLS;
    for model in [
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Smagorinsky { cs: 0.18 },
        SgsModel::Vreman { cv: 0.07 },
    ] {
        let low = Source {
            state: Ok(fluid(0.4)),
            reads: Cell::new(0),
        };
        let high = Source {
            state: Ok(fluid(1.6)),
            reads: Cell::new(0),
        };
        let source = |cell: [usize; 3]| if cell[0] < 4 { &low } else { &high };
        let mut live = seeded(&ACTIVE);
        let mut frozen = live.clone();
        let bound = scaling
            .bind_sparse_sgs_materials(&live, count, &gate, |p| Ok(source(p)))
            .unwrap();
        let mut fw = frozen.sgs_workspace(&gate).unwrap();
        let mut work = scaling
            .sparse_rheology_sgs_workspace(&live, count, &gate)
            .unwrap();
        let capacity = work.workspace_bytes();
        let collisions = work.flow.collisions().as_ptr();
        for _ in 0..3 {
            let expected = bound
                .step_serial(&mut frozen, &gate, settings(model), &mut fw)
                .unwrap();
            let report = scaling
                .step_sparse_with_material_rheology_and_sgs(
                    &mut live,
                    &gate,
                    settings(model),
                    &mut work,
                    |p| Ok(source(p.cell)),
                )
                .unwrap();
            assert_eq!(live.state_bits(), frozen.state_bits());
            assert_eq!(report, expected);
            assert_eq!(work.last_report(), Some(report));
            assert_eq!(work.workspace_bytes(), capacity);
            assert_eq!(work.flow.collisions().as_ptr(), collisions);
            for (applied, held) in work.applied().iter().zip(bound.samples()) {
                assert!(core::ptr::eq(applied.material, source(applied.point.cell)));
                assert_eq!(applied.transport, held.transport);
                assert_eq!(
                    applied.shear_relaxation_time,
                    held.molecular_relaxation_time
                );
            }
        }
        assert_eq!(low.reads.get() + high.reads.get(), 4 * count);
    }
}

#[test]
fn resolved_shear_uses_physical_units_half_force_walls_and_live_source_feedback() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, 0.01, 0.01, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    let mut grid = SparseGrid3::new(4, 8, 4, 0.8, [2e-5, -1e-5, 0.0]).unwrap();
    grid.activate_tiles(&[(0, 0, 0), (0, 1, 0)]).unwrap();
    for slot in 0..grid.keys.len() {
        let (_, ty, _) = demorton3(grid.keys[slot]);
        for lane in 0..TILE_CELLS {
            let y = (ty as usize * TILE + lane / TILE % TILE) as f64 + 0.5;
            let u = [1e-3 * y * (8.0 - y), 0.0, 0.0];
            grid.set_populations(
                slot,
                lane,
                equilibrium3(1.0, core::array::from_fn(|a| u[a] - 0.5 * grid.force[a])),
            );
        }
    }
    let config = settings(SgsModel::Smagorinsky { cs: 0.18 });
    let mut work = scaling
        .sparse_rheology_sgs_workspace(&grid, 128, &gate)
        .unwrap();
    let resolve = |p: MaterialImmersedPoint3| Ok(fluid(1.0 / (1.0 + p.shear_rate_per_s)));
    scaling
        .step_sparse_with_material_rheology_and_sgs(&mut grid, &gate, config, &mut work, resolve)
        .unwrap();
    for y in 0..8 {
        let (i, s) = work
            .applied()
            .iter()
            .enumerate()
            .find(|(_, s)| s.point.cell == [1, y, 1])
            .unwrap();
        let derivative = (1e-3 * (8.0 - 2.0 * (y as f64 + 0.5))).abs();
        let rate = derivative / 0.01;
        assert!((s.point.shear_rate_per_s - rate).abs() < 2e-14);
        let mu = 1.0 / (1.0 + rate);
        assert!((s.transport.dynamic_viscosity_pa_s() - mu).abs() < 2e-14);
        let molecular_nu = mu / 1000.0 * 0.01 / (0.01 * 0.01);
        let total = work.flow.collisions()[i].kinematic_viscosity();
        assert!((total - molecular_nu - 0.18 * 0.18 * derivative).abs() < 3e-15);
    }
    let previous = epoch(&work);
    let mut held_flow = grid.clone();
    let held = scaling
        .bind_sparse_sgs_materials(&held_flow, 128, &gate, |cell| {
            Ok(previous.iter().find(|(p, _, _)| p.cell == cell).unwrap().1)
        })
        .unwrap();
    let mut held_work = held_flow.sgs_workspace(&gate).unwrap();
    held.step_serial(&mut held_flow, &gate, config, &mut held_work)
        .unwrap();
    scaling
        .step_sparse_with_material_rheology_and_sgs(&mut grid, &gate, config, &mut work, resolve)
        .unwrap();
    assert_ne!(
        epoch(&work),
        previous,
        "source was not reevaluated at the new strain"
    );
    assert_ne!(
        grid.state_bits(),
        held_flow.state_bits(),
        "source feedback did not reach collision"
    );

    // Native full-gradient queries cancel rigid rotation and retain diagonal
    // strain. Check an interior cell so stationary walls do not add strain.
    for dilation in [0.0, 0.002] {
        let mut spin = seeded(&[(0, 0, 0)]);
        for lane in 0..TILE_CELLS {
            let [x, y, z] =
                [lane % TILE, lane / TILE % TILE, lane / (TILE * TILE)].map(|p| p as f64 + 0.5);
            let u = [
                dilation * x - 0.003 * y,
                0.003 * x + dilation * y,
                dilation * z,
            ];
            spin.set_populations(
                0,
                lane,
                equilibrium3(1.0, core::array::from_fn(|a| u[a] - 0.5 * spin.force[a])),
            );
        }
        let mut sw = scaling
            .sparse_rheology_sgs_workspace(&spin, 64, &gate)
            .unwrap();
        scaling
            .step_sparse_with_material_rheology_and_sgs(
                &mut spin,
                &gate,
                settings(SgsModel::Wale { cw: 0.0 }),
                &mut sw,
                |_| Ok(reference),
            )
            .unwrap();
        let s = sw
            .applied()
            .iter()
            .find(|s| s.point.cell == [1, 1, 1])
            .unwrap();
        assert!((s.point.shear_rate_per_s - 6.0_f64.sqrt() * dilation / 0.01).abs() < 2e-13);
    }
}

#[test]
fn source_ratio_and_cancellation_refusals_keep_the_applied_epoch() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, 0.01, 0.01, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    let good = Source {
        state: Ok(reference),
        reads: Cell::new(0),
    };
    let source_error = MaterialError::State {
        what: "original source domain refusal".into(),
    };
    let bad = Source {
        state: Err(source_error.clone()),
        reads: Cell::new(0),
    };
    let mut grid = seeded(&ACTIVE);
    let mut work = scaling
        .sparse_rheology_sgs_workspace(&grid, ACTIVE.len() * TILE_CELLS, &gate)
        .unwrap();
    let config = settings(SgsModel::Smagorinsky { cs: 0.18 });
    scaling
        .step_sparse_with_material_rheology_and_sgs(&mut grid, &gate, config, &mut work, |_| {
            Ok(&good)
        })
        .unwrap();
    let before = grid.state_bits();
    let accepted = epoch(&work);
    let report = work.last_report();
    let last = accepted.last().unwrap().0.cell;
    for transport_failure in [false, true] {
        let refusal = scaling.step_sparse_with_material_rheology_and_sgs(
            &mut grid,
            &gate,
            config,
            &mut work,
            |p| {
                if p.cell != last {
                    Ok(&good)
                } else if transport_failure {
                    Ok(&bad)
                } else {
                    Err(source_error.clone())
                }
            },
        );
        assert_eq!(
            refusal,
            Err(SparseMaterialSgsError3::Source {
                cell: last,
                source: MaterialLbmError::Transport(source_error.clone())
            })
        );
        assert_eq!(grid.state_bits(), before);
        assert_eq!(epoch(&work), accepted);
        assert_eq!(work.last_report(), report);
    }
    let mut limited = config;
    limited.maximum_eddy_to_molecular = 0.0;
    assert!(matches!(
        scaling.step_sparse_with_material_rheology_and_sgs(
            &mut grid,
            &gate,
            limited,
            &mut work,
            |_| Ok(&good)
        ),
        Err(SparseMaterialSgsError3::Sgs(SparseSgsError3::Closure(
            SgsError::ViscosityRatioExceeded { .. }
        )))
    ));
    let cancelled = CancelGate::new_clock_free();
    assert!(matches!(
        scaling.step_sparse_with_material_rheology_and_sgs(
            &mut grid,
            &cancelled,
            config,
            &mut work,
            |p| {
                if p.cell == last {
                    cancelled.request();
                }
                Ok(&good)
            }
        ),
        Err(SparseMaterialSgsError3::Sgs(SparseSgsError3::Flow(
            SparseError3::Cancelled
        )))
    ));
    assert_eq!(grid.state_bits(), before);
    assert_eq!(grid.steps(), 1);
    assert_eq!(epoch(&work), accepted);
    assert_eq!(work.last_report(), report);
    let other_scaling = plan_scaling_from_material(&reference, 0.01, 0.02, 0.05).unwrap();
    let reads = good.reads.get();
    assert!(matches!(
        other_scaling.step_sparse_with_material_rheology_and_sgs(
            &mut grid,
            &gate,
            config,
            &mut work,
            |_| Ok(&good)
        ),
        Err(SparseMaterialSgsError3::Invalid(_))
    ));
    grid.deactivate_tiles(&[(0, 0, 0)]).unwrap();
    assert!(matches!(
        scaling.step_sparse_with_material_rheology_and_sgs(
            &mut grid,
            &gate,
            config,
            &mut work,
            |_| Ok(&good)
        ),
        Err(SparseMaterialSgsError3::Sgs(
            SparseSgsError3::TopologyChanged
        ))
    ));
    assert_eq!(good.reads.get(), reads);
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
        let result = self.pool.run_with_gate(kernel, gate);
        if self.passes.fetch_add(1, Ordering::Relaxed) == 1 {
            gate.request();
        }
        result
    }
}

#[test]
fn live_rheology_pool_matches_serial_and_native_late_cancellation_is_atomic() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, 0.01, 0.01, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    let tiles: Vec<_> = (0..3)
        .flat_map(|x| (0..2).flat_map(move |y| (0..2).map(move |z| (x, y, z))))
        .filter(|p| *p != (2, 1, 1))
        .collect();
    let mut serial = seeded(&tiles);
    let mut pooled = serial.clone();
    let mut sw = scaling
        .sparse_rheology_sgs_workspace(&serial, tiles.len() * 64, &gate)
        .unwrap();
    let mut pw = scaling
        .sparse_rheology_sgs_workspace(&pooled, tiles.len() * 64, &gate)
        .unwrap();
    let pool = TilePool::for_host(3, 977);
    let config = settings(SgsModel::Wale { cw: 0.5 });
    let resolve = |p: MaterialImmersedPoint3| Ok(fluid(0.4 + 0.1 * p.shear_rate_per_s));
    for _ in 0..2 {
        let sr = scaling
            .step_sparse_with_material_rheology_and_sgs(
                &mut serial,
                &gate,
                config,
                &mut sw,
                resolve,
            )
            .unwrap();
        let pr = scaling
            .step_sparse_pooled_with_material_rheology_and_sgs(
                &mut pooled,
                &pool,
                &gate,
                config,
                &mut pw,
                resolve,
            )
            .unwrap();
        assert_eq!(sr, pr);
        assert_eq!(serial.state_bits(), pooled.state_bits());
        assert_eq!(epoch(&sw), epoch(&pw));
    }
    let before = pooled.state_bits();
    let accepted = epoch(&pw);
    let report = pw.last_report();
    let cancelled = CancelGate::new_clock_free();
    let stop = CancelAfterStream {
        pool,
        passes: AtomicUsize::new(0),
    };
    assert!(matches!(
        scaling.step_sparse_pooled_with_material_rheology_and_sgs(
            &mut pooled,
            &stop,
            &cancelled,
            config,
            &mut pw,
            resolve
        ),
        Err(SparseMaterialSgsError3::Sgs(SparseSgsError3::Flow(
            SparseError3::Cancelled
        )))
    ));
    assert_eq!(pooled.state_bits(), before);
    assert_eq!(pooled.steps(), 2);
    assert_eq!(epoch(&pw), accepted);
    assert_eq!(pw.last_report(), report);
    scaling
        .step_sparse_with_material_rheology_and_sgs(&mut serial, &gate, config, &mut sw, resolve)
        .unwrap();
    scaling
        .step_sparse_pooled_with_material_rheology_and_sgs(
            &mut pooled,
            &stop,
            &gate,
            config,
            &mut pw,
            resolve,
        )
        .unwrap();
    assert_eq!(serial.state_bits(), pooled.state_bits());
    assert_eq!(epoch(&sw), epoch(&pw));
}
