use super::*;
use crate::d3q19::{BoundarySpec3, E3, Q3, equilibrium3};
use crate::sgs3::{SgsModel, SgsUse};

fn fixture() -> BoundaryGrid3 {
    let mut grid = BoundaryGrid3::with_collision_model(
        8,
        8,
        8,
        CollisionModel3::CentralMoment {
            second_order_rate: 1.3,
            higher_order_rate: 1.6,
        },
        [0.0; 3],
        BoundarySpec3::periodic(),
    );
    grid.voxelize_sdf(|p| {
        if (3.0..5.0).contains(&p[0]) && (3.0..5.0).contains(&p[1]) {
            -1.0
        } else {
            1.0
        }
    });
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if !grid.is_solid(x, y, z) {
                    let u = [
                        0.012 + 0.001 * y as f64,
                        0.003 * (x as f64).sin(),
                        0.005 * (z as f64).cos(),
                    ];
                    grid.set_populations(x, y, z, &equilibrium3(1.0, u));
                }
            }
        }
    }
    grid
}

fn probe(grid: &BoundaryGrid3) -> WallPatchProbe3 {
    let links = grid.boundary_links();
    let ids: Vec<_> = links.iter().map(|link| link.cell[2] / 4).collect();
    grid.wall_patch_probe(&links, &ids, &[[0.0; 3]; 2], 1_000_000)
        .unwrap()
}

fn bits(grid: &BoundaryGrid3) -> Vec<[u64; Q3]> {
    let mut result = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if !grid.is_solid(x, y, z) {
                    result.push(grid.populations(x, y, z).map(f64::to_bits));
                }
            }
        }
    }
    result
}

fn momentum(grid: &BoundaryGrid3) -> [f64; 3] {
    let mut sum = [0.0; 3];
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if !grid.is_solid(x, y, z) {
                    for (f, e) in grid.populations(x, y, z).into_iter().zip(E3) {
                        for (a, c) in [e.0, e.1, e.2].into_iter().enumerate() {
                            sum[a] += f * f64::from(c);
                        }
                    }
                }
            }
        }
    }
    sum
}

fn closure(coefficient: f64) -> SgsConfig {
    SgsConfig {
        model: SgsModel::Wale { cw: coefficient },
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

#[test]
fn g1_patch_forces_close_all_three_components_and_sum_to_native_loads() {
    let gate = CancelGate::new_clock_free();
    for sgs in [None, Some(closure(0.4))] {
        let mut distributed = fixture();
        let mut aggregate = fixture();
        let mut patches = probe(&distributed);
        let links = aggregate.boundary_links();
        let mut a_work = aggregate.sgs_workspace().unwrap();
        let mut d_work = distributed.sgs_workspace().unwrap();
        let before = momentum(&distributed);
        let total = if let Some(config) = sgs {
            distributed
                .step_with_sgs_and_wall_patches(&gate, config, &mut d_work, &mut patches)
                .unwrap();
            aggregate
                .step_with_sgs(&gate, config, &mut a_work, Some((&links, [0.0; 3])))
                .unwrap()
                .wall_momentum
                .unwrap()
        } else {
            distributed
                .step_with_wall_patches(&gate, &mut patches)
                .unwrap();
            aggregate.step_with_wall_momentum(&links, [0.0; 3])
        };
        assert_eq!(bits(&distributed), bits(&aggregate));
        let after = momentum(&distributed);
        let exchanges = patches.exchanges().unwrap();
        assert_eq!(
            exchanges.iter().map(|x| x.measured_links).sum::<usize>(),
            links.len()
        );
        assert!(
            exchanges
                .iter()
                .any(|x| x.fluid_population_impulse[2].abs() > 1e-5)
        );
        for a in 0..3 {
            let fluid: f64 = exchanges
                .iter()
                .map(|x| x.fluid_population_impulse[a])
                .sum();
            let wall: f64 = exchanges.iter().map(|x| x.wall_impulse[a]).sum();
            let moment: f64 = exchanges.iter().map(|x| x.wall_angular_impulse[a]).sum();
            assert!((fluid - (after[a] - before[a])).abs() < 3e-12);
            assert!((wall - total.wall_impulse[a]).abs() < 1e-12);
            assert!((moment - total.wall_angular_impulse[a]).abs() < 1e-11);
            assert!((wall + fluid).abs() < 1e-14);
        }
    }
}

#[test]
fn g3_zero_sgs_retains_patch_bits_and_g4_late_cancel_keeps_last_sample() {
    let gate = CancelGate::new_clock_free();
    let mut plain = fixture();
    let mut zero = fixture();
    let mut p = probe(&plain);
    let mut z = probe(&zero);
    let mut workspace = zero.sgs_workspace().unwrap();
    for _ in 0..3 {
        plain.step_with_wall_patches(&gate, &mut p).unwrap();
        zero.step_with_sgs_and_wall_patches(&gate, closure(0.0), &mut workspace, &mut z)
            .unwrap();
        assert_eq!(p.exchanges(), z.exchanges());
        assert_eq!(bits(&plain), bits(&zero));
    }
    let accepted = z.exchanges().unwrap().to_vec();
    let populations = bits(&zero);
    zero.prepare_sgs(&gate, closure(0.4), &mut workspace)
        .unwrap();
    let result = zero.step_sgs_prepared_measured(&gate, &workspace, |grid| {
        z.stage(grid, &gate)?;
        gate.request(); // After measured post-collision loads, before publication.
        Ok(())
    });
    assert!(matches!(
        result,
        Err(BoundarySgsError3::Boundary(BoundaryStepError3::Cancelled))
    ));
    assert_eq!(z.exchanges().unwrap(), accepted);
    assert_eq!(bits(&zero), populations);
    let retry = CancelGate::new_clock_free();
    zero.step_with_sgs_and_wall_patches(&retry, closure(0.4), &mut workspace, &mut z)
        .unwrap();
    assert_ne!(bits(&zero), populations);
}

#[test]
fn g4_duplicate_empty_partition_and_memory_refuse_before_stepping() {
    let grid = fixture();
    let links = grid.boundary_links();
    let ids = vec![0; links.len()];
    assert!(
        grid.wall_patch_probe(&links, &ids, &[[0.0; 3]; 2], 1_000_000)
            .is_err()
    );
    assert!(grid.wall_patch_probe(&links, &ids, &[[0.0; 3]], 0).is_err());
    assert!(
        grid.wall_patch_probe(&[links[0], links[0]], &[0, 0], &[[0.0; 3]], 1_000_000)
            .is_err()
    );
    let good = grid
        .wall_patch_probe(&links, &ids, &[[0.0; 3]], 1_000_000)
        .unwrap();
    assert!(
        grid.wall_patch_probe(&links, &ids, &[[0.0; 3]], good.workspace_bytes() - 1)
            .is_err()
    );
    assert!(good.exchanges().is_none());
}
