use super::*;
use crate::d3q19::force_corrected_stress3;
use crate::d3q19::freesurface3::{ContactModel3, PlanarPiston3, PlanarWall3};

fn controls() -> ImpulseControls {
    ImpulseControls {
        max_cells: 128,
        maximum_speed: 0.12,
        maximum_momentum_defect: 1e-12,
        maximum_work_defect: 1e-12,
    }
}
fn surface() -> FreeSurface3 {
    FreeSurface3::new(
        4,
        4,
        4,
        0.9,
        [1e-5, -2e-5, 3e-5],
        0.0,
        ContactModel3::Neutral,
        |x, _, _| x < 3,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn compact_readback_keeps_partial_inventory_guo_frame_and_stress() {
    let mut liquid = surface();
    let ids = [1, 2];
    let gate = CancelGate::new_clock_free();
    for x in [1, 2] {
        let (slot, lane) = liquid.grid.resolve_source(x, 0, 0).unwrap();
        let mut f = liquid.grid.populations(slot, lane);
        let positive = E3.iter().position(|&e| e == (1, 0, 0)).unwrap();
        let negative = E3.iter().position(|&e| e == (-1, 0, 0)).unwrap();
        f[positive] += 2e-5;
        f[negative] += 2e-5;
        f[0] -= 4e-5;
        liquid.grid.set_populations(slot, lane, f);
    }
    let before_bits = liquid.grid.state_bits();
    let ledger = liquid.ledger_mass();
    let epoch = liquid.geometry_epoch();
    let stage = liquid
        .prepare_impulse_stencil(&ids, controls(), &gate)
        .unwrap();
    let samples: Vec<_> = stage.cells().collect();
    near(samples[0].mass, 1.0, 3e-16);
    near(samples[1].mass, 0.5, 3e-16);
    let old: Vec<_> = stage
        .cells
        .iter()
        .map(|c| {
            let density = c.populations.iter().sum::<f64>();
            (
                density,
                force_corrected_stress3(&c.populations, density, c.sample.velocity, c.force),
                c.force,
            )
        })
        .collect();
    let impulses = [[0.003, -0.001, 0.002], [-0.002, 0.001, -0.003]];
    let candidate = stage.with_impulses(&impulses, &gate).unwrap();
    assert_eq!(candidate.owner.grid.state_bits(), before_bits);
    assert_eq!(candidate.report().affected_cells, 2);
    for (i, (read, cell)) in candidate.readback().zip(&candidate.cells).enumerate() {
        for axis in 0..3 {
            near(
                read.velocity_after[axis] - read.before.velocity[axis],
                impulses[i][axis] / read.before.mass,
                1e-16,
            );
            near(read.momentum_change[axis], impulses[i][axis], 1e-16);
        }
        assert_eq!(cell.force, old[i].2);
        let rho = cell.populations.iter().sum::<f64>();
        near(rho, old[i].0, 3e-16);
        let stress =
            force_corrected_stress3(&cell.populations, rho, read.velocity_after, cell.force);
        for (a, b) in stress.into_iter().zip(old[i].1) {
            near(a, b, 2e-16);
        }
    }
    let report = candidate.commit(&gate).unwrap();
    assert!(report.absolute_momentum_defect < 1e-15);
    assert_eq!(liquid.grid.steps(), 0);
    assert_eq!(liquid.geometry_epoch(), epoch);
    near(liquid.ledger_mass(), ledger, 3e-15);
    let after = liquid
        .prepare_impulse_stencil(&ids, controls(), &gate)
        .unwrap();
    for (a, b) in after.cells().zip(samples) {
        assert_eq!(a.mass.to_bits(), b.mass.to_bits());
    }
}

#[test]
fn dropping_and_final_cancellation_leave_original_population_bits() {
    let mut liquid = surface();
    let original = liquid.grid.state_bits();
    let gate = CancelGate::new_clock_free();
    {
        let stage = liquid
            .prepare_impulse_stencil(&[1, 2], controls(), &gate)
            .unwrap();
        let _candidate = stage.with_impulses(&[[0.001, 0.0, 0.0]; 2], &gate).unwrap();
    }
    assert_eq!(liquid.grid.state_bits(), original);
    let stage = liquid
        .prepare_impulse_stencil(&[1, 2], controls(), &gate)
        .unwrap();
    let candidate = stage.with_impulses(&[[0.001, 0.0, 0.0]; 2], &gate).unwrap();
    gate.request();
    assert!(matches!(
        candidate.commit(&gate),
        Err(ImpulseError::Cancelled)
    ));
    assert_eq!(liquid.grid.state_bits(), original);
    let gate = CancelGate::new_clock_free();
    let candidate = liquid
        .prepare_impulse_stencil(&[1, 2], controls(), &gate)
        .unwrap()
        .with_impulses(&[[0.0; 3]; 2], &gate)
        .unwrap();
    assert_eq!(candidate.commit(&gate).unwrap().affected_cells, 0);
    assert_eq!(liquid.grid.state_bits(), original);
}

#[test]
fn selection_rejects_gas_solids_pistons_duplicates_and_unlocated_inventory() {
    let gate = CancelGate::new_clock_free();
    let mut liquid = surface();
    let bits = liquid.grid.state_bits();
    for ids in [&[1, 1][..], &[2, 1][..], &[3][..], &[64][..], &[][..]] {
        assert!(
            liquid
                .prepare_impulse_stencil(ids, controls(), &gate)
                .is_err()
        );
        assert_eq!(liquid.grid.state_bits(), bits);
    }
    liquid.planar_piston = Some(PlanarPiston3 {
        wall: PlanarWall3::XMin,
        depth_cells: 2,
    });
    assert!(
        liquid
            .prepare_impulse_stencil(&[1], controls(), &gate)
            .is_err()
    );
    liquid.planar_piston = None;
    let (slot, lane) = liquid.grid.resolve_source(2, 0, 0).unwrap();
    let key = liquid.grid.active_keys()[slot];
    for mass in [0.0, -0.1, 1.2, f64::NAN] {
        liquid.tiles.get_mut(&key).unwrap().mass[lane] = mass;
        assert!(
            liquid
                .prepare_impulse_stencil(&[2], controls(), &gate)
                .is_err()
        );
        assert_eq!(liquid.grid.state_bits(), bits);
    }
    let mut solid = FreeSurface3::with_solid_tiles(
        8,
        4,
        4,
        0.9,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        &[(1, 0, 0)],
        |x, _, _| x < 3,
    )
    .unwrap();
    let bits = solid.grid.state_bits();
    assert!(solid.grid.resolve_source(4, 0, 0).is_none());
    assert!(
        solid
            .prepare_impulse_stencil(&[4], controls(), &gate)
            .is_err()
    );
    assert_eq!(solid.grid.state_bits(), bits);
}
