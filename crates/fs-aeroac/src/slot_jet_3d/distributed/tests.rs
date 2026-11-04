use super::*;

fn config() -> SlotJet3dConfig {
    let mut cfg = ladder_config(1.3);
    cfg.nx = 64;
    cfg.ny = 32;
    cfg.nz = 8;
    cfg.plate_length = 3;
    cfg.fringe_width = 8;
    cfg.steps_settle = 1;
    cfg.steps_record = 64;
    cfg
}

fn budget() -> SlotJet3dPatchBudget {
    SlotJet3dPatchBudget {
        maximum_patches: 32,
        maximum_record_bytes: 1_000_000,
        maximum_probe_bytes: 1_000_000,
    }
}

#[test]
fn g1_actual_plate_links_have_bounded_geometry_and_keep_spanwise_loads() {
    let cfg = config();
    let mut distributed = Rig3::build(&cfg);
    let mut aggregate = Rig3::build(&cfg);
    let layout = SlotJet3dPatchLayout {
        cells_x: 2,
        cells_z: 3,
    };
    let (count, _) = admit(&cfg, layout, budget()).unwrap();
    let (patches, mut probe) = patch_probe(
        &cfg,
        &distributed,
        layout,
        count,
        budget().maximum_probe_bytes,
    )
    .unwrap();
    assert_eq!(count, 6);
    assert_eq!(
        patches.iter().map(|p| p.measured_links).sum::<usize>(),
        distributed.plate_links.len()
    );
    assert!(
        patches
            .iter()
            .all(|p| p.radius_lu.is_finite() && p.radius_lu > 0.0 && p.radius_lu < 3.0)
    );
    // Every actual link point, including a diagonal crossing z=0, belongs to
    // its solid-cell patch without creating a box-spanning bounding sphere.
    for link in &distributed.plate_links {
        let e = E3[link.direction];
        let e = [e.0, e.1, e.2];
        let dims = [cfg.nx, cfg.ny, cfg.nz];
        let solid: [usize; 3] = core::array::from_fn(|a| {
            (link.cell[a] as isize - e[a] as isize).rem_euclid(dims[a] as isize) as usize
        });
        let p = solid[2] / 3 * 2 + (solid[0] - cfg.nozzle_thickness - cfg.edge_distance) / 2;
        let distance = det::sqrt(
            (0..3)
                .map(|a| {
                    let d = solid[a] as f64 + 0.5 + 0.5 * f64::from(e[a]) - patches[p].center_lu[a];
                    d * d
                })
                .sum(),
        );
        assert!(distance <= patches[p].radius_lu + 1e-14);
    }
    distributed
        .grid
        .step_with_wall_patches(&CancelGate::new_clock_free(), &mut probe)
        .unwrap();
    let total = aggregate
        .grid
        .step_with_wall_momentum(&aggregate.plate_links, [0.0; 3]);
    for a in 0..3 {
        let resolved: f64 = probe
            .exchanges()
            .unwrap()
            .iter()
            .map(|p| p.fluid_population_impulse[a])
            .sum();
        assert!((resolved + total.wall_impulse[a]).abs() < 1e-12);
    }
}

#[test]
fn g4_source_budgets_and_cancel_refuse_before_a_run_is_published() {
    let cfg = config();
    let layout = SlotJet3dPatchLayout {
        cells_x: 1,
        cells_z: 1,
    };
    let mut small = budget();
    small.maximum_record_bytes = 1;
    assert!(
        run_slot_jet_3d_patches(&cfg, layout, small, None, &CancelGate::new_clock_free()).is_err()
    );
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        run_slot_jet_3d_patches(&cfg, layout, budget(), None, &gate),
        Err(SlotJet3dSgsError::Source(AeroacError::Cancelled { .. }))
    ));
}

#[test]
fn g1_complete_spatial_run_keeps_the_legacy_aggregate_force_sign() {
    let cfg = config();
    let spatial = run_slot_jet_3d_patches(
        &cfg,
        SlotJet3dPatchLayout {
            cells_x: 2,
            cells_z: 3,
        },
        budget(),
        None,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let total = run_slot_jet_3d(&cfg).unwrap();
    assert_eq!(spatial.force_on_fluid_lu.len(), 6 * cfg.steps_record);
    for (sample, aggregate) in spatial
        .force_on_fluid_lu
        .chunks_exact(6)
        .zip(&total.force_series)
    {
        for a in 0..2 {
            assert!((sample.iter().map(|p| p[a]).sum::<f64>() + aggregate[a]).abs() < 2e-11);
        }
    }
    assert_eq!(
        spatial.diagnostics.mach_max_lattice,
        total.diagnostics.mach_max_lattice
    );
}
