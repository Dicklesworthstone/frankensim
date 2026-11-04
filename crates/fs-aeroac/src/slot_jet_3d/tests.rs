use super::*;
use fs_lbm::d3q19::{OPP3, collide_cell3};
use fs_lbm::sgs3::{SgsError, SgsModel, SgsUse};

fn collision() -> CollisionModel3 {
    CollisionModel3::CentralMoment {
        second_order_rate: 1.3,
        higher_order_rate: 1.6,
    }
}

fn closure(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn fixture() -> Rig3 {
    let mut grid = BoundaryGrid3::with_collision_model(
        8,
        8,
        4,
        collision(),
        [0.0; 3],
        BoundarySpec3::periodic(),
    );
    grid.voxelize_sdf(|[x, y, _]| {
        if (3.0..4.0).contains(&x) && (3.0..5.0).contains(&y) {
            -1.0
        } else {
            1.0
        }
    });
    for z in 0..4 {
        for y in 0..8 {
            for x in 0..8 {
                if grid.is_solid(x, y, z) {
                    continue;
                }
                let velocity = [
                    0.018 + 0.002 * (0.7 * y as f64).sin(),
                    0.007 * (0.5 * x as f64).cos(),
                    0.003 * (0.9 * z as f64).sin(),
                ];
                let mut f = equilibrium3(1.0 + 0.001 * ((x + y + z) as f64).cos(), velocity);
                // A mass- and momentum-free non-equilibrium normal stress.
                let stress = 0.0007 * (x + 1) as f64 / 8.0;
                f[0] -= 2.0 * stress;
                f[1] += stress;
                f[OPP3[1]] += stress;
                grid.set_populations(x, y, z, &f);
            }
        }
    }
    let plate_links = grid.boundary_links();
    Rig3 {
        grid,
        plate_links,
        // Disable only the external fringe source in this closed-box momentum
        // oracle. Collision, streaming, voxel walls and loads are production.
        fringe: Fringe3 {
            width: 1,
            sigma_max: 0.0,
            targets: vec![equilibrium3(1.0, [0.0; 3]); 8],
        },
        x_plate_plane: 1,
        x_fringe_plane: 6,
    }
}

fn populations(grid: &BoundaryGrid3) -> Vec<[u64; Q3]> {
    let [nx, ny, nz] = grid.dimensions();
    let mut values = Vec::new();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if !grid.is_solid(x, y, z) {
                    values.push(grid.populations(x, y, z).map(f64::to_bits));
                }
            }
        }
    }
    values
}

fn momentum(grid: &BoundaryGrid3) -> [f64; 3] {
    let [nx, ny, nz] = grid.dimensions();
    let mut result = [0.0; 3];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if grid.is_solid(x, y, z) {
                    continue;
                }
                for (f, e) in grid.populations(x, y, z).into_iter().zip(E3) {
                    result[0] += f * f64::from(e.0);
                    result[1] += f * f64::from(e.1);
                    result[2] += f * f64::from(e.2);
                }
            }
        }
    }
    result
}

#[test]
fn g1_plate_record_uses_the_collided_impulse_and_closes_fluid_momentum() {
    let mut rig = fixture();
    let before = momentum(&rig.grid);
    let mut expected = [0.0; 2];
    for link in &rig.plate_links {
        let [x, y, z] = link.cell;
        let post = collide_cell3(rig.grid.populations(x, y, z), collision(), [0.0; 3]).unwrap();
        let incoming = E3[link.direction];
        expected[0] -= 2.0 * post[OPP3[link.direction]] * f64::from(incoming.0);
        expected[1] -= 2.0 * post[OPP3[link.direction]] * f64::from(incoming.1);
    }
    let actual = rig.step_plate_impulse();
    let after = momentum(&rig.grid);
    for axis in 0..2 {
        assert!((actual[axis] - expected[axis]).abs() < 2e-14);
        assert!((after[axis] - before[axis] + actual[axis]).abs() < 3e-13);
    }
    let mut old_proxy = [0.0; 2];
    for link in &rig.plate_links {
        let [x, y, z] = link.cell;
        let f = rig.grid.populations(x, y, z)[OPP3[link.direction]];
        let incoming = E3[link.direction];
        old_proxy[0] -= 2.0 * f * f64::from(incoming.0);
        old_proxy[1] -= 2.0 * f * f64::from(incoming.1);
    }
    assert!((old_proxy[0] - actual[0]).abs() + (old_proxy[1] - actual[1]).abs() > 1e-8);
}

#[test]
fn g0_g3_native_sgs_preserves_zero_bits_and_uses_the_same_collided_loads() {
    let gate = CancelGate::new_clock_free();
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let mut ordinary = fixture();
        let mut sgs = fixture();
        let mut workspace = sgs.grid.sgs_workspace().unwrap();
        for _ in 0..3 {
            let expected = ordinary.step_plate_impulse();
            let report = sgs
                .step_sgs(&gate, closure(model), &mut workspace, true)
                .unwrap();
            let force = report.wall_momentum.unwrap().wall_impulse;
            assert_eq!(
                expected.map(f64::to_bits),
                [force[0], force[1]].map(f64::to_bits)
            );
            assert_eq!(populations(&ordinary.grid), populations(&sgs.grid));
            assert_eq!(report.sgs.sgs_dissipation_lu, 0.0);
        }
    }
    for model in [
        SgsModel::Smagorinsky { cs: 0.12 },
        SgsModel::Wale { cw: 0.35 },
        SgsModel::Vreman { cv: 0.07 },
    ] {
        let mut rig = fixture();
        let before = momentum(&rig.grid);
        let mut workspace = rig.grid.sgs_workspace().unwrap();
        let report = rig
            .step_sgs(&gate, closure(model), &mut workspace, true)
            .unwrap();
        assert!(report.sgs.maximum_eddy_viscosity_lu > 0.0);
        assert!(report.sgs.sgs_dissipation_lu > 0.0);
        let after = momentum(&rig.grid);
        let exchange = report.wall_momentum.unwrap();
        for axis in 0..3 {
            assert!((after[axis] - before[axis] + exchange.wall_impulse[axis]).abs() < 3e-13);
            assert!(
                (exchange.wall_impulse[axis] + exchange.fluid_population_impulse[axis]).abs()
                    < 1e-14
            );
        }
    }
}

fn config() -> SlotJet3dConfig {
    SlotJet3dConfig {
        nx: 64,
        ny: 32,
        nz: 8,
        slot_half: 2.5,
        u_jet: 0.04,
        collision: collision(),
        nozzle_thickness: 1,
        edge_distance: 8,
        plate_length: 2,
        fringe_width: 8,
        fringe_sigma: 0.4,
        seed_amplitude: 0.02,
        steps_settle: 3,
        steps_record: 64,
    }
}

#[test]
fn g3_chunk_partitions_record_the_same_plate_links_and_step_forces() {
    let cfg = config();
    let mut ordinary = Rig3::build(&cfg);
    assert!(!ordinary.plate_links.is_empty());
    for link in &ordinary.plate_links {
        let e = E3[link.direction];
        let wall_x = (link.cell[0] as isize - e.0 as isize).rem_euclid(cfg.nx as isize) as usize;
        let wall_y = (link.cell[1] as isize - e.1 as isize).rem_euclid(cfg.ny as isize) as usize;
        assert!((9..11).contains(&wall_x));
        assert!(wall_y == 15 || wall_y == 16);
    }
    ordinary.settle(3).unwrap();
    let expected = ordinary.record(5).unwrap();
    let mut chunked = Rig3::build(&cfg);
    let mut state = ChunkState {
        settle_done: 0,
        record_done: 0,
        force_series: Vec::new(),
    };
    for budget in [2, 3, 3] {
        step_chunk_loop(&cfg, &mut chunked, &mut state, budget).unwrap();
    }
    assert_eq!((state.settle_done, state.record_done), (3, 5));
    assert_eq!(expected, state.force_series);
    assert_eq!(populations(&ordinary.grid), populations(&chunked.grid));
}

#[test]
fn g0_closure_refusals_and_cancellation_do_not_publish_a_source_step() {
    let gate = CancelGate::new_clock_free();
    let mut invalid = closure(SgsModel::Smagorinsky { cs: 0.1 });
    invalid.usage = SgsUse::AcceptanceRender;
    assert!(matches!(
        run_slot_jet_3d_with_sgs(&config(), invalid, &gate),
        Err(SlotJet3dSgsError::Flow(BoundarySgsError3::Closure(
            SgsError::AcceptanceRenderForbidden
        )))
    ));
    gate.request();
    assert!(matches!(
        run_slot_jet_3d_with_sgs(&config(), closure(SgsModel::Wale { cw: 0.2 }), &gate),
        Err(SlotJet3dSgsError::Source(AeroacError::Cancelled { .. }))
    ));
    let mut rig = fixture();
    let before = populations(&rig.grid);
    let mut workspace = rig.grid.sgs_workspace().unwrap();
    let mut bounded = closure(SgsModel::Smagorinsky { cs: 0.1 });
    bounded.maximum_eddy_to_molecular = 0.0;
    assert!(matches!(
        rig.step_sgs(&CancelGate::new_clock_free(), bounded, &mut workspace, true),
        Err(BoundarySgsError3::Closure(
            SgsError::ViscosityRatioExceeded { .. }
        ))
    ));
    assert_eq!(before, populations(&rig.grid));
}

#[test]
fn g1_material_sgs_source_preserves_physical_force_and_molecular_scaling() {
    use crate::material::{JetPhysicalState, MaterialSlotJet3d, SlotJet3dSetup};
    use fs_material::gas::GasState;

    let cfg = config();
    let setup = SlotJet3dSetup {
        nx: cfg.nx,
        ny: cfg.ny,
        nz: cfg.nz,
        slot_half: cfg.slot_half,
        higher_order_rate: 1.6,
        nozzle_thickness: cfg.nozzle_thickness,
        edge_distance: cfg.edge_distance,
        plate_length: cfg.plate_length,
        fringe_width: cfg.fringe_width,
        fringe_sigma: cfg.fringe_sigma,
        seed_amplitude: cfg.seed_amplitude,
        steps_settle: 1,
        steps_record: 64,
    };
    let physical = JetPhysicalState {
        cell_size_m: 1e-4,
        timestep_s: 1e-5,
        jet_speed_m_s: 0.4,
        max_speed_m_s: 1.2,
        max_cells: 64 * 32 * 8,
    };
    let gas = GasState::try_new_moist_air(300.0, 100_000.0, 0.0).unwrap();
    let source = MaterialSlotJet3d::try_new(&gas, &setup, physical).unwrap();
    let molecular = source.lattice_config().collision;
    let model = closure(SgsModel::Smagorinsky { cs: 0.12 });
    let gate = CancelGate::new_clock_free();
    let (lattice, diagnostics) =
        run_slot_jet_3d_with_sgs(source.lattice_config(), model, &gate).unwrap();
    let (run, material_diagnostics) = source.run_with_sgs(model, &gate).unwrap();
    assert_eq!(diagnostics, material_diagnostics);
    assert!(diagnostics.recorded.sgs_dissipation_lu > 0.0);
    assert_eq!(molecular, source.lattice_config().collision);
    let force_unit = gas.density * physical.cell_size_m.powi(4) / physical.timestep_s.powi(2);
    for (actual, native) in run.force_series_n().iter().zip(&lattice.force_series) {
        for axis in 0..2 {
            let expected = force_unit * native[axis];
            assert!((actual[axis] - expected).abs() < 2e-14 * expected.abs().max(1e-20));
        }
    }
    assert_eq!(run.timestep_s(), physical.timestep_s);
}
