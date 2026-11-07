use super::*;

#[test]
fn interval_keeps_each_native_flux_history_and_all_step_energy_accounts() {
    let steps = 4;
    let walls = vec![
        face(1, false, NsfWallThermal::Isothermal(310.0)),
        face(
            1,
            true,
            NsfWallThermal::Film {
                coefficient_w_m2_k: 1000.0,
                ambient_temperature_k: 290.0,
            },
        ),
        face(2, true, NsfWallThermal::HeatFlux(125.0)),
    ];
    let mut direct = grid();
    direct
        .set_thermal_walls_recorded(walls.clone(), || false)
        .unwrap();
    let native: Vec<_> = (0..steps).map(|_| direct.step().unwrap()).collect();
    let mut fluid = grid();
    assert_eq!(fluid.cell_count(), 64);
    let trial = fluid
        .preview_thermal_wall_interval(walls.clone(), steps, 64, 64 * steps, || false)
        .unwrap();
    assert_eq!(trial.grid().state_hash(), direct.state_hash());
    assert_eq!(trial.grid().div, direct.div);
    assert_eq!(trial.grid().cumulative_ledger(), direct.cumulative_ledger());
    assert_eq!(trial.grid().steps(), steps);
    let ledger = *trial.energy_ledger();
    let expected = EnergyLedger {
        internal_j: native.last().unwrap().internal_j,
        kinetic_j: native.last().unwrap().kinetic_j,
        viscous_heating_j: native.iter().map(|x| x.viscous_heating_j).sum(),
        unresolved_internal_exchange_j: native
            .iter()
            .map(|x| x.unresolved_internal_exchange_j)
            .sum(),
        pressure_work_j: native.iter().map(|x| x.pressure_work_j).sum(),
        wall_heat_j: native.iter().map(|x| x.wall_heat_j).sum(),
        wall_work_j: native.iter().map(|x| x.wall_work_j).sum(),
        open_boundary_energy_j: native.iter().map(|x| x.open_boundary_energy_j).sum(),
        open_energy_j: native.iter().map(|x| x.open_energy_j).sum(),
        source_energy_j: native.iter().map(|x| x.source_energy_j).sum(),
        gravity_work_j: native.iter().map(|x| x.gravity_work_j).sum(),
        body_force_work_j: native.iter().map(|x| x.body_force_work_j).sum(),
        residual_rel: ledger.residual_rel,
    };
    assert_eq!(ledger, expected);
    near(
        trial
            .wall_exchange()
            .iter()
            .map(|x| x.heat_into_fluid_j)
            .sum(),
        ledger.wall_heat_j,
        1e-21,
    );
    let prescribed = trial.wall_exchange().iter().find(|x| x.axis == 2).unwrap();
    let expected_heat = 125.0 * direct.law().dx_m.powi(2) * direct.law().dt_s * steps as f64;
    near(
        prescribed.heat_into_fluid_j,
        expected_heat,
        5e-14 * expected_heat.abs(),
    );
    assert!(ledger.residual_rel < 1e-11);
    trial.accept();
    assert_eq!(fluid.state_hash(), direct.state_hash());
    // Resetting the boundary history at EVERY substep would silently replace
    // the corrector's retained flux with a different accepted-temperature query.
    let mut restarted = grid();
    for _ in 0..steps {
        restarted
            .set_thermal_walls_recorded(walls.clone(), || false)
            .unwrap();
        restarted.step().unwrap();
    }
    assert_ne!(restarted.state_hash(), direct.state_hash());
}

#[test]
fn interval_caps_count_solid_cells_and_refuse_zero_or_overflow_before_copying() {
    let mut fluid =
        NsfGrid3::<D3q27>::new(spec(), |cell| cell == [1, 1, 1], |_| (1.0, [0.0; 3], 300.0))
            .unwrap();
    assert_eq!(fluid.cell_count(), 64);
    let before = (
        fluid.state_hash(),
        fluid.cumulative_ledger(),
        fluid.steps(),
        fluid.div.clone(),
    );
    let walls = vec![face(1, false, NsfWallThermal::HeatFlux(250.0))];
    for (steps, cells, work) in [
        (0, 64, 64),
        (2, 63, 128),
        (2, 64, 127),
        (usize::MAX, 64, usize::MAX),
    ] {
        assert!(
            fluid
                .preview_thermal_wall_interval(walls.clone(), steps, cells, work, || false)
                .is_err()
        );
        assert_eq!(
            (
                fluid.state_hash(),
                fluid.cumulative_ledger(),
                fluid.steps(),
                fluid.div.clone()
            ),
            before
        );
        assert!(fluid.thermal_walls().is_empty());
    }
    let trial = fluid
        .preview_thermal_wall_interval(walls.clone(), 2, 64, 128, || false)
        .unwrap();
    assert_eq!(trial.grid().steps(), 2);
    drop(trial);
    assert_eq!(
        (
            fluid.state_hash(),
            fluid.cumulative_ledger(),
            fluid.steps(),
            fluid.div.clone()
        ),
        before
    );
    fluid.steps = usize::MAX;
    assert!(matches!(
        fluid.preview_thermal_wall_interval(walls, 1, 64, 64, || false),
        Err(NsfError::Input("thermal wall trial clock overflow"))
    ));
    assert_eq!(fluid.steps(), usize::MAX);
}
