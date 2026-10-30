//! Independent physical quadrature and original cached-work interval checks.
use super::*;
use crate::nsf3::NsfBoundaryDrive;

#[test]
fn g1_fv_reservoir_pressure_work_retains_spatial_product_covariance() {
    for axis in 1..3 {
        for high in [false, true] {
            let mut grid = aperture_grid(0, 1)
                .with_reservoir_aperture(aperture(axis, high, 293.15)).unwrap();
            let reference = grid.characteristic_reference_gas(false).unwrap();
            let mut state = grid.reservoir_aperture_state_at(reference.pressure + 5.0, 293.15).unwrap();
            state.velocity_m_s[axis] = -0.3;
            grid.set_reservoir_aperture_state(state).unwrap();
            let (mut rho, momentum) = grid.macro_fields(&grid.f);
            let mut velocity = grid.velocity_field(&rho, &momentum, &grid.force);
            let theta = grid.theta.clone();
            let mut energy = grid.e.clone();
            let speed_unit = grid.law.dx_m / grid.law.dt_s;
            let area = grid.law.dx_m.powi(2);
            let sign = if high { 1.0 } else { -1.0 };
            let mut pressure_sum = 0.0;
            let mut acoustic_sum = 0.0;
            let mut pressure_power = 0.0;
            let mut acoustic_power = 0.0;
            let mut midpoint_acoustic_power = 0.0;
            let bound = grid.reservoir_aperture.as_ref().unwrap();
            for (n, &index) in bound.indices.iter().enumerate() {
                let perturbation = 2.0 * (n + 1) as f64;
                let pressure = reference.pressure + perturbation;
                let speed = 1.0 + 0.1 * n as f64;
                rho[index] = pressure / (reference.specific_gas_constant * 293.15 * grid.law.rho_ref_kg_m3);
                velocity[index][axis] = speed / speed_unit;
                energy[index] = rho[index] * (theta[index] / (grid.law.gamma - 1.0)
                    + 0.5 * velocity[index].iter().map(|v| v * v).sum::<f64>());
                pressure_sum += 0.5 * (pressure + reference.pressure + 5.0);
                acoustic_sum += 0.5 * (perturbation + 5.0);
                pressure_power += sign * area * 0.5
                    * (pressure * speed + (reference.pressure + 5.0) * -0.3);
                acoustic_power += sign * area * 0.5 * (perturbation * speed + 5.0 * -0.3);
                midpoint_acoustic_power += sign * area * 0.25
                    * (perturbation + 5.0) * (speed - 0.3);
            }
            let before = (grid.state_hash(), grid.total_energy_j(), grid.last_reservoir_exchange);
            let (_, _, flux) = grid.flux_divergence(&energy, &rho, &velocity, &theta,
                &vec![[0.0; 6]; grid.n], &[], grid.open_faces.as_ref()).unwrap();
            assert!((flux.reservoir_pressure_sum_pa - pressure_sum).abs() < 1e-9);
            assert!((flux.reservoir_acoustic_pressure_sum_pa - acoustic_sum).abs() < 1e-9);
            assert!((flux.reservoir_pressure_power_w - pressure_power).abs() < 2e-14 * pressure_power.abs());
            assert!((flux.reservoir_acoustic_pressure_power_w - acoustic_power).abs() < 2e-11 * acoustic_power.abs());
            assert!((acoustic_power - midpoint_acoustic_power).abs() > 0.1 * acoustic_power.abs(),
                "the manufactured spatial covariance must distinguish the two quadratures");
            assert_eq!((grid.state_hash(), grid.total_energy_j(), grid.last_reservoir_exchange), before);
        }
    }
}

#[test]
fn g3_reservoir_work_uses_the_same_cached_predictors_as_total_energy() {
    let mut grid = aperture_grid(0, 1)
        .with_reservoir_aperture(aperture(1, true, 293.15)).unwrap();
    assert!(grid.last_reservoir_exchange().is_none());
    let mut source = grid.reservoir_aperture_state_at(101_345.0, 295.0).unwrap();
    source.velocity_m_s[1] = 0.3;
    grid.set_reservoir_aperture_state(source).unwrap();
    for _ in 0..4 {
        let old_flux = grid.open_out;
        let before = (grid.total_mass_kg(), grid.total_energy_j());
        let ledger = grid.step().unwrap();
        let after_flux = grid.open_out;
        let exchange = *grid.last_reservoir_exchange().unwrap();
        let dt = grid.law.dt_s;
        assert_eq!(exchange.reference_pressure_pa, 101_325.0);
        assert_eq!(exchange.area_m2, 9.0 * grid.law.dx_m.powi(2));
        assert_eq!(exchange.pressure_work_out_j, 0.5 * dt
            * (old_flux.reservoir_pressure_power_w + after_flux.reservoir_pressure_power_w));
        assert_eq!(exchange.acoustic_pressure_work_out_j, 0.5 * dt
            * (old_flux.reservoir_acoustic_pressure_power_w + after_flux.reservoir_acoustic_pressure_power_w));
        assert_eq!(exchange.mean_pressure_pa,
            0.5 * (old_flux.reservoir_pressure_sum_pa + after_flux.reservoir_pressure_sum_pa) / 9.0);
        assert_eq!(exchange.normal_volume_flow_m3_s,
            grid.last_reservoir_flows().iter().map(|flow| flow.outward_volume_flow_m3_s).sum::<f64>());
        assert_eq!(exchange.projection_covariance_j, exchange.acoustic_pressure_work_out_j
            - exchange.mean_acoustic_pressure_pa * exchange.normal_volume_flow_m3_s * dt);
        let mass_out = grid.last_reservoir_mass_out_kg().unwrap()
            + grid.last_boundary_exchange().unwrap().iter().map(|face| face.mass_out_kg).sum::<f64>();
        assert!((grid.total_mass_kg() - before.0 + mass_out).abs() < 5e-12 * before.0);
        assert!((grid.total_energy_j() - before.1 - ledger.wall_heat_j
            - ledger.open_boundary_energy_j - ledger.open_energy_j).abs() < 5e-12 * before.1);
    }
    let previous = *grid.last_reservoir_exchange().unwrap();
    let physical = (grid.total_mass_kg(), grid.total_energy_j());
    let source = grid.reservoir_aperture_state_at(101_335.0, 294.0).unwrap();
    grid.set_reservoir_aperture_state(source).unwrap();
    assert_eq!(grid.last_reservoir_exchange(), Some(&previous));
    assert!((grid.total_mass_kg() - physical.0).abs() < 2e-14 * physical.0);
    assert_eq!(grid.total_energy_j(), physical.1);
}

#[test]
fn g4_g5_reservoir_work_intervals_rollback_and_replay_across_workers() {
    let build = |workers| aperture_grid(0, workers)
        .with_reservoir_aperture(aperture(1, true, 293.15)).unwrap();
    let mut grid = build(1);
    let mut replay = build(3);
    let drives = [NsfBoundaryDrive::PressureOffset(2.0), NsfBoundaryDrive::PressureOffset(0.0)];
    let source = grid.reservoir_aperture_state_at(101_345.0, 293.35).unwrap();
    let before = (grid.state_hash(), grid.last_reservoir_exchange, grid.last_reservoir_flows.clone(),
        grid.cumulative_ledger());
    let work = 3 * grid.cell_count();
    let trial = grid.preview_boundary_interval_with_reservoir(drives, None, source, 3, work, || false).unwrap();
    let expected = (*trial.reservoir_exchange().unwrap(), trial.grid().state_hash(), *trial.energy_ledger());
    let frame = expected.0;
    assert_eq!(frame.projection_covariance_j, frame.acoustic_pressure_work_out_j
        - frame.mean_acoustic_pressure_pa * frame.normal_volume_flow_m3_s * (3.0 * trial.grid().law.dt_s));
    let original_flows: f64 = trial.reservoir_flows().iter().map(|flow| flow.outward_volume_flow_m3_s).sum();
    assert!((frame.normal_volume_flow_m3_s - original_flows).abs() < 2e-13 * original_flows.abs().max(1e-30));
    let mass_out = trial.reservoir_mass_out_kg()
        + trial.boundary_exchange().iter().map(|face| face.mass_out_kg).sum::<f64>();
    assert!((trial.grid().total_mass_kg() - trial.previous_mass_kg() + mass_out).abs()
        < 5e-12 * trial.previous_mass_kg());
    assert!(trial.energy_ledger().residual_rel < 5e-12);
    drop(trial);
    assert_eq!((grid.state_hash(), grid.last_reservoir_exchange, grid.last_reservoir_flows.clone(),
        grid.cumulative_ledger()), before);
    let calls = std::cell::Cell::new(0);
    assert!(grid.preview_boundary_interval_with_reservoir(drives, None, source, 3, work, || {
        calls.set(calls.get() + 1); calls.get() == 3
    }).is_err());
    assert_eq!((grid.state_hash(), grid.last_reservoir_exchange, grid.last_reservoir_flows.clone(),
        grid.cumulative_ledger()), before);
    let retry = grid.preview_boundary_interval_with_reservoir(drives, None, source, 3, work, || false).unwrap();
    assert_eq!((*retry.reservoir_exchange().unwrap(), retry.grid().state_hash(), *retry.energy_ledger()), expected);
    retry.accept();
    let other = replay.preview_boundary_interval_with_reservoir(drives, None, source, 3, work, || false).unwrap();
    assert_eq!((*other.reservoir_exchange().unwrap(), other.grid().state_hash(), *other.energy_ledger()), expected);
    other.accept();
    assert_eq!(grid.last_reservoir_exchange(), replay.last_reservoir_exchange());
}
