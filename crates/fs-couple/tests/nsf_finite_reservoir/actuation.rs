//! Heat and moving volume compose with the original finite flow exchange.
use super::*;
use fs_lbm::nsf3::NsfGasMedium;

fn drive(volume_rate_m3_s: f64, heat_power_w: f64) -> NsfDgDrive {
    NsfDgDrive::FiniteReservoirActuation {
        aperture: 0,
        volume_rate_m3_s,
        heat_power_w,
        maximum_relative_volume_change: 0.05,
    }
}

#[test]
fn actuation_reuses_original_midpoint_gas_after_identical_flow_and_feeds_next_interval() {
    let gate = CancelGate::new_clock_free();
    for (pressure, upstream, substeps) in [
        (100_002.0, false, 1), (100_002.0, true, 3), (100_000.0, true, 3),
    ] {
        for (rate, power) in [(-2e-7, 0.0), (2e-7, 0.0), (0.0, 0.02), (-2e-7, 0.02)] {
            let attach = || bare(pressure, upstream, substeps)
                .with_finite_reservoir(0, VOLUME, source_limits(), &gate).unwrap();
            let mut passive = attach();
            let mut actuated = attach();
            let dt = actuated.load().propagator().dt();
            let NsfGasMedium::Thermochemical(medium) = actuated.zone().medium() else {
                panic!("original thermochemical source required");
            };
            let solve = medium.energy_solve();
            let flow = passive.step_with_drive(
                NsfDgDrive::FiniteReservoir { aperture: 0 }, &gate).unwrap();
            let flowed = passive.finite_reservoir().unwrap();
            let expected = flowed.midpoint_heat_volume(
                flowed.volume_m3() + rate * dt, power * dt, solve, || false).unwrap();
            let frame = actuated.step_with_drive(drive(rate, power), &gate).unwrap();

            // The first interval holds the same starting gas: motion and heat
            // act on the original post-flow parcel, without changing its flux.
            assert_eq!(actuated.zone().state_hash(), passive.zone().state_hash());
            assert_eq!(actuated.load().state(), passive.load().state());
            assert_eq!(actuated.upstream_load().map(|load| load.state()),
                passive.upstream_load().map(|load| load.state()));
            assert_eq!(frame.reservoir_exchanges, flow.reservoir_exchanges);
            assert_eq!(actuated.finite_reservoir().unwrap(), expected.endpoint());
            let account = frame.finite_reservoir.as_ref().unwrap();
            let motion = account.actuation.as_ref().unwrap();
            assert_eq!(motion.step, *expected.step());
            assert_eq!(motion.midpoint_pressure_pa, expected.midpoint().state().pressure);
            assert_eq!(account.volume_m3, VOLUME);
            assert_eq!(account.final_volume_m3, expected.endpoint().volume_m3());
            assert_eq!(account.final_mass_kg, flowed.mass_kg());
            assert!((account.sensible_energy_change_j
                - flow.finite_reservoir.as_ref().unwrap().sensible_energy_change_j
                - motion.step.stored_energy_change_j).abs() <= 1e-9);
            assert!((account.sensible_energy_change_j - account.source_exchange.total_energy_out_j
                - motion.step.supplied_heat_j + motion.step.work_by_gas_j).abs() <= 1e-9);
            assert!(account.energy_residual_j.abs() <= 1e-9);

            let before_next = *expected.endpoint().state().state();
            let flowed_state = flowed.state().state();
            if rate <= 0.0 {
                assert!(before_next.temperature > flowed_state.temperature);
                assert!(before_next.pressure > flowed_state.pressure);
            } else {
                assert!(before_next.temperature < flowed_state.temperature);
                assert!(before_next.pressure < flowed_state.pressure);
            }
            let next = actuated.step_with_drive(
                NsfDgDrive::FiniteReservoir { aperture: 0 }, &gate).unwrap();
            let next_source = next.finite_reservoir.unwrap();
            assert_eq!(next_source.initial_pressure_pa, before_next.pressure);
            assert_eq!(next_source.initial_temperature_k, before_next.temperature);
            assert_eq!(next_source.volume_m3, expected.endpoint().volume_m3());
            let held = actuated.zone().reservoir_aperture_at(0).unwrap().state;
            assert_eq!(held.temperature_k, before_next.temperature);
            assert!((held.density_kg_m3 - before_next.density).abs()
                <= 32.0 * f64::EPSILON * before_next.density);
        }
    }
}

#[test]
fn explicit_zero_actuation_preserves_the_passive_physical_state_exactly() {
    let gate = CancelGate::new_clock_free();
    let attach = || bare(100_002.0, true, 3)
        .with_finite_reservoir(0, VOLUME, source_limits(), &gate).unwrap();
    let mut passive = attach();
    let mut actuated = attach();
    passive.step_with_drive(NsfDgDrive::FiniteReservoir { aperture: 0 }, &gate).unwrap();
    let frame = actuated.step_with_drive(drive(0.0, 0.0), &gate).unwrap();
    assert_eq!(snapshot!(actuated), snapshot!(passive));
    let account = frame.finite_reservoir.unwrap();
    assert_eq!(account.volume_m3, account.final_volume_m3);
    let motion = account.actuation.unwrap();
    assert_eq!(motion.step.supplied_heat_j, 0.0);
    assert_eq!(motion.step.work_by_gas_j, 0.0);
    assert_eq!(motion.step.stored_energy_change_j, 0.0);
    assert_eq!(motion.step.iterations, 0);
}

#[test]
fn cancelled_or_invalid_actuation_preserves_the_accepted_source_fluid_and_both_ducts() {
    let gate = CancelGate::new_clock_free();
    let mut owner = bare(100_002.0, true, 3)
        .with_finite_reservoir(0, VOLUME, source_limits(), &gate).unwrap();
    owner.step_with_drive(NsfDgDrive::FiniteReservoir { aperture: 0 }, &gate).unwrap();
    let before = snapshot!(owner);
    let dt = owner.load().propagator().dt();
    for (aperture, rate, power, bound) in [
        (0, -0.2 * VOLUME / dt, 0.0, 0.01),
        (0, -2.0 * VOLUME / dt, 0.0, 0.5),
        (0, 0.0, 1e9, 0.05),
        (0, 0.0, 0.0, 0.0),
        (0, 0.0, 0.0, 1.0),
        (0, 0.0, 0.0, f64::NAN),
        (0, f64::INFINITY, 0.0, 0.05),
        (0, 0.0, f64::NAN, 0.05),
        (1, 0.0, 0.0, 0.05),
    ] {
        assert!(owner.step_with_drive(NsfDgDrive::FiniteReservoirActuation {
            aperture, volume_rate_m3_s: rate, heat_power_w: power,
            maximum_relative_volume_change: bound,
        }, &gate).is_err());
        assert_eq!(snapshot!(owner), before);
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(owner.step_with_drive(drive(-2e-7, 0.02), &cancelled).is_err());
    assert_eq!(snapshot!(owner), before);
    owner.step_with_drive(drive(-2e-7, 0.02), &gate).unwrap();
    assert_eq!(owner.accepted_steps(), 2);
}

#[test]
fn opposing_heat_and_work_cannot_hide_excessive_source_energy_movement() {
    let gate = CancelGate::new_clock_free();
    let limits = NsfFiniteReservoirLimits {
        maximum_relative_energy_exchange: 1e-5, ..source_limits()
    };
    let relative_budget = limits.maximum_relative_energy_exchange;
    let mut owner = bare(100_000.0, true, 3)
        .with_finite_reservoir(0, VOLUME, limits, &gate).unwrap();
    owner.step_with_drive(NsfDgDrive::FiniteReservoir { aperture: 0 }, &gate).unwrap();
    let before = snapshot!(owner);
    let inventory = owner.finite_reservoir_inventory().unwrap();
    let budget = inventory.sensible_energy_j.abs() * relative_budget;
    let volume_change = 4.0 * budget / inventory.gas.pressure;
    let heat = inventory.gas.pressure * volume_change;
    let NsfGasMedium::Thermochemical(medium) = owner.zone().medium() else {
        panic!("original thermochemical source required");
    };
    let expected = owner.finite_reservoir().unwrap().midpoint_heat_volume(
        inventory.volume_m3 + volume_change, heat, medium.energy_solve(), || false).unwrap();
    assert!(expected.step().supplied_heat_j > budget);
    assert!(expected.step().work_by_gas_j > budget);
    assert!(expected.step().stored_energy_change_j.abs() < budget);
    let dt = owner.load().propagator().dt();
    assert!(owner.step_with_drive(drive(volume_change / dt, heat / dt), &gate).is_err());
    assert_eq!(snapshot!(owner), before);
}
