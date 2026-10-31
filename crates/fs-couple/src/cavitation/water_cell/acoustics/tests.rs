use super::super::{EquilibriumWaterCell, WaterCellControls};
use super::*;
use fs_duct::fluid::segment_wave;
use fs_material::phase::water::isochoric::{IsochoricWater, IsochoricWaterControls};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::MassSpecificVolume;
use fs_thermochem::water::if97::{If97Region, If97SaturationPair, If97State};

fn controls(minimum_temperature_k: f64, maximum_temperature_k: f64) -> WaterCellControls {
    WaterCellControls {
        equilibrium: IsochoricWaterControls {
            minimum_temperature_k,
            maximum_temperature_k,
            minimum_pressure_pa: 100.0,
            maximum_pressure_pa: 100e6,
            volume_tolerance_m3_kg: 1e-12,
            internal_energy_tolerance_j_kg: 1e-6,
            maximum_pressure_iterations: 96,
            maximum_temperature_iterations: 96,
        },
        maximum_relative_mass_exchange: 0.25,
        maximum_inflow_parcels: 16,
        mass_tolerance_kg: 1e-12,
        energy_tolerance_j: 1e-6,
        volume_tolerance_m3: 1e-10,
    }
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-20),
        "{actual} != {expected}"
    );
}

#[test]
fn finite_phase_inventory_drives_independent_wood_compliance_and_native_duct() {
    let cell = EquilibriumWaterCell::at_temperature(
        0.1,
        1.0,
        Temperature::new(400.0),
        0.0,
        controls(380.0, 430.0),
        &mut || false,
    )
    .unwrap();
    let accepted = cell.state();
    let acoustic = FrozenWaterCellAcoustics::new(&accepted).unwrap();
    assert!(core::ptr::eq(acoustic.source(), &accepted));
    let phase = accepted.source_state();
    let pair = phase.saturation_pair().unwrap();
    let x = phase.vapor_mass_fraction();
    assert!(x > 0.0 && x < 1.0);
    let liquid_volume = (1.0 - x) * pair.liquid().specific_volume().value();
    let vapor_volume = x * pair.vapor().specific_volume().value();
    let alpha = vapor_volume / (liquid_volume + vapor_volume);
    let rho = 1.0 / (liquid_volume + vapor_volume);
    let rho_l = pair.liquid().density().value();
    let rho_v = pair.vapor().density().value();
    let c_l = pair.liquid().sound_speed().value();
    let c_v = pair.vapor().sound_speed().value();
    let compliance = (1.0 - alpha) / (rho_l * c_l * c_l) + alpha / (rho_v * c_v * c_v);
    let sound_speed = (1.0 / (rho * compliance)).sqrt();
    near(acoustic.state().density_kg_m3(), rho, 5e-13);
    near(acoustic.state().sound_speed_m_s(), sound_speed, 5e-13);
    near(acoustic.mixture().unwrap().void_fraction(), alpha, 5e-13);
    near(accepted.vapor_volume_fraction(), alpha, 5e-13);
    assert!((alpha - x).abs() > 0.1); // Volume fraction cannot be mass quality.

    let omega = 2.0 * core::f64::consts::PI * 100.0;
    let radius: f64 = 0.005;
    let wave = segment_wave(&acoustic, radius, omega).unwrap();
    near(wave.wavenumber.re, omega / sound_speed, 5e-13);
    near(
        wave.characteristic_impedance.re,
        rho * sound_speed / (core::f64::consts::PI * radius * radius),
        5e-13,
    );
    assert_eq!(wave.wavenumber.im, 0.0);
    assert_eq!(wave.characteristic_impedance.im, 0.0);
    assert_eq!(
        cell.state().internal_energy_j(),
        accepted.internal_energy_j()
    );
    assert_eq!(cell.state().mass_kg(), accepted.mass_kg());
}

#[test]
fn liquid_and_steam_observations_retain_exact_original_single_phase_response() {
    for (region, pressure_pa, temperature_k) in [
        (If97Region::Region1, 10e6, 300.0),
        (If97Region::Region2, 1e5, 500.0),
    ] {
        let original = If97State::at_pressure_temperature(
            Pressure::new(pressure_pa),
            Temperature::new(temperature_k),
            region,
        )
        .unwrap();
        let cell = EquilibriumWaterCell::at_temperature(
            original.specific_volume().value(),
            1.0,
            Temperature::new(temperature_k),
            2.0,
            controls(temperature_k - 1.0, temperature_k + 1.0),
            &mut || false,
        )
        .unwrap();
        let accepted = cell.state();
        let native = accepted.source_state().single_phase_source().unwrap();
        let acoustic = FrozenWaterCellAcoustics::new(&accepted).unwrap();
        assert!(acoustic.mixture().is_none());
        let equilibrium = EquilibriumWaterCellAcoustics::new(&accepted).unwrap();
        assert!(equilibrium.equilibrium().is_none());
        assert_eq!(equilibrium.state(), acoustic.state());
        assert_eq!(
            acoustic.acoustic_state().unwrap(),
            native.acoustic_state().unwrap()
        );
        near(
            acoustic.state().density_kg_m3(),
            original.density().value(),
            1e-8,
        );
        near(
            acoustic.state().sound_speed_m_s(),
            original.sound_speed().value(),
            1e-8,
        );
        assert_eq!(acoustic.source().time_s(), 2.0);
    }
}

#[test]
fn explicit_equilibrium_phase_exchange_changes_native_duct_response_at_same_inventory() {
    let mut cell = EquilibriumWaterCell::at_temperature(
        0.1,
        1.0,
        Temperature::new(400.0),
        0.0,
        controls(380.0, 430.0),
        || false,
    )
    .unwrap();
    let before = cell.state();
    let frozen = FrozenWaterCellAcoustics::new(&before).unwrap();
    let equilibrium = EquilibriumWaterCellAcoustics::new(&before).unwrap();
    assert!(core::ptr::eq(equilibrium.source(), &before));
    assert_eq!(
        equilibrium.equilibrium().unwrap().source(),
        &before.source_state().saturation_pair().unwrap()
    );
    near(
        equilibrium.state().density_kg_m3(),
        frozen.state().density_kg_m3(),
        1e-13,
    );
    assert!(equilibrium.state().sound_speed_m_s() < frozen.state().sound_speed_m_s());
    let omega = 2.0 * core::f64::consts::PI * 10.0;
    let wave = segment_wave(&equilibrium, 0.005, omega).unwrap();
    near(
        wave.wavenumber.re,
        omega / equilibrium.state().sound_speed_m_s(),
        1e-13,
    );
    assert!(wave.wavenumber.re > segment_wave(&frozen, 0.005, omega).unwrap().wavenumber.re);
    assert_eq!(cell.state(), before);

    let step = cell.advance(&[], 0.0, 10_000.0, 0.1, || false).unwrap();
    let after = cell.state();
    let changed = EquilibriumWaterCellAcoustics::new(&after).unwrap();
    assert_eq!(changed.source().time_s(), 0.1);
    assert_ne!(changed.state(), equilibrium.state());
    assert!(after.phase_masses_kg()[1] > before.phase_masses_kg()[1]);
    near(
        after.internal_energy_j() - before.internal_energy_j(),
        10_000.0,
        1e-12,
    );
    assert!(step.energy_residual_j.abs() <= cell.controls().energy_tolerance_j);
    assert_eq!(equilibrium.source().time_s(), 0.0);
}

#[test]
fn accepted_phase_change_updates_propagation_without_mutating_an_earlier_observer() {
    let controls = controls(380.0, 430.0);
    let mut cell = EquilibriumWaterCell::at_temperature(
        0.1,
        1.0,
        Temperature::new(400.0),
        0.0,
        controls,
        &mut || false,
    )
    .unwrap();
    let before = cell.state();
    let old = FrozenWaterCellAcoustics::new(&before).unwrap();
    let original_law =
        IsochoricWater::new(MassSpecificVolume::new(0.1), controls.equilibrium, || false).unwrap();
    let target = original_law
        .state_at_temperature(Temperature::new(410.0), || false)
        .unwrap();
    let heat =
        before.mass_kg() * target.specific_internal_energy_j_kg() - before.internal_energy_j();
    cell.advance(&[], 0.0, heat, 0.25, &mut || false).unwrap();
    let after = cell.state();
    let updated = FrozenWaterCellAcoustics::new(&after).unwrap();
    let pair =
        If97SaturationPair::at_temperature(Temperature::new(after.source_state().temperature_k()))
            .unwrap();
    let expected = FrozenHomogeneousMixture::from_mass_fraction(
        pair.liquid().acoustic_state().unwrap(),
        pair.vapor().acoustic_state().unwrap(),
        after.source_state().vapor_mass_fraction(),
    )
    .unwrap();
    assert_eq!(updated.state(), expected.state());
    assert_eq!(old.source().time_s(), 0.0);
    assert_eq!(updated.source().time_s(), 0.25);
    near(after.mass_kg(), before.mass_kg(), 1e-14);
    near(
        after.internal_energy_j() - before.internal_energy_j(),
        heat,
        1e-12,
    );
    assert!(after.phase_masses_kg()[1] > before.phase_masses_kg()[1]);
    assert!((updated.state().sound_speed_m_s() - old.state().sound_speed_m_s()).abs() > 1.0);
    let omega = 2.0 * core::f64::consts::PI * 100.0;
    let initial_wave = segment_wave(&old, 0.005, omega).unwrap();
    let final_wave = segment_wave(&updated, 0.005, omega).unwrap();
    assert_ne!(initial_wave.wavenumber, final_wave.wavenumber);
    assert_eq!(
        old.state(),
        FrozenWaterCellAcoustics::new(&before).unwrap().state()
    );
}
