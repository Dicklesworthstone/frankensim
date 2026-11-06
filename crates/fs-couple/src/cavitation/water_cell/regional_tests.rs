use super::acoustics::{EquilibriumWaterCellAcoustics, FrozenWaterCellAcoustics};
use super::volume::{WaterVolumeControls, WaterVolumeError};
use super::*;
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;
use fs_thermochem::water::if97::regions::If97FluidSaturationPair;

fn dense_controls() -> DensePressureSolve {
    DensePressureSolve {
        pressure_tolerance_pa: 1e-4,
        density_tolerance_kg_m3: 1e-8,
        maximum_iterations: 128,
    }
}

fn controls() -> WaterCellControls {
    WaterCellControls {
        equilibrium: IsochoricWaterControls {
            minimum_temperature_k: 638.0,
            maximum_temperature_k: 642.0,
            minimum_pressure_pa: 1e5,
            maximum_pressure_pa: 22.064e6,
            volume_tolerance_m3_kg: 1e-14,
            internal_energy_tolerance_j_kg: 1e-4,
            maximum_pressure_iterations: 96,
            maximum_temperature_iterations: 96,
        },
        maximum_relative_mass_exchange: 0.1,
        maximum_inflow_parcels: 8,
        mass_tolerance_kg: 1e-13,
        energy_tolerance_j: 5e-5,
        volume_tolerance_m3: 1e-13,
    }
}

fn pair() -> If97FluidSaturationPair {
    If97FluidSaturationPair::at_temperature(Temperature::new(640.0), dense_controls(), || false)
        .unwrap()
}

fn cell(quality: f64, mass: f64) -> EquilibriumWaterCell {
    let pair = pair();
    let vl = pair.liquid().specific_volume().value();
    let vv = pair.vapor().specific_volume().value();
    let chart = IsochoricWater::new_regional_coexistence(
        MassSpecificVolume::new(vl + quality * (vv - vl)),
        controls().source_controls(mass),
        dense_controls(),
        || false,
    )
    .unwrap();
    let source = chart
        .state_at_temperature(Temperature::new(640.0), || false)
        .unwrap();
    EquilibriumWaterCell::from_isochoric_state(mass, source, 0.0, controls(), || false).unwrap()
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance={tolerance:e}"
    );
}

fn assert_inventory(state: WaterCellState) {
    let source = state.source_state();
    let pair = source.regional_saturation_pair().unwrap();
    assert!(source.law().regional_coexistence_controls().is_some());
    let masses = state.phase_masses_kg();
    let volumes = state.phase_volumes_m3();
    near(masses[0] + masses[1], state.mass_kg(), 1e-13);
    near(volumes[0] + volumes[1], state.volume_m3(), 1e-13);
    near(
        volumes[0],
        masses[0] * pair.liquid().specific_volume().value(),
        1e-16,
    );
    near(
        volumes[1],
        masses[1] * pair.vapor().specific_volume().value(),
        1e-16,
    );
    near(
        source.bulk_density_kg_m3() * source.specific_volume_m3_kg(),
        1.0,
        3e-16,
    );
    assert!(state.source_energy_residual_j().abs() <= controls().energy_tolerance_j);
}

#[test]
fn high_pressure_heat_updates_phase_inventory_and_both_acoustic_closures() {
    let mut water = cell(0.4, 0.2);
    let initial = water.state();
    let original = initial.source_state();
    assert!(original.saturation_pair().is_none());
    assert!(original.temperature_k() > 623.15);
    let target = original
        .law()
        .state_at_temperature(Temperature::new(640.2), || false)
        .unwrap();
    let heat =
        initial.mass_kg() * target.specific_internal_energy_j_kg() - initial.internal_energy_j();
    let step = water.advance(&[], 0.0, heat, 0.1, || false).unwrap();
    let final_state = water.state();
    assert_eq!(final_state.mass_kg(), initial.mass_kg());
    assert_eq!(final_state.volume_m3(), initial.volume_m3());
    near(final_state.source_state().temperature_k(), 640.2, 2e-7);
    near(step.internal_energy_change_j, heat, 5e-5);
    assert_ne!(final_state.phase_masses_kg(), initial.phase_masses_kg());
    assert_inventory(initial);
    assert_inventory(final_state);

    for state in [&initial, &final_state] {
        let frozen = FrozenWaterCellAcoustics::new(state).unwrap();
        let equilibrium = EquilibriumWaterCellAcoustics::new(state).unwrap();
        assert!(core::ptr::eq(frozen.source(), state));
        assert!(core::ptr::eq(equilibrium.source(), state));
        assert!(equilibrium.equilibrium().is_none());
        assert_eq!(
            equilibrium.regional_equilibrium().unwrap().source(),
            &state.source_state().regional_saturation_pair().unwrap(),
        );
        for acoustic in [
            frozen.acoustic_state().unwrap(),
            equilibrium.acoustic_state().unwrap(),
        ] {
            assert_eq!(
                acoustic.density_kg_m3().to_bits(),
                state.source_state().bulk_density_kg_m3().to_bits(),
            );
        }
        assert!(equilibrium.state().sound_speed_m_s() < frozen.state().sound_speed_m_s());
        let selected = equilibrium.regional_equilibrium().unwrap().state();
        near(
            equilibrium.state().adiabatic_bulk_modulus_pa(),
            selected.adiabatic_bulk_modulus_pa(),
            2e-7,
        );
    }
}

#[test]
fn high_pressure_two_cell_transfer_retains_source_chart_and_conserves_total_inventory() {
    let mut donor = cell(0.25, 0.2);
    let mut receiver = cell(0.65, 0.15);
    let before = [donor.state(), receiver.state()];
    let transfer =
        transfer_water(&mut donor, &mut receiver, 0.0005, [0.0, 0.0], 0.1, || false).unwrap();
    assert_eq!(
        transfer.parcel.source(),
        WaterParcelSource::Equilibrium(before[0].source_state())
    );
    assert_eq!(
        transfer.parcel.enthalpy_j(),
        transfer.parcel.mass_kg() * before[0].source_state().specific_enthalpy_j_kg()
    );
    near(
        donor.state().mass_kg() + receiver.state().mass_kg(),
        before[0].mass_kg() + before[1].mass_kg(),
        1e-13,
    );
    near(
        donor.state().internal_energy_j() + receiver.state().internal_energy_j(),
        before[0].internal_energy_j() + before[1].internal_energy_j(),
        1e-8,
    );
    for (after, initial) in [donor.state(), receiver.state()].into_iter().zip(before) {
        assert_inventory(after);
        assert_eq!(after.volume_m3(), initial.volume_m3());
        assert_ne!(
            after.source_state().bulk_density_kg_m3(),
            initial.source_state().bulk_density_kg_m3()
        );
        assert_eq!(
            after.source_state().law().regional_coexistence_controls(),
            Some(dense_controls())
        );
    }
}

fn volume_controls() -> WaterVolumeControls {
    WaterVolumeControls {
        maximum_step_s: 1.0,
        maximum_relative_volume_change: 0.1,
        maximum_temperature_increment_k: 2.0,
        energy_tolerance_j: 5e-5,
        maximum_iterations: 96,
    }
}

#[test]
fn high_pressure_phase_expansion_and_return_use_raw_energy_and_reciprocal_pressure_work() {
    let mut water = cell(0.4, 0.2);
    let initial = water.state();
    let pair = initial.source_state().regional_saturation_pair().unwrap();
    let quality_change = 0.002;
    let delta_volume = initial.mass_kg()
        * quality_change
        * (pair.vapor().specific_volume().value() - pair.liquid().specific_volume().value());
    let delta_u = initial.mass_kg()
        * quality_change
        * (pair.vapor().specific_internal_energy().value()
            - pair.liquid().specific_internal_energy().value());
    let work = pair.pressure().value() * delta_volume;
    // The original dense roots need not have exactly the declared saturation
    // pressure. Use raw Delta U plus boundary work, never silently replace it
    // with an independently fitted latent enthalpy.
    let heat = delta_u + work;
    let step = water
        .advance_volume(
            initial.volume_m3() + delta_volume,
            heat,
            0.5,
            volume_controls(),
            || false,
        )
        .unwrap();
    near(step.water.state.source_state().temperature_k(), 640.0, 2e-7);
    near(step.water.boundary_work_j, work, 5e-5);
    near(step.water.internal_energy_change_j, delta_u, 5e-5);
    near(
        step.water.state.phase_masses_kg()[1],
        initial.phase_masses_kg()[1] + initial.mass_kg() * quality_change,
        5e-11,
    );
    near(step.water.energy_residual_j, 0.0, 5e-5);
    assert_inventory(water.state());
    let returned = water
        .advance_volume(initial.volume_m3(), -heat, 1.0, volume_controls(), || false)
        .unwrap();
    near(returned.water.boundary_work_j, -work, 5e-5);
    near(
        water.state().internal_energy_j(),
        initial.internal_energy_j(),
        5e-5,
    );
    near(water.state().source_state().temperature_k(), 640.0, 2e-7);
    assert_inventory(water.state());
}

#[test]
fn original_source_controls_domain_and_late_cancellation_preserve_the_cell() {
    let original = cell(0.4, 0.2);
    let mut incompatible = controls();
    incompatible.equilibrium.maximum_temperature_k = 641.0;
    assert!(
        EquilibriumWaterCell::from_isochoric_state(
            original.state().mass_kg(),
            original.state().source_state(),
            0.0,
            incompatible,
            || false,
        )
        .is_err()
    );

    let mut probe = original;
    let mut count = 0;
    probe
        .advance_volume(
            original.state().volume_m3() * 1.001,
            0.0,
            0.1,
            volume_controls(),
            || {
                count += 1;
                false
            },
        )
        .unwrap();
    let mut cancelled = original;
    let mut polls = 0;
    let result = cancelled.advance_volume(
        original.state().volume_m3() * 1.001,
        0.0,
        0.1,
        volume_controls(),
        || {
            polls += 1;
            polls == count
        },
    );
    assert!(matches!(
        result,
        Err(WaterVolumeError::Water(WaterCellError::Cancelled))
    ));
    assert_eq!(cancelled, original);
    assert!(cancelled.advance(&[], 0.0, 1e9, 0.1, || false).is_err());
    assert_eq!(cancelled, original);
}

#[test]
fn high_pressure_finite_wall_heat_uses_the_same_regional_inventory() {
    use super::thermal::{HeatExchangeSolve, WaterWallControls, advance_water_wall};
    use fs_blake3::ContentHash;
    use fs_material::phase::{
        EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, UniformEnthalpyStepInput,
    };
    use fs_material::thermal::conduction::{
        ConductivityCurve, LayerConductionConfig, LayerStackConduction,
    };

    let mut water = cell(0.4, 0.2);
    let initial = water.state();
    let target = initial
        .source_state()
        .law()
        .state_at_temperature(Temperature::new(640.02), || false)
        .unwrap();
    let heat = initial.mass_kg()
        * (target.specific_internal_energy_j_kg()
            - initial.source_state().specific_internal_energy_j_kg());
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xd3; 32]),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2700.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 50_000.0,
                temperature_k: 700.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2700.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 150_000.0,
                temperature_k: 700.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 2700.0,
            },
        ],
    )
    .unwrap();
    let mut wall = curve
        .state_at_specific_enthalpy(500.0 * (645.0 - 600.0))
        .unwrap();
    let initial_wall = wall;
    let final_wall_temperature = 645.0 - heat / (0.4 * 500.0);
    let duration = heat / (10.0 * (final_wall_temperature - 640.02));
    assert!(duration > 0.0 && duration < 10.0);
    let conductivity = [ConductivityCurve::new(vec![[600.0, 2.0], [700.0, 2.0]]).unwrap()];
    let lengths = [0.001];
    let layers = LayerStackConduction::new(&conductivity, &lengths, &[], 0.005).unwrap();
    let report = advance_water_wall(
        &mut water,
        &mut wall,
        UniformEnthalpyStepInput {
            curve: &curve,
            initial: initial_wall,
            mass_kg: 0.4,
            volume_m3: 0.4 / 2700.0,
            surface_area_m2: 0.005,
            internal_heat_j: 0.0,
            duration_s: duration,
        },
        layers,
        0.0,
        duration,
        WaterWallControls {
            exchange: HeatExchangeSolve {
                heat_bracket_j: [0.0, 2.0 * heat],
                source_tolerance_j: 5e-5,
                energy_tolerance_j: 1e-4,
                maximum_iterations: 96,
            },
            conduction: LayerConductionConfig {
                max_iterations: 96,
                absolute_heat_tolerance_w: 1e-10,
                relative_heat_tolerance: 1e-12,
            },
            maximum_step_s: 10.0,
            maximum_temperature_increment_k: [5.0, 1.0],
        },
        || false,
    )
    .unwrap();
    near(report.wall_to_water_heat_j, heat, 1e-4);
    near(water.state().source_state().temperature_k(), 640.02, 2e-7);
    near(wall.temperature_k(), final_wall_temperature, 1e-6);
    near(report.combined_energy_residual_j, 0.0, 1e-4);
    assert_eq!(water.state().mass_kg(), initial.mass_kg());
    assert_eq!(water.state().volume_m3(), initial.volume_m3());
    assert_inventory(water.state());
}
