//! The same reciprocal Coulomb source drives FEEC current and thermal exchange.
//! Field-derived pool powers are explicitly held during one frozen-population
//! energy interval; the final voltage solve is the subsequent field feedback.

mod common;
#[path = "../../fs-material/tests/common/ionized_gas.rs"]
mod fixture;

use common::{config, ends, layered, limits};
use fs_em::conduction::{ConductionControls, PlanarConduction};
use fs_gas::chemistry::energy::{ElectronHeavyEnergySolve, advance_frozen_energy_with_coulomb};
use fs_material::ionized_gas::{CarrierMomentumRelaxation, IonizedGasConductivity};
use fs_qty::{Power, Temperature, Time};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::coulomb::CoulombModelV1;
use fs_thermochem::electron_collision::{
    ElectronCrossSectionV1, MaxwellianRateControlsV1, MaxwellianRateWeightingV1,
};
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;

fn rate_controls() -> MaxwellianRateControlsV1 {
    MaxwellianRateControlsV1 {
        relative_tolerance: 1e-9,
        maximum_subdivisions: 2048,
    }
}

fn material(
    state: &FrozenIonizedGasStateV1,
    neutral: &ElectronCrossSectionV1,
    coulomb: &CoulombModelV1,
) -> IonizedGasConductivity {
    let evaluated = neutral
        .evaluate(
            state.electron_temperature(),
            MaxwellianRateWeightingV1::WeakDriftMomentum,
            rate_controls(),
            || false,
        )
        .unwrap();
    let electron = CarrierMomentumRelaxation::from_electron_neutral_cross_sections(
        state,
        vec![evaluated],
        0.1,
        || false,
    )
    .unwrap();
    let ion = fixture::source_two_temperature(
        "H+",
        state.heavy_temperature().value(),
        state.electron_temperature().value(),
        state.pressure().value(),
        1e6,
    );
    IonizedGasConductivity::from_state(
        state.clone(),
        vec![electron, ion],
        fixture::limits(),
        || false,
    )
    .unwrap()
    .with_coulomb_drag(coulomb.clone(), || false)
    .unwrap()
}

fn field_controls() -> ConductionControls {
    ConductionControls {
        linear: config(),
        current_tolerance_a: 1e-6,
        power_tolerance_w: 1e-6,
    }
}

#[test]
fn g1_reciprocal_coulomb_current_heating_and_thermal_exchange_close_one_energy_budget() {
    let initial = fixture::frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let neutral = fixture::momentum_cross_section("H", 1e-22);
    let coulomb = fixture::coulomb_model(&initial);
    let before = material(&initial, &neutral, &coulomb);
    let sigma_before = before.conductivity_s_m();
    let electron = before.carriers().iter().find(|c| c.is_electron()).unwrap();
    let background_drag =
        electron.number_density_per_m3() * electron.mass_kg() * electron.relaxation().rate_per_s();
    assert!(before.coulomb_pairs()[0].friction_coefficient_kg_m3_s() > 10.0 * background_drag);
    assert!(before.admits_electric_field([1.0, 0.0]));
    let original = PlanarConduction::new(
        layered(0.5, 1.0, 1.0),
        vec![before; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = original
        .solve(&ends(1.0), field_controls(), &mut || false)
        .unwrap();
    let resistance_before = original
        .resistance(&[2, 5], &[0, 3], field_controls(), &mut || false)
        .unwrap()
        .resistance_ohm();
    fixture::relative(field.current_density_a_m2()[0][0], -sigma_before, 1e-12);
    fixture::relative(field.power_w(), sigma_before, 1e-12);
    let mut powers = [0.0, 0.0];
    let mut correction_to_electrons = 0.0;
    let mut volume = 0.0;
    for (cell, actual_field) in field.electric_field_v_m().iter().enumerate() {
        let heat = original.conductors()[cell]
            .joule_heating_w_m3(*actual_field)
            .unwrap();
        let cell_volume = original.mesh().face_measure(cell).unwrap();
        powers[0] += cell_volume * heat.heavy_power_density_w_m3();
        powers[1] += cell_volume * heat.electron_power_density_w_m3();
        correction_to_electrons +=
            cell_volume * heat.thermal_exchange_correction_to_electrons_w_m3();
        volume += cell_volume;
        fixture::relative(
            cell_volume * heat.total_power_density_w_m3(),
            field.cell_power_w()[cell],
            2e-12,
        );
    }
    assert!(correction_to_electrons > 0.0);
    fixture::relative(volume, initial.volume().value(), 1e-14);
    fixture::relative(powers[0] + powers[1], field.power_w(), 2e-12);
    let mut controls = ElectronHeavyEnergySolve::new(
        Temperature::new(5000.0),
        Temperature::new(16000.0),
        0.2,
        0.1,
        rate_controls(),
    );
    controls.relative_tolerance = 1e-6;
    controls.initial_step_s = 1e-10;
    controls.maximum_step_s = 1e-9;
    let duration = 1e-8;
    // Powers already include the signed drift-exchange correction. The energy
    // owner separately evaluates zero-drift Coulomb AND neutral exchange;
    // neither total E.J nor the correction is supplied a second time.
    let result = advance_frozen_energy_with_coulomb(
        &initial,
        std::slice::from_ref(&neutral),
        &coulomb,
        Power::new(powers[1]),
        Power::new(powers[0]),
        Time::new(duration),
        controls,
        || false,
    )
    .unwrap();
    let next = result.state();
    assert_eq!(next.amounts(), initial.amounts());
    assert_eq!(next.volume(), initial.volume());
    assert!(next.electron_temperature().value() < initial.electron_temperature().value());
    assert!(next.heavy_temperature().value() > initial.heavy_temperature().value());
    assert_ne!(next.electron_temperature(), next.heavy_temperature());
    let balance = result.balance();
    assert!(balance.electron_to_heavy_neutral_energy_j > 0.0);
    assert!(
        balance.electron_to_heavy_coulomb_energy_j > balance.electron_to_heavy_neutral_energy_j
    );
    fixture::relative(
        balance.electron_to_heavy_energy_j,
        balance.electron_to_heavy_coulomb_energy_j + balance.electron_to_heavy_neutral_energy_j,
        1e-13,
    );
    let electron_energy = 1.5
        * R
        * initial.amounts()[2].value()
        * (next.electron_temperature().value() - initial.electron_temperature().value());
    let heavy_energy =
        1.5 * R * (next.heavy_temperature().value() - initial.heavy_temperature().value());
    assert!(
        (electron_energy - powers[1] * duration + balance.electron_to_heavy_energy_j).abs() < 1e-8
    );
    assert!(
        (heavy_energy - powers[0] * duration - balance.electron_to_heavy_energy_j).abs() < 1e-8
    );
    assert!((electron_energy + heavy_energy - field.power_w() * duration).abs() < 2e-8);
    assert!(balance.electron_energy_residual_j.abs() <= controls.energy_tolerance_j);
    assert!(balance.heavy_energy_residual_j.abs() <= controls.energy_tolerance_j);
    assert_eq!(result.coulomb_model(), Some(&coulomb));

    let retained = result.coulomb_model().unwrap();
    assert!(retained.is_bound_to_populations(next));
    let after = material(next, &neutral, retained);
    assert!(after.conductivity_s_m() < sigma_before);
    let sigma_after = after.conductivity_s_m();
    let updated = PlanarConduction::new(
        layered(0.5, 1.0, 1.0),
        vec![after; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let new_field = updated
        .solve(&ends(1.0), field_controls(), &mut || false)
        .unwrap();
    fixture::relative(new_field.current_density_a_m2()[0][0], -sigma_after, 1e-12);
    fixture::relative(new_field.power_w(), sigma_after, 1e-12);
    assert!(new_field.power_w() < field.power_w());
    let resistance_after = updated
        .resistance(&[2, 5], &[0, 3], field_controls(), &mut || false)
        .unwrap()
        .resistance_ohm();
    assert!(resistance_after > resistance_before);
}
