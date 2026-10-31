//! Actual nonequilibrium populations and separate temperatures reach transport.
use super::*;
use fs_material::ionized_gas::{CarrierMomentumRelaxation, CarrierMomentumSource};
use fs_qty::{Amount, Temperature, Volume};
use fs_thermochem::electron_collision::{
    MaxwellianRateControlsV1, MaxwellianRateWeightingV1, nist,
};
use fs_thermochem::equilibrium::composition::state::FrozenIonizedGasStateV1;

#[test]
fn g1_frozen_carriers_preserve_actual_populations_and_separate_temperatures() {
    let state = frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let gas =
        IonizedGasConductivity::from_state(state.clone(), frozen_rates(&state), limits(), || false)
            .unwrap();
    assert!(gas.equilibrium().is_none());
    assert_eq!(gas.state(), &state);
    for carrier in gas.carriers() {
        relative(carrier.number_density_per_m3(), 0.001 * NA, 1e-14);
        assert_eq!(
            carrier.temperature_k(),
            if carrier.species() == &id("e-") {
                12000.0
            } else {
                6000.0
            }
        );
        relative(
            carrier.conductivity_s_m(),
            0.001 * NA * E * E / (carrier.mass_kg() * carrier.relaxation().rate_per_s()),
            1e-14,
        );
    }
    let heat = gas.joule_heating_w_m3([-2.0, 3.0]).unwrap();
    relative(
        heat.total_power_density_w_m3(),
        13.0 * gas.conductivity_s_m(),
        1e-14,
    );
    for carrier in gas.carriers() {
        let power = if carrier.is_electron() {
            heat.electron_power_density_w_m3()
        } else {
            heat.heavy_power_density_w_m3()
        };
        relative(power, 13.0 * carrier.conductivity_s_m(), 1e-14);
    }
    assert!(
        gas.joule_heating_w_m3([2.0 * gas.maximum_electric_field_v_m(), 0.0])
            .is_err()
    );
    let doubled = frozen_state(0.002, 6000.0, 12000.0, 1.0);
    let more = IonizedGasConductivity::from_state(
        doubled.clone(),
        frozen_rates(&doubled),
        limits(),
        || false,
    )
    .unwrap();
    relative(more.conductivity_s_m(), 2.0 * gas.conductivity_s_m(), 1e-14);

    let common = frozen_state(0.001, 6000.0, 6000.0, 1.0);
    let cold =
        IonizedGasConductivity::from_state(common.clone(), frozen_rates(&common), limits(), || {
            false
        })
        .unwrap();
    relative(cold.conductivity_s_m(), gas.conductivity_s_m(), 1e-14);
    let electron = |m: &IonizedGasConductivity| {
        m.carriers()
            .iter()
            .find(|c| c.species() == &id("e-"))
            .unwrap()
            .maximum_electric_field_v_m()
    };
    relative(electron(&gas), 2.0_f64.sqrt() * electron(&cold), 1e-14);
    relative(
        gas.maximum_thermal_phase_space_density(),
        cold.maximum_thermal_phase_space_density() / 2.0_f64.powf(1.5),
        1e-13,
    );
}

#[test]
fn g0_two_temperature_transport_refuses_stale_or_single_temperature_collision_sources() {
    let state = frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let p = state.pressure().value();
    for rates in [
        rates(6000.0, p, 1.0),
        vec![
            source_two_temperature("e-", 6000.0, 10000.0, p, 1e10),
            source_two_temperature("H+", 6000.0, 12000.0, p, 1e9),
        ],
        vec![
            source_two_temperature("e-", 6000.0, 12000.0, 2.0 * p, 1e10),
            source_two_temperature("H+", 6000.0, 12000.0, p, 1e9),
        ],
    ] {
        assert!(
            IonizedGasConductivity::from_state(state.clone(), rates, limits(), || false).is_err()
        );
    }
    assert_eq!(
        IonizedGasConductivity::from_state(state.clone(), frozen_rates(&state), limits(), || true),
        Err(IonizedGasError::Cancelled)
    );
}

#[test]
fn g1_cross_sections_supply_momentum_drag_at_actual_neutral_density() {
    let state = frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let rate = electron_relaxation(&state, 1e-20);
    let thermal_mean_speed = (8.0 * KB * 12000.0 / (core::f64::consts::PI * ME)).sqrt();
    let expected = (4.0 / 3.0) * (1.0 - 0.001) * NA * 1e-20 * thermal_mean_speed;
    relative(rate.rate_per_s(), expected, 2e-9);
    let CarrierMomentumSource::ElectronNeutral {
        state: retained,
        contributions,
        ..
    } = rate.source()
    else {
        panic!("original cross sections must remain attached")
    };
    assert_eq!(retained, &state);
    assert_eq!(contributions.len(), 1);
    assert_eq!(contributions[0].evaluation().target_species(), &id("H"));
    relative(contributions[0].rate_per_s(), expected, 2e-9);
    let ion = source_two_temperature("H+", 6000.0, 12000.0, state.pressure().value(), 1e9);
    let gas =
        IonizedGasConductivity::from_state(state, vec![rate, ion], limits(), || false).unwrap();
    let electron = gas
        .carriers()
        .iter()
        .find(|c| c.species() == &id("e-"))
        .unwrap();
    relative(
        electron.conductivity_s_m(),
        0.001 * NA * E * E / (ME * expected),
        2e-9,
    );
}

#[test]
fn g0_cross_section_drag_refuses_event_rates_missing_targets_and_stale_populations() {
    let state = frozen_state(0.001, 6000.0, 12000.0, 1.0);
    let evaluated = momentum_cross_section("H", 1e-20)
        .evaluate(
            state.electron_temperature(),
            MaxwellianRateWeightingV1::WeakDriftMomentum,
            MaxwellianRateControlsV1 {
                relative_tolerance: 1e-10,
                maximum_subdivisions: 1024,
            },
            || false,
        )
        .unwrap();
    for coverage in [vec![], vec![evaluated.clone(), evaluated.clone()]] {
        assert!(
            CarrierMomentumRelaxation::from_electron_neutral_cross_sections(
                &state,
                coverage,
                0.1,
                || false
            )
            .is_err()
        );
    }
    assert!(
        CarrierMomentumRelaxation::from_electron_neutral_cross_sections(
            &state,
            vec![evaluated],
            1e-4,
            || false
        )
        .is_err()
    );
    let event = nist::hydrogen(id("H"))
        .evaluate(
            state.electron_temperature(),
            MaxwellianRateWeightingV1::Event,
            MaxwellianRateControlsV1 {
                relative_tolerance: 1e-8,
                maximum_subdivisions: 1024,
            },
            || false,
        )
        .unwrap();
    assert!(
        CarrierMomentumRelaxation::from_electron_neutral_cross_sections(
            &state,
            vec![event],
            0.1,
            || false
        )
        .is_err()
    );

    // Both specimens have one total mole at the same T and V, but different
    // carrier/background fractions. An already resolved electron rate is stale.
    let owner = chemistry();
    let make = |x| {
        FrozenIonizedGasStateV1::new(
            &owner,
            &[Amount::new(1.0 - 2.0 * x), Amount::new(x), Amount::new(x)],
            Volume::new(1.0),
            Temperature::new(8000.0),
            Temperature::new(8000.0),
            1e-12,
            || false,
        )
        .unwrap()
    };
    let initial = make(0.001);
    let changed = make(0.002);
    let electron = electron_relaxation(&initial, 1e-20);
    let ion = source_two_temperature("H+", 8000.0, 8000.0, changed.pressure().value(), 1e9);
    assert!(
        IonizedGasConductivity::from_state(changed, vec![electron, ion], limits(), || false)
            .is_err()
    );
}

#[test]
fn g0_cold_light_target_cannot_bypass_the_massive_background_requirement() {
    use fs_qty::MolarMass;
    use fs_thermochem::Nasa9StandardStateModelV1;
    use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
    let original = chemistry();
    let mut models = original.models().to_vec();
    let neutral = &models[0];
    models[0] = Nasa9StandardStateModelV1::new(
        neutral.species().clone(),
        MolarMass::new(ME * NA),
        neutral.convention().clone(),
        neutral.regions().to_vec(),
    )
    .unwrap();
    let owner = IdealGasEquilibriumV1::new(
        original.elemental_matrix().clone(),
        original.charge_vector().clone(),
        models,
    )
    .unwrap();
    let state = FrozenIonizedGasStateV1::new(
        &owner,
        &[Amount::new(0.999), Amount::new(0.001), Amount::new(0.001)],
        Volume::new(1.0),
        Temperature::new(2000.0),
        Temperature::new(20000.0),
        1e-12,
        || false,
    )
    .unwrap();
    let evaluated = momentum_cross_section("H", 1e-20)
        .evaluate(
            state.electron_temperature(),
            MaxwellianRateWeightingV1::WeakDriftMomentum,
            MaxwellianRateControlsV1 {
                relative_tolerance: 1e-10,
                maximum_subdivisions: 1024,
            },
            || false,
        )
        .unwrap();
    // Thermal-speed ratio=sqrt(0.1)<0.5, but sqrt(me/M)=1 is not small.
    assert!((0.1_f64).sqrt() < 0.5);
    assert!(
        CarrierMomentumRelaxation::from_electron_neutral_cross_sections(
            &state,
            vec![evaluated],
            0.5,
            || false,
        )
        .is_err()
    );
}
