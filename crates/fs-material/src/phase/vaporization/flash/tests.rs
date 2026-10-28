use super::*;
use fs_qty::Temperature;

fn law(temperature: f64) -> IsobaricWaterVaporization {
    IsobaricWaterVaporization::from_saturation(
        SaturatedWaterState::at_temperature(Temperature::new(temperature)).unwrap(),
    )
    .unwrap()
}

fn controls() -> WaterFlashControls {
    WaterFlashControls {
        saturation: WaterSaturationSolve {
            pressure_tolerance_pa: 1e-6,
            maximum_iterations: 128,
        },
        energy_tolerance_j_kg: 1e-8,
        entropy_tolerance_j_kg_k: 1e-10,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn saturated_hot_liquid_flashes_without_added_heat_or_shaft_work() {
    let inlet_law = law(450.0);
    let inlet = inlet_law
        .state_at_specific_enthalpy(inlet_law.saturation().liquid_specific_enthalpy().value())
        .unwrap();
    let flash = inlet
        .flash_to_pressure(Pressure::new(101_325.0), controls(), || false)
        .unwrap();
    assert_eq!(flash.inlet, inlet);
    assert_eq!(
        flash.outlet.specific_enthalpy_j_kg(),
        inlet.specific_enthalpy_j_kg()
    );
    assert_eq!(inlet.vapor_mass_fraction(), 0.0);
    assert!(flash.outlet.vapor_mass_fraction() > 0.14);
    assert!(flash.outlet.vapor_mass_fraction() < 0.15);
    assert!(flash.outlet.bulk_density_kg_m3() < 0.01 * inlet.bulk_density_kg_m3());
    assert!(flash.internal_energy_change_j_kg < 0.0);
    assert!(flash.flow_work_change_j_kg > 0.0);
    assert!(flash.entropy_generation_j_kg_k > 0.0);
    near(
        flash.first_law_residual_j_kg,
        0.0,
        controls().energy_tolerance_j_kg,
    );
    near(
        flash.outlet_pressure_residual_pa,
        0.0,
        controls().saturation.pressure_tolerance_pa,
    );
    let source =
        SaturatedWaterState::at_temperature(Temperature::new(flash.outlet.temperature_k()))
            .unwrap();
    let quality = (inlet.specific_enthalpy_j_kg() - source.liquid_specific_enthalpy().value())
        / (source.vapor_specific_enthalpy().value() - source.liquid_specific_enthalpy().value());
    near(flash.outlet.vapor_mass_fraction(), quality, 1e-15);
    let mixture_u = (1.0 - quality) * source.liquid_specific_internal_energy().value()
        + quality * source.vapor_specific_internal_energy().value();
    let mixture_v = (1.0 - quality) / source.liquid_density().value()
        + quality / source.vapor_density().value();
    near(
        mixture_u + source.pressure().value() * mixture_v,
        inlet.specific_enthalpy_j_kg(),
        1e-8,
    );
    let mixture_s = (1.0 - quality) * source.liquid_specific_entropy().value()
        + quality * source.vapor_specific_entropy().value();
    near(
        flash.entropy_generation_j_kg_k,
        mixture_s - inlet.specific_entropy_j_kg_k(),
        1e-10,
    );
}

#[test]
fn heating_then_pressure_flash_then_cooling_reuses_original_material_owners() {
    let inlet_law = law(450.0);
    let mass = 2.0;
    let liquid = inlet_law
        .state_at_specific_enthalpy(inlet_law.saturation().liquid_specific_enthalpy().value())
        .unwrap();
    let heating = inlet_law
        .advance_heat(liquid, mass, 50_000.0, 1e-7, || false)
        .unwrap();
    let flash = heating
        .state
        .flash_to_pressure(Pressure::new(101_325.0), controls(), || false)
        .unwrap();
    assert!(flash.outlet.vapor_mass_fraction() > heating.state.vapor_mass_fraction());
    near(
        mass * (flash.outlet.specific_enthalpy_j_kg() - liquid.specific_enthalpy_j_kg()),
        50_000.0,
        1e-7,
    );
    let outlet_law = flash.outlet.law();
    let cooling = outlet_law
        .advance_heat(flash.outlet, mass, -25_000.0, 1e-7, || false)
        .unwrap();
    assert!(cooling.state.vapor_mass_fraction() < flash.outlet.vapor_mass_fraction());
    assert_eq!(cooling.state.temperature_k(), flash.outlet.temperature_k());
    assert!(cooling.boundary_work_j < 0.0);
    near(
        mass * (cooling.state.specific_enthalpy_j_kg() - liquid.specific_enthalpy_j_kg()),
        25_000.0,
        1e-7,
    );
    near(
        heating.internal_energy_change_j
            + heating.boundary_work_j
            + mass * (flash.internal_energy_change_j_kg + flash.flow_work_change_j_kg)
            + cooling.internal_energy_change_j
            + cooling.boundary_work_j,
        25_000.0,
        1e-7,
    );
}

#[test]
fn unsupported_phase_domain_pressure_direction_accuracy_and_cancellation_refuse() {
    let law = law(450.0);
    let liquid = law
        .state_at_specific_enthalpy(law.saturation().liquid_specific_enthalpy().value())
        .unwrap();
    let vapor = law
        .state_at_specific_enthalpy(law.saturation().vapor_specific_enthalpy().value())
        .unwrap();
    for pressure in [
        liquid.pressure_pa(),
        liquid.pressure_pa() + 1.0,
        0.0,
        f64::NAN,
    ] {
        assert_eq!(
            liquid.flash_to_pressure(Pressure::new(pressure), controls(), || false),
            Err(WaterFlashError::PressureDropRequired)
        );
    }
    assert!(matches!(
        vapor.flash_to_pressure(Pressure::new(101_325.0), controls(), || false),
        Err(WaterFlashError::Coexistence(
            WaterVaporizationError::OutsideCoexistence
        ))
    ));
    assert!(matches!(
        liquid.flash_to_pressure(Pressure::new(100.0), controls(), || false),
        Err(WaterFlashError::Saturation(
            WaterSaturationError::PressureOutOfRange { .. }
        ))
    ));
    let mut strict = controls();
    strict.energy_tolerance_j_kg = 0.0;
    assert!(matches!(
        liquid.flash_to_pressure(Pressure::new(101_325.0), strict, || false),
        Err(WaterFlashError::EnergyBalance { .. })
    ));
    assert_eq!(
        liquid.flash_to_pressure(Pressure::new(101_325.0), controls(), || true),
        Err(WaterFlashError::Cancelled)
    );
    let mut polls = 0;
    assert!(matches!(
        liquid.flash_to_pressure(Pressure::new(101_325.0), controls(), || {
            polls += 1;
            polls == 2
        }),
        Err(WaterFlashError::Saturation(WaterSaturationError::Cancelled))
    ));
    assert_eq!(liquid.vapor_mass_fraction(), 0.0);
    assert_eq!(liquid.temperature_k(), 450.0);
}
