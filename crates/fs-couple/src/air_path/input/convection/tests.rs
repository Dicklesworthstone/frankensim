use super::*;
use fs_material::gas::{GasSpec, GasState};
use fs_material::gas_species::input::read_gas_record;

fn air(pressure: f64) -> IsobaricGasLaw {
    IsobaricGasLaw::try_from_record(GasRecord::Legacy(
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, pressure).unwrap(),
    )).unwrap()
}
fn wall(geometry: ForcedGasGeometry) -> ForcedGasBoundary {
    ForcedGasBoundary { geometry, characteristic_length_m: 0.1,
        exposed_area_m2: 0.01, speed_m_s: 1.0 }
}
fn near(a: f64, b: f64) { assert!((a-b).abs() <= 2e-12*a.abs().max(b.abs()).max(1e-12), "{a} != {b}"); }

#[test]
fn original_film_properties_and_full_temperature_drop_match_independent_formula() {
    let law = air(100000.0);
    for geometry in [ForcedGasGeometry::CylinderCrossflow, ForcedGasGeometry::LaminarFlatPlate] {
        let boundary = wall(geometry);
        let sample = boundary.evaluate(&law, 340.0).unwrap();
        let gas = sample.property_gas.state();
        assert!(sample.surface_gas.is_none());
        assert_eq!(gas.temperature, 320.0);
        assert_eq!(gas.pressure, 100000.0);
        let re = gas.density * boundary.speed_m_s * boundary.characteristic_length_m / gas.dynamic_viscosity;
        let pr = gas.dynamic_viscosity * gas.specific_heat_cp / gas.thermal_conductivity;
        let nu = match geometry {
            ForcedGasGeometry::LaminarFlatPlate => 0.664*re.sqrt()*pr.powf(1.0/3.0),
            ForcedGasGeometry::CylinderCrossflow => 0.3 + 0.62*re.sqrt()*pr.powf(1.0/3.0)
                / (1.0+(0.4/pr).powf(2.0/3.0)).powf(0.25)
                * (1.0+(re/282000.0).powf(5.0/8.0)).powf(4.0/5.0),
            ForcedGasGeometry::Sphere => unreachable!("sphere uses bulk/surface properties"),
        };
        near(sample.nusselt.evidence().value, nu);
        near(sample.conductance_w_k, nu*gas.thermal_conductivity/0.1*0.01);
        near(sample.heat_into_wall_w, -40.0*sample.conductance_w_k);
        near(sample.entropy_generation_w_k, sample.heat_into_wall_w*(1.0/340.0-1.0/300.0));
        assert_eq!(law.reference().state().temperature, 300.0);
    }
}

#[test]
fn pressure_and_composition_substitutions_change_the_same_boundary() {
    let boundary = wall(ForcedGasGeometry::LaminarFlatPlate);
    let a = boundary.evaluate(&air(100000.0), 320.0).unwrap();
    let b = boundary.evaluate(&air(200000.0), 320.0).unwrap();
    near(b.conductance_w_k/a.conductance_w_k, 2.0_f64.sqrt());
    let gas = read_gas_record(&mut "mixture 300 100000 mole 1 He 1".split_ascii_whitespace(), None, None).unwrap();
    let helium = IsobaricGasLaw::try_from_record(gas).unwrap();
    let c = boundary.evaluate(&helium, 320.0).unwrap();
    assert!((c.conductance_w_k/a.conductance_w_k-1.0).abs() > 0.1);
    assert_eq!(c.property_gas.state().specific_gas_constant, helium.reference().state().specific_gas_constant);
}

#[test]
fn humid_composition_and_passive_hot_cold_signs_are_retained() {
    let reference = GasState::try_new_moist_air(300.0, 100000.0, 0.25).unwrap();
    let law = IsobaricGasLaw::try_from_record(GasRecord::Legacy(reference)).unwrap();
    let boundary = wall(ForcedGasGeometry::CylinderCrossflow);
    for temperature in [290.0, 300.0, 340.0] {
        let sample = boundary.evaluate(&law, temperature).unwrap();
        near(sample.property_gas.state().water_mole_fraction, reference.water_mole_fraction);
        assert!(sample.heat_into_wall_w*(300.0-temperature) >= 0.0);
        assert!(sample.entropy_generation_w_k >= 0.0);
        if temperature == 300.0 { assert_eq!(sample.heat_into_wall_w, 0.0); }
    }
}

#[test]
fn malformed_geometry_source_domain_and_correlation_regime_refuse() {
    let law = air(100000.0);
    let boundary = wall(ForcedGasGeometry::LaminarFlatPlate);
    for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(ForcedGasBoundary { exposed_area_m2: value, ..boundary }.evaluate(&law, 320.0).is_err());
        assert!(ForcedGasBoundary { speed_m_s: value, ..boundary }.evaluate(&law, 320.0).is_err());
        assert!(ForcedGasBoundary { characteristic_length_m: value, ..boundary }.evaluate(&law, 320.0).is_err());
        assert!(boundary.evaluate(&law, value).is_err());
    }
    assert!(ForcedGasBoundary { speed_m_s: 1e8, ..boundary }.evaluate(&law, 320.0).is_err());
    assert!(boundary.evaluate(&law, 1e9).is_err());
    let snapshot = *law.reference().state();
    assert_eq!(law.reference().state(), &snapshot);
}

fn nitrogen() -> IsobaricGasLaw {
    let gas = read_gas_record(
        &mut "mixture 300 100000 mole 1 N2 1".split_ascii_whitespace(), None, None,
    ).unwrap();
    IsobaricGasLaw::try_from_record(gas).unwrap()
}

#[test]
fn g1_sphere_uses_bulk_transport_and_source_owned_surface_viscosity() {
    let law = nitrogen();
    let boundary = ForcedGasBoundary {
        exposed_area_m2: core::f64::consts::PI * 0.1 * 0.1,
        ..wall(ForcedGasGeometry::Sphere)
    };
    let sample = boundary.evaluate(&law, 280.0).unwrap();
    let bulk = sample.property_gas.state();
    let surface = sample.surface_gas.as_ref().unwrap().state();
    assert_eq!(&sample.property_gas, law.reference());
    assert_eq!(surface.temperature, 280.0);
    assert_eq!(surface.pressure, bulk.pressure);
    assert_eq!(surface.specific_gas_constant, bulk.specific_gas_constant);
    assert_eq!(sample.nusselt.card().id, CorrelationId::WhitakerSphere);
    let re = bulk.density * boundary.speed_m_s * boundary.characteristic_length_m
        / bulk.dynamic_viscosity;
    let pr = bulk.dynamic_viscosity * bulk.specific_heat_cp / bulk.thermal_conductivity;
    let viscosity_ratio = bulk.dynamic_viscosity / surface.dynamic_viscosity;
    let nu = 2.0 + (0.4 * re.sqrt() + 0.06 * re.powf(2.0 / 3.0))
        * pr.powf(0.4) * viscosity_ratio.powf(0.25);
    near(sample.nusselt.groups()["Re"], re);
    near(sample.nusselt.groups()["mu_bulk_over_surface"], viscosity_ratio);
    near(sample.nusselt.evidence().value, nu);
    near(sample.conductance_w_k, nu * bulk.thermal_conductivity
        / boundary.characteristic_length_m * boundary.exposed_area_m2);
    near(sample.heat_into_wall_w, 20.0 * sample.conductance_w_k);
    near(sample.entropy_generation_w_k,
        sample.heat_into_wall_w * (1.0 / 280.0 - 1.0 / 300.0));
    assert!(sample.heat_into_wall_w > 0.0);
    assert!(sample.entropy_generation_w_k > 0.0);

    // The surface viscosity is queried again at each proposed wall state;
    // moving it cannot change the bulk properties or rewrite the source law.
    let warmer = boundary.evaluate(&law, 290.0).unwrap();
    assert_eq!(warmer.property_gas, sample.property_gas);
    assert!(warmer.coefficient.value.value() < sample.coefficient.value.value());
    assert_eq!(law.reference().state().temperature, 300.0);
}

#[test]
fn g3_sphere_equilibrium_and_unsupported_viscosity_ordering_are_explicit() {
    let law = nitrogen();
    let boundary = wall(ForcedGasGeometry::Sphere);
    let equilibrium = boundary.evaluate(&law, 300.0).unwrap();
    assert_eq!(equilibrium.nusselt.groups()["mu_bulk_over_surface"], 1.0);
    assert_eq!(equilibrium.heat_into_wall_w, 0.0);
    assert_eq!(equilibrium.entropy_generation_w_k, 0.0);
    assert!(matches!(boundary.evaluate(&law, 320.0),
        Err(GasConvectionError::Correlation(_))));
    // The USSA dry-air fit has Pr < 0.71 at 300 K. The source card owns that
    // boundary; this adapter must not replace its Pr with a nominal 0.71.
    assert!(matches!(boundary.evaluate(&air(100000.0), 280.0),
        Err(GasConvectionError::Correlation(_))));
    assert!(ForcedGasBoundary { speed_m_s: 1e-8, ..boundary }
        .evaluate(&law, 280.0).is_err());
    assert!(ForcedGasBoundary { exposed_area_m2: 0.04, ..boundary }
        .evaluate(&law, 280.0).is_err());
    assert_eq!(law.reference().state().temperature, 300.0);
}
