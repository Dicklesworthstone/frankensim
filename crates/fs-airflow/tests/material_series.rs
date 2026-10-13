//! Numerical wiring regressions; these are not calibrated coolant measurements.
use fs_airflow::graph::thermal::medium::{
    FlowHeatGeometry, FlowHeatSection, FrozenFlowMedium, MAX_MATERIAL_FLOW_SECTIONS,
};
use fs_convection::{CorrelationId, ThermalDirection};
use fs_material::fluid_transport::FluidTransport;
use fs_material::gas::{ConductivityModel, GasSpec, GasState};
use fs_qty::{Area, Length, Temperature, VolumetricFlowRate};

fn gas(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, pressure).unwrap()
}
fn section(region: &str) -> FlowHeatSection<'_> {
    FlowHeatSection {
        region,
        geometry: FlowHeatGeometry {
            flow_area: Area::new(1e-4), characteristic_length: Length::new(0.01),
            heated_length: Length::new(2.0), wetted_area: Area::new(0.002), aspect_ratio: None,
        },
        correlation: CorrelationId::CircularDuctLaminarCwt,
        direction: ThermalDirection::HeatingFluid,
    }
}
fn close(a: f64, b: f64) {
    assert!((a-b).abs() < 2e-11 * a.abs().max(b.abs()).max(1.0), "{a} != {b}");
}

#[test]
fn single_section_is_the_same_public_channel() {
    let gas = gas(100000.0);
    let medium = FrozenFlowMedium::DeclaredGas(&gas);
    let input = section("wall");
    let single = medium.channel(input.region, VolumetricFlowRate::new(1e-4), Temperature::new(300.0),
        input.geometry, input.correlation, input.direction).unwrap();
    let series = medium.channel_series(VolumetricFlowRate::new(1e-4), Temperature::new(300.0), &[input]).unwrap();
    assert_eq!(series.path(), single.path());
    assert_eq!(series.exchanges()[0].coefficient(), single.coefficient());
    assert_eq!(series.exchanges()[0].correlation(), single.correlation());
    assert_eq!(series.medium(), medium);
    assert_eq!(series.path().march(&[340.0]).unwrap(), single.path().march(&[340.0]).unwrap());
}

#[test]
fn downstream_regions_receive_upstream_heat_and_robin_rows_close_watts() {
    let gas = gas(100000.0);
    let medium = FrozenFlowMedium::DeclaredGas(&gas);
    let series = medium.channel_series(VolumetricFlowRate::new(1e-4), Temperature::new(300.0),
        &[section("upstream"), section("downstream")]).unwrap();
    let walls = [340.0, 320.0];
    let march = series.path().march(&walls).unwrap();
    let capacity = series.path().capacity_rate_w_per_k();
    let mut expected = 300.0;
    for (region, wall) in series.path().segments().iter().zip(walls) {
        expected = wall + (expected-wall) * (-region.htc_w_per_m2_k() * region.area_m2() / capacity).exp();
    }
    close(march.outlet_temperature_k, expected);
    close(march.total_heat_rate_w, capacity * (expected-300.0));
    let reset_outlet = 320.0 + (300.0-320.0) * (-series.path().segments()[1].ntu(capacity)).exp();
    assert!(march.outlet_temperature_k > reset_outlet);
    let references = march.reference_temperatures_k();
    let mut robin_heat = 0.0;
    for (i, exchange) in series.exchanges().iter().enumerate() {
        let boundary = exchange.robin_boundary(Temperature::new(references[i])).unwrap();
        close(boundary.coefficient().value.value(), series.path().segments()[i].htc_w_per_m2_k());
        robin_heat += boundary.coefficient().value.value() * exchange.geometry().wetted_area.value()
            * (walls[i] - references[i]);
    }
    close(robin_heat, march.total_heat_rate_w);
}

#[test]
fn local_geometry_changes_reynolds_and_htc_not_the_conserved_capacity() {
    let gas = gas(100000.0);
    let a = section("first");
    let mut b = section("second");
    b.geometry.flow_area = Area::new(4e-4);
    b.geometry.characteristic_length = Length::new(0.02);
    let series = FrozenFlowMedium::DeclaredGas(&gas).channel_series(VolumetricFlowRate::new(1e-4),
        Temperature::new(300.0), &[a,b]).unwrap();
    let first = &series.exchanges()[0];
    let second = &series.exchanges()[1];
    close(second.correlation().groups()["Re"], 0.5 * first.correlation().groups()["Re"]);
    close(second.coefficient().value.value(), 0.5 * first.coefficient().value.value());
    close(series.path().capacity_rate_w_per_k(), gas.density * 1e-4 * gas.specific_heat_cp);
    assert_eq!(*first.transport_state(), gas.transport_state().unwrap());
    assert_eq!(first.transport_state(), second.transport_state());
}

#[test]
fn constant_wall_refinement_preserves_outlet_and_total_heat() {
    let gas = gas(100000.0);
    let medium = FrozenFlowMedium::DeclaredGas(&gas);
    let mut base = section("whole");
    base.geometry.heated_length = Length::new(64.0);
    let whole = medium.channel_series(VolumetricFlowRate::new(1e-4), Temperature::new(300.0), &[base]).unwrap();
    let expected = whole.path().march(&[350.0]).unwrap();
    for count in [2, 8, 32] {
        let names: Vec<String> = (0..count).map(|i| format!("wall-{i}")).collect();
        let inputs: Vec<_> = names.iter().map(|name| {
            let mut part = base;
            part.region = name;
            part.geometry.heated_length = Length::new(64.0 / count as f64);
            part.geometry.wetted_area = Area::new(base.geometry.wetted_area.value() / count as f64);
            part
        }).collect();
        let split = medium.channel_series(VolumetricFlowRate::new(1e-4), Temperature::new(300.0), &inputs).unwrap();
        let actual = split.path().march(&vec![350.0; count]).unwrap();
        close(actual.outlet_temperature_k, expected.outlet_temperature_k);
        close(actual.total_heat_rate_w, expected.total_heat_rate_w);
    }
}

#[test]
fn changing_gas_data_updates_every_region_and_the_shared_capacity() {
    let air = gas(100000.0);
    let other = GasState::try_new(&GasSpec {
        molar_mass: 0.004, gamma: 5.0/3.0, sutherland_beta: 1e-6,
        sutherland_s: 80.0, conductivity: ConductivityModel::Eucken,
    }, 300.0, 100000.0).unwrap();
    let parts = [section("one"), section("two")];
    let first = FrozenFlowMedium::DeclaredGas(&air).channel_series(VolumetricFlowRate::new(1e-4), Temperature::new(300.0), &parts).unwrap();
    let second = FrozenFlowMedium::DeclaredGas(&other).channel_series(VolumetricFlowRate::new(1e-4), Temperature::new(300.0), &parts).unwrap();
    for exchange in second.exchanges() {
        assert_eq!(exchange.medium(), FrozenFlowMedium::DeclaredGas(&other));
        close(exchange.coefficient().value.value(), 3.66 * other.thermal_conductivity / 0.01);
    }
    assert_ne!(first.path().capacity_rate_w_per_k(), second.path().capacity_rate_w_per_k());
    assert_ne!(first.path().march(&[350.0,340.0]).unwrap(), second.path().march(&[350.0,340.0]).unwrap());
}

#[test]
fn malformed_or_late_refused_region_never_publishes_a_partial_path() {
    let gas = gas(100000.0);
    let medium = FrozenFlowMedium::DeclaredGas(&gas);
    let build = |parts: &[FlowHeatSection<'_>]| medium.channel_series(
        VolumetricFlowRate::new(1e-4), Temperature::new(300.0), parts);
    assert!(build(&[]).is_err());
    assert!(build(&[section("duplicate"),section("duplicate")]).is_err());
    assert!(build(&vec![section("x"); MAX_MATERIAL_FLOW_SECTIONS + 1]).is_err());
    let good = [section("one"),section("two")];
    let before = build(&good).unwrap();
    let mut bad = good;
    bad[1].correlation = CorrelationId::Gnielinski;
    assert!(build(&bad).is_err());
    bad[1] = section("two");
    bad[1].geometry.wetted_area = Area::new(f64::NAN);
    assert!(build(&bad).is_err());
    assert_eq!(build(&good).unwrap(), before);
    assert!(medium.channel_series(VolumetricFlowRate::new(1e-4), Temperature::new(0.0), &good).is_err());
}

#[test]
fn stale_derived_gas_fields_do_not_change_the_shared_correlation_inputs() {
    let mut gas = gas(100000.0);
    gas.prandtl = 123.0;
    let medium = FrozenFlowMedium::DeclaredGas(&gas);
    let series = medium.channel_series(VolumetricFlowRate::new(1e-4), Temperature::new(300.0), &[section("one")]).unwrap();
    assert_eq!(series.exchanges()[0].correlation().groups()["Pr"], gas.transport_state().unwrap().prandtl());
    assert_eq!(medium.transport_state().unwrap(), gas.transport_state().unwrap());
}
