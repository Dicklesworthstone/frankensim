//! Actual native SP/UV source-reference integration, independent of material/runtime.
use super::*;
use fs_qty::Volume;
use fs_thermochem::equilibrium::composition::{
    isentrope::EquilibriumIsentropeV1,
    volume::{EquilibriumInternalEnergySolveV1, EquilibriumVolumeSolveV1,
        energy_path::EquilibriumEnergyVolumePathV1},
};

fn uv_controls() -> EquilibriumInternalEnergySolveV1 {
    let sp = controls().isentropic;
    EquilibriumInternalEnergySolveV1 {
        volume: EquilibriumVolumeSolveV1 {
            equilibrium: sp.equilibrium,
            relative_volume_tolerance: 1e-12,
            maximum_pressure_iterations: 80,
        },
        minimum_temperature_k: sp.minimum_temperature_k,
        maximum_temperature_k: sp.maximum_temperature_k,
        energy_tolerance_j: 1e-7,
        maximum_temperature_iterations: 80,
    }
}

#[test]
fn g1_anchored_uv_recovers_frozen_energy_volume_and_direct_entropy() {
    let chemistry = single();
    let feed = [(id("A"), Amount::new(1.0))];
    let source = EquilibriumIsentropeV1::new(
        &chemistry, Temperature::new(T0), Pressure::new(P0),
        &feed, controls().isentropic, || false,
    ).unwrap();
    let anchor = source.caloric_anchor(|| false).unwrap();
    let v0 = R * T0 / P0;
    let mut path = EquilibriumEnergyVolumePathV1::new(
        &chemistry, Temperature::new(T0), Volume::new(v0),
        &feed, uv_controls(), || false,
    ).unwrap();
    let birth = path.reference_state().clone();
    for (temperature, volume) in [(1200.0, 0.7 * v0), (800.0, 1.3 * v0)] {
        let du = 3.0 * R * (temperature - T0);
        let result = path.sample_from_caloric_anchor(
            &anchor, Volume::new(volume), Energy::new(du), || false,
        ).unwrap();
        close(result.sample.temperature.value(), temperature, 1e-7);
        close(result.sample.pressure.value(), R * temperature / volume, 2e-5);
        close(result.internal_energy_change_j, du, 1e-7);
        let ds = 3.0 * R * (temperature / T0).ln() + R * (volume / v0).ln();
        close(result.entropy_change_j_k, ds, 2e-8);
        let returned = path.publish_current(|| false).unwrap();
        let actual = EquilibriumIsentropeV1::new(
            &chemistry, returned.equilibrium().temperature(),
            returned.equilibrium().pressure(), &feed, controls().isentropic, || false,
        ).unwrap();
        close(actual.specific_entropy_change_from(&anchor, || false).unwrap() * 0.01,
            ds, 2e-8);
        assert_eq!(path.reference_state(), &birth);
    }
}

#[test]
fn g3_small_received_heat_survives_large_formation_datum_and_different_uv_birth() {
    let chemistry = owner(vec![model("A", 0.01, 4.0, 1e12, 0.0)], vec![1]);
    let feed = [(id("A"), Amount::new(1.0))];
    let source = EquilibriumIsentropeV1::new(
        &chemistry, Temperature::new(T0), Pressure::new(P0),
        &feed, controls().isentropic, || false,
    ).unwrap();
    let anchor = source.caloric_anchor(|| false).unwrap();
    let v0 = R * T0 / P0;
    let mut path = EquilibriumEnergyVolumePathV1::new(
        &chemistry, Temperature::new(1100.0), Volume::new(v0),
        &feed, uv_controls(), || false,
    ).unwrap();
    let birth = path.reference_state().clone();
    let heat = 3.0 * R * 1e-3;
    let result = path.sample_from_caloric_anchor(
        &anchor, Volume::new(v0 * 0.9), Energy::new(heat), || false,
    ).unwrap();
    close(result.sample.temperature.value(), T0 + 1e-3, 1e-8);
    close(result.internal_energy_change_j, heat, 1e-7);
    assert!(result.sample.internal_energy_change_j < -2000.0);
    assert_eq!(path.reference_state(), &birth);
}

#[test]
fn g4_uv_anchor_source_inventory_and_cancellation_refuse_stale_publication() {
    let chemistry = single();
    let feed = [(id("A"), Amount::new(1.0))];
    let v0 = R * T0 / P0;
    let mut path = EquilibriumEnergyVolumePathV1::new(
        &chemistry, Temperature::new(T0), Volume::new(v0),
        &feed, uv_controls(), || false,
    ).unwrap();
    let birth = path.reference_state().clone();
    let changed = owner(vec![model("A", 0.01, 4.1, 0.0, 0.0)], vec![1]);
    for (owner, amount) in [(&chemistry, 2.0), (&changed, 1.0)] {
        let other = EquilibriumIsentropeV1::new(
            owner, Temperature::new(T0), Pressure::new(P0),
            &[(id("A"), Amount::new(amount))], controls().isentropic, || false,
        ).unwrap();
        assert!(path.sample_from_caloric_anchor(
            &other.caloric_anchor(|| false).unwrap(), Volume::new(v0),
            Energy::new(1.0), || false,
        ).is_err());
        assert!(path.publish_current(|| false).is_err());
        assert_eq!(path.reference_state(), &birth);
    }
    let correct = EquilibriumIsentropeV1::new(
        &chemistry, Temperature::new(T0), Pressure::new(P0),
        &feed, controls().isentropic, || false,
    ).unwrap().caloric_anchor(|| false).unwrap();
    assert!(path.sample_from_caloric_anchor(
        &correct, Volume::new(v0), Energy::new(1.0), || true,
    ).is_err());
    assert!(path.publish_current(|| false).is_err());
    assert_eq!(path.reference_state(), &birth);
}

#[test]
fn g4_anchor_bridge_requires_original_source_continuity_across_separate_ranges() {
    let split = |gap: f64| {
        let original = model("A", 0.01, 4.0, 0.0, 0.0);
        let mut lower = original.regions()[0].card().clone();
        lower.validity = ValidityDomain::unconstrained().with("T", 200.0, 1000.0);
        let mut upper = lower.clone();
        upper.validity = ValidityDomain::unconstrained().with("T", 1000.0 + gap, 6000.0);
        upper.parameters.get_mut("a7").unwrap().value += 10000.0 / R;
        Nasa9StandardStateModelV1::new(
            id("A"), MolarMass::new(0.01), original.convention().clone(),
            vec![Nasa9RegionV1::from_card(lower).unwrap(),
                Nasa9RegionV1::from_card(upper).unwrap()],
        ).unwrap()
    };
    let feed = [(id("A"), Amount::new(1.0))];
    let mut low = controls().isentropic;
    low.minimum_temperature_k = 300.0;
    low.maximum_temperature_k = 900.0;
    let mut high = controls().isentropic;
    high.minimum_temperature_k = 1200.0;
    high.maximum_temperature_k = 2000.0;
    let mut uv = uv_controls();
    uv.minimum_temperature_k = high.minimum_temperature_k;
    uv.maximum_temperature_k = high.maximum_temperature_k;
    let volume = Volume::new(R * 1500.0 / P0);
    for gap in [0.0, 100.0] {
        let chemistry = owner(vec![split(gap)], vec![1]);
        let source = EquilibriumIsentropeV1::new(
            &chemistry, Temperature::new(800.0), Pressure::new(P0),
            &feed, low, || false,
        ).unwrap();
        let anchor = source.caloric_anchor(|| false).unwrap();
        let target = EquilibriumIsentropeV1::new(
            &chemistry, Temperature::new(1500.0), Pressure::new(P0),
            &feed, high, || false,
        ).unwrap();
        assert!(matches!(target.specific_enthalpy_change_from(&anchor, || false),
            Err(EquilibriumExpansionErrorV1::Invalid(_))));
        assert!(matches!(target.specific_entropy_change_from(&anchor, || false),
            Err(EquilibriumExpansionErrorV1::Invalid(_))));
        let mut path = EquilibriumEnergyVolumePathV1::new(
            &chemistry, Temperature::new(1500.0), volume, &feed, uv, || false,
        ).unwrap();
        assert!(matches!(path.sample_from_caloric_anchor(
            &anchor, volume, Energy::new(3.0 * R * 700.0), || false,
        ), Err(fs_thermochem::equilibrium::composition::volume::EquilibriumVolumeErrorV1::Response(
            EquilibriumExpansionErrorV1::Invalid(_)))));
        assert!(path.publish_current(|| false).is_err());
    }
    let corrected = split(0.0).with_caloric_continuity(
        Temperature::new(500.0),
        fs_thermochem::MolarEnergyQuantityV1::new(10001.0),
        fs_thermochem::MolarThermalQuantityV1::new(1.0),
    ).unwrap();
    let chemistry = owner(vec![corrected], vec![1]);
    let source = EquilibriumIsentropeV1::new(
        &chemistry, Temperature::new(800.0), Pressure::new(P0), &feed, low, || false,
    ).unwrap();
    let mut path = EquilibriumEnergyVolumePathV1::new(
        &chemistry, Temperature::new(1500.0), volume, &feed, uv, || false,
    ).unwrap();
    let result = path.sample_from_caloric_anchor(
        &source.caloric_anchor(|| false).unwrap(), volume,
        Energy::new(3.0 * R * 700.0), || false,
    ).unwrap();
    close(result.sample.temperature.value(), 1500.0, 1e-7);
}
