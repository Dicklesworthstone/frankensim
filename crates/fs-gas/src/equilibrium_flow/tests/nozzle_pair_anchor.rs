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
