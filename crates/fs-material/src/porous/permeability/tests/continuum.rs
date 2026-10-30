//! Original tensor/source geometry, continuum admission and exact closures.
use super::*;

#[test]
fn continuum_uses_the_original_layer_series_and_tensor_frame() {
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let tensor = IntrinsicPermeability::orthotropic(
        [1e-12, 9e-12, 0.0],
        [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap();
    let first = PorousGasLayer::try_new(
        tensor,
        0.01,
        0.0,
        "synthetic rotated continuum permeability",
    )
    .unwrap();
    let original =
        DarcyGasBarrier::try_new(vec![first, layer(4e-12, 0.02)], [1.0, 0.0, 0.0], support())
            .unwrap();
    let geometry = original.continuum_geometry().unwrap();
    assert!(core::ptr::eq(geometry.barrier(), &original));
    near(
        geometry.geometric_resistance_per_m().unwrap(),
        0.01 / 5e-12 + 0.02 / 4e-12,
    );
    assert_eq!(geometry.support(), original.support());
    // A continuum coefficient uses the same geometry as the original no-slip
    // gas limit, while the consuming fluid supplies its own viscosity.
    let air = gas();
    let area = 0.02;
    near(
        original
            .linear_volume_resistance(&air, area)
            .unwrap()
            .unwrap(),
        air.dynamic_viscosity * geometry.geometric_resistance_per_m().unwrap() / area,
    );
    geometry
        .admit_flow(300.0, 160000.0, 140000.0, 0.02)
        .unwrap();
    geometry
        .admit_flow(300.0, 140000.0, 160000.0, -0.02)
        .unwrap();
}

#[test]
fn continuum_keeps_original_source_receipts_and_actual_state_limits() {
    let source = material(INTRINSIC_PERMEABILITY_DIMS);
    let layer = PorousGasLayer::from_material_isotropic(
        source.clone(),
        "intrinsic_permeability",
        0.01,
        None,
    )
    .unwrap();
    let original = DarcyGasBarrier::try_new(vec![layer], [1.0, 0.0, 0.0], support()).unwrap();
    let geometry = original.continuum_geometry().unwrap();
    assert_eq!(geometry.barrier().layers()[0].source(), Some(&source));
    assert_eq!(geometry.support().temperature_k, (280.0, 350.0));
    assert_eq!(geometry.support().pressure_pa, (100000.0, 200000.0));
    geometry
        .admit_flow(280.0, 100000.0, 200000.0, 20.0)
        .unwrap();
    geometry
        .admit_flow(350.0, 200000.0, 100000.0, -20.0)
        .unwrap();
    geometry.admit_temperature_speed(350.0, -20.0).unwrap();
    assert!(geometry.admit_temperature_speed(351.0, 0.0).is_err());
    assert!(geometry.admit_temperature_speed(300.0, 20.01).is_err());
    for (temperature, left, right, speed) in [
        (279.0, 150000.0, 150000.0, 0.0),
        (351.0, 150000.0, 150000.0, 0.0),
        (300.0, 99999.0, 150000.0, 0.0),
        (300.0, 150000.0, 200001.0, 0.0),
        (300.0, 150000.0, 150000.0, 20.01),
        (300.0, 150000.0, 150000.0, -20.01),
        (f64::NAN, 150000.0, 150000.0, 0.0),
        (300.0, f64::NAN, 150000.0, 0.0),
        (300.0, 150000.0, f64::INFINITY, 0.0),
        (300.0, 150000.0, 150000.0, f64::NAN),
    ] {
        assert!(
            geometry
                .admit_flow(temperature, left, right, speed)
                .is_err()
        );
    }
}

#[test]
fn a_sealed_continuum_layer_never_becomes_a_small_open_path() {
    let original = DarcyGasBarrier::try_new(
        vec![layer(1e-12, 0.01), layer(0.0, 0.001)],
        [1.0, 0.0, 0.0],
        support(),
    )
    .unwrap();
    let geometry = original.continuum_geometry().unwrap();
    assert_eq!(geometry.geometric_resistance_per_m(), None);
    geometry.admit_flow(300.0, 200000.0, 100000.0, 0.0).unwrap();
    assert!(
        geometry
            .admit_flow(300.0, 200000.0, 100000.0, f64::MIN_POSITIVE)
            .is_err()
    );
    assert!(
        geometry
            .admit_flow(300.0, 100000.0, 200000.0, -f64::MIN_POSITIVE)
            .is_err()
    );
}

#[test]
fn gas_calibrations_cannot_silently_become_liquid_mobility() {
    let model = model(vec![0.75, 0.25], 0.0);
    let original = DarcyGasBarrier::try_new_for_gas_model(
        vec![slip_layer(1e-12, 0.01, 50000.0)],
        [1.0, 0.0, 0.0],
        support(),
        &model,
    )
    .unwrap();
    assert!(original.continuum_geometry().is_err());

    use crate::porous::wet::{GasDiffusivityLaw, WetPoreGasLaw, WetPoreSupport};
    let wet = WetPoreGasLaw::try_new(
        0.0,
        0.0,
        0.5,
        GasDiffusivityLaw::ConstantTransportTortuosity { factor: 1.0 },
        WetPoreSupport {
            liquid_saturation: (0.0, 1.0),
            total_porosity: (0.5, 0.5),
            mean_gas_pore_radius_m: (1e-5, 1e-5),
        },
        "synthetic gas-only wet mobility",
    )
    .unwrap();
    for saturation in [0.0, 0.5, 1.0] {
        let state = wet.sample(saturation, 0.5, 1e-5).unwrap();
        let original = barrier(1e-12, 0.01).with_wet_states(&[state]).unwrap();
        assert!(original.continuum_geometry().is_err());
    }
}
