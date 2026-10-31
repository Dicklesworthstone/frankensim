use super::*;
use crate::porous::wet::liquid::{WetPoreLiquidLaw, WetPoreLiquidState, WetPoreLiquidSupport};

fn liquid_state(saturation: f64) -> WetPoreLiquidState {
    WetPoreLiquidLaw::try_new(
        0.0,
        0.0,
        0.5,
        0.5,
        1000.0,
        WetPoreLiquidSupport {
            liquid_saturation: (0.0, 1.0),
            temperature_k: (290.0, 330.0),
            liquid_pressure_pa: (120000.0, 180000.0),
            maximum_superficial_speed_m_s: 0.25,
        },
        "synthetic liquid retention calibration",
    )
    .unwrap()
    .sample(saturation)
    .unwrap()
}

#[test]
fn liquid_mobility_preserves_original_rotated_layers_and_rebinding() {
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let tensor = IntrinsicPermeability::orthotropic(
        [1e-12, 9e-12, 0.0],
        [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap();
    let first = PorousGasLayer::try_new(tensor, 0.01, 0.0, "synthetic anisotropic solid").unwrap();
    let original =
        DarcyGasBarrier::try_new(vec![first, layer(4e-12, 0.02)], [1.0, 0.0, 0.0], support())
            .unwrap();
    let states = [liquid_state(0.6), liquid_state(0.8)];
    let geometry = original
        .continuum_geometry()
        .unwrap()
        .with_liquid_states(&states)
        .unwrap();
    let kr0 = 0.6_f64.sqrt() * 0.04;
    let kr1 = 0.8_f64.sqrt() * 0.16;
    near(
        geometry.geometric_resistance_per_m().unwrap(),
        0.01 / 5e-12 / kr0 + 0.02 / 4e-12 / kr1,
    );
    assert!(core::ptr::eq(geometry.barrier(), &original));
    assert_eq!(geometry.liquid_states(), Some(states.as_slice()));
    let confined = original
        .with_flow_constraint(DarcyFlowConstraint::ConfinedNormalFlow)
        .unwrap();
    let confined_geometry = confined
        .continuum_geometry()
        .unwrap()
        .with_liquid_states(&states)
        .unwrap();
    near(
        confined_geometry.geometric_resistance_per_m().unwrap(),
        0.01 / 1.8e-12 / kr0 + 0.02 / 4e-12 / kr1,
    );
    let flooded = [liquid_state(1.0), liquid_state(1.0)];
    let rebound = geometry.with_liquid_states(&flooded).unwrap();
    assert_eq!(
        rebound.geometric_resistance_per_m(),
        original.geometric_resistance_per_m()
    );
    assert_eq!(geometry.with_liquid_states(&states).unwrap(), geometry);
}

#[test]
fn liquid_support_and_exact_dry_seal_keep_original_source() {
    let source = material(INTRINSIC_PERMEABILITY_DIMS);
    let original = DarcyGasBarrier::try_new(
        vec![
            PorousGasLayer::from_material_isotropic(
                source.clone(),
                "intrinsic_permeability",
                0.01,
                None,
            )
            .unwrap(),
        ],
        [1.0, 0.0, 0.0],
        support(),
    )
    .unwrap();
    let states = [liquid_state(0.6)];
    let geometry = original
        .continuum_geometry()
        .unwrap()
        .with_liquid_states(&states)
        .unwrap();
    assert_eq!(geometry.barrier().layers()[0].source(), Some(&source));
    assert_eq!(geometry.support().temperature_k, (290.0, 330.0));
    assert_eq!(geometry.support().pressure_pa, (120000.0, 180000.0));
    geometry
        .admit_flow(300.0, 150000.0, 140000.0, -0.2)
        .unwrap();
    for (t, left, right, speed) in [
        (289.0, 150000.0, 140000.0, 0.1),
        (300.0, 119999.0, 140000.0, 0.1),
        (300.0, 150000.0, 180001.0, 0.1),
        (300.0, 150000.0, 140000.0, 0.251),
    ] {
        assert!(geometry.admit_flow(t, left, right, speed).is_err());
    }
    let dry = [liquid_state(0.0)];
    let sealed = geometry.with_liquid_states(&dry).unwrap();
    assert_eq!(sealed.geometric_resistance_per_m(), None);
    sealed.admit_flow(300.0, 150000.0, 140000.0, 0.0).unwrap();
    assert!(
        sealed
            .admit_flow(300.0, 150000.0, 140000.0, f64::MIN_POSITIVE)
            .is_err()
    );
    assert!(geometry.with_liquid_states(&[]).is_err());
    let intrinsic_seal = barrier(0.0, 0.01);
    assert_eq!(
        intrinsic_seal
            .continuum_geometry()
            .unwrap()
            .with_liquid_states(&states)
            .unwrap()
            .geometric_resistance_per_m(),
        None
    );
}

#[test]
fn saturation_does_not_infer_an_inertial_calibration() {
    let inertial = ForchheimerCalibration::try_new(
        2e6,
        [1.0, 0.0, 0.0],
        DarcyFlowConstraint::NormalPressureGradient,
        support(),
        "synthetic saturated inertial drag",
    )
    .unwrap();
    let layer = layer(1e-12, 0.01)
        .with_forchheimer_calibration(inertial)
        .unwrap();
    let original = DarcyGasBarrier::try_new(vec![layer], [1.0, 0.0, 0.0], support()).unwrap();
    for saturation in [0.0, 0.5, 1.0] {
        let states = [liquid_state(saturation)];
        assert!(
            original
                .continuum_geometry()
                .unwrap()
                .with_liquid_states(&states)
                .is_err()
        );
    }
}
