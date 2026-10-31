//! G0/G1/G3: independent tensor inversion and original shared flow consumers.
use super::*;

fn rotated(k: f64) -> IntrinsicPermeability {
    let s = core::f64::consts::FRAC_1_SQRT_2;
    IntrinsicPermeability::orthotropic(
        [k, 9.0 * k, 0.0],
        [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap()
}

#[test]
fn confinement_uses_inverse_tensor_and_keeps_transverse_velocity_zero() {
    let k = rotated(1e-12);
    let normal = [1.0, 0.0, 0.0];
    near(k.normal_permeability(normal).unwrap(), 5e-12);
    near(k.confined_permeability(normal).unwrap(), 1.8e-12);
    // Independent 2x2 tensor: [[5,-4],[-4,5]]*1e-12. A confined
    // pressure gradient [-10,-8] gives zero transverse Darcy velocity.
    let mu = 2e-5;
    let velocity = [
        -(5e-12 * -10.0 - 4e-12 * -8.0) / mu,
        -(-4e-12 * -10.0 + 5e-12 * -8.0) / mu,
    ];
    near(
        velocity[0],
        k.confined_permeability(normal).unwrap() * 10.0 / mu,
    );
    assert_eq!(velocity[1], 0.0);
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let unrotated = IntrinsicPermeability::orthotropic(
        [1e-12, 9e-12, 0.0],
        [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    )
    .unwrap();
    near(
        k.confined_permeability(normal).unwrap(),
        unrotated.confined_permeability([s, -s, 0.0]).unwrap(),
    );
    near(k.confined_permeability([s, s, 0.0]).unwrap(), 1e-12);
    assert_eq!(
        k.confined_permeability(normal).unwrap(),
        k.confined_permeability([-1.0, 0.0, 0.0]).unwrap()
    );
}

#[test]
fn a_sealed_material_axis_blocks_confined_flow_without_a_permeability_floor() {
    let k = rotated(1e-12);
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let normal = [s, 0.0, s];
    assert!(k.normal_permeability(normal).unwrap() > 0.0);
    assert_eq!(k.confined_permeability(normal).unwrap(), 0.0);
    let original = DarcyGasBarrier::try_new(
        vec![PorousGasLayer::try_new(k, 0.01, 0.0, "synthetic rank-two tensor").unwrap()],
        normal,
        support(),
    )
    .unwrap();
    let sealed = original
        .with_flow_constraint(DarcyFlowConstraint::ConfinedNormalFlow)
        .unwrap();
    assert_eq!(sealed.geometric_resistance_per_m(), None);
    assert_eq!(
        sealed
            .mass_flow(&gas(), 200000.0, 100000.0, 0.02)
            .unwrap()
            .mass_rate_left_to_right_kg_s,
        0.0
    );
    assert!(
        sealed
            .continuum_geometry()
            .unwrap()
            .admit_flow(300.0, 200000.0, 100000.0, 1e-100)
            .is_err()
    );
    for normal in [[0., 0., 0.], [2., 0., 0.], [f64::NAN, 0., 0.]] {
        assert!(k.confined_permeability(normal).is_err());
    }
}

#[test]
fn confined_projection_does_not_overflow_reciprocal_small_permeability() {
    for scale in [1e-310, 1e-12, 1e300] {
        let projected = rotated(scale)
            .confined_permeability([1.0, 0.0, 0.0])
            .unwrap();
        assert!((projected / scale - 1.8).abs() < 1e-12);
    }
    for value in [f64::from_bits(1), f64::MIN_POSITIVE, 1.0, 1e300] {
        let projected = IntrinsicPermeability::isotropic(value)
            .unwrap()
            .confined_permeability([1.0, 0.0, 0.0])
            .unwrap();
        assert!((projected / value - 1.0).abs() < 1e-12);
    }
}

#[test]
fn the_same_confined_series_drives_compressible_gas_and_continuum_flow() {
    let first = PorousGasLayer::try_new(
        rotated(1e-12),
        0.01,
        0.0,
        "synthetic fully developed confined laminate",
    )
    .unwrap();
    let wide =
        DarcyGasBarrier::try_new(vec![first, layer(4e-12, 0.02)], [1.0, 0.0, 0.0], support())
            .unwrap();
    let confined = wide
        .with_flow_constraint(DarcyFlowConstraint::ConfinedNormalFlow)
        .unwrap();
    assert_eq!(confined.layers(), wide.layers());
    assert_eq!(
        wide.flow_constraint(),
        DarcyFlowConstraint::NormalPressureGradient
    );
    let resistance = 0.01 / 1.8e-12 + 0.02 / 4e-12;
    near(confined.geometric_resistance_per_m().unwrap(), resistance);
    let gas = gas();
    let profile = confined
        .mass_flow_profile(&gas, 200000.0, 100000.0, 0.02)
        .unwrap();
    let expected = 0.02 * 100000.0 * 150000.0
        / (resistance * gas.dynamic_viscosity * gas.specific_gas_constant * gas.temperature);
    near(profile.flow.mass_rate_left_to_right_kg_s, expected);
    near(
        profile.flow.right_volume_rate_m3_s,
        2.0 * profile.flow.left_volume_rate_m3_s,
    );
    let p = profile.interface_pressures_pa.unwrap();
    for (i, r) in [0.01 / 1.8e-12, 0.02 / 4e-12].into_iter().enumerate() {
        near(
            (p[i] - p[i + 1]) * (0.5 * p[i] + 0.5 * p[i + 1]) / r,
            expected * gas.dynamic_viscosity * gas.specific_gas_constant * gas.temperature / 0.02,
        );
    }
    let continuum = confined.continuum_geometry().unwrap();
    near(continuum.geometric_resistance_per_m().unwrap(), resistance);
    near(
        confined
            .linear_volume_resistance(&gas, 0.02)
            .unwrap()
            .unwrap(),
        gas.dynamic_viscosity * resistance / 0.02,
    );
    assert_eq!(
        continuum.barrier().flow_constraint(),
        DarcyFlowConstraint::ConfinedNormalFlow
    );
    assert_eq!(
        confined
            .with_flow_constraint(DarcyFlowConstraint::NormalPressureGradient)
            .unwrap(),
        wide
    );
}

#[test]
fn source_domains_and_wet_updates_preserve_the_chosen_constraint() {
    use crate::porous::wet::{GasDiffusivityLaw, WetPoreGasLaw, WetPoreSupport};
    let source = material(INTRINSIC_PERMEABILITY_DIMS);
    let sourced = PorousGasLayer::from_material_isotropic(
        source.clone(),
        "intrinsic_permeability",
        0.01,
        None,
    )
    .unwrap();
    let original = DarcyGasBarrier::try_new(vec![sourced], [1., 0., 0.], support()).unwrap();
    let confined = original
        .with_flow_constraint(DarcyFlowConstraint::ConfinedNormalFlow)
        .unwrap();
    assert_eq!(confined.layers()[0].source(), Some(&source));
    assert_eq!(confined.support(), original.support());
    assert!(
        confined
            .mass_flow(&gas(), 210000.0, 150000.0, 0.01)
            .is_err()
    );

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
        "synthetic gas mobility",
    )
    .unwrap();
    let state = wet.sample(0.5, 0.5, 1e-5).unwrap();
    let response = confined.with_wet_states(&[state.clone()]).unwrap();
    assert_eq!(
        response.flow_constraint(),
        DarcyFlowConstraint::ConfinedNormalFlow
    );
    near(
        response.geometric_resistance_per_m().unwrap(),
        confined.geometric_resistance_per_m().unwrap() / state.gas_relative_permeability(),
    );
    assert_eq!(response.with_wet_states(&[state]).unwrap(), response);
    assert!(response.continuum_geometry().is_err());
}
