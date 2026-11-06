//! G0/G1/G3 combined section-stress oracles with explicitly synthetic sources.
//! Extreme arithmetic cases test representation, not a physical strain envelope.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, IntegratedIsotropicThermalExpansion, MaterialPropertySelection,
    ScalarAdmissibility, ScalarPropertyRequirement, integrate_isotropic_thermal_expansion,
    resolve_material_state_point,
};
use fs_qty::{Density, Pressure, Temperature};

const LIMIT: &str = "declared_tensile_strength";
const SHEAR: &str = "shear_modulus";

fn card(name: &str, young: f64, shear: f64, expansion: f64, limit: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, 8000.0),
        ("young_modulus", Pressure::DIMS, young),
        (SHEAR, Pressure::DIMS, shear),
        (
            "linear_thermal_expansion_coefficient",
            INVERSE_TEMPERATURE_DIMS,
            expansion,
        ),
        (LIMIT, Pressure::DIMS, limit),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if name == LIMIT {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, value), (500.0, 0.5 * value)],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0),
                interpolation: if name == LIMIT {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic combined section strength oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic circular layers".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn layer(
    card: &MaterialCard,
    radius: f64,
    temperature: f64,
) -> (ConcentricMaterialLayer, IntegratedIsotropicThermalExpansion) {
    let point = QueryPoint::new().with("T", temperature).unwrap();
    let requirements = [
        ("density", Density::DIMS),
        ("young_modulus", Pressure::DIMS),
        (SHEAR, Pressure::DIMS),
        (LIMIT, Pressure::DIMS),
    ]
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    });
    let material = resolve_material_state_point(
        card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let path = integrate_isotropic_thermal_expansion(
        card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (
        ConcentricMaterialLayer {
            outer_radius_m: radius,
            material,
        },
        path,
    )
}

fn homogeneous() -> ConcentricShearSection {
    let source = card("homogeneous annulus", 200e9, 80e9, 1e-5, 300e6);
    ConcentricElasticSection::try_new(0.0005, vec![layer(&source, 0.002, 300.0).0])
        .unwrap()
        .with_shear_modulus(SHEAR, 0.8)
        .unwrap()
}

fn thermal(temperature: f64) -> (ConcentricThermoelasticSection, ConcentricShearSection) {
    let (a, pa) = layer(&card("core", 200e9, 80e9, 1e-5, 200e6), 0.001, temperature);
    let (b, pb) = layer(&card("sleeve", 70e9, 20e9, 2e-5, 100e6), 0.002, temperature);
    let section = ConcentricElasticSection::try_new(0.0, vec![a, b]).unwrap();
    let shear = section.clone().with_shear_modulus(SHEAR, 0.8).unwrap();
    (
        ConcentricThermoelasticSection::try_new(section, vec![pa, pb]).unwrap(),
        shear,
    )
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-13 * expected.abs().max(1e-300),
        "actual {actual:e}, expected {expected:e}"
    );
}

fn ordinary_principal(normal: f64, shear: f64) -> f64 {
    0.5 * normal + (0.25 * normal * normal + shear * shear).sqrt()
}

#[test]
fn g1_combined_annular_loads_include_torsional_principal_tension_and_sign_invariance() {
    let shear = homogeneous();
    let section = shear.elastic();
    let capacity = section.axial_tensile_capacity(LIMIT).unwrap();
    let force = 60e6 * section.area_m2();
    let moment = 90e6 * section.second_moment_m4() / section.outer_radius_m();
    let torque = 240e6 * (2.0 * section.second_moment_m4()) / section.outer_radius_m();
    let axial = capacity.assess_axial_bending(force, [moment, 0.0]).unwrap();
    assert!(!axial.exceeds_limit());
    for (bending, torsion) in [
        ([0.6 * moment, 0.8 * moment], torque),
        ([-moment, 0.0], -torque),
    ] {
        let result = capacity
            .assess_axial_bending_torsion(&shear, force, bending, torsion)
            .unwrap();
        let layer = result.layers()[0];
        near(layer.minimum_axial_stress_pa, -30e6);
        near(layer.maximum_axial_stress_pa, 150e6);
        near(layer.maximum_torsional_shear_pa, 240e6);
        near(
            layer.principal_tensile_stress_pa,
            ordinary_principal(150e6, 240e6),
        );
        near(
            result.utilization(),
            ordinary_principal(150e6, 240e6) / 300e6,
        );
        assert!(result.reaches_limit());
        assert!(result.exceeds_limit());
    }
    let zero_torque = capacity
        .assess_axial_bending_torsion(&shear, force, [moment, 0.0], 0.0)
        .unwrap();
    near(zero_torque.utilization(), axial.utilization());
    near(
        zero_torque.layers()[0].principal_tensile_stress_pa,
        axial.layers()[0].maximum_axial_stress_pa,
    );
    for sign in [-1.0, 1.0] {
        let pure = capacity
            .assess_axial_bending_torsion(&shear, 0.0, [0.0; 2], sign * torque)
            .unwrap();
        near(pure.layers()[0].principal_tensile_stress_pa, 240e6);
        near(pure.utilization(), 0.8);
    }
}

#[test]
fn g0_dominant_compression_retains_small_torsional_tension_without_dimensional_squares() {
    let shear = homogeneous();
    let section = shear.elastic();
    let capacity = section.axial_tensile_capacity(LIMIT).unwrap();
    for (normal, tau, expected) in [
        (-1e200, 1e90, 1e-20),
        (1e200, 1e200, (0.5 + 1.25_f64.sqrt()) * 1e200),
        (-1e-100, 1e-100, (-0.5 + 1.25_f64.sqrt()) * 1e-100),
    ] {
        let force = normal * section.area_m2();
        let torque = tau * (2.0 * section.second_moment_m4()) / section.outer_radius_m();
        let result = capacity
            .assess_axial_bending_torsion(&shear, force, [0.0; 2], torque)
            .unwrap();
        near(result.layers()[0].principal_tensile_stress_pa, expected);
        assert!(result.utilization() > 0.0);
    }
    let compressed_force = -1e200 * section.area_m2();
    let no_torsion = capacity
        .assess_axial_bending_torsion(&shear, compressed_force, [0.0; 2], 0.0)
        .unwrap();
    assert_eq!(no_torsion.utilization(), 0.0);
    let tiny_torque = 1e-100 * (2.0 * section.second_moment_m4()) / section.outer_radius_m();
    // True positive principal tension is about 1e-400 Pa: do not report zero.
    assert!(
        capacity
            .assess_axial_bending_torsion(&shear, compressed_force, [0.0; 2], tiny_torque)
            .is_err()
    );
}

#[test]
fn g1_each_annulus_keeps_its_own_shear_modulus_radius_and_strength() {
    let core = card("weak core", 200e9, 80e9, 1e-5, 50e6);
    let sleeve = card("sleeve", 70e9, 20e9, 1e-5, 80e6);
    let section = ConcentricElasticSection::try_new(
        0.0,
        vec![layer(&core, 0.001, 300.0).0, layer(&sleeve, 0.002, 300.0).0],
    )
    .unwrap();
    let shear = section.clone().with_shear_modulus(SHEAR, 0.4).unwrap();
    let capacity = section.axial_tensile_capacity(LIMIT).unwrap();
    let force = 0.0001 * section.axial_stiffness_n();
    let moment = 0.1 * section.bending_stiffness_n_m2();
    let torque = shear.torsional_stiffness_n_m2(); // one radian per meter
    let result = capacity
        .assess_axial_bending_torsion(&shear, force, [moment, 0.0], torque)
        .unwrap();
    for ((layer, normal), tau) in result.layers().iter().zip([40e6, 21e6]).zip([80e6, 40e6]) {
        near(layer.maximum_axial_stress_pa, normal);
        near(layer.maximum_torsional_shear_pa, tau);
        near(
            layer.principal_tensile_stress_pa,
            ordinary_principal(normal, tau),
        );
    }
    assert_eq!(result.governing_layer(), 0); // internal core surface governs
    assert!(result.exceeds_limit());
    near(
        result.layers()[0].tensile_utilization,
        ordinary_principal(40e6, 80e6) / 50e6,
    );
    near(
        result.layers()[1].tensile_utilization,
        ordinary_principal(21e6, 40e6) / 80e6,
    );
    let another_correction = section.with_shear_modulus(SHEAR, 1.0).unwrap();
    assert_eq!(
        result,
        capacity
            .assess_axial_bending_torsion(&another_correction, force, [moment, 0.0], torque)
            .unwrap()
    );
}

#[test]
fn g1_zero_resultant_thermal_mismatch_combines_with_bending_and_torsion() {
    let (thermal, shear) = thermal(400.0);
    let capacity = thermal.section().axial_tensile_capacity(LIMIT).unwrap();
    let strain = thermal.zero_resultant_axial_strain();
    let native = thermal.response_at_axial_strain(strain, 0.01).unwrap();
    assert!(native.force_n().abs() < 1e-9);
    let moment = 0.05 * thermal.section().bending_stiffness_n_m2();
    let torque = shear.torsional_stiffness_n_m2();
    let result = capacity
        .assess_thermoelastic_axial_bending_torsion(
            &thermal,
            &shear,
            strain,
            [moment, 0.0],
            torque,
            0.01,
        )
        .unwrap();
    assert_eq!(capacity.layer_limits_pa(), &[150e6, 75e6]);
    for (i, (modulus, radius, tau)) in [(200e9, 0.001, 80e6), (70e9, 0.002, 40e6)]
        .into_iter()
        .enumerate()
    {
        let maximum = native.layers()[i].stress_pa() + modulus * radius * 0.05;
        near(result.layers()[i].maximum_axial_stress_pa, maximum);
        near(
            result.layers()[i].principal_tensile_stress_pa,
            ordinary_principal(maximum, tau),
        );
        assert_eq!(
            &capacity.layer_limit_sources()[i],
            thermal.section().layers()[i]
                .material
                .property(LIMIT)
                .unwrap()
        );
    }
    let normal_only = capacity
        .assess_thermoelastic_axial_bending(&thermal, strain, [moment, 0.0], 0.01)
        .unwrap();
    assert!(!normal_only.exceeds_limit());
    assert!(result.exceeds_limit());
    assert_eq!(result.governing_layer(), 0);
    assert!(result.layers()[1].maximum_axial_stress_pa < 0.0);
    assert!(result.layers()[1].principal_tensile_stress_pa > 0.0);
}

#[test]
fn g0_foreign_states_nonfinite_loads_and_original_thermal_strain_limit_refuse() {
    let (cold, shear) = thermal(300.0);
    let (hot, hot_shear) = thermal(400.0);
    let capacity = cold.section().axial_tensile_capacity(LIMIT).unwrap();
    assert!(
        capacity
            .assess_axial_bending_torsion(&hot_shear, 0.0, [0.0; 2], 0.0)
            .is_err()
    );
    assert!(
        capacity
            .assess_thermoelastic_axial_bending_torsion(&hot, &shear, 0.0, [0.0; 2], 0.0, 0.01)
            .is_err()
    );
    assert!(
        capacity
            .assess_thermoelastic_axial_bending_torsion(&cold, &hot_shear, 0.0, [0.0; 2], 0.0, 0.01)
            .is_err()
    );
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            capacity
                .assess_axial_bending_torsion(&shear, bad, [0.0; 2], 0.0)
                .is_err()
        );
        assert!(
            capacity
                .assess_axial_bending_torsion(&shear, 0.0, [bad, 0.0], 0.0)
                .is_err()
        );
        assert!(
            capacity
                .assess_axial_bending_torsion(&shear, 0.0, [0.0; 2], bad)
                .is_err()
        );
    }
    let unit_curvature_moment = cold.section().bending_stiffness_n_m2();
    assert!(cold.response_at_axial_strain(0.009, 0.01).is_ok());
    assert!(
        capacity
            .assess_thermoelastic_axial_bending_torsion(
                &cold,
                &shear,
                0.009,
                [unit_curvature_moment, 0.0],
                0.0,
                0.01
            )
            .is_err()
    );
    for (strain, admitted) in [(0.009, true), (0.011, false)] {
        let torque = (strain / cold.section().outer_radius_m()) * shear.torsional_stiffness_n_m2();
        let result = capacity
            .assess_thermoelastic_axial_bending_torsion(&cold, &shear, 0.0, [0.0; 2], torque, 0.01);
        assert_eq!(result.is_ok(), admitted);
    }
    let zero = capacity
        .assess_thermoelastic_axial_bending_torsion(&cold, &shear, 0.0, [0.0; 2], 0.0, 0.01)
        .unwrap();
    assert_eq!(zero.governing_layer(), 0);
    assert_eq!(zero.utilization(), 0.0);
}
