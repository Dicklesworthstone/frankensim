//! G0/G1/G3 checks of original source roles across actual section consumers.
use super::*;
use crate::conductor::{ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY, OhmicDrive};
use crate::section_temperature::{
    ConcentricSectionTemperatureLaw,
    thermoelastic::{ConcentricThermoelasticLaw, rod::ConcentricRodThermoelasticLaw},
};
use crate::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA,
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY as CP,
    ScalarAdmissibility, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY as K,
    resolve_material_state_point,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{
    Dims, Temperature, Time,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn requirement(name: &str, dims: Dims) -> ScalarPropertyRequirement {
    ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
}
fn roles(pair: u8) -> LayerIsotropicParameters {
    LayerIsotropicParameters {
        density: requirement("rho-source", Density::DIMS),
        parameters: match pair {
            0 => IsotropicElasticParameters::YoungPoisson {
                young_modulus: requirement("E-source", Pressure::DIMS),
                poisson_ratio: requirement("nu-source", Dims::NONE),
            },
            1 => IsotropicElasticParameters::YoungShear {
                young_modulus: requirement("E-source", Pressure::DIMS),
                shear_modulus: requirement("G-source", Pressure::DIMS),
            },
            _ => IsotropicElasticParameters::BulkShear {
                bulk_modulus: requirement("K-source", Pressure::DIMS),
                shear_modulus: requirement("G-source", Pressure::DIMS),
            },
        },
    }
}
fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn fixture(pair: u8, varying: bool) -> (MaterialCard, ResolvedMaterialStatePoint) {
    let mut values = vec![
        ("rho-source", Density::DIMS, [8000.0, 8000.0]),
        ("eta-b", DynViscosity::DIMS, [1e7, 1e7]),
        ("eta-a", DynViscosity::DIMS, [2e7, 2e7]),
        ("eta-s", DynViscosity::DIMS, [3e7, 3e7]),
        ("loss", Dims::NONE, [0.01, 0.01]),
        ("strength", Pressure::DIMS, [200e6, 200e6]),
        ("delta", Pressure::DIMS, [1e9, 1e9]),
        ("tau", Time::DIMS, [0.001, 0.001]),
        (ALPHA, INVERSE_TEMPERATURE_DIMS, [1e-5, 1e-5]),
        (CP, SPECIFIC_HEAT_CAPACITY_DIMS, [400.0, 800.0]),
        (K, THERMAL_CONDUCTIVITY_DIMS, [20.0, 40.0]),
        (
            ELECTRICAL_RESISTIVITY_PROPERTY,
            ELECTRICAL_RESISTIVITY_DIMS,
            [1e-7, 2e-7],
        ),
    ];
    values.extend(match pair {
        0 => vec![
            ("E-source", Pressure::DIMS, [200e9, 100e9]),
            ("nu-source", Dims::NONE, [0.25, 0.25]),
        ],
        1 => vec![
            ("E-source", Pressure::DIMS, [200e9, 100e9]),
            ("G-source", Pressure::DIMS, [80e9, 40e9]),
        ],
        _ => vec![
            ("K-source", Pressure::DIMS, [100e9, 80e9]),
            ("G-source", Pressure::DIMS, [40e9, 10e9]),
        ],
    });
    // Equal at the reference temperature but with a different thermal law:
    // this auxiliary observation cannot replace the selected shear source.
    let unrelated_shear = if pair == 2 { 40e9 } else { 80e9 };
    values.push((
        "unselected-G",
        Pressure::DIMS,
        [unrelated_shear, unrelated_shear],
    ));
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, endpoints) in values {
        let curve = varying && endpoints[0] != endpoints[1];
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if curve {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, endpoints[0]), (500.0, endpoints[1])],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar {
                        value: endpoints[0],
                        dims,
                    }
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", absolute(), 300.0, 500.0)
                    .with("omega", 1.0, 20000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if curve {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured section parameter oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(requirement(name, dims));
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: format!("parameter-pair-{pair}"),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("T", absolute(), 300.0)
        .unwrap()
        .with("omega", 500.0)
        .unwrap();
    let state = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (card, state)
}
fn section(material: ResolvedMaterialStatePoint, pair: u8) -> ConcentricElasticSection {
    ConcentricElasticSection::try_new_with_isotropic_parameters(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: 0.002,
            material,
        }],
        &[roles(pair)],
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-15),
        "{a} != {b}"
    );
}

#[test]
fn source_parameterized_layers_share_mass_stiffness_heat_and_losses() {
    let (_, kg) = fixture(2, false);
    let (_, enu) = fixture(0, false);
    let elastic = ConcentricElasticSection::try_new_with_isotropic_parameters(
        0.0,
        vec![
            ConcentricMaterialLayer {
                outer_radius_m: 0.001,
                material: kg.clone(),
            },
            ConcentricMaterialLayer {
                outer_radius_m: 0.002,
                material: enu.clone(),
            },
        ],
        &[roles(2), roles(0)],
    )
    .unwrap();
    let e = [9.0 * 100e9 * 40e9 / (3.0 * 100e9 + 40e9), 200e9];
    let a = [core::f64::consts::PI * 1e-6, core::f64::consts::PI * 3e-6];
    let i = [
        core::f64::consts::PI * 1e-12 / 4.0,
        core::f64::consts::PI * 15e-12 / 4.0,
    ];
    close(elastic.line_density_kg_m(), 8000.0 * (a[0] + a[1]), 2e-14);
    close(
        elastic.axial_stiffness_n(),
        e[0] * a[0] + e[1] * a[1],
        2e-14,
    );
    close(
        elastic.bending_stiffness_n_m2(),
        e[0] * i[0] + e[1] * i[1],
        2e-14,
    );
    assert_eq!(elastic.layers()[0].material, kg);
    assert_eq!(elastic.layers()[1].material, enu);
    assert!(kg.property(YOUNG_MODULUS_PROPERTY).is_none());
    assert!(kg.property(DENSITY_PROPERTY).is_none());
    let shear = elastic.clone().with_isotropic_shear(0.8).unwrap();
    assert_eq!(shear.shear_property(), None);
    close(
        shear.torsional_stiffness_n_m2(),
        2.0 * (40e9 * i[0] + 80e9 * i[1]),
        2e-14,
    );
    let torsion = shear.with_torsional_viscosity("eta-s").unwrap();
    close(
        torsion.rotary_inertia_kg_m(),
        16000.0 * (i[0] + i[1]),
        2e-14,
    );
    close(
        elastic
            .clone()
            .with_kelvin_voigt_bending("eta-b")
            .unwrap()
            .viscous_stiffness_n_m2_s(),
        1e7 * (i[0] + i[1]),
        2e-14,
    );
    close(
        elastic
            .clone()
            .with_axial_viscosity("eta-a")
            .unwrap()
            .viscous_rigidity_n_s(),
        2e7 * (a[0] + a[1]),
        2e-14,
    );
    close(
        elastic
            .clone()
            .with_loss_angle_bending("loss")
            .unwrap()
            .hysteretic_stiffness_n_m2(),
        0.01 * elastic.bending_stiffness_n_m2(),
        2e-14,
    );
    let capacity = elastic.axial_tensile_capacity("strength").unwrap();
    close(
        capacity.capacity_n(),
        200e6 * elastic.axial_stiffness_n() / e[1],
        2e-14,
    );
    let transport = elastic
        .thermal_transport()
        .unwrap()
        .with_ohmic_transport(1.0)
        .unwrap();
    let response = transport.evaluate(OhmicDrive::Voltage(1.0), 10.0).unwrap();
    close(
        response.layers.iter().map(|x| x.heat_capacity_j_k).sum(),
        400.0 * elastic.line_density_kg_m(),
        2e-14,
    );
    close(response.joule_power_w, (a[0] + a[1]) / 1e-7, 2e-14);
}

#[test]
fn nonlinear_kg_temperature_refresh_preserves_roles_and_auxiliary_physics() {
    let (card, source) = fixture(2, true);
    let original = section(source, 2);
    let law = ConcentricSectionTemperatureLaw::try_new(original.clone(), &[&card], &["T"]).unwrap();
    let middle = law.at_temperatures_k(&[400.0]).unwrap();
    let hot = law.at_temperatures_k(&[500.0]).unwrap();
    let expected_e = 9.0 * 90e9 * 25e9 / (3.0 * 90e9 + 25e9);
    close(
        middle.layer_integrals()[0].young_modulus_pa(),
        expected_e,
        2e-14,
    );
    assert!(
        (expected_e
            - 0.5
                * (original.layer_integrals()[0].young_modulus_pa()
                    + hot.layer_integrals()[0].young_modulus_pa()))
        .abs()
            > 1e9
    );
    assert_eq!(hot.isotropic_parameters(), original.isotropic_parameters());
    assert_eq!(hot.line_density_kg_m(), original.line_density_kg_m());
    assert!(
        hot.layers()[0]
            .material
            .property(YOUNG_MODULUS_PROPERTY)
            .is_none()
    );
    let old_electric = original
        .thermal_transport()
        .unwrap()
        .with_ohmic_transport(1.0)
        .unwrap()
        .evaluate(OhmicDrive::Voltage(1.0), 0.0)
        .unwrap();
    let hot_electric = hot
        .thermal_transport()
        .unwrap()
        .with_ohmic_transport(1.0)
        .unwrap()
        .evaluate(OhmicDrive::Voltage(1.0), 0.0)
        .unwrap();
    close(
        hot_electric.joule_power_w,
        old_electric.joule_power_w / 2.0,
        2e-14,
    );
    close(
        hot_electric.layers[0].heat_capacity_j_k,
        2.0 * old_electric.layers[0].heat_capacity_j_k,
        2e-14,
    );
    let thermal = law.thermoelastic_at_temperatures_k(&[500.0]).unwrap();
    let clamped = thermal.response_at_axial_strain(0.0, 0.1).unwrap();
    close(clamped.force_n(), -hot.axial_stiffness_n() * 0.002, 2e-14);
    assert!(
        ConcentricThermoelasticLaw::try_new(law).is_err(),
        "nonlinear derived E must not enter an affine thermal potential"
    );
}

#[test]
fn derived_shear_thermodynamics_matches_original_potential_derivatives() {
    for pair in [0, 1, 2] {
        let (card, source) = fixture(pair, pair != 2);
        let original = section(source, pair);
        let temperature =
            ConcentricSectionTemperatureLaw::try_new(original, &[&card], &["T"]).unwrap();
        let law =
            ConcentricRodThermoelasticLaw::try_new_with_isotropic_shear(temperature, 0.8).unwrap();
        let evaluate = |t| {
            law.at_strain(&[t], [0.004, 0.003, 0.002], [0.1, 0.2, 0.3], 0.1)
                .unwrap()
        };
        let state = evaluate(400.0);
        let plus = evaluate(400.01);
        let minus = evaluate(399.99);
        let layer = &state.layers[0];
        close(
            layer.entropy_j_m_k,
            -(plus.layers[0].free_energy_j_m - minus.layers[0].free_energy_j_m) / 0.02,
            2e-8,
        );
        close(
            layer.shear_modulus_derivative_pa_k,
            if pair == 2 { 0.0 } else { -2e8 },
            2e-14,
        );
        assert_eq!(state.section.shear_property(), None);
    }
}

#[test]
fn named_shear_cannot_replace_original_source_support_or_thermal_derivatives() {
    for pair in [0, 1, 2] {
        let (card, source) = fixture(pair, true);
        let original = section(source, pair);
        assert_eq!(
            original.layers()[0]
                .material
                .property("unselected-G")
                .unwrap()
                .value_si(),
            original.layer_integrals()[0]
                .isotropic_shear_modulus_pa()
                .unwrap()
        );
        assert!(
            original
                .clone()
                .with_shear_modulus("unselected-G", 0.8)
                .is_err(),
            "equal sampled G does not establish the selected source or its derivative"
        );
        let temperature =
            ConcentricSectionTemperatureLaw::try_new(original.clone(), &[&card], &["T"]).unwrap();
        assert!(
            ConcentricRodThermoelasticLaw::try_new(temperature.clone(), "unselected-G", 0.8)
                .is_err()
        );
        if pair != 0 {
            let selected = original.with_shear_modulus("G-source", 0.8).unwrap();
            assert_eq!(selected.shear_property(), Some("G-source"));
            assert_eq!(
                selected.source_properties(0).unwrap()[1]
                    .requirement()
                    .name(),
                "G-source"
            );
            if pair == 1 {
                let law =
                    ConcentricRodThermoelasticLaw::try_new(temperature, "G-source", 0.8).unwrap();
                let state = law
                    .at_strain(&[400.0], [0.004, 0.003, 0.002], [0.1, 0.2, 0.3], 0.1)
                    .unwrap();
                close(state.layers[0].shear_modulus_derivative_pa_k, -2e8, 2e-14);
            }
        }
    }
}

#[test]
fn parameterized_sections_support_source_spectra_and_thermoelastic_loss() {
    let (_, source) = fixture(1, false);
    let elastic = section(source, 1);
    let prony = elastic
        .clone()
        .with_prony_bending(&[relaxation::LayerRelaxationProperties {
            equilibrium_modulus: "E-source".into(),
            branches: vec![("delta".into(), "tau".into())],
        }])
        .unwrap();
    close(
        prony.branches()[0].relaxing_stiffness_n_m2,
        1e9 * elastic.second_moment_m4(),
        2e-14,
    );
    let torsion = elastic
        .clone()
        .with_isotropic_shear(1.0)
        .unwrap()
        .with_torsional_viscosity("eta-s")
        .unwrap()
        .with_prony_torsion(&[relaxation::LayerRelaxationProperties {
            equilibrium_modulus: "G-source".into(),
            branches: vec![("delta".into(), "tau".into())],
        }])
        .unwrap();
    close(
        torsion.branches()[0].relaxing_stiffness_n_m2,
        2e9 * elastic.second_moment_m4(),
        2e-14,
    );
    let thermal = elastic.with_thermoelastic_bending().unwrap();
    close(thermal.law().e, 200e9, 2e-14);
    close(thermal.law().rho, 8000.0, 2e-14);
}

#[test]
fn parameterized_section_refresh_never_invents_sources_or_reinterprets_geometry() {
    let (_, source) = fixture(2, false);
    let original = section(source.clone(), 2);
    assert!(
        original
            .clone()
            .with_shear_modulus("K-source", 1.0)
            .is_err()
    );
    assert_eq!(
        original
            .clone()
            .with_shear_modulus("G-source", 1.0)
            .unwrap()
            .layer_shear_modulus_pa(0)
            .unwrap(),
        40e9
    );
    assert!(ConcentricElasticSection::try_new(0.0, original.layers().to_vec()).is_err());
    assert!(
        ConcentricElasticSection::try_new_with_isotropic_parameters(
            0.0,
            original.layers().to_vec(),
            &[]
        )
        .is_err()
    );
    let mut changed = original.layers().to_vec();
    changed[0].outer_radius_m *= 2.0;
    assert!(original.with_material_layers(changed).is_err());
    assert!(original.elastic_source_properties(1).is_err());
    let mut wrong = roles(2);
    wrong.density = requirement("rho-source", Pressure::DIMS);
    assert!(
        ConcentricElasticSection::try_new_with_isotropic_parameters(
            0.0,
            original.layers().to_vec(),
            &[wrong]
        )
        .is_err()
    );
    let mut swapped = roles(2);
    swapped.parameters = IsotropicElasticParameters::BulkShear {
        bulk_modulus: requirement("G-source", Pressure::DIMS),
        shear_modulus: requirement("K-source", Pressure::DIMS),
    };
    let mixed = ConcentricElasticSection::try_new_with_isotropic_parameters(
        0.0,
        vec![
            ConcentricMaterialLayer {
                outer_radius_m: 0.001,
                material: source.clone(),
            },
            ConcentricMaterialLayer {
                outer_radius_m: 0.002,
                material: source,
            },
        ],
        &[roles(2), swapped],
    )
    .unwrap();
    assert!(
        mixed.with_thermoelastic_bending().is_err(),
        "same raw card does not imply the same physical role binding"
    );
    assert_eq!(
        original,
        original
            .with_material_layers(original.layers().to_vec())
            .unwrap()
    );
}
