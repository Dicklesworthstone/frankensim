//! G0/G1/G3 analytical section checks, using explicitly synthetic sources.
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS, IntegratedIsotropicThermalExpansion,
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    integrate_isotropic_thermal_expansion, resolve_material_state_point};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_qty::semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm};

const LIMIT: &str = "tensile_limit";

fn card(name: &str, young: f64, expansion: f64, limit: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (key, dims, value) in [
        ("density", Density::DIMS, 8000.0),
        ("young_modulus", Pressure::DIMS, young),
        ("linear_thermal_expansion_coefficient", INVERSE_TEMPERATURE_DIMS, expansion),
        (LIMIT, Pressure::DIMS, limit),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims),
            value: if key == LIMIT {
                PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, value), (500.0, 0.5 * value)], dims }
            } else { PropertyValue::Scalar { value, dims } },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if key == LIMIT { InterpolationPolicy::LinearInside }
                else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![],
            provenance: Provenance { source: "synthetic tensile section oracle".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: name.into(), phase: "solid".into(),
        process: "synthetic".into(), revision: 0 }, claims, vec![]).unwrap()
}

fn layer(card: &MaterialCard, radius: f64, t: f64)
    -> (ConcentricMaterialLayer, IntegratedIsotropicThermalExpansion) {
    let point = QueryPoint::new().with("T", t).unwrap();
    let requirements = [("density", Density::DIMS), ("young_modulus", Pressure::DIMS),
        (LIMIT, Pressure::DIMS)].map(|(name, dims)|
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap());
    let material = resolve_material_state_point(card, &point, &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let path = integrate_isotropic_thermal_expansion(card,
        &QueryPoint::new().with("T", 300.0).unwrap(), &point,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    (ConcentricMaterialLayer { outer_radius_m: radius, material }, path)
}

fn thermal(t: f64) -> ConcentricThermoelasticSection {
    let (a, pa) = layer(&card("core", 200e9, 1e-5, 100e6), 0.001, t);
    let (b, pb) = layer(&card("sleeve", 70e9, 2e-5, 100e6), 0.002, t);
    ConcentricThermoelasticSection::try_new(
        ConcentricElasticSection::try_new(0.0, vec![a, b]).unwrap(), vec![pa, pb],
    ).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 3e-13 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}");
}

#[test]
fn g1_combined_force_and_biaxial_bending_match_annular_extrema() {
    let (a, _) = layer(&card("homogeneous", 200e9, 1e-5, 300e6), 0.002, 300.0);
    let section = ConcentricElasticSection::try_new(0.0, vec![a]).unwrap();
    let capacity = section.axial_tensile_capacity(LIMIT).unwrap();
    assert_eq!(capacity.section(), &section);
    let force = 0.2 * 300e6 * section.area_m2();
    let moment = 0.9 * 300e6 * section.second_moment_m4() / section.outer_radius_m();
    let result = capacity.assess_axial_bending(force, [0.6 * moment, 0.8 * moment]).unwrap();
    close(result.layers()[0].minimum_axial_stress_pa, -0.7 * 300e6);
    close(result.layers()[0].maximum_axial_stress_pa, 1.1 * 300e6);
    close(result.utilization(), 1.1);
    assert!(result.exceeds_limit());
    // Circular geometry makes the extrema invariant to bending-axis rotation.
    let rotated = capacity.assess_axial_bending(force, [-moment, 0.0]).unwrap();
    close(rotated.utilization(), result.utilization());
    // Compression and bending can still leave tensile fibers.
    let compressed = capacity.assess_axial_bending(-force, [moment, 0.0]).unwrap();
    close(compressed.utilization(), 0.7);
    let all_compressed = capacity.assess_axial_bending(-10.0 * force, [moment, 0.0]).unwrap();
    assert_eq!(all_compressed.utilization(), 0.0);
    assert!(all_compressed.layers()[0].maximum_axial_stress_pa < 0.0);
    assert!(!all_compressed.exceeds_limit());
}

#[test]
fn g1_each_layer_keeps_its_modulus_radius_and_limit() {
    let (a, _) = layer(&card("core", 200e9, 1e-5, 200e6), 0.001, 300.0);
    let (b, _) = layer(&card("sleeve", 70e9, 1e-5, 100e6), 0.002, 300.0);
    let section = ConcentricElasticSection::try_new(0.0, vec![a, b]).unwrap();
    let capacity = section.axial_tensile_capacity(LIMIT).unwrap();
    let axial = capacity.assess_axial_bending(100.0, [0.0; 2]).unwrap();
    assert_eq!(axial.governing_layer(), capacity.governing_layer());
    close(axial.utilization(), capacity.utilization(100.0).unwrap());
    assert_eq!(axial.governing_layer(), 0);
    let moment = 0.5 * section.bending_stiffness_n_m2();
    let bent = capacity.assess_axial_bending(0.0, [moment, 0.0]).unwrap();
    close(bent.layers()[0].maximum_axial_stress_pa, 200e9 * 0.001 * 0.5);
    close(bent.layers()[1].maximum_axial_stress_pa, 70e9 * 0.002 * 0.5);
    assert_eq!(bent.governing_layer(), 1);
    close(bent.utilization(), 0.7);
}

#[test]
fn g1_zero_resultant_thermal_stress_can_exceed_tensile_limit() {
    let thermal = thermal(400.0);
    let capacity = thermal.section().axial_tensile_capacity(LIMIT).unwrap();
    let strain = (200.0 * 0.001 + 210.0 * 0.002) / 410.0;
    close(thermal.zero_resultant_axial_strain(), strain);
    let native = thermal.response_at_axial_strain(strain, 0.01).unwrap();
    assert!(native.force_n().abs() < 1e-9);
    let assessment = capacity.assess_thermoelastic_axial_bending(&thermal, strain, [0.0; 2], 0.01).unwrap();
    // Strength is queried at 400 K, not copied from the 300 K reference.
    assert_eq!(capacity.layer_limits_pa(), &[75e6; 2]);
    close(assessment.layers()[0].maximum_axial_stress_pa, 200e9 * (strain - 0.001));
    close(assessment.utilization(), 200e9 * (strain - 0.001) / 75e6);
    assert!(assessment.exceeds_limit());
    assert_eq!(assessment.governing_layer(), 0);
    assert_eq!(assessment.layers()[1].tensile_utilization, 0.0);
    for (native, result) in native.layers().iter().zip(assessment.layers()) {
        assert_eq!(result.minimum_axial_stress_pa, native.stress_pa());
        assert_eq!(result.maximum_axial_stress_pa, native.stress_pa());
    }
    let curvature = 0.1;
    let moment = curvature * thermal.section().bending_stiffness_n_m2();
    let bent = capacity.assess_thermoelastic_axial_bending(&thermal, strain, [moment, 0.0], 0.01).unwrap();
    close(bent.layers()[0].maximum_axial_stress_pa,
        native.layers()[0].stress_pa() + 200e9 * 0.001 * curvature);
}

#[test]
fn g0_thermal_assessment_refuses_foreign_or_stale_states_and_bending_strain_escape() {
    let cold = thermal(300.0);
    let hot = thermal(400.0);
    let capacity = cold.section().axial_tensile_capacity(LIMIT).unwrap();
    assert!(capacity.assess_thermoelastic_axial_bending(&hot, 0.0, [0.0; 2], 0.01).is_err());
    let mut foreign_layers = cold.section().layers().to_vec();
    foreign_layers.swap(0, 1);
    foreign_layers[0].outer_radius_m = 0.001;
    foreign_layers[1].outer_radius_m = 0.002;
    let foreign_section = ConcentricElasticSection::try_new(0.0, foreign_layers).unwrap();
    let foreign = ConcentricThermoelasticSection::try_new(foreign_section,
        cold.expansion_paths().iter().rev().cloned().collect()).unwrap();
    assert!(capacity.assess_thermoelastic_axial_bending(&foreign, 0.0, [0.0; 2], 0.01).is_err());
    let moment = cold.section().bending_stiffness_n_m2();
    // Mean strain alone is admissible; outer bending fibers are not.
    assert!(cold.response_at_axial_strain(0.009, 0.01).is_ok());
    assert!(capacity.assess_thermoelastic_axial_bending(&cold, 0.009, [moment, 0.0], 0.01).is_err());
    assert!(capacity.assess_thermoelastic_axial_bending(&cold, 0.007, [moment, 0.0], 0.01).is_ok());
}

#[test]
fn g3_bending_tensile_utilization_is_invariant_to_layer_subdivision() {
    let source = card("homogeneous", 200e9, 1e-5, 300e6);
    let whole = ConcentricElasticSection::try_new(0.0,
        vec![layer(&source, 0.002, 300.0).0]).unwrap();
    let split = ConcentricElasticSection::try_new(0.0,
        [0.0005, 0.001, 0.002].map(|r| layer(&source, r, 300.0).0).to_vec()).unwrap();
    let a = whole.axial_tensile_capacity(LIMIT).unwrap().assess_axial_bending(100.0, [0.3, 0.4]).unwrap();
    let b = split.axial_tensile_capacity(LIMIT).unwrap().assess_axial_bending(100.0, [0.3, 0.4]).unwrap();
    close(a.utilization(), b.utilization());
    close(a.layers()[0].maximum_axial_stress_pa, b.layers()[2].maximum_axial_stress_pa);
    assert_eq!(b.governing_layer(), 2);
}

fn contextual_section(axis: &str, quantity: QuantitySpec, curve: bool) -> ConcentricElasticSection {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [("density", Density::DIMS, 8000.0),
        ("young_modulus", Pressure::DIMS, 200e9), (LIMIT, Pressure::DIMS, 300e6)] {
        let bounded = name == LIMIT;
        let validity = ValidityDomain::unconstrained().with("T", 300.0, 500.0);
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: if bounded && curve {
                PropertyValue::Curve { abscissa: axis.into(), abscissa_dims: quantity.dims(),
                    knots: vec![(0.0, value), (0.2, 0.5 * value)], dims }
            } else { PropertyValue::Scalar { value, dims } },
            validity: if bounded { validity.with_quantity(axis, quantity, 0.0, 0.2) } else { validity },
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if bounded && curve { InterpolationPolicy::LinearInside }
                else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![], provenance: Provenance { source: "synthetic strength context".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(), phase: "solid".into(),
        process: "held state".into(), revision: 0 }, claims, vec![]).unwrap();
    let point = QueryPoint::new().with("T", 300.0).unwrap().with_quantity(axis, quantity, 0.1).unwrap();
    let material = resolve_material_state_point(&card, &point, &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer { outer_radius_m: 0.002, material }]).unwrap()
}

#[test]
fn g0_fixed_moisture_strength_is_admitted_but_unowned_axes_are_not() {
    let moisture = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Moisture(MoistureBasis::DryMass), ValueForm::Static));
    let humid = contextual_section("moisture", moisture, true);
    let capacity = humid.axial_tensile_capacity(LIMIT).unwrap();
    assert_eq!(capacity.layer_limits_pa(), &[225e6]);
    assert_eq!(capacity.layer_limit_sources()[0], *humid.layers()[0].material.property(LIMIT).unwrap());
    let external = humid.axial_tensile_capacity_from_resolved(
        LIMIT, capacity.layer_limit_sources(),
    ).unwrap();
    assert_eq!(external, capacity);
    for (axis, quantity) in [("moisture", QuantitySpec::dimensional(Dims::NONE)),
        ("strain", QuantitySpec::dimensional(Dims::NONE)),
        ("strain_rate", QuantitySpec::dimensional(fs_qty::Frequency::DIMS)),
        ("stress", QuantitySpec::dimensional(Pressure::DIMS)),
        ("omega", QuantitySpec::dimensional(fs_qty::Frequency::DIMS))] {
        for curve in [false, true] {
            let unsupported = contextual_section(axis, quantity, curve);
            assert!(unsupported.axial_tensile_capacity(LIMIT).is_err(),
                "unowned axis {axis}, curve={curve}");
            let source = unsupported.layers()[0].material.property(LIMIT).unwrap();
            assert!(unsupported.axial_tensile_capacity_from_resolved(
                LIMIT, core::slice::from_ref(source),
            ).is_err(), "resolved source must not freeze unowned axis {axis}, curve={curve}");
        }
    }
}

#[test]
fn g0_resolved_drawn_wire_diameter_remains_an_explicit_fixed_source_coordinate() {
    let length = QuantitySpec::dimensional(fs_qty::Length::DIMS);
    for curve in [false, true] {
        let section = contextual_section("wire_diameter", length, curve);
        let source = section.layers()[0].material.property(LIMIT).unwrap();
        // Only the explicit external-source path owns the drawn-wire diameter.
        assert!(section.axial_tensile_capacity(LIMIT).is_err());
        let capacity = section.axial_tensile_capacity_from_resolved(
            LIMIT, core::slice::from_ref(source),
        ).unwrap();
        assert_eq!(capacity.layer_limit_sources(), core::slice::from_ref(source));
        assert_eq!(capacity.layer_limits_pa(), &[if curve {225e6} else {300e6}]);
        let untyped = contextual_section("wire_diameter", QuantitySpec::dimensional(Dims::NONE), curve);
        assert!(untyped.axial_tensile_capacity_from_resolved(
            LIMIT, core::slice::from_ref(untyped.layers()[0].material.property(LIMIT).unwrap()),
        ).is_err());
    }
}

#[test]
fn g0_nonfinite_loads_and_utilization_overflow_refuse() {
    let section = thermal(300.0);
    let capacity = section.section().axial_tensile_capacity(LIMIT).unwrap();
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(capacity.assess_axial_bending(bad, [0.0; 2]).is_err());
        assert!(capacity.assess_axial_bending(0.0, [0.0, bad]).is_err());
    }
    // Ordinary geometry with a very low positive sourced limit gives a finite
    // capacity. A finite force can still have an unrepresentable utilization.
    let (a, _) = layer(&card("weak", 200e9, 1e-5, 1e-290), 0.002, 300.0);
    let weak = ConcentricElasticSection::try_new(0.0, vec![a]).unwrap()
        .axial_tensile_capacity(LIMIT).unwrap();
    assert!(weak.capacity_n().is_finite() && weak.capacity_n() > 0.0);
    assert!(weak.utilization(1e100).is_err());
    assert!(weak.assess_axial_bending(1e100, [0.0; 2]).is_err());
    assert_eq!(weak.utilization(0.0).unwrap(), 0.0);
}
