//! Synthetic analytical fixtures; these constants are not empirical metal data.
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::state_point::{DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS,
    IntegratedIsotropicThermalExpansion, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    YOUNG_MODULUS_PROPERTY, integrate_isotropic_thermal_expansion, resolve_material_state_point};
use fs_qty::{Density, Pressure, Temperature};

fn card(name: &str, young: f64, alpha: f64, hot_modulus_ratio: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (key, dims, cold, hot) in [
        (DENSITY_PROPERTY, Density::DIMS, 8000.0, 8000.0),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, young, young * hot_modulus_ratio),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS, alpha, 2.0 * alpha),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, cold), (400.0, hot)], dims },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0).with("p", 90000.0, 200000.0),
            uncertainty: UncertaintyModel::RelativeHalfWidth { fraction: 0.01, confidence: 0.95 },
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance { source: "synthetic layered thermoelastic fixture".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: name.into(), phase: "solid".into(),
        process: "synthetic".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap().with("p", p).unwrap()
}
fn layer(card: &MaterialCard, radius: f64, t: f64) -> (ConcentricMaterialLayer, IntegratedIsotropicThermalExpansion) {
    let current = point(t, 100000.0);
    let material = resolve_material_state_point(card, &current, &[
        ScalarPropertyRequirement::try_new(DENSITY_PROPERTY, Density::DIMS, ScalarAdmissibility::StrictlyPositive).unwrap(),
        ScalarPropertyRequirement::try_new(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, ScalarAdmissibility::StrictlyPositive).unwrap(),
    ], MaterialPropertySelection::SingleClaimOnly).unwrap();
    let path = integrate_isotropic_thermal_expansion(card, &point(300.0, 100000.0),
        &current, MaterialPropertySelection::SingleClaimOnly).unwrap();
    (ConcentricMaterialLayer { outer_radius_m: radius, material }, path)
}
fn specimen(t: f64) -> ConcentricThermoelasticSection {
    let (a, pa) = layer(&card("core", 200e9, 1e-5, 1.0), 0.001, t);
    let (b, pb) = layer(&card("sleeve", 70e9, 2e-5, 1.0), 0.002, t);
    ConcentricThermoelasticSection::try_new(
        ConcentricElasticSection::try_new(0.0, vec![a, b]).unwrap(), vec![pa, pb],
    ).unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-11 * b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn restrained_layer_stresses_and_force_match_independent_integrals() {
    let s = specimen(350.0);
    // Integral of alpha0*(1+(T-300)/100) from 300 to 350 is 62.5*alpha0,
    // not the endpoint coefficient multiplied by 50 K.
    let response = s.response_at_axial_strain(0.0, 0.01).unwrap();
    close(response.layers()[0].stress_pa(), -200e9 * 0.000625);
    close(response.layers()[1].stress_pa(), -70e9 * 0.00125);
    let a = std::f64::consts::PI * 1e-6;
    close(response.force_n(), -a * (200e9 * 0.000625 + 3.0 * 70e9 * 0.00125));
    close(response.force_n(), -s.thermal_force_n());
    assert!(response.strain_energy_per_length_j_m() > 0.0);
}

#[test]
fn zero_resultant_preserves_self_equilibrated_layer_stress() {
    let s = specimen(350.0);
    let strain = (200.0 * 0.000625 + 210.0 * 0.00125) / 410.0;
    close(s.zero_resultant_axial_strain(), strain);
    let response = s.response_at_axial_strain(strain, 0.01).unwrap();
    assert!(response.force_n().abs() < 1e-9);
    assert!(response.layers()[0].stress_pa() > 0.0);
    assert!(response.layers()[1].stress_pa() < 0.0);
    assert!(response.strain_energy_per_length_j_m() > 0.0);
}

#[test]
fn reference_temperature_recovers_elastic_section_and_energy_derivative() {
    let s = specimen(300.0);
    let response = s.response_at_axial_strain(0.002, 0.01).unwrap();
    close(response.force_n(), s.section().axial_stiffness_n() * 0.002);
    close(response.strain_energy_per_length_j_m(), 0.5 * s.section().axial_stiffness_n() * 0.002 * 0.002);
    let s = specimen(350.0);
    let plus = s.response_at_axial_strain(0.002 + 1e-7, 0.01).unwrap();
    let minus = s.response_at_axial_strain(0.002 - 1e-7, 0.01).unwrap();
    let derivative = (plus.strain_energy_per_length_j_m() - minus.strain_energy_per_length_j_m()) / 2e-7;
    assert!((derivative / s.response_at_axial_strain(0.002, 0.01).unwrap().force_n() - 1.0).abs() < 1e-9);
}

#[test]
fn current_temperature_modulus_and_layer_temperature_both_propagate() {
    let a = card("temperature dependent core", 200e9, 1e-5, 0.5);
    let b = card("different-temperature sleeve", 70e9, 2e-5, 1.0);
    let (la, pa) = layer(&a, 0.001, 350.0);
    let (lb, pb) = layer(&b, 0.002, 325.0);
    let s = ConcentricThermoelasticSection::try_new(
        ConcentricElasticSection::try_new(0.0, vec![la, lb]).unwrap(), vec![pa, pb],
    ).unwrap();
    let r = s.response_at_axial_strain(0.0, 0.01).unwrap();
    close(r.layers()[0].stress_pa(), -150e9 * 0.000625);
    close(r.layers()[1].stress_pa(), -70e9 * (28.125 * 2e-5));
}

#[test]
fn identical_material_layer_splitting_preserves_response() {
    let source = card("homogeneous", 200e9, 1e-5, 1.0);
    let (a, pa) = layer(&source, 0.001, 350.0);
    let (b, pb) = layer(&source, 0.002, 350.0);
    let split = ConcentricThermoelasticSection::try_new(
        ConcentricElasticSection::try_new(0.0, vec![a, b.clone()]).unwrap(), vec![pa, pb.clone()],
    ).unwrap();
    let whole = ConcentricThermoelasticSection::try_new(
        ConcentricElasticSection::try_new(0.0, vec![b]).unwrap(), vec![pb],
    ).unwrap();
    let r = split.response_at_axial_strain(0.002, 0.01).unwrap();
    let q = whole.response_at_axial_strain(0.002, 0.01).unwrap();
    close(r.force_n(), q.force_n());
    close(r.strain_energy_per_length_j_m(), q.strain_energy_per_length_j_m());
    close(split.section().bending_stiffness_n_m2(), whole.section().bending_stiffness_n_m2());
}

#[test]
fn foreign_stale_reordered_and_inconsistent_reference_paths_refuse() {
    let a = card("a", 200e9, 1e-5, 1.0);
    let b = card("b", 70e9, 2e-5, 1.0);
    let (la, pa) = layer(&a, 0.001, 350.0);
    let (lb, pb) = layer(&b, 0.002, 350.0);
    let section = ConcentricElasticSection::try_new(0.0, vec![la, lb]).unwrap();
    for paths in [vec![], vec![pa.clone()], vec![pb.clone(), pa.clone()]] {
        assert!(ConcentricThermoelasticSection::try_new(section.clone(), paths).is_err());
    }
    for (reference, current) in [
        (point(300.0, 100000.0), point(325.0, 100000.0)),
        (point(300.0, 120000.0), point(350.0, 120000.0)),
        (point(310.0, 100000.0), point(350.0, 100000.0)),
    ] {
        let wrong = integrate_isotropic_thermal_expansion(&b, &reference, &current,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        assert!(ConcentricThermoelasticSection::try_new(section.clone(), vec![pa.clone(), wrong]).is_err());
    }
}

#[test]
fn strain_limits_apply_to_each_layer_and_failures_do_not_mutate_sources() {
    let s = specimen(350.0);
    let before = s.clone();
    for (strain, limit) in [(f64::NAN, 0.01), (0.0, 0.0), (0.0, 1.0),
        (0.0, f64::INFINITY), (0.0, 0.001), (0.002, 0.001), (-0.009, 0.01)] {
        assert!(s.response_at_axial_strain(strain, limit).is_err());
    }
    assert_eq!(s, before);
    assert_eq!(s.response_at_axial_strain(0.002, 0.01).unwrap(),
        s.response_at_axial_strain(0.002, 0.01).unwrap());
}
