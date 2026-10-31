//! Synthetic-data wiring regressions, not measured wound-string validation.
use fs_couple::acoustic_realize::string_mode_omega;
use fs_couple::material_assembly::AcousticMaterialBindings;
use fs_couple::material_assembly::concentric::{
    compile_concentric_material_assembly, with_concentric_section,
};
use fs_couple::string_specimen::StringPrestress;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Pressure};
use fs_scenario::{AcousticAssembly, AmbientGas, Listener, Pluck, PrestressedString};

fn material(name: &str, rho: f64, young: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [("density", rho, Density::DIMS), ("young_modulus", young, Pressure::DIMS)] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::RelativeHalfWidth { fraction: 0.01, confidence: 0.95 },
            interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
            provenance: Provenance { source: "synthetic elastic section fixture".into(), license: "test-only".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: name.into(), phase: "solid".into(),
        process: "synthetic elastic specimen".into(), revision: 0 }, claims, vec![]).unwrap();
    resolve_material_state_point(&card, &QueryPoint::new(), &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn section(sleeve: ResolvedMaterialStatePoint) -> ConcentricElasticSection {
    ConcentricElasticSection::try_new(0.0, vec![
        ConcentricMaterialLayer { outer_radius_m: 0.00025, material: material("core", 8000.0, 2e11) },
        ConcentricMaterialLayer { outer_radius_m: 0.0005, material: sleeve },
    ]).unwrap()
}
fn string() -> PrestressedString {
    PrestressedString { length_m: 0.65, tension_n: 10.0, lin_density_kg_m: 1.0,
        axial_stiffness_n: 0.0, width_m: 0.01, n_modes: 4, damping_ratio: 0.0,
        rayleigh: None, bending_stiffness_n_m2: 0.0, kelvin_voigt_bending: None,
        relaxation_bending: None, polarization_detune: 0.0, moving_end: false,
        tension_loss_factor: None, bending_loss_factor: None }
}
fn assembly() -> AcousticAssembly {
    AcousticAssembly { ambient: AmbientGas::sea_level(), string: Some(string()), duct: None,
        pluck: Some(Pluck { station_frac: 0.25, height_m: 1e-5 }), bow: None, blow: None,
        reed: None, soundboard: None, body_modes: vec![], plate: None, cavity: None,
        obstacles: vec![], contact_texture: None, listener: Listener { distance_m: 1.0 },
        // The tension-modulated string steps midpoint: its four-mode pole rate
        // (about 880 rad/s) needs dt <= 4.5e-5 s for the one-cent band.
        sample_rate_hz: 32000, duration_s: 0.008 }
}
fn near(a: f64, b: f64) {
    assert!((a-b).abs() <= 2e-12 * a.abs().max(b.abs()), "{a} != {b}");
}

#[test]
fn homogeneous_subdivision_and_hollow_limit_preserve_the_section_equations() {
    let source = material("uniform", 8000.0, 1e11);
    for inner in [0.0, 0.0001] {
        let whole = ConcentricElasticSection::try_new(inner, vec![
            ConcentricMaterialLayer { outer_radius_m: 0.0005, material: source.clone() },
        ]).unwrap();
        let split = ConcentricElasticSection::try_new(inner, vec![
            ConcentricMaterialLayer { outer_radius_m: 0.00025, material: source.clone() },
            ConcentricMaterialLayer { outer_radius_m: 0.0005, material: source.clone() },
        ]).unwrap();
        let area = core::f64::consts::PI * (0.0005_f64.powi(2) - inner.powi(2));
        let inertia = core::f64::consts::PI / 4.0 * (0.0005_f64.powi(4) - inner.powi(4));
        near(whole.area_m2(), area);
        near(whole.second_moment_m4(), inertia);
        near(whole.line_density_kg_m(), 8000.0 * area);
        near(whole.axial_stiffness_n(), 1e11 * area);
        near(whole.bending_stiffness_n_m2(), 1e11 * inertia);
        near(split.line_density_kg_m(), whole.line_density_kg_m());
        near(split.axial_stiffness_n(), whole.axial_stiffness_n());
        near(split.bending_stiffness_n_m2(), whole.bending_stiffness_n_m2());
        assert_eq!(split.layers()[0].material, source);
    }
}

#[test]
fn replacing_only_sleeve_data_changes_mass_stiffness_modes_and_actual_waveform() {
    let a = section(material("sleeve-a", 8000.0, 16e10));
    let b = section(material("sleeve-b", 16000.0, 8e10));
    let template = assembly();
    let snapshot = template.clone();
    let bindings = AcousticMaterialBindings::default();
    let ca = compile_concentric_material_assembly(&template, a, StringPrestress::FixedTension(10.0), &bindings).unwrap();
    let cb = compile_concentric_material_assembly(&template, b, StringPrestress::FixedTension(10.0), &bindings).unwrap();
    let sa = ca.compiled().assembly().string.as_ref().unwrap();
    let sb = cb.compiled().assembly().string.as_ref().unwrap();
    assert_eq!(template, snapshot);
    assert_eq!(sa.width_m, 0.001);
    assert_eq!(sa.width_m, sb.width_m);
    assert!(sb.lin_density_kg_m > sa.lin_density_kg_m);
    assert!(sb.axial_stiffness_n < sa.axial_stiffness_n);
    assert!(sb.bending_stiffness_n_m2 < sa.bending_stiffness_n_m2);
    for mode in 1..=4 { assert!(string_mode_omega(sb, mode) < string_mode_omega(sa, mode)); }
    let ra = ca.realize().unwrap();
    let rb = cb.realize().unwrap();
    assert_eq!(ra.gas, rb.gas);
    assert_eq!(ra.pressure_pa.len(), 256);
    assert!(ra.pressure_pa.iter().chain(&rb.pressure_pa).all(|v| v.is_finite()));
    assert!(ra.pressure_pa.iter().any(|v| *v != 0.0));
    assert_ne!(ra.pressure_pa, rb.pressure_pa);
}

#[test]
fn support_prescription_is_independent_of_material_and_retuning_is_not_implicit() {
    let section = section(material("sleeve", 10000.0, 8e10));
    let rest = 0.65 / 1.0002;
    let extended = with_concentric_section(string(), &section,
        StringPrestress::FixedExtension { stress_free_length_m: rest, linear_strain_limit: 0.001 }).unwrap();
    near(extended.tension_n, section.axial_stiffness_n() * ((0.65 - rest) / rest));
    let tuned = with_concentric_section(string(), &section, StringPrestress::TargetFundamentalHz(120.0)).unwrap();
    near(string_mode_omega(&tuned, 1), 2.0 * core::f64::consts::PI * 120.0);
    assert!(with_concentric_section(string(), &section,
        StringPrestress::FixedExtension { stress_free_length_m: 0.7, linear_strain_limit: 0.01 }).is_err());
    let mut lossy = string();
    lossy.damping_ratio = 0.01;
    assert!(with_concentric_section(lossy, &section, StringPrestress::FixedTension(10.0)).is_err());
    let mut moving = string();
    moving.moving_end = true;
    assert!(with_concentric_section(moving, &section, StringPrestress::TargetFundamentalHz(120.0)).is_err());
}

#[test]
fn malformed_layer_assignments_refuse_instead_of_sorting_or_filling_gaps() {
    let source = material("uniform", 8000.0, 1e11);
    assert!(ConcentricElasticSection::try_new(0.0, vec![]).is_err());
    for radius in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(ConcentricElasticSection::try_new(0.0, vec![
            ConcentricMaterialLayer { outer_radius_m: radius, material: source.clone() },
        ]).is_err());
    }
    assert!(ConcentricElasticSection::try_new(0.001, vec![
        ConcentricMaterialLayer { outer_radius_m: 0.0005, material: source.clone() },
    ]).is_err());
    assert!(ConcentricElasticSection::try_new(0.0, vec![
        ConcentricMaterialLayer { outer_radius_m: 0.0005, material: source.clone() },
        ConcentricMaterialLayer { outer_radius_m: 0.00025, material: source },
    ]).is_err());
}
