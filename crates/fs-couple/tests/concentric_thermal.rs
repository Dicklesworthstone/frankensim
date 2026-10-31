//! Synthetic thermal/material substitutions through the real acoustic solver.
use fs_couple::acoustic_realize::string_mode_omega;
use fs_couple::material_assembly::AcousticMaterialBindings;
use fs_couple::material_assembly::concentric::with_concentric_section;
use fs_couple::material_assembly::concentric::thermal::{
    compile_concentric_thermal_material_assembly, with_concentric_thermal_extension,
};
use fs_couple::string_specimen::StringPrestress;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection,
    ScalarAdmissibility, ScalarPropertyRequirement, integrate_isotropic_thermal_expansion,
    resolve_material_state_point};
use fs_qty::{Density, Pressure};
use fs_scenario::{AcousticAssembly, AmbientGas, Listener, Pluck, PrestressedString};

fn section(t: f64, sleeve_alpha: f64) -> ConcentricThermoelasticSection {
    let mut layers = Vec::new();
    let mut paths = Vec::new();
    for (name, radius, rho, young, alpha) in [
        ("core", 0.00025, 8000.0, 2e11, 1e-5),
        ("sleeve", 0.0005, 10000.0, 8e10, sleeve_alpha),
    ] {
        let mut claims = ClaimSet::new();
        for (key, value, dims) in [("density", rho, Density::DIMS),
            ("young_modulus", young, Pressure::DIMS),
            ("linear_thermal_expansion_coefficient", alpha, INVERSE_TEMPERATURE_DIMS)]
        {
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
                uncertainty: UncertaintyModel::RelativeHalfWidth { fraction: 0.01, confidence: 0.95 },
                interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
                provenance: Provenance { source: "synthetic thermal section fixture".into(),
                    license: "CC0-1.0".into(), artifact: None },
            }).unwrap();
        }
        let card = MaterialCard::assemble(MaterialStateId { chemistry: name.into(), phase: "solid".into(),
            process: "synthetic".into(), revision: 0 }, claims, vec![]).unwrap();
        let point = QueryPoint::new().with("T", t).unwrap();
        let material = resolve_material_state_point(&card, &point, &[
            ScalarPropertyRequirement::try_new("density", Density::DIMS, ScalarAdmissibility::StrictlyPositive).unwrap(),
            ScalarPropertyRequirement::try_new("young_modulus", Pressure::DIMS, ScalarAdmissibility::StrictlyPositive).unwrap(),
        ], MaterialPropertySelection::SingleClaimOnly).unwrap();
        paths.push(integrate_isotropic_thermal_expansion(&card,
            &QueryPoint::new().with("T", 300.0).unwrap(), &point,
            MaterialPropertySelection::SingleClaimOnly).unwrap());
        layers.push(ConcentricMaterialLayer { outer_radius_m: radius, material });
    }
    ConcentricThermoelasticSection::try_new(ConcentricElasticSection::try_new(0.0, layers).unwrap(), paths).unwrap()
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
        // (about 2100 rad/s) needs dt <= 1.9e-5 s for the one-cent band.
        sample_rate_hz: 64000, duration_s: 0.008 }
}
fn near(a: f64, b: f64) {
    assert!((a-b).abs() <= 2e-12 * a.abs().max(b.abs()).max(1.0), "{a} != {b}");
}

#[test]
fn heating_and_expansion_data_change_tension_modes_and_actual_waveform() {
    let template = assembly();
    let before = template.clone();
    let bindings = AcousticMaterialBindings::default();
    let rest = 0.65 / 1.001;
    let compile = |t, alpha| compile_concentric_thermal_material_assembly(&template,
        section(t, alpha), rest, 0.01, &bindings).unwrap();
    let cold = compile(300.0, 2e-5);
    let warm = compile(325.0, 2e-5);
    let changed = compile(325.0, 3e-5);
    let a = cold.compiled().assembly().string.as_ref().unwrap();
    let b = warm.compiled().assembly().string.as_ref().unwrap();
    let c = changed.compiled().assembly().string.as_ref().unwrap();
    assert_eq!(template, before);
    assert_eq!(a.lin_density_kg_m, b.lin_density_kg_m);
    assert_eq!(a.bending_stiffness_n_m2, b.bending_stiffness_n_m2);
    assert_eq!(a.width_m, b.width_m);
    near(b.tension_n, warm.axial_response().force_n());
    near(a.tension_n - b.tension_n, warm.section().thermal_force_n());
    assert!(c.tension_n < b.tension_n && b.tension_n < a.tension_n);
    for mode in 1..=4 {
        assert!(string_mode_omega(c, mode) < string_mode_omega(b, mode));
        assert!(string_mode_omega(b, mode) < string_mode_omega(a, mode));
    }
    let ra = cold.realize().unwrap();
    let rb = warm.realize().unwrap();
    let rc = changed.realize().unwrap();
    assert_eq!(ra.gas, rb.gas);
    assert_eq!(rb.gas, rc.gas);
    assert_eq!(rb.pressure_pa.len(), 512);
    assert!(rb.pressure_pa.iter().chain(&rc.pressure_pa).all(|v| v.is_finite()));
    assert!(rb.pressure_pa.iter().any(|v| *v != 0.0));
    assert_ne!(ra.pressure_pa, rb.pressure_pa);
    assert_ne!(rb.pressure_pa, rc.pressure_pa);
}

#[test]
fn zero_temperature_change_recovers_existing_fixed_extension_binding() {
    let section = section(300.0, 2e-5);
    let rest = 0.65 / 1.001;
    let ordinary = with_concentric_section(string(), section.section(),
        StringPrestress::FixedExtension { stress_free_length_m: rest, linear_strain_limit: 0.01 }).unwrap();
    let (thermal, response) = with_concentric_thermal_extension(string(), &section, rest, 0.01).unwrap();
    near(ordinary.tension_n, thermal.tension_n);
    assert_eq!(ordinary.lin_density_kg_m, thermal.lin_density_kg_m);
    assert_eq!(ordinary.axial_stiffness_n, thermal.axial_stiffness_n);
    near(response.force_n(), thermal.tension_n);
    near(string_mode_omega(&ordinary, 1), string_mode_omega(&thermal, 1));
}

#[test]
fn thermal_binding_refuses_slack_moving_supports_inherited_loss_and_bad_lengths() {
    let rest = 0.65 / 1.001;
    assert!(with_concentric_thermal_extension(string(), &section(400.0, 2e-5), rest, 0.01).is_err());
    let s = section(325.0, 2e-5);
    for length in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(with_concentric_thermal_extension(string(), &s, length, 0.01).is_err());
    }
    let mut moving = string();
    moving.moving_end = true;
    assert!(with_concentric_thermal_extension(moving, &s, rest, 0.01).is_err());
    let mut lossy = string();
    lossy.damping_ratio = 0.01;
    assert!(with_concentric_thermal_extension(lossy, &s, rest, 0.01).is_err());
    assert!(with_concentric_thermal_extension(string(), &s, rest, 0.0001).is_err());
    let mut template = assembly();
    template.string = None;
    assert!(compile_concentric_thermal_material_assembly(&template, s, rest, 0.01,
        &AcousticMaterialBindings::default()).is_err());
}
