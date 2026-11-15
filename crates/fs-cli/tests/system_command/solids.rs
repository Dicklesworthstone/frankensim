//! `rigid_network` bodies declared as solids: a shape with its dimensions
//! and a material card, and no typed mass or inertia.
//!
//! examples/system/steel_flywheel.json is a steel disc of radius 0.1 m and
//! length 0.02 m. Its card is the AISI 1018 supplier-guide card of
//! data/matdb/seed-v1, compiled into examples/system/steel1018/. The disc is
//! spun up from rest by 1 N m about its axis for 1 s. The claims, declared
//! before the first run: the density the run reports is the card's printed
//! 0.284 lb/in^3; the mass and the principal inertia are those of the
//! cylinder at that density; the spin is `omega = tau t / I` and the torque's
//! work `tau^2 t^2 / (2 I)`, all of it stored. The same body as a tube, a
//! cuboid, a sphere and an ellipsoid gives each shape's closed form. Then the
//! refusals, among them a typed density or mass beside a solid.
//!
//! examples/system/glass_flywheel_hot.json is the same disc in soda-lime
//! silica glass, the NIST SRM 1826b card compiled into
//! examples/system/glass1826b/, declared at the card's 293.15 K and given a
//! `temperature_k` of 373.15 K. The claims, declared before the first run:
//! the density is the certificate's 2548.668 kg/m^3 at 293.15 K; the free
//! linear strain is the card's constant 9.6e-6 /K times the 80 K rise; the
//! mass is bit for bit that of the same disc without a temperature, its
//! principal inertia that disc's times (1 + strain)^2 and its volume times
//! (1 + strain)^3; the spin is `tau t / I` of the expanded disc.

use std::f64::consts::PI;

use super::*;

/// The card's printed density, 0.284 lb/in^3, in kg/m^3.
fn printed_density_kg_m3() -> f64 {
    0.284 * 0.453_592_37 / (0.0254 * 0.0254 * 0.0254)
}

/// The three numbers of the array `key` that follows `anchor` in a report.
fn triple_after(json: &str, anchor: &str, key: &str) -> [f64; 3] {
    let tail = &json[json
        .find(anchor)
        .unwrap_or_else(|| panic!("{anchor} not in report"))..];
    let field = format!("\"{key}\":[");
    let start = tail
        .find(&field)
        .unwrap_or_else(|| panic!("{key} not after {anchor}"))
        + field.len();
    let numbers: Vec<f64> = tail[start..]
        .split(']')
        .next()
        .unwrap()
        .split(',')
        .map(|number| number.parse().unwrap())
        .collect();
    numbers.try_into().unwrap()
}

fn relative(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

const CYLINDER: &str =
    "\"shape\": \"cylinder\",\n            \"radius_m\": 0.1,\n            \"length_m\": 0.02";

/// The example with its shape replaced, written beside its card pack.
fn reshaped(name: &str, shape: &str) -> PathBuf {
    let scene = std::fs::read_to_string(example("steel_flywheel.json")).unwrap();
    let changed = scene.replace(CYLINDER, shape);
    assert_ne!(changed, scene);
    let path = scratch(name, &changed);
    let pack = path.parent().unwrap().join("steel1018");
    std::fs::create_dir_all(&pack).unwrap();
    std::fs::copy(
        example("steel1018/steel-1018.fsmcdpk"),
        pack.join("steel-1018.fsmcdpk"),
    )
    .unwrap();
    path
}

#[test]
fn a_steel_flywheel_takes_its_mass_and_inertia_from_its_card_and_its_shape() {
    let (code, json, stderr) = run_binary(&["--json", "system", &example("steel_flywheel.json")]);
    assert_eq!(code, 0, "{stderr}");
    let rho = value(&json, "\"solids\"", "density_kg_m3");
    assert!(relative(rho, printed_density_kg_m3()) <= 1.0e-12, "{rho}");
    let (radius, length, torque, t) = (0.1_f64, 0.02_f64, 1.0_f64, 1.0_f64);
    let volume = PI * radius * radius * length;
    let mass = rho * volume;
    let transverse = mass * (3.0 * radius * radius + length * length) / 12.0;
    let axial = mass * radius * radius / 2.0;
    assert!(relative(value(&json, "\"solids\"", "volume_m3"), volume) <= 1.0e-15);
    assert!(relative(value(&json, "\"solids\"", "mass_kg"), mass) <= 1.0e-15);
    let inertia = triple_after(&json, "\"solids\"", "principal_inertia_kg_m2");
    for (got, want) in inertia.into_iter().zip([transverse, transverse, axial]) {
        assert!(relative(got, want) <= 1.0e-15, "{got} {want}");
    }
    for key in [
        "pack_hash",
        "card_hash",
        "density_receipt_hash",
        "selected_claim",
    ] {
        assert_eq!(reported_material_hash(&json, key).len(), 64);
    }
    // Spun up from rest about the axis: omega = tau t / I, the work
    // tau^2 t^2 / (2 I), all of it stored, and no motion off the axis.
    let flywheel = element("flywheel");
    let spin = value(&json, &flywheel, "body0.angular_velocity_z_rad_s");
    assert!(relative(spin, torque * t / axial) <= 1.0e-13, "{spin}");
    assert!(value(&json, &flywheel, "body0.angular_velocity_x_rad_s").abs() <= 0.0);
    assert!(value(&json, &flywheel, "body0.velocity_x_m_s").abs() <= 0.0);
    let work = torque * torque * t * t / (2.0 * axial);
    assert!(relative(value(&json, "flywheel.torque0", "energy_in_j"), work) <= 1.0e-13);
    assert!(relative(value(&json, &flywheel, "stored_change_j"), work) <= 1.0e-12);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-11);
    assert!(
        json.contains("a solid body's mass and principal inertia are"),
        "{json}"
    );
    assert!(json.contains("no material-card resolution besides solid bodies' densities"));
    // The text report names the body's mass.
    let (code, text, stderr) = run_binary(&["system", &example("steel_flywheel.json")]);
    assert_eq!(code, 0, "{stderr}");
    let line = text
        .lines()
        .find(|line| line.starts_with("flywheel.body0: "))
        .unwrap();
    let said: f64 = line["flywheel.body0: ".len()..]
        .split(' ')
        .next()
        .unwrap()
        .parse()
        .unwrap();
    assert!(relative(said, mass) <= 1.0e-12, "{line}");
}

/// The glass example changed by `edit`, written beside its card pack.
fn glass(name: &str, edit: impl Fn(&str) -> String) -> PathBuf {
    let scene = std::fs::read_to_string(example("glass_flywheel_hot.json")).unwrap();
    let changed = edit(&scene);
    assert_ne!(changed, scene, "{name}");
    let path = scratch(name, &changed);
    let pack = path.parent().unwrap().join("glass1826b");
    std::fs::create_dir_all(&pack).unwrap();
    std::fs::copy(
        example("glass1826b/glass-srm-1826b.fsmcdpk"),
        pack.join("glass-srm-1826b.fsmcdpk"),
    )
    .unwrap();
    path
}

const HOT: &str = "\"temperature_k\": 373.15,";

#[test]
fn a_glass_flywheel_given_a_temperature_is_its_cold_disc_freely_expanded_with_its_mass_kept() {
    let (code, json, stderr) =
        run_binary(&["--json", "system", &example("glass_flywheel_hot.json")]);
    assert_eq!(code, 0, "{stderr}");
    let cold = glass("glass_cold.json", |scene| scene.replace(HOT, ""));
    let (code, cold, stderr) = run_binary(&["--json", "system", cold.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let rho = value(&json, "\"solids\"", "density_kg_m3");
    assert!(relative(rho, 2548.668) <= 1.0e-15, "{rho}");
    for (key, kelvin) in [
        ("temperature_k", 373.15_f64),
        ("reference_temperature_k", 293.15),
    ] {
        assert_eq!(value(&json, "\"solids\"", key).to_bits(), kelvin.to_bits());
    }
    let strain = value(&json, "\"solids\"", "free_linear_strain");
    let want = 9.6e-6 * (373.15 - 293.15);
    assert!(relative(strain, want) <= 1.0e-12, "{strain} {want}");
    let stretch = 1.0 + strain;
    // Mass is kept bit for bit; the shape's moments scale as s^2, s^3.
    let mass = value(&json, "\"solids\"", "mass_kg");
    assert_eq!(
        mass.to_bits(),
        value(&cold, "\"solids\"", "mass_kg").to_bits()
    );
    let volume = PI * 0.1 * 0.1 * 0.02;
    assert!(relative(mass, rho * volume) <= 1.0e-15, "{mass}");
    assert_eq!(
        value(&json, "\"solids\"", "volume_m3").to_bits(),
        value(&cold, "\"solids\"", "volume_m3").to_bits()
    );
    let expanded = value(&json, "\"solids\"", "expanded_volume_m3");
    assert!(
        relative(expanded, volume * stretch.powi(3)) <= 1.0e-14,
        "{expanded}"
    );
    let hot = triple_after(&json, "\"solids\"", "principal_inertia_kg_m2");
    let reference = triple_after(&cold, "\"solids\"", "principal_inertia_kg_m2");
    for (got, was) in hot.into_iter().zip(reference) {
        assert!(
            relative(got, was * stretch * stretch) <= 1.0e-14,
            "{got} {was}"
        );
        assert!(got > was);
    }
    for key in ["expansion_hash", "expansion_claim", "density_receipt_hash"] {
        assert_eq!(reported_material_hash(&json, key).len(), 64);
    }
    assert!(!cold.contains("\"expansion_hash\""));
    // Spun by 1 N m for 1 s: the expanded disc's axial inertia sets the spin.
    let spin = value(
        &json,
        &element("flywheel"),
        "body0.angular_velocity_z_rad_s",
    );
    assert!(relative(spin, 1.0 / hot[2]) <= 1.0e-13, "{spin}");
    assert!(json.contains("freely and uniformly expanded"), "{json}");
    assert!(cold.contains("with no thermal expansion"), "{cold}");
}

#[test]
fn a_solid_temperature_refuses_a_state_without_a_temperature_and_a_coefficient_off_its_card() {
    let steel = std::fs::read_to_string(example("steel_flywheel.json")).unwrap();
    let changed = steel.replace("\"shape\"", "\"temperature_k\": 373.15, \"shape\"");
    assert_ne!(changed, steel);
    let path = scratch("steel_hot.json", &changed);
    let pack = path.parent().unwrap().join("steel1018");
    std::fs::create_dir_all(&pack).unwrap();
    std::fs::copy(
        example("steel1018/steel-1018.fsmcdpk"),
        pack.join("steel-1018.fsmcdpk"),
    )
    .unwrap();
    let mut cases = vec![(
        "steel_hot.json",
        path,
        "solid.temperature_k: the state has no absolute-temperature axis to expand along",
    )];
    for (name, from, to, said) in [
        (
            "glass_negative.json",
            HOT,
            "\"temperature_k\": -373.15,",
            "solid: temperature_k must be positive and finite",
        ),
        (
            "glass_range_known.json",
            "{\"name\": \"source-temperature-range-known\", \"value\": \"0\"",
            "{\"name\": \"source-temperature-range-known\", \"value\": \"1\"",
            "solid.temperature_k: ",
        ),
    ] {
        cases.push((name, glass(name, |scene| scene.replace(from, to)), said));
    }
    for (name, path, said) in cases {
        let (code, stdout, stderr) = run_binary(&["system", path.to_str().unwrap()]);
        assert!(stdout.is_empty(), "{name}: {stdout}");
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(
            stderr.contains("elements.flywheel.bodies[0]"),
            "{name}: {stderr}"
        );
        assert!(stderr.contains(said), "{name}: {stderr}");
    }
}

#[test]
fn each_shape_of_a_solid_gives_its_closed_form_mass_and_axial_inertia() {
    let rho = printed_density_kg_m3();
    let third = 4.0 / 3.0 * PI;
    for (name, shape, volume, axial_m2) in [
        (
            "tube.json",
            "\"shape\": \"cylinder\", \"radius_m\": 0.1, \"inner_radius_m\": 0.05, \"length_m\": 0.02",
            PI * (0.01 - 0.0025) * 0.02,
            (0.01 + 0.0025) / 2.0,
        ),
        (
            "cuboid.json",
            "\"shape\": \"cuboid\", \"sides_m\": [0.2, 0.1, 0.05]",
            0.2 * 0.1 * 0.05,
            (0.04 + 0.01) / 12.0,
        ),
        (
            "sphere.json",
            "\"shape\": \"sphere\", \"radius_m\": 0.05",
            third * 1.25e-4,
            0.001,
        ),
        (
            "ellipsoid.json",
            "\"shape\": \"ellipsoid\", \"semiaxes_m\": [0.1, 0.05, 0.02]",
            third * 1.0e-4,
            (0.01 + 0.0025) / 5.0,
        ),
    ] {
        let path = reshaped(name, shape);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{name}: {stderr}");
        let mass = rho * volume;
        let reported = value(&json, "\"solids\"", "mass_kg");
        assert!(
            relative(reported, mass) <= 1.0e-12,
            "{name}: {reported} {mass}"
        );
        let axial = mass * axial_m2;
        let inertia = triple_after(&json, "\"solids\"", "principal_inertia_kg_m2");
        assert!(
            relative(inertia[2], axial) <= 1.0e-12,
            "{name}: {inertia:?} {axial}"
        );
        let spin = value(
            &json,
            &element("flywheel"),
            "body0.angular_velocity_z_rad_s",
        );
        assert!(relative(spin, 1.0 / axial) <= 1.0e-12, "{name}: {spin}");
    }
}

#[test]
fn a_solid_refuses_typed_mass_or_density_other_shapes_dimensions_and_states_off_its_card() {
    let scene = std::fs::read_to_string(example("steel_flywheel.json")).unwrap();
    for (name, from, to, said) in [
        (
            "typed_mass.json",
            "\"position_m\"",
            "\"mass_kg\": 4.9, \"position_m\"",
            "mass_kg and principal_inertia_kg_m2 are not given with it",
        ),
        (
            "typed_inertia.json",
            "\"position_m\"",
            "\"principal_inertia_kg_m2\": [1, 1, 1], \"position_m\"",
            "mass_kg and principal_inertia_kg_m2 are not given with it",
        ),
        (
            "typed_density.json",
            "\"shape\"",
            "\"density_kg_m3\": 7800, \"shape\"",
            "unsupported field \"density_kg_m3\"",
        ),
        (
            "cone.json",
            "\"shape\": \"cylinder\"",
            "\"shape\": \"cone\"",
            "solid.shape: expected sphere, cuboid, ellipsoid or cylinder",
        ),
        (
            "sides.json",
            "\"radius_m\": 0.1,",
            "\"radius_m\": 0.1, \"sides_m\": [1, 1, 1],",
            "solid.sides_m: not a dimension of a cylinder",
        ),
        (
            "inner.json",
            "\"radius_m\": 0.1,",
            "\"radius_m\": 0.1, \"inner_radius_m\": 0.1,",
            "0 <= inner < outer",
        ),
        (
            "no_length.json",
            ",\n            \"length_m\": 0.02",
            "",
            "missing length_m",
        ),
        (
            "no_pack.json",
            "steel1018/steel-1018.fsmcdpk",
            "steel1018/absent.fsmcdpk",
            "solid.pack: ",
        ),
        (
            "off_card.json",
            "{\"name\": \"source-cold-finished\", \"value\": \"1\"",
            "{\"name\": \"source-cold-finished\", \"value\": \"0\"",
            "material inertia source",
        ),
        (
            "bad_kind.json",
            "\"kind\": \"dimensional\"}",
            "\"kind\": \"no-such-kind\"}",
            "solid.state: unsupported quantity kind",
        ),
    ] {
        let changed = scene.replace(from, to);
        assert_ne!(changed, scene, "{name}");
        let path = scratch(name, &changed);
        let pack = path.parent().unwrap().join("steel1018");
        std::fs::create_dir_all(&pack).unwrap();
        std::fs::copy(
            example("steel1018/steel-1018.fsmcdpk"),
            pack.join("steel-1018.fsmcdpk"),
        )
        .unwrap();
        let (code, stdout, stderr) = run_binary(&["system", path.to_str().unwrap()]);
        assert!(stdout.is_empty(), "{name}: {stdout}");
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(
            stderr.contains("elements.flywheel.bodies[0]"),
            "{name}: {stderr}"
        );
        assert!(stderr.contains(said), "{name}: {stderr}");
    }
}
