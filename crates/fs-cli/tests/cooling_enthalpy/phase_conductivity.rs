//! G1/G3 actual-binary phase-dependent conduction with explicit material laws.
//! Synthetic one-element algebra checks the accepted nonlinear endpoint.
use super::*;

const INITIAL: [f64; 4] = [500.0, 2000.0, 2100.0, 2200.0];
const DT: f64 = 0.02;
const SOURCE: &str = "Synthetic liquid-mass-fraction multiplier; no volume-mixture claim";

fn law(solid: f64, liquid: f64) -> J {
    J::Object(vec![
        ("law".into(), J::Str("linear-liquid-mass-fraction".into())),
        ("solid_multiplier".into(), number(solid)),
        ("liquid_multiplier".into(), number(liquid)),
        ("source".into(), J::Str(SOURCE.into())),
    ])
}

fn single_step(request: &mut J) {
    remove(request, "radiation");
    let schedule = member(request, "transient");
    put(schedule, "max_step_s", number(DT));
    put(
        schedule,
        "intervals",
        J::parse(r#"[{"duration_s":0.02,"power_scale":0}]"#).unwrap(),
    );
}

fn request() -> J {
    let mut request = J::parse(FIXTURE).unwrap();
    single_step(&mut request);
    let solid = member(&mut request, "solid");
    remove(solid, "conductivity_w_m_k");
    put(solid, "materials", J::parse(r#"[{"name":"base-conduction","conductivity_curve":{"temperature_k":[250,500],"conductivity_w_m_k":[8,18]},"source":"Synthetic base k(T), independently assigned"}]"#).unwrap());
    put(
        solid,
        "element_materials",
        J::parse(r#"["base-conduction"]"#).unwrap(),
    );
    let storage = member(member(&mut request, "transient"), "enthalpy");
    remove(storage, "initial_specific_enthalpy_j_kg");
    put(
        storage,
        "initial_specific_enthalpies_j_kg",
        J::Array(INITIAL.into_iter().map(number).collect()),
    );
    put(storage, "phase_conductivity", law(0.5, 4.0));
    request
}

fn make_named(request: &mut J) {
    let storage = member(member(request, "transient"), "enthalpy");
    let mut material = storage.clone();
    for key in ["initial_specific_enthalpies_j_kg", "newton"] {
        remove(&mut material, key);
    }
    put(&mut material, "name", J::Str("phase-storage".into()));
    for key in [
        "material_card_identity",
        "source",
        "reference_density_kg_m3",
        "knots",
        "phase",
        "phase_conductivity",
    ] {
        remove(storage, key);
    }
    put(storage, "materials", J::Array(vec![material]));
    put(
        storage,
        "element_materials",
        J::parse(r#"["phase-storage"]"#).unwrap(),
    );
}

fn check_echo(material: &J, solid: f64, liquid: f64) {
    let declaration = material.get("phase_conductivity").unwrap();
    assert_eq!(
        declaration.str_field("law"),
        Some("linear-liquid-mass-fraction")
    );
    near(n(declaration, "solid_multiplier"), solid, 0.0);
    near(n(declaration, "liquid_multiplier"), liquid, 0.0);
    assert_eq!(declaration.str_field("source"), Some(SOURCE));
}

fn check_endpoint(result: &J) {
    let h = values(result, "solid_specific_enthalpies_j_kg");
    let temperature = values(result, "solid_temperatures_k");
    let fractions = values(result, "solid_liquid_mass_fractions");
    assert!((0.0..1000.0).contains(&h[0]));
    near(temperature[0], 250.0 + h[0] / 10.0, 1e-9);
    near(fractions[0], 0.0, 0.0);
    for i in 1..4 {
        assert!((1000.0..3000.0).contains(&h[i]));
        near(temperature[i], 350.0, 1e-9);
        near(fractions[i], (h[i] - 1000.0) / 2000.0, 1e-10);
    }

    // Right unit tetrahedron: V=1/6, constant P1 gradients. Both k(T) and
    // s(f) use the accepted element means, not initial/frozen material data.
    let mean_t = temperature.iter().sum::<f64>() / 4.0;
    let scale = 0.5 + 3.5 * fractions.iter().sum::<f64>() / 4.0;
    let conductivity = scale * (8.0 + 0.04 * (mean_t - 250.0));
    let gradient = [
        temperature[1] - temperature[0],
        temperature[2] - temperature[0],
        temperature[3] - temperature[0],
    ];
    let gradients = [
        [-1.0, -1.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let faces = [
        ([1, 2, 3], 3.0_f64.sqrt() / 2.0),
        ([0, 2, 3], 0.5),
        ([0, 1, 3], 0.5),
        ([0, 1, 2], 0.5),
    ];
    let area = (3.0 + 3.0_f64.sqrt()) / 2.0;
    let wall_mean = faces
        .iter()
        .map(|(face, area)| area / 3.0 * face.iter().map(|&i| temperature[i]).sum::<f64>())
        .sum::<f64>()
        / area;
    let air_heat = (wall_mean - 300.0) * (1.0 - (-2.0 * area).exp());
    let air_reference = wall_mean - air_heat / (2.0 * area);
    for i in 0..4 {
        let conduction =
            conductivity / 6.0 * (0..3).map(|d| gradients[i][d] * gradient[d]).sum::<f64>();
        let mut convection = 0.0;
        for (face, area) in faces {
            if face.contains(&i) {
                for j in face {
                    convection += 2.0 * area / 12.0
                        * if i == j { 2.0 } else { 1.0 }
                        * (temperature[j] - air_reference);
                }
            }
        }
        near(
            RHO / 24.0 * (h[i] - INITIAL[i]) + DT * (conduction + convection),
            0.0,
            2e-7,
        );
    }
    let transient = result.get("transient").unwrap();
    near(n(transient, "stored_energy_change_j"), -DT * air_heat, 2e-7);
    near(n(transient, "energy_residual_j"), 0.0, 2e-7);
}

#[test]
fn uniform_and_named_phase_laws_reach_the_nonlinear_endpoint_and_preserve_declarations() {
    let mut input = request();
    let uniform = run(&input);
    check_echo(uniform.path(&["transient", "enthalpy"]).unwrap(), 0.5, 4.0);
    check_endpoint(&uniform);

    let mut base_only = input.clone();
    remove(
        member(member(&mut base_only, "transient"), "enthalpy"),
        "phase_conductivity",
    );
    let base_only = run(&base_only);
    assert!(
        base_only
            .path(&["transient", "enthalpy", "phase_conductivity"])
            .is_none()
    );
    assert!(
        values(&uniform, "solid_specific_enthalpies_j_kg")[0]
            > values(&base_only, "solid_specific_enthalpies_j_kg")[0] + 1.0,
        "the declared phase multiplier must change physical conduction"
    );

    make_named(&mut input);
    let named = run(&input);
    let materials = named
        .path(&["transient", "enthalpy", "materials"])
        .unwrap()
        .as_array()
        .unwrap();
    check_echo(&materials[0], 0.5, 4.0);
    check_endpoint(&named);
    for key in ["solid_specific_enthalpies_j_kg", "solid_temperatures_k"] {
        assert_eq!(values(&uniform, key), values(&named, key));
    }
}

#[test]
fn named_material_assignment_keeps_omitted_laws_at_identity() {
    let mut input = J::parse(CONTACT_FIXTURE).unwrap();
    single_step(&mut input);
    let storage = member(member(&mut input, "transient"), "enthalpy");
    put(
        storage,
        "initial_specific_enthalpies_j_kg",
        J::parse("[500,2000,2100,2200,500,3000,3100,3200]").unwrap(),
    );
    let J::Array(materials) = member(storage, "materials") else {
        panic!()
    };
    put(&mut materials[0], "phase_conductivity", law(0.5, 4.0));
    let implicit = run(&input);
    let declarations = implicit
        .path(&["transient", "enthalpy", "materials"])
        .unwrap()
        .as_array()
        .unwrap();
    check_echo(&declarations[0], 0.5, 4.0);
    assert!(declarations[1].get("phase_conductivity").is_none());

    let mut explicit = input.clone();
    let J::Array(materials) = member(
        member(member(&mut explicit, "transient"), "enthalpy"),
        "materials",
    ) else {
        panic!()
    };
    put(&mut materials[1], "phase_conductivity", law(1.0, 1.0));
    let explicit = run(&explicit);
    let J::Array(materials) = member(
        member(member(&mut input, "transient"), "enthalpy"),
        "materials",
    ) else {
        panic!()
    };
    materials.reverse();
    let reordered = run(&input);
    for key in ["solid_specific_enthalpies_j_kg", "solid_temperatures_k"] {
        assert_eq!(values(&implicit, key), values(&explicit, key));
        assert_eq!(values(&implicit, key), values(&reordered, key));
    }
    near(
        n(implicit.get("transient").unwrap(), "energy_residual_j"),
        0.0,
        2e-7,
    );
}

#[test]
fn malformed_phase_laws_and_conflicting_storage_policies_refuse_before_output() {
    let input = request();
    for (key, value) in [
        ("law", J::Str("linear-liquid-volume-fraction".into())),
        ("solid_multiplier", number(0.0)),
        ("liquid_multiplier", number(-1.0)),
        ("source", J::Str(String::new())),
    ] {
        let mut invalid = input.clone();
        let declaration = member(
            member(member(&mut invalid, "transient"), "enthalpy"),
            "phase_conductivity",
        );
        put(declaration, key, value);
        assert!(refuses(&invalid).contains("phase_conductivity"));
    }
    let mut missing_source = input.clone();
    remove(
        member(
            member(member(&mut missing_source, "transient"), "enthalpy"),
            "phase_conductivity",
        ),
        "source",
    );
    assert!(refuses(&missing_source).contains("source"));
    let mut nonfinite = input.clone();
    put(
        member(
            member(member(&mut nonfinite, "transient"), "enthalpy"),
            "phase_conductivity",
        ),
        "solid_multiplier",
        J::Number {
            value: f64::INFINITY,
            raw: "1e999".into(),
        },
    );
    refuses(&nonfinite);

    let mut named = input.clone();
    make_named(&mut named);
    put(
        member(member(&mut named, "transient"), "enthalpy"),
        "phase_conductivity",
        law(1.0, 2.0),
    );
    assert!(refuses(&named).contains("choose uniform enthalpy chart fields"));
    let mut duplicate_policy = input;
    put(
        member(&mut duplicate_policy, "transient"),
        "nonlinear",
        J::Object(Vec::new()),
    );
    assert!(refuses(&duplicate_policy).contains("transient.nonlinear"));
}
