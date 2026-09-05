//! G0/G1/G3: material resolution → circular specimen → existing pressure solver.
//! Synthetic material data check the implementation, not measured wire fidelity.

use fs_couple::acoustic_realize::{realize_assembly, string_mode_omega};
use fs_couple::string_specimen::{
    StringGeometryConstraint, StringPrestress, with_uniform_circular_material_and_constraints,
    with_uniform_circular_material_and_prestress, with_uniform_circular_material_state,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{
    Density, Dims, Pressure, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_scenario::{
    AcousticAssembly, AmbientGas, Listener, Pluck, PrestressedString, RayleighParams,
};

fn state(name: &str, properties: &[(&str, Dims, f64)]) -> ResolvedMaterialStatePoint {
    let quantities: Vec<_> = properties
        .iter()
        .map(|&(key, dims, value)| (key, QuantitySpec::dimensional(dims), value))
        .collect();
    state_with_quantities(name, &quantities)
}

fn state_with_quantities(
    name: &str,
    properties: &[(&str, QuantitySpec, f64)],
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(key, quantity, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(key, quantity),
                value: PropertyValue::Scalar {
                    value,
                    dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained().with("T", 290.0, 300.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic uniform string".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_with_quantity(
                key,
                quantity,
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new().with("T", 293.15).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn elastic(name: &str, rho: f64, young: f64) -> ResolvedMaterialStatePoint {
    state(
        name,
        &[
            ("density", Density::DIMS, rho),
            ("young_modulus", Pressure::DIMS, young),
        ],
    )
}

fn template() -> PrestressedString {
    PrestressedString {
        length_m: 0.5,
        tension_n: 20.0,
        // These independent legacy inputs must be replaced by the binding.
        lin_density_kg_m: 1.0,
        axial_stiffness_n: 1.0,
        bending_stiffness_n_m2: 1.0,
        width_m: 1.0,
        n_modes: 1,
        damping_ratio: 0.0,
        rayleigh: Some(RayleighParams {
            alpha_per_s: 2.0,
            beta_s: 0.0,
        }),
        polarization_detune: 0.0,
        moving_end: false,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual / expected - 1.0).abs() < 2.0e-13,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_circular_mass_and_stiffness_match_diameter_formulas_and_beam_modes() {
    let material = elastic("specimen", 7800.0, 200.0e9);
    let resolved = with_uniform_circular_material_state(template(), 0.0005, &material).unwrap();
    let string = resolved.string();
    let diameter: f64 = 0.001;
    let area = core::f64::consts::PI * diameter.powi(2) / 4.0;
    let moment = core::f64::consts::PI * diameter.powi(4) / 64.0;
    close(resolved.area_m2(), area);
    close(resolved.second_moment_m4(), moment);
    close(resolved.mass_kg(), 7800.0 * area * 0.5);
    close(string.lin_density_kg_m, 7800.0 * area);
    close(string.axial_stiffness_n, 200.0e9 * area);
    close(string.bending_stiffness_n_m2, 200.0e9 * moment);
    assert_eq!(string.width_m.to_bits(), diameter.to_bits());
    assert_eq!(string.tension_n.to_bits(), template().tension_n.to_bits());
    assert_eq!(string.rayleigh, template().rayleigh);
    assert_eq!(resolved.material(), &material);
    for n in 1..=4 {
        let k = f64::from(n) * core::f64::consts::PI / 0.5;
        let omega = ((20.0 * k.powi(2) + 200.0e9 * moment * k.powi(4)) / (7800.0 * area)).sqrt();
        close(string_mode_omega(string, n as usize), omega);
    }
}

#[test]
fn g3_geometric_scaling_and_material_replacement_move_all_derived_inputs() {
    let material = elastic("first", 1000.0, 2.0e9);
    let first = with_uniform_circular_material_state(template(), 0.0005, &material).unwrap();
    let twice = with_uniform_circular_material_state(
        PrestressedString {
            length_m: 1.0,
            ..template()
        },
        0.001,
        &material,
    )
    .unwrap();
    close(twice.mass_kg() / first.mass_kg(), 8.0);
    close(
        twice.string().lin_density_kg_m / first.string().lin_density_kg_m,
        4.0,
    );
    close(
        twice.string().axial_stiffness_n / first.string().axial_stiffness_n,
        4.0,
    );
    close(
        twice.string().bending_stiffness_n_m2 / first.string().bending_stiffness_n_m2,
        16.0,
    );
    let changed = with_uniform_circular_material_state(
        template(),
        0.0005,
        &elastic("replacement", 4000.0, 6.0e9),
    )
    .unwrap();
    close(changed.mass_kg() / first.mass_kg(), 4.0);
    close(
        changed.string().axial_stiffness_n / first.string().axial_stiffness_n,
        3.0,
    );
    close(
        changed.string().bending_stiffness_n_m2 / first.string().bending_stiffness_n_m2,
        3.0,
    );
    assert_ne!(changed.specimen_identity(), first.specimen_identity());
    let same_values = with_uniform_circular_material_state(
        template(),
        0.0005,
        &elastic("different-lot", 1000.0, 2.0e9),
    )
    .unwrap();
    assert_eq!(same_values.string(), first.string());
    assert_ne!(same_values.specimen_identity(), first.specimen_identity());
}

#[test]
fn g1_fixed_mass_derives_geometry_stiffness_and_independent_prestress() {
    let mass = 0.001;
    let length = template().length_m;
    let young = 200.0e9;
    let constraint = StringGeometryConstraint::FixedMass(mass);
    for density in [1000.0, 4000.0] {
        let material = elastic("fixed-mass", density, young);
        for prestress in [
            StringPrestress::FixedTension(20.0),
            StringPrestress::FixedExtension {
                stress_free_length_m: 0.498,
                linear_strain_limit: 0.005,
            },
            StringPrestress::TargetFundamentalHz(160.0),
        ] {
            let specimen = with_uniform_circular_material_and_constraints(
                template(),
                constraint,
                &material,
                prestress,
            )
            .unwrap();
            let string = specimen.string();
            let area = mass / (density * length);
            let moment = area.powi(2) / (4.0 * core::f64::consts::PI);
            close(specimen.mass_kg(), mass);
            close(specimen.area_m2(), area);
            close(specimen.second_moment_m4(), moment);
            close(string.lin_density_kg_m, mass / length);
            close(string.axial_stiffness_n, young * area);
            close(string.bending_stiffness_n_m2, young * moment);
            close(string.width_m, (4.0 * area / core::f64::consts::PI).sqrt());
            let expected_tension = match prestress {
                StringPrestress::FixedTension(force) => force,
                StringPrestress::FixedExtension {
                    stress_free_length_m,
                    ..
                } => young * area * (length / stress_free_length_m - 1.0),
                StringPrestress::TargetFundamentalHz(hz) => {
                    4.0 * mass * length * hz.powi(2)
                        - young * moment * (core::f64::consts::PI / length).powi(2)
                }
            };
            close(string.tension_n, expected_tension);
            assert_eq!(specimen.geometry_constraint(), constraint);
            assert_eq!(specimen.prestress(), prestress);
            assert_eq!(specimen.material(), &material);
            assert_eq!(string.rayleigh, template().rayleigh);
            let same_geometry = with_uniform_circular_material_and_prestress(
                template(),
                string.width_m / 2.0,
                &material,
                prestress,
            )
            .unwrap();
            assert_eq!(same_geometry.string(), string);
            assert_eq!(
                same_geometry.specimen_identity(),
                specimen.specimen_identity()
            );
            assert_ne!(same_geometry.geometry_constraint(), constraint);
        }
    }
}

#[test]
fn g3_fixed_mass_material_swap_reaches_pressure_with_bending_change() {
    let bind = |density| {
        with_uniform_circular_material_and_constraints(
            template(),
            StringGeometryConstraint::FixedMass(0.001),
            &elastic("fixed-mass-pressure", density, 200.0e9),
            StringPrestress::FixedTension(20.0),
        )
        .unwrap()
    };
    let light = bind(1000.0);
    let heavy = bind(4000.0);
    close(heavy.mass_kg(), light.mass_kg());
    close(heavy.string().width_m / light.string().width_m, 0.5);
    close(
        heavy.string().axial_stiffness_n / light.string().axial_stiffness_n,
        0.25,
    );
    close(
        heavy.string().bending_stiffness_n_m2 / light.string().bending_stiffness_n_m2,
        1.0 / 16.0,
    );
    let mut pitches = Vec::new();
    for (specimen, density) in [(&light, 1000.0_f64), (&heavy, 4000.0_f64)] {
        // Eliminate A and I independently: f_n² = n² T/(4mL)
        // + n⁴ E pi m/(16 rho² L⁵). Fixed mass preserves only the tension term.
        let frequency = |n: f64| {
            (n.powi(2) * 20.0 / (4.0 * 0.001 * 0.5)
                + n.powi(4) * 200.0e9 * core::f64::consts::PI * 0.001
                    / (16.0 * density.powi(2) * 0.5_f64.powi(5)))
            .sqrt()
        };
        for n in 1..=4 {
            close(
                string_mode_omega(specimen.string(), n) / core::f64::consts::TAU,
                frequency(n as f64),
            );
        }
        let waveform = pressure(specimen.string());
        assert!(waveform.iter().all(|p| p.is_finite()));
        assert_eq!(waveform, pressure(specimen.string()));
        let measured = measured_hz(&waveform);
        assert!((measured / frequency(1.0) - 1.0).abs() < 0.01);
        pitches.push(measured);
    }
    assert!(
        pitches[0] > 1.04 * pitches[1],
        "bending must change pitch at fixed mass"
    );
    eprintln!(
        "G3 fixed mass, density x4: pressure pitch {:.6} -> {:.6} Hz",
        pitches[0], pitches[1]
    );
}

fn pressure(string: PrestressedString) -> Vec<f64> {
    realize_assembly(&AcousticAssembly {
        ambient: AmbientGas::sea_level(),
        string: Some(string),
        duct: None,
        pluck: Some(Pluck {
            station_frac: 0.5,
            height_m: 1.0e-6,
        }),
        bow: None,
        blow: None,
        reed: None,
        soundboard: None,
        body_modes: vec![],
        plate: None,
        cavity: None,
        obstacles: vec![],
        contact_texture: None,
        listener: Listener { distance_m: 1.0 },
        sample_rate_hz: 8000,
        duration_s: 0.12,
    })
    .unwrap()
    .pressure_pa
}

fn measured_hz(pressure: &[f64]) -> f64 {
    let crossings: Vec<_> = pressure
        .windows(2)
        .enumerate()
        .skip(160)
        .filter(|(_, p)| p[0] > 0.0 && p[1] <= 0.0)
        .map(|(i, p)| i as f64 + p[0] / (p[0] - p[1]))
        .collect();
    assert!(
        crossings.len() >= 5,
        "a live decaying pressure wave is required"
    );
    8000.0 * (crossings.len() - 1) as f64 / (crossings.last().unwrap() - crossings[0])
}

#[test]
fn g3_material_density_changes_pitch_in_the_existing_pressure_simulation() {
    let light =
        with_uniform_circular_material_state(template(), 0.0005, &elastic("light", 1000.0, 2.0e9))
            .unwrap();
    let heavy =
        with_uniform_circular_material_state(template(), 0.0005, &elastic("heavy", 4000.0, 2.0e9))
            .unwrap();
    let a = pressure(light.string());
    let b = pressure(heavy.string());
    assert!(a.iter().chain(&b).all(|p| p.is_finite()));
    assert_eq!(a, pressure(light.string()));
    let f_a = measured_hz(&a);
    let f_b = measured_hz(&b);
    eprintln!("G3 material density 1000 -> 4000 kg/m3: pressure pitch {f_a:.6} -> {f_b:.6} Hz");
    assert!(
        (f_a / f_b - 2.0).abs() < 0.015,
        "density ×4 must halve pitch: {f_a}, {f_b}"
    );
    let diameter: f64 = 0.001;
    let area = core::f64::consts::PI * diameter.powi(2) / 4.0;
    let moment = core::f64::consts::PI * diameter.powi(4) / 64.0;
    let k = core::f64::consts::PI / 0.5;
    let expected = ((20.0 * k.powi(2) + 2.0e9 * moment * k.powi(4)) / (1000.0 * area)).sqrt()
        / (2.0 * core::f64::consts::PI);
    assert!(
        (f_a / expected - 1.0).abs() < 0.01,
        "{f_a} vs beam reference {expected}"
    );
}

#[test]
fn g0_binding_refuses_missing_wrong_dimension_and_unrepresentable_inputs() {
    let material = elastic("valid", 1000.0, 2.0e9);
    for radius in [
        0.0,
        -0.001,
        f64::NAN,
        f64::INFINITY,
        f64::MAX,
        f64::MIN_POSITIVE,
    ] {
        assert!(with_uniform_circular_material_state(template(), radius, &material).is_err());
    }
    for bad in [
        state("missing", &[("density", Density::DIMS, 1000.0)]),
        state(
            "wrong-dimension",
            &[
                ("density", Dims::NONE, 1000.0),
                ("young_modulus", Pressure::DIMS, 2.0e9),
            ],
        ),
        elastic("negative", -1000.0, 2.0e9),
    ] {
        assert!(with_uniform_circular_material_state(template(), 0.0005, &bad).is_err());
    }
    let original = template();
    assert!(
        with_uniform_circular_material_state(
            PrestressedString {
                tension_n: f64::INFINITY,
                ..original
            },
            0.0005,
            &material,
        )
        .is_err()
    );
    assert_eq!(original, template());
}

#[test]
fn g0_binding_preserves_quantity_kinds_and_value_forms_at_the_solver_boundary() {
    let density = QuantitySpec::dimensional(Density::DIMS);
    let modulus = QuantitySpec::dimensional(Pressure::DIMS);
    for (property, quantity) in [
        (
            "density",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::MassConcentration,
                ValueForm::Static,
            )),
        ),
        (
            "young_modulus",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AcousticPressure,
                ValueForm::Rms,
            )),
        ),
    ] {
        let mut properties = [
            ("density", density, 1000.0),
            ("young_modulus", modulus, 2.0e9),
        ];
        let entry = properties.iter_mut().find(|p| p.0 == property).unwrap();
        assert_eq!(entry.1.dims(), quantity.dims());
        entry.1 = quantity;
        let material = state_with_quantities("explicit-semantic-alias", &properties);
        assert_eq!(
            material
                .property(property)
                .unwrap()
                .requirement()
                .quantity(),
            quantity
        );
        assert!(with_uniform_circular_material_state(template(), 0.0005, &material).is_err());
        assert!(
            with_uniform_circular_material_and_constraints(
                template(),
                StringGeometryConstraint::FixedMass(0.001),
                &material,
                StringPrestress::FixedTension(20.0),
            )
            .is_err()
        );
    }
}

#[test]
fn g0_fixed_mass_refuses_invalid_and_unrepresentable_geometry() {
    let material = elastic("mass-admission", 1000.0, 2.0e9);
    for mass in [
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::MAX,
        f64::MIN_POSITIVE,
    ] {
        assert!(
            with_uniform_circular_material_and_constraints(
                template(),
                StringGeometryConstraint::FixedMass(mass),
                &material,
                StringPrestress::FixedTension(20.0),
            )
            .is_err()
        );
    }
    for length in [
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::MAX,
        f64::MIN_POSITIVE,
    ] {
        assert!(
            with_uniform_circular_material_and_constraints(
                PrestressedString {
                    length_m: length,
                    ..template()
                },
                StringGeometryConstraint::FixedMass(0.001),
                &material,
                StringPrestress::FixedTension(20.0),
            )
            .is_err()
        );
    }
}

#[test]
fn g1_prestress_prescriptions_preserve_material_authority_and_solve_beam_tension() {
    let material = elastic("same-specimen", 1000.0, 2.0e9);
    let fixed = with_uniform_circular_material_state(template(), 0.0005, &material).unwrap();
    let extension = StringPrestress::FixedExtension {
        stress_free_length_m: 0.498,
        linear_strain_limit: 0.005,
    };
    let extended = with_uniform_circular_material_and_prestress(
        PrestressedString {
            tension_n: f64::NAN,
            ..template()
        },
        0.0005,
        &material,
        extension,
    )
    .unwrap();
    let area = core::f64::consts::PI * 0.001_f64.powi(2) / 4.0;
    close(extended.string().tension_n, 2.0e9 * area * 0.002 / 0.498);
    assert_eq!(extended.prestress(), extension);
    assert_eq!(extended.material(), &material);
    assert_eq!(extended.specimen_identity(), fixed.specimen_identity());
    assert_eq!(extended.mass_kg().to_bits(), fixed.mass_kg().to_bits());
    assert_eq!(extended.string().rayleigh, fixed.string().rayleigh);
    assert_ne!(
        extended.string().tension_n.to_bits(),
        fixed.string().tension_n.to_bits()
    );

    let target = StringPrestress::TargetFundamentalHz(160.0);
    let tuned = with_uniform_circular_material_and_prestress(template(), 0.0005, &material, target)
        .unwrap();
    let moment = core::f64::consts::PI * 0.001_f64.powi(4) / 64.0;
    let expected_tension = 4.0 * 1000.0 * area * 0.5_f64.powi(2) * 160.0_f64.powi(2)
        - core::f64::consts::PI.powi(2) * 2.0e9 * moment / 0.5_f64.powi(2);
    close(tuned.string().tension_n, expected_tension);
    close(
        string_mode_omega(tuned.string(), 1) / core::f64::consts::TAU,
        160.0,
    );
    assert_eq!(tuned.prestress(), target);
    assert_eq!(tuned.specimen_identity(), fixed.specimen_identity());
    assert_eq!(fixed.prestress(), StringPrestress::FixedTension(20.0));
}

#[test]
fn g3_constraint_choice_changes_material_swap_pressure_response() {
    let extension = StringPrestress::FixedExtension {
        stress_free_length_m: 0.498,
        linear_strain_limit: 0.005,
    };
    let compliant = with_uniform_circular_material_and_prestress(
        template(),
        0.0005,
        &elastic("compliant", 1000.0, 2.0e9),
        extension,
    )
    .unwrap();
    let stiff = with_uniform_circular_material_and_prestress(
        template(),
        0.0005,
        &elastic("stiff", 1000.0, 8.0e9),
        extension,
    )
    .unwrap();
    close(stiff.string().tension_n / compliant.string().tension_n, 4.0);
    let compliant_hz = measured_hz(&pressure(compliant.string()));
    let stiff_hz = measured_hz(&pressure(stiff.string()));
    assert!(
        (stiff_hz / compliant_hz - 2.0).abs() < 0.015,
        "fixed extension and E x4 must double small-amplitude pitch: {compliant_hz}, {stiff_hz}"
    );

    let target = StringPrestress::TargetFundamentalHz(160.0);
    let light = with_uniform_circular_material_and_prestress(
        template(),
        0.0005,
        &elastic("light-tuned", 1000.0, 2.0e9),
        target,
    )
    .unwrap();
    let heavy = with_uniform_circular_material_and_prestress(
        template(),
        0.0005,
        &elastic("heavy-tuned", 4000.0, 2.0e9),
        target,
    )
    .unwrap();
    assert!(heavy.string().tension_n > 4.0 * light.string().tension_n);
    let light_hz = measured_hz(&pressure(light.string()));
    let heavy_hz = measured_hz(&pressure(heavy.string()));
    for actual in [light_hz, heavy_hz] {
        assert!(
            (actual / 160.0 - 1.0).abs() < 0.01,
            "target pressure pitch: {actual}"
        );
    }
    eprintln!(
        "G3 fixed extension E x4: {compliant_hz:.6} -> {stiff_hz:.6} Hz; target pitch density x4: {light_hz:.6} -> {heavy_hz:.6} Hz"
    );
}

#[test]
fn g0_prestress_refuses_slack_overstrain_invalid_targets_and_moving_end_mismatch() {
    let material = elastic("valid", 1000.0, 2.0e9);
    let original = template();
    let bind = |prestress| {
        with_uniform_circular_material_and_prestress(original, 0.0005, &material, prestress)
    };
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(bind(StringPrestress::FixedTension(invalid)).is_err());
        assert!(bind(StringPrestress::TargetFundamentalHz(invalid)).is_err());
        assert!(
            bind(StringPrestress::FixedExtension {
                stress_free_length_m: invalid,
                linear_strain_limit: 0.005,
            })
            .is_err()
        );
        assert!(
            bind(StringPrestress::FixedExtension {
                stress_free_length_m: 0.498,
                linear_strain_limit: invalid,
            })
            .is_err()
        );
    }
    for stress_free_length_m in [0.5, 0.6, f64::MIN_POSITIVE] {
        assert!(
            bind(StringPrestress::FixedExtension {
                stress_free_length_m,
                linear_strain_limit: 0.005,
            })
            .is_err()
        );
    }
    assert!(
        bind(StringPrestress::FixedExtension {
            stress_free_length_m: 0.498,
            linear_strain_limit: 0.001,
        })
        .is_err()
    );
    for hz in [1.0e-10, f64::MAX] {
        assert!(bind(StringPrestress::TargetFundamentalHz(hz)).is_err());
    }
    let moving = PrestressedString {
        moving_end: true,
        ..original
    };
    assert!(
        with_uniform_circular_material_and_prestress(
            moving,
            0.0005,
            &material,
            StringPrestress::FixedTension(20.0),
        )
        .is_ok()
    );
    for prestress in [
        StringPrestress::TargetFundamentalHz(160.0),
        StringPrestress::FixedExtension {
            stress_free_length_m: 0.498,
            linear_strain_limit: 0.005,
        },
    ] {
        assert!(
            with_uniform_circular_material_and_prestress(moving, 0.0005, &material, prestress,)
                .is_err()
        );
    }
    assert_eq!(original, template());
}
