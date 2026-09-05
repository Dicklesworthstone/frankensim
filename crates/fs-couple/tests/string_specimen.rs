//! G0/G1/G3: material resolution → circular specimen → existing pressure solver.
//! Synthetic material data check the implementation, not measured wire fidelity.

use fs_couple::acoustic_realize::{realize_assembly, string_mode_omega};
use fs_couple::string_specimen::{
    BENDING_RELAXATION_TIME_PROPERTY, EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, RELAXING_BENDING_MODULUS_PROPERTY,
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
    Density, Dims, DynViscosity, Pressure, QuantitySpec, Time,
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
    state_with_domain(
        name,
        properties,
        ValidityDomain::unconstrained().with("T", 290.0, 300.0),
        QueryPoint::new().with("T", 293.15).unwrap(),
        None,
    )
}

fn state_with_domain(
    name: &str,
    properties: &[(&str, QuantitySpec, f64)],
    validity: ValidityDomain,
    point: QueryPoint,
    curve_property: Option<&str>,
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(key, quantity, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(key, quantity),
                value: if curve_property == Some(key) {
                    PropertyValue::Curve {
                        abscissa: "omega".into(),
                        abscissa_dims: Dims([0, 0, -1, 0, 0, 0]),
                        knots: vec![(1.0, value), (20_000.0, 2.0 * value)],
                        dims: quantity.dims(),
                    }
                } else {
                    PropertyValue::Scalar {
                        value,
                        dims: quantity.dims(),
                    }
                },
                validity: validity.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if curve_property == Some(key) {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
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
        &point,
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
        kelvin_voigt_bending: None,
        relaxation_bending: None,
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
    realize_assembly(&assembly(string)).unwrap().pressure_pa
}

fn assembly(string: PrestressedString) -> AcousticAssembly {
    AcousticAssembly {
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
    }
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

fn viscous_material(
    eta: f64,
    quantity: QuantitySpec,
    band: (f64, f64),
) -> ResolvedMaterialStatePoint {
    state_with_domain(
        "synthetic Kelvin-Voigt solid",
        &[
            ("density", QuantitySpec::dimensional(Density::DIMS), 1000.0),
            (
                "young_modulus",
                QuantitySpec::dimensional(Pressure::DIMS),
                2.0e9,
            ),
            (KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, quantity, eta),
        ],
        ValidityDomain::unconstrained()
            .with("T", 290.0, 300.0)
            .with("omega", band.0, band.1),
        QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", band.0)
            .unwrap(),
        None,
    )
}

fn loss_template() -> PrestressedString {
    PrestressedString {
        rayleigh: None,
        ..template()
    }
}

fn relaxing_material(
    delta: f64,
    tau: f64,
    band: (f64, f64),
    curve: Option<&str>,
) -> ResolvedMaterialStatePoint {
    state_with_domain(
        "synthetic standard-linear-solid",
        &[
            ("density", QuantitySpec::dimensional(Density::DIMS), 1000.0),
            (
                "young_modulus",
                QuantitySpec::dimensional(Pressure::DIMS),
                2.0e9,
            ),
            (
                EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
                QuantitySpec::dimensional(Pressure::DIMS),
                2.0e9,
            ),
            (
                RELAXING_BENDING_MODULUS_PROPERTY,
                QuantitySpec::dimensional(Pressure::DIMS),
                delta,
            ),
            (
                BENDING_RELAXATION_TIME_PROPERTY,
                QuantitySpec::dimensional(Time::DIMS),
                tau,
            ),
        ],
        ValidityDomain::unconstrained().with("omega", band.0, band.1),
        QueryPoint::new().with("omega", band.0).unwrap(),
        curve,
    )
}

#[test]
fn g1_sls_material_binding_recomputes_geometry_and_preserves_sources() {
    let material = relaxing_material(1.0e10, 0.005, (1.0, 20_000.0), None);
    let specimen = with_uniform_circular_material_state(loss_template(), 0.002, &material)
        .unwrap()
        .with_standard_linear_solid_bending_loss()
        .unwrap();
    let law = specimen.string().relaxation_bending.unwrap();
    close(
        law.relaxing_stiffness_n_m2,
        1.0e10 * core::f64::consts::PI * 0.004_f64.powi(4) / 64.0,
    );
    close(law.relaxation_time_s, 0.005);
    assert_eq!(law.material_state_identity, Some(material.identity()));
    assert_eq!(specimen.material(), &material);
    assert!(
        !material
            .property(RELAXING_BENDING_MODULUS_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .observation_backed
    );
    let replacement = relaxing_material(2.0e10, 0.01, (1.0, 30_000.0), None);
    let rebound =
        with_uniform_circular_material_state(specimen.string(), 0.004, &replacement).unwrap();
    let new_law = rebound.string().relaxation_bending.unwrap();
    close(
        new_law.relaxing_stiffness_n_m2 / law.relaxing_stiffness_n_m2,
        32.0,
    );
    close(new_law.relaxation_time_s, 0.01);
    assert_eq!(
        new_law.material_state_identity,
        Some(replacement.identity())
    );
    assert_eq!(new_law.omega_band_rad_s, (1.0, 30_000.0));
    assert!(
        with_uniform_circular_material_state(
            specimen.string(),
            0.004,
            &elastic("missing relaxation", 1000.0, 2.0e9)
        )
        .is_err()
    );
}

#[test]
fn g0_sls_material_refuses_missing_equilibrium_wrong_time_curves_and_duplicate_loss() {
    let bind = |material: &ResolvedMaterialStatePoint| {
        with_uniform_circular_material_state(loss_template(), 0.002, material)
            .unwrap()
            .with_standard_linear_solid_bending_loss()
    };
    assert!(bind(&elastic("missing", 1000.0, 2.0e9)).is_err());
    for (equilibrium, tau_dims, band) in [
        (3.0e9, Time::DIMS, true),
        (2.0e9, Pressure::DIMS, true),
        (2.0e9, Time::DIMS, false),
    ] {
        let material = state_with_domain(
            "invalid SLS",
            &[
                ("density", QuantitySpec::dimensional(Density::DIMS), 1000.0),
                (
                    "young_modulus",
                    QuantitySpec::dimensional(Pressure::DIMS),
                    2.0e9,
                ),
                (
                    EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
                    QuantitySpec::dimensional(Pressure::DIMS),
                    equilibrium,
                ),
                (
                    RELAXING_BENDING_MODULUS_PROPERTY,
                    QuantitySpec::dimensional(Pressure::DIMS),
                    1.0e10,
                ),
                (
                    BENDING_RELAXATION_TIME_PROPERTY,
                    QuantitySpec::dimensional(tau_dims),
                    0.005,
                ),
            ],
            if band {
                ValidityDomain::unconstrained().with("omega", 1.0, 20_000.0)
            } else {
                ValidityDomain::unconstrained()
            },
            QueryPoint::new().with("omega", 1.0).unwrap(),
            None,
        );
        assert!(bind(&material).is_err());
    }
    for key in [
        "density",
        "young_modulus",
        EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
        RELAXING_BENDING_MODULUS_PROPERTY,
        BENDING_RELAXATION_TIME_PROPERTY,
    ] {
        assert!(
            bind(&relaxing_material(
                1.0e10,
                0.005,
                (1.0, 20_000.0),
                Some(key)
            ))
            .is_err()
        );
    }
    for (delta, tau) in [(-1.0, 0.005), (1.0e10, 0.0)] {
        assert!(bind(&relaxing_material(delta, tau, (1.0, 20_000.0), None)).is_err());
    }
    let material = relaxing_material(1.0e10, 0.005, (1.0, 20_000.0), None);
    assert!(
        with_uniform_circular_material_state(template(), 0.002, &material)
            .unwrap()
            .with_standard_linear_solid_bending_loss()
            .is_err()
    );
    let specimen = bind(&material).unwrap();
    assert!(specimen.clone().with_kelvin_voigt_bending_loss().is_err());
    let mut duplicate = specimen.string();
    duplicate.damping_ratio = 0.01;
    assert!(realize_assembly(&assembly(duplicate)).is_err());
}

#[test]
fn g0_sls_realization_checks_instantaneous_band_nyquist_and_polarizations() {
    let material = relaxing_material(1.0e10, 0.005, (1.0, 20_000.0), None);
    let string = with_uniform_circular_material_state(loss_template(), 0.002, &material)
        .unwrap()
        .with_standard_linear_solid_bending_loss()
        .unwrap()
        .string();
    let w = string_mode_omega(string, 1);
    let mut narrow = string;
    narrow
        .relaxation_bending
        .as_mut()
        .unwrap()
        .omega_band_rad_s
        .1 = 1.01 * w;
    assert!(
        realize_assembly(&assembly(narrow)).is_err(),
        "instantaneous stiffness exceeds admitted band"
    );
    let mut valid = string;
    valid.relaxation_bending.as_mut().unwrap().omega_band_rad_s = (0.9 * w, 1.2 * w);
    assert!(realize_assembly(&assembly(valid)).is_ok());
    for modified in [
        PrestressedString {
            n_modes: 2,
            ..valid
        },
        PrestressedString {
            polarization_detune: 0.5,
            ..valid
        },
        PrestressedString {
            moving_end: true,
            ..valid
        },
        PrestressedString {
            moving_end: true,
            polarization_detune: 0.01,
            ..string
        },
    ] {
        assert!(realize_assembly(&assembly(modified)).is_err());
    }
    let mut too_stiff = string;
    too_stiff
        .relaxation_bending
        .as_mut()
        .unwrap()
        .relaxing_stiffness_n_m2 *= 1.0e6;
    too_stiff
        .relaxation_bending
        .as_mut()
        .unwrap()
        .omega_band_rad_s
        .1 = 1.0e9;
    assert!(
        realize_assembly(&assembly(too_stiff)).is_err(),
        "material band cannot license aliased modes"
    );
    let mut fast = string;
    fast.relaxation_bending.as_mut().unwrap().relaxation_time_s = 1.0e-6;
    assert!(
        realize_assembly(&assembly(fast))
            .unwrap_err()
            .to_string()
            .contains("dt/tau")
    );
    let mut coarse = assembly(string);
    coarse.sample_rate_hz = 1000; // below Nyquist, but not an accurate midpoint phase
    coarse
        .string
        .as_mut()
        .unwrap()
        .relaxation_bending
        .as_mut()
        .unwrap()
        .relaxation_time_s = 0.05;
    assert!(
        realize_assembly(&coarse)
            .unwrap_err()
            .to_string()
            .contains("reference phase")
    );
}

#[test]
fn g1_sls_pressure_pitch_and_decay_match_independent_characteristic_roots() {
    for (nonlinear, moving_end) in [(false, false), (true, false), (true, true)] {
        let mut pitches = Vec::new();
        for delta in [0.0, 1.0e10] {
            let material = relaxing_material(delta, 0.005, (1.0, 20_000.0), None);
            let mut string =
                with_uniform_circular_material_state(loss_template(), 0.002, &material)
                    .unwrap()
                    .with_standard_linear_solid_bending_loss()
                    .unwrap()
                    .string();
            string.moving_end = moving_end;
            if !nonlinear {
                string.axial_stiffness_n = 0.0;
            }
            let mut scene = assembly(string);
            scene.sample_rate_hz = 16_000;
            scene.duration_s = 0.4;
            let output = realize_assembly(&scene).unwrap();
            let pi = core::f64::consts::PI;
            let k = if moving_end { 0.5 } else { 1.0 } * pi / string.length_m;
            // Independent diameter-form reduction, not the runtime branch builder.
            let w2 = string.tension_n * k * k / string.lin_density_kg_m
                + 2.0e9 * string.width_m.powi(2) * k.powi(4) / (16.0 * 1000.0);
            let a = delta * string.width_m.powi(2) * k.powi(4) / (16.0 * 1000.0);
            let gas = output.gas;
            let c = (2.0 * pi * gas.dynamic_viscosity
                + 2.0
                    * pi
                    * string.width_m
                    * (gas.dynamic_viscosity * gas.density * w2.sqrt() / 2.0).sqrt())
                / string.lin_density_kg_m;
            // q'' + c q' + w0² q + a(q-v)=0, v'=(q-v)/tau.
            // The real root lies in [-1/tau,0]; divide the cubic to get the pair.
            let rate = 200.0;
            let (mut lo, mut hi) = (-rate, 0.0);
            for _ in 0..80 {
                let s = 0.5 * (lo + hi);
                if (s * s + c * s + w2) * (s + rate) + a * s > 0.0 {
                    hi = s;
                } else {
                    lo = s;
                }
            }
            let real = 0.5 * (lo + hi);
            let decay = 0.5 * (c + rate + real);
            let hz = (-w2 * rate / real - decay * decay).sqrt() / (2.0 * pi);
            let peaks: Vec<_> = output
                .pressure_pa
                .windows(3)
                .enumerate()
                .filter(|(i, p)| {
                    (960..5600).contains(i) && p[1] > 0.0 && p[1] > p[0] && p[1] >= p[2]
                })
                .map(|(i, p)| ((i + 1) as f64 / 16_000.0, p[1]))
                .collect();
            assert!(peaks.len() >= 4, "live oscillation required");
            let (t0, p0) = peaks[0];
            let (t1, p1) = *peaks.last().unwrap();
            let measured_hz = (peaks.len() - 1) as f64 / (t1 - t0);
            let measured_decay = (p0 / p1).ln() / (t1 - t0);
            assert!(
                (measured_hz / hz - 1.0).abs() < 0.005,
                "{measured_hz} vs {hz} Hz"
            );
            if delta > 0.0 {
                assert!(
                    (measured_decay / decay - 1.0).abs() < 0.03,
                    "nonlinear={nonlinear}, moving={moving_end}: {measured_decay} vs {decay} /s"
                );
            }
            pitches.push(measured_hz);
            eprintln!(
                "G1 SLS deltaE={delta}, nonlinear={nonlinear}, moving={moving_end}: {measured_hz} Hz / {measured_decay} per s; oracle {hz} / {decay}"
            );
        }
        assert!(
            pitches[1] > 1.005 * pitches[0],
            "relaxation must change storage, not just attenuate PCM"
        );
    }
}

#[test]
fn g1_material_bending_viscosity_follows_geometry_and_retains_sources() {
    let material = viscous_material(
        3.0e7,
        QuantitySpec::dimensional(DynViscosity::DIMS),
        (1.0, 20_000.0),
    );
    let specimen = with_uniform_circular_material_state(loss_template(), 0.0005, &material)
        .unwrap()
        .with_kelvin_voigt_bending_loss()
        .unwrap();
    let loss = specimen.string().kelvin_voigt_bending.unwrap();
    close(
        loss.viscous_stiffness_n_m2_s,
        3.0e7 * core::f64::consts::PI * 0.001_f64.powi(4) / 64.0,
    );
    assert_eq!(loss.material_state_identity, Some(material.identity()));
    assert_eq!(specimen.material(), &material);
    assert!(
        !specimen
            .material()
            .property(KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .observation_backed
    );
    let replacement = viscous_material(
        6.0e7,
        QuantitySpec::dimensional(DynViscosity::DIMS),
        (1.0, 30_000.0),
    );
    let rebound =
        with_uniform_circular_material_state(specimen.string(), 0.001, &replacement).unwrap();
    let new_loss = rebound.string().kelvin_voigt_bending.unwrap();
    close(
        new_loss.viscous_stiffness_n_m2_s / loss.viscous_stiffness_n_m2_s,
        32.0,
    );
    assert_eq!(
        new_loss.material_state_identity,
        Some(replacement.identity())
    );
    assert_eq!(new_loss.omega_band_rad_s.0.to_bits(), 1.0_f64.to_bits());
    assert_eq!(
        new_loss.omega_band_rad_s.1.to_bits(),
        30_000.0_f64.to_bits()
    );
    assert_ne!(rebound.specimen_identity(), specimen.specimen_identity());
    assert!(
        with_uniform_circular_material_state(
            specimen.string(),
            0.001,
            &elastic("missing viscosity", 1000.0, 2.0e9)
        )
        .is_err()
    );
}

#[test]
fn g3_material_viscosity_changes_pressure_decay_without_changing_elastic_modes() {
    let pi = core::f64::consts::PI;
    // The same material loss enters linear, nonlinear KC, and moving-end modes.
    for (nonlinear, moving_end) in [(false, false), (true, false), (true, true)] {
        let mut decays = Vec::new();
        let mut reference_frequency = None;
        for eta in [0.0, 3.0e7, 1.2e8] {
            let material = viscous_material(
                eta,
                QuantitySpec::dimensional(DynViscosity::DIMS),
                (1.0, 20_000.0),
            );
            let specimen = with_uniform_circular_material_state(loss_template(), 0.0005, &material)
                .unwrap()
                .with_kelvin_voigt_bending_loss()
                .unwrap();
            let mut string = specimen.string();
            string.moving_end = moving_end;
            if !nonlinear {
                string.axial_stiffness_n = 0.0;
            }
            let wave_number = if moving_end { 0.5 } else { 1.0 } * pi / string.length_m;
            let omega = ((string.tension_n * wave_number.powi(2)
                + string.bending_stiffness_n_m2 * wave_number.powi(4))
                / string.lin_density_kg_m)
                .sqrt();
            if let Some(first) = reference_frequency {
                assert_eq!(omega.to_bits(), first);
            }
            reference_frequency = Some(omega.to_bits());
            let mut scene = assembly(string);
            scene.sample_rate_hz = 16_000;
            scene.duration_s = 0.4;
            let output = realize_assembly(&scene).unwrap();
            assert!(output.pressure_pa.iter().all(|p| p.is_finite()));
            // Independent diameter-form projection: I/A = d²/16. The static
            // tensile energy must not be counted as dissipative bending energy.
            let bending_decay =
                eta * string.width_m.powi(2) * wave_number.powi(4) / (32.0 * 1000.0);
            let gas = output.gas;
            let resistance = 2.0 * pi * gas.dynamic_viscosity
                + 2.0
                    * pi
                    * string.width_m
                    * (gas.dynamic_viscosity * gas.density * omega / 2.0).sqrt();
            let expected = bending_decay + resistance / (2.0 * string.lin_density_kg_m);
            let peaks: Vec<_> = output
                .pressure_pa
                .windows(3)
                .enumerate()
                .filter(|(i, p)| {
                    (800..5600).contains(i) && p[1] > 0.0 && p[1] > p[0] && p[1] >= p[2]
                })
                .map(|(i, p)| ((i + 1) as f64 / 16_000.0, p[1]))
                .collect();
            assert!(peaks.len() > 20, "live pressure oscillation required");
            let (t0, p0) = peaks[0];
            let (t1, p1) = *peaks.last().unwrap();
            let measured = (p0 / p1).ln() / (t1 - t0);
            assert!(
                (measured / expected - 1.0).abs() < 0.02,
                "nonlinear={nonlinear}, moving={moving_end}, eta={eta}: {measured} vs {expected} /s"
            );
            decays.push(measured);
        }
        let ratio = (decays[2] - decays[0]) / (decays[1] - decays[0]);
        assert!(
            (ratio - 4.0).abs() < 0.12,
            "material-only decay ratio {ratio}"
        );
        eprintln!(
            "G3 Kelvin-Voigt viscosity 0, 3e7, 1.2e8 Pa s; nonlinear={nonlinear}, moving={moving_end}: {decays:?} /s"
        );
    }
}

#[test]
fn g0_material_bending_loss_refuses_missing_data_aliases_and_duplicate_losses() {
    let bind = |material: &ResolvedMaterialStatePoint| {
        with_uniform_circular_material_state(loss_template(), 0.0005, material)
            .unwrap()
            .with_kelvin_voigt_bending_loss()
    };
    assert!(bind(&elastic("missing", 1000.0, 2.0e9)).is_err());
    assert!(
        bind(&state(
            "no band",
            &[
                ("density", Density::DIMS, 1000.0),
                ("young_modulus", Pressure::DIMS, 2.0e9),
                (
                    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
                    DynViscosity::DIMS,
                    3.0e7
                )
            ]
        ))
        .is_err()
    );
    for (eta, quantity) in [
        (-1.0, QuantitySpec::dimensional(DynViscosity::DIMS)),
        (3.0e7, QuantitySpec::dimensional(Pressure::DIMS)),
    ] {
        assert!(bind(&viscous_material(eta, quantity, (1.0, 20_000.0))).is_err());
    }
    let material = viscous_material(
        3.0e7,
        QuantitySpec::dimensional(DynViscosity::DIMS),
        (1.0, 20_000.0),
    );
    assert!(
        with_uniform_circular_material_state(template(), 0.0005, &material)
            .unwrap()
            .with_kelvin_voigt_bending_loss()
            .is_err()
    );
    let specimen = bind(&material).unwrap();
    let mut duplicate = specimen.string();
    duplicate.damping_ratio = 0.01;
    assert!(realize_assembly(&assembly(duplicate)).is_err());
    duplicate = specimen.string();
    duplicate.rayleigh = template().rayleigh;
    assert!(realize_assembly(&assembly(duplicate)).is_err());
}

#[test]
fn g0_material_bending_loss_does_not_freeze_a_sampled_frequency_curve() {
    for key in [
        "density",
        "young_modulus",
        KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
    ] {
        let material = state_with_domain(
            "frequency-dependent synthetic solid",
            &[
                ("density", QuantitySpec::dimensional(Density::DIMS), 1000.0),
                (
                    "young_modulus",
                    QuantitySpec::dimensional(Pressure::DIMS),
                    2.0e9,
                ),
                (
                    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
                    QuantitySpec::dimensional(DynViscosity::DIMS),
                    3.0e7,
                ),
            ],
            ValidityDomain::unconstrained().with("omega", 1.0, 20_000.0),
            QueryPoint::new().with("omega", 1000.0).unwrap(),
            Some(key),
        );
        assert!(matches!(
            material.property(key).unwrap().answer().receipt.decision,
            fs_matdb::EvaluationDecision::LinearInside { .. }
        ));
        let sampled =
            with_uniform_circular_material_state(loss_template(), 0.0005, &material).unwrap();
        let error = sampled.with_kelvin_voigt_bending_loss().unwrap_err();
        assert!(
            error.to_string().contains("validity-wide scalar constants"),
            "{error}"
        );
    }
}

#[test]
fn g0_material_bending_loss_checks_every_retained_frequency_and_polarization() {
    let omega = string_mode_omega(
        with_uniform_circular_material_state(
            loss_template(),
            0.0005,
            &elastic("base", 1000.0, 2.0e9),
        )
        .unwrap()
        .string(),
        1,
    );
    let material = viscous_material(
        3.0e7,
        QuantitySpec::dimensional(DynViscosity::DIMS),
        (omega * 0.9, omega * 1.1),
    );
    let string = with_uniform_circular_material_state(loss_template(), 0.0005, &material)
        .unwrap()
        .with_kelvin_voigt_bending_loss()
        .unwrap()
        .string();
    assert!(realize_assembly(&assembly(string)).is_ok());
    for modified in [
        PrestressedString {
            n_modes: 2,
            ..string
        },
        PrestressedString {
            polarization_detune: 0.5,
            ..string
        },
        PrestressedString {
            moving_end: true,
            ..string
        },
    ] {
        assert!(realize_assembly(&assembly(modified)).is_err());
    }
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
