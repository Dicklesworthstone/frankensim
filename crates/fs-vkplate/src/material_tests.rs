//! G1/G3: canonical source data reaches actual kit operators and live pressure.
//! Manufactured coefficients establish causal wiring, not a measured instrument.

use super::*;
use crate::materials::{
    EDGE_LOSS_STATED_UNKNOWN_RECEIPT, FILM_LOSS_REFUSAL, FILM_LOSS_STATED_UNKNOWN_RECEIPT,
    HERTZ_HALF_SPACE_CAVEAT, KitMaterials, TIP_ANISOTROPY, axial_indentation_modulus,
};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::{
    Density, Dims, DynViscosity, Frequency, Pressure, QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};

fn pack(name: &str, rho: f64, young: f64, poisson: Option<f64>) -> NormalizedMaterialCardPack {
    pack_with_viscosity(name, rho, young, poisson, None)
}

fn pack_with_viscosity(
    name: &str,
    rho: f64,
    young: f64,
    poisson: Option<f64>,
    viscosity: Option<([f64; 6], (f64, f64), bool)>,
) -> NormalizedMaterialCardPack {
    manufactured_pack(name, rho, young, poisson, viscosity, &[])
}

/// `extra` adds dimensionless scalar claims with the same validity.
fn manufactured_pack(
    name: &str,
    rho: f64,
    young: f64,
    poisson: Option<f64>,
    viscosity: Option<([f64; 6], (f64, f64), bool)>,
    extra: &[(&str, f64)],
) -> NormalizedMaterialCardPack {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut claims = ClaimSet::new();
    // Retain the actual manufactured inputs, including the frequency-curve
    // prescription. These are numerical fixtures, not measured material data.
    let mut raw = format!(
        "{name}\ndensity={rho:?}\nyoung_modulus={young:?}\npoisson_ratio={poisson:?}\nviscosity={viscosity:?}\n"
    );
    if !extra.is_empty() {
        raw.push_str(&format!("extra={extra:?}\n"));
    }
    let artifact = hash_domain("frankensim.kit.material-test", raw.as_bytes());
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: name.into(),
            method: "manufactured SI coefficients supplied by this test".into(),
            artifact,
            caveats: raw,
            provenance: Provenance {
                source: format!("manufactured kit test {name}"),
                license: "CC0-1.0".into(),
                artifact: Some(artifact),
            },
        })
        .unwrap();
    for (key, value, dims) in [
        ("density", Some(rho), Density::DIMS),
        ("young_modulus", Some(young), Pressure::DIMS),
        ("poisson_ratio", poisson, Dims::NONE),
    ]
    .into_iter()
    .chain(
        extra
            .iter()
            .map(|&(key, value)| (key, Some(value), Dims::NONE)),
    ) {
        if let Some(value) = value {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(key, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained().with_quantity(
                        "T",
                        temperature,
                        273.15,
                        323.15,
                    ),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![observation],
                    provenance: Provenance {
                        source: format!("manufactured kit test {name}/{key}"),
                        license: "CC0-1.0".into(),
                        artifact: Some(artifact),
                    },
                })
                .unwrap();
        }
    }
    if let Some((values, (lo, hi), frequency_curve)) = viscosity {
        let frequency = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Frequency(FrequencyConvention::Cyclic),
            ValueForm::Static,
        ));
        for (i, (key, value)) in
            fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES
                .into_iter()
                .zip(values)
                .enumerate()
        {
            let curve = frequency_curve && i == 0;
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(key, DynViscosity::DIMS),
                    value: if curve {
                        PropertyValue::Curve {
                            abscissa: "frequency".into(),
                            abscissa_dims: Frequency::DIMS,
                            knots: vec![(lo, value), (hi, value * 2.0)],
                            dims: DynViscosity::DIMS,
                        }
                    } else {
                        PropertyValue::Scalar {
                            value,
                            dims: DynViscosity::DIMS,
                        }
                    },
                    validity: ValidityDomain::unconstrained()
                        .with_quantity("T", temperature, 273.15, 323.15)
                        .with_quantity("frequency", frequency, lo, hi),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: if curve {
                        InterpolationPolicy::LinearInside
                    } else {
                        InterpolationPolicy::ConstantWithinValidity
                    },
                    observations: vec![observation],
                    provenance: Provenance {
                        source: format!("manufactured {name}/{key}"),
                        license: "CC0-1.0".into(),
                        artifact: Some(artifact),
                    },
                })
                .unwrap();
        }
    }
    let claims = NormalizedPack::new(
        name,
        "manufactured-kit-g1-g3",
        artifact,
        "CC0-1.0 synthetic test data",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "manufactured numerical test".into(),
            revision: 0,
        },
        claims,
    )
    .unwrap()
}

/// These manufactured cards make no loss claim, so the binding declares the
/// head film loss a stated unknown; without that a drum refuses to build.
fn binding(
    primary: &NormalizedMaterialCardPack,
    secondary: &NormalizedMaterialCardPack,
    striker: &NormalizedMaterialCardPack,
) -> KitMaterials {
    KitMaterials::from_packs_at_temperature(
        &primary.to_bytes(),
        &secondary.to_bytes(),
        &striker.to_bytes(),
        "T",
        293.15,
        9e-6,
        0.005,
        0.05,
    )
    .and_then(KitMaterials::with_stated_unknown_film_loss)
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-11 * expected.abs().max(1e-20),
        "{actual} != {expected}"
    );
}

fn impact_k(ins: &Instrument) -> f64 {
    let Law::Power { k, .. } = ins.contacts.last().unwrap().law;
    k
}

#[test]
fn g1_g3_one_material_swap_changes_mass_radiation_and_pair_contact_at_fixed_geometry() {
    let film = pack("film", 1390.0, 4e9, Some(0.38));
    let dense = pack("dense film", 2780.0, 4e9, Some(0.38));
    let stiff = pack("stiff film", 1390.0, 8e9, Some(0.38));
    let stick = pack("elastic impactor", 800.0, 12e9, Some(0.3));
    let heavy_stick = pack("denser impactor", 1600.0, 12e9, Some(0.3));
    let base = binding(&film, &film, &stick);
    let mut spec = kit::snare_wires_off();
    spec.grid *= 6.0;
    let build = |m: &KitMaterials| {
        kit::drum_with_materials(&spec, &[], standard_gas(), m)
            .unwrap()
            .instrument
    };
    let a = build(&base);
    let b = build(&binding(&dense, &film, &stick));
    let c = build(&binding(&stiff, &film, &stick));
    let d = build(&binding(&film, &film, &heavy_stick));
    assert_eq!(a.bodies[0].xy, b.bodies[0].xy);
    assert_eq!(a.bodies[0].tension, b.bodies[0].tension);
    for (&mass, &changed) in a.bodies[0].mass.iter().zip(&b.bodies[0].mass) {
        near(changed, 2.0 * mass);
    }
    assert_eq!(a.bodies[1].mass, b.bodies[1].mass);
    near(b.rads[0].resistance, 0.5 * a.rads[0].resistance);
    near(impact_k(&b), impact_k(&a)); // Density is not Hertz elasticity.
    let e_star = 1.0 / ((1.0 - 0.3_f64.powi(2)) / 12e9 + (1.0 - 0.38_f64.powi(2)) / 8e9);
    near(impact_k(&c), 4.0 / 3.0 * e_star * 0.005_f64.sqrt());
    assert!(impact_k(&c) > impact_k(&a));
    near(a.bodies.last().unwrap().mass[0], 800.0 * 9e-6);
    near(
        d.bodies.last().unwrap().mass[0],
        2.0 * a.bodies.last().unwrap().mass[0],
    );
    near(impact_k(&d), impact_k(&a));
    assert_eq!(
        a.materials().unwrap().primary().resolved().card_identity(),
        film.card().content_hash()
    );
    assert_ne!(
        a.materials().unwrap().primary().resolved().identity(),
        b.materials().unwrap().primary().resolved().identity()
    );
}

#[test]
fn g1_plate_bending_frequency_follows_elasticity_and_density_from_the_same_card() {
    let film = pack("base plate", 1390.0, 4e9, Some(0.3));
    let stiff = pack("twice modulus", 1390.0, 8e9, Some(0.3));
    let dense = pack("twice density", 2780.0, 4e9, Some(0.3));
    let stick = pack("impactor", 800.0, 12e9, Some(0.3));
    let mut spec = kit::kick();
    spec.grid *= 8.0;
    spec.batter.tension = 0.0;
    spec.reso.tension = 0.0;
    spec.muffle = None;
    let build = |primary| {
        kit::drum_with_materials(&spec, &[], standard_gas(), &binding(primary, &film, &stick))
            .unwrap()
            .instrument
    };
    let mut a = build(&film);
    let mut b = build(&stiff);
    let mut c = build(&dense);
    let lambda = a.bodies[0].max_rate2();
    assert!(lambda > 0.0);
    near(b.bodies[0].max_rate2(), 2.0 * lambda);
    near(c.bodies[0].max_rate2(), 0.5 * lambda);
    // Image export retains the same source bundle used by downstream steppers.
    let image = b.image();
    assert_eq!(image.materials.as_deref(), b.materials());
}

#[test]
fn g3_imported_cards_reach_live_pressure_clones_and_overlapping_impacts() {
    let film = pack("live film", 1390.0, 4e9, Some(0.38));
    let stiff = pack("live alternative", 1390.0, 8e9, Some(0.38));
    let stick = pack("live striker", 800.0, 12e9, Some(0.3));
    let material = binding(&film, &film, &stick);
    assert_eq!(
        material.source_packs()[0].content_hash(),
        film.content_hash()
    );
    let mut a = live::Live::with_materials_record("kick", 6.0, "293.15 101325 0", material.clone())
        .unwrap();
    let mut b = live::Live::with_materials_record(
        "kick",
        6.0,
        "293.15 101325 0",
        binding(&stiff, &film, &stick),
    )
    .unwrap();
    for player in [&mut a, &mut b] {
        let before = player.ins.bodies.len();
        player.strike(0.04, 0.01, 1.0);
        player.strike(0.06, 0.02, 0.7);
        assert_eq!(player.ins.bodies.len(), before + 1);
        near(
            player.ins.bodies.last().unwrap().mass[0],
            material.striker_mass_kg(),
        );
    }
    let mut copy = a.clone();
    let (mut first, mut changed, mut cloned) = ([0.0f32; 1024], [0.0f32; 1024], [0.0f32; 1024]);
    a.render(&mut first, 1.0);
    b.render(&mut changed, 1.0);
    copy.render(&mut cloned, 1.0);
    assert!(first.iter().any(|x| x.abs() > 1e-9));
    assert!(first.iter().chain(&changed).all(|x| x.is_finite()));
    assert_eq!(first, cloned);
    assert_ne!(first, changed);
    assert_eq!(a.ins.materials(), copy.ins.materials());
    assert_eq!(a.ins.materials().unwrap().source_packs().len(), 3);
}

#[test]
fn g3_hihat_declared_striker_overrides_default_tip_and_survives_overlapping_impacts() {
    let bronze = pack("hi-hat bronze", 8800.0, 110e9, Some(0.34));
    let stick = pack("hi-hat striker", 800.0, 12e9, Some(0.3));
    // The existing binding declares a 5 mm striker, distinct from the 4 mm default.
    let material = binding(&bronze, &bronze, &stick);
    let coefficients = |law| {
        let Law::Power { k, alpha, mu } = law;
        (k, alpha, mu)
    };
    let expected = coefficients(material.impact_law());
    assert_ne!(expected, coefficients(kit::hertz_stick(0.004, 110e9, 0.34)));
    for name in ["hihat", "hihat-open", "hihat-closed"] {
        let mut player = live::Live::with_materials_record(
            name, 7.0, "293.15 101325 0", material.clone(),
        ).unwrap();
        assert_eq!(coefficients(player.ins.contacts.last().unwrap().law), expected, "{name}");
        let before = player.ins.bodies.len();
        player.strike(0.85 * player.radius, 0.0, 2.5);
        player.strike(0.8 * player.radius, 0.0, 1.5);
        assert_eq!(player.ins.bodies.len(), before + 1, "{name}");
        assert_eq!(coefficients(player.ins.contacts.last().unwrap().law), expected, "{name}");
        near(player.ins.bodies.last().unwrap().mass[0], material.striker_mass_kg());
        let copy = player.clone();
        assert_eq!(coefficients(copy.ins.contacts.last().unwrap().law), expected, "{name}");
        assert_eq!(copy.ins.materials(), player.ins.materials());
    }
}

#[test]
fn g0_material_pack_admission_refuses_incomplete_out_of_domain_and_malformed_sources() {
    let film = pack("complete", 1390.0, 4e9, Some(0.38)).to_bytes();
    let missing = pack("no poisson claim", 1390.0, 4e9, None).to_bytes();
    let admit = |primary: &[u8], axis, temperature| {
        KitMaterials::from_packs_at_temperature(
            primary,
            &[],
            &film,
            axis,
            temperature,
            9e-6,
            0.005,
            0.05,
        )
    };
    assert_eq!(
        admit(&film, "T", 293.15).unwrap().primary(),
        admit(&film, "T", 293.15).unwrap().secondary()
    );
    assert!(admit(&missing, "T", 293.15).is_err());
    assert!(admit(&film, "wrong_axis", 293.15).is_err());
    for temperature in [0.0, -1.0, f64::NAN, 600.0] {
        assert!(admit(&film, "T", temperature).is_err());
    }
    assert!(admit(b"not a material card", "T", 293.15).is_err());
    let mut trailing = film.clone();
    trailing.push(0);
    assert!(admit(&trailing, "T", 293.15).is_err());
    assert!(
        KitMaterials::from_packs_at_temperature(&film, &[], &film, "T", 293.15, 0.0, 0.005, 0.05)
            .is_err()
    );
}

fn viscous_binding(source: &NormalizedMaterialCardPack, angle: f64) -> KitMaterials {
    let striker = pack("viscous-test striker", 800.0, 12e9, Some(0.3)).to_bytes();
    KitMaterials::from_packs_with_bending_loss(
        &source.to_bytes(),
        &[],
        &striker,
        "T",
        293.15,
        9e-6,
        0.005,
        0.05,
        "frequency",
        1000.0,
        [angle, 0.0],
    )
    .unwrap()
}

fn viscous_plate(materials: &KitMaterials, edge: Edge) -> Body {
    let (law, angle) = materials.primary_bending_loss().unwrap();
    Body::plate_with_viscosity(
        "material-loss plate",
        &PlateSpec {
            radius: 0.06,
            hole: 0.0,
            grid: 0.012,
            edge,
            thickness: &|_| 0.001,
            rest_height: &|_| 0.0,
            young: materials.primary().young_modulus_pa(),
            poisson: materials.primary().poisson_ratio(),
            density: materials.primary().density_kg_m3(),
            tension: 0.0,
            bending: true,
            nonlinear: Nonlinear::None,
            t60: [(500.0, 0.2), (10000.0, 0.1)],
            support: None,
            muffle: None,
        },
        law,
        angle,
    )
    .unwrap()
}

#[test]
fn g1_bending_viscosity_has_the_sourced_tensor_power_and_physical_orientation() {
    let source = pack_with_viscosity(
        "anisotropic reduced viscosity",
        1390.0,
        4e9,
        Some(0.3),
        Some(([4e4, 7e4, 3e4, 1e4, 5e3, -4e3], (0.0, 1e6), false)),
    );
    for angle in [0.0, PI / 2.0] {
        let materials = viscous_binding(&source, angle);
        let body = viscous_plate(&materials, Edge::Free);
        assert_eq!(
            body.loss(),
            (0.0, 0.0),
            "source loss replaces the T60 prescription"
        );
        let v: Vec<f64> = body
            .xy
            .iter()
            .map(|[x, y]| 0.2 * x * x - 0.35 * y * y + 0.6 * x * y)
            .collect();
        let law = &body.viscosity.as_ref().unwrap().source;
        let local_curvature = if angle == 0.0 {
            [0.4, -0.7, 1.2]
        } else {
            [-0.7, 0.4, -1.2]
        };
        let (_, power_density) = law.stress_and_power(local_curvature).unwrap();
        let root = &body.viscosity.as_ref().unwrap().root;
        let points = root.ptr.windows(2).filter(|p| p[0] != p[1]).count() / BENDING_ROOT_ROWS_PER_NODE;
        let expected = power_density * 0.001_f64.powi(3) / 12.0 * points as f64 * body.h.powi(2);
        let mut rates = vec![0.0; root.n];
        root.mul_add(&v, 1.0, &mut rates);
        near(dot(&rates, &rates), expected);
        near(body.viscosity.as_ref().unwrap().c.quad(&v), expected);
    }
}

#[test]
fn g0_centered_material_dissipation_closes_the_actual_staggered_energy() {
    let source = pack_with_viscosity(
        "passive viscosity",
        1390.0,
        4e9,
        Some(0.3),
        Some(([8e4, 6e4, 3e4, 1e4, 4e3, -3e3], (0.0, 1e6), false)),
    );
    let mut body = viscous_plate(&viscous_binding(&source, 0.23), Edge::Clamped);
    for (i, p) in body.xy.iter().enumerate() {
        body.w[i] = 1e-5 * (1.0 - (p[0] * p[0] + p[1] * p[1]) / 0.06_f64.powi(2)).powi(2);
        body.wp[i] = body.w[i];
    }
    body.sync();
    let mut ins = Instrument::new(384000.0);
    ins.add(body);
    ins.validate_bending_loss_clock().unwrap();
    assert!(ins.stability().0 < 0.6);
    let initial = ins.bodies[0].discrete_energy(ins.k());
    let mut previous = initial;
    for _ in 0..512 {
        ins.step(&mut []);
        let body = &ins.bodies[0];
        let energy = body.discrete_energy(ins.k());
        assert!(body.viscosity.as_ref().unwrap().last_j >= 0.0);
        assert!(energy <= previous + 2e-12 * initial);
        assert!((energy + body.bending_loss_dissipated_j() - initial).abs() < 2e-10 * initial);
        previous = energy;
    }
    assert!(ins.bodies[0].bending_loss_dissipated_j() > initial * 0.01);
    let mut copy = ins.clone();
    for _ in 0..64 {
        ins.step(&mut []);
        copy.step(&mut []);
    }
    assert_eq!(ins.bodies[0].w, copy.bodies[0].w);
    assert_eq!(
        ins.bodies[0].bending_loss_dissipated_j(),
        copy.bodies[0].bending_loss_dissipated_j()
    );
}

#[test]
fn g1_kelvin_voigt_bending_converges_to_the_analytic_damped_mode() {
    let (young, nu, tau) = (4e9, 0.3, 2e-5);
    let eta = tau * young / (1.0 - nu * nu);
    let source = pack_with_viscosity(
        "proportional manufactured viscosity",
        1390.0,
        young,
        Some(nu),
        Some((
            [eta, eta, eta * (1.0 - nu) / 2.0, eta * nu, 0.0, 0.0],
            (0.0, 1e7),
            false,
        )),
    );
    let body = viscous_plate(&viscous_binding(&source, 0.0), Edge::Clamped);
    // A reference eigenvector of the actual elastic grid; the independent
    // time solution is q'' + tau*lambda*q' + lambda*q=0, not another stepper.
    let mut mode: Vec<f64> = (0..body.n)
        .map(|i| ((i * 17 % 13) as f64 - 6.0) / 6.0)
        .collect();
    for _ in 0..700 {
        let mut next = vec![0.0; body.n];
        body.stiff.mul_add(&mode, 1.0, &mut next);
        for (v, m) in next.iter_mut().zip(&body.mass) {
            *v /= m;
        }
        let norm = dot(&next, &next).sqrt();
        for v in &mut next {
            *v /= norm;
        }
        mode = next;
    }
    let modal_mass = body
        .mass
        .iter()
        .zip(&mode)
        .map(|(m, v)| m * v * v)
        .sum::<f64>();
    let lambda = body.stiff.quad(&mode) / modal_mass;
    near(
        body.viscosity.as_ref().unwrap().c.quad(&mode) / modal_mass,
        tau * lambda,
    );
    let gamma = tau * lambda / 2.0;
    let omega = (lambda - gamma * gamma).sqrt();
    assert!(omega.is_finite() && omega > 0.0);
    let exact =
        |t: f64| (-gamma * t).exp() * ((omega * t).cos() + gamma / omega * (omega * t).sin());
    let end = PI / (2.0 * lambda.sqrt());
    let mut errors = Vec::new();
    for steps in [40, 80, 160] {
        let k = end / steps as f64;
        let mut ins = Instrument::new(1.0 / k);
        let mut b = body.clone();
        for (i, v) in mode.iter().enumerate() {
            b.w[i] = *v;
            b.wp[i] = *v * exact(-k);
        }
        b.sync();
        ins.add(b);
        ins.validate_bending_loss_clock().unwrap();
        assert!(ins.stability().0 < 0.6);
        for _ in 0..steps {
            ins.step(&mut []);
        }
        let amplitude = ins.bodies[0]
            .w
            .iter()
            .zip(&mode)
            .zip(&body.mass)
            .map(|((w, v), m)| w * v * m)
            .sum::<f64>()
            / modal_mass;
        errors.push((amplitude - exact(end)).abs());
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (1.8..2.2).contains(&ratio),
            "backward viscous velocity is first-order: {errors:?}"
        );
    }
}

#[test]
fn g3_source_loss_changes_live_pressure_without_preset_decay_or_hidden_fallback() {
    let values = [1e6, 1e6, 4e5, 2e5, 0.0, 0.0];
    let source = pack_with_viscosity(
        "live loss",
        1390.0,
        4e9,
        Some(0.38),
        Some((values, (0.0, 1e6), false)),
    );
    let changed = pack_with_viscosity(
        "live doubled loss",
        1390.0,
        4e9,
        Some(0.38),
        Some((values.map(|v| 2.0 * v), (0.0, 1e6), false)),
    );
    let make = |source| {
        live::Live::with_materials_record(
            "kick",
            6.0,
            "293.15 101325 0",
            viscous_binding(source, 0.0),
        )
        .unwrap()
    };
    let mut a = make(&source);
    let mut b = make(&changed);
    assert_eq!(a.ins.bodies[0].mass, b.ins.bodies[0].mass);
    assert_eq!(a.ins.bodies[0].stiff.val, b.ins.bodies[0].stiff.val);
    assert_eq!(a.ins.bodies[0].s1, 0.0);
    assert_eq!(b.ins.bodies[0].s1, 0.0);
    // Air loading does not depend on the film-loss source: the sourced route
    // carries the head radiation too.
    assert!(a.ins.bodies[0].head_radiation().is_some());
    for player in [&mut a, &mut b] {
        player.strike(0.03, 0.01, 1.0);
    }
    let mut copy = a.clone();
    let (mut pa, mut pb, mut pc) = ([0.0; 2048], [0.0; 2048], [0.0; 2048]);
    a.render(&mut pa, 1.0);
    b.render(&mut pb, 1.0);
    copy.render(&mut pc, 1.0);
    assert!(pa.iter().chain(&pb).all(|x| x.is_finite()));
    assert!(pa.iter().any(|x| x.abs() > 1e-9));
    assert_ne!(pa, pb);
    assert_eq!(pa, pc);
    assert!(a.ins.bodies[0].bending_loss_dissipated_j() > 0.0);
    assert_eq!(
        a.ins.bodies[0].bending_loss_dissipated_j(),
        copy.ins.bodies[0].bending_loss_dissipated_j()
    );

    let zero = pack_with_viscosity(
        "explicit zero loss",
        1390.0,
        4e9,
        Some(0.3),
        Some(([0.0; 6], (0.0, 1e6), false)),
    );
    assert_eq!(
        viscous_plate(&viscous_binding(&zero, 0.0), Edge::Clamped).loss(),
        (0.0, 0.0)
    );
    for (band, curve) in [
        ((100.0, 1e6), false),
        ((0.0, 2000.0), false),
        ((0.0, 1e6), true),
    ] {
        let unsupported = pack_with_viscosity(
            "unsupported loss domain",
            1390.0,
            4e9,
            Some(0.3),
            Some((values, band, curve)),
        );
        let bytes = unsupported.to_bytes();
        let result = KitMaterials::from_packs_with_bending_loss(
            &bytes,
            &[],
            &bytes,
            "T",
            293.15,
            9e-6,
            0.005,
            0.05,
            "frequency",
            1000.0,
            [0.0; 2],
        )
        .and_then(|m| live::Live::with_materials_record("kick", 6.0, "293.15 101325 0", m));
        assert!(result.is_err());
    }
}

#[test]
fn g0_source_viscosity_participates_in_the_real_kit_clock_admission() {
    let source = |name, values| {
        pack_with_viscosity(
            name,
            1390.0,
            4e9,
            Some(0.3),
            Some((values, (0.0, 1e7), false)),
        )
    };
    let zero = source("clock zero loss", [0.0; 6]);
    let too_stiff = source(
        "clock excessive viscosity",
        [1e15, 1e15, 4e14, 2e14, 0.0, 0.0],
    );
    let mut spec = kit::kick();
    spec.grid *= 8.0;
    let admitted =
        kit::drum_with_materials(&spec, &[], standard_gas(), &viscous_binding(&zero, 0.0)).unwrap();
    assert_eq!(admitted.oversample, spec.oversample);
    // Identical inertia, elasticity, geometry and source band. The enormous
    // viscosity alone exceeds every supported clock; silently retaining the
    // old T60 or checking elastic CFL alone would incorrectly admit this kit.
    let result = kit::drum_with_materials(
        &spec,
        &[],
        standard_gas(),
        &viscous_binding(&too_stiff, 0.0),
    );
    let error = result.err().expect("excessive viscosity must refuse every supported clock");
    assert!(error.contains("material bending viscosity"), "{error}");
    assert!(error.contains("32x oversampling"), "{error}");
}

#[test]
fn g0_imported_head_films_without_a_loss_claim_refuse_by_name_instead_of_keeping_the_fitted_t60() {
    let film = pack("film card without loss", 1390.0, 4e9, Some(0.38)).to_bytes();
    let stick = pack("striker card", 800.0, 12e9, Some(0.3)).to_bytes();
    let undeclared =
        KitMaterials::from_packs_at_temperature(&film, &[], &stick, "T", 293.15, 9e-6, 0.005, 0.05)
            .unwrap();
    assert!(!undeclared.film_loss_is_stated_unknown());
    assert!(
        undeclared
            .info_json()
            .contains("\"filmLossStatedUnknown\":null")
    );
    let mut spec = kit::kick();
    spec.grid *= 6.0;
    assert_eq!(
        kit::drum_with_materials(&spec, &[], standard_gas(), &undeclared)
            .err()
            .as_deref(),
        Some(FILM_LOSS_REFUSAL)
    );
    for piece in ["kick", "snare", "snare-off", "tom-12", "tom-14"] {
        let live =
            live::Live::with_materials_record(piece, 8.0, "293.15 101325 0", undeclared.clone());
        assert_eq!(live.err().as_deref(), Some(FILM_LOSS_REFUSAL), "{piece}");
    }
    // Cymbal and hi-hat loss is a separate bead: the same binding still builds them.
    assert!(live::Live::with_materials_record("crash", 8.0, "293.15 101325 0", undeclared).is_ok());

    // A sourced bending law and the declaration exclude each other in either order.
    let source = pack_with_viscosity(
        "film card with loss",
        1390.0,
        4e9,
        Some(0.38),
        Some(([1e4, 1e4, 4e3, 2e3, 0.0, 0.0], (0.0, 1e7), false)),
    );
    let viscous = viscous_binding(&source, 0.0);
    assert!(viscous.clone().with_stated_unknown_film_loss().is_err());
    let (law, angle) = viscous.primary_bending_loss().unwrap();
    let declared = KitMaterials::from_packs_at_temperature(
        &source.to_bytes(),
        &[],
        &stick,
        "T",
        293.15,
        9e-6,
        0.005,
        0.05,
    )
    .and_then(KitMaterials::with_stated_unknown_film_loss)
    .unwrap();
    assert!(matches!(
        declared.with_bending_viscosity(law.clone(), law.clone(), [angle; 2]),
        Err(error) if error.contains("exclusive")
    ));
}

#[test]
fn g1_a_stated_unknown_film_loss_is_zero_surface_loss_with_its_receipt_and_the_declared_muffle() {
    let film = pack("declared film", 1390.0, 4e9, Some(0.38));
    let stick = pack("declared striker", 800.0, 12e9, Some(0.3));
    let declared = binding(&film, &film, &stick);
    assert!(declared.film_loss_is_stated_unknown());
    assert!(declared.info_json().contains(&format!(
        "\"filmLossStatedUnknown\":\"{FILM_LOSS_STATED_UNKNOWN_RECEIPT}\""
    )));
    // The kick declares a muffle on both heads; its dashpot is not film loss.
    let mut spec = kit::kick();
    spec.grid *= 6.0;
    let piece = kit::drum_with_materials(&spec, &[], standard_gas(), &declared).unwrap();
    assert!(
        piece
            .notes
            .iter()
            .any(|note| note == FILM_LOSS_STATED_UNKNOWN_RECEIPT)
    );
    // The preset's Mylar has this card's density and elasticity, so the two
    // builds differ only in the recording-fitted loss the declaration withholds.
    let preset = kit::drum_with_gas(&spec, &[], standard_gas()).unwrap();
    assert!(
        preset
            .notes
            .iter()
            .all(|note| note != FILM_LOSS_STATED_UNKNOWN_RECEIPT)
    );
    for head in 0..2 {
        let (body, fitted) = (
            &piece.instrument.bodies[head],
            &preset.instrument.bodies[head],
        );
        assert_eq!(body.mass, fitted.mass);
        assert_eq!(body.stiff.val, fitted.stiff.val);
        assert_eq!(body.loss(), (0.0, 0.0));
        assert_ne!(fitted.loss(), (0.0, 0.0));
        let floor = fitted.loss().0;
        for (s, f) in body.s0.iter().zip(&fitted.s0) {
            assert!(
                (s - (f - floor)).abs() <= 1e-12 * f,
                "{s} against {f} - {floor}"
            );
        }
        assert!(
            body.s0.iter().any(|&s| s > 0.0),
            "the muffle dashpot remains"
        );
        // The withheld fit is replaced by the baked exterior load, not by nothing.
        assert!(body.head_radiation().is_some());
        assert!(fitted.head_radiation().is_none());
    }
    for label in ["batter head radiation:", "resonant head radiation:"] {
        assert!(
            piece.notes.iter().any(|note| note.starts_with(label)),
            "{label}"
        );
    }
    let mut live =
        live::Live::with_materials_record("tom-12", 8.0, "293.15 101325 0", declared).unwrap();
    live.strike(0.03, 0.01, 1.0);
    let mut out = [0.0f32; 1024];
    live.render(&mut out, 1.0);
    assert!(out.iter().all(|x| x.is_finite()) && out.iter().any(|x| x.abs() > 1e-9));
}

/// The head-film card compiled from data/matdb/seed-v1/head-film-pet-mylar-dupont-matsuzawa1959.
const HEAD_FILM: &[u8] =
    include_bytes!("../packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk");

#[test]
fn g1_the_sourced_head_film_card_supplies_density_and_modulus_and_claims_no_loss() {
    let stick = pack("film-card striker", 800.0, 12e9, Some(0.3)).to_bytes();
    let at = |temperature| {
        KitMaterials::from_packs_at_temperature(
            HEAD_FILM,
            &[],
            &stick,
            "T",
            temperature,
            9e-6,
            0.005,
            0.05,
        )
    };
    let film = at(293.15).unwrap();
    near(film.primary().density_kg_m3(), 1390.0);
    near(film.primary().young_modulus_pa(), 4903.325e6);
    // The isotropic reduction and its directional spread are stated, then
    // the (isotropic) striker tip with the Hertz half-space caveat.
    let notes = film.isotropic_reduction_notes();
    assert_eq!(
        notes[..2],
        ["batter", "resonant"].map(|label| format!(
            "isotropic {label} head: E = mean(MD, TD) = 4903.3 MPa; MD/TD = 4805.3 / 5001.4 MPa (+/- 2.0 %)"
        ))
    );
    assert_eq!(
        notes[2..],
        [format!(
            "striker tip: isotropic, E / (1 - nu^2) = 13.187 GPa, E* on the batter film 3.995 GPa; {HERTZ_HALF_SPACE_CAVEAT}"
        )]
    );
    near(film.primary().poisson_ratio(), 0.38);
    // The datasheet states no test temperature; the card admits 20 to 23 C only.
    assert!(at(296.15).is_ok());
    for temperature in [293.0, 296.2, 288.15] {
        assert!(at(temperature).is_err(), "{temperature} K");
    }
    let loss = KitMaterials::from_packs_with_bending_loss(
        HEAD_FILM,
        &[],
        &stick,
        "T",
        293.15,
        9e-6,
        0.005,
        0.05,
        "frequency",
        1000.0,
        [0.0; 2],
    );
    assert!(loss.is_err(), "the card makes no bending-viscosity claim");
    let mut spec = kit::tom(12.0);
    spec.grid *= 6.0;
    assert_eq!(
        kit::drum_with_materials(&spec, &[], standard_gas(), &film)
            .err()
            .as_deref(),
        Some(FILM_LOSS_REFUSAL)
    );
    let declared = film.with_stated_unknown_film_loss().unwrap();
    let piece = kit::drum_with_materials(&spec, &[], standard_gas(), &declared).unwrap();
    assert!(
        piece
            .notes
            .iter()
            .any(|note| note == FILM_LOSS_STATED_UNKNOWN_RECEIPT)
    );
    assert!(
        piece
            .notes
            .iter()
            .any(|note| note.starts_with("isotropic batter head: E = mean(MD, TD)"))
    );
    assert_eq!(piece.instrument.bodies[0].loss(), (0.0, 0.0));
}

/// The stick card compiled from data/matdb/seed-v1/stick-hickory-shagbark-fpl-gtr282.
const HICKORY: &[u8] = include_bytes!("../packs/stick-hickory-shagbark-fpl-gtr282.fsmatpk");

#[test]
fn g1_the_hickory_stick_card_sets_the_striker_mass_and_its_hertz_pair_with_the_film() {
    let at = |primary: &[u8], temperature| {
        KitMaterials::from_packs_at_temperature(
            primary,
            &[],
            HICKORY,
            "T",
            temperature,
            9e-6,
            0.005,
            0.05,
        )
    };
    let sticks = at(HEAD_FILM, 293.15).unwrap();
    // Table 5-3a, shagbark at 12%: SG 0.72 -> 806.4 kg/m^3; E_L = 1.1 x 14,900 MPa.
    near(sticks.striker().density_kg_m3(), 806.4);
    near(sticks.striker().young_modulus_pa(), 16_390e6);
    near(sticks.striker().poisson_ratio(), 0.371);
    near(sticks.striker_mass_kg(), 806.4 * 9e-6);
    // Ruling 6441: the card's composition-proxy elastic ratios (hardwood
    // means of Tables 5-1 and 5-2) make the tip transversely isotropic about
    // the grain, and its Hertz modulus is the axial indentation modulus.
    let e_l = 16_390e6;
    let tip = closed_form_axial_modulus(e_l, 0.05693 * e_l, 0.06657 * e_l, 0.4668, 0.518);
    near(sticks.striker_tip_modulus_pa(), tip);
    let e_star = 1.0 / (1.0 / tip + (1.0 - 0.38_f64.powi(2)) / 4_903.325e6);
    let Law::Power { k, .. } = sticks.impact_law();
    near(k, 4.0 / 3.0 * e_star * 0.005_f64.sqrt());
    // The receipt prints both ends of the bracket and the half-space caveat.
    assert_eq!(
        sticks.isotropic_reduction_notes().last().map(String::as_str),
        Some(
            format!(
                "striker tip: axial indentation modulus (Delafargue-Ulm 2004, transversely isotropic about the grain from E_L and the card's E_T/E_L, G_LT/E_L, nu_LT, nu_P) 7.245 GPa, E* on the batter film 3.200 GPa; isotropic-E_L end 19.006 GPa, E* 4.403 GPa; contact time 1.136 x the isotropic end; {HERTZ_HALF_SPACE_CAVEAT}"
            )
            .as_str()
        )
    );
    // The stick card alone admits 15 to 25 C (Table 5-15 bounds the MOE change).
    let wide = pack("wide-band film", 1390.0, 4e9, Some(0.38)).to_bytes();
    for temperature in [288.15, 298.15] {
        assert!(at(&wide, temperature).is_ok(), "{temperature} K");
    }
    for temperature in [287.9, 298.4] {
        assert!(at(&wide, temperature).is_err(), "{temperature} K");
    }
}

/// Delafargue and Ulm's M3 from the closed-form stiffnesses of a solid
/// transversely isotropic about axis 3: an oracle independent of the
/// compliance inversion in `axial_indentation_modulus`.
fn closed_form_axial_modulus(e_l: f64, e_p: f64, g_lp: f64, nu_lp: f64, nu_p: f64) -> f64 {
    let d = 1.0 - nu_p - 2.0 * nu_lp * nu_lp * e_p / e_l;
    let c11 = e_p * (1.0 - nu_lp * nu_lp * e_p / e_l) / ((1.0 + nu_p) * d);
    let (c13, c33) = (e_p * nu_lp / d, e_l * (1.0 - nu_p) / d);
    2.0 * ((c11 * c33 - c13 * c13) / c11 / (1.0 / g_lp + 2.0 / ((c11 * c33).sqrt() + c13))).sqrt()
}

#[test]
fn g1_the_axial_indentation_modulus_reduces_to_e_over_one_minus_nu_squared_when_isotropic() {
    let e = 12e9;
    for nu in [0.0, 0.25, 0.371, 0.45] {
        let m = axial_indentation_modulus(e, e, e / (2.0 * (1.0 + nu)), nu, nu).unwrap();
        near(m, e / (1.0 - nu * nu));
    }
    let e_l = 16_390e6;
    let hardwood = (0.05693 * e_l, 0.06657 * e_l, 0.4668, 0.518);
    near(
        axial_indentation_modulus(e_l, hardwood.0, hardwood.1, hardwood.2, hardwood.3).unwrap(),
        closed_form_axial_modulus(e_l, hardwood.0, hardwood.1, hardwood.2, hardwood.3),
    );
    // Constants that are not positive definite refuse.
    for (nu_lp, nu_p) in [(hardwood.2, 1.0), (5.0, hardwood.3)] {
        assert!(axial_indentation_modulus(e_l, hardwood.0, hardwood.1, nu_lp, nu_p).is_err());
    }
}

/// A striker card carries all of TIP_ANISOTROPY or none. Its isotropic values
/// (E_T = E_L, G = E / 2(1 + nu), nu_LT = nu_P = nu) reproduce the isotropic
/// Hertz pair; a partial set refuses instead of leaving the tip isotropic.
#[test]
fn g1_a_striker_tip_takes_all_four_anisotropy_claims_or_none() {
    let film = pack("tip film", 1390.0, 4.9e9, Some(0.38));
    let nu = 0.3;
    let striker = |claims: &[(&str, f64)]| {
        manufactured_pack("tip striker", 800.0, 12e9, Some(nu), None, claims)
    };
    let isotropic = [1.0, 1.0 / (2.0 * (1.0 + nu)), nu, nu];
    let claims: Vec<_> = TIP_ANISOTROPY.into_iter().zip(isotropic).collect();
    let plain = binding(&film, &film, &striker(&[]));
    let limit = binding(&film, &film, &striker(&claims));
    let Law::Power { k: k_plain, .. } = plain.impact_law();
    let Law::Power { k: k_limit, .. } = limit.impact_law();
    near(k_limit, k_plain);
    near(limit.striker_tip_modulus_pa(), 12e9 / (1.0 - nu * nu));
    assert!(plain.isotropic_reduction_notes()[0].starts_with("striker tip: isotropic"));
    assert!(
        limit.isotropic_reduction_notes()[0].starts_with("striker tip: axial indentation modulus")
    );
    let bind = |claims: &[(&str, f64)]| {
        KitMaterials::from_packs_at_temperature(
            &film.to_bytes(),
            &[],
            &striker(claims).to_bytes(),
            "T",
            293.15,
            9e-6,
            0.005,
            0.05,
        )
    };
    for missing in 0..4 {
        let mut partial = claims.clone();
        partial.remove(missing);
        let refused = bind(&partial).err();
        assert!(
            refused
                .as_deref()
                .is_some_and(|e| e.contains("only part of")),
            "{missing}: {refused:?}"
        );
    }
    let mut singular = claims.clone();
    singular[3].1 = 1.0;
    assert!(bind(&singular).is_err());
}

/// Ruling 6433 (b): on the card route the kick keeps its preset felt beater,
/// whose receipt replaces the card striker's tip note; the same spec through
/// `drum_with_materials` strikes with the card's hickory sphere.
#[test]
fn g1_the_card_built_kick_keeps_its_preset_beater_and_says_so() {
    let cards = KitMaterials::from_packs_at_temperature(
        HEAD_FILM,
        &[],
        HICKORY,
        "T",
        293.15,
        9e-6,
        0.005,
        0.05,
    )
    .and_then(KitMaterials::with_stated_unknown_film_loss)
    .unwrap();
    let mut spec = kit::kick();
    spec.grid *= 4.0;
    let striker = |piece: &kit::Piece| {
        let body = piece
            .instrument
            .bodies
            .iter()
            .find(|b| b.name == "striker")
            .unwrap();
        (body.mass[0], impact_k(&piece.instrument))
    };
    let kick = kit::kick_with_materials(&spec, &[], standard_gas(), &cards).unwrap();
    let (mass, Law::Power { k, .. }, _) = spec.striker;
    assert_eq!(striker(&kick), (mass, k));
    assert!(
        kick.notes
            .iter()
            .any(|note| note == kit::KICK_BEATER_PRESET_RECEIPT)
    );
    assert!(
        !kick
            .notes
            .iter()
            .any(|note| note.starts_with("striker tip:"))
    );
    assert!(
        kick.notes
            .iter()
            .any(|note| note == FILM_LOSS_STATED_UNKNOWN_RECEIPT)
    );
    let sphere = kit::drum_with_materials(&spec, &[], standard_gas(), &cards).unwrap();
    let Law::Power { k: card_k, .. } = cards.impact_law();
    assert_eq!(striker(&sphere), (cards.striker_mass_kg(), card_k));
    assert!(
        sphere
            .notes
            .iter()
            .any(|note| note == &cards.striker_tip_note())
    );
    assert!(
        !sphere
            .notes
            .iter()
            .any(|note| note == kit::KICK_BEATER_PRESET_RECEIPT)
    );
}

/// .370 DONE-WHEN 4: a data-only film swap (7.5 mil -> 10 mil resonant head of
/// the 12" tom, same card, same tension) changes the computed decay in the
/// predicted direction, inside a band derived before the run.
///
/// The (1,1) shape loses energy only by radiation here (stated-unknown film
/// loss; no uniform-cavity coupling). With T60 = 13.8 (M + Ma) / R, a fixed
/// stiffness K and omega^2 = K / (M + Ma), Ma about constant at small ka and
/// R ~ omega^n: T60_thick / T60_thin = q^(1 + n/2), q = (r M1 + Ma) / (M1 + Ma).
/// r = 2.54e-4 / 1.9e-4 = 1.3368; the thin head's receipt (5d6f99ed7) gives
/// M1 = 1.769e-2 kg and Ma = 4.052e-3 kg, so q = 1.2741. n is 4 in the dipole
/// limit and bounded to [3.4, 4] at ka 0.61; bending changes K by under 0.5 %.
/// Prediction [1.923, 2.068]; with 10 % for the time-domain measurement the
/// declared band is [1.73, 2.28]. It excludes "no change" (1) and the
/// mass-only effect (1.34), and asserts the direction: thicker rings longer.
#[test]
fn g1_a_thicker_film_lengthens_the_radiation_decay_inside_its_predicted_band() {
    let materials = KitMaterials::from_packs_at_temperature(
        HEAD_FILM,
        &[],
        HICKORY,
        "T",
        293.15,
        9e-6,
        0.005,
        0.05,
    )
    .and_then(KitMaterials::with_stated_unknown_film_loss)
    .unwrap();
    let shape = crate::radiation::ModeShape::new(1, 1);
    let t60 = |thickness: f64| {
        let mut spec = kit::tom(12.0);
        spec.grid *= 4.0;
        spec.reso.thickness = thickness;
        let mut ins = kit::drum_with_materials(&spec, &[], standard_gas(), &materials)
            .unwrap()
            .instrument;
        let reso = &mut ins.bodies[1];
        for i in 0..reso.n {
            let [x, y] = reso.xy[i];
            reso.w[i] = 1e-5 * shape.value(x / spec.radius, y / spec.radius, false);
            reso.wp[i] = reso.w[i];
        }
        reso.sync();
        let k = ins.k();
        let mut out = vec![0.0; ins.mics.len()];
        let initial = ins.bodies[1].discrete_energy(k);
        let steps = (0.25 / k) as usize;
        for _ in 0..steps {
            ins.step(&mut out);
        }
        // The batter is not driven: the shape has no net volume on the grid.
        assert!(ins.bodies[0].discrete_energy(k) < 1e-9 * initial);
        let rate = -(ins.bodies[1].discrete_energy(k) / initial).ln() / (steps as f64 * k);
        2.0 * 6.907755 / rate
    };
    let (thin, thick) = (t60(1.9e-4), t60(2.54e-4));
    let ratio = thick / thin;
    assert!(
        (1.73..=2.28).contains(&ratio),
        "T60 thin {thin:.3} s, thick {thick:.3} s: ratio {ratio:.4} outside the predicted [1.73, 2.28]"
    );
}

/// The snare wires' cards as the card route ships them (crate README), each at its own state.
fn snare_wire_binding(m: KitMaterials) -> Result<KitMaterials, String> {
    let steel = include_bytes!("../packs/steel-music-wire-string-heptonstall-ligo-p1400229.fsmatpk");
    let shear = include_bytes!("../packs/music-wire-spring-moduli-fuchs-1968.fsmatpk");
    let shear_point = fs_matdb::QueryPoint::new().with("temperature", 293.0).unwrap();
    m.with_snare_wire(
        steel,
        &crate::materials::temperature_frequency_point("T", 293.15, "frequency", 1000.0).unwrap(),
        shear,
        &shear_point,
    )
}

#[test]
fn g1_card_snare_wires_are_coils_of_the_sourced_steel_and_say_what_is_declared() {
    use crate::coil::{CANOPUS_VINTAGE, COIL_ENDS_NOTE, CoilWire};
    let film = pack("declared film", 1390.0, 4e9, Some(0.38));
    let stick = pack("declared striker", 800.0, 12e9, Some(0.3));
    let cards = snare_wire_binding(binding(&film, &film, &stick)).unwrap();
    let steel = cards.snare_wire().unwrap();
    // Independent arithmetic from the card values (Heptonstall rho 7800 kg/m^3, E 212 GPa;
    // Fuchs G 79.29 GPa) and the maker's wire 0.5 mm, outer 1.35 mm (mean 0.85), pitch 3.4 mm.
    let (d, dm, p) = (0.5e-3_f64, 0.85e-3_f64, 3.4e-3_f64);
    let alpha = (p / (PI * dm)).atan();
    let (sin, cos2) = (alpha.sin(), alpha.cos().powi(2));
    let i = PI * d.powi(4) / 64.0;
    let (bend, twist) = ((1.0 + sin * sin) / (212e9 * i), cos2 / (79.29e9 * 2.0 * i));
    let wire = CoilWire::new(&CANOPUS_VINTAGE, steel, &standard_gas()).unwrap();
    near(wire.bending_n_m2(), 2.0 * sin / (bend + twist));
    near(wire.mass_per_m(), 7800.0 * PI * d * d / 4.0 / sin);
    near(wire.wire_bending_share(), bend / (bend + twist));
    // Material loss: fs-material's exact thermoelastic series of the same state, plus the card's
    // 1.4e-4 excess only inside its 8-3000 Hz band, on the wire-bending share of the mode's
    // bending energy; a pure-tension mode has none.
    let section = fs_material::section::ConcentricElasticSection::try_new(
        0.0,
        vec![fs_material::section::ConcentricMaterialLayer {
            outer_radius_m: d / 2.0,
            material: steel.steel().clone(),
        }],
    )
    .unwrap()
    .with_thermoelastic_bending()
    .unwrap();
    let (inside, outside) = (2.0 * PI * 1000.0, 2.0 * PI * 5000.0);
    let s_w = wire.wire_bending_share();
    near(
        wire.material_rate(inside, 0.4),
        0.5 * inside * 0.4 * s_w * (section.loss_tangent(inside) + 1.4e-4),
    );
    near(
        wire.material_rate(outside, 0.4),
        0.5 * outside * 0.4 * s_w * section.loss_tangent(outside),
    );
    assert_eq!(wire.material_rate(inside, 0.0), 0.0);
    let gas = wire.gas_rate(inside).unwrap();
    assert!(gas > 0.0);
    near(wire.decay_rate(inside, 0.4).unwrap(), wire.material_rate(inside, 0.4) + gas);

    let mut spec = kit::snare();
    spec.grid *= 6.0;
    let piece = kit::drum_with_materials(&spec, &[], standard_gas(), &cards).unwrap();
    let wires: Vec<&Body> = piece
        .instrument
        .bodies
        .iter()
        .filter(|b| b.name.starts_with("wire "))
        .collect();
    assert!(!wires.is_empty());
    for w in &wires {
        assert!(w.modal.is_some(), "{} is not a coil", w.name);
        assert!(w.s0.iter().all(|&s| s == 0.0));
        for &m in &w.mass {
            near(m, wire.mass_per_m() * w.h);
        }
    }
    assert!(piece.notes.iter().any(|n| n.starts_with("snare wires: sourced coils")));
    assert!(piece.notes.iter().any(|n| n == COIL_ENDS_NOTE));
    assert!(!piece.notes.iter().any(|n| n == kit::FITTED_WIRES_RECEIPT));
    // The wire's contact (ruling 7520 L12 (a)): elliptical Hertz of a turn's surface on the film,
    // the wire's radius across it and the helix's radius of curvature plus that along it (steel
    // nu 0.3 assumed), grid / pitch turns per node at one approach; the spec's loss kept.
    let (rm, lead) = (dm / 2.0, p / (2.0 * PI));
    let turn = crate::hertz::EllipticHertz::isotropic((212e9, 0.3), (4e9, 0.38), d / 2.0, (rm * rm + lead * lead) / rm + d / 2.0).unwrap();
    let Law::Power { mu: declared_loss, .. } = spec.wires.unwrap().contact;
    let on_wires = |c: &&Contact| c.b.iter().any(|&(b, _, _)| piece.instrument.bodies[b].name.starts_with("wire "));
    let contacts: Vec<&Contact> = piece.instrument.contacts.iter().filter(on_wires).collect();
    assert!(!contacts.is_empty());
    for c in contacts {
        let Law::Power { k, alpha, mu } = c.law;
        near(k, turn.k() * (spec.grid / p));
        assert_eq!((alpha, mu), (1.5, declared_loss));
    }
    assert!(piece.notes.iter().any(|n| n.starts_with("snare wire contact: computed from the cards (Hertz, elliptical, helix)")));

    // A wire card without the thermal and loss data refuses; absent loss is not a lossless wire.
    let refused = binding(&film, &film, &stick)
        .with_snare_wire(
            &film.to_bytes(),
            &crate::materials::temperature_frequency_point("T", 293.15, "frequency", 1000.0)
                .unwrap(),
            include_bytes!("../packs/music-wire-spring-moduli-fuchs-1968.fsmatpk"),
            &fs_matdb::QueryPoint::new().with("temperature", 293.0).unwrap(),
        )
        .err()
        .unwrap();
    assert!(refused.starts_with("snare wire steel:"), "{refused}");
}

#[test]
fn g1_snare_wires_without_a_sourced_steel_stay_the_fitted_strings_and_say_so() {
    let film = pack("declared film", 1390.0, 4e9, Some(0.38));
    let stick = pack("declared striker", 800.0, 12e9, Some(0.3));
    let mut spec = kit::snare();
    spec.grid *= 6.0;
    // the preset (and so live and GPU), and a card drum bound without the wire cards
    for piece in [
        kit::drum(&spec, &[]).unwrap(),
        kit::drum_with_materials(&spec, &[], standard_gas(), &binding(&film, &film, &stick))
            .unwrap(),
    ] {
        let wires: Vec<&Body> = piece
            .instrument
            .bodies
            .iter()
            .filter(|b| b.name.starts_with("wire "))
            .collect();
        assert!(!wires.is_empty());
        for w in &wires {
            assert!(w.modal.is_none());
            assert!(w.s0.iter().all(|&s| s == 6.907755 / 0.35));
            assert!(w.mass.iter().all(|&m| m == 0.004 * w.h));
        }
        assert!(piece.notes.iter().any(|n| n == kit::FITTED_WIRES_RECEIPT));
        assert!(!piece.notes.iter().any(|n| n.starts_with("snare wires: sourced coils")));
        // the spec's contact law, not a computed one
        let on_wires = |c: &&Contact| c.b.iter().any(|&(b, _, _)| piece.instrument.bodies[b].name.starts_with("wire "));
        let Law::Power { k: k0, alpha: a0, mu: m0 } = spec.wires.unwrap().contact;
        assert!(piece.instrument.contacts.iter().filter(on_wires).all(|c| {
            let Law::Power { k, alpha, mu } = c.law;
            (k, alpha, mu) == (k0, a0, m0)
        }));
        assert!(!piece.notes.iter().any(|n| n.starts_with("snare wire contact:")));
    }
}

#[test]
fn g2_card_drums_keep_no_value_graded_against_a_recording_and_presets_say_theirs_were() {
    // Ruling 7486: 88cccbb78's head ratios, strike points, wires and kick beater; the presets (so
    // live and GPU) keep 96ccf8f1b's and carry FITTED_DRUM_RECEIPT.
    let ratio = |s: &kit::DrumSpec| s.reso.tension / s.batter.tension;
    let declared = [
        (kit::snare(), 2080.0 / 1330.0, (0.06, 0.02)),
        (kit::snare_wires_off(), 2080.0 / 1330.0, (0.06, 0.02)),
        (kit::kick(), 253.0 / 358.0, (0.0, 0.02)),
        (kit::tom(12.0), 950.0 / 1270.0, (0.35 * 0.146, 0.0)),
        (kit::tom(14.0), 175.0 / 330.0, (0.35 * 0.1703, 0.0)),
    ];
    for (preset, head_ratio, strike_xy) in declared {
        assert_eq!(preset.provenance, kit::FITTED_DRUM_RECEIPT);
        let card = kit::card_drum(&preset).unwrap();
        assert!(card.provenance.starts_with("drum (card route, ruling 7486): tensions tuned to the recording's pitch"), "{}", preset.name);
        near(ratio(&card), head_ratio);
        assert_eq!(card.strike_xy, strike_xy, "{}", preset.name);
        assert_eq!((card.radius, card.depth, card.grid, card.mic, card.floor), (preset.radius, preset.depth, preset.grid, preset.mic, preset.floor));
        assert_eq!(card.wires.is_some(), preset.wires.is_some());
        if let Some(w) = card.wires {
            let Law::Power { k, alpha, mu } = w.contact;
            assert_eq!((w.tension, w.preload_per_m, k, alpha, mu), (2.0, 0.25, 2.0e7, 1.5, 2.0));
        }
        if preset.name == "kick 18x18" {
            let (m, Law::Power { k, alpha, mu }, foot) = card.striker;
            assert_eq!((m, k, alpha, mu, foot), (0.12, 1.0e8, 2.5, 0.3, 0.03));
            assert_eq!(card.muffle.map(|m| m.3), Some(25.0));
        }
    }
    assert!(kit::card_drum(&kit::DrumSpec { name: "tom 16", ..kit::tom(14.0) }).is_err());

    let film = pack("declared film", 1390.0, 4e9, Some(0.38));
    let stick = pack("declared striker", 800.0, 12e9, Some(0.3));
    let mut spec = kit::tom(12.0);
    spec.grid *= 4.0;
    let preset = kit::drum(&spec, &[]).unwrap();
    assert!(preset.notes.iter().any(|n| n == kit::FITTED_DRUM_RECEIPT));
    let card = kit::card_drum(&spec).unwrap();
    let piece = kit::drum_with_materials(&card, &[], standard_gas(), &binding(&film, &film, &stick)).unwrap();
    assert!(piece.notes.iter().any(|n| n == card.provenance));
    assert!(!piece.notes.iter().any(|n| n == kit::FITTED_DRUM_RECEIPT));
}

#[test]
fn g2_card_cymbals_and_hats_take_the_declared_interface_and_presets_say_theirs_was_fitted() {
    let coarse = |c: kit::CymbalSpec| kit::CymbalSpec { grid: c.grid * 3.0, ..c };
    let stick_mu = |p: &kit::Piece| {
        let c = p.instrument.contacts.iter().find(|c| c.b.iter().any(|&(b, _, _)| p.instrument.bodies[b].name == "stick")).unwrap();
        let Law::Power { mu, .. } = c.law;
        mu
    };
    let hit = [(0.0, 1.0, 0.6, 0.004)];
    let fitted = kit::cymbal(&coarse(kit::ride()), &hit).unwrap();
    assert!(fitted.notes.iter().any(|n| n == kit::FITTED_CYMBAL_INTERFACE.cymbal_receipt));
    assert!(fitted.notes.iter().any(|n| n.starts_with("ride: fitted to a recording (96ccf8f1b)")));
    assert_eq!(stick_mu(&fitted), 0.05);
    let ride = kit::card_ride();
    assert_eq!((ride.bow_thickness, ride.bell_thickness, ride.t60[1], ride.grid), (1.6e-3, 2.2e-3, (10000.0, 3.0), kit::ride().grid));
    let card = kit::cymbal_with_interface(&coarse(ride), &hit, &kit::DECLARED_CYMBAL_INTERFACE).unwrap();
    assert!(card.notes.iter().any(|n| n == kit::DECLARED_CYMBAL_INTERFACE.cymbal_receipt));
    assert!(!card.notes.iter().any(|n| n.contains("fitted to a recording")));
    assert_eq!(stick_mu(&card), kit::DECLARED_STICK_LOSS_S_M);

    // The hi-hat's rim gaps: the presets' two-lobe waviness, the card route's flat rims.
    let rim_gaps = |p: &kit::Piece| p.instrument.contacts.iter()
        .filter(|c| !c.b.iter().any(|&(b, _, _)| p.instrument.bodies[b].name == "stick"))
        .map(|c| c.gap).collect::<Vec<_>>();
    let spread = |g: &[f64]| g.iter().cloned().fold(f64::MIN, f64::max) - g.iter().cloned().fold(f64::MAX, f64::min);
    let fitted = kit::hihat(0.002, 0.0, &[(0.0, 1.0)]).unwrap();
    assert!(fitted.notes.iter().any(|n| n == kit::FITTED_CYMBAL_INTERFACE.hat_receipt));
    assert_eq!(stick_mu(&fitted), 0.05);
    let card = kit::hihat_with_interface(0.002, 0.0, &[(0.0, 1.0)], kit::DECLARED_HIHAT_RIM_LAW, &kit::DECLARED_CYMBAL_INTERFACE).unwrap();
    assert!(card.notes.iter().any(|n| n == kit::DECLARED_CYMBAL_INTERFACE.hat_receipt));
    assert_eq!(stick_mu(&card), kit::DECLARED_STICK_LOSS_S_M);
    let (fitted_gaps, card_gaps) = (rim_gaps(&fitted), rim_gaps(&card));
    assert!(!card_gaps.is_empty() && card_gaps.len() == fitted_gaps.len());
    // the shells' rest shapes vary the gap a little either way; the waviness adds up to 0.15 mm
    assert!(spread(&fitted_gaps) > spread(&card_gaps) + 1.0e-4, "{} {}", spread(&fitted_gaps), spread(&card_gaps));
    for c in card.instrument.contacts.iter().filter(|c| !c.b.iter().any(|&(b, _, _)| card.instrument.bodies[b].name == "stick")) {
        let Law::Power { mu, .. } = c.law;
        assert_eq!(mu, 3.0);
    }
    // Rim friction (ruling 7520 G1): the presets' graded RIM_FRICTION, the card route's measured
    // tin bronze on itself; both receipts say which.
    let rubs = |p: &kit::Piece| p.instrument.frictions.iter().map(|f| f.mu).collect::<Vec<_>>();
    assert_eq!((kit::RIM_FRICTION, kit::BRONZE_RIM_FRICTION), (0.5, 0.6));
    assert!(rubs(&fitted).len() == fitted_gaps.len() && rubs(&fitted).iter().all(|&m| m == kit::RIM_FRICTION));
    assert!(rubs(&card).len() == card_gaps.len() && rubs(&card).iter().all(|&m| m == kit::BRONZE_RIM_FRICTION));
    assert!(kit::FITTED_CYMBAL_INTERFACE.hat_receipt.contains("rim friction RIM_FRICTION 0.5, picked by grade"));
    assert!(kit::DECLARED_CYMBAL_INTERFACE.hat_receipt.contains("rim friction 0.6 (Buckley, Cu-10Sn self-mated; ruling 7520)"));
}

/// The value between `before` and the next `after` in a receipt line.
fn receipt_number(line: &str, before: &str, after: &str) -> f64 {
    let rest = &line[line.find(before).unwrap_or_else(|| panic!("{before:?} in {line}")) + before.len()..];
    rest[..rest.find(after).unwrap()].parse().unwrap()
}

#[test]
fn g1_the_film_card_curve_gives_each_head_mode_the_loss_tangent_at_its_own_frequency() {
    let cards = || {
        KitMaterials::from_packs_at_temperature(
            HEAD_FILM, &[], HICKORY, "T", 293.15, 9e-6, 0.005, 0.05,
        )
    };
    let film = cards()
        .and_then(|m| m.with_film_loss_curve("T", "frequency", 293.15))
        .unwrap();
    // The card's knots, linear between them, and nothing outside 20 Hz to 20 kHz.
    for head in 0..2 {
        for (hz, knot) in [(20.0, 0.0158), (1000.0, 0.0254), (20000.0, 0.0372)] {
            near(film.film_tan_delta(head, hz).unwrap().unwrap(), knot);
        }
        near(film.film_tan_delta(head, 1125.0).unwrap().unwrap(), 0.5 * (0.0254 + 0.0266));
        for hz in [19.9, 20000.5] {
            assert_eq!(film.film_tan_delta(head, hz), Ok(None), "{hz} Hz");
        }
    }
    assert!(film.film_tan_delta(0, f64::NAN).is_err());
    // The curve is derived at 20 C only; a film loss is one declaration.
    assert!(cards().and_then(|m| m.with_film_loss_curve("T", "frequency", 296.15)).is_err());
    assert!(film.clone().with_stated_unknown_film_loss().is_err());
    let declared = cards().and_then(KitMaterials::with_stated_unknown_film_loss).unwrap();
    assert!(declared.clone().with_film_loss_curve("T", "frequency", 293.15).is_err());
    let plain = pack("declared film", 1390.0, 4e9, Some(0.38)).to_bytes();
    let refused = KitMaterials::from_packs_at_temperature(&plain, &[], HICKORY, "T", 293.15, 9e-6, 0.005, 0.05)
        .and_then(|m| m.with_film_loss_curve("T", "frequency", 293.15))
        .err()
        .unwrap();
    assert!(refused.contains("exactly one loss_tangent_tensile claim, not 0"), "{refused}");
    assert!(film.info_json().contains("\"filmLossCurve\":{\"property\":\"loss_tangent_tensile\",\"temperatureK\":293.15,\"bandHz\":[[20,20000],[20,20000]]"));

    let mut spec = kit::tom(12.0);
    spec.grid *= 4.0;
    let hits = [(0.002, 3.5)];
    let piece = kit::drum_with_materials(&spec, &hits, standard_gas(), &film).unwrap();
    let notes = &piece.notes;
    assert!(notes.contains(&film.film_loss_note().unwrap()));
    assert!(notes.iter().any(|n| n == EDGE_LOSS_STATED_UNKNOWN_RECEIPT));
    assert!(!notes.iter().any(|n| n == FILM_LOSS_STATED_UNKNOWN_RECEIPT));
    assert_eq!(piece.instrument.bodies[0].loss(), (0.0, 0.0));
    // Each retained mode's tan delta is the card's at that mode's loaded
    // frequency [Hz], and its eta is that times its bending share.
    let mut modes = 0;
    for line in notes.iter().filter(|n| n.starts_with("batter (") && n.contains("; film: tan delta")) {
        let hz = receipt_number(line, "loaded ", " Hz");
        let tan_delta = receipt_number(line, "film: tan delta ", ",");
        let share = receipt_number(line, "bending share ", ",");
        let eta = receipt_number(line, "eta ", ",");
        let card = film.film_tan_delta(0, hz).unwrap().unwrap();
        assert!((tan_delta - card).abs() <= 5.1e-5, "{line}: card {card}");
        assert!((eta - tan_delta * share).abs() <= 4e-3 * eta, "{line}");
        modes += 1;
    }
    assert!(modes >= 4, "{modes} batter modes carry a film loss");
    assert!(notes.iter().any(|n| n.starts_with("batter head film loss: tan delta 0.0") && n.contains("at each mode's own loaded frequency")));
    // The film's damping is in the stepped operator: it removes energy, part
    // of what the head radiation's modal damping removes.
    let mut ins = piece.instrument;
    let mut out = vec![0.0; ins.mics.len()];
    for _ in 0..(0.03 / ins.k()) as usize {
        ins.step(&mut out);
    }
    let load = ins.bodies[0].head_radiation().unwrap();
    assert!(load.film_dissipated_j() > 0.0);
    assert!(load.film_dissipated_j() < ins.bodies[0].head_radiated_j());
    let mut stated = kit::drum_with_materials(&spec, &hits, standard_gas(), &declared).unwrap().instrument;
    for _ in 0..(0.03 / stated.k()) as usize {
        stated.step(&mut out);
    }
    assert_eq!(stated.bodies[0].head_radiation().unwrap().film_dissipated_j(), 0.0);

    // The card kick says that no pillow is modelled.
    let mut kick = kit::kick();
    kick.grid *= 6.0;
    let kick = kit::kick_with_materials(&kick, &[], standard_gas(), &film).unwrap();
    assert!(kick.notes.iter().any(|n| n.starts_with(kit::KICK_PILLOW_RECEIPT) && n.contains("muffle dashpot (22 N s/m^3 over a 70 mm radius")));
}
