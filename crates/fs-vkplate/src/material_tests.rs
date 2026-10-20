//! G1/G3: canonical source data reaches actual kit operators and live pressure.
//! Manufactured coefficients establish causal wiring, not a measured instrument.

use super::*;
use crate::materials::KitMaterials;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::{
    Density, Dims, Pressure, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn pack(name: &str, rho: f64, young: f64, poisson: Option<f64>) -> NormalizedMaterialCardPack {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut claims = ClaimSet::new();
    for (key, value, dims) in [
        ("density", Some(rho), Density::DIMS),
        ("young_modulus", Some(young), Pressure::DIMS),
        ("poisson_ratio", poisson, Dims::NONE),
    ] {
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
                    observations: vec![],
                    provenance: Provenance {
                        source: format!("manufactured kit test {name}/{key}"),
                        license: "CC0-1.0".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    let claims = NormalizedPack::new(
        name,
        "manufactured-kit-g1-g3",
        hash_domain("frankensim.kit.material-test", name.as_bytes()),
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
