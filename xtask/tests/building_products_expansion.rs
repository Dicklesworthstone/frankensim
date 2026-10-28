//! G0/G3 source -> real compiler -> typed query tests for building products.
//! Product specifications and nominal estimates are not experimental validation.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::path::Path;
use std::process::Command;

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);

fn compile(slug: &str, count: usize) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let time = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-building-products-{}-{time}-{slug}",
        std::process::id()
    ));
    std::fs::create_dir(&dir).unwrap();
    let out = dir.join("material.fsmatpk");
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&out)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let bytes = std::fs::read(out).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let verified = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(verified.to_bytes(), bytes);
    assert_eq!(verified.pack_id(), slug);
    assert_eq!(verified.claims().claims_ordered().count(), count);
    verified
}

fn flags(values: &[(&str, f64)]) -> QueryPoint {
    values.iter().fold(QueryPoint::new(), |p, (axis, value)| {
        p.with_quantity(*axis, QuantitySpec::dimensional(Dims::NONE), *value)
            .unwrap()
    })
}

fn unknown(product: &str) -> QueryPoint {
    flags(&[
        (product, 1.0),
        ("source-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
    ])
}

fn quantity(p: QueryPoint, name: &str, dims: Dims, value: f64) -> QueryPoint {
    p.with_quantity(name, QuantitySpec::dimensional(dims), value)
        .unwrap()
}

fn temperature(p: QueryPoint, kelvin: f64) -> QueryPoint {
    p.with_quantity(
        "temperature",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        )),
        kelvin,
    )
    .unwrap()
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    let spec = if name == "thermal-conductivity" {
        assert_eq!(dims, K);
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ))
    } else {
        QuantitySpec::dimensional(dims)
    };
    PropertyKey::with_quantity(name, spec)
}

fn sample(pack: &NormalizedPack, property: &str, dims: Dims, p: &QueryPoint, expected: f64) {
    let answer = pack
        .claims()
        .query_typed(&key(property, dims), p, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|e| panic!("{} {property}: {e:?}", pack.pack_id()));
    assert!((answer.evidence.value.value - expected).abs() <= 1e-11 * expected.abs().max(1.0));
    assert_eq!(answer.evidence.value.dims, dims);
    assert!(matches!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    ));
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
}

fn refuses(pack: &NormalizedPack, property: &str, dims: Dims, p: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(&key(property, dims), p, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "{} must refuse {property} at {p:?}",
        pack.pack_id()
    );
}

#[test]
fn g0_g3_mdf_hdf_and_particleboard_preserve_specification_and_thickness() {
    for (slug, count, product, thickness, modulus, density) in [
        (
            "finsa-fibrapan-bio-18mm-2026",
            6,
            "source-product-fibrapan-bio-2026",
            0.018,
            2.4e9,
            None,
        ),
        (
            "finsa-fibrapan-model-30mm-2026",
            5,
            "source-product-fibrapan-model-2026",
            0.030,
            3.6e9,
            Some(900.0),
        ),
        (
            "finsa-fimapan-af-38mm-2026",
            5,
            "source-product-fimapan-af-2026",
            0.038,
            1.25e9,
            Some(600.0),
        ),
    ] {
        let pack = compile(slug, count);
        let p = quantity(unknown(product), "nominal-thickness", LENGTH, thickness);
        sample(&pack, "bending-moe-specification", PA, &p, modulus);
        refuses(&pack, "young-modulus", PA, &p);
        refuses(
            &pack,
            "bending-moe-specification",
            PA,
            &p.clone()
                .with("nominal-thickness", thickness + 0.001)
                .unwrap(),
        );
        refuses(
            &pack,
            "bending-moe-specification",
            PA,
            &p.clone().with(product, 0.0).unwrap(),
        );
        refuses(
            &pack,
            "bending-moe-specification",
            PA,
            &p.clone()
                .with("source-moisture-condition-known", 1.0)
                .unwrap(),
        );
        if let Some(density) = density {
            sample(&pack, "density", RHO, &p, density);
        } else {
            sample(&pack, "density-guide-low", RHO, &p, 700.0);
            sample(&pack, "density-guide-high", RHO, &p, 720.0);
            refuses(&pack, "density", RHO, &p);
        }
    }
}

#[test]
fn g0_g3_melamine_face_does_not_turn_core_specifications_into_laminate_properties() {
    let pack = compile("kronospan-melamine-faced-p2-18mm-2024", 4);
    let p = quantity(
        unknown("source-product-kronospan-mf-p2-2024"),
        "nominal-thickness",
        LENGTH,
        0.018,
    );
    sample(&pack, "core-bending-moe-specification", PA, &p, 1.6e9);
    sample(&pack, "core-bending-mor-specification", PA, &p, 11e6);
    refuses(&pack, "young-modulus", PA, &p);
    refuses(&pack, "bending-moe", PA, &p);
    let p = quantity(p, "pearl-effect-decor", Dims::NONE, 0.0);
    let force = Dims([1, 1, -2, 0, 0, 0]);
    sample(&pack, "surface-scratch-force-lower-limit", force, &p, 1.5);
    refuses(
        &pack,
        "surface-scratch-force-lower-limit",
        force,
        &p.with("pearl-effect-decor", 1.0).unwrap(),
    );
}

#[test]
fn g0_g3_gypsum_and_cement_boards_keep_thermal_test_conditions() {
    let gypsum = compile("fermacell-gypsum-fibreboard-reference", 6);
    let p = unknown("source-product-fermacell-gypsum-standard");
    sample(&gypsum, "density", RHO, &p, 1150.0);
    sample(&gypsum, "density-declared-tolerance", RHO, &p, 50.0);
    sample(&gypsum, "thermal-conductivity", K, &p, 0.32);
    sample(&gypsum, "specific-heat", CP, &p, 1100.0);
    sample(
        &gypsum,
        "thermal-expansion",
        Dims([0, 0, 0, -1, 0, 0]),
        &p,
        1e-5,
    );
    assert_eq!(
        gypsum.to_bytes(),
        compile("fermacell-gypsum-fibreboard-reference", 6).to_bytes()
    );

    let cement = compile("fermacell-powerpanel-h2o-reference", 6);
    let product = "source-product-fermacell-powerpanel-h2o";
    let p = quantity(unknown(product), "nominal-thickness", LENGTH, 0.0125);
    sample(&cement, "density", RHO, &p, 1000.0);
    sample(&cement, "specific-heat", CP, &p, 1000.0);
    sample(&cement, "source-elastic-modulus", PA, &p, 4.2e9);
    refuses(&cement, "young-modulus", PA, &p);
    refuses(&cement, "thermal-conductivity", K, &p);
    let dry = temperature(
        quantity(
            flags(&[
                (product, 1.0),
                ("source-temperature-known", 1.0),
                ("source-moisture-condition-known", 1.0),
                ("dry-condition", 1.0),
            ]),
            "nominal-thickness",
            LENGTH,
            0.0125,
        ),
        283.15,
    );
    sample(&cement, "thermal-conductivity", K, &dry, 0.17);
    refuses(
        &cement,
        "thermal-conductivity",
        K,
        &dry.clone().with("temperature", 293.15).unwrap(),
    );
    refuses(
        &cement,
        "thermal-conductivity",
        K,
        &dry.with("dry-condition", 0.0).unwrap(),
    );
}

#[test]
fn g0_g3_mortar_curing_is_discrete_and_fresh_density_cannot_supply_hardened_mass() {
    let pack = compile("sika-monotop-412n-nz-2026", 9);
    let base = unknown("source-product-monotop-412n-nz-2026");
    let fresh = quantity(base.clone(), "fresh-mortar", Dims::NONE, 1.0);
    sample(&pack, "fresh-mortar-density", RHO, &fresh, 2100.0);
    refuses(&pack, "density", RHO, &fresh);
    let cured = quantity(
        quantity(base, "fresh-mortar", Dims::NONE, 0.0),
        "source-cure-age-known",
        Dims::NONE,
        1.0,
    );
    for (days, compression, bending) in [(1.0, 17e6, 4e6), (7.0, 40e6, 6e6), (28.0, 55e6, 8e6)] {
        let p = quantity(cured.clone(), "cure-age", TIME, days * 86400.0);
        sample(&pack, "compressive-strength", PA, &p, compression);
        sample(&pack, "flexural-strength", PA, &p, bending);
        refuses(&pack, "fresh-mortar-density", RHO, &p);
    }
    refuses(
        &pack,
        "compressive-strength",
        PA,
        &quantity(cured, "cure-age", TIME, 14.0 * 86400.0),
    );
}

#[test]
fn g0_g3_aerated_concretes_preserve_formulation_and_dry_ten_degree_conductivity() {
    for (slug, product, rho, tolerance, k, strength) in [
        (
            "bauroc-classic-aac-2020",
            "source-product-bauroc-classic-2020",
            425.0,
            30.0,
            0.10,
            3e6,
        ),
        (
            "bauroc-ecoterm-plus-aac-2020",
            "source-product-bauroc-ecoterm-plus-2020",
            300.0,
            25.0,
            0.072,
            1.8e6,
        ),
    ] {
        let pack = compile(slug, 4);
        sample(
            &pack,
            "compressive-strength",
            PA,
            &unknown(product),
            strength,
        );
        let dry = flags(&[
            (product, 1.0),
            ("source-temperature-known", 0.0),
            ("source-moisture-condition-known", 1.0),
            ("dry-condition", 1.0),
        ]);
        sample(&pack, "density", RHO, &dry, rho);
        sample(&pack, "density-declared-tolerance", RHO, &dry, tolerance);
        let thermal = temperature(dry.with("source-temperature-known", 1.0).unwrap(), 283.15);
        sample(&pack, "thermal-conductivity", K, &thermal, k);
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &thermal.clone().with("temperature", 284.15).unwrap(),
        );
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &thermal.with("dry-condition", 0.0).unwrap(),
        );
    }
}
