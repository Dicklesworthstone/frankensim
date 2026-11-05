//! G0/G3: source-backed flotation products and explicit dosing/phase bases.
#![deny(unsafe_code)]

use fs_matdb::{
    MatDbError, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const ID: &str = "mibc-arkema-flotation-reference";
const PAX: &str = "pax-ixom-solid-reference";
const PAX22: &str = "pax22-ixom-solution-reference";
const SIBX: &str = "sibx-ixom-solid-reference";
const COPPER_SULFATE: &str = "copper-sulfate-pentahydrate-ixom-reference";
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);

fn compile() -> NormalizedPack {
    compile_pack(ID, 5)
}

fn compile_pack(id: &str, claims: usize) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let out = std::env::temp_dir().join(format!("{id}-{}-{nonce}.fsmatpk", std::process::id()));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{id}/manifest.tsv")))
        .arg("--out")
        .arg(&out)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(out).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), id);
    assert_eq!(pack.claims().claims_ordered().count(), claims);
    pack
}

fn product() -> QueryPoint {
    product_for("source-arkema-mibc-product")
}

fn product_for(source: &str) -> QueryPoint {
    [
        (source, 1.0),
        ("engineering-reference-model", 1.0),
        ("source-lot-qualified", 0.0),
    ]
    .into_iter()
    .fold(QueryPoint::new(), |p, (name, value)| {
        p.with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap()
    })
}

fn flag(point: QueryPoint, axis: &str, value: f64) -> QueryPoint {
    point
        .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
        .unwrap()
}

fn sibx_solubility(temperature: f64) -> QueryPoint {
    let point = product_for("source-ixom-sibx-solid-product")
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature,
        )
        .unwrap();
    [
        ("source-test-pressure-known", 0.0),
        ("source-solubility-volume-basis-known", 0.0),
        ("source-approximate-value", 1.0),
    ]
    .into_iter()
    .fold(point, |p, (axis, value)| flag(p, axis, value))
}

fn copper_product() -> QueryPoint {
    flag(
        product_for("source-ixom-copper-sulfate-pentahydrate-product"),
        "source-hydration-water-moles-per-formula-unit",
        5.0,
    )
}

fn copper_relative_density() -> QueryPoint {
    [
        "source-test-temperature-known",
        "source-test-pressure-known",
        "source-relative-density-reference-known",
    ]
    .into_iter()
    .fold(copper_product(), |p, axis| flag(p, axis, 0.0))
}

fn liquid(temperature: f64) -> QueryPoint {
    product()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "source-test-pressure-known",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap()
}

fn value(pack: &NormalizedPack, name: &str, dims: Dims, point: &QueryPoint) -> f64 {
    let result = pack
        .claims()
        .query_typed(
            &PropertyKey::new(name, dims),
            point,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert_eq!(
        result.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(result.receipt.observation_backed);
    assert!(!result.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&result.receipt).unwrap();
    result.evidence.value.value
}

#[test]
fn source_density_purity_and_directional_solubilities_keep_original_bases() {
    let pack = compile();
    let p = liquid(293.15);
    for (name, dims, expected) in [
        ("density-reported-range-minimum", DENSITY, 806.0),
        ("density-reported-range-maximum", DENSITY, 808.0),
        (
            "water-rich-phase-mibc-saturation-mass-fraction",
            Dims::NONE,
            0.017,
        ),
        (
            "mibc-rich-phase-water-saturation-mass-fraction",
            Dims::NONE,
            0.058,
        ),
    ] {
        let actual = value(&pack, name, dims, &p);
        assert!((actual - expected).abs() <= 1e-13 * expected.abs().max(1.0));
    }
    assert_eq!(
        value(
            &pack,
            "mibc-product-purity-mass-fraction-minimum",
            Dims::NONE,
            &product()
        ),
        0.985
    );
    // An authored volume flow of 1 L/min and 60 dry tonnes/h: product kg/t,
    // not a recommended frother dose or a forecast of mineral recovery.
    let volume_m3_s = 0.001 / 60.0;
    let dry_kg_s = 60000.0 / 3600.0;
    for (name, expected_kg_t) in [
        ("density-reported-range-minimum", 0.806),
        ("density-reported-range-maximum", 0.808),
    ] {
        let product_kg_t = volume_m3_s * value(&pack, name, DENSITY, &p) / dry_kg_s * 1000.0;
        assert!((product_kg_t - expected_kg_t).abs() < 1e-14);
    }
}

#[test]
fn frother_reference_refuses_invented_temperature_lot_and_recovery_properties() {
    let pack = compile();
    let key = PropertyKey::new("density-reported-range-minimum", DENSITY);
    for temperature in [293.14, 293.16, 298.15] {
        assert!(matches!(
            pack.claims()
                .query_typed(&key, &liquid(temperature), SelectionPolicy::SingleClaimOnly),
            Err(MatDbError::NoClaimInDomain { .. })
        ));
    }
    for (name, wrong) in [
        ("source-arkema-mibc-product", 0.0),
        ("engineering-reference-model", 0.0),
        ("source-lot-qualified", 1.0),
        ("source-test-pressure-known", 1.0),
    ] {
        let p = liquid(293.15)
            .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), wrong)
            .unwrap();
        assert!(matches!(
            pack.claims()
                .query_typed(&key, &p, SelectionPolicy::SingleClaimOnly),
            Err(MatDbError::NoClaimInDomain { .. })
        ));
    }
    assert!(
        pack.claims()
            .query_typed(&key, &product(), SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
    for name in [
        "density",
        "flotation-rate-constant",
        "surface-tension",
        "standard-molar-formation-enthalpy",
    ] {
        assert!(matches!(
            pack.claims().query_typed(
                &PropertyKey::new(name, Dims::NONE),
                &liquid(293.15),
                SelectionPolicy::SingleClaimOnly
            ),
            Err(MatDbError::UnknownProperty { .. })
        ));
    }
    assert!(
        pack.claims()
            .query_typed(
                &PropertyKey::new("water-rich-phase-mibc-saturation-mass-fraction", DENSITY),
                &liquid(293.15),
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
}

#[test]
fn frother_pack_is_deterministic_and_discoverable_by_ordinary_catalog() {
    assert_eq!(compile().to_bytes(), compile().to_bytes());
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let args = [
        "matdb-catalog",
        "--pack",
        ID,
        "--property",
        "density-reported-range-minimum",
        "--at",
        "temperature=293.15",
        "--at",
        "source-arkema-mibc-product=1",
        "--at",
        "engineering-reference-model=1",
        "--at",
        "source-lot-qualified=0",
        "--at",
        "source-test-pressure-known=0",
    ];
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(args)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let output = String::from_utf8(run.stdout).unwrap();
    assert!(output.contains(ID));
    assert!(output.contains("density-reported-range-minimum"));
}

#[test]
fn collector_and_activator_products_preserve_composition_bases() {
    let pax = compile_pack(PAX, 2);
    let pax22 = compile_pack(PAX22, 1);
    let sibx = compile_pack(SIBX, 3);
    let copper = compile_pack(COPPER_SULFATE, 2);
    let solid = product_for("source-ixom-pax-solid-product");
    let solution = flag(
        product_for("source-ixom-pax22-solution-product"),
        "source-product-label-nominal",
        1.0,
    );
    let sibx_product = product_for("source-ixom-sibx-solid-product");
    let solid_pax = value(
        &pax,
        "pax-product-component-mass-fraction-reported",
        Dims::NONE,
        &solid,
    );
    let liquid_pax = value(
        &pax22,
        "pax-solution-component-mass-fraction-nominal",
        Dims::NONE,
        &solution,
    );
    assert_eq!(solid_pax, 0.9);
    assert_eq!(liquid_pax, 0.22);
    assert_eq!(
        value(
            &pax,
            "potassium-hydroxide-product-component-mass-fraction-reported",
            Dims::NONE,
            &solid,
        ),
        0.015
    );
    for (property, expected) in [
        ("sibx-product-component-mass-fraction-minimum", 0.9),
        (
            "sodium-hydroxide-product-component-mass-fraction-maximum",
            0.01,
        ),
    ] {
        assert_eq!(value(&sibx, property, Dims::NONE, &sibx_product), expected);
    }
    assert_eq!(
        value(
            &copper,
            "copper-sulfate-pentahydrate-product-component-mass-fraction-minimum",
            Dims::NONE,
            &copper_product(),
        ),
        0.98
    );

    // Authored 20 kg product/h at 100 dry tonnes/h; this is reference mass
    // accounting, not an operating recommendation or a recovery prediction.
    let product_kg_per_t = 20.0 / 100.0;
    assert!((product_kg_per_t * solid_pax - 0.18).abs() < 1e-14);
    assert!((product_kg_per_t * liquid_pax - 0.044).abs() < 1e-14);
    // Replacing a mass of solid product with the same mass of 22% solution
    // cannot preserve active collector mass. A volume substitution needs
    // independent density evidence, deliberately absent from these packs.
    assert!((solid_pax / liquid_pax - 4.090909090909091).abs() < 1e-14);
    for (pack, point) in [(&pax, &solid), (&pax22, &solution)] {
        assert!(matches!(
            pack.claims().query_typed(
                &PropertyKey::new("density", DENSITY),
                point,
                SelectionPolicy::SingleClaimOnly,
            ),
            Err(MatDbError::UnknownProperty { .. })
        ));
    }
}

#[test]
fn sibx_solubility_retains_approximation_temperature_and_unknown_volume_basis() {
    let pack = compile_pack(SIBX, 3);
    let name = "sibx-water-solubility-reported-approximate-mass-per-volume";
    let key = PropertyKey::new(name, DENSITY);
    assert_eq!(value(&pack, name, DENSITY, &sibx_solubility(293.15)), 500.0);
    for temperature in [293.14, 293.16, 298.15] {
        assert!(matches!(
            pack.claims().query_typed(
                &key,
                &sibx_solubility(temperature),
                SelectionPolicy::SingleClaimOnly,
            ),
            Err(MatDbError::NoClaimInDomain { .. })
        ));
    }
    for (axis, wrong) in [
        ("source-solubility-volume-basis-known", 1.0),
        ("source-approximate-value", 0.0),
        ("source-test-pressure-known", 1.0),
        ("source-lot-qualified", 1.0),
    ] {
        assert!(matches!(
            pack.claims().query_typed(
                &key,
                &flag(sibx_solubility(293.15), axis, wrong),
                SelectionPolicy::SingleClaimOnly,
            ),
            Err(MatDbError::NoClaimInDomain { .. })
        ));
    }
    assert!(
        pack.claims()
            .query_typed(
                &key,
                &product_for("source-ixom-sibx-solid-product"),
                SelectionPolicy::SingleClaimOnly,
            )
            .is_err()
    );
    assert!(
        pack.claims()
            .query_typed(
                &PropertyKey::new(name, Dims::NONE),
                &sibx_solubility(293.15),
                SelectionPolicy::SingleClaimOnly,
            )
            .is_err()
    );
}

#[test]
fn copper_reference_keeps_hydrate_and_relative_density_separate_from_mass_density() {
    let pack = compile_pack(COPPER_SULFATE, 2);
    let relative = "relative-density-reported";
    assert_eq!(
        value(&pack, relative, Dims::NONE, &copper_relative_density()),
        2.284
    );
    for (axis, wrong) in [
        ("source-hydration-water-moles-per-formula-unit", 0.0),
        ("source-hydration-water-moles-per-formula-unit", 1.0),
        ("source-test-temperature-known", 1.0),
        ("source-test-pressure-known", 1.0),
        ("source-relative-density-reference-known", 1.0),
    ] {
        assert!(matches!(
            pack.claims().query_typed(
                &PropertyKey::new(relative, Dims::NONE),
                &flag(copper_relative_density(), axis, wrong),
                SelectionPolicy::SingleClaimOnly,
            ),
            Err(MatDbError::NoClaimInDomain { .. })
        ));
    }
    assert!(
        pack.claims()
            .query_typed(
                &PropertyKey::new(relative, DENSITY),
                &copper_relative_density(),
                SelectionPolicy::SingleClaimOnly,
            )
            .is_err()
    );
    for name in [
        "density",
        "bulk-density",
        "copper-product-component-mass-fraction",
        "standard-molar-formation-enthalpy",
        "flotation-activation-rate-constant",
    ] {
        assert!(matches!(
            pack.claims().query_typed(
                &PropertyKey::new(name, Dims::NONE),
                &copper_relative_density(),
                SelectionPolicy::SingleClaimOnly,
            ),
            Err(MatDbError::UnknownProperty { .. })
        ));
    }
}

#[test]
fn collector_products_refuse_other_products_lot_certification_and_nominal_upgrades() {
    for (id, count, source, property) in [
        (
            PAX,
            2,
            "source-ixom-pax-solid-product",
            "pax-product-component-mass-fraction-reported",
        ),
        (
            SIBX,
            3,
            "source-ixom-sibx-solid-product",
            "sibx-product-component-mass-fraction-minimum",
        ),
    ] {
        let pack = compile_pack(id, count);
        for (axis, wrong) in [
            (source, 0.0),
            ("source-lot-qualified", 1.0),
            ("engineering-reference-model", 0.0),
        ] {
            assert!(matches!(
                pack.claims().query_typed(
                    &PropertyKey::new(property, Dims::NONE),
                    &flag(product_for(source), axis, wrong),
                    SelectionPolicy::SingleClaimOnly,
                ),
                Err(MatDbError::NoClaimInDomain { .. })
            ));
        }
    }
    let pack = compile_pack(PAX22, 1);
    let key = PropertyKey::new("pax-solution-component-mass-fraction-nominal", Dims::NONE);
    assert!(
        pack.claims()
            .query_typed(
                &key,
                &product_for("source-ixom-pax-solid-product"),
                SelectionPolicy::SingleClaimOnly,
            )
            .is_err()
    );
    assert!(matches!(
        pack.claims().query_typed(
            &key,
            &flag(
                product_for("source-ixom-pax22-solution-product"),
                "source-product-label-nominal",
                0.0,
            ),
            SelectionPolicy::SingleClaimOnly,
        ),
        Err(MatDbError::NoClaimInDomain { .. })
    ));
}

#[test]
fn collector_and_activator_packs_are_deterministic_and_catalog_discoverable() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    for (id, count, source, property, extra) in [
        (
            PAX,
            2,
            "source-ixom-pax-solid-product",
            "pax-product-component-mass-fraction-reported",
            None,
        ),
        (
            PAX22,
            1,
            "source-ixom-pax22-solution-product",
            "pax-solution-component-mass-fraction-nominal",
            Some("source-product-label-nominal=1"),
        ),
        (
            SIBX,
            3,
            "source-ixom-sibx-solid-product",
            "sibx-product-component-mass-fraction-minimum",
            None,
        ),
        (
            COPPER_SULFATE,
            2,
            "source-ixom-copper-sulfate-pentahydrate-product",
            "copper-sulfate-pentahydrate-product-component-mass-fraction-minimum",
            Some("source-hydration-water-moles-per-formula-unit=5"),
        ),
    ] {
        assert_eq!(
            compile_pack(id, count).to_bytes(),
            compile_pack(id, count).to_bytes()
        );
        let mut command = Command::new(env!("CARGO_BIN_EXE_xtask"));
        command
            .args([
                "matdb-catalog",
                "--pack",
                id,
                "--property",
                property,
                "--at",
            ])
            .arg(format!("{source}=1"))
            .args([
                "--at",
                "engineering-reference-model=1",
                "--at",
                "source-lot-qualified=0",
            ]);
        if let Some(extra) = extra {
            command.args(["--at", extra]);
        }
        let run = command
            .env("CARGO_WORKSPACE_DIR", root)
            .current_dir(root)
            .output()
            .unwrap();
        assert!(
            run.status.success(),
            "{id}: {}",
            String::from_utf8_lossy(&run.stderr)
        );
        let output = String::from_utf8(run.stdout).unwrap();
        assert!(output.contains(id));
        assert!(output.contains(property));
    }
}
