//! G0/G3/G5: actual manufacturer media packs, typed queries and inventory bases.
//! Source specifications do not become measured ore hardness or plant wear laws.
#![deny(unsafe_code)]

use fs_matdb::{
    MatDbError, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{HardnessScale, QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const MEDIA: &str = "gerdau-forged-grinding-balls-high-hardness";
const PACKING: &str = "molycop-sag-charge-packing-reference";
const MASS: Dims = Dims([0, 1, 0, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const FLAGS: [(&str, f64); 4] = [
    ("source-gerdau-high-hardness-forged-steel", 1.0),
    ("source-product-specification", 1.0),
    ("source-test-temperature-known", 0.0),
    ("source-lot-qualified", 0.0),
];

fn compile(id: &str, claims: usize) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!("{id}-{}-{nonce}.fsmatpk", std::process::id()));
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{id}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(pack.pack_id(), id);
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.claims().claims_ordered().count(), claims);
    pack
}

fn flags(values: &[(&str, f64)]) -> QueryPoint {
    values
        .iter()
        .fold(QueryPoint::new(), |point, &(name, value)| {
            point
                .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
                .unwrap()
        })
}

fn point(size_mm: f64) -> QueryPoint {
    flags(&FLAGS)
        .with_quantity(
            "source-nominal-ball-size-label",
            QuantitySpec::dimensional(LENGTH),
            size_mm / 1000.0,
        )
        .unwrap()
}

fn value(pack: &NormalizedPack, name: &str, quantity: QuantitySpec, point: &QueryPoint) -> f64 {
    let answer = pack
        .claims()
        .query_typed(
            &PropertyKey::with_quantity(name, quantity),
            point,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    answer.evidence.value.value
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-14 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn manufacturer_rows_retain_mass_geometry_and_distinct_hardness_locations() {
    let pack = compile(MEDIA, 7);
    // Independently transcribed metric columns of Gerdau brochure page 2.
    for (size, mass, low_cm, high_cm, out_mm, avh) in [
        (25.0, 0.067, 2.568, 2.634, 1.016, 62.0),
        (32.0, 0.131, 3.206, 3.277, 1.016, 62.0),
        (38.0, 0.226, 3.815, 3.889, 1.016, 62.0),
        (50.0, 0.536, 5.108, 5.202, 1.524, 62.0),
        (64.0, 1.046, 6.431, 6.535, 1.651, 62.0),
        (76.0, 1.809, 7.742, 7.892, 2.032, 61.0),
        (89.0, 2.885, 9.002, 9.187, 2.540, 62.0),
    ] {
        let p = point(size);
        for (name, dims, expected) in [
            ("grinding-ball-nominal-mass", MASS, mass),
            (
                "grinding-ball-average-diameter-specification-minimum",
                LENGTH,
                low_cm / 100.0,
            ),
            (
                "grinding-ball-average-diameter-specification-maximum",
                LENGTH,
                high_cm / 100.0,
            ),
            (
                "grinding-ball-out-of-round-specification-maximum",
                LENGTH,
                out_mm / 1000.0,
            ),
            (
                "grinding-ball-rockwell-c-specification-tolerance",
                Dims::NONE,
                4.0,
            ),
        ] {
            close(
                value(&pack, name, QuantitySpec::dimensional(dims), &p),
                expected,
            );
        }
        close(
            value(
                &pack,
                "grinding-ball-surface-rockwell-c-nominal-specification",
                QuantitySpec::dimensional(Dims::NONE),
                &p,
            ),
            62.0,
        );
        close(
            value(
                &pack,
                "grinding-ball-average-volumetric-rockwell-c-nominal-specification",
                QuantitySpec::dimensional(Dims::NONE),
                &p,
            ),
            avh,
        );
    }
}

#[test]
fn discrete_product_labels_never_interpolate_or_silently_convert_inch_designations() {
    let pack = compile(MEDIA, 7);
    let key = PropertyKey::new("grinding-ball-nominal-mass", MASS);
    for size in [25.4, 40.0, 50.8, 90.0] {
        assert!(
            pack.claims()
                .query_typed(&key, &point(size), SelectionPolicy::SingleClaimOnly)
                .is_err(),
            "unsupported nominal size {size}"
        );
    }
    for (axis, original) in FLAGS {
        let p = point(50.0)
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 1.0 - original)
            .unwrap();
        assert!(matches!(
            pack.claims()
                .query_typed(&key, &p, SelectionPolicy::SingleClaimOnly),
            Err(MatDbError::NoClaimInDomain { .. })
        ));
    }
    assert!(
        pack.claims()
            .query_typed(&key, &QueryPoint::new(), SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
    let wrong_unit = point(50.0)
        .with_quantity(
            "source-nominal-ball-size-label",
            QuantitySpec::dimensional(MASS),
            0.05,
        )
        .unwrap();
    assert!(
        pack.claims()
            .query_typed(&key, &wrong_unit, SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
    // A source specification cannot answer either scale's measured-hardness API.
    for scale in [HardnessScale::RockwellC, HardnessScale::Vickers] {
        let measured = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Hardness(scale),
            ValueForm::Static,
        ));
        assert!(
            pack.claims()
                .query_typed(
                    &PropertyKey::with_quantity(
                        "grinding-ball-surface-rockwell-c-nominal-specification",
                        measured
                    ),
                    &point(50.0),
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
    for name in [
        "density",
        "bond-ball-mill-work-index",
        "archard-wear-coefficient",
    ] {
        assert!(matches!(
            pack.claims().query_typed(
                &PropertyKey::new(name, Dims::NONE),
                &point(50.0),
                SelectionPolicy::SingleClaimOnly
            ),
            Err(MatDbError::UnknownProperty { .. })
        ));
    }
}

#[test]
fn media_inventory_uses_nominal_ball_mass_and_packed_charge_density_on_separate_bases() {
    let balls = compile(MEDIA, 7);
    let count = 2000.0; // Authored inventory, not a measured plant charge.
    let single_mass = value(
        &balls,
        "grinding-ball-nominal-mass",
        QuantitySpec::dimensional(MASS),
        &point(50.0),
    );
    close(count * single_mass, 1072.0);

    let packing = compile(PACKING, 1);
    let p = flags(&[
        ("source-molycop-sag-ball-charge", 1.0),
        ("source-theoretical-packing-reference", 1.0),
        ("source-site-packing-qualified", 0.0),
        ("source-test-temperature-known", 0.0),
    ]);
    let density = value(
        &packing,
        "grinding-media-theoretical-apparent-packing-density",
        QuantitySpec::dimensional(DENSITY),
        &p,
    );
    close(density, 4520.0);
    // Independently declared bulk charge volume includes voids; no grain-density substitution.
    close(30.0 * density, 135600.0);
    assert!(matches!(
        packing.claims().query_typed(
            &PropertyKey::new("density", DENSITY),
            &p,
            SelectionPolicy::SingleClaimOnly
        ),
        Err(MatDbError::UnknownProperty { .. })
    ));
    let site = p
        .with_quantity(
            "source-site-packing-qualified",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap();
    assert!(matches!(
        packing.claims().query_typed(
            &PropertyKey::new(
                "grinding-media-theoretical-apparent-packing-density",
                DENSITY
            ),
            &site,
            SelectionPolicy::SingleClaimOnly
        ),
        Err(MatDbError::NoClaimInDomain { .. })
    ));
}

#[test]
fn media_compilation_is_deterministic_and_catalog_exposes_usable_mass() {
    assert_eq!(compile(MEDIA, 7).to_bytes(), compile(MEDIA, 7).to_bytes());
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let mut args = vec![
        "matdb-catalog".into(),
        "--pack".into(),
        MEDIA.into(),
        "--property".into(),
        "grinding-ball-nominal-mass".into(),
        "--at".into(),
        "source-nominal-ball-size-label=0.05".into(),
    ];
    for (name, number) in FLAGS {
        args.extend(["--at".to_owned(), format!("{name}={number}")]);
    }
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
    assert!(output.contains(MEDIA));
    assert!(output.contains("grinding-ball-nominal-mass"));
}
