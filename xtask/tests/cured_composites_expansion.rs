//! G0/G3 source -> actual compiler -> sealed pack -> typed composite queries.
//! Producer coupon values are not a full laminate tensor or a design allowable.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

fn compile(slug: &str, count: usize) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-cured-composites-{}-{nonce}-{slug}",
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
    let bytes = std::fs::read(&out).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let verified = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(verified.to_bytes(), bytes);
    assert_eq!(verified.pack_id(), slug);
    assert_eq!(verified.claims().claims_ordered().count(), count);
    verified
}

fn flag(point: QueryPoint, axis: &str, value: f64) -> QueryPoint {
    point
        .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
        .unwrap()
}

fn base(slug: &str, resin: f64) -> QueryPoint {
    [
        (format!("source-grade-{slug}"), 1.0),
        ("nominal-prepreg-resin-mass-fraction".into(), resin),
        ("source-cure-schedule-known".into(), 0.0),
        ("source-specimens-matched".into(), 0.0),
    ]
    .into_iter()
    .fold(QueryPoint::new(), |p, (axis, value)| flag(p, &axis, value))
}

fn temperature(point: QueryPoint, kelvin: f64) -> QueryPoint {
    point
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            kelvin,
        )
        .unwrap()
}

fn dry_mechanical(point: QueryPoint) -> QueryPoint {
    [
        ("source-dry", 1.0),
        ("source-pressure-known", 0.0),
        ("source-test-rate-known", 0.0),
        ("load-along-0deg", 1.0),
    ]
    .into_iter()
    .fold(point, |p, (axis, value)| flag(p, axis, value))
}

fn key(property: &str, dims: Dims) -> PropertyKey {
    PropertyKey::with_quantity(property, QuantitySpec::dimensional(dims))
}

fn sample(pack: &NormalizedPack, property: &str, dims: Dims, p: &QueryPoint, expected: f64) {
    let answer = pack
        .claims()
        .query_typed(&key(property, dims), p, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|e| panic!("{} {property}: {e:?}", pack.pack_id()));
    let v = &answer.evidence.value;
    assert!((v.value - expected).abs() <= 1e-12 * expected.abs().max(1e-6));
    assert_eq!(v.dims, dims);
    assert!(matches!(v.uncertainty, UncertaintyModel::Unstated));
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
fn g0_g3_carbon_laminates_preserve_fiber_resin_weave_and_direction() {
    // Literal metric entries from Hexcel 8552 EU (2023), pages 2--4.
    for (slug, count, resin, rho, ply_mm, vf, e0, e90, xt, yt, ec, xc, ilss, shear) in [
        (
            "hexcel-8552-as4-ud-2023",
            11,
            0.35,
            1580.0,
            0.130,
            0.5742,
            141.0,
            10.0,
            2207.0,
            81.0,
            Some(128.0),
            Some(1531.0),
            128.0,
            Some(114.0),
        ),
        (
            "hexcel-8552-im7-ud-2023",
            11,
            0.35,
            1570.0,
            0.131,
            0.5770,
            164.0,
            12.0,
            2724.0,
            64.0,
            Some(150.0),
            Some(1690.0),
            137.0,
            Some(120.0),
        ),
        (
            "hexcel-8552-as4-agp193-pw-2023",
            10,
            0.37,
            1570.0,
            0.195,
            0.5529,
            68.0,
            66.0,
            828.0,
            793.0,
            Some(60.0),
            Some(883.0),
            84.0,
            None,
        ),
        (
            "hexcel-8552-as4-agp280-5h-2023",
            10,
            0.37,
            1570.0,
            0.289,
            0.5529,
            67.0,
            66.0,
            876.0,
            800.0,
            Some(64.0),
            Some(924.0),
            79.0,
            None,
        ),
        (
            "hexcel-8552-im7-spg196-p-2023",
            8,
            0.37,
            1560.0,
            0.199,
            0.5557,
            85.0,
            80.0,
            1090.0,
            945.0,
            None,
            None,
            88.0,
            None,
        ),
        (
            "hexcel-8552-im7-spg370-8h-2023",
            8,
            0.37,
            1560.0,
            0.380,
            0.5557,
            86.0,
            81.0,
            1014.0,
            959.0,
            None,
            None,
            90.0,
            None,
        ),
    ] {
        let pack = compile(slug, count);
        let physical = flag(
            flag(base(slug, resin), "source-temperature-known", 0.0),
            "source-moisture-condition-known",
            0.0,
        );
        sample(&pack, "density", RHO, &physical, rho);
        sample(
            &pack,
            "cured-ply-thickness",
            LENGTH,
            &physical,
            ply_mm * 1e-3,
        );
        sample(
            &pack,
            "nominal-fiber-volume-fraction",
            Dims::NONE,
            &physical,
            vf,
        );
        refuses(
            &pack,
            "density",
            RHO,
            &flag(physical.clone(), "source-temperature-known", 1.0),
        );
        refuses(
            &pack,
            "density",
            RHO,
            &flag(physical.clone(), "source-moisture-condition-known", 1.0),
        );

        let p = dry_mechanical(temperature(base(slug, resin), 298.15));
        sample(&pack, "composite-tensile-modulus", PA, &p, e0 * 1e9);
        sample(&pack, "composite-tensile-strength", PA, &p, xt * 1e6);
        sample(
            &pack,
            "interlaminar-short-beam-strength",
            PA,
            &p,
            ilss * 1e6,
        );
        let transverse = flag(p.clone(), "load-along-0deg", 0.0);
        sample(
            &pack,
            "composite-tensile-modulus",
            PA,
            &transverse,
            e90 * 1e9,
        );
        sample(
            &pack,
            "composite-tensile-strength",
            PA,
            &transverse,
            yt * 1e6,
        );
        refuses(&pack, "interlaminar-short-beam-strength", PA, &transverse);
        for (name, value, scale) in [
            ("composite-compression-modulus", ec, 1e9),
            ("composite-compression-strength", xc, 1e6),
            ("in-plane-shear-strength", shear, 1e6),
        ] {
            if let Some(value) = value {
                sample(&pack, name, PA, &p, value * scale);
            } else {
                refuses(&pack, name, PA, &p);
            }
        }
        refuses(&pack, "composite-compression-modulus", PA, &transverse);
        refuses(&pack, "young-modulus", PA, &p);
        refuses(&pack, "shear-modulus", PA, &p);
        refuses(&pack, "poisson-ratio", Dims::NONE, &p);
        refuses(&pack, "composite-tensile-modulus", PA, &QueryPoint::new());
        for bad in [
            temperature(p.clone(), 293.15),
            flag(p.clone(), "source-dry", 0.0),
            flag(p.clone(), "nominal-prepreg-resin-mass-fraction", 0.40),
            flag(p.clone(), "source-specimens-matched", 1.0),
            flag(p.clone(), "source-cure-schedule-known", 1.0),
            flag(p.clone(), &format!("source-grade-{slug}"), 0.0),
        ] {
            refuses(&pack, "composite-tensile-modulus", PA, &bad);
        }
    }
}

#[test]
fn g0_g3_fiberglass_keeps_room_temperature_unknown_and_loading_roles_separate() {
    // Hexcel 913 US FTA-317-MY20, page 4 metric column values.
    for (slug, count, resin, rho, et, xt, ec, xc, ilss, ef, xf) in [
        (
            "hexcel-913-120gl-2020",
            6,
            0.37,
            1830.0,
            21.0,
            489.0,
            None,
            588.0,
            74.0,
            None,
            Some(714.0),
        ),
        (
            "hexcel-913-7781gl-2020",
            8,
            0.37,
            1830.0,
            22.0,
            450.0,
            Some(28.0),
            460.0,
            65.0,
            Some(23.0),
            Some(610.0),
        ),
        (
            "hexcel-913-s2gl-2020",
            6,
            0.33,
            1800.0,
            44.0,
            1400.0,
            Some(46.0),
            1100.0,
            82.0,
            None,
            None,
        ),
    ] {
        let pack = compile(slug, count);
        let physical = flag(
            flag(base(slug, resin), "source-temperature-known", 0.0),
            "source-moisture-condition-known",
            0.0,
        );
        sample(&pack, "density", RHO, &physical, rho);
        refuses(
            &pack,
            "density",
            RHO,
            &flag(physical.clone(), "source-moisture-condition-known", 1.0),
        );
        let p = flag(
            dry_mechanical(flag(base(slug, resin), "source-temperature-known", 0.0)),
            "source-room-temperature",
            1.0,
        );
        for (name, value, scale) in [
            ("composite-tensile-modulus", Some(et), 1e9),
            ("composite-tensile-strength", Some(xt), 1e6),
            ("composite-compression-modulus", ec, 1e9),
            ("composite-compression-strength", Some(xc), 1e6),
            ("interlaminar-short-beam-strength", Some(ilss), 1e6),
            ("composite-flexural-modulus", ef, 1e9),
            ("composite-flexural-strength", xf, 1e6),
        ] {
            if let Some(value) = value {
                sample(&pack, name, PA, &p, value * scale);
            } else {
                refuses(&pack, name, PA, &p);
            }
        }
        refuses(&pack, "young-modulus", PA, &p);
        refuses(&pack, "shear-modulus", PA, &p);
        for bad in [
            flag(p.clone(), "source-temperature-known", 1.0),
            flag(p.clone(), "source-room-temperature", 0.0),
            flag(p.clone(), "source-dry", 0.0),
            flag(p.clone(), "load-along-0deg", 0.0),
            flag(p.clone(), "nominal-prepreg-resin-mass-fraction", 0.40),
            flag(p.clone(), &format!("source-grade-{slug}"), 0.0),
        ] {
            refuses(&pack, "composite-tensile-modulus", PA, &bad);
        }
    }
}
