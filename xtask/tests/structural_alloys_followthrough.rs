//! Source-level structural steel facts: product minima, measured rate groups,
//! and separately identified reinforcing bars. No nonlinear law is inferred.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{Dims, QuantitySpec};
use std::{path::Path, process::Command};

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);

fn compile(slug: &str, count: usize) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "frankensim-structural-{}-{stamp}-{slug}.fsmatpk",
        std::process::id()
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), count);
    pack
}

fn quantity(point: QueryPoint, name: &str, dims: Dims, value: f64) -> QueryPoint {
    point
        .with_quantity(name, QuantitySpec::dimensional(dims), value)
        .unwrap()
}

fn flags(items: &[(&str, f64)]) -> QueryPoint {
    items.iter().fold(QueryPoint::new(), |p, (name, value)| {
        quantity(p, name, Dims::NONE, *value)
    })
}

fn sample(pack: &NormalizedPack, name: &str, dims: Dims, at: &QueryPoint, expected: f64) -> f64 {
    let result = pack
        .claims()
        .query_typed(
            &PropertyKey::with_quantity(name, QuantitySpec::dimensional(dims)),
            at,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap_or_else(|error| panic!("{} {name}: {error:?}", pack.pack_id()));
    let actual = result.evidence.value.value;
    assert!(
        (actual - expected).abs() <= (expected.abs() * 1e-12).max(1e-20),
        "{name}: {actual} != {expected}"
    );
    assert_eq!(result.evidence.value.dims, dims);
    assert_eq!(
        result.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(result.receipt.observation_backed);
    assert!(!result.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&result.receipt).unwrap();
    actual
}

fn refuses(pack: &NormalizedPack, name: &str, dims: Dims, at: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(
                &PropertyKey::with_quantity(name, QuantitySpec::dimensional(dims)),
                at,
                SelectionPolicy::SingleClaimOnly,
            )
            .is_err(),
        "{} {name} unexpectedly admitted {at:?}",
        pack.pack_id()
    );
}

fn ksi() -> f64 {
    1000.0 * 0.45359237 * 9.80665 / 0.0254_f64.powi(2)
}

#[test]
fn structural_product_minima_keep_geometry_revision_and_elongation_basis() {
    let plate = compile("steel-a572-grade50-millplate-ssab-reference", 4);
    let base = quantity(
        quantity(
            flags(&[
                ("source-specification-minimum", 1.0),
                ("source-mill-plate", 1.0),
                ("source-transverse-tension", 1.0),
                ("source-test-temperature-known", 0.0),
                ("source-additional-elongation-adjustments", 0.0),
            ]),
            "plate-thickness",
            LENGTH,
            0.0127,
        ),
        "plate-width",
        LENGTH,
        1.524,
    );
    sample(
        &plate,
        "tensile-yield-strength-specified-minimum",
        PA,
        &base,
        50.0 * ksi(),
    );
    sample(
        &plate,
        "tensile-strength-specified-minimum",
        PA,
        &base,
        65.0 * ksi(),
    );
    // All admitted mill plates are over 24 in wide: apply the source's two-point
    // reduction and retain separate 2 in and 8 in gauge requirements.
    for (gauge, expected) in [(0.0508, 0.19), (0.2032, 0.16)] {
        let at = quantity(base.clone(), "gauge-length", LENGTH, gauge);
        sample(
            &plate,
            "elongation-at-break-specified-minimum",
            Dims::NONE,
            &at,
            expected,
        );
        refuses(
            &plate,
            "elongation-at-break-specified-minimum",
            Dims::NONE,
            &quantity(
                at,
                "source-additional-elongation-adjustments",
                Dims::NONE,
                1.0,
            ),
        );
    }
    refuses(
        &plate,
        "elongation-at-break-specified-minimum",
        Dims::NONE,
        &base,
    );
    for (axis, dims, value) in [
        ("plate-thickness", LENGTH, 0.102),
        ("plate-width", LENGTH, 0.5),
        ("source-mill-plate", Dims::NONE, 0.0),
        ("source-specification-minimum", Dims::NONE, 0.0),
        ("source-test-temperature-known", Dims::NONE, 1.0),
    ] {
        refuses(
            &plate,
            "tensile-yield-strength-specified-minimum",
            PA,
            &quantity(base.clone(), axis, dims, value),
        );
    }
    refuses(&plate, "tensile-yield-strength-mean", PA, &base);
    refuses(&plate, "young-modulus", PA, &base);

    let hss = compile("steel-a500-grade-c-hss-sti-2021", 2);
    let at = flags(&[
        ("source-specification-minimum", 1.0),
        ("source-a500-grade-c", 1.0),
        ("source-cold-formed-hss", 1.0),
        ("source-astm-revision-year", 2021.0),
        ("source-test-temperature-known", 0.0),
    ]);
    // In the described 2021 revision, 50 ksi also applies to ROUND Grade C HSS.
    sample(
        &hss,
        "tensile-yield-strength-specified-minimum",
        PA,
        &at,
        50.0 * ksi(),
    );
    sample(
        &hss,
        "tensile-strength-specified-minimum",
        PA,
        &at,
        62.0 * ksi(),
    );
    for (axis, value) in [
        ("source-astm-revision-year", 2018.0),
        ("source-a500-grade-c", 0.0),
        ("source-specification-minimum", 0.0),
    ] {
        refuses(
            &hss,
            "tensile-yield-strength-specified-minimum",
            PA,
            &quantity(at.clone(), axis, Dims::NONE, value),
        );
    }
    refuses(&hss, "young-modulus", PA, &at);
}

#[test]
fn measured_plate_and_web_rate_groups_preserve_primary_scatter_and_upper_lower_yield() {
    for (
        slug,
        cohort,
        upper,
        ultimate,
        upper_sd,
        ultimate_sd,
        rate_sd,
        count,
        lower,
        elong,
        area,
    ) in [
        (
            "steel-a572-grade50-rate-murray2020",
            "source-murray2020-a572-cohort",
            [61.5, 65.7, 69.8, 71.7, 83.5],
            [81.9, 85.4, 87.5, 89.8, 97.9],
            [0.27, 1.37, 1.20, 1.24, 0.85],
            [1.43, 0.36, 0.42, 0.55, 0.90],
            [2.3e-6, 9.4e-5, 5.9e-3, 2.1e-2, 1.8e-1],
            [8.0, 5.0, 8.0, 5.0, 6.0],
            58.7,
            0.278,
            0.651,
        ),
        (
            "steel-a992-web-rate-murray2020",
            "source-murray2020-a992-web-cohort",
            [54.0, 60.3, 64.3, 67.6, 78.1],
            [68.9, 71.6, 73.9, 76.2, 83.4],
            [1.38, 1.13, 1.20, 0.70, 1.00],
            [0.48, 0.53, 0.27, 0.18, 0.39],
            [5.3e-6, 1.9e-4, 3.9e-3, 1.6e-2, 3.1e-1],
            [14.0, 5.0, 7.0, 5.0, 10.0],
            51.0,
            0.320,
            0.691,
        ),
    ] {
        let pack = compile(slug, 9);
        let base = flags(&[
            (cohort, 1.0),
            ("source-test-temperature-known", 0.0),
            ("source-room-temperature", 1.0),
        ]);
        for (i, rate) in [0.00002, 0.002, 0.05, 0.2, 2.0].into_iter().enumerate() {
            let at = quantity(base.clone(), "source-target-strain-rate", RATE, rate);
            for (name, expected) in [
                ("tensile-upper-yield-strength-mean", upper[i]),
                ("tensile-strength-mean", ultimate[i]),
                (
                    "tensile-upper-yield-strength-reported-standard-deviation",
                    upper_sd[i],
                ),
                (
                    "tensile-strength-reported-standard-deviation",
                    ultimate_sd[i],
                ),
            ] {
                sample(&pack, name, PA, &at, expected * ksi());
            }
            sample(
                &pack,
                "strain-rate-reported-standard-deviation",
                RATE,
                &at,
                rate_sd[i],
            );
            sample(&pack, "source-test-count", Dims::NONE, &at, count[i]);
            refuses(&pack, "tensile-yield-strength-specified-minimum", PA, &at);
            refuses(
                &pack,
                "tensile-upper-yield-strength-mean",
                PA,
                &quantity(at.clone(), cohort, Dims::NONE, 0.0),
            );
            if i == 0 {
                sample(
                    &pack,
                    "tensile-lower-yield-strength-mean",
                    PA,
                    &at,
                    lower * ksi(),
                );
                sample(
                    &pack,
                    "elongation-after-fracture-mean",
                    Dims::NONE,
                    &at,
                    elong,
                );
                sample(&pack, "reduction-of-area-mean", Dims::NONE, &at, area);
                assert!(upper[i] > lower);
            } else {
                refuses(&pack, "tensile-lower-yield-strength-mean", PA, &at);
                refuses(&pack, "elongation-after-fracture-mean", Dims::NONE, &at);
            }
        }
        for rate in [0.00001, 0.01, 3.0] {
            refuses(
                &pack,
                "tensile-upper-yield-strength-mean",
                PA,
                &quantity(base.clone(), "source-target-strain-rate", RATE, rate),
            );
        }
        // A bare number cannot erase the physical dimensions of strain rate.
        refuses(
            &pack,
            "tensile-upper-yield-strength-mean",
            PA,
            &quantity(base.clone(), "source-target-strain-rate", Dims::NONE, 0.002),
        );
        refuses(
            &pack,
            "young-modulus",
            PA,
            &quantity(base, "source-target-strain-rate", RATE, 0.00002),
        );
    }
}

#[test]
fn reinforcing_bar_groups_keep_manufacturer_size_and_distinct_elongations() {
    for (slug, manufacturer, yield_ksi, ultimate_ksi, uniform, fracture, ratio) in [
        (
            "rebar-a615-grade60-zhao2016-m1",
            1.0,
            [63.0, 63.7, 63.2],
            [97.9, 101.3, 104.0],
            [0.112, 0.107, 0.112],
            [0.165, 0.162, 0.179],
            [1.55, 1.59, 1.64],
        ),
        (
            "rebar-a615-grade60-zhao2016-m4",
            4.0,
            [80.2, 68.1, 63.6],
            [102.9, 95.8, 90.7],
            [0.103, 0.120, 0.121],
            [0.150, 0.183, 0.171],
            [1.28, 1.41, 1.43],
        ),
    ] {
        let pack = compile(slug, 5);
        let base = quantity(
            flags(&[
                ("source-manufacturer-code", manufacturer),
                ("source-zhao2016-a615-grade60-cohort", 1.0),
                ("source-straight-unbent-bar", 1.0),
                ("source-experimental-strain-aging", 0.0),
                ("source-test-temperature-known", 0.0),
                ("source-test-rate-known", 0.0),
            ]),
            "gauge-length",
            LENGTH,
            0.2032,
        );
        for (i, bar) in [5.0, 8.0, 11.0].into_iter().enumerate() {
            let at = quantity(base.clone(), "source-us-bar-number", Dims::NONE, bar);
            sample(
                &pack,
                "tensile-yield-strength-mean",
                PA,
                &at,
                yield_ksi[i] * ksi(),
            );
            sample(
                &pack,
                "tensile-strength-mean",
                PA,
                &at,
                ultimate_ksi[i] * ksi(),
            );
            sample(
                &pack,
                "uniform-elongation-mean",
                Dims::NONE,
                &at,
                uniform[i],
            );
            sample(
                &pack,
                "elongation-at-fracture-mean",
                Dims::NONE,
                &at,
                fracture[i],
            );
            sample(
                &pack,
                "reported-tensile-to-yield-strength-ratio",
                Dims::NONE,
                &at,
                ratio[i],
            );
            assert!(uniform[i] < fracture[i]);
            for (axis, value) in [
                ("source-manufacturer-code", 2.0),
                ("source-straight-unbent-bar", 0.0),
                ("source-experimental-strain-aging", 1.0),
                ("source-test-temperature-known", 1.0),
            ] {
                refuses(
                    &pack,
                    "tensile-yield-strength-mean",
                    PA,
                    &quantity(at.clone(), axis, Dims::NONE, value),
                );
            }
            refuses(
                &pack,
                "elongation-at-fracture-mean",
                Dims::NONE,
                &quantity(at.clone(), "gauge-length", LENGTH, 0.0508),
            );
            refuses(&pack, "tensile-yield-strength-specified-minimum", PA, &at);
            refuses(&pack, "young-modulus", PA, &at);
        }
        refuses(
            &pack,
            "tensile-yield-strength-mean",
            PA,
            &quantity(base, "source-us-bar-number", Dims::NONE, 6.0),
        );
    }
}
