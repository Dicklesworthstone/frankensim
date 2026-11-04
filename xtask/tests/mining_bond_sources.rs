//! Actual primary-source test rows -> compiler -> typed resolver -> Bond physics.
//! Source measurements remain historical specimen tests; energy and throughput
//! below are explicit model calculations, not measured plant performance.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::comminution::{
    BondSizeDomain, WorkIndexMassBasis, source::resolve_bond_ball_mill,
};
use fs_qty::{Dimensionless, Dims, Length, Power, QuantitySpec, Temperature};
use std::{path::Path, process::Command};

const PACK: &str = "joanna-hosco-sgs-bond-2011";
const FLAGS: [(&str, f64); 5] = [
    ("source-test-temperature-known", 0.0),
    ("source-test-pressure-known", 0.0),
    ("source-preparation-detail-known", 0.0),
    ("operating-ore-qualified", 0.0),
    ("historical-testwork-reference", 1.0),
];

fn compile() -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output =
        std::env::temp_dir().join(format!("{PACK}-{}-{nonce}.fsmatpk", std::process::id(),));
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{PACK}/manifest.tsv")))
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
    assert_eq!(pack.pack_id(), PACK);
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.claims().claim_count(), 20);
    assert_eq!(pack.claims().observation_ids().count(), 5);
    pack
}

fn point(row: f64) -> QueryPoint {
    FLAGS
        .into_iter()
        .chain([("source-joanna-composite", row)])
        .fold(QueryPoint::new(), |point, (axis, value)| {
            point
                .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
                .unwrap()
        })
}

fn domain(feed_um: f64, product_um: f64) -> BondSizeDomain {
    BondSizeDomain::new(
        [Length::new(feed_um * 1e-6); 2],
        [Length::new(product_um * 1e-6); 2],
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-13 * expected.abs().max(1.0),
        "{actual} != {expected}",
    );
}

#[test]
fn five_published_composites_keep_their_own_test_points_and_drive_real_bond_cards() {
    let pack = compile();
    // Independent transcription: SGS Joanna report, effective 2011-12-31,
    // section 13.1.3, Table 13.2 p.90. The column and section 2.3 specify metric tonnes.
    for (row, specimen, feed_um, product_um, wi) in [
        (1.0, "Joanna Hosco Comp 1 (S)", 2166.0_f64, 114.0_f64, 12.6),
        (2.0, "Joanna Hosco Comp 2 (S)", 2271.0, 113.0, 14.3),
        (3.0, "Joanna Hosco Comp 3 (N)", 2239.0, 112.0, 12.5),
        (4.0, "Joanna Hosco Comp 4 (S)", 1998.0, 112.0, 12.2),
        (5.0, "Joanna Hosco Comp 5 (S)", 2114.0, 111.0, 12.5),
    ] {
        let card = resolve_bond_ball_mill(
            &pack,
            &point(row),
            WorkIndexMassBasis::MetricTonne,
            SelectionPolicy::SingleClaimOnly,
            domain(feed_um, product_um),
        )
        .unwrap();
        let source = card.calibration_source().unwrap();
        assert_eq!(card.test_context().specimen, specimen);
        assert_eq!(source.pack_hash(), pack.content_hash());
        assert_eq!(source.query(), &point(row));
        assert!(source.observation().caveats.contains("Historical PFS"));
        assert!(
            source
                .observation()
                .provenance
                .source
                .contains("Table 13.2")
        );
        close(card.work_index().source_value(), wi);
        assert_eq!(
            card.work_index().source_basis(),
            WorkIndexMassBasis::MetricTonne
        );
        close(card.test_context().closing_screen.value(), 150e-6);
        let [feed, product] = source.test_sizes();
        close(feed.value(), feed_um * 1e-6);
        close(product.value(), product_um * 1e-6);
        for (claim, answer) in source.claims().iter().zip(source.answers()) {
            assert_eq!(pack.claims().claim(answer.receipt.selected), Some(claim));
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed);
            pack.claims().verify_receipt(&answer.receipt).unwrap();
        }
        // Direct textbook expression is independent of the stabilized SI kernel.
        let expected_kwh_t = 10.0 * wi * (1.0 / product_um.sqrt() - 1.0 / feed_um.sqrt());
        close(
            card.specific_energy(feed, product).unwrap().value(),
            expected_kwh_t * 3600.0,
        );
        close(
            card.dry_throughput(feed, product, Power::new(1e6), Dimensionless::new(0.9))
                .unwrap()
                .value()
                * 3.6,
            900.0 / expected_kwh_t,
        );
        // The exact-point model declaration cannot become a broad operating range.
        assert!(card.specific_energy(feed, Length::new(100e-6)).is_err());
        assert_eq!(card.clone().calibration_source(), Some(source));
    }
    assert!(
        pack.normalizations()
            .iter()
            .any(|n| n.source_basis() == "kW*h/Mg")
    );
}

#[test]
fn missing_context_fractional_specimens_unit_swaps_and_plant_qualification_refuse() {
    let pack = compile();
    let window = domain(2166.0, 114.0);
    let evaluate = |query: &QueryPoint, basis| {
        resolve_bond_ball_mill(
            &pack,
            query,
            basis,
            SelectionPolicy::SingleClaimOnly,
            window,
        )
    };
    assert!(evaluate(&QueryPoint::new(), WorkIndexMassBasis::MetricTonne).is_err());
    for row in [0.0, 1.5, 6.0] {
        assert!(evaluate(&point(row), WorkIndexMassBasis::MetricTonne).is_err());
    }
    for (axis, _) in FLAGS {
        let value = if axis == "historical-testwork-reference" {
            0.0
        } else {
            1.0
        };
        let wrong = point(1.0)
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
        assert!(evaluate(&wrong, WorkIndexMassBasis::MetricTonne).is_err());
    }
    let unsupported_temperature = point(1.0)
        .with_quantity(
            "temperature",
            QuantitySpec::dimensional(Temperature::DIMS),
            298.15,
        )
        .unwrap();
    assert!(evaluate(&unsupported_temperature, WorkIndexMassBasis::MetricTonne).is_err());
    assert!(evaluate(&point(1.0), WorkIndexMassBasis::ShortTon).is_err());
    for (name, dims) in [
        ("bond-ball-mill-work-index-metric-tonne", Dims::NONE),
        ("bond-test-closing-screen", Temperature::DIMS),
        ("ore-recovery", Dims::NONE),
    ] {
        assert!(
            pack.claims()
                .query_typed(
                    &PropertyKey::new(name, dims),
                    &point(1.0),
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
}
