//! Published JKDWT fits -> actual pack compiler -> typed source -> impact model.
#![deny(unsafe_code)]
use fs_matdb::{NormalizedPack, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::comminution::{
    SpecificGrindingEnergy,
    drop_weight::{ImpactEnergyDomain, resolve_jk_drop_weight},
};
use fs_qty::{Dims, QuantitySpec};
use std::{path::Path, process::Command};

fn point(ore: f64, qualified: f64) -> QueryPoint {
    [
        ("source-bueno-ore", ore),
        ("source-test-temperature-known", 0.0),
        ("source-test-pressure-known", 0.0),
        ("operating-ore-qualified", qualified),
        ("historical-testwork-reference", 1.0),
    ]
    .into_iter()
    .fold(QueryPoint::new(), |point, (axis, value)| {
        point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap()
    })
}
fn compile() -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "bueno-jkdwt-{}-{nonce}.fsmatpk",
        std::process::id()
    ));
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join("data/matdb/seed-v1/bueno-jkdwt-2021/manifest.tsv"))
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
    NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= expected.abs().max(1.0) * 1e-12,
        "{actual} != {expected}"
    );
}

#[test]
fn five_published_ores_drive_jkdwt_energy_response_without_relabelling_as_plant_measurements() {
    let pack = compile();
    assert_eq!(pack.claims().claim_count(), 10);
    assert_eq!(pack.claims().observation_ids().count(), 5);
    // This explicitly declared calculation interval equals the source's overall
    // energy envelope. It does not assert every size/energy pair was tested.
    let domain = ImpactEnergyDomain::new([
        SpecificGrindingEnergy::new(360.0),
        SpecificGrindingEnergy::new(9000.0),
    ])
    .unwrap();
    // Independently transcribed Table 7 JKDWT coefficients and Table 8 apparent
    // calculated t10 (NOT Geopyora measurements), Bueno et al. (2021), DOI
    // 10.1016/j.mineng.2021.106905. Printed coefficients/energies are rounded;
    // agreement to 0.2 percentage points checks published curve scale only.
    for (ore, a_percent, b_per_kwh_t, e_kwh_t, t10_percent) in [
        (1.0, 57.9, 0.55, 0.95, 23.5),
        (2.0, 71.2, 0.49, 1.21, 31.9),
        (3.0, 78.3, 0.40, 1.18, 29.4),
        (4.0, 67.6, 0.89, 0.71, 31.7),
        (5.0, 67.8, 1.31, 1.15, 52.8),
    ] {
        let query = point(ore, 0.0);
        let card = resolve_jk_drop_weight(&pack, &query, SelectionPolicy::SingleClaimOnly, domain)
            .unwrap();
        close(card.asymptotic_t10().value(), a_percent / 100.0);
        close(card.energy_sensitivity().value(), b_per_kwh_t / 3600.0);
        let response = card
            .response(SpecificGrindingEnergy::new(e_kwh_t * 3600.0))
            .unwrap();
        assert!((response.passing_one_tenth().value() * 100.0 - t10_percent).abs() < 0.2);
        close(
            response.passing_one_tenth().value() + response.retained_one_tenth().value(),
            1.0,
        );
        close(
            card.energy_for_t10(response.passing_one_tenth())
                .unwrap()
                .value(),
            e_kwh_t * 3600.0,
        );
        let source = card.calibration_source();
        assert_eq!(source.pack_id(), "bueno-jkdwt-2021");
        assert_eq!(source.pack_hash(), pack.content_hash());
        assert_eq!(source.query(), &query);
        assert_eq!(
            source.claims()[0].observations,
            source.claims()[1].observations
        );
        for answer in source.answers() {
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
        }
        assert!(card.response(SpecificGrindingEnergy::new(100.0)).is_err());
        assert!(
            resolve_jk_drop_weight(
                &pack,
                &point(ore, 1.0),
                SelectionPolicy::SingleClaimOnly,
                domain
            )
            .is_err()
        );
    }
    assert!(
        resolve_jk_drop_weight(
            &pack,
            &point(6.0, 0.0),
            SelectionPolicy::SingleClaimOnly,
            domain
        )
        .is_err()
    );
}
