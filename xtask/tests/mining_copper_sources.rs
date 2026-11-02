//! G0/G3: copper redox source facts survive compilation and refuse other states.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_evidence::NumericalKind;
use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn compile(slug: &str) -> NormalizedPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "frankensim-mining-copper-{}-{nonce}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&directory).unwrap();
    let output = directory.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace().join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", workspace())
        .current_dir(workspace())
        .output()
        .expect("run the actual material compiler");
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes, "{slug}: canonical round trip");
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claim_count(), 3);
    pack
}

fn point(gas: bool) -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            298.15,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap();
    for (axis, value) in [
        ("source-standard-state", 1.0),
        ("source-elemental-reference-janaf", 1.0),
        ("source-pure-species", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
        (
            if gas {
                "source-phase-gas"
            } else {
                "source-phase-solid"
            },
            1.0,
        ),
    ] {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

#[test]
fn g0_g3_copper_redox_references_keep_standard_state_and_source_precision() {
    // Hf and S are literal 298.15 K rows, https://janaf.nist.gov/tables/.
    // Molar masses retain WebBook display precision, including oxide rounding.
    for (slug, gas, molar_mass, hf, entropy) in [
        (
            "copper-ii-oxide-solid-janaf-mining", // Cu-015
            false,
            0.079545,
            -156_063.0,
            42.594,
        ),
        (
            "copper-i-oxide-solid-janaf-mining", // Cu-019
            false,
            0.143091,
            -170_707.0,
            92.360,
        ),
        (
            "copper-solid-janaf-mining", // Cu-002: Hf zero, absolute S nonzero
            false,
            0.063546,
            0.0,
            33.164,
        ),
        (
            "carbon-monoxide-gas-janaf-mining", // C-093
            true,
            0.0280101,
            -110_527.0,
            197.653,
        ),
    ] {
        let pack = compile(slug);
        let query = point(gas);
        for (property, dims, expected) in [
            ("molar-mass", Dims([0, 1, 0, 0, 0, -1]), molar_mass),
            (
                "standard-molar-formation-enthalpy",
                Dims([2, 1, -2, 0, 0, -1]),
                hf,
            ),
            (
                "standard-molar-entropy",
                Dims([2, 1, -2, -1, 0, -1]),
                entropy,
            ),
        ] {
            let key = PropertyKey::new(property, dims);
            let answer = pack
                .claims()
                .query_typed(&key, &query, SelectionPolicy::SingleClaimOnly)
                .unwrap();
            assert!(
                (answer.evidence.value.value - expected).abs() <= 1e-12 * expected.abs().max(1.0),
                "{slug} {property}: {} != {expected}",
                answer.evidence.value.value,
            );
            assert_eq!(answer.evidence.value.dims, dims);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert_eq!(answer.evidence.numerical.kind, NumericalKind::NoClaim);
            assert!(answer.receipt.observation_backed);
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            let claim = pack.claims().claim(answer.receipt.selected).unwrap();
            for (axis, value) in query.axes() {
                assert_eq!(claim.validity.bound(axis), Some((*value, *value)));
                assert_eq!(
                    claim.validity.axis_quantities().get(axis),
                    query.axis_quantities().get(axis),
                );
            }
            for wrong in [
                QueryPoint::new(),
                point(!gas),
                query.clone().with("temperature", 1000.0).unwrap(),
                query.clone().with("pressure", 101_325.0).unwrap(),
                query.clone().with("source-standard-state", 0.0).unwrap(),
                query
                    .clone()
                    .with("source-elemental-reference-janaf", 0.0)
                    .unwrap(),
                query.clone().with("source-pure-species", 0.9).unwrap(),
                query
                    .clone()
                    .with("engineering-reference-model", 0.0)
                    .unwrap(),
                query
                    .clone()
                    .with("source-specimen-qualified", 1.0)
                    .unwrap(),
            ] {
                assert!(
                    pack.claims()
                        .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                        .is_err(),
                    "{slug} {property}: admitted wrong state {wrong:?}",
                );
            }
        }
    }
}
