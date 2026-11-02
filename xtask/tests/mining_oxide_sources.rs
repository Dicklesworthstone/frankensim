//! G0/G3 checks of JANAF mineral/gas standard points through the real compiler.
//! These verify source transport and admissible states, not process conversion.
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

const MOLAR_MASS: Dims = Dims([0, 1, 0, 0, 0, -1]);
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const MOLAR_ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const HF: &str = "standard-molar-formation-enthalpy";

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
        "frankensim-mining-oxide-{}-{nonce}-{}",
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

fn point(gas: bool, kelvin: f64) -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            kelvin,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap();
    for (name, value) in [
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
            .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-11 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_g3_janaf_reference_points_preserve_units_phases_and_common_reference() {
    // Hf and S are literal 298.15 K JANAF rows, not Shomate constants or
    // rounded WebBook thermochemistry. Molecular weights come from WebBook.
    for (slug, gas, mass, hf, entropy) in [
        (
            "magnesite-solid-janaf-mining", // C-079
            false,
            0.0843139,
            -1_111_689.0,
            65.854,
        ),
        (
            "magnesium-oxide-solid-janaf-mining", // Mg-008
            false,
            0.0403044,
            -601_241.0,
            26.924,
        ),
        (
            "carbon-dioxide-gas-janaf-mining", // C-095
            true,
            0.0440095,
            -393_522.0,
            213.795,
        ),
        (
            "magnetite-solid-janaf-mining", // Fe-032
            false,
            0.231533,
            -1_120_894.0,
            145.266,
        ),
        (
            "hematite-solid-janaf-mining", // Fe-030
            false,
            0.159688,
            -825_503.0,
            87.400,
        ),
        (
            "oxygen-gas-janaf-mining", // O-029, including reference zero Hf
            true,
            0.0319988,
            0.0,
            205.147,
        ),
    ] {
        let pack = compile(slug);
        let query = point(gas, 298.15);
        for (name, dims, expected) in [
            ("molar-mass", MOLAR_MASS, mass),
            (HF, MOLAR_ENERGY, hf),
            ("standard-molar-entropy", MOLAR_ENTROPY, entropy),
        ] {
            let key = PropertyKey::new(name, dims);
            let answer = pack
                .claims()
                .query_typed(&key, &query, SelectionPolicy::SingleClaimOnly)
                .unwrap_or_else(|error| panic!("{slug} {name}: {error:?}"));
            close(answer.evidence.value.value, expected);
            assert_eq!(answer.evidence.value.dims, dims);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert_eq!(answer.evidence.numerical.kind, NumericalKind::NoClaim);
            assert!(answer.receipt.observation_backed);
            pack.claims().verify_receipt(&answer.receipt).unwrap();

            // Required meaning must be declared by the selected source claim,
            // rather than merely added as otherwise unconstrained query axes.
            let claim = pack.claims().claim(answer.receipt.selected).unwrap();
            for (axis, coordinate) in query.axes() {
                assert_eq!(
                    claim.validity.bound(axis),
                    Some((*coordinate, *coordinate)),
                    "{slug} {name}: {axis}"
                );
                assert_eq!(
                    claim.validity.axis_quantities().get(axis),
                    query.axis_quantities().get(axis)
                );
            }

            for wrong in [
                QueryPoint::new(),
                point(!gas, 298.15),
                point(gas, 300.0),
                query.clone().with("pressure", 101_325.0).unwrap(),
                query.clone().with("source-standard-state", 0.0).unwrap(),
                query
                    .clone()
                    .with("source-elemental-reference-janaf", 0.0)
                    .unwrap(),
                query.clone().with("source-pure-species", 0.98).unwrap(),
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
                    "{slug} {name}: must refuse {wrong:?}"
                );
            }
        }
        // Molar formation energy cannot be requested as specific energy.
        assert!(
            pack.claims()
                .query_typed(
                    &PropertyKey::new(HF, Dims([2, 0, -2, 0, 0, 0])),
                    &query,
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
}
