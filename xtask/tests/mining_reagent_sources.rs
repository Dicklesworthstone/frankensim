//! G0/G1/G3 checks of mining reagent source transport, phase and reference state.
//! These check the real compiler/query path, not plant chemistry validation.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_evidence::NumericalKind;
use fs_matdb::{
    MaterialAnswer, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const MOLAR_MASS: Dims = Dims([0, 1, 0, 0, 0, -1]);
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const MOLAR_ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const SPECIFIC_HEAT: Dims = Dims([2, 0, -2, -1, 0, 0]);
const HF: &str = "standard-molar-formation-enthalpy";
const CP: &str = "specific-heat-capacity";

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
        "frankensim-mining-reagent-{}-{nonce}-{}",
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
    assert_eq!(pack.claims().claims_ordered().count(), 4);
    pack
}

fn point(liquid: bool, kelvin: f64) -> QueryPoint {
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
            if liquid {
                "source-phase-liquid"
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

fn sample(pack: &NormalizedPack, name: &str, dims: Dims, query: &QueryPoint) -> MaterialAnswer {
    let answer = pack
        .claims()
        .query_typed(
            &PropertyKey::new(name, dims),
            query,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap_or_else(|error| panic!("{} {name}: {error:?}", pack.pack_id()));
    assert_eq!(answer.evidence.value.dims, dims);
    // Unstated source uncertainty cannot acquire a numerical accuracy claim.
    assert_eq!(answer.evidence.numerical.kind, NumericalKind::NoClaim);
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    answer
}

fn refuses(pack: &NormalizedPack, name: &str, dims: Dims, query: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(
                &PropertyKey::new(name, dims),
                query,
                SelectionPolicy::SingleClaimOnly,
            )
            .is_err(),
        "{} {name}: must refuse {query:?}",
        pack.pack_id()
    );
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-11 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_g3_standard_points_preserve_sources_and_refuse_wrong_phases_and_references() {
    let mut formation = Vec::new();
    for (slug, liquid, mass, hf, entropy) in [
        (
            "calcium-oxide-solid-nist-mining",
            false,
            0.056077,
            -635_090.0,
            38.19,
        ),
        (
            "calcium-hydroxide-solid-nist-mining",
            false,
            0.074093,
            -986_090.0,
            83.36,
        ),
        (
            "water-liquid-janaf-mining",
            true,
            0.0180153,
            -285_830.0,
            69.950,
        ),
        (
            "sodium-hydroxide-solid-nist-mining",
            false,
            0.0399971,
            -425_930.0,
            64.46,
        ),
        (
            "sodium-cyanide-solid-nist-mining",
            false,
            0.0490072,
            -90_710.0,
            118.49,
        ),
        (
            "sulfuric-acid-liquid-janaf-mining",
            true,
            0.098078,
            -813_989.0,
            156.895,
        ),
    ] {
        let pack = compile(slug);
        let query = point(liquid, 298.15);
        close(
            sample(&pack, "molar-mass", MOLAR_MASS, &query)
                .evidence
                .value
                .value,
            mass,
        );
        let enthalpy = sample(&pack, HF, MOLAR_ENERGY, &query).evidence.value.value;
        close(enthalpy, hf);
        formation.push(enthalpy);
        close(
            sample(&pack, "standard-molar-entropy", MOLAR_ENTROPY, &query)
                .evidence
                .value
                .value,
            entropy,
        );
        for wrong in [
            QueryPoint::new(),
            point(!liquid, 298.15),
            point(liquid, 300.0),
            query.clone().with("pressure", 101_325.0).unwrap(),
            query.clone().with("source-standard-state", 0.0).unwrap(),
            query
                .clone()
                .with("source-elemental-reference-janaf", 0.0)
                .unwrap(),
            query.clone().with("source-pure-species", 0.98).unwrap(),
            query
                .clone()
                .with("source-specimen-qualified", 1.0)
                .unwrap(),
        ] {
            refuses(&pack, HF, MOLAR_ENERGY, &wrong);
        }
        refuses(&pack, HF, Dims([2, 0, -2, 0, 0, 0]), &query);
        refuses(
            &pack,
            "aqueous-ion-formation-enthalpy",
            MOLAR_ENERGY,
            &query,
        );
    }
    // Independent Hess-law result from the three named source values. This
    // says nothing about slurry concentration, conversion rate or plant heat loss.
    close(formation[1] - formation[0] - formation[2], -65_170.0);
}

#[test]
fn g1_g3_sensible_heat_retains_shomate_branches_and_literal_liquid_points() {
    for (slug, mass, coefficients, upper) in [
        (
            "calcium-oxide-solid-nist-mining",
            0.056077,
            [49.95403, 4.887916, -0.352056, 0.046187, -0.825097],
            1000.0,
        ),
        (
            "calcium-hydroxide-solid-nist-mining",
            0.074093,
            [130.8253, -82.69216, 122.7690, -50.39210, -2.513146],
            500.0,
        ),
        (
            "sodium-hydroxide-solid-nist-mining",
            0.0399971,
            [419.4837, -1717.754, 2953.573, -1597.221, -6.046884],
            500.0,
        ),
        (
            "sodium-cyanide-solid-nist-mining",
            0.0490072,
            [68.14732, 1.091197, 0.010951, -0.010567, 0.013628],
            500.0,
        ),
    ] {
        let pack = compile(slug);
        let [a, b, c, d, e] = coefficients;
        let equation = |kelvin: f64| {
            let t = kelvin / 1000.0;
            (a + b * t + c * t.powi(2) + d * t.powi(3) + e / t.powi(2)) / mass
        };
        for temperature in [298.15, 300.0, 400.0, upper] {
            close(
                sample(&pack, CP, SPECIFIC_HEAT, &point(false, temperature))
                    .evidence
                    .value
                    .value,
                equation(temperature),
            );
        }
        // Interior probes test the committed interpolation rather than merely
        // checking copied knots. This numerical discrepancy is not source uncertainty.
        let mut temperature = 312.5;
        while temperature < upper {
            let actual = sample(&pack, CP, SPECIFIC_HEAT, &point(false, temperature))
                .evidence
                .value
                .value;
            assert!((actual / equation(temperature) - 1.0).abs() < 0.0013);
            temperature += 25.0;
        }
        refuses(&pack, CP, SPECIFIC_HEAT, &point(false, 298.0));
        refuses(&pack, CP, SPECIFIC_HEAT, &point(false, upper + 0.01));
    }
    for (slug, mass, temperatures, molar_cp, upper) in [
        (
            "water-liquid-janaf-mining",
            0.0180153,
            [298.15, 300.0, 360.0],
            [75.351, 75.349, 75.679],
            360.0,
        ),
        (
            "sulfuric-acid-liquid-janaf-mining",
            0.098078,
            [298.15, 300.0, 400.0],
            [138.584, 138.940, 158.205],
            400.0,
        ),
    ] {
        let pack = compile(slug);
        for (temperature, cp) in temperatures.into_iter().zip(molar_cp) {
            close(
                sample(&pack, CP, SPECIFIC_HEAT, &point(true, temperature))
                    .evidence
                    .value
                    .value,
                cp / mass,
            );
        }
        refuses(&pack, CP, SPECIFIC_HEAT, &point(true, upper + 0.01));
        refuses(&pack, CP, SPECIFIC_HEAT, &point(false, 298.15));
    }
}
