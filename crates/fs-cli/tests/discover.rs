//! G0/G3 CLI discovery checks with real canonical packs and FrankenSQLite.
//! Authored density values below are synthetic software fixtures.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::Dims;
use fs_qty::semantic::{QuantityKind, QuantitySpec, SemanticType, ValueForm};

static NEXT: AtomicU64 = AtomicU64::new(0);

fn scratch() -> PathBuf {
    let path = std::env::temp_dir().join(format!(
        "fs-cli-discovery-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&path).unwrap();
    path
}

fn pack(id: &str, typed: bool) -> NormalizedMaterialCardPack {
    let provenance = Provenance {
        source: "synthetic CLI density".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "authored coupon".into(),
            method: "software test".into(),
            artifact: hash_bytes(id.as_bytes()),
            caveats: "not a measured material".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let validity = if typed {
        ValidityDomain::unconstrained().with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            200.0,
            400.0,
        )
    } else {
        ValidityDomain::unconstrained().with("temperature", 200.0, 400.0)
    };
    let dims = Dims([-3, 1, 0, 0, 0, 0]);
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("density", dims),
            value: PropertyValue::Scalar {
                value: 1000.0,
                dims,
            },
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    let claims = NormalizedPack::new(
        id,
        "synthetic-discovery-v1",
        hash_bytes(id.as_bytes()),
        "CC0-1.0",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: id.into(),
            phase: "synthetic solid".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
    )
    .unwrap()
}

fn request(domain: &str, properties: &str) -> String {
    format!(
        r#"{{"schema":"frankensim.discovery.v1","target":"materials","properties":{properties},"models":[],"domain":{domain},"selection":"single-claim-only"}}"#
    )
}

const DENSITY: &str = r#"[{"name":"density","unit":"kg/m3","kind":"dimensional"}]"#;
const LOCAL: &str =
    r#"{"mode":"local-state","axes":[{"name":"temperature","kind":"legacy","value":"300 K"}]}"#;

fn run(dir: &Path, source: &str, packs: &[PathBuf], json: bool) -> fs_cli::CommandOutput {
    let path = dir.join("request.json");
    std::fs::write(&path, source).unwrap();
    let mut args = vec!["discover".into(), path.to_str().unwrap().into()];
    args.extend(packs.iter().map(|p| p.to_str().unwrap().to_owned()));
    if json {
        args.push("--json".into());
    }
    fs_cli::run(args)
}

#[test]
fn g0_discover_binary_returns_complete_named_candidates_deterministically() {
    let dir = scratch();
    let first = pack("first", false);
    let second = pack("second", false);
    let paths = [dir.join("first.pack"), dir.join("second.pack")];
    std::fs::write(&paths[0], first.to_bytes()).unwrap();
    std::fs::write(&paths[1], second.to_bytes()).unwrap();
    let source = request(LOCAL, DENSITY);
    let output = run(&dir, &source, &paths, true);
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS, "{}", output.stderr);
    assert_eq!(output.stdout.matches("\"status\":\"complete\"").count(), 2);
    assert!(output.stdout.contains("\"lower_si\":1000"));
    assert!(output.stdout.contains(&first.content_hash().to_hex()));
    assert!(
        output
            .stdout
            .contains("first/synthetic solid/authored rev 0")
    );
    assert!(output.stdout.contains("\"unknown_properties\":[]"));
    let reversed = run(
        &dir,
        &source,
        &[paths[1].clone(), paths[0].clone(), paths[0].clone()],
        true,
    );
    assert_eq!(
        output, reversed,
        "input order and exact duplicates cannot alter discovery"
    );
    let binary = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .arg("--json")
        .arg("discover")
        .arg(dir.join("request.json"))
        .args(&paths)
        .output()
        .unwrap();
    assert!(
        binary.status.success(),
        "{}",
        String::from_utf8_lossy(&binary.stderr)
    );
    assert_eq!(binary.stdout, output.stdout.as_bytes());
    assert!(binary.stderr.is_empty());
    let ratio = run(
        &dir,
        &request(
            LOCAL,
            r#"[{"name":"ratio","unit":"1","kind":"dimensional"}]"#,
        ),
        &paths,
        true,
    );
    assert_eq!(ratio.exit_code, fs_cli::exit::SUCCESS, "{}", ratio.stderr);
    assert!(ratio.stdout.contains("\"unknown_properties\":[\"ratio\"]"));
}

#[test]
fn g0_discover_reports_partial_local_support_and_full_envelope_gaps() {
    let dir = scratch();
    let path = dir.join("card.pack");
    std::fs::write(&path, pack("authored", false).to_bytes()).unwrap();
    let properties = r#"[{"name":"density","unit":"kg/m3","kind":"dimensional"},{"name":"specific-heat-capacity","unit":"J/kg/K","kind":"dimensional"}]"#;
    let output = run(&dir, &request(LOCAL, properties), &[path.clone()], true);
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS, "{}", output.stderr);
    assert!(output.stdout.contains("\"status\":\"partial\""));
    assert!(
        output
            .stdout
            .contains("\"unknown_properties\":[\"specific-heat-capacity\"]")
    );
    let envelope = r#"{"mode":"envelope","axes":[{"name":"temperature","kind":"legacy","lower":"300 K","upper":"450 K"}]}"#;
    let output = run(&dir, &request(envelope, properties), &[path.clone()], true);
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS);
    assert!(output.stdout.contains("\"status\":\"unavailable\""));
    assert!(output.stdout.contains("450"));
    assert!(
        output
            .stdout
            .contains(r#"NoClaimInDomain { property: \"density\", considered: 1 }"#),
        "{}",
        output.stdout
    );
    assert!(
        output
            .stdout
            .contains(r#"point: QueryPoint { axes: {\"temperature\": 450.0}"#)
    );
    let text = run(&dir, &request(LOCAL, properties), &[path], false);
    assert!(text.stdout.contains("status=partial"));
    assert!(
        text.stdout
            .contains("property=density supported lower_si=1000")
    );
    assert!(
        text.stdout
            .contains("unknown_property=specific-heat-capacity")
    );
}

#[test]
fn g0_discover_typed_axes_preserve_conventions_and_convert_declared_units() {
    let dir = scratch();
    let path = dir.join("typed.pack");
    std::fs::write(&path, pack("typed", true).to_bytes()).unwrap();
    let typed = r#"{"mode":"local-state","axes":[{"name":"temperature","kind":"absolute-temperature","value":"26.85 degC"}]}"#;
    let output = run(&dir, &request(typed, DENSITY), &[path.clone()], true);
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS, "{}", output.stderr);
    assert!(output.stdout.contains("\"status\":\"complete\""));
    let legacy = run(&dir, &request(LOCAL, DENSITY), &[path], true);
    assert_eq!(legacy.exit_code, fs_cli::exit::SUCCESS);
    assert!(legacy.stdout.contains("\"status\":\"unavailable\""));
    assert!(legacy.stdout.contains("Quantity"));
}

#[test]
fn g0_discover_request_refuses_unsupported_or_inconsistent_semantics() {
    let dir = scratch();
    let path = dir.join("unused.pack");
    // Request validation precedes file admission; the absent pack must not
    // obscure the exact malformed field in each case.
    for (source, expected) in [
        (
            request(
                LOCAL,
                r#"[{"name":"density","unit":"kg/m3","kind":"dimensional","tensor":{}}]"#,
            ),
            "tensor",
        ),
        (
            request(
                LOCAL,
                r#"[{"name":"density","unit":"kg/m3","kind":"pressure"}]"#,
            ),
            "disagrees",
        ),
        (
            request(
                r#"{"mode":"local-state","axes":[{"name":"T","kind":"legacy","value":"1 K","upper":"2 K"}]}"#,
                DENSITY,
            ),
            "upper",
        ),
        (
            request(
                r#"{"mode":"envelope","axes":[{"name":"T","kind":"legacy","lower":"400 K","upper":"300 K"}]}"#,
                DENSITY,
            ),
            "lower <= upper",
        ),
        (
            request(
                r#"{"mode":"local-state","axes":[{"name":"T","kind":"legacy","value":"1 K"},{"name":"T","kind":"legacy","value":"2 K"}]}"#,
                DENSITY,
            ),
            "duplicate domain axis",
        ),
        (
            request(
                r#"{"mode":"local-state","axes":[{"name":"T","kind":"absolute-temperature","value":"-1 K"}]}"#,
                DENSITY,
            ),
            "FiniteNonnegative",
        ),
    ] {
        let output = run(&dir, &source, &[path.clone()], true);
        assert_eq!(output.exit_code, fs_cli::exit::REFUSED, "{}", output.stderr);
        assert!(output.stderr.contains("cli-discover-request"));
        assert!(
            output.stderr.contains(expected),
            "expected {expected:?}: {}",
            output.stderr
        );
    }
}

#[test]
fn g0_discover_refuses_invalid_packs_and_does_not_infer_material_membership() {
    let dir = scratch();
    let path = dir.join("input.pack");
    let card = pack("unbound", false);
    std::fs::write(&path, card.claims_pack().to_bytes()).unwrap();
    let output = run(&dir, &request(LOCAL, DENSITY), &[path.clone()], true);
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS);
    assert!(output.stdout.contains("\"candidates\":[]"));
    assert!(output.stdout.contains("\"unknown_properties\":[]"));
    std::fs::write(&path, b"not a canonical pack").unwrap();
    let invalid = run(&dir, &request(LOCAL, DENSITY), &[path.clone()], true);
    assert_eq!(invalid.exit_code, fs_cli::exit::REFUSED);
    assert!(invalid.stderr.contains("cli-discover-pack"));
    let too_many = run(&dir, &request(LOCAL, DENSITY), &vec![path; 33], true);
    assert_eq!(too_many.exit_code, fs_cli::exit::INPUT);
    assert!(too_many.stderr.contains("cli-discover-pack-count"));
}

#[test]
fn g0_discover_preserves_model_requirements_and_refuses_malformed_pins() {
    let dir = scratch();
    let path = dir.join("card.pack");
    std::fs::write(&path, pack("no-associated-law", false).to_bytes()).unwrap();
    let source = request(LOCAL, DENSITY).replace(
        "\"models\":[]",
        r#""models":[{"law":"required-law","version":7}]"#,
    );
    let output = run(&dir, &source, &[path.clone()], true);
    assert_eq!(output.exit_code, fs_cli::exit::SUCCESS, "{}", output.stderr);
    assert!(output.stdout.contains("\"status\":\"partial\""));
    assert!(
        output
            .stdout
            .contains("\"unknown_models\":[\"required-law@7\"]")
    );
    assert!(output.stdout.contains("Missing { available_versions: [] }"));
    assert!(output.stdout.contains("\"version\":7"));
    let malformed = source.replace("\"version\":7", "\"version\":7,\"pin\":\"not-a-hash\"");
    let output = run(&dir, &malformed, &[path], true);
    assert_eq!(output.exit_code, fs_cli::exit::REFUSED);
    assert!(
        output
            .stderr
            .contains("pin: expected a 64-digit content hash")
    );
}
