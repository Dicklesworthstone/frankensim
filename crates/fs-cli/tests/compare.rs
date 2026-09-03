//! Integration tests for `frankensim compare` (bead frankensim-rc-root-q61wp.47).
//!
//! Every comparison here runs the real verbs against the tracked reference
//! project: import, `run` (all seven stages), then `compare` on the retained
//! runs. Nothing is stubbed; a compare row that could not be traced to a
//! retained receipt would be a fabricated change.

use fs_cli::{CardPackKind, CardPackSet, RawCardPack, exit, run};

/// The aluminium card the tracked reference project binds its region to.
const REFERENCE_CARD: &str = "2117f2e3d70c07676e8776654aee88a45d4a36cc42852f17bdbb20672609b48e";

const CARD_FIXTURE_DOMAIN: &str = "org.frankensim.fs-cli.tests.compare-card-fixture.v1";
const CONDUCTIVITY_DIMS: fs_qty::Dims = fs_qty::Dims([1, 1, -3, -1, 0, 0]);

fn args(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).to_string()).collect()
}

fn scratch(tag: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("fs-cli-compare-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("scratch directory");
    dir
}

/// A material card pack under the reference project's chemistry key with a
/// caller-chosen conductivity: 167 W/mK is the aluminium twin, 0.04 W/mK the
/// hostile foam twin (Journey A falsifier, bead q61wp.14 item 2).
fn material_pack_bytes(pack_id: &str, conductivity_w_mk: f64) -> Vec<u8> {
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
        ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance,
        UncertaintyModel,
    };
    let provenance = || Provenance {
        source: "compare fixture conductivity table".to_string(),
        license: "CC-BY-4.0; redistribution permitted with attribution".to_string(),
        artifact: Some(fs_blake3::hash_domain(
            CARD_FIXTURE_DOMAIN,
            b"fixture-table",
        )),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "compare fixture coupon".to_string(),
            method: "compare fixture campaign".to_string(),
            artifact: fs_blake3::hash_domain(CARD_FIXTURE_DOMAIN, b"raw-observation"),
            caveats: "fixture value; not a seed-dataset authority".to_string(),
            provenance: provenance(),
        })
        .expect("licensed observation inserts");
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal-conductivity", CONDUCTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: conductivity_w_mk,
                dims: CONDUCTIVITY_DIMS,
            },
            validity: fs_evidence::ValidityDomain::unconstrained().with("T", 200.0, 450.0),
            uncertainty: UncertaintyModel::HalfWidth {
                half_width: 3.0,
                confidence: 0.95,
            },
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance(),
        })
        .expect("conductivity claim inserts");
    let claims_pack = NormalizedPack::new(
        pack_id,
        "frankensim-material-card-pack-compiler-v1",
        fs_blake3::hash_domain(CARD_FIXTURE_DOMAIN, b"source-envelope"),
        "CC-BY-4.0: redistribution permitted with attribution",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .expect("claim pack admits");
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "AA6061".to_string(),
            phase: "wrought".to_string(),
            process: "T6".to_string(),
            revision: 0,
        },
        claims_pack,
    )
    .expect("material-card pack admits")
    .to_bytes()
}

struct Reference {
    fsim: std::path::PathBuf,
    stl: std::path::PathBuf,
    aluminium: std::path::PathBuf,
}

fn reference() -> Reference {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    Reference {
        fsim: root.join("data/reference-project/cooling-reference.fsim"),
        stl: root.join("data/reference-project/plate.stl"),
        aluminium: root.join("data/reference-project/aa6061.fsmcdpk"),
    }
}

fn import(fsim: &std::path::Path, stl: &std::path::Path, ledger: &std::path::Path) {
    let imported = run(args(&[
        "--json",
        "import",
        fsim.to_string_lossy().as_ref(),
        stl.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]));
    assert_eq!(
        imported.exit_code,
        exit::SUCCESS,
        "stderr: {}",
        imported.stderr
    );
}

/// Run all seven stages and return the retained run id.
fn complete_run(
    fsim: &std::path::Path,
    ledger: &std::path::Path,
    pack: &std::path::Path,
) -> String {
    let output = run(args(&[
        "--json",
        "run",
        fsim.to_string_lossy().as_ref(),
        ledger.to_string_lossy().as_ref(),
        "--materials",
        pack.to_string_lossy().as_ref(),
    ]));
    assert_eq!(
        output.exit_code,
        exit::SUCCESS,
        "stdout: {} / stderr: {}",
        output.stdout,
        output.stderr
    );
    assert!(output.stdout.contains("\"stages_completed\":7"));
    let run_id = output
        .stdout
        .split("\"run\":\"")
        .nth(1)
        .and_then(|rest| rest.split('"').next())
        .expect("run id in the result")
        .to_string();
    assert_eq!(run_id.len(), 64);
    run_id
}

fn compare_json(left: &str, right: &str, ledger: &std::path::Path) -> fs_cli::CommandOutput {
    run(args(&[
        "--json",
        "compare",
        left,
        right,
        ledger.to_string_lossy().as_ref(),
    ]))
}

fn json_field<'a>(text: &'a str, key: &str) -> &'a str {
    let marker = format!("\"{key}\":");
    let rest = text
        .split(&marker)
        .nth(1)
        .unwrap_or_else(|| panic!("field `{key}` in {text}"));
    let rest = rest.strip_prefix('"').unwrap_or(rest);
    rest.split(['"', ',', '}', ']'])
        .next()
        .expect("field value")
}

#[test]
fn compare_without_a_ledger_refuses_before_reading_anything() {
    let output = run(args(&["compare", "run_base", "run_opt"]));
    assert_eq!(output.exit_code, exit::INPUT);
    assert!(
        output.stdout.contains("status=refused"),
        "{}",
        output.stdout
    );
    assert!(
        output.stderr.contains("cli-export-ledger-required"),
        "{}",
        output.stderr
    );
    assert!(!output.stdout.contains("changed="));
}

#[test]
fn a_missing_ledger_cannot_mint_a_comparison() {
    let output = run(args(&[
        "--json",
        "compare",
        "same_run",
        "same_run",
        "/definitely/missing/frankensim-ledger.db",
    ]));
    assert_eq!(output.exit_code, exit::INPUT);
    assert!(output.stdout.contains("\"command\":\"compare\""));
    assert!(output.stdout.contains("\"status\":\"refused\""));
    assert!(
        output
            .stderr
            .contains("\"code\":\"cli-export-ledger-missing\""),
        "{}",
        output.stderr
    );
    assert!(!output.stdout.contains("\"changed\":"));
}

#[test]
fn an_unknown_run_refuses_through_the_solve_loader() {
    let reference = reference();
    let dir = scratch("unknown-run");
    let ledger = dir.join("ledger.db");
    import(&reference.fsim, &reference.stl, &ledger);
    let run_id = complete_run(&reference.fsim, &ledger, &reference.aluminium);
    let bogus = "f".repeat(64);
    let output = compare_json(&run_id, &bogus, &ledger);
    assert_eq!(output.exit_code, exit::REFUSED, "{}", output.stderr);
    assert!(output.stdout.contains("\"status\":\"refused\""));
    assert!(!output.stdout.contains("\"changed\":"));
    assert!(output.stderr.contains(&bogus), "{}", output.stderr);
}

/// Self-compare is empty: every stage receipt is the same receipt, every QoI
/// delta is exactly zero, and the summary says so.
#[test]
fn comparing_a_run_with_itself_reports_no_change_anywhere() {
    let reference = reference();
    let dir = scratch("self");
    let ledger = dir.join("ledger.db");
    import(&reference.fsim, &reference.stl, &ledger);
    let run_id = complete_run(&reference.fsim, &ledger, &reference.aluminium);

    let output = compare_json(&run_id, &run_id, &ledger);
    assert_eq!(
        output.exit_code,
        exit::SUCCESS,
        "stdout: {} / stderr: {}",
        output.stdout,
        output.stderr
    );
    let out = &output.stdout;
    assert!(
        out.contains("\"command\":\"compare\",\"status\":\"ok\""),
        "{out}"
    );
    assert!(out.contains("\"changed\":false"), "{out}");
    assert!(
        out.contains("\"summary\":\"identical runs: no differences in any retained receipt\""),
        "{out}"
    );
    assert_eq!(json_field(out, "qoi_count"), "1");
    assert!(out.contains("\"delta\":0,\"rel_delta\":0"), "{out}");
    assert!(out.contains("\"classification\":\"same\""), "{out}");
    assert!(!out.contains("\"classification\":\"changed\""), "{out}");
    assert!(out.contains("\"verdict_changed\":false"), "{out}");
    assert_eq!(
        out.matches("\"status\":\"unchanged (same receipt)\"")
            .count(),
        7,
        "all seven stages are the same receipt: {out}"
    );
    assert!(!out.contains("\"status\":\"changed\""), "{out}");
    assert!(out.contains("\"authority\":\"projection-of-retained-receipts\""));
    assert!(
        out.contains("\"verification\":\"sealed-evidence\""),
        "{out}"
    );

    let text = run(args(&[
        "compare",
        &run_id,
        &run_id,
        ledger.to_string_lossy().as_ref(),
    ]));
    assert_eq!(text.exit_code, exit::SUCCESS);
    assert!(text.stdout.contains("changed=false"), "{}", text.stdout);
    assert!(
        text.stdout
            .contains("stage.report=unchanged (same receipt)"),
        "{}",
        text.stdout
    );
}

/// Write a twin of the reference project whose region binding names `card`
/// instead of the aluminium card. The binding is part of the canonical
/// project, so the twin is a different project hash with the same geometry,
/// requirement, and QoI.
fn rebound_project(dir: &std::path::Path, reference: &Reference, card: &str) -> std::path::PathBuf {
    let source = std::fs::read_to_string(&reference.fsim).expect("reference project reads");
    assert_eq!(
        source.matches(REFERENCE_CARD).count(),
        1,
        "the reference project binds one region to the aluminium card: {source}"
    );
    let twin_source = source.replacen(REFERENCE_CARD, card, 1);
    let twin = dir.join("twin.fsim");
    std::fs::write(&twin, twin_source).expect("twin project writes");
    twin
}

/// The hostile foam twin (Journey A falsifier 2): the reference project
/// rebound to a 0.04 W/mK card under the same chemistry key. Compare must
/// show exactly the material identity change and its consequences — a hotter
/// maximum, a smaller nominal margin — while the verdict stays
/// `indeterminate` on both sides (seven budget terms are NO-DATA, so the
/// composition rule cannot move it; see bead q61wp.14). Stages whose inputs
/// did not change (import-verify, assign) must be reported as unchanged in
/// their inputs even though their receipts name different runs, projects,
/// and import ops.
#[test]
fn the_foam_twin_shows_the_material_identity_change_and_its_consequences() {
    let reference = reference();
    let dir = scratch("foam");
    let ledger = dir.join("ledger.db");
    let foam_bytes = material_pack_bytes("compare-foam", 0.04);
    let foam_card = CardPackSet::admit(vec![RawCardPack {
        kind: CardPackKind::Material,
        source: "foam.fsmcdpk".to_string(),
        bytes: foam_bytes.clone(),
        expect: None,
    }])
    .expect("the foam pack admits")
    .iter()
    .next()
    .expect("one pack")
    .card()
    .to_hex();
    assert_ne!(foam_card, REFERENCE_CARD);
    let foam = dir.join("foam.fsmcdpk");
    std::fs::write(&foam, foam_bytes).expect("foam pack writes");
    let twin = rebound_project(&dir, &reference, &foam_card);

    import(&reference.fsim, &reference.stl, &ledger);
    let aluminium = complete_run(&reference.fsim, &ledger, &reference.aluminium);
    import(&twin, &reference.stl, &ledger);
    let foam_run = complete_run(&twin, &ledger, &foam);
    assert_ne!(
        aluminium, foam_run,
        "the card-pack-set root is bound into run identity"
    );

    let output = compare_json(&aluminium, &foam_run, &ledger);
    assert_eq!(
        output.exit_code,
        exit::SUCCESS,
        "stdout: {} / stderr: {}",
        output.stdout,
        output.stderr
    );
    let out = &output.stdout;
    assert!(out.contains("\"changed\":true"), "{out}");
    assert!(out.contains("\"same_project\":false"), "{out}");
    assert_ne!(
        json_field(out, "project_hash_left"),
        json_field(out, "project_hash_right")
    );

    // Material identity: the one material pack differs by card and identity.
    let materials = out.split("\"materials\":").nth(1).expect("materials block");
    assert!(materials.contains("\"changed\":true"), "{out}");
    let pack = materials.split("\"packs\":[").nth(1).expect("packs");
    assert!(pack.starts_with("{\"kind\":\"material\""), "{out}");
    let card_left = json_field(pack, "card_left");
    let card_right = json_field(pack, "card_right");
    assert_ne!(card_left, card_right, "{out}");
    assert_eq!(card_left.len(), 64);
    assert_eq!(card_right.len(), 64);
    assert!(pack.contains("\"changed\":true"), "{out}");
    assert_eq!(
        materials.matches("\"kind\":").count(),
        1,
        "one pack kind: {out}"
    );

    // QoI: the foam runs hotter, so the delta is strictly positive and the
    // margin shrinks; both sides keep the estimate-only colour.
    let qoi = out.split("\"qoi_diffs\":[").nth(1).expect("qoi_diffs");
    assert!(qoi.starts_with("{\"name\":\"temperature-max\""), "{out}");
    let delta: f64 = json_field(qoi, "delta").parse().expect("delta parses");
    assert!(delta > 0.0, "the foam must run hotter: {out}");
    let nominal_left: f64 = json_field(qoi, "nominal_left").parse().expect("left");
    let nominal_right: f64 = json_field(qoi, "nominal_right").parse().expect("right");
    assert!(nominal_right > nominal_left, "{out}");
    assert!(
        ((nominal_right - nominal_left) - delta).abs() <= 1e-9 * nominal_left.abs(),
        "delta is right minus left: {out}"
    );
    assert!(qoi.contains("\"color_left\":\"estimated\",\"color_right\":\"estimated\",\"color_evolution\":\"same\""), "{out}");
    assert!(
        qoi.contains("\"identity_same\":false,\"classification\":\"changed\""),
        "{out}"
    );

    // Requirement: the verdict does not flip (both indeterminate) but the
    // nominal margin moves with the maximum.
    let requirement = out
        .split("\"requirements\":[")
        .nth(1)
        .expect("requirements");
    assert!(
        requirement.contains("\"outcome_left\":\"indeterminate\",\"outcome_right\":\"indeterminate\",\"verdict_changed\":false"),
        "{out}"
    );
    let margin_left: f64 = json_field(requirement, "nominal_margin_left")
        .parse()
        .expect("l");
    let margin_right: f64 = json_field(requirement, "nominal_margin_right")
        .parse()
        .expect("r");
    assert!(margin_right < margin_left, "the margin shrinks: {out}");
    assert!(
        ((margin_left - margin_right) - delta).abs() <= 1e-9 * margin_left.abs(),
        "the margin moves by exactly the maximum's delta: {out}"
    );

    // Budget terms: the default fidelity measures none, so all eight stay
    // NO-DATA on both sides and none is reported as changed.
    let terms = out
        .split("\"budget_terms\":[")
        .nth(1)
        .expect("budget_terms");
    let terms = terms.split("],\"materials\"").next().expect("terms block");
    assert_eq!(terms.matches("\"kind\":").count(), 8, "{out}");
    assert_eq!(
        terms
            .matches("\"state_left\":\"no-data\",\"state_right\":\"no-data\"")
            .count(),
        8,
        "{out}"
    );
    assert!(!terms.contains("\"changed\":true"), "{out}");

    // Stages: geometry import, assignment, and the flow network saw the same
    // inputs (their receipts differ only in binding keys: the run, the
    // project hash, and the twin's own import op). The flow network is
    // material-independent — the air path sees geometry and duty, not the
    // solid's conductivity — and the comparison reports that from the
    // receipt bytes. Every stage from material-resolve on that consumes the
    // card changed.
    let stages = out.split("\"stages\":[").nth(1).expect("stages");
    let stage_row = |name: &str| -> &str {
        let rest = stages
            .split(&format!("\"stage\":\"{name}\""))
            .nth(1)
            .expect("stage row");
        // Each row ends with its `differing_keys` array; keep that closing
        // bracket so the key list can be matched exactly.
        let end = rest.find("]}").expect("row end");
        &rest[..=end]
    };
    for (unchanged, expected_keys) in [
        ("import-verify", "[\"run\",\"project_hash\",\"import_op\"]"),
        ("assign", "[\"run\"]"),
        ("flow-network", "[\"run\"]"),
    ] {
        let row = stage_row(unchanged);
        assert!(
            row.contains("\"status\":\"unchanged (same inputs; differs only by binding keys)\""),
            "{unchanged}: {row}"
        );
        assert!(
            row.contains(&format!("\"differing_keys\":{expected_keys}")),
            "{unchanged} differing keys: {row}"
        );
    }
    for changed in ["material-resolve", "conduction", "qoi", "report"] {
        let row = stage_row(changed);
        assert!(row.contains("\"status\":\"changed\""), "{changed}: {row}");
        assert!(
            !row.contains("\"differing_keys\":[\"run\"]"),
            "{changed} must differ in more than its run binding: {row}"
        );
    }
    // The summary embeds the requirement's JSON identity (escaped quotes),
    // so it is matched inside the raw output rather than through
    // `json_field`, which stops at the first quote.
    let summary = out
        .split("\"summary\":\"")
        .nth(1)
        .expect("summary")
        .split("\",\"qoi_count\"")
        .next()
        .expect("summary end");
    assert!(
        summary.starts_with("4 of 7 stages changed; project hash "),
        "{summary}"
    );
    assert!(summary.contains("; temperature-max "), "{summary}");
    // The requirement id is its JSON identity document, so the verdict
    // fragment names it as `verdict {...} stays indeterminate (margin ...)`.
    assert!(summary.contains("verdict {"), "{summary}");
    assert!(
        summary.contains("stays indeterminate (margin "),
        "{summary}"
    );
    assert!(summary.contains("material card "), "{summary}");

    // The verb is symmetric in what it reports: swapping the operands negates
    // the delta and swaps the card columns.
    let reverse = compare_json(&foam_run, &aluminium, &ledger);
    assert_eq!(reverse.exit_code, exit::SUCCESS);
    let reverse_delta: f64 = json_field(
        reverse.stdout.split("\"qoi_diffs\":[").nth(1).expect("qoi"),
        "delta",
    )
    .parse()
    .expect("delta");
    assert_eq!(reverse_delta, -delta, "{}", reverse.stdout);
    let reverse_pack = reverse.stdout.split("\"packs\":[").nth(1).expect("packs");
    assert_eq!(json_field(reverse_pack, "card_left"), card_right);
    assert_eq!(json_field(reverse_pack, "card_right"), card_left);
}

/// A renamed project with the same geometry, card, and requirement is a
/// different canonical project but a comparable one: compare reports the
/// project-hash change and shows that nothing physical moved (bit-identical
/// QoI value, same card, no verdict or budget change). The refusal for runs
/// with no common requirement is pinned at the unit level in
/// `compare::tests`, because every project the driver can complete today
/// carries the same single QoI and requirement identity.
#[test]
fn a_renamed_project_compares_and_shows_only_the_identity_change() {
    let reference = reference();
    let dir = scratch("renamed");
    let ledger = dir.join("ledger.db");
    let source = std::fs::read_to_string(&reference.fsim).expect("reference project reads");
    assert!(source.contains(":name \"solve-reference\""), "{source}");
    let twin_source = source.replacen(":name \"solve-reference\"", ":name \"solve-twin\"", 1);
    let twin = dir.join("twin.fsim");
    std::fs::write(&twin, twin_source).expect("twin project writes");

    import(&reference.fsim, &reference.stl, &ledger);
    let left = complete_run(&reference.fsim, &ledger, &reference.aluminium);
    import(&twin, &reference.stl, &ledger);
    let right = complete_run(&twin, &ledger, &reference.aluminium);
    assert_ne!(left, right);

    let output = compare_json(&left, &right, &ledger);
    assert_eq!(
        output.exit_code,
        exit::SUCCESS,
        "stdout: {} / stderr: {}",
        output.stdout,
        output.stderr
    );
    let out = &output.stdout;
    assert!(out.contains("\"same_project\":false"), "{out}");
    assert!(out.contains("\"changed\":true"), "{out}");
    let summary = json_field(out, "summary");
    assert!(summary.contains("project hash "), "{summary}");
    assert!(!summary.contains("temperature-max "), "{summary}");
    assert!(!summary.contains("card "), "{summary}");

    // Same physics: the value is bit-identical, so the row is `same` even
    // though its identity (which carries the project binding) differs.
    let qoi = out.split("\"qoi_diffs\":[").nth(1).expect("qoi_diffs");
    assert_eq!(
        json_field(qoi, "nominal_left"),
        json_field(qoi, "nominal_right")
    );
    assert!(qoi.contains("\"delta\":0,\"rel_delta\":0"), "{out}");
    assert!(
        qoi.contains("\"identity_same\":false,\"classification\":\"same\""),
        "{out}"
    );
    assert!(out.contains("\"verdict_changed\":false"), "{out}");

    // Same card: the material block is unchanged.
    let materials = out.split("\"materials\":").nth(1).expect("materials block");
    let pack = materials.split("\"packs\":[").nth(1).expect("packs");
    assert_eq!(
        json_field(pack, "card_left"),
        json_field(pack, "card_right")
    );
    assert_eq!(json_field(pack, "card_left"), REFERENCE_CARD);
    assert!(pack.contains("\"changed\":false"), "{out}");
}
