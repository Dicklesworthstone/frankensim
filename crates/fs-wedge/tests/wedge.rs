//! Battery for go-to-market wedge selection (addendum Proposal 7). Verifies
//! historical-score supersession, evidence-complete measured inputs, workspace
//! evidence drift, candidate rankings, and the cycle-time kill criterion.

use fs_blake3::hash_domain;
use fs_wedge::historical::{
    COMPARISON_EVIDENCE_BUNDLE_IDENTITY_DOMAIN, COMPARISON_EVIDENCE_MANIFEST_IDENTITY_DOMAIN,
};
use fs_wedge::{
    BaselineProvenance, CHT_BASELINE, COMPARISON_EVIDENCE_SNAPSHOT, ComparisonCandidate,
    DEFAULT_FACTOR_WEIGHTS, EvidenceKind, HistoricalEvidenceError, HistoricalEvidenceSnapshot,
    IncumbentStep, InputAxis, KillCriterionError, KillVerdict, Measurement,
    RETIRED_PLACEHOLDER_BASELINE, Readiness, STRONG_THRESHOLD, ScoreUse, ScoringError,
    ScoringFactor, WEDGE_DOCTRINE, WedgeCriterion, audit, chosen_wedge, comparison_candidates,
    comparison_evidence_bundle, comparison_evidence_manifest, default_recommendation,
    four_criteria, measured_inputs_for, measured_wedge_inputs, render_comparison_report,
    score_candidates, to_json, verify_comparison_evidence, verify_default_comparison_evidence,
    verticals,
};
use std::path::Path;

#[test]
fn the_historical_beachhead_is_conjugate_heat_transfer() {
    let w = chosen_wedge();
    assert_eq!(w.name, "conjugate-heat-transfer");
    assert_eq!(w.rank, 1);
    // it exercises incremental re-solve (2), adjoints (1), the ladder (3),
    // and the evidence package (12).
    assert!(w.exercises.contains(&"2") && w.exercises.contains(&"3"));
}

#[test]
fn historical_scores_are_preserved_but_superseded_for_decisions() {
    let w = chosen_wedge();
    assert_eq!(w.score_use, ScoreUse::SupersededForDecisionUse);
    for c in four_criteria() {
        // Replay retains the plan's values; the decision API refuses them.
        assert!(
            w.score(c) >= STRONG_THRESHOLD,
            "historical {} score changed on {}",
            w.name,
            c.label()
        );
        assert_eq!(w.decision_score(c), None);
    }
    assert!(w.weakest_criterion_score() >= STRONG_THRESHOLD);
    assert!(
        verticals()
            .iter()
            .all(|vertical| !vertical.score_use.permits_decision())
    );
}

#[test]
fn every_candidate_has_complete_measured_inputs_on_all_four_axes() {
    let inputs = measured_wedge_inputs();
    assert_eq!(inputs.len(), verticals().len());
    for vertical in verticals() {
        let measured = measured_inputs_for(vertical.name)
            .unwrap_or_else(|| panic!("missing measured inputs for {}", vertical.name));
        assert!(measured.is_complete(), "incomplete: {}", measured.vertical);
        assert!(!measured.kernels.is_empty());
        assert!(!measured.validation_data.is_empty());
        assert!(!measured.cad_burden.is_empty());
        assert!(!measured.compute_cost.is_empty());
        for measurement in measured.measurements() {
            assert!(measurement.is_complete(), "{measurement:?}");
            assert!(!measurement.evidence.is_empty());
            assert!(
                measurement
                    .evidence
                    .iter()
                    .all(|pointer| pointer.is_complete())
            );
        }
    }
}

#[test]
fn absent_inputs_cannot_carry_strong_scores() {
    for inputs in measured_wedge_inputs() {
        for measurement in inputs.measurements() {
            assert!(
                measurement.score <= measurement.readiness.score_ceiling(),
                "{} has score {} above {:?} ceiling {}",
                inputs.vertical,
                measurement.score,
                measurement.readiness,
                measurement.readiness.score_ceiling()
            );
            if measurement.readiness == Readiness::Absent {
                assert!(
                    measurement.score < STRONG_THRESHOLD,
                    "absent {} input scored {}",
                    inputs.vertical,
                    measurement.score
                );
            }
        }
    }
}

fn check_workspace_measurement(
    root: &Path,
    vertical: &str,
    axis: &str,
    label: &str,
    measurement: Measurement,
) -> Vec<String> {
    let mut failures = Vec::new();
    for pointer in measurement
        .evidence
        .iter()
        .filter(|pointer| pointer.kind == EvidenceKind::WorkspacePath)
    {
        let path = root.join(pointer.reference);
        let result = std::fs::read_to_string(&path);
        let (passed, detail) = match result {
            Ok(contents) if contents.contains(pointer.locator) => {
                (true, "marker-found".to_string())
            }
            Ok(_) => (false, format!("missing marker {:?}", pointer.locator)),
            Err(error) => (false, format!("read failed: {error}")),
        };
        eprintln!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            if passed { "PASS" } else { "FAIL" },
            vertical,
            axis,
            label,
            pointer.reference,
            detail
        );
        if !passed {
            failures.push(format!(
                "{} {} {}: {} ({detail})",
                vertical, axis, label, pointer.reference
            ));
        }
    }
    failures
}

#[derive(Debug, Clone, Copy)]
struct KernelProbe {
    path: &'static str,
    marker: &'static str,
}

const fn kernel_probe(path: &'static str, marker: &'static str) -> KernelProbe {
    KernelProbe { path, marker }
}

#[allow(clippy::too_many_lines)] // One exhaustive capability-to-probe table is easiest to audit.
fn kernel_probes(vertical: &str, capability: &str) -> Option<Vec<KernelProbe>> {
    match (vertical, capability) {
        ("conjugate-heat-transfer", "steady-conduction-fem") => Some(vec![kernel_probe(
            "crates/fs-conduction/src/solve.rs",
            "pub struct ConductionSolver",
        )]),
        ("conjugate-heat-transfer", "thermal-natural-convection-lbm") => Some(vec![kernel_probe(
            "crates/fs-lbm/src/thermal.rs",
            "pub struct ThermalLbm",
        )]),
        ("conjugate-heat-transfer", "forced-convection-correlations-and-fan-curve") => Some(vec![
            kernel_probe(
                "crates/fs-convection/src/lib.rs",
                "pub struct CorrelationInputs",
            ),
            kernel_probe("crates/fs-airflow/src/lib.rs", "pub struct FanCurve"),
        ]),
        ("conjugate-heat-transfer", "time-dependent-heat-adjoint") => Some(vec![
            kernel_probe("crates/fs-adjoint/src/timedep.rs", "pub struct HeatAdjoint"),
            kernel_probe(
                "crates/fs-adjoint/src/timedep.rs",
                "pub struct CoupledChtAdjoint",
            ),
        ]),
        ("conjugate-heat-transfer", "temperature-dependent-material-properties") => Some(vec![
            kernel_probe("crates/fs-matdb/src/lib.rs", "conductivity(T)"),
            kernel_probe(
                "crates/fs-matdb/src/lib.rs",
                "pub struct ElectronicsThermalMaterialSet",
            ),
        ]),
        ("conjugate-heat-transfer", "solid-fluid-thermal-coupling-and-contact-resistance") => {
            Some(vec![
                kernel_probe(
                    "crates/fs-conduction/src/interface.rs",
                    "pub struct InterfaceResistance",
                ),
                kernel_probe(
                    "crates/fs-couple/src/lib.rs",
                    "pub struct ThermalFieldTransfer",
                ),
            ])
        }
        ("aeroelastic-screening", "wing-shell-structure-and-modes") => Some(vec![
            kernel_probe(
                "crates/fs-solid/src/stability.rs",
                "pub struct BucklingResult",
            ),
            kernel_probe("crates/fs-solid/src/lib.rs", "pub struct ShellElement3d"),
        ]),
        ("aeroelastic-screening", "unsteady-aerodynamic-loads") => Some(vec![
            kernel_probe("crates/fs-vpm/src/lib.rs", "pub struct VortexParticle"),
            kernel_probe("crates/fs-vpm/src/lib.rs", "pub struct VortexFilament3d"),
        ]),
        ("aeroelastic-screening", "nonlinear-field-fsi") => Some(vec![
            kernel_probe(
                "crates/fs-couple/src/lib.rs",
                "pub struct AeroStructureFieldTransfer",
            ),
            kernel_probe(
                "crates/fs-couple/src/lib.rs",
                "pub struct NonlinearFsiSolver",
            ),
        ]),
        ("aeroelastic-screening", "coupled-flutter-gradient") => Some(vec![
            kernel_probe("crates/fs-flutter-e2e/src/lib.rs", "pub fn run_campaign"),
            kernel_probe(
                "crates/fs-flutter-e2e/src/lib.rs",
                "pub struct CoupledFlutterGradient",
            ),
        ]),
        ("additive-manufacturing-distortion", "moving-heat-source-and-phase-change") => Some(vec![
            kernel_probe(
                "crates/fs-conduction/src/lib.rs",
                "pub struct MovingHeatSource",
            ),
            kernel_probe(
                "crates/fs-conduction/src/lib.rs",
                "pub struct PhaseChangeModel",
            ),
        ]),
        ("additive-manufacturing-distortion", "three-dimensional-inelastic-distortion") => {
            Some(vec![
                kernel_probe(
                    "crates/fs-solid/src/lib.rs",
                    "pub struct J2ContinuumElement3d",
                ),
                kernel_probe(
                    "crates/fs-solid/src/lib.rs",
                    "pub struct AdditiveDistortionSolve",
                ),
            ])
        }
        ("additive-manufacturing-distortion", "layer-activation-time-sequencing") => Some(vec![
            kernel_probe("crates/fs-time/src/slabs.rs", "pub fn activation_report"),
            kernel_probe(
                "crates/fs-time/src/slabs.rs",
                "pub struct AmLayerActivation",
            ),
        ]),
        ("additive-manufacturing-distortion", "manufacturing-constraint-screen") => Some(vec![
            kernel_probe("crates/fs-fab/src/lib.rs", "pub fn overhang_angle"),
            kernel_probe(
                "crates/fs-fab/src/lib.rs",
                "pub struct GeometryOverhangEvaluator",
            ),
        ]),
        ("additive-manufacturing-distortion", "as-built-registration") => Some(vec![
            kernel_probe("crates/fs-asbuilt/src/lib.rs", "pub struct Point2"),
            kernel_probe("crates/fs-asbuilt/src/rigid3.rs", "pub struct Point3"),
            kernel_probe("crates/fs-asbuilt/src/lib.rs", "pub struct IcpSolver"),
        ]),
        _ => None,
    }
}

fn derive_kernel_readiness(root: &Path, probes: &[KernelProbe]) -> (Readiness, Vec<String>) {
    assert!(!probes.is_empty(), "kernel probe set must not be empty");
    let mut present = 0_usize;
    let mut evidence = Vec::with_capacity(probes.len());
    for probe in probes {
        let path = root.join(probe.path);
        let found =
            std::fs::read_to_string(&path).is_ok_and(|contents| contents.contains(probe.marker));
        if found {
            present += 1;
        }
        evidence.push(format!(
            "{}:{}={}",
            probe.path,
            probe.marker,
            if found { "found" } else { "missing" }
        ));
    }
    let readiness = if present == probes.len() {
        Readiness::Present
    } else if present == 0 {
        Readiness::Absent
    } else {
        Readiness::Partial
    };
    (readiness, evidence)
}

#[test]
fn kernel_readiness_matrix_matches_independent_workspace_probes() {
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    let root = manifest
        .parent()
        .and_then(Path::parent)
        .expect("fs-wedge lives at <workspace>/crates/fs-wedge");
    let mut failures = Vec::new();
    let mut checked = 0_usize;

    eprintln!("RESULT\tVERTICAL\tCAPABILITY\tMATRIX\tOBSERVED\tPROBES");
    for inputs in measured_wedge_inputs() {
        for entry in inputs.kernels {
            checked += 1;
            let Some(probes) = kernel_probes(inputs.vertical, entry.capability) else {
                failures.push(format!(
                    "{} {}: no independent probe specification",
                    inputs.vertical, entry.capability
                ));
                continue;
            };
            let (observed, evidence) = derive_kernel_readiness(root, &probes);
            let passed = observed == entry.measurement.readiness;
            eprintln!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                if passed { "PASS" } else { "FAIL" },
                inputs.vertical,
                entry.capability,
                entry.measurement.readiness.label(),
                observed.label(),
                evidence.join("; ")
            );
            if !passed {
                failures.push(format!(
                    "{} {}: matrix={} workspace={} [{}]",
                    inputs.vertical,
                    entry.capability,
                    entry.measurement.readiness.label(),
                    observed.label(),
                    evidence.join("; ")
                ));
            }
        }
    }
    assert_eq!(checked, 15, "kernel probe inventory must remain exhaustive");
    assert!(
        failures.is_empty(),
        "kernel readiness drift:\n{}",
        failures.join("\n")
    );
}

#[test]
fn current_workspace_evidence_paths_and_markers_have_not_drifted() {
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    let root = manifest
        .parent()
        .and_then(Path::parent)
        .expect("fs-wedge lives at <workspace>/crates/fs-wedge");
    let mut failures = Vec::new();

    eprintln!("RESULT\tVERTICAL\tAXIS\tENTRY\tPATH\tDETAIL");
    for inputs in measured_wedge_inputs() {
        for entry in inputs.kernels {
            failures.extend(check_workspace_measurement(
                root,
                inputs.vertical,
                InputAxis::KernelReadiness.label(),
                entry.capability,
                entry.measurement,
            ));
        }
        for entry in inputs.validation_data {
            failures.extend(check_workspace_measurement(
                root,
                inputs.vertical,
                InputAxis::ValidationDataAccess.label(),
                entry.dataset,
                entry.measurement,
            ));
        }
        for entry in inputs.cad_burden {
            failures.extend(check_workspace_measurement(
                root,
                inputs.vertical,
                InputAxis::CadBurden.label(),
                entry.required_geometry,
                entry.measurement,
            ));
        }
        for entry in inputs.compute_cost {
            failures.extend(check_workspace_measurement(
                root,
                inputs.vertical,
                InputAxis::ComputeCost.label(),
                entry.rung,
                entry.measurement,
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "evidence drift:\n{}",
        failures.join("\n")
    );
}

fn authenticated_snapshot(manifest: &str, bundle: &[u8]) -> HistoricalEvidenceSnapshot {
    let manifest_identity = Box::leak(
        hash_domain(
            COMPARISON_EVIDENCE_MANIFEST_IDENTITY_DOMAIN,
            manifest.as_bytes(),
        )
        .to_hex()
        .into_boxed_str(),
    );
    let bundle_identity = Box::leak(
        hash_domain(COMPARISON_EVIDENCE_BUNDLE_IDENTITY_DOMAIN, bundle)
            .to_hex()
            .into_boxed_str(),
    );
    HistoricalEvidenceSnapshot {
        manifest_identity_blake3: manifest_identity,
        bundle_identity_blake3: bundle_identity,
        ..COMPARISON_EVIDENCE_SNAPSHOT
    }
}

fn find_subslice(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack
        .windows(needle.len())
        .position(|window| window == needle)
}

fn test_tar_size(header: &[u8]) -> usize {
    let field = &header[124..136];
    let end = field
        .iter()
        .position(|byte| *byte == 0)
        .unwrap_or(field.len());
    usize::from_str_radix(
        std::str::from_utf8(&field[..end])
            .expect("test TAR size is ASCII")
            .trim(),
        8,
    )
    .expect("test TAR size is octal")
}

fn test_tar_header_offset(bundle: &[u8], path: &str) -> usize {
    let mut offset = 0usize;
    loop {
        let header = &bundle[offset..offset + 512];
        assert!(
            !header.iter().all(|byte| *byte == 0),
            "missing test TAR path {path}"
        );
        let end = header[..100]
            .iter()
            .position(|byte| *byte == 0)
            .unwrap_or(100);
        if &header[..end] == path.as_bytes() {
            return offset;
        }
        let size = test_tar_size(header);
        offset += 512 + size.div_ceil(512) * 512;
    }
}

fn rename_test_tar_file(bundle: &mut [u8], old: &str, new: &str) {
    assert_eq!(old.len(), new.len(), "test rename preserves header layout");
    let offset = test_tar_header_offset(bundle, old);
    let header = &mut bundle[offset..offset + 512];
    header[..new.len()].copy_from_slice(new.as_bytes());
    header[148..156].fill(b' ');
    let checksum: usize = header.iter().map(|byte| usize::from(*byte)).sum();
    let encoded = format!("{checksum:06o}\0 ");
    assert_eq!(encoded.len(), 8);
    header[148..156].copy_from_slice(encoded.as_bytes());
}

#[test]
fn historical_comparison_replays_from_retained_exact_source_bytes() {
    let first = verify_default_comparison_evidence().expect("historical snapshot replays");
    let second = verify_default_comparison_evidence().expect("historical replay is deterministic");
    assert_eq!(first, second);
    assert_eq!(
        first.inventory_revision(),
        COMPARISON_EVIDENCE_SNAPSHOT.inventory_revision
    );
    assert_eq!(first.source_count(), 13);
    assert_eq!(first.pointer_count(), 31);
    assert_eq!(
        first.manifest_identity_blake3(),
        COMPARISON_EVIDENCE_SNAPSHOT.manifest_identity_blake3
    );
    assert_eq!(
        first.bundle_identity_blake3(),
        COMPARISON_EVIDENCE_SNAPSHOT.bundle_identity_blake3
    );

    // The historical marker is deliberately absent from today's contract.
    // Replay succeeds because it reads the retained b3b bytes, never HEAD.
    let live_contract = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../fs-ladder/CONTRACT.md"),
    )
    .expect("live fs-ladder contract is readable");
    assert!(!live_contract.contains("The ladder does not run solves"));
    assert!(live_contract.contains("The `RANS` and `LES` rungs remain declarations"));
    assert!(
        comparison_evidence_manifest().contains(
            "full-electronics-cooling-cht\tlow-compute-cost\tcrates/fs-ladder/CONTRACT.md\tThe ladder does not run solves"
        )
    );
}

#[test]
fn historical_comparison_replay_fails_closed_on_artifact_and_revision_mutation() {
    let mut tampered_bundle = comparison_evidence_bundle().to_vec();
    tampered_bundle[1024] ^= 1;
    assert!(matches!(
        verify_comparison_evidence(
            COMPARISON_EVIDENCE_SNAPSHOT,
            comparison_evidence_manifest(),
            &tampered_bundle,
            comparison_candidates(),
        ),
        Err(HistoricalEvidenceError::IdentityMismatch {
            artifact: "bundle",
            ..
        })
    ));

    let tampered_manifest = comparison_evidence_manifest().replacen(
        "The ladder does not run solves",
        "The ladder does run solves",
        1,
    );
    assert!(matches!(
        verify_comparison_evidence(
            COMPARISON_EVIDENCE_SNAPSHOT,
            &tampered_manifest,
            comparison_evidence_bundle(),
            comparison_candidates(),
        ),
        Err(HistoricalEvidenceError::IdentityMismatch {
            artifact: "manifest",
            ..
        })
    ));

    let mut candidates = comparison_candidates().to_vec();
    candidates[0].inventory_revision = "0000000000000000000000000000000000000000";
    assert!(matches!(
        verify_comparison_evidence(
            COMPARISON_EVIDENCE_SNAPSHOT,
            comparison_evidence_manifest(),
            comparison_evidence_bundle(),
            &candidates,
        ),
        Err(HistoricalEvidenceError::CandidateRevisionMismatch {
            candidate: "full-electronics-cooling-cht",
            ..
        })
    ));
}

#[test]
fn historical_comparison_replay_refuses_missing_and_unavailable_sources() {
    let mut missing_source_bundle = comparison_evidence_bundle().to_vec();
    rename_test_tar_file(
        &mut missing_source_bundle,
        "crates/fs-adjoint/CONTRACT.md",
        "crates/fs-adjoinx/CONTRACT.md",
    );
    let missing_source_snapshot =
        authenticated_snapshot(comparison_evidence_manifest(), &missing_source_bundle);
    let error = verify_comparison_evidence(
        missing_source_snapshot,
        comparison_evidence_manifest(),
        &missing_source_bundle,
        comparison_candidates(),
    )
    .expect_err("renamed retained source must fail closed");
    match &error {
        HistoricalEvidenceError::SourceMismatch {
            inventory_revision,
            bundle_identity_blake3,
            reference,
            expected_bytes: Some(_),
            observed_bytes: None,
        } => {
            assert_eq!(
                *inventory_revision,
                COMPARISON_EVIDENCE_SNAPSHOT.inventory_revision
            );
            assert_eq!(
                *bundle_identity_blake3,
                missing_source_snapshot.bundle_identity_blake3
            );
            assert_eq!(reference, "crates/fs-adjoint/CONTRACT.md");
        }
        other => panic!("unexpected missing-source diagnostic: {other:?}"),
    }
    let diagnostic = error.to_string();
    assert!(diagnostic.contains(COMPARISON_EVIDENCE_SNAPSHOT.inventory_revision));
    assert!(diagnostic.contains("crates/fs-adjoint/CONTRACT.md"));

    let empty_snapshot = authenticated_snapshot(comparison_evidence_manifest(), &[]);
    let unavailable = verify_comparison_evidence(
        empty_snapshot,
        comparison_evidence_manifest(),
        &[],
        comparison_candidates(),
    )
    .expect_err("empty retained bundle must fail closed");
    assert!(matches!(
        &unavailable,
        HistoricalEvidenceError::MalformedBundle { offset: 0, .. }
    ));
    assert!(unavailable.to_string().contains("TAR length"));
}

#[test]
#[allow(clippy::too_many_lines)] // One two-stage mutation proves identity refusal precedes marker refusal.
fn historical_comparison_replay_refuses_reauthenticated_content_and_marker_mutation() {
    let mut content_mutation = comparison_evidence_bundle().to_vec();
    let non_marker = b"The fidelity-ladder registry";
    let content_offset =
        find_subslice(&content_mutation, non_marker).expect("non-marker source text is retained");
    content_mutation[content_offset] = b't';
    let content_snapshot =
        authenticated_snapshot(comparison_evidence_manifest(), &content_mutation);
    let content_error = verify_comparison_evidence(
        content_snapshot,
        comparison_evidence_manifest(),
        &content_mutation,
        comparison_candidates(),
    )
    .expect_err("per-source identity must reject reauthenticated content mutation");
    match &content_error {
        HistoricalEvidenceError::SourceIdentityMismatch {
            inventory_revision,
            candidate,
            factor,
            reference,
            locator,
            expected_blake3,
            observed_blake3,
        } => {
            assert_eq!(
                *inventory_revision,
                COMPARISON_EVIDENCE_SNAPSHOT.inventory_revision
            );
            assert_eq!(candidate.as_ref(), "full-electronics-cooling-cht");
            assert_eq!(factor.as_ref(), "low-compute-cost");
            assert_eq!(reference.as_ref(), "crates/fs-ladder/CONTRACT.md");
            assert_eq!(locator.as_ref(), "The ladder does not run solves");
            assert_ne!(expected_blake3, observed_blake3);
        }
        other => panic!("unexpected content-mutation diagnostic: {other:?}"),
    }
    let diagnostic = content_error.to_string();
    assert!(diagnostic.contains("full-electronics-cooling-cht/low-compute-cost"));
    assert!(diagnostic.contains("expected BLAKE3"));
    assert!(diagnostic.contains("observed"));

    let marker = b"The ladder does not run solves";
    let mut markerless_bundle = comparison_evidence_bundle().to_vec();
    let marker_offset =
        find_subslice(&markerless_bundle, marker).expect("historical marker is retained");
    markerless_bundle[marker_offset] = b't';
    let preliminary_snapshot =
        authenticated_snapshot(comparison_evidence_manifest(), &markerless_bundle);
    let observed_markerless_identity = match verify_comparison_evidence(
        preliminary_snapshot,
        comparison_evidence_manifest(),
        &markerless_bundle,
        comparison_candidates(),
    ) {
        Err(HistoricalEvidenceError::SourceIdentityMismatch {
            observed_blake3, ..
        }) => observed_blake3,
        other => panic!("expected source-identity refusal before marker replay: {other:?}"),
    };
    let original_ladder_identity = comparison_evidence_manifest()
        .lines()
        .find(|line| line.starts_with("SOURCE\tcrates/fs-ladder/CONTRACT.md\t"))
        .and_then(|line| line.split('\t').nth(5))
        .expect("ladder source identity is present");
    let markerless_manifest = comparison_evidence_manifest().replacen(
        original_ladder_identity,
        &observed_markerless_identity,
        1,
    );
    let markerless_snapshot = authenticated_snapshot(&markerless_manifest, &markerless_bundle);
    let marker_error = verify_comparison_evidence(
        markerless_snapshot,
        &markerless_manifest,
        &markerless_bundle,
        comparison_candidates(),
    )
    .expect_err("authenticated source without locator must fail closed");
    match &marker_error {
        HistoricalEvidenceError::MarkerMissing {
            inventory_revision,
            source_identity_blake3,
            candidate,
            factor,
            reference,
            locator,
        } => {
            assert_eq!(
                *inventory_revision,
                COMPARISON_EVIDENCE_SNAPSHOT.inventory_revision
            );
            assert_eq!(source_identity_blake3, &observed_markerless_identity);
            assert_eq!(candidate.as_ref(), "full-electronics-cooling-cht");
            assert_eq!(factor.as_ref(), "low-compute-cost");
            assert_eq!(reference.as_ref(), "crates/fs-ladder/CONTRACT.md");
            assert_eq!(locator.as_ref(), "The ladder does not run solves");
        }
        other => panic!("unexpected marker diagnostic: {other:?}"),
    }
    assert!(
        marker_error
            .to_string()
            .contains("after verifying source BLAKE3")
    );
}

#[test]
fn historical_comparison_replay_refuses_authenticated_pointer_path_mutation() {
    let canonical = "POINTER\tfull-electronics-cooling-cht\tlow-compute-cost\tcrates/fs-ladder/CONTRACT.md\tThe ladder does not run solves";
    let mutated = "POINTER\tfull-electronics-cooling-cht\tlow-compute-cost\tcrates/fs-laddeq/CONTRACT.md\tThe ladder does not run solves";
    let path_mutation = comparison_evidence_manifest().replacen(canonical, mutated, 1);
    let snapshot = authenticated_snapshot(&path_mutation, comparison_evidence_bundle());
    let error = verify_comparison_evidence(
        snapshot,
        &path_mutation,
        comparison_evidence_bundle(),
        comparison_candidates(),
    )
    .expect_err("authenticated pointer-path mutation must fail closed");
    match &error {
        HistoricalEvidenceError::PointerMismatch {
            inventory_revision,
            manifest_identity_blake3,
            expected,
            observed,
            ..
        } => {
            assert_eq!(
                *inventory_revision,
                COMPARISON_EVIDENCE_SNAPSHOT.inventory_revision
            );
            assert_eq!(*manifest_identity_blake3, snapshot.manifest_identity_blake3);
            assert!(expected.contains("crates/fs-ladder/CONTRACT.md"));
            assert!(observed.contains("crates/fs-laddeq/CONTRACT.md"));
        }
        other => panic!("unexpected pointer-path diagnostic: {other:?}"),
    }
    let diagnostic = error.to_string();
    assert!(diagnostic.contains(COMPARISON_EVIDENCE_SNAPSHOT.inventory_revision));
    assert!(diagnostic.contains("manifest BLAKE3"));
}

#[test]
fn historical_pointer_sequence_comparison_is_structural() {
    let prefix = "POINTER\tfull-electronics-cooling-cht\t";
    let collision_manifest = comparison_evidence_manifest()
        .lines()
        .map(|line| {
            line.strip_prefix(prefix)
                .map_or_else(|| line.to_string(), |rest| format!("{prefix}foo/{rest}"))
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let snapshot = authenticated_snapshot(&collision_manifest, comparison_evidence_bundle());
    let mut candidates = comparison_candidates().to_vec();
    candidates[0].name = "full-electronics-cooling-cht/foo";
    let error = verify_comparison_evidence(
        snapshot,
        &collision_manifest,
        comparison_evidence_bundle(),
        &candidates,
    )
    .expect_err("delimiter-colliding pointer tuples must remain unequal");
    assert!(matches!(
        error,
        HistoricalEvidenceError::PointerMismatch { index: 0, .. }
    ));
}

#[test]
fn explicit_comparison_is_evidence_complete_and_ranked() {
    let candidates = comparison_candidates();
    assert_eq!(candidates.len(), 3);
    for candidate in candidates {
        assert_eq!(candidate.factors.len(), ScoringFactor::ALL.len());
        for factor in ScoringFactor::ALL {
            let input = candidate
                .factors
                .iter()
                .find(|input| input.factor == factor)
                .expect("every comparison factor is present");
            assert!(input.is_complete(), "{candidate:?} {input:?}");
            assert!(!input.measurement.evidence.is_empty());
        }
    }

    let record = default_recommendation().expect("default comparison is valid");
    assert_eq!(record.recommended, "thermal-design-assurance");
    assert_eq!(record.runner_up, "sdf-structural-topology-assurance");
    assert!(record.minority_report.contains("lowest technical-risk"));
    assert_eq!(record.ranked[0].weighted_total, 638);
    assert_eq!(record.ranked[1].weighted_total, 623);
    assert_eq!(record.ranked[2].weighted_total, 502);
}

#[test]
fn scoring_refuses_bad_weights_and_is_weight_order_invariant() {
    let baseline = score_candidates(&DEFAULT_FACTOR_WEIGHTS, comparison_candidates())
        .expect("default scoring succeeds");
    let mut reversed = DEFAULT_FACTOR_WEIGHTS;
    reversed.reverse();
    assert_eq!(
        baseline,
        score_candidates(&reversed, comparison_candidates()).expect("reordered weights succeed")
    );

    let mut wrong_sum = DEFAULT_FACTOR_WEIGHTS;
    wrong_sum[0].weight += 1;
    assert_eq!(
        score_candidates(&wrong_sum, comparison_candidates()),
        Err(ScoringError::WeightsNotNormalized { sum: 101 })
    );

    let mut duplicate = DEFAULT_FACTOR_WEIGHTS;
    duplicate[0].factor = duplicate[1].factor;
    assert_eq!(
        score_candidates(&duplicate, comparison_candidates()),
        Err(ScoringError::DuplicateWeight {
            factor: ScoringFactor::KernelReadiness
        })
    );
}

#[test]
fn every_factor_is_monotone_under_a_positive_weight() {
    let source = comparison_candidates()[1];
    let baseline = score_candidates(&DEFAULT_FACTOR_WEIGHTS, &[source])
        .expect("one-candidate score succeeds")[0]
        .weighted_total;
    for factor in ScoringFactor::ALL {
        let mut factors: [fs_wedge::FactorRating; 9] = source
            .factors
            .try_into()
            .expect("comparison has exactly nine factors");
        let input = factors
            .iter_mut()
            .find(|input| input.factor == factor)
            .expect("factor exists");
        input.rating += 1;
        let improved = ComparisonCandidate {
            factors: Box::leak(Box::new(factors)),
            ..source
        };
        let improved_total = score_candidates(&DEFAULT_FACTOR_WEIGHTS, &[improved])
            .expect("improved candidate remains valid")[0]
            .weighted_total;
        assert!(
            improved_total > baseline,
            "{} was not monotone",
            factor.label()
        );
    }
}

#[test]
fn candidate_permutation_and_tie_breaking_are_deterministic() {
    let candidates = comparison_candidates();
    let baseline = score_candidates(&DEFAULT_FACTOR_WEIGHTS, candidates).unwrap();
    let permuted = [candidates[2], candidates[0], candidates[1]];
    assert_eq!(
        baseline,
        score_candidates(&DEFAULT_FACTOR_WEIGHTS, &permuted).unwrap()
    );

    let alpha = ComparisonCandidate {
        name: "alpha",
        display: "Alpha",
        ..candidates[0]
    };
    let beta = ComparisonCandidate {
        name: "beta",
        display: "Beta",
        ..candidates[0]
    };
    let tied = score_candidates(&DEFAULT_FACTOR_WEIGHTS, &[beta, alpha]).unwrap();
    assert_eq!(tied[0].candidate, "alpha");
    assert_eq!(tied[0].weighted_total, tied[1].weighted_total);
}

#[test]
fn sensitivity_tables_expose_flips_and_degenerate_ties() {
    let record = default_recommendation().unwrap();
    let expected = 2 * ScoringFactor::ALL.len();
    assert_eq!(record.rating_sensitivities.len(), expected);
    assert_eq!(record.weight_sensitivities.len(), expected);
    assert!(record.rating_sensitivities.iter().any(|row| {
        row.challenger == record.runner_up
            && row.factor == ScoringFactor::KernelReadiness
            && row.required_rating.is_some()
    }));
    assert!(record.weight_sensitivities.iter().any(|row| {
        row.challenger == record.runner_up
            && row.factor == ScoringFactor::KernelReadiness
            && row.required_weight.is_some()
    }));
    assert!(record.rating_sensitivities.iter().all(|row| {
        row.challenger != "full-electronics-cooling-cht" || row.required_rating.is_none()
    }));
    for row in record
        .weight_sensitivities
        .iter()
        .filter(|row| row.challenger == "full-electronics-cooling-cht")
    {
        let ties_thermal_at_full_weight = matches!(
            row.factor,
            ScoringFactor::CustomerPain | ScoringFactor::DataAccess | ScoringFactor::RegulatoryRisk
        );
        assert_eq!(
            row.required_weight,
            ties_thermal_at_full_weight.then_some(100),
            "unexpected full-CHT weight sensitivity for {}",
            row.factor.label()
        );
    }
}

#[test]
fn verbose_comparison_report_is_deterministic() {
    let first = render_comparison_report().expect("comparison report renders");
    let second = render_comparison_report().expect("comparison report replays");
    assert_eq!(first, second);
    assert!(first.starts_with("FS-WEDGE-COMPARISON\tv2\n"));
    assert_eq!(first.matches("HISTORICAL_EVIDENCE\t").count(), 1);
    assert_eq!(first.matches("FACTOR\t").count(), 27);
    assert_eq!(first.matches("RATING_FLIP\t").count(), 18);
    assert_eq!(first.matches("WEIGHT_FLIP\t").count(), 18);
    assert!(first.contains("RECOMMENDED\tthermal-design-assurance"));
    assert!(first.contains("MINORITY_REPORT\tSDF structural assurance"));
    eprintln!("{first}");
}

#[test]
fn three_verticals_are_ranked_with_proposal_mappings() {
    let vs = verticals();
    assert_eq!(vs.len(), 3);
    let mut ranks: Vec<u8> = vs.iter().map(|v| v.rank).collect();
    ranks.sort_unstable();
    assert_eq!(ranks, vec![1, 2, 3]);
    // second vertical exercises Proposal 1; third exercises 11 and 4.
    let aero = vs
        .iter()
        .find(|v| v.name == "aeroelastic-screening")
        .unwrap();
    assert_eq!(aero.rank, 2);
    assert!(aero.exercises.contains(&"1"));
    let am = vs
        .iter()
        .find(|v| v.name == "additive-manufacturing-distortion")
        .unwrap();
    assert_eq!(am.rank, 3);
    assert!(am.exercises.contains(&"11") && am.exercises.contains(&"4"));
    // every vertical names at least one exercised proposal.
    assert!(vs.iter().all(|v| !v.exercises.is_empty()));
}

#[test]
fn the_cycle_time_kill_criterion_is_measurable_and_conservative() {
    assert!((CHT_BASELINE.target_reduction - 3.0).abs() < 1e-12);
    assert_eq!(CHT_BASELINE.kill_within_quarters, 2);

    // Met: even the LOW bound (0.5625 days) clears 3x against 0.15 days.
    let met = CHT_BASELINE
        .evaluate_kill_criterion(0.15)
        .expect("measurable");
    assert_eq!(met.verdict, KillVerdict::Met);
    assert!(met.reduction_low >= 3.0);

    // NotMet: even the HIGH bound (16 days) cannot clear 3x against 6 days.
    let not_met = CHT_BASELINE
        .evaluate_kill_criterion(6.0)
        .expect("measurable");
    assert_eq!(not_met.verdict, KillVerdict::NotMet);
    assert!(not_met.reduction_high < 3.0);

    // Indeterminate: the envelope straddles the target against 2 days, so
    // marketing may not claim "met".
    let straddle = CHT_BASELINE
        .evaluate_kill_criterion(2.0)
        .expect("measurable");
    assert_eq!(straddle.verdict, KillVerdict::Indeterminate);
    assert!(straddle.reduction_low < 3.0 && straddle.reduction_high >= 3.0);

    // Non-measurable cycle times are refused, not coerced.
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let refusal = CHT_BASELINE.evaluate_kill_criterion(bad);
        assert!(
            matches!(
                refusal,
                Err(KillCriterionError::NonMeasurableCycleTime { .. })
            ),
            "accepted non-measurable cycle time {bad}"
        );
    }
}

#[test]
fn the_audit_is_complete() {
    let a = audit();
    assert!(a.ok(), "gaps: {:?}", a.gaps);
    assert!(a.passed("historic-scores-superseded"));
    assert!(a.passed("measured-inputs-complete"));
    assert!(a.passed("no-absent-strong-scores"));
    assert!(a.passed("comparison-inputs-complete"));
    assert!(a.passed("default-weights-normalized"));
    assert!(a.passed("comparison-history-bound"));
    assert!(a.passed("comparison-ranking-complete"));
    assert!(a.passed("comparison-sensitivity-complete"));
    assert!(a.passed("ranks-complete"));
    assert!(a.passed("all-exercise-proposals"));
    assert!(a.passed("kill-criterion-measurable"));
    assert!(a.passed("cycle-time-baseline-measured"));
    assert!(a.passed("placeholder-baseline-refused"));
    assert_eq!(a.checks.len(), 13);
}

#[test]
fn the_negative_doctrine_is_stated() {
    // the load-bearing anti-pattern: don't sell against peak single-physics.
    assert!(
        WEDGE_DOCTRINE
            .to_lowercase()
            .contains("peak single-physics")
    );
    // criterion labels are unique.
    let labels: Vec<&str> = WedgeCriterion::ALL.iter().map(|c| c.label()).collect();
    let mut sorted = labels.clone();
    sorted.sort_unstable();
    sorted.dedup();
    assert_eq!(sorted.len(), labels.len());
}

#[test]
fn json_is_well_formed_and_deterministic() {
    let j = to_json();
    assert_eq!(j, to_json());
    assert!(j.starts_with('{') && j.ends_with('}'));
    assert!(j.contains("conjugate-heat-transfer"));
    assert!(j.contains("\"score_use\":\"superseded-for-decision-use\""));
    assert!(j.contains("\"measured_inputs\":"));
    assert!(j.contains("\"validation_data\":"));
    assert!(j.contains("NIST Additive Manufacturing Benchmark Test Series"));
    assert!(j.contains("\"target_reduction\":3"));
    assert_eq!(j.matches("\"rank\":").count(), 3);
}

#[test]
fn the_measured_baseline_is_complete_and_published_source_derived() {
    let baseline = CHT_BASELINE;
    assert!(
        baseline.is_complete(),
        "measured baseline record incomplete"
    );
    assert_eq!(
        baseline.provenance,
        BaselineProvenance::PublishedSourceDerived
    );
    assert_eq!(baseline.steps.len(), IncumbentStep::ALL.len());
    for (estimate, step) in baseline.steps.iter().zip(IncumbentStep::ALL) {
        assert_eq!(estimate.step, step, "steps out of pipeline order");
        assert!(
            estimate
                .sources
                .iter()
                .any(|source| source.class.is_load_bearing()),
            "step {} rests on vendor material alone",
            estimate.step.label()
        );
        for source in estimate.sources {
            assert!(
                source.is_complete(),
                "incomplete source in {}",
                estimate.step.label()
            );
            assert!(
                source.url.starts_with("https://"),
                "source URL is not https in {}",
                estimate.step.label()
            );
        }
    }
}

#[test]
fn the_baseline_dossiers_exist_and_carry_their_markers() {
    // Both the live record and the retired placeholder point outside the Rust
    // declaration that embeds the pointer, so neither marker can satisfy its
    // own `contains` check.
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let repo_root = Path::new(manifest_dir)
        .ancestors()
        .nth(2)
        .expect("workspace root above crates/fs-wedge");
    for pointer in [CHT_BASELINE.dossier, RETIRED_PLACEHOLDER_BASELINE.dossier] {
        assert_ne!(
            pointer.reference, "crates/fs-wedge/src/lib.rs",
            "baseline dossier pointers must not self-verify against their declaring source"
        );
        assert_eq!(
            pointer.reference, "crates/fs-wedge/data/cycle-time-baseline-dossier.md",
            "baseline records must resolve through the retained external dossier"
        );
        let path = repo_root.join(pointer.reference);
        let contents = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("unreadable dossier {}: {error}", path.display()));
        assert!(
            contents.contains(pointer.locator),
            "dossier {} lost its marker {}",
            path.display(),
            pointer.locator
        );
    }

    let dossier_path = repo_root.join(CHT_BASELINE.dossier.reference);
    let dossier = std::fs::read_to_string(&dossier_path)
        .unwrap_or_else(|error| panic!("unreadable dossier {}: {error}", dossier_path.display()));
    // Every load-bearing URL in the record must appear in the dossier, so the
    // Rust constants cannot drift away from the provenance document.
    for step in CHT_BASELINE.steps {
        for source in step.sources {
            assert!(
                dossier.contains(source.url),
                "dossier does not cite {} (step {})",
                source.url,
                step.step.label()
            );
        }
    }
}

#[test]
fn the_baseline_envelope_matches_the_dossier_derivation() {
    // Dossier section 3: low 4.5 h, high 128 h; 8 h/working day.
    let baseline = CHT_BASELINE;
    assert!((baseline.low_hours() - 4.5).abs() < 1e-12);
    assert!((baseline.high_hours() - 128.0).abs() < 1e-12);
    assert!((baseline.baseline_days_low() - 0.5625).abs() < 1e-12);
    assert!((baseline.baseline_days_high() - 16.0).abs() < 1e-12);
    // The retired placeholder's flat figure lies inside the envelope; the
    // range, not the point, is the citable object.
    assert!(baseline.baseline_days_low() < 5.0 && 5.0 < baseline.baseline_days_high());
}

#[test]
fn the_placeholder_baseline_is_refused_with_a_typed_error() {
    let refusal = RETIRED_PLACEHOLDER_BASELINE.evaluate_kill_criterion(1.0);
    assert_eq!(
        refusal,
        Err(KillCriterionError::PlaceholderBaseline {
            vertical: "conjugate-heat-transfer",
        })
    );
}

#[test]
fn the_kill_evaluation_prints_its_full_derivation() {
    // The e2e auditability requirement: recompute the criterion from the
    // stored record and log the complete derivation.
    let evaluation = CHT_BASELINE
        .evaluate_kill_criterion(2.0)
        .expect("measurable");
    let derivation = &evaluation.derivation;
    println!("{derivation}");
    // Record identity: which baseline record was used.
    assert!(derivation.contains("conjugate-heat-transfer"));
    assert!(derivation.contains("2026-07-22"));
    assert!(derivation.contains("published-source-derived"));
    assert!(derivation.contains("cycle-time-baseline-dossier.md"));
    // Baseline breakdown: every step appears with its bounds.
    for step in IncumbentStep::ALL {
        assert!(
            derivation.contains(step.label()),
            "derivation omits step {}",
            step.label()
        );
    }
    // Envelope, conversion, target, measured value, verdict.
    assert!(derivation.contains("0.56..16.00 working days"));
    assert!(derivation.contains("8 h/day"));
    assert!(derivation.contains("target 3.0x"));
    assert!(derivation.contains("2.000 days"));
    assert!(derivation.contains("verdict: indeterminate"));
}

#[test]
fn the_manifest_round_trips_the_baseline_envelope() {
    // The manifest renders the measured envelope; parse the numbers back out
    // and compare them with the record so the two spellings cannot drift.
    let json = to_json();
    assert!(json.contains("\"cycle_time_baseline\":{"));
    assert!(json.contains("\"provenance\":\"published-source-derived\""));
    let extract = |key: &str| -> f64 {
        let marker = format!("\"{key}\":");
        let start = json
            .find(&marker)
            .unwrap_or_else(|| panic!("missing {key}"))
            + marker.len();
        let rest = &json[start..];
        let end = rest
            .find([',', '}'])
            .unwrap_or_else(|| panic!("unterminated {key}"));
        rest[..end]
            .parse::<f64>()
            .unwrap_or_else(|error| panic!("unparsable {key}: {error}"))
    };
    assert_eq!(
        extract("baseline_days_low").to_bits(),
        CHT_BASELINE.baseline_days_low().to_bits()
    );
    assert_eq!(
        extract("baseline_days_high").to_bits(),
        CHT_BASELINE.baseline_days_high().to_bits()
    );
    assert_eq!(
        extract("target_reduction").to_bits(),
        CHT_BASELINE.target_reduction.to_bits()
    );
    assert_eq!(
        extract("hours_per_working_day").to_bits(),
        CHT_BASELINE.hours_per_working_day.to_bits()
    );
}

#[test]
fn the_audit_covers_the_measured_baseline_and_the_refusal_drill() {
    let report = audit();
    assert!(report.passed("cycle-time-baseline-measured"));
    assert!(report.passed("placeholder-baseline-refused"));
    assert!(report.passed("kill-criterion-measurable"));
}
