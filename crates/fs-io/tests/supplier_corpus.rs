//! G0/G2/G5 retained real-supplier corpus coverage
//! (`frankensim-extreal-program-f85xj.11.6`).

use fs_blake3::hash_bytes;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_io::{
    CorpusSourceKind, SupplierCadFormat, SupplierCorpusPolicy, parse_supplier_corpus_manifest,
    run_supplier_corpus,
};
use std::path::{Path, PathBuf};

fn corpus_root() -> PathBuf {
    if let Some(root) = std::env::var_os("FRANKENSIM_SUPPLIER_CORPUS_ROOT") {
        return PathBuf::from(root);
    }
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/cad-import-corpus")
}

fn read_manifest() -> fs_io::CorpusManifest {
    let source = std::fs::read_to_string(corpus_root().join("corpus-v1.tsv"))
        .expect("retained manifest must be readable");
    parse_supplier_corpus_manifest(&source).expect("retained manifest must satisfy the v1 schema")
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x11_06_2026,
                kernel_id: 11,
                tile: 6,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

#[test]
fn g0_retained_sources_match_manifest_byte_identities() {
    let manifest = read_manifest();
    let mut mismatches = Vec::new();
    for case in &manifest.cases {
        let bytes = std::fs::read(corpus_root().join(&case.relative_path))
            .unwrap_or_else(|error| panic!("{} must be retained: {error}", case.case_id));
        let actual = hash_bytes(&bytes);
        if actual != case.content_blake3 {
            mismatches.push(format!(
                "{}\tmanifest={}\tactual={}",
                case.case_id, case.content_blake3, actual
            ));
        }
    }
    assert!(
        mismatches.is_empty(),
        "retained byte identities disagree:\n{}",
        mismatches.join("\n")
    );
}

#[test]
fn g2_g5_real_supplier_sweep_is_population_complete_and_deterministic() {
    let manifest = read_manifest();
    assert_eq!(manifest.cases.len(), 21);
    assert!(manifest.meets_minimum_population());
    assert_eq!(
        manifest
            .cases
            .iter()
            .filter(|case| case.format == SupplierCadFormat::Step)
            .count(),
        10
    );
    assert!(manifest.covers_required_formats());
    assert!(manifest.covers_required_quality_tiers());
    assert_eq!(
        manifest
            .cases
            .iter()
            .filter(|case| case.format == SupplierCadFormat::Stl)
            .count(),
        10
    );
    assert_eq!(
        manifest
            .cases
            .iter()
            .filter(|case| case.format == SupplierCadFormat::Ply)
            .count(),
        1
    );
    assert_eq!(
        manifest
            .cases
            .iter()
            .filter(|case| case.source_kind == CorpusSourceKind::HttpSnapshot)
            .count(),
        1
    );
    assert!(
        !manifest.annotations_locked(),
        "agent-proposed baselines must not impersonate human review"
    );

    let policy = SupplierCorpusPolicy::try_standing_lane().expect("standing policy is valid");
    let run = || {
        with_cx(|cx| {
            run_supplier_corpus(&manifest, policy, cx, |case| {
                std::fs::read(corpus_root().join(&case.relative_path))
                    .map_err(|error| error.to_string())
            })
        })
        .expect("real corpus sweep must not be cancelled")
    };
    let first = run();
    let second = run();
    println!("supplier_corpus_scorecard={}", first.to_json());
    assert_eq!(first, second);
    assert_eq!(first.artifact_identity(), second.artifact_identity());
    assert_eq!(first.to_json(), second.to_json());
    assert_eq!(first.rows.len(), 21);
    assert_eq!(first.unreviewed, 21);
    assert_eq!(first.mismatches, 0);
    assert_eq!(first.proposed_mismatches, 0);
    assert_eq!(first.clean, 3);
    assert_eq!(first.repaired, 0);
    assert_eq!(first.refused, 18);
    let metrics = first.import_metrics();
    assert_eq!(metrics.total(), 21);
    assert_eq!(metrics.reviewed(), 0);
    assert_eq!(metrics.clean(), 0);
    assert_eq!(metrics.repaired(), 0);
    assert_eq!(metrics.refused(), 0);
    assert_eq!(metrics.annotation_mismatches(), 0);
    assert_eq!(first.summary_json(), second.summary_json());
    let tracked_summary = std::fs::read_to_string(corpus_root().join("scorecard-summary-v1.json"))
        .expect("tracked compact scorecard projection must be readable");
    assert_eq!(
        tracked_summary.trim_end_matches(['\r', '\n']),
        first.summary_json(),
        "tracked compact scorecard projection must equal the real full sweep"
    );
    assert!(first.summary_json().contains(
        "\"reviewed\":{\"total\":0,\"clean\":0,\"repaired\":0,\"refused\":0,\
         \"annotation_mismatch\":0}"
    ));
    assert!(first.covers_required_formats());
    assert!(first.covers_required_quality_tiers());
    assert!(!first.lane_passes());
    assert!(
        first
            .to_json()
            .contains("\"authority\":\"retained-population-metric-not-universal-import-rate\"")
    );
}
