//! Conformance battery for the acoustic corpus catalog (music bead
//! `frankensim-music-v8-root-3ez8g.1.1`): structural invariants, licensing
//! completeness, and manifest coherence — the tracked TSV is a projection
//! of the code catalog and must byte-equal its rendering, so neither can
//! drift from the other.

use fs_vvreg::acoustic::{
    ACOUSTIC_MANIFEST_LOCATOR, AcousticAcceptance, AcousticLevel, acoustic_absences,
    acoustic_cases, render_acoustic_manifest,
};
use std::collections::BTreeSet;

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("workspace root")
        .to_path_buf()
}

#[test]
fn catalog_ids_are_unique_and_prefixed() {
    let mut seen = BTreeSet::new();
    for case in acoustic_cases() {
        assert!(
            case.id.starts_with("acoustic-"),
            "case id {} must carry the acoustic- prefix",
            case.id
        );
        assert!(seen.insert(case.id), "duplicate case id {}", case.id);
    }
    for absence in acoustic_absences() {
        assert!(
            absence.id.starts_with("acoustic-"),
            "absence id {} must carry the acoustic- prefix",
            absence.id
        );
        assert!(seen.insert(absence.id), "duplicate id {}", absence.id);
    }
    assert_eq!(
        seen.len(),
        acoustic_cases().len() + acoustic_absences().len()
    );
}

#[test]
fn every_case_is_structurally_complete() {
    for case in acoustic_cases() {
        assert!(
            case.reference_value_si.is_finite() && case.reference_value_si > 0.0,
            "{}: reference must be finite and positive",
            case.id
        );
        for field in [
            case.title,
            case.metric,
            case.formula,
            case.source,
            case.license,
            case.retention,
            case.no_claim_reason,
        ] {
            assert!(!field.trim().is_empty(), "{}: empty field", case.id);
            assert!(
                !field.contains('\t') && !field.contains('\n'),
                "{}: manifest fields must be single-line, tab-free",
                case.id
            );
        }
        match case.acceptance {
            AcousticAcceptance::CentsEnvelope { outer, inner } => {
                assert!(
                    outer > inner && inner > 0.0,
                    "{}: cents envelope must nest (outer > inner > 0)",
                    case.id
                );
            }
            AcousticAcceptance::Relative { rtol } => {
                assert!(
                    rtol > 0.0 && rtol < 1.0,
                    "{}: rtol must be in (0, 1)",
                    case.id
                );
            }
            AcousticAcceptance::AbsoluteTolerance { atol } => {
                assert!(
                    atol.is_finite() && atol > 0.0,
                    "{}: atol must be finite and positive",
                    case.id
                );
            }
            AcousticAcceptance::ClassBand { lo_hz, hi_hz } => {
                assert!(
                    lo_hz.is_finite() && hi_hz.is_finite() && hi_hz > lo_hz,
                    "{}: class band must be finite and ordered (lo < hi)",
                    case.id
                );
                assert!(
                    case.reference_value_si >= lo_hz && case.reference_value_si <= hi_hz,
                    "{}: reference value must land inside its own class band",
                    case.id
                );
            }
        }
        for context in case.context {
            assert!(
                context.lo.is_finite() && context.hi >= context.lo,
                "{}: context {} must be a finite ordered range",
                case.id,
                context.name
            );
        }
        // Licensing law: published-experiment rows must state their license
        // decision and retention explicitly; analytic rows must say why no
        // transcription trust is needed.
        match case.level {
            AcousticLevel::PublishedExperiment | AcousticLevel::PublishedModelData => {
                assert!(
                    case.license.contains("CC-BY")
                        || case.license.contains("two-source")
                        || case.license.contains("government"),
                    "{}: published-experiment rows need a named license basis",
                    case.id
                );
            }
            AcousticLevel::AnalyticDefinition => {
                assert!(
                    case.license.contains("analytic")
                        || case.license.contains("factual")
                        || case.license.contains("government"),
                    "{}: analytic rows state why transcription trust is not needed",
                    case.id
                );
            }
        }
    }
}

#[test]
fn absences_carry_their_hunt_rationale() {
    for absence in acoustic_absences() {
        for field in [absence.title, absence.what, absence.why, absence.unblocks] {
            assert!(!field.trim().is_empty(), "{}: empty field", absence.id);
            assert!(
                !field.contains('\t') && !field.contains('\n'),
                "{}: fields must be single-line, tab-free",
                absence.id
            );
        }
    }
    // The openwind refusal is a deliberate decision, not a gap; it must
    // exist and say so (the licensing DONE-WHEN of bead 3ez8g.1.1).
    let refusal = acoustic_absences()
        .iter()
        .find(|absence| absence.id == "acoustic-refused-openwind-curves")
        .expect("the openwind retention refusal must stay recorded");
    assert!(refusal.why.contains("GPLv3"));
}

#[test]
fn free_free_bar_ratios_are_rederivable() {
    // Self-verified pin: recompute the cosh*cos = 1 roots by bisection and
    // confirm the catalog ratios — the analytic rows never rest on
    // transcription trust (music-plan doctrine: boundary conditions do not
    // travel with function names, so the characteristic equation itself is
    // the authority).
    fn root(lo: f64, hi: f64) -> f64 {
        let f = |x: f64| x.cosh() * x.cos() - 1.0;
        let (mut lo, mut hi) = (lo, hi);
        assert!(f(lo) * f(hi) < 0.0, "root must be bracketed");
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if f(lo) * f(mid) <= 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        0.5 * (lo + hi)
    }
    let beta1 = root(4.6, 4.8);
    let beta2 = root(7.8, 7.9);
    let beta3 = root(10.9, 11.1);
    let f2f1 = (beta2 / beta1).powi(2);
    let f3f1 = (beta3 / beta1).powi(2);
    let cases = acoustic_cases();
    let cat_f2f1 = cases
        .iter()
        .find(|case| case.id == "acoustic-bar-free-free-f2f1")
        .expect("f2f1 row");
    let cat_f3f1 = cases
        .iter()
        .find(|case| case.id == "acoustic-bar-free-free-f3f1")
        .expect("f3f1 row");
    assert!(
        (f2f1 - cat_f2f1.reference_value_si).abs() < 1.0e-9,
        "f2/f1 rederived {f2f1} vs catalog {}",
        cat_f2f1.reference_value_si
    );
    assert!(
        (f3f1 - cat_f3f1.reference_value_si).abs() < 1.0e-9,
        "f3/f1 rederived {f3f1} vs catalog {}",
        cat_f3f1.reference_value_si
    );
    // Falsifier direction: the PINNED tensioned-string family would give
    // f2/f1 = 4 exactly; the free-free ratio must be far from it, so the
    // gate provably discriminates boundary conditions.
    assert!((cat_f2f1.reference_value_si - 4.0).abs() > 1.0);
}

#[test]
fn tracked_manifest_matches_the_catalog_rendering() {
    let rendered = render_acoustic_manifest();
    let path = repo_root().join(ACOUSTIC_MANIFEST_LOCATOR);
    let tracked = std::fs::read_to_string(&path).unwrap_or_else(|error| {
        panic!(
            "{ACOUSTIC_MANIFEST_LOCATOR} unreadable ({error}); if this is the first run, \
             write the rendering below to that path:\n{rendered}"
        )
    });
    assert_eq!(
        tracked, rendered,
        "tracked acoustic manifest differs from the catalog rendering; regenerate the TSV \
         from render_acoustic_manifest() in the same commit as the catalog change"
    );
    // Sanity on the projection itself: every id appears exactly once.
    for case in acoustic_cases() {
        assert_eq!(
            tracked.matches(case.id).count(),
            1,
            "{} must appear exactly once in the manifest",
            case.id
        );
    }
}
