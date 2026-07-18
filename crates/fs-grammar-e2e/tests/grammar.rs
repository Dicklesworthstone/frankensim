//! End-to-end battery: a diverse fabricable program family is illuminated, the
//! best matches the target, and every simplification's certificate is re-verified.

use fs_evidence::Color;
use fs_grammar_e2e::{
    SimplificationCheckStatus, SimplificationSummary, assess_simplification, build_program,
    run_campaign, target,
};
use fs_shapeprog::{Geom, SimplifyRefusal, max_sdf_discrepancy};

#[test]
fn a_fabricable_program_family_is_illuminated_and_simplified_soundly() {
    let report = run_campaign(0.2, 0.03);
    // ILLUMINATION: a diverse archive, not one model.
    assert!(
        report.num_elites >= 5,
        "too few niches: {}",
        report.num_elites
    );
    assert!(report.coverage > 0.0 && report.qd_score > 0.0); // fitness = 1/(1+discrepancy)
    // the best program genuinely matches the peanut target.
    assert!(
        report.best_discrepancy < 0.2,
        "best discrepancy {}",
        report.best_discrepancy
    );
    // CERTIFICATE-PRESERVING SIMPLIFICATION: some programs shrank, and EVERY
    // simplification's re-measured error stayed within its certified bound.
    assert!(
        report.simplification.simplified_count() > 0,
        "nothing simplified"
    );
    assert!(
        report.simplification.size_after() < report.simplification.size_before(),
        "no size reduction"
    );
    assert!(
        report.simplification.is_sound(),
        "a rewrite certificate was refused or its independent check failed: {:#?}",
        report.simplification
    );
    assert!(
        report
            .simplification
            .is_complete_and_sound(report.num_elites),
        "every elite must have one sound assessment"
    );
    assert_eq!(
        report.simplification.radius_threshold().to_bits(),
        0.03_f64.to_bits()
    );
    assert_eq!(
        report.simplification.max_certified_error().to_bits(),
        0.04_f64.to_bits(),
        "the 0.02 offset admitted by a 0.03 radius threshold has a 0.04 context-free certificate"
    );
    assert!(
        report.simplification.max_certified_error() > report.simplification.radius_threshold(),
        "the local admission threshold must not be mislabeled as a global error budget"
    );
    assert!(
        report.simplification.max_sampled_discrepancy()
            <= report.simplification.max_certified_error(),
        "conservative outward sampled check exceeds the compositional certificate: {:#?}",
        report.simplification
    );
    assert_eq!(
        report.simplification.status(),
        SimplificationCheckStatus::Certified
    );
    assert_eq!(report.simplification.simplifier_refusals(), 0);
    assert_eq!(report.simplification.non_finite_certificates(), 0);
    assert_eq!(report.simplification.negative_certificates(), 0);
    assert_eq!(report.simplification.discrepancy_evidence_refusals(), 0);
    assert_eq!(report.simplification.structural_empty_agreements(), 0);
    assert_eq!(report.simplification.certificate_check_exceedances(), 0);
    assert_eq!(report.simplification.threshold_mismatches(), 0);
    assert_eq!(
        report.simplified_count,
        report.simplification.simplified_count()
    );
    assert_eq!(report.size_before, report.simplification.size_before());
    assert_eq!(report.size_after, report.simplification.size_after());
    assert_eq!(
        report.max_certified_error.to_bits(),
        report.simplification.max_certified_error().to_bits()
    );
    assert_eq!(
        report.simplification_sound,
        report
            .simplification
            .is_complete_and_sound(report.num_elites)
    );
    // FABRICABILITY: the constraint genuinely discriminates — some elites pass
    // the minimum feature size and some (the thin-sphere ones) do not.
    assert!(report.fab_satisfied > 0);
    assert!(
        report.fab_satisfied < report.num_elites,
        "fab did not discriminate: {}/{}",
        report.fab_satisfied,
        report.num_elites
    );
    // the headline claim is Verified (matches + fabricable + sound).
    assert!(matches!(report.headline_color, Color::Verified { .. }));
    println!(
        "{{\"campaign\":\"grammarforge\",\"niches\":{},\"coverage\":{:.3},\"best_disc\":{:.4},\
         \"best_params\":{:?},\"simplified\":{},\"size\":{}->{},\"radius_threshold\":{:.4},\
         \"max_cert_err\":{:.4},\"max_sampled\":{:.4},\"status\":{},\"assessments\":{},\
         \"sound\":{},\"complete_sound\":{},\"fab_ok\":{}}}",
        report.num_elites,
        report.coverage,
        report.best_discrepancy,
        report.best_params,
        report.simplification.simplified_count(),
        report.simplification.size_before(),
        report.simplification.size_after(),
        report.simplification.radius_threshold(),
        report.simplification.max_certified_error(),
        report.simplification.max_sampled_discrepancy(),
        report.simplification.status().wire_code(),
        report.simplification.assessments(),
        report.simplification.is_sound(),
        report.simplification_sound,
        report.fab_satisfied,
    );
}

#[test]
fn the_local_threshold_and_compositional_envelope_are_distinct() {
    let prog = build_program(1.0, 1.0, 0.8, 0.02);
    let samples: Vec<[f64; 3]> = (0..5)
        .flat_map(|i| (0..5).map(move |j| [-2.0 + f64::from(i), (-1.0 + f64::from(j) * 0.5), 0.0]))
        .collect();
    let assessment = assess_simplification(&prog, 0.03, &samples);

    assert_eq!(assessment.radius_threshold().to_bits(), 0.03_f64.to_bits());
    assert_eq!(assessment.certified_error(), Some(0.04));
    assert_eq!(assessment.status(), SimplificationCheckStatus::Certified);
    assert_eq!(assessment.size_after() + 1, assessment.size_before());
    assert!(assessment.refusal().is_none());
    assert!(
        assessment.sampled_discrepancy().expect("finite witness") <= 0.04,
        "sampled witness exceeds the exact 2|radius| envelope: {assessment:#?}"
    );
}

#[test]
fn zero_drop_keep_and_sequential_threshold_cases_are_typed() {
    let samples = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let zero = assess_simplification(&Geom::sphere(1.0).offset(0.0), 0.03, &samples);
    let drop = assess_simplification(&Geom::sphere(1.0).offset(0.02), 0.03, &samples);
    let keep = assess_simplification(&Geom::sphere(1.0).offset(0.03), 0.03, &samples);
    let sequential =
        assess_simplification(&Geom::sphere(1.0).offset(0.02).offset(0.02), 0.03, &samples);

    assert_eq!(zero.certified_error(), Some(0.0));
    assert!(zero.size_after() < zero.size_before());
    assert_eq!(drop.certified_error(), Some(0.04));
    assert!(drop.size_after() < drop.size_before());
    assert_eq!(keep.certified_error(), Some(0.0));
    assert_eq!(keep.size_after(), keep.size_before(), "threshold is strict");
    assert_eq!(sequential.certified_error(), Some(0.08));
    for assessment in [&zero, &drop, &keep, &sequential] {
        assert_eq!(assessment.status(), SimplificationCheckStatus::Certified);
        assert!(
            assessment.sampled_discrepancy().expect("finite witness")
                <= assessment.certified_error().expect("finite certificate"),
            "{assessment:#?}"
        );
    }
}

#[test]
fn assessment_is_covariant_under_translation_sign_and_constructor_nesting() {
    let base_samples = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let translated_samples = [[3.0, -2.0, 1.0], [4.0, -2.0, 1.0]];
    let positive = assess_simplification(&Geom::sphere(1.0).offset(0.02), 0.03, &base_samples);
    let negative = assess_simplification(&Geom::sphere(1.0).offset(-0.02), 0.03, &base_samples);
    let translated = assess_simplification(
        &Geom::sphere(1.0).offset(0.02).translate([3.0, -2.0, 1.0]),
        0.03,
        &translated_samples,
    );
    let union_left = assess_simplification(
        &Geom::Empty.union(Geom::sphere(1.0).offset(0.02)),
        0.03,
        &base_samples,
    );
    let union_right = assess_simplification(
        &Geom::sphere(1.0).offset(0.02).union(Geom::Empty),
        0.03,
        &base_samples,
    );

    assert_eq!(positive.certified_error(), negative.certified_error());
    assert_eq!(positive.certified_error(), translated.certified_error());
    assert_eq!(positive.certified_error(), union_left.certified_error());
    assert_eq!(positive.certified_error(), union_right.certified_error());
    assert_eq!(
        positive.sampled_discrepancy(),
        negative.sampled_discrepancy()
    );
    assert_eq!(
        positive.sampled_discrepancy(),
        translated.sampled_discrepancy()
    );
    assert_eq!(
        positive.sampled_discrepancy(),
        union_left.sampled_discrepancy()
    );
    assert_eq!(
        positive.sampled_discrepancy(),
        union_right.sampled_discrepancy()
    );
    for assessment in [&positive, &negative, &translated, &union_left, &union_right] {
        assert_eq!(assessment.status(), SimplificationCheckStatus::Certified);
    }

    // Preserve the core P0 witness: both 0.006 offsets are still admitted at a
    // 0.01 local radius threshold even though their composed envelope is 0.024.
    let separated = Geom::sphere(1.0)
        .offset(0.006)
        .union(Geom::sphere(1.0).translate([100.0, 0.0, 0.0]))
        .offset(0.006);
    let nested = assess_simplification(&separated, 0.01, &[[0.0, 0.0, 0.0]]);
    assert_eq!(nested.certified_error(), Some(0.024));
    assert_eq!(nested.status(), SimplificationCheckStatus::Certified);
    assert!(nested.sampled_discrepancy().expect("finite witness") > 0.006);
}

#[test]
fn refusal_evidence_and_structural_empty_states_are_not_laundered() {
    let samples = [[0.0, 0.0, 0.0]];
    let refused = assess_simplification(&Geom::sphere(1.0), f64::NAN, &samples);
    assert_eq!(
        refused.status(),
        SimplificationCheckStatus::SimplifierRefused
    );
    assert!(matches!(
        refused.refusal(),
        Some(SimplifyRefusal::InvalidTolerance { tolerance_bits })
            if *tolerance_bits == f64::NAN.to_bits()
    ));
    assert!(refused.certified_error().is_none());
    assert!(refused.sampled_discrepancy().is_none());
    assert_eq!(
        refused,
        assess_simplification(&Geom::sphere(1.0), f64::NAN, &samples),
        "an identical typed refusal must replay bit-exactly"
    );

    let non_finite_program = assess_simplification(&Geom::sphere(f64::NAN), 0.03, &samples);
    assert_eq!(
        non_finite_program.status(),
        SimplificationCheckStatus::SimplifierRefused
    );
    assert_eq!(
        non_finite_program.refusal(),
        Some(&SimplifyRefusal::NonFiniteProgram)
    );

    let missing_evidence = assess_simplification(&Geom::sphere(1.0), 0.03, &[]);
    assert_eq!(
        missing_evidence.status(),
        SimplificationCheckStatus::DiscrepancyEvidenceRefused
    );
    assert_eq!(missing_evidence.certified_error(), Some(0.0));
    assert!(missing_evidence.sampled_discrepancy().is_none());
    let non_finite_evidence =
        assess_simplification(&Geom::sphere(1.0), 0.03, &[[f64::MAX, f64::MAX, f64::MAX]]);
    assert_eq!(
        non_finite_evidence.status(),
        SimplificationCheckStatus::DiscrepancyEvidenceRefused
    );

    let structural = assess_simplification(&Geom::Empty.offset(0.02), 0.03, &samples);
    assert_eq!(
        structural.status(),
        SimplificationCheckStatus::StructuralEmptyAgreement
    );
    assert_eq!(structural.certified_error(), Some(0.0));
    assert_eq!(structural.sampled_discrepancy(), Some(0.0));
    assert!(
        structural.size_after() < structural.size_before(),
        "the structural-empty rewrite should be exact and simplifying"
    );

    let mut structural_only = SimplificationSummary::new(0.03);
    structural_only.observe(&structural);
    assert_eq!(
        structural_only.status(),
        SimplificationCheckStatus::StructuralEmptyAgreement
    );
    assert!(structural_only.is_sound());
    assert!(structural_only.is_complete_and_sound(1));
    assert!(!structural_only.is_complete_and_sound(2));
    let finite = assess_simplification(&Geom::sphere(1.0), 0.03, &samples);
    structural_only.observe(&finite);
    assert_eq!(
        structural_only.status(),
        SimplificationCheckStatus::Certified,
        "mixed finite and structural success is the generic certified state"
    );
    assert_eq!(structural_only.structural_empty_agreements(), 1);
    assert!(structural_only.is_sound());
    assert!(structural_only.is_complete_and_sound(2));

    let mut summary = SimplificationSummary::new(0.03);
    summary.observe(&refused);
    summary.observe(&non_finite_program);
    summary.observe(&missing_evidence);
    summary.observe(&non_finite_evidence);
    summary.observe(&structural);
    assert!(!summary.is_sound());
    assert_eq!(summary.simplifier_refusals(), 2);
    assert_eq!(summary.discrepancy_evidence_refusals(), 2);
    assert_eq!(summary.structural_empty_agreements(), 1);

    let mut mismatched = SimplificationSummary::new(0.03);
    let other_threshold = assess_simplification(&Geom::sphere(1.0), 0.02, &samples);
    mismatched.observe(&other_threshold);
    assert_eq!(mismatched.threshold_mismatches(), 1);
    assert_eq!(
        mismatched.status(),
        SimplificationCheckStatus::ThresholdMismatch
    );
    assert!(!mismatched.is_sound());
}

#[test]
fn simplification_status_wire_codes_are_stable_and_unambiguous() {
    let statuses = [
        SimplificationCheckStatus::Certified,
        SimplificationCheckStatus::StructuralEmptyAgreement,
        SimplificationCheckStatus::SimplifierRefused,
        SimplificationCheckStatus::NonFiniteCertificate,
        SimplificationCheckStatus::NegativeCertificate,
        SimplificationCheckStatus::DiscrepancyEvidenceRefused,
        SimplificationCheckStatus::CertificateCheckExceeded,
        SimplificationCheckStatus::ThresholdMismatch,
    ];
    let codes: Vec<_> = statuses.iter().map(|status| status.wire_code()).collect();
    assert_eq!(codes, vec![0, 1, 2, 3, 4, 5, 6, 7]);
    assert!(SimplificationCheckStatus::Certified.is_sound());
    assert!(SimplificationCheckStatus::StructuralEmptyAgreement.is_sound());
    assert!(statuses[2..].iter().all(|status| !status.is_sound()));
    assert!(
        !SimplificationSummary::new(0.03).is_sound(),
        "an empty campaign must not be vacuously certified"
    );
    assert!(
        !SimplificationSummary::new(0.03).is_complete_and_sound(0),
        "zero expected records must not authorize a vacuous claim"
    );
    assert_ne!(
        SimplificationSummary::new(0.0),
        SimplificationSummary::new(-0.0),
        "signed-zero threshold provenance is bit-exact"
    );
}

#[test]
fn the_target_matches_itself_exactly() {
    let t = target();
    let samples: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.8, 0.0, 0.0], [1.5, 0.0, 0.0]];
    assert!(max_sdf_discrepancy(&t, &t, &samples).abs() < 1e-12);
}

#[test]
fn the_campaign_is_deterministic() {
    let a = run_campaign(0.2, 0.03);
    let b = run_campaign(0.2, 0.03);
    assert_eq!(a.num_elites, b.num_elites);
    assert_eq!(a.best_discrepancy.to_bits(), b.best_discrepancy.to_bits());
    assert_eq!(a.simplification, b.simplification);
}
