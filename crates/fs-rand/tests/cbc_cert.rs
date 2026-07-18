//! Battery for the exact-CBC per-prefix certificates and their independent
//! checker (bead 6ys.20, certificate tranche): certified runs stay
//! byte-identical to uncertified ones, every emitted certificate passes
//! both checker modes, certificates are invariant under tile shape and
//! pause/resume, structural mirror ties are really captured, and seeded
//! tampering of every bound field fails closed in the named error class.

use fs_rand::cbc::{CbcAdmission, CbcBudget, CbcExecutionMode, CbcProblem};
use fs_rand::cbc_cert::{
    ADMISSIBLE_RULE_UNITS, CBC_CERTIFICATE_SCHEMA_VERSION, CbcCertError, CbcPrefixCertificate,
    TIE_RULE_LOWEST_CANDIDATE, audit_minimality, audit_minimality_admitted, verify_consistency,
    verify_consistency_admitted,
};
use fs_rand::cbc_exec::{CbcControl, CbcExecutor, CbcRunStatus, CbcTileShape};
use fs_rand::qmc::Lattice;

const CASES: [(u32, usize); 4] = [(5, 3), (8, 3), (16, 4), (127, 4)];

fn certified_admission(n: u32, dim: usize) -> CbcAdmission {
    let problem = CbcProblem::new(n, dim).expect("battery cases are structurally valid");
    problem
        .admit_for(CbcExecutionMode::Certified, CbcBudget::UNBOUNDED)
        .expect("battery cases admit under the unbounded budget")
}

fn certified_run(n: u32, dim: usize, tile: CbcTileShape) -> CbcExecutor {
    let admission = certified_admission(n, dim);
    let mut executor =
        CbcExecutor::new(admission).expect("fresh certified admission matches authority");
    executor
        .enable_certificates()
        .expect("certificates enabled before any work");
    let mut keep_going = || CbcControl::Continue;
    let status = executor
        .run(&mut keep_going, tile, u128::MAX)
        .expect("certified runs must stay within the admitted schedule");
    assert_eq!(status, CbcRunStatus::Completed);
    assert_eq!(
        executor.work_spent(),
        admission.estimate().work_units(),
        "certified execution must consume its exact admitted schedule"
    );
    executor
}

#[test]
fn crt_001_certified_runs_match_the_authority_and_both_checker_modes() {
    for (n, dim) in CASES {
        let reference = Lattice::cbc(n, dim);
        let checker_admission = certified_admission(n, dim);
        let executor = certified_run(n, dim, CbcTileShape::new(3, 7).expect("static tile"));
        let certificates = executor.certificates().to_vec();
        assert_eq!(
            certificates.len(),
            dim - 1,
            "n={n} dim={dim}: one certificate per scanned component"
        );
        for (index, certificate) in certificates.iter().enumerate() {
            assert_eq!(certificate.point_count, n);
            assert_eq!(
                certificate.prefix,
                reference.z[..index + 2].to_vec(),
                "n={n}: certificate {index} binds the authority prefix"
            );
            assert_eq!(certificate.chosen(), reference.z[index + 1]);
            verify_consistency(certificate).unwrap_or_else(|error| {
                panic!("n={n} certificate {index} failed consistency: {error:?}")
            });
            audit_minimality(certificate).unwrap_or_else(|error| {
                panic!("n={n} certificate {index} failed the full audit: {error:?}")
            });
            verify_consistency_admitted(checker_admission, certificate).unwrap_or_else(|error| {
                panic!("n={n} certificate {index} failed admitted consistency: {error:?}")
            });
            audit_minimality_admitted(checker_admission, certificate).unwrap_or_else(|error| {
                panic!("n={n} certificate {index} failed the admitted audit: {error:?}")
            });
        }
        assert_eq!(
            executor
                .into_lattice()
                .expect("completed construction yields the lattice")
                .z,
            reference.z,
            "n={n} dim={dim}: certifying changed the constructed bytes"
        );
    }
    assert_eq!(CBC_CERTIFICATE_SCHEMA_VERSION, 2);
}

#[test]
fn crt_002_certificates_are_invariant_under_tiling_and_pause_resume() {
    let reference = certified_run(16, 4, CbcTileShape::new(1, 1).expect("static tile"))
        .certificates()
        .to_vec();

    for (candidate_block, point_block) in [(2, 4), (64, 64)] {
        let executor = certified_run(
            16,
            4,
            CbcTileShape::new(candidate_block, point_block).expect("nonzero tile"),
        );
        assert_eq!(
            executor.certificates(),
            reference.as_slice(),
            "tile ({candidate_block},{point_block}) changed certificate bytes"
        );
    }

    // Sliced allowances (pause/resume) must also reproduce them exactly.
    let problem = CbcProblem::new(16, 4).expect("structurally valid");
    let admission = problem
        .admit_for(CbcExecutionMode::Certified, CbcBudget::UNBOUNDED)
        .expect("admits unbounded");
    let mut executor =
        CbcExecutor::new(admission).expect("fresh certified admission matches authority");
    executor
        .enable_certificates()
        .expect("certificates enabled before any work");
    let tile = CbcTileShape::new(2, 4).expect("static tile");
    let mut keep_going = || CbcControl::Continue;
    let mut guard = 0_u32;
    loop {
        guard += 1;
        assert!(guard < 1_000_000, "allowance loop failed to converge");
        match executor
            .run(&mut keep_going, tile, 50_000)
            .expect("sliced certified runs stay within the schedule")
        {
            CbcRunStatus::Completed => break,
            CbcRunStatus::AllowanceExhausted(_) => {}
            CbcRunStatus::Cancelled(_) => panic!("nothing requested cancellation"),
        }
    }
    assert_eq!(
        executor.certificates(),
        reference.as_slice(),
        "pause/resume changed certificate bytes"
    );
}

#[test]
fn crt_003_mirror_ties_are_structural_and_lowest_candidate_wins() {
    // Mirror symmetry: candidate c and n−c share a residue multiset, so
    // every certificate's tie class pairs c with n−c (unless c is
    // self-mirrored). The chosen component is always the class minimum.
    for (n, dim) in CASES {
        let executor = certified_run(n, dim, CbcTileShape::new(3, 7).expect("static tile"));
        for certificate in executor.certificates() {
            let chosen = certificate.chosen();
            assert_eq!(
                certificate.tie_class.first(),
                Some(&chosen),
                "n={n}: the chosen candidate must be the tie-class minimum"
            );
            for &member in &certificate.tie_class {
                let mirror = n - member;
                assert!(
                    certificate.tie_class.contains(&mirror),
                    "n={n}: tie class {:?} lacks the mirror {mirror} of {member}",
                    certificate.tie_class
                );
            }
            assert!(
                certificate.tie_class.len() >= 2 || chosen * 2 == n,
                "n={n}: a lone tie class is only possible for a self-mirrored candidate"
            );
        }
    }
}

fn baseline_certificate() -> CbcPrefixCertificate {
    certified_run(16, 4, CbcTileShape::new(3, 7).expect("static tile")).certificates()[0].clone()
}

#[test]
fn crt_004_tampering_fails_closed_in_the_named_class() {
    let good = baseline_certificate();
    verify_consistency(&good).expect("the baseline certificate is green");
    audit_minimality(&good).expect("the baseline certificate audits green");

    // A different chosen candidate: prefix binding breaks.
    let mut tampered = good.clone();
    let position = tampered.prefix.len() - 1;
    tampered.prefix[position] = 3; // coprime to 16, not the winner
    assert!(
        matches!(
            verify_consistency(&tampered),
            Err(CbcCertError::MalformedTieClass | CbcCertError::TieClassScoreMismatch { .. })
        ),
        "swapping the chosen candidate must refuse"
    );

    // A corrupted winning limb.
    let mut tampered = good.clone();
    tampered.winning_score_limbs[0] ^= 1;
    assert!(matches!(
        verify_consistency(&tampered),
        Err(CbcCertError::TieClassScoreMismatch { .. })
    ));
    assert!(matches!(
        audit_minimality(&tampered),
        Err(CbcCertError::NotMinimal { .. } | CbcCertError::TieClassIncomplete)
    ));

    // A non-coprime tie-class member.
    let mut tampered = good.clone();
    tampered.tie_class.push(good.point_count / 2);
    tampered.tie_class.sort_unstable();
    assert!(matches!(
        verify_consistency(&tampered),
        Err(CbcCertError::MalformedTieClass)
    ));

    // A dropped tie-class member: consistency may pass (the remaining
    // members really do score the winning value) but the full audit must
    // catch the incompleteness.
    let mut tampered = good.clone();
    if tampered.tie_class.len() > 1 {
        tampered.tie_class.pop();
        assert!(matches!(
            audit_minimality(&tampered),
            Err(CbcCertError::TieClassIncomplete)
        ));
    }

    // A forged runner-up score.
    let mut tampered = good.clone();
    if let Some((limbs, _)) = &mut tampered.runner_up {
        limbs[0] ^= 1;
        assert!(matches!(
            verify_consistency(&tampered),
            Err(CbcCertError::RunnerUpMismatch)
        ));
    }

    // A wrong denominator derivation.
    let mut tampered = good.clone();
    tampered.denominator_exponent += 1;
    assert_eq!(
        verify_consistency(&tampered),
        Err(CbcCertError::DenominatorMismatch)
    );

    // A foreign rule token.
    let mut tampered = good.clone();
    tampered.tie_rule = "highest-candidate-wins";
    assert_eq!(
        verify_consistency(&tampered),
        Err(CbcCertError::UnknownRule)
    );

    // A theorem-component certificate (prefix too short).
    let mut tampered = good.clone();
    tampered.prefix.truncate(1);
    tampered.denominator_exponent = 1;
    assert_eq!(
        verify_consistency(&tampered),
        Err(CbcCertError::MalformedPrefix)
    );

    // Rule tokens are the schema's declared ones.
    assert_eq!(good.tie_rule, TIE_RULE_LOWEST_CANDIDATE);
    assert_eq!(good.admissible_rule, ADMISSIBLE_RULE_UNITS);
}

#[test]
fn crt_005_uncertified_runs_emit_nothing_and_late_enabling_refuses() {
    let problem = CbcProblem::new(8, 3).expect("structurally valid");
    let admission = problem
        .admit(CbcBudget::UNBOUNDED)
        .expect("admits unbounded");
    let mut executor =
        CbcExecutor::new(admission).expect("fresh construction admission matches authority");
    let mut keep_going = || CbcControl::Continue;
    let tile = CbcTileShape::new(2, 4).expect("static tile");
    let status = executor
        .run(&mut keep_going, tile, u128::MAX)
        .expect("uncertified run completes");
    assert_eq!(status, CbcRunStatus::Completed);
    assert!(
        executor.certificates().is_empty(),
        "uncertified runs must not allocate certificates"
    );
    assert!(
        executor.enable_certificates()
            == Err(fs_rand::cbc_exec::CbcExecError::CertificatesNotAdmitted),
        "enabling after work must refuse (certificates cover all or none)"
    );
}

#[test]
fn crt_006_admitted_checkers_refuse_wrong_authority_and_hostile_lengths_first() {
    let good = baseline_certificate();
    let admission = certified_admission(good.point_count, 4);
    verify_consistency_admitted(admission, &good).expect("current receipt covers consistency");
    audit_minimality_admitted(admission, &good).expect("current receipt covers the full audit");

    let construction_only = CbcProblem::new(good.point_count, 4)
        .expect("fixture is structural")
        .admit(CbcBudget::UNBOUNDED)
        .expect("fixture admits");
    assert_eq!(
        verify_consistency_admitted(construction_only, &good),
        Err(CbcCertError::AdmissionMismatch),
        "a construction receipt cannot authorize certificate checking"
    );

    let wrong_point_count = certified_admission(good.point_count + 1, 4);
    assert_eq!(
        verify_consistency_admitted(wrong_point_count, &good),
        Err(CbcCertError::AdmissionMismatch)
    );

    let short_problem = certified_admission(good.point_count, good.prefix.len() - 1);
    assert_eq!(
        verify_consistency_admitted(short_problem, &good),
        Err(CbcCertError::AdmissionEnvelopeExceeded {
            storage_class: "prefix words",
            required: u128::try_from(good.prefix.len()).expect("usize fits u128"),
            admitted: u128::try_from(good.prefix.len() - 1).expect("usize fits u128"),
        }),
        "the receipt must cover the complete declared prefix"
    );

    let score_limit = usize::try_from(admission.estimate().score_capacity_limbs())
        .expect("the admitted target bound fits usize");
    let mut oversized_winner = good.clone();
    oversized_winner
        .winning_score_limbs
        .resize(score_limit + 1, 1);
    assert_eq!(
        verify_consistency_admitted(admission, &oversized_winner),
        Err(CbcCertError::AdmissionEnvelopeExceeded {
            storage_class: "winning score limbs",
            required: u128::try_from(score_limit + 1).expect("usize fits u128"),
            admitted: admission.estimate().score_capacity_limbs(),
        })
    );

    let tie_limit = usize::try_from(admission.estimate().admissible_candidates_per_prefix())
        .expect("the admitted target bound fits usize");
    let mut oversized_ties = good.clone();
    oversized_ties.tie_class.resize(tie_limit + 1, 1);
    assert_eq!(
        audit_minimality_admitted(admission, &oversized_ties),
        Err(CbcCertError::AdmissionEnvelopeExceeded {
            storage_class: "tie-class words",
            required: u128::try_from(tie_limit + 1).expect("usize fits u128"),
            admitted: admission.estimate().admissible_candidates_per_prefix(),
        }),
        "envelope refusal must precede semantic tie-class parsing"
    );

    let mut oversized_runner = good.clone();
    oversized_runner.runner_up = Some((vec![1; score_limit + 1], 1));
    assert_eq!(
        verify_consistency_admitted(admission, &oversized_runner),
        Err(CbcCertError::AdmissionEnvelopeExceeded {
            storage_class: "runner-up score limbs",
            required: u128::try_from(score_limit + 1).expect("usize fits u128"),
            admitted: admission.estimate().score_capacity_limbs(),
        })
    );
}
