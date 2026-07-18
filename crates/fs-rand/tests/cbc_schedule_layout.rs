//! G0/G3 ratchets for CBC schema v4's single execution schedule and
//! requested-capacity storage authority (bead 6ys.20.4).

#![deny(unsafe_code)]

use fs_rand::cbc::{
    CBC_ADMISSION_SCHEMA_VERSION, CbcAdmissionError, CbcBudget, CbcExecutionMode, CbcProblem,
};
use fs_rand::cbc_cert::CbcPrefixCertificate;
use fs_rand::cbc_exec::{CbcControl, CbcExecError, CbcExecutor, CbcRunStatus, CbcTileShape};

#[test]
fn csl_001_exact_prime_composite_and_dimension_one_schedule_kats() {
    let dimension_one = CbcProblem::new(3, 1)
        .expect("structural fixture")
        .estimate()
        .expect("finite fixture");
    assert_eq!(dimension_one.admissible_candidates_per_prefix(), 2);
    assert_eq!(dimension_one.candidate_count(), 0);
    assert_eq!(dimension_one.candidate_visits(), 0);
    assert_eq!(dimension_one.limb_work_units(), 30);
    assert_eq!(dimension_one.scalar_work_units(), 99);
    assert_eq!(dimension_one.work_units(), 129);

    let prime = CbcProblem::new(5, 6)
        .expect("structural fixture")
        .estimate()
        .expect("finite fixture");
    assert_eq!(prime.admissible_candidates_per_prefix(), 4);
    assert_eq!(prime.candidate_count(), 20);
    assert_eq!(prime.admissible_candidate_count(), 20);
    assert_eq!(prime.candidate_visits(), 100);
    assert_eq!(prime.work_units(), 6_309);

    let composite = CbcProblem::new(8, 3)
        .expect("structural fixture")
        .estimate()
        .expect("finite fixture");
    assert_eq!(composite.admissible_candidates_per_prefix(), 4);
    assert_eq!(composite.candidate_count(), 14);
    assert_eq!(composite.admissible_candidate_count(), 8);
    assert_eq!(composite.candidate_visits(), 64);
    assert_eq!(composite.lattice_visits(), 88);
    assert_eq!(composite.limb_work_units(), 594);
    assert_eq!(composite.scalar_work_units(), 3_105);
    assert_eq!(composite.work_units(), 3_699);
}

#[test]
fn csl_002_certified_mode_charges_its_full_retained_graph_and_work() {
    let problem = CbcProblem::new(8, 3).expect("structural fixture");
    let plain = problem.estimate().expect("plain estimate");
    let certified = problem
        .estimate_for(CbcExecutionMode::Certified)
        .expect("certified estimate");

    let certificate_owners = 2_u128
        * u128::try_from(core::mem::size_of::<CbcPrefixCertificate>())
            .expect("layout size fits u128");
    assert_eq!(plain.certificate_work_units(), 0);
    assert_eq!(plain.certificate_retained_bytes(), 0);
    assert_eq!(certified.certificate_work_units(), 43);
    assert_eq!(certified.work_units(), 3_742);
    assert_eq!(certified.certificate_owner_bytes(), certificate_owners);
    assert_eq!(certified.certificate_prefix_payload_bytes(), 20);
    assert_eq!(certified.certificate_score_payload_bytes(), 48);
    assert_eq!(certified.certificate_tie_payload_bytes(), 32);
    assert_eq!(
        certified.certificate_retained_bytes(),
        certificate_owners + 100
    );
    assert_eq!(
        certified.candidate_phase_bytes(),
        plain.candidate_phase_bytes() + certified.certificate_retained_bytes() + 12,
        "certified scan adds the retained graph and third live score"
    );
}

#[test]
fn csl_003_mode_and_budget_boundaries_fail_closed() {
    let problem = CbcProblem::new(8, 3).expect("structural fixture");
    let estimate = problem
        .estimate_for(CbcExecutionMode::Certified)
        .expect("certified estimate");
    let exact = CbcBudget::new(estimate.work_units(), estimate.logical_state_bytes());
    let admission = problem
        .admit_for(CbcExecutionMode::Certified, exact)
        .expect("exact certified envelope admits");
    assert_eq!(CBC_ADMISSION_SCHEMA_VERSION, 4);
    assert_eq!(admission.mode(), CbcExecutionMode::Certified);

    assert_eq!(
        problem.admit_for(
            CbcExecutionMode::Certified,
            CbcBudget::new(estimate.work_units() - 1, estimate.logical_state_bytes()),
        ),
        Err(CbcAdmissionError::WorkBudgetExceeded {
            required: estimate.work_units(),
            available: estimate.work_units() - 1,
        })
    );
    assert_eq!(
        problem.admit_for(
            CbcExecutionMode::Certified,
            CbcBudget::new(estimate.work_units(), estimate.logical_state_bytes() - 1),
        ),
        Err(CbcAdmissionError::MemoryBudgetExceeded {
            required: estimate.logical_state_bytes(),
            available: estimate.logical_state_bytes() - 1,
        })
    );

    let plain_admission = problem
        .admit(CbcBudget::UNBOUNDED)
        .expect("plain mode admits");
    let mut plain = CbcExecutor::new(plain_admission).expect("fresh authority");
    assert_eq!(
        plain.enable_certificates(),
        Err(CbcExecError::CertificatesNotAdmitted)
    );

    let mut certified = CbcExecutor::new(admission).expect("fresh certified authority");
    certified
        .enable_certificates()
        .expect("certificate capability was admitted");
    let mut keep_going = || CbcControl::Continue;
    assert_eq!(
        certified
            .run(
                &mut keep_going,
                CbcTileShape::new(2, 3).expect("static tile"),
                u128::MAX,
            )
            .expect("certified execution follows its sealed schedule"),
        CbcRunStatus::Completed
    );
    assert_eq!(certified.work_spent(), estimate.work_units());
    assert_eq!(certified.certificates().len(), 2);
    let observed = certified.storage_observation();
    assert_eq!(observed.requested_certificate_records(), 2);
    assert_eq!(observed.retained_certificate_records(), 2);
    assert!(observed.observed_certificate_record_capacity() >= 2);
    assert_eq!(
        observed.retained_certificate_prefix_words(),
        usize::try_from(estimate.certificate_prefix_payload_bytes() / 4)
            .expect("fixture payload fits usize")
    );
    assert!(
        observed.observed_certificate_prefix_capacity_words()
            >= observed.retained_certificate_prefix_words()
    );
    assert!(
        observed.retained_certificate_score_limbs()
            <= usize::try_from(estimate.certificate_score_payload_bytes() / 4)
                .expect("fixture payload fits usize")
    );
    assert!(
        observed.observed_certificate_score_capacity_limbs()
            >= observed.retained_certificate_score_limbs()
    );
    assert!(
        observed.retained_certificate_tie_words()
            <= usize::try_from(estimate.certificate_tie_payload_bytes() / 4)
                .expect("fixture payload fits usize")
    );
    assert!(
        observed.observed_certificate_tie_capacity_words()
            >= observed.retained_certificate_tie_words()
    );
}

#[test]
fn csl_004_certified_mode_cannot_start_certificate_tracking_late() {
    let problem = CbcProblem::new(5, 3).expect("structural fixture");
    let admission = problem
        .admit_for(CbcExecutionMode::Certified, CbcBudget::UNBOUNDED)
        .expect("certified mode admits");
    let mut executor = CbcExecutor::new(admission).expect("fresh authority");
    let mut keep_going = || CbcControl::Continue;
    assert!(matches!(
        executor
            .run(
                &mut keep_going,
                CbcTileShape::new(1, 1).expect("static tile"),
                1,
            )
            .expect("a partial run remains replayable"),
        CbcRunStatus::AllowanceExhausted(_)
    ));
    assert!(executor.work_spent() > 0);
    assert_eq!(
        executor.enable_certificates(),
        Err(CbcExecError::CertificatesAfterStart)
    );
}
