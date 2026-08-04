//! Outcome, censoring, refinement, and calibration-admission regressions.

use fs_euler_disc_e2e::convergence::{
    CalibrationEvidenceKind, CalibrationReadinessError, CalibrationReadinessInput,
    CensorAwareDurationOrdering, CensorAwareRankingRefusal, ConvergenceScales, DeclaredEvidence,
    HorizonContinuationPolicy, ObservedDurationOrdering, ObservedOrder, OrderUnavailableReason,
    RefinementMode, RunOutcome, ThreeRungConvergence, admit_calibration_readiness,
    analyse_three_rung_convergence, classify_outcome, compare_censor_aware_durations,
    compare_observed_durations,
};
use fs_euler_disc_e2e::coupled_runner::{
    CoupledControls, CoupledFactors, CoupledInitialState, CoupledNumericalRefusalReason,
    run_closed_reduced,
};

fn factors() -> CoupledFactors {
    let radius_m = 0.038;
    let thickness_m = 0.006;
    let mass_kg = std::f64::consts::PI * radius_m * radius_m * thickness_m * 2_680.0;
    CoupledFactors {
        mass_kg,
        radius_m,
        thickness_m,
        transverse_inertia_kg_m2: mass_kg * (3.0 * radius_m * radius_m + thickness_m.powi(2))
            / 12.0,
        axial_inertia_kg_m2: 0.5 * mass_kg * radius_m * radius_m,
        gravity_m_per_s2: 9.806_65,
        sliding_friction_coefficient: 0.42,
        rolling_resistance_m: 4.0e-5,
        contact_stiffness_n_per_m: 8.0e4,
        contact_damping_n_s_per_m: 3.0,
        base_effective_mass_kg: 0.25,
        base_stiffness_n_per_m: 4.0e4,
        base_damping_n_s_per_m: 4.0,
        gas_rotational_damping_n_m_s: 2.0e-7,
        gas_translation_damping_n_s_per_m: 4.0e-4,
    }
}

#[test]
fn manual_non_finite_or_negative_terminal_times_refuse_duration_comparison() {
    let valid = RunOutcome::PhysicalTerminal {
        kind: fs_euler_disc_e2e::convergence::PhysicalTerminal::InclinationThreshold,
        event_time_s: 1.0,
    };
    let non_finite = RunOutcome::PhysicalTerminal {
        kind: fs_euler_disc_e2e::convergence::PhysicalTerminal::InclinationThreshold,
        event_time_s: f64::NAN,
    };
    let negative = RunOutcome::PhysicalTerminal {
        kind: fs_euler_disc_e2e::convergence::PhysicalTerminal::InclinationThreshold,
        event_time_s: -1.0,
    };
    assert_eq!(
        compare_observed_durations(non_finite, valid),
        Err(fs_euler_disc_e2e::convergence::RankingRefusal::InvalidLeftObservedTerminalTime)
    );
    assert_eq!(
        compare_observed_durations(valid, negative),
        Err(fs_euler_disc_e2e::convergence::RankingRefusal::InvalidRightObservedTerminalTime)
    );
}

fn observed(time_s: f64) -> RunOutcome {
    RunOutcome::PhysicalTerminal {
        kind: fs_euler_disc_e2e::convergence::PhysicalTerminal::InclinationThreshold,
        event_time_s: time_s,
    }
}

fn censored(lower_bound_s: f64) -> RunOutcome {
    RunOutcome::RightCensored {
        censor_time_s: lower_bound_s,
    }
}

#[test]
fn censor_aware_ordering_proves_only_logically_established_duration_relations() {
    assert_eq!(
        compare_censor_aware_durations(observed(1.0), observed(2.0)),
        Ok(CensorAwareDurationOrdering::ProvenLeftShorter)
    );
    assert_eq!(
        compare_censor_aware_durations(observed(2.0), observed(1.0)),
        Ok(CensorAwareDurationOrdering::ProvenLeftLonger)
    );
    assert_eq!(
        compare_censor_aware_durations(observed(1.0), observed(1.0)),
        Ok(CensorAwareDurationOrdering::EqualObserved)
    );
    assert_eq!(
        compare_censor_aware_durations(observed(0.8), censored(1.0)),
        Ok(CensorAwareDurationOrdering::ProvenLeftShorter)
    );
    assert_eq!(
        compare_censor_aware_durations(censored(1.2), observed(1.0)),
        Ok(CensorAwareDurationOrdering::ProvenLeftLonger)
    );
}

#[test]
fn censor_aware_ordering_retains_overlapping_bounds_as_indeterminate() {
    assert_eq!(
        compare_censor_aware_durations(observed(1.0), censored(1.0)),
        Ok(CensorAwareDurationOrdering::Indeterminate)
    );
    assert_eq!(
        compare_censor_aware_durations(censored(1.0), observed(1.0)),
        Ok(CensorAwareDurationOrdering::Indeterminate)
    );
    assert_eq!(
        compare_censor_aware_durations(censored(2.0), censored(1.0)),
        Ok(CensorAwareDurationOrdering::Indeterminate)
    );
}

#[test]
fn censor_aware_ordering_refuses_numerical_failures() {
    let refusal = RunOutcome::NumericalRefusal {
        last_valid_time_s: 0.5,
        reason: CoupledNumericalRefusalReason::ReimpactLimitExceeded,
    };
    assert_eq!(
        compare_censor_aware_durations(refusal, observed(1.0)),
        Err(CensorAwareRankingRefusal::LeftNumericalRefusal)
    );
    assert_eq!(
        compare_censor_aware_durations(observed(1.0), refusal),
        Err(CensorAwareRankingRefusal::RightNumericalRefusal)
    );
}

fn initial() -> CoupledInitialState {
    CoupledInitialState {
        inclination_rad: 0.08,
        precession_rad_per_s: 16.0,
        spin_rad_per_s: 120.0,
    }
}

fn run(timestep_s: f64, maximum_steps: u32) -> fs_euler_disc_e2e::coupled_runner::CoupledRun {
    run_closed_reduced(
        factors(),
        CoupledControls {
            timestep_s,
            maximum_steps,
            terminal_inclination_rad: 0.002,
            reimpact_limit: 128,
        },
        initial(),
        None,
    )
    .expect("reduced run")
}

fn scales() -> ConvergenceScales {
    ConvergenceScales {
        inclination_rad: 1.0,
        precession_rad_per_s: 200.0,
        spin_rad_per_s: 200.0,
        work_j: 1.0,
        energy_j: 1.0,
    }
}

#[test]
fn horizon_is_right_censored_and_cannot_be_duration_ranked() {
    let censored = run(2.0e-5, 8);
    let outcome = classify_outcome(&censored).expect("outcome");
    assert!(matches!(outcome, RunOutcome::RightCensored { .. }));
    assert_eq!(
        compare_observed_durations(outcome, outcome),
        Err(
            fs_euler_disc_e2e::convergence::RankingRefusal::LeftNotObservedPhysicalTerminal {
                class: fs_euler_disc_e2e::convergence::OutcomeClass::RightCensored,
            }
        )
    );
    let policy = HorizonContinuationPolicy {
        initial_horizon_s: 0.1,
        maximum_horizon_s: 0.4,
        multiplier: 2.0,
        maximum_extensions: 2,
    };
    let policy_outcome = RunOutcome::RightCensored { censor_time_s: 0.1 };
    assert_eq!(
        policy
            .next_horizon_s(policy_outcome, 0.1, 0)
            .expect("next horizon"),
        Some(0.2)
    );
    let final_policy_outcome = RunOutcome::RightCensored { censor_time_s: 0.4 };
    assert_eq!(
        policy
            .next_horizon_s(final_policy_outcome, 0.4, 1)
            .expect("final horizon"),
        None
    );
}

#[test]
fn continuation_refuses_mismatched_or_non_finite_censor_times() {
    let policy = HorizonContinuationPolicy {
        initial_horizon_s: 0.1,
        maximum_horizon_s: 0.4,
        multiplier: 2.0,
        maximum_extensions: 2,
    };
    assert_eq!(
        policy.next_horizon_s(
            RunOutcome::RightCensored {
                censor_time_s: 0.15,
            },
            0.1,
            0,
        ),
        Err(fs_euler_disc_e2e::convergence::HorizonPolicyError::CensorTimeDoesNotMatchCurrentHorizon)
    );
    assert_eq!(
        policy.next_horizon_s(
            RunOutcome::RightCensored {
                censor_time_s: f64::INFINITY,
            },
            0.1,
            0,
        ),
        Err(fs_euler_disc_e2e::convergence::HorizonPolicyError::InvalidCensorTime)
    );
}

#[test]
fn observed_physical_terminals_are_comparable_only_when_both_observed() {
    let outcome = RunOutcome::PhysicalTerminal {
        kind: fs_euler_disc_e2e::convergence::PhysicalTerminal::InclinationThreshold,
        event_time_s: 0.25,
    };
    assert_eq!(
        compare_observed_durations(outcome, outcome),
        Ok(ObservedDurationOrdering::Equal)
    );
}

#[test]
fn three_rung_analysis_retains_deltas_but_withholds_order_for_eventful_mode() {
    let coarse = run(4.0e-5, 8);
    let fine = run(2.0e-5, 16);
    let reference = run(1.0e-5, 32);
    let receipt = analyse_three_rung_convergence(ThreeRungConvergence {
        coarse: &coarse,
        fine: &fine,
        reference: &reference,
        coarse_timestep_s: 4.0e-5,
        fine_timestep_s: 2.0e-5,
        reference_timestep_s: 1.0e-5,
        mode: RefinementMode::Eventful {
            reason: "contact/reimpact classification requires the event lane".to_owned(),
        },
        scales: scales(),
    })
    .expect("three rung analysis");
    assert!(receipt.terminal_class_agreement);
    assert!(receipt.coarse_fine_qoi.inclination.is_finite());
    assert!(receipt.fine_reference_work_energy.energy_defect.is_finite());
    assert_eq!(
        receipt.observed_order,
        ObservedOrder::NotApplicable {
            reason: OrderUnavailableReason::NonSmoothOrUnresolvedMode,
        }
    );
}

#[test]
fn non_finite_normalized_delta_refuses_instead_of_publishing_a_receipt() {
    let mut coarse = run(4.0e-5, 8);
    let mut fine = run(2.0e-5, 16);
    let reference = run(1.0e-5, 32);
    coarse.samples.last_mut().expect("sample").spin_rad_per_s = f64::MAX;
    fine.samples.last_mut().expect("sample").spin_rad_per_s = -f64::MAX;
    assert_eq!(
        analyse_three_rung_convergence(ThreeRungConvergence {
            coarse: &coarse,
            fine: &fine,
            reference: &reference,
            coarse_timestep_s: 4.0e-5,
            fine_timestep_s: 2.0e-5,
            reference_timestep_s: 1.0e-5,
            mode: RefinementMode::Eventful {
                reason: "overflow adversary".to_owned(),
            },
            scales: scales(),
        }),
        Err(
            fs_euler_disc_e2e::convergence::ConvergenceError::NonFiniteNormalizedDelta {
                field: "spin",
            }
        )
    );
}

#[test]
fn calibration_readiness_refuses_missing_or_aliased_evidence_without_fitting() {
    let missing = CalibrationReadinessInput {
        specimen: DeclaredEvidence::Present {
            identity: "specimen-v1".to_owned(),
        },
        rig: DeclaredEvidence::Present {
            identity: "rig-v1".to_owned(),
        },
        instrument: DeclaredEvidence::Missing,
        raw_observations: DeclaredEvidence::Missing,
        observation_covariance: DeclaredEvidence::Missing,
        calibration_partition: DeclaredEvidence::Missing,
        blind_holdout: DeclaredEvidence::Missing,
    };
    assert_eq!(
        admit_calibration_readiness(missing),
        Err(CalibrationReadinessError::MissingEvidence {
            kind: CalibrationEvidenceKind::Instrument,
        })
    );

    let present = || DeclaredEvidence::Present {
        identity: "declared-artifact-v1".to_owned(),
    };
    let aliased = CalibrationReadinessInput {
        specimen: present(),
        rig: present(),
        instrument: present(),
        raw_observations: present(),
        observation_covariance: present(),
        calibration_partition: DeclaredEvidence::Present {
            identity: "partition-same".to_owned(),
        },
        blind_holdout: DeclaredEvidence::Present {
            identity: "partition-same".to_owned(),
        },
    };
    assert_eq!(
        admit_calibration_readiness(aliased),
        Err(CalibrationReadinessError::PartitionAlias)
    );
}

#[test]
fn calibration_readiness_accepts_distinct_declared_partitions_without_promoting_them() {
    let evidence = |identity: &str| DeclaredEvidence::Present {
        identity: identity.to_owned(),
    };
    let receipt = admit_calibration_readiness(CalibrationReadinessInput {
        specimen: evidence("specimen-v1"),
        rig: evidence("rig-v1"),
        instrument: evidence("instrument-calibration-v1"),
        raw_observations: evidence("raw-stream-v1"),
        observation_covariance: evidence("covariance-v1"),
        calibration_partition: evidence("calibration-partition-v1"),
        blind_holdout: evidence("blind-holdout-v1"),
    })
    .expect("structural readiness");
    assert_eq!(
        receipt.calibration_partition_identity,
        "calibration-partition-v1"
    );
    assert_eq!(receipt.blind_holdout_identity, "blind-holdout-v1");
}
