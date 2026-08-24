//! Horizon trigger 13b battery (bead `frankensim-epic-addendum-xpck.5.5`):
//! boundary, mutant, refusal, and receipt tests for Proposal 13b symmetry
//! prevalence and isotypic solver gate.

use fs_govern::horizon_symmetry::{
    evaluate_trigger_13b, mint_trigger_13b_receipt, SymmetryDisposition,
    Trigger13bReceipt, Trigger13bRefusal, Trigger13bVerdict, WorkloadSymmetryAssessment,
    SYMMETRY_PREVALENCE_MIN,
};

fn symmetric_workload(id: &str) -> WorkloadSymmetryAssessment {
    WorkloadSymmetryAssessment {
        workload_id: id.into(),
        group_name: "C2v".into(),
        group_order: 4,
        asymmetry_residual: 1e-6,
        full_solve_falsifier_available: true,
    }
}

fn asymmetric_workload(id: &str) -> WorkloadSymmetryAssessment {
    WorkloadSymmetryAssessment {
        workload_id: id.into(),
        group_name: "C1".into(),
        group_order: 1,
        asymmetry_residual: 0.0,
        full_solve_falsifier_available: true,
    }
}

#[test]
fn gate_activates_at_and_above_15_percent_prevalence() {
    // 2 out of 10 = 20% >= 15%: activates!
    let population = vec![
        symmetric_workload("wl-1"),
        symmetric_workload("wl-2"),
        asymmetric_workload("wl-3"),
        asymmetric_workload("wl-4"),
        asymmetric_workload("wl-5"),
        asymmetric_workload("wl-6"),
        asymmetric_workload("wl-7"),
        asymmetric_workload("wl-8"),
        asymmetric_workload("wl-9"),
        asymmetric_workload("wl-10"),
    ];
    assert_eq!(
        evaluate_trigger_13b(&population),
        Ok(Trigger13bVerdict::ActivateIsotypicSolver)
    );

    let receipt = mint_trigger_13b_receipt(Some(&population));
    assert_eq!(receipt.disposition, SymmetryDisposition::Activate);
    assert_eq!(receipt.verdict, Trigger13bVerdict::ActivateIsotypicSolver);
    assert_eq!(receipt.qualifying_count, 2);
    assert_eq!(receipt.total_count, 10);
    assert!((receipt.prevalence - 0.20).abs() < 1e-9);
}

#[test]
fn boundary_exactly_at_15_percent_activates() {
    // 3 out of 20 = 15%: activates!
    let mut population = Vec::new();
    for i in 0..3 {
        population.push(symmetric_workload(&format!("sym-{i}")));
    }
    for i in 3..20 {
        population.push(asymmetric_workload(&format!("asym-{i}")));
    }
    assert_eq!(
        evaluate_trigger_13b(&population),
        Ok(Trigger13bVerdict::ActivateIsotypicSolver)
    );
}

#[test]
fn boundary_one_below_15_percent_defers_to_detection_only() {
    // 1 out of 10 = 10% < 15%: detection only.
    let mut population = vec![symmetric_workload("wl-1")];
    for i in 2..=10 {
        population.push(asymmetric_workload(&format!("wl-{i}")));
    }
    assert_eq!(
        evaluate_trigger_13b(&population),
        Ok(Trigger13bVerdict::DetectionOnly)
    );

    let receipt = mint_trigger_13b_receipt(Some(&population));
    assert_eq!(receipt.disposition, SymmetryDisposition::DetectionOnly);
    assert_eq!(receipt.verdict, Trigger13bVerdict::DetectionOnly);
    assert_eq!(receipt.qualifying_count, 1);
    assert_eq!(receipt.total_count, 10);
    assert!((receipt.prevalence - 0.10).abs() < 1e-9);
}

#[test]
fn uncertified_wide_residual_disqualifies_workload() {
    let mut sym = symmetric_workload("wl-1");
    sym.asymmetry_residual = 0.05; // 5% residual >> 1e-4 cap
    assert!(!sym.is_qualifying());

    let population = vec![sym, asymmetric_workload("wl-2")];
    assert_eq!(
        evaluate_trigger_13b(&population),
        Ok(Trigger13bVerdict::DetectionOnly)
    );
}

#[test]
fn missing_falsifier_disqualifies_workload() {
    let mut sym = symmetric_workload("wl-1");
    sym.full_solve_falsifier_available = false; // No full-solve falsifier
    assert!(!sym.is_qualifying());
}

#[test]
fn malformed_workloads_refuse_by_name() {
    let empty: Vec<WorkloadSymmetryAssessment> = Vec::new();
    assert_eq!(
        evaluate_trigger_13b(&empty),
        Err(Trigger13bRefusal::EmptyPopulation)
    );

    let zero_group = vec![WorkloadSymmetryAssessment {
        workload_id: "wl-bad".into(),
        group_name: "invalid".into(),
        group_order: 0,
        asymmetry_residual: 0.0,
        full_solve_falsifier_available: true,
    }];
    assert!(matches!(
        evaluate_trigger_13b(&zero_group),
        Err(Trigger13bRefusal::ZeroGroupOrder { .. })
    ));
}

#[test]
fn mint_receipt_returns_nodata_when_population_absent() {
    let receipt: Trigger13bReceipt = mint_trigger_13b_receipt(None);
    assert_eq!(receipt.disposition, SymmetryDisposition::NoData);
    assert_eq!(receipt.verdict, Trigger13bVerdict::DetectionOnly);
    assert_eq!(receipt.qualifying_count, 0);
    assert_eq!(receipt.total_count, 0);
    assert_eq!(receipt.prevalence, 0.0);
}
