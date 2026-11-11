use super::*;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-affine-grid-refinement.fsim"
));

fn source() -> String {
    // The curved cut domain changes its complete affine response on enrichment.
    // Zero change tolerances force a new-grid experiment even when every
    // original physical cap is comfortably satisfied on the source grid.
    FIXTURE
        .replacen(":max-updates 80", ":max-updates 1", 1)
        .replace(":max-updates 20", ":max-updates 2")
        .replace(
            ":compliance-change-tolerance 0.02",
            ":compliance-change-tolerance 0.0",
        )
        .replace(
            ":response-change-tolerance 0.02",
            ":response-change-tolerance 0.0",
        )
        .replace(
            ":stress-change-tolerance 0.02",
            ":stress-change-tolerance 0.0",
        )
        .replace(
            ":inventory-change-tolerance 0.0001",
            ":inventory-change-tolerance 0.0",
        )
}

fn field<'a>(value: &'a JsonValue, key: &str) -> &'a JsonValue {
    value.get(key).unwrap_or_else(|| panic!("missing {key}"))
}

fn document(bytes: &[u8]) -> JsonValue {
    JsonValue::parse(std::str::from_utf8(bytes).unwrap()).unwrap()
}

fn report(ledger: &Ledger, receipt: &Loaded) -> JsonValue {
    document(&linked(ledger, &receipt.value, "report_json", "study-report-json").unwrap())
}

fn phase(receipt: &Loaded) -> &str {
    field(&receipt.value, "grid_refinement")
        .str_field("phase")
        .unwrap()
}

fn source_receipt(ledger: &Ledger, receipt: &Loaded) -> Loaded {
    let context = field(&receipt.value, "grid_refinement");
    load(
        ledger,
        &format!("study-{}", context.str_field("source_receipt").unwrap()),
    )
    .unwrap()
}

fn assert_native_state_equal(left: &JsonValue, right: &JsonValue) {
    for key in [
        "binding",
        "point",
        "best_feasible_density",
        "sample",
        "multiplier",
        "penalty",
        "spectral_step",
        "inner_tolerance",
        "previous_violation",
    ] {
        assert_eq!(field(left, key), field(right, key), "native {key}");
    }
}

fn assert_source_selected(ledger: &Ledger, source: &Loaded, current: &Loaded) {
    assert_eq!(
        linked(ledger, &current.value, "design", "study-design").unwrap(),
        linked(ledger, &source.value, "design", "study-design").unwrap(),
    );
    let original = report(ledger, source);
    let active = report(ledger, current);
    for key in [
        "selected",
        "selected_design",
        "selected_feasible",
        "grid_assessment",
    ] {
        assert_eq!(
            field(&active, key),
            field(&original, key),
            "installed {key}"
        );
    }
    let refinement = field(&active, "grid_refinement");
    assert_eq!(integer(refinement, "selected_generation").unwrap(), 0);
    assert_eq!(
        integer(field(refinement, "source_grid"), "iterations_completed").unwrap(),
        1
    );
    let checkpoint = document(&linked(ledger, &current.value, "checkpoint", resume::KIND).unwrap());
    assert_eq!(
        checkpoint.str_field("design"),
        current.value.str_field("working_design")
    );
}

fn assert_charged_after(spec: &Spec, before: &Loaded, after: &Loaded) {
    let before = Spent::read(&before.value).unwrap();
    let after = Spent::read(&after.value).unwrap();
    after.validate(spec).unwrap();
    assert!(after.wall_s >= before.wall_s);
    assert!(after.linear.linear_solves > before.linear.linear_solves);
    assert!(after.linear.linear_iterations > before.linear.linear_iterations);
    assert!(
        after.linear.preconditioner_operator_applications
            >= before.linear.preconditioner_operator_applications
    );
    assert!(
        after.linear.preconditioner_galerkin_products
            >= before.linear.preconditioner_galerkin_products
    );
    assert!(after.geometry.boxes > before.geometry.boxes);
    assert!(after.geometry.points > before.geometry.points);
    assert!(after.geometry.field_evaluations > before.geometry.field_evaluations);
}

#[test]
fn an_exhausted_stress_point_field_allowance_stops_the_final_assessment_as_a_budget() {
    // 100000 clears the source and candidate grids (10692 and 49248 point-fields) and not the finer grid of the final assessment (171072): the allowance runs out there, and that is a budget,
    // named by its stage, not a numerical failure.
    let spec = spec::parse(&source().replace(":max-stress-point-fields 200000", ":max-stress-point-fields 100000")).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let first = drive(&spec, &ledger, Some(1), &gate).unwrap();
    let preparing = load(&ledger, &first.pointer).unwrap();
    let second = resume(&ledger, &preparing, Some(1), &gate).unwrap();
    let candidate = load(&ledger, &second.pointer).unwrap();
    let third = resume(&ledger, &candidate, Some(1), &gate).unwrap();
    let promoted = load(&ledger, &third.pointer).unwrap();
    assert_eq!(phase(&promoted), "promoted");
    let promoted_report = report(&ledger, &promoted);
    let assessment = field(&promoted_report, "grid_assessment");
    assert_eq!(assessment.str_field("status"), Some("budget-exhausted"), "{:?}", assessment.str_field("reason"));
    assert!(
        assessment.str_field("reason").is_some_and(|reason| reason.contains("sdf3-affine-stress-point-budget")),
        "{:?}",
        assessment.str_field("reason")
    );
}

#[test]
fn g3_g5_affine_grid_refinement_resumes_candidate_without_replaying_or_installing_it_early() {
    let spec = spec::parse(&source()).unwrap();
    assert_eq!(spec.updates, 1);
    assert_eq!(spec.affine_grid_refinement.unwrap().max_updates, 2);
    assert_eq!(
        spec.affine_grid_assessment.unwrap().compliance_tolerance,
        0.0
    );
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();

    // The first invocation spends its single update on the source. It must
    // durably stop before constructing or optimizing the new grid.
    let first = drive(&spec, &ledger, Some(1), &gate).unwrap();
    let preparing = load(&ledger, &first.pointer).unwrap();
    assert_eq!(phase(&preparing), "preparing");
    assert_eq!(
        integer(&preparing.value, "iterations_completed").unwrap(),
        1
    );
    let original = source_receipt(&ledger, &preparing);
    assert_source_selected(&ledger, &original, &preparing);
    // The new artifact binds the Preparing context, while the underlying
    // source optimizer state and retained physical result are unchanged.
    let preparing_checkpoint =
        document(&linked(&ledger, &preparing.value, "checkpoint", resume::KIND).unwrap());
    let source_checkpoint =
        document(&linked(&ledger, &original.value, "checkpoint", resume::KIND).unwrap());
    assert_native_state_equal(&preparing_checkpoint, &source_checkpoint);
    for key in ["design", "iterations", "work", "restoration_evaluations"] {
        assert_eq!(
            field(&preparing_checkpoint, key),
            field(&source_checkpoint, key),
            "source {key}"
        );
    }
    for key in [
        "iterations",
        "evaluations",
        "multiplier_updates",
        "rejected_trials",
        "restoration_evaluations",
    ] {
        assert_eq!(
            field(field(&preparing.value, "optimizer_work"), key),
            field(field(&original.value, "optimizer_work"), key),
            "source {key}"
        );
    }
    let idle = resume(&ledger, &preparing, Some(0), &gate).unwrap();
    assert_eq!(idle.pointer, first.pointer);
    assert_eq!(idle.receipt, preparing.bytes);

    let second = resume(&ledger, &preparing, Some(1), &gate).unwrap();
    let candidate = load(&ledger, &second.pointer).unwrap();
    assert_eq!(phase(&candidate), "candidate");
    assert_eq!(
        integer(&candidate.value, "iterations_completed").unwrap(),
        1
    );
    assert_eq!(integer(&candidate.value, "target_iterations").unwrap(), 2);
    assert_source_selected(&ledger, &original, &candidate);
    assert_charged_after(&spec, &preparing, &candidate);
    let candidate_report = report(&ledger, &candidate);
    let candidate_refinement = field(&candidate_report, "grid_refinement");
    let source_cells = integer(field(candidate_refinement, "source_grid"), "cells").unwrap();
    let candidate_cells = integer(field(candidate_refinement, "working_grid"), "cells").unwrap();
    assert!(candidate_cells > source_cells);
    assert_eq!(
        field(candidate_refinement, "fresh_recheck").str_field("status"),
        Some("pending")
    );
    let working =
        document(&linked(&ledger, &candidate.value, "working_design", "study-design").unwrap());
    assert_eq!(
        field(field(&working, "selected"), "cells")
            .as_array()
            .unwrap()
            .len(),
        candidate_cells
    );
    assert_ne!(
        candidate.value.str_field("working_design"),
        candidate.value.str_field("design")
    );

    // Precancellation cannot silently install the candidate or erase the
    // immutable source and independently retained working-grid checkpoint.
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        resume(&ledger, &candidate, Some(1), &cancelled)
            .unwrap_err()
            .exit,
        exit::CANCELLED
    );
    assert_eq!(
        load(&ledger, &second.pointer).unwrap().bytes,
        candidate.bytes
    );

    let third = resume(&ledger, &candidate, Some(1), &gate).unwrap();
    let promoted = load(&ledger, &third.pointer).unwrap();
    assert_eq!(phase(&promoted), "promoted");
    assert_eq!(integer(&promoted.value, "iterations_completed").unwrap(), 2);
    assert_eq!(
        promoted.value.str_field("design"),
        promoted.value.str_field("working_design")
    );
    assert_charged_after(&spec, &candidate, &promoted);
    let promoted_report = report(&ledger, &promoted);
    let refinement = field(&promoted_report, "grid_refinement");
    assert_eq!(integer(refinement, "selected_generation").unwrap(), 1);
    assert_eq!(integer(refinement, "attempt").unwrap(), 1);
    let recheck = field(refinement, "fresh_recheck");
    assert_eq!(recheck.str_field("status"), Some("checked"));
    assert_eq!(field(recheck, "feasible"), &JsonValue::Bool(true));
    assert_eq!(
        field(field(recheck, "measurement"), "feasible"),
        &JsonValue::Bool(true)
    );
    let assessment = field(&promoted_report, "grid_assessment");
    assert_eq!(assessment.str_field("status"), Some("checked"));
    assert_eq!(
        integer(field(assessment, "coarse_grid"), "cells").unwrap(),
        candidate_cells
    );
    assert!(integer(field(assessment, "fine_grid"), "cells").unwrap() > candidate_cells);
    assert_eq!(field(assessment, "passed"), &JsonValue::Bool(false));
    assert_ne!(promoted.value.str_field("status"), Some("completed"));

    // A run without interruption uses the same candidate initialization,
    // accepted prefix and physical fields. Only endpoint restoration may add
    // optimizer evaluations; the original gradient experiment is retained.
    let full = drive(&spec, &ledger, Some(3), &gate).unwrap();
    let full = load(&ledger, &full.pointer).unwrap();
    assert_eq!(phase(&full), "promoted");
    for (key, kind) in [
        ("design", "study-design"),
        ("working_design", "study-design"),
        ("iterations", "study-iterations"),
    ] {
        assert_eq!(
            linked(&ledger, &promoted.value, key, kind).unwrap(),
            linked(&ledger, &full.value, key, kind).unwrap(),
            "{key}"
        );
    }
    let full_report = report(&ledger, &full);
    assert_eq!(
        field(&promoted_report, "gradient_check"),
        field(&candidate_report, "gradient_check")
    );
    assert_eq!(
        field(&promoted_report, "gradient_check"),
        field(&full_report, "gradient_check")
    );
    for key in ["selected", "last_accepted", "stop"] {
        assert_eq!(
            field(&promoted_report, key),
            field(&full_report, key),
            "{key}"
        );
    }
    let split_work = field(&promoted.value, "optimizer_work");
    let full_work = field(&full.value, "optimizer_work");
    let extra = integer(split_work, "restoration_evaluations").unwrap();
    assert!((1..=2).contains(&extra));
    assert_eq!(integer(full_work, "restoration_evaluations").unwrap(), 0);
    assert_eq!(
        integer(split_work, "evaluations").unwrap(),
        integer(full_work, "evaluations").unwrap() + extra
    );
    for key in ["iterations", "rejected_trials", "multiplier_updates"] {
        assert_eq!(field(split_work, key), field(full_work, key), "{key}");
    }
    let split_checkpoint =
        document(&linked(&ledger, &promoted.value, "checkpoint", resume::KIND).unwrap());
    let full_checkpoint =
        document(&linked(&ledger, &full.value, "checkpoint", resume::KIND).unwrap());
    assert_native_state_equal(&split_checkpoint, &full_checkpoint);
    let split_spent = Spent::read(&promoted.value).unwrap();
    let full_spent = Spent::read(&full.value).unwrap();
    assert!(split_spent.linear.linear_solves > full_spent.linear.linear_solves);
    assert!(split_spent.geometry.points > full_spent.geometry.points);
    assert_eq!(
        field(refinement, "source_grid"),
        field(field(&candidate_report, "grid_refinement"), "source_grid")
    );
    let finished = resume(&ledger, &promoted, Some(1), &gate).unwrap();
    assert_eq!(finished.pointer, third.pointer);
    assert_eq!(finished.receipt, promoted.bytes);
    for saved in [&original, &preparing, &candidate] {
        assert_eq!(
            load(&ledger, &format!("study-{}", saved.hash.to_hex()))
                .unwrap()
                .bytes,
            saved.bytes
        );
    }
}

#[test]
fn g4_affine_grid_refinement_keeps_source_when_candidate_evaluations_run_out() {
    let text = source().replace(
        ":max-optimizer-evaluations 500",
        ":max-optimizer-evaluations 2",
    );
    let spec = spec::parse(&text).unwrap();
    assert_eq!(spec.updates, 1);
    assert_eq!(spec.affine_grid_refinement.unwrap().max_evaluations, 2);
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let first = drive(&spec, &ledger, Some(1), &gate).unwrap();
    let preparing = load(&ledger, &first.pointer).unwrap();
    assert_eq!(phase(&preparing), "preparing");
    let original = source_receipt(&ledger, &preparing);
    let second = resume(&ledger, &preparing, Some(2), &gate).unwrap();
    let stopped = load(&ledger, &second.pointer).unwrap();
    assert_eq!(phase(&stopped), "candidate");
    assert_eq!(stopped.value.str_field("status"), Some("budget-exhausted"));
    assert!(integer(&stopped.value, "iterations_completed").unwrap() < 2);
    assert_eq!(
        integer(field(&stopped.value, "optimizer_work"), "evaluations").unwrap(),
        2
    );
    assert_eq!(
        field(&stopped.value, "stop").str_field("reason"),
        Some("EvaluationLimit")
    );
    assert_source_selected(&ledger, &original, &stopped);
    assert_charged_after(&spec, &preparing, &stopped);
    let stopped_report = report(&ledger, &stopped);
    let refinement = field(&stopped_report, "grid_refinement");
    assert_eq!(
        field(refinement, "fresh_recheck").str_field("status"),
        Some("pending")
    );
    assert!(
        integer(field(refinement, "working_grid"), "cells").unwrap()
            > integer(field(refinement, "source_grid"), "cells").unwrap()
    );

    // No fresh allowance may be substituted for the candidate's authored cap:
    // even its retained endpoint requires a charged native reconstruction.
    let error = resume(&ledger, &stopped, Some(100), &gate).unwrap_err();
    assert_eq!(error.code, "cli-study-sdf3-resume-budget");
    assert_eq!(error.exit, exit::BUDGET);
    assert_eq!(load(&ledger, &second.pointer).unwrap().bytes, stopped.bytes);
    assert_eq!(
        load(&ledger, &first.pointer).unwrap().bytes,
        preparing.bytes
    );
    assert_source_selected(&ledger, &original, &stopped);
}

#[test]
fn g0_affine_grid_refinement_requires_assessment_and_room_for_its_final_probe() {
    let spec = spec::parse(FIXTURE).unwrap();
    assert_eq!(spec.level, 1);
    assert_eq!(spec.max_level, 3);
    assert_eq!(spec.leaves, 512);
    let (head, tail) = FIXTURE.split_once("  (affine-grid-assessment").unwrap();
    let (assessment_tail, refinement_tail) = tail.split_once("  (affine-grid-refinement").unwrap();
    let assessment = format!("  (affine-grid-assessment{assessment_tail}");
    let refinement = format!("  (affine-grid-refinement{refinement_tail}");
    let refinement = refinement.trim_end().strip_suffix(')').unwrap();
    for (reason, source) in [
        (
            "missing original assessment",
            format!("{head}{refinement})\n"),
        ),
        (
            "no final-probe level",
            FIXTURE.replace(":maximum-level 3", ":maximum-level 2"),
        ),
        (
            "no final-probe leaf quota",
            FIXTURE.replace(":maximum-leaves 512", ":maximum-leaves 511"),
        ),
        (
            "duplicate refinement",
            format!("{head}{assessment}{refinement}\n{refinement})\n"),
        ),
        (
            "refinement before assessment",
            format!("{head}{refinement}\n{assessment})\n"),
        ),
        (
            "zero candidate updates",
            FIXTURE.replace(":max-updates 20", ":max-updates 0"),
        ),
        (
            "no candidate trial allowance",
            FIXTURE.replace(
                ":max-optimizer-evaluations 500",
                ":max-optimizer-evaluations 1",
            ),
        ),
    ] {
        assert!(spec::parse(&source).is_err(), "admitted {reason}");
    }
}
