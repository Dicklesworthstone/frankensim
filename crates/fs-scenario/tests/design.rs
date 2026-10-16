//! G0 battery for design variables, constraints, objectives, and failure modes.
//!
//! The load-bearing tests are the refusals. An optimization request with no
//! manufacturing constraints must fail closed, because a design change nobody
//! can build is not an optimum; a robustness posture must be unrepresentable
//! as "unstated" rather than merely discouraged; and every declaration must
//! carry provenance, mirroring the crate's existing "no unsourced tolerance"
//! rule. The sidecar tests pin that the design block is bound to a specific
//! scenario and that its identity moves when any declared field moves.

use fs_qty::{Dims, QtyAny};
use fs_scenario::design::{
    DEFAULT_DESIGN_BUDGET, DesignBlock, DesignError, DesignSeverity, DesignSource, DesignVariable,
    DesignVariableKind, FailureMode, ManufacturingConstraint, ManufacturingConstraintKind,
    ObjectiveSense, OptimizationObjective, OptimizationRefusal, RobustnessPosture,
    ScenarioDesignExtension, VariableDomain, design_identity,
};
use fs_scenario::entity::{EntityCatalog, EntityDeclaration, EntityId, EntityRef, KindExpectation};
use fs_scenario::scenario::{Environment, Scenario};

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);

fn metres(value: f64) -> QtyAny {
    QtyAny::new(value, LENGTH)
}

fn kelvin(value: f64) -> QtyAny {
    QtyAny::new(value, TEMPERATURE)
}

fn standard(clause: &str) -> DesignSource {
    DesignSource::Standard {
        clause: clause.to_string(),
    }
}

struct Fixture {
    catalog: EntityCatalog,
    plate: EntityId,
    fin: EntityId,
}

fn fixture() -> Fixture {
    let mut catalog = EntityCatalog::new();
    let assembly = catalog
        .declare(EntityDeclaration::assembly("cold-plate-stack"))
        .expect("assembly");
    let plate = catalog
        .declare(EntityDeclaration::part(assembly, "cold-plate"))
        .expect("plate");
    let fin = catalog
        .declare(EntityDeclaration::part(assembly, "fin-array"))
        .expect("fin");
    Fixture {
        catalog,
        plate,
        fin,
    }
}

fn part_ref(id: EntityId) -> EntityRef {
    EntityRef::new(id, KindExpectation::Domain)
}

/// A worked design block: what may change, what may be built, what is being
/// optimized, and what counts as failure.
fn reference_block(fx: &Fixture) -> DesignBlock {
    let mut block = DesignBlock::new();
    block
        .declare_variable(DesignVariable::new(
            "fin-thickness",
            part_ref(fx.fin),
            DesignVariableKind::GeometryParameter,
            VariableDomain::Continuous {
                lo: metres(0.000_5),
                hi: metres(0.003),
            },
            standard("internal fin parameterization rev C"),
        ))
        .expect("variable admits");
    block
        .declare_constraint(ManufacturingConstraint::new(
            part_ref(fx.fin),
            ManufacturingConstraintKind::MinimumFeatureSpacing,
            metres(0.001),
            DesignSource::ProcessCapability {
                supplier: "extrusion-vendor-a".to_string(),
                document: "cap-2026-03 rev B".to_string(),
            },
        ))
        .expect("constraint admits");
    block
        .declare_objective(OptimizationObjective::new(
            "junction-temperature",
            ObjectiveSense::Minimize,
            RobustnessPosture::ConditionalValueAtRisk { beta: 0.95 },
            DesignSource::InternalPolicy {
                policy: "thermal-design-policy-2026".to_string(),
            },
        ))
        .expect("objective admits");
    block
        .declare_failure_mode(FailureMode::new(
            "thermal-runaway",
            "power-die",
            DesignSeverity::SafetyCritical,
            kelvin(423.15),
            DesignSource::Datasheet {
                part: "SiC-1200-40".to_string(),
                revision: "rev D".to_string(),
            },
        ))
        .expect("failure mode admits");
    block
}

/// The reference block is clean, and every declaration is reachable.
#[test]
fn the_reference_design_block_validates() {
    let fx = fixture();
    let block = reference_block(&fx);
    let findings = block.validate(&fx.catalog);
    assert!(
        findings.is_empty(),
        "reference block should be clean: {findings:?}"
    );

    assert_eq!(block.variables().len(), 1);
    assert_eq!(block.constraints().len(), 1);
    assert_eq!(block.objectives().len(), 1);
    assert_eq!(block.failure_modes().len(), 1);
    // Rows are assigned by the block, not the caller.
    assert_eq!(block.variables()[0].row(), 0);
    assert_eq!(block.failure_modes()[0].row(), 0);
    assert!(block.failure_mode("thermal-runaway", "power-die").is_some());
    assert!(block.failure_mode("thermal-runaway", "connector").is_none());
}

/// THE central refusal: an optimization request without manufacturing
/// constraints fails closed. A design change nobody can build is not an
/// optimum, so the absence of constraints is an unanswered question rather
/// than "none apply".
#[test]
fn an_optimization_request_without_manufacturing_constraints_refuses() {
    let fx = fixture();

    // Everything except constraints.
    let mut block = DesignBlock::new();
    block
        .declare_variable(DesignVariable::new(
            "fin-thickness",
            part_ref(fx.fin),
            DesignVariableKind::GeometryParameter,
            VariableDomain::Continuous {
                lo: metres(0.000_5),
                hi: metres(0.003),
            },
            standard("clause"),
        ))
        .expect("variable admits");
    block
        .declare_objective(OptimizationObjective::new(
            "junction-temperature",
            ObjectiveSense::Minimize,
            RobustnessPosture::Nominal,
            standard("clause"),
        ))
        .expect("objective admits");

    let refusals = block
        .optimization_request()
        .expect_err("an unconstrained request must refuse");
    assert_eq!(
        refusals,
        vec![OptimizationRefusal::NoManufacturingConstraint]
    );
    let refusal = &refusals[0];
    assert_eq!(refusal.code(), "design-optimization-unconstrained");
    // The message must be actionable, and must say why rather than only what.
    assert!(refusal.to_string().contains("nobody can build"));
    assert!(!refusal.fix().trim().is_empty());
    assert!(refusal.fix().contains("DesignSource"));

    // The reference block, which declares constraints, is admitted.
    let full = reference_block(&fx);
    let request = full
        .optimization_request()
        .expect("complete request admits");
    assert_eq!(request.objectives().len(), 1);
    assert_eq!(request.variables().len(), 1);
    assert_eq!(request.constraints().len(), 1);
}

/// Refusal is total: an empty block reports every missing piece at once, so a
/// caller repairs the whole set rather than one problem per attempt.
#[test]
fn refusals_report_every_missing_piece_at_once() {
    let block = DesignBlock::new();
    let refusals = block
        .optimization_request()
        .expect_err("an empty block cannot optimize");
    assert_eq!(refusals.len(), 3, "{refusals:?}");
    assert!(refusals.contains(&OptimizationRefusal::NoObjective));
    assert!(refusals.contains(&OptimizationRefusal::NoDesignVariable));
    assert!(refusals.contains(&OptimizationRefusal::NoManufacturingConstraint));

    let mut codes: Vec<&str> = refusals.iter().map(OptimizationRefusal::code).collect();
    codes.sort_unstable();
    codes.dedup();
    assert_eq!(codes.len(), 3, "codes are distinct");
    for refusal in &refusals {
        assert!(!refusal.fix().trim().is_empty(), "{refusal:?} has no fix");
    }
}

/// A robustness posture is a required constructor argument, so "unstated" is
/// unrepresentable rather than defaulted. Nominal is a DECLARED choice that
/// makes no robustness claim — not the absence of one.
#[test]
fn robustness_posture_is_declared_never_inferred() {
    let nominal = OptimizationObjective::new(
        "t-max",
        ObjectiveSense::Minimize,
        RobustnessPosture::Nominal,
        standard("clause"),
    );
    assert_eq!(nominal.posture(), RobustnessPosture::Nominal);
    assert!(
        !nominal.posture().is_robust(),
        "nominal declares no robustness claim"
    );

    let robust = OptimizationObjective::new(
        "t-max",
        ObjectiveSense::Minimize,
        RobustnessPosture::ConditionalValueAtRisk { beta: 0.9 },
        standard("clause"),
    );
    assert!(robust.posture().is_robust());
    assert_eq!(robust.posture().label(), "cvar");

    // A request is only "fully robust" when EVERY objective makes a claim.
    let fx = fixture();
    let mut block = reference_block(&fx);
    let request = block.optimization_request().expect("admits");
    assert!(request.fully_robust());
    block
        .declare_objective(nominal)
        .expect("second objective admits");
    let request = block.optimization_request().expect("admits");
    assert!(
        !request.fully_robust(),
        "one nominal objective removes the robustness claim for the request"
    );
}

/// Every declaration carries provenance, and an empty source field is named.
#[test]
fn there_is_no_unsourced_declaration() {
    let fx = fixture();
    let mut block = DesignBlock::new();
    block
        .declare_constraint(ManufacturingConstraint::new(
            part_ref(fx.plate),
            ManufacturingConstraintKind::MinimumWallThickness,
            metres(0.002),
            DesignSource::ProcessCapability {
                supplier: "   ".to_string(),
                document: String::new(),
            },
        ))
        .expect("admits structurally");

    let findings = block.validate(&fx.catalog);
    let empties: Vec<&str> = findings
        .iter()
        .filter(|violation| violation.code == "design-source-empty")
        .map(|violation| violation.what.as_str())
        .collect();
    assert_eq!(
        empties.len(),
        2,
        "both source fields are reported: {findings:?}"
    );
    assert!(empties.iter().any(|what| what.contains("supplier")));
    assert!(empties.iter().any(|what| what.contains("document")));
    for violation in &findings {
        assert!(!violation.fix.trim().is_empty(), "{violation:?} has no fix");
    }

    // Every source variant labels itself.
    for source in [
        standard("c"),
        DesignSource::ProcessCapability {
            supplier: "s".to_string(),
            document: "d".to_string(),
        },
        DesignSource::Datasheet {
            part: "p".to_string(),
            revision: "r".to_string(),
        },
        DesignSource::InternalPolicy {
            policy: "p".to_string(),
        },
        DesignSource::Assumed {
            rationale: "r".to_string(),
        },
    ] {
        assert!(!source.label().is_empty());
    }
}

/// Bounds carry units, are ordered, and are finite. A unit mismatch is a named
/// violation, never a silent reinterpretation.
#[test]
fn variable_bounds_are_dimensioned_ordered_and_finite() {
    let fx = fixture();
    let mut block = DesignBlock::new();
    // Dimensionless bounds on a geometry parameter, reversed, one non-finite.
    block
        .declare_variable(DesignVariable::new(
            "bad",
            part_ref(fx.fin),
            DesignVariableKind::GeometryParameter,
            VariableDomain::Continuous {
                lo: QtyAny::new(5.0, Dims([0, 0, 0, 0, 0, 0])),
                hi: QtyAny::new(f64::NAN, Dims([0, 0, 0, 0, 0, 0])),
            },
            standard("clause"),
        ))
        .expect("admits structurally");

    let codes: Vec<&str> = block
        .validate(&fx.catalog)
        .iter()
        .map(|violation| violation.code)
        .collect();
    assert!(codes.contains(&"design-variable-dims"), "{codes:?}");
    assert!(codes.contains(&"design-variable-nonfinite"), "{codes:?}");

    // Reversed but well-dimensioned bounds are an ordering fault.
    let mut ordered = DesignBlock::new();
    ordered
        .declare_variable(DesignVariable::new(
            "reversed",
            part_ref(fx.fin),
            DesignVariableKind::GeometryParameter,
            VariableDomain::Continuous {
                lo: metres(0.01),
                hi: metres(0.001),
            },
            standard("clause"),
        ))
        .expect("admits");
    let codes: Vec<&str> = ordered
        .validate(&fx.catalog)
        .iter()
        .map(|violation| violation.code)
        .collect();
    assert_eq!(codes, vec!["design-variable-bound-order"], "{codes:?}");
}

/// A discrete variable needs a nonempty, duplicate-free choice set, and a
/// continuous/discrete mismatch is named rather than coerced.
#[test]
fn discrete_domains_are_checked_against_their_kind() {
    let fx = fixture();
    let mut block = DesignBlock::new();
    // A material selection with an EMPTY choice set.
    block
        .declare_variable(DesignVariable::new(
            "plate-material",
            part_ref(fx.plate),
            DesignVariableKind::MaterialSelection,
            VariableDomain::Discrete { choices: vec![] },
            standard("clause"),
        ))
        .expect("admits");
    // A geometry parameter given a discrete domain.
    block
        .declare_variable(DesignVariable::new(
            "fin-thickness",
            part_ref(fx.fin),
            DesignVariableKind::GeometryParameter,
            VariableDomain::Discrete {
                choices: vec![
                    fs_blake3::ContentHash([1; 32]),
                    fs_blake3::ContentHash([1; 32]),
                ],
            },
            standard("clause"),
        ))
        .expect("admits");

    let codes: Vec<&str> = block
        .validate(&fx.catalog)
        .iter()
        .map(|violation| violation.code)
        .collect();
    assert!(codes.contains(&"design-variable-domain-empty"), "{codes:?}");
    assert!(codes.contains(&"design-variable-domain-kind"), "{codes:?}");
    assert!(
        codes.contains(&"design-variable-choice-duplicate"),
        "a repeated citation is reported: {codes:?}"
    );
}

/// A CVaR tail level must be a real probability, and a dangling entity target
/// is caught rather than silently accepted.
#[test]
fn objectives_and_targets_fail_closed() {
    let fx = fixture();
    let mut block = DesignBlock::new();
    block
        .declare_objective(OptimizationObjective::new(
            "",
            ObjectiveSense::Maximize,
            RobustnessPosture::ConditionalValueAtRisk { beta: 1.5 },
            standard("clause"),
        ))
        .expect("admits");

    let codes: Vec<&str> = block
        .validate(&fx.catalog)
        .iter()
        .map(|violation| violation.code)
        .collect();
    assert!(codes.contains(&"design-objective-qoi-empty"), "{codes:?}");
    assert!(codes.contains(&"design-objective-cvar-beta"), "{codes:?}");

    // A target outside the catalog is dangling, not assumed live.
    let orphan = EntityDeclaration::assembly("never-declared").identity();
    let mut dangling = DesignBlock::new();
    dangling
        .declare_constraint(ManufacturingConstraint::new(
            part_ref(orphan),
            ManufacturingConstraintKind::MinimumWallThickness,
            metres(0.002),
            standard("clause"),
        ))
        .expect("admits");
    let codes: Vec<&str> = dangling
        .validate(&fx.catalog)
        .iter()
        .map(|violation| violation.code)
        .collect();
    assert!(
        codes.contains(&"design-constraint-dangling-subject"),
        "{codes:?}"
    );
}

/// Failure-mode severity orders consequences, and thresholds are checked.
#[test]
fn failure_modes_are_ordered_and_dimensioned() {
    assert!(DesignSeverity::ReliabilityDerating < DesignSeverity::DamageLimit);
    assert!(DesignSeverity::DamageLimit < DesignSeverity::SafetyCritical);

    let fx = fixture();
    let mut block = DesignBlock::new();
    block
        .declare_failure_mode(FailureMode::new(
            "",
            "",
            DesignSeverity::DamageLimit,
            metres(1.0),
            standard("clause"),
        ))
        .expect("admits");
    let codes: Vec<&str> = block
        .validate(&fx.catalog)
        .iter()
        .map(|violation| violation.code)
        .collect();
    assert!(
        codes.contains(&"design-failure-mode-name-empty"),
        "{codes:?}"
    );
    assert!(
        codes.contains(&"design-failure-mode-class-empty"),
        "{codes:?}"
    );
    assert!(
        codes.contains(&"design-failure-mode-dims"),
        "a length is not a thermal threshold: {codes:?}"
    );
}

/// Declarations are budgeted, and the bound is a refusal rather than a panic.
#[test]
fn declarations_are_budgeted() {
    let fx = fixture();
    let mut block = DesignBlock::with_budget(fs_scenario::design::DesignBudget {
        max_objectives: 1,
        ..DEFAULT_DESIGN_BUDGET
    });
    let objective = || {
        OptimizationObjective::new(
            "t-max",
            ObjectiveSense::Minimize,
            RobustnessPosture::Nominal,
            standard("clause"),
        )
    };
    assert_eq!(block.declare_objective(objective()).expect("first"), 0);
    let error = block
        .declare_objective(objective())
        .expect_err("second exceeds the budget");
    assert_eq!(
        error,
        DesignError::CapacityExceeded {
            resource: "objectives",
            requested: 2,
            limit: 1,
        }
    );
    assert_eq!(error.code(), "design-capacity-exceeded");
    assert!(block.validate(&fx.catalog).is_empty());
}

/// The sidecar binds to ONE scenario, and its identity moves when any declared
/// field moves. A sidecar that does not verify is not evidence about that
/// scenario.
#[test]
fn the_design_sidecar_is_bound_and_content_addressed() {
    let fx = fixture();
    let block = reference_block(&fx);
    let scenario = Scenario::new("cold-plate", 7, Environment::earth_lab());
    let other = Scenario::new("other-study", 9, Environment::earth_lab());

    let extension = ScenarioDesignExtension::new(&scenario, block.clone());
    assert!(extension.verifies_scenario(&scenario));
    assert!(
        !extension.verifies_scenario(&other),
        "a sidecar must not verify a scenario it was not built from"
    );
    assert_eq!(extension.block(), &block);

    // Rebuilding from the same inputs is identical.
    let again = ScenarioDesignExtension::new(&scenario, block.clone());
    assert_eq!(extension.identity(), again.identity());
    assert_eq!(
        extension.identity(),
        design_identity(extension.scenario_ir_hash(), &block)
    );

    // A different scenario yields a different identity for the same block.
    let moved = ScenarioDesignExtension::new(&other, block.clone());
    assert_ne!(extension.identity(), moved.identity());

    // Moving ANY declared field moves the identity.
    let mut altered = block.clone();
    altered
        .declare_failure_mode(FailureMode::new(
            "solder-fatigue",
            "bga",
            DesignSeverity::ReliabilityDerating,
            kelvin(398.15),
            standard("IPC-9701"),
        ))
        .expect("admits");
    assert_ne!(
        extension.identity(),
        ScenarioDesignExtension::new(&scenario, altered).identity()
    );

    // A changed source is a changed declaration, even at the same magnitude.
    let mut resourced = DesignBlock::new();
    resourced
        .declare_constraint(ManufacturingConstraint::new(
            part_ref(fx.plate),
            ManufacturingConstraintKind::MinimumWallThickness,
            metres(0.002),
            standard("clause-a"),
        ))
        .expect("admits");
    let mut other_source = DesignBlock::new();
    other_source
        .declare_constraint(ManufacturingConstraint::new(
            part_ref(fx.plate),
            ManufacturingConstraintKind::MinimumWallThickness,
            metres(0.002),
            standard("clause-b"),
        ))
        .expect("admits");
    assert_ne!(
        design_identity(extension.scenario_ir_hash(), &resourced),
        design_identity(extension.scenario_ir_hash(), &other_source),
        "provenance is part of the declaration, not decoration"
    );
}
