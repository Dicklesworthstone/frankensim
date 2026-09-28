//! Explicit final-design DWR assessment. The two bounded solves run only after
//! the actual optimizer endpoint has been durably retained. This includes a
//! feasible projected baseline when the search accepts no update. A stop at
//! either phase boundary retains the same geometry with assessment pending.

use super::*;
use fs_topols::{ComplianceDwrAssessment, ComplianceDwrStage};
use std::ops::ControlFlow;

// This file is itself included through #[path], so a bare `mod multi_load;`
// would resolve beside it (study_elasticity/multi_load.rs), not in assessment/.
#[path = "assessment/multi_load.rs"]
mod multi_load;
pub(super) use multi_load::run as run_multi_load;

pub(super) enum FinalAssessment {
    Estimated(ComplianceDwrAssessment),
    Weighted(fs_topols::WeightedComplianceDwrAssessment),
    Refused { reason: String, max_solves: usize },
}

const SCOPE: &str = "Estimated compliance goal error on the unchanged final bilinear level set. Signed DWR residual and absolute cell-indicator sum are not certified continuum-error bounds; the latter is marking mass, not an interval radius. Recomputed Euclidean solver residuals establish discrete solve accuracy only. Two bounded solves and residual integration are indivisible between cancellation checks. Each interrupted attempt consumes wall time and may be repeated on resume; a completed assessment is reused. Memory is an admitted allowance, not measured peak RSS. No physical validation or optimum is certified.";

fn solve_limit(projected: Option<&continuation::ProjectedControls>) -> Result<usize> {
    match projected {
        Some(continuation::ProjectedControls::Stress(policy)) => policy.dwr_case_count()
            .map(|count| 2 * count).ok_or_else(|| fail("cli-study-elasticity-assessment",
                "stress-mode DWR requires an explicit weighted-sum load family; worst-case objectives are not differentiable compliance sums")),
        _ => Ok(2),
    }
}

pub(super) fn max_solves(spec: &ElasticitySpec) -> usize {
    solve_limit(spec.projected.as_ref()).expect("admitted assessment mode")
}

pub(super) fn parse(root: &Node, projected: Option<&continuation::ProjectedControls>) -> Result<bool> {
    let Some(fields) = list(root, "study root")?.iter().find_map(|node| {
        let NodeKind::List(fields) = &node.kind else { return None };
        fields.first().is_some_and(|node| matches!(&node.kind,
            NodeKind::Symbol(name) if name == "assessment")).then_some(fields)
    }) else { return Ok(false) };
    let max_solves = solve_limit(projected)?;
    if !matches!(&field(fields, "type")?.kind, NodeKind::Symbol(kind) if kind == "elasticity-dwr")
        || integer_node(field(fields, "max-solves-per-attempt")?, "assessment.max-solves-per-attempt")? != max_solves
    {
        return Err(fail("cli-study-elasticity-assessment",
            format!("final assessment requires (assessment :type elasticity-dwr :max-solves-per-attempt {max_solves}); each independent load requires two solves")));
    }
    Ok(true)
}

pub(super) fn validate(spec: &ElasticitySpec) -> Result<()> {
    if spec.final_dwr && spec.memory_bytes < 256 * 1024 * 1024
    {
        return Err(fail("cli-study-elasticity-assessment",
            "final elasticity DWR requires at least 256 MiB admitted memory; independent cases are assessed sequentially, retaining indicators but not displacement states"));
    }
    Ok(())
}

/// Independently measured mechanics of the geometry being assessed. A feasible
/// projected baseline is valid even when the optimizer accepted no updates.
pub(super) struct Design {
    pub(super) snapshot: u64,
    pub(super) compliance: f64,
    pub(super) volume: f64,
}

impl Design {
    pub(super) fn from_report(report: &OptimizeReport) -> Result<Self> {
        match (report.snapshots.last(), report.compliance.last(), report.volume.last()) {
            (Some(&snapshot), Some(&compliance), Some(&volume)) => Ok(Self { snapshot, compliance, volume }),
            _ => Err(fail("cli-study-elasticity-dwr-design", "final assessment requires an evaluated accepted design")),
        }
    }
}

/// Shared numerical assessment for plain and same-area projected endpoints.
/// The caller durably retains its optimizer terminal before entering here.
pub(super) fn run(
    spec: &ElasticitySpec,
    phi: &GridSdf,
    expected: Design,
    terminal: &'static str,
    mut control: impl FnMut(ComplianceDwrStage) -> ControlFlow<&'static str>,
) -> (&'static str, Option<FinalAssessment>) {
    let assessed = fs_topols::assess_compliance_dwr_controlled(
        phi, fixture(spec), settings(spec, spec.steps), &mut control,
    ).map_err(|error| error.to_string()).and_then(|outcome| match outcome {
        ControlFlow::Continue(assessment) => {
            if expected.snapshot != assessment.snapshot
                || expected.compliance.to_bits() != assessment.estimate.j_primal.to_bits()
                || expected.volume.to_bits() != assessment.volume.to_bits()
            {
                Err("DWR coarse solve does not reproduce the accepted geometry, compliance and material area".into())
            } else {
                Ok(ControlFlow::Continue(assessment))
            }
        }
        ControlFlow::Break(stop) => Ok(ControlFlow::Break(stop)),
    });
    match assessed {
        Ok(ControlFlow::Continue(assessment)) => (terminal, Some(FinalAssessment::Estimated(assessment))),
        Ok(ControlFlow::Break(stop)) => (stop, None),
        Err(reason) => {
            let status = match control(ComplianceDwrStage::BeforePublish) {
                ControlFlow::Break(stop) => stop,
                ControlFlow::Continue(()) => "numerical-failure",
            };
            (status, Some(FinalAssessment::Refused { reason, max_solves: 2 }))
        }
    }
}

pub(super) fn json(assessment: Option<&FinalAssessment>, requested_solves: Option<usize>) -> String {
    let Some(a) = assessment else {
        return if let Some(max_solves) = requested_solves {
            format!(",\"goal_error_assessment\":{{\"status\":\"pending\",\"method\":\"elasticity-compliance-dwr\",\"max_solves_per_attempt\":{max_solves}}}")
        } else { String::new() };
    };
    let a = match a {
        FinalAssessment::Estimated(a) => a,
        FinalAssessment::Weighted(a) => return multi_load::json(a),
        FinalAssessment::Refused { reason, max_solves } => return format!(
            ",\"goal_error_assessment\":{{\"status\":\"refused\",\"method\":\"elasticity-compliance-dwr\",\"reason\":{},\"max_solves_per_attempt\":{max_solves}}}",
            quoted(reason),
        ),
    };
    let e = &a.estimate;
    format!(
        ",\"goal_error_assessment\":{{\"status\":\"estimated\",\"method\":\"elasticity-compliance-dwr\",\"snapshot\":\"{:#018x}\",\"material_area_m2\":{:.17e},\"coarse_level\":{},\"enriched_level\":{},\"eta_signed_j\":{:.17e},\"absolute_indicator_sum_j\":{:.17e},\"coarse_compliance_j\":{:.17e},\"enriched_compliance_j\":{:.17e},\"coarse_dofs\":{},\"enriched_dofs\":{},\"residual_terms_j\":{{\"bulk\":{:.17e},\"nitsche\":{:.17e},\"outer_traction\":{:.17e},\"ghost\":{:.17e}}},\"ghost_method\":{:?},\"solver\":{{\"relative_residual_kind\":\"recomputed-euclidean\",\"coarse_relative_residual\":{:.17e},\"enriched_relative_residual\":{:.17e},\"coarse_iterations\":{},\"enriched_iterations\":{},\"solves\":2,\"max_solves_per_attempt\":2,\"max_iterations_per_solve\":60000,\"relative_tolerance\":1e-12}},\"no_claim\":{}}}",
        a.snapshot, a.volume, a.level, a.level + 1, e.eta_signed, e.eta_abs,
        e.j_primal, e.j_enriched, e.dofs, e.enriched_dofs,
        e.terms.bulk, e.terms.nitsche, e.terms.outer_traction, e.terms.ghost,
        e.ghost_method.as_str(), e.primal_relative_residual, e.enriched_relative_residual,
        e.primal_iterations, e.enriched_iterations, quoted(SCOPE),
    )
}

pub(super) fn html(assessment: Option<&FinalAssessment>, requested: bool) -> String {
    let Some(a) = assessment else {
        return if requested { "<p>Final compliance goal-error assessment: pending.</p>".into() }
            else { String::new() };
    };
    if let FinalAssessment::Weighted(a) = a { return multi_load::html(a); }
    let FinalAssessment::Estimated(a) = a else {
        return "<p>Final compliance goal-error assessment: refused. The retained JSON reports the numerical reason; accepted geometry is unchanged.</p>".into();
    };
    let e = &a.estimate;
    format!("<h2>Final compliance goal-error estimate</h2><p>Signed DWR estimate {:.6e} J; absolute cell-indicator sum {:.6e} J. Coarse compliance {:.6e} J (level {}), enriched compliance {:.6e} J (level {}). Recomputed Euclidean relative residuals: {:.3e} and {:.3e}.</p><p>{SCOPE}</p>",
        e.eta_signed, e.eta_abs, e.j_primal, a.level, e.j_enriched, a.level + 1,
        e.primal_relative_residual, e.enriched_relative_residual)
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/marquee/bracket-2d.fsim"));
    const REQUEST: &str = "  (assessment :type elasticity-dwr :max-solves-per-attempt 2)\n";

    fn with_assessment(source: &str) -> String {
        format!("{}\n{REQUEST})\n", source.trim_end().strip_suffix(')').unwrap())
    }

    #[test]
    fn final_assessment_requires_explicit_bounded_work_and_supported_mode() {
        let source = with_assessment(FIXTURE);
        let spec = super::super::parse(&source).unwrap();
        assert!(spec.final_dwr);
        assert_eq!(super::super::parse(&spec.canonical).unwrap().id, spec.id);
        assert!(!super::super::parse(FIXTURE).unwrap().final_dwr);
        assert!(super::super::parse(&source.replace(":max-solves-per-attempt 2",
            ":max-solves-per-attempt 1")).is_err());
        assert!(super::super::parse(&source.replace(":memory 536870912 B",
            ":memory 134217728 B")).is_err());
        let projected = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
            "/../../examples/marquee/bracket-projected-volume-2d.fsim"));
        assert!(super::super::parse(&with_assessment(projected)).unwrap().final_dwr);
        let stress = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
            "/../../examples/marquee/bracket-projected-stress-2d.fsim"));
        assert!(super::super::parse(&with_assessment(stress)).is_err());
    }
}
