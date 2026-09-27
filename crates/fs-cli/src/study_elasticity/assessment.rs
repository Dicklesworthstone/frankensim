//! Explicit final-design DWR assessment. The two bounded solves run only after
//! the last requested update has been durably retained. A stop at either phase
//! boundary publishes the same accepted geometry with assessment pending.

use super::*;
use fs_topols::ComplianceDwrAssessment;

pub(super) enum FinalAssessment {
    Estimated(ComplianceDwrAssessment),
    Refused(String),
}

const SCOPE: &str = "Estimated compliance goal error on the unchanged final bilinear level set. Signed DWR residual and absolute cell-indicator sum are not certified continuum-error bounds; the latter is marking mass, not an interval radius. Recomputed Euclidean solver residuals establish discrete solve accuracy only. Two bounded solves and residual integration are indivisible between cancellation checks. Each interrupted attempt consumes wall time and may be repeated on resume; a completed assessment is reused. Memory is an admitted allowance, not measured peak RSS. No physical validation or optimum is certified.";

pub(super) fn parse(root: &Node) -> Result<bool> {
    let Some(fields) = list(root, "study root")?.iter().find_map(|node| {
        let NodeKind::List(fields) = &node.kind else { return None };
        fields.first().is_some_and(|node| matches!(&node.kind,
            NodeKind::Symbol(name) if name == "assessment")).then_some(fields)
    }) else { return Ok(false) };
    if !matches!(&field(fields, "type")?.kind, NodeKind::Symbol(kind) if kind == "elasticity-dwr")
        || integer_node(field(fields, "max-solves-per-attempt")?, "assessment.max-solves-per-attempt")? != 2
    {
        return Err(fail("cli-study-elasticity-assessment",
            "final assessment requires (assessment :type elasticity-dwr :max-solves-per-attempt 2)"));
    }
    Ok(true)
}

pub(super) fn validate(spec: &ElasticitySpec) -> Result<()> {
    if spec.final_dwr && (spec.projected.is_some() || spec.memory_bytes < 256 * 1024 * 1024) {
        return Err(fail("cli-study-elasticity-assessment",
            "final elasticity DWR currently requires the plain single-load mode and at least 256 MiB admitted memory for coarse and enriched states"));
    }
    Ok(())
}

pub(super) fn require_same_design(assessment: &ComplianceDwrAssessment, report: &OptimizeReport) -> Result<()> {
    if report.snapshots.last() != Some(&assessment.snapshot)
        || report.compliance.last().map(|value| value.to_bits()) != Some(assessment.estimate.j_primal.to_bits())
        || report.volume.last().map(|value| value.to_bits()) != Some(assessment.volume.to_bits())
    {
        return Err(fail("cli-study-elasticity-dwr-design",
            "DWR coarse solve does not reproduce the accepted geometry, compliance and material area"));
    }
    Ok(())
}

pub(super) fn json(assessment: Option<&FinalAssessment>, requested: bool) -> String {
    let Some(a) = assessment else {
        return if requested {
            ",\"goal_error_assessment\":{\"status\":\"pending\",\"method\":\"elasticity-compliance-dwr\",\"max_solves_per_attempt\":2}".into()
        } else { String::new() };
    };
    let a = match a {
        FinalAssessment::Estimated(a) => a,
        FinalAssessment::Refused(reason) => return format!(
            ",\"goal_error_assessment\":{{\"status\":\"refused\",\"method\":\"elasticity-compliance-dwr\",\"reason\":{},\"max_solves_per_attempt\":2}}",
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
        assert!(super::super::parse(&with_assessment(projected)).is_err());
    }
}
