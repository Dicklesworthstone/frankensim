//! Weighted-sum compliance error, with independent load evidence. Sampled
//! stress feasibility remains the optimizer's separate admission decision.
use super::*;
use fs_topols::{SampledStressEvaluation, WeightedComplianceDwrAssessment, WeightedComplianceDwrStage};

const SCOPE: &str = "Estimated error in the declared weighted-sum compliance objective on one unchanged bilinear level set. Each traction is solved independently; weights are not probabilities and zero-weight cases retain their own evidence. Signed residual estimates and weighted absolute indicator mass are not certified continuum-error bounds. This is not a stress error estimator or a stress certificate. Each two-solve case is indivisible between cancellation checks; interrupted attempts may repeat and consume lifetime wall time. No optimization convergence or physical validation is claimed.";

pub(in super::super) fn run(
    spec: &ElasticitySpec, phi: &GridSdf, expected: Design,
    expected_cases: &[SampledStressEvaluation], terminal: &'static str,
    mut control: impl FnMut(WeightedComplianceDwrStage) -> ControlFlow<&'static str>,
) -> (&'static str, Option<FinalAssessment>) {
    let assessed = (|| {
        let Some(continuation::ProjectedControls::Stress(policy)) = &spec.projected
            else { return Err("missing weighted-sum load family".into()) };
        let cases = policy.dwr_cases(spec).map_err(|error| error.message)?;
        if cases.len() != expected_cases.len() || expected_cases.iter().any(|case|
            case.snapshot != expected.snapshot || case.volume.to_bits() != expected.volume.to_bits()) {
            return Err("retained independent cases do not describe one complete design".into());
        }
        match fs_topols::assess_weighted_compliance_dwr_controlled(
            phi, &cases, settings(spec, spec.steps), &mut control,
        ).map_err(|error| error.to_string())? {
            ControlFlow::Break(stop) => Ok(ControlFlow::Break(stop)),
            ControlFlow::Continue(a) => {
                if a.snapshot != expected.snapshot || a.volume.to_bits() != expected.volume.to_bits()
                    || a.coarse_compliance.to_bits() != expected.compliance.to_bits()
                    || a.cases.iter().zip(expected_cases).any(|(actual, expected)|
                        actual.estimate.j_primal.to_bits() != expected.compliance.to_bits()) {
                    return Err("DWR coarse case solves do not reproduce the retained geometry, independent compliances and material area".into());
                }
                Ok(ControlFlow::Continue(a))
            }
        }
    })();
    match assessed {
        Ok(ControlFlow::Continue(a)) => (terminal, Some(FinalAssessment::Weighted(a))),
        Ok(ControlFlow::Break(stop)) => (stop, None),
        Err(reason) => {
            let status = match control(WeightedComplianceDwrStage::BeforePublish) {
                ControlFlow::Break(stop) => stop,
                ControlFlow::Continue(()) => "numerical-failure",
            };
            (status, Some(FinalAssessment::Refused { reason, max_solves: max_solves(spec) }))
        }
    }
}

pub(super) fn json(a: &WeightedComplianceDwrAssessment) -> String {
    let cases = a.cases.iter().enumerate().map(|(index, case)| {
        let e = &case.estimate;
        let [start, end] = case.load.interval();
        let [x, y] = case.load.traction();
        format!(concat!("{{\"case\":{index},\"weight\":{:.17e},\"edge\":\"right\",",
            "\"band\":[{start:.17e},{end:.17e}],\"traction_pa\":[{x:.17e},{y:.17e}],",
            "\"coarse_compliance_j\":{:.17e},\"enriched_compliance_j\":{:.17e},",
            "\"eta_signed_j\":{:.17e},\"absolute_indicator_sum_j\":{:.17e},",
            "\"residual_terms_j\":{{\"bulk\":{:.17e},\"nitsche\":{:.17e},\"outer_traction\":{:.17e},\"ghost\":{:.17e}}},",
            "\"ghost_method\":{:?},\"coarse_dofs\":{},\"enriched_dofs\":{},",
            "\"coarse_relative_residual\":{:.17e},\"enriched_relative_residual\":{:.17e},",
            "\"coarse_iterations\":{},\"enriched_iterations\":{}}}"),
            case.load.weight(), e.j_primal, e.j_enriched, e.eta_signed, e.eta_abs,
            e.terms.bulk, e.terms.nitsche, e.terms.outer_traction, e.terms.ghost,
            e.ghost_method.as_str(), e.dofs, e.enriched_dofs, e.primal_relative_residual,
            e.enriched_relative_residual, e.primal_iterations, e.enriched_iterations,
            index = index, start = start, end = end, x = x, y = y)
    }).collect::<Vec<_>>().join(",");
    let solves = 2 * a.cases.len();
    format!(concat!(",\"goal_error_assessment\":{{\"status\":\"estimated\",",
        "\"method\":\"elasticity-compliance-dwr\",\"aggregate\":\"weighted-sum\",",
        "\"snapshot\":\"{:#018x}\",\"material_area_m2\":{:.17e},",
        "\"coarse_level\":{},\"enriched_level\":{},",
        "\"coarse_compliance_j\":{:.17e},\"enriched_compliance_j\":{:.17e},",
        "\"eta_signed_j\":{:.17e},\"absolute_indicator_sum_j\":{:.17e},",
        "\"residual_terms_j\":{{\"bulk\":{:.17e},\"nitsche\":{:.17e},\"outer_traction\":{:.17e},\"ghost\":{:.17e}}},",
        "\"cases\":[{cases}],\"solver\":{{\"relative_residual_kind\":\"recomputed-euclidean\",",
        "\"solves\":{solves},\"max_solves_per_attempt\":{solves},\"max_iterations_per_solve\":60000,",
        "\"relative_tolerance\":1e-12}},\"no_claim\":{}}}"),
        a.snapshot, a.volume, a.level, a.level + 1, a.coarse_compliance, a.enriched_compliance,
        a.eta_signed, a.absolute_indicator_sum, a.terms.bulk, a.terms.nitsche,
        a.terms.outer_traction, a.terms.ghost, quoted(SCOPE), cases = cases, solves = solves)
}

pub(super) fn html(a: &WeightedComplianceDwrAssessment) -> String {
    format!("<h2>Final weighted-sum compliance goal-error estimate</h2><p>{} independent load cases. Coarse compliance {:.6e} J (level {}); enriched compliance {:.6e} J (level {}). Signed DWR estimate {:.6e} J; weighted absolute cell-indicator mass {:.6e} J. Per-case loads, residual breakdown and actual solver residuals are retained in JSON.</p><p>{SCOPE}</p>",
        a.cases.len(), a.coarse_compliance, a.level, a.enriched_compliance, a.level + 1,
        a.eta_signed, a.absolute_indicator_sum)
}
