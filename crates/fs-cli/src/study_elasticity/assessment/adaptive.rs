//! Final-design local refinement through the existing adaptive CutFEM owner.
//! The optimizer endpoint is already durable. This module neither evolves the
//! material design nor substitutes adaptive-grid compliance into its history.
use super::*;
use fs_topols::adaptive_dwr::{AdaptiveDwrPolicy, AdaptiveDwrReport, AdaptiveDwrStage,
    AdaptiveDwrStop, assess_adaptive_dwr_controlled};
use fs_topols::{DesignBoxEdge, RobustLoadCase};

pub(in crate::study::elasticity) struct Assessment {
    report: AdaptiveDwrReport,
    policy: AdaptiveDwrPolicy,
    cases: Vec<RobustLoadCase>,
}

pub(super) fn parse(fields: &[Node]) -> Result<AdaptiveDwrPolicy> {
    let count = |key| integer_node(field(fields, key)?, key);
    let real = |key| number(field(fields, key)?, key);
    let policy = AdaptiveDwrPolicy {
        max_level: u32::try_from(count("max-level")?)
            .map_err(|_| fail("cli-study-adaptive-assessment", "max-level exceeds u32"))?,
        max_rounds: count("max-rounds")?, max_cells: count("max-cells")?,
        absolute_tolerance: real("absolute-tolerance-j")?,
        relative_tolerance: real("relative-tolerance")?, mark_fraction: real("mark-fraction")?,
    };
    policy.validate(1, 1).map_err(|e| fail("cli-study-adaptive-assessment", e.to_string()))?;
    Ok(policy)
}

pub(super) fn canonical(policy: AdaptiveDwrPolicy, max_solves: usize, out: &mut String) {
    let _ = writeln!(out, "  (assessment :type elasticity-adaptive-dwr");
    let _ = writeln!(out, "    :max-solves-per-attempt {max_solves}");
    let _ = writeln!(out, "    :max-level {}", policy.max_level);
    let _ = writeln!(out, "    :max-rounds {}", policy.max_rounds);
    let _ = writeln!(out, "    :max-cells {}", policy.max_cells);
    let _ = writeln!(out, "    :absolute-tolerance-j {}", canonical_float(policy.absolute_tolerance));
    let _ = writeln!(out, "    :relative-tolerance {}", canonical_float(policy.relative_tolerance));
    let _ = writeln!(out, "    :mark-fraction {})", canonical_float(policy.mark_fraction));
}

pub(super) fn run_single(spec: &ElasticitySpec, phi: &GridSdf, expected: Design,
    terminal: &'static str, mut control: impl FnMut(ComplianceDwrStage) -> ControlFlow<&'static str>)
    -> (&'static str, Option<FinalAssessment>)
{
    let fixture = fixture(spec);
    let case = RobustLoadCase::new(DesignBoxEdge::Right,
        0.5 - fixture.band, 0.5 + fixture.band, [0.0, -fixture.load], 1.0);
    let case = match case {
        Ok(case) => case,
        Err(e) => return ("numerical-failure", Some(FinalAssessment::Refused {
            reason: e.to_string(), max_solves: max_solves(spec),
        })),
    };
    let compliance = expected.compliance;
    run(spec, phi, expected, &[case], &[compliance], terminal, |stage| {
        // Keep existing native phase callbacks. They repeat for each adaptive
        // round/case; each before-work and after-work boundary remains polled.
        control(match stage {
            AdaptiveDwrStage::BeforeCase { .. } | AdaptiveDwrStage::BeforeRefine(_) =>
                ComplianceDwrStage::BeforeEstimate,
            AdaptiveDwrStage::RoundComplete(_) | AdaptiveDwrStage::BeforePublish =>
                ComplianceDwrStage::BeforePublish,
        })
    })
}

/// Assessment-only retry uses the identical retained endpoint, never optimizer
/// replay. A stop publishes no partial analysis. Per-attempt solve allowance and
/// cumulative wall charging retain the existing final-assessment semantics.
pub(super) fn run(spec: &ElasticitySpec, phi: &GridSdf, expected: Design,
    cases: &[RobustLoadCase], expected_cases: &[f64], terminal: &'static str,
    mut control: impl FnMut(AdaptiveDwrStage) -> ControlFlow<&'static str>)
    -> (&'static str, Option<FinalAssessment>)
{
    let assessed = (|| -> std::result::Result<ControlFlow<&'static str, Assessment>, String> {
        let policy = spec.adaptive_dwr.ok_or("missing adaptive assessment policy")?;
        if cases.is_empty() || cases.len() != expected_cases.len()
            || expected_cases.iter().any(|c| !c.is_finite() || *c < 0.0)
            || snapshot(phi) != expected.snapshot
        { return Err("adaptive assessment requires a complete retained design/load family".into()); }
        let report = match assess_adaptive_dwr_controlled(phi, cases, settings(spec, spec.steps),
            policy, &mut control).map_err(|e| e.to_string())? {
            ControlFlow::Break(stop) => return Ok(ControlFlow::Break(stop)),
            ControlFlow::Continue(report) => report,
        };
        let first = report.rows.first().ok_or("adaptive assessment has no initial solve")?;
        if report.snapshot != expected.snapshot || report.original_area.to_bits() != expected.volume.to_bits()
            || first.coarse_compliance.to_bits() != expected.compliance.to_bits()
            || report.initial_case_compliances.len() != expected_cases.len()
            || report.initial_case_compliances.iter().zip(expected_cases)
                .any(|(actual, expected)| actual.to_bits() != expected.to_bits())
            || report.solves > max_solves(spec)
        { return Err("adaptive initial solves do not reproduce the retained geometry, area and independent compliances".into()); }
        Ok(ControlFlow::Continue(Assessment { report, policy, cases: cases.to_vec() }))
    })();
    match assessed {
        Ok(ControlFlow::Continue(a)) => (terminal, Some(FinalAssessment::Adaptive(Box::new(a)))),
        Ok(ControlFlow::Break(stop)) => (stop, None),
        Err(reason) => {
            let status = match control(AdaptiveDwrStage::BeforePublish) {
                ControlFlow::Break(stop) => stop, ControlFlow::Continue(()) => "numerical-failure",
            };
            (status, Some(FinalAssessment::Refused { reason, max_solves: max_solves(spec) }))
        }
    }
}

const SCOPE: &str = "Estimated goal-driven local refinement of the analysis grid on the unchanged retained bilinear design. Adaptive compliance and area do not replace the original optimizer metrics. Tolerance refers to numerical signed-DWR corrections and observed two-grid changes, not a continuum bound, stress certificate, optimum or physical validation. The optimizer status and assessment stop are separate; completed optimization does not imply the requested assessment tolerance was met. Each case's two solves and residual integration are indivisible. Interrupted assessment retries start again on the original analysis grid under the declared per-attempt allowance and cumulative wall budget.";

pub(super) fn json(a: &Assessment) -> String {
    let r = &a.report;
    let p = a.policy;
    let rows = r.rows.iter().map(|row| format!(concat!(
        "{{\"cells\":{},\"max_cell_level\":{},\"coarse_compliance_j\":{:.17e},",
        "\"enriched_compliance_j\":{:.17e},\"eta_signed_j\":{:.17e},",
        "\"correction_magnitude_j\":{:.17e},\"observed_change_j\":{:.17e},",
        "\"tolerance_j\":{:.17e},\"coarse_dofs\":{},\"enriched_dofs\":{},",
        "\"relative_residual\":{:.17e}}}"), row.cells, row.max_cell_level,
        row.coarse_compliance, row.enriched_compliance, row.eta_signed, row.correction_magnitude,
        row.observed_change, row.tolerance, row.coarse_dofs, row.enriched_dofs, row.relative_residual))
        .collect::<Vec<_>>().join(",");
    let cells = r.grid.leaves().map(|(level, i, j)| format!("[{level},{i},{j}]"))
        .collect::<Vec<_>>().join(",");
    let cases = a.cases.iter().zip(&r.estimates).enumerate().map(|(i, (load, e))| {
        let [start, end] = load.interval(); let [x, y] = load.traction();
        let indicators = e.indicators.iter().map(|(&(level, x, y), eta)|
            format!("[{level},{x},{y},{eta:.17e}]")).collect::<Vec<_>>().join(",");
        format!(concat!("{{\"case\":{},\"weight\":{:.17e},\"edge\":\"right\",",
            "\"band\":[{:.17e},{:.17e}],\"traction_pa\":[{:.17e},{:.17e}],",
            "\"initial_compliance_j\":{:.17e},\"coarse_compliance_j\":{:.17e},",
            "\"enriched_compliance_j\":{:.17e},\"eta_signed_j\":{:.17e},",
            "\"residual_terms_j\":{{\"bulk\":{:.17e},\"nitsche\":{:.17e},",
            "\"outer_traction\":{:.17e},\"ghost\":{:.17e}}},\"ghost_method\":{},",
            "\"coarse_relative_residual\":{:.17e},\"enriched_relative_residual\":{:.17e},",
            "\"coarse_iterations\":{},\"enriched_iterations\":{},\"indicators\":[{}]}}"),
            i, load.weight(), start, end, x, y, r.initial_case_compliances[i], e.j_primal,
            e.j_enriched, e.eta_signed, e.terms.bulk, e.terms.nitsche, e.terms.outer_traction,
            e.terms.ghost, quoted(e.ghost_method.as_str()), e.primal_relative_residual,
            e.enriched_relative_residual, e.primal_iterations, e.enriched_iterations, indicators)
    }).collect::<Vec<_>>().join(",");
    format!(concat!(",\"goal_error_assessment\":{{\"status\":\"estimated\",",
        "\"method\":\"adaptive-elasticity-compliance-dwr\",\"aggregate\":\"weighted-sum\",",
        "\"stop\":{},\"tolerance_satisfied\":{},\"snapshot\":\"{:#018x}\",",
        "\"original_material_area_m2\":{:.17e},\"assessed_material_area_m2\":{:.17e},",
        "\"policy\":{{\"max_level\":{},\"max_rounds\":{},\"max_cells\":{},",
        "\"absolute_tolerance_j\":{:.17e},\"relative_tolerance\":{:.17e},\"mark_fraction\":{:.17e}}},",
        "\"solves\":{},\"max_solves_per_attempt\":{},\"rounds\":[{}],",
        "\"grid_cells\":[{}],\"cases\":[{}],\"no_claim\":{}}}"),
        quoted(r.stop.as_str()), r.stop == AdaptiveDwrStop::EstimatedTolerance, r.snapshot,
        r.original_area, r.assessed_area, p.max_level, p.max_rounds, p.max_cells,
        p.absolute_tolerance, p.relative_tolerance, p.mark_fraction, r.solves,
        p.max_solves(a.cases.len()).expect("admitted work"), rows, cells, cases, quoted(SCOPE))
}

pub(super) fn html(a: &Assessment) -> String {
    let mut out = format!("<h2>Adaptive compliance assessment</h2><p>Analysis stop: {}. Requested estimated tolerance met: {}. {} independent case(s), {} solves. The optimized material design is unchanged.</p><table><tr><th>Round</th><th>Cells</th><th>Coarse J</th><th>Enriched J</th><th>Correction magnitude J</th><th>Observed change J</th><th>Tolerance J</th></tr>",
        a.report.stop.as_str(), a.report.stop == AdaptiveDwrStop::EstimatedTolerance,
        a.cases.len(), a.report.solves);
    for (i, row) in a.report.rows.iter().enumerate() {
        let _ = write!(out, "<tr><td>{}</td><td>{}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td></tr>",
            i, row.cells, row.coarse_compliance, row.enriched_compliance,
            row.correction_magnitude, row.observed_change, row.tolerance);
    }
    let _ = write!(out, "</table><p>{SCOPE}</p>");
    out
}

#[cfg(test)]
#[path = "adaptive/tests.rs"]
mod tests;
