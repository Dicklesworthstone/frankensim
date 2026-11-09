//! A separate physical check of the selected design on two actual grids.
//! Native owners retain every load direction and every original constraint.
use super::super::grid_assessment::{Outcome as Assessment, Policy};
use super::*;
use fs_topopt::sdf3::robust::envelope::affine_load_space::{
    AffineLoadSpaceGrid3, refinement::AffineLoadSpaceRefinement3,
    stress::refinement::AffineStressChange3,
};

const SCOPE: &str = "Numerical two-grid assessment of the retained selected design. The fine grid inherits eroded physical occupancy and stiffness and dilated inventory from source parents; it does not refilter or optimize the design. Both grids independently maximize the complete affine load family and check the original compliance or volume cap, signed-response and sampled-strength requirements. For a compliance objective the reference compliance normalizes the objective and full-model changes; physical feasibility uses the dilated-volume cap. Mean displacement uses each grid's own reference cut-volume normalization on the same authored patch. Full-model change metrics include directions other than the compliance maximizers. Frozen-direction DWR terms and load-switch corrections are separate diagnostics. Stress aggregates use each grid's actual quadrature; matched-site changes use both complete field families on the same fine sites. Changing sample counts can change the unweighted stress aggregate. These finite-sample, floating-point observations are not continuum or interval error bounds. The optimizer checkpoint and selected coarse field artifact remain separate from this assessment.";

fn invalid(message: impl Into<String>) -> Failure {
    fail("cli-study-sdf3-grid-assessment-output", message)
}

pub(super) fn policy(p: Policy) -> String {
    format!(
        "{{\"max_transfer_terms\":{},\"compliance_tolerance\":{:.17e},\"response_tolerance\":{:.17e},\"stress_tolerance\":{:.17e},\"inventory_tolerance\":{:.17e}}}",
        p.max_transfer_terms,
        p.compliance_tolerance,
        p.response_tolerance,
        p.stress_tolerance,
        p.inventory_tolerance,
    )
}

fn grid(g: &AffineLoadSpaceGrid3, config: Config, volume: f64, joint: f64) -> String {
    let o = &g.oracle;
    let requirement = match config.goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit } => {
            let feasible = volume.is_finite() && volume >= 0.0 && volume <= volume_limit;
            format!(
                "\"compliance_limit_j\":null,\"compliance_feasible\":null,\"volume_limit_fraction\":{volume_limit:.17e},\"volume_feasible\":{feasible}"
            )
        }
        AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => {
            let feasible = g.maximum.is_finite()
                && g.worst_compliance.is_finite()
                && (0.0..=1.0).contains(&g.maximum)
                && g.worst_compliance >= 0.0
                && g.worst_compliance <= config.cap_j;
            format!(
                "\"compliance_limit_j\":{:.17e},\"compliance_feasible\":{feasible},\"volume_limit_fraction\":null,\"volume_feasible\":null",
                config.cap_j,
            )
        }
    };
    format!(
        "{{{requirement},\"compliance_reference_j\":{:.17e},\"dilated_volume_fraction\":{volume:.17e},\"worst_compliance_j\":{:.17e},\"maximum_ratio\":{:.17e},\"joint_constraint\":{joint:.17e},\"worst_coefficients\":{:?},\"normalized_gram_row_major\":{:?},\"normalized_linear\":{:?},\"normalized_offset\":{:.17e},\"reciprocity_defect\":{:.17e},\"witness_relative_defect\":{:.17e},\"witness_relative_residual\":{:.17e},\"oracle\":{{\"branch\":{:?},\"gradient_admissible\":{},\"value_without_offset\":{:.17e},\"dual_value_without_offset\":{:.17e},\"relative_duality_gap\":{:.17e},\"stationarity_residual\":{:.17e},\"feasibility_residual\":{:.17e},\"complementarity_residual\":{:.17e},\"eigenpair_residual\":{:.17e},\"iterations\":{}}}}}",
        config.cap_j,
        g.worst_compliance,
        g.maximum,
        g.worst_coefficients,
        g.gram,
        g.linear,
        g.offset,
        g.reciprocity_defect,
        g.witness_relative_defect,
        g.witness_relative_residual,
        format!("{:?}", o.branch),
        o.gradient_admissible(),
        o.value,
        o.dual_value,
        o.relative_duality_gap,
        o.stationarity_residual,
        o.feasibility_residual,
        o.complementarity_residual,
        o.eigenpair_residual,
        o.iterations,
    )
}

fn goal(r: &AffineLoadSpaceRefinement3, fine_maximizer: bool) -> String {
    let g = if fine_maximizer {
        &r.fine_witness
    } else {
        &r.coarse_witness
    };
    format!(
        "{{\"coarse_value_j\":{:.17e},\"fine_value_j\":{:.17e},\"dwr_j\":{:.17e},\"coarse_space_j\":{:.17e},\"goal_transfer_j\":{:.17e},\"algebraic_j\":{:.17e},\"correction_j\":{:.17e},\"identity_relative_defect\":{:.17e},\"field_residuals\":{:?},\"diagnostic_cells\":{}}}",
        g.coarse_value,
        g.fine_value,
        g.dwr,
        g.coarse_space,
        g.goal_transfer,
        g.algebraic,
        g.correction(),
        g.identity_relative_defect,
        g.field_residuals,
        g.cells.len(),
    )
}

fn volumes(values: &[Option<f64>]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|v| v.map_or_else(|| "null".into(), |v| format!("{v:.17e}")))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn responses(
    spec: &Spec,
    a: &Assessment,
    r: &AffineLoadSpaceRefinement3,
    p: Policy,
) -> Result<(String, String)> {
    let (family, response) = match (&spec.affine_responses, &r.responses) {
        (None, None) if a.coarse_mean_volumes.is_empty() && a.fine_mean_volumes.is_empty() => {
            return Ok(("null".into(), String::new()));
        }
        (Some(family), Some(response)) => (family, response),
        _ => {
            return Err(invalid(
                "assessment response family differs from its source",
            ));
        }
    };
    if response.normalized_model_changes.len() != family.limits.len()
        || a.coarse_mean_volumes.len() != family.limits.len()
        || a.fine_mean_volumes.len() != family.limits.len()
    {
        return Err(invalid(
            "assessment does not cover every original response limit",
        ));
    }
    let coarse = response_output::measurement(family, &response.coarse)?;
    let fine = response_output::measurement(family, &response.fine)?;
    let mut rows = String::new();
    for (index, limit) in family.limits.iter().enumerate() {
        let coarse_volume = a.coarse_mean_volumes[index];
        let fine_volume = a.fine_mean_volumes[index];
        let mean = matches!(
            limit.observation,
            spec::responses::Observation::MeanDisplacement { .. }
        );
        if [coarse_volume, fine_volume].iter().any(|v| match v {
            Some(v) => !mean || !v.is_finite() || *v <= 0.0,
            None => mean,
        }) {
            return Err(invalid(
                "assessment mean-response normalization is missing or invalid",
            ));
        }
        let volume_text = |v: Option<f64>| v.map_or_else(|| "—".into(), |v| format!("{v:.8e}"));
        let _ = write!(
            rows,
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td><td>{}</td><td>{}</td></tr>",
            response_output::html_text(&limit.name),
            limit.observation.unit(),
            response_output::sense(limit.sense),
            limit.bound,
            response.coarse.witnesses[index].worst_value,
            response.fine.witnesses[index].worst_value,
            response.normalized_model_changes[index],
            volume_text(coarse_volume),
            volume_text(fine_volume),
        );
    }
    let json = format!(
        "{{\"coarse\":{coarse},\"fine\":{fine},\"normalized_model_changes\":{:?},\"maximum_model_change\":{:.17e},\"tolerance\":{:.17e},\"stable_and_feasible\":{},\"coarse_mean_volumes_m3\":{},\"fine_mean_volumes_m3\":{}}}",
        response.normalized_model_changes,
        response.maximum_model_change,
        p.response_tolerance,
        response.stable_and_feasible(p.response_tolerance),
        volumes(&a.coarse_mean_volumes),
        volumes(&a.fine_mean_volumes),
    );
    let html = format!(
        "<h3>Every original signed response on both grids</h3><p>Maximum normalized whole-response-model change: {:.8e}; tolerance: {:.8e}; strict bounds and model-change gate pass: {}. Each grid retains its own adverse coefficients and independently measured witness. Patch volumes below normalize the same fixed physical mean using each grid's native reference quadrature.</p><table><tr><th>Limit</th><th>Unit</th><th>Sense</th><th>Bound</th><th>Coarse adverse value</th><th>Fine adverse value</th><th>Normalized model change</th><th>Coarse patch volume (m³)</th><th>Fine patch volume (m³)</th></tr>{rows}</table>",
        response.maximum_model_change,
        p.response_tolerance,
        response.stable_and_feasible(p.response_tolerance),
    );
    Ok((json, html))
}

fn stress_change(c: &AffineStressChange3) -> String {
    format!(
        "{{\"cell\":{},\"position_m\":{:?},\"quadrature_weight_m3\":{:.17e},\"normalized_reference_change\":{:.17e},\"normalized_relaxed_change\":{:.17e},\"normalized_physical_change\":{:.17e}}}",
        c.cell,
        c.position,
        c.weight,
        c.normalized_reference_change,
        c.normalized_relaxed_change,
        c.normalized_physical_change,
    )
}

fn strength(config: Config, r: &AffineLoadSpaceRefinement3, p: Policy) -> Result<(String, String)> {
    let s = match (config.stress, &r.stress) {
        (None, None) => return Ok(("null".into(), String::new())),
        (Some(limit), Some(s)) if s.coarse.limit == limit && s.fine.limit == limit => s,
        _ => {
            return Err(invalid(
                "assessment strength policy differs from its source",
            ));
        }
    };
    let first = s
        .changes
        .first()
        .ok_or_else(|| invalid("assessment strength has no matched fine-site changes"))?;
    let governing = |measure: fn(&AffineStressChange3) -> f64| {
        s.changes.iter().skip(1).fold(first, |best, current| {
            if measure(current) > measure(best) {
                current
            } else {
                best
            }
        })
    };
    let physical = governing(|v| v.normalized_physical_change);
    let relaxed = governing(|v| v.normalized_relaxed_change);
    let json = format!(
        "{{\"coarse\":{},\"fine\":{},\"normalized_aggregate_change\":{:.17e},\"maximum_model_change\":{:.17e},\"maximum_reference_change\":{:.17e},\"maximum_relaxed_change\":{:.17e},\"maximum_physical_change\":{:.17e},\"matched_points\":{},\"governing_physical_change\":{},\"governing_relaxed_change\":{},\"tolerance\":{:.17e},\"stable_and_feasible\":{}}}",
        strength::measurement(&s.coarse)?,
        strength::measurement(&s.fine)?,
        s.normalized_aggregate_change,
        s.maximum_model_change,
        s.maximum_reference_change,
        s.maximum_relaxed_change,
        s.maximum_physical_change,
        s.changes.len(),
        stress_change(physical),
        stress_change(relaxed),
        p.stress_tolerance,
        s.stable_and_feasible(p.stress_tolerance),
    );
    let html = format!(
        "<h3>Strength on both grids and matched fine sites</h3><p>The original stress cap is {:.8e} Pa. Coarse sampled aggregate: {:.8e} Pa at {} points (feasible: {}); fine aggregate: {:.8e} Pa at {} points (feasible: {}). Normalized aggregate change: {:.8e}. The maximum whole-family physical and relaxed field changes at {} matched fine sites are {:.8e} and {:.8e}; tolerance: {:.8e}. Strict caps, aggregate change and matched field-change gate pass: {}.</p>",
        s.coarse.limit.limit,
        s.coarse.aggregate,
        s.coarse.point_count,
        s.coarse.physically_feasible(),
        s.fine.aggregate,
        s.fine.point_count,
        s.fine.physically_feasible(),
        s.normalized_aggregate_change,
        s.changes.len(),
        s.maximum_physical_change,
        s.maximum_relaxed_change,
        p.stress_tolerance,
        s.stable_and_feasible(p.stress_tolerance),
    );
    Ok((json, html))
}

fn result(
    spec: &Spec,
    state: &State,
    a: &Assessment,
    r: &AffineLoadSpaceRefinement3,
    p: Policy,
) -> Result<(String, String)> {
    let config = spec
        .affine_load_space
        .ok_or_else(|| invalid("missing affine design policy"))?;
    if r.rho != state.selected().rho
        || r.policy != config.policy(spec.loads.len())?
        || r.robustness != config.options.robustness
        || r.goal != Some(config.goal)
        || r.constraints != config.constraint_policy()
    {
        return Err(invalid(
            "assessment differs from the selected design or original physical problem",
        ));
    }
    let joint = r.joint_constraints().map_err(|e| invalid(e.to_string()))?;
    let inventory = r
        .inventory
        .as_ref()
        .ok_or_else(|| invalid("assessment lacks dilated inventory"))?;
    let objective_change = r
        .objective_change()
        .ok_or_else(|| invalid("invalid assessed objective change"))?;
    let inventory_stable = inventory.volume_change.is_finite()
        && inventory.volume_change >= 0.0
        && inventory.volume_change <= p.inventory_tolerance;
    let objective_tolerance = p.objective_tolerance(config.goal);
    let objective_stable = objective_change <= objective_tolerance;
    let (
        objective_name, objective_unit, coarse_objective, fine_objective,
        normalized_coarse, normalized_fine, inventory_role, requirement,
    ) = match config.goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit } => (
            "compliance", "J", r.coarse.worst_compliance, r.fine.worst_compliance,
            r.coarse.maximum, r.fine.maximum, "volume-constraint",
            format!("The original dilated-volume cap is {volume_limit:.8e}. The compliance reference {:.8e} J supplies objective normalization.", config.cap_j),
        ),
        AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => (
            "dilated-volume-fraction", "1", inventory.coarse_volume_fraction,
            inventory.fine_volume_fraction, inventory.coarse_objective,
            inventory.fine_objective, "design-objective",
            format!("The original compliance cap is {:.8e} J.", config.cap_j),
        ),
    };
    let native_pass = r.stable_and_feasible_with_constraint_tolerances(
        p.compliance_tolerance,
        p.response_tolerance,
        p.stress_tolerance,
    );
    let (responses, response_html) = responses(spec, a, r, p)?;
    let (strength, strength_html) = strength(config, r, p)?;
    let json = format!(
        "{{\"coarse\":{},\"fine\":{},\"normalized_model_change\":{:.17e},\"normalized_maximum_change\":{:.17e},\"normalized_witness_change\":{:.17e},\"native_stable_and_feasible\":{native_pass},\"design_objective\":{{\"kind\":{objective_name:?},\"unit\":{objective_unit:?},\"coarse_value\":{coarse_objective:.17e},\"fine_value\":{fine_objective:.17e},\"normalized_coarse\":{normalized_coarse:.17e},\"normalized_fine\":{normalized_fine:.17e},\"normalized_change\":{objective_change:.17e},\"tolerance\":{objective_tolerance:.17e},\"stable\":{objective_stable}}},\"inventory\":{{\"objective\":\"dilated-volume-fraction\",\"role\":{inventory_role:?},\"volume_scale\":{:.17e},\"coarse_volume_fraction\":{:.17e},\"fine_volume_fraction\":{:.17e},\"volume_change\":{:.17e},\"coarse_objective\":{:.17e},\"fine_objective\":{:.17e},\"objective_change\":{:.17e},\"tolerance\":{:.17e},\"tolerance_applies_to\":\"absolute-volume-fraction-change\",\"stable\":{inventory_stable}}},\"responses\":{responses},\"strength\":{strength},\"dwr\":{{\"coarse_maximizer\":{},\"fine_maximizer\":{},\"fine_switch_correction\":{:.17e},\"coarse_switch_correction\":{:.17e},\"normalized_reconstructed_change\":{:.17e},\"identity_relative_defect\":{:.17e},\"diagnostic_marking_cells\":{}}}}}",
        grid(&r.coarse, config, inventory.coarse_volume_fraction, joint[0]),
        grid(&r.fine, config, inventory.fine_volume_fraction, joint[1]),
        r.normalized_model_change,
        r.normalized_maximum_change,
        r.normalized_witness_change,
        inventory.volume_scale,
        inventory.coarse_volume_fraction,
        inventory.fine_volume_fraction,
        inventory.volume_change,
        inventory.coarse_objective,
        inventory.fine_objective,
        inventory.objective_change,
        p.inventory_tolerance,
        goal(r, false),
        goal(r, true),
        r.fine_switch_correction,
        r.coarse_switch_correction,
        r.correction(),
        r.identity_relative_defect,
        r.marking_mass.len(),
    );
    let html = format!(
        "<p>Worst compliance on the coarse and fine grids: {:.8e} J and {:.8e} J. {requirement} Their complete joint constraints are {:.8e} and {:.8e}. The normalized whole-affine-model change is {:.8e}, against tolerance {:.8e}; the signed change of independently maximized values is {:.8e}. Native physical constraints and full-family change gates pass: {native_pass}.</p><p>The optimized {objective_name} changes from {coarse_objective:.8e} {objective_unit} to {fine_objective:.8e} {objective_unit}. Its normalized objective change is {objective_change:.8e}, against tolerance {objective_tolerance:.8e}; objective gate passes: {objective_stable}.</p><p>Dilated volume fraction changes from {:.8e} to {:.8e}. The separate absolute volume-fraction change is {:.8e}, against tolerance {:.8e}; inventory gate passes: {inventory_stable}. The two-grid identity defect is {:.8e}; the fine and coarse load-switch corrections are {:.8e} and {:.8e}. The retained JSON keeps both frozen-direction DWR decompositions and numerical residuals.</p>{response_html}{strength_html}",
        r.coarse.worst_compliance,
        r.fine.worst_compliance,
        joint[0],
        joint[1],
        r.normalized_model_change,
        p.compliance_tolerance,
        r.normalized_maximum_change,
        inventory.coarse_volume_fraction,
        inventory.fine_volume_fraction,
        inventory.volume_change,
        p.inventory_tolerance,
        r.identity_relative_defect,
        r.fine_switch_correction,
        r.coarse_switch_correction,
    );
    Ok((json, html))
}

pub(super) fn render(spec: &Spec, state: &State) -> Result<Option<(String, String)>> {
    let Some(p) = spec.affine_grid_assessment else {
        return if state.assessment.is_none() {
            Ok(None)
        } else {
            Err(invalid("unrequested grid assessment cannot be published"))
        };
    };
    let config = spec
        .affine_load_space
        .ok_or_else(|| invalid("missing affine design policy"))?;
    let selection = if state.best.is_some() {
        match config.goal {
            AffineLoadSpaceDesignGoal3::WorstCompliance { .. } => "least-compliance-feasible",
            AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => "least-material-feasible",
        }
    } else {
        "last-accepted"
    };
    let declared = policy(p);
    let Some(a) = &state.assessment else {
        return Ok(Some((
            format!(
                "{{\"schema\":\"sdf3-affine-grid-assessment.v1\",\"status\":\"pending\",\"passed\":false,\"assessed_design\":{selection:?},\"policy\":{declared},\"coarse_grid\":null,\"fine_grid\":null,\"reason\":null,\"work\":null,\"result\":null,\"no_claim\":{}}}",
                quoted(SCOPE),
            ),
            format!(
                "<h2>Enriched-grid assessment</h2><p>The requested assessment of the selected design is pending. No two-grid result has been produced for this retained point.</p><p>{SCOPE}</p>"
            ),
        )));
    };
    if !matches!(
        a.status,
        "checked" | "budget-exhausted" | "cancelled" | "numerical-failure"
    ) || (a.status == "checked") != a.report.is_some()
    {
        return Err(invalid("inconsistent grid-assessment completion state"));
    }
    let (result, details) = a
        .report
        .as_ref()
        .map(|r| result(spec, state, a, r, p))
        .transpose()?
        .unwrap_or_else(|| ("null".into(), String::new()));
    let work = Spent {
        wall_s: 0.0,
        linear: a.work,
        geometry: a.geometry,
    }
    .json();
    let reason = a
        .reason
        .as_ref()
        .map_or_else(|| "null".into(), |r| quoted(r));
    let reason_html = a.reason.as_ref().map_or_else(String::new, |r| {
        format!("<p>Assessment stop: {}.</p>", response_output::html_text(r))
    });
    let passed = a.passed(p, config.goal);
    let json = format!(
        "{{\"schema\":\"sdf3-affine-grid-assessment.v1\",\"status\":{:?},\"passed\":{passed},\"assessed_design\":{selection:?},\"policy\":{declared},\"coarse_grid\":{{\"cells\":{},\"dofs\":{}}},\"fine_grid\":{{\"cells\":{},\"dofs\":{}}},\"reason\":{reason},\"work\":{work},\"result\":{result},\"no_claim\":{}}}",
        a.status,
        a.coarse_cells,
        a.coarse_dofs,
        a.fine_cells,
        a.fine_dofs,
        quoted(SCOPE),
    );
    let html = format!(
        "<h2>Enriched-grid assessment</h2><p>Assessment status: {}. All strict physical and declared grid-change gates pass: {passed}. The assessed optimizer selection is {selection}. Coarse grid: {} cells and {} DOFs; fine grid: {} cells and {} DOFs.</p>{reason_html}{details}<p>{SCOPE}</p>",
        a.status, a.coarse_cells, a.coarse_dofs, a.fine_cells, a.fine_dofs,
    );
    Ok(Some((json, html)))
}
