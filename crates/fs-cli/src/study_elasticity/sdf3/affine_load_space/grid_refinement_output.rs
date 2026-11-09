//! Keep the installed design distinct from a resumable finer-grid candidate.
use super::super::grid_refinement::Policy;
use super::*;
use fs_solver::LinearOp;

const SCOPE: &str = "One bounded native refit inherits raw design controls and protected regions, rebuilds self-weight and reintegrates the original fixed mechanical preload on the finer grid, and starts a separate optimizer with its own declared update and evaluation allowances. Geometry, Krylov and wall work remain cumulative. Until a completed candidate passes a fresh complete-load-family physical recheck, the selected design and its field export remain the immutable source design. Promotion does not establish cross-grid objective descent or continuum feasibility; the promoted design receives its own enriched-grid assessment. Iteration histories, gradient audits and stopping diagnostics belong to their respective working grids.";

fn invalid(message: impl Into<String>) -> Failure {
    fail("cli-study-sdf3-grid-refinement-output", message)
}

fn document(bytes: &[u8]) -> Result<JsonValue> {
    JsonValue::parse(std::str::from_utf8(bytes).map_err(|e| invalid(e.to_string()))?)
        .map_err(|e| invalid(e.to_string()))
}

fn field<'a>(value: &'a JsonValue, key: &str) -> Result<&'a JsonValue> {
    value
        .get(key)
        .ok_or_else(|| invalid(format!("missing retained source {key}")))
}

fn real(value: &JsonValue, key: &str) -> Result<f64> {
    field(value, key)?
        .as_f64()
        .filter(|x| x.is_finite())
        .ok_or_else(|| invalid(format!("nonfinite retained source {key}")))
}

fn boolean(value: &JsonValue, key: &str) -> Result<bool> {
    match field(value, key)? {
        JsonValue::Bool(value) => Ok(*value),
        _ => Err(invalid(format!("invalid retained source {key}"))),
    }
}

// Retain the exact spelling of stored numerical values. The bounded JSON
// reader already admits nesting and refuses duplicate object members.
fn json(value: &JsonValue) -> Result<String> {
    Ok(match value {
        JsonValue::Null => "null".into(),
        JsonValue::Bool(value) => value.to_string(),
        JsonValue::Number { value, raw } if value.is_finite() => raw.clone(),
        JsonValue::Number { .. } => return Err(invalid("nonfinite retained source value")),
        JsonValue::Str(value) => quoted(value),
        JsonValue::Array(values) => format!(
            "[{}]",
            values
                .iter()
                .map(json)
                .collect::<Result<Vec<_>>>()?
                .join(",")
        ),
        JsonValue::Object(values) => format!(
            "{{{}}}",
            values
                .iter()
                .map(|(key, value)| Ok(format!("{}:{}", quoted(key), json(value)?)))
                .collect::<Result<Vec<_>>>()?
                .join(",")
        ),
    })
}

pub(super) fn policy(p: Policy) -> String {
    format!(
        "{{\"max_refits\":1,\"candidate_max_updates\":{},\"candidate_max_optimizer_evaluations\":{}}}",
        p.max_updates, p.max_evaluations
    )
}

pub(super) struct Source {
    pub design: String,
    pub selection: String,
    pub feasible: bool,
    pub measurement: String,
    pub volume: f64,
    pub compliance: f64,
    pub assessment: String,
    receipt: JsonValue,
    report: JsonValue,
    cells: usize,
    nodes: usize,
}

pub(super) fn source(
    spec: &Spec,
    ledger: &Ledger,
    state: &State,
    producer: ContentHash,
) -> Result<Option<Source>> {
    let Some(context) = &state.refinement else {
        return Ok(None);
    };
    let config = spec.affine_load_space
        .ok_or_else(|| invalid("missing source affine design goal"))?;
    let expected_goal = document(goal_model(config).as_bytes())?;
    let (objective, objective_unit) = objective_metadata(config.goal);
    let old = load(
        ledger,
        &format!("study-{}", context.source_receipt.to_hex()),
    )?;
    if old.value.str_field("driver") != Some(AFFINE_LOAD_SPACE3_DRIVER)
        || old.value.str_field("study_id") != Some(spec.id.to_hex().as_str())
        || old.value.str_field("producer") != Some(producer.to_hex().as_str())
        || old.value.str_field("source")
            != Some(hash_bytes(spec.canonical.as_bytes()).to_hex().as_str())
        || old.value.get("grid_refinement").is_some()
    {
        return Err(invalid(
            "refinement source must be this study's original sealed grid endpoint",
        ));
    }
    let report = document(&linked(
        ledger,
        &old.value,
        "report_json",
        "study-report-json",
    )?)?;
    let design = String::from_utf8(linked(ledger, &old.value, "design", "study-design")?)
        .map_err(|e| invalid(e.to_string()))?;
    let fields = document(design.as_bytes())?;
    let selection = report
        .str_field("selected_design")
        .filter(|s| *s == feasible_selection(config.goal) || *s == "last-accepted")
        .ok_or_else(|| invalid("missing retained source selection"))?
        .to_owned();
    let feasible = boolean(&report, "selected_feasible")?;
    let selected = field(&report, "selected")?;
    let assessment = field(&report, "grid_assessment")?;
    if report.str_field("driver") != Some(AFFINE_LOAD_SPACE3_DRIVER)
        || report.str_field("study_id") != old.value.str_field("study_id")
        || report.str_field("producer") != old.value.str_field("producer")
        || report.str_field("objective") != Some(objective)
        || report.str_field("objective_unit") != Some(objective_unit)
        || report.get("goal") != Some(&expected_goal)
        || fields.str_field("schema") != Some("sdf3-affine-load-space-design.v1")
        || fields.get("goal") != Some(&expected_goal)
        || fields.str_field("selected_design") != Some(selection.as_str())
        || old.value.str_field("selected_design") != Some(selection.as_str())
        || boolean(&fields, "selected_feasible")? != feasible
        || boolean(&old.value, "selected_feasible")? != feasible
        || boolean(selected, "feasible")? != feasible
        || !feasible
        || assessment.str_field("status") != Some("checked")
        || boolean(assessment, "passed")?
    {
        return Err(invalid(
            "refinement source fields, selected result or unresolved assessment disagree",
        ));
    }
    let selected_fields = field(&fields, "selected")?;
    let count = |key: &str| -> Result<usize> {
        field(selected_fields, key)?
            .as_array()
            .map(|values| values.len())
            .filter(|n| *n > 0)
            .ok_or_else(|| invalid(format!("empty retained source {key}")))
    };
    Ok(Some(Source {
        cells: count("cells")?,
        nodes: count("physical_nodes_m")?,
        measurement: json(selected)?,
        volume: real(selected, "dilated_volume_fraction")?,
        compliance: real(selected, "worst_compliance_j")?,
        assessment: json(assessment)?,
        design,
        selection,
        feasible,
        receipt: old.value,
        report,
    }))
}

impl Source {
    pub(super) fn selected_constraints_html(&self) -> Result<String> {
        let selected = field(&self.report, "selected")?;
        let mut html = String::new();
        if let Some(strength) = selected
            .get("strength")
            .filter(|v| !matches!(v, JsonValue::Null))
        {
            let _ = write!(
                html,
                "<h2>Retained source strength</h2><p>The selected source design has sampled affine strength aggregate {:.8e} Pa against limit {:.8e} Pa; cap-feasible: {}. The complete joint constraint is {:.8e}.</p>",
                real(strength, "aggregate_pa")?,
                real(strength, "stress_limit_pa")?,
                boolean(strength, "feasible")?,
                real(selected, "joint_constraint")?
            );
        }
        if let Some(responses) = selected.get("responses") {
            let limits = field(responses, "limits")?
                .as_array()
                .ok_or_else(|| invalid("missing source signed response rows"))?;
            let _ = write!(
                html,
                "<h2>Retained source signed responses</h2><p>Physical bounds and response aggregate pass: {}. Each row retains its own adverse load over the original complete affine family.</p><table><tr><th>Limit</th><th>Unit</th><th>Sense</th><th>Bound</th><th>Adverse value</th><th>Measured witness</th><th>Normalized residual</th></tr>",
                boolean(responses, "feasible")?
            );
            for row in limits {
                let text = |key| {
                    row.str_field(key)
                        .map(response_output::html_text)
                        .ok_or_else(|| invalid(format!("missing source response {key}")))
                };
                let _ = write!(
                    html,
                    "<tr><td>{}</td><td>{}</td><td>{}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td></tr>",
                    text("name")?,
                    text("unit")?,
                    text("sense")?,
                    real(row, "bound")?,
                    real(row, "worst_value")?,
                    real(row, "witness_value")?,
                    real(row, "residual")?
                );
            }
            html.push_str("</table>");
        }
        Ok(html)
    }

    pub(super) fn assessment_html(&self, config: Config) -> Result<String> {
        let assessment = field(&self.report, "grid_assessment")?;
        let result = field(assessment, "result")?;
        let coarse = field(result, "coarse")?;
        let fine = field(result, "fine")?;
        let requirement = match config.goal {
            AffineLoadSpaceDesignGoal3::WorstCompliance { .. } => format!(
                "The dilated volume cap is {:.8e}; coarse and fine volume fractions are {:.8e} and {:.8e}. The {:.8e} J compliance reference is a numerical normalization scale, not a physical upper bound.",
                real(coarse, "volume_limit_fraction")?,
                real(coarse, "dilated_volume_fraction")?,
                real(fine, "dilated_volume_fraction")?,
                real(coarse, "compliance_reference_j")?,
            ),
            AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => format!(
                "The original physical compliance cap is {:.8e} J.",
                real(coarse, "compliance_limit_j")?,
            ),
        };
        Ok(format!(
            "<h2>Selected source grid assessment</h2><p>The original complete-family assessment is checked and unresolved. Its coarse and fine worst compliance values are {:.8e} J and {:.8e} J; joint constraints are {:.8e} and {:.8e}. {requirement} The whole-family normalized compliance-model change is {:.8e}. This unchanged assessment belongs to the retained source design; the working candidate has not replaced it. All original signed-bound rows, sampled-strength evidence, objective and inventory changes remain in the retained report JSON.</p>",
            real(coarse, "worst_compliance_j")?,
            real(fine, "worst_compliance_j")?,
            real(coarse, "joint_constraint")?,
            real(fine, "joint_constraint")?,
            real(result, "normalized_model_change")?,
        ))
    }
}

pub(super) fn render(
    spec: &Spec,
    study: &CutDensityStudy3<AdaptiveSolveSpace3>,
    state: &State,
    source: &Source,
    working_selected: &str,
) -> Result<(String, String)> {
    let context = state
        .refinement
        .as_ref()
        .ok_or_else(|| invalid("missing refinement context"))?;
    let p = spec
        .affine_grid_refinement
        .ok_or_else(|| invalid("unrequested grid refinement"))?;
    let config = spec
        .affine_load_space
        .ok_or_else(|| invalid("missing affine refinement optimizer"))?;
    if context.attempt != 1 {
        return Err(invalid("refinement exceeds its single candidate identity"));
    }
    let phase = context.phase.as_str();
    let generation = context.working_generation();
    let retained = context.uses_source_selection();
    let selected_generation = if retained { 0 } else { generation };
    let recheck = match &state.refinement_recheck {
        None => "{\"status\":\"pending\",\"feasible\":false,\"measurement\":null,\"reason\":null}"
            .into(),
        Some(Ok(point)) => {
            if point.rho != state.selected().rho {
                return Err(invalid(
                    "fresh recheck differs from the working selected design",
                ));
            }
            format!(
                "{{\"status\":\"checked\",\"feasible\":{},\"measurement\":{},\"reason\":null}}",
                config.feasible(point),
                measurement(point, config, spec)?
            )
        }
        Some(Err(error)) => format!(
            "{{\"status\":\"failed\",\"feasible\":false,\"measurement\":null,\"reason\":{}}}",
            quoted(&error.to_string())
        ),
    };
    if !retained
        && !state
            .refinement_recheck
            .as_ref()
            .is_some_and(|r| r.as_ref().is_ok_and(|point| config.feasible(point)))
    {
        return Err(invalid(
            "promoted design lacks its successful complete-family physical recheck",
        ));
    }
    let source_selected = field(&source.report, "selected")?;
    let goal = goal_model(config);
    let (objective, objective_unit) = objective_metadata(config.goal);
    let source_objective = match config.goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { .. } => source.compliance,
        AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => source.volume,
    };
    let source_grid = format!(
        "{{\"objective\":{objective:?},\"objective_unit\":{objective_unit:?},\"objective_value\":{source_objective:.17e},\"cells\":{},\"physical_nodes\":{},\"iterations_completed\":{},\"target_iterations\":{},\"optimization_status\":{},\"optimizer_work\":{},\"iterations_artifact\":{},\"selected_feasible\":{},\"dilated_volume_fraction\":{:.17e},\"worst_compliance_j\":{:.17e},\"joint_constraint\":{:.17e},\"assessment_status\":\"checked\",\"assessment_passed\":false}}",
        source.cells,
        source.nodes,
        integer(&source.receipt, "iterations_completed")?,
        integer(&source.receipt, "target_iterations")?,
        json(field(&source.report, "optimization_status")?)?,
        json(field(&source.report, "optimizer_work")?)?,
        json(field(&source.receipt, "iterations")?)?,
        source.feasible,
        source.volume,
        source.compliance,
        real(source_selected, "joint_constraint")?,
    );
    let operator = study.operator().elasticity();
    let working_grid = format!(
        "{{\"generation\":{generation},\"cells\":{},\"physical_nodes\":{},\"dofs\":{},\"iterations_completed\":{},\"target_iterations\":{},\"max_optimizer_evaluations\":{}}}",
        study.cells(),
        operator.physical_nodes().len(),
        operator.n(),
        state.optimizer.iterations,
        spec.updates,
        config.options.optimizer.max_evaluations,
    );
    let base = context.json();
    let base = base
        .strip_suffix('}')
        .ok_or_else(|| invalid("invalid refinement context object"))?;
    let failure = state
        .refinement_failure
        .as_ref()
        .map_or_else(|| "null".into(), |reason| quoted(reason));
    let json = format!(
        "{base},\"goal\":{goal},\"policy\":{},\"selected_generation\":{selected_generation},\"last_accepted_generation\":{generation},\"gradient_check_generation\":{generation},\"source_grid\":{source_grid},\"working_grid\":{working_grid},\"working_selected\":{working_selected},\"fresh_recheck\":{recheck},\"failure\":{failure},\"history_scope\":\"working-grid-only\",\"physical_budget_scope\":\"cumulative-across-source-candidate-assessments-and-restoration\",\"no_claim\":{}}}",
        policy(p),
        quoted(SCOPE),
    );
    let selection = if retained {
        "The published design and VTU remain exactly the source selection. The working design is retained separately and is not installed."
    } else {
        "The selected design and VTU now use the promoted grid and its actual physical fields."
    };
    let recheck_html = match &state.refinement_recheck {
        None => "Fresh complete-family candidate recheck: pending.".into(),
        Some(Ok(point)) => format!(
            "Fresh complete-family candidate recheck: checked; every original physical cap and joint constraint passes: {}.",
            config.feasible(point)
        ),
        Some(Err(error)) => format!(
            "Fresh complete-family candidate recheck stopped: {}.",
            response_output::html_text(&error.to_string())
        ),
    };
    let failure_html = state
        .refinement_failure
        .as_ref()
        .map_or_else(String::new, |reason| {
            format!(
                "<p>Refinement stopped: {}.</p>",
                response_output::html_text(reason)
            )
        });
    let html = format!(
        "<h2>Automatic grid refinement</h2><p>Phase: {phase}. Attempt: 1 of 1. {selection}</p><p>Source grid: {} active cells, {} physical nodes and {} retained updates of {}. Working grid generation {generation}: {} active cells, {} physical nodes and {} updates of {}; its optimizer evaluation allowance is {}. Gradient, optimizer stopping/KKT values and the iteration table below apply only to the working grid. Original source history remains linked separately.</p><p>{recheck_html}</p>{failure_html}<p>{SCOPE}</p>",
        source.cells,
        source.nodes,
        integer(&source.receipt, "iterations_completed")?,
        integer(&source.receipt, "target_iterations")?,
        study.cells(),
        operator.physical_nodes().len(),
        state.optimizer.iterations,
        spec.updates,
        config.options.optimizer.max_evaluations,
    );
    Ok((json, html))
}
