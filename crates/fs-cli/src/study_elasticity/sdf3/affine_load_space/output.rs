//! Retain the actual affine witness and selected physical design in the ledger.
use super::*;

#[path = "grid_assessment_output.rs"]
mod grid_assessment_output;
#[path = "grid_refinement_output.rs"]
mod grid_refinement_output;
#[path = "response_output.rs"]
mod response_output;
#[path = "strength_output.rs"]
mod strength;

fn objective_metadata(goal: AffineLoadSpaceDesignGoal3) -> (&'static str, &'static str) {
    match goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { .. } => ("compliance", "J"),
        AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => ("dilated-volume-fraction", "1"),
    }
}

fn feasible_selection(goal: AffineLoadSpaceDesignGoal3) -> &'static str {
    match goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { .. } => "least-compliance-feasible",
        AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => "least-material-feasible",
    }
}

/// Keep the numerical compliance reference separate from the authored cap.
/// In stiffness design the physical base constraint is dilated inventory.
fn goal_model(config: Config) -> String {
    let (objective, unit) = objective_metadata(config.goal);
    let (kind, compliance_limit, volume_limit) = match config.goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit } => (
            "worst-compliance", "null".into(), format!("{volume_limit:.17e}"),
        ),
        AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => (
            "minimum-material", format!("{:.17e}", config.cap_j), "null".into(),
        ),
    };
    format!(
        "{{\"kind\":{kind:?},\"objective\":{objective:?},\"objective_unit\":{unit:?},\"compliance_reference_j\":{:.17e},\"compliance_limit_j\":{compliance_limit},\"volume_limit_fraction\":{volume_limit}}}",
        config.cap_j,
    )
}

fn has_dead_load(spec: &Spec) -> bool {
    spec.affine_dead_load.is_some_and(|law|
        law.body_n_m3.iter().any(|value| *value != 0.0) || law.surface.is_some())
}

/// The authored reference preload is independent of occupancy and uncertainty.
/// Its displacement is combined with native material self-weight in the common
/// offset; the uncertain basis columns remain separate native solves.
fn dead_load_model(spec: &Spec) -> String {
    let Some(law) = spec.affine_dead_load else { return "null".into(); };
    let surface = law.surface.map_or_else(|| "null".into(), |surface| {
        let (kind, pa) = match surface.force {
            loading::Force::Pressure(value) => ("pressure", format!("{value:.17e}")),
            loading::Force::Traction(value) => ("traction", format!("{value:?}")),
        };
        format!("{{\"kind\":{kind:?},\"pa\":{pa},\"x_fraction\":{:?}}}", surface.x_fraction)
    });
    format!("{{\"body_n_m3\":{:?},\"surface\":{surface}}}", law.body_n_m3)
}

pub(super) fn measurement(
    point: &AffineLoadSpaceDesignPoint3,
    config: Config,
    spec: &Spec,
) -> Result<String> {
    let e = &point.evaluation;
    let p = &e.physical;
    let o = &e.oracle;
    if e.stress.as_ref().map(|s| s.limit) != config.stress {
        return Err(fail(
            "cli-study-sdf3-strength-output",
            "retained strength policy differs from its source",
        ));
    }
    let constraint = point
        .constraint(config.goal)
        .map_err(|e| fail("cli-study-sdf3-strength-output", e.to_string()))?;
    let strength = e
        .stress
        .as_ref()
        .map(strength::measurement)
        .transpose()?
        .unwrap_or_else(|| "null".into());
    let responses = match (&spec.affine_responses, &e.responses) {
        (None, None) => String::new(),
        (Some(family), Some(evaluation)) => format!(
            ",\"responses\":{}",
            response_output::measurement(family, evaluation)?
        ),
        _ => {
            return Err(fail(
                "cli-study-sdf3-response-output",
                "retained response family differs from its source",
            ));
        }
    };
    let (_, objective_unit) = objective_metadata(config.goal);
    let (objective_value, normalized_objective, cap_fields) = match config.goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit } => (
            e.worst_compliance, e.maximum,
            format!("\"compliance_reference_j\":{:.17e},\"compliance_limit_j\":null,\"compliance_feasible\":null,\"volume_limit_fraction\":{volume_limit:.17e},\"volume_feasible\":{}",
                config.cap_j, p.dilated_volume_fraction.is_finite() && p.dilated_volume_fraction <= volume_limit),
        ),
        AffineLoadSpaceDesignGoal3::MinimumMaterial { volume_scale, .. } => (
            p.dilated_volume_fraction, p.dilated_volume_fraction / volume_scale,
            format!("\"compliance_reference_j\":{:.17e},\"compliance_limit_j\":{:.17e},\"compliance_feasible\":{},\"volume_limit_fraction\":null,\"volume_feasible\":null",
                config.cap_j, config.cap_j, e.worst_compliance.is_finite()
                    && e.worst_compliance >= 0.0 && e.worst_compliance <= config.cap_j
                    && (0.0..=1.0).contains(&e.maximum)),
        ),
    };
    Ok(format!(
        "{{\"objective_value\":{objective_value:.17e},\"objective_unit\":{objective_unit:?},\"normalized_objective\":{normalized_objective:.17e},{cap_fields},\"dilated_volume_fraction\":{:.17e},\"nominal_volume_fraction\":{:.17e},\"eroded_volume_fraction\":{:.17e},\"worst_compliance_j\":{:.17e},\"maximum_ratio\":{:.17e},\"feasible\":{},\"joint_constraint\":{constraint:.17e},\"strength\":{strength}{responses},\"worst_coefficients\":{:?},\"normalized_gram_row_major\":{:?},\"normalized_linear\":{:?},\"normalized_offset\":{:.17e},\"reciprocity_defect\":{:.17e},\"witness_relative_defect\":{:.17e},\"witness_relative_residual\":{:.17e},\"oracle\":{{\"branch\":{:?},\"gradient_admissible\":{},\"value_without_offset\":{:.17e},\"dual_value_without_offset\":{:.17e},\"duality_gap\":{:.17e},\"relative_duality_gap\":{:.17e},\"relative_separation\":{:.17e},\"stationarity_residual\":{:.17e},\"feasibility_residual\":{:.17e},\"complementarity_residual\":{:.17e},\"eigenpair_residual\":{:.17e},\"iterations\":{}}}}}",
        p.dilated_volume_fraction,
        p.nominal_volume_fraction,
        p.eroded_volume_fraction,
        e.worst_compliance,
        e.maximum,
        config.feasible(point),
        e.worst_coefficients,
        e.gram,
        e.linear,
        e.offset,
        e.reciprocity_defect,
        e.witness_relative_defect,
        e.witness_relative_residual,
        format!("{:?}", o.branch),
        o.gradient_admissible(),
        o.value,
        o.dual_value,
        o.duality_gap,
        o.relative_duality_gap,
        o.relative_separation,
        o.stationarity_residual,
        o.feasibility_residual,
        o.complementarity_residual,
        o.eigenpair_residual,
        o.iterations,
    ))
}

fn fields(
    study: &CutDensityStudy3<AdaptiveSolveSpace3>,
    point: &AffineLoadSpaceDesignPoint3,
    spec: &Spec,
) -> Result<String> {
    let e = &point.evaluation;
    let operator = study.operator().elasticity();
    let volumes = operator.volumes();
    let mut cells = String::new();
    for (i, cell) in operator.leaves().iter().enumerate() {
        if i != 0 {
            cells.push(',');
        }
        let region = study
            .physical_regions()
            .map_or(fs_topopt::sdf3::PhysicalRegion3::Design, |r| r[i])
            .label();
        let _ = write!(
            cells,
            "{{\"level\":{},\"index\":{:?},\"physical_region\":{region:?},\"raw_density\":{:.17e},\"eroded_density\":{:.17e},\"nominal_density\":{:.17e},\"dilated_density\":{:.17e},\"cut_volume_m3\":{:.17e}}}",
            cell.level(),
            cell.index(),
            point.rho[i],
            e.physical.eroded_rho[i],
            e.physical.nominal_rho[i],
            e.physical.dilated_rho[i],
            volumes[i]
        );
    }
    let physical = |field: &[f64]| {
        operator
            .physical_displacements(field)
            .map_err(|e| fail("cli-study-sdf3-field", e.to_string()))
    };
    let basis = e
        .basis_displacements
        .iter()
        .map(|u| physical(u))
        .collect::<Result<Vec<_>>>()?;
    let offset = physical(&e.offset_displacement)?;
    let worst = physical(&e.worst_displacement)?;
    let strength = e
        .stress
        .as_ref()
        .map(|stress| strength::fields(stress, study.cells()))
        .transpose()?
        .map_or_else(String::new, |fields| format!(",\"strength\":{fields}"));
    let responses = match (&spec.affine_responses, &e.responses) {
        (None, None) => String::new(),
        (Some(family), Some(evaluation)) => format!(
            ",\"responses\":{}",
            response_output::fields(family, evaluation)?
        ),
        _ => {
            return Err(fail(
                "cli-study-sdf3-response-output",
                "retained response fields differ from their source",
            ));
        }
    };
    Ok(format!(
        "{{\"cells\":[{cells}],\"physical_nodes_m\":{:?},\"cell_nodes\":{:?},\"displacement_layout\":\"basis column then node-major xyz, on eroded material\",\"offset_displacement_scope\":\"complete common load f_dead+B*rho_eroded on eroded stiffness\",\"basis_displacements_m\":{basis:?},\"offset_displacement_m\":{offset:?},\"worst_displacement_m\":{worst:?}{strength}{responses}}}",
        operator.physical_nodes(),
        operator.cell_nodes(),
    ))
}

fn stopping(state: &State) -> String {
    match &state.report {
        Some(r) => format!(
            "{{\"reason\":{:?},\"kkt_scope\":\"last-accepted\",\"stationarity\":{:.17e},\"feasibility\":{:.17e},\"dual_feasibility\":{:.17e},\"complementarity\":{:.17e}}}",
            format!("{:?}", r.stop),
            r.kkt.stationarity,
            r.kkt.feasibility,
            r.kkt.dual_feasibility,
            r.kkt.complementarity
        ),
        None => format!(
            "{{\"reason\":\"EvaluationStopped\",\"error\":{},\"kkt\":null}}",
            state
                .error
                .as_ref()
                .map_or_else(|| "null".into(), |e| quoted(&e.to_string()))
        ),
    }
}

pub(super) fn persist(
    spec: &Spec,
    ledger: &Ledger,
    study: &CutDensityStudy3<AdaptiveSolveSpace3>,
    state: &State,
    producer: ContentHash,
    previous: Option<ContentHash>,
    limit: usize,
) -> Result<Outcome> {
    let config = spec.affine_load_space.ok_or_else(|| {
        fail(
            "cli-study-sdf3-input",
            "missing affine optimizer during publication",
        )
    })?;
    let law = spec.material_acceleration.ok_or_else(|| {
        fail(
            "cli-study-sdf3-input",
            "missing affine self-weight during publication",
        )
    })?;
    if state.gradient.goal != config.goal {
        return Err(fail("cli-study-sdf3-goal-output", "retained gradient experiment belongs to another design goal"));
    }
    let goal = goal_model(config);
    let (objective, objective_unit) = objective_metadata(config.goal);
    let source = grid_refinement_output::source(spec, ledger, state, producer)?;
    let retain_source = state
        .refinement
        .as_ref()
        .is_some_and(|context| context.uses_source_selection());
    let working_selection = if state.best.is_some() {
        feasible_selection(config.goal)
    } else {
        "last-accepted"
    };
    let working_feasible = config.feasible(state.selected());
    let working_selected = measurement(state.selected(), config, spec)?;
    let source_selection = source.as_ref().filter(|_| retain_source);
    let selection = source_selection.map_or(working_selection, |s| s.selection.as_str());
    let feasible = source_selection.map_or(working_feasible, |s| s.feasible);
    let selected = source_selection.map_or(working_selected.as_str(), |s| s.measurement.as_str());
    let accepted = measurement(&state.accepted, config, spec)?;
    let strength_marker = if config.stress.is_some() {
        ",\"strength_constrained\":true"
    } else {
        ""
    };
    let response_marker = if spec.affine_responses.is_some() {
        ",\"response_constrained\":true"
    } else {
        ""
    };
    let working_design = format!(
        "{{\"schema\":\"sdf3-affine-load-space-design.v1\"{strength_marker}{response_marker},\"goal\":{goal},\"selected_design\":{working_selection:?},\"selected_feasible\":{working_feasible},\"selected\":{},\"last_accepted\":{}}}",
        fields(study, state.selected(), spec)?,
        fields(study, &state.accepted, spec)?
    );
    // Never serialize a source point against the candidate operator. The
    // installed design remains the immutable source artifact until promotion.
    let design = source_selection.map_or(working_design.as_str(), |s| s.design.as_str());
    let mut history = String::new();
    let mut table = String::new();
    for (i, row) in state.history.iter().enumerate() {
        if i != 0 {
            history.push(',');
        }
        let _ = write!(
            history,
            "{{\"iteration\":{},\"objective\":{:.17e},\"constraint\":{:.17e},\"dilated_volume_fraction\":{:.17e},\"worst_compliance_j\":{:.17e},\"maximum_ratio\":{:.17e},\"feasible\":{}}}",
            row.iteration,
            row.objective,
            row.constraint,
            row.dilated_volume_fraction,
            row.worst_compliance,
            row.maximum,
            row.feasible
        );
        let _ = write!(
            table,
            "<tr><td>{}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td><td>{}</td></tr>",
            row.iteration, row.objective, row.dilated_volume_fraction, row.worst_compliance,
            row.maximum, row.constraint, row.feasible
        );
    }
    let rows =
        format!("{{\"schema\":\"sdf3-affine-load-space-iterations.v1\",\"history\":[{history}]}}");
    let work = state.spent.json();
    let stop = stopping(state);
    let gradient = &state.gradient;
    let base_constraint = match config.goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { .. } => "volume-only",
        AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => "compliance-only",
    };
    let constraint_policy = gradient.constraints.stress.map_or_else(
        || {
            format!(
                "{{\"kind\":{base_constraint:?},\"joint_sharpness\":{:.17e}}}",
                gradient.constraints.sharpness
            )
        },
        |limit| strength::policy(limit, gradient.constraints.sharpness),
    );
    let constraint_policy = match gradient.constraints.responses {
        None => constraint_policy,
        Some(options) => {
            let strength = gradient.constraints.stress.map_or_else(
                || "null".into(),
                |limit| strength::policy(limit, gradient.constraints.sharpness),
            );
            format!(
                "{{\"kind\":\"joint-affine-constraints\",\"joint_sharpness\":{:.17e},\"responses\":{},\"strength\":{strength}}}",
                gradient.constraints.sharpness,
                response_output::options(options)
            )
        }
    };
    let check = format!(
        "{{\"goal\":{goal},\"passed\":{},\"directions\":{},\"step\":{:.17e},\"relative_tolerance\":{:.17e},\"baseline_objective\":{:.17e},\"baseline_constraint\":{:.17e},\"objective_relative_error\":{:.17e},\"constraint_relative_error\":{:.17e},\"linear_iterations\":{},\"linear_solves\":{},\"preconditioner_operator_applications\":{},\"preconditioner_galerkin_products\":{},\"constraint_policy\":{constraint_policy}}}",
        gradient.passed(),
        gradient.directions,
        gradient.options.step,
        gradient.options.relative_tolerance,
        gradient.baseline_objective,
        gradient.baseline_constraint,
        gradient.objective_relative_error,
        gradient.constraint_relative_error,
        gradient.work.linear_iterations,
        gradient.work.linear_solves,
        gradient.work.preconditioner_operator_applications,
        gradient.work.preconditioner_galerkin_products
    );
    let evaluation_scope = if state.refinement.is_some() {
        "working-grid optimizer and restored endpoint evaluations; each grid has its declared evaluation allowance, while gradient probes, fresh rechecks and all grids share cumulative physical allowances"
    } else {
        "optimizer and restored endpoint evaluations; initial gradient probes share the original Krylov and wall allowances"
    };
    let optimizer = format!(
        "{{\"iterations\":{},\"evaluations\":{},\"multiplier_updates\":{},\"rejected_trials\":{},\"restoration_evaluations\":{},\"evaluation_budget_scope\":{evaluation_scope:?}}}",
        state.optimizer.iterations,
        state.optimizer.evaluations,
        state.optimizer.multiplier_updates,
        state.optimizer.rejected_trials,
        state.checkpoint.restoration_evaluations
    );
    let policy = config.policy(spec.loads.len())?;
    let oracle = policy.oracle_options();
    let strength_model = config.stress.map_or_else(
        || "null".into(),
        |limit| strength::policy(limit, config.joint_sharpness),
    );
    let response_model = spec
        .affine_responses
        .as_ref()
        .map(response_output::model)
        .transpose()?
        .map_or_else(String::new, |model| {
            format!(",\"response_constraints\":{model}")
        });
    let assessment_model = spec.affine_grid_assessment.map_or_else(String::new, |p| {
        format!(",\"grid_assessment\":{}", grid_assessment_output::policy(p))
    });
    let refinement_model = spec.affine_grid_refinement.map_or_else(String::new, |p| {
        format!(",\"grid_refinement\":{}", grid_refinement_output::policy(p))
    });
    let assessment = match source_selection {
        Some(source) => Some((source.assessment.clone(), source.assessment_html(config)?)),
        None => grid_assessment_output::render(spec, state)?,
    };
    let assessment_json = assessment.as_ref().map_or_else(String::new, |(json, _)| {
        format!(",\"grid_assessment\":{json}")
    });
    let assessment_html = assessment.as_ref().map_or("", |(_, html)| html.as_str());
    let refinement = source
        .as_ref()
        .map(|source| grid_refinement_output::render(spec, study, state, source, &working_selected))
        .transpose()?;
    let refinement_json = refinement.as_ref().map_or_else(String::new, |(json, _)| {
        format!(",\"grid_refinement\":{json}")
    });
    let refinement_html = refinement.as_ref().map_or("", |(_, html)| html.as_str());
    let (compliance_limit, volume_limit, objective_material, constraint_material) = match config.goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit } => (
            "null".into(), format!("{volume_limit:.17e}"), "eroded", "dilated-inventory-and-eroded-physics",
        ),
        AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => (
            format!("{:.17e}", config.cap_j), "null".into(), "dilated", "eroded",
        ),
    };
    let dead_load = dead_load_model(spec);
    let (force_set, dead_offset) = if has_dead_load(spec) {
        ("f=f_dead+B*rho_eroded+Fq, ||q||_2<=1", "authored-reference-load")
    } else {
        ("f=B*rho_eroded+Fq, ||q||_2<=1", "zero")
    };
    let model = format!(
        "{{\"goal\":{goal},\"force_set\":{force_set:?},\"dead_offset\":{dead_offset:?},\"affine_dead_load\":{dead_load},\"basis_columns\":{},\"material_density_kg_m3\":{},\"acceleration_m_s2\":{:?},\"compliance_reference_j\":{},\"compliance_limit_j\":{compliance_limit},\"volume_limit_fraction\":{volume_limit},\"threshold_offset\":{},\"nominal_threshold\":{},\"penal\":{},\"beta\":{},\"constraint_material\":{constraint_material:?},\"objective_material\":{objective_material:?},\"oracle_max_iterations\":{},\"oracle_tolerance\":{},\"oracle_separation_tolerance\":{},\"strength_constraint\":{strength_model},\"joint_sharpness\":{:.17e}{response_model}{assessment_model}{refinement_model}}}",
        spec.loads.len(),
        law.density_kg_m3,
        law.acceleration_m_s2,
        config.cap_j,
        config.options.robustness.threshold_offset,
        spec.schedule[0].eta,
        spec.schedule[0].penal,
        spec.schedule[0].beta,
        oracle.max_iterations,
        oracle.tolerance,
        oracle.separation_tolerance,
        config.joint_sharpness
    );
    let summary = format!(
        "{{\"driver\":{AFFINE_LOAD_SPACE3_DRIVER:?},\"study_id\":{:?},\"status\":{:?},\"optimization_status\":{:?},\"goal\":{goal},\"objective\":{objective:?},\"objective_unit\":{objective_unit:?},\"iterations_completed\":{},\"target_iterations\":{},\"selected_design\":{selection:?},\"selected_feasible\":{feasible},\"selected\":{selected},\"last_accepted\":{accepted},\"model\":{model},\"gradient_check\":{check},\"optimizer_work\":{optimizer},\"work\":{work},\"stop\":{stop}{assessment_json}{refinement_json},\"resume_supported\":true,\"resume_mode\":{:?},\"producer\":{:?},\"authority\":\"Estimated\",\"no_claim\":{}}}",
        spec.id.to_hex(),
        state.status,
        state.optimization_status,
        state.optimizer.iterations,
        spec.updates,
        resume::MODE,
        producer.to_hex(),
        quoted(SCOPE)
    );
    let chosen = &state.selected().evaluation;
    let latest = &state.accepted.evaluation;
    let strength_html = match (&chosen.stress, &latest.stress) {
        (Some(selected), Some(accepted)) => strength::html(selected, accepted),
        _ => String::new(),
    };
    let response_html = match (&spec.affine_responses, &chosen.responses, &latest.responses) {
        (Some(family), Some(selected), Some(accepted)) => {
            response_output::html(family, selected, accepted)?
        }
        _ => String::new(),
    };
    let constraints_html = match source_selection {
        Some(source) => format!(
            "{}<section><h2>Working grid constraint results</h2><p>The following selected and last-accepted values describe the working grid's own optimizer archive, separately from the installed source design.</p>{strength_html}{response_html}</section>",
            source.selected_constraints_html()?
        ),
        None => format!("{strength_html}{response_html}"),
    };
    let selected_volume =
        source_selection.map_or(chosen.physical.dilated_volume_fraction, |s| s.volume);
    let selected_compliance = source_selection.map_or(chosen.worst_compliance, |s| s.compliance);
    let (title, requirement, ratio_label) = match config.goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit } => (
            "3-D minimum worst compliance with self-weight and uncertain loads",
            format!("Dilated volume fraction must be at most {volume_limit:.8e}. Worst compliance is minimized in joules; {:.8e} J is its optimizer normalization reference and is not a physical compliance limit. Every declared signed-response and sampled-strength bound and the joint constraint remain active.", config.cap_j),
            "Worst compliance / reference",
        ),
        AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => (
            "3-D minimum material with self-weight and uncertain loads",
            format!("Dilated material is minimized while worst compliance must be at most {:.8e} J. Every declared signed-response and sampled-strength bound and the joint constraint remain active.", config.cap_j),
            "Worst compliance / limit",
        ),
    };
    let preload_html = spec.affine_dead_load.map_or_else(String::new, |law| {
        let surface = law.surface.map_or_else(|| "none".into(), |surface| match surface.force {
            loading::Force::Pressure(value) => format!(
                "inward-positive pressure {value:.8e} Pa on x-fraction {:?}", surface.x_fraction),
            loading::Force::Traction(value) => format!(
                "traction {value:?} Pa on x-fraction {:?}", surface.x_fraction),
        });
        format!("<p>Common fixed mechanical preload: body force density {:?} N/m³; implicit-surface load: {surface}. This reference-configuration load is added once to the density-dependent self-weight and remains fixed across the uncertain force combinations. Pressure follows the reference normal; follower-load mechanics are not included. The exported offset displacement responds to the complete common preload plus self-weight.</p>", law.body_n_m3)
    });
    let html = format!(
        "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>{title}</title><body><h1>{title}</h1><p>Status: {}. Optimization status: {}. Estimated numerical results.</p><p>Selected design: {selection}; physically feasible for the declared design goal over the complete affine load set: {feasible}.</p><p>The selected design has dilated volume fraction {:.8e} and worst compliance {:.8e} J. The last accepted iterate has dilated volume fraction {:.8e} and worst compliance {:.8e} J; optimizer stopping and KKT diagnostics apply to that iterate.</p><p>{requirement} The common self-weight changes with eroded occupancy; all simultaneous external-force combinations in the declared unit ball are considered.</p>{preload_html}{refinement_html}{constraints_html}{assessment_html}<p>{SCOPE}</p><table><tr><th>Update</th><th>Normalized optimizer objective</th><th>Dilated volume fraction</th><th>Worst compliance (J)</th><th>{ratio_label}</th><th>Joint constraint</th><th>Feasible</th></tr>{table}</table></body></html>",
        state.status,
        state.optimization_status,
        selected_volume,
        selected_compliance,
        latest.physical.dilated_volume_fraction,
        latest.worst_compliance,
    );
    let constellation = hash_bytes(include_bytes!("../../../../../../constellation.lock")).to_hex();
    let mut package = EvidencePackage::new(Provenance::new(
        format!(
            "fs-cli/{}+{AFFINE_LOAD_SPACE3_DRIVER}",
            env!("CARGO_PKG_VERSION")
        ),
        &constellation,
    ));
    let mut claims = vec![
        ("source", spec.canonical.as_str()),
        ("report_json", summary.as_str()),
        ("design", design),
        ("iterations", rows.as_str()),
    ];
    if state.refinement.is_some() {
        claims.push(("working_design", working_design.as_str()));
    }
    for (name, payload) in claims {
        package = package.with_claim(fs_package::Claim::estimated(
            format!("study.sdf3-affine-load-space.{name}"),
            format!(
                "{{\"kind\":{name:?},\"hash\":{:?},\"contents\":{}}}",
                hash_bytes(payload.as_bytes()).to_hex(),
                quoted(payload)
            ),
            AFFINE_LOAD_SPACE3_DRIVER,
            f64::INFINITY,
        ));
    }
    if !fs_checker::check(&package).passed() {
        return Err(fail(
            "cli-study-sdf3-package",
            "structural package check failed",
        ));
    }
    let package = package
        .to_json()
        .map_err(|e| fail("cli-study-sdf3-package", e.to_string()))?;
    let checkpoint = resume::encode(spec, state, producer, &working_design, &rows);
    let mut artifacts = vec![
        ("iterations", "study-iterations", rows.as_bytes()),
        ("design", "study-design", design.as_bytes()),
        ("report_html", "study-report-html", html.as_bytes()),
        ("report_json", "study-report-json", summary.as_bytes()),
        ("package", "study-package", package.as_bytes()),
        ("checkpoint", resume::KIND, checkpoint.as_bytes()),
    ];
    if state.refinement.is_some() {
        artifacts.push(("working_design", "study-design", working_design.as_bytes()));
    }
    if artifacts
        .iter()
        .any(|(_, _, bytes)| bytes.len() as u64 > MAX_ARTIFACT_BYTES)
    {
        return Err(fail(
            "cli-study-sdf3-artifact",
            "affine fields/report exceed the 16 MiB artifact envelope",
        ));
    }
    let versions = format!(
        "{{\"driver\":{AFFINE_LOAD_SPACE3_DRIVER:?},\"crate\":{:?},\"constellation_lock\":{constellation:?},\"producer\":{:?}}}",
        env!("CARGO_PKG_VERSION"),
        producer.to_hex()
    );
    let budgets = format!(
        "{{\"wall_s\":{},\"memory_bytes\":{},\"linear_iterations\":{},\"per_solve_iterations\":{},\"quadrature_boxes\":{},\"quadrature_points\":{},\"max_optimizer_evaluations\":{},\"consumed_wall_s\":{},\"work\":{work},\"optimizer_work\":{optimizer}}}",
        spec.wall_s,
        spec.memory,
        spec.linear,
        spec.per_solve,
        spec.boxes,
        spec.points,
        config.options.optimizer.max_evaluations,
        state.spent.wall_s
    );
    let trace = hash_domain(
        "org.frankensim.cli.sdf3-affine-load-space.trace.v1",
        rows.as_bytes(),
    );
    let predecessor = previous.map_or_else(|| "null".into(), |p| quoted(&p.to_hex()));
    let refinement_context = state
        .refinement
        .as_ref()
        .map_or_else(String::new, |context| {
            format!(",\"grid_refinement\":{}", context.json())
        });
    ledger.begin()?;
    let result = (|| -> Result<Outcome> {
        let op = ledger.begin_op(Some(spec.id.as_bytes()),
            &ir_for(AFFINE_LOAD_SPACE3_DRIVER, spec.id, state.optimizer.iterations),
            &FiveExplicits { seed: &spec.seed.to_le_bytes(), versions: &versions, budget: &budgets,
                capability: "{\"ops\":[\"optimization.marquee-topopt\",\"geometry.sdf\",\"physics.cutfem\"]}" }, 0)?;
        if let Some(previous) = previous {
            ledger.link(op, &previous, EdgeRole::In)?;
        }
        if let Some(context) = &state.refinement {
            if previous != Some(context.source_receipt) {
                ledger.link(op, &context.source_receipt, EdgeRole::In)?;
            }
        }
        let source = ledger.put_artifact("study-source", spec.canonical.as_bytes(), None)?;
        ledger.link(op, &source.hash, EdgeRole::In)?;
        let mut refs = String::new();
        let mut linked_outputs = Vec::with_capacity(artifacts.len());
        for (name, kind, bytes) in artifacts {
            let artifact = ledger.put_artifact(kind, bytes, None)?;
            if !linked_outputs.contains(&artifact.hash) {
                ledger.link(op, &artifact.hash, EdgeRole::Out)?;
                linked_outputs.push(artifact.hash);
            }
            let _ = write!(refs, ",{name:?}:\"{}\"", artifact.hash.to_hex());
        }
        let receipt = format!(
            "{{\"schema\":{STUDY_RUN_RECEIPT_SCHEMA:?},\"driver\":{AFFINE_LOAD_SPACE3_DRIVER:?},\"study_id\":{:?},\"status\":{:?},\"optimization_status\":{:?},\"source\":{:?},\"iterations_completed\":{},\"target_iterations\":{},\"iteration_limit_this_invocation\":{limit},\"trace_hash\":{:?},\"consumed_wall_s\":{},\"work\":{work},\"optimizer_work\":{optimizer},\"selected_design\":{selection:?},\"selected_feasible\":{feasible},\"stop\":{stop},\"predecessor\":{predecessor},\"resume_supported\":true,\"resume_mode\":{:?},\"producer\":{:?}{refinement_context}{refs}}}",
            spec.id.to_hex(),
            state.status,
            state.optimization_status,
            source.hash.to_hex(),
            state.optimizer.iterations,
            spec.updates,
            trace.to_hex(),
            state.spent.wall_s,
            resume::MODE,
            producer.to_hex()
        );
        let stored = ledger.put_artifact(RECEIPT_KIND, receipt.as_bytes(), None)?;
        ledger.link(op, &stored.hash, EdgeRole::Out)?;
        ledger.seal_artifact_output(&stored.hash, op)?;
        ledger.finish_op(op, OpOutcome::Ok, None, 1)?;
        Ok(Outcome {
            pointer: format!("study-{}", stored.hash.to_hex()),
            receipt,
            status: state.status,
        })
    })();
    match result {
        Ok(value) => match ledger.commit() {
            Ok(()) => Ok(value),
            Err(e) => {
                ledger.rollback()?;
                Err(e.into())
            }
        },
        Err(e) => {
            ledger.rollback()?;
            Err(e)
        }
    }
}
