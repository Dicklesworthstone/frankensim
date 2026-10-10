//! Sealed exports of the accepted and best-feasible numerical designs. No solve.
use super::*;

fn measurement(value: Option<&StressEvaluation3>, options: StressDesignOptions3) -> String {
    value.map_or_else(|| "null".into(), |e| format!(
        "{{\"volume_fraction\":{:.17e},\"stress_aggregate_pa\":{:.17e},\"sampled_relaxed_max_pa\":{:.17e},\"sampled_physical_max_pa\":{:.17e},\"constraint_violation\":{:.17e},\"feasible\":{},\"case_compliances_j\":{:?},\"case_relaxed_max_pa\":{:?},\"case_physical_max_pa\":{:?},\"normalized_load_weights\":{:?},\"stress_points\":{}}}",
        e.volume_fraction,e.aggregate,e.sampled_relaxed_max,e.sampled_physical_max,
        (e.aggregate/options.stress_limit-1.0).max(0.0),
        e.aggregate/options.stress_limit-1.0 <= options.optimizer.tolerance,
        e.case_compliances,e.case_relaxed_max,e.case_physical_max,e.normalized_load_weights,e.point_count))
}

fn fields(
    study: &CutDensityStudy3<AdaptiveSolveSpace3>,
    value: Option<&StressEvaluation3>,
) -> Result<String> {
    let Some(value) = value else {
        return Ok("null".into());
    };
    if value
        .rho
        .iter()
        .chain(&value.projected_rho)
        .chain(value.displacements.iter().flatten())
        .any(|v| !v.is_finite())
    {
        return Err(fail(
            "cli-study-sdf3-field",
            "nonfinite retained stress field",
        ));
    }
    let op = study.operator().elasticity();
    let mut cells = String::new();
    for (i, cell) in op.leaves().iter().enumerate() {
        if i != 0 {
            cells.push(',');
        }
        let _ = write!(
            cells,
            "{{\"level\":{},\"index\":{:?},\"raw_density\":{:.17e},\"projected_density\":{:.17e},\"cut_volume_m3\":{:.17e}}}",
            cell.level(),
            cell.index(),
            value.rho[i],
            value.projected_rho[i],
            op.volumes()[i]
        );
    }
    let mut displacements = String::new();
    for (i, field) in value.displacements.iter().enumerate() {
        if i != 0 {
            displacements.push(',');
        }
        let physical = op
            .physical_displacements(field)
            .map_err(|e| fail("cli-study-sdf3-field", e.to_string()))?;
        let _ = write!(displacements, "{physical:?}");
    }
    Ok(format!(
        "{{\"cells\":[{cells}],\"physical_nodes_m\":{:?},\"cell_nodes\":{:?},\"displacement_layout\":\"load-case then node-major xyz\",\"displacements_m\":[{displacements}],\"cell_relaxed_max_pa\":{:?}}}",
        op.physical_nodes(),
        op.cell_nodes(),
        value.cell_relaxed_max
    ))
}

fn rows(state: &State) -> (String, String) {
    let mut rows = String::new();
    let mut html = String::new();
    for (i, row) in state.history.iter().enumerate() {
        if i != 0 {
            rows.push(',');
        }
        let _ = write!(
            rows,
            "{{\"iteration\":{},\"volume_fraction\":{:.17e},\"stress_aggregate_pa\":{:.17e},\"sampled_relaxed_max_pa\":{:.17e},\"sampled_physical_max_pa\":{:.17e},\"constraint_violation\":{:.17e},\"feasible\":{}}}",
            row.iteration,
            row.volume_fraction,
            row.stress_aggregate,
            row.sampled_relaxed_max,
            row.sampled_physical_max,
            row.constraint_violation,
            row.feasible
        );
        let _ = write!(
            html,
            "<tr><td>{}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td><td>{}</td></tr>",
            row.iteration,
            row.volume_fraction,
            row.stress_aggregate,
            row.constraint_violation,
            row.feasible
        );
    }
    (
        format!("{{\"schema\":\"sdf3-stress-iterations.v1\",\"history\":[{rows}]}}"),
        html,
    )
}

fn probes(state: &State) -> String {
    let mut rows = String::new();
    for (i, p) in state.audit.probes.iter().enumerate() {
        if i != 0 {
            rows.push(',');
        }
        let _ = write!(
            rows,
            "{{\"direction\":{:?},\"active_densities\":{},\"stress_analytic\":{:.17e},\"stress_difference\":{:.17e},\"stress_relative_error\":{:.17e},\"volume_analytic\":{:.17e},\"volume_difference\":{:.17e},\"volume_relative_error\":{:.17e}}}",
            p.direction,
            p.active,
            p.stress_analytic,
            p.stress_difference,
            p.stress_relative_error,
            p.volume_analytic,
            p.volume_difference,
            p.volume_relative_error
        );
    }
    format!(
        "{{\"passed\":{},\"step\":1e-4,\"relative_tolerance\":5e-4,\"evaluations\":{},\"probes\":[{rows}]}}",
        state.audit.passed, state.audit.evaluations
    )
}

fn stopping(state: &State) -> String {
    if let Some(report) = &state.report {
        format!(
            "{{\"reason\":{:?},\"kkt_scope\":\"last-accepted\",\"kkt\":{{\"stationarity\":{:.17e},\"feasibility\":{:.17e},\"dual_feasibility\":{:.17e},\"complementarity\":{:.17e}}},\"multiplier\":{:.17e},\"penalty\":{:.17e}}}",
            format!("{:?}", report.stop),
            report.kkt.stationarity,
            report.kkt.feasibility,
            report.kkt.dual_feasibility,
            report.kkt.complementarity,
            report.multiplier,
            report.penalty
        )
    } else {
        format!(
            "{{\"reason\":\"EvaluationStopped\",\"error\":{},\"kkt\":null}}",
            state
                .error
                .as_ref()
                .map_or_else(|| "null".into(), |e| quoted(&e.to_string()))
        )
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
    let options = spec.stress.expect("stress producer admission");
    let selection = if state.best.is_some() {
        "best-feasible"
    } else if state.accepted.is_some() {
        "last-accepted"
    } else {
        "none"
    };
    let selected = measurement(state.selected(), options);
    let last = measurement(state.accepted.as_ref(), options);
    let final_volume = state
        .selected()
        .map_or_else(|| "null".into(), |e| format!("{:.17e}", e.volume_fraction));
    let feasible = state
        .selected()
        .is_some_and(|e| e.aggregate / options.stress_limit - 1.0 <= options.optimizer.tolerance);
    let design = format!(
        "{{\"schema\":\"sdf3-stress-design.v1\",\"state\":{:?},\"selected_design\":{selection:?},\"selected_feasible\":{feasible},\"penal\":{},\"beta\":{},\"selected\":{},\"last_accepted\":{}}}",
        if state.accepted.is_some() {
            "solved"
        } else {
            "no-solved-baseline"
        },
        spec.schedule[0].penal,
        spec.schedule[0].beta,
        fields(study, state.selected())?,
        fields(study, state.accepted.as_ref())?
    );
    let (rows, table) = rows(state);
    let checkpoint = resume::encode(spec, state, producer, &design, &rows);
    let resumable = checkpoint.is_some();
    let restoration_evaluations = state.checkpoint.as_ref().map_or(0, |c| c.restoration_evaluations);
    let trace = hash_domain("org.frankensim.cli.sdf3-stress.trace.v1", rows.as_bytes());
    let stop = stopping(state);
    let gradient = probes(state);
    let work = state.spent.json();
    let evaluations = state.audit.evaluations + state.optimizer_work.evaluations;
    let optimizer = format!(
        "{{\"iterations\":{},\"evaluations\":{},\"multiplier_updates\":{},\"rejected_trials\":{},\"restoration_evaluations\":{restoration_evaluations},\"total_evaluations_including_gradient_gate\":{evaluations}}}",
        state.iterations(),
        state.optimizer_work.evaluations,
        state.optimizer_work.multiplier_updates,
        state.optimizer_work.rejected_trials
    );
    let measure = format!(
        "{{\"kind\":\"normalized-volume-and-load-weighted-qp-von-mises\",\"stress_limit_pa\":{:.17e},\"relaxation_power\":{},\"aggregation_power\":{},\"maximum_stress_is_constrained\":false}}",
        options.stress_limit, options.stress.relaxation_power, options.stress.aggregation_power
    );
    let summary = format!(
        "{{\"driver\":{STRESS3_DRIVER:?},\"study_id\":{:?},\"status\":{:?},\"objective\":\"volume-fraction\",\"objective_unit\":\"1\",\"final_volume_fraction\":{final_volume},\"iterations_completed\":{},\"target_iterations\":{},\"selected_design\":{selection:?},\"selected_feasible\":{feasible},\"selected\":{selected},\"last_accepted\":{last},\"stress_measure\":{measure},\"gradient_check\":{gradient},\"optimizer_work\":{optimizer},\"work\":{work},\"stop\":{stop},\"resume_supported\":{resumable},\"resume_mode\":{:?},\"authority\":\"Estimated\",\"no_claim\":{}}}",
        spec.id.to_hex(),
        state.status,
        state.iterations(),
        spec.updates,
        resume::MODE,
        quoted(SCOPE)
    );
    let html = format!(
        "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>3-D minimum-volume stress study</title><body><h1>3-D minimum-volume stress study</h1><p>Status: {}. Estimated numerical evidence.</p><p>Selected design: {selection}; aggregate-feasible: {feasible}; projected material fraction: {final_volume}.</p><p>Constraint: normalized qp von Mises aggregate ≤ {:.8e} Pa, q = {}, p = {}. This is not a cap on the sampled or continuum maximum stress.</p><p>{SCOPE}</p><table><tr><th>Update</th><th>Material fraction</th><th>Stress aggregate Pa</th><th>Relative violation</th><th>Feasible</th></tr>{table}</table></body></html>",
        state.status,
        options.stress_limit,
        options.stress.relaxation_power,
        options.stress.aggregation_power
    );
    let constellation = hash_bytes(include_bytes!("../../../../../../constellation.lock")).to_hex();
    let mut package = EvidencePackage::new(Provenance::new(
        format!("fs-cli/{}+{STRESS3_DRIVER}", env!("CARGO_PKG_VERSION")),
        &constellation,
    ));
    // Preserve the actual numerical results in the standalone package, not
    // merely its provenance. Estimated source statements bind and carry the
    // exact artifact bytes; structural checking does not certify mechanics.
    for (name, payload) in [
        ("source", spec.canonical.as_str()),
        ("report_json", summary.as_str()),
        ("design", design.as_str()),
        ("iterations", rows.as_str()),
    ] {
        package = package.with_claim(fs_package::Claim::estimated(
            format!("study.sdf3-stress.{name}"),
            format!(
                "{{\"schema\":\"sdf3-stress-package-artifact.v1\",\"kind\":{name:?},\"hash\":{:?},\"contents\":{}}}",
                hash_bytes(payload.as_bytes()).to_hex(),
                quoted(payload)
            ),
            STRESS3_DRIVER,
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
    let versions = format!(
        "{{\"driver\":{STRESS3_DRIVER:?},\"crate\":{:?},\"constellation_lock\":{constellation:?},\"producer\":{:?},\"policy\":\"PHR projected AL v1 defaults except explicit input; initial-gradient-step=1e-4; gradient-rtol=5e-4; quadrature-depth=2; fixed-background; accepted optimizer-state restoration v1\"}}",
        env!("CARGO_PKG_VERSION"), producer.to_hex()
    );
    let budgets = format!(
        "{{\"wall_s\":{},\"consumed_wall_s\":{},\"memory_bytes\":{},\"linear_iterations\":{},\"per_solve_iterations\":{},\"quadrature_boxes\":{},\"quadrature_points\":{},\"max_stress_points_per_evaluation\":{},\"max_evaluations_including_gradient_gate\":{},\"work\":{work},\"optimizer_work\":{optimizer}}}",
        spec.wall_s,
        state.spent.wall_s,
        spec.memory,
        spec.linear,
        spec.per_solve,
        spec.boxes,
        spec.points,
        options.stress.max_points,
        options.optimizer.max_evaluations
    );
    let mut artifacts = vec![
        ("iterations", "study-iterations", rows.as_bytes()),
        ("design", "study-design", design.as_bytes()),
        ("report_html", "study-report-html", html.as_bytes()),
        ("report_json", "study-report-json", summary.as_bytes()),
        ("package", "study-package", package.as_bytes()),
    ];
    if let Some(checkpoint) = &checkpoint {
        artifacts.push(("checkpoint", resume::KIND, checkpoint.as_bytes()));
    }
    if artifacts
        .iter()
        .any(|(_, _, bytes)| bytes.len() as u64 > MAX_ARTIFACT_BYTES)
    {
        return Err(fail(
            "cli-study-sdf3-artifact",
            "stress design/report exceeds the 16 MiB artifact envelope",
        ));
    }
    let predecessor = previous.map_or_else(|| "null".into(), |p| quoted(&p.to_hex()));
    ledger.begin()?;
    let result = (|| -> Result<Outcome> {
        let op = ledger.begin_op(Some(spec.id.as_bytes()),&ir_for(STRESS3_DRIVER,spec.id,state.iterations()),&FiveExplicits {
            seed: &spec.seed.to_le_bytes(), versions: &versions,budget: &budgets,
            capability: "{\"ops\":[\"optimization.marquee-topopt\",\"geometry.sdf\",\"physics.cutfem\"]}" },0)?;
        if let Some(previous) = previous {
            ledger.link(op, &previous, EdgeRole::In)?;
        }
        let source = ledger.put_artifact("study-source", spec.canonical.as_bytes(), None)?;
        ledger.link(op, &source.hash, EdgeRole::In)?;
        let mut refs = String::new();
        for (name, kind, bytes) in artifacts {
            let artifact = ledger.put_artifact(kind, bytes, None)?;
            ledger.link(op, &artifact.hash, EdgeRole::Out)?;
            let _ = write!(refs, ",{name:?}:\"{}\"", artifact.hash.to_hex());
        }
        if !resumable { refs.push_str(",\"checkpoint\":null"); }
        let receipt = format!(
            "{{\"schema\":{STUDY_RUN_RECEIPT_SCHEMA:?},\"driver\":{STRESS3_DRIVER:?},\"study_id\":{:?},\"status\":{:?},\"source\":{:?},\"iterations_completed\":{},\"target_iterations\":{},\"iteration_limit_this_invocation\":{limit},\"trace_hash\":{:?},\"consumed_wall_s\":{},\"work\":{work},\"optimizer_work\":{optimizer},\"selected_design\":{selection:?},\"selected_feasible\":{feasible},\"stop\":{stop},\"predecessor\":{predecessor},\"producer\":{:?},\"resume_supported\":{resumable},\"resume_mode\":{:?}{refs}}}",
            spec.id.to_hex(),
            state.status,
            source.hash.to_hex(),
            state.iterations(),
            spec.updates,
            trace.to_hex(),
            state.spent.wall_s,
            producer.to_hex(), resume::MODE,
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
        Ok(result) => match ledger.commit() {
            Ok(()) => Ok(result),
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
