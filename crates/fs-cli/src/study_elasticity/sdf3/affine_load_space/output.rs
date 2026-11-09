//! Retain the actual affine witness and selected physical design in the ledger.
use super::*;

pub(super) fn measurement(point: &AffineLoadSpaceDesignPoint3, config: Config) -> String {
    let e = &point.evaluation;
    let p = &e.physical;
    let o = &e.oracle;
    format!(
        "{{\"dilated_volume_fraction\":{:.17e},\"nominal_volume_fraction\":{:.17e},\"eroded_volume_fraction\":{:.17e},\"worst_compliance_j\":{:.17e},\"maximum_ratio\":{:.17e},\"feasible\":{},\"worst_coefficients\":{:?},\"normalized_gram_row_major\":{:?},\"normalized_linear\":{:?},\"normalized_offset\":{:.17e},\"reciprocity_defect\":{:.17e},\"witness_relative_defect\":{:.17e},\"witness_relative_residual\":{:.17e},\"oracle\":{{\"branch\":{:?},\"gradient_admissible\":{},\"value_without_offset\":{:.17e},\"dual_value_without_offset\":{:.17e},\"duality_gap\":{:.17e},\"relative_duality_gap\":{:.17e},\"relative_separation\":{:.17e},\"stationarity_residual\":{:.17e},\"feasibility_residual\":{:.17e},\"complementarity_residual\":{:.17e},\"eigenpair_residual\":{:.17e},\"iterations\":{}}}}}",
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
    )
}

fn fields(
    study: &CutDensityStudy3<AdaptiveSolveSpace3>,
    point: &AffineLoadSpaceDesignPoint3,
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
    Ok(format!(
        "{{\"cells\":[{cells}],\"physical_nodes_m\":{:?},\"cell_nodes\":{:?},\"displacement_layout\":\"basis column then node-major xyz, on eroded material\",\"basis_displacements_m\":{basis:?},\"offset_displacement_m\":{offset:?},\"worst_displacement_m\":{worst:?}}}",
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
    let selection = if state.best.is_some() {
        "least-material-feasible"
    } else {
        "last-accepted"
    };
    let feasible = config.feasible(state.selected());
    let selected = measurement(state.selected(), config);
    let accepted = measurement(&state.accepted, config);
    let design = format!(
        "{{\"schema\":\"sdf3-affine-load-space-design.v1\",\"selected_design\":{selection:?},\"selected_feasible\":{feasible},\"selected\":{},\"last_accepted\":{}}}",
        fields(study, state.selected())?,
        fields(study, &state.accepted)?
    );
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
            "<tr><td>{}</td><td>{:.8e}</td><td>{:.8e}</td><td>{}</td></tr>",
            row.iteration, row.dilated_volume_fraction, row.maximum, row.feasible
        );
    }
    let rows =
        format!("{{\"schema\":\"sdf3-affine-load-space-iterations.v1\",\"history\":[{history}]}}");
    let work = state.spent.json();
    let stop = stopping(state);
    let gradient = &state.gradient;
    let check = format!(
        "{{\"passed\":{},\"directions\":{},\"step\":{:.17e},\"relative_tolerance\":{:.17e},\"baseline_objective\":{:.17e},\"baseline_constraint\":{:.17e},\"objective_relative_error\":{:.17e},\"constraint_relative_error\":{:.17e},\"linear_iterations\":{},\"linear_solves\":{},\"preconditioner_operator_applications\":{},\"preconditioner_galerkin_products\":{}}}",
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
    let optimizer = format!(
        "{{\"iterations\":{},\"evaluations\":{},\"multiplier_updates\":{},\"rejected_trials\":{},\"restoration_evaluations\":{},\"evaluation_budget_scope\":\"optimizer and restored endpoint evaluations; initial gradient probes share the original Krylov and wall allowances\"}}",
        state.optimizer.iterations,
        state.optimizer.evaluations,
        state.optimizer.multiplier_updates,
        state.optimizer.rejected_trials,
        state.checkpoint.restoration_evaluations
    );
    let policy = config.policy(spec.loads.len())?;
    let oracle = policy.oracle_options();
    let model = format!(
        "{{\"force_set\":\"f=B*rho_eroded+Fq, ||q||_2<=1\",\"dead_offset\":\"zero\",\"basis_columns\":{},\"material_density_kg_m3\":{},\"acceleration_m_s2\":{:?},\"compliance_limit_j\":{},\"threshold_offset\":{},\"nominal_threshold\":{},\"penal\":{},\"beta\":{},\"constraint_material\":\"eroded\",\"objective_material\":\"dilated\",\"oracle_max_iterations\":{},\"oracle_tolerance\":{},\"oracle_separation_tolerance\":{}}}",
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
        oracle.separation_tolerance
    );
    let summary = format!(
        "{{\"driver\":{AFFINE_LOAD_SPACE3_DRIVER:?},\"study_id\":{:?},\"status\":{:?},\"objective\":\"dilated-volume-fraction\",\"objective_unit\":\"1\",\"iterations_completed\":{},\"target_iterations\":{},\"selected_design\":{selection:?},\"selected_feasible\":{feasible},\"selected\":{selected},\"last_accepted\":{accepted},\"model\":{model},\"gradient_check\":{check},\"optimizer_work\":{optimizer},\"work\":{work},\"stop\":{stop},\"resume_supported\":true,\"resume_mode\":{:?},\"producer\":{:?},\"authority\":\"Estimated\",\"no_claim\":{}}}",
        spec.id.to_hex(),
        state.status,
        state.optimizer.iterations,
        spec.updates,
        resume::MODE,
        producer.to_hex(),
        quoted(SCOPE)
    );
    let chosen = &state.selected().evaluation;
    let latest = &state.accepted.evaluation;
    let html = format!(
        "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>3-D minimum material with self-weight and uncertain loads</title><body><h1>3-D minimum material with self-weight and uncertain loads</h1><p>Status: {}. Estimated numerical results.</p><p>Selected design: {selection}; cap-feasible over the declared affine load set: {feasible}.</p><p>The selected design has dilated volume fraction {:.8e} and worst compliance {:.8e} J. The last accepted iterate has dilated volume fraction {:.8e} and worst compliance {:.8e} J; optimizer stopping and KKT diagnostics apply to that iterate.</p><p>Worst compliance must be at most {:.8e} J. The common self-weight changes with eroded occupancy; all simultaneous external-force combinations in the declared unit ball are considered.</p><p>{SCOPE}</p><table><tr><th>Update</th><th>Dilated volume fraction</th><th>Worst compliance / limit</th><th>Feasible</th></tr>{table}</table></body></html>",
        state.status,
        chosen.physical.dilated_volume_fraction,
        chosen.worst_compliance,
        latest.physical.dilated_volume_fraction,
        latest.worst_compliance,
        config.cap_j
    );
    let constellation = hash_bytes(include_bytes!("../../../../../../constellation.lock")).to_hex();
    let mut package = EvidencePackage::new(Provenance::new(
        format!(
            "fs-cli/{}+{AFFINE_LOAD_SPACE3_DRIVER}",
            env!("CARGO_PKG_VERSION")
        ),
        &constellation,
    ));
    for (name, payload) in [
        ("source", spec.canonical.as_str()),
        ("report_json", summary.as_str()),
        ("design", design.as_str()),
        ("iterations", rows.as_str()),
    ] {
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
    let checkpoint = resume::encode(spec, state, producer, &design, &rows);
    let artifacts = [
        ("iterations", "study-iterations", rows.as_bytes()),
        ("design", "study-design", design.as_bytes()),
        ("report_html", "study-report-html", html.as_bytes()),
        ("report_json", "study-report-json", summary.as_bytes()),
        ("package", "study-package", package.as_bytes()),
        ("checkpoint", resume::KIND, checkpoint.as_bytes()),
    ];
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
    ledger.begin()?;
    let result = (|| -> Result<Outcome> {
        let op = ledger.begin_op(Some(spec.id.as_bytes()),
            &ir_for(AFFINE_LOAD_SPACE3_DRIVER, spec.id, state.optimizer.iterations),
            &FiveExplicits { seed: &spec.seed.to_le_bytes(), versions: &versions, budget: &budgets,
                capability: "{\"ops\":[\"optimization.marquee-topopt\",\"geometry.sdf\",\"physics.cutfem\"]}" }, 0)?;
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
        let receipt = format!(
            "{{\"schema\":{STUDY_RUN_RECEIPT_SCHEMA:?},\"driver\":{AFFINE_LOAD_SPACE3_DRIVER:?},\"study_id\":{:?},\"status\":{:?},\"source\":{:?},\"iterations_completed\":{},\"target_iterations\":{},\"iteration_limit_this_invocation\":{limit},\"trace_hash\":{:?},\"consumed_wall_s\":{},\"work\":{work},\"optimizer_work\":{optimizer},\"selected_design\":{selection:?},\"selected_feasible\":{feasible},\"stop\":{stop},\"predecessor\":{predecessor},\"resume_supported\":true,\"resume_mode\":{:?},\"producer\":{:?}{refs}}}",
            spec.id.to_hex(),
            state.status,
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
