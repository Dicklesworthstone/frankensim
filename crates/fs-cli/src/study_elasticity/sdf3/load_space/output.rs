//! Persist real accepted load-space results through the existing study artifacts.
use super::*;

pub(super) fn measurement(point: &LoadSpaceDesignPoint3, config: Config) -> Result<String> {
    let e = &point.evaluation;
    let feasible = config
        .cap(e.policy.dimension())?
        .admits(e, config.options.optimizer.tolerance);
    Ok(format!(
        "{{\"dilated_volume_fraction\":{:.17e},\"nominal_volume_fraction\":{:.17e},\"eroded_volume_fraction\":{:.17e},\"worst_compliance_j\":{:.17e},\"maximum_ratio\":{:.17e},\"smooth_maximum_ratio\":{:.17e},\"feasible\":{feasible},\"basis_compliances_j\":{:?},\"worst_coefficients\":{:?},\"normalized_gram_row_major\":{:?},\"normalized_eigenvalues\":{:?},\"reciprocity_defect\":{:.17e},\"eigenpair_residual\":{:.17e},\"orthogonality_defect\":{:.17e}}}",
        e.physical.dilated_volume_fraction,
        e.physical.nominal_volume_fraction,
        e.physical.eroded.volume_fraction,
        e.worst_compliance,
        e.spectrum.maximum,
        e.spectrum.smooth_maximum,
        e.physical.eroded.objective.case_compliances,
        e.worst_coefficients,
        e.spectrum.gram,
        e.spectrum.eigenvalues,
        e.spectrum.reciprocity_defect,
        e.spectrum.eigenpair_residual,
        e.spectrum.orthogonality_defect
    ))
}

fn fields(
    study: &CutDensityStudy3<AdaptiveSolveSpace3>,
    point: &LoadSpaceDesignPoint3,
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
            e.physical.eroded.projected_rho[i],
            e.physical.nominal_rho[i],
            e.physical.dilated_rho[i],
            volumes[i]
        );
    }
    let basis = e
        .physical
        .eroded
        .objective
        .displacements
        .iter()
        .map(|field| {
            operator
                .physical_displacements(field)
                .map_err(|e| fail("cli-study-sdf3-field", e.to_string()))
        })
        .collect::<Result<Vec<_>>>()?;
    let worst = operator
        .physical_displacements(&e.worst_displacement)
        .map_err(|e| fail("cli-study-sdf3-field", e.to_string()))?;
    Ok(format!(
        "{{\"cells\":[{cells}],\"physical_nodes_m\":{:?},\"cell_nodes\":{:?},\"displacement_layout\":\"basis column then node-major xyz, on eroded material\",\"basis_displacements_m\":{basis:?},\"worst_displacement_m\":{worst:?}}}",
        operator.physical_nodes(),
        operator.cell_nodes()
    ))
}

fn gradient(state: &State) -> String {
    let check = &state.gradient.numerical;
    let rows = check.probes.iter().map(|p| format!(
        "{{\"direction\":{:?},\"active_densities\":{},\"spectral_analytic\":{:.17e},\"spectral_difference\":{:.17e},\"spectral_relative_error\":{:.17e},\"material_analytic\":{:.17e},\"material_difference\":{:.17e},\"material_relative_error\":{:.17e}}}",
        format!("{:?}", p.direction), p.active_densities, p.compliance_analytic,
        p.compliance_difference, p.compliance_relative_error, p.volume_analytic,
        p.volume_difference, p.volume_relative_error)).collect::<Vec<_>>().join(",");
    format!(
        "{{\"passed\":{},\"step\":{},\"relative_tolerance\":{},\"evaluations\":{},\"baseline_spectral\":{:.17e},\"baseline_material\":{:.17e},\"work\":{{\"linear_iterations\":{},\"linear_solves\":{},\"preconditioner_operator_applications\":{},\"preconditioner_galerkin_products\":{}}},\"probes\":[{rows}]}}",
        check.passed(),
        check.options.step,
        check.options.relative_tolerance,
        1 + 2 * check.probes.len(),
        check.baseline_compliance, check.baseline_volume_fraction,
        check.work.linear_iterations, check.work.linear_solves,
        check.work.preconditioner_operator_applications, check.work.preconditioner_galerkin_products,
    )
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
    let config = spec.load_space.expect("load-space producer admission");
    let cap = config.cap(spec.loads.len())?;
    let selection = if state.best.is_some() {
        "least-material-feasible"
    } else {
        "last-accepted"
    };
    let feasible = cap.admits(
        &state.selected().evaluation,
        config.options.optimizer.tolerance,
    );
    let selected = measurement(state.selected(), config)?;
    let accepted = measurement(&state.accepted, config)?;
    let design = format!(
        "{{\"schema\":\"sdf3-load-space-design.v1\",\"selected_design\":{selection:?},\"selected_feasible\":{feasible},\"selected\":{},\"last_accepted\":{}}}",
        fields(study, state.selected())?,
        fields(study, &state.accepted)?
    );
    let mut table = String::new();
    let mut history = String::new();
    for (i, row) in state.history.iter().enumerate() {
        if i != 0 {
            history.push(',');
        }
        let _ = write!(
            history,
            "{{\"iteration\":{},\"dilated_volume_fraction\":{:.17e},\"maximum_ratio\":{:.17e},\"smooth_maximum_ratio\":{:.17e},\"feasible\":{}}}",
            row.iteration, row.dilated_volume, row.maximum, row.smooth_maximum, row.feasible
        );
        let _ = write!(
            table,
            "<tr><td>{}</td><td>{:.8e}</td><td>{:.8e}</td><td>{}</td></tr>",
            row.iteration, row.dilated_volume, row.maximum, row.feasible
        );
    }
    let rows = format!("{{\"schema\":\"sdf3-load-space-iterations.v1\",\"history\":[{history}]}}");
    let work = state.spent.json();
    let stop = stopping(state);
    let check = gradient(state);
    let optimizer = format!(
        "{{\"iterations\":{},\"evaluations\":{},\"multiplier_updates\":{},\"rejected_trials\":{},\"gradient_probe_evaluations\":{},\"restoration_evaluations\":{},\"evaluation_budget_scope\":\"optimizer and restoration evaluations; gradient probes share the original Krylov and wall budgets\"}}",
        state.optimizer.iterations,
        state.optimizer.evaluations,
        state.optimizer.multiplier_updates,
        state.optimizer.rejected_trials,
        1 + 2 * state.gradient.numerical.probes.len(),
        state.checkpoint.restoration_evaluations,
    );
    let model = format!(
        "{{\"force_set\":\"f=Fq, ||q||_2<=1\",\"basis_columns\":{},\"compliance_limit_j\":{},\"spectral_temperature\":{},\"smooth_bias_bound\":{},\"threshold_offset\":{},\"nominal_threshold\":{},\"penal\":{},\"beta\":{},\"constraint_material\":\"eroded\",\"objective_material\":\"dilated\"}}",
        spec.loads.len(),
        config.cap_j,
        config.temperature,
        cap.policy().bias_bound(),
        config.options.robustness.threshold_offset,
        spec.schedule[0].eta,
        spec.schedule[0].penal,
        spec.schedule[0].beta
    );
    let summary = format!(
        "{{\"driver\":{LOAD_SPACE3_DRIVER:?},\"study_id\":{:?},\"status\":{:?},\"objective\":\"dilated-volume-fraction\",\"objective_unit\":\"1\",\"iterations_completed\":{},\"target_iterations\":{},\"selected_design\":{selection:?},\"selected_feasible\":{feasible},\"selected\":{selected},\"last_accepted\":{accepted},\"model\":{model},\"gradient_check\":{check},\"optimizer_work\":{optimizer},\"work\":{work},\"stop\":{stop},\"resume_supported\":true,\"resume_mode\":{:?},\"producer\":{:?},\"authority\":\"Estimated\",\"no_claim\":{}}}",
        spec.id.to_hex(),
        state.status,
        state.optimizer.iterations,
        spec.updates,
        resume::MODE, producer.to_hex(),
        quoted(SCOPE)
    );
    let html = format!(
        "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>3-D robust minimum-material study</title><body><h1>3-D robust minimum-material study</h1><p>Status: {}. Estimated numerical results.</p><p>Selected design: {selection}; compliant over the declared load set: {feasible}.</p><p>Every simultaneous combination f=Fq with ||q||₂ ≤ 1 must satisfy compliance ≤ {:.8e} J on eroded material. The minimized material fraction uses the dilated field.</p><p>{SCOPE}</p><table><tr><th>Update</th><th>Dilated material fraction</th><th>Worst compliance / limit</th><th>Feasible</th></tr>{table}</table></body></html>",
        state.status, config.cap_j
    );
    let constellation = hash_bytes(include_bytes!("../../../../../../constellation.lock")).to_hex();
    let mut package = EvidencePackage::new(Provenance::new(
        format!("fs-cli/{}+{LOAD_SPACE3_DRIVER}", env!("CARGO_PKG_VERSION")),
        &constellation,
    ));
    for (name, payload) in [
        ("source", spec.canonical.as_str()),
        ("report_json", summary.as_str()),
        ("design", design.as_str()),
        ("iterations", rows.as_str()),
    ] {
        package = package.with_claim(fs_package::Claim::estimated(
            format!("study.sdf3-load-space.{name}"),
            format!(
                "{{\"kind\":{name:?},\"hash\":{:?},\"contents\":{}}}",
                hash_bytes(payload.as_bytes()).to_hex(),
                quoted(payload)
            ),
            LOAD_SPACE3_DRIVER,
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
            "load-space fields/report exceed the 16 MiB artifact envelope",
        ));
    }
    let versions = format!(
        "{{\"driver\":{LOAD_SPACE3_DRIVER:?},\"crate\":{:?},\"constellation_lock\":{constellation:?},\"producer\":{:?}}}",
        env!("CARGO_PKG_VERSION"), producer.to_hex()
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
        "org.frankensim.cli.sdf3-load-space.trace.v1",
        rows.as_bytes(),
    );
    let predecessor = previous.map_or_else(|| "null".into(), |p| quoted(&p.to_hex()));
    ledger.begin()?;
    let result = (|| -> Result<Outcome> {
        let op = ledger.begin_op(Some(spec.id.as_bytes()), &ir_for(LOAD_SPACE3_DRIVER, spec.id, state.optimizer.iterations),
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
            "{{\"schema\":{STUDY_RUN_RECEIPT_SCHEMA:?},\"driver\":{LOAD_SPACE3_DRIVER:?},\"study_id\":{:?},\"status\":{:?},\"source\":{:?},\"iterations_completed\":{},\"target_iterations\":{},\"iteration_limit_this_invocation\":{limit},\"trace_hash\":{:?},\"consumed_wall_s\":{},\"work\":{work},\"optimizer_work\":{optimizer},\"selected_design\":{selection:?},\"selected_feasible\":{feasible},\"stop\":{stop},\"predecessor\":{predecessor},\"resume_supported\":true,\"resume_mode\":{:?},\"producer\":{:?}{refs}}}",
            spec.id.to_hex(),
            state.status,
            source.hash.to_hex(),
            state.optimizer.iterations,
            spec.updates,
            trace.to_hex(),
            state.spent.wall_s,
            resume::MODE, producer.to_hex()
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
