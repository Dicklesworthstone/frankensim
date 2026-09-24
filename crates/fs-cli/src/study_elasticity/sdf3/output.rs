//! Retained numerical fields and stage-local evidence; export performs no solve.
use super::*;

fn design(run: &Computation) -> Result<String> {
    let Some(last) = &run.report.continuation.last else {
        return Ok("{\"schema\":\"sdf3-study-design.v1\",\"state\":\"no-solved-baseline\",\"cells\":[],\"displacements_m\":[]}".into());
    };
    let operator = run.study.operator().elasticity();
    let mut cells = String::new();
    let volumes = operator.volumes();
    for (index, cell) in operator.leaves().iter().enumerate() {
        if index != 0 {
            cells.push(',');
        }
        let _ = write!(
            cells,
            "{{\"level\":{},\"index\":{:?},\"raw_density\":{:.17e},\"projected_density\":{:.17e},\"cut_volume_m3\":{:.17e}}}",
            cell.level(),
            cell.index(),
            last.rho[index],
            last.projected_rho[index],
            volumes[index]
        );
    }
    let mut fields = String::new();
    for (index, field) in last.displacements.iter().enumerate() {
        if index != 0 {
            fields.push(',');
        }
        let physical = operator
            .physical_displacements(field)
            .map_err(|e| fail("cli-study-sdf3-field", e.to_string()))?;
        let _ = write!(fields, "{physical:?}");
    }
    Ok(format!(
        "{{\"schema\":\"sdf3-study-design.v1\",\"state\":\"solved\",\"penal\":{},\"beta\":{},\"cells\":[{cells}],\"physical_nodes_m\":{:?},\"cell_nodes\":{:?},\"displacement_layout\":\"load-case then node-major xyz\",\"displacements_m\":[{fields}]}}",
        run.report.continuation.params.penal,
        run.report.continuation.params.beta,
        operator.physical_nodes(),
        operator.cell_nodes()
    ))
}

fn stages(run: &Computation) -> (String, String, usize) {
    let mut output = String::new();
    let mut html = String::new();
    let mut updates = 0;
    for (ordinal, stage) in run.report.continuation.stages.iter().enumerate() {
        if ordinal != 0 {
            output.push(',');
        }
        let mut rows = String::new();
        for (index, row) in stage.history.iter().enumerate() {
            if index != 0 {
                rows.push(',');
            }
            let _ = write!(
                rows,
                "{{\"iteration\":{},\"compliance_j\":{:.17e},\"case_compliances_j\":{:?},\"volume_fraction\":{:.17e},\"max_density_change\":{:.17e}}}",
                row.iteration,
                row.compliance,
                row.case_compliances,
                row.volume_fraction,
                row.max_change
            );
            let _ = write!(
                html,
                "<tr><td>{}</td><td>{}</td><td>{:.8e}</td><td>{:.8e}</td></tr>",
                stage.stage, row.iteration, row.compliance, row.volume_fraction
            );
        }
        updates += stage.history.len().saturating_sub(1);
        let mut probes = String::new();
        let passed = stage
            .gradient_check
            .as_ref()
            .is_some_and(|check| check.passed());
        if let Some(check) = &stage.gradient_check {
            for (index, probe) in check.probes.iter().enumerate() {
                if index != 0 {
                    probes.push(',');
                }
                let _ = write!(
                    probes,
                    "{{\"direction\":\"{:?}\",\"active_densities\":{},\"compliance_analytic\":{:.17e},\"compliance_difference\":{:.17e},\"compliance_relative_error\":{:.17e},\"volume_analytic\":{:.17e},\"volume_difference\":{:.17e},\"volume_relative_error\":{:.17e}}}",
                    probe.direction,
                    probe.active_densities,
                    probe.compliance_analytic,
                    probe.compliance_difference,
                    probe.compliance_relative_error,
                    probe.volume_analytic,
                    probe.volume_difference,
                    probe.volume_relative_error
                );
            }
        }
        let _ = write!(
            output,
            "{{\"stage\":{},\"penal\":{},\"beta\":{},\"incoming_volume_fraction\":{:.17e},\"restoration_scale\":{:.17e},\"gradient_check_passed\":{passed},\"gradient_probes\":[{probes}],\"termination\":\"{:?}\",\"history\":[{rows}]}}",
            stage.stage,
            stage.params.penal,
            stage.params.beta,
            stage.incoming_volume_fraction,
            stage.restoration_scale,
            stage.termination
        );
    }
    (
        format!("{{\"schema\":\"sdf3-study-stages.v1\",\"stages\":[{output}]}}"),
        html,
        updates,
    )
}

fn refinements(run: &Computation) -> String {
    let mut output = String::new();
    for (index, step) in run.report.refinements.iter().enumerate() {
        if index != 0 {
            output.push(',');
        }
        let mut marked = String::new();
        for (index, cell) in step.marking.marked.iter().enumerate() {
            if index != 0 {
                marked.push(',');
            }
            let _ = write!(
                marked,
                "{{\"level\":{},\"index\":{:?}}}",
                cell.level(),
                cell.index()
            );
        }
        let mut cases = String::new();
        for (index, case) in step.estimate.cases.iter().enumerate() {
            if index != 0 {
                cases.push(',');
            }
            let _ = write!(
                cases,
                "{{\"dwr\":{:.17e},\"coarse_space\":{:.17e},\"goal_transfer\":{:.17e},\"algebraic\":{:.17e},\"identity_relative_defect\":{:.17e},\"field_residuals\":{:?}}}",
                case.dwr,
                case.coarse_space,
                case.goal_transfer,
                case.algebraic,
                case.identity_relative_defect,
                case.field_residuals
            );
        }
        let _ = write!(
            output,
            "{{\"destination_stage\":{},\"source_active_cells\":{},\"target_active_cells\":{},\"source_background_cells\":{},\"target_background_cells\":{},\"installed\":{},\"marked\":[{marked}],\"marking_fraction\":{:.17e},\"marking_target_met\":{},\"coarse_compliance_j\":{:.17e},\"enriched_compliance_j\":{:.17e},\"two_grid_correction_j\":{:.17e},\"cases\":[{cases}]}}",
            step.destination_stage,
            step.source_active_cells,
            step.target_active_cells,
            step.source_background_cells,
            step.target_background_cells,
            step.installed,
            step.marking.achieved_fraction,
            step.marking.target_met,
            step.estimate.coarse_value,
            step.estimate.fine_value,
            step.estimate.correction
        );
    }
    format!("[{output}]")
}

fn rejected_gradient(run: &Computation) -> String {
    let Some(check) = &run.report.continuation.rejected_gradient_check else {
        return "null".into();
    };
    let mut probes = String::new();
    for (index, probe) in check.probes.iter().enumerate() {
        if index != 0 {
            probes.push(',');
        }
        let _ = write!(
            probes,
            "{{\"direction\":\"{:?}\",\"compliance_relative_error\":{:.17e},\"volume_relative_error\":{:.17e}}}",
            probe.direction, probe.compliance_relative_error, probe.volume_relative_error
        );
    }
    format!(
        "{{\"penal\":{},\"beta\":{},\"step\":{},\"relative_tolerance\":{},\"probes\":[{probes}]}}",
        check.params.penal, check.params.beta, check.options.step, check.options.relative_tolerance
    )
}

pub(super) fn persist(spec: &Spec, ledger: &Ledger, run: &Computation) -> Result<Outcome> {
    let (rows, table, updates) = stages(run);
    let design = design(run)?;
    let trace = hash_domain("org.frankensim.cli.sdf3-study.trace.v1", rows.as_bytes());
    let report = &run.report.continuation;
    let final_row = report.last.as_ref().and_then(|last| last.history.last());
    let compliance = final_row.map_or("null".into(), |row| format!("{:.17e}", row.compliance));
    let volume = final_row.map_or("null".into(), |row| format!("{:.17e}", row.volume_fraction));
    let work = format!(
        "{{\"linear_iterations\":{},\"linear_solves\":{},\"preconditioner_operator_applications\":{},\"geometry_boxes\":{},\"geometry_points\":{},\"geometry_field_evaluations\":{}}}",
        report.work.linear_iterations,
        report.work.linear_solves,
        report.work.preconditioner_operator_applications,
        run.geometry.boxes,
        run.geometry.points,
        run.geometry.field_evaluations
    );
    let refinement = refinements(run);
    let stop = format!(
        "{{\"termination\":\"{:?}\",\"stopped_stage\":{},\"evaluation\":{},\"refinement\":{},\"rejected_gradient_check\":{}}}",
        report.termination,
        report
            .stopped_stage
            .map_or("null".into(), |stage| stage.to_string()),
        report
            .evaluation_stop
            .as_ref()
            .map_or("null".into(), |reason| quoted(&format!("{reason:?}"))),
        run.report
            .refinement_error
            .as_ref()
            .map_or("null".into(), |reason| quoted(&format!("{reason:?}"))),
        rejected_gradient(run)
    );
    let summary = format!(
        "{{\"driver\":{SDF3_DRIVER:?},\"study_id\":\"{}\",\"status\":{:?},\"iterations_completed\":{updates},\"stage_models_retained\":{},\"final_compliance_j\":{compliance},\"final_volume_fraction\":{volume},\"active_cells\":{},\"background_cells\":{},\"work\":{work},\"stop\":{stop},\"refinements\":{refinement},\"stage_evidence\":{rows},\"authority\":\"Estimated\",\"no_claim\":{}}}",
        spec.id.to_hex(),
        run.status,
        report.stages.len(),
        run.study.cells(),
        run.background_cells,
        quoted(SCOPE)
    );
    let html = format!(
        "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Adaptive 3-D topology study</title><body><h1>Adaptive 3-D topology study</h1><p>Status: {}. Estimated numerical evidence.</p><p>{SCOPE}</p><p>Final compliance: {compliance} J; projected volume fraction: {volume}.</p><p>Retained active cells: {}. Installed background refinements: {}.</p><table><tr><th>Stage</th><th>Update</th><th>Compliance J</th><th>Volume fraction</th></tr>{table}</table><p>The JSON report retains gradient discrepancies, every two-grid correction and stop cause. The design artifact contains raw/projected cell densities and physical nodal displacements for every independent load.</p></body></html>",
        run.status,
        run.study.cells(),
        run.report
            .refinements
            .iter()
            .filter(|step| step.installed)
            .count()
    );
    let constellation = hash_bytes(include_bytes!("../../../../../constellation.lock")).to_hex();
    let package = EvidencePackage::new(Provenance::new(
        format!("fs-cli/{}+{SDF3_DRIVER}", env!("CARGO_PKG_VERSION")),
        &constellation,
    ));
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
        "{{\"driver\":{SDF3_DRIVER:?},\"crate\":{:?},\"constellation_lock\":{constellation:?},\"policy\":\"gradient-step=1e-4; gradient-rtol=5e-4; volume-atol=1e-8; quadrature-depth=2; two-level-default-bounds\"}}",
        env!("CARGO_PKG_VERSION")
    );
    let budgets = format!(
        "{{\"wall_s\":{},\"consumed_wall_s\":{},\"memory_bytes\":{},\"linear_iterations\":{},\"per_solve_iterations\":{},\"quadrature_boxes\":{},\"quadrature_points\":{},\"work\":{work}}}",
        spec.wall_s, run.wall_s, spec.memory, spec.linear, spec.per_solve, spec.boxes, spec.points
    );
    let artifacts = [
        ("iterations", "study-iterations", rows.as_bytes()),
        ("design", "study-design", design.as_bytes()),
        ("report_html", "study-report-html", html.as_bytes()),
        ("report_json", "study-report-json", summary.as_bytes()),
        ("package", "study-package", package.as_bytes()),
    ];
    if artifacts
        .iter()
        .any(|(_, _, bytes)| bytes.len() as u64 > MAX_ARTIFACT_BYTES)
    {
        return Err(fail(
            "cli-study-sdf3-artifact",
            "retained field/report exceeds the 16 MiB artifact envelope",
        ));
    }
    ledger.begin()?;
    let result = (|| -> Result<Outcome> {
        let op = ledger.begin_op(Some(spec.id.as_bytes()), &ir_for(SDF3_DRIVER, spec.id, updates), &FiveExplicits {
            seed: &spec.seed.to_le_bytes(), versions: &versions, budget: &budgets,
            capability: "{\"ops\":[\"optimization.marquee-topopt\",\"geometry.sdf\",\"physics.cutfem\"]}",
        }, 0)?;
        let source = ledger.put_artifact("study-source", spec.canonical.as_bytes(), None)?;
        ledger.link(op, &source.hash, EdgeRole::In)?;
        let mut refs = String::new();
        for (name, kind, bytes) in artifacts {
            let artifact = ledger.put_artifact(kind, bytes, None)?;
            ledger.link(op, &artifact.hash, EdgeRole::Out)?;
            let _ = write!(refs, ",{name:?}:\"{}\"", artifact.hash.to_hex());
        }
        let receipt = format!(
            "{{\"schema\":{STUDY_RUN_RECEIPT_SCHEMA:?},\"driver\":{SDF3_DRIVER:?},\"study_id\":\"{}\",\"status\":{:?},\"source\":\"{}\",\"iterations_completed\":{updates},\"target_iterations\":{},\"trace_hash\":\"{}\",\"consumed_wall_s\":{},\"work\":{work},\"stop\":{stop},\"predecessor\":null,\"resume_supported\":false{refs}}}",
            spec.id.to_hex(),
            run.status,
            source.hash.to_hex(),
            spec.updates * spec.schedule.len(),
            trace.to_hex(),
            run.wall_s
        );
        let artifact = ledger.put_artifact(RECEIPT_KIND, receipt.as_bytes(), None)?;
        ledger.link(op, &artifact.hash, EdgeRole::Out)?;
        ledger.seal_artifact_output(&artifact.hash, op)?;
        ledger.finish_op(op, OpOutcome::Ok, None, 1)?;
        Ok(Outcome {
            pointer: format!("study-{}", artifact.hash.to_hex()),
            receipt,
            status: run.status,
        })
    })();
    match result {
        Ok(result) => match ledger.commit() {
            Ok(()) => Ok(result),
            Err(error) => {
                ledger.rollback()?;
                Err(error.into())
            }
        },
        Err(error) => {
            ledger.rollback()?;
            Err(error)
        }
    }
}
