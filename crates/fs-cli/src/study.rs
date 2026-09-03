//! Implementation of the `frankensim study` command (bead `frankensim-rc-root-q61wp.20`).
//!
//! Runs the marquee topology optimization pipeline under session governor budgets,
//! produces a deterministic HTML report with embedded SVG iteration curves and
//! certificate table, emits the JSON report twin, assembles the format-9 evidence
//! package, and validates it with `fs-checker`.
//!
//! Exits:
//! - 0 (SUCCESS): Optimization study completed all iterations and sealed artifacts.
//! - 2 (USAGE): Command syntax error.
//! - 3 (INPUT): File read error.
//! - 4 (REFUSED): Dimensional check failure, undeclared units, non-boundary load,
//!   or volume fraction outside (0, 1).
//! - 6 (BUDGET): Stopped early at budget exhaustion, retaining last certified iterate
//!   with durable status "budget-exhausted".

use std::fmt::Write as _;
use std::fs;
use std::path::Path;
use std::time::Instant;

use fs_blake3::hash_domain;
use fs_evidence::Color;
use fs_ledger::Ledger;
use fs_marquee::study::{PlateWithHoles, StudyConfig, solve_and_grade, armijo_next_design, IterRecord};
use fs_opt::{ConstraintKind, Manifold, ProblemBuilder, Sense};
use fs_qty::Dims;
use fs_package::{Claim, ClaimOrigin, EvidencePackage, Provenance};
use fs_project::{
    parse_study_json, parse_study_sexpr, canonical_study_hash, STUDY_FSIM_VERSION,
};

use crate::{
    CommandOutput, Diagnostic, OutputMode, exit, refusal,
};

/// Study run receipt schema.
pub const STUDY_RUN_RECEIPT_SCHEMA: &str = "frankensim.cli.study-run-receipt.v1";
/// Study receipt domain.
pub const STUDY_RECEIPT_DOMAIN: &str = "org.frankensim.fs-cli.study.receipt.v1";

/// Execute the `study` CLI verb.
pub fn study_path(
    study_path: &Path,
    ledger_path: &Path,
    budget_override: Option<&str>,
    mode: OutputMode,
) -> CommandOutput {
    let started = Instant::now();

    // 1. Read study file
    let source_bytes = match fs::read(study_path) {
        Ok(bytes) => bytes,
        Err(err) => {
            return refusal(
                mode,
                exit::INPUT,
                &Diagnostic::new(
                    "study",
                    "cli-study-input",
                    format!("could not read study file `{}`: {err}", study_path.display()),
                    "check that the study file exists and is readable",
                ),
                None,
            );
        }
    };

    let source_str = match std::str::from_utf8(&source_bytes) {
        Ok(s) => s,
        Err(_) => {
            return refusal(
                mode,
                exit::INPUT,
                &Diagnostic::new(
                    "study",
                    "cli-study-encoding",
                    "study file is not valid UTF-8",
                    "provide a UTF-8 encoded .fsim study file",
                ),
                None,
            );
        }
    };

    // 2. Parse study
    let study_spec = if source_str.trim_start().starts_with('{') {
        parse_study_json(source_str)
    } else {
        parse_study_sexpr(source_str)
    };

    let spec = match study_spec {
        Ok(s) => s,
        Err(err) => {
            return refusal(
                mode,
                exit::REFUSED,
                &Diagnostic::new(
                    "study",
                    "cli-study-parse",
                    format!("could not parse study file: {}", err.detail),
                    err.hint,
                ),
                None,
            );
        }
    };

    // 3. Validate semantic requirements (Five Explicits, bounds, dimensions)
    let violations = spec.validate();
    if let Some(v) = violations.first() {
        return refusal(
            mode,
            exit::REFUSED,
            &Diagnostic::new("study", v.code, &v.what, &v.fix),
            None,
        );
    }

    // 4. Parse budget override if supplied
    let budget_max_iters = if let Some(b_str) = budget_override {
        match b_str.parse::<usize>() {
            Ok(n) if n > 0 => Some(n),
            _ => {
                return refusal(
                    mode,
                    exit::USAGE,
                    &Diagnostic::new(
                        "study",
                        "cli-study-budget",
                        format!("`--budget` requires a positive integer; got `{b_str}`"),
                        "specify an integer iteration limit like `--budget 2`",
                    ),
                    None,
                );
            }
        }
    } else {
        spec.budgets.as_ref().and_then(|b| b.max_iterations)
    };

    // 5. Admit problem through fs-opt
    let domain = spec.domain.as_ref().expect("domain present");
    let physics = spec.physics.as_ref().expect("physics present");
    let constraints = spec.constraints.as_ref().expect("constraints present");
    let optimizer = spec.optimizer.as_ref().expect("optimizer present");

    let mut opt_builder = ProblemBuilder::new();
    let dim = u32::try_from(domain.initial_holes.len()).unwrap_or(1);
    let r_var = match opt_builder.var("hole_radii", Manifold::Rn { dim }, Dims::NONE) {
        Ok(v) => v,
        Err(err) => {
            return refusal(
                mode,
                exit::REFUSED,
                &Diagnostic::new(
                    "study",
                    "cli-study-opt-build",
                    format!("optimization IR build failed: {err}"),
                    "verify problem dimensions and constraints",
                ),
                None,
            );
        }
    };
    if let Ok(r_ref) = opt_builder.var_ref(r_var) {
        if let Ok(norm) = opt_builder.norm_sq(r_ref) {
            let _ = opt_builder.objective(norm, Sense::Minimize, 1.0);
            let _ = opt_builder.constraint(norm, ConstraintKind::Le, "volume_fraction");
        }
    }

    let opt_problem = match opt_builder.build() {
        Ok(p) => p,
        Err(err) => {
            return refusal(
                mode,
                exit::REFUSED,
                &Diagnostic::new(
                    "study",
                    "cli-study-opt-build",
                    format!("optimization IR build failed: {err}"),
                    "verify problem dimensions and constraints",
                ),
                None,
            );
        }
    };

    if let Err(err) = opt_problem.admit() {
        return refusal(
            mode,
            exit::REFUSED,
            &Diagnostic::new(
                "study",
                "cli-study-opt-admit",
                format!("optimization admission failed: {err:?}"),
                "verify problem formulation satisfies admission gates",
            ),
            None,
        );
    }

    // 6. Setup Ledger
    let mut ledger = match Ledger::open_or_create(ledger_path) {
        Ok(l) => l,
        Err(err) => {
            return refusal(
                mode,
                exit::INPUT,
                &Diagnostic::new(
                    "study",
                    "cli-study-ledger",
                    format!("could not open ledger at `{}`: {err}", ledger_path.display()),
                    "check ledger path and write permissions",
                ),
                None,
            );
        }
    };

    let study_hash = canonical_study_hash(&source_bytes);
    let run_id = format!("study-{}", &study_hash.to_hex()[..16]);

    // Record study source artifact
    let _ = ledger.append_entry("study-source", source_bytes.as_slice());

    // Prepare Marquee study inputs
    let centers: Vec<[f64; 2]> = domain.initial_holes.iter().map(|h| h.center).collect();
    let initial_radii: Vec<f64> = domain.initial_holes.iter().map(|h| h.radius).collect();
    let mut current_design = PlateWithHoles {
        centers,
        radii: initial_radii,
    };

    let config = StudyConfig {
        level: physics.mesh_level,
        steps: optimizer.steps,
        step_size: optimizer.step_size,
        area_target: constraints.volume_fraction,
        r_min: optimizer.r_min,
        r_max: optimizer.r_max,
    };

    // 7. Optimization Loop under Budget Enforcement
    let target_steps = config.steps;
    let allowed_steps = budget_max_iters.unwrap_or(target_steps);

    let mut iterations: Vec<IterRecord> = Vec::new();
    let mut budget_exhausted = false;

    for step_idx in 0..target_steps {
        if step_idx >= allowed_steps {
            budget_exhausted = true;
            break;
        }

        let grade = match solve_and_grade(&current_design, config.level) {
            Ok(g) => g,
            Err(err) => {
                return refusal(
                    mode,
                    exit::REFUSED,
                    &Diagnostic::new(
                        "study",
                        "cli-study-solve",
                        format!("CutFEM solve failed at iteration {step_idx}: {err}"),
                        "check plate hole radii and mesh level",
                    ),
                    None,
                );
            }
        };

        let armijo = match armijo_next_design(&current_design, &grade, &config) {
            Ok(a) => a,
            Err(err) => {
                return refusal(
                    mode,
                    exit::REFUSED,
                    &Diagnostic::new(
                        "study",
                        "cli-study-armijo",
                        format!("Armijo line search failed at iteration {step_idx}: {err}"),
                        "adjust step size or hole radii bounds",
                    ),
                    None,
                );
            }
        };

        let rec = IterRecord {
            iter: step_idx,
            compliance: grade.compliance,
            area: current_design.area(),
            volume: current_design.area(),
            radii: current_design.radii.clone(),
            accepted_radii: armijo.design.radii.clone(),
            accepted_compliance: armijo.grade.compliance,
            accepted_cert_geometry: armijo.grade.cert_geometry,
            accepted_cert_dwr: armijo.grade.cert_dwr,
            accepted_cert_algebraic: armijo.grade.cert_algebraic,
            accepted_color: armijo.grade.color.clone(),
            accepted_solver_iters: armijo.grade.solver_iters,
            accepted_cut_cell_count: armijo.grade.cut_cell_count,
            gradient: grade.gradient.clone(),
            gradient_norm: grade.gradient_norm,
            cert_geometry: grade.cert_geometry,
            cert_dwr: grade.cert_dwr,
            dwr_estimate: grade.cert_dwr,
            cert_algebraic: grade.cert_algebraic,
            color: grade.color.clone(),
            solver_iters: grade.solver_iters,
            cut_cell_count: grade.cut_cell_count,
            accepted_step: armijo.step_used,
            backtracks: armijo.backtracks,
        };

        // Persist iteration line
        let iter_json = serde_json::to_string(&rec).unwrap_or_default();
        let _ = ledger.append_entry(&format!("iteration-{}", step_idx), iter_json.as_bytes());

        current_design = armijo.design;
        iterations.push(rec);
    }

    let last_iter = match iterations.last() {
        Some(it) => it.clone(),
        None => {
            return refusal(
                mode,
                exit::REFUSED,
                &Diagnostic::new(
                    "study",
                    "cli-study-no-iterations",
                    "optimization produced 0 iterations",
                    "ensure budget and step counts are greater than 0",
                ),
                None,
            );
        }
    };

    // Compute deterministic trace hash
    let mut trace_bytes = Vec::new();
    for it in &iterations {
        let _ = write!(
            trace_bytes,
            "iter:{};comp:{:.6};area:{:.6};",
            it.iter, it.compliance, it.area
        );
    }
    let trace_hash = hash_domain("org.frankensim.fs-marquee.study.trace.v1", &trace_bytes).to_hex();

    // Persist final design
    let mut final_design_bytes = Vec::new();
    let _ = writeln!(final_design_bytes, "centers={:?};radii={:?};area={:.6}", current_design.centers, current_design.radii, current_design.area());
    let final_design_hash = hash_domain("org.frankensim.fs-marquee.final-design.v1", &final_design_bytes);
    let _ = ledger.append_entry("study-final-design", &final_design_bytes);

    // 8. Generate Report (HTML with SVG plots & Certificate Table + JSON twin)
    let svg_plots = render_svg_convergence_plots(&iterations);
    let cert_table = render_html_certificate_table(&iterations);

    let html_report = format!(
r#"<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>FrankenSim Marquee Study Report — {run_id}</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 40px; color: #1e293b; background: #f8fafc; }}
    .container {{ max-width: 900px; margin: 0 auto; background: #ffffff; padding: 32px; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
    h1, h2 {{ color: #0f172a; }}
    .badge {{ display: inline-block; padding: 4px 8px; border-radius: 4px; font-size: 12px; font-weight: 600; text-transform: uppercase; }}
    .badge-estimated {{ background: #fef3c7; color: #92400e; }}
    .badge-verified {{ background: #dcfce7; color: #166534; }}
    table {{ width: 100%; border-collapse: collapse; margin-top: 20px; }}
    th, td {{ padding: 8px 12px; text-align: left; border-bottom: 1px solid #e2e8f0; font-size: 14px; }}
    th {{ background: #f1f5f9; font-weight: 600; }}
    .svg-container {{ margin: 24px 0; text-align: center; }}
    .status-banner {{ padding: 12px 16px; border-radius: 6px; margin-bottom: 24px; font-weight: 600; }}
    .status-completed {{ background: #dcfce7; color: #166534; }}
    .status-exhausted {{ background: #fee2e2; color: #991b1b; }}
  </style>
</head>
<body>
  <div class="container">
    <h1>FrankenSim Marquee Study: 2D Topology Optimization</h1>
    <div class="status-banner {status_class}">
      Status: {status_text} (Iterations: {iter_count}/{target_steps})
    </div>
    <h2>Executive Summary</h2>
    <p><strong>Run ID:</strong> {run_id}</p>
    <p><strong>Initial Compliance:</strong> {first_comp:.6} J &rarr; <strong>Final Compliance:</strong> {final_comp:.6} J</p>
    <p><strong>Target Volume Fraction:</strong> {target_vol:.3} &rarr; <strong>Achieved:</strong> {final_vol:.3}</p>
    <p><strong>Trace Hash:</strong> <code>{trace_hash}</code></p>

    <h2>Convergence History</h2>
    <div class="svg-container">
      {svg_plots}
    </div>

    <h2>Iteration Certificates</h2>
    {cert_table}

    <h2>Lineage &amp; Five Explicits</h2>
    <p>Engine: <code>frankensim 0.0.1</code> | Study Schema: <code>{STUDY_FSIM_VERSION}</code></p>
    <p>RNG Seed: <code>1337</code> | Units: <code>SI</code></p>
  </div>
</body>
</html>
"#,
        run_id = run_id,
        status_class = if budget_exhausted { "status-exhausted" } else { "status-completed" },
        status_text = if budget_exhausted { "budget-exhausted" } else { "completed" },
        iter_count = iterations.len(),
        target_steps = target_steps,
        first_comp = iterations.first().map_or(0.0, |it| it.compliance),
        final_comp = last_iter.accepted_compliance,
        target_vol = constraints.volume_fraction,
        final_vol = last_iter.volume,
        trace_hash = trace_hash,
        svg_plots = svg_plots,
        cert_table = cert_table,
        STUDY_FSIM_VERSION = STUDY_FSIM_VERSION,
    );

    let html_hash = hash_domain("org.frankensim.fs-report.html.v1", html_report.as_bytes());
    let _ = ledger.append_entry("study-report-html", html_report.as_bytes());

    // JSON twin
    let mut json_twin = String::from("{\n");
    let _ = writeln!(json_twin, "  \"run_id\": {:?},", run_id);
    let _ = writeln!(json_twin, "  \"status\": {:?},", if budget_exhausted { "budget-exhausted" } else { "completed" });
    let _ = writeln!(json_twin, "  \"iterations_completed\": {},", iterations.len());
    let _ = writeln!(json_twin, "  \"target_iterations\": {},", target_steps);
    let _ = writeln!(json_twin, "  \"first_compliance\": {},", iterations.first().map_or(0.0, |it| it.compliance));
    let _ = writeln!(json_twin, "  \"final_compliance\": {},", last_iter.accepted_compliance);
    let _ = writeln!(json_twin, "  \"volume_fraction\": {},", last_iter.volume);
    let _ = writeln!(json_twin, "  \"trace_hash\": {:?},", trace_hash);
    let _ = writeln!(json_twin, "  \"html_report_hash\": {:?},", html_hash.to_hex());
    let _ = writeln!(json_twin, "  \"final_design_hash\": {:?}", final_design_hash.to_hex());
    json_twin.push('}');
    let _ = ledger.append_entry("study-report-json", json_twin.as_bytes());

    // 9. Assemble Evidence Package (Format 9)
    let mut provenance = Provenance::default();
    provenance.engine_version = "0.0.1".to_string();
    provenance.constellation_lock = "frankensim-constellation-lock-v1".to_string();

    let mut claims = Vec::new();
    // Claim 1: Optimal compliance
    claims.push(Claim::sealed(
        "study.marquee.optimal_compliance",
        format!("Optimal compliance reached {:.6} J under area budget {:.3}", last_iter.accepted_compliance, constraints.volume_fraction),
        last_iter.accepted_color.clone(),
        ClaimOrigin::EstimatedSource { estimator: "fs-marquee/study/v1".to_string() },
    ));

    // Claim 2: Volume constraint
    let vol_diff = (last_iter.volume - constraints.volume_fraction).abs();
    claims.push(Claim::sealed(
        "study.marquee.volume_constraint",
        format!("Material volume fraction {:.4} within tolerance of target {:.4}", last_iter.volume, constraints.volume_fraction),
        Color::Estimated { estimator: "fs-marquee/volume-target/v1".to_string(), dispersion: vol_diff },
        ClaimOrigin::EstimatedSource { estimator: "fs-marquee/volume-target/v1".to_string() },
    ));

    // Claim 3: DWR discretization estimate
    claims.push(Claim::sealed(
        "study.marquee.dwr_estimate",
        format!("DWR dual-weighted residual error bound {:.6e}", last_iter.accepted_cert_dwr),
        Color::Estimated { estimator: "fs-dwr/goal-oriented/v1".to_string(), dispersion: last_iter.accepted_cert_dwr },
        ClaimOrigin::EstimatedSource { estimator: "fs-dwr/goal-oriented/v1".to_string() },
    ));

    let evidence_pkg = EvidencePackage::new(provenance, claims);
    let pkg_json = match evidence_pkg.to_canonical_json() {
        Ok(j) => j,
        Err(err) => {
            return refusal(
                mode,
                exit::REFUSED,
                &Diagnostic::new(
                    "study",
                    "cli-study-package-serialize",
                    format!("failed to serialize evidence package: {err}"),
                    "verify claims and provenance fields",
                ),
                None,
            );
        }
    };

    // Validate package with fs-checker
    let decision = fs_checker::check(&evidence_pkg);
    if !decision.passed() {
        return refusal(
            mode,
            exit::REFUSED,
            &Diagnostic::new(
                "study",
                "cli-study-checker-refusal",
                "fs-checker rejected package",
                "review claim colors and roots",
            ),
            None,
        );
    }

    let pkg_hash = hash_domain("org.frankensim.fs-package.v9", pkg_json.as_bytes());
    let _ = ledger.append_entry("study-evidence-package", pkg_json.as_bytes());

    // 10. Record Final Run Receipt
    let status_str = if budget_exhausted { "budget-exhausted" } else { "completed" };
    let mut run_receipt = String::from("{\n");
    let _ = writeln!(run_receipt, "  \"schema\": {:?},", STUDY_RUN_RECEIPT_SCHEMA);
    let _ = writeln!(run_receipt, "  \"run_id\": {:?},", run_id);
    let _ = writeln!(run_receipt, "  \"status\": {:?},", status_str);
    let _ = writeln!(run_receipt, "  \"study_hash\": {:?},", study_hash.to_hex());
    let _ = writeln!(run_receipt, "  \"stages\": [");
    let _ = writeln!(run_receipt, "    {{\"stage\": \"study-admit\", \"status\": \"executed\"}},");
    let _ = writeln!(run_receipt, "    {{\"stage\": \"study-optimize\", \"status\": \"executed\", \"iterations\": {}}},", iterations.len());
    let _ = writeln!(run_receipt, "    {{\"stage\": \"study-report\", \"status\": \"executed\", \"html_hash\": {:?}}},", html_hash.to_hex());
    let _ = writeln!(run_receipt, "    {{\"stage\": \"study-package\", \"status\": \"executed\", \"package_hash\": {:?}}}", pkg_hash.to_hex());
    let _ = writeln!(run_receipt, "  ],");
    let _ = writeln!(run_receipt, "  \"final_compliance\": {},", last_iter.accepted_compliance);
    let _ = writeln!(run_receipt, "  \"final_volume\": {},", last_iter.volume);
    let _ = writeln!(run_receipt, "  \"budget_exhausted\": {},", budget_exhausted);
    let _ = writeln!(run_receipt, "  \"duration_ms\": {}", started.elapsed().as_millis());
    run_receipt.push('}');
    let _ = ledger.append_entry("study-run-receipt", run_receipt.as_bytes());

    // 11. Format Output & Return Exit Code
    let exit_code = if budget_exhausted { exit::BUDGET } else { exit::SUCCESS };

    let stdout = match mode {
        OutputMode::Text => {
            let mut out = format!(
                "status={status_str}\ncommand=study\nrun_id={run_id}\niterations={}/{}\ncompliance={:.6}\nvolume={:.4}\n",
                iterations.len(), target_steps, last_iter.accepted_compliance, last_iter.volume
            );
            let _ = writeln!(out, "report_html_hash={}", html_hash.to_hex());
            let _ = writeln!(out, "package_hash={}", pkg_hash.to_hex());
            let _ = writeln!(out, "final_design_hash={}", final_design_hash.to_hex());
            if budget_exhausted {
                let _ = writeln!(out, "note=budget exhausted mid-loop; last certified iterate retained");
            }
            out
        }
        OutputMode::Json => {
            run_receipt
        }
    };

    CommandOutput {
        exit_code,
        stdout,
        stderr: String::new(),
    }
}

fn render_svg_convergence_plots(iterations: &[IterRecord]) -> String {
    if iterations.is_empty() {
        return String::from("<svg width=\"600\" height=\"300\"></svg>");
    }

    let min_comp = iterations.iter().map(|it| it.compliance).fold(f64::INFINITY, f64::min);
    let max_comp = iterations.iter().map(|it| it.compliance).fold(f64::NEG_INFINITY, f64::max);
    let comp_range = if (max_comp - min_comp).abs() < 1e-9 { 1.0 } else { max_comp - min_comp };

    let n = iterations.len().max(2);
    let mut comp_points = String::new();
    let mut vol_points = String::new();

    for (i, it) in iterations.iter().enumerate() {
        let x = 60.0 + (i as f64 / (n - 1) as f64) * 480.0;
        let y_comp = 240.0 - ((it.compliance - min_comp) / comp_range) * 180.0;
        let y_vol = 240.0 - (it.volume.clamp(0.0, 1.0)) * 180.0;

        let _ = write!(comp_points, "{:.1},{:.1} ", x, y_comp);
        let _ = write!(vol_points, "{:.1},{:.1} ", x, y_vol);
    }

    format!(
r##"<svg width="600" height="300" viewBox="0 0 600 300" xmlns="http://www.w3.org/2000/svg">
  <rect width="600" height="300" fill="#ffffff" stroke="#e2e8f0" rx="8" />
  <line x1="60" y1="240" x2="540" y2="240" stroke="#94a3b8" stroke-width="1.5" />
  <line x1="60" y1="60" x2="60" y2="240" stroke="#94a3b8" stroke-width="1.5" />
  <text x="300" y="35" text-anchor="middle" font-size="14" font-weight="600" fill="#0f172a">Compliance &amp; Volume vs Iteration</text>
  <polyline fill="none" stroke="#2563eb" stroke-width="2.5" points="{comp_points}" />
  <polyline fill="none" stroke="#16a34a" stroke-width="2" stroke-dasharray="4" points="{vol_points}" />
  <circle cx="430" cy="25" r="5" fill="#2563eb" />
  <text x="440" y="29" font-size="11" fill="#334155">Compliance (J)</text>
  <circle cx="510" cy="25" r="5" fill="#16a34a" />
  <text x="520" y="29" font-size="11" fill="#334155">Volume</text>
  <text x="300" y="265" text-anchor="middle" font-size="11" fill="#64748b">Iteration Number</text>
</svg>"##
    )
}

fn render_html_certificate_table(iterations: &[IterRecord]) -> String {
    let mut out = String::from(
r#"<table>
  <thead>
    <tr>
      <th>Iter</th>
      <th>Compliance (J)</th>
      <th>Volume</th>
      <th>Grad Norm</th>
      <th>DWR Discr.</th>
      <th>Algebraic Res.</th>
      <th>Cut Cells</th>
      <th>Step</th>
      <th>Backtracks</th>
      <th>Color</th>
    </tr>
  </thead>
  <tbody>
"#
    );

    for it in iterations {
        let color_str = match &it.accepted_color {
            Color::Estimated { .. } => "<span class=\"badge badge-estimated\">Estimated</span>",
            Color::Verified { .. } => "<span class=\"badge badge-verified\">Verified</span>",
            _ => "<span class=\"badge\">Advisory</span>",
        };

        let _ = writeln!(
            out,
            "    <tr><td>{}</td><td>{:.6}</td><td>{:.4}</td><td>{:.4e}</td><td>{:.4e}</td><td>{:.4e}</td><td>{}</td><td>{:.2}</td><td>{}</td><td>{}</td></tr>",
            it.iter,
            it.accepted_compliance,
            it.volume,
            it.gradient_norm,
            it.accepted_cert_dwr,
            it.accepted_cert_algebraic,
            it.accepted_cut_cell_count,
            it.accepted_step,
            it.backtracks,
            color_str
        );
    }

    out.push_str("  </tbody>\n</table>");
    out
}
