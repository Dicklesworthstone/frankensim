//! File-driven physical inverse design. The library owns parsing, physics,
//! adjoints, SQP steps and KKT checks; this executable only connects those owners.
#[path = "equilibrium_fit/scenarios.rs"]
mod scenarios;

use fs_ascent::{EquilibriumStudy, SqpRunReport};
use fs_couple::render::schedule::force::file::MAX_MODAL_PERFORMANCE_BYTES;
use fs_couple::render::schedule::force::file::design::{
    EquilibriumDesignFile, MAX_EQUILIBRIUM_DESIGN_BYTES,
};
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::{
    DesignControl, DesignEvaluation, DesignWork,
    constraints::{ConstraintSense, ResponseQuantity},
};
use fs_exec::CancelGate;
use std::fmt::Write as _;
use std::io::Read;

const USAGE: &str = "equilibrium_fit MODEL.performance DESIGN.fit [--scenarios TOLERANCES.txt] [--iterations N] [--evaluations N] [--tolerance T] [--max-kkt-dimension N]";

#[derive(Clone, Copy, Debug, PartialEq)]
struct Limits { iterations: usize, evaluations: usize, tolerance: f64, kkt_dimension: usize }
impl Default for Limits {
    fn default() -> Self { Self { iterations:128, evaluations:256, tolerance:1e-8, kkt_dimension:384 } }
}
struct Options { model: String, design: String, scenarios: Option<String>, limits: Limits }
fn options(args: &[String]) -> Result<Options, String> {
    let mut paths = Vec::new();
    let mut limits = Limits::default();
    let mut scenario_file = None;
    let mut seen = Vec::new();
    let mut args = args.iter();
    while let Some(arg) = args.next() {
        if !arg.starts_with('-') { paths.push(arg.clone()); continue; }
        if seen.contains(arg) { return Err(format!("duplicate option {arg}")); }
        seen.push(arg.clone());
        let value = args.next().ok_or_else(|| format!("missing value for {arg}"))?;
        match arg.as_str() {
            "--scenarios" => scenario_file = Some(value.clone()),
            "--iterations" => limits.iterations = value.parse::<usize>().ok().filter(|v| *v <= 512)
                .ok_or("--iterations must be in 0..=512")?,
            "--evaluations" => limits.evaluations = value.parse::<usize>().ok().filter(|v| (2..=4096).contains(v))
                .ok_or("--evaluations must be in 2..=4096, including final re-solve")?,
            "--tolerance" => limits.tolerance = value.parse::<f64>().ok().filter(|v| v.is_finite() && *v > 0.0)
                .ok_or("--tolerance must be finite and positive")?,
            "--max-kkt-dimension" => limits.kkt_dimension = value.parse::<usize>().ok().filter(|v| (1..=384).contains(v))
                .ok_or("--max-kkt-dimension must be in 1..=384")?,
            _ => return Err(format!("unknown option {arg}")),
        }
    }
    let [model, design] = paths.as_slice() else { return Err(USAGE.into()); };
    Ok(Options { model:model.clone(), design:design.clone(), scenarios:scenario_file, limits })
}
fn read_bounded(path: &str, limit: usize) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)?.take((limit+1) as u64).read_to_end(&mut bytes)?;
    if bytes.len() > limit { return Err(format!("input exceeds {limit} bytes: {path}").into()); }
    Ok(bytes)
}

struct Outcome { initial: f64, report: SqpRunReport, audited: DesignEvaluation, work: DesignWork }
fn solve(loaded: &EquilibriumDesignFile, limits: Limits, gate: &CancelGate)
    -> Result<Outcome, Box<dyn std::error::Error>>
{
    let problem = loaded.problem();
    let case_limit = limits.evaluations.checked_mul(problem.load_cases().len()).ok_or("case budget overflow")?;
    let mut control = DesignControl::new(limits.evaluations, case_limit);
    let (initial, report, cached) = {
        let mut study = EquilibriumStudy::new(problem, &vec![0.0;problem.variables().len()],
            &mut control, limits.kkt_dimension, gate)?;
        let initial = study.accepted().value;
        // Reserve one complete physical/adjoint evaluation within the TOTAL
        // allowance. The SQP engine, not an infinity barrier, owns the bounds.
        let report = study.run(limits.tolerance, limits.iterations, limits.evaluations-1, gate)?;
        (initial, report, study.accepted().clone())
    };
    let audited = problem.evaluate(&report.solution.x, &mut control, gate)?;
    // The immutable problem and deterministic owner must reproduce the complete
    // accepted sample. Do not attach cached multipliers to a changed objective.
    if audited != cached { return Err("final physical re-solve differs from the accepted objective/derivatives".into()); }
    Ok(Outcome { initial, report, audited, work:control.work() })
}

fn json_string(text: &str) -> String {
    let mut out = String::from("\"");
    for ch in text.chars() {
        match ch {
            '"' => out.push_str("\\\""), '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"), '\r' => out.push_str("\\r"), '\t' => out.push_str("\\t"),
            ch if ch <= '\u{001f}' => { write!(&mut out, "\\u{:04x}", u32::from(ch)).expect("String write"); }
            ch => out.push(ch),
        }
    }
    out.push('"'); out
}
fn output(loaded: &EquilibriumDesignFile, result: &Outcome) -> String {
    let s = &result.report.solution;
    let kkt = &s.kkt;
    let mut out = format!("{{\"schema\":\"frankensim-equilibrium-fit-v1\",\"scope\":\"local-static-authored-model\",\"model_blake3\":\"{}\",\"design_blake3\":\"{}\",\"stop\":\"{:?}\",\"converged\":{},\"initial_objective\":{:.17e},\"objective\":{:.17e},\"iterations\":{},\"evaluations_including_audit\":{},\"case_solves\":{},\"kkt\":{{\"stationarity\":{:.17e},\"feasibility\":{:.17e},\"dual_feasibility\":{:.17e},\"complementarity\":{:.17e}}},\"parameters\":[",
        loaded.model_info().input_hash.to_hex(), loaded.design_hash().to_hex(), result.report.stop, s.converged,
        result.initial, result.audited.value, s.iters, result.work.evaluations, result.work.case_solves,
        kkt.stationarity, kkt.feasibility, kkt.dual_feasibility, kkt.complementarity);
    for (i, (variable, value)) in loaded.problem().variables().iter().zip(&result.audited.physical_parameters).enumerate() {
        if i != 0 { out.push(','); }
        write!(&mut out, "{{\"name\":{},\"value\":{value:.17e},\"decision\":{:.17e},\"lower_multiplier_decision\":{:.17e},\"upper_multiplier_decision\":{:.17e}}}",
            json_string(&variable.name), s.x[i], s.nu[2*i], s.nu[2*i+1]).expect("String write");
    }
    out.push(']');
    append_physics(&mut out, loaded, &result.audited, &s.lambda, &s.nu[2*loaded.problem().variables().len()..]);
    out.push('}'); out
}
fn append_physics(out: &mut String, loaded: &EquilibriumDesignFile, evaluation: &DesignEvaluation,
    equality_multipliers: &[f64], inequality_multipliers: &[f64])
{
    out.push_str(",\"cases\":[");
    for (i, (case, actual)) in loaded.problem().load_cases().iter().zip(&evaluation.cases).enumerate() {
        if i != 0 { out.push(','); }
        write!(out, "{{\"name\":{},\"objective\":{:.17e},\"active_contacts\":{},\"adjoint_relative_residual\":{:.17e},\"observations\":[",
            json_string(&case.name), actual.value, actual.equilibrium.active_contacts, actual.adjoint_relative_residual).expect("String write");
        for (j, (target, observed)) in case.targets.iter().zip(&actual.observations_m).enumerate() {
            if j != 0 { out.push(','); }
            write!(out, "{{\"predicted_m\":{observed:.17e},\"target_m\":{:.17e},\"scale_m\":{:.17e},\"weight\":{:.17e}}}",
                target.target_m, target.scale_m, target.weight).expect("String write");
        }
        out.push_str("]}");
    }
    out.push(']');
    if !loaded.problem().constraints().is_empty() {
        out.push_str(",\"constraints\":[");
        let mut equality = 0;
        let mut inequality = 0;
        for (i, (constraint, row)) in loaded.problem().constraints().iter().zip(&evaluation.constraints).enumerate() {
            if i != 0 { out.push(','); }
            let (sense, multiplier, violation) = match constraint.sense {
                ConstraintSense::Equal => { let dual = equality_multipliers[equality]; equality += 1; ("equal", dual, row.residual.abs()) }
                other => {
                    let dual = inequality_multipliers[inequality]; inequality += 1;
                    (if other == ConstraintSense::AtMost { "at-most" } else { "at-least" }, dual, row.residual.max(0.0))
                }
            };
            let (quantity, unit) = match &constraint.quantity {
                ResponseQuantity::Displacement(_) => ("displacement", "m"),
                ResponseQuantity::SpringForce(_) => ("spring-force", "N"),
                ResponseQuantity::ContactForce(_) => ("contact-force", "N"),
                ResponseQuantity::ContactPenetration(_) => ("contact-penetration", "m"),
            };
            write!(out, "{{\"name\":{},\"case\":{},\"quantity\":\"{quantity}\",\"unit\":\"{unit}\",\"sense\":\"{sense}\",\"value\":{:.17e},\"bound\":{:.17e},\"scale\":{:.17e},\"residual\":{:.17e},\"violation\":{violation:.17e},\"multiplier_normalized\":{multiplier:.17e},\"adjoint_relative_residual\":{:.17e}}}",
                json_string(&constraint.name), constraint.case, row.value, constraint.bound, constraint.scale,
                row.residual, row.adjoint_relative_residual).expect("String write");
        }
        out.push(']');
    }
}
fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() == 1 && args[0] == "--help" { println!("{USAGE}"); return Ok(()); }
    let options = options(&args)?;
    let model = read_bounded(&options.model, MAX_MODAL_PERFORMANCE_BYTES)?;
    let design = read_bounded(&options.design, MAX_EQUILIBRIUM_DESIGN_BYTES)?;
    let gate = CancelGate::new();
    let loaded = EquilibriumDesignFile::from_bytes(&model, &design, &gate)?;
    if let Some(path) = options.scenarios {
        let bytes = read_bounded(&path, scenarios::MAX_SCENARIO_BYTES)?;
        println!("{}", scenarios::run(&loaded, &bytes, options.limits, &gate)?);
    } else {
        let result = solve(&loaded, options.limits, &gate)?;
        println!("{}", output(&loaded, &result));
    }
    Ok(())
}
fn main() {
    if let Err(error) = run() { eprintln!("equilibrium_fit refused: {error}"); std::process::exit(1); }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn malformed_or_duplicate_execution_limits_are_refused() {
        for tail in [vec!["--evaluations","1"], vec!["--iterations","513"], vec!["--tolerance","NaN"],
            vec!["--max-kkt-dimension","385"], vec!["--unknown","x"], vec!["--iterations"],
            vec!["--evaluations","4","--evaluations","5"]] {
            let args: Vec<String> = ["model","design"].into_iter().chain(tail).map(str::to_owned).collect();
            assert!(options(&args).is_err());
        }
    }
    #[test]
    fn user_supplied_names_are_escaped_without_losing_their_values() {
        assert_eq!(json_string("a\"b\\c\n\u{0001}"), "\"a\\\"b\\\\c\\n\\u0001\"");
        assert_eq!(json_string("force-α"), "\"force-α\"");
    }
}
