//! Execute uncertainty plans against caller-authored stationary mechanical models.
//! Parsing, primal equilibrium, contact admission and adjoints remain fs-couple's
//! responsibility; fs-uq owns sampling and statistics. No substitute physics.
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::{
    DesignControl, DesignError, DesignWork, EquilibriumDesign,
};
use fs_couple::render::schedule::force::file::MAX_MODAL_PERFORMANCE_BYTES;
use fs_couple::render::schedule::force::file::design::{
    EquilibriumDesignFile, MAX_EQUILIBRIUM_DESIGN_BYTES,
};
use fs_exec::CancelGate;
use fs_uq::{
    CorrelationModel, ParameterUncertainty, PropagationMethod, QmcConfig,
    QmcExecution, UqExecution, UqPlan, UqStatus,
};
use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::io::Read;

const USAGE: &str = "equilibrium_uq MODEL.performance DESIGN.fit --method mc|rqmc [--replicates R] --samples N --seed U64 --case NAME --target INDEX --limit-m METRES --independent (--uniform-x VARIABLE LOWER UPPER | --fixed-x VARIABLE VALUE)...";
const NO_CLAIM: &str = "Estimated under the declared independent bounded input law and fixed numerical model; no physical-validation, fit-error, finite-grid or inverse-normal bias, discretization, confidence-interval or optional-stopping guarantee";

#[derive(Clone, Debug, PartialEq)]
struct Law { name: String, lo: f64, hi: f64 }
struct Options {
    model: String,
    design: String,
    method: PropagationMethod,
    replicates: Option<usize>,
    samples: usize,
    seed: u64,
    case: String,
    target: usize,
    limit_m: f64,
    laws: Vec<Law>,
}

fn finite(value: &str) -> Result<f64, String> {
    value.parse::<f64>().ok().filter(|v| v.is_finite())
        .ok_or_else(|| format!("expected a finite number, received {value}"))
}
fn next<'a>(args: &mut std::slice::Iter<'a, String>, flag: &str) -> Result<&'a str, String> {
    args.next().map(String::as_str).ok_or_else(|| format!("missing value for {flag}"))
}
fn options(args: &[String]) -> Result<Options, String> {
    if args.len() < 2 || args.len() > 1024 { return Err(USAGE.into()); }
    let mut args = args.iter();
    let model = next(&mut args, "model")?.to_owned();
    let design = next(&mut args, "design")?.to_owned();
    let mut method = None;
    let mut replicates = None;
    let mut samples = None;
    let mut seed = None;
    let mut case = None;
    let mut target = None;
    let mut limit_m = None;
    let mut independent = false;
    let mut laws = Vec::new();
    let mut seen = BTreeSet::new();
    while let Some(flag) = args.next() {
        if flag == "--uniform-x" || flag == "--fixed-x" {
            let name = next(&mut args, flag)?.to_owned();
            let lo = finite(next(&mut args, flag)?)?;
            let hi = if flag == "--uniform-x" { finite(next(&mut args, flag)?)? } else { lo };
            if name.is_empty() || lo > hi || laws.iter().any(|law: &Law| law.name == name) {
                return Err("each variable needs one named, finite, ordered uncertainty declaration".into());
            }
            laws.push(Law { name, lo, hi });
            continue;
        }
        if !seen.insert(flag.as_str()) { return Err(format!("duplicate option {flag}")); }
        if flag == "--independent" { independent = true; continue; }
        let value = next(&mut args, flag)?;
        match flag.as_str() {
            "--method" if value == "mc" => method = Some(PropagationMethod::MonteCarlo),
            "--method" if value == "rqmc" => method = Some(PropagationMethod::QuasiMonteCarlo),
            "--method" => return Err("--method must be mc or rqmc".into()),
            "--replicates" => replicates = Some(value.parse::<usize>().ok()
                .filter(|v| (2..=256).contains(v)).ok_or("--replicates must be in 2..=256")?),
            "--samples" => samples = Some(value.parse::<usize>().ok()
                .filter(|v| (2..=1_000_000).contains(v)).ok_or("--samples must be in 2..=1000000")?),
            "--seed" => seed = Some(value.parse::<u64>().map_err(|_| "--seed must be a decimal u64")?),
            "--case" => case = Some(value.to_owned()),
            "--target" => target = Some(value.parse::<usize>().map_err(|_| "--target must be a zero-based index")?),
            "--limit-m" => limit_m = Some(finite(value)?),
            _ => return Err(format!("unknown option {flag}")),
        }
    }
    if !independent { return Err("declare --independent explicitly; bounds are not a joint probability law".into()); }
    let method = method.ok_or("--method is required")?;
    let samples = samples.ok_or("--samples is required")?;
    match (method, replicates) {
        (PropagationMethod::MonteCarlo, None) => {}
        (PropagationMethod::QuasiMonteCarlo, Some(r))
            if samples % r == 0 && samples / r >= 2 && (samples / r).is_power_of_two() => {}
        _ => return Err("mc forbids --replicates; rqmc needs explicit replicates times a complete power-of-two net of at least two points".into()),
    }
    Ok(Options {
        model, design, method, replicates, samples,
        seed: seed.ok_or("--seed is required")?, case: case.ok_or("--case is required")?,
        target: target.ok_or("--target is required")?, limit_m: limit_m.ok_or("--limit-m is required")?, laws,
    })
}

// Uncertainty is specified in the FILE'S dimensionless design coordinates:
// p = reference + scale*x. Thus units and all sharing come from the admitted
// physical bindings; the sampler neither guesses units nor independently varies
// fields that the design declared to be shared.
fn plan(problem: &EquilibriumDesign, options: &Options) -> Result<(UqPlan, usize), String> {
    let case = problem.load_cases().iter().position(|c| c.name == options.case)
        .ok_or_else(|| format!("unknown load case {}", options.case))?;
    if options.target >= problem.load_cases()[case].targets.len() {
        return Err("target index is outside the selected case's displacement observations".into());
    }
    if options.laws.len() != problem.variables().len() {
        return Err("declare --uniform-x or --fixed-x for EVERY variable; none may be left unstated".into());
    }
    let mut plan = UqPlan::new("selected-equilibrium-displacement-m", options.method, options.samples)
        .with_correlation(CorrelationModel::Independent).with_compliance_threshold(options.limit_m);
    plan.seed = options.seed;
    let mut lower = Vec::new();
    let mut upper = Vec::new();
    for variable in problem.variables() {
        let law = options.laws.iter().find(|law| law.name == variable.name)
            .ok_or_else(|| format!("missing uncertainty declaration for {}", variable.name))?;
        lower.push(law.lo);
        upper.push(law.hi);
        plan.parameters.push(ParameterUncertainty::uniform(&law.name, law.lo, law.hi, "1"));
    }
    // This is domain admission, NOT interval propagation of displacement.
    // The affine decoding is monotone because every admitted scale is positive.
    problem.physical_parameters(&lower).map_err(|e| format!("lower uncertainty endpoint: {e}"))?;
    problem.physical_parameters(&upper).map_err(|e| format!("upper uncertainty endpoint: {e}"))?;
    Ok((plan, case))
}

struct Outcome {
    mean_m: f64,
    mean_standard_error_m: Option<f64>,
    displacement_std_dev_m: Option<f64>,
    compliance_probability: f64,
    compliance_standard_error: Option<f64>,
    completed_replicates: Option<usize>,
    samples: usize,
    work: DesignWork,
}

fn solve(loaded: &EquilibriumDesignFile, options: &Options, gate: &CancelGate)
    -> Result<Outcome, Box<dyn std::error::Error>>
{
    let problem = loaded.problem();
    let (plan, case) = plan(problem, options)?;
    let max_cases = options.samples.checked_mul(problem.load_cases().len()).ok_or("case budget overflow")?;
    let mut work = DesignControl::new(options.samples, max_cases);
    // Each draw rebuilds its own independent equilibria. Original primal,
    // contact-activity and derivative residual refusals remain terminal; no
    // clipping, resampling, partial-case objective or failed-sample deletion.
    match options.method {
        PropagationMethod::MonteCarlo => {
            let mut execution = UqExecution::new(&plan)?;
            let result = execution.advance(options.samples, || gate.is_requested(), |x| {
                let solved = problem.evaluate(x, &mut work, gate)?;
                Ok::<_, DesignError>(solved.cases[case].observations_m[options.target])
            });
            if result.status != UqStatus::Complete {
                return Err(incomplete(result.status, result.samples_evaluated, work.work(),
                    result.rejection_reason.as_deref()).into());
            }
            Ok(Outcome {
                mean_m: result.mean.ok_or("completed sampling omitted its mean")?,
                mean_standard_error_m: result.std_dev.map(|_| result.sampling_error),
                displacement_std_dev_m: result.std_dev,
                compliance_probability: result.probability_of_compliance.ok_or("missing compliance estimate")?,
                compliance_standard_error: None, completed_replicates: None,
                samples: result.samples_evaluated, work: work.work(),
            })
        }
        PropagationMethod::QuasiMonteCarlo => {
            let replicates = options.replicates.ok_or("missing QMC replicate layout")?;
            let config = QmcConfig { replicates, samples_per_replicate: options.samples / replicates };
            let mut execution = QmcExecution::new(&plan, config)?;
            let result = execution.advance(options.samples, || gate.is_requested(), |x| {
                let solved = problem.evaluate(x, &mut work, gate)?;
                Ok::<_, DesignError>(solved.cases[case].observations_m[options.target])
            });
            if result.status != UqStatus::Complete {
                return Err(incomplete(result.status, result.samples_evaluated, work.work(),
                    result.rejection_reason.as_deref()).into());
            }
            let estimate = result.estimate.ok_or("missing complete-net displacement estimate")?;
            let compliance = result.compliance.ok_or("missing complete-net compliance estimate")?;
            // The upstream QMC owner computes these errors across whole nets.
            // Never apply the MC iid confidence sequence or pointwise sd/sqrt(n).
            Ok(Outcome {
                mean_m: estimate.mean, mean_standard_error_m: estimate.standard_error,
                displacement_std_dev_m: None,
                compliance_probability: compliance.mean,
                compliance_standard_error: compliance.standard_error,
                completed_replicates: Some(result.completed_replicates),
                samples: result.samples_evaluated, work: work.work(),
            })
        }
        _ => Err("unsupported uncertainty method".into()),
    }
}

fn incomplete(status: UqStatus, samples: usize, work: DesignWork, reason: Option<&str>) -> String {
    format!("uncertainty execution {status:?} after {samples} calls ({} case solves): {}",
        work.case_solves, reason.unwrap_or("execution interrupted"))
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
fn number(value: Option<f64>) -> String {
    value.map_or_else(|| "null".into(), |v| format!("{v:.17e}"))
}
fn output(loaded: &EquilibriumDesignFile, options: &Options, result: &Outcome) -> String {
    let mut out = format!("{{\"schema\":\"frankensim-equilibrium-uq-v1\",\"method\":\"{}\",\"status\":\"complete\",\"evidence\":\"Estimated\",\"no_claim\":{},\"model_blake3\":\"{}\",\"design_blake3\":\"{}\",\"seed\":\"{}\",\"case\":{},\"target\":{},\"unit\":\"m\",\"compliance_event\":\"displacement_m <= limit_m\",\"limit_m\":{:.17e},\"samples\":{},\"case_solves\":{},\"mean_m\":{:.17e},\"mean_standard_error_m\":{},\"displacement_std_dev_m\":{},\"compliance_probability\":{:.17e},\"uncertainty_coordinates\":\"dimensionless x; physical p = reference + scale*x\",\"dependence\":\"independent\",\"variables\":[",
        if options.method == PropagationMethod::MonteCarlo { "mc" } else { "rqmc" },
        json_string(NO_CLAIM), loaded.model_info().input_hash.to_hex(), loaded.design_hash().to_hex(), options.seed,
        json_string(&options.case), options.target, options.limit_m, result.samples, result.work.case_solves,
        result.mean_m, number(result.mean_standard_error_m), number(result.displacement_std_dev_m), result.compliance_probability);
    for (i, variable) in loaded.problem().variables().iter().enumerate() {
        if i != 0 { out.push(','); }
        // Admission already found exactly one law for every physical variable.
        let law = options.laws.iter().find(|law| law.name == variable.name).expect("admitted law");
        write!(&mut out, "{{\"name\":{},\"law\":\"{}\",\"lower_x\":{:.17e},\"upper_x\":{:.17e},\"physical_reference\":{:.17e},\"physical_scale\":{:.17e}}}",
            json_string(&variable.name), if law.lo == law.hi { "fixed" } else { "uniform" },
            law.lo, law.hi, variable.reference, variable.scale).expect("String write");
    }
    write!(&mut out, "],\"completed_replicates\":{},\"compliance_standard_error\":{},\"standard_error_basis\":\"{}\"}}",
        result.completed_replicates.map_or_else(|| "null".into(), |n| n.to_string()),
        number(result.compliance_standard_error),
        if options.method == PropagationMethod::MonteCarlo { "individual-monte-carlo-observations" }
        else { "complete-independent-scramble-means" }).expect("String write");
    out
}
fn read_bounded(path: &str, cap: usize) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)?.take((cap + 1) as u64).read_to_end(&mut bytes)?;
    if bytes.len() > cap { return Err(format!("input exceeds {cap} bytes: {path}").into()); }
    Ok(bytes)
}
fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() == 1 && args[0] == "--help" { println!("{USAGE}"); return Ok(()); }
    let options = options(&args)?;
    let model = read_bounded(&options.model, MAX_MODAL_PERFORMANCE_BYTES)?;
    let design = read_bounded(&options.design, MAX_EQUILIBRIUM_DESIGN_BYTES)?;
    let gate = CancelGate::new();
    let loaded = EquilibriumDesignFile::from_bytes(&model, &design, &gate)?;
    let outcome = solve(&loaded, &options, &gate)?;
    println!("{}", output(&loaded, &options, &outcome));
    Ok(())
}
fn main() {
    if let Err(error) = run() { eprintln!("equilibrium_uq refused: {error}"); std::process::exit(1); }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args() -> Vec<String> {
        "m d --method mc --samples 64 --seed 73 --case load --target 0 --limit-m 0.1 --independent --uniform-x k -1 1"
            .split_whitespace().map(str::to_owned).collect()
    }
    #[test]
    fn required_laws_controls_and_numeric_inputs_are_not_silently_defaulted() {
        assert!(options(&args()).is_ok());
        for flag in ["--method", "--samples", "--seed", "--case", "--target", "--limit-m"] {
            let mut bad = args();
            let start = bad.iter().position(|a| a == flag).unwrap();
            bad.drain(start..start + 2);
            assert!(options(&bad).is_err(), "{flag}");
        }
        let mut bad = args(); bad.retain(|v| v != "--independent"); assert!(options(&bad).is_err());
        for tail in ["--seed 99", "--uniform-x k 0 1", "--fixed-x k NaN", "--unexpected 1"] {
            let mut bad = args(); bad.extend(tail.split_whitespace().map(str::to_owned));
            assert!(options(&bad).is_err());
        }
    }
    #[test]
    fn randomized_qmc_requires_a_complete_explicit_net_layout() {
        let mut qmc = args();
        qmc[3] = "rqmc".into();
        assert!(options(&qmc).is_err());
        qmc.extend(["--replicates".into(), "4".into()]);
        assert!(options(&qmc).is_ok());
        for count in ["1", "3", "257"] {
            let mut bad = qmc.clone(); *bad.last_mut().unwrap() = count.into();
            assert!(options(&bad).is_err());
        }
        let mut bad = args(); bad.extend(["--replicates".into(), "4".into()]);
        assert!(options(&bad).is_err());
    }
    #[test]
    fn json_names_preserve_quotes_backslashes_and_controls() {
        assert_eq!(json_string("a\"\\\n\u{0001}"), "\"a\\\"\\\\\\n\\u0001\"");
    }
}
