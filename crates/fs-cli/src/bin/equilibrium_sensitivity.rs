//! Authored model/experiment files -> observation Jacobian and local information.
//! The physical/adjoint and dense algebra owners remain in their existing crates.
use fs_couple::render::schedule::force::file::{MAX_MODAL_PERFORMANCE_BYTES, design::{EquilibriumDesignFile, MAX_EQUILIBRIUM_DESIGN_BYTES}};
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::{
    DesignControl, observations::{ObservationControl, ObservationEvaluation, information::{ObservationInformation, admit_information}},
};
use fs_exec::CancelGate;
use std::{collections::BTreeMap, fmt::Write as _, io::Read};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const USAGE: &str = "equilibrium_sensitivity MODEL.performance DESIGN.fit --rank-relative-tolerance R --max-observations N --max-adjoints N --point-x NAME VALUE...";
const SCOPE: &str = "Local numerical information of weighted displacement observations in declared dimensionless design coordinates; no global/structural identifiability, calibrated noise, covariance, posterior, physical validation or active-constraint identifiability certificate. Gram analysis squares conditioning; ambiguous numerical ranks remain null";
struct Options {
    model: String, design: String, relative: f64, maximum_rows: usize, maximum_adjoints: usize,
    point: BTreeMap<String, f64>,
}
fn finite(s: &str) -> Result<f64> {
    Ok(s.parse::<f64>().ok().filter(|x| x.is_finite()).ok_or("expected a finite number")?)
}
fn options(args: &[String]) -> Result<Options> {
    if args.len() < 2 || args.len() > 1024 { return Err(USAGE.into()); }
    let mut point = BTreeMap::new(); let mut relative = None; let mut rows = None; let mut adjoints = None;
    let mut i = 2;
    while i < args.len() {
        let flag = args[i].as_str(); i += 1;
        let value = args.get(i).ok_or("option requires a value")?; i += 1;
        match flag {
            "--point-x" => {
                let x = finite(args.get(i).ok_or("--point-x requires NAME VALUE")?)?; i += 1;
                if value.is_empty() || point.insert(value.clone(), x).is_some() {
                    return Err("each variable must have exactly one named point coordinate".into());
                }
            }
            "--rank-relative-tolerance" if relative.is_none() => relative = Some(finite(value)?),
            "--max-observations" if rows.is_none() => rows = Some(value.parse::<usize>()?),
            "--max-adjoints" if adjoints.is_none() => adjoints = Some(value.parse::<usize>()?),
            _ => return Err(format!("unknown or duplicate option {flag}").into()),
        }
    }
    let relative = relative.ok_or("--rank-relative-tolerance is required")?;
    // Reject unsupported accuracy and missing budgets before reading files.
    admit_information(1, relative)?;
    let maximum_rows = rows.ok_or("--max-observations is required")?;
    let maximum_adjoints = adjoints.ok_or("--max-adjoints is required")?;
    ObservationControl::new(maximum_rows, maximum_adjoints)?;
    Ok(Options { model: args[0].clone(), design: args[1].clone(), relative, maximum_rows, maximum_adjoints, point })
}
fn read(path: &str, cap: usize) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)?.take((cap+1) as u64).read_to_end(&mut bytes)?;
    if bytes.len() > cap { return Err(format!("input exceeds {cap} bytes: {path}").into()); }
    Ok(bytes)
}
fn quoted(text: &str) -> String {
    let mut result = String::from("\"");
    for ch in text.chars() {
        match ch {
            '"' => result.push_str("\\\""), '\\' => result.push_str("\\\\"),
            ch if ch <= '\u{001f}' => { write!(&mut result, "\\u{:04x}", u32::from(ch)).expect("String write"); }
            ch => result.push(ch),
        }
    }
    result.push('"'); result
}
fn array(values: &[f64]) -> String {
    let mut result = String::from("[");
    for (i, v) in values.iter().enumerate() {
        if i != 0 { result.push(','); }
        write!(&mut result, "{v:.17e}").expect("String write");
    }
    result.push(']'); result
}
fn optional(value: Option<f64>) -> String { value.map_or_else(|| "null".into(), |v| format!("{v:.17e}")) }
fn output(loaded: &EquilibriumDesignFile, data: &ObservationEvaluation, information: &ObservationInformation,
    work: &DesignControl, observations: &ObservationControl) -> String
{
    let p = loaded.problem();
    let mut out = format!("{{\"schema\":\"frankensim-equilibrium-sensitivity-v1\",\"status\":\"complete\",\"evidence\":\"Estimated\",\"scope\":{},\"model_blake3\":\"{}\",\"design_blake3\":\"{}\",\"evaluations\":{},\"case_solves\":{},\"observation_adjoints\":{},\"unassessed_response_constraints\":{},\"objective\":{:.17e},\"objective_gradient\":{},\"residual_definition\":\"sqrt(weight)*(displacement-target_m)/scale_m\",\"variables\":[",
        quoted(SCOPE), loaded.model_info().input_hash.to_hex(), loaded.design_hash().to_hex(),
        work.work().evaluations, work.work().case_solves, observations.adjoints_attempted(),
        data.unassessed_response_constraints(), data.objective(), array(data.gradient()));
    for (i, variable) in p.variables().iter().enumerate() {
        if i != 0 { out.push(','); }
        write!(&mut out, "{{\"name\":{},\"x\":{:.17e},\"physical_value\":{:.17e},\"physical_reference\":{:.17e},\"physical_scale\":{:.17e}}}",
            quoted(&variable.name), data.point()[i], data.physical_parameters()[i], variable.reference, variable.scale).expect("String write");
    }
    out.push_str("],\"observations\":[");
    for (i, row) in data.rows().iter().enumerate() {
        if i != 0 { out.push(','); }
        let case = &p.load_cases()[row.case]; let target = &case.targets[row.target];
        write!(&mut out, "{{\"case\":{},\"case_index\":{},\"target\":{},\"value_m\":{:.17e},\"target_m\":{:.17e},\"scale_m\":{:.17e},\"weight\":{:.17e},\"d_displacement_d_x\":{},\"residual\":{:.17e},\"d_residual_d_x\":{},\"adjoint_relative_residual\":{:.17e}}}",
            quoted(&case.name), row.case, row.target, row.value_m, target.target_m, target.scale_m, target.weight,
            array(&row.gradient_m), row.residual, array(&row.residual_gradient), row.adjoint_relative_residual).expect("String write");
    }
    out.push_str("],\"equilibria\":[");
    for (i, report) in data.equilibria().iter().enumerate() {
        if i != 0 { out.push(','); }
        write!(&mut out, "{{\"case_index\":{},\"primal_residual_to_tolerance\":{:.17e},\"minimum_contact_margin_m\":{},\"active_contacts\":{}}}",
            i, report.primal_residual_to_tolerance, optional(report.minimum_contact_margin_m), report.active_contacts).expect("String write");
    }
    write!(&mut out, "],\"information\":{{\"columns\":{},\"numerical_rank\":{},\"relative_threshold\":{:.17e},\"jacobian_scale\":{:.17e},\"relative_singular_values\":{},\"condition_number\":{},\"eigen_residual_relative\":{:.17e},\"relative_gram_guard\":{:.17e},\"parameter_directions\":[",
        data.point().len(), information.numerical_rank.map_or_else(|| "null".into(), |n| n.to_string()),
        information.relative_threshold, information.jacobian_scale, array(&information.relative_singular_values),
        optional(information.condition_number), information.eigen_residual_relative, information.relative_gram_guard).expect("String write");
    for (i, direction) in information.directions.iter().enumerate() {
        if i != 0 { out.push(','); } out.push_str(&array(direction));
    }
    out.push_str("]}}"); out
}
fn run() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() == 1 && args[0] == "--help" { println!("{USAGE}"); return Ok(()); }
    let options = options(&args)?;
    let model = read(&options.model, MAX_MODAL_PERFORMANCE_BYTES)?;
    let design = read(&options.design, MAX_EQUILIBRIUM_DESIGN_BYTES)?;
    let gate = CancelGate::new();
    let loaded = EquilibriumDesignFile::from_bytes(&model, &design, &gate)?;
    let p = loaded.problem();
    admit_information(p.variables().len(), options.relative)?;
    if options.point.len() != p.variables().len() { return Err("declare --point-x for EVERY design variable".into()); }
    let point: Vec<f64> = p.variables().iter().map(|v| options.point.get(&v.name).copied()
        .ok_or_else(|| format!("missing point coordinate {}", v.name))).collect::<std::result::Result<_,_>>()?;
    let mut work = DesignControl::new(1, p.load_cases().len());
    let mut observations = ObservationControl::new(options.maximum_rows, options.maximum_adjoints)?;
    let data = p.evaluate_observations(&point, &mut work, &mut observations, &gate)
        .map_err(|e| format!("{e}; attempted {} cases and {} observation adjoints", work.work().case_solves, observations.adjoints_attempted()))?;
    let information = data.information(options.relative, &gate)?;
    println!("{}", output(&loaded, &data, &information, &work, &observations));
    Ok(())
}
fn main() { if let Err(error) = run() { eprintln!("equilibrium_sensitivity refused: {error}"); std::process::exit(1); } }
