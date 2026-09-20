//! Optional finite-tolerance minimax execution, reusing the existing command's
//! model/design loading, physical reporting and execution-limit interpretation.
use super::{append_physics, json_string, ConstraintSense, DesignControl, EquilibriumDesignFile, Limits};
use fs_ascent::equilibrium::scenarios::{EquilibriumScenario, ScenarioEquilibriumStudy, ScenarioProblem};
use fs_exec::CancelGate;
use std::fmt::Write as _;

pub(super) const MAX_SCENARIO_BYTES: usize = 64 * 1024;
const SCHEMA: &str = "frankensim-equilibrium-scenarios-v1";

fn row<'a>(lines: &mut std::str::Lines<'a>, number: &mut usize) -> Result<Vec<&'a str>, String> {
    *number += 1;
    lines.next().map(|r| r.split_ascii_whitespace().collect())
        .ok_or_else(|| format!("missing scenario record on line {number}"))
}
fn parse(bytes: &[u8], names: &[&str]) -> Result<Vec<EquilibriumScenario>, String> {
    if bytes.len() > MAX_SCENARIO_BYTES || names.is_empty() || names.len() > 128 {
        return Err("scenario input exceeds 64 KiB or its variable family is invalid".into());
    }
    let text = std::str::from_utf8(bytes).map_err(|_| "scenario input must be UTF-8")?;
    let mut lines = text.lines(); let mut number = 0;
    if row(&mut lines, &mut number)? != [SCHEMA] { return Err("unsupported scenario input schema".into()); }
    let fields = row(&mut lines, &mut number)?;
    let ["scenarios", count] = fields.as_slice() else { return Err("expected scenarios COUNT".into()); };
    let count = count.parse::<usize>().ok().filter(|c| (1..=32).contains(c))
        .ok_or("scenario count must be in 1..=32")?;
    let mut result = Vec::with_capacity(count);
    for _ in 0..count {
        let fields = row(&mut lines, &mut number)?;
        let ["scenario", name] = fields.as_slice() else { return Err(format!("expected scenario NAME on line {number}")); };
        if name.len() > 128 { return Err("scenario name exceeds 128 bytes".into()); }
        let mut scenario = EquilibriumScenario { name: (*name).to_owned(), physical_offsets: vec![0.0; names.len()] };
        let mut seen = vec![false; names.len()];
        for _ in 0..names.len() {
            let fields = row(&mut lines, &mut number)?;
            let ["offset", name, value] = fields.as_slice() else {
                return Err(format!("expected offset VARIABLE PHYSICAL_VALUE on line {number}"));
            };
            let index = names.iter().position(|n| n == name).ok_or_else(|| format!("unknown design variable {name}"))?;
            if seen[index] { return Err(format!("duplicate offset for {name}")); }
            let value = value.parse::<f64>().ok().filter(|v| v.is_finite()).ok_or("scenario offsets must be finite")?;
            scenario.physical_offsets[index] = value; seen[index] = true;
        }
        result.push(scenario);
    }
    if lines.next().is_some() { return Err(format!("unexpected trailing scenario record on line {}", number+1)); }
    Ok(result)
}

pub(super) fn run(loaded: &EquilibriumDesignFile, bytes: &[u8], limits: Limits, gate: &CancelGate)
    -> Result<String, Box<dyn std::error::Error>>
{
    let problem = loaded.problem();
    let names: Vec<&str> = problem.variables().iter().map(|v| v.name.as_str()).collect();
    let scenarios = parse(bytes, &names)?;
    let count = scenarios.len(); let n = names.len();
    // Keep --evaluations as a TOTAL physical-evaluation allowance, rather than
    // multiplying the old budget silently by the number of scenarios.
    if limits.evaluations / count < 2 {
        return Err("--evaluations must admit initialization and final re-solve of EVERY scenario".into());
    }
    let callback_limit = limits.evaluations/count-1;
    let case_limit = limits.evaluations.checked_mul(problem.load_cases().len()).ok_or("scenario case budget overflow")?;
    let ensemble = ScenarioProblem::new(problem, scenarios, 32, limits.kkt_dimension, gate)?;
    let mut control = DesignControl::new(limits.evaluations, case_limit);
    let mut study = ScenarioEquilibriumStudy::new(ensemble, &vec![0.0;n], &mut control, gate)?;
    let initial = study.accepted().worst_objective;
    let report = study.run(limits.tolerance, limits.iterations, callback_limit, gate)?;
    let audit = study.recheck(gate)?;
    let work = study.work();
    let s = &report.solution; let k = &s.kkt;
    let epigraph = s.x[n];
    let mut out = format!("{{\"schema\":\"frankensim-equilibrium-scenario-fit-v1\",\"scope\":\"local-static-finite-scenario-minimax\",\"model_blake3\":\"{}\",\"design_blake3\":\"{}\",\"stop\":\"{:?}\",\"converged\":{},\"scenario_count\":{count},\"initial_worst_objective\":{initial:.17e},\"worst_objective\":{:.17e},\"epigraph\":{epigraph:.17e},\"epigraph_violation\":{:.17e},\"iterations\":{},\"ensemble_evaluations_including_audit\":{},\"physical_evaluations_including_audit\":{},\"case_solves\":{},\"kkt\":{{\"stationarity\":{:.17e},\"feasibility\":{:.17e},\"dual_feasibility\":{:.17e},\"complementarity\":{:.17e}}},\"parameters\":[",
        loaded.model_info().input_hash.to_hex(), loaded.design_hash().to_hex(), report.stop, s.converged,
        audit.worst_objective, (audit.worst_objective-epigraph).max(0.0), s.iters, s.evals+1,
        work.evaluations, work.case_solves, k.stationarity, k.feasibility, k.dual_feasibility, k.complementarity);
    let (lower, upper) = study.ensemble().decision_bounds();
    for (i, (v, value)) in problem.variables().iter().zip(&audit.nominal_parameters).enumerate() {
        if i != 0 { out.push(','); }
        write!(&mut out, "{{\"name\":{},\"value\":{value:.17e},\"decision\":{:.17e},\"common_lower_decision\":{:.17e},\"common_upper_decision\":{:.17e},\"lower_multiplier_decision\":{:.17e},\"upper_multiplier_decision\":{:.17e}}}",
            json_string(&v.name), s.x[i], lower[i], upper[i], s.nu[2*i], s.nu[2*i+1]).expect("String write");
    }
    out.push_str("],\"scenarios\":[");
    let equalities = problem.constraints().iter().filter(|c| c.sense == ConstraintSense::Equal).count();
    let inequalities = problem.constraints().len()-equalities;
    for (i, (scenario, actual)) in study.ensemble().scenarios().iter().zip(&audit.scenarios).enumerate() {
        if i != 0 { out.push(','); }
        let risk_index = 2*n+i*(1+inequalities);
        write!(&mut out, "{{\"name\":{},\"objective\":{:.17e},\"epigraph_residual\":{:.17e},\"epigraph_multiplier\":{:.17e},\"realized_parameters\":[",
            json_string(&scenario.name), actual.value, actual.value-epigraph, s.nu[risk_index]).expect("String write");
        for (j, (v, value)) in problem.variables().iter().zip(&actual.physical_parameters).enumerate() {
            if j != 0 { out.push(','); }
            write!(&mut out, "{{\"name\":{},\"physical_offset\":{:.17e},\"value\":{value:.17e}}}",
                json_string(&v.name), scenario.physical_offsets[j]).expect("String write");
        }
        out.push(']');
        append_physics(&mut out, loaded, actual, &s.lambda[i*equalities..(i+1)*equalities],
            &s.nu[risk_index+1..risk_index+1+inequalities]);
        out.push('}');
    }
    out.push_str("]}");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn named_offsets_bind_by_identity_not_column_position() {
        let input=b"frankensim-equilibrium-scenarios-v1\nscenarios 2\nscenario low\noffset stiffness -20\noffset gap 0.00001\nscenario high\noffset gap -0.00002\noffset stiffness 30\n";
        let values=parse(input,&["gap","stiffness"]).unwrap();
        assert_eq!(values[0].physical_offsets,[0.00001,-20.0]);
        assert_eq!(values[1].physical_offsets,[-0.00002,30.0]);
    }
    #[test]
    fn malformed_or_incomplete_scenarios_are_never_silently_dropped() {
        let valid="frankensim-equilibrium-scenarios-v1\nscenarios 1\nscenario test\noffset load-N 0.1\n";
        for input in [valid.replace("scenarios 1","scenarios 0"),valid.replace("scenarios 1","scenarios 33"),
            valid.replace("load-N","missing"),valid.replace("0.1","NaN"),valid.replace("0.1","0.1 extra"),
            valid.replace("offset load-N 0.1\n",""),format!("{valid}ignored\n")] {
            assert!(parse(input.as_bytes(),&["load-N"]).is_err());
        }
        assert!(parse(&vec![b' ';MAX_SCENARIO_BYTES+1],&["load-N"]).is_err());
        assert!(parse(b"\xff",&["load-N"]).is_err());
        let duplicate="frankensim-equilibrium-scenarios-v1\nscenarios 1\nscenario test\noffset a 0\noffset a 1\n";
        assert!(parse(duplicate.as_bytes(),&["a","b"]).is_err());
    }
}
