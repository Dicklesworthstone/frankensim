//! Optional finite-tolerance minimax/CVaR execution, reusing the existing command's
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
fn parse(bytes: &[u8], names: &[&str]) -> Result<(Vec<EquilibriumScenario>, Option<f64>), String> {
    if bytes.len() > MAX_SCENARIO_BYTES || names.is_empty() || names.len() > 128 {
        return Err("scenario input exceeds 64 KiB or its variable family is invalid".into());
    }
    let text = std::str::from_utf8(bytes).map_err(|_| "scenario input must be UTF-8")?;
    let mut lines = text.lines(); let mut number = 0;
    if row(&mut lines, &mut number)? != [SCHEMA] { return Err("unsupported scenario input schema".into()); }
    let fields = row(&mut lines, &mut number)?;
    // Optional, explicit risk selection. Old files still mean worst case.
    // Repeated realizations retain empirical mass; no deduplication is valid here.
    let (fields, cvar_alpha) = if fields.first() == Some(&"risk") {
        let ["risk", "empirical-cvar", alpha] = fields.as_slice() else {
            return Err("expected risk empirical-cvar ALPHA".into());
        };
        let alpha = alpha.parse::<f64>().ok().filter(|a| a.is_finite() && *a > 0.0 && *a < 1.0)
            .ok_or("empirical CVaR alpha must be strictly between zero and one")?;
        (row(&mut lines, &mut number)?, Some(alpha))
    } else { (fields, None) };
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
    Ok((result, cvar_alpha))
}

pub(super) fn run(loaded: &EquilibriumDesignFile, bytes: &[u8], limits: Limits, gate: &CancelGate)
    -> Result<String, Box<dyn std::error::Error>>
{
    let problem = loaded.problem();
    let names: Vec<&str> = problem.variables().iter().map(|v| v.name.as_str()).collect();
    let (scenarios, cvar_alpha) = parse(bytes, &names)?;
    let count = scenarios.len(); let n = names.len();
    // Keep --evaluations as a TOTAL physical-evaluation allowance, rather than
    // multiplying the old budget silently by the number of scenarios.
    if limits.evaluations / count < 2 {
        return Err("--evaluations must admit initialization and final re-solve of EVERY scenario".into());
    }
    let callback_limit = limits.evaluations/count-1;
    let case_limit = limits.evaluations.checked_mul(problem.load_cases().len()).ok_or("scenario case budget overflow")?;
    let ensemble = ScenarioProblem::new(problem, scenarios, 32, limits.kkt_dimension, gate)?;
    let ensemble = match cvar_alpha { Some(alpha) => ensemble.with_cvar(alpha)?, None => ensemble };
    let mut control = DesignControl::new(limits.evaluations, case_limit);
    let mut study = ScenarioEquilibriumStudy::new(ensemble, &vec![0.0;n], &mut control, gate)?;
    let initial = study.accepted().worst_objective;
    let report = study.run(limits.tolerance, limits.iterations, callback_limit, gate)?;
    let audit = study.recheck(gate)?;
    let work = study.work();
    let s = &report.solution; let k = &s.kkt;
    let epigraph = s.x[n];
    let mut out = if let Some(alpha) = cvar_alpha {
        let score = study.ensemble().cvar_upper_bound(&audit, epigraph)?;
        let mass = 1.0/(count as f64);
        let tail_mass = (count as f64)*(1.0-alpha);
        let mut tail_violation = 0.0_f64;
        for (i, actual) in audit.scenarios.iter().enumerate() {
            let slack = s.x[n+1+i];
            tail_violation = tail_violation.max(actual.value-epigraph-slack).max(-slack);
        }
        format!("{{\"schema\":\"frankensim-equilibrium-cvar-fit-v1\",\"scope\":\"local-static-equal-mass-empirical-cvar\",\"model_blake3\":\"{}\",\"design_blake3\":\"{}\",\"stop\":\"{:?}\",\"converged\":{},\"scenario_count\":{count},\"alpha\":{alpha:.17e},\"scenario_mass\":{mass:.17e},\"tail_mass\":{tail_mass:.17e},\"initial_cvar_upper_bound\":{initial:.17e},\"objective\":{:.17e},\"cvar_upper_bound\":{score:.17e},\"threshold\":{epigraph:.17e},\"tail_violation\":{tail_violation:.17e},\"worst_objective\":{:.17e},\"iterations\":{},\"ensemble_evaluations_including_audit\":{},\"physical_evaluations_including_audit\":{},\"case_solves\":{},\"kkt\":{{\"stationarity\":{:.17e},\"feasibility\":{:.17e},\"dual_feasibility\":{:.17e},\"complementarity\":{:.17e}}},\"parameters\":[",
            loaded.model_info().input_hash.to_hex(), loaded.design_hash().to_hex(), report.stop, s.converged,
            s.f, audit.worst_objective, s.iters, s.evals+1, work.evaluations, work.case_solves,
            k.stationarity, k.feasibility, k.dual_feasibility, k.complementarity)
    } else {
        format!("{{\"schema\":\"frankensim-equilibrium-scenario-fit-v1\",\"scope\":\"local-static-finite-scenario-minimax\",\"model_blake3\":\"{}\",\"design_blake3\":\"{}\",\"stop\":\"{:?}\",\"converged\":{},\"scenario_count\":{count},\"initial_worst_objective\":{initial:.17e},\"worst_objective\":{:.17e},\"epigraph\":{epigraph:.17e},\"epigraph_violation\":{:.17e},\"iterations\":{},\"ensemble_evaluations_including_audit\":{},\"physical_evaluations_including_audit\":{},\"case_solves\":{},\"kkt\":{{\"stationarity\":{:.17e},\"feasibility\":{:.17e},\"dual_feasibility\":{:.17e},\"complementarity\":{:.17e}}},\"parameters\":[",
        loaded.model_info().input_hash.to_hex(), loaded.design_hash().to_hex(), report.stop, s.converged,
        audit.worst_objective, (audit.worst_objective-epigraph).max(0.0), s.iters, s.evals+1,
        work.evaluations, work.case_solves, k.stationarity, k.feasibility, k.dual_feasibility, k.complementarity)
    };
    let (lower, upper) = study.ensemble().decision_bounds();
    for (i, (v, value)) in problem.variables().iter().zip(&audit.nominal_parameters).enumerate() {
        if i != 0 { out.push(','); }
        write!(&mut out, "{{\"name\":{},\"value\":{value:.17e},\"decision\":{:.17e},\"common_lower_decision\":{:.17e},\"common_upper_decision\":{:.17e},\"lower_multiplier_decision\":{:.17e},\"upper_multiplier_decision\":{:.17e}}}",
            json_string(&v.name), s.x[i], lower[i], upper[i], s.nu[2*i], s.nu[2*i+1]).expect("String write");
    }
    out.push_str("],\"scenarios\":[");
    let equalities = problem.constraints().iter().filter(|c| c.sense == ConstraintSense::Equal).count();
    let inequalities = problem.constraints().len()-equalities;
    let risk_rows = if cvar_alpha.is_some() { 2 } else { 1 };
    for (i, (scenario, actual)) in study.ensemble().scenarios().iter().zip(&audit.scenarios).enumerate() {
        if i != 0 { out.push(','); }
        let risk_index = 2*n+i*(risk_rows+inequalities);
        if cvar_alpha.is_some() {
            let excess = s.x[n+1+i];
            write!(&mut out, "{{\"name\":{},\"objective\":{:.17e},\"tail_excess\":{excess:.17e},\"tail_residual\":{:.17e},\"tail_multiplier\":{:.17e},\"excess_nonnegative_residual\":{:.17e},\"excess_multiplier\":{:.17e},\"realized_parameters\":[",
                json_string(&scenario.name), actual.value, actual.value-epigraph-excess,
                s.nu[risk_index], -excess, s.nu[risk_index+1]).expect("String write");
        } else {
            write!(&mut out, "{{\"name\":{},\"objective\":{:.17e},\"epigraph_residual\":{:.17e},\"epigraph_multiplier\":{:.17e},\"realized_parameters\":[",
                json_string(&scenario.name), actual.value, actual.value-epigraph, s.nu[risk_index]).expect("String write");
        }
        for (j, (v, value)) in problem.variables().iter().zip(&actual.physical_parameters).enumerate() {
            if j != 0 { out.push(','); }
            write!(&mut out, "{{\"name\":{},\"physical_offset\":{:.17e},\"value\":{value:.17e}}}",
                json_string(&v.name), scenario.physical_offsets[j]).expect("String write");
        }
        out.push(']');
        append_physics(&mut out, loaded, actual, &s.lambda[i*equalities..(i+1)*equalities],
            &s.nu[risk_index+risk_rows..risk_index+risk_rows+inequalities]);
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
        let (values, risk)=parse(input,&["gap","stiffness"]).unwrap();
        assert_eq!(risk,None);
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
    #[test]
    fn empirical_cvar_requires_explicit_finite_risk_selection() {
        let input="frankensim-equilibrium-scenarios-v1\nrisk empirical-cvar 0.4\nscenarios 1\nscenario one\noffset load-N 0\n";
        let (entries,risk)=parse(input.as_bytes(),&["load-N"]).unwrap();
        assert_eq!(entries.len(),1);assert_eq!(risk,Some(0.4));
        for value in ["0","1","-1","NaN","inf","0.4 extra"] {
            assert!(parse(input.replace("0.4",value).as_bytes(),&["load-N"]).is_err());
        }
        assert!(parse(input.replace("empirical-cvar","mean").as_bytes(),&["load-N"]).is_err());
        assert!(parse(input.replace("scenarios 1","risk empirical-cvar 0.7\nscenarios 1").as_bytes(),&["load-N"]).is_err());
    }

}
