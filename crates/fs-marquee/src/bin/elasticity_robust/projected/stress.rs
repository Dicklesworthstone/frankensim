//! Explicit sampled-stress input and output for the projected CLI.
use super::numbers;
use fs_topols::{RobustSampledStressEvaluation, SampledStressLimit};
use std::error::Error;

pub(super) fn options(args: &[String]) -> Result<(Vec<String>, Option<SampledStressLimit>), Box<dyn Error>> {
    if args.len() > 13 { return Err("too many projected study arguments".into()); }
    let mut positional = Vec::new();
    let mut maximum = None;
    let mut tolerance = None;
    let mut index = 0;
    while index < args.len() {
        let slot = match args[index].as_str() {
            "--stress-limit" => Some(&mut maximum),
            "--stress-tolerance" => Some(&mut tolerance),
            value if value.starts_with("--") => return Err(format!("unknown projected option {value}").into()),
            _ => None,
        };
        if let Some(slot) = slot {
            if slot.is_some() { return Err(format!("duplicate {}", args[index]).into()); }
            let value = args.get(index + 1).ok_or_else(|| format!("{} requires a value", args[index]))?;
            *slot = Some(value.parse::<f64>()?);
            index += 2;
        } else {
            positional.push(args[index].clone());
            index += 1;
        }
    }
    let limit = match maximum {
        Some(maximum) => Some(SampledStressLimit::new(maximum, tolerance.unwrap_or(0.0))?),
        None if tolerance.is_some() => return Err("--stress-tolerance requires --stress-limit".into()),
        None => None,
    };
    Ok((positional, limit))
}

pub(super) fn json(limit: SampledStressLimit, report: Option<&RobustSampledStressEvaluation>) -> String {
    let prefix = format!(
        "\"scope\":\"q1-positive-volume-material-cell-probes-v1\",\"units\":\"normalized_stress\",\"limit\":{:.17e},\"absolute_tolerance\":{:.17e},\"admitted_max\":{:.17e},\"continuous_maximum_certificate\":false,\"physical_allowable_validated\":false",
        limit.max_von_mises, limit.absolute_tolerance, limit.admitted_max(),
    );
    let Some(report) = report else {
        return format!("{{{prefix},\"status\":\"unavailable\"}}");
    };
    let status = if report.worst_sampled_von_mises <= limit.admitted_max() {
        "sampled_feasible"
    } else { "sampled_limit_exceeded" };
    let counts = report.case_sample_counts.iter().map(usize::to_string).collect::<Vec<_>>().join(",");
    let locations = report.case_max_locations.iter().map(|point| numbers(point)).collect::<Vec<_>>().join(",");
    format!(
        "{{{prefix},\"status\":\"{status}\",\"case_maxima\":{},\"case_sample_counts\":[{counts}],\"case_max_locations\":[{locations}],\"worst_von_mises\":{:.17e},\"worst_case\":{},\"snapshot\":\"{:#018x}\"}}",
        numbers(&report.case_sampled_max_von_mises), report.worst_sampled_von_mises,
        report.worst_stress_case, report.snapshot,
    )
}

pub(super) fn fields(
    limit: Option<SampledStressLimit>, report: Option<&RobustSampledStressEvaluation>,
) -> String {
    limit.map_or_else(String::new, |limit| format!(",\"sampled_stress\":{}", json(limit, report)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(values: &[&str]) -> Vec<String> { values.iter().map(|value| (*value).into()).collect() }

    #[test]
    fn stress_options_preserve_positionals_and_validate_declared_limits() {
        let (paths, limit) = options(&args(&["out", "loads.csv", "--stress-limit", "4", "3", "--stress-tolerance", "0.1"])).unwrap();
        assert_eq!(paths, ["out", "loads.csv", "3"]);
        assert_eq!(limit.unwrap(), SampledStressLimit::new(4.0, 0.1).unwrap());
        for invalid in [
            vec!["--stress-limit"], vec!["--stress-limit", "NaN"],
            vec!["--stress-limit", "-1"], vec!["--stress-limit", "inf"],
            vec!["--stress-tolerance", "0.1"],
            vec!["--stress-limit", "4", "--stress-limit", "5"],
            vec!["--stress-limit", "4", "--stress-tolerance", "-1"],
            vec!["--stress-limit", "4", "--unknown", "0"],
        ] {
            assert!(options(&args(&invalid)).is_err(), "{invalid:?}");
        }
        assert!(options(&args(&["out", "loads.csv"])).unwrap().1.is_none());
    }

    #[test]
    fn missing_or_disabled_measurements_never_print_a_feasible_constraint() {
        let limit = SampledStressLimit::new(2.0, 0.0).unwrap();
        assert!(json(limit, None).contains("\"status\":\"unavailable\""));
        assert!(!json(limit, None).contains("sampled_feasible"));
        assert_eq!(fields(None, None), "");
    }
}
