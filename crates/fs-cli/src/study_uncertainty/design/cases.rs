//! Fixed operating-case projection over the existing design driver. Every
//! physical child uses Model.sample; only a fully observed group can supply a
//! feasible incumbent. Partial groups remain ordinary retained observations.
use std::borrow::Cow;
use fs_project::uncertainty::UncertaintyStudy;
use super::{DesignPolicy, Trial, quoted};

pub(super) fn count(study: &UncertaintyStudy) -> usize {
    study.operating_cases().len().max(1)
}

pub(super) fn point(study: &UncertaintyStudy, coordinate: f64, completed: usize) -> Vec<f64> {
    let mut values = Vec::with_capacity(study.parameters().len());
    values.push(coordinate);
    if !study.operating_cases().is_empty() {
        values.extend_from_slice(&study.operating_cases()[completed % count(study)].values);
    }
    values
}

/// The first exact maximum is the deterministic governing operating case.
/// Ordinary one-case studies borrow their existing rows without allocation.
pub(super) fn completed(rows: &[Trial], case_count: usize) -> Cow<'_, [Trial]> {
    if case_count <= 1 { return Cow::Borrowed(rows); }
    Cow::Owned(rows.chunks_exact(case_count).map(|group| {
        group.iter().skip(1).fold(&group[0], |worst, row|
            if row.temperature_k > worst.temperature_k { row } else { worst }).clone()
    }).collect())
}

pub(super) fn scope(study: &UncertaintyStudy, original: &str) -> String {
    if study.operating_cases().is_empty() { return original.into(); }
    format!("{original} Each design coordinate is evaluated under every authored finite operating case. The monotonicity assumption and bracket apply to the maximum temperature across this fixed case set. Only complete case groups establish feasibility. Budgets count physical case evaluations; a partially evaluated candidate is never a feasible incumbent. The selected project is the actually evaluated governing-case project. These cases carry no probabilities and establish no interval, unlisted-condition or distributional robustness.")
}

pub(super) fn json(study: &UncertaintyStudy, rows: &[Trial], policy: DesignPolicy,
    chosen: Option<&Trial>) -> String
{
    if study.operating_cases().is_empty() { return String::new(); }
    let count = count(study);
    let declarations = study.operating_cases().iter().map(|case| {
        let inputs = study.parameters()[1..].iter().zip(&case.values).map(|(p, value)|
            format!("{{\"coordinate\":{},\"entity\":{},\"unit\":{},\"value\":{value}}}",
                quoted(&p.name), quoted(&p.entity), quoted(p.target.unit())))
            .collect::<Vec<_>>().join(",");
        format!("{{\"name\":{},\"inputs\":[{inputs}]}}", quoted(&case.name))
    }).collect::<Vec<_>>().join(",");
    let evaluations = rows.iter().enumerate().map(|(i, row)|
        format!("{{\"candidate\":{},\"operating_case\":{},\"observation\":{}}}",
            i / count, quoted(&study.operating_cases()[i % count].name), row.json(policy)))
        .collect::<Vec<_>>().join(",");
    let remainder = rows.len() % count;
    let pending = if remainder == 0 { "null".into() } else {
        format!("{{\"coordinate\":{},\"cases_completed\":{remainder},\"cases_required\":{count}}}",
            rows.last().expect("nonempty partial group").coordinate())
    };
    let governing = chosen.and_then(|chosen| rows.iter().position(|row| row == chosen))
        .map_or_else(|| "null".into(), |i| quoted(&study.operating_cases()[i % count].name));
    format!(",\"operating_cases\":[{declarations}],\"cases_per_candidate\":{count},\"candidates_completed\":{},\"case_evaluations\":[{evaluations}],\"pending_candidate\":{pending},\"selected_operating_case\":{governing},\"selected_project_scope\":\"evaluated-governing-operating-case\",\"robustness_scope\":\"authored-finite-cases-only\"",
        rows.len() / count)
}

pub(super) fn html(study: &UncertaintyStudy, rows: &[Trial]) -> String {
    if study.operating_cases().is_empty() { return String::new(); }
    let count = count(study);
    let mut html = format!("<h2>Finite operating cases</h2><p>{} complete candidates, {} physical case evaluations, {} cases in the pending candidate. Only complete candidates may be selected. The exported project is the evaluated governing case.</p><table><tr><th>Candidate</th><th>Operating case</th><th>Decision K</th></tr>",
        rows.len() / count, rows.len(), rows.len() % count);
    for (i, row) in rows.iter().enumerate() {
        let name = study.operating_cases()[i % count].name.replace('&', "&amp;")
            .replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
            .replace('\'', "&#39;");
        use std::fmt::Write as _;
        let _ = write!(html, "<tr><td>{}</td><td>{name}</td><td>{}</td></tr>",
            i / count, row.temperature_k);
    }
    html.push_str("</table>");
    html
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{Sample, DesignControl, Next, UniformParameter, next, selected};
    use fs_project::uncertainty::Target;
    use fs_blake3::hash_bytes;

    fn trial(coordinate: f64, case: f64, temperature: f64) -> Trial {
        Trial { sample: Sample { run: format!("{coordinate}/{case}"), project_hash: String::new(),
            qoi_receipt: hash_bytes(b"qoi"), value_k: temperature,
            parameters: vec![coordinate, case] }, temperature_k: temperature,
            time_s: None, state: None, conduction: hash_bytes(b"conduction") }
    }

    #[test]
    fn partial_case_groups_never_supply_a_feasible_candidate() {
        let rows = vec![trial(1.2, 0.0, 300.0), trial(1.2, 1.0, 310.0),
            trial(0.8, 0.0, 299.0)];
        let full = completed(&rows, 2);
        assert_eq!(full.len(), 1);
        assert_eq!(full[0], rows[1]);
        let policy = DesignPolicy { control: DesignControl::FanSpeedRatio { tolerance: 0.01 },
            temperature_tolerance_k: 0.01 };
        assert_eq!(selected(&full, policy, 315.0).unwrap().coordinate(), 1.2);
        assert!(completed(&rows[..1], 2).is_empty());
        assert_eq!(completed(&rows, 1).as_ref(), rows.as_slice());
    }

    #[test]
    fn governing_case_switches_and_exact_ties_keep_authored_order() {
        let rows = vec![trial(1.2, 0.0, 300.0), trial(1.2, 1.0, 310.0),
            trial(0.8, 0.0, 330.0), trial(0.8, 1.0, 320.0),
            trial(1.0, 0.0, 315.0), trial(1.0, 1.0, 315.0)];
        let full = completed(&rows, 2);
        assert_eq!(full.as_ref(), &[rows[1].clone(), rows[2].clone(), rows[4].clone()]);
    }

    #[test]
    fn worst_case_bracket_selects_the_feasible_design_not_a_passing_partial() {
        let policy = DesignPolicy { control: DesignControl::FanSpeedRatio { tolerance: 0.001 },
            temperature_tolerance_k: 0.001 };
        let control = UniformParameter { name:"fan".into(), target:Target::FanSpeedRatio,
            entity:"bank".into(), low:0.8, high:1.2 };
        let mut rows = Vec::new();
        for _ in 0..64 {
            let full = completed(&rows, 2);
            let coordinate = match next(&control, policy, 315.0, &full) {
                Next::Evaluate(value) => value,
                Next::Complete("target-bracketed") => break,
                other => panic!("unexpected decision {other:?}"),
            };
            rows.push(trial(coordinate, 0.0, 285.0 + 20.0 / coordinate));
            let before = completed(&rows, 2);
            assert_eq!(before.len(), (rows.len() - 1) / 2);
            rows.push(trial(coordinate, 1.0, 295.0 + 20.0 / coordinate));
        }
        let full = completed(&rows, 2);
        assert_eq!(next(&control, policy, 315.0, &full), Next::Complete("target-bracketed"));
        let chosen = selected(&full, policy, 315.0).unwrap();
        assert!((1.0..=1.001).contains(&chosen.coordinate()));
        assert_eq!(chosen.sample.parameters[1], 1.0);
    }
}
