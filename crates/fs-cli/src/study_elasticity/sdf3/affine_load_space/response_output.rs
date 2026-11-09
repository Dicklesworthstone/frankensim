//! Original signed response limits and their independent native adversaries.
use super::super::super::spec::responses::{Family, Limit, Observation};
use super::*;
use fs_topopt::sdf3::response::limits::ResponseLimitSense3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::response::{
    AffineResponseEvaluation3, AffineResponseOptions3,
};

pub(super) const SCOPE: &str = "Each signed response inequality applies to the complete affine load ball and retains its own adverse coefficients. Mean displacement uses the fixed reference cut-volume average on the declared patch. Force and moment are the native embedded numerical support traction tested by a Cartesian translation or rotation mode. Physical extrema and independently evaluated combined-field witnesses are reported with the original signed bounds and same-unit scales. Response and joint aggregates are conservative in exact arithmetic; these are estimated discrete linear-model results with numerical identity diagnostics.";

fn invalid(message: impl Into<String>) -> Failure {
    fail("cli-study-sdf3-response-output", message)
}

fn sense(value: ResponseLimitSense3) -> &'static str {
    match value {
        ResponseLimitSense3::AtMost => "at-most",
        ResponseLimitSense3::AtLeast => "at-least",
    }
}

fn declaration(limit: &Limit) -> Result<String> {
    let axis = match limit.observation.axis() {
        0 => "x",
        1 => "y",
        2 => "z",
        _ => return Err(invalid("invalid response observation axis")),
    };
    let detail = match limit.observation {
        Observation::MeanDisplacement { x_fraction, .. } => format!(
            ",\"x_fraction\":{x_fraction:?},\"normalization\":\"reference-cut-volume-average\"",
        ),
        Observation::SupportForce { .. } => String::new(),
        Observation::SupportMoment { about_m, .. } => format!(",\"about_m\":{about_m:?}"),
    };
    Ok(format!(
        "\"name\":{},\"kind\":{:?},\"unit\":{:?},\"axis\":{axis:?}{detail},\"sense\":{:?},\"bound\":{:.17e},\"scale\":{:.17e}",
        quoted(&limit.name),
        limit.observation.kind(),
        limit.observation.unit(),
        sense(limit.sense),
        limit.bound,
        limit.scale,
    ))
}

pub(super) fn options(options: AffineResponseOptions3) -> String {
    format!(
        "{{\"sharpness\":{:.17e},\"max_limits\":{}}}",
        options.sharpness, options.max_limits
    )
}

pub(super) fn model(family: &Family) -> Result<String> {
    let limits = family
        .limits
        .iter()
        .map(|limit| declaration(limit).map(|row| format!("{{{row}}}")))
        .collect::<Result<Vec<_>>>()?
        .join(",");
    Ok(format!(
        "{{\"options\":{},\"limits\":[{limits}]}}",
        options(family.options)
    ))
}

fn admit(family: &Family, evaluation: &AffineResponseEvaluation3) -> Result<()> {
    if family.options != evaluation.options
        || family.limits.is_empty()
        || family.limits.len() != evaluation.witnesses.len()
        || family
            .limits
            .iter()
            .zip(&evaluation.witnesses)
            .any(|(limit, witness)| {
                limit.sense != witness.sense
                    || limit.bound != witness.bound
                    || limit.scale != witness.scale
            })
    {
        return Err(invalid(
            "retained response witnesses differ from their ordered source limits",
        ));
    }
    Ok(())
}

pub(super) fn measurement(family: &Family, e: &AffineResponseEvaluation3) -> Result<String> {
    admit(family, e)?;
    let rows = family.limits.iter().zip(&e.witnesses).enumerate()
        .map(|(index, (limit, w))| -> Result<String> {
            Ok(format!(
                "{{\"index\":{index},{},\"baseline\":{:.17e},\"basis_responses\":{:?},\"radius\":{:.17e},\"coefficients\":{:?},\"worst_value\":{:.17e},\"witness_value\":{:.17e},\"residual\":{:.17e},\"witness_residual\":{:.17e},\"relative_defect\":{:.17e}}}",
                declaration(limit)?, w.baseline, w.basis_responses, w.radius,
                w.coefficients, w.worst_value, w.witness_value, w.residual,
                w.witness_residual, w.relative_defect,
            ))
        }).collect::<Result<Vec<_>>>()?.join(",");
    Ok(format!(
        "{{\"constraint\":{:.17e},\"worst_residual\":{:.17e},\"feasible\":{},\"gradient_available\":{},\"options\":{},\"adjoint_layout\":{:?},\"adjoint_count\":{},\"limits\":[{rows}]}}",
        e.constraint,
        e.worst_residual,
        e.physically_feasible(),
        e.gradient.is_some(),
        options(e.options),
        format!("{:?}", e.adjoint_layout),
        e.adjoints.len(),
    ))
}

/// The already retained offset and basis fields reconstruct each physical
/// response witness; the report keeps its scalar measurements and bounds.
pub(super) fn fields(family: &Family, e: &AffineResponseEvaluation3) -> Result<String> {
    admit(family, e)?;
    let rows = family.limits.iter().zip(&e.witnesses).enumerate()
        .map(|(index, (limit, w))| format!(
            "{{\"index\":{index},\"name\":{},\"kind\":{:?},\"unit\":{:?},\"sense\":{:?},\"coefficients\":{:?}}}",
            quoted(&limit.name), limit.observation.kind(), limit.observation.unit(),
            sense(limit.sense), w.coefficients,
        )).collect::<Vec<_>>().join(",");
    Ok(format!(
        "{{\"schema\":\"sdf3-affine-response-fields.v1\",\"witnesses\":[{rows}]}}"
    ))
}

fn html_text(value: &str) -> String {
    let mut escaped = String::new();
    for ch in value.chars() {
        match ch {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            '\'' => escaped.push_str("&#39;"),
            _ => escaped.push(ch),
        }
    }
    escaped
}

pub(super) fn html(
    family: &Family,
    selected: &AffineResponseEvaluation3,
    accepted: &AffineResponseEvaluation3,
) -> Result<String> {
    admit(family, selected)?;
    admit(family, accepted)?;
    let rows = family.limits.iter().zip(&selected.witnesses).zip(&accepted.witnesses)
        .map(|((limit, chosen), current)| format!(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td></tr>",
            html_text(&limit.name), limit.observation.unit(), sense(limit.sense), limit.bound,
            chosen.worst_value, chosen.witness_value, chosen.residual, current.worst_value,
        )).collect::<Vec<_>>().join("");
    Ok(format!(
        "<h2>Signed response limits over the affine load family</h2><p>Selected response aggregate: {:.8e}; physical limits and aggregate pass: {}. Last accepted response aggregate: {:.8e}; physical limits and aggregate pass: {}. Nonpositive normalized residuals satisfy their original signed inequalities; the full joint constraint remains separate.</p><table><tr><th>Limit</th><th>Unit</th><th>Sense</th><th>Bound</th><th>Selected adverse value</th><th>Selected measured witness</th><th>Selected normalized residual</th><th>Last accepted adverse value</th></tr>{rows}</table><p>{SCOPE}</p>",
        selected.constraint,
        selected.physically_feasible(),
        accepted.constraint,
        accepted.physically_feasible(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_topopt::sdf3::robust::envelope::affine_load_space::response::{
        AffineResponseAdjointLayout3, AffineResponseWitness3,
    };

    #[test]
    fn retained_response_rows_keep_signed_units_order_and_independent_coefficients() {
        let family = Family {
            options: AffineResponseOptions3::default(),
            limits: vec![
                Limit {
                    name: "motion-lower".into(),
                    observation: Observation::MeanDisplacement {
                        axis: 2,
                        x_fraction: [0.5, 1.0],
                    },
                    sense: ResponseLimitSense3::AtLeast,
                    bound: -2.0,
                    scale: 2.0,
                },
                Limit {
                    name: "support-upper".into(),
                    observation: Observation::SupportMoment {
                        axis: 1,
                        about_m: [1.0, 2.0, 3.0],
                    },
                    sense: ResponseLimitSense3::AtMost,
                    bound: 2.0,
                    scale: 2.0,
                },
            ],
        };
        // Exercise the serialization boundary with explicitly authored native
        // records; the solver tests establish these response values themselves.
        let witnesses = family
            .limits
            .iter()
            .enumerate()
            .map(|(index, limit)| {
                let adverse = if index == 0 { -1.0 } else { 1.0 };
                AffineResponseWitness3 {
                    baseline: 0.0,
                    basis_responses: vec![1.0, 0.0],
                    radius: 1.0,
                    coefficients: vec![adverse, 0.0],
                    worst_value: adverse,
                    witness_value: adverse,
                    bound: limit.bound,
                    sense: limit.sense,
                    scale: limit.scale,
                    residual: -0.5,
                    witness_residual: -0.5,
                    relative_defect: 0.0,
                }
            })
            .collect();
        let mut evaluation = AffineResponseEvaluation3 {
            options: family.options,
            constraint: -0.4,
            gradient: None,
            worst_residual: -0.5,
            witnesses,
            adjoint_layout: AffineResponseAdjointLayout3::NotRequested,
            adjoints: Vec::new(),
            work: SolveWork::default(),
        };
        let report = JsonValue::parse(&measurement(&family, &evaluation).unwrap()).unwrap();
        let rows = report.get("limits").unwrap().as_array().unwrap();
        assert_eq!(rows[0].str_field("unit"), Some("m"));
        assert_eq!(rows[0].str_field("sense"), Some("at-least"));
        assert_eq!(rows[0].get("bound").unwrap().as_f64(), Some(-2.0));
        assert_eq!(rows[0].get("worst_value").unwrap().as_f64(), Some(-1.0));
        assert_eq!(rows[1].str_field("unit"), Some("N*m"));
        assert_eq!(rows[1].str_field("axis"), Some("y"));
        assert_eq!(rows[1].get("about_m").unwrap().as_array().unwrap().len(), 3);
        assert_ne!(rows[0].get("coefficients"), rows[1].get("coefficients"));
        let design = JsonValue::parse(&fields(&family, &evaluation).unwrap()).unwrap();
        let retained = design.get("witnesses").unwrap().as_array().unwrap();
        for (row, retained) in rows.iter().zip(retained) {
            assert_eq!(retained.get("coefficients"), row.get("coefficients"));
            assert_eq!(retained.get("name"), row.get("name"));
        }
        evaluation.witnesses[0].bound = 2.0;
        assert!(measurement(&family, &evaluation).is_err());
        assert!(fields(&family, &evaluation).is_err());
    }
}
