//! `cooling-cht` design studies: a scene with a `study` block is evaluated
//! over the Cartesian product of declared parameter values (each a JSON
//! path into the scene and a list of values), every variant solved as an
//! ordinary scene under its own wall budget, and the variants ranked by a
//! declared objective subject to optional bounds on other result
//! quantities. A variant the solver refuses (not steady, a fin the grid
//! cannot represent, ...) is recorded with its refusal and never ranked.
//!
//! No-claims: a grid search over the declared values (no optimality claim
//! between grid points); every evaluation inherits the `cooling-cht`
//! no-claims; no sensitivity or uncertainty is attached to the ranking.

use std::fmt::Write as _;
use std::time::Instant;

use super::json::JsonValue as J;
use super::{NO_CLAIM, Result, Scene, bad, execute_within_budget, num, quote};

const STUDY_SCHEMA: &str = "frankensim.cooling-cht.study.v1";
const MAX_EVALUATIONS: usize = 64;

/// One step of a parameter path: an object key or an array index.
enum Step {
    Key(String),
    Index(usize),
}

struct Parameter {
    name: String,
    path: Vec<Step>,
    values: Vec<J>,
}

/// A bound on a result quantity.
struct Bound {
    quantity: String,
    min: Option<f64>,
    max: Option<f64>,
}

fn parse_path(item: &J, at: &str) -> Result<Vec<Step>> {
    let steps = item
        .get("path")
        .and_then(J::as_array)
        .filter(|steps| !steps.is_empty())
        .ok_or_else(|| {
            bad(format!(
                "{at}.path must be a non-empty array of keys and indices"
            ))
        })?;
    steps
        .iter()
        .map(|step| match step {
            J::Str(key) => Ok(Step::Key(key.clone())),
            J::Number { value, .. } if *value >= 0.0 && value.fract() == 0.0 => {
                Ok(Step::Index(*value as usize))
            }
            _ => Err(bad(format!(
                "{at}.path entries must be object keys or array indices"
            ))),
        })
        .collect()
}

/// Replace the value at `path` (which must exist) by `value`.
fn set_at(node: &mut J, path: &[Step], value: J, at: &str) -> Result<()> {
    let Some((step, rest)) = path.split_first() else {
        *node = value;
        return Ok(());
    };
    let child = match (step, node) {
        (Step::Key(key), J::Object(entries)) => {
            entries.iter_mut().find(|(k, _)| k == key).map(|(_, v)| v)
        }
        (Step::Index(i), J::Array(items)) => items.get_mut(*i),
        _ => None,
    };
    let child = child.ok_or_else(|| bad(format!("{at}: the path does not exist in the scene")))?;
    set_at(child, rest, value, at)
}

/// A named result quantity: `max_solid_temperature_k`, `fan_flow_m3_s`,
/// `inflow_m3_s`, `source:<name>` (its peak temperature),
/// `component:<name>` (its junction temperature), or
/// `internal_fan:<name>` (its flow).
fn quantity(result: &J, name: &str) -> Option<f64> {
    let named = |list: &str, key: &str, field: &str| {
        result
            .get(list)
            .and_then(J::as_array)?
            .iter()
            .find(|item| item.str_field("name") == Some(key))?
            .get(field)
            .and_then(J::as_f64)
    };
    if let Some(source) = name.strip_prefix("source:") {
        return named("sources", source, "max_temperature_k");
    }
    if let Some(part) = name.strip_prefix("component:") {
        return named("components", part, "junction_temperature_k");
    }
    if let Some(fan) = name.strip_prefix("internal_fan:") {
        return result
            .path(&["flow", "internal_fans"])
            .and_then(J::as_array)?
            .iter()
            .find(|item| item.str_field("name") == Some(fan))?
            .get("flow_m3_s")
            .and_then(J::as_f64);
    }
    match name {
        "max_solid_temperature_k" => result.get(name).and_then(J::as_f64),
        "fan_flow_m3_s" | "inflow_m3_s" => result.path(&["flow", name]).and_then(J::as_f64),
        _ => None,
    }
}

fn known_quantity(name: &str) -> bool {
    ["source:", "component:", "internal_fan:"]
        .iter()
        .any(|prefix| {
            name.strip_prefix(prefix)
                .is_some_and(|rest| !rest.is_empty())
        })
        || matches!(
            name,
            "max_solid_temperature_k" | "fan_flow_m3_s" | "inflow_m3_s"
        )
}

/// One evaluated variant.
struct Evaluation {
    values: Vec<J>,
    outcome: std::result::Result<(f64, Vec<Option<f64>>), (String, String)>,
}

impl Evaluation {
    fn feasible(&self, bounds: &[Bound]) -> bool {
        match &self.outcome {
            Ok((_, measured)) => bounds.iter().zip(measured).all(|(bound, value)| {
                value.is_some_and(|v| {
                    bound.min.is_none_or(|m| v >= m) && bound.max.is_none_or(|m| v <= m)
                })
            }),
            Err(_) => false,
        }
    }
}

fn json_text(value: &J) -> String {
    match value {
        J::Null => "null".into(),
        J::Bool(flag) => flag.to_string(),
        J::Number { raw, .. } => raw.clone(),
        J::Str(text) => quote(text),
        J::Array(items) => format!(
            "[{}]",
            items.iter().map(json_text).collect::<Vec<_>>().join(",")
        ),
        J::Object(entries) => format!(
            "{{{}}}",
            entries
                .iter()
                .map(|(k, v)| format!("{}:{}", quote(k), json_text(v)))
                .collect::<Vec<_>>()
                .join(",")
        ),
    }
}

/// Run the study declared in `root`.
#[allow(clippy::too_many_lines)] // admission, loop, one report
pub(super) fn run(root: &J, base: &std::path::Path, json_mode: bool) -> Result<String> {
    let study = root.get("study").expect("checked by the caller");
    let mut parameters = Vec::new();
    for (i, item) in study
        .get("parameters")
        .and_then(J::as_array)
        .filter(|items| !items.is_empty())
        .ok_or_else(|| bad("study.parameters must be a non-empty array"))?
        .iter()
        .enumerate()
    {
        let at = format!("study.parameters[{i}]");
        let values = item
            .get("values")
            .and_then(J::as_array)
            .filter(|values| !values.is_empty())
            .ok_or_else(|| bad(format!("{at}.values must be a non-empty array")))?
            .to_vec();
        parameters.push(Parameter {
            name: item.str_field("name").unwrap_or("parameter").to_string(),
            path: parse_path(item, &at)?,
            values,
        });
    }
    let evaluations: usize = parameters.iter().map(|p| p.values.len()).product();
    if evaluations > MAX_EVALUATIONS {
        return Err(bad(format!(
            "study declares {evaluations} evaluations; at most {MAX_EVALUATIONS}"
        )));
    }
    let objective = study
        .get("objective")
        .and_then(|o| o.str_field("minimize"))
        .filter(|name| known_quantity(name))
        .ok_or_else(|| {
            bad("study.objective.minimize must name a result quantity (max_solid_temperature_k, source:<name>, component:<name>, internal_fan:<name>, fan_flow_m3_s, inflow_m3_s)")
        })?
        .to_string();
    let mut bounds = Vec::new();
    if let Some(items) = study.get("constraints") {
        for (i, item) in items
            .as_array()
            .ok_or_else(|| bad("study.constraints must be an array"))?
            .iter()
            .enumerate()
        {
            let at = format!("study.constraints[{i}]");
            let name = item
                .str_field("quantity")
                .filter(|name| known_quantity(name))
                .ok_or_else(|| bad(format!("{at}.quantity must name a result quantity")))?;
            let (min, max) = (item.f64_field("min"), item.f64_field("max"));
            if min.is_none() && max.is_none() {
                return Err(bad(format!("{at} needs min or max")));
            }
            bounds.push(Bound {
                quantity: name.to_string(),
                min,
                max,
            });
        }
    }
    // The base scene without its study block.
    let mut template = root.clone();
    if let J::Object(entries) = &mut template {
        entries.retain(|(key, _)| key != "study");
    }
    let started = Instant::now();
    let mut results: Vec<Evaluation> = Vec::with_capacity(evaluations);
    for index in 0..evaluations {
        // Mixed-radix digits, the first parameter slowest.
        let mut rest = index;
        let mut picks = vec![0usize; parameters.len()];
        for (slot, parameter) in picks.iter_mut().zip(&parameters).rev() {
            *slot = rest % parameter.values.len();
            rest /= parameter.values.len();
        }
        let mut variant = template.clone();
        let values: Vec<J> = parameters
            .iter()
            .zip(&picks)
            .map(|(p, &k)| p.values[k].clone())
            .collect();
        for (parameter, value) in parameters.iter().zip(&values) {
            set_at(
                &mut variant,
                &parameter.path,
                value.clone(),
                &format!("study.parameters.{}", parameter.name),
            )?;
        }
        let outcome = Scene::from_root(&variant, base)
            .and_then(|scene| execute_within_budget(&scene, true))
            .and_then(|text| {
                J::parse(text.trim()).map_err(|e| bad(format!("unreadable result: {e:?}")))
            })
            .map_err(|failure| (failure.code.to_string(), failure.message))
            .and_then(|result| {
                let value = quantity(&result, &objective).ok_or_else(|| {
                    (
                        "cooling-cht-study".to_string(),
                        format!("the result has no {objective}"),
                    )
                })?;
                Ok((
                    value,
                    bounds
                        .iter()
                        .map(|b| quantity(&result, &b.quantity))
                        .collect(),
                ))
            });
        results.push(Evaluation { values, outcome });
    }
    let best = results
        .iter()
        .enumerate()
        .filter(|(_, e)| e.feasible(&bounds))
        .filter_map(|(i, e)| e.outcome.as_ref().ok().map(|(v, _)| (i, *v)))
        .min_by(|a, b| a.1.total_cmp(&b.1));
    let wall_s = started.elapsed().as_secs_f64();
    let names: Vec<&str> = parameters.iter().map(|p| p.name.as_str()).collect();
    if json_mode {
        let mut out = format!(
            "{{\"schema\":{},\"status\":\"completed\",\"objective\":{{\"minimize\":{}}},\"parameters\":[{}],\"evaluations\":[",
            quote(STUDY_SCHEMA),
            quote(&objective),
            names.iter().map(|n| quote(n)).collect::<Vec<_>>().join(",")
        );
        for (i, evaluation) in results.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let values = evaluation
                .values
                .iter()
                .map(json_text)
                .collect::<Vec<_>>()
                .join(",");
            match &evaluation.outcome {
                Ok((value, measured)) => {
                    let measured = measured
                        .iter()
                        .map(|v| v.map_or_else(|| "null".to_string(), |v| v.to_string()))
                        .collect::<Vec<_>>()
                        .join(",");
                    let _ = write!(
                        out,
                        "{{\"index\":{i},\"values\":[{values}],\"status\":\"completed\",\"objective\":{},\"constraints\":[{measured}],\"feasible\":{}}}",
                        num(*value)?,
                        evaluation.feasible(&bounds)
                    );
                }
                Err((code, message)) => {
                    let _ = write!(
                        out,
                        "{{\"index\":{i},\"values\":[{values}],\"status\":\"refused\",\"code\":{},\"message\":{},\"feasible\":false}}",
                        quote(code),
                        quote(message)
                    );
                }
            }
        }
        out.push(']');
        match best {
            Some((i, value)) => {
                let _ = write!(
                    out,
                    ",\"best\":{{\"index\":{i},\"values\":[{}],\"objective\":{}}}",
                    results[i]
                        .values
                        .iter()
                        .map(json_text)
                        .collect::<Vec<_>>()
                        .join(","),
                    num(value)?
                );
            }
            None => out.push_str(",\"best\":null"),
        }
        let _ = writeln!(
            out,
            ",\"wall_s\":{},\"evidence\":\"Estimated\",\"no_claim\":{}}}",
            num(wall_s)?,
            quote(&format!(
                "grid search over the declared values (no optimality claim between grid points); each evaluation: {NO_CLAIM}"
            ))
        );
        Ok(out)
    } else {
        let mut out = format!(
            "status=completed\nobjective=minimize {objective}\nevaluations={}\n",
            results.len()
        );
        for (i, evaluation) in results.iter().enumerate() {
            let values = names
                .iter()
                .zip(&evaluation.values)
                .map(|(n, v)| format!("{n}={}", json_text(v)))
                .collect::<Vec<_>>()
                .join(" ");
            match &evaluation.outcome {
                Ok((value, _)) => {
                    let _ = writeln!(
                        out,
                        "evaluation={i} {values} objective={value} feasible={}",
                        evaluation.feasible(&bounds)
                    );
                }
                Err((code, message)) => {
                    let _ = writeln!(out, "evaluation={i} {values} refused={code}: {message}");
                }
            }
        }
        match best {
            Some((i, value)) => {
                let _ = writeln!(out, "best={i} objective={value}");
            }
            None => out.push_str("best=none\n"),
        }
        let _ = writeln!(out, "wall_s={wall_s:.3}\nevidence=Estimated");
        Ok(out)
    }
}
