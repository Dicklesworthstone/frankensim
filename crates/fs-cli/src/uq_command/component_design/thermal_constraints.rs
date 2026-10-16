//! Independent temperature limits on one fixed physical experiment.
//! An envelope compares T_i - limit_i, NOT the largest absolute temperature.
//! It is a sampled-model decision and never a continuous-time certificate.
use super::*;
use std::collections::BTreeSet;

pub(super) const MAX_CONSTRAINTS: usize = 16;

#[derive(Debug, Clone)]
pub(super) struct Constraint {
    pub name: String,
    pub limit: f64,
    pub objective: J,
}

impl Constraint {
    pub fn primary(base: &J, limit: f64) -> Result<Self> {
        Ok(Self { name: "primary".into(), limit, objective: field(base, "objective")?.clone() })
    }

    /// Only the observation and threshold change. The physical experiment,
    /// fixed grid, source controls and requested adjoint policy stay identical.
    pub fn request(&self, base: &J) -> Result<J> {
        let mut request = base.clone();
        input::put(&mut request, "objective", self.objective.clone())?;
        input::put(input::member_mut(&mut request, "transient")?,
            "temperature_limit_k", number_value(self.limit)?)?;
        Ok(request)
    }
}

/// Extra limits are optional. The original objective/limit remains mandatory
/// and is always first. Sixteen includes that original primary constraint.
pub(super) fn parse(base: &J, value: Option<&J>) -> Result<Vec<Constraint>> {
    let Some(value) = value else { return Ok(Vec::new()); };
    let rows = array(value, "thermal_constraints", MAX_CONSTRAINTS - 1)?;
    if rows.is_empty() { return Err(bad("thermal_constraints must be nonempty when supplied")); }
    let mut names = BTreeSet::from(["primary".to_string()]);
    let mut constraints = Vec::with_capacity(rows.len());
    for row in rows {
        object(row, &["name", "temperature_limit_k", "component", "objective"], "thermal constraint")?;
        let name = string(field(row, "name")?, "thermal constraint name")?;
        if !names.insert(name.clone()) { return Err(bad("thermal constraint names must be unique; primary is reserved")); }
        let limit = positive(field(row, "temperature_limit_k")?, "thermal constraint temperature_limit_k")?;
        let objective = match (row.get("component"), row.get("objective")) {
            (Some(component), None) => component_objective(base, string(component, "thermal constraint component")?.as_str())?,
            (None, Some(objective)) => admitted_objective(base, objective)?,
            _ => return Err(bad("each thermal constraint requires exactly one component or objective selector")),
        };
        constraints.push(Constraint { name, limit, objective });
    }
    Ok(constraints)
}

fn vertex_count(base: &J) -> Result<usize> {
    let vertices = field(field(base, "solid")?, "vertices_m")?.as_array()
        .ok_or_else(|| bad("temperature constraints require the declared solid vertices"))?;
    if vertices.is_empty() { return Err(bad("temperature constraints cannot select an empty mesh")); }
    Ok(vertices.len())
}

fn component_objective(base: &J, name: &str) -> Result<J> {
    let rows = field(field(field(base, "solid")?, "component_power")?, "components")?
        .as_array().ok_or_else(|| bad("thermal component limit requires the original component map"))?;
    let matches: Vec<_> = rows.iter().filter(|row| row.str_field("name") == Some(name)).collect();
    if matches.len() != 1 { return Err(bad("thermal constraint names an unknown or ambiguous component")); }
    let count = vertex_count(base)?;
    let mut selected = BTreeSet::new();
    for vertex in array(field(matches[0], "vertices")?, "component footprint", count)? {
        selected.insert(integer(vertex, "component temperature vertex", count - 1)?);
    }
    // PowerMap also sorts/deduplicates a footprint. This does not invent a
    // junction temperature or include unrelated nodes outside that footprint.
    let selector = J::Object(vec![("max_vertices".into(),
        J::Array(selected.into_iter().map(usize_value).collect()))]);
    admitted_objective(base, &selector)
}

fn admitted_objective(base: &J, selector: &J) -> Result<J> {
    const KEYS: [&str; 4] = ["mean_wall_region", "max_wall_region", "max_solid_temperature", "max_vertices"];
    object(selector, &KEYS, "thermal constraint objective")?;
    let keys: Vec<_> = KEYS.into_iter().filter(|key| selector.get(key).is_some()).collect();
    if keys.len() != 1 { return Err(bad("thermal constraint objective requires exactly one temperature selector")); }
    let key = keys[0];
    match key {
        "max_solid_temperature" => {
            if field(selector, key)? != &J::Bool(true) { return Err(bad("max_solid_temperature must be true")); }
        }
        "max_vertices" => {
            let count = vertex_count(base)?;
            let rows = array(field(selector, key)?, "constraint vertices", count)?;
            if rows.is_empty() { return Err(bad("a thermal vertex limit requires a nonempty observation set")); }
            let mut seen = BTreeSet::new();
            for row in rows {
                if !seen.insert(integer(row, "thermal constraint vertex", count - 1)?) {
                    return Err(bad("thermal constraint repeats an observation vertex"));
                }
            }
        }
        _ => {
            let name = string(field(selector, key)?, "thermal constraint surface")?;
            let rows = field(field(base, "solid")?, "surfaces")?.as_array()
                .ok_or_else(|| bad("thermal wall limit requires declared surfaces"))?;
            if rows.iter().filter(|row| row.str_field("name") == Some(name.as_str())).count() != 1 {
                return Err(bad("thermal constraint names an unknown or ambiguous cooling surface"));
            }
        }
    }
    let mut objective = selector.clone();
    input::put(&mut objective, "gradient", J::Bool(false))?;
    Ok(objective)
}

pub(super) struct Assessment {
    pub constraint: Constraint,
    pub peak: f64,
    pub time: f64,
    pub slopes: Vec<Option<f64>>,
}
impl Assessment {
    fn margin(&self) -> Result<f64> { checked(self.peak - self.constraint.limit) }
    fn render(&self) -> Result<J> {
        Ok(J::Object(vec![
            ("name".into(), J::Str(self.constraint.name.clone())),
            ("objective".into(), self.constraint.objective.clone()),
            ("temperature_limit_k".into(), number_value(self.constraint.limit)?),
            ("sampled_peak_k".into(), number_value(self.peak)?),
            ("sampled_peak_time_s".into(), number_value(self.time)?),
            ("temperature_excess_k".into(), number_value(self.margin()?)?),
            ("passing".into(), J::Bool(self.peak <= self.constraint.limit)),
            ("dpeak_dcontrolled_power_w_k_per_w".into(), J::Array(self.slopes.iter()
                .map(|value| value.map(number_value).transpose().map(|v| v.unwrap_or(J::Null)))
                .collect::<Result<Vec<_>>>()?)),
        ]))
    }
}

pub(super) struct Envelope {
    rows: Vec<Assessment>,
    active: usize,
    pub margin: f64,
    pub ties: usize,
}
impl Envelope {
    pub fn new(rows: Vec<Assessment>) -> Result<Self> {
        if rows.is_empty() || rows.len() > MAX_CONSTRAINTS { return Err(model_failure("incomplete thermal constraint envelope")); }
        let controls = rows[0].slopes.len();
        let mut active = 0;
        let mut margin = f64::NEG_INFINITY;
        let mut ties = 0;
        for (index, row) in rows.iter().enumerate() {
            if !(row.peak.is_finite() && row.peak > 0.0 && row.time.is_finite() && row.time >= 0.0
                && row.constraint.limit.is_finite() && row.constraint.limit > 0.0)
                || row.slopes.len() != controls
                || row.slopes.iter().flatten().any(|v| !v.is_finite()) {
                return Err(model_failure("invalid complete-trajectory thermal constraint result"));
            }
            let next = row.margin()?;
            if next > margin { margin = next; active = index; ties = 1; }
            else if next == margin { ties += 1; }
        }
        Ok(Self { rows, active, margin, ties })
    }
    pub fn active_name(&self) -> &str { &self.rows[self.active].constraint.name }
    pub fn slopes(&self) -> Vec<Option<f64>> {
        if self.ties == 1 { self.rows[self.active].slopes.clone() }
        else { vec![None; self.rows[self.active].slopes.len()] }
    }
    pub fn report(&self) -> Result<J> {
        Ok(J::Object(vec![
            ("all_passed".into(), J::Bool(self.margin <= 0.0)),
            ("maximum_temperature_excess_k".into(), number_value(self.margin)?),
            ("active_constraint".into(), J::Str(self.active_name().into())),
            ("exact_active_count".into(), usize_value(self.ties)),
            ("rows".into(), J::Array(self.rows.iter().map(Assessment::render).collect::<Result<Vec<_>>>()?)),
            ("scope".into(), J::Str("all initial/accepted-endpoint peaks obey their own limits; each criterion observes the same fixed complete trajectory; the largest absolute temperature need not be the limiting constraint; active-margin derivatives only guide verified trials, exact margin ties use bisection; component selectors observe original footprint vertices, not certified physical junction temperatures; no continuous-time or physical-compliance certificate".into())),
        ]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn row(name: &str, peak: f64, limit: f64, slope: f64) -> Assessment {
        Assessment { constraint: Constraint { name: name.into(), limit, objective: J::Null },
            peak, time: 2.0, slopes: vec![Some(slope)] }
    }
    #[test]
    fn colder_component_can_be_the_binding_temperature_limit() {
        let e = Envelope::new(vec![row("primary", 330.0, 350.0, 10.0), row("memory", 305.0, 304.0, 2.0)]).unwrap();
        assert_eq!(e.margin, 1.0); assert_eq!(e.active_name(), "memory");
        assert_eq!(e.slopes(), vec![Some(2.0)]);
        let e = Envelope::new(vec![row("primary", 330.0, 350.0, 10.0), row("memory", 303.0, 304.0, 2.0)]).unwrap();
        assert_eq!(e.margin, -1.0);
    }
    #[test]
    fn active_switches_and_exact_ties_never_reuse_a_stale_constraint_slope() {
        let a = Envelope::new(vec![row("a", 301.0, 300.0, 1.0), row("b", 301.5, 300.0, 5.0)]).unwrap();
        assert_eq!(a.slopes(), vec![Some(5.0)]);
        let b = Envelope::new(vec![row("a", 302.0, 300.0, 1.0), row("b", 301.5, 300.0, 5.0)]).unwrap();
        assert_eq!(b.slopes(), vec![Some(1.0)]);
        let tied = Envelope::new(vec![row("a", 301.0, 300.0, 1.0), row("b", 305.0, 304.0, 5.0)]).unwrap();
        assert_eq!(tied.ties, 2); assert_eq!(tied.slopes(), vec![None]);
        assert!(Envelope::new(vec![row("bad", f64::NAN, 300.0, 1.0)]).is_err());
    }
    #[test]
    fn component_selectors_preserve_physics_and_refuse_ambiguous_declarations() {
        let base = J::parse(r#"{"solid":{"vertices_m":[[0,0,0],[1,0,0]],"component_power":{"components":[{"name":"memory","vertices":[1,0,1]}]},"surfaces":[{"name":"wall"}]},"objective":{"max_solid_temperature":true,"gradient":false},"transient":{"temperature_limit_k":350}}"#).unwrap();
        let spec = J::parse(r#"[{"name":"memory-limit","component":"memory","temperature_limit_k":304}]"#).unwrap();
        // A footprint may repeat IDs just as PowerMap permits; resolve the
        // sorted observation set without changing the original component map.
        let mut expanded = base.clone();
        let components = input::member_mut(input::member_mut(input::member_mut(&mut expanded,"solid").unwrap(),"component_power").unwrap(),"components").unwrap();
        if let J::Array(rows) = components { input::put(&mut rows[0],"vertices",J::parse("[1,0]").unwrap()).unwrap(); }
        let c = parse(&expanded, Some(&spec)).unwrap();
        let request = c[0].request(&expanded).unwrap();
        assert_eq!(request.get("solid"), expanded.get("solid"));
        assert_eq!(request.path(&["objective","max_vertices"]), Some(&J::parse("[0,1]").unwrap()));
        for text in [r#"[]"#, r#"[{"name":"primary","component":"memory","temperature_limit_k":304}]"#,
            r#"[{"name":"x","component":"missing","temperature_limit_k":304}]"#,
            r#"[{"name":"x","objective":{"max_vertices":[2]},"temperature_limit_k":304}]"#,
            r#"[{"name":"x","objective":{"max_vertices":[1,1]},"temperature_limit_k":304}]"#,
            r#"[{"name":"x","objective":{"max_solid_temperature":true,"gradient":true},"temperature_limit_k":304}]"#] {
            assert!(parse(&expanded, Some(&J::parse(text).unwrap())).is_err(), "{text}");
        }
    }
}
