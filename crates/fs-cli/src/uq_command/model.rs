use super::*;
use std::collections::BTreeSet;

mod material;

#[derive(Debug, Clone)]
pub(super) enum Target {
    AirDensity,
    AirSpecificHeat,
    InletTemperature(usize),
    FanSpeedRatio,
    SurfaceHtc(String),
    ComponentPower(String),
    Material(material::Target),
}
impl Target {
    fn name(&self) -> String {
        match self {
            Self::AirDensity => "air.density_kg_m3".into(),
            Self::AirSpecificHeat => "air.specific_heat_j_kg_k".into(),
            Self::InletTemperature(index) => format!("inlet[{index}].temperature_k"),
            Self::FanSpeedRatio => "hydraulics.fan.speed_ratio".into(),
            Self::SurfaceHtc(name) => format!("surface[{name}].htc_w_m2_k"),
            Self::ComponentPower(name) => format!("component[{name}].power_w"),
            Self::Material(target) => target.name(),
        }
    }
    fn unit(&self) -> &'static str {
        match self {
            Self::AirDensity => "kg/m3",
            Self::AirSpecificHeat => "J/(kg K)",
            Self::InletTemperature(_) => "K",
            Self::FanSpeedRatio => "1",
            Self::SurfaceHtc(_) => "W/(m2 K)",
            Self::ComponentPower(_) => "W",
            Self::Material(_) => "W/(m K)",
        }
    }
    fn allows_zero(&self) -> bool { matches!(self, Self::ComponentPower(_)) }
    fn render(&self) -> String {
        match self {
            Self::AirDensity => "{\"kind\":\"air-density\"}".into(),
            Self::AirSpecificHeat => "{\"kind\":\"air-specific-heat\"}".into(),
            Self::InletTemperature(index) => format!("{{\"kind\":\"inlet-temperature\",\"index\":{index}}}"),
            Self::FanSpeedRatio => "{\"kind\":\"fan-speed-ratio\"}".into(),
            Self::SurfaceHtc(name) => format!("{{\"kind\":\"surface-htc\",\"surface\":{}}}", quote(name)),
            Self::ComponentPower(name) => format!("{{\"kind\":\"component-power\",\"component\":{}}}", quote(name)),
            Self::Material(target) => target.render(),
        }
    }
}

#[derive(Debug, Clone)]
enum Distribution {
    Gaussian { mean: f64, std_dev: f64 },
    Uniform { lo: f64, hi: f64 },
}
impl Distribution {
    fn parameter(&self, target: &Target) -> ParameterUncertainty {
        match *self {
            Self::Gaussian { mean, std_dev } => ParameterUncertainty::gaussian(target.name(), mean, std_dev, target.unit()),
            Self::Uniform { lo, hi } => ParameterUncertainty::uniform(target.name(), lo, hi, target.unit()),
        }
    }
    fn validate_support(&self, target: &Target) -> Result<()> {
        if let Self::Uniform { lo, hi } = *self {
            let lower_ok = if target.allows_zero() { lo >= 0.0 } else { lo > 0.0 };
            if !lower_ok || !hi.is_finite() || lo > hi {
                return Err(bad(format!("uniform support for {} leaves its admitted physical domain", target.name())));
            }
        }
        Ok(())
    }
    fn render(&self) -> Result<String> {
        match *self {
            Self::Gaussian { mean, std_dev } => Ok(format!("{{\"kind\":\"gaussian\",\"mean\":{},\"std_dev\":{}}}", number_json(mean)?, number_json(std_dev)?)),
            Self::Uniform { lo, hi } => Ok(format!("{{\"kind\":\"uniform\",\"lo\":{},\"hi\":{}}}", number_json(lo)?, number_json(hi)?)),
        }
    }
}

#[derive(Debug, Clone)]
struct ParameterSpec { target: Target, distribution: Distribution }

#[derive(Debug)]
pub(super) struct Config {
    pub(super) seed: u64,
    pub(super) samples: usize,
    pub(super) wall_seconds: f64,
    pub(super) threshold_k: Option<f64>,
    correlation: CorrelationModel,
    pub(super) correlation_label: &'static str,
    parameters: Vec<ParameterSpec>,
}

impl Config {
    pub(super) fn parse(text: &str, base: &J) -> Result<Self> {
        if text.len() as u64 > MAX_UQ_BYTES { return Err(bad("UQ request exceeds 4 MiB")); }
        let root = J::parse(text).map_err(|error| bad(error.to_string()))?;
        object(&root, &["schema", "seed", "samples", "wall_seconds", "temperature_limit_k", "correlation", "parameters"], "UQ request")?;
        if field(&root, "schema")?.as_str() != Some(SCHEMA) { return Err(bad("expected schema frankensim.cooling-network-uq.v1")); }
        let seed = string(field(&root, "seed")?, "seed")?.parse::<u64>().map_err(|_| bad("seed must be a decimal u64 string"))?;
        let samples = integer(field(&root, "samples")?, "samples", MAX_PRODUCT_SAMPLES)?;
        if samples < 2 { return Err(bad("samples must be at least two")); }
        let wall_seconds = positive(field(&root, "wall_seconds")?, "wall_seconds")?;
        if wall_seconds > 86_400.0 { return Err(bad("wall_seconds must be <= 86400")); }
        let threshold_k = root.get("temperature_limit_k").map(|value| positive(value, "temperature_limit_k")).transpose()?;
        let rows = array(field(&root, "parameters")?, "parameters", 256)?;
        if rows.is_empty() { return Err(bad("at least one uncertain parameter is required")); }
        validate_base(base)?;
        let mut names = BTreeSet::new();
        let mut parameters = Vec::with_capacity(rows.len());
        for row in rows {
            object(row, &["target", "distribution"], "parameter")?;
            let target = parse_target(field(row, "target")?)?;
            if !names.insert(target.name()) { return Err(bad(format!("duplicate uncertain target {}", target.name()))); }
            let distribution = parse_distribution(field(row, "distribution")?)?;
            distribution.validate_support(&target)?;
            validate_target(base, &target)?;
            parameters.push(ParameterSpec { target, distribution });
        }
        let (correlation, correlation_label) = parse_correlation(field(&root, "correlation")?, parameters.len())?;
        Ok(Self { seed, samples, wall_seconds, threshold_k, correlation, correlation_label, parameters })
    }

    pub(super) fn plan(&self) -> UqPlan {
        let mut plan = UqPlan::new("cooling-objective-temperature-k", PropagationMethod::MonteCarlo, self.samples).with_correlation(self.correlation.clone());
        plan.seed = self.seed;
        if let Some(threshold) = self.threshold_k { plan = plan.with_compliance_threshold(threshold); }
        for parameter in &self.parameters { plan = plan.with_parameter(parameter.distribution.parameter(&parameter.target)); }
        plan
    }

    pub(super) fn sample_request(&self, base: &J, values: &[f64]) -> Result<String> {
        if values.len() != self.parameters.len() { return Err(bad("sample parameter arity mismatch")); }
        let mut request = base.clone();
        for (parameter, &value) in self.parameters.iter().zip(values) {
            if !value.is_finite() { return Err(model_failure("sampled parameter is nonfinite")); }
            let positive_ok = if parameter.target.allows_zero() { value >= 0.0 } else { value > 0.0 };
            if !positive_ok { return Err(model_failure(format!("sampled {}={} leaves its admitted physical domain", parameter.target.name(), value))); }
            apply_target(&mut request, &parameter.target, value)?;
        }
        if self.parameters.iter().any(|parameter| matches!(parameter.target, Target::ComponentPower(_))) {
            recompute_component_total(&mut request)?;
        }
        render_json(&request)
    }

    pub(super) fn render_parameters(&self) -> Result<String> {
        self.parameters.iter().map(|parameter| Ok(format!(
            "{{\"target\":{},\"name\":{},\"unit\":{},\"distribution\":{}}}",
            parameter.target.render(), quote(&parameter.target.name()), quote(parameter.target.unit()), parameter.distribution.render()?
        ))).collect::<Result<Vec<_>>>().map(|rows| rows.join(","))
    }
}

fn validate_base(base: &J) -> Result<()> {
    if field(base, "schema")?.as_str() != Some("frankensim.cooling-network.v1") { return Err(bad("base request must use frankensim.cooling-network.v1")); }
    if base.get("transient").is_some() || base.get("design").is_some() || base.get("fan_speed_design").is_some() {
        return Err(bad("cooling-network-uq currently requires a steady base request without design controls"));
    }
    if field(field(base, "objective")?, "gradient")? != &J::Bool(false) { return Err(bad("set base objective.gradient=false for UQ sampling")); }
    Ok(())
}

fn validate_target(base: &J, target: &Target) -> Result<()> {
    match target {
        Target::AirDensity => { field(field(base, "air")?, "density_kg_m3")?; }
        Target::AirSpecificHeat => { field(field(base, "air")?, "specific_heat_j_kg_k")?; }
        Target::InletTemperature(index) => inlet_location(base, *index)?,
        Target::FanSpeedRatio => {
            let fan = field(base, "hydraulics")?.get("fan").ok_or_else(|| bad("fan-speed-ratio uncertainty requires hydraulics.fan"))?;
            field(fan, "speed_ratio")?;
        }
        Target::SurfaceHtc(name) => {
            let surfaces = array(field(field(base, "solid")?, "surfaces")?, "solid.surfaces", 4096)?;
            let surface = surfaces.iter().find(|surface| surface.str_field("name") == Some(name.as_str())).ok_or_else(|| bad(format!("unknown surface {name}")))?;
            if surface.get("convection").is_some() || surface.get("htc_w_m2_k").is_none() {
                return Err(bad(format!("surface {name} does not have a declared scalar h; do not override a correlation-derived coefficient")));
            }
        }
        Target::ComponentPower(name) => component_location(base, name)?,
        Target::Material(target) => target.validate(base)?,
    }
    Ok(())
}

fn parse_target(value: &J) -> Result<Target> {
    if let Some(target) = material::Target::parse(value)? { return Ok(Target::Material(target)); }
    object(value, &["kind", "index", "surface", "component"], "target")?;
    match field(value, "kind")?.as_str() {
        Some("air-density") => Ok(Target::AirDensity),
        Some("air-specific-heat") => Ok(Target::AirSpecificHeat),
        Some("inlet-temperature") => Ok(Target::InletTemperature(integer_raw(field(value, "index")?, "inlet index")?)),
        Some("fan-speed-ratio") => Ok(Target::FanSpeedRatio),
        Some("surface-htc") => Ok(Target::SurfaceHtc(string(field(value, "surface")?, "surface")?)),
        Some("component-power") => Ok(Target::ComponentPower(string(field(value, "component")?, "component")?)),
        _ => Err(bad("unknown UQ target kind")),
    }
}

fn parse_distribution(value: &J) -> Result<Distribution> {
    object(value, &["kind", "mean", "std_dev", "lo", "hi"], "distribution")?;
    match field(value, "kind")?.as_str() {
        Some("gaussian") => {
            let mean = number(field(value, "mean")?, "Gaussian mean")?;
            let std_dev = number(field(value, "std_dev")?, "Gaussian std_dev")?;
            if std_dev < 0.0 { return Err(bad("Gaussian std_dev must be nonnegative")); }
            Ok(Distribution::Gaussian { mean, std_dev })
        }
        Some("uniform") => {
            let lo = number(field(value, "lo")?, "uniform lo")?;
            let hi = number(field(value, "hi")?, "uniform hi")?;
            if lo > hi { return Err(bad("uniform lo must not exceed hi")); }
            Ok(Distribution::Uniform { lo, hi })
        }
        _ => Err(bad("distribution kind must be gaussian or uniform")),
    }
}

fn parse_correlation(value: &J, dimensions: usize) -> Result<(CorrelationModel, &'static str)> {
    object(value, &["kind", "matrix"], "correlation")?;
    match field(value, "kind")?.as_str() {
        Some("independent") => Ok((CorrelationModel::Independent, "independent")),
        Some("unknown") => Ok((CorrelationModel::Unknown, "unknown")),
        Some("joint-gaussian") => {
            let rows = array(field(value, "matrix")?, "correlation matrix", dimensions)?;
            if rows.len() != dimensions { return Err(bad("correlation matrix row count must match parameter count")); }
            let mut matrix = Vec::with_capacity(dimensions);
            for row in rows {
                let entries = array(row, "correlation row", dimensions)?;
                if entries.len() != dimensions { return Err(bad("correlation matrix must be square")); }
                matrix.push(entries.iter().map(|entry| number(entry, "correlation coefficient")).collect::<Result<Vec<_>>>()?);
            }
            Ok((CorrelationModel::JointGaussian { matrix }, "joint-gaussian"))
        }
        _ => Err(bad("correlation kind must be independent, unknown, or joint-gaussian")),
    }
}

fn apply_target(root: &mut J, target: &Target, value: f64) -> Result<()> {
    match target {
        Target::AirDensity => set_path_number(root, &["air", "density_kg_m3"], value),
        Target::AirSpecificHeat => set_path_number(root, &["air", "specific_heat_j_kg_k"], value),
        Target::InletTemperature(index) => set_inlet_temperature(root, *index, value),
        Target::FanSpeedRatio => set_path_number(root, &["hydraulics", "fan", "speed_ratio"], value),
        Target::SurfaceHtc(name) => {
            let surfaces = array_mut_path(root, &["solid", "surfaces"])?;
            let surface = surfaces.iter_mut().find(|surface| surface.str_field("name") == Some(name.as_str())).ok_or_else(|| bad("sample target surface disappeared"))?;
            set_member_number(surface, "htc_w_m2_k", value)
        }
        Target::ComponentPower(name) => {
            let components = array_mut_path(root, &["solid", "component_power", "components"])?;
            let component = components.iter_mut().find(|component| component.str_field("name") == Some(name.as_str())).ok_or_else(|| bad("sample target component disappeared"))?;
            set_member_number(component, "watts", value)
        }
        Target::Material(target) => target.apply(root, value),
    }
}

fn recompute_component_total(root: &mut J) -> Result<()> {
    let components = array_path(root, &["solid", "component_power", "components"])?;
    let mut total = 0.0_f64;
    for component in components {
        total += number(field(component, "watts")?, "component watts")?;
        if !total.is_finite() { return Err(model_failure("sampled component total power overflowed")); }
    }
    set_path_number(root, &["solid", "component_power", "total_w"], total)
}

fn inlet_location(base: &J, index: usize) -> Result<()> {
    let hydraulics = field(base, "hydraulics")?;
    if let Some(fan) = hydraulics.get("fan") {
        if index != 0 { return Err(bad("fan-driven cooling has exactly one external inlet temperature")); }
        field(fan, "temperature_k")?;
        return Ok(());
    }
    let boundaries = array(field(hydraulics, "boundaries")?, "hydraulic boundaries", 4096)?;
    let count = boundaries.iter().filter(|row| row.get("temperature_k").is_some()).count();
    if index >= count { Err(bad(format!("inlet temperature index {index} is out of range"))) } else { Ok(()) }
}

fn set_inlet_temperature(base: &mut J, index: usize, value: f64) -> Result<()> {
    let hydraulics = member_mut(base, "hydraulics")?;
    if let Some(fan) = member_opt_mut(hydraulics, "fan")? {
        if index != 0 { return Err(bad("fan-driven cooling has exactly one external inlet temperature")); }
        return set_member_number(fan, "temperature_k", value);
    }
    let boundaries = member_mut(hydraulics, "boundaries")?;
    let J::Array(rows) = boundaries else { return Err(bad("hydraulics.boundaries must be an array")); };
    let mut ordinal = 0usize;
    for row in rows {
        if row.get("temperature_k").is_some() {
            if ordinal == index { return set_member_number(row, "temperature_k", value); }
            ordinal += 1;
        }
    }
    Err(bad("sample inlet index disappeared"))
}

fn component_location(base: &J, name: &str) -> Result<()> {
    let components = array_path(base, &["solid", "component_power", "components"])?;
    if components.iter().any(|component| component.str_field("name") == Some(name)) { Ok(()) } else { Err(bad(format!("unknown component {name}"))) }
}

fn set_path_number(root: &mut J, path: &[&str], value: f64) -> Result<()> {
    let (last, parents) = path.split_last().ok_or_else(|| bad("empty JSON path"))?;
    let mut current = root;
    for key in parents { current = member_mut(current, key)?; }
    set_member_number(current, last, value)
}
fn set_member_number(value: &mut J, key: &str, number: f64) -> Result<()> {
    if !number.is_finite() { return Err(model_failure("sampled JSON number is nonfinite")); }
    let J::Object(members) = value else { return Err(bad(format!("expected object while setting {key}"))); };
    let slot = members.iter_mut().find_map(|(name, value)| (name == key).then_some(value)).ok_or_else(|| bad(format!("missing base field {key}")))?;
    *slot = J::Number { value: number, raw: number.to_string() };
    Ok(())
}
fn member_mut<'a>(value: &'a mut J, key: &str) -> Result<&'a mut J> {
    let J::Object(members) = value else { return Err(bad(format!("expected object while looking up {key}"))); };
    members.iter_mut().find_map(|(name, value)| (name == key).then_some(value)).ok_or_else(|| bad(format!("missing base field {key}")))
}
fn member_opt_mut<'a>(value: &'a mut J, key: &str) -> Result<Option<&'a mut J>> {
    let J::Object(members) = value else { return Err(bad(format!("expected object while looking up {key}"))); };
    Ok(members.iter_mut().find_map(|(name, value)| (name == key).then_some(value)))
}
fn array_path<'a>(value: &'a J, path: &[&str]) -> Result<&'a [J]> {
    let mut current = value;
    for key in path { current = field(current, key)?; }
    current.as_array().ok_or_else(|| bad("expected JSON array"))
}
fn array_mut_path<'a>(value: &'a mut J, path: &[&str]) -> Result<&'a mut Vec<J>> {
    let mut current = value;
    for key in path { current = member_mut(current, key)?; }
    match current { J::Array(items) => Ok(items), _ => Err(bad("expected JSON array")) }
}

fn render_json(value: &J) -> Result<String> {
    let mut output = String::new();
    write_json(value, &mut output)?;
    output.push('\n');
    Ok(output)
}
fn write_json(value: &J, output: &mut String) -> Result<()> {
    use std::fmt::Write as _;
    match value {
        J::Null => output.push_str("null"),
        J::Bool(flag) => output.push_str(if *flag { "true" } else { "false" }),
        J::Number { value, raw } => {
            if !value.is_finite() { return Err(model_failure("cannot serialize a nonfinite sample request")); }
            output.push_str(raw);
        }
        J::Str(text) => output.push_str(&quote(text)),
        J::Array(items) => {
            output.push('[');
            for (index, item) in items.iter().enumerate() { if index > 0 { output.push(','); } write_json(item, output)?; }
            output.push(']');
        }
        J::Object(members) => {
            output.push('{');
            for (index, (key, item)) in members.iter().enumerate() {
                if index > 0 { output.push(','); }
                let _ = write!(output, "{}:", quote(key));
                write_json(item, output)?;
            }
            output.push('}');
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/fan-correlated-hotspot.json"));

    #[test]
    fn sampling_mutates_only_targets_and_recomputes_power_total() {
        let base = J::parse(BASE).unwrap();
        let plan = r#"{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":4,"wall_seconds":60,"temperature_limit_k":303,"correlation":{"kind":"independent"},"parameters":[{"target":{"kind":"air-density"},"distribution":{"kind":"uniform","lo":1.18,"hi":1.22}},{"target":{"kind":"component-power","component":"chip"},"distribution":{"kind":"uniform","lo":0.9,"hi":1.1}}]}"#;
        let config = Config::parse(plan, &base).unwrap();
        let sampled = J::parse(&config.sample_request(&base, &[1.21, 1.07]).unwrap()).unwrap();
        assert_eq!(sampled.path(&["air", "density_kg_m3"]).and_then(J::as_f64), Some(1.21));
        assert_eq!(sampled.path(&["solid", "component_power", "total_w"]).and_then(J::as_f64), Some(1.07));
        assert_eq!(sampled.path(&["solid", "component_power", "components"]).unwrap().as_array().unwrap()[0].f64_field("watts"), Some(1.07));
    }

    #[test]
    fn dependence_admission_and_correlation_derived_h_are_fail_closed() {
        let base = J::parse(BASE).unwrap();
        let unknown = r#"{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":4,"wall_seconds":60,"correlation":{"kind":"unknown"},"parameters":[{"target":{"kind":"air-density"},"distribution":{"kind":"uniform","lo":1.18,"hi":1.22}},{"target":{"kind":"air-specific-heat"},"distribution":{"kind":"uniform","lo":1000,"hi":1010}}]}"#;
        let config = Config::parse(unknown, &base).unwrap();
        assert!(UqExecution::new(&config.plan()).is_err());
        let surface = unknown.replace("{\"kind\":\"air-density\"}", "{\"kind\":\"surface-htc\",\"surface\":\"first-face\"}");
        assert!(Config::parse(&surface, &base).is_err());
    }
}
