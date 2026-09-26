//! A probability study of an existing native cooling project, not another
//! physical-model grammar. Every sample changes only explicit project inputs;
//! geometry, material identities, solver policy and requirements stay intact.
//!
//! Version 1 admits fixed-count independent uniform inputs. An engineering
//! interval or card tolerance is never silently interpreted as a probability
//! law. The study inherits units, capabilities, physics seed, versions, memory
//! and per-solve budgets from its validated base project; its own sampling seed
//! and total wall/sample allowances are mandatory. Statistics remain Estimated.

use std::collections::{BTreeMap, BTreeSet};

use fs_ir::ast::{Node, NodeKind};
use fs_qty::{Dims, QtyAny};

use crate::{ConsequenceClass, DecisionGate, ProjectError, ProjectSpec,
    RequirementDirection, ThermalBoundaryCondition};

/// Independent schema version; existing project and optimization-study bytes
/// do not change merely because probability studies are now executable.
pub const VERSION: u32 = 1;
/// Bounded native-study source size, before parsing.
pub const MAX_SOURCE_BYTES: usize = 65_536;
/// Full native import/solve pipelines per study, not scalar callback evaluations.
pub const MAX_SAMPLES: usize = 256;

/// Exact project field addressed by a random parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Target {
    /// Declared dissipation watts; the original duty factor is retained.
    Power,
    /// A declared convection coefficient, never a correlation-derived one.
    ConvectionCoefficient,
    /// A declared convection reservoir temperature.
    ConvectionTemperature,
    /// All segment inlet declarations belonging to one named air branch.
    AirInletTemperature,
    /// Prescribed outward flux; negative values still mean inward heating.
    HeatFlux,
}

impl Target {
    /// Coherent parameter unit used by the sampler and report.
    #[must_use]
    pub const fn unit(self) -> &'static str {
        match self {
            Self::Power => "W",
            Self::ConvectionCoefficient => "W/m^2/K",
            Self::ConvectionTemperature | Self::AirInletTemperature => "K",
            Self::HeatFlux => "W/m^2",
        }
    }
    fn dims(self) -> Dims {
        use crate::spec::dims;
        match self {
            Self::Power => dims::POWER,
            Self::ConvectionCoefficient => dims::HEAT_TRANSFER_COEFFICIENT,
            Self::ConvectionTemperature | Self::AirInletTemperature => dims::TEMPERATURE,
            Self::HeatFlux => dims::HEAT_FLUX,
        }
    }
}

/// An explicit probability law in coherent units, in sampler declaration order.
#[derive(Debug, Clone, PartialEq)]
pub struct UniformParameter {
    /// Unique statistical parameter name.
    pub name: String,
    /// The native field to vary.
    pub target: Target,
    /// Region, boundary target, or air branch according to `target`.
    pub entity: String,
    /// Closed distribution support in coherent SI units; equality is deterministic.
    pub low: f64,
    /// Upper endpoint, not a confidence interval or certified bound on outputs.
    pub high: f64,
}

/// A source to pass through the ordinary native mesh-import quarantine.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MeshSource {
    /// Exact geometry-artifact role in the base project.
    pub role: String,
    /// Path relative to the uncertainty study, not the current working directory.
    pub path: String,
    /// Explicit source-coordinate unit.
    pub unit: String,
    /// Existing mesh promotion policy; no implicit repair allowance.
    pub max_hole_edges: usize,
}

/// Strictly parsed study. Private fields prevent edits that detach semantics
/// from the retained canonical declaration.
#[derive(Debug, Clone, PartialEq)]
pub struct UncertaintyStudy {
    canonical: String,
    project: String,
    samples: usize,
    seed: u64,
    wall_seconds: f64,
    qoi: String,
    geometry: Vec<MeshSource>,
    materials: Vec<String>,
    interfaces: Vec<String>,
    parameters: Vec<UniformParameter>,
}

/// A study bound to one admitted base project. Samples are fresh copies, so
/// parameters never accumulate and failed proposals never alter the base.
#[derive(Debug, Clone)]
pub struct BoundStudy {
    study: UncertaintyStudy,
    base: ProjectSpec,
}

type Result<T> = std::result::Result<T, ProjectError>;
fn error(detail: impl Into<String>) -> ProjectError {
    ProjectError { code: "project-uncertainty", detail: detail.into(),
        hint: "declare a version-1 fsim-uncertainty-study with explicit independent uniform inputs on the existing cooling project".into() }
}
fn list(node: &Node) -> Result<&[Node]> {
    match &node.kind { NodeKind::List(values) => Ok(values), _ => Err(error("expected a list")) }
}
fn fields<'a>(nodes: &'a [Node], keys: &[&str]) -> Result<BTreeMap<&'a str, &'a Node>> {
    if nodes.len() % 2 != 0 { return Err(error("every field requires a value")); }
    let mut result = BTreeMap::new();
    for pair in nodes.chunks_exact(2) {
        let NodeKind::Keyword(key) = &pair[0].kind else { return Err(error("expected a keyword")); };
        if !keys.contains(&key.as_str()) || result.insert(key.as_str(), &pair[1]).is_some() {
            return Err(error(format!("unknown or repeated field :{key}")));
        }
    }
    if result.len() != keys.len() { return Err(error(format!("required fields: {}", keys.join(", ")))); }
    Ok(result)
}
fn text(node: &Node) -> Result<String> {
    match &node.kind {
        NodeKind::Str(value) if !value.is_empty() && value.len() <= 1024
            && !value.chars().any(char::is_control) => Ok(value.clone()),
        _ => Err(error("expected a nonempty bounded string without control characters")),
    }
}
fn symbol(node: &Node, expected: &str) -> Result<()> {
    if matches!(&node.kind, NodeKind::Symbol(value) if value == expected) { Ok(()) }
    else { Err(error(format!("expected {expected}"))) }
}
fn integer(node: &Node) -> Result<u64> {
    match node.kind {
        NodeKind::Int(value) if value >= 0 => u64::try_from(value).map_err(|_| error("integer overflow")),
        _ => Err(error("expected a nonnegative exact integer")),
    }
}
fn quantity(node: &Node, dims: Dims) -> Result<f64> {
    match &node.kind {
        NodeKind::Qty { value, dims: found, .. } if *found == dims && value.is_finite() => Ok(*value),
        _ => Err(error(format!("expected an explicit finite {} quantity", dims.unit_string()))),
    }
}
fn paths(node: &Node) -> Result<Vec<String>> {
    let nodes = list(node)?;
    if nodes.len() > 32 { return Err(error("at most 32 paths per asset family")); }
    nodes.iter().map(text).collect()
}

impl UncertaintyStudy {
    /// Parse using the existing typed FrankenScript AST. Unknown/repeated fields,
    /// inferred distributions, implicit units, and undeclared independence refuse.
    ///
    /// # Errors
    /// Returns a project diagnostic before any file or numerical work.
    pub fn parse(source: &str) -> Result<Self> {
        if source.len() > MAX_SOURCE_BYTES { return Err(error("study source exceeds 64 KiB")); }
        let root = fs_ir::sexpr::parse(source).map_err(|e| error(e.to_string()))?;
        let nodes = list(&root)?;
        symbol(nodes.first().ok_or_else(|| error("empty study"))?, "fsim-uncertainty-study")?;
        let f = fields(&nodes[1..], &["version", "project", "samples", "seed", "wall-time",
            "method", "correlation", "qoi", "geometry", "materials", "interfaces", "parameters"])?;
        if integer(f["version"])? != u64::from(VERSION) { return Err(error("unsupported study version")); }
        symbol(f["method"], "monte-carlo")?;
        symbol(f["correlation"], "independent")?;
        let samples = usize::try_from(integer(f["samples"])?).map_err(|_| error("sample count overflow"))?;
        if !(2..=MAX_SAMPLES).contains(&samples) { return Err(error("samples must be in 2..=256")); }
        let wall_seconds = quantity(f["wall-time"], crate::spec::dims::TIME)?;
        if !(wall_seconds > 0.0 && wall_seconds <= 86_400.0) { return Err(error("wall-time must be in (0, 86400] seconds")); }
        let geometry_nodes = list(f["geometry"])?;
        if geometry_nodes.is_empty() || geometry_nodes.len() > 32 { return Err(error("declare 1..=32 mesh sources")); }
        let mut geometry = Vec::new();
        let mut roles = BTreeSet::new();
        for node in geometry_nodes {
            let n = list(node)?;
            symbol(n.first().ok_or_else(|| error("empty mesh source"))?, "mesh")?;
            let m = fields(&n[1..], &["role", "path", "unit", "max-hole-edges"])?;
            let role = text(m["role"])?;
            if !roles.insert(role.clone()) { return Err(error("duplicate geometry role")); }
            let max_hole_edges = usize::try_from(integer(m["max-hole-edges"])?).map_err(|_| error("repair count overflow"))?;
            if max_hole_edges > 10_000 { return Err(error("repair edge budget exceeds 10000")); }
            geometry.push(MeshSource { role, path: text(m["path"])?, unit: text(m["unit"])?, max_hole_edges });
        }
        let parameter_nodes = list(f["parameters"])?;
        if parameter_nodes.is_empty() || parameter_nodes.len() > 32 { return Err(error("declare 1..=32 parameters")); }
        let mut parameters = Vec::new();
        let mut names = BTreeSet::new();
        let mut targets = BTreeSet::new();
        for node in parameter_nodes {
            let n = list(node)?;
            symbol(n.first().ok_or_else(|| error("empty probability law"))?, "uniform")?;
            let p = fields(&n[1..], &["name", "target", "entity", "low", "high"])?;
            let target = match &p["target"].kind {
                NodeKind::Symbol(s) => match s.as_str() {
                    "power" => Target::Power,
                    "convection-coefficient" => Target::ConvectionCoefficient,
                    "convection-temperature" => Target::ConvectionTemperature,
                    "air-inlet-temperature" => Target::AirInletTemperature,
                    "heat-flux" => Target::HeatFlux,
                    _ => return Err(error("unsupported random project field")),
                },
                _ => return Err(error("parameter target must be a symbol")),
            };
            let name = text(p["name"])?;
            let entity = text(p["entity"])?;
            if !names.insert(name.clone()) || !targets.insert((target, entity.clone())) {
                return Err(error("duplicate parameter name or physical target"));
            }
            let low = quantity(p["low"], target.dims())?;
            let high = quantity(p["high"], target.dims())?;
            if low > high || (target == Target::Power && low < 0.0)
                || (matches!(target, Target::ConvectionCoefficient | Target::ConvectionTemperature
                    | Target::AirInletTemperature) && low <= 0.0) {
                return Err(error("invalid probability support for the physical target"));
            }
            parameters.push(UniformParameter { name, target, entity, low, high });
        }
        let qoi = text(f["qoi"])?;
        if qoi != "temperature-max" { return Err(error("this native lane requires temperature-max")); }
        Ok(Self { canonical: fs_ir::sexpr::print(&root).map_err(|e| error(e.to_string()))?,
            project: text(f["project"])?, samples, seed: integer(f["seed"])?, wall_seconds, qoi,
            geometry, materials: paths(f["materials"])?, interfaces: paths(f["interfaces"])?, parameters })
    }
    /// Canonical source, including all explicit path and parameter declarations.
    #[must_use] pub fn canonical(&self) -> &str { &self.canonical }
    /// Referenced native project path.
    #[must_use] pub fn project_path(&self) -> &str { &self.project }
    /// Original fixed sample budget.
    #[must_use] pub const fn samples(&self) -> usize { self.samples }
    /// Statistical sampler seed, separate from the base project's physics seed.
    #[must_use] pub const fn seed(&self) -> u64 { self.seed }
    /// Shared numerical-work wall allowance.
    #[must_use] pub const fn wall_seconds(&self) -> f64 { self.wall_seconds }
    /// Exact native observable.
    #[must_use] pub fn qoi(&self) -> &str { &self.qoi }
    /// Ordered probability laws.
    #[must_use] pub fn parameters(&self) -> &[UniformParameter] { &self.parameters }
    /// Geometry sources, matched by role, not by path order.
    #[must_use] pub fn geometry(&self) -> &[MeshSource] { &self.geometry }
    /// Material pack sources.
    #[must_use] pub fn materials(&self) -> &[String] { &self.materials }
    /// Interface pack sources.
    #[must_use] pub fn interfaces(&self) -> &[String] { &self.interfaces }

    /// Bind semantic targets and check endpoint project admission before any
    /// model evaluation. Each actual sample is independently revalidated too.
    ///
    /// # Errors
    /// Refuses invalid bases, non-advisory decisions, ambiguous/missing fields,
    /// geometry role mismatches, unsupported QoIs, or support outside the base's
    /// declared envelope. This never expands a material or operating domain.
    pub fn bind(self, base: &ProjectSpec) -> Result<BoundStudy> {
        validate_project(base)?;
        let metadata = base.metadata.as_ref().ok_or_else(|| error("missing metadata"))?;
        if metadata.decision_gate != DecisionGate::ScopingEstimate || metadata.consequence != ConsequenceClass::Advisory {
            return Err(error("empirical native UQ is advisory scoping, not a compliance signoff"));
        }
        let requirements = base.requirements.as_deref().unwrap_or(&[]);
        if requirements.len() != 1 || requirements[0].qoi != self.qoi
            || requirements[0].direction != RequirementDirection::AtMost {
            return Err(error("declare one matching at-most temperature-max requirement"));
        }
        let envelope = base.envelope.as_ref().ok_or_else(|| error("missing operating envelope"))?;
        for parameter in &self.parameters {
            if matches!(parameter.target, Target::ConvectionTemperature | Target::AirInletTemperature)
                && (parameter.low < envelope.ambient_lo.value || parameter.high > envelope.ambient_hi.value) {
                return Err(error("temperature support exceeds the unchanged operating envelope"));
            }
        }
        let threshold = requirements[0].limit.value - requirements[0].margin.value;
        if !threshold.is_finite() { return Err(error("effective compliance threshold is not finite")); }
        let geometry = base.geometry.as_deref().unwrap_or(&[]);
        if geometry.len() != self.geometry.len() || geometry.iter().any(|g|
            !self.geometry.iter().any(|source| source.role == g.role)
            || !matches!(g.format.as_str(), "stl" | "obj" | "ply")) {
            return Err(error("every native mesh artifact needs exactly one source; STEP is not admitted in this lane"));
        }
        let bound = BoundStudy { study: self, base: base.clone() };
        for high in [false, true] {
            let values = bound.study.parameters.iter().map(|p| if high { p.high } else { p.low }).collect::<Vec<_>>();
            bound.sample_project(&values)?;
        }
        Ok(bound)
    }
}

impl BoundStudy {
    /// Admitted immutable probability declaration.
    #[must_use] pub const fn study(&self) -> &UncertaintyStudy { &self.study }
    /// Unmodified native base.
    #[must_use] pub const fn base(&self) -> &ProjectSpec { &self.base }
    /// Numerical pass threshold with the original explicit margin applied.
    /// This is not the native engineering uncertainty verdict.
    #[must_use] pub fn threshold_k(&self) -> f64 {
        let requirement = &self.base.requirements.as_ref().expect("bound requirement")[0];
        requirement.limit.value - requirement.margin.value
    }
    /// Produce a fresh validated native sample without changing any unrelated
    /// physical field or overwriting immutable material-card data.
    ///
    /// # Errors
    /// Refuses incorrect arity, out-of-support/nonfinite values, missing target
    /// fields and native project violations. No clipping or sample skipping.
    pub fn sample_project(&self, values: &[f64]) -> Result<ProjectSpec> {
        if values.len() != self.study.parameters.len() { return Err(error("sample arity differs")); }
        let mut project = self.base.clone();
        for (parameter, &value) in self.study.parameters.iter().zip(values) {
            if !value.is_finite() || value < parameter.low || value > parameter.high { return Err(error("sample outside declared probability support")); }
            apply(&mut project, parameter, value)?;
        }
        validate_project(&project)?;
        Ok(project)
    }
}
fn validate_project(project: &ProjectSpec) -> Result<()> {
    if let Some(finding) = project.validate().first() { Err(error(format!("{}: {}", finding.code, finding.what))) }
    else { Ok(()) }
}
fn apply(project: &mut ProjectSpec, parameter: &UniformParameter, value: f64) -> Result<()> {
    let mut matches = 0;
    if parameter.target == Target::Power {
        for row in project.power.as_mut().ok_or_else(|| error("missing power map"))? {
            if row.region == parameter.entity { row.watts = QtyAny::new(value, parameter.target.dims()); matches += 1; }
        }
    } else {
        let setup = project.cooling.as_mut().and_then(|c| c.conduction.as_mut())
            .ok_or_else(|| error("native conduction setup is required"))?;
        for row in &mut setup.boundaries {
            if parameter.target == Target::AirInletTemperature {
                if let ThermalBoundaryCondition::AirflowConvection { branch, inlet_temperature, .. } = &mut row.condition {
                    if branch == &parameter.entity { inlet_temperature.value = value; matches += 1; }
                }
            } else if row.target == parameter.entity {
                match (&mut row.condition, parameter.target) {
                    (ThermalBoundaryCondition::Convection { coefficient, .. }, Target::ConvectionCoefficient) => coefficient.value = value,
                    (ThermalBoundaryCondition::Convection { reference_temperature, .. }, Target::ConvectionTemperature) => reference_temperature.value = value,
                    (ThermalBoundaryCondition::HeatFlux { outward_flux }, Target::HeatFlux) => outward_flux.value = value,
                    _ => return Err(error("random field does not match the declared boundary law; derived coefficients cannot be overwritten")),
                }
                matches += 1;
            }
        }
    }
    if matches == 0 || (parameter.target != Target::AirInletTemperature && matches != 1) {
        return Err(error(format!("random parameter {} has a missing or ambiguous target {}", parameter.name, parameter.entity)));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
