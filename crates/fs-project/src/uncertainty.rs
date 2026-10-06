//! A probability study of an existing native cooling project, not another
//! physical-model grammar. Every sample changes only explicit project inputs;
//! geometry, material identities, solver policy and requirements stay intact.
//!
//! Uniform inputs require explicit independence or a declared Gaussian copula.
//! Version 1 admits fixed-count propagation with Monte Carlo
//! or explicitly replicated randomized Sobol quadrature. Version 2 requires
//! an explicit Bernoulli-mixture policy for sequential Monte Carlo decisions.
//! Version 3 adds fixed-count mean controls from support secants or one nominal adjoint.
//! An engineering
//! interval or card tolerance is never silently interpreted as a probability
//! law. The study inherits units, capabilities, physics seed, versions, memory
//! and per-solve budgets from its validated base project; its own sampling seed
//! and total wall/sample allowances are mandatory. Statistics remain Estimated.

use std::collections::{BTreeMap, BTreeSet};

use fs_ir::ast::{Node, NodeKind};
use fs_qty::{Dims, QtyAny};

use crate::{ConsequenceClass, DecisionGate, ProjectError, ProjectSpec,
    RequirementDirection, ThermalBoundaryCondition};

pub mod mean_control;
pub use mean_control::MeanControlPolicy;

/// Latest native probability-study schema. Versions 1 and 2 preserve their
/// original semantics; version 3 requires an explicit mean-control solve cap.
pub const VERSION: u32 = 3;
/// Bounded native-study source size, before parsing.
pub const MAX_SOURCE_BYTES: usize = 65_536;
/// Full native import/solve sample pipelines, excluding explicitly capped probes.
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
    /// Absolute speed ratio of one explicitly named native fan-system bank.
    /// The source curve and its admitted speed domain remain unchanged.
    FanSpeedRatio,
    /// Ambient temperature of one declared natural-convection boundary.
    NaturalConvectionAmbient,
    /// Reservoir temperature of one named radiating surface, not its wall.
    RadiationReservoirTemperature,
    /// Absolute prescribed temperature of one uniform Dirichlet boundary.
    /// This is a solid-boundary value, not an ambient-fluid temperature.
    FixedTemperature,
}

impl Target {
    /// Coherent parameter unit used by the sampler and report.
    #[must_use]
    pub const fn unit(self) -> &'static str {
        match self {
            Self::Power => "W",
            Self::ConvectionCoefficient => "W/m^2/K",
            Self::ConvectionTemperature | Self::AirInletTemperature
            | Self::NaturalConvectionAmbient | Self::RadiationReservoirTemperature
            | Self::FixedTemperature => "K",
            Self::HeatFlux => "W/m^2",
            Self::FanSpeedRatio => "1",
        }
    }
    fn dims(self) -> Dims {
        use crate::spec::dims;
        match self {
            Self::Power => dims::POWER,
            Self::ConvectionCoefficient => dims::HEAT_TRANSFER_COEFFICIENT,
            Self::ConvectionTemperature | Self::AirInletTemperature
            | Self::NaturalConvectionAmbient | Self::RadiationReservoirTemperature
            | Self::FixedTemperature => dims::TEMPERATURE,
            Self::HeatFlux => dims::HEAT_FLUX,
            Self::FanSpeedRatio => Dims::NONE,
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
    /// Region, boundary target, air branch, fan bank, or radiating-surface name
    /// according to `target`; a radiating surface is not its geometric target.
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

/// Predeclared stopping policy for the probability of the native numerical
/// event `temperature-max <= requirement limit - margin`. This does not
/// change the project's engineering verdict or confer physical validation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CompliancePolicy {
    /// Probability target strictly inside (0, 1).
    pub required_probability: f64,
    /// Error level strictly inside (0, 1), fixed before any samples.
    pub alpha: f64,
    /// Earliest decision ordinal, between two and the original sample cap.
    pub min_samples: usize,
}

/// Fixed randomized-quadrature layout, declared before any native solves.
/// Replicates use independent Owen scramble keys; points inside a net are
/// dependent and cannot enter a Bernoulli-iid confidence sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QmcLayout {
    /// Independently scrambled nets, at least two.
    pub replicates: usize,
    /// Power-of-two point count per net, at least two.
    pub samples_per_replicate: usize,
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
    latent_correlation: Option<Vec<Vec<f64>>>,
    compliance: Option<CompliancePolicy>,
    qmc: Option<QmcLayout>,
    mean_control: Option<MeanControlPolicy>,
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
        hint: "declare uniform inputs with explicit dependence; version 1 is fixed-count, version 2 requires Bernoulli compliance, version 3 requires coordinate-secant or nominal-adjoint mean controls with max-solves".into() }
}
fn list(node: &Node) -> Result<&[Node]> {
    match &node.kind { NodeKind::List(values) => Ok(values), _ => Err(error("expected a list")) }
}
fn fields<'a>(nodes: &'a [Node], keys: &[&str]) -> Result<BTreeMap<&'a str, &'a Node>> {
    fields_with_optional(nodes, keys, &[])
}
fn fields_with_optional<'a>(nodes: &'a [Node], keys: &[&str], optional: &[&str]) -> Result<BTreeMap<&'a str, &'a Node>> {
    if nodes.len() % 2 != 0 { return Err(error("every field requires a value")); }
    let mut result = BTreeMap::new();
    for pair in nodes.chunks_exact(2) {
        let NodeKind::Keyword(key) = &pair[0].kind else { return Err(error("expected a keyword")); };
        if (!keys.contains(&key.as_str()) && !optional.contains(&key.as_str()))
            || result.insert(key.as_str(), &pair[1]).is_some() {
            return Err(error(format!("unknown or repeated field :{key}")));
        }
    }
    if keys.iter().any(|key| !result.contains_key(key)) { return Err(error(format!("required fields: {}", keys.join(", ")))); }
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
        NodeKind::Float(value) if dims == Dims::NONE && value.is_finite() => Ok(*value),
        NodeKind::Int(value) if dims == Dims::NONE => Ok(*value as f64),
        NodeKind::Qty { value, dims: found, .. } if *found == dims && value.is_finite() => Ok(*value),
        _ => Err(error(format!("expected an explicit finite {} quantity", dims.unit_string()))),
    }
}
fn probability(node: &Node) -> Result<f64> {
    match node.kind {
        NodeKind::Float(value) if value.is_finite() && value > 0.0 && value < 1.0 => Ok(value),
        _ => Err(error("probability and alpha require dimensionless numbers strictly inside (0, 1)")),
    }
}
fn compliance_policy(node: &Node, samples: usize) -> Result<CompliancePolicy> {
    let nodes = list(node)?;
    symbol(nodes.first().ok_or_else(|| error("empty compliance policy"))?, "bernoulli-mixture")?;
    let f = fields(&nodes[1..], &["required-probability", "alpha", "min-samples"])?;
    let required_probability = probability(f["required-probability"])?;
    let alpha = probability(f["alpha"])?;
    if !alpha.recip().is_finite() { return Err(error("compliance alpha reciprocal exceeds the finite range")); }
    let min_samples = usize::try_from(integer(f["min-samples"])?).map_err(|_| error("minimum sample count overflow"))?;
    if !(2..=samples).contains(&min_samples) { return Err(error("compliance min-samples must be in 2..=samples")); }
    Ok(CompliancePolicy { required_probability, alpha, min_samples })
}
fn qmc_layout(node: &Node, samples: usize) -> Result<QmcLayout> {
    let nodes = list(node)?;
    symbol(nodes.first().ok_or_else(|| error("empty QMC layout"))?, "owen-scrambled-sobol")?;
    let f = fields(&nodes[1..], &["replicates", "samples-per-replicate"])?;
    let replicates = usize::try_from(integer(f["replicates"])?)
        .map_err(|_| error("QMC replicate count overflow"))?;
    let samples_per_replicate = usize::try_from(integer(f["samples-per-replicate"])?)
        .map_err(|_| error("QMC point count overflow"))?;
    if !(2..=256).contains(&replicates) || samples_per_replicate < 2
        || !samples_per_replicate.is_power_of_two()
        || replicates.checked_mul(samples_per_replicate) != Some(samples) {
        return Err(error("QMC requires at least two replicates times a power-of-two point count >=2, exactly matching :samples"));
    }
    Ok(QmcLayout { replicates, samples_per_replicate })
}
fn paths(node: &Node) -> Result<Vec<String>> {
    let nodes = list(node)?;
    if nodes.len() > 32 { return Err(error("at most 32 paths per asset family")); }
    nodes.iter().map(text).collect()
}

fn latent_correlation(node: &Node, dimension: usize) -> Result<Option<Vec<Vec<f64>>>> {
    if matches!(&node.kind, NodeKind::Symbol(value) if value == "independent") {
        return Ok(None);
    }
    let nodes = list(node)?;
    symbol(nodes.first().ok_or_else(|| error("empty dependence declaration"))?, "gaussian-copula")?;
    let f = fields(&nodes[1..], &["latent-correlation"])?;
    let rows = list(f["latent-correlation"])?;
    if rows.len() != dimension {
        return Err(error("latent correlation matrix must match parameter declaration order, including constant marginals"));
    }
    let matrix = rows.iter().map(|row| {
        let entries = list(row)?;
        if entries.len() != dimension { return Err(error("latent correlation matrix must be square")); }
        entries.iter().map(|entry| match entry.kind {
            NodeKind::Int(value) => Ok(value as f64),
            NodeKind::Float(value) if value.is_finite() => Ok(value),
            _ => Err(error("latent correlations must be finite dimensionless numbers")),
        }).collect::<Result<Vec<_>>>()
    }).collect::<Result<Vec<_>>>()?;
    // Numerical symmetry, range, unit-diagonal and PSD admission belong to
    // fs-uq's copula executor, before any native ledger or physical solve.
    Ok(Some(matrix))
}

impl UncertaintyStudy {
    /// Parse using the existing typed FrankenScript AST. Unknown/repeated fields,
    /// inferred distributions, implicit units, and undeclared dependence refuse.
    ///
    /// # Errors
    /// Returns a project diagnostic before any file or numerical work.
    pub fn parse(source: &str) -> Result<Self> {
        if source.len() > MAX_SOURCE_BYTES { return Err(error("study source exceeds 64 KiB")); }
        let root = fs_ir::sexpr::parse(source).map_err(|e| error(e.to_string()))?;
        let nodes = list(&root)?;
        symbol(nodes.first().ok_or_else(|| error("empty study"))?, "fsim-uncertainty-study")?;
        let f = fields_with_optional(&nodes[1..], &["version", "project", "samples", "seed", "wall-time",
            "method", "correlation", "qoi", "geometry", "materials", "interfaces", "parameters"], &["compliance", "qmc", "mean-control"])?;
        let version = integer(f["version"])?;
        if !(1..=u64::from(VERSION)).contains(&version) { return Err(error("unsupported study version")); }
        if (version == 2) != f.contains_key("compliance") {
            return Err(error("only version 2 requires and admits a compliance policy"));
        }
        if (version == 3) != f.contains_key("mean-control") {
            return Err(error("only version 3 requires and admits a mean-control policy"));
        }
        let randomized_qmc = match &f["method"].kind {
            NodeKind::Symbol(method) if method == "monte-carlo" => false,
            NodeKind::Symbol(method) if method == "quasi-monte-carlo" => true,
            _ => return Err(error("method must be monte-carlo or quasi-monte-carlo")),
        };
        if randomized_qmc != f.contains_key("qmc") {
            return Err(error("quasi-monte-carlo requires an explicit :qmc layout; monte-carlo forbids it"));
        }
        if randomized_qmc && version == 2 {
            return Err(error("dependent QMC points cannot use the version 2 Bernoulli-iid compliance policy"));
        }
        let samples = usize::try_from(integer(f["samples"])?).map_err(|_| error("sample count overflow"))?;
        if !(2..=samples).contains(&samples) || samples > MAX_SAMPLES { return Err(error("samples must be in 2..=256")); }
        let compliance = f.get("compliance").map(|node| compliance_policy(node, samples)).transpose()?;
        let qmc = f.get("qmc").map(|node| qmc_layout(node, samples)).transpose()?;
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
        if qmc.is_some() && parameter_nodes.len() > 10 {
            return Err(error("replicated Sobol QMC supports at most 10 parameters; no Monte Carlo tail fallback"));
        }
        let latent_correlation = latent_correlation(f["correlation"], parameter_nodes.len())?;
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
                    "fan-speed-ratio" => Target::FanSpeedRatio,
                    "natural-convection-ambient" => Target::NaturalConvectionAmbient,
                    "radiation-reservoir-temperature" => Target::RadiationReservoirTemperature,
                    "fixed-temperature" => Target::FixedTemperature,
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
                    | Target::AirInletTemperature | Target::FanSpeedRatio
                    | Target::NaturalConvectionAmbient | Target::RadiationReservoirTemperature
                    | Target::FixedTemperature) && low <= 0.0) {
                return Err(error("invalid probability support for the physical target"));
            }
            parameters.push(UniformParameter { name, target, entity, low, high });
        }
        let mean_control = f.get("mean-control").map(|node| mean_control::parse(node, &parameters)).transpose()?;
        let qoi = text(f["qoi"])?;
        if qoi != "temperature-max" { return Err(error("this native lane requires temperature-max")); }
        Ok(Self { canonical: fs_ir::sexpr::print(&root).map_err(|e| error(e.to_string()))?,
            project: text(f["project"])?, samples, seed: integer(f["seed"])?, wall_seconds, qoi,
            geometry, materials: paths(f["materials"])?, interfaces: paths(f["interfaces"])?, parameters,
            latent_correlation, compliance, qmc, mean_control })
    }
    /// Canonical source, including all explicit path and parameter declarations.
    #[must_use] pub fn canonical(&self) -> &str { &self.canonical }
    /// Referenced native project path.
    #[must_use] pub fn project_path(&self) -> &str { &self.project }
    /// Original lifetime sample budget, including sequential decision studies.
    #[must_use] pub const fn samples(&self) -> usize { self.samples }
    /// Statistical sampler seed, separate from the base project's physics seed.
    #[must_use] pub const fn seed(&self) -> u64 { self.seed }
    /// Shared numerical-work wall allowance.
    #[must_use] pub const fn wall_seconds(&self) -> f64 { self.wall_seconds }
    /// Exact native observable.
    #[must_use] pub fn qoi(&self) -> &str { &self.qoi }
    /// Ordered probability laws.
    #[must_use] pub fn parameters(&self) -> &[UniformParameter] { &self.parameters }
    /// Gaussian-copula correlation of latent standard normals, in parameter
    /// declaration order; not Pearson correlation of physical uniform inputs.
    /// Absence means the source explicitly declared independent marginals.
    /// The statistical executor owns numerical matrix/PSD admission.
    #[must_use] pub fn latent_correlation(&self) -> Option<&[Vec<f64>]> {
        self.latent_correlation.as_deref()
    }
    /// Predeclared Bernoulli stopping policy; absent for fixed-count studies.
    #[must_use] pub fn compliance(&self) -> Option<&CompliancePolicy> { self.compliance.as_ref() }
    /// Explicit randomized Sobol layout; absent for the Monte Carlo method.
    #[must_use] pub const fn qmc(&self) -> Option<QmcLayout> { self.qmc }
    /// Explicit whole-model mean calibration; absent in versions 1 and 2.
    #[must_use] pub const fn mean_control(&self) -> Option<MeanControlPolicy> { self.mean_control }
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
            if matches!(parameter.target, Target::ConvectionTemperature | Target::AirInletTemperature
                | Target::NaturalConvectionAmbient)
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
        if let Some(policy) = bound.study.mean_control {
            for ordinal in 0..policy.probe_count(&bound.study.parameters) {
                bound.sample_project(&policy.probe(&bound.study.parameters, ordinal)?)?;
            }
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
    } else if parameter.target == Target::FanSpeedRatio {
        let system = project.cooling.as_mut().and_then(|c| c.fan_system.as_mut())
            .ok_or_else(|| error("fan-speed-ratio requires a declared native fan system"))?;
        // ProjectSpec's structural validation does not replace this family's
        // admission. Check identities, topology, curves and domains before a
        // random input can select a bank or repair an invalid base declaration.
        system.validate().map_err(|e| error(format!("{}: {}", e.code, e.detail)))?;
        for bank in &mut system.banks {
            if bank.bank_id == parameter.entity {
                let (low, high) = bank.speed_ratio_domain;
                if value < low || value > high {
                    return Err(error(format!(
                        "fan speed support for {} exceeds the unchanged declared domain [{low}, {high}]",
                        parameter.entity
                    )));
                }
                // Absolute ratio relative to the retained source curve, not a
                // multiplier of the base operating speed or a previous sample.
                // The ordinary native flow producer owns fan affinity laws.
                bank.speed_ratio = value;
                matches += 1;
            }
        }
    } else if parameter.target == Target::RadiationReservoirTemperature {
        // A radiative reservoir is independent of the ambient-fluid envelope.
        // Change only its explicit temperature, never emissivity, card query,
        // convection ambient, surface ownership, or a derived effective Robin row.
        let radiation = project.cooling.as_mut().and_then(|c| c.conduction.as_mut())
            .and_then(|c| c.radiation.as_mut()).ok_or_else(|| error("declared radiation is required"))?;
        for surface in &mut radiation.surfaces {
            if surface.name == parameter.entity {
                surface.reservoir_temperature.value = value;
                matches += 1;
            }
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
                    // Preserve the original Dirichlet law and every unrelated
                    // value. The solid/material solver admits each actual field;
                    // the fluid ambient envelope does not bound a fixed wall.
                    (ThermalBoundaryCondition::FixedTemperature { temperature }, Target::FixedTemperature) => temperature.value = value,
                    (ThermalBoundaryCondition::NaturalConvection { ambient_temperature, .. }, Target::NaturalConvectionAmbient) => ambient_temperature.value = value,
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
#[cfg(test)]
mod boundary_temperature_tests;

#[cfg(test)]
mod fixed_temperature_tests;
