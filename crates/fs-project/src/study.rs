//! The versioned `.fsim` study schema (bead `frankensim-rc-root-q61wp.20`):
//! user-facing contract for topology optimization and design study workflows.
//!
//! A study defines the optimization problem over the Five Explicits:
//! design domain (SDF plate with holes or grid), physics (CutFEM elasticity),
//! boundary conditions and loads, compliance objective, volume fraction constraint,
//! and optimizer settings.
//!
//! Every omission is a named [`Violation`], undeclared units refuse, volume
//! fractions outside (0, 1) refuse, loads on non-boundary regions refuse, and
//! wrong objective dimensions refuse at admission.

use std::fmt::Write as _;

use fs_blake3::{ContentHash, hash_domain};
use fs_ir::ast::{Node, NodeKind};
use fs_qty::{Dims, QtyAny};
use fs_scenario::Violation;

use crate::spec::{ConsequenceClass, DecisionGate, Metadata, Seeds, UnitsDoctrine, Versions};

/// Current `.fsim` study schema version.
pub const STUDY_FSIM_VERSION: u32 = 1;

/// Domain for canonical `.fsim` study hashing.
pub const STUDY_CANONICAL_DOMAIN: &str = "org.frankensim.fs-project.study.canonical.v1";

/// Hash canonical `.fsim` study bytes under the schema's domain.
#[must_use]
pub fn canonical_study_hash(bytes: &[u8]) -> ContentHash {
    hash_domain(STUDY_CANONICAL_DOMAIN, bytes)
}

/// The complete study specification.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct StudySpec {
    /// Study metadata.
    pub metadata: Option<Metadata>,
    /// Five Explicits: versions.
    pub versions: Option<Versions>,
    /// Five Explicits: seeds.
    pub seeds: Option<Seeds>,
    /// Five Explicits: budgets.
    pub budgets: Option<StudyBudgets>,
    /// Five Explicits: capabilities.
    pub capabilities: Option<Vec<String>>,
    /// Five Explicits: units doctrine.
    pub units: Option<UnitsDoctrine>,
    /// Design domain declaration.
    pub domain: Option<StudyDomain>,
    /// Physics model settings.
    pub physics: Option<StudyPhysics>,
    /// Boundary conditions and loads.
    pub scenario: Option<StudyScenario>,
    /// Optimization objective.
    pub objective: Option<StudyObjective>,
    /// Design constraints.
    pub constraints: Option<StudyConstraints>,
    /// Optimizer algorithm and hyperparameters.
    pub optimizer: Option<StudyOptimizer>,
}

/// Study budgets.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct StudyBudgets {
    /// Wall time budget.
    pub wall_time: Option<QtyAny>,
    /// Maximum resident memory.
    pub memory_bytes: Option<u64>,
    /// Maximum optimization iterations.
    pub max_iterations: Option<usize>,
}

/// Design domain declaration.
#[derive(Debug, Clone, PartialEq)]
pub struct StudyDomain {
    /// Domain type (e.g. "sdf-plate-with-holes").
    pub domain_type: String,
    /// Bounding box: ([x_min, y_min], [x_max, y_max]).
    pub bounds: ([f64; 2], [f64; 2]),
    /// Initial circular holes.
    pub initial_holes: Vec<StudyHole>,
}

impl Default for StudyDomain {
    fn default() -> Self {
        Self {
            domain_type: "sdf-plate-with-holes".to_string(),
            bounds: ([0.0, 0.0], [1.0, 1.0]),
            initial_holes: Vec::new(),
        }
    }
}

/// Circular hole definition in the plate.
#[derive(Debug, Clone, PartialEq)]
pub struct StudyHole {
    /// Hole center [x, y].
    pub center: [f64; 2],
    /// Hole radius.
    pub radius: f64,
}

/// Physics configuration.
#[derive(Debug, Clone, PartialEq)]
pub struct StudyPhysics {
    /// Physics type (e.g. "elasticity-2d").
    pub physics_type: String,
    /// Background quadtree mesh refinement level.
    pub mesh_level: u32,
    /// Young's modulus with units (default 1.0 Pa).
    pub youngs_modulus: Option<QtyAny>,
    /// Poisson's ratio (dimensionless).
    pub poissons_ratio: Option<f64>,
}

impl Default for StudyPhysics {
    fn default() -> Self {
        Self {
            physics_type: "elasticity-2d".to_string(),
            mesh_level: 4,
            youngs_modulus: None,
            poissons_ratio: Some(0.3),
        }
    }
}

/// Scenario boundary conditions and applied loads.
#[derive(Debug, Clone, PartialEq)]
pub struct StudyScenario {
    /// Boundary segment held fixed (Dirichlet), e.g. "left".
    pub fixed_boundary: String,
    /// Region or boundary where load is applied, e.g. "right".
    pub load_region: String,
    /// Applied force magnitude.
    pub load_force: Option<QtyAny>,
}

impl Default for StudyScenario {
    fn default() -> Self {
        Self {
            fixed_boundary: "left".to_string(),
            load_region: "right".to_string(),
            load_force: None,
        }
    }
}

/// Optimization objective.
#[derive(Debug, Clone, PartialEq)]
pub struct StudyObjective {
    /// Type of objective, e.g. "compliance".
    pub objective_type: String,
    /// Optimization sense, e.g. "minimize".
    pub sense: String,
    /// Unit expression declared for the objective, e.g. "J".
    pub unit: String,
}

impl Default for StudyObjective {
    fn default() -> Self {
        Self {
            objective_type: "compliance".to_string(),
            sense: "minimize".to_string(),
            unit: "J".to_string(),
        }
    }
}

/// Design constraints.
#[derive(Debug, Clone, PartialEq)]
pub struct StudyConstraints {
    /// Material volume/area fraction target in (0, 1).
    pub volume_fraction: f64,
}

impl Default for StudyConstraints {
    fn default() -> Self {
        Self {
            volume_fraction: 0.5,
        }
    }
}

/// Optimizer settings.
#[derive(Debug, Clone, PartialEq)]
pub struct StudyOptimizer {
    /// Optimizer algorithm name, e.g. "projected-gradient".
    pub optimizer_type: String,
    /// Gradient step size on radii.
    pub step_size: f64,
    /// Minimum allowed hole radius.
    pub r_min: f64,
    /// Maximum allowed hole radius.
    pub r_max: f64,
    /// Number of optimization steps.
    pub steps: usize,
}

impl Default for StudyOptimizer {
    fn default() -> Self {
        Self {
            optimizer_type: "projected-gradient".to_string(),
            step_size: 1.0,
            r_min: 0.08,
            r_max: 0.20,
            steps: 8,
        }
    }
}

fn violation(code: &'static str, what: impl Into<String>, fix: impl Into<String>) -> Violation {
    Violation {
        code,
        what: what.into(),
        fix: fix.into(),
    }
}

impl StudySpec {
    /// Validate all semantic requirements and Five Explicits.
    #[must_use]
    pub fn validate(&self) -> Vec<Violation> {
        let mut violations = Vec::new();

        // 1. Metadata
        if self.metadata.is_none() {
            violations.push(violation(
                "project-missing-section",
                "missing mandatory `metadata` section",
                "declare project metadata (name, context-of-use, decision-gate)",
            ));
        }

        // 2. Five Explicits: versions
        if self.versions.is_none() {
            violations.push(violation(
                "project-missing-section",
                "missing mandatory `versions` section",
                "declare engine and schema versions",
            ));
        }

        // 3. Five Explicits: seeds
        if self.seeds.is_none() {
            violations.push(violation(
                "project-missing-section",
                "missing mandatory `seeds` section",
                "declare pseudo-random generator seeds",
            ));
        }

        // 4. Five Explicits: budgets
        if self.budgets.is_none() {
            violations.push(violation(
                "project-missing-section",
                "missing mandatory `budgets` section",
                "declare explicit execution budgets (wall-time, memory, max-iterations)",
            ));
        }

        // 5. Five Explicits: capabilities
        match &self.capabilities {
            None => {
                violations.push(violation(
                    "project-missing-section",
                    "missing mandatory `capabilities` section",
                    "declare required capability verbs",
                ));
            }
            Some(caps) if caps.is_empty() => {
                violations.push(violation(
                    "project-empty-capabilities",
                    "capabilities section declared no capability verbs",
                    "declare at least one required capability verb",
                ));
            }
            _ => {}
        }

        // 6. Five Explicits: units doctrine
        match &self.units {
            None => {
                violations.push(violation(
                    "project-undeclared-units",
                    "units doctrine must be explicitly declared (Five Explicits P10)",
                    "declare (units :system \"SI\")",
                ));
            }
            Some(u) if u.storage.trim().is_empty() => {
                violations.push(violation(
                    "project-undeclared-units",
                    "units storage convention was declared empty",
                    "declare (units :system \"SI\")",
                ));
            }
            _ => {}
        }

        // 7. Domain
        match &self.domain {
            None => {
                violations.push(violation(
                    "project-missing-section",
                    "missing mandatory `domain` section",
                    "declare the design domain and initial geometry",
                ));
            }
            Some(domain) => {
                if domain.initial_holes.is_empty() {
                    violations.push(violation(
                        "study-empty-initial-holes",
                        "domain declared zero initial holes",
                        "declare at least one initial hole in the design domain",
                    ));
                }
                for (idx, hole) in domain.initial_holes.iter().enumerate() {
                    if hole.radius <= 0.0 || !hole.radius.is_finite() {
                        violations.push(violation(
                            "study-invalid-hole-radius",
                            format!("hole {idx} has invalid radius {}", hole.radius),
                            "declare a positive finite hole radius",
                        ));
                    }
                    if !hole.center[0].is_finite() || !hole.center[1].is_finite() {
                        violations.push(violation(
                            "study-invalid-hole-center",
                            format!("hole {idx} has non-finite center coordinates"),
                            "declare finite coordinates within domain bounds",
                        ));
                    }
                }
            }
        }

        // 8. Physics
        match &self.physics {
            None => {
                violations.push(violation(
                    "project-missing-section",
                    "missing mandatory `physics` section",
                    "declare the physics formulation and discretization settings",
                ));
            }
            Some(physics) => {
                if physics.mesh_level < 1 || physics.mesh_level > 12 {
                    violations.push(violation(
                        "study-invalid-mesh-level",
                        format!("physics mesh level {} is outside [1, 12]", physics.mesh_level),
                        "choose a background grid refinement level between 1 and 12",
                    ));
                }
            }
        }

        // 9. Scenario & Boundary Conditions
        match &self.scenario {
            None => {
                violations.push(violation(
                    "project-missing-section",
                    "missing mandatory `scenario` section",
                    "declare supports and applied loads",
                ));
            }
            Some(scenario) => {
                let valid_boundaries = [
                    "left", "right", "top", "bottom",
                    "boundary/left", "boundary/right", "boundary/top", "boundary/bottom",
                ];
                let load_region = scenario.load_region.trim().to_ascii_lowercase();
                if !valid_boundaries.contains(&load_region.as_str()) {
                    violations.push(violation(
                        "study-load-non-boundary",
                        format!(
                            "load region `{}` is not a boundary region; loads must be applied on a boundary region",
                            scenario.load_region
                        ),
                        "place the applied load on a recognized boundary region (e.g. `left`, `right`, `top`, `bottom`)",
                    ));
                }
            }
        }

        // 10. Objective & Dimensional Checks
        match &self.objective {
            None => {
                violations.push(violation(
                    "project-missing-section",
                    "missing mandatory `objective` section",
                    "declare the optimization objective, sense, and units",
                ));
            }
            Some(objective) => {
                // Compliance functional J(u) = ∫ f·u represents energy/work with dimensions [m^2 kg s^-2].
                let energy_dims = Dims([2, 1, -2, 0, 0, 0]);
                let test_qty_str = format!("1.0 {}", objective.unit.trim());
                match fs_qty::parse::parse_qty(&test_qty_str) {
                    Ok(qty) => {
                        if qty.dims != energy_dims {
                            violations.push(violation(
                                "study-objective-dimension-mismatch",
                                format!(
                                    "objective `{}` requires energy dimensions [m^2 kg s^-2] (unit J or N*m); got `{}` with dimensions {}",
                                    objective.objective_type,
                                    objective.unit,
                                    qty.dims.unit_string()
                                ),
                                "declare the compliance objective in energy units such as `J` or `N*m`",
                            ));
                        }
                    }
                    Err(_) => {
                        violations.push(violation(
                            "study-objective-invalid-unit",
                            format!(
                                "objective `{}` declared unparseable unit `{}`",
                                objective.objective_type, objective.unit
                            ),
                            "provide a valid SI unit symbol like `J`",
                        ));
                    }
                }
            }
        }

        // 11. Constraints
        match &self.constraints {
            None => {
                violations.push(violation(
                    "project-missing-section",
                    "missing mandatory `constraints` section",
                    "declare design constraints (volume-fraction)",
                ));
            }
            Some(constraints) => {
                if constraints.volume_fraction <= 0.0 || constraints.volume_fraction >= 1.0 || !constraints.volume_fraction.is_finite() {
                    violations.push(violation(
                        "study-volume-fraction-out-of-bounds",
                        format!(
                            "volume fraction must be strictly between 0 and 1; got {}",
                            constraints.volume_fraction
                        ),
                        "set the material volume fraction target strictly inside (0, 1)",
                    ));
                }
            }
        }

        // 12. Optimizer
        match &self.optimizer {
            None => {
                violations.push(violation(
                    "project-missing-section",
                    "missing mandatory `optimizer` section",
                    "declare optimizer hyperparameters (step-size, r-min, r-max, steps)",
                ));
            }
            Some(opt) => {
                if opt.step_size <= 0.0 || !opt.step_size.is_finite() {
                    violations.push(violation(
                        "study-invalid-step-size",
                        format!("optimizer step size {} is not positive finite", opt.step_size),
                        "set a positive finite gradient step size",
                    ));
                }
                if opt.r_min <= 0.0 || opt.r_min > opt.r_max || !opt.r_min.is_finite() || !opt.r_max.is_finite() {
                    violations.push(violation(
                        "study-invalid-radius-bounds",
                        format!("radius bounds [{}, {}] violate 0 < r_min <= r_max", opt.r_min, opt.r_max),
                        "set valid positive radius bounds with r_min <= r_max",
                    ));
                }
                if opt.steps == 0 {
                    violations.push(violation(
                        "study-invalid-step-count",
                        "optimizer step count must be at least 1",
                        "declare at least one optimization step",
                    ));
                }
            }
        }

        violations
    }
}

/// Print canonical s-expression syntax.
#[must_use]
pub fn print_study_sexpr(study: &StudySpec) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "(fsim-study :version {}", STUDY_FSIM_VERSION);

    if let Some(m) = &study.metadata {
        let _ = writeln!(out, "  (metadata");
        let _ = writeln!(out, "    :name {:?}", m.name);
        let _ = writeln!(out, "    :created {:?}", m.created);
        let _ = writeln!(out, "    :context-of-use {:?}", m.context_of_use);
        let _ = writeln!(out, "    :intended-decision {:?}", m.intended_decision);
        let _ = writeln!(out, "    :decision-gate {}", match m.decision_gate {
            DecisionGate::ScopingEstimate => "scoping-estimate",
            DecisionGate::DesignSelection => "design-selection",
            DecisionGate::ComplianceSignoff => "compliance-signoff",
        });
        let _ = writeln!(out, "    :consequence {})", match m.consequence {
            ConsequenceClass::Advisory => "advisory",
            ConsequenceClass::Reliability => "reliability",
            ConsequenceClass::SafetyCritical => "safety-critical",
        });
    }

    if let Some(v) = &study.versions {
        let _ = writeln!(out, "  (versions");
        let _ = writeln!(out, "    :schema {})", v.schema);
    }

    if let Some(s) = &study.seeds {
        let _ = writeln!(out, "  (seeds");
        let _ = writeln!(out, "    :root {})", s.root);
    }

    if let Some(b) = &study.budgets {
        let _ = writeln!(out, "  (budgets");
        if let Some(wt) = &b.wall_time {
            let _ = writeln!(out, "    :wall-time {} s", wt.value);
        }
        if let Some(mem) = b.memory_bytes {
            let _ = writeln!(out, "    :memory {} B", mem);
        }
        if let Some(max_iter) = b.max_iterations {
            let _ = writeln!(out, "    :max-iterations {})", max_iter);
        } else {
            let _ = writeln!(out, "  )");
        }
    }

    if let Some(caps) = &study.capabilities {
        let _ = writeln!(out, "  (capabilities");
        for cap in caps {
            let _ = writeln!(out, "    {:?}", cap);
        }
        let _ = writeln!(out, "  )");
    }

    if let Some(u) = &study.units {
        let _ = writeln!(out, "  (units");
        let _ = writeln!(out, "    :storage {:?})", u.storage);
    }

    if let Some(d) = &study.domain {
        let _ = writeln!(out, "  (domain");
        let _ = writeln!(out, "    :type {}", d.domain_type);
        let _ = writeln!(out, "    :bounds (({} {}) ({} {}))", d.bounds.0[0], d.bounds.0[1], d.bounds.1[0], d.bounds.1[1]);
        let _ = writeln!(out, "    :initial-holes (");
        for hole in &d.initial_holes {
            let _ = writeln!(out, "      (hole :center ({} {}) :radius {})", hole.center[0], hole.center[1], hole.radius);
        }
        let _ = writeln!(out, "    ))");
    }

    if let Some(p) = &study.physics {
        let _ = writeln!(out, "  (physics");
        let _ = writeln!(out, "    :type {}", p.physics_type);
        let _ = writeln!(out, "    :mesh-level {})", p.mesh_level);
    }

    if let Some(sc) = &study.scenario {
        let _ = writeln!(out, "  (scenario");
        let _ = writeln!(out, "    :fixed-boundary {:?}", sc.fixed_boundary);
        let _ = writeln!(out, "    :load-region {:?}", sc.load_region);
        if let Some(force) = &sc.load_force {
            let _ = writeln!(out, "    :load-force {} N)", force.value);
        } else {
            let _ = writeln!(out, "  )");
        }
    }

    if let Some(obj) = &study.objective {
        let _ = writeln!(out, "  (objective");
        let _ = writeln!(out, "    :type {}", obj.objective_type);
        let _ = writeln!(out, "    :sense {}", obj.sense);
        let _ = writeln!(out, "    :unit {:?})", obj.unit);
    }

    if let Some(c) = &study.constraints {
        let _ = writeln!(out, "  (constraints");
        let _ = writeln!(out, "    :volume-fraction {})", c.volume_fraction);
    }

    if let Some(opt) = &study.optimizer {
        let _ = writeln!(out, "  (optimizer");
        let _ = writeln!(out, "    :type {}", opt.optimizer_type);
        let _ = writeln!(out, "    :step-size {}", opt.step_size);
        let _ = writeln!(out, "    :r-min {}", opt.r_min);
        let _ = writeln!(out, "    :r-max {}", opt.r_max);
        let _ = writeln!(out, "    :steps {})", opt.steps);
    }

    let _ = writeln!(out, ")");
    out
}

/// Print JSON representation.
#[must_use]
pub fn print_study_json(study: &StudySpec) -> String {
    let sexpr = print_study_sexpr(study);
    if let Ok(node) = fs_ir::sexpr::parse(&sexpr) {
        if let Ok(json_str) = fs_ir::json::print(&node) {
            return json_str;
        }
    }
    String::new()
}

/// Parse study from s-expression source string.
///
/// # Errors
/// Returns [`crate::ProjectError`] on syntax or envelope structure errors.
pub fn parse_study_sexpr(source: &str) -> Result<StudySpec, crate::ProjectError> {
    let node = fs_ir::sexpr::parse(source).map_err(|e| crate::ProjectError {
        code: "study-syntax",
        detail: format!("s-expression parse error: {e}"),
        hint: "check parentheses and literal formatting".to_string(),
    })?;
    parse_study_node(&node)
}

/// Parse study from JSON source string.
///
/// # Errors
/// Returns [`crate::ProjectError`] on JSON syntax errors.
pub fn parse_study_json(source: &str) -> Result<StudySpec, crate::ProjectError> {
    let node = fs_ir::json::parse(source).map_err(|e| crate::ProjectError {
        code: "study-json-syntax",
        detail: format!("JSON parse error: {e}"),
        hint: "check valid JSON syntax".to_string(),
    })?;
    parse_study_node(&node)
}

fn parse_study_node(root: &Node) -> Result<StudySpec, crate::ProjectError> {
    let NodeKind::List(items) = &root.kind else {
        return Err(crate::ProjectError {
            code: "study-not-a-list",
            detail: "study root must be a list".to_string(),
            hint: "wrap document in (fsim-study :version 1 ...)".to_string(),
        });
    };

    let Some(first) = items.first() else {
        return Err(crate::ProjectError {
            code: "study-empty-list",
            detail: "empty study root list".to_string(),
            hint: "document must begin with `fsim-study`".to_string(),
        });
    };

    match &first.kind {
        NodeKind::Symbol(sym) if sym == "fsim-study" || sym == "study" => {}
        _ => {
            return Err(crate::ProjectError {
                code: "study-wrong-root",
                detail: "expected `fsim-study` root symbol".to_string(),
                hint: "document must begin with `fsim-study`".to_string(),
            });
        }
    }

    // Check version: either `:version <int>`
    let mut version = 1u32;
    let mut start_idx = 1;
    if items.len() >= 3 {
        if let (NodeKind::Keyword(k), NodeKind::Int(v)) = (&items[1].kind, &items[2].kind) {
            if k == "version" {
                version = *v as u32;
                start_idx = 3;
            }
        }
    }

    if version != STUDY_FSIM_VERSION {
        return Err(crate::ProjectError {
            code: "study-unsupported-version",
            detail: format!("study declares version {version}; this reader admits only {STUDY_FSIM_VERSION}"),
            hint: "use version 1".to_string(),
        });
    }

    let mut spec = StudySpec::default();

    for item in &items[start_idx..] {
        let NodeKind::List(section_items) = &item.kind else {
            continue;
        };
        let Some(first_node) = section_items.first() else {
            continue;
        };
        let section_name = match &first_node.kind {
            NodeKind::Symbol(s) => s.as_str(),
            _ => continue,
        };

        match section_name {
            "metadata" => {
                let mut name = "bracket-2d".to_string();
                let mut created = "2026-09-01".to_string();
                let mut context_of_use = "marquee study".to_string();
                let mut intended_decision = "material layout".to_string();
                let mut decision_gate = DecisionGate::ScopingEstimate;
                let mut consequence = ConsequenceClass::Advisory;

                let mut idx = 1;
                while idx < section_items.len() {
                    if let NodeKind::Keyword(k) = &section_items[idx].kind {
                        if let Some(val_node) = section_items.get(idx + 1) {
                            match (k.as_str(), &val_node.kind) {
                                ("name", NodeKind::Str(s)) => name = s.clone(),
                                ("created", NodeKind::Str(s)) => created = s.clone(),
                                ("context-of-use", NodeKind::Str(s)) => context_of_use = s.clone(),
                                ("intended-decision", NodeKind::Str(s)) => intended_decision = s.clone(),
                                ("decision-gate", NodeKind::Symbol(s)) => {
                                    if s == "design-selection" {
                                        decision_gate = DecisionGate::DesignSelection;
                                    } else if s == "compliance-signoff" {
                                        decision_gate = DecisionGate::ComplianceSignoff;
                                    }
                                }
                                ("consequence", NodeKind::Symbol(s)) => {
                                    if s == "reliability" {
                                        consequence = ConsequenceClass::Reliability;
                                    } else if s == "safety-critical" {
                                        consequence = ConsequenceClass::SafetyCritical;
                                    }
                                }
                                _ => {}
                            }
                            idx += 2;
                            continue;
                        }
                    }
                    idx += 1;
                }

                spec.metadata = Some(Metadata {
                    name,
                    created,
                    context_of_use,
                    intended_decision,
                    decision_gate,
                    consequence,
                });
            }
            "versions" => {
                let mut schema = 1u32;
                let mut idx = 1;
                while idx < section_items.len() {
                    if let NodeKind::Keyword(k) = &section_items[idx].kind {
                        if let Some(val_node) = section_items.get(idx + 1) {
                            if k == "schema" {
                                if let NodeKind::Int(i) = &val_node.kind {
                                    schema = *i as u32;
                                }
                            }
                            idx += 2;
                            continue;
                        }
                    }
                    idx += 1;
                }
                spec.versions = Some(Versions {
                    schema,
                    constellation: "frankensim-constellation-lock-v1".to_string(),
                    workspace: "frankensim-workspace-v1".to_string(),
                });
            }
            "seeds" => {
                let mut root = 1337u64;
                let mut idx = 1;
                while idx < section_items.len() {
                    if let NodeKind::Keyword(k) = &section_items[idx].kind {
                        if let Some(val_node) = section_items.get(idx + 1) {
                            if k == "rng" || k == "root" {
                                if let NodeKind::Int(i) = &val_node.kind {
                                    root = *i as u64;
                                }
                            }
                            idx += 2;
                            continue;
                        }
                    }
                    idx += 1;
                }
                spec.seeds = Some(Seeds { root });
            }
            "budgets" => {
                let mut wall_time = Some(QtyAny::new(60.0, fs_qty::Dims([0, 0, 1, 0, 0, 0])));
                let mut memory_bytes = Some(1024 * 1024 * 1024);
                let mut max_iterations = Some(8);

                let mut idx = 1;
                while idx < section_items.len() {
                    if let NodeKind::Keyword(k) = &section_items[idx].kind {
                        if let Some(val_node) = section_items.get(idx + 1) {
                            match (k.as_str(), &val_node.kind) {
                                ("wall-time", NodeKind::Int(i)) => {
                                    wall_time = Some(QtyAny::new(*i as f64, fs_qty::Dims([0, 0, 1, 0, 0, 0])));
                                }
                                ("wall-time", NodeKind::Float(f)) => {
                                    wall_time = Some(QtyAny::new(*f, fs_qty::Dims([0, 0, 1, 0, 0, 0])));
                                }
                                ("memory", NodeKind::Int(i)) => {
                                    memory_bytes = Some(*i as u64);
                                }
                                ("max-iterations", NodeKind::Int(i)) => {
                                    max_iterations = Some(*i as usize);
                                }
                                _ => {}
                            }
                            idx += 2;
                            continue;
                        }
                    }
                    idx += 1;
                }
                spec.budgets = Some(StudyBudgets {
                    wall_time,
                    memory_bytes,
                    max_iterations,
                });
            }
            "capabilities" => {
                let mut caps = Vec::new();
                for node in &section_items[1..] {
                    if let NodeKind::Str(s) = &node.kind {
                        caps.push(s.clone());
                    } else if let NodeKind::Symbol(s) = &node.kind {
                        caps.push(s.clone());
                    }
                }
                spec.capabilities = Some(caps);
            }
            "units" => {
                let mut storage = "si-base".to_string();
                let mut idx = 1;
                while idx < section_items.len() {
                    if let NodeKind::Keyword(k) = &section_items[idx].kind {
                        if let Some(val_node) = section_items.get(idx + 1) {
                            if k == "system" || k == "storage" {
                                if let NodeKind::Str(s) = &val_node.kind {
                                    storage = s.clone();
                                }
                            }
                            idx += 2;
                            continue;
                        }
                    }
                    idx += 1;
                }
                spec.units = Some(UnitsDoctrine { storage: storage.clone(), display: storage });
            }
            "domain" => {
                let mut domain_type = "sdf-plate-with-holes".to_string();
                let bounds = ([0.0, 0.0], [1.0, 1.0]);
                let mut initial_holes = Vec::new();

                let mut idx = 1;
                while idx < section_items.len() {
                    if let NodeKind::Keyword(k) = &section_items[idx].kind {
                        if let Some(val_node) = section_items.get(idx + 1) {
                            match (k.as_str(), &val_node.kind) {
                                ("type", NodeKind::Symbol(s)) => domain_type = s.clone(),
                                ("initial-holes", NodeKind::List(holes_list)) => {
                                    for h_item in holes_list {
                                        if let NodeKind::List(h_fields) = &h_item.kind {
                                            let mut center = [0.5, 0.5];
                                            let mut radius = 0.1;
                                            let mut h_idx = 0;
                                            while h_idx < h_fields.len() {
                                                if let NodeKind::Keyword(hk) = &h_fields[h_idx].kind {
                                                    if let Some(h_val) = h_fields.get(h_idx + 1) {
                                                        if hk == "center" {
                                                            if let NodeKind::List(coords) = &h_val.kind {
                                                                if coords.len() >= 2 {
                                                                    let x = match coords[0].kind {
                                                                        NodeKind::Float(f) => f,
                                                                        NodeKind::Int(i) => i as f64,
                                                                        _ => 0.0,
                                                                    };
                                                                    let y = match coords[1].kind {
                                                                        NodeKind::Float(f) => f,
                                                                        NodeKind::Int(i) => i as f64,
                                                                        _ => 0.0,
                                                                    };
                                                                    center = [x, y];
                                                                }
                                                            }
                                                        } else if hk == "radius" {
                                                            radius = match h_val.kind {
                                                                NodeKind::Float(f) => f,
                                                                NodeKind::Int(i) => i as f64,
                                                                _ => 0.1,
                                                            };
                                                        }
                                                        h_idx += 2;
                                                        continue;
                                                    }
                                                }
                                                h_idx += 1;
                                            }
                                            initial_holes.push(StudyHole { center, radius });
                                        }
                                    }
                                }
                                _ => {}
                            }
                            idx += 2;
                            continue;
                        }
                    }
                    idx += 1;
                }
                spec.domain = Some(StudyDomain {
                    domain_type,
                    bounds,
                    initial_holes,
                });
            }
            "physics" => {
                let mut physics_type = "elasticity-2d".to_string();
                let mut mesh_level = 4u32;
                let mut idx = 1;
                while idx < section_items.len() {
                    if let NodeKind::Keyword(k) = &section_items[idx].kind {
                        if let Some(val_node) = section_items.get(idx + 1) {
                            match (k.as_str(), &val_node.kind) {
                                ("type", NodeKind::Symbol(s)) => physics_type = s.clone(),
                                ("mesh-level", NodeKind::Int(i)) => mesh_level = *i as u32,
                                _ => {}
                            }
                            idx += 2;
                            continue;
                        }
                    }
                    idx += 1;
                }
                spec.physics = Some(StudyPhysics {
                    physics_type,
                    mesh_level,
                    youngs_modulus: None,
                    poissons_ratio: Some(0.3),
                });
            }
            "scenario" => {
                let mut fixed_boundary = "left".to_string();
                let mut load_region = "right".to_string();
                let mut idx = 1;
                while idx < section_items.len() {
                    if let NodeKind::Keyword(k) = &section_items[idx].kind {
                        if let Some(val_node) = section_items.get(idx + 1) {
                            match (k.as_str(), &val_node.kind) {
                                ("fixed-boundary", NodeKind::Str(s)) => fixed_boundary = s.clone(),
                                ("fixed-boundary", NodeKind::Symbol(s)) => fixed_boundary = s.clone(),
                                ("load-region", NodeKind::Str(s)) => load_region = s.clone(),
                                ("load-region", NodeKind::Symbol(s)) => load_region = s.clone(),
                                _ => {}
                            }
                            idx += 2;
                            continue;
                        }
                    }
                    idx += 1;
                }
                spec.scenario = Some(StudyScenario {
                    fixed_boundary,
                    load_region,
                    load_force: None,
                });
            }
            "objective" => {
                let mut objective_type = "compliance".to_string();
                let mut sense = "minimize".to_string();
                let mut unit = "J".to_string();
                let mut idx = 1;
                while idx < section_items.len() {
                    if let NodeKind::Keyword(k) = &section_items[idx].kind {
                        if let Some(val_node) = section_items.get(idx + 1) {
                            match (k.as_str(), &val_node.kind) {
                                ("type", NodeKind::Symbol(s)) => objective_type = s.clone(),
                                ("sense", NodeKind::Symbol(s)) => sense = s.clone(),
                                ("unit", NodeKind::Str(s)) => unit = s.clone(),
                                _ => {}
                            }
                            idx += 2;
                            continue;
                        }
                    }
                    idx += 1;
                }
                spec.objective = Some(StudyObjective {
                    objective_type,
                    sense,
                    unit,
                });
            }
            "constraints" => {
                let mut volume_fraction = 0.5;
                let mut idx = 1;
                while idx < section_items.len() {
                    if let NodeKind::Keyword(k) = &section_items[idx].kind {
                        if let Some(val_node) = section_items.get(idx + 1) {
                            if k == "volume-fraction" {
                                volume_fraction = match val_node.kind {
                                    NodeKind::Float(f) => f,
                                    NodeKind::Int(i) => i as f64,
                                    _ => 0.5,
                                };
                            }
                            idx += 2;
                            continue;
                        }
                    }
                    idx += 1;
                }
                spec.constraints = Some(StudyConstraints { volume_fraction });
            }
            "optimizer" => {
                let mut optimizer_type = "projected-gradient".to_string();
                let mut step_size = 1.0;
                let mut r_min = 0.08;
                let mut r_max = 0.20;
                let mut steps = 8;
                let mut idx = 1;
                while idx < section_items.len() {
                    if let NodeKind::Keyword(k) = &section_items[idx].kind {
                        if let Some(val_node) = section_items.get(idx + 1) {
                            match (k.as_str(), &val_node.kind) {
                                ("type", NodeKind::Symbol(s)) => optimizer_type = s.clone(),
                                ("step-size", NodeKind::Float(f)) => step_size = *f,
                                ("step-size", NodeKind::Int(i)) => step_size = *i as f64,
                                ("r-min", NodeKind::Float(f)) => r_min = *f,
                                ("r-min", NodeKind::Int(i)) => r_min = *i as f64,
                                ("r-max", NodeKind::Float(f)) => r_max = *f,
                                ("r-max", NodeKind::Int(i)) => r_max = *i as f64,
                                ("steps", NodeKind::Int(i)) => steps = *i as usize,
                                _ => {}
                            }
                            idx += 2;
                            continue;
                        }
                    }
                    idx += 1;
                }
                spec.optimizer = Some(StudyOptimizer {
                    optimizer_type,
                    step_size,
                    r_min,
                    r_max,
                    steps,
                });
            }
            _ => {}
        }
    }

    Ok(spec)
}
