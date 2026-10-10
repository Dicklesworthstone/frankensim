//! Explicit, bounded 3-D input. The 2-D study schema is deliberately separate.
use super::*;
use fs_ascent::projected_al::ProjectedAlOptions;
use fs_ir::ast::{Node, NodeKind};
use fs_topopt::sdf3::design::StressDesignOptions3;
use fs_topopt::sdf3::stress::StressOptions3;

#[derive(Clone, Debug)]
pub(super) struct Spec {
    pub canonical: String,
    pub id: ContentHash,
    pub seed: u64,
    pub wall_s: f64,
    pub memory: usize,
    pub linear: usize,
    pub per_solve: usize,
    pub boxes: usize,
    pub points: usize,
    pub physical: bool,
    pub bounds: ([f64; 3], [f64; 3]),
    pub fixed: FixedFace,
    pub height: f64,
    pub curvature: f64,
    pub level: u32,
    pub max_level: u32,
    pub leaves: usize,
    pub youngs: f64,
    pub poisson: f64,
    pub loads: Vec<([f64; 3], f64)>,
    pub surfaces: Vec<Option<loading::SurfaceSpec>>,
    pub density: f64,
    pub volume: f64,
    pub radius: f64,
    pub updates: usize,
    pub move_limit: f64,
    pub marking: f64,
    pub max_marks: usize,
    pub schedule: Vec<SimpParams>,
    /// A fixed-background minimum-volume problem; adaptive controls apply only
    /// to the compliance mode when this is absent.
    pub stress: Option<StressDesignOptions3>,
}

fn invalid(message: impl Into<String>) -> Failure {
    fail("cli-study-sdf3-input", message)
}

fn fields<'a>(node: &'a Node, name: &str, keys: &[&str]) -> Result<Vec<&'a Node>> {
    let items = super::super::list(node, "3-D study section")?;
    if items.len() != 1 + 2 * keys.len()
        || !matches!(&items[0].kind, NodeKind::Symbol(value) if value == name)
    {
        return Err(invalid(format!(
            "{name} requires exactly the declared fields, in order"
        )));
    }
    keys.iter()
        .enumerate()
        .map(|(index, key)| {
            if !matches!(&items[1 + 2 * index].kind, NodeKind::Keyword(value) if value == key) {
                return Err(invalid(format!(
                    "{name} requires :{key} at position {}",
                    index + 1
                )));
            }
            Ok(&items[2 + 2 * index])
        })
        .collect()
}

fn word(node: &Node, expected: &str) -> Result<()> {
    if matches!(&node.kind, NodeKind::Symbol(value) | NodeKind::Str(value) if value == expected) {
        Ok(())
    } else {
        Err(invalid(format!("expected {expected}")))
    }
}

fn scalar(node: &Node) -> Result<f64> {
    super::super::number(node, "3-D scalar")
}
fn count(node: &Node) -> Result<usize> {
    super::super::integer_node(node, "3-D count")
}
fn vector(node: &Node) -> Result<[f64; 3]> {
    let items = super::super::list(node, "3-D vector")?;
    if items.len() != 3 {
        return Err(invalid("expected three vector components"));
    }
    Ok([scalar(&items[0])?, scalar(&items[1])?, scalar(&items[2])?])
}

fn surface(node: &Node) -> Result<Option<loading::SurfaceSpec>> {
    if matches!(&node.kind, NodeKind::Symbol(s) if s == "none") {
        return Ok(None);
    }
    let items = super::super::list(node, "reference surface load")?;
    let (name, pressure) = match items.first().map(|n| &n.kind) {
        Some(NodeKind::Symbol(s)) if s == "pressure" => ("pressure", true),
        Some(NodeKind::Symbol(s)) if s == "traction" => ("traction", false),
        _ => {
            return Err(invalid(
                "surface requires none, (pressure ...) or (traction ...)",
            ));
        }
    };
    let values = fields(node, name, &["pa", "x-fraction"])?;
    let force = if pressure {
        loading::Force::Pressure(scalar(values[0])?)
    } else {
        loading::Force::Traction(vector(values[0])?)
    };
    Ok(Some(loading::SurfaceSpec {
        force,
        x_fraction: super::super::pair(values[1], "surface x-fraction")?,
    }))
}

pub(super) fn parse(source: &str) -> Result<Spec> {
    let root = fs_ir::sexpr::parse(source).map_err(|e| invalid(e.to_string()))?;
    let items = super::super::list(&root, "3-D study")?;
    if items.len() != 13
        || !matches!(&items[0].kind, NodeKind::Symbol(s) if s == "fsim-sdf3-study")
        || !matches!(&items[1].kind, NodeKind::Keyword(s) if s == "version")
        || !matches!(items[2].kind, NodeKind::Int(1))
    {
        return Err(invalid(
            "expected the complete fsim-sdf3-study version 1 schema",
        ));
    }
    let meta = fields(
        &items[3],
        "metadata",
        &[
            "name",
            "created",
            "context-of-use",
            "intended-decision",
            "decision-gate",
            "consequence",
        ],
    )?;
    for value in &meta[..4] {
        if !matches!(&value.kind, NodeKind::Str(s) if !s.is_empty() && s.len() <= 4096) {
            return Err(invalid("metadata strings must contain 1..=4096 bytes"));
        }
    }
    word(meta[4], "scoping-estimate")?;
    word(meta[5], "advisory")?;
    let seed = count(fields(&items[4], "seeds", &["root"])?[0])? as u64;
    let budgets = fields(
        &items[5],
        "budgets",
        &[
            "wall-seconds",
            "memory-bytes",
            "linear-iterations",
            "per-solve-iterations",
            "quadrature-boxes",
            "quadrature-points",
        ],
    )?;
    let capabilities = super::super::list(&items[6], "capabilities")?;
    if capabilities.len() != 4 {
        return Err(invalid("three explicit capabilities are required"));
    }
    for (value, expected) in capabilities.iter().zip([
        "capabilities",
        "optimization.marquee-topopt",
        "geometry.sdf",
        "physics.cutfem",
    ]) {
        word(value, expected)?;
    }
    word(fields(&items[7], "units", &["storage"])?[0], "SI")?;
    let domain = fields(
        &items[8],
        "domain",
        &["type", "bounds", "height-m", "curvature-per-m"],
    )?;
    let physical = matches!(&domain[0].kind,
        NodeKind::Symbol(s) | NodeKind::Str(s) if s == "physical-curved-height-sdf");
    if !physical {
        word(domain[0], "curved-height-sdf")?;
    }
    let bounds = super::super::list(domain[1], "bounds")?;
    if bounds.len() != 2 {
        return Err(invalid("bounds require lower and upper SI vectors"));
    }
    let bounds = (vector(&bounds[0])?, vector(&bounds[1])?);
    let physics = fields(
        &items[9],
        "physics",
        &[
            "type",
            "initial-level",
            "maximum-level",
            "maximum-leaves",
            "youngs-modulus-pa",
            "poissons-ratio",
        ],
    )?;
    word(physics[0], "elasticity-3d")?;
    let scenario_items = super::super::list(&items[10], "scenario")?;
    let mixed = scenario_items
        .get(3)
        .is_some_and(|node| matches!(&node.kind, NodeKind::Keyword(s) if s == "loads"));
    let scenario = fields(
        &items[10],
        "scenario",
        &["fixed-boundary", if mixed { "loads" } else { "body-loads" }],
    )?;
    let fixed = match &scenario[0].kind {
        NodeKind::Symbol(s) | NodeKind::Str(s) => match s.as_str() {
            "left" => FixedFace::Left,
            "right" => FixedFace::Right,
            "front" => FixedFace::Front,
            "back" => FixedFace::Back,
            "bottom" => FixedFace::Bottom,
            _ => {
                return Err(invalid(
                    "fixed-boundary requires left, right, front, back or bottom",
                ));
            }
        },
        _ => return Err(invalid("fixed-boundary must name a box face")),
    };
    let load_nodes = super::super::list(scenario[1], "independent loads")?;
    if load_nodes.is_empty() || load_nodes.len() > 4 {
        return Err(invalid("declare 1..=4 independent reference loads"));
    }
    let mut loads = Vec::with_capacity(load_nodes.len());
    let mut surfaces = Vec::with_capacity(load_nodes.len());
    for load in load_nodes {
        if mixed {
            let values = fields(load, "load", &["body-n-m3", "surface", "weight"])?;
            loads.push((vector(values[0])?, scalar(values[2])?));
            surfaces.push(surface(values[1])?);
        } else {
            let values = fields(load, "load", &["density-n-m3", "weight"])?;
            loads.push((vector(values[0])?, scalar(values[1])?));
            surfaces.push(None);
        }
    }
    let objective = fields(&items[11], "objective", &["type", "sense", "unit"])?;
    let stress_mode = matches!(&objective[0].kind,
        NodeKind::Symbol(s) | NodeKind::Str(s) if s == "volume-fraction");
    let objective_words = if stress_mode {
        ["volume-fraction", "minimize", "1"]
    } else {
        ["compliance", "minimize", "J"]
    };
    for (value, expected) in objective.iter().zip(objective_words) {
        word(value, expected)?;
    }
    let (density, volume, radius, updates, move_limit, marking, max_marks, schedule, stress) =
        if stress_mode {
            let opt = fields(
                &items[12],
                "optimizer",
                &[
                    "type",
                    "initial-density",
                    "filter-radius-m",
                    "max-updates",
                    "density-floor",
                    "stress-limit-pa",
                    "relaxation-power",
                    "aggregation-power",
                    "max-stress-points",
                    "max-evaluations",
                    "max-backtracks",
                    "tolerance",
                    "penal",
                    "beta",
                ],
            )?;
            word(opt[0], "stress-limited-simp")?;
            let options = StressDesignOptions3 {
                stress: StressOptions3 {
                    relaxation_power: scalar(opt[6])?,
                    aggregation_power: scalar(opt[7])?,
                    max_cases: loads.len(),
                    max_points: count(opt[8])?,
                },
                stress_limit: scalar(opt[5])?,
                density_floor: scalar(opt[4])?,
                optimizer: ProjectedAlOptions {
                    max_dimension: count(physics[3])?,
                    max_evaluations: count(opt[9])?,
                    max_backtracks: count(opt[10])?,
                    tolerance: scalar(opt[11])?,
                    ..Default::default()
                },
            };
            let params = SimpParams {
                penal: scalar(opt[12])?,
                beta: scalar(opt[13])?,
                ..Default::default()
            };
            // These neutral values are unused by the explicitly separate stress
            // driver; there is no hidden volume target or refinement policy.
            (
                scalar(opt[1])?,
                1.0,
                scalar(opt[2])?,
                count(opt[3])?,
                0.0,
                0.0,
                0,
                vec![params],
                Some(options),
            )
        } else {
            let opt = fields(
                &items[12],
                "optimizer",
                &[
                    "type",
                    "initial-density",
                    "volume-fraction",
                    "filter-radius-m",
                    "updates-per-stage",
                    "move-limit",
                    "marking-fraction",
                    "max-marks",
                    "schedule",
                ],
            )?;
            word(opt[0], "adaptive-simp")?;
            let stage_nodes = super::super::list(opt[8], "SIMP schedule")?;
            if stage_nodes.is_empty() || stage_nodes.len() > 4 {
                return Err(invalid("declare 1..=4 SIMP stages"));
            }
            let mut schedule = Vec::with_capacity(stage_nodes.len());
            for stage in stage_nodes {
                let pair = super::super::pair(stage, "(penal beta)")?;
                schedule.push(SimpParams {
                    penal: pair[0],
                    beta: pair[1],
                    ..Default::default()
                });
            }
            (
                scalar(opt[1])?,
                scalar(opt[2])?,
                scalar(opt[3])?,
                count(opt[4])?,
                scalar(opt[5])?,
                scalar(opt[6])?,
                count(opt[7])?,
                schedule,
                None,
            )
        };
    let canonical = fs_ir::sexpr::print(&root).map_err(|e| invalid(e.to_string()))?;
    let id = hash_domain(
        "org.frankensim.cli.sdf3-study.v1",
        format!(
            "{}\n{}\n{canonical}",
            if stress_mode {
                STRESS3_DRIVER
            } else {
                SDF3_DRIVER
            },
            hash_bytes(include_bytes!("../../../../../constellation.lock")).to_hex()
        )
        .as_bytes(),
    );
    let spec = Spec {
        canonical,
        id,
        seed,
        wall_s: scalar(budgets[0])?,
        memory: count(budgets[1])?,
        linear: count(budgets[2])?,
        per_solve: count(budgets[3])?,
        boxes: count(budgets[4])?,
        points: count(budgets[5])?,
        physical,
        bounds,
        fixed,
        height: scalar(domain[2])?,
        curvature: scalar(domain[3])?,
        level: u32::try_from(count(physics[1])?).map_err(|_| invalid("initial-level overflow"))?,
        max_level: u32::try_from(count(physics[2])?)
            .map_err(|_| invalid("maximum-level overflow"))?,
        leaves: count(physics[3])?,
        youngs: scalar(physics[4])?,
        poisson: scalar(physics[5])?,
        loads,
        surfaces,
        density,
        volume,
        radius,
        updates,
        move_limit,
        marking,
        max_marks,
        schedule,
        stress,
    };
    spec.validate()?;
    Ok(spec)
}

impl Spec {
    fn validate(&self) -> Result<()> {
        if !(0.0 < self.wall_s && self.wall_s <= 86_400.0)
            || !(128 * 1024 * 1024..=8 * 1024 * 1024 * 1024).contains(&self.memory)
            || !(1..=2_000_000).contains(&self.linear)
            || !(1..=100_000).contains(&self.per_solve)
            || !(1..=2_000_000).contains(&self.boxes)
            || !(1..=4_000_000).contains(&self.points)
        {
            return Err(invalid(
                "wall, memory, Krylov or quadrature budget exceeds the admitted envelope",
            ));
        }
        if !(1..=2).contains(&self.level)
            || !(self.level..=5).contains(&self.max_level)
            || !(8usize.pow(self.level)..=2048).contains(&self.leaves)
            || (self.stress.is_none() && !(1..=8).contains(&self.max_marks))
        {
            return Err(invalid(
                "initial octree level must be 1..=2, maximum 5, at most 2048 leaves and 1..=8 marks",
            ));
        }
        // Admission envelope, not a measured or certified peak-RSS bound. Covers
        // retained quadrature, simultaneous proposal fields and bounded setup.
        let stress_memory = self
            .stress
            .map_or(0, |s| 1024usize.saturating_mul(s.stress.max_points));
        let admitted_memory = (64 * 1024 * 1024 + 256 * self.points + 65_536 * self.leaves)
            .saturating_add(stress_memory);
        if admitted_memory > self.memory {
            return Err(invalid(
                "declared memory is too small for the requested geometry/field envelope",
            ));
        }
        geometry::validate(self)?;
        if !(1e-6..=1e15).contains(&self.youngs) || !(0.0..=1.0 / 3.0).contains(&self.poisson) {
            return Err(invalid(
                "unsupported finite graph-domain or isotropic material parameters",
            ));
        }
        if self.surfaces.len() != self.loads.len() {
            return Err(invalid("surface and body load families must align"));
        }
        for ((f, w), surface) in self.loads.iter().zip(&self.surfaces) {
            if !(0.0 < *w && *w <= 1.0)
                || f.iter().any(|v| !v.is_finite() || v.abs() > 1e12)
                || (f.iter().all(|v| *v == 0.0) && surface.is_none())
            {
                return Err(invalid(
                    "each case requires a nonzero reference load and weight in (0,1]",
                ));
            }
            if let Some(surface) = surface {
                surface.validate(self.level)?;
            }
        }
        if self
            .schedule
            .iter()
            .any(|s| !(1.0..=5.0).contains(&s.penal) || !(0.0..=16.0).contains(&s.beta))
        {
            return Err(invalid(
                "SIMP penal must be 1..=5 and projection beta 0..=16",
            ));
        }
        if let Some(stress) = self.stress {
            stress
                .optimizer
                .validate::<fs_topopt::sdf3::stress::StressError3>(self.leaves)
                .map_err(|e| invalid(e.to_string()))?;
            let params = self.schedule[0];
            let q = stress.stress.relaxation_power;
            let floor = fs_topopt::heaviside(stress.density_floor, params.beta, params.eta);
            let turnover = q * params.e_min / ((params.penal - q) * (1.0 - params.e_min));
            if self.max_level != self.level
                || !(1..=1024).contains(&self.updates)
                || !(0.01..=0.9).contains(&stress.density_floor)
                || !(stress.density_floor..=1.0).contains(&self.density)
                || !(1e-12..=1e15).contains(&stress.stress_limit)
                || !q.is_finite()
                || q < 1.0
                || q >= params.penal
                || !(2.0..=64.0).contains(&stress.stress.aggregation_power)
                || !(1..=1_000_000).contains(&stress.stress.max_points)
                || !(1..=10_000).contains(&stress.optimizer.max_evaluations)
                || !(1e-10..=1e-2).contains(&stress.optimizer.tolerance)
                || !floor.is_finite()
                || floor <= 0.0
                || !(fs_math::det::pow(floor, params.penal) > turnover)
            {
                return Err(invalid(
                    "stress-limited-simp requires a fixed initial octree, 1..=1024 updates, finite positive aggregate cap, 1<=q<penal, bounded stress/evaluation work, and a density floor above ersatz stress foldback",
                ));
            }
        } else if !(0.01..=0.99).contains(&self.density)
            || !(0.01..=0.99).contains(&self.volume)
            || !(1..=16).contains(&self.updates)
            || !(0.001..=0.5).contains(&self.move_limit)
            || !(0.0 < self.marking && self.marking <= 1.0)
        {
            return Err(invalid(
                "invalid raw density, volume, filter, update or marking control",
            ));
        }
        Ok(())
    }
}
