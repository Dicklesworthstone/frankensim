//! Canonical `.fsim` free-boundary elasticity/topology study producer.
//!
//! This driver shares the `study` command, ledger artifact kinds, report/package
//! surface, and receipt envelope with the retained thermal study driver, while
//! keeping a distinct producer identity and explicit numerical no-claim boundary.
//! Resume restores exact geometry, multiplier and global iteration state under
//! the same executable. Legacy receipts receive one verified prefix replay.

use std::fmt::Write as _;
use std::io::Read;
use std::path::Path;
use std::time::Instant;

use fs_blake3::{ContentHash, hash_bytes, hash_domain};
use fs_exec::CancelGate;
use fs_ir::ast::{Node, NodeKind};
use fs_ledger::{EdgeRole, FiveExplicits, Ledger, LedgerError, OpOutcome};
use fs_package::{EvidencePackage, Provenance};
use fs_project::{ConsequenceClass, DecisionGate};
use fs_topols::{Cantilever, GridSdf, OptimizeReport, OptimizeSettings, optimize_compliance};

use crate::json_read::JsonValue;
use crate::{
    CommandOutput, Diagnostic, MAX_PROJECT_BYTES, OutputMode, exit, push_json_string, refusal,
};
use super::STUDY_RUN_RECEIPT_SCHEMA;

#[path = "study_elasticity/continuation.rs"]
mod continuation;
use continuation::drive;

const DRIVER: &str = "free-boundary-elasticity-study-v1";
const RECEIPT_KIND: &str = "study-run-receipt";
const MAX_ARTIFACT_BYTES: u64 = 16 * 1024 * 1024;
const TRACE_DOMAIN: &str = "org.frankensim.cli.elasticity-study.trace.v1";
const ID_DOMAIN: &str = "org.frankensim.cli.free-boundary-elasticity-study.v1";
const NO_CLAIM: &str = "2-D plane-strain CutFEM compliance on an evolving bilinear level set. Each retained trajectory row binds post-evolution compliance, cut-quadrature material area and geometry snapshot to the same canonically re-solved design. This is Estimated numerical design evidence: no physical validation, continuous stress certificate, KKT/global-optimum claim, 3-D claim, or guaranteed discretization-error bound. Cancellation and wall enforcement occur between accepted-state updates, not inside one CutFEM solve. Each accepted update is durably retained before more physics.";

type Result<T> = std::result::Result<T, Failure>;

#[derive(Debug)]
struct Failure {
    code: &'static str,
    message: String,
    exit: u8,
}

impl From<LedgerError> for Failure {
    fn from(error: LedgerError) -> Self {
        fail("cli-study-elasticity-ledger", error.to_string())
    }
}

fn fail(code: &'static str, message: impl Into<String>) -> Failure {
    Failure { code, message: message.into(), exit: exit::REFUSED }
}

fn quoted(value: &str) -> String {
    let mut text = String::new();
    push_json_string(&mut text, value);
    text
}

fn failure_output(command: &'static str, mode: OutputMode, error: Failure) -> CommandOutput {
    refusal(
        mode,
        error.exit,
        &Diagnostic::new(
            command,
            error.code,
            error.message,
            "use examples/marquee/bracket-2d.fsim and preserve every explicit elasticity-study field",
        ),
        None,
    )
}

#[derive(Debug, Clone)]
struct ElasticitySpec {
    base: fs_project::study::StudySpec,
    canonical: String,
    youngs_pa: f64,
    poisson: f64,
    load_traction_pa: f64,
    load_direction: [f64; 2],
    load_band: [f64; 2],
    move_cells: f64,
    band_cells: f64,
    ell0: f64,
    mu_al: f64,
    sobolev_alpha: f64,
    nucleation_period: usize,
    hole_radius_cells: f64,
    steps: usize,
    wall_s: f64,
    memory_bytes: u64,
    max_iterations: usize,
    id: ContentHash,
}

fn list<'a>(node: &'a Node, what: &'static str) -> Result<&'a [Node]> {
    match &node.kind {
        NodeKind::List(items) => Ok(items),
        _ => Err(fail("cli-study-elasticity-shape", format!("{what} must be a list"))),
    }
}

fn section<'a>(root: &'a Node, name: &'static str) -> Result<&'a [Node]> {
    let items = list(root, "study root")?;
    for item in items {
        let NodeKind::List(fields) = &item.kind else { continue };
        if fields.first().is_some_and(
            |first| matches!(&first.kind, NodeKind::Symbol(value) if value == name),
        ) {
            return Ok(fields);
        }
    }
    Err(fail("cli-study-elasticity-section", format!("missing ({name} ...) section")))
}

fn field<'a>(fields: &'a [Node], key: &'static str) -> Result<&'a Node> {
    for pair in fields.windows(2) {
        if matches!(&pair[0].kind, NodeKind::Keyword(value) if value == key) {
            return Ok(&pair[1]);
        }
    }
    Err(fail("cli-study-elasticity-field", format!("missing explicit :{key}")))
}

fn number(node: &Node, what: &'static str) -> Result<f64> {
    let value = match node.kind {
        NodeKind::Float(value) => value,
        NodeKind::Int(value) => value as f64,
        _ => return Err(fail("cli-study-elasticity-number", format!("{what} must be numeric"))),
    };
    if !value.is_finite() {
        return Err(fail("cli-study-elasticity-number", format!("{what} must be finite")));
    }
    Ok(value)
}

fn integer_node(node: &Node, what: &'static str) -> Result<usize> {
    match node.kind {
        NodeKind::Int(value) if value >= 0 => usize::try_from(value)
            .map_err(|_| fail("cli-study-elasticity-integer", format!("{what} is too large"))),
        _ => Err(fail("cli-study-elasticity-integer", format!("{what} must be a nonnegative integer"))),
    }
}

fn pair(node: &Node, what: &'static str) -> Result<[f64; 2]> {
    let values = list(node, what)?;
    if values.len() != 2 {
        return Err(fail("cli-study-elasticity-pair", format!("{what} requires exactly two numbers")));
    }
    Ok([number(&values[0], what)?, number(&values[1], what)?])
}

fn canonical_float(value: f64) -> String {
    format!("{value:.17e}")
}

fn canonical(spec: &ElasticitySpec) -> String {
    let base = &spec.base;
    let metadata = base.metadata.as_ref().expect("admitted metadata");
    let versions = base.versions.as_ref().expect("admitted versions");
    let seeds = base.seeds.as_ref().expect("admitted seeds");
    let caps = base.capabilities.as_ref().expect("admitted capabilities");
    let units = base.units.as_ref().expect("admitted units");
    let domain = base.domain.as_ref().expect("admitted domain");
    let objective = base.objective.as_ref().expect("admitted objective");
    let constraints = base.constraints.as_ref().expect("admitted constraints");
    let physics = base.physics.as_ref().expect("admitted physics");
    let mut out = String::new();
    let _ = writeln!(out, "(fsim-study :version 1");
    let _ = writeln!(out, "  (metadata");
    let _ = writeln!(out, "    :name {:?}", metadata.name);
    let _ = writeln!(out, "    :created {:?}", metadata.created);
    let _ = writeln!(out, "    :context-of-use {:?}", metadata.context_of_use);
    let _ = writeln!(out, "    :intended-decision {:?}", metadata.intended_decision);
    let _ = writeln!(out, "    :decision-gate {}", metadata.decision_gate.slug());
    let _ = writeln!(out, "    :consequence {})", metadata.consequence.slug());
    let _ = writeln!(out, "  (versions :schema {})", versions.schema);
    let _ = writeln!(out, "  (seeds :root {})", seeds.root);
    let _ = writeln!(out, "  (budgets");
    let _ = writeln!(out, "    :wall-time {} s", canonical_float(spec.wall_s));
    let _ = writeln!(out, "    :memory {} B", spec.memory_bytes);
    let _ = writeln!(out, "    :max-iterations {})", spec.max_iterations);
    let _ = writeln!(out, "  (capabilities");
    for cap in caps {
        let _ = writeln!(out, "    {:?}", cap);
    }
    let _ = writeln!(out, "  )");
    let _ = writeln!(out, "  (units :storage {:?})", units.storage);
    let _ = writeln!(out, "  (domain");
    let _ = writeln!(out, "    :type {}", domain.domain_type);
    let _ = writeln!(out, "    :bounds ((0.0 0.0) (1.0 1.0))");
    let _ = writeln!(out, "    :initial-holes (");
    for hole in &domain.initial_holes {
        let _ = writeln!(
            out,
            "      (hole :center ({} {}) :radius {})",
            canonical_float(hole.center[0]),
            canonical_float(hole.center[1]),
            canonical_float(hole.radius)
        );
    }
    let _ = writeln!(out, "    ))");
    let _ = writeln!(out, "  (physics");
    let _ = writeln!(out, "    :type {}", physics.physics_type);
    let _ = writeln!(out, "    :mesh-level {}", physics.mesh_level);
    let _ = writeln!(out, "    :youngs-modulus-pa {}", canonical_float(spec.youngs_pa));
    let _ = writeln!(out, "    :poissons-ratio {})", canonical_float(spec.poisson));
    let _ = writeln!(out, "  (scenario");
    let _ = writeln!(out, "    :fixed-boundary left");
    let _ = writeln!(out, "    :load-region right");
    let _ = writeln!(out, "    :load-traction-pa {}", canonical_float(spec.load_traction_pa));
    let _ = writeln!(
        out,
        "    :load-direction ({} {})",
        canonical_float(spec.load_direction[0]),
        canonical_float(spec.load_direction[1])
    );
    let _ = writeln!(
        out,
        "    :load-band ({} {}))",
        canonical_float(spec.load_band[0]),
        canonical_float(spec.load_band[1])
    );
    let _ = writeln!(
        out,
        "  (objective :type {} :sense {} :unit {:?})",
        objective.objective_type,
        objective.sense,
        objective.unit
    );
    let _ = writeln!(
        out,
        "  (constraints :volume-fraction {})",
        canonical_float(constraints.volume_fraction)
    );
    let _ = writeln!(out, "  (optimizer");
    let _ = writeln!(out, "    :type level-set-compliance");
    let _ = writeln!(out, "    :move-cells {}", canonical_float(spec.move_cells));
    let _ = writeln!(out, "    :band-cells {}", canonical_float(spec.band_cells));
    let _ = writeln!(out, "    :ell0 {}", canonical_float(spec.ell0));
    let _ = writeln!(out, "    :mu-al {}", canonical_float(spec.mu_al));
    let _ = writeln!(out, "    :sobolev-alpha {}", canonical_float(spec.sobolev_alpha));
    let _ = writeln!(out, "    :nucleation-period {}", spec.nucleation_period);
    let _ = writeln!(
        out,
        "    :hole-radius-cells {}",
        canonical_float(spec.hole_radius_cells)
    );
    let _ = writeln!(out, "    :steps {})", spec.steps);
    let _ = writeln!(out, ")");
    out
}

fn parse(source: &str) -> Result<ElasticitySpec> {
    let base = fs_project::study::parse_study_sexpr(source)
        .map_err(|error| fail(error.code, error.detail))?;
    if let Some(violation) = base.validate().first() {
        return Err(fail(violation.code, &violation.what));
    }
    let root = fs_ir::sexpr::parse(source)
        .map_err(|error| fail("cli-study-elasticity-syntax", error.to_string()))?;
    let root_items = list(&root, "study root")?;
    if root_items.len() < 3
        || !matches!(&root_items[0].kind, NodeKind::Symbol(value) if value == "fsim-study")
        || !matches!(&root_items[1].kind, NodeKind::Keyword(value) if value == "version")
        || !matches!(root_items[2].kind, NodeKind::Int(1))
    {
        return Err(fail(
            "cli-study-elasticity-version",
            "elasticity studies require `(fsim-study :version 1 ...)`",
        ));
    }
    let physics_fields = section(&root, "physics")?;
    let scenario_fields = section(&root, "scenario")?;
    let optimizer_fields = section(&root, "optimizer")?;
    let youngs_pa = number(field(physics_fields, "youngs-modulus-pa")?, "youngs-modulus-pa")?;
    let poisson = number(field(physics_fields, "poissons-ratio")?, "poissons-ratio")?;
    let load_traction_pa = number(field(scenario_fields, "load-traction-pa")?, "load-traction-pa")?;
    let load_direction = pair(field(scenario_fields, "load-direction")?, "load-direction")?;
    let load_band = pair(field(scenario_fields, "load-band")?, "load-band")?;
    let move_cells = number(field(optimizer_fields, "move-cells")?, "move-cells")?;
    let band_cells = number(field(optimizer_fields, "band-cells")?, "band-cells")?;
    let ell0 = number(field(optimizer_fields, "ell0")?, "ell0")?;
    let mu_al = number(field(optimizer_fields, "mu-al")?, "mu-al")?;
    let sobolev_alpha = number(field(optimizer_fields, "sobolev-alpha")?, "sobolev-alpha")?;
    let nucleation_period = integer_node(
        field(optimizer_fields, "nucleation-period")?,
        "nucleation-period",
    )?;
    let hole_radius_cells = number(
        field(optimizer_fields, "hole-radius-cells")?,
        "hole-radius-cells",
    )?;
    let steps = integer_node(field(optimizer_fields, "steps")?, "steps")?;
    let budgets = base.budgets.as_ref().expect("validated budgets");
    let wall_s = budgets
        .wall_time
        .as_ref()
        .ok_or_else(|| fail("cli-study-elasticity-budget", "wall-time is required"))?
        .value;
    let memory_bytes = budgets
        .memory_bytes
        .ok_or_else(|| fail("cli-study-elasticity-budget", "memory is required"))?;
    let max_iterations = budgets
        .max_iterations
        .ok_or_else(|| fail("cli-study-elasticity-budget", "max-iterations is required"))?;
    let mut parsed = ElasticitySpec {
        base,
        canonical: String::new(),
        youngs_pa,
        poisson,
        load_traction_pa,
        load_direction,
        load_band,
        move_cells,
        band_cells,
        ell0,
        mu_al,
        sobolev_alpha,
        nucleation_period,
        hole_radius_cells,
        steps,
        wall_s,
        memory_bytes,
        max_iterations,
        id: ContentHash([0; 32]),
    };
    validate_model(&parsed)?;
    parsed.canonical = canonical(&parsed);
    let canonical_node = fs_ir::sexpr::parse(&parsed.canonical)
        .map_err(|error| fail("cli-study-elasticity-canonical", error.to_string()))?;
    let source_printed = fs_ir::sexpr::print(&root)
        .map_err(|error| fail("cli-study-elasticity-canonical", error.to_string()))?;
    let canonical_printed = fs_ir::sexpr::print(&canonical_node)
        .map_err(|error| fail("cli-study-elasticity-canonical", error.to_string()))?;
    if source_printed != canonical_printed {
        return Err(fail(
            "cli-study-elasticity-noncanonical",
            "elasticity study contains a missing, unsupported, reordered, aliased or defaulted field",
        ));
    }
    let identity = format!(
        "{DRIVER}\n{}\n{}",
        hash_bytes(include_bytes!("../../../constellation.lock")).to_hex(),
        parsed.canonical
    );
    parsed.id = hash_domain(ID_DOMAIN, identity.as_bytes());
    Ok(parsed)
}

fn validate_model(spec: &ElasticitySpec) -> Result<()> {
    let base = &spec.base;
    let metadata = base.metadata.as_ref().expect("validated metadata");
    if metadata.decision_gate != DecisionGate::ScopingEstimate
        || metadata.consequence != ConsequenceClass::Advisory
    {
        return Err(fail(
            "cli-study-elasticity-decision",
            "elasticity study v1 supports advisory scoping estimates only",
        ));
    }
    if base.versions.as_ref().expect("versions").schema != fs_project::STUDY_FSIM_VERSION {
        return Err(fail(
            "cli-study-elasticity-version",
            "versions.schema must equal the current study schema",
        ));
    }
    let caps = base.capabilities.as_ref().expect("capabilities");
    let required = ["optimization.marquee-topopt", "geometry.sdf", "physics.cutfem"];
    if caps.len() != required.len()
        || !required.iter().all(|needed| caps.iter().any(|value| value == needed))
    {
        return Err(fail(
            "cli-study-elasticity-capability",
            "required capabilities are optimization.marquee-topopt, geometry.sdf and physics.cutfem",
        ));
    }
    if base.units.as_ref().expect("units").storage != "SI" {
        return Err(fail(
            "cli-study-elasticity-units",
            "elasticity study v1 requires `(units :storage \"SI\")`",
        ));
    }
    let domain = base.domain.as_ref().expect("domain");
    if domain.domain_type != "sdf-plate-with-holes"
        || domain.bounds != ([0.0, 0.0], [1.0, 1.0])
    {
        return Err(fail(
            "cli-study-elasticity-domain",
            "elasticity study v1 requires the 1 m x 1 m SDF plate domain",
        ));
    }
    for (index, hole) in domain.initial_holes.iter().enumerate() {
        let edge = hole.center[0]
            .min(1.0 - hole.center[0])
            .min(hole.center[1])
            .min(1.0 - hole.center[1]);
        if !(hole.radius.is_finite() && hole.radius > 0.0 && edge.is_finite() && hole.radius < edge)
        {
            return Err(fail(
                "cli-study-elasticity-hole",
                format!("initial hole {index} must lie strictly inside the unit plate"),
            ));
        }
    }
    let physics = base.physics.as_ref().expect("physics");
    if physics.physics_type != "elasticity-2d" || !(2..=5).contains(&physics.mesh_level) {
        return Err(fail(
            "cli-study-elasticity-physics",
            "supported elasticity envelope is elasticity-2d at mesh level 2..=5",
        ));
    }
    if !(spec.youngs_pa.is_finite() && spec.youngs_pa > 0.0) {
        return Err(fail(
            "cli-study-elasticity-material",
            "youngs-modulus-pa must be finite and positive",
        ));
    }
    if !(spec.poisson.is_finite() && spec.poisson > -1.0 && spec.poisson <= 1.0 / 3.0) {
        return Err(fail(
            "cli-study-elasticity-material",
            "poissons-ratio must lie in (-1, 1/3] for the admitted plane-strain regime",
        ));
    }
    let scenario = base.scenario.as_ref().expect("scenario");
    if scenario.fixed_boundary != "left" || scenario.load_region != "right" {
        return Err(fail(
            "cli-study-elasticity-scenario",
            "v1 elasticity study requires left clamp and right-edge traction",
        ));
    }
    if !(spec.load_traction_pa.is_finite() && spec.load_traction_pa > 0.0)
        || spec.load_direction != [0.0, -1.0]
        || !(spec.load_band[0].is_finite()
            && spec.load_band[1].is_finite()
            && 0.0 <= spec.load_band[0]
            && spec.load_band[0] < spec.load_band[1]
            && spec.load_band[1] <= 1.0
            && (spec.load_band[0] + spec.load_band[1] - 1.0).abs() <= 1e-12)
    {
        return Err(fail(
            "cli-study-elasticity-load",
            "v1 requires finite positive downward traction and a right-edge load band symmetric about y=0.5",
        ));
    }
    let objective = base.objective.as_ref().expect("objective");
    if objective.objective_type != "compliance"
        || objective.sense != "minimize"
        || objective.unit != "J"
    {
        return Err(fail(
            "cli-study-elasticity-objective",
            "v1 requires compliance minimization with objective unit J",
        ));
    }
    if !(spec.move_cells.is_finite()
        && spec.move_cells > 0.0
        && spec.band_cells.is_finite()
        && spec.band_cells > 0.0
        && spec.move_cells <= 0.5 * spec.band_cells
        && spec.ell0.is_finite()
        && spec.ell0 >= 0.0
        && spec.mu_al.is_finite()
        && spec.mu_al > 0.0
        && spec.sobolev_alpha.is_finite()
        && spec.sobolev_alpha >= 0.0
        && spec.hole_radius_cells.is_finite()
        && spec.hole_radius_cells > 0.0)
    {
        return Err(fail(
            "cli-study-elasticity-optimizer",
            "invalid finite level-set move/band/multiplier/smoothing/hole controls",
        ));
    }
    if spec.steps == 0
        || spec.steps > 32
        || spec.steps > spec.max_iterations
        || spec.max_iterations > 32
    {
        return Err(fail(
            "cli-study-elasticity-budget",
            "elasticity v1 admits 1..=32 steps and max-iterations must cover the declared steps",
        ));
    }
    if !spec.wall_s.is_finite()
        || spec.wall_s <= 0.0
        || spec.memory_bytes < 128 * 1024 * 1024
    {
        return Err(fail(
            "cli-study-elasticity-budget",
            "elasticity v1 requires positive wall seconds and at least 128 MiB admitted memory",
        ));
    }
    Ok(())
}

fn initial_phi(spec: &ElasticitySpec) -> GridSdf {
    let n = 1usize << spec.base.physics.as_ref().expect("physics").mesh_level;
    let holes = spec.base.domain.as_ref().expect("domain").initial_holes.clone();
    GridSdf::from_fn(n, &move |x, y| {
        holes
            .iter()
            .map(|hole| hole.radius - (x - hole.center[0]).hypot(y - hole.center[1]))
            .fold(-1.0, f64::max)
    })
}

fn settings(spec: &ElasticitySpec, iterations: usize) -> OptimizeSettings {
    let physics = spec.base.physics.as_ref().expect("physics");
    OptimizeSettings {
        level: physics.mesh_level,
        volfrac: spec.base.constraints.as_ref().expect("constraints").volume_fraction,
        iterations,
        band_cells: spec.band_cells,
        move_cells: spec.move_cells,
        ell0: spec.ell0,
        mu_al: spec.mu_al,
        sobolev_alpha: spec.sobolev_alpha,
        nucleation_period: spec.nucleation_period,
        hole_radius_cells: spec.hole_radius_cells,
        youngs: spec.youngs_pa,
        poisson: spec.poisson,
    }
}

fn fixture(spec: &ElasticitySpec) -> Cantilever {
    Cantilever {
        load: spec.load_traction_pa,
        band: 0.5 * (spec.load_band[1] - spec.load_band[0]),
    }
}

fn trace_hash(rows: &[String]) -> ContentHash {
    let mut bytes = Vec::new();
    for row in rows {
        bytes.extend_from_slice(row.as_bytes());
        bytes.push(b'\n');
    }
    hash_domain(TRACE_DOMAIN, &bytes)
}

fn budget(text: Option<&str>) -> Result<Option<usize>> {
    text.map(|value| {
        value
            .parse::<usize>()
            .ok()
            .filter(|n| *n > 0 && *n <= 32)
            .ok_or_else(|| {
                fail(
                    "cli-study-elasticity-budget-override",
                    "--budget must be an integer in 1..=32",
                )
            })
    })
    .transpose()
}

fn ir(id: ContentHash, ordinal: usize) -> String {
    format!(
        "{{\"driver\":{DRIVER:?},\"study_id\":\"{}\",\"ordinal\":{ordinal},\"units\":\"SI\",\"objective\":\"J\"}}",
        id.to_hex()
    )
}

#[derive(Debug)]
struct Outcome {
    pointer: String,
    receipt: String,
    status: &'static str,
}

fn render(out: Outcome, mode: OutputMode) -> CommandOutput {
    let exit_code = match out.status {
        "completed" => exit::SUCCESS,
        "cancelled" => exit::CANCELLED,
        _ => exit::BUDGET,
    };
    let stdout = match mode {
        OutputMode::Json => format!(
            "{{\"command\":\"study\",\"status\":{:?},\"run_id\":{:?},\"run\":{:?},\"receipt\":{}}}\n",
            out.status, out.pointer, out.pointer, out.receipt
        ),
        OutputMode::Text => format!(
            "command=study\nstatus={}\nrun={}\nauthority=estimated-elasticity\n",
            out.status, out.pointer
        ),
    };
    CommandOutput { exit_code, stdout, stderr: String::new() }
}

fn snapshot(phi: &GridSdf) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for value in phi.nodes() {
        for byte in value.to_bits().to_le_bytes() {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    hash
}

fn design_json(phi: &GridSdf) -> String {
    let mut out = format!(
        "{{\"n\":{},\"snapshot\":\"{:#018x}\",\"phi_bits\":[",
        phi.n(),
        snapshot(phi)
    );
    for (index, value) in phi.nodes().iter().enumerate() {
        if index > 0 { out.push(','); }
        let _ = write!(out, "\"{:016x}\"", value.to_bits());
    }
    out.push_str("]}");
    out
}

fn geometry_svg(phi: &GridSdf) -> String {
    let n = phi.n();
    let mut rects = String::new();
    for j in 0..n {
        for i in 0..n {
            let x = (i as f64 + 0.5) / n as f64;
            let y = (j as f64 + 0.5) / n as f64;
            if phi.value_at([x, y]) <= 0.0 {
                let _ = write!(
                    rects,
                    "<rect x=\"{}\" y=\"{}\" width=\"1\" height=\"1\"/>",
                    i,
                    n - 1 - j
                );
            }
        }
    }
    format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {n} {n}\" width=\"420\" role=\"img\" aria-label=\"Final material layout\"><g fill=\"#111827\">{rects}</g></svg>"
    )
}

#[allow(clippy::too_many_arguments)]
fn persist(
    spec: &ElasticitySpec,
    ledger: &Ledger,
    phi: &GridSdf,
    report: &OptimizeReport,
    status: &'static str,
    wall_s: f64,
    predecessor: Option<ContentHash>,
    evidence: &continuation::Evidence,
) -> Result<Outcome> {
    let continuation = evidence.json();
    let count = report.rows.len();
    let trace = trace_hash(&report.rows);
    let mut rows = format!(
        "{{\"schema\":\"elasticity-study-iterations-v1\",\"study_id\":\"{}\",\"iterations\":[",
        spec.id.to_hex()
    );
    for (index, row) in report.rows.iter().enumerate() {
        if index > 0 { rows.push(','); }
        rows.push_str(row);
    }
    rows.push_str("]}");
    let design = design_json(phi);
    let final_compliance = report
        .compliance
        .last()
        .map_or("null".to_string(), |value| format!("{value:.17e}"));
    let final_volume = report.volume.last().copied();
    let final_snapshot = report.snapshots.last().copied().unwrap_or_else(|| snapshot(phi));
    let summary = format!(
        "{{\"driver\":{DRIVER:?},\"study_id\":\"{}\",\"status\":{status:?},\"iterations_completed\":{count},\"target_iterations\":{},\"final_compliance_j\":{final_compliance},\"final_material_area_m2\":{},\"snapshot\":\"{final_snapshot:#018x}\",\"trace_hash\":\"{}\",\"authority\":\"Estimated\",\"no_claim\":{}}}",
        spec.id.to_hex(),
        spec.steps,
        final_volume.map_or("null".to_string(), |value| format!("{value:.17e}")),
        trace.to_hex(),
        quoted(NO_CLAIM)
    );
    let mut table = String::new();
    for (index, ((compliance, volume), snap)) in report
        .compliance
        .iter()
        .zip(&report.volume)
        .zip(&report.snapshots)
        .enumerate()
    {
        let _ = writeln!(
            table,
            "<tr><td>{}</td><td>{:.8e}</td><td>{:.8e}</td><td>{snap:#018x}</td></tr>",
            index + 1,
            compliance,
            volume
        );
    }
    let html = format!(
        "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Elasticity topology study</title><body><h1>Free-boundary 2-D elasticity topology study</h1><p>Status: {status}. {count}/{} iterations. Estimated.</p><p>{NO_CLAIM}</p><p>Final discrete compliance: {final_compliance} J.</p>{}<table><tr><th>Iteration</th><th>Compliance J</th><th>Material area m²</th><th>Snapshot</th></tr>{table}</table><p>Trace: {}</p></body></html>",
        spec.steps,
        geometry_svg(phi),
        trace.to_hex()
    );
    let package = EvidencePackage::new(Provenance::new(
        format!("fs-cli/{}+{DRIVER}", env!("CARGO_PKG_VERSION")),
        hash_bytes(include_bytes!("../../../constellation.lock")).to_hex(),
    ));
    let package_json = package
        .to_json()
        .map_err(|error| fail("cli-study-elasticity-package", error.to_string()))?;
    if !fs_checker::check(&package).passed() {
        return Err(fail(
            "cli-study-elasticity-package",
            "checker refused the structural evidence package",
        ));
    }
    let versions = format!(
        "{{\"driver\":{DRIVER:?},\"crate\":{:?},\"fs_topols\":{:?},\"constellation_lock\":\"{}\"}}",
        env!("CARGO_PKG_VERSION"),
        fs_topols::VERSION,
        hash_bytes(include_bytes!("../../../constellation.lock")).to_hex()
    );
    let budget_json = format!(
        "{{\"wall_s\":{},\"memory_bytes\":{},\"max_iterations\":{},\"consumed_wall_s\":{wall_s}}}",
        spec.wall_s, spec.memory_bytes, spec.max_iterations
    );
    let seed = spec.base.seeds.as_ref().expect("seed").root.to_le_bytes();
    ledger.begin()?;
    let result = (|| -> Result<Outcome> {
        let op = ledger.begin_op(
            Some(spec.id.as_bytes()),
            &ir(spec.id, count),
            &FiveExplicits {
                seed: &seed,
                versions: &versions,
                budget: &budget_json,
                capability: "{\"ops\":[\"optimization.marquee-topopt\",\"geometry.sdf\",\"physics.cutfem\"]}",
            },
            0,
        )?;
        if let Some(previous) = predecessor {
            ledger.link(op, &previous, EdgeRole::In)?;
        }
        let source = ledger.put_artifact("study-source", spec.canonical.as_bytes(), None)?;
        ledger.link(op, &source.hash, EdgeRole::In)?;
        let mut refs = String::new();
        for (name, kind, bytes) in [
            ("iterations", "study-iterations", rows.as_bytes()),
            ("design", "study-design", design.as_bytes()),
            ("report_html", "study-report-html", html.as_bytes()),
            ("report_json", "study-report-json", summary.as_bytes()),
            ("package", "study-package", package_json.as_bytes()),
        ] {
            let artifact = ledger.put_artifact(kind, bytes, None)?;
            ledger.link(op, &artifact.hash, EdgeRole::Out)?;
            let _ = write!(refs, ",{name:?}:\"{}\"", artifact.hash.to_hex());
        }
        let previous_json = predecessor.map_or("null".to_string(), |hash| quoted(&hash.to_hex()));
        let receipt = format!(
            "{{\"schema\":{STUDY_RUN_RECEIPT_SCHEMA:?},\"driver\":{DRIVER:?},\"study_id\":\"{}\",\"status\":{status:?},\"source\":\"{}\",\"iterations_completed\":{count},\"target_iterations\":{},\"trace_hash\":\"{}\",\"consumed_wall_s\":{wall_s},\"predecessor\":{previous_json},\"continuation\":{continuation}{refs}}}",
            spec.id.to_hex(),
            source.hash.to_hex(),
            spec.steps,
            trace.to_hex()
        );
        let artifact = ledger.put_artifact(RECEIPT_KIND, receipt.as_bytes(), None)?;
        ledger.link(op, &artifact.hash, EdgeRole::Out)?;
        if ledger.artifact_output_seal(&artifact.hash)?.is_none() {
            ledger.seal_artifact_output(&artifact.hash, op)?;
        }
        ledger.finish_op(op, OpOutcome::Ok, None, 1)?;
        Ok(Outcome {
            pointer: format!("study-{}", artifact.hash.to_hex()),
            receipt,
            status,
        })
    })();
    match result {
        Ok(out) => {
            if let Err(error) = ledger.commit() {
                return match ledger.rollback() {
                    Ok(()) => Err(error.into()),
                    Err(rollback) => Err(fail(
                        "cli-study-elasticity-ledger",
                        format!("commit failed: {error}; rollback also failed: {rollback}"),
                    )),
                };
            }
            Ok(out)
        }
        Err(error) => match ledger.rollback() {
            Ok(()) => Err(error),
            Err(rollback) => Err(fail(
                "cli-study-elasticity-ledger",
                format!("{}; rollback also failed: {rollback}", error.message),
            )),
        },
    }
}

#[derive(Debug)]
struct Loaded {
    hash: ContentHash,
    value: JsonValue,
    bytes: String,
}

fn artifact(ledger: &Ledger, hash: ContentHash, kind: &str) -> Result<Vec<u8>> {
    let info = ledger
        .artifact_info(&hash)?
        .ok_or_else(|| fail("cli-study-elasticity-artifact", "missing study artifact"))?;
    if info.kind != kind {
        return Err(fail(
            "cli-study-elasticity-artifact",
            format!("expected {kind}, found {}", info.kind),
        ));
    }
    ledger
        .get_artifact_bounded(&hash, MAX_ARTIFACT_BYTES)?
        .ok_or_else(|| fail("cli-study-elasticity-artifact", "missing study bytes"))
}

fn linked(ledger: &Ledger, value: &JsonValue, key: &str, kind: &str) -> Result<Vec<u8>> {
    let hash = value
        .str_field(key)
        .and_then(ContentHash::from_hex)
        .ok_or_else(|| fail("cli-study-elasticity-receipt", format!("missing {key} hash")))?;
    artifact(ledger, hash, kind)
}

fn integer(value: &JsonValue, key: &str) -> Result<usize> {
    value
        .get(key)
        .and_then(JsonValue::number_raw)
        .and_then(|text| text.parse().ok())
        .ok_or_else(|| fail("cli-study-elasticity-receipt", format!("invalid {key}")))
}

fn load(ledger: &Ledger, pointer: &str) -> Result<Loaded> {
    let hash = pointer
        .strip_prefix("study-")
        .and_then(ContentHash::from_hex)
        .ok_or_else(|| {
            fail(
                "cli-study-elasticity-run-id",
                "expected study- followed by a 64-digit receipt hash",
            )
        })?;
    let bytes = artifact(ledger, hash, RECEIPT_KIND)?;
    let text = std::str::from_utf8(&bytes)
        .map_err(|error| fail("cli-study-elasticity-receipt", error.to_string()))?
        .to_string();
    let value = JsonValue::parse(&text)
        .map_err(|error| fail("cli-study-elasticity-receipt", error.to_string()))?;
    if value.str_field("schema") != Some(STUDY_RUN_RECEIPT_SCHEMA)
        || value.str_field("driver") != Some(DRIVER)
    {
        return Err(fail(
            "cli-study-elasticity-receipt",
            "receipt is not an elasticity study",
        ));
    }
    let study_id = value
        .str_field("study_id")
        .and_then(ContentHash::from_hex)
        .ok_or_else(|| fail("cli-study-elasticity-receipt", "missing study identity"))?;
    let ordinal = integer(&value, "iterations_completed")?;
    let producer = ledger
        .artifact_output_seal(&hash)?
        .ok_or_else(|| fail("cli-study-elasticity-receipt", "unsealed study receipt"))?;
    let op = ledger
        .op(producer)?
        .ok_or_else(|| fail("cli-study-elasticity-receipt", "missing producing operation"))?;
    if op.session.as_deref() != Some(study_id.as_bytes().as_slice())
        || op.ir != ir(study_id, ordinal)
        || op.outcome.as_deref() != Some("ok")
        || !ledger.edge_exists(producer, &hash, EdgeRole::Out)?
    {
        return Err(fail(
            "cli-study-elasticity-receipt",
            "receipt is not bound to its completed elasticity operation",
        ));
    }
    for (key, kind, role) in [
        ("source", "study-source", EdgeRole::In),
        ("iterations", "study-iterations", EdgeRole::Out),
        ("design", "study-design", EdgeRole::Out),
        ("report_html", "study-report-html", EdgeRole::Out),
        ("report_json", "study-report-json", EdgeRole::Out),
        ("package", "study-package", EdgeRole::Out),
    ] {
        let linked_hash = value
            .str_field(key)
            .and_then(ContentHash::from_hex)
            .ok_or_else(|| fail("cli-study-elasticity-receipt", format!("missing {key}")))?;
        if !ledger.edge_exists(producer, &linked_hash, role)? {
            return Err(fail(
                "cli-study-elasticity-receipt",
                format!("missing {key} lineage"),
            ));
        }
        artifact(ledger, linked_hash, kind)?;
    }
    Ok(Loaded { hash, value, bytes: text })
}

fn run_prefix(spec: &ElasticitySpec, count: usize) -> Result<(GridSdf, OptimizeReport)> {
    let mut phi = initial_phi(spec);
    if count == 0 {
        return Ok((phi, OptimizeReport::default()));
    }
    let report = optimize_compliance(&mut phi, fixture(spec), settings(spec, count))
        .map_err(|error| fail("cli-study-elasticity-solve", format!("{error:?}")))?;
    if report.rows.len() != count || report.snapshots.len() != count {
        return Err(fail(
            "cli-study-elasticity-solve",
            "optimizer returned an incomplete trajectory",
        ));
    }
    Ok((phi, report))
}

pub(crate) fn study_path(
    path: &Path,
    ledger_path: &Path,
    override_text: Option<&str>,
    mode: OutputMode,
) -> CommandOutput {
    let result = (|| {
        if path.extension().and_then(|value| value.to_str()) != Some("fsim") {
            return Err(fail(
                "cli-study-elasticity-format",
                "elasticity studies currently require canonical .fsim s-expression input",
            ));
        }
        let cap = budget(override_text)?;
        let mut bytes = Vec::new();
        let file = std::fs::File::open(path).map_err(|error| Failure {
            code: "cli-study-elasticity-read",
            message: error.to_string(),
            exit: exit::INPUT,
        })?;
        file.take(MAX_PROJECT_BYTES + 1)
            .read_to_end(&mut bytes)
            .map_err(|error| Failure {
                code: "cli-study-elasticity-read",
                message: error.to_string(),
                exit: exit::INPUT,
            })?;
        if bytes.len() as u64 > MAX_PROJECT_BYTES {
            return Err(Failure {
                code: "cli-study-elasticity-size",
                message: "study exceeds 16 MiB".into(),
                exit: exit::INPUT,
            });
        }
        let source = std::str::from_utf8(&bytes).map_err(|error| Failure {
            code: "cli-study-elasticity-utf8",
            message: error.to_string(),
            exit: exit::INPUT,
        })?;
        let spec = parse(source)?;
        let ledger = Ledger::open(
            ledger_path
                .to_str()
                .ok_or_else(|| fail("cli-study-elasticity-ledger-path", "ledger path is not UTF-8"))?,
        )?;
        drive(&spec, &ledger, cap, &CancelGate::new(), None)
    })();
    match result {
        Ok(out) => render(out, mode),
        Err(error) => failure_output("study", mode, error),
    }
}

pub(crate) fn resume_path(
    pointer: &str,
    path: &Path,
    override_text: Option<&str>,
    mode: OutputMode,
) -> CommandOutput {
    let result = (|| {
        let cap = budget(override_text)?;
        if !path.is_file() {
            return Err(fail(
                "cli-study-elasticity-ledger-missing",
                "resume requires an existing ledger",
            ));
        }
        let ledger = Ledger::open(
            path.to_str()
                .ok_or_else(|| fail("cli-study-elasticity-ledger-path", "ledger path is not UTF-8"))?,
        )?;
        let old = load(&ledger, pointer)?;
        let source = linked(&ledger, &old.value, "source", "study-source")?;
        let source = std::str::from_utf8(&source)
            .map_err(|error| fail("cli-study-elasticity-receipt", error.to_string()))?;
        let spec = parse(source)?;
        drive(&spec, &ledger, cap, &CancelGate::new(), Some(&old))
    })();
    match result {
        Ok(out) => render(out, mode),
        Err(error) => failure_output("study", mode, error),
    }
}

pub(crate) fn owns_run(pointer: &str, path: &Path) -> bool {
    if !path.is_file() {
        return false;
    }
    let Ok(ledger) = Ledger::open(path.to_str().unwrap_or("")) else {
        return false;
    };
    let Some(hash) = pointer.strip_prefix("study-").and_then(ContentHash::from_hex) else {
        return false;
    };
    let Ok(Some(bytes)) = ledger.get_artifact_bounded(&hash, MAX_ARTIFACT_BYTES) else {
        return false;
    };
    let Ok(text) = std::str::from_utf8(&bytes) else {
        return false;
    };
    let Ok(value) = JsonValue::parse(text) else {
        return false;
    };
    value.str_field("schema") == Some(STUDY_RUN_RECEIPT_SCHEMA)
        && value.str_field("driver") == Some(DRIVER)
}

pub(crate) fn looks_like(path: &Path) -> bool {
    if path.extension().and_then(|value| value.to_str()) != Some("fsim") {
        return false;
    }
    let Ok(file) = std::fs::File::open(path) else {
        return false;
    };
    let mut bytes = Vec::new();
    if file
        .take(MAX_PROJECT_BYTES + 1)
        .read_to_end(&mut bytes)
        .is_err()
        || bytes.len() as u64 > MAX_PROJECT_BYTES
    {
        return false;
    }
    std::str::from_utf8(&bytes).is_ok_and(|text| text.contains("elasticity-2d"))
}

pub(crate) fn export(
    command: &'static str,
    pointer: &str,
    path: Option<&Path>,
    mode: OutputMode,
) -> CommandOutput {
    let result = (|| -> Result<String> {
        let path = path.filter(|value| value.is_file()).ok_or_else(|| {
            fail(
                "cli-study-elasticity-ledger-missing",
                "study export requires an existing ledger operand",
            )
        })?;
        let ledger = Ledger::open(
            path.to_str()
                .ok_or_else(|| fail("cli-study-elasticity-ledger-path", "ledger path is not UTF-8"))?,
        )?;
        let loaded = load(&ledger, pointer)?;
        let dir = path.parent().unwrap_or_else(|| Path::new("."));
        let fields: &[(&str, &str, &str)] = if command == "package" {
            &[("package", "study-package", "fspkg")]
        } else {
            &[
                ("report_html", "study-report-html", "html"),
                ("report_json", "study-report-json", "json"),
            ]
        };
        let mut paths = String::new();
        for &(key, kind, extension) in fields {
            let bytes = linked(&ledger, &loaded.value, key, kind)?;
            if key == "package" {
                let text = std::str::from_utf8(&bytes)
                    .map_err(|error| fail("cli-study-elasticity-package", error.to_string()))?;
                let package = EvidencePackage::from_json(text)
                    .map_err(|error| fail("cli-study-elasticity-package", error.to_string()))?;
                if !fs_checker::check(&package).passed() {
                    return Err(fail(
                        "cli-study-elasticity-package",
                        "retained package failed structural verification",
                    ));
                }
            }
            let dest = dir.join(format!("{pointer}.{extension}"));
            crate::report::write_retained(&dest, &bytes)
                .map_err(|error| fail("cli-study-elasticity-export", error))?;
            let _ = write!(paths, ",{key:?}:{}", quoted(&dest.to_string_lossy()));
        }
        Ok(format!(
            "{{\"command\":{command:?},\"status\":\"ok\",\"run\":{pointer:?},\"study_status\":{},\"authority\":\"projection-of-retained-estimates\",\"verification\":\"sealed-evidence\"{paths}}}\n",
            quoted(loaded.value.str_field("status").unwrap_or("unknown"))
        ))
    })();
    match result {
        Ok(mut stdout) => {
            if matches!(mode, OutputMode::Text) {
                stdout = format!(
                    "command={command}\nstatus=ok\nrun={pointer}\nauthority=projection-of-retained-estimates\n"
                );
            }
            CommandOutput { exit_code: exit::SUCCESS, stdout, stderr: String::new() }
        }
        Err(error) => failure_output(command, mode, error),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIXTURE: &str = include_str!("../../../examples/marquee/bracket-2d.fsim");

    #[test]
    fn canonical_fixture_is_explicit_and_admitted() {
        let spec = parse(FIXTURE).expect("canonical elasticity fixture");
        assert_eq!(spec.steps, 8);
        assert_eq!(spec.load_direction, [0.0, -1.0]);
        assert_eq!(spec.base.physics.as_ref().unwrap().mesh_level, 4);
    }

    #[test]
    fn missing_material_or_optimizer_field_refuses() {
        for field in [
            "    :youngs-modulus-pa 70000000000.0\n",
            "    :move-cells 0.35\n",
            "    :load-band (0.375 0.625))\n",
        ] {
            let bad = FIXTURE.replace(field, "");
            assert!(parse(&bad).is_err(), "must refuse missing {field:?}");
        }
    }

    #[test]
    fn deterministic_prefix_replay_extends_without_changing_retained_rows() {
        let mut spec = parse(FIXTURE).expect("fixture");
        spec.base.physics.as_mut().unwrap().mesh_level = 2;
        spec.steps = 2;
        spec.max_iterations = 2;
        let first = run_prefix(&spec, 1).expect("first prefix");
        let second = run_prefix(&spec, 2).expect("second prefix");
        assert_eq!(first.1.rows[0], second.1.rows[0]);
        assert_eq!(first.1.snapshots[0], second.1.snapshots[0]);
        assert_eq!(trace_hash(&first.1.rows), trace_hash(&second.1.rows[..1]));
    }
}
