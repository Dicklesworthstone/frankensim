//! Assembly `.fsim` (s-expression & isomorphic `fs-ir` JSON) parser, Five
//! Explicits admission, and `CardPackSet -> fs_system::CardStore` bridge.

use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{ContentHash, hash_bytes};
use fs_evidence::uncertainty::RequirementRelation;
use fs_ir::ast::{CountUnit, Node, NodeKind};
use fs_matdb::{NormalizedInterfacePack, NormalizedMaterialCardPack};
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftAssemblyInput, DraftEnvironment, DraftGeometrySource, DraftInitialState,
    DraftInterfaceDecl, DraftInterfaceKind, DraftPartDecl, DraftPartKind, FidelityErrorBudget,
    FidelityPolicy, FidelityRung, JointDecl, JointKind, MicrophonePattern, PhysicsDomain,
    PrimitiveShape, ProbeDecl, ProbeKind, SourceDecl, SourceTarget, StepPolicy, TimeHistory,
    TimePlan, TimingEdge,
};
use fs_scenario::assembly::{MaterialRef, ParameterKey, ReferenceState};
use fs_scenario::{EntityCatalog, EntityDeclaration, EntityId, FrameId, InterfacePair, Scenario};
use fs_system::CardStore;

use crate::cards::{CardPackKind, CardPackSet};
use crate::solve::SolveRefusal;
use crate::stage_driver::staged_refusal;

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const FREQUENCY: Dims = Dims([0, 0, -1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const AREA: Dims = Dims([2, 0, 0, 0, 0, 0]);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ProbeStat {
    Max,
    Min,
    Final,
}

impl ProbeStat {
    pub(crate) const fn as_str(self) -> &'static str {
        match self {
            Self::Max => "max",
            Self::Min => "min",
            Self::Final => "final",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AssemblyRequirement {
    pub(crate) id: String,
    pub(crate) probe: String,
    pub(crate) stat: ProbeStat,
    pub(crate) relation: RequirementRelation,
    pub(crate) limit: QtyAny,
}

/// Admitted assembly document, Five Explicits, and resolved card store.
pub(crate) struct AdmittedAssembly {
    pub(crate) canonical_source: Vec<u8>,
    pub(crate) source_hash: ContentHash,
    pub(crate) name: String,
    pub(crate) seed: u64,
    pub(crate) constellation: String,
    pub(crate) workspace: String,
    pub(crate) wall_budget_s: f64,
    pub(crate) mem_budget_bytes: u64,
    pub(crate) cores: u64,
    pub(crate) granted_ops: Vec<String>,
    pub(crate) draft: DraftAssembly,
    pub(crate) card_store: CardStore,
    pub(crate) entity_names: BTreeMap<EntityId, String>,
    pub(crate) harmonic_frequencies_hz: Vec<f64>,
    pub(crate) modal_max_modes: usize,
    pub(crate) requirements: Vec<AssemblyRequirement>,
}

fn err(what: impl Into<String>) -> SolveRefusal {
    staged_refusal(
        "cli-assembly-admit",
        "admit",
        what,
        "fix the assembly declaration",
    )
}

fn as_list<'a>(n: &'a Node, w: &str) -> Result<&'a [Node], SolveRefusal> {
    match &n.kind {
        NodeKind::List(v) if !v.is_empty() => Ok(v),
        _ => Err(err(format!("expected non-empty list for {w}"))),
    }
}
fn as_sym<'a>(n: &'a Node, w: &str) -> Result<&'a str, SolveRefusal> {
    match &n.kind {
        NodeKind::Symbol(s) => Ok(s.as_str()),
        _ => Err(err(format!("expected symbol for {w}"))),
    }
}
fn as_str<'a>(n: &'a Node, w: &str) -> Result<&'a str, SolveRefusal> {
    match &n.kind {
        NodeKind::Str(s) => Ok(s.as_str()),
        _ => Err(err(format!("expected string for {w}"))),
    }
}
fn as_f64(n: &Node, w: &str) -> Result<f64, SolveRefusal> {
    match &n.kind {
        NodeKind::Float(v) if v.is_finite() => Ok(*v),
        #[allow(clippy::cast_precision_loss)]
        NodeKind::Int(v) => Ok(*v as f64),
        _ => Err(err(format!("expected finite number for {w}"))),
    }
}
fn as_qty(n: &Node, w: &str) -> Result<QtyAny, SolveRefusal> {
    match &n.kind {
        NodeKind::Qty { value, dims, .. } if value.is_finite() => Ok(QtyAny::new(*value, *dims)),
        _ => Err(err(format!("expected dimensioned quantity for {w}"))),
    }
}
fn as_qty_d(n: &Node, d: Dims, w: &str) -> Result<QtyAny, SolveRefusal> {
    let q = as_qty(n, w)?;
    if q.dims != d {
        return Err(err(format!("{w} has wrong SI dimensions")));
    }
    Ok(q)
}
fn as_bytes(n: &Node, w: &str) -> Result<u64, SolveRefusal> {
    match &n.kind {
        NodeKind::Count { value, unit } => value
            .integral_bytes(*unit)
            .filter(|&b| b > 0)
            .ok_or_else(|| err(format!("{w} must be a positive byte count"))),
        _ => Err(err(format!("expected byte count for {w}"))),
    }
}
fn as_vec3(n: &Node, w: &str) -> Result<[f64; 3], SolveRefusal> {
    let v = as_list(n, w)?;
    if v.len() != 3 {
        return Err(err(format!("{w} must have 3 elements")));
    }
    Ok([as_f64(&v[0], w)?, as_f64(&v[1], w)?, as_f64(&v[2], w)?])
}
fn as_qty3(n: &Node, w: &str) -> Result<[QtyAny; 3], SolveRefusal> {
    let v = as_list(n, w)?;
    if v.len() != 3 {
        return Err(err(format!("{w} must have 3 elements")));
    }
    Ok([as_qty(&v[0], w)?, as_qty(&v[1], w)?, as_qty(&v[2], w)?])
}
fn as_quat4(n: &Node, w: &str) -> Result<[f64; 4], SolveRefusal> {
    let v = as_list(n, w)?;
    if v.len() != 4 {
        return Err(err(format!("{w} must have 4 elements")));
    }
    Ok([
        as_f64(&v[0], w)?,
        as_f64(&v[1], w)?,
        as_f64(&v[2], w)?,
        as_f64(&v[3], w)?,
    ])
}
fn parse_initial(n: &Node, part_scope: &str) -> Result<DraftInitialState, SolveRefusal> {
    let scope = format!("{part_scope} :initial");
    let kw = kwargs(
        as_list(n, &scope)?,
        &[
            "position",
            "orientation",
            "linear-velocity",
            "angular-velocity",
            "temperature",
        ],
        &scope,
    )?;
    Ok(DraftInitialState {
        position: as_qty3(req(&kw, "position", &scope)?, "initial.position")?,
        orientation: as_quat4(req(&kw, "orientation", &scope)?, "initial.orientation")?,
        linear_velocity: as_qty3(
            req(&kw, "linear-velocity", &scope)?,
            "initial.linear-velocity",
        )?,
        angular_velocity: as_qty3(
            req(&kw, "angular-velocity", &scope)?,
            "initial.angular-velocity",
        )?,
        temperature: kw
            .get("temperature")
            .map(|t| as_qty(t, "initial.temperature"))
            .transpose()?,
    })
}

fn kwargs<'a>(
    nodes: &'a [Node],
    ok: &[&str],
    scope: &str,
) -> Result<BTreeMap<&'a str, &'a Node>, SolveRefusal> {
    if !nodes.len().is_multiple_of(2) {
        return Err(err(format!("{scope} has unpaired keyword")));
    }
    let mut map = BTreeMap::new();
    let (pairs, _) = nodes.as_chunks::<2>();
    for pair in pairs {
        let NodeKind::Keyword(k) = &pair[0].kind else {
            return Err(err(format!("{scope} expected :keyword")));
        };
        if !ok.contains(&k.as_str()) {
            return Err(err(format!("{scope} rejects unknown keyword `:{k}`")));
        }
        if map.insert(k.as_str(), &pair[1]).is_some() {
            return Err(err(format!("{scope} repeats `:{k}`")));
        }
    }
    Ok(map)
}
fn req<'a>(kw: &BTreeMap<&str, &'a Node>, k: &str, scope: &str) -> Result<&'a Node, SolveRefusal> {
    kw.get(k)
        .copied()
        .ok_or_else(|| err(format!("{scope} missing `:{k}`")))
}

fn parse_geom(n: &Node) -> Result<DraftGeometrySource, SolveRefusal> {
    let v = as_list(n, "geometry")?;
    Ok(DraftGeometrySource::Primitive(
        match (as_sym(&v[0], "geometry")?, v.len()) {
            ("box", 4) => PrimitiveShape::Box {
                lx: as_qty(&v[1], "lx")?,
                ly: as_qty(&v[2], "ly")?,
                lz: as_qty(&v[3], "lz")?,
            },
            ("cylinder", 3) => PrimitiveShape::Cylinder {
                radius: as_qty(&v[1], "r")?,
                height: as_qty(&v[2], "h")?,
            },
            ("sphere", 2) => PrimitiveShape::Sphere {
                radius: as_qty(&v[1], "r")?,
            },
            ("tube", 4) => PrimitiveShape::Tube {
                inner_radius: as_qty(&v[1], "ri")?,
                outer_radius: as_qty(&v[2], "ro")?,
                length: as_qty(&v[3], "l")?,
            },
            (o, _) => return Err(err(format!("unsupported geometry `{o}`"))),
        },
    ))
}

fn parse_history(n: &Node) -> Result<TimeHistory, SolveRefusal> {
    let v = as_list(n, "history")?;
    match as_sym(&v[0], "history")? {
        "constant" if v.len() == 2 => Ok(TimeHistory::Constant(as_qty(&v[1], "constant")?)),
        "harmonic" => {
            let kw = kwargs(
                &v[1..],
                &["amplitude", "frequency", "phase", "offset"],
                "harmonic",
            )?;
            let phase = kw.get("phase").map_or(Ok(0.0), |x| as_f64(x, "phase"))?;
            Ok(TimeHistory::Harmonic {
                amplitude: as_qty(req(&kw, "amplitude", "harmonic")?, "amplitude")?,
                frequency: as_qty_d(req(&kw, "frequency", "harmonic")?, FREQUENCY, "frequency")?,
                phase: QtyAny::dimensionless(phase),
                offset: as_qty(req(&kw, "offset", "harmonic")?, "offset")?,
            })
        }
        "step" => {
            let kw = kwargs(&v[1..], &["t", "before", "after"], "step")?;
            Ok(TimeHistory::Step {
                t_switch: as_qty_d(req(&kw, "t", "step")?, TIME, "t")?,
                before: as_qty(req(&kw, "before", "step")?, "before")?,
                after: as_qty(req(&kw, "after", "step")?, "after")?,
            })
        }
        o => Err(err(format!("unsupported history `{o}`"))),
    }
}

fn resolve_card(
    packs: &CardPackSet,
    kind: CardPackKind,
    sel: &str,
    scope: &str,
) -> Result<ContentHash, SolveRefusal> {
    let slice = match kind {
        CardPackKind::Material => packs.materials(),
        CardPackKind::Interface => packs.interfaces(),
    };
    let hits: Vec<ContentHash> = slice
        .iter()
        .filter(|p| {
            p.card().to_hex() == sel
                || p.identity() == sel
                || p.identity().starts_with(&format!("{sel}/"))
        })
        .map(crate::cards::AdmittedCardPack::card)
        .collect();
    match hits.as_slice() {
        [only] => Ok(*only),
        [] => Err(staged_refusal(
            "cli-assembly-card-missing",
            "admit",
            format!("{scope} missing {} card `{sel}`", kind.label()),
            format!("supply the card pack via {}", kind.flag()),
        )),
        _ => Err(staged_refusal(
            "cli-assembly-card-missing",
            "admit",
            format!("{scope} card selector `{sel}` is ambiguous"),
            "specify full identity or 64-hex hash",
        )),
    }
}

fn build_card_store(
    packs: &CardPackSet,
    families: &BTreeMap<ContentHash, String>,
) -> Result<CardStore, SolveRefusal> {
    let mut store = CardStore::new();
    for p in packs.materials() {
        let dec = NormalizedMaterialCardPack::from_bytes(p.bytes())
            .map_err(|e| err(format!("material pack decode failed: {e}")))?;
        store.insert_material(dec.card().clone());
    }
    for p in packs.interfaces() {
        let dec = NormalizedInterfacePack::from_bytes(p.bytes())
            .map_err(|e| err(format!("interface pack decode failed: {e}")))?;
        let fam = families
            .get(&p.card())
            .cloned()
            .unwrap_or_else(|| "thermal-conductance".into());
        store
            .insert_interface(fam, dec.card().clone())
            .map_err(|e| err(format!("interface store error: {e}")))?;
    }
    Ok(store)
}

fn parse_probe(
    items: &[Node],
    parts: &BTreeMap<String, EntityId>,
) -> Result<ProbeDecl, SolveRefusal> {
    let name = as_str(
        items.get(1).ok_or_else(|| err("probe missing name"))?,
        "probe.name",
    )?
    .to_string();
    let scope = format!("probe `{name}`");
    let kw = kwargs(
        &items[2..],
        &[
            "kind",
            "target",
            "sample-rate",
            "time-constant",
            "threshold",
            "edge",
            "quantity",
            "dims",
            "direction",
        ],
        &scope,
    )?;
    let sample_rate = as_qty_d(req(&kw, "sample-rate", &scope)?, FREQUENCY, "sample-rate")?;
    let zero = [QtyAny::new(0.0, Dims([1, 0, 0, 0, 0, 0])); 3];
    let target = || -> Result<EntityId, SolveRefusal> {
        let tn = as_str(req(&kw, "target", &scope)?, "target")?;
        parts
            .get(tn)
            .copied()
            .ok_or_else(|| err(format!("{scope} unknown target `{tn}`")))
    };
    let kind = match as_sym(req(&kw, "kind", &scope)?, "probe.kind")? {
        "thermocouple" => ProbeKind::Thermocouple {
            location: target()?,
            position: zero,
            time_constant: kw
                .get("time-constant")
                .map(|n| as_qty_d(n, TIME, "time-constant"))
                .transpose()?,
        },
        "timing" => ProbeKind::Timing {
            target: target()?,
            threshold: as_qty(req(&kw, "threshold", &scope)?, "threshold")?,
            edge: match kw.get("edge").map(|n| as_sym(n, "edge")).transpose()? {
                None | Some("rising") => TimingEdge::Rising,
                Some("falling") => TimingEdge::Falling,
                Some("both") => TimingEdge::Both,
                Some(o) => return Err(err(format!("unknown timing edge `{o}`"))),
            },
        },
        "field-quantity" => {
            let qn = as_str(req(&kw, "quantity", &scope)?, "quantity")?.to_string();
            let dn = as_list(req(&kw, "dims", &scope)?, "dims")?;
            if dn.len() != 6 {
                return Err(err("probe :dims requires 6 integers"));
            }
            let mut exps = [0i8; 6];
            for (slot, node) in exps.iter_mut().zip(dn) {
                let NodeKind::Int(i) = node.kind else {
                    return Err(err("dims must be ints"));
                };
                *slot = i8::try_from(i).map_err(|_| err("dims exponent out of i8 range"))?;
            }
            ProbeKind::FieldQuantity {
                target: target()?,
                position: zero,
                quantity: ParameterKey::new(qn, Dims(exps)),
            }
        }
        "strain-gauge" => ProbeKind::StrainGauge {
            target: target()?,
            position: zero,
            gauge_direction: kw
                .get("direction")
                .map_or(Ok([1.0, 0.0, 0.0]), |n| as_vec3(n, "direction"))?,
        },
        "microphone" => ProbeKind::Microphone {
            position: zero,
            aim_direction: kw
                .get("direction")
                .map_or(Ok([1.0, 0.0, 0.0]), |n| as_vec3(n, "direction"))?,
            pattern: MicrophonePattern::Omnidirectional,
        },
        o => return Err(err(format!("unsupported probe kind `{o}`"))),
    };
    Ok(ProbeDecl {
        name,
        frame: FrameId(0),
        sample_rate,
        kind,
    })
}

#[allow(clippy::too_many_lines)]
pub(crate) fn admit_assembly_source(
    source: &[u8],
    packs: &CardPackSet,
) -> Result<AdmittedAssembly, SolveRefusal> {
    let text = std::str::from_utf8(source).map_err(|e| err(format!("non-UTF-8 source: {e}")))?;
    let ast = if text.trim_start().starts_with(['{', '[']) {
        fs_ir::json::parse(text)
    } else {
        fs_ir::sexpr::parse(text)
    }
    .map_err(|e| err(format!("assembly parse failed: {e}")))?;
    let canonical_source = fs_ir::sexpr::print_checked(&ast)
        .map_err(|e| err(format!("canonicalization failed: {e}")))?
        .into_bytes();
    let source_hash = hash_bytes(&canonical_source);

    let root_items = as_list(&ast, "assembly")?;
    if as_sym(&root_items[0], "root")? != "assembly" || root_items.len() < 3 {
        return Err(err("top-level form must be (assembly \"name\" ...)"));
    }
    let name = as_str(&root_items[1], "assembly.name")?.to_string();
    let (mut version, mut units_ok, mut seed, mut versions, mut budget, mut capability) =
        (None, false, None, None, None, None);
    let (mut ref_state, mut env, mut time_plan) = (None, None, None);
    let (mut harm_hz, mut modal_max) = (Vec::new(), 8usize);

    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly(&name))
        .map_err(|e| err(format!("catalog error: {e}")))?;
    let mut entity_names = BTreeMap::from([(root, name.clone())]);
    let (mut by_name, mut families) = (BTreeMap::new(), BTreeMap::new());
    let (mut parts, mut joints, mut interfaces, mut sources, mut probes, mut requirements) = (
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    );
    let mut req_ids = BTreeSet::new();

    for clause in &root_items[2..] {
        let items = as_list(clause, "clause")?;
        match as_sym(&items[0], "clause")? {
            "version" if items.len() == 2 => match items[1].kind {
                NodeKind::Int(v) => version = Some(v),
                _ => return Err(err("version must be integer 1")),
            },
            "units" if items.len() == 2 => {
                if as_str(&items[1], "units")? != "SI" {
                    return Err(err("only SI units are admitted"));
                }
                units_ok = true;
            }
            "seed" if items.len() == 2 => match items[1].kind {
                NodeKind::Seed(s) => seed = Some(s),
                _ => return Err(err("seed must be hex literal 0x...")),
            },
            "versions" => {
                let (mut c, mut w) = (None, None);
                for sub in &items[1..] {
                    let s = as_list(sub, "versions")?;
                    match (as_sym(&s[0], "versions")?, s.len()) {
                        ("constellation", 2) => {
                            c = Some(as_str(&s[1], "constellation")?.to_string());
                        }
                        ("workspace", 2) => {
                            w = Some(as_str(&s[1], "workspace")?.to_string());
                        }
                        (o, _) => return Err(err(format!("unknown versions field `{o}`"))),
                    }
                }
                versions = Some((
                    c.ok_or_else(|| err("missing constellation"))?,
                    w.ok_or_else(|| err("missing workspace"))?,
                ));
            }
            "budget" => {
                let (mut w, mut m) = (None, None);
                for sub in &items[1..] {
                    let s = as_list(sub, "budget")?;
                    match (as_sym(&s[0], "budget")?, s.len()) {
                        ("wall", 2) => w = Some(as_qty_d(&s[1], TIME, "budget.wall")?.value),
                        ("mem", 2) => m = Some(as_bytes(&s[1], "budget.mem")?),
                        (o, _) => return Err(err(format!("unknown budget field `{o}`"))),
                    }
                }
                budget = Some((
                    w.filter(|&v| v > 0.0)
                        .ok_or_else(|| err("positive budget wall required"))?,
                    m.ok_or_else(|| err("budget mem required"))?,
                ));
            }
            "capability" => {
                let kw = kwargs(&items[1..], &["cores", "mem", "wall", "ops"], "capability")?;
                let cores = match &req(&kw, "cores", "capability")?.kind {
                    NodeKind::Count {
                        value,
                        unit: CountUnit::Cores,
                    } => value
                        .integral_count()
                        .filter(|&c| c > 0)
                        .ok_or_else(|| err("positive cores required"))?,
                    NodeKind::Int(v) if *v > 0 => u64::try_from(*v).unwrap_or(1),
                    _ => return Err(err("capability :cores must be positive")),
                };
                let ops = as_list(req(&kw, "ops", "capability")?, "ops")?
                    .iter()
                    .map(|n| as_sym(n, "op").map(str::to_string))
                    .collect::<Result<Vec<_>, _>>()?;
                capability = Some((
                    cores,
                    as_qty_d(req(&kw, "wall", "capability")?, TIME, "capability.wall")?.value,
                    as_bytes(req(&kw, "mem", "capability")?, "capability.mem")?,
                    ops,
                ));
            }
            "reference-state" if items.len() == 2 => {
                ref_state = Some(match as_sym(&items[1], "reference-state")? {
                    "unloaded" => ReferenceState::Unloaded,
                    "strung" => ReferenceState::Strung,
                    "assembled" => ReferenceState::Assembled,
                    "played" => ReferenceState::Played,
                    o => return Err(err(format!("unknown reference-state `{o}`"))),
                });
            }
            "environment" => {
                let mut e = DraftEnvironment::earth_lab_air();
                if !(items.len() == 2
                    && matches!(&items[1].kind, NodeKind::Symbol(s) if s == "earth-lab-air"))
                {
                    let kw = kwargs(
                        &items[1..],
                        &["ambient-temp", "ambient-pressure"],
                        "environment",
                    )?;
                    if let Some(t) = kw.get("ambient-temp") {
                        e.base.ambient_temperature = as_qty_d(t, TEMPERATURE, "ambient-temp")?;
                    }
                    if let Some(p) = kw.get("ambient-pressure") {
                        e.base.ambient_pressure = as_qty_d(p, PRESSURE, "ambient-pressure")?;
                    }
                }
                env = Some(e);
            }
            "time-plan" => {
                let kw = kwargs(&items[1..], &["horizon", "dt", "sample-rate"], "time-plan")?;
                time_plan = Some(TimePlan {
                    horizon: as_qty_d(req(&kw, "horizon", "time-plan")?, TIME, "horizon")?,
                    step_policy: StepPolicy::Fixed {
                        dt: as_qty_d(req(&kw, "dt", "time-plan")?, TIME, "dt")?,
                    },
                    default_sample_rate: as_qty_d(
                        req(&kw, "sample-rate", "time-plan")?,
                        FREQUENCY,
                        "sample-rate",
                    )?,
                });
            }
            "part" => {
                let pname = as_str(
                    items.get(1).ok_or_else(|| err("part missing name"))?,
                    "part.name",
                )?
                .to_string();
                let scope = format!("part `{pname}`");
                let kw = kwargs(
                    &items[2..],
                    &[
                        "kind",
                        "geometry",
                        "material",
                        "query-temp",
                        "physics",
                        "fidelity",
                        "initial",
                    ],
                    &scope,
                )?;
                let entity = catalog
                    .declare(EntityDeclaration::part(root, &pname))
                    .map_err(|e| err(format!("{scope}: {e}")))?;
                by_name.insert(pname.clone(), entity);
                entity_names.insert(entity, pname.clone());
                let kind = match as_sym(req(&kw, "kind", &scope)?, "part.kind")? {
                    "rigid-body" => DraftPartKind::RigidBody,
                    "solid" => DraftPartKind::Solid,
                    "spring-beam" => DraftPartKind::SpringBeam,
                    "driven-body" => DraftPartKind::DrivenBody,
                    o => return Err(err(format!("{scope} unsupported kind `{o}`"))),
                };
                let mut materials = Vec::new();
                if let Some(mn) = kw.get("material") {
                    let card = resolve_card(
                        packs,
                        CardPackKind::Material,
                        as_str(mn, "material")?,
                        &scope,
                    )?;
                    let qt = kw
                        .get("query-temp")
                        .map_or(Ok(QtyAny::new(293.15, TEMPERATURE)), |n| {
                            as_qty_d(n, TEMPERATURE, "query-temp")
                        })?;
                    materials.push(MaterialRef {
                        card,
                        query_point: vec![(ParameterKey::new("T", TEMPERATURE), qt)],
                        frame: FrameId(0),
                        layer: None,
                        unknown: Vec::new(),
                    });
                }
                let mut physics = BTreeSet::new();
                for d in as_list(req(&kw, "physics", &scope)?, "physics")? {
                    physics.insert(match as_sym(d, "domain")? {
                        "mechanics" => PhysicsDomain::Mechanics,
                        "thermal" => PhysicsDomain::Thermal,
                        "electromagnetics" => PhysicsDomain::Electromagnetics,
                        "fluid" => PhysicsDomain::Fluid,
                        "acoustics" => PhysicsDomain::Acoustics,
                        "tribology" => PhysicsDomain::Tribology,
                        o => return Err(err(format!("unknown domain `{o}`"))),
                    });
                }
                if as_sym(req(&kw, "fidelity", &scope)?, "fidelity")? != "lumped" {
                    return Err(err(format!("{scope} only :fidelity lumped is supported")));
                }
                let initial = kw
                    .get("initial")
                    .map(|n| parse_initial(n, &scope))
                    .transpose()?;
                parts.push(DraftPartDecl {
                    name: pname,
                    entity,
                    kind,
                    geometry: parse_geom(req(&kw, "geometry", &scope)?)?,
                    materials,
                    physics,
                    fidelity: FidelityPolicy {
                        rungs: vec![FidelityRung::Lumped],
                        error_budget: FidelityErrorBudget {
                            relative_tol: 1.0e-5,
                            absolute_tol: QtyAny::dimensionless(1.0e-12),
                            max_work_units: Some(1_000_000),
                        },
                    },
                    uncertainty: Vec::new(),
                    initial,
                });
            }
            "joint" => {
                let jname = as_str(
                    items.get(1).ok_or_else(|| err("joint missing name"))?,
                    "joint.name",
                )?
                .to_string();
                let scope = format!("joint `{jname}`");
                let kw = kwargs(
                    &items[2..],
                    &["kind", "parent", "child", "axis", "damping", "stiffness"],
                    &scope,
                )?;
                let kind = match as_sym(req(&kw, "kind", &scope)?, "joint.kind")? {
                    "fixed" => JointKind::Fixed,
                    "revolute" => JointKind::Revolute,
                    "prismatic" => JointKind::Prismatic,
                    o => return Err(err(format!("{scope} unsupported kind `{o}`"))),
                };
                let parent = match &req(&kw, "parent", &scope)?.kind {
                    NodeKind::Symbol(s) if s == "root" => root,
                    NodeKind::Str(s) => *by_name
                        .get(s)
                        .ok_or_else(|| err(format!("{scope} unknown parent `{s}`")))?,
                    _ => return Err(err(format!("{scope} invalid parent"))),
                };
                let cn = as_str(req(&kw, "child", &scope)?, "child")?;
                let child = *by_name
                    .get(cn)
                    .ok_or_else(|| err(format!("{scope} unknown child `{cn}`")))?;
                joints.push(JointDecl {
                    name: jname,
                    parent,
                    child,
                    frame: FrameId(0),
                    kind,
                    axis: as_vec3(req(&kw, "axis", &scope)?, "axis")?,
                    pitch: None,
                    limits: None,
                    damping: kw
                        .get("damping")
                        .map(|n| as_qty(n, "damping"))
                        .transpose()?,
                    stiffness: kw
                        .get("stiffness")
                        .map(|n| as_qty(n, "stiffness"))
                        .transpose()?,
                });
            }
            "interface" => {
                let iname = as_str(
                    items.get(1).ok_or_else(|| err("interface missing name"))?,
                    "interface.name",
                )?
                .to_string();
                let scope = format!("interface `{iname}`");
                let kw = kwargs(
                    &items[2..],
                    &[
                        "kind",
                        "a",
                        "b",
                        "card",
                        "family",
                        "area",
                        "normal-pressure",
                        "thickness",
                        "fraction",
                    ],
                    &scope,
                )?;
                let ea = *by_name
                    .get(as_str(req(&kw, "a", &scope)?, "a")?)
                    .ok_or_else(|| err(format!("{scope} unknown part a")))?;
                let eb = *by_name
                    .get(as_str(req(&kw, "b", &scope)?, "b")?)
                    .ok_or_else(|| err(format!("{scope} unknown part b")))?;
                let sa = catalog
                    .declare(EntityDeclaration::surface(ea, &format!("{iname}_a")))
                    .map_err(|e| err(format!("{scope}: {e}")))?;
                let sb = catalog
                    .declare(EntityDeclaration::surface(eb, &format!("{iname}_b")))
                    .map_err(|e| err(format!("{scope}: {e}")))?;
                let entity = catalog
                    .declare(EntityDeclaration::interface(
                        root,
                        &iname,
                        InterfacePair::unordered(sa, sb),
                    ))
                    .map_err(|e| err(format!("{scope}: {e}")))?;
                entity_names.insert(entity, iname.clone());
                let mut card = || -> Result<ContentHash, SolveRefusal> {
                    let card = resolve_card(
                        packs,
                        CardPackKind::Interface,
                        as_str(req(&kw, "card", &scope)?, "card")?,
                        &scope,
                    )?;
                    let family = as_str(req(&kw, "family", &scope)?, "family")?;
                    families.insert(card, family.to_string());
                    Ok(card)
                };
                let mut state = Vec::new();
                for (flag, key, dims) in [
                    ("area", "area", AREA),
                    ("normal-pressure", "normal_pressure", PRESSURE),
                    ("thickness", "thickness", LENGTH),
                ] {
                    if let Some(n) = kw.get(flag) {
                        state.push((ParameterKey::new(key, dims), as_qty_d(n, dims, flag)?));
                    }
                }
                if let Some(n) = kw.get("fraction") {
                    let fraction = QtyAny::dimensionless(as_f64(n, "fraction")?);
                    state.push((ParameterKey::new("fraction", Dims::NONE), fraction));
                }
                let kind = match as_sym(req(&kw, "kind", &scope)?, "interface.kind")? {
                    // A heat route is bookkeeping, not a material interface: no card.
                    "heat-route" if kw.contains_key("card") || kw.contains_key("family") => {
                        return Err(err(format!("{scope}: a heat route takes no card")));
                    }
                    "heat-route" => DraftInterfaceKind::HeatRoute { state },
                    "thermal-contact" => DraftInterfaceKind::ThermalContact {
                        law_card: card()?,
                        third_body: None,
                        state,
                    },
                    "mechanical-contact" => DraftInterfaceKind::MechanicalContact {
                        law_card: card()?,
                        third_body: None,
                        state,
                    },
                    "glue" => DraftInterfaceKind::Glue { law_card: card()? },
                    o => return Err(err(format!("{scope} unsupported kind `{o}`"))),
                };
                interfaces.push(DraftInterfaceDecl {
                    name: iname,
                    entity,
                    feature_a: sa,
                    feature_b: sb,
                    kind,
                });
            }
            "source" => {
                let sname = as_str(
                    items.get(1).ok_or_else(|| err("source missing name"))?,
                    "source.name",
                )?
                .to_string();
                let scope = format!("source `{sname}`");
                let kw = kwargs(&items[2..], &["part", "port", "history"], &scope)?;
                let entity = *by_name
                    .get(as_str(req(&kw, "part", &scope)?, "part")?)
                    .ok_or_else(|| err(format!("{scope} unknown part")))?;
                sources.push(SourceDecl {
                    name: sname,
                    target: SourceTarget::Port {
                        entity,
                        port: as_str(req(&kw, "port", &scope)?, "port")?.to_string(),
                    },
                    history: parse_history(req(&kw, "history", &scope)?)?,
                });
            }
            "probe" => probes.push(parse_probe(items, &by_name)?),
            "requirement" => {
                let id = as_str(
                    items.get(1).ok_or_else(|| err("requirement missing id"))?,
                    "requirement.id",
                )?
                .to_string();
                if id.trim().is_empty() || !req_ids.insert(id.clone()) {
                    return Err(err(format!("duplicate or empty requirement `{id}`")));
                }
                let sc = format!("requirement `{id}`");
                let kw = kwargs(&items[2..], &["probe", "stat", "relation", "limit"], &sc)?;
                let probe = as_str(req(&kw, "probe", &sc)?, "probe")?.to_string();
                let stat = match as_sym(req(&kw, "stat", &sc)?, "stat")? {
                    "max" => ProbeStat::Max,
                    "min" => ProbeStat::Min,
                    "final" => ProbeStat::Final,
                    o => return Err(err(format!("{sc} unsupported stat `{o}`"))),
                };
                let relation = match as_sym(req(&kw, "relation", &sc)?, "relation")? {
                    "at-most" => RequirementRelation::AtMost,
                    "at-least" => RequirementRelation::AtLeast,
                    o => return Err(err(format!("{sc} unsupported relation `{o}`"))),
                };
                let limit = as_qty(req(&kw, "limit", &sc)?, "limit")?;
                requirements.push(AssemblyRequirement {
                    id,
                    probe,
                    stat,
                    relation,
                    limit,
                });
            }
            "analysis-options" => {
                let kw = kwargs(
                    &items[1..],
                    &["harmonic-frequencies", "modal-max"],
                    "analysis-options",
                )?;
                if let Some(fs) = kw.get("harmonic-frequencies") {
                    for f in as_list(fs, "harmonic-frequencies")? {
                        let hz = as_f64(f, "freq")?;
                        if hz <= 0.0 {
                            return Err(err("harmonic frequencies must be positive"));
                        }
                        harm_hz.push(hz);
                    }
                }
                if let Some(NodeKind::Int(v)) = kw.get("modal-max").map(|n| &n.kind)
                    && *v > 0
                {
                    modal_max = usize::try_from(*v).unwrap_or(8);
                }
            }
            o => return Err(err(format!("unknown assembly clause `{o}`"))),
        }
    }

    if version != Some(1) || !units_ok {
        return Err(err("assembly requires (version 1) and (units \"SI\")"));
    }
    let seed = seed.ok_or_else(|| err("missing required (seed 0x...)"))?;
    let (constellation, workspace) =
        versions.ok_or_else(|| err("missing required (versions ...)"))?;
    let (wall_budget_s, mem_budget_bytes) =
        budget.ok_or_else(|| err("missing required (budget ...)"))?;
    let (cores, cap_w, cap_m, granted_ops) =
        capability.ok_or_else(|| err("missing required (capability ...)"))?;
    if cap_w < wall_budget_s || cap_m < mem_budget_bytes {
        return Err(err("capability grant does not cover budget"));
    }
    let environment = env.ok_or_else(|| err("missing required (environment ...)"))?;
    let time_plan = time_plan.ok_or_else(|| err("missing required (time-plan ...)"))?;
    let scenario = Scenario::new(&name, seed, environment.base.clone());
    let draft = DraftAssembly::admit(
        DraftAssemblyInput {
            name: name.clone(),
            root,
            catalog,
            reference_state: ref_state,
            environment,
            time_plan,
            parts,
            joints,
            gear_meshes: Vec::new(),
            interfaces,
            sources,
            controllers: Vec::new(),
            probes,
        },
        &scenario,
    )
    .map_err(|e| err(format!("structural admission refused: {e}")))?;
    for r in &requirements {
        if !draft.input().probes.iter().any(|p| p.name == r.probe) {
            return Err(err(format!(
                "requirement `{}` references unknown probe `{}`",
                r.id, r.probe
            )));
        }
    }
    let card_store = build_card_store(packs, &families)?;
    Ok(AdmittedAssembly {
        canonical_source,
        source_hash,
        name,
        seed,
        constellation,
        workspace,
        wall_budget_s,
        mem_budget_bytes,
        cores,
        granted_ops,
        draft,
        card_store,
        entity_names,
        harmonic_frequencies_hz: harm_hz,
        modal_max_modes: modal_max,
        requirements,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);

    fn admit(attrs: &str) -> Result<AdmittedAssembly, SolveRefusal> {
        let src = format!(
            "(assembly \"init_check\" (version 1) (units \"SI\") (seed 0x1) \
             (versions (constellation \"2026-07\") (workspace \"pkg-1c\")) \
             (budget (wall 60s) (mem 64MiB)) \
             (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*)) \
             (reference-state unloaded) (environment earth-lab-air) \
             (time-plan :horizon 0.05s :dt 0.001s :sample-rate 1000Hz) \
             (part \"rotor\" :kind rigid-body :geometry (box 0.1m 0.05m 0.05m) \
             :physics (mechanics thermal) :fidelity lumped {attrs}))\n"
        );
        let packs = CardPackSet::admit(Vec::new()).unwrap();
        admit_assembly_source(src.as_bytes(), &packs)
    }

    #[test]
    fn part_initial_round_trips_sexpr_and_json_with_rad_per_s_and_optional_temp() {
        let q = QtyAny::new;
        let full = ":initial (:position (0.1m -0.2m 0.05m) :orientation (0.6 0.0 0.8 0.0) \
                    :linear-velocity (1.5m/s 0m/s -0.5m/s) \
                    :angular-velocity (0rad/s -2.5rad/s 120rad/s) :temperature 320K)";
        let expected = DraftInitialState {
            position: [q(0.1, LENGTH), q(-0.2, LENGTH), q(0.05, LENGTH)],
            orientation: [0.6, 0.0, 0.8, 0.0],
            linear_velocity: [q(1.5, VELOCITY), q(0.0, VELOCITY), q(-0.5, VELOCITY)],
            angular_velocity: [q(0.0, FREQUENCY), q(-2.5, FREQUENCY), q(120.0, FREQUENCY)],
            temperature: Some(q(320.0, TEMPERATURE)),
        };
        let packs = CardPackSet::admit(Vec::new()).unwrap();
        let a1 = admit(full).expect("full :initial admits");
        let init1 = a1.draft.input().parts[0].initial.as_ref().unwrap();
        assert_eq!(init1, &expected);
        assert_eq!(init1.angular_velocity[2].dims, FREQUENCY);

        let a2 = admit_assembly_source(&a1.canonical_source, &packs).unwrap();
        assert_eq!(a2.draft.input().parts[0].initial.as_ref(), Some(&expected));

        let ast = fs_ir::sexpr::parse(std::str::from_utf8(&a1.canonical_source).unwrap()).unwrap();
        let json = fs_ir::json::print_checked(&ast).unwrap();
        let a3 = admit_assembly_source(json.as_bytes(), &packs).unwrap();
        assert_eq!(a3.draft.input().parts[0].initial.as_ref(), Some(&expected));

        let no_temp = full.replace(" :temperature 320K", "");
        let a4 = admit(&no_temp).expect("omitted :temperature admits");
        let want_no_temp = DraftInitialState {
            temperature: None,
            ..expected
        };
        assert_eq!(
            a4.draft.input().parts[0].initial.as_ref(),
            Some(&want_no_temp)
        );
        assert_eq!(admit("").unwrap().draft.input().parts[0].initial, None);
    }

    #[test]
    fn part_initial_refusals_report_shape_and_admission_item_names() {
        let base = ":position (0m 0m 0m) :orientation (1 0 0 0) \
                    :linear-velocity (0m/s 0m/s 0m/s) :angular-velocity (0rad/s 0rad/s 0rad/s)";
        let rep = |old: &str, new: &str| format!(":initial ({})", base.replace(old, new));
        let cases = [
            (
                rep("(1 0 0 0)", "(2.0 0.0 0.0 0.0)"),
                "part:rotor.initial.orientation",
            ),
            (
                format!(":initial ({base} :temperature 0K)"),
                "part:rotor.initial.temperature",
            ),
            (
                rep("(0m 0m 0m)", "(1K 0m 0m)"),
                "part:rotor.initial.position[0]",
            ),
            (
                rep("(0m/s 0m/s 0m/s)", "(1m 0m/s 0m/s)"),
                "part:rotor.initial.linear_velocity[0]",
            ),
            (
                rep("(0rad/s 0rad/s 0rad/s)", "(0rad/s 1K 0rad/s)"),
                "part:rotor.initial.angular_velocity[1]",
            ),
            (
                format!(":initial ({base} :temperature 1m)"),
                "part:rotor.initial.temperature",
            ),
            (
                rep("(0m 0m 0m)", "(0m 0m)"),
                "initial.position must have 3 elements",
            ),
            (
                rep("(1 0 0 0)", "(1 0 0)"),
                "initial.orientation must have 4 elements",
            ),
            (
                rep("(1 0 0 0)", "(1 \"x\" 0 0)"),
                "expected finite number for initial.orientation",
            ),
            (
                rep(":linear-velocity (0m/s 0m/s 0m/s) ", ""),
                "missing `:linear-velocity`",
            ),
            (
                format!(":initial ({base} :bogus 1m)"),
                "rejects unknown keyword `:bogus`",
            ),
            (
                format!(":initial ({base}) :initial ({base})"),
                "repeats `:initial`",
            ),
        ];
        for (attrs, needle) in cases {
            let Err(r) = admit(&attrs) else {
                panic!("expected refusal for `{attrs}`");
            };
            assert_eq!((r.code, r.stage), ("cli-assembly-admit", Some("admit")));
            assert!(
                r.what.contains(needle),
                "got `{}` without `{needle}`",
                r.what
            );
        }
    }

    #[test]
    fn requirement_clause_parses_and_refuses_bad_probe_stat_relation_and_dup() {
        let with = |extra: &str| {
            let src = format!(
                "(assembly \"req_check\" (version 1) (units \"SI\") (seed 0x1) \
                 (versions (constellation \"2026-07\") (workspace \"t\")) \
                 (budget (wall 60s) (mem 64MiB)) \
                 (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*)) \
                 (reference-state unloaded) (environment earth-lab-air) \
                 (time-plan :horizon 0.05s :dt 0.001s :sample-rate 1000Hz) \
                 (part \"p\" :kind rigid-body :geometry (box 0.1m 0.05m 0.05m) :physics (thermal) :fidelity lumped) \
                 (probe \"tp\" :kind thermocouple :target \"p\" :sample-rate 1000Hz) {extra})\n"
            );
            admit_assembly_source(src.as_bytes(), &CardPackSet::admit(Vec::new()).unwrap())
        };
        let ok = with("(requirement \"r1\" :probe \"tp\" :stat max :relation at-most :limit 310K)")
            .unwrap();
        assert_eq!(
            (
                ok.requirements.len(),
                ok.requirements[0].id.as_str(),
                ok.requirements[0].stat
            ),
            (1, "r1", ProbeStat::Max)
        );
        for (bad, needle) in [
            (
                "(requirement \"r1\" :probe \"missing\" :stat max :relation at-most :limit 310K)",
                "unknown probe",
            ),
            (
                "(requirement \"r1\" :probe \"tp\" :stat mean :relation at-most :limit 310K)",
                "unsupported stat",
            ),
            (
                "(requirement \"r1\" :probe \"tp\" :stat min :relation eq :limit 310K)",
                "unsupported relation",
            ),
            (
                "(requirement \"r1\" :probe \"tp\" :stat final :relation at-least :limit 290K) (requirement \"r1\" :probe \"tp\" :stat max :relation at-most :limit 310K)",
                "duplicate",
            ),
        ] {
            let Err(e) = with(bad) else {
                panic!("expected refusal for {bad}")
            };
            assert!(e.what.contains(needle), "`{}` missing `{needle}`", e.what);
        }
    }

    #[test]
    fn heat_route_takes_a_fraction_and_no_card() {
        // `admit` wraps its text in the rotor's part form; this closes it.
        let route = |extra: &str| {
            admit(&format!(
                ") (part \"cell\" :kind solid :geometry (box 0.1m 0.1m 0.1m) \
                 :physics (thermal) :fidelity lumped) (interface \"route\" :kind heat-route \
                 :a \"rotor\" :b \"cell\" :fraction 0.7{extra}"
            ))
        };
        let a = route("").expect("a heat route admits without a card");
        let key = ParameterKey::new("fraction", Dims::NONE);
        let state = vec![(key, QtyAny::dimensionless(0.7))];
        let kind = &a.draft.input().interfaces[0].kind;
        assert_eq!(kind, &DraftInterfaceKind::HeatRoute { state });
        let r = route(" :card \"x\" :family \"y\"").err().unwrap();
        assert!(r.what.contains("a heat route takes no card"), "{}", r.what);
    }
}
