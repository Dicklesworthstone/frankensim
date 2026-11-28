//! Assembly `.fsim` (s-expression & isomorphic JSON) parser, Five Explicits
//! admission, and `PackCatalog -> fs_system::CardStore` bridge.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{ContentHash, hash_bytes};
use fs_evidence::uncertainty::RequirementRelation;
use fs_matdb::{NormalizedInterfacePack, NormalizedMaterialCardPack, PackError};
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::draft::{
    ControllerDecl, ControllerLaw, DraftAssembly, DraftAssemblyInput, DraftEnvironment,
    DraftGeometrySource, DraftInitialState, DraftInterfaceDecl, DraftInterfaceKind, DraftPartDecl,
    DraftPartKind, FidelityErrorBudget, FidelityPolicy, FidelityRung, JointDecl, JointKind,
    MicrophonePattern, PhysicsDomain, PrimitiveShape, ProbeDecl, ProbeKind, SourceDecl,
    SourceTarget, StepPolicy, TensionedDiskGenerator, ThirdBodyDecl, TimeHistory, TimePlan,
    TimeScheme, TimingEdge,
};
use fs_scenario::assembly::{DimensionBound, MaterialRef, ParameterKey, ReferenceState};
use fs_scenario::entity::ToleranceSource;
use fs_scenario::{EntityCatalog, EntityDeclaration, EntityId, FrameId, InterfacePair, Scenario};
use fs_system::CardStore;

use crate::ast::{CountUnit, Node, NodeKind};

/// Selectable analysis mode for `--analysis transient|stationary|harmonic|modal`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyAnalysisMode {
    /// Coupled transient stepping via `System::step`.
    Transient,
    /// Stationary equilibrium via `System::run_analysis(Stationary)`.
    Stationary,
    /// Linearised frequency response via `System::run_analysis(Harmonic)`.
    Harmonic,
    /// Linearised eigenmodes via `System::run_analysis(Modal)`.
    Modal,
}

impl AssemblyAnalysisMode {
    /// Canonical lowercase mode identifier.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Transient => "transient",
            Self::Stationary => "stationary",
            Self::Harmonic => "harmonic",
            Self::Modal => "modal",
        }
    }

    /// Parse a canonical lowercase analysis mode label.
    #[must_use]
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "transient" => Some(Self::Transient),
            "stationary" => Some(Self::Stationary),
            "harmonic" => Some(Self::Harmonic),
            "modal" => Some(Self::Modal),
            _ => None,
        }
    }
}

/// Which normalized `fs-matdb` pack envelope a card belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackKind {
    /// `FSMCDPK\0` material card pack.
    Material,
    /// `FSINTPK\0` interface card pack.
    Interface,
}

impl PackKind {
    /// Lowercase diagnostic label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Material => "material",
            Self::Interface => "interface",
        }
    }

    /// CLI flag supplying this pack kind.
    #[must_use]
    pub const fn flag(self) -> &'static str {
        match self {
            Self::Material => "--materials",
            Self::Interface => "--interfaces",
        }
    }
}

/// One admitted normalized card pack entry.
#[derive(Debug, Clone)]
pub struct PackEntry {
    /// Content hash of the reconstructed card.
    pub card: ContentHash,
    /// Human-readable card identity selector.
    pub identity: String,
    /// Exact canonical pack bytes.
    pub bytes: Vec<u8>,
}

/// Maximum distinct card packs admitted by [`PackCatalog::from_raw_bytes`]
/// (matches `fs-cli::cards::MAX_CARD_PACKS`).
pub const MAX_PACKS: usize = 128;

/// Maximum byte length of a single card pack accepted by [`PackCatalog::from_raw_bytes`]
/// (matches `fs-cli::cards::MAX_CARD_PACK_BYTES`, 4 MiB).
pub const MAX_PACK_BYTES: usize = 4 * 1024 * 1024;

/// Admitted material and interface card packs passed to [`admit_assembly_source`].
#[derive(Debug, Clone, Default)]
pub struct PackCatalog {
    materials: Vec<PackEntry>,
    interfaces: Vec<PackEntry>,
}

fn pack_refusal(code: &'static str, what: impl Into<String>, fix: &'static str) -> CodecRefusal {
    CodecRefusal {
        code,
        what: what.into(),
        fix: fix.to_owned(),
    }
}

fn map_pack_err(kind: &str, e: &PackError) -> CodecRefusal {
    let (code, fix) = match e {
        PackError::IdentityMismatch { .. } => (
            "cli-solve-card-pack-identity",
            "re-export the card pack from its canonical source",
        ),
        PackError::ResourceLimit { .. } => (
            "cli-solve-card-pack-size",
            "supply a normalized card pack within the 4 MiB ceiling",
        ),
        _ => (
            "cli-solve-card-pack-decode",
            "supply a valid canonical FSMCDPK\0 or FSINTPK\0 pack",
        ),
    };
    pack_refusal(code, format!("{kind} pack decode failed: {e}"), fix)
}

fn push_unique(
    out: &mut Vec<(ContentHash, PackEntry)>,
    total_distinct: usize,
    kind: &str,
    root: ContentHash,
    entry: PackEntry,
) -> Result<(), CodecRefusal> {
    if out.iter().any(|(r, _)| *r == root) {
        return Ok(());
    }
    if out.iter().any(|(_, e)| e.card == entry.card) {
        return Err(pack_refusal(
            "cli-solve-card-pack-conflict",
            format!("conflicting {kind} packs for card {}", entry.card.to_hex()),
            "supply only one pack per material or interface card",
        ));
    }
    if total_distinct >= MAX_PACKS {
        return Err(pack_refusal(
            "cli-solve-card-pack-count",
            format!("distinct card pack count exceeds limit {MAX_PACKS}"),
            "supply at most 128 distinct card packs",
        ));
    }
    out.push((root, entry));
    Ok(())
}

impl PackCatalog {
    /// Build a catalog from already-admitted material and interface entries.
    #[must_use]
    pub fn from_parts(materials: Vec<PackEntry>, interfaces: Vec<PackEntry>) -> Self {
        Self {
            materials,
            interfaces,
        }
    }

    /// Decode canonical `FSMCDPK\0` and `FSINTPK\0` byte buffers into a catalog.
    ///
    /// Matches `fs-cli::cards::CardPackSetBuilder`: enforces `MAX_PACK_BYTES` (4 MiB)
    /// per pack before decoding, collapses byte-identical repeats before enforcing
    /// `MAX_PACKS` (128) distinct packs, and returns `cli-solve-card-pack-*` codes.
    pub fn from_raw_bytes(packs: &[&[u8]]) -> Result<Self, CodecRefusal> {
        let (mut mats, mut ints, mut seen_bytes) = (Vec::new(), Vec::new(), BTreeSet::new());
        for &raw in packs {
            if raw.len() > MAX_PACK_BYTES {
                return Err(pack_refusal(
                    "cli-solve-card-pack-size",
                    format!(
                        "card pack is {} bytes, above the {MAX_PACK_BYTES}-byte ceiling",
                        raw.len()
                    ),
                    "supply a normalized card pack within the 4 MiB ceiling",
                ));
            }
            if !seen_bytes.insert(hash_bytes(raw)) {
                continue;
            }
            let total = mats.len() + ints.len();
            match raw.get(..8) {
                Some(b"FSMCDPK\0") => {
                    let dec = NormalizedMaterialCardPack::from_bytes(raw)
                        .map_err(|e| map_pack_err("material", &e))?;
                    let e = PackEntry {
                        card: dec.card().content_hash(),
                        identity: dec.card().id().to_string(),
                        bytes: raw.to_vec(),
                    };
                    push_unique(&mut mats, total, "material", dec.content_hash(), e)?;
                }
                Some(b"FSINTPK\0") => {
                    let dec = NormalizedInterfacePack::from_bytes(raw)
                        .map_err(|e| map_pack_err("interface", &e))?;
                    let id = format!(
                        "{} | {} | {}",
                        dec.card().surface_a().material,
                        dec.card().surface_b().material,
                        dec.card().medium()
                    );
                    let e = PackEntry {
                        card: dec.card().content_hash(),
                        identity: id,
                        bytes: raw.to_vec(),
                    };
                    push_unique(&mut ints, total, "interface", dec.content_hash(), e)?;
                }
                _ => {
                    return Err(pack_refusal(
                        "cli-solve-card-pack-decode",
                        "unknown card pack magic (expected FSMCDPK\0 or FSINTPK\0)",
                        "supply a valid canonical FSMCDPK\0 or FSINTPK\0 pack",
                    ));
                }
            }
        }
        mats.sort_by(|a, b| a.0.as_bytes().cmp(b.0.as_bytes()));
        ints.sort_by(|a, b| a.0.as_bytes().cmp(b.0.as_bytes()));
        Ok(Self {
            materials: mats.into_iter().map(|(_, e)| e).collect(),
            interfaces: ints.into_iter().map(|(_, e)| e).collect(),
        })
    }

    /// Admitted material packs.
    #[must_use]
    pub fn materials(&self) -> &[PackEntry] {
        &self.materials
    }

    /// Admitted interface packs.
    #[must_use]
    pub fn interfaces(&self) -> &[PackEntry] {
        &self.interfaces
    }
}

/// Structured refusal surfaced by `.fsim` assembly admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CodecRefusal {
    /// Stable machine refusal code.
    pub code: &'static str,
    /// Human-readable description of what failed.
    pub what: String,
    /// Actionable fix hint.
    pub fix: String,
}

impl fmt::Display for CodecRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}; fix: {}", self.code, self.what, self.fix)
    }
}

impl std::error::Error for CodecRefusal {}

type SolveRefusal = CodecRefusal;
type CardPackSet = PackCatalog;
type CardPackKind = PackKind;

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const FREQUENCY: Dims = Dims([0, 0, -1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const AREA: Dims = Dims([2, 0, 0, 0, 0, 0]);
const FORCE: Dims = Dims([1, 1, -2, 0, 0, 0]);
const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);
const SURFACE_TENSION: Dims = Dims([0, 1, -2, 0, 0, 0]);

/// Scalar summary statistic applied to a probe series by a requirement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProbeStat {
    /// Maximum finite value.
    Max,
    /// Minimum finite value.
    Min,
    /// Final sample.
    Final,
}

impl ProbeStat {
    /// Canonical lowercase label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Max => "max",
            Self::Min => "min",
            Self::Final => "final",
        }
    }
}

/// Declared requirement linking a probe summary to a threshold.
#[derive(Debug, Clone, PartialEq)]
pub struct AssemblyRequirement {
    /// Unique requirement id.
    pub id: String,
    /// Optional analysis mode override.
    pub analysis: Option<AssemblyAnalysisMode>,
    /// Target probe name.
    pub probe: String,
    /// Reduction statistic.
    pub stat: ProbeStat,
    /// Bound direction.
    pub relation: RequirementRelation,
    /// Dimensioned limit.
    pub limit: QtyAny,
}

/// Canonical box dimension axis labels.
pub const BOX_AXES: [&str; 3] = ["lx", "ly", "lz"];

/// Declared geometric tolerance envelope on one box axis.
#[derive(Debug, Clone, PartialEq)]
pub struct GeomEnvelope {
    /// Part name.
    pub part: String,
    /// Box axis index (`0..3`).
    pub axis: usize,
    /// Positive half-width [m].
    pub half_width_m: f64,
    /// Cited tolerance source.
    pub source: String,
}

/// Declared uncertainty envelope on a boundary condition or probe.
#[derive(Debug, Clone, PartialEq)]
pub struct TargetEnvelope {
    /// Target source or probe name.
    pub target: String,
    /// Optional half-width (`None` when covered-by-source).
    pub half_width: Option<QtyAny>,
    /// Cited source string.
    pub source: String,
}

/// Uncertainty envelopes declared in an assembly document.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DeclaredEnvelopes {
    /// Geometry tolerance envelopes.
    pub geometry: Vec<GeomEnvelope>,
    /// Boundary-condition envelopes.
    pub boundary_conditions: Vec<TargetEnvelope>,
    /// Measurement envelopes.
    pub measurement: Vec<TargetEnvelope>,
    /// Model-form envelopes.
    pub model_form: Vec<TargetEnvelope>,
}

/// Admitted assembly document, Five Explicits, and resolved card store.
pub struct AdmittedAssembly {
    /// Canonical s-expression bytes.
    pub canonical_source: Vec<u8>,
    /// BLAKE3 hash of `canonical_source`.
    pub source_hash: ContentHash,
    /// Assembly name.
    pub name: String,
    /// Deterministic seed.
    pub seed: u64,
    /// Constellation version.
    pub constellation: String,
    /// Workspace version.
    pub workspace: String,
    /// Wall-clock budget [s].
    pub wall_budget_s: f64,
    /// Memory budget [bytes].
    pub mem_budget_bytes: u64,
    /// Granted core count.
    pub cores: u64,
    /// Granted operation patterns.
    pub granted_ops: Vec<String>,
    /// Admitted structural draft assembly.
    pub draft: DraftAssembly,
    /// Resolved material and interface card store.
    pub card_store: CardStore,
    /// Entity id to declared name map.
    pub entity_names: BTreeMap<EntityId, String>,
    /// Harmonic sweep frequencies [Hz].
    pub harmonic_frequencies_hz: Vec<f64>,
    /// Maximum modal count requested.
    pub modal_max_modes: usize,
    /// Declared requirements.
    pub requirements: Vec<AssemblyRequirement>,
    /// Declared uncertainty envelopes.
    pub envelopes: DeclaredEnvelopes,
}

fn err(what: impl Into<String>) -> CodecRefusal {
    CodecRefusal {
        code: "cli-assembly-admit",
        what: what.into(),
        fix: "fix the assembly declaration".into(),
    }
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

fn parse_query_pairs(
    n: &Node,
    scope: &str,
    flag: &str,
) -> Result<Vec<(ParameterKey, QtyAny)>, SolveRefusal> {
    let at = as_list(n, flag)?;
    if !at.len().is_multiple_of(2) {
        return Err(err(format!("{scope}: :{flag} pairs keys with values")));
    }
    at.as_chunks::<2>()
        .0
        .iter()
        .map(|pair| {
            let key = as_sym(&pair[0], &format!("{flag} key"))?;
            let q = as_qty(&pair[1], key)
                .or_else(|_| as_f64(&pair[1], key).map(QtyAny::dimensionless))?;
            Ok((ParameterKey::new(key, q.dims), q))
        })
        .collect()
}

/// A relaxation bushing's branch memory, `:memory (displacement-<k>-<axis> value ...)`: each
/// pair is the state key `displacement_<k>_<axis>`, in m on x, y and z and in rad on rx, ry and rz,
/// `<k>` a branch's number or `eq`, the equilibrium spring.
fn parse_memory(n: &Node, scope: &str) -> Result<Vec<(ParameterKey, QtyAny)>, SolveRefusal> {
    let pairs = as_list(n, "memory")?;
    if !pairs.len().is_multiple_of(2) {
        return Err(err(format!("{scope}: :memory pairs keys with values")));
    }
    let mut state: Vec<(ParameterKey, QtyAny)> = Vec::new();
    for pair in pairs.as_chunks::<2>().0 {
        let key = as_sym(&pair[0], "memory key")?;
        let rest = key.strip_prefix("displacement-");
        let split = rest.and_then(|rest| rest.split_once('-'));
        let branch = |k: &str| k == "eq" || k.parse::<usize>().is_ok_and(|n| n.to_string() == k);
        let dims = match split {
            Some((k, "x" | "y" | "z")) if branch(k) => LENGTH,
            Some((k, "rx" | "ry" | "rz")) if branch(k) => Dims::NONE,
            _ => return Err(err(format!("{scope}: no :memory key `{key}`"))),
        };
        let name = key.replace('-', "_");
        if state.iter().any(|(have, _)| have.name == name) {
            return Err(err(format!("{scope} repeats `{key}` in :memory")));
        }
        let value = as_qty_d(&pair[1], dims, key)?;
        state.push((ParameterKey::new(name, dims), value));
    }
    Ok(state)
}

fn parse_fidelity(n: &Node, scope: &str) -> Result<FidelityRung, SolveRefusal> {
    match &n.kind {
        NodeKind::Symbol(s) if s == "lumped" => Ok(FidelityRung::Lumped),
        NodeKind::List(v) if !v.is_empty() => match as_sym(&v[0], "fidelity")? {
            "one-d-field" if v.len() == 2 => match v[1].kind {
                NodeKind::Int(1) => Ok(FidelityRung::OneDField { cells: 1 }),
                NodeKind::Int(n) => Err(err(format!(
                    "{scope}: :fidelity (one-d-field {n}) unsupported"
                ))),
                _ => Err(err(format!(
                    "{scope}: :fidelity (one-d-field ...) needs integer cells"
                ))),
            },
            "modal" if v.len() == 2 => match v[1].kind {
                NodeKind::Int(n) if n > 0 => u32::try_from(n)
                    .map(|mode_count| FidelityRung::Modal { mode_count })
                    .map_err(|_| err(format!("{scope}: :fidelity (modal {n}) exceeds u32"))),
                NodeKind::Int(n) => Err(err(format!(
                    "{scope}: :fidelity (modal {n}) requires positive mode_count"
                ))),
                _ => Err(err(format!(
                    "{scope}: :fidelity (modal ...) requires integer mode_count"
                ))),
            },
            "modal" => Err(err(format!(
                "{scope}: :fidelity (modal <mode_count>) requires exactly one mode_count"
            ))),
            _ => Err(err(format!("{scope}: unsupported :fidelity"))),
        },
        NodeKind::Symbol(s) => Err(err(format!("{scope}: unsupported :fidelity `{s}`"))),
        _ => Err(err(format!("{scope}: unsupported :fidelity"))),
    }
}

fn parse_geom(n: &Node) -> Result<DraftGeometrySource, SolveRefusal> {
    let v = as_list(n, "geometry")?;
    let head = v.first().ok_or_else(|| err("empty geometry"))?;
    if as_sym(head, "geometry")? == "tensioned-disk" {
        let scope = "geometry tensioned-disk";
        let kw = kwargs(
            &v[1..],
            &[
                "radius",
                "thickness",
                "tension",
                "radial-intervals",
                "azimuths",
            ],
            scope,
        )?;
        let pos_u32 = |k: &str| -> Result<u32, SolveRefusal> {
            match req(&kw, k, scope)?.kind {
                NodeKind::Int(i) if i > 0 => {
                    u32::try_from(i).map_err(|_| err(format!("{scope}: :{k} `{i}` exceeds u32")))
                }
                _ => Err(err(format!("{scope}: :{k} requires a positive integer"))),
            }
        };
        return Ok(DraftGeometrySource::TensionedDisk(TensionedDiskGenerator {
            radius: as_qty_d(req(&kw, "radius", scope)?, LENGTH, "tensioned-disk.radius")?,
            thickness: as_qty_d(
                req(&kw, "thickness", scope)?,
                LENGTH,
                "tensioned-disk.thickness",
            )?,
            tension: as_qty_d(
                req(&kw, "tension", scope)?,
                SURFACE_TENSION,
                "tensioned-disk.tension",
            )?,
            radial_intervals: pos_u32("radial-intervals")?,
            azimuths: pos_u32("azimuths")?,
        }));
    }
    Ok(DraftGeometrySource::Primitive(
        match (as_sym(head, "geometry")?, v.len()) {
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
        "table" if v.len() >= 3 => {
            let mut knots = Vec::with_capacity(v.len() - 1);
            for row in &v[1..] {
                let pair = as_list(row, "table knot")?;
                if pair.len() != 2 {
                    return Err(err("table knot must be (time value)"));
                }
                knots.push((
                    as_qty_d(&pair[0], TIME, "table.time")?,
                    as_qty(&pair[1], "table.value")?,
                ));
            }
            Ok(TimeHistory::PiecewiseLinear { knots })
        }
        "table" => Err(err(
            "table history requires at least two (time value) knots",
        )),
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
            p.card.to_hex() == sel
                || p.identity == sel
                || p.identity.starts_with(&format!("{sel}/"))
        })
        .map(|p| p.card)
        .collect();
    match hits.as_slice() {
        [only] => Ok(*only),
        [] => Err(CodecRefusal {
            code: "cli-assembly-card-missing",
            what: format!("{scope} missing {} card `{sel}`", kind.label()),
            fix: format!("supply the card pack via {}", kind.flag()),
        }),
        _ => Err(CodecRefusal {
            code: "cli-assembly-card-missing",
            what: format!("{scope} card selector `{sel}` is ambiguous"),
            fix: "specify full identity or 64-hex hash".into(),
        }),
    }
}

fn build_card_store(
    packs: &CardPackSet,
    families: &BTreeMap<ContentHash, String>,
) -> Result<CardStore, SolveRefusal> {
    let mut store = CardStore::new();
    for p in packs.materials() {
        let dec = NormalizedMaterialCardPack::from_bytes(&p.bytes)
            .map_err(|e| err(format!("material pack decode failed: {e}")))?;
        store.insert_material(dec.card().clone());
    }
    for p in packs.interfaces() {
        let dec = NormalizedInterfacePack::from_bytes(&p.bytes)
            .map_err(|e| err(format!("interface pack decode failed: {e}")))?;
        let fam = families
            .get(&p.card)
            .cloned()
            .unwrap_or_else(|| "thermal-conductance".into());
        store
            .insert_interface(fam, dec.card().clone())
            .map_err(|e| err(format!("interface store error: {e}")))?;
    }
    Ok(store)
}

/// `(controller "name" :sensor ("part" "signal") :actuator ("part" "port")
/// :sample-rate <Hz> :law (pid :kp .. :ki .. :kd .. :setpoint ..))`. The gains'
/// units are checked against the sensor and actuator when the assembly lowers.
fn parse_controller(
    items: &[Node],
    parts: &BTreeMap<String, EntityId>,
) -> Result<ControllerDecl, SolveRefusal> {
    let name = as_str(
        items.get(1).ok_or_else(|| err("controller missing name"))?,
        "controller.name",
    )?
    .to_string();
    let scope = format!("controller `{name}`");
    let kw = kwargs(
        &items[2..],
        &["sensor", "actuator", "sample-rate", "law"],
        &scope,
    )?;
    let port = |key: &str| -> Result<(EntityId, String), SolveRefusal> {
        let [part, signal] = as_list(req(&kw, key, &scope)?, key)? else {
            return Err(err(format!("{scope} :{key} is (\"part\" \"port\")")));
        };
        let part = as_str(part, key)?;
        let entity = parts.get(part).copied();
        let entity = entity.ok_or_else(|| err(format!("{scope} unknown part `{part}`")))?;
        Ok((entity, as_str(signal, key)?.to_string()))
    };
    let law = as_list(req(&kw, "law", &scope)?, "law")?;
    let law = match as_sym(&law[0], "law")? {
        "pid" => {
            let gains = kwargs(&law[1..], &["kp", "ki", "kd", "setpoint"], &scope)?;
            let q = |k: &str| as_qty(req(&gains, k, &scope)?, k);
            let (kp, ki, kd, setpoint) = (q("kp")?, q("ki")?, q("kd")?, q("setpoint")?);
            ControllerLaw::Pid {
                kp,
                ki,
                kd,
                setpoint,
            }
        }
        o => {
            return Err(err(format!(
                "{scope} law `{o}` does not parse; declare (pid ...)"
            )));
        }
    };
    Ok(ControllerDecl {
        name,
        sensor_ports: vec![port("sensor")?],
        actuator_port: port("actuator")?,
        sample_rate: as_qty_d(req(&kw, "sample-rate", &scope)?, FREQUENCY, "sample-rate")?,
        law,
    })
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
            "position",
            "pattern",
        ],
        &scope,
    )?;
    let sample_rate = as_qty_d(req(&kw, "sample-rate", &scope)?, FREQUENCY, "sample-rate")?;
    let zero = [QtyAny::new(0.0, Dims([1, 0, 0, 0, 0, 0])); 3];
    let pos = kw
        .get("position")
        .map(|n| as_qty3(n, "position"))
        .transpose()?
        .unwrap_or(zero);
    let target = || -> Result<EntityId, SolveRefusal> {
        let tn = as_str(req(&kw, "target", &scope)?, "target")?;
        parts
            .get(tn)
            .copied()
            .ok_or_else(|| err(format!("{scope} unknown target `{tn}`")))
    };
    let kind_sym = as_sym(req(&kw, "kind", &scope)?, "probe.kind")?;
    if kind_sym != "microphone" && kw.contains_key("pattern") {
        return Err(err(format!("{scope} :pattern requires :kind microphone")));
    }
    let kind = match kind_sym {
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
        "field-quantity" => ProbeKind::FieldQuantity {
            target: target()?,
            position: zero,
            quantity: ParameterKey::new(
                as_str(req(&kw, "quantity", &scope)?, "quantity")?.to_string(),
                parse_dims(req(&kw, "dims", &scope)?)?,
            ),
        },
        "strain-gauge" => ProbeKind::StrainGauge {
            target: target()?,
            position: zero,
            gauge_direction: kw
                .get("direction")
                .map_or(Ok([1.0, 0.0, 0.0]), |n| as_vec3(n, "direction"))?,
        },
        "microphone" => ProbeKind::Microphone {
            target: target()?,
            position: pos,
            aim_direction: kw
                .get("direction")
                .map_or(Ok([0.0, 0.0, -1.0]), |n| as_vec3(n, "direction"))?,
            pattern: parse_mic_pattern(kw.get("pattern").copied())?,
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

fn parse_dims(node: &Node) -> Result<Dims, SolveRefusal> {
    let dn = as_list(node, "dims")?;
    if dn.len() != 6 {
        return Err(err("probe :dims requires 6 integers"));
    }
    let mut exps = [0i8; 6];
    for (slot, n) in exps.iter_mut().zip(dn) {
        let NodeKind::Int(i) = n.kind else {
            return Err(err("dims must be ints"));
        };
        *slot = i8::try_from(i).map_err(|_| err("dims exponent out of i8 range"))?;
    }
    Ok(Dims(exps))
}

fn parse_mic_pattern(node: Option<&Node>) -> Result<MicrophonePattern, SolveRefusal> {
    match node.map(|n| as_sym(n, "pattern")).transpose()? {
        None | Some("omnidirectional") => Ok(MicrophonePattern::Omnidirectional),
        Some("cardioid") => Ok(MicrophonePattern::Cardioid),
        Some("figure-eight") => Ok(MicrophonePattern::FigureEight),
        Some("supercardioid") => Ok(MicrophonePattern::Supercardioid),
        Some(o) => Err(err(format!("unknown microphone pattern `{o}`"))),
    }
}

fn parse_envelope(items: &[Node], out: &mut DeclaredEnvelopes) -> Result<(), SolveRefusal> {
    let k = as_sym(
        items.get(1).ok_or_else(|| err("envelope missing kind"))?,
        "envelope.kind",
    )?;
    let sc = format!("envelope `{k}`");
    let src = |kw: &BTreeMap<&str, &Node>| {
        let s = as_str(req(kw, "source", &sc)?, "source")?.trim();
        (!s.is_empty())
            .then(|| s.to_string())
            .ok_or_else(|| err(format!("{sc} requires non-empty :source")))
    };
    if k == "geometry" {
        let kw = kwargs(&items[2..], &["part", "axis", "half-width", "source"], &sc)?;
        let part = as_str(req(&kw, "part", &sc)?, "part")?.to_string();
        let ax_s = as_sym(req(&kw, "axis", &sc)?, "axis")?;
        let axis = BOX_AXES
            .iter()
            .position(|&a| a == ax_s)
            .ok_or_else(|| err(format!("{sc} unsupported box axis `{ax_s}`")))?;
        let hw = as_qty_d(req(&kw, "half-width", &sc)?, LENGTH, "half-width")?.value;
        if hw <= 0.0
            || out
                .geometry
                .iter()
                .any(|g| g.part == part && g.axis == axis)
        {
            return Err(err(format!("{sc} invalid or duplicate ({part}, {ax_s})")));
        }
        out.geometry.push(GeomEnvelope {
            part,
            axis,
            half_width_m: hw,
            source: src(&kw)?,
        });
        return Ok(());
    }
    let (key, list, zero_ok, cov_ok) = match k {
        "boundary-conditions" => ("target", &mut out.boundary_conditions, false, false),
        "measurement" => ("probe", &mut out.measurement, true, false),
        "model-form" => ("probe", &mut out.model_form, false, true),
        o => return Err(err(format!("unsupported envelope kind `{o}`"))),
    };
    let ok_kw: &[&str] = if cov_ok {
        &[key, "half-width", "source", "covered-by-source"]
    } else {
        &[key, "half-width", "source"]
    };
    let kw = kwargs(&items[2..], ok_kw, &sc)?;
    let target = as_str(req(&kw, key, &sc)?, key)?.to_string();
    if list.iter().any(|d| d.target == target) {
        return Err(err(format!("{sc} duplicate `{target}`")));
    }
    if let Some(cov) = kw.get("covered-by-source") {
        if kw.contains_key("half-width") || kw.contains_key("source") {
            return Err(err(format!(
                "{sc}: :covered-by-source cannot be combined with :half-width or :source"
            )));
        }
        let cite = as_str(cov, "covered-by-source")?.trim();
        if cite.is_empty() {
            return Err(err(format!("{sc} requires non-empty :covered-by-source")));
        }
        list.push(TargetEnvelope {
            target,
            half_width: None,
            source: cite.to_string(),
        });
        return Ok(());
    }
    let hw = as_qty(req(&kw, "half-width", &sc)?, "half-width")?;
    if hw.value < 0.0 || (!zero_ok && hw.value == 0.0) {
        return Err(err(format!("{sc} invalid :half-width for `{target}`")));
    }
    list.push(TargetEnvelope {
        target,
        half_width: Some(hw),
        source: src(&kw)?,
    });
    Ok(())
}

fn bind_envelopes(
    envelopes: &DeclaredEnvelopes,
    parts: &mut [DraftPartDecl],
    sources: &[SourceDecl],
    probes: &[ProbeDecl],
) -> Result<(), SolveRefusal> {
    for g in &envelopes.geometry {
        let p = parts
            .iter_mut()
            .find(|p| p.name == g.part)
            .ok_or_else(|| err(format!("geometry envelope unknown part `{}`", g.part)))?;
        let dim = match &p.geometry {
            DraftGeometrySource::Primitive(PrimitiveShape::Box { lx, ly, lz }) => {
                [lx, ly, lz][g.axis].value
            }
            other => {
                let sh = match other {
                    DraftGeometrySource::Primitive(PrimitiveShape::Cylinder { .. }) => "cylinder",
                    DraftGeometrySource::Primitive(PrimitiveShape::Sphere { .. }) => "sphere",
                    DraftGeometrySource::Primitive(PrimitiveShape::Tube { .. }) => "tube",
                    DraftGeometrySource::TensionedDisk(_) => "tensioned-disk",
                    _ => "non-box",
                };
                return Err(err(format!(
                    "geometry envelope `{}`: box only, not `{sh}`",
                    g.part
                )));
            }
        };
        let ax = BOX_AXES[g.axis];
        if g.half_width_m >= dim {
            return Err(err(format!(
                "geometry envelope `{}` {ax} half-width must be < nominal",
                g.part
            )));
        }
        p.uncertainty.push(DimensionBound {
            property: ParameterKey::new(ax, LENGTH),
            lower: QtyAny::new(dim - g.half_width_m, LENGTH),
            upper: QtyAny::new(dim + g.half_width_m, LENGTH),
            source: ToleranceSource::Standard {
                clause: g.source.clone(),
            },
        });
    }
    for b in &envelopes.boundary_conditions {
        let s = sources
            .iter()
            .find(|s| s.name == b.target)
            .ok_or_else(|| err(format!("bc envelope unknown source `{}`", b.target)))?;
        let (TimeHistory::Constant(q) | TimeHistory::Harmonic { offset: q, .. }) = &s.history
        else {
            return Err(err(format!(
                "bc envelope `{}`: constant/harmonic sources only",
                b.target
            )));
        };
        if b.half_width.is_none_or(|hw| hw.dims != q.dims) {
            return Err(err(format!(
                "bc envelope `{}` dims do not match source",
                b.target
            )));
        }
    }
    for d in envelopes.measurement.iter().chain(&envelopes.model_form) {
        if !probes.iter().any(|p| p.name == d.target) {
            return Err(err(format!("envelope unknown probe `{}`", d.target)));
        }
    }
    Ok(())
}

/// Parse, admit, and resolve card packs for a `.fsim` s-expression or JSON source.
#[allow(clippy::too_many_lines)]
pub fn admit_assembly_source(
    source: &[u8],
    packs: &CardPackSet,
) -> Result<AdmittedAssembly, SolveRefusal> {
    let text = std::str::from_utf8(source).map_err(|e| err(format!("non-UTF-8 source: {e}")))?;
    let ast = if text.trim_start().starts_with(['{', '[']) {
        crate::json::parse(text)
    } else {
        crate::sexpr::parse(text)
    }
    .map_err(|e| err(format!("assembly parse failed: {e}")))?;
    let canonical_source = crate::sexpr::print_checked(&ast)
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
    let mut by_interface = BTreeMap::new();
    let mut controllers = Vec::new();
    let (mut parts, mut joints, mut interfaces, mut sources, mut probes, mut requirements) = (
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    );
    let mut req_ids = BTreeSet::new();
    let mut envelopes = DeclaredEnvelopes::default();

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
                        &["ambient-temp", "ambient-pressure", "relative-humidity"],
                        "environment",
                    )?;
                    if let Some(t) = kw.get("ambient-temp") {
                        e.base.ambient_temperature = as_qty_d(t, TEMPERATURE, "ambient-temp")?;
                    }
                    if let Some(p) = kw.get("ambient-pressure") {
                        e.base.ambient_pressure = as_qty_d(p, PRESSURE, "ambient-pressure")?;
                    }
                    if let Some(rh) = kw.get("relative-humidity") {
                        e.relative_humidity = as_f64(rh, "relative-humidity")?;
                    }
                }
                env = Some(e);
            }
            "time-plan" => {
                let keys = ["horizon", "dt", "sample-rate", "scheme"];
                let kw = kwargs(&items[1..], &keys, "time-plan")?;
                let scheme = kw.get("scheme").map(|s| as_sym(s, "scheme")).transpose()?;
                let scheme = match scheme {
                    None | Some("explicit") => TimeScheme::Explicit,
                    Some("backward-euler") => TimeScheme::BackwardEuler,
                    Some("sdirk2") => TimeScheme::Sdirk2,
                    Some(o) => return Err(err(format!("unknown scheme `{o}`"))),
                };
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
                    scheme,
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
                        "query-at",
                        "physics",
                        "fidelity",
                        "initial",
                        "max-work-units",
                    ],
                    &scope,
                )?;
                if kw.contains_key("query-at") {
                    if kw.contains_key("query-temp") {
                        return Err(err(format!(
                            "{scope}: use :query-at or :query-temp, not both"
                        )));
                    }
                    if !kw.contains_key("material") {
                        return Err(err(format!("{scope}: :query-at requires :material")));
                    }
                }
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
                    "cavity" => DraftPartKind::Cavity,
                    "shell" => DraftPartKind::Shell,
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
                    let query_point = if let Some(at) = kw.get("query-at") {
                        parse_query_pairs(at, &scope, "query-at")?
                    } else {
                        let qt = kw
                            .get("query-temp")
                            .map_or(Ok(QtyAny::new(293.15, TEMPERATURE)), |n| {
                                as_qty_d(n, TEMPERATURE, "query-temp")
                            })?;
                        vec![(ParameterKey::new("T", TEMPERATURE), qt)]
                    };
                    materials.push(MaterialRef {
                        card,
                        query_point,
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
                let fidelity = parse_fidelity(req(&kw, "fidelity", &scope)?, &scope)?;
                // The work a part's owner may do in a step: a network's dense preparation.
                let positive = ":max-work-units is a positive integer";
                let work = match kw.get("max-work-units").map(|n| &n.kind) {
                    None => 1_000_000,
                    Some(NodeKind::Int(n)) if *n > 0 => n.unsigned_abs(),
                    Some(_) => return Err(err(format!("{scope}: {positive}"))),
                };
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
                        rungs: vec![fidelity],
                        error_budget: FidelityErrorBudget {
                            relative_tol: 1.0e-5,
                            absolute_tol: QtyAny::dimensionless(1.0e-12),
                            max_work_units: Some(work),
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
                        "fluid",
                        "gap",
                        "fluid-at",
                        "memory",
                        "sign",
                        "normal-force",
                        "heat-partition-a",
                        "heat-partition-effusivity",
                        "initial-slip",
                        "query-at",
                        "maximum-force",
                        "maximum-penetration",
                        "initial-closure",
                        "approach-speed",
                        "restitution",
                        "station-x",
                        "station-y",
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
                by_interface.insert(iname.clone(), entity);
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
                // A shear film's fluid, its gap `h` and the whole point `mu` is read at.
                let film = || -> Result<ThirdBodyDecl, SolveRefusal> {
                    let (fluid, gap) = (req(&kw, "fluid", &scope)?, req(&kw, "gap", &scope)?);
                    let query_point =
                        parse_query_pairs(req(&kw, "fluid-at", &scope)?, &scope, "fluid-at")?;
                    let fluid = as_str(fluid, "fluid")?;
                    Ok(ThirdBodyDecl {
                        material_card: resolve_card(packs, CardPackKind::Material, fluid, &scope)?,
                        film_thickness: as_qty_d(gap, LENGTH, "gap")?,
                        query_point,
                    })
                };
                let mut state = Vec::new();
                for (flag, key, dims) in [
                    ("area", "area", AREA),
                    ("normal-pressure", "normal_pressure", PRESSURE),
                    ("thickness", "thickness", LENGTH),
                    ("normal-force", "normal_force", FORCE),
                    ("initial-slip", "initial_slip", VELOCITY),
                ] {
                    if let Some(n) = kw.get(flag) {
                        state.push((ParameterKey::new(key, dims), as_qty_d(n, dims, flag)?));
                    }
                }
                for (flag, key) in [
                    ("fraction", "fraction"),
                    ("sign", "sign"),
                    ("heat-partition-a", "heat_partition_a"),
                    ("heat-partition-effusivity", "heat_partition_effusivity"),
                ] {
                    if let Some(n) = kw.get(flag) {
                        let val = QtyAny::dimensionless(as_f64(n, flag)?);
                        state.push((ParameterKey::new(key, Dims::NONE), val));
                    }
                }
                if let Some(n) = kw.get("memory") {
                    state.extend(parse_memory(n, &scope)?);
                }
                if let Some(n) = kw.get("query-at") {
                    for (k, q) in parse_query_pairs(n, &scope, "query-at")? {
                        if state.iter().any(|(have, _)| have.name == k.name) {
                            return Err(err(format!("{scope} repeats `{}` in :query-at", k.name)));
                        }
                        state.push((k, q));
                    }
                }
                // A Hertz contact's ceilings, gap, initial closure and Hunt-Crossley loss; on
                // every other kind these are refused by name, but `:gap` is a shear film's own.
                let kind_sym = as_sym(req(&kw, "kind", &scope)?, "interface.kind")?;
                let hertz = kind_sym == "hertz-contact";
                for (flag, key, dims) in [
                    ("maximum-force", "maximum_force", FORCE),
                    ("maximum-penetration", "maximum_penetration", LENGTH),
                    ("gap", "gap", LENGTH),
                    ("initial-closure", "initial_closure", LENGTH),
                    ("approach-speed", "approach_speed", VELOCITY),
                    ("restitution", "restitution", Dims::NONE),
                    ("station-x", "station_x", LENGTH),
                    ("station-y", "station_y", LENGTH),
                ] {
                    let Some(n) = kw.get(flag) else {
                        continue;
                    };
                    if !hertz && flag != "gap" {
                        return Err(err(format!("{scope}: :{flag} needs a hertz contact")));
                    }
                    if hertz {
                        let q = if dims == Dims::NONE {
                            QtyAny::dimensionless(as_f64(n, flag)?)
                        } else {
                            as_qty_d(n, dims, flag)?
                        };
                        state.push((ParameterKey::new(key, dims), q));
                    }
                }
                let bare = "a temperature source takes no card or state";
                let kind = match as_sym(req(&kw, "kind", &scope)?, "interface.kind")? {
                    // A heat route is bookkeeping, not a material interface: no card.
                    "heat-route" if kw.contains_key("card") || kw.contains_key("family") => {
                        return Err(err(format!("{scope}: a heat route takes no card")));
                    }
                    "heat-route" => DraftInterfaceKind::HeatRoute { state },
                    "piston-face" if kw.contains_key("card") || kw.contains_key("family") => {
                        return Err(err(format!("{scope}: a piston face takes no card")));
                    }
                    "piston-face" => DraftInterfaceKind::PistonFace { state },
                    "shear-film" if kw.contains_key("card") || kw.contains_key("family") => {
                        return Err(err(format!("{scope}: a film's one card is its fluid")));
                    }
                    "shear-film" => DraftInterfaceKind::ShearFilm {
                        third_body: film()?,
                        state,
                    },
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
                    "hertz-contact" if kw.contains_key("card") || kw.contains_key("family") => {
                        return Err(err(format!("{scope}: a Hertz contact takes no card")));
                    }
                    "hertz-contact" => DraftInterfaceKind::HertzContact { state },
                    "glue" => DraftInterfaceKind::Glue { law_card: card()? },
                    // Its bushing joins the body `:a`; it reads the cell `:b`.
                    "temperature-source" if kw.len() > 3 => {
                        return Err(err(format!("{scope}: {bare}")));
                    }
                    "temperature-source" => DraftInterfaceKind::TemperatureSource,
                    o => return Err(err(format!("{scope} unsupported kind `{o}`"))),
                };
                let contact = matches!(kind, DraftInterfaceKind::MechanicalContact { .. });
                for flag in [
                    "memory",
                    "normal-force",
                    "heat-partition-a",
                    "heat-partition-effusivity",
                    "initial-slip",
                    "query-at",
                ] {
                    if kw.contains_key(flag) && !contact {
                        return Err(err(format!("{scope}: :{flag} needs a mechanical contact")));
                    }
                }
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
                let ok = ["part", "interface", "port", "history"];
                let kw = kwargs(&items[2..], &ok, &scope)?;
                // A part's port, or an interface's, named earlier in the file.
                let (names, what, named) = match (kw.get("part"), kw.get("interface")) {
                    (Some(named), None) => (&by_name, "part", named),
                    (None, Some(named)) => (&by_interface, "interface", named),
                    _ => return Err(err(format!("{scope} names one :part or one :interface"))),
                };
                let entity = names.get(as_str(named, what)?);
                let entity = *entity.ok_or_else(|| err(format!("{scope} unknown {what}")))?;
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
            "controller" => controllers.push(parse_controller(items, &by_name)?),
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
                let kw = kwargs(
                    &items[2..],
                    &["analysis", "probe", "stat", "relation", "limit"],
                    &sc,
                )?;
                let analysis = kw
                    .get("analysis")
                    .map(|n| {
                        let s = as_sym(n, "analysis")?;
                        AssemblyAnalysisMode::parse(s)
                            .ok_or_else(|| err(format!("{sc} unsupported analysis `{s}`")))
                    })
                    .transpose()?;
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
                let lim_n = req(&kw, "limit", &sc)?;
                let limit = as_qty(lim_n, "limit")
                    .or_else(|_| as_f64(lim_n, "limit").map(QtyAny::dimensionless))?;
                requirements.push(AssemblyRequirement {
                    id,
                    analysis,
                    probe,
                    stat,
                    relation,
                    limit,
                });
            }
            "envelope" => parse_envelope(items, &mut envelopes)?,
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
    bind_envelopes(&envelopes, &mut parts, &sources, &probes)?;
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
            controllers,
            probes,
        },
        &scenario,
    )
    .map_err(|e| err(format!("structural admission refused: {e}")))?;
    for r in &requirements {
        let ok = match r.analysis {
            Some(AssemblyAnalysisMode::Modal) => matches!(
                r.probe.as_str(),
                "modal.natural_rad_s" | "modal.frequency_hz" | "modal.damping_ratio"
            ),
            Some(AssemblyAnalysisMode::Harmonic) => r.probe == "harmonic.transfer_magnitude",
            Some(AssemblyAnalysisMode::Stationary) => {
                r.probe == "stationary.residual_norm"
                    || r.probe
                        .strip_prefix("stationary.state.")
                        .is_some_and(|s| !s.is_empty())
                    || draft.input().probes.iter().any(|p| p.name == r.probe)
            }
            None | Some(AssemblyAnalysisMode::Transient) => {
                draft.input().probes.iter().any(|p| p.name == r.probe)
            }
        };
        if !ok {
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
        envelopes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);

    fn admit(attrs: &str) -> Result<AdmittedAssembly, SolveRefusal> {
        admit_with(attrs, &[])
    }

    fn admit_with(attrs: &str, raw: &[&[u8]]) -> Result<AdmittedAssembly, SolveRefusal> {
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
        let packs = PackCatalog::from_raw_bytes(raw).unwrap();
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
        let packs = PackCatalog::default();
        let a1 = admit(full).expect("full :initial admits");
        let init1 = a1.draft.input().parts[0].initial.as_ref().unwrap();
        assert_eq!(init1, &expected);
        assert_eq!(init1.angular_velocity[2].dims, FREQUENCY);

        let a2 = admit_assembly_source(&a1.canonical_source, &packs).unwrap();
        assert_eq!(a2.draft.input().parts[0].initial.as_ref(), Some(&expected));

        let ast = crate::sexpr::parse(std::str::from_utf8(&a1.canonical_source).unwrap()).unwrap();
        let json = crate::json::print_checked(&ast).unwrap();
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
    fn a_time_plan_is_explicit_unless_it_names_an_implicit_scheme() {
        let packs = PackCatalog::default();
        let plain = admit("").unwrap();
        assert_eq!(plain.draft.input().time_plan.scheme, TimeScheme::Explicit);
        let source = std::str::from_utf8(&plain.canonical_source).unwrap();
        let scheme = |name: &str| {
            let named = format!(":scheme {name} :sample-rate");
            let src = source.replace(":sample-rate", &named);
            let admitted = admit_assembly_source(src.as_bytes(), &packs);
            admitted.map(|a| a.draft.input().time_plan.scheme)
        };
        assert_eq!(scheme("backward-euler").unwrap(), TimeScheme::BackwardEuler);
        assert_eq!(scheme("sdirk2").unwrap(), TimeScheme::Sdirk2);
        assert!(scheme("rk4").is_err());
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
            assert_eq!(r.code, "cli-assembly-admit");
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
            admit_assembly_source(src.as_bytes(), &PackCatalog::default())
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

    /// `admit_with`'s text closing the rotor's part, a thermal cell, and
    /// `forms`, whose last parenthesis `admit_with` writes.
    fn with_cell(forms: &str, raw: &[&[u8]]) -> Result<AdmittedAssembly, SolveRefusal> {
        let forms = forms.strip_suffix(')').expect("closed forms");
        let cell = "(part \"cell\" :kind solid :geometry (box 0.1m 0.1m 0.1m) \
                    :physics (thermal) :fidelity lumped)";
        admit_with(&format!(") {cell} {forms}"), raw)
    }

    #[test]
    fn memory_a_temperature_source_and_a_source_on_an_interface_admit_and_round_trip() {
        // Any interface card shows the memory; admission reads no bushing law.
        let bytes: &[u8] = include_bytes!("../../../examples/contact-pair/cold-hot.fsintpk");
        let packs = PackCatalog::from_raw_bytes(&[bytes]).unwrap();
        let law_card = packs.interfaces()[0].card;
        let pad = "(part \"pad\" :kind rigid-body :geometry (box 0.1m 0.1m 0.1m) \
                   :physics (mechanics) :fidelity lumped)";
        let memory = "displacement-0-x 0.4m displacement-0-y 0m displacement-0-z 0m \
                      displacement-0-rx 0rad displacement-0-ry 0rad displacement-1-rz 0.1rad \
                      displacement-eq-x 0.2m";
        let mount = format!(
            "(interface \"mount\" :kind mechanical-contact :a \"rotor\" :b \"pad\" \
             :card \"{}\" :family \"relaxation-bushing\" :memory ({memory}))",
            law_card.to_hex()
        );
        let warm = "(interface \"warm\" :kind temperature-source :a \"rotor\" :b \"cell\")";
        // Admission knows no ports: the lowering owns them.
        let hold = "(source \"hold\" :interface \"mount\" :port \"drive_x\" \
                    :history (constant 0.7N))";
        let forms = format!("{pad} {mount} {warm} {hold}");
        let a = with_cell(&forms, &[bytes]).expect("all three admit");
        let input = a.draft.input();
        let key = |name, value, dims| (ParameterKey::new(name, dims), QtyAny::new(value, dims));
        let state = vec![
            key("displacement_0_x", 0.4, LENGTH),
            key("displacement_0_y", 0.0, LENGTH),
            key("displacement_0_z", 0.0, LENGTH),
            key("displacement_0_rx", 0.0, Dims::NONE),
            key("displacement_0_ry", 0.0, Dims::NONE),
            key("displacement_1_rz", 0.1, Dims::NONE),
            key("displacement_eq_x", 0.2, LENGTH),
        ];
        let third_body = None;
        let contact = DraftInterfaceKind::MechanicalContact {
            law_card,
            third_body,
            state,
        };
        assert_eq!(input.interfaces[0].kind, contact);
        let source = &DraftInterfaceKind::TemperatureSource;
        assert_eq!(&input.interfaces[1].kind, source);
        let (entity, port) = (input.interfaces[0].entity, "drive_x".to_owned());
        assert_eq!(input.sources[0].target, SourceTarget::Port { entity, port });
        let again = admit_assembly_source(&a.canonical_source, &packs).unwrap();
        let ast = crate::sexpr::parse(std::str::from_utf8(&a.canonical_source).unwrap()).unwrap();
        let json = crate::json::print_checked(&ast).unwrap();
        let from_json = admit_assembly_source(json.as_bytes(), &packs).unwrap();
        for b in [again, from_json] {
            assert_eq!(b.draft.input().interfaces, input.interfaces);
            assert_eq!(b.draft.input().sources, input.sources);
        }
    }

    #[test]
    fn a_part_states_its_work_cap_as_a_positive_integer() {
        let a = admit(" :max-work-units 340023768").expect("a positive cap admits");
        let budget = &a.draft.input().parts[0].fidelity.error_budget;
        assert_eq!(budget.max_work_units, Some(340_023_768));
        for bad in ["0", "-1", "1.5e9"] {
            let why = admit(&format!(" :max-work-units {bad}"))
                .err()
                .unwrap()
                .what;
            assert!(
                why.contains(":max-work-units is a positive integer"),
                "{bad}: {why}"
            );
        }
    }

    #[test]
    fn a_pid_controller_admits_round_trips_and_refuses_by_name() {
        let law = "(pid :kp 72N/m :ki 36N/m/s :kd 38.4N*s/m :setpoint 0.05m)";
        let form = |sensor: &str, law: &str| {
            format!(
                "(controller \"servo\" :sensor ({sensor}) :actuator (\"rotor\" \"mechanics\") \
                 :sample-rate 16Hz :law {law})"
            )
        };
        let a1 = with_cell(&form("\"rotor\" \"position\"", law), &[]).unwrap();
        let input = a1.draft.input();
        let rotor = input
            .parts
            .iter()
            .find(|p| p.name == "rotor")
            .unwrap()
            .entity;
        let q = QtyAny::new;
        let expected = ControllerDecl {
            name: "servo".to_owned(),
            sensor_ports: vec![(rotor, "position".to_owned())],
            actuator_port: (rotor, "mechanics".to_owned()),
            sample_rate: q(16.0, FREQUENCY),
            law: ControllerLaw::Pid {
                kp: q(72.0, Dims([0, 1, -2, 0, 0, 0])),
                ki: q(36.0, Dims([0, 1, -3, 0, 0, 0])),
                kd: q(38.4, Dims([0, 1, -1, 0, 0, 0])),
                setpoint: q(0.05, LENGTH),
            },
        };
        assert_eq!(input.controllers, std::slice::from_ref(&expected));
        let a2 = admit_assembly_source(&a1.canonical_source, &PackCatalog::default()).unwrap();
        assert_eq!(a2.draft.input().controllers, [expected]);

        let refused = |forms: &str| with_cell(forms, &[]).err().unwrap().what;
        for (sensor, law, named) in [
            (
                "\"rotor\" \"position\"",
                "(table)",
                "law `table` does not parse",
            ),
            ("\"nowhere\" \"position\"", law, "unknown part `nowhere`"),
            ("\"rotor\"", law, ":sensor is (\"part\" \"port\")"),
            (
                "\"rotor\" \"position\"",
                "(pid :kp 72N/m :ki 36N/m/s :kd 38.4N*s/m)",
                "missing `:setpoint`",
            ),
        ] {
            let why = refused(&form(sensor, law));
            assert!(why.contains(named), "{named}: {why}");
        }
    }

    #[test]
    fn memory_keys_temperature_sources_and_interface_sources_refuse_by_name() {
        let refused = |forms: &str| with_cell(forms, &[]).err().unwrap().what;
        let route = "(interface \"route\" :kind heat-route :a \"rotor\" :b \"cell\" :fraction 0.7";
        for row in [
            "displacement-0-w 0m | no :memory key `displacement-0-w`",
            "displacement-01-x 0m | no :memory key `displacement-01-x`",
            "displacement-Eq-x 0m | no :memory key `displacement-Eq-x`",
            "displacement-0-x 0.1s | displacement-0-x has wrong SI dimensions",
            "displacement-0-x 0m displacement-0-x 0m | repeats `displacement-0-x`",
            "displacement-0-x 0m | :memory needs a mechanical contact",
        ] {
            let (memory, named) = row.split_once(" | ").unwrap();
            let why = refused(&format!("{route} :memory ({memory}))"));
            assert!(why.contains(named), "{memory}: {why}");
        }
        let warm =
            "(interface \"warm\" :kind temperature-source :a \"rotor\" :b \"cell\" :card \"x\")";
        let (why, bare) = (refused(warm), "a temperature source takes no card or state");
        assert!(why.contains(bare), "{why}");
        for row in [
            ":part \"rotor\" :interface \"x\" | names one :part or one :interface",
            " | names one :part or one :interface",
            ":interface \"nowhere\" | unknown interface",
        ] {
            let (n, named) = row.split_once(" | ").unwrap();
            let source = format!("(source \"s\" {n} :port \"p\" :history (constant 1N))");
            let why = refused(&source);
            assert!(why.contains(named), "{n}: {why}");
        }
    }

    #[test]
    fn shear_film_names_its_fluid_gap_and_query_point_and_no_card() {
        // Any material card shows the binding; the lowering reads its viscosity.
        let bytes: &[u8] = include_bytes!(
            "../../../examples/assemblies/packs/aluminum-pure-thermal-nbs-c447.fsmcdpk"
        );
        let film = |extra: &str| {
            let pad = "(part \"pad\" :kind rigid-body :geometry (box 0.1m 0.1m 0.1m) \
                 :physics (mechanics) :fidelity lumped)";
            let film = "(interface \"film\" :kind shear-film :a \"rotor\" :b \"pad\" \
                 :fluid \"aluminum-pure-nbs-c447\" :gap 0.0001m";
            let at = ":fluid-at (T 313.15K phase-liquid 1)";
            admit_with(
                &format!(") {pad} {film} {at} :area 0.01m^2{extra}"),
                &[bytes],
            )
        };
        let a = film("").expect("a shear film admits with its fluid");
        let kind = &a.draft.input().interfaces[0].kind;
        let DraftInterfaceKind::ShearFilm { third_body, state } = kind else {
            panic!("{kind:?}")
        };
        let fluid = a.card_store.material(&third_body.material_card).unwrap();
        assert_eq!(fluid.id().chemistry, "aluminum-pure-nbs-c447");
        assert_eq!(third_body.film_thickness, QtyAny::new(0.0001, LENGTH));
        let key = |name, value, dims| (ParameterKey::new(name, dims), QtyAny::new(value, dims));
        let t = key("T", 313.15, TEMPERATURE);
        let liquid = key("phase-liquid", 1.0, Dims::NONE);
        assert_eq!(third_body.query_point, vec![t, liquid]);
        assert_eq!(state, &vec![key("area", 0.01, AREA)]);
        let why = film(" :card \"x\"").err().unwrap().what;
        assert!(why.contains("a film's one card is its fluid"), "{why}");
    }

    #[test]
    fn envelope_clauses_bind_uncertainty_and_refuse_non_box_or_bad_targets() {
        let with = |extra: &str| {
            admit(&format!(
                ") (part \"c\" :kind solid :geometry (cylinder 0.02m 0.05m) :physics (thermal) :fidelity lumped) \
                 (source \"h\" :part \"rotor\" :port \"thermal\" :history (constant 10W)) {extra} \
                 (probe \"tp\" :kind thermocouple :target \"rotor\" :sample-rate 1000Hz"
            ))
        };
        let ok = with(
            "(envelope geometry :part \"rotor\" :axis lx :half-width 0.001m :source \"ISO 2768-m\") \
             (envelope boundary-conditions :target \"h\" :half-width 0.5W :source \"PSU 5%\") \
             (envelope measurement :probe \"tp\" :half-width 0.05K :source \"SPRT\") \
             (envelope model-form :probe \"tp\" :half-width 0.1K :source \"Biot<0.01\")",
        )
        .unwrap();
        assert_eq!(
            (
                ok.draft.input().parts[0].uncertainty.len(),
                ok.envelopes.boundary_conditions.len()
            ),
            (1, 1)
        );
        let chk = |bad: &str, needle: &str| {
            let e = with(bad).err().unwrap();
            assert!(e.what.contains(needle), "`{}` missing `{needle}`", e.what);
        };
        chk(
            "(envelope geometry :part \"c\" :axis lx :half-width 0.001m :source \"s\")",
            "cylinder",
        );
        chk(
            "(envelope geometry :part \"rotor\" :axis lx :half-width 0.1m :source \"s\")",
            "< nominal",
        );
        chk(
            "(envelope boundary-conditions :target \"h\" :half-width 1K :source \"s\")",
            "dims do not match",
        );
        chk(
            "(envelope model-form :probe \"tp\" :half-width 0K :source \"s\")",
            "invalid :half-width",
        );
        chk(
            "(envelope measurement :probe \"missing\" :half-width 0.1K :source \"s\")",
            "unknown probe",
        );
        chk(
            "(envelope measurement :probe \"tp\" :half-width 0.1K :source \" \")",
            "non-empty :source",
        );
        let cov =
            with("(envelope model-form :probe \"tp\" :covered-by-source \"NIST SP 260-175\")")
                .unwrap();
        assert!(cov.envelopes.model_form[0].half_width.is_none());
        chk(
            "(envelope model-form :probe \"tp\" :covered-by-source \"s\" :half-width 0.1K)",
            "cannot be combined",
        );
        chk(
            "(envelope model-form :probe \"tp\" :covered-by-source \" \")",
            "non-empty :covered-by-source",
        );
        chk(
            "(envelope measurement :probe \"tp\" :covered-by-source \"s\")",
            "rejects unknown keyword `:covered-by-source`",
        );
    }

    #[test]
    fn part_query_at_and_one_d_field_fidelity_admit_and_refuse() {
        let b: &[u8] =
            include_bytes!("../../../examples/assemblies/packs/eps-nist-srm1453-2013.fsmcdpk");
        let slab = |e: &str, f: &str| {
            let s = format!(
                ") (part \"s\" :kind solid :geometry (box 0.0254m 0.25m 0.25m) :material \"eps-nist-srm1453-2013\" {e} :physics (thermal) :fidelity {f}"
            );
            admit_with(&s, &[b])
        };
        let ok = slab(
            ":query-at (mean-temperature 297.15K bulk-density 41kg/m^3 source-lot-srm1453 1)",
            "(one-d-field 1)",
        )
        .unwrap();
        let p = &ok.draft.input().parts[1];
        assert_eq!(p.fidelity.rungs, [FidelityRung::OneDField { cells: 1 }]);
        assert_eq!(p.materials[0].query_point.len(), 3);
        for (e, f, n) in [
            (":query-at (T 300K) :query-temp 300K", "lumped", "not both"),
            (":query-at (T)", "lumped", "pairs keys"),
            (":query-at (T 300K)", "(one-d-field 2)", "(one-d-field 2)"),
            (":query-at (T 300K)", "modal", "`modal`"),
        ] {
            assert!(slab(e, f).err().unwrap().what.contains(n));
        }
        let reqs = admit(
            ") (requirement \"m\" :analysis modal :probe \"modal.frequency_hz\" :stat min :relation at-least :limit 1Hz) \
               (requirement \"h\" :analysis harmonic :probe \"harmonic.transfer_magnitude\" :stat max :relation at-most :limit 0.5s) \
               (requirement \"s\" :analysis stationary :probe \"stationary.state.x\" :stat final :relation at-least :limit -1m",
        )
        .unwrap();
        assert_eq!(reqs.requirements.len(), 3);
        assert!(admit(") (requirement \"b\" :analysis bogus :probe \"modal.frequency_hz\" :stat min :relation at-least :limit 1Hz").is_err());
        assert!(admit(") (requirement \"b\" :analysis modal :probe \"mass_pos\" :stat min :relation at-least :limit 1Hz").is_err());
    }

    const VALID_TENSIONED_DISK: &str = "(tensioned-disk :radius 0.1778m :thickness 0.00025m \
:tension 2500N/m :radial-intervals 8 :azimuths 24)";

    fn admit_drum_head(geom: &str, fid: &str) -> Result<AdmittedAssembly, SolveRefusal> {
        let steel: &[u8] = include_bytes!(
            "../../../examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk"
        );
        let packs = PackCatalog::from_raw_bytes(&[steel]).unwrap();
        let src = format!(
            "(assembly \"drum_head\" (version 1) (units \"SI\") (seed 0x1) \
             (versions (constellation \"2026-07\") (workspace \"test\")) \
             (budget (wall 60s) (mem 64MiB)) \
             (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*)) \
             (reference-state assembled) (environment earth-lab-air) \
             (time-plan :horizon 0.01s :dt 0.0001s :sample-rate 10000Hz) \
             (part \"head\" :kind shell :geometry {geom} \
              :material \"steel-music-wire-ligo-p1400229\" :query-temp 293.15K \
              :physics (mechanics) :fidelity {fid}))"
        );
        admit_assembly_source(src.as_bytes(), &packs)
    }

    #[test]
    fn shell_modal_tensioned_disk_admits_roundtrips_and_pins_missing_lowering() {
        let steel: &[u8] = include_bytes!(
            "../../../examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk"
        );
        let packs = PackCatalog::from_raw_bytes(&[steel]).unwrap();
        let ok = admit_drum_head(VALID_TENSIONED_DISK, "(modal 12)").unwrap();
        let part = &ok.draft.input().parts[0];
        assert_eq!(part.kind, DraftPartKind::Shell);
        assert_eq!(
            part.fidelity.rungs,
            [FidelityRung::Modal { mode_count: 12 }]
        );
        assert_eq!(
            part.geometry,
            DraftGeometrySource::TensionedDisk(TensionedDiskGenerator {
                radius: QtyAny::new(0.1778, LENGTH),
                thickness: QtyAny::new(0.00025, LENGTH),
                tension: QtyAny::new(2500.0, SURFACE_TENSION),
                radial_intervals: 8,
                azimuths: 24,
            })
        );
        // Round-trip through canonical s-expr and isomorphic JSON preserves geometry & modal rung.
        let ast = crate::sexpr::parse(std::str::from_utf8(&ok.canonical_source).unwrap()).unwrap();
        let json = crate::json::print_checked(&ast).unwrap();
        for src in [ok.canonical_source.as_slice(), json.as_bytes()] {
            let again = admit_assembly_source(src, &packs).unwrap();
            assert_eq!(again.draft.input().parts[0].geometry, part.geometry);
            assert_eq!(
                again.draft.input().parts[0].fidelity.rungs,
                part.fidelity.rungs
            );
        }
        // Condition 2: lowering refuses by name until .11.9.2 registers (Shell, {Mechanics}, Modal).
        let low_err = fs_system::LoweringRegistry::with_kernel_defaults()
            .lower(&ok.draft, &ok.card_store)
            .err()
            .expect("expected MissingPartConstructor refusal");
        let msg = low_err.to_string();
        assert!(
            matches!(
                low_err,
                fs_system::LoweringError::MissingPartConstructor {
                    ref part,
                    kind: DraftPartKind::Shell,
                    ref tried_rungs,
                    ..
                } if part == "head"
                    && tried_rungs == &[fs_system::lowering::FidelityRungKind::Modal]
            ),
            "unexpected lowering error: {msg}"
        );
        for tok in ["head", "shell", "mechanics", "modal"] {
            assert!(msg.contains(tok), "missing `{tok}` in `{msg}`");
        }
    }

    #[test]
    fn shell_modal_tensioned_disk_refuses_invalid_inputs() {
        for (g, f, needle) in [
            (
                "(tensioned-disk :radius 0.17m :thickness 0.00025m :tension 100N :radial-intervals 8 :azimuths 24)",
                "(modal 12)",
                "tensioned-disk.tension has wrong SI dimensions",
            ),
            (
                "(tensioned-disk :radius 0.17m :thickness 0.00025m :tension 100Pa :radial-intervals 8 :azimuths 24)",
                "(modal 12)",
                "tensioned-disk.tension has wrong SI dimensions",
            ),
            (
                "(tensioned-disk :thickness 0.00025m :tension 2500N/m :radial-intervals 8 :azimuths 24)",
                "(modal 12)",
                "missing `:radius`",
            ),
            (
                "(tensioned-disk :radius 0.17m :tension 2500N/m :radial-intervals 8 :azimuths 24)",
                "(modal 12)",
                "missing `:thickness`",
            ),
            (
                "(tensioned-disk :radius 0.17m :thickness 0.00025m :radial-intervals 8 :azimuths 24)",
                "(modal 12)",
                "missing `:tension`",
            ),
            (
                "(tensioned-disk :radius 0.17m :thickness 0.00025m :tension 2500N/m :azimuths 24)",
                "(modal 12)",
                "missing `:radial-intervals`",
            ),
            (
                "(tensioned-disk :radius 0.17m :thickness 0.00025m :tension 2500N/m :radial-intervals 8)",
                "(modal 12)",
                "missing `:azimuths`",
            ),
            (
                "(tensioned-disk :radius 0.17m :thickness 0.00025m :tension 0N/m :radial-intervals 8 :azimuths 24)",
                "(modal 12)",
                "tensioned_disk.tension",
            ),
            (
                "(tensioned-disk :radius 0.17m :thickness 0.00025m :tension 2500N/m :radial-intervals 0 :azimuths 24)",
                "(modal 12)",
                ":radial-intervals requires a positive integer",
            ),
            (
                "(tensioned-disk :radius 0.17m :thickness 0.00025m :tension 2500N/m :radial-intervals 8 :azimuths 2)",
                "(modal 12)",
                "azimuths must be at least 3",
            ),
            (VALID_TENSIONED_DISK, "(modal)", "mode_count"),
            (VALID_TENSIONED_DISK, "(modal 0)", "positive mode_count"),
            (VALID_TENSIONED_DISK, "(modal -3)", "positive mode_count"),
            (VALID_TENSIONED_DISK, "(modal 1.5)", "integer mode_count"),
            (VALID_TENSIONED_DISK, "modal", "`modal`"),
        ] {
            let err = admit_drum_head(g, f)
                .err()
                .unwrap_or_else(|| panic!("expected err for {g} / {f}"));
            assert!(
                err.what.contains(needle),
                "got `{}` expected `{needle}`",
                err.what
            );
        }
    }

    #[test]
    fn hertz_contact_station_admits_roundtrips_and_refuses() {
        let s: &[u8] = include_bytes!(
            "../../../examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk"
        );
        let packs = PackCatalog::from_raw_bytes(&[s]).unwrap();
        let m =
            ":material \"steel-music-wire-ligo-p1400229\" :query-temp 293.15K :physics (mechanics)";
        let run = |k: &str, ex: &str| {
            let src = format!(
                "(assembly \"s\" (version 1) (units \"SI\") (seed 0x1) \
                 (versions (constellation \"2026-07\") (workspace \"t\")) \
                 (budget (wall 60s) (mem 64MiB)) \
                 (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*)) \
                 (reference-state assembled) (environment earth-lab-air) \
                 (time-plan :horizon 0.01s :dt 0.0001s :sample-rate 10000Hz) \
                 (part \"st\" :kind rigid-body :geometry (sphere 0.005m) {m} :fidelity lumped) \
                 (part \"hd\" :kind shell :geometry {VALID_TENSIONED_DISK} {m} :fidelity (modal 3)) \
                 (interface \"hit\" :kind {k} :a \"st\" :b \"hd\" {ex}))"
            );
            admit_assembly_source(src.as_bytes(), &packs)
        };
        let f = ":maximum-force 10000N :maximum-penetration 0.01m";
        let st = format!("{f} :station-x 0.06m :station-y -0.02m");
        let ok = run("hertz-contact", &st).unwrap();
        let txt = std::str::from_utf8(&ok.canonical_source).unwrap();
        let ast = crate::sexpr::parse(txt).unwrap();
        let js = crate::json::print_checked(&ast).unwrap();
        for b in [ok.canonical_source.as_slice(), js.as_bytes()] {
            let a = admit_assembly_source(b, &packs).unwrap();
            let DraftInterfaceKind::HertzContact { state } = &a.draft.input().interfaces[0].kind
            else {
                panic!()
            };
            let v = |k| state.iter().find(|(p, _)| p.name == k).unwrap().1;
            assert_eq!(v("station_x"), QtyAny::new(0.06, LENGTH));
            assert_eq!(v("station_y"), QtyAny::new(-0.02, LENGTH));
        }
        for (k, ex, w) in [
            ("hertz-contact", ":station-x 0.06m", "pair required"),
            ("hertz-contact", ":station-y 0.02m", "pair required"),
            ("hertz-contact", ":station-x 1N :station-y 0m", "wrong SI"),
            ("piston-face", ":station-x 0.06m", "needs a hertz contact"),
        ] {
            let arg = if k == "piston-face" {
                ex.into()
            } else {
                format!("{f} {ex}")
            };
            assert!(run(k, &arg).err().unwrap().what.contains(w));
        }
    }

    #[test]
    fn microphone_probe_admits_roundtrips_and_refuses() {
        let s: &[u8] = include_bytes!(
            "../../../examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk"
        );
        let packs = PackCatalog::from_raw_bytes(&[s]).unwrap();
        let run = |pr: &str| {
            let src = format!(
                "(assembly \"m\" (version 1) (units \"SI\") (seed 0x1) \
                 (versions (constellation \"2026-07\") (workspace \"t\")) \
                 (budget (wall 60s) (mem 64MiB)) \
                 (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*)) \
                 (reference-state assembled) \
                 (environment :ambient-temp 293.15K :ambient-pressure 101325Pa :relative-humidity 0) \
                 (time-plan :horizon 0.01s :dt 0.0001s :sample-rate 10000Hz) \
                 (part \"hd\" :kind shell :geometry {VALID_TENSIONED_DISK} \
                   :material \"steel-music-wire-ligo-p1400229\" :query-temp 293.15K :physics (mechanics) :fidelity (modal 3)) \
                 {pr})"
            );
            admit_assembly_source(src.as_bytes(), &packs)
        };
        let ok = run("(probe \"mic\" :kind microphone :target \"hd\" :position (0.08m 0.05m 0.35m) :pattern omnidirectional :sample-rate 48000Hz)").unwrap();
        let ast = crate::sexpr::parse(std::str::from_utf8(&ok.canonical_source).unwrap()).unwrap();
        let js = crate::json::print_checked(&ast).unwrap();
        for b in [ok.canonical_source.as_slice(), js.as_bytes()] {
            let a = admit_assembly_source(b, &packs).unwrap();
            assert_eq!(
                a.draft.input().environment.relative_humidity.to_bits(),
                0.0_f64.to_bits()
            );
            let ProbeKind::Microphone {
                target,
                position,
                pattern,
                ..
            } = &a.draft.input().probes[0].kind
            else {
                panic!()
            };
            assert_eq!(*target, a.draft.input().parts[0].entity);
            assert_eq!(position[2], QtyAny::new(0.35, LENGTH));
            assert_eq!(*pattern, MicrophonePattern::Omnidirectional);
        }
        for (pr, w) in [
            (
                "(probe \"m\" :kind microphone :position (0m 0m 0.5m) :sample-rate 48000Hz)",
                "missing `:target`",
            ),
            (
                "(probe \"m\" :kind microphone :target \"nope\" :position (0m 0m 0.5m) :sample-rate 48000Hz)",
                "unknown target",
            ),
            (
                "(probe \"m\" :kind microphone :target \"hd\" :position (0m 0m 0.5m) :pattern laser :sample-rate 48000Hz)",
                "unknown microphone pattern",
            ),
            (
                "(probe \"m\" :kind microphone :target \"hd\" :position (0m 0m 1N) :sample-rate 48000Hz)",
                "dimension",
            ),
            (
                "(probe \"t\" :kind thermocouple :target \"hd\" :pattern cardioid :sample-rate 100Hz)",
                ":pattern requires :kind microphone",
            ),
        ] {
            assert!(run(pr).err().unwrap().what.contains(w));
        }
    }
}
