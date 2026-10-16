//! Score commands, independent of the physical interaction state.
//!
//! SMF notes select named assembly configurations. They never synthesize a
//! frequency. Velocity bins, controller meanings, tempo before the first tempo
//! event, dynamics, command bindings and observation prior are explicit inputs.
//! This v1 importer describes one monophonic performer and refuses overlap.
//! Timing uses the ceiling of the complete rational tempo integral (G0).
//! The frozen sidecar retains its MIDI bytes and hashes of all interpretation
//! inputs; decoding requires those same inputs, as with an assembly sidecar.
use crate::assembly::{AssemblyDescription, DEFAULT_ASSEMBLY_BUDGET, ParameterPath};
use crate::gesture::{GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue};
use fs_blake3::{ContentHash, hash_bytes};
use fs_qty::{Dims, QtyAny};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

pub mod smf;
/// Frozen wire version; interpretation inputs are content addressed separately.
pub const PERFORMANCE_SCHEMA_VERSION: u32 = 1;

#[derive(Clone, Copy, Debug)]
pub struct SmfLimits {
    pub max_input_bytes: usize,
    pub max_tracks: usize,
    /// Includes metadata and end-of-track events.
    pub max_events: usize,
    pub max_control_ticks: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum PerformanceError {
    Invalid { what: &'static str },
    Truncated { offset: usize },
    Limit { resource: &'static str, cap: u64 },
    UnsupportedFormat { format: u16 },
    SmpteDivision { division: u16 },
    OverlappingNote { channel: u8, note: u8 },
    UnmatchedNoteOff { channel: u8, note: u8 },
    UnterminatedNote { channel: u8, note: u8 },
    ZeroLengthNote { channel: u8, note: u8 },
    TempoOutsideConductor { track: u16 },
    UnsupportedMeta { kind: u8 },
    UnsupportedMessage { track: u16, status: u8 },
    MissingTempo,
    UnknownNote { channel: u8, note: u8 },
    UnknownController { channel: u8, controller: u8 },
    UnknownConfiguration { name: String },
    MissingVelocityBin { velocity: u8 },
    MissingDynamicRow { marking: String },
    MissingBinding { path: ParameterPath },
    PhysicalStateCommand { path: ParameterPath },
    UnresolvedControl { index: usize },
    ClockRepresentation { tick: u64 },
    HashMismatch { input: &'static str },
    ConflictingCommands { track: String, tick: u64 },
}
impl fmt::Display for PerformanceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "performance refusal: {self:?}")
    }
}
impl std::error::Error for PerformanceError {}
pub(crate) fn invalid(what: &'static str) -> PerformanceError {
    PerformanceError::Invalid { what }
}
pub(crate) fn cap(n: usize, max: usize, resource: &'static str) -> Result<(), PerformanceError> {
    if n > max {
        Err(PerformanceError::Limit {
            resource,
            cap: max as u64,
        })
    } else {
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum DeclaredObjective {
    /// Declaration only: no MIDI-number to frequency conversion occurs here.
    EqualTemperamentPitch {
        a4_hz: f64,
        ratio_source: Option<String>,
    },
    SlotCentre,
    MinimumThresholdPressure,
}
#[derive(Clone, Debug, PartialEq)]
pub struct PerNoteControl {
    pub note_index: usize,
    pub target: ParameterPath,
    pub objective: DeclaredObjective,
    pub domain: (QtyAny, QtyAny),
    pub resolved: Option<QtyAny>,
    pub receipt: Option<ContentHash>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct NoteEvent {
    pub on_tick: u64,
    pub off_tick: u64,
    pub configuration: String,
    pub dynamic: String,
    pub velocity: u8,
}
#[derive(Clone, Debug)]
pub struct NoteChartEntry {
    pub channel: u8,
    pub note: u8,
    pub configuration: String,
}
#[derive(Clone, Debug)]
pub struct NoteChart {
    pub source: ContentHash,
    pub entries: Vec<NoteChartEntry>,
}
#[derive(Clone, Debug)]
pub struct VelocityBin {
    pub first: u8,
    pub last: u8,
    pub marking: String,
}
#[derive(Clone, Debug)]
pub struct VelocityBins {
    pub source: ContentHash,
    pub bins: Vec<VelocityBin>,
}
#[derive(Clone, Debug)]
pub enum ControlChoice {
    /// A measured value or an owner-authored performance; source is mandatory.
    Fixed {
        value: QtyAny,
        source: ContentHash,
    },
    PerNote {
        objective: DeclaredObjective,
    },
}
#[derive(Clone, Debug)]
pub struct DynamicControl {
    pub target: ParameterPath,
    pub domain: (QtyAny, QtyAny),
    pub choice: ControlChoice,
}
#[derive(Clone, Debug)]
pub struct DynamicsRow {
    pub marking: String,
    pub controls: Vec<DynamicControl>,
}
#[derive(Clone, Debug)]
pub struct DynamicsTable {
    pub source: ContentHash,
    pub rows: Vec<DynamicsRow>,
}

/// Only declared root-level player commands may be prescribed. Interaction
/// states on parts (airway/contact/tissue/holding) belong to coupled physics.
#[derive(Clone, Debug)]
pub struct CommandBinding {
    pub path: ParameterPath,
    pub track_id: String,
    pub target: GestureTarget,
    /// Also declares the scalar gesture unit and the note-release command.
    pub initial: GestureValue,
}
#[derive(Clone, Debug)]
pub struct ControllerRule {
    pub channel: u8,
    pub controller: u8,
    pub target: ParameterPath,
    /// Linear interpolation over MIDI values 0..127, explicitly declared.
    pub range: (QtyAny, QtyAny),
}
#[derive(Clone, Debug)]
pub struct ImportPolicy {
    pub chart: NoteChart,
    pub velocities: VelocityBins,
    pub bindings: Vec<CommandBinding>,
    pub controllers: Vec<ControllerRule>,
    pub control_rate_hz: u32,
    /// None requires an SMF Set Tempo at tick zero; no implicit 120 BPM.
    pub initial_tempo_us: Option<u32>,
    /// A prior declared before grading, retained as a content-addressed input.
    pub observation_prior: ContentHash,
}

#[derive(Clone, Debug)]
pub struct Performance {
    description: ContentHash,
    policy: ContentHash,
    dynamics: ContentHash,
    midi: Vec<u8>,
    notes: Vec<NoteEvent>,
    controllers: Vec<(u64, smf::MidiController)>,
    per_note: Vec<PerNoteControl>,
    control_rate_hz: u32,
    observation_prior: ContentHash,
}
#[derive(Clone, Debug)]
pub struct PerformanceProvenance {
    pub performance: ContentHash,
    pub description: ContentHash,
    pub midi: ContentHash,
    pub chart: ContentHash,
    pub dynamics: ContentHash,
    pub observation_prior: ContentHash,
    pub per_note: Vec<PerNoteControl>,
}
#[derive(Clone, Debug)]
pub struct LoweredPerformance {
    pub schedule: GestureSchedule,
    pub provenance: PerformanceProvenance,
}

// Length-prefix every variable field; float bit patterns are preserved.
fn text(b: &mut Vec<u8>, s: &str) {
    b.extend_from_slice(&(s.len() as u64).to_be_bytes());
    b.extend_from_slice(s.as_bytes());
}
fn qty(b: &mut Vec<u8>, q: QtyAny) {
    b.extend(q.dims.0.iter().map(|v| *v as u8));
    b.extend_from_slice(&q.value.to_bits().to_be_bytes());
}
fn path(b: &mut Vec<u8>, p: &ParameterPath) {
    text(b, &p.entity.token());
    text(b, &p.parameter.name);
    b.extend(p.parameter.dims.0.iter().map(|v| *v as u8));
}
fn objective(b: &mut Vec<u8>, o: &DeclaredObjective) {
    match o {
        DeclaredObjective::EqualTemperamentPitch {
            a4_hz,
            ratio_source,
        } => {
            b.push(0);
            b.extend_from_slice(&a4_hz.to_bits().to_be_bytes());
            b.push(u8::from(ratio_source.is_some()));
            if let Some(s) = ratio_source {
                text(b, s);
            }
        }
        DeclaredObjective::SlotCentre => b.push(1),
        DeclaredObjective::MinimumThresholdPressure => b.push(2),
    }
}
impl NoteChart {
    pub fn content_hash(&self) -> ContentHash {
        let mut b = b"fs-note-chart-v1".to_vec();
        b.extend_from_slice(self.source.as_bytes());
        for e in &self.entries {
            b.push(e.channel);
            b.push(e.note);
            text(&mut b, &e.configuration);
        }
        hash_bytes(&b)
    }
}
impl DynamicsTable {
    pub fn content_hash(&self) -> ContentHash {
        let mut b = b"fs-dynamics-v1".to_vec();
        b.extend_from_slice(self.source.as_bytes());
        for row in &self.rows {
            text(&mut b, &row.marking);
            b.extend_from_slice(&(row.controls.len() as u64).to_be_bytes());
            for c in &row.controls {
                path(&mut b, &c.target);
                qty(&mut b, c.domain.0);
                qty(&mut b, c.domain.1);
                match &c.choice {
                    ControlChoice::Fixed { value, source } => {
                        b.push(0);
                        qty(&mut b, *value);
                        b.extend_from_slice(source.as_bytes());
                    }
                    ControlChoice::PerNote { objective: o } => {
                        b.push(1);
                        objective(&mut b, o);
                    }
                }
            }
        }
        hash_bytes(&b)
    }
}
impl ImportPolicy {
    pub fn content_hash(&self) -> Result<ContentHash, PerformanceError> {
        let mut b = b"fs-performance-policy-v1".to_vec();
        b.extend_from_slice(self.chart.content_hash().as_bytes());
        b.extend_from_slice(self.velocities.source.as_bytes());
        b.extend_from_slice(&(self.velocities.bins.len() as u64).to_be_bytes());
        for v in &self.velocities.bins {
            b.extend_from_slice(&[v.first, v.last]);
            text(&mut b, &v.marking);
        }
        b.extend_from_slice(&self.control_rate_hz.to_be_bytes());
        b.push(u8::from(self.initial_tempo_us.is_some()));
        b.extend_from_slice(&self.initial_tempo_us.unwrap_or(0).to_be_bytes());
        b.extend_from_slice(self.observation_prior.as_bytes());
        b.extend_from_slice(&(self.bindings.len() as u64).to_be_bytes());
        for bind in &self.bindings {
            path(&mut b, &bind.path);
            let s = GestureSchedule::try_new(
                self.control_rate_hz,
                vec![GestureTrack {
                    id: bind.track_id.clone(),
                    target: bind.target.clone(),
                    initial: bind.initial.clone(),
                    events: vec![],
                }],
            )
            .map_err(|_| invalid("command binding gesture"))?;
            let wire = s.to_canonical_bytes();
            b.extend_from_slice(&(wire.len() as u64).to_be_bytes());
            b.extend(wire);
        }
        for c in &self.controllers {
            b.extend_from_slice(&[c.channel, c.controller]);
            path(&mut b, &c.target);
            qty(&mut b, c.range.0);
            qty(&mut b, c.range.1);
        }
        Ok(hash_bytes(&b))
    }
}

fn check_domain(p: &ParameterPath, d: (QtyAny, QtyAny)) -> Result<(), PerformanceError> {
    if d.0.dims != p.parameter.dims
        || d.1.dims != d.0.dims
        || !d.0.value.is_finite()
        || !d.1.value.is_finite()
        || d.0.value > d.1.value
    {
        return Err(invalid("control domain units or bounds"));
    }
    Ok(())
}
fn in_domain(v: QtyAny, d: (QtyAny, QtyAny)) -> bool {
    v.dims == d.0.dims && v.value.is_finite() && v.value >= d.0.value && v.value <= d.1.value
}
fn value(template: &GestureValue, q: QtyAny) -> Result<GestureValue, PerformanceError> {
    let (dims, result) = match template {
        GestureValue::PressurePa(_) => (
            Dims([-1, 1, -2, 0, 0, 0]),
            GestureValue::PressurePa(q.value),
        ),
        GestureValue::LengthM(_) => (Dims([1, 0, 0, 0, 0, 0]), GestureValue::LengthM(q.value)),
        GestureValue::Fraction(_) => (Dims::NONE, GestureValue::Fraction(q.value)),
        GestureValue::Multiplier(_) => (Dims::NONE, GestureValue::Multiplier(q.value)),
        GestureValue::TensionN(_) => (Dims([1, 1, -2, 0, 0, 0]), GestureValue::TensionN(q.value)),
        GestureValue::VelocityMPerS(_) => (
            Dims([1, 0, -1, 0, 0, 0]),
            GestureValue::VelocityMPerS(q.value),
        ),
        GestureValue::AngleRad(_) => (Dims::NONE, GestureValue::AngleRad(q.value)),
        _ => {
            return Err(invalid(
                "structured gesture needs a declared articulation adapter",
            ));
        }
    };
    if q.dims != dims || !q.value.is_finite() {
        return Err(invalid("command units or nonfinite value"));
    }
    Ok(result)
}
fn binding<'a>(
    policy: &'a ImportPolicy,
    p: &ParameterPath,
) -> Result<&'a CommandBinding, PerformanceError> {
    policy
        .bindings
        .iter()
        .find(|b| &b.path == p)
        .ok_or_else(|| PerformanceError::MissingBinding { path: p.clone() })
}
fn validate(
    desc: &AssemblyDescription,
    policy: &ImportPolicy,
    table: &DynamicsTable,
) -> Result<(), PerformanceError> {
    if policy.control_rate_hz == 0 || policy.initial_tempo_us == Some(0) {
        return Err(invalid("control rate or initial tempo"));
    }
    let mut paths = BTreeSet::new();
    let mut tracks = BTreeSet::new();
    let mut targets = Vec::new();
    for b in &policy.bindings {
        if b.path.entity != desc.input().root
            || matches!(
                b.target,
                GestureTarget::PreStress
                    | GestureTarget::RestAperture
                    | GestureTarget::StringTension { .. }
            )
        {
            return Err(PerformanceError::PhysicalStateCommand {
                path: b.path.clone(),
            });
        }
        if !paths.insert(b.path.clone())
            || !tracks.insert(b.track_id.clone())
            || targets.contains(&&b.target)
        {
            return Err(invalid("duplicate command path, track or gesture target"));
        }
        targets.push(&b.target);
        // Check declared units even when a track has no score events.
        value(&b.initial, QtyAny::new(0.0, b.path.parameter.dims))?;
    }
    let mut notes = BTreeSet::new();
    for e in &policy.chart.entries {
        if e.channel > 15 || e.note > 127 || !notes.insert((e.channel, e.note)) {
            return Err(invalid("duplicate or invalid note chart key"));
        }
        if !desc
            .input()
            .configurations
            .iter()
            .any(|c| c.name == e.configuration)
        {
            return Err(PerformanceError::UnknownConfiguration {
                name: e.configuration.clone(),
            });
        }
    }
    let mut velocities = BTreeSet::new();
    for v in &policy.velocities.bins {
        if v.first == 0 || v.last > 127 || v.first > v.last || v.marking.is_empty() {
            return Err(invalid("velocity bin"));
        }
        for n in v.first..=v.last {
            if !velocities.insert(n) {
                return Err(invalid("overlapping velocity bins"));
            }
        }
    }
    let mut markings = BTreeSet::new();
    for r in &table.rows {
        if r.marking.is_empty() || !markings.insert(&r.marking) || r.controls.is_empty() {
            return Err(invalid("duplicate or empty dynamics row"));
        }
        let mut row_paths = BTreeSet::new();
        for c in &r.controls {
            if !row_paths.insert(c.target.clone()) {
                return Err(invalid("duplicate dynamics control"));
            }
            let b = binding(policy, &c.target)?;
            check_domain(&c.target, c.domain)?;
            value(&b.initial, c.domain.0)?;
            value(&b.initial, c.domain.1)?;
            match &c.choice {
                ControlChoice::Fixed { value: v, .. } if !in_domain(*v, c.domain) => {
                    return Err(invalid("fixed control outside domain"));
                }
                ControlChoice::PerNote {
                    objective:
                        DeclaredObjective::EqualTemperamentPitch {
                            a4_hz,
                            ratio_source,
                        },
                } if !a4_hz.is_finite()
                    || *a4_hz <= 0.0
                    || ratio_source.as_ref().is_some_and(|s| s.is_empty()) =>
                {
                    return Err(invalid("pitch objective declaration"));
                }
                _ => {}
            }
        }
    }
    let mut cc = BTreeSet::new();
    for c in &policy.controllers {
        if c.channel > 15 || c.controller > 127 || !cc.insert((c.channel, c.controller)) {
            return Err(invalid("duplicate or invalid controller rule"));
        }
        let b = binding(policy, &c.target)?;
        check_domain(&c.target, c.range)?;
        value(&b.initial, c.range.0)?;
        value(&b.initial, c.range.1)?;
        // Independent CC ownership avoids silently overwriting a score command.
        if table
            .rows
            .iter()
            .any(|r| r.controls.iter().any(|v| v.target == c.target))
            || desc
                .input()
                .configurations
                .iter()
                .any(|s| s.assignments.iter().any(|(p, _)| p == &c.target))
        {
            return Err(invalid("controller and note command share ownership"));
        }
    }
    policy.content_hash()?;
    Ok(())
}

/// Import one performer. Overlap, unmapped events and unowned controls refuse.
pub fn import_smf(
    bytes: &[u8],
    desc: &AssemblyDescription,
    policy: &ImportPolicy,
    table: &DynamicsTable,
    limits: SmfLimits,
) -> Result<Performance, PerformanceError> {
    validate(desc, policy, table)?;
    let midi = smf::decode(bytes, limits)?;
    let clock = |t| {
        midi.control_tick(
            t,
            policy.control_rate_hz,
            policy.initial_tempo_us,
            limits.max_control_ticks,
        )
    };
    clock(midi.end_tick)?;
    let mut notes = Vec::new();
    let mut per_note = Vec::new();
    for n in &midi.notes {
        let chart = policy
            .chart
            .entries
            .iter()
            .find(|e| e.channel == n.channel && e.note == n.note)
            .ok_or(PerformanceError::UnknownNote {
                channel: n.channel,
                note: n.note,
            })?;
        let dynamic = &policy
            .velocities
            .bins
            .iter()
            .find(|v| n.velocity >= v.first && n.velocity <= v.last)
            .ok_or(PerformanceError::MissingVelocityBin {
                velocity: n.velocity,
            })?
            .marking;
        let row = table
            .rows
            .iter()
            .find(|r| &r.marking == dynamic)
            .ok_or_else(|| PerformanceError::MissingDynamicRow {
                marking: dynamic.clone(),
            })?;
        let on_tick = clock(n.on)?;
        let off_tick = clock(n.off)?;
        if on_tick == off_tick {
            return Err(PerformanceError::ZeroLengthNote {
                channel: n.channel,
                note: n.note,
            });
        }
        for c in &row.controls {
            if let ControlChoice::PerNote { objective } = &c.choice {
                per_note.push(PerNoteControl {
                    note_index: notes.len(),
                    target: c.target.clone(),
                    objective: objective.clone(),
                    domain: c.domain,
                    resolved: None,
                    receipt: None,
                });
            }
        }
        notes.push(NoteEvent {
            on_tick,
            off_tick,
            configuration: chart.configuration.clone(),
            dynamic: dynamic.clone(),
            velocity: n.velocity,
        });
    }
    let mut controllers = Vec::new();
    for c in &midi.controllers {
        if !policy
            .controllers
            .iter()
            .any(|r| r.channel == c.channel && r.controller == c.controller)
        {
            return Err(PerformanceError::UnknownController {
                channel: c.channel,
                controller: c.controller,
            });
        }
        controllers.push((clock(c.tick)?, c.clone()));
    }
    Ok(Performance {
        description: desc
            .content_hash(DEFAULT_ASSEMBLY_BUDGET)
            .map_err(|_| invalid("assembly hash"))?,
        policy: policy.content_hash()?,
        dynamics: table.content_hash(),
        midi: bytes.to_vec(),
        notes,
        controllers,
        per_note,
        control_rate_hz: policy.control_rate_hz,
        observation_prior: policy.observation_prior,
    })
}

impl Performance {
    pub fn notes(&self) -> &[NoteEvent] {
        &self.notes
    }
    pub fn per_note(&self) -> &[PerNoteControl] {
        &self.per_note
    }
    /// Only a bounded, unit-correct result with an evaluation receipt is stored.
    pub fn resolve(
        &mut self,
        index: usize,
        result: QtyAny,
        receipt: ContentHash,
    ) -> Result<(), PerformanceError> {
        let c = self
            .per_note
            .get_mut(index)
            .ok_or_else(|| invalid("per-note index"))?;
        if !in_domain(result, c.domain) {
            return Err(invalid("resolution outside declared domain"));
        }
        c.resolved = Some(result);
        c.receipt = Some(receipt);
        Ok(())
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut b = b"FSPERF\0".to_vec();
        b.extend_from_slice(&PERFORMANCE_SCHEMA_VERSION.to_be_bytes());
        for h in [self.description, self.policy, self.dynamics] {
            b.extend_from_slice(h.as_bytes());
        }
        b.extend_from_slice(&(self.midi.len() as u64).to_be_bytes());
        b.extend_from_slice(&self.midi);
        b.extend_from_slice(&(self.per_note.len() as u64).to_be_bytes());
        for c in &self.per_note {
            b.push(u8::from(c.resolved.is_some()));
            if let (Some(v), Some(h)) = (c.resolved, c.receipt) {
                b.extend_from_slice(&v.value.to_bits().to_be_bytes());
                b.extend_from_slice(h.as_bytes());
            }
        }
        b
    }
    pub fn content_hash(&self) -> ContentHash {
        hash_bytes(&self.to_bytes())
    }
    /// References are mandatory and verified before any interpretation is used.
    pub fn from_bytes(
        bytes: &[u8],
        desc: &AssemblyDescription,
        policy: &ImportPolicy,
        table: &DynamicsTable,
        limits: SmfLimits,
    ) -> Result<Self, PerformanceError> {
        let mut r = smf::Reader { bytes, at: 0 };
        if r.take(7)? != b"FSPERF\0" || r.u32()? != PERFORMANCE_SCHEMA_VERSION {
            return Err(invalid("performance schema"));
        }
        for (name, hash) in [
            (
                "description",
                desc.content_hash(DEFAULT_ASSEMBLY_BUDGET)
                    .map_err(|_| invalid("assembly hash"))?,
            ),
            ("policy", policy.content_hash()?),
            ("dynamics", table.content_hash()),
        ] {
            if r.take(32)? != hash.as_bytes() {
                return Err(PerformanceError::HashMismatch { input: name });
            }
        }
        let len = usize::try_from(r.u64()?).map_err(|_| invalid("MIDI extent"))?;
        cap(len, limits.max_input_bytes, "MIDI input bytes")?;
        let mut p = import_smf(r.take(len)?, desc, policy, table, limits)?;
        if r.u64()? != p.per_note.len() as u64 {
            return Err(invalid("resolution count"));
        }
        for i in 0..p.per_note.len() {
            match r.byte()? {
                0 => {}
                1 => {
                    let v = f64::from_bits(r.u64()?);
                    let h = ContentHash::from_slice(r.take(32)?)
                        .ok_or_else(|| invalid("receipt hash"))?;
                    p.resolve(i, QtyAny::new(v, p.per_note[i].domain.0.dims), h)?;
                }
                _ => return Err(invalid("resolution flag")),
            }
        }
        if r.at != bytes.len() {
            return Err(invalid("trailing performance bytes"));
        }
        Ok(p)
    }
}

fn seconds(tick: u64, rate: u32) -> Result<f64, PerformanceError> {
    let t = tick as f64 / f64::from(rate);
    if tick >= (1u64 << 53) || t * f64::from(rate) != tick as f64 {
        return Err(PerformanceError::ClockRepresentation { tick });
    }
    Ok(t)
}
/// Lower only after all free controls have model-only resolution receipts.
pub fn lower(
    perf: &Performance,
    desc: &AssemblyDescription,
    policy: &ImportPolicy,
    table: &DynamicsTable,
) -> Result<LoweredPerformance, PerformanceError> {
    validate(desc, policy, table)?;
    for (name, a, b) in [
        (
            "description",
            perf.description,
            desc.content_hash(DEFAULT_ASSEMBLY_BUDGET)
                .map_err(|_| invalid("assembly hash"))?,
        ),
        ("policy", perf.policy, policy.content_hash()?),
        ("dynamics", perf.dynamics, table.content_hash()),
    ] {
        if a != b {
            return Err(PerformanceError::HashMismatch { input: name });
        }
    }
    for (i, c) in perf.per_note.iter().enumerate() {
        if c.resolved.is_none() || c.receipt.is_none() {
            return Err(PerformanceError::UnresolvedControl { index: i });
        }
    }
    let mut events: BTreeMap<ParameterPath, BTreeMap<u64, GestureValue>> = policy
        .bindings
        .iter()
        .map(|b| (b.path.clone(), BTreeMap::new()))
        .collect();
    // Process notes in score order: a same-tick release precedes the next attack.
    for (index, n) in perf.notes.iter().enumerate() {
        let config = desc
            .input()
            .configurations
            .iter()
            .find(|c| c.name == n.configuration)
            .ok_or_else(|| PerformanceError::UnknownConfiguration {
                name: n.configuration.clone(),
            })?;
        let row = table
            .rows
            .iter()
            .find(|r| r.marking == n.dynamic)
            .ok_or_else(|| PerformanceError::MissingDynamicRow {
                marking: n.dynamic.clone(),
            })?;
        let mut commands = BTreeMap::new();
        for (p, q) in &config.assignments {
            commands.insert(p.clone(), *q);
        }
        for c in &row.controls {
            let q = match &c.choice {
                ControlChoice::Fixed { value, .. } => *value,
                ControlChoice::PerNote { .. } => perf
                    .per_note
                    .iter()
                    .find(|r| r.note_index == index && r.target == c.target)
                    .and_then(|r| r.resolved)
                    .ok_or_else(|| invalid("missing resolved control"))?,
            };
            if commands.insert(c.target.clone(), q).is_some() {
                return Err(invalid(
                    "configuration and dynamics share command ownership",
                ));
            }
        }
        for (p, q) in commands {
            let b = binding(policy, &p)?;
            let v = value(&b.initial, q)?;
            let track = events
                .get_mut(&p)
                .ok_or_else(|| PerformanceError::MissingBinding { path: p.clone() })?;
            track.insert(n.on_tick, v);
            track.insert(n.off_tick, b.initial.clone());
        }
    }
    for (tick, c) in &perf.controllers {
        let rule = policy
            .controllers
            .iter()
            .find(|r| r.channel == c.channel && r.controller == c.controller)
            .ok_or(PerformanceError::UnknownController {
                channel: c.channel,
                controller: c.controller,
            })?;
        let b = binding(policy, &rule.target)?;
        let q = QtyAny::new(
            rule.range.0.value
                + (rule.range.1.value - rule.range.0.value) * (f64::from(c.value) / 127.0),
            rule.range.0.dims,
        );
        let v = value(&b.initial, q)?;
        let track = events
            .get_mut(&rule.target)
            .ok_or_else(|| invalid("controller binding"))?;
        if track.get(tick).is_some_and(|old| old != &v) {
            return Err(PerformanceError::ConflictingCommands {
                track: b.track_id.clone(),
                tick: *tick,
            });
        }
        track.insert(*tick, v);
    }
    let mut tracks = Vec::new();
    for b in &policy.bindings {
        let es = events
            .remove(&b.path)
            .ok_or_else(|| invalid("command track"))?;
        let mut converted = Vec::new();
        for (tick, v) in es {
            converted.push(GestureEvent {
                time_s: seconds(tick, perf.control_rate_hz)?,
                transition_s: 0.0,
                value: v,
            });
        }
        tracks.push(GestureTrack {
            id: b.track_id.clone(),
            target: b.target.clone(),
            initial: b.initial.clone(),
            events: converted,
        });
    }
    let schedule = GestureSchedule::try_new(perf.control_rate_hz, tracks)
        .map_err(|_| invalid("lowered gesture schedule"))?;
    Ok(LoweredPerformance {
        schedule,
        provenance: PerformanceProvenance {
            performance: perf.content_hash(),
            description: perf.description,
            midi: hash_bytes(&perf.midi),
            chart: policy.chart.content_hash(),
            dynamics: perf.dynamics,
            observation_prior: perf.observation_prior,
            per_note: perf.per_note.clone(),
        },
    })
}
