//! Canonical scenario IR: a deterministic s-expression encoding with
//! LOSSLESS round-trip (floats print in shortest-round-trip form; dims
//! travel as explicit SI exponent vectors, so no unit-string parsing is
//! involved). Versioned typed boundary payloads are embedded as lowercase hex
//! of their canonical bounded payload envelope; historical scenario v1 refuses
//! that form rather than inventing a five-base crosswalk. `write_ir` output is
//! byte-stable — the ledger stores it as a content-addressed artifact. Parsing
//! inverts it exactly when the caller's explicit resource budget admits the
//! complete artifact; the convenience [`parse_ir`] retains a conservative
//! 16 MiB total-input ceiling.

use crate::bc::{BcKind, BcValue, BoundaryCondition, Compat, Physics};
use crate::ensemble::{SpectrumModel, StochasticEnsemble};
use crate::frame::{Frame, FrameId, FrameMotion, FrameTree};
use crate::payload::{
    MAX_PAYLOAD_WIRE_BYTES, PAYLOAD_WIRE_VERSION, Payload, PayloadDecodeLimits,
    canonical_payload_bytes, decode_payload_with_limits,
};
use crate::scenario::{
    Combination, ContactLaw, ContactModel, Environment, LoadCase, Scenario, Violation,
};
use crate::signal::{ChebProfile, Interp, TimeSignal};
use crate::{IrSourceSpan, ScenarioError};
use fs_blake3::{ContentHash, hash_bytes};
use fs_cheb::Cheb1;
use fs_ga::{Quat, Vec3};
use fs_qty::{Dims, QtyAny};
use std::fmt::Write as _;

/// Canonical scenario-IR version written by this build.
pub const SCENARIO_IR_VERSION: u32 = 2;
/// Historical unversioned scenario-IR semantics: five SI base exponents.
pub const LEGACY_SCENARIO_IR_VERSION: u32 = 1;

/// The only admitted semantic rule for historical five-base dimensions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FiveToSixRule {
    /// Preserve the five exponents and append an exact zero mole exponent.
    AppendMoleZero,
}

/// Immutable evidence binding exact historical bytes to exact canonical v2.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DimensionCrosswalkReceipt {
    /// Source scenario-IR version.
    source_version: u32,
    /// Canonical scenario-IR version.
    target_version: u32,
    /// BLAKE3 hash of the exact supplied historical bytes.
    old_hash: ContentHash,
    /// BLAKE3 hash of the exact canonical v2 bytes.
    new_hash: ContentHash,
    /// Number of source dimension exponents.
    source_width: u8,
    /// Number of canonical dimension exponents.
    target_width: u8,
    /// Semantic migration rule applied.
    rule: FiveToSixRule,
}

impl DimensionCrosswalkReceipt {
    /// Historical source version.
    #[must_use]
    pub const fn source_version(&self) -> u32 {
        self.source_version
    }

    /// Canonical target version.
    #[must_use]
    pub const fn target_version(&self) -> u32 {
        self.target_version
    }

    /// BLAKE3 hash of the exact source bytes.
    #[must_use]
    pub const fn old_hash(&self) -> ContentHash {
        self.old_hash
    }

    /// BLAKE3 hash of the exact canonical target bytes.
    #[must_use]
    pub const fn new_hash(&self) -> ContentHash {
        self.new_hash
    }

    /// Source dimension-vector width.
    #[must_use]
    pub const fn source_width(&self) -> u8 {
        self.source_width
    }

    /// Target dimension-vector width.
    #[must_use]
    pub const fn target_width(&self) -> u8 {
        self.target_width
    }

    /// Semantic rule used by the crosswalk.
    #[must_use]
    pub const fn rule(&self) -> FiveToSixRule {
        self.rule
    }

    /// Verify the receipt against exact preserved source and canonical bytes.
    #[must_use]
    pub fn verifies(&self, old_bytes: &[u8], new_bytes: &[u8]) -> bool {
        self.source_version == LEGACY_SCENARIO_IR_VERSION
            && self.target_version == SCENARIO_IR_VERSION
            && self.source_width == 5
            && self.target_width == 6
            && self.rule == FiveToSixRule::AppendMoleZero
            && hash_bytes(old_bytes) == self.old_hash
            && hash_bytes(new_bytes) == self.new_hash
    }
}

/// The semantic rule witnessed by a current-version source receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceCanonicalizationRule {
    /// Decode admitted v2 source text and re-emit it through [`write_ir`].
    CanonicalV2Reemission,
}

impl SourceCanonicalizationRule {
    const fn for_source_version(version: u32) -> Option<Self> {
        match version {
            2 => Some(Self::CanonicalV2Reemission),
            _ => None,
        }
    }

    const fn source_version(self) -> u32 {
        match self {
            Self::CanonicalV2Reemission => 2,
        }
    }
}

/// Immutable evidence binding accepted noncanonical v2 bytes to canonical v2.
///
/// Canonical writer output does not need this receipt. It is present only when
/// an accepted current-version source uses an alternate spelling or layout,
/// so callers never lose the exact authority identity that was supplied.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceCanonicalizationReceipt {
    /// Version declared by the supplied source.
    source_version: u32,
    /// Version emitted by the canonical writer.
    canonical_version: u32,
    /// BLAKE3 hash of the exact supplied source bytes.
    source_hash: ContentHash,
    /// BLAKE3 hash of the exact canonical re-emission.
    canonical_hash: ContentHash,
    /// Canonicalization rule applied.
    rule: SourceCanonicalizationRule,
}

impl SourceCanonicalizationReceipt {
    /// Version declared by the supplied source.
    #[must_use]
    pub const fn source_version(&self) -> u32 {
        self.source_version
    }

    /// Version emitted by the canonical writer.
    #[must_use]
    pub const fn canonical_version(&self) -> u32 {
        self.canonical_version
    }

    /// BLAKE3 hash of the exact supplied source bytes.
    #[must_use]
    pub const fn source_hash(&self) -> ContentHash {
        self.source_hash
    }

    /// BLAKE3 hash of the exact canonical re-emission.
    #[must_use]
    pub const fn canonical_hash(&self) -> ContentHash {
        self.canonical_hash
    }

    /// Semantic canonicalization rule used by this receipt.
    #[must_use]
    pub const fn rule(&self) -> SourceCanonicalizationRule {
        self.rule
    }

    /// Verify the parser-issued receipt against exact preserved byte strings.
    ///
    /// This authenticates the immutable source/canonical hash pair created by
    /// the parser. It does not independently parse arbitrary target bytes to
    /// establish canonicality.
    #[must_use]
    pub fn verifies(&self, source_bytes: &[u8], canonical_bytes: &[u8]) -> bool {
        self.source_version == self.rule.source_version()
            && self.canonical_version == self.rule.source_version()
            && source_bytes != canonical_bytes
            && self.source_hash != self.canonical_hash
            && hash_bytes(source_bytes) == self.source_hash
            && hash_bytes(canonical_bytes) == self.canonical_hash
    }
}

/// A decoded scenario together with mutually exclusive wire-source evidence.
///
/// Exact canonical v2 has neither receipt, accepted noncanonical v2 has only a
/// [`SourceCanonicalizationReceipt`], and legacy v1 has only a
/// [`DimensionCrosswalkReceipt`]. `PartialEq` includes that source evidence;
/// compare [`DecodedScenario::scenario`] values when only semantic equality is
/// intended.
#[derive(Debug, Clone, PartialEq)]
pub struct DecodedScenario {
    /// Decoded scenario value.
    scenario: Scenario,
    /// Version found on the wire; historical unversioned forms are v1.
    source_version: u32,
    /// Present exactly when a legacy five-base form was crossed into six-base memory.
    dimension_crosswalk: Option<DimensionCrosswalkReceipt>,
    /// Present exactly when admitted v2 source bytes were not canonical writer output.
    source_canonicalization: Option<SourceCanonicalizationReceipt>,
}

impl DecodedScenario {
    /// Decoded six-base scenario.
    #[must_use]
    pub const fn scenario(&self) -> &Scenario {
        &self.scenario
    }

    /// Version found on the supplied wire bytes.
    #[must_use]
    pub const fn source_version(&self) -> u32 {
        self.source_version
    }

    /// Mandatory receipt for v1 input; absent for canonical v2.
    #[must_use]
    pub const fn migration(&self) -> Option<&DimensionCrosswalkReceipt> {
        self.dimension_crosswalk.as_ref()
    }

    /// Receipt binding accepted noncanonical v2 bytes to canonical v2 output.
    #[must_use]
    pub const fn canonicalization(&self) -> Option<&SourceCanonicalizationReceipt> {
        self.source_canonicalization.as_ref()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DimensionWire {
    LegacyFive,
    CanonicalSix,
}

// ---------------------------------------------------------------- writer

fn canonical_float(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn w_qty(out: &mut String, q: &QtyAny) {
    let _ = write!(out, "(qty {}", canonical_float(q.value));
    for d in q.dims.0 {
        let _ = write!(out, " {d}");
    }
    out.push(')');
}

fn w_dims(out: &mut String, d: Dims) {
    out.push_str("(dims");
    for e in d.0 {
        let _ = write!(out, " {e}");
    }
    out.push(')');
}

fn w_str(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        if c == '"' || c == '\\' {
            out.push('\\');
        }
        out.push(c);
    }
    out.push('"');
}

fn w_payload_hex(out: &mut String, payload: &Payload) {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let bytes = canonical_payload_bytes(payload);
    let _ = write!(out, "(typed :version {PAYLOAD_WIRE_VERSION} \"");
    for byte in bytes {
        out.push(char::from(HEX[usize::from(byte >> 4)]));
        out.push(char::from(HEX[usize::from(byte & 0x0f)]));
    }
    out.push_str("\")");
}

fn w_vec3(out: &mut String, v: Vec3) {
    let _ = write!(
        out,
        "(vec {} {} {})",
        canonical_float(v.x),
        canonical_float(v.y),
        canonical_float(v.z)
    );
}

fn w_floats(out: &mut String, head: &str, vs: &[f64]) {
    let _ = write!(out, "({head}");
    for v in vs {
        let _ = write!(out, " {}", canonical_float(*v));
    }
    out.push(')');
}

fn w_signal(out: &mut String, s: &TimeSignal) {
    match s {
        TimeSignal::Constant(q) => {
            out.push_str("(constant ");
            w_qty(out, q);
            out.push(')');
        }
        TimeSignal::Ramp {
            t_start,
            t_end,
            from,
            to,
        } => {
            let _ = write!(
                out,
                "(ramp {} {} ",
                canonical_float(*t_start),
                canonical_float(*t_end)
            );
            w_qty(out, from);
            out.push(' ');
            w_qty(out, to);
            out.push(')');
        }
        TimeSignal::Table {
            times,
            values,
            dims,
            interp,
        } => {
            let tag = match interp {
                Interp::Linear => "linear",
                Interp::Hold => "hold",
            };
            let _ = write!(out, "(table {tag} ");
            w_dims(out, *dims);
            out.push(' ');
            w_floats(out, "times", times);
            out.push(' ');
            w_floats(out, "values", values);
            out.push(')');
        }
        TimeSignal::Chebfun(p) => {
            w_profile(out, "chebfun", p);
        }
    }
}

fn w_profile(out: &mut String, head: &str, p: &ChebProfile) {
    let (a, b) = p.cheb.domain();
    let _ = write!(out, "({head} ");
    w_dims(out, p.dims);
    let _ = write!(out, " {} {} ", canonical_float(a), canonical_float(b));
    w_floats(out, "coeffs", p.cheb.coeffs());
    out.push(')');
}

fn physics_tag(p: Physics) -> &'static str {
    match p {
        Physics::IncompressibleFlow => "incompressible-flow",
        Physics::Thermal => "thermal",
        Physics::Elasticity => "elasticity",
        Physics::Magnetics => "magnetics",
        Physics::Electrics => "electrics",
        Physics::GasExchange => "gas-exchange",
    }
}

fn kind_tag(k: BcKind) -> &'static str {
    match k {
        BcKind::Dirichlet => "dirichlet",
        BcKind::Neumann => "neumann",
        BcKind::Robin => "robin",
        BcKind::MassFlowInlet => "mass-flow-inlet",
        BcKind::PressureOutlet => "pressure-outlet",
        BcKind::WallNoSlip => "wall-no-slip",
        BcKind::WallSlip => "wall-slip",
        BcKind::Traction => "traction",
        BcKind::MagneticVectorPotential => "magnetic-vector-potential",
        BcKind::NormalMagneticFluxDensity => "normal-magnetic-flux-density",
        BcKind::ElectricPotential => "electric-potential",
        BcKind::NormalCurrentDensity => "normal-current-density",
        BcKind::SpeciesAmountFlux => "species-amount-flux",
        BcKind::SpeciesMassFlux => "species-mass-flux",
        BcKind::GasCharacteristicInlet => "gas-characteristic-inlet",
        BcKind::GasCharacteristicOutlet => "gas-characteristic-outlet",
    }
}

fn w_bc(out: &mut String, bc: &BoundaryCondition) {
    out.push_str("(bc ");
    w_str(out, &bc.region);
    let _ = write!(
        out,
        " {} {} {} ",
        physics_tag(bc.physics),
        kind_tag(bc.kind),
        bc.frame
    );
    match &bc.value {
        None => out.push_str("none"),
        Some(BcValue::Uniform(q)) => {
            out.push_str("(uniform ");
            w_qty(out, q);
            out.push(')');
        }
        Some(BcValue::Signal(s)) => {
            out.push_str("(signal ");
            w_signal(out, s);
            out.push(')');
        }
        Some(BcValue::Profile(p)) => w_profile(out, "profile", p),
        Some(BcValue::Typed(payload)) => w_payload_hex(out, payload),
    }
    out.push(' ');
    match bc.compatibility {
        None => out.push_str("none"),
        Some(Compat::Incompressible) => out.push_str("incompressible"),
    }
    out.push(')');
}

fn w_frame(out: &mut String, f: &Frame) {
    let _ = write!(out, "(frame {} ", f.id.0);
    w_str(out, &f.name);
    let _ = write!(out, " {} ", f.parent.0);
    match &f.motion {
        FrameMotion::Fixed {
            orientation,
            translation,
        } => {
            let _ = write!(
                out,
                "(fixed (quat {} {} {} {}) ",
                canonical_float(orientation.w),
                canonical_float(orientation.x),
                canonical_float(orientation.y),
                canonical_float(orientation.z)
            );
            w_vec3(out, *translation);
            out.push(')');
        }
        FrameMotion::Rotating { axis, center, rate } => {
            let _ = write!(
                out,
                "(rotating (vec {} {} {}) ",
                canonical_float(axis[0]),
                canonical_float(axis[1]),
                canonical_float(axis[2])
            );
            w_vec3(out, *center);
            out.push(' ');
            w_qty(out, rate);
            out.push(')');
        }
        FrameMotion::Tilt {
            axis,
            center,
            angle,
        } => {
            let _ = write!(
                out,
                "(tilt (vec {} {} {}) ",
                canonical_float(axis[0]),
                canonical_float(axis[1]),
                canonical_float(axis[2])
            );
            w_vec3(out, *center);
            out.push(' ');
            w_signal(out, angle);
            out.push(')');
        }
    }
    out.push(')');
}

fn w_ensemble(out: &mut String, e: &StochasticEnsemble) {
    out.push_str("(ensemble ");
    w_str(out, &e.name);
    let _ = write!(out, " {} {} ", e.seed, e.members);
    w_qty(out, &e.duration);
    out.push(' ');
    w_qty(out, &e.dt);
    out.push(' ');
    match &e.model {
        SpectrumModel::Dryden {
            sigma,
            length_scale,
            mean_speed,
        } => {
            out.push_str("(dryden ");
            w_qty(out, sigma);
            out.push(' ');
            w_qty(out, length_scale);
            out.push(' ');
            w_qty(out, mean_speed);
            out.push(')');
        }
        SpectrumModel::KanaiTajimi {
            s0,
            omega_g,
            zeta_g,
        } => {
            let _ = write!(out, "(kanai-tajimi {} ", canonical_float(*s0));
            w_qty(out, omega_g);
            let _ = write!(out, " {})", canonical_float(*zeta_g));
        }
        SpectrumModel::CarreauBand {
            eta_zero,
            eta_inf,
            lambda,
            n,
        } => {
            out.push_str("(carreau");
            for q in [
                &eta_zero[0],
                &eta_zero[1],
                &eta_inf[0],
                &eta_inf[1],
                &lambda[0],
                &lambda[1],
            ] {
                out.push(' ');
                w_qty(out, q);
            }
            let _ = write!(out, " {} {})", canonical_float(n[0]), canonical_float(n[1]));
        }
    }
    out.push(')');
}

/// Serialize a scenario to canonical, byte-stable IR text.
#[must_use]
pub fn write_ir(s: &Scenario) -> String {
    let mut out = String::with_capacity(1024);
    let _ = write!(out, "(scenario :version {SCENARIO_IR_VERSION} ");
    w_str(&mut out, &s.name);
    let _ = write!(out, " {} (environment ", s.seed);
    for g in &s.environment.gravity {
        w_qty(&mut out, g);
        out.push(' ');
    }
    w_qty(&mut out, &s.environment.ambient_temperature);
    out.push(' ');
    w_qty(&mut out, &s.environment.ambient_pressure);
    out.push_str(") (frames");
    for f in &s.frames.frames {
        out.push(' ');
        w_frame(&mut out, f);
    }
    out.push_str(") (bcs");
    for bc in &s.base_bcs {
        out.push(' ');
        w_bc(&mut out, bc);
    }
    out.push_str(") (cases");
    for case in &s.cases {
        out.push_str(" (case ");
        w_str(&mut out, &case.name);
        for bc in &case.bcs {
            out.push(' ');
            w_bc(&mut out, bc);
        }
        out.push(')');
    }
    out.push_str(") (combos");
    for combo in &s.combinations {
        out.push_str(" (combo ");
        w_str(&mut out, &combo.name);
        for (case, factor) in &combo.terms {
            out.push_str(" (term ");
            w_str(&mut out, case);
            let _ = write!(out, " {})", canonical_float(*factor));
        }
        out.push(')');
    }
    out.push_str(") (ensembles");
    for e in &s.ensembles {
        out.push(' ');
        w_ensemble(&mut out, e);
    }
    out.push_str(") (contacts");
    for c in &s.contacts {
        out.push_str(" (contact ");
        w_str(&mut out, &c.region_a);
        out.push(' ');
        w_str(&mut out, &c.region_b);
        out.push(' ');
        match c.model {
            ContactModel::Frictionless => out.push_str("frictionless"),
            ContactModel::Tied => out.push_str("tied"),
            ContactModel::Coulomb { mu_s, mu_k } => {
                let _ = write!(
                    out,
                    "(coulomb {} {})",
                    canonical_float(mu_s),
                    canonical_float(mu_k)
                );
            }
        }
        out.push(')');
    }
    out.push_str("))");
    out
}

// ---------------------------------------------------------------- parser

#[derive(Debug, Clone, PartialEq)]
struct Sx {
    span: IrSourceSpan,
    kind: SxKind,
}

#[derive(Debug, Clone, PartialEq)]
enum SxKind {
    Atom(String),
    Str(String),
    List(Vec<Sx>),
}

/// Explicit resource budget for scenario-IR parsing.
///
/// Limits are checked before recursion or growth. `max_depth` counts the root
/// form as depth one; `max_list_items` applies independently to every list;
/// `max_atom_bytes` applies to both atoms and decoded string contents.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IrParseBudget {
    /// Maximum input size in bytes.
    pub max_bytes: usize,
    /// Maximum recursive form depth, including the root.
    pub max_depth: usize,
    /// Maximum total atom/string/list nodes.
    pub max_nodes: usize,
    /// Maximum bytes in one atom or decoded string.
    pub max_atom_bytes: usize,
    /// Maximum direct children in one list.
    pub max_list_items: usize,
}

/// Conservative default used by [`parse_ir`].
pub const DEFAULT_IR_PARSE_BUDGET: IrParseBudget = IrParseBudget {
    max_bytes: 16 * 1024 * 1024,
    max_depth: 128,
    max_nodes: 1_000_000,
    // A typed payload is one hex string. Keep the per-atom ceiling equal to
    // the already-bounded total input ceiling so canonical payloads do not hit
    // an unrelated smaller limit after `write_ir`; callers may still tighten
    // either dimension explicitly.
    max_atom_bytes: 16 * 1024 * 1024,
    max_list_items: 100_000,
};

/// Absolute recursion ceiling for the recursive-descent implementation.
/// Callers may tighten but not raise this safety boundary.
pub const MAX_IR_PARSE_DEPTH: usize = 256;

impl Default for IrParseBudget {
    fn default() -> Self {
        DEFAULT_IR_PARSE_BUDGET
    }
}

fn err(at: usize, what: &str) -> ScenarioError {
    ScenarioError::Parse {
        span: IrSourceSpan { start: at, end: at },
        path: String::new(),
        what: what.to_string(),
    }
}

fn err_span(span: IrSourceSpan, what: &str) -> ScenarioError {
    ScenarioError::Parse {
        span,
        path: String::new(),
        what: what.to_string(),
    }
}

fn err_node(node: &Sx, what: &str) -> ScenarioError {
    err_span(node.span, what)
}

impl ScenarioError {
    fn prepend_path(mut self, prefix: &str) -> Self {
        match &mut self {
            Self::Parse { path, .. } | Self::ReservedBoundaryRole { path, .. } => {
                path.insert_str(0, prefix);
            }
            Self::Dimensions { .. } | Self::Frame { .. } | Self::Evaluate { .. } => {}
        }
        self
    }

    fn rooted(mut self) -> Self {
        match &mut self {
            Self::Parse { path, .. } | Self::ReservedBoundaryRole { path, .. } => {
                path.insert(0, '$');
            }
            Self::Dimensions { .. } | Self::Frame { .. } | Self::Evaluate { .. } => {}
        }
        self
    }
}

trait ParsePathExt<T> {
    fn field(self, field: &str) -> Result<T, ScenarioError>;
    fn index(self, index: usize) -> Result<T, ScenarioError>;
}

impl<T> ParsePathExt<T> for Result<T, ScenarioError> {
    fn field(self, field: &str) -> Result<T, ScenarioError> {
        self.map_err(|error| error.prepend_path(&format!(".{field}")))
    }

    fn index(self, index: usize) -> Result<T, ScenarioError> {
        self.map_err(|error| error.prepend_path(&format!("[{index}]")))
    }
}

struct SxParser<'a> {
    bytes: &'a [u8],
    pos: usize,
    budget: IrParseBudget,
    nodes: usize,
}

impl<'a> SxParser<'a> {
    fn new(text: &'a str, budget: IrParseBudget) -> Result<Self, ScenarioError> {
        if budget.max_depth > MAX_IR_PARSE_DEPTH {
            return Err(err(
                0,
                &format!(
                    "configured IR depth budget {} exceeds hard safety limit {MAX_IR_PARSE_DEPTH}",
                    budget.max_depth
                ),
            ));
        }
        if text.len() > budget.max_bytes {
            return Err(err_span(
                IrSourceSpan {
                    start: budget.max_bytes.min(text.len()),
                    end: text.len(),
                },
                &format!(
                    "IR byte budget exceeded: {} bytes > {}",
                    text.len(),
                    budget.max_bytes
                ),
            ));
        }
        Ok(Self {
            bytes: text.as_bytes(),
            pos: 0,
            budget,
            nodes: 0,
        })
    }

    fn parse(mut self) -> Result<Sx, ScenarioError> {
        let root = self.parse_one(1)?;
        self.skip_ws();
        if self.pos != self.bytes.len() {
            return Err(err_span(
                IrSourceSpan {
                    start: self.pos,
                    end: self.bytes.len(),
                },
                "trailing bytes after the scenario form",
            ));
        }
        Ok(root)
    }

    fn skip_ws(&mut self) {
        while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }

    fn admit_node(&mut self, depth: usize) -> Result<(), ScenarioError> {
        if depth > self.budget.max_depth {
            return Err(err(
                self.pos,
                &format!(
                    "IR depth budget exceeded: depth {depth} > {}",
                    self.budget.max_depth
                ),
            ));
        }
        self.nodes = self
            .nodes
            .checked_add(1)
            .ok_or_else(|| err(self.pos, "IR node counter overflowed"))?;
        if self.nodes > self.budget.max_nodes {
            return Err(err(
                self.pos,
                &format!(
                    "IR node budget exceeded: node {} > {}",
                    self.nodes, self.budget.max_nodes
                ),
            ));
        }
        Ok(())
    }

    fn reserve<T>(&self, values: &mut Vec<T>, additional: usize) -> Result<(), ScenarioError> {
        values.try_reserve(additional).map_err(|allocation_error| {
            err(
                self.pos,
                &format!("IR parser allocation refused: {allocation_error}"),
            )
        })
    }

    fn parse_one(&mut self, depth: usize) -> Result<Sx, ScenarioError> {
        self.skip_ws();
        let start = self.pos;
        self.admit_node(depth)?;
        match self.bytes.get(self.pos) {
            None => Err(err(self.pos, "unexpected end of input")),
            Some(b'(') => {
                self.pos += 1;
                let mut items = Vec::new();
                loop {
                    self.skip_ws();
                    match self.bytes.get(self.pos) {
                        None => return Err(err(self.pos, "unclosed list")),
                        Some(b')') => {
                            self.pos += 1;
                            return Ok(Sx {
                                span: IrSourceSpan {
                                    start,
                                    end: self.pos,
                                },
                                kind: SxKind::List(items),
                            });
                        }
                        _ => {
                            if items.len() >= self.budget.max_list_items {
                                return Err(err(
                                    self.pos,
                                    &format!(
                                        "IR list-item budget exceeded: more than {} children",
                                        self.budget.max_list_items
                                    ),
                                ));
                            }
                            self.reserve(&mut items, 1)?;
                            items.push(self.parse_one(depth + 1)?);
                        }
                    }
                }
            }
            Some(b')') => Err(err_span(
                IrSourceSpan {
                    start: self.pos,
                    end: self.pos + 1,
                },
                "unexpected ')'",
            )),
            Some(b'"') => {
                self.pos += 1;
                // Accumulate RAW BYTES and decode UTF-8 once at the end. The writer
                // (`w_str`) emits each `char` UTF-8-encoded and escapes only the
                // ASCII bytes `"`/`\`; pushing `byte as char` here would Latin-1-
                // decode, splitting every multi-byte code point (e.g. `é` → `Ã©`)
                // and breaking round-trip losslessness for non-ASCII names. Escape
                // bytes and the closing quote never collide with a UTF-8
                // continuation byte (those are all ≥ 0x80, `"`/`\` are ASCII).
                let mut buf: Vec<u8> = Vec::new();
                loop {
                    match self.bytes.get(self.pos) {
                        None => return Err(err(self.pos, "unterminated string")),
                        Some(b'\\') => {
                            let c = *self.bytes.get(self.pos + 1).ok_or_else(|| {
                                err_span(
                                    IrSourceSpan {
                                        start: self.pos,
                                        end: self.pos + 1,
                                    },
                                    "bad escape",
                                )
                            })?;
                            if !matches!(c, b'"' | b'\\') {
                                return Err(err_span(
                                    IrSourceSpan {
                                        start: self.pos,
                                        end: self.pos + 2,
                                    },
                                    "unsupported string escape; only quote and backslash escapes are canonical",
                                ));
                            }
                            if buf.len() >= self.budget.max_atom_bytes {
                                return Err(err_span(
                                    IrSourceSpan {
                                        start: self.pos,
                                        end: self.pos + 1,
                                    },
                                    "IR string exceeds atom-byte budget",
                                ));
                            }
                            self.reserve(&mut buf, 1)?;
                            buf.push(c);
                            self.pos += 2;
                        }
                        Some(b'"') => {
                            self.pos += 1;
                            let s = String::from_utf8(buf).map_err(|_| {
                                err_span(
                                    IrSourceSpan {
                                        start,
                                        end: self.pos,
                                    },
                                    "string is not valid UTF-8",
                                )
                            })?;
                            return Ok(Sx {
                                span: IrSourceSpan {
                                    start,
                                    end: self.pos,
                                },
                                kind: SxKind::Str(s),
                            });
                        }
                        Some(&c) => {
                            if buf.len() >= self.budget.max_atom_bytes {
                                return Err(err_span(
                                    IrSourceSpan {
                                        start: self.pos,
                                        end: self.pos + 1,
                                    },
                                    "IR string exceeds atom-byte budget",
                                ));
                            }
                            self.reserve(&mut buf, 1)?;
                            buf.push(c);
                            self.pos += 1;
                        }
                    }
                }
            }
            Some(_) => {
                while self.pos < self.bytes.len()
                    && !self.bytes[self.pos].is_ascii_whitespace()
                    && self.bytes[self.pos] != b'('
                    && self.bytes[self.pos] != b')'
                {
                    self.pos += 1;
                    if self.pos - start > self.budget.max_atom_bytes {
                        return Err(err_span(
                            IrSourceSpan {
                                start,
                                end: self.pos,
                            },
                            "IR atom exceeds atom-byte budget",
                        ));
                    }
                }
                let atom = std::str::from_utf8(&self.bytes[start..self.pos]).map_err(|_| {
                    err_span(
                        IrSourceSpan {
                            start,
                            end: self.pos,
                        },
                        "atom is not valid UTF-8",
                    )
                })?;
                let mut owned = String::new();
                owned
                    .try_reserve_exact(atom.len())
                    .map_err(|allocation_error| {
                        err_span(
                            IrSourceSpan {
                                start,
                                end: self.pos,
                            },
                            &format!("IR atom allocation refused: {allocation_error}"),
                        )
                    })?;
                owned.push_str(atom);
                Ok(Sx {
                    span: IrSourceSpan {
                        start,
                        end: self.pos,
                    },
                    kind: SxKind::Atom(owned),
                })
            }
        }
    }
}

fn parse_sx(text: &str, budget: IrParseBudget) -> Result<Sx, ScenarioError> {
    SxParser::new(text, budget)
        .and_then(SxParser::parse)
        .map_err(ScenarioError::rooted)
}

// -------------------------------------------------------------- decoding

fn as_list<'a>(sx: &'a Sx, head: &str) -> Result<&'a [Sx], ScenarioError> {
    if let SxKind::List(items) = &sx.kind
        && let Some(Sx {
            kind: SxKind::Atom(a),
            ..
        }) = items.first()
        && a == head
    {
        return Ok(&items[1..]);
    }
    Err(err_node(sx, &format!("expected ({head} ...) form")))
}

fn as_f64(sx: &Sx) -> Result<f64, ScenarioError> {
    if let SxKind::Atom(a) = &sx.kind
        && let Ok(v) = a.parse::<f64>()
        && v.is_finite()
    {
        return Ok(canonical_float(v));
    }
    Err(err_node(sx, "expected a finite number"))
}

fn reserve_decoded_string(
    value: &mut String,
    additional: usize,
    span: IrSourceSpan,
) -> Result<(), ScenarioError> {
    value
        .try_reserve_exact(additional)
        .map_err(|allocation_error| {
            err_span(
                span,
                &format!(
                    "IR decoded string allocation for {additional} bytes was refused: {allocation_error}"
                ),
            )
        })
}

fn as_str(sx: &Sx) -> Result<String, ScenarioError> {
    if let SxKind::Str(s) = &sx.kind {
        let mut owned = String::new();
        reserve_decoded_string(&mut owned, s.len(), sx.span)?;
        owned.push_str(s);
        return Ok(owned);
    }
    Err(err_node(sx, "expected a string"))
}

fn as_str_ref(sx: &Sx) -> Result<&str, ScenarioError> {
    if let SxKind::Str(value) = &sx.kind {
        Ok(value)
    } else {
        Err(err_node(sx, "expected a string"))
    }
}

fn lower_hex_nibble(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        _ => None,
    }
}

fn decode_payload_hex(source: &Sx) -> Result<Vec<u8>, ScenarioError> {
    let text = as_str_ref(source)?;
    if text.len() % 2 != 0 {
        return Err(err_node(
            source,
            "typed payload hex must contain complete byte pairs",
        ));
    }
    let byte_count = text.len() / 2;
    if byte_count > MAX_PAYLOAD_WIRE_BYTES {
        return Err(err_node(
            source,
            &format!(
                "typed payload requests {byte_count} bytes; hard V1 limit is {MAX_PAYLOAD_WIRE_BYTES}"
            ),
        ));
    }
    let mut decoded = Vec::new();
    decoded
        .try_reserve_exact(byte_count)
        .map_err(|allocation_error| {
            err_node(
                source,
                &format!(
                    "typed payload allocation for {byte_count} bytes was refused: {allocation_error}"
                ),
            )
        })?;
    for (index, pair) in text.as_bytes().chunks_exact(2).enumerate() {
        let high = lower_hex_nibble(pair[0]).ok_or_else(|| {
            err_span(
                IrSourceSpan {
                    start: source.span.start + 1 + index * 2,
                    end: source.span.start + 2 + index * 2,
                },
                "typed payload hex must use canonical lowercase digits",
            )
        })?;
        let low = lower_hex_nibble(pair[1]).ok_or_else(|| {
            err_span(
                IrSourceSpan {
                    start: source.span.start + 2 + index * 2,
                    end: source.span.start + 3 + index * 2,
                },
                "typed payload hex must use canonical lowercase digits",
            )
        })?;
        decoded.push((high << 4) | low);
    }
    Ok(decoded)
}

fn as_atom(sx: &Sx) -> Result<&str, ScenarioError> {
    if let SxKind::Atom(a) = &sx.kind {
        return Ok(a);
    }
    Err(err_node(sx, "expected an atom"))
}

fn as_u64(sx: &Sx) -> Result<u64, ScenarioError> {
    if let SxKind::Atom(a) = &sx.kind
        && let Ok(v) = a.parse::<u64>()
    {
        return Ok(v);
    }
    Err(err_node(sx, "expected an unsigned integer"))
}

fn as_u32(sx: &Sx) -> Result<u32, ScenarioError> {
    if let SxKind::Atom(a) = &sx.kind
        && let Ok(v) = a.parse::<u32>()
    {
        return Ok(v);
    }
    Err(err_node(sx, "expected a u32"))
}

fn as_i8(sx: &Sx) -> Result<i8, ScenarioError> {
    if let SxKind::Atom(a) = &sx.kind
        && let Ok(v) = a.parse::<i8>()
    {
        return Ok(v);
    }
    Err(err_node(sx, "expected a small integer exponent"))
}

fn as_dims_at(container: &Sx, items: &[Sx], wire: DimensionWire) -> Result<Dims, ScenarioError> {
    let expected = match wire {
        DimensionWire::LegacyFive => 5,
        DimensionWire::CanonicalSix => 6,
    };
    if items.len() != expected {
        return Err(err_node(
            container,
            &format!("scenario IR needs exactly {expected} dimension exponents"),
        ));
    }
    let mut d = [0i8; 6];
    for (index, (slot, sx)) in d.iter_mut().zip(items).enumerate() {
        *slot = as_i8(sx).index(index)?;
    }
    Ok(Dims(d))
}

fn as_qty(sx: &Sx, wire: DimensionWire) -> Result<QtyAny, ScenarioError> {
    let items = as_list(sx, "qty")?;
    let expected = match wire {
        DimensionWire::LegacyFive => 6,
        DimensionWire::CanonicalSix => 7,
    };
    if items.len() != expected {
        return Err(err_node(
            sx,
            &format!("qty needs a value plus {} exponents", expected - 1),
        ));
    }
    Ok(QtyAny::new(
        as_f64(&items[0]).field("value")?,
        as_dims_at(sx, &items[1..], wire).field("dims")?,
    ))
}

fn as_dims(sx: &Sx, wire: DimensionWire) -> Result<Dims, ScenarioError> {
    as_dims_at(sx, as_list(sx, "dims")?, wire)
}

fn reserve_decoded<T>(
    values: &mut Vec<T>,
    additional: usize,
    resource: &str,
    span: IrSourceSpan,
) -> Result<(), ScenarioError> {
    values
        .try_reserve_exact(additional)
        .map_err(|allocation_error| {
            err_span(
                span,
                &format!(
                    "IR decoded {resource} allocation for {additional} elements was refused: {allocation_error}"
                ),
            )
        })
}

fn decode_items<T>(
    container: &Sx,
    items: &[Sx],
    resource: &str,
    mut decode: impl FnMut(&Sx) -> Result<T, ScenarioError>,
) -> Result<Vec<T>, ScenarioError> {
    let mut values = Vec::new();
    reserve_decoded(&mut values, items.len(), resource, container.span)?;
    for (index, item) in items.iter().enumerate() {
        values.push(decode(item).index(index)?);
    }
    Ok(values)
}

fn as_floats(sx: &Sx, head: &str) -> Result<Vec<f64>, ScenarioError> {
    decode_items(sx, as_list(sx, head)?, head, as_f64)
}

fn as_vec3(sx: &Sx) -> Result<Vec3, ScenarioError> {
    let items = as_list(sx, "vec")?;
    if items.len() != 3 {
        return Err(err_node(sx, "vec needs three components"));
    }
    Ok(Vec3::new(
        as_f64(&items[0]).field("x")?,
        as_f64(&items[1]).field("y")?,
        as_f64(&items[2]).field("z")?,
    ))
}

fn as_profile(sx: &Sx, head: &str, wire: DimensionWire) -> Result<ChebProfile, ScenarioError> {
    let items = as_list(sx, head)?;
    if items.len() != 4 {
        return Err(err_node(sx, "profile needs dims, domain, coeffs"));
    }
    let dims = as_dims(&items[0], wire).field("dims")?;
    let a = as_f64(&items[1]).field("domain_start")?;
    let b = as_f64(&items[2]).field("domain_end")?;
    let coeffs = as_floats(&items[3], "coeffs").field("coeffs")?;
    if !(a < b) {
        return Err(err_span(
            IrSourceSpan {
                start: items[1].span.start,
                end: items[2].span.end,
            },
            "profile domain must satisfy finite a < b",
        )
        .prepend_path(".domain"));
    }
    if coeffs.is_empty() {
        return Err(
            err_node(&items[3], "profile needs at least one finite coefficient")
                .prepend_path(".coeffs"),
        );
    }
    // `as_f64` already rejects non-finite coefficients. Check every public
    // constructor precondition before calling `fs_cheb`, whose infallible
    // constructor deliberately asserts them.
    Ok(ChebProfile {
        cheb: Cheb1::from_coeffs(a, b, coeffs),
        dims,
    })
}

fn as_signal(sx: &Sx, wire: DimensionWire) -> Result<TimeSignal, ScenarioError> {
    let SxKind::List(items) = &sx.kind else {
        return Err(err_node(sx, "expected a signal form"));
    };
    let head_node = items
        .first()
        .ok_or_else(|| err_node(sx, "empty signal"))
        .field("kind")?;
    let head = as_atom(head_node).field("kind")?;
    let rest = &items[1..];
    match head {
        "constant" => {
            if rest.len() != 1 {
                return Err(err_node(sx, "constant needs a value"));
            }
            Ok(TimeSignal::Constant(as_qty(&rest[0], wire).field("value")?))
        }
        "ramp" => {
            if rest.len() != 4 {
                return Err(err_node(sx, "ramp needs t0 t1 from to"));
            }
            Ok(TimeSignal::Ramp {
                t_start: as_f64(&rest[0]).field("t_start")?,
                t_end: as_f64(&rest[1]).field("t_end")?,
                from: as_qty(&rest[2], wire).field("from")?,
                to: as_qty(&rest[3], wire).field("to")?,
            })
        }
        "table" => {
            if rest.len() != 4 {
                return Err(err_node(sx, "table needs interp dims times values"));
            }
            let interp = match as_atom(&rest[0]).field("interp")? {
                "linear" => Interp::Linear,
                "hold" => Interp::Hold,
                other => {
                    return Err(err_node(&rest[0], &format!("unknown interp {other:?}"))
                        .prepend_path(".interp"));
                }
            };
            Ok(TimeSignal::Table {
                interp,
                dims: as_dims(&rest[1], wire).field("dims")?,
                times: as_floats(&rest[2], "times").field("times")?,
                values: as_floats(&rest[3], "values").field("values")?,
            })
        }
        "chebfun" => Ok(TimeSignal::Chebfun(as_profile(sx, "chebfun", wire)?)),
        other => {
            Err(err_node(head_node, &format!("unknown signal {other:?}")).prepend_path(".kind"))
        }
    }
}

fn as_physics(sx: &Sx) -> Result<Physics, ScenarioError> {
    let a = as_atom(sx)?;
    match a {
        "incompressible-flow" => Ok(Physics::IncompressibleFlow),
        "thermal" => Ok(Physics::Thermal),
        "elasticity" => Ok(Physics::Elasticity),
        "magnetics" => Ok(Physics::Magnetics),
        "electrics" => Ok(Physics::Electrics),
        "gas-exchange" => Ok(Physics::GasExchange),
        other => Err(err_node(sx, &format!("unknown physics {other:?}"))),
    }
}

fn reserved_machine_role(a: &str) -> Option<&'static str> {
    match a {
        "joint" => Some("joint"),
        "terminal" => Some("terminal"),
        "controller" => Some("controller"),
        "reset" => Some("reset"),
        _ => None,
    }
}

fn as_kind(sx: &Sx) -> Result<BcKind, ScenarioError> {
    let a = as_atom(sx)?;
    if let Some(role) = reserved_machine_role(a) {
        return Err(ScenarioError::ReservedBoundaryRole {
            role,
            span: sx.span,
            path: String::new(),
        });
    }
    match a {
        "dirichlet" => Ok(BcKind::Dirichlet),
        "neumann" => Ok(BcKind::Neumann),
        "robin" => Ok(BcKind::Robin),
        "mass-flow-inlet" => Ok(BcKind::MassFlowInlet),
        "pressure-outlet" => Ok(BcKind::PressureOutlet),
        "wall-no-slip" => Ok(BcKind::WallNoSlip),
        "wall-slip" => Ok(BcKind::WallSlip),
        "traction" => Ok(BcKind::Traction),
        "magnetic-vector-potential" => Ok(BcKind::MagneticVectorPotential),
        "normal-magnetic-flux-density" => Ok(BcKind::NormalMagneticFluxDensity),
        "electric-potential" => Ok(BcKind::ElectricPotential),
        "normal-current-density" => Ok(BcKind::NormalCurrentDensity),
        "species-amount-flux" => Ok(BcKind::SpeciesAmountFlux),
        "species-mass-flux" => Ok(BcKind::SpeciesMassFlux),
        "gas-characteristic-inlet" => Ok(BcKind::GasCharacteristicInlet),
        "gas-characteristic-outlet" => Ok(BcKind::GasCharacteristicOutlet),
        other => Err(err_node(sx, &format!("unknown bc kind {other:?}"))),
    }
}

fn as_typed_payload(sx: &Sx, wire: DimensionWire) -> Result<Payload, ScenarioError> {
    if wire != DimensionWire::CanonicalSix {
        return Err(err_node(sx, "typed payloads require scenario IR version 2"));
    }
    let inner = as_list(sx, "typed")?;
    if inner.len() != 3 {
        return Err(err_node(
            sx,
            "typed payload needs :version, version number, and canonical hex bytes",
        ));
    }
    if as_atom(&inner[0]).field("version_tag")? != ":version" {
        return Err(
            err_node(&inner[0], "typed payload version tag must be :version")
                .prepend_path(".version_tag"),
        );
    }
    let version = as_u32(&inner[1]).field("version")?;
    if version != u32::from(PAYLOAD_WIRE_VERSION) {
        return Err(err_node(
            &inner[1],
            &format!(
                "unsupported typed payload version {version}; supported version is {PAYLOAD_WIRE_VERSION}"
            ),
        )
        .prepend_path(".version"));
    }
    let bytes = decode_payload_hex(&inner[2]).field("bytes")?;
    let limits = PayloadDecodeLimits {
        max_bytes: bytes.len(),
        ..PayloadDecodeLimits::DEFAULT
    };
    decode_payload_with_limits(&bytes, limits).map_err(|error| {
        err_node(&inner[2], &format!("typed payload refused: {error}")).prepend_path(".bytes")
    })
}

fn as_bc(sx: &Sx, wire: DimensionWire) -> Result<BoundaryCondition, ScenarioError> {
    let items = as_list(sx, "bc")?;
    if items.len() != 6 {
        return Err(err_node(
            sx,
            "bc needs region physics kind frame value compat",
        ));
    }
    let value = match &items[4].kind {
        SxKind::Atom(a) if a == "none" => None,
        SxKind::List(inner) => {
            let head_node = inner
                .first()
                .ok_or_else(|| err_node(&items[4], "empty bc value"))
                .field("value.kind")?;
            let head = as_atom(head_node).field("value.kind")?;
            match head {
                "uniform" => {
                    if inner.len() != 2 {
                        return Err(err_node(&items[4], "uniform bc value needs a quantity")
                            .prepend_path(".value"));
                    }
                    Some(BcValue::Uniform(
                        as_qty(&inner[1], wire).field("value.uniform")?,
                    ))
                }
                "signal" => {
                    if inner.len() != 2 {
                        return Err(err_node(&items[4], "signal bc value needs a signal form")
                            .prepend_path(".value"));
                    }
                    Some(BcValue::Signal(
                        as_signal(&inner[1], wire).field("value.signal")?,
                    ))
                }
                "profile" => Some(BcValue::Profile(
                    as_profile(&items[4], "profile", wire).field("value.profile")?,
                )),
                "typed" => Some(BcValue::Typed(
                    as_typed_payload(&items[4], wire).field("value.typed")?,
                )),
                other => {
                    return Err(err_node(head_node, &format!("unknown bc value {other:?}"))
                        .prepend_path(".value.kind"));
                }
            }
        }
        _ => return Err(err_node(&items[4], "bad bc value").prepend_path(".value")),
    };
    let compatibility = match as_atom(&items[5]).field("compatibility")? {
        "none" => None,
        "incompressible" => Some(Compat::Incompressible),
        other => {
            return Err(err_node(&items[5], &format!("unknown compat {other:?}"))
                .prepend_path(".compatibility"));
        }
    };
    Ok(BoundaryCondition {
        region: as_str(&items[0]).field("region")?,
        physics: as_physics(&items[1]).field("physics")?,
        kind: as_kind(&items[2]).field("kind")?,
        frame: as_u32(&items[3]).field("frame")?,
        value,
        compatibility,
    })
}

fn as_case(sx: &Sx, wire: DimensionWire) -> Result<LoadCase, ScenarioError> {
    let items = as_list(sx, "case")?;
    let name_node = items
        .first()
        .ok_or_else(|| err_node(sx, "case needs a name"))
        .field("name")?;
    let name = as_str(name_node).field("name")?;
    let bcs = decode_items(sx, &items[1..], "case boundary conditions", |bc| {
        as_bc(bc, wire)
    })
    .field("bcs")?;
    Ok(LoadCase { name, bcs })
}

fn as_combination(sx: &Sx) -> Result<Combination, ScenarioError> {
    let items = as_list(sx, "combo")?;
    let name_node = items
        .first()
        .ok_or_else(|| err_node(sx, "combo needs a name"))
        .field("name")?;
    let name = as_str(name_node).field("name")?;
    let terms = decode_items(sx, &items[1..], "combination terms", |term| {
        let term_items = as_list(term, "term")?;
        if term_items.len() != 2 {
            return Err(err_node(term, "term needs case + factor"));
        }
        Ok((
            as_str(&term_items[0]).field("case")?,
            as_f64(&term_items[1]).field("factor")?,
        ))
    })
    .field("terms")?;
    Ok(Combination { name, terms })
}

fn as_frame(sx: &Sx, wire: DimensionWire) -> Result<Frame, ScenarioError> {
    let items = as_list(sx, "frame")?;
    if items.len() != 4 {
        return Err(err_node(sx, "frame needs id name parent motion"));
    }
    let SxKind::List(motion_items) = &items[3].kind else {
        return Err(err_node(&items[3], "bad frame motion").prepend_path(".motion"));
    };
    let head_node = motion_items
        .first()
        .ok_or_else(|| err_node(&items[3], "empty motion"))
        .field("motion.kind")?;
    let head = as_atom(head_node).field("motion.kind")?;
    let rest = &motion_items[1..];
    let motion = match head {
        "fixed" => {
            if rest.len() != 2 {
                return Err(
                    err_node(&items[3], "fixed motion needs orientation and translation")
                        .prepend_path(".motion"),
                );
            }
            let q_items = as_list(&rest[0], "quat").field("motion.orientation")?;
            if q_items.len() != 4 {
                return Err(err_node(&rest[0], "quat needs four components")
                    .prepend_path(".motion.orientation"));
            }
            FrameMotion::Fixed {
                orientation: Quat {
                    w: as_f64(&q_items[0]).field("motion.orientation.w")?,
                    x: as_f64(&q_items[1]).field("motion.orientation.x")?,
                    y: as_f64(&q_items[2]).field("motion.orientation.y")?,
                    z: as_f64(&q_items[3]).field("motion.orientation.z")?,
                },
                translation: as_vec3(&rest[1]).field("motion.translation")?,
            }
        }
        "rotating" => {
            if rest.len() != 3 {
                return Err(
                    err_node(&items[3], "rotating motion needs axis center rate")
                        .prepend_path(".motion"),
                );
            }
            let axis = as_vec3(&rest[0]).field("motion.axis")?;
            FrameMotion::Rotating {
                axis: [axis.x, axis.y, axis.z],
                center: as_vec3(&rest[1]).field("motion.center")?,
                rate: as_qty(&rest[2], wire).field("motion.rate")?,
            }
        }
        "tilt" => {
            if rest.len() != 3 {
                return Err(err_node(&items[3], "tilt motion needs axis center angle")
                    .prepend_path(".motion"));
            }
            let axis = as_vec3(&rest[0]).field("motion.axis")?;
            FrameMotion::Tilt {
                axis: [axis.x, axis.y, axis.z],
                center: as_vec3(&rest[1]).field("motion.center")?,
                angle: as_signal(&rest[2], wire).field("motion.angle")?,
            }
        }
        other => {
            return Err(err_node(head_node, &format!("unknown motion {other:?}"))
                .prepend_path(".motion.kind"));
        }
    };
    Ok(Frame {
        id: FrameId(as_u32(&items[0]).field("id")?),
        name: as_str(&items[1]).field("name")?,
        parent: FrameId(as_u32(&items[2]).field("parent")?),
        motion,
    })
}

fn as_model(sx: &Sx, wire: DimensionWire) -> Result<SpectrumModel, ScenarioError> {
    let SxKind::List(items) = &sx.kind else {
        return Err(err_node(sx, "expected a spectrum model form"));
    };
    let head_node = items
        .first()
        .ok_or_else(|| err_node(sx, "empty model"))
        .field("kind")?;
    let head = as_atom(head_node).field("kind")?;
    let rest = &items[1..];
    match head {
        "dryden" => {
            if rest.len() != 3 {
                return Err(err_node(sx, "dryden needs sigma length_scale mean_speed"));
            }
            Ok(SpectrumModel::Dryden {
                sigma: as_qty(&rest[0], wire).field("sigma")?,
                length_scale: as_qty(&rest[1], wire).field("length_scale")?,
                mean_speed: as_qty(&rest[2], wire).field("mean_speed")?,
            })
        }
        "kanai-tajimi" => {
            if rest.len() != 3 {
                return Err(err_node(sx, "kanai-tajimi needs s0 omega_g zeta_g"));
            }
            Ok(SpectrumModel::KanaiTajimi {
                s0: as_f64(&rest[0]).field("s0")?,
                omega_g: as_qty(&rest[1], wire).field("omega_g")?,
                zeta_g: as_f64(&rest[2]).field("zeta_g")?,
            })
        }
        "carreau" => {
            if rest.len() != 8 {
                return Err(err_node(sx, "carreau needs six qty bounds + two n bounds"));
            }
            Ok(SpectrumModel::CarreauBand {
                eta_zero: [
                    as_qty(&rest[0], wire).field("eta_zero[0]")?,
                    as_qty(&rest[1], wire).field("eta_zero[1]")?,
                ],
                eta_inf: [
                    as_qty(&rest[2], wire).field("eta_inf[0]")?,
                    as_qty(&rest[3], wire).field("eta_inf[1]")?,
                ],
                lambda: [
                    as_qty(&rest[4], wire).field("lambda[0]")?,
                    as_qty(&rest[5], wire).field("lambda[1]")?,
                ],
                n: [
                    as_f64(&rest[6]).field("n[0]")?,
                    as_f64(&rest[7]).field("n[1]")?,
                ],
            })
        }
        other => {
            Err(err_node(head_node, &format!("unknown model {other:?}")).prepend_path(".kind"))
        }
    }
}

fn as_ensemble(sx: &Sx, wire: DimensionWire) -> Result<StochasticEnsemble, ScenarioError> {
    let items = as_list(sx, "ensemble")?;
    if items.len() != 6 {
        return Err(err_node(
            sx,
            "ensemble needs name seed members duration dt model",
        ));
    }
    Ok(StochasticEnsemble {
        name: as_str(&items[0]).field("name")?,
        seed: as_u64(&items[1]).field("seed")?,
        members: as_u32(&items[2]).field("members")?,
        duration: as_qty(&items[3], wire).field("duration")?,
        dt: as_qty(&items[4], wire).field("dt")?,
        model: as_model(&items[5], wire).field("model")?,
    })
}

fn as_contact(sx: &Sx) -> Result<ContactLaw, ScenarioError> {
    let contact_items = as_list(sx, "contact")?;
    if contact_items.len() != 3 {
        return Err(err_node(sx, "contact needs two regions + model"));
    }
    let model = match &contact_items[2] {
        Sx {
            kind: SxKind::Atom(a),
            ..
        } if a == "frictionless" => ContactModel::Frictionless,
        Sx {
            kind: SxKind::Atom(a),
            ..
        } if a == "tied" => ContactModel::Tied,
        other @ Sx {
            kind: SxKind::List(_),
            ..
        } => {
            let coulomb = as_list(other, "coulomb").field("model")?;
            if coulomb.len() != 2 {
                return Err(err_node(other, "coulomb needs mu_s + mu_k").prepend_path(".model"));
            }
            ContactModel::Coulomb {
                mu_s: as_f64(&coulomb[0]).field("model.mu_static")?,
                mu_k: as_f64(&coulomb[1]).field("model.mu_kinetic")?,
            }
        }
        _ => {
            return Err(err_node(&contact_items[2], "bad contact model").prepend_path(".model"));
        }
    };
    Ok(ContactLaw {
        region_a: as_str(&contact_items[0]).field("region_a")?,
        region_b: as_str(&contact_items[1]).field("region_b")?,
        model,
    })
}

fn as_environment(sx: &Sx, wire: DimensionWire) -> Result<Environment, ScenarioError> {
    let items = as_list(sx, "environment")?;
    if items.len() != 5 {
        return Err(err_node(
            sx,
            "environment needs gravity x3 + temperature + pressure",
        ));
    }
    Ok(Environment {
        gravity: [
            as_qty(&items[0], wire).field("gravity[0]")?,
            as_qty(&items[1], wire).field("gravity[1]")?,
            as_qty(&items[2], wire).field("gravity[2]")?,
        ],
        ambient_temperature: as_qty(&items[3], wire).field("ambient_temperature")?,
        ambient_pressure: as_qty(&items[4], wire).field("ambient_pressure")?,
    })
}

fn scenario_wire_header<'a>(
    root: &Sx,
    root_items: &'a [Sx],
) -> Result<(u32, DimensionWire, &'a [Sx]), ScenarioError> {
    if matches!(root_items.first(), Some(Sx { kind: SxKind::Atom(key), .. }) if key == ":version") {
        let version = root_items
            .get(1)
            .ok_or_else(|| err_node(root, "scenario :version needs a value"))
            .and_then(as_u32)
            .field("version")?;
        let wire = match version {
            LEGACY_SCENARIO_IR_VERSION => DimensionWire::LegacyFive,
            SCENARIO_IR_VERSION => DimensionWire::CanonicalSix,
            _ => {
                return Err(err_node(
                    &root_items[1],
                    &format!(
                        "unsupported scenario IR version {version}; supported versions are 1 and {SCENARIO_IR_VERSION}"
                    ),
                )
                .prepend_path(".version"));
            }
        };
        Ok((version, wire, &root_items[2..]))
    } else {
        Ok((
            LEGACY_SCENARIO_IR_VERSION,
            DimensionWire::LegacyFive,
            root_items,
        ))
    }
}

/// Parse scenario IR under [`DEFAULT_IR_PARSE_BUDGET`].
///
/// The decoded value retains its wire version and any dimension crosswalk or
/// current-version source-canonicalization receipt. Canonical v2 uses six
/// exponents and may embed versioned typed payloads; explicit v1 and the
/// historical unversioned form use five, append `mol = 0`, and refuse typed
/// payload forms.
///
/// # Errors
/// [`ScenarioError::Parse`] for malformed, non-finite, or over-budget input,
/// or [`ScenarioError::ReservedBoundaryRole`] when Machine-IR graph semantics
/// are presented in a boundary-kind slot.
pub fn parse_ir(text: &str) -> Result<DecodedScenario, ScenarioError> {
    parse_ir_with_budget(text, IrParseBudget::default())
}

/// Parse scenario IR under an explicit byte/depth/node/atom/list budget.
///
/// Syntactic resource admission happens before recursive descent or syntax-tree
/// growth; no over-budget input is partially decoded into a [`Scenario`]. The
/// resulting scenario still requires [`Scenario::validate`] for semantic
/// admission.
///
/// # Errors
/// [`ScenarioError::Parse`] for malformed, non-finite, or over-budget input,
/// or [`ScenarioError::ReservedBoundaryRole`] when Machine-IR graph semantics
/// are presented in a boundary-kind slot.
pub fn parse_ir_with_budget(
    text: &str,
    budget: IrParseBudget,
) -> Result<DecodedScenario, ScenarioError> {
    let root = parse_sx(text, budget)?;
    (|| {
        let root_items = as_list(&root, "scenario")?;
        let (source_version, wire, items) = scenario_wire_header(&root, root_items)?;
        if items.len() != 9 {
            return Err(err_node(
                &root,
                "scenario needs name seed + seven sections",
            ));
        }
        let environment = as_environment(&items[2], wire).field("environment")?;
        let frames = FrameTree {
            frames: decode_items(
                &items[3],
                as_list(&items[3], "frames").field("frames")?,
                "frames",
                |frame| as_frame(frame, wire),
            )
            .field("frames")?,
        };
        let base_bcs = decode_items(
            &items[4],
            as_list(&items[4], "bcs").field("base_bcs")?,
            "base boundary conditions",
            |bc| as_bc(bc, wire),
        )
        .field("base_bcs")?;
        let cases = decode_items(
            &items[5],
            as_list(&items[5], "cases").field("cases")?,
            "load cases",
            |case| as_case(case, wire),
        )
        .field("cases")?;
        let combinations = decode_items(
            &items[6],
            as_list(&items[6], "combos").field("combinations")?,
            "load combinations",
            as_combination,
        )
        .field("combinations")?;
        let ensembles = decode_items(
            &items[7],
            as_list(&items[7], "ensembles").field("ensembles")?,
            "stochastic ensembles",
            |ensemble| as_ensemble(ensemble, wire),
        )
        .field("ensembles")?;
        let contacts = decode_items(
            &items[8],
            as_list(&items[8], "contacts").field("contacts")?,
            "contact laws",
            as_contact,
        )
        .field("contacts")?;
        let scenario = Scenario {
            name: as_str(&items[0]).field("name")?,
            seed: as_u64(&items[1]).field("seed")?,
            frames,
            base_bcs,
            cases,
            combinations,
            ensembles,
            contacts,
            environment,
        };
        let canonical = write_ir(&scenario);
        let (dimension_crosswalk, source_canonicalization) =
            if wire == DimensionWire::LegacyFive {
                (
                    Some(DimensionCrosswalkReceipt {
                        source_version,
                        target_version: SCENARIO_IR_VERSION,
                        old_hash: hash_bytes(text.as_bytes()),
                        new_hash: hash_bytes(canonical.as_bytes()),
                        source_width: 5,
                        target_width: 6,
                        rule: FiveToSixRule::AppendMoleZero,
                    }),
                    None,
                )
            } else if text.as_bytes() == canonical.as_bytes() {
                (None, None)
            } else {
                let rule = SourceCanonicalizationRule::for_source_version(source_version)
                    .ok_or_else(|| {
                        err_node(
                            &root,
                            &format!(
                                "scenario IR v{source_version} has no registered source-canonicalization rule"
                            ),
                        )
                        .prepend_path(".version")
                    })?;
                (
                    None,
                    Some(SourceCanonicalizationReceipt {
                        source_version,
                        canonical_version: SCENARIO_IR_VERSION,
                        source_hash: hash_bytes(text.as_bytes()),
                        canonical_hash: hash_bytes(canonical.as_bytes()),
                        rule,
                    }),
                )
            };
        Ok(DecodedScenario {
            scenario,
            source_version,
            dimension_crosswalk,
            source_canonicalization,
        })
    })()
    .map_err(ScenarioError::rooted)
}

/// Round-trip helper for lints/tests: violations if reparse ≠ original.
///
/// The exact already-materialized writer output supplies its own byte and atom
/// authority for this semantic check; depth, node, and list ceilings remain the
/// conservative parser defaults. This keeps resource policy distinct from
/// canonical inversion without granting open-ended syntax-tree growth.
pub fn check_round_trip(s: &Scenario, out: &mut Vec<Violation>) {
    let text = write_ir(s);
    let budget = IrParseBudget {
        max_bytes: text.len(),
        max_atom_bytes: text.len(),
        ..IrParseBudget::default()
    };
    match parse_ir_with_budget(&text, budget) {
        Ok(back) if back.scenario() == s => {}
        Ok(_) => out.push(Violation {
            code: "ir-round-trip-drift",
            what: format!("scenario {:?} reparses to a different value", s.name),
            fix: "report this as an fs-scenario IR bug (canonical form must be lossless)"
                .to_string(),
        }),
        Err(e) => out.push(Violation {
            code: "ir-round-trip-parse",
            what: format!("scenario {:?} canonical IR failed to reparse: {e}", s.name),
            fix: "report this as an fs-scenario IR bug".to_string(),
        }),
    }
}

#[cfg(test)]
mod allocation_internal_tests {
    use super::{reserve_decoded, reserve_decoded_string};
    use crate::{IrSourceSpan, ScenarioError};

    #[test]
    fn decoded_collection_allocation_refusal_is_typed() {
        let mut values = Vec::<f64>::new();
        let error = reserve_decoded(
            &mut values,
            usize::MAX,
            "collection test",
            IrSourceSpan { start: 0, end: 0 },
        )
        .expect_err("impossible decoded capacity must be refused");

        assert!(matches!(
            error,
            ScenarioError::Parse {
                span: IrSourceSpan { start: 0, end: 0 },
                path,
                what,
            }
                if path.is_empty()
                    && what.contains("IR decoded collection test allocation")
                    && what.contains("was refused")
        ));
        assert!(values.is_empty());

        let mut value = String::new();
        let error =
            reserve_decoded_string(&mut value, usize::MAX, IrSourceSpan { start: 0, end: 0 })
                .expect_err("impossible decoded string capacity must be refused");
        assert!(matches!(
            error,
            ScenarioError::Parse {
                span: IrSourceSpan { start: 0, end: 0 },
                path,
                what,
            }
                if path.is_empty()
                    && what.contains("IR decoded string allocation")
                    && what.contains("was refused")
        ));
        assert!(value.is_empty());
    }
}
