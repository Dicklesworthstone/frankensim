//! Canonical scenario IR: a deterministic s-expression encoding with
//! LOSSLESS round-trip (floats print in shortest-round-trip form; dims
//! travel as explicit SI exponent vectors, so no unit-string parsing is
//! involved). `write_ir` output is byte-stable — the ledger stores it as
//! a content-addressed artifact; `parse_ir` inverts it exactly.

use crate::ScenarioError;
use crate::bc::{BcKind, BcValue, BoundaryCondition, Compat, Physics};
use crate::ensemble::{SpectrumModel, StochasticEnsemble};
use crate::frame::{Frame, FrameId, FrameMotion, FrameTree};
use crate::scenario::{
    Combination, ContactLaw, ContactModel, Environment, LoadCase, Scenario, Violation,
};
use crate::signal::{ChebProfile, Interp, TimeSignal};
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

/// A decoded scenario together with its wire-version migration context.
#[derive(Debug, Clone, PartialEq)]
pub struct DecodedScenario {
    /// Decoded scenario value.
    scenario: Scenario,
    /// Version found on the wire; historical unversioned forms are v1.
    source_version: u32,
    /// Present exactly when a legacy five-base form was crossed into six-base memory.
    dimension_crosswalk: Option<DimensionCrosswalkReceipt>,
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
enum Sx {
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
    max_atom_bytes: 1024 * 1024,
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
        at,
        what: what.to_string(),
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
            return Err(err(
                budget.max_bytes.min(text.len()),
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
            return Err(err(self.pos, "trailing bytes after the scenario form"));
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
                            return Ok(Sx::List(items));
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
            Some(b')') => Err(err(self.pos, "unexpected ')'")),
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
                            let c = *self
                                .bytes
                                .get(self.pos + 1)
                                .ok_or_else(|| err(self.pos, "bad escape"))?;
                            if !matches!(c, b'"' | b'\\') {
                                return Err(err(
                                    self.pos,
                                    "unsupported string escape; only quote and backslash escapes are canonical",
                                ));
                            }
                            if buf.len() >= self.budget.max_atom_bytes {
                                return Err(err(self.pos, "IR string exceeds atom-byte budget"));
                            }
                            self.reserve(&mut buf, 1)?;
                            buf.push(c);
                            self.pos += 2;
                        }
                        Some(b'"') => {
                            self.pos += 1;
                            let s = String::from_utf8(buf)
                                .map_err(|_| err(self.pos, "string is not valid UTF-8"))?;
                            return Ok(Sx::Str(s));
                        }
                        Some(&c) => {
                            if buf.len() >= self.budget.max_atom_bytes {
                                return Err(err(self.pos, "IR string exceeds atom-byte budget"));
                            }
                            self.reserve(&mut buf, 1)?;
                            buf.push(c);
                            self.pos += 1;
                        }
                    }
                }
            }
            Some(_) => {
                let start = self.pos;
                while self.pos < self.bytes.len()
                    && !self.bytes[self.pos].is_ascii_whitespace()
                    && self.bytes[self.pos] != b'('
                    && self.bytes[self.pos] != b')'
                {
                    self.pos += 1;
                    if self.pos - start > self.budget.max_atom_bytes {
                        return Err(err(start, "IR atom exceeds atom-byte budget"));
                    }
                }
                let atom = std::str::from_utf8(&self.bytes[start..self.pos])
                    .map_err(|_| err(start, "atom is not valid UTF-8"))?;
                let mut owned = String::new();
                owned
                    .try_reserve_exact(atom.len())
                    .map_err(|allocation_error| {
                        err(
                            start,
                            &format!("IR atom allocation refused: {allocation_error}"),
                        )
                    })?;
                owned.push_str(atom);
                Ok(Sx::Atom(owned))
            }
        }
    }
}

fn parse_sx(text: &str, budget: IrParseBudget) -> Result<Sx, ScenarioError> {
    SxParser::new(text, budget)?.parse()
}

// -------------------------------------------------------------- decoding

fn as_list<'a>(sx: &'a Sx, head: &str) -> Result<&'a [Sx], ScenarioError> {
    if let Sx::List(items) = sx
        && let Some(Sx::Atom(a)) = items.first()
        && a == head
    {
        return Ok(&items[1..]);
    }
    Err(err(0, &format!("expected ({head} ...) form")))
}

fn as_f64(sx: &Sx) -> Result<f64, ScenarioError> {
    if let Sx::Atom(a) = sx
        && let Ok(v) = a.parse::<f64>()
        && v.is_finite()
    {
        return Ok(v);
    }
    Err(err(0, "expected a finite number"))
}

fn as_str(sx: &Sx) -> Result<String, ScenarioError> {
    if let Sx::Str(s) = sx {
        return Ok(s.clone());
    }
    Err(err(0, "expected a string"))
}

fn as_atom(sx: &Sx) -> Result<&str, ScenarioError> {
    if let Sx::Atom(a) = sx {
        return Ok(a);
    }
    Err(err(0, "expected an atom"))
}

fn as_u64(sx: &Sx) -> Result<u64, ScenarioError> {
    if let Sx::Atom(a) = sx
        && let Ok(v) = a.parse::<u64>()
    {
        return Ok(v);
    }
    Err(err(0, "expected an unsigned integer"))
}

fn as_u32(sx: &Sx) -> Result<u32, ScenarioError> {
    if let Sx::Atom(a) = sx
        && let Ok(v) = a.parse::<u32>()
    {
        return Ok(v);
    }
    Err(err(0, "expected a u32"))
}

fn as_i8(sx: &Sx) -> Result<i8, ScenarioError> {
    if let Sx::Atom(a) = sx
        && let Ok(v) = a.parse::<i8>()
    {
        return Ok(v);
    }
    Err(err(0, "expected a small integer exponent"))
}

fn as_dims_at(items: &[Sx], wire: DimensionWire) -> Result<Dims, ScenarioError> {
    let expected = match wire {
        DimensionWire::LegacyFive => 5,
        DimensionWire::CanonicalSix => 6,
    };
    if items.len() != expected {
        return Err(err(
            0,
            &format!("scenario IR needs exactly {expected} dimension exponents"),
        ));
    }
    let mut d = [0i8; 6];
    for (slot, sx) in d.iter_mut().zip(items) {
        *slot = as_i8(sx)?;
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
        return Err(err(
            0,
            &format!("qty needs a value plus {} exponents", expected - 1),
        ));
    }
    Ok(QtyAny::new(
        as_f64(&items[0])?,
        as_dims_at(&items[1..], wire)?,
    ))
}

fn as_dims(sx: &Sx, wire: DimensionWire) -> Result<Dims, ScenarioError> {
    as_dims_at(as_list(sx, "dims")?, wire)
}

fn as_floats(sx: &Sx, head: &str) -> Result<Vec<f64>, ScenarioError> {
    as_list(sx, head)?.iter().map(as_f64).collect()
}

fn as_vec3(sx: &Sx) -> Result<Vec3, ScenarioError> {
    let items = as_list(sx, "vec")?;
    if items.len() != 3 {
        return Err(err(0, "vec needs three components"));
    }
    Ok(Vec3::new(
        as_f64(&items[0])?,
        as_f64(&items[1])?,
        as_f64(&items[2])?,
    ))
}

fn as_profile(items: &[Sx], wire: DimensionWire) -> Result<ChebProfile, ScenarioError> {
    if items.len() != 4 {
        return Err(err(0, "profile needs dims, domain, coeffs"));
    }
    let dims = as_dims(&items[0], wire)?;
    let a = as_f64(&items[1])?;
    let b = as_f64(&items[2])?;
    let coeffs = as_floats(&items[3], "coeffs")?;
    if !(a < b) {
        return Err(err(0, "profile domain must satisfy finite a < b"));
    }
    if coeffs.is_empty() {
        return Err(err(0, "profile needs at least one finite coefficient"));
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
    let Sx::List(items) = sx else {
        return Err(err(0, "expected a signal form"));
    };
    let head = as_atom(items.first().ok_or_else(|| err(0, "empty signal"))?)?;
    let rest = &items[1..];
    match head {
        "constant" => {
            if rest.len() != 1 {
                return Err(err(0, "constant needs a value"));
            }
            Ok(TimeSignal::Constant(as_qty(&rest[0], wire)?))
        }
        "ramp" => {
            if rest.len() != 4 {
                return Err(err(0, "ramp needs t0 t1 from to"));
            }
            Ok(TimeSignal::Ramp {
                t_start: as_f64(&rest[0])?,
                t_end: as_f64(&rest[1])?,
                from: as_qty(&rest[2], wire)?,
                to: as_qty(&rest[3], wire)?,
            })
        }
        "table" => {
            if rest.len() != 4 {
                return Err(err(0, "table needs interp dims times values"));
            }
            let interp = match as_atom(&rest[0])? {
                "linear" => Interp::Linear,
                "hold" => Interp::Hold,
                other => return Err(err(0, &format!("unknown interp {other:?}"))),
            };
            Ok(TimeSignal::Table {
                interp,
                dims: as_dims(&rest[1], wire)?,
                times: as_floats(&rest[2], "times")?,
                values: as_floats(&rest[3], "values")?,
            })
        }
        "chebfun" => Ok(TimeSignal::Chebfun(as_profile(rest, wire)?)),
        other => Err(err(0, &format!("unknown signal {other:?}"))),
    }
}

fn as_physics(a: &str) -> Result<Physics, ScenarioError> {
    match a {
        "incompressible-flow" => Ok(Physics::IncompressibleFlow),
        "thermal" => Ok(Physics::Thermal),
        "elasticity" => Ok(Physics::Elasticity),
        other => Err(err(0, &format!("unknown physics {other:?}"))),
    }
}

fn as_kind(a: &str) -> Result<BcKind, ScenarioError> {
    match a {
        "dirichlet" => Ok(BcKind::Dirichlet),
        "neumann" => Ok(BcKind::Neumann),
        "robin" => Ok(BcKind::Robin),
        "mass-flow-inlet" => Ok(BcKind::MassFlowInlet),
        "pressure-outlet" => Ok(BcKind::PressureOutlet),
        "wall-no-slip" => Ok(BcKind::WallNoSlip),
        "wall-slip" => Ok(BcKind::WallSlip),
        "traction" => Ok(BcKind::Traction),
        other => Err(err(0, &format!("unknown bc kind {other:?}"))),
    }
}

fn as_bc(sx: &Sx, wire: DimensionWire) -> Result<BoundaryCondition, ScenarioError> {
    let items = as_list(sx, "bc")?;
    if items.len() != 6 {
        return Err(err(0, "bc needs region physics kind frame value compat"));
    }
    let value = match &items[4] {
        Sx::Atom(a) if a == "none" => None,
        Sx::List(inner) => {
            let head = as_atom(inner.first().ok_or_else(|| err(0, "empty bc value"))?)?;
            match head {
                "uniform" => {
                    if inner.len() != 2 {
                        return Err(err(0, "uniform bc value needs a quantity"));
                    }
                    Some(BcValue::Uniform(as_qty(&inner[1], wire)?))
                }
                "signal" => {
                    if inner.len() != 2 {
                        return Err(err(0, "signal bc value needs a signal form"));
                    }
                    Some(BcValue::Signal(as_signal(&inner[1], wire)?))
                }
                "profile" => Some(BcValue::Profile(as_profile(&inner[1..], wire)?)),
                other => return Err(err(0, &format!("unknown bc value {other:?}"))),
            }
        }
        _ => return Err(err(0, "bad bc value")),
    };
    let compatibility = match as_atom(&items[5])? {
        "none" => None,
        "incompressible" => Some(Compat::Incompressible),
        other => return Err(err(0, &format!("unknown compat {other:?}"))),
    };
    Ok(BoundaryCondition {
        region: as_str(&items[0])?,
        physics: as_physics(as_atom(&items[1])?)?,
        kind: as_kind(as_atom(&items[2])?)?,
        frame: as_u32(&items[3])?,
        value,
        compatibility,
    })
}

fn as_frame(sx: &Sx, wire: DimensionWire) -> Result<Frame, ScenarioError> {
    let items = as_list(sx, "frame")?;
    if items.len() != 4 {
        return Err(err(0, "frame needs id name parent motion"));
    }
    let Sx::List(motion_items) = &items[3] else {
        return Err(err(0, "bad frame motion"));
    };
    let head = as_atom(motion_items.first().ok_or_else(|| err(0, "empty motion"))?)?;
    let rest = &motion_items[1..];
    let motion = match head {
        "fixed" => {
            if rest.len() != 2 {
                return Err(err(0, "fixed motion needs orientation and translation"));
            }
            let q_items = as_list(&rest[0], "quat")?;
            if q_items.len() != 4 {
                return Err(err(0, "quat needs four components"));
            }
            FrameMotion::Fixed {
                orientation: Quat {
                    w: as_f64(&q_items[0])?,
                    x: as_f64(&q_items[1])?,
                    y: as_f64(&q_items[2])?,
                    z: as_f64(&q_items[3])?,
                },
                translation: as_vec3(&rest[1])?,
            }
        }
        "rotating" => {
            if rest.len() != 3 {
                return Err(err(0, "rotating motion needs axis center rate"));
            }
            let axis = as_vec3(&rest[0])?;
            FrameMotion::Rotating {
                axis: [axis.x, axis.y, axis.z],
                center: as_vec3(&rest[1])?,
                rate: as_qty(&rest[2], wire)?,
            }
        }
        "tilt" => {
            if rest.len() != 3 {
                return Err(err(0, "tilt motion needs axis center angle"));
            }
            let axis = as_vec3(&rest[0])?;
            FrameMotion::Tilt {
                axis: [axis.x, axis.y, axis.z],
                center: as_vec3(&rest[1])?,
                angle: as_signal(&rest[2], wire)?,
            }
        }
        other => return Err(err(0, &format!("unknown motion {other:?}"))),
    };
    Ok(Frame {
        id: FrameId(as_u32(&items[0])?),
        name: as_str(&items[1])?,
        parent: FrameId(as_u32(&items[2])?),
        motion,
    })
}

fn as_model(sx: &Sx, wire: DimensionWire) -> Result<SpectrumModel, ScenarioError> {
    let Sx::List(items) = sx else {
        return Err(err(0, "expected a spectrum model form"));
    };
    let head = as_atom(items.first().ok_or_else(|| err(0, "empty model"))?)?;
    let rest = &items[1..];
    match head {
        "dryden" => {
            if rest.len() != 3 {
                return Err(err(0, "dryden needs sigma length_scale mean_speed"));
            }
            Ok(SpectrumModel::Dryden {
                sigma: as_qty(&rest[0], wire)?,
                length_scale: as_qty(&rest[1], wire)?,
                mean_speed: as_qty(&rest[2], wire)?,
            })
        }
        "kanai-tajimi" => {
            if rest.len() != 3 {
                return Err(err(0, "kanai-tajimi needs s0 omega_g zeta_g"));
            }
            Ok(SpectrumModel::KanaiTajimi {
                s0: as_f64(&rest[0])?,
                omega_g: as_qty(&rest[1], wire)?,
                zeta_g: as_f64(&rest[2])?,
            })
        }
        "carreau" => {
            if rest.len() != 8 {
                return Err(err(0, "carreau needs six qty bounds + two n bounds"));
            }
            Ok(SpectrumModel::CarreauBand {
                eta_zero: [as_qty(&rest[0], wire)?, as_qty(&rest[1], wire)?],
                eta_inf: [as_qty(&rest[2], wire)?, as_qty(&rest[3], wire)?],
                lambda: [as_qty(&rest[4], wire)?, as_qty(&rest[5], wire)?],
                n: [as_f64(&rest[6])?, as_f64(&rest[7])?],
            })
        }
        other => Err(err(0, &format!("unknown model {other:?}"))),
    }
}

fn as_ensemble(sx: &Sx, wire: DimensionWire) -> Result<StochasticEnsemble, ScenarioError> {
    let items = as_list(sx, "ensemble")?;
    if items.len() != 6 {
        return Err(err(0, "ensemble needs name seed members duration dt model"));
    }
    Ok(StochasticEnsemble {
        name: as_str(&items[0])?,
        seed: as_u64(&items[1])?,
        members: as_u32(&items[2])?,
        duration: as_qty(&items[3], wire)?,
        dt: as_qty(&items[4], wire)?,
        model: as_model(&items[5], wire)?,
    })
}

fn as_contact(sx: &Sx) -> Result<ContactLaw, ScenarioError> {
    let contact_items = as_list(sx, "contact")?;
    if contact_items.len() != 3 {
        return Err(err(0, "contact needs two regions + model"));
    }
    let model = match &contact_items[2] {
        Sx::Atom(a) if a == "frictionless" => ContactModel::Frictionless,
        Sx::Atom(a) if a == "tied" => ContactModel::Tied,
        other @ Sx::List(_) => {
            let coulomb = as_list(other, "coulomb")?;
            if coulomb.len() != 2 {
                return Err(err(0, "coulomb needs mu_s + mu_k"));
            }
            ContactModel::Coulomb {
                mu_s: as_f64(&coulomb[0])?,
                mu_k: as_f64(&coulomb[1])?,
            }
        }
        _ => return Err(err(0, "bad contact model")),
    };
    Ok(ContactLaw {
        region_a: as_str(&contact_items[0])?,
        region_b: as_str(&contact_items[1])?,
        model,
    })
}

fn as_environment(sx: &Sx, wire: DimensionWire) -> Result<Environment, ScenarioError> {
    let items = as_list(sx, "environment")?;
    if items.len() != 5 {
        return Err(err(
            0,
            "environment needs gravity x3 + temperature + pressure",
        ));
    }
    Ok(Environment {
        gravity: [
            as_qty(&items[0], wire)?,
            as_qty(&items[1], wire)?,
            as_qty(&items[2], wire)?,
        ],
        ambient_temperature: as_qty(&items[3], wire)?,
        ambient_pressure: as_qty(&items[4], wire)?,
    })
}

fn scenario_wire_header(root_items: &[Sx]) -> Result<(u32, DimensionWire, &[Sx]), ScenarioError> {
    if matches!(root_items.first(), Some(Sx::Atom(key)) if key == ":version") {
        let version = root_items
            .get(1)
            .ok_or_else(|| err(0, "scenario :version needs a value"))
            .and_then(as_u32)?;
        let wire = match version {
            LEGACY_SCENARIO_IR_VERSION => DimensionWire::LegacyFive,
            SCENARIO_IR_VERSION => DimensionWire::CanonicalSix,
            _ => {
                return Err(err(
                    0,
                    &format!(
                        "unsupported scenario IR version {version}; supported versions are 1 and {SCENARIO_IR_VERSION}"
                    ),
                ));
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
/// The decoded value retains its wire version and any dimension crosswalk
/// receipt. Canonical v2 uses six exponents; explicit v1 and the historical
/// unversioned form use five and append `mol = 0`.
///
/// # Errors
/// [`ScenarioError::Parse`] for malformed, non-finite, or over-budget input.
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
/// [`ScenarioError::Parse`] for malformed, non-finite, or over-budget input.
pub fn parse_ir_with_budget(
    text: &str,
    budget: IrParseBudget,
) -> Result<DecodedScenario, ScenarioError> {
    let root = parse_sx(text, budget)?;
    let root_items = as_list(&root, "scenario")?;
    let (source_version, wire, items) = scenario_wire_header(root_items)?;
    if items.len() != 9 {
        return Err(err(0, "scenario needs name seed + seven sections"));
    }
    let environment = as_environment(&items[2], wire)?;
    let mut frames = FrameTree::new();
    for f in as_list(&items[3], "frames")? {
        frames.add(as_frame(f, wire)?);
    }
    let base_bcs = as_list(&items[4], "bcs")?
        .iter()
        .map(|bc| as_bc(bc, wire))
        .collect::<Result<Vec<_>, _>>()?;
    let mut cases = Vec::new();
    for c in as_list(&items[5], "cases")? {
        let case_items = as_list(c, "case")?;
        let name = as_str(
            case_items
                .first()
                .ok_or_else(|| err(0, "case needs a name"))?,
        )?;
        let bcs = case_items[1..]
            .iter()
            .map(|bc| as_bc(bc, wire))
            .collect::<Result<Vec<_>, _>>()?;
        cases.push(LoadCase { name, bcs });
    }
    let mut combinations = Vec::new();
    for c in as_list(&items[6], "combos")? {
        let combo_items = as_list(c, "combo")?;
        let name = as_str(
            combo_items
                .first()
                .ok_or_else(|| err(0, "combo needs a name"))?,
        )?;
        let mut terms = Vec::new();
        for t in &combo_items[1..] {
            let term_items = as_list(t, "term")?;
            if term_items.len() != 2 {
                return Err(err(0, "term needs case + factor"));
            }
            terms.push((as_str(&term_items[0])?, as_f64(&term_items[1])?));
        }
        combinations.push(Combination { name, terms });
    }
    let ensembles = as_list(&items[7], "ensembles")?
        .iter()
        .map(|ensemble| as_ensemble(ensemble, wire))
        .collect::<Result<Vec<_>, _>>()?;
    let contacts = as_list(&items[8], "contacts")?
        .iter()
        .map(as_contact)
        .collect::<Result<Vec<_>, _>>()?;
    let scenario = Scenario {
        name: as_str(&items[0])?,
        seed: as_u64(&items[1])?,
        frames,
        base_bcs,
        cases,
        combinations,
        ensembles,
        contacts,
        environment,
    };
    let dimension_crosswalk = if wire == DimensionWire::LegacyFive {
        let canonical = write_ir(&scenario);
        Some(DimensionCrosswalkReceipt {
            source_version,
            target_version: SCENARIO_IR_VERSION,
            old_hash: hash_bytes(text.as_bytes()),
            new_hash: hash_bytes(canonical.as_bytes()),
            source_width: 5,
            target_width: 6,
            rule: FiveToSixRule::AppendMoleZero,
        })
    } else {
        None
    };
    Ok(DecodedScenario {
        scenario,
        source_version,
        dimension_crosswalk,
    })
}

/// Round-trip helper for lints/tests: violations if reparse ≠ original.
pub fn check_round_trip(s: &Scenario, out: &mut Vec<Violation>) {
    let text = write_ir(s);
    match parse_ir(&text) {
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
