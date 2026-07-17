//! Local execution-containment model (i94v.7.3.2).
//!
//! Within one Attempt, every execution event has exactly ONE primary
//! containment parent, forming a deterministic tree rooted at the propagated
//! [`AttemptId`]. DSR-hosting, campaign-selection, and shard-membership are
//! typed CONTEXTUAL edges carried alongside — structurally separate from
//! parentage so a convenient local tree can never alias or flatten the real
//! multi-parent durable provenance owned by V.3.8 (submissions, jobs,
//! attempts, checkpoints, artifacts, forks, replays, idempotency).
//!
//! The tree is honest about what it does not know: reordered delivery is
//! buffered, duplicate delivery is idempotent, and a parent that never
//! arrives seals as an explicit [`IntegrityGap`] rather than silently
//! re-rooting the orphan.

use core::fmt;

/// Logical containment-model version. v1 projects through the typed
/// [`crate::EventKind::ContainmentNode`] / [`crate::EventKind::ContainmentGap`]
/// kinds (identity-bound per field); bump on any containment-semantics
/// change.
pub const CONTAINMENT_WIRE_VERSION: u32 = 1;

fn valid_id_chars(s: &str) -> bool {
    !s.is_empty()
        && s.len() <= 256
        && s.chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '-' | '_' | '.' | ':'))
}

/// Refusal for one malformed identity component.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IdError {
    /// Which typed identity refused.
    pub role: &'static str,
    /// The offending raw text.
    pub raw: String,
}

impl fmt::Display for IdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} must be 1..=256 chars of [A-Za-z0-9._:-]; got {:?}",
            self.role, self.raw
        )
    }
}

impl core::error::Error for IdError {}

/// One propagated Attempt identity. OPAQUE HERE: V.3.8 owns attempt
/// semantics (durable submission, retries, idempotency); this module only
/// roots the local tree at whatever token was propagated to the process.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AttemptId(String);

/// One DSR invocation.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct DsrRunId(String);

/// One canonical manifest execution.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CampaignRunId(String);

/// One deterministic partition of a campaign.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ShardId(String);

/// One logical journey definition (not an execution).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct JourneyId(String);

/// One logical case definition (not an execution).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CaseId(String);

/// One user-meaningful operation inside an Attempt.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ExecutionOpId(String);

/// One execution scope inside an operation (mirrors the emitter scope tree).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ExecutionScopeId(String);

/// One tile of work inside a scope (the same logical tile identity that
/// keys RNG streams and deterministic reductions).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct TileId(String);

macro_rules! id_impl {
    ($ty:ident, $role:literal) => {
        impl $ty {
            /// Validate and adopt one raw identity token.
            ///
            /// # Errors
            /// [`IdError`] when the token is empty, oversized, or carries
            /// characters outside `[A-Za-z0-9._:-]`.
            pub fn new(raw: impl Into<String>) -> Result<Self, IdError> {
                let raw = raw.into();
                if valid_id_chars(&raw) {
                    Ok(Self(raw))
                } else {
                    Err(IdError { role: $role, raw })
                }
            }

            /// The validated raw token.
            #[must_use]
            pub fn as_str(&self) -> &str {
                &self.0
            }
        }
    };
}

id_impl!(AttemptId, "AttemptId");
id_impl!(DsrRunId, "DsrRunId");
id_impl!(CampaignRunId, "CampaignRunId");
id_impl!(ShardId, "ShardId");
id_impl!(JourneyId, "JourneyId");
id_impl!(CaseId, "CaseId");
id_impl!(ExecutionOpId, "ExecutionOpId");
id_impl!(ExecutionScopeId, "ExecutionScopeId");
id_impl!(TileId, "TileId");

/// Identity of one local tree node. The role tag is part of identity:
/// an op and a tile sharing raw text are still different nodes, so ID text
/// can never confuse roles.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum LocalNodeId {
    /// A user-meaningful operation.
    Op(ExecutionOpId),
    /// An execution scope.
    Scope(ExecutionScopeId),
    /// A tile of work.
    Tile(TileId),
}

impl LocalNodeId {
    /// Stable role name for wire projection and diagnostics.
    #[must_use]
    pub fn role(&self) -> &'static str {
        match self {
            Self::Op(_) => "op",
            Self::Scope(_) => "scope",
            Self::Tile(_) => "tile",
        }
    }

    /// The raw identity token inside the role.
    #[must_use]
    pub fn raw(&self) -> &str {
        match self {
            Self::Op(id) => id.as_str(),
            Self::Scope(id) => id.as_str(),
            Self::Tile(id) => id.as_str(),
        }
    }
}

/// The one primary containment parent of a local node. Exactly one; the
/// root is always the propagated Attempt.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum LocalParent {
    /// Directly contained by the Attempt root.
    AttemptRoot,
    /// Contained by another local node.
    Node(LocalNodeId),
}

/// Typed contextual edges for one node: WHERE the work sits in the
/// campaign geometry. Deliberately not a parent: hosting/selection/
/// membership never alias containment.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct ContainmentContext {
    /// Hosting DSR invocation, when hosted by a DSR.
    pub dsr_run: Option<DsrRunId>,
    /// Selecting campaign execution, when campaign-selected.
    pub campaign_run: Option<CampaignRunId>,
    /// Deterministic shard membership.
    pub shard: Option<ShardId>,
    /// Logical journey definition this work realizes.
    pub journey: Option<JourneyId>,
    /// Logical case definition this work realizes.
    pub case: Option<CaseId>,
}

/// One local execution-containment record: a node, its single primary
/// parent, its deterministic sibling sequence, and its typed context.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContainmentRecord {
    /// Node identity.
    pub node: LocalNodeId,
    /// The one primary containment parent.
    pub parent: LocalParent,
    /// Deterministic sequence among siblings (logical order, never wall
    /// clock).
    pub seq: u64,
    /// Typed contextual edges.
    pub context: ContainmentContext,
}

impl ContainmentRecord {
    fn fingerprint(&self) -> u64 {
        let mut bytes = Vec::with_capacity(96);
        let mut push = |tag: &str, val: &str| {
            bytes.extend_from_slice(tag.as_bytes());
            bytes.extend_from_slice(&(val.len() as u64).to_le_bytes());
            bytes.extend_from_slice(val.as_bytes());
        };
        push("role", self.node.role());
        push("node", self.node.raw());
        match &self.parent {
            LocalParent::AttemptRoot => push("parent", "\u{1}attempt-root"),
            LocalParent::Node(id) => {
                push("parent-role", id.role());
                push("parent", id.raw());
            }
        }
        push("seq", &self.seq.to_string());
        let ctx = &self.context;
        push("dsr", ctx.dsr_run.as_ref().map_or("", |v| v.as_str()));
        push(
            "campaign",
            ctx.campaign_run.as_ref().map_or("", |v| v.as_str()),
        );
        push("shard", ctx.shard.as_ref().map_or("", |v| v.as_str()));
        push("journey", ctx.journey.as_ref().map_or("", |v| v.as_str()));
        push("case", ctx.case.as_ref().map_or("", |v| v.as_str()));
        crate::fnv1a64(&bytes)
    }
}

/// Fail-closed refusals for containment ingest and seal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ContainmentError {
    /// The same node identity arrived again with DIFFERENT content: this is
    /// aliasing or corruption, never idempotent redelivery.
    ConflictingRedelivery {
        /// Node role.
        role: &'static str,
        /// Node raw identity.
        node: String,
    },
    /// A node named itself as its own parent.
    SelfParent {
        /// Node role.
        role: &'static str,
        /// Node raw identity.
        node: String,
    },
    /// Admitting this record would close a containment cycle.
    Cycle {
        /// Node role at which the cycle was detected.
        role: &'static str,
        /// Node raw identity.
        node: String,
    },
    /// The tree was asked to embed under a different Attempt root than the
    /// one it was propagated.
    WrongAttemptRoot {
        /// Root this tree was built under.
        expected: String,
        /// Root the caller supplied.
        found: String,
    },
}

impl fmt::Display for ContainmentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ConflictingRedelivery { role, node } => write!(
                f,
                "{role} node {node:?} redelivered with different content; \
                 refusing aliased identity"
            ),
            Self::SelfParent { role, node } => {
                write!(f, "{role} node {node:?} cannot contain itself")
            }
            Self::Cycle { role, node } => write!(
                f,
                "admitting {role} node {node:?} would close a containment cycle"
            ),
            Self::WrongAttemptRoot { expected, found } => write!(
                f,
                "local tree is rooted at attempt {expected:?}; refusing embedding \
                 under attempt {found:?}"
            ),
        }
    }
}

impl core::error::Error for ContainmentError {}

/// One explicit integrity gap in a sealed tree: evidence that lineage is
/// incomplete, never silently repaired.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntegrityGap {
    /// Node whose lineage is incomplete.
    pub node_role: &'static str,
    /// Raw identity of that node.
    pub node: String,
    /// Role of the parent that never arrived.
    pub missing_parent_role: &'static str,
    /// Raw identity of the parent that never arrived.
    pub missing_parent: String,
}

/// A deterministic local containment tree under one propagated Attempt.
#[derive(Debug)]
pub struct AttemptTree {
    root: AttemptId,
    // Admitted nodes: id -> (record, fingerprint).
    admitted: Vec<(ContainmentRecord, u64)>,
    // Records whose parent has not arrived yet.
    pending: Vec<(ContainmentRecord, u64)>,
}

/// Outcome of one ingest call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ingest {
    /// Newly admitted into the tree.
    Admitted,
    /// Byte-identical redelivery; nothing changed.
    Duplicate,
    /// Parent not yet known; buffered until it arrives.
    Buffered,
}

impl AttemptTree {
    /// Start an empty tree rooted at the propagated Attempt.
    #[must_use]
    pub fn new(root: AttemptId) -> Self {
        Self {
            root,
            admitted: Vec::new(),
            pending: Vec::new(),
        }
    }

    /// The propagated Attempt this tree is rooted at.
    #[must_use]
    pub fn root(&self) -> &AttemptId {
        &self.root
    }

    fn find(&self, node: &LocalNodeId) -> Option<&(ContainmentRecord, u64)> {
        self.admitted
            .iter()
            .chain(self.pending.iter())
            .find(|(r, _)| &r.node == node)
    }

    fn parent_known(&self, parent: &LocalParent) -> bool {
        match parent {
            LocalParent::AttemptRoot => true,
            LocalParent::Node(id) => self.admitted.iter().any(|(r, _)| &r.node == id),
        }
    }

    fn would_cycle(&self, record: &ContainmentRecord) -> bool {
        // Walk ancestors of the proposed parent through admitted AND pending
        // records; if the new node appears, admitting it closes a cycle.
        let mut cursor = record.parent.clone();
        let mut steps = 0usize;
        while let LocalParent::Node(id) = cursor {
            if id == record.node {
                return true;
            }
            steps += 1;
            if steps > self.admitted.len() + self.pending.len() + 1 {
                return true; // defensive: malformed chains terminate
            }
            match self.find(&id) {
                Some((r, _)) => cursor = r.parent.clone(),
                None => return false,
            }
        }
        false
    }

    /// Ingest one record: idempotent for byte-identical redelivery, buffered
    /// for reordered delivery, refused for aliasing, self-parentage, and
    /// cycles.
    ///
    /// # Errors
    /// [`ContainmentError`] on conflicting redelivery, self-parentage, or a
    /// containment cycle.
    pub fn ingest(&mut self, record: ContainmentRecord) -> Result<Ingest, ContainmentError> {
        if let LocalParent::Node(parent_id) = &record.parent
            && parent_id == &record.node
        {
            return Err(ContainmentError::SelfParent {
                role: record.node.role(),
                node: record.node.raw().to_string(),
            });
        }
        let fp = record.fingerprint();
        if let Some((_, existing_fp)) = self.find(&record.node) {
            return if *existing_fp == fp {
                Ok(Ingest::Duplicate)
            } else {
                Err(ContainmentError::ConflictingRedelivery {
                    role: record.node.role(),
                    node: record.node.raw().to_string(),
                })
            };
        }
        if self.would_cycle(&record) {
            return Err(ContainmentError::Cycle {
                role: record.node.role(),
                node: record.node.raw().to_string(),
            });
        }
        if self.parent_known(&record.parent) {
            self.admitted.push((record, fp));
            self.drain_pending();
            Ok(Ingest::Admitted)
        } else {
            self.pending.push((record, fp));
            Ok(Ingest::Buffered)
        }
    }

    fn drain_pending(&mut self) {
        loop {
            let Some(pos) = self
                .pending
                .iter()
                .position(|(r, _)| self.parent_known(&r.parent))
            else {
                return;
            };
            let entry = self.pending.remove(pos);
            self.admitted.push(entry);
        }
    }

    /// Seal the tree: deterministic parent-major, seq-minor order over the
    /// admitted nodes, plus an explicit gap ledger for every record whose
    /// parent never arrived. Gaps are evidence, not errors: closure is
    /// incomplete and downstream adjudication must see that.
    #[must_use]
    pub fn seal(mut self) -> SealedAttemptTree {
        self.admitted.sort_by(|(a, _), (b, _)| {
            let key = |r: &ContainmentRecord| {
                (
                    match &r.parent {
                        LocalParent::AttemptRoot => (0u8, String::new(), String::new()),
                        LocalParent::Node(id) => (1u8, id.role().to_string(), id.raw().to_string()),
                    },
                    r.seq,
                    r.node.role().to_string(),
                    r.node.raw().to_string(),
                )
            };
            key(a).cmp(&key(b))
        });
        let gaps = self
            .pending
            .iter()
            .map(|(r, _)| {
                let (role, parent) = match &r.parent {
                    LocalParent::Node(id) => (id.role(), id.raw().to_string()),
                    LocalParent::AttemptRoot => unreachable!("root parents are always known"),
                };
                IntegrityGap {
                    node_role: r.node.role(),
                    node: r.node.raw().to_string(),
                    missing_parent_role: role,
                    missing_parent: parent,
                }
            })
            .collect();
        SealedAttemptTree {
            root: self.root,
            nodes: self.admitted.into_iter().map(|(r, _)| r).collect(),
            gaps,
        }
    }
}

/// An immutable, deterministically ordered local tree plus its explicit
/// gap ledger — the exact shape handed to the V.3.8 embedder.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SealedAttemptTree {
    root: AttemptId,
    nodes: Vec<ContainmentRecord>,
    gaps: Vec<IntegrityGap>,
}

impl SealedAttemptTree {
    /// The propagated Attempt root.
    #[must_use]
    pub fn root(&self) -> &AttemptId {
        &self.root
    }

    /// Admitted nodes in deterministic parent-major, seq-minor order.
    #[must_use]
    pub fn nodes(&self) -> &[ContainmentRecord] {
        &self.nodes
    }

    /// Explicit lineage gaps. Non-empty gaps mean closure is incomplete.
    #[must_use]
    pub fn gaps(&self) -> &[IntegrityGap] {
        &self.gaps
    }

    /// Project the sealed tree into fs-obs JSONL events under the typed
    /// [`crate::EventKind::ContainmentNode`] /
    /// [`crate::EventKind::ContainmentGap`] kinds.
    #[must_use]
    pub fn to_events(&self, session: &str, scope: &str) -> Vec<crate::Event> {
        let mut em = crate::Emitter::new(session, scope);
        let mut events = Vec::with_capacity(self.nodes.len() + self.gaps.len());
        let opt = |v: Option<&str>| v.unwrap_or("").to_string();
        for record in &self.nodes {
            let (parent_role, parent) = match &record.parent {
                LocalParent::AttemptRoot => ("attempt", self.root.as_str()),
                LocalParent::Node(id) => (id.role(), id.raw()),
            };
            let ctx = &record.context;
            events.push(em.emit(
                crate::Severity::Trace,
                crate::EventKind::ContainmentNode {
                    attempt: self.root.as_str().to_string(),
                    role: record.node.role().to_string(),
                    node: record.node.raw().to_string(),
                    parent_role: parent_role.to_string(),
                    parent: parent.to_string(),
                    seq: record.seq,
                    dsr_run: opt(ctx.dsr_run.as_ref().map(DsrRunId::as_str)),
                    campaign_run: opt(ctx.campaign_run.as_ref().map(CampaignRunId::as_str)),
                    shard: opt(ctx.shard.as_ref().map(ShardId::as_str)),
                    journey: opt(ctx.journey.as_ref().map(JourneyId::as_str)),
                    case: opt(ctx.case.as_ref().map(CaseId::as_str)),
                },
                None,
            ));
        }
        for gap in &self.gaps {
            events.push(em.emit(
                crate::Severity::Warn,
                crate::EventKind::ContainmentGap {
                    attempt: self.root.as_str().to_string(),
                    node_role: gap.node_role.to_string(),
                    node: gap.node.clone(),
                    missing_parent_role: gap.missing_parent_role.to_string(),
                    missing_parent: gap.missing_parent.clone(),
                },
                None,
            ));
        }
        events
    }
}

/// The V.3.8 embedding seam: the durable-provenance owner places each
/// sealed LOCAL tree into the authoritative global causal DAG. Local trees
/// never carry cross-attempt parentage; anything multi-parent (artifacts,
/// checkpoints, forks, replays) belongs to the embedder's side of this
/// boundary.
pub trait GlobalDagEmbedder {
    /// Refusal type for embeddings the global DAG cannot admit.
    type Error;

    /// Embed one sealed local tree under its propagated Attempt.
    ///
    /// Implementations MUST verify the tree's root matches the attempt
    /// node they are embedding under (see
    /// [`ContainmentError::WrongAttemptRoot`]).
    ///
    /// # Errors
    /// Implementation-defined refusal when the global DAG cannot admit the
    /// tree.
    fn embed(&mut self, tree: &SealedAttemptTree) -> Result<(), Self::Error>;
}

/// Guard for embedder implementations: refuse a tree whose root is not the
/// attempt being embedded under.
///
/// # Errors
/// [`ContainmentError::WrongAttemptRoot`] on mismatch.
pub fn check_embedding_root(
    tree: &SealedAttemptTree,
    under: &AttemptId,
) -> Result<(), ContainmentError> {
    if tree.root() == under {
        Ok(())
    } else {
        Err(ContainmentError::WrongAttemptRoot {
            expected: tree.root().as_str().to_string(),
            found: under.as_str().to_string(),
        })
    }
}

/// One containment line parsed back off the wire (i94v.7.3.2
/// cross-process propagation): either a tree node or a gap-ledger entry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParsedContainmentLine {
    /// A `containment_node` line: the record plus its attempt root.
    Node {
        /// Attempt root named on the line.
        attempt: AttemptId,
        /// Reconstructed typed record.
        record: ContainmentRecord,
    },
    /// A `containment_gap` line.
    Gap {
        /// Attempt root named on the line.
        attempt: AttemptId,
        /// Reconstructed gap entry.
        gap: IntegrityGap,
    },
}

/// Refusals when reading containment lines off the wire.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ContainmentWireError {
    /// The line's wire schema version is not the one this build reads.
    UnsupportedWireVersion {
        /// Version declared on the line.
        declared: u64,
    },
    /// The line is not a containment kind at all.
    NotAContainmentLine,
    /// A required field is missing or malformed.
    MalformedField {
        /// Field key that failed.
        key: &'static str,
    },
}

impl fmt::Display for ContainmentWireError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedWireVersion { declared } => write!(
                f,
                "containment reader accepts wire schema v{} exactly; line declares v{declared}",
                crate::SCHEMA_VERSION
            ),
            Self::NotAContainmentLine => {
                write!(f, "line is not a containment_node/containment_gap event")
            }
            Self::MalformedField { key } => {
                write!(f, "containment line field {key:?} is missing or malformed")
            }
        }
    }
}

impl core::error::Error for ContainmentWireError {}

fn wire_str(line: &str, key: &str) -> Option<String> {
    let tag = format!("\"{key}\":\"");
    let start = line.find(&tag)? + tag.len();
    let rest = &line[start..];
    let mut out = String::new();
    let mut chars = rest.chars();
    while let Some(c) = chars.next() {
        match c {
            '"' => return Some(out),
            '\\' => match chars.next()? {
                '"' => out.push('"'),
                '\\' => out.push('\\'),
                'n' => out.push('\n'),
                't' => out.push('\t'),
                'r' => out.push('\r'),
                _ => return None,
            },
            c => out.push(c),
        }
    }
    None
}

fn wire_u64(line: &str, key: &str) -> Option<u64> {
    let tag = format!("\"{key}\":");
    let start = line.find(&tag)? + tag.len();
    let digits: String = line[start..]
        .chars()
        .take_while(char::is_ascii_digit)
        .collect();
    digits.parse().ok()
}

fn node_id(role: &str, raw: &str) -> Result<LocalNodeId, ContainmentWireError> {
    let malformed = |key| ContainmentWireError::MalformedField { key };
    match role {
        "op" => Ok(LocalNodeId::Op(
            ExecutionOpId::new(raw).map_err(|_| malformed("node"))?,
        )),
        "scope" => Ok(LocalNodeId::Scope(
            ExecutionScopeId::new(raw).map_err(|_| malformed("node"))?,
        )),
        "tile" => Ok(LocalNodeId::Tile(
            TileId::new(raw).map_err(|_| malformed("node"))?,
        )),
        _ => Err(malformed("role")),
    }
}

/// Parse one JSONL line produced by [`SealedAttemptTree::to_events`] back
/// into its typed form. Fail-closed: unknown wire versions and malformed
/// identity fields refuse; non-containment kinds refuse with
/// [`ContainmentWireError::NotAContainmentLine`] so callers can skip
/// interleaved event traffic explicitly.
///
/// # Errors
/// [`ContainmentWireError`] as described above.
pub fn parse_containment_line(line: &str) -> Result<ParsedContainmentLine, ContainmentWireError> {
    let malformed = |key| ContainmentWireError::MalformedField { key };
    let declared = wire_u64(line, "v").ok_or(ContainmentWireError::MalformedField { key: "v" })?;
    if declared != u64::from(crate::SCHEMA_VERSION) {
        return Err(ContainmentWireError::UnsupportedWireVersion { declared });
    }
    let kind = wire_str(line, "kind").ok_or(malformed("kind"))?;
    // Payload keys can collide with envelope keys (the envelope carries its
    // own "seq"), so payload fields are extracted from the payload object
    // only.
    let payload_start = line.find("\"payload\":{").ok_or(malformed("payload"))?;
    let line = &line[payload_start..];
    match kind.as_str() {
        "containment_node" => {
            let attempt = AttemptId::new(wire_str(line, "attempt").ok_or(malformed("attempt"))?)
                .map_err(|_| malformed("attempt"))?;
            let role = wire_str(line, "role").ok_or(malformed("role"))?;
            let node_raw = wire_str(line, "node").ok_or(malformed("node"))?;
            let node = node_id(&role, &node_raw)?;
            let parent_role = wire_str(line, "parent_role").ok_or(malformed("parent_role"))?;
            let parent_raw = wire_str(line, "parent").ok_or(malformed("parent"))?;
            let parent = if parent_role == "attempt" {
                LocalParent::AttemptRoot
            } else {
                LocalParent::Node(node_id(&parent_role, &parent_raw)?)
            };
            let seq = wire_u64(line, "seq").ok_or(malformed("seq"))?;
            let opt = |key: &'static str| -> Result<Option<String>, ContainmentWireError> {
                let raw =
                    wire_str(line, key).ok_or(ContainmentWireError::MalformedField { key })?;
                Ok(if raw.is_empty() { None } else { Some(raw) })
            };
            let context = ContainmentContext {
                dsr_run: opt("dsr_run")?
                    .map(DsrRunId::new)
                    .transpose()
                    .map_err(|_| malformed("dsr_run"))?,
                campaign_run: opt("campaign_run")?
                    .map(CampaignRunId::new)
                    .transpose()
                    .map_err(|_| malformed("campaign_run"))?,
                shard: opt("shard")?
                    .map(ShardId::new)
                    .transpose()
                    .map_err(|_| malformed("shard"))?,
                journey: opt("journey")?
                    .map(JourneyId::new)
                    .transpose()
                    .map_err(|_| malformed("journey"))?,
                case: opt("case")?
                    .map(CaseId::new)
                    .transpose()
                    .map_err(|_| malformed("case"))?,
            };
            Ok(ParsedContainmentLine::Node {
                attempt,
                record: ContainmentRecord {
                    node,
                    parent,
                    seq,
                    context,
                },
            })
        }
        "containment_gap" => {
            let attempt = AttemptId::new(wire_str(line, "attempt").ok_or(malformed("attempt"))?)
                .map_err(|_| malformed("attempt"))?;
            let node_role = wire_str(line, "node_role").ok_or(malformed("node_role"))?;
            let missing_parent_role =
                wire_str(line, "missing_parent_role").ok_or(malformed("missing_parent_role"))?;
            let role_static = |r: &str| -> Result<&'static str, ContainmentWireError> {
                match r {
                    "op" => Ok("op"),
                    "scope" => Ok("scope"),
                    "tile" => Ok("tile"),
                    _ => Err(malformed("node_role")),
                }
            };
            Ok(ParsedContainmentLine::Gap {
                attempt,
                gap: IntegrityGap {
                    node_role: role_static(&node_role)?,
                    node: wire_str(line, "node").ok_or(malformed("node"))?,
                    missing_parent_role: role_static(&missing_parent_role)?,
                    missing_parent: wire_str(line, "missing_parent")
                        .ok_or(malformed("missing_parent"))?,
                },
            })
        }
        _ => Err(ContainmentWireError::NotAContainmentLine),
    }
}
