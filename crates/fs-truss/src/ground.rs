//! Ground structures: a node grid plus every candidate member that
//! passes the fabrication rules -- length bounds and an allowed angle
//! set. The candidate graph is a FrankenNetworkx [`Graph`]. Generation
//! is reproducible, explicitly bounded, and cancellation-aware; the
//! stats row is the compact ledger evidence.

use fnx_classes::Graph;
use fnx_runtime::{CompatibilityMode, RuntimePolicy};
use fs_exec::Cx;
use std::fmt::Write as _;
use std::mem::size_of;

/// Hard ceiling on nodes admitted by one ground-structure construction.
pub const HARD_MAX_GROUND_NODES: usize = 1_000_000;
/// Hard ceiling on candidate pairs examined by grid construction.
pub const HARD_MAX_CANDIDATE_PAIRS: usize = 50_000_000;
/// Hard ceiling on point/member triplets examined for through-node removal.
pub const HARD_MAX_THROUGH_NODE_CHECKS: usize = 250_000_000;
/// Hard ceiling on members retained by one ground structure.
pub const HARD_MAX_GROUND_MEMBERS: usize = 10_000_000;
/// Hard ceiling on conservatively estimated retained ground/graph storage.
pub const HARD_MAX_GROUND_BYTES: usize = 1_073_741_824;
/// Conservative retained FrankenNetworkx bookkeeping estimate per node.
pub const ESTIMATED_GRAPH_BYTES_PER_NODE: usize = 512;
/// Conservative retained FrankenNetworkx bookkeeping estimate per member.
pub const ESTIMATED_GRAPH_BYTES_PER_MEMBER: usize = 512;

/// Smallest grid extent admitted without subnormal-input underflow ambiguity.
pub const MIN_SAFE_GROUND_EXTENT: f64 = f64::MIN_POSITIVE;
/// Conservative coordinate-magnitude limit for two-dimensional predicates.
///
/// This is one quarter of `sqrt(f64::MAX)`. Thus even two coordinates at
/// opposite extremes have differences whose squared products, dot products,
/// and two-term cross products remain finite. The bound is deliberately
/// smaller than the range supported by `hypot`: the through-node predicates
/// also multiply coordinate differences.
pub const MAX_SAFE_GROUND_COORDINATE: f64 = 3.351_951_982_485_649e153;

const CANCELLATION_STRIDE: usize = 256;
const LENGTH_REL_TOLERANCE: f64 = 1e-12;

/// Structured refusal from truss ground or LP construction.
///
/// The variants are shared by the ground-structure and LP admission layers so
/// callers can handle malformed data, resource refusal, allocation pressure,
/// and cancellation without parsing text or observing partial state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrussConstructionError {
    /// A scalar or collection violates its admitted numerical domain.
    InvalidInput {
        /// Stable public field name.
        field: &'static str,
        /// Stable description of the admitted domain.
        requirement: &'static str,
    },
    /// Parallel input vectors do not have the required shape.
    VectorLength {
        /// Stable public field name.
        field: &'static str,
        /// Required element count.
        expected: usize,
        /// Supplied element count.
        actual: usize,
    },
    /// A checked work or retained-resource budget was exceeded.
    WorkBudget {
        /// Stable resource name.
        resource: &'static str,
        /// Admitted maximum.
        limit: usize,
        /// Requested or conservatively estimated amount.
        observed: usize,
    },
    /// A fallible vector reservation failed.
    AllocationFailed {
        /// Stable resource name.
        resource: &'static str,
        /// Number of elements requested for that resource.
        requested: usize,
    },
    /// The caller cancelled construction at a deterministic polling stage.
    Cancelled {
        /// Stable construction stage.
        stage: &'static str,
    },
}

impl core::fmt::Display for TrussConstructionError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput { field, requirement } => {
                write!(formatter, "truss input {field} {requirement}")
            }
            Self::VectorLength {
                field,
                expected,
                actual,
            } => write!(
                formatter,
                "truss input {field} length {actual}; expected {expected}"
            ),
            Self::WorkBudget {
                resource,
                limit,
                observed,
            } => write!(
                formatter,
                "truss construction {resource} budget {limit} exceeded by request {observed}"
            ),
            Self::AllocationFailed {
                resource,
                requested,
            } => write!(
                formatter,
                "truss construction could not reserve {requested} {resource}"
            ),
            Self::Cancelled { stage } => {
                write!(formatter, "truss construction cancelled during {stage}")
            }
        }
    }
}

impl std::error::Error for TrussConstructionError {}

/// Fabrication rules for candidate members.
///
/// Values can only be constructed through [`GroundRules::try_new`] or
/// [`Default`], so every instance has finite positive length bounds, a finite
/// nonnegative angular tolerance, and canonical strictly increasing angles.
#[derive(Debug, Clone)]
pub struct GroundRules {
    min_len: f64,
    max_len: f64,
    angles: Vec<f64>,
    angle_tol: f64,
}

impl GroundRules {
    /// Admit fabrication rules.
    ///
    /// Angles are degrees in `[0, 180)` and must already be in strictly
    /// increasing canonical order. An empty angle collection admits every
    /// direction.
    pub fn try_new(
        min_len: f64,
        max_len: f64,
        angles: Vec<f64>,
        angle_tol: f64,
    ) -> Result<Self, TrussConstructionError> {
        if !min_len.is_finite() || min_len <= 0.0 {
            return Err(TrussConstructionError::InvalidInput {
                field: "min_len",
                requirement: "must be finite and positive",
            });
        }
        if !max_len.is_finite() || max_len < min_len {
            return Err(TrussConstructionError::InvalidInput {
                field: "max_len",
                requirement: "must be finite and at least min_len",
            });
        }
        if !angle_tol.is_finite() || angle_tol < 0.0 {
            return Err(TrussConstructionError::InvalidInput {
                field: "angle_tol",
                requirement: "must be finite and nonnegative",
            });
        }
        let mut previous = None;
        for &angle in &angles {
            if !angle.is_finite()
                || !(0.0..180.0).contains(&angle)
                || (angle == 0.0 && angle.is_sign_negative())
            {
                return Err(TrussConstructionError::InvalidInput {
                    field: "angles",
                    requirement: "must contain only canonical finite angles in [0, 180), including +0",
                });
            }
            if previous.is_some_and(|prior| angle <= prior) {
                return Err(TrussConstructionError::InvalidInput {
                    field: "angles",
                    requirement: "must be in strictly increasing canonical order",
                });
            }
            previous = Some(angle);
        }
        Ok(Self {
            min_len,
            max_len,
            angles,
            angle_tol,
        })
    }

    /// Minimum admitted member length.
    #[must_use]
    pub fn min_len(&self) -> f64 {
        self.min_len
    }

    /// Maximum admitted member length.
    #[must_use]
    pub fn max_len(&self) -> f64 {
        self.max_len
    }

    /// Canonical allowed direction angles in degrees.
    #[must_use]
    pub fn angles(&self) -> &[f64] {
        &self.angles
    }

    /// Angular matching tolerance in degrees.
    #[must_use]
    pub fn angle_tol(&self) -> f64 {
        self.angle_tol
    }
}

impl Default for GroundRules {
    fn default() -> Self {
        Self {
            min_len: 1e-9,
            max_len: f64::MAX,
            angles: Vec::new(),
            angle_tol: 1e-6,
        }
    }
}

/// Caller-selected construction budgets, bounded by crate-wide hard ceilings.
///
/// Small custom limits are useful both for resource-constrained callers and
/// for exact-cap/cap-plus-one conformance tests.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(clippy::struct_field_names)]
pub struct GroundLimits {
    max_nodes: usize,
    max_candidate_pairs: usize,
    max_through_node_checks: usize,
    max_members: usize,
    max_retained_bytes: usize,
}

impl GroundLimits {
    /// Admit a set of positive construction limits.
    pub fn try_new(
        max_nodes: usize,
        max_candidate_pairs: usize,
        max_through_node_checks: usize,
        max_members: usize,
        max_retained_bytes: usize,
    ) -> Result<Self, TrussConstructionError> {
        validate_limit("nodes", max_nodes, HARD_MAX_GROUND_NODES)?;
        validate_limit(
            "candidate_pairs",
            max_candidate_pairs,
            HARD_MAX_CANDIDATE_PAIRS,
        )?;
        validate_limit(
            "through_node_checks",
            max_through_node_checks,
            HARD_MAX_THROUGH_NODE_CHECKS,
        )?;
        validate_limit("members", max_members, HARD_MAX_GROUND_MEMBERS)?;
        validate_limit("retained_bytes", max_retained_bytes, HARD_MAX_GROUND_BYTES)?;
        Ok(Self {
            max_nodes,
            max_candidate_pairs,
            max_through_node_checks,
            max_members,
            max_retained_bytes,
        })
    }

    /// Maximum admitted nodes.
    #[must_use]
    pub fn max_nodes(&self) -> usize {
        self.max_nodes
    }

    /// Maximum candidate node pairs examined.
    #[must_use]
    pub fn max_candidate_pairs(&self) -> usize {
        self.max_candidate_pairs
    }

    /// Maximum point/member triplets examined for through-node removal.
    #[must_use]
    pub fn max_through_node_checks(&self) -> usize {
        self.max_through_node_checks
    }

    /// Maximum members retained in the authoritative output.
    #[must_use]
    pub fn max_members(&self) -> usize {
        self.max_members
    }

    /// Maximum conservative retained ground/graph storage estimate in bytes.
    #[must_use]
    pub fn max_retained_bytes(&self) -> usize {
        self.max_retained_bytes
    }
}

impl Default for GroundLimits {
    fn default() -> Self {
        Self {
            max_nodes: HARD_MAX_GROUND_NODES,
            max_candidate_pairs: HARD_MAX_CANDIDATE_PAIRS,
            max_through_node_checks: HARD_MAX_THROUGH_NODE_CHECKS,
            max_members: HARD_MAX_GROUND_MEMBERS,
            max_retained_bytes: HARD_MAX_GROUND_BYTES,
        }
    }
}

/// An admitted ground structure.
///
/// The parallel member and length arrays are externally immutable. Values are
/// published only after complete validation, graph construction, and the final
/// cancellation checkpoint.
pub struct GroundStructure {
    /// Admitted node positions.
    nodes: Vec<[f64; 2]>,
    /// Members as canonical node index pairs.
    members: Vec<(usize, usize)>,
    /// Positive member lengths parallel to `members`.
    lengths: Vec<f64>,
    /// FrankenNetworkx candidate graph.
    graph: Graph,
}

impl GroundStructure {
    /// Construct an `nx x ny` node grid on `[0, w] x [0, h]` and filter all
    /// pairs by the fabrication rules.
    ///
    /// Admission conservatively checks `nodes choose 2` pair work and
    /// `pairs * (nodes - 2)` through-node triplet work before allocation. It
    /// also bounds conservative retained vector/graph storage and polls `cx` every 256
    /// deterministic loop iterations. No authoritative value is returned
    /// unless the complete graph is built successfully.
    #[allow(clippy::cast_precision_loss, clippy::too_many_lines)]
    pub fn try_grid(
        nx: usize,
        ny: usize,
        w: f64,
        h: f64,
        rules: &GroundRules,
        limits: GroundLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, TrussConstructionError> {
        if nx < 2 {
            return Err(TrussConstructionError::InvalidInput {
                field: "nx",
                requirement: "must be at least 2",
            });
        }
        if ny < 2 {
            return Err(TrussConstructionError::InvalidInput {
                field: "ny",
                requirement: "must be at least 2",
            });
        }
        validate_extent("w", w)?;
        validate_extent("h", h)?;

        let node_count = nx
            .checked_mul(ny)
            .ok_or(TrussConstructionError::WorkBudget {
                resource: "nodes",
                limit: limits.max_nodes,
                observed: usize::MAX,
            })?;
        check_budget("nodes", node_count, limits.max_nodes)?;
        let pair_count = checked_pair_count(node_count, limits.max_candidate_pairs)?;
        check_budget("candidate_pairs", pair_count, limits.max_candidate_pairs)?;
        let through_count = pair_count.checked_mul(node_count.saturating_sub(2)).ok_or(
            TrussConstructionError::WorkBudget {
                resource: "through_node_checks",
                limit: limits.max_through_node_checks,
                observed: usize::MAX,
            },
        )?;
        check_budget(
            "through_node_checks",
            through_count,
            limits.max_through_node_checks,
        )?;
        let member_capacity = pair_count.min(limits.max_members);
        let retained_bytes = retained_vector_bytes(node_count, member_capacity, limits)?;
        check_budget("retained_bytes", retained_bytes, limits.max_retained_bytes)?;
        checkpoint(cx, "grid admission")?;

        let mut nodes = Vec::new();
        reserve_exact(&mut nodes, node_count, "nodes")?;
        let mut node_iteration = 0usize;
        for j in 0..ny {
            for i in 0..nx {
                poll_stride(cx, node_iteration, "grid nodes")?;
                nodes.push([
                    w * i as f64 / (nx - 1) as f64,
                    h * j as f64 / (ny - 1) as f64,
                ]);
                node_iteration += 1;
            }
        }

        let mut members = Vec::new();
        reserve_exact(&mut members, member_capacity, "members")?;
        let mut lengths = Vec::new();
        reserve_exact(&mut lengths, member_capacity, "lengths")?;
        let mut pair_iteration = 0usize;
        let mut through_iteration = 0usize;
        for a in 0..node_count {
            for b in (a + 1)..node_count {
                poll_stride(cx, pair_iteration, "candidate pairs")?;
                pair_iteration += 1;

                let dx = nodes[b][0] - nodes[a][0];
                let dy = nodes[b][1] - nodes[a][1];
                let len = dx.hypot(dy);
                if len < rules.min_len || len > rules.max_len {
                    continue;
                }
                if !rules.angles.is_empty() {
                    let angle = dy.atan2(dx).to_degrees().rem_euclid(180.0);
                    if !angle_is_allowed(angle, &rules.angles, rules.angle_tol) {
                        continue;
                    }
                }

                // A longer collinear member that passes exactly through a
                // grid node adds no independent statics and is omitted.
                let mut through = false;
                for (c, node) in nodes.iter().enumerate() {
                    if c == a || c == b {
                        continue;
                    }
                    poll_stride(cx, through_iteration, "through-node checks")?;
                    through_iteration += 1;
                    let cx_offset = node[0] - nodes[a][0];
                    let cy_offset = node[1] - nodes[a][1];
                    let cross = cx_offset * dy - cy_offset * dx;
                    let dot = cx_offset * dx + cy_offset * dy;
                    if cross.abs() < 1e-9 * len && dot > 1e-12 && dot < len * len - 1e-12 {
                        through = true;
                        break;
                    }
                }
                if through {
                    continue;
                }
                if members.len() == limits.max_members {
                    return Err(TrussConstructionError::WorkBudget {
                        resource: "members",
                        limit: limits.max_members,
                        observed: members.len().saturating_add(1),
                    });
                }
                members.push((a, b));
                lengths.push(len);
            }
        }

        let graph = build_graph(node_count, &members, cx)?;
        checkpoint(cx, "ground publication")?;
        Ok(Self {
            nodes,
            members,
            lengths,
            graph,
        })
    }

    /// Construct an admitted ground structure from borrowed hand-fixture data.
    ///
    /// Nodes must be nonempty and finite within
    /// [`MAX_SAFE_GROUND_COORDINATE`]. Members must be in strict lexicographic
    /// order with canonical endpoints `a < b`, and `lengths` must contain the
    /// corresponding positive finite Euclidean lengths within a relative
    /// tolerance of `1e-12`. Published lengths are always the recomputed
    /// Euclidean values, so within-tolerance caller representations cannot
    /// create distinct authoritative identities.
    #[allow(clippy::too_many_lines)]
    pub fn try_from_parts(
        nodes: &[[f64; 2]],
        members: &[(usize, usize)],
        lengths: &[f64],
        limits: GroundLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, TrussConstructionError> {
        if nodes.is_empty() {
            return Err(TrussConstructionError::InvalidInput {
                field: "nodes",
                requirement: "must be nonempty",
            });
        }
        if lengths.len() != members.len() {
            return Err(TrussConstructionError::VectorLength {
                field: "lengths",
                expected: members.len(),
                actual: lengths.len(),
            });
        }
        check_budget("nodes", nodes.len(), limits.max_nodes)?;
        check_budget("candidate_pairs", members.len(), limits.max_candidate_pairs)?;
        check_budget("members", members.len(), limits.max_members)?;
        let retained_bytes = retained_vector_bytes(nodes.len(), members.len(), limits)?;
        check_budget("retained_bytes", retained_bytes, limits.max_retained_bytes)?;
        checkpoint(cx, "parts admission")?;

        for (index, point) in nodes.iter().enumerate() {
            poll_stride(cx, index, "parts nodes")?;
            if point.iter().any(|coordinate| {
                !coordinate.is_finite() || coordinate.abs() > MAX_SAFE_GROUND_COORDINATE
            }) {
                return Err(TrussConstructionError::InvalidInput {
                    field: "nodes",
                    requirement: "must contain only finite coordinates within the safe bound",
                });
            }
        }

        let mut canonical_lengths = Vec::new();
        reserve_exact(&mut canonical_lengths, lengths.len(), "lengths")?;
        let mut previous = None;
        for (index, (&(a, b), &stored_length)) in members.iter().zip(lengths.iter()).enumerate() {
            poll_stride(cx, index, "parts members")?;
            if a >= b {
                return Err(TrussConstructionError::InvalidInput {
                    field: "members",
                    requirement: "must have canonical endpoints a < b",
                });
            }
            if b >= nodes.len() {
                return Err(TrussConstructionError::InvalidInput {
                    field: "members",
                    requirement: "must contain only in-range endpoints",
                });
            }
            let member = (a, b);
            if previous.is_some_and(|prior| member <= prior) {
                return Err(TrussConstructionError::InvalidInput {
                    field: "members",
                    requirement: "must be in strictly increasing canonical order without duplicates",
                });
            }
            previous = Some(member);

            if !stored_length.is_finite() || stored_length <= 0.0 {
                return Err(TrussConstructionError::InvalidInput {
                    field: "lengths",
                    requirement: "must contain only finite positive values",
                });
            }
            let dx = nodes[b][0] - nodes[a][0];
            let dy = nodes[b][1] - nodes[a][1];
            let recomputed = dx.hypot(dy);
            if !recomputed.is_finite() || recomputed <= 0.0 {
                return Err(TrussConstructionError::InvalidInput {
                    field: "members",
                    requirement: "must connect geometrically distinct safe nodes",
                });
            }
            let tolerance = LENGTH_REL_TOLERANCE * stored_length.max(recomputed);
            if (stored_length - recomputed).abs() > tolerance {
                return Err(TrussConstructionError::InvalidInput {
                    field: "lengths",
                    requirement: "must agree with recomputed Euclidean member lengths",
                });
            }
            canonical_lengths.push(recomputed);
        }

        let mut owned_nodes = Vec::new();
        reserve_exact(&mut owned_nodes, nodes.len(), "nodes")?;
        owned_nodes.extend_from_slice(nodes);
        let mut owned_members = Vec::new();
        reserve_exact(&mut owned_members, members.len(), "members")?;
        owned_members.extend_from_slice(members);
        let graph = build_graph(owned_nodes.len(), &owned_members, cx)?;
        checkpoint(cx, "ground publication")?;
        Ok(Self {
            nodes: owned_nodes,
            members: owned_members,
            lengths: canonical_lengths,
            graph,
        })
    }

    /// Admitted node coordinates.
    #[must_use]
    pub fn nodes(&self) -> &[[f64; 2]] {
        &self.nodes
    }

    /// Canonical admitted member endpoints.
    #[must_use]
    pub fn members(&self) -> &[(usize, usize)] {
        &self.members
    }

    /// Positive member lengths parallel to [`Self::members`].
    #[must_use]
    pub fn lengths(&self) -> &[f64] {
        &self.lengths
    }

    /// FrankenNetworkx candidate graph.
    #[must_use]
    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    /// Ledger stats row (counts plus FNV hash of the member list).
    #[must_use]
    pub fn stats(&self) -> String {
        let mut hash = 0xcbf2_9ce4_8422_2325u64;
        let mut mix = |value: u64| {
            for byte in value.to_le_bytes() {
                hash ^= u64::from(byte);
                hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
            }
        };
        for &(a, b) in &self.members {
            mix(a as u64);
            mix(b as u64);
        }
        let mut output = String::new();
        let _ = write!(
            output,
            "{{\"nodes\":{},\"members\":{},\"hash\":\"{hash:#018x}\"}}",
            self.nodes.len(),
            self.members.len()
        );
        output
    }
}

fn validate_limit(
    resource: &'static str,
    value: usize,
    hard_limit: usize,
) -> Result<(), TrussConstructionError> {
    if value == 0 {
        return Err(TrussConstructionError::InvalidInput {
            field: resource,
            requirement: "limit must be positive",
        });
    }
    check_budget(resource, value, hard_limit)
}

fn validate_extent(field: &'static str, value: f64) -> Result<(), TrussConstructionError> {
    if !value.is_finite() || value < MIN_SAFE_GROUND_EXTENT || value > MAX_SAFE_GROUND_COORDINATE {
        return Err(TrussConstructionError::InvalidInput {
            field,
            requirement: "must be finite and within the safe positive extent bounds",
        });
    }
    Ok(())
}

fn angle_is_allowed(angle: f64, allowed: &[f64], tolerance: f64) -> bool {
    debug_assert!(!allowed.is_empty());
    let insertion = allowed.partition_point(|&candidate| candidate < angle);
    let matches = |candidate: f64| {
        let distance = (angle - candidate).abs();
        distance.min(180.0 - distance) <= tolerance
    };
    (insertion < allowed.len() && matches(allowed[insertion]))
        || (insertion > 0 && matches(allowed[insertion - 1]))
        || matches(allowed[0])
        || matches(allowed[allowed.len() - 1])
}

fn check_budget(
    resource: &'static str,
    observed: usize,
    limit: usize,
) -> Result<(), TrussConstructionError> {
    if observed > limit {
        Err(TrussConstructionError::WorkBudget {
            resource,
            limit,
            observed,
        })
    } else {
        Ok(())
    }
}

fn checked_pair_count(nodes: usize, limit: usize) -> Result<usize, TrussConstructionError> {
    nodes
        .checked_mul(nodes.saturating_sub(1))
        .and_then(|product| product.checked_div(2))
        .ok_or(TrussConstructionError::WorkBudget {
            resource: "candidate_pairs",
            limit,
            observed: usize::MAX,
        })
}

fn retained_vector_bytes(
    nodes: usize,
    members: usize,
    limits: GroundLimits,
) -> Result<usize, TrussConstructionError> {
    let node_bytes = nodes.checked_mul(
        size_of::<[f64; 2]>()
            .checked_add(ESTIMATED_GRAPH_BYTES_PER_NODE)
            .ok_or(TrussConstructionError::WorkBudget {
                resource: "retained_bytes",
                limit: limits.max_retained_bytes,
                observed: usize::MAX,
            })?,
    );
    let member_bytes = members.checked_mul(
        size_of::<(usize, usize)>()
            .checked_add(size_of::<f64>())
            .and_then(|bytes| bytes.checked_add(ESTIMATED_GRAPH_BYTES_PER_MEMBER))
            .ok_or(TrussConstructionError::WorkBudget {
                resource: "retained_bytes",
                limit: limits.max_retained_bytes,
                observed: usize::MAX,
            })?,
    );
    node_bytes
        .and_then(|total| member_bytes.and_then(|bytes| total.checked_add(bytes)))
        .ok_or(TrussConstructionError::WorkBudget {
            resource: "retained_bytes",
            limit: limits.max_retained_bytes,
            observed: usize::MAX,
        })
}

fn reserve_exact<T>(
    vector: &mut Vec<T>,
    requested: usize,
    resource: &'static str,
) -> Result<(), TrussConstructionError> {
    vector
        .try_reserve_exact(requested)
        .map_err(|_| TrussConstructionError::AllocationFailed {
            resource,
            requested,
        })
}

fn checkpoint(cx: &Cx<'_>, stage: &'static str) -> Result<(), TrussConstructionError> {
    cx.checkpoint()
        .map_err(|_| TrussConstructionError::Cancelled { stage })
}

fn poll_stride(
    cx: &Cx<'_>,
    iteration: usize,
    stage: &'static str,
) -> Result<(), TrussConstructionError> {
    if iteration.is_multiple_of(CANCELLATION_STRIDE) {
        checkpoint(cx, stage)?;
    }
    Ok(())
}

fn build_graph(
    node_count: usize,
    members: &[(usize, usize)],
    cx: &Cx<'_>,
) -> Result<Graph, TrussConstructionError> {
    let mut graph = Graph::new(CompatibilityMode::Strict);
    for node_start in (0..node_count).step_by(CANCELLATION_STRIDE) {
        checkpoint(cx, "graph nodes")?;
        let node_end = node_start
            .saturating_add(CANCELLATION_STRIDE)
            .min(node_count);
        let inserted =
            graph.extend_nodes_unrecorded((node_start..node_end).map(|node| format!("n{node}")));
        if inserted != node_end - node_start {
            return Err(TrussConstructionError::InvalidInput {
                field: "candidate graph",
                requirement: "must accept every canonical node identity exactly once",
            });
        }
    }
    for member_chunk in members.chunks(CANCELLATION_STRIDE) {
        checkpoint(cx, "graph members")?;
        let inserted = graph.extend_existing_index_edges_unrecorded(member_chunk.iter().copied());
        if inserted != member_chunk.len() {
            return Err(TrussConstructionError::InvalidInput {
                field: "candidate graph",
                requirement: "must accept every canonical member identity exactly once",
            });
        }
    }
    // Graph mutation records wall-clock compatibility evidence internally.
    // Construction evidence is not part of the mathematical candidate graph,
    // so clear it before publication to make returned graph state replayable.
    graph.set_runtime_policy(RuntimePolicy::new(CompatibilityMode::Strict));
    Ok(graph)
}
