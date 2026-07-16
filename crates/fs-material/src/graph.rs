//! The executable ConstitutiveGraph and law-node protocol (bead kagp,
//! increment 1).
//!
//! Matter is a typed constitutive graph, not a bag of scalars. The
//! seven-role decomposition is load-bearing: (1) topology/balance and
//! (2) bulk storage/metric are OWNED ELSEWHERE (fs-feec/fs-rep-mesh) —
//! this graph can reference them as external roles but refuses to
//! execute them; (3) bulk transport/dissipation, (4) reversible/
//! cross-coupled blocks, (5) interface laws, (6) reaction/source, and
//! (7) internal memory are executable law nodes here.
//!
//! Every node DECLARES its ports (name + dims), state schema,
//! calibration domain, differentiability class, energy/dissipation
//! behavior, and whether it claims a consistent tangent — and admission
//! refuses nodes whose declarations are incomplete or whose claims
//! probe false, with typed diagnostics naming the node, the law, and
//! the failed obligation. Implementations are keyed by the immutable
//! L1 fs-matdb (LawId, LawVersion) identity and instantiated from
//! validated `ConstitutiveModelCard`s — fs-material CONSUMES the card
//! metadata, never redefines it.

use std::collections::BTreeMap;
use std::fmt;

use fs_matdb::{ConstitutiveModelCard, LawId};
use fs_qty::Dims;

/// The seven constitutive roles. The first two are structural anchors
/// owned by other crates; nodes claiming them are declarable (so a
/// graph can NAME its balance/storage context) but never executable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeRole {
    /// Exact discrete d / balance topology (fs-feec/fs-rep-mesh owns).
    TopologyBalance,
    /// Weighted mass / Hodge storage metric (fs-feec owns).
    BulkStorage,
    /// Monotone flux operators / dissipation potentials (Fourier, Ohm,
    /// Darcy, viscosity).
    BulkTransport,
    /// Reversible and cross-coupled blocks (Hall, gyroscopic, piezo,
    /// thermoelectric) with declared skew part and reciprocity class.
    ReversibleCoupling,
    /// Oriented trace relations at interfaces (friction, wetting,
    /// contact conductance, radiation exchange).
    InterfaceLaw,
    /// Stoichiometric / source operators (combustion, electrochemistry,
    /// impressed current) under conservation constraints.
    ReactionSource,
    /// State-transition laws with memory (plasticity, damage,
    /// hysteresis, viscoelasticity, fatigue, wear).
    InternalMemory,
}

impl NodeRole {
    /// Roles owned by other crates: declarable, never executable here.
    #[must_use]
    pub fn externally_owned(self) -> bool {
        matches!(self, NodeRole::TopologyBalance | NodeRole::BulkStorage)
    }

    /// Stable diagnostic tag.
    #[must_use]
    pub fn tag(self) -> &'static str {
        match self {
            NodeRole::TopologyBalance => "topology-balance",
            NodeRole::BulkStorage => "bulk-storage",
            NodeRole::BulkTransport => "bulk-transport",
            NodeRole::ReversibleCoupling => "reversible-coupling",
            NodeRole::InterfaceLaw => "interface-law",
            NodeRole::ReactionSource => "reaction-source",
            NodeRole::InternalMemory => "internal-memory",
        }
    }
}

/// Differentiability class a node declares for its maps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Differentiability {
    /// C¹ or better on the calibration domain.
    Smooth,
    /// Piecewise smooth with declared kink sets (return mapping).
    PiecewiseSmooth,
    /// No derivative claim (lookup/hysteresis operators).
    NonSmooth,
}

/// Declared energy/dissipation behavior. This is a CONTRACT surface:
/// the runtime audits dissipation non-negativity when declared, and a
/// node with `Empirical` says out loud that no admissibility claim is
/// made.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyBehavior {
    /// Derives from a free-energy potential; storage is conservative.
    FreeEnergyStorage,
    /// Purely dissipative; reported rate must be non-negative.
    NonNegativeDissipation,
    /// Mixed storage + dissipation, both audited.
    StorageAndDissipation,
    /// EXPLICIT empirical no-claim: no thermodynamic admissibility is
    /// asserted (the honest label for fitted hysteresis).
    Empirical,
}

/// Time-reversal parity of a port's quantity (Onsager–Casimir): heat
/// flux and temperature gradient are Even; magnetization rates,
/// velocities, and Hall-type variables are Odd. Required on
/// reversible-coupling nodes, where the reciprocity gate reads it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeParity {
    /// Invariant under time reversal.
    Even,
    /// Sign-flips under time reversal.
    Odd,
}

/// One typed port: a named quantity with dimensions and time parity.
#[derive(Debug, Clone, PartialEq)]
pub struct Port {
    /// Port name (unique within its direction on a node).
    pub name: String,
    /// The quantity's dimensions.
    pub dims: Dims,
    /// Time-reversal parity (read by the Onsager–Casimir gate on
    /// reversible-coupling nodes; informational elsewhere).
    pub parity: TimeParity,
}

/// A node's complete declaration. Admission validates every field; an
/// admitted declaration is the node's public contract.
#[derive(Debug, Clone, PartialEq)]
pub struct NodeDeclaration {
    /// The law identity this node implements (fs-matdb keyed).
    pub law: LawId,
    /// The law's semantic version.
    pub law_version: u32,
    /// The constitutive role.
    pub role: NodeRole,
    /// Input ports (nonempty for executable roles).
    pub inputs: Vec<Port>,
    /// Output ports (nonempty for executable roles).
    pub outputs: Vec<Port>,
    /// Named internal-state slots, in codec order (empty = stateless).
    pub state_slots: Vec<String>,
    /// The state schema version the slots follow.
    pub state_schema_version: u32,
    /// Where the node's calibration is valid.
    pub calibration: fs_evidence::ValidityDomain,
    /// Declared differentiability class.
    pub differentiability: Differentiability,
    /// Declared energy/dissipation behavior.
    pub energy: EnergyBehavior,
    /// Whether the node claims a consistent tangent (d outputs / d
    /// inputs, algorithmically consistent with `evaluate`).
    pub tangent_claimed: bool,
}

/// One evaluation's result.
#[derive(Debug, Clone, PartialEq)]
pub struct NodeOutput {
    /// Output port values, in declaration order.
    pub outputs: Vec<f64>,
    /// Updated state, in slot order (same length as `state_slots`).
    pub next_state: Vec<f64>,
    /// Reported dissipation rate (required non-negative when the node
    /// declares dissipative behavior; `None` = not reported).
    pub dissipation_rate: Option<f64>,
}

/// The law-node protocol: a declaration plus the executable maps.
pub trait LawNode {
    /// The node's declaration (stable across calls).
    fn declaration(&self) -> &NodeDeclaration;

    /// Evaluate the node: `inputs` in input-port order, `state` in
    /// slot order.
    ///
    /// # Errors
    /// [`GraphError`] for arity mismatches or non-finite payloads (the
    /// node's own numerical refusals).
    fn evaluate(&self, state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError>;

    /// The consistent tangent d(outputs)/d(inputs) at (state, inputs),
    /// row-major `[outputs.len() × inputs.len()]`. MUST be provided
    /// exactly when the declaration claims it.
    fn tangent(&self, _state: &[f64], _inputs: &[f64]) -> Option<Vec<f64>> {
        None
    }

    /// Free energy ψ(state, inputs) for nodes whose declared energy
    /// behavior claims storage (`FreeEnergyStorage` /
    /// `StorageAndDissipation`). MUST be provided exactly when storage
    /// is claimed; the Hessian gates differentiate it.
    fn free_energy(&self, _state: &[f64], _inputs: &[f64]) -> Option<f64> {
        None
    }
}

/// Everything the graph boundary can refuse, always naming the node,
/// the law, and the failed obligation.
#[derive(Debug, Clone, PartialEq)]
pub enum GraphError {
    /// A declaration field is structurally incomplete.
    IncompleteDeclaration {
        /// The node id (registry/graph name).
        node: String,
        /// The law id.
        law: String,
        /// The failed obligation.
        obligation: &'static str,
    },
    /// An externally-owned role was asked to execute.
    ExternallyOwnedRole {
        /// The node id.
        node: String,
        /// The role.
        role: &'static str,
    },
    /// The node claims a tangent but its probe returned none or a
    /// wrongly-sized matrix.
    TangentClaimUnmet {
        /// The node id.
        node: String,
        /// The law id.
        law: String,
    },
    /// The node declares free-energy storage but provides no ψ
    /// evaluator.
    EnergyClaimUnmet {
        /// The node id.
        node: String,
        /// The law id.
        law: String,
    },
    /// Evaluation arity mismatch (state or input slice length).
    ArityMismatch {
        /// The node id.
        node: String,
        /// Which slice.
        which: &'static str,
        /// Expected length.
        expected: usize,
        /// Offered length.
        found: usize,
    },
    /// A payload value is non-finite.
    NonFinite {
        /// The node id.
        node: String,
        /// Which field.
        field: &'static str,
        /// The offending bits.
        bits: u64,
    },
    /// A declared-dissipative node reported a negative rate.
    NegativeDissipation {
        /// The node id.
        node: String,
        /// The reported rate.
        rate: f64,
    },
    /// A card's identity does not match the registered factory or the
    /// node it built.
    CardMismatch {
        /// The law id offered.
        law: String,
        /// What disagreed.
        obligation: &'static str,
    },
    /// No factory is registered for the (law, version) identity.
    UnknownLaw {
        /// The law id.
        law: String,
        /// The version.
        version: u32,
    },
    /// A required card parameter is missing.
    MissingParameter {
        /// The law id.
        law: String,
        /// The parameter name.
        parameter: String,
    },
    /// An edge connects dimensionally incompatible ports.
    EdgeDimsMismatch {
        /// Source node id.
        from: String,
        /// Source port.
        from_port: String,
        /// Target node id.
        to: String,
        /// Target port.
        to_port: String,
    },
    /// An edge references a node or port that does not exist.
    UnknownEndpoint {
        /// The dangling reference.
        reference: String,
    },
    /// The graph contains a cycle (increment 1 executes DAGs only;
    /// implicit/history coupling belongs to the solver loop, not the
    /// graph's single pass).
    CycleDetected {
        /// A node on the cycle.
        node: String,
    },
    /// Aggregate state codec refusals.
    StateSchemaMismatch {
        /// What disagreed.
        obligation: &'static str,
        /// Expected value.
        expected: u64,
        /// Found value.
        found: u64,
    },
    /// Underlying card validation failed.
    Card(fs_matdb::MatDbError),
}

impl fmt::Display for GraphError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GraphError::IncompleteDeclaration {
                node,
                law,
                obligation,
            } => write!(
                f,
                "node '{node}' (law '{law}'): incomplete declaration: {obligation}"
            ),
            GraphError::ExternallyOwnedRole { node, role } => write!(
                f,
                "node '{node}': role '{role}' is owned by fs-feec/fs-rep-mesh and cannot execute \
                 here"
            ),
            GraphError::TangentClaimUnmet { node, law } => write!(
                f,
                "node '{node}' (law '{law}') claims a consistent tangent but does not provide one"
            ),
            GraphError::EnergyClaimUnmet { node, law } => write!(
                f,
                "node '{node}' (law '{law}') declares free-energy storage but provides no \
                 evaluator"
            ),
            GraphError::ArityMismatch {
                node,
                which,
                expected,
                found,
            } => write!(
                f,
                "node '{node}': {which} arity {found} (expected {expected})"
            ),
            GraphError::NonFinite { node, field, bits } => {
                write!(
                    f,
                    "node '{node}': {field} is non-finite (bits {bits:#018x})"
                )
            }
            GraphError::NegativeDissipation { node, rate } => write!(
                f,
                "node '{node}': declared-dissipative law reported negative rate {rate}"
            ),
            GraphError::CardMismatch { law, obligation } => {
                write!(f, "card for law '{law}' mismatches: {obligation}")
            }
            GraphError::UnknownLaw { law, version } => {
                write!(f, "no factory registered for law '{law}' v{version}")
            }
            GraphError::MissingParameter { law, parameter } => {
                write!(
                    f,
                    "law '{law}': card lacks required parameter '{parameter}'"
                )
            }
            GraphError::EdgeDimsMismatch {
                from,
                from_port,
                to,
                to_port,
            } => write!(
                f,
                "edge {from}.{from_port} -> {to}.{to_port} connects dimensionally incompatible \
                 ports"
            ),
            GraphError::UnknownEndpoint { reference } => {
                write!(f, "edge references unknown endpoint '{reference}'")
            }
            GraphError::CycleDetected { node } => write!(
                f,
                "graph cycle through '{node}': single-pass execution is DAG-only; implicit \
                 coupling belongs to the solver loop"
            ),
            GraphError::StateSchemaMismatch {
                obligation,
                expected,
                found,
            } => write!(
                f,
                "aggregate state codec: {obligation} (expected {expected}, found {found})"
            ),
            GraphError::Card(inner) => write!(f, "card validation failed: {inner}"),
        }
    }
}

impl std::error::Error for GraphError {}

impl From<fs_matdb::MatDbError> for GraphError {
    fn from(inner: fs_matdb::MatDbError) -> GraphError {
        GraphError::Card(inner)
    }
}

/// Admit a node under a graph-local id: validate the declaration and
/// probe the tangent claim. Every refusal names node, law, and the
/// failed obligation.
///
/// # Errors
/// Typed [`GraphError`] refusals; nothing partial.
pub fn admit_node(node_id: &str, node: &dyn LawNode) -> Result<(), GraphError> {
    let declaration = node.declaration();
    let law = declaration.law.0.clone();
    let incomplete = |obligation: &'static str| GraphError::IncompleteDeclaration {
        node: node_id.to_string(),
        law: law.clone(),
        obligation,
    };
    if declaration.role.externally_owned() {
        return Err(GraphError::ExternallyOwnedRole {
            node: node_id.to_string(),
            role: declaration.role.tag(),
        });
    }
    if declaration.law.0.trim().is_empty() {
        return Err(incomplete("law id is blank"));
    }
    if declaration.inputs.is_empty() {
        return Err(incomplete(
            "an executable node needs at least one input port",
        ));
    }
    if declaration.outputs.is_empty() {
        return Err(incomplete(
            "an executable node needs at least one output port",
        ));
    }
    let mut seen = BTreeMap::new();
    for (direction, ports) in [
        ("input", &declaration.inputs),
        ("output", &declaration.outputs),
    ] {
        seen.clear();
        for port in ports {
            if port.name.trim().is_empty() {
                return Err(incomplete("a port has a blank name"));
            }
            if seen.insert(port.name.clone(), ()).is_some() {
                return Err(match direction {
                    "input" => incomplete("duplicate input port name"),
                    _ => incomplete("duplicate output port name"),
                });
            }
        }
    }
    seen.clear();
    for slot in &declaration.state_slots {
        if slot.trim().is_empty() {
            return Err(incomplete("a state slot has a blank name"));
        }
        if seen.insert(slot.clone(), ()).is_some() {
            return Err(incomplete("duplicate state slot name"));
        }
    }
    for (axis, &(lo, hi)) in declaration.calibration.bounds() {
        if lo.is_nan() || hi.is_nan() {
            let _ = axis;
            return Err(incomplete("calibration domain has NaN bounds"));
        }
    }
    if declaration.tangent_claimed {
        let state = vec![0.0; declaration.state_slots.len()];
        let inputs = vec![0.0; declaration.inputs.len()];
        match node.tangent(&state, &inputs) {
            Some(matrix)
                if matrix.len() == declaration.outputs.len() * declaration.inputs.len() => {}
            _ => {
                return Err(GraphError::TangentClaimUnmet {
                    node: node_id.to_string(),
                    law,
                });
            }
        }
    }
    if matches!(
        declaration.energy,
        EnergyBehavior::FreeEnergyStorage | EnergyBehavior::StorageAndDissipation
    ) {
        let state = vec![0.0; declaration.state_slots.len()];
        let inputs = vec![0.0; declaration.inputs.len()];
        if node.free_energy(&state, &inputs).is_none() {
            return Err(GraphError::EnergyClaimUnmet {
                node: node_id.to_string(),
                law,
            });
        }
    }
    Ok(())
}

/// Gate a storage-claiming node's energy consistency at a point:
/// (a) the outputs are the conjugate forces ∂ψ/∂inputs (central FD),
/// and (b) the tangent — which for such a node IS the Hessian of ψ —
/// is symmetric within `atol` (Maxwell reciprocity).
///
/// # Errors
/// [`GraphError::EnergyClaimUnmet`] when ψ is not provided;
/// [`GraphError::TangentClaimUnmet`] when the tangent is missing;
/// [`GraphError::IncompleteDeclaration`] naming the failed obligation
/// (gradient mismatch or Hessian asymmetry).
pub fn check_free_energy_consistency(
    node_id: &str,
    node: &dyn LawNode,
    state: &[f64],
    inputs: &[f64],
    atol: f64,
) -> Result<(), GraphError> {
    let declaration = node.declaration();
    let law = declaration.law.0.clone();
    if node.free_energy(state, inputs).is_none() {
        return Err(GraphError::EnergyClaimUnmet {
            node: node_id.to_string(),
            law,
        });
    }
    let n = declaration.inputs.len();
    let m = declaration.outputs.len();
    let h = 1e-6;
    let mut probe = inputs.to_vec();
    let value = node.evaluate(state, inputs)?;
    for j in 0..n.min(m) {
        let base = probe[j];
        probe[j] = base + h;
        let plus = node
            .free_energy(state, &probe)
            .ok_or_else(|| GraphError::EnergyClaimUnmet {
                node: node_id.to_string(),
                law: law.clone(),
            })?;
        probe[j] = base - h;
        let minus =
            node.free_energy(state, &probe)
                .ok_or_else(|| GraphError::EnergyClaimUnmet {
                    node: node_id.to_string(),
                    law: law.clone(),
                })?;
        probe[j] = base;
        let gradient = (plus - minus) / (2.0 * h);
        if (gradient - value.outputs[j]).abs() > atol {
            return Err(GraphError::IncompleteDeclaration {
                node: node_id.to_string(),
                law: law.clone(),
                obligation: "outputs are not the conjugate forces of the declared free energy",
            });
        }
    }
    let Some(tangent) = node.tangent(state, inputs) else {
        return Err(GraphError::TangentClaimUnmet {
            node: node_id.to_string(),
            law,
        });
    };
    check_symmetry(
        node_id,
        &law,
        &tangent,
        m,
        n,
        atol,
        "free-energy Hessian is asymmetric",
    )
}

/// Gate the PSD symmetric part of a node's tangent: the second-law
/// obligation for dissipative/transport blocks — the symmetric part of
/// the conductivity/coupling matrix must be positive semidefinite.
/// Checked by Sylvester's criterion on the symmetric part (leading
/// principal minors `>= -atol`); fixture-scale `n` only.
///
/// # Errors
/// [`GraphError::TangentClaimUnmet`] when no tangent is provided;
/// [`GraphError::IncompleteDeclaration`] when a leading principal minor
/// is negative beyond `atol`.
pub fn check_psd_symmetric_part(
    node_id: &str,
    node: &dyn LawNode,
    state: &[f64],
    inputs: &[f64],
    atol: f64,
) -> Result<(), GraphError> {
    let declaration = node.declaration();
    let law = declaration.law.0.clone();
    let n = declaration.inputs.len();
    let m = declaration.outputs.len();
    let Some(tangent) = node.tangent(state, inputs) else {
        return Err(GraphError::TangentClaimUnmet {
            node: node_id.to_string(),
            law,
        });
    };
    if m != n {
        return Err(GraphError::IncompleteDeclaration {
            node: node_id.to_string(),
            law,
            obligation: "PSD audit needs a square tangent (matched input/output ports)",
        });
    }
    // Symmetric part, then leading principal minors by Gaussian
    // elimination without pivoting (adequate at fixture scale).
    let mut s = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            s[i * n + j] = f64::midpoint(tangent[i * n + j], tangent[j * n + i]);
        }
    }
    for k in 1..=n {
        if minor_determinant(&s, n, k) < -atol {
            return Err(GraphError::IncompleteDeclaration {
                node: node_id.to_string(),
                law,
                obligation: "symmetric dissipative part is not positive semidefinite",
            });
        }
    }
    Ok(())
}

/// Gate Onsager–Casimir reciprocity on a reversible-coupling block at
/// zero applied field: `L[i][j] = e_i * e_j * L[j][i]` with `e` the
/// declared time parities (+1 Even, −1 Odd). Even–even pairs must be
/// symmetric (Onsager); mixed-parity pairs antisymmetric (Casimir).
///
/// # Errors
/// [`GraphError::TangentClaimUnmet`] when no tangent is provided;
/// [`GraphError::IncompleteDeclaration`] when the block is not square
/// or a pair violates its reciprocity class beyond `atol`.
pub fn check_onsager_casimir(
    node_id: &str,
    node: &dyn LawNode,
    state: &[f64],
    inputs: &[f64],
    atol: f64,
) -> Result<(), GraphError> {
    let declaration = node.declaration();
    let law = declaration.law.0.clone();
    let n = declaration.inputs.len();
    let m = declaration.outputs.len();
    if m != n {
        return Err(GraphError::IncompleteDeclaration {
            node: node_id.to_string(),
            law,
            obligation: "reciprocity audit needs a square coupling block",
        });
    }
    let Some(tangent) = node.tangent(state, inputs) else {
        return Err(GraphError::TangentClaimUnmet {
            node: node_id.to_string(),
            law,
        });
    };
    let sign = |parity: TimeParity| -> f64 {
        match parity {
            TimeParity::Even => 1.0,
            TimeParity::Odd => -1.0,
        }
    };
    for i in 0..n {
        for j in 0..n {
            let expected = sign(declaration.outputs[i].parity)
                * sign(declaration.outputs[j].parity)
                * tangent[j * n + i];
            if (tangent[i * n + j] - expected).abs() > atol {
                return Err(GraphError::IncompleteDeclaration {
                    node: node_id.to_string(),
                    law,
                    obligation: "coupling block violates Onsager-Casimir reciprocity",
                });
            }
        }
    }
    Ok(())
}

fn check_symmetry(
    node_id: &str,
    law: &str,
    tangent: &[f64],
    m: usize,
    n: usize,
    atol: f64,
    obligation: &'static str,
) -> Result<(), GraphError> {
    if m != n {
        return Err(GraphError::IncompleteDeclaration {
            node: node_id.to_string(),
            law: law.to_string(),
            obligation: "symmetry audit needs a square tangent",
        });
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if (tangent[i * n + j] - tangent[j * n + i]).abs() > atol {
                return Err(GraphError::IncompleteDeclaration {
                    node: node_id.to_string(),
                    law: law.to_string(),
                    obligation,
                });
            }
        }
    }
    Ok(())
}

/// Leading principal minor determinant of order `k` (Laplace at fixture
/// scale).
fn minor_determinant(s: &[f64], n: usize, k: usize) -> f64 {
    let mut a = vec![0.0f64; k * k];
    for i in 0..k {
        for j in 0..k {
            a[i * k + j] = s[i * n + j];
        }
    }
    det(&mut a, k)
}

fn det(a: &mut [f64], k: usize) -> f64 {
    let mut result = 1.0;
    for col in 0..k {
        // Partial pivot for stability at fixture scale.
        let mut pivot = col;
        for row in (col + 1)..k {
            if a[row * k + col].abs() > a[pivot * k + col].abs() {
                pivot = row;
            }
        }
        if a[pivot * k + col] == 0.0 {
            return 0.0;
        }
        if pivot != col {
            for j in 0..k {
                a.swap(col * k + j, pivot * k + j);
            }
            result = -result;
        }
        result *= a[col * k + col];
        for row in (col + 1)..k {
            let factor = a[row * k + col] / a[col * k + col];
            for j in col..k {
                a[row * k + j] -= factor * a[col * k + j];
            }
        }
    }
    result
}

/// Gate a node's tangent against central finite differences at a point:
/// the consistent-tangent-vs-FD obligation from the bead. Absolute
/// tolerance; entries are compared one by one and the first divergence
/// is named.
///
/// # Errors
/// [`GraphError::TangentClaimUnmet`] when no tangent is provided;
/// [`GraphError::NonFinite`] for a non-finite entry;
/// [`GraphError::IncompleteDeclaration`] (obligation names the entry)
/// when an entry diverges from the FD probe beyond `atol`.
pub fn check_consistent_tangent(
    node_id: &str,
    node: &dyn LawNode,
    state: &[f64],
    inputs: &[f64],
    atol: f64,
) -> Result<(), GraphError> {
    let declaration = node.declaration();
    let law = declaration.law.0.clone();
    let Some(tangent) = node.tangent(state, inputs) else {
        return Err(GraphError::TangentClaimUnmet {
            node: node_id.to_string(),
            law,
        });
    };
    let m = declaration.outputs.len();
    let n = declaration.inputs.len();
    if tangent.len() != m * n {
        return Err(GraphError::TangentClaimUnmet {
            node: node_id.to_string(),
            law,
        });
    }
    let h = 1e-6;
    let mut probe = inputs.to_vec();
    for j in 0..n {
        let base = probe[j];
        probe[j] = base + h;
        let plus = node.evaluate(state, &probe)?;
        probe[j] = base - h;
        let minus = node.evaluate(state, &probe)?;
        probe[j] = base;
        for i in 0..m {
            let fd = (plus.outputs[i] - minus.outputs[i]) / (2.0 * h);
            let analytic = tangent[i * n + j];
            if !analytic.is_finite() {
                return Err(GraphError::NonFinite {
                    node: node_id.to_string(),
                    field: "tangent entry",
                    bits: analytic.to_bits(),
                });
            }
            if (analytic - fd).abs() > atol {
                return Err(GraphError::IncompleteDeclaration {
                    node: node_id.to_string(),
                    law: law.clone(),
                    obligation: "consistent tangent diverges from finite-difference probe",
                });
            }
        }
    }
    Ok(())
}

/// The aggregate runtime-state schema when several laws coexist: the
/// concatenation of each node's state slots, in graph order, under one
/// schema version. The codec is exact and fail-closed — a buffer from
/// another schema version or with the wrong length refuses instead of
/// misaligning silently.
#[derive(Debug, Clone, PartialEq)]
pub struct AggregateStateSchema {
    entries: Vec<(String, usize)>,
    total: usize,
    version: u64,
}

impl AggregateStateSchema {
    /// Build from (node id, declaration) pairs in graph order. The
    /// schema version is an FNV-1a fold over node ids, slot names, and
    /// per-node state schema versions, so ANY layout change moves it.
    #[must_use]
    pub fn assemble(nodes: &[(&str, &NodeDeclaration)]) -> AggregateStateSchema {
        let mut version: u64 = 0xcbf2_9ce4_8422_2325;
        let mut fold = |bytes: &[u8]| {
            for &b in bytes {
                version ^= u64::from(b);
                version = version.wrapping_mul(0x0000_0100_0000_01b3);
            }
        };
        let mut entries = Vec::new();
        let mut total = 0;
        for (node_id, declaration) in nodes {
            fold(node_id.as_bytes());
            fold(&declaration.state_schema_version.to_le_bytes());
            for slot in &declaration.state_slots {
                fold(slot.as_bytes());
            }
            entries.push(((*node_id).to_string(), declaration.state_slots.len()));
            total += declaration.state_slots.len();
        }
        AggregateStateSchema {
            entries,
            total,
            version,
        }
    }

    /// The schema version (layout-sensitive).
    #[must_use]
    pub fn version(&self) -> u64 {
        self.version
    }

    /// Total state slots across all nodes.
    #[must_use]
    pub fn total_slots(&self) -> usize {
        self.total
    }

    /// Encode per-node state slices into one aggregate buffer, prefixed
    /// by the schema version (as exact bits in slot 0..1 is WRONG for a
    /// float buffer — the version travels beside the buffer, so encode
    /// returns both).
    ///
    /// # Errors
    /// [`GraphError::StateSchemaMismatch`] on wrong node count or any
    /// per-node arity mismatch.
    pub fn encode(&self, per_node: &[&[f64]]) -> Result<(u64, Vec<f64>), GraphError> {
        if per_node.len() != self.entries.len() {
            return Err(GraphError::StateSchemaMismatch {
                obligation: "node count differs",
                expected: self.entries.len() as u64,
                found: per_node.len() as u64,
            });
        }
        let mut buffer = Vec::with_capacity(self.total);
        for (slice, (_, len)) in per_node.iter().zip(&self.entries) {
            if slice.len() != *len {
                return Err(GraphError::StateSchemaMismatch {
                    obligation: "per-node slot count differs",
                    expected: *len as u64,
                    found: slice.len() as u64,
                });
            }
            buffer.extend_from_slice(slice);
        }
        Ok((self.version, buffer))
    }

    /// Decode an aggregate buffer back into per-node vectors, refusing
    /// version drift and length mismatch.
    ///
    /// # Errors
    /// [`GraphError::StateSchemaMismatch`] naming the failed check.
    pub fn decode(&self, version: u64, buffer: &[f64]) -> Result<Vec<Vec<f64>>, GraphError> {
        if version != self.version {
            return Err(GraphError::StateSchemaMismatch {
                obligation: "schema version differs",
                expected: self.version,
                found: version,
            });
        }
        if buffer.len() != self.total {
            return Err(GraphError::StateSchemaMismatch {
                obligation: "buffer length differs",
                expected: self.total as u64,
                found: buffer.len() as u64,
            });
        }
        let mut out = Vec::with_capacity(self.entries.len());
        let mut offset = 0;
        for (_, len) in &self.entries {
            out.push(buffer[offset..offset + len].to_vec());
            offset += len;
        }
        Ok(out)
    }
}

/// A factory building a law node from a validated card.
pub type NodeFactory = fn(&ConstitutiveModelCard) -> Result<Box<dyn LawNode>, GraphError>;

/// The implementation registry, keyed by the immutable fs-matdb
/// (LawId, LawVersion) identity. fs-material consumes the card
/// metadata — the registry validates the card, checks identity and
/// state-schema agreement, and admits the built node before returning
/// it.
#[derive(Default)]
pub struct LawRegistry {
    factories: BTreeMap<(String, u32), NodeFactory>,
}

impl LawRegistry {
    /// An empty registry.
    #[must_use]
    pub fn new() -> LawRegistry {
        LawRegistry::default()
    }

    /// Register a factory for a (law, version) identity.
    pub fn register(&mut self, law: &LawId, version: u32, factory: NodeFactory) {
        self.factories.insert((law.0.clone(), version), factory);
    }

    /// Instantiate a node from a card: validate the card, find the
    /// factory, build, verify the built declaration matches the card's
    /// identity and state schema, and admit it.
    ///
    /// # Errors
    /// [`GraphError::UnknownLaw`] for an unregistered identity; card
    /// validation refusals; [`GraphError::CardMismatch`] when the built
    /// node disagrees with the card; admission refusals.
    pub fn instantiate(
        &self,
        node_id: &str,
        card: &ConstitutiveModelCard,
    ) -> Result<Box<dyn LawNode>, GraphError> {
        card.validate()?;
        let factory = self
            .factories
            .get(&(card.law.0.clone(), card.law_version))
            .ok_or_else(|| GraphError::UnknownLaw {
                law: card.law.0.clone(),
                version: card.law_version,
            })?;
        let node = factory(card)?;
        let declaration = node.declaration();
        if declaration.law != card.law || declaration.law_version != card.law_version {
            return Err(GraphError::CardMismatch {
                law: card.law.0.clone(),
                obligation: "built node declares a different (law, version) identity",
            });
        }
        if declaration.state_schema_version != card.state_schema_version {
            return Err(GraphError::CardMismatch {
                law: card.law.0.clone(),
                obligation: "built node's state schema version differs from the card's",
            });
        }
        admit_node(node_id, node.as_ref())?;
        Ok(node)
    }
}
