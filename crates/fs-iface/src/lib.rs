//! fs-iface — interface types + the coupling-graph static checker (plan
//! addendum, Proposal 13). Layer: L3 (FLUX-adjacent; the FEEC periodic table
//! is a FLUX concept, but this crate is pure static analysis with no
//! numerical dependencies).
//!
//! Type every interface by its FUNCTION-SPACE role so that ill-posed
//! couplings are COMPILE-TIME errors. Arnold's FEEC periodic table is
//! literally the type lattice: `H(grad)`, `H(curl)`, `H(div)`, `L²` — the de
//! Rham complex `H(grad) --grad--> H(curl) --curl--> H(div) --div--> L²`. A
//! static [`check`] over a [`CouplingGraph`] validates every coupling against
//! an inf-sup-certified [`PairingRegistry`] BEFORE anything runs. Coupling a
//! pressure (`L²`) to a displacement trace (`H(grad)`) in an illegal pairing
//! becomes a localized rejection, not a silently garbage solve — which
//! converts the largest class of silent agent-generated nonsense into
//! compile-time errors (humans catch these by experience; agents have none).
//!
//! Two coupling ROLES are distinguished:
//! - [`CouplingRole::Continuity`] — matching the SAME field across an
//!   interface (trace continuity). Legal iff both sides live in the same
//!   trace space (an `H(div)` normal trace cannot match an `H(curl)`
//!   tangential trace).
//! - [`CouplingRole::Saddle`] — a mixed / saddle-point block pairing a trial
//!   space with a test (multiplier) space. Legality is inf-sup (LBB)
//!   stability, read from the [`PairingRegistry`].
//!
//! CONSERVATIVE BY DEFAULT (the load-bearing soundness rule): a saddle pairing
//! that is neither certified nor known-unstable is [`PairingVerdict::Unknown`]
//! and is REJECTED as illegal-until-certified — a checker that silently admits
//! pairings it does not recognize is unsound in practice.
//!
//! The registry is a DECLARATIVE literature table (LBB/inf-sup results). The
//! checker guarantees coupling-graph LEGALITY against the registry; it does
//! NOT require the corresponding saddle-point element families or solvers to
//! be implemented. Determinism: couplings are processed in order and findings
//! are emitted in order, so a replayed check reproduces the same report.

use std::collections::BTreeMap;

/// The FEEC periodic table: the de Rham complex spaces, in exact-sequence
/// order. This IS the interface type lattice.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SpaceType {
    /// `H(grad) = H¹`: scalar potentials (0-forms). Full trace continuity.
    HGrad,
    /// `H(curl)`: fields with square-integrable curl (1-forms). Tangential
    /// trace continuity.
    HCurl,
    /// `H(div)`: fields with square-integrable divergence (2-forms). Normal
    /// trace continuity. Fluxes live here.
    HDiv,
    /// `L²`: square-integrable (3-forms). No trace continuity — the natural
    /// home of pressures and multipliers.
    L2,
}

impl SpaceType {
    /// Stable machine-readable name.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            SpaceType::HGrad => "H(grad)",
            SpaceType::HCurl => "H(curl)",
            SpaceType::HDiv => "H(div)",
            SpaceType::L2 => "L2",
        }
    }

    /// The differential-form degree in 3D (`H(grad)`→0 … `L²`→3): the
    /// position in the de Rham complex.
    #[must_use]
    pub fn form_degree(self) -> u8 {
        match self {
            SpaceType::HGrad => 0,
            SpaceType::HCurl => 1,
            SpaceType::HDiv => 2,
            SpaceType::L2 => 3,
        }
    }
}

/// The role a coupling plays.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CouplingRole {
    /// Trace continuity of the SAME field across an interface.
    Continuity,
    /// A mixed / saddle-point block pairing a trial and a test space.
    Saddle,
}

/// A field living on an interface, typed by its function space.
#[derive(Debug, Clone, PartialEq)]
pub struct InterfaceField {
    /// A stable identifier (referenced by couplings).
    pub id: String,
    /// The function-space role.
    pub space: SpaceType,
}

/// A declared coupling between two interface fields.
#[derive(Debug, Clone, PartialEq)]
pub struct Coupling {
    /// A stable identifier (referenced by findings).
    pub id: String,
    /// The trial field's id.
    pub trial: String,
    /// The test / multiplier field's id.
    pub test: String,
    /// What kind of coupling this is.
    pub role: CouplingRole,
}

/// A coupling graph: the interface fields and the couplings among them.
#[derive(Debug, Clone, Default)]
pub struct CouplingGraph {
    fields: BTreeMap<String, SpaceType>,
    couplings: Vec<Coupling>,
}

impl CouplingGraph {
    /// An empty graph.
    #[must_use]
    pub fn new() -> CouplingGraph {
        CouplingGraph {
            fields: BTreeMap::new(),
            couplings: Vec::new(),
        }
    }

    /// Declare an interface field (builder style).
    #[must_use]
    pub fn field(mut self, id: &str, space: SpaceType) -> CouplingGraph {
        self.fields.insert(id.to_string(), space);
        self
    }

    /// Declare a coupling (builder style).
    #[must_use]
    pub fn couple(
        mut self,
        id: &str,
        trial: &str,
        test: &str,
        role: CouplingRole,
    ) -> CouplingGraph {
        self.couplings.push(Coupling {
            id: id.to_string(),
            trial: trial.to_string(),
            test: test.to_string(),
            role,
        });
        self
    }

    /// The declared field ids, sorted (deterministic teaching output).
    #[must_use]
    pub fn field_ids(&self) -> Vec<&str> {
        self.fields.keys().map(String::as_str).collect()
    }
}

/// The inf-sup classification of a saddle pairing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairingVerdict {
    /// A known inf-sup-stable pairing, with a literature citation.
    Certified {
        /// The result it rests on.
        cite: &'static str,
    },
    /// A known inf-sup-UNSTABLE pairing, with the reason.
    Unstable {
        /// Why it is unstable.
        reason: &'static str,
    },
    /// Not in the registry — treated as illegal-until-certified.
    Unknown,
}

/// A declarative registry of inf-sup (LBB) saddle-point pairings. This is a
/// LITERATURE TABLE, not an assertion that the element families are built.
#[derive(Debug, Clone)]
pub struct PairingRegistry {
    saddle: BTreeMap<(SpaceType, SpaceType), PairingVerdict>,
}

impl PairingRegistry {
    /// An empty registry (every saddle pairing is `Unknown`).
    #[must_use]
    pub fn empty() -> PairingRegistry {
        PairingRegistry {
            saddle: BTreeMap::new(),
        }
    }

    /// The standard table of well-established inf-sup results.
    #[must_use]
    pub fn standard() -> PairingRegistry {
        let mut r = PairingRegistry::empty();
        // Certified stable mixed methods (trial space, test/multiplier space):
        r.set_certified(
            SpaceType::HDiv,
            SpaceType::L2,
            "Brezzi–Fortin: RT/BDM mixed Poisson / Darcy (H(div)–L² is LBB-stable)",
        );
        r.set_certified(
            SpaceType::HGrad,
            SpaceType::L2,
            "Taylor–Hood / Boffi–Brezzi–Fortin: H¹ velocity with L² pressure (LBB-stable)",
        );
        // Known-UNSTABLE equal-order pairings (the classic LBB violations):
        r.set_unstable(
            SpaceType::L2,
            SpaceType::L2,
            "no inf-sup control between equal L² spaces — spurious/checkerboard modes",
        );
        r.set_unstable(
            SpaceType::HGrad,
            SpaceType::HGrad,
            "equal-order H¹ velocity–pressure violates LBB (needs Taylor–Hood or stabilization)",
        );
        r
    }

    /// Register a certified pairing.
    pub fn set_certified(&mut self, trial: SpaceType, test: SpaceType, cite: &'static str) {
        self.saddle
            .insert((trial, test), PairingVerdict::Certified { cite });
    }

    /// Register a known-unstable pairing.
    pub fn set_unstable(&mut self, trial: SpaceType, test: SpaceType, reason: &'static str) {
        self.saddle
            .insert((trial, test), PairingVerdict::Unstable { reason });
    }

    /// Classify a saddle pairing (`Unknown` if not registered).
    #[must_use]
    pub fn classify_saddle(&self, trial: SpaceType, test: SpaceType) -> PairingVerdict {
        self.saddle
            .get(&(trial, test))
            .copied()
            .unwrap_or(PairingVerdict::Unknown)
    }

    /// The certified test partners for a trial space, sorted — used to teach
    /// the operator a legal alternative when their pairing is rejected.
    #[must_use]
    pub fn certified_partners(&self, trial: SpaceType) -> Vec<SpaceType> {
        self.saddle
            .iter()
            .filter(|((t, _), v)| *t == trial && matches!(v, PairingVerdict::Certified { .. }))
            .map(|((_, test), _)| *test)
            .collect()
    }
}

impl Default for PairingRegistry {
    fn default() -> Self {
        PairingRegistry::standard()
    }
}

/// Finding severity. `Reject` blocks the coupling graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    /// Advisory.
    Warn,
    /// Blocks the graph.
    Reject,
}

/// One localized finding against a coupling.
#[derive(Debug, Clone, PartialEq)]
pub struct CheckFinding {
    /// The offending coupling's id (localization).
    pub coupling: String,
    /// Which check produced it.
    pub check: &'static str,
    /// Reject or Warn.
    pub severity: Severity,
    /// The diagnosis.
    pub what: String,
    /// How to fix it.
    pub fix: String,
}

/// The verdict of checking a coupling graph.
#[derive(Debug, Clone, PartialEq)]
pub struct CheckReport {
    /// True iff no finding is a `Reject`.
    pub admitted: bool,
    /// Findings, in coupling order (deterministic).
    pub findings: Vec<CheckFinding>,
}

impl CheckReport {
    /// A one-line structured diagnosis for logging (never printed to stdout by
    /// library code — the caller decides).
    #[must_use]
    pub fn diagnosis(&self) -> String {
        let rejects = self
            .findings
            .iter()
            .filter(|f| f.severity == Severity::Reject)
            .count();
        if self.admitted {
            format!(
                "coupling graph admitted ({} coupling checks)",
                self.findings.len()
            )
        } else {
            format!("coupling graph rejected ({rejects} illegal coupling(s))")
        }
    }
}

/// Statically check a coupling graph against an inf-sup pairing registry.
/// Pure and deterministic. Every illegal or unrecognized coupling yields a
/// localized [`CheckFinding`] with a teaching fix; unknown saddle pairings are
/// rejected conservatively (illegal-until-certified). An empty graph is
/// vacuously legal.
#[must_use]
pub fn check(graph: &CouplingGraph, registry: &PairingRegistry) -> CheckReport {
    let mut findings = Vec::new();
    for c in &graph.couplings {
        if let Some(f) = check_coupling(c, graph, registry) {
            findings.push(f);
        }
    }
    let admitted = !findings.iter().any(|f| f.severity == Severity::Reject);
    CheckReport { admitted, findings }
}

/// Check one coupling, returning a finding if it is illegal.
fn check_coupling(
    c: &Coupling,
    graph: &CouplingGraph,
    registry: &PairingRegistry,
) -> Option<CheckFinding> {
    // Resolve both endpoints; a dangling reference is a malformed graph.
    let (Some(&trial_space), Some(&test_space)) =
        (graph.fields.get(&c.trial), graph.fields.get(&c.test))
    else {
        let missing = if graph.fields.contains_key(&c.trial) {
            &c.test
        } else {
            &c.trial
        };
        return Some(CheckFinding {
            coupling: c.id.clone(),
            check: "graph.field",
            severity: Severity::Reject,
            what: format!(
                "coupling '{}' references undeclared field '{missing}'",
                c.id
            ),
            fix: format!("declare '{missing}' with its function space, or fix the reference"),
        });
    };
    // A field coupled to itself is a degenerate (self-loop) coupling.
    if c.trial == c.test {
        return Some(CheckFinding {
            coupling: c.id.clone(),
            check: "graph.self-loop",
            severity: Severity::Reject,
            what: format!("coupling '{}' couples field '{}' to itself", c.id, c.trial),
            fix: "couple two distinct interface fields".to_string(),
        });
    }
    match c.role {
        CouplingRole::Continuity => check_continuity(c, trial_space, test_space),
        CouplingRole::Saddle => check_saddle(c, trial_space, test_space, registry),
    }
}

/// Continuity is legal iff both sides share a trace space.
fn check_continuity(c: &Coupling, trial: SpaceType, test: SpaceType) -> Option<CheckFinding> {
    if trial == test {
        return None;
    }
    Some(CheckFinding {
        coupling: c.id.clone(),
        check: "coupling.continuity",
        severity: Severity::Reject,
        what: format!(
            "incompatible traces: {} cannot continuity-couple to {} (different trace spaces)",
            trial.name(),
            test.name()
        ),
        fix: format!(
            "match the trace spaces ({} to {}), or use a saddle coupling with a certified pairing",
            trial.name(),
            trial.name()
        ),
    })
}

/// A saddle pairing is legal iff the registry certifies it.
fn check_saddle(
    c: &Coupling,
    trial: SpaceType,
    test: SpaceType,
    registry: &PairingRegistry,
) -> Option<CheckFinding> {
    match registry.classify_saddle(trial, test) {
        PairingVerdict::Certified { .. } => None,
        PairingVerdict::Unstable { reason } => Some(CheckFinding {
            coupling: c.id.clone(),
            check: "coupling.infsup",
            severity: Severity::Reject,
            what: format!(
                "inf-sup-UNSTABLE saddle pairing ({}, {}): {reason}",
                trial.name(),
                test.name()
            ),
            fix: certified_fix(trial, registry),
        }),
        PairingVerdict::Unknown => Some(CheckFinding {
            coupling: c.id.clone(),
            check: "coupling.infsup",
            severity: Severity::Reject,
            what: format!(
                "unrecognized saddle pairing ({}, {}): illegal until certified in the registry",
                trial.name(),
                test.name()
            ),
            fix: certified_fix(trial, registry),
        }),
    }
}

/// Teach a legal test partner for a trial space, if one exists.
fn certified_fix(trial: SpaceType, registry: &PairingRegistry) -> String {
    let partners = registry.certified_partners(trial);
    if partners.is_empty() {
        format!(
            "no certified test partner for {} in this registry — register the inf-sup result first",
            trial.name()
        )
    } else {
        let names: Vec<&str> = partners.iter().map(|s| s.name()).collect();
        format!(
            "pair {} with a certified test space: {}",
            trial.name(),
            names.join(", ")
        )
    }
}
