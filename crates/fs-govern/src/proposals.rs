//! The nineteen addendum proposals as machine-readable governance data (plan
//! addendum, Part II + Part III + Part IV). Each carries its composite score,
//! phase, KILL METRIC, and OWNING bead, and [`governance_audit`] enforces
//! Governance Rule 2 / design principle P8: a proposal with no instrumented
//! kill measurement counts as killed, so every proposal must at least DECLARE
//! a kill metric and an owner.

use crate::json_escape;
use core::fmt::Write as _;

/// One addendum proposal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Proposal {
    /// The stable id (`"9"`, `"E"`, `"D"`, …).
    pub id: &'static str,
    /// The proposal name.
    pub name: &'static str,
    /// The build phase (`spine` / `flywheel` / `leverage` / `horizon` /
    /// `commercial`).
    pub phase: &'static str,
    /// The composite Mean score (0–1000).
    pub mean: u16,
    /// The measurable kill criterion (P8).
    pub kill_metric: &'static str,
    /// The owning bead.
    pub owning_bead: &'static str,
    /// Is the kill metric actually instrumented on a dashboard yet?
    pub instrumented: bool,
}

/// The nineteen proposals, in composite (Mean) order.
const PROPOSALS: [Proposal; 19] = [
    Proposal {
        id: "9",
        name: "Certified speculation",
        phase: "flywheel",
        mean: 850,
        kill_metric: "accept-rate >30% AND median warm-start >=1.5x at customer tolerances (6-month checkpoint)",
        owning_bead: "frankensim-epic-flywheel-lmp4.1",
        instrumented: false,
    },
    Proposal {
        id: "2",
        name: "Tolerance-aware incremental recomputation",
        phase: "flywheel",
        mean: 840,
        kill_metric: "certified skip-yield >=2x median wall-clock vs plain hash memoization",
        owning_bead: "frankensim-epic-flywheel-lmp4.8",
        instrumented: false,
    },
    Proposal {
        id: "10",
        name: "Version control for physics",
        phase: "flywheel",
        mean: 820,
        kill_metric: "<25% of realistic merges surface harmonic (structural) conflicts",
        owning_bead: "frankensim-epic-flywheel-lmp4.12",
        instrumented: false,
    },
    Proposal {
        id: "8",
        name: "Declarative queries against physics",
        phase: "leverage",
        mean: 810,
        kill_metric: "greedy planner beats a fixed baseline by >=2x cost at equal certified accuracy (else ship the interface anyway)",
        owning_bead: "frankensim-epic-flywheel-lmp4.16",
        instrumented: false,
    },
    Proposal {
        id: "E",
        name: "Compounding swarm memory (contracts + tombstones)",
        phase: "spine",
        mean: 810,
        kill_metric: "re-exploration rate falls materially; envelope containment is satisfiable on real assemblies",
        owning_bead: "frankensim-epic-flywheel-lmp4.13",
        instrumented: false,
    },
    Proposal {
        id: "3",
        name: "The three-color ledger",
        phase: "spine",
        mean: 810,
        kill_metric: "probe-derived model-form maps actually change downstream decisions (probe compute capped)",
        owning_bead: "frankensim-epic-epistype-qmao.1",
        instrumented: false,
    },
    Proposal {
        id: "6",
        name: "Falsifier pairing",
        phase: "spine",
        mean: 790,
        kill_metric: "falsifier yield (true catches per compute) per class stays positive",
        owning_bead: "frankensim-epic-epistype-qmao.4",
        instrumented: false,
    },
    Proposal {
        id: "F",
        name: "Objective epistemics",
        phase: "leverage",
        mean: 790,
        kill_metric: "robust optima not consistently dominated by nominal-optimum-plus-safety-factor on realized cost",
        owning_bead: "frankensim-epic-epistype-qmao.3",
        instrumented: false,
    },
    Proposal {
        id: "A",
        name: "Certified abstraction ladder",
        phase: "horizon",
        mean: 780,
        kill_metric: "RB-certified regions cover >=20% of wedge-vertical query volume",
        owning_bead: "frankensim-epic-selfknow-knh1.4",
        instrumented: false,
    },
    Proposal {
        id: "C",
        name: "Value-of-information queries",
        phase: "horizon",
        mean: 780,
        kill_metric: "VoI-recommended purchases outperform agent-chosen at matched cost on realized decision-changes",
        owning_bead: "frankensim-epic-selfknow-knh1.6",
        instrumented: false,
    },
    Proposal {
        id: "1",
        name: "End-to-end adjoints",
        phase: "leverage",
        mean: 770,
        kill_metric: "adjoint-driven optimization beats the best derivative-free baseline on >=70% of benchmark design tasks",
        owning_bead: "frankensim-epic-coupling-bk0o.1",
        instrumented: false,
    },
    Proposal {
        id: "B",
        name: "Explanation objects",
        phase: "leverage",
        mean: 770,
        kill_metric: "attributed channels + residual reconcile to the observed change within bounds on >=90% of cases",
        owning_bead: "frankensim-epic-selfknow-knh1.5",
        instrumented: false,
    },
    Proposal {
        id: "D",
        name: "The Goodhart guard",
        phase: "spine",
        mean: 770,
        kill_metric: "guard endpoint catch-rate exceeds its catch-rate on random non-endpoint designs",
        owning_bead: "frankensim-epic-epistype-qmao.5",
        instrumented: false,
    },
    Proposal {
        id: "11",
        name: "Reality is just another chart",
        phase: "horizon",
        mean: 760,
        kill_metric: "registration uncertainty stays below the geometric deviations being certified",
        owning_bead: "frankensim-epic-coupling-bk0o.4",
        instrumented: false,
    },
    Proposal {
        id: "13",
        name: "Interface types + symmetry harvesting",
        phase: "spine",
        mean: 720,
        kill_metric: "type checker ships (no kill); symmetry: >=15% of workloads present exploitable symmetry",
        owning_bead: "frankensim-epic-selfknow-knh1.1",
        instrumented: false,
    },
    Proposal {
        id: "12",
        name: "Evidence packages",
        phase: "leverage",
        mean: 720,
        kill_metric: "an auditor/certification body engages the machine-checkable format as at least supplementary evidence",
        owning_bead: "frankensim-epic-epistype-qmao.6",
        instrumented: false,
    },
    Proposal {
        id: "5",
        name: "Spectral health monitoring",
        phase: "flywheel",
        mean: 660,
        kill_metric: "gap collapse observed outside synthetic cases at volume (else demote to sampled spot checks)",
        owning_bead: "frankensim-epic-selfknow-knh1.3",
        instrumented: false,
    },
    Proposal {
        id: "7",
        name: "The wedge and the plugin surface",
        phase: "commercial",
        mean: 640,
        kill_metric: "referenceable customer with measured cycle-time reduction >=3x within two quarters of GA",
        owning_bead: "frankensim-epic-gtm-jwq8.1",
        instrumented: false,
    },
    Proposal {
        id: "4",
        name: "Extend the complex into time",
        phase: "horizon",
        mean: 590,
        kill_metric: "a paying workload's error budget is dominated by splitting error (>=20% of budget)",
        owning_bead: "frankensim-epic-coupling-bk0o.7",
        instrumented: false,
    },
];

/// The nineteen proposals, in composite (Mean) order.
#[must_use]
pub fn proposals() -> &'static [Proposal] {
    &PROPOSALS
}

/// The result of auditing the proposals for governance completeness.
#[derive(Debug, Clone, PartialEq)]
pub struct GovernanceAudit {
    /// Total proposals.
    pub total: usize,
    /// Proposals that DECLARE both a kill metric and an owning bead.
    pub with_kill_metric_and_owner: usize,
    /// Proposals whose kill metric is actually instrumented.
    pub instrumented: usize,
    /// `(proposal id, reason)` for every incomplete entry.
    pub gaps: Vec<(&'static str, &'static str)>,
}

impl GovernanceAudit {
    /// Does every proposal declare a kill metric and an owner?
    #[must_use]
    pub fn ok(&self) -> bool {
        self.gaps.is_empty()
    }
}

/// Audit the proposals: every one must declare a non-empty kill metric AND an
/// owning bead (Governance Rule 2 — a proposal whose kill measurement was
/// never instrumented counts as killed). Also counts how many kill metrics are
/// actually instrumented.
#[must_use]
pub fn governance_audit() -> GovernanceAudit {
    let mut gaps = Vec::new();
    let mut complete = 0usize;
    let mut instrumented = 0usize;
    for p in &PROPOSALS {
        let mut ok = true;
        if p.kill_metric.trim().is_empty() {
            gaps.push((p.id, "missing kill metric"));
            ok = false;
        }
        if p.owning_bead.trim().is_empty() {
            gaps.push((p.id, "missing owning bead"));
            ok = false;
        }
        if ok {
            complete += 1;
        }
        if p.instrumented {
            instrumented += 1;
        }
    }
    GovernanceAudit {
        total: PROPOSALS.len(),
        with_kill_metric_and_owner: complete,
        instrumented,
        gaps,
    }
}

/// Emit the proposals as a machine-readable JSON array (id, name, phase, mean,
/// kill_metric, owning_bead, instrumented). Deterministic.
#[must_use]
pub fn proposals_json() -> String {
    let mut out = String::from("[");
    for (i, p) in PROPOSALS.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        write!(
            out,
            "{{\"id\":\"{}\",\"name\":\"{}\",\"phase\":\"{}\",\"mean\":{},\"kill_metric\":\"{}\",\"owning_bead\":\"{}\",\"instrumented\":{}}}",
            json_escape(p.id),
            json_escape(p.name),
            json_escape(p.phase),
            p.mean,
            json_escape(p.kill_metric),
            json_escape(p.owning_bead),
            p.instrumented,
        )
        .expect("writing to a String is infallible");
    }
    out.push(']');
    out
}
