//! The addendum's DOCTRINE (plan addendum, Part I.5 + Part IV.3): the eight
//! design principles P1–P8 and the four governance rules, as machine-readable
//! data. Every addendum bead is meant to honor these; encoding them here lets
//! a review or CI gate cite the exact principle a change respects or violates.

/// One design principle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Principle {
    /// The id (`"P1"` … `"P8"`).
    pub id: &'static str,
    /// A short name.
    pub name: &'static str,
    /// The principle, stated.
    pub statement: &'static str,
}

/// The eight design principles.
pub const PRINCIPLES: [Principle; 8] = [
    Principle {
        id: "P1",
        name: "Buy correctness at the level of algebra",
        statement: "Prefer formulations where the invariant holds identically (δ∘δ=0, discrete de Rham exactness, type-checked couplings) over formulations patched per-bug.",
    },
    Principle {
        id: "P2",
        name: "Certificates must compose",
        statement: "Any accuracy, sensitivity, or consistency statement that cannot propagate through the ledger DAG is a side channel, and side channels rot.",
    },
    Principle {
        id: "P3",
        name: "Every certificate ships its cheapest independent falsifier",
        statement: "Certificates prove what the model claims; falsifiers probe whether the claims connect to reality. The gap between them is where systems silently rot.",
    },
    Principle {
        id: "P4",
        name: "Type the epistemics",
        statement: "Verified, validated, and estimated are different kinds of knowledge with different composition rules; the type system must refuse to launder one into another.",
    },
    Principle {
        id: "P5",
        name: "Design for the swarm's deficits",
        statement: "No intuition -> certified abstraction ladders + explanation objects. No memory -> contracts + tombstones. Miscalibration -> falsification budgets. Goodhart tendency -> distrust optima by default.",
    },
    Principle {
        id: "P6",
        name: "Exploit the check/produce asymmetry",
        statement: "Verifying a candidate solution is vastly cheaper than producing one. This asymmetry is the economic engine: speculation, cheap re-verification, and third-party audit all run on it.",
    },
    Principle {
        id: "P7",
        name: "Concentrate research risk",
        statement: "One research-grade bet at a time, chosen so its hard problem unlocks the most downstream value; everything else ships as engineering over known results.",
    },
    Principle {
        id: "P8",
        name: "The plan itself must be falsifiable",
        statement: "Every proposal carries a kill criterion with a measurement and a deadline. A proposal that survives only because nobody defined its failure condition is dead weight.",
    },
];

/// The eight design principles.
#[must_use]
pub fn principles() -> &'static [Principle] {
    &PRINCIPLES
}

/// One governance rule (Part IV.3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GovernanceRule {
    /// The rule number (1–4).
    pub number: u8,
    /// A short name.
    pub name: &'static str,
    /// The rule, stated.
    pub statement: &'static str,
}

/// The four governance rules.
pub const RULES: [GovernanceRule; 4] = [
    GovernanceRule {
        number: 1,
        name: "One research bet at a time",
        statement: "Phase 1's bet is Proposal 9, chosen because its hard problem is partially purchasable and unlocks the widest downstream cone. Proposal 8's planner and Proposal A's general certification do NOT open as research programs.",
    },
    GovernanceRule {
        number: 2,
        name: "Kill criteria enforced quarterly, in writing",
        statement: "Against the named measurements. A proposal whose kill measurement was never instrumented counts as killed — unmeasured survival is not survival.",
    },
    GovernanceRule {
        number: 3,
        name: "Build in composite order, sell in Ex-A order",
        statement: "Engineering priority follows the Mean column; commercial narrative and partnership investment follow the Ex-A column (12 and 7 lead the sell-side despite mid-table composite ranks).",
    },
    GovernanceRule {
        number: 4,
        name: "The circularity discount stands until lifted by evidence",
        statement: "Proposals A, B, C, D's top-tier standing is conditional on agents being the dominant operator; if the first vertical ships human-driven, they are formally deferred to Phase 3 regardless of score.",
    },
];

/// The four governance rules.
#[must_use]
pub fn rules() -> &'static [GovernanceRule] {
    &RULES
}
