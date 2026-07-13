//! The net-new crates the addendum introduced (plan addendum, Proposal 7's
//! contract discipline + the crate-contracts governance bead). Each addendum
//! crate must ship a `CONTRACT.md` with no-claim boundaries before it becomes
//! a dependency target (AGENTS.md); this module is the canonical registry of
//! those crates + their key no-claim boundary, so a governance review can see
//! the full inventory at a glance. Actual `CONTRACT.md` PRESENCE is enforced
//! by `xtask check-contracts`; this registry tracks the obligations.

/// One net-new addendum crate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AddendumCrate {
    /// The crate name.
    pub name: &'static str,
    /// What it does.
    pub purpose: &'static str,
    /// The addendum proposal it realizes (or `"governance"`).
    pub owning_proposal: &'static str,
    /// Its architecture layer.
    pub layer: &'static str,
    /// Its key no-claim boundary (the honest edge of what it does NOT do).
    pub no_claim: &'static str,
}

/// The net-new addendum crates.
const ADDENDUM_CRATES: [AddendumCrate; 7] = [
    AddendumCrate {
        name: "fs-iface",
        purpose: "coupling-graph static checker over the FEEC periodic-table type lattice",
        owning_proposal: "13",
        layer: "L3",
        no_claim: "the inf-sup pairing registry is a literature table; it does not require the element families to be built",
    },
    AddendumCrate {
        name: "fs-ladder",
        purpose: "fidelity-ladder registry: per-kernel rungs + typed prolongation/restriction",
        owning_proposal: "3",
        layer: "L3",
        no_claim: "owns rung declarations + adjacency; numerical transfers are consumer-supplied (Refine1d is a demonstrator)",
    },
    AddendumCrate {
        name: "fs-probe",
        purpose: "adjacent-rung discrepancy probes + the budget pie",
        owning_proposal: "3",
        layer: "L3",
        no_claim: "consumes fs-ladder + fs-evidence, does not run solves; the discrepancy is estimated-color, not a certified bound",
    },
    AddendumCrate {
        name: "fs-spececo",
        purpose: "certified-speculation accept/reject economics + telemetry + drift",
        owning_proposal: "9",
        layer: "L6",
        no_claim: "owns decision + telemetry + drift; solve-node schema persistence is fs-ledger's; the speculative race is the executor's",
    },
    AddendumCrate {
        name: "fs-verify",
        purpose: "certified-speculation verifier: equilibrated-flux a-posteriori error bounds, interval-evaluated",
        owning_proposal: "9",
        layer: "L3",
        no_claim: "see fs-verify/CONTRACT.md for its stated boundaries",
    },
    AddendumCrate {
        name: "fs-recompute",
        purpose: "tolerance-aware incremental-recompute store: content-addressed Merkle DAG with first-class slack",
        owning_proposal: "2",
        layer: "L6",
        no_claim: "see fs-recompute/CONTRACT.md for its stated boundaries",
    },
    AddendumCrate {
        name: "fs-govern",
        purpose: "governance as data: doctrine (P1-P8 + rules), the 19 proposals, the risk register, and this crate inventory",
        owning_proposal: "governance",
        layer: "UTIL",
        no_claim: "encodes governance data; does not measure metrics, read the beads DB, or check the filesystem",
    },
];

/// The net-new addendum crates.
#[must_use]
pub fn addendum_crates() -> &'static [AddendumCrate] {
    &ADDENDUM_CRATES
}

/// The result of auditing the crate registry for governance completeness.
#[derive(Debug, Clone, PartialEq)]
pub struct CrateAudit {
    /// Total crates.
    pub total: usize,
    /// Crates with a purpose, an owning proposal, and a no-claim boundary.
    pub complete: usize,
    /// `(crate name, reason)` for every incomplete entry.
    pub gaps: Vec<(&'static str, &'static str)>,
}

impl CrateAudit {
    /// Does every crate declare a purpose, an owner, and a no-claim boundary?
    /// Fails closed on an empty scope and requires the complete count to equal
    /// the nonzero total, mirroring [`crate::RiskAudit`]: a zero-row audit is a
    /// coverage gap, never a green.
    #[must_use]
    pub fn ok(&self) -> bool {
        self.total > 0 && self.complete == self.total && self.gaps.is_empty()
    }
}

/// Audit the crate registry: every addendum crate must declare a purpose, an
/// owning proposal, and a no-claim boundary (the AGENTS.md contract discipline
/// made governance-legible). `CONTRACT.md` file PRESENCE is enforced
/// separately by `xtask check-contracts`.
#[must_use]
pub fn crate_audit() -> CrateAudit {
    let mut gaps = Vec::new();
    let mut complete = 0usize;
    for c in &ADDENDUM_CRATES {
        let mut ok = true;
        if c.purpose.trim().is_empty() {
            gaps.push((c.name, "missing purpose"));
            ok = false;
        }
        if c.owning_proposal.trim().is_empty() {
            gaps.push((c.name, "missing owning proposal"));
            ok = false;
        }
        if c.no_claim.trim().is_empty() {
            gaps.push((c.name, "missing no-claim boundary"));
            ok = false;
        }
        if ok {
            complete += 1;
        }
    }
    CrateAudit {
        total: ADDENDUM_CRATES.len(),
        complete,
        gaps,
    }
}

/// Emit the crate registry as a deterministic machine-readable JSON array.
#[must_use]
pub fn crates_json() -> String {
    use crate::json_escape;
    use core::fmt::Write as _;
    let mut out = String::from("[");
    for (i, c) in ADDENDUM_CRATES.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        write!(
            out,
            "{{\"name\":\"{}\",\"purpose\":\"{}\",\"owning_proposal\":\"{}\",\"layer\":\"{}\",\"no_claim\":\"{}\"}}",
            json_escape(c.name),
            json_escape(c.purpose),
            json_escape(c.owning_proposal),
            json_escape(c.layer),
            json_escape(c.no_claim),
        )
        .expect("writing to a String is infallible");
    }
    out.push(']');
    out
}
