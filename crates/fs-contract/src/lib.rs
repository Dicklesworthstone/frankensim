//! fs-contract — assume-guarantee component contracts (plan addendum,
//! Proposal E, the ambitious half of compounding swarm memory). Layer: L3.
//!
//! A contract is `(environment envelope ⇒ guaranteed behavior, certificate)`:
//! a certified design motif ("this stiffener topology dominates in this load
//! regime, envelope attached") that a thousand agent instances can REUSE
//! without re-derivation. Contracts compose into SYSTEM claims by the
//! deliberately-primitive v1 rule — **envelope containment**: component
//! contracts compose iff each component's operating conditions, as computed by
//! the system model, provably land INSIDE its envelope (an interval-box
//! containment check, hence cheap and verified-color).
//!
//! The load-bearing soundness invariant (the Gauntlet contract-composition
//! test): a composed system claim is NEVER tighter than its weakest member's
//! certificate permits — so the composed certificate takes the WEAKEST
//! member's epistemic color ([`fs_evidence::ColorRank`]). Contracts are not
//! exempt from the type system: a NONLINEAR contract cannot carry a
//! verified-color certificate until validated per regime.
//!
//! Envelope quantities are typed by their interface function-space role
//! ([`fs_iface::SpaceType`]). Everything here is pure and deterministic.

/// Typed floating-design to deployed-target refinement contracts.
pub mod deployment;

use fs_evidence::{Color, ColorRank};
use fs_iface::SpaceType;
use std::collections::BTreeMap;

/// A closed interval `[lo, hi]` (finite, `lo <= hi`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    /// Lower bound.
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
}

impl Interval {
    /// A finite, ordered interval.
    ///
    /// # Errors
    /// [`ContractError::BadInterval`] if `lo`/`hi` are non-finite or `lo > hi`.
    pub fn new(lo: f64, hi: f64) -> Result<Interval, ContractError> {
        if lo.is_finite() && hi.is_finite() && lo <= hi {
            Ok(Interval { lo, hi })
        } else {
            Err(ContractError::BadInterval { lo, hi })
        }
    }

    /// Does this (closed) interval CONTAIN `other`? Inclusive on the boundary;
    /// a point just outside is not contained.
    #[must_use]
    pub fn contains(&self, other: &Interval) -> bool {
        self.lo <= other.lo && other.hi <= self.hi
    }
}

/// An operating envelope: the interval box over named, typed interface
/// quantities within which a contract's guarantee holds.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Envelope {
    bounds: BTreeMap<String, (SpaceType, Interval)>,
}

impl Envelope {
    /// An empty envelope (no constraints).
    #[must_use]
    pub fn new() -> Envelope {
        Envelope {
            bounds: BTreeMap::new(),
        }
    }

    /// Constrain a typed quantity to an interval (builder style).
    #[must_use]
    pub fn with(mut self, quantity: &str, space: SpaceType, interval: Interval) -> Envelope {
        self.bounds.insert(quantity.to_string(), (space, interval));
        self
    }

    /// The constrained quantity names, sorted.
    #[must_use]
    pub fn quantities(&self) -> Vec<&str> {
        self.bounds.keys().map(String::as_str).collect()
    }
}

/// The system model's computed operating conditions: an interval per quantity.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct OperatingConditions {
    conds: BTreeMap<String, Interval>,
}

impl OperatingConditions {
    /// Empty conditions.
    #[must_use]
    pub fn new() -> OperatingConditions {
        OperatingConditions {
            conds: BTreeMap::new(),
        }
    }

    /// Record a quantity's operating interval (builder style).
    #[must_use]
    pub fn with(mut self, quantity: &str, interval: Interval) -> OperatingConditions {
        self.conds.insert(quantity.to_string(), interval);
        self
    }
}

/// An assume-guarantee contract.
#[derive(Debug, Clone, PartialEq)]
pub struct Contract {
    /// A stable name (referenced by composition + the library).
    pub name: String,
    /// The interface function-space this contract governs.
    pub interface: SpaceType,
    /// Is the contract's regime LINEAR? A nonlinear contract cannot carry a
    /// verified-color certificate until validated per regime.
    pub linear: bool,
    /// The operating envelope.
    pub envelope: Envelope,
    /// The guaranteed behavior, stated.
    pub guarantee: String,
    /// The certificate backing the guarantee.
    pub certificate: Color,
    /// Names of sub-contracts this contract composes over.
    pub requires: Vec<String>,
}

impl Contract {
    /// Is the certificate color legal for this contract's regime? A nonlinear
    /// contract may NOT be verified-color (Proposal 3's rule applies —
    /// contracts are not exempt from the type system).
    #[must_use]
    pub fn color_ok(&self) -> bool {
        self.linear || self.certificate.rank() != ColorRank::Verified
    }
}

/// A certified-motif catalog: the swarm's reusable contract vocabulary.
#[derive(Debug, Clone, Default)]
pub struct ContractLibrary {
    contracts: BTreeMap<String, Contract>,
}

impl ContractLibrary {
    /// An empty library.
    #[must_use]
    pub fn new() -> ContractLibrary {
        ContractLibrary {
            contracts: BTreeMap::new(),
        }
    }

    /// Add (or replace) a contract.
    pub fn insert(&mut self, contract: Contract) {
        self.contracts.insert(contract.name.clone(), contract);
    }

    /// A contract by name.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&Contract> {
        self.contracts.get(name)
    }
}

/// A composed system claim.
#[derive(Debug, Clone, PartialEq)]
pub struct SystemClaim {
    /// The member contracts used (root + transitive requires), sorted.
    pub members: Vec<String>,
    /// The composed certificate — the WEAKEST member's color (never tighter).
    pub certificate: Color,
    /// A summary of the composed guarantee.
    pub guarantee: String,
}

/// A structured contract-composition failure.
#[derive(Debug, Clone, PartialEq)]
pub enum ContractError {
    /// A non-finite or inverted interval.
    BadInterval {
        /// The offending lower bound.
        lo: f64,
        /// The offending upper bound.
        hi: f64,
    },
    /// A referenced contract is not in the library.
    UnknownContract {
        /// The missing name.
        name: String,
    },
    /// The operating conditions do not constrain a quantity the envelope needs.
    MissingCondition {
        /// The contract.
        contract: String,
        /// The quantity with no operating condition.
        quantity: String,
    },
    /// A quantity's operating conditions fall OUTSIDE the contract's envelope.
    OutsideEnvelope {
        /// The contract.
        contract: String,
        /// The out-of-envelope quantity.
        quantity: String,
    },
    /// A nonlinear contract carries a verified-color certificate (illegal).
    ColorDiscipline {
        /// The offending contract.
        contract: String,
    },
    /// The requires-graph has a cycle.
    CircularDependency {
        /// The contract at which the cycle was detected.
        at: String,
    },
}

/// Compose the system claim for `root`, resolving its transitive `requires`.
/// Every member's operating conditions (from `ops`) must land inside its
/// envelope; every member must satisfy color discipline; the requires-graph
/// must be acyclic. The composed certificate is the WEAKEST member's color.
///
/// # Errors
/// See [`ContractError`].
pub fn compose(
    lib: &ContractLibrary,
    root: &str,
    ops: &OperatingConditions,
) -> Result<SystemClaim, ContractError> {
    let mut members: Vec<String> = Vec::new();
    let mut on_stack: Vec<String> = Vec::new();
    resolve(lib, root, &mut members, &mut on_stack)?;
    members.sort();
    members.dedup();

    // weakest certificate starts at the strongest and is lowered to any weaker
    // member; the composed claim can never outrank its weakest member.
    let mut weakest: Option<Color> = None;
    for name in &members {
        let c = lib.get(name).expect("resolved members exist");
        if !c.color_ok() {
            return Err(ContractError::ColorDiscipline {
                contract: name.clone(),
            });
        }
        check_envelope(c, ops)?;
        weakest = Some(match weakest {
            None => c.certificate.clone(),
            Some(w) if c.certificate.rank() < w.rank() => c.certificate.clone(),
            Some(w) => w,
        });
    }
    let certificate = weakest.unwrap_or(Color::Estimated {
        estimator: "empty-composition".to_string(),
        dispersion: f64::INFINITY,
    });
    Ok(SystemClaim {
        guarantee: format!("system claim over {} contract(s)", members.len()),
        members,
        certificate,
    })
}

/// DFS-resolve `name` and its requires into `members`, detecting cycles via
/// the `on_stack` path.
fn resolve(
    lib: &ContractLibrary,
    name: &str,
    members: &mut Vec<String>,
    on_stack: &mut Vec<String>,
) -> Result<(), ContractError> {
    let c = lib
        .get(name)
        .ok_or_else(|| ContractError::UnknownContract {
            name: name.to_string(),
        })?;
    if on_stack.iter().any(|n| n == name) {
        return Err(ContractError::CircularDependency {
            at: name.to_string(),
        });
    }
    if members.iter().any(|n| n == name) {
        return Ok(()); // already resolved (shared sub-contract, not a cycle)
    }
    on_stack.push(name.to_string());
    for req in &c.requires {
        resolve(lib, req, members, on_stack)?;
    }
    on_stack.pop();
    members.push(name.to_string());
    Ok(())
}

/// Every envelope quantity must have an operating condition contained in its
/// interval.
fn check_envelope(c: &Contract, ops: &OperatingConditions) -> Result<(), ContractError> {
    for (quantity, (_space, allowed)) in &c.envelope.bounds {
        match ops.conds.get(quantity) {
            None => {
                return Err(ContractError::MissingCondition {
                    contract: c.name.clone(),
                    quantity: quantity.clone(),
                });
            }
            Some(op) if !allowed.contains(op) => {
                return Err(ContractError::OutsideEnvelope {
                    contract: c.name.clone(),
                    quantity: quantity.clone(),
                });
            }
            Some(_) => {}
        }
    }
    Ok(())
}
