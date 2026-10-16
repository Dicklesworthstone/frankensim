//! Cross-crate IR contracts (bead frankensim-epic-gauntlet-6nb.8,
//! slice 3): integration points between crates expressed at the IR level
//! and checked WITHOUT compiling either crate.
//!
//! A provider crate declares what its operator's output CERTIFIES (named
//! certified queries: "signed-distance", "gradient", …); a consumer crate
//! declares what one of its operator's arguments REQUIRES. Given only an
//! IR program, [`check_program`] lowers it through the real path, resolves
//! the let-binding dataflow to find which operator actually produces each
//! consumed value, and verifies requirement ⊆ certification — so
//! "fs-lbm consumes an SDF chart satisfying these certified queries" is a
//! machine-checked statement about programs, negotiated against the
//! contract rather than against either implementation. Where the operator
//! catalog knows both heads, the typed signature (produces/consumes
//! [`SigType`]) is cross-checked too; where it cannot be positionally
//! resolved, the check says so instead of pretending.
//!
//! Fail-closed: an unresolvable producer, a producer with no declared
//! contract, or a missing certified query each refuse with the exact gap
//! named — an integration seam is never silently assumed sound.
//!
//! No-claims: contracts govern PROGRAM shape; whether the provider's
//! implementation actually honors its certified queries is that crate's
//! own conformance suite's burden (this module is the negotiation layer,
//! not a proof of the physics).

use std::collections::{BTreeMap, BTreeSet};

use crate::ast::{Node, NodeKind};
use crate::catalog::{Catalog, CatalogQuery};
use crate::lower;

/// What a provider operator's output certifies.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProviderContract {
    /// The producing operator head (catalog name).
    pub head: String,
    /// The named certified queries the output supports.
    pub certifies: BTreeSet<String>,
}

/// What a consumer operator requires of one argument.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConsumerContract {
    /// The consuming operator head (catalog name).
    pub head: String,
    /// Positional index into the form's arguments (0 = first argument
    /// after the head).
    pub arg_index: usize,
    /// The certified queries the argument must support.
    pub requires: BTreeSet<String>,
}

/// The registered contract surface for one check run.
#[derive(Debug, Clone, Default)]
pub struct ContractRegistry {
    /// Provider contracts by operator head.
    pub providers: BTreeMap<String, ProviderContract>,
    /// Consumer contracts by operator head.
    pub consumers: BTreeMap<String, ConsumerContract>,
}

/// One checked integration seam.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeamCheck {
    /// The consuming operator head.
    pub consumer: String,
    /// The resolved producing operator head, when resolution succeeded.
    pub provider: Option<String>,
    /// The verdict.
    pub verdict: SeamVerdict,
    /// Whether the catalog's typed signature was cross-checked
    /// ("checked", "kind-mismatch", or a reason it could not be).
    pub signature_note: String,
}

/// The seam verdict.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum SeamVerdict {
    /// Every required certified query is provided.
    Satisfied {
        /// The queries that discharged the requirement.
        queries: BTreeSet<String>,
    },
    /// The producer could not be resolved from the program dataflow.
    UnresolvedProducer {
        /// Why.
        reason: String,
    },
    /// The producer has no declared provider contract (fail closed).
    UndeclaredProvider,
    /// Declared but insufficient: the named queries are missing.
    MissingQueries {
        /// The exact gap.
        missing: BTreeSet<String>,
    },
}

/// The whole-program contract report.
#[derive(Debug, Clone)]
pub struct ContractReport {
    /// Every consumer application found, in deterministic traversal
    /// order.
    pub seams: Vec<SeamCheck>,
}

impl ContractReport {
    /// Whether every seam is satisfied AND at least one seam was
    /// checked (a program with no seams proves nothing about them).
    #[must_use]
    pub fn clean(&self) -> bool {
        !self.seams.is_empty()
            && self
                .seams
                .iter()
                .all(|s| matches!(s.verdict, SeamVerdict::Satisfied { .. }))
    }
}

const MAX_BINDING_HOPS: usize = 64;

fn binding_env(root: &Node) -> BTreeMap<String, Node> {
    // Iterative walk (explicit-stack doctrine): collect (let name expr).
    let mut env = BTreeMap::new();
    let mut stack = vec![root];
    while let Some(node) = stack.pop() {
        if let Some(items) = node.items() {
            if node.head() == Some("let")
                && items.len() >= 3
                && let NodeKind::Symbol(name) = &items[1].kind
            {
                env.insert(name.clone(), items[2].clone());
            }
            for child in items {
                stack.push(child);
            }
        }
    }
    env
}

fn producing_head(argument: &Node, env: &BTreeMap<String, Node>) -> Result<String, String> {
    let mut current = argument.clone();
    for _ in 0..MAX_BINDING_HOPS {
        match &current.kind {
            NodeKind::Symbol(name) => match env.get(name) {
                Some(bound) => current = bound.clone(),
                None => {
                    return Err(format!("unbound symbol {name:?} — no producer in scope"));
                }
            },
            NodeKind::List(_) => {
                return current
                    .head()
                    .map(str::to_owned)
                    .ok_or_else(|| "headless form produces the consumed value".to_owned());
            }
            other => {
                return Err(format!("literal {other:?} cannot certify queries"));
            }
        }
    }
    Err("binding chain exceeded the hop cap".to_owned())
}

fn signature_note(catalog: &Catalog, provider: &str, consumer: &str, arg_index: usize) -> String {
    let find = |name: &str| {
        let query = CatalogQuery {
            name: Some(name.to_owned()),
            ..CatalogQuery::default()
        };
        catalog.query(&query).into_iter().next().cloned()
    };
    match (find(provider), find(consumer)) {
        (Some(p), Some(c)) => {
            if let Some(expected) = c.consumes.get(arg_index) {
                if *expected == p.produces {
                    format!("checked: {} feeds {}", p.produces.name(), expected.name())
                } else {
                    format!(
                        "kind-mismatch: provider produces {}, consumer expects {}",
                        p.produces.name(),
                        expected.name()
                    )
                }
            } else {
                "unchecked: argument index beyond the catalog's positional signature".to_owned()
            }
        }
        _ => "unchecked: one head is not in the operator catalog".to_owned(),
    }
}

/// Check every registered integration seam in an IR program: lower
/// through the real path, resolve the dataflow, and verify certified
/// queries — without compiling either crate.
///
/// # Errors
/// A string diagnostic when the program itself refuses to parse/lower
/// (contract checking needs a well-formed program to talk about).
pub fn check_program(
    program_sexpr: &str,
    registry: &ContractRegistry,
    catalog: &Catalog,
) -> Result<ContractReport, String> {
    let node = crate::sexpr::parse(program_sexpr).map_err(|e| e.to_string())?;
    let lowered = lower::lower(&node).map_err(|e| e.to_string())?;
    let env = binding_env(&lowered.node);

    let mut seams = Vec::new();
    // Deterministic traversal: depth-first, children in order.
    let mut stack = vec![&lowered.node];
    let mut ordered = Vec::new();
    while let Some(current) = stack.pop() {
        ordered.push(current);
        if let Some(items) = current.items() {
            for child in items.iter().rev() {
                stack.push(child);
            }
        }
    }
    for current in ordered {
        let Some(head) = current.head() else {
            continue;
        };
        let Some(contract) = registry.consumers.get(head) else {
            continue;
        };
        let items = current.items().unwrap_or_default();
        let argument = items.get(1 + contract.arg_index);
        let (provider, verdict) = match argument {
            None => (
                None,
                SeamVerdict::UnresolvedProducer {
                    reason: format!(
                        "argument {} absent from the {head} form",
                        contract.arg_index
                    ),
                },
            ),
            Some(argument) => match producing_head(argument, &env) {
                Err(reason) => (None, SeamVerdict::UnresolvedProducer { reason }),
                Ok(provider_head) => match registry.providers.get(&provider_head) {
                    None => (Some(provider_head), SeamVerdict::UndeclaredProvider),
                    Some(provider) => {
                        let missing: BTreeSet<String> = contract
                            .requires
                            .difference(&provider.certifies)
                            .cloned()
                            .collect();
                        let verdict = if missing.is_empty() {
                            SeamVerdict::Satisfied {
                                queries: contract.requires.clone(),
                            }
                        } else {
                            SeamVerdict::MissingQueries { missing }
                        };
                        (Some(provider_head), verdict)
                    }
                },
            },
        };
        let signature = provider.as_ref().map_or_else(
            || "unchecked: no resolved provider".to_owned(),
            |p| signature_note(catalog, p, head, contract.arg_index),
        );
        seams.push(SeamCheck {
            consumer: head.to_owned(),
            provider,
            verdict,
            signature_note: signature,
        });
    }
    Ok(ContractReport { seams })
}
