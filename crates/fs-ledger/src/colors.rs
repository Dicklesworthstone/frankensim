//! WRITE-TIME enforcement of the three-color schema (Proposal 3,
//! bead qmao.1): the [`ColorGraph`] accepts only writes whose claimed
//! color is consistent with what the composition algebra derives from
//! the parents — an estimated result CANNOT be written as verified
//! (the laundering refusal), validated claims are re-checked against
//! the CURRENT execution state and AUTO-DEMOTE on regime exit, and the
//! only override is a SIGNED WAIVER that participates in the node's
//! provenance hash (it cannot be quietly dropped later).
//!
//! The color enum and pairwise algebra live in fs-evidence (usable by
//! every layer); this module is the HELM-side gatekeeper over
//! already-colored values. Rows are canonical JSON lines ready for the
//! event stream; a dedicated schema table is a CONTRACT no-claim.

use crate::hash::{ContentHash, hash_bytes};
use fs_evidence::{Color, ColorRank, Demotion, IntervalOp, check_regime, compose};
use std::collections::BTreeMap;

/// A signed waiver: the only path past a laundering refusal, and it
/// travels IN the provenance hash.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Waiver {
    /// Waiver identifier (ticket, memo).
    pub id: String,
    /// The human who signed it.
    pub signer: String,
    /// Why.
    pub reason: String,
}

fn json_f64(value: f64) -> String {
    if value.is_finite() {
        value.to_string()
    } else {
        format!("\"non-finite:{value}\"")
    }
}

fn json_string(value: &str) -> String {
    use core::fmt::Write as _;
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if u32::from(c) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// One colored ledger node.
#[derive(Debug, Clone)]
pub struct ColorNode {
    /// Node id (write order).
    pub id: u64,
    /// Human name.
    pub name: String,
    /// The color as WRITTEN (post demotion, post waiver).
    pub color: Color,
    /// Parent node ids.
    pub parents: Vec<u64>,
    /// Demotion flag, when the regime check fired.
    pub demotion: Option<Demotion>,
    /// The waiver, when one authorized an upgrade.
    pub waiver: Option<Waiver>,
    /// Provenance hash (name, payload, parent hashes, waiver).
    pub hash: ContentHash,
}

/// Teaching errors at the write gate.
#[derive(Debug, Clone, PartialEq)]
pub enum ColorWriteError {
    /// The claimed color outranks what the parents support.
    LaunderingRefused {
        /// The claimed rank.
        claimed: ColorRank,
        /// The rank the composition algebra derived.
        derived: ColorRank,
        /// The parents that cap the rank.
        offending_parents: Vec<u64>,
    },
    /// A referenced parent does not exist.
    UnknownParent {
        /// The offending id.
        id: u64,
    },
    /// Derivations need at least one parent.
    NoParents,
}

impl core::fmt::Display for ColorWriteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ColorWriteError::LaunderingRefused {
                claimed,
                derived,
                offending_parents,
            } => write!(
                f,
                "laundering refused: the write claims {claimed:?} but the parents \
                 support at most {derived:?} (capped by nodes {offending_parents:?}); \
                 estimates cannot become certificates by assertion — attach a signed \
                 waiver if a human accepts responsibility, and it will travel in \
                 provenance"
            ),
            ColorWriteError::UnknownParent { id } => {
                write!(f, "parent node {id} does not exist in this color graph")
            }
            ColorWriteError::NoParents => {
                write!(f, "derived nodes need parents; use `source` for leaves")
            }
        }
    }
}

impl std::error::Error for ColorWriteError {}

/// The write-time color gatekeeper (append-only, deterministic).
#[derive(Debug, Default)]
pub struct ColorGraph {
    nodes: Vec<ColorNode>,
    rows: Vec<String>,
}

impl ColorGraph {
    /// Empty graph.
    #[must_use]
    pub fn new() -> Self {
        ColorGraph::default()
    }

    /// The nodes written so far.
    #[must_use]
    pub fn nodes(&self) -> &[ColorNode] {
        &self.nodes
    }

    /// The canonical JSON rows (one per write, plus demotion events).
    #[must_use]
    pub fn rows(&self) -> &[String] {
        &self.rows
    }

    fn node_hash(
        &self,
        name: &str,
        color: &Color,
        parents: &[u64],
        waiver: Option<&Waiver>,
    ) -> ContentHash {
        let mut buf = String::new();
        buf.push_str(name);
        buf.push('\n');
        buf.push_str(color.name());
        buf.push('\n');
        buf.push_str(&color.payload_json());
        buf.push('\n');
        for &p in parents {
            buf.push_str(&self.nodes[p as usize].hash.to_hex());
            buf.push('\n');
        }
        if let Some(w) = waiver {
            use core::fmt::Write as _;
            let _ = writeln!(buf, "waiver:{}:{}:{}", w.id, w.signer, w.reason);
        }
        hash_bytes(buf.as_bytes())
    }

    fn push_node(
        &mut self,
        name: &str,
        color: Color,
        parents: Vec<u64>,
        demotion: Option<Demotion>,
        waiver: Option<Waiver>,
    ) -> u64 {
        let id = self.nodes.len() as u64;
        let hash = self.node_hash(name, &color, &parents, waiver.as_ref());
        if let Some(d) = &demotion {
            self.rows.push(format!(
                "{{\"event\":\"demotion\",\"node\":{id},\"dataset\":{},\
                 \"axis\":{},\"value\":{}}}",
                json_string(&d.dataset),
                json_string(&d.axis),
                json_f64(d.value)
            ));
        }
        let waiver_json = waiver.as_ref().map_or("null".to_string(), |w| {
            format!(
                "{{\"id\":{},\"signer\":{},\"reason\":{}}}",
                json_string(&w.id),
                json_string(&w.signer),
                json_string(&w.reason)
            )
        });
        self.rows.push(format!(
            "{{\"event\":\"color-write\",\"node\":{id},\"name\":{},\
             \"color\":\"{}\",\"payload\":{},\"parents\":{:?},\"waiver\":{},\
             \"hash\":\"{}\"}}",
            json_string(name),
            color.name(),
            color.payload_json(),
            parents,
            waiver_json,
            hash.to_hex()
        ));
        self.nodes.push(ColorNode {
            id,
            name: name.to_string(),
            color,
            parents,
            demotion,
            waiver,
            hash,
        });
        id
    }

    /// Write a colored LEAF (a measurement, a certified input, an
    /// estimator output). Leaves state their color; derivations must
    /// EARN theirs.
    pub fn source(&mut self, name: &str, color: Color) -> u64 {
        self.push_node(name, color, Vec::new(), None, None)
    }

    /// Write a DERIVED node: the composition algebra folds the parent
    /// colors (with regime re-checks against `state`, auto-demoting on
    /// exit), and the claimed color must not outrank the derivation —
    /// unless a signed waiver authorizes it, in which case the waiver
    /// travels in the provenance hash.
    ///
    /// # Errors
    /// [`ColorWriteError`] teaching errors; the laundering refusal
    /// names the capping parents.
    pub fn derive(
        &mut self,
        name: &str,
        parents: &[u64],
        op: IntervalOp,
        claimed: Option<Color>,
        state: &BTreeMap<String, f64>,
        waiver: Option<Waiver>,
    ) -> Result<u64, ColorWriteError> {
        if parents.is_empty() {
            return Err(ColorWriteError::NoParents);
        }
        for &p in parents {
            if p as usize >= self.nodes.len() {
                return Err(ColorWriteError::UnknownParent { id: p });
            }
        }
        // Regime re-checks per parent (validated is REGIONAL).
        let mut demotion = None;
        let mut effective: Vec<Color> = Vec::with_capacity(parents.len());
        for &p in parents {
            let (c, d) = check_regime(&self.nodes[p as usize].color, state);
            if demotion.is_none() {
                demotion = d;
            }
            effective.push(c);
        }
        // Fold the composition algebra.
        let mut derived = effective[0].clone();
        for c in &effective[1..] {
            derived = compose(&derived, c, op);
        }
        // The claim gate.
        let written = match claimed {
            None => derived,
            Some(c) if c.rank() <= derived.rank() => c,
            Some(c) => {
                if waiver.is_none() {
                    let cap = derived.rank();
                    let offending: Vec<u64> = parents
                        .iter()
                        .copied()
                        .filter(|&p| {
                            let (eff, _) = check_regime(&self.nodes[p as usize].color, state);
                            eff.rank() <= cap
                        })
                        .collect();
                    return Err(ColorWriteError::LaunderingRefused {
                        claimed: c.rank(),
                        derived: cap,
                        offending_parents: offending,
                    });
                }
                c // waived upgrade: recorded WITH the waiver in provenance
            }
        };
        Ok(self.push_node(name, written, parents.to_vec(), demotion, waiver))
    }

    /// The node by id.
    #[must_use]
    pub fn node(&self, id: u64) -> &ColorNode {
        &self.nodes[id as usize]
    }
}
