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

/// A human ANNOTATION (ticket, memo, name, rationale). It travels in
/// provenance but AUTHORIZES NOTHING (bead qmao.1.1): presence of
/// caller-created strings is not proof. The only path past a
/// laundering refusal is an authenticated [`WaiverGrant`] through
/// [`ColorGraph::derive_waived`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Waiver {
    /// Waiver identifier (ticket, memo).
    pub id: String,
    /// The human who accepts responsibility.
    pub signer: String,
    /// Why.
    pub reason: String,
}

/// The canonical scope string a color-upgrade grant must carry.
pub const WAIVER_SCOPE_COLOR_UPGRADE: &str = "color-upgrade";

/// An AUTHENTICATED waiver: a versioned, length-prefixed payload bound
/// to the exact node identity, evidence lineage, claimed color, scope,
/// signer key, and expiry — plus signature bytes over that payload.
/// Verification happens through a caller-supplied [`WaiverVerifier`]
/// capability; the grant travels whole in the provenance hash.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WaiverGrant {
    /// The human annotation riding along (never itself authorizing).
    pub annotation: Waiver,
    /// Issuer key identity the verifier resolves.
    pub key_id: String,
    /// Must equal [`WAIVER_SCOPE_COLOR_UPGRADE`] for color upgrades.
    pub scope: String,
    /// The node name this grant is bound to.
    pub node_name: String,
    /// The color name (`Color::name`) being authorized.
    pub claimed_color: String,
    /// The exact parent provenance hashes, in write order — binds the
    /// grant to one evidence lineage (replay to another node fails).
    pub parent_hashes: Vec<ContentHash>,
    /// Last day the grant is valid (days since 2026-01-01).
    pub expires_day: u32,
    /// Signature bytes over [`WaiverGrant::signing_payload`].
    pub signature: Vec<u8>,
}

impl WaiverGrant {
    /// Canonical signing payload, VERSIONED and LENGTH-PREFIXED (no
    /// delimiters, so adversarial text cannot collide structurally):
    /// version byte 1, then each field as u32-LE length + bytes, then
    /// parent count + 32-byte hashes, then expiry as u32 LE. The
    /// signature is NOT part of its own payload.
    #[must_use]
    pub fn signing_payload(&self) -> Vec<u8> {
        let mut out = vec![1u8];
        for field in [
            self.key_id.as_str(),
            self.scope.as_str(),
            self.node_name.as_str(),
            self.claimed_color.as_str(),
            self.annotation.id.as_str(),
            self.annotation.signer.as_str(),
            self.annotation.reason.as_str(),
        ] {
            out.extend_from_slice(&(field.len() as u32).to_le_bytes());
            out.extend_from_slice(field.as_bytes());
        }
        out.extend_from_slice(&(self.parent_hashes.len() as u32).to_le_bytes());
        for h in &self.parent_hashes {
            out.extend_from_slice(h.to_hex().as_bytes());
        }
        out.extend_from_slice(&self.expires_day.to_le_bytes());
        out
    }
}

/// The signature-verification CAPABILITY (injected; this crate ships
/// no cryptography). Implementations resolve `key_id` and check
/// `signature` over `payload`.
pub trait WaiverVerifier {
    /// True iff `signature` authenticates `payload` under `key_id`.
    fn verify(&self, key_id: &str, payload: &[u8], signature: &[u8]) -> bool;
}

/// The in-tree default: NO verifier exists, so NOTHING authenticates
/// (the no-crypto no-claim — fail closed until a Franken-compliant
/// signature capability is wired in).
#[derive(Debug, Default, Clone, Copy)]
pub struct NoWaiverVerifier;

impl WaiverVerifier for NoWaiverVerifier {
    fn verify(&self, _key_id: &str, _payload: &[u8], _signature: &[u8]) -> bool {
        false
    }
}

/// Why a grant failed to authorize (structured, teaching).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WaiverRejection {
    /// Scope is not [`WAIVER_SCOPE_COLOR_UPGRADE`].
    ScopeMismatch,
    /// The grant names a different node.
    NodeMismatch,
    /// The grant authorizes a different color than claimed.
    ColorMismatch,
    /// The grant's parent hashes differ from the actual lineage
    /// (replay to another node / tampered evidence).
    LineageMismatch,
    /// Expired as of the supplied day.
    Expired,
    /// The verifier refused the signature (wrong key, tampered
    /// payload, rotated-out key, or no verifier capability at all).
    BadSignature,
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
    /// The human annotation, when one was recorded (never authorizing).
    pub waiver: Option<Waiver>,
    /// The authenticated grant, when one authorized an upgrade.
    pub grant: Option<WaiverGrant>,
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
    /// A waiver grant failed authentication or binding checks; the
    /// promotion is refused (fail closed).
    WaiverRefused {
        /// The structured reason.
        rejection: WaiverRejection,
    },
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
                 estimates cannot become certificates by assertion — an authenticated \
                 WaiverGrant via derive_waived is the only path past this refusal, and \
                 it travels whole in provenance"
            ),
            ColorWriteError::UnknownParent { id } => {
                write!(f, "parent node {id} does not exist in this color graph")
            }
            ColorWriteError::NoParents => {
                write!(f, "derived nodes need parents; use `source` for leaves")
            }
            ColorWriteError::WaiverRefused { rejection } => write!(
                f,
                "waiver refused ({rejection:?}): promotion requires an authenticated \
                 grant bound to this node, lineage, color, and scope, unexpired, with \
                 a signature the verifier capability accepts — fail closed otherwise"
            ),
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

    /// Provenance hash over a VERSIONED, LENGTH-PREFIXED encoding
    /// (bead qmao.1.1): the former newline/colon-delimited encoding let
    /// adversarial text collide structurally (a name containing a
    /// newline could impersonate a parent hash line).
    fn node_hash(
        &self,
        name: &str,
        color: &Color,
        parents: &[u64],
        waiver: Option<&Waiver>,
        grant: Option<&WaiverGrant>,
    ) -> ContentHash {
        let mut buf = vec![2u8]; // encoding version
        let mut field = |b: &mut Vec<u8>, s: &str| {
            b.extend_from_slice(&(s.len() as u32).to_le_bytes());
            b.extend_from_slice(s.as_bytes());
        };
        field(&mut buf, name);
        field(&mut buf, color.name());
        field(&mut buf, &color.payload_json());
        buf.extend_from_slice(&(parents.len() as u32).to_le_bytes());
        for &p in parents {
            let hex = self.nodes[p as usize].hash.to_hex();
            field(&mut buf, &hex);
        }
        match waiver {
            Some(w) => {
                buf.push(1);
                field(&mut buf, &w.id);
                field(&mut buf, &w.signer);
                field(&mut buf, &w.reason);
            }
            None => buf.push(0),
        }
        match grant {
            Some(g) => {
                buf.push(1);
                let payload = g.signing_payload();
                buf.extend_from_slice(&(payload.len() as u32).to_le_bytes());
                buf.extend_from_slice(&payload);
                buf.extend_from_slice(&(g.signature.len() as u32).to_le_bytes());
                buf.extend_from_slice(&g.signature);
            }
            None => buf.push(0),
        }
        hash_bytes(&buf)
    }

    #[allow(clippy::needless_pass_by_value)]
    fn push_node(
        &mut self,
        name: &str,
        color: Color,
        parents: Vec<u64>,
        demotion: Option<Demotion>,
        waiver: Option<Waiver>,
        grant: Option<WaiverGrant>,
    ) -> u64 {
        let id = self.nodes.len() as u64;
        let hash = self.node_hash(name, &color, &parents, waiver.as_ref(), grant.as_ref());
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
        let grant_json = grant.as_ref().map_or("null".to_string(), |g| {
            format!(
                "{{\"key_id\":{},\"scope\":{},\"expires_day\":{},\"authorized\":true}}",
                json_string(&g.key_id),
                json_string(&g.scope),
                g.expires_day
            )
        });
        self.rows.push(format!(
            "{{\"event\":\"color-write\",\"node\":{id},\"name\":{},\
             \"color\":\"{}\",\"payload\":{},\"parents\":{:?},\"waiver\":{},\
             \"grant\":{},\"hash\":\"{}\"}}",
            json_string(name),
            color.name(),
            color.payload_json(),
            parents,
            waiver_json,
            grant_json,
            hash.to_hex()
        ));
        self.nodes.push(ColorNode {
            id,
            name: name.to_string(),
            color,
            parents,
            demotion,
            waiver,
            grant,
            hash,
        });
        id
    }

    /// Write a colored LEAF (a measurement, a certified input, an
    /// estimator output). Leaves state their color; derivations must
    /// EARN theirs.
    pub fn source(&mut self, name: &str, color: Color) -> u64 {
        self.push_node(name, color, Vec::new(), None, None, None)
    }

    /// Regime re-checks + composition fold shared by the derive paths.
    fn fold_parents(
        &self,
        parents: &[u64],
        op: IntervalOp,
        state: &BTreeMap<String, f64>,
    ) -> Result<(Color, Option<Demotion>), ColorWriteError> {
        if parents.is_empty() {
            return Err(ColorWriteError::NoParents);
        }
        for &p in parents {
            if p as usize >= self.nodes.len() {
                return Err(ColorWriteError::UnknownParent { id: p });
            }
        }
        let mut demotion = None;
        let mut effective: Vec<Color> = Vec::with_capacity(parents.len());
        for &p in parents {
            let (c, d) = check_regime(&self.nodes[p as usize].color, state);
            if demotion.is_none() {
                demotion = d;
            }
            effective.push(c);
        }
        let mut derived = effective[0].clone();
        for c in &effective[1..] {
            derived = compose(&derived, c, op);
        }
        Ok((derived, demotion))
    }

    fn laundering_error(
        &self,
        parents: &[u64],
        state: &BTreeMap<String, f64>,
        claimed: ColorRank,
        cap: ColorRank,
    ) -> ColorWriteError {
        let offending: Vec<u64> = parents
            .iter()
            .copied()
            .filter(|&p| {
                let (eff, _) = check_regime(&self.nodes[p as usize].color, state);
                eff.rank() <= cap
            })
            .collect();
        ColorWriteError::LaunderingRefused {
            claimed,
            derived: cap,
            offending_parents: offending,
        }
    }

    /// Write a DERIVED node: the composition algebra folds the parent
    /// colors (with regime re-checks against `state`, auto-demoting on
    /// exit), and the claimed color must not outrank the derivation.
    /// The `waiver` argument is a HUMAN ANNOTATION only (bead
    /// qmao.1.1): it is recorded and hashed but authorizes NOTHING —
    /// an upgrade claim is refused here regardless. The authorized
    /// path is [`ColorGraph::derive_waived`].
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
        let (derived, demotion) = self.fold_parents(parents, op, state)?;
        let written = match claimed {
            None => derived,
            Some(c) if c.rank() <= derived.rank() => c,
            Some(c) => {
                return Err(self.laundering_error(parents, state, c.rank(), derived.rank()));
            }
        };
        Ok(self.push_node(name, written, parents.to_vec(), demotion, waiver, None))
    }

    /// Write a DERIVED node whose upgrade past the composition cap is
    /// authorized by an AUTHENTICATED [`WaiverGrant`] (bead qmao.1.1):
    /// the grant must carry the color-upgrade scope, name THIS node,
    /// authorize exactly the claimed color, bind the exact parent
    /// provenance hashes (replay to another node fails), be unexpired
    /// as of `today_day`, and carry a signature the `verifier`
    /// capability accepts over the canonical length-prefixed payload.
    /// Any failure refuses the write (fail closed) — with the in-tree
    /// [`NoWaiverVerifier`] every promotion is refused (the no-crypto
    /// no-claim).
    ///
    /// # Errors
    /// [`ColorWriteError::WaiverRefused`] with the structured
    /// rejection, plus the ordinary derive errors.
    #[allow(clippy::too_many_arguments)] // the authorization surface is the point
    pub fn derive_waived(
        &mut self,
        name: &str,
        parents: &[u64],
        op: IntervalOp,
        claimed: Color,
        state: &BTreeMap<String, f64>,
        grant: WaiverGrant,
        verifier: &dyn WaiverVerifier,
        today_day: u32,
    ) -> Result<u64, ColorWriteError> {
        let (derived, demotion) = self.fold_parents(parents, op, state)?;
        if claimed.rank() <= derived.rank() {
            // No upgrade needed; record the annotation, drop nothing.
            return Ok(self.push_node(
                name,
                claimed,
                parents.to_vec(),
                demotion,
                Some(grant.annotation.clone()),
                Some(grant),
            ));
        }
        let refuse = |rejection| Err(ColorWriteError::WaiverRefused { rejection });
        if grant.scope != WAIVER_SCOPE_COLOR_UPGRADE {
            return refuse(WaiverRejection::ScopeMismatch);
        }
        if grant.node_name != name {
            return refuse(WaiverRejection::NodeMismatch);
        }
        if grant.claimed_color != claimed.name() {
            return refuse(WaiverRejection::ColorMismatch);
        }
        let lineage: Vec<ContentHash> = parents
            .iter()
            .map(|&p| self.nodes[p as usize].hash)
            .collect();
        if grant.parent_hashes != lineage {
            return refuse(WaiverRejection::LineageMismatch);
        }
        if today_day > grant.expires_day {
            return refuse(WaiverRejection::Expired);
        }
        if !verifier.verify(&grant.key_id, &grant.signing_payload(), &grant.signature) {
            return refuse(WaiverRejection::BadSignature);
        }
        Ok(self.push_node(
            name,
            claimed,
            parents.to_vec(),
            demotion,
            Some(grant.annotation.clone()),
            Some(grant),
        ))
    }

    /// The node by id.
    #[must_use]
    pub fn node(&self, id: u64) -> &ColorNode {
        &self.nodes[id as usize]
    }
}
