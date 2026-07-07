//! TOMBSTONE LEDGER (addendum Proposal E): compounding swarm memory's
//! cheap half. A tombstone is an indexed, retrievable record of a
//! falsified hypothesis, so failed explorations are never silently
//! re-run — the swarm's negative results are just rows it declined to
//! write, so it writes them.
//!
//! Two indexes: the PRIMARY, domain-native π-space signature
//! (dimensionless groups via fs-regime, so "aluminum bracket at
//! Re 2×10⁵" and "steel bracket at Re 2.1×10⁵" collide as the SAME
//! death even though raw parameters differ) and a deterministic
//! feature-vector cosine embedding (hand-built over tokens +
//! π-coordinates — no external embedding model, per the Franken-only
//! rule). The orchestrator gate refuses to fund an exploration inside a
//! tombstone's neighborhood unless a VALIDATED distinguishing feature is
//! cited — and cited distinguishers are themselves logged, so they
//! accumulate.

use crate::{EventRow, Ledger, LedgerError};
use fs_qty::QtyAny;
use fs_regime::{Input, pi_groups};
use std::collections::BTreeMap;
use std::fmt::Write as _;

/// Neighborhood radius in log10 (decades) summed over π-coordinates:
/// two deaths within this distance collide.
pub const PI_RADIUS_DECADES: f64 = 0.1;

/// Cosine threshold for the embedding index.
pub const EMBED_MIN_COSINE: f64 = 0.95;

/// A distinguisher must move the named parameter by at least this many
/// decades to count as a genuine difference.
pub const DISTINGUISHER_MIN_DECADES: f64 = 0.05;

/// Embedding dimensionality (deterministic hashed features).
const EMBED_DIM: usize = 64;

/// The hypothesis/design-region descriptor: a name plus dimensioned
/// parameters (positive values — π-space is multiplicative).
#[derive(Debug, Clone, PartialEq)]
pub struct Descriptor {
    /// Human/agent-readable hypothesis name.
    pub name: String,
    /// Named, dimensioned parameters (BTreeMap: deterministic order).
    pub params: BTreeMap<String, QtyAny>,
}

/// The π-space signature: exponent vectors over the (name-sorted)
/// parameters plus the group's log10 value.
#[derive(Debug, Clone, PartialEq)]
pub struct PiSignature {
    /// Parameter names in the order the exponents index.
    pub basis: Vec<String>,
    /// (integer exponents, log10 of the group value) per group.
    pub groups: Vec<(Vec<i64>, f64)>,
}

fn log10(x: f64) -> f64 {
    fs_math_ln(x) / core::f64::consts::LN_10
}

// fs-ledger already links fs-evidence (UTIL); avoid a new fs-math dep by
// the std path — determinism here is index-internal, not a P2 artifact.
fn fs_math_ln(x: f64) -> f64 {
    x.ln()
}

impl Descriptor {
    /// Compute the π-space signature via fs-regime's exact Buckingham
    /// machinery.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] when parameters are empty/non-positive
    /// (π-space is multiplicative — a signed quantity has no signature).
    pub fn pi_signature(&self) -> Result<PiSignature, LedgerError> {
        if self.params.is_empty() {
            return Err(LedgerError::Invalid {
                field: "descriptor.params".to_string(),
                problem: "descriptor has no parameters".to_string(),
            });
        }
        let inputs: Vec<Input> = self
            .params
            .iter()
            .map(|(name, qty)| Input {
                name: name.clone(),
                qty: *qty,
            })
            .collect();
        let basis = pi_groups(&inputs).map_err(|e| LedgerError::Invalid {
            field: "descriptor.params".to_string(),
            problem: format!("pi signature failed: {e}"),
        })?;
        let mut groups: Vec<(Vec<i64>, f64)> = basis
            .groups
            .iter()
            .map(|g| (g.exponents.clone(), log10(g.value)))
            .collect();
        groups.sort_by(|a, b| a.0.cmp(&b.0));
        Ok(PiSignature {
            basis: self.params.keys().cloned().collect(),
            groups,
        })
    }

    /// The deterministic feature-vector embedding: hashed name tokens +
    /// hashed parameter names + π log-values folded into fixed buckets.
    #[must_use]
    pub fn embedding(&self) -> [f64; EMBED_DIM] {
        let mut v = [0.0f64; EMBED_DIM];
        let bucket = |text: &str| -> usize {
            (fs_obs_fnv(text.as_bytes()) as usize) % EMBED_DIM
        };
        for token in self.name.split(|c: char| !c.is_alphanumeric()) {
            if !token.is_empty() {
                v[bucket(&token.to_ascii_lowercase())] += 1.0;
            }
        }
        for (name, qty) in &self.params {
            v[bucket(name)] += 1.0;
            // Fold the magnitude class (decade) in, so scale matters but
            // small perturbations do not.
            if qty.value > 0.0 {
                let decade = log10(qty.value).round();
                v[bucket(&format!("{name}:{decade}"))] += 1.0;
            }
        }
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-300);
        for x in &mut v {
            *x /= norm;
        }
        v
    }
}

fn fs_obs_fnv(bytes: &[u8]) -> u64 {
    // FNV-1a, inline to keep this module self-contained.
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    h
}

/// π-space distance: sum over matched exponent vectors of |Δlog10|;
/// `None` when the group structures differ (different physics — no
/// collision claim either way).
#[must_use]
pub fn pi_distance(a: &PiSignature, b: &PiSignature) -> Option<f64> {
    if a.groups.len() != b.groups.len() || a.groups.is_empty() {
        return None;
    }
    let mut d = 0.0f64;
    for ((ea, va), (eb, vb)) in a.groups.iter().zip(&b.groups) {
        if ea != eb {
            return None;
        }
        d += (va - vb).abs();
    }
    Some(d)
}

fn cosine(a: &[f64; EMBED_DIM], b: &[f64; EMBED_DIM]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// One tombstone: a falsified hypothesis with its evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct TombstoneRecord {
    /// The dead hypothesis.
    pub descriptor: Descriptor,
    /// Evidence against (falsifier detail, branch metrics…).
    pub evidence: String,
    /// Certificate colors involved (from the three-color schema).
    pub colors: Vec<String>,
    /// Compute spent before death (seconds).
    pub compute_spent_s: f64,
    /// Logical date (caller-supplied; determinism).
    pub date: String,
    /// Authoring context (agent/session identity).
    pub context: String,
    /// Distinguishers cited against this tombstone (they accumulate).
    pub distinguishers: Vec<String>,
}

impl TombstoneRecord {
    /// Canonical JSON (the ledger event payload).
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = format!(
            "{{\"kind\":\"tombstone\",\"name\":{:?},\"evidence\":{:?},\"colors\":[",
            self.descriptor.name, self.evidence
        );
        for (i, c) in self.colors.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(s, "{c:?}");
        }
        let _ = write!(
            s,
            "],\"compute_s\":{},\"date\":{:?},\"context\":{:?},\"params\":{{",
            self.compute_spent_s, self.date, self.context
        );
        for (i, (k, v)) in self.descriptor.params.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(s, "{k:?}:{}", v.value);
        }
        s.push_str("},\"distinguishers\":[");
        for (i, d) in self.distinguishers.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(s, "{d:?}");
        }
        s.push_str("]}");
        s
    }
}

/// The orchestrator's pre-exploration verdict.
#[derive(Debug, Clone, PartialEq)]
pub enum ExplorationVerdict {
    /// No tombstone nearby: fund it.
    Clear,
    /// Inside a tombstone's neighborhood: SKIP, or cite a validated
    /// distinguishing feature.
    Blocked {
        /// Indexes of the colliding tombstones (π-space first).
        neighbors: Vec<usize>,
        /// Which index fired ("pi-space" / "embedding").
        via: &'static str,
    },
}

/// A refused distinguisher (free text is not accepted).
#[derive(Debug, Clone, PartialEq)]
pub struct DistinguisherRefusal {
    /// Why.
    pub what: String,
}

/// The in-memory tombstone index (rebuildable from ledger rows).
#[derive(Debug, Default)]
pub struct TombstoneIndex {
    records: Vec<TombstoneRecord>,
    signatures: Vec<Option<PiSignature>>,
    embeddings: Vec<[f64; EMBED_DIM]>,
    /// Metrics: (funded clear, blocked hits, funded via distinguisher).
    stats: (u64, u64, u64),
}

impl TombstoneIndex {
    /// An empty index.
    #[must_use]
    pub fn new() -> Self {
        TombstoneIndex::default()
    }

    /// Number of tombstones.
    #[must_use]
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// True when no tombstones exist.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// A tombstone by index.
    #[must_use]
    pub fn get(&self, i: usize) -> Option<&TombstoneRecord> {
        self.records.get(i)
    }

    fn append(&mut self, record: TombstoneRecord) -> usize {
        self.signatures.push(record.descriptor.pi_signature().ok());
        self.embeddings.push(record.descriptor.embedding());
        self.records.push(record);
        self.records.len() - 1
    }

    /// AUTOMATIC APPEND on a falsification kill (Proposal 6 wiring): the
    /// falsifier's tombstone JSON is carried as evidence verbatim.
    pub fn record_falsification_kill(
        &mut self,
        descriptor: Descriptor,
        falsifier_tombstone_json: &str,
        colors: Vec<String>,
        compute_spent_s: f64,
        date: &str,
        context: &str,
    ) -> usize {
        self.append(TombstoneRecord {
            descriptor,
            evidence: falsifier_tombstone_json.to_string(),
            colors,
            compute_spent_s,
            date: date.to_string(),
            context: context.to_string(),
            distinguishers: Vec::new(),
        })
    }

    /// AUTOMATIC APPEND on an abandoned optimization branch — but only
    /// above the cost threshold (cheap failures are noise, not memory).
    /// Returns the index, or `None` below threshold.
    pub fn record_abandoned_branch(
        &mut self,
        descriptor: Descriptor,
        best_objective: f64,
        compute_spent_s: f64,
        cost_threshold_s: f64,
        date: &str,
        context: &str,
    ) -> Option<usize> {
        if compute_spent_s < cost_threshold_s {
            return None;
        }
        Some(self.append(TombstoneRecord {
            descriptor,
            evidence: format!(
                "optimization branch abandoned at objective {best_objective} after \
                 {compute_spent_s}s"
            ),
            colors: vec!["estimated".to_string()],
            compute_spent_s,
            date: date.to_string(),
            context: context.to_string(),
            distinguishers: Vec::new(),
        }))
    }

    /// π-space neighbors within the radius (the PRIMARY index).
    #[must_use]
    pub fn pi_neighbors(&self, descriptor: &Descriptor) -> Vec<usize> {
        let Ok(sig) = descriptor.pi_signature() else {
            return Vec::new();
        };
        self.signatures
            .iter()
            .enumerate()
            .filter_map(|(i, s)| {
                s.as_ref()
                    .and_then(|s| pi_distance(&sig, s))
                    .filter(|&d| d <= PI_RADIUS_DECADES)
                    .map(|_| i)
            })
            .collect()
    }

    /// Embedding neighbors above the cosine threshold (secondary index).
    #[must_use]
    pub fn embed_neighbors(&self, descriptor: &Descriptor) -> Vec<usize> {
        let e = descriptor.embedding();
        self.embeddings
            .iter()
            .enumerate()
            .filter(|(_, x)| cosine(&e, x) >= EMBED_MIN_COSINE)
            .map(|(i, _)| i)
            .collect()
    }

    /// THE ORCHESTRATOR GATE: query before funding. π-space is primary;
    /// the embedding index catches descriptor-similar deaths whose
    /// π-signatures were unavailable.
    pub fn pre_exploration_check(&mut self, descriptor: &Descriptor) -> ExplorationVerdict {
        let pi = self.pi_neighbors(descriptor);
        if !pi.is_empty() {
            self.stats.1 += 1;
            return ExplorationVerdict::Blocked {
                neighbors: pi,
                via: "pi-space",
            };
        }
        let em = self.embed_neighbors(descriptor);
        if !em.is_empty() {
            self.stats.1 += 1;
            return ExplorationVerdict::Blocked {
                neighbors: em,
                via: "embedding",
            };
        }
        self.stats.0 += 1;
        ExplorationVerdict::Clear
    }

    /// Fund a BLOCKED exploration by citing a distinguishing feature.
    /// The distinguisher is VALIDATED, not accepted as free text: it must
    /// name a parameter present in the new descriptor whose value differs
    /// from the tombstone's by at least [`DISTINGUISHER_MIN_DECADES`]
    /// (or be absent from the tombstone entirely). Accepted
    /// distinguishers are LOGGED on the tombstone, so they accumulate.
    ///
    /// # Errors
    /// [`DistinguisherRefusal`] naming what failed.
    pub fn fund_with_distinguisher(
        &mut self,
        descriptor: &Descriptor,
        neighbor: usize,
        distinguisher_param: &str,
    ) -> Result<(), DistinguisherRefusal> {
        let Some(new_value) = descriptor.params.get(distinguisher_param) else {
            return Err(DistinguisherRefusal {
                what: format!(
                    "distinguisher {distinguisher_param:?} is not a parameter of the \
                     proposed exploration — cite a real feature"
                ),
            });
        };
        let Some(tomb) = self.records.get_mut(neighbor) else {
            return Err(DistinguisherRefusal {
                what: format!("tombstone index {neighbor} does not exist"),
            });
        };
        if let Some(old_value) = tomb.descriptor.params.get(distinguisher_param) {
            if new_value.dims != old_value.dims {
                // Different dimensions for the same name: genuinely distinct.
            } else if new_value.value > 0.0 && old_value.value > 0.0 {
                let delta = (log10(new_value.value) - log10(old_value.value)).abs();
                if delta < DISTINGUISHER_MIN_DECADES {
                    return Err(DistinguisherRefusal {
                        what: format!(
                            "distinguisher {distinguisher_param:?} differs by only \
                             {delta:.3} decades (< {DISTINGUISHER_MIN_DECADES}) from the \
                             tombstone — that is the same death"
                        ),
                    });
                }
            }
        }
        tomb.distinguishers
            .push(format!("{distinguisher_param}={}", new_value.value));
        self.stats.2 += 1;
        Ok(())
    }

    /// The RE-EXPLORATION-RATE metric (the proposal's kill criterion
    /// input): `(clear, blocked, funded_via_distinguisher, rate)` where
    /// rate = blocked / (clear + blocked).
    #[must_use]
    pub fn re_exploration_rate(&self) -> (u64, u64, u64, f64) {
        let (clear, blocked, funded) = self.stats;
        let total = clear + blocked;
        #[allow(clippy::cast_precision_loss)]
        let rate = if total == 0 {
            0.0
        } else {
            blocked as f64 / total as f64
        };
        (clear, blocked, funded, rate)
    }

    /// Persist every tombstone as a ledger event (kind "tombstone").
    ///
    /// # Errors
    /// Propagates ledger write failures.
    pub fn flush_to_ledger(&self, ledger: &Ledger) -> Result<(), LedgerError> {
        for (i, r) in self.records.iter().enumerate() {
            let payload = r.to_json();
            ledger.append_event(&EventRow {
                session: None,
                t: i64::try_from(i).unwrap_or(i64::MAX),
                kind: "tombstone",
                payload: Some(&payload),
            })?;
        }
        Ok(())
    }
}
