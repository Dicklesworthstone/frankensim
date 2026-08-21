//! E10.2b — the distributed exact H-campaign EXECUTOR (bead
//! wf-root-guzez.11.4, de-circularized per Round-5): consume a
//! HistoricalCampaignIntentManifestV1, execute deterministic
//! RESUMABLE member×flight shards through the REAL lifecycle engine,
//! merge by member index into a HistoricalCampaignReceiptV1.
//!
//! Laws:
//!   - the executor SCORES NOTHING (the receipt carries outcomes and
//!     digests; scoring is E10.2c's separate owner);
//!   - refusal accounting is exact: a member×flight whose physics
//!     refuses is a RECORDED row, never a dropped one;
//!   - surrogate members ALLOCATE (they size the exact-run budget in
//!     the intent) and never execute — a surrogate result appearing
//!     in the finals is the hostile twin and the merge refuses;
//!   - the merge digest is a pure function of the result SET: any
//!     shard partitioning of the same intent merges to the same
//!     digest, which is what makes shards resumable.

use crate::simloop::{Phase, PilotMode, ScenarioSpec, SimLoop, TerminalEvent};
use crate::{Refusal, refuse};
use fs_blake3::hash_domain;

/// Member cap per campaign.
pub const MAX_MEMBERS: usize = 64;
/// Flight cap per campaign.
pub const MAX_FLIGHTS: usize = 16;

/// One campaign member.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MemberSpec {
    /// Member index (identity within the campaign).
    pub member: u32,
    /// Study seed for this member.
    pub seed: u64,
    /// Surrogate members size the allocation; they NEVER execute.
    pub surrogate: bool,
}

/// One historical flight condition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlightSpec {
    /// Flight index (Dec-17: 0..=3).
    pub flight: u32,
    /// Headwind [m/s].
    pub headwind_mps: f64,
    /// Air density [kg/m³].
    pub rho_kg_m3: f64,
    /// Tick budget for the exact run.
    pub max_ticks: u64,
}

/// The campaign intent (pre-execution; its digest is the campaign id
/// ingredient).
#[derive(Clone, Debug, PartialEq)]
pub struct HistoricalCampaignIntentManifestV1 {
    /// Campaign label.
    pub label: &'static str,
    /// Members (exact + surrogate).
    pub members: Vec<MemberSpec>,
    /// Flights.
    pub flights: Vec<FlightSpec>,
}

impl HistoricalCampaignIntentManifestV1 {
    /// Admit + digest the intent.
    ///
    /// # Errors
    /// `campaign-intent-invalid` (caps at cap AND cap+1; duplicate
    /// member/flight indices; non-finite conditions; zero exact
    /// members).
    pub fn admit(&self) -> Result<String, Refusal> {
        let nm = self.members.len();
        let nf = self.flights.len();
        let mut member_ids: Vec<u32> = self.members.iter().map(|m| m.member).collect();
        member_ids.sort_unstable();
        member_ids.dedup();
        let mut flight_ids: Vec<u32> = self.flights.iter().map(|f| f.flight).collect();
        flight_ids.sort_unstable();
        flight_ids.dedup();
        let exact = self.members.iter().filter(|m| !m.surrogate).count();
        if nm == 0
            || nm > MAX_MEMBERS
            || nf == 0
            || nf > MAX_FLIGHTS
            || member_ids.len() != nm
            || flight_ids.len() != nf
            || exact == 0
            || self
                .flights
                .iter()
                .any(|f| !(f.headwind_mps.is_finite() && f.rho_kg_m3.is_finite()))
        {
            return Err(refuse(
                "campaign-intent-invalid",
                format!("{nm} members ({exact} exact), {nf} flights"),
                "1..=64 unique members (>=1 exact), 1..=16 unique flights, finite conditions",
            ));
        }
        let mut b = self.label.as_bytes().to_vec();
        for m in &self.members {
            b.extend_from_slice(&m.member.to_le_bytes());
            b.extend_from_slice(&m.seed.to_le_bytes());
            b.push(u8::from(m.surrogate));
        }
        for f in &self.flights {
            b.extend_from_slice(&f.flight.to_le_bytes());
            b.extend_from_slice(&f.headwind_mps.to_bits().to_le_bytes());
            b.extend_from_slice(&f.rho_kg_m3.to_bits().to_le_bytes());
            b.extend_from_slice(&f.max_ticks.to_le_bytes());
        }
        Ok(hash_domain("org.frankensim.wf.h-campaign-intent.v1", &b).to_hex())
    }

    /// The exact (executable) member×flight work list, deterministic
    /// order (member-major).
    #[must_use]
    pub fn work_units(&self) -> Vec<(MemberSpec, FlightSpec)> {
        let mut units = Vec::new();
        for m in self.members.iter().filter(|m| !m.surrogate) {
            for f in &self.flights {
                units.push((*m, *f));
            }
        }
        units
    }
}

/// One executed member×flight outcome.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RunOutcome {
    /// The run reached a terminal phase; the lifecycle digest is the
    /// exact identity of what happened.
    Completed {
        /// Terminal phase word.
        terminal: &'static str,
        /// Lifecycle digest (hex).
        digest: String,
    },
    /// The physics refused (typed) — RECORDED, never dropped.
    Refused {
        /// Refusal code.
        code: String,
    },
}

/// One receipt row.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CampaignRow {
    /// Member index.
    pub member: u32,
    /// Flight index.
    pub flight: u32,
    /// Outcome.
    pub outcome: RunOutcome,
}

/// A shard: a contiguous slice of the deterministic work list.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ShardResult {
    /// Intent digest the shard executed under.
    pub intent_digest: String,
    /// First work-unit index (inclusive).
    pub from: usize,
    /// Last work-unit index (exclusive).
    pub to: usize,
    /// Rows, work-list order.
    pub rows: Vec<CampaignRow>,
}

/// Execute one shard of the intent's exact work list.
///
/// # Errors
/// Intent refusals; `campaign-shard-invalid` (empty or out-of-range
/// slice — AT the end admits, one past refuses).
pub fn execute_shard(
    intent: &HistoricalCampaignIntentManifestV1,
    from: usize,
    to: usize,
) -> Result<ShardResult, Refusal> {
    let intent_digest = intent.admit()?;
    let units = intent.work_units();
    if from >= to || to > units.len() {
        return Err(refuse(
            "campaign-shard-invalid",
            format!("[{from}, {to}) of {} units", units.len()),
            "shards are non-empty in-range slices of the work list",
        ));
    }
    let mut rows = Vec::with_capacity(to - from);
    for (m, f) in &units[from..to] {
        let spec = ScenarioSpec {
            seed: m.seed,
            rho_kg_m3: f.rho_kg_m3,
            headwind_mps: f.headwind_mps,
            pilot_mode: PilotMode::Historical(0),
            assist: None,
            catapult: None,
            rail_length_m: 18.3,
            max_ticks: f.max_ticks,
        };
        let outcome = match SimLoop::init(spec) {
            Err(e) => RunOutcome::Refused {
                code: e.code.to_string(),
            },
            Ok(mut sim) => {
                let mut terminal = "max-ticks";
                loop {
                    match sim.step(None) {
                        Err(e) => {
                            terminal = "step-refused";
                            let _ = e;
                            break;
                        }
                        Ok(out) => match out.phase {
                            Phase::Ended(TerminalEvent::GroundContact) => {
                                terminal = "ground-contact";
                                break;
                            }
                            Phase::Ended(TerminalEvent::RailEndWithoutLift) => {
                                terminal = "rail-end-without-lift";
                                break;
                            }
                            Phase::Ended(TerminalEvent::MaxTicks) => {
                                terminal = "max-ticks";
                                break;
                            }
                            Phase::Ended(TerminalEvent::EnvelopeExceeded) => {
                                terminal = "envelope-exceeded";
                                break;
                            }
                            Phase::Ended(TerminalEvent::DamageModelUnavailable) => {
                                terminal = "damage-model-unavailable";
                                break;
                            }
                            _ => {}
                        },
                    }
                }
                RunOutcome::Completed {
                    terminal,
                    digest: sim.digest_hex(),
                }
            }
        };
        rows.push(CampaignRow {
            member: m.member,
            flight: f.flight,
            outcome,
        });
    }
    Ok(ShardResult {
        intent_digest,
        from,
        to,
        rows,
    })
}

/// The merged campaign receipt (NO scoring fields exist).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HistoricalCampaignReceiptV1 {
    /// Schema id.
    pub schema: &'static str,
    /// Intent digest.
    pub intent_digest: String,
    /// Rows, member-major (the member-index merge).
    pub rows: Vec<CampaignRow>,
    /// Completed count.
    pub completed: usize,
    /// Refused count (accounting law: completed + refused = rows).
    pub refused: usize,
    /// Merge digest (partition-independent).
    pub merge_digest: String,
}

/// Merge shards into the receipt. Partition-independent: any shard
/// split of the same intent merges to the same digest.
///
/// # Errors
/// `campaign-shards-incomplete` (gaps or overlaps in coverage — the
/// resume path is to execute the missing slice, never to pad);
/// `campaign-shards-mismatched` (a shard from a different intent);
/// `surrogate-in-finals` (a result row for a surrogate member — the
/// hostile twin).
pub fn merge_shards(
    intent: &HistoricalCampaignIntentManifestV1,
    shards: &[ShardResult],
) -> Result<HistoricalCampaignReceiptV1, Refusal> {
    let intent_digest = intent.admit()?;
    let units = intent.work_units();
    for s in shards {
        if s.intent_digest != intent_digest {
            return Err(refuse(
                "campaign-shards-mismatched",
                format!("shard [{}, {}) from another intent", s.from, s.to),
                "shards bind to exactly one intent digest",
            ));
        }
    }
    // Coverage check: every unit exactly once.
    let mut seen = vec![0u32; units.len()];
    for s in shards {
        for i in s.from..s.to {
            seen[i] += 1;
        }
    }
    if seen.iter().any(|c| *c != 1) {
        return Err(refuse(
            "campaign-shards-incomplete",
            format!(
                "coverage {:?} (gaps or overlaps)",
                seen.iter().filter(|c| **c != 1).count()
            ),
            "execute the missing slices; never pad or double-count",
        ));
    }
    // Surrogate hostile twin: no finals row may belong to a surrogate.
    let surrogates: Vec<u32> = intent
        .members
        .iter()
        .filter(|m| m.surrogate)
        .map(|m| m.member)
        .collect();
    let mut rows: Vec<CampaignRow> = shards.iter().flat_map(|s| s.rows.iter().cloned()).collect();
    if rows.iter().any(|r| surrogates.contains(&r.member)) {
        return Err(refuse(
            "surrogate-in-finals",
            "a surrogate member's result reached the finals".into(),
            "surrogates allocate; only exact members' runs are evidence",
        ));
    }
    // Member-index merge (member-major, then flight).
    rows.sort_by_key(|r| (r.member, r.flight));
    let completed = rows
        .iter()
        .filter(|r| matches!(r.outcome, RunOutcome::Completed { .. }))
        .count();
    let refused = rows.len() - completed;
    let mut b = intent_digest.as_bytes().to_vec();
    for r in &rows {
        b.extend_from_slice(&r.member.to_le_bytes());
        b.extend_from_slice(&r.flight.to_le_bytes());
        match &r.outcome {
            RunOutcome::Completed { terminal, digest } => {
                b.push(1);
                b.extend_from_slice(terminal.as_bytes());
                b.extend_from_slice(digest.as_bytes());
            }
            RunOutcome::Refused { code } => {
                b.push(0);
                b.extend_from_slice(code.as_bytes());
            }
        }
    }
    Ok(HistoricalCampaignReceiptV1 {
        schema: "org.frankensim.wf.h-campaign-receipt.v1",
        intent_digest,
        rows,
        completed,
        refused,
        merge_digest: hash_domain("org.frankensim.wf.h-campaign-merge.v1", &b).to_hex(),
    })
}
