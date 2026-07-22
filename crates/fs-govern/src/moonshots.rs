//! Fixed-size moonshot portfolio policy (bead
//! `frankensim-extreal-program-f85xj.16.3`).
//!
//! This module specializes [`crate::lanes::PortfolioLedger`] for `[M]`
//! research. A declaration must name an owner, falsifier, effort/calendar
//! budget, quarterly review, boring fallback, and critical-path-disjointness
//! boundary. New work enters only through a replacement that terminalizes one
//! named active lane as completed, falsified, or shelved with retained state.
//!
//! The state machine is deliberately ignorant of Beads and product dependency
//! graphs. The repository adapter in `xtask` binds declarations to those live
//! authorities and rejects direct or transitive critical-path blocking.

use crate::lanes::{
    FinalizationReceipt, IdempotencyKey, LaneCharter, LaneError, MechanismId, PortfolioLedger,
    PortfolioPolicy, ResourceEnvelope, TerminalKind,
};
use std::collections::{BTreeMap, BTreeSet};

/// Version of the fixed-size moonshot portfolio policy.
pub const MOONSHOT_POLICY_VERSION: u32 = 1;

/// Immutable maximum active count captured by the v1 status-quo inventory.
///
/// A live portfolio may shrink below this ceiling, but cannot grow back into a
/// slot released without a simultaneous replacement.
pub const MOONSHOT_V1_INITIAL_CAP: u32 = 6;

/// Maximum UTF-8 bytes in one moonshot policy text field.
pub const MAX_MOONSHOT_FIELD_BYTES: usize = 4096;

/// Domain for disposition evidence passed to the underlying lane ledger.
pub const MOONSHOT_DISPOSITION_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.moonshot-disposition.v1";

/// Structured refusal from the moonshot portfolio policy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MoonshotError {
    /// A required text field was empty after trimming.
    EmptyField {
        /// Stable field name.
        field: &'static str,
    },
    /// A text field exceeded [`MAX_MOONSHOT_FIELD_BYTES`].
    FieldTooLarge {
        /// Stable field name.
        field: &'static str,
        /// Observed byte length.
        len: usize,
        /// Maximum admitted byte length.
        cap: usize,
    },
    /// The effort budget was zero.
    ZeroEffortBudget,
    /// A calendar field was not a real `YYYY-MM-DD` date.
    InvalidDate {
        /// Stable field name.
        field: &'static str,
    },
    /// The next falsifier review was later than the lane deadline.
    ReviewAfterDeadline,
    /// The declared active cap exceeded the immutable v1 ceiling.
    CapExceeded {
        /// Requested cap.
        cap: u32,
        /// Immutable ceiling.
        ceiling: u32,
    },
    /// The cap and number of active declarations differed.
    ActiveCountMismatch {
        /// Declared cap.
        cap: u32,
        /// Number of declarations.
        declarations: usize,
    },
    /// Adding the active effort budgets overflowed `u64`.
    WorkBudgetOverflow,
    /// Two active declarations named the same Bead.
    DuplicateActiveBead {
        /// Duplicate Bead id.
        bead_id: String,
    },
    /// A replacement candidate was incorrectly marked as baseline work.
    CandidateIsLegacyBaseline,
    /// The candidate Bead was already active.
    CandidateAlreadyActive {
        /// Active candidate Bead id.
        bead_id: String,
    },
    /// A candidate attempted to displace itself.
    SelfDisplacement,
    /// The named displaced Bead was not active.
    DisplacedBeadNotActive {
        /// Missing active Bead id.
        bead_id: String,
    },
    /// A terminalized Bead attempted to return as a replacement.
    TerminalBeadRevival {
        /// Terminal Bead id.
        bead_id: String,
    },
    /// The underlying proof-lane admission ledger refused the operation.
    Lane(LaneError),
    /// A supposedly one-for-one replacement changed active WIP.
    ReplacementChangedActiveCount {
        /// Policy cap.
        expected: u32,
        /// Ledger count after the attempted replacement.
        actual: u32,
    },
}

impl core::fmt::Display for MoonshotError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyField { field } => write!(f, "{field} must be non-empty"),
            Self::FieldTooLarge { field, len, cap } => {
                write!(f, "{field} is {len} bytes; cap is {cap}")
            }
            Self::ZeroEffortBudget => f.write_str("budget.effort_minutes must be positive"),
            Self::InvalidDate { field } => {
                write!(f, "{field} must be a real YYYY-MM-DD date")
            }
            Self::ReviewAfterDeadline => f.write_str(
                "quarterly_review.next_review must not be later than budget.calendar_deadline",
            ),
            Self::CapExceeded { cap, ceiling } => {
                write!(
                    f,
                    "moonshot cap {cap} exceeds immutable v1 ceiling {ceiling}"
                )
            }
            Self::ActiveCountMismatch { cap, declarations } => write!(
                f,
                "moonshot cap must equal declared active WIP: cap {cap}, declarations {declarations}"
            ),
            Self::WorkBudgetOverflow => f.write_str("moonshot work budget sum overflowed"),
            Self::DuplicateActiveBead { bead_id } => {
                write!(f, "duplicate active moonshot declaration {bead_id:?}")
            }
            Self::CandidateIsLegacyBaseline => {
                f.write_str("replacement candidate cannot be a legacy-baseline declaration")
            }
            Self::CandidateAlreadyActive { bead_id } => {
                write!(f, "replacement candidate {bead_id:?} is already active")
            }
            Self::SelfDisplacement => f.write_str("a moonshot cannot displace itself"),
            Self::DisplacedBeadNotActive { bead_id } => {
                write!(f, "replacement must name active moonshot {bead_id:?}")
            }
            Self::TerminalBeadRevival { bead_id } => {
                write!(f, "terminal moonshot {bead_id:?} cannot be revived")
            }
            Self::Lane(error) => write!(f, "proof-lane ledger refused moonshot policy: {error}"),
            Self::ReplacementChangedActiveCount { expected, actual } => write!(
                f,
                "replacement changed active WIP: expected {expected}, found {actual}"
            ),
        }
    }
}

impl std::error::Error for MoonshotError {}

impl From<LaneError> for MoonshotError {
    fn from(value: LaneError) -> Self {
        Self::Lane(value)
    }
}

fn validated_text(field: &'static str, raw: &str) -> Result<String, MoonshotError> {
    if raw.len() > MAX_MOONSHOT_FIELD_BYTES {
        return Err(MoonshotError::FieldTooLarge {
            field,
            len: raw.len(),
            cap: MAX_MOONSHOT_FIELD_BYTES,
        });
    }
    let value = raw.trim();
    if value.is_empty() {
        return Err(MoonshotError::EmptyField { field });
    }
    Ok(value.to_string())
}

fn date_key(value: &str) -> Option<(u16, u8, u8)> {
    if value.len() != 10
        || value.as_bytes().get(4) != Some(&b'-')
        || value.as_bytes().get(7) != Some(&b'-')
        || !value
            .bytes()
            .enumerate()
            .all(|(index, byte)| index == 4 || index == 7 || byte.is_ascii_digit())
    {
        return None;
    }
    let year = value.get(..4)?.parse().ok()?;
    let month: u8 = value.get(5..7)?.parse().ok()?;
    let day: u8 = value.get(8..)?.parse().ok()?;
    let max_day = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 if year % 400 == 0 || (year % 4 == 0 && year % 100 != 0) => 29,
        2 => 28,
        _ => return None,
    };
    (day != 0 && day <= max_day).then_some((year, month, day))
}

/// Named observation and decision rule capable of killing a moonshot.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NamedFalsifier {
    id: String,
    observation: String,
    decision_rule: String,
}

impl NamedFalsifier {
    /// Construct a bounded, non-empty falsifier declaration.
    pub fn new(id: &str, observation: &str, decision_rule: &str) -> Result<Self, MoonshotError> {
        Ok(Self {
            id: validated_text("falsifier.id", id)?,
            observation: validated_text("falsifier.observation", observation)?,
            decision_rule: validated_text("falsifier.decision_rule", decision_rule)?,
        })
    }

    /// Stable falsifier id.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Observation that triggers falsification.
    #[must_use]
    pub fn observation(&self) -> &str {
        &self.observation
    }

    /// Required policy action after the observation.
    #[must_use]
    pub fn decision_rule(&self) -> &str {
        &self.decision_rule
    }
}

/// Effort and calendar cap for one moonshot.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoonshotBudget {
    effort_minutes: u64,
    calendar_deadline: String,
}

impl MoonshotBudget {
    /// Construct a positive effort budget with a real calendar deadline.
    pub fn new(effort_minutes: u64, calendar_deadline: &str) -> Result<Self, MoonshotError> {
        if effort_minutes == 0 {
            return Err(MoonshotError::ZeroEffortBudget);
        }
        if date_key(calendar_deadline).is_none() {
            return Err(MoonshotError::InvalidDate {
                field: "budget.calendar_deadline",
            });
        }
        Ok(Self {
            effort_minutes,
            calendar_deadline: calendar_deadline.to_string(),
        })
    }

    /// Maximum charged effort in minutes.
    #[must_use]
    pub const fn effort_minutes(&self) -> u64 {
        self.effort_minutes
    }

    /// Calendar deadline in canonical `YYYY-MM-DD` form.
    #[must_use]
    pub fn calendar_deadline(&self) -> &str {
        &self.calendar_deadline
    }
}

/// Scheduled falsifier review and the evidence it must inspect.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuarterlyReview {
    next_review: String,
    status: String,
    required_evidence: String,
}

impl QuarterlyReview {
    /// Construct a bounded review declaration.
    pub fn new(
        next_review: &str,
        status: &str,
        required_evidence: &str,
    ) -> Result<Self, MoonshotError> {
        if date_key(next_review).is_none() {
            return Err(MoonshotError::InvalidDate {
                field: "quarterly_review.next_review",
            });
        }
        Ok(Self {
            next_review: next_review.to_string(),
            status: validated_text("quarterly_review.status", status)?,
            required_evidence: validated_text(
                "quarterly_review.required_evidence",
                required_evidence,
            )?,
        })
    }

    /// Next review date in canonical `YYYY-MM-DD` form.
    #[must_use]
    pub fn next_review(&self) -> &str {
        &self.next_review
    }

    /// Review lifecycle status.
    #[must_use]
    pub fn status(&self) -> &str {
        &self.status
    }

    /// Evidence required at the review.
    #[must_use]
    pub fn required_evidence(&self) -> &str {
        &self.required_evidence
    }
}

/// Complete declaration for one active or proposed moonshot lane.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoonshotDeclaration {
    bead_id: String,
    lane: String,
    owner: String,
    legacy_baseline: bool,
    falsifier: NamedFalsifier,
    budget: MoonshotBudget,
    critical_path_disjointness: String,
    review: QuarterlyReview,
    charter: LaneCharter,
    mechanism: MechanismId,
    reservation: ResourceEnvelope,
}

impl MoonshotDeclaration {
    /// Construct a declaration and its underlying proof-lane authority.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bead_id: &str,
        lane: &str,
        owner: &str,
        legacy_baseline: bool,
        falsifier: NamedFalsifier,
        budget: MoonshotBudget,
        critical_path_disjointness: &str,
        review: QuarterlyReview,
    ) -> Result<Self, MoonshotError> {
        let bead_id = validated_text("bead_id", bead_id)?;
        let lane = validated_text("lane", lane)?;
        let owner = validated_text("owner", owner)?;
        let critical_path_disjointness =
            validated_text("critical_path_disjointness", critical_path_disjointness)?;
        if date_key(review.next_review()) > date_key(budget.calendar_deadline()) {
            return Err(MoonshotError::ReviewAfterDeadline);
        }
        let statement = format!("moonshot lane {lane}: {bead_id}");
        let baseline = format!("fallback required by rule: {}", falsifier.decision_rule());
        let falsifier_family = format!("{}: {}", falsifier.id(), falsifier.observation());
        let assumptions = [
            format!("owner is {owner}"),
            format!("calendar deadline is {}", budget.calendar_deadline()),
            format!(
                "review {} is {} and requires {}",
                review.next_review(),
                review.status(),
                review.required_evidence()
            ),
        ];
        let assumption_refs: Vec<_> = assumptions.iter().map(String::as_str).collect();
        let charter = LaneCharter::new(
            &statement,
            &critical_path_disjointness,
            &assumption_refs,
            "research-only descriptive authority",
            &baseline,
            &falsifier_family,
            &lane,
        )?;
        let mechanism = charter.mechanism_id(&bead_id, MOONSHOT_POLICY_VERSION)?;
        let reservation = ResourceEnvelope {
            work_units: budget.effort_minutes(),
            memory_bytes: 0,
            reviewer_slots: 1,
            falsification_capacity: 1,
        };
        Ok(Self {
            bead_id,
            lane,
            owner,
            legacy_baseline,
            falsifier,
            budget,
            critical_path_disjointness,
            review,
            charter,
            mechanism,
            reservation,
        })
    }

    /// Authoritative Bead id for this declaration.
    #[must_use]
    pub fn bead_id(&self) -> &str {
        &self.bead_id
    }

    /// Independently falsifiable lane name.
    #[must_use]
    pub fn lane(&self) -> &str {
        &self.lane
    }

    /// Named owner.
    #[must_use]
    pub fn owner(&self) -> &str {
        &self.owner
    }

    /// Whether this declaration belongs to the immutable v1 baseline.
    #[must_use]
    pub const fn is_legacy_baseline(&self) -> bool {
        self.legacy_baseline
    }

    /// Named falsifier.
    #[must_use]
    pub const fn falsifier(&self) -> &NamedFalsifier {
        &self.falsifier
    }

    /// Effort and calendar budget.
    #[must_use]
    pub const fn budget(&self) -> &MoonshotBudget {
        &self.budget
    }

    /// Declared product-path disjointness boundary.
    #[must_use]
    pub fn critical_path_disjointness(&self) -> &str {
        &self.critical_path_disjointness
    }

    /// Scheduled falsifier review.
    #[must_use]
    pub const fn review(&self) -> &QuarterlyReview {
        &self.review
    }

    /// Canonical proof-lane charter used by the admission ledger.
    #[must_use]
    pub const fn charter(&self) -> &LaneCharter {
        &self.charter
    }

    /// Mechanism identity minted by the charter.
    #[must_use]
    pub const fn mechanism(&self) -> MechanismId {
        self.mechanism
    }

    /// Four-axis reservation projected from this declaration.
    #[must_use]
    pub const fn reservation(&self) -> ResourceEnvelope {
        self.reservation
    }
}

/// Terminal disposition required before a moonshot slot can move.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MoonshotDisposition {
    /// The declared work completed within its no-claim boundary.
    Completed,
    /// The named falsifier fired.
    Falsified,
    /// Work stopped with enough retained state for later review.
    ShelvedWithState,
}

impl MoonshotDisposition {
    fn terminal_kind(self) -> TerminalKind {
        match self {
            Self::Completed | Self::ShelvedWithState => TerminalKind::Withdrawn,
            Self::Falsified => TerminalKind::Refuted,
        }
    }

    fn code(self) -> &'static str {
        match self {
            Self::Completed => "completed",
            Self::Falsified => "falsified",
            Self::ShelvedWithState => "shelved-with-state",
        }
    }
}

/// Retained liquidation record for a named active lane.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DisplacementRecord {
    bead_id: String,
    disposition: MoonshotDisposition,
    state_artifact: String,
    reason: String,
}

impl DisplacementRecord {
    /// Construct a complete terminal disposition.
    pub fn new(
        bead_id: &str,
        disposition: MoonshotDisposition,
        state_artifact: &str,
        reason: &str,
    ) -> Result<Self, MoonshotError> {
        Ok(Self {
            bead_id: validated_text("displaces.bead_id", bead_id)?,
            disposition,
            state_artifact: validated_text("displaces.state_artifact", state_artifact)?,
            reason: validated_text("displaces.reason", reason)?,
        })
    }

    /// Bead leaving the active portfolio.
    #[must_use]
    pub fn bead_id(&self) -> &str {
        &self.bead_id
    }

    /// Exact terminal disposition.
    #[must_use]
    pub const fn disposition(&self) -> MoonshotDisposition {
        self.disposition
    }

    /// Retained evidence or resumable-state locator.
    #[must_use]
    pub fn state_artifact(&self) -> &str {
        &self.state_artifact
    }

    /// Human-auditable liquidation reason.
    #[must_use]
    pub fn reason(&self) -> &str {
        &self.reason
    }
}

/// One-for-one proposal to replace an active moonshot.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReplacementAdmission {
    candidate: MoonshotDeclaration,
    displacement: DisplacementRecord,
}

impl ReplacementAdmission {
    /// Construct a replacement request and reject obvious self/legacy cases.
    pub fn new(
        candidate: MoonshotDeclaration,
        displacement: DisplacementRecord,
    ) -> Result<Self, MoonshotError> {
        if candidate.is_legacy_baseline() {
            return Err(MoonshotError::CandidateIsLegacyBaseline);
        }
        if candidate.bead_id() == displacement.bead_id() {
            return Err(MoonshotError::SelfDisplacement);
        }
        Ok(Self {
            candidate,
            displacement,
        })
    }

    /// Proposed candidate.
    #[must_use]
    pub const fn candidate(&self) -> &MoonshotDeclaration {
        &self.candidate
    }

    /// Named lane disposition that releases the slot.
    #[must_use]
    pub const fn displacement(&self) -> &DisplacementRecord {
        &self.displacement
    }
}

#[derive(Debug)]
struct ActiveMoonshot {
    mechanism: MechanismId,
}

fn finalization_receipt(
    mechanism: MechanismId,
    displacement: &DisplacementRecord,
) -> Result<FinalizationReceipt, MoonshotError> {
    let evidence = fs_blake3::hash_domain(
        MOONSHOT_DISPOSITION_IDENTITY_DOMAIN,
        format!(
            "{}\0{}\0{}\0{}",
            displacement.disposition().code(),
            displacement.state_artifact(),
            displacement.reason(),
            displacement.bead_id()
        )
        .as_bytes(),
    );
    FinalizationReceipt::new(
        mechanism,
        displacement.disposition().terminal_kind(),
        None,
        evidence,
    )
    .map_err(Into::into)
}

/// Non-clone fixed-size moonshot authority over the generic proof-lane ledger.
///
/// Replacement consumes `self`. A failed second step drops the partially
/// evaluated candidate ledger, so callers can never observe a released slot
/// without the admitted successor.
#[derive(Debug)]
pub struct MoonshotPortfolio {
    cap: u32,
    ledger: PortfolioLedger,
    active: BTreeMap<String, ActiveMoonshot>,
    terminal_beads: BTreeSet<String>,
    displacements: Vec<DisplacementRecord>,
}

impl MoonshotPortfolio {
    /// Admit the complete current portfolio under a fixed cap.
    pub fn new(cap: u32, declarations: Vec<MoonshotDeclaration>) -> Result<Self, MoonshotError> {
        if cap > MOONSHOT_V1_INITIAL_CAP {
            return Err(MoonshotError::CapExceeded {
                cap,
                ceiling: MOONSHOT_V1_INITIAL_CAP,
            });
        }
        if usize::try_from(cap).ok() != Some(declarations.len()) {
            return Err(MoonshotError::ActiveCountMismatch {
                cap,
                declarations: declarations.len(),
            });
        }
        let total_work = declarations.iter().try_fold(0_u64, |total, declaration| {
            total.checked_add(declaration.reservation().work_units)
        });
        let Some(total_work) = total_work else {
            return Err(MoonshotError::WorkBudgetOverflow);
        };
        let policy = PortfolioPolicy {
            global: ResourceEnvelope {
                work_units: total_work,
                memory_bytes: 0,
                reviewer_slots: u64::from(cap),
                falsification_capacity: u64::from(cap),
            },
            max_active_mechanisms: cap,
        };
        let mut ledger = PortfolioLedger::new(policy);
        let mut active = BTreeMap::new();
        for declaration in declarations {
            if active.contains_key(declaration.bead_id()) {
                return Err(MoonshotError::DuplicateActiveBead {
                    bead_id: declaration.bead_id().to_string(),
                });
            }
            ledger.admit(
                declaration.charter(),
                declaration.mechanism(),
                declaration.reservation(),
                IdempotencyKey::derive(&format!(
                    "moonshot-current-admit:{}",
                    declaration.bead_id()
                )),
            )?;
            active.insert(
                declaration.bead_id().to_string(),
                ActiveMoonshot {
                    mechanism: declaration.mechanism(),
                },
            );
        }
        Ok(Self {
            cap,
            ledger,
            active,
            terminal_beads: BTreeSet::new(),
            displacements: Vec::new(),
        })
    }

    /// Terminalize one active moonshot without a successor and permanently
    /// shrink the available WIP cap.
    pub fn liquidate(mut self, displacement: &DisplacementRecord) -> Result<Self, MoonshotError> {
        let displaced = self.active.get(displacement.bead_id()).ok_or_else(|| {
            MoonshotError::DisplacedBeadNotActive {
                bead_id: displacement.bead_id().to_string(),
            }
        })?;
        let receipt = finalization_receipt(displaced.mechanism, displacement)?;
        self.ledger.finalize(
            &receipt,
            IdempotencyKey::derive(&format!("moonshot-liquidate:{}", displacement.bead_id())),
        )?;
        self.active.remove(displacement.bead_id());
        self.terminal_beads
            .insert(displacement.bead_id().to_string());
        self.displacements.push(displacement.clone());
        self.cap = self
            .cap
            .checked_sub(1)
            .ok_or(MoonshotError::ReplacementChangedActiveCount {
                expected: 0,
                actual: self.ledger.active_count(),
            })?;
        if self.ledger.active_count() != self.cap {
            return Err(MoonshotError::ReplacementChangedActiveCount {
                expected: self.cap,
                actual: self.ledger.active_count(),
            });
        }
        Ok(self)
    }

    /// Evaluate and apply a one-for-one replacement.
    ///
    /// The portfolio is consumed so a failed candidate admission cannot expose
    /// the internally finalized predecessor as usable state.
    pub fn assess_replacement(
        mut self,
        admission: &ReplacementAdmission,
    ) -> Result<Self, MoonshotError> {
        let candidate = admission.candidate();
        let displacement = admission.displacement();
        if self.active.contains_key(candidate.bead_id()) {
            return Err(MoonshotError::CandidateAlreadyActive {
                bead_id: candidate.bead_id().to_string(),
            });
        }
        if self.terminal_beads.contains(candidate.bead_id()) {
            return Err(MoonshotError::TerminalBeadRevival {
                bead_id: candidate.bead_id().to_string(),
            });
        }
        let displaced = self.active.get(displacement.bead_id()).ok_or_else(|| {
            MoonshotError::DisplacedBeadNotActive {
                bead_id: displacement.bead_id().to_string(),
            }
        })?;
        let receipt = finalization_receipt(displaced.mechanism, displacement)?;
        self.ledger.finalize(
            &receipt,
            IdempotencyKey::derive(&format!("moonshot-displace:{}", displacement.bead_id())),
        )?;
        self.ledger.admit(
            candidate.charter(),
            candidate.mechanism(),
            candidate.reservation(),
            IdempotencyKey::derive(&format!("moonshot-admit:{}", candidate.bead_id())),
        )?;
        if self.ledger.active_count() != self.cap {
            return Err(MoonshotError::ReplacementChangedActiveCount {
                expected: self.cap,
                actual: self.ledger.active_count(),
            });
        }
        self.active.remove(displacement.bead_id());
        self.terminal_beads
            .insert(displacement.bead_id().to_string());
        self.active.insert(
            candidate.bead_id().to_string(),
            ActiveMoonshot {
                mechanism: candidate.mechanism(),
            },
        );
        self.displacements.push(displacement.clone());
        Ok(self)
    }

    /// Fixed active WIP cap.
    #[must_use]
    pub const fn cap(&self) -> u32 {
        self.cap
    }

    /// Current active count.
    #[must_use]
    pub fn active_count(&self) -> u32 {
        self.ledger.active_count()
    }

    /// Whether a Bead currently owns a slot.
    #[must_use]
    pub fn is_active(&self, bead_id: &str) -> bool {
        self.active.contains_key(bead_id)
    }

    /// Whether a Bead has permanently left this portfolio instance.
    #[must_use]
    pub fn is_terminal(&self, bead_id: &str) -> bool {
        self.terminal_beads.contains(bead_id)
    }

    /// Exact retained displacement history for this portfolio instance.
    #[must_use]
    pub fn displacements(&self) -> &[DisplacementRecord] {
        &self.displacements
    }

    /// Underlying deterministic decision log.
    #[must_use]
    pub const fn ledger(&self) -> &PortfolioLedger {
        &self.ledger
    }
}
