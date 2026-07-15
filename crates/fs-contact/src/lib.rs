//! fs-contact — capability-routed body-to-body contact (bead tqag,
//! Stage 1 / increment 1).
//!
//! Layer: L3. Blocker B3's opening move: bodies are geometry bound to
//! CERTIFIED motion ([`fs_motion::CertifiedMotorTube`]), the broad
//! phase prunes with conservative spacetime supports (an AABB
//! containing the body's image over the WHOLE query window — no
//! sampled instants, no tunneling through the gaps between samples),
//! and the narrow phase is capability-routed: a pairing either
//! carries the query theorem it needs (convex support maps, exact
//! SDFs) or REFUSES with a typed error naming the missing capability.
//! Candidate management runs under an explicit budget (program risk
//! #2): exhaustion returns the unresolved pairs, never a truncated
//! set presented as complete.
//!
//! Increment boundary (recorded in CONTRACT no-claims): certified CCD
//! (feature-pair conservative advancement over tubes — deliberately
//! NOT a simple-root guard on a global separation function),
//! EPA-class penetration certificates, Rep Router error inflation,
//! and the tube-source-agnostic trajectory interface land in later
//! increments of this bead's staging plan.

use fs_exec::Cx;
use fs_geom::Aabb;
use fs_ivl::Interval;
use fs_motion::{CertifiedMotorTube, MotionError};
use fs_query::{ConvexSeparation, ConvexSupportMap, QueryError, convex_separation};

/// Hard bound on bodies per broad-phase call.
pub const MAX_CONTACT_BODIES: usize = 1 << 16;

/// Teaching errors (P10): every refusal names the violated assumption.
#[derive(Debug)]
pub enum ContactError {
    /// A motion enclosure refused (domain, chart transition, budget…).
    Motion(MotionError),
    /// A geometry query refused (capability, evidence, budget…).
    Query(QueryError),
    /// Too many bodies for the deterministic broad phase.
    TooManyBodies {
        /// Bodies supplied.
        bodies: usize,
        /// Public ceiling.
        max: usize,
    },
    /// A body's body-frame support box is non-finite or inverted.
    InvalidSupport {
        /// The offending body index.
        body: usize,
    },
    /// The query window is empty, inverted, or non-finite.
    InvalidWindow {
        /// Window low endpoint.
        lo: f64,
        /// Window high endpoint.
        hi: f64,
    },
    /// The candidate budget is exhausted. The resolved prefix is
    /// sound but INCOMPLETE; the unresolved pairs are listed so the
    /// caller can split the window or raise the budget — a truncated
    /// set is never presented as the full candidate set.
    CandidateBudgetExhausted {
        /// The caller's cap.
        max_pairs: usize,
        /// Overlapping pairs that no longer fit the budget.
        unresolved: Vec<(usize, usize)>,
    },
    /// The pairing's narrow-phase route needs a capability neither
    /// side declared. Refusal, not guessing, per the routing doctrine.
    MissingCapability {
        /// First body index of the pair.
        body_a: usize,
        /// Second body index of the pair.
        body_b: usize,
        /// The absent capability, by stable name.
        capability: &'static str,
    },
    /// Cancelled mid-scan.
    Cancelled,
}

impl core::fmt::Display for ContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ContactError::Motion(e) => write!(f, "contact motion enclosure refused: {e}"),
            ContactError::Query(e) => write!(f, "contact geometry query refused: {e}"),
            ContactError::TooManyBodies { bodies, max } => write!(
                f,
                "{bodies} bodies exceed the deterministic {max}-body broad-phase ceiling"
            ),
            ContactError::InvalidSupport { body } => write!(
                f,
                "body {body} has a non-finite or inverted body-frame support box"
            ),
            ContactError::InvalidWindow { lo, hi } => {
                write!(f, "contact window [{lo}, {hi}] must be finite with lo < hi")
            }
            ContactError::CandidateBudgetExhausted {
                max_pairs,
                unresolved,
            } => write!(
                f,
                "candidate budget {max_pairs} exhausted with {} unresolved pairs; split \
                 the window or raise the budget — the resolved prefix is not the full set",
                unresolved.len()
            ),
            ContactError::MissingCapability {
                body_a,
                body_b,
                capability,
            } => write!(
                f,
                "pair ({body_a}, {body_b}) needs the {capability:?} capability that neither \
                 side declares; refusing rather than guessing"
            ),
            ContactError::Cancelled => write!(f, "cancelled mid-contact-query"),
        }
    }
}

impl core::error::Error for ContactError {}

impl From<MotionError> for ContactError {
    fn from(e: MotionError) -> Self {
        ContactError::Motion(e)
    }
}

impl From<QueryError> for ContactError {
    fn from(e: QueryError) -> Self {
        ContactError::Query(e)
    }
}

/// One body in the broad phase: a body-frame support box bound to a
/// certified body-to-world motor tube.
pub struct SpacetimeBody<'a> {
    support: Aabb,
    tube: &'a CertifiedMotorTube,
}

impl<'a> SpacetimeBody<'a> {
    /// Bind a finite body-frame support box to a tube.
    ///
    /// # Errors
    /// [`ContactError::InvalidSupport`] (reported with body index 0;
    /// the broad phase re-reports with the true index) for non-finite
    /// or inverted boxes.
    pub fn new(support: Aabb, tube: &'a CertifiedMotorTube) -> Result<Self, ContactError> {
        let finite = [
            support.min.x,
            support.min.y,
            support.min.z,
            support.max.x,
            support.max.y,
            support.max.z,
        ]
        .iter()
        .all(|v| v.is_finite())
            && support.min.x <= support.max.x
            && support.min.y <= support.max.y
            && support.min.z <= support.max.z;
        if finite {
            Ok(SpacetimeBody { support, tube })
        } else {
            Err(ContactError::InvalidSupport { body: 0 })
        }
    }

    /// The body-frame support box.
    #[must_use]
    pub fn support(&self) -> &Aabb {
        &self.support
    }

    /// The bound tube.
    #[must_use]
    pub fn tube(&self) -> &CertifiedMotorTube {
        &self.tube
    }
}

/// Broad-phase output: candidate pairs plus honest statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct BroadPhaseReport {
    /// Lexicographically sorted candidate pairs `(i, j)` with `i < j`.
    pub pairs: Vec<(usize, usize)>,
    /// Pairs whose windowed boxes were compared.
    pub checked_pairs: usize,
    /// Pairs pruned by the sweep before any box comparison.
    pub pruned_pairs: usize,
    /// Worst versor-defect bound among the windowed enclosures — the
    /// motion-side quality signal consumers must carry forward.
    pub max_defect: f64,
}

/// Conservative spacetime broad phase over one shared window.
///
/// Each body's windowed world box comes from
/// [`CertifiedMotorTube::box_action_over`] — an enclosure of the
/// body's image for EVERY `t` in the window, so a pair whose boxes do
/// not overlap provably cannot touch inside the window (no sampling,
/// no tunneling). Sweep axis: world `x`, sorted by `total_cmp` with
/// index tie-breaks — output is a pure function of the inputs.
///
/// # Errors
/// Window/support/body-count refusals; motion enclosure refusals;
/// [`ContactError::CandidateBudgetExhausted`] listing every
/// unresolved overlapping pair beyond `max_pairs`;
/// [`ContactError::Cancelled`] (checked per body and per sweep
/// stride).
pub fn spacetime_candidates(
    bodies: &[SpacetimeBody<'_>],
    window: Interval,
    max_pairs: usize,
    cx: &Cx<'_>,
) -> Result<BroadPhaseReport, ContactError> {
    if bodies.len() > MAX_CONTACT_BODIES {
        return Err(ContactError::TooManyBodies {
            bodies: bodies.len(),
            max: MAX_CONTACT_BODIES,
        });
    }
    if !(window.lo().is_finite() && window.hi().is_finite() && window.lo() < window.hi()) {
        return Err(ContactError::InvalidWindow {
            lo: window.lo(),
            hi: window.hi(),
        });
    }
    let mut boxes = Vec::with_capacity(bodies.len());
    let mut max_defect = 0.0f64;
    for body in bodies {
        if cx.checkpoint().is_err() {
            return Err(ContactError::Cancelled);
        }
        let enclosure = body.tube.box_action_over(&body.support, window, cx)?;
        max_defect = max_defect.max(enclosure.defect);
        boxes.push(enclosure.bounds);
    }
    // Sweep and prune on x: sort by min.x, scan forward while the
    // next candidate's min.x is at most the current box's max.x.
    let mut order: Vec<usize> = (0..boxes.len()).collect();
    order.sort_unstable_by(|&a, &b| {
        boxes[a]
            .min
            .x
            .total_cmp(&boxes[b].min.x)
            .then_with(|| a.cmp(&b))
    });
    let mut pairs = Vec::new();
    let mut unresolved = Vec::new();
    let mut checked_pairs = 0usize;
    let mut pruned_pairs = 0usize;
    for (rank, &i) in order.iter().enumerate() {
        if cx.checkpoint().is_err() {
            return Err(ContactError::Cancelled);
        }
        let mut row_checked = 0usize;
        for &j in order.iter().skip(rank + 1) {
            if boxes[j].min.x > boxes[i].max.x {
                // Sorted by min.x: every later body starts even
                // further right; the rest of this row is pruned.
                pruned_pairs += order.len() - rank - 1 - row_checked;
                break;
            }
            row_checked += 1;
            checked_pairs += 1;
            let overlap_yz = boxes[i].min.y <= boxes[j].max.y
                && boxes[j].min.y <= boxes[i].max.y
                && boxes[i].min.z <= boxes[j].max.z
                && boxes[j].min.z <= boxes[i].max.z;
            if overlap_yz {
                let pair = (i.min(j), i.max(j));
                if pairs.len() < max_pairs {
                    pairs.push(pair);
                } else {
                    unresolved.push(pair);
                }
            }
        }
    }
    if !unresolved.is_empty() {
        unresolved.sort_unstable();
        return Err(ContactError::CandidateBudgetExhausted {
            max_pairs,
            unresolved,
        });
    }
    pairs.sort_unstable();
    Ok(BroadPhaseReport {
        pairs,
        checked_pairs,
        pruned_pairs,
        max_defect,
    })
}

/// The narrow-phase capability a body declares for a frozen-time
/// query. Routing NEVER guesses: a pairing without a compatible
/// declared route refuses.
pub enum NarrowRoute<'a> {
    /// A world-frame convex support map valid at the query time.
    Convex(&'a dyn ConvexSupportMap),
    /// No narrow-phase capability declared.
    Undeclared,
}

/// One routed narrow-phase verdict.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NarrowVerdict {
    /// Certified convex separation enclosure (fs-query semantics:
    /// `separation_proven ⇔ lo > 0`; overlap is never claimed).
    Convex(ConvexSeparation),
}

/// Route one candidate pair through its declared capabilities at a
/// frozen query time.
///
/// # Errors
/// [`ContactError::MissingCapability`] when either side is
/// [`NarrowRoute::Undeclared`] (stable name `"convex-support-map"` for
/// the only Stage-1 route); refusals from the underlying certified
/// query pass through typed.
pub fn narrow_phase(
    pair: (usize, usize),
    route_a: &NarrowRoute<'_>,
    route_b: &NarrowRoute<'_>,
    max_iterations: u32,
    cx: &Cx<'_>,
) -> Result<NarrowVerdict, ContactError> {
    match (route_a, route_b) {
        (NarrowRoute::Convex(a), NarrowRoute::Convex(b)) => {
            let separation = convex_separation(*a, *b, max_iterations, cx)?;
            Ok(NarrowVerdict::Convex(separation))
        }
        _ => Err(ContactError::MissingCapability {
            body_a: pair.0,
            body_b: pair.1,
            capability: "convex-support-map",
        }),
    }
}
