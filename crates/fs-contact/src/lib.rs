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
    /// The CCD subdivision budget is exhausted. The partial state is
    /// sound but INCOMPLETE: `pending` windows were never examined and
    /// `possible` windows were confirmed unresolved — neither is a
    /// clear verdict, and the caller must widen the budget, loosen the
    /// tolerance, or split the query window.
    CcdBudgetExhausted {
        /// The caller's cap on examined subwindows.
        max_windows: usize,
        /// Subwindows actually examined before exhaustion.
        examined: usize,
        /// Unexamined subwindows, ascending in time.
        pending: Vec<Interval>,
        /// Sub-tolerance windows already confirmed unresolved.
        possible: Vec<Interval>,
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
            ContactError::CcdBudgetExhausted {
                max_windows,
                examined,
                pending,
                possible,
            } => write!(
                f,
                "CCD window budget {max_windows} exhausted after {examined} enclosures with \
                 {} pending and {} unresolved windows; widen the budget, loosen the time \
                 tolerance, or split the query window — this is not a clear verdict",
                pending.len(),
                possible.len()
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
        self.tube
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
        for (row_checked, &j) in order.iter().skip(rank + 1).enumerate() {
            if boxes[j].min.x > boxes[i].max.x {
                // Sorted by min.x: every later body starts even
                // further right; the rest of this row is pruned.
                pruned_pairs += order.len() - rank - 1 - row_checked;
                break;
            }
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

// ── Certified continuous collision detection (bead tqag, increment 2) ──────

/// One pair's certified CCD verdict over a query window.
///
/// Soundness (the Sev-0 no-tunneling claim): every per-window test uses
/// [`CertifiedMotorTube::box_action_over`], an enclosure of the body's
/// image at EVERY instant of the subwindow. A subwindow is declared
/// clear only when the two enclosures are disjoint along a coordinate
/// axis — a proof that no contact exists anywhere inside it. Everything
/// not proven clear is subdivided down to the caller's time tolerance
/// and reported as a possible-contact window, so the union of reported
/// windows CONTAINS every true contact instant. Contact itself is never
/// claimed: box overlap is necessary, not sufficient. This is
/// deliberately NOT a sign-change root guard on a global separation
/// function — persistent or grazing contact has no sign change to find,
/// and the window report stays honest exactly there.
#[derive(Debug, Clone, PartialEq)]
pub enum CcdVerdict {
    /// PROVEN: no contact anywhere in the window.
    ClearWindow {
        /// Certified lower bound on the pair's axis-aligned gap over
        /// the whole window (the smallest separating-axis gap among
        /// the cleared subwindows).
        min_gap: f64,
    },
    /// Time-ordered, disjoint, tolerance-width windows that together
    /// contain every instant at which the pair COULD touch.
    PossibleContact {
        /// The unresolved windows, merged where adjacent.
        windows: Vec<Interval>,
    },
}

/// Certified CCD output: verdict plus honest work statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct CcdReport {
    /// The verdict.
    pub verdict: CcdVerdict,
    /// Subwindows whose enclosures were computed.
    pub examined_windows: usize,
    /// Worst versor-defect bound among every enclosure consulted.
    pub max_defect: f64,
}

/// Certified pairwise CCD by conservative window bisection.
///
/// # Errors
/// Window/tolerance refusals ([`ContactError::InvalidWindow`]); motion
/// enclosure refusals; [`ContactError::CcdBudgetExhausted`] when more
/// than `max_windows` subwindows would be examined (the partial state
/// is returned, never presented as complete);
/// [`ContactError::Cancelled`] (checked per subwindow).
pub fn certified_ccd(
    a: &SpacetimeBody<'_>,
    b: &SpacetimeBody<'_>,
    window: Interval,
    time_tolerance: f64,
    max_windows: usize,
    cx: &Cx<'_>,
) -> Result<CcdReport, ContactError> {
    if !(window.lo().is_finite() && window.hi().is_finite() && window.lo() < window.hi()) {
        return Err(ContactError::InvalidWindow {
            lo: window.lo(),
            hi: window.hi(),
        });
    }
    if !(time_tolerance.is_finite() && time_tolerance > 0.0) {
        return Err(ContactError::InvalidWindow {
            lo: time_tolerance,
            hi: time_tolerance,
        });
    }
    // LIFO with the later half pushed first, so subwindows are examined
    // (and possible windows emitted) in ascending time order.
    let mut stack = vec![window];
    let mut possible: Vec<Interval> = Vec::new();
    let mut examined = 0usize;
    let mut max_defect = 0.0f64;
    let mut min_gap = f64::INFINITY;
    while let Some(w) = stack.pop() {
        if cx.checkpoint().is_err() {
            return Err(ContactError::Cancelled);
        }
        if examined >= max_windows {
            stack.push(w);
            stack.reverse(); // back to ascending time order for the report
            return Err(ContactError::CcdBudgetExhausted {
                max_windows,
                examined,
                pending: stack,
                possible,
            });
        }
        examined += 1;
        let ea = a.tube.box_action_over(&a.support, w, cx)?;
        let eb = b.tube.box_action_over(&b.support, w, cx)?;
        max_defect = max_defect.max(ea.defect).max(eb.defect);
        let (ba, bb) = (&ea.bounds, &eb.bounds);
        let axis_gaps = [
            (bb.min.x - ba.max.x).max(ba.min.x - bb.max.x),
            (bb.min.y - ba.max.y).max(ba.min.y - bb.max.y),
            (bb.min.z - ba.max.z).max(ba.min.z - bb.max.z),
        ];
        let gap = axis_gaps[0].max(axis_gaps[1]).max(axis_gaps[2]);
        if gap > 0.0 {
            // Disjoint along some axis for EVERY instant of `w`: proven
            // clear; the true distance in `w` is at least `gap`.
            min_gap = min_gap.min(gap);
            continue;
        }
        let mid = w.midpoint();
        if w.width() <= time_tolerance || !(mid > w.lo() && mid < w.hi()) {
            // Tolerance reached (or the window can no longer split in
            // f64): report, never claim.
            possible.push(w);
            continue;
        }
        stack.push(Interval::new(mid, w.hi()));
        stack.push(Interval::new(w.lo(), mid));
    }
    if possible.is_empty() {
        return Ok(CcdReport {
            verdict: CcdVerdict::ClearWindow { min_gap },
            examined_windows: examined,
            max_defect,
        });
    }
    // Merge windows that share an endpoint (they arrive time-ordered).
    let mut merged: Vec<Interval> = Vec::with_capacity(possible.len());
    for w in possible {
        match merged.last_mut() {
            Some(last) if last.hi() >= w.lo() => *last = Interval::new(last.lo(), w.hi()),
            _ => merged.push(w),
        }
    }
    Ok(CcdReport {
        verdict: CcdVerdict::PossibleContact { windows: merged },
        examined_windows: examined,
        max_defect,
    })
}

// ── Swept-vertex-hull refinement (bead tqag, increment 3) ──────────────────

use fs_geom::Vec3;

/// Hard bound on vertices per polytope body in the refinement route
/// (each vertex contributes eight trajectory-box corners per window).
pub const MAX_CCD_VERTICES: usize = 1 << 10;

/// The convex hull of a finite corner set, presented as a support map.
/// Corner selection is the trait's documented exact case
/// (`support_slack` = 0): the returned point is always a member corner.
struct SweptVertexHull {
    corners: Vec<fs_geom::Point3>,
}

impl ConvexSupportMap for SweptVertexHull {
    fn support_point(&self, direction: Vec3) -> fs_geom::Point3 {
        // First-strict-max selection: deterministic under permutation of
        // equal dots because corner order is itself deterministic
        // (vertex order × fixed corner order).
        let mut best = self.corners[0];
        let mut best_dot = Vec3::new(best.x, best.y, best.z).dot(direction);
        for corner in &self.corners[1..] {
            let dot = Vec3::new(corner.x, corner.y, corner.z).dot(direction);
            if dot > best_dot {
                best = *corner;
                best_dot = dot;
            }
        }
        best
    }

    fn interior_point(&self) -> fs_geom::Point3 {
        self.corners[0]
    }

    fn support_slack(&self) -> f64 {
        0.0
    }

    fn contained_ball_radius(&self, _center: fs_geom::Point3) -> Option<f64> {
        // A corner cloud proves no interior ball; refusal per the trait.
        None
    }

    fn name(&self) -> &'static str {
        "swept-vertex-hull"
    }
}

/// Build the swept hull of `vertices` under `tube` over `window`: every
/// vertex trajectory is enclosed by [`CertifiedMotorTube::point_action_over`]
/// and the hull of all box corners contains the body's image at every
/// instant (the image of a convex hull under a rigid motion is the hull
/// of the vertex images).
fn swept_vertex_hull(
    vertices: &[fs_geom::Point3],
    tube: &CertifiedMotorTube,
    window: Interval,
    cx: &Cx<'_>,
) -> Result<(SweptVertexHull, f64), ContactError> {
    let mut corners = Vec::with_capacity(vertices.len() * 8);
    let mut defect = 0.0f64;
    for v in vertices {
        if cx.checkpoint().is_err() {
            return Err(ContactError::Cancelled);
        }
        let enc = tube.point_action_over(*v, window, cx)?;
        defect = defect.max(enc.defect);
        let [x, y, z] = enc.coords;
        for &cx_ in &[x.lo(), x.hi()] {
            for &cy in &[y.lo(), y.hi()] {
                for &cz in &[z.lo(), z.hi()] {
                    corners.push(fs_geom::Point3::new(cx_, cy, cz));
                }
            }
        }
    }
    Ok((SweptVertexHull { corners }, defect))
}

/// One refined possible-contact window.
#[derive(Debug, Clone, PartialEq)]
pub enum RefinedWindow {
    /// The swept hulls are PROVEN disjoint over the whole window: the
    /// box verdict was a false alarm here.
    Pruned {
        /// The window.
        window: Interval,
        /// Certified lower bound on the pair's distance over it.
        gap: f64,
    },
    /// The tighter test could not clear the window either; it remains
    /// a possible-contact window.
    Retained {
        /// The window.
        window: Interval,
    },
}

/// Refinement output over a set of possible-contact windows.
#[derive(Debug, Clone, PartialEq)]
pub struct CcdRefinement {
    /// One entry per input window, in input order.
    pub windows: Vec<RefinedWindow>,
    /// Worst versor-defect bound among every enclosure consulted.
    pub max_defect: f64,
}

/// Feature-pair refinement of [`CcdVerdict::PossibleContact`] windows
/// for POLYTOPE bodies (finite vertex sets in body frame): each window
/// is re-tested with the certified separation of the two swept vertex
/// hulls — tight where per-instant axis-aligned boxes are structurally
/// loose (rotated or diagonally-moving bodies), so windows the box
/// verdict could NEVER clear are pruned with a certified gap.
///
/// Soundness: pruning uses `separation_proven` (a certified positive
/// lower bound between supersets of the two swept bodies); a window
/// containing a true contact can therefore never be pruned. Retention
/// claims nothing, exactly like the box verdict.
///
/// # Errors
/// [`ContactError::TooManyBodies`] reusing the vertex ceiling
/// ([`MAX_CCD_VERTICES`]); empty vertex sets refuse as
/// [`ContactError::InvalidSupport`]; motion/query refusals pass
/// through typed; [`ContactError::Cancelled`] per vertex and window.
pub fn refine_possible_windows(
    a_vertices: &[fs_geom::Point3],
    a_tube: &CertifiedMotorTube,
    b_vertices: &[fs_geom::Point3],
    b_tube: &CertifiedMotorTube,
    windows: &[Interval],
    max_iterations: u32,
    cx: &Cx<'_>,
) -> Result<CcdRefinement, ContactError> {
    for (index, vertices) in [(0usize, a_vertices), (1usize, b_vertices)] {
        if vertices.is_empty() {
            return Err(ContactError::InvalidSupport { body: index });
        }
        if vertices.len() > MAX_CCD_VERTICES {
            return Err(ContactError::TooManyBodies {
                bodies: vertices.len(),
                max: MAX_CCD_VERTICES,
            });
        }
    }
    let mut out = Vec::with_capacity(windows.len());
    let mut max_defect = 0.0f64;
    for &window in windows {
        if cx.checkpoint().is_err() {
            return Err(ContactError::Cancelled);
        }
        let (hull_a, defect_a) = swept_vertex_hull(a_vertices, a_tube, window, cx)?;
        let (hull_b, defect_b) = swept_vertex_hull(b_vertices, b_tube, window, cx)?;
        max_defect = max_defect.max(defect_a).max(defect_b);
        let separation = convex_separation(&hull_a, &hull_b, max_iterations, cx)?;
        if separation.separation_proven {
            out.push(RefinedWindow::Pruned {
                window,
                gap: separation.lo,
            });
        } else {
            out.push(RefinedWindow::Retained { window });
        }
    }
    Ok(CcdRefinement {
        windows: out,
        max_defect,
    })
}

// ── SDF-obstacle window pruning (bead tqag, SDF route) ─────────────────────

use fs_geom::{Chart, TraceStepClaim};

/// Refine possible-contact windows for a POLYTOPE body against a STATIC
/// exact-distance obstacle chart (the SDF route of this bead's staging
/// plan).
///
/// Soundness: an exact Euclidean signed-distance function is 1-Lipschitz
/// (the [`TraceStepClaim::ExactDistance`] theorem carried by the chart —
/// charts with weaker claims refuse at entry). The swept body over a
/// window is contained in the ball around any center `c` with radius
/// `r = max |corner − c|` over the swept-vertex-hull corners, so
/// `φ(q) ≥ φ_lo(c) − r` for every swept point `q`; when that bound is
/// positive the whole window is PROVEN clear of the obstacle with a
/// certified gap. The radius computation is inflated multiplicatively
/// (documented in-line) to cover its own f64 rounding; the center choice
/// affects tightness only, never soundness.
///
/// # Errors
/// [`ContactError::MissingCapability`] (stable name
/// `"exact-distance-chart"`) when the obstacle's trace claim is weaker;
/// [`ContactError::InvalidSupport`]/[`ContactError::TooManyBodies`] on
/// vertex-set refusals; motion refusals pass through typed;
/// [`ContactError::Query`] wrapping an unusable chart enclosure;
/// [`ContactError::Cancelled`] per vertex and window.
#[allow(clippy::too_many_lines)] // One conservative bisection transaction, mirroring certified_ccd.
pub fn refine_windows_against_sdf(
    vertices: &[fs_geom::Point3],
    tube: &CertifiedMotorTube,
    obstacle: &dyn Chart,
    windows: &[Interval],
    time_tolerance: f64,
    max_windows: usize,
    cx: &Cx<'_>,
) -> Result<CcdRefinement, ContactError> {
    if obstacle.trace_step_claim() != TraceStepClaim::ExactDistance {
        return Err(ContactError::MissingCapability {
            body_a: 0,
            body_b: 1,
            capability: "exact-distance-chart",
        });
    }
    if vertices.is_empty() {
        return Err(ContactError::InvalidSupport { body: 0 });
    }
    if vertices.len() > MAX_CCD_VERTICES {
        return Err(ContactError::TooManyBodies {
            bodies: vertices.len(),
            max: MAX_CCD_VERTICES,
        });
    }
    if !(time_tolerance.is_finite() && time_tolerance > 0.0) {
        return Err(ContactError::InvalidWindow {
            lo: time_tolerance,
            hi: time_tolerance,
        });
    }
    // The ball around a long sweep is hopelessly loose (the input windows
    // arrive MERGED from certified_ccd), so the route bisects internally:
    // LIFO with the later half pushed first keeps emission time-ordered.
    let mut stack: Vec<Interval> = windows.iter().rev().copied().collect();
    let mut out = Vec::new();
    let mut examined = 0usize;
    let mut max_defect = 0.0f64;
    while let Some(window) = stack.pop() {
        if cx.checkpoint().is_err() {
            return Err(ContactError::Cancelled);
        }
        if examined >= max_windows {
            stack.push(window);
            stack.reverse();
            let possible = out
                .iter()
                .filter_map(|w| match w {
                    RefinedWindow::Retained { window } => Some(*window),
                    RefinedWindow::Pruned { .. } => None,
                })
                .collect();
            return Err(ContactError::CcdBudgetExhausted {
                max_windows,
                examined,
                pending: stack,
                possible,
            });
        }
        examined += 1;
        let (hull, defect) = swept_vertex_hull(vertices, tube, window, cx)?;
        max_defect = max_defect.max(defect);
        // Center: component-wise corner mean. Any finite center is sound;
        // the radius below is measured from THIS center.
        #[allow(clippy::cast_precision_loss)]
        let n = hull.corners.len() as f64;
        let (mut sx, mut sy, mut sz) = (0.0f64, 0.0f64, 0.0f64);
        for corner in &hull.corners {
            sx += corner.x;
            sy += corner.y;
            sz += corner.z;
        }
        let center = fs_geom::Point3::new(sx / n, sy / n, sz / n);
        let mut radius = 0.0f64;
        for corner in &hull.corners {
            let dx = corner.x - center.x;
            let dy = corner.y - center.y;
            let dz = corner.z - center.z;
            radius = radius.max((dx * dx + dy * dy + dz * dz).sqrt());
        }
        // ≤ 6 rounded f64 ops per corner distance at ≤ 1 ulp each: an
        // 8-eps multiplicative inflation strictly covers the rounding.
        let radius = (radius * (1.0 + 8.0 * f64::EPSILON)).next_up();

        let sample = obstacle.eval(center, cx);
        if cx.checkpoint().is_err() {
            return Err(ContactError::Cancelled);
        }
        let enclosure = obstacle.trace_value_enclosure(center, &sample, cx);
        let usable = matches!(
            enclosure.kind,
            fs_evidence::NumericalKind::Exact | fs_evidence::NumericalKind::Enclosure
        ) && enclosure.lo.is_finite()
            && enclosure.hi.is_finite()
            && enclosure.lo <= enclosure.hi;
        if !usable {
            return Err(ContactError::Query(QueryError::InvalidPointSample {
                at: [center.x, center.y, center.z],
            }));
        }
        let gap = (enclosure.lo.next_down() - radius).next_down();
        if gap > 0.0 {
            out.push(RefinedWindow::Pruned { window, gap });
            continue;
        }
        let mid = window.midpoint();
        if window.width() <= time_tolerance || !(mid > window.lo() && mid < window.hi()) {
            out.push(RefinedWindow::Retained { window });
            continue;
        }
        stack.push(Interval::new(mid, window.hi()));
        stack.push(Interval::new(window.lo(), mid));
    }
    Ok(CcdRefinement {
        windows: out,
        max_defect,
    })
}
