//! Validated events on prescribed analytic motor paths (bead
//! `frankensim-ext-events-validated-prescribed-6b8h`).
//!
//! Doctrine D3: a discrete transition carries a class-valid
//! certificate, a typed possible-event, a set-valued outcome, an
//! explicit refusal, or an Estimated baseline label — a silent event
//! miss is Sev-0. This module implements the PRESCRIBED-PATH rung:
//! for analytic motor tubes, Taylor-model evaluation encloses the
//! exact path directly (no flowpipe), so complete finite root
//! accounting is achievable with fs-ivl primitives.
//!
//! The claim boundary is the guard class:
//!
//! - [`GuardFamily`] (Taylor pair): the caller supplies, per tube
//!   segment, a Taylor model of the guard AND a Taylor model whose
//!   band rigorously encloses the guard's true time derivative. The
//!   built-in [`plane_crossing_guard`] constructs such a pair from a
//!   tube + rate tube. For this class the scan proves complete
//!   accounting: every subinterval of the span ends classified as
//!   root-free (band excludes zero, or certified-monotone with equal
//!   endpoint signs), as exactly one certified crossing
//!   (certified-monotone with opposite endpoint signs), or as an
//!   explicit [`PossibleEvent`] window. Nothing is dropped.
//! - Anything weaker (black-box guards) must use
//!   [`estimated_scan`], whose outputs are LABELED Estimated and are
//!   also the falsifier lane: scans falsify, they never prove.
//!
//! Grazing (guard and derivative bands both containing zero at the
//! resolution floor) yields a [`PossibleEvent`] with the joint
//! enclosure — Unknown IS the correct verdict there; manufacturing a
//! finite certificate from a zero-containing interval alone would be
//! the Sev-0 failure this module exists to prevent.

use crate::MotionError;
use crate::algebra::{
    POINT_BLADE_NX, POINT_BLADE_NZ, POINT_BLADE_W, POINT_BLADE_Y, TmMv, point_to_mv,
};
use crate::analytic::MotorRateTube;
use crate::tube::CertifiedMotorTube;
use fs_exec::Cx;
use fs_geom::Point3;
use fs_ivl::{Interval, TaylorModel1};

/// One segment's guard pair. Class invariant (the caller's
/// obligation, discharged by construction for the built-in guards):
/// `g` encloses one fixed real function on the segment domain and
/// `gdot` encloses that function's true derivative.
#[derive(Debug, Clone)]
pub struct GuardModel {
    g: TaylorModel1,
    gdot: TaylorModel1,
}

impl GuardModel {
    /// Assemble a pair from caller-built models (domains must match).
    pub fn new(g: TaylorModel1, gdot: TaylorModel1) -> Result<GuardModel, MotionError> {
        if g.domain() != gdot.domain() {
            return Err(MotionError::InvalidConfiguration {
                what: "guard and derivative models must share one domain",
            });
        }
        Ok(GuardModel { g, gdot })
    }

    /// The guard model.
    #[must_use]
    pub fn guard(&self) -> &TaylorModel1 {
        &self.g
    }

    /// The derivative model.
    #[must_use]
    pub fn derivative(&self) -> &TaylorModel1 {
        &self.gdot
    }
}

/// A guard family aligned with a tube's segments.
#[derive(Debug, Clone)]
pub struct GuardFamily {
    segments: Vec<GuardModel>,
}

impl GuardFamily {
    /// Assemble from per-segment pairs (at least one).
    pub fn new(segments: Vec<GuardModel>) -> Result<GuardFamily, MotionError> {
        if segments.is_empty() {
            return Err(MotionError::EmptyTimeDomain);
        }
        Ok(GuardFamily { segments })
    }

    /// The per-segment pairs.
    #[must_use]
    pub fn segments(&self) -> &[GuardModel] {
        &self.segments
    }
}

/// Crossing direction of a certified event (sign of the certified
/// derivative band over the event window).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossingDirection {
    /// Guard goes negative → positive.
    Rising,
    /// Guard goes positive → negative.
    Falling,
}

/// A certified crossing: EXACTLY ONE root of the true guard lies in
/// `window` (certified-monotone derivative band + rigorous opposite
/// endpoint signs).
#[derive(Debug, Clone, Copy)]
pub struct CertifiedEvent {
    /// Time window containing exactly one root.
    pub window: Interval,
    /// Crossing direction.
    pub direction: CrossingDirection,
    /// Guard band over the window (context for consumers).
    pub guard_band: Interval,
}

/// Why a window stayed possible instead of certified.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PossibleReason {
    /// Guard and derivative bands both contain zero at the resolution
    /// floor: tangency/grazing cannot be excluded.
    Grazing,
    /// An endpoint sign could not be certified at the resolution
    /// floor (root at or near the sample point).
    ResolutionFloor,
    /// The subdivision budget ran out before this window was
    /// classified.
    BudgetExhausted,
}

/// A window where an event can be neither certified nor excluded.
/// Unknown is a first-class verdict, never a silent drop.
#[derive(Debug, Clone, Copy)]
pub struct PossibleEvent {
    /// The unresolved time window.
    pub window: Interval,
    /// Guard band over the window.
    pub guard_band: Interval,
    /// Derivative band over the window.
    pub derivative_band: Interval,
    /// Why it stayed possible.
    pub reason: PossibleReason,
}

/// Overall verdict of a scan.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScanVerdict {
    /// Every subinterval classified; the confirmed count is THE count.
    Complete,
    /// All work finished but possible-event windows remain. `confirmed`
    /// is a certified lower bound; without additional evidence, an
    /// unresolved window has no finite certified root-count upper bound.
    IncompleteUnknownWindows,
    /// The subdivision budget ran out; unclassified remainders were
    /// converted to possible-event windows.
    SubdivisionBudgetExhausted,
    /// The certified-event budget was reached while unclassified work
    /// remained; the untouched tail is surfaced as possible windows.
    ZenoBudgetExceeded,
}

/// Root accounting with a certified count when complete, or a certified
/// lower bound plus explicit unresolved windows otherwise.
#[derive(Debug, Clone, Copy)]
pub struct RootCountCertificate {
    /// Certified crossings (each window holds exactly one root).
    pub confirmed: usize,
    /// Possible-event windows (each holds zero or more roots).
    pub possible_windows: usize,
    /// The verdict; `Complete` iff `possible_windows == 0` and no
    /// budget tripped.
    pub verdict: ScanVerdict,
}

/// PO-4 discipline: what the scan actually did.
#[derive(Debug, Clone, Copy)]
pub struct ScanReceipt {
    /// Intervals examined.
    pub intervals_examined: usize,
    /// Deepest bisection level reached.
    pub max_depth: usize,
    /// Leaves proven root-free.
    pub excluded_leaves: usize,
    /// Widest guard band observed at a leaf (enclosure-width logging).
    pub widest_leaf_guard_band: f64,
}

/// A completed scan over one guard family and span.
#[derive(Debug, Clone)]
pub struct EventScan {
    /// Certified crossings in ascending window order.
    pub certified: Vec<CertifiedEvent>,
    /// Possible-event windows in ascending window order.
    pub possible: Vec<PossibleEvent>,
    /// Certified root accounting and unresolved-window count.
    pub count: RootCountCertificate,
    /// Work receipt.
    pub receipt: ScanReceipt,
}

/// Scan budgets. Every exhaustion is a typed, visible outcome.
#[derive(Debug, Clone, Copy)]
pub struct EventScanConfig {
    /// Windows narrower than this stop subdividing and become
    /// possible events when unresolved.
    pub min_width: f64,
    /// Certified windows are bisection-refined (certificate-preserving)
    /// down to roughly this width, budget permitting.
    pub refine_width: f64,
    /// Total interval-examination budget.
    pub max_subdivisions: usize,
    /// Zeno guard: maximum certified events before the scan refuses
    /// to continue.
    pub max_certified_events: usize,
}

impl Default for EventScanConfig {
    fn default() -> Self {
        EventScanConfig {
            min_width: 1e-9,
            refine_width: 1e-6,
            max_subdivisions: 65_536,
            max_certified_events: 4_096,
        }
    }
}

fn contains_zero(iv: Interval) -> bool {
    iv.lo() <= 0.0 && iv.hi() >= 0.0
}

enum LeafClass {
    Excluded,
    Certified(CertifiedEvent),
    Possible(PossibleEvent),
    Split,
}

fn validate_guard_alignment(
    tube: &CertifiedMotorTube,
    family: &GuardFamily,
) -> Result<(), MotionError> {
    if family.segments.len() != tube.segments().len() {
        return Err(MotionError::InvalidConfiguration {
            what: "guard family must align with the tube's segments",
        });
    }
    if tube
        .segments()
        .iter()
        .zip(&family.segments)
        .any(|(segment, model)| segment.domain() != model.g.domain())
    {
        return Err(MotionError::InvalidConfiguration {
            what: "each guard model domain must equal its tube segment domain",
        });
    }
    Ok(())
}

fn classify_leaf(model: &GuardModel, window: Interval, min_width: f64) -> LeafClass {
    let guard_band = model.g.eval_interval(window);
    if !contains_zero(guard_band) {
        return LeafClass::Excluded;
    }
    let derivative_band = model.gdot.eval_interval(window);
    let at_floor = window.width() <= min_width;
    if contains_zero(derivative_band) {
        // Possible tangency; only the resolution floor stops us.
        if at_floor {
            return LeafClass::Possible(PossibleEvent {
                window,
                guard_band,
                derivative_band,
                reason: PossibleReason::Grazing,
            });
        }
        return LeafClass::Split;
    }
    // Certified-monotone window: decide by endpoint signs.
    let at_lo = model.g.eval_interval(Interval::point(window.lo()));
    let at_hi = model.g.eval_interval(Interval::point(window.hi()));
    if contains_zero(at_lo) || contains_zero(at_hi) {
        if at_floor {
            return LeafClass::Possible(PossibleEvent {
                window,
                guard_band,
                derivative_band,
                reason: PossibleReason::ResolutionFloor,
            });
        }
        return LeafClass::Split;
    }
    let lo_negative = at_lo.hi() < 0.0;
    let hi_negative = at_hi.hi() < 0.0;
    if lo_negative == hi_negative {
        // Strictly monotone with equal endpoint signs: root-free.
        return LeafClass::Excluded;
    }
    let direction = if derivative_band.lo() > 0.0 {
        CrossingDirection::Rising
    } else {
        CrossingDirection::Falling
    };
    LeafClass::Certified(CertifiedEvent {
        window,
        direction,
        guard_band,
    })
}

/// Certificate-preserving bisection refinement of a certified window:
/// the derivative band stays sign-definite on every subwindow of the
/// parent, and each kept half retains rigorously opposite endpoint
/// signs. Stops at the refine target, an undecidable midpoint sign, a
/// bisection floor, or budget exhaustion — the window simply stays
/// wider; the certificate is never weakened.
fn refine_certified(
    model: &GuardModel,
    event: CertifiedEvent,
    config: &EventScanConfig,
    receipt: &mut ScanReceipt,
) -> CertifiedEvent {
    let mut window = event.window;
    while window.width() > config.refine_width
        && receipt.intervals_examined < config.max_subdivisions
    {
        let mid = window.midpoint();
        if !(mid > window.lo() && mid < window.hi()) {
            break;
        }
        receipt.intervals_examined += 1;
        let at_mid = model.g.eval_interval(Interval::point(mid));
        if contains_zero(at_mid) {
            break;
        }
        let mid_negative = at_mid.hi() < 0.0;
        let root_in_left = match event.direction {
            // Rising: negative before the root, positive after.
            CrossingDirection::Rising => !mid_negative,
            // Falling: positive before the root, negative after.
            CrossingDirection::Falling => mid_negative,
        };
        window = if root_in_left {
            Interval::new(window.lo(), mid)
        } else {
            Interval::new(mid, window.hi())
        };
    }
    CertifiedEvent {
        window,
        direction: event.direction,
        guard_band: model.g.eval_interval(window),
    }
}

/// Scan one guard family over `span` with complete accounting for the
/// Taylor-pair class. Deterministic: fixed left-first bisection order,
/// no scheduler dependence. Polls cancellation at every interval pop.
pub fn scan_events(
    tube: &CertifiedMotorTube,
    family: &GuardFamily,
    span: Interval,
    config: &EventScanConfig,
    cx: &Cx<'_>,
) -> Result<EventScan, MotionError> {
    if !(config.min_width > 0.0 && config.min_width.is_finite())
        || !(config.refine_width > 0.0 && config.refine_width.is_finite())
    {
        return Err(MotionError::InvalidConfiguration {
            what: "min_width and refine_width must be positive and finite",
        });
    }
    if config.max_subdivisions == 0 || config.max_certified_events == 0 {
        return Err(MotionError::InvalidConfiguration {
            what: "subdivision and event budgets must be positive",
        });
    }
    validate_guard_alignment(tube, family)?;
    let domain = tube.domain();
    if !domain.encloses(span) {
        return Err(MotionError::OutOfDomain {
            lo: span.lo(),
            hi: span.hi(),
            domain_lo: domain.lo(),
            domain_hi: domain.hi(),
        });
    }
    let mut certified = Vec::new();
    let mut possible = Vec::new();
    let mut receipt = ScanReceipt {
        intervals_examined: 0,
        max_depth: 0,
        excluded_leaves: 0,
        widest_leaf_guard_band: 0.0,
    };
    let mut verdict = ScanVerdict::Complete;

    let mut zeno_budget_hit = false;
    'segments: for (segment_index, (segment, model)) in
        tube.segments().iter().zip(&family.segments).enumerate()
    {
        let d = segment.domain();
        let lo = span.lo().max(d.lo());
        let hi = span.hi().min(d.hi());
        if lo >= hi {
            continue;
        }
        if zeno_budget_hit {
            cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
            let window = Interval::new(lo, hi);
            possible.push(PossibleEvent {
                window,
                guard_band: model.g.eval_interval(window),
                derivative_band: model.gdot.eval_interval(window),
                reason: PossibleReason::BudgetExhausted,
            });
            continue;
        }
        // LIFO stack, right pushed first so the left half is examined
        // first: ascending-time deterministic order.
        let mut stack: Vec<(Interval, usize)> = vec![(Interval::new(lo, hi), 0)];
        while let Some((window, depth)) = stack.pop() {
            cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
            if receipt.intervals_examined >= config.max_subdivisions {
                // Budget exhausted: every unclassified window (this one
                // and the rest of the stack) becomes a possible event.
                verdict = ScanVerdict::SubdivisionBudgetExhausted;
                let mut leftovers = vec![window];
                leftovers.extend(stack.drain(..).map(|(w, _)| w));
                for w in leftovers {
                    possible.push(PossibleEvent {
                        window: w,
                        guard_band: model.g.eval_interval(w),
                        derivative_band: model.gdot.eval_interval(w),
                        reason: PossibleReason::BudgetExhausted,
                    });
                }
                continue 'segments;
            }
            receipt.intervals_examined += 1;
            receipt.max_depth = receipt.max_depth.max(depth);
            match classify_leaf(model, window, config.min_width) {
                LeafClass::Excluded => {
                    receipt.excluded_leaves += 1;
                    let width = model.g.eval_interval(window).width();
                    receipt.widest_leaf_guard_band = receipt.widest_leaf_guard_band.max(width);
                }
                LeafClass::Certified(event) => {
                    receipt.widest_leaf_guard_band =
                        receipt.widest_leaf_guard_band.max(event.guard_band.width());
                    certified.push(refine_certified(model, event, config, &mut receipt));
                    if certified.len() >= config.max_certified_events {
                        let future_segment_intersects_span =
                            tube.segments()[segment_index + 1..].iter().any(|future| {
                                let future_domain = future.domain();
                                span.lo().max(future_domain.lo())
                                    < span.hi().min(future_domain.hi())
                            });
                        if !stack.is_empty() || future_segment_intersects_span {
                            // Zeno guard: stop and surface the untouched
                            // current-segment remainder plus every later
                            // segment as possible windows.
                            verdict = ScanVerdict::ZenoBudgetExceeded;
                            for (w, _) in stack.drain(..) {
                                possible.push(PossibleEvent {
                                    window: w,
                                    guard_band: model.g.eval_interval(w),
                                    derivative_band: model.gdot.eval_interval(w),
                                    reason: PossibleReason::BudgetExhausted,
                                });
                            }
                            zeno_budget_hit = true;
                            continue 'segments;
                        }
                    }
                }
                LeafClass::Possible(event) => {
                    receipt.widest_leaf_guard_band =
                        receipt.widest_leaf_guard_band.max(event.guard_band.width());
                    possible.push(event);
                }
                LeafClass::Split => {
                    let mid = window.midpoint();
                    if !(mid > window.lo() && mid < window.hi()) {
                        // Cannot bisect further in f64: resolution floor.
                        let event = PossibleEvent {
                            window,
                            guard_band: model.g.eval_interval(window),
                            derivative_band: model.gdot.eval_interval(window),
                            reason: PossibleReason::ResolutionFloor,
                        };
                        receipt.widest_leaf_guard_band =
                            receipt.widest_leaf_guard_band.max(event.guard_band.width());
                        possible.push(event);
                        continue;
                    }
                    stack.push((Interval::new(mid, window.hi()), depth + 1));
                    stack.push((Interval::new(window.lo(), mid), depth + 1));
                }
            }
        }
    }
    certified.sort_by(|a, b| a.window.lo().total_cmp(&b.window.lo()));
    possible.sort_by(|a, b| a.window.lo().total_cmp(&b.window.lo()));
    if verdict == ScanVerdict::Complete && !possible.is_empty() {
        verdict = ScanVerdict::IncompleteUnknownWindows;
    }
    Ok(EventScan {
        count: RootCountCertificate {
            confirmed: certified.len(),
            possible_windows: possible.len(),
            verdict,
        },
        certified,
        possible,
        receipt,
    })
}

/// Build a plane-crossing guard family for a moving point: the guard
/// is `⟨M(t)·x, n⟩ − offset`, represented (roots unchanged) by the
/// weight-cleared numerator `n·p_num(t) − offset·w(t)` after
/// certifying the homogeneous weight band strictly positive on every
/// segment. The derivative model is built from the rate tube by the
/// sandwich product rule, so its band rigorously encloses the true
/// derivative of the SAME real function.
pub fn plane_crossing_guard(
    tube: &CertifiedMotorTube,
    rate: &MotorRateTube,
    point: Point3,
    normal: [f64; 3],
    offset: f64,
) -> Result<GuardFamily, MotionError> {
    if !(point.x.is_finite() && point.y.is_finite() && point.z.is_finite()) {
        return Err(MotionError::NonFiniteInput { what: "point" });
    }
    if !normal.iter().all(|v| v.is_finite()) || !offset.is_finite() {
        return Err(MotionError::NonFiniteInput { what: "plane" });
    }
    if tube.segments().len() != rate.segments().len() {
        return Err(MotionError::InvalidConfiguration {
            what: "rate tube must align with the primal tube's segments",
        });
    }
    let mut models = Vec::with_capacity(tube.segments().len());
    for (segment, rate_mv) in tube.segments().iter().zip(rate.segments()) {
        let mv = segment.components();
        let p = TmMv::constant(
            &point_to_mv(point.x, point.y, point.z),
            mv.domain(),
            mv.order(),
        )?;
        let rev = mv.reverse()?;
        let sandwich = mv.gp(&p)?.gp(&rev)?;
        // Weight positivity makes numerator roots equal to the real
        // guard's roots.
        let w_band = sandwich.component(POINT_BLADE_W).bound();
        if w_band.lo() <= 0.0 {
            return Err(MotionError::DegenerateWeight {
                lo: w_band.lo(),
                hi: w_band.hi(),
            });
        }
        // d(M P M̃)/dt = M' P M̃ + M P M̃'  (P constant in t).
        let rate_rev = rate_mv.reverse()?;
        let d_sandwich = rate_mv
            .gp(&p)?
            .gp(&rev)?
            .add_componentwise(&mv.gp(&p)?.gp(&rate_rev)?)?;
        let combine = |s: &TmMv| -> Result<TaylorModel1, MotionError> {
            // x = −S[13]/w, y = S[11]/w, z = −S[7]/w ⇒ the numerator of
            // ⟨p, n⟩ − offset (times w) is
            //   −n_x·S[13] + n_y·S[11] − n_z·S[7] − offset·S[14].
            let mut acc = s.component(POINT_BLADE_NX).scale(-normal[0])?;
            acc = acc.try_add(&s.component(POINT_BLADE_Y).scale(normal[1])?)?;
            acc = acc.try_add(&s.component(POINT_BLADE_NZ).scale(-normal[2])?)?;
            acc = acc.try_add(&s.component(POINT_BLADE_W).scale(-offset)?)?;
            Ok(acc)
        };
        models.push(GuardModel::new(combine(&sandwich)?, combine(&d_sandwich)?)?);
    }
    GuardFamily::new(models)
}

/// An Estimated event from the classical dense-scan lane.
#[derive(Debug, Clone, Copy)]
pub struct EstimatedEvent {
    /// Bracketing sample times (sign change observed between them).
    pub bracket: Interval,
    /// Direction implied by the two samples.
    pub direction: CrossingDirection,
}

/// The classical dense-output lane: pointwise midpoint samples of the
/// guard at uniform resolution. Outputs are ESTIMATED — a sampling
/// scan can miss events between samples and proves nothing. Its role
/// is (a) a fast baseline and (b) the independent falsifier for the
/// certified lane: every estimated event must land inside a certified
/// or possible window.
pub fn estimated_scan(
    tube: &CertifiedMotorTube,
    family: &GuardFamily,
    span: Interval,
    samples_per_segment: usize,
    cx: &Cx<'_>,
) -> Result<Vec<EstimatedEvent>, MotionError> {
    if samples_per_segment < 2 {
        return Err(MotionError::InvalidConfiguration {
            what: "estimated scan needs at least two samples per segment",
        });
    }
    validate_guard_alignment(tube, family)?;
    let mut events = Vec::new();
    for (segment, model) in tube.segments().iter().zip(&family.segments) {
        let d = segment.domain();
        let lo = span.lo().max(d.lo());
        let hi = span.hi().min(d.hi());
        if lo >= hi {
            continue;
        }
        let n = samples_per_segment;
        let mut prev_t = lo;
        let mut prev_g = model.g.eval_interval(Interval::point(lo)).midpoint();
        for k in 1..=n {
            cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
            let t = lo + (hi - lo) * (k as f64 / n as f64);
            let g = model.g.eval_interval(Interval::point(t)).midpoint();
            if prev_g != 0.0 && g != 0.0 && (prev_g < 0.0) != (g < 0.0) {
                events.push(EstimatedEvent {
                    bracket: Interval::new(prev_t, t),
                    direction: if g > prev_g {
                        CrossingDirection::Rising
                    } else {
                        CrossingDirection::Falling
                    },
                });
            }
            prev_t = t;
            prev_g = g;
        }
    }
    Ok(events)
}

/// A group of certified events (from different guards) whose windows
/// overlap: the true event order inside the group is undetermined, so
/// admissible orders are enumerated instead of picked.
#[derive(Debug, Clone)]
pub struct SimultaneousGroup {
    /// (guard index, event) members, ascending by window start.
    pub members: Vec<(usize, CertifiedEvent)>,
    /// Admissible orderings (permutations of member indices) when the
    /// group is small enough to enumerate; `None` when enumeration was
    /// refused (too many members) — the group is then explicitly
    /// set-valued and unordered.
    pub admissible_orders: Option<Vec<Vec<usize>>>,
}

/// Maximum group size for explicit order enumeration (4! = 24).
pub const MAX_ENUMERATED_GROUP: usize = 4;

fn permutations(n: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut items: Vec<usize> = (0..n).collect();
    fn heap(k: usize, items: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if k <= 1 {
            out.push(items.clone());
            return;
        }
        for i in 0..k {
            heap(k - 1, items, out);
            if k % 2 == 0 {
                items.swap(i, k - 1);
            } else {
                items.swap(0, k - 1);
            }
        }
    }
    heap(items.len(), &mut items, &mut out);
    out
}

fn order_respects_certified_precedence(
    order: &[usize],
    members: &[(usize, CertifiedEvent)],
) -> bool {
    let mut position = vec![0usize; order.len()];
    for (ordinal, &member) in order.iter().enumerate() {
        position[member] = ordinal;
    }
    for earlier in 0..members.len() {
        for later in 0..members.len() {
            if members[earlier].1.window.hi() < members[later].1.window.lo()
                && position[earlier] >= position[later]
            {
                return false;
            }
        }
    }
    true
}

/// Group certified events from multiple guard scans by window overlap
/// and enumerate admissible orders. Events whose windows are disjoint
/// are totally ordered by time; overlapping windows form set-valued
/// groups.
#[must_use]
pub fn enumerate_simultaneous(scans: &[EventScan]) -> Vec<SimultaneousGroup> {
    let mut all: Vec<(usize, CertifiedEvent)> = Vec::new();
    for (gi, scan) in scans.iter().enumerate() {
        for e in &scan.certified {
            all.push((gi, *e));
        }
    }
    all.sort_by(|a, b| {
        a.1.window
            .lo()
            .total_cmp(&b.1.window.lo())
            .then(a.0.cmp(&b.0))
    });
    let mut groups: Vec<SimultaneousGroup> = Vec::new();
    for (gi, e) in all {
        let overlaps = groups.last().is_some_and(|g| {
            g.members
                .iter()
                .any(|(_, m)| e.window.lo() <= m.window.hi() && m.window.lo() <= e.window.hi())
        });
        if overlaps {
            let group = groups.last_mut().expect("non-empty by overlap check");
            group.members.push((gi, e));
        } else {
            groups.push(SimultaneousGroup {
                members: vec![(gi, e)],
                admissible_orders: None,
            });
        }
    }
    for group in &mut groups {
        if group.members.len() <= MAX_ENUMERATED_GROUP {
            group.admissible_orders = Some(if group.members.len() == 1 {
                vec![vec![0]]
            } else {
                permutations(group.members.len())
                    .into_iter()
                    .filter(|order| order_respects_certified_precedence(order, &group.members))
                    .collect()
            });
        }
    }
    groups
}
