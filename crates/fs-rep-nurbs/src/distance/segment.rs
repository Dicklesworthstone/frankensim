//! Continuous segment clearance from outward point-distance queries.

use super::{
    CLOSEST_MAX_SPLITS, MAX_BYTES, MAX_WORK, Meter, SearchRun, SurfaceDistanceBudget, TrimDomain,
    checked, distance_to_box, finite, refusal, search_with_poll,
};
use crate::{AdmittedNurbsSurface, AdmittedTrimmedPatch, NurbsError, NurbsSurface, TrimmedPatch};
use fs_exec::Cx;
use fs_ivl::Interval;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::mem::size_of;

/// Combined limits for the segment frontier and its sequential surface queries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SegmentDistanceBudget {
    /// Binary refinements of the segment parameter interval [0, 1].
    pub max_splits: u32,
    /// Per-point limits, further capped by the remaining combined allowance.
    pub point_budget: SurfaceDistanceBudget,
    /// Aggregate conservative work, capped at the surface query's ceiling.
    pub max_work_units: u128,
    /// Peak payload including the segment frontier and one nested query.
    /// The requested segment frontier is reserved before the first query.
    pub max_retained_bytes: usize,
}

/// Why complete segment coverage did not close the requested distance width.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SegmentDistanceStop {
    /// The segment refinement allowance ended.
    SplitLimit,
    /// The next query cannot fit in the combined work allowance.
    WorkLimit,
    /// A segment cell has no strictly interior binary midpoint.
    ParameterResolution,
    /// A nested query or interval refinement refused; the parent stays covered.
    /// The original diagnostic is retained by `SegmentDistance::refusal`.
    RefinementRefused,
}

/// Work for this invocation, including unsuccessful nested query attempts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SegmentDistanceWork {
    /// Surface point-distance queries attempted.
    pub point_queries: u64,
    /// Committed binary segment refinements; both children are covered first.
    pub splits: u32,
    /// Segment cells remaining in the exhaustive frontier.
    pub remaining_cells: usize,
    /// Cells excluded by distance bounds or an empty-target proof.
    pub excluded_cells: usize,
    /// Aggregate work charge; a refused nested call retains its full allowance.
    pub charged_work_units: u128,
    /// Conservative peak live payload; a refused call retains its reservation.
    /// Allocator bookkeeping and bounded stack storage are excluded.
    pub retained_payload_bytes: usize,
}

/// An actual segment parameter paired with an actual surface parameter.
/// These copied values carry their binding through the enclosing result.
#[derive(Debug, Clone, Copy)]
pub struct SegmentDistanceWitness {
    /// Parameter of the real affine segment, in [0, 1].
    pub segment_parameter: f64,
    /// Outward evaluation of that affine point, including interpolation rounding.
    pub segment_position: [Interval; 3],
    /// Actual surface parameter, proved kept when trimming is requested.
    pub surface_parameter: [f64; 2],
    /// Outward evaluation of the actual surface point.
    pub surface_position: [Interval; 3],
}

/// A strict clearance decision for a ball swept along the entire segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SweptBallClearance {
    /// Every segment point is strictly farther than the supplied radius.
    Separated,
    /// An actual retained pair is strictly closer than the supplied radius.
    WithinRadius,
    /// The bounds do not decide the strict comparison; touching is included here.
    Undetermined,
    /// The trim kept region was proved empty.
    EmptyTarget,
}

/// Bounds on minimum unsigned clearance between a real affine segment and the
/// exact source surface, optionally restricted to the trim's kept-image closure.
#[derive(Debug, Clone)]
pub struct SegmentDistance<'s, 't> {
    source: AdmittedNurbsSurface<'s, f64>,
    trim: Option<AdmittedTrimmedPatch<'t>>,
    endpoints: [[f64; 3]; 2],
    tolerance: f64,
    bounds: Option<Interval>,
    witness: Option<SegmentDistanceWitness>,
    work: SegmentDistanceWork,
    refusal: Option<NurbsError>,
}

impl<'s, 't> SegmentDistance<'s, 't> {
    /// Original immutable geometric source, with stored homogeneous semantics.
    #[must_use]
    pub const fn source(&self) -> &'s NurbsSurface<f64> {
        self.source.source()
    }
    /// Original immutable trim source, if requested.
    #[must_use]
    pub fn trim(&self) -> Option<&'t TrimmedPatch> {
        self.trim.map(|trim| trim.source())
    }
    /// Exact stored binary endpoints defining P(t) = (1-t) A + t B.
    #[must_use]
    pub const fn endpoints(&self) -> [[f64; 3]; 2] {
        self.endpoints
    }
    /// Requested absolute clearance-interval width, in source Cartesian units.
    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }
    /// Outward unsigned-distance bounds; absent only for a proved empty target.
    #[must_use]
    pub const fn bounds(&self) -> Option<Interval> {
        self.bounds
    }
    /// Retained actual pair; absent when no exact kept surface point was found.
    #[must_use]
    pub const fn witness(&self) -> Option<SegmentDistanceWitness> {
        self.witness
    }
    /// Exhaustive coverage and combined resource accounting.
    #[must_use]
    pub const fn work(&self) -> SegmentDistanceWork {
        self.work
    }
    /// First refinement refusal, including inherited trim-classifier diagnostics.
    #[must_use]
    pub fn refusal(&self) -> Option<&NurbsError> {
        self.refusal.as_ref()
    }
    /// Compare a nonnegative finite ball radius with the rigorous clearance.
    /// This is unsigned separation or proximity, not signed solid overlap, an
    /// exact intersection, a contact time, or a first-ray-hit claim.
    ///
    /// # Errors
    /// Refuses a negative or nonfinite radius.
    pub fn swept_ball_clearance(&self, radius: f64) -> Result<SweptBallClearance, NurbsError> {
        if !radius.is_finite() || radius < 0.0 {
            return Err(refusal(
                "segment clearance needs a finite nonnegative radius",
            ));
        }
        Ok(match self.bounds {
            None => SweptBallClearance::EmptyTarget,
            Some(bounds) if bounds.lo() > radius => SweptBallClearance::Separated,
            Some(bounds) if self.witness.is_some() && bounds.hi() < radius => {
                SweptBallClearance::WithinRadius
            }
            Some(_) => SweptBallClearance::Undetermined,
        })
    }
}

/// Every published nonempty result covers the entire segment. Cancellation
/// publishes no bound, witness, partial frontier, or clearance decision.
#[derive(Debug, Clone)]
#[must_use]
pub enum SegmentDistanceRun<'s, 't> {
    /// A finite interval with an actual pair meets the requested width.
    Complete {
        /// Source-, trim-, and endpoints-bound result.
        distance: SegmentDistance<'s, 't>,
    },
    /// Valid complete coverage; the upper bound may be +infinity without a pair.
    Unresolved {
        /// Complete coverage, including the parent of any unfinished split.
        distance: SegmentDistance<'s, 't>,
        /// First binding refinement limit or refusal.
        stop: SegmentDistanceStop,
    },
    /// The surface's entire kept region was proved empty.
    Empty {
        /// Source binding and work without a distance or witness.
        distance: SegmentDistance<'s, 't>,
    },
    /// No partial result is published.
    Cancelled,
}

impl<'s> AdmittedNurbsSurface<'s, f64> {
    /// Enclose clearance over the entire real affine segment between finite
    /// binary endpoints, optionally against an exact nonzero-winding trim.
    ///
    /// Each segment cell uses an outward point-distance lower bound minus an
    /// outward radius enclosing the whole cell about the actual rounded query.
    /// Upper bounds are re-evaluated from an interval affine segment point and
    /// a retained actual surface witness; the rounded query is never the pair.
    /// A zero-length segment and zero refinement budgets are supported.
    ///
    /// Nested queries consume one combined allowance. A complete parent cell is
    /// replaced only after both child queries succeed. Point-query refinement
    /// limits can still yield useful segment bounds; a classifier refusal keeps
    /// its diagnostic and prevents a completion claim. No signed-solid, hit,
    /// normal, uniqueness, infinite-ray-tail, or executor-finalization claim.
    ///
    /// # Errors
    /// Refuses invalid input, excessive split limits, allocation failure, or
    /// resources/arithmetic insufficient to establish the initial full coverage.
    /// Later refusals preserve the parent and return `Unresolved`.
    ///
    /// # Panics
    /// Optional trim queries inherit the exact trim kernel's documented i128
    /// arithmetic envelope. A named Rat overflow publishes no result.
    pub fn segment_distance_enclosure_with_cx<'t>(
        &self,
        endpoints: [[f64; 3]; 2],
        trim: Option<AdmittedTrimmedPatch<'t>>,
        tolerance: f64,
        budget: SegmentDistanceBudget,
        cx: &Cx<'_>,
    ) -> Result<SegmentDistanceRun<'s, 't>, NurbsError> {
        segment_with_poll(*self, endpoints, trim, tolerance, budget, &mut || {
            cx.checkpoint().is_err()
        })
    }
}

#[derive(Clone, Copy)]
struct Cell {
    min: f64,
    max: f64,
    lower: f64,
    order: u64,
}
impl PartialEq for Cell {
    fn eq(&self, rhs: &Self) -> bool {
        self.order == rhs.order && self.lower.to_bits() == rhs.lower.to_bits()
    }
}
impl Eq for Cell {}
impl PartialOrd for Cell {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}
impl Ord for Cell {
    fn cmp(&self, rhs: &Self) -> Ordering {
        rhs.lower
            .total_cmp(&self.lower)
            .then_with(|| rhs.order.cmp(&self.order))
    }
}

fn affine(endpoints: [[f64; 3]; 2], t: Interval) -> Result<[Interval; 3], NurbsError> {
    let mut position = [Interval::point(0.0); 3];
    for k in 0..3 {
        let a = Interval::point(endpoints[0][k]);
        let b = Interval::point(endpoints[1][k]);
        position[k] = if endpoints[0][k] == endpoints[1][k] {
            a
        } else {
            finite(
                ((Interval::point(1.0) - t) * a + t * b)
                    .intersect(a.hull(b))
                    .ok_or_else(|| refusal("segment interval lost its affine convex hull"))?,
            )?
        };
    }
    Ok(position)
}

enum Evaluation {
    Covered,
    Empty,
    WorkLimit,
    Cancelled,
}

struct Search<'s, 't> {
    source: AdmittedNurbsSurface<'s, f64>,
    trim: Option<AdmittedTrimmedPatch<'t>>,
    endpoints: [[f64; 3]; 2],
    tolerance: f64,
    budget: SegmentDistanceBudget,
    meter: Meter,
    cell_work: u128,
    frontier_bytes: usize,
    retained: usize,
    point_queries: u64,
    upper: f64,
    witness: Option<SegmentDistanceWitness>,
    refusal: Option<NurbsError>,
}

impl Search<'_, '_> {
    fn evaluate(
        &mut self,
        cell: &mut Cell,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<Evaluation, NurbsError> {
        if cancel() {
            return Ok(Evaluation::Cancelled);
        }
        if !self.meter.charge(self.cell_work) {
            return Ok(Evaluation::WorkLimit);
        }
        let t = f64::midpoint(cell.min, cell.max);
        let point = affine(self.endpoints, Interval::point(t))?;
        let query = point.map(|x| f64::midpoint(x.lo(), x.hi()));
        let cell_box = affine(self.endpoints, Interval::new(cell.min, cell.max))?;
        let radius = distance_to_box(query, cell_box)?.hi();
        let work = self
            .budget
            .point_budget
            .max_work_units
            .min(MAX_WORK)
            .min(self.meter.limit - self.meter.used);
        if work == 0 {
            return Ok(Evaluation::WorkLimit);
        }
        // Prepay before entering the nested query: even an Err without a work
        // receipt stays bounded. Only successful receipts refund unused work.
        self.meter.used += work;
        let bytes = self
            .budget
            .point_budget
            .max_retained_bytes
            .min(self.budget.max_retained_bytes.min(MAX_BYTES) - self.frontier_bytes);
        self.point_queries += 1;
        let mut domain = self.trim.map(TrimDomain::new);
        let nested = search_with_poll(
            self.source,
            query,
            self.tolerance,
            SurfaceDistanceBudget {
                max_splits: self.budget.point_budget.max_splits,
                max_work_units: work,
                max_retained_bytes: bytes,
            },
            domain.as_mut(),
            cancel,
        );
        let result = match nested {
            Ok(SearchRun::Complete(result)) => result,
            Ok(SearchRun::Cancelled) => return Ok(Evaluation::Cancelled),
            Err(error) => {
                self.retained = self.retained.max(self.frontier_bytes + bytes);
                return Err(error);
            }
        };
        self.meter.used -= work - result.work.charged_work_units;
        self.retained = self
            .retained
            .max(self.frontier_bytes + result.work.retained_payload_bytes);
        if let Some(error) = domain.and_then(|domain| domain.refusal) {
            self.refusal.get_or_insert(error);
        }
        let Some(bounds) = result.bounds else {
            return Ok(Evaluation::Empty);
        };
        // Distance to any fixed nonempty set (or its closure) is 1-Lipschitz.
        let lower = (Interval::point(bounds.lo()) - Interval::point(radius))
            .lo()
            .max(0.0);
        cell.lower = cell.lower.max(lower);
        if let Some(witness) = result.witness {
            let difference = std::array::from_fn(|k| witness.position[k] - point[k]);
            let upper = distance_to_box([0.0; 3], difference)?.hi();
            if upper < self.upper {
                self.upper = upper;
                self.witness = Some(SegmentDistanceWitness {
                    segment_parameter: t,
                    segment_position: point,
                    surface_parameter: witness.param,
                    surface_position: witness.position,
                });
            }
        }
        Ok(Evaluation::Covered)
    }
}

fn segment_with_poll<'s, 't>(
    source: AdmittedNurbsSurface<'s, f64>,
    endpoints: [[f64; 3]; 2],
    trim: Option<AdmittedTrimmedPatch<'t>>,
    tolerance: f64,
    budget: SegmentDistanceBudget,
    cancel: &mut impl FnMut() -> bool,
) -> Result<SegmentDistanceRun<'s, 't>, NurbsError> {
    if endpoints.iter().flatten().any(|x| !x.is_finite())
        || !tolerance.is_finite()
        || tolerance <= 0.0
        || budget.max_splits > CLOSEST_MAX_SPLITS
        || budget.point_budget.max_splits > CLOSEST_MAX_SPLITS
    {
        return Err(refusal(
            "segment-distance needs finite endpoints, positive finite tolerance, and bounded splits",
        ));
    }
    if cancel() {
        return Ok(SegmentDistanceRun::Cancelled);
    }
    let capacity = checked((budget.max_splits as usize).checked_add(1))?;
    let bytes_limit = budget.max_retained_bytes.min(MAX_BYTES);
    if checked(capacity.checked_mul(size_of::<Cell>()))? > bytes_limit {
        return Err(refusal(
            "segment-distance frontier exceeds combined retained payload",
        ));
    }
    let mut frontier = BinaryHeap::new();
    frontier
        .try_reserve_exact(capacity)
        .map_err(|_| refusal("segment-distance frontier allocation refused"))?;
    let frontier_bytes = checked(frontier.capacity().checked_mul(size_of::<Cell>()))?;
    if frontier_bytes > bytes_limit {
        return Err(refusal(
            "segment-distance allocator capacity exceeds retained payload",
        ));
    }
    let mut search = Search {
        source,
        trim,
        endpoints,
        tolerance,
        budget,
        meter: Meter {
            used: 0,
            limit: budget.max_work_units.min(MAX_WORK),
        },
        cell_work: 1024 + 16 * (usize::BITS - capacity.leading_zeros()) as u128,
        frontier_bytes,
        retained: frontier_bytes,
        point_queries: 0,
        upper: f64::INFINITY,
        witness: None,
        refusal: None,
    };
    let mut seed = Cell {
        min: 0.0,
        max: 1.0,
        lower: 0.0,
        order: 0,
    };
    let mut empty = match search.evaluate(&mut seed, cancel)? {
        Evaluation::Covered => {
            frontier.push(seed);
            false
        }
        Evaluation::Empty => true,
        Evaluation::Cancelled => return Ok(SegmentDistanceRun::Cancelled),
        Evaluation::WorkLimit => {
            return Err(refusal("segment-distance work cannot cover the segment"));
        }
    };
    let mut splits = 0;
    let mut excluded = usize::from(empty);
    let mut order = 1;
    let stop = 'refine: loop {
        if cancel() {
            return Ok(SegmentDistanceRun::Cancelled);
        }
        if search.refusal.is_some() {
            break Some(SegmentDistanceStop::RefinementRefused);
        }
        let Some(cell) = frontier.peek().copied() else {
            break None;
        };
        if cell.lower >= search.upper {
            excluded += frontier.len();
            frontier.clear();
            break None;
        }
        if search.witness.is_some() && Interval::new(cell.lower, search.upper).width() <= tolerance
        {
            break None;
        }
        if splits == budget.max_splits {
            break Some(SegmentDistanceStop::SplitLimit);
        }
        let mid = f64::midpoint(cell.min, cell.max);
        if mid <= cell.min || mid >= cell.max {
            break Some(SegmentDistanceStop::ParameterResolution);
        }
        let mut children = [cell; 2];
        children[0].max = mid;
        children[1].min = mid;
        // Keep the complete parent live until BOTH child calls have succeeded.
        for child in &mut children {
            child.order = order;
            order += 1;
            match search.evaluate(child, cancel) {
                Ok(Evaluation::Covered) => {}
                Ok(Evaluation::Empty) => {
                    empty = true;
                    search.witness = None;
                    excluded += frontier.len();
                    frontier.clear();
                    break 'refine None;
                }
                Ok(Evaluation::Cancelled) => return Ok(SegmentDistanceRun::Cancelled),
                Ok(Evaluation::WorkLimit) => break 'refine Some(SegmentDistanceStop::WorkLimit),
                Err(error) => {
                    search.refusal.get_or_insert(error);
                    break 'refine Some(SegmentDistanceStop::RefinementRefused);
                }
            }
            if search.refusal.is_some() {
                break 'refine Some(SegmentDistanceStop::RefinementRefused);
            }
        }
        frontier.pop();
        for child in children {
            if child.lower >= search.upper {
                excluded += 1;
            } else {
                frontier.push(child);
            }
        }
        splits += 1;
    };
    let bounds = (!empty).then(|| {
        let lower = frontier
            .peek()
            .map_or(search.upper, |cell| cell.lower)
            .min(search.upper);
        Interval::new(lower, search.upper)
    });
    let distance = SegmentDistance {
        source,
        trim,
        endpoints,
        tolerance,
        bounds,
        witness: search.witness,
        work: SegmentDistanceWork {
            point_queries: search.point_queries,
            splits,
            remaining_cells: frontier.len(),
            excluded_cells: excluded,
            charged_work_units: search.meter.used,
            retained_payload_bytes: search.retained,
        },
        refusal: search.refusal,
    };
    drop(frontier);
    if cancel() {
        return Ok(SegmentDistanceRun::Cancelled);
    }
    Ok(match stop {
        Some(stop) => SegmentDistanceRun::Unresolved { distance, stop },
        None if empty => SegmentDistanceRun::Empty { distance },
        None => SegmentDistanceRun::Complete { distance },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::KnotVector;

    fn fixture() -> (NurbsSurface<f64>, SegmentDistanceBudget, [[f64; 3]; 2]) {
        let k = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
        let surface = NurbsSurface::new(
            k.clone(),
            k,
            &[
                vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
            ],
            &[vec![1.0; 2], vec![1.0; 2]],
        )
        .unwrap();
        let budget = SegmentDistanceBudget {
            max_splits: 4,
            point_budget: SurfaceDistanceBudget {
                max_splits: 0,
                max_work_units: 1_000_000,
                max_retained_bytes: 1_000_000,
            },
            max_work_units: 16_000_000,
            max_retained_bytes: 2_000_000,
        };
        (surface, budget, [[0.25, 0.5, 1.0], [0.75, 0.5, 1.0]])
    }

    #[test]
    fn g4_nested_refusal_keeps_parent_coverage_and_charges_attempted_work() {
        let (surface, mut budget, endpoints) = fixture();
        let source = surface.admit().unwrap();
        let SegmentDistanceRun::Complete { distance: seed } =
            segment_with_poll(source, endpoints, None, 10.0, budget, &mut || false).unwrap()
        else {
            panic!("loose seed width");
        };
        let cell_work = 1024 + 16 * (usize::BITS - 5usize.leading_zeros()) as u128;
        budget.max_work_units = seed.work().charged_work_units + cell_work + 1;
        let SegmentDistanceRun::Unresolved {
            distance,
            stop: SegmentDistanceStop::RefinementRefused,
        } = segment_with_poll(source, endpoints, None, 1e-9, budget, &mut || false).unwrap()
        else {
            panic!("nested source admission should refuse its one work unit");
        };
        assert!(distance.refusal().is_some());
        assert_eq!(distance.work().point_queries, 2);
        assert_eq!(distance.work().splits, 0);
        assert_eq!(distance.work().remaining_cells, 1);
        assert_eq!(distance.work().charged_work_units, budget.max_work_units);
        assert!(distance.bounds().unwrap().contains(1.0));
    }

    #[test]
    fn g4_nested_and_final_cancellation_publish_no_segment_result() {
        let (surface, budget, endpoints) = fixture();
        let source = surface.admit().unwrap();
        let mut polls = 0;
        let baseline = segment_with_poll(source, endpoints, None, 1e-9, budget, &mut || {
            polls += 1;
            false
        })
        .unwrap();
        assert!(!matches!(baseline, SegmentDistanceRun::Cancelled));
        assert!(polls > 20);
        for boundary in [1, 2, polls / 3, polls / 2, polls - 1, polls] {
            let mut visited = 0;
            let run = segment_with_poll(source, endpoints, None, 1e-9, budget, &mut || {
                visited += 1;
                visited == boundary
            })
            .unwrap();
            assert!(
                matches!(run, SegmentDistanceRun::Cancelled),
                "checkpoint {boundary}"
            );
            assert_eq!(visited, boundary);
        }
    }
}
