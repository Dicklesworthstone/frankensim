//! Outward unsigned-distance bounds for one admitted, untrimmed rational
//! B-spline surface. This lane does not use the measured Bezier conversion.
//!
//! Each original nonzero knot-span rectangle is covered. Interval tensor de
//! Boor encloses its homogeneous image, positive weights permit projection,
//! and a Cartesian control hull independently bounds that projection. A
//! deterministic frontier refines the smallest lower bound. Only an interval
//! evaluation at an actual parameter supplies a witness upper bound.
//!
//! The mathematical object is the real spline defined by the STORED binary
//! homogeneous controls and knots, in the same Cartesian units as the query.
//! Construction roundoff relative to earlier Cartesian input is not undone.
//! Discontinuous knot spans contribute their image infimum, including limits;
//! the retained witness is always evaluated on the actual owned side.
//!
//! This is a source-bound arithmetic primitive, not a portable authority
//! artifact, signed field, trim certificate, or executor lifecycle owner.

use crate::NurbsError;
use crate::basis::AdmittedKnotVector;
use crate::closest::CLOSEST_MAX_SPLITS;
use crate::surface::{AdmittedNurbsSurface, NurbsSurface};
use fs_exec::Cx;
use fs_ivl::Interval;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::mem::size_of;

const MAX_WORK: u128 = 1_073_741_824;
const MAX_BYTES: usize = 256 * 1024 * 1024;

type Homogeneous = [Interval; 4];

/// Caller limits for the complete admitted query, including original-span
/// coverage. Work units conservatively price interval blends and source visits;
/// they are not machine instructions or elapsed time. Retained payload includes
/// the borrowed source, frontier, span indices and reusable scratch. Allocator
/// bookkeeping and bounded stack storage are outside that payload count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SurfaceDistanceBudget {
    /// Number of binary parameter-cell refinements after complete seeding.
    pub max_splits: u32,
    /// One allowance shared by planning, hull construction and all evaluations.
    pub max_work_units: u128,
    /// Simultaneously retained payload limit; a 256 MiB ceiling also applies.
    pub max_retained_bytes: usize,
}

/// Why a valid distance enclosure did not meet the requested width.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceDistanceStop {
    /// No further binary refinements were admitted.
    SplitLimit,
    /// Refining the next complete cell would exceed the shared work allowance.
    WorkLimit,
    /// The next complete cell expansion would exceed the retained payload cap.
    MemoryLimit,
    /// The most relevant cell has no strictly interior representable midpoint.
    ParameterResolution,
}

/// Measured work and retained coverage information for this invocation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SurfaceDistanceWork {
    /// Every nonzero original span rectangle was seeded before publication.
    pub original_span_cells: usize,
    /// Actual binary refinements.
    pub splits: u32,
    /// Cells still represented by the final frontier.
    pub remaining_cells: usize,
    /// Cells removed only after an outward lower bound excluded improvement.
    pub excluded_cells: usize,
    /// Charged conservative interval/source/heap work units.
    pub charged_work_units: u128,
    /// Admitted simultaneously live payload, including the borrowed source.
    pub retained_payload_bytes: usize,
}

/// Unforgeable through this API: the interval, witness and original source are
/// published together after complete span coverage and a final checkpoint.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceDistanceEnclosure<'a> {
    source: AdmittedNurbsSurface<'a, f64>,
    query: [f64; 3],
    tolerance: f64,
    bounds: Interval,
    witness_param: [f64; 2],
    witness_position: [Interval; 3],
    work: SurfaceDistanceWork,
}

impl<'a> SurfaceDistanceEnclosure<'a> {
    /// Exact immutable generation to which every bound applies.
    #[must_use]
    pub const fn source(&self) -> &'a NurbsSurface<f64> {
        self.source.source()
    }
    /// Original finite query in the source's Cartesian units.
    #[must_use]
    pub const fn query(&self) -> [f64; 3] {
        self.query
    }
    /// Requested absolute Euclidean-distance interval width, in source units.
    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }
    /// Outward bounds on the unsigned distance infimum over the whole domain.
    #[must_use]
    pub const fn bounds(&self) -> Interval {
        self.bounds
    }
    /// Actual in-domain parameter whose evaluated point supplies the upper bound.
    #[must_use]
    pub const fn witness_param(&self) -> [f64; 2] {
        self.witness_param
    }
    /// Outward enclosure of that witness's exact projected Cartesian position.
    #[must_use]
    pub const fn witness_position(&self) -> [Interval; 3] {
        self.witness_position
    }
    /// Coverage and budget accounting for the completed invocation.
    #[must_use]
    pub const fn work(&self) -> SurfaceDistanceWork {
        self.work
    }
}

/// Cancellation publishes no enclosure. Both other variants have complete
/// original-span coverage; only `Complete` meets the requested outward width.
#[derive(Debug, Clone, Copy)]
#[must_use]
pub enum SurfaceDistanceRun<'a> {
    /// The outward distance interval meets the requested absolute width.
    Complete {
        /// Bound to the admitted source and original query.
        enclosure: SurfaceDistanceEnclosure<'a>,
    },
    /// A valid, possibly useful enclosure without requested-width closure.
    Unresolved {
        /// Complete coverage with a wider interval than requested.
        enclosure: SurfaceDistanceEnclosure<'a>,
        /// First binding refinement limit.
        stop: SurfaceDistanceStop,
    },
    /// No partial bound, point, or resumable frontier is published.
    Cancelled,
}

impl<'a> AdmittedNurbsSurface<'a, f64> {
    /// Enclose unsigned distance to this untrimmed surface under one shared
    /// work/memory allowance and cancellation gate. Structural admission is
    /// already represented by `self`; this method does not repeat it.
    ///
    /// Norms and coordinates use the source's Cartesian units. Zero splits
    /// still cover all original spans and can return a useful enclosure.
    /// The caller owns executor scope drain/finalization. This primitive does
    /// not consume `Cx`'s budget or establish sign, topology, trim membership,
    /// a unique minimizer, or geometric accuracy of the authored input.
    ///
    /// # Errors
    /// Refuses nonfinite input, nonpositive tolerance, excessive split limits,
    /// insufficient resources to cover ALL original spans, allocation failure,
    /// or arithmetic that cannot retain finite outward bounds.
    pub fn distance_enclosure_with_cx(
        &self,
        query: [f64; 3],
        tolerance: f64,
        budget: SurfaceDistanceBudget,
        cx: &Cx<'_>,
    ) -> Result<SurfaceDistanceRun<'a>, NurbsError> {
        distance_with_poll(*self, query, tolerance, budget, &mut || {
            cx.checkpoint().is_err()
        })
    }
}

#[derive(Debug, Clone, Copy)]
struct Cell {
    span: [usize; 2],
    min: [f64; 2],
    max: [f64; 2],
    lower: f64,
    order: u64,
    depth: [u32; 2],
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

struct Meter {
    used: u128,
    limit: u128,
}
impl Meter {
    fn charge(&mut self, amount: u128) -> bool {
        if amount > self.limit - self.used {
            return false;
        }
        self.used += amount;
        true
    }
}

fn refusal(what: &str) -> NurbsError {
    NurbsError::Domain {
        what: what.to_string(),
    }
}
fn checked(value: Option<usize>) -> Result<usize, NurbsError> {
    value.ok_or_else(|| refusal("surface-distance size accounting overflow"))
}
fn alloc<T>(count: usize) -> Result<Vec<T>, NurbsError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| refusal("surface-distance allocation refused"))?;
    Ok(result)
}
fn tick(steps: &mut usize, cancel: &mut impl FnMut() -> bool) -> bool {
    *steps += 1;
    if *steps == 64 {
        *steps = 0;
        cancel()
    } else {
        false
    }
}
fn finite(value: Interval) -> Result<Interval, NurbsError> {
    if value.lo().is_finite() && value.hi().is_finite() {
        Ok(value)
    } else {
        Err(refusal(
            "surface-distance outward arithmetic exceeds the finite range",
        ))
    }
}

fn source_payload(
    source: AdmittedNurbsSurface<'_, f64>,
    meter: &mut Meter,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Option<usize>, NurbsError> {
    let raw = source.source();
    let mut bytes = checked(
        raw.cpw
            .capacity()
            .checked_mul(size_of::<Vec<[f64; 4]>>())
            .and_then(|n| n.checked_add(size_of::<NurbsSurface<f64>>()))
            .and_then(|n| {
                n.checked_add(raw.knots_u.knots.capacity().checked_mul(size_of::<f64>())?)
            })
            .and_then(|n| {
                n.checked_add(raw.knots_v.knots.capacity().checked_mul(size_of::<f64>())?)
            }),
    )?;
    let mut steps = 0;
    for row in &raw.cpw {
        if !meter.charge(1) {
            return Err(refusal("surface-distance work cannot admit source storage"));
        }
        if tick(&mut steps, cancel) {
            return Ok(None);
        }
        bytes = checked(
            row.capacity()
                .checked_mul(size_of::<[f64; 4]>())
                .and_then(|n| bytes.checked_add(n)),
        )?;
    }
    if cancel() { Ok(None) } else { Ok(Some(bytes)) }
}

fn span_count(
    knots: AdmittedKnotVector<'_, f64>,
    meter: &mut Meter,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Option<usize>, NurbsError> {
    let mut count = 0;
    let mut steps = 0;
    for i in knots.degree()..knots.control_count() {
        if !meter.charge(1) {
            return Err(refusal(
                "surface-distance work cannot cover original knot spans",
            ));
        }
        if tick(&mut steps, cancel) {
            return Ok(None);
        }
        count += usize::from(knots.knots()[i] < knots.knots()[i + 1]);
    }
    if cancel() { Ok(None) } else { Ok(Some(count)) }
}

fn span_indices(
    knots: AdmittedKnotVector<'_, f64>,
    out: &mut Vec<usize>,
    cancel: &mut impl FnMut() -> bool,
) -> bool {
    let mut steps = 0;
    for i in knots.degree()..knots.control_count() {
        if tick(&mut steps, cancel) {
            return false;
        }
        if knots.knots()[i] < knots.knots()[i + 1] {
            out.push(i);
        }
    }
    !cancel()
}

fn control_hull(
    source: AdmittedNurbsSurface<'_, f64>,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Option<[Interval; 3]>, NurbsError> {
    let mut hull: Option<[Interval; 3]> = None;
    let mut steps = 0;
    for row in source.homogeneous_control_net() {
        for h in row {
            if tick(&mut steps, cancel) {
                return Ok(None);
            }
            let mut p = [Interval::point(0.0); 3];
            for k in 0..3 {
                p[k] = finite(Interval::point(h[k]) / Interval::point(h[3]))?;
            }
            hull = Some(match hull {
                None => p,
                Some(old) => std::array::from_fn(|k| old[k].hull(p[k])),
            });
        }
    }
    if cancel() { Ok(None) } else { Ok(hull) }
}

/// Stable outward Euclidean norm, including subnormal magnitudes. Scaling is
/// itself interval arithmetic; the selected finite scale is an exact binary
/// constant, not an assumed exact norm or an unrounded difference.
fn distance_to_box(query: [f64; 3], position: [Interval; 3]) -> Result<Interval, NurbsError> {
    let mut delta = [Interval::point(0.0); 3];
    let mut scale = 0.0f64;
    for k in 0..3 {
        delta[k] = finite(position[k] - Interval::point(query[k]))?.abs();
        scale = scale.max(delta[k].hi());
    }
    if scale == 0.0 {
        return Ok(Interval::point(0.0));
    }
    let mut square = Interval::point(0.0);
    for d in delta {
        let unit = d / Interval::point(scale);
        square = square + unit * unit;
    }
    let result = finite(square.sqrt() * Interval::point(scale))?;
    Ok(Interval::new(result.lo().max(0.0), result.hi()))
}

fn de_boor(
    values: &mut [Homogeneous],
    knots: &[f64],
    span: usize,
    parameter: Interval,
    steps: &mut usize,
    cancel: &mut impl FnMut() -> bool,
) -> Result<bool, NurbsError> {
    let degree = values.len() - 1;
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            if tick(steps, cancel) {
                return Ok(false);
            }
            let index = span - degree + j;
            // Every parameter in the current original span lies between these
            // two knots. Intersecting alpha with [0,1] remains sound even if
            // outward subtraction of extremely close knots contains zero.
            let alpha = ((parameter - Interval::point(knots[index]))
                / (Interval::point(knots[index + degree + 1 - r]) - Interval::point(knots[index])))
            .intersect(Interval::new(0.0, 1.0))
            .ok_or_else(|| refusal("surface-distance de Boor parameter escaped its span"))?;
            let left = values[j - 1];
            let right = values[j];
            for k in 0..4 {
                values[j][k] = ((Interval::point(1.0) - alpha) * left[k] + alpha * right[k])
                    .intersect(left[k].hull(right[k]))
                    .ok_or_else(|| {
                        refusal("surface-distance interval blend lost its convex hull")
                    })?;
            }
        }
    }
    Ok(true)
}

struct Scratch {
    row: Vec<Homogeneous>,
    column: Vec<Homogeneous>,
}
fn enclose_cell(
    source: AdmittedNurbsSurface<'_, f64>,
    span: [usize; 2],
    min: [f64; 2],
    max: [f64; 2],
    hull: [Interval; 3],
    scratch: &mut Scratch,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Option<[Interval; 3]>, NurbsError> {
    let ku = source.knots_u();
    let kv = source.knots_v();
    let mut steps = 0;
    for j in 0..scratch.column.len() {
        for i in 0..scratch.row.len() {
            if tick(&mut steps, cancel) {
                return Ok(None);
            }
            let control = source.homogeneous_control_net()[span[0] - ku.degree() + i]
                [span[1] - kv.degree() + j];
            scratch.row[i] = control.map(Interval::point);
        }
        if !de_boor(
            &mut scratch.row,
            ku.knots(),
            span[0],
            Interval::new(min[0], max[0]),
            &mut steps,
            cancel,
        )? {
            return Ok(None);
        }
        scratch.column[j] = scratch.row[ku.degree()];
    }
    if !de_boor(
        &mut scratch.column,
        kv.knots(),
        span[1],
        Interval::new(min[1], max[1]),
        &mut steps,
        cancel,
    )? {
        return Ok(None);
    }
    let h = scratch.column[kv.degree()];
    if h[3].lo() <= 0.0 {
        return Err(refusal(
            "surface-distance projection lacks a strictly positive weight enclosure",
        ));
    }
    let mut point = [Interval::point(0.0); 3];
    for k in 0..3 {
        point[k] = finite((h[k] / h[3]).intersect(hull[k]).ok_or_else(|| {
            refusal("surface-distance projection lost the Cartesian control hull")
        })?)?;
    }
    if cancel() { Ok(None) } else { Ok(Some(point)) }
}

#[derive(Clone, Copy)]
struct Witness {
    upper: f64,
    param: [f64; 2],
    position: [Interval; 3],
}
fn evaluate_cell(
    source: AdmittedNurbsSurface<'_, f64>,
    query: [f64; 3],
    cell: &mut Cell,
    hull: [Interval; 3],
    scratch: &mut Scratch,
    witness: &mut Option<Witness>,
    cancel: &mut impl FnMut() -> bool,
) -> Result<bool, NurbsError> {
    let Some(position) =
        enclose_cell(source, cell.span, cell.min, cell.max, hull, scratch, cancel)?
    else {
        return Ok(false);
    };
    cell.lower = cell.lower.max(distance_to_box(query, position)?.lo());
    let mut param = [
        f64::midpoint(cell.min[0], cell.max[0]),
        f64::midpoint(cell.min[1], cell.max[1]),
    ];
    // At full-multiplicity knots the actual evaluator owns the right span.
    // A rounded midpoint must not become a left-limit witness for that knot.
    for (axis, knots) in [source.knots_u(), source.knots_v()].into_iter().enumerate() {
        if param[axis] == knots.knots()[cell.span[axis] + 1] {
            param[axis] = cell.min[axis];
        }
    }
    let Some(position) = enclose_cell(source, cell.span, param, param, hull, scratch, cancel)?
    else {
        return Ok(false);
    };
    let upper = distance_to_box(query, position)?.hi();
    if witness.is_none_or(|old| upper < old.upper) {
        *witness = Some(Witness {
            upper,
            param,
            position,
        });
    }
    Ok(!cancel())
}

fn distance_with_poll<'a>(
    source: AdmittedNurbsSurface<'a, f64>,
    query: [f64; 3],
    tolerance: f64,
    budget: SurfaceDistanceBudget,
    cancel: &mut impl FnMut() -> bool,
) -> Result<SurfaceDistanceRun<'a>, NurbsError> {
    if query.iter().any(|x| !x.is_finite()) || !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(refusal(
            "surface-distance needs a finite point and positive finite absolute tolerance",
        ));
    }
    if budget.max_splits > CLOSEST_MAX_SPLITS {
        return Err(refusal(
            "surface-distance split request exceeds the defensive ceiling",
        ));
    }
    if cancel() {
        return Ok(SurfaceDistanceRun::Cancelled);
    }
    let ku = source.knots_u();
    let kv = source.knots_v();
    let (pu, pv, nu, nv) = (
        ku.degree(),
        kv.degree(),
        ku.control_count(),
        kv.control_count(),
    );
    let mut meter = Meter {
        used: 0,
        limit: budget.max_work_units.min(MAX_WORK),
    };
    let Some(source_bytes) = source_payload(source, &mut meter, cancel)? else {
        return Ok(SurfaceDistanceRun::Cancelled);
    };
    let Some(su) = span_count(ku, &mut meter, cancel)? else {
        return Ok(SurfaceDistanceRun::Cancelled);
    };
    let Some(sv) = span_count(kv, &mut meter, cancel)? else {
        return Ok(SurfaceDistanceRun::Cancelled);
    };
    let seeds = checked(su.checked_mul(sv))?;
    if seeds == 0 {
        return Err(refusal(
            "surface-distance has no nonzero original knot span",
        ));
    }
    let scratch_entries = checked(pu.checked_add(pv).and_then(|n| n.checked_add(2)))?;
    let scratch_bytes = checked(scratch_entries.checked_mul(size_of::<Homogeneous>()))?;
    let span_bytes = checked(
        su.checked_add(sv)
            .and_then(|n| n.checked_mul(size_of::<usize>())),
    )?;
    let fixed_bytes = checked(
        source_bytes
            .checked_add(scratch_bytes)
            .and_then(|n| n.checked_add(span_bytes)),
    )?;
    let bytes_limit = budget.max_retained_bytes.min(MAX_BYTES);
    let available = bytes_limit.checked_sub(fixed_bytes).ok_or_else(|| {
        refusal("surface-distance retained payload cannot hold source and scratch")
    })?;
    let desired = checked(seeds.checked_add(budget.max_splits as usize))?;
    let capacity = desired.min(available / size_of::<Cell>());
    if capacity < seeds {
        return Err(refusal(
            "surface-distance retained payload cannot cover all original spans",
        ));
    }
    let blend_count =
        (pu as u128 * (pu as u128 + 1) / 2) * (pv as u128 + 1) + pv as u128 * (pv as u128 + 1) / 2;
    let evaluation_work = 64 * blend_count + 4 * (pu as u128 + 1) * (pv as u128 + 1) + 256;
    let heap_work = 8 * (usize::BITS - capacity.leading_zeros()) as u128 + 8;
    let cell_work = 2 * evaluation_work + heap_work;
    let seed_work = seeds as u128 * cell_work + nu as u128 * nv as u128 * 64 + (nu + nv) as u128;
    if !meter.charge(seed_work) {
        return Err(refusal(
            "surface-distance work cannot establish complete original-span coverage",
        ));
    }
    let mut spans_u = alloc(su)?;
    let mut spans_v = alloc(sv)?;
    let mut scratch = Scratch {
        row: alloc(pu + 1)?,
        column: alloc(pv + 1)?,
    };
    scratch.row.resize(pu + 1, [Interval::point(0.0); 4]);
    scratch.column.resize(pv + 1, [Interval::point(0.0); 4]);
    let mut frontier = BinaryHeap::new();
    frontier
        .try_reserve_exact(capacity)
        .map_err(|_| refusal("surface-distance frontier allocation refused"))?;
    let actual_spans = checked(
        spans_u
            .capacity()
            .checked_add(spans_v.capacity())
            .and_then(|n| n.checked_mul(size_of::<usize>())),
    )?;
    let actual_scratch = checked(
        scratch
            .row
            .capacity()
            .checked_add(scratch.column.capacity())
            .and_then(|n| n.checked_mul(size_of::<Homogeneous>())),
    )?;
    let retained = checked(
        source_bytes
            .checked_add(actual_spans)
            .and_then(|n| n.checked_add(actual_scratch))
            .and_then(|n| n.checked_add(frontier.capacity().checked_mul(size_of::<Cell>())?)),
    )?;
    if retained > bytes_limit {
        return Err(refusal(
            "surface-distance allocator capacity exceeds retained payload admission",
        ));
    }
    if !span_indices(ku, &mut spans_u, cancel) || !span_indices(kv, &mut spans_v, cancel) {
        return Ok(SurfaceDistanceRun::Cancelled);
    }
    let Some(hull) = control_hull(source, cancel)? else {
        return Ok(SurfaceDistanceRun::Cancelled);
    };
    let mut witness = None;
    let mut order = 0u64;
    for &u in &spans_u {
        for &v in &spans_v {
            let mut cell = Cell {
                span: [u, v],
                min: [ku.knots()[u], kv.knots()[v]],
                max: [ku.knots()[u + 1], kv.knots()[v + 1]],
                lower: 0.0,
                order,
                depth: [0, 0],
            };
            order += 1;
            if !evaluate_cell(
                source,
                query,
                &mut cell,
                hull,
                &mut scratch,
                &mut witness,
                cancel,
            )? {
                return Ok(SurfaceDistanceRun::Cancelled);
            }
            frontier.push(cell);
        }
    }
    let mut splits = 0u32;
    let mut excluded = 0usize;
    let stop = loop {
        if cancel() {
            return Ok(SurfaceDistanceRun::Cancelled);
        }
        let upper = witness
            .expect("all original span cells were evaluated")
            .upper;
        let Some(cell) = frontier.peek().copied() else {
            break None;
        };
        if cell.lower >= upper {
            excluded += frontier.len();
            frontier.clear();
            break None;
        }
        if Interval::new(cell.lower, upper).width() <= tolerance {
            break None;
        }
        if splits == budget.max_splits {
            break Some(SurfaceDistanceStop::SplitLimit);
        }
        if frontier.len() == capacity {
            break Some(SurfaceDistanceStop::MemoryLimit);
        }
        let mids = [
            f64::midpoint(cell.min[0], cell.max[0]),
            f64::midpoint(cell.min[1], cell.max[1]),
        ];
        let possible = [
            mids[0] > cell.min[0] && mids[0] < cell.max[0],
            mids[1] > cell.min[1] && mids[1] < cell.max[1],
        ];
        if !possible[0] && !possible[1] {
            break Some(SurfaceDistanceStop::ParameterResolution);
        }
        // Balance halvings in the original span coordinates. Logical depths
        // avoid overflow/underflow and dependence on unrelated U/V knot units.
        let axis = if possible[0] && (!possible[1] || cell.depth[0] <= cell.depth[1]) {
            0
        } else {
            1
        };
        if !meter.charge(2 * cell_work + heap_work) {
            break Some(SurfaceDistanceStop::WorkLimit);
        }
        frontier.pop();
        for half in 0..2 {
            let mut child = cell;
            child.order = order;
            order += 1;
            child.depth[axis] += 1;
            if half == 0 {
                child.max[axis] = mids[axis];
            } else {
                child.min[axis] = mids[axis];
            }
            if !evaluate_cell(
                source,
                query,
                &mut child,
                hull,
                &mut scratch,
                &mut witness,
                cancel,
            )? {
                return Ok(SurfaceDistanceRun::Cancelled);
            }
            if child.lower >= witness.expect("retained witness").upper {
                excluded += 1;
            } else {
                frontier.push(child);
            }
        }
        splits += 1;
    };
    let witness = witness.expect("complete initial coverage supplies a witness");
    let lower = frontier
        .peek()
        .map_or(witness.upper, |c| c.lower)
        .min(witness.upper);
    let bounds = Interval::new(lower, witness.upper);
    let work = SurfaceDistanceWork {
        original_span_cells: seeds,
        splits,
        remaining_cells: frontier.len(),
        excluded_cells: excluded,
        charged_work_units: meter.used,
        retained_payload_bytes: retained,
    };
    drop(frontier);
    drop(scratch);
    drop(spans_u);
    drop(spans_v);
    if cancel() {
        return Ok(SurfaceDistanceRun::Cancelled);
    }
    let enclosure = SurfaceDistanceEnclosure {
        source,
        query,
        tolerance,
        bounds,
        witness_param: witness.param,
        witness_position: witness.position,
        work,
    };
    Ok(match stop {
        None => SurfaceDistanceRun::Complete { enclosure },
        Some(stop) => SurfaceDistanceRun::Unresolved { enclosure, stop },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::KnotVector;

    #[test]
    fn g4_cancel_during_enclosure_refinement_or_final_publication_returns_no_bound() {
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
        let source = surface.admit().unwrap();
        let budget = SurfaceDistanceBudget {
            max_splits: 32,
            max_work_units: 1_000_000,
            max_retained_bytes: 1_000_000,
        };
        let mut polls = 0;
        let run = distance_with_poll(source, [0.13, 0.22, 1.0], 1e-12, budget, &mut || {
            polls += 1;
            false
        })
        .unwrap();
        assert!(!matches!(run, SurfaceDistanceRun::Cancelled));
        assert!(polls > 20);
        for boundary in [1, 2, polls / 4, polls / 2, polls - 1, polls] {
            let mut visited = 0;
            let run = distance_with_poll(source, [0.13, 0.22, 1.0], 1e-12, budget, &mut || {
                visited += 1;
                visited == boundary
            })
            .unwrap();
            assert!(
                matches!(run, SurfaceDistanceRun::Cancelled),
                "checkpoint {boundary}"
            );
            assert_eq!(visited, boundary, "no work after observed cancellation");
        }
    }
}
