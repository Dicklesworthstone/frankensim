//! Global, bounded queries on six-node Lagrange triangle patches.
//!
//! The geometry is the quadratic patch, including its curved edges, rather
//! than a tessellation. Directed Bernstein hulls and supporting planes bound
//! distance; subdivision covers the complete closed parameter triangle.
//! Local optimization supplies feasible upper witnesses only. Directed
//! first-order necessary conditions can exclude cells containing no global
//! minimizer, including during positional localization of affine patches.
//!
//! A successful receipt encloses the global minimum distance AND every
//! possible minimizing point near the returned witness. It does not prove
//! uniqueness. Well-separated equal minima cannot meet a smaller authored
//! witness tolerance and cause an explicit bounded-work refusal.

use super::{TriangleFeature, TriangleQueryError, closest_triangle_point, cross, dot, norm, poll};
use fs_exec::Cx;
use fs_ivl::Interval;
use std::collections::BinaryHeap;

mod surface;
pub use surface::{
    PreparedQuadraticTriangleSurface, QuadraticTrianglePointQuery, RefitQuadraticTriangleSurface,
    SweptQuadraticTriangleSurface,
};

/// Conservative admission of declared linear sphere/quadratic-patch paths.
pub mod sweep;

const MAXIMUM_CELLS: usize = 4_000_000;
const MAXIMUM_ACTIVE_CELLS: usize = 262_144;
const MAXIMUM_DEPTH: u8 = 48;
const CORNERS: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const EDGE_PAIRS: [(usize, usize); 3] = [(0, 1), (1, 2), (2, 0)];

/// Absolute distance and position tolerances, in the caller's coordinate unit.
/// A cell has a fixed bounded number of candidate/derivative evaluations.
/// `maximum_cells` counts every constructed cell, including discarded ones;
/// `maximum_active_cells` bounds stored cells, plus four fixed child scratch
/// cells. Surface queries share one cell budget across all examined patches.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticQueryControls {
    pub distance_tolerance: f64,
    pub witness_tolerance: f64,
    pub maximum_cells: usize,
    pub maximum_active_cells: usize,
}

impl Default for QuadraticQueryControls {
    fn default() -> Self {
        Self {
            distance_tolerance: 1.0e-9,
            witness_tolerance: 1.0e-7,
            maximum_cells: 32_768,
            maximum_active_cells: 4_096,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QuadraticQueryError {
    Geometry(TriangleQueryError),
    Invalid(&'static str),
    DegeneratePatch,
    NumericalRange,
    ParameterResolution,
    CellBudgetExhausted { maximum_cells: usize },
    ActiveCellBudgetExhausted { maximum_active_cells: usize },
}

impl From<TriangleQueryError> for QuadraticQueryError {
    fn from(value: TriangleQueryError) -> Self {
        Self::Geometry(value)
    }
}
impl core::fmt::Display for QuadraticQueryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "quadratic triangle query refused: {self:?}")
    }
}
impl std::error::Error for QuadraticQueryError {}

/// Evaluation in original material coordinates. Coordinates are interpreted
/// as normalized nonnegative ratios; the six shape functions remain signed.
/// The normal is the actual local quadratic tangent cross product and follows
/// corner winding. This is floating geometry, not an orientation certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticTriangleEvaluation {
    pub point: [f64; 3],
    pub barycentric: [f64; 3],
    pub unit_normal: [f64; 3],
    pub feature: TriangleFeature,
}

/// One globally admissible witness, with enclosures of every global minimizer.
/// `distance` is the nominal relative-coordinate distance; subtraction of
/// rounded absolute `point` coordinates can differ at large translations.
/// `witness_error_bound` encloses distance from the returned floating point to
/// EVERY exact global minimizing point, including coordinate rounding.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticTrianglePoint {
    pub point: [f64; 3],
    pub barycentric: [f64; 3],
    pub distance: f64,
    pub distance_lower: f64,
    pub distance_upper: f64,
    pub unit_normal: [f64; 3],
    pub feature: TriangleFeature,
    pub witness_error_bound: f64,
    pub minimizer_bounds: [[f64; 2]; 3],
    pub minimizer_world_bounds: [[f64; 2]; 3],
    /// Enclosure of `(query-Q) dot (dQ/du cross dQ/dv)` at EVERY possible
    /// global minimizer, with the original winding (coordinate unit cubed).
    pub oriented_side_bounds: [f64; 2],
    pub examined_cells: usize,
}

/// Evaluate corners 0,1,2 then midside nodes 01,12,20. No corner-plane or
/// tessellation normal is substituted for the quadratic patch's normal.
pub fn evaluate_quadratic_triangle(
    nodes: [[f64; 3]; 6],
    barycentric: [f64; 3],
) -> Result<QuadraticTriangleEvaluation, QuadraticQueryError> {
    if !barycentric.iter().all(|x| x.is_finite() && *x >= 0.0)
        || (barycentric.iter().sum::<f64>() - 1.0).abs() > 64.0 * f64::EPSILON
    {
        return Err(QuadraticQueryError::Invalid(
            "closed-triangle barycentrics required",
        ));
    }
    let barycentric = normalized_weights(barycentric)?;
    let controls = relative_bernstein_controls(nodes, nodes[0])?;
    let approximate = controls.map(|p| p.map(Interval::midpoint));
    let point = nominal_point(nodes[0], approximate, barycentric)?;
    let [du, dv] = derivatives(approximate, barycentric);
    let scale = du.into_iter().chain(dv).map(f64::abs).fold(0.0, f64::max);
    if scale == 0.0 || !scale.is_finite() {
        return Err(QuadraticQueryError::DegeneratePatch);
    }
    let n = cross(du.map(|x| x / scale), dv.map(|x| x / scale));
    let area = norm(n)?;
    if area <= 64.0 * f64::EPSILON {
        return Err(QuadraticQueryError::DegeneratePatch);
    }
    Ok(QuadraticTriangleEvaluation {
        point,
        barycentric,
        unit_normal: n.map(|x| x / area),
        feature: feature(barycentric),
    })
}

/// Globally bound the closest point on the complete curved closed triangle.
/// Work/storage exhaustion and insufficient parameter resolution return no
/// partial nearest-point claim. Set a large finite witness tolerance for a
/// distance-only proof that does not require localized material coordinates.
pub fn closest_quadratic_triangle_point(
    point: [f64; 3],
    nodes: [[f64; 3]; 6],
    controls: QuadraticQueryControls,
    cx: &Cx<'_>,
) -> Result<QuadraticTrianglePoint, QuadraticQueryError> {
    validate_controls(controls)?;
    let mut budget = CellBudget {
        used: 0,
        maximum: controls.maximum_cells,
    };
    solve(point, nodes, controls, &mut budget, cx)
}

pub(super) fn validate_controls(
    controls: QuadraticQueryControls,
) -> Result<(), QuadraticQueryError> {
    if !controls.distance_tolerance.is_finite()
        || controls.distance_tolerance <= 0.0
        || !controls.witness_tolerance.is_finite()
        || controls.witness_tolerance <= 0.0
        || controls.maximum_cells == 0
        || controls.maximum_cells > MAXIMUM_CELLS
        || controls.maximum_active_cells == 0
        || controls.maximum_active_cells > MAXIMUM_ACTIVE_CELLS
        || controls.maximum_active_cells > controls.maximum_cells
    {
        return Err(QuadraticQueryError::Invalid(
            "positive finite tolerances and bounded cell/storage budgets required",
        ));
    }
    Ok(())
}

/// Exact Lagrange-to-Bernstein conversion enclosed relative to `anchor`.
/// Controls have corner/edge order 0,1,2,01,12,20 and enclose node minus anchor.
pub(super) fn relative_bernstein_controls(
    nodes: [[f64; 3]; 6],
    anchor: [f64; 3],
) -> Result<[[Interval; 3]; 6], QuadraticQueryError> {
    super::vector(anchor)?;
    for node in nodes {
        super::vector(node)?;
    }
    let relative = nodes.map(|node| {
        core::array::from_fn(|i| Interval::point(node[i]) - Interval::point(anchor[i]))
    });
    let mut controls = relative;
    for (edge, (a, b)) in EDGE_PAIRS.into_iter().enumerate() {
        for axis in 0..3 {
            controls[3 + edge][axis] = relative[3 + edge][axis] * Interval::point(2.0)
                - (relative[a][axis] + relative[b][axis]) * Interval::point(0.5);
        }
    }
    for c in controls {
        for x in c {
            finite_interval(x)?;
        }
    }
    Ok(controls)
}

pub(super) fn enclosed_norm_upper(v: [Interval; 3]) -> Result<f64, QuadraticQueryError> {
    let mut scale = 0.0_f64;
    for x in v {
        finite_interval(x)?;
        scale = scale.max(x.abs().hi());
    }
    if scale == 0.0 {
        return Ok(0.0);
    }
    let mut squared = Interval::point(0.0);
    for x in v {
        let y = (x / Interval::point(scale)).abs();
        squared = squared + y * y;
    }
    finite_interval(squared)?;
    let result = (squared.sqrt() * Interval::point(scale)).hi();
    if result.is_finite() {
        Ok(result)
    } else {
        Err(QuadraticQueryError::NumericalRange)
    }
}

fn finite_interval(x: Interval) -> Result<(), QuadraticQueryError> {
    if x.lo().is_finite() && x.hi().is_finite() {
        Ok(())
    } else {
        Err(QuadraticQueryError::NumericalRange)
    }
}

struct CellBudget {
    used: usize,
    maximum: usize,
}
impl CellBudget {
    fn consume(&mut self) -> Result<(), QuadraticQueryError> {
        if self.used >= self.maximum {
            return Err(QuadraticQueryError::CellBudgetExhausted {
                maximum_cells: self.maximum,
            });
        }
        self.used += 1;
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct Candidate {
    barycentric: [f64; 3],
    point: [f64; 3],
    relative: [Interval; 3],
    distance: f64,
    upper: f64,
}

#[derive(Clone)]
struct Cell {
    domain: [[f64; 3]; 3],
    controls: [[Interval; 3]; 6],
    lower: f64,
    id: usize,
    depth: u8,
}
impl PartialEq for Cell {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}
impl Eq for Cell {}
impl PartialOrd for Cell {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Cell {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        other
            .lower
            .total_cmp(&self.lower)
            .then_with(|| other.id.cmp(&self.id))
    }
}

struct Patch {
    nodes: [[f64; 3]; 6],
    point: [f64; 3],
    controls: [[Interval; 3]; 6],
    geometry: [[f64; 3]; 6],
}

fn solve(
    point: [f64; 3],
    nodes: [[f64; 3]; 6],
    controls: QuadraticQueryControls,
    budget: &mut CellBudget,
    cx: &Cx<'_>,
) -> Result<QuadraticTrianglePoint, QuadraticQueryError> {
    poll(cx)?;
    let before = budget.used;
    let patch = Patch {
        nodes,
        point,
        controls: relative_bernstein_controls(nodes, point)?,
        geometry: relative_bernstein_controls(nodes, nodes[0])?.map(|p| p.map(Interval::midpoint)),
    };
    let mut best = patch.candidate([1.0, 0.0, 0.0])?;
    let root = make_cell(
        &patch,
        CORNERS,
        patch.controls,
        0,
        0.0,
        budget,
        &mut best,
        cx,
    )?
    .ok_or(QuadraticQueryError::NumericalRange)?;
    let mut heap = BinaryHeap::with_capacity(controls.maximum_active_cells.min(16));
    heap.push(root);
    loop {
        poll(cx)?;
        while heap.peek().is_some_and(|cell| cell.lower > best.upper) {
            heap.pop();
        }
        let minimum = heap
            .peek()
            .ok_or(QuadraticQueryError::NumericalRange)?
            .lower;
        if gap(best.upper, minimum)? <= controls.distance_tolerance {
            break;
        }
        let cell = heap.pop().ok_or(QuadraticQueryError::NumericalRange)?;
        for child in split(&patch, &cell, budget, &mut best, cx)?
            .into_iter()
            .flatten()
        {
            if child.lower <= best.upper {
                if heap.len() == controls.maximum_active_cells {
                    return Err(QuadraticQueryError::ActiveCellBudgetExhausted {
                        maximum_active_cells: controls.maximum_active_cells,
                    });
                }
                if heap.len() == heap.capacity() {
                    heap.reserve_exact(
                        (heap.capacity().max(1)).min(controls.maximum_active_cells - heap.len()),
                    );
                }
                heap.push(child);
            }
        }
    }
    // After the distance proof, round-based refinement localizes ALL surviving
    // possible minimizers. No repeated global scan for every heap insertion,
    // and no tie-biased refinement of only one of several distant minima.
    let mut active = heap.into_vec();
    loop {
        poll(cx)?;
        active.retain(|cell| cell.lower <= best.upper);
        let enclosure = minimizer_enclosure(&patch, &active, best, cx)?;
        if enclosure.error <= controls.witness_tolerance {
            let evaluation = evaluate_quadratic_triangle(nodes, best.barycentric)?;
            return Ok(QuadraticTrianglePoint {
                point: best.point,
                barycentric: best.barycentric,
                distance: best.distance,
                distance_lower: enclosure.lower,
                distance_upper: best.upper,
                unit_normal: evaluation.unit_normal,
                feature: feature(best.barycentric),
                witness_error_bound: enclosure.error,
                minimizer_bounds: enclosure.parameter,
                minimizer_world_bounds: enclosure.world,
                oriented_side_bounds: enclosure.side,
                examined_cells: budget.used - before,
            });
        }
        let cutoff = budget.used;
        let mut index = 0;
        while index < active.len() {
            poll(cx)?;
            if active[index].id >= cutoff {
                index += 1;
                continue;
            }
            if active[index].lower > best.upper {
                active.swap_remove(index);
                continue;
            }
            if cell_witness_error(&active[index], point, best.point)? <= controls.witness_tolerance
            {
                index += 1;
                continue;
            }
            let cell = active.swap_remove(index);
            for child in split(&patch, &cell, budget, &mut best, cx)?
                .into_iter()
                .flatten()
            {
                if child.lower <= best.upper {
                    if active.len() == controls.maximum_active_cells {
                        return Err(QuadraticQueryError::ActiveCellBudgetExhausted {
                            maximum_active_cells: controls.maximum_active_cells,
                        });
                    }
                    if active.len() == active.capacity() {
                        active.reserve_exact(
                            (active.capacity().max(1))
                                .min(controls.maximum_active_cells - active.len()),
                        );
                    }
                    active.push(child);
                }
            }
        }
    }
}

impl Patch {
    fn candidate(&self, weights: [f64; 3]) -> Result<Candidate, QuadraticQueryError> {
        let barycentric = normalized_weights(weights)?;
        let enclosed = enclosed_weights(barycentric)?;
        let relative = blossom(self.controls, enclosed, enclosed);
        let upper = enclosed_norm_upper(relative)?;
        let point = nominal_point(self.nodes[0], self.geometry, barycentric)?;
        let distance = norm(relative.map(Interval::midpoint))?;
        Ok(Candidate {
            barycentric,
            point,
            relative,
            distance,
            upper,
        })
    }
}

fn make_cell(
    patch: &Patch,
    domain: [[f64; 3]; 3],
    controls: [[Interval; 3]; 6],
    depth: u8,
    inherited_lower: f64,
    budget: &mut CellBudget,
    best: &mut Candidate,
    cx: &Cx<'_>,
) -> Result<Option<Cell>, QuadraticQueryError> {
    poll(cx)?;
    let id = budget.used;
    budget.consume()?;
    for control in controls {
        for component in control {
            finite_interval(component)?;
        }
    }
    let approximate = controls.map(|p| p.map(Interval::midpoint));
    let mut local = [1.0 / 3.0; 3];
    let mut local_distance = f64::INFINITY;
    for seed in [
        CORNERS[0],
        CORNERS[1],
        CORNERS[2],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        local,
    ] {
        let d = norm(evaluate_nominal(approximate, seed))?;
        if d < local_distance {
            local = seed;
            local_distance = d;
        }
    }
    let linear =
        closest_triangle_point([0.0; 3], [approximate[0], approximate[1], approximate[2]]).ok();
    if let Some(linear) = linear {
        let d = norm(evaluate_nominal(approximate, linear.barycentric))?;
        if d < local_distance {
            local = linear.barycentric;
        }
    }
    local = improve_candidate(approximate, local, cx)?;
    let barycentric = core::array::from_fn(|axis| dot(local, domain.map(|v| v[axis])));
    let candidate = patch.candidate(barycentric)?;
    if candidate.upper < best.upper {
        *best = candidate;
    }
    let mut lower = hull_distance_lower(controls).max(inherited_lower);
    lower = lower.max(plane_lower(
        controls,
        candidate.relative.map(Interval::midpoint),
    )?);
    if let Some(linear) = linear {
        lower = lower.max(plane_lower(controls, linear.point)?);
    }
    if lower > best.upper || excludes_stationary_minimum(patch.controls, domain, controls)? {
        return Ok(None);
    }
    Ok(Some(Cell {
        domain,
        controls,
        lower,
        id,
        depth,
    }))
}

fn split(
    patch: &Patch,
    cell: &Cell,
    budget: &mut CellBudget,
    best: &mut Candidate,
    cx: &Cx<'_>,
) -> Result<[Option<Cell>; 4], QuadraticQueryError> {
    if cell.depth >= MAXIMUM_DEPTH {
        return Err(QuadraticQueryError::ParameterResolution);
    }
    let a = [0.5, 0.5, 0.0];
    let b = [0.0, 0.5, 0.5];
    let c = [0.5, 0.0, 0.5];
    let domains = [
        [CORNERS[0], a, c],
        [a, CORNERS[1], b],
        [c, b, CORNERS[2]],
        [a, b, c],
    ];
    let mut children: [Option<Cell>; 4] = core::array::from_fn(|_| None);
    for (slot, child) in children.iter_mut().zip(domains) {
        let domain = child
            .map(|weights| core::array::from_fn(|axis| dot(weights, cell.domain.map(|v| v[axis]))));
        let weights = child.map(|v| v.map(Interval::point));
        let mut controls = [blossom(cell.controls, weights[0], weights[0]); 6];
        for i in 0..3 {
            controls[i] = blossom(cell.controls, weights[i], weights[i]);
        }
        for (i, (a, b)) in EDGE_PAIRS.into_iter().enumerate() {
            controls[3 + i] = blossom(cell.controls, weights[a], weights[b]);
        }
        *slot = make_cell(
            patch,
            domain,
            controls,
            cell.depth + 1,
            cell.lower,
            budget,
            best,
            cx,
        )?;
    }
    Ok(children)
}

// A true global minimizer on a convex parameter triangle must have a
// nonnegative directional derivative toward each original parameter corner.
// This exclusion is valid even at original curved edges and vertices; no
// assumption of convex distance, local uniqueness, or a Newton basin is used.
fn excludes_stationary_minimum(
    original: [[Interval; 3]; 6],
    domain: [[f64; 3]; 3],
    controls: [[Interval; 3]; 6],
) -> Result<bool, QuadraticQueryError> {
    let residual = hull(controls);
    let [tangent_u, tangent_v] = tangent_hulls(original, domain);
    let rs = residual
        .map(|x| x.abs().hi())
        .into_iter()
        .fold(0.0, f64::max);
    let ts = tangent_u
        .into_iter()
        .chain(tangent_v)
        .map(|x| x.abs().hi())
        .fold(0.0, f64::max);
    if !rs.is_finite() || !ts.is_finite() {
        return Err(QuadraticQueryError::NumericalRange);
    }
    if rs == 0.0 || ts == 0.0 {
        return Ok(false);
    }
    let mut gu = Interval::point(0.0);
    let mut gv = Interval::point(0.0);
    for axis in 0..3 {
        let r = residual[axis] / Interval::point(rs);
        gu = gu + r * (tangent_u[axis] / Interval::point(ts));
        gv = gv + r * (tangent_v[axis] / Interval::point(ts));
    }
    finite_interval(gu)?;
    finite_interval(gv)?;
    let u = parameter_hull(domain, 1);
    let v = parameter_hull(domain, 2);
    for corner in CORNERS {
        let derivative =
            gu * (Interval::point(corner[1]) - u) + gv * (Interval::point(corner[2]) - v);
        if derivative.hi() < 0.0 {
            return Ok(true);
        }
    }
    Ok(false)
}

fn tangent_hulls(original: [[Interval; 3]; 6], domain: [[f64; 3]; 3]) -> [[Interval; 3]; 2] {
    let mut tangent_u = [Interval::point(0.0); 3];
    let mut tangent_v = tangent_u;
    for (i, vertex) in domain.into_iter().enumerate() {
        for axis in 0..3 {
            let u = ((original[3][axis] - original[0][axis]) * Interval::point(vertex[0])
                + (original[1][axis] - original[3][axis]) * Interval::point(vertex[1])
                + (original[4][axis] - original[5][axis]) * Interval::point(vertex[2]))
                * Interval::point(2.0);
            let v = ((original[5][axis] - original[0][axis]) * Interval::point(vertex[0])
                + (original[4][axis] - original[3][axis]) * Interval::point(vertex[1])
                + (original[2][axis] - original[5][axis]) * Interval::point(vertex[2]))
                * Interval::point(2.0);
            tangent_u[axis] = if i == 0 { u } else { tangent_u[axis].hull(u) };
            tangent_v[axis] = if i == 0 { v } else { tangent_v[axis].hull(v) };
        }
    }
    [tangent_u, tangent_v]
}

fn plane_lower(
    controls: [[Interval; 3]; 6],
    direction: [f64; 3],
) -> Result<f64, QuadraticQueryError> {
    let scale = direction.into_iter().map(f64::abs).fold(0.0, f64::max);
    if scale == 0.0 {
        return Ok(0.0);
    }
    let direction = direction.map(|x| Interval::point(x / scale));
    let length = enclosed_norm_upper(direction)?;
    let mut projection = f64::INFINITY;
    for point in controls {
        let mut p = Interval::point(0.0);
        for axis in 0..3 {
            p = p + point[axis] * direction[axis];
        }
        finite_interval(p)?;
        projection = projection.min(p.lo());
    }
    if projection <= 0.0 {
        return Ok(0.0);
    }
    Ok((Interval::point(projection) / Interval::point(length))
        .lo()
        .max(0.0))
}

fn hull_distance_lower(controls: [[Interval; 3]; 6]) -> f64 {
    hull(controls)
        .into_iter()
        .map(|x| {
            if x.lo() > 0.0 {
                x.lo()
            } else if x.hi() < 0.0 {
                -x.hi()
            } else {
                0.0
            }
        })
        .fold(0.0, f64::max)
}

fn hull(controls: [[Interval; 3]; 6]) -> [Interval; 3] {
    core::array::from_fn(|axis| {
        controls
            .iter()
            .skip(1)
            .fold(controls[0][axis], |h, p| h.hull(p[axis]))
    })
}

fn parameter_hull(domain: [[f64; 3]; 3], axis: usize) -> Interval {
    Interval::new(
        domain.iter().map(|x| x[axis]).fold(f64::INFINITY, f64::min),
        domain
            .iter()
            .map(|x| x[axis])
            .fold(f64::NEG_INFINITY, f64::max),
    )
}

fn blossom(controls: [[Interval; 3]; 6], a: [Interval; 3], b: [Interval; 3]) -> [Interval; 3] {
    core::array::from_fn(|axis| {
        let mut p = Interval::point(0.0);
        for i in 0..3 {
            p = p + controls[i][axis] * (a[i] * b[i]);
        }
        for (edge, (i, j)) in EDGE_PAIRS.into_iter().enumerate() {
            p = p + controls[3 + edge][axis] * (a[i] * b[j] + a[j] * b[i]);
        }
        p
    })
}

fn normalized_weights(weights: [f64; 3]) -> Result<[f64; 3], QuadraticQueryError> {
    let sum = weights.iter().sum::<f64>();
    if !sum.is_finite() || sum <= 0.0 || weights.iter().any(|x| !x.is_finite() || *x < 0.0) {
        return Err(QuadraticQueryError::NumericalRange);
    }
    Ok(weights.map(|x| x / sum))
}

fn enclosed_weights(weights: [f64; 3]) -> Result<[Interval; 3], QuadraticQueryError> {
    let sum = weights
        .into_iter()
        .fold(Interval::point(0.0), |s, x| s + Interval::point(x));
    if sum.lo() <= 0.0 {
        return Err(QuadraticQueryError::NumericalRange);
    }
    let mut result = [Interval::point(0.0); 3];
    for i in 0..3 {
        result[i] = (Interval::point(weights[i]) / sum)
            .intersect(Interval::new(0.0, 1.0))
            .ok_or(QuadraticQueryError::NumericalRange)?;
    }
    Ok(result)
}

fn evaluate_nominal(controls: [[f64; 3]; 6], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|axis| {
        let mut p = 0.0;
        for i in 0..3 {
            p = (b[i] * b[i]).mul_add(controls[i][axis], p);
        }
        for (edge, (i, j)) in EDGE_PAIRS.into_iter().enumerate() {
            p = (2.0 * b[i] * b[j]).mul_add(controls[3 + edge][axis], p);
        }
        p
    })
}

fn nominal_point(
    anchor: [f64; 3],
    controls: [[f64; 3]; 6],
    b: [f64; 3],
) -> Result<[f64; 3], QuadraticQueryError> {
    let relative = evaluate_nominal(controls, b);
    Ok(super::vector(core::array::from_fn(|axis| {
        anchor[axis] + relative[axis]
    }))?)
}

fn derivatives(c: [[f64; 3]; 6], b: [f64; 3]) -> [[f64; 3]; 2] {
    [
        core::array::from_fn(|i| {
            2.0 * ((c[3][i] - c[0][i]) * b[0]
                + (c[1][i] - c[3][i]) * b[1]
                + (c[4][i] - c[5][i]) * b[2])
        }),
        core::array::from_fn(|i| {
            2.0 * ((c[5][i] - c[0][i]) * b[0]
                + (c[4][i] - c[3][i]) * b[1]
                + (c[2][i] - c[5][i]) * b[2])
        }),
    ]
}

// Fixed bounded projected Gauss-Newton iterations; these merely improve a
// feasible upper witness. Neither stationarity nor convergence admits a result.
fn improve_candidate(
    controls: [[f64; 3]; 6],
    mut b: [f64; 3],
    cx: &Cx<'_>,
) -> Result<[f64; 3], QuadraticQueryError> {
    let scale = controls
        .into_iter()
        .flatten()
        .map(f64::abs)
        .fold(0.0, f64::max);
    if scale == 0.0 {
        return Ok(b);
    }
    let c = controls.map(|p| p.map(|x| x / scale));
    for _ in 0..12 {
        poll(cx)?;
        let r = evaluate_nominal(c, b);
        let [u, v] = derivatives(c, b);
        let a = dot(u, u);
        let d = dot(v, v);
        let e = dot(u, v);
        let gu = dot(r, u);
        let gv = dot(r, v);
        let determinant = a.mul_add(d, -e * e);
        let step = if determinant > 64.0 * f64::EPSILON * a * d {
            [
                (d * gu - e * gv) / determinant,
                (a * gv - e * gu) / determinant,
            ]
        } else if a.max(d) > 0.0 {
            [gu / a.max(d), gv / a.max(d)]
        } else {
            break;
        };
        if step.iter().any(|x| !x.is_finite()) {
            break;
        }
        let distance = norm(r)?;
        let mut improved = false;
        for exponent in 0..8 {
            let alpha = 0.5_f64.powi(exponent);
            let next = project_parameter(b[1] - alpha * step[0], b[2] - alpha * step[1]);
            if norm(evaluate_nominal(c, next))? < distance {
                b = next;
                improved = true;
                break;
            }
        }
        if !improved {
            break;
        }
    }
    Ok(b)
}

fn project_parameter(u: f64, v: f64) -> [f64; 3] {
    let mut u = u.max(0.0);
    let mut v = v.max(0.0);
    if u + v > 1.0 {
        u = (0.5 * (u - v + 1.0)).clamp(0.0, 1.0);
        v = 1.0 - u;
    }
    [(1.0 - u - v).max(0.0), u, v]
}

fn feature(b: [f64; 3]) -> TriangleFeature {
    for i in 0..3 {
        if b[(i + 1) % 3] == 0.0 && b[(i + 2) % 3] == 0.0 {
            return TriangleFeature::Vertex(i);
        }
    }
    if b[2] == 0.0 {
        TriangleFeature::Edge(0, 1)
    } else if b[1] == 0.0 {
        TriangleFeature::Edge(0, 2)
    } else if b[0] == 0.0 {
        TriangleFeature::Edge(1, 2)
    } else {
        TriangleFeature::Face
    }
}

fn gap(upper: f64, lower: f64) -> Result<f64, QuadraticQueryError> {
    let width = Interval::point(upper) - Interval::point(lower);
    finite_interval(width)?;
    Ok(width.hi().max(0.0))
}

fn cell_witness_error(
    cell: &Cell,
    query: [f64; 3],
    witness: [f64; 3],
) -> Result<f64, QuadraticQueryError> {
    let relative: [Interval; 3] =
        core::array::from_fn(|i| Interval::point(witness[i]) - Interval::point(query[i]));
    let mut radius = 0.0_f64;
    for control in cell.controls {
        radius = radius.max(enclosed_norm_upper(core::array::from_fn(|i| {
            control[i] - relative[i]
        }))?);
    }
    Ok(radius)
}

struct MinimizerEnclosure {
    lower: f64,
    error: f64,
    parameter: [[f64; 2]; 3],
    world: [[f64; 2]; 3],
    side: [f64; 2],
}
fn minimizer_enclosure(
    patch: &Patch,
    cells: &[Cell],
    best: Candidate,
    cx: &Cx<'_>,
) -> Result<MinimizerEnclosure, QuadraticQueryError> {
    if cells.is_empty() {
        return Err(QuadraticQueryError::NumericalRange);
    }
    let mut lower = best.upper;
    let mut error = 0.0_f64;
    let mut parameter = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
    let mut world = parameter;
    let mut side = [f64::INFINITY, f64::NEG_INFINITY];
    for cell in cells {
        poll(cx)?;
        lower = lower.min(cell.lower);
        error = error.max(cell_witness_error(cell, patch.point, best.point)?);
        let bounds = hull(cell.controls);
        let [du, dv] = tangent_hulls(patch.controls, cell.domain);
        let normal: [Interval; 3] = core::array::from_fn(|i| {
            du[(i + 1) % 3] * dv[(i + 2) % 3] - du[(i + 2) % 3] * dv[(i + 1) % 3]
        });
        let signed = (0..3).fold(Interval::point(0.0), |sum, i| sum - bounds[i] * normal[i]);
        finite_interval(signed)?;
        side[0] = side[0].min(signed.lo());
        side[1] = side[1].max(signed.hi());
        for axis in 0..3 {
            let p = parameter_hull(cell.domain, axis);
            parameter[axis][0] = parameter[axis][0].min(p.lo());
            parameter[axis][1] = parameter[axis][1].max(p.hi());
            let w = bounds[axis] + Interval::point(patch.point[axis]);
            finite_interval(w)?;
            world[axis][0] = world[axis][0].min(w.lo());
            world[axis][1] = world[axis][1].max(w.hi());
        }
    }
    Ok(MinimizerEnclosure {
        lower,
        error,
        parameter,
        world,
        side,
    })
}

#[cfg(test)]
mod tests;
