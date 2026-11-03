//! Positive-Jacobian admission over a closed TET10 simplex and affine nodal path.
//!
//! Corners are 0,1,2,3; midsides follow 01,12,02,03,13,23. The complete
//! quadratic geometry is bounded, including faces, edges, and vertices. An
//! interval determinant is tried first; cubic simplex × cubic time Bernstein
//! coefficients retain dependencies when that enclosure is insufficient.
//! Bounded subdivision refines the material tetrahedron and fractional time.
//! Point evaluations can refuse a path but never admit one.
//!
//! This proves local positive orientation along the declared affine motion of
//! all ten nodes. It does not prove global injectivity, absence of intersections
//! between distinct material points, or the path of an underlying time solver.

use fs_exec::Cx;
use fs_ivl::Interval;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

const MAXIMUM_CELLS: usize = 1 << 22;
const MAXIMUM_ACTIVE_CELLS: usize = 1 << 18;
const MAXIMUM_DEPTH: u8 = 52;
const EDGES: [[usize; 2]; 6] = [[0, 1], [1, 2], [0, 2], [0, 3], [1, 3], [2, 3]];
const SIMPLEX: [[f64; 4]; 4] = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
];
type Matrix = [[Interval; 3]; 3];
type JacobianControls = [[Matrix; 2]; 4];

/// Work, storage, fractional-time, material-coordinate, and quality limits.
/// `maximum_cells` counts evaluated cells, including cells subsequently split.
/// `maximum_active_cells` bounds the pending frontier; scratch is fixed size.
/// `parameter_tolerance` is the maximum barycentric-coordinate span used as a
/// spatial refinement floor. Unresolved cells can remain spatially broad when
/// temporal variation already reaches its fractional-time floor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticTetPathControls {
    pub maximum_cells: usize,
    pub maximum_active_cells: usize,
    pub fraction_tolerance: f64,
    pub parameter_tolerance: f64,
    /// Strictly positive admission requires a lower bound greater than this
    /// nonnegative threshold, in the normalization reported by the query.
    pub minimum_scaled_determinant: f64,
}

impl Default for QuadraticTetPathControls {
    fn default() -> Self {
        Self {
            maximum_cells: 1 << 18,
            maximum_active_cells: 4096,
            fraction_tolerance: 1.0e-6,
            parameter_tolerance: 1.0e-5,
            minimum_scaled_determinant: 1.0e-12,
        }
    }
}

impl QuadraticTetPathControls {
    pub fn validate(self) -> Result<(), QuadraticTetPathError> {
        if self.maximum_cells == 0
            || self.maximum_cells > MAXIMUM_CELLS
            || self.maximum_active_cells == 0
            || self.maximum_active_cells > MAXIMUM_ACTIVE_CELLS
            || !self.fraction_tolerance.is_finite()
            || self.fraction_tolerance <= 0.0
            || self.fraction_tolerance >= 1.0
            || !self.parameter_tolerance.is_finite()
            || self.parameter_tolerance <= 0.0
            || self.parameter_tolerance >= 1.0
            || !self.minimum_scaled_determinant.is_finite()
            || self.minimum_scaled_determinant < 0.0
        {
            return Err(QuadraticTetPathError::Invalid(
                "bounded positive work/storage, tolerances in (0,1), and finite nonnegative determinant threshold required",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QuadraticTetPathError {
    Invalid(&'static str),
    Cancelled,
    DegenerateGeometry,
    NumericalRange,
    ParameterResolution,
    CellBudgetExhausted { maximum_cells: usize },
    ActiveCellBudgetExhausted { maximum_active_cells: usize },
}

impl core::fmt::Display for QuadraticTetPathError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "quadratic tetrahedron path query refused: {self:?}")
    }
}
impl std::error::Error for QuadraticTetPathError {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum QuadraticTetPathVerdict {
    /// Every closed-domain material point at every fraction in [0,1] has a
    /// scaled determinant at least this lower bound, strictly above the limit.
    Positive { minimum_scaled_determinant: f64 },
    /// The earliest uncertified window. Fractions strictly before its start
    /// have been admitted over the entire material simplex. This is neither a
    /// collision time nor a claim that an actual fold exists. Bounds refer to
    /// one unresolved material subdomain (or a demonstrated point refusal).
    Unresolved {
        fraction_start: f64,
        fraction_end: f64,
        barycentric_bounds: [[f64; 2]; 4],
        scaled_determinant_bounds: [f64; 2],
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticTetPathReport {
    pub verdict: QuadraticTetPathVerdict,
    pub examined_cells: usize,
    /// Let S be this positive finite power of two. At each path endpoint,
    /// subtract that endpoint's node 0 with directed arithmetic, take the
    /// largest absolute coordinate bound over all ten nodes and both ends,
    /// and round upward to a power of two. Reported determinant bounds enclose
    /// det(dX/d(lambda1,lambda2,lambda3)) / S^3. S^3 is never formed in f64.
    /// This is a per-element/path scale, not a corner-edge quality metric.
    pub normalization_length: f64,
}

/// Admit only a strictly positive Jacobian over the complete closed material
/// tetrahedron and every fraction of the declared affine ten-node path.
/// Inputs are finite floating coordinates interpreted as exact real numbers.
/// Endpoint-relative normalization removes common translation before any
/// determinant arithmetic. An unresolved report and all typed errors refuse
/// full-path admission; work or storage exhaustion never becomes `Positive`.
pub fn quadratic_tetrahedron_orientation(
    start: [[f64; 3]; 10],
    end: [[f64; 3]; 10],
    controls: QuadraticTetPathControls,
    cx: &Cx<'_>,
) -> Result<QuadraticTetPathReport, QuadraticTetPathError> {
    controls.validate()?;
    poll(cx)?;
    let (jacobians, normalization_length) = normalized_jacobians(start, end)?;
    let static_path = start == end;
    let mut frontier = BinaryHeap::with_capacity(controls.maximum_active_cells.min(16));
    frontier.push(Cell {
        vertices: SIMPLEX,
        fraction: [0.0, 1.0],
        spatial_depth: 0,
        time_depth: 0,
        serial: 0,
    });
    let mut examined_cells = 0;
    let mut next_serial = 1;
    let mut minimum = f64::INFINITY;
    while let Some(cell) = frontier.pop() {
        poll(cx)?;
        if examined_cells == controls.maximum_cells {
            return Err(QuadraticTetPathError::CellBudgetExhausted {
                maximum_cells: controls.maximum_cells,
            });
        }
        examined_cells += 1;
        let restricted = restrict(&jacobians, &cell);
        let direct = finite_interval(determinant(&matrix_hull(&restricted)))?;
        if direct.lo() > controls.minimum_scaled_determinant {
            minimum = minimum.min(direct.lo());
            continue;
        }
        // A rigorous point counterexample can terminate a proof search. It
        // cannot certify the rest of the element or move past an earlier cell.
        for (vertex, matrices) in restricted.iter().enumerate() {
            let at_start = finite_interval(determinant(&matrices[0]))?;
            if at_start.hi() <= controls.minimum_scaled_determinant {
                let mut point_cell = cell;
                point_cell.vertices = [cell.vertices[vertex]; 4];
                point_cell.fraction[1] = point_cell.fraction[0];
                return Ok(unresolved(
                    point_cell,
                    at_start,
                    examined_cells,
                    normalization_length,
                ));
            }
        }
        let bernstein = determinant_bernstein_hull(&restricted, cx)?;
        let bound = direct
            .intersect(bernstein)
            .ok_or(QuadraticTetPathError::NumericalRange)?;
        if bound.lo() > controls.minimum_scaled_determinant {
            minimum = minimum.min(bound.lo());
            continue;
        }
        if bound.hi() <= controls.minimum_scaled_determinant {
            let mut first = cell;
            first.fraction[1] = first.fraction[0];
            return Ok(unresolved(
                first,
                bound,
                examined_cells,
                normalization_length,
            ));
        }
        let (spatial_variation, time_variation) = variations(&restricted);
        let spatial_small = cell.parameter_span() <= controls.parameter_tolerance;
        let time_small =
            static_path || cell.fraction[1] - cell.fraction[0] <= controls.fraction_tolerance;
        if time_small && (spatial_small || (!static_path && time_variation >= spatial_variation)) {
            let mut first = cell;
            if static_path {
                first.fraction[1] = first.fraction[0];
            }
            return Ok(unresolved(
                first,
                bound,
                examined_cells,
                normalization_length,
            ));
        }
        let split_time =
            !static_path && !time_small && (spatial_small || time_variation >= spatial_variation);
        let children = cell.split(split_time, next_serial)?;
        next_serial += 2;
        if frontier.len() + 2 > controls.maximum_active_cells {
            return Err(QuadraticTetPathError::ActiveCellBudgetExhausted {
                maximum_active_cells: controls.maximum_active_cells,
            });
        }
        if frontier.len() + 2 > frontier.capacity() {
            let capacity = (2 * frontier.capacity())
                .max(frontier.len() + 2)
                .min(controls.maximum_active_cells);
            frontier.reserve_exact(capacity - frontier.len());
        }
        frontier.extend(children);
    }
    Ok(QuadraticTetPathReport {
        verdict: QuadraticTetPathVerdict::Positive {
            minimum_scaled_determinant: minimum,
        },
        examined_cells,
        normalization_length,
    })
}

fn poll(cx: &Cx<'_>) -> Result<(), QuadraticTetPathError> {
    cx.checkpoint()
        .map_err(|_| QuadraticTetPathError::Cancelled)
}

fn finite_interval(value: Interval) -> Result<Interval, QuadraticTetPathError> {
    if value.lo().is_finite() && value.hi().is_finite() {
        Ok(value)
    } else {
        Err(QuadraticTetPathError::NumericalRange)
    }
}

fn power_of_two_ceiling(value: f64) -> Result<f64, QuadraticTetPathError> {
    if value == 0.0 {
        return Err(QuadraticTetPathError::DegenerateGeometry);
    }
    if !value.is_finite() || value < 0.0 {
        return Err(QuadraticTetPathError::NumericalRange);
    }
    let bits = value.to_bits();
    let exponent = bits >> 52;
    let fraction = bits & ((1_u64 << 52) - 1);
    let rounded = if exponent == 0 {
        // Subnormal powers of two are powers of two in the fraction bits.
        let shift = 63 - fraction.leading_zeros();
        let lower = 1_u64 << shift;
        if fraction == lower { lower } else { lower << 1 }
    } else {
        (exponent + u64::from(fraction != 0)) << 52
    };
    let scale = f64::from_bits(rounded);
    if scale.is_finite() {
        Ok(scale)
    } else {
        Err(QuadraticTetPathError::NumericalRange)
    }
}

fn normalized_jacobians(
    start: [[f64; 3]; 10],
    end: [[f64; 3]; 10],
) -> Result<(JacobianControls, f64), QuadraticTetPathError> {
    let zero = Interval::point(0.0);
    let mut relative = [[[zero; 3]; 10]; 2];
    let mut largest = 0.0_f64;
    for (time, nodes) in [start, end].iter().enumerate() {
        for (node, position) in nodes.iter().enumerate() {
            for axis in 0..3 {
                if !position[axis].is_finite() || !nodes[0][axis].is_finite() {
                    return Err(QuadraticTetPathError::Invalid(
                        "finite TET10 nodes required",
                    ));
                }
                let offset = if position[axis] == nodes[0][axis] {
                    zero
                } else {
                    finite_interval(
                        Interval::point(position[axis]) - Interval::point(nodes[0][axis]),
                    )?
                };
                largest = largest.max(offset.lo().abs()).max(offset.hi().abs());
                relative[time][node][axis] = offset;
            }
        }
    }
    let length = power_of_two_ceiling(largest)?;
    let scale = Interval::point(length);
    let mut jacobians = [[[[zero; 3]; 3]; 2]; 4];
    for (time, nodes) in relative.iter_mut().enumerate() {
        for position in nodes.iter_mut() {
            for coordinate in position {
                *coordinate = finite_interval(*coordinate / scale)?;
            }
        }
        let mut coefficients = [[[zero; 3]; 4]; 4];
        for corner in 0..4 {
            coefficients[corner][corner] = nodes[corner];
        }
        for (edge, [a, b]) in EDGES.iter().copied().enumerate() {
            for axis in 0..3 {
                let value = Interval::point(2.0) * nodes[4 + edge][axis]
                    - Interval::point(0.5) * (nodes[a][axis] + nodes[b][axis]);
                coefficients[a][b][axis] = value;
                coefficients[b][a][axis] = value;
            }
        }
        for vertex in 0..4 {
            for axis in 0..3 {
                for column in 0..3 {
                    jacobians[vertex][time][axis][column] = finite_interval(
                        Interval::point(2.0)
                            * (coefficients[vertex][column + 1][axis]
                                - coefficients[vertex][0][axis]),
                    )?;
                }
            }
        }
    }
    Ok((jacobians, length))
}

#[derive(Clone, Copy, Debug)]
struct Cell {
    vertices: [[f64; 4]; 4],
    fraction: [f64; 2],
    spatial_depth: u8,
    time_depth: u8,
    serial: usize,
}

impl PartialEq for Cell {
    fn eq(&self, other: &Self) -> bool {
        self.serial == other.serial
    }
}
impl Eq for Cell {}
impl PartialOrd for Cell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Cell {
    fn cmp(&self, other: &Self) -> Ordering {
        other.fraction[0]
            .total_cmp(&self.fraction[0])
            .then_with(|| other.fraction[1].total_cmp(&self.fraction[1]))
            .then_with(|| self.spatial_depth.cmp(&other.spatial_depth))
            .then_with(|| other.serial.cmp(&self.serial))
    }
}

impl Cell {
    fn parameter_bounds(self) -> [[f64; 2]; 4] {
        core::array::from_fn(|axis| {
            let mut bounds = [f64::INFINITY, f64::NEG_INFINITY];
            for vertex in self.vertices {
                bounds[0] = bounds[0].min(vertex[axis]);
                bounds[1] = bounds[1].max(vertex[axis]);
            }
            bounds
        })
    }

    fn parameter_span(self) -> f64 {
        self.parameter_bounds()
            .iter()
            .map(|bounds| bounds[1] - bounds[0])
            .fold(0.0, f64::max)
    }

    fn split(self, time: bool, serial: usize) -> Result<[Self; 2], QuadraticTetPathError> {
        let mut children = [self; 2];
        children[0].serial = serial;
        children[1].serial = serial + 1;
        if time {
            if self.time_depth == MAXIMUM_DEPTH {
                return Err(QuadraticTetPathError::ParameterResolution);
            }
            let midpoint = self.fraction[0].midpoint(self.fraction[1]);
            if midpoint <= self.fraction[0] || midpoint >= self.fraction[1] {
                return Err(QuadraticTetPathError::ParameterResolution);
            }
            children[0].fraction[1] = midpoint;
            children[1].fraction[0] = midpoint;
            for child in &mut children {
                child.time_depth += 1;
            }
        } else {
            if self.spatial_depth == MAXIMUM_DEPTH {
                return Err(QuadraticTetPathError::ParameterResolution);
            }
            let mut longest = -1.0;
            let mut edge = EDGES[0];
            for [a, b] in EDGES {
                let squared = (0..4)
                    .map(|axis| (self.vertices[a][axis] - self.vertices[b][axis]).powi(2))
                    .sum::<f64>();
                if squared > longest {
                    longest = squared;
                    edge = [a, b];
                }
            }
            let [a, b] = edge;
            // At depth <=52, these [0,1] dyadic barycentrics and their
            // midpoints are exact f64 values. The children cover the closed
            // parent tetrahedron without an omitted rounding sliver.
            let midpoint = core::array::from_fn(|axis| {
                self.vertices[a][axis].midpoint(self.vertices[b][axis])
            });
            children[0].vertices[a] = midpoint;
            children[1].vertices[b] = midpoint;
            for child in &mut children {
                child.spatial_depth += 1;
            }
        }
        Ok(children)
    }
}

fn restrict(jacobians: &JacobianControls, cell: &Cell) -> JacobianControls {
    let zero = Interval::point(0.0);
    let mut output = [[[[zero; 3]; 3]; 2]; 4];
    for (vertex, weights) in cell.vertices.iter().enumerate() {
        for row in 0..3 {
            for column in 0..3 {
                let mut endpoints = [zero; 2];
                for time in 0..2 {
                    for original in 0..4 {
                        if weights[original] != 0.0 {
                            let value = jacobians[original][time][row][column];
                            endpoints[time] = endpoints[time]
                                + if weights[original] == 1.0 {
                                    value
                                } else {
                                    Interval::point(weights[original]) * value
                                };
                        }
                    }
                }
                for (time, fraction) in cell.fraction.iter().copied().enumerate() {
                    output[vertex][time][row][column] = if fraction == 0.0 {
                        endpoints[0]
                    } else if fraction == 1.0 {
                        endpoints[1]
                    } else {
                        (Interval::point(1.0) - Interval::point(fraction)) * endpoints[0]
                            + Interval::point(fraction) * endpoints[1]
                    };
                }
            }
        }
    }
    output
}

fn matrix_hull(jacobians: &JacobianControls) -> Matrix {
    let mut hull = jacobians[0][0];
    for vertex in jacobians {
        for matrix in vertex {
            for row in 0..3 {
                for column in 0..3 {
                    hull[row][column] = hull[row][column].hull(matrix[row][column]);
                }
            }
        }
    }
    hull
}

fn determinant(matrix: &Matrix) -> Interval {
    determinant_columns(matrix, matrix, matrix)
}

fn determinant_columns(a: &Matrix, b: &Matrix, c: &Matrix) -> Interval {
    a[0][0] * (b[1][1] * c[2][2] - b[2][1] * c[1][2])
        - a[1][0] * (b[0][1] * c[2][2] - b[2][1] * c[0][2])
        + a[2][0] * (b[0][1] * c[1][2] - b[1][1] * c[0][2])
}

fn coefficient_index(mut triple: [usize; 3]) -> (usize, usize) {
    if triple[0] > triple[1] {
        triple.swap(0, 1);
    }
    if triple[1] > triple[2] {
        triple.swap(1, 2);
    }
    if triple[0] > triple[1] {
        triple.swap(0, 1);
    }
    let [a, b, c] = triple;
    let index = [0, 10, 16, 19][a] + (a..b).map(|i| 4 - i).sum::<usize>() + c - b;
    let multiplicity = if a == c {
        1
    } else if a == b || b == c {
        3
    } else {
        6
    };
    (index, multiplicity)
}

fn determinant_bernstein_hull(
    jacobians: &JacobianControls,
    cx: &Cx<'_>,
) -> Result<Interval, QuadraticTetPathError> {
    let mut coefficients = [[Interval::point(0.0); 4]; 20];
    let mut multiplicities = [0; 20];
    // The determinant is multilinear in its three columns. Group their
    // 4^3 × 2^3 products by simplex multiindex and time degree, then divide
    // by the exact multinomial/binomial integers. Inexact 1/3 or 1/6
    // constants are never interpreted as exact Bernstein weights.
    for a in 0..4 {
        poll(cx)?;
        for b in 0..4 {
            for c in 0..4 {
                let (index, multiplicity) = coefficient_index([a, b, c]);
                multiplicities[index] = multiplicity;
                for time_mask in 0_usize..8 {
                    let degree = time_mask.count_ones() as usize;
                    coefficients[index][degree] = coefficients[index][degree]
                        + determinant_columns(
                            &jacobians[a][time_mask & 1],
                            &jacobians[b][(time_mask >> 1) & 1],
                            &jacobians[c][time_mask >> 2],
                        );
                }
            }
        }
    }
    let mut hull: Option<Interval> = None;
    for (index, row) in coefficients.iter().enumerate() {
        for (time_degree, value) in row.iter().copied().enumerate() {
            let divisor = multiplicities[index] * [1, 3, 3, 1][time_degree];
            let bound = finite_interval(value / Interval::point(divisor as f64))?;
            hull = Some(hull.map_or(bound, |previous| previous.hull(bound)));
        }
    }
    hull.ok_or(QuadraticTetPathError::NumericalRange)
}

fn variations(jacobians: &JacobianControls) -> (f64, f64) {
    // These nominal variations choose a subdivision axis only. Neither
    // their accuracy nor their sign enters any admission decision.
    let mut spatial = 0.0_f64;
    let mut temporal = 0.0_f64;
    for vertex in 0..4 {
        for row in 0..3 {
            for column in 0..3 {
                temporal = temporal.max(
                    (jacobians[vertex][0][row][column].midpoint()
                        - jacobians[vertex][1][row][column].midpoint())
                    .abs(),
                );
                for time in 0..2 {
                    for other in 0..vertex {
                        spatial = spatial.max(
                            (jacobians[vertex][time][row][column].midpoint()
                                - jacobians[other][time][row][column].midpoint())
                            .abs(),
                        );
                    }
                }
            }
        }
    }
    (spatial, temporal)
}

fn unresolved(
    cell: Cell,
    bound: Interval,
    examined_cells: usize,
    normalization_length: f64,
) -> QuadraticTetPathReport {
    QuadraticTetPathReport {
        verdict: QuadraticTetPathVerdict::Unresolved {
            fraction_start: cell.fraction[0],
            fraction_end: cell.fraction[1],
            barycentric_bounds: cell.parameter_bounds(),
            scaled_determinant_bounds: [bound.lo(), bound.hi()],
        },
        examined_cells,
        normalization_length,
    }
}

#[cfg(test)]
mod tests;
