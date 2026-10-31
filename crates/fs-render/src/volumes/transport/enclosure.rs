//! Outward enclosures of the original piecewise-constant physical grid.
//!
//! A segment is parameterized by its two retained binary64 endpoints, not an
//! already-rounded unit direction. Cell planes are the real expressions
//! `origin + integer_index * cell_size`. An implicit binary hierarchy prunes
//! disjoint index boxes; only possible leaf cells contribute column mass.
//! No density copy, second medium, precomputed voxel tree, or ray allocation is
//! needed. Interval clipping bounds every cell length independently, so ties
//! need not be resolved by an uncertified DDA ordering.
//!
//! Cells are half-open on their upper faces. A ray possibly lying along a face
//! includes both ownership possibilities: uncertain cells contribute only to
//! the upper column-mass bound. This can be loose, but never invents a clear
//! path. Calibration, field discretization and model discrepancy are excluded.

use super::{AbsorbingGrid, MassAbsorptionSpectrum, VolumeTransportError, checkpoint};
use fs_exec::Cx;

mod prepared;
pub use prepared::{MAX_PREPARED_EXTINCTION_VALUES, PreparedGridExtinction};
use fs_geom::Point3;
use fs_ivl::{Interval, decay::exp_neg};

/// Version of the retained cell model, clipping and attenuation enclosure.
pub const GRID_EXTINCTION_ENCLOSURE_SEMANTICS_VERSION: u32 = 1;

/// Additional bound on optional interval traversal work. The original grid's
/// `max_steps_per_segment` also bounds the number of possible leaf cells.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GridEnclosureBudget {
    /// Maximum index boxes tested, including boxes rejected by clipping.
    pub max_nodes_per_segment: usize,
}

impl Default for GridEnclosureBudget {
    fn default() -> Self {
        Self { max_nodes_per_segment: 1_048_576 }
    }
}

/// Column density of the retained numerical field, not of an unknown specimen.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ColumnDensityEnclosure {
    /// Enclosure of `integral rho ds`, in kg/m². Outside the grid is vacuum.
    pub kg_per_m2: Interval,
    /// Number of tested hierarchy boxes; not a wall-clock performance claim.
    pub visited_nodes: usize,
    /// Cells with possibly positive intersection length, including empty cells.
    pub possible_cells: u32,
    /// Possible cells whose parallel-face ownership was not decided.
    pub uncertain_parallel_cells: u32,
}

impl MassAbsorptionSpectrum {
    /// Enclose the retained affine spectral law, never extrapolating its domain.
    pub(crate) fn coefficient_enclosure(&self, wavelength: f64) -> Result<Interval, VolumeTransportError> {
        let (lower, upper) = self.wavelength_bounds_nm();
        if !wavelength.is_finite() || !(lower..=upper).contains(&wavelength) {
            return Err(VolumeTransportError::Wavelength);
        }
        let knots = &self.knots[..self.count];
        let upper = knots.partition_point(|&(lambda, _)| lambda < wavelength);
        if upper < knots.len() && knots[upper].0 == wavelength {
            return Ok(Interval::point(knots[upper].1));
        }
        let (left, a) = knots[upper - 1];
        let (right, b) = knots[upper];
        if a == b { return Ok(Interval::point(a)); }
        let t = (Interval::point(wavelength) - Interval::point(left))
            / (Interval::point(right) - Interval::point(left));
        let t = Interval::new(t.lo().max(0.0), t.hi().min(1.0));
        let u = Interval::point(1.0) - t;
        let u = Interval::new(u.lo().max(0.0), u.hi().min(1.0));
        let value = Interval::point(a) * u + Interval::point(b) * t;
        Ok(Interval::new(value.lo().max(a.min(b)), value.hi().min(a.max(b))))
    }
}

impl AbsorbingGrid<'_> {
    /// Enclose the retained absorption interpolant at the requested wavelengths.
    /// Indexed optical palettes return `Palette`; use `prepare_extinction_enclosure`
    /// for their cell-selected laws. Preparation may be shared across a whole grid
    /// of receivers. Missing source wavelengths refuse; no extrapolation occurs.
    pub fn extinction_coefficients_enclosure<const N: usize>(
        &self,
        wavelengths_nm: &[f64; N],
    ) -> Result<[Interval; N], VolumeTransportError> {
        if N == 0 { return Err(VolumeTransportError::Sample); }
        if self.absorption_palette().is_some() {
            // An indexed field has no representative global spectrum.
            return Err(VolumeTransportError::Palette);
        }
        let mut result = [Interval::point(0.0); N];
        for (lane, wavelength) in wavelengths_nm.iter().enumerate() {
            result[lane] = self.absorption.coefficient_enclosure(*wavelength)?;
        }
        Ok(result)
    }

    /// Bound column density along a finite endpoint-to-endpoint segment.
    ///
    /// The original field remains borrowed and unchanged. The fixed stack has
    /// room for every balanced split of a `usize`-addressable cell product.
    /// Cancellation, exhausted work and unrepresentable arithmetic return no
    /// partial enclosure. A zero-length segment has exact zero column mass.
    pub fn column_density_enclosure(
        &self,
        cx: &Cx<'_>,
        start: Point3,
        end: Point3,
        budget: GridEnclosureBudget,
    ) -> Result<ColumnDensityEnclosure, VolumeTransportError> {
        self.visit_column(cx, start, end, budget, |_, _| Ok(()))
    }

    // One traversal owner for column density and cell-specific optical depth.
    // The callback sees an outward physical column contribution for each cell.
    fn visit_column(
        &self,
        cx: &Cx<'_>,
        start: Point3,
        end: Point3,
        budget: GridEnclosureBudget,
        mut visit: impl FnMut(usize, Interval) -> Result<(), VolumeTransportError>,
    ) -> Result<ColumnDensityEnclosure, VolumeTransportError> {
        checkpoint(cx)?;
        if budget.max_nodes_per_segment == 0 {
            return Err(VolumeTransportError::Budget);
        }
        let a = [start.x, start.y, start.z];
        let b = [end.x, end.y, end.z];
        if a.iter().chain(&b).any(|v| !v.is_finite()) {
            return Err(VolumeTransportError::Ray);
        }
        // Index-to-binary64 conversion must be exact before interval arithmetic.
        if self.grid.dims.iter().any(|&n| n as u128 > (1_u128 << 53)) {
            return Err(VolumeTransportError::Unrepresentable);
        }
        let mut result = ColumnDensityEnclosure {
            kg_per_m2: Interval::point(0.0), visited_nodes: 0,
            possible_cells: 0, uncertain_parallel_cells: 0,
        };
        if a == b { return Ok(result); }
        let segment = Segment::new(a, b)?;
        let mut stack = [IndexBox::EMPTY; STACK_CAPACITY];
        stack[0] = IndexBox { lo: [0; 3], hi: self.grid.dims };
        let mut pending = 1;
        let mut mass_per_parameter = Interval::point(0.0);
        while pending > 0 {
            checkpoint(cx)?;
            if result.visited_nodes == budget.max_nodes_per_segment {
                return Err(VolumeTransportError::Budget);
            }
            result.visited_nodes += 1;
            pending -= 1;
            let indices = stack[pending];
            let lower = std::array::from_fn(|axis| plane(
                self.grid.origin[axis], self.grid.cell[axis], indices.lo[axis],
            ));
            let upper = std::array::from_fn(|axis| plane(
                self.grid.origin[axis], self.grid.cell[axis], indices.hi[axis],
            ));
            if lower.iter().chain(&upper).any(|&v| !finite(v)) {
                return Err(VolumeTransportError::Unrepresentable);
            }
            let Some((fraction, parallel_uncertain)) = segment.overlap(lower, upper) else {
                continue;
            };
            let spans = std::array::from_fn::<_, 3, _>(|axis| indices.hi[axis] - indices.lo[axis]);
            // Largest index extent; ties retain x, then y, then z ordering.
            let mut axis = 0;
            for candidate in 1..3 {
                if spans[candidate] > spans[axis] { axis = candidate; }
            }
            if spans[axis] == 1 {
                if result.possible_cells == self.max_steps {
                    return Err(VolumeTransportError::Budget);
                }
                result.possible_cells += 1;
                result.uncertain_parallel_cells += u32::from(parallel_uncertain);
                let [x, y, z] = indices.lo;
                let index = (z * self.grid.dims[1] + y) * self.grid.dims[0] + x;
                let rho = self.grid.data[index];
                if rho != 0.0 {
                    let fraction_mass = nonnegative(fraction * Interval::point(rho));
                    let column = nonnegative(fraction_mass * segment.length);
                    if !finite(column) { return Err(VolumeTransportError::Unrepresentable); }
                    visit(index, column)?;
                    mass_per_parameter = nonnegative(mass_per_parameter + fraction_mass);
                    if !finite(mass_per_parameter) {
                        return Err(VolumeTransportError::Unrepresentable);
                    }
                }
            } else {
                if pending + 2 > stack.len() { return Err(VolumeTransportError::Budget); }
                let split = indices.lo[axis] + spans[axis] / 2;
                let mut left = indices;
                let mut right = indices;
                left.hi[axis] = split;
                right.lo[axis] = split;
                stack[pending] = right;
                stack[pending + 1] = left;
                pending += 2;
            }
        }
        result.kg_per_m2 = if mass_per_parameter.hi() == 0.0 {
            Interval::point(0.0)
        } else {
            nonnegative(mass_per_parameter * segment.length)
        };
        if !finite(result.kg_per_m2) { return Err(VolumeTransportError::Unrepresentable); }
        checkpoint(cx)?;
        Ok(result)
    }

    /// Bound direct, unscattered packet transmittance. LTE emission is not part
    /// of transmittance; neither it nor in-scattering is silently added here.
    /// Geometry visibility remains the responsibility of the calling transport.
    /// Coefficients are prepared within the default retained-value cap; repeated
    /// rays should use `prepare_extinction_enclosure` to amortize that work.
    pub fn transmittance_enclosure<const N: usize>(
        &self,
        cx: &Cx<'_>,
        start: Point3,
        end: Point3,
        wavelengths_nm: &[f64; N],
        budget: GridEnclosureBudget,
    ) -> Result<[Interval; N], VolumeTransportError> {
        checkpoint(cx)?;
        self.prepare_extinction_enclosure(cx, wavelengths_nm, MAX_PREPARED_EXTINCTION_VALUES)?
            .transmittance(cx, start, end, budget)
    }
}

/// Apply the same column integral to independently prepared spectral lanes.
pub(crate) fn attenuate<const N: usize>(
    column: Interval,
    coefficients: &[Interval; N],
) -> Result<[Interval; N], VolumeTransportError> {
    let mut result = [Interval::point(1.0); N];
    for (lane, coefficient) in coefficients.iter().enumerate() {
        if column.hi() == 0.0 || coefficient.hi() == 0.0 { continue; }
        let depth = nonnegative(column * *coefficient);
        result[lane] = exp_neg(depth).ok_or(VolumeTransportError::Unrepresentable)?;
    }
    Ok(result)
}

// At most sum(ceil(log2(dims))) <= usize::BITS + 2 splits on one
// root-to-leaf path, since the product of dimensions fits in usize.
const STACK_CAPACITY: usize = usize::BITS as usize + 4;

#[derive(Clone, Copy)]
struct IndexBox { lo: [usize; 3], hi: [usize; 3] }
impl IndexBox { const EMPTY: Self = Self { lo: [0; 3], hi: [0; 3] }; }

struct Segment {
    start: [f64; 3],
    end: [f64; 3],
    delta: [Interval; 3],
    length: Interval,
}

impl Segment {
    fn new(start: [f64; 3], end: [f64; 3]) -> Result<Self, VolumeTransportError> {
        let delta = std::array::from_fn(|axis| {
            if start[axis] == end[axis] { Interval::point(0.0) }
            else { Interval::point(end[axis]) - Interval::point(start[axis]) }
        });
        if delta.iter().any(|&v| !finite(v)) {
            return Err(VolumeTransportError::Unrepresentable);
        }
        for axis in 0..3 {
            if start[axis] != end[axis] && delta[axis].lo() <= 0.0 && delta[axis].hi() >= 0.0 {
                return Err(VolumeTransportError::Unrepresentable);
            }
        }
        let scale = delta.iter().map(|v| v.lo().abs().max(v.hi().abs())).fold(0.0_f64, f64::max);
        if scale == 0.0 { return Err(VolumeTransportError::Unrepresentable); }
        let scaled = delta.map(|v| v / Interval::point(scale));
        let squares = scaled.map(|v| nonnegative(v * v));
        let length = nonnegative(squares[0] + squares[1] + squares[2]).sqrt() * Interval::point(scale);
        if !finite(length) || length.lo() <= 0.0 {
            return Err(VolumeTransportError::Unrepresentable);
        }
        Ok(Self { start, end, delta, length })
    }

    fn overlap(&self, lower: [Interval; 3], upper: [Interval; 3]) -> Option<(Interval, bool)> {
        let mut enter = Interval::point(0.0);
        let mut exit = Interval::point(1.0);
        let mut uncertain_parallel = false;
        for axis in 0..3 {
            if self.start[axis] == self.end[axis] {
                let value = self.start[axis];
                if value < lower[axis].lo() || value >= upper[axis].hi() { return None; }
                uncertain_parallel |= value < lower[axis].hi() || value >= upper[axis].lo();
                continue;
            }
            let first = (lower[axis] - Interval::point(self.start[axis])) / self.delta[axis];
            let second = (upper[axis] - Interval::point(self.start[axis])) / self.delta[axis];
            let (near, far) = if self.start[axis] < self.end[axis] { (first, second) } else { (second, first) };
            enter = Interval::new(enter.lo().max(near.lo()), enter.hi().max(near.hi()));
            exit = Interval::new(exit.lo().min(far.lo()), exit.hi().min(far.hi()));
            if enter.lo() >= exit.hi() { return None; }
        }
        let width = nonnegative(exit - enter);
        let upper = width.hi().min(1.0);
        if upper == 0.0 { return None; }
        Some((Interval::new(if uncertain_parallel { 0.0 } else { width.lo() }, upper), uncertain_parallel))
    }
}

fn plane(origin: f64, cell: f64, index: usize) -> Interval {
    if index == 0 { Interval::point(origin) }
    else { Interval::point(origin) + Interval::point(cell) * Interval::point(index as f64) }
}
fn nonnegative(value: Interval) -> Interval { Interval::new(value.lo().max(0.0), value.hi().max(0.0)) }
fn finite(value: Interval) -> bool { value.lo().is_finite() && value.hi().is_finite() }

#[cfg(test)]
mod tests;
