//! Facing-surface overlap of two independently deformed deposited roads.
//!
//! Merge the two world-x partitions, invert BOTH surface maps, and measure
//! each quadratic gap on their common interval. Reference areas on the two
//! bodies can differ; neither is replaced by a shared undeformed area. The
//! shared projected area is not curved area or a thermal/adhesion model.

use super::{Bridge, Error, finite, fraction_below, poll, positive};

/// Parallel local x/z frames translated into one common world frame, in metres.
/// Lower uses its outermost top face; upper uses its outermost bottom face.
#[derive(Debug, Clone, Copy)]
pub struct Query {
    pub lower_origin_m: [f64; 2],
    pub upper_origin_m: [f64; 2],
    /// Centred overlap width; must fit both selected faces. No lateral offsets.
    pub width_m: f64,
    /// Explicit geometric capture distance; not inferred contact or weld strength.
    pub capture_gap_m: f64,
    /// Maximum admitted overlap in metres; larger penetration refuses all rows.
    pub penetration_tolerance_m: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Budget {
    /// Total cells across both input roads, counted even for no overlap.
    pub max_cells: usize,
    /// All common-partition rows, including those with zero captured area.
    pub max_rows: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Row {
    pub lower_cell_id: u64,
    pub upper_cell_id: u64,
    pub world_interval_m: [f64; 2],
    /// Full overlap intervals in each road's own material-reference x coordinate.
    pub lower_reference_interval_m: [f64; 2],
    pub upper_reference_interval_m: [f64; 2],
    /// Extrema include the interior stationary point of the relative gap.
    pub minimum_gap_m: f64,
    pub maximum_gap_m: f64,
    /// Captured material areas can differ because the two surface stretches differ.
    pub lower_reference_area_m2: f64,
    pub upper_reference_area_m2: f64,
    /// Captured world-x projected area, not either reference or curved surface area.
    pub projected_area_m2: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Footprint {
    pub rows: Vec<Row>,
    /// Entire projected length where these two faces overlap, captured or not.
    pub overlap_projected_length_m: f64,
    /// None means these faces do not overlap; it says nothing about other bodies.
    pub minimum_gap_m: Option<f64>,
    pub cells_visited: usize,
}

struct Segment {
    id: u64,
    start: f64,
    length: f64,
    left: f64,
    right: f64,
    stretch: f64,
    height: f64,
    slope: f64,
    curvature: f64,
}
impl Segment {
    // Coefficients of height at x=a+(b-a)*u, with u in [0,1]. Keep operations
    // local rather than subtracting large squared world coordinates.
    fn polynomial(&self, a: f64, b: f64) -> Result<([f64; 3], f64, f64), Error> {
        let s = finite((a-self.left)/self.stretch)?;
        let ds = finite((b-a)/self.stretch)?;
        Ok(([finite(self.height+s*(self.slope+0.5*s*self.curvature))?,
            finite(ds*(self.slope+s*self.curvature))?,
            finite(0.5*self.curvature*ds*ds)?], s, ds))
    }
}

impl Bridge<'_> {
    /// Geometric capture between this lower road and an independently deformed
    /// upper road at the SAME accepted time. Meshes, lengths, thermal histories
    /// and translations may differ. No history is changed and no forces are solved.
    ///
    /// O(n_lower+n_upper+rows) work and storage, with two bounded surface lists.
    /// Each pair of cell intervals is visited once, not a quadratic all-pairs
    /// search. Cancellation/invalid input/penetration/budget failures return no
    /// partial footprint. Results are nominal floating-point geometry.
    ///
    /// A reference-station neighbour-force solve does not prove that these world
    /// surfaces are nonpenetrating between stations. This independent query can
    /// refuse such a state. Capture distance must be supplied separately from
    /// mechanics; captured area alone creates neither force nor a healed bond.
    pub fn neighbour_contact(&self, upper: &Bridge<'_>, query: Query, budget: Budget,
        cancelled: &mut impl FnMut() -> bool) -> Result<Footprint, Error> {
        poll(cancelled)?;
        if self.time_s() != upper.time_s() {
            return Err(Error::Invalid("neighbour surfaces require the same accepted time"));
        }
        let count = self.cells.len().checked_add(upper.cells.len()).ok_or(Error::Budget)?;
        if count > budget.max_cells { return Err(Error::Budget); }
        positive(query.width_m)?;
        for v in query.lower_origin_m.into_iter().chain(query.upper_origin_m) { finite(v)?; }
        if !query.capture_gap_m.is_finite() || query.capture_gap_m < 0.0
            || !query.penetration_tolerance_m.is_finite() || query.penetration_tolerance_m < 0.0 {
            return Err(Error::Invalid("invalid pair capture/penetration allowance"));
        }
        let lower = surface(self, true, query.lower_origin_m, query.width_m, cancelled)?;
        let upper = surface(upper, false, query.upper_origin_m, query.width_m, cancelled)?;
        let (mut i, mut j) = (0, 0);
        let mut rows = Vec::new();
        let mut covered = 0.0;
        let mut minimum: Option<f64> = None;
        while i < lower.len() && j < upper.len() {
            poll(cancelled)?;
            let l = &lower[i]; let u = &upper[j];
            let a = l.left.max(u.left); let b = l.right.min(u.right);
            if b > a {
                if rows.len() >= budget.max_rows { return Err(Error::Budget); }
                let (lp, ls, ld) = l.polynomial(a,b)?;
                let (up, us, ud) = u.polynomial(a,b)?;
                let c = [finite(up[0]-lp[0])?, finite(up[1]-lp[1])?, finite(up[2]-lp[2])?];
                let end = finite(c[0]+c[1]+c[2])?;
                let mut lo = c[0].min(end); let mut hi = c[0].max(end);
                if c[2] != 0.0 {
                    let vertex = -0.5*c[1]/c[2];
                    if vertex > 0.0 && vertex < 1.0 {
                        let value = finite(c[0]+vertex*(c[1]+vertex*c[2]))?;
                        lo = lo.min(value); hi = hi.max(value);
                    }
                }
                if lo < -query.penetration_tolerance_m {
                    return Err(Error::Invalid("deformed neighbour faces penetrate between stations"));
                }
                let f = fraction_below(c[2],c[1],finite(c[0]-query.capture_gap_m)?)?;
                rows.push(Row { lower_cell_id: l.id, upper_cell_id: u.id,
                    world_interval_m: [a,b],
                    lower_reference_interval_m: reference_interval(l,ls,ld)?,
                    upper_reference_interval_m: reference_interval(u,us,ud)?,
                    minimum_gap_m: lo, maximum_gap_m: hi,
                    lower_reference_area_m2: finite(query.width_m*(ld*f))?,
                    upper_reference_area_m2: finite(query.width_m*(ud*f))?,
                    projected_area_m2: finite(query.width_m*((b-a)*f))? });
                minimum = Some(minimum.map_or(lo, |old| old.min(lo)));
                covered = finite(covered+(b-a))?;
            }
            // Advance both at a shared boundary, avoiding duplicate zero-width rows.
            if l.right <= u.right { i += 1; }
            if u.right <= l.right { j += 1; }
        }
        poll(cancelled)?;
        Ok(Footprint { rows, overlap_projected_length_m: covered,
            minimum_gap_m: minimum, cells_visited: count })
    }
}

fn reference_interval(segment: &Segment, s: f64, ds: f64) -> Result<[f64; 2], Error> {
    // Exact endpoints use the retained reference partition, not re-derived sums.
    let a = if s == 0.0 { segment.start } else { finite(segment.start+s)? };
    let end = if s+ds == segment.length { finite(segment.start+segment.length)? }
        else { finite(segment.start+s+ds)? };
    Ok([a,end])
}

fn surface(beam: &Bridge<'_>, top: bool, origin: [f64; 2], width: f64,
    cancelled: &mut impl FnMut() -> bool) -> Result<Vec<Segment>, Error> {
    let level = |g: &crate::deposition::LayerGeometry| g.centroid_z_m
        + if top { 0.5*g.thickness_m } else { -0.5*g.thickness_m };
    let layer = beam.cells[0].section.layers().values().max_by(|a,b| {
        if top { level(&a.geometry).total_cmp(&level(&b.geometry)) }
        else { level(&b.geometry).total_cmp(&level(&a.geometry)) }
    }).ok_or(Error::Invalid("unborn road has no neighbour surface"))?;
    let z = finite(level(&layer.geometry))?;
    if width > layer.geometry.area_m2/layer.geometry.thickness_m {
        return Err(Error::Invalid("pair width exceeds a facing layer"));
    }
    let mut result = Vec::with_capacity(beam.cells.len());
    for (k, cell) in beam.cells.iter().enumerate() {
        poll(cancelled)?;
        let left = finite(origin[0]+cell.start_m+cell.u_m-z*cell.slope)?;
        let right = if let Some(next) = beam.cells.get(k+1) {
            finite(origin[0]+next.start_m+next.u_m-z*next.slope)?
        } else { finite(origin[0]+beam.deformed_position_m(beam.length_m,z)?[0])? };
        let stretch = finite((right-left)/cell.length_m())?;
        positive(stretch)?;
        result.push(Segment { id: cell.id, start: cell.start_m, length: cell.length_m(),
            left, right, stretch, height: finite(origin[1]+z+cell.w_m)?,
            slope: cell.slope, curvature: cell.section.curvature_per_m() });
    }
    Ok(result)
}

#[cfg(test)]
mod tests;
