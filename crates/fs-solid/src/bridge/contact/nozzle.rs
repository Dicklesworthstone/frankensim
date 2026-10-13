//! Continuous clearance of a finite nozzle over a frozen deformed top surface.
//!
//! The supplied axisymmetric lower envelope has a flat tip and conical flank.
//! Its axis stays on the beam's transverse centreline, so the smallest clearance
//! occurs in the x/z section. Motion is affine in x/z during one sweep. This
//! checks the represented TOP surface only: it is not a whole-scene collision,
//! side-wall/overhang, rotating/laterally offset nozzle, or impact-force model.

use super::{Bridge, Error, finite, poll, positive};

/// Explicit lower envelope, not an inferred printer geometry. The material
/// above this envelope within `outer_radius_m` is treated as occupied.
#[derive(Debug, Clone, Copy)]
pub struct Profile {
    pub flat_radius_m: f64,
    pub outer_radius_m: f64,
    /// Rise in underside height per radial metre beyond the flat tip; >= 0.
    pub flank_slope: f64,
}
#[derive(Debug, Clone, Copy)]
pub struct Sweep {
    /// [axis world x, flat-tip underside world z], metres.
    pub start_m: [f64; 2],
    pub end_m: [f64; 2],
    pub profile: Profile,
    /// A nominal gap <= this nonnegative margin is reported as interference.
    pub required_clearance_m: f64,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Witness {
    pub cell_id: u64,
    pub reference_x_m: f64,
    pub world_x_m: f64,
    /// Interpolates the supplied motion, not the exact first time of impact.
    pub sweep_fraction: f64,
    /// Underside minus top-surface height, m; negative means penetration.
    pub gap_m: f64,
}
#[derive(Debug, Clone, PartialEq)]
pub struct Report {
    /// None means no swept overlap with THIS top surface, not scene clearance.
    pub minimum: Option<Witness>,
    pub interferes: bool,
    pub checked_cells: usize,
}

impl Bridge<'_> {
    /// Minimize clearance over the entire sweep and finite footprint. Neither
    /// endpoint-only motion checks nor centreline-only tool checks are used.
    /// Returns penetration as an observation rather than refusing its geometry.
    ///
    /// On each nozzle flank and beam cell the feasible (position,time) set is
    /// a clipped square and the gap is quadratic in position, affine in time.
    /// Its minimum is attained on the polygon boundary (also for a time-flat
    /// objective); endpoint and stationary edge values exhaust that boundary.
    /// This is nominal floating-point geometry, NOT an interval certificate.
    /// Axial strain and fibre rotation enter the same world-x inversion as the
    /// deposition footprint. Freeze mechanical history for this sweep; subdivide
    /// the physical timeline when that approximation is not adequate.
    ///
    /// `max_cells` bounds the whole traversal. Cancellation or invalid arithmetic
    /// returns no partial clearance claim. O(cells) time, O(1) working storage.
    pub fn nozzle_clearance(&self, sweep: Sweep, max_cells: usize,
        cancelled: &mut impl FnMut() -> bool) -> Result<Report, Error> {
        poll(cancelled)?;
        if self.cells.len() > max_cells { return Err(Error::Budget); }
        for v in sweep.start_m.into_iter().chain(sweep.end_m) { finite(v)?; }
        let p = sweep.profile;
        positive(p.outer_radius_m)?;
        if !p.flat_radius_m.is_finite() || p.flat_radius_m < 0.0 || p.flat_radius_m > p.outer_radius_m
            || !p.flank_slope.is_finite() || p.flank_slope < 0.0
            || !sweep.required_clearance_m.is_finite() || sweep.required_clearance_m < 0.0 {
            return Err(Error::Invalid("invalid nozzle envelope or clearance margin"));
        }
        let dx = finite(sweep.end_m[0]-sweep.start_m[0])?;
        let dz = finite(sweep.end_m[1]-sweep.start_m[1])?;
        let top = self.cells[0].section.layers().values().map(|l|
            l.geometry.centroid_z_m+0.5*l.geometry.thickness_m)
            .max_by(f64::total_cmp).ok_or(Error::Invalid("an unborn span has no nozzle-check surface"))?;
        let mut minimum: Option<Witness> = None;
        for (index, cell) in self.cells.iter().enumerate() {
            poll(cancelled)?;
            let h = cell.length_m();
            let left = finite(cell.start_m+cell.u_m-top*cell.slope)?;
            let right = if let Some(next) = self.cells.get(index+1) {
                finite(next.start_m+next.u_m-top*next.slope)?
            } else { self.deformed_position_m(self.length_m,top)?[0] };
            let width = finite(right-left)?; positive(width)?;
            let relative = finite(left-sweep.start_m[0])?;
            for (lo,hi,slope) in [(-p.outer_radius_m,-p.flat_radius_m,-p.flank_slope),
                (-p.flat_radius_m,p.flat_radius_m,0.0),(p.flat_radius_m,p.outer_radius_m,p.flank_slope)] {
                // u is normalized reference-cell position; t is motion fraction.
                let polygon = vec![[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0]];
                let polygon = clip(polygon, [width,-dx,finite(relative-hi)?])?;
                let polygon = clip(polygon, [-width,dx,finite(lo-relative)?])?;
                if polygon.is_empty() { continue; }
                let rise = if slope == 0.0 { 0.0 } else { -p.flank_slope*p.flat_radius_m };
                let c = [finite(sweep.start_m[1]+rise+slope*relative-top-cell.w_m)?,
                    finite(slope*width-h*cell.slope)?,
                    finite(-0.5*cell.section.curvature_per_m()*h*h)?, finite(dz-slope*dx)?];
                for i in 0..polygon.len() {
                    let a = polygon[i]; let b = polygon[(i+1)%polygon.len()];
                    let du = b[0]-a[0]; let dt = b[1]-a[1];
                    let quadratic = finite(c[2]*du*du)?;
                    let linear = finite((c[1]+2.0*c[2]*a[0])*du+c[3]*dt)?;
                    let interior = if quadratic > 0.0 { -0.5*linear/quadratic } else { -1.0 };
                    for v in [0.0,1.0,interior] {
                        if !(0.0..=1.0).contains(&v) { continue; }
                        let u = a[0]+v*du; let t = a[1]+v*dt;
                        let gap = finite(c[0]+u*(c[1]+u*c[2])+c[3]*t)?;
                        let witness = Witness { cell_id: cell.id,
                            reference_x_m: finite(cell.start_m+h*u)?, world_x_m: finite(left+width*u)?,
                            sweep_fraction: t, gap_m: gap };
                        if minimum.is_none_or(|old| gap < old.gap_m || (gap == old.gap_m
                            && (t,cell.id,witness.reference_x_m) < (old.sweep_fraction,old.cell_id,old.reference_x_m))) {
                            minimum = Some(witness);
                        }
                    }
                }
            }
        }
        poll(cancelled)?;
        Ok(Report { interferes: minimum.is_some_and(|w| w.gap_m <= sweep.required_clearance_m),
            minimum, checked_cells: self.cells.len() })
    }
}

// Clip a bounded polygon by a*u+b*t+c <= 0. Normalize before crossing
// interpolation so finite but large input coefficients do not overflow f0-f1.
fn clip(polygon: Vec<[f64;2]>, plane: [f64;3]) -> Result<Vec<[f64;2]>,Error> {
    let scale = plane.iter().fold(0.0_f64, |a,b| a.max(b.abs()));
    if scale == 0.0 { return Ok(polygon); }
    positive(scale)?;
    let [a,b,c] = plane.map(|v| v/scale);
    let side = |p: [f64;2]| a*p[0]+b*p[1]+c;
    let mut out = Vec::with_capacity(8);
    for i in 0..polygon.len() {
        let p = polygon[i]; let q = polygon[(i+1)%polygon.len()];
        let fp = side(p); let fq = side(q);
        if fp <= 0.0 { out.push(p); }
        if (fp <= 0.0) != (fq <= 0.0) {
            let t = finite(fp/(fp-fq))?;
            out.push([p[0]+t*(q[0]-p[0]),p[1]+t*(q[1]-p[1])]);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests;
