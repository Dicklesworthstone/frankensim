//! Straight later-pass footprint against the accepted small-deformation surface.
//! This is a geometric capture model, not a contact-force or adhesion law.

use super::{Bridge, Error, finite, poll, positive};

/// Finite-nozzle sweeps against the frozen, actually deformed top surface.
pub mod nozzle;
/// Facing-surface contact footprints between two independently deformed roads.
pub mod pair;

/// A centered straight footprint in the same local WORLD x-z frame as the beam.
#[derive(Debug, Clone, Copy)]
pub struct Pass {
    pub start_x_m: f64,
    pub end_x_m: f64,
    pub start_lower_z_m: f64,
    pub end_lower_z_m: f64,
    /// Centered overlap width; cannot exceed the current top layer's width.
    pub width_m: f64,
    /// Explicit geometric capture distance, not inferred roughness or strength.
    pub capture_gap_m: f64,
    /// Allowed signed-overlap roundoff/model tolerance; larger penetration refuses.
    pub penetration_tolerance_m: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ContactRow {
    pub cell_id: u64,
    pub world_interval_m: [f64; 2],
    /// Extrema include an interior quadratic vertex, not only the endpoints.
    pub minimum_gap_m: f64,
    pub maximum_gap_m: f64,
    /// Original-material area, for reference-area thermal/contact consumers.
    pub reference_area_m2: f64,
    /// Current horizontal projected area; NOT curved surface area.
    pub projected_area_m2: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ContactFootprint {
    pub rows: Vec<ContactRow>,
    /// Footprint length with no surface from THIS bridge. Other bodies are unknown.
    pub uncovered_projected_length_m: f64,
}

impl Bridge<'_> {
    /// Intersect the actual deformed top surface with a later straight pass.
    /// Its x map is inverted cellwise, including axial contraction and the
    /// top fiber's rotation. The gap is quadratic; both extrema and capture
    /// length are solved on that polynomial rather than sampled at cell centres.
    ///
    /// Areas are distinct reference/current-projected quantities. Conversion to
    /// conductance or strength needs an explicit independent law. Excessive
    /// penetration returns no footprint, not fabricated support or impact force.
    /// Fixed frames, rectangular centered widths and single-valued top surface
    /// are required. Results are floating-point geometry, not interval certificates.
    pub fn next_pass_contact(&self, pass: Pass, max_rows: usize,
        cancelled: &mut impl FnMut() -> bool) -> Result<ContactFootprint, Error> {
        poll(cancelled)?;
        finite(pass.start_x_m)?; finite(pass.end_x_m)?;
        finite(pass.start_lower_z_m)?; finite(pass.end_lower_z_m)?;
        let length = finite(pass.end_x_m-pass.start_x_m)?;
        positive(length)?; positive(pass.width_m)?;
        if !pass.capture_gap_m.is_finite() || pass.capture_gap_m < 0.0
            || !pass.penetration_tolerance_m.is_finite() || pass.penetration_tolerance_m < 0.0 {
            return Err(Error::Invalid("invalid capture or penetration tolerance"));
        }
        let top = self.cells[0].section.layers().values().max_by(|a,b| {
            (a.geometry.centroid_z_m+0.5*a.geometry.thickness_m)
                .total_cmp(&(b.geometry.centroid_z_m+0.5*b.geometry.thickness_m))
        }).ok_or(Error::Invalid("an unborn span has no contact surface"))?;
        let z = top.geometry.centroid_z_m+0.5*top.geometry.thickness_m;
        if pass.width_m > top.geometry.area_m2/top.geometry.thickness_m {
            return Err(Error::Invalid("pass width exceeds the centered top-face overlap"));
        }
        let path_slope = finite((pass.end_lower_z_m-pass.start_lower_z_m)/length)?;
        let mut rows = Vec::new();
        let mut covered = 0.0;
        for (index, cell) in self.cells.iter().enumerate() {
            poll(cancelled)?;
            let h = cell.length_m();
            let left = finite(cell.start_m+cell.u_m-z*cell.slope)?;
            // Shared endpoints keep the inverted world partition continuous.
            let right = if let Some(next) = self.cells.get(index+1) {
                finite(next.start_m+next.u_m-z*next.slope)?
            } else {
                self.deformed_position_m(self.length_m,z)?[0]
            };
            let stretch = finite((right-left)/h)?;
            positive(stretch)?;
            let a = left.max(pass.start_x_m);
            let b = right.min(pass.end_x_m);
            if b <= a { continue; }
            if rows.len() >= max_rows { return Err(Error::Budget); }
            let s0 = (a-left)/stretch;
            let ds = (b-a)/stretch;
            let c0 = finite(pass.start_lower_z_m+path_slope*(a-pass.start_x_m)-z-cell.height(s0))?;
            let c1 = finite(path_slope*(b-a)-ds*(cell.slope+cell.section.curvature_per_m()*s0))?;
            let c2 = finite(-0.5*cell.section.curvature_per_m()*ds*ds)?;
            let gap = |u: f64| c0+u*(c1+u*c2);
            let end = finite(gap(1.0))?;
            let mut minimum = c0.min(end);
            let mut maximum = c0.max(end);
            if c2 != 0.0 {
                let critical = -0.5*c1/c2;
                if critical > 0.0 && critical < 1.0 {
                    let value = finite(gap(critical))?;
                    minimum = minimum.min(value); maximum = maximum.max(value);
                }
            }
            if minimum < -pass.penetration_tolerance_m {
                return Err(Error::Invalid("later pass penetrates the accepted bridge surface"));
            }
            let fraction = fraction_below(c2,c1,finite(c0-pass.capture_gap_m)?)?;
            let reference_area = finite(pass.width_m*(ds*fraction))?;
            let projected_area = finite(pass.width_m*((b-a)*fraction))?;
            rows.push(ContactRow { cell_id: cell.id, world_interval_m: [a,b],
                minimum_gap_m: minimum, maximum_gap_m: maximum,
                reference_area_m2: reference_area, projected_area_m2: projected_area });
            covered = finite(covered+(b-a))?;
        }
        poll(cancelled)?;
        Ok(ContactFootprint { rows, uncovered_projected_length_m: (length-covered).max(0.0) })
    }
}

// Measure of {u in [0,1]: a*u^2+b*u+c <= 0}. Normalization and the conjugate
// quadratic root avoid scale overflow and catastrophic cancellation of one root.
fn fraction_below(a: f64, b: f64, c: f64) -> Result<f64, Error> {
    let scale = a.abs().max(b.abs()).max(c.abs());
    if scale == 0.0 { return Ok(1.0); }
    positive(scale)?;
    let (a,b,c) = (a/scale,b/scale,c/scale);
    let mut cuts = vec![0.0,1.0];
    let mut add = |root: f64| { if root > 0.0 && root < 1.0 { cuts.push(root); } };
    if a == 0.0 {
        if b != 0.0 { add(-c/b); }
    } else {
        let discriminant = b*b-4.0*a*c;
        if discriminant >= 0.0 {
            let root = fs_math::det::sqrt(discriminant);
            let q = -0.5*(b+root.copysign(b));
            if q == 0.0 { add(-0.5*b/a); }
            else { add(q/a); add(c/q); }
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    let mut fraction = 0.0;
    for pair in cuts.windows(2) {
        let middle = 0.5*(pair[0]+pair[1]);
        if c+middle*(b+middle*a) <= 0.0 { fraction += pair[1]-pair[0]; }
    }
    finite(fraction)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bridge::tests::{birth, material, temperatures, config};
    use crate::bridge::Support;
    fn near(a: f64,b: f64) { assert!((a-b).abs()<1e-11,"{a:e} != {b:e}"); }
    fn pass(start: f64,end: f64,z: f64,gap: f64) -> Pass {
        Pass { start_x_m: start,end_x_m: end,start_lower_z_m: z,end_lower_z_m: z,
            width_m: 0.001,capture_gap_m: gap,penetration_tolerance_m: 1e-12 }
    }

    #[test]
    fn quadratic_capture_resolves_interior_gaps_and_is_scale_invariant() {
        for scale in [1e-200,1.0,1e200] {
            near(fraction_below(-scale,scale,-scale/16.0).unwrap(),1.0-(0.75_f64).sqrt());
            near(fraction_below(scale,-scale,3.0*scale/16.0).unwrap(),0.5);
            near(fraction_below(0.0,scale,-scale/4.0).unwrap(),0.25);
        }
        near(fraction_below(0.0,0.0,0.0).unwrap(),1.0);
        near(fraction_below(0.0,0.0,1.0).unwrap(),0.0);
        near(fraction_below(1.0,-1.0,0.25).unwrap(),0.0);
    }

    #[test]
    fn solved_sag_changes_contact_despite_identical_endpoint_heights() {
        let law = material(1e8,vec![]);
        let state = birth(1,Support::SimplySupported,&law);
        let (sagged,report) = state.advance_controlled(1.0,&temperatures(&state,300.0),9.81,
            &[0.0],config(),&mut || false).unwrap();
        let left = sagged.deformed_position_m(0.0,0.0005).unwrap()[0];
        let right = sagged.deformed_position_m(sagged.length_m(),0.0005).unwrap()[0];
        let query = pass(left,right,0.0005,report.maximum_sag_m/4.0);
        let contact = sagged.next_pass_contact(query,1,&mut || false).unwrap();
        near(contact.rows[0].reference_area_m2/(0.001*0.02),1.0-(0.75_f64).sqrt());
        near(contact.rows[0].maximum_gap_m,report.maximum_sag_m);
        assert!(contact.rows[0].projected_area_m2 < contact.rows[0].reference_area_m2);
        let flat = state.next_pass_contact(query,1,&mut || false).unwrap();
        assert!(flat.rows[0].reference_area_m2 > 5.0*contact.rows[0].reference_area_m2);
    }

    #[test]
    fn world_footprints_include_contraction_and_report_uncovered_length() {
        let law = material(1e8,vec![]);
        let state = birth(4,Support::SimplySupported,&law);
        let (expanded,_) = state.advance_controlled(1.0,&temperatures(&state,500.0),0.0,
            &[0.0;4],config(),&mut || false).unwrap();
        let end = expanded.deformed_position_m(expanded.length_m(),0.0005).unwrap()[0];
        let contact = expanded.next_pass_contact(pass(0.0,end+0.001,0.0005,0.0),4,&mut || false).unwrap();
        near(contact.uncovered_projected_length_m,0.001);
        near(contact.rows.iter().map(|r| r.reference_area_m2).sum(),0.001*0.02);
        near(contact.rows.iter().map(|r| r.projected_area_m2).sum(),0.001*end);
    }

    #[test]
    fn penetration_budget_and_cancellation_return_no_partial_footprint() {
        let law = material(1e8,vec![]);
        let state = birth(4,Support::SimplySupported,&law);
        let query = pass(0.0,0.02,0.0005,1e-6);
        assert!(state.next_pass_contact(query,3,&mut || false).is_err());
        assert!(state.next_pass_contact(Pass { start_lower_z_m: 0.0,..query },4,&mut || false).is_err());
        let mut calls = 0;
        assert!(state.next_pass_contact(query,4,&mut || { calls+=1; calls>3 }).is_err());
        assert_eq!(state.time_s(),0.0);
    }
}
