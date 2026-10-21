//! Nodal clearance and spatial slit area on the retained plate shape.
//!
//! Contact uses area-lumped triangle quadrature of the existing power-law
//! potential. Slit flow integrates the positive part of the P1 gap EXACTLY
//! along each declared boundary edge, including partially closed edges. The
//! structural basis stays fixed; this is not a multimode or large-strain model.
//! An opposed pair sums both blades' gap motion but has ONE contact interface:
//! neither its contact quadrature areas nor its slit width are doubled.
/// Matching-interface geometry for two independently moving blade modes.
pub mod paired;

use super::{AcousticRealizeError, PlateApertureReduction, invalid};
use fs_dcontact::Obstacle;
use crate::unilateral_contact::distributed::MAX_APERTURE_CONTACT_POINTS;
use crate::bernoulli_aperture::viscous::{ViscousProfileFlow, ViscousSlit};

/// Complete explicit clearance/contact description on the original plate mesh.
#[derive(Debug, Clone)]
pub struct PlateClosureSpec {
    /// Signed TOTAL rest clearances at EVERY mesh node [m]. Negative means
    /// installed interference. For a symmetric pair these are blade-to-blade
    /// gaps, not each blade's distance from the symmetry plane.
    pub nodal_rest_gap_m: Vec<f64>,
    /// Unique zero-based triangles of the common compliant contact interface.
    /// A stationary lay for one blade; blade-to-blade contact for an opposed pair.
    pub lay_triangles: Vec<usize>,
    /// Normal pressure coefficient [Pa / m^alpha]. Triangle areas become
    /// quadrature weights [m²]; no dimensionless weight normalization occurs.
    pub stiffness_pa_per_m_alpha: f64,
    /// Power-law exponent, at least one.
    pub alpha: f64,
    /// Explicit Hunt--Crossley coefficient [s/m], nonnegative.
    pub internal_loss_s_per_m: f64,
    /// Original source or explicit authored-model label for THIS contact pair.
    /// A blade/lay law is not automatically a blade/blade law.
    pub provenance: String,
    /// Maximum admitted local contact penetration [m], strictly positive.
    pub max_penetration_m: f64,
}

/// Physical geometric observation at one retained scalar opening.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateClosureProbe {
    /// Sum of positive-gap edge areas [m²], used by the Bernoulli junction.
    pub open_area_m2: f64,
    /// Contact quadrature nodes at positive penetration; not exact contact area.
    pub active_lay_points: usize,
    /// Sum of active nodal area weights [m²]; quadrature estimate only.
    pub active_lay_area_m2: f64,
    /// Maximum nodal contact penetration [m].
    pub max_penetration_m: f64,
}

/// Immutable mapping bound by the runtime constructor to its exact plate.
/// Source clearances and law remain separately inspectable. No observer changes
/// the model, and no gap is changed to force an authored closing time.
#[derive(Debug, Clone)]
pub struct PlateClosure {
    spec: PlateClosureSpec,
    law: Obstacle,
    // Coordinate y is the reduction's mean TOTAL opening; a supplied
    // nonuniform rest profile may have a different physical mean clearance.
    rest_coordinate_m: f64,
    shape: Vec<f64>,
    edges: Vec<([usize;2],f64)>,
    nodes: Vec<usize>,
}
impl PlateApertureReduction {
    /// Compile the supplied interface against the actual triangles and gap mode.
    /// This is a geometry adapter over fs-dcontact, not a new collision law.
    ///
    /// # Errors
    /// Incomplete/nonfinite profile, invalid/duplicate patch triangles,
    /// quadrature budget, contact-law admission or unrepresentable geometry.
    pub fn compile_closure(&self, spec: PlateClosureSpec) -> Result<PlateClosure, AcousticRealizeError> {
        let mesh = &self.chart().mesh;
        let (weights, nodes) = contact_quadrature(mesh, &spec)?;
        let blades = self.moving_blade_count() as f64;
        let shape: Vec<_> = self.shape_per_opening().iter().map(|v| blades * v[0]).collect();
        let h = self.options().rest_opening_m;
        // Local TOTAL gap = g_i + shape_i*(y-H). Obstacle penetration b_i*y-c_i
        // therefore uses b_i=-shape_i and c_i=g_i-shape_i*H. Each physical
        // interface point contributes once, including for two moving blades.
        let rows = nodes.iter().map(|&i| -shape[i]).collect();
        let gaps = nodes.iter().map(|&i| spec.nodal_rest_gap_m[i]-shape[i]*h).collect();
        let areas = nodes.iter().map(|&i| weights[i]).collect();
        let law = Obstacle::new(rows,nodes.len(),1,gaps,areas,
            spec.stiffness_pa_per_m_alpha,spec.alpha,spec.provenance.clone())
            .and_then(|law| law.with_internal_loss(spec.internal_loss_s_per_m))
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let edges = self.options().slit_edges.iter().map(|&edge| {
            let [a,b] = edge.map(|i| mesh.nodes[i]);
            (edge,(a.0-b.0).hypot(a.1-b.1))
        }).collect();
        Ok(PlateClosure { spec,law,rest_coordinate_m:h,shape,edges,nodes })
    }
}

// One physical area quadrature for scalar, paired and multimode interfaces.
// Modal adapters supply their original displacement columns separately.
pub(crate) fn contact_quadrature(mesh: &fs_plate::PlateMesh, spec: &PlateClosureSpec)
    -> Result<(Vec<f64>, Vec<usize>), AcousticRealizeError>
{
        if spec.nodal_rest_gap_m.len() != mesh.node_count()
            || spec.nodal_rest_gap_m.iter().any(|x| !x.is_finite())
            || spec.lay_triangles.is_empty() || spec.lay_triangles.len() > mesh.tris.len()
            || !spec.max_penetration_m.is_finite() || spec.max_penetration_m <= 0.0
            || spec.provenance.trim().is_empty() {
            return Err(invalid("plate closure requires complete finite nodal clearances, a lay patch and explicit penetration allowance"));
        }
        let mut weights = vec![0.0; mesh.node_count()];
        let mut seen = std::collections::BTreeSet::new();
        for &i in &spec.lay_triangles {
            let tri = mesh.tris.get(i).ok_or_else(|| invalid("lay triangle is outside the plate mesh"))?;
            if !seen.insert(i) { return Err(invalid("lay patch repeats a triangle")); }
            let [a,b,c] = tri.map(|n| mesh.nodes[n]);
            let area = 0.5*((b.0-a.0)*(c.1-a.1)-(b.1-a.1)*(c.0-a.0));
            if !area.is_finite() || area <= 0.0 { return Err(invalid("lay triangle has no finite positive area")); }
            for &node in tri { weights[node] += area/3.0; }
        }
        let nodes: Vec<_> = weights.iter().enumerate().filter(|(_,w)| **w > 0.0).map(|(i,_)| i).collect();
        if nodes.is_empty() || nodes.len() > MAX_APERTURE_CONTACT_POINTS
            || weights.iter().any(|w| !w.is_finite()) {
            return Err(invalid("plate lay exceeds its 4096-point contact budget or area overflowed"));
        }
        Ok((weights, nodes))
}
impl PlateClosure {
    /// Original unmodified clearance/contact input, with no normalized weights.
    #[must_use]
    pub const fn spec(&self) -> &PlateClosureSpec { &self.spec }
    /// Compiled scalar contact on original mesh coordinates and physical areas.
    #[must_use]
    pub const fn contact_law(&self) -> &Obstacle { &self.law }
    /// Mesh node for each row in the compiled contact obstacle.
    #[must_use]
    pub fn contact_nodes(&self) -> &[usize] { &self.nodes }
    /// Physical gap at a supplied mesh node. Negative means penetration.
    ///
    /// # Errors
    /// Absent node, nonfinite coordinate or arithmetic overflow.
    pub fn nodal_gap_m(&self, node: usize, opening_m: f64) -> Result<f64, AcousticRealizeError> {
        let s = self.shape.get(node).ok_or_else(|| invalid("closure node is outside the plate mesh"))?;
        let gap = self.spec.nodal_rest_gap_m[node] + s*(opening_m-self.rest_coordinate_m);
        if !gap.is_finite() { return Err(invalid("plate closure gap overflowed")); }
        Ok(gap)
    }
    /// Actual slit area at one midpoint. Integrate max(gap(s),0), NOT the
    /// positive part of a mean gap and NOT the trapezoid of clipped endpoints.
    ///
    /// # Errors
    /// Nonfinite coordinate or arithmetic overflow. This geometric trial query
    /// does not enforce final-state slope/penetration limits during root search.
    pub fn open_area_m2(&self, opening_m: f64) -> Result<f64, AcousticRealizeError> {
        if !opening_m.is_finite() { return Err(invalid("plate slit needs a finite coordinate")); }
        let mut area = 0.0;
        for &(edge,length) in &self.edges {
            let a = self.nodal_gap_m(edge[0],opening_m)?;
            let b = self.nodal_gap_m(edge[1],opening_m)?;
            area += length*positive_linear_mean(a,b);
        }
        if !area.is_finite() || area < 0.0 { return Err(invalid("plate slit area overflowed")); }
        Ok(area)
    }

    /// Integrate the channel law over the actual positive P1 gaps, edge by edge.
    /// A material change can alter the retained blade shape and therefore this
    /// flow even when its mean aperture is unchanged. No effective-height
    /// substitution or doubling for a symmetric opposed pair is performed.
    pub fn viscous_flow(&self, law: ViscousSlit, opening_m: f64,
        density_kg_m3: f64, pressure_drop_pa: f64)
        -> Result<ViscousProfileFlow, AcousticRealizeError>
    {
        let mut total = ViscousProfileFlow::default();
        for &(edge, length) in &self.edges {
            let gaps = [self.nodal_gap_m(edge[0], opening_m)?, self.nodal_gap_m(edge[1], opening_m)?];
            let local = law.linear_gap_flow(length, gaps, density_kg_m3, pressure_drop_pa)?;
            total.flow_m3_s += local.flow_m3_s;
            total.flow_per_pressure_m3_s_pa += local.flow_per_pressure_m3_s_pa;
        }
        if !total.flow_m3_s.is_finite() || !total.flow_per_pressure_m3_s_pa.is_finite() {
            return Err(invalid("spatial viscous slit flow overflowed"));
        }
        Ok(total)
    }
    /// Inspect area, nodal activity and penetration without moving physical state.
    ///
    /// # Errors
    /// Nonfinite coordinate or arithmetic overflow.
    pub fn probe(&self, opening_m: f64) -> Result<PlateClosureProbe, AcousticRealizeError> {
        let mut out = PlateClosureProbe { open_area_m2:self.open_area_m2(opening_m)?,
            active_lay_points:0,active_lay_area_m2:0.0,max_penetration_m:0.0 };
        for (i,&node) in self.nodes.iter().enumerate() {
            let gap = self.nodal_gap_m(node,opening_m)?;
            if gap < 0.0 {
                out.active_lay_points += 1;
                out.active_lay_area_m2 += self.law.weights()[i];
                out.max_penetration_m = out.max_penetration_m.max(-gap);
            }
        }
        if !out.active_lay_area_m2.is_finite() { return Err(invalid("active lay area overflowed")); }
        Ok(out)
    }
    /// Enforce the supplied penetration budget on an initial/candidate state.
    ///
    /// # Errors
    /// Nonfinite geometry or excessive local contact penetration.
    pub fn validate_opening(&self, opening_m: f64) -> Result<(), AcousticRealizeError> {
        if self.probe(opening_m)?.max_penetration_m > self.spec.max_penetration_m {
            return Err(invalid("plate closure exceeds its declared penetration allowance"));
        }
        Ok(())
    }
}

// Positive-part integral on a unit edge. Scaling avoids overflow in a-b for
// opposite large finite endpoints. At a zero crossing the open part is a
// triangle, not a trapezoid formed by clipping both endpoint samples first.
pub(crate) fn positive_linear_mean(a: f64,b: f64) -> f64 {
    if a >= 0.0 && b >= 0.0 { return 0.5*a+0.5*b; }
    if a <= 0.0 && b <= 0.0 { return 0.0; }
    let p = a.max(b);
    let n = -a.min(b);
    let scale = p.max(n);
    0.5*p*((p/scale)/(p/scale+n/scale))
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn partly_closed_edge_keeps_its_true_positive_area() {
        assert_eq!(positive_linear_mean(2.0,-2.0),0.5);
        assert_eq!(positive_linear_mean(-2.0,2.0),0.5);
        assert_eq!(positive_linear_mean(2.0,0.0),1.0);
        assert_eq!(positive_linear_mean(-2.0,0.0),0.0);
        assert_eq!(positive_linear_mean(2.0,4.0),3.0);
        let huge = positive_linear_mean(f64::MAX,-f64::MAX);
        assert!(huge.is_finite()); assert_eq!(huge,0.25*f64::MAX);
    }
    #[test]
    fn edge_subdivision_preserves_the_piecewise_linear_open_area() {
        for (a,b) in [(3.0,-1.0),(-3.0,1.0),(1.0,3.0),(-1.0,-3.0)] {
            let midpoint = f64::midpoint(a,b);
            let split = 0.5*positive_linear_mean(a,midpoint)+0.5*positive_linear_mean(midpoint,b);
            assert!((split-positive_linear_mean(a,b)).abs() <= 1e-14);
        }
    }
}
