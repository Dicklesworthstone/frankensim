//! Nodal lay clearance and spatial slit area on the retained plate shape.
//!
//! Contact uses area-lumped triangle quadrature of the existing power-law
//! potential. Slit flow integrates the positive part of the P1 gap EXACTLY
//! along each declared boundary edge, including partially closed edges. The
//! structural basis stays fixed; this is not a multimode or large-strain model.
use super::{AcousticRealizeError, PlateApertureReduction, invalid};
use fs_dcontact::Obstacle;
use crate::unilateral_contact::distributed::MAX_APERTURE_CONTACT_POINTS;

/// Complete explicit clearance/contact description on the original plate mesh.
#[derive(Debug, Clone)]
pub struct PlateClosureSpec {
    /// Signed rest clearances at EVERY mesh node [m]. Negative means installed
    /// interference, not an inferred correction. Used at the slit and the lay.
    pub nodal_rest_gap_m: Vec<f64>,
    /// Unique zero-based triangle indices covered by the compliant lay.
    pub lay_triangles: Vec<usize>,
    /// Normal pressure coefficient [Pa / m^alpha]. Triangle areas become
    /// quadrature weights [m²]; no dimensionless weight normalization occurs.
    pub stiffness_pa_per_m_alpha: f64,
    /// Power-law exponent, at least one.
    pub alpha: f64,
    /// Explicit Hunt--Crossley coefficient [s/m], nonnegative.
    pub internal_loss_s_per_m: f64,
    /// Original source or explicit authored-model label, never inferred.
    pub provenance: String,
    /// Maximum admitted local lay penetration [m], strictly positive.
    pub max_penetration_m: f64,
}

/// Physical geometric observation at one retained scalar opening.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateClosureProbe {
    /// Sum of positive-gap edge areas [m²], used by the Bernoulli junction.
    pub open_area_m2: f64,
    /// Lay quadrature nodes at positive penetration; not exact contact area.
    pub active_lay_points: usize,
    /// Sum of active nodal area weights [m²]; quadrature estimate only.
    pub active_lay_area_m2: f64,
    /// Maximum nodal lay penetration [m].
    pub max_penetration_m: f64,
}

/// Immutable mapping bound by the runtime constructor to its exact plate.
/// Source clearances and law remain separately inspectable. No observer changes
/// the model, and no gap is changed to force an authored closing time.
#[derive(Debug, Clone)]
pub struct PlateClosure {
    spec: PlateClosureSpec,
    law: Obstacle,
    // Coordinate y is the earlier reduction's mean-slit opening; a supplied
    // nonuniform rest profile may have a different physical mean clearance.
    rest_coordinate_m: f64,
    shape: Vec<f64>,
    edges: Vec<([usize;2],f64)>,
    nodes: Vec<usize>,
}
impl PlateApertureReduction {
    /// Compile the supplied lay against this plate's actual triangles and mode.
    /// This is a geometry adapter over fs-dcontact, not a new collision law.
    ///
    /// # Errors
    /// Incomplete/nonfinite profile, invalid/duplicate patch triangles,
    /// quadrature budget, contact-law admission or unrepresentable geometry.
    pub fn compile_closure(&self, spec: PlateClosureSpec) -> Result<PlateClosure, AcousticRealizeError> {
        let mesh = &self.chart().mesh;
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
        let shape: Vec<_> = self.shape_per_opening().iter().map(|v| v[0]).collect();
        let h = self.options().rest_opening_m;
        // Local gap = g_i + shape_i*(y-H). Obstacle penetration b_i*y-c_i
        // therefore uses b_i=-shape_i and c_i=g_i-shape_i*H.
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
    /// Nonfinite geometry or excessive local lay penetration.
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
fn positive_linear_mean(a: f64,b: f64) -> f64 {
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
