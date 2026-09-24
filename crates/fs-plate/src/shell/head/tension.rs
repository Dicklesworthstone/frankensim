//! Equilibrated affine membrane prestress on the existing P1 head geometry.
//!
//! This is installed force per length, not an eigenfrequency adjustment. With
//! base tension T, the full tensor is
//! Nxx=T+xx+c*x+d*y, Nyy=T+yy+e*x+f*y, Nxy=xy-f*x-c*y.
//! The two in-plane divergence equations vanish identically. Rim tractions are
//! N*n; no unbalanced interior body force or fictitious tuning-lug force is
//! inserted. This family is the stress from a cubic Airy potential, not a full
//! model of hoops, lugs, bearing-edge friction or a measured tension map.
use crate::{PlateError, PlateMesh, PlateModel};
use fs_sparse::Coo;

/// Variation about a disk's explicitly supplied uniform installed tension.
/// Coordinates are the original head's Cartesian metres, centered on the disk.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TensionVariation {
    /// Constant additions [N/m] in the order xx, yy, xy (tensor shear).
    pub constant_n_m: [f64; 3],
    /// Spatial derivatives [N/m²] in the order c, d, e, f in this module's law.
    pub gradient_n_m2: [f64; 4],
}
fn bad(what: &'static str) -> PlateError { PlateError::BadSection { what } }

impl TensionVariation {
    /// Admit finite coefficients; tensile admissibility additionally depends
    /// on the original disk's base tension and actual meshed domain.
    /// # Errors
    /// Nonfinite coefficients.
    pub fn validate(self) -> Result<(), PlateError> {
        if self.constant_n_m.iter().chain(&self.gradient_n_m2).any(|v| !v.is_finite()) {
            return Err(bad("head tension variation needs finite SI coefficients"));
        }
        Ok(())
    }
    /// The original uniform image, including signed-zero input.
    #[must_use]
    pub fn is_zero(self) -> bool {
        self.constant_n_m.iter().chain(&self.gradient_n_m2).all(|v| *v == 0.0)
    }
    /// Full installed resultant [N/m] in xx, yy, xy order, not stress [Pa].
    /// Nonfinite results are rejected by preparation, never clipped.
    #[must_use]
    pub fn resultant(self, base_n_m: f64, point_m: [f64; 2]) -> [f64; 3] {
        let [x,y]=point_m; let [xx,yy,xy]=self.constant_n_m;
        let [c,d,e,f]=self.gradient_n_m2;
        [base_n_m+xx+c*x+d*y, base_n_m+yy+e*x+f*y, xy-f*x-c*y]
    }
    pub(super) fn validate_mesh(self, base: f64, mesh: &PlateMesh) -> Result<(), PlateError> {
        self.validate()?;
        if !base.is_finite() || base<=0.0 { return Err(bad("head base tension must be finite and positive")); }
        // An affine symmetric tensor at an interior point is a convex
        // combination of the vertex tensors. Positive definiteness at EVERY
        // vertex therefore admits the ENTIRE polygonal mesh, not just samples
        // at the integration points. No claim is made outside that mesh.
        for &(x,y) in &mesh.nodes {
            let [xx,yy,xy]=self.resultant(base,[x,y]);
            if [xx,yy,xy].iter().any(|v| !v.is_finite()) || xx<=0.0 || yy<=0.0 {
                return Err(bad("installed head tensor must remain finite and tensile across the mesh"));
            }
            let scale=xx.max(yy).max(xy.abs());
            if (xx/scale)*(yy/scale)<=(xy/scale).powi(2) {
                return Err(bad("head tension admits slack/compression or an unresolved positive-definite margin"));
            }
        }
        Ok(())
    }
    /// Add ONLY tensor geometric stiffness to a zero-prestress DKT pencil.
    /// P1 gradients are constant: the centroid value integrates an affine
    /// stress field exactly on each triangle. Mass, slopes and support map are
    /// untouched. The disk constructor owns this cold, one-time operation.
    pub(super) fn add_to(self, base: f64, mesh: &PlateMesh, model: &mut PlateModel) -> Result<(), PlateError> {
        self.validate_mesh(base,mesh)?;
        let n=model.free;
        let mut k=Coo::new(n,n);
        for r in 0..n {
            let (columns,values)=model.k.row(r);
            for (&c,&value) in columns.iter().zip(values) { k.push(r,c,value); }
        }
        for (element,tri) in mesh.tris.iter().enumerate() {
            let [(x0,y0),(x1,y1),(x2,y2)]=tri.map(|i|mesh.nodes[i]);
            let twice_area=(x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
            if !twice_area.is_finite() || twice_area<=0.0 {
                return Err(PlateError::DegenerateElement {element,twice_area});
            }
            let grad=[[(y1-y2)/twice_area,(x2-x1)/twice_area],
                [(y2-y0)/twice_area,(x0-x2)/twice_area],
                [(y0-y1)/twice_area,(x1-x0)/twice_area]];
            let [xx,yy,xy]=self.resultant(base,[x0/3.0+x1/3.0+x2/3.0,y0/3.0+y1/3.0+y2/3.0]);
            for i in 0..3 { for j in i..3 {
                let [ix,iy]=grad[i]; let [jx,jy]=grad[j];
                let value=0.5*twice_area*(xx*ix*jx+xy*(ix*jy+iy*jx)+yy*iy*jy);
                if !value.is_finite() { return Err(bad("head tensor geometric stiffness overflow")); }
                if let (Some(r),Some(c))=(model.dof_map[3*tri[i]],model.dof_map[3*tri[j]]) {
                    k.push(r,c,value); if r!=c { k.push(c,r,value); }
                }
            }}
        }
        let candidate=k.assemble();
        for r in 0..n {
            if candidate.row(r).1.iter().any(|v|!v.is_finite()) {
                return Err(bad("assembled head prestress overflows the stiffness pencil"));
            }
        }
        model.k=candidate;
        Ok(())
    }
}

#[cfg(test)]
#[path="tension_tests.rs"]
mod tests;
