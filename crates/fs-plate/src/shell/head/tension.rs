//! Equilibrated installed membrane prestress on the existing P1 head geometry.
//!
//! The original affine field remains exact and bit-for-bit unchanged when no
//! harmonic is supplied. Optional order-2..12 circumferential terms come from
//! polynomial Airy stress potentials, so each term is divergence-free in the
//! head interior rather than an authored per-mode pitch adjustment.
//!
//! Harmonic coefficients are normalized installed-resultant amplitudes [N/m].
//! For order n and coefficient pair (a,b), the Airy image is proportional to
//! r^(n+2)*(a*cos(n theta)+b*sin(n theta)); normalization makes hypot(a,b)
//! a conservative Frobenius-norm bound on that term everywhere inside the rim.
//! This supports resolved multi-lug patterns while remaining a prestress field,
//! not a hoop/lug/bearing-edge dynamics model or a measured tension map.
use crate::{PlateError, PlateMesh, PlateModel};
use fs_sparse::Coo;

pub const MIN_TENSION_HARMONIC_ORDER: usize = 2;
pub const MAX_TENSION_HARMONIC_ORDER: usize = 12;
pub const TENSION_HARMONIC_COUNT: usize =
    MAX_TENSION_HARMONIC_ORDER - MIN_TENSION_HARMONIC_ORDER + 1;

/// Variation about a disk's explicitly supplied uniform installed tension.
/// Coordinates are the original head's Cartesian metres, centered on the disk.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TensionVariation {
    /// Constant additions [N/m] in the order xx, yy, xy (tensor shear).
    pub constant_n_m: [f64; 3],
    /// Spatial derivatives [N/m²] in the order c, d, e, f in this module's law.
    pub gradient_n_m2: [f64; 4],
    /// Normalized Airy harmonic coefficients [cos,sin] in N/m for orders 2..=12.
    pub harmonics_n_m: [[f64; 2]; TENSION_HARMONIC_COUNT],
}
fn bad(what: &'static str) -> PlateError { PlateError::BadSection { what } }

const GL8_X: [f64;8] = [
    0.0198550717512319,0.101666761293187,0.237233795041836,0.408282678752175,
    0.591717321247825,0.762766204958164,0.898333238706813,0.980144928248768,
];
const GL8_W: [f64;8] = [
    0.0506142681451881,0.111190517226687,0.156853322938944,0.181341891689181,
    0.181341891689181,0.156853322938944,0.111190517226687,0.0506142681451881,
];

impl TensionVariation {
    /// Admit finite coefficients. Tensile admissibility additionally depends on
    /// the base tension and actual disk radius.
    pub fn validate(self) -> Result<(), PlateError> {
        if self.constant_n_m.iter().chain(&self.gradient_n_m2)
            .chain(self.harmonics_n_m.iter().flatten()).any(|v| !v.is_finite()) {
            return Err(bad("head tension variation needs finite SI coefficients"));
        }
        Ok(())
    }
    /// The original uniform image, including signed-zero input.
    #[must_use]
    pub fn is_zero(self) -> bool {
        self.constant_n_m.iter().chain(&self.gradient_n_m2)
            .chain(self.harmonics_n_m.iter().flatten()).all(|v| *v == 0.0)
    }
    #[must_use]
    pub fn has_harmonics(self) -> bool {
        self.harmonics_n_m.iter().flatten().any(|v| *v != 0.0)
    }
    /// Affine installed resultant [N/m] in xx, yy, xy order.
    #[must_use]
    pub fn resultant(self, base_n_m: f64, point_m: [f64; 2]) -> [f64; 3] {
        let [x,y]=point_m; let [xx,yy,xy]=self.constant_n_m;
        let [c,d,e,f]=self.gradient_n_m2;
        [base_n_m+xx+c*x+d*y, base_n_m+yy+e*x+f*y, xy-f*x-c*y]
    }
    /// Full affine + Airy-harmonic resultant on a disk of the declared radius.
    #[must_use]
    pub fn resultant_on_disk(self, base_n_m:f64, radius_m:f64, point_m:[f64;2])->[f64;3] {
        let mut out=self.resultant(base_n_m,point_m);
        if !self.has_harmonics() { return out; }
        let [x,y]=point_m; let r=x.hypot(y);
        if r==0.0 { return out; }
        let theta=y.atan2(x); let rho=r/radius_m;
        let (st,ct)=theta.sin_cos();
        for (slot,pair) in self.harmonics_n_m.iter().enumerate() {
            let [ac,as_]=*pair;
            if ac==0.0 && as_==0.0 {continue;}
            let n=(slot+MIN_TENSION_HARMONIC_ORDER) as f64;
            let angle=n*theta; let (sn,cn)=angle.sin_cos();
            let aa=n+2.0-n*n;
            let bb=(n+2.0)*(n+1.0);
            let cc=n*(n+1.0);
            let norm=(aa*aa+bb*bb+2.0*cc*cc).sqrt();
            let radial=rho.powi(n as i32)/norm;
            let s=radial*(ac*cn+as_*sn);
            let t=radial*(ac*sn-as_*cn);
            let rr=aa*s; let tt=bb*s; let rt=cc*t;
            let c2=ct*ct; let s2=st*st; let cs=ct*st;
            out[0]+=c2*rr+s2*tt-2.0*cs*rt;
            out[1]+=s2*rr+c2*tt+2.0*cs*rt;
            out[2]+=cs*(rr-tt)+(c2-s2)*rt;
        }
        out
    }
    fn radius(mesh:&PlateMesh)->Result<f64,PlateError> {
        let radius=mesh.nodes.iter().map(|&(x,y)|x.hypot(y)).fold(0.0_f64,f64::max);
        if !radius.is_finite() || radius<=0.0 {
            return Err(bad("head tension harmonics need a finite positive disk radius"));
        }
        Ok(radius)
    }
    fn min_eigenvalue([xx,yy,xy]:[f64;3])->f64 {
        let scale=xx.abs().max(yy.abs()).max(xy.abs());
        if scale==0.0 {return 0.0;}
        scale*(0.5*(xx/scale+yy/scale)
            - (0.25*(xx/scale-yy/scale).powi(2)+(xy/scale).powi(2)).sqrt())
    }
    fn harmonic_bound(self)->f64 {
        self.harmonics_n_m.iter().map(|v|v[0].hypot(v[1])).sum()
    }
    pub(super) fn validate_mesh(self, base: f64, mesh: &PlateMesh) -> Result<(), PlateError> {
        self.validate()?;
        if !base.is_finite() || base<=0.0 { return Err(bad("head base tension must be finite and positive")); }
        let radius=Self::radius(mesh)?;
        let mut affine_margin=f64::INFINITY;
        for &(x,y) in &mesh.nodes {
            if x.hypot(y)>radius*(1.0+64.0*f64::EPSILON) {
                return Err(bad("head tension harmonic point lies outside its disk radius"));
            }
            let affine=self.resultant(base,[x,y]);
            if affine.iter().any(|v|!v.is_finite()) {
                return Err(bad("installed head tensor must remain finite across the mesh"));
            }
            affine_margin=affine_margin.min(Self::min_eigenvalue(affine));
        }
        // lambda_min is concave on symmetric matrices, so the minimum affine
        // margin is bounded by the vertices. Each normalized harmonic has
        // spectral norm <= hypot(cos,sin); the sum therefore gives a rigorous
        // disk-wide positivity guard without pretending sampled points certify it.
        let bound=self.harmonic_bound();
        if !affine_margin.is_finite() || affine_margin<=0.0 || !bound.is_finite()
            || bound>=affine_margin {
            return Err(bad("head tension harmonics exceed the certified tensile margin"));
        }
        Ok(())
    }
    fn integrated_resultant(self,base:f64,radius:f64,p:[[f64;2];3],twice_area:f64)->[f64;3] {
        if !self.has_harmonics() {
            let point=[(p[0][0]+p[1][0]+p[2][0])/3.0,(p[0][1]+p[1][1]+p[2][1])/3.0];
            return self.resultant(base,point).map(|v|0.5*twice_area*v);
        }
        // Duffy transform of 8x8 Gauss-Legendre. The Airy stresses are
        // Cartesian polynomials of degree <=12; after the transform and its
        // Jacobian the degree is <=13 in either coordinate, so GL8 integrates
        // this family exactly up to floating-point evaluation error.
        let mut out=[0.0;3];
        for (iu,&u) in GL8_X.iter().enumerate() {
            for (it,&t) in GL8_X.iter().enumerate() {
                let b1=u; let b2=(1.0-u)*t; let b0=1.0-b1-b2;
                let point=[b0*p[0][0]+b1*p[1][0]+b2*p[2][0],
                    b0*p[0][1]+b1*p[1][1]+b2*p[2][1]];
                let n=self.resultant_on_disk(base,radius,point);
                let w=GL8_W[iu]*GL8_W[it]*(1.0-u)*twice_area;
                for j in 0..3 {out[j]+=w*n[j];}
            }
        }
        out
    }
    /// Add ONLY tensor geometric stiffness to a zero-prestress DKT pencil.
    /// Affine-only input retains the original centroid arithmetic. Harmonics use
    /// exact-order cold quadrature for the declared polynomial Airy family.
    pub(super) fn add_to(self, base: f64, mesh: &PlateMesh, model: &mut PlateModel) -> Result<(), PlateError> {
        self.validate_mesh(base,mesh)?;
        let radius=Self::radius(mesh)?;
        let n=model.free;
        let mut k=Coo::new(n,n);
        for r in 0..n {
            let (columns,values)=model.k.row(r);
            for (&c,&value) in columns.iter().zip(values) { k.push(r,c,value); }
        }
        for (element,tri) in mesh.tris.iter().enumerate() {
            let p=tri.map(|i|{let (x,y)=mesh.nodes[i];[x,y]});
            let [[x0,y0],[x1,y1],[x2,y2]]=p;
            let twice_area=(x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
            if !twice_area.is_finite() || twice_area<=0.0 {
                return Err(PlateError::DegenerateElement {element,twice_area});
            }
            let grad=[[(y1-y2)/twice_area,(x2-x1)/twice_area],
                [(y2-y0)/twice_area,(x0-x2)/twice_area],
                [(y0-y1)/twice_area,(x1-x0)/twice_area]];
            let [xx,yy,xy]=self.integrated_resultant(base,radius,p,twice_area);
            for i in 0..3 { for j in i..3 {
                let [ix,iy]=grad[i]; let [jx,jy]=grad[j];
                let value=xx*ix*jx+xy*(ix*jy+iy*jx)+yy*iy*jy;
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
