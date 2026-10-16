//! Separable, source-selectable nonlinear membrane response.
//!
//! Follows FrankenTextiles' explicit `tabulated_green_strain/1` piecewise-linear
//! image: second-Piola resultants versus [E11,E22,2E12]. This is NOT a nominal
//! force/engineering-extension curve. Work is integrated analytically from the
//! stress-free origin; no secant stiffness substitutes for the actual tangent.
//! Cubic interpolation, hysteresis and extrapolation are separate capabilities.

use super::{SheetError, SheetLaw, SheetResponse};

/// Signed piecewise-linear response with exact work and one-sided tangents.
/// At interior knots the right-hand tangent is selected; at the final endpoint
/// it is the left-hand tangent. Energy is C1, not C2 at slope changes.
#[derive(Clone, Debug, PartialEq)]
pub struct StrainCurve {
    points: Vec<(f64, f64)>,
    slopes: Vec<f64>,
    work: Vec<f64>,
    minimum: f64,
    support: [f64;2],
}
impl StrainCurve {
    /// At most 128 ordered (Green strain, second-Piola resultant [N/m]) pairs.
    /// Both signs and an explicit (0,0) point are required. Strict positive
    /// segment slopes are a constitutive admission rule, not material accuracy.
    pub fn new(points: &[(f64, f64)]) -> Result<Self, SheetError> {
        if !(3..=128).contains(&points.len())
            || points.iter().any(|(x,y)| !x.is_finite() || !y.is_finite()) {
            return Err(SheetError::Invalid("strain curve needs 3..128 finite pairs"));
        }
        let origin = points.iter().position(|p| *p == (0.0,0.0))
            .ok_or(SheetError::Invalid("strain curve requires an explicit stress-free origin"))?;
        if origin == 0 || origin+1 == points.len() {
            return Err(SheetError::Invalid("strain curve must characterize compression and tension"));
        }
        let mut slopes = Vec::with_capacity(points.len()-1);
        for w in points.windows(2) {
            let dx=w[1].0-w[0].0; let dy=w[1].1-w[0].1; let slope=dy/dx;
            if !dx.is_finite() || dx<=0.0 || !dy.is_finite() || dy<=0.0
                || !slope.is_finite() || slope<=0.0 {
                return Err(SheetError::Invalid("strain and resultant knots must strictly increase with representable slopes"));
            }
            slopes.push(slope);
        }
        // Integrate OUTWARD from zero on both sides. This avoids subtracting a
        // large arbitrary energy datum when evaluating near the unstressed state.
        let mut work=vec![0.0;points.len()];
        for i in origin..points.len()-1 {
            work[i+1]=work[i]+(points[i+1].0-points[i].0)*(0.5*points[i].1+0.5*points[i+1].1);
        }
        for i in (0..origin).rev() {
            work[i]=work[i+1]-(points[i+1].0-points[i].0)*(0.5*points[i].1+0.5*points[i+1].1);
        }
        if work.iter().any(|w| !w.is_finite() || *w<0.0) {
            return Err(SheetError::Invalid("unrepresentable strain-curve work"));
        }
        let minimum=slopes.iter().copied().fold(f64::INFINITY,f64::min);
        Ok(Self {points:points.to_vec(),slopes,work,minimum,
            support:[points[0].0,points[points.len()-1].0]})
    }
    /// Original knots, retained without resampling or fitting.
    pub fn points(&self) -> &[(f64,f64)] {&self.points}
    /// Closed supported strain interval.
    pub const fn domain(&self) -> [f64;2] {self.support}
    /// Minimum segment tangent [N/m].
    pub const fn minimum_tangent(&self) -> f64 {self.minimum}
    /// Require a two-sided local tangent for small-vibration extraction.
    /// Newton may still use evaluate() at a kink. A conservative 8-epsilon
    /// strain neighborhood also refuses numerically unresolved knot proximity;
    /// this numerical screen is not an interval material-accuracy certificate.
    pub fn admit_linearization(&self, strain:f64) -> Result<(),SheetError> {
        self.evaluate(strain)?;
        let guard=8.0*f64::EPSILON*strain.abs().max(1.0);
        let [lo,hi]=self.domain();
        if strain-lo<=guard || hi-strain<=guard {
            return Err(SheetError::Invalid("small vibration needs two-sided source strain support"));
        }
        for i in 1..self.points.len()-1 {
            if (strain-self.points[i].0).abs()<=guard && self.slopes[i-1]!=self.slopes[i] {
                return Err(SheetError::Invalid("no unique small-vibration tangent at a slope-changing source knot"));
            }
        }
        Ok(())
    }
    /// [work J/m2, resultant N/m, tangent N/m] at the SAME strain.
    pub fn evaluate(&self, strain:f64) -> Result<[f64;3],SheetError> {
        let [lo,hi]=self.domain();
        if !strain.is_finite() || strain<lo || strain>hi {
            return Err(SheetError::Invalid("Green strain outside the source curve; no extrapolation"));
        }
        let i=self.points.partition_point(|p|p.0<=strain).saturating_sub(1).min(self.slopes.len()-1);
        let anchor=if strain<0.0 {i+1} else {i};
        let dx=strain-self.points[anchor].0; let y=self.points[anchor].1; let slope=self.slopes[i];
        let out=[self.work[anchor]+dx*(y+0.5*slope*dx), y+slope*dx,slope];
        if out.iter().any(|v|!v.is_finite()) || out[0]<0.0 {
            return Err(SheetError::Invalid("unrepresentable strain-curve response"));
        }
        Ok(out)
    }
}

/// Independent nonlinear warp, weft and shear curves with constant cross
/// derivative. No scalar membrane stiffness is added to these curves.
#[derive(Clone, Debug, PartialEq)]
pub struct TabulatedMembrane {
    curves: [StrainCurve;3],
    coupling_n_m: f64,
}
impl TabulatedMembrane {
    /// Cross-coupling must keep the tangent positive definite throughout the
    /// whole rectangular strain domain, not merely at the current sample.
    /// An eight-epsilon relative margin conservatively excludes unresolved
    /// proximity to singular coupling; it is not an empirical error bound.
    pub fn new(curves:[StrainCurve;3], coupling_n_m:f64) -> Result<Self,SheetError> {
        let ratio=coupling_n_m/curves[0].minimum_tangent().sqrt()/curves[1].minimum_tangent().sqrt();
        // A singular coupling can round just below one after two divisions.
        // Reserve a small numerical margin rather than admitting that boundary.
        if !coupling_n_m.is_finite() || !ratio.is_finite()
            || ratio.abs()>=1.0-8.0*f64::EPSILON {
            return Err(SheetError::Invalid("nonlinear membrane coupling is not strictly stable over its source domain"));
        }
        Ok(Self {curves,coupling_n_m})
    }
    /// Require each source direction to admit a two-sided local tangent.
    pub fn admit_linearization(&self, strain:[f64;3]) -> Result<(),SheetError> {
        for (curve,e) in self.curves.iter().zip(strain) {curve.admit_linearization(e)?;}
        Ok(())
    }
    /// Original constitutive curves in material-axis order.
    pub const fn curves(&self) -> &[StrainCurve;3] {&self.curves}
    /// Constant dS11/dE22 = dS22/dE11 [N/m].
    pub const fn coupling_n_m(&self) -> f64 {self.coupling_n_m}
}
impl SheetLaw for TabulatedMembrane {
    fn response(&self, strain:[f64;3]) -> Result<SheetResponse,SheetError> {
        let a=self.curves[0].evaluate(strain[0])?;
        let b=self.curves[1].evaluate(strain[1])?;
        let s=self.curves[2].evaluate(strain[2])?;
        let c=self.coupling_n_m;
        let energy=a[0]+b[0]+s[0]+c*strain[0]*strain[1];
        let stress=[a[1]+c*strain[1],b[1]+c*strain[0],s[1]];
        if !energy.is_finite() || energy<0.0 || stress.iter().any(|x|!x.is_finite()) {
            return Err(SheetError::Invalid("unrepresentable nonlinear membrane work or resultant"));
        }
        Ok(SheetResponse {energy_j_m2:energy,resultants_n_m:stress,
            tangent_n_m:[a[2],c,0.0,c,b[2],0.0,0.0,0.0,s[2]]})
    }
}

mod source;
pub use source::{STRAIN_AXIS, CURVE_PROPERTIES, COUPLING_PROPERTY};
