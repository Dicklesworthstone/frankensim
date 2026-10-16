//! Observe the very same P1 double layer/P0 single layer used in boundary work.
use super::*;

/// One admitted Galerkin solve. Its source geometry is borrowed, not inferred
/// from an independently supplied panel-pressure array.
pub struct Solution<'p,'s> {
    pub(super) owner:&'p Prepared<'s>, pub(super) coefficients:Vec<C64>,
    /// Wavenumber [1/m].
    pub k:f64,
    /// Actual homogeneous exterior gas.
    pub medium:Medium,
    /// Exact area mean of the P1 pressure trace on every original panel [Pa].
    /// Use these means for work against prescribed P0 normal velocity, NOT as
    /// replacement P0 densities for a different observer formulation.
    pub pressure:Vec<C64>,
    /// Original prescribed P0 outward velocity, not the auxiliary density.
    pub velocity:Vec<C64>,
    /// Complex boundary work 1/2 integral p conj(v) dS [W], left unclipped.
    pub surface_power:C64,
    /// Nonnegative flux from the same angular kernel quadrature [W]. This is
    /// not a bound on the continuum radiation or quadrature truncation error.
    pub outgoing_power_w:f64,
    /// Absolute flux/work defect divided by apparent-power scale.
    pub flux_balance:f64,
    /// Componentwise-scaled infinity residual of the solved mixed operator.
    pub backward_residual:f64,
    /// Same area-based wavelength diagnostic as the collocation image.
    pub panels_per_wavelength:f64,
}
impl Solution<'_, '_> {
    /// Angular amplitude F with p ~ F exp(ikr)/r. Direction is normalized here.
    /// # Errors
    /// Nonfinite/zero direction or unrepresentable field arithmetic.
    pub fn far_field(&self,direction:Point)->Result<C64,HelmholtzError> {
        let radius=norm(direction);
        if direction.iter().any(|v|!v.is_finite()) || !radius.is_finite() || radius<=0.0 {
            return Err(bad("invalid reciprocal far-field direction"));
        }
        let coefficients=self.owner.wave_coefficients(self.k,direction.map(|v|v/radius))?;
        let value=coefficients.iter().zip(&self.coefficients)
            .fold(C64::ZERO,|v,(a,x)|v+*a * *x).scale(1.0/(4.0*PI));
        if !finite(value) {return Err(bad("reciprocal far-field overflow"));}Ok(value)
    }
    /// Pressure at a genuinely exterior point [Pa], including physical phase
    /// and spreading. Reuse the existing near-field geometry admission, then
    /// integrate the P1/P0 representation (analytic static moments plus 4x4
    /// regular remainder). This fixed rule is NOT the adaptive near-mic image.
    /// # Errors
    /// Invalid, interior, surface-touching or arithmetically unresolved receiver.
    pub fn pressure_at(&self,position:Point)->Result<C64,HelmholtzError> {
        let points=[position];
        let _admitted=crate::near_field::Geometry::new(self.owner.surface,&points,
            64.0*f64::EPSILON*self.owner.length).map_err(|_|bad("reciprocal receiver must be a resolved exterior point"))?;
        let mut value=C64::ZERO;
        for (i,f) in self.owner.facets.iter().enumerate() {
            let (single,_,double)=static_moments(position,f,false);
            let lambda=self.coefficients[self.owner.vertices+i];
            value=value+lambda.scale(single);
            for a in 0..3 {value=value+self.coefficients[f.ids[a]].scale(double[a]);}
            for sample in f.samples {
                let d=sub(position,sample.point);let r=norm(d);let z=self.k*r;
                let (re_single,re_derivative)=regular_real(self.k,r);
                let im_derivative=if z.abs()<1e-3 {
                    let z2=z*z;z*z2*(1.0/3.0+z2*(-1.0/30.0+z2/840.0))
                }else{det::sin(z)-z*det::cos(z)};
                let g=C64::new(re_single,det::sin(z)/(4.0*PI*r));
                let dg=C64::new(re_derivative,im_derivative/(4.0*PI*r*r*r)).scale(dot(f.normal,d));
                let p=(0..3).fold(C64::ZERO,|sum,a|sum+self.coefficients[f.ids[a]].scale(sample.shape[a]));
                value=value+(g*lambda+dg*p).scale(sample.weight);
            }
        }
        if !finite(value) {return Err(bad("reciprocal point-pressure overflow"));}Ok(value)
    }
}
