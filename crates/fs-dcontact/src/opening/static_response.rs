//! Differential of the existing one-point normal law at zero opening speed.
//! No new potential, preload algorithm or force clipping is introduced.
use super::*;

/// Partial derivatives at a fixed strictly active or strictly separated opening.
/// These are derivatives of the supplied law, not uncertainty or material data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StaticContactDifferential {
    /// Positive penetration, or negative separation [m].
    pub signed_penetration_m: f64,
    /// Existing conservative normal reaction [N].
    pub force_n: f64,
    /// Derivative with respect to increasing closure, opposite opening [N/m].
    pub closure_stiffness_n_m: f64,
    /// Derivative with respect to the original stiffness coefficient.
    pub force_per_stiffness: f64,
    /// Derivative with respect to the original quadrature weight.
    pub force_per_weight: f64,
}

impl OpeningContactStep {
    /// Differentiate the stationary reaction at this step's initial opening.
    /// A caller-supplied exclusion distance [m] refuses contact switching points,
    /// including exactly touching. It does NOT certify that a parameter step
    /// stays in this activity region. Zero stiffness/weight retain their proper
    /// one-sided coefficient partials; no division by those coefficients occurs.
    /// Internal loss contributes zero at rest and has no static derivative here.
    pub fn static_differential(&self, minimum_margin_m: f64)
        -> Result<StaticContactDifferential, DContactError>
    {
        if !minimum_margin_m.is_finite() || minimum_margin_m < 0.0 { return Err(bad()); }
        let penetration = -self.before - self.gap;
        if !penetration.is_finite() || penetration.abs() <= minimum_margin_m {
            return Err(DContactError::Parameter { what: "static contact derivative is inside the activity-switch exclusion distance" });
        }
        let mut out = StaticContactDifferential {
            signed_penetration_m: penetration, force_n: 0.0,
            closure_stiffness_n_m: 0.0, force_per_stiffness: 0.0, force_per_weight: 0.0,
        };
        if penetration > 0.0 {
            // Evaluate force through the existing stable scalar owner. The
            // derivative uses the same deterministic power at exponent alpha-1.
            out.force_n = self.coefficients(self.before)?.0;
            let power = det::pow(penetration, self.alpha);
            out.closure_stiffness_n_m = self.weight * self.stiffness
                * self.alpha * det::pow(penetration, self.alpha - 1.0);
            out.force_per_stiffness = self.weight * power;
            out.force_per_weight = self.stiffness * power;
        }
        if [out.force_n, out.closure_stiffness_n_m, out.force_per_stiffness,
            out.force_per_weight].iter().any(|v| !v.is_finite() || *v < 0.0) { return Err(bad()); }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn point(k:f64,w:f64,gap:f64,alpha:f64,opening:f64)->OpeningContactStep {
        let law=Obstacle::new(vec![-1.0],1,1,vec![gap],vec![w],k,alpha,"differential test".into()).unwrap();
        OpeningContactStep::new(&law,opening).unwrap()
    }
    #[test]
    fn stationary_partials_match_independent_polynomial_and_difference_checks() {
        for alpha in [1.0,1.5,2.0,3.7] {
            let p=point(3e6,0.7,0.0001,alpha,-0.002);
            let d=p.static_differential(1e-8).unwrap();
            let h=1e-8;
            let force=|k,w,gap,opening|point(k,w,gap,alpha,opening).coefficients(opening).unwrap().0;
            let closure=(force(3e6,0.7,0.0001,-0.002-h)-force(3e6,0.7,0.0001,-0.002+h))/(2.0*h);
            assert!((d.closure_stiffness_n_m-closure).abs()<1e-7*d.closure_stiffness_n_m);
            assert!((d.force_per_stiffness*3e6-d.force_n).abs()<1e-13*d.force_n);
            assert!((d.force_per_weight*0.7-d.force_n).abs()<1e-13*d.force_n);
        }
    }
    #[test]
    fn coefficient_zeros_keep_nonzero_coefficient_partials_and_separation_is_zero() {
        let d=point(0.0,2.0,0.0,2.0,-0.01).static_differential(0.0).unwrap();
        assert_eq!(d.force_n,0.0);assert_eq!(d.closure_stiffness_n_m,0.0);
        assert!((d.force_per_stiffness-0.0002).abs()<1e-16);
        assert!(point(10.0,0.0,0.0,2.0,-0.01).static_differential(0.0).unwrap().force_per_weight>0.0);
        let d=point(10.0,2.0,0.0,2.0,0.01).static_differential(0.0).unwrap();
        assert_eq!((d.force_n,d.closure_stiffness_n_m,d.force_per_stiffness,d.force_per_weight),(0.0,0.0,0.0,0.0));
    }
    #[test]
    fn switching_and_invalid_margins_refuse() {
        assert!(point(1.0,1.0,0.0,1.0,0.0).static_differential(0.0).is_err());
        assert!(point(1.0,1.0,0.0,2.0,-1e-9).static_differential(1e-8).is_err());
        for margin in [-1.0,f64::NAN,f64::INFINITY] {
            assert!(point(1.0,1.0,0.0,2.0,-0.01).static_differential(margin).is_err());
        }
    }
}
