//! Stable one-coordinate secant of the existing unilateral power-law potential.
//! This is the scalar discrete gradient, evaluated without subtracting almost
//! equal energies. It changes rounding, not the law or the time integrator.
use crate::{DContactError, Obstacle};
use fs_math::det;

/// One fixed initial opening and its admitted fs-dcontact coefficients.
/// The obstacle must use the one-point unit opening map [-1]. Evaluation is
/// allocation-free; no state, energy, or reaction is changed by a trial.
#[derive(Clone, Copy, Debug)]
pub struct OpeningContactStep {
    before: f64,
    gap: f64,
    stiffness: f64,
    weight: f64,
    alpha: f64,
    loss: f64,
}
impl OpeningContactStep {
    /// Copy the original law, refusing malformed raw-parts obstacles before any
    /// indexing. This numerical adapter makes no source-validity assertion.
    pub fn new(law: &Obstacle, opening_before_m: f64) -> Result<Self, DContactError> {
        if law.n_points() != 1 || law.collocation() != [-1.0]
            || law.gaps().len() != 1 || law.weights().len() != 1 {
            return Err(DContactError::Shape { what: "opening secant requires one unit-opening contact" });
        }
        if !opening_before_m.is_finite() || !law.gaps()[0].is_finite()
            || [law.stiffness(), law.weights()[0], law.internal_loss()].iter().any(|v| !v.is_finite() || *v < 0.0)
            || !law.alpha().is_finite() || law.alpha() < 1.0 {
            return Err(bad());
        }
        Ok(Self { before: opening_before_m, gap: law.gaps()[0], stiffness: law.stiffness(),
            weight: law.weights()[0], alpha: law.alpha(), loss: law.internal_loss() })
    }

    /// Positive conservative reaction and its Hunt-Crossley damping coefficient.
    /// At opening velocity v use max(elastic - damping*v, 0), as in the existing
    /// nonadhesive law. Entry and exit use the full opening displacement, not
    /// only the distance spent inside contact. Coincident endpoints use the
    /// original potential derivative. Refuses nonfinite arithmetic, never clips.
    pub fn coefficients(&self, after: f64) -> Result<(f64, f64), DContactError> {
        if !after.is_finite() { return Err(bad()); }
        let lo = -self.before.max(after) - self.gap;
        let hi = -self.before.min(after) - self.gap;
        let delta = (after-self.before).abs();
        if !lo.is_finite() || !hi.is_finite() || !delta.is_finite() { return Err(bad()); }
        if hi <= 0.0 || self.stiffness == 0.0 || self.weight == 0.0 { return Ok((0.0, 0.0)); }
        let amplitude = self.weight*self.stiffness*det::pow(hi, self.alpha);
        let beta = self.alpha+1.0;
        let elastic = if delta == 0.0 {
            amplitude
        } else if lo <= 0.0 {
            amplitude*(hi/delta)/beta
        } else {
            let r = delta/hi;
            if r == 0.0 { amplitude } else {
                // log(1-r) = -2 atanh(r/(2-r)). Near coincidence this fixed
                // eight-term series avoids forming a ratio rounded to one.
                // For r<0.01 its analytic tail is below 2e-40; floating-point
                // rounding remains, and no interval certificate is asserted.
                let log_ratio = if r < 0.01 {
                    let z = r/(2.0-r);
                    let mut term = z;
                    let mut sum = z;
                    for n in 1..8 { term *= z*z; sum += term/(2*n+1) as f64; }
                    -2.0*sum
                } else { det::ln(lo/hi) };
                amplitude*(-det::expm1(beta*log_ratio))/(beta*r)
            }
        };
        let damping = elastic*self.loss;
        if !elastic.is_finite() || elastic < 0.0 || !damping.is_finite() || damping < 0.0 { return Err(bad()); }
        Ok((elastic, damping))
    }
}
fn bad() -> DContactError { DContactError::Parameter { what: "opening contact secant requires finite physical coefficients and openings" } }

#[cfg(test)]
mod tests {
    use super::*;
    fn law(alpha:f64,gap:f64)->Obstacle {
        Obstacle::new(vec![-1.0],1,1,vec![gap],vec![1.0],1e8,alpha,"authored secant test".into())
            .unwrap().with_internal_loss(0.2).unwrap()
    }
    #[test]
    fn tiny_motion_matches_polynomial_without_subtracting_energies() {
        let gap=0.0001;let a=0.000569113194841233_f64;
        for relative in [0.0,1e-15,1e-12,1e-9,1e-6,0.001] {
            let b=a*(1.0-relative);let p=a-gap;let q=b-gap;
            let expected=1e8*(p*p+p*q+q*q)/3.0;
            let (force,damping)=OpeningContactStep::new(&law(2.0,gap),-a).unwrap().coefficients(-b).unwrap();
            assert!((force-expected).abs()<1e-13*expected,"{relative}: {force} != {expected}");
            assert_eq!(damping,force*0.2);
        }
    }
    #[test]
    fn arbitrary_exponents_are_symmetric_and_approach_the_original_derivative() {
        for alpha in [1.0,1.01,1.5,2.5,3.7,8.0,20.0] {
            let ob=law(alpha,0.0);let a=-0.01;let b=a*(1.0-1e-12);
            let forward=OpeningContactStep::new(&ob,a).unwrap().coefficients(b).unwrap();
            let reverse=OpeningContactStep::new(&ob,b).unwrap().coefficients(a).unwrap();
            assert_eq!(forward,reverse);
            let derivative=1e8*det::pow(-a,alpha);
            assert!((forward.0-derivative).abs()<3e-11*derivative);
        }
    }
    #[test]
    fn crossing_accounts_for_the_entire_opening_and_noncontact_is_zero() {
        let ob=law(2.0,0.0);
        let f=OpeningContactStep::new(&ob,0.002).unwrap().coefficients(-0.001).unwrap().0;
        assert!((f-1e8*0.001_f64.powi(3)/(3.0*0.003)).abs()<1e-12);
        assert_eq!(OpeningContactStep::new(&ob,0.002).unwrap().coefficients(0.001).unwrap(),(0.0,0.0));
        assert_eq!(OpeningContactStep::new(&ob,0.0).unwrap().coefficients(0.0).unwrap(),(0.0,0.0));
    }
    #[test]
    fn bad_raw_shapes_and_nonfinite_trials_refuse_without_indexing_or_mutation() {
        let bad_shape=Obstacle::from_raw_parts(vec![-1.0],1,vec![],vec![1.0],1e8,2.0,"bad".into());
        assert!(matches!(OpeningContactStep::new(&bad_shape,0.0),Err(DContactError::Shape {..})));
        assert!(OpeningContactStep::new(&law(2.0,0.0),f64::NAN).is_err());
        let step=OpeningContactStep::new(&law(2.0,0.0),-0.01).unwrap();
        let expected=step.coefficients(-0.01).unwrap();
        assert!(step.coefficients(f64::INFINITY).is_err());assert_eq!(step.coefficients(-0.01).unwrap(),expected);
    }
}
