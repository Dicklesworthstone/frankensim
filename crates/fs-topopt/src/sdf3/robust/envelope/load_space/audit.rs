//! Cancellation-resistant values for the existing directional gradient audit.
use super::{failure, EvaluationStop, LoadSpace3, LoadSpaceSpectrum3};

/// J - temperature*ln(n), computed without first rounding away tiny eigenvalues.
/// Its derivative is exactly J's derivative. This value is ONLY used by the
/// numerical audit: the optimizer, spectrum, physical caps and weights do not
/// change, and no absolute error floor is introduced.
pub(super) fn centered_maximum(s: &LoadSpaceSpectrum3, p: &LoadSpace3)
    -> Result<f64, EvaluationStop> {
    let magnitude = s.eigenvalues.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let value = if magnitude > 0.25*p.temperature() {
        // Outside the offset-dominated regime, the original stable LSE is
        // already resolved at the scale of its centered value.
        s.smooth_maximum-p.bias_bound()
    } else {
        // With m=max(lambda), x=mean(expm1((lambda-m)/t)),
        // J-t*ln(n) = m+t*log1p(x). Evaluate t*x in eigenvalue units:
        // delta*expm1(delta/t)/(delta/t) also survives delta/t underflow.
        let n = p.dimension() as f64;
        let mut x = 0.0;
        let mut tx = 0.0;
        for &lambda in &s.eigenvalues {
            let delta = lambda-s.maximum;
            let z = delta/p.temperature();
            let e = fs_math::det::expm1(z);
            x += e/n;
            tx += (delta/n)*if z == 0.0 { 1.0 } else { e/z };
        }
        // Correct for rounding in 1+x using ln(1+x)/(rounded(1+x)-1).
        // At rounded(1+x)==1 the analytic limiting ratio is one. Only the
        // existing deterministic math owner is used; no platform log1p.
        let u = 1.0+x;
        let ratio = if u == 1.0 { 1.0 } else { fs_math::det::ln(u)/(u-1.0) };
        s.maximum+tx*ratio
    };
    if value.is_finite() { Ok(value) }
    else { Err(failure("sdf3-load-space-centered-audit")) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{spectrum, controlled_load_space_gradient_check, RobustOptions3};
    use crate::gradient_check::{audit_design_map, GradientCheckOptions, GradientSample};
    use crate::{SimpParams, SolveBudget, SolveControl};
    use std::ops::ControlFlow;

    fn model(x: f64, scale: f64, p: &LoadSpace3) -> (LoadSpaceSpectrum3, f64) {
        let n = p.dimension();
        let mut gram = vec![0.0; n*n];
        for i in 0..n { gram[i*n+i] = scale*(1.0+(i+1) as f64/n as f64)/(1.0+x); }
        let s = spectrum(gram, p).unwrap();
        let derivative = s.weights.iter().zip(&s.eigenvalues).map(|(w,v)| -w*v/(1.0+x)).sum();
        (s, derivative)
    }

    #[test]
    fn g0_centered_audit_resolves_small_sensitivities_without_changing_the_objective() {
        for n in [1,2,16] { for scale in [1e-200,1e-20,1e-12,1e-8,1.0,1e100,1e200] {
            for temperature in [0.02,1e200] {
                let p = LoadSpace3::new(n,1.0,temperature).unwrap();
                let (x,h) = (0.6,1e-4);
                let (s,analytic) = model(x,scale,&p);
                let before = (s.smooth_maximum,s.maximum,s.weights.clone());
                let base = centered_maximum(&s,&p).unwrap();
                let first = centered_maximum(&model(x+h,scale,&p).0,&p).unwrap();
                let second = centered_maximum(&model(x+2.0*h,scale,&p).0,&p).unwrap();
                let difference = (2.0*(first-base)-0.5*(second-base))/h;
                let relative = (difference/analytic-1.0).abs();
                assert!(relative < GradientCheckOptions::default().relative_tolerance,
                    "n={n} scale={scale} temperature={temperature} discrepancy={relative}");
                assert_eq!(before,(s.smooth_maximum,s.maximum,s.weights.clone()));
            }
        }}
    }

    #[test]
    fn g0_centering_does_not_hide_zero_or_wrong_sign_derivatives() {
        let policy = LoadSpace3::new(2,1.0,0.02).unwrap();
        for factor in [0.0,-1.0,1.0] {
            let mut poll = |_| ControlFlow::Continue(());
            let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
            let check = audit_design_map(SimpParams::default(),&[0.6],GradientCheckOptions::default(),
                &mut control,|point,_,_| {
                    let (s,derivative) = model(point[0],1e-100,&policy);
                    Ok(GradientSample { compliance:centered_maximum(&s,&policy)?,
                        compliance_gradient:vec![factor*derivative],volume:point[0],volume_gradient:vec![1.0] })
                }).unwrap();
            assert_eq!(check.passed(),factor == 1.0,"factor={factor}: {check:?}");
        }
        let zero = spectrum(vec![0.0;4],&policy).unwrap();
        assert_eq!(centered_maximum(&zero,&policy).unwrap(),0.0);
    }

    #[test]
    fn g1_loose_caps_keep_actual_cutfem_gradient_checks_and_original_baseline_values() {
        use crate::sdf3::robust::envelope::tests::fixture;
        let (mut study,y,z) = fixture(2.0);
        let rho = vec![0.75;study.cells()];
        let original = study.operator().scales().to_vec();
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
        for reference in [1e6,1e12,1e100] {
            let policy = LoadSpace3::new(2,reference,0.02).unwrap();
            let baseline = study.evaluate_robust_load_space(&rho,&[&y,&z],RobustOptions3::default(),
                &policy,&mut control).unwrap();
            let check = controlled_load_space_gradient_check(&mut study,&rho,&[&y,&z],RobustOptions3::default(),
                &policy,GradientCheckOptions::default(),&mut control).unwrap();
            assert!(check.numerical.passed(),"reference={reference}: {check:?}");
            assert_eq!(check.numerical.baseline_compliance,baseline.spectrum.smooth_maximum);
            assert_eq!(study.operator().scales(),original);
        }
    }
}
