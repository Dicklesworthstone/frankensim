//! Shared vacancy balance for backward-Euler competitive Langmuir kinetics.
//! A diagonal-plus-rank-one solve costs O(species), not O(species^3).
use super::*;

/// Allocation-free one-species limit and its analytic rate derivative.
pub mod scalar;
pub use scalar::{ScalarSiteError, ScalarSiteStep, backward_euler_site};
/// Allocation-free common-site evaluation for spatial transport consumers.
mod prepared;
pub use prepared::{CompetitiveSiteBuffer, backward_euler_sites_into};

#[derive(Clone, Debug, PartialEq)]
pub struct CompetitiveSiteStep {
    pub initial_moles: Vec<f64>,
    pub final_moles: Vec<f64>,
    pub gross_adsorption_mol: Vec<f64>,
    pub gross_desorption_mol: Vec<f64>,
    pub vacancy_moles: [f64; 2],
    pub maximum_balance_residual_mol: f64,
}

/// All components consume ONE vacancy inventory. Rates are per vacant/occupied
/// site [1/s], respectively; gross streams use the implicit endpoint rates.
/// This is first-order backward Euler with frozen physical forcing, not an
/// exact time integral or a finite-bath equilibrium solve. Zero rates are inert.
pub fn backward_euler_sites(
    initial: &[f64], capacity: f64, adsorption: &[f64], desorption: &[f64],
    duration: f64, mut cancelled: impl FnMut() -> bool,
) -> Result<CompetitiveSiteStep> {
    let mut work = CompetitiveSiteBuffer::default();
    backward_euler_sites_into(initial, capacity, adsorption, desorption, duration,
        &mut work, &mut cancelled)?;
    Ok(CompetitiveSiteStep { initial_moles: initial.to_vec(),
        final_moles: work.final_moles().to_vec(),
        gross_adsorption_mol: work.gross_adsorption_mol().to_vec(),
        gross_desorption_mol: work.gross_desorption_mol().to_vec(),
        vacancy_moles: work.vacancy_moles(),
        maximum_balance_residual_mol: work.maximum_balance_residual_mol() })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn near(a: f64, b: f64) { assert!((a-b).abs() <= 2e-13 * (1.0+a.abs()+b.abs()), "{a} != {b}"); }
    #[test]
    fn shared_vacancies_follow_independent_two_by_two_inverse() {
        let s = backward_euler_sites(&[0.2,0.3], 1.0, &[4.0,2.0], &[1.0,3.0], 0.1, || false).unwrap();
        // [1.5,.4; .2,1.5] n1 = [.6,.5]; determinant=2.17.
        near(s.final_moles[0], (0.6*1.5-0.4*0.5)/2.17);
        near(s.final_moles[1], (1.5*0.5-0.2*0.6)/2.17);
        near(s.vacancy_moles[1], 1.0-s.final_moles.iter().sum::<f64>());
        for i in 0..2 { near(s.final_moles[i]-s.initial_moles[i],s.gross_adsorption_mol[i]-s.gross_desorption_mol[i]); }
    }
    #[test]
    fn equilibrium_is_shared_and_gross_turnover_does_not_disappear() {
        let s = backward_euler_sites(&[0.5,0.25],1.0,&[2.0,4.0],&[1.0,4.0],0.25,||false).unwrap();
        near(s.final_moles[0],0.5); near(s.final_moles[1],0.25);
        for i in 0..2 { assert!(s.gross_adsorption_mol[i]>0.0); near(s.gross_adsorption_mol[i],s.gross_desorption_mol[i]); }
    }
    #[test]
    fn competition_displaces_a_weak_adsorbate_without_overfilling() {
        let s = backward_euler_sites(&[0.8,0.0],1.0,&[1.0,100.0],&[1.0,1.0],1.0,||false).unwrap();
        assert!(s.final_moles[0]<0.8 && s.final_moles[1]>0.0);
        assert!(s.final_moles.iter().sum::<f64>()<1.0);
        let full = backward_euler_sites(&[1.0,0.0],1.0,&[0.0,10.0],&[1.0,1.0],0.1,||false).unwrap();
        assert!(full.vacancy_moles[1]>0.0 && full.final_moles[1]>0.0);
    }
    #[test]
    fn scalar_limit_refines_at_first_order_toward_exact_relaxation() {
        let exact = 0.75 + (0.2-0.75)*(-4.0_f64).exp();
        let run = |count: usize| {
            let mut n = 0.2;
            for _ in 0..count { n=backward_euler_sites(&[n],1.0,&[3.0],&[1.0],1.0/count as f64,||false).unwrap().final_moles[0]; }
            (n-exact).abs()
        };
        let ratio=run(128)/run(256); assert!(ratio>1.9 && ratio<2.1);
    }
    #[test]
    fn zero_rates_closed_limits_order_and_cancellation() {
        let s=backward_euler_sites(&[0.25,0.75],1.0,&[0.0,0.0],&[0.0,0.0],1e20,||false).unwrap();
        assert_eq!(s.initial_moles,s.final_moles);
        let a=backward_euler_sites(&[0.2,0.3],1.0,&[4.0,2.0],&[1.0,3.0],0.1,||false).unwrap();
        let b=backward_euler_sites(&[0.3,0.2],1.0,&[2.0,4.0],&[3.0,1.0],0.1,||false).unwrap();
        near(a.final_moles[0],b.final_moles[1]);
        assert!(backward_euler_sites(&[0.8,0.8],1.0,&[1.0,1.0],&[1.0,1.0],0.1,||false).is_err());
        assert!(matches!(backward_euler_sites(&[0.1],1.0,&[1.0],&[1.0],0.1,||true),Err(AdsorptionError::Cancelled)));
    }
}
