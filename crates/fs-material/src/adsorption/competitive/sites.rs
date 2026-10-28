//! Shared vacancy balance for backward-Euler competitive Langmuir kinetics.
//! A diagonal-plus-rank-one solve costs O(species), not O(species^3).
use super::*;

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
    poll(&mut cancelled)?;
    let n = initial.len();
    if n == 0 || n > MAX_COMPETITIVE_SPECIES || adsorption.len() != n
        || desorption.len() != n || !positive(capacity) || !nonnegative(duration)
        || initial.iter().chain(adsorption).chain(desorption).any(|&x| !nonnegative(x)) {
        return Err(invalid("bounded matched competitive inventories and finite nonnegative rates required"));
    }
    let vacancy0 = finite(capacity - sum(initial.iter().copied())?)?;
    if vacancy0 < 0.0 { return Err(invalid("competitive loadings exceed their ONE shared site capacity")); }
    let mut result = CompetitiveSiteStep { initial_moles: initial.to_vec(),
        final_moles: initial.to_vec(), gross_adsorption_mol: vec![0.0; n],
        gross_desorption_mol: vec![0.0; n], vacancy_moles: [vacancy0; 2],
        maximum_balance_residual_mol: 0.0 };
    if duration == 0.0 { poll(&mut cancelled)?; return Ok(result); }
    let mut qa = [0.0; MAX_COMPETITIVE_SPECIES];
    let mut qd = [0.0; MAX_COMPETITIVE_SPECIES];
    let mut retained = [0.0; MAX_COMPETITIVE_SPECIES];
    let mut released = [0.0; MAX_COMPETITIVE_SPECIES];
    for i in 0..n {
        poll(&mut cancelled)?;
        qa[i] = finite(duration * adsorption[i])?;
        qd[i] = finite(duration * desorption[i])?;
        if (adsorption[i] > 0.0 && qa[i] == 0.0) || (desorption[i] > 0.0 && qd[i] == 0.0) {
            return Err(invalid("nonzero competitive kinetic interval underflow"));
        }
        let denominator = finite(1.0 + qd[i])?;
        retained[i] = qa[i] / denominator;
        released[i] = initial[i] * (qd[i] / denominator);
        if (qa[i] > 0.0 && retained[i] == 0.0)
            || (initial[i] > 0.0 && qd[i] > 0.0 && released[i] == 0.0) {
            return Err(invalid("competitive vacancy solve lost a nonzero contribution"));
        }
    }
    // Derive vacancy from old vacancies PLUS desorbed old sites, avoiding
    // capacity - sum(theta/(1+dt*d)) cancellation near complete saturation.
    let numerator = sum(std::iter::once(vacancy0).chain(released[..n].iter().copied()))?;
    let denominator = finite(1.0 + sum(retained[..n].iter().copied())?)?;
    let vacancy1 = numerator / denominator;
    if !nonnegative(vacancy1) || (numerator > 0.0 && vacancy1 == 0.0) {
        return Err(invalid("competitive vacancy inventory is unrepresentable"));
    }
    for i in 0..n {
        poll(&mut cancelled)?;
        let incoming = finite(qa[i] * vacancy1)?;
        // Rational positive form for large unloading; incremental form for
        // tiny changes. Neither clips occupancy nor renormalizes a species.
        let delta = finite((incoming - qd[i] * initial[i]) / (1.0 + qd[i]))?;
        let next = if delta.abs() < 0.5 * initial[i] { initial[i] + delta }
            else { (initial[i] + incoming) / (1.0 + qd[i]) };
        let outgoing = finite(qd[i] * next)?;
        if !nonnegative(next) || next > capacity
            || ((initial[i] > 0.0 || incoming > 0.0) && next == 0.0)
            || (qa[i] > 0.0 && vacancy1 > 0.0 && incoming == 0.0)
            || (qd[i] > 0.0 && next > 0.0 && outgoing == 0.0) {
            return Err(invalid("competitive occupied sites or gross streams are unrepresentable"));
        }
        result.final_moles[i] = next;
        result.gross_adsorption_mol[i] = incoming;
        result.gross_desorption_mol[i] = outgoing;
        let residual = finite((next - initial[i]) - (incoming - outgoing))?;
        result.maximum_balance_residual_mol = result.maximum_balance_residual_mol.max(residual.abs());
    }
    let residual = finite(sum(result.final_moles.iter().copied())? + vacancy1 - capacity)?;
    result.maximum_balance_residual_mol = result.maximum_balance_residual_mol.max(residual.abs());
    let roundoff = 128.0 * f64::EPSILON * capacity;
    if result.maximum_balance_residual_mol > roundoff {
        return Err(invalid("original competitive site equations failed their roundoff balance"));
    }
    result.vacancy_moles[1] = vacancy1;
    poll(&mut cancelled)?;
    Ok(result)
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
