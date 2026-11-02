//! Allocation-free scalar limit of the shared-vacancy backward-Euler site law.
//!
//! `dn/dt = a (C-n) - d n`. Rates are independently nonnegative; a Langmuir
//! consumer supplies `a = k_d b p`, while other callers retain their original
//! rate convention. The derivative is with respect to `a`, holding `d` fixed.

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScalarSiteError(&'static str);
impl ScalarSiteError {
    pub const fn message(self) -> &'static str { self.0 }
}
impl core::fmt::Display for ScalarSiteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result { f.write_str(self.0) }
}
impl std::error::Error for ScalarSiteError {}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ScalarSiteStep {
    pub occupied_moles: f64,
    pub vacancy_moles: f64,
    pub occupied_change_mol: f64,
    pub gross_adsorption_mol: f64,
    pub gross_desorption_mol: f64,
    /// d(n_new)/d(a) [mol s], needed when pressure and storage are coupled.
    pub adsorption_rate_derivative_mol_s: f64,
    pub maximum_balance_residual_mol: f64,
}

pub fn backward_euler_site(initial: f64, capacity: f64, adsorption_rate: f64,
    desorption_rate: f64, duration_s: f64) -> Result<ScalarSiteStep, ScalarSiteError> {
    let nonnegative = |x: f64| x.is_finite() && x >= 0.0;
    if !capacity.is_finite() || capacity <= 0.0 || !nonnegative(initial) || initial > capacity
        || !nonnegative(adsorption_rate) || !nonnegative(desorption_rate) || !nonnegative(duration_s) {
        return Err(ScalarSiteError("finite scalar site capacity, bounded loading and nonnegative rates/time required"));
    }
    let qa = duration_s * adsorption_rate;
    let qd = duration_s * desorption_rate;
    if !nonnegative(qa) || !nonnegative(qd)
        || (duration_s > 0.0 && adsorption_rate > 0.0 && qa == 0.0)
        || (duration_s > 0.0 && desorption_rate > 0.0 && qd == 0.0) {
        return Err(ScalarSiteError("scalar integrated site rate is unrepresentable"));
    }
    // Scale before summing: 1+qa+qd may overflow even when all fractions fit.
    let scale = 1.0_f64.max(qa).max(qd);
    let denominator = 1.0 / scale + qa / scale + qd / scale;
    let inverse = (1.0 / scale) / denominator;
    let wa = (qa / scale) / denominator;
    let wd = (qd / scale) / denominator;
    let vacant0 = capacity - initial;
    let change = wa * vacant0 - wd * initial;
    let occupied = if change.abs() < 0.5 * initial { initial + change }
        else { inverse * initial + wa * capacity };
    // An empty population with no integrated adsorption remains exactly empty.
    // Evaluating the two complementary weights separately can otherwise round
    // vacancy one ULP above capacity. This analytical limit also keeps the
    // nonzero derivative with respect to adsorption at zero loading.
    let vacancy = if initial == 0.0 && qa == 0.0 { capacity }
        else { inverse * vacant0 + wd * capacity };
    let incoming = qa * vacancy;
    let outgoing = qd * occupied;
    let derivative = (duration_s * inverse) * vacancy;
    if ![occupied, vacancy, incoming, outgoing, derivative].into_iter().all(nonnegative)
        || occupied > capacity || vacancy > capacity
        || ((initial > 0.0 || qa > 0.0) && occupied == 0.0)
        || ((vacant0 > 0.0 || qd > 0.0) && vacancy == 0.0)
        || (qa > 0.0 && vacancy > 0.0 && incoming == 0.0)
        || (qd > 0.0 && occupied > 0.0 && outgoing == 0.0) {
        return Err(ScalarSiteError("scalar occupied/vacant sites or gross streams are unrepresentable"));
    }
    let change = occupied - initial;
    let residual = (change - (incoming - outgoing)).abs()
        .max((occupied + vacancy - capacity).abs());
    if !residual.is_finite() || residual > 128.0 * f64::EPSILON * capacity {
        return Err(ScalarSiteError("original scalar site equations failed their roundoff balance"));
    }
    Ok(ScalarSiteStep { occupied_moles: occupied, vacancy_moles: vacancy,
        occupied_change_mol: change, gross_adsorption_mol: incoming,
        gross_desorption_mol: outgoing, adsorption_rate_derivative_mol_s: derivative,
        maximum_balance_residual_mol: residual })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn scalar_site_amount_derivative_and_gross_streams_match_the_exact_rational_law() {
        let (n0,c,a,d,dt)=(0.2,1.3,2.5,0.7,0.4);
        let step=backward_euler_site(n0,c,a,d,dt).unwrap();
        let denominator=1.0+dt*(a+d);
        let exact=(n0+dt*a*c)/denominator;
        let derivative=dt*((1.0+dt*d)*c-n0)/(denominator*denominator);
        assert!((step.occupied_moles-exact).abs()<1e-15);
        assert!((step.adsorption_rate_derivative_mol_s-derivative).abs()<1e-15);
        assert!((step.gross_adsorption_mol-dt*a*(c-exact)).abs()<1e-15);
        assert!((step.gross_desorption_mol-dt*d*exact).abs()<1e-15);
    }
    #[test]
    fn inert_equilibrium_and_empty_limits_keep_finite_sites_exactly() {
        let inert=backward_euler_site(0.75,1.0,0.0,0.0,1e200).unwrap();
        assert_eq!(inert.occupied_moles,0.75);
        let equilibrium=backward_euler_site(0.5,1.0,2.0,2.0,1e8).unwrap();
        assert_eq!(equilibrium.occupied_moles,0.5);
        assert_eq!(equilibrium.occupied_change_mol,0.0);
        assert_eq!(equilibrium.gross_adsorption_mol,equilibrium.gross_desorption_mol);
        assert_eq!(backward_euler_site(0.0,1.0,0.0,2.0,1.0).unwrap().occupied_moles,0.0);
        assert!(backward_euler_site(1.1,1.0,1.0,1.0,1.0).is_err());
    }
    #[test]
    fn empty_sites_preserve_inventory_at_small_steps_and_retain_the_uptake_jacobian() {
        let capacity = 1e-6;
        let desorption = 0.5;
        for duration in [0.0625, 0.03125, 1e-12, 1e6] {
            let step = backward_euler_site(0.0,capacity,0.0,desorption,duration).unwrap();
            assert_eq!(step.occupied_moles,0.0);
            assert_eq!(step.vacancy_moles,capacity);
            assert_eq!(step.occupied_change_mol,0.0);
            assert_eq!(step.gross_adsorption_mol,0.0);
            assert_eq!(step.gross_desorption_mol,0.0);
            assert_eq!(step.maximum_balance_residual_mol,0.0);
            let derivative = duration/(1.0+duration*desorption)*capacity;
            assert!(step.adsorption_rate_derivative_mol_s>0.0);
            assert!((step.adsorption_rate_derivative_mol_s-derivative).abs()
                <=4.0*f64::EPSILON*derivative);
        }
        // A small positive adsorption perturbation must create the amount
        // predicted by the zero-loading Jacobian used in the coupled Newton solve.
        let empty = backward_euler_site(0.0,capacity,0.0,desorption,0.0625).unwrap();
        let perturbed = backward_euler_site(0.0,capacity,1e-6,desorption,0.0625).unwrap();
        assert!((perturbed.occupied_moles/1e-6-empty.adsorption_rate_derivative_mol_s).abs()
            <=1e-7*empty.adsorption_rate_derivative_mol_s);
        assert_eq!(perturbed.occupied_moles+perturbed.vacancy_moles,capacity);
    }
}
