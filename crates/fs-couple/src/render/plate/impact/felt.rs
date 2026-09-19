//! Compression-only felt with the existing WoolFelt envelope/crush history.
//! The frozen-history primitive has exactly the owner's stress derivative.
//! History changes occur after a mechanical solve and remove explicit crush
//! energy. Kelvin creep branches add rate-dependent recovery in the SAME pHS
//! solve; no lagged force, invented damping curve or piano preset is selected.
use fs_material::fiber::{Uniaxial, WoolFelt, WoolFeltState};
use super::{ImpactError, invalid};

/// One series Kelvin spring/dashpot deformation, not a second time stepper.
#[derive(Debug, Clone, Copy)]
pub struct KelvinBranch {
    /// Stiffness [N/m].
    pub stiffness_n_m: f64,
    /// Dashpot viscosity [N s/m].
    pub viscosity_n_s_m: f64,
}
/// A physical compressed pad. Positive weights dot q increases compression.
/// Use opposite weights for the top and bottom washers; distributed supports
/// use separate patches whose areas sum to the actual loaded felt area.
#[derive(Debug, Clone)]
pub struct FeltPad {
    /// Loaded contact area [m^2], not the stand's outer bounding box.
    pub area_m2: f64,
    /// Uncompressed felt thickness [m].
    pub thickness_m: f64,
    /// Geometric preload compression at q=0 and zero creep [m].
    pub precompression_m: f64,
    /// Reciprocal displacement/force weights [1/sqrt(kg)], all mechanical modes.
    pub weights: Vec<f64>,
    /// Existing physical law; stand-specific parameters require identification.
    pub law: WoolFelt,
    /// Largest compression strain before this run (conditioning history).
    pub prior_maximum_strain: f64,
    /// Series viscoelastic recovery elements, at most four; initially relaxed.
    pub creep: Vec<KelvinBranch>,
}
impl FeltPad {
    pub(super) fn validate(&self, modes: usize) -> Result<(), ImpactError> {
        let l=&self.law;
        if ![self.area_m2,self.thickness_m,l.f_ref,l.eps_ref,l.p,l.q,l.crush_fraction,l.eps_densify]
            .iter().all(|v|v.is_finite() && *v>0.0)
            || !self.precompression_m.is_finite() || !self.prior_maximum_strain.is_finite()
            || self.prior_maximum_strain<0.0 || self.prior_maximum_strain>l.eps_densify
            || l.p<=1.0 || l.q<l.p || l.crush_fraction>=1.0 || l.eps_densify<=l.eps_ref
            || self.weights.len()!=modes || self.weights.iter().any(|v|!v.is_finite())
            || self.creep.len()>4 || self.creep.iter().any(|b|
                !b.stiffness_n_m.is_finite() || b.stiffness_n_m<=0.0
                || !b.viscosity_n_s_m.is_finite() || b.viscosity_n_s_m<=0.0
                || !(b.stiffness_n_m/b.viscosity_n_s_m).is_finite()) {
            return Err(invalid("felt needs explicit finite geometry, admissible WoolFelt history and passive creep"));
        }
        Ok(())
    }
    pub(super) fn recovered(&self, strain: f64, state:&WoolFeltState)->f64 {
        if state.eps_max<=0.0 {return 0.0;}
        let residual=self.law.eps_residual(state);
        let span=state.eps_max-residual;
        let x=((strain-residual)/span).max(0.0);
        self.area_m2*self.thickness_m*state.sig_max*span/(self.law.q+1.0)*x.powf(self.law.q+1.0)
    }
    pub(super) fn path_energy(&self,strain:f64,state:&WoolFeltState)->f64 {
        if strain<=state.eps_max {return self.recovered(strain,state);}
        let l=&self.law;
        // Positive difference computed without cancellation near a prior maximum.
        let power=l.p+1.0;
        let high=(strain/l.eps_ref).powf(power);
        let difference=if state.eps_max==0.0 {high}else{
            high * (-(-power*((strain-state.eps_max)/state.eps_max).ln_1p()).exp_m1())
        };
        self.recovered(state.eps_max,state)+self.area_m2*self.thickness_m*l.f_ref*l.eps_ref/power*difference
    }
    pub(super) fn force(&self,strain:f64,state:&WoolFeltState)->f64 {
        self.area_m2*self.law.stress(strain,state)
    }
    pub(super) fn history_at(&self,strain:f64)->WoolFeltState {
        self.law.update_state(strain.max(0.0),&self.law.initial_state())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn pad()->FeltPad {FeltPad{area_m2:0.001,thickness_m:0.006,precompression_m:0.0,
        weights:vec![1.0],law:WoolFelt::new(30000.0,0.2,2.2,3.0,0.15,0.7).unwrap(),
        prior_maximum_strain:0.0,creep:vec![]}}
    #[test]
    fn history_primitive_has_the_existing_material_stress_derivative() {
        let p=pad();let h=p.history_at(0.3);
        for strain in [0.01,0.08,0.2,0.299,0.301,0.45] {
            let delta=1e-7;
            let fd=(p.path_energy(strain+delta,&h)-p.path_energy(strain-delta,&h))/(2.0*delta*p.thickness_m);
            let force=p.force(strain,&h);
            assert!((fd-force).abs()<2e-7*force.abs().max(1.0));
        }
    }
    #[test]
    fn permanent_crush_dissipates_while_conditioned_subloops_are_elastic() {
        let p=pad();let virgin=p.history_at(0.0);let h=p.history_at(0.4);
        let loss=p.path_energy(0.4,&virgin)-p.recovered(0.4,&h);
        assert!(loss>0.0);
        assert_eq!(p.force(0.0,&h),0.0);
        assert_eq!(p.recovered(0.0,&h),0.0);
        assert!(p.law.eps_residual(&h)>0.0);
        let lower=0.2;let reloaded=p.law.update_state(lower,&h);
        assert_eq!(h,reloaded);
        assert_eq!(p.path_energy(lower,&h),p.recovered(lower,&reloaded));
        // Increase to a new maximum: only irreversible energy leaves storage.
        let next=p.history_at(0.5);
        assert!(p.path_energy(0.5,&h)-p.recovered(0.5,&next)>0.0);
    }
}
