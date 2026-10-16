//! Force-driven realization of the EXISTING fs-material GeneralizedMaxwell card.
//!
//! A Prony solid has creep compliance J(s) = 1/E_instant + sum w_j/(1+s*tau_j).
//! We use fs-modal to compute that equivalent positive retardation spectrum.
//! Replacing its instantaneous spring by the existing unilateral WoolFelt law
//! adds rate-dependent relaxation WITHOUT a tensile-force clamp or a new felt
//! envelope. This is a series composition, not stress addition followed by max(0).
//!
//! Each retained branch, in mechanical coordinates, obeys c*z_dot+k*z=F.
//! A held force gives z1=a*z0+(1-a)*F/k. Its EXACT loss is
//! k*(1-a*a)*(z0-F/k)^2/2 >= 0. Thus F*(z1-z0)=U1-U0+loss.
//! Adding its displacement compliance to the implicit contact solve closes the
//! work of the coupled island. Memory relaxes even when no string is touched.
//! No allocations or transcendental functions occur in the stepping methods.
use fs_material::visco::GeneralizedMaxwell;
use fs_math::det;

pub const MAX_TERMS: usize = 8;
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Memory(pub [f64; MAX_TERMS]);
#[derive(Clone, Copy, Debug, Default)]
struct Mode { k: f64, decay: f64, response: f64, loss_fraction: f64 }
#[derive(Clone, Debug)]
pub struct Prepared {
    modes: [Mode; MAX_TERMS], count: usize, compliance: f64,
}
#[derive(Clone, Debug)]
pub struct Spectrum { terms: Vec<(f64, f64)> }

/// Explicitly estimated felt time scales. Instantaneous reference modulus is
/// 5 MPa, matching the virgin tangent of demonstration_law at strain 0.2.
/// These are NOT a measured Steinway hammer coupon or a claimed Stulov fit.
pub fn demonstration_prony() -> GeneralizedMaxwell {
    GeneralizedMaxwell::new(2.5e6, vec![(2.0e6, 0.0002), (0.5e6, 0.004)])
        .expect("admissible authored Prony constants")
}

impl Spectrum {
    /// Cold algebraic realization, NOT a fresh material fit. Eliminating the
    /// instantaneous strain gives H=diag(Ej)-E E^T/Einst, R=diag(Ej*tauj).
    /// R*z_dot+H*z=(E/Einst)*stress. R-orthonormal eigenvectors of (H,R)
    /// diagonalize this to positive Kelvin branches; no polynomial root finder.
    pub fn from_prony(card: &GeneralizedMaxwell) -> Result<Self,String> {
        GeneralizedMaxwell::new(card.e_inf,card.terms.clone()).map_err(|e|e.to_string())?;
        let terms:Vec<_>=card.terms.iter().copied().filter(|t|t.0>0.0).collect();
        let n=terms.len();
        if n>MAX_TERMS {return Err("felt creep supports at most eight nonzero Prony terms".into());}
        if n==0 {return Ok(Self{terms:Vec::new()});}
        let instant=card.e_inf+terms.iter().map(|t|t.0).sum::<f64>();
        if !instant.is_finite() {return Err("instantaneous felt modulus overflow".into());}
        let mut h=vec![0.0;n*n];let mut r=vec![0.0;n*n];
        for i in 0..n {
            r[i*n+i]=terms[i].0*terms[i].1;
            for j in 0..n {h[i*n+j]=if i==j {terms[i].0}else{0.0};h[i*n+j]-=terms[i].0/instant*terms[j].0;}
        }
        let eigen=fs_modal::eigh_gen_dense(&h,&r,n).map_err(|e|e.to_string())?;
        let mut creep=Vec::new();
        for m in eigen {
            if !m.lambda.is_finite()||m.lambda<=0.0||m.residual>1e-6*m.lambda {
                return Err("Prony creep spectrum not resolved positive definite".into());
            }
            let coupling=m.phi.iter().zip(&terms).map(|(v,t)|v*t.0/instant).sum::<f64>();
            let weight=coupling*coupling/m.lambda; // equilibrium compliance [1/Pa]
            // Exactly unobservable equal-time branches can produce roundoff
            // weights. Only machine-scale compliance is discarded.
            if weight*instant<=1e-14 {continue;}
            let k=1.0/weight;let tau=1.0/m.lambda;
            if !k.is_finite()||!tau.is_finite()||k<=0.0||tau<=0.0 {
                return Err("nonfinite creep realization".into());
            }
            creep.push((k,tau));
        }
        Ok(Self{terms:creep})
    }
    pub fn prepare(&self,area_m2:f64,thickness_m:f64,dt:f64)->Result<Prepared,String> {
        if [area_m2,thickness_m,dt].iter().any(|x|!x.is_finite()||*x<=0.0) {
            return Err("creep preparation needs positive SI area/thickness/dt".into());
        }
        let mut result=Prepared{modes:[Mode::default();MAX_TERMS],count:self.terms.len(),compliance:0.0};
        for (i,&(modulus,tau)) in self.terms.iter().enumerate() {
            let k=modulus*area_m2/thickness_m;
            let decay=det::exp(-dt/tau);
            let response=-det::expm1(-dt/tau)/k;
            let loss_fraction=-det::expm1(-2.0*dt/tau);
            if !k.is_finite()||k<=0.0||!response.is_finite() {
                return Err("creep mechanical coefficients overflow".into());
            }
            result.modes[i]=Mode{k,decay,response,loss_fraction};result.compliance+=response;
        }
        if !result.compliance.is_finite() {return Err("creep compliance overflow".into());}
        Ok(result)
    }
}
impl Prepared {
    pub fn compliance(&self)->f64 {self.compliance}
    pub fn deformation(&self,m:&Memory)->f64 {m.0[..self.count].iter().sum()}
    pub fn free_deformation(&self,m:&Memory)->f64 {
        self.modes[..self.count].iter().zip(&m.0).map(|(b,z)|b.decay*z).sum()
    }
    pub fn stored(&self,m:&Memory)->f64 {
        self.modes[..self.count].iter().zip(&m.0).map(|(b,z)|0.5*b.k*z*z).sum()
    }
    /// Trial only. The caller commits this memory after its coupled step passes.
    /// The force must be the same held force used in the hammer/string solve.
    pub fn advance(&self,m:&Memory,force:f64)->(Memory,f64) {
        let mut next=Memory::default();let mut loss=0.0;
        for (i,b) in self.modes[..self.count].iter().enumerate() {
            next.0[i]=b.decay*m.0[i]+b.response*force;
            loss+=0.5*b.k*b.loss_fraction*(m.0[i]-force/b.k).powi(2);
        }
        (next,loss)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn creep_realization_matches_its_existing_prony_owner() {
        let gm=demonstration_prony();let spectrum=Spectrum::from_prony(&gm).unwrap();
        let instant=gm.e_inf+gm.terms.iter().map(|t|t.0).sum::<f64>();
        for i in 0..=80 {
            let w=10.0f64.powf(-1.0+8.0*i as f64/80.0);
            let(mut jr,mut ji)=(1.0/instant,0.0);
            for &(k,tau) in &spectrum.terms {let x=w*tau;jr+=1.0/(k*(1.0+x*x));ji-=x/(k*(1.0+x*x));}
            let (er,ei)=gm.modulus(w);let den=er*er+ei*ei;
            assert!(((jr-er/den).abs()+(ji+ei/den).abs())*instant<1e-8);
        }
    }
    #[test]
    fn held_force_and_release_close_work_without_negative_loss() {
        let s=Spectrum::from_prony(&demonstration_prony()).unwrap();
        for dt in [1e-8,1.0/192_000.0,1.0/8_000.0,0.1] {
            let p=s.prepare(1e-4,0.008,dt).unwrap();let mut m=Memory::default();
            let mut work=0.0;let mut dissipated=0.0;
            for i in 0..2000 {
                let f=if i<1000 {30.0*(i as f64*0.03).sin().abs()}else{0.0};
                let (next,loss)=p.advance(&m,f);
                work+=f*(p.deformation(&next)-p.deformation(&m));dissipated+=loss;
                assert!(loss>=0.0);assert!(next.0.iter().all(|z|*z>=0.0));
                assert!((work-p.stored(&next)-dissipated).abs()<1e-9);
                m=next;
            }
        }
    }
    #[test]
    fn elastic_card_has_no_relaxation_memory_and_zero_terms_are_legal() {
        let gm=GeneralizedMaxwell::new(1e6,vec![(0.0,0.1)]).unwrap();
        let p=Spectrum::from_prony(&gm).unwrap().prepare(1.,1.,0.01).unwrap();
        assert_eq!(p.compliance(),0.0);
        assert_eq!(p.advance(&Memory::default(),3.0),(Memory::default(),0.0));
        assert!(Spectrum::from_prony(&GeneralizedMaxwell{e_inf:f64::NAN,terms:vec![]}).is_err());
    }
}
