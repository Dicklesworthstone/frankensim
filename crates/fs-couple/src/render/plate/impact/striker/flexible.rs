//! Geometry-derived planar shaft flexure in the existing impact time owner.
//! The exact rigid rotation and elastic modes share one consistent Rayleigh
//! mass pencil. Splitting their storage is only an address layout: tip/hand
//! rows still act on EVERY mode. No second time integrator or sampled attack.
use super::{ImpactError, RadiusStation, invalid};
use super::super::{BodyPotential, ImpactBody};
use crate::modal_acoustic_time::ModalAcousticState;
use fs_plate::shell::stiffened::beam::{RoundBeamModes, RoundBeamSpec, RoundStation};

/// Cold, geometry-owned pin-supported shaft basis. A pin fixes transverse
/// translation but not rotation. Planar small-deflection bending only; no
/// inferred hand impedance, wood anisotropy, shear deformation or tip law.
#[derive(Debug, Clone)]
pub struct FlexibleStriker {
    beam: RoundBeamModes,
    damping_ratio: f64,
}
/// Explicit original rigid address and appended elastic range. These refer to
/// mechanical q/p pairs, never Kelvin/material/acoustic-history scalar slots.
#[derive(Debug, Clone)]
pub struct StrikerPorts {
    rigid: usize,
    elastic_start: usize,
    omega: Vec<f64>,
    tip: Vec<f64>,
    hand: Vec<f64>,
}
/// Endpoint shaft observations; the energy is already in the joint ledger.
#[derive(Debug, Clone, Copy)]
pub struct StrikerObservation {
    pub tip_displacement_m: f64,
    pub tip_velocity_m_s: f64,
    pub hand_displacement_m: f64,
    pub hand_velocity_m_s: f64,
    pub flexural_energy_j: f64,
}
impl FlexibleStriker {
    /// Use the original tapered-beam FEM and retain ALL modes in the supplied
    /// frequency window. Damping is an explicit modal research input, not
    /// identified from a wood name. The rigid mode is never damped or tethered.
    pub fn new(profile: &[RadiusStation], spec: RoundBeamSpec, damping_ratio: f64)
        -> Result<Self, ImpactError>
    {
        if !damping_ratio.is_finite() || !(0.0..1.0).contains(&damping_ratio) {
            return Err(invalid("flexible striker needs finite modal damping in [0,1)"));
        }
        let stations: Vec<_> = profile.iter().map(|p| RoundStation {
            x_m: p.position_m, radius_m: p.radius_m,
        }).collect();
        let beam = RoundBeamModes::new(&stations, spec)
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        if beam.omega.iter().any(|w| !w.powi(2).is_finite())
            || beam.tip[0] <= 0.0 || !beam.tip[0].is_finite()
            || beam.omega.iter().any(|w| !(2.0*damping_ratio*w).is_finite()) {
            return Err(invalid("flexible striker modal coefficients overflow"));
        }
        Ok(Self { beam, damping_ratio })
    }
    pub fn omegas(&self) -> &[f64] { &self.beam.omega }
    pub fn tip_weights(&self) -> &[f64] { &self.beam.tip }
    pub fn hand_weights(&self) -> &[f64] { &self.beam.hand }
    pub fn pivot_inertia_kg_m2(&self) -> f64 { self.beam.pivot_inertia_kg_m2 }
    pub fn elastic_modes(&self) -> usize { self.beam.omega.len()-1 }

    /// Place initial displacement and velocity in the true rigid mode ONLY.
    /// The entire stick is launched; no initial bending or arbitrary elastic
    /// phase is manufactured. Returns (rigid storage, elastic storage) so an
    /// instrument can keep its original head/shell/wire/cavity addresses.
    pub fn split_bodies(&self, position_m: f64, velocity_m_s: f64)
        -> Result<(ImpactBody, ImpactBody), ImpactError>
    {
        let q = position_m/self.beam.tip[0];
        let v = velocity_m_s/self.beam.tip[0];
        if !q.is_finite() || !v.is_finite() || !(q*q+v*v).is_finite() {
            return Err(invalid("flexible striker initial rigid motion overflow"));
        }
        let rigid = ImpactBody { potential: BodyPotential::Linear(vec![0.0]),
            initial: vec![ModalAcousticState { displacement_m_sqrt_kg: q,
                velocity_m_sqrt_kg_per_s: v }], damping_per_s: vec![0.0] };
        let elastic = ImpactBody { potential: BodyPotential::Linear(self.beam.omega[1..].to_vec()),
            initial: vec![ModalAcousticState::default(); self.elastic_modes()],
            damping_per_s: self.beam.omega[1..].iter().map(|w|2.0*self.damping_ratio*w).collect() };
        Ok((rigid, elastic))
    }
    /// Bind the separated storage to explicit, nonoverlapping mechanical slots.
    /// The clock guard is the SAME one used by ImpactSystem::new.
    pub fn ports(&self, rigid: usize, elastic_start: usize, modes: usize, dt_s: f64)
        -> Result<StrikerPorts, ImpactError>
    {
        let end = elastic_start.checked_add(self.elastic_modes())
            .ok_or_else(||invalid("flexible striker layout overflow"))?;
        if modes > super::super::MAX_IMPACT_MODES || rigid >= modes || end > modes
            || (elastic_start..end).contains(&rigid) || !dt_s.is_finite() || dt_s <= 0.0
            || self.beam.omega.iter().any(|w|w*dt_s >= 0.9*std::f64::consts::PI) {
            return Err(invalid("flexible striker layout or mechanical Nyquist guard failed"));
        }
        Ok(StrikerPorts { rigid, elastic_start, omega: self.beam.omega.clone(),
            tip: self.beam.tip.clone(), hand: self.beam.hand.clone() })
    }
}
impl StrikerPorts {
    pub fn rigid_coordinate(&self) -> usize { self.rigid }
    pub fn elastic_start(&self) -> usize { self.elastic_start }
    pub fn elastic_modes(&self) -> usize { self.omega.len()-1 }
    fn index(&self, mode: usize) -> usize {
        if mode == 0 { self.rigid } else { self.elastic_start+mode-1 }
    }
    fn row(&self, local: &[f64], modes: usize) -> Result<Vec<f64>, ImpactError> {
        if modes > super::super::MAX_IMPACT_MODES || self.rigid >= modes
            || self.elastic_start+self.elastic_modes() > modes {
            return Err(invalid("flexible striker port exceeds mechanical layout"));
        }
        let mut row = vec![0.0; modes];
        for (i,&b) in local.iter().enumerate() { row[self.index(i)] = b; }
        Ok(row)
    }
    /// Physical tip motion / generalized motion, or its work-conjugate force
    /// transpose. Subtract a real surface row to form unilateral contact.
    pub fn tip_row(&self, modes: usize) -> Result<Vec<f64>, ImpactError> { self.row(&self.tip,modes) }
    /// Hand force must use the declared hand station, NOT the tip row. Signed
    /// elastic participation is retained, including arbitrary eigenvector signs.
    pub fn hand_row(&self, modes: usize) -> Result<Vec<f64>, ImpactError> { self.row(&self.hand,modes) }
    pub fn observe(&self, state: &[f64]) -> Result<StrikerObservation, ImpactError> {
        let required = 2*(self.rigid+1).max(self.elastic_start+self.elastic_modes());
        if state.len() < required { return Err(invalid("short flexible striker state")); }
        let mut out = StrikerObservation { tip_displacement_m: 0.0, tip_velocity_m_s: 0.0,
            hand_displacement_m: 0.0, hand_velocity_m_s: 0.0, flexural_energy_j: 0.0 };
        for i in 0..self.omega.len() {
            let k=2*self.index(i); let q=state[k]; let v=state[k+1];
            if !q.is_finite() || !v.is_finite() {return Err(invalid("nonfinite flexible striker state"));}
            out.tip_displacement_m += self.tip[i]*q; out.tip_velocity_m_s += self.tip[i]*v;
            out.hand_displacement_m += self.hand[i]*q; out.hand_velocity_m_s += self.hand[i]*v;
            if i>0 {out.flexural_energy_j += 0.5*((self.omega[i]*q).powi(2)+v*v);}
        }
        if [out.tip_displacement_m,out.tip_velocity_m_s,out.hand_displacement_m,
            out.hand_velocity_m_s,out.flexural_energy_j].iter().any(|v|!v.is_finite()) {
            return Err(invalid("flexible striker observation overflow"));
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::StrikerProperties;
    use super::super::super::{ImpactConfig,ImpactSystem};
    use fs_dcontact::Obstacle;
    use fs_exec::CancelGate;
    fn profile() -> [RadiusStation;2] {
        [(0.0,0.005),(0.4,0.005)].map(|(position_m,radius_m)|RadiusStation{position_m,radius_m})
    }
    fn spec() -> RoundBeamSpec { RoundBeamSpec {young_pa:12e9,density_kg_m3:800.0,
        pivot_m:0.1,contact_m:0.39,hand_m:0.16,subdivisions:8,maximum_hz:3000.0,maximum_modes:17} }
    fn shaft() -> FlexibleStriker {FlexibleStriker::new(&profile(),spec(),0.001).unwrap()}
    fn fixture() -> (ImpactSystem,StrikerPorts) {
        let s=shaft(); let n=2+s.elastic_modes(); let p=s.ports(0,2,n,2e-6).unwrap();
        let (rigid,elastic)=s.split_bodies(-0.00002,0.5).unwrap();
        let (target,_)=ImpactBody::free_mass(0.1,0.0,0.0).unwrap();
        let mut contact=p.tip_row(n).unwrap();contact[1]=-1.0/0.1_f64.sqrt();
        let ob=Obstacle::new(contact,1,n,vec![0.0],vec![1.0],1e7,1.5,"declared elastic test target".into()).unwrap();
        let system=ImpactSystem::new(vec![rigid,target,elastic],vec![ob],vec![],vec![],
            ImpactConfig{dt_s:2e-6,max_steps:512,maximum_energy_j:1.0,
                energy_absolute_tolerance_j:1e-9,energy_relative_tolerance:1e-6,
                maximum_generalized_force:1e6}).unwrap();
        (system,p)
    }
    #[test]
    fn rigid_launch_matches_original_profile_inertia_without_initial_bending() {
        let s=shaft();let old=StrikerProperties::from_profile(&profile(),800.0,0.1,0.39).unwrap();
        assert!((s.pivot_inertia_kg_m2()/old.pivot_inertia_kg_m2-1.0).abs()<1e-12);
        assert!((s.tip_weights()[0]*old.contact_effective_mass_kg.sqrt()-1.0).abs()<1e-12);
        let (rigid,elastic)=s.split_bodies(-0.0002,0.8).unwrap();
        assert_eq!(rigid.damping_per_s,[0.0]);assert!(elastic.damping_per_s.iter().all(|d|*d>0.0));
        assert!(elastic.initial.iter().all(|q|q.displacement_m_sqrt_kg==0.0 && q.velocity_m_sqrt_kg_per_s==0.0));
        let (a,p)=fixture();let o=p.observe(a.state()).unwrap();
        assert!((o.tip_velocity_m_s-0.5).abs()<1e-14);assert_eq!(o.flexural_energy_j,0.0);
    }
    #[test]
    fn contacts_excite_real_flexure_and_close_shared_energy_with_analytic_retry() {
        let (a,p)=fixture();let initial=a.stored_energy_j();let mut a=a.prepare_analytic().unwrap();
        let (mut b,_)=fixture();let force=vec![0.0;2+p.elastic_modes()];
        let gate=CancelGate::new_clock_free();let mut loss=0.0;let mut bending=0.0_f64;
        for tick in 0..256 {
            if tick==80 {
                let before=a.state().to_vec();let mut invalid=force.clone();invalid[0]=f64::NAN;
                assert!(a.step(&invalid,&gate).is_err());assert_eq!(a.state(),before);
            }
            let f=a.step(&force,&gate).unwrap();b.step(&force,&gate).unwrap();
            loss+=f.dissipated_energy_j;
            assert!((f.stored_energy_j+loss-initial).abs()<2e-7);
            bending=bending.max(p.observe(a.state()).unwrap().flexural_energy_j);
        }
        assert!(bending>1e-10,"contact must excite elastic motion, not merely add unused modes");
        assert!(a.state()[3].abs()>1e-7,"equal-and-opposite contact must accelerate the target");
        assert!((p.observe(a.state()).unwrap().tip_velocity_m_s-p.observe(b.state()).unwrap().tip_velocity_m_s).abs()<1e-4);
    }
    #[test]
    fn hand_force_is_work_conjugate_at_its_own_station_and_layout_refuses_aliases() {
        let s=shaft();let p=s.ports(3,5,5+s.elastic_modes(),2e-6).unwrap();
        let tip=p.tip_row(5+s.elastic_modes()).unwrap();let hand=p.hand_row(tip.len()).unwrap();
        assert_ne!(tip,hand);assert_eq!(&tip[..3],&[0.0;3]);assert_eq!(tip[4],0.0);
        let mut x=vec![0.0;2*tip.len()];for i in 0..tip.len(){x[2*i+1]=0.01*(i+1) as f64;}
        let o=p.observe(&x).unwrap();
        let work=hand.iter().enumerate().map(|(i,b)|2.5*b*x[2*i+1]).sum::<f64>();
        assert!((work-2.5*o.hand_velocity_m_s).abs()<1e-12);
        assert!(s.ports(3,3,30,2e-6).is_err());assert!(s.ports(0,2,2,2e-6).is_err());
        assert!(s.ports(0,2,30,1.0).is_err());assert!(p.observe(&[]).is_err());
        assert!(FlexibleStriker::new(&profile(),spec(),f64::NAN).is_err());
    }
}
