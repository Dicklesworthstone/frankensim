//! Two-coordinate, geometry-derived hammer/shank reduction, not a point mass.
//!
//! Chabassier & Durufle, JSV 333 (2014), Tables 1-3 and Fig. 2:
//! https://www.math.u-bordeaux.fr/~durufle/article/HammerJSVV2.pdf
//! The published L,A,I,rho,E,G,kappa,H and jack station are used below.
//! This is OUR linearized Ritz image, NOT their full nonlinear rotating beam.
//! One rigid rotation and one shear-corrected bending shape are retained.
//! Gravity is linearized about the horizontal shank. The head follows the
//! paper's offset point-mass kinematics; crown geometry, head rotary inertia,
//! roller-felt compliance and the rest of the action linkage are not inferred.
//!
//! All modes come from fs-modal. The existing exact-ZOH owner advances bending
//! through its allocation-free workspace. Contact and jack use reciprocal
//! force/displacement projections. No oscillator or constitutive is duplicated.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeModel, ModalAcousticWorkspace};
use fs_math::{c64::C64, det};

const GRAVITY: f64 = 9.80665;

#[derive(Clone, Copy, Debug)]
pub struct Geometry {
    pub length_m: f64, pub area_m2: f64, pub inertia_m4: f64,
    pub density_kg_m3: f64, pub young_pa: f64, pub shear_pa: f64,
    pub shear_factor: f64, pub head_offset_m: f64, pub jack_station_m: f64,
    /// An explicit estimate; the cited table does not provide shank damping.
    pub damping_ratio: f64,
}
impl Geometry {
    pub fn published() -> Self {
        Self { length_m:0.086, area_m2:32.38e-6, inertia_m4:83.44e-12,
            density_kg_m3:560.0, young_pa:10.18e9, shear_pa:0.64e9,
            shear_factor:0.85, head_offset_m:0.04, jack_station_m:0.0155,
            damping_ratio:0.01 }
    }
    fn validate(self) -> Result<(),String> {
        if [self.length_m,self.area_m2,self.inertia_m4,self.density_kg_m3,
            self.young_pa,self.shear_pa,self.shear_factor,self.jack_station_m]
            .iter().any(|x|!x.is_finite()||*x<=0.0)
            || self.jack_station_m>=self.length_m || !self.head_offset_m.is_finite()
            || self.head_offset_m<0.0 || !self.damping_ratio.is_finite()
            || !(0.0..=1.0).contains(&self.damping_ratio) {
            return Err("invalid SI shank geometry/material or jack station".into());
        }
        Ok(())
    }
}

/// Authoritative state lives here, not in the prepared owner's scratch model.
/// q[1] is measured from static bending equilibrium under gravity.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct State { pub q:[f64;2], pub v:[f64;2] }

pub struct Prepared {
    dt:f64, head:[f64;2], jack:[f64;2], gravity_rigid:f64, omega2:f64,
    compliance:f64, rest_deflection:f64,
    oscillator:Option<(ModalAcousticTimeModel,ModalAcousticWorkspace)>,
}
impl Prepared {
    /// Retain the old point-mass image exactly for direct coupon comparisons.
    pub fn point_mass(mass:f64,rate:u32) -> Result<Self,String> {
        if !mass.is_finite()||mass<=0.0||rate<8_000 {return Err("invalid hammer mass/rate".into());}
        let dt=1.0/f64::from(rate);let root=det::sqrt(mass);
        Ok(Self {dt,head:[1.0/root,0.0],jack:[0.0;2],gravity_rigid:-GRAVITY*root,
            omega2:0.0,compliance:0.5*dt*dt/mass,rest_deflection:0.0,oscillator:None})
    }

    pub fn from_geometry(g:Geometry,head_mass:f64,rate:u32) -> Result<Self,String> {
        g.validate()?;
        if !head_mass.is_finite()||head_mass<=0.0||rate<8_000 {return Err("invalid hammer mass/rate".into());}
        let l=g.length_m;let mr=g.density_kg_m3*g.area_m2*l;
        let rotary=g.density_kg_m3*g.inertia_m4/l;
        let bending=l.powi(3)/(3.0*g.young_pa*g.inertia_m4);
        let shear=l/(g.shear_factor*g.shear_pa*g.area_m2);
        let a=bending/(bending+shear);let b=1.0-a;
        // Coordinates u=L*theta, q=tip bending. The static Timoshenko shape
        // is psi(t)=a*(3t^2-t^3)/2+b*t, section rotation=a*(3t-1.5t^2)/L.
        // These analytic integrals include distributed translational AND rotary
        // inertia. The offset head contributes H^2 only to rigid rotation,
        // consistently with xi=L*ur+(w-H)*utheta in the cited linearization.
        let m00=head_mass*(1.0+(g.head_offset_m/l).powi(2))+mr/3.0+rotary;
        let m01=head_mass+mr*(a*11.0/40.0+b/3.0)+rotary*a;
        let m11=head_mass+mr*(a*a*33.0/140.0+2.0*a*b*11.0/40.0+b*b/3.0)
            +rotary*a*a*6.0/5.0;
        let stiffness=1.0/(bending+shear);
        let modes=fs_modal::eigh_gen_dense(&[0.0,0.0,0.0,stiffness],
            &[m00,m01,m01,m11],2).map_err(|e|e.to_string())?;
        if modes.len()!=2 || modes[1].lambda<=0.0 || !modes[1].lambda.is_finite()
            || modes[0].lambda.abs()>1e-8*modes[1].lambda {
            return Err("shank must resolve one rigid and one positive bending mode".into());
        }
        let mut head=[0.0;2];let mut jack=[0.0;2];let mut gravity=[0.0;2];
        let t=g.jack_station_m/l;let psi=a*(1.5*t*t-0.5*t*t*t)+b*t;
        for (i,m) in modes.iter().enumerate() {
            if m.residual>1e-7*modes[1].lambda {return Err("unresolved shank mode".into());}
            let sign=if (if i==0 {m.phi[0]}else{m.phi[1]})<0.0 {-1.0}else{1.0};
            head[i]=sign*(m.phi[0]+m.phi[1]);
            jack[i]=sign*(t*m.phi[0]+psi*m.phi[1]);
            gravity[i]=-GRAVITY*sign*((head_mass+mr/2.0)*m.phi[0]
                +(head_mass+mr*(a*3.0/8.0+b/2.0))*m.phi[1]);
        }
        if head[0]<=0.0||gravity[0]>=0.0 {return Err("invalid rigid shank projection".into());}
        let omega2=modes[1].lambda;
        let mut model=ModalAcousticTimeModel::try_new(rate,vec![ModalAcousticMode {
            angular_frequency_rad_s:det::sqrt(omega2), damping_ratio:g.damping_ratio,
            pressure_per_modal_velocity:C64::new(0.0,0.0),
        }],ModalAcousticTimeBudget::audible_reference()).map_err(|e|e.to_string())?;
        let mut workspace=ModalAcousticWorkspace::new(&model);
        // One cold impulse-response probe of the existing owner, not another
        // ZOH formula. The gravity equilibrium is absorbed into the state origin.
        model.step_into(&[1.0],&mut workspace).map_err(|e|e.to_string())?;
        let bq=model.states()[0].displacement_m_sqrt_kg;
        let dt=1.0/f64::from(rate);
        let compliance=0.5*dt*dt*head[0]*head[0]+head[1]*head[1]*bq;
        if !compliance.is_finite()||compliance<=0.0 {return Err("invalid shank contact compliance".into());}
        let rest_deflection=-gravity[0]*head[1]/(head[0]*omega2);
        Ok(Self {dt,head,jack,gravity_rigid:gravity[0],omega2,compliance,rest_deflection,
            oscillator:Some((model,workspace))})
    }
    pub fn is_flexible(&self)->bool {self.oscillator.is_some()}
    pub fn compliance(&self)->f64 {self.compliance}
    pub fn position(&self,s:&State)->f64 {self.head[0]*s.q[0]+self.head[1]*s.q[1]}
    pub fn velocity(&self,s:&State)->f64 {self.head[0]*s.v[0]+self.head[1]*s.v[1]}
    pub fn jack_position(&self,s:&State)->f64 {self.jack[0]*s.q[0]+self.jack[1]*s.q[1]}
    pub fn launch(&self,y:f64,v:f64)->State {State{q:[y/self.head[0],0.0],v:[v/self.head[0],0.0]}}
    pub fn rest(&self,distance:f64)->State {
        State{q:[(-distance-self.head[1]*self.rest_deflection)/self.head[0],self.rest_deflection],v:[0.0;2]}
    }
    pub fn energy(&self,s:&State,rest_distance:f64)->f64 {
        0.5*(s.v[0]*s.v[0]+s.v[1]*s.v[1]+self.omega2*s.q[1]*s.q[1])
            -self.gravity_rigid*(s.q[0]+rest_distance/self.head[0])
            +0.5*self.omega2*self.rest_deflection*self.rest_deflection
    }
    /// Trial: restore the caller's accepted state on EVERY call. The caller may
    /// freely discard this result or roll back the entire audio sample.
    pub fn advance(&mut self,s:&State,head_force:f64,jack_force:f64)->Result<(State,f64),&'static str> {
        let acceleration=self.gravity_rigid+self.head[0]*head_force+self.jack[0]*jack_force;
        let mut next=*s;let mut loss=0.0;
        next.q[0]+=self.dt*s.v[0]+0.5*self.dt*self.dt*acceleration;
        next.v[0]+=self.dt*acceleration;
        if let Some((model,workspace))=&mut self.oscillator {
            model.restore_states(&[ModalAcousticState {
                displacement_m_sqrt_kg:s.q[1],velocity_m_sqrt_kg_per_s:s.v[1],
            }]).map_err(|_|"shank modal checkpoint refused")?;
            let frame=model.step_into(&[self.head[1]*head_force+self.jack[1]*jack_force],workspace)
                .map_err(|_|"shank modal transition refused")?;
            loss=frame.viscous_dissipation_j;
            next.q[1]=model.states()[0].displacement_m_sqrt_kg;
            next.v[1]=model.states()[0].velocity_m_sqrt_kg_per_s;
        }
        if next.q.iter().chain(&next.v).any(|x|!x.is_finite()) {return Err("nonfinite shank state");}
        Ok((next,loss))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn point_mass_image_recovers_existing_flight_and_energy() {
        let mut p=Prepared::point_mass(0.01,48_000).unwrap();let s=p.launch(-0.002,2.0);
        let (next,loss)=p.advance(&s,-10.0,0.0).unwrap();let dt=1.0/48_000.0;
        assert!((p.position(&next)-(-0.002+2.0*dt-0.5*(1000.0+GRAVITY)*dt*dt)).abs()<1e-15);
        assert!((p.velocity(&next)-(2.0-(1000.0+GRAVITY)*dt)).abs()<1e-12);
        assert_eq!(loss,0.0);
    }
    #[test]
    fn geometry_changes_bending_not_an_authored_frequency() {
        let g=Geometry::published();let p=Prepared::from_geometry(g,0.008,192_000).unwrap();
        let q=Prepared::from_geometry(Geometry{young_pa:2.0*g.young_pa,shear_pa:2.0*g.shear_pa,..g},0.008,192_000).unwrap();
        assert!((q.omega2/p.omega2-2.0).abs()<1e-8);
        assert!((200.0..400.0).contains(&(det::sqrt(p.omega2)/std::f64::consts::TAU)));
    }
    #[test]
    fn reciprocal_head_and_jack_ports_close_work_and_superposition() {
        let mut p=Prepared::from_geometry(Geometry::published(),0.008,192_000).unwrap();
        let mut s=p.launch(-0.015,0.0);let mut work=0.0;let mut loss=0.0;
        let e0=p.energy(&s,0.02);
        for i in 0..3000 {
            let jack=if i<1000 {30.0}else{0.0};
            let force=if i>1500&&i<1800 {-10.0}else{0.0};
            let (free,_)=p.advance(&s,0.0,jack).unwrap();
            let (next,d)=p.advance(&s,force,jack).unwrap();
            assert!((p.position(&next)-p.position(&free)-force*p.compliance()).abs()<1e-12);
            work+=force*(p.position(&next)-p.position(&s))+jack*(p.jack_position(&next)-p.jack_position(&s));
            loss+=d;
            assert!((p.energy(&next,0.02)-e0+loss-work).abs()<1e-8);
            s=next;
        }
        assert!(loss>0.0);
    }
    #[test]
    fn invalid_geometry_is_not_silently_repaired() {
        let g=Geometry::published();
        for bad in [Geometry{area_m2:0.0,..g},Geometry{young_pa:f64::NAN,..g},
            Geometry{jack_station_m:g.length_m,..g}] {
            assert!(Prepared::from_geometry(bad,0.008,48_000).is_err());
        }
    }
}
