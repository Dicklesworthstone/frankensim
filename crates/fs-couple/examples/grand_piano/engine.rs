//! Allocation-free prepared composition: moving-boundary modal strings,
//! one shared board, separate felt contact patches, free-flight hammers and
//! pedal-controlled dissipative ports. No samples, envelopes or reverb presets.
//!
//! note_on supplies a POST-ESCAPEMENT hammer velocity. The catch is idealized;
//! this is not a claim to reconstruct the complete Steinway grand action.
use fs_material::{Uniaxial, WoolFelt};
use super::{felt,geometry::Course,linear::{Bank,BoardMode}};

const GRAVITY: f64 = 9.80665;
const CATCH_DISTANCE: f64 = 0.020; // authored ideal backcheck, not a factory dimension

#[derive(Debug,Clone,Copy)]
pub enum Error {
    InvalidControl, UnknownKey, NotRearmed, NoConvergence,
    Contact(&'static str), NonFinite, Energy { defect_j:f64 },
}
impl std::fmt::Display for Error {
    fn fmt(&self,f:&mut std::fmt::Formatter<'_>)->std::fmt::Result {
        match self {
            Self::InvalidControl=>write!(f,"control outside its finite admitted range"),
            Self::UnknownKey=>write!(f,"key is not present in the admitted scale"),
            Self::NotRearmed=>write!(f,"hammer/key has not returned for another strike"),
            Self::NoConvergence=>write!(f,"coupled contact solve exhausted 32 sweeps"),
            Self::Contact(s)=>write!(f,"{s}"),Self::NonFinite=>write!(f,"nonfinite or excessive mechanical state"),
            Self::Energy{defect_j}=>write!(f,"mechanical work/energy defect {defect_j:e} J"),
        }
    }
}
impl std::error::Error for Error {}

#[derive(Clone,Copy,Debug)]
struct Hammer { y:f64,v:f64,active:bool,held:bool,latched:bool }
impl Default for Hammer {
    fn default()->Self {Self{y:-CATCH_DISTANCE,v:0.0,active:false,held:false,latched:false}}
}
#[derive(Clone,Debug)]
struct Contact { state:felt::State,overlap:f64,force:f64,enabled:bool }

#[derive(Clone,Copy,Debug,Default)]
pub struct Accounting {
    pub input_work_j:f64,
    pub felt_loss_j:f64,
    pub modal_loss_j:f64,
    pub damper_loss_j:f64,
    pub catch_loss_j:f64,
    pub max_balance_error_j:f64,
}
impl Accounting {
    pub fn dissipated_j(self)->f64 {self.felt_loss_j+self.modal_loss_j+self.damper_loss_j+self.catch_loss_j}
}

pub struct Instrument {
    pub bank:Bank,
    courses:Vec<Course>,law:WoolFelt,hammers:Vec<Hammer>,contacts:Vec<Contact>,
    output_rate:u32,substeps:usize,sustain:f64,sostenuto:bool,una_corda:bool,
    /// Explicitly authored upper damper break; not a verified Steinway D value.
    pub last_damped_midi:u8,
    pub damper_drag_ns_m:f64,
    pub accounting:Accounting,
    contact_h:Vec<f64>,force:Vec<f64>,gap:Vec<f64>,active:Vec<usize>,
    hammer_free:Vec<f64>,hammer_next:Vec<Hammer>,
    saved_q:Vec<f64>,saved_v:Vec<f64>,saved_hammers:Vec<Hammer>,saved_contacts:Vec<Contact>,
}

impl Instrument {
    pub fn new(courses:Vec<Course>,board:&[BoardMode],rate:u32,substeps:usize,
        modes_per_string:usize,damping:bool)->Result<Self,String> {
        if !(8_000..=192_000).contains(&rate)||!(1..=16).contains(&substeps){return Err("invalid rate/substep budget".into());}
        let mechanics_rate=rate.checked_mul(substeps as u32).ok_or("mechanics rate overflow")?;
        let bank=Bank::new(&courses,board,mechanics_rate,0.45*f64::from(rate),modes_per_string,damping)?;
        let law=felt::demonstration_law()?;
        let nc=bank.contact_strings.len();let dt=1.0/f64::from(mechanics_rate);
        let hammers=vec![Hammer::default();courses.len()];
        let contacts=vec![Contact{state:law.initial_state(),overlap:-CATCH_DISTANCE,force:0.0,enabled:true};nc];
        let mut contact_h=bank.contact_compliance.clone();
        for i in 0..nc {for j in 0..nc {
            let ci=bank.strings[bank.contact_strings[i]].course;
            let cj=bank.strings[bank.contact_strings[j]].course;
            if ci==cj {contact_h[i*nc+j]+=0.5*dt*dt/courses[ci].hammer_mass_kg;}
        }}
        Ok(Self {saved_q:bank.q.clone(),saved_v:bank.v.clone(),saved_hammers:hammers.clone(),
            saved_contacts:contacts.clone(),hammer_next:hammers.clone(),hammer_free:vec![0.0;courses.len()],
            bank,courses,law,hammers,contacts,output_rate:rate,substeps,sustain:0.0,
            sostenuto:false,una_corda:false,last_damped_midi:88,damper_drag_ns_m:0.4,
            accounting:Accounting::default(),contact_h,force:vec![0.0;nc],gap:vec![0.0;nc],
            active:Vec::with_capacity(nc)})
    }

    pub fn sample_rate(&self)->u32{self.output_rate}
    pub fn set_sustain(&mut self,value:f64)->Result<(),Error>{
        if !value.is_finite()||!(0.0..=1.0).contains(&value){return Err(Error::InvalidControl);}
        self.sustain=value;Ok(())
    }
    pub fn set_una_corda(&mut self,on:bool){self.una_corda=on;}
    pub fn set_sostenuto(&mut self,on:bool){
        if on&&!self.sostenuto {for h in &mut self.hammers {h.latched=h.held;}}
        if !on {for h in &mut self.hammers {h.latched=false;}}
        self.sostenuto=on;
    }
    fn key_index(&self,midi:u8)->Result<usize,Error>{self.courses.iter().position(|c|c.midi==midi).ok_or(Error::UnknownKey)}

    /// Velocity is SI m/s, not an arbitrary MIDI brightness parameter. All
    /// force/spectral changes emerge from the contact law and coupled motion.
    pub fn note_on(&mut self,midi:u8,hammer_velocity_m_s:f64)->Result<(),Error>{
        if !hammer_velocity_m_s.is_finite()||!(0.0..=8.0).contains(&hammer_velocity_m_s)||hammer_velocity_m_s==0.0 {
            return Err(Error::InvalidControl);
        }
        let ci=self.key_index(midi)?;
        if self.hammers[ci].active||self.hammers[ci].held{return Err(Error::NotRearmed);}
        let before=self.energy_j();let mut launch:f64=-0.002;let mut member=0;
        for c in 0..self.contacts.len(){
            let s=&self.bank.strings[self.bank.contact_strings[c]];
            if s.course!=ci {continue;}
            // Same area per physical string; una corda does not duplicate the
            // remaining patch's force or replace the instrument's material.
            let count=if self.una_corda{self.courses[ci].unison.saturating_sub(1).max(1)}else{self.courses[ci].unison};
            self.contacts[c].enabled=member<count;member+=1;
            if self.contacts[c].enabled {
                let x=self.bank.contact_position(c,&self.bank.q);
                let free=self.law.eps_residual(&self.contacts[c].state)*self.courses[ci].felt_thickness_m;
                launch=launch.min(x+free-0.002);
            }
        }
        self.hammers[ci].y=launch;self.hammers[ci].v=hammer_velocity_m_s;
        self.hammers[ci].held=true;self.hammers[ci].active=true;
        for c in 0..self.contacts.len(){if self.bank.strings[self.bank.contact_strings[c]].course==ci {
            self.contacts[c].overlap=launch-self.bank.contact_position(c,&self.bank.q);self.contacts[c].force=0.0;
        }}
        self.accounting.input_work_j+=self.energy_j()-before;Ok(())
    }
    pub fn note_off(&mut self,midi:u8)->Result<(),Error>{let i=self.key_index(midi)?;self.hammers[i].held=false;Ok(())}

    pub fn energy_j(&self)->f64 {
        let mut e=self.bank.energy();
        for (h,c) in self.hammers.iter().zip(&self.courses){
            e+=0.5*c.hammer_mass_kg*h.v*h.v+c.hammer_mass_kg*GRAVITY*(h.y+CATCH_DISTANCE);
        }
        for (i,p) in self.contacts.iter().enumerate(){if p.enabled {
            let c=&self.courses[self.bank.strings[self.bank.contact_strings[i]].course];
            e+=c.felt_area_m2/c.unison as f64*c.felt_thickness_m
                *felt::stored(&self.law,p.overlap/c.felt_thickness_m,&p.state);
        }}
        e
    }

    fn damp(&mut self,dt:f64)->Result<f64,Error>{
        if !self.damper_drag_ns_m.is_finite()||self.damper_drag_ns_m<0.0{return Err(Error::InvalidControl);}
        let mut loss=0.0;
        for si in 0..self.bank.strings.len(){
            let ci=self.bank.strings[si].course;let h=self.hammers[ci];
            if self.bank.strings[si].contact.is_some()&&!h.held&&!h.latched&&self.courses[ci].midi<=self.last_damped_midi {
                let drag=self.damper_drag_ns_m*(1.0-self.sustain).powi(2);
                loss+=self.bank.damp_string(si,drag,dt);
            }
        }
        Ok(loss)
    }

    fn mechanics_step(&mut self)->Result<(),Error>{
        let dt=1.0/f64::from(self.bank.rate);let nc=self.contacts.len();
        let before=self.energy_j();let mut damper_loss=self.damp(0.5*dt)?;
        self.bank.predict();self.active.clear();
        for (i,h) in self.hammers.iter().enumerate(){
            self.hammer_free[i]=if h.active{h.y+dt*h.v-0.5*dt*dt*GRAVITY}else{h.y};
        }
        self.force.fill(0.0);
        for c in 0..nc {
            let ci=self.bank.strings[self.bank.contact_strings[c]].course;
            self.gap[c]=self.hammer_free[ci]-self.bank.free_contact[c];
            if self.hammers[ci].active&&self.contacts[c].enabled {
                self.active.push(c);self.force[c]=self.contacts[c].force;
            }
        }
        for &i in &self.active {for &j in &self.active{self.gap[i]-=self.contact_h[i*nc+j]*self.force[j];}}
        let mut converged=self.active.is_empty();
        for _ in 0..32 {
            if converged{break;}
            for index in 0..self.active.len(){
                let i=self.active[index];let ci=self.bank.strings[self.bank.contact_strings[i]].course;
                let c=self.courses[ci];let diagonal=self.contact_h[i*nc+i];
                let free=self.gap[i]+diagonal*self.force[i];
                let next=felt::solve(&self.law,&self.contacts[i].state,self.contacts[i].overlap,
                    free,diagonal,c.felt_thickness_m,c.felt_area_m2/c.unison as f64).map_err(Error::Contact)?;
                let change=next-self.force[i];self.force[i]=next;
                for &j in &self.active{self.gap[j]-=self.contact_h[j*nc+i]*change;}
            }
            converged=true;
            for &i in &self.active {
                let c=self.courses[self.bank.strings[self.bank.contact_strings[i]].course];
                let expected=felt::average(&self.law,&self.contacts[i].state,self.contacts[i].overlap,
                    self.gap[i],c.felt_thickness_m,c.felt_area_m2/c.unison as f64).0;
                if !expected.is_finite()||(self.force[i]-expected).abs()>1e-5+1e-8*expected.abs(){converged=false;}
            }
        }
        if !converged{return Err(Error::NoConvergence);}
        self.bank.finish(&self.force);
        if self.bank.next_q.iter().chain(&self.bank.next_v).any(|x|!x.is_finite()||x.abs()>1e5){return Err(Error::NonFinite);}
        self.hammer_next.copy_from_slice(&self.hammers);
        for (i,h) in self.hammer_next.iter_mut().enumerate(){if h.active {
            let force=(0..nc).filter(|&c|self.bank.strings[self.bank.contact_strings[c]].course==i).map(|c|self.force[c]).sum::<f64>();
            let accel=GRAVITY+force/self.courses[i].hammer_mass_kg;
            h.y=self.hammers[i].y+dt*self.hammers[i].v-0.5*dt*dt*accel;
            h.v=self.hammers[i].v-dt*accel;
        }}
        let mut felt_loss=0.0;
        for i in 0..nc {
            let ci=self.bank.strings[self.bank.contact_strings[i]].course;let c=self.courses[ci];
            let end=self.hammer_next[ci].y-self.bank.contact_position(i,&self.bank.next_q);
            if self.contacts[i].enabled {
                if end/c.felt_thickness_m>self.law.eps_densify+1e-10{return Err(Error::Contact("felt densification bound exceeded"));}
                let old=&self.contacts[i];let state=self.law.update_state(end/c.felt_thickness_m,&old.state);
                let volume=c.felt_area_m2/c.unison as f64*c.felt_thickness_m;
                let delta=volume*(felt::stored(&self.law,end/c.felt_thickness_m,&state)
                    -felt::stored(&self.law,old.overlap/c.felt_thickness_m,&old.state));
                felt_loss+=self.force[i]*(end-old.overlap)-delta;
                self.contacts[i].state=state;
            }
            self.contacts[i].overlap=end;self.contacts[i].force=self.force[i];
        }
        self.hammers.copy_from_slice(&self.hammer_next);self.bank.commit();
        let mut catch_loss=0.0;
        for i in 0..self.hammers.len(){
            let h=self.hammers[i];
            if h.active&&h.y< -CATCH_DISTANCE&&h.v<0.0 {
                let clear=(0..nc).filter(|&c|self.bank.strings[self.bank.contact_strings[c]].course==i)
                    .all(|c|self.contacts[c].overlap<=0.0);
                if clear {
                    let m=self.courses[i].hammer_mass_kg;
                    catch_loss+=0.5*m*h.v*h.v+m*GRAVITY*(h.y+CATCH_DISTANCE);
                    self.hammers[i].y=-CATCH_DISTANCE;self.hammers[i].v=0.0;self.hammers[i].active=false;
                    for c in 0..nc {if self.bank.strings[self.bank.contact_strings[c]].course==i {
                        self.contacts[c].overlap=-CATCH_DISTANCE-self.bank.contact_position(c,&self.bank.q);self.contacts[c].force=0.0;
                    }}
                }
            }
        }
        damper_loss+=self.damp(0.5*dt)?;
        let modal_loss=self.bank.last_modal_loss_j;
        let after=self.energy_j();let balance=after-before+felt_loss+modal_loss+damper_loss+catch_loss;
        let tolerance=1e-9+1e-8*before.abs().max(after.abs());
        if !after.is_finite()||after>100.0{return Err(Error::NonFinite);}
        if balance.abs()>tolerance||felt_loss< -tolerance||modal_loss< -tolerance||catch_loss< -tolerance {
            return Err(Error::Energy{defect_j:balance});
        }
        self.accounting.felt_loss_j+=felt_loss;self.accounting.modal_loss_j+=modal_loss;
        self.accounting.damper_loss_j+=damper_loss;self.accounting.catch_loss_j+=catch_loss;
        self.accounting.max_balance_error_j=self.accounting.max_balance_error_j.max(balance.abs());
        Ok(())
    }

    /// Transactional audio sample. Every scratch/history buffer is allocated
    /// at construction; a refused substep restores the ENTIRE sample, including
    /// prior substeps, felt maxima, hammers and all component loss accounts.
    /// Returns a surface volume-velocity diagnostic, NOT calibrated pressure.
    pub fn step(&mut self)->Result<f64,Error>{
        self.saved_q.copy_from_slice(&self.bank.q);self.saved_v.copy_from_slice(&self.bank.v);
        self.saved_hammers.copy_from_slice(&self.hammers);self.saved_contacts.clone_from_slice(&self.contacts);
        let saved=self.accounting;let modal=self.bank.last_modal_loss_j;
        let mut average=0.0;
        for _ in 0..self.substeps {
            if let Err(error)=self.mechanics_step(){
                self.bank.q.copy_from_slice(&self.saved_q);self.bank.v.copy_from_slice(&self.saved_v);
                self.hammers.copy_from_slice(&self.saved_hammers);self.contacts.clone_from_slice(&self.saved_contacts);
                self.accounting=saved;self.bank.last_modal_loss_j=modal;return Err(error);
            }
            average+=self.bank.volume_velocity();
        }
        // Box averaging is a diagnostic decimator, not a certified audio
        // antialias filter. The renderer must disclose this approximation.
        Ok(average/self.substeps as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn instrument()->Instrument{
        let scale=super::super::geometry::demonstration_scale().unwrap();
        Instrument::new(vec![scale[48]],&super::super::board::demonstration(),48_000,4,12,true).unwrap()
    }
    #[test]
    fn real_contact_is_passive_audible_and_replayable(){
        let mut a=instrument();let mut b=instrument();
        a.note_on(69,3.0).unwrap();b.note_on(69,3.0).unwrap();
        let mut peak:f64=0.0;
        for _ in 0..4800{let x=a.step().unwrap();let y=b.step().unwrap();assert_eq!(x.to_bits(),y.to_bits());peak=peak.max(x.abs());}
        assert!(peak>1e-10);assert!(a.accounting.felt_loss_j>0.0);
        let balance=a.accounting.input_work_j-a.accounting.dissipated_j()-a.energy_j();
        assert!(balance.abs()<1e-7,"{balance:e}");
    }
    #[test]
    fn invalid_controls_and_failed_samples_preserve_state(){
        let mut p=instrument();assert!(p.note_on(69,f64::NAN).is_err());assert!(p.set_sustain(1.1).is_err());
        p.note_on(69,2.0).unwrap();let before=p.energy_j();let q=p.bank.q.clone();
        p.damper_drag_ns_m=f64::NAN;assert!(p.step().is_err());assert_eq!(p.energy_j(),before);assert_eq!(p.bank.q,q);
    }
    #[test]
    fn sostenuto_captures_only_keys_held_on_its_rising_edge(){
        let mut p=instrument();p.set_sostenuto(true);p.note_on(69,1.0).unwrap();
        assert!(!p.hammers[0].latched);p.set_sostenuto(false);p.set_sostenuto(true);
        assert!(p.hammers[0].latched);p.note_off(69).unwrap();assert!(p.hammers[0].latched);
        p.set_sostenuto(false);assert!(!p.hammers[0].latched);
    }
}
