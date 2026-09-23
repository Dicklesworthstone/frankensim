//! Prepared reciprocal strings/board, hysteretic/viscoelastic felt and hammers.
//! Both the point-mass coupon image and geometry-derived flexible shanks remain
//! available. note_on supplies a post-escapement velocity; jack_on instead drives
//! the shank through its physical jack port and disengages at let-off.
//! The backcheck/rest stop is idealized, not a full grand-action reconstruction.
use fs_material::{Uniaxial, WoolFelt};
use fs_material::visco::GeneralizedMaxwell;
use super::{felt,geometry::Course,linear::{Bank,BoardMode,dampers,hammer_footprint}};
#[path = "felt_relaxation.rs"]
mod relaxation;
#[path = "hammer_shank.rs"]
mod shank;
pub use shank::Geometry as ShankGeometry;
#[path = "radiation.rs"]
pub mod radiation;
#[path = "hammer_contact_engine.rs"]
mod contact_solver;

const CATCH_DISTANCE: f64 = 0.020; // idealized backcheck/rest gap
const LET_OFF_DISTANCE: f64 = 0.0015; // Chabassier/Durufle JSV 2014 Table 3

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

#[derive(Clone,Copy,Debug,Default)]
struct Jack { peak_n:f64, duration_s:f64, elapsed_s:f64 }
#[derive(Clone,Copy,Debug)]
struct Hammer { motion:shank::State, active:bool, held:bool, latched:bool, on_rest:bool, jack:Jack }
#[derive(Clone,Debug)]
struct Contact { state:felt::State,memory:relaxation::Memory,overlap:f64,force:f64,enabled:bool }

#[derive(Clone,Copy,Debug,Default)]
pub struct Accounting {
    pub input_work_j:f64,
    /// Total felt dissipation: permanent crush plus reversible-branch viscosity.
    pub felt_loss_j:f64,
    /// The viscous subset of felt_loss_j; NOT added a second time to the total.
    pub felt_relaxation_loss_j:f64,
    pub shank_loss_j:f64,
    pub modal_loss_j:f64,
    /// Dissipation in the passive acoustic realization, not wood/felt loss.
    pub radiation_loss_j:f64,
    pub damper_loss_j:f64,
    pub catch_loss_j:f64,
    pub max_balance_error_j:f64,
}
impl Accounting {
    pub fn dissipated_j(self)->f64 {self.felt_loss_j+self.shank_loss_j+self.modal_loss_j+self.radiation_loss_j+self.damper_loss_j+self.catch_loss_j}
}

pub struct Instrument {
    pub bank:Bank,
    courses:Vec<Course>,laws:Vec<WoolFelt>,hammers:Vec<Hammer>,contacts:Vec<Contact>,
    hammer_models:Vec<shank::Prepared>,
    creep:Vec<relaxation::Prepared>,
    // Actual per-site areas: unison allocation times longitudinal quadrature.
    contact_areas:Vec<f64>,
    contact_solver:Option<contact_solver::Prepared>,
    spatial_dampers:Option<dampers::Prepared>,
    radiation:Option<radiation::Prepared>,
    output_rate:u32,substeps:usize,sustain:f64,sostenuto:bool,una_corda:bool,
    /// Point-image controls only. A spatial specification owns its pad/free
    /// break and individual drag values instead. Neither is verified Model D data.
    pub last_damped_midi:u8,
    pub damper_drag_ns_m:f64,
    pub accounting:Accounting,
    contact_h:Vec<f64>,force:Vec<f64>,gap:Vec<f64>,active:Vec<usize>,
    hammer_free:Vec<f64>,hammer_next:Vec<Hammer>,jack_force:Vec<f64>,rest_force:Vec<f64>,
    saved_q:Vec<f64>,saved_v:Vec<f64>,saved_hammers:Vec<Hammer>,saved_contacts:Vec<Contact>,
}

impl Instrument {
    pub fn new(courses:Vec<Course>,board:&[BoardMode],rate:u32,substeps:usize,
        modes_per_string:usize,damping:bool)->Result<Self,String> {
        Self::new_with_felt(courses,board,rate,substeps,modes_per_string,damping,
            felt::demonstration_law()?,&relaxation::demonstration_prony())
    }

    /// Uniform-material, point-mass front door retained for coupon comparisons.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_felt(courses:Vec<Course>,board:&[BoardMode],rate:u32,substeps:usize,
        modes_per_string:usize,damping:bool,law:WoolFelt,prony:&GeneralizedMaxwell)->Result<Self,String> {
        let materials=(0..courses.len()).map(|_|(law.clone(),GeneralizedMaxwell {
            e_inf:prony.e_inf,terms:prony.terms.clone(),
        })).collect();
        Self::new_with_course_felts(courses,board,rate,substeps,modes_per_string,damping,materials)
    }

    /// One material card per course, with the original point-mass hammer image.
    /// Unison members share a material, never their internal contact histories.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_course_felts(courses:Vec<Course>,board:&[BoardMode],rate:u32,substeps:usize,
        modes_per_string:usize,damping:bool,materials:Vec<(WoolFelt,GeneralizedMaxwell)>)->Result<Self,String> {
        Self::build(courses,board,rate,substeps,modes_per_string,damping,materials,None)
    }

    /// Geometry-derived rigid rotation plus Timoshenko static bending reduction.
    /// Head masses and felt cards are per-course; the shank geometry is supplied.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_course_shanks(courses:Vec<Course>,board:&[BoardMode],rate:u32,substeps:usize,
        modes_per_string:usize,damping:bool,materials:Vec<(WoolFelt,GeneralizedMaxwell)>,
        geometry:ShankGeometry)->Result<Self,String> {
        Self::build(courses,board,rate,substeps,modes_per_string,damping,materials,Some(geometry))
    }

    #[allow(clippy::too_many_arguments)]
    fn build(courses:Vec<Course>,board:&[BoardMode],rate:u32,substeps:usize,
        modes_per_string:usize,damping:bool,materials:Vec<(WoolFelt,GeneralizedMaxwell)>,
        shank_geometry:Option<ShankGeometry>)->Result<Self,String> {
        Self::new_with_contact_geometry(courses,board,rate,substeps,modes_per_string,damping,
            materials,shank_geometry,None)
    }

    /// Demonstration felt/Prony material with explicit longitudinal contact
    /// geometry. This changes neither the scale's total felt area nor its mass.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_footprints(courses:Vec<Course>,board:&[BoardMode],rate:u32,substeps:usize,
        modes_per_string:usize,damping:bool,footprints:&hammer_footprint::Specification)->Result<Self,String> {
        let law=felt::demonstration_law()?;let prony=relaxation::demonstration_prony();
        let materials=(0..courses.len()).map(|_|(law.clone(),GeneralizedMaxwell {
            e_inf:prony.e_inf,terms:prony.terms.clone(),
        })).collect();
        Self::new_with_contact_geometry(courses,board,rate,substeps,modes_per_string,damping,
            materials,None,Some(footprints))
    }

    /// Complete cold physical construction. Every footprint site has separate
    /// existing felt/crush/Prony history, but all sites of one key share its ONE
    /// original hammer inertia and shank. No active material/geometry replacement.
    /// None retains the original contact image and arithmetic.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_contact_geometry(courses:Vec<Course>,board:&[BoardMode],rate:u32,substeps:usize,
        modes_per_string:usize,damping:bool,materials:Vec<(WoolFelt,GeneralizedMaxwell)>,
        shank_geometry:Option<ShankGeometry>,footprints:Option<&hammer_footprint::Specification>)->Result<Self,String> {
        if !(8_000..=192_000).contains(&rate)||!(1..=16).contains(&substeps){return Err("invalid rate/substep budget".into());}
        if materials.len()!=courses.len() {return Err("one felt/Prony card is required for every course".into());}
        let mut laws=Vec::with_capacity(materials.len());
        let mut spectra=Vec::with_capacity(materials.len());
        for (law,prony) in materials {
            if [law.f_ref,law.eps_ref,law.p,law.q,law.crush_fraction,law.eps_densify].iter().any(|x|!x.is_finite()) {
                return Err("nonfinite felt material card".into());
            }
            WoolFelt::new(law.f_ref,law.eps_ref,law.p,law.q,law.crush_fraction,law.eps_densify).map_err(|e|e.to_string())?;
            spectra.push(relaxation::Spectrum::from_prony(&prony)?);laws.push(law);
        }
        let mechanics_rate=rate.checked_mul(substeps as u32).ok_or("mechanics rate overflow")?;
        let bank=match footprints {
            Some(spec)=>Bank::new_with_hammer_footprints(&courses,board,mechanics_rate,
                0.45*f64::from(rate),modes_per_string,damping,spec)?,
            None=>Bank::new(&courses,board,mechanics_rate,0.45*f64::from(rate),modes_per_string,damping)?,
        };
        let contact_solver=contact_solver::Prepared::new(&bank,&courses)?;
        let nc=bank.contact_strings.len();let dt=1.0/f64::from(mechanics_rate);
        let contact_areas:Vec<f64>=bank.contact_strings.iter().enumerate().map(|(i,&si)| {
            let c=&courses[bank.strings[si].course];
            (c.felt_area_m2/c.unison as f64)*bank.contact_area_fraction(i)
        }).collect();
        if contact_areas.iter().any(|a|!a.is_finite()||*a<=0.0) {
            return Err("hammer contact area is not finite positive".into());
        }
        let creep=bank.contact_strings.iter().enumerate().map(|(i,&si)| {
            let ci=bank.strings[si].course;let c=&courses[ci];
            spectra[ci].prepare(contact_areas[i],c.felt_thickness_m,dt)
        }).collect::<Result<Vec<_>,_>>()?;
        let hammer_models=courses.iter().map(|c| match shank_geometry {
            Some(g)=>shank::Prepared::from_geometry(g,c.hammer_mass_kg,mechanics_rate),
            None=>shank::Prepared::point_mass(c.hammer_mass_kg,mechanics_rate),
        }).collect::<Result<Vec<_>,_>>()?;
        let hammers:Vec<_>=hammer_models.iter().map(|p|Hammer {
            motion:p.rest(CATCH_DISTANCE),active:false,held:false,latched:false,on_rest:true,jack:Jack::default(),
        }).collect();
        let contacts:Vec<_>=bank.contact_strings.iter().map(|&si| Contact {
            state:laws[bank.strings[si].course].initial_state(),memory:relaxation::Memory::default(),
            overlap:-CATCH_DISTANCE,force:0.0,enabled:true,
        }).collect();
        let mut contact_h=bank.contact_compliance.clone();
        for i in 0..nc {for j in 0..nc {
            let ci=bank.strings[bank.contact_strings[i]].course;
            let cj=bank.strings[bank.contact_strings[j]].course;
            if ci==cj {contact_h[i*nc+j]+=hammer_models[ci].compliance();}
        }}
        Ok(Self {saved_q:bank.q.clone(),saved_v:bank.v.clone(),saved_hammers:hammers.clone(),
            saved_contacts:contacts.clone(),hammer_next:hammers.clone(),hammer_free:vec![0.0;courses.len()],
            jack_force:vec![0.0;courses.len()],rest_force:vec![0.0;courses.len()],
            bank,courses,laws,hammers,contacts,hammer_models,creep,contact_areas,contact_solver,spatial_dampers:None,radiation:None,output_rate:rate,substeps,sustain:0.0,
            sostenuto:false,una_corda:false,last_damped_midi:88,damper_drag_ns_m:0.4,
            accounting:Accounting::default(),contact_h,force:vec![0.0;nc],gap:vec![0.0;nc],
            active:Vec::with_capacity(nc)})
    }

    pub fn sample_rate(&self)->u32{self.output_rate}
    /// Attach before any excitation. Rows must already be in this bank's
    /// complete mass-loaded basis. No state-reset/replacement while playing.
    pub fn configure_radiation(&mut self,model:&radiation::Model)->Result<(),String>{
        if self.radiation.is_some() || self.accounting.input_work_j!=0.
            || self.bank.q.iter().chain(&self.bank.v).any(|v|*v!=0.)
            || self.hammers.iter().any(|h|h.active||h.held) {
            return Err("radiation must be prepared once, before piano excitation".into());
        }
        let prepared=radiation::Prepared::new(model,self.bank.rate,self.bank.board_count)?;
        self.radiation=Some(prepared);Ok(())
    }
    pub fn radiation_energy_j(&self)->f64 {
        self.radiation.as_ref().map_or(0.,radiation::Prepared::energy)
    }
    pub fn has_radiation(&self)->bool {self.radiation.is_some()}
    /// Number of independent hammer/felt sites, not number of strings or voices.
    pub fn hammer_contact_count(&self)->usize{self.contacts.len()}
    /// Cold preparation in the existing loaded basis. Publish only on complete
    /// admission; configuring a viscous law neither stores energy nor resets
    /// ongoing string, board, hammer or felt history. This is not a hot control.
    pub fn configure_dampers(&mut self,spec:&dampers::Specification)->Result<(),String>{
        let prepared=dampers::Prepared::new(spec,&self.courses,&self.bank)?;
        self.spatial_dampers=Some(prepared);Ok(())
    }
    /// None is the original point image; Some counts actual retained string
    /// pads and quadrature stations, including a valid all-free specification.
    pub fn damper_resolution(&self)->Option<(usize,usize)>{
        self.spatial_dampers.as_ref().map(|d|(d.string_count(),d.cell_count()))
    }
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
    fn arm(&mut self,ci:usize)->f64 {
        let mut launch:f64=-0.002;let mut member=0;let mut previous_string=None;
        for c in 0..self.contacts.len(){
            let si=self.bank.contact_strings[c];
            if self.bank.strings[si].course!=ci {continue;}
            // Una corda selects physical strings, not a fraction of one face's
            // quadrature sites. All sites on the same string move together.
            if previous_string.is_some_and(|previous|previous!=si){member+=1;}
            previous_string=Some(si);
            let count=if self.una_corda{self.courses[ci].unison.saturating_sub(1).max(1)}else{self.courses[ci].unison};
            self.contacts[c].enabled=member<count;
            if self.contacts[c].enabled {
                let x=self.bank.contact_position(c,&self.bank.q);
                let free=self.laws[ci].eps_residual(&self.contacts[c].state)*self.courses[ci].felt_thickness_m
                    +self.creep[c].deformation(&self.contacts[c].memory);
                launch=launch.min(x+free-0.002);
            }
        }
        launch
    }
    fn reset_overlaps(&mut self,ci:usize) {
        let y=self.hammer_models[ci].position(&self.hammers[ci].motion);
        for c in 0..self.contacts.len(){if self.bank.strings[self.bank.contact_strings[c]].course==ci {
            self.contacts[c].overlap=y-self.bank.contact_position(c,&self.bank.q);self.contacts[c].force=0.0;
        }}
    }

    /// Post-escapement SI crown velocity. No preset spectrum or gain envelope.
    pub fn note_on(&mut self,midi:u8,hammer_velocity_m_s:f64)->Result<(),Error>{
        if !hammer_velocity_m_s.is_finite()||hammer_velocity_m_s<=0.0||hammer_velocity_m_s>8.0 {
            return Err(Error::InvalidControl);
        }
        let ci=self.key_index(midi)?;
        if self.hammers[ci].active||self.hammers[ci].held{return Err(Error::NotRearmed);}
        let before=self.energy_j();let launch=self.arm(ci);
        self.hammers[ci].motion=self.hammer_models[ci].launch(launch,hammer_velocity_m_s);
        self.hammers[ci].held=true;self.hammers[ci].active=true;self.hammers[ci].on_rest=false;
        self.hammers[ci].jack=Jack::default();self.reset_overlaps(ci);
        self.accounting.input_work_j+=self.energy_j()-before;Ok(())
    }

    /// Push the published jack station with F=A*sin(pi*t/duration)^2, starting
    /// from rest. Table 3 examples: (70 N,7 ms) and (30 N,100 ms). Force ceases
    /// at the earlier of pulse end or 1.5 mm let-off, resolved to one mechanical
    /// substep. This is a force-driven action fragment, not a key-velocity alias.
    pub fn jack_on(&mut self,midi:u8,peak_n:f64,duration_s:f64)->Result<(),Error> {
        if !peak_n.is_finite()||!(0.0..=200.0).contains(&peak_n)||peak_n==0.0
            || !duration_s.is_finite()||!(0.001..=0.2).contains(&duration_s) {return Err(Error::InvalidControl);}
        let ci=self.key_index(midi)?;
        if !self.hammer_models[ci].is_flexible() {return Err(Error::Contact("jack drive needs a geometry-derived shank"));}
        if self.hammers[ci].active||self.hammers[ci].held {return Err(Error::NotRearmed);}
        let before=self.energy_j();self.arm(ci);
        self.hammers[ci].motion=self.hammer_models[ci].rest(CATCH_DISTANCE);
        self.hammers[ci].held=true;self.hammers[ci].active=true;self.hammers[ci].on_rest=true;
        self.hammers[ci].jack=Jack{peak_n,duration_s,elapsed_s:0.0};self.reset_overlaps(ci);
        self.accounting.input_work_j+=self.energy_j()-before;Ok(())
    }
    pub fn note_off(&mut self,midi:u8)->Result<(),Error>{
        let i=self.key_index(midi)?;self.hammers[i].held=false;self.hammers[i].jack.peak_n=0.0;Ok(())
    }

    pub fn energy_j(&self)->f64 {
        let mut e=self.bank.energy()+self.radiation_energy_j();
        for (h,p) in self.hammers.iter().zip(&self.hammer_models) {e+=p.energy(&h.motion,CATCH_DISTANCE);}
        for (i,p) in self.contacts.iter().enumerate(){
            e+=self.creep[i].stored(&p.memory);
            if p.enabled {
                let ci=self.bank.strings[self.bank.contact_strings[i]].course;let c=&self.courses[ci];
                let elastic=p.overlap-self.creep[i].deformation(&p.memory);
                e+=self.contact_areas[i]*c.felt_thickness_m
                    *felt::stored(&self.laws[ci],elastic/c.felt_thickness_m,&p.state);
            }
        }
        e
    }

    fn damp(&mut self,dt:f64)->Result<f64,Error>{
        if let Some(dampers)=&self.spatial_dampers {
            let hammers=&self.hammers;
            return dampers.apply(&mut self.bank.v,dt,self.sustain,
                |ci|hammers[ci].held||hammers[ci].latched).map_err(Error::Contact);
        }
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
        let before=self.energy_j();
        let mut radiation_loss=if let Some(air)=&mut self.radiation {
            air.before(&mut self.bank.v[self.bank.modes.len()..])
        } else {0.};
        let mut damper_loss=self.damp(0.5*dt)?;
        self.bank.predict();self.active.clear();
        self.jack_force.fill(0.0);self.rest_force.fill(0.0);
        for i in 0..self.hammers.len() {
            let h=&mut self.hammers[i];let p=&mut self.hammer_models[i];
            let y=p.position(&h.motion);self.hammer_free[i]=y;
            if !h.active {continue;}
            if y>=-LET_OFF_DISTANCE||h.jack.elapsed_s>=h.jack.duration_s {h.jack.peak_n=0.0;}
            if h.jack.peak_n>0.0 {
                let phase=std::f64::consts::PI*(h.jack.elapsed_s+0.5*dt).min(h.jack.duration_s)/h.jack.duration_s;
                self.jack_force[i]=h.jack.peak_n*fs_math::det::sin(phase).powi(2);
                h.jack.elapsed_s+=dt;
            }
            let (free,_)=p.advance(&h.motion,0.0,self.jack_force[i]).map_err(Error::Contact)?;
            self.hammer_free[i]=p.position(&free);
            if h.on_rest {
                if self.hammer_free[i]< -CATCH_DISTANCE {
                    // Work-conjugate unilateral support: the crown stays on
                    // the rest, so its reaction does zero work, not a clamp.
                    self.rest_force[i]=(-CATCH_DISTANCE-self.hammer_free[i])/p.compliance();
                    self.hammer_free[i]=-CATCH_DISTANCE;
                } else {h.on_rest=false;}
            }
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
            self.contact_sweep()?;
            converged=true;
            for &i in &self.active {
                let ci=self.bank.strings[self.bank.contact_strings[i]].course;let c=self.courses[ci];
                let material=&self.creep[i];let old=&self.contacts[i];
                let start=old.overlap-material.deformation(&old.memory);
                let end=self.gap[i]-material.free_deformation(&old.memory)-material.compliance()*self.force[i];
                let expected=felt::average(&self.laws[ci],&old.state,start,end,
                    c.felt_thickness_m,self.contact_areas[i]).0;
                if !expected.is_finite()||(self.force[i]-expected).abs()>1e-5+1e-8*expected.abs(){converged=false;}
            }
        }
        if !converged{return Err(Error::NoConvergence);}
        self.bank.finish(&self.force);
        if self.bank.next_q.iter().chain(&self.bank.next_v).any(|x|!x.is_finite()||x.abs()>1e5){return Err(Error::NonFinite);}
        self.hammer_next.copy_from_slice(&self.hammers);
        let mut shank_loss=0.0;let mut jack_work=0.0;
        for (i,h) in self.hammer_next.iter_mut().enumerate(){if h.active {
            let force=(0..nc).filter(|&c|self.bank.strings[self.bank.contact_strings[c]].course==i).map(|c|self.force[c]).sum::<f64>();
            let model=&mut self.hammer_models[i];
            let (motion,loss)=model.advance(&self.hammers[i].motion,self.rest_force[i]-force,
                self.jack_force[i]).map_err(Error::Contact)?;
            jack_work+=self.jack_force[i]*(model.jack_position(&motion)-model.jack_position(&self.hammers[i].motion));
            shank_loss+=loss;h.motion=motion;
        }}
        let mut felt_loss=0.0;let mut relaxation_loss=0.0;
        for i in 0..nc {
            let ci=self.bank.strings[self.bank.contact_strings[i]].course;let c=self.courses[ci];
            let end=self.hammer_models[ci].position(&self.hammer_next[ci].motion)-self.bank.contact_position(i,&self.bank.next_q);
            let material=&self.creep[i];let old=&self.contacts[i];
            let (memory,viscous)=material.advance(&old.memory,self.force[i]);
            relaxation_loss+=viscous;felt_loss+=viscous;
            if old.enabled {
                let start_elastic=old.overlap-material.deformation(&old.memory);
                let end_elastic=end-material.deformation(&memory);
                if end/c.felt_thickness_m>self.laws[ci].eps_densify+1e-10 {
                    return Err(Error::Contact("total felt densification bound exceeded"));
                }
                let state=self.laws[ci].update_state(end_elastic/c.felt_thickness_m,&old.state);
                let volume=self.contact_areas[i]*c.felt_thickness_m;
                let delta=volume*(felt::stored(&self.laws[ci],end_elastic/c.felt_thickness_m,&state)
                    -felt::stored(&self.laws[ci],start_elastic/c.felt_thickness_m,&old.state));
                felt_loss+=self.force[i]*(end_elastic-start_elastic)-delta;
                self.contacts[i].state=state;
            }
            self.contacts[i].memory=memory;
            self.contacts[i].overlap=end;self.contacts[i].force=self.force[i];
        }
        self.hammers.copy_from_slice(&self.hammer_next);self.bank.commit();
        let mut catch_loss=0.0;
        for i in 0..self.hammers.len(){
            let h=self.hammers[i];let model=&self.hammer_models[i];
            let landed=model.position(&h.motion)< -CATCH_DISTANCE-1e-10&&model.velocity(&h.motion)<0.0;
            if h.active&&(landed||(h.on_rest&&h.jack.peak_n==0.0)) {
                let clear=(0..nc).filter(|&c|self.bank.strings[self.bank.contact_strings[c]].course==i)
                    .all(|c|self.contacts[c].overlap<=0.0);
                if clear {
                    let rest=model.rest(CATCH_DISTANCE);
                    catch_loss+=model.energy(&h.motion,CATCH_DISTANCE)-model.energy(&rest,CATCH_DISTANCE);
                    self.hammers[i].motion=rest;self.hammers[i].on_rest=true;
                    self.hammers[i].active=h.jack.peak_n>0.0;self.reset_overlaps(i);
                }
            }
        }
        damper_loss+=self.damp(0.5*dt)?;
        if let Some(air)=&mut self.radiation {
            radiation_loss+=air.after(&mut self.bank.v[self.bank.modes.len()..]);
        }
        let modal_loss=self.bank.last_modal_loss_j;
        let after=self.energy_j();let balance=after-before+felt_loss+shank_loss+modal_loss+radiation_loss+damper_loss+catch_loss-jack_work;
        let tolerance=1e-9+1e-8*before.abs().max(after.abs());
        if !after.is_finite()||after>100.0{return Err(Error::NonFinite);}
        if balance.abs()>tolerance||felt_loss< -tolerance||shank_loss< -tolerance||modal_loss< -tolerance||radiation_loss< -tolerance||catch_loss< -tolerance {
            return Err(Error::Energy{defect_j:balance});
        }
        self.accounting.input_work_j+=jack_work;self.accounting.shank_loss_j+=shank_loss;
        self.accounting.felt_loss_j+=felt_loss;self.accounting.felt_relaxation_loss_j+=relaxation_loss;
        self.accounting.modal_loss_j+=modal_loss;
        self.accounting.radiation_loss_j+=radiation_loss;
        self.accounting.damper_loss_j+=damper_loss;self.accounting.catch_loss_j+=catch_loss;
        self.accounting.max_balance_error_j=self.accounting.max_balance_error_j.max(balance.abs());
        Ok(())
    }

    /// Transactional audio sample, including jack timing and bending state.
    /// Prepared modal scratch is always restored from authoritative motion on
    /// each trial. No allocation is required to predict, accept or roll back.
    pub fn step(&mut self)->Result<f64,Error>{self.step_observed(None)}

    pub fn board_trace_len(&self)->usize {self.substeps*self.bank.board_count}

    /// Capture every substep's loaded-board end velocity, interleaved by mode.
    /// Only consume the trace on Ok: a refused step may overwrite this caller's
    /// scratch buffer, but restores all authoritative mechanical/material state.
    pub fn step_with_board_trace(&mut self,trace:&mut[f64])->Result<f64,Error>{
        if trace.len()!=self.board_trace_len(){return Err(Error::InvalidControl);}
        self.step_observed(Some(trace))
    }
    fn step_observed(&mut self,mut trace:Option<&mut[f64]>)->Result<f64,Error>{
        self.saved_q.copy_from_slice(&self.bank.q);self.saved_v.copy_from_slice(&self.bank.v);
        self.saved_hammers.copy_from_slice(&self.hammers);self.saved_contacts.clone_from_slice(&self.contacts);
        let saved=self.accounting;let modal=self.bank.last_modal_loss_j;
        if let Some(air)=&mut self.radiation {air.checkpoint();}
        let mut average=0.0;
        for substep in 0..self.substeps {
            if let Err(error)=self.mechanics_step(){
                self.bank.q.copy_from_slice(&self.saved_q);self.bank.v.copy_from_slice(&self.saved_v);
                self.hammers.copy_from_slice(&self.saved_hammers);self.contacts.clone_from_slice(&self.saved_contacts);
                if let Some(air)=&mut self.radiation {air.restore();}
                self.accounting=saved;self.bank.last_modal_loss_j=modal;return Err(error);
            }
            if let Some(buffer)=trace.as_deref_mut(){
                let r=self.bank.board_count;
                buffer[substep*r..(substep+1)*r].copy_from_slice(&self.bank.v[self.bank.modes.len()..]);
            }
            average+=self.bank.volume_velocity();
        }
        // Diagnostic box decimation; physical pressure consumes the full trace.
        Ok(average/self.substeps as f64)
    }
}

#[cfg(test)]
#[path = "hammer_footprint_engine_tests.rs"]
mod footprint_tests;

#[cfg(test)]
#[path = "radiation_engine_tests.rs"]
mod radiation_tests;

#[cfg(test)]
#[path = "damper_engine_tests.rs"]
mod damper_tests;

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
        assert!(a.accounting.felt_relaxation_loss_j>0.0);
        assert!(a.accounting.felt_loss_j>=a.accounting.felt_relaxation_loss_j);
        let balance=a.accounting.input_work_j-a.accounting.dissipated_j()-a.energy_j();
        assert!(balance.abs()<1e-7,"{balance:e}");
    }
    #[test]
    fn invalid_controls_and_failed_samples_preserve_state(){
        let mut p=instrument();assert!(p.note_on(69,f64::NAN).is_err());assert!(p.set_sustain(1.1).is_err());
        p.note_on(69,2.0).unwrap();
        for _ in 0..200 {p.step().unwrap();}
        let before=p.energy_j();let q=p.bank.q.clone();let memory=p.contacts[0].memory;
        p.damper_drag_ns_m=f64::NAN;assert!(p.step().is_err());assert_eq!(p.energy_j(),before);
        assert_eq!(p.bank.q,q);assert_eq!(p.contacts[0].memory,memory);
    }
    #[test]
    fn sostenuto_captures_only_keys_held_on_its_rising_edge(){
        let mut p=instrument();p.set_sostenuto(true);p.note_on(69,1.0).unwrap();
        assert!(!p.hammers[0].latched);p.set_sostenuto(false);p.set_sostenuto(true);
        assert!(p.hammers[0].latched);p.note_off(69).unwrap();assert!(p.hammers[0].latched);
        p.set_sostenuto(false);assert!(!p.hammers[0].latched);
    }
    #[test]
    fn airborne_felt_recovers_without_erasing_permanent_crush(){
        let mut p=instrument();p.note_on(69,3.0).unwrap();let mut peak:f64=0.0;
        for _ in 0..4800 {
            p.step().unwrap();peak=peak.max(p.creep[0].deformation(&p.contacts[0].memory));
        }
        assert!(peak>1e-7);
        assert!(p.creep[0].deformation(&p.contacts[0].memory)<peak*1e-4);
        assert!(p.contacts[0].state.eps_max>0.0);
    }
    #[test]
    fn empty_prony_keeps_the_rate_independent_image_available(){
        let scale=super::super::geometry::demonstration_scale().unwrap();
        let card=GeneralizedMaxwell::new(5e6,vec![]).unwrap();
        let mut p=Instrument::new_with_felt(vec![scale[48]],&super::super::board::demonstration(),
            48_000,4,12,true,felt::demonstration_law().unwrap(),&card).unwrap();
        p.note_on(69,2.0).unwrap();for _ in 0..1000 {p.step().unwrap();}
        assert!(p.accounting.felt_loss_j>0.0);assert_eq!(p.accounting.felt_relaxation_loss_j,0.0);
    }
    #[test]
    fn independent_course_materials_keep_contact_history_and_close_energy(){
        let scale=super::super::geometry::demonstration_scale().unwrap();
        let courses=vec![scale[27],scale[48]];
        let first=felt::demonstration_law().unwrap();let mut second=first.clone();second.f_ref*=2.0;
        let mut p=Instrument::new_with_course_felts(courses.clone(),&super::super::board::demonstration(),
            48_000,4,12,true,vec![(first,relaxation::demonstration_prony()),
                (second,relaxation::demonstration_prony())]).unwrap();
        assert_eq!(p.laws[1].f_ref,2.0*p.laws[0].f_ref);
        p.note_on(courses[0].midi,2.0).unwrap();p.note_on(courses[1].midi,2.0).unwrap();
        for _ in 0..2400 {p.step().unwrap();}
        for course in 0..2 {
            assert!(p.bank.contact_strings.iter().enumerate().any(|(i,&s)|
                p.bank.strings[s].course==course&&p.contacts[i].state.eps_max>0.0));
        }
        assert!((p.accounting.input_work_j-p.accounting.dissipated_j()-p.energy_j()).abs()<1e-7);
        assert!(Instrument::new_with_course_felts(courses,&super::super::board::demonstration(),
            48_000,4,12,true,vec![]).is_err());
    }
    #[test]
    fn jack_touch_drives_contact_disengages_and_preserves_energy_and_rollback() {
        let c=super::super::geometry::demonstration_scale().unwrap()[48];
        for (peak,duration) in [(70.0,0.007),(30.0,0.1)] {
            let mut p=Instrument::new_with_course_shanks(vec![c],&super::super::board::demonstration(),
                48_000,4,12,true,vec![(felt::demonstration_law().unwrap(),relaxation::demonstration_prony())],
                ShankGeometry::published()).unwrap();
            p.jack_on(69,peak,duration).unwrap();
            for _ in 0..120 {p.step().unwrap();}
            let motion=p.hammers[0].motion;let elapsed=p.hammers[0].jack.elapsed_s;let input=p.accounting.input_work_j;
            p.damper_drag_ns_m=f64::NAN;assert!(p.step().is_err());
            assert_eq!(p.hammers[0].motion,motion);assert_eq!(p.hammers[0].jack.elapsed_s,elapsed);
            assert_eq!(p.accounting.input_work_j,input);p.damper_drag_ns_m=0.4;
            for _ in 0..7200 {p.step().unwrap();}
            assert!(p.contacts.iter().any(|c|c.state.eps_max>0.0));
            assert_eq!(p.hammers[0].jack.peak_n,0.0);
            assert!(p.accounting.input_work_j>0.0);assert!(p.accounting.shank_loss_j>0.0);
            assert!((p.accounting.input_work_j-p.accounting.dissipated_j()-p.energy_j()).abs()<1e-7);
        }
        assert!(instrument().jack_on(69,70.0,0.007).is_err());
    }
    #[test]
    fn observed_step_preserves_dynamics_and_every_substep_velocity() {
        let mut a=instrument();let mut b=instrument();
        a.note_on(69,2.0).unwrap();b.note_on(69,2.0).unwrap();
        let r=a.bank.board_count;let mut trace=vec![0.0;a.board_trace_len()];
        for _ in 0..1000 {
            let volume=a.step_with_board_trace(&mut trace).unwrap();let mut average=0.0;
            for substep in 0..b.substeps {
                b.mechanics_step().unwrap();average+=b.bank.volume_velocity();
                assert_eq!(&trace[substep*r..(substep+1)*r],&b.bank.v[b.bank.modes.len()..]);
            }
            assert_eq!(volume,average/b.substeps as f64);
            assert_eq!(a.bank.q,b.bank.q);assert_eq!(a.bank.v,b.bank.v);
            assert_eq!(a.accounting.felt_loss_j,b.accounting.felt_loss_j);
        }
    }
    #[test]
    fn observed_step_rejects_bad_shape_and_keeps_mechanical_rollback() {
        let mut p=instrument();p.note_on(69,2.0).unwrap();
        let q=p.bank.q.clone();let v=p.bank.v.clone();let energy=p.energy_j();
        assert!(p.step_with_board_trace(&mut[]).is_err());
        assert_eq!(p.bank.q,q);assert_eq!(p.bank.v,v);assert_eq!(p.energy_j(),energy);
        let mut trace=vec![0.0;p.board_trace_len()];p.damper_drag_ns_m=f64::NAN;
        assert!(p.step_with_board_trace(&mut trace).is_err());
        assert_eq!(p.bank.q,q);assert_eq!(p.bank.v,v);assert_eq!(p.energy_j(),energy);
    }
}
