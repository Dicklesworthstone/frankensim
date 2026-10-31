//! Kirchhoff--Carrier extension with the piano's actual moving bridge.
//!
//! fs-nlmodal supplies the same fixed-interface stress channel used by wires.
//! For y=sum(phi*q_relative)+x/L*b, the extra chord slope contributes 2*b²/L²
//! to that channel. Pull forces back through q_relative=z-beta*b. Omitting
//! either term would make bridge work inconsistent with the string geometry.
//! The two transverse polarizations of ONE physical segment share the SUM of
//! their quadratic strains, hence one tension and one quartic storage term.
//! No new time propagator, retuned frequencies, or endpoint force lag.
use std::{collections::BTreeMap, io::Read};
use fs_nlmodal::{KcStringParams, kirchhoff_carrier_string};
use fs_math::det;
use super::{Bank, StringMode, StringPort};
use super::super::geometry::Course;

#[path = "longitudinal_string.rs"]
pub mod longitudinal;

pub const HEADER:&str="frankensim-piano-string-stretching-v1";
const MAX_BYTES:u64=64*1024;
/// Bounded nonlinear force iteration; this is not a real-time guarantee.
pub const MAX_TRIALS:usize=64;

#[derive(Clone,Copy,Debug)]
pub struct Material {
    /// Axial rigidity EA [N], NOT installed tension or bending rigidity EI.
    pub axial_rigidity_n:f64,
    /// Continuous-span moderate-slope bound. Never a displacement clamp.
    pub maximum_slope:f64,
}
/// Explicit retained longitudinal inertia and its material/domain budget.
#[derive(Clone,Copy,Debug)]
pub struct WaveMaterial {
    pub bandwidth_hz:f64,
    pub mode_cap:usize,
    pub damping_ratio:f64,
    pub maximum_strain:f64,
}
#[derive(Clone,Debug)]
pub struct Specification {
    pub courses:BTreeMap<u8,Option<Material>>,
    pub waves:BTreeMap<u8,WaveMaterial>,
    /// Axial participation row of each `waves` key, in the BARE board's mass-normalized coordinates
    /// (length = the board's mode count): the string termination's displacement along the string axis,
    /// positive when it extends the string, per unit modal coordinate. It is not in the file: the board
    /// supplies it (`BoardGeometry::bridge_attachments` with `MotionSurface::project_at`), and a `waves`
    /// key without its row refuses when the bank is configured. There is no default coupling.
    pub axial_rows:BTreeMap<u8,Vec<f64>>,
}
impl Specification {
    /// Attach the axial participation rows of the `waves` keys (see `axial_rows`). Every row must belong to
    /// a `waves` key, be finite and nonempty; a key may be given once.
    pub fn set_bridge_axial_rows(&mut self,rows:BTreeMap<u8,Vec<f64>>)->Result<(),String> {
        for (key,row) in &rows {
            if !self.waves.contains_key(key) {return Err(format!("key {key}: an axial participation row needs a waves row"));}
            if row.is_empty()||row.iter().any(|a|!a.is_finite()) {
                return Err(format!("key {key}: the axial participation row must be finite and nonempty"));
            }
        }
        self.axial_rows=rows;Ok(())
    }
    pub fn read(text:&str,courses:&[Course])->Result<Self,String> {
        if text.len() as u64>MAX_BYTES {return Err("string stretching file exceeds 64 KiB".into());}
        let mut rows=BTreeMap::new();let mut waves=BTreeMap::new();let mut header=false;
        for (line,raw) in text.lines().enumerate() {
            let row=raw.split('#').next().unwrap_or("").trim();if row.is_empty(){continue;}
            let bad=||format!("string stretching line {}: expected linear,key; stretch,key,EA_N,slope_bound; or waves,key,EA_N,slope_bound,band_hz,mode_cap,zeta,strain_bound",line+1);
            if !header {if row!=HEADER{return Err(bad());}header=true;continue;}
            let fields:Vec<_>=row.split(',').map(str::trim).collect();
            if !matches!(fields.len(),2|4|8){return Err(bad());}
            let key=fields[1].parse::<u8>().map_err(|_|bad())?;
            let axial = |token: &str| -> Result<f64, String> {
                let course = courses.iter().find(|course| course.midi == key)
                    .ok_or_else(|| format!("string stretching key {key} is absent from the scale"))?;
                match (token, course.bonded_section) {
                    ("section", Some(section)) => Ok(section.axial_rigidity_n()),
                    ("section", None) => Err(format!("key {key}: section EA requires a resolved bonded material section")),
                    (_, Some(_)) => Err(format!("key {key}: use section for EA; a separate number would duplicate the material input")),
                    (_, None) => token.parse().map_err(|_| bad()),
                }
            };
            let material=match (fields[0],fields.len()) {
                ("linear",2)=>None,
                ("stretch",4)=>Some(Material {axial_rigidity_n:axial(fields[2])?,
                    maximum_slope:fields[3].parse().map_err(|_|bad())?}),
                ("waves",8)=>{
                    waves.insert(key,WaveMaterial {
                        bandwidth_hz:fields[4].parse().map_err(|_|bad())?,
                        mode_cap:fields[5].parse().map_err(|_|bad())?,
                        damping_ratio:fields[6].parse().map_err(|_|bad())?,
                        maximum_strain:fields[7].parse().map_err(|_|bad())?,
                    });
                    Some(Material {axial_rigidity_n:axial(fields[2])?,
                        maximum_slope:fields[3].parse().map_err(|_|bad())?})
                },
                _=>return Err(bad()),
            };
            if rows.insert(key,material).is_some(){return Err(format!("duplicate string stretching key {key}"));}
        }
        if !header{return Err("missing string stretching header".into());}
        let s=Self{courses:rows,waves,axial_rows:BTreeMap::new()};s.validate(courses)?;Ok(s)
    }
    pub fn load(path:&str,courses:&[Course])->Result<Self,String> {
        let mut text=String::new();std::fs::File::open(path).map_err(|e|format!("{path}: {e}"))?
            .take(MAX_BYTES+1).read_to_string(&mut text).map_err(|e|format!("{path}: {e}"))?;
        Self::read(&text,courses)
    }
    fn validate(&self,courses:&[Course])->Result<(),String> {
        if courses.is_empty()||courses.len()>88||self.courses.len()!=courses.len()
            || courses.iter().enumerate().any(|(i,c)|courses[..i].iter().any(|p|p.midi==c.midi)) {
            return Err("string stretching must cover the complete unique scale".into());
        }
        for c in courses {
            c.validate()?;
            let p=self.courses.get(&c.midi).ok_or_else(||format!("missing string stretching key {}",c.midi))?;
            if let Some(p)=p {
                if c.bonded_section.is_some_and(|section| p.axial_rigidity_n != section.axial_rigidity_n()) {
                    return Err(format!("key {}: stretching EA is stale relative to the bound material section; resolve the specification again", c.midi));
                }
                if !p.axial_rigidity_n.is_finite()||p.axial_rigidity_n<=0.
                    || !p.maximum_slope.is_finite()||p.maximum_slope<=0.||p.maximum_slope>0.3 {
                    return Err("string stretching requires finite positive EA and slope bound in (0,0.3]".into());
                }
            }
        }
        for (&key,w) in &self.waves {
            if !self.courses.get(&key).is_some_and(Option::is_some)
                || !w.bandwidth_hz.is_finite() || w.bandwidth_hz<=0.
                || !(1..=32).contains(&w.mode_cap)
                || !w.damping_ratio.is_finite() || !(0.0..=1.0).contains(&w.damping_ratio)
                || !w.maximum_strain.is_finite() || w.maximum_strain<=0. || w.maximum_strain>0.1 {
                return Err("longitudinal waves require a stretching material, positive band, 1..32 mode cap, zeta in [0,1], and strain bound in (0,0.1]".into());
            }
        }
        Ok(())
    }
}

#[derive(Clone,Copy,Debug,Default)]
pub struct Observation {
    pub additional_strain:f64,
    pub tension_n:f64,
    pub slope_bound:f64,
    /// ONE physical segment's total, already included in bank energy. Querying
    /// either polarization returns this same quantity; do not sum it twice.
    pub stretching_energy_j:f64,
    /// Retained axial coordinates for this ONE physical segment.
    pub longitudinal_modes:usize,
    pub longitudinal_kinetic_energy_j:f64,
    /// Conservative continuous-span bounds, not collocation samples.
    pub minimum_tension_n:f64,
    pub maximum_tension_n:f64,
}
struct Channel {
    string:usize,
    second:Option<usize>,
    coefficient:f64,
    diagonal:Vec<f64>,
    chord:f64,
    length:f64,
    rest_tension:f64,
    material:Material,
    waves:Option<longitudinal::Prepared>,
}
impl Channel {
    fn plane_values(&self,si:usize,q:&[f64],strings:&[StringPort],modes:&[StringMode])->(f64,f64,f64) {
        let s=&strings[si];let b=s.bridge.iter().zip(&q[modes.len()..]).map(|(g,x)|g*x).sum::<f64>();
        let mut strain= self.chord*b*b;let mut slope=b.abs()/self.length;
        for (k,e) in s.modes.clone().zip(&self.diagonal) {
            let relative=q[k]-modes[k].beta*b;
            strain+=e*relative*relative;slope+=det::sqrt(*e)*relative.abs();
        }
        (strain,b,slope)
    }
    fn values(&self,q:&[f64],strings:&[StringPort],modes:&[StringMode])->(f64,f64) {
        let (mut strain,_,mut slope)=self.plane_values(self.string,q,strings,modes);
        if let Some(si)=self.second {
            let (extra,_,bound)=self.plane_values(si,q,strings,modes);
            strain+=extra;slope=slope.hypot(bound);
        }
        (strain,slope)
    }
    fn observe(&self,q:&[f64],strings:&[StringPort],modes:&[StringMode])->Observation {
        let (s,slope_bound)=self.values(q,strings,modes);
        let tension_n=self.rest_tension+0.25*self.material.axial_rigidity_n*s;
        let mut result=Observation {additional_strain:0.25*s,tension_n,
            slope_bound,stretching_energy_j:0.25*self.coefficient*s*s,
            minimum_tension_n:tension_n,maximum_tension_n:tension_n,..Observation::default()};
        if let Some(w)=&self.waves {w.observe(q,strings,modes,&mut result);}
        result
    }
}

pub(super) struct Prepared {
    channels:Vec<Channel>,
    pub force:Vec<f64>,
    required:Vec<f64>,
    residual:f64,
    relaxation:f64,
}
impl Prepared {
    fn new(bank:&Bank,courses:&[Course],spec:&Specification)->Result<Self,String> {
        spec.validate(courses)?;let mut channels=Vec::new();
        let mut wave_coordinates=0usize;let mut wave_terms=0usize;
        for (si,s) in bank.strings.iter().enumerate() {
            if s.polarization!=0 {continue;} // this segment's second direction shares its first channel
            let c=courses.get(s.course).ok_or("missing original string course")?;
            let Some(material)=spec.courses[&c.midi] else {continue;};
            let length=if s.duplex{c.duplex_length_m}else{c.length_m};
            let cents=(s.member as f64-0.5*(c.unison-1) as f64)*c.detune_cents;
            let tension=c.tension_at_cents(cents)?;
            let second=bank.strings.iter().position(|p|p.course==s.course && p.member==s.member
                && p.duplex==s.duplex && p.polarization==1);
            if let Some(j)=second {
                let p=&bank.strings[j];
                if p.modes.len()!=s.modes.len() || s.modes.clone().zip(p.modes.clone()).any(|(a,b)|
                    bank.modes[a].omega!=bank.modes[b].omega || bank.modes[a].beta!=bank.modes[b].beta) {
                    return Err("physical string polarizations have incompatible span/material coordinates".into());
                }
            }
            // A zero-retained-mode duplex still has chord strain and endpoint
            // mass. Obtain the same n-independent coefficient without adding
            // an oscillator or deleting the segment from the physical model.
            let count=s.modes.len().max(1);
            let owner=kirchhoff_carrier_string(&KcStringParams {length,tension,
                lin_density:c.linear_density_kg_m,ea:material.axial_rigidity_n},count)
                .map_err(|e|e.to_string())?;
            if owner.channels.len()!=1 || owner.channels[0].coupling.len()!=count*count {
                return Err("string owner did not return its single complete stress channel".into());
            }
            let ch=&owner.channels[0];
            let diagonal:Vec<_>=(0..s.modes.len()).map(|i|ch.coupling[i*count+i]).collect();
            let chord=2./(length*length);
            if !ch.coefficient.is_finite()||ch.coefficient<=0.||!chord.is_finite()
                || diagonal.iter().any(|x|!x.is_finite()||*x<=0.) {
                return Err("derived string extension channel is unrepresentable".into());
            }
            let mut channel=Channel {string:si,second,coefficient:ch.coefficient,diagonal,chord,length,
                rest_tension:tension,material,waves:None};
            if let Some(w)=spec.waves.get(&c.midi) {
                // The termination's axial motion per board mode is geometry (arm, string axis, the board's
                // rotations), not a material constant: refuse a key that has no row rather than assume one.
                let row=spec.axial_rows.get(&c.midi).ok_or_else(||format!(
                    "key {}: waves needs the board's axial participation row (the string termination's motion along its axis per board mode); none was supplied and there is no default coupling",c.midi))?;
                if row.len()!=bank.board_count {
                    return Err(format!("key {}: the axial participation row has {} entries for {} board modes",c.midi,row.len(),bank.board_count));
                }
                let axial=bank.project_board_shape(row)?;
                let prepared=longitudinal::Prepared::new(&channel,c.linear_density_kg_m,bank.rate,*w,axial)?;
                wave_coordinates+=prepared.mode_count();
                wave_terms+=prepared.coupling_terms();
                if wave_coordinates>8192 || wave_terms>4_000_000 {
                    return Err("longitudinal string bank exceeds 8192 axial coordinates or four million coupling terms".into());
                }
                // Zero modes is an explicit band-limited quasistatic channel.
                // Keep its strain/domain admission even on a short duplex.
                channel.waves=Some(prepared);
            }
            channels.push(channel);
        }
        Ok(Self {channels,force:vec![0.;bank.q.len()],required:vec![0.;bank.q.len()],
            residual:f64::INFINITY,relaxation:1.})
    }
    pub(super) fn energy(&self,q:&[f64],strings:&[StringPort],modes:&[StringMode])->f64 {
        self.channels.iter().map(|c|{
            let o=c.observe(q,strings,modes);o.stretching_energy_j+o.longitudinal_kinetic_energy_j
        }).sum()
    }
    pub(super) fn candidate_energy(&self,q:&[f64],strings:&[StringPort],modes:&[StringMode])->f64 {
        self.channels.iter().map(|c|match &c.waves {
            Some(w)=>w.energy(q,strings,modes,true),
            None=>c.observe(q,strings,modes).stretching_energy_j,
        }).sum()
    }
    fn begin(&mut self) {
        self.force.fill(0.);self.residual=f64::INFINITY;self.relaxation=1.;
        for c in &mut self.channels {if let Some(w)=&mut c.waves {w.begin();}}
    }
    pub(super) fn predict(&mut self) {
        for c in &mut self.channels {if let Some(w)=&mut c.waves {w.predict();}}
    }
    pub(super) fn modal_loss(&self)->f64 {
        self.channels.iter().filter_map(|c|c.waves.as_ref()).map(|w|w.loss_j).sum()
    }
    pub(super) fn commit(&mut self) {
        for c in &mut self.channels {if let Some(w)=&mut c.waves {w.commit();}}
    }
    /// Exact discrete gradient of c/4*(q^T E q)^2, pulled through BOTH moving-
    /// boundary directions. Its scalar tension uses the total geometric strain.
    fn evaluate(&mut self,a:&[f64],b:&[f64],strings:&[StringPort],modes:&[StringMode])->Result<(),&'static str> {
        if a.len()!=self.force.len()||b.len()!=a.len()||a.iter().chain(b).any(|v|!v.is_finite()) {
            return Err("string stretching requires complete finite mechanical coordinates");
        }
        self.required.fill(0.);
        for c in &mut self.channels {
            if let Some(w)=&mut c.waves {
                w.evaluate(a,b,strings,modes,&mut self.required)?;
                continue;
            }
            let (sa,_) =c.values(a,strings,modes);let (sb,_) =c.values(b,strings,modes);
            let scalar=c.coefficient*f64::midpoint(sa,sb);
            for si in std::iter::once(c.string).chain(c.second) {
                let s=&strings[si];
                let (_,ba,_)=c.plane_values(si,a,strings,modes);let (_,bb,_)=c.plane_values(si,b,strings,modes);
                let bridge=f64::midpoint(ba,bb);
                let mut reaction=scalar*c.chord*bridge;
                for (k,e) in s.modes.clone().zip(&c.diagonal) {
                    let relative=f64::midpoint(a[k]-modes[k].beta*ba,b[k]-modes[k].beta*bb);
                    let g=scalar*e*relative;self.required[k]-=g;reaction-=modes[k].beta*g;
                }
                for (g,out) in s.bridge.iter().zip(&mut self.required[modes.len()..]) {*out-=g*reaction;}
            }
        }
        if self.required.iter().any(|x|!x.is_finite()){return Err("string stretching force overflow");}
        Ok(())
    }
    fn validate_state(&self,q:&[f64],strings:&[StringPort],modes:&[StringMode])->Result<(),&'static str> {
        for c in &self.channels {
            let o=c.observe(q,strings,modes);
            if ![o.slope_bound,o.tension_n,o.stretching_energy_j].iter().all(|x|x.is_finite())
                || o.slope_bound>c.material.maximum_slope {
                return Err("piano string exceeds its finite moderate-slope domain");
            }
        }
        Ok(())
    }
    fn correct(&mut self,a:&[f64],b:&[f64],strings:&[StringPort],modes:&[StringMode])->Result<bool,&'static str> {
        self.validate_state(a,strings,modes)?;
        self.evaluate(a,b,strings,modes)?;
        let mut residual=0.0_f64;let mut scale=1.0_f64;
        for (&want,&got) in self.required.iter().zip(&self.force) {
            residual=residual.max((want-got).abs());scale=scale.max(want.abs());
        }
        for c in &self.channels {if let Some(w)=&c.waves {
            w.validate(a,strings,modes,false)?;
            w.residual(&mut residual,&mut scale);
        }}
        if residual<=1e-10*scale {
            self.validate_state(b,strings,modes)?;
            for c in &self.channels {if let Some(w)=&c.waves {w.validate(b,strings,modes,true)?;}}
            return Ok(true); // force must still match the candidate it produced
        }
        // Underrelax only NUMERICAL iterates, never the accepted force law.
        // A difficult contraction still refuses at MAX_TRIALS; no frozen-tension
        // fallback, force repair, physical tolerance change or mode deletion.
        if residual>self.residual {self.relaxation*=0.5;}
        for (got,&want) in self.force.iter_mut().zip(&self.required) {
            *got+=self.relaxation*(want-*got);
        }
        for c in &mut self.channels {if let Some(w)=&mut c.waves {w.correct(self.relaxation);}}
        self.residual=residual;Ok(false)
    }
}

impl Bank {
    /// Cold selection using this bank's original scale/order. Call before any
    /// excitation; all-linear selection leaves the old image bit for bit.
    pub(crate) fn configure_string_stretching(&mut self,courses:&[Course],spec:&Specification)->Result<(),String> {
        if self.stretching.is_some()||self.q.iter().chain(&self.v).any(|q|*q!=0.) {
            return Err("string stretching requires a fresh unconfigured bank".into());
        }
        let p=Prepared::new(self,courses,spec)?;
        if !p.channels.is_empty(){self.stretching=Some(p);}Ok(())
    }
    pub fn has_string_stretching(&self)->bool {self.stretching.is_some()}
    pub fn longitudinal_mode_count(&self)->usize {
        self.stretching.as_ref().map_or(0,|p|p.channels.iter()
            .filter_map(|c|c.waves.as_ref()).map(longitudinal::Prepared::mode_count).sum())
    }
    /// Complete strain potential plus axial kinetic energy of wave-selected
    /// segments. Included once in energy(); not a separately additive ledger.
    pub fn longitudinal_energy_j(&self)->f64 {
        self.stretching.as_ref().map_or(0.,|p|p.channels.iter().filter_map(|c|c.waves.as_ref())
            .map(|w|w.energy(&self.q,&self.strings,&self.modes,false)).sum())
    }
    pub fn string_longitudinal_state(&self,string:usize)->Option<(&[f64],&[f64])> {
        let c=self.stretching.as_ref()?.channels.iter()
            .find(|c|c.string==string||c.second==Some(string))?;
        let w=c.waves.as_ref()?;Some((&w.q,&w.v))
    }
    pub(crate) fn checkpoint_longitudinal(&mut self) {
        if let Some(p)=&mut self.stretching {
            for c in &mut p.channels {if let Some(w)=&mut c.waves {w.checkpoint();}}
        }
    }
    pub(crate) fn restore_longitudinal(&mut self) {
        if let Some(p)=&mut self.stretching {
            for c in &mut p.channels {if let Some(w)=&mut c.waves {w.restore();}}
        }
    }
    pub fn string_stretching_observation(&self,string:usize)->Option<Observation> {
        let p=self.stretching.as_ref()?;
        p.channels.iter().find(|c|c.string==string || c.second==Some(string))
            .map(|c|c.observe(&self.q,&self.strings,&self.modes))
    }
    /// Reset trial scratch once per mechanical tick, including after refusal.
    pub fn begin_string_stretching_step(&mut self) {
        if let Some(p)=&mut self.stretching {p.begin();}
    }
    /// After predict/contact/finish, check the same-tick nonlinear force equation.
    /// False requires another predict/contact/finish, NOT a mechanical commit.
    /// Only a true result followed by the original energy gate can be published.
    pub fn correct_string_stretching_step(&mut self)->Result<bool,&'static str> {
        match &mut self.stretching {
            Some(p)=>p.correct(&self.q,&self.next_q,&self.strings,&self.modes),None=>Ok(true),
        }
    }
}

/// TEST FIXTURE ONLY (never compiled into a render): the axial coupling of the longitudinal mechanics tests,
/// a fixed fraction of the vertical bridge row. Those tests exercise energy closure, rollback and band
/// admission of the retained axial modes, not the physics of the coupling, so any nonzero row serves; this
/// one is the value they were written against. A render takes its row from the board
/// (`BoardGeometry::bridge_attachments`).
#[cfg(test)]
pub(crate) const FIXTURE_AXIAL_RATIO:f64=0.06;
#[cfg(test)]
impl Specification {
    /// Give every `waves` key the fixture row `FIXTURE_AXIAL_RATIO` x the board's vertical bridge row.
    pub(crate) fn with_fixture_axial_rows(mut self,board:&[super::BoardMode])->Self {
        let rows=self.waves.keys().map(|&key|(key,board.iter()
            .map(|m|FIXTURE_AXIAL_RATIO*m.bridge[usize::from(key-21)]).collect())).collect();
        self.set_bridge_axial_rows(rows).unwrap();self
    }
}

#[cfg(test)]
#[path="string_stretching_tests.rs"]
mod tests;

#[cfg(test)]
#[path="longitudinal_bank_tests.rs"]
mod longitudinal_tests;
