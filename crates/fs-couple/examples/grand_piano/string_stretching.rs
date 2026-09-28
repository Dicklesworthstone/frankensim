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
#[derive(Clone,Debug)]
pub struct Specification {pub courses:BTreeMap<u8,Option<Material>>}
impl Specification {
    pub fn read(text:&str,courses:&[Course])->Result<Self,String> {
        if text.len() as u64>MAX_BYTES {return Err("string stretching file exceeds 64 KiB".into());}
        let mut rows=BTreeMap::new();let mut header=false;
        for (line,raw) in text.lines().enumerate() {
            let row=raw.split('#').next().unwrap_or("").trim();if row.is_empty(){continue;}
            let bad=||format!("string stretching line {}: expected linear,key or stretch,key,EA_N,slope_bound",line+1);
            if !header {if row!=HEADER{return Err(bad());}header=true;continue;}
            let fields:Vec<_>=row.split(',').map(str::trim).collect();
            if !matches!(fields.len(),2|4){return Err(bad());}
            let key=fields[1].parse::<u8>().map_err(|_|bad())?;
            let material=match (fields[0],fields.len()) {
                ("linear",2)=>None,
                ("stretch",4)=>Some(Material {axial_rigidity_n:fields[2].parse().map_err(|_|bad())?,
                    maximum_slope:fields[3].parse().map_err(|_|bad())?}),
                _=>return Err(bad()),
            };
            if rows.insert(key,material).is_some(){return Err(format!("duplicate string stretching key {key}"));}
        }
        if !header{return Err("missing string stretching header".into());}
        let s=Self{courses:rows};s.validate(courses)?;Ok(s)
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
                if !p.axial_rigidity_n.is_finite()||p.axial_rigidity_n<=0.
                    || !p.maximum_slope.is_finite()||p.maximum_slope<=0.||p.maximum_slope>0.3 {
                    return Err("string stretching requires finite positive EA and slope bound in (0,0.3]".into());
                }
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
        Observation {additional_strain:0.25*s,tension_n:self.rest_tension+0.25*self.material.axial_rigidity_n*s,
            slope_bound,stretching_energy_j:0.25*self.coefficient*s*s}
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
            channels.push(Channel {string:si,second,coefficient:ch.coefficient,diagonal,chord,length,
                rest_tension:tension,material});
        }
        Ok(Self {channels,force:vec![0.;bank.q.len()],required:vec![0.;bank.q.len()],
            residual:f64::INFINITY,relaxation:1.})
    }
    pub(super) fn energy(&self,q:&[f64],strings:&[StringPort],modes:&[StringMode])->f64 {
        self.channels.iter().map(|c|c.observe(q,strings,modes).stretching_energy_j).sum()
    }
    fn begin(&mut self) {self.force.fill(0.);self.residual=f64::INFINITY;self.relaxation=1.;}
    /// Exact discrete gradient of c/4*(q^T E q)^2, pulled through BOTH moving-
    /// boundary directions. Its scalar tension uses the total geometric strain.
    fn evaluate(&mut self,a:&[f64],b:&[f64],strings:&[StringPort],modes:&[StringMode])->Result<(),&'static str> {
        if a.len()!=self.force.len()||b.len()!=a.len()||a.iter().chain(b).any(|v|!v.is_finite()) {
            return Err("string stretching requires complete finite mechanical coordinates");
        }
        self.required.fill(0.);
        for c in &self.channels {
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
        if residual<=1e-10*scale {
            self.validate_state(b,strings,modes)?;
            return Ok(true); // force must still match the candidate it produced
        }
        // Underrelax only NUMERICAL iterates, never the accepted force law.
        // A difficult contraction still refuses at MAX_TRIALS; no frozen-tension
        // fallback, force repair, physical tolerance change or mode deletion.
        if residual>self.residual {self.relaxation*=0.5;}
        for (got,&want) in self.force.iter_mut().zip(&self.required) {
            *got+=self.relaxation*(want-*got);
        }
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

#[cfg(test)]
#[path="string_stretching_tests.rs"]
mod tests;
