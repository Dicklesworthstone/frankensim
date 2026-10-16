//! Exact sample-addressed physical edge forces. No frequency or state reset.
use super::Error;
use fs_couple::render::plate::impact::dilation::DilationConfig;

#[derive(Clone,Debug)]
pub struct Tensioner {
    pub config:DilationConfig,
    pub drag_n_s_m:f64,
    events:Vec<(u64,f64)>,
}
impl Tensioner {
    /// Held physical force before the named mechanical sample. No mutable
    /// cursor: a rejected solve reads the same assignment on exact retry.
    pub fn force_n(&self,sample:u64)->f64 {
        let end=self.events.partition_point(|e|e.0<=sample);
        if end==0 {0.0}else{self.events[end-1].1}
    }
    pub fn admit(&self,steps:u64)->Result<(),Error> {
        let c=self.config;
        if ![c.effective_mass_kg,c.hardware_stiffness_n_m,c.maximum_abs_strain,self.drag_n_s_m].iter().all(|v|v.is_finite())
            || c.effective_mass_kg<=0.0 || c.hardware_stiffness_n_m<0.0 || self.drag_n_s_m<0.0
            || !(0.0<c.maximum_abs_strain && c.maximum_abs_strain<0.05)
            || !(self.drag_n_s_m/c.effective_mass_kg).is_finite()
            || self.events.len()>65_536 || self.events.windows(2).any(|e|e[1].0<=e[0].0)
            || self.events.iter().any(|e|e.0>=steps || !e.1.is_finite()
                || !(e.1/c.effective_mass_kg.sqrt()).is_finite() || (e.1/c.effective_mass_kg.sqrt()).abs()>1e6) {
            return Err("edge actuation needs positive finite total mass, passive stiffness/drag, a small strain cap and ordered in-budget physical force assignments".into());
        }
        Ok(())
    }
}
/// Remove only these options after complete parse; malformed input is unchanged.
/// All times are exact MECHANICAL sample indices, including in WAV mode.
pub fn option(args:&mut Vec<String>,maximum_slope:f64)->Result<Option<Tensioner>,Error> {
    let mut config=None;let mut events=Vec::new();let mut rest=Vec::new();let mut i=0;
    while i<args.len() {
        match args[i].as_str() {
            "--edge-tensioner"=>{
                if config.is_some() || args.len()-i<5 {return Err("one --edge-tensioner MASS_KG EXTRA_STIFFNESS_N_M DRAG_NS_M MAX_ABS_STRAIN is required".into());}
                let mut v=[0.0;4];for (v,s) in v.iter_mut().zip(&args[i+1..i+5]) {*v=s.parse::<f64>()?;}
                config=Some((DilationConfig {effective_mass_kg:v[0],hardware_stiffness_n_m:v[1],
                    maximum_abs_strain:v[3],maximum_slope},v[2]));i+=5;
            },
            "--edge-force"=>{
                if args.len()-i<3 || events.len()==65_536 {return Err("--edge-force requires SAMPLE FORCE_N within 65536 assignments".into());}
                events.push((args[i+1].parse::<u64>()?,args[i+2].parse::<f64>()?));i+=3;
            },
            _=>{rest.push(args[i].clone());i+=1;},
        }
    }
    let Some((config,drag_n_s_m))=config else {
        if !events.is_empty() {return Err("--edge-force requires an explicit --edge-tensioner".into());}
        return Ok(None);
    };
    let program=Tensioner {config,drag_n_s_m,events};program.admit(u64::MAX)?;*args=rest;Ok(Some(program))
}
#[cfg(test)]
mod tests;
