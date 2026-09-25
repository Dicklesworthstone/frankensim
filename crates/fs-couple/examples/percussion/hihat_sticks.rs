//! Explicit shaft geometry/material input for the existing flexible-striker owner.
//! Tip contact is still the separate declared percussion Hertz model.
use super::Error;
use fs_couple::render::plate::impact::striker::{RadiusStation, flexible::FlexibleStriker};
use fs_plate::shell::stiffened::beam::RoundBeamSpec;
use std::io::Read;

const MAX_BYTES: usize = 8192;
struct Spec {
    profile: Vec<RadiusStation>, beam: RoundBeamSpec, damping: f64,
}
impl Spec {
    fn parse(text: &str) -> Result<Self, Error> {
        if text.len()>MAX_BYTES {return Err("flexible stick exceeds 8 KiB".into());}
        let (mut header,mut material,mut support,mut basis)=(false,None,None,None);
        let mut profile=Vec::new();
        for (line,raw) in text.lines().enumerate() {
            let row=raw.split('#').next().unwrap_or("").trim();if row.is_empty(){continue;}
            let bad=||->Error{format!("flexible stick line {}: invalid or duplicate record",line+1).into()};
            if !header {if row!="frankensim-flexible-stick-v1"{return Err(bad());}header=true;continue;}
            let f:Vec<_>=row.split(',').map(str::trim).collect();
            let number=|i:usize|->Result<f64,Error>{let n:f64=f[i].parse().map_err(|_|bad())?;
                if !n.is_finite(){return Err(bad());}Ok(n)};
            match f[0] {
                "material" if f.len()==4 && material.is_none()=>material=Some([number(1)?,number(2)?,number(3)?]),
                "support" if f.len()==4 && support.is_none()=>support=Some([number(1)?,number(2)?,number(3)?]),
                "basis" if f.len()==4 && basis.is_none()=>{
                    let subdivisions:usize=f[1].parse().map_err(|_|bad())?;
                    let maximum_modes:usize=f[3].parse().map_err(|_|bad())?;
                    basis=Some((subdivisions,number(2)?,maximum_modes));
                }
                "station" if f.len()==3 && profile.len()<33=>profile.push(RadiusStation{position_m:number(1)?,radius_m:number(2)?}),
                _=>return Err(bad()),
            }
        }
        let [young_pa,density_kg_m3,damping]=material.ok_or("missing shaft material")?;
        let [pivot_m,contact_m,hand_m]=support.ok_or("missing shaft support/tip/hand stations")?;
        let (subdivisions,maximum_hz,maximum_modes)=basis.ok_or("missing shaft basis limits")?;
        if profile.len()<2 || young_pa<=0.0 || density_kg_m3<=0.0 || !(0.0..1.0).contains(&damping)
            || subdivisions==0 || subdivisions>16 || !(2..=17).contains(&maximum_modes) || maximum_hz<=0.0
            || (profile.len()-1)*subdivisions+2>66
            || profile.iter().any(|s|s.radius_m<0.0)
            || profile.windows(2).any(|w|w[0].position_m>=w[1].position_m || (w[0].radius_m==0.0 && w[1].radius_m==0.0)) {
            return Err("invalid bounded shaft geometry, material or complete frequency slice".into());
        }
        let start=profile[0].position_m;let end=profile[profile.len()-1].position_m;
        if !(start..end).contains(&pivot_m) || contact_m<=pivot_m || contact_m>end || hand_m<=pivot_m || hand_m>end {
            return Err("pivot, tip and hand must be on the supplied shaft with positive force levers".into());
        }
        Ok(Self{profile,beam:RoundBeamSpec{young_pa,density_kg_m3,pivot_m,contact_m,hand_m,
            subdivisions,maximum_hz,maximum_modes},damping})
    }
    fn prepare(self)->Result<FlexibleStriker,Error>{Ok(FlexibleStriker::new(&self.profile,self.beam,self.damping)?)}
}

pub(super) fn option(args:&mut Vec<String>,flag:&str)->Result<Option<FlexibleStriker>,Error> {
    let positions:Vec<_>=args.iter().enumerate().filter_map(|(i,s)|(s==flag).then_some(i)).collect();
    if positions.len()>1{return Err(format!("{flag} may be supplied only once").into());}
    let Some(&i)=positions.first() else{return Ok(None);};
    let path=args.get(i+1).filter(|s|!s.is_empty()&&!s.starts_with("--"))
        .ok_or_else(||format!("{flag} requires a shaft file; see FLEXIBLE_STICKS.md"))?;
    let mut text=String::new();std::fs::File::open(path)?.take((MAX_BYTES+1) as u64).read_to_string(&mut text)?;
    let shaft=Spec::parse(&text)?.prepare()?;args.drain(i..i+2);Ok(Some(shaft))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn complete_physical_input_has_no_inferred_stiffness_or_force_station() {
        let text=include_str!("estimated-flexible-stick.fst");assert!(Spec::parse(text).is_ok());
        for bad in [text.replace("material,12000000000,800,0.001","material,0,800,0.001"),
            text.replace("basis,8,3000,17","basis,100,3000,17"),
            text.replace("basis,8,3000,17","basis,8,3000,1"),
            text.replace("support,0.1,0.39,0.16","support,0.1,0.39,0.09"),
            text.replace("station,0.4,0.005","station,0,0.005"),
            text.replace("station,0.4,0.005","station,0.4,NaN"),
            text.replace("support,0.1,0.39,0.16", ""),format!("{text}\nunknown,0")]
        {assert!(Spec::parse(&bad).is_err(),"{bad}");}
        let mut args=vec!["hihat".into()];let before=args.clone();
        assert!(option(&mut args,"--flexible-stick").unwrap().is_none());assert_eq!(args,before);
        for mut args in [vec!["--flexible-stick".into()],
            vec!["--flexible-stick".into(),"--analytic-newton".into()],
            vec!["--flexible-stick".into(),"x".into(),"--flexible-stick".into(),"y".into()]] {
            let before=args.clone();assert!(option(&mut args,"--flexible-stick").is_err());assert_eq!(args,before);
        }
    }
}

pub(super) struct Launch {
    pub body: super::ImpactBody,
    pub weight: f64,
    pub elastic: Option<super::ImpactBody>,
    pub ports: Option<fs_couple::render::plate::impact::striker::flexible::StrikerPorts>,
}
impl Launch {
    pub(super) fn new(shaft:Option<&FlexibleStriker>,rigid:usize,elastic_start:usize,
        total:usize,dt:f64,speed:f64)->Result<Self,Error> {
        if !speed.is_finite() || !(0.0..=20.0).contains(&speed) {
            return Err("strike speed must be finite in 0..=20 m/s".into());
        }
        match shaft {
            Some(shaft)=>{
                let ports=shaft.ports(rigid,elastic_start,total,dt)?;
                let (body,elastic)=shaft.split_bodies(-0.0002,speed)?;
                Ok(Self{body,weight:shaft.tip_weights()[0],elastic:Some(elastic),ports:Some(ports)})
            }
            None=>{let (body,weight)=super::stick_with_speed(speed)?;
                Ok(Self{body,weight,elastic:None,ports:None})}
        }
    }
    pub(super) fn tip_row(&self,rigid:usize,total:usize)->Result<Vec<f64>,Error> {
        match &self.ports {
            Some(p)=>Ok(p.tip_row(total)?),
            None=>{if rigid>=total{return Err("rigid stick lies outside the force basis".into());}
                let mut row=vec![0.0;total];row[rigid]=self.weight;Ok(row)}
        }
    }
}
