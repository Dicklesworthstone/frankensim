//! Complete bounded SI input for two separate shells and one axial mechanism.
use super::*;
use std::path::{Path,PathBuf};
use std::io::Read;
use std::fmt::Write as _;
const MAX_BYTES:usize=16384;

pub(super) struct Mount {
    pub radius:f64,pub area:f64,pub thickness:f64,pub precompression:f64,
    pub law:WoolFelt,pub prior:f64,pub k:f64,pub eta:f64,
}
pub(super) struct Source {pub mesh:bool,pub path:PathBuf}
pub(super) struct Spec {
    pub sources:[Source;2],pub separation:f64,pub carriage:[f64;3],pub damping:[f64;2],
    pub mounts:[Mount;2],pub contact:[f64;3],pub sites:Vec<(f64,f64,f64)>,
    pub strike:[f64;2],pub pedal:String,
}
impl Spec {
    pub fn load(path:&Path)->Result<Self,Error> {
        let mut s=String::new();std::fs::File::open(path)?.take((MAX_BYTES+1) as u64).read_to_string(&mut s)?;
        let mut out=Self::parse(&s)?;
        let base=path.parent().unwrap_or_else(||Path::new("."));
        for source in &mut out.sources {source.path=base.join(&source.path);}
        Ok(out)
    }
    pub fn shells(&self)->Result<[specimen::Specimen;2],Error> {
        let load=|s:&Source|if s.mesh {specimen::load_selection(None,Some(&s.path))?.ok_or_else(||"missing explicit hi-hat mesh".into())}else{specimen::Specimen::load(&s.path)};
        Ok([load(&self.sources[0])?,load(&self.sources[1])?])
    }
    pub fn parse(text:&str)->Result<Self,Error> {
        if text.len()>MAX_BYTES {return Err("hi-hat input exceeds 16 KiB".into());}
        let mut header=false;let mut sources:[Option<Source>;2]=[None,None];
        let mut mounts:[Option<Mount>;2]=[None,None];
        let (mut gap,mut carriage,mut damping,mut law,mut strike)=(None,None,None,None,None);
        let mut sites=Vec::new();let mut pedal=String::new();
        for (line,raw) in text.lines().enumerate() {
            let row=raw.split('#').next().unwrap_or("").trim();if row.is_empty(){continue;}
            let bad=||->Error{format!("hi-hat line {}: malformed, missing or duplicate record",line+1).into()};
            if !header {if row!="frankensim-hihat-v1"{return Err(bad());}header=true;continue;}
            let f:Vec<_>=row.split(',').map(str::trim).collect();
            let number=|i:usize|->Result<f64,Error>{
                let x:f64=f[i].parse().map_err(|_|bad())?;if !x.is_finite(){return Err(bad());}Ok(x)
            };
            match f[0] {
                "shell" if f.len()==4=>{
                    let i=side(f[1])?;if sources[i].is_some()||f[3].is_empty(){return Err(bad());}
                    let mesh=match f[2]{"profile"=>false,"mesh"=>true,_=>return Err(bad())};
                    sources[i]=Some(Source{mesh,path:f[3].into()});
                }
                "separation_m" if f.len()==2 && gap.is_none()=>gap=Some(number(1)?),
                "carriage" if f.len()==4 && carriage.is_none()=>carriage=Some([number(1)?,number(2)?,number(3)?]),
                "damping_ratio" if f.len()==3 && damping.is_none()=>damping=Some([number(1)?,number(2)?]),
                "contact" if f.len()==4 && law.is_none()=>law=Some([number(1)?,number(2)?,number(3)?]),
                "strike_m" if f.len()==3 && strike.is_none()=>strike=Some([number(1)?,number(2)?]),
                "site" if f.len()==4 && sites.len()<64=>sites.push((number(1)?,number(2)?,number(3)?)),
                "pedal" if f.len()==3=>{writeln!(pedal,"{},{}",number(1)?,number(2)?)?;}
                "mount" if f.len()==15=>{
                    let i=side(f[1])?;if mounts[i].is_some(){return Err(bad());}
                    mounts[i]=Some(Mount{radius:number(2)?,area:number(3)?,thickness:number(4)?,
                        precompression:number(5)?,law:WoolFelt::new(number(6)?,number(7)?,number(8)?,number(9)?,number(10)?,number(11)?)?,
                        prior:number(12)?,k:number(13)?,eta:number(14)?});
                }
                _=>return Err(bad()),
            }
        }
        let [u,l]=sources;let [um,lm]=mounts;
        let out=Self{sources:[u.ok_or("missing upper shell")?,l.ok_or("missing lower shell")?],
            mounts:[um.ok_or("missing upper mount")?,lm.ok_or("missing lower mount")?],
            separation:gap.ok_or("missing separation")?,carriage:carriage.ok_or("missing carriage")?,
            damping:damping.ok_or("missing damping")?,contact:law.ok_or("missing contact law")?,
            strike:strike.ok_or("missing default strike")?,sites,pedal};
        out.validate()?;Ok(out)
    }
    pub fn validate(&self)->Result<(),Error> {
        if !self.separation.is_finite()||self.separation<=0.||self.sites.is_empty()||self.sites.len()>64
            ||self.carriage.iter().any(|x|!x.is_finite())||self.carriage[0]<=0.||self.carriage[1]<=0.||self.carriage[2]<0.
            ||!(self.carriage[1]/self.carriage[0]).is_finite()||!(self.carriage[2]/self.carriage[0]).is_finite()
            ||self.damping.iter().any(|x|!x.is_finite()||*x<0.)||self.strike.iter().any(|x|!x.is_finite())
            ||self.sites.iter().any(|&(x,y,w)|!x.is_finite()||!y.is_finite()||!w.is_finite()||w<=0.)
            ||(self.sites.iter().map(|s|s.2).sum::<f64>()-1.).abs()>1e-12 {
            return Err("hi-hat requires positive separation/mass/return stiffness and positive unit-sum contact weights".into());
        }
        for (i,a) in self.sites.iter().enumerate(){if self.sites[..i].iter().any(|b|a.0==b.0&&a.1==b.1){return Err("duplicate hi-hat collision station".into());}}
        // The contact owner, not this adapter, admits alpha/chi and K.
        Obstacle::new(vec![1.],1,1,vec![0.],vec![1.],self.contact[0],self.contact[1],"supplied paired contact".into())?
            .with_internal_loss(self.contact[2])?;
        for m in &self.mounts {
            if [m.radius,m.area,m.thickness,m.k,m.eta].iter().any(|v|!v.is_finite()||*v<=0.)
                ||!m.precompression.is_finite()||m.precompression<0.||m.precompression/m.thickness>=m.law.eps_densify
                ||!m.prior.is_finite()||m.prior<0.||m.prior>m.law.eps_densify||!(m.k/m.eta).is_finite() {
                return Err("invalid hi-hat washer geometry, prior strain or Kelvin branch".into());
            }
        }
        mechanics::drive::Program::parse(&self.pedal)?;Ok(())
    }
}
fn side(s:&str)->Result<usize,Error>{match s{"upper"=>Ok(0),"lower"=>Ok(1),_=>Err("expected upper or lower".into())}}
