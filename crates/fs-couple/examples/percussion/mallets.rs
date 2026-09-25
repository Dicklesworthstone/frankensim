//! Finite circular felt faces on the real head or curved-shell displacement.
//! The existing WoolFelt/Kelvin/impact owners supply all material history and
//! time evolution. No fitted strike pulse, parallel Hertz tip, or output gain.
use super::{Error, ModePair, Stroke, TensionedDisk};
use fs_couple::render::plate::impact::compliant::{CompliantJaw, PadSide, PadSite};
use fs_couple::render::plate::impact::compliant::striker::FeltStriker;
use fs_couple::render::plate::impact::felt::KelvinBranch;
use fs_couple::render::plate::impact::linear::wire::film_shapes;
use fs_material::fiber::WoolFelt;
use std::{collections::BTreeMap, io::Read};

#[path="mallet_shell.rs"]
mod shell;
#[path="mallet_shaft.rs"]
mod shaft;

const HEADER:&str="frankensim-felt-mallet-v1";
const MAX_BYTES:u64=65_536;
#[derive(Clone,Debug)]
pub struct Spec { radius_m:f64, jaw:CompliantJaw, head:Option<Head> }
/// V2 uses actual head mass (geometry), rotary inertia and the shaft's XY
/// azimuth. V1's effective mass is never reinterpreted as a physical head.
#[derive(Clone,Copy,Debug)]
struct Head { rotary_kg_m2:f64, azimuth_rad:f64 }
impl Spec {
    pub fn parse(text:&str)->Result<Self,Error> {
        if text.len() as u64>MAX_BYTES {return Err("mallet specification exceeds 64 KiB".into());}
        let mut rows=text.lines().enumerate().filter_map(|(i,l)| {
            let l=l.split('#').next().unwrap_or("").trim();(!l.is_empty()).then_some((i+1,l))
        });
        let physical_head=match rows.next().map(|(_,r)|r) {
            Some(HEADER)=>false,
            Some("frankensim-felt-mallet-v2")=>true,
            _=>return Err("expected frankensim-felt-mallet-v1 or -v2".into()),
        };
        let mut geometry=None;let mut law=None;let mut conditioning=None;let mut creep=Vec::new();let mut head=None;
        for (line,row) in rows {
            let f:Vec<_>=row.split(',').map(str::trim).collect();
            let bad=||format!("invalid felt-mallet record at line {line}");
            let value=|i:usize|->Result<f64,Error> {
                let v=f.get(i).ok_or_else(bad)?.parse::<f64>()?;
                if !v.is_finite() {return Err(bad().into());}Ok(v)
            };
            match (f[0],f.len()) {
                ("geometry",5) if geometry.is_none()=>geometry=Some([value(1)?,value(2)?,value(3)?,value(4)?]),
                ("attachment",3) if physical_head && head.is_none()=>head=Some(Head {
                    rotary_kg_m2:value(1)?,azimuth_rad:value(2)?,
                }),
                ("felt",7) if law.is_none()=>law=Some(WoolFelt::new(value(1)?,value(2)?,value(3)?,value(4)?,value(5)?,value(6)?)?),
                ("conditioning",2) if conditioning.is_none()=>conditioning=Some(value(1)?),
                ("creep",3) if creep.len()<4=>creep.push(KelvinBranch {stiffness_n_m:value(1)?,viscosity_n_s_m:value(2)?}),
                _=>return Err(bad().into()),
            }
        }
        let [mass_kg,radius_m,thickness_m,initial_gap_m]=geometry.ok_or("mallet geometry record is required")?;
        if radius_m<=0.0 || !(std::f64::consts::PI*radius_m*radius_m).is_finite()
            || std::f64::consts::PI*radius_m*radius_m/4.0==0.0 {
            return Err("mallet face radius/area must be positive and representable".into());
        }
        if physical_head && head.is_none() {return Err("v2 mallet requires explicit attachment inertia and shaft azimuth".into());}
        if head.is_some_and(|h| h.rotary_kg_m2<=0.0 || h.azimuth_rad.abs()>std::f64::consts::PI) {
            return Err("v2 head needs positive rotary inertia and shaft azimuth in [-pi,pi] radians".into());
        }
        let result=Self {radius_m,head,jaw:CompliantJaw {side:PadSide::Negative,mass_kg,drag_n_s_m:0.0,
            initial_gap_m,initial_velocity_m_s:0.0,thickness_m,
            law:law.ok_or("mallet felt record is required")?,
            prior_maximum_strain:conditioning.ok_or("mallet conditioning record is required")?,creep}};
        // Existing owner admits all material/history/mass/creep fields before FEM.
        result.from_rows(&vec![vec![1.0];4],0,2,0.0)?;
        Ok(result)
    }
    pub fn load(path:&str)->Result<Self,Error> {
        let mut text=String::new();std::fs::File::open(path)?.take(MAX_BYTES+1).read_to_string(&mut text)?;
        Self::parse(&text)
    }
    fn points(&self,center:[f64;2])->[[f64;2];4] {
        let a=self.radius_m/std::f64::consts::SQRT_2;
        [[center[0]+a,center[1]],[center[0],center[1]+a],
         [center[0]-a,center[1]],[center[0],center[1]-a]]
    }
    fn from_rows(&self,rows:&[Vec<f64>],coordinate:usize,total:usize,speed:f64)->Result<FeltStriker,Error> {
        if rows.len()!=4 || total>fs_couple::render::plate::impact::MAX_IMPACT_MODES || !speed.is_finite() || !(0.0..=20.0).contains(&speed) {
            return Err("felt mallet requires four sites and an admitted physical launch speed".into());
        }
        let area=std::f64::consts::PI*self.radius_m*self.radius_m/4.0;
        let mut sites=Vec::with_capacity(4);
        for row in rows {
            if row.is_empty() || row.len()+1>total || coordinate>=total
                || (coordinate>0 && coordinate<=row.len()) {
                return Err("mallet coordinate overlaps the batter basis or exceeds the layout".into());
            }
            let mut weights=vec![0.0;total];weights[1..1+row.len()].copy_from_slice(row);
            sites.push(PadSite {weights,area_m2:area});
        }
        let mut jaw=self.jaw.clone();jaw.initial_velocity_m_s=speed;
        Ok(FeltStriker::new(total,coordinate,&sites,&jaw)?)
    }
    pub fn compile(&self,film:&TensionedDisk,modes:&[ModePair],stroke:Stroke,
        coordinate:usize,total:usize)->Result<FeltStriker,Error> {
        let center=stroke.position_m.ok_or("felt mallets require an explicit strike XY position")?;
        if !center.iter().all(|v|v.is_finite()) {return Err("nonfinite mallet position".into());}
        // The full disk must fit, not just its quadrature nodes. TensionedDisk
        // has one convex polygonal boundary; test distance to every actual edge.
        let mut edges=BTreeMap::<(usize,usize),(usize,usize,usize)>::new();
        for tri in &film.mesh.tris {for j in 0..3 {
            let (a,b)=(tri[j],tri[(j+1)%3]);
            let entry=edges.entry((a.min(b),a.max(b))).or_insert((a,b,0));entry.2+=1;
        }}
        for &(a,b,count) in edges.values() {
            if count!=1 {continue;}
            let (a,b)=(film.mesh.nodes[a],film.mesh.nodes[b]);
            let (dx,dy)=(b.0-a.0,b.1-a.1);let length=dx.hypot(dy);
            // Use the known interior origin to orient even a reversed edge.
            let sign=(dx*(-a.1)-dy*(-a.0)).signum();
            let distance=sign*(dx*(center[1]-a.1)-dy*(center[0]-a.0))/length;
            if !distance.is_finite() || distance<=self.radius_m {
                return Err("entire mallet footprint must lie inside the actual head boundary; no rim contact or clipping".into());
            }
        }
        let rows=film_shapes(film,modes,&self.points(center))?;
        self.from_rows(&rows,coordinate,total,stroke.speed_m_s)
    }
}

#[derive(Default)]
pub struct Selection {pub first:Option<Spec>,pub second:Option<Spec>}
impl Selection {
    pub fn enabled(&self)->bool {self.first.is_some() || self.second.is_some()}
    pub fn load(paths:[Option<String>;2])->Result<Self,Error> {
        Ok(Self {first:paths[0].as_deref().map(Spec::load).transpose()?,
            second:paths[1].as_deref().map(Spec::load).transpose()?})
    }
    pub fn admit(&self,command:&str,first:Stroke,second:Option<Stroke>)->Result<(),Error> {
        if !self.enabled() {return Ok(());}
        if !matches!(command,"drum"|"drum-wav"|"drum-mic"|"drum-stretch"|"drum-stretch-wav"|"drum-stretch-mic"|
            "snare"|"snare-wav"|"snare-mic"|"snare-off"|"snare-off-wav"|"snare-off-mic"|
            "splash"|"splash-wav"|"splash-mic") {
            return Err("felt mallets require nonlinear-capable drum, snare or splash mechanics; modal-only and paired-hi-hat images are not selected here".into());
        }
        if self.first.is_some() && first.position_m.is_none() {
            return Err("--mallet-spec requires --strike-position-m X Y".into());
        }
        if self.second.is_some() && second.is_none() {
            return Err("--second-mallet-spec requires --second-stick-position-m X Y".into());
        }
        Ok(())
    }
}
pub fn options(args:&mut Vec<String>)->Result<[Option<String>;2],Error> {
    let mut paths=[None,None];let mut i=0;
    while i<args.len() {
        let index=match args[i].as_str() {"--mallet-spec"=>0,"--second-mallet-spec"=>1,_=>{i+=1;continue;}};
        if paths[index].is_some() {return Err("duplicate mallet specification".into());}
        let path=args.get(i+1).ok_or("mallet option needs a specification path")?;
        if path.starts_with("--") || path.is_empty() {return Err("mallet option needs a path, not another option".into());}
        paths[index]=Some(path.clone());args.drain(i..i+2);
    }
    Ok(paths)
}

#[cfg(test)]
#[path="mallets_tests.rs"]
mod tests;
