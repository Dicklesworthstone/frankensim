//! Authored moving compliant mutes. Geometry and contact force act on the actual
//! resonator; neither force knots nor jaw travel prescribe its microphone sound.
use super::{Error, Mechanics, muffling::{self, Muffler, Surface}, mechanics::drive};
use fs_couple::render::plate::impact::compliant::{CompliantJaw, JawPort, MovingPads, PadSide, PadSite};
use fs_couple::render::plate::impact::felt::KelvinBranch;
use fs_material::fiber::WoolFelt;
use std::io::{Read, Write};

const MAX_BYTES: u64 = 65_536;
/// Supplied bound on the whole jaw's XY projection, including its felt.
/// Sites remain contact quadrature; this disk is not an averaged modal port.
#[derive(Debug, Clone, Copy)]
struct Envelope { center_m: [f64;2], radius_m: f64 }
#[derive(Debug)]
pub struct Spec {
    pub surface: Surface,
    sites: Vec<([f64;2], f64)>,
    envelope: Option<Envelope>,
    jaws: Vec<CompliantJaw>,
    programs: Vec<drive::Program>,
}
fn side(text: &str) -> Result<PadSide, Error> {
    // All example surface rows point DOWN, including both drumhead bases.
    match text { "above"=>Ok(PadSide::Negative), "below"=>Ok(PadSide::Positive),
        _=>Err("compliant mute jaw side must be above or below".into()) }
}
fn side_index(s: PadSide) -> usize { usize::from(s == PadSide::Positive) }
fn side_name(s: PadSide) -> &'static str { if s == PadSide::Positive {"below"}else{"above"} }
impl Spec {
    pub fn parse(text: &str) -> Result<Self, Error> {
        if text.len() as u64 > MAX_BYTES { return Err("compliant mute exceeds 64 KiB".into()); }
        let mut rows=text.lines().enumerate().filter_map(|(i,l)| {
            let l=l.split('#').next().unwrap_or("").trim(); (!l.is_empty()).then_some((i+1,l))
        });
        if rows.next().map(|(_,s)|s)!=Some("frankensim-compliant-mute-v1") {
            return Err("compliant mute needs frankensim-compliant-mute-v1 header".into());
        }
        let mut surface=None;let mut sites=Vec::new();let mut envelope=None;
        let mut jaws:Vec<CompliantJaw>=Vec::new();
        let mut force=[String::new(),String::new()];let mut creep=[Vec::new(),Vec::new()];
        for (line,text) in rows {
            let f:Vec<_>=text.split(',').map(str::trim).collect();
            let bad=||format!("invalid compliant mute record at line {line}");
            let number=|i:usize| -> Result<f64,Error> {
                let value=f.get(i).ok_or_else(bad)?.parse::<f64>()?;
                if !value.is_finite() {return Err(bad().into());} Ok(value)
            };
            match (f[0],f.len()) {
                ("surface",2) if surface.is_none()=>surface=Some(match f[1] {
                    "shell"=>Surface::Shell,"batter"=>Surface::Batter,"resonant"=>Surface::Resonant,
                    _=>return Err(bad().into()),
                }),
                ("site",4) if sites.len()<4=> {
                    let area=number(3)?;if area<=0.0 {return Err(bad().into());}
                    sites.push(([number(1)?,number(2)?],area));
                }
                ("envelope",4) if envelope.is_none()=> {
                    envelope=Some(Envelope {center_m:[number(1)?,number(2)?],radius_m:number(3)?});
                }
                ("jaw",13) if jaws.len()<2=> {
                    let s=side(f[1])?;
                    if jaws.iter().any(|j|j.side==s) {return Err("duplicate compliant mute jaw".into());}
                    jaws.push(CompliantJaw {side:s,mass_kg:number(2)?,drag_n_s_m:number(3)?,
                        initial_gap_m:number(4)?,initial_velocity_m_s:0.0,thickness_m:number(5)?,
                        law:WoolFelt::new(number(6)?,number(7)?,number(8)?,number(9)?,number(10)?,number(11)?)?,
                        prior_maximum_strain:number(12)?,creep:Vec::new()});
                }
                ("creep",4)=> {
                    let i=side_index(side(f[1])?);
                    if creep[i].len()==4 {return Err("compliant mute exceeds four Kelvin elements per jaw".into());}
                    creep[i].push(KelvinBranch {stiffness_n_m:number(2)?,viscosity_n_s_m:number(3)?});
                }
                ("force",4)=> {
                    let i=side_index(side(f[1])?);
                    // Reuse the existing time/force parser, interval integration,
                    // force limits and accepted-tick cursor. No second scheduler.
                    writeln!(&mut force[i],"{},{}",f[2],f[3])?;
                }
                _=>return Err(bad().into()),
            }
        }
        let surface=surface.ok_or("compliant mute needs a surface")?;
        let mut programs=Vec::new();
        for j in &mut jaws {
            let i=side_index(j.side);j.creep=std::mem::take(&mut creep[i]);
            programs.push(drive::Program::parse(&std::mem::take(&mut force[i]))?);
        }
        if force.iter().any(|s|!s.is_empty()) || creep.iter().any(|v|!v.is_empty()) {
            return Err("compliant mute force or creep references an absent jaw".into());
        }
        // Validate all physical cards before any instrument eigensolve. This
        // one-coordinate row checks cards only; real geometry is lowered later.
        let checks:Vec<_>=sites.iter().map(|(_,area)|PadSite {weights:vec![1.0],area_m2:*area}).collect();
        MovingPads::new(1,&checks,&jaws)?;
        let result=Self {surface,sites,envelope,jaws,programs};
        result.admit_envelope()?;
        result.admit_command(if surface==Surface::Shell {"splash"}else{"drum"})?;
        Ok(result)
    }
    pub fn load(path: &str) -> Result<Self, Error> {
        let mut text=String::new();std::fs::File::open(path)?.take(MAX_BYTES+1).read_to_string(&mut text)?;
        Self::parse(&text)
    }
    pub fn jaw_count(&self) -> usize {self.jaws.len()}
    pub fn admit_command(&self, command: &str) -> Result<(),Error> {
        let shell=matches!(command,"splash"|"splash-wav"|"splash-mic");
        let snare=matches!(command,"snare"|"snare-wav"|"snare-mic"|"snare-off"|"snare-off-wav"|"snare-off-mic");
        let drum=snare || matches!(command,"drum"|"drum-wav"|"drum-mic"|"drum-stretch"|"drum-stretch-wav"|"drum-stretch-mic");
        if !(shell || drum) || (self.surface==Surface::Shell)!=shell {
            return Err("compliant mute requires the matching nonlinear-capable splash, drum or snare image".into());
        }
        if snare && self.surface==Surface::Resonant && self.envelope.is_none() {
            return Err("a resonant-head snare mute needs a whole-jaw envelope separated from the wire bank".into());
        }
        if !shell && (self.jaws.len()!=1 || self.jaws[0].side!=
            if self.surface==Surface::Batter {PadSide::Negative}else{PadSide::Positive}) {
            return Err("drum mutes must approach the exterior: batter above or resonant below; no interior pad displacing unmodelled cavity air".into());
        }
        Ok(())
    }
    fn admit_envelope(&self) -> Result<(),Error> {
        let Some(e)=self.envelope else {return Ok(());};
        let area=std::f64::consts::PI*e.radius_m*e.radius_m;
        if self.surface==Surface::Shell || !e.radius_m.is_finite() || e.radius_m<=0.0
            || e.center_m.iter().any(|x|!x.is_finite()) || !area.is_finite() || area<=0.0
            || self.sites.iter().map(|s|s.1).sum::<f64>()>area
            || self.sites.iter().any(|(p,area)|
                (p[0]-e.center_m[0]).hypot(p[1]-e.center_m[1])
                    +(area/std::f64::consts::PI).sqrt()>e.radius_m) {
            return Err("head-mute envelope needs a finite positive radius containing every circular contact cell and its total pad area".into());
        }
        Ok(())
    }
    /// Cold geometric exclusion before either head is assembled. All admitted
    /// wire, carrier and jaw motion is axial, so their disjoint XY projections
    /// stay disjoint throughout the performance; no jaw/wire collision is hidden.
    pub fn admit_drum_geometry(&self, radius_m:f64, wires:Option<super::snare::SnareSet>)
        -> Result<(),Error>
    {
        self.admit_command(if wires.is_some() {"snare"}else{"drum"})?;
        self.admit_envelope()?;
        let Some(e)=self.envelope else {return Ok(());};
        let extent=e.center_m[0].hypot(e.center_m[1])+e.radius_m;
        let margin=64.0*f64::EPSILON*(radius_m+extent);
        if !radius_m.is_finite() || radius_m<=0.0 || !extent.is_finite() || !margin.is_finite()
            || radius_m-extent<=margin {
            return Err("the entire head-mute envelope must lie strictly inside the clear head radius; rim contact is not represented".into());
        }
        if self.surface==Surface::Resonant {
            if let Some(wires)=wires {
                wires.mode_count()?;
                wires.coil.linear_density_kg_m()?;
                // Enclose every full strand, including the coil/wire radius.
                // The bank is a line when there is only one central strand.
                let collar=wires.coil.coil_radius_m+wires.coil.wire_radius_m;
                let half=[0.5*wires.length_m+collar,
                    if wires.strands==1 {collar}else{0.5*wires.width_m+collar}];
                let dx=(e.center_m[0].abs()-half[0]).max(0.0);
                let dy=(e.center_m[1].abs()-half[1]).max(0.0);
                let separation=dx.hypot(dy)-e.radius_m;
                let margin=64.0*f64::EPSILON*(extent+half[0]+half[1]);
                if half.iter().any(|x|!x.is_finite()) || !separation.is_finite()
                    || !margin.is_finite() || separation<=margin {
                    return Err("resonant-head mute envelope intersects or touches the complete wire bank; jaw/wire collision is not represented".into());
                }
            }
        }
        Ok(())
    }
    fn locations(&self) -> Vec<Muffler> {
        self.sites.iter().map(|(position_m,_)|Muffler {surface:self.surface,
            position_m:*position_m,resistance_n_s_m:0.0}).collect()
    }
    fn compile(&self, rows: Vec<fs_couple::render::plate::impact::damping::ViscousDamper>,
        structural: usize) -> Result<MovingPads, Error> {
        if rows.len()!=self.sites.len() || rows.iter().any(|r|r.weights.len()>structural) {
            return Err("compliant mute cannot truncate its geometric surface rows".into());
        }
        let sites:Vec<_>=rows.into_iter().zip(&self.sites).map(|(mut row,(_,area))| {
            row.weights.resize(structural,0.0);PadSite {weights:row.weights,area_m2:*area}
        }).collect();
        Ok(MovingPads::new(structural,&sites,&self.jaws)?)
    }
    pub fn shell(&self, reduction:&super::ShellReduction,nodes:&[[f64;3]],triangles:&[[usize;3]],
        structural:usize) -> Result<MovingPads,Error> {
        self.admit_command("splash")?;
        self.compile(muffling::shell_ports(&self.locations(),reduction,nodes,triangles)?,structural)
    }
    pub fn head(&self, films:&[super::TensionedDisk], modes:&[Vec<super::ModePair>],
        structural:usize) -> Result<MovingPads,Error> {
        self.admit_command("drum")?;
        self.compile(muffling::head_ports(&self.locations(),films,modes,structural)?,structural)
    }
    pub fn into_inputs(self, attachment:&Attachment) -> Result<Vec<drive::Input>,Error> {
        if self.programs.len()!=attachment.ports.len() {return Err("missing physical mute jaw ports".into());}
        Ok(self.programs.into_iter().zip(&attachment.ports).map(|(program,p)|drive::Input {
            program,coordinate:p.coordinate,tip_weight:p.inverse_sqrt_mass,
        }).collect())
    }
    pub fn observation(&self, ports:Vec<JawPort>, first_pad:usize) -> Attachment {
        Attachment {ports,first_pad,sites_per_jaw:self.sites.len()}
    }
}
// std::fmt::Write supplies String's write_fmt; std::io::Write supplies CSV output.
use std::fmt::Write as _;

pub fn option(args:&mut Vec<String>) -> Result<Option<Spec>,Error> {
    let mut found=None;let mut i=0;
    while i<args.len() {
        if args[i]!="--compliant-mute" {i+=1;continue;}
        if found.is_some() {return Err("--compliant-mute may be supplied only once".into());}
        let path=args.get(i+1).ok_or("--compliant-mute needs a specification path")?;
        if path.starts_with("--") {return Err("--compliant-mute needs a path, not an option".into());}
        found=Some(Spec::load(path)?);args.drain(i..i+2);
    }
    Ok(found)
}

pub struct Attachment {pub ports:Vec<JawPort>,first_pad:usize,sites_per_jaw:usize}
fn felt(system:&Mechanics, index:usize) -> Option<(f64,f64)> {
    match system {
        Mechanics::Reference(s)=>s.felt_observation(index),
        Mechanics::Nonlinear(s)=>s.felt_observation(index),
        Mechanics::Substepped(s)=>s.felt_observation(index),
        Mechanics::Driven {inner,..}=>felt(inner,index),Mechanics::Prepared(_)=>None,
    }
}
impl Attachment {
    pub fn header(&self,out:&mut impl Write) -> Result<(),Error> {
        for p in &self.ports {
            let side=side_name(p.side);
            write!(out,",mute_{side}_inward_m,mute_{side}_velocity_m_s,mute_{side}_contact_n")?;
        } Ok(())
    }
    pub fn row(&self,system:&Mechanics,out:&mut impl Write) -> Result<(),Error> {
        for (jaw,p) in self.ports.iter().enumerate() {
            let force=(0..self.sites_per_jaw).try_fold(0.0,|sum,i| -> Result<f64,Error> {
                Ok(sum+felt(system,self.first_pad+jaw*self.sites_per_jaw+i)
                    .ok_or("missing accepted mute contact state")?.1)
            })?;
            let x=system.state();write!(out,",{:.17e},{:.17e},{force:.17e}",
                x[2*p.coordinate]*p.inverse_sqrt_mass,x[2*p.coordinate+1]*p.inverse_sqrt_mass)?;
        } Ok(())
    }
}

#[cfg(test)]
mod tests;
