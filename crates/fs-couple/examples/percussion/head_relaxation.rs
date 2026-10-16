//! Hereditary bending from supplied material moduli and the actual head FEM.
//! Installed tension and geometric stretching are NOT viscoelastic on this path.
use super::{Error, Mechanics, ModePair, TensionedDisk, drum_spec};
use fs_couple::render::plate::impact::{ImpactSystem, relaxation::{InitialMemory,RelaxationObservation}};
use fs_material::visco::GeneralizedMaxwell;
use fs_phs::RelaxationBranch;
use std::io::Read;

const HEADER:&str="frankensim-head-relaxation-v1";
const MAX_BYTES:usize=65536;
const MAX_MEMORY:usize=256;

#[derive(Clone,Debug)]
pub struct Spec {
    heads:[GeneralizedMaxwell;2],
    initial:InitialMemory,
    band_hz:[f64;2],
}
impl Spec {
    pub fn read(text:&str)->Result<Self,Error> {
        if text.len()>MAX_BYTES {return Err("head relaxation exceeds 64 KiB".into());}
        let (mut header,mut initial,mut band)=(false,None,None);
        let mut heads:[Option<GeneralizedMaxwell>;2]=[None,None];
        for (index,raw) in text.lines().enumerate() {
            let row=raw.split('#').next().unwrap_or("").trim();if row.is_empty() {continue;}
            let bad=|why:&str|->Error {format!("head relaxation line {}: {why}",index+1).into()};
            if !header {if row!=HEADER {return Err(bad("expected frankensim-head-relaxation-v1"));}header=true;continue;}
            let f:Vec<_>=row.split(',').map(str::trim).collect();
            let number=|i:usize|->Result<f64,Error> {
                let x:f64=f[i].parse().map_err(|_|bad("invalid scalar"))?;
                if !x.is_finite() {return Err(bad("nonfinite scalar"));}Ok(x)
            };
            match f[0] {
                "initial" if f.len()==2 && initial.is_none()=>initial=Some(match f[1] {
                    "relaxed"=>InitialMemory::Relaxed,"unrelaxed"=>InitialMemory::Unrelaxed,
                    _=>return Err(bad("initial must be relaxed or unrelaxed")),
                }),
                "band_hz" if f.len()==3 && band.is_none()=> {
                    let b=[number(1)?,number(2)?];
                    if b[0]<0.0 || b[1]<=b[0] {return Err(bad("band must be nonnegative and increasing"));}band=Some(b);
                }
                "head"|"branch"=> {
                    if f.len()!=(if f[0]=="head" {3}else{4}) {return Err(bad("wrong material field count"));}
                    let i=match f[1] {"batter"=>0,"resonant"=>1,_=>return Err(bad("unknown head"))};
                    if f[0]=="head" {
                        if heads[i].is_some() {return Err(bad("duplicate head"));}
                        heads[i]=Some(GeneralizedMaxwell::new(number(2)?,vec![])?);
                    } else {
                        let law=heads[i].as_mut().ok_or_else(||bad("branch must follow its head"))?;
                        if law.terms.len()==8 {return Err(bad("at most eight branches per head"));}
                        let mut terms=law.terms.clone();terms.push((number(2)?,number(3)?));
                        *law=GeneralizedMaxwell::new(law.e_inf,terms)?;
                    }
                }
                _=>return Err(bad("unknown, duplicate or malformed record")),
            }
        }
        let [a,b]=heads;
        let result=Self {heads:[a.ok_or("missing batter material")?,b.ok_or("missing resonant material")?],
            initial:initial.ok_or("missing initial material state")?,band_hz:band.ok_or("missing material band")?};
        for law in &result.heads {
            if !(law.e_inf+law.terms.iter().map(|t|t.0).sum::<f64>()).is_finite() {
                return Err("instantaneous head modulus overflows".into());
            }
        }
        Ok(result)
    }
    pub fn load(path:&str)->Result<Self,Error> {
        let mut text=String::new();std::fs::File::open(path)?.take((MAX_BYTES+1) as u64).read_to_string(&mut text)?;
        Self::read(&text)
    }
    pub fn admit(&self,drum:&drum_spec::Spec)->Result<(),Error> {
        drum.validate()?;
        for (head,law) in drum.heads.iter().zip(&self.heads) {
            if head.young_pa!=law.e_inf {
                return Err("head Young modulus must equal its supplied equilibrium Maxwell modulus; no retuning is inferred".into());
            }
            if head.damping_ratio!=0.0 {
                return Err("head relaxation requires zero head damping ratios in --drum-spec; material loss is not added twice".into());
            }
        }
        if drum.band_hz[0]<self.band_hz[0] || drum.band_hz[1]>self.band_hz[1] {
            return Err("head mode window exceeds the declared material band".into());
        }
        Ok(())
    }
    pub fn attach(&self,system:ImpactSystem,films:&[TensionedDisk],modes:&[Vec<ModePair>],dt:f64)->Result<ImpactSystem,Error> {
        if films.len()!=2 || modes.len()!=2 || !dt.is_finite() || dt<=0.0 {
            return Err("head relaxation requires the two original finite modal pencils and clock".into());
        }
        let count:usize=modes.iter().zip(&self.heads).map(|(m,h)|m.len()*h.terms.iter().filter(|t|t.0>0.0).count()).sum();
        if count>MAX_MEMORY {return Err("head relaxation exceeds 256 material states; no branches or head modes are dropped".into());}
        let mut branches=Vec::with_capacity(count);let mut start=1;
        for ((film,modes),law) in films.iter().zip(modes).zip(&self.heads) {
            if film.spec.young_pa!=law.e_inf || modes.iter().any(|m|!m.lambda.is_finite() || m.lambda<=0.0
                || m.lambda.sqrt()/std::f64::consts::TAU<self.band_hz[0]
                || m.lambda.sqrt()/std::f64::consts::TAU>self.band_hz[1]) {
                return Err("head material does not match the actual equilibrium pencil or frequency band".into());
            }
            let n=modes.len();
            if start+n>system.state().len()/2 {return Err("head material projection exceeds the mechanical basis".into());}
            if law.terms.iter().any(|t|t.0>0.0) {
                let (matrix,factor)=bending_factor(film,modes)?;
                let ratio=law.terms.iter().map(|t|t.0/law.e_inf).sum::<f64>();
                // Gershgorin bounds the instantaneous stiffness, including every
                // off-diagonal bending term. Do not rely on relaxed frequencies
                // alone when adding large instantaneous material stiffness.
                let maximum=(0..n).map(|i|modes[i].lambda+ratio*matrix[i*n..(i+1)*n].iter().map(|v|v.abs()).sum::<f64>())
                    .fold(0.0,f64::max).sqrt();
                if !maximum.is_finite() || maximum*dt>=0.9*std::f64::consts::PI
                    || maximum/std::f64::consts::TAU>self.band_hz[1] {
                    return Err("instantaneous head stiffness exceeds the clock or declared material band".into());
                }
                for &(modulus,tau) in &law.terms {
                    if modulus==0.0 {continue;}
                    for i in 0..n {
                        let mut projection=vec![0.0;system.state().len()];
                        for j in i..n {projection[2*(start+j)]=factor.l(j,i);}
                        branches.push(RelaxationBranch {projection,stiffness:modulus/law.e_inf,relaxation_time_s:tau});
                    }
                }
            }
            start+=n;
        }
        Ok(system.with_relaxation_branches(branches,self.initial,MAX_MEMORY)?)
    }
}

// K_b is assembled with ZERO prestress, not extracted by subtracting nearly
// equal tension-dominated matrices. Keep every modal cross term. L^T q is the
// energy projection because K_b=L L^T. No diagonal damping approximation.
fn bending_factor(film:&TensionedDisk,modes:&[ModePair])->Result<(Vec<f64>,fs_la::factor::Cholesky),Error> {
    if modes.is_empty() || modes.len()>63 || modes.iter().any(|m|m.phi.len()!=film.model.free
        || m.phi.iter().any(|v|!v.is_finite())) {return Err("invalid head bending subspace".into());}
    let rim:Vec<_>=(0..film.mesh.nodes.len()).filter(|&i|film.model.dof_map[3*i].is_none()).collect();
    let bending=fs_plate::assemble(&film.mesh,&film.section,&rim,&[],&fs_plate::AssemblyOptions {
        pretension:0.0,support:fs_plate::EdgeSupport::SimplySupported})?;
    if bending.dof_map!=film.model.dof_map {return Err("bending and prestressed head coordinates differ".into());}
    let n=modes.len();let mut matrix=vec![0.0;n*n];let mut work=vec![0.0;film.model.free];
    for (j,mode) in modes.iter().enumerate() {
        bending.k.spmv(&mode.phi,&mut work);
        for (i,mode) in modes.iter().enumerate() {matrix[i*n+j]=mode.phi.iter().zip(&work).map(|(a,b)|a*b).sum();}
    }
    let scale=matrix.iter().map(|v|v.abs()).fold(0.0,f64::max);
    if !scale.is_finite() || scale<=0.0 || matrix.iter().any(|v|!v.is_finite()) {return Err("invalid projected bending energy".into());}
    for i in 0..n {for j in 0..i {
        if (matrix[i*n+j]-matrix[j*n+i]).abs()>1e-10*scale {return Err("projected bending lost reciprocity".into());}
        let value=f64::midpoint(matrix[i*n+j],matrix[j*n+i]);matrix[i*n+j]=value;matrix[j*n+i]=value;
    }}
    let factor=fs_la::factor::cholesky(&matrix,n).map_err(|e|format!("head bending subspace refused: {e}"))?;
    Ok((matrix,factor))
}

pub fn option(args:&mut Vec<String>)->Result<Option<String>,Error> {
    let mut path=None;let mut i=0;
    while i<args.len() {
        if args[i]!="--head-relaxation" {i+=1;continue;}
        if path.is_some() {return Err("--head-relaxation may be supplied once".into());}
        let value=args.get(i+1).filter(|s|!s.starts_with("--")).ok_or("--head-relaxation needs a material file")?.clone();
        args.drain(i..i+2);path=Some(value);
    }
    Ok(path)
}
pub fn admit_command(selected:bool,command:&str)->Result<(),Error> {
    if selected && !matches!(command,"drum"|"drum-wav"|"drum-mic"|"drum-stretch"|"drum-stretch-wav"|"drum-stretch-mic"|
        "snare"|"snare-wav"|"snare-mic"|"snare-off"|"snare-off-wav"|"snare-off-mic") {
        return Err("head relaxation needs drum[-stretch] or snare[-off], not the linear-only drum-modal image or a cymbal".into());
    }
    Ok(())
}
pub fn observation(system:&Mechanics)->RelaxationObservation {
    match system {
        Mechanics::Reference(s)=>s.relaxation_observation().unwrap_or_default(),
        Mechanics::Nonlinear(s)=>s.relaxation_observation().unwrap_or_default(),
        Mechanics::Substepped(s)=>s.relaxation_observation().unwrap_or_default(),
        Mechanics::Driven {inner,..}=>observation(inner),
        Mechanics::Prepared(_)=>RelaxationObservation::default(),
    }
}

#[cfg(test)]
#[path = "head_relaxation_tests.rs"]
mod tests;
