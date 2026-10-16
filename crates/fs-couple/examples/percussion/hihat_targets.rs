//! Exterior hand targets in the pair's world frame. Both source meshes retain
//! their original, outward-positive local coordinates.
use super::{Error, Stroke, Shell, mallets};
use fs_couple::render::plate::impact::compliant::striker::FeltStriker;

#[derive(Clone,Copy,Debug,Default,Eq,PartialEq)]
pub(super) enum Target { #[default] Upper, Lower }
impl Target {
    pub fn local_point(self,world:[f64;2])->[f64;2] {
        match self {Self::Upper=>world,Self::Lower=>[world[0],-world[1]]}
    }
    pub fn local_stroke(self,stroke:Stroke)->Stroke {
        Stroke{position_m:stroke.position_m.map(|p|self.local_point(p)),..stroke}
    }
    pub fn label(self)->&'static str {
        match self {Self::Upper=>"upper",Self::Lower=>"lower"}
    }
    #[allow(clippy::too_many_arguments)]
    pub fn mallet(self,spec:&mallets::Spec,shell:&Shell,stroke:Stroke,coordinate:usize,
        first_shell:usize,total:usize)->Result<FeltStriker,Error>
    {
        let mut tip=spec.compile_shell(&shell.reduction,&shell.mesh,self.local_stroke(stroke),
            coordinate,first_shell,total)?;
        // The proper lower-body rotation reverses world Y. Restore the world
        // +X,+Y,-X,-Y site order before bind_shaft applies its supplied world-XY
        // azimuth and reciprocal footprint moments. Each gap stays with its row.
        if self==Self::Lower {tip.pads.swap(1,3);}
        Ok(tip)
    }
}

#[derive(Clone,Copy,Debug,Default,Eq,PartialEq)]
pub(super) struct Selection {pub first:Target,pub second:Target}
impl Selection {
    pub fn options(args:&mut Vec<String>,second:Option<Stroke>)->Result<Self,Error> {
        let mut selected=[None,None];let mut i=0;
        while i<args.len() {
            let hand=match args[i].as_str() {
                "--strike-target"=>0,"--second-stick-target"=>1,_=>{i+=1;continue;}
            };
            if selected[hand].is_some() {return Err("duplicate hi-hat hand target".into());}
            let target=match args.get(i+1).map(String::as_str) {
                Some("upper")=>Target::Upper,Some("lower")=>Target::Lower,
                _=>return Err("hi-hat strike target must be upper or lower; only the exterior skin is supported".into()),
            };
            selected[hand]=Some(target);args.drain(i..i+2);
        }
        if selected[1].is_some() && second.is_none() {
            return Err("--second-stick-target requires --second-stick-position-m X Y".into());
        }
        Ok(Self{first:selected[0].unwrap_or_default(),second:selected[1].unwrap_or_default()})
    }
    pub fn admit(self,second:Option<Stroke>)->Result<(),Error> {
        if self.second==Target::Lower && second.is_none() {
            return Err("a lower second-hand target requires its own strike position".into());
        }
        Ok(())
    }
}
