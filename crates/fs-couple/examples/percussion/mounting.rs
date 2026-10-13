//! File-selected, finite-inertia backing for the shell's distributed felt.
//! Cold geometry composition only; the existing impact/material/time owners
//! retain every patch, support inertia, reaction and accepted-step history.
use super::{Error, Obstacle, playing};
use fs_couple::render::plate::impact::felt::{FeltPad, KelvinBranch};
use fs_couple::render::plate::impact::felt::mount::{Carrier, CarrierPort, MountSite, MountedFelt};
use fs_material::fiber::WoolFelt;
use fs_plate::shell::{survey::MeshShell, reduction::ShellReduction};
use fs_plate::shell::reduction::radiation::ShellFace;
use std::io::{Read, Write};

const MAX_BYTES:u64=65_536;
#[derive(Clone, Copy, Debug)]
struct Backing { z_m:f64, thickness_m:f64, prior_strain:f64 }
#[derive(Clone, Copy, Debug)]
struct Site { xy_m:[f64;2], area_m2:f64 }
#[derive(Clone, Debug)]
pub struct Spec {
    carrier:Carrier,
    upper:Backing,
    lower:Backing,
    law:WoolFelt,
    /// Whole-washer coefficients, scaled by each declared loaded area fraction.
    creep:Vec<KelvinBranch>,
    sites:Vec<Site>,
}

impl Spec {
    pub fn parse(text:&str)->Result<Self,Error> {
        if text.len() as u64>MAX_BYTES {return Err("felt mount specification exceeds 64 KiB".into());}
        let mut rows=text.lines().enumerate().filter_map(|(i,line)| {
            let line=line.split('#').next().unwrap_or("").trim();
            (!line.is_empty()).then_some((i+1,line))
        });
        if rows.next().map(|(_,line)|line)!=Some("frankensim-felt-mount-v1") {
            return Err("expected frankensim-felt-mount-v1".into());
        }
        let (mut carrier,mut upper,mut lower,mut law)=(None,None,None,None);
        let mut creep=Vec::new();let mut sites=Vec::new();
        for (line,row) in rows {
            let fields:Vec<_>=row.split(',').map(str::trim).collect();
            let bad=||format!("invalid felt-mount record at line {line}");
            let value=|i:usize|->Result<f64,Error> {
                let v=fields.get(i).ok_or_else(bad)?.parse::<f64>()?;
                if !v.is_finite() {return Err(bad().into());}Ok(v)
            };
            match (fields[0],fields.len()) {
                ("carrier",6) if carrier.is_none()=>carrier=Some(Carrier {
                    mass_kg:value(1)?,stiffness_n_m:value(2)?,resistance_n_s_m:value(3)?,
                    initial_position_m:value(4)?,initial_velocity_m_s:value(5)?,
                }),
                ("upper",4) if upper.is_none()=>upper=Some(Backing {
                    z_m:value(1)?,thickness_m:value(2)?,prior_strain:value(3)?,
                }),
                ("lower",4) if lower.is_none()=>lower=Some(Backing {
                    z_m:value(1)?,thickness_m:value(2)?,prior_strain:value(3)?,
                }),
                ("felt",7) if law.is_none()=>law=Some(WoolFelt::new(value(1)?,value(2)?,
                    value(3)?,value(4)?,value(5)?,value(6)?)?),
                ("creep",3) if creep.len()<4=>creep.push(KelvinBranch {
                    stiffness_n_m:value(1)?,viscosity_n_s_m:value(2)?,
                }),
                ("site",4) if sites.len()<8=>sites.push(Site {xy_m:[value(1)?,value(2)?],area_m2:value(3)?}),
                _=>return Err(bad().into()),
            }
        }
        let result=Self {
            carrier:carrier.ok_or("felt mount needs explicit carrier mass, spring, resistance and initial motion")?,
            upper:upper.ok_or("felt mount needs an upper backing plane, thickness and conditioning")?,
            lower:lower.ok_or("felt mount needs a lower backing plane, thickness and conditioning")?,
            law:law.ok_or("felt mount needs an explicit material law")?,creep,sites,
        };
        if result.upper.z_m<=result.lower.z_m || result.sites.is_empty()
            || result.sites.iter().any(|s|s.area_m2<=0.0) {
            return Err("felt mount needs ordered backings and 1..=8 positive-area sites per face".into());
        }
        for (i,site) in result.sites.iter().enumerate() {
            if result.sites[..i].iter().any(|previous|previous.xy_m==site.xy_m) {
                return Err("felt mount sites must be distinct; do not count a loaded area twice".into());
            }
        }
        let area=result.area()?;
        // Reuse the existing material/carrier admission before any modal solve.
        // These are validation rows, not fabricated physical geometry.
        let mut sites=Vec::new();
        for site in &result.sites {
            for (backing,sign) in [(result.upper,1.0),(result.lower,-1.0)] {
                sites.push(MountSite {pad:result.pad(*site,backing,0.0,vec![sign],area)?,
                    compression_per_carrier_m:-sign});
            }
        }
        MountedFelt::new(1,result.carrier,&sites)?;
        Ok(result)
    }
    pub fn load(path:&str)->Result<Self,Error> {
        let mut text=String::new();std::fs::File::open(path)?.take(MAX_BYTES+1).read_to_string(&mut text)?;
        Self::parse(&text)
    }
    fn area(&self)->Result<f64,Error> {
        let area=self.sites.iter().map(|s|s.area_m2).sum::<f64>();
        if !area.is_finite() || area<=0.0 {return Err("felt loaded area is not representable".into());}
        Ok(area)
    }
    fn pad(&self,site:Site,backing:Backing,precompression_m:f64,weights:Vec<f64>,area:f64)
        ->Result<FeltPad,Error>
    {
        let fraction=site.area_m2/area;
        if !fraction.is_finite() || fraction<=0.0 {return Err("felt area fraction is not representable".into());}
        Ok(FeltPad {area_m2:site.area_m2,thickness_m:backing.thickness_m,precompression_m,weights,
            law:self.law.clone(),prior_maximum_strain:backing.prior_strain,
            creep:self.creep.iter().map(|b|KelvinBranch {
                stiffness_n_m:b.stiffness_n_m*fraction,viscosity_n_s_m:b.viscosity_n_s_m*fraction,
            }).collect()})
    }
    /// Append the carrier after the ENTIRE original structural layout. The shell
    /// keeps its existing prefix and all other bodies are zero in the felt rows.
    /// Each site is vertical-area quadrature, not a claimed continuous footprint
    /// or a surveyed pressure map. Neither point positions nor preload are snapped.
    pub fn compile(&self,shell:&MeshShell,reduction:&ShellReduction,shell_offset:usize,
        structural_modes:usize)->Result<MountedFelt,Error>
    {
        let shell_end=shell_offset.checked_add(reduction.mode_count()).ok_or("mount shell address overflow")?;
        if shell_end>structural_modes || structural_modes==0
            || structural_modes>=fs_couple::render::plate::impact::MAX_IMPACT_MODES {
            return Err("mount shell does not fit a bounded original body layout".into());
        }
        let area=self.area()?;
        let positive=reduction.surface_positions(&shell.nodal_thickness_m,ShellFace::Positive)?;
        let negative=reduction.surface_positions(&shell.nodal_thickness_m,ShellFace::Negative)?;
        let mut sites=Vec::with_capacity(2*self.sites.len());
        for site in &self.sites {
            let skin=|nodes:&[[f64;3]],side:ShellFace|->Result<_,Error> {
                let (face,bary)=playing::shell_location(nodes,&shell.mesh.tris,site.xy_m)?;
                Ok(reduction.surface_point_port(&shell.nodal_thickness_m,face,bary,side,[0.,0.,1.])?)
            };
            let p=skin(&positive,ShellFace::Positive)?;let n=skin(&negative,ShellFace::Negative)?;
            // Surface orientation is not an instruction to put the upper washer
            // under the shell. Order actual skin heights at the same XY station.
            let (top,bottom)=if p.position_m[2]>n.position_m[2] {(p,n)}else{(n,p)};
            if !(top.position_m[2]>bottom.position_m[2]) {
                return Err("mount station needs two distinct finite shell skins".into());
            }
            for (port,backing,sign) in [(top,self.upper,1.0),(bottom,self.lower,-1.0)] {
                let gap=sign*(backing.z_m-port.position_m[2]);
                let initial_gap=gap+sign*self.carrier.initial_position_m;
                let compression=backing.thickness_m-gap;
                let initial_compression=backing.thickness_m-initial_gap;
                if !gap.is_finite() || !initial_gap.is_finite() || initial_gap<=0.0
                    || !compression.is_finite() || !initial_compression.is_finite()
                    || initial_compression/backing.thickness_m>self.law.eps_densify {
                    return Err("mount backing intersects the shell or initial felt exceeds densification; geometry is not clamped".into());
                }
                let mut weights=vec![0.0;structural_modes];
                for (target,source) in weights[shell_offset..shell_end].iter_mut().zip(&port.weights) {
                    *target=sign*source;
                }
                sites.push(MountSite {pad:self.pad(*site,backing,compression,weights,area)?,
                    compression_per_carrier_m:-sign});
            }
        }
        let mount=MountedFelt::new(structural_modes,self.carrier,&sites)?;
        eprintln!("felt mount: carrier_coordinate={}, mass_kg={}, ground_K_N_m={}, ground_R_Ns_m={}, loaded_area_per_face_m2={area}, sites_per_face={}; actual finite-thickness skin gaps; input is not measured/calibrated authority; no gravity equilibrium or full stand flexure",
            mount.port.coordinate,self.carrier.mass_kg,self.carrier.stiffness_n_m,
            self.carrier.resistance_n_s_m,self.sites.len());
        eprintln!("felt mount zero-coordinate precompressions_m={:?}",mount.pads.iter().map(|p|p.precompression_m).collect::<Vec<_>>());
        Ok(mount)
    }
    /// The host appends the carrier last, BEFORE private Kelvin/material memory
    /// and later acoustic states. Capture this port before exterior preparation.
    pub fn port(&self,structural_modes:usize)->Result<CarrierPort,Error> {
        Ok(CarrierPort {coordinate:structural_modes.checked_sub(1).ok_or("missing mount coordinate")?,
            inverse_sqrt_mass:1.0/self.carrier.mass_kg.sqrt()})
    }
}

/// Zero-extend each contact row, not the flattened matrix's tail. Preserve the
/// existing gaps, quadrature, force law, provenance and internal-loss coefficient.
/// This is cold coordinate embedding, not another contact law.
pub fn extend_contact(contact:&Obstacle,old:usize,new:usize)->Result<Obstacle,Error> {
    if old==0 || new<old || new>fs_couple::render::plate::impact::MAX_IMPACT_MODES
        || contact.n_points()>4096
        || old.checked_mul(contact.n_points())!=Some(contact.collocation().len()) {
        return Err("invalid mounted contact embedding".into());
    }
    let size=new.checked_mul(contact.n_points()).ok_or("mount contact extent overflow")?;
    let mut rows=vec![0.0;size];
    for (source,target) in contact.collocation().chunks_exact(old).zip(rows.chunks_exact_mut(new)) {
        target[..old].copy_from_slice(source);
    }
    Ok(Obstacle::new(rows,contact.n_points(),new,contact.gaps().to_vec(),contact.weights().to_vec(),
        contact.stiffness(),contact.alpha(),contact.provenance().to_string())?
        .with_internal_loss(contact.internal_loss())?)
}

pub fn option(args:&mut Vec<String>)->Result<Option<String>,Error> {
    let mut positions=args.iter().enumerate().filter(|(_,a)|a.as_str()=="--mount-spec");
    let Some((at,_))=positions.next() else {return Ok(None);};
    if positions.next().is_some() {return Err("--mount-spec may be supplied once".into());}
    let path=args.get(at+1).ok_or("--mount-spec needs a file path")?;
    if path.is_empty() || path.starts_with("--") {return Err("--mount-spec needs a file path, not another option".into());}
    let path=path.clone();args.drain(at..at+2);Ok(Some(path))
}
pub fn admit_command(selected:bool,command:&str)->Result<(),Error> {
    if selected && !matches!(command,"splash"|"splash-wav"|"splash-mic") {
        return Err("--mount-spec currently requires a single supplied/reference shell (splash[-wav|-mic]); paired hi-hat and drum supports are not inferred".into());
    }
    Ok(())
}
pub fn header(out:&mut impl Write)->Result<(),Error> {
    write!(out,",mount_position_m,mount_velocity_m_s")?;Ok(())
}
pub fn row(port:CarrierPort,state:&[f64],out:&mut impl Write)->Result<(),Error> {
    let [x,v]=port.motion(state)?;write!(out,",{x:.17e},{v:.17e}")?;Ok(())
}

#[cfg(test)]
#[path="mounting_tests.rs"]
mod tests;
