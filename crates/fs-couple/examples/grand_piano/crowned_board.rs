//! Measured/authored crowned midsurface -> existing stiffened CST/DKT shell
//! -> mass-normalized piano bridge ports. No fitted frequencies or new stepper.
//!
//! Crown is the supplied REFERENCE shape. Explicit downbearing rows solve a
//! nonlinear static equilibrium before computing small-vibration modes.
//! The board is a shallow graph above the XY chart, with clamped rim. Ribs and
//! bridges retain their existing E,G,A,Iy,J,e,rho; Iz follows the explicitly
//! declared rectangular section reconstruction from A and Iy. Perfect bonds.
//! Grain directions are projected into each actual facet before assembly.
//!
//! Acoustics remains an explicitly PROJECTED flat-baffle approximation:
//! signed normal displacement * true area is preserved at each quadrature
//! station but its source position is moved to z=0. Not exterior BEM. This
//! does not flatten the STRUCTURAL mesh or suppress in-plane shell motion.
//! See Mamou-Mani et al., JASA 123 (2008), doi:10.1121/1.2836787 for the
//! distinction between initial crown and a downbearing/prestress calculation.
use super::board_geometry::{BoardGeometry, PreparedBoard, SurfaceSample};
use super::linear::{BoardMode, MAX_BOARD_MODES};
use fs_plate::{PlateSection, ShellMesh, ShellModel, ShellSupport};
use fs_plate::shell::stiffened::{BeamSection, ShellBeam, assemble_stiffened_shell};
use std::{collections::{BTreeMap, BTreeSet}, f64::consts::TAU, fmt::Write};

#[path = "downbearing.rs"]
mod downbearing;

pub const HEADER: &str = "frankensim-crowned-board-si-v1";
const SHAPE: &str = "beam-section,rectangular-from-area-inertia";
const MAX_BYTES: usize = 8*1024*1024;
const MAX_HEIGHT_M: f64 = 0.05;

fn number(x:&str)->Result<f64,String> {
    x.parse::<f64>().ok().filter(|x|x.is_finite()).ok_or_else(||"invalid finite shell scalar".into())
}
fn index(x:&str)->Result<usize,String> {x.parse().map_err(|_|"invalid shell index".into())}
fn dot(a:[f64;3],b:[f64;3])->f64 {a.iter().zip(b).map(|(a,b)|a*b).sum()}
fn cross(a:[f64;3],b:[f64;3])->[f64;3] {
    [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
}
fn unit(a:[f64;3])->Result<[f64;3],String> {
    let l=dot(a,a).sqrt();
    if !l.is_finite() || l<1e-12 {return Err("unresolved shell direction".into());}
    Ok(a.map(|x|x/l))
}
fn meaningful(text:&str)->impl Iterator<Item=&str> {
    text.lines().map(str::trim).filter(|s|!s.is_empty()&&!s.starts_with('#'))
}
pub fn is_crowned(text:&str)->bool {meaningful(text).next()==Some(HEADER)}

/// Retain a native board's full sections, topology, supports and bridge map;
/// replace ONLY its reference midsurface heights. All heights are SI and are
/// required, not interpolated or filled. The shape reconstruction is explicit.
pub fn elevate(flat:&str, heights:&[f64], source:&str)->Result<String,String> {
    if flat.len()>MAX_BYTES || source.len()>8192 || source.contains(['\n','\r']) {
        return Err("crown input/source exceeds the format budget".into());
    }
    BoardGeometry::read(flat)?;
    let mut out=format!("{HEADER}\n{SHAPE}\n");
    let mut count=0;
    for row in meaningful(flat).skip(1) {
        let f:Vec<_>=row.split(',').map(str::trim).collect();
        if f[0]=="node" {
            let id=index(f[1])?;
            let z=*heights.get(id).ok_or("crown height missing for a panel node")?;
            if !z.is_finite() {return Err("nonfinite crown height".into());}
            writeln!(out,"{row},{z:.17e}").unwrap();count+=1;
        } else if f[0]=="source" {
            writeln!(out,"source,mixed,{}; crown: {source}",f[1..].join(",")).unwrap();
        } else {writeln!(out,"{row}").unwrap();}
    }
    if count!=heights.len() {return Err("crown heights must cover exactly the panel nodes".into());}
    CrownedBoard::read(&out)?;
    Ok(out)
}

#[derive(Clone)]
struct Site {key:u8, tri:usize, weights:[f64;3], arm:[f64;3]}
pub struct CrownedBoard {
    mesh:ShellMesh, sections:Vec<PlateSection>, beams:Vec<ShellBeam>, fixed:Vec<usize>,
    sites:Vec<Site>, damping:f64, source:String,
    preload:Option<downbearing::Specification>,
    pub max_height_m:f64,
    pub area_m2:f64,
    pub mass_kg:f64,
}
impl CrownedBoard {
    /// Native FSB row semantics, except node,id,x,y,z, the explicit rectangular
    /// beam declaration, and optional bridge_arm,key,dx,dy,dz [m]. An arm moves
    /// the vertical force port from the panel midsurface to its bearing point.
    /// No arm means the original midsurface port, NOT an inferred bridge top.
    /// Clamped supports and zero authored membrane pretension only. Optional
    /// downbearing rows require an unloaded reference and complete course loads;
    /// the full-coordinate static solve precedes tangent-modal analysis.
    pub fn read(text:&str)->Result<Self,String> {
        if text.len()>MAX_BYTES || !is_crowned(text) {return Err(format!("expected bounded {HEADER}"));}
        let preload=downbearing::Specification::read(text)?;
        let mut flat=String::from("frankensim-board-geometry-si-v1\n");
        let mut nodes=Vec::new();let mut arms=BTreeMap::new();let mut shape=false;
        for row in meaningful(text).skip(1) {
            let f:Vec<_>=row.split(',').map(str::trim).collect();
            match f[0] {
                "downbearing" | "downbearing-source" | "preload-reference" => {},
                "node" => {
                    if f.len()!=5 || nodes.len()>=20_000 || index(f[1])?!=nodes.len() {
                        return Err("crowned node rows must be contiguous bounded x,y,z coordinates".into());
                    }
                    let p=[number(f[2])?,number(f[3])?,number(f[4])?];
                    if p[2].abs()>MAX_HEIGHT_M {return Err("crown exceeds the 50 mm shallow-board height budget".into());}
                    nodes.push(p);writeln!(flat,"{}",f[..4].join(",")).unwrap();
                }
                "beam-section" => {
                    if row!=SHAPE || shape {return Err("declare rectangular beam sections exactly once".into());}
                    shape=true;
                }
                "bridge_arm" => {
                    if f.len()!=5 {return Err("bridge_arm needs key,dx,dy,dz in metres".into());}
                    let key:u8=f[1].parse().map_err(|_|"invalid bridge arm key")?;
                    let arm=[number(f[2])?,number(f[3])?,number(f[4])?];
                    if !(21..=108).contains(&key) || arm.iter().any(|x|x.abs()>0.2)
                        || arms.insert(key,arm).is_some() {return Err("invalid/duplicate bridge arm".into());}
                }
                _=>{writeln!(flat,"{row}").unwrap();}
            }
        }
        if !shape {return Err("missing rectangular beam-section declaration".into());}
        // The existing front door owns all projected topology, section, index,
        // source, support, beam-duplication and bridge-barycentric admission.
        BoardGeometry::read(&flat)?;
        let mut triangles=Vec::new();let mut cards=Vec::new();let mut fixed=Vec::new();
        let mut raw_beams=Vec::new();let mut sites=Vec::new();let mut source=String::new();
        let mut damping=0.;
        for row in meaningful(&flat).skip(1) {
            let f:Vec<_>=row.split(',').map(str::trim).collect();
            match f[0] {
                "triangle"=>{
                    triangles.push([index(f[2])?,index(f[3])?,index(f[4])?]);
                    let mut card=[0.;7];for i in 0..7 {card[i]=number(f[5+i])?;}cards.push(card);
                }
                "fixed"=>fixed.push(index(f[1])?),
                "support" if f[1]!="clamped"=>return Err("crowned shells require explicit clamped supports; plate simply-supported is not a shell pin".into()),
                "pretension" if number(f[1])?!=0.=>return Err("authored membrane pretension is not downbearing; supply explicit bridge loads instead".into()),
                "damping"=>damping=number(f[1])?,
                "source"=>source=f[1..].join(","),
                "stiffener"=>{
                    let mut card=[0.;7];for i in 0..7 {card[i]=number(f[i+1])?;}
                    let path=f[8..].iter().map(|x|index(x)).collect::<Result<Vec<_>,_>>()?;
                    raw_beams.push((card,path));
                }
                "bridge"=>{
                    let key=f[1].parse::<u8>().map_err(|_|"invalid bridge key")?;
                    sites.push(Site {key,tri:index(f[2])?,weights:[number(f[3])?,number(f[4])?,number(f[5])?],
                        arm:arms.remove(&key).unwrap_or([0.;3])});
                }
                _=>{}
            }
        }
        if !arms.is_empty() {return Err("bridge arm has no matching bridge station".into());}
        if let Some(spec)=&preload {spec.check_sites(&sites)?;}
        let mesh=ShellMesh::new(nodes,triangles).map_err(|e|e.to_string())?;
        let mut normals=vec![[0.;3];mesh.nodes.len()];let mut sections=Vec::new();
        let mut area_m2=0.;let mut mass_kg=0.;
        for (i,card) in cards.iter().enumerate() {
            let g=mesh.facet(i).map_err(|e|e.to_string())?;
            if g.frame[2][2]<0.95 {return Err("soundboard must be an upward shallow graph (normal.z >= 0.95)".into());}
            for &node in &mesh.tris[i] {for d in 0..3 {normals[node][d]+=g.area_m2*g.frame[2][d];}}
            // The authored angle is in the board XY chart, not in whichever
            // direction the first edge of this triangle happens to point.
            let grain=[card[6].cos(),card[6].sin(),0.];
            let angle=dot(grain,g.frame[1]).atan2(dot(grain,g.frame[0]));
            let s=PlateSection::orthotropic_plane_stress_at_angle(card[2],card[3],card[4],card[5],
                card[0],card[1],angle).map_err(|e|e.to_string())?;
            area_m2+=g.area_m2;mass_kg+=s.density*s.thickness*g.area_m2;sections.push(s);
        }
        for normal in &mut normals {*normal=unit(*normal)?;}
        let mut beams=Vec::new();
        for (c,path) in raw_beams {
            // c = E,G,A,Iy,J,e,rho. A and Iy determine the explicitly declared
            // rectangle, so the second bending plane gets Iz, not a copy of Iy.
            let height=(12.*c[3]/c[2]).sqrt();let width=c[2]/height;
            let iz=c[2]*width*width/12.;
            if !height.is_finite() || height<=0. || !iz.is_finite() || iz<=0. || c[4]<=0. {
                return Err("rectangular shell beams require positive finite Iy,Iz,J".into());
            }
            for pair in path.windows(2) {
                let offsets=[normals[pair[0]].map(|x|c[5]*x),normals[pair[1]].map(|x|c[5]*x)];
                let direction=unit(std::array::from_fn(|i|normals[pair[0]][i]+normals[pair[1]][i]))?;
                let chord:[f64;3]=std::array::from_fn(|i|mesh.nodes[pair[1]][i]+offsets[1][i]
                    -mesh.nodes[pair[0]][i]-offsets[0][i]);
                mass_kg+=c[6]*c[2]*dot(chord,chord).sqrt();
                beams.push(ShellBeam {nodes:[pair[0],pair[1]],offsets_m:offsets,section_z:direction,
                    section:BeamSection {young_pa:c[0],shear_pa:c[1],density_kg_m3:c[6],area_m2:c[2],
                        iy_m4:c[3],iz_m4:iz,torsion_m4:c[4]}});
            }
        }
        if !mass_kg.is_finite() || mass_kg<=0. {return Err("crowned board mass overflow".into());}
        let max_height_m=mesh.nodes.iter().map(|p|p[2].abs()).fold(0.,f64::max);
        Ok(Self {mesh,sections,beams,fixed,sites,damping,source,preload,max_height_m,area_m2,mass_kg})
    }
    pub fn prepare(&self,keys:&[u8],upper_hz:f64)->Result<PreparedBoard,String> {
        if !upper_hz.is_finite() || upper_hz<=0. || upper_hz>80_000. || keys.is_empty() {
            return Err("invalid crowned soundboard frequency/key budget".into());
        }
        let mut seen=BTreeSet::new();
        for &key in keys {
            if !seen.insert(key) || !self.sites.iter().any(|s|s.key==key) {
                return Err(format!("missing or duplicated crowned bridge key {key}"));
            }
        }
        let (model,acoustic_mesh,equilibrium)=downbearing::prepare(self)?;
        let report=fs_plate::modes_shell(&model,(0.,(TAU*upper_hz).powi(2)),&fs_plate::SliceOptions::default())
            .map_err(|e|e.to_string())?;
        if report.below_low!=0 || report.modes.is_empty() || report.modes.len()>MAX_BOARD_MODES {
            return Err(format!("crowned shell has {} modes in band and {} below zero; admitted complete stable slice is 1..={MAX_BOARD_MODES}",report.modes.len(),report.below_low));
        }
        let mut modes=Vec::new();let mut intervals=Vec::new();let mut mphi=vec![0.;model.free];
        for (i,pair) in report.modes.iter().enumerate() {
            if pair.phi.len()!=model.free || pair.phi.iter().any(|x|!x.is_finite())
                || !pair.lambda.is_finite() || pair.lambda<=0. || pair.interval.0<=0.
                || !pair.interval.0.is_finite() || !pair.interval.1.is_finite()
                || pair.interval.1<pair.interval.0 {
                return Err("unresolved crowned shell eigenpair".into());
            }
            model.m.spmv(&pair.phi,&mut mphi);
            let product=|phi:&[f64]|phi.iter().zip(&mphi).map(|(a,b)|a*b).sum::<f64>();
            if (product(&pair.phi)-1.).abs()>1e-7 || !product(&pair.phi).is_finite()
                || report.modes[..i].iter().any(|p| {let v=product(&p.phi);!v.is_finite()||v.abs()>1e-7}) {
                return Err("crowned modes are not mass orthonormal".into());
            }
            let mut bridge=[0.;88];
            for site in &self.sites {
                let mut displacement=0.;
                for (n,&node) in self.mesh.tris[site.tri].iter().enumerate() {
                    let u=nodal(&model,&pair.phi,node,0);let theta=nodal(&model,&pair.phi,node,3);
                    displacement+=site.weights[n]*(u[2]+cross(theta,site.arm)[2]);
                }
                bridge[usize::from(site.key-21)]=displacement;
            }
            if bridge.iter().any(|x|!x.is_finite()) {return Err("crowned bridge projection overflow".into());}
            modes.push(BoardMode {frequency_hz:pair.lambda.sqrt()/TAU,damping_ratio:self.damping,bridge,volume:0.});
            intervals.push((pair.interval.0.sqrt()/TAU,pair.interval.1.sqrt()/TAU));
        }
        let mut surface=Vec::with_capacity(self.mesh.tris.len()*3);
        let mut surface_area=0.;
        for (element,tri) in acoustic_mesh.tris.iter().enumerate() {
            let g=acoustic_mesh.facet(element).map_err(|e|e.to_string())?;
            surface_area+=g.area_m2;
            let n=g.frame[2];let projected=g.area_m2*n[2]/3.;
            for weights in [[2./3.,1./6.,1./6.],[1./6.,2./3.,1./6.],[1./6.,1./6.,2./3.]] {
                let position=[(0..3).map(|i|weights[i]*acoustic_mesh.nodes[tri[i]][0]).sum(),
                    (0..3).map(|i|weights[i]*acoustic_mesh.nodes[tri[i]][1]).sum(),0.];
                let shape:Vec<f64>=report.modes.iter().map(|pair|(0..3).map(|i|
                    weights[i]*dot(n,nodal(&model,&pair.phi,tri[i],0))/n[2]).sum()).collect();
                if shape.iter().any(|x|!x.is_finite()) {return Err("crowned radiation projection overflow".into());}
                for (mode,&value) in modes.iter_mut().zip(&shape) {mode.volume+=projected*value;}
                surface.push(SurfaceSample {position_m:position,area_m2:projected,mode_shape:shape});
            }
        }
        Ok(PreparedBoard {modes,surface,area_m2:surface_area,mass_kg:self.mass_kg,
            provenance:format!("{}; 3-D CST/DKT crowned shell, {} eccentric rectangular beam segments; reference max |z|={} m; projected flat-baffle radiation; {}",self.source,self.beams.len(),self.max_height_m,equilibrium),
            frequency_intervals_hz:intervals,free_dofs:model.free})
    }
}
fn nodal(model:&ShellModel,phi:&[f64],node:usize,start:usize)->[f64;3] {
    std::array::from_fn(|i|model.dof_map[6*node+start+i].map_or(0.,|d|phi[d]))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture(height:f64)->String {
        let mut f=String::from("frankensim-board-geometry-si-v1\nsource,estimated,regression only\nsupport,clamped\npretension,0\ndamping,0.01\n");
        for (i,p) in [[0.,0.],[1.,0.],[1.,1.],[0.,1.],[0.5,0.5]].iter().enumerate() {
            writeln!(f,"node,{i},{},{}",p[0],p[1]).unwrap();
        }
        for (i,t) in [[0,1,4],[1,2,4],[2,3,4],[3,0,4]].iter().enumerate() {
            writeln!(f,"triangle,{i},{},{},{},0.008,450,1e10,8e8,0.3,6e8,0.27",t[0],t[1],t[2]).unwrap();
        }
        f.push_str("fixed,0\nfixed,1\nfixed,2\nfixed,3\nbridge,69,0,0,0,1\n");
        elevate(&f,&[0.,0.,0.,0.,height],"authored test crown").unwrap()
    }
    #[test]
    fn crown_changes_the_computed_structural_operator_not_an_audio_filter() {
        let a=CrownedBoard::read(&fixture(0.)).unwrap();
        let b=CrownedBoard::read(&fixture(0.015)).unwrap();
        let assemble=|b:&CrownedBoard|assemble_stiffened_shell(&b.mesh,&b.sections,&b.fixed,ShellSupport::Clamped,&b.beams).unwrap();
        let flat=assemble(&a);let curved=assemble(&b);
        assert!(b.mass_kg>a.mass_kg);assert_ne!(flat.k,curved.k);
        let z=flat.dof_map[6*4+2].unwrap();
        assert!(curved.k.get(z,z)>flat.k.get(z,z));
        let flat=a.prepare(&[69],400.).unwrap();let curved=b.prepare(&[69],400.).unwrap();
        assert!((curved.modes[0].frequency_hz-flat.modes[0].frequency_hz).abs()>0.1);
        assert!(curved.modes[0].bridge[48].abs()>1e-10);
        for (i,m) in curved.modes.iter().enumerate() {
            let volume:f64=curved.surface.iter().map(|p|p.area_m2*p.mode_shape[i]).sum();
            assert_eq!(volume,m.volume);
        }
    }
    #[test]
    fn real_model_d_sections_and_all_ribs_bridges_survive_crown_lifting() {
        let preset=super::super::steinway_d::build(4).unwrap();
        let heights:Vec<_>=meaningful(&preset.geometry).filter(|l|l.starts_with("node,")).map(|_|0.005).collect();
        let text=elevate(&preset.geometry,&heights,"5 mm offset test, NOT measured crown").unwrap();
        let b=CrownedBoard::read(&text).unwrap();
        assert_eq!(b.sites.len(),88);assert!(b.beams.len()>100);
        assert_eq!(b.sections.len(),b.mesh.tris.len());
        assert!(b.beams.iter().any(|b|b.offsets_m[0][2]<0.));
        assert!(b.beams.iter().any(|b|b.offsets_m[0][2]>0.));
    }
    #[test]
    fn unsupported_crown_preload_and_missing_data_refuse() {
        let good=fixture(0.01);
        for bad in [good.replace(SHAPE,""),good.replace("pretension,0","pretension,100"),
            good.replace("support,clamped","support,simply_supported"),
            good.replace("node,4,0.5,0.5,","node,4,0.5,0.5,NaN"),
            format!("{good}bridge_arm,60,0,0,0.03\n")] {assert!(CrownedBoard::read(&bad).is_err());}
        let board=CrownedBoard::read(&good).unwrap();assert!(board.prepare(&[60],400.).is_err());
    }
}
