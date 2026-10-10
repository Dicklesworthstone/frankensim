//! Supplied closed acoustic skins -> existing Helmholtz BEM -> receiver Pa.
//! Exact OBJ labels select moving skin versus acoustically rigid scatterers.
//! No baffle, invented box/lid, material inference, or change to string physics.
//! Multiple disjoint outward closed components are allowed. Self-intersection,
//! component overlap and cavity accessibility remain input responsibilities.
use super::{board_geometry::motion::MotionSurface, linear::Bank};
use fs_bem::{helmholtz::{self, Medium}, panel3d::SpherePanels,
    radiation_policy::GeometryPolicy, near_field::FirstOrder};
use fs_math::c64::C64;
use std::{collections::{BTreeMap,BTreeSet},f64::consts::TAU,
    sync::atomic::{AtomicUsize,Ordering}};

#[path = "exterior_section_skin.rs"]
mod section_skin;

#[path = "exterior_rigid.rs"]
pub mod rigid;

#[path = "exterior_receivers.rs"]
mod receivers;
pub use receivers::ReceiverSet;

pub const RATE:u32=48_000;
pub const MAX_OBJ_BYTES:usize=32*1024*1024;
pub const MAX_SPEC_BYTES:usize=64*1024;
pub const MAX_PANELS:usize=2048;
// Match the existing fs-bem near-field receiver bound; audio remains mono/stereo.
pub const MAX_HARMONIC_RECEIVERS:usize=64;
const MAX_PROJECTION_WORK:usize=3_000_000;
// Each Helmholtz worker owns two dense complex matrices and one LU image.
// Bound their aggregate working set independently of the panel-count gate.
const MAX_FREQUENCY_WORKERS:usize=4;
const MAX_FREQUENCY_DENSE_BYTES:usize=768*1024*1024;

struct FrequencySample {
    pressure:Vec<Vec<C64>>, // input-major, receiver-major
    minimum_ppw:f64,
    maximum_condition_lower_bound:f64,
}

#[derive(Debug)]
pub struct Specification {
    pub source:String,
    pub scale_m:f64,
    pub origin_obj:[f64;3],
    pub offset_m:f64,
    pub receivers:Vec<[f64;3]>,
    /// Explicit triangle-integrated close receivers; omission retains centroid observation.
    pub near_field_receivers:bool,
    /// Optional first-order patterns, indexed by receiver row (zero-based).
    /// Missing entries retain the original omnidirectional pressure receiver.
    pub receiver_patterns:BTreeMap<usize,FirstOrder>,
    pub medium:Medium,
    pub band_hz:(f64,f64),
    pub frequencies:usize,
    pub board_band_hz:f64,
    pub fit_order:usize,
    pub min_ppw:f64,
    pub full_scale_pa:f64,
    rules:BTreeMap<String,bool>,
}
fn number(text:&str)->Result<f64,String> {
    text.parse::<f64>().ok().filter(|x|x.is_finite()).ok_or_else(||"expected finite acoustic scalar".into())
}
impl Specification {
    pub fn read(text:&str)->Result<Self,String> {
        if text.len()>MAX_SPEC_BYTES {return Err("exterior specification exceeds 64 KiB".into());}
        let mut rows=text.lines().map(str::trim).filter(|l|!l.is_empty()&&!l.starts_with('#'));
        if rows.next()!=Some("frankensim-piano-exterior-si-v1") {return Err("expected frankensim-piano-exterior-si-v1".into());}
        let mut fields=BTreeMap::<String,Vec<String>>::new();let mut receivers=Vec::new();let mut rules=BTreeMap::new();
        let mut receiver_patterns=BTreeMap::new();
        for row in rows {
            let f:Vec<_>=row.split(',').map(str::trim).collect();
            match f[0] {
                "moving"|"rigid" if f.len()==2 && !f[1].is_empty()=>{
                    if rules.len()>=256 || rules.insert(f[1].into(),f[0]=="moving").is_some() {return Err("duplicate or excessive acoustic part mapping".into());}
                }
                "receiver-pattern" if f.len()==6=>{
                    let index:usize=f[1].parse().map_err(|_|"receiver-pattern needs a zero-based receiver index")?;
                    let pattern=FirstOrder::new(number(f[2])?,[number(f[3])?,number(f[4])?,number(f[5])?])
                        .map_err(|e|e.to_string())?;
                    if index>=MAX_HARMONIC_RECEIVERS || receiver_patterns.insert(index,pattern).is_some() {
                        return Err(format!("receiver-pattern index must be unique and in 0..{}",MAX_HARMONIC_RECEIVERS));
                    }
                }
                "receiver-m" if f.len()==4=>{
                    if receivers.len()>=MAX_HARMONIC_RECEIVERS {
                        return Err(format!("at most {MAX_HARMONIC_RECEIVERS} harmonic receivers"));
                    }
                    receivers.push([number(f[1])?,number(f[2])?,number(f[3])?]);
                }
                "receiver-evaluation"|"source"|"obj-scale-m"|"obj-origin"|"max-skin-offset-m"|"medium"|
                "band-hz"|"board-band-hz"|"fit-order"|"min-panels-per-wavelength"|"full-scale-pa"=>{
                    if fields.insert(f[0].into(),f[1..].iter().map(|s|(*s).to_owned()).collect()).is_some() {return Err(format!("duplicate acoustic row {}",f[0]));}
                }
                _=>return Err(format!("unknown acoustic row or wrong field count: {}",f[0])),
            }
        }
        let get=|key:&str,n:usize|->Result<&[String],String>{
            fields.get(key).filter(|v|v.len()==n).map(Vec::as_slice).ok_or_else(||format!("missing/wrong-sized {key} row"))
        };
        let scalar=|key:&str|number(&get(key,1)?[0]);
        let source=fields.get("source").ok_or("missing acoustic source attribution")?;
        if source.len()<2 || !["estimated","mixed","published","measured"].contains(&source[0].as_str())
            || source[1..].iter().all(String::is_empty) {return Err("source requires authority and attribution".into());}
        let origin=get("obj-origin",3)?;let gas=get("medium",2)?;let band=get("band-hz",3)?;
        let near_field_receivers=match fields.get("receiver-evaluation").map(Vec::as_slice) {
            None=>false,
            Some([value]) if value=="centroid"=>false,
            Some([value]) if value=="near-field"=>true,
            _=>return Err("receiver-evaluation must be centroid or near-field".into()),
        };
        let mut out=Self {source:source.join(","),scale_m:scalar("obj-scale-m")?,
            origin_obj:[number(&origin[0])?,number(&origin[1])?,number(&origin[2])?],
            offset_m:scalar("max-skin-offset-m")?,receivers,near_field_receivers,receiver_patterns,
            medium:Medium {density:number(&gas[0])?,sound_speed:number(&gas[1])?},
            band_hz:(number(&band[0])?,number(&band[1])?),frequencies:band[2].parse().map_err(|_|"invalid frequency count")?,
            board_band_hz:scalar("board-band-hz")?,fit_order:get("fit-order",1)?[0].parse().map_err(|_|"invalid fit order")?,
            min_ppw:scalar("min-panels-per-wavelength")?,full_scale_pa:scalar("full-scale-pa")?,rules};
        if out.receiver_patterns.keys().any(|i|*i>=out.receivers.len())
            || (!near_field_receivers && out.receiver_patterns.values().any(|p|p.pressure_fraction()<1.)) {
            return Err("directional receivers require receiver-evaluation,near-field and an existing receiver index".into());
        }
        if !(1..=MAX_HARMONIC_RECEIVERS).contains(&out.receivers.len()) || !out.rules.values().any(|moving|*moving)
            || !(1e-9..=1e6).contains(&out.scale_m) || !(0.0..=0.05).contains(&out.offset_m)
            || out.medium.density<=0. || out.medium.sound_speed<=0.
            || out.band_hz.0<10. || out.band_hz.1<=out.band_hz.0 || out.band_hz.1>=0.45*f64::from(RATE)
            || !(17..=257).contains(&out.frequencies) || out.frequencies%2==0
            || out.board_band_hz<=0. || out.board_band_hz>out.band_hz.1
            || !(2..=32).contains(&out.fit_order) || 2*out.fit_order>(out.frequencies+1)/2
            || out.min_ppw<6. || out.full_scale_pa<=0. {
            return Err("invalid explicit exterior units, geometry, medium, frequency/fit or receiver budget".into());
        }
        if near_field_receivers {out.source.push_str("; near-field triangle-integrated receivers, surface-distance propagation");}
        for (index,pattern) in &out.receiver_patterns {
            out.source.push_str(&format!("; receiver {index}: ideal first-order alpha={}, front={:?}; Pa-equivalent on-axis plane-wave sensitivity",
                pattern.pressure_fraction(),pattern.front_axis()));
        }
        Ok(out)
    }
    /// Harmonic arrays never implicitly select a multichannel audio image.
    pub fn require_audio_receivers(&self)->Result<(),String> {
        if !(1..=2).contains(&self.receivers.len()) {
            return Err("played exterior audio requires one or two receivers; use response or admittance for harmonic arrays".into());
        }
        Ok(())
    }
    pub fn omega(&self)->Vec<f64> {(0..self.frequencies).map(|i|
        TAU*(self.band_hz.0+(self.band_hz.1-self.band_hz.0)*i as f64/(self.frequencies-1) as f64)).collect()}
}
fn sub(a:[f64;3],b:[f64;3])->[f64;3] {std::array::from_fn(|i|a[i]-b[i])}
fn dot(a:[f64;3],b:[f64;3])->f64 {a.iter().zip(b).map(|(a,b)|a*b).sum()}
fn cross(a:[f64;3],b:[f64;3])->[f64;3] {[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]}
fn norm(a:[f64;3])->f64 {a.iter().fold(0.0_f64,|r,v|r.hypot(*v))}

pub struct Boundary {
    pub surface:SpherePanels,
    /// Input-major normal velocity per unit modal velocity, in the selected basis.
    pub weights:Vec<Vec<f64>>,
    pub center:[f64;3],
    pub radius:f64,
    pub components:usize,
}
impl Boundary {
    pub fn from_obj(text:&str,spec:&Specification,motion:&MotionSurface)->Result<Self,String> {
        if text.len()>MAX_OBJ_BYTES {return Err("acoustic OBJ exceeds its byte budget".into());}
        let obj=fs_io::obj::read_obj_document(text).map_err(|e|e.to_string())?;
        let count=obj.soup.triangles.len();
        if count>MAX_PANELS {return Err("acoustic mesh exceeds 2048 panels".into());}
        let mut moving=vec![false;count];let mut used=BTreeSet::new();
        for region in &obj.regions {
            let matches:Vec<_>=spec.rules.iter().filter(|(label,_)|region.has_label(label)).collect();
            if matches.len()!=1 {return Err("every OBJ region must match exactly one explicit moving/rigid label".into());}
            let (label,flag)=matches[0];used.insert(label.clone());
            moving[region.triangles.clone()].fill(*flag);
        }
        if used.len()!=spec.rules.len() {return Err("an acoustic part label matches no faces".into());}
        let triangles:Vec<[[f64;3];3]>=obj.soup.triangles.iter().map(|t|t.map(|i| {
            let p=obj.soup.positions[i as usize];let raw=[p.x,p.y,p.z];
            std::array::from_fn(|c|spec.scale_m*(raw[c]-spec.origin_obj[c]))
        })).collect();
        if triangles.iter().flatten().flatten().any(|v|!v.is_finite() || v.abs()>100.) {
            return Err("acoustic skin must be finite and within the 100 metre coordinate budget".into());
        }
        let components=closed_components(&triangles)?;
        let center=std::array::from_fn(|c| {
            let lo=triangles.iter().flatten().map(|p|p[c]).fold(f64::INFINITY,f64::min);
            let hi=triangles.iter().flatten().map(|p|p[c]).fold(f64::NEG_INFINITY,f64::max);
            0.5*(lo+hi)
        });
        let radius=triangles.iter().flatten().map(|&p|norm(sub(p,center))).fold(0.0_f64,f64::max);
        let surface=SpherePanels::from_triangles(triangles.clone()).map_err(|e|e.to_string())?;
        let mut weights=vec![vec![0.;count];motion.shapes.len()];
        let mut projection_work=0usize;
        for (face,tri) in triangles.iter().enumerate() {
            if !moving[face] {continue;}
            let normal=surface.normals()[face];
            // Admit all corners, not only a centroid that could span a hole.
            // Input nonoverlap and topology remain the source mesh's contract.
            for &point in tri {
                let (_, work)=motion.normal_weights_with_work(point,normal,spec.offset_m)?;
                projection_work=projection_work.saturating_add(work);
                if projection_work>MAX_PROJECTION_WORK {return Err("acoustic skin exceeds bounded projection-work budget".into());}
            }
            for bary in [[2./3.,1./6.,1./6.],[1./6.,2./3.,1./6.],[1./6.,1./6.,2./3.]] {
                let point=std::array::from_fn(|c|(0..3).map(|i|bary[i]*tri[i][c]).sum());
                let (row,work)=motion.normal_weights_with_work(point,normal,spec.offset_m)?;
                projection_work=projection_work.saturating_add(work);
                if projection_work>MAX_PROJECTION_WORK {return Err("acoustic skin exceeds bounded projection-work budget".into());}
                for (out,value) in weights.iter_mut().zip(row) {out[face]+=value/3.;}
            }
        }
        Ok(Self {surface,weights,center,radius,components})
    }
    /// Move the skin rows into the very same mass-loaded basis as the bank.
    /// Consume the bare boundary; do not repeat this coordinate transformation.
    pub fn loaded(mut self,bank:&Bank)->Result<Self,String> {
        let panels=self.surface.areas().len();
        let mut transformed=vec![vec![0.;panels];bank.board_count];
        for face in 0..panels {
            let row:Vec<_>=self.weights.iter().map(|m|m[face]).collect();
            for (output,value) in transformed.iter_mut().zip(bank.project_board_shape(&row)?) {
                output[face]=value;
            }
        }
        self.weights=transformed;Ok(self)
    }
    pub fn sample(&self,spec:&Specification)->Result<Samples,String> {
        let plan=ReceiverSet::for_spec(self,spec)?;
        self.sample_grid_plan(&spec.omega(),&spec.receivers,spec.medium,spec.min_ppw,plan)
    }
    fn sample_grid(&self,omega:&[f64],receivers:&[[f64;3]],medium:Medium,min_ppw:f64)
        ->Result<Samples,String> {
        self.sample_grid_mode(omega,receivers,medium,min_ppw,false)
    }
    fn sample_grid_mode(&self,omega:&[f64],receivers:&[[f64;3]],medium:Medium,min_ppw:f64,near:bool)
        ->Result<Samples,String> {
        let plan=ReceiverSet::new(self,receivers,medium,near)?;
        self.sample_grid_plan(omega,receivers,medium,min_ppw,plan)
    }
    fn sample_grid_plan(&self,omega:&[f64],receivers:&[[f64;3]],medium:Medium,min_ppw:f64,plan:ReceiverSet<'_>)
        ->Result<Samples,String> {
        let panels=self.surface.areas().len();
        let dense_bytes=panels.checked_mul(panels).and_then(|n|n.checked_mul(3*std::mem::size_of::<C64>()))
            .ok_or("exterior dense-work size overflow")?;
        let memory_workers=(MAX_FREQUENCY_DENSE_BYTES/dense_bytes.max(1)).max(1);
        let workers=if panels>=1024 {
            std::thread::available_parallelism().map_or(1,|n|n.get())
                .min(MAX_FREQUENCY_WORKERS).min(memory_workers).min(omega.len().max(1))
        } else {1};
        self.sample_grid_plan_with_workers(omega,receivers,medium,min_ppw,plan,workers)
    }
    fn sample_grid_plan_with_workers(&self,omega:&[f64],receivers:&[[f64;3]],medium:Medium,
        min_ppw:f64,plan:ReceiverSet<'_>,workers:usize)->Result<Samples,String> {
        let count=self.weights.len();let panels=self.surface.areas().len();
        if workers==0 || workers>MAX_FREQUENCY_WORKERS || count==0
            || count>super::linear::MAX_BOARD_MODES || panels>MAX_PANELS
            || !(1..=MAX_HARMONIC_RECEIVERS).contains(&receivers.len()) || omega.is_empty() || omega.len()>257
            || omega.iter().enumerate().any(|(i,w)|!w.is_finite() || *w<=0.
                || (i>0 && *w<=omega[i-1]))
            || !medium.density.is_finite() || medium.density<=0.
            || !medium.sound_speed.is_finite() || medium.sound_speed<=0.
            || !min_ppw.is_finite() || min_ppw<6.
            || self.weights.iter().any(|r|r.len()!=panels || r.iter().any(|v|!v.is_finite())) {
            return Err("invalid complete exterior response grid, basis or medium".into());
        }
        // Select from actual closed-component spectra, not the radius of the
        // empty space between the soundboard and disjoint rigid scatterers.
        let policy=GeometryPolicy::new(&self.surface).map_err(|e|e.to_string())?;
        let delays_s=plan.delays_s().to_vec();
        let mut values=vec![vec![Vec::with_capacity(omega.len());count];receivers.len()];
        let mut minimum_ppw=f64::INFINITY;let mut maximum_condition_lower_bound=0.0_f64;
        let sample_at=|w:f64|->Result<FrequencySample,String> {
            let k=w/medium.sound_speed;
            let evaluation=plan.prepare(k)?;
            let fields:Vec<Vec<C64>>=self.weights.iter().map(|r|r.iter().map(|b|C64::new(0.,b/w)).collect()).collect();
            let refs:Vec<&[C64]>=fields.iter().map(Vec::as_slice).collect();
            let formulation=policy.formulation(k).map_err(|e|e.to_string())?;
            let solutions=policy.solve_batch(k,medium,&refs)
                .map_err(|e|format!("exterior solve at {} Hz ({formulation:?}): {e}",w/TAU))?;
            let mut pressure_rows=Vec::with_capacity(count);
            let mut frequency_ppw=f64::INFINITY;
            let mut frequency_condition=0.0_f64;
            for solution in &solutions {
                if !solution.panels_per_wavelength.is_finite() || solution.panels_per_wavelength<min_ppw
                    || !solution.condition_lower_bound.is_finite()
                    || !solution.radiated_power_roundoff_interval.1.is_finite()
                    || solution.radiated_power_roundoff_interval.1<0. {
                    return Err(format!("exterior solve at {} Hz ({formulation:?}) is under-resolved or has inadmissible radiation power/conditioning: ppw={}, condition lower bound={}, power interval={:?}",
                        w/TAU,solution.panels_per_wavelength,solution.condition_lower_bound,
                        solution.radiated_power_roundoff_interval));
                }
                frequency_ppw=frequency_ppw.min(solution.panels_per_wavelength);
                frequency_condition=frequency_condition.max(solution.condition_lower_bound);
                let pressure=evaluation.pressure(solution)?;
                if pressure.len()!=receivers.len() {return Err("exterior receiver pressure count differs".into());}
                for &p in &pressure {
                    if !p.re.is_finite() || !p.im.is_finite() {return Err("exterior receiver pressure is nonfinite".into());}
                }
                pressure_rows.push(pressure);
            }
            if pressure_rows.len()!=count {return Err("exterior modal response count differs".into());}
            Ok(FrequencySample {pressure:pressure_rows,minimum_ppw:frequency_ppw,
                maximum_condition_lower_bound:frequency_condition})
        };
        let samples=if workers==1 {
            omega.iter().copied().map(sample_at).collect::<Result<Vec<_>,_>>()?
        } else {
            let earliest_error=AtomicUsize::new(usize::MAX);
            std::thread::scope(|scope| {
                let mut tasks=Vec::with_capacity(workers);
                for worker in 0..workers {
                    let solve=&sample_at;
                    let earliest_error=&earliest_error;
                    tasks.push(scope.spawn(move || {
                        let mut rows=Vec::new();
                        for i in (worker..omega.len()).step_by(workers) {
                            if i>earliest_error.load(Ordering::Relaxed) {break;}
                            let row=solve(omega[i]);
                            if row.is_err() {earliest_error.fetch_min(i,Ordering::Relaxed);}
                            rows.push((i,row));
                        }
                        rows
                    }));
                }
                let mut ordered:Vec<Option<Result<FrequencySample,String>>>=
                    (0..omega.len()).map(|_|None).collect();
                let mut panicked=false;
                for task in tasks {
                    match task.join() {
                        Ok(rows)=>for (i,row) in rows {ordered[i]=Some(row);},
                        Err(_)=>panicked=true,
                    }
                }
                if panicked {return Err("exterior frequency worker panicked".into());}
                ordered.into_iter().map(|row|row.unwrap_or_else(||Err("exterior frequency result missing".to_string())))
                    .collect::<Result<Vec<_>,String>>()
            })?
        };
        for sample in samples {
            minimum_ppw=minimum_ppw.min(sample.minimum_ppw);
            maximum_condition_lower_bound=maximum_condition_lower_bound.max(sample.maximum_condition_lower_bound);
            for (input,pressure) in sample.pressure.into_iter().enumerate() {
                for (channel,p) in pressure.into_iter().enumerate() {values[channel][input].push(p);}
            }
        }
        Ok(Samples {omega:omega.to_vec(),values,delays_s,minimum_ppw,maximum_condition_lower_bound})
    }
}

/// Receiver-major, mode-major complex pressure per unit generalized
/// acceleration (Pa-equivalent for selected first-order receivers). Time convention exp(-i omega t), no fitted transfer yet.
/// Minimum panels/wavelength is a resolution check, NOT a convergence proof.
pub struct Samples {
    pub omega:Vec<f64>,
    pub values:Vec<Vec<Vec<C64>>>,
    pub delays_s:Vec<f64>,
    pub minimum_ppw:f64,
    pub maximum_condition_lower_bound:f64,
}

// Weld EXACT shared positions only (including +0/-0). Render-export seam
// vertices need not share OBJ indices. No tolerance repair or flipped faces.
fn closed_components(triangles:&[[[f64;3];3]])->Result<usize,String> {
    let mut points=BTreeMap::<[u64;3],usize>::new();let mut indices=Vec::new();
    let mut faces=BTreeSet::new();
    for tri in triangles {
        let mut ids=[0;3];
        for (slot,p) in ids.iter_mut().zip(tri) {
            let key=p.map(|v|if v==0. {0}else{v.to_bits()});let next=points.len();
            *slot=*points.entry(key).or_insert(next);
        }
        let mut sorted=ids;sorted.sort_unstable();
        if sorted[0]==sorted[1] || sorted[1]==sorted[2] || !faces.insert(sorted) {
            return Err("duplicate or degenerate acoustic face".into());
        }
        indices.push(ids);
    }
    let mut edges=BTreeMap::<(usize,usize),Vec<(usize,usize,usize)>>::new();
    for (face,tri) in indices.iter().enumerate() {for i in 0..3 {
        let (a,b)=(tri[i],tri[(i+1)%3]);edges.entry((a.min(b),a.max(b))).or_default().push((face,a,b));
    }}
    let mut adjacency=vec![Vec::new();triangles.len()];
    for uses in edges.values() {
        if uses.len()!=2 || uses[0].1!=uses[1].2 || uses[0].2!=uses[1].1 {
            return Err("acoustic skins must be consistently closed: every edge has two opposite uses".into());
        }
        adjacency[uses[0].0].push(uses[1].0);adjacency[uses[1].0].push(uses[0].0);
    }
    let mut visited=vec![false;triangles.len()];let mut components=0;
    for first in 0..triangles.len() {
        if visited[first] {continue;}
        components+=1;let origin=triangles[first][0];let mut volume6=0.;
        let mut pending=vec![first];visited[first]=true;
        while let Some(face)=pending.pop() {
            let [a,b,c]=triangles[face].map(|p|sub(p,origin));volume6+=dot(a,cross(b,c));
            for &next in &adjacency[face] {if !visited[next] {visited[next]=true;pending.push(next);}}
        }
        if !volume6.is_finite() || volume6<=1e-15 {return Err("every acoustic component must enclose positive outward-oriented volume".into());}
    }
    Ok(components)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    pub fn specification()->String {
        "frankensim-piano-exterior-si-v1\nsource,estimated,authored regression not Steinway\nobj-scale-m,1\nobj-origin,0,0,0\nmax-skin-offset-m,0.02\nmedium,1.2,340\nband-hz,40,400,17\nboard-band-hz,400\nfit-order,4\nmin-panels-per-wavelength,6\nfull-scale-pa,2\nreceiver-m,0.05,0.05,1\nmoving,skin\n".into()
    }
    /// Authored cuboid, not measured instrument geometry. Exact outward faces.
    pub fn box_obj(label:&str,shift:[f64;3],size:[f64;3])->String {
        let mut text=format!("o {label}\n");
        for p in [[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],
                  [0.,0.,1.],[1.,0.,1.],[1.,1.,1.],[0.,1.,1.]] {
            let q:[f64;3]=std::array::from_fn(|c|shift[c]+size[c]*p[c]);
            text.push_str(&format!("v {} {} {}\n",q[0],q[1],q[2]));
        }
        for t in [[1,3,2],[1,4,3],[5,6,7],[5,7,8],[1,2,6],[1,6,5],
                  [2,3,7],[2,7,6],[3,4,8],[3,8,7],[4,1,5],[4,5,8]] {
            // Relative indices let several components be concatenated exactly.
            text.push_str(&format!("f {} {} {}\n",t[0]-9,t[1]-9,t[2]-9));
        }
        text
    }
    fn motion()->MotionSurface {
        let mesh=fs_plate::ShellMesh::new(vec![[0.,0.,0.],[0.1,0.,0.],[0.1,0.1,0.],[0.,0.1,0.]],
            vec![[0,1,2],[0,2,3]]).unwrap();
        MotionSurface::new(mesh,vec![vec![[0.,0.,1.,0.,0.,0.];4]]).unwrap()
    }
    #[test]
    fn both_faces_keep_opposite_normal_motion_and_rigid_lid_is_not_another_source() {
        let obj=box_obj("skin",[0.,0.,-0.01],[0.1,0.1,0.02]);
        let spec=Specification::read(&specification()).unwrap();
        let b=Boundary::from_obj(&obj,&spec,&motion()).unwrap();
        assert_eq!(b.components,1);assert_eq!(b.weights[0][0],-1.);assert_eq!(b.weights[0][2],1.);
        assert!(b.weights[0][4..].iter().all(|v|*v==0.));
        let flux:f64=b.weights[0].iter().zip(b.surface.areas()).map(|(a,b)|a*b).sum();assert!(flux.abs()<1e-15);
        let joined=format!("{}{}",obj,box_obj("lid",[0.,0.,0.04],[0.1,0.1,0.02]));
        let spec=Specification::read(&format!("{}rigid,lid\n",specification())).unwrap();
        let b=Boundary::from_obj(&joined,&spec,&motion()).unwrap();
        assert_eq!(b.components,2);assert!(b.weights[0][12..].iter().all(|v|*v==0.));
    }
    #[test]
    fn open_inward_unmapped_and_remote_moving_geometry_refuse() {
        let obj=box_obj("skin",[0.,0.,-0.01],[0.1,0.1,0.02]);
        let spec=Specification::read(&specification()).unwrap();
        assert!(Boundary::from_obj(&obj.replace("f -8 -6 -7\n",""),&spec,&motion()).is_err());
        assert!(Boundary::from_obj(&obj.replace("o skin","o unknown"),&spec,&motion()).is_err());
        assert!(Boundary::from_obj(&box_obj("skin",[2.,0.,0.],[0.1,0.1,0.02]),&spec,&motion()).is_err());
        let inward=obj.lines().map(|l| {
            if let Some(face)=l.strip_prefix("f ") {let f:Vec<_>=face.split_whitespace().collect();format!("f {} {} {}",f[0],f[2],f[1])}
            else {l.into()}
        }).collect::<Vec<_>>().join("\n");
        assert!(Boundary::from_obj(&inward,&spec,&motion()).is_err());
        assert!(Specification::read(&specification().replace("min-panels-per-wavelength,6","min-panels-per-wavelength,1")).is_err());
    }
    #[test]
    fn harmonic_array_cardinality_and_pattern_indices_are_bounded() {
        let header=specification().replace("receiver-m,0.05,0.05,1\n","");
        let rows=|count|format!("{header}{}","receiver-m,0.05,0.05,1\n".repeat(count));
        for count in [1,2,3,MAX_HARMONIC_RECEIVERS] {
            let text=format!("{}receiver-pattern,{},1,0,0,-1\n",rows(count),count-1);
            let spec=Specification::read(&text).unwrap();
            assert_eq!(spec.receivers.len(),count);
            assert!(spec.receiver_patterns.contains_key(&(count-1)));
            assert_eq!(spec.require_audio_receivers().is_ok(),count<=2);
        }
        assert!(Specification::read(&rows(0)).is_err());
        assert!(Specification::read(&rows(MAX_HARMONIC_RECEIVERS+1)).is_err());
        for suffix in ["receiver-pattern,64,1,0,0,1\n","receiver-m,NaN,0,1\n"] {
            assert!(Specification::read(&format!("{}{suffix}",rows(MAX_HARMONIC_RECEIVERS-1))).is_err());
        }
        assert!(Specification::read(&format!("{}receiver-pattern,3,1,0,0,1\n",rows(3))).is_err());
    }
    #[test]
    fn harmonic_arrays_match_individual_and_permuted_scalar_receivers() {
        let spec=Specification::read(&specification()).unwrap();
        let boundary=Boundary::from_obj(&box_obj("skin",[0.,0.,-0.01],[0.1,0.1,0.02]),&spec,&motion()).unwrap();
        let mut points:Vec<_>=(0..MAX_HARMONIC_RECEIVERS).map(|i|
            [0.02+0.001*i as f64,0.04,if i%2==0 {1.}else {-1.}]).collect();
        points[MAX_HARMONIC_RECEIVERS-1]=points[0]; // coincident observers remain independent rows
        let omega=[TAU*100.,TAU*200.];
        for near in [false,true] {
            let batch=boundary.sample_grid_mode(&omega,&points,spec.medium,6.,near).unwrap();
            assert_eq!(batch.values.len(),MAX_HARMONIC_RECEIVERS);
            let three=boundary.sample_grid_mode(&omega,&points[..3],spec.medium,6.,near).unwrap();
            assert_eq!(three.values,batch.values[..3]);
            for (i,&point) in points.iter().enumerate() {
                let single=boundary.sample_grid_mode(&omega,&[point],spec.medium,6.,near).unwrap();
                assert_eq!(batch.values[i],single.values[0]);
                assert_eq!(batch.delays_s[i],single.delays_s[0]);
                assert_eq!(batch.minimum_ppw,single.minimum_ppw);
                assert_eq!(batch.maximum_condition_lower_bound,single.maximum_condition_lower_bound);
            }
            let reversed:Vec<_>=points.iter().copied().rev().collect();
            let reverse=boundary.sample_grid_mode(&omega,&reversed,spec.medium,6.,near).unwrap();
            for i in 0..points.len() {assert_eq!(batch.values[i],reverse.values[points.len()-1-i]);}
            for invalid in [Vec::new(),vec![points[0];MAX_HARMONIC_RECEIVERS+1],vec![[f64::NAN,0.,1.]]] {
                assert!(boundary.sample_grid_mode(&omega,&invalid,spec.medium,6.,near).is_err());
            }
        }
    }
    #[test]
    fn harmonic_array_near_field_work_budget_is_aggregate() {
        let mut spec=Specification::read(&format!("{}rigid,case\n",specification())).unwrap();
        spec.near_field_receivers=true;
        spec.receivers=vec![[0.05,0.05,2.];MAX_HARMONIC_RECEIVERS];
        let mut obj=box_obj("skin",[0.,0.,-0.01],[0.1,0.1,0.02]);
        for i in 1..33 {obj.push_str(&box_obj("case",[0.2*f64::from(i),0.,-0.01],[0.1,0.1,0.02]));}
        let boundary=Boundary::from_obj(&obj,&spec,&motion()).unwrap();
        // 396 panels * 64 receivers * at least 80 kernel evaluations exceeds
        // the existing 2M cap. A single receiver on this same surface fits.
        let one=ReceiverSet::new(&boundary,&spec.receivers[..1],spec.medium,true).unwrap();
        assert!(one.prepare(1.).is_ok());
        let all=ReceiverSet::for_spec(&boundary,&spec).unwrap();
        let error=all.prepare(1.).err().unwrap();
        assert!(error.contains("kernel-evaluation budget exhausted"),"{error}");
    }
    #[test]
    fn real_bem_observes_both_sides_and_rigid_scattering_changes_pressure() {
        let obj=box_obj("skin",[0.,0.,-0.01],[0.1,0.1,0.02]);
        let spec=Specification::read(&specification()).unwrap();
        let b=Boundary::from_obj(&obj,&spec,&motion()).unwrap();
        let receivers=[[0.05,0.05,1.],[0.05,0.05,-1.]];
        let a=b.sample_grid(&[TAU*100.,TAU*200.],&receivers,spec.medium,6.).unwrap();
        assert!(a.values.iter().flatten().flatten().any(|v|v.abs()>1e-10));
        // Uniform translation's two receivers have opposite, not equal, pressure.
        let p=a.values[0][0][0];let q=a.values[1][0][0];assert!((p+q).abs()<0.02*p.abs());
        let joined=format!("{}{}",obj,box_obj("lid",[0.,0.,0.06],[0.1,0.1,0.02]));
        let spec=Specification::read(&format!("{}rigid,lid\n",specification())).unwrap();
        let lid=Boundary::from_obj(&joined,&spec,&motion()).unwrap();
        let c=lid.sample_grid(&[TAU*100.,TAU*200.],&receivers,spec.medium,6.).unwrap();
        assert!((c.values[0][0][1]-a.values[0][0][1]).abs()>1e-6*a.values[0][0][1].abs());
        assert!(b.sample_grid(&[TAU*100.],&[[0.05,0.05,0.]],spec.medium,6.).is_err());
    }
    #[test]
    fn parallel_frequency_sweep_matches_serial_output_exactly() {
        let spec=Specification::read(&specification()).unwrap();
        let b=Boundary::from_obj(&box_obj("skin",[0.,0.,-0.01],[0.1,0.1,0.02]),&spec,&motion()).unwrap();
        let receivers=[[0.05,0.05,1.],[0.05,0.05,-1.]];
        let omega=[TAU*100.,TAU*200.,TAU*300.];
        let run=|workers| b.sample_grid_plan_with_workers(&omega,&receivers,spec.medium,6.,
            ReceiverSet::new(&b,&receivers,spec.medium,false).unwrap(),workers).unwrap();
        let serial=run(1);let parallel=run(2);
        assert_eq!(serial.values,parallel.values);
        assert_eq!(serial.minimum_ppw,parallel.minimum_ppw);
        assert_eq!(serial.maximum_condition_lower_bound,parallel.maximum_condition_lower_bound);
        let unresolved=[TAU*100.,TAU*50_000.];
        let error=|workers| b.sample_grid_plan_with_workers(&unresolved,&receivers,spec.medium,6.,
            ReceiverSet::new(&b,&receivers,spec.medium,false).unwrap(),workers).err().unwrap();
        assert_eq!(error(1),error(2));
    }
}

#[cfg(test)]
#[path="near_receiver_tests.rs"]
mod near_receiver_tests;

#[cfg(test)]
#[path="directional_receiver_tests.rs"]
mod directional_receiver_tests;
