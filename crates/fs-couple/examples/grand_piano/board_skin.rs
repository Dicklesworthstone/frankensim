//! Closed acoustic boundary of the ADMITTED piecewise-section board.
//!
//! A shallow XY graph is thickened into vertical columns. For each facet the
//! column height is h / n_z, so its two faces are h apart in the facet-normal
//! direction. Thickness jumps retain their exposed step walls; holes and the
//! rim are closed, never capped across the XY plane. No ribs/cabinet/lid are
//! invented. This is the FE section-column image, not a scanned smooth finish.
//!
//! Heights use a declared 1 nm grid solely to avoid sub-roundoff step facets.
//! The actual maximum thickness change and enclosed-volume error are returned.
//! Folded graphs, nonmanifold steps and unresolved panels refuse, not repair.
use super::MotionSurface;
use fs_plate::ShellMesh;
use std::{collections::{BTreeMap, BTreeSet}, fmt::Write};

pub const HEIGHT_QUANTUM_M: f64 = 1e-9;
pub const LABEL: &str = "soundboard_skin";
const MAX_SOURCE_BYTES: usize = 8 * 1024 * 1024;
const MAX_EXPORT_PANELS: usize = 250_000;
const QUADRATURE: [[f64; 3]; 3] = [
    [2. / 3., 1. / 6., 1. / 6.], [1. / 6., 2. / 3., 1. / 6.], [1. / 6., 1. / 6., 2. / 3.],
];
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { std::array::from_fn(|i| a[i] - b[i]) }
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }
fn unit(i: usize) -> [f64; 3] { std::array::from_fn(|j| if i == j { 1. } else { 0. }) }

#[derive(Clone, Debug)]
struct Embedding {
    element: usize,
    // Each boundary vertex's source-facet coordinates and vertical offset.
    bary: [[f64; 3]; 3],
    z: [f64; 3],
}
#[derive(Debug)]
pub struct Skin {
    pub vertices: Vec<[f64; 3]>,
    pub triangles: Vec<[usize; 3]>,
    embeddings: Vec<Embedding>,
    source_triangles: Vec<[usize; 3]>,
    source_nodes: Vec<[f64; 3]>,
    pub section_volume_m3: f64,
    pub volume_m3: f64,
    pub maximum_thickness_change_m: f64,
    pub maximum_offset_m: f64,
}
struct Builder<'a> {
    mesh: &'a ShellMesh,
    skin: Skin,
    ids: BTreeMap<(usize, i64), usize>,
    levels: Vec<BTreeSet<i64>>,
    cap: usize,
}
impl Builder<'_> {
    fn vertex(&mut self, node: usize, level: i64) -> usize {
        *self.ids.entry((node, level)).or_insert_with(|| {
            let mut p = self.mesh.nodes[node]; p[2] += level as f64 * HEIGHT_QUANTUM_M;
            let id = self.skin.vertices.len(); self.skin.vertices.push(p); id
        })
    }
    fn face(&mut self, ids: [usize; 3], embedding: Embedding) -> Result<(), String> {
        if self.skin.triangles.len() == self.cap {
            return Err(format!("section skin exceeds its explicit {}-panel budget; no coarsening or section loss was performed", self.cap));
        }
        let [a,b,c] = ids.map(|i| self.skin.vertices[i]);
        let area2 = dot(cross(sub(b,a),sub(c,a)),cross(sub(b,a),sub(c,a))).sqrt();
        if !area2.is_finite() || area2 < 2e-14 {
            return Err("section skin contains an unresolved tiny step/panel; supply a resolved acoustic mesh instead".into());
        }
        self.skin.triangles.push(ids); self.skin.embeddings.push(embedding); Ok(())
    }
    // Edge a->b follows the outward column's CCW source triangle. The wall's
    // normal points to its right. Split BOTH vertical edges at every incident
    // height, including third-facet levels, to eliminate T-junctions.
    fn wall(&mut self, e: usize, a: usize, b: usize, lo: i64, hi: i64) -> Result<(), String> {
        if lo == hi { return Ok(()); }
        let tri = self.mesh.tris[e];
        let ia = tri.iter().position(|&n| n == a).ok_or("skin edge absent from its element")?;
        let ib = tri.iter().position(|&n| n == b).ok_or("skin edge absent from its element")?;
        let mut ring = vec![(self.vertex(a, lo), 0., lo as f64 * HEIGHT_QUANTUM_M)];
        let up: Vec<_> = self.levels[b].range(lo..=hi).copied().collect();
        for z in up { ring.push((self.vertex(b,z), 1., z as f64 * HEIGHT_QUANTUM_M)); }
        let down: Vec<_> = self.levels[a].range(lo..=hi).rev().copied().collect();
        for z in down { if z != lo { ring.push((self.vertex(a,z), 0., z as f64 * HEIGHT_QUANTUM_M)); } }
        let zc = 0.5 * (lo as f64 + hi as f64) * HEIGHT_QUANTUM_M;
        let center = std::array::from_fn(|c| 0.5 * (self.mesh.nodes[a][c] + self.mesh.nodes[b][c]) + if c == 2 {zc} else {0.});
        let id = self.skin.vertices.len(); self.skin.vertices.push(center);
        let weights = |t: f64| { let mut w = [0.;3]; w[ia]=1.-t; w[ib]=t; w };
        for i in 0..ring.len() {
            let p = ring[i]; let q = ring[(i+1)%ring.len()];
            self.face([id,p.0,q.0], Embedding {element:e,
                bary:[weights(0.5),weights(p.1),weights(q.1)], z:[zc,p.2,q.2]})?;
        }
        Ok(())
    }
}
impl Skin {
    /// Section rows must come from the SAME already-admitted native board used
    /// to prepare `motion`. Exact element IDs/topology are cross-checked; the
    /// loaded coordinates may differ from its unloaded reference. No material
    /// constants, node normals or thicknesses are inferred from a visual asset.
    pub fn from_source(motion: &MotionSurface, source: &str, offset_limit_m: f64, cap: usize) -> Result<Self,String> {
        if source.len() > MAX_SOURCE_BYTES { return Err("section skin source exceeds 8 MiB".into()); }
        let mut rows = source.lines().map(str::trim).filter(|l| !l.is_empty() && !l.starts_with('#'));
        if !matches!(rows.next(), Some("frankensim-board-geometry-si-v1" | "frankensim-crowned-board-si-v1")) {
            return Err("section skin needs its admitted native flat/crowned board source".into());
        }
        let mut thickness = vec![None; motion.mesh.tris.len()];
        for row in rows {
            let f: Vec<_> = row.split(',').map(str::trim).collect();
            if f[0] != "triangle" { continue; }
            if f.len() != 12 { return Err("invalid native section row for skin".into()); }
            let e: usize = f[1].parse().map_err(|_| "invalid skin section ID")?;
            let nodes = [f[2],f[3],f[4]].map(str::parse::<usize>);
            let [a,b,c] = nodes;
            let tri = [a.map_err(|_| "invalid section node")?,b.map_err(|_| "invalid section node")?,c.map_err(|_| "invalid section node")?];
            if motion.mesh.tris.get(e) != Some(&tri) || thickness[e].is_some() {
                return Err("section skin source does not match the prepared board's exact element topology".into());
            }
            thickness[e] = Some(f[5].parse::<f64>().map_err(|_| "invalid skin section thickness")?);
        }
        let thickness = thickness.into_iter().collect::<Option<Vec<_>>>().ok_or("missing skin section thickness")?;
        Self::build(&motion.mesh, &thickness, offset_limit_m, cap)
    }
    /// Construct a column-union boundary. 1 nm height quantization is explicit;
    /// source thickness and structural operators are NEVER changed. No source
    /// facet is removed to meet the caller's panel budget.
    pub fn build(mesh: &ShellMesh, thickness: &[f64], offset_limit_m: f64, cap: usize) -> Result<Self,String> {
        if mesh.nodes.len()>20_000 || mesh.tris.len()>40_000 || mesh.tris.is_empty()
            || thickness.len()!=mesh.tris.len() || !offset_limit_m.is_finite()
            || offset_limit_m<=0. || offset_limit_m>0.05 || !(1..=MAX_EXPORT_PANELS).contains(&cap)
            || mesh.nodes.iter().flatten().any(|v|!v.is_finite() || v.abs()>100.) {
            return Err("invalid complete shallow-board skin, SI offset or panel budget".into());
        }
        let mut levels=vec![BTreeSet::new();mesh.nodes.len()];
        let mut edges=BTreeMap::<(usize,usize),Vec<(usize,usize,usize)>>::new();
        let mut faces=BTreeSet::new(); let mut heights=Vec::new();
        let mut source_volume=0.;let mut rounded_volume=0.;let mut maximum_change=0.0_f64;let mut maximum_offset=0.0_f64;
        for (e,(&h,tri)) in thickness.iter().zip(&mesh.tris).enumerate() {
            let g=mesh.facet(e).map_err(|e|e.to_string())?;
            if !h.is_finite() || h<=0. || g.frame[2][2]<0.95 {
                return Err("section skin requires positive finite thickness and an upward shallow graph".into());
            }
            let ideal=0.5*h/g.frame[2][2];
            if !ideal.is_finite() || ideal>offset_limit_m {return Err("section skin exceeds the declared offset limit".into());}
            let level=(ideal/HEIGHT_QUANTUM_M).round() as i64;
            let height=level as f64*HEIGHT_QUANTUM_M;
            if level<=0 || height>offset_limit_m {return Err("section thickness/offset is unresolved on the 1 nm skin grid".into());}
            heights.push(level);maximum_offset=maximum_offset.max(height);
            maximum_change=maximum_change.max((2.*height*g.frame[2][2]-h).abs());
            source_volume+=g.area_m2*h;rounded_volume+=g.area_m2*2.*height*g.frame[2][2];
            let mut sorted=*tri;sorted.sort_unstable();
            if !faces.insert(sorted) {return Err("duplicate source facet in section skin".into());}
            for i in 0..3 {
                levels[tri[i]].extend([-level,level]);
                let (a,b)=(tri[i],tri[(i+1)%3]);edges.entry((a.min(b),a.max(b))).or_default().push((e,a,b));
            }
        }
        if levels.iter().any(BTreeSet::is_empty) {return Err("unused source node in section skin".into());}
        let mut boundary_degree=vec![0;mesh.nodes.len()];
        for uses in edges.values() {
            if uses.len()==1 {boundary_degree[uses[0].1]+=1;boundary_degree[uses[0].2]+=1;}
            else if uses.len()!=2 || uses[0].1!=uses[1].2 || uses[0].2!=uses[1].1 {
                return Err("nonmanifold/inconsistently oriented source skin edge".into());
            }
        }
        if boundary_degree.iter().any(|&n|n!=0 && n!=2) {return Err("pinched source boundary in section skin".into());}
        let skin=Self {vertices:Vec::new(),triangles:Vec::new(),embeddings:Vec::new(),source_triangles:mesh.tris.clone(),source_nodes:mesh.nodes.clone(),
            section_volume_m3:source_volume,volume_m3:0.,maximum_thickness_change_m:maximum_change,maximum_offset_m:maximum_offset};
        let mut b=Builder {mesh,skin,ids:BTreeMap::new(),levels,cap};
        for (e,&tri) in mesh.tris.iter().enumerate() {
            for sign in [1_i64,-1] {
                let z=sign*heights[e];let order=if sign==1 {[0,1,2]} else {[0,2,1]};
                let ids=order.map(|i|b.vertex(tri[i],z));
                b.face(ids,Embedding {element:e,bary:order.map(unit),z:[z as f64*HEIGHT_QUANTUM_M;3]})?;
            }
        }
        for uses in edges.values() {
            let (e,a,c)=uses[0];
            if uses.len()==1 {b.wall(e,a,c,-heights[e],heights[e])?;}
            else {
                let (f,u,v)=uses[1];
                let (e,a,c,high,low)=if heights[e]>=heights[f] {(e,a,c,heights[e],heights[f])}
                    else {(f,u,v,heights[f],heights[e])};
                b.wall(e,a,c,low,high)?;b.wall(e,a,c,-high,-low)?;
            }
        }
        // Prove exact indexed closure after splitting all level junctions.
        // Four incident walls on a vertical edge (a pinched step) must refuse.
        let mut closure=BTreeMap::<(usize,usize),(usize,i32)>::new();
        for t in &b.skin.triangles {for i in 0..3 {
            let (a,c)=(t[i],t[(i+1)%3]);let use_count=closure.entry((a.min(c),a.max(c))).or_default();
            use_count.0+=1;use_count.1+=if a<c {1}else{-1};
        }}
        if closure.values().any(|v|*v!=(2,0)) {return Err("generated section steps are nonmanifold; no boundary repair was performed".into());}
        let origin=mesh.nodes[0];
        let volume=b.skin.triangles.iter().map(|t| {
            let [a,c,d]=t.map(|i|sub(b.skin.vertices[i],origin));dot(a,cross(c,d))/6.
        }).sum::<f64>();
        if !volume.is_finite() || !source_volume.is_finite() || source_volume<=0.
            || (volume-rounded_volume).abs()>1e-10*rounded_volume+1e-15 {
            return Err("section skin did not close the column volume".into());
        }
        b.skin.volume_m3=volume;Ok(b.skin)
    }
    pub fn panel_triangles(&self)->Vec<[[f64;3];3]> {
        self.triangles.iter().map(|t|t.map(|i|self.vertices[i])).collect()
    }
    /// Exact positive three-point quadrature for P1 translations plus the
    /// bilinear rotation/offset field. Known section sites are retained; no
    /// nearest-facet search, vertex snap or extrapolation across a hole occurs.
    pub fn normal_weights(&self,motion:&MotionSurface,normals:&[[f64;3]])->Result<Vec<Vec<f64>>,String> {
        if motion.mesh.tris!=self.source_triangles || motion.mesh.nodes!=self.source_nodes || normals.len()!=self.triangles.len()
            || normals.iter().any(|n|n.iter().any(|x|!x.is_finite()) || (dot(*n,*n)-1.).abs()>1e-8) {
            return Err("skin normal/motion topology mismatch".into());
        }
        let mut out=vec![vec![0.;self.triangles.len()];motion.shapes.len()];
        for (f,site) in self.embeddings.iter().enumerate() {
            let tri=motion.mesh.tris[site.element];
            for q in QUADRATURE {
                let bary:[f64;3]=std::array::from_fn(|i|(0..3).map(|j|q[j]*site.bary[j][i]).sum());
                let z=(0..3).map(|j|q[j]*site.z[j]).sum::<f64>();
                for (result,shape) in out.iter_mut().zip(&motion.shapes) {
                    let state:[f64;6]=std::array::from_fn(|c|(0..3).map(|i|bary[i]*shape[tri[i]][c]).sum());
                    let velocity=[state[0]+z*state[4],state[1]-z*state[3],state[2]];
                    result[f]+=dot(normals[f],velocity)/3.;
                }
            }
        }
        if out.iter().flatten().any(|v|!v.is_finite()) {return Err("section-skin velocity overflow".into());}
        Ok(out)
    }
    pub fn obj(&self)->String {
        let mut out=format!("# Admitted section-column acoustic image, SI metres; not a scanned cabinet.\n# height quantum {} m; maximum section thickness change {} m\n# section volume {} m3; skin volume {} m3\no {LABEL}\n",HEIGHT_QUANTUM_M,self.maximum_thickness_change_m,self.section_volume_m3,self.volume_m3);
        for p in &self.vertices {writeln!(out,"v {:.17e} {:.17e} {:.17e}",p[0],p[1],p[2]).unwrap();}
        for t in &self.triangles {writeln!(out,"f {} {} {}",t[0]+1,t[1]+1,t[2]+1).unwrap();}
        out
    }
}

#[cfg(test)]
#[path="board_skin_tests.rs"]
mod tests;
