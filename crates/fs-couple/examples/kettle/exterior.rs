//! Closed, oriented zero-thickness rigid bowl plus its actual moving film.
//! The meridian is the supplied inner wall; using it outside is an EXPLICIT
//! zero-wall-thickness approximation, not an inferred copper skin or thickness.
use super::Error;
use fs_couple::render::plate::impact::cavity::meridian::Station;
use fs_exec::CancelGate;
use fs_plate::{ModePair, shell::head::TensionedDisk};
use std::collections::BTreeMap;

pub struct Boundary {
    /// Translation only: all returned triangles are relative to this point.
    pub origin_m: [f64;3],
    pub triangles: Vec<[[f64;3];3]>,
    /// Outward normal velocity per unit positive-DOWNWARD modal velocity.
    /// P0 panel means of the original P1 film. The bowl wall has exact zeros.
    pub weights: Vec<Vec<f64>>,
    pub head_panels: usize,
    pub radius_m: f64,
}
impl Boundary {
    pub fn new(film: &TensionedDisk, modes: &[ModePair], profile: &[Station],
        max_panels: usize, max_terms: usize, gate: &CancelGate) -> Result<Self,Error>
    {
        if gate.is_requested() { return Err("exterior preparation cancelled".into()); }
        if profile.len()<3 || profile.len()>64 || profile[0].radius_m!=0.0
            || profile.iter().any(|s|!s.z_m.is_finite() || !s.radius_m.is_finite() || s.radius_m<0.0)
            || profile[1..].iter().any(|s|s.radius_m<=0.0)
            || profile.windows(2).any(|s|s[0].z_m>=s[1].z_m)
            || !film.spec.radius_m.is_finite() || film.spec.radius_m<=0.0
            || profile[profile.len()-1].radius_m!=film.spec.radius_m
            || modes.is_empty() || modes.len()>32 || film.mesh.tris.is_empty()
            || film.mesh.tris.len()>max_panels
            || film.mesh.nodes.len().checked_mul(3)!=Some(film.model.dof_map.len())
            || film.model.dof_map.iter().flatten().any(|&i|i>=film.model.free)
            || film.mesh.nodes.iter().any(|&(x,y)|!x.is_finite() || !y.is_finite())
            || film.mesh.tris.iter().flatten().any(|&i|i>=film.mesh.nodes.len())
            || modes.iter().any(|m|m.phi.len()!=film.model.free || m.phi.iter().any(|v|!v.is_finite()))
        { return Err("exterior requires the matching finite film and a sealed axis-bottom wall".into()); }
        let mut edges = BTreeMap::<(usize,usize),Vec<(usize,usize)>>::new();
        for tri in &film.mesh.tris { for i in 0..3 {
            let (a,b)=(tri[i],tri[(i+1)%3]);
            edges.entry((a.min(b),a.max(b))).or_default().push((a,b));
        }}
        let mut successor=BTreeMap::new();
        for occurrences in edges.values() {
            match occurrences.as_slice() {
                &[(a,b)] => {if successor.insert(a,b).is_some() {return Err("branched head rim".into());}},
                &[(a,b),(c,d)] if a==d && b==c => {},
                _ => return Err("head is not an oriented manifold".into()),
            }
        }
        let start=*successor.keys().next().ok_or("head has no boundary")?;
        let mut ring=Vec::with_capacity(successor.len()); let mut at=start;
        loop {
            if ring.contains(&at) {return Err("head boundary is not a single simple cycle".into());}
            ring.push(at);
            at=*successor.get(&at).ok_or("open head boundary")?;
            if at==start {break;}
        }
        if ring.len()!=successor.len() || ring.len()<3 {return Err("head has several boundary components".into());}
        let count=film.mesh.tris.len().checked_add(ring.len()*(2*profile.len()-3))
            .ok_or("exterior panel count overflow")?;
        if count>max_panels || count.checked_mul(modes.len()).is_none_or(|n|n>max_terms) {
            return Err("exterior panel/shape budget exceeded; no head modes or wall stations removed".into());
        }
        for &i in &ring {
            let (x,y)=film.mesh.nodes[i];
            if (x.hypot(y)-film.spec.radius_m).abs()>64.0*f64::EPSILON*film.spec.radius_m {
                return Err("actual head rim does not match the circular wall".into());
            }
        }
        let top=profile[profile.len()-1].z_m;
        let middle=profile[0].z_m/2.0+top/2.0;
        let origin_m=[0.0,0.0,middle];
        let mut vertices:Vec<_>=film.mesh.nodes.iter().map(|&(x,y)|[x,y,top-middle]).collect();
        let mut indexed=film.mesh.tris.clone();
        let axis=vertices.len();vertices.push([0.0,0.0,profile[0].z_m-middle]);
        let mut lower=vec![axis];
        for (station_index,station) in profile.iter().enumerate().skip(1) {
            if gate.is_requested() {return Err("exterior preparation cancelled".into());}
            let upper=if station_index+1==profile.len() {ring.clone()} else {
                ring.iter().map(|&i| {
                    let (x,y)=film.mesh.nodes[i];let id=vertices.len();
                    vertices.push([x/film.spec.radius_m*station.radius_m,
                        y/film.spec.radius_m*station.radius_m,station.z_m-middle]);id
                }).collect::<Vec<_>>()
            };
            for i in 0..ring.len() {
                let next=(i+1)%ring.len();
                if lower.len()==1 {indexed.push([lower[0],upper[next],upper[i]]);}
                else {indexed.extend([[lower[i],lower[next],upper[next]],[lower[i],upper[next],upper[i]]]);}
            }
            lower=upper;
        }
        // Check exact shared identities, not approximate proximity welding.
        let mut incidence=BTreeMap::<(usize,usize),(usize,i32)>::new();
        for tri in &indexed {for i in 0..3 {
            let (a,b)=(tri[i],tri[(i+1)%3]);let e=incidence.entry((a.min(b),a.max(b))).or_default();
            e.0+=1;e.1+=if a<b {1}else{-1};
        }}
        if incidence.values().any(|&(n,s)|n!=2 || s!=0) {return Err("exterior is not closed and consistently oriented".into());}
        let triangles:Vec<_>=indexed.iter().map(|t|t.map(|i|vertices[i])).collect();
        let head_panels=film.mesh.tris.len();
        let mut weights=vec![vec![0.0;count];modes.len()];
        for (p,tri) in film.mesh.tris.iter().enumerate() {
            let t=tri.map(|i|film.mesh.nodes[i]);
            let area2=(t[1].0-t[0].0)*(t[2].1-t[0].1)-(t[2].0-t[0].0)*(t[1].1-t[0].1);
            if !area2.is_finite() || area2<=0.0 {return Err("reversed or degenerate moving head".into());}
            for (j,m) in modes.iter().enumerate() {
                weights[j][p]=-tri.iter().map(|&i|film.model.dof_map[3*i].map_or(0.0,|d|m.phi[d])/3.0).sum::<f64>();
            }
        }
        let radius_m=vertices.iter().map(|p|p[0].hypot(p[1]).hypot(p[2])).fold(0.0_f64,f64::max).next_up();
        if !radius_m.is_finite() || radius_m<=0.0 || triangles.iter().flatten().flatten().any(|v|!v.is_finite())
            || weights.iter().flatten().any(|v|!v.is_finite()) {return Err("exterior geometry or velocity overflow".into());}
        Ok(Self {origin_m,triangles,weights,head_panels,radius_m})
    }
}

#[cfg(test)]
#[path="exterior_tests.rs"]
mod tests;
