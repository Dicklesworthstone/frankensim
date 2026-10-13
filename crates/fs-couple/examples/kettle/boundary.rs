//! Actual head, rigid bearing annulus and independently supplied closed outer wall.
//! The pressure basis keeps its INNER wall; this mesh is the EXTERIOR boundary.
use super::{Error, model::rim, spec::Spec};
use fs_plate::{ModePair, shell::head::TensionedDisk};
use std::collections::BTreeMap;

pub struct Boundary {
    pub triangles: Vec<[[f64;3];3]>,
    /// Positive OUTWARD velocity per mass-normalized head coordinate.
    pub weights: Vec<Vec<f64>>,
    pub radius_m: f64,
    pub minimum_box_width_m: f64,
}
impl Boundary {
    pub fn new(spec: &Spec, film: &TensionedDisk, modes: &[ModePair]) -> Result<Self, Error> {
        spec.validate_exterior()?;
        if film.spec.radius_m != spec.head.radius_m || modes.is_empty() || modes.len() > 32
            || modes.iter().any(|m| m.phi.len()!=film.model.free || m.phi.iter().any(|p| !p.is_finite()))
        { return Err("exterior requires the same complete head geometry and basis".into()); }
        let edges = rim(film);
        let mut rim_nodes: Vec<_> = edges.iter().flat_map(|&(a,b)| [a,b]).collect();
        rim_nodes.sort_unstable(); rim_nodes.dedup();
        if edges.len()!=rim_nodes.len() || edges.len()<8 { return Err("expected one circular head boundary".into()); }
        let panels = film.mesh.tris.len()+2*edges.len()
            +(2*(spec.exterior.len()-2)+1)*edges.len();
        if panels > 2048 { return Err("closed kettle exterior exceeds 2048 panels; no wall station is dropped".into()); }
        let z = spec.profile.last().ok_or("missing head plane")?.z_m;
        let mut vertices: Vec<_> = film.mesh.nodes.iter().map(|&(x,y)| [x,y,z]).collect();
        let mut rings = Vec::new();
        for station in spec.exterior[1..].iter().rev() {
            let mut ring = BTreeMap::new();
            for &i in &rim_nodes {
                let (x,y) = film.mesh.nodes[i]; let r = x.hypot(y);
                if (r-film.spec.radius_m).abs()>1e-10*film.spec.radius_m { return Err("outer wall requires the actual circular head seam".into()); }
                ring.insert(i,vertices.len());
                vertices.push([station.radius_m*x/r,station.radius_m*y/r,station.z_m]);
            }
            rings.push(ring);
        }
        let bottom = vertices.len(); vertices.push([0.0,0.0,spec.exterior[0].z_m]);
        let mut indices = film.mesh.tris.clone();
        for &(a,b) in &edges {
            let (ta,tb)=(rings[0][&a],rings[0][&b]);
            indices.extend([[a,ta,tb],[a,tb,b]]); // rigid bearing annulus
            for pair in rings.windows(2) {
                let (ta,tb,ba,bb)=(pair[0][&a],pair[0][&b],pair[1][&a],pair[1][&b]);
                indices.extend([[ta,ba,bb],[ta,bb,tb]]);
            }
            let last=rings.last().ok_or("missing bottom ring")?;
            indices.push([last[&a],bottom,last[&b]]);
        }
        let mut uses=BTreeMap::<(usize,usize),Vec<(usize,usize)>>::new();
        for t in &indices { for j in 0..3 {
            let (a,b)=(t[j],t[(j+1)%3]);
            uses.entry((a.min(b),a.max(b))).or_default().push((a,b));
        }}
        if indices.len()!=panels || uses.values().any(|e|e.len()!=2 || e[0]!=(e[1].1,e[1].0)) {
            return Err("outer kettle surface is not an oriented closed manifold".into());
        }
        let triangles: Vec<_> = indices.iter().map(|t|t.map(|i|vertices[i])).collect();
        let mut low=[f64::INFINITY;3]; let mut high=[f64::NEG_INFINITY;3];
        let mut radius_m=0.0_f64;
        for p in &vertices {
            for j in 0..3 { low[j]=low[j].min(p[j]); high[j]=high[j].max(p[j]); }
            radius_m=radius_m.max(p[0].hypot(p[1]).hypot(p[2]));
        }
        let minimum_box_width_m=(0..3).map(|j|(high[j]-low[j]).next_up()).fold(f64::INFINITY,f64::min);
        if !radius_m.is_finite() || radius_m<=0.0 || !minimum_box_width_m.is_finite() || minimum_box_width_m<=0.0 {
            return Err("unrepresentable acoustic bounds".into());
        }
        let mut weights=vec![vec![0.0;panels];modes.len()];
        for (row,mode) in weights.iter_mut().zip(modes) {
            for (i,tri) in film.mesh.tris.iter().enumerate() {
                row[i]=-tri.iter().map(|&node|film.model.dof_map[3*node].map_or(0.0,|j|mode.phi[j])/3.0).sum::<f64>();
            }
        }
        Ok(Self {triangles,weights,radius_m:radius_m.next_up(),minimum_box_width_m})
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_bem::panel3d::SpherePanels;
    use fs_plate::shell::profile::ProfileBudget;
    fn fixture() -> (Spec,TensionedDisk,Vec<ModePair>) {
        let mut s=Spec::read(include_str!("illustrative-parabolic.fsk")).unwrap();
        s.head.radial_intervals=2; s.head.azimuths=16;
        let f=TensionedDisk::new(s.head,ProfileBudget {
            max_nodes:100,max_triangles:200,max_feature_evaluations:0 }).unwrap();
        let mut phi=vec![0.0;f.model.free];
        for (i,&(x,y)) in f.mesh.nodes.iter().enumerate() {
            if let Some(j)=f.model.dof_map[3*i] {phi[j]=1.0-(x*x+y*y)/s.head.radius_m.powi(2);}
        }
        // Synthetic kinematic fixture, NOT an eigenpair or acoustic measurement.
        let mode=ModePair {lambda:1.0,phi,residual:0.0,interval:(1.0,1.0)};
        (s,f,vec![mode])
    }
    #[test]
    fn closed_outer_volume_and_head_flux_use_their_correct_surfaces() {
        let (s,f,m)=fixture();let b=Boundary::new(&s,&f,&m).unwrap();
        let surface=SpherePanels::from_triangles(b.triangles).unwrap();
        let volume:f64=surface.centroids().iter().zip(surface.normals()).zip(surface.areas())
            .map(|((p,n),a)|p.iter().zip(n).map(|(p,n)|p*n).sum::<f64>()*a/3.0).sum();
        let frustum:f64=s.exterior.windows(2).map(|p|core::f64::consts::PI/3.0*(p[1].z_m-p[0].z_m)
            *(p[0].radius_m.powi(2)+p[0].radius_m*p[1].radius_m+p[1].radius_m.powi(2))).sum();
        let angle=core::f64::consts::TAU/s.head.azimuths as f64;
        assert!((volume-frustum*angle.sin()/angle).abs()<1e-12);
        let flux:f64=b.weights[0].iter().zip(surface.areas()).map(|(b,a)|b*a).sum();
        assert!((flux+f.modal_area(&m[0].phi).unwrap()).abs()<1e-12);
        assert!(b.weights[0][f.mesh.tris.len()..].iter().all(|v|*v==0.0));
    }
    #[test]
    fn outer_wall_may_not_cut_inside_at_its_own_unsampled_inner_knot() {
        let (mut s,_,_)=fixture();
        use fs_couple::render::plate::impact::cavity::meridian::Station;
        s.exterior.insert(5,Station {z_m:0.09,radius_m:0.001});
        assert!(s.validate_exterior().is_err());
    }
    #[test]
    fn omitted_outer_wall_is_legal_for_mechanics_but_not_audio() {
        let (mut s,f,m)=fixture();s.exterior.clear();
        s.validate().unwrap();
        assert!(Boundary::new(&s,&f,&m).is_err());
    }
}
