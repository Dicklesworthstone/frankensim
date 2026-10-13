//! Same closed BEM owner, but with the declared outer profile, not a cylinder.
use super::{Boundary, Error, MAX_INPUTS, MAX_PANELS, check_closed};
use crate::{chamber, ModePair, TensionedDisk};
use std::collections::BTreeMap;

impl Boundary {
    pub fn profiled_heads(films: &[TensionedDisk], modes: &[Vec<ModePair>],
        profile: &chamber::Spec) -> Result<Self, Error> {
        profile.validate()?;
        if films.len()!=2 || modes.len()!=2 || films[0].mesh.nodes!=films[1].mesh.nodes
            || films[0].mesh.tris!=films[1].mesh.tris {
            return Err("profiled exterior needs matching head meshes and mode pencils".into());
        }
        let mesh=&films[0].mesh; let n=mesh.nodes.len();
        let last=profile.stations.len()-1; let depth=profile.stations[last].z_m;
        let count=modes.iter().map(Vec::len).sum::<usize>();
        if count==0 || count>MAX_INPUTS || modes.iter().zip(films).any(|(ms,f)|
            ms.is_empty() || ms.iter().any(|m|m.phi.len()!=f.model.free || m.phi.iter().any(|v|!v.is_finite())))
            || [0,last].into_iter().enumerate().any(|(head,i)|profile.stations[i].radius_m!=films[head].spec.radius_m) {
            return Err("profiled exterior needs actual matching head radii and finite modes within its source budget".into());
        }
        let mut uses=BTreeMap::<(usize,usize),Vec<(usize,usize)>>::new();
        for t in &mesh.tris {for i in 0..3 {
            let (a,b)=(t[i],t[(i+1)%3]);uses.entry((a.min(b),a.max(b))).or_default().push((a,b));
        }}
        let edges:Vec<_>=uses.values().filter(|e|e.len()==1).map(|e|e[0]).collect();
        if edges.len()<3 {return Err("profiled head has no boundary loop".into());}
        let mut rim:Vec<_>=edges.iter().flat_map(|&(a,b)|[a,b]).collect();rim.sort_unstable();rim.dedup();
        let panels=2*mesh.tris.len()+edges.len()*(4+2*last);
        if panels>MAX_PANELS {return Err("profiled exterior exceeds its panel budget; no profile stations are discarded".into());}
        let mut points:Vec<_>=mesh.nodes.iter().map(|&(x,y)|[x,y,0.5*depth])
            .chain(mesh.nodes.iter().map(|&(x,y)|[x,y,-0.5*depth])).collect();
        let mut rings=Vec::with_capacity(last+1);
        for (station,&radius) in profile.stations.iter().zip(&profile.outer_radii_m) {
            let mut ring=BTreeMap::new();
            for &i in &rim {
                let (x,y)=mesh.nodes[i];let r=x.hypot(y);
                if (r-films[0].spec.radius_m).abs()>1e-10*films[0].spec.radius_m {
                    return Err("profiled exterior requires one circular head rim".into());
                }
                ring.insert(i,points.len());points.push([radius*x/r,radius*y/r,0.5*depth-station.z_m]);
            }
            rings.push(ring);
        }
        let mut indices=Vec::with_capacity(panels);
        indices.extend(mesh.tris.iter().copied());
        indices.extend(mesh.tris.iter().map(|t|[t[2]+n,t[1]+n,t[0]+n]));
        for &(a,b) in &edges {
            let (ta,tb)=(rings[0][&a],rings[0][&b]);
            indices.extend([[a,ta,tb],[a,tb,b]]);
            let (ba,bb)=(rings[last][&a],rings[last][&b]);
            indices.extend([[a+n,bb,ba],[a+n,b+n,bb]]);
            for level in 0..last {
                let (ta,tb,ba,bb)=(rings[level][&a],rings[level][&b],rings[level+1][&a],rings[level+1][&b]);
                indices.extend([[ta,ba,bb],[ta,bb,tb]]);
            }
        }
        check_closed(&indices)?;
        let triangles=indices.iter().map(|t|t.map(|i|points[i])).collect();
        let mut weights=vec![vec![0.0;panels];count];let mut at=0;
        for (head,(film,ms)) in films.iter().zip(modes).enumerate() {
            for mode in ms {
                for (f,t) in mesh.tris.iter().enumerate() {
                    let shape=t.iter().map(|&i|film.model.dof_map[3*i].map_or(0.0,|d|mode.phi[d])/3.0).sum::<f64>();
                    weights[at][head*mesh.tris.len()+f]=if head==0 {-shape}else{shape};
                }
                at+=1;
            }
        }
        Ok(Self {triangles,weights,state_modes:(1..=count).collect()})
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_profile_closes_outward_and_only_the_original_heads_radiate() {
        let drum=crate::drum_spec::Spec {radial_intervals:2,azimuths:16,..crate::drum_spec::Spec::reference()};
        let (films,modes)=drum.prepare(2e-6,false).unwrap();
        let profile=chamber::Spec::read(include_str!("../estimated-waisted-chamber.fscp")).unwrap();
        let boundary=Boundary::profiled_heads(&films,&modes,&profile).unwrap();
        let mut volume=0.0;
        for t in &boundary.triangles {
            let cross=[t[1][1]*t[2][2]-t[1][2]*t[2][1],t[1][2]*t[2][0]-t[1][0]*t[2][2],t[1][0]*t[2][1]-t[1][1]*t[2][0]];
            volume+=t[0].iter().zip(cross).map(|(a,b)|a*b).sum::<f64>()/6.0;
        }
        let polygon_factor=16.0*(core::f64::consts::TAU/16.0).sin()/core::f64::consts::TAU;
        let exact: f64=profile.stations.windows(2).zip(profile.outer_radii_m.windows(2))
            .map(|(s,r)|core::f64::consts::PI/3.0*(s[1].z_m-s[0].z_m)*(r[0]*r[0]+r[0]*r[1]+r[1]*r[1])).sum();
        assert!(volume>0.0 && (volume/(exact*polygon_factor)-1.0).abs()<1e-12);
        assert_eq!(boundary.state_modes,(1..=modes.iter().map(Vec::len).sum::<usize>()).collect::<Vec<_>>());
        let head_faces=films[0].mesh.tris.len();
        for (i,row) in boundary.weights.iter().enumerate() {
            let head=usize::from(i>=modes[0].len());
            assert!(row[2*head_faces..].iter().all(|v|*v==0.0),"rigid walls never become invented sound sources");
            assert!(row[(1-head)*head_faces..(2-head)*head_faces].iter().all(|v|*v==0.0));
        }
    }
}
