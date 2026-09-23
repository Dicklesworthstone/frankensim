//! Both physical cymbal skins in ONE stationary exterior solve.
use super::*;
impl Boundary {
    pub(crate) fn shell_pair(upper:&ShellRadiationSurface,upper_start:usize,
        lower:&ShellRadiationSurface,lower_start:usize,separation:f64)->Result<Self,Error> {
        let ua=upper.normal_velocity_weights().len();let lb=lower.normal_velocity_weights().len();
        let up=upper.triangles().len();let lp=lower.triangles().len();
        let count=ua.checked_add(lb).ok_or("paired acoustic input overflow")?;
        let panels=up.checked_add(lp).ok_or("paired acoustic panel overflow")?;
        let ue=upper_start.checked_add(ua).ok_or("paired upper address overflow")?;
        let le=lower_start.checked_add(lb).ok_or("paired lower address overflow")?;
        if ua==0||lb==0||count>MAX_INPUTS||panels>MAX_PANELS||!separation.is_finite()||separation<=0.
            ||upper_start<le && lower_start<ue {
            return Err("paired boundary exceeds original budgets or overlaps source addresses".into());
        }
        let mut triangles=Vec::with_capacity(panels);
        triangles.extend(upper.triangles().iter().map(|t|t.map(|p|[p[0],p[1],p[2]+0.5*separation])));
        // Proper rigid rotation, determinant +1: normal and velocity rotate
        // together, so their scalar product has NO additional sign change.
        triangles.extend(lower.triangles().iter().map(|t|t.map(|p|[p[0],-p[1],-p[2]-0.5*separation])));
        if triangles.iter().flatten().flatten().any(|x|!x.is_finite()) {return Err("paired placement overflows".into());}
        let upper_min=triangles[..up].iter().flatten().map(|p|p[2]).fold(f64::INFINITY,f64::min);
        let lower_max=triangles[up..].iter().flatten().map(|p|p[2]).fold(f64::NEG_INFINITY,f64::max);
        if upper_min<=lower_max {return Err("paired acoustic skins require a separating reference plane".into());}
        let mut weights=Vec::with_capacity(count);
        for row in upper.normal_velocity_weights(){let mut r=row.clone();r.resize(panels,0.);weights.push(r);}
        for row in lower.normal_velocity_weights(){let mut r=vec![0.;up];r.extend(row);weights.push(r);}
        let state_modes=(upper_start..ue).chain(lower_start..le).collect();
        Ok(Self{triangles,weights,state_modes})
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn paired_scene_retains_both_skin_orientations_and_source_flux_without_cross_wiring() {
        use crate::{shell_prepare,specimen};
        use fs_plate::shell::reduction::radiation::RadiationSurfaceBudget;
        let mut s=specimen::Specimen::reference();s.azimuths=8;s.band_hz=[10.,11.];
        let (mesh,r)=shell_prepare::prepare(&s,2e-6).unwrap();
        let skin=r.radiation_surface(&mesh.nodal_thickness_m,
            RadiationSurfaceBudget{max_panels:2048,max_panel_modes:65536}).unwrap();
        let n=r.mode_count();let p=skin.triangles().len();
        let b=Boundary::shell_pair(&skin,1,&skin,1+n,0.002).unwrap();
        assert_eq!(b.triangles.len(),2*p);assert_eq!(b.weights.len(),2*n);
        for k in 0..n {
            assert_eq!(&b.weights[k][..p],&skin.normal_velocity_weights()[k]);
            assert!(b.weights[k][p..].iter().all(|x|*x==0.));
            assert!(b.weights[n+k][..p].iter().all(|x|*x==0.));
            assert_eq!(&b.weights[n+k][p..],&skin.normal_velocity_weights()[k]);
        }
        let actual=SpherePanels::from_triangles(b.triangles.clone()).unwrap();
        for (i,t) in skin.triangles().iter().enumerate(){
            for a in 0..3 {assert_eq!(b.triangles[p+i][a],[t[a][0],-t[a][1],-t[a][2]-0.001]);}
            // A proper rotation preserves material-side normals and flux.
            assert!((actual.areas()[i]-actual.areas()[p+i]).abs()<1e-15);
            assert!((actual.normals()[i][2]+actual.normals()[p+i][2]).abs()<1e-12);
        }
        for row in &b.weights {
            assert!(row.iter().zip(actual.areas()).map(|(v,a)|v*a).sum::<f64>().abs()<1e-12);
        }
        assert!(Boundary::shell_pair(&skin,1,&skin,1,0.002).is_err());
        assert!(Boundary::shell_pair(&skin,1,&skin,1+n,0.00001).is_err());
        assert!(Boundary::shell_pair(&skin,usize::MAX,&skin,1,0.002).is_err());
    }
}
