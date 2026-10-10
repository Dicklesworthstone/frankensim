//! One source-geometry preparation for a complete harmonic frequency grid.
//! Cache exact static CBIE triangle integrals; dynamic quadrature, formulation
//! selection and the shared modal LU remain independent at every frequency.
use super::{Boundary,MAX_PANELS,MAX_FREQUENCY_WORKERS,MAX_FREQUENCY_DENSE_BYTES};
use fs_bem::{helmholtz::{Formulation,HelmholtzError,Medium,PreparedCbieGeometry,RadiationSolution},
    radiation_policy::GeometryPolicy};
use fs_math::c64::C64;

pub struct SourceSweep<'a> {
    policy:GeometryPolicy<'a>,
    geometry:Option<PreparedCbieGeometry<'a>>,
    workers:usize,
}

// Reserve the one shared static table BEFORE choosing the number of dense
// frequency operators. At 2048 panels, cached sweeps admit three workers,
// while an uncached sweep can still use the original four-worker image.
fn storage_plan(panels:usize,requested:usize,cache:bool)->Result<(usize,usize),String> {
    if !(1..=MAX_PANELS).contains(&panels) || !(1..=MAX_FREQUENCY_WORKERS).contains(&requested) {
        return Err("invalid bounded exterior source preparation".into());
    }
    let pairs=panels*panels;
    let cache_bytes=if cache {pairs*std::mem::size_of::<(f64,f64)>()} else {0};
    let dense_bytes=pairs*3*std::mem::size_of::<C64>();
    let workers=requested.min((MAX_FREQUENCY_DENSE_BYTES-cache_bytes)/dense_bytes);
    if workers==0 {return Err("exterior source exceeds aggregate dense and static-geometry memory budget".into());}
    Ok((workers,cache_bytes))
}

impl<'a> SourceSweep<'a> {
    pub fn new(boundary:&'a Boundary,omega:&[f64],medium:Medium,requested_workers:usize)
        ->Result<Self,String> {
        if omega.is_empty() || omega.len()>257
            || omega.iter().enumerate().any(|(i,w)|!w.is_finite() || *w<=0.
                || (i>0 && *w<=omega[i-1]))
            || !medium.density.is_finite() || medium.density<=0.
            || !medium.sound_speed.is_finite() || medium.sound_speed<=0. {
            return Err("invalid exterior source frequency grid or medium".into());
        }
        storage_plan(boundary.surface.areas().len(),requested_workers,false)?;
        let policy=GeometryPolicy::new(&boundary.surface).map_err(|e|e.to_string())?;
        let mut plain_count=0;let mut highest_plain=0.;
        for &w in omega {
            let k=w/medium.sound_speed;
            if policy.formulation(k).map_err(|e|e.to_string())?==Formulation::PlainCbie {
                plain_count+=1;highest_plain=k;
            }
        }
        let (workers,cache_bytes)=storage_plan(boundary.surface.areas().len(),requested_workers,plain_count>1)?;
        let geometry=if cache_bytes>0 {
            Some(PreparedCbieGeometry::new(&boundary.surface,highest_plain,cache_bytes)
                .map_err(|e|format!("exterior static geometry through {} Hz: {e}",
                    highest_plain*medium.sound_speed/std::f64::consts::TAU))?)
        } else {None};
        Ok(Self {policy,geometry,workers})
    }
    pub fn workers(&self)->usize {self.workers}
    pub fn formulation(&self,k:f64)->Result<Formulation,HelmholtzError> {self.policy.formulation(k)}
    pub fn solve_batch(&self,k:f64,medium:Medium,fields:&[&[C64]])
        ->Result<Vec<RadiationSolution>,HelmholtzError> {
        if self.policy.formulation(k)?==Formulation::PlainCbie {
            if let Some(geometry)=&self.geometry {return geometry.solve_batch(k,medium,fields);}
        }
        self.policy.solve_batch(k,medium,fields)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{tests as fixture,Specification};

    fn boundary()->(Boundary,Specification) {
        let spec=Specification::read(&fixture::specification()).unwrap();
        let boundary=Boundary::from_obj(&fixture::box_obj("skin",[0.,0.,-0.01],[0.1,0.1,0.02]),
            &spec,&fixture::motion()).unwrap();
        (boundary,spec)
    }

    #[test]
    fn harmonic_source_cache_matches_original_fields_and_power_at_every_frequency() {
        let (boundary,spec)=boundary();let omega=spec.omega();
        let prepared=SourceSweep::new(&boundary,&omega,spec.medium,2).unwrap();
        assert_eq!(prepared.geometry.as_ref().unwrap().cache_bytes(),16*12*12);
        let original=GeometryPolicy::new(&boundary.surface).unwrap();
        let a:Vec<_>=boundary.surface.normals().iter().map(|n|C64::new(n[2],0.)).collect();
        let b:Vec<_>=boundary.surface.normals().iter().map(|n|C64::new(0.,n[0])).collect();
        for w in omega {
            let k=w/spec.medium.sound_speed;
            let cached=prepared.solve_batch(k,spec.medium,&[&a,&b]).unwrap();
            let direct=original.solve_batch(k,spec.medium,&[&a,&b]).unwrap();
            for (cached,direct) in cached.iter().zip(direct) {
                assert_eq!(cached.pressure,direct.pressure);
                assert_eq!(cached.radiated_power_roundoff_interval,direct.radiated_power_roundoff_interval);
                assert_eq!(cached.condition_lower_bound,direct.condition_lower_bound);
                assert_eq!(cached.panels_per_wavelength,direct.panels_per_wavelength);
            }
        }
    }

    #[test]
    fn harmonic_source_cache_preserves_frequency_formulation_and_one_shot_admission() {
        let (boundary,spec)=boundary();
        let limit=GeometryPolicy::new(&boundary.surface).unwrap().plain_cbie_limit();
        let omega=[0.01*limit*spec.medium.sound_speed,0.02*limit*spec.medium.sound_speed,
            1.01*limit*spec.medium.sound_speed];
        let mixed=SourceSweep::new(&boundary,&omega,spec.medium,1).unwrap();
        assert!(mixed.geometry.is_some());
        assert_eq!(mixed.formulation(omega[0]/spec.medium.sound_speed).unwrap(),Formulation::PlainCbie);
        assert_eq!(mixed.formulation(omega[2]/spec.medium.sound_speed).unwrap(),Formulation::BurtonMiller);
        let single=SourceSweep::new(&boundary,&omega[..1],spec.medium,1).unwrap();
        assert!(single.geometry.is_none());
        let high=SourceSweep::new(&boundary,&omega[2..],spec.medium,1).unwrap();
        assert!(high.geometry.is_none());
        for invalid in [vec![],vec![0.],vec![f64::NAN],vec![2.,1.],vec![1.;258]] {
            assert!(SourceSweep::new(&boundary,&invalid,spec.medium,1).is_err());
        }
    }

    #[test]
    fn harmonic_source_cache_shares_the_existing_dense_memory_budget() {
        assert_eq!(storage_plan(2048,4,false).unwrap(),(4,0));
        assert_eq!(storage_plan(2048,4,true).unwrap(),(3,64*1024*1024));
        for panels in [12,1024,1536,2048] {for requested in 1..=4 {
            let (workers,cache)=storage_plan(panels,requested,true).unwrap();
            assert!(workers>0 && workers<=requested);
            assert!(workers*panels*panels*3*std::mem::size_of::<C64>()+cache<=MAX_FREQUENCY_DENSE_BYTES);
        }}
        for (panels,workers) in [(0,1),(2049,1),(12,0),(12,5)] {
            assert!(storage_plan(panels,workers,true).is_err());
        }
    }
}
