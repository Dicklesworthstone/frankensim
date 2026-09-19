//! Raw-design initialization after local refinement. No solved fields transfer.
use std::collections::BTreeMap;
use fs_cutfem::octree3::Octant3;
use crate::{EvaluationStop,SolveControl};
use super::{CutDensityStudy3,Sdf3Elasticity};

/// Copy each target leaf's closest source-ancestor RAW density. This is an
/// initialization operation, not reuse of a solved displacement, projected
/// material volume or certificate. Caller owns geometric/problem identity.
/// New filter/projection/volume and load solves must be evaluated afterward.
/// Missing ancestors and duplicates refuse rather than zero-filling new cells.
pub fn inherit_raw_densities3(source:&[Octant3],rho:&[f64],target:&[Octant3],control:&mut SolveControl<'_>)
    ->Result<Vec<f64>,EvaluationStop> {
    control.checkpoint("sdf3-density-transfer")?;
    assert_eq!(source.len(),rho.len(),"one source density per leaf required");
    assert!(rho.iter().all(|r|r.is_finite()&&(1e-3..=1.0).contains(r)),"invalid source densities");
    let mut by_key=BTreeMap::new();
    for (&leaf,&value) in source.iter().zip(rho) {
        control.checkpoint("sdf3-density-transfer")?;
        if by_key.insert((leaf.level(),leaf.index()),value).is_some() {
            return Err(EvaluationStop::Breakdown{stage:"sdf3-duplicate-source-leaf"});
        }
    }
    let mut result=Vec::with_capacity(target.len());
    for &leaf in target {
        control.checkpoint("sdf3-density-transfer")?;
        let value=(0..=leaf.level()).rev().find_map(|level| {
            let index=leaf.index().map(|i|i>>(leaf.level()-level));
            by_key.get(&(level,index)).copied()
        }).ok_or(EvaluationStop::Breakdown{stage:"sdf3-missing-source-ancestor"})?;
        result.push(value);
    }
    control.checkpoint("sdf3-density-transfer-publish")?;Ok(result)
}
impl<O:Sdf3Elasticity> CutDensityStudy3<O> {
    /// Restore a starting design's actual projected-volume feasibility under
    /// this study's OWN new filter, geometry measures and projection. At most
    /// 64 filter evaluations bisect toward the raw density floor. No operator
    /// stiffness or accepted field is mutated. A failed floor test is failure
    /// of this restoration path, not a proof of global infeasibility.
    pub fn feasible_start(&self,incoming:&[f64],cap:f64,tolerance:f64,control:&mut SolveControl<'_>)
        ->Result<Vec<f64>,EvaluationStop> {
        assert!(cap.is_finite()&&cap>0.0&&cap<=1.0&&tolerance.is_finite()&&tolerance>=0.0&&tolerance<cap,"invalid volume restoration policy");
        assert!(incoming.iter().all(|r|r.is_finite()&&(1e-3..=1.0).contains(r)),"invalid incoming raw design");
        let limit=cap+tolerance;
        if self.design(incoming,control)?.volume<=limit {return Ok(incoming.to_vec());}
        let mut feasible=vec![1e-3;incoming.len()];
        if self.design(&feasible,control)?.volume>limit {return Err(EvaluationStop::Breakdown{stage:"sdf3-volume-restoration"});}
        let (mut low,mut high)=(0.0,1.0);
        for _ in 0..64 {
            control.checkpoint("sdf3-volume-restoration")?;
            let mid=f64::midpoint(low,high);if mid<=low||mid>=high {break;}
            let trial:Vec<f64>=incoming.iter().map(|r|(1e-3+mid*(r-1e-3)).clamp(1e-3,1.0)).collect();
            if self.design(&trial,control)?.volume<=limit {low=mid;feasible=trial;}else{high=mid;}
        }
        control.checkpoint("sdf3-volume-restoration-publish")?;Ok(feasible)
    }
}
