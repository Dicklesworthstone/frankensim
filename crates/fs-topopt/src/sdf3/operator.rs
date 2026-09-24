//! Geometry-backed operator interface shared by uniform and adaptive studies.
use std::collections::BTreeMap;
use std::ops::ControlFlow;
use fs_cutfem::elastic3::{CutElasticity3,ElasticityError3};
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::{AdaptivePrepared3, AdaptivePreconditionError3, AdaptiveSolveSpace3, AdaptiveSetupWork3};
use fs_solver::op::{LinearOp, two_level::TwoLevelError, multilevel::MultilevelError};
use fs_sparse::precond::{IdentityPrecond, Precond};
use crate::{EvaluationStop, SolveControl};
mod sealed {pub trait Sealed {}}
impl sealed::Sealed for CutElasticity3 {}
impl sealed::Sealed for AdaptiveElasticity3 {}
impl sealed::Sealed for AdaptiveSolveSpace3 {}
/// Geometry-owned 3-D elasticity backends admitted by the density study.
/// Every index below refers to the SAME retained active-cell/unknown ordering.
pub trait Sdf3Elasticity: LinearOp + sealed::Sealed {
    /// Fixed linear SPD action, borrowing exactly the current density state.
    type Prepared<'a>: Precond where Self: 'a;
    /// Prepare once per design evaluation, before any independent load solve.
    /// Existing bare operators return identity without setup or new checkpoints.
    fn prepare_elasticity(&self, control: &mut SolveControl<'_>) -> Result<Self::Prepared<'_>, EvaluationStop>;
    /// Number of independent density variables.
    fn cells(&self)->usize;
    /// Numerical active-cell measures.
    fn volumes(&self)->Vec<f64>;
    /// Current reference-material multipliers.
    fn scales(&self)->&[f64];
    /// Transactional reference-material update.
    fn set_scales(&mut self,values:&[f64])->Result<(),ElasticityError3>;
    /// One flag per independent displacement node.
    fn fixed(&self)->&[bool];
    /// Bulk-plus-ghost contractions with any hanging-node map applied.
    fn scale_quadratic_forms(&self,u:&[f64])->Result<Vec<f64>,ElasticityError3>;
    /// Each shared active face once, with normal center separation.
    fn filter_edges(&self)->Vec<(usize,usize,f64)>;
}
impl Sdf3Elasticity for CutElasticity3 {
    type Prepared<'a> = IdentityPrecond;
    fn prepare_elasticity(&self, _: &mut SolveControl<'_>) -> Result<IdentityPrecond, EvaluationStop> { Ok(IdentityPrecond) }
    fn cells(&self)->usize {CutElasticity3::cells(self)}
    fn volumes(&self)->Vec<f64> {CutElasticity3::volumes(self)}
    fn scales(&self)->&[f64] {CutElasticity3::scales(self)}
    fn set_scales(&mut self,v:&[f64])->Result<(),ElasticityError3> {CutElasticity3::set_scales(self,v)}
    fn fixed(&self)->&[bool] {CutElasticity3::fixed(self)}
    fn scale_quadratic_forms(&self,u:&[f64])->Result<Vec<f64>,ElasticityError3> {CutElasticity3::scale_quadratic_forms(self,u)}
    fn filter_edges(&self)->Vec<(usize,usize,f64)> {
        // Preserve the original Cartesian neighbor and insertion order.
        let keys=self.cell_keys();
        let ids:BTreeMap<_,_>=keys.iter().enumerate().map(|(i,&k)|(k,i)).collect();
        let centers:Vec<[f64;3]>=self.cell_nodes().iter().map(|nodes| {
            std::array::from_fn(|a|f64::midpoint(self.nodes()[nodes[0]][a],self.nodes()[nodes[7]][a]))
        }).collect();
        let mut edges=Vec::new();
        for (i,&key) in keys.iter().enumerate() {for axis in 0..3 {
            let mut next=key;next[axis]+=1;
            if let Some(&j)=ids.get(&next) {edges.push((i,j,centers[j][axis]-centers[i][axis]));}
        }}edges
    }
}
impl Sdf3Elasticity for AdaptiveElasticity3 {
    type Prepared<'a> = IdentityPrecond;
    fn prepare_elasticity(&self, _: &mut SolveControl<'_>) -> Result<IdentityPrecond, EvaluationStop> { Ok(IdentityPrecond) }
    fn cells(&self)->usize {AdaptiveElasticity3::cells(self)}
    fn volumes(&self)->Vec<f64> {AdaptiveElasticity3::volumes(self)}
    fn scales(&self)->&[f64] {AdaptiveElasticity3::scales(self)}
    fn set_scales(&mut self,v:&[f64])->Result<(),ElasticityError3> {AdaptiveElasticity3::set_scales(self,v)}
    fn fixed(&self)->&[bool] {AdaptiveElasticity3::fixed(self)}
    fn scale_quadratic_forms(&self,u:&[f64])->Result<Vec<f64>,ElasticityError3> {AdaptiveElasticity3::scale_quadratic_forms(self,u)}
    fn filter_edges(&self)->Vec<(usize,usize,f64)> {AdaptiveElasticity3::filter_edges(self).to_vec()}
}
impl Sdf3Elasticity for AdaptiveSolveSpace3 {
    type Prepared<'a> = AdaptivePrepared3<'a>;
    fn prepare_elasticity(&self, control: &mut SolveControl<'_>) -> Result<AdaptivePrepared3<'_>, EvaluationStop> {
        control.checkpoint("sdf3-preconditioner-start")?;
        let mut recorded = AdaptiveSetupWork3::default();
        let mut setup_stop = None;
        let prepared = self.prepare_with_work(|work| {
            let additional = work.operator_applications.checked_sub(recorded.operator_applications)
                .zip(work.galerkin_products.checked_sub(recorded.galerkin_products));
            let Some((applications, products)) = additional else {
                setup_stop = Some(EvaluationStop::Breakdown { stage: "preconditioner-accounting" });
                return ControlFlow::Break(());
            };
            recorded = work;
            match control.record_preconditioner_setup(applications, products) {
                Ok(()) => ControlFlow::Continue(()),
                Err(stop) => { setup_stop = Some(stop); ControlFlow::Break(()) }
            }
        });
        if let Some(stop) = setup_stop { return Err(stop); }
        let prepared = prepared.map_err(|e| match e {
            AdaptivePreconditionError3::Physics(ElasticityError3::Cancelled)
                | AdaptivePreconditionError3::Coarse(TwoLevelError::Cancelled)
                | AdaptivePreconditionError3::Hierarchy(MultilevelError::Cancelled) => EvaluationStop::Cancelled,
            AdaptivePreconditionError3::Coarse(TwoLevelError::Budget(_))
                | AdaptivePreconditionError3::Hierarchy(MultilevelError::Budget(_)) =>
                EvaluationStop::TotalBudget { stage: "sdf3-preconditioner-setup" },
            _ => EvaluationStop::Breakdown { stage: "sdf3-preconditioner" },
        })?;
        control.checkpoint("sdf3-preconditioner-ready")?;
        Ok(prepared)
    }
    fn cells(&self)->usize {self.elasticity().cells()}
    fn volumes(&self)->Vec<f64> {self.elasticity().volumes()}
    fn scales(&self)->&[f64] {self.elasticity().scales()}
    fn set_scales(&mut self,v:&[f64])->Result<(),ElasticityError3> {AdaptiveSolveSpace3::set_scales(self,v)}
    fn fixed(&self)->&[bool] {self.elasticity().fixed()}
    fn scale_quadratic_forms(&self,u:&[f64])->Result<Vec<f64>,ElasticityError3> {self.elasticity().scale_quadratic_forms(u)}
    fn filter_edges(&self)->Vec<(usize,usize,f64)> {self.elasticity().filter_edges().to_vec()}
}

/// Adaptive geometry behind either a bare or explicitly preconditioned study.
/// Goal estimation consumes exactly these same scales, fields and constraints.
pub trait AdaptiveSdf3Elasticity: Sdf3Elasticity {
    /// Read-only adaptive operator, not a reassembled approximation.
    fn adaptive(&self) -> &AdaptiveElasticity3;
    /// Consume the solve wrapper while retaining its exact geometry and scales.
    /// An enrichment estimator chooses and accounts for its own solve policy.
    fn into_adaptive(self) -> AdaptiveElasticity3;
}
impl AdaptiveSdf3Elasticity for AdaptiveElasticity3 {
    fn adaptive(&self) -> &AdaptiveElasticity3 { self }
    fn into_adaptive(self) -> AdaptiveElasticity3 { self }
}
impl AdaptiveSdf3Elasticity for AdaptiveSolveSpace3 {
    fn adaptive(&self) -> &AdaptiveElasticity3 { self.elasticity() }
    fn into_adaptive(self) -> AdaptiveElasticity3 { self.into_elasticity() }
}
