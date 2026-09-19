//! Geometry-backed operator interface shared by uniform and adaptive studies.
use std::collections::BTreeMap;
use fs_cutfem::elastic3::{CutElasticity3,ElasticityError3};
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_solver::op::LinearOp;
mod sealed {pub trait Sealed {}}
impl sealed::Sealed for CutElasticity3 {}
impl sealed::Sealed for AdaptiveElasticity3 {}
/// The two geometry-owned 3-D elasticity backends admitted by the density study.
/// Every index below refers to the SAME retained active-cell/unknown ordering.
pub trait Sdf3Elasticity: LinearOp + sealed::Sealed {
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
    fn cells(&self)->usize {AdaptiveElasticity3::cells(self)}
    fn volumes(&self)->Vec<f64> {AdaptiveElasticity3::volumes(self)}
    fn scales(&self)->&[f64] {AdaptiveElasticity3::scales(self)}
    fn set_scales(&mut self,v:&[f64])->Result<(),ElasticityError3> {AdaptiveElasticity3::set_scales(self,v)}
    fn fixed(&self)->&[bool] {AdaptiveElasticity3::fixed(self)}
    fn scale_quadratic_forms(&self,u:&[f64])->Result<Vec<f64>,ElasticityError3> {AdaptiveElasticity3::scale_quadratic_forms(self,u)}
    fn filter_edges(&self)->Vec<(usize,usize,f64)> {AdaptiveElasticity3::filter_edges(self).to_vec()}
}
