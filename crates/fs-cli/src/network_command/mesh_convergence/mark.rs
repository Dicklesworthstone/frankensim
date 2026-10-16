//! Goal-weighted gradient-recovery MARKING, not an error certificate.
//! The primal and total coupled-adjoint P1 gradients are recovered on each
//! vertex/material patch. A cell score is volume times the product of their
//! recovery defects. Global normalizations only scale all scores together.
//! Neither these scores nor their sum are used as a temperature error bound.
use super::*;

#[derive(Debug)]
pub(super) struct Marking {
    pub(super) cells:Vec<usize>,
    pub(super) total:f64,
    pub(super) captured_fraction:f64,
}
fn finite(value:f64)->Result<f64> {
    if value.is_finite() { Ok(value) } else { Err(producer("nonfinite adaptive recovery arithmetic")) }
}

pub(super) fn evaluate(cx:&Cx<'_>,request:&Request,primal:&[f64],dual:&[f64],fraction:f64)->Result<Marking> {
    let mesh=&request.mesh; let tets=&mesh.complex().tets;
    if primal.len()!=mesh.vertex_count() || dual.len()!=mesh.vertex_count()
        || primal.iter().chain(dual).any(|x|!x.is_finite()) {
        return Err(bad("adaptive marking requires complete finite primal/adjoint fields"));
    }
    let mut gradients=Vec::with_capacity(tets.len());
    let mut scales=[0.0_f64;2]; let mut max_volume=0.0_f64;
    for (index,tet) in tets.iter().enumerate() {
        if index%256==0 { poll(cx)?; }
        let mut pair=[[0.0;3];2];
        for (which,field) in [primal,dual].into_iter().enumerate() {
            // Subtract the local constant before differentiating. A uniform
            // temperature/inlet offset must not alter the marked cells.
            for local in 1..4 {
                let difference=finite(field[tet[local] as usize]-field[tet[0] as usize])?;
                for axis in 0..3 {
                    pair[which][axis]=finite(pair[which][axis]+difference*mesh.geometry().grads[index][local][axis])?;
                }
            }
            for value in pair[which] { scales[which]=scales[which].max(value.abs()); }
        }
        max_volume=max_volume.max(mesh.element_volume(index));
        gradients.push(pair);
    }
    if !(max_volume.is_finite() && max_volume>0.0) { return Err(producer("invalid adaptive cell volumes")); }
    for pair in &mut gradients {
        for which in 0..2 {
            if scales[which]>0.0 { for value in &mut pair[which] { *value/=scales[which]; } }
        }
    }
    let label=|i:usize|request.solid_data.element_materials.as_ref()
        .map_or(0,|materials|materials.of_element()[i].0);
    // Distinct material labels keep real conductivity-interface gradient jumps
    // out of a single recovery patch. Contact traces already have distinct IDs.
    let mut patches:BTreeMap<(u32,u32),([[f64;3];2],f64)>=BTreeMap::new();
    for (index,tet) in tets.iter().enumerate() {
        if index%256==0 { poll(cx)?; }
        let volume=mesh.element_volume(index)/max_volume;
        if volume<=0.0 { return Err(producer("adaptive cell-volume scaling underflow")); }
        for &vertex in tet {
            let (sums,weight)=patches.entry((vertex,label(index))).or_insert(([[0.0;3];2],0.0));
            *weight=finite(*weight+volume)?;
            for which in 0..2 { for axis in 0..3 {
                sums[which][axis]=finite(sums[which][axis]+volume*gradients[index][which][axis])?;
            } }
        }
    }
    let mut scores=Vec::with_capacity(tets.len());
    for (index,tet) in tets.iter().enumerate() {
        if index%256==0 { poll(cx)?; }
        let mut squared=[0.0;2];
        for &vertex in tet {
            let (sums,weight)=&patches[&(vertex,label(index))];
            for which in 0..2 { for axis in 0..3 {
                let difference=gradients[index][which][axis]-sums[which][axis]/weight;
                squared[which]=finite(squared[which]+difference*difference)?;
            } }
        }
        scores.push(finite((mesh.element_volume(index)/max_volume)
            *fs_math::det::sqrt(squared[0]/4.0)*fs_math::det::sqrt(squared[1]/4.0))?);
    }
    select(&scores,fraction)
}

fn select(scores:&[f64],fraction:f64)->Result<Marking> {
    if scores.is_empty() || !fraction.is_finite() || fraction<=0.0 || fraction>1.0
        || scores.iter().any(|x|!x.is_finite() || *x<0.0) {
        return Err(bad("adaptive marking requires nonnegative finite scores and a fraction in (0,1]"));
    }
    let total=scores.iter().try_fold(0.0,|sum,&value|finite(sum+value))?;
    if total==0.0 { return Ok(Marking{cells:Vec::new(),total,captured_fraction:0.0}); }
    let mut order:Vec<_>=(0..scores.len()).collect();
    order.sort_by(|&a,&b|scores[b].total_cmp(&scores[a]).then(a.cmp(&b)));
    let mut cells=Vec::new(); let mut captured=0.0;
    for index in order {
        if scores[index]==0.0 { break; }
        captured=finite(captured+scores[index])?; cells.push(index);
        if captured>=fraction*total { break; }
    }
    cells.sort_unstable();
    Ok(Marking{cells,total,captured_fraction:(captured/total).min(1.0)})
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bulk_marking_retains_the_requested_score_fraction_and_deterministic_ties() {
        let result=select(&[1.0,4.0,4.0,1.0],0.6).unwrap();
        assert_eq!(result.cells,vec![1,2]);
        assert_eq!(result.total,10.0); assert_eq!(result.captured_fraction,0.8);
        assert_eq!(select(&[2.0,2.0,2.0],0.3).unwrap().cells,vec![0]);
        assert_eq!(select(&[1.0,0.0,2.0],1.0).unwrap().cells,vec![0,2]);
    }
    #[test]
    fn no_recovery_signal_requests_a_global_probe_not_success() {
        assert!(select(&[0.0;4],0.5).unwrap().cells.is_empty());
        for scores in [vec![f64::NAN],vec![-1.0],vec![f64::MAX;2],vec![]] {
            assert!(select(&scores,0.5).is_err());
        }
        for fraction in [0.0,-1.0,1.01,f64::NAN] { assert!(select(&[1.0],fraction).is_err()); }
    }
}
