//! Regional material-law controls, using the accepted complete native dual.
use std::collections::BTreeMap;
use fs_conduction::{ConductionProblem, adjoint::RobinResponse};
use super::{Cx, ProjectSpec, RungSolved, SolveRefusal, bad, finite, json_string,
    lower, poll, row};

/// A shared scale on a region's effective K(T) is the sum of its element
/// scale partials. Never divide by volume, replace anisotropy with a scalar,
/// merge regions sharing a card, or discard fixed-node temperature values.
pub(super) fn rows(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved,
    ids: &BTreeMap<String, u32>, lambda: &[f64], max_elements: usize,
) -> Result<Vec<String>, SolveRefusal> {
    poll(cx)?;
    let bindings = spec.materials.as_deref().unwrap_or(&[]);
    if bindings.is_empty() { return Ok(Vec::new()); }
    let mut sums = BTreeMap::<u32,(f64,usize)>::new();
    for binding in bindings {
        poll(cx)?;
        let id = *ids.get(&binding.region).ok_or_else(|| bad("material control has no native region label"))?;
        if sums.insert(id,(0.0,0)).is_some() {
            return Err(bad("material controls require one binding per native region"));
        }
    }
    if solved.labels.len() != solved.mesh.element_count() {
        return Err(bad("material controls require the exact retained element labels"));
    }
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("material control lost the accepted material operator"))?;
    let problem = ConductionProblem {mesh:&solved.mesh,boundary:&data.boundary,
        material:&data.fallback,element_materials:Some(&data.materials),source:&data.source};
    let gradient = RobinResponse::conductivity_scale_pullback(cx,problem,
        &solved.solution.temperature,lambda,max_elements).map_err(lower)?;
    for (e,(&id,value)) in solved.labels.iter().zip(gradient).enumerate() {
        if e % 512 == 0 { poll(cx)?; }
        let sum = sums.get_mut(&id).ok_or_else(|| bad("a retained element has no declared material control owner"))?;
        sum.0 = finite(sum.0+value)?;
        sum.1 += 1;
    }
    let mut result = Vec::with_capacity(bindings.len());
    for (ordinal,binding) in bindings.iter().enumerate() {
        poll(cx)?;
        let (value,count) = sums[&ids[&binding.region]];
        if count == 0 { return Err(bad("declared material control has no retained volume elements")); }
        let mut entry = row("conductivity-multiplier",&binding.region,ordinal,"1",value)?;
        // Extend the common parameter row with the exact normalized coordinate
        // and selected material identity. No card or physical input is changed.
        entry.pop();
        entry.push_str(&format!(",\"reference_value\":1,\"material_card\":{},\"material_state\":{},\"parameterization\":\"K(T) -> s K(T), evaluated at s=1\"}}",
            json_string(&binding.card),json_string(&binding.state)));
        result.push(entry);
    }
    poll(cx)?;
    Ok(result)
}
