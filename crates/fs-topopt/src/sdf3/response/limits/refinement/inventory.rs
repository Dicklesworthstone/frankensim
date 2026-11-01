//! Direct inventory changes on the existing source-bound parent map.
use super::*;
use crate::sdf3::Sdf3Elasticity;

pub(super) fn normalized(value: f64, scale: f64) -> Result<f64, GoalRefinementError3> {
    let result = value/scale;
    if !value.is_finite() || value < 0.0 || !result.is_finite()
        || (value != 0.0 && result == 0.0) {
        return Err(invalid("unrepresentable normalized refinement inventory"));
    }
    Ok(result)
}
fn scale(options: ResponseLimitRefinementOptions3) -> Result<f64, GoalRefinementError3> {
    let scale = match options.objective {
        MinimumMaterialObjective3::VolumeFraction => options.volume_scale,
        MinimumMaterialObjective3::MassKg { scale_kg } => scale_kg,
    };
    if [options.volume_scale, scale].iter()
        .any(|v| !v.is_finite() || *v <= 0.0 || !(1.0/v).is_finite()) {
        return Err(invalid("refinement inventory requires positive invertible objective scales"));
    }
    Ok(scale)
}
fn weights(op: &impl Sdf3Elasticity) -> Result<&[f64], GoalRefinementError3> {
    crate::sdf3::inventory::mass_weights(op).map_err(Into::into)
}
pub(super) fn admit(coarse: &impl Sdf3Elasticity, fine: &impl Sdf3Elasticity,
    accepted: &ResponseLimitEvaluation3, options: ResponseLimitRefinementOptions3)
    -> Result<(), GoalRefinementError3> {
    scale(options)?;
    let mass = matches!(options.objective, MinimumMaterialObjective3::MassKg { .. });
    if mass != accepted.mass_kg.is_some() {
        return Err(invalid("refinement objective disagrees with the retained material inventory"));
    }
    if mass { weights(coarse)?; weights(fine)?; }
    Ok(())
}

pub(super) fn measure(op: &impl Sdf3Elasticity, projected: &[f64], volume: f64,
    options: ResponseLimitRefinementOptions3, control: &mut SolveControl<'_>)
    -> Result<(Option<f64>, f64), GoalRefinementError3> {
    let scale = scale(options)?;
    if matches!(options.objective, MinimumMaterialObjective3::VolumeFraction) {
        return Ok((None, normalized(volume, scale)?));
    }
    let mass = crate::sdf3::inventory::weighted_occupancy(weights(op)?, projected, control)?;
    Ok((Some(mass), normalized(mass, scale)?))
}

pub(super) fn cell_corrections(coarse: &impl AdaptiveSdf3Elasticity, fine: &impl Sdf3Elasticity,
    projected: &[f64], parents: &[usize], options: ResponseLimitRefinementOptions3,
    control: &mut SolveControl<'_>) -> Result<BTreeMap<Octant3, f64>, GoalRefinementError3> {
    let mut result = BTreeMap::new();
    if matches!(options.objective, MinimumMaterialObjective3::VolumeFraction) { return Ok(result); }
    let coarse_weights = weights(coarse)?; let fine_weights = weights(fine)?;
    if parents.len() != fine_weights.len() || projected.len() != coarse_weights.len() {
        return Err(invalid("refinement inventory parent-map length mismatch"));
    }
    let mut children = vec![0.0; coarse_weights.len()];
    for (i, (&weight, &parent)) in fine_weights.iter().zip(parents).enumerate() {
        if i%1024 == 0 { control.checkpoint("response-limit-goal-mass-cells")?; }
        let total = children.get_mut(parent).ok_or(invalid("invalid refinement inventory parent"))?;
        *total = number(*total+weight)?;
    }
    for (i, &leaf) in coarse.adaptive().leaves().iter().enumerate() {
        if i%1024 == 0 { control.checkpoint("response-limit-goal-mass-cells")?; }
        result.insert(leaf, number(number(children[i]-coarse_weights[i])?*projected[i])?);
    }
    Ok(result)
}
