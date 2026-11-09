//! Bind authored physical observations to the actual fixed cut-cell grid.
//! Displacement means use reference cut volume, including supported material;
//! support forces/moments remain live native material-dependent observations.
use super::*;
use fs_cutfem::elastic3::{
    ElasticityError3, adaptive::AdaptiveElasticity3, surface::ReferenceLoad3,
};
use fs_topopt::sdf3::response::limits::{BoundedResponse3, ResponseLimit3};
use spec::responses::{Family, Observation};

pub(super) struct Prepared<'a> {
    family: &'a Family,
    observations: Vec<Option<Vec<f64>>>,
}

pub(super) fn assemble<'a>(
    spec: &'a Spec,
    operator: &AdaptiveElasticity3,
    mut checkpoint: impl FnMut() -> ControlFlow<()>,
) -> std::result::Result<Option<Prepared<'a>>, ElasticityError3> {
    let Some(family) = &spec.affine_responses else {
        return Ok(None);
    };
    let volumes = operator.volumes();
    let mut observations = Vec::with_capacity(family.limits.len());
    for limit in &family.limits {
        if checkpoint().is_break() {
            return Err(ElasticityError3::Cancelled);
        }
        let observation = match limit.observation {
            Observation::MeanDisplacement { axis, x_fraction } => {
                let side = f64::from(1_u32 << spec.level);
                let planes = x_fraction.map(|fraction| (side * fraction) as u32);
                let mut volume = 0.0;
                for (leaf, &measure) in operator.leaves().iter().zip(&volumes) {
                    if checkpoint().is_break() {
                        return Err(ElasticityError3::Cancelled);
                    }
                    let level = u32::from(leaf.level());
                    if level < spec.level {
                        return Err(ElasticityError3::Invalid(
                            "response patch cannot coarsen the initial lattice",
                        ));
                    }
                    let ancestor = leaf.index()[0] >> (level - spec.level);
                    if planes[0] <= ancestor && ancestor < planes[1] {
                        volume += measure;
                    }
                }
                if !volume.is_finite() || volume <= 0.0 || !(1.0 / volume).is_finite() {
                    return Err(ElasticityError3::Invalid(
                        "mean-displacement patch has no finite positive reference cut volume",
                    ));
                }
                let lower = spec.bounds.0[0];
                let upper = spec.bounds.1[0];
                let extent = upper - lower;
                let x = x_fraction.map(|fraction| {
                    if fraction == 1.0 {
                        upper
                    } else {
                        lower + extent * fraction
                    }
                });
                let density = |p: [f64; 3]| {
                    let mut weight = [0.0; 3];
                    if x[0] <= p[0] && p[0] <= x[1] {
                        weight[axis] = 1.0 / volume;
                    }
                    weight
                };
                Some(operator.reference_load(ReferenceLoad3::body(&density), &mut checkpoint)?)
            }
            Observation::SupportForce { .. } | Observation::SupportMoment { .. } => None,
        };
        observations.push(observation);
    }
    if checkpoint().is_break() {
        return Err(ElasticityError3::Cancelled);
    }
    Ok(Some(Prepared {
        family,
        observations,
    }))
}

/// The vectors and pure reaction modes outlive every evaluation, audit and
/// restored endpoint in `apply`; no density-dependent reaction is frozen as q.
pub(super) fn with_limits<T>(
    prepared: Option<&Prepared<'_>>,
    apply: impl FnOnce(Option<AffineResponseLimits3<'_>>) -> T,
) -> T {
    let Some(prepared) = prepared else {
        return apply(None);
    };
    let modes: Vec<_> = prepared
        .family
        .limits
        .iter()
        .map(|limit| {
            let observation = limit.observation;
            move |p: [f64; 3], _normal: [f64; 3]| match observation {
                Observation::SupportForce { axis } => {
                    let mut direction = [0.0; 3];
                    direction[axis] = 1.0;
                    direction
                }
                Observation::SupportMoment { axis, about_m } => {
                    let mut direction = [0.0; 3];
                    direction[axis] = 1.0;
                    let r: [f64; 3] = std::array::from_fn(|i| p[i] - about_m[i]);
                    [
                        direction[1] * r[2] - direction[2] * r[1],
                        direction[2] * r[0] - direction[0] * r[2],
                        direction[0] * r[1] - direction[1] * r[0],
                    ]
                }
                Observation::MeanDisplacement { .. } => [0.0; 3],
            }
        })
        .collect();
    let limits: Vec<_> = prepared
        .family
        .limits
        .iter()
        .zip(&prepared.observations)
        .zip(&modes)
        .map(|((limit, observation), mode)| ResponseLimit3 {
            response: match observation {
                Some(q) => BoundedResponse3::Displacement(q),
                None => BoundedResponse3::Reaction(mode),
            },
            sense: limit.sense,
            bound: limit.bound,
            scale: limit.scale,
        })
        .collect();
    apply(Some(AffineResponseLimits3 {
        limits: &limits,
        options: prepared.family.options,
    }))
}
