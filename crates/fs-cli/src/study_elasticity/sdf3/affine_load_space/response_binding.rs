//! Bind authored physical observations to the actual fixed cut-cell grid.
//! Displacement means use reference cut volume, including supported material;
//! support forces/moments remain live native material-dependent observations.
use super::*;
use fs_cutfem::elastic3::{
    ElasticityError3, adaptive::AdaptiveElasticity3, surface::ReferenceLoad3,
};
use fs_topopt::sdf3::response::limits::refinement::{
    ReferenceBoundedResponse3, ReferenceResponseLimit3,
};
use fs_topopt::sdf3::response::limits::{BoundedResponse3, ResponseLimit3};
use fs_topopt::sdf3::robust::envelope::affine_load_space::response::reference::ReferenceAffineResponseLimits3;
use spec::responses::{Family, Observation};

pub(super) struct Prepared<'a> {
    family: &'a Family,
    observations: Vec<Option<Vec<f64>>>,
}

#[derive(Clone, Copy)]
struct Mean {
    axis: usize,
    x: [f64; 2],
    volume: f64,
}
impl Mean {
    fn density(self, p: [f64; 3]) -> [f64; 3] {
        let mut weight = [0.0; 3];
        if self.x[0] <= p[0] && p[0] <= self.x[1] {
            weight[self.axis] = 1.0 / self.volume;
        }
        weight
    }
}

pub(super) struct ReferencePrepared<'a> {
    family: &'a Family,
    means: Vec<Option<Mean>>,
}
impl ReferencePrepared<'_> {
    pub(super) fn mean_volumes(&self) -> Vec<Option<f64>> {
        self.means
            .iter()
            .map(|mean| mean.map(|m| m.volume))
            .collect()
    }
}

pub(super) fn assemble<'a>(
    spec: &'a Spec,
    operator: &AdaptiveElasticity3,
    mut checkpoint: impl FnMut() -> ControlFlow<()>,
) -> std::result::Result<Option<Prepared<'a>>, ElasticityError3> {
    let Some(reference) = prepare_reference(spec, operator, &mut checkpoint)? else {
        return Ok(None);
    };
    let mut observations = Vec::with_capacity(reference.means.len());
    for mean in &reference.means {
        if checkpoint().is_break() {
            return Err(ElasticityError3::Cancelled);
        }
        observations.push(match mean {
            Some(mean) => Some(
                operator
                    .reference_load(ReferenceLoad3::body(&|p| mean.density(p)), &mut checkpoint)?,
            ),
            None => None,
        });
    }
    Ok(Some(Prepared {
        family: reference.family,
        observations,
    }))
}

/// Bind each mean's geometry independently before borrowing either operator
/// for native two-grid physics. The original physical patch stays fixed.
pub(super) fn prepare_reference<'a>(
    spec: &'a Spec,
    operator: &AdaptiveElasticity3,
    mut checkpoint: impl FnMut() -> ControlFlow<()>,
) -> std::result::Result<Option<ReferencePrepared<'a>>, ElasticityError3> {
    let Some(family) = &spec.affine_responses else {
        return Ok(None);
    };
    let volumes = operator.volumes();
    let mut means = Vec::with_capacity(family.limits.len());
    for limit in &family.limits {
        if checkpoint().is_break() {
            return Err(ElasticityError3::Cancelled);
        }
        let mean = match limit.observation {
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
                Some(Mean { axis, x, volume })
            }
            Observation::SupportForce { .. } | Observation::SupportMoment { .. } => None,
        };
        means.push(mean);
    }
    if checkpoint().is_break() {
        return Err(ElasticityError3::Cancelled);
    }
    Ok(Some(ReferencePrepared { family, means }))
}

fn motion(observation: Observation, p: [f64; 3]) -> [f64; 3] {
    match observation {
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
            move |p: [f64; 3], _normal: [f64; 3]| motion(observation, p)
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

pub(super) fn with_reference_limits<T>(
    prepared: Option<&ReferencePrepared<'_>>,
    apply: impl FnOnce(Option<ReferenceAffineResponseLimits3<'_>>) -> T,
) -> T {
    let Some(prepared) = prepared else {
        return apply(None);
    };
    let bodies: Vec<_> = prepared
        .means
        .iter()
        .map(|&mean| move |p| mean.map_or([0.0; 3], |mean| mean.density(p)))
        .collect();
    let modes: Vec<_> = prepared
        .family
        .limits
        .iter()
        .map(|limit| {
            let observation = limit.observation;
            move |p: [f64; 3], _normal: [f64; 3]| motion(observation, p)
        })
        .collect();
    let limits: Vec<_> = prepared
        .family
        .limits
        .iter()
        .enumerate()
        .map(|(i, limit)| ReferenceResponseLimit3 {
            response: if prepared.means[i].is_some() {
                ReferenceBoundedResponse3::Displacement(ReferenceLoad3::body(&bodies[i]))
            } else {
                ReferenceBoundedResponse3::Reaction(&modes[i])
            },
            sense: limit.sense,
            bound: limit.bound,
            scale: limit.scale,
        })
        .collect();
    apply(Some(ReferenceAffineResponseLimits3 {
        limits: &limits,
        options: prepared.family.options,
    }))
}
