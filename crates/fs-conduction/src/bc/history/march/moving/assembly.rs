//! Compose original moving-source projections without changing their power law.

use super::motion::ScanMotionLimit;
use super::{
    ConductionError, ConductionMesh, MovingHeatSample, MovingHeatSource, ScalarField, invalid, poll,
};
use fs_exec::Cx;

pub(in crate::bc::history::march) fn admit(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    sources: &[MovingHeatSource],
    max_sources: usize,
) -> Result<(), ConductionError> {
    poll(cx, 0)?;
    if sources.is_empty() || sources.len() > max_sources {
        return Err(invalid(
            "moving heating requires a nonempty source set within its explicit cap",
        ));
    }
    sources
        .len()
        .checked_mul(mesh.vertex_count())
        .and_then(|v| v.checked_add(sources.len().checked_mul(mesh.element_count())?))
        .ok_or_else(|| invalid("combined moving-source work count overflows"))?;
    Ok(())
}

/// The observer sees each independently checked projection before its field is
/// released. Owners may retain deposition metadata without retaining N copies
/// of the mesh field; a native transport owner needs only the summed density.
pub(in crate::bc::history::march) fn project(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    background: &ScalarField,
    sources: &[MovingHeatSource],
    time_s: f64,
    before: bool,
    mut observe: impl FnMut(&MovingHeatSample),
) -> Result<ScalarField, ConductionError> {
    poll(cx, 0)?;
    let mut field = Vec::new();
    field
        .try_reserve_exact(mesh.vertex_count())
        .map_err(|_| invalid("moving heat field allocation refused"))?;
    for i in 0..mesh.vertex_count() {
        if i % 256 == 0 {
            poll(cx, i)?;
        }
        field.push(background.at(i));
    }
    for (index, source) in sources.iter().enumerate() {
        poll(cx, index)?;
        let sample = if before {
            source.before(cx, mesh, time_s)
        } else {
            source.at(cx, mesh, time_s)
        }?;
        for (i, value) in field.iter_mut().enumerate() {
            if i % 256 == 0 {
                poll(cx, i)?;
            }
            *value += sample.source.at(i);
            crate::require_finite("combined moving heat source", *value)?;
        }
        observe(&sample);
    }
    poll(cx, field.len())?;
    Ok(ScalarField::Nodal(field))
}

pub(in crate::bc::history::march) fn admit_motion(
    cx: &Cx<'_>,
    sources: &[MovingHeatSource],
    policy: ScanMotionLimit,
) -> Result<(), ConductionError> {
    poll(cx, 0)?;
    policy.validate()?;
    for source in sources {
        source.admit_motion_limit(cx, policy)?;
    }
    poll(cx, sources.len())
}

pub(in crate::bc::history::march) fn motion_endpoint(
    cx: &Cx<'_>,
    sources: &[MovingHeatSource],
    policy: Option<ScanMotionLimit>,
    start_s: f64,
    mut end_s: f64,
) -> Result<f64, ConductionError> {
    if let Some(policy) = policy {
        for source in sources {
            end_s = source.motion_limited_endpoint(cx, start_s, end_s, policy)?;
        }
    }
    Ok(end_s)
}
