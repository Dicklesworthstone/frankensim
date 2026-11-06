//! Conservative nesting admission; no generic intersection or repair engine.
use super::*;

pub(super) fn admit(
    outer: &SpherePanels,
    core: &SpherePanels,
    minimum_clearance_m: f64,
    options: Options,
    gate: &CancelGate,
) -> Result<usize, TransmissionError> {
    if !minimum_clearance_m.is_finite()
        || minimum_clearance_m <= 0.0
        || options.maximum_panels == 0
        || options.maximum_panels > MAX_TRANSMISSION_PANELS
    {
        return Err(bad("invalid core-shell clearance or panel budget").into());
    }
    let count = outer
        .areas()
        .len()
        .checked_add(core.areas().len())
        .ok_or_else(|| bad("core-shell panel count overflow"))?;
    if count > options.maximum_panels {
        return Err(HelmholtzError::WorkCap { panels: count }.into());
    }
    let mut scale = 0.0_f64;
    for surface in [outer, core] {
        poll(gate)?;
        surface
            .validate()
            .map_err(|_| bad("invalid core-shell surface"))?;
        let triangles = surface
            .triangles()
            .ok_or_else(|| bad("core-shell interfaces require retained triangles"))?;
        if triangles.is_empty() || near_field::components(triangles)?.len() != 1 {
            return Err(bad("each core-shell interface needs one closed outward component").into());
        }
        for vertex in triangles.iter().flatten() {
            if vertex
                .iter()
                .any(|value| !value.is_finite() || value.abs() > 1e9)
            {
                return Err(bad("core-shell coordinates exceed the finite receiver domain").into());
            }
            scale = vertex
                .iter()
                .fold(scale, |largest, value| largest.max(value.abs()));
        }
    }
    let roundoff = 128.0 * f64::EPSILON * scale;
    let clearance = minimum_clearance_m + roundoff;
    if !clearance.is_finite() {
        return Err(bad("core-shell clearance plus numerical guard is unrepresentable").into());
    }
    let outer_triangles = outer.triangles().expect("admitted triangles");
    let core_triangles = core.triangles().expect("admitted triangles");
    for (face, triangle) in outer_triangles.iter().enumerate() {
        poll(gate)?;
        let normal = outer.normals()[face];
        let origin = triangle[0];
        // Every retained outer vertex lies behind every outward face. This
        // intentionally refuses concave shells rather than infer nesting from
        // vertices or a bounding sphere in an unsupported geometry domain.
        // The roundoff tolerance is a numerical admission allowance, not an
        // exact-arithmetic convexity or global nonintersection certificate.
        for (vertex, &point) in outer_triangles.iter().flatten().enumerate() {
            if vertex % 256 == 0 {
                poll(gate)?;
            }
            let signed = dot(normal, sub(point, origin));
            if !signed.is_finite() || signed > roundoff {
                return Err(bad("outer core-shell interface must be convex").into());
            }
        }
        // Entire core triangles lie in these convex half-spaces whenever
        // their vertices do. Plane distance is a conservative lower bound on
        // distance to each finite outer face, including edges and corners.
        for (vertex, &point) in core_triangles.iter().flatten().enumerate() {
            if vertex % 256 == 0 {
                poll(gate)?;
            }
            let signed = dot(normal, sub(point, origin));
            if !signed.is_finite() || signed > -clearance {
                return Err(
                    bad("core leaves the convex shell or violates interface clearance").into(),
                );
            }
        }
    }
    let m = 2 * count;
    let workers = DENSE_WORKERS.with(core::cell::Cell::get);
    let scratch = lu_complex_blocked_workspace_bytes(m, workers)
        .ok_or_else(|| bad("core-shell LU workspace byte count overflow"))?;
    // During preparation A, forcing and the LU copy coexist. Temporary local
    // contrast blocks are no larger than the LU copy and are dropped first.
    // Retained matrices, original LU descriptors, incident/solution/trace
    // vectors and representation copies fit this conservative payload bound.
    // Geometry maps, allocator/runtime metadata and thread stacks are excluded.
    let entries = m
        .checked_mul(m)
        .and_then(|n| n.checked_mul(3))
        .and_then(|n| n.checked_add(14 * m))
        .ok_or_else(|| bad("core-shell numerical payload count overflow"))?;
    let required_bytes = entries
        .checked_mul(core::mem::size_of::<C64>())
        .and_then(|n| n.checked_add(m * core::mem::size_of::<usize>()))
        .and_then(|n| n.checked_add(scratch))
        .ok_or_else(|| bad("core-shell numerical payload byte count overflow"))?;
    if required_bytes > options.maximum_bytes {
        return Err(HelmholtzError::CacheBudget {
            required_bytes,
            maximum_bytes: options.maximum_bytes,
        }
        .into());
    }
    poll(gate)?;
    Ok(workers)
}
