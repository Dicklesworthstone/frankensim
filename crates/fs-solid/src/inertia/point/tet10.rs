//! Signed quadratic material traces sharing the original mass response.

use super::{
    Cx, InertiaError, ReferenceMass, SolidImpulseControls, SolidPointKinematics,
    SolidPointMobility, SolidPointMobilityMatrix, SolveReport, barycentric, checkpoint, invalid,
    mapped_distributions, mapped_mobility, mapped_mobility_with, response_mobility,
    traced_kinematics, unit_direction,
};
use crate::linear3::tet10::tet10_shapes;

/// A material point in an original quadratic tetrahedron, including its boundary.
/// These are four nonnegative, unit-sum natural coordinates, not the ten signed
/// displacement shape weights. All corner and midside nodes contribute to the
/// current isoparametric position, velocity and work-conjugate impulse.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tet10ImpulsePoint {
    pub tetrahedron: usize,
    pub barycentric: [f64; 4],
}

fn validate(
    point: Tet10ImpulsePoint,
    tetrahedra: &[[usize; 10]],
    node_count: usize,
) -> Result<([usize; 10], [f64; 10]), InertiaError> {
    let lambda = barycentric(point.barycentric)?;
    let nodes = *tetrahedra
        .get(point.tetrahedron)
        .ok_or_else(|| invalid("material point tetrahedron is outside the original mesh"))?;
    if nodes.iter().any(|&node| node >= node_count) {
        return Err(invalid("material point node is outside the original mesh"));
    }
    Ok((nodes, tet10_shapes(lambda).0))
}

pub(crate) fn kinematics(
    point: Tet10ImpulsePoint,
    tetrahedra: &[[usize; 10]],
    nodes_m: &[[f64; 3]],
    displacement_m: &[[f64; 3]],
    velocity_mps: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<SolidPointKinematics, InertiaError> {
    checkpoint(cx)?;
    let trace = validate(point, tetrahedra, nodes_m.len())?;
    traced_kinematics(trace, nodes_m, displacement_m, velocity_mps, cx)
}

/// B^T J retains negative corner weights and the complete curved position map.
pub(crate) fn distribution(
    point: Tet10ImpulsePoint,
    tetrahedra: &[[usize; 10]],
    node_count: usize,
    impulse_world: [f64; 3],
    cx: &Cx<'_>,
) -> Result<[(usize, f64); 30], InertiaError> {
    checkpoint(cx)?;
    if impulse_world.iter().any(|value| !value.is_finite()) {
        return Err(invalid("finite world material-point impulse required"));
    }
    let (nodes, weights) = validate(point, tetrahedra, node_count)?;
    Ok(std::array::from_fn(|i| {
        (
            3 * nodes[i / 3] + i % 3,
            weights[i / 3] * impulse_world[i % 3],
        )
    }))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn directional_mobility(
    mass: &ReferenceMass,
    tetrahedra: &[[usize; 10]],
    held: &[usize],
    point: Tet10ImpulsePoint,
    direction: [f64; 3],
    controls: SolidImpulseControls,
    cx: &Cx<'_>,
) -> Result<SolidPointMobility, InertiaError> {
    checkpoint(cx)?;
    unit_direction(direction)?;
    let impulse = distribution(point, tetrahedra, mass.ndofs() / 3, direction, cx)?;
    response_mobility(mass, held, &impulse, controls, cx)
}

/// Full cross-point response through the same signed quadratic traces.
/// The shared owner retains one full velocity response per column at a time.
#[allow(clippy::too_many_arguments)]
pub(crate) fn mobility_matrix(
    mass: &ReferenceMass,
    tetrahedra: &[[usize; 10]],
    held: &[usize],
    rows: &[(Tet10ImpulsePoint, [f64; 3])],
    maximum_rows: usize,
    controls: SolidImpulseControls,
    cx: &Cx<'_>,
) -> Result<SolidPointMobilityMatrix, InertiaError> {
    mapped_mobility(
        mass,
        held,
        rows,
        maximum_rows,
        controls,
        |point| validate(point, tetrahedra, mass.ndofs() / 3),
        cx,
    )
}

/// Share the original quadratic traces with a prepared constrained response.
#[allow(clippy::too_many_arguments)]
pub(crate) fn mobility_matrix_with<E>(
    node_count: usize,
    tetrahedra: &[[usize; 10]],
    rows: &[(Tet10ImpulsePoint, [f64; 3])],
    maximum_rows: usize,
    controls: SolidImpulseControls,
    respond: impl FnMut(&[(usize, f64)]) -> Result<(Vec<[f64; 3]>, SolveReport), E>,
    map_error: impl Fn(InertiaError) -> E,
    cx: &Cx<'_>,
) -> Result<SolidPointMobilityMatrix, E> {
    let ndofs = node_count
        .checked_mul(3)
        .ok_or_else(|| map_error(invalid("material-point DOF count overflow")))?;
    mapped_mobility_with(
        ndofs,
        rows,
        maximum_rows,
        controls,
        |point| validate(point, tetrahedra, node_count),
        respond,
        map_error,
        cx,
    )
}

/// Keep contributions separate until native compensated nodal accumulation.
pub(crate) fn distributions(
    tetrahedra: &[[usize; 10]],
    node_count: usize,
    impulses: &[(Tet10ImpulsePoint, [f64; 3])],
    controls: SolidImpulseControls,
    cx: &Cx<'_>,
) -> Result<Vec<(usize, f64)>, InertiaError> {
    mapped_distributions(
        node_count,
        impulses,
        controls,
        |point| validate(point, tetrahedra, node_count),
        cx,
    )
}

#[cfg(test)]
mod tests;
