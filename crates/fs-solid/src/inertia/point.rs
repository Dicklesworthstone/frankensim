//! One material trace and its transpose for supported finite-solid impacts.
use fs_solver::SolveReport;

use super::{Cx, InertiaError, ReferenceMass, SolidImpulseControls, checkpoint, finite, invalid};

/// A material point in an original P1 tetrahedron, including faces and vertices.
/// Each use checks the cell index and nonnegative unit-sum weights. A sum defect
/// within 64 machine epsilons is normalized identically for every trace/load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TetImpulsePoint {
    pub tetrahedron: usize,
    pub weights: [f64; 4],
}

/// Current physical point coordinates and its work-conjugate velocity trace.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolidPointKinematics {
    pub position_m: [f64; 3],
    pub velocity_mps: [f64; 3],
}

/// Directional mobility of the supported original consistent mass.
#[derive(Clone, Debug)]
pub struct SolidPointMobility {
    /// n^T B M_ff^-1 B^T n [1/kg]; zero for a fully restrained trace.
    pub inverse_mass_per_kg: f64,
    /// Native checked-CG evidence, including the recomputed true residual.
    pub solve: SolveReport,
}

fn validate(
    point: TetImpulsePoint,
    tetrahedra: &[[usize; 4]],
    node_count: usize,
) -> Result<([usize; 4], [f64; 4]), InertiaError> {
    let sum: f64 = point.weights.iter().sum();
    if point
        .weights
        .iter()
        .any(|w| !w.is_finite() || !(0.0..=1.0).contains(w))
        || !sum.is_finite()
        || (sum - 1.0).abs() > 64.0 * f64::EPSILON
    {
        return Err(invalid(
            "material point requires nonnegative unit-sum barycentric weights",
        ));
    }
    let nodes = *tetrahedra
        .get(point.tetrahedron)
        .ok_or_else(|| invalid("material point tetrahedron is outside the original mesh"))?;
    if nodes.iter().any(|&node| node >= node_count) {
        return Err(invalid("material point node is outside the original mesh"));
    }
    Ok((nodes, point.weights.map(|w| w / sum)))
}

pub(crate) fn kinematics(
    point: TetImpulsePoint,
    tetrahedra: &[[usize; 4]],
    nodes_m: &[[f64; 3]],
    displacement_m: &[[f64; 3]],
    velocity_mps: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<SolidPointKinematics, InertiaError> {
    checkpoint(cx)?;
    if nodes_m.len() != displacement_m.len() || nodes_m.len() != velocity_mps.len() {
        return Err(invalid(
            "complete material-point displacement and velocity required",
        ));
    }
    let (nodes, weights) = validate(point, tetrahedra, nodes_m.len())?;
    let mut position_m = [0.0; 3];
    let mut velocity = [0.0; 3];
    for (&node, &weight) in nodes.iter().zip(&weights) {
        for axis in 0..3 {
            let x = finite(nodes_m[node][axis] + displacement_m[node][axis])?;
            position_m[axis] = finite(weight.mul_add(x, position_m[axis]))?;
            velocity[axis] = finite(weight.mul_add(velocity_mps[node][axis], velocity[axis]))?;
        }
    }
    checkpoint(cx)?;
    Ok(SolidPointKinematics {
        position_m,
        velocity_mps: velocity,
    })
}

/// B^T J uses exactly the same admitted shape weights as the point trace.
pub(crate) fn distribution(
    point: TetImpulsePoint,
    tetrahedra: &[[usize; 4]],
    node_count: usize,
    impulse_world: [f64; 3],
    cx: &Cx<'_>,
) -> Result<[(usize, f64); 12], InertiaError> {
    checkpoint(cx)?;
    if impulse_world.iter().any(|value| !value.is_finite()) {
        return Err(invalid("finite world material-point impulse required"));
    }
    let (nodes, weights) = validate(point, tetrahedra, node_count)?;
    let entries = std::array::from_fn(|i| {
        (
            3 * nodes[i / 3] + i % 3,
            weights[i / 3] * impulse_world[i % 3],
        )
    });
    Ok(entries)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn directional_mobility(
    mass: &ReferenceMass,
    tetrahedra: &[[usize; 4]],
    held: &[usize],
    point: TetImpulsePoint,
    direction: [f64; 3],
    controls: SolidImpulseControls,
    cx: &Cx<'_>,
) -> Result<SolidPointMobility, InertiaError> {
    checkpoint(cx)?;
    let norm = direction[0].hypot(direction[1]).hypot(direction[2]);
    if !norm.is_finite() || (norm - 1.0).abs() > 64.0 * f64::EPSILON {
        return Err(invalid(
            "material-point mobility requires a finite unit world direction",
        ));
    }
    let impulse = distribution(point, tetrahedra, mass.ndofs() / 3, direction, cx)?;
    let (response, solve) = mass.free_velocity_response(held, &impulse, controls, cx)?;
    let inverse_mass_per_kg = finite(impulse.iter().fold(0.0, |sum, &(dof, load)| {
        load.mul_add(response[dof / 3][dof % 3], sum)
    }))?;
    if inverse_mass_per_kg < 0.0 {
        return Err(invalid(
            "negative material-point mobility from original positive mass",
        ));
    }
    checkpoint(cx)?;
    Ok(SolidPointMobility {
        inverse_mass_per_kg,
        solve,
    })
}
