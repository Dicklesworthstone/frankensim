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

/// A material point in an original trilinear Hex8 cell, including its boundary.
/// Coordinates use the original `[-1,1]^3` chart and physical Q1 interpolation,
/// independently of any assumed mean-dilatation constitutive field.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HexImpulsePoint {
    pub hexahedron: usize,
    pub reference_coordinates: [f64; 3],
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

/// Complete cross-point/direction response of the supported original mass.
#[derive(Clone, Debug)]
pub struct SolidPointMobilityMatrix {
    pub dimension: usize,
    /// Row-major B M_ff^-1 B^T [1/kg]. Columns retain their actual independently
    /// checked solves; entries are never symmetrized or clipped. Consumers own
    /// any additional symmetry, rank or positive-energy admission they require.
    pub matrix_per_kg: Vec<f64>,
    /// One original native mass solve per column, in input order.
    pub solves: Vec<SolveReport>,
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
    let trace = validate(point, tetrahedra, nodes_m.len())?;
    traced_kinematics(trace, nodes_m, displacement_m, velocity_mps, cx)
}

fn traced_kinematics<const N: usize>(
    (nodes, weights): ([usize; N], [f64; N]),
    nodes_m: &[[f64; 3]],
    displacement_m: &[[f64; 3]],
    velocity_mps: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<SolidPointKinematics, InertiaError> {
    if nodes_m.len() != displacement_m.len() || nodes_m.len() != velocity_mps.len() {
        return Err(invalid(
            "complete material-point displacement and velocity required",
        ));
    }
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
    unit_direction(direction)?;
    let impulse = distribution(point, tetrahedra, mass.ndofs() / 3, direction, cx)?;
    response_mobility(mass, held, &impulse, controls, cx)
}

fn response_mobility(
    mass: &ReferenceMass,
    held: &[usize],
    impulse: &[(usize, f64)],
    controls: SolidImpulseControls,
    cx: &Cx<'_>,
) -> Result<SolidPointMobility, InertiaError> {
    let (response, solve) = mass.free_velocity_response(held, impulse, controls, cx)?;
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

fn unit_direction(direction: [f64; 3]) -> Result<(), InertiaError> {
    let norm = direction[0].hypot(direction[1]).hypot(direction[2]);
    if !norm.is_finite() || (norm - 1.0).abs() > 64.0 * f64::EPSILON {
        return Err(invalid(
            "material-point mobility requires a finite unit world direction",
        ));
    }
    Ok(())
}

fn validate_hex(
    point: HexImpulsePoint,
    hexahedra: &[[usize; 8]],
    node_count: usize,
) -> Result<([usize; 8], [f64; 8]), InertiaError> {
    if point
        .reference_coordinates
        .iter()
        .any(|x| !x.is_finite() || !(-1.0..=1.0).contains(x))
    {
        return Err(invalid(
            "Hex8 material-point coordinates must lie in [-1,1]^3",
        ));
    }
    let nodes = *hexahedra
        .get(point.hexahedron)
        .ok_or_else(|| invalid("material point hexahedron is outside the original mesh"))?;
    if nodes.iter().any(|&node| node >= node_count) {
        return Err(invalid("material point node is outside the original mesh"));
    }
    let (weights, _) = crate::linear3::hex::hex8_shapes(point.reference_coordinates);
    let sum: f64 = weights.iter().sum();
    Ok((nodes, weights.map(|w| w / sum)))
}

pub(crate) fn hex_kinematics(
    point: HexImpulsePoint,
    hexahedra: &[[usize; 8]],
    nodes_m: &[[f64; 3]],
    displacement_m: &[[f64; 3]],
    velocity_mps: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<SolidPointKinematics, InertiaError> {
    checkpoint(cx)?;
    traced_kinematics(
        validate_hex(point, hexahedra, nodes_m.len())?,
        nodes_m,
        displacement_m,
        velocity_mps,
        cx,
    )
}

pub(crate) fn hex_distribution(
    point: HexImpulsePoint,
    hexahedra: &[[usize; 8]],
    node_count: usize,
    impulse_world: [f64; 3],
    cx: &Cx<'_>,
) -> Result<[(usize, f64); 24], InertiaError> {
    checkpoint(cx)?;
    if impulse_world.iter().any(|value| !value.is_finite()) {
        return Err(invalid("finite world material-point impulse required"));
    }
    let (nodes, weights) = validate_hex(point, hexahedra, node_count)?;
    Ok(std::array::from_fn(|i| {
        (
            3 * nodes[i / 3] + i % 3,
            weights[i / 3] * impulse_world[i % 3],
        )
    }))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn hex_directional_mobility(
    mass: &ReferenceMass,
    hexahedra: &[[usize; 8]],
    held: &[usize],
    point: HexImpulsePoint,
    direction: [f64; 3],
    controls: SolidImpulseControls,
    cx: &Cx<'_>,
) -> Result<SolidPointMobility, InertiaError> {
    checkpoint(cx)?;
    unit_direction(direction)?;
    let impulse = hex_distribution(point, hexahedra, mass.ndofs() / 3, direction, cx)?;
    response_mobility(mass, held, &impulse, controls, cx)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn mobility_matrix(
    mass: &ReferenceMass,
    tetrahedra: &[[usize; 4]],
    held: &[usize],
    rows: &[(TetImpulsePoint, [f64; 3])],
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

#[allow(clippy::too_many_arguments)]
pub(crate) fn hex_mobility_matrix(
    mass: &ReferenceMass,
    hexahedra: &[[usize; 8]],
    held: &[usize],
    rows: &[(HexImpulsePoint, [f64; 3])],
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
        |point| validate_hex(point, hexahedra, mass.ndofs() / 3),
        cx,
    )
}

/// Store only sparse point rows and one full velocity response at a time.
/// Repeated/same-point rows remain physical cross terms; no diagonal shortcut.
#[allow(clippy::too_many_arguments)]
fn mapped_mobility<P: Copy, const N: usize>(
    mass: &ReferenceMass,
    held: &[usize],
    rows: &[(P, [f64; 3])],
    maximum_rows: usize,
    controls: SolidImpulseControls,
    map: impl Fn(P) -> Result<([usize; N], [f64; N]), InertiaError>,
    cx: &Cx<'_>,
) -> Result<SolidPointMobilityMatrix, InertiaError> {
    checkpoint(cx)?;
    let dimension = rows.len();
    if dimension > maximum_rows {
        return Err(invalid("material-point mobility exceeds its row cap"));
    }
    let entries = dimension
        .checked_mul(dimension)
        .ok_or_else(|| invalid("material-point mobility size overflow"))?;
    let extra = controls
        .solve
        .max_iterations
        .checked_add(1)
        .and_then(|h| h.checked_mul(4))
        .and_then(|h| h.checked_add(6 * N + 64))
        .and_then(|per_row| per_row.checked_mul(dimension))
        .and_then(|rows| rows.checked_add(entries))
        .ok_or_else(|| invalid("material-point mobility workspace overflow"))?;
    admit_point_workspace(mass.ndofs(), extra, controls)?;
    if dimension == 0 {
        super::admit_reference_impulses(mass.ndofs(), held, &[], controls, cx)?;
        return Ok(SolidPointMobilityMatrix {
            dimension,
            matrix_per_kg: Vec::new(),
            solves: Vec::new(),
        });
    }
    let mut loads = Vec::with_capacity(dimension);
    for &(point, direction) in rows {
        checkpoint(cx)?;
        unit_direction(direction)?;
        let (nodes, weights) = map(point)?;
        let mut load = Vec::with_capacity(3 * N);
        for (&node, &weight) in nodes.iter().zip(&weights) {
            for (axis, &component) in direction.iter().enumerate() {
                load.push((3 * node + axis, finite(weight * component)?));
            }
        }
        loads.push(load);
    }
    let mut matrix_per_kg = vec![0.0; entries];
    let mut solves = Vec::with_capacity(dimension);
    for (column, load) in loads.iter().enumerate() {
        checkpoint(cx)?;
        let (response, solve) = mass.free_velocity_response(held, load, controls, cx)?;
        for (row, trace) in loads.iter().enumerate() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            let value = finite(trace.iter().fold(0.0, |sum, &(dof, weight)| {
                weight.mul_add(response[dof / 3][dof % 3], sum)
            }))?;
            if row == column && value < 0.0 {
                return Err(invalid("negative material-point diagonal mobility"));
            }
            matrix_per_kg[row * dimension + column] = value;
        }
        solves.push(solve);
    }
    checkpoint(cx)?;
    Ok(SolidPointMobilityMatrix {
        dimension,
        matrix_per_kg,
        solves,
    })
}

pub(crate) fn distributions(
    tetrahedra: &[[usize; 4]],
    node_count: usize,
    impulses: &[(TetImpulsePoint, [f64; 3])],
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

pub(crate) fn hex_distributions(
    hexahedra: &[[usize; 8]],
    node_count: usize,
    impulses: &[(HexImpulsePoint, [f64; 3])],
    controls: SolidImpulseControls,
    cx: &Cx<'_>,
) -> Result<Vec<(usize, f64)>, InertiaError> {
    mapped_distributions(
        node_count,
        impulses,
        controls,
        |point| validate_hex(point, hexahedra, node_count),
        cx,
    )
}

fn mapped_distributions<P: Copy, const N: usize>(
    node_count: usize,
    impulses: &[(P, [f64; 3])],
    controls: SolidImpulseControls,
    map: impl Fn(P) -> Result<([usize; N], [f64; N]), InertiaError>,
    cx: &Cx<'_>,
) -> Result<Vec<(usize, f64)>, InertiaError> {
    checkpoint(cx)?;
    let ndofs = node_count
        .checked_mul(3)
        .ok_or_else(|| invalid("material-point DOF count overflow"))?;
    let count = impulses
        .len()
        .checked_mul(3 * N)
        .ok_or_else(|| invalid("material-point load count overflow"))?;
    let extra = count
        .checked_mul(2)
        .ok_or_else(|| invalid("material-point load workspace overflow"))?;
    admit_point_workspace(ndofs, extra, controls)?;
    let mut output = Vec::with_capacity(count);
    for &(point, impulse) in impulses {
        checkpoint(cx)?;
        if impulse.iter().any(|v| !v.is_finite()) {
            return Err(invalid("finite world material-point impulses required"));
        }
        let (nodes, weights) = map(point)?;
        for (&node, &weight) in nodes.iter().zip(&weights) {
            for (axis, &component) in impulse.iter().enumerate() {
                output.push((3 * node + axis, finite(weight * component)?));
            }
        }
    }
    // Preserve individual contributions until the original mass owner forms
    // each simultaneous resultant with compensated accumulation.
    checkpoint(cx)?;
    Ok(output)
}

fn admit_point_workspace(
    n: usize,
    extra: usize,
    controls: SolidImpulseControls,
) -> Result<(), InertiaError> {
    let required = super::impulse::workspace_components(n, controls)?
        .checked_add(extra)
        .ok_or_else(|| invalid("material-point workspace overflow"))?;
    if required > controls.maximum_workspace_components {
        return Err(invalid(
            "material-point operation exceeds its workspace component cap",
        ));
    }
    Ok(())
}
