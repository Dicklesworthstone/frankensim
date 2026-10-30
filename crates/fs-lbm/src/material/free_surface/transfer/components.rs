//! No impulse, torque or numerical-energy cancellation between detached pieces.

use super::*;

/// One original face-connected tetrahedral piece and its native ingress.
#[derive(Clone, Debug, PartialEq)]
pub struct TetLiquidComponentTransfer {
    /// Indices into the unchanged source tetrahedron list, in input order.
    /// Components are ordered by their lowest original tetrahedron index.
    pub source_tetrahedra: Vec<usize>,
    /// Occupied original VOF cell coordinates belonging to this piece at ingress.
    /// A cell cannot alias another disconnected source piece.
    pub target_cells: Vec<[usize; 3]>,
    pub source: LiquidTransferMoments,
    pub target: LiquidTransferMoments,
    pub centroid_shift_m: [f64; 3],
    pub kinetic_energy_change_j: f64,
    pub gravitational_energy_change_j: f64,
}

#[derive(Default)]
pub(super) struct Component {
    pub(super) source_tetrahedra: Vec<usize>,
    pub(super) cell_indices: Vec<usize>,
    pub(super) source: LiquidTransferMoments,
}

fn root(parent: &mut [usize], mut node: usize) -> usize {
    while parent[node] != node {
        parent[node] = parent[parent[node]];
        node = parent[node];
    }
    node
}

pub(super) fn partition(
    source: TetrahedralLiquidState<'_>,
    cx: &Cx<'_>,
) -> Result<(Vec<usize>, Vec<Component>), TetLiquidTransferError> {
    let mut parent: Vec<_> = (0..source.tetrahedra.len()).collect();
    let mut faces = BTreeMap::new();
    for (index, tet) in source.tetrahedra.iter().enumerate() {
        checkpoint(cx)?;
        if tet
            .iter()
            .any(|&node| node as usize >= source.positions_m.len())
        {
            return Err(invalid("tetrahedron index lies outside the source field"));
        }
        let mut distinct = *tet;
        distinct.sort_unstable();
        if distinct.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(invalid("a source tetrahedron repeats a vertex"));
        }
        for excluded in 0..4 {
            let mut face = [0; 3];
            for (slot, node) in tet
                .iter()
                .enumerate()
                .filter(|(slot, _)| *slot != excluded)
                .map(|(_, node)| node)
                .enumerate()
            {
                face[slot] = *node;
            }
            face.sort_unstable();
            let (other, shared) = faces.entry(face).or_insert((index, false));
            if *other != index {
                if *shared {
                    return Err(invalid(
                        "source tetrahedral face has more than two incident cells",
                    ));
                }
                *shared = true;
                let a = root(&mut parent, index);
                let b = root(&mut parent, *other);
                parent[a.max(b)] = a.min(b);
            }
        }
    }
    let mut number = BTreeMap::new();
    let mut components = Vec::<Component>::new();
    let mut assignment = Vec::with_capacity(parent.len());
    for index in 0..parent.len() {
        checkpoint(cx)?;
        let root = root(&mut parent, index);
        let next = number.len();
        let component = *number.entry(root).or_insert(next);
        if component == components.len() {
            components.push(Component::default());
        }
        components[component].source_tetrahedra.push(index);
        assignment.push(component);
    }
    Ok((assignment, components))
}

pub(super) fn correct(
    initial: &mut [MaterialFreeSurfaceInitialCell3],
    components: &[Component],
    mass_unit: f64,
    dx: f64,
    limits: TetLiquidTransferLimits,
    cx: &Cx<'_>,
) -> Result<(), TetLiquidTransferError> {
    for component in components {
        checkpoint(cx)?;
        let mut cells: Vec<_> = component
            .cell_indices
            .iter()
            .map(|&index| initial[index])
            .collect();
        correct_momenta(&mut cells, mass_unit, dx, &component.source, limits, cx)?;
        for (&index, cell) in component.cell_indices.iter().zip(cells) {
            initial[index] = cell;
        }
    }
    Ok(())
}

pub(super) fn reports(
    components: Vec<Component>,
    actual: Vec<LiquidTransferMoments>,
    initial: &[MaterialFreeSurfaceInitialCell3],
    gravity: [f64; 3],
    limits: TetLiquidTransferLimits,
    cx: &Cx<'_>,
) -> Result<Vec<TetLiquidComponentTransfer>, TetLiquidTransferError> {
    let mut reports = Vec::with_capacity(components.len());
    let (mut mass, mut momentum, mut angular, mut energy) = (0.0, 0.0, 0.0, 0.0);
    for (component, target) in components.into_iter().zip(actual) {
        checkpoint(cx)?;
        let (centroid_shift_m, kinetic_energy_change_j, gravitational_energy_change_j) =
            check_projection(&component.source, &target, gravity, limits)?;
        mass += (target.mass_kg - component.source.mass_kg).abs();
        momentum += norm(sub(
            target.momentum_kg_m_s,
            component.source.momentum_kg_m_s,
        ));
        angular += norm(sub(
            target.angular_momentum_kg_m2_s,
            component.source.angular_momentum_kg_m2_s,
        ));
        energy += kinetic_energy_change_j.abs() + gravitational_energy_change_j.abs();
        reports.push(TetLiquidComponentTransfer {
            source_tetrahedra: component.source_tetrahedra,
            target_cells: component
                .cell_indices
                .iter()
                .map(|&index| initial[index].coordinates)
                .collect(),
            source: component.source,
            target,
            centroid_shift_m,
            kinetic_energy_change_j,
            gravitational_energy_change_j,
        });
    }
    // These are ONE original work/error budget. Opposite discrepancies cannot
    // cancel between physically detached pieces, or each spend the full budget.
    if ![mass, momentum, angular, energy]
        .into_iter()
        .all(f64::is_finite)
        || mass > limits.mass_tolerance_kg
        || momentum > limits.momentum_tolerance_kg_m_s
        || angular > limits.angular_momentum_tolerance_kg_m2_s
        || energy > limits.maximum_projection_energy_change_j
    {
        return Err(invalid(
            "component transfer exceeds the shared absolute conservation/error budget",
        ));
    }
    Ok(reports)
}
