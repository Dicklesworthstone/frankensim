//! Exact affine elimination of mechanical constraints and conducting electrodes.

use super::{ElectrodeCondition, PiezoelectricElectrode, Result, invalid, poll};
use fs_exec::Cx;
use std::collections::BTreeMap;

#[derive(Clone)]
pub(super) struct BoundaryMap {
    pub index: Vec<Option<usize>>,
    pub prescribed: Vec<f64>,
    pub rhs: Vec<f64>,
    pub voltage_reference: Vec<f64>,
    /// First prescribed-voltage electrode in each connected electric component.
    /// Retaining its identity lets transient sources change their gauge without
    /// rebuilding the mesh connectivity or confusing equal voltage values.
    pub voltage_reference_electrode: Vec<usize>,
}

fn root(parent: &mut [usize], mut node: usize) -> usize {
    while parent[node] != node {
        parent[node] = parent[parent[node]];
        node = parent[node];
    }
    node
}
fn join(parent: &mut [usize], a: usize, b: usize) {
    let a = root(parent, a);
    let b = root(parent, b);
    if a != b {
        parent[b] = a;
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare_boundary(
    nodes: &[[f64; 3]],
    tetrahedra: &[[usize; 4]],
    prescribed_displacements: &[(usize, f64)],
    electrodes: &[PiezoelectricElectrode<'_>],
    forces: &[[f64; 3]],
    charges: &[f64],
    cx: &Cx<'_>,
) -> Result<BoundaryMap> {
    let n = nodes.len();
    let full_dofs = n
        .checked_mul(4)
        .ok_or_else(|| invalid("coupled DOF count overflows"))?;
    if prescribed_displacements.len() > 3 * n || electrodes.len() > n {
        return Err(invalid("boundary data exceeds original coordinate count"));
    }
    let mut prescribed = vec![0.0; full_dofs];
    let mut fixed = vec![false; 3 * n];
    for &(dof, value) in prescribed_displacements {
        poll(cx)?;
        if dof >= 3 * n || fixed[dof] || !value.is_finite() {
            return Err(invalid(
                "finite unique displacement constraints in original DOF range required",
            ));
        }
        fixed[dof] = true;
        prescribed[dof] = value;
    }
    let mut parent = (0..n).collect::<Vec<_>>();
    let mut face_owners = BTreeMap::<[usize; 3], usize>::new();
    for tet in tetrahedra {
        poll(cx)?;
        for &node in &tet[1..] {
            join(&mut parent, tet[0], node);
        }
        for skip in 0..4 {
            let mut face = [0; 3];
            let mut index = 0;
            for (local, &node) in tet.iter().enumerate() {
                if local != skip {
                    face[index] = node;
                    index += 1;
                }
            }
            face.sort_unstable();
            let count = face_owners.entry(face).or_insert(0);
            *count += 1;
            if *count > 2 {
                return Err(invalid(
                    "nonmanifold tetrahedral face has more than two owners",
                ));
            }
        }
    }
    let exterior = face_owners
        .into_iter()
        .filter_map(|(face, count)| (count == 1).then_some(face))
        .collect::<Vec<_>>();
    let mut electrode_of = vec![None; n];
    for (id, electrode) in electrodes.iter().enumerate() {
        poll(cx)?;
        let value = match electrode.condition {
            ElectrodeCondition::Voltage(v) | ElectrodeCondition::Charge(v) => v,
        };
        if electrode.nodes.is_empty() || electrode.nodes.len() > n || !value.is_finite() {
            return Err(invalid(
                "each electrode requires a nonempty bounded face trace and finite condition",
            ));
        }
        for &node in electrode.nodes {
            if node >= n || electrode_of[node].is_some() {
                return Err(invalid(
                    "electrode vertices must be in range, unique and disjoint",
                ));
            }
            electrode_of[node] = Some(id);
        }
        let mut supported = vec![false; n];
        for face in &exterior {
            poll(cx)?;
            if face.iter().all(|&node| electrode_of[node] == Some(id)) {
                for &node in face {
                    supported[node] = true;
                }
            }
        }
        if electrode.nodes.iter().any(|&node| !supported[node]) {
            return Err(invalid(
                "electrode must be a union of complete exterior faces, not point contacts",
            ));
        }
        // A physical wire joining disconnected patches also identifies their
        // potential gauges, so count voltage references on this joined graph.
        for &node in &electrode.nodes[1..] {
            join(&mut parent, electrode.nodes[0], node);
        }
    }
    let components = (0..n)
        .map(|node| root(&mut parent, node))
        .collect::<Vec<_>>();
    let mut references = vec![None; n];
    for (id, electrode) in electrodes.iter().enumerate() {
        if let ElectrodeCondition::Voltage(voltage) = electrode.condition {
            references[components[electrode.nodes[0]]].get_or_insert((voltage, id));
        }
    }
    let mut voltage_reference = Vec::with_capacity(n);
    let mut voltage_reference_electrode = Vec::with_capacity(n);
    for &component in &components {
        let (voltage, id) = references[component].ok_or_else(|| {
            invalid(
                "each connected dielectric/electrode component needs an explicit voltage reference",
            )
        })?;
        voltage_reference.push(voltage);
        voltage_reference_electrode.push(id);
    }
    let mut index = vec![None; full_dofs];
    let mut next = 0;
    for dof in 0..3 * n {
        if !fixed[dof] {
            index[dof] = Some(next);
            next += 1;
        }
    }
    let mut electrode_index = vec![None; electrodes.len()];
    for node in 0..n {
        if let Some(id) = electrode_of[node] {
            match electrodes[id].condition {
                ElectrodeCondition::Voltage(voltage) => {
                    let relative = voltage - voltage_reference[node];
                    if !relative.is_finite() {
                        return Err(invalid("electrode voltage difference overflows"));
                    }
                    prescribed[3 * n + node] = relative;
                }
                ElectrodeCondition::Charge(_) => {
                    let coordinate = *electrode_index[id].get_or_insert_with(|| {
                        let coordinate = next;
                        next += 1;
                        coordinate
                    });
                    index[3 * n + node] = Some(coordinate);
                }
            }
        } else {
            index[3 * n + node] = Some(next);
            next += 1;
        }
    }
    let mut rhs = vec![0.0; next];
    for node in 0..n {
        poll(cx)?;
        for axis in 0..3 {
            if let Some(row) = index[3 * node + axis] {
                rhs[row] += forces[node][axis];
            }
        }
        if let Some(row) = index[3 * n + node] {
            rhs[row] -= charges[node];
        }
    }
    for (id, electrode) in electrodes.iter().enumerate() {
        if let ElectrodeCondition::Charge(charge) = electrode.condition {
            rhs[electrode_index[id].expect("charge electrode has coordinate")] -= charge;
        }
    }
    if rhs.iter().any(|v| !v.is_finite()) {
        return Err(invalid("aggregated nodal load overflows"));
    }
    Ok(BoundaryMap {
        index,
        prescribed,
        rhs,
        voltage_reference,
        voltage_reference_electrode,
    })
}
