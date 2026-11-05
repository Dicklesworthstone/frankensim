//! Boundary-rooted circulation basis for exact incompressible junctions.
//! Integer through-flow/cycle columns satisfy incidence before any floating
//! arithmetic. Orthonormalizing sqrt(L)*F preserves the ORIGINAL kinetic metric.
//! No pressure storage, penalty, fitted resistance or discarded loop mode.
use fs_math::det;

pub(super) const ROUNDOFF: f64 = 4096.0 * f64::EPSILON;

pub(super) fn norm(values: &[f64]) -> f64 {
    values.iter().fold(0.0, |n, x| det::hypot(n, *x))
}

pub(super) struct ContinuityBasis {
    pub columns: usize,
    /// Row-major branch-mean energy map x_mean = U*y, with U^T U = I.
    pub energy_rows: Vec<f64>,
    pub boundary_slot: Vec<Option<usize>>,
    /// One pressure gauge per physical connected component.
    pub pressure_reference: Vec<usize>,
}

fn root(parent: &[usize], mut node: usize) -> usize {
    while parent[node] != node {
        node = parent[node];
    }
    node
}

pub(super) fn compile(
    node_count: usize,
    ends: &[[usize; 2]],
    pressure_nodes: &[usize],
    mean_scales: &[f64],
) -> Result<ContinuityBasis, &'static str> {
    let edges = ends.len();
    if node_count == 0
        || node_count > 256
        || edges == 0
        || edges > 128
        || pressure_nodes.is_empty()
        || pressure_nodes.len() > node_count
        || mean_scales.len() != edges
        || mean_scales
            .iter()
            .any(|s| !s.is_finite() || *s <= 0.0 || !(1.0 / s).is_finite())
        || ends
            .iter()
            .any(|e| e[0] >= node_count || e[1] >= node_count || e[0] == e[1])
    {
        return Err(
            "bounded nonempty topology, distinct branch ends and finite inertia scales required",
        );
    }
    let mut boundary_slot = vec![None; node_count];
    for (slot, &node) in pressure_nodes.iter().enumerate() {
        if node >= node_count || boundary_slot[node].is_some() {
            return Err("pressure nodes must be distinct original topology nodes");
        }
        boundary_slot[node] = Some(slot);
    }
    let mut components: Vec<_> = (0..node_count).collect();
    for &[a, b] in ends {
        let (a, b) = (root(&components, a), root(&components, b));
        if a != b {
            components[b] = a;
        }
    }
    let mut pressure_reference = Vec::with_capacity(pressure_nodes.len());
    for &node in pressure_nodes {
        let component = root(&components, node);
        let first = pressure_nodes
            .iter()
            .position(|&other| root(&components, other) == component)
            .ok_or("missing component pressure gauge")?;
        pressure_reference.push(first);
    }
    for node in 0..node_count {
        let component = root(&components, node);
        if !pressure_nodes
            .iter()
            .any(|&p| root(&components, p) == component)
        {
            return Err(
                "every physical component needs a pressure node; one may be a closed-loop gauge",
            );
        }
    }
    // Merge pressure roots only for spanning-forest selection. They retain
    // independent physical pressures and are NEVER identified in the port map.
    let mut forest: Vec<_> = (0..node_count).collect();
    for &node in &pressure_nodes[1..] {
        forest[node] = pressure_nodes[0];
    }
    let mut tree = vec![false; edges];
    let mut adjacency = vec![Vec::new(); node_count];
    let mut tree_count = 0;
    for (edge, &[a, b]) in ends.iter().enumerate() {
        let (ra, rb) = (root(&forest, a), root(&forest, b));
        if ra != rb {
            forest[rb] = ra;
            tree[edge] = true;
            tree_count += 1;
            adjacency[a].push((b, edge));
            adjacency[b].push((a, edge));
        }
    }
    if tree_count != node_count - pressure_nodes.len() {
        return Err("internal continuity rows have unresolved topological rank");
    }
    let mut parent_edge = vec![None; node_count];
    let mut parent_node = vec![0; node_count];
    let mut visited = vec![false; node_count];
    let mut queue = pressure_nodes.to_vec();
    for &node in pressure_nodes {
        visited[node] = true;
    }
    let mut cursor = 0;
    while cursor < queue.len() {
        let node = queue[cursor];
        cursor += 1;
        for &(next, edge) in &adjacency[node] {
            if !visited[next] {
                visited[next] = true;
                parent_node[next] = node;
                parent_edge[next] = Some(edge);
                queue.push(next);
            }
        }
    }
    if visited.iter().any(|v| !*v) {
        return Err("boundary-rooted physical spanning forest is incomplete");
    }
    let columns = edges - tree_count;
    let mut candidates = Vec::with_capacity(columns);
    for chord in 0..edges {
        if tree[chord] {
            continue;
        }
        let mut flow = vec![0.0; edges];
        flow[chord] = 1.0;
        for (mut node, amount) in [(ends[chord][0], 1.0), (ends[chord][1], -1.0)] {
            while let Some(edge) = parent_edge[node] {
                let parent = parent_node[node];
                let direction = if ends[edge][0] == parent { 1.0 } else { -1.0 };
                flow[edge] += amount * direction;
                node = parent;
            }
        }
        for (value, scale) in flow.iter_mut().zip(mean_scales) {
            *value *= scale;
        }
        let scale = norm(&flow);
        if !scale.is_finite() || scale <= 0.0 {
            return Err("circulation energy scale is unrepresentable");
        }
        for value in &mut flow {
            *value /= scale;
        }
        candidates.push(flow);
    }
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(columns);
    // Deterministic pivoted, twice reorthogonalized modified Gram-Schmidt.
    // A tiny remaining column is a refusal, never permission to drop a mode.
    for k in 0..columns {
        let mut pivot = k;
        let mut largest = norm(&candidates[k]);
        for (j, column) in candidates.iter().enumerate().skip(k + 1) {
            let length = norm(column);
            if length > largest {
                pivot = j;
                largest = length;
            }
        }
        candidates.swap(k, pivot);
        let mut column = candidates[k].clone();
        for _ in 0..2 {
            for previous in &basis {
                let dot: f64 = previous.iter().zip(&column).map(|(a, b)| a * b).sum();
                for (value, prior) in column.iter_mut().zip(previous) {
                    *value -= dot * prior;
                }
            }
        }
        let length = norm(&column);
        if !length.is_finite() || length <= ROUNDOFF * edges as f64 {
            return Err(
                "complete circulation basis is numerically unresolved; no mode was discarded",
            );
        }
        for value in &mut column {
            *value /= length;
        }
        for remaining in candidates.iter_mut().skip(k + 1) {
            for _ in 0..2 {
                let dot: f64 = column
                    .iter()
                    .zip(remaining.iter())
                    .map(|(a, b)| a * b)
                    .sum();
                for (value, q) in remaining.iter_mut().zip(&column) {
                    *value -= dot * q;
                }
            }
        }
        basis.push(column);
    }
    let mut energy_rows = vec![0.0; edges * columns];
    for (j, column) in basis.iter().enumerate() {
        for (edge, value) in column.iter().enumerate() {
            energy_rows[edge * columns + j] = *value;
        }
        for (k, other) in basis.iter().enumerate() {
            let dot: f64 = column.iter().zip(other).map(|(a, b)| a * b).sum();
            if (dot - if j == k { 1.0 } else { 0.0 }).abs() > ROUNDOFF * edges as f64 {
                return Err("circulation basis does not preserve original kinetic storage");
            }
        }
        let mut outflow = vec![0.0; node_count];
        let mut scale = vec![0.0; node_count];
        for (edge, &[a, b]) in ends.iter().enumerate() {
            let flow = column[edge] / mean_scales[edge];
            if !flow.is_finite() {
                return Err("circulation flow observation is unrepresentable");
            }
            outflow[a] += flow;
            outflow[b] -= flow;
            scale[a] += flow.abs();
            scale[b] += flow.abs();
        }
        for node in 0..node_count {
            if boundary_slot[node].is_none()
                && outflow[node].abs() > ROUNDOFF * edges as f64 * scale[node]
            {
                return Err("energy projection does not retain internal incompressibility");
            }
        }
    }
    Ok(ContinuityBasis {
        columns,
        energy_rows,
        boundary_slot,
        pressure_reference,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn verify(nodes: usize, ends: &[[usize; 2]], ports: &[usize], scales: &[f64], expected: usize) {
        let b = compile(nodes, ends, ports, scales).unwrap();
        assert_eq!(b.columns, expected);
        for j in 0..b.columns {
            let mut balance = vec![0.0; nodes];
            for (e, &[a, z]) in ends.iter().enumerate() {
                let q = b.energy_rows[e * b.columns + j] / scales[e];
                balance[a] += q;
                balance[z] -= q;
            }
            for node in 0..nodes {
                if !ports.contains(&node) {
                    assert!(balance[node].abs() < 2e-13);
                }
            }
            for k in 0..b.columns {
                let dot: f64 = (0..ends.len())
                    .map(|e| b.energy_rows[e * b.columns + j] * b.energy_rows[e * b.columns + k])
                    .sum();
                assert!((dot - if j == k { 1.0 } else { 0.0 }).abs() < 3e-14);
            }
        }
    }

    #[test]
    fn series_parallel_fork_and_mesh_keep_all_physical_modes() {
        verify(3, &[[0, 1], [1, 2]], &[0, 2], &[2.0, 5.0], 1);
        verify(2, &[[0, 1], [0, 1]], &[0, 1], &[2.0, 5.0], 2);
        verify(
            4,
            &[[0, 1], [1, 2], [1, 3]],
            &[0, 2, 3],
            &[2.0, 5.0, 0.7],
            2,
        );
        verify(
            4,
            &[[0, 1], [1, 2], [0, 3], [3, 2], [1, 3]],
            &[0, 2],
            &[2.0, 5.0, 0.7, 1.1, 3.0],
            3,
        );
    }

    #[test]
    fn closed_loop_has_one_gauge_and_retains_circulation() {
        verify(3, &[[0, 1], [1, 2], [2, 0]], &[0], &[2.0, 5.0, 3.0], 1);
        let b = compile(3, &[[0, 1], [1, 2], [2, 0]], &[0], &[2.0, 5.0, 3.0]).unwrap();
        let flow = b.energy_rows[0] / 2.0;
        assert!((b.energy_rows[1] / 5.0 - flow).abs() < 1e-15);
        assert!((b.energy_rows[2] / 3.0 - flow).abs() < 1e-15);
    }

    #[test]
    fn disconnected_components_keep_separate_gauges_and_dead_end_has_no_mean_mode() {
        let b = compile(4, &[[0, 1], [2, 3]], &[0, 1, 2, 3], &[2.0, 3.0]).unwrap();
        assert_eq!(b.pressure_reference, vec![0, 0, 2, 2]);
        verify(3, &[[0, 1], [1, 2]], &[0], &[2.0, 3.0], 0);
        assert!(compile(4, &[[0, 1], [2, 3]], &[0], &[2.0, 3.0]).is_err());
    }

    #[test]
    fn reversed_orientations_and_wide_inertia_scales_preserve_weighted_continuity() {
        verify(
            4,
            &[[1, 0], [2, 1], [3, 0], [2, 3], [3, 1]],
            &[0, 2],
            &[0.001, 3.0, 70.0, 0.1, 2.0],
            3,
        );
        assert!(compile(2, &[[0, 1]], &[0, 0], &[1.0]).is_err());
        assert!(compile(2, &[[0, 1]], &[2], &[1.0]).is_err());
        assert!(compile(2, &[[0, 1]], &[0], &[f64::NAN]).is_err());
        assert!(compile(2, &[[0, 0]], &[0], &[1.0]).is_err());
    }
}
