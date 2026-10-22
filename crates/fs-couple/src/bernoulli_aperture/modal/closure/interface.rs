//! Admission and original P1 traces for independently meshed opposing blades.
//! No registration, extrapolation, remeshing or modal renormalization occurs.

use super::{AcousticRealizeError, invalid, nonlinear};
use fs_plate::{
    PlateMesh,
    loading::{PlateLoadBudget, PlateTraceStencil},
};
use std::collections::{BTreeMap, BTreeSet};

type Segment = [[f64; 2]; 2];

/// Both meshes must cover the same simple planar blade in the same frame. The
/// original first mesh remains the explicitly authored integration surface.
pub(super) fn traces(
    interface: &PlateMesh,
    slit: &[[usize; 2]],
    source: &PlateMesh,
    source_slit: &[[usize; 2]],
    budget: PlateLoadBudget,
    max_terms: usize,
) -> Result<Vec<PlateTraceStencil>, AcousticRealizeError> {
    let terms = interface
        .nodes
        .len()
        .checked_mul(source.tris.len())
        .ok_or_else(|| invalid("modal interface interpolation work overflow"))?;
    if terms > max_terms {
        return Err(invalid("modal interface exceeds its source projection budget"));
    }
    let first = boundary(
        interface,
        PlateLoadBudget {
            max_nodes: 512,
            max_triangles: 2048,
        },
    )?;
    let second = boundary(source, budget)?;
    same_cover(&first, &second, true)?;
    let slit_segments = segments(interface, slit)?;
    let source_segments = segments(source, source_slit)?;
    same_cover(&slit_segments, &source_segments, false)?;
    interface
        .nodes
        .iter()
        .map(|&(x, y)| {
            PlateTraceStencil::locate(source, [x, y], budget)
                .map_err(|e| nonlinear(format!("modal interface trace: {e}")))
        })
        .collect()
}

fn segments(mesh: &PlateMesh, edges: &[[usize; 2]]) -> Result<Vec<Segment>, AcousticRealizeError> {
    edges
        .iter()
        .map(|&edge| {
            if edge.iter().any(|&node| node >= mesh.nodes.len()) || edge[0] == edge[1] {
                return Err(invalid("modal interface edge is outside its original mesh"));
            }
            Ok(edge.map(|node| [mesh.nodes[node].0, mesh.nodes[node].1]))
        })
        .collect()
}

fn boundary(mesh: &PlateMesh, budget: PlateLoadBudget) -> Result<Vec<Segment>, AcousticRealizeError> {
    let &(x, y) = mesh
        .nodes
        .first()
        .ok_or_else(|| invalid("empty modal blade mesh"))?;
    // Shared fs-plate admission visits EVERY triangle, including those not used
    // by the first trace. Degenerate/reversed cells never enter an interface.
    PlateTraceStencil::locate(mesh, [x, y], budget)
        .map_err(|e| nonlinear(format!("modal interface geometry: {e}")))?;
    let mut incidence = BTreeMap::<(usize, usize), Vec<[usize; 2]>>::new();
    let mut triangles = BTreeSet::new();
    let mut used = vec![false; mesh.nodes.len()];
    for &tri in &mesh.tris {
        let mut key = tri;
        key.sort_unstable();
        if !triangles.insert(key) {
            return Err(invalid("modal interface contains duplicate triangles"));
        }
        for node in tri {
            used[node] = true;
        }
        for edge in [[tri[0], tri[1]], [tri[1], tri[2]], [tri[2], tri[0]]] {
            incidence
                .entry((edge[0].min(edge[1]), edge[0].max(edge[1])))
                .or_default()
                .push(edge);
        }
    }
    if used.iter().any(|&present| !present) {
        return Err(invalid("modal interface contains unused nodes"));
    }
    for (i, p) in mesh.nodes.iter().enumerate() {
        if mesh.nodes[..i].contains(p) {
            return Err(invalid("modal interface contains coincident independent nodes"));
        }
    }
    let mut next = BTreeMap::new();
    let mut incoming = BTreeSet::new();
    for edges in incidence.values() {
        match edges.as_slice() {
            &[edge] => {
                if next.insert(edge[0], edge[1]).is_some() || !incoming.insert(edge[1]) {
                    return Err(invalid("modal blade boundary branches or touches itself"));
                }
            }
            &[a, b] if a == [b[1], b[0]] => {}
            _ => return Err(invalid("modal blade interface is not an oriented manifold")),
        }
    }
    let (&start, _) = next
        .first_key_value()
        .ok_or_else(|| invalid("modal blade has no boundary"))?;
    let mut node = start;
    let mut edges = Vec::with_capacity(next.len());
    loop {
        let end = *next
            .get(&node)
            .ok_or_else(|| invalid("modal blade boundary is open"))?;
        edges.push([node, end]);
        node = end;
        if node == start {
            break;
        }
        if edges.len() >= next.len() {
            return Err(invalid("modal blade boundary does not close"));
        }
    }
    if edges.len() != next.len() {
        return Err(invalid(
            "independent modal blades require one simple boundary without holes",
        ));
    }
    let boundary = segments(mesh, &edges)?;
    let tolerance = tolerance(&boundary)?;
    for (i, a) in boundary.iter().enumerate() {
        for (j, b) in boundary[..i].iter().enumerate() {
            if j + 1 != i
                && !(j == 0 && i + 1 == boundary.len())
                && intersects(*a, *b, tolerance)
            {
                return Err(invalid("modal blade boundary crosses or touches itself"));
            }
        }
    }
    // Positive cells plus a simple oriented boundary form a single covered
    // planar sheet. Matching its complete boundary rules out holes, hidden
    // overlapping sheets and a mere sample-point approximation to coverage.
    Ok(boundary)
}

fn tolerance(segments: &[Segment]) -> Result<f64, AcousticRealizeError> {
    let first = segments
        .first()
        .ok_or_else(|| invalid("empty modal interface boundary"))?[0];
    let scale = segments.iter().flatten().fold(0.0_f64, |s, p| {
        s.max((p[0] - first[0]).abs())
            .max((p[1] - first[1]).abs())
    });
    if !scale.is_finite() || scale <= 0.0 {
        return Err(invalid("modal interface geometry has no finite physical extent"));
    }
    Ok(128.0 * f64::EPSILON * scale)
}

fn delta(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn intersects(a: Segment, b: Segment, tolerance: f64) -> bool {
    for coordinate in 0..2 {
        if a[0][coordinate].min(a[1][coordinate])
            > b[0][coordinate].max(b[1][coordinate]) + tolerance
            || b[0][coordinate].min(b[1][coordinate])
                > a[0][coordinate].max(a[1][coordinate]) + tolerance
        {
            return false;
        }
    }
    let straddles = |line: Segment, points: Segment| {
        let d = delta(line[1], line[0]);
        let length = d[0].hypot(d[1]);
        let unit = [d[0] / length, d[1] / length];
        let p = cross(unit, delta(points[0], line[0]));
        let q = cross(unit, delta(points[1], line[0]));
        p.min(q) <= tolerance && p.max(q) >= -tolerance
    };
    straddles(a, b) && straddles(b, a)
}

fn same_cover(a: &[Segment], b: &[Segment], oriented: bool) -> Result<(), AcousticRealizeError> {
    let tolerance = tolerance(a)?.max(tolerance(b)?);
    for (reference, candidates) in [(a, b), (b, a)] {
        for &[start, end] in reference {
            let d = delta(end, start);
            let length = d[0].hypot(d[1]);
            if length <= 2.0 * tolerance {
                return Err(invalid("modal interface edge is unresolved at geometry roundoff"));
            }
            let unit = [d[0] / length, d[1] / length];
            let mut intervals = Vec::new();
            for &[p, q] in candidates {
                let (p, q) = (delta(p, start), delta(q, start));
                if cross(unit, p).abs() > tolerance || cross(unit, q).abs() > tolerance {
                    continue;
                }
                let (x, y) = (dot(unit, p), dot(unit, q));
                let (lo, hi) = (x.min(y).max(0.0), x.max(y).min(length));
                if hi - lo > tolerance {
                    if oriented && y <= x {
                        return Err(invalid("modal blade interfaces have opposite orientation"));
                    }
                    intervals.push((lo, hi));
                }
            }
            intervals.sort_by(|a, b| a.0.total_cmp(&b.0));
            let mut covered = 0.0;
            for (lo, hi) in intervals {
                if (lo - covered).abs() > tolerance {
                    return Err(invalid(
                        "modal blades require the same nonoverlapping physical surface and slit",
                    ));
                }
                covered = hi;
            }
            if (covered - length).abs() > tolerance {
                return Err(invalid(
                    "modal blades require complete physical surface and slit coverage",
                ));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn budget() -> PlateLoadBudget {
        PlateLoadBudget {
            max_nodes: 512,
            max_triangles: 2048,
        }
    }
    fn slit(mesh: &PlateMesh, x: f64) -> Vec<[usize; 2]> {
        mesh.boundary_edges()
            .into_iter()
            .filter(|&(a, b)| mesh.nodes[a].0 == x && mesh.nodes[b].0 == x)
            .map(|(a, b)| [a, b])
            .collect()
    }

    #[test]
    fn g0_independent_mesh_trace_preserves_rigid_and_affine_motion() {
        let a = PlateMesh::rectangle(1.0, 1.0, 3, 2);
        let b = PlateMesh::rectangle(1.0, 1.0, 2, 3);
        for (target, source) in [(&a, &b), (&b, &a)] {
            let mapped = traces(
                target, &slit(target, 1.0), source, &slit(source, 1.0),
                budget(), 100_000,
            ).unwrap();
            for (&(x, y), trace) in target.nodes.iter().zip(mapped) {
                let actual: f64 = trace.nodes().into_iter().zip(trace.weights())
                    .map(|(node, weight)| {
                        let (u, v) = source.nodes[node];
                        weight * (0.7 + 0.2 * u - 0.3 * v)
                    }).sum();
                assert!((actual - (0.7 + 0.2 * x - 0.3 * y)).abs() < 2e-14);
                assert!((trace.weights().iter().sum::<f64>() - 1.0).abs() < 2e-14);
            }
        }
    }

    #[test]
    fn g0_independent_interfaces_refuse_ambiguous_or_noncovering_geometry() {
        let a = PlateMesh::rectangle(1.0, 1.0, 3, 2);
        let source = || PlateMesh::rectangle(1.0, 1.0, 2, 3);
        let transfer = |b: &PlateMesh, edges: &[[usize; 2]], max_terms| {
            traces(&a, &slit(&a, 1.0), b, edges, budget(), max_terms)
        };
        let b = source();
        assert!(transfer(&b, &slit(&b, 1.0), 100_000).is_ok());
        assert!(transfer(&b, &slit(&b, 0.0), 100_000).is_err());
        assert!(transfer(&b, &slit(&b, 1.0), 1).is_err());

        let mut shifted = source();
        for point in &mut shifted.nodes { point.0 += 0.001; }
        let mut duplicate = source();
        duplicate.tris.push(duplicate.tris[0]);
        let mut missing = source();
        missing.tris.remove(2);
        let mut degenerate = source();
        degenerate.nodes[1] = degenerate.nodes[0];
        for bad in [shifted, duplicate, missing, degenerate] {
            assert!(transfer(&bad, &[[2, 5]], 100_000).is_err());
        }

        // All fan cells are positively oriented, but their star boundary
        // crosses itself. Positive element areas alone do not admit a sheet.
        let mut crossed = source();
        crossed.nodes = vec![(0.0, 0.0)];
        for vertex in [0, 2, 4, 1, 3] {
            let angle = std::f64::consts::TAU * f64::from(vertex) / 5.0;
            crossed.nodes.push((angle.cos(), angle.sin()));
        }
        crossed.tris = (0..5).map(|i| [0, i + 1, (i + 1) % 5 + 1]).collect();
        assert!(transfer(&crossed, &[[1, 2]], 100_000).is_err());
    }
}
