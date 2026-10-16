//! Nested h-refinement with a fixed polygonal geometry axis.
use super::*;
impl PlanarMesh {
    /// Split every triangle into four similar children. Input boundary chains
    /// are retained exactly; unrepresentable boundary midpoints refuse.
    /// This preserves element shapes, not a new Ruppert quality guarantee:
    /// descendants of acute-corner exempt elements need a fresh quality audit.
    pub fn refine_uniform(&self, budget: PlanarBudget, cx: &Cx<'_>) -> Result<Self, PlanarError> {
        checkpoint(cx)?;
        let edges = triangulation::edges(&self.triangles);
        let boundary: BTreeSet<_> = self
            .receipt
            .boundary_loops
            .iter()
            .flat_map(|ring| (0..ring.len()).map(|i| edge(ring[i], ring[(i + 1) % ring.len()])))
            .collect();
        let nodes = self
            .vertices
            .len()
            .checked_add(edges.len())
            .ok_or(PlanarError::InvalidInput("refinement vertex overflow"))?;
        let faces = self
            .triangles
            .len()
            .checked_mul(4)
            .ok_or(PlanarError::InvalidInput("refinement triangle overflow"))?;
        cap("vertices", nodes, budget.max_vertices)?;
        cap("triangles", faces, budget.max_triangles)?;
        cap("refinement steps", edges.len(), budget.max_refinement_steps)?;
        let mut vertices = self.vertices.clone();
        let mut midpoint = BTreeMap::new();
        for &(a, b) in edges.keys() {
            checkpoint(cx)?;
            let p = [
                f64::midpoint(vertices[a][0], vertices[b][0]),
                f64::midpoint(vertices[a][1], vertices[b][1]),
            ];
            if p == vertices[a]
                || p == vertices[b]
                || (boundary.contains(&(a, b)) && !on_segment(vertices[a], vertices[b], p))
            {
                return Err(PlanarError::Precision(
                    "uniform refinement needs a representable collinear midpoint",
                ));
            }
            midpoint.insert((a, b), vertices.len());
            vertices.push(p);
        }
        let mut triangles = Vec::with_capacity(faces);
        for &[a, b, c] in &self.triangles {
            checkpoint(cx)?;
            let ab = midpoint[&edge(a, b)];
            let bc = midpoint[&edge(b, c)];
            let ca = midpoint[&edge(c, a)];
            triangles.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]);
        }
        if triangles
            .iter()
            .any(|t| orient2d(vertices[t[0]], vertices[t[1]], vertices[t[2]]) != Sign::Positive)
        {
            return Err(PlanarError::Precision(
                "uniform refinement inverted a represented child",
            ));
        }
        let chains = self
            .receipt
            .source_segment_chains
            .iter()
            .map(|chain| {
                let mut out = Vec::with_capacity(2 * chain.len() - 1);
                for pair in chain.windows(2) {
                    out.extend([pair[0], midpoint[&edge(pair[0], pair[1])]]);
                }
                out.push(*chain.last().unwrap());
                out
            })
            .collect();
        let loops = self
            .receipt
            .boundary_loops
            .iter()
            .map(|ring| {
                let mut out = Vec::with_capacity(2 * ring.len());
                for i in 0..ring.len() {
                    out.extend([
                        ring[i],
                        midpoint[&edge(ring[i], ring[(i + 1) % ring.len()])],
                    ]);
                }
                out
            })
            .collect();
        let mut receipt = self.receipt.clone();
        receipt.source_segment_chains = chains;
        receipt.boundary_loops = loops;
        receipt.exempt_triangles = triangles
            .iter()
            .enumerate()
            .filter(|(_, t)| t.iter().any(|i| receipt.acute_input_vertices.contains(i)))
            .map(|(i, _)| i)
            .collect();
        receipt.minimum_angle_deg = triangles
            .iter()
            .map(|t| minimum_angle(t.map(|i| vertices[i])))
            .fold(f64::INFINITY, f64::min);
        receipt.mesh_area_m2 = triangles
            .iter()
            .map(|t| area(&t.map(|i| vertices[i])))
            .sum();
        receipt.refinement_insertions += edges.len();
        if (receipt.mesh_area_m2 / receipt.polygon_area_m2 - 1.0).abs() > 1e-12 {
            return Err(PlanarError::Topology(
                "uniform refinement failed area conservation",
            ));
        }
        Ok(Self {
            vertices,
            triangles,
            receipt,
        })
    }
}
