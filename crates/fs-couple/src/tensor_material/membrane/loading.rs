//! Consistent mechanical loads on the SAME immutable P1 surface and supports.
use super::*;
use std::collections::BTreeMap;

/// Vector line traction on one actual boundary edge, linear between endpoints.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MembraneEdgeTraction {
    /// Original vertex indices; reversing them requires reversing the values.
    pub nodes: [usize; 2],
    /// Cartesian force per reference edge length [N/m] at the two endpoints.
    /// Already integrated through thickness; thickness is not applied twice.
    pub values_n_m: [[f64; 3]; 2],
}

impl TensorMembrane {
    /// Integrate surface and boundary tractions and uniform body acceleration.
    /// `traction_pa[e][vertex][axis]` is a per-facet P1 Cartesian surface load
    /// [N/m^2], permitting discontinuities between adjacent material facets.
    /// Explicit zeros describe unloaded faces. Acceleration is in m/s^2.
    ///
    /// Exact consistent P1 weights are A/12*(2*t_i+t_j+t_k) on triangles and
    /// L/6*(2*t_i+t_j) on edges. Consequently nodal work equals the physical
    /// traction work for any represented P1 displacement/velocity field.
    /// The unprojected output includes loads on supports. No force is moved
    /// from a constrained node to a free node. Loads use reference geometry,
    /// not a nonlinear follower-pressure law or a computed support reaction.
    ///
    /// Every input and assembled value must be finite. Interior/nonexistent
    /// edges, duplicate edges, and cancellation return no partial load vector.
    pub fn distributed_nodal_forces(&self, traction_pa: &[[[f64;3];3]],
        edges: &[MembraneEdgeTraction], acceleration_m_s2: [f64;3], gate: &CancelGate)
        -> Result<Vec<[f64;3]>>
    {
        poll(gate)?;
        if traction_pa.len() != self.mesh.tris.len()
            || traction_pa.iter().flatten().flatten().any(|v| !v.is_finite())
            || acceleration_m_s2.iter().any(|v| !v.is_finite())
            || edges.len() > 3*self.mesh.tris.len()
            || edges.iter().any(|e| e.values_n_m.iter().flatten().any(|v| !v.is_finite())) {
            return Err(bad("complete finite surface tractions, bounded edge loads and finite acceleration required"));
        }
        let mut incidence = BTreeMap::new();
        for tri in &self.mesh.tris {
            poll(gate)?;
            for mut edge in [[tri[0],tri[1]],[tri[1],tri[2]],[tri[2],tri[0]]] {
                edge.sort_unstable(); *incidence.entry(edge).or_insert(0_usize) += 1;
            }
        }
        let mut loaded = BTreeSet::new();
        for edge in edges {
            poll(gate)?;
            let mut key = edge.nodes; key.sort_unstable();
            if incidence.get(&key) != Some(&1) { return Err(bad("line traction must name an actual boundary edge")); }
            if !loaded.insert(key) { return Err(bad("duplicate boundary-edge traction")); }
        }
        let mut force: Vec<[f64;3]> = self.nodal_mass.iter().map(|mass|
            acceleration_m_s2.map(|a| mass*a)).collect();
        for ((tri, facet), traction) in self.mesh.tris.iter().zip(&self.facets).zip(traction_pa) {
            poll(gate)?;
            for i in 0..3 {for axis in 0..3 {
                force[tri[i]][axis] += facet.area_m2/12.0*
                    (2.0*traction[i][axis]+traction[(i+1)%3][axis]+traction[(i+2)%3][axis]);
            }}
        }
        for edge in edges {
            poll(gate)?;
            let [a,b] = edge.nodes;
            let delta: [f64;3] = core::array::from_fn(|i| self.mesh.nodes[b][i]-self.mesh.nodes[a][i]);
            let length = fs_math::det::sqrt(delta.iter().map(|v|v*v).sum());
            if !length.is_finite() || length <= 0.0 { return Err(bad("unresolved loaded edge length")); }
            for i in 0..2 {for axis in 0..3 {
                force[edge.nodes[i]][axis] += length/6.0*
                    (2.0*edge.values_n_m[i][axis]+edge.values_n_m[1-i][axis]);
            }}
        }
        self.check_nodal(&force)?; poll(gate)?;
        Ok(force)
    }

    /// Apply a signed pressure jump at each triangle's three vertices [Pa].
    /// Positive pressure acts AGAINST its original oriented ez normal; negative
    /// pressure reverses that load. Uses the same consistent P1 work map.
    /// This supplies a load, not bending/prestress or an acoustic field solver.
    pub fn pressure_nodal_forces(&self, pressure_pa: &[[f64;3]], gate: &CancelGate)
        -> Result<Vec<[f64;3]>>
    {
        poll(gate)?;
        if pressure_pa.len() != self.mesh.tris.len()
            || pressure_pa.iter().flatten().any(|v| !v.is_finite()) {
            return Err(bad("complete finite per-facet nodal pressure field required"));
        }
        let mut traction = Vec::with_capacity(pressure_pa.len());
        for (p, facet) in pressure_pa.iter().zip(&self.facets) {
            poll(gate)?;
            traction.push(p.map(|p| facet.frame[2].map(|n| -p*n)));
        }
        self.distributed_nodal_forces(&traction, &[], [0.0;3], gate)
    }

    /// Advance THIS mesh's model-bound history under physical nodal forces.
    /// The existing exact integrator owns work, loss, state and publication.
    /// A foreign history, invalid force, clock or cancellation leaves it intact.
    pub fn advance_nodal(&self, history: &mut TensorMaterialHistory<'_>, force_n: &[[f64;3]],
        dt_s: f64, gate: &CancelGate) -> Result<fs_phs::StepRecord>
    {
        poll(gate)?;
        if !std::ptr::eq(history.body, &self.body) { return Err(bad("history belongs to a different membrane body")); }
        let force = self.project_nodal_forces(force_n)?;
        history.advance(&force, dt_s, gate).map_err(MembraneError::Dynamics)
    }
}

#[cfg(test)]
mod tests;
