//! Current-area pressure on the original oriented P1 membrane triangles.
//!
//! A uniform signed pressure jump acts along the current reference-oriented
//! normal: each triangle contributes `p (x_b-x_a) cross (x_c-x_a) / 6` to
//! each of its three nodes. Differentiating the area vector gives the exact
//! follower-load tangent; freezing the reference normal omits this term.

use super::{MembraneError, MembranePatch, finite_vectors, poll};
use std::collections::BTreeMap;

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

impl MembranePatch {
    fn admit_pressure_area(&self, element: usize, area: [f64; 3]) -> Result<(), MembraneError> {
        let ratio = area[0].hypot(area[1]).hypot(area[2]) / (2.0 * self.elements[element].area);
        if !ratio.is_finite() || ratio < self.domain.minimum_area_ratio {
            return Err(MembraneError::State {
                element,
                reason: "pressure surface collapse or unrepresentable area",
            });
        }
        Ok(())
    }

    fn pressure_edges(&self, nodes: [usize; 3], u: &[[f64; 3]]) -> [[f64; 3]; 2] {
        std::array::from_fn(|edge| {
            std::array::from_fn(|axis| {
                let a = nodes[0];
                let b = nodes[edge + 1];
                let reference = if axis < 2 {
                    self.mesh.nodes[b][axis] - self.mesh.nodes[a][axis]
                } else {
                    0.0
                };
                reference + (u[b][axis] - u[a][axis])
            })
        })
    }

    /// Integrate a uniform signed pressure jump [Pa] on CURRENT triangle areas.
    /// Positive pressure acts along the oriented normal (initially +z), and
    /// negative pressure acts oppositely. This is a pressure difference, not
    /// absolute thermodynamic pressure. Loads at prescribed nodes are retained.
    /// The original geometric envelope applies; material response is separate.
    pub fn pressure_load_n(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, MembraneError> {
        poll(&mut cancelled)?;
        let n = self.mesh.nodes.len();
        if !pressure_pa.is_finite() || !finite_vectors(displacement_m, n) {
            return Err(MembraneError::Invalid(
                "finite pressure jump and per-node displacement required",
            ));
        }
        let mut force = vec![[0.0; 3]; n];
        for (index, e) in self.elements.iter().enumerate() {
            poll(&mut cancelled)?;
            let [a, b] = self.pressure_edges(e.nodes, displacement_m);
            let area = cross(a, b);
            self.admit_pressure_area(index, area)?;
            for &node in &e.nodes {
                for axis in 0..3 {
                    force[node][axis] += (pressure_pa / 6.0) * area[axis];
                }
            }
        }
        if !finite_vectors(&force, n) {
            return Err(MembraneError::Invalid("unrepresentable pressure load"));
        }
        poll(&mut cancelled)?;
        Ok(force)
    }

    /// Exact matrix-free derivative of `pressure_load_n` with respect to the
    /// displacement, in the supplied direction [m]. This is the EXTERNAL force
    /// tangent; force balance subtracts it from the material/geometric tangent.
    /// A free-edge pressure load need not have a symmetric tangent.
    pub fn pressure_tangent_apply(
        &self,
        displacement_m: &[[f64; 3]],
        direction_m: &[[f64; 3]],
        pressure_pa: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, MembraneError> {
        poll(&mut cancelled)?;
        let n = self.mesh.nodes.len();
        if !pressure_pa.is_finite()
            || !finite_vectors(displacement_m, n)
            || !finite_vectors(direction_m, n)
        {
            return Err(MembraneError::Invalid(
                "finite pressure jump, displacement and direction required",
            ));
        }
        let mut force = vec![[0.0; 3]; n];
        for (index, e) in self.elements.iter().enumerate() {
            poll(&mut cancelled)?;
            let [a, b] = self.pressure_edges(e.nodes, displacement_m);
            self.admit_pressure_area(index, cross(a, b))?;
            let da =
                std::array::from_fn(|i| direction_m[e.nodes[1]][i] - direction_m[e.nodes[0]][i]);
            let db =
                std::array::from_fn(|i| direction_m[e.nodes[2]][i] - direction_m[e.nodes[0]][i]);
            let left = cross(da, b);
            let right = cross(a, db);
            for &node in &e.nodes {
                for axis in 0..3 {
                    force[node][axis] += (pressure_pa / 6.0) * (left[axis] + right[axis]);
                }
            }
        }
        if !finite_vectors(&force, n) {
            return Err(MembraneError::Invalid("unrepresentable pressure tangent"));
        }
        poll(&mut cancelled)?;
        Ok(force)
    }

    /// Exact pressure work [J] along a STRAIGHT nodal-displacement segment at
    /// held pressure. The area vector is quadratic along the segment, so its
    /// mean is integrated analytically rather than approximated at a midpoint.
    /// With the perimeter fixed this is `p * swept_volume`, a conservative
    /// potential difference on the free coordinates. With moving edges it is
    /// path-dependent; no closed cavity or fluid inventory is inferred.
    pub fn pressure_work_j(
        &self,
        before_m: &[[f64; 3]],
        after_m: &[[f64; 3]],
        pressure_pa: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, MembraneError> {
        poll(&mut cancelled)?;
        let n = self.mesh.nodes.len();
        if !pressure_pa.is_finite() || !finite_vectors(before_m, n) || !finite_vectors(after_m, n) {
            return Err(MembraneError::Invalid(
                "finite pressure jump and segment endpoints required",
            ));
        }
        let mut work = 0.0;
        for e in &self.elements {
            poll(&mut cancelled)?;
            let [a, b] = self.pressure_edges(e.nodes, before_m);
            let change: [[f64; 3]; 3] = e
                .nodes
                .map(|node| std::array::from_fn(|i| after_m[node][i] - before_m[node][i]));
            let da = std::array::from_fn(|i| change[1][i] - change[0][i]);
            let db = std::array::from_fn(|i| change[2][i] - change[0][i]);
            let c0 = cross(a, b);
            let c1a = cross(da, b);
            let c1b = cross(a, db);
            let c2 = cross(da, db);
            for i in 0..3 {
                let mean_area = c0[i] + 0.5 * (c1a[i] + c1b[i]) + c2[i] / 3.0;
                work +=
                    (pressure_pa / 6.0) * mean_area * (change[0][i] + change[1][i] + change[2][i]);
            }
        }
        if !work.is_finite() {
            return Err(MembraneError::Invalid("unrepresentable pressure work"));
        }
        poll(&mut cancelled)?;
        Ok(work)
    }

    /// Admit the fixed perimeter required for a conservative pressure-volume
    /// coupling. `fixed` is a flat, three-coordinate-per-node displacement map:
    /// `Some(value)` prescribes that coordinate [m], and `None` leaves it free.
    /// Every open edge must have all coordinates prescribed at both endpoints;
    /// interior edges must have consistent orientation and at most two incident
    /// triangles. This checks actual triangle incidence, not named patches.
    ///
    /// Callers must preserve these prescribed values in every configuration
    /// used with [`Self::pressure_work_j`] to obtain a path-independent volume
    /// change. This admission does not establish an absolute enclosed volume.
    pub fn admit_pressure_perimeter(
        &self,
        fixed: &[Option<f64>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), MembraneError> {
        poll(&mut cancelled)?;
        if fixed.len() != 3 * self.mesh.nodes.len()
            || fixed.iter().flatten().any(|value| !value.is_finite())
        {
            return Err(MembraneError::Invalid(
                "finite prescribed displacement map with three coordinates per node required",
            ));
        }
        let mut edges = BTreeMap::<[usize; 2], (u8, i8)>::new();
        for e in &self.elements {
            poll(&mut cancelled)?;
            for (a, b) in [
                (e.nodes[0], e.nodes[1]),
                (e.nodes[1], e.nodes[2]),
                (e.nodes[2], e.nodes[0]),
            ] {
                let (key, sign) = if a < b { ([a, b], 1) } else { ([b, a], -1) };
                let entry = edges.entry(key).or_insert((0, 0));
                if entry.0 == 2 {
                    return Err(MembraneError::Invalid("nonmanifold pressure membrane edge"));
                }
                entry.0 += 1;
                entry.1 += sign;
            }
        }
        for (edge, (count, sign)) in edges {
            poll(&mut cancelled)?;
            if count == 2 && sign != 0 {
                return Err(MembraneError::Invalid(
                    "inconsistent pressure membrane orientation",
                ));
            }
            if count == 1
                && edge
                    .iter()
                    .any(|&node| (0..3).any(|i| fixed[3 * node + i].is_none()))
            {
                return Err(MembraneError::Invalid(
                    "uniform-pressure equilibrium requires a fully prescribed perimeter",
                ));
            }
        }
        poll(&mut cancelled)
    }
}
