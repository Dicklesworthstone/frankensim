//! Current-area pressure on the original oriented P1 membrane triangles.
//!
//! Signed pressure jumps act along the current reference-oriented normal.
//! The P1 surface mass integrates material pressure samples against the nodal
//! displacement basis. Differentiating the area vector gives the exact follower
//! tangent; integrating its quadratic drift gives the exact held-profile work.

use super::{MembraneError, MembranePatch, finite_vectors, poll};
use std::collections::BTreeMap;

/// Signed pressure difference [Pa] acting on the original membrane triangles.
///
/// Positive values act along the current oriented normal, initially +z.
/// `TriangleLinear` holds three P1 samples per original triangle, in that
/// triangle's node order. Adjacent triangles may have different traces. Samples
/// are held in MATERIAL/barycentric coordinates while geometry changes: an
/// Eulerian hydrostatic or fluid field must be resampled by its source owner.
///
/// A nonuniform profile need not be conservative even with a fixed perimeter.
/// It supplies no absolute gas state, pressure derivative, fluid inertia or
/// implicit interpolation from another mesh.
#[derive(Clone, Copy, Debug)]
pub enum MembranePressure<'a> {
    Uniform(f64),
    TriangleLinear(&'a [[f64; 3]]),
}

impl From<f64> for MembranePressure<'_> {
    fn from(pressure_pa: f64) -> Self {
        Self::Uniform(pressure_pa)
    }
}

impl MembranePressure<'_> {
    pub(super) fn validate(
        self,
        triangles: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MembraneError> {
        match self {
            Self::Uniform(pressure) if pressure.is_finite() => Ok(()),
            Self::TriangleLinear(values) if values.len() == triangles => {
                for samples in values {
                    poll(cancelled)?;
                    if samples.iter().any(|p| !p.is_finite()) {
                        return Err(MembraneError::Invalid(
                            "finite P1 pressure samples required",
                        ));
                    }
                }
                Ok(())
            }
            _ => Err(MembraneError::Invalid(
                "finite uniform pressure or three finite samples per original triangle required",
            )),
        }
    }

    pub(super) fn is_zero(self) -> bool {
        matches!(self, Self::Uniform(pressure) if pressure == 0.0)
    }

    fn samples(self, element: usize) -> [f64; 3] {
        match self {
            Self::Uniform(pressure) => [pressure; 3],
            Self::TriangleLinear(values) => values[element],
        }
    }

    // Integral Ni*p dA = cross(a,b) * (2*pi + pj + pk)/24.
    // Divide before adding to avoid overflow in otherwise finite averages.
    fn weights(self, element: usize) -> [f64; 3] {
        let p = self.samples(element);
        if p[0] == p[1] && p[0] == p[2] {
            [p[0] / 6.0; 3]
        } else {
            std::array::from_fn(|i| p[i] / 12.0 + p[(i + 1) % 3] / 24.0 + p[(i + 2) % 3] / 24.0)
        }
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn mean_area(a: [f64; 3], b: [f64; 3], da: [f64; 3], db: [f64; 3]) -> [f64; 3] {
    let c0 = cross(a, b);
    let c1a = cross(da, b);
    let c1b = cross(a, db);
    let c2 = cross(da, db);
    std::array::from_fn(|i| c0[i] + 0.5 * (c1a[i] + c1b[i]) + c2[i] / 3.0)
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

    /// Integrate uniform or material P1 pressure [Pa] on CURRENT triangle areas.
    /// Positive pressure acts along the oriented normal (initially +z), and
    /// negative pressure acts oppositely. This is a pressure difference, not
    /// absolute thermodynamic pressure. Loads at prescribed nodes are retained.
    /// The original geometric envelope applies; material response is separate.
    pub fn pressure_load_n<'p>(
        &self,
        displacement_m: &[[f64; 3]],
        pressure: impl Into<MembranePressure<'p>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, MembraneError> {
        let pressure = pressure.into();
        pressure.validate(self.elements.len(), &mut cancelled)?;
        let mut force = vec![[0.0; 3]; self.mesh.nodes.len()];
        self.pressure_load_into_n(displacement_m, pressure, &mut force, &mut cancelled)?;
        Ok(force)
    }

    /// Shared current-area assembly into the time owner's retained workspace.
    /// The immutable pressure profile has already passed shape/value admission.
    pub(super) fn pressure_load_into_n(
        &self,
        displacement_m: &[[f64; 3]],
        pressure: MembranePressure<'_>,
        force: &mut [[f64; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), MembraneError> {
        poll(cancelled)?;
        let n = self.mesh.nodes.len();
        if !finite_vectors(displacement_m, n) || force.len() != n {
            return Err(MembraneError::Invalid(
                "finite pressure jump, per-node displacement and complete force workspace required",
            ));
        }
        force.fill([0.0; 3]);
        for (index, e) in self.elements.iter().enumerate() {
            poll(cancelled)?;
            let [a, b] = self.pressure_edges(e.nodes, displacement_m);
            let area = cross(a, b);
            self.admit_pressure_area(index, area)?;
            let weights = pressure.weights(index);
            for (local, &node) in e.nodes.iter().enumerate() {
                for axis in 0..3 {
                    force[node][axis] += weights[local] * area[axis];
                }
            }
        }
        if !finite_vectors(force, n) {
            return Err(MembraneError::Invalid("unrepresentable pressure load"));
        }
        poll(cancelled)?;
        Ok(())
    }

    /// Exact matrix-free derivative of `pressure_load_n` with respect to the
    /// displacement, in the supplied direction [m]. This is the EXTERNAL force
    /// tangent; force balance subtracts it from the material/geometric tangent.
    /// Material samples are held fixed. Nonuniform pressure, or uniform pressure
    /// on a free perimeter, need not have a symmetric tangent.
    pub fn pressure_tangent_apply<'p>(
        &self,
        displacement_m: &[[f64; 3]],
        direction_m: &[[f64; 3]],
        pressure: impl Into<MembranePressure<'p>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, MembraneError> {
        poll(&mut cancelled)?;
        let n = self.mesh.nodes.len();
        let pressure = pressure.into();
        pressure.validate(self.elements.len(), &mut cancelled)?;
        if !finite_vectors(displacement_m, n) || !finite_vectors(direction_m, n) {
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
            let weights = pressure.weights(index);
            for (local, &node) in e.nodes.iter().enumerate() {
                for axis in 0..3 {
                    force[node][axis] += weights[local] * (left[axis] + right[axis]);
                }
            }
        }
        if !finite_vectors(&force, n) {
            return Err(MembraneError::Invalid("unrepresentable pressure tangent"));
        }
        poll(&mut cancelled)?;
        Ok(force)
    }

    /// Mean current-area pressure force along a STRAIGHT displacement segment.
    ///
    /// Pressure samples stay fixed in material coordinates. The current area
    /// vector is quadratic on this path, so this mean is integrated exactly.
    /// Its nodal dot product with `after_m-before_m` equals the original
    /// `pressure_work_j` up to floating-point summation. For uniform pressure
    /// and a fixed perimeter, unit pressure supplies a discrete volume gradient.
    /// This is a path average, not the instantaneous force at the midpoint.
    ///
    /// Both endpoint areas are admitted. This operator alone does not certify
    /// intermediate material states or make a nonuniform profile conservative.
    pub fn pressure_segment_load_n<'p>(
        &self,
        before_m: &[[f64; 3]],
        after_m: &[[f64; 3]],
        pressure: impl Into<MembranePressure<'p>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, MembraneError> {
        self.pressure_segment_action(before_m, after_m, None, pressure.into(), &mut cancelled)
    }

    /// Exact directional derivative of `pressure_segment_load_n` with respect
    /// to its AFTER endpoint, holding the BEFORE endpoint and pressure fixed.
    /// The direction is an endpoint displacement [m], without a hidden time or
    /// midpoint factor. A midpoint coordinate with `q1=2*q_mid-q0` therefore
    /// contributes twice this action by the chain rule.
    pub fn pressure_segment_tangent_apply<'p>(
        &self,
        before_m: &[[f64; 3]],
        after_m: &[[f64; 3]],
        direction_m: &[[f64; 3]],
        pressure: impl Into<MembranePressure<'p>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, MembraneError> {
        self.pressure_segment_action(
            before_m,
            after_m,
            Some(direction_m),
            pressure.into(),
            &mut cancelled,
        )
    }

    fn pressure_segment_action(
        &self,
        before_m: &[[f64; 3]],
        after_m: &[[f64; 3]],
        direction_m: Option<&[[f64; 3]]>,
        pressure: MembranePressure<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, MembraneError> {
        poll(cancelled)?;
        let n = self.mesh.nodes.len();
        pressure.validate(self.elements.len(), cancelled)?;
        if !finite_vectors(before_m, n)
            || !finite_vectors(after_m, n)
            || direction_m.is_some_and(|direction| !finite_vectors(direction, n))
        {
            return Err(MembraneError::Invalid(
                "finite complete pressure-segment endpoints and tangent direction required",
            ));
        }
        let mut force = vec![[0.0; 3]; n];
        for (index, e) in self.elements.iter().enumerate() {
            poll(cancelled)?;
            let [a, b] = self.pressure_edges(e.nodes, before_m);
            let [end_a, end_b] = self.pressure_edges(e.nodes, after_m);
            self.admit_pressure_area(index, cross(a, b))?;
            self.admit_pressure_area(index, cross(end_a, end_b))?;
            let change = e.nodes.map(|node| {
                std::array::from_fn::<_, 3, _>(|i| after_m[node][i] - before_m[node][i])
            });
            let da = std::array::from_fn(|i| change[1][i] - change[0][i]);
            let db = std::array::from_fn(|i| change[2][i] - change[0][i]);
            let area = if let Some(direction) = direction_m {
                let ha =
                    std::array::from_fn(|i| direction[e.nodes[1]][i] - direction[e.nodes[0]][i]);
                let hb =
                    std::array::from_fn(|i| direction[e.nodes[2]][i] - direction[e.nodes[0]][i]);
                let first = cross(ha, b);
                let second = cross(a, hb);
                let quadratic_a = cross(ha, db);
                let quadratic_b = cross(da, hb);
                std::array::from_fn(|i| {
                    0.5 * (first[i] + second[i]) + (quadratic_a[i] + quadratic_b[i]) / 3.0
                })
            } else {
                mean_area(a, b, da, db)
            };
            let weights = pressure.weights(index);
            for (local, &node) in e.nodes.iter().enumerate() {
                for axis in 0..3 {
                    force[node][axis] += weights[local] * area[axis];
                }
            }
        }
        if !finite_vectors(&force, n) {
            return Err(MembraneError::Invalid(
                "unrepresentable pressure-segment action",
            ));
        }
        poll(cancelled)?;
        Ok(force)
    }

    /// Exact pressure work [J] along a STRAIGHT nodal-displacement segment at
    /// held material pressure. The area vector is quadratic along the segment,
    /// so its mean is integrated analytically. Uniform pressure with a fixed
    /// perimeter gives the conservative `p * swept_volume`; nonuniform samples
    /// or moving edges generally give path-dependent work. No closed cavity or
    /// fluid inventory is inferred.
    pub fn pressure_work_j<'p>(
        &self,
        before_m: &[[f64; 3]],
        after_m: &[[f64; 3]],
        pressure: impl Into<MembranePressure<'p>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, MembraneError> {
        poll(&mut cancelled)?;
        let n = self.mesh.nodes.len();
        let pressure = pressure.into();
        pressure.validate(self.elements.len(), &mut cancelled)?;
        if !finite_vectors(before_m, n) || !finite_vectors(after_m, n) {
            return Err(MembraneError::Invalid(
                "finite pressure jump and segment endpoints required",
            ));
        }
        self.pressure_segment_work_j(
            before_m,
            pressure,
            |node| std::array::from_fn(|i| after_m[node][i] - before_m[node][i]),
            &mut cancelled,
        )
    }

    /// Work over the original resolved drift, even if adding that increment to
    /// a translated absolute displacement would round it away. The pressure
    /// owner still integrates the same changing triangle area analytically.
    /// The immutable pressure profile is admitted once by the time owner.
    pub(super) fn pressure_work_from_increment_j(
        &self,
        before_m: &[[f64; 3]],
        increment_m: &[[f64; 3]],
        pressure: MembranePressure<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, MembraneError> {
        poll(cancelled)?;
        let n = self.mesh.nodes.len();
        if !finite_vectors(before_m, n) || !finite_vectors(increment_m, n) {
            return Err(MembraneError::Invalid(
                "finite pressure jump, initial displacement and drift required",
            ));
        }
        self.pressure_segment_work_j(before_m, pressure, |node| increment_m[node], cancelled)
    }

    fn pressure_segment_work_j(
        &self,
        before_m: &[[f64; 3]],
        pressure: MembranePressure<'_>,
        change_at: impl Fn(usize) -> [f64; 3],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, MembraneError> {
        let mut work = 0.0;
        for (index, e) in self.elements.iter().enumerate() {
            poll(cancelled)?;
            let [a, b] = self.pressure_edges(e.nodes, before_m);
            let change = e.nodes.map(&change_at);
            let da = std::array::from_fn(|i| change[1][i] - change[0][i]);
            let db = std::array::from_fn(|i| change[2][i] - change[0][i]);
            let area = mean_area(a, b, da, db);
            let samples = pressure.samples(index);
            let weights = pressure.weights(index);
            for i in 0..3 {
                let mean_area = area[i];
                if samples[0] == samples[1] && samples[0] == samples[2] {
                    // Retain the original uniform-pressure arithmetic exactly.
                    work += (samples[0] / 6.0)
                        * mean_area
                        * (change[0][i] + change[1][i] + change[2][i]);
                } else {
                    for local in 0..3 {
                        work += weights[local] * mean_area * change[local][i];
                    }
                }
            }
        }
        if !work.is_finite() {
            return Err(MembraneError::Invalid("unrepresentable pressure work"));
        }
        poll(cancelled)?;
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
    /// used with UNIFORM pressure in [`Self::pressure_work_j`] to obtain a
    /// path-independent volume change. Nonuniform samples retain path-dependent
    /// work. This admission does not establish an absolute enclosed volume.
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
