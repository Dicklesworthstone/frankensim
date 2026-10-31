//! Follower pressure on the oriented exterior faces of a tet mesh.
//!
//! Positive pressure acts inward on CURRENT triangle areas. Uniform pressure
//! on a closed boundary has force -p grad(V), tangent -p Hessian(V), and
//! mechanical work over a held-pressure interval is -p (V_after - V_before).
//! The last expression is the endpoint volume difference, not midpoint force
//! dotted with a displacement increment. P1 tetrahedral volume is cubic in
//! position, so those two work expressions need not agree.
//!
//! Preparation reuses the native load-geometry incidence/quality admission,
//! then requires exactly two oppositely oriented faces at every boundary edge.
//! A per-face pressure field can load selected patches (zero on other faces)
//! or separate pressure regions. Its tangent is generally nonsymmetric and
//! its work is path-dependent. Held pressures have exact work along a straight
//! nodal-displacement segment, because current face area is quadratic along
//! that segment. This does not assert that a dynamic trajectory is straight.
//! There is no global injectivity certificate, fluid equation of state, or
//! time integrator in this module.

use std::collections::BTreeMap;

use fs_exec::Cx;
use fs_sparse::{Coo, Csr};

use super::{HyperTetError, checkpoint};
use crate::linear3::TetAssemblyBudget;
use crate::loads3::{LoadError3, TetLoadGeometry3};

/// Enclosed void volume and work-conjugate pressure geometry for gas cavities.
pub mod cavity;

/// Signed differential pressure in Pa, held in material-face coordinates.
/// Positive pressure acts inward; negative pressure prescribes outward traction.
#[derive(Debug, Clone, Copy)]
pub enum TetPressure<'a> {
    /// One pressure on the complete closed boundary, with potential `p V`.
    Uniform(f64),
    /// One pressure per canonical [`ClosedTetPressureSurface::faces`] entry.
    /// Zero leaves a face unloaded. No conservative potential is assumed.
    PerFace(&'a [f64]),
}

impl From<f64> for TetPressure<'_> {
    fn from(pressure_pa: f64) -> Self {
        Self::Uniform(pressure_pa)
    }
}

impl TetPressure<'_> {
    fn at(self, face: usize) -> f64 {
        match self {
            Self::Uniform(value) => value,
            Self::PerFace(values) => values[face],
        }
    }
}

/// Prepared complete exterior of one immutable reference tetrahedral mesh.
///
/// All indices and reference coordinates are owned. A consumer must bind this
/// surface to the same reference mesh, rather than just a matching node count.
/// Either reference tet winding is admitted; current tets must retain it.
#[derive(Debug, Clone)]
pub struct ClosedTetPressureSurface {
    nodes_m: Vec<[f64; 3]>,
    tetrahedra: Vec<[usize; 4]>,
    orientations: Vec<f64>,
    faces: Vec<[usize; 3]>,
}

impl ClosedTetPressureSurface {
    /// Derive all exterior faces and admit their closed, outward orientation.
    ///
    /// The face budget bounds the native incidence map (four entries per tet),
    /// the boundary-edge map (at most three entries per face), and a pressure
    /// tangent (at most 81 scalar contributions per exterior face). All such
    /// products are checked before allocation. Cancellation is polled per
    /// node, tetrahedron, face and edge operation.
    pub fn new(
        nodes_m: &[[f64; 3]],
        tetrahedra: &[[usize; 4]],
        budget: TetAssemblyBudget,
        maximum_face_entries: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, HyperTetError> {
        checkpoint(cx)?;
        let entries = tetrahedra
            .len()
            .checked_mul(4)
            .ok_or_else(|| invalid("pressure face count overflow"))?;
        entries
            .checked_mul(81)
            .ok_or_else(|| invalid("pressure tangent work count overflow"))?;
        let geometry = TetLoadGeometry3::new(nodes_m, tetrahedra, budget, maximum_face_entries, cx)
            .map_err(|error| match error {
                LoadError3::Cancelled => HyperTetError::Cancelled,
                other => invalid(&other.to_string()),
            })?;
        if geometry.faces().is_empty() {
            return Err(invalid("pressure requires a nonempty complete exterior"));
        }
        let mut faces = Vec::with_capacity(geometry.faces().len());
        let mut edges: BTreeMap<[usize; 2], (u8, i8)> = BTreeMap::new();
        for face in geometry.faces() {
            checkpoint(cx)?;
            let mut ids = face.nodes;
            let [a, b, c] = ids.map(|node| nodes_m[node]);
            let e1 = difference(b, a)?;
            let e2 = difference(c, a)?;
            let scale = maximum_component(&[e1, e2]);
            let normal = cross(e1.map(|x| x / scale), e2.map(|x| x / scale));
            let alignment = dot(normal, face.outward_normal);
            if !alignment.is_finite() || alignment == 0.0 {
                return Err(invalid(
                    "pressure boundary orientation is not representable",
                ));
            }
            if alignment < 0.0 {
                ids.swap(1, 2);
            }
            for (a, b) in [(ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])] {
                checkpoint(cx)?;
                let (edge, orientation) = if a < b { ([a, b], 1) } else { ([b, a], -1) };
                let incidence = edges.entry(edge).or_insert((0, 0));
                if incidence.0 == 2 {
                    return Err(invalid("pressure exterior has a nonmanifold edge"));
                }
                incidence.0 += 1;
                incidence.1 += orientation;
            }
            faces.push(ids);
        }
        for &(count, orientation) in edges.values() {
            checkpoint(cx)?;
            if count != 2 || orientation != 0 {
                return Err(invalid(
                    "pressure exterior is not closed and consistently oriented",
                ));
            }
        }
        let mut orientations = Vec::with_capacity(tetrahedra.len());
        for tet in tetrahedra {
            checkpoint(cx)?;
            orientations.push(signed_volume(tet.map(|node| nodes_m[node]))?.signum());
        }
        let prepared = Self {
            nodes_m: nodes_m.to_vec(),
            tetrahedra: tetrahedra.to_vec(),
            orientations,
            faces,
        };
        checkpoint(cx)?;
        Ok(prepared)
    }

    /// Canonically ordered exterior triangles with outward reference winding.
    pub fn faces(&self) -> &[[usize; 3]] {
        &self.faces
    }

    /// Admit finite pressures and exact face coverage before a solve begins.
    /// A sparse patch is represented by explicit zeros on unselected faces.
    /// Returns whether all faces are unloaded, without a second unbounded scan.
    pub fn validate_pressure(
        &self,
        pressure: TetPressure<'_>,
        cx: &Cx<'_>,
    ) -> Result<bool, HyperTetError> {
        checkpoint(cx)?;
        let unloaded = match pressure {
            TetPressure::Uniform(value) => {
                finite(value)?;
                value == 0.0
            }
            TetPressure::PerFace(values) => {
                if values.len() != self.faces.len() {
                    return Err(invalid(
                        "pressure values must cover the canonical exterior faces exactly",
                    ));
                }
                let mut unloaded = true;
                for (index, &value) in values.iter().enumerate() {
                    if index % 1024 == 0 {
                        checkpoint(cx)?;
                    }
                    finite(value)?;
                    unloaded &= value == 0.0;
                }
                unloaded
            }
        };
        checkpoint(cx)?;
        Ok(unloaded)
    }

    /// Exact reference-coordinate bits and ordered connectivity match.
    pub fn matches_reference(&self, nodes_m: &[[f64; 3]], tetrahedra: &[[usize; 4]]) -> bool {
        self.tetrahedra == tetrahedra
            && self.nodes_m.len() == nodes_m.len()
            && self
                .nodes_m
                .iter()
                .flatten()
                .zip(nodes_m.iter().flatten())
                .all(|(a, b)| a.to_bits() == b.to_bits())
    }

    /// Current total material volume, m^3; element inversion/collapse refuses.
    pub fn volume_m3(&self, displacement_m: &[f64], cx: &Cx<'_>) -> Result<f64, HyperTetError> {
        let (_, volume) = self.positions_and_volume(displacement_m, cx)?;
        Ok(volume)
    }

    /// Current pressure forces, N, in full node-major Cartesian order.
    ///
    /// Pressure is a finite signed differential pressure; positive values
    /// compress the solid. A negative value prescribes outward traction, not
    /// a claim about the tensile capability of any particular fluid.
    pub fn force_n<'p>(
        &self,
        displacement_m: &[f64],
        pressure_pa: impl Into<TetPressure<'p>>,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, HyperTetError> {
        self.force_and_volume_loading(displacement_m, pressure_pa.into(), cx)
            .map(|(force, _)| force)
    }

    /// Shared load evaluation retains the volume for conservative static work.
    pub(super) fn force_and_volume<'p>(
        &self,
        displacement_m: &[f64],
        pressure_pa: impl Into<TetPressure<'p>>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), HyperTetError> {
        self.force_and_volume_loading(displacement_m, pressure_pa.into(), cx)
    }

    fn force_and_volume_loading(
        &self,
        displacement_m: &[f64],
        pressure: TetPressure<'_>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), HyperTetError> {
        self.validate_pressure(pressure, cx)?;
        let (positions, volume) = self.positions_and_volume(displacement_m, cx)?;
        let mut force = vec![0.0; displacement_m.len()];
        for (index, face) in self.faces.iter().enumerate() {
            checkpoint(cx)?;
            let pressure_pa = pressure.at(index);
            let [a, b, c] = face.map(|node| positions[node]);
            let area = area_vector(a, b, c)?;
            for &node in face {
                for axis in 0..3 {
                    let load = finite((-pressure_pa / 3.0) * area[axis])?;
                    force[3 * node + axis] = finite(force[3 * node + axis] + load)?;
                }
            }
        }
        checkpoint(cx)?;
        Ok((force, volume))
    }

    /// Analytic derivative of the external pressure force with respect to u.
    ///
    /// For an oriented face (a,b,c), each vertex receives -p/6 times
    /// (b-a) cross (c-a). Its variations are h cross (b-c), h cross
    /// (c-a), and h cross (a-b) for perturbations of a, b, and c.
    /// Uniform closed-boundary assembly is symmetric to rounding. Per-face
    /// assembly is generally nonsymmetric and must not be symmetrized. This
    /// is the EXTERNAL force Jacobian, which an internal residual subtracts.
    pub fn force_tangent<'p>(
        &self,
        displacement_m: &[f64],
        pressure_pa: impl Into<TetPressure<'p>>,
        cx: &Cx<'_>,
    ) -> Result<Csr, HyperTetError> {
        let pressure = pressure_pa.into();
        self.validate_pressure(pressure, cx)?;
        let (positions, _) = self.positions_and_volume(displacement_m, cx)?;
        let mut entries = Coo::new(displacement_m.len(), displacement_m.len());
        for (index, face) in self.faces.iter().enumerate() {
            checkpoint(cx)?;
            let pressure_pa = pressure.at(index);
            let [a, b, c] = face.map(|node| positions[node]);
            let variations = [difference(b, c)?, difference(c, a)?, difference(a, b)?];
            for (&column_node, edge) in face.iter().zip(variations) {
                for column_axis in 0..3 {
                    let mut basis = [0.0; 3];
                    basis[column_axis] = 1.0;
                    let derivative = cross(basis, edge);
                    for &row_node in face {
                        for (row_axis, value) in derivative.into_iter().enumerate() {
                            let value = finite((-pressure_pa / 6.0) * value)?;
                            if value != 0.0 {
                                entries.push(
                                    3 * row_node + row_axis,
                                    3 * column_node + column_axis,
                                    value,
                                );
                            }
                        }
                    }
                }
            }
        }
        let tangent = entries.assemble();
        for row in 0..tangent.nrows() {
            checkpoint(cx)?;
            for &value in tangent.row(row).1 {
                finite(value)?;
            }
        }
        Ok(tangent)
    }

    /// Work of held pressure over the declared displacement segment.
    ///
    /// Uniform pressure uses exact endpoint-volume work, independently of
    /// the path. Per-face pressure integrates force along the straight segment
    /// with Simpson's rule, exact for the quadratic P1 area vector to rounding.
    /// It retains common translation work on open patches. Subdividing a
    /// curved trajectory is necessary to resolve its physical path work.
    pub fn work_j<'p>(
        &self,
        before_m: &[f64],
        after_m: &[f64],
        pressure_pa: impl Into<TetPressure<'p>>,
        cx: &Cx<'_>,
    ) -> Result<f64, HyperTetError> {
        let pressure = pressure_pa.into();
        self.validate_pressure(pressure, cx)?;
        if let TetPressure::Uniform(value) = pressure {
            let before = self.volume_m3(before_m, cx)?;
            let after = self.volume_m3(after_m, cx)?;
            return finite(-value * finite(after - before)?);
        }
        let before_force = self.force_n(before_m, pressure, cx)?;
        let after_force = self.force_n(after_m, pressure, cx)?;
        let mut midpoint = Vec::with_capacity(before_m.len());
        for (index, (&before, &after)) in before_m.iter().zip(after_m).enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            midpoint.push(finite(0.5 * before + 0.5 * after)?);
        }
        let middle_force = self.force_n(&midpoint, pressure, cx)?;
        let mut work = 0.0;
        let mut compensation = 0.0;
        for index in 0..before_m.len() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            let average = finite(
                before_force[index] / 6.0
                    + (2.0 / 3.0) * middle_force[index]
                    + after_force[index] / 6.0,
            )?;
            let increment = finite(average * finite(after_m[index] - before_m[index])?)?;
            let corrected = finite(increment - compensation)?;
            let next = finite(work + corrected)?;
            compensation = (next - work) - corrected;
            work = next;
        }
        checkpoint(cx)?;
        finite(work)
    }

    fn positions_and_volume(
        &self,
        displacement_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 3]>, f64), HyperTetError> {
        self.positions_and_volume_about(displacement_m, 0, cx)
    }

    fn positions_and_volume_about(
        &self,
        displacement_m: &[f64],
        origin_node: usize,
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 3]>, f64), HyperTetError> {
        checkpoint(cx)?;
        if displacement_m.len() != 3 * self.nodes_m.len() {
            return Err(invalid(
                "pressure displacement dimension does not match the reference mesh",
            ));
        }
        // Pressure geometry depends only on relative positions. Remove the
        // common motion before adding it to reference coordinates, otherwise
        // a large rigid translation can round a small face or tet to zero.
        let reference_origin = self.nodes_m[origin_node];
        let displacement_origin = [
            finite(displacement_m[3 * origin_node])?,
            finite(displacement_m[3 * origin_node + 1])?,
            finite(displacement_m[3 * origin_node + 2])?,
        ];
        let mut positions = Vec::with_capacity(self.nodes_m.len());
        for (node, reference) in self.nodes_m.iter().enumerate() {
            checkpoint(cx)?;
            let mut position = [0.0; 3];
            for axis in 0..3 {
                let reference_offset = finite(reference[axis] - reference_origin[axis])?;
                let displacement_offset =
                    finite(finite(displacement_m[3 * node + axis])? - displacement_origin[axis])?;
                position[axis] = finite(reference_offset + displacement_offset)?;
            }
            positions.push(position);
        }
        let mut volume = 0.0;
        for (tet, &orientation) in self.tetrahedra.iter().zip(&self.orientations) {
            checkpoint(cx)?;
            let current = orientation * signed_volume(tet.map(|node| positions[node]))?;
            if current <= 0.0 {
                return Err(invalid(
                    "pressure volume contains an inverted or collapsed tetrahedron",
                ));
            }
            volume = finite(volume + current)?;
        }
        if volume <= 0.0 {
            return Err(invalid("pressure volume is not representable"));
        }
        checkpoint(cx)?;
        Ok((positions, volume))
    }
}

fn invalid(what: &str) -> HyperTetError {
    HyperTetError::InvalidInput {
        what: what.to_owned(),
    }
}

fn finite(value: f64) -> Result<f64, HyperTetError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("nonfinite pressure geometry or load arithmetic"))
    }
}

fn difference(a: [f64; 3], b: [f64; 3]) -> Result<[f64; 3], HyperTetError> {
    Ok([
        finite(a[0] - b[0])?,
        finite(a[1] - b[1])?,
        finite(a[2] - b[2])?,
    ])
}

fn maximum_component(vectors: &[[f64; 3]]) -> f64 {
    vectors
        .iter()
        .flatten()
        .fold(0.0, |largest, value| largest.max(value.abs()))
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn signed_volume(points: [[f64; 3]; 4]) -> Result<f64, HyperTetError> {
    let edges = [
        difference(points[1], points[0])?,
        difference(points[2], points[0])?,
        difference(points[3], points[0])?,
    ];
    let scale = maximum_component(&edges);
    if scale == 0.0 {
        return Err(invalid("zero-size pressure tetrahedron"));
    }
    let [a, b, c] = edges.map(|edge| edge.map(|value| value / scale));
    let volume = finite(((dot(a, cross(b, c)) / 6.0 * scale) * scale) * scale)?;
    if volume == 0.0 {
        return Err(invalid("collapsed or unrepresentable pressure tetrahedron"));
    }
    Ok(volume)
}

fn area_vector(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Result<[f64; 3], HyperTetError> {
    let edges = [difference(b, a)?, difference(c, a)?];
    let scale = maximum_component(&edges);
    if scale == 0.0 {
        return Err(invalid("zero-size pressure face"));
    }
    let [a, b] = edges.map(|edge| edge.map(|value| value / scale));
    let normal = cross(a, b);
    let mut area = [0.0; 3];
    for axis in 0..3 {
        area[axis] = finite((0.5 * normal[axis] * scale) * scale)?;
    }
    if area.iter().all(|value| *value == 0.0) {
        return Err(invalid("collapsed or unrepresentable pressure face"));
    }
    Ok(area)
}
