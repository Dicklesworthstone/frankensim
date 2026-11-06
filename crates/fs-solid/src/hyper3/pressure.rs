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
//! An affine Eulerian field instead follows each CURRENT point and is integrated
//! exactly against the P1 surface basis. Its closed-surface potential is
//! `integral_volume p(x) dV`; the pressure gradient therefore supplies buoyancy
//! without replacing the physical surface traction by a separate body force.
//! A prescribed planar free surface clips this field at zero gauge pressure,
//! integrating wet-face traction and clipped-volume work exactly. Prepared
//! stratified reservoirs extend this same potential to depth-dependent density,
//! retaining every layer interface in forces, tangents and endpoint work.
//! There is no global injectivity certificate, fluid equation of state, or
//! time integrator in this module.

use std::collections::BTreeMap;

use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_sparse::{Coo, Csr};

use super::{HyperTetError, checkpoint};
use crate::linear3::TetAssemblyBudget;
use crate::loads3::{LoadError3, TetLoadGeometry3};

pub(super) mod algorithmic;
/// Enclosed void volume and work-conjugate pressure geometry for gas cavities.
pub mod cavity;
mod free_surface;
mod stratified;

pub use stratified::{HydrostaticDensityLayer, StratifiedHydrostaticPressure};

use free_surface::FaceMoments;

/// Shared exact wet-face integration for reference and follower pressure loads.
/// Values are UNCLIPPED P1 gauge-pressure coefficients; clipping is spatial.
pub(crate) fn free_surface_face_pressure_weights(
    values: [f64; 3],
) -> Result<[f64; 3], HyperTetError> {
    for value in values {
        finite(value)?;
    }
    Ok(free_surface::face_moments::<false>(values)?.weights_pa)
}

/// Signed differential pressure in Pa, held over the declared load interval.
/// Positive pressure acts inward; negative pressure prescribes outward traction.
#[derive(Debug, Clone, Copy)]
pub enum TetPressure<'a> {
    /// One pressure on the complete closed boundary, with potential `p V`.
    Uniform(f64),
    /// One pressure per canonical [`ClosedTetPressureSurface::faces`] entry.
    /// Zero leaves a face unloaded. No conservative potential is assumed.
    PerFace(&'a [f64]),
    /// Spatially fixed field `p(x) = pressure_pa + gradient_pa_m dot (x-origin_m)`.
    /// Hydrostatics in a uniform fluid uses `gradient_pa_m = rho_fluid * gravity`.
    /// This loads the complete closed boundary: no wet/dry clipping, free surface,
    /// cavitation, or equation of state is inferred. Signed values are traction
    /// data, so the caller remains responsible for the field's physical domain.
    Affine {
        /// World position at which `pressure_pa` is prescribed [m].
        origin_m: [f64; 3],
        /// Signed differential pressure at `origin_m` [Pa].
        pressure_pa: f64,
        /// World pressure gradient [Pa/m].
        gradient_pa_m: [f64; 3],
    },
    /// Prescribed planar free surface with gauge pressure
    /// `max(0, pressure_pa + gradient_pa_m dot (x-origin_m))`.
    ///
    /// Uniform-fluid hydrostatics uses `gradient_pa_m = rho_fluid * gravity`;
    /// the zero plane is the prescribed liquid surface and the other side is
    /// unloaded. Atmospheric/reference pressure is already subtracted. This
    /// does not evolve the fluid, model wetting/surface tension, or implement
    /// cavitation. Zero gradient reduces to uniform nonnegative pressure.
    FreeSurface {
        /// World reference position [m].
        origin_m: [f64; 3],
        /// Unclipped gauge pressure at `origin_m` [Pa].
        pressure_pa: f64,
        /// World hydrostatic pressure gradient [Pa/m].
        gradient_pa_m: [f64; 3],
    },
    /// Prescribed stratified reservoir, with continuous hydrostatic gauge
    /// pressure and explicit density interfaces. Nonnegative load scaling
    /// multiplies pressure AFTER clipping and never moves those interfaces.
    Stratified {
        /// Prepared physical density profile and world-frame gravity.
        profile: &'a StratifiedHydrostaticPressure,
        /// Finite nonnegative load multiplier; physical loading uses one.
        load_factor: f64,
    },
}

impl From<f64> for TetPressure<'_> {
    fn from(pressure_pa: f64) -> Self {
        Self::Uniform(pressure_pa)
    }
}

impl TetPressure<'_> {
    fn constant_at(self, face: usize) -> Option<f64> {
        match self {
            Self::Uniform(value) => Some(value),
            Self::PerFace(values) => Some(values[face]),
            Self::Affine {
                pressure_pa,
                gradient_pa_m,
                ..
            } if gradient_pa_m == [0.0; 3] => Some(pressure_pa),
            Self::FreeSurface {
                pressure_pa,
                gradient_pa_m,
                ..
            } if gradient_pa_m == [0.0; 3] => Some(pressure_pa.max(0.0)),
            Self::Stratified {
                profile,
                load_factor,
            } if load_factor == 0.0 || profile.is_zero() => Some(0.0),
            Self::Affine { .. } | Self::FreeSurface { .. } | Self::Stratified { .. } => None,
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
            TetPressure::Stratified {
                profile,
                load_factor,
            } => {
                finite(load_factor)?;
                if load_factor < 0.0 {
                    return Err(invalid(
                        "stratified pressure load factor must be nonnegative",
                    ));
                }
                load_factor == 0.0 || profile.is_zero()
            }
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
            TetPressure::Affine {
                origin_m,
                pressure_pa,
                gradient_pa_m,
            }
            | TetPressure::FreeSurface {
                origin_m,
                pressure_pa,
                gradient_pa_m,
            } => {
                finite(pressure_pa)?;
                for value in origin_m.into_iter().chain(gradient_pa_m) {
                    finite(value)?;
                }
                gradient_pa_m == [0.0; 3]
                    && (pressure_pa == 0.0
                        || matches!(pressure, TetPressure::FreeSurface { .. }) && pressure_pa < 0.0)
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
            let weights = self
                .face_pressure_moments::<false>(displacement_m, index, pressure, cx)?
                .weights_pa;
            let [a, b, c] = face.map(|node| positions[node]);
            let area = area_vector(a, b, c)?;
            for (&node, weight) in face.iter().zip(weights) {
                for axis in 0..3 {
                    let load = finite(-weight * area[axis])?;
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
    /// assembly is generally nonsymmetric and must not be symmetrized. Affine
    /// fields include both area and spatial-pressure derivatives; their complete
    /// closed-boundary assembly is also symmetric. This is the EXTERNAL force
    /// Jacobian, which an internal residual subtracts.
    ///
    /// A free surface integrates both terms over the current wet polygon. Its
    /// moving-edge term vanishes at zero pressure. A whole face exactly on the
    /// plane has no unique two-sided tangent; the dry-side generalized tangent
    /// is selected. Isolated zero-pressure edges/vertices have zero measure.
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
            let moments =
                self.face_pressure_moments::<true>(displacement_m, index, pressure, cx)?;
            let constant_factor = pressure.constant_at(index).map(|value| -value / 6.0);
            let gradient = match pressure {
                TetPressure::Affine { gradient_pa_m, .. }
                | TetPressure::FreeSurface { gradient_pa_m, .. } => gradient_pa_m,
                TetPressure::Stratified {
                    profile,
                    load_factor,
                } if load_factor != 0.0 => profile.depth_gradient(),
                _ => [0.0; 3],
            };
            let [a, b, c] = face.map(|node| positions[node]);
            let area = area_vector(a, b, c)?;
            let variations = [difference(b, c)?, difference(c, a)?, difference(a, b)?];
            for (column_local, (&column_node, edge)) in face.iter().zip(variations).enumerate() {
                for column_axis in 0..3 {
                    let mut basis = [0.0; 3];
                    basis[column_axis] = 1.0;
                    let derivative = cross(basis, edge);
                    for (row_local, (&row_node, weight)) in
                        face.iter().zip(moments.weights_pa).enumerate()
                    {
                        let factor = constant_factor.unwrap_or(-0.5 * weight);
                        let pressure_slope =
                            gradient[column_axis] * moments.basis[row_local][column_local];
                        for (row_axis, value) in derivative.into_iter().enumerate() {
                            // The basis moment is over the actual wet portion,
                            // normalized by the full current triangle area.
                            let value = finite(factor * value - pressure_slope * area[row_axis])?;
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
    /// Uniform pressure uses exact endpoint-volume work. Affine, free-surface
    /// and stratified fields use differences of their original volume moments,
    /// retaining small volume changes before absolute potentials are rounded.
    /// Free surfaces use `integral max(0,p(x)) dV`, including changes in immersion.
    /// A held per-face pressure integrates its quadratic force along the
    /// straight segment with Simpson's rule, retaining common translation work
    /// on loaded patches. Subdividing a curved trajectory is necessary for that
    /// path-dependent face field. No fixed time quadrature is assumed exact
    /// across a wet/dry transition.
    pub fn work_j<'p>(
        &self,
        before_m: &[f64],
        after_m: &[f64],
        pressure_pa: impl Into<TetPressure<'p>>,
        cx: &Cx<'_>,
    ) -> Result<f64, HyperTetError> {
        let pressure = pressure_pa.into();
        self.validate_pressure(pressure, cx)?;
        let uniform = match pressure {
            TetPressure::Stratified {
                profile,
                load_factor,
            } if load_factor == 0.0 || profile.is_zero() => Some(0.0),
            TetPressure::Uniform(value) => Some(value),
            TetPressure::Affine {
                pressure_pa,
                gradient_pa_m,
                ..
            } if gradient_pa_m == [0.0; 3] => Some(pressure_pa),
            TetPressure::FreeSurface {
                pressure_pa,
                gradient_pa_m,
                ..
            } if gradient_pa_m == [0.0; 3] => Some(pressure_pa.max(0.0)),
            _ => None,
        };
        if let Some(value) = uniform {
            // Keep the original endpoint geometry/domain checks, but form
            // each physical volume increment before a large unrelated mesh
            // volume can round it away in either global diagnostic total.
            let (before, _) = self.positions_and_volume(before_m, cx)?;
            let (after, _) = self.positions_and_volume(after_m, cx)?;
            let mut change = Dd::ZERO;
            for (tet, &orientation) in self.tetrahedra.iter().zip(&self.orientations) {
                checkpoint(cx)?;
                change = change
                    + Dd::from_f64(orientation)
                        * signed_volume_increment(
                            tet.map(|node| before[node]),
                            tet.map(|node| after[node]),
                        )?;
            }
            checkpoint(cx)?;
            return finite((-Dd::from_f64(value) * change).to_f64());
        }
        if let TetPressure::Stratified {
            profile,
            load_factor,
        } = pressure
        {
            return profile.work_j(self, before_m, after_m, load_factor, cx);
        }
        if let TetPressure::Affine {
            origin_m,
            pressure_pa,
            gradient_pa_m,
        }
        | TetPressure::FreeSurface {
            origin_m,
            pressure_pa,
            gradient_pa_m,
        } = pressure
        {
            return self.spatial_potential_work_j(
                before_m,
                after_m,
                origin_m,
                pressure_pa,
                gradient_pa_m,
                matches!(pressure, TetPressure::FreeSurface { .. }),
                cx,
            );
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

    /// Exact pressure and, when requested, wet basis-product moments.
    /// The field is sampled in world coordinates; translated relative geometry
    /// from `positions_and_volume` is deliberately not used as the field origin.
    fn face_pressure_moments<const TANGENT: bool>(
        &self,
        displacement_m: &[f64],
        face_index: usize,
        pressure: TetPressure<'_>,
        cx: &Cx<'_>,
    ) -> Result<FaceMoments, HyperTetError> {
        if let Some(value) = pressure.constant_at(face_index) {
            return Ok(FaceMoments::constant(value));
        }
        if let TetPressure::Stratified {
            profile,
            load_factor,
        } = pressure
        {
            return profile.face_moments::<TANGENT>(
                self,
                displacement_m,
                face_index,
                load_factor,
                cx,
            );
        }
        let (origin_m, pressure_pa, gradient_pa_m) = match pressure {
            TetPressure::Affine {
                origin_m,
                pressure_pa,
                gradient_pa_m,
            }
            | TetPressure::FreeSurface {
                origin_m,
                pressure_pa,
                gradient_pa_m,
            } => (origin_m, pressure_pa, gradient_pa_m),
            _ => unreachable!("constant pressure handled above"),
        };
        let mut values = [0.0; 3];
        for (local, &node) in self.faces[face_index].iter().enumerate() {
            values[local] = self.spatial_pressure_at(
                displacement_m,
                node,
                origin_m,
                pressure_pa,
                gradient_pa_m,
            )?;
        }
        if matches!(pressure, TetPressure::FreeSurface { .. }) {
            free_surface::face_moments::<TANGENT>(values)
        } else {
            FaceMoments::affine(values)
        }
    }

    pub(super) fn spatial_pressure_at(
        &self,
        displacement_m: &[f64],
        node: usize,
        origin_m: [f64; 3],
        pressure_pa: f64,
        gradient_pa_m: [f64; 3],
    ) -> Result<f64, HyperTetError> {
        let mut value = pressure_pa;
        for axis in 0..3 {
            let offset = finite(
                finite(self.nodes_m[node][axis] - origin_m[axis])?
                    + displacement_m[3 * node + axis],
            )?;
            value = finite(gradient_pa_m[axis].mul_add(offset, value))?;
        }
        Ok(value)
    }

    /// Conservative endpoint work from affine or positive clipped-volume moments.
    /// Summing element increments avoids subtracting large global potentials;
    /// the symmetric product difference avoids large endpoint energy products.
    #[allow(clippy::too_many_arguments)]
    fn spatial_potential_work_j(
        &self,
        before_m: &[f64],
        after_m: &[f64],
        origin_m: [f64; 3],
        pressure_pa: f64,
        gradient_pa_m: [f64; 3],
        clipped: bool,
        cx: &Cx<'_>,
    ) -> Result<f64, HyperTetError> {
        let (before, _) = self.positions_and_volume(before_m, cx)?;
        let (after, _) = self.positions_and_volume(after_m, cx)?;
        let mut before_pressure = Vec::with_capacity(self.nodes_m.len());
        let mut after_pressure = Vec::with_capacity(self.nodes_m.len());
        for node in 0..self.nodes_m.len() {
            checkpoint(cx)?;
            before_pressure.push(self.spatial_pressure_at(
                before_m,
                node,
                origin_m,
                pressure_pa,
                gradient_pa_m,
            )?);
            after_pressure.push(self.spatial_pressure_at(
                after_m,
                node,
                origin_m,
                pressure_pa,
                gradient_pa_m,
            )?);
        }
        let mut work = 0.0;
        let mut compensation = 0.0;
        for (tet, &orientation) in self.tetrahedra.iter().zip(&self.orientations) {
            checkpoint(cx)?;
            let before_values = tet.map(|node| before_pressure[node]);
            let after_values = tet.map(|node| after_pressure[node]);
            let mean = |values: [f64; 4]| {
                if clipped {
                    free_surface::tet_mean(values)
                } else {
                    finite(values.iter().fold(Dd::ZERO, |sum, &value| {
                        sum + Dd::from_f64(0.25) * Dd::from_f64(value)
                    }).to_f64())
                }
            };
            let mean_before = mean(before_values)?;
            let mean_after = mean(after_values)?;
            if clipped && mean_before == 0.0 && mean_after == 0.0 {
                continue;
            }
            let volume_before = orientation * signed_volume(tet.map(|node| before[node]))?;
            let volume_after = orientation * signed_volume(tet.map(|node| after[node]))?;
            let mean_change = if !clipped || before_values
                .iter()
                .chain(&after_values)
                .all(|&value| value >= 0.0)
            {
                // In a fully wet element, the datum cancels analytically. This
                // retains small translation work even at large ambient depth.
                let mut change = 0.0;
                for &node in tet {
                    for axis in 0..3 {
                        change = finite((gradient_pa_m[axis] / 4.0).mul_add(
                            finite(after_m[3 * node + axis] - before_m[3 * node + axis])?,
                            change,
                        ))?;
                    }
                }
                change
            } else {
                finite(mean_after - mean_before)?
            };
            let mean_volume = finite(0.5 * volume_before + 0.5 * volume_after)?;
            let mean_pressure = finite(0.5 * mean_before + 0.5 * mean_after)?;
            let volume_change = Dd::from_f64(orientation)
                * signed_volume_increment(
                    tet.map(|node| before[node]),
                    tet.map(|node| after[node]),
                )?;
            let increment = finite(
                (-Dd::from_f64(mean_pressure) * volume_change
                    - Dd::from_f64(mean_volume) * Dd::from_f64(mean_change))
                .to_f64(),
            )?;
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

/// Difference of the original cubic volume potential at represented endpoints.
/// Keep the midpoint and endpoint subtraction tails before evaluating the
/// polarized determinant; an independently rounded midpoint loses those terms.
pub(super) fn determinant_increment(before: &[f64; 9], after: &[f64; 9]) -> Dd {
    determinant_increment_dd(&before.map(Dd::from_f64), &after.map(Dd::from_f64))
}

fn determinant_increment_dd(before: &[Dd; 9], after: &[Dd; 9]) -> Dd {
    let half = Dd::from_f64(0.5);
    let midpoint = core::array::from_fn::<_, 9, _>(|i| half * before[i] + half * after[i]);
    let delta = core::array::from_fn::<_, 9, _>(|i| after[i] - before[i]);
    let cofactor = |f: &[Dd; 9], ij: usize| {
        let (a, b) = ((ij / 3 + 1) % 3, (ij / 3 + 2) % 3);
        let (c, d) = ((ij % 3 + 1) % 3, (ij % 3 + 2) % 3);
        f[3 * a + c] * f[3 * b + d] - f[3 * a + d] * f[3 * b + c]
    };
    let linear = delta.iter().enumerate().fold(Dd::ZERO, |sum, (i, &value)| {
        sum + cofactor(&midpoint, i) * value
    });
    let cubic = (0..3).fold(Dd::ZERO, |sum, i| sum + delta[i] * cofactor(&delta, i));
    linear + Dd::from_f64(0.25) * cubic
}

fn signed_volume_increment(
    before: [[f64; 3]; 4],
    after: [[f64; 3]; 4],
) -> Result<Dd, HyperTetError> {
    // Columns are physical edges. Scaling in compensated precision preserves
    // the original endpoint differences while avoiding cubic overflow.
    let edges = |points: &[[f64; 3]; 4]| {
        core::array::from_fn::<_, 9, _>(|ij| {
            Dd::from_f64(points[ij % 3 + 1][ij / 3]) - Dd::from_f64(points[0][ij / 3])
        })
    };
    let before = edges(&before);
    let after = edges(&after);
    let scale = before
        .iter()
        .chain(&after)
        .fold(0.0_f64, |largest, edge| largest.max(edge.hi.abs()));
    if !scale.is_finite() || scale == 0.0 {
        return Err(invalid(
            "unrepresentable pressure volume increment geometry",
        ));
    }
    let scale = Dd::from_f64(scale);
    let change = determinant_increment_dd(
        &before.map(|edge| edge / scale),
        &after.map(|edge| edge / scale),
    );
    let change = (((change / Dd::from_f64(6.0)) * scale) * scale) * scale;
    finite(change.to_f64())?;
    Ok(change)
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

#[cfg(test)]
mod potential_increment_tests {
    use super::{determinant_increment, signed_volume, signed_volume_increment};

    const UNIT: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];

    #[test]
    fn cubic_increment_retains_the_represented_midpoint_tail() {
        let before = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let mut after = before;
        let epsilon = f64::EPSILON;
        after[0] += epsilon;
        after[4] -= epsilon;
        // (1+epsilon)(1-epsilon)-1 = -epsilon². Both absolute
        // determinants round to one, and an f64 midpoint loses half the work.
        assert_eq!(after[0] * after[4], 1.0);
        assert_eq!(
            determinant_increment(&before, &after).to_f64(),
            -epsilon * epsilon
        );
        assert_eq!(
            determinant_increment(&after, &before).to_f64(),
            epsilon * epsilon
        );
        assert_eq!(determinant_increment(&after, &after).to_f64(), 0.0);
    }

    #[test]
    fn physical_volume_increment_retains_quadratic_work_below_endpoint_roundoff() {
        let shear = 2.0_f64.powi(-28);
        let mut after = UNIT;
        after[1][1] = shear;
        after[2][0] = shear;
        assert_eq!(signed_volume(UNIT).unwrap(), signed_volume(after).unwrap());
        let expected = -shear * shear / 6.0;
        let change = signed_volume_increment(UNIT, after).unwrap().to_f64();
        assert!((change / expected - 1.0).abs() < 4.0 * f64::EPSILON);
        assert_eq!(
            signed_volume_increment(after, UNIT).unwrap().to_f64(),
            -change
        );
    }

    #[test]
    fn physical_volume_increment_preserves_scale_and_orientation() {
        for scale in [1e-90, 1.0, 1e90] {
            let before = UNIT.map(|point| point.map(|coordinate| coordinate * scale));
            let mut after = before;
            after[1][0] *= 1.125;
            // Independent axial extrusion: area * the represented edge travel / 6.
            let expected = (((after[1][0] - before[1][0]) / 6.0) * scale) * scale;
            let change = signed_volume_increment(before, after).unwrap().to_f64();
            assert!((change / expected - 1.0).abs() < 16.0 * f64::EPSILON);
            let reversed_before = [before[0], before[2], before[1], before[3]];
            let reversed_after = [after[0], after[2], after[1], after[3]];
            let reversed = signed_volume_increment(reversed_before, reversed_after)
                .unwrap()
                .to_f64();
            assert!((reversed / expected + 1.0).abs() < 16.0 * f64::EPSILON);
        }
    }
}
