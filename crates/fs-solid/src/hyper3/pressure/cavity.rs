//! Sealed gas-cavity geometry on an admitted tetrahedral solid boundary.
//!
//! A complete solid exterior can contain both an outer surface and inner void
//! surfaces. Its material volume is NOT the gas volume. This owner selects one
//! connected, closed inner boundary, retains canonical solid-face indices for
//! pressure loading, and reverses the orientation only for gas-volume geometry.
//! Positive differential gas pressure therefore has force +p grad(V_gas).
//!
//! An evaluated geometry epoch exposes compact boundary-node coordinates,
//! dV/du, and an analytic Hessian-vector action. These can be consumed by
//! fs-material's GasMechanicalLoad without a fixed-area or fixed-gamma closure.
//! The material owner supplies p and dp/dV; this module supplies geometry only.
//! Geometry, force and tangent actions support the admitted mesh boundary size.
//! The gas tangent retains sparse surface geometry and one acoustic rank-one
//! factor; it neither imposes the small local-coordinate limit nor materializes
//! the dense compressibility coupling. Its reusable-buffer action allocates no
//! heap memory and polls cancellation in bounded coordinate/nonzero tiles.
//!
//! Every evaluation checks the complete solid mesh for local tetrahedron
//! inversion/collapse and requires positive enclosed cavity volume. Closure,
//! edge incidence and connectedness are checked at preparation. These checks
//! do not establish global injectivity, absence of self-intersection, contact,
//! or geometric validity of an arbitrary imported mesh. Pressure evolution,
//! gas heat exchange and time integration remain with their existing owners.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::thermochemical_gas::energy::entropy::{GasEntropyError, GasEntropySolve};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_math::dd::Dd;

mod tangent;
pub use tangent::TetCavityGasTangent;
pub(crate) mod equilibrium;
pub use equilibrium::AdiabaticCavity;

use super::{
    ClosedTetPressureSurface, HyperTetError, area_vector, checkpoint, cross, difference, dot,
    finite, invalid, maximum_component,
};

/// One closed inner boundary of an immutable reference tetrahedral solid.
///
/// The complete exterior is shared, so a generated face-pressure field can be
/// passed directly to the existing native static or dynamic pressure owner.
/// Other boundary faces receive zero; additional loads remain explicit.
#[derive(Debug, Clone)]
pub struct ClosedTetCavity {
    surface: Arc<ClosedTetPressureSurface>,
    solid_face_indices: Vec<usize>,
    nodes: Vec<usize>,
    // Gas-outward triangles, indexing nodes instead of the entire mesh.
    local_faces: Vec<[usize; 3]>,
    reference_positions_m: Vec<[f64; 3]>,
    reference_volume_m3: f64,
}

/// Immutable cavity geometry at one admitted solid displacement.
///
/// Coordinates are Cartesian displacements of [`ClosedTetCavity::nodes`] in
/// that exact node-major x,y,z order. A solid global DOF is 3*nodes[i]+axis.
/// Nonboundary coordinates have identically zero volume derivatives and are
/// deliberately omitted. Evaluation owns current positions, so later changes
/// to the input displacement cannot change this geometry epoch.
#[derive(Debug)]
pub struct TetCavityVolume<'a> {
    cavity: &'a ClosedTetCavity,
    positions_m: Vec<[f64; 3]>,
    volume_m3: f64,
    gradient_m2: Vec<f64>,
}

/// One source-resolved adiabatic gas load pulled back to the full solid chart.
/// The endpoint gas is reconstructed from the original parcel at the actual
/// cavity volume. The tangent includes both geometry and acoustic dp/dV terms.
#[derive(Debug)]
pub struct TetCavityGasLoad {
    /// Gas endpoint at exactly `volume_m3`, retaining the original source model.
    pub gas: ClosedGasParcel,
    /// Actual enclosed gas volume [m^3].
    pub volume_m3: f64,
    /// Gas pressure minus the fixed exterior pressure [Pa].
    pub net_pressure_pa: f64,
    /// Full node-major force on the surrounding solid [N].
    pub force_n: Vec<f64>,
    /// Full force derivative: sparse geometry plus acoustic rank-one coupling.
    /// Use apply_into with reusable buffers for large admitted cavity meshes.
    pub force_tangent: TetCavityGasTangent,
    /// Delta U_gas + p_exterior Delta V [J] from the reference cavity state.
    pub potential_change_j: f64,
    // Retain the already evaluated compact geometry for original-volume work
    // between Newton trials. No second mesh admission or position copy is needed.
    positions_m: Vec<[f64; 3]>,
}

impl ClosedTetCavity {
    /// Select one complete inner boundary by canonical solid-face indices.
    ///
    /// The input order is immaterial; duplicate or out-of-range indices refuse.
    /// Every selected edge must have exactly two opposite incidences and all
    /// faces must be connected through edges. Reversing the solid-outward
    /// orientation must enclose positive volume: selecting the outer body shell
    /// therefore refuses instead of silently treating material as gas.
    ///
    /// Work and allocations are bounded by the already admitted complete
    /// surface; no second mesh-sized face/tangent budget is introduced.
    pub fn new(
        surface: Arc<ClosedTetPressureSurface>,
        face_indices: &[usize],
        cx: &Cx<'_>,
    ) -> Result<Self, HyperTetError> {
        checkpoint(cx)?;
        if face_indices.is_empty() || face_indices.len() > surface.faces.len() {
            return Err(invalid(
                "a cavity needs a nonempty subset of boundary faces",
            ));
        }
        let mut selected = BTreeSet::new();
        for &index in face_indices {
            checkpoint(cx)?;
            if index >= surface.faces.len() || !selected.insert(index) {
                return Err(invalid("cavity face indices are repeated or out of range"));
            }
        }
        let solid_face_indices: Vec<_> = selected.into_iter().collect();
        let mut nodes = BTreeSet::new();
        let mut edges: BTreeMap<[usize; 2], (usize, usize, i8)> = BTreeMap::new();
        let mut neighbors = vec![[usize::MAX; 3]; solid_face_indices.len()];
        for (local, &index) in solid_face_indices.iter().enumerate() {
            checkpoint(cx)?;
            let [a, b, c] = surface.faces[index];
            for node in [a, b, c] {
                nodes.insert(node);
            }
            for (slot, (a, b)) in [(a, b), (b, c), (c, a)].into_iter().enumerate() {
                checkpoint(cx)?;
                let (edge, orientation) = if a < b {
                    ([a, b], 1_i8)
                } else {
                    ([b, a], -1_i8)
                };
                if let Some((other, other_slot, other_orientation)) = edges.remove(&edge) {
                    if orientation == other_orientation {
                        return Err(invalid("cavity edges are not consistently oriented"));
                    }
                    neighbors[local][slot] = other;
                    neighbors[other][other_slot] = local;
                } else {
                    edges.insert(edge, (local, slot, orientation));
                }
            }
        }
        if !edges.is_empty() {
            return Err(invalid("selected cavity boundary is open"));
        }
        let mut visited = vec![false; solid_face_indices.len()];
        let mut pending = vec![0_usize];
        visited[0] = true;
        let mut count = 0;
        while let Some(face) = pending.pop() {
            checkpoint(cx)?;
            count += 1;
            for &other in &neighbors[face] {
                if other == usize::MAX {
                    return Err(invalid("selected cavity boundary has incomplete incidence"));
                }
                if !visited[other] {
                    visited[other] = true;
                    pending.push(other);
                }
            }
        }
        if count != solid_face_indices.len() {
            return Err(invalid("a cavity must have one connected boundary"));
        }

        let nodes: Vec<_> = nodes.into_iter().collect();
        let mut local_faces = Vec::with_capacity(solid_face_indices.len());
        for &index in &solid_face_indices {
            checkpoint(cx)?;
            let [a, b, c] = surface.faces[index];
            let local = |node| {
                nodes
                    .binary_search(&node)
                    .map_err(|_| invalid("cavity face is missing a boundary node"))
            };
            // The solid outward normal on an inner void points INTO the gas.
            local_faces.push([local(a)?, local(c)?, local(b)?]);
        }
        let origin = surface.nodes_m[nodes[0]];
        let mut positions = Vec::with_capacity(nodes.len());
        for &node in &nodes {
            checkpoint(cx)?;
            positions.push(difference(surface.nodes_m[node], origin)?);
        }
        let reference_volume_m3 = enclosed_volume(&positions, &local_faces, cx)?;
        let cavity = Self {
            surface,
            solid_face_indices,
            nodes,
            local_faces,
            reference_positions_m: positions,
            reference_volume_m3,
        };
        checkpoint(cx)?;
        Ok(cavity)
    }

    /// Complete prepared exterior, suitable for the existing pressure solvers.
    pub fn surface(&self) -> &ClosedTetPressureSurface {
        &self.surface
    }

    /// Selected canonical indices into the complete solid exterior's faces.
    pub fn solid_face_indices(&self) -> &[usize] {
        &self.solid_face_indices
    }

    /// Sorted solid node IDs; derivatives use their compact x,y,z DOF order.
    pub fn nodes(&self) -> &[usize] {
        &self.nodes
    }

    /// Undeformed gas volume, excluding all solid material, m^3.
    pub const fn reference_volume_m3(&self) -> f64 {
        self.reference_volume_m3
    }

    /// Evaluate actual enclosed volume and its exact first derivative.
    ///
    /// Displacements cover the COMPLETE reference solid in node-major order.
    /// Reference and displacement offsets are formed about a cavity vertex
    /// before addition, preserving geometry under a large common translation.
    /// All material tetrahedra retain their reference orientation, including
    /// those away from this cavity. A zero or negative current gas volume refuses.
    pub fn evaluate(
        &self,
        displacement_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<TetCavityVolume<'_>, HyperTetError> {
        let (positions_m, volume_m3) = self.positions_and_volume(displacement_m, cx)?;
        let mut gradient_m2 = vec![0.0; 3 * self.nodes.len()];
        for face in &self.local_faces {
            checkpoint(cx)?;
            let [a, b, c] = face.map(|node| positions_m[node]);
            let area = area_vector(a, b, c)?;
            for &node in face {
                for axis in 0..3 {
                    let index = 3 * node + axis;
                    gradient_m2[index] = finite(gradient_m2[index] + area[axis] / 3.0)?;
                }
            }
        }
        checkpoint(cx)?;
        Ok(TetCavityVolume {
            cavity: self,
            positions_m,
            volume_m3,
            gradient_m2,
        })
    }

    fn positions_and_volume(
        &self, displacement_m: &[f64], cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 3]>, f64), HyperTetError> {
        let (positions, _) = self.surface
            .positions_and_volume_about(displacement_m, self.nodes[0], cx)?;
        let mut compact = Vec::with_capacity(self.nodes.len());
        for &node in &self.nodes {
            checkpoint(cx)?;
            compact.push(positions[node]);
        }
        let volume = enclosed_volume(&compact, &self.local_faces, cx)?;
        Ok((compact, volume))
    }

    /// Original gas-volume change, retaining changes hidden by rounded totals.
    ///
    /// Both complete solid configurations and positive cavity volumes are
    /// admitted first. Each oriented boundary cone then uses the same polarized
    /// cubic determinant as the pressure owner. This is endpoint volume work,
    /// not a sampled pressure-force integral or a material-volume difference.
    pub fn volume_change_m3(
        &self, before_m: &[f64], after_m: &[f64], cx: &Cx<'_>,
    ) -> Result<f64, HyperTetError> {
        let (before, _) = self.positions_and_volume(before_m, cx)?;
        let (after, _) = self.positions_and_volume(after_m, cx)?;
        finite(self.volume_increment(&before, &after, cx)?.to_f64())
    }

    fn volume_increment(
        &self, before: &[[f64; 3]], after: &[[f64; 3]], cx: &Cx<'_>,
    ) -> Result<Dd, HyperTetError> {
        checkpoint(cx)?;
        if before.len() != self.nodes.len() || after.len() != self.nodes.len() {
            return Err(invalid("cavity volume increment has the wrong coordinate chart"));
        }
        let mut change = Dd::ZERO;
        for &[a, b, c] in &self.local_faces {
            checkpoint(cx)?;
            change = change + super::signed_volume_increment(
                [[0.0; 3], before[a], before[b], before[c]],
                [[0.0; 3], after[a], after[b], after[c]],
            )?;
        }
        finite(change.to_f64())?;
        checkpoint(cx)?;
        Ok(change)
    }

    /// Pull the original closed gas through this cavity's exact nonlinear volume.
    ///
    /// The supplied parcel must be at the undeformed cavity volume. Every call
    /// starts from that immutable thermodynamic state, follows its source-owned
    /// reversible adiabatic path to the CURRENT cavity volume, and then uses the
    /// common acoustic mechanical-load owner for pressure and dp/dV. This avoids
    /// a fixed-gamma or frozen-pressure surrogate and keeps source, force,
    /// acoustic stiffness and solid geometry on one state.
    pub fn adiabatic_gas_load(
        &self,
        displacement_m: &[f64],
        gas: &ClosedGasParcel,
        exterior_pressure_pa: f64,
        controls: GasEntropySolve,
        cx: &Cx<'_>,
    ) -> Result<TetCavityGasLoad, HyperTetError> {
        checkpoint(cx)?;
        if gas.volume_m3().to_bits() != self.reference_volume_m3.to_bits()
            || !exterior_pressure_pa.is_finite()
            || exterior_pressure_pa < 0.0
        {
            return Err(invalid(
                "adiabatic cavity gas must start at the reference volume with finite nonnegative exterior pressure",
            ));
        }
        let geometry = self.evaluate(displacement_m, cx)?;
        let mut endpoint = gas.clone();
        let mut was_cancelled = false;
        let volume_step = endpoint.set_volume_isentropic(
            geometry.volume_m3,
            controls,
            || {
                if cx.checkpoint().is_err() {
                    was_cancelled = true;
                    true
                } else {
                    false
                }
            },
        );
        let volume_step = match volume_step {
            Ok(step) => step,
            Err(_) if was_cancelled => return Err(HyperTetError::Cancelled),
            Err(error) => return Err(invalid(&format!("adiabatic cavity gas refused state: {error}"))),
        };
        let map_error = |error: GasEntropyError| match error {
            GasEntropyError::Cancelled => HyperTetError::Cancelled,
            other => invalid(&format!("cavity gas mechanical load refused: {other}")),
        };
        let load = endpoint
            .mechanical_load_with_limit(
                geometry.volume_m3,
                geometry.gradient_m2(),
                exterior_pressure_pa,
                geometry.gradient_m2().len(),
                || cx.checkpoint().is_err(),
            )
            .map_err(map_error)?;
        let compact_force = load
            .generalized_forces_cancellable(|| cx.checkpoint().is_err())
            .map_err(map_error)?;
        let mut force_n = vec![0.0; displacement_m.len()];
        let mut volume_gradient = vec![0.0; displacement_m.len()];
        for (local, &value) in compact_force.iter().enumerate() {
            if local % 256 == 0 { checkpoint(cx)?; }
            let node = self.nodes[local / 3];
            let global = 3 * node + local % 3;
            force_n[global] = value;
            volume_gradient[global] = geometry.gradient_m2()[local];
        }
        let face_pressure = self.face_pressures_pa(load.net_pressure_pa(), cx)?;
        let geometric = self.surface.force_tangent(
            displacement_m, super::TetPressure::PerFace(&face_pressure), cx,
        )?;
        let tangent = TetCavityGasTangent::new(
            geometric,
            volume_gradient,
            load.pressure_response().pressure_volume_derivative_pa_m3,
            cx,
        )?;
        let exterior_work = if exterior_pressure_pa == 0.0 { Dd::ZERO } else {
            Dd::from_f64(exterior_pressure_pa) * self.volume_increment(
                &self.reference_positions_m, &geometry.positions_m, cx,
            )?
        };
        let potential_change_j = finite((Dd::from_f64(volume_step.stored_energy_change_j)
            + exterior_work).to_f64())?;
        let net_pressure_pa = load.net_pressure_pa();
        checkpoint(cx)?;
        Ok(TetCavityGasLoad {
            gas: endpoint,
            volume_m3: geometry.volume_m3,
            net_pressure_pa,
            force_n,
            force_tangent: tangent,
            potential_change_j,
            positions_m: geometry.positions_m,
        })
    }
    /// Uniform differential gas pressure in canonical SOLID face order.
    ///
    /// Positive p expands the cavity: it acts opposite the solid-outward normal
    /// on the inner boundary, exactly the existing positive-pressure convention.
    /// All other faces are unloaded. A negative value represents a larger
    /// counterpressure, not an assertion about a gas's tensile capability.
    /// Use absolute gas pressure if the exterior pressure acts on a different
    /// solid boundary; apply that spatially separate exterior load explicitly.
    pub fn face_pressures_pa(
        &self,
        differential_pressure_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, HyperTetError> {
        checkpoint(cx)?;
        finite(differential_pressure_pa)?;
        let mut values = vec![0.0; self.surface.faces.len()];
        for &index in &self.solid_face_indices {
            checkpoint(cx)?;
            values[index] = differential_pressure_pa;
        }
        checkpoint(cx)?;
        Ok(values)
    }

    /// Exact work ON the solid from held differential pressure: p*(V_after-V_before).
    ///
    /// This is endpoint-volume work, independent of the intervening admissible
    /// path. It is positive for gas expansion. Changing gas pressure requires
    /// the material energy/entropy owner; multiplying one sampled p by delta V
    /// is not a finite-step adiabatic gas law.
    pub fn held_pressure_work_j(
        &self,
        before_m: &[f64],
        after_m: &[f64],
        differential_pressure_pa: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, HyperTetError> {
        checkpoint(cx)?;
        finite(differential_pressure_pa)?;
        let (before, _) = self.positions_and_volume(before_m, cx)?;
        let (after, _) = self.positions_and_volume(after_m, cx)?;
        let change = self.volume_increment(&before, &after, cx)?;
        finite((Dd::from_f64(differential_pressure_pa) * change).to_f64())
    }
}

impl TetCavityVolume<'_> {
    /// Current enclosed gas volume, m^3.
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }

    /// dV/du, m^2, in the compact boundary-node chart.
    ///
    /// Multiplying by differential pressure gives nodal force on the solid.
    /// Dotting with boundary-node velocity gives the actual gas-volume rate.
    pub fn gradient_m2(&self) -> &[f64] {
        &self.gradient_m2
    }

    /// Sorted solid node IDs for this geometry epoch's compact derivative chart.
    pub fn nodes(&self) -> &[usize] {
        &self.cavity.nodes
    }

    /// Exact Hess(V)*direction, without building a sparse or dense matrix.
    ///
    /// The direction covers precisely the compact boundary-node coordinates.
    /// For a displacement direction in m the result is in m^2; for a velocity
    /// direction in m/s it is in m^2/s. Every triangle contributes the derivative
    /// of its current area vector, in O(number of boundary faces) work.
    /// Use the material pressure derivative as well as this geometric term
    /// when linearizing a sealed, compressible gas load.
    pub fn hessian_action(
        &self,
        direction: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, HyperTetError> {
        checkpoint(cx)?;
        if direction.len() != self.gradient_m2.len() {
            return Err(invalid(
                "cavity Hessian direction has the wrong coordinate count",
            ));
        }
        for (index, &value) in direction.iter().enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            finite(value)?;
        }
        let mut action = vec![0.0; direction.len()];
        for face in &self.cavity.local_faces {
            checkpoint(cx)?;
            let [a, b, c] = face.map(|node| self.positions_m[node]);
            let [va, vb, vc] = face.map(|node| {
                [
                    direction[3 * node],
                    direction[3 * node + 1],
                    direction[3 * node + 2],
                ]
            });
            let first = cross(difference(vb, va)?, difference(c, a)?);
            let second = cross(difference(b, a)?, difference(vc, va)?);
            for &node in face {
                for axis in 0..3 {
                    let derivative =
                        finite(finite(first[axis])? / 6.0 + finite(second[axis])? / 6.0)?;
                    let index = 3 * node + axis;
                    action[index] = finite(action[index] + derivative)?;
                }
            }
        }
        checkpoint(cx)?;
        Ok(action)
    }
}

fn enclosed_volume(
    positions: &[[f64; 3]],
    faces: &[[usize; 3]],
    cx: &Cx<'_>,
) -> Result<f64, HyperTetError> {
    let mut volume = 0.0;
    let mut compensation = 0.0;
    for face in faces {
        checkpoint(cx)?;
        let points = face.map(|node| positions[node]);
        let scale = maximum_component(&points);
        if scale == 0.0 {
            return Err(invalid("zero-size cavity triangle"));
        }
        let [a, b, c] = points.map(|point| point.map(|value| value / scale));
        // Cones through the origin legitimately contribute zero. A concave
        // cavity can have negative cones; only the oriented total must be > 0.
        let contribution = finite(((dot(a, cross(b, c)) / 6.0 * scale) * scale) * scale)?;
        let corrected = finite(contribution - compensation)?;
        let next = finite(volume + corrected)?;
        compensation = (next - volume) - corrected;
        volume = next;
    }
    if volume <= 0.0 {
        return Err(invalid(
            "cavity must enclose positive gas volume with gas-outward orientation",
        ));
    }
    checkpoint(cx)?;
    Ok(volume)
}
