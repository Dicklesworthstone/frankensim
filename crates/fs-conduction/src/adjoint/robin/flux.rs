//! Original P1 prescribed-flux controls, shared by steady and transient adjoints.

use std::collections::BTreeMap;
use super::{ConductionError, Cx, ScalarField, ThermalBc, checked, invalid, poll, vector};
use crate::{ConductionMesh, ThermalBoundary};

/// Derivative in the original boundary's flux coordinate, in objective units
/// per W/m². A world vector is never collapsed to an averaged face normal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PrescribedFluxGradient {
    /// Uniform additive outward-normal flux, including over a fixed nodal profile.
    Normal(f64),
    /// Independent world x/y/z components of the original uniform vector.
    WorldVector([f64; 3]),
}

struct Face {
    control: usize,
    vertices: [usize; 3],
    /// P1 face integral of each active coordinate, before the outward-load minus.
    weights: [f64; 3],
}

/// Sparse original-face basis for independently controlled Neumann regions.
///
/// `load_direction` feeds the existing [`super::RobinResponse::apply`] with
/// nodal watts; `pullback` consumes the SAME complete nodal-load adjoint from
/// steady, capacity, enthalpy, contact or radiation owners. State feedback and
/// the time-step convention belong to that owner; no extra solve, heat source,
/// storage factor or frozen-feedback approximation is introduced here.
/// Geometry, the selected partition and prescribed nodes stay fixed. Scenario
/// frame/schedule parameters require their own chain rule; these coordinates
/// are the native sampled WORLD fluxes. Explicit objective dependence on flux
/// must be added separately. This is a discrete derivative, not an error bound.
///
/// Storage is linear in selected faces, controls and prescribed vertices;
/// there is no per-control full nodal vector or dense feedback matrix.
pub struct PrescribedFluxControls {
    vertices: usize,
    vectors: Vec<bool>,
    faces: Vec<Face>,
    fixed: Vec<usize>,
}

impl PrescribedFluxControls {
    /// Bind nonempty original scalar/vector Neumann regions in caller order.
    /// `max_boundary_faces` bounds the full face traversal and retained basis.
    /// Wrong kinds, duplicate/empty regions, mesh-size mismatch, invalid sizes,
    /// nonfinite arithmetic, allocation refusal and cancellation are errors.
    pub fn new(cx: &Cx<'_>, mesh: &ConductionMesh, boundary: &ThermalBoundary,
        regions: &[&str], max_boundary_faces: usize) -> Result<Self, ConductionError>
    {
        poll(cx, 0)?;
        if boundary.face_count() != mesh.boundary().len()
            || mesh.boundary().len() > max_boundary_faces || regions.len() > max_boundary_faces
        { return Err(invalid("prescribed-flux basis exceeds its face allowance or mesh partition")); }
        let mut selected = BTreeMap::new();
        let mut vectors = Vec::new();
        vectors.try_reserve_exact(regions.len()).map_err(|_| invalid("flux control allocation refused"))?;
        for (control, &name) in regions.iter().enumerate() {
            poll(cx, control)?;
            let region = boundary.region_names().iter().position(|v| v == name)
                .ok_or_else(|| invalid("unknown prescribed-flux control region"))?;
            if selected.insert(region, control).is_some() {
                return Err(invalid("duplicate prescribed-flux control region"));
            }
            vectors.push(match boundary.conditions()[region] {
                ThermalBc::Neumann { .. } => false,
                ThermalBc::NeumannVector { .. } => true,
                _ => return Err(invalid("prescribed-flux controls require original Neumann regions")),
            });
        }
        let mut faces = Vec::new();
        let mut counts = vec![0usize; regions.len()];
        for (slot, face) in mesh.boundary().iter().enumerate() {
            poll(cx, slot)?;
            let Some(control) = boundary.region_for(slot).and_then(|r| selected.get(&r)).copied()
                else { continue; };
            let weights = if vectors[control] {
                let mut w = [0.0; 3];
                for (i, value) in w.iter_mut().enumerate() {
                    *value = checked((face.area / 3.0) * face.outward_normal[i])?;
                }
                w
            } else { [checked(face.area / 3.0)?, 0.0, 0.0] };
            faces.try_reserve(1).map_err(|_| invalid("flux face allocation refused"))?;
            faces.push(Face { control, vertices: face.vertices.map(|v| v as usize), weights });
            counts[control] += 1;
        }
        if counts.contains(&0) { return Err(invalid("prescribed-flux control has no original face")); }
        let mut fixed = Vec::new();
        fixed.try_reserve_exact(boundary.dirichlet().len())
            .map_err(|_| invalid("flux prescribed-node allocation refused"))?;
        for &(vertex, _) in boundary.dirichlet() {
            poll(cx, vertex)?;
            if vertex >= mesh.vertex_count() { return Err(invalid("flux prescribed node is outside the mesh")); }
            fixed.push(vertex);
        }
        poll(cx, mesh.boundary().len())?;
        Ok(Self { vertices: mesh.vertex_count(), vectors, faces, fixed })
    }

    /// Number of original mesh vertices, including prescribed equations.
    #[must_use]
    pub const fn vertex_count(&self) -> usize { self.vertices }

    /// Assemble a full nodal-watt perturbation from original native carriers.
    /// Each selected scalar region accepts a uniform normal-flux direction;
    /// each vector region accepts a WORLD vector direction. Other kinds refuse.
    /// Assign this to `RobinDirection::nodal_load_w` before the native tangent.
    /// Prescribed equations remain zero; their reaction is not differentiated.
    pub fn load_direction(&self, cx: &Cx<'_>, directions: &[ThermalBc])
        -> Result<Vec<f64>, ConductionError>
    {
        poll(cx, 0)?;
        if directions.len() != self.vectors.len() { return Err(invalid("flux directions lost their control order")); }
        let mut coordinates = Vec::new();
        coordinates.try_reserve_exact(directions.len()).map_err(|_| invalid("flux direction allocation refused"))?;
        for (i, direction) in directions.iter().enumerate() {
            poll(cx, i)?;
            let value = match (self.vectors[i], direction) {
                (false, ThermalBc::Neumann { outward_flux: ScalarField::Uniform(q) }) => [*q, 0.0, 0.0],
                (true, ThermalBc::NeumannVector { flux }) => *flux,
                _ => return Err(invalid("flux direction must retain the original scalar/vector coordinate")),
            };
            vector(cx, &value, 3)?;
            coordinates.push(value);
        }
        let mut load = Vec::new();
        load.try_reserve_exact(self.vertices).map_err(|_| invalid("flux tangent allocation refused"))?;
        load.resize(self.vertices, 0.0);
        for (slot, face) in self.faces.iter().enumerate() {
            poll(cx, slot)?;
            let mut flux = 0.0;
            for (w, q) in face.weights.iter().zip(coordinates[face.control]) {
                flux = checked(flux + checked(w * q)?)?;
            }
            for &v in &face.vertices { load[v] = checked(load[v] - flux)?; }
        }
        for (slot, &v) in self.fixed.iter().enumerate() {
            if slot % 512 == 0 { poll(cx, slot)?; }
            load[v] = 0.0;
        }
        poll(cx, self.faces.len())?;
        Ok(load)
    }

    /// Contract the complete nodal-load multiplier with the original P1 traces.
    /// Fixed-node multipliers must already be zero. Output follows the declared
    /// control order and retains scalar versus world-vector coordinates.
    pub fn pullback(&self, cx: &Cx<'_>, nodal_load_adjoint: &[f64])
        -> Result<Vec<PrescribedFluxGradient>, ConductionError>
    {
        vector(cx, nodal_load_adjoint, self.vertices)?;
        for (slot, &v) in self.fixed.iter().enumerate() {
            if slot % 512 == 0 { poll(cx, slot)?; }
            if nodal_load_adjoint[v] != 0.0 {
                return Err(invalid("flux pullback requires zero prescribed-node multipliers"));
            }
        }
        let mut result = Vec::new();
        result.try_reserve_exact(self.vectors.len()).map_err(|_| invalid("flux gradient allocation refused"))?;
        for (slot, &vector) in self.vectors.iter().enumerate() {
            if slot % 512 == 0 { poll(cx, slot)?; }
            result.push(if vector { PrescribedFluxGradient::WorldVector([0.0; 3]) }
                else { PrescribedFluxGradient::Normal(0.0) });
        }
        for (slot, face) in self.faces.iter().enumerate() {
            poll(cx, slot)?;
            for &vertex in &face.vertices {
                let lambda = nodal_load_adjoint[vertex];
                match &mut result[face.control] {
                    PrescribedFluxGradient::Normal(value) => *value = checked(*value - checked(face.weights[0] * lambda)?)?,
                    PrescribedFluxGradient::WorldVector(values) => for (value, weight) in values.iter_mut().zip(face.weights) {
                        *value = checked(*value - checked(weight * lambda)?)?;
                    },
                }
            }
        }
        poll(cx, self.faces.len())?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
