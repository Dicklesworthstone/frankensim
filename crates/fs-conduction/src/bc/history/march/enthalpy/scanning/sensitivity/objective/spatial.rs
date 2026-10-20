//! Mesh-bound spatial readouts of the original nodal thermal fields.
//!
//! A readout is a nonnegative normalized sparse row. Interior P1 probes use
//! barycentric weights; region means use actual tetrahedral volumes; declared
//! sensor or reference-mass weights use `weighted_mean`. Temperature and phase
//! fraction are sampled BEFORE averaging, never inferred from mean enthalpy.
//! Geometry and weights are fixed; sensor dynamics and optical response are
//! not modeled. Volume means of mass fraction are NOT bulk mass fractions:
//! use reference-mass weights for that quantity.

use super::{ScanSensitivityError, finite, invalid, poll, zeros};
use crate::ConductionMesh;
use fs_exec::Cx;
use std::collections::BTreeMap;

/// One immutable spatial mean tied to the exact mesh used for construction.
/// Cloning copies only its bounded sparse row, not the physical mesh.
#[derive(Debug, Clone)]
pub struct ScanReadout<'m> {
    mesh: &'m ConductionMesh,
    weights: Vec<(usize, f64)>,
}

impl PartialEq for ScanReadout<'_> {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.mesh, other.mesh) && self.weights == other.weights
    }
}

impl<'m> ScanReadout<'m> {
    /// Exact original vertex readout, without interpolation or normalization drift.
    pub fn vertex(cx: &Cx<'_>, mesh: &'m ConductionMesh, vertex: usize)
        -> Result<Self, ScanSensitivityError>
    {
        Self::weighted_mean(cx, mesh, &[(vertex, 1.0)], 1)
    }

    /// A P1 point in one explicit tetrahedron. Nonnegative homogeneous
    /// barycentric weights are normalized; `[1;4]` means its centroid. Negative
    /// coordinates (extrapolation) and the all-zero tuple refuse. No nearest-node
    /// snapping or search across coincident material interfaces is performed.
    pub fn in_tetrahedron(cx: &Cx<'_>, mesh: &'m ConductionMesh, element: usize,
        barycentric: [f64; 4]) -> Result<Self, ScanSensitivityError>
    {
        poll(cx, 0)?;
        let tet = mesh.complex().tets.get(element)
            .ok_or_else(|| invalid("spatial readout names an absent tetrahedron"))?;
        let mut weights: Vec<_> = tet.iter().zip(barycentric)
            .map(|(&v, w)| (v as usize, w)).collect();
        weights.sort_by_key(|&(v, _)| v);
        Self::weighted_mean(cx, mesh, &weights, 4)
    }

    /// Volume mean of a P1 field over selected complete tetrahedra. This is the
    /// exact element integral of the represented linear field, evaluated in
    /// binary64, not a node-count average or an error bound on the continuum.
    /// Elements must be nonempty, strictly increasing, and fit `max_elements`.
    pub fn volume_mean(cx: &Cx<'_>, mesh: &'m ConductionMesh,
        elements: &[usize], max_elements: usize) -> Result<Self, ScanSensitivityError>
    {
        poll(cx, 0)?;
        if elements.is_empty() || elements.len() > max_elements {
            return Err(invalid("volume readout needs a nonempty bounded element selection"));
        }
        let mut weights = BTreeMap::<usize, f64>::new();
        let mut previous = None;
        for (i, &element) in elements.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            if element >= mesh.element_count() || previous.is_some_and(|old| old >= element) {
                return Err(invalid("volume readout elements must be distinct increasing IDs inside the mesh"));
            }
            previous = Some(element);
            let weight = finite(mesh.element_volume(element) / 4.0)?;
            if weight <= 0.0 { return Err(invalid("volume readout element measure is not representable")); }
            for &vertex in &mesh.complex().tets[element] {
                let entry = weights.entry(vertex as usize).or_default();
                *entry = finite(*entry + weight)?;
            }
        }
        let weights: Vec<_> = weights.into_iter().collect();
        Self::weighted_mean(cx, mesh, &weights, mesh.vertex_count())
    }

    /// A declared finite sensor footprint or reference-mass mean. Entries must
    /// have distinct increasing vertex IDs, finite nonnegative weights and at
    /// least one positive weight. Common weight units cancel on normalization.
    /// Zero entries are removed, not allowed to seed unrelated material corners.
    ///
    /// Admission precedes allocation. Normalization scales by the largest weight
    /// to avoid an overflowing raw sum. A positive coefficient that underflows
    /// to zero refuses instead of silently discarding a sensor contribution.
    pub fn weighted_mean(cx: &Cx<'_>, mesh: &'m ConductionMesh,
        weights: &[(usize, f64)], max_weights: usize) -> Result<Self, ScanSensitivityError>
    {
        poll(cx, 0)?;
        if weights.is_empty() || weights.len() > max_weights || weights.len() > mesh.vertex_count() {
            return Err(invalid("spatial readout needs a nonempty bounded sparse row"));
        }
        let mut maximum = 0.0_f64;
        let mut previous = None;
        for (i, &(vertex, weight)) in weights.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            if vertex >= mesh.vertex_count() || previous.is_some_and(|old| old >= vertex)
                || !weight.is_finite() || weight < 0.0
            { return Err(invalid("readout weights need distinct increasing mesh vertices and finite nonnegative values")); }
            previous = Some(vertex);
            maximum = maximum.max(weight);
        }
        if maximum == 0.0 { return Err(invalid("spatial readout has zero total weight")); }
        let mut sum = 0.0;
        for (i, &(_, weight)) in weights.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            let scaled = weight / maximum;
            if weight > 0.0 && scaled == 0.0 {
                return Err(invalid("positive readout weight underflows during normalization"));
            }
            sum = finite(sum + scaled)?;
        }
        let mut normalized = Vec::new();
        normalized.try_reserve_exact(weights.len())
            .map_err(|_| invalid("spatial readout allocation refused"))?;
        for (i, &(vertex, weight)) in weights.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            if weight == 0.0 { continue; }
            let coefficient = finite((weight / maximum) / sum)?;
            if coefficient == 0.0 { return Err(invalid("positive normalized readout coefficient underflows")); }
            normalized.push((vertex, coefficient));
        }
        poll(cx, weights.len())?;
        Ok(Self { mesh, weights: normalized })
    }

    /// Actual normalized coefficients in ascending original vertex order.
    pub fn weights(&self) -> &[(usize, f64)] { &self.weights }

    /// Immutable mesh identity. Equal vertex counts do not establish identity.
    pub fn mesh(&self) -> &'m ConductionMesh { self.mesh }

    /// Sample a field already expressed in its measured physical coordinate.
    /// Unsampled vertices are not read. A refusal returns no partial scalar.
    pub fn evaluate(&self, cx: &Cx<'_>, field: &[f64]) -> Result<f64, ScanSensitivityError> {
        poll(cx, 0)?;
        if field.len() != self.mesh.vertex_count() { return Err(invalid("readout field must match its exact mesh shape")); }
        let mut value = 0.0;
        for (i, &(vertex, weight)) in self.weights.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            value = finite(value + finite(weight * finite(field[vertex])?)?)?;
        }
        poll(cx, self.weights.len())?;
        Ok(value)
    }

    /// Transpose the SAME spatial row into the original nodal coordinates.
    /// Thermal-chart conversion is a separate step, after this scatter.
    pub fn pullback(&self, cx: &Cx<'_>, seed: f64) -> Result<Vec<f64>, ScanSensitivityError> {
        poll(cx, 0)?;
        finite(seed)?;
        let mut out = zeros(self.mesh.vertex_count())?;
        self.add_seed(cx, &mut out, seed)?;
        Ok(out)
    }

    fn add_seed(&self, cx: &Cx<'_>, out: &mut [f64], seed: f64) -> Result<(), ScanSensitivityError> {
        for (i, &(vertex, weight)) in self.weights.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            out[vertex] = finite(out[vertex] + finite(weight * seed)?)?;
        }
        poll(cx, self.weights.len())
    }
}
