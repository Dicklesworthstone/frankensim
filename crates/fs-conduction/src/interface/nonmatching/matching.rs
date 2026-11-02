//! Exact matching triangles using the existing pointwise contact stencil.
//! Geometry remains owned by the matching admission; no global plane,
//! intersection tolerance, nearest-point search or second contact operator.

use super::{Bound, Face, ResistanceField, Stencil, error, finite, poll};
use crate::interface::{InterfaceSurface, matching::matching_vertices};
use crate::{ConductionError, ConductionMesh, ThermalBoundary};
use fs_exec::Cx;

impl Bound {
    /// Called only after the complete matching owner has admitted all pairs,
    /// opposing traces and exclusive boundary ownership.
    pub(in crate::interface) fn matching_quadrature(
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        boundary: &ThermalBoundary,
        surface: &InterfaceSurface,
    ) -> Result<Self, ConductionError> {
        let name = surface.name();
        poll(cx, 0)?;
        let resistance = surface.resistance().ok_or_else(|| {
            error(
                name,
                "pointwise source binding requires one original reference law",
            )
        })?;
        let mut pairs = surface.face_pairs().to_vec();
        pairs.sort_by_key(|pair| (pair.side_a.min(pair.side_b), pair.side_a.max(pair.side_b)));
        let mut faces = Vec::new();
        let face_count = pairs
            .len()
            .checked_mul(2)
            .ok_or_else(|| error(name, "matching quadrature face count overflows"))?;
        faces
            .try_reserve_exact(face_count)
            .map_err(|_| error(name, "matching quadrature face allocation refused"))?;
        let mut stencils = Vec::new();
        stencils
            .try_reserve_exact(pairs.len())
            .map_err(|_| error(name, "matching quadrature stencil allocation refused"))?;
        for (index, pair) in pairs.into_iter().enumerate() {
            poll(cx, index)?;
            let a = Face::read(mesh, boundary, pair.side_a, name)?;
            let b = Face::read(mesh, boundary, pair.side_b, name)?;
            let mapping = matching_vertices(mesh, pair)?;
            // Canonical physical vertex order makes the integration points
            // independent of which original trace is named A or B.
            let mut order = [0, 1, 2];
            order.sort_by(|&i, &j| {
                a.points[i][0]
                    .total_cmp(&a.points[j][0])
                    .then(a.points[i][1].total_cmp(&a.points[j][1]))
                    .then(a.points[i][2].total_cmp(&a.points[j][2]))
            });
            let va = order.map(|i| a.vertices[i]);
            let vb = order.map(|i| b.vertices[mapping[i]]);
            let mut basis = [[0.0; 6]; 3];
            for (q, row) in basis.iter_mut().enumerate() {
                for i in 0..3 {
                    row[i] = if i == q { 2.0 / 3.0 } else { 1.0 / 6.0 };
                    row[i + 3] = -row[i];
                }
            }
            let weight_m2 = finite(a.area / 3.0, name)?;
            let weighted_conductance = finite(weight_m2 / resistance.value_m2_k_per_w(), name)?;
            if weight_m2 <= 0.0 || weighted_conductance <= 0.0 {
                return Err(error(
                    name,
                    "matching quadrature weight or conductance underflow",
                ));
            }
            stencils.push(Stencil {
                pair,
                vertices: [va[0], va[1], va[2], vb[0], vb[1], vb[2]],
                resistance_index: 0,
                basis,
                weight_m2,
                delegated: true,
            });
            faces.extend([a, b]);
        }
        Ok(Self {
            name: name.into(),
            matching_only: true,
            resistances: ResistanceField::Uniform(resistance.clone()),
            faces: faces.into(),
            stencils: stencils.into(),
            quadrature_resistances: None,
        })
    }

    /// Guard the actual paired P1 corner coordinates, not an invented
    /// Cartesian product of unrelated extrema from the two sides. These are
    /// support checks, not extra integration points or averaged source claims.
    pub(in crate::interface) fn matching_source_support(
        &self,
        cx: &Cx<'_>,
        temperature: &[f64],
        mut query: impl FnMut([f64; 2]) -> Result<(), ConductionError>,
    ) -> Result<(), ConductionError> {
        self.matching_state_support(cx,temperature,None,|_,traces,_| query(traces))
    }

    /// Admit each actual paired corner with its companion P1 source coordinate.
    pub(in crate::interface) fn matching_state_support(
        &self, cx: &Cx<'_>, temperature: &[f64], companion: Option<&[f64]>,
        mut query: impl FnMut([usize; 2], [f64; 2], Option<[f64; 2]>) -> Result<(), ConductionError>,
    ) -> Result<(), ConductionError> {
        if !self.matching_only {
            return Ok(());
        }
        self.check_field(temperature)?;
        if let Some(field) = companion { self.check_field(field)?; }
        for (index, stencil) in self.stencils.iter().enumerate() {
            poll(cx, index)?;
            for i in 0..3 {
                query([stencil.vertices[i], stencil.vertices[i + 3]], [
                    temperature[stencil.vertices[i]],
                    temperature[stencil.vertices[i + 3]],
                ], companion.map(|field| [field[stencil.vertices[i]], field[stencil.vertices[i + 3]]]))?;
            }
        }
        Ok(())
    }
}
