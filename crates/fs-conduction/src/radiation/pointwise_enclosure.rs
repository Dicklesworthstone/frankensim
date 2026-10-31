//! Local P1 emission with reflected exchange through an admitted patch system.
//!
//! Each patch has uniform irradiation G_i, but the surface flux is
//! epsilon_i * (sigma T_h(x)^4 - G_i). The radiosity RHS uses the integral
//! of T_h^4, never the fourth power of a mean temperature. Geometry, emissivity
//! and view factors stay fixed. This resolves local emission and absorption
//! under the supplied patch irradiation model, not a spatial visibility kernel.

use fs_blake3::ContentHash;
use fs_exec::Cx;
use std::{collections::BTreeMap, sync::Arc};

use super::{
    GrayDiffuseEnclosure, OpenDiffuseEnclosure, PreparedRadiosity, RadiationSurface,
    ViewFactorMatrix, STEFAN_BOLTZMANN_W_M2_K4 as SIGMA,
    SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, balance, pointwise::TRIANGLE_DEGREE_FIVE,
    radiation_error, pointwise::PointwiseRadiationParameterGradient,
};
use crate::{ConductionError, ConductionMesh, ThermalBc, ThermalBoundary, assemble::ASSEMBLY_TILE};

/// Independent watt gates for the radiosity equations and closed heat sum.
#[derive(Debug, Clone, Copy)]
pub struct PointwiseEnclosureTolerance {
    pub absolute_heat_w: f64,
    pub relative_heat: f64,
}

impl Default for PointwiseEnclosureTolerance {
    fn default() -> Self {
        Self {
            absolute_heat_w: 1e-9,
            relative_heat: 1e-10,
        }
    }
}

/// Source declaration and physical integrals on one radiating patch.
#[derive(Debug, Clone, PartialEq)]
pub struct PointwiseEnclosurePatchReport {
    pub name: String,
    pub emissivity: SurfaceEmissivity,
    pub area_m2: f64,
    pub mean_temperature_k: f64,
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    /// Area mean of sigma T_h^4, W/m2; not sigma mean(T_h)^4.
    pub mean_blackbody_emission_w_m2: f64,
    pub radiosity_w_m2: f64,
    pub irradiation_w_m2: f64,
    pub outward_heat_w: f64,
}

/// Explicit black bath, independently evaluated from its own reciprocal row.
/// This is an external reservoir, never a solid thermal degree of freedom.
#[derive(Debug, Clone, PartialEq)]
pub struct PointwiseBlackBathReport {
    pub temperature_k: f64,
    /// Sum of physical surface area times its ambient view factor, m².
    pub effective_area_m2: f64,
    pub radiosity_w_m2: f64,
    pub irradiation_w_m2: f64,
    /// Positive heat leaves the bath; negative heat is absorbed by the bath.
    pub outward_heat_w: f64,
}

/// Radiation accounting, separate from external convection/prescribed flux.
#[derive(Debug, Clone, PartialEq)]
pub struct PointwiseEnclosureReport {
    pub patches: Vec<PointwiseEnclosurePatchReport>,
    /// Actual signed PHYSICAL SURFACE sum, included once in solid energy.
    /// Open cavities generally lose or gain net heat through the black bath.
    pub outward_heat_w: f64,
    pub absolute_patch_heat_w: f64,
    pub linear_residual_max_w_m2: f64,
    pub view_factor_identity: ContentHash,
    /// Absent for a closed enclosure; otherwise the single external black bath.
    pub ambient: Option<PointwiseBlackBathReport>,
    /// Sum of surface and bath powers; neither term is forced to close.
    pub system_energy_closure_w: f64,
}

#[derive(Debug)]
struct PatchGeometry {
    vertices: Vec<usize>,
    /// Exact P1 integral of N_i, m2, in local vertex order.
    weights: Vec<f64>,
    emissivity: f64,
}

#[derive(Debug)]
struct FaceGeometry {
    vertices: [usize; 3],
    local: [usize; 3],
    patch: usize,
    area: f64,
}

#[derive(Debug)]
struct Geometry {
    vertices: usize,
    patches: Vec<PatchGeometry>,
    faces: Vec<FaceGeometry>,
    /// Sparse memberships permit exact diagonal assembly at shared edges.
    memberships: Vec<(usize, Vec<(usize, usize)>)>,
    /// F [I-(1-epsilon)F]^-1 diag(epsilon), prepared once.
    irradiation_response: Vec<Vec<f64>>,
    /// Same reflected response to unit BLACK BATH emission, when present.
    ambient_irradiation_response: Option<Vec<f64>>,
}

#[derive(Debug, Clone)]
struct FaceTangent {
    vertices: [usize; 3],
    patch: usize,
    matrix: [[f64; 3]; 3],
    emission_delta_load: [f64; 3],
}

/// Local blocks plus reflected patch feedback. This tangent is generally
/// NONSYMMETRIC: the source column carries the local T_h^3 derivative.
#[derive(Debug, Clone)]
pub struct PointwiseEnclosureAssembly {
    pub outward_load_w: Vec<f64>,
    /// Exact diagonal of the complete reflected temperature tangent. This can
    /// be negative on a nonuniform patch; it is not a positive preconditioner.
    pub temperature_diagonal_w_k: Vec<f64>,
    pub report: PointwiseEnclosureReport,
    geometry: Arc<Geometry>,
    local: Vec<FaceTangent>,
    emission_slopes: Vec<Vec<f64>>,
    emission_diagonal: Vec<f64>,
    irradiation_corrections: Vec<f64>,
    emission_minus_irradiation: Vec<f64>,
}

impl PointwiseEnclosureAssembly {
    pub(crate) fn emission_diagonal(&self) -> &[f64] {
        &self.emission_diagonal
    }

    /// Replace output by the complete nonsymmetric temperature tangent action.
    pub fn tangent_apply(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
        output: &mut [f64],
    ) -> Result<(), ConductionError> {
        self.check_fields(direction.len(), output.len())?;
        output.fill(0.0);
        self.add_scaled_tangent(cx, direction, None, 1.0, |_| true, output)
    }

    pub(crate) fn add_scaled_tangent(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
        column_scale: Option<&[f64]>,
        scale: f64,
        include_row: impl Fn(usize) -> bool,
        output: &mut [f64],
    ) -> Result<(), ConductionError> {
        poll(cx, 0)?;
        self.check_fields(direction.len(), output.len())?;
        if column_scale.is_some_and(|s| s.len() != direction.len()) {
            return Err(invalid(
                "enclosure tangent column scaling has the wrong length",
            ));
        }
        let value = |vertex| finite(direction[vertex] * column_scale.map_or(1.0, |s| s[vertex]));
        for (slot, face) in self.local.iter().enumerate() {
            if slot % ASSEMBLY_TILE == 0 {
                poll(cx, slot)?;
            }
            let local = [
                value(face.vertices[0])?,
                value(face.vertices[1])?,
                value(face.vertices[2])?,
            ];
            for (i, &vertex) in face.vertices.iter().enumerate() {
                if include_row(vertex) {
                    let applied = face.matrix[i][0].mul_add(
                        local[0],
                        face.matrix[i][1].mul_add(local[1], face.matrix[i][2] * local[2]),
                    );
                    output[vertex] = finite(scale.mul_add(applied, output[vertex]))?;
                }
            }
        }
        let mut emission = vec![0.0; self.geometry.patches.len()];
        for (patch, (geometry, slopes)) in self
            .geometry
            .patches
            .iter()
            .zip(&self.emission_slopes)
            .enumerate()
        {
            poll(cx, patch)?;
            for (slot, (&vertex, &slope)) in geometry.vertices.iter().zip(slopes).enumerate() {
                if slot % ASSEMBLY_TILE == 0 {
                    poll(cx, slot)?;
                }
                emission[patch] = finite(slope.mul_add(value(vertex)?, emission[patch]))?;
            }
        }
        for (patch, geometry) in self.geometry.patches.iter().enumerate() {
            poll(cx, patch)?;
            let irradiation = dot(&self.geometry.irradiation_response[patch], &emission)?;
            let coefficient = finite(-scale * geometry.emissivity * irradiation)?;
            for (slot, (&vertex, &weight)) in
                geometry.vertices.iter().zip(&geometry.weights).enumerate()
            {
                if slot % ASSEMBLY_TILE == 0 {
                    poll(cx, slot)?;
                }
                if include_row(vertex) {
                    output[vertex] = finite(coefficient.mul_add(weight, output[vertex]))?;
                }
            }
        }
        poll(cx, self.local.len())
    }

    /// Apply the genuine transpose, including nonlocal reflected irradiation.
    pub fn tangent_apply_transpose(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
        output: &mut [f64],
    ) -> Result<(), ConductionError> {
        self.check_fields(direction.len(), output.len())?;
        output.fill(0.0);
        self.add_scaled_tangent_transpose(cx, direction, None, 1.0, |_| true, output)
    }

    // Transpose of scale * P * (D - W diag(epsilon) B S) * diag(dT/dh).
    // The prescribed-row mask acts on input, never on the resulting columns.
    pub(crate) fn add_scaled_tangent_transpose(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
        column_scale: Option<&[f64]>,
        scale: f64,
        include_row: impl Fn(usize) -> bool,
        output: &mut [f64],
    ) -> Result<(), ConductionError> {
        poll(cx, 0)?;
        self.check_fields(direction.len(), output.len())?;
        if column_scale.is_some_and(|s| s.len() != direction.len()) {
            return Err(invalid("enclosure transpose scaling has the wrong length"));
        }
        let value = |vertex| {
            if include_row(vertex) { finite(direction[vertex]) } else { Ok(0.0) }
        };
        for (slot, face) in self.local.iter().enumerate() {
            if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
            let local = [
                value(face.vertices[0])?,
                value(face.vertices[1])?,
                value(face.vertices[2])?,
            ];
            for (j, &vertex) in face.vertices.iter().enumerate() {
                let applied = finite(face.matrix[0][j].mul_add(
                    local[0],
                    face.matrix[1][j].mul_add(local[1], face.matrix[2][j] * local[2]),
                ))?;
                let coefficient = finite(scale * column_scale.map_or(1.0, |s| s[vertex]))?;
                output[vertex] = finite(coefficient.mul_add(applied, output[vertex]))?;
            }
        }
        let p = self.geometry.patches.len();
        let mut irradiation_bar = vec![0.0; p];
        for (patch, geometry) in self.geometry.patches.iter().enumerate() {
            poll(cx, patch)?;
            let mut integral = 0.0;
            for (slot, (&vertex, &weight)) in
                geometry.vertices.iter().zip(&geometry.weights).enumerate()
            {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                integral = finite(weight.mul_add(value(vertex)?, integral))?;
            }
            irradiation_bar[patch] = finite(geometry.emissivity * integral)?;
        }
        for (column, (geometry, slopes)) in self.geometry.patches.iter()
            .zip(&self.emission_slopes).enumerate()
        {
            poll(cx, column)?;
            let mut emission_bar = 0.0;
            for (row, &bar) in irradiation_bar.iter().enumerate() {
                emission_bar = finite(
                    self.geometry.irradiation_response[row][column]
                        .mul_add(bar, emission_bar),
                )?;
            }
            for (slot, (&vertex, &slope)) in geometry.vertices.iter().zip(slopes).enumerate() {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                let coefficient = finite(
                    -scale * column_scale.map_or(1.0, |s| s[vertex]) * emission_bar,
                )?;
                output[vertex] = finite(coefficient.mul_add(slope, output[vertex]))?;
            }
        }
        poll(cx, self.local.len())
    }

    /// Contract absolute emissivity derivatives of the actual spatial load.
    /// The prepared irradiation response includes every reflected path, so
    /// this contraction needs no additional radiosity or thermal solve.
    pub fn parameter_pullback(
        &self,
        cx: &Cx<'_>,
        outward_load_weights: &[f64],
    ) -> Result<PointwiseRadiationParameterGradient, ConductionError> {
        self.check_fields(outward_load_weights.len(), outward_load_weights.len())?;
        self.admit_adjoint(cx)?;
        for (vertex, &weight) in outward_load_weights.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 { poll(cx, vertex)?; }
            finite(weight)?;
        }
        let p = self.geometry.patches.len();
        let mut area_weights = vec![0.0; p];
        let mut emissivities = vec![0.0; p];
        for (patch, geometry) in self.geometry.patches.iter().enumerate() {
            poll(cx, patch)?;
            for (slot, (&vertex, &area)) in
                geometry.vertices.iter().zip(&geometry.weights).enumerate()
            {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                area_weights[patch] = finite(
                    area.mul_add(outward_load_weights[vertex], area_weights[patch]),
                )?;
            }
            emissivities[patch] = finite(
                -self.irradiation_corrections[patch] * area_weights[patch],
            )?;
        }
        // Retain the spatial integral of v_h sigma(T_h^4-reference^4).
        // A patch mean heat multiplied by a mean cotangent loses this term.
        for (slot, face) in self.local.iter().enumerate() {
            if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
            for (i, &vertex) in face.vertices.iter().enumerate() {
                emissivities[face.patch] = finite(
                    face.emission_delta_load[i]
                        .mul_add(outward_load_weights[vertex], emissivities[face.patch]),
                )?;
            }
        }
        for (column, derivative) in emissivities.iter_mut().enumerate() {
            poll(cx, column)?;
            let epsilon = self.geometry.patches[column].emissivity;
            let mut irradiation_bar = 0.0;
            for (row, &area_weight) in area_weights.iter().enumerate() {
                // B = F A^-1 diag(epsilon); divide its column before the
                // contraction to avoid forming (E-G)/epsilon at small epsilon.
                let response = finite(
                    self.geometry.irradiation_response[row][column] / epsilon,
                )?;
                irradiation_bar = finite(
                    response.mul_add(
                        finite(self.geometry.patches[row].emissivity * area_weight)?,
                        irradiation_bar,
                    ),
                )?;
            }
            *derivative = finite(
                (-self.emission_minus_irradiation[column])
                    .mul_add(irradiation_bar, *derivative),
            )?;
        }
        let reservoir_temperatures = if let Some(ambient) = &self.report.ambient {
            let response = self.geometry.ambient_irradiation_response.as_ref()
                .ok_or_else(|| invalid("open radiation lost its prepared ambient response"))?;
            let t = ambient.temperature_k;
            let derivative = finite(4.0 * SIGMA * t * t * t)?;
            let mut contraction = 0.0;
            for (patch, &weight) in area_weights.iter().enumerate() {
                let coefficient = finite(-self.geometry.patches[patch].emissivity
                    * response[patch] * derivative)?;
                contraction = finite(coefficient.mul_add(weight, contraction))?;
            }
            vec![contraction]
        } else { Vec::new() };
        poll(cx, p)?;
        Ok(PointwiseRadiationParameterGradient {
            radiation_regions: self.report.patches.iter().map(|p| p.name.clone()).collect(),
            emissivities,
            reservoir_temperatures,
        })
    }

    pub(crate) fn admit_adjoint(&self, cx: &Cx<'_>) -> Result<(), ConductionError> {
        for (patch, report) in self.report.patches.iter().enumerate() {
            poll(cx, patch)?;
            if let Some((low, high)) = report.emissivity.temperature_validity_k() {
                for value in [report.minimum_temperature_k, report.maximum_temperature_k] {
                    if value <= low || value >= high {
                        return Err(ConductionError::MaterialQuery {
                            property: SURFACE_EMISSIVITY_PROPERTY.into(),
                            temperature: value,
                            upstream: format!(
                                "surface {:?} adjoint needs every node strictly inside retained emissivity support ({low}, {high}) K",
                                report.name,
                            ),
                        });
                    }
                }
            }
        }
        Ok(())
    }

    fn check_fields(&self, direction: usize, output: usize) -> Result<(), ConductionError> {
        if direction != self.geometry.vertices || output != self.geometry.vertices {
            return Err(invalid("enclosure tangent requires full nodal fields"));
        }
        Ok(())
    }
}

#[derive(Debug)]
enum PatchSystem<'a> {
    Closed { enclosure: &'a GrayDiffuseEnclosure, factors: PreparedRadiosity<'a> },
    // Own an immutable snapshot so a caller may lower a local open declaration
    // and return this spatial owner. Bath-temperature changes clone the Arc.
    Open(Arc<OpenDiffuseEnclosure>),
}

impl PatchSystem<'_> {
    fn surfaces(&self) -> &[RadiationSurface] {
        match self { Self::Closed { enclosure, .. } => enclosure.surfaces(),
            Self::Open(enclosure) => enclosure.surfaces() }
    }

    fn view_factors(&self) -> &ViewFactorMatrix {
        match self { Self::Closed { enclosure, .. } => enclosure.view_factors(),
            Self::Open(enclosure) => enclosure.view_factors() }
    }

    fn solve_rhs(&self, cx: &Cx<'_>, rhs: &[f64]) -> Result<Vec<f64>, ConductionError> {
        match self { Self::Closed { factors, .. } => factors.solve_rhs(cx, rhs),
            Self::Open(enclosure) => enclosure.solve_rhs(cx, rhs) }
    }
}

/// Prepared fixed-geometry enclosure, with sparse trace geometry and one small
/// radiosity factorization. At most 64 patches; trace storage is linear in
/// boundary faces/vertices and the reflected matrix is O(patches^2). An open
/// owner includes one external black bath without any corresponding solid node.
#[derive(Debug)]
pub struct PointwiseEnclosureRadiation<'a> {
    mesh: &'a ConductionMesh,
    boundary: &'a ThermalBoundary,
    system: PatchSystem<'a>,
    ambient_temperature_k: Option<f64>,
    geometry: Arc<Geometry>,
    tolerance: PointwiseEnclosureTolerance,
}

impl<'a> PointwiseEnclosureRadiation<'a> {
    pub fn new(
        cx: &Cx<'_>,
        mesh: &'a ConductionMesh,
        boundary: &'a ThermalBoundary,
        enclosure: &'a GrayDiffuseEnclosure,
        tolerance: PointwiseEnclosureTolerance,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        admit_configuration(enclosure.surfaces().len(), tolerance)?;
        let factors = enclosure.prepare_radiosity(cx)?;
        Self::prepare(cx, mesh, boundary,
            PatchSystem::Closed { enclosure, factors }, None, tolerance)
    }

    /// Integrate actual P1 emission with reflected exchange and one explicitly
    /// black ambient reservoir. The completed open matrix remains the original
    /// radiosity owner; a bath is never assigned a fictitious solid node.
    ///
    /// Retains an owned immutable snapshot of the supplied model, cloning its
    /// prepared factors without another factorization. The returned lifetime
    /// depends only on mesh/boundary, so a locally lowered model is supported.
    /// All existing pointwise enthalpy, adaptive and accepted-adjoint consumers
    /// use this same load and complete tangent. Storage supplies a transient
    /// thermal anchor even when all original boundaries are adiabatic.
    pub fn new_open(
        cx: &Cx<'_>,
        mesh: &'a ConductionMesh,
        boundary: &'a ThermalBoundary,
        enclosure: &OpenDiffuseEnclosure,
        ambient_temperature_k: f64,
        tolerance: PointwiseEnclosureTolerance,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        admit_configuration(enclosure.surfaces().len(), tolerance)?;
        super::require_temperature("open-ambient", "black ambient temperature", ambient_temperature_k)?;
        Self::prepare(cx, mesh, boundary, PatchSystem::Open(Arc::new(enclosure.clone())),
            Some(ambient_temperature_k), tolerance)
    }

    /// Change a prescribed black-bath temperature without rebuilding geometry,
    /// refactoring radiosity or advancing thermal history. Closed owners refuse.
    /// The caller chooses the endpoint bath temperature for each time interval.
    pub fn with_ambient_temperature_k(&self, ambient_temperature_k: f64)
        -> Result<Self, ConductionError>
    {
        super::require_temperature("open-ambient", "black ambient temperature", ambient_temperature_k)?;
        let PatchSystem::Open(enclosure) = &self.system else {
            return Err(invalid("a closed enclosure has no ambient temperature to replace"));
        };
        Ok(Self { mesh: self.mesh, boundary: self.boundary,
            system: PatchSystem::Open(Arc::clone(enclosure)),
            ambient_temperature_k: Some(ambient_temperature_k),
            geometry: Arc::clone(&self.geometry), tolerance: self.tolerance })
    }

    fn prepare(
        cx: &Cx<'_>,
        mesh: &'a ConductionMesh,
        boundary: &'a ThermalBoundary,
        system: PatchSystem<'a>,
        ambient_temperature_k: Option<f64>,
        tolerance: PointwiseEnclosureTolerance,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        let count = system.surfaces().len();
        admit_configuration(count, tolerance)?;
        mesh.boundary()
            .len()
            .checked_mul(std::mem::size_of::<FaceGeometry>() + std::mem::size_of::<FaceTangent>())
            .filter(|&bytes| bytes <= isize::MAX as usize)
            .ok_or_else(|| invalid("pointwise enclosure face storage overflows"))?;
        let identity = super::radiation_mesh_identity_with_poll(mesh, |at| poll(cx, at))?;
        let mut owners = vec![None; mesh.boundary().len()];
        let mut region_face_counts = vec![0usize; boundary.region_names().len()];
        for slot in 0..mesh.boundary().len() {
            if slot % ASSEMBLY_TILE == 0 {
                poll(cx, slot)?;
            }
            if let Some(region) = boundary.region_for(slot) {
                region_face_counts[region] += 1;
            }
        }
        let mut patches = Vec::with_capacity(count);
        let mut memberships: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
        for (patch, surface) in system.surfaces().iter().enumerate() {
            poll(cx, patch)?;
            if surface.mesh_identity != identity {
                return Err(invalid("pointwise enclosure belongs to a different mesh"));
            }
            let region = boundary
                .region_names()
                .iter()
                .position(|name| name == surface.name())
                .ok_or_else(|| invalid("pointwise enclosure requires existing named regions"))?;
            if matches!(boundary.conditions()[region], ThermalBc::Dirichlet { .. }) {
                return Err(invalid(
                    "pointwise enclosure requires natural boundary traces",
                ));
            }
            boundary.conditions()[region].validate(mesh.vertex_count())?;
            let mut weights = BTreeMap::<usize, f64>::new();
            for (ordinal, &slot) in surface.face_slots().iter().enumerate() {
                if ordinal % ASSEMBLY_TILE == 0 {
                    poll(cx, ordinal)?;
                }
                if boundary.region_for(slot) != Some(region)
                    || owners[slot].replace(patch).is_some()
                {
                    return Err(invalid(
                        "enclosure trace overlaps or differs from its named region",
                    ));
                }
                for &vertex in &mesh.boundary()[slot].vertices {
                    let weight = weights.entry(vertex as usize).or_default();
                    *weight = finite(*weight + mesh.boundary()[slot].area / 3.0)?;
                }
            }
            if region_face_counts[region] != surface.face_slots().len() {
                return Err(invalid("enclosure trace must cover its whole named region"));
            }
            let (vertices, weights): (Vec<_>, Vec<_>) = weights.into_iter().unzip();
            for (local, &vertex) in vertices.iter().enumerate() {
                memberships.entry(vertex).or_default().push((patch, local));
            }
            patches.push(PatchGeometry {
                vertices,
                weights,
                emissivity: surface.emissivity().value(),
            });
        }
        let mut faces = Vec::new();
        for (slot, owner) in owners.into_iter().enumerate() {
            if slot % ASSEMBLY_TILE == 0 {
                poll(cx, slot)?;
            }
            if let Some(patch) = owner {
                let face = &mesh.boundary()[slot];
                let vertices = face.vertices.map(|v| v as usize);
                let local = vertices.map(|v| {
                    patches[patch]
                        .vertices
                        .binary_search(&v)
                        .expect("prepared trace vertex")
                });
                faces.push(FaceGeometry {
                    vertices,
                    local,
                    patch,
                    area: face.area,
                });
            }
        }
        let mut irradiation_response = vec![vec![0.0; count]; count];
        let mut ambient_irradiation_response = ambient_temperature_k.map(|_| vec![0.0; count]);
        let equations = system.view_factors().len();
        for column in 0..equations {
            poll(cx, column)?;
            let mut rhs = vec![0.0; equations];
            rhs[column] = if column < count { patches[column].emissivity } else { 1.0 };
            let response = system.solve_rhs(cx, &rhs)?;
            for (row, values) in irradiation_response.iter_mut().enumerate() {
                let value = dot(&system.view_factors().factors()[row], &response)?;
                if column < count { values[column] = value; }
                else if let Some(ambient) = &mut ambient_irradiation_response { ambient[row] = value; }
            }
        }
        Ok(Self {
            mesh,
            boundary,
            system,
            ambient_temperature_k,
            tolerance,
            geometry: Arc::new(Geometry {
                vertices: mesh.vertex_count(),
                patches,
                faces,
                memberships: memberships.into_iter().collect(),
                irradiation_response,
                ambient_irradiation_response,
            }),
        })
    }

    pub(crate) fn patch_count(&self) -> usize { self.geometry.patches.len() }

    pub(crate) fn matches(&self, mesh: &ConductionMesh, boundary: &ThermalBoundary) -> bool {
        std::ptr::eq(self.mesh, mesh) && std::ptr::eq(self.boundary, boundary)
    }

    /// Integrate local emission and solve reflected irradiation, with the exact
    /// temperature tangent. Use the nonsymmetric Krylov lane for this operator.
    pub fn assemble(
        &self,
        cx: &Cx<'_>,
        temperature: &[f64],
    ) -> Result<PointwiseEnclosureAssembly, ConductionError> {
        self.integrate::<true>(cx, temperature)
    }

    pub(crate) fn load_and_report(
        &self,
        cx: &Cx<'_>,
        temperature: &[f64],
    ) -> Result<(Vec<f64>, PointwiseEnclosureReport), ConductionError> {
        let assembly = self.integrate::<false>(cx, temperature)?;
        Ok((assembly.outward_load_w, assembly.report))
    }

    fn integrate<const TANGENT: bool>(
        &self,
        cx: &Cx<'_>,
        temperature: &[f64],
    ) -> Result<PointwiseEnclosureAssembly, ConductionError> {
        poll(cx, 0)?;
        let n = self.geometry.vertices;
        let p = self.geometry.patches.len();
        if temperature.len() != n {
            return Err(ConductionError::FieldLength {
                field: "pointwise enclosure temperature",
                expected: n,
                found: temperature.len(),
            });
        }
        let mut minimum = vec![f64::INFINITY; p];
        let mut maximum = vec![f64::NEG_INFINITY; p];
        for (patch, geometry) in self.geometry.patches.iter().enumerate() {
            let surface = &self.system.surfaces()[patch];
            for (slot, &vertex) in geometry.vertices.iter().enumerate() {
                if slot % ASSEMBLY_TILE == 0 {
                    poll(cx, slot)?;
                }
                let t = temperature[vertex];
                super::require_temperature(surface.name(), "radiating surface node", t)?;
                if let Some((lo, hi)) = surface.emissivity().temperature_validity_k()
                    && !(lo..=hi).contains(&t)
                {
                    return Err(ConductionError::MaterialQuery {
                        property: SURFACE_EMISSIVITY_PROPERTY.into(),
                        temperature: t,
                        upstream: format!(
                            "surface {:?} node is outside emissivity support [{lo}, {hi}] K",
                            surface.name()
                        ),
                    });
                }
                minimum[patch] = minimum[patch].min(t);
                maximum[patch] = maximum[patch].max(t);
            }
        }
        let surface_reference = minimum.iter().copied().fold(f64::INFINITY, f64::min);
        let reference = self.ambient_temperature_k
            .map_or(surface_reference, |ambient| ambient.min(surface_reference));
        let background = finite(SIGMA * reference * reference * reference * reference)?;
        let mut emission_delta = vec![0.0; p];
        let mut load = vec![0.0; n];
        let mut diagonal = if TANGENT { vec![0.0; n] } else { Vec::new() };
        let mut local = if TANGENT {
            Vec::with_capacity(self.geometry.faces.len())
        } else {
            Vec::new()
        };
        let mut slopes: Vec<Vec<f64>> = if TANGENT {
            self.geometry
                .patches
                .iter()
                .map(|g| vec![0.0; g.vertices.len()])
                .collect()
        } else {
            Vec::new()
        };
        for (slot, face) in self.geometry.faces.iter().enumerate() {
            if slot % ASSEMBLY_TILE == 0 {
                poll(cx, slot)?;
            }
            let t = face.vertices.map(|vertex| temperature[vertex]);
            let epsilon = self.geometry.patches[face.patch].emissivity;
            let mut matrix = [[0.0; 3]; 3];
            let mut emission_delta_load = [0.0; 3];
            for (bary, weight) in TRIANGLE_DEGREE_FIVE {
                let tq = finite(t[0] + (t[1] - t[0]).mul_add(bary[1], (t[2] - t[0]) * bary[2]))?;
                let area_weight = finite(face.area * weight)?;
                let emission =
                    finite(area_weight * balance::emission_difference(SIGMA, tq, reference))?;
                emission_delta[face.patch] = finite(emission_delta[face.patch] + emission)?;
                for i in 0..3 {
                    load[face.vertices[i]] =
                        finite(load[face.vertices[i]] + epsilon * bary[i] * emission)?;
                }
                if TANGENT {
                    let derivative = finite(area_weight * 4.0 * SIGMA * tq * tq * tq)?;
                    for i in 0..3 {
                        emission_delta_load[i] = finite(
                            bary[i].mul_add(emission, emission_delta_load[i]),
                        )?;
                        let value = &mut slopes[face.patch][face.local[i]];
                        *value = finite(*value + derivative * bary[i])?;
                        for j in 0..=i {
                            matrix[i][j] =
                                finite(matrix[i][j] + epsilon * derivative * bary[i] * bary[j])?;
                        }
                    }
                }
            }
            if TANGENT {
                for i in 0..3 {
                    diagonal[face.vertices[i]] = finite(diagonal[face.vertices[i]] + matrix[i][i])?;
                    for j in 0..i {
                        matrix[j][i] = matrix[i][j];
                    }
                }
                local.push(FaceTangent {
                    vertices: face.vertices,
                    patch: face.patch,
                    matrix,
                    emission_delta_load,
                });
            }
        }
        let factors = self.system.view_factors().factors();
        let mut row_sums = Vec::with_capacity(p);
        let mut rhs = Vec::with_capacity(factors.len());
        for patch in 0..p {
            poll(cx, patch)?;
            let area = self.system.surfaces()[patch].area_m2();
            emission_delta[patch] = finite(emission_delta[patch] / area)?;
            if TANGENT {
                for value in &mut slopes[patch] {
                    *value = finite(*value / area)?;
                }
            }
            let sum = finite(factors[patch].iter().sum())?;
            row_sums.push(sum);
            rhs.push(finite(balance::shifted_load(
                self.geometry.patches[patch].emissivity,
                emission_delta[patch],
                background,
                sum,
            ))?);
        }
        // A black bath's radiosity is its Planck integral independently of
        // incoming light. It contributes a source equation, never a solid row.
        if let Some(ambient) = self.ambient_temperature_k {
            rhs.push(finite(balance::emission_difference(SIGMA, ambient, reference))?);
        }
        let radiosity_delta = self.system.solve_rhs(cx, &rhs)?;
        let mut reports = Vec::with_capacity(p);
        let mut irradiation_corrections = Vec::with_capacity(if TANGENT { p } else { 0 });
        let mut emission_minus_irradiation = Vec::with_capacity(if TANGENT { p } else { 0 });
        let mut total = 0.0;
        let mut absolute = 0.0;
        let mut residual_max = 0.0_f64;
        let mut max_area = 0.0_f64;
        for (patch, geometry) in self.geometry.patches.iter().enumerate() {
            poll(cx, patch)?;
            let surface = &self.system.surfaces()[patch];
            let irradiation_delta = dot(&factors[patch], &radiosity_delta)?;
            let residual = finite(
                radiosity_delta[patch]
                    - (1.0 - geometry.emissivity) * irradiation_delta
                    - rhs[patch],
            )?;
            residual_max = residual_max.max(residual.abs());
            max_area = max_area.max(surface.area_m2());
            let correction = finite(background.mul_add(row_sums[patch] - 1.0, irradiation_delta))?;
            if TANGENT {
                irradiation_corrections.push(correction);
                emission_minus_irradiation.push(finite(emission_delta[patch] - correction)?);
            }
            let mut mean_delta = 0.0;
            for (slot, (&vertex, &weight)) in
                geometry.vertices.iter().zip(&geometry.weights).enumerate()
            {
                if slot % ASSEMBLY_TILE == 0 {
                    poll(cx, slot)?;
                }
                load[vertex] = finite(load[vertex] - geometry.emissivity * weight * correction)?;
                mean_delta = finite(weight.mul_add(temperature[vertex] - reference, mean_delta))?;
            }
            let heat = finite(
                surface.area_m2()
                    * balance::outward_flux(
                        geometry.emissivity,
                        emission_delta[patch],
                        irradiation_delta,
                        background,
                        row_sums[patch],
                    ),
            )?;
            total = finite(total + heat)?;
            absolute = finite(absolute + heat.abs())?;
            let radiosity = finite(background + radiosity_delta[patch])?;
            let irradiation = finite(background.mul_add(row_sums[patch], irradiation_delta))?;
            if radiosity < 0.0 || irradiation < 0.0 {
                return Err(invalid(
                    "pointwise enclosure produced negative radiosity or irradiation",
                ));
            }
            reports.push(PointwiseEnclosurePatchReport {
                name: surface.name().to_owned(),
                emissivity: surface.emissivity().clone(),
                area_m2: surface.area_m2(),
                mean_temperature_k: finite(reference + mean_delta / surface.area_m2())?,
                minimum_temperature_k: minimum[patch],
                maximum_temperature_k: maximum[patch],
                mean_blackbody_emission_w_m2: finite(background + emission_delta[patch])?,
                radiosity_w_m2: radiosity,
                irradiation_w_m2: irradiation,
                outward_heat_w: heat,
            });
        }
        let ambient = if let Some(temperature_k) = self.ambient_temperature_k {
            let sum = finite(factors[p].iter().sum())?;
            let irradiation_delta = dot(&factors[p], &radiosity_delta)?;
            let effective_area_m2 = self.system.view_factors().areas_m2()[p];
            let outward_heat_w = finite(effective_area_m2 * balance::outward_flux(
                1.0, rhs[p], irradiation_delta, background, sum,
            ))?;
            residual_max = residual_max.max(finite(radiosity_delta[p] - rhs[p])?.abs());
            max_area = max_area.max(effective_area_m2);
            let radiosity_w_m2 = finite(background + radiosity_delta[p])?;
            let irradiation_w_m2 = finite(background.mul_add(sum, irradiation_delta))?;
            if radiosity_w_m2 < 0.0 || irradiation_w_m2 < 0.0 {
                return Err(invalid("pointwise black bath produced negative radiosity or irradiation"));
            }
            Some(PointwiseBlackBathReport { temperature_k, effective_area_m2,
                radiosity_w_m2, irradiation_w_m2, outward_heat_w })
        } else { None };
        let ambient_heat = ambient.as_ref().map_or(0.0, |bath| bath.outward_heat_w);
        let system_energy_closure_w = finite(total + ambient_heat)?;
        let total_absolute = finite(absolute + ambient_heat.abs())?;
        let tolerance = finite(
            self.tolerance
                .relative_heat
                .mul_add(total_absolute, self.tolerance.absolute_heat_w),
        )?;
        if system_energy_closure_w.abs() > tolerance || finite(residual_max * max_area)? > tolerance {
            return Err(invalid(
                "pointwise radiosity equations or complete-system heat sum exceed the watt budget",
            ));
        }
        let emission_diagonal = diagonal.clone();
        if TANGENT {
            for (ordinal, &(vertex, ref members)) in self.geometry.memberships.iter().enumerate() {
                if ordinal % ASSEMBLY_TILE == 0 {
                    poll(cx, ordinal)?;
                }
                for &(row, a) in members {
                    let geometry = &self.geometry.patches[row];
                    for &(column, b) in members {
                        diagonal[vertex] = finite(
                            diagonal[vertex]
                                - geometry.emissivity
                                    * geometry.weights[a]
                                    * self.geometry.irradiation_response[row][column]
                                    * slopes[column][b],
                        )?;
                    }
                }
            }
        }
        poll(cx, self.geometry.faces.len())?;
        Ok(PointwiseEnclosureAssembly {
            outward_load_w: load,
            temperature_diagonal_w_k: diagonal,
            report: PointwiseEnclosureReport {
                patches: reports,
                outward_heat_w: total,
                absolute_patch_heat_w: absolute,
                linear_residual_max_w_m2: residual_max,
                view_factor_identity: self.system.view_factors().identity(),
                ambient,
                system_energy_closure_w,
            },
            geometry: Arc::clone(&self.geometry),
            local,
            emission_slopes: slopes,
            emission_diagonal,
            irradiation_corrections,
            emission_minus_irradiation,
        })
    }
}

fn admit_configuration(count: usize, tolerance: PointwiseEnclosureTolerance) -> Result<(), ConductionError> {
    if !(1..=64).contains(&count) || !tolerance.absolute_heat_w.is_finite()
        || tolerance.absolute_heat_w < 0.0 || !tolerance.relative_heat.is_finite()
        || tolerance.relative_heat < 0.0
    {
        Err(invalid("enclosure requires 1..=64 patches and finite nonnegative watt tolerances"))
    } else { Ok(()) }
}

fn dot(a: &[f64], b: &[f64]) -> Result<f64, ConductionError> {
    let mut sum = 0.0;
    for (&x, &y) in a.iter().zip(b) {
        sum = finite(x.mul_add(y, sum))?;
    }
    Ok(sum)
}
fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "pointwise-enclosure",
        at,
    })
}
fn invalid(reason: &str) -> ConductionError {
    radiation_error(
        "pointwise-enclosure",
        reason,
        "supply supported traces, view factors and explicit watt budgets",
    )
}
fn finite(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("nonfinite pointwise enclosure arithmetic"))
    }
}

