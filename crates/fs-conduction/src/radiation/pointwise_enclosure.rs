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
    GrayDiffuseEnclosure, PreparedRadiosity, STEFAN_BOLTZMANN_W_M2_K4 as SIGMA,
    SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, balance, pointwise::TRIANGLE_DEGREE_FIVE,
    radiation_error,
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

/// Closed-system radiation accounting, separate from external convection/flux.
#[derive(Debug, Clone, PartialEq)]
pub struct PointwiseEnclosureReport {
    pub patches: Vec<PointwiseEnclosurePatchReport>,
    /// Actual signed sum; never projected to zero.
    pub outward_heat_w: f64,
    pub absolute_patch_heat_w: f64,
    pub linear_residual_max_w_m2: f64,
    pub view_factor_identity: ContentHash,
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
}

#[derive(Debug, Clone)]
struct FaceTangent {
    vertices: [usize; 3],
    matrix: [[f64; 3]; 3],
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

    fn check_fields(&self, direction: usize, output: usize) -> Result<(), ConductionError> {
        if direction != self.geometry.vertices || output != self.geometry.vertices {
            return Err(invalid("enclosure tangent requires full nodal fields"));
        }
        Ok(())
    }
}

/// Prepared fixed-geometry enclosure, with sparse trace geometry and one small
/// radiosity factorization. At most 64 patches; trace storage is linear in
/// boundary faces/vertices and the reflected matrix is O(patches^2).
#[derive(Debug)]
pub struct PointwiseEnclosureRadiation<'a> {
    mesh: &'a ConductionMesh,
    boundary: &'a ThermalBoundary,
    enclosure: &'a GrayDiffuseEnclosure,
    factors: PreparedRadiosity<'a>,
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
        let count = enclosure.surfaces().len();
        if !(2..=64).contains(&count)
            || !tolerance.absolute_heat_w.is_finite()
            || tolerance.absolute_heat_w < 0.0
            || !tolerance.relative_heat.is_finite()
            || tolerance.relative_heat < 0.0
        {
            return Err(invalid(
                "enclosure requires 2..=64 patches and finite nonnegative watt tolerances",
            ));
        }
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
        for (patch, surface) in enclosure.surfaces().iter().enumerate() {
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
        let factors = enclosure.prepare_radiosity(cx)?;
        let mut irradiation_response = vec![vec![0.0; count]; count];
        for column in 0..count {
            poll(cx, column)?;
            let mut rhs = vec![0.0; count];
            rhs[column] = patches[column].emissivity;
            let response = factors.solve_rhs(cx, &rhs)?;
            for (row, values) in irradiation_response.iter_mut().enumerate() {
                values[column] = dot(&enclosure.view_factors().factors()[row], &response)?;
            }
        }
        Ok(Self {
            mesh,
            boundary,
            enclosure,
            factors,
            tolerance,
            geometry: Arc::new(Geometry {
                vertices: mesh.vertex_count(),
                patches,
                faces,
                memberships: memberships.into_iter().collect(),
                irradiation_response,
            }),
        })
    }

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
            let surface = &self.enclosure.surfaces()[patch];
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
        let reference = minimum.iter().copied().fold(f64::INFINITY, f64::min);
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
                    matrix,
                });
            }
        }
        let factors = self.enclosure.view_factors().factors();
        let mut row_sums = Vec::with_capacity(p);
        let mut rhs = Vec::with_capacity(p);
        for patch in 0..p {
            poll(cx, patch)?;
            let area = self.enclosure.surfaces()[patch].area_m2();
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
        let radiosity_delta = self.factors.solve_rhs(cx, &rhs)?;
        let mut reports = Vec::with_capacity(p);
        let mut total = 0.0;
        let mut absolute = 0.0;
        let mut residual_max = 0.0_f64;
        let mut max_area = 0.0_f64;
        for (patch, geometry) in self.geometry.patches.iter().enumerate() {
            poll(cx, patch)?;
            let surface = &self.enclosure.surfaces()[patch];
            let irradiation_delta = dot(&factors[patch], &radiosity_delta)?;
            let residual = finite(
                radiosity_delta[patch]
                    - (1.0 - geometry.emissivity) * irradiation_delta
                    - rhs[patch],
            )?;
            residual_max = residual_max.max(residual.abs());
            max_area = max_area.max(surface.area_m2());
            let correction = finite(background.mul_add(row_sums[patch] - 1.0, irradiation_delta))?;
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
        let tolerance = finite(
            self.tolerance
                .relative_heat
                .mul_add(absolute, self.tolerance.absolute_heat_w),
        )?;
        if total.abs() > tolerance || finite(residual_max * max_area)? > tolerance {
            return Err(invalid(
                "pointwise radiosity equations or closed heat sum exceed the watt budget",
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
                view_factor_identity: self.enclosure.view_factors().identity(),
            },
            geometry: Arc::clone(&self.geometry),
            local,
            emission_slopes: slopes,
            emission_diagonal,
        })
    }
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

