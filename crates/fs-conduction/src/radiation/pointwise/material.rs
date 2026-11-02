//! Source-backed nodal emissivity on the original ambient P1 radiation trace.
//!
//! The declared discretization interpolates nodal temperature AND emissivity.
//! Original caloric states supply each node's T and liquid mass fraction. This
//! is not a mixture rule, quadrature-point inversion of interpolated enthalpy,
//! spectral radiation, or a replacement thermal/nonlinear solver.

use std::collections::BTreeSet;

use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint};
use fs_material::thermal::source::ThermalMaterialPoint;

use super::*;
use crate::SurfaceEmissivity;

pub(in crate::radiation) mod source;
pub mod adjoint;
pub use source::MaterialRadiationSmoothnessError;
use source::{DerivativeMode, MaterialEmissivitySource};

/// One homogeneous material trace and its independently declared black bath.
/// Different material charts require different named traces. Other conditions
/// such as pressure and finish remain on the original material binding.
#[derive(Debug)]
pub struct MaterialAmbientRadiationPatch<'a> {
    pub region: &'a str,
    pub ambient_temperature_k: f64,
    /// The SAME caloric object used by every participating storage node.
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    pub material: ThermalMaterialPoint<'a>,
    /// Resolved once at preparation; that exact claim stays pinned thereafter.
    pub selection: MaterialPropertySelection,
}

/// Actual material/phase sample defining one nodal emissivity coefficient.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialRadiationNodeState {
    pub vertex: usize,
    pub phase: EquilibriumPhaseState,
    pub emissivity: SurfaceEmissivity,
    pub material: ResolvedMaterialStatePoint,
}

#[derive(Debug)]
struct PreparedPatch<'a> {
    snapshot: AmbientRadiationPatch,
    source: MaterialEmissivitySource<'a>,
    vertices: Vec<usize>,
}

/// Immutable prepared trace and original source laws. Every physical query
/// refreshes emissivity at its actual nodal h, including latent-phase changes.
/// Supported sources are static scalars and original linear curves in T or
/// liquid mass fraction. No source switching, extrapolation or phase averaging.
#[derive(Debug)]
pub struct MaterialPointwiseAmbientRadiation<'a> {
    mesh: &'a ConductionMesh,
    boundary: &'a ThermalBoundary,
    patches: Vec<PreparedPatch<'a>>,
    faces: Vec<(usize, usize)>,
    face_samples: Vec<[usize; 3]>,
    areas: Vec<f64>,
    curves: Vec<Option<&'a EquilibriumEnthalpyPhaseCurve>>,
}

impl<'a> MaterialPointwiseAmbientRadiation<'a> {
    pub fn new(
        cx: &Cx<'_>, mesh: &'a ConductionMesh, boundary: &'a ThermalBoundary,
        patches: &[MaterialAmbientRadiationPatch<'a>],
    ) -> Result<Self, ConductionError> {
        let names: Vec<_> = patches.iter().map(|patch| patch.region).collect();
        let (faces, areas) = prepare_geometry(cx, mesh, boundary, &names)?;
        let mut vertices = vec![BTreeSet::new(); patches.len()];
        for (index, &(face, patch)) in faces.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 { poll(cx, index)?; }
            for vertex in mesh.boundary()[face].vertices { vertices[patch].insert(vertex as usize); }
        }
        let mut prepared = Vec::with_capacity(patches.len());
        let mut curves: Vec<Option<&'a EquilibriumEnthalpyPhaseCurve>> = vec![None; mesh.vertex_count()];
        for (index, patch) in patches.iter().enumerate() {
            poll(cx, index)?;
            let source = MaterialEmissivitySource::prepare(
                patch.region, patch.curve, patch.material, patch.selection.clone(),
            )?;
            let local: Vec<_> = vertices[index].iter().copied().collect();
            for (index, &vertex) in local.iter().enumerate() {
                if index % ASSEMBLY_TILE == 0 { poll(cx, index)?; }
                if curves[vertex].is_some_and(|curve| !std::ptr::eq(curve, patch.curve)) {
                    return Err(invalid("shared surface vertices require the same original caloric object"));
                }
                curves[vertex] = Some(patch.curve);
            }
            prepared.push(PreparedPatch {
                snapshot: AmbientRadiationPatch::new(patch.region, source.snapshot().clone(), patch.ambient_temperature_k)?,
                source, vertices: local,
            });
        }
        let mut face_samples = Vec::with_capacity(faces.len());
        for (index, &(face, patch)) in faces.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 { poll(cx, index)?; }
            face_samples.push(mesh.boundary()[face].vertices.map(|vertex|
                prepared[patch].vertices.binary_search(&(vertex as usize)).expect("prepared trace vertex")));
        }
        poll(cx, faces.len())?;
        Ok(Self { mesh, boundary, patches: prepared, faces, face_samples, areas, curves })
    }

    pub(crate) fn matches(&self, mesh: &ConductionMesh, boundary: &ThermalBoundary) -> bool {
        std::ptr::eq(self.mesh, mesh) && std::ptr::eq(self.boundary, boundary)
    }

    pub(crate) fn matches_curves<'c>(
        &self, cx: &Cx<'_>,
        mut curve_at: impl FnMut(usize) -> &'c EquilibriumEnthalpyPhaseCurve,
    ) -> Result<bool, ConductionError> {
        for (vertex, curve) in self.curves.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 { poll(cx, vertex)?; }
            if curve.is_some_and(|curve| !std::ptr::eq(curve, curve_at(vertex))) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Assemble the physical nodal load and its complete h-coordinate tangent.
    /// Call `enthalpy_tangent_apply` on the result. The original positive
    /// emission-temperature diagonal remains available for preconditioning.
    pub fn assemble(&self, cx: &Cx<'_>, specific_h: &[f64])
        -> Result<PointwiseRadiationAssembly, ConductionError>
    {
        self.integrate::<true>(cx, specific_h, None)
    }

    pub(crate) fn assemble_at_temperature(&self, cx: &Cx<'_>, specific_h: &[f64], temperature: &[f64])
        -> Result<PointwiseRadiationAssembly, ConductionError>
    {
        self.integrate::<true>(cx, specific_h, Some(temperature))
    }

    pub(crate) fn load_and_report(&self, cx: &Cx<'_>, specific_h: &[f64], temperature: &[f64])
        -> Result<(Vec<f64>, PointwiseRadiationReport), ConductionError>
    {
        let assembly = self.integrate::<false>(cx, specific_h, Some(temperature))?;
        Ok((assembly.outward_load_w, assembly.report))
    }

    fn integrate<const WITH_TANGENT: bool>(
        &self, cx: &Cx<'_>, specific_h: &[f64], supplied_temperature: Option<&[f64]>,
    ) -> Result<PointwiseRadiationAssembly, ConductionError> {
        poll(cx, 0)?;
        let n = self.mesh.vertex_count();
        if specific_h.len() != n || supplied_temperature.is_some_and(|temperature| temperature.len() != n) {
            return Err(invalid("material radiation requires complete nodal enthalpy and temperature fields"));
        }
        let mut reports = Vec::with_capacity(self.patches.len());
        let mut coefficients = Vec::with_capacity(self.patches.len());
        let mut temperature_slopes = if WITH_TANGENT { vec![0.0; n] } else { Vec::new() };
        for (index, patch) in self.patches.iter().enumerate() {
            let mut states = Vec::with_capacity(patch.vertices.len());
            let mut values = Vec::with_capacity(patch.vertices.len());
            let mut minimum = f64::INFINITY;
            let mut maximum = f64::NEG_INFINITY;
            for &vertex in &patch.vertices {
                poll(cx, vertex)?;
                let response = patch.source.resolve(vertex, specific_h[vertex],
                    supplied_temperature.map(|temperature| temperature[vertex]),
                    if WITH_TANGENT { DerivativeMode::Primal } else { DerivativeMode::None })?;
                let temperature = response.state.phase.temperature_k();
                if WITH_TANGENT { temperature_slopes[vertex] = response.dt_dh; }
                minimum = minimum.min(temperature); maximum = maximum.max(temperature);
                values.push((temperature, response.state.emissivity.value(), response.depsilon_dh));
                states.push(response.state);
            }
            coefficients.push(values);
            reports.push(PointwiseRadiationPatchReport {
                patch: patch.snapshot.clone(), area_m2: self.areas[index],
                mean_surface_temperature_k: 0.0, minimum_surface_temperature_k: minimum,
                maximum_surface_temperature_k: maximum, outward_heat_w: 0.0, material_states: states,
            });
        }
        let mut load = vec![0.0; n];
        let mut diagonal = if WITH_TANGENT { vec![0.0; n] } else { Vec::new() };
        let mut tangent = Vec::new(); let mut enthalpy_feedback = Vec::new();
        for (index, &(slot, patch)) in self.faces.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 { poll(cx, index)?; }
            let face = &self.mesh.boundary()[slot];
            let vertices = face.vertices.map(|vertex| vertex as usize);
            let local = self.face_samples[index].map(|sample| coefficients[patch][sample]);
            let t = local.map(|value| value.0); let epsilon = local.map(|value| value.1);
            let epsilon_h = local.map(|value| value.2);
            let ambient = self.patches[patch].snapshot.ambient_temperature_k();
            reports[patch].mean_surface_temperature_k = finite(reports[patch].mean_surface_temperature_k
                + face.area * (t[0] / 3.0 + t[1] / 3.0 + t[2] / 3.0))?;
            let mut emission = [[0.0; 3]; 3]; let mut feedback = [[0.0; 3]; 3];
            for (bary, weight) in TRIANGLE_DEGREE_SIX {
                let local_t = finite(t[0] + (t[1] - t[0]).mul_add(bary[1], (t[2] - t[0]) * bary[2]))?;
                let local_epsilon = finite(epsilon[0] + (epsilon[1] - epsilon[0]).mul_add(
                    bary[1], (epsilon[2] - epsilon[0]) * bary[2]))?;
                let measure = finite(face.area * weight * STEFAN_BOLTZMANN_W_M2_K4)?;
                let difference = finite((local_t - ambient) * (local_t + ambient)
                    * local_t.mul_add(local_t, ambient * ambient))?;
                let heat = finite(measure * local_epsilon * difference)?;
                reports[patch].outward_heat_w = finite(reports[patch].outward_heat_w + heat)?;
                for i in 0..3 {
                    load[vertices[i]] = finite(bary[i].mul_add(heat, load[vertices[i]]))?;
                    if WITH_TANGENT {
                        for j in 0..3 {
                            let common = finite(measure * bary[i] * bary[j])?;
                            emission[i][j] = finite(emission[i][j]
                                + common * 4.0 * local_epsilon * local_t * local_t * local_t)?;
                            feedback[i][j] = finite(feedback[i][j] + common * difference * epsilon_h[j])?;
                        }
                    }
                }
            }
            if WITH_TANGENT {
                for i in 0..3 { diagonal[vertices[i]] = finite(diagonal[vertices[i]] + emission[i][i])?; }
                let block = |matrix| FaceTangent { vertices, patch, matrix,
                    emissivity_load: [0.0; 3], reservoir_load: [0.0; 3] };
                tangent.push(block(emission)); enthalpy_feedback.push(block(feedback));
            }
        }
        let mut outward_heat_w = 0.0;
        for report in &mut reports {
            report.mean_surface_temperature_k = finite(report.mean_surface_temperature_k / report.area_m2)?;
            outward_heat_w = finite(outward_heat_w + report.outward_heat_w)?;
        }
        poll(cx, self.faces.len())?;
        Ok(PointwiseRadiationAssembly { outward_load_w: load, temperature_diagonal_w_k: diagonal,
            report: PointwiseRadiationReport { patches: reports, outward_heat_w }, tangent,
            enthalpy_feedback, material_temperature_slopes: WITH_TANGENT.then_some(temperature_slopes),
            material_parameters: None })
    }
}

// Four-point Gauss-Legendre on [0,1], with the Duffy map
// N=[u,(1-u)v,(1-u)(1-v)]. Degree-six barycentric polynomials become
// degree <=7 in u (including the Jacobian) and <=6 in v. All 16 points
// and weights are positive; weights below integrate a unit-normalized area.
const GAUSS_FOUR: [(f64, f64); 4] = [
    (0.06943184420297371, 0.17392742256872692),
    (0.33000947820757187, 0.32607257743127307),
    (0.6699905217924281, 0.32607257743127307),
    (0.9305681557970262, 0.17392742256872692),
];
pub(in crate::radiation) const TRIANGLE_DEGREE_SIX: [([f64; 3], f64); 16] = {
    let mut result = [([0.0; 3], 0.0); 16];
    let mut i = 0;
    while i < 16 {
        let (u, wu) = GAUSS_FOUR[i / 4]; let (v, wv) = GAUSS_FOUR[i % 4];
        result[i] = ([u, (1.0 - u) * v, (1.0 - u) * (1.0 - v)], 2.0 * wu * wv * (1.0 - u));
        i += 1;
    }
    result
};

#[cfg(test)]
mod tests {
    use super::TRIANGLE_DEGREE_SIX;

    #[test]
    fn positive_duffy_rule_matches_all_exact_simplex_moments_through_degree_six() {
        let factorial = |n: i32| (1..=n).map(f64::from).product::<f64>();
        for degree in 0..=6 {
            for a in 0..=degree {
                for b in 0..=degree - a {
                    let c = degree - a - b;
                    let exact = 2.0 * factorial(a) * factorial(b) * factorial(c) / factorial(degree + 2);
                    let actual: f64 = TRIANGLE_DEGREE_SIX.iter().map(|(n, w)| {
                        assert!(*w > 0.0 && n.iter().all(|v| *v > 0.0));
                        w * n[0].powi(a) * n[1].powi(b) * n[2].powi(c)
                    }).sum();
                    assert!((actual - exact).abs() <= 3e-15 * exact.max(1e-3), "{a},{b},{c}: {actual} != {exact}");
                }
            }
        }
    }
}
