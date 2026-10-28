//! Local gray-body radiation on the actual P1 surface temperature field.
//!
//! On each triangle, integrate `N_i epsilon sigma (T_h^4 - T_a^4)` and its
//! exact temperature tangent `4 N_i N_j epsilon sigma T_h^3`. Both integrands
//! have degree five, so the positive seven-point degree-five triangle rule is
//! exact up to floating-point roundoff. No patch-average temperature replaces
//! a hot spot. The frozen, card-backed emissivity and black reservoir are the
//! same declarations used by the explicitly reduced area-mean radiation lane.
//!
//! Geometry and the trace partition are prepared once. The tangent is stored
//! as symmetric 3x3 face blocks and applied without a dense matrix, finite
//! differences, or a global sparse assembly. This is ambient exchange with
//! view factor one, not occlusion, participating media or enclosure radiosity.

use fs_exec::Cx;

use super::{
    AmbientRadiationPatch, STEFAN_BOLTZMANN_W_M2_K4, SURFACE_EMISSIVITY_PROPERTY, radiation_error,
    require_temperature,
};
use crate::{ConductionError, ConductionMesh, ThermalBc, ThermalBoundary, assemble::ASSEMBLY_TILE};

/// Nonlinear heat and material evidence on one spatially resolved trace.
#[derive(Debug, Clone, PartialEq)]
pub struct PointwiseRadiationPatchReport {
    /// Exact emissivity declaration, source receipt and external reservoir.
    pub patch: AmbientRadiationPatch,
    /// Sum of the selected triangle areas, m2.
    pub area_m2: f64,
    /// Exact P1 area-mean temperature, retained only as an observation, K.
    pub mean_surface_temperature_k: f64,
    /// Lowest nodal surface temperature, K.
    pub minimum_surface_temperature_k: f64,
    /// Highest nodal surface temperature, K.
    pub maximum_surface_temperature_k: f64,
    /// Integral of epsilon sigma (T_h^4 - T_a^4), W outward.
    /// Negative values are absorbed heat from a hotter reservoir.
    pub outward_heat_w: f64,
}

/// Direct integrated radiation balance, separate from convective heat.
#[derive(Debug, Clone, PartialEq)]
pub struct PointwiseRadiationReport {
    /// Trace records in the supplied declaration order.
    pub patches: Vec<PointwiseRadiationPatchReport>,
    /// Sum of the signed radiative trace powers, W outward.
    pub outward_heat_w: f64,
}

#[derive(Debug, Clone)]
struct FaceTangent {
    vertices: [usize; 3],
    matrix: [[f64; 3]; 3],
}

/// Exact discrete radiation load and its local analytic temperature tangent.
#[derive(Debug, Clone)]
pub struct PointwiseRadiationAssembly {
    /// Full nodal outward load `integral N_i q_rad dA`, W.
    pub outward_load_w: Vec<f64>,
    /// Diagonal of the temperature tangent, W/K, including shared faces.
    pub temperature_diagonal_w_k: Vec<f64>,
    /// Independently integrated trace powers and retained material declarations.
    pub report: PointwiseRadiationReport,
    tangent: Vec<FaceTangent>,
}

impl PointwiseRadiationAssembly {
    /// Apply the symmetric temperature tangent to a full nodal direction.
    /// The output is replaced, not accumulated; no state or coefficients change.
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

    // The enthalpy owner supplies its exact dT/dh and suppresses prescribed
    // rows. Applying local blocks avoids rebuilding a global matrix per step.
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
        if column_scale.is_some_and(|values| values.len() != direction.len()) {
            return Err(invalid(
                "radiation tangent scaling has the wrong field length",
            ));
        }
        for (slot, face) in self.tangent.iter().enumerate() {
            if slot % ASSEMBLY_TILE == 0 {
                poll(cx, slot)?;
            }
            let mut local = [0.0; 3];
            for (i, &vertex) in face.vertices.iter().enumerate() {
                local[i] = finite(direction[vertex] * column_scale.map_or(1.0, |s| s[vertex]))?;
            }
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
        poll(cx, self.tangent.len())
    }

    fn check_fields(&self, direction: usize, output: usize) -> Result<(), ConductionError> {
        let n = self.outward_load_w.len();
        if direction != n || output != n {
            return Err(invalid(
                "radiation tangent requires full nodal input and output fields",
            ));
        }
        Ok(())
    }
}

/// Prepared local radiation on immutable geometry and a named boundary partition.
///
/// Borrowing the originals makes the binding exact without hashing the entire
/// mesh in each nonlinear residual. Prepare again if geometry or partition
/// changes. Constant emissivity is the retained snapshot; temperature-dependent
/// emissivity, enclosure reflection and moving boundaries need their own model.
#[derive(Debug)]
pub struct PointwiseAmbientRadiation<'a> {
    mesh: &'a ConductionMesh,
    boundary: &'a ThermalBoundary,
    patches: &'a [AmbientRadiationPatch],
    // Mesh-face order fixes reductions regardless of patch declaration order.
    faces: Vec<(usize, usize)>,
    areas: Vec<f64>,
}

impl<'a> PointwiseAmbientRadiation<'a> {
    /// Bind nonempty distinct patches to existing natural boundary regions.
    /// Radiation adds to the original Robin/Neumann law; it never replaces
    /// convection or relabels radiative watts as heat carried by a coolant.
    /// Fixed nodes shared with adjacent Dirichlet faces retain their reactions.
    pub fn new(
        cx: &Cx<'_>,
        mesh: &'a ConductionMesh,
        boundary: &'a ThermalBoundary,
        patches: &'a [AmbientRadiationPatch],
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        if patches.is_empty() || patches.len() > boundary.region_names().len() {
            return Err(invalid(
                "nonempty distinct pointwise radiation regions required",
            ));
        }
        mesh.boundary()
            .len()
            .checked_mul(9 * std::mem::size_of::<f64>())
            .ok_or_else(|| invalid("radiative face storage size overflow"))?;
        let mut owners = vec![None; boundary.region_names().len()];
        for (index, patch) in patches.iter().enumerate() {
            poll(cx, index)?;
            let region = boundary
                .region_names()
                .iter()
                .position(|name| name == patch.region())
                .ok_or_else(|| invalid("unknown pointwise radiation region"))?;
            if owners[region].replace(index).is_some() {
                return Err(invalid("duplicate pointwise radiation region"));
            }
            if matches!(boundary.conditions()[region], ThermalBc::Dirichlet { .. }) {
                return Err(invalid(
                    "pointwise radiation requires a natural boundary trace",
                ));
            }
        }
        let mut faces = Vec::new();
        let mut areas = vec![0.0; patches.len()];
        for (slot, face) in mesh.boundary().iter().enumerate() {
            if slot % ASSEMBLY_TILE == 0 {
                poll(cx, slot)?;
            }
            if let Some(patch) = boundary.region_for(slot).and_then(|region| owners[region]) {
                faces.push((slot, patch));
                areas[patch] = finite(areas[patch] + face.area)?;
            }
        }
        if areas.iter().any(|&area| area <= 0.0) {
            return Err(invalid("pointwise radiation region has no exterior area"));
        }
        poll(cx, faces.len())?;
        Ok(Self {
            mesh,
            boundary,
            patches,
            faces,
            areas,
        })
    }

    pub(crate) fn matches(&self, mesh: &ConductionMesh, boundary: &ThermalBoundary) -> bool {
        std::ptr::eq(self.mesh, mesh) && std::ptr::eq(self.boundary, boundary)
    }

    /// Integrate the actual P1 temperatures with a degree-five positive rule.
    /// Each surface node must lie within the retained emissivity claim's support;
    /// checking only a patch average would conceal unsupported hot/cold spots.
    pub fn assemble(
        &self,
        cx: &Cx<'_>,
        temperature: &[f64],
    ) -> Result<PointwiseRadiationAssembly, ConductionError> {
        self.integrate::<true>(cx, temperature)
    }

    // Line-search trials need the nonlinear load, not another face tangent.
    // A const specialization removes matrix arithmetic and storage from that
    // hot path while retaining exactly the same quadrature and domain checks.
    pub(crate) fn load_and_report(
        &self,
        cx: &Cx<'_>,
        temperature: &[f64],
    ) -> Result<(Vec<f64>, PointwiseRadiationReport), ConductionError> {
        let assembly = self.integrate::<false>(cx, temperature)?;
        Ok((assembly.outward_load_w, assembly.report))
    }

    fn integrate<const WITH_TANGENT: bool>(
        &self,
        cx: &Cx<'_>,
        temperature: &[f64],
    ) -> Result<PointwiseRadiationAssembly, ConductionError> {
        poll(cx, 0)?;
        let n = self.mesh.vertex_count();
        if temperature.len() != n {
            return Err(ConductionError::FieldLength {
                field: "pointwise radiation temperature",
                expected: n,
                found: temperature.len(),
            });
        }
        let mut load = vec![0.0; n];
        let mut diagonal = if WITH_TANGENT {
            vec![0.0; n]
        } else {
            Vec::new()
        };
        let mut tangent = if WITH_TANGENT {
            Vec::with_capacity(self.faces.len())
        } else {
            Vec::new()
        };
        let mut reports: Vec<_> = self
            .patches
            .iter()
            .zip(&self.areas)
            .map(|(patch, &area)| PointwiseRadiationPatchReport {
                patch: patch.clone(),
                area_m2: area,
                mean_surface_temperature_k: 0.0,
                minimum_surface_temperature_k: f64::INFINITY,
                maximum_surface_temperature_k: f64::NEG_INFINITY,
                outward_heat_w: 0.0,
            })
            .collect();
        for (index, &(slot, patch_index)) in self.faces.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            let face = &self.mesh.boundary()[slot];
            let vertices = face.vertices.map(|v| v as usize);
            let t = vertices.map(|v| temperature[v]);
            let patch = &self.patches[patch_index];
            let report = &mut reports[patch_index];
            for value in t {
                require_temperature(patch.region(), "radiating surface node", value)?;
                if let Some((low, high)) = patch.emissivity().temperature_validity_k()
                    && !(low..=high).contains(&value)
                {
                    return Err(ConductionError::MaterialQuery {
                        property: SURFACE_EMISSIVITY_PROPERTY.into(),
                        temperature: value,
                        upstream: format!(
                            "surface {:?} node is outside retained emissivity support [{low}, {high}] K",
                            patch.region()
                        ),
                    });
                }
                report.minimum_surface_temperature_k =
                    report.minimum_surface_temperature_k.min(value);
                report.maximum_surface_temperature_k =
                    report.maximum_surface_temperature_k.max(value);
            }
            let mean = finite(t[0] / 3.0 + t[1] / 3.0 + t[2] / 3.0)?;
            report.mean_surface_temperature_k =
                finite(report.mean_surface_temperature_k + face.area * mean)?;
            let ambient = patch.ambient_temperature_k();
            let coefficient = patch.emissivity().value() * STEFAN_BOLTZMANN_W_M2_K4;
            let mut matrix = [[0.0; 3]; 3];
            for (bary, weight) in TRIANGLE_DEGREE_FIVE {
                // Affine evaluation preserves an exactly uniform field.
                let local_t =
                    finite(t[0] + (t[1] - t[0]).mul_add(bary[1], (t[2] - t[0]) * bary[2]))?;
                let area_weight = finite(face.area * weight)?;
                // Factoring the difference avoids subtracting two nearly equal
                // absolute fourth powers at thermal equilibrium.
                let flux = finite(
                    coefficient
                        * (local_t - ambient)
                        * (local_t + ambient)
                        * local_t.mul_add(local_t, ambient * ambient),
                )?;
                let heat = finite(area_weight * flux)?;
                report.outward_heat_w = finite(report.outward_heat_w + heat)?;
                for i in 0..3 {
                    load[vertices[i]] = finite(load[vertices[i]] + bary[i] * heat)?;
                }
                if WITH_TANGENT {
                    let derivative =
                        finite(area_weight * (4.0 * coefficient * local_t * local_t * local_t))?;
                    for i in 0..3 {
                        for j in 0..=i {
                            matrix[i][j] = finite(matrix[i][j] + derivative * bary[i] * bary[j])?;
                        }
                    }
                }
            }
            if WITH_TANGENT {
                for i in 0..3 {
                    diagonal[vertices[i]] = finite(diagonal[vertices[i]] + matrix[i][i])?;
                    for j in 0..i {
                        matrix[j][i] = matrix[i][j];
                    }
                }
                tangent.push(FaceTangent { vertices, matrix });
            }
        }
        let mut outward_heat_w = 0.0;
        for report in &mut reports {
            report.mean_surface_temperature_k =
                finite(report.mean_surface_temperature_k / report.area_m2)?;
            outward_heat_w = finite(outward_heat_w + report.outward_heat_w)?;
        }
        poll(cx, self.faces.len())?;
        Ok(PointwiseRadiationAssembly {
            outward_load_w: load,
            temperature_diagonal_w_k: diagonal,
            report: PointwiseRadiationReport {
                patches: reports,
                outward_heat_w,
            },
            tangent,
        })
    }
}

// Dunavant's degree-five rule, normalized so weights sum to one. The exact
// noncentroid coordinates are (6 +/- sqrt(15))/21 and their complements;
// their weights are (155 +/- sqrt(15))/1200. Positive interior samples keep
// the passive tangent positive semidefinite. Polynomial exactness is checked
// independently against barycentric simplex moments in the regression tests.
pub(super) const TRIANGLE_DEGREE_FIVE: [([f64; 3], f64); 7] = [
    ([1.0 / 3.0; 3], 9.0 / 40.0),
    (
        [
            0.05971587178976982,
            0.47014206410511509,
            0.47014206410511509,
        ],
        0.13239415278850618,
    ),
    (
        [
            0.47014206410511509,
            0.05971587178976982,
            0.47014206410511509,
        ],
        0.13239415278850618,
    ),
    (
        [
            0.47014206410511509,
            0.47014206410511509,
            0.05971587178976982,
        ],
        0.13239415278850618,
    ),
    (
        [0.7974269853530872, 0.10128650732345634, 0.10128650732345634],
        0.12593918054482715,
    ),
    (
        [0.10128650732345634, 0.7974269853530872, 0.10128650732345634],
        0.12593918054482715,
    ),
    (
        [0.10128650732345634, 0.10128650732345634, 0.7974269853530872],
        0.12593918054482715,
    ),
];

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "pointwise-radiation",
        at,
    })
}
fn invalid(reason: &str) -> ConductionError {
    radiation_error(
        "pointwise-ambient",
        reason,
        "supply supported surface temperatures and a matching natural trace",
    )
}
fn finite(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("nonfinite pointwise radiation arithmetic"))
    }
}
