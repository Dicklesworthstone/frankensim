//! Original scalar/curve-ordinate controls of an admitted material trace.
//! Geometry, knot abscissae, caloric chart and fixed source conditions stay fixed.

use super::source::{Coordinate, Response};
use super::*;
use fs_blake3::ContentHash;
use fs_matdb::ClaimId;

/// Derivatives with respect to original emissivity values, not nodal snapshots.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialEmissivityGradient {
    Scalar {
        derivative: f64,
    },
    Curve {
        axis: String,
        knot_abscissae: Vec<f64>,
        knot_value_derivatives: Vec<f64>,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct MaterialRadiationPatchGradient {
    pub region: String,
    pub material_identity: ContentHash,
    pub claim: ClaimId,
    pub emissivity: MaterialEmissivityGradient,
    /// Derivative with respect to this patch's independent black-bath K.
    pub reservoir_temperature: f64,
}

/// Original trace order. If a source value is shared across patches, sum its
/// entries with the same material/claim identity; geometry is not a control.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialRadiationParameterGradient {
    pub patches: Vec<MaterialRadiationPatchGradient>,
}

impl MaterialEmissivitySource<'_> {
    pub(in crate::radiation) fn control_count(&self) -> usize {
        match self.response {
            Response::Constant => 1,
            Response::Curve { knots, .. } => knots.len(),
        }
    }

    pub(in crate::radiation) fn control_template(
        &self,
        cx: &Cx<'_>,
    ) -> Result<MaterialEmissivityGradient, ConductionError> {
        Ok(match self.response {
            Response::Constant => MaterialEmissivityGradient::Scalar { derivative: 0.0 },
            Response::Curve { coordinate, knots } => {
                let mut knot_abscissae = Vec::with_capacity(knots.len());
                for (index, knot) in knots.iter().enumerate() {
                    if index % ASSEMBLY_TILE == 0 { poll(cx, index)?; }
                    knot_abscissae.push(knot.0);
                }
                MaterialEmissivityGradient::Curve {
                    axis: match coordinate {
                        Coordinate::Temperature => self.material.temperature_axis,
                        Coordinate::LiquidFraction => self.material.liquid_fraction_axis
                            .expect("prepared fraction axis"),
                    }.into(),
                    knot_abscissae,
                    knot_value_derivatives: vec![0.0; knots.len()],
                }
            }
        })
    }

    /// Original ordinate basis at an already resolved, strictly admitted state.
    pub(in crate::radiation) fn control_basis(
        &self,
        phase: EquilibriumPhaseState,
    ) -> Result<[(usize, f64); 2], ConductionError> {
        Ok(match self.response {
            Response::Constant => [(0, 1.0), (0, 0.0)],
            Response::Curve { coordinate, knots } => {
                let value = match coordinate {
                    Coordinate::Temperature => phase.temperature_k(),
                    Coordinate::LiquidFraction => phase.liquid_mass_fraction(),
                };
                let upper = knots.partition_point(|knot| knot.0 <= value).min(knots.len() - 1);
                let fraction = finite((value - knots[upper - 1].0) / (knots[upper].0 - knots[upper - 1].0))?;
                [(upper - 1, 1.0 - fraction), (upper, fraction)]
            }
        })
    }
}

#[derive(Debug, Clone)]
struct FaceParameters {
    vertices: [usize; 3],
    patch: usize,
    emissivity: [[f64; 3]; 3],
    reservoir: [f64; 3],
    basis: [[(usize, f64); 2]; 3],
}

#[derive(Debug, Clone)]
pub(crate) struct MaterialRadiationParameterData {
    template: MaterialRadiationParameterGradient,
    faces: Vec<FaceParameters>,
}

impl MaterialPointwiseAmbientRadiation<'_> {
    /// Conservative numerical/index-word cap for the retained full h action,
    /// source-control data and contraction scratch. Source strings and immutable
    /// cards are excluded, as in the existing pointwise adjoint allowance.
    pub fn adjoint_workspace_entries(&self) -> Result<usize, ConductionError> {
        let fail = || invalid("material radiation adjoint workspace overflow");
        if self.patches.is_empty() || self.patches.len() > 64 {
            return Err(invalid(
                "material radiation adjoint requires one to 64 patches",
            ));
        }
        let parameters = self.patches.iter().try_fold(0usize, |count, patch| {
            count
                .checked_add(patch.source.control_count())
                .ok_or_else(fail)
        })?;
        self.mesh
            .vertex_count()
            .checked_mul(8)
            .and_then(|n| {
                self.faces
                    .len()
                    .checked_mul(128)
                    .and_then(|f| n.checked_add(f))
            })
            .and_then(|n| parameters.checked_mul(4).and_then(|p| n.checked_add(p)))
            .and_then(|n| {
                self.patches
                    .len()
                    .checked_mul(64)
                    .and_then(|p| n.checked_add(p))
            })
            .filter(|&n| n <= isize::MAX as usize / std::mem::size_of::<usize>())
            .ok_or_else(fail)
    }

    pub(crate) fn bind_adjoint(
        &self,
        cx: &Cx<'_>,
        h: &[f64],
        assembly: &mut PointwiseRadiationAssembly,
        max_entries: usize,
    ) -> Result<(), ConductionError> {
        poll(cx, 0)?;
        if self.adjoint_workspace_entries()? > max_entries {
            return Err(invalid(
                "material radiation adjoint exceeds its declared entry allowance",
            ));
        }
        if h.len() != self.mesh.vertex_count()
            || assembly.report.patches.len() != self.patches.len()
            || assembly.material_temperature_slopes.is_none()
        {
            return Err(invalid(
                "material radiation adjoint needs the original complete h assembly",
            ));
        }
        let mut template = MaterialRadiationParameterGradient {
            patches: Vec::with_capacity(self.patches.len()),
        };
        let mut bases = Vec::with_capacity(self.patches.len());
        for (index, patch) in self.patches.iter().enumerate() {
            poll(cx, index)?;
            let source = &patch.source;
            let emissivity = source.control_template(cx)?;
            template.patches.push(MaterialRadiationPatchGradient {
                region: source.region.clone(),
                material_identity: source.material.material.content_hash(),
                claim: source.claim,
                emissivity,
                reservoir_temperature: 0.0,
            });
            let mut basis = Vec::with_capacity(patch.vertices.len());
            let states = &assembly.report.patches[index].material_states;
            if states.len() != patch.vertices.len() {
                return Err(invalid(
                    "material adjoint nodal source receipts have changed",
                ));
            }
            for (sample, &vertex) in patch.vertices.iter().enumerate() {
                poll(cx, vertex)?;
                let response = source.resolve(
                    vertex,
                    h[vertex],
                    Some(states[sample].phase.temperature_k()),
                    DerivativeMode::Strict,
                )?;
                if response.state != states[sample] {
                    return Err(invalid(
                        "material adjoint source state differs from the accepted endpoint",
                    ));
                }
                basis.push(source.control_basis(response.state.phase)?);
            }
            bases.push(basis);
        }
        let mut faces = Vec::with_capacity(self.faces.len());
        for (index, &(slot, patch)) in self.faces.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            let face = &self.mesh.boundary()[slot];
            let vertices = face.vertices.map(|v| v as usize);
            let samples = self.face_samples[index];
            let states = &assembly.report.patches[patch].material_states;
            let t = samples.map(|sample| states[sample].phase.temperature_k());
            let epsilon = samples.map(|sample| states[sample].emissivity.value());
            let bath = self.patches[patch].snapshot.ambient_temperature_k();
            let mut emissivity = [[0.0; 3]; 3];
            let mut reservoir = [0.0; 3];
            for (bary, weight) in TRIANGLE_DEGREE_SIX {
                let local_t =
                    finite(t[0] + (t[1] - t[0]).mul_add(bary[1], (t[2] - t[0]) * bary[2]))?;
                let local_epsilon = finite(
                    epsilon[0]
                        + (epsilon[1] - epsilon[0])
                            .mul_add(bary[1], (epsilon[2] - epsilon[0]) * bary[2]),
                )?;
                let measure = finite(face.area * weight * STEFAN_BOLTZMANN_W_M2_K4)?;
                let difference = finite(
                    (local_t - bath) * (local_t + bath) * local_t.mul_add(local_t, bath * bath),
                )?;
                let bath_derivative = finite(-4.0 * local_epsilon * bath * bath * bath)?;
                for i in 0..3 {
                    reservoir[i] = finite(reservoir[i] + measure * bary[i] * bath_derivative)?;
                    for j in 0..3 {
                        emissivity[i][j] =
                            finite(emissivity[i][j] + measure * bary[i] * bary[j] * difference)?;
                    }
                }
            }
            faces.push(FaceParameters {
                vertices,
                patch,
                emissivity,
                reservoir,
                basis: samples.map(|sample| bases[patch][sample]),
            });
        }
        poll(cx, self.faces.len())?;
        assembly.material_parameters = Some(MaterialRadiationParameterData { template, faces });
        Ok(())
    }
}

impl PointwiseRadiationAssembly {
    /// Contract original scalar/curve values and bath temperature at a strictly
    /// admitted endpoint. These are outward-load weights; a solved thermal
    /// objective supplies `-dt*lambda` on free rows, not an incoming heat weight.
    pub fn material_parameter_pullback(
        &self,
        cx: &Cx<'_>,
        outward_load_weights: &[f64],
    ) -> Result<MaterialRadiationParameterGradient, ConductionError> {
        self.check_fields(outward_load_weights.len(), outward_load_weights.len())?;
        let data = self.material_parameters.as_ref()
            .ok_or_else(|| invalid("material emissivity controls require an explicitly admitted original-source adjoint"))?;
        for (vertex, &weight) in outward_load_weights.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            finite(weight)?;
        }
        let mut result = data.template.clone();
        for (index, face) in data.faces.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            let output = &mut result.patches[face.patch];
            for i in 0..3 {
                output.reservoir_temperature = finite(face.reservoir[i].mul_add(
                    outward_load_weights[face.vertices[i]],
                    output.reservoir_temperature,
                ))?;
            }
            for j in 0..3 {
                let mut nodal = 0.0;
                for i in 0..3 {
                    nodal = finite(
                        face.emissivity[i][j]
                            .mul_add(outward_load_weights[face.vertices[i]], nodal),
                    )?;
                }
                for (parameter, weight) in face.basis[j] {
                    let target = match &mut output.emissivity {
                        MaterialEmissivityGradient::Scalar { derivative } => derivative,
                        MaterialEmissivityGradient::Curve {
                            knot_value_derivatives,
                            ..
                        } => &mut knot_value_derivatives[parameter],
                    };
                    *target = finite(nodal.mul_add(weight, *target))?;
                }
            }
        }
        poll(cx, data.faces.len())?;
        Ok(result)
    }
}
