//! Original emissivity controls through the current diffuse-gray reflection.
//! Geometry, view factors, chart and knot abscissae are held fixed.

use super::*;
use crate::radiation::pointwise::material::adjoint::MaterialEmissivityGradient;
use fs_blake3::ContentHash;
use fs_matdb::ClaimId;

#[derive(Debug, Clone, PartialEq)]
pub struct MaterialEnclosurePatchGradient {
    pub region: String,
    pub material_identity: ContentHash,
    pub claim: ClaimId,
    pub emissivity: MaterialEmissivityGradient,
}

/// Original patch/source order, with one shared black-bath control for an open
/// enclosure. Sum entries with the same material/claim when a source is shared.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialEnclosureParameterGradient {
    pub patches: Vec<MaterialEnclosurePatchGradient>,
    /// Derivative with respect to the open black bath's temperature, in K.
    pub reservoir_temperature: Option<f64>,
}

#[derive(Debug, Clone)]
struct FaceParameters {
    vertices: [usize; 3],
    patch: usize,
    /// Integral of Na Nb (sigma T_h^4 - G_patch), in W.
    emissivity: [[f64; 3]; 3],
    basis: [[(usize, f64); 2]; 3],
}

#[derive(Debug, Clone)]
pub(super) struct MaterialEnclosureParameterData {
    template: MaterialEnclosureParameterGradient,
    faces: Vec<FaceParameters>,
}

impl MaterialPointwiseEnclosureRadiation<'_> {
    /// Conservative numerical/index-word cap for the retained h action,
    /// original source controls, current patch factors and contraction scratch.
    /// Immutable cards and source strings are excluded.
    pub fn adjoint_workspace_entries(&self) -> Result<usize, ConductionError> {
        let fail = || invalid("material enclosure adjoint workspace overflow");
        if self.sources.is_empty() || self.sources.len() > 64 {
            return Err(invalid("material enclosure adjoint requires one to 64 patches"));
        }
        let controls = self.sources.iter().try_fold(0usize, |n, source| {
            n.checked_add(source.control_count()).ok_or_else(fail)
        })?;
        let equations = self.view_factors.len();
        self.base.geometry.vertices.checked_mul(8)
            .and_then(|n| self.base.geometry.faces.len().checked_mul(128).and_then(|f| n.checked_add(f)))
            .and_then(|n| equations.checked_mul(equations).and_then(|p| p.checked_mul(8)).and_then(|p| n.checked_add(p)))
            .and_then(|n| equations.checked_mul(64).and_then(|p| n.checked_add(p)))
            .and_then(|n| controls.checked_mul(4).and_then(|p| n.checked_add(p)))
            .filter(|&n| n <= isize::MAX as usize / std::mem::size_of::<usize>())
            .ok_or_else(fail)
    }

    pub(crate) fn bind_adjoint(
        &self,
        cx: &Cx<'_>,
        h: &[f64],
        assembly: &mut MaterialPointwiseEnclosureAssembly,
        max_entries: usize,
    ) -> Result<(), ConductionError> {
        poll(cx, 0)?;
        if self.adjoint_workspace_entries()? > max_entries {
            return Err(invalid("material enclosure adjoint exceeds its declared entry allowance"));
        }
        if h.len() != self.base.geometry.vertices
            || !Arc::ptr_eq(&assembly.geometry, &self.base.geometry)
            || !Arc::ptr_eq(&assembly.view_factors, &self.view_factors)
            || assembly.report.patches.len() != self.sources.len()
            || assembly.patch_slopes.len() != self.sources.len()
            || assembly.local.len() != self.base.geometry.faces.len()
        {
            return Err(invalid("material enclosure adjoint needs its original complete h assembly"));
        }
        let mut template = MaterialEnclosureParameterGradient {
            patches: Vec::with_capacity(self.sources.len()),
            reservoir_temperature: self.base.ambient_temperature_k.map(|_| 0.0),
        };
        let mut bases = Vec::with_capacity(self.sources.len());
        for (patch, source) in self.sources.iter().enumerate() {
            poll(cx, patch)?;
            let geometry = &self.base.geometry.patches[patch];
            let states = &assembly.report.patches[patch].material_states;
            if states.len() != geometry.vertices.len() {
                return Err(invalid("material enclosure nodal source receipts have changed"));
            }
            template.patches.push(MaterialEnclosurePatchGradient {
                region: self.base.system.surfaces()[patch].name().into(),
                material_identity: source.snapshot().card_identity(),
                claim: source.claim(),
                emissivity: source.control_template(cx)?,
            });
            let mut basis = Vec::with_capacity(states.len());
            for (slot, &vertex) in geometry.vertices.iter().enumerate() {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                let response = source.resolve(vertex, h[vertex],
                    Some(states[slot].phase.temperature_k()), DerivativeMode::Strict)?;
                if response.state != states[slot] {
                    return Err(invalid("material enclosure source state differs from the accepted endpoint"));
                }
                basis.push(source.control_basis(response.state.phase)?);
            }
            bases.push(basis);
        }
        let mut faces = Vec::with_capacity(self.base.geometry.faces.len());
        for (slot, face) in self.base.geometry.faces.iter().enumerate() {
            if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
            let states = &assembly.report.patches[face.patch].material_states;
            let temperature = face.local.map(|i| states[i].phase.temperature_k());
            let mut emissivity = [[0.0; 3]; 3];
            for (bary, weight) in TRIANGLE_DEGREE_SIX {
                let t = finite(temperature[0] + (temperature[1] - temperature[0])
                    .mul_add(bary[1], (temperature[2] - temperature[0]) * bary[2]))?;
                let net = finite(balance::emission_difference(SIGMA, t, assembly.reference_temperature_k)
                    - assembly.irradiation_corrections_w_m2[face.patch])?;
                let measure = finite(face.area * weight * net)?;
                for i in 0..3 { for j in 0..3 {
                    emissivity[i][j] = finite((bary[i] * bary[j]).mul_add(measure, emissivity[i][j]))?;
                } }
            }
            faces.push(FaceParameters { vertices: face.vertices, patch: face.patch, emissivity,
                basis: face.local.map(|i| bases[face.patch][i]) });
        }
        poll(cx, faces.len())?;
        assembly.material_parameters = Some(MaterialEnclosureParameterData { template, faces });
        Ok(())
    }
}

impl MaterialPointwiseEnclosureAssembly {
    /// Contract original scalar/curve ordinates and the open black bath using
    /// outward-load weights. The thermal owner supplies -dt*lambda on free rows.
    ///
    /// With K=I-diag(1-<epsilon>)F, solve K^T z=F^T a once using the retained
    /// current factors, where a_i=sum_a w_a integral(Na epsilon). Each original
    /// nodal emissivity contributes sum_a (w_a-z_i/A_i) integral(Na Nb (B-G_i)).
    /// The G_i term differentiates reflection as well as emitted power.
    pub fn material_parameter_pullback(
        &self,
        cx: &Cx<'_>,
        outward_load_weights: &[f64],
    ) -> Result<MaterialEnclosureParameterGradient, ConductionError> {
        poll(cx, 0)?;
        self.check_fields(outward_load_weights, outward_load_weights)?;
        let data = self.material_parameters.as_ref().ok_or_else(||
            invalid("material enclosure controls require an explicitly admitted original-source adjoint"))?;
        for (vertex, &weight) in outward_load_weights.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 { poll(cx, vertex)?; }
            finite(weight)?;
        }
        let mut rhs = vec![0.0; self.view_factors.len()];
        for (patch, geometry) in self.geometry.patches.iter().enumerate() {
            poll(cx, patch)?;
            let mut absorption = 0.0;
            for (slot, (&vertex, &weight)) in geometry.vertices.iter()
                .zip(&self.absorption_weights[patch]).enumerate()
            {
                if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
                absorption = finite(weight.mul_add(outward_load_weights[vertex], absorption))?;
            }
            for (column, value) in rhs.iter_mut().enumerate() {
                *value = finite(self.view_factors.factors()[patch][column].mul_add(absorption, *value))?;
            }
        }
        let reflected = self.factors.solve_transpose(cx, &rhs)?;
        let mut result = data.template.clone();
        for (slot, face) in data.faces.iter().enumerate() {
            if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
            let reflection_weight = finite(reflected[face.patch] / self.view_factors.areas_m2()[face.patch])?;
            for j in 0..3 {
                let mut derivative = 0.0;
                for i in 0..3 {
                    derivative = finite(face.emissivity[i][j].mul_add(
                        outward_load_weights[face.vertices[i]] - reflection_weight, derivative))?;
                }
                match &mut result.patches[face.patch].emissivity {
                    MaterialEmissivityGradient::Scalar { derivative: value } => {
                        *value = finite(*value + derivative)?;
                    }
                    MaterialEmissivityGradient::Curve { knot_value_derivatives, .. } => {
                        for &(knot, basis) in &face.basis[j] {
                            knot_value_derivatives[knot] = finite(basis.mul_add(derivative, knot_value_derivatives[knot]))?;
                        }
                    }
                }
            }
        }
        if let Some(bath) = &self.report.ambient {
            // The completed enclosure has a black final row: dJ_bath=4 sigma T^3 dT.
            result.reservoir_temperature = Some(finite(-4.0 * SIGMA * bath.temperature_k
                * bath.temperature_k * bath.temperature_k * reflected[self.geometry.patches.len()])?);
        }
        poll(cx, data.faces.len())?;
        Ok(result)
    }
}
