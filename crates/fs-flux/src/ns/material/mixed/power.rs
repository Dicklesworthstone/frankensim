//! One retained reference-density conversion for native mechanical power.

use super::{
    MaterialFluxError, MaterialStokesRheologySolution, MaterialStokesSolution, StokesSolution,
};

/// Mechanical powers ON the fluid per metre of declared extrusion [W/m].
/// Actual prescribed traces and numerical/physical reaction functionals remain
/// native. These are separate source exchanges, not a complete energy balance,
/// an inferred heat source, or a thermodynamic absolute-pressure convention.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialBoundaryPower {
    pub prescribed_velocity_w_per_m: f64,
    pub imposed_traction_w_per_m: f64,
    /// Nonpositive selected backflow stress work; never automatic heating.
    pub backflow_resistance_w_per_m: f64,
    /// Original selected component datums, with mean-zero offsets in Pa.
    pub pressure_offsets_pa: Vec<(usize, Option<f64>)>,
}

pub(super) fn physical_boundary_power(
    flow: &StokesSolution,
    density: f64,
    edges: &[usize],
    pressure_offsets_pa: &[Option<f64>],
) -> Result<MaterialBoundaryPower, MaterialFluxError> {
    let offsets: Vec<_> = pressure_offsets_pa
        .iter()
        .map(|value| value.map(|pressure| pressure / density))
        .collect();
    let power = flow
        .reactions
        .power(edges, &offsets)
        .map_err(MaterialFluxError::Boundary)?;
    let values = [
        power.prescribed_velocity,
        power.imposed_traction,
        power.backflow_resistance,
    ]
    .map(|value| density * value);
    if values.iter().any(|value| !value.is_finite()) {
        return Err(MaterialFluxError::Invalid(
            "unrepresentable physical boundary power",
        ));
    }
    Ok(MaterialBoundaryPower {
        prescribed_velocity_w_per_m: values[0],
        imposed_traction_w_per_m: values[1],
        backflow_resistance_w_per_m: values[2],
        pressure_offsets_pa: power
            .pressure_offsets
            .iter()
            .map(|(component, _)| (*component, pressure_offsets_pa[*component]))
            .collect(),
    })
}

impl MaterialStokesSolution {
    /// Original prescribed/slip, supplied-traction and backflow powers at the
    /// accepted endpoint. The SAME reference rho converts pressure offsets
    /// and powers; no boundary or material callback is repeated. Mean-zero
    /// components require explicit Some(offset_pa), traction datums None.
    ///
    /// # Errors
    /// Preserves native edge/datum refusal and rejects SI conversion overflow.
    pub fn boundary_power(
        &self,
        edges: &[usize],
        pressure_offsets_pa: &[Option<f64>],
    ) -> Result<MaterialBoundaryPower, MaterialFluxError> {
        physical_boundary_power(
            &self.flow,
            self.transport.density_kg_m3(),
            edges,
            pressure_offsets_pa,
        )
    }
}

impl<G> MaterialStokesRheologySolution<G> {
    /// Source-work powers from the accepted nonlinear constitutive endpoint.
    /// Uses the original final sample's reference density without resampling.
    ///
    /// # Errors
    /// Requires an endpoint source and native compatible edges/pressure datums.
    pub fn boundary_power(
        &self,
        edges: &[usize],
        pressure_offsets_pa: &[Option<f64>],
    ) -> Result<MaterialBoundaryPower, MaterialFluxError> {
        let density = self
            .samples
            .first()
            .ok_or(MaterialFluxError::Invalid("missing endpoint material"))?
            .transport()
            .density_kg_m3();
        physical_boundary_power(&self.flow.solution, density, edges, pressure_offsets_pa)
    }
}
