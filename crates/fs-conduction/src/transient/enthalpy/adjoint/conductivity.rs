//! Conductivity scales of the accepted enthalpy endpoint.
//! The retained local loads include its actual phase multiplier, but no
//! boundary, contact or source terms and no second state derivative.

use crate::assemble::{element_stiffness, element_temperature};

use super::{
    ASSEMBLY_TILE, ConductionProblem, Cx, EnthalpyAdjointError,
    EnthalpyStepLinearization, finite, poll, vector,
};

/// Reuse the production sampling and stiffness scaling order. Phase scales
/// come only from StepContext's typed laws at the accepted nodal enthalpies.
pub(super) fn element_loads(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    temperature: &[f64],
    phase_scales: Option<&[f64]>,
) -> Result<Vec<[f64; 4]>, EnthalpyAdjointError> {
    poll(cx, 0)?;
    let mesh = problem.mesh;
    let elements = mesh.element_count();
    if phase_scales.is_some_and(|scales| scales.len() != elements)
        || elements.checked_mul(std::mem::size_of::<[f64; 4]>())
            .is_none_or(|bytes| isize::try_from(bytes).is_err())
    {
        return Err(EnthalpyAdjointError::InvalidInput(
            "enthalpy conductivity load shape or allocation size differs",
        ));
    }
    vector(cx, temperature, mesh.vertex_count())?;
    let mut loads = Vec::new();
    loads.try_reserve_exact(elements).map_err(|_| {
        EnthalpyAdjointError::InvalidInput("enthalpy conductivity load allocation refused")
    })?;
    for (element, tet) in mesh.complex().tets.iter().enumerate() {
        if element % ASSEMBLY_TILE == 0 { poll(cx, element)?; }
        let model = match problem.element_materials {
            Some(materials) => materials.model_for(element)?,
            None => problem.material,
        };
        let sample = finite(element_temperature(mesh, element, temperature))?;
        let tensor = model.tensor_at(sample)?;
        let mut stiffness = element_stiffness(mesh, element, &tensor);
        if let Some(scales) = phase_scales {
            let scale = finite(scales[element])?;
            if scale <= 0.0 {
                return Err(EnthalpyAdjointError::InvalidInput(
                    "enthalpy conductivity requires positive phase scales",
                ));
            }
            for row in &mut stiffness {
                for entry in row { *entry = finite(*entry * scale)?; }
            }
        }
        let mut local = [0.0; 4];
        for (a, row) in stiffness.iter().enumerate() {
            for (b, &vertex) in tet.iter().enumerate() {
                local[a] = finite(finite(row[b])?.mul_add(
                    temperature[vertex as usize], local[a],
                ))?;
            }
        }
        loads.push(local);
    }
    poll(cx, elements)?;
    Ok(loads)
}

impl EnthalpyStepLinearization<'_> {
    /// Pull the same complete joule-residual multiplier back to element scales.
    ///
    /// Element e has control K_eff,e(T,f) -> s_e K_eff,e(T,f), with s_e=1.
    /// Return dJ/ds_e at s_e=1, equivalently dJ/dln(s_e) at ln(s_e)=0,
    /// in mesh element order.
    /// K_eff includes the assigned full tensor/temperature curve and the
    /// accepted typed phase multiplier. Sum element entries for a shared
    /// regional scale without dividing by volume.
    ///
    /// The caller supplies the unchanged complete multiplier from the checked
    /// enthalpy transpose, including any natural, air and radiation feedback.
    /// The result is -dt * lambda_e^T * (K_eff,e*T_e). Neither k'(T) nor the
    /// phase-state chain rule is inserted again; both belong to that transpose.
    /// Constitutive card identity, table knots, phase-law settings, enthalpy
    /// chart, storage mass, contact, source, geometry and time grid stay fixed.
    /// This is a full-law scale, not a derivative of material selection or an
    /// individual tensor entry. Explicit objective dependence must be added
    /// separately. No new state, transpose solve or material query is performed.
    ///
    /// The owner retains 32 bytes per element. This call allocates only one
    /// 8-byte-per-element result vector, after checking max_elements and size.
    ///
    /// # Errors
    /// Wrong full-dual length, nonfinite arithmetic, exhausted element or
    /// allocation allowance, and cancellation publish no partial gradient.
    pub fn conductivity_scale_pullback(
        &self,
        cx: &Cx<'_>,
        full_joule_multiplier: &[f64],
        max_elements: usize,
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        poll(cx, 0)?;
        let elements = self.mesh.element_count();
        if elements > max_elements || self.conductive_loads.len() != elements
            || elements.checked_mul(std::mem::size_of::<f64>())
                .is_none_or(|bytes| isize::try_from(bytes).is_err())
        {
            return Err(EnthalpyAdjointError::InvalidInput(
                "enthalpy conductivity pullback exceeds its element allowance",
            ));
        }
        vector(cx, full_joule_multiplier, self.mesh.vertex_count())?;
        let mut result = Vec::new();
        result.try_reserve_exact(elements).map_err(|_| {
            EnthalpyAdjointError::InvalidInput("enthalpy conductivity gradient allocation refused")
        })?;
        for (element, (tet, loads)) in self.mesh.complex().tets.iter()
            .zip(&self.conductive_loads).enumerate()
        {
            if element % ASSEMBLY_TILE == 0 { poll(cx, element)?; }
            let mut contraction = 0.0;
            for (&vertex, &load) in tet.iter().zip(loads) {
                if self.is_prescribed(vertex as usize) { continue; }
                contraction = finite(
                    full_joule_multiplier[vertex as usize].mul_add(load, contraction),
                )?;
            }
            result.push(finite(-self.dt * contraction)?);
        }
        poll(cx, elements)?;
        Ok(result)
    }
}
