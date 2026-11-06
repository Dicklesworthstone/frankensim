//! Original coefficient paths in the existing finite acoustic preparation.

use super::{
    CaloricReference, ConservativeThermalField, FiniteSolidAdiabaticState, FiniteSolidMaterial,
    FiniteThermoelasticSourceError, ReferenceDeformation, ResolvedHeatCapacityCurve,
    ResolvedNeoHookeanThermalCoefficients, invalid,
};

impl FiniteSolidMaterial {
    /// Prepare entropy-constrained waves with the original Neo-Hookean E(T)
    /// and nu(T) derivatives as well as the actual natural-map derivatives.
    ///
    /// The coefficient source must describe this exact card, selected model,
    /// fixed context and current parameters. Its original source/model owner
    /// is retained alongside the caloric chart. The shared native potential
    /// adds both elastic-coefficient curvature and mixed expansion terms to
    /// c_F exactly once, then applies the original adiabatic acoustic tangent.
    /// No finite-temperature update of the accepted mechanical state is made.
    #[allow(clippy::too_many_arguments)]
    pub fn adiabatic_acoustic_state_with_coefficients<'a>(
        &'a self,
        total_f: &[f64; 9],
        before_thermal: ReferenceDeformation,
        coefficients: &'a ResolvedNeoHookeanThermalCoefficients,
        caloric: &'a ResolvedHeatCapacityCurve,
        specific_h_j_kg: f64,
        caloric_reference: CaloricReference,
    ) -> Result<FiniteSolidAdiabaticState<'a>, FiniteThermoelasticSourceError> {
        self.adiabatic_acoustic_state_inner(
            total_f,
            before_thermal,
            caloric,
            specific_h_j_kg,
            caloric_reference,
            Some(coefficients),
        )
    }
}

impl ConservativeThermalField {
    /// Same original element, manufacturing map and coefficient evidence in
    /// the full compressible local response. This does not infer a constrained
    /// mixed-body acoustic projection or change the retained thermal field.
    #[allow(clippy::too_many_arguments)]
    pub fn adiabatic_acoustic_state_with_coefficients<'a>(
        &'a self,
        element: usize,
        total_f: &[f64; 9],
        coefficients: &'a ResolvedNeoHookeanThermalCoefficients,
        caloric: &'a ResolvedHeatCapacityCurve,
        specific_h_j_kg: f64,
        caloric_reference: CaloricReference,
    ) -> Result<FiniteSolidAdiabaticState<'a>, FiniteThermoelasticSourceError> {
        let material = self
            .material_at(element)
            .ok_or_else(|| invalid("thermoelastic source element is outside the original field"))?;
        let before = self
            .reference_before_thermal_at(element)
            .ok_or_else(|| invalid("the original manufacturing map is absent"))?;
        material.adiabatic_acoustic_state_with_coefficients(
            total_f,
            *before,
            coefficients,
            caloric,
            specific_h_j_kg,
            caloric_reference,
        )
    }
}
