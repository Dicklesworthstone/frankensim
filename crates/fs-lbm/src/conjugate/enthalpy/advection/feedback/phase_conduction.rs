//! The original h chart drives conductivity on the owned native fluid faces.

use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::phase::conductivity::LiquidMassFractionConductivity;

use super::*;

#[derive(Clone, Copy)]
pub(super) struct PhaseConduction<'a> {
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    law: LiquidMassFractionConductivity,
}

impl<'a, F: FluidTransport + Clone, R, O: DepositionObserver> MaterialFlowFeedback<'a, F, R, O>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Scale owned fluid-face conduction with the shared declared phase law.
    ///
    /// The exact registered h chart supplies liquid mass fraction at each
    /// applied source epoch, including along an isothermal latent plateau.
    /// Each half-cell's existing source-backed k(T) is scaled before the
    /// original series-resistance reduction. The same resulting conductance
    /// drives heat exchange and the original conductive/advective CFL.
    /// Authored contacts, boundary exchange, viscosity and caloric storage
    /// remain with their existing owners.
    ///
    /// This is an explicit phenomenological mass-fraction multiplier, not an
    /// inferred volume-mixture rule or measured solid conductivity. The caller
    /// must admit its base source and multipliers for the declared composition,
    /// pressure and phase regime. Constant reference density is required; no
    /// shrinkage source, solid locking, thermal solver or latent heat is added.
    /// Attach `Self::with_mushy_mobility` separately for phase-dependent drag.
    ///
    /// This consumed builder replaces any previous multiplier assignment.
    /// Source/endpoint refusals use the existing all-or-nothing advance.
    pub fn with_phase_conductivity(
        mut self,
        curve: &'a EquilibriumEnthalpyPhaseCurve,
        law: LiquidMassFractionConductivity,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(self.gate, &mut || Ok(()))?;
        if !core::ptr::addr_eq(
            self.input.caloric as *const dyn fs_conduction::transient::deposition::ThermalLaw,
            curve as *const EquilibriumEnthalpyPhaseCurve,
        ) {
            return Err(invalid("phase conduction must use the original caloric object").into());
        }
        for knot in curve.knots() {
            check_gate(self.gate, &mut || Ok(()))?;
            if knot.bulk_density_kg_m3.to_bits() != self.reference.density_kg_m3().to_bits() {
                return Err(invalid(
                    "phase conduction requires the existing fixed reference density",
                )
                .into());
            }
        }
        self.phase_conduction = Some(PhaseConduction { curve, law });
        self.check_phase_conduction(&self.endpoint, &mut || Ok(()))?;
        Ok(self)
    }

    /// Declared multiplier actually used with the retained source samples.
    #[must_use]
    pub fn phase_conductivity_law(&self) -> Option<LiquidMassFractionConductivity> {
        self.phase_conduction.map(|binding| binding.law)
    }

    pub(super) fn effective_conductivity(
        &self,
        sample: &FlowMaterialSample<F>,
    ) -> Result<f64, MaterialFlowError<O::Error>> {
        let base = sample.transport.thermal_conductivity_w_m_k();
        let Some(binding) = self.phase_conduction else {
            return Ok(base);
        };
        let material_error = |what| MaterialFlowError::Material {
            cell: sample.cell,
            source: MaterialError::State { what },
        };
        let state = binding
            .curve
            .state_at_specific_enthalpy(sample.specific_enthalpy_j_kg)
            .map_err(|error| material_error(error.to_string()))?;
        if state.temperature_k().to_bits() != sample.temperature_k.to_bits() {
            return Err(invalid(
                "phase conductivity temperature differs from the original caloric owner",
            )
            .into());
        }
        let multiplier = binding
            .law
            .multiplier_at(state.liquid_mass_fraction())
            .map_err(|error| material_error(error.to_string()))?;
        let effective = base * multiplier;
        if !effective.is_finite() || effective <= 0.0 {
            return Err(invalid("phase-scaled fluid conductivity is not representable").into());
        }
        Ok(effective)
    }

    pub(super) fn check_phase_conduction(
        &self,
        samples: &[FlowMaterialSample<F>],
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), MaterialFlowError<O::Error>> {
        if self.phase_conduction.is_some() {
            for sample in samples {
                check_gate(self.gate, check)?;
                self.effective_conductivity(sample)?;
            }
        }
        Ok(())
    }
}
