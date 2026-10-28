//! Original two-solid source compliance in the finite-caloric normal owner.
//! Both bulk states and caloric charts retain the same immutable cards, fixed
//! external context and exact selected density source. Temperatures are never
//! replaced by a representative interface temperature.

use super::{NormalModulus, NormalModulusState, NormalThermalError, modulus_and_damping};
use fs_contact::normal_patch::NormalPatchRequest;
use fs_material::elastic::thermal::contact::{ThermalHertzError, ThermalHertzPair};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, ResolvedHeatCapacityCurve};
use fs_material::state_point::DENSITY_PROPERTY;

/// Source-backed fixed-geometry Hertz/Hunt--Crossley free energy. The source
/// pair computes the nonlinear reduced modulus and its reciprocal derivatives.
/// Density does not replace the accepted reference geometry or mass.
#[derive(Clone)]
pub struct MaterialNormalModulus<'a> {
    reference: NormalPatchRequest,
    source: ThermalHertzPair,
    caloric: [&'a ResolvedHeatCapacityCurve; 2],
}
fn source_error(error: ThermalHertzError) -> NormalThermalError {
    NormalThermalError::MaterialModel(Box::new(error))
}
impl<'a> MaterialNormalModulus<'a> {
    pub fn new(
        reference: NormalPatchRequest,
        source: ThermalHertzPair,
        caloric: [&'a ResolvedHeatCapacityCurve; 2],
    ) -> Result<Self, NormalThermalError> {
        let invalid = NormalThermalError::Invalid;
        let native = source
            .at_temperatures(source.reference_temperatures_k())
            .map_err(source_error)?;
        let (modulus, _) = modulus_and_damping(&reference.law)?;
        if modulus.to_bits() != native.reduced_modulus.value_pa.to_bits()
            || reference.applicability.yield_strength_pa.to_bits()
                != source.limiting_yield_stress_pa().to_bits()
        {
            return Err(invalid(
                "normal request must use the original pair modulus and limiting yield stress",
            ));
        }
        for i in 0..2 {
            let range = source.temperature_ranges_k()[i];
            if range[0] < reference.limits.min_temperature_k
                || range[1] > reference.limits.max_temperature_k
                || caloric[i].curve().material_card_identity() != source.card_identities()[i]
            {
                return Err(invalid(
                    "normal source domain or ordered caloric material mismatch",
                ));
            }
            let original = native.surfaces[i].resolved();
            let density = original
                .property(DENSITY_PROPERTY)
                .ok_or(invalid("original solid density source is missing"))?;
            for receipt in caloric[i].receipts() {
                if receipt.axis_quantities != *original.axis_quantities()
                    || receipt
                        .query_point
                        .iter()
                        .filter(|(axis, _)| axis != "T")
                        .ne(original
                            .query_point()
                            .iter()
                            .filter(|(axis, _)| axis != "T"))
                    || (receipt.property == DENSITY_PROPERTY
                        && receipt.selected != density.answer().receipt.selected)
                {
                    return Err(invalid(
                        "caloric and elastic sources must share complete external conditions and original density claim",
                    ));
                }
            }
        }
        Ok(Self {
            reference,
            source,
            caloric,
        })
    }
    pub const fn source(&self) -> &ThermalHertzPair {
        &self.source
    }
}
impl NormalModulus for MaterialNormalModulus<'_> {
    fn reference_request(&self) -> &NormalPatchRequest {
        &self.reference
    }
    fn evaluate(&self, temperatures_k: [f64; 2]) -> Result<NormalModulusState, NormalThermalError> {
        let result = self
            .source
            .at_temperatures(temperatures_k)
            .map_err(source_error)?;
        Ok(NormalModulusState {
            modulus_pa: result.reduced_modulus.value_pa,
            gradient_pa_k: result.gradient_pa_k,
            hessian_pa_k2: result.hessian_pa_k2,
        })
    }
    fn expected_caloric_curves(&self) -> Option<[&EquilibriumEnthalpyPhaseCurve; 2]> {
        Some(self.caloric.map(|source| source.curve()))
    }
}
