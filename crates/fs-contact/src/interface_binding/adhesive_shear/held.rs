//! Retain both original source owners around a frozen normal patch.
use super::*;
use crate::normal_patch::maugis::MaugisResponse;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};
use fs_tribo::adhesive::AdhesiveTangentialContact;
use fs_tribo::{DryTangentialState, TangentialLaw, TangentialStep, TangentialStepInput};

/// A frozen area/load with the original ordered normal and shear source states.
/// The whole source query remains fixed. This does not infer a thermal refresh,
/// normal evolution, shear-induced peeling or a new material calibration.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceAdhesiveTangential {
    normal: BoundNormalInterfaceState,
    shear: SourceAdhesiveShear,
    tangential: AdhesiveTangentialContact,
    temperature_kelvin: Option<f64>,
}
impl SourceAdhesiveShear {
    /// Freeze actual JKR contact geometry, including finite tensile-load area.
    pub fn bind_jkr(
        &self,
        normal: &BoundNormalInterfaceState,
        response: &JkrResponse,
    ) -> Result<SourceAdhesiveTangential, AdhesiveShearError> {
        self.bind_patch(
            normal,
            response.normal_state_identity(),
            response.friction_contact_state()?,
        )
    }
    /// Freeze actual DMT solid area; noncontact attraction carries no shear.
    pub fn bind_dmt(
        &self,
        normal: &BoundNormalInterfaceState,
        response: &DmtResponse,
    ) -> Result<SourceAdhesiveTangential, AdhesiveShearError> {
        self.bind_patch(
            normal,
            response.normal_state_identity(),
            response.friction_contact_state()?,
        )
    }
    /// Freeze Maugis real contact area, excluding its noncontact cohesive annulus.
    pub fn bind_maugis(
        &self,
        normal: &BoundNormalInterfaceState,
        response: &MaugisResponse,
    ) -> Result<SourceAdhesiveTangential, AdhesiveShearError> {
        self.bind_patch(
            normal,
            response.normal_state_identity(),
            response.friction_contact_state()?,
        )
    }
    fn bind_patch(
        &self,
        normal: &BoundNormalInterfaceState,
        response_identity: ContentHash,
        contact: AdhesiveContactState,
    ) -> Result<SourceAdhesiveTangential, AdhesiveShearError> {
        if normal.identity() != self.normal_identity || response_identity != self.normal_identity {
            return Err(AdhesiveShearError::Invalid {
                reason: "held adhesive patch must retain the original normal and shear source state",
            });
        }
        let mut temperature_kelvin = None;
        for (axis, value) in self.resolved.query_point() {
            let quantity = self.resolved.axis_quantities().get(axis).copied();
            let temperature = quantity
                == Some(QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )))
                || axis == "T"
                    && (quantity.is_none()
                        || quantity == Some(QuantitySpec::dimensional(Temperature::DIMS)));
            if temperature {
                if temperature_kelvin.is_some_and(|old| old != *value) {
                    return Err(AdhesiveShearError::Invalid {
                        reason: "held adhesive interface needs one common fixed contact temperature",
                    });
                }
                temperature_kelvin = Some(*value);
            }
        }
        Ok(SourceAdhesiveTangential {
            normal: normal.clone(),
            shear: self.clone(),
            tangential: AdhesiveTangentialContact::new(self.law.clone(), contact)?,
            temperature_kelvin,
        })
    }
    /// Evaluate the instantaneous Maugis solid-contact junction shear.
    pub fn evaluate_maugis<'a>(
        &self,
        normal: &'a BoundNormalInterfaceState,
        response: &MaugisResponse,
        slip: TangentialSlip,
    ) -> Result<AdhesiveFrictionResponse<'a>, AdhesiveShearError> {
        if response.normal_state_identity() != self.normal_identity {
            return Err(AdhesiveShearError::Invalid {
                reason: "Maugis response belongs to a different normal source state",
            });
        }
        self.evaluate(normal, response.friction_contact_state()?, slip)
    }
}
impl SourceAdhesiveTangential {
    /// Complete original normal source and its exact frozen external query.
    #[must_use]
    pub const fn normal_state(&self) -> &BoundNormalInterfaceState {
        &self.normal
    }
    /// Complete original shear claims and their usage receipts.
    #[must_use]
    pub const fn shear_source(&self) -> &SourceAdhesiveShear {
        &self.shear
    }
    /// Constitutive leaf with real area and repulsion, not signed net load.
    #[must_use]
    pub const fn tangential(&self) -> &AdhesiveTangentialContact {
        &self.tangential
    }
    /// Original common contact temperature, preserving absence when unprovided.
    #[must_use]
    pub const fn temperature_kelvin(&self) -> Option<f64> {
        self.temperature_kelvin
    }
}
impl TangentialLaw for SourceAdhesiveTangential {
    type State = DryTangentialState;
    fn model_id(&self) -> &str {
        "source-adhesive-junction-held-effort"
    }
    fn initial_state(&self) -> Self::State {
        self.tangential.initial_state()
    }
    fn stored_energy_j(&self, s: &Self::State) -> f64 {
        self.tangential.stored_energy_j(s)
    }
    fn advance(
        &self,
        s: &Self::State,
        x: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        if x.contact_temperature_k != self.temperature_kelvin {
            return Err(TriboError::InvalidInput {
                field: "held adhesive temperature must match the original frozen source query",
            });
        }
        if x.nominal_pressure_pa.is_some() {
            return Err(TriboError::InvalidInput {
                field: "held adhesive patch does not own a nominal-pressure source refresh",
            });
        }
        self.tangential.advance(s, x)
    }
}
