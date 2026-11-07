//! Original H(T, depth), elastic and wear-coefficient temperature continuation.
//!
//! Independent body temperatures evaluate each original complete depth profile;
//! the interface temperature evaluates its original k_A/k_B. The current exposed
//! depths are retained, so heating/cooling cannot restore removed material or
//! reset the relative plastic dent. Native normal re-equilibration owns the
//! endpoint work account at held approach, force or finite loading spring.
//!
//! Temperatures are prescribed by the caller's thermal owner. Parameter work is
//! not heat, and this endpoint update does not add thermal expansion, conduction,
//! annealing kinetics or certify an intermediate temperature path.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_matdb::MaterialCard;
use fs_material::hardness::depth::DepthHardnessError;
use fs_material::state_point::MaterialStatePointError;

use super::{
    DepthHardnessWearError, SourceDepthHardnessWear, Sources, admit_depth_sources, checkpoint,
    depth_source_identity,
};
use crate::normal_patch::adhesive::temperature::{AdhesiveTemperatureError, AdhesiveTemperatures};
use crate::normal_patch::{
    FiniteGapPlasticContinuation, FiniteGapPlasticContinuationReceipt, FiniteGapPlasticDrive,
    FiniteGapPlasticDriveError, FiniteGapPlasticDriveReceipt, FiniteGapPlasticError,
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadError, FiniteGapPlasticLoadReceipt,
    FiniteGapPlasticModel,
};

/// Original source or physical continuation refusal; no partial owner is returned.
#[derive(Clone, Debug, PartialEq)]
pub enum DepthHardnessTemperatureError {
    /// Original bulk/interface temperature source or binding refused.
    NormalSource(AdhesiveTemperatureError),
    /// Original selected interface wear coefficients refused their temperature.
    CoefficientSource(MaterialStatePointError),
    /// Original depth support, native boundary, applicability or inventory refused.
    Wear(DepthHardnessWearError),
    /// A source coordinate is outside this prescribed-temperature adapter.
    Invalid(&'static str),
}
impl From<AdhesiveTemperatureError> for DepthHardnessTemperatureError {
    fn from(error: AdhesiveTemperatureError) -> Self {
        Self::NormalSource(error)
    }
}
impl From<MaterialStatePointError> for DepthHardnessTemperatureError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::CoefficientSource(error)
    }
}
impl From<DepthHardnessWearError> for DepthHardnessTemperatureError {
    fn from(error: DepthHardnessWearError) -> Self {
        Self::Wear(error)
    }
}
impl From<DepthHardnessError> for DepthHardnessTemperatureError {
    fn from(error: DepthHardnessError) -> Self {
        Self::Wear(error.into())
    }
}
impl From<FiniteGapPlasticError> for DepthHardnessTemperatureError {
    fn from(error: FiniteGapPlasticError) -> Self {
        Self::Wear(error.into())
    }
}
impl From<FiniteGapPlasticLoadError> for DepthHardnessTemperatureError {
    fn from(error: FiniteGapPlasticLoadError) -> Self {
        Self::Wear(error.into())
    }
}
impl From<FiniteGapPlasticDriveError> for DepthHardnessTemperatureError {
    fn from(error: FiniteGapPlasticDriveError) -> Self {
        Self::Wear(error.into())
    }
}
impl core::fmt::Display for DepthHardnessTemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "depth-hardness temperature continuation: {self:?}")
    }
}
impl std::error::Error for DepthHardnessTemperatureError {}

/// Complete original-source thermal endpoint and native normal-work account.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthHardnessThermalRefresh {
    /// Original accepted source/configuration identity before temperature change.
    pub original_source_identity: ContentHash,
    /// Complete successor with both original body inventories and exposed depths.
    pub candidate: SourceDepthHardnessWear,
    /// Native storage, irreversible indentation and signed parameter work.
    pub receipt: FiniteGapPlasticContinuationReceipt,
}

/// Thermal continuation at a maintained compressive resultant.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthHardnessThermalLoadRefresh {
    /// Complete source/contact/wear successor and native work account.
    pub thermal: DepthHardnessThermalRefresh,
    /// Original native force equilibrium, approach and search receipt.
    pub load: FiniteGapPlasticLoadReceipt,
}

/// Thermal continuation at a fixed finite loading spring and position.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthHardnessThermalDriveRefresh {
    /// Complete source/contact/wear successor and native work account.
    pub thermal: DepthHardnessThermalRefresh,
    /// Original native spring equilibrium and reversible spring storage.
    pub drive: FiniteGapPlasticDriveReceipt,
}

struct PreparedTemperature {
    sources: Arc<Sources>,
    model: FiniteGapPlasticModel,
}

impl SourceDepthHardnessWear {
    /// Refresh original E(T), nu(T), H(T,z) and k(T) at held approach.
    ///
    /// Each body keeps its own prescribed temperature and actual exposed depth;
    /// interface coefficients follow the independent interface temperature.
    /// Complete prepared depth support is required, including currently unexposed
    /// material. Heating may add indentation; cooling cannot heal a retained dent.
    pub fn refresh_at_temperatures(
        &self,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessThermalRefresh, DepthHardnessTemperatureError> {
        let prepared = self.prepare_temperature_sources(
            material_cards,
            bulk_temperatures_k,
            interface_temperature_k,
            cx,
        )?;
        let next = prepared
            .model
            .refresh_from(&self.model, &self.contact, cx)?;
        self.finish_temperature(prepared, next, cx)
    }

    /// Refresh all original sources directly at the maintained force. The
    /// predecessor must already support that force; a newly exposed/hot profile
    /// unable to support it refuses without accepting an intermediate dent.
    #[allow(clippy::too_many_arguments)]
    pub fn refresh_at_temperatures_at_load(
        &self,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessThermalLoadRefresh, DepthHardnessTemperatureError> {
        self.admit_control_bounds(controls, cx)?;
        let prepared = self.prepare_temperature_sources(
            material_cards,
            bulk_temperatures_k,
            interface_temperature_k,
            cx,
        )?;
        let next = prepared.model.refresh_from_at_load(
            &self.model,
            &self.contact,
            force_n,
            controls,
            cx,
        )?;
        Ok(DepthHardnessThermalLoadRefresh {
            thermal: self.finish_temperature(prepared, next.continuation, cx)?,
            load: next.load,
        })
    }

    /// Refresh original sources while holding the loading spring and position.
    /// Native spring storage supplies normal work; signed parameter work remains
    /// a separate endpoint account, without an invented heat or actuator input.
    #[allow(clippy::too_many_arguments)]
    pub fn refresh_at_temperatures_with_drive(
        &self,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessThermalDriveRefresh, DepthHardnessTemperatureError> {
        self.admit_control_bounds(controls, cx)?;
        let prepared = self.prepare_temperature_sources(
            material_cards,
            bulk_temperatures_k,
            interface_temperature_k,
            cx,
        )?;
        let next = prepared.model.refresh_from_with_drive(
            &self.model,
            &self.contact,
            drive,
            controls,
            cx,
        )?;
        Ok(DepthHardnessThermalDriveRefresh {
            thermal: self.finish_temperature(prepared, next.continuation, cx)?,
            drive: next.drive,
        })
    }

    fn prepare_temperature_sources(
        &self,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        cx: &Cx<'_>,
    ) -> Result<PreparedTemperature, DepthHardnessTemperatureError> {
        checkpoint(cx)?;
        if self
            .sources
            .coefficients
            .properties()
            .iter()
            .any(|property| !property.has_only_fixed_external_conditions(Some("T"), false))
        {
            return Err(DepthHardnessTemperatureError::Invalid(
                "wear temperature continuation requires one original T and otherwise fixed source coordinates",
            ));
        }
        // The original normal owner verifies both immutable cards even for a
        // no-op, retains exact source roles, and changes only canonical T.
        let normal = self.sources.normal.at_temperatures(
            material_cards,
            AdhesiveTemperatures {
                surface_a_k: bulk_temperatures_k[0],
                surface_b_k: bulk_temperatures_k[1],
                interface_k: interface_temperature_k,
            },
        )?;
        let hardness = [
            self.sources.hardness[0].at_temperature(material_cards[0], bulk_temperatures_k[0])?,
            self.sources.hardness[1].at_temperature(material_cards[1], bulk_temperatures_k[1])?,
        ];
        let query = normal.interface_state().query();
        let coefficients = if self.sources.coefficients.query() == query {
            self.sources.coefficients.clone()
        } else {
            self.sources
                .coefficients
                .at_query_point(normal.elastic().interface().card(), &query)?
        };
        checkpoint(cx)?;
        if normal.identity() == self.sources.normal.identity()
            && coefficients.identity() == self.sources.coefficients.identity()
            && hardness
                .iter()
                .zip(&self.sources.hardness)
                .all(|(next, old)| next.identity() == old.identity())
        {
            return Ok(PreparedTemperature {
                sources: Arc::clone(&self.sources),
                model: self.model.clone(),
            });
        }
        let (k, selected) = admit_depth_sources(
            &normal,
            &coefficients,
            &self.sources.coefficient_keys,
            &hardness,
            self.sources.envelope,
            self.sources.maximum_recession_m,
        )?;
        let mut caps = Vec::with_capacity(self.contact.pressure_pa().len());
        for cell in 0..self.contact.pressure_pa().len() {
            if cell % 64 == 0 {
                checkpoint(cx)?;
            }
            let a = hardness[0].hardness_pa_at(self.recession_m[0][cell])?;
            let b = hardness[1].hardness_pa_at(self.recession_m[1][cell])?;
            caps.push(a.min(b));
        }
        let model =
            self.model
                .with_material_properties(normal.elastic().reduced_modulus_pa(), caps, cx)?;
        let identity = depth_source_identity(
            &normal,
            &coefficients,
            &hardness,
            &selected,
            &model,
            self.sources.envelope,
            self.sources.maximum_recession_m,
        );
        let sources = Arc::new(Sources {
            normal,
            coefficients,
            hardness,
            k,
            identity,
            coefficient_keys: self.sources.coefficient_keys.clone(),
            envelope: self.sources.envelope,
            maximum_recession_m: self.sources.maximum_recession_m,
        });
        checkpoint(cx)?;
        Ok(PreparedTemperature { sources, model })
    }

    fn finish_temperature(
        &self,
        prepared: PreparedTemperature,
        next: FiniteGapPlasticContinuation,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessThermalRefresh, DepthHardnessTemperatureError> {
        self.admit(&prepared.model, &next.candidate, &self.recession_m, cx)?;
        checkpoint(cx)?;
        Ok(DepthHardnessThermalRefresh {
            original_source_identity: self.source_identity(),
            candidate: Self {
                sources: prepared.sources,
                model: prepared.model,
                contact: next.candidate,
                volumes: self.volumes,
                recession_m: self.recession_m.clone(),
            },
            receipt: next.receipt,
        })
    }
}
