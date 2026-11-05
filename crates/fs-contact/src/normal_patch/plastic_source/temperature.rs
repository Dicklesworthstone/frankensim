//! Prescribed material temperatures with retained indentation and wear history.
//!
//! The chosen normal boundary, original sampled geometry, and applicability
//! envelope stay fixed. Original selected material claims supply E(T), nu(T), and H(T)
//! independently for the two bodies. The pressure field re-equilibrates and a
//! reduced hardness may produce additional permanent indentation. Cooling does
//! not heal an existing dent or replace removed material.
//!
//! At held approach the signed parameter supply is the elastic-energy change
//! plus incremental plastic work. At maintained force it additionally subtracts
//! normal work `F * delta(approach)`. This is an endpoint account for the prescribed
//! condition change, not heat, entropy, thermal expansion, or a continuous
//! thermal-path integral.

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_matdb::MaterialCard;
use fs_material::state_point::MaterialStatePointError;

use crate::normal_patch::adhesive::temperature::{AdhesiveTemperatureError, AdhesiveTemperatures};
use crate::normal_patch::{
    FiniteGapContactError, FiniteGapPlasticContinuationReceipt, FiniteGapPlasticError,
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadError, FiniteGapPlasticLoadReceipt,
};

use super::{
    PlasticContactControls, SourcePlasticContact, SourcePlasticContactError,
    SourcePlasticContactState,
};

/// Exact material, interface, contact, or applicability refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum SourcePlasticTemperatureError {
    /// Original bulk/interface temperature resolution or binding refused.
    NormalSource(AdhesiveTemperatureError),
    /// Original selected hardness bundle refused its temperature query.
    HardnessSource(MaterialStatePointError),
    /// Native history transfer or the original physical envelope refused.
    Contact(SourcePlasticContactError),
    /// The maintained normal force or native load search refused.
    Load(FiniteGapPlasticLoadError),
}

impl From<AdhesiveTemperatureError> for SourcePlasticTemperatureError {
    fn from(error: AdhesiveTemperatureError) -> Self {
        Self::NormalSource(error)
    }
}

impl From<MaterialStatePointError> for SourcePlasticTemperatureError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::HardnessSource(error)
    }
}

impl From<SourcePlasticContactError> for SourcePlasticTemperatureError {
    fn from(error: SourcePlasticContactError) -> Self {
        Self::Contact(error)
    }
}

impl From<FiniteGapPlasticError> for SourcePlasticTemperatureError {
    fn from(error: FiniteGapPlasticError) -> Self {
        Self::Contact(error.into())
    }
}

impl From<FiniteGapPlasticLoadError> for SourcePlasticTemperatureError {
    fn from(error: FiniteGapPlasticLoadError) -> Self {
        Self::Load(error)
    }
}

impl core::fmt::Display for SourcePlasticTemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "source plastic-contact temperature continuation: {self:?}"
        )
    }
}

impl std::error::Error for SourcePlasticTemperatureError {}

/// Refreshed original sources and their immutable continued physical history.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticThermalRefresh {
    original_source_identity: ContentHash,
    source: SourcePlasticContact,
    candidate: SourcePlasticContactState,
    receipt: FiniteGapPlasticContinuationReceipt,
}

impl SourcePlasticThermalRefresh {
    /// Complete original source owner, before the prescribed temperature change.
    pub const fn original_source_identity(&self) -> ContentHash {
        self.original_source_identity
    }

    /// Refreshed original material/interface sources and numerical owner.
    pub const fn source(&self) -> &SourcePlasticContact {
        &self.source
    }

    /// Retained indentation and recession, equilibrated at the held approach.
    pub const fn candidate(&self) -> &SourcePlasticContactState {
        &self.candidate
    }

    /// Original native endpoint response and signed parameter-work account.
    pub const fn receipt(&self) -> &FiniteGapPlasticContinuationReceipt {
        &self.receipt
    }

    /// Move the complete accepted continuation without copying its field arrays.
    pub fn into_parts(
        self,
    ) -> (
        SourcePlasticContact,
        SourcePlasticContactState,
        FiniteGapPlasticContinuationReceipt,
    ) {
        (self.source, self.candidate, self.receipt)
    }
}

/// Original-source temperature continuation with a maintained normal resultant.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticThermalLoadRefresh {
    thermal: SourcePlasticThermalRefresh,
    load: FiniteGapPlasticLoadReceipt,
}

impl SourcePlasticThermalLoadRefresh {
    /// Complete original source owner before changing prescribed temperatures.
    pub const fn original_source_identity(&self) -> ContentHash {
        self.thermal.original_source_identity()
    }

    /// Original selected sources, evaluated at the new independent temperatures.
    pub const fn source(&self) -> &SourcePlasticContact {
        self.thermal.source()
    }

    /// Direct maintained-force equilibrium with retained original physical history.
    pub const fn candidate(&self) -> &SourcePlasticContactState {
        self.thermal.candidate()
    }

    /// Storage, irreversible indentation, normal work, and parameter work.
    pub const fn receipt(&self) -> &FiniteGapPlasticContinuationReceipt {
        self.thermal.receipt()
    }

    /// Native force residual, approach bracket, and bounded work count.
    pub const fn load(&self) -> &FiniteGapPlasticLoadReceipt {
        &self.load
    }

    /// Move the complete original-source continuation without copying fields.
    pub fn into_parts(
        self,
    ) -> (
        SourcePlasticContact,
        SourcePlasticContactState,
        FiniteGapPlasticContinuationReceipt,
        FiniteGapPlasticLoadReceipt,
    ) {
        let (source, candidate, receipt) = self.thermal.into_parts();
        (source, candidate, receipt, self.load)
    }
}

impl SourcePlasticContact {
    /// Change prescribed body and interface temperatures at fixed rigid approach.
    ///
    /// The two immutable material cards must be the original ordered cards,
    /// including for a no-op temperature change. Only the canonical absolute T
    /// coordinate changes. Every selected claim, hardness test, and other fixed
    /// condition remains original. All thermal properties are queried before a
    /// candidate is produced, and the original geometry/depth limits still apply.
    ///
    /// This endpoint return mapping may add plastic work when hardness falls.
    /// Its parameter-work account requires refinement along changing conditions;
    /// it does not provide a caloric or autonomous temperature evolution law.
    pub fn refresh_at_temperatures(
        &self,
        before: &SourcePlasticContactState,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticThermalRefresh, SourcePlasticTemperatureError> {
        let source = self.resolve_temperature_sources(
            before,
            material_cards,
            bulk_temperatures_k,
            interface_temperature_k,
            cx,
        )?;
        let continuation = source
            .model
            .refresh_from(&self.model, &before.contact, cx)?;
        source.admit_candidate(&continuation.candidate, cx)?;
        Ok(SourcePlasticThermalRefresh {
            original_source_identity: self.identity,
            candidate: SourcePlasticContactState {
                owner: source.identity,
                contact: continuation.candidate,
            },
            source,
            receipt: continuation.receipt,
        })
    }

    /// Change original prescribed temperatures while retaining the applied force.
    ///
    /// The accepted predecessor must already carry `target_force_n` within the
    /// declared force tolerance. Temperature-dependent E, nu, and hardness are
    /// resolved once, then the native solver finds the new approach directly
    /// from the ORIGINAL indentation/recession history. A held-approach thermal
    /// return map is never used as an intermediate; it could create a spurious
    /// irreversible dent before the force is restored.
    ///
    /// The whole approach bracket must fit the original physical envelope. A
    /// softened surface that cannot sustain the force, or a force on an
    /// unresolved plastic plateau, is explicitly refused. The returned normal
    /// work is `F * delta(approach)`; signed parameter work is an endpoint account
    /// requiring path refinement, not heat or an autonomous temperature law.
    pub fn refresh_at_temperatures_at_load(
        &self,
        before: &SourcePlasticContactState,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticThermalLoadRefresh, SourcePlasticTemperatureError> {
        self.admit_load_controls(controls)?;
        let source = self.resolve_temperature_sources(
            before,
            material_cards,
            bulk_temperatures_k,
            interface_temperature_k,
            cx,
        )?;
        let next = source.model.refresh_from_at_load(
            &self.model,
            &before.contact,
            target_force_n,
            controls,
            cx,
        )?;
        source.admit_candidate(&next.continuation.candidate, cx)?;
        Ok(SourcePlasticThermalLoadRefresh {
            thermal: SourcePlasticThermalRefresh {
                original_source_identity: self.identity,
                candidate: SourcePlasticContactState {
                    owner: source.identity,
                    contact: next.continuation.candidate,
                },
                source,
                receipt: next.continuation.receipt,
            },
            load: next.load,
        })
    }

    // Source selection is independent of the normal boundary condition. Both
    // continuations retain these exact claims, test contexts, and fixed axes.
    fn resolve_temperature_sources(
        &self,
        before: &SourcePlasticContactState,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, SourcePlasticTemperatureError> {
        cx.checkpoint()
            .map_err(|_| FiniteGapPlasticError::Contact(FiniteGapContactError::Cancelled))?;
        if before.owner != self.identity {
            return Err(SourcePlasticContactError::Invalid(
                "thermal continuation requires the original source-bound history",
            )
            .into());
        }
        if self.hardness.iter().any(|source| {
            source
                .resolved()
                .properties()
                .iter()
                .any(|property| !property.has_only_fixed_external_conditions(Some("T"), false))
        }) {
            return Err(SourcePlasticContactError::Invalid(
                "hardness refresh requires one owned temperature and otherwise fixed conditions",
            )
            .into());
        }
        let normal = self.normal.at_temperatures(
            material_cards,
            AdhesiveTemperatures {
                surface_a_k: bulk_temperatures_k[0],
                surface_b_k: bulk_temperatures_k[1],
                interface_k: interface_temperature_k,
            },
        )?;
        let points = [
            normal.elastic().surface_a_state().resolved().query(),
            normal.elastic().surface_b_state().resolved().query(),
        ];
        let mut hardness = self.hardness.clone();
        for side in 0..2 {
            if self.hardness[side].resolved().state_coordinate("T")
                != Some(bulk_temperatures_k[side])
            {
                hardness[side] =
                    self.hardness[side].at_query_point(material_cards[side], &points[side])?;
            }
        }
        let unchanged = normal.identity() == self.normal.identity()
            && hardness[0].identity() == self.hardness[0].identity()
            && hardness[1].identity() == self.hardness[1].identity();
        if unchanged {
            Ok(self.clone())
        } else {
            let config = self.model.config();
            Ok(Self::new(
                normal,
                hardness,
                config.grid,
                config.undeformed_gap_m.clone(),
                PlasticContactControls {
                    maximum_active_set_iterations: config.maximum_active_set_iterations,
                    complementarity_tolerance_m: config.complementarity_tolerance_m,
                    boundary_clearance_cells: config.boundary_clearance_cells,
                },
                self.envelope,
                cx,
            )?)
        }
    }
}
