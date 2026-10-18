//! Source-bound electrical and thermal response on the shared enthalpy timeline.
//!
//! A binding uses the timeline's actual caloric law and immutable material card.
//! Both transport properties are resolved together at the evolving temperature
//! and liquid mass fraction, including changes along an isothermal melt plateau.
//! The existing explicit and implicit drivers remain the current/heat owners.

use fs_conduction::transient::deposition::{
    CellId, DepositedCell, DepositionSimulation, MaterialId, ThermalLaw,
};
use fs_matdb::PropertyKey;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

use super::{EnthalpyTransport, JouleEnthalpyError};
use crate::conjugate::ChtError;

#[derive(Debug, Clone, PartialEq)]
enum Source {
    Conducting(ResolvedConductorState),
    Insulating(ResolvedMaterialStatePoint),
}

/// One immutable phase/material evaluation used by both transport systems.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialTransportState {
    phase: EquilibriumPhaseState,
    source: Source,
}

impl MaterialTransportState {
    /// Original caloric-chart state, including phase and both identities.
    #[must_use]
    pub const fn phase(&self) -> EquilibriumPhaseState {
        self.phase
    }

    /// Complete original query, source claims and usage receipts. A derived
    /// electrical resistivity is not inserted as a fictitious measured property.
    #[must_use]
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        match &self.source {
            Source::Conducting(source) => source.material(),
            Source::Insulating(source) => source,
        }
    }

    /// Same resolved bundle lowered by the common scalar Ohmic owner.
    /// None is explicitly declared electrical insulation, never missing data.
    #[must_use]
    pub const fn conductor(&self) -> Option<&ResolvedConductorState> {
        match &self.source {
            Source::Conducting(source) => Some(source),
            Source::Insulating(_) => None,
        }
    }
}

/// An actual thermal cell bound to its original material and caloric law.
///
/// Bind once, then pass `|_, cell, temperature| binding.resolve(cell, temperature)`
/// to the existing enthalpy-transport drivers. For a heterogeneous body, select
/// the binding by stable cell ID. The same callback works for explicit substeps,
/// implicit residual trials, endpoint checks and independent coupled observers.
///
/// Geometry, mass and storage remain owned by the thermal timeline. This is
/// stationary isotropic DC/conduction coupling; it adds no density remapping,
/// pressure work, moving melt, magnetic law or material validity extension.
pub struct MaterialCell<'a> {
    cell: CellId,
    material_id: MaterialId,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    point: ThermalMaterialPoint<'a>,
    conductivity_property: String,
    electrical: Option<ElectricalPropertyBasis>,
    requirements: Vec<ScalarPropertyRequirement>,
    selection: MaterialPropertySelection,
}

impl<'a> MaterialCell<'a> {
    /// Bind the timeline's actual cell/law to one original material card.
    ///
    /// `conductivity` must explicitly identify a scalar static thermal
    /// conductivity. Zero is a supported thermal insulator. `electrical=None`
    /// explicitly declares electrical insulation; otherwise the selected
    /// resistivity or conductivity source is required, even for electrode cells.
    /// `selection` applies to the complete initial bundle, so explicit claim
    /// pins must include both required properties. The selected claims are then
    /// pinned for the binding's lifetime. Crossing a source's support requires
    /// an explicitly admitted new binding, not an implicit source-table switch.
    ///
    /// # Errors
    /// Refuses a missing cell, foreign caloric object/card, wrong conductivity
    /// quantity, or incompatible temperature/phase coordinates. The in-process
    /// law-identity check is not experimental material validation.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        curve: &'a EquilibriumEnthalpyPhaseCurve,
        point: ThermalMaterialPoint<'a>,
        conductivity: &PropertyKey,
        electrical: Option<ElectricalPropertyBasis>,
        selection: MaterialPropertySelection,
    ) -> Result<Self, JouleEnthalpyError> {
        let active = thermal
            .cells()
            .get(&cell)
            .ok_or_else(|| refusal("material binding requires an active thermal cell"))?;
        if !thermal.uses_material_law(active.material_id, curve as &dyn ThermalLaw) {
            return Err(refusal(
                "material binding must use the cell's actual borrowed caloric law",
            ));
        }
        let quantity = conductivity.quantity();
        if quantity != QuantitySpec::dimensional(THERMAL_CONDUCTIVITY_DIMS)
            && quantity
                != QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::ThermalConductivity,
                    ValueForm::Static,
                ))
        {
            return Err(refusal("scalar static thermal conductivity key required"));
        }
        let phase = curve
            .state_at_specific_enthalpy(active.specific_h_j_kg)
            .map_err(refusal)?;
        let initial_point = point.at(phase).map_err(refusal)?;
        let mut requirements = vec![
            ScalarPropertyRequirement::try_with_key(conductivity, ScalarAdmissibility::NonNegative)
                .map_err(refusal)?,
        ];
        if let Some(basis) = electrical {
            requirements.push(basis.requirement().map_err(refusal)?);
        }
        let initial =
            resolve_material_state_point(point.material, &initial_point, &requirements, selection)
                .map_err(refusal)?;
        for property in initial.properties() {
            point.admit_static(property).map_err(refusal)?;
        }
        let selection = MaterialPropertySelection::PinnedByProperty(
            initial
                .properties()
                .iter()
                .map(|property| {
                    (
                        property.requirement().name().to_owned(),
                        property.answer().receipt.selected,
                    )
                })
                .collect(),
        );
        if let Some(basis) = electrical {
            // Admission includes any reciprocal required by the common Ohmic
            // owner. The initial state must not defer an invalid source law
            // until after a coupled advance has already been requested.
            ResolvedConductorState::try_from_material_with_basis(initial, basis)
                .map_err(refusal)?;
        }
        Ok(Self {
            cell,
            material_id: active.material_id,
            curve,
            point,
            conductivity_property: conductivity.name().into(),
            electrical,
            requirements,
            selection,
        })
    }

    /// Original stable thermal-cell assignment.
    #[must_use]
    pub const fn cell_id(&self) -> CellId {
        self.cell
    }

    /// Resolve both transport laws at the SAME current enthalpy-derived state.
    /// External coordinates such as pressure and composition are retained from
    /// the original `ThermalMaterialPoint`. No separate conductivity temperature,
    /// electrical temperature or constant-temperature phase cache exists. Both
    /// selected source identities stay fixed even when another same-name claim
    /// becomes available at a later temperature or phase fraction.
    ///
    /// # Errors
    /// Refuses changed cell/material assignments, foreign temperature samples,
    /// absent/out-of-domain properties or unrepresentable derived electrical
    /// response. The existing driver's transactional acceptance propagates these
    /// failures without publishing a partially heated state.
    pub fn resolve(
        &self,
        cell: &DepositedCell,
        temperature_k: f64,
    ) -> Result<EnthalpyTransport<MaterialTransportState>, JouleEnthalpyError> {
        if cell.id != self.cell || cell.material_id != self.material_id {
            return Err(refusal(
                "material binding does not match the receiving thermal cell",
            ));
        }
        let phase = self
            .curve
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(refusal)?;
        if phase.temperature_k().to_bits() != temperature_k.to_bits() {
            return Err(refusal(
                "transport temperature must come from the same cell and caloric law",
            ));
        }
        let point = self.point.at(phase).map_err(refusal)?;
        let material = resolve_material_state_point(
            self.point.material,
            &point,
            &self.requirements,
            self.selection.clone(),
        )
        .map_err(refusal)?;
        for property in material.properties() {
            self.point.admit_static(property).map_err(refusal)?;
        }
        let conductivity = material
            .property(&self.conductivity_property)
            .ok_or_else(|| refusal("resolved thermal conductivity is absent"))?
            .value_si();
        let (source, resistivity) = if let Some(basis) = self.electrical {
            let conductor = ResolvedConductorState::try_from_material_with_basis(material, basis)
                .map_err(refusal)?;
            let resistivity = conductor.resistivity_ohm_m();
            (Source::Conducting(conductor), Some(resistivity))
        } else {
            (Source::Insulating(material), None)
        };
        EnthalpyTransport::try_new(
            MaterialTransportState { phase, source },
            resistivity,
            conductivity,
        )
    }
}

fn refusal(error: impl core::fmt::Display) -> JouleEnthalpyError {
    ChtError::InvalidInput {
        field: "enthalpy.material",
        reason: error.to_string(),
    }
    .into()
}
