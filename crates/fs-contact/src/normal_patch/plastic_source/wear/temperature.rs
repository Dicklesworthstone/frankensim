//! Prescribed material temperatures with original loading and wear histories.
//!
//! Normal contact is solved directly under the selected boundary condition.
//! The same original-source refresh retains separate Archard coefficients,
//! hardness tests, per-body removal profiles, and inventories for held-approach,
//! maintained-force, and fixed-device temperature continuation.

use fs_exec::Cx;
use fs_matdb::MaterialCard;
use fs_qty::Temperature;

use crate::normal_patch::{
    FiniteGapPlasticContinuationReceipt, FiniteGapPlasticDrive, FiniteGapPlasticDriveReceipt,
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadReceipt,
};

use super::{
    SourcePlasticArchardWear, SourcePlasticContact, SourcePlasticContactState,
    SourcePlasticWearError, SourcePlasticWearState, SourcePlasticWearThermalRefresh, checkpoint,
    invalid,
};

/// Temperature-dependent contact and wear sources at maintained normal force.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticWearThermalLoadRefresh {
    thermal: SourcePlasticWearThermalRefresh,
    load: FiniteGapPlasticLoadReceipt,
}

impl SourcePlasticWearThermalLoadRefresh {
    /// Original normal, coefficient, and hardness sources at their new temperatures.
    pub const fn source(&self) -> &SourcePlasticArchardWear {
        self.thermal.source()
    }

    /// Continued equilibrium with original plastic history and both wear inventories.
    pub const fn candidate(&self) -> &SourcePlasticWearState {
        self.thermal.candidate()
    }

    /// Elastic storage, plastic work, normal work, and signed parameter work.
    pub const fn receipt(&self) -> &FiniteGapPlasticContinuationReceipt {
        self.thermal.receipt()
    }

    /// Original native load-search residual, bracket, and work count.
    pub const fn load(&self) -> &FiniteGapPlasticLoadReceipt {
        &self.load
    }

    /// Adopt the complete continuation without cloning cell arrays.
    pub fn into_parts(
        self,
    ) -> (
        SourcePlasticArchardWear,
        SourcePlasticWearState,
        FiniteGapPlasticContinuationReceipt,
        FiniteGapPlasticLoadReceipt,
    ) {
        let (source, candidate, receipt) = self.thermal.into_parts();
        (source, candidate, receipt, self.load)
    }
}

/// Temperature-dependent original contact and wear sources at a held device.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticWearThermalDriveRefresh {
    thermal: SourcePlasticWearThermalRefresh,
    drive: FiniteGapPlasticDriveReceipt,
}

impl SourcePlasticWearThermalDriveRefresh {
    /// Original normal, coefficient, and hardness sources at their new temperatures.
    pub const fn source(&self) -> &SourcePlasticArchardWear {
        self.thermal.source()
    }

    /// Continued equilibrium with original indentation and both wear inventories.
    pub const fn candidate(&self) -> &SourcePlasticWearState {
        self.thermal.candidate()
    }

    /// Contact storage, plastic work, spring-to-contact work, and parameter work.
    pub const fn receipt(&self) -> &FiniteGapPlasticContinuationReceipt {
        self.thermal.receipt()
    }

    /// Fixed device, final spring storage, and the native force-balance diagnostics.
    pub const fn drive(&self) -> &FiniteGapPlasticDriveReceipt {
        &self.drive
    }

    /// Adopt the complete original-source continuation without cloning fields.
    pub fn into_parts(
        self,
    ) -> (
        SourcePlasticArchardWear,
        SourcePlasticWearState,
        FiniteGapPlasticContinuationReceipt,
        FiniteGapPlasticDriveReceipt,
    ) {
        let (source, candidate, receipt) = self.thermal.into_parts();
        (source, candidate, receipt, self.drive)
    }
}

impl SourcePlasticArchardWear {
    /// Refresh original temperatures while retaining the accepted normal load.
    ///
    /// Each body's indentation/wear hardness follows that body's temperature;
    /// the original wear coefficients follow the separate interface temperature.
    /// Neither material is replaced, and no removed volume or permanent dent is
    /// reset. The candidate must still fit the remaining homogeneous depth of
    /// each original body after all accepted recession.
    ///
    /// The native solve uses original history directly at maintained force;
    /// there is no held-approach thermal step that could create an artificial
    /// dent before restoring the load. This operation removes no material by
    /// itself and supplies no temperature evolution or heat-partition law.
    pub fn refresh_at_temperatures_at_load(
        &self,
        before: &SourcePlasticWearState,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearThermalLoadRefresh, SourcePlasticWearError> {
        let (thermal, load) =
            self.refresh_temperature_with(before, material_cards, bulk_temperatures_k, cx, || {
                self.contact
                    .refresh_at_temperatures_at_load(
                        &before.contact,
                        material_cards,
                        bulk_temperatures_k,
                        interface_temperature_k,
                        target_force_n,
                        controls,
                        cx,
                    )
                    .map(|refresh| refresh.into_parts())
                    .map_err(SourcePlasticWearError::from)
            })?;
        Ok(SourcePlasticWearThermalLoadRefresh { thermal, load })
    }

    /// Refresh original temperatures while holding the loading device fixed.
    ///
    /// Both original removal profiles and inventories survive, while changed
    /// contact hardness and elasticity can alter force, approach, and permanent
    /// indentation through the device compliance. Each original wear hardness
    /// follows its body's temperature; original coefficients follow the separate
    /// interface temperature. Subsequent wear uses those refreshed sources.
    ///
    /// No material is removed by this temperature operation. The resulting
    /// contact must still fit both remaining homogeneous depths. Normal work
    /// comes from loading-spring storage; the held actuator performs no work.
    /// Thermal evolution, heat partition, and unresolved material paths remain
    /// the enclosing physical owner's responsibility.
    pub fn refresh_at_temperatures_with_drive(
        &self,
        before: &SourcePlasticWearState,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticWearThermalDriveRefresh, SourcePlasticWearError> {
        let (thermal, drive) =
            self.refresh_temperature_with(before, material_cards, bulk_temperatures_k, cx, || {
                self.contact
                    .refresh_at_temperatures_with_drive(
                        &before.contact,
                        material_cards,
                        bulk_temperatures_k,
                        interface_temperature_k,
                        drive,
                        controls,
                        cx,
                    )
                    .map(|refresh| refresh.into_parts())
                    .map_err(SourcePlasticWearError::from)
            })?;
        Ok(SourcePlasticWearThermalDriveRefresh { thermal, drive })
    }

    // Reuse one original-source and inventory transfer for every normal boundary.
    // The callback's candidate is provisional until all refreshed wear claims
    // and remaining-material checks have succeeded.
    pub(super) fn refresh_temperature_with<T>(
        &self,
        before: &SourcePlasticWearState,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        cx: &Cx<'_>,
        refresh_contact: impl FnOnce() -> Result<
            (
                SourcePlasticContact,
                SourcePlasticContactState,
                FiniteGapPlasticContinuationReceipt,
                T,
            ),
            SourcePlasticWearError,
        >,
    ) -> Result<(SourcePlasticWearThermalRefresh, T), SourcePlasticWearError> {
        self.admit(before)?;
        checkpoint(cx)?;
        let owned = self
            .local_axes
            .as_ref()
            .map_or_else(Vec::new, |axes| axes.owned());
        if self.coefficients.properties().iter().any(|property| {
            if self.local_axes.is_none() {
                return !property.has_only_fixed_external_conditions(Some("T"), false);
            }
            !property.has_only_owned_coordinates_and_fixed_external_conditions(&owned)
                || property
                    .source_axes()
                    .iter()
                    .any(|(axis, dims)| axis != "T" && *dims == Temperature::DIMS)
                || property
                    .answer()
                    .evidence
                    .model
                    .validity
                    .bounds()
                    .keys()
                    .any(|axis| {
                        axis != "T"
                            && self
                                .coefficients
                                .axis_quantities()
                                .get(axis)
                                .is_some_and(|quantity| quantity.dims() == Temperature::DIMS)
                    })
        }) || self.hardness.iter().any(|source| {
            source
                .resolved()
                .properties()
                .iter()
                .any(|property| !property.has_only_fixed_external_conditions(Some("T"), false))
        }) {
            return Err(invalid(
                "wear refresh requires owned local coordinates/temperature and otherwise fixed conditions",
            ));
        }
        let (contact, contact_state, receipt, boundary) = refresh_contact()?;
        let normal = contact.normal_source();
        let mut coefficient_query = self.coefficients.query();
        for (name, value) in normal.interface_state().query_point() {
            coefficient_query = coefficient_query
                .with(name, *value)
                .map_err(|_| invalid("finite original wear temperature query"))?;
        }
        let coefficients = if self.coefficients.query() == coefficient_query {
            self.coefficients.clone()
        } else {
            self.coefficients
                .at_query_point(normal.elastic().interface().card(), &coefficient_query)?
        };
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
        let source = Self::build(
            contact,
            coefficients,
            self.coefficient_keys.clone(),
            hardness,
            self.maximum_recession_m,
            self.local_axes.clone(),
        )?;
        source.admit_receded_profile(&contact_state, &before.recession_m, cx)?;
        let candidate = SourcePlasticWearState {
            owner: source.identity,
            contact: contact_state,
            volumes: before.volumes,
            recession_m: before.recession_m.clone(),
        };
        checkpoint(cx)?;
        Ok((
            SourcePlasticWearThermalRefresh {
                source,
                candidate,
                receipt,
            },
            boundary,
        ))
    }
}
