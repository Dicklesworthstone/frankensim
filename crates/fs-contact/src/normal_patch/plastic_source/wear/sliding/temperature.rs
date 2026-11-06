//! Independent material/interface temperatures with original sliding history.
//!
//! The original wear owner re-equilibrates directly under the chosen normal
//! boundary. Original friction claims are then checked at the new interface
//! temperature and each actual loaded pressure. No tangential step is taken:
//! signed slip, removed volume, and both recession fields retain their accepted
//! histories. Thermal evolution and heat partition belong to the enclosing
//! physical owner; normal parameter work is not silently converted into heat.

use std::collections::BTreeMap;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_matdb::MaterialCard;
use fs_material::state_point::ResolvedInterfaceStatePoint;

use crate::normal_patch::{
    FiniteGapPlasticContinuationReceipt, FiniteGapPlasticDrive, FiniteGapPlasticDriveReceipt,
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadReceipt,
};

use super::{
    SourcePlasticArchardWear, SourcePlasticSlidingError, SourcePlasticSlidingWear,
    SourcePlasticSlidingWearState, SourcePlasticWearError, SourcePlasticWearState, checkpoint,
};

/// Original normal response and friction applicability at prescribed temperatures.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticSlidingThermalReceipt<B = ()> {
    /// Complete original normal, wear, friction and registration owner.
    pub original_source_identity: ContentHash,
    /// Continued owner, including the new body and interface source queries.
    pub source_identity: ContentHash,
    /// Contact storage, permanent indentation, normal work and parameter work.
    /// These endpoint accounts imply neither friction work nor a heat partition.
    pub normal: FiniteGapPlasticContinuationReceipt,
    /// Original pinned friction claims at the actual successor loaded pressures.
    /// Identical complete temperature/pressure queries share one sample.
    pub friction_sources: Vec<ResolvedInterfaceStatePoint>,
    /// Per-cell index into `friction_sources`. Detached cells have no sample
    /// and make no claim of friction-source coverage at zero pressure.
    pub cell_source_indices: Vec<Option<usize>>,
    /// Unit at held approach, or the original native load/drive receipt.
    pub boundary: B,
}

/// Complete temperature continuation, published only after every source accepts.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticSlidingThermalRefresh<B = ()> {
    source: SourcePlasticSlidingWear,
    candidate: SourcePlasticSlidingWearState,
    receipt: SourcePlasticSlidingThermalReceipt<B>,
}

impl<B> SourcePlasticSlidingThermalRefresh<B> {
    /// Original selected sources at the independent prescribed temperatures.
    pub const fn source(&self) -> &SourcePlasticSlidingWear {
        &self.source
    }

    /// Continued normal equilibrium and exactly retained native shear histories.
    pub const fn candidate(&self) -> &SourcePlasticSlidingWearState {
        &self.candidate
    }

    /// Native normal accounts and exact current friction-source queries.
    pub const fn receipt(&self) -> &SourcePlasticSlidingThermalReceipt<B> {
        &self.receipt
    }

    /// Adopt the complete continuation without copying the physical field arrays.
    pub fn into_parts(
        self,
    ) -> (
        SourcePlasticSlidingWear,
        SourcePlasticSlidingWearState,
        SourcePlasticSlidingThermalReceipt<B>,
    ) {
        (self.source, self.candidate, self.receipt)
    }
}

impl SourcePlasticSlidingWear {
    /// Continue the original sources at a held normal approach without a shear step.
    ///
    /// Each body's elasticity and indentation/wear hardness follow its own
    /// prescribed temperature. Friction and the two original wear coefficients
    /// follow the independent interface temperature. Original material cards,
    /// claim selections, fixed coordinates, tangent registration, numerical
    /// controls and physical envelopes remain bound to their original owners.
    ///
    /// Reduced hardness may add permanent indentation during normal equilibrium.
    /// No material is removed and no shear distance or friction work is invented.
    /// Native tangential states retain the temperature of their last genuine
    /// shear step; `interface_temperature_k()` reports this new held temperature.
    /// A detached cell does not sample friction at an unsupported zero pressure;
    /// subsequent loaded motion must resolve its actual current source query.
    pub fn refresh_at_temperatures(
        &self,
        before: &SourcePlasticSlidingWearState,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticSlidingThermalRefresh, SourcePlasticSlidingError> {
        self.refresh_temperature_with(before, interface_temperature_k, cx, || {
            self.wear
                .refresh_at_temperatures(
                    &before.wear,
                    material_cards,
                    bulk_temperatures_k,
                    interface_temperature_k,
                    cx,
                )
                .map(|refresh| {
                    let (source, candidate, receipt) = refresh.into_parts();
                    (source, candidate, receipt, ())
                })
        })
    }

    /// Continue temperatures directly at maintained normal force, retaining shear.
    ///
    /// No intermediate held-approach solve can create an artificial dent before
    /// restoring the requested force. Friction is checked at the accepted loaded
    /// pressures, and failure leaves all original normal/wear/shear history intact.
    pub fn refresh_at_temperatures_at_load(
        &self,
        before: &SourcePlasticSlidingWearState,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<
        SourcePlasticSlidingThermalRefresh<FiniteGapPlasticLoadReceipt>,
        SourcePlasticSlidingError,
    > {
        self.refresh_temperature_with(before, interface_temperature_k, cx, || {
            self.wear
                .refresh_at_temperatures_at_load(
                    &before.wear,
                    material_cards,
                    bulk_temperatures_k,
                    interface_temperature_k,
                    target_force_n,
                    controls,
                    cx,
                )
                .map(|refresh| refresh.into_parts())
        })
    }

    /// Continue temperatures against a held loading device with original shear.
    ///
    /// Changed hardness and elasticity can alter approach, pressure and permanent
    /// indentation through the spring. The original native spring/contact account
    /// is retained separately from friction; the held actuator performs no work.
    pub fn refresh_at_temperatures_with_drive(
        &self,
        before: &SourcePlasticSlidingWearState,
        material_cards: [&MaterialCard; 2],
        bulk_temperatures_k: [f64; 2],
        interface_temperature_k: f64,
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<
        SourcePlasticSlidingThermalRefresh<FiniteGapPlasticDriveReceipt>,
        SourcePlasticSlidingError,
    > {
        self.refresh_temperature_with(before, interface_temperature_k, cx, || {
            self.wear
                .refresh_at_temperatures_with_drive(
                    &before.wear,
                    material_cards,
                    bulk_temperatures_k,
                    interface_temperature_k,
                    drive,
                    controls,
                    cx,
                )
                .map(|refresh| refresh.into_parts())
        })
    }

    // Every boundary shares the same original-source admission and atomic
    // history transfer. The normal/wear result stays provisional until all
    // currently loaded cells admit their new friction query.
    fn refresh_temperature_with<B>(
        &self,
        before: &SourcePlasticSlidingWearState,
        interface_temperature_k: f64,
        cx: &Cx<'_>,
        refresh: impl FnOnce() -> Result<
            (
                SourcePlasticArchardWear,
                SourcePlasticWearState,
                FiniteGapPlasticContinuationReceipt,
                B,
            ),
            SourcePlasticWearError,
        >,
    ) -> Result<SourcePlasticSlidingThermalRefresh<B>, SourcePlasticSlidingError> {
        checkpoint(cx)?;
        self.admit(before)?;
        let query = self
            .query_at_temperature(interface_temperature_k)
            .map_err(SourcePlasticSlidingError::Source)?;
        let (wear, wear_state, normal, boundary) = refresh()?;
        Self::admit_query(&wear, &self.friction, &query)?;
        let source = Self {
            identity: Self::source_identity(&wear, &self.friction),
            wear,
            friction: self.friction.clone(),
            initial_law: self.initial_law.clone(),
            temperature_k: interface_temperature_k,
        };
        let pressure = wear_state.contact_state().contact_state().pressure_pa();
        let grid = source.wear.contact().model().config().grid;
        let area = grid.cell_width_m * grid.cell_depth_m;
        let mut friction_sources = Vec::new();
        let mut cell_source_indices = Vec::with_capacity(pressure.len());
        let mut sample_indices = BTreeMap::new();
        for (index, &p) in pressure.iter().enumerate() {
            checkpoint(cx)?;
            let error = |error| SourcePlasticSlidingError::Cell { index, error };
            Self::cell_normal_force(p, area).map_err(error)?;
            let slot = if p == 0.0 {
                None
            } else {
                let key = if source.friction.axes().nominal_pressure.is_some() {
                    p.to_bits()
                } else {
                    0
                };
                let slot = if let Some(&slot) = sample_indices.get(&key) {
                    slot
                } else {
                    let (sample, _) = source.coefficients_at_pressure(p).map_err(error)?;
                    let slot = friction_sources.len();
                    friction_sources.push(sample);
                    sample_indices.insert(key, slot);
                    slot
                };
                Some(slot)
            };
            cell_source_indices.push(slot);
        }
        let receipt = SourcePlasticSlidingThermalReceipt {
            original_source_identity: self.identity,
            source_identity: source.identity,
            normal,
            friction_sources,
            cell_source_indices,
            boundary,
        };
        let candidate = SourcePlasticSlidingWearState {
            owner: source.identity,
            wear: wear_state,
            tangential: before.tangential.clone(),
        };
        checkpoint(cx)?;
        Ok(SourcePlasticSlidingThermalRefresh {
            source,
            candidate,
            receipt,
        })
    }
}
