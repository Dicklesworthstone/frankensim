//! Explicit native destination for the original passive-load acoustic work.
//! The load's impedance already specifies its physical area and mode. This
//! binding supplies only the separate thermal control volume; it infers no
//! material, heat fraction, geometry, or temperature dependence of impedance.
//! Binding declares that the load's full real-impedance work thermalizes in
//! that body. A radiation load's outgoing energy must remain exported unless
//! a separate receiving absorber actually justifies this thermal destination.

use super::*;

#[derive(Clone, Copy)]
pub(super) struct LoadReceiver {
    pub(super) cell: CellId,
    geometry: CellGeometry,
    material: MaterialId,
    mass_kg: f64,
    body_id: u64,
    segment_id: u64,
    storage_basis: ThermalStorageBasis,
    caloric_identity: Option<[u8; 32]>,
}

impl LoadReceiver {
    pub(super) fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), PorousHeatingError> {
        checkpoint()?;
        let cell = thermal
            .cells()
            .get(&self.cell)
            .ok_or(PorousHeatingError::Invalid(
                "porous load heat receiver is no longer active",
            ))?;
        if cell.geometry != self.geometry
            || cell.material_id != self.material
            || cell.mass_kg != self.mass_kg
            || cell.body_id != self.body_id
            || cell.segment_id != self.segment_id
            || thermal.material_storage_basis(cell.material_id)? != self.storage_basis
            || thermal.material_caloric_identity(cell.material_id)? != self.caloric_identity
        {
            return Err(PorousHeatingError::Invalid(
                "porous load receiver changed geometry, mass or original thermal material",
            ));
        }
        let temperature = thermal.temperature_k(self.cell)?;
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(DepositionError::MaterialDomain.into());
        }
        Ok(())
    }
}

impl PreparedHeat {
    pub(super) fn bind_load_receiver(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        cell: CellId,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), PorousHeatingError> {
        checkpoint()?;
        self.validate(thermal, checkpoint)?;
        if self.load_receiver.is_some() || self.sources_w.contains_key(&cell) {
            return Err(PorousHeatingError::Invalid(
                "porous load heat needs one distinct receiver outside the porous cells",
            ));
        }
        let native = thermal
            .cells()
            .get(&cell)
            .ok_or(PorousHeatingError::Invalid(
                "porous load heat receiver is not active",
            ))?;
        let receiver = LoadReceiver {
            cell,
            geometry: native.geometry,
            material: native.material_id,
            mass_kg: native.mass_kg,
            body_id: native.body_id,
            segment_id: native.segment_id,
            storage_basis: thermal.material_storage_basis(native.material_id)?,
            caloric_identity: thermal.material_caloric_identity(native.material_id)?,
        };
        receiver.validate(thermal, checkpoint)?;
        self.load_receiver = Some(receiver);
        Ok(())
    }

    // One native source map is used for all modes/frequencies together. The
    // distinct receiver prevents its work overwriting or duplicating pore heat.
    pub(super) fn sources(&self) -> BTreeMap<CellId, f64> {
        let mut sources = self.sources_w.clone();
        if let Some(receiver) = self.load_receiver {
            sources.insert(receiver.cell, self.load_absorbed_power_w);
        }
        sources
    }

    pub(super) fn load_deposited_power_w(&self) -> f64 {
        if self.load_receiver.is_some() {
            self.load_absorbed_power_w
        } else {
            0.0
        }
    }

    pub(super) fn load_exported_power_w(&self) -> f64 {
        if self.load_receiver.is_some() {
            0.0
        } else {
            self.load_absorbed_power_w
        }
    }
}
