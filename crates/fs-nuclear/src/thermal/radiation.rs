//! Bind an externally resolved radiation deposition response to original cells.
//!
//! The response describes absorbed ENERGY fractions, not optical extinction or
//! photon-interaction probabilities. It is fixed over subsequent intervals and
//! applies to the emitter's complete integrated outward-radiation channel. Its
//! spectrum, geometry and material applicability remain the caller's data claim.
//! No thermal emissivity is substituted for ionizing-radiation absorption.

use super::*;
use crate::radiation::{RadiationDeposition, RadiationError, RadiationResponse};

/// Expected original material at a physical radiation receiver. A receiver can
/// also be a nuclear source; local deposition and radiation return add once.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RadiationReceiver {
    /// Existing active thermal cell; its ID is the response's `receiver_id`.
    pub cell_id: CellId,
    /// Material for which this externally supplied response was authored.
    pub material_id: MaterialId,
}

/// Immutable correspondence is shared across trial clones and accepted steps.
#[derive(Clone, Debug)]
pub(super) struct BoundRadiation {
    response: Arc<RadiationResponse>,
    receivers: BTreeMap<CellId, ThermalBinding>,
}

impl<S: NuclearHeatSource> NuclearHeating<S> {
    /// Deposit emitted radiation through an explicit, fixed transport response.
    ///
    /// Response source IDs must cover this complete source collection exactly;
    /// receiver bindings must cover every response destination exactly. Both use
    /// the original `CellId.0`. Receiving cells need no nuclear inventory. The
    /// response retains its explicit escape fraction and never changes the
    /// source's emission history or its nonrecoverable energy channel.
    ///
    /// The response must describe absorption of the full emitted energy mix over
    /// its intended temperature/composition/geometry domain. Changing spectra,
    /// scattering, photon flight times and response recalculation are not solved
    /// here. Binding snapshots the original material, calorics, mass, geometry
    /// and storage basis; subsequent changes refuse instead of reusing stale
    /// cell correspondence. A caller may replace the response for the same
    /// still-valid cell bindings at an accepted clock boundary, without replaying
    /// past radiation.
    pub fn with_radiation_response(
        mut self,
        thermal: &DepositionSimulation<'_>,
        response: Arc<RadiationResponse>,
        receivers: Vec<RadiationReceiver>,
    ) -> Result<Self, FissionThermalError> {
        self.validate(thermal, self.time_s, &mut || Ok(()))?;
        if self.pending.is_some() {
            return Err(FissionThermalError::Invalid(
                "cannot change radiation response during a pending interval",
            ));
        }
        if !response
            .source_ids()
            .eq(self.sources.keys().map(|cell| cell.0))
        {
            return Err(FissionThermalError::Invalid(
                "radiation response must cover the original source cohort exactly",
            ));
        }
        let bindings = bind_receivers(thermal, receivers, response.receiver_ids(), &mut || Ok(()))?;
        self.spectral_radiation = None;
        self.radiation = Some(Arc::new(BoundRadiation {
            response,
            receivers: bindings,
        }));
        Ok(self)
    }

    /// The immutable response used for future intervals, if one was supplied.
    #[must_use]
    pub fn radiation_response(&self) -> Option<&Arc<RadiationResponse>> {
        self.radiation.as_ref().map(|binding| &binding.response)
    }
}

impl BoundRadiation {
    pub(super) fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), FissionThermalError> {
        for (&id, binding) in &self.receivers {
            checkpoint()?;
            binding.validate(thermal, id)?;
            thermal.temperature_k(id)?;
        }
        Ok(())
    }

    pub(super) fn deposit(
        &self,
        emitted: &[(u64, f64)],
        tolerance_j: f64,
        heat_j: &mut BTreeMap<CellId, Sum>,
        account: &mut Account,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), FissionThermalError> {
        let mut checkpoint_error = None;
        let deposition = self.response.apply(emitted, tolerance_j, || {
            checkpoint().map_err(|error| {
                checkpoint_error = Some(error);
                RadiationError::Cancelled
            })
        });
        // Preserve the caller's exact native cancellation/error value; all
        // source histories and receiver heat are still unpublished candidates.
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        credit_deposition(deposition?, heat_j, account, checkpoint)
    }
}

pub(super) fn bind_receivers(
    thermal: &DepositionSimulation<'_>,
    receivers: Vec<RadiationReceiver>,
    expected: impl Iterator<Item = u64>,
    checkpoint: &mut dyn FnMut() -> Result<(), FissionThermalError>,
) -> Result<BTreeMap<CellId, ThermalBinding>, FissionThermalError> {
    if receivers.len() > thermal.cells().len() {
        return Err(FissionThermalError::Invalid(
            "radiation receivers need distinct active cells",
        ));
    }
    let mut bindings = BTreeMap::new();
    for receiver in receivers {
        checkpoint()?;
        let binding = ThermalBinding::new(thermal, receiver.cell_id, receiver.material_id)?;
        thermal.temperature_k(receiver.cell_id)?;
        if bindings.insert(receiver.cell_id, binding).is_some() {
            return Err(FissionThermalError::Invalid("duplicate radiation receiver"));
        }
    }
    if !expected.eq(bindings.keys().map(|cell| cell.0)) {
        return Err(FissionThermalError::Invalid(
            "radiation receivers must cover every response destination exactly",
        ));
    }
    Ok(bindings)
}

pub(super) fn credit_deposition(
    deposition: RadiationDeposition,
    heat_j: &mut BTreeMap<CellId, Sum>,
    account: &mut Account,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), FissionThermalError> {
    for (id, energy_j) in deposition.receiver_energy_j {
        checkpoint()?;
        heat_j.entry(CellId(id)).or_default().add(energy_j)?;
    }
    account.radiation_deposited.add(deposition.deposited_j)?;
    account.radiation_escaped.add(deposition.escaped_j)?;
    account.radiation_roundoff.add(deposition.roundoff_j)?;
    Ok(())
}
