//! Original finite-cell heat transport shared by angular adhesion loss owners.
//!
//! Each owner regenerates its physical loss, pins its actual receivers and
//! admits the positive interval before constructing this observer. The observer
//! distributes that loss at interval-average power, preserving the original
//! native time integration and independent per-cell transfer checks.

use super::JkrSnapHeatError;
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
};
use fs_tribo::DissipationStep;
use std::collections::BTreeMap;

/// Immutable receiver admission plus candidate-local integrated source heat.
/// The admission callback has no mutable side effects and retains the owner's
/// original source/material/cell identity checks on every native interval.
pub(super) struct SurfacePower<'a, A> {
    admit: &'a A,
    cells: [CellId; 2],
    start: f64,
    end: f64,
    powers: [f64; 2],
    injected: [Dd; 2],
}

// Only mutable physical accounting is copied; the original source/binding
// stays immutably borrowed. No Clone implementation is required of a callback.
impl<A> Clone for SurfacePower<'_, A> {
    fn clone(&self) -> Self {
        Self {
            admit: self.admit,
            cells: self.cells,
            start: self.start,
            end: self.end,
            powers: self.powers,
            injected: self.injected,
        }
    }
}

impl<'a, A> SurfacePower<'a, A> {
    /// Reconstruct average power over the owner's admitted mechanical interval.
    /// Actual native durations, including clock rounding, own the integral;
    /// `finish` refuses a mismatch instead of applying an endpoint correction.
    pub(super) fn new(
        cells: [CellId; 2],
        start: f64,
        end: f64,
        dt_s: f64,
        loss: DissipationStep,
        admit: &'a A,
    ) -> Result<Self, JkrSnapHeatError> {
        Ok(Self {
            admit,
            cells,
            start,
            end,
            powers: [
                finite(loss.surface_a_heat_j() / dt_s)?,
                finite(loss.surface_b_heat_j() / dt_s)?,
            ],
            injected: [Dd::ZERO; 2],
        })
    }

    pub(super) fn finish(
        self,
        expected: DissipationStep,
        tolerance: f64,
    ) -> Result<([Dd; 2], f64, Dd), JkrSnapHeatError> {
        let a = self.injected[0] - Dd::from_f64(expected.surface_a_heat_j());
        let b = self.injected[1] - Dd::from_f64(expected.surface_b_heat_j());
        let partition = Dd::from_f64(expected.surface_a_heat_j())
            + Dd::from_f64(expected.surface_b_heat_j())
            + Dd::from_f64(expected.other_work_j())
            - Dd::from_f64(expected.total_work_j());
        let residual = finite((a + b + partition).to_f64())?;
        let absolute = a.abs() + b.abs() + partition.abs();
        if residual.abs() > tolerance || finite(absolute.to_f64())? > tolerance {
            return Err(JkrSnapHeatError::TransferBalance(
                residual.abs().max(absolute.to_f64()),
            ));
        }
        Ok((self.injected, residual, absolute))
    }
}

impl<A> DepositionObserver for SurfacePower<'_, A>
where
    A: Fn(&DepositionSimulation<'_>) -> Result<(), JkrSnapHeatError>,
{
    type Error = JkrSnapHeatError;

    fn requires_step_interval(&self) -> bool {
        self.powers.iter().any(|p| *p != 0.0)
    }

    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        (self.admit)(thermal)?;
        if start_s < self.start || end_s > self.end || end_s <= start_s || dt_s <= 0.0 {
            return Err(JkrSnapHeatError::Invalid(
                "contact source interval differs from native thermal interval",
            ));
        }
        let mut sources = BTreeMap::new();
        for (side, power) in self.powers.into_iter().enumerate() {
            if power != 0.0 {
                self.injected[side] =
                    self.injected[side] + Dd::from_f64(power) * Dd::from_f64(dt_s);
                finite(self.injected[side].to_f64())?;
                sources.insert(self.cells[side], power);
            }
        }
        Ok(sources)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        (self.admit)(thermal)?;
        Ok(())
    }

    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        (self.admit)(thermal)?;
        Ok(())
    }
}

/// Preserve original compensated relative storage, heat and incoming-energy
/// audits without recomputing residuals from rounded aggregate energy totals.
pub(super) fn sum_audits(audits: [EnergyAudit; 3]) -> Result<EnergyAudit, JkrSnapHeatError> {
    let mut values = [Dd::ZERO; 5];
    let mut substeps = 0usize;
    let mut events = 0usize;
    for a in audits {
        for (sum, value) in values.iter_mut().zip([
            a.incoming_mass_kg,
            a.incoming_energy_j,
            a.external_heat_j,
            a.stored_energy_change_j,
            a.residual_j,
        ]) {
            *sum = *sum + Dd::from_f64(finite(value)?);
        }
        substeps = substeps
            .checked_add(a.substeps)
            .ok_or(DepositionError::Arithmetic)?;
        events = events
            .checked_add(a.events)
            .ok_or(DepositionError::Arithmetic)?;
    }
    Ok(EnergyAudit {
        incoming_mass_kg: finite(values[0].to_f64())?,
        incoming_energy_j: finite(values[1].to_f64())?,
        external_heat_j: finite(values[2].to_f64())?,
        stored_energy_change_j: finite(values[3].to_f64())?,
        residual_j: finite(values[4].to_f64())?,
        substeps,
        events,
    })
}

fn finite(value: f64) -> Result<f64, JkrSnapHeatError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
