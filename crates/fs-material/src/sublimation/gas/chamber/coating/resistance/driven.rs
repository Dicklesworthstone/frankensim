//! External heat and nonlinear contact in the original resistive chamber.
//!
//! Prescribed support/deposit powers are integrated over each ACTUAL physical
//! trial interval. The shared adaptive owner counts only accepted half-step
//! heat, keeps both storage states and resolves molecular feedback on subsequent
//! intervals. This is first-order coupled evolution, not an imposed-temperature
//! thermostat, feedback controller or unbounded heat reservoir.

use super::advance::{
    GasResistiveAdvanceCandidate, GasResistiveAdvanceError, GasResistiveAdvanceSolve,
};
use super::{
    ClosedGasParcel, GasResistiveDrive, GasResistiveSurface, GasSurfaceError,
    ResistiveInterfaceLaw, ThermalAccommodation, ThermochemicalGasModel, invalid, poll,
};
use std::collections::{BTreeMap, BTreeSet};

/// Prescribed external powers on the two original condensed bodies [W].
/// Positive values heat the named body; negative values cool it. An optional
/// passive interface replaces the surface's scalar conductance, independently
/// of these external sources and any gas-film mass resistance.
#[derive(Clone, Copy, Default)]
pub struct GasResistivePower<'a> {
    pub substrate_power_w: f64,
    pub deposit_power_w: f64,
    pub interface: Option<&'a dyn ResistiveInterfaceLaw>,
}
impl<'a> GasResistivePower<'a> {
    /// Integrate each prescribed power over one actual interval. A zero-time
    /// admission supplies zero heat; nonfinite or underflowed active work refuses.
    pub fn for_interval(self, duration_s: f64) -> Result<GasResistiveDrive<'a>, GasSurfaceError> {
        if !duration_s.is_finite()
            || duration_s < 0.0
            || !self.substrate_power_w.is_finite()
            || !self.deposit_power_w.is_finite()
        {
            return Err(invalid(
                "finite prescribed powers and nonnegative duration required",
            ));
        }
        let substrate_heat_j = self.substrate_power_w * duration_s;
        let deposit_heat_j = self.deposit_power_w * duration_s;
        if !substrate_heat_j.is_finite()
            || !deposit_heat_j.is_finite()
            || !(substrate_heat_j + deposit_heat_j).is_finite()
            || (duration_s > 0.0 && self.substrate_power_w != 0.0 && substrate_heat_j == 0.0)
            || (duration_s > 0.0 && self.deposit_power_w != 0.0 && deposit_heat_j == 0.0)
        {
            return Err(invalid(
                "prescribed resistive heating/cooling is unrepresentable",
            ));
        }
        Ok(GasResistiveDrive {
            substrate_heat_j,
            deposit_heat_j,
            interface: self.interface,
        })
    }
}

/// Check the accepted source history against the independently integrated
/// prescribed powers. Every body is checked separately, including net-zero pairs.
pub(super) fn admit_power_history(
    powers: &BTreeMap<u64, GasResistivePower<'_>>,
    history: &GasResistiveAdvanceCandidate<'_>,
    duration_s: f64,
    tolerance_j: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), GasResistiveAdvanceError> {
    for row in &history.audit().thermal {
        poll(cancelled)?;
        let expected = powers
            .get(&row.id)
            .copied()
            .unwrap_or_default()
            .for_interval(duration_s)?;
        for residual_j in [
            row.supplied_substrate_heat_j - expected.substrate_heat_j,
            row.supplied_deposit_heat_j - expected.deposit_heat_j,
        ] {
            if !residual_j.is_finite() || residual_j.abs() > tolerance_j {
                return Err(GasSurfaceError::EnergyBalance {
                    residual_j,
                    tolerance_j,
                }
                .into());
            }
        }
    }
    Ok(())
}

/// Evolve prescribed support/coating heating or cooling with one finite gas.
///
/// The sparse power map must reference existing surface IDs. Other surfaces
/// retain zero external heat and their original contact. Source models are
/// borrowed; one reusable drive map is prepared before the trial loop. Full
/// and half intervals receive `P * dt`, so rejected trials cannot contribute
/// source energy or change its prescribed integral. Both local source heats
/// participate in the shared physical error norm and accepted energy accounts.
/// Initial positive films, source domains, work budgets and rollback semantics
/// retain the original adaptive owner's meaning.
#[allow(clippy::too_many_arguments)]
pub fn advance_resistive_chamber_adaptive_driven<'catalog>(
    surfaces: &[GasResistiveSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    films: Option<&BTreeMap<u64, f64>>,
    powers: &BTreeMap<u64, GasResistivePower<'_>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasResistiveAdvanceSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasResistiveAdvanceCandidate<'catalog>, GasResistiveAdvanceError> {
    let mut films = films;
    advance_with_films(
        surfaces,
        accommodation,
        &mut films,
        powers,
        gas,
        catalog,
        duration_s,
        controls,
        cancelled,
    )
}

/// The same power/history owner with interval-resolved transport inputs.
#[allow(clippy::too_many_arguments)]
pub(super) fn advance_with_films<'catalog>(
    surfaces: &[GasResistiveSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    films: &mut dyn super::film::TrialFilms,
    powers: &BTreeMap<u64, GasResistivePower<'_>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasResistiveAdvanceSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasResistiveAdvanceCandidate<'catalog>, GasResistiveAdvanceError> {
    poll(&mut cancelled)?;
    if surfaces.len() > controls.step.chamber.maximum_surfaces || powers.len() > surfaces.len() {
        return Err(GasResistiveAdvanceError::Invalid(
            "resistive power roster exceeds the surface work bound",
        ));
    }
    let ids: BTreeSet<_> = surfaces.iter().map(|row| row.id).collect();
    let mut drives = BTreeMap::new();
    for (&id, &power) in powers {
        poll(&mut cancelled)?;
        if !ids.contains(&id) {
            return Err(GasResistiveAdvanceError::Invalid(
                "power references an unknown resistive surface",
            ));
        }
        drives.insert(id, power.for_interval(0.0)?);
    }
    let history = super::advance::advance_with(
        surfaces,
        gas,
        catalog,
        duration_s,
        controls,
        |rows, gas, dt, cancelled| {
            for ((_, power), (_, drive)) in powers.iter().zip(drives.iter_mut()) {
                poll(cancelled)?;
                *drive = power.for_interval(dt)?;
            }
            let films = films.resolve(rows, gas, cancelled)?;
            super::advance(
                rows,
                accommodation,
                films,
                Some(&drives),
                gas,
                catalog,
                dt,
                controls.step,
                false,
                cancelled,
            )
            .map_err(Into::into)
        },
        |_, _, _, _, _| Ok(0.0),
        &mut cancelled,
    )?;
    admit_power_history(
        powers,
        &history,
        duration_s,
        controls.step.chamber.surface.gas.energy_tolerance_j,
        &mut cancelled,
    )?;
    poll(&mut cancelled)?;
    Ok(history)
}

#[cfg(test)]
mod tests;
