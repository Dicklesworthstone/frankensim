//! Adaptive evolution of two-temperature coatings and one finite gas inventory.
//!
//! Each trial compares one original chamber interval with two half intervals.
//! Only the conservative fine history is accepted. The support and deposit keep
//! separate caloric coordinates; deposited mass is never scaled by support mass.
//! Their enthalpy increments, both temperatures, every gas species, original-law
//! gas sensible energy and integrated interface heat participate in the norm.
//!
//! This is local first-order temporal error control, not a global error bound or
//! a morphology/nucleation model. Both condensed inventories must start positive.
//! Only original inventory-limit refusals request a smaller step. Constitutive
//! domain failures, cancellation and exhausted work publish no partial history.

use super::super::super::{
    Total,
    advance::{GasChamberAdvanceAudit, GasSurfaceHistory},
};
use super::{
    ClosedGasParcel, GasResistiveChamberCandidate, GasResistiveSurface, GasResistiveSurfaceStep,
    GasSupportedChamberSolve, GasSurfaceError, SublimationError, ThermalAccommodation,
    ThermochemicalGasModel,
};
use crate::{MaterialError, thermochemical_gas::energy::GasEnergyError};
use fs_thermal::enthalpy::DepositionError;
use std::collections::BTreeMap;

/// Physical local-error scales and one trial allowance for the complete window.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasResistiveAdvanceSolve {
    pub step: GasSupportedChamberSolve,
    /// Applied independently to the deposit and EACH gas species [kg].
    pub absolute_mass_tolerance_kg: f64,
    /// Applied independently to support, deposit and gas temperatures [K].
    pub absolute_temperature_tolerance_k: f64,
    /// Applied separately to the original support/deposit caloric charts [J/kg].
    /// Relative scales use changes from the current state, not datum magnitudes.
    pub absolute_specific_enthalpy_tolerance_j_kg: f64,
    /// Applied to gas sensible energy and accepted interface-heat integrals [J].
    pub absolute_heat_tolerance_j: f64,
    pub relative_tolerance: f64,
    pub minimum_dt_s: f64,
    pub maximum_dt_s: f64,
    /// Includes rejected/refused trials. Each trial invokes at most THREE
    /// original kernels, each bounded by the original surface/storage controls.
    /// A zero-duration call invokes one original admission kernel.
    pub maximum_trials: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum GasResistiveAdvanceError {
    Invalid(&'static str),
    Surface(GasSurfaceError),
    /// Original source/geometry thermal failure in a composing interval owner.
    Thermal(DepositionError),
    TrialBudget {
        maximum_trials: usize,
    },
    MinimumStep {
        requested_dt_s: f64,
        minimum_dt_s: f64,
    },
}
impl From<GasSurfaceError> for GasResistiveAdvanceError {
    fn from(error: GasSurfaceError) -> Self {
        Self::Surface(error)
    }
}
impl From<SublimationError> for GasResistiveAdvanceError {
    fn from(error: SublimationError) -> Self {
        Self::Surface(error.into())
    }
}
impl From<GasEnergyError> for GasResistiveAdvanceError {
    fn from(error: GasEnergyError) -> Self {
        Self::Surface(error.into())
    }
}
impl From<DepositionError> for GasResistiveAdvanceError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl core::fmt::Display for GasResistiveAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive resistive deposition: {self:?}")
    }
}
impl std::error::Error for GasResistiveAdvanceError {}

/// Separate condensed storage along ONLY the accepted two-half-step histories.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasResistiveThermalHistory {
    pub id: u64,
    pub substrate_temperature_k: [f64; 2],
    /// Positive support -> deposit; internal to the combined chamber balance.
    pub interface_heat_j: f64,
    pub substrate_energy_change_j: f64,
    /// Includes the existing deposit-to-gas enthalpy datum conversion.
    pub deposit_energy_change_j: f64,
    pub substrate_energy_residual_j: f64,
    pub deposit_energy_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct GasResistiveAdvanceAudit {
    /// Existing complete molecular/species/energy and attempted-work accounts.
    /// Its surface temperatures are the deposit temperatures.
    pub chamber: GasChamberAdvanceAudit,
    /// The same canonical surface order, with independently retained support heat.
    pub thermal: Vec<GasResistiveThermalHistory>,
}

/// The entire proposed window is transactional; inputs are never mutated.
#[derive(Clone, Debug, PartialEq)]
pub struct GasResistiveAdvanceCandidate<'catalog> {
    endpoint: GasResistiveChamberCandidate<'catalog>,
    audit: GasResistiveAdvanceAudit,
}
impl<'catalog> GasResistiveAdvanceCandidate<'catalog> {
    /// Its step receipt describes only the last physical half step.
    #[must_use]
    pub const fn endpoint(&self) -> &GasResistiveChamberCandidate<'catalog> {
        &self.endpoint
    }
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.endpoint.gas()
    }
    #[must_use]
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.endpoint.catalog()
    }
    /// Accounts covering the complete accepted window, excluding rejected work.
    #[must_use]
    pub const fn audit(&self) -> &GasResistiveAdvanceAudit {
        &self.audit
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        GasResistiveChamberCandidate<'catalog>,
        GasResistiveAdvanceAudit,
    ) {
        (self.endpoint, self.audit)
    }
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasResistiveAdvanceError> {
    if cancelled() {
        Err(SublimationError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn finite(x: f64) -> Result<f64, GasResistiveAdvanceError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(GasResistiveAdvanceError::Invalid(
            "unrepresentable adaptive physical account",
        ))
    }
}
fn ratio(a: f64, b: f64, scale: f64) -> Result<f64, GasResistiveAdvanceError> {
    let value = (a - b).abs() / scale;
    if positive(scale) && value.is_finite() {
        Ok(value)
    } else {
        Err(GasResistiveAdvanceError::Invalid(
            "unrepresentable local physical error scale",
        ))
    }
}
fn state_ratio(
    a: f64,
    b: f64,
    start: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, GasResistiveAdvanceError> {
    ratio(
        a,
        b,
        absolute + relative * start.abs().max(a.abs()).max(b.abs()),
    )
}
fn change_ratio(
    a: f64,
    b: f64,
    start: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, GasResistiveAdvanceError> {
    ratio(
        a,
        b,
        absolute + relative * (a - start).abs().max((b - start).abs()),
    )
}

// Reuse the active gas composition's ORIGINAL integrated caloric law. Cancel
// species formation constants before arithmetic; never subtract absolute U.
fn sensible_energy(
    gas: &ClosedGasParcel,
    start_temperature_k: f64,
) -> Result<f64, GasResistiveAdvanceError> {
    let change = gas
        .model()
        .thermodynamics()
        .specific_internal_energy_change(
            fs_qty::Temperature::new(start_temperature_k),
            fs_qty::Temperature::new(gas.state().state().temperature),
        )
        .map_err(|error| {
            GasEnergyError::Material(MaterialError::Parameters {
                what: format!("adaptive resistive gas caloric increment: {error}"),
            })
        })?
        .value();
    finite(gas.mass_kg() * change)
}

fn compare(
    before: &[GasResistiveSurface<'_>],
    gas: &ClosedGasParcel,
    coarse: &GasResistiveChamberCandidate<'_>,
    first: &GasResistiveChamberCandidate<'_>,
    fine: &GasResistiveChamberCandidate<'_>,
    controls: GasResistiveAdvanceSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, GasResistiveAdvanceError> {
    let t0 = gas.state().state().temperature;
    let mut maximum = state_ratio(
        coarse.gas().state().state().temperature,
        fine.gas().state().state().temperature,
        t0,
        controls.absolute_temperature_tolerance_k,
        controls.relative_tolerance,
    )?;
    maximum = maximum.max(change_ratio(
        sensible_energy(coarse.gas(), t0)?,
        sensible_energy(fine.gas(), t0)?,
        0.0,
        controls.absolute_heat_tolerance_j,
        controls.relative_tolerance,
    )?);
    for (((original, a), half), b) in before
        .iter()
        .zip(coarse.surfaces())
        .zip(first.surfaces())
        .zip(fine.surfaces())
    {
        poll(cancelled)?;
        maximum = maximum.max(state_ratio(
            a.state.deposit.mass_kg,
            b.state.deposit.mass_kg,
            original.state.deposit.mass_kg,
            controls.absolute_mass_tolerance_kg,
            controls.relative_tolerance,
        )?);
        for (ha, hb, h0) in [
            (
                a.state.substrate_specific_enthalpy_j_kg,
                b.state.substrate_specific_enthalpy_j_kg,
                original.state.substrate_specific_enthalpy_j_kg,
            ),
            (
                a.state.deposit.specific_enthalpy_j_kg,
                b.state.deposit.specific_enthalpy_j_kg,
                original.state.deposit.specific_enthalpy_j_kg,
            ),
        ] {
            maximum = maximum.max(change_ratio(
                ha,
                hb,
                h0,
                controls.absolute_specific_enthalpy_tolerance_j_kg,
                controls.relative_tolerance,
            )?);
        }
        for (ta, tb, t0) in [
            (
                a.substrate_temperature_k[1],
                b.substrate_temperature_k[1],
                a.substrate_temperature_k[0],
            ),
            (
                a.exchange.surface_temperature_k[1],
                b.exchange.surface_temperature_k[1],
                a.exchange.surface_temperature_k[0],
            ),
        ] {
            maximum = maximum.max(state_ratio(
                ta,
                tb,
                t0,
                controls.absolute_temperature_tolerance_k,
                controls.relative_tolerance,
            )?);
        }
        maximum = maximum.max(change_ratio(
            a.interface_heat_j,
            finite(half.interface_heat_j + b.interface_heat_j)?,
            0.0,
            controls.absolute_heat_tolerance_j,
            controls.relative_tolerance,
        )?);
    }
    for species in fine.catalog().thermodynamics().components() {
        poll(cancelled)?;
        let id = species.species();
        maximum = maximum.max(state_ratio(
            coarse.gas().species_mass_kg(id)?,
            fine.gas().species_mass_kg(id)?,
            gas.species_mass_kg(id)?,
            controls.absolute_mass_tolerance_kg,
            controls.relative_tolerance,
        )?);
    }
    Ok(maximum)
}

fn smaller(requested: f64, previous: f64, minimum: f64) -> Result<f64, GasResistiveAdvanceError> {
    let next = requested.min(0.5 * previous);
    if !positive(next) || next < minimum {
        Err(GasResistiveAdvanceError::MinimumStep {
            requested_dt_s: next,
            minimum_dt_s: minimum,
        })
    } else {
        Ok(next)
    }
}

#[derive(Default)]
struct SurfaceTotals {
    outgoing: Total,
    incoming: Total,
    transferred: Total,
    stream: Total,
    reflected: Total,
    stored: Total,
    residual: Total,
    recession: Total,
    support_stored: Total,
    deposit_stored: Total,
    interface: Total,
    support_residual: Total,
    deposit_residual: Total,
    initial_support_temperature: Option<f64>,
    initial_deposit_temperature: Option<f64>,
}
impl SurfaceTotals {
    fn add(&mut self, row: &GasResistiveSurfaceStep) -> Result<(), GasSurfaceError> {
        let exchange = &row.exchange;
        self.initial_support_temperature
            .get_or_insert(row.substrate_temperature_k[0]);
        self.initial_deposit_temperature
            .get_or_insert(exchange.surface_temperature_k[0]);
        self.outgoing.add(exchange.outgoing_mass_kg)?;
        self.incoming.add(exchange.incoming_mass_kg)?;
        self.transferred.add(exchange.transferred_mass_kg)?;
        self.stream.add(exchange.energy_into_gas_j)?;
        self.reflected.add(row.reflected_heat_into_gas_j)?;
        self.stored.add(exchange.solid_energy_change_j)?;
        self.residual.add(exchange.energy_residual_j)?;
        self.recession.add(exchange.recession_m)?;
        self.support_stored.add(row.substrate_energy_change_j)?;
        self.deposit_stored.add(row.deposit_energy_change_j)?;
        self.interface.add(row.interface_heat_j)?;
        self.support_residual.add(row.substrate_energy_residual_j)?;
        self.deposit_residual.add(row.deposit_energy_residual_j)
    }
}

/// Evolve separate coating/support temperatures against one finite gas.
///
/// Areas, contact conductances and coefficients remain authored and fixed here.
/// Original molecular and thermal laws are reevaluated on every physical half
/// step. Source/caloric refusal is not treated as a timestep error. A final
/// remainder below minimum_dt_s is attempted once; any needed reduction refuses.
#[allow(clippy::too_many_arguments)]
pub fn advance_resistive_chamber_adaptive<'catalog>(
    surfaces: &[GasResistiveSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasResistiveAdvanceSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasResistiveAdvanceCandidate<'catalog>, GasResistiveAdvanceError> {
    advance_with(
        surfaces,
        gas,
        catalog,
        duration_s,
        controls,
        |rows, gas, dt, cancel| {
            super::advance(
                rows,
                accommodation,
                None,
                gas,
                catalog,
                dt,
                controls.step,
                false,
                cancel,
            )
            .map_err(Into::into)
        },
        |_, _, _, _, _| Ok(0.0),
        cancelled,
    )
}

/// The same adaptive owner with independent per-surface gas-film resistance.
/// Missing IDs use the free-molecular limit; zero explicitly seals mass transfer.
#[allow(clippy::too_many_arguments)]
pub fn advance_resistive_chamber_adaptive_with_films<'catalog>(
    surfaces: &[GasResistiveSurface<'_>],
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    films: &BTreeMap<u64, f64>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasResistiveAdvanceSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasResistiveAdvanceCandidate<'catalog>, GasResistiveAdvanceError> {
    advance_with(
        surfaces,
        gas,
        catalog,
        duration_s,
        controls,
        |rows, gas, dt, cancel| {
            super::advance(
                rows,
                accommodation,
                Some(films),
                gas,
                catalog,
                dt,
                controls.step,
                false,
                cancel,
            )
            .map_err(Into::into)
        },
        |_, _, _, _, _| Ok(0.0),
        cancelled,
    )
}

// One controller carries BOTH condensed states. Composing physical geometry
// owners may rebind each interval and add an error coordinate; they do not own
// another retry loop, gas reconstruction, heat integral or acceptance decision.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(super) fn advance_with<'catalog, C: FnMut() -> bool>(
    surfaces: &[GasResistiveSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasResistiveAdvanceSolve,
    mut evaluate_interval: impl FnMut(
        &[GasResistiveSurface<'_>],
        &ClosedGasParcel,
        f64,
        &mut C,
    ) -> Result<
        GasResistiveChamberCandidate<'catalog>,
        GasResistiveAdvanceError,
    >,
    mut extra_error: impl FnMut(
        &[GasResistiveSurface<'_>],
        &GasResistiveChamberCandidate<'_>,
        &GasResistiveChamberCandidate<'_>,
        f64,
        &mut C,
    ) -> Result<f64, GasResistiveAdvanceError>,
    mut cancelled: C,
) -> Result<GasResistiveAdvanceCandidate<'catalog>, GasResistiveAdvanceError> {
    poll(&mut cancelled)?;
    if !duration_s.is_finite()
        || duration_s < 0.0
        || !positive(controls.absolute_mass_tolerance_kg)
        || !positive(controls.absolute_temperature_tolerance_k)
        || !positive(controls.absolute_specific_enthalpy_tolerance_j_kg)
        || !positive(controls.absolute_heat_tolerance_j)
        || !controls.relative_tolerance.is_finite()
        || controls.relative_tolerance < 0.0
        || !positive(controls.minimum_dt_s)
        || !positive(controls.maximum_dt_s)
        || controls.minimum_dt_s > controls.maximum_dt_s
        || controls.maximum_trials == 0
        || surfaces.len() > controls.step.chamber.maximum_surfaces
    {
        return Err(GasResistiveAdvanceError::Invalid(
            "finite forward duration and positive bounded physical tolerances/work required",
        ));
    }
    let mut current = surfaces.to_vec();
    current.sort_by_key(|row| row.id);
    if current.windows(2).any(|pair| pair[0].id == pair[1].id) {
        return Err(GasResistiveAdvanceError::Invalid(
            "duplicate resistive surface identity",
        ));
    }
    for row in &current {
        poll(&mut cancelled)?;
        if !positive(row.state.substrate_mass_kg) || !positive(row.state.deposit.mass_kg) {
            return Err(GasResistiveAdvanceError::Invalid(
                "adaptive resistance needs existing positive support and deposit masses",
            ));
        }
    }
    let original = current.clone();
    let mut middle = current.clone();
    let mut sums: Vec<_> = (0..current.len())
        .map(|_| SurfaceTotals::default())
        .collect();
    let mut gas_energy = Total::default();
    let mut gas_residual = Total::default();
    let mut endpoint: Option<GasResistiveChamberCandidate<'catalog>> = None;
    let mut elapsed = 0.0;
    let mut elapsed_sum = Total::default();
    let mut dt = controls.maximum_dt_s.min(duration_s);
    let (mut trials, mut rejected, mut accepted, mut evaluations) =
        (0usize, 0usize, 0usize, 0usize);
    let mut maximum_error: f64 = 0.0;
    while elapsed < duration_s {
        poll(&mut cancelled)?;
        if trials == controls.maximum_trials {
            return Err(GasResistiveAdvanceError::TrialBudget {
                maximum_trials: controls.maximum_trials,
            });
        }
        trials += 1;
        let remaining = duration_s - elapsed;
        dt = dt.min(remaining);
        if !positive(dt) || elapsed + dt == elapsed || dt * 0.5 == 0.0 {
            return Err(GasResistiveAdvanceError::MinimumStep {
                requested_dt_s: dt,
                minimum_dt_s: controls.minimum_dt_s,
            });
        }
        let current_gas = endpoint.as_ref().map_or(gas, |value| value.gas());
        let mut evaluate = |rows: &[GasResistiveSurface<'_>], gas: &ClosedGasParcel, duration| {
            evaluations = evaluations
                .checked_add(1)
                .ok_or(GasResistiveAdvanceError::Invalid(
                    "kernel-work count overflow",
                ))?;
            let result = evaluate_interval(rows, gas, duration, &mut cancelled)?;
            if rows.len() != result.surfaces().len()
                || rows.iter().zip(result.surfaces()).any(|(old, next)| {
                    old.id != next.id || old.state.substrate_mass_kg != next.state.substrate_mass_kg
                })
            {
                return Err(GasResistiveAdvanceError::Invalid(
                    "interval changed the original support roster",
                ));
            }
            Ok(result)
        };
        let trial = (|| {
            let coarse = evaluate(&current, current_gas, dt)?;
            let first = evaluate(&current, current_gas, 0.5 * dt)?;
            for (row, next) in middle.iter_mut().zip(first.surfaces()) {
                row.state = next.state;
            }
            let second = evaluate(&middle, first.gas(), 0.5 * dt)?;
            Ok::<_, GasResistiveAdvanceError>((coarse, first, second))
        })();
        let (coarse, first, second) = match trial {
            Ok(value) => value,
            Err(GasResistiveAdvanceError::Surface(GasSurfaceError::InventoryLimit {
                maximum_duration_s,
            })) => {
                rejected += 1;
                dt = smaller(0.8 * maximum_duration_s, dt, controls.minimum_dt_s)?;
                continue;
            }
            Err(error) => return Err(error),
        };
        let mut error = compare(
            &current,
            current_gas,
            &coarse,
            &first,
            &second,
            controls,
            &mut cancelled,
        )?;
        let extra = extra_error(
            &current,
            &coarse,
            &second,
            controls.relative_tolerance,
            &mut cancelled,
        )?;
        if !extra.is_finite() || extra < 0.0 {
            return Err(GasResistiveAdvanceError::Invalid(
                "invalid additional physical error ratio",
            ));
        }
        error = error.max(extra);
        if error > 1.0 {
            rejected += 1;
            dt = smaller(
                0.8 * dt / fs_math::det::sqrt(error),
                dt,
                controls.minimum_dt_s,
            )?;
            continue;
        }
        for ((sum, a), b) in sums.iter_mut().zip(first.surfaces()).zip(second.surfaces()) {
            poll(&mut cancelled)?;
            sum.add(a)?;
            sum.add(b)?;
        }
        for part in [&first, &second] {
            gas_energy.add(part.step().gas_energy_change_j)?;
            gas_residual.add(part.gas_exchange().step().energy_residual_j)?;
        }
        for (row, next) in current.iter_mut().zip(second.surfaces()) {
            row.state = next.state;
        }
        endpoint = Some(second);
        accepted += 1;
        maximum_error = maximum_error.max(error);
        elapsed_sum.add(dt)?;
        elapsed = if dt == remaining {
            duration_s
        } else {
            elapsed_sum.value()?
        };
        let factor = if error == 0.0 {
            2.0
        } else {
            (0.8 / fs_math::det::sqrt(error)).clamp(0.5, 2.0)
        };
        dt = (dt * factor)
            .min(controls.maximum_dt_s)
            .max(controls.minimum_dt_s);
    }
    let endpoint = match endpoint {
        Some(value) => value,
        None => {
            evaluations += 1;
            evaluate_interval(&current, gas, 0.0, &mut cancelled)?
        }
    };
    let mut solid_energy = Total::default();
    let mut solid_residual = Total::default();
    let mut solid_mass = Total::default();
    let mut stream_energy = Total::default();
    let mut reflected_energy = Total::default();
    let mut species_accounts = BTreeMap::<_, (Total, Total)>::new();
    let mut histories = Vec::with_capacity(original.len());
    let mut thermal = Vec::with_capacity(original.len());
    let mass_tolerance = controls.step.chamber.surface.gas.mass_tolerance_kg;
    let energy_tolerance = controls.step.chamber.surface.gas.energy_tolerance_j;
    for ((initial, end), sum) in original.iter().zip(endpoint.surfaces()).zip(sums) {
        poll(&mut cancelled)?;
        let mass_change = finite(end.state.deposit.mass_kg - initial.state.deposit.mass_kg)?;
        let transferred = sum.transferred.value()?;
        if (mass_change + transferred).abs() > mass_tolerance {
            return Err(GasResistiveAdvanceError::Invalid(
                "whole-history coating mass balance failed",
            ));
        }
        for residual in [
            sum.residual.value()?,
            sum.support_residual.value()?,
            sum.deposit_residual.value()?,
        ] {
            if residual.abs() > energy_tolerance {
                return Err(GasSurfaceError::EnergyBalance {
                    residual_j: residual,
                    tolerance_j: energy_tolerance,
                }
                .into());
            }
        }
        solid_energy.add(sum.stored.value()?)?;
        solid_residual.add(sum.residual.value()?)?;
        solid_mass.add(mass_change)?;
        stream_energy.add(sum.stream.value()?)?;
        reflected_energy.add(sum.reflected.value()?)?;
        let species = species_accounts
            .entry(initial.binding.species.clone())
            .or_default();
        species.0.add(transferred)?;
        species.1.add(mass_change)?;
        histories.push(GasSurfaceHistory {
            id: initial.id,
            outgoing_mass_kg: sum.outgoing.value()?,
            incoming_mass_kg: sum.incoming.value()?,
            transferred_mass_kg: transferred,
            energy_into_gas_j: sum.stream.value()?,
            reflected_heat_into_gas_j: sum.reflected.value()?,
            recession_m: sum.recession.value()?,
            surface_temperature_k: [
                sum.initial_deposit_temperature
                    .unwrap_or(end.exchange.surface_temperature_k[0]),
                end.exchange.surface_temperature_k[1],
            ],
        });
        thermal.push(GasResistiveThermalHistory {
            id: initial.id,
            substrate_temperature_k: [
                sum.initial_support_temperature
                    .unwrap_or(end.substrate_temperature_k[0]),
                end.substrate_temperature_k[1],
            ],
            interface_heat_j: sum.interface.value()?,
            substrate_energy_change_j: sum.support_stored.value()?,
            deposit_energy_change_j: sum.deposit_stored.value()?,
            substrate_energy_residual_j: sum.support_residual.value()?,
            deposit_energy_residual_j: sum.deposit_residual.value()?,
        });
    }
    let mut maximum_species_mass_residual: f64 = 0.0;
    for species in catalog.thermodynamics().components() {
        poll(&mut cancelled)?;
        let id = species.species();
        let change = finite(endpoint.gas().species_mass_kg(id)? - gas.species_mass_kg(id)?)?;
        let (transported, condensed) = match species_accounts.get(id) {
            Some((stream, stored)) => (stream.value()?, stored.value()?),
            None => (0.0, 0.0),
        };
        maximum_species_mass_residual = maximum_species_mass_residual
            .max(finite(change - transported)?.abs())
            .max(finite(change + condensed)?.abs());
    }
    let gas_stream_residual = gas_residual.value()?;
    let solid_stream_residual = solid_residual.value()?;
    let energy_residual = finite(gas_stream_residual + solid_stream_residual)?;
    let worst_energy = energy_residual
        .abs()
        .max(gas_stream_residual.abs())
        .max(solid_stream_residual.abs());
    if worst_energy > energy_tolerance {
        return Err(GasSurfaceError::EnergyBalance {
            residual_j: worst_energy,
            tolerance_j: energy_tolerance,
        }
        .into());
    }
    let mass_residual = finite(solid_mass.value()? + (endpoint.gas().mass_kg() - gas.mass_kg()))?;
    if mass_residual.abs() > mass_tolerance || maximum_species_mass_residual > mass_tolerance {
        return Err(GasResistiveAdvanceError::Invalid(
            "whole-history solid/gas species balance failed",
        ));
    }
    poll(&mut cancelled)?;
    Ok(GasResistiveAdvanceCandidate {
        endpoint,
        audit: GasResistiveAdvanceAudit {
            chamber: GasChamberAdvanceAudit {
                elapsed_s: elapsed,
                trials,
                rejected_trials: rejected,
                accepted_intervals: accepted,
                kernel_evaluations: evaluations,
                maximum_accepted_error_ratio: maximum_error,
                surfaces: histories,
                energy_into_gas_j: stream_energy.value()?,
                reflected_heat_into_gas_j: reflected_energy.value()?,
                solid_energy_change_j: solid_energy.value()?,
                gas_energy_change_j: gas_energy.value()?,
                energy_residual_j: energy_residual,
                gas_stream_energy_residual_j: gas_stream_residual,
                solid_stream_energy_residual_j: solid_stream_residual,
                mass_residual_kg: mass_residual,
                maximum_species_mass_residual_kg: maximum_species_mass_residual,
            },
            thermal,
        },
    })
}

#[cfg(test)]
mod tests;
