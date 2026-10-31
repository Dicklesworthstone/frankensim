//! Complete finite-mixture breakthrough windows using original implicit steps.
//!
//! A full step estimates the error of two half steps. Only the latter advance
//! inventories and integrated streams. Every film cell and every actual gas
//! species is compared independently; a small net transfer cannot hide opposing
//! permeation. This first-order local estimate is not a global error certificate.
//! All original source, linear, gas, pressure and caloric refusals stay hard.
//! Isothermal gas energy is audited in the original source-relative coordinate;
//! no unknown sorbed caloric store or heat of sorption is manufactured.

use super::{
    Error as NativeError, IsothermalGasMembraneMixture, IsothermalGasMembraneMixtureStep,
    IsothermalGasMembraneMixtureWorkspace, R, SpeciesId, Sum, inventory,
};
use fs_qty::Temperature;

/// Temporal accuracy and fixed whole-window physical budgets, all explicit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MembraneMixtureAdvanceControls {
    pub initial_step_s: f64,
    /// Smallest comparison interval; accepted physical halves may be smaller.
    /// A shorter final remainder is attempted without relaxing its accuracy.
    pub minimum_step_s: f64,
    /// Nominal comparison ceiling. A terminal clock one ULP beyond the
    /// proposed endpoint is included before solving, so the actual final span
    /// may exceed this ceiling by that clock roundoff only.
    pub maximum_step_s: f64,
    /// Original complete-mixture calls, including rejected probes; 3..=100000.
    pub maximum_evaluations: usize,
    /// All retained film cells over all species, admitted before cloning states.
    pub maximum_total_cells: usize,
    /// Per-cell Henry pressure and per-gas species partial-pressure error [Pa].
    pub absolute_pressure_tolerance_pa: f64,
    /// Each actual gas species inventory receives its own amount comparison.
    pub absolute_amount_tolerance_mol: f64,
    /// Multiplies each coordinate's change from its last accepted value; [0,1).
    pub relative_tolerance: f64,
    /// One fixed amount budget per species, per gas boundary and film account.
    pub species_balance_tolerance_mol: f64,
    pub mass_balance_tolerance_kg: f64,
    /// Each original gas's source-relative complete-window thermal residual.
    pub gas_energy_tolerance_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MembraneMixtureAdvanceError {
    Invalid(&'static str),
    Native(NativeError),
    Accuracy {
        time_s: f64,
        step_s: f64,
        error_ratio: f64,
    },
    WorkLimit {
        time_s: f64,
        evaluations: usize,
    },
    TimeResolution {
        time_s: f64,
        target_s: f64,
    },
    SpeciesBalance {
        species: SpeciesId,
        account: MembraneMixtureBalanceAccount,
        residual_mol: f64,
        tolerance_mol: f64,
    },
    MassBalance {
        residual_kg: f64,
        tolerance_kg: f64,
    },
    GasEnergyBalance {
        reservoir: usize,
        residual_j: f64,
        tolerance_j: f64,
    },
    Cancelled,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MembraneMixtureBalanceAccount {
    LeftGas,
    RightGas,
    FilmBoundary,
    ClosedSpecies,
}
impl core::fmt::Display for MembraneMixtureAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive finite membrane mixture: {self:?}")
    }
}
impl std::error::Error for MembraneMixtureAdvanceError {}
impl From<NativeError> for MembraneMixtureAdvanceError {
    fn from(error: NativeError) -> Self {
        if error == NativeError::Cancelled {
            Self::Cancelled
        } else {
            Self::Native(error)
        }
    }
}
type Result<T> = core::result::Result<T, MembraneMixtureAdvanceError>;
use MembraneMixtureAdvanceError as Error;

/// Separate original full/half-step factor caches. Repeated accepted intervals
/// reuse both native tridiagonal factorizations; coefficients are never replaced.
/// Failed proposals may alter these numerical caches, never physical state.
#[derive(Debug)]
pub struct MembraneMixtureAdvanceWorkspace {
    full: IsothermalGasMembraneMixtureWorkspace,
    half: IsothermalGasMembraneMixtureWorkspace,
}

/// Accepted-only streams in canonical species order. Sides are [left,right];
/// positive signed boundary amounts/energies enter the gas, not the film.
#[derive(Clone, Debug, PartialEq)]
pub struct MembraneMixtureIntegratedSpecies {
    pub species: SpeciesId,
    pub amount_into_gas_mol: [f64; 2],
    pub gross_amount_into_gas_mol: [f64; 2],
    pub gross_amount_out_of_gas_mol: [f64; 2],
    /// Original partial-enthalpy streams; their formation reference is retained.
    pub enthalpy_into_gas_j: [f64; 2],
    pub gross_absolute_enthalpy_j: [f64; 2],
    pub thermostat_heat_into_gas_j: [f64; 2],
    pub gross_absolute_thermostat_heat_j: [f64; 2],
    pub actual_gas_amount_changes_mol: [f64; 2],
    /// Computed from original-to-final CELL changes, not summed step residuals.
    pub actual_sorbed_amount_change_mol: f64,
    pub gas_boundary_residuals_mol: [f64; 2],
    pub film_boundary_residual_mol: f64,
    pub closed_species_residual_mol: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MembraneMixtureAdvanceReport {
    pub interval_s: [f64; 2],
    pub accepted_intervals: usize,
    pub accepted_substeps: usize,
    pub rejected_intervals: usize,
    pub evaluations: usize,
    pub minimum_accepted_substep_s: f64,
    pub maximum_accepted_substep_s: f64,
    pub maximum_accepted_error_ratio: f64,
    pub species: Vec<MembraneMixtureIntegratedSpecies>,
    pub actual_gas_mass_changes_kg: [f64; 2],
    pub actual_sorbed_mass_change_kg: f64,
    pub total_mass_residual_kg: f64,
    pub direct_total_mass_residual_kg: f64,
    pub enthalpy_into_gas_j: [f64; 2],
    pub thermostat_heat_into_gas_j: [f64; 2],
    /// Original final-mixture sensible storage minus accepted relative input.
    /// Species-reference transport cancels analytically; inventory rounding is
    /// separately bounded. Absolute enthalpy subtraction is not this residual.
    pub source_relative_gas_energy_residual_j: [f64; 2],
    /// Original native receipts for the last accepted half step, not history.
    pub last_substep: IsothermalGasMembraneMixtureStep,
}
#[derive(Clone, Debug)]
pub struct IsothermalGasMembraneMixtureAdvanceCandidate<'a> {
    state: IsothermalGasMembraneMixture<'a>,
    report: MembraneMixtureAdvanceReport,
}
impl<'a> IsothermalGasMembraneMixtureAdvanceCandidate<'a> {
    pub const fn state(&self) -> &IsothermalGasMembraneMixture<'a> {
        &self.state
    }
    pub const fn report(&self) -> &MembraneMixtureAdvanceReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        IsothermalGasMembraneMixture<'a>,
        MembraneMixtureAdvanceReport,
    ) {
        (self.state, self.report)
    }
}

#[derive(Default)]
struct Ledger {
    amount: [Sum; 2],
    incoming: [Sum; 2],
    outgoing: [Sum; 2],
    enthalpy: [Sum; 2],
    absolute_enthalpy: [Sum; 2],
    heat: [Sum; 2],
    absolute_heat: [Sum; 2],
    relative: [Sum; 2],
}
impl Ledger {
    fn add_step(
        &mut self,
        row: &super::MembraneMixtureSpeciesStep,
        mass: f64,
        temperature: f64,
    ) -> Result<()> {
        let amounts = [
            -row.diffusion.membrane.left_reservoir_to_membrane_mol,
            row.diffusion.membrane.membrane_to_right_reservoir_mol,
        ];
        for side in 0..2 {
            let amount = amounts[side];
            add(&mut self.amount[side], amount)?;
            add(&mut self.incoming[side], amount.max(0.0))?;
            add(&mut self.outgoing[side], (-amount).max(0.0))?;
            add(
                &mut self.enthalpy[side],
                row.stream_enthalpy_into_gas_j[side],
            )?;
            add(
                &mut self.absolute_enthalpy[side],
                row.stream_enthalpy_into_gas_j[side].abs(),
            )?;
            add(&mut self.heat[side], row.thermostat_heat_into_gas_j[side])?;
            add(
                &mut self.absolute_heat[side],
                row.thermostat_heat_into_gas_j[side].abs(),
            )?;
            // EXACT grouping of the original isothermal_stream/boundary owner.
            // Keep the paired relative thermal input before summing streams.
            let flow_work = finite((amount * mass) * (R * temperature / mass))?;
            add(
                &mut self.relative[side],
                finite(flow_work + row.thermostat_heat_into_gas_j[side])?,
            )?;
        }
        Ok(())
    }
}

impl<'a> IsothermalGasMembraneMixture<'a> {
    pub fn adaptive_workspace(&self) -> Result<MembraneMixtureAdvanceWorkspace> {
        Ok(MembraneMixtureAdvanceWorkspace {
            full: self.workspace()?,
            half: self.workspace()?,
        })
    }

    /// Complete a positive duration using the SAME original finite-gas/film
    /// steps. Every rejected full/fine probe spends work but contributes no
    /// accepted inventory, boundary stream or thermostat heat. Actual source
    /// support is rechecked by the native owner at every probe; no coefficient,
    /// phase storage or source-selection fallback is invented.
    ///
    /// Prepared native workspaces are reused, but native immutable candidates,
    /// gas reconstruction and original source receipts allocate. Complete cell
    /// count and probe work are bounded explicitly. One-species catalogs use
    /// this same API. Failure leaves this original owner and clock unchanged.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        controls: MembraneMixtureAdvanceControls,
        work: &mut MembraneMixtureAdvanceWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsothermalGasMembraneMixtureAdvanceCandidate<'a>> {
        poll(&mut cancelled)?;
        validate(self, duration_s, controls)?;
        let start = self.elapsed_s();
        let target = start + duration_s;
        let mut state = self.clone();
        let mut ledgers: Vec<Ledger> = self.channels.iter().map(|_| Ledger::default()).collect();
        let (mut accepted, mut rejected, mut evaluations) = (0, 0, 0);
        let (mut minimum, mut maximum, mut maximum_error) = (f64::INFINITY, 0.0_f64, 0.0_f64);
        let mut step = controls.initial_step_s;
        let mut last = None;
        while state.elapsed_s() < target {
            poll(&mut cancelled)?;
            let time = state.elapsed_s();
            let planned_end = if step >= target - time {
                target
            } else {
                time + step
            };
            // Do not leave a one-ULP tail which cannot contain two physical
            // half steps. Include it in THIS proposal and recompute both
            // durations before any native solve; never snap an evolved clock.
            let end = if planned_end < target && planned_end.next_up() >= target {
                target
            } else {
                planned_end
            };
            let span = end - time;
            let middle = time + 0.5 * span;
            if !(time < middle && middle < end) {
                return Err(Error::TimeResolution {
                    time_s: time,
                    target_s: end,
                });
            }
            let (h1, h2) = (middle - time, end - middle);
            let coarse = probe(
                &state,
                span,
                &mut work.full,
                &mut evaluations,
                controls,
                &mut cancelled,
            )?
            .into_parts()
            .0;
            let (half, first) = probe(
                &state,
                h1,
                &mut work.half,
                &mut evaluations,
                controls,
                &mut cancelled,
            )?
            .into_parts();
            check_clock(&half, middle)?;
            let (fine, second) = probe(
                &half,
                h2,
                &mut work.half,
                &mut evaluations,
                controls,
                &mut cancelled,
            )?
            .into_parts();
            check_clock(&coarse, end)?;
            check_clock(&fine, end)?;
            let (u, v) = (h1 / span, h2 / span);
            let factor = (u * u + v * v) / (2.0 * u * v);
            let error = factor * distance(&state, &coarse, &fine, controls, &mut cancelled)?;
            if !error.is_finite() || error > 1.0 {
                let smaller = (0.5 * span).max(controls.minimum_step_s);
                if smaller >= span {
                    return Err(Error::Accuracy {
                        time_s: time,
                        step_s: span,
                        error_ratio: error,
                    });
                }
                rejected += 1;
                step = smaller;
                continue;
            }
            for (i, ledger) in ledgers.iter_mut().enumerate() {
                poll(&mut cancelled)?;
                ledger.add_step(
                    &first.species[i],
                    self.channels[i].molar_mass,
                    self.temperature_k(),
                )?;
                ledger.add_step(
                    &second.species[i],
                    self.channels[i].molar_mass,
                    self.temperature_k(),
                )?;
            }
            accepted += 1;
            minimum = minimum.min(h1.min(h2));
            maximum = maximum.max(h1.max(h2));
            maximum_error = maximum_error.max(error);
            state = fine;
            last = Some(second);
            step = if error <= 0.25 {
                (2.0 * span).min(controls.maximum_step_s)
            } else {
                span
            };
            step = step.max(controls.minimum_step_s);
        }
        check_clock(&state, target)?;
        let mut report = audit(
            self,
            &state,
            ledgers,
            controls,
            last.expect("positive completed interval contains native half steps"),
            &mut cancelled,
        )?;
        report.interval_s = [start, target];
        report.accepted_intervals = accepted;
        report.accepted_substeps = 2 * accepted;
        report.rejected_intervals = rejected;
        report.evaluations = evaluations;
        report.minimum_accepted_substep_s = minimum;
        report.maximum_accepted_substep_s = maximum;
        report.maximum_accepted_error_ratio = maximum_error;
        poll(&mut cancelled)?;
        Ok(IsothermalGasMembraneMixtureAdvanceCandidate { state, report })
    }
}

fn validate(
    state: &IsothermalGasMembraneMixture<'_>,
    duration: f64,
    c: MembraneMixtureAdvanceControls,
) -> Result<()> {
    let target = state.elapsed_s() + duration;
    if !positive(duration)
        || !target.is_finite()
        || target <= state.elapsed_s()
        || !positive(c.minimum_step_s)
        || !positive(c.maximum_step_s)
        || !positive(c.initial_step_s)
        || c.minimum_step_s > c.initial_step_s
        || c.initial_step_s > c.maximum_step_s
        || !(3..=100_000).contains(&c.maximum_evaluations)
        || !positive(c.absolute_pressure_tolerance_pa)
        || !positive(c.absolute_amount_tolerance_mol)
        || !c.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&c.relative_tolerance)
        || !positive(c.species_balance_tolerance_mol)
        || !positive(c.mass_balance_tolerance_kg)
        || !positive(c.gas_energy_tolerance_j)
    {
        return Err(Error::Invalid(
            "positive forward duration, explicit accuracy/balance limits and bounded native work required",
        ));
    }
    let cells = state.channels.iter().try_fold(0usize, |n, channel| {
        n.checked_add(channel.finite.membrane().model().cell_count())
    });
    let maximum =
        super::super::MAX_TRANSIENT_MEMBRANE_CELLS * super::MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1;
    if c.maximum_total_cells == 0
        || c.maximum_total_cells > maximum
        || cells.is_none_or(|n| n > c.maximum_total_cells)
    {
        return Err(Error::Invalid(
            "complete shared film profiles exceed the admitted total cell budget",
        ));
    }
    check_clock(state, state.elapsed_s())
}
fn check_clock(state: &IsothermalGasMembraneMixture<'_>, target: f64) -> Result<()> {
    if state
        .channels
        .iter()
        .any(|channel| channel.finite.membrane().elapsed_s().to_bits() != target.to_bits())
    {
        Err(Error::TimeResolution {
            time_s: state.elapsed_s(),
            target_s: target,
        })
    } else {
        Ok(())
    }
}
fn probe<'a>(
    state: &IsothermalGasMembraneMixture<'a>,
    dt: f64,
    work: &mut IsothermalGasMembraneMixtureWorkspace,
    calls: &mut usize,
    controls: MembraneMixtureAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<super::IsothermalGasMembraneMixtureCandidate<'a>> {
    poll(cancelled)?;
    if *calls >= controls.maximum_evaluations {
        return Err(Error::WorkLimit {
            time_s: state.elapsed_s(),
            evaluations: *calls,
        });
    }
    *calls += 1;
    Ok(state.step_candidate(dt, work, cancelled)?)
}
fn gas_amounts(state: &IsothermalGasMembraneMixture<'_>) -> Result<Vec<[f64; 2]>> {
    let mass = [
        inventory(state.gases[0].model(), state.gases[0].mass_kg()).map_err(NativeError::from)?,
        inventory(state.gases[1].model(), state.gases[1].mass_kg()).map_err(NativeError::from)?,
    ];
    state
        .channels
        .iter()
        .map(|channel| {
            Ok([
                finite(
                    mass[0].get(channel.species()).copied().unwrap_or(0.0) / channel.molar_mass,
                )?,
                finite(
                    mass[1].get(channel.species()).copied().unwrap_or(0.0) / channel.molar_mass,
                )?,
            ])
        })
        .collect()
}
fn scaled(old: f64, a: f64, b: f64, absolute: f64, relative: f64) -> Result<f64> {
    let allowance = finite(absolute + relative * (a - old).abs().max((b - old).abs()))?;
    Ok((a - b).abs() / allowance)
}
fn distance(
    old: &IsothermalGasMembraneMixture<'_>,
    a: &IsothermalGasMembraneMixture<'_>,
    b: &IsothermalGasMembraneMixture<'_>,
    c: MembraneMixtureAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let (old_gas, a_gas, b_gas) = (gas_amounts(old)?, gas_amounts(a)?, gas_amounts(b)?);
    let mut error = 0.0_f64;
    for i in 0..old.channels.len() {
        poll(cancelled)?;
        let (old, a, b) = (
            &old.channels[i].finite,
            &a.channels[i].finite,
            &b.channels[i].finite,
        );
        for side in 0..2 {
            error = error.max(scaled(
                old_gas[i][side],
                a_gas[i][side],
                b_gas[i][side],
                c.absolute_amount_tolerance_mol,
                c.relative_tolerance,
            )?);
            error = error.max(scaled(
                old.reservoir_pressures_pa()[side],
                a.reservoir_pressures_pa()[side],
                b.reservoir_pressures_pa()[side],
                c.absolute_pressure_tolerance_pa,
                c.relative_tolerance,
            )?);
        }
        for (j, ((&old, &a), &b)) in old
            .membrane()
            .partial_pressures_pa()
            .iter()
            .zip(a.membrane().partial_pressures_pa())
            .zip(b.membrane().partial_pressures_pa())
            .enumerate()
        {
            if j % 128 == 0 {
                poll(cancelled)?;
            }
            error = error.max(scaled(
                old,
                a,
                b,
                c.absolute_pressure_tolerance_pa,
                c.relative_tolerance,
            )?);
        }
    }
    Ok(error)
}

fn audit(
    old: &IsothermalGasMembraneMixture<'_>,
    next: &IsothermalGasMembraneMixture<'_>,
    ledgers: Vec<Ledger>,
    c: MembraneMixtureAdvanceControls,
    last_substep: IsothermalGasMembraneMixtureStep,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MembraneMixtureAdvanceReport> {
    let mut rows = Vec::with_capacity(old.channels.len());
    let mut storage_mass = Sum::default();
    let mut enthalpy = [Sum::default(), Sum::default()];
    let mut heat = [Sum::default(), Sum::default()];
    let mut relative = [Sum::default(), Sum::default()];
    for (i, ledger) in ledgers.into_iter().enumerate() {
        poll(cancelled)?;
        let channel = &old.channels[i];
        let id = channel.species();
        let model = channel.finite.membrane().model();
        let before = channel.finite.membrane().partial_pressures_pa();
        let after = next.channels[i].finite.membrane().partial_pressures_pa();
        let mut storage = Sum::default();
        for (layer, range) in model.layers().iter().zip(model.layer_cell_ranges()) {
            // Same original S*dx then area grouping as the native cell owner.
            let capacity = finite(
                (layer.solubility.value_si() * (layer.thickness_m / layer.cells as f64))
                    * model.area_m2(),
            )?;
            for j in range.clone() {
                if j % 128 == 0 {
                    poll(cancelled)?;
                }
                add(&mut storage, finite(capacity * (after[j] - before[j]))?)?;
            }
        }
        let mut changes = [0.0; 2];
        let mut gas_residual = [0.0; 2];
        for side in 0..2 {
            let a = old.gases[side]
                .species_mass_kg(id)
                .map_err(NativeError::from)?;
            let b = next.gases[side]
                .species_mass_kg(id)
                .map_err(NativeError::from)?;
            changes[side] = finite((b - a) / channel.molar_mass)?;
            gas_residual[side] = finite(changes[side] - ledger.amount[side].value)?;
            species_gate(
                id,
                if side == 0 {
                    MembraneMixtureBalanceAccount::LeftGas
                } else {
                    MembraneMixtureBalanceAccount::RightGas
                },
                gas_residual[side],
                c.species_balance_tolerance_mol,
            )?;
            add(&mut enthalpy[side], ledger.enthalpy[side].value)?;
            add(&mut heat[side], ledger.heat[side].value)?;
            add(&mut relative[side], ledger.relative[side].value)?;
        }
        let mut film = Sum::default();
        add(&mut film, storage.value)?;
        add(&mut film, ledger.amount[0].value)?;
        add(&mut film, ledger.amount[1].value)?;
        species_gate(
            id,
            MembraneMixtureBalanceAccount::FilmBoundary,
            film.value,
            c.species_balance_tolerance_mol,
        )?;
        let mut closed = Sum::default();
        add(&mut closed, storage.value)?;
        add(&mut closed, changes[0])?;
        add(&mut closed, changes[1])?;
        species_gate(
            id,
            MembraneMixtureBalanceAccount::ClosedSpecies,
            closed.value,
            c.species_balance_tolerance_mol,
        )?;
        add(
            &mut storage_mass,
            finite(storage.value * channel.molar_mass)?,
        )?;
        rows.push(MembraneMixtureIntegratedSpecies {
            species: id.clone(),
            amount_into_gas_mol: ledger.amount.map(|x| x.value),
            gross_amount_into_gas_mol: ledger.incoming.map(|x| x.value),
            gross_amount_out_of_gas_mol: ledger.outgoing.map(|x| x.value),
            enthalpy_into_gas_j: ledger.enthalpy.map(|x| x.value),
            gross_absolute_enthalpy_j: ledger.absolute_enthalpy.map(|x| x.value),
            thermostat_heat_into_gas_j: ledger.heat.map(|x| x.value),
            gross_absolute_thermostat_heat_j: ledger.absolute_heat.map(|x| x.value),
            actual_gas_amount_changes_mol: changes,
            actual_sorbed_amount_change_mol: storage.value,
            gas_boundary_residuals_mol: gas_residual,
            film_boundary_residual_mol: film.value,
            closed_species_residual_mol: closed.value,
        });
    }
    let gas_changes = [
        next.gases[0].mass_kg() - old.gases[0].mass_kg(),
        next.gases[1].mass_kg() - old.gases[1].mass_kg(),
    ];
    let mut total = Sum::default();
    add(&mut total, storage_mass.value)?;
    add(&mut total, gas_changes[0])?;
    add(&mut total, gas_changes[1])?;
    let direct = finite(next.total_mass_kg()? - old.total_mass_kg()?)?;
    for residual in [total.value, direct] {
        if residual.abs() > c.mass_balance_tolerance_kg {
            return Err(Error::MassBalance {
                residual_kg: residual,
                tolerance_kg: c.mass_balance_tolerance_kg,
            });
        }
    }
    let mut energy = [0.0; 2];
    for side in 0..2 {
        poll(cancelled)?;
        let actual = next.gases[side]
            .model()
            .thermodynamics()
            .specific_internal_energy_change(
                Temperature::new(old.gases[side].state().state().temperature),
                Temperature::new(next.gases[side].state().state().temperature),
            )
            .map_err(super::super::source_error)
            .map_err(NativeError::from)?
            .value();
        energy[side] = finite(next.gases[side].mass_kg() * actual - relative[side].value)?;
        if energy[side].abs() > c.gas_energy_tolerance_j {
            return Err(Error::GasEnergyBalance {
                reservoir: side,
                residual_j: energy[side],
                tolerance_j: c.gas_energy_tolerance_j,
            });
        }
    }
    Ok(MembraneMixtureAdvanceReport {
        interval_s: [0.0; 2],
        accepted_intervals: 0,
        accepted_substeps: 0,
        rejected_intervals: 0,
        evaluations: 0,
        minimum_accepted_substep_s: 0.0,
        maximum_accepted_substep_s: 0.0,
        maximum_accepted_error_ratio: 0.0,
        species: rows,
        actual_gas_mass_changes_kg: gas_changes,
        actual_sorbed_mass_change_kg: storage_mass.value,
        total_mass_residual_kg: total.value,
        direct_total_mass_residual_kg: direct,
        enthalpy_into_gas_j: enthalpy.map(|x| x.value),
        thermostat_heat_into_gas_j: heat.map(|x| x.value),
        source_relative_gas_energy_residual_j: energy,
        last_substep,
    })
}
fn species_gate(
    id: &SpeciesId,
    account: MembraneMixtureBalanceAccount,
    residual: f64,
    tolerance: f64,
) -> Result<()> {
    if !residual.is_finite() || residual.abs() > tolerance {
        Err(Error::SpeciesBalance {
            species: id.clone(),
            account,
            residual_mol: residual,
            tolerance_mol: tolerance,
        })
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Invalid(
            "unrepresentable adaptive membrane arithmetic",
        ))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn add(sum: &mut Sum, value: f64) -> Result<()> {
    Ok(sum.add(value).map_err(NativeError::from)?)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
