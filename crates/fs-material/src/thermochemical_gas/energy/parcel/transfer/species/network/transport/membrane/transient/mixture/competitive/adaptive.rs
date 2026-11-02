//! Accepted competitive breakthrough histories with original finite gas owners.
//!
//! The joint competitive time driver invokes the source-qualified gas candidate
//! on every probe. Original gas reconstruction, mobile/site applicability and
//! all native numerical budgets remain hard requirements. Cloned proposals
//! carry small stream ledgers; only the accepted half-step path survives.
//! Final accounts independently compare actual original-to-final gas species,
//! mobile/trapped stores, total mass and source-relative gas thermal change.
//! These are gas-side energy accounts at externally held temperature, not a
//! sorbed caloric store, heat-of-sorption model or whole-system energy closure.

use super::*;
use super::super::adaptive::{
    MembraneMixtureAdvanceControls, MembraneMixtureBalanceAccount,
    MembraneMixtureIntegratedSpecies,
};
use super::super::super::slab::{
    adaptive::{MembraneDiffusionAdvanceControls, MembraneDiffusionAdvanceError, MembraneDiffusionAdvanceWorkspace},
    competitive::adaptive::{self as driver, CompetitiveMembraneAdvanceControls,
        CompetitiveMembraneAdvanceReport, CompetitiveNative},
};
use fs_qty::Temperature;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

type NativeResult<T> = core::result::Result<T, IsothermalCompetitiveGasError>;
type Result<T> = core::result::Result<T, IsothermalCompetitiveGasAdvanceError>;

#[derive(Clone, Debug, PartialEq)]
pub enum IsothermalCompetitiveGasAdvanceError {
    Invalid(&'static str),
    Native(IsothermalCompetitiveGasError),
    Diffusion(MembraneDiffusionAdvanceError<IsothermalCompetitiveGasError>),
    SpeciesBalance { species: SpeciesId, account: MembraneMixtureBalanceAccount,
        residual_mol: f64, tolerance_mol: f64 },
    MassBalance { residual_kg: f64, tolerance_kg: f64 },
    GasEnergyBalance { reservoir: usize, residual_j: f64, tolerance_j: f64 },
    Cancelled,
}
use IsothermalCompetitiveGasAdvanceError as Error;
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive source-qualified competitive gas membrane: {self:?}")
    }
}
impl std::error::Error for Error {}
impl From<IsothermalCompetitiveGasError> for Error {
    fn from(error: IsothermalCompetitiveGasError) -> Self {
        if error == IsothermalCompetitiveGasError::Cancelled { Self::Cancelled } else { Self::Native(error) }
    }
}
impl From<IsothermalGasMembraneError> for Error {
    fn from(error: IsothermalGasMembraneError) -> Self { IsothermalCompetitiveGasError::from(error).into() }
}
impl From<MembraneDiffusionAdvanceError<IsothermalCompetitiveGasError>> for Error {
    fn from(error: MembraneDiffusionAdvanceError<IsothermalCompetitiveGasError>) -> Self {
        match error {
            MembraneDiffusionAdvanceError::Cancelled => Self::Cancelled,
            MembraneDiffusionAdvanceError::Native(error) => Self::from(error),
            other => Self::Diffusion(other),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct IsothermalCompetitiveGasAdvanceReport {
    /// Original coupled diffusion timeline, per-species mobile/trapped accounts
    /// and accepted-only signed/gross face transfers.
    pub diffusion: CompetitiveMembraneAdvanceReport,
    /// Canonical species order; `actual_sorbed_amount_change_mol` contains BOTH
    /// mobile Henry storage and finite-site storage, separately retained above.
    pub species: Vec<MembraneMixtureIntegratedSpecies>,
    pub actual_gas_mass_changes_kg: [f64; 2],
    pub actual_mobile_mass_change_kg: f64,
    pub actual_trapped_mass_change_kg: f64,
    pub total_mass_residual_kg: f64,
    pub direct_total_mass_residual_kg: f64,
    pub enthalpy_into_gas_j: [f64; 2],
    pub thermostat_heat_into_gas_j: [f64; 2],
    pub source_relative_gas_energy_residual_j: [f64; 2],
    pub last_substep: IsothermalCompetitiveGasStep,
}

#[derive(Clone, Debug)]
pub struct IsothermalCompetitiveGasAdvanceCandidate<'a> {
    state: IsothermalCompetitiveGasMembrane<'a>,
    report: IsothermalCompetitiveGasAdvanceReport,
}
impl<'a> IsothermalCompetitiveGasAdvanceCandidate<'a> {
    pub const fn state(&self) -> &IsothermalCompetitiveGasMembrane<'a> { &self.state }
    pub const fn report(&self) -> &IsothermalCompetitiveGasAdvanceReport { &self.report }
    pub fn into_parts(self) -> (IsothermalCompetitiveGasMembrane<'a>, IsothermalCompetitiveGasAdvanceReport) {
        (self.state, self.report)
    }
}

#[derive(Default)]
struct Ledger {
    amount: [Sum; 2], incoming: [Sum; 2], outgoing: [Sum; 2],
    enthalpy: [Sum; 2], absolute_enthalpy: [Sum; 2],
    heat: [Sum; 2], absolute_heat: [Sum; 2], relative: [Sum; 2],
}
impl Clone for Ledger {
    fn clone(&self) -> Self {
        let copy = |s: &[Sum; 2]| core::array::from_fn(|i| Sum { value: s[i].value, correction: s[i].correction });
        Self { amount: copy(&self.amount), incoming: copy(&self.incoming), outgoing: copy(&self.outgoing),
            enthalpy: copy(&self.enthalpy), absolute_enthalpy: copy(&self.absolute_enthalpy),
            heat: copy(&self.heat), absolute_heat: copy(&self.absolute_heat), relative: copy(&self.relative) }
    }
}
impl Ledger {
    fn add(&mut self, amount: [f64; 2], row: &CompetitiveGasSpeciesAccount,
        molar_mass: f64, temperature: f64,
    ) -> NativeResult<()> {
        for side in 0..2 {
            self.amount[side].add(amount[side])?;
            self.incoming[side].add(amount[side].max(0.0))?;
            self.outgoing[side].add((-amount[side]).max(0.0))?;
            self.enthalpy[side].add(row.stream_enthalpy_into_gas_j[side])?;
            self.absolute_enthalpy[side].add(row.stream_enthalpy_into_gas_j[side].abs())?;
            self.heat[side].add(row.thermostat_heat_into_gas_j[side])?;
            self.absolute_heat[side].add(row.thermostat_heat_into_gas_j[side].abs())?;
            // Original isothermal_stream grouping. Pair thermal work and its
            // thermostat offset BEFORE summing across time/species or adding
            // any arbitrarily large formation-reference enthalpies.
            let flow_work = finite((amount[side] * molar_mass) * (R * temperature / molar_mass))?;
            self.relative[side].add(finite(flow_work + row.thermostat_heat_into_gas_j[side])?)?;
        }
        Ok(())
    }
}

#[derive(Clone)]
struct Tracked<'a> {
    state: IsothermalCompetitiveGasMembrane<'a>,
    ledgers: Vec<Ledger>,
    controls: MembraneMixtureAdvanceControls,
    window: f64,
    minimum_molar_mass: f64,
    last: Option<IsothermalCompetitiveGasStep>,
}
impl CompetitiveNative for Tracked<'_> {
    type Workspace = CompetitiveMembraneWorkspace;
    type Error = IsothermalCompetitiveGasError;
    fn competition(&self) -> &CompetitiveMembraneDiffusionState { self.state.membrane() }
    fn finite_reservoirs(&self) -> Option<([f64; 2], &[[f64; 2]])> {
        Some((self.state.finite.capacities().map(|c| c.value_si()), self.state.finite.reservoir_partial_pressures_pa()))
    }
    fn is_cancelled(error: &Self::Error) -> bool { *error == IsothermalCompetitiveGasError::Cancelled }
    fn additional_distance(
        &self,
        coarse: &Self,
        fine: &Self,
        controls: MembraneDiffusionAdvanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<f64, MembraneDiffusionAdvanceError<Self::Error>> {
        let native = |e: GasEnergyError| MembraneDiffusionAdvanceError::Native(IsothermalCompetitiveGasError::from(e));
        let mut error = 0.0_f64;
        for side in 0..2 {
            if cancelled() { return Err(MembraneDiffusionAdvanceError::Cancelled); }
            // Each original composition is reconstructed once for this error
            // comparison. Never substitute a pressure-capacity inventory for
            // the actual gas species account, including absent species.
            let old = inventory(self.state.gases[side].model(), self.state.gases[side].mass_kg()).map_err(native)?;
            let a = inventory(coarse.state.gases[side].model(), coarse.state.gases[side].mass_kg()).map_err(native)?;
            let b = inventory(fine.state.gases[side].model(), fine.state.gases[side].mass_kg()).map_err(native)?;
            for channel in self.state.channels.iter() {
                let id = &channel.bindings.species;
                let old = old.get(id).copied().unwrap_or(0.0);
                let a = a.get(id).copied().unwrap_or(0.0);
                let b = b.get(id).copied().unwrap_or(0.0);
                let difference = (a-b) / channel.molar_mass;
                let change = (a-old).abs().max((b-old).abs()) / channel.molar_mass;
                let allowance = controls.absolute_amount_tolerance_mol + controls.relative_tolerance*change;
                if !difference.is_finite() || !positive(allowance) {
                    return Err(MembraneDiffusionAdvanceError::Invalid("actual gas temporal amount scale is unrepresentable"));
                }
                error = error.max(difference.abs()/allowance);
            }
        }
        Ok(error)
    }
    fn competitive_step(&self, _: &[[f64; 2]], duration: f64, requested: CompetitiveMembraneSolve,
        work: &mut Self::Workspace, cancelled: &mut impl FnMut() -> bool,
    ) -> NativeResult<(Self, CompetitiveMembraneStep)> {
        poll(cancelled)?;
        let mut probe = self.state.clone();
        let original = self.state.controls;
        // Every original native bound is a hard ceiling. Tightened controls
        // belong only to this proposal, never to the returned model definition.
        probe.controls = CompetitiveMembraneSolve {
            balance: super::super::super::slab::MembraneDiffusionSolve {
                absolute_amount_tolerance_mol: requested.balance.absolute_amount_tolerance_mol.min(original.balance.absolute_amount_tolerance_mol),
                relative_amount_tolerance: requested.balance.relative_amount_tolerance.min(original.balance.relative_amount_tolerance),
            },
            maximum_iterations: requested.maximum_iterations.min(original.maximum_iterations),
            maximum_backtracks: requested.maximum_backtracks.min(original.maximum_backtracks),
            maximum_linear_sweeps: requested.maximum_linear_sweeps.min(original.maximum_linear_sweeps),
            linear_relative_tolerance: requested.linear_relative_tolerance.min(original.linear_relative_tolerance),
        };
        // Reserve a duration-weighted part of fixed whole-window gas/species
        // budgets, beneath temporal amount accuracy and original gas limits.
        let share = duration / self.window;
        let mass = ((0.01 * self.controls.mass_balance_tolerance_kg) * share / (2.0 * self.state.channels.len() as f64))
            .min((0.01 * self.controls.species_balance_tolerance_mol) * share * self.minimum_molar_mass)
            .min((0.01 * self.controls.absolute_amount_tolerance_mol) * self.minimum_molar_mass);
        let energy = (0.01 * self.controls.gas_energy_tolerance_j) * share;
        if !positive(mass) || !positive(energy) { return Err(invalid("adaptive gas algebraic budget underflowed")); }
        probe.gas_controls.mass_tolerance_kg = probe.gas_controls.mass_tolerance_kg.min(mass);
        for gas in &mut probe.gas_controls.gas {
            gas.mass_tolerance_kg = gas.mass_tolerance_kg.min(mass);
            gas.energy_tolerance_j = gas.energy_tolerance_j.min(energy);
        }
        let (mut state, step) = probe.step_candidate(duration, work, &mut *cancelled)?.into_parts();
        state.controls = self.state.controls;
        state.gas_controls = self.state.gas_controls;
        let mut ledgers = self.ledgers.clone();
        for (s, row) in step.species.iter().enumerate() {
            poll(cancelled)?;
            let report = &step.diffusion.species[s].1;
            ledgers[s].add([-report.boundary_transfers_mol[0], report.boundary_transfers_mol[1]],
                row, self.state.channels[s].molar_mass, self.state.temperature_k())?;
        }
        let diffusion = step.diffusion.clone();
        Ok((Self { state, ledgers, controls: self.controls, window: self.window,
            minimum_molar_mass: self.minimum_molar_mass, last: Some(step) }, diffusion))
    }
}

fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn add(sum: &mut Sum, value: f64) -> Result<()> { sum.add(value).map_err(|e| Error::Native(e.into())) }
fn checked(value: f64) -> Result<f64> { Ok(finite(value)?) }
fn species_gate(id: &SpeciesId, account: MembraneMixtureBalanceAccount, residual: f64, tolerance: f64) -> Result<()> {
    if !residual.is_finite() || residual.abs() > tolerance {
        Err(Error::SpeciesBalance { species: id.clone(), account, residual_mol: residual, tolerance_mol: tolerance })
    } else { Ok(()) }
}

fn audit(old: &IsothermalCompetitiveGasMembrane<'_>, next: &Tracked<'_>,
    diffusion: CompetitiveMembraneAdvanceReport, cancelled: &mut impl FnMut() -> bool,
) -> Result<IsothermalCompetitiveGasAdvanceReport> {
    let c = next.controls;
    let before_inventory = [
        inventory(old.gases[0].model(), old.gases[0].mass_kg()).map_err(IsothermalCompetitiveGasError::from)?,
        inventory(old.gases[1].model(), old.gases[1].mass_kg()).map_err(IsothermalCompetitiveGasError::from)?,
    ];
    let after_inventory = [
        inventory(next.state.gases[0].model(), next.state.gases[0].mass_kg()).map_err(IsothermalCompetitiveGasError::from)?,
        inventory(next.state.gases[1].model(), next.state.gases[1].mass_kg()).map_err(IsothermalCompetitiveGasError::from)?,
    ];
    let mut rows = Vec::with_capacity(old.channels.len());
    let (mut mobile_mass, mut trapped_mass) = (Sum::default(), Sum::default());
    let (mut enthalpy, mut heat, mut relative) = ([Sum::default(), Sum::default()], [Sum::default(), Sum::default()], [Sum::default(), Sum::default()]);
    for (s, ledger) in next.ledgers.iter().enumerate() {
        poll(cancelled)?;
        let channel = &old.channels[s];
        let id = &channel.bindings.species;
        let record = &diffusion.species[s].1;
        let storage = checked(record.mobile_change_mol + record.trapped_change_mol)?;
        let (mut changes, mut residuals) = ([0.0; 2], [0.0; 2]);
        let expected_stream = [-record.transport.face_amounts_mol[0], *record.transport.face_amounts_mol.last().expect("membrane exterior face")];
        for side in 0..2 {
            let account = if side == 0 { MembraneMixtureBalanceAccount::LeftGas } else { MembraneMixtureBalanceAccount::RightGas };
            let a = before_inventory[side].get(id).copied().unwrap_or(0.0);
            let b = after_inventory[side].get(id).copied().unwrap_or(0.0);
            changes[side] = checked((b-a) / channel.molar_mass)?;
            residuals[side] = checked(changes[side] - ledger.amount[side].value)?;
            species_gate(id, account, residuals[side], c.species_balance_tolerance_mol)?;
            species_gate(id, account, checked(ledger.amount[side].value - expected_stream[side])?, c.species_balance_tolerance_mol)?;
            add(&mut enthalpy[side], ledger.enthalpy[side].value)?;
            add(&mut heat[side], ledger.heat[side].value)?;
            add(&mut relative[side], ledger.relative[side].value)?;
        }
        let (mut film, mut closed) = (Sum::default(), Sum::default());
        for value in [record.mobile_change_mol, record.trapped_change_mol, ledger.amount[0].value, ledger.amount[1].value] { add(&mut film, value)?; }
        for value in [record.mobile_change_mol, record.trapped_change_mol, changes[0], changes[1]] { add(&mut closed, value)?; }
        species_gate(id, MembraneMixtureBalanceAccount::FilmBoundary, film.value, c.species_balance_tolerance_mol)?;
        species_gate(id, MembraneMixtureBalanceAccount::ClosedSpecies, closed.value, c.species_balance_tolerance_mol)?;
        add(&mut mobile_mass, checked(record.mobile_change_mol * channel.molar_mass)?)?;
        add(&mut trapped_mass, checked(record.trapped_change_mol * channel.molar_mass)?)?;
        rows.push(MembraneMixtureIntegratedSpecies {
            species: id.clone(), amount_into_gas_mol: ledger.amount.each_ref().map(|s| s.value),
            gross_amount_into_gas_mol: ledger.incoming.each_ref().map(|s| s.value),
            gross_amount_out_of_gas_mol: ledger.outgoing.each_ref().map(|s| s.value),
            enthalpy_into_gas_j: ledger.enthalpy.each_ref().map(|s| s.value),
            gross_absolute_enthalpy_j: ledger.absolute_enthalpy.each_ref().map(|s| s.value),
            thermostat_heat_into_gas_j: ledger.heat.each_ref().map(|s| s.value),
            gross_absolute_thermostat_heat_j: ledger.absolute_heat.each_ref().map(|s| s.value),
            actual_gas_amount_changes_mol: changes, actual_sorbed_amount_change_mol: storage,
            gas_boundary_residuals_mol: residuals, film_boundary_residual_mol: film.value,
            closed_species_residual_mol: closed.value,
        });
    }
    let changes = core::array::from_fn(|side| next.state.gases[side].mass_kg()-old.gases[side].mass_kg());
    let mut total = Sum::default();
    for value in [mobile_mass.value, trapped_mass.value, changes[0], changes[1]] { add(&mut total, value)?; }
    let direct = checked(next.state.total_mass_kg()? - old.total_mass_kg()?)?;
    for residual in [total.value, direct] {
        if residual.abs() > c.mass_balance_tolerance_kg {
            return Err(Error::MassBalance { residual_kg: residual, tolerance_kg: c.mass_balance_tolerance_kg });
        }
    }
    let mut energy = [0.0; 2];
    for side in 0..2 {
        poll(cancelled)?;
        let change = next.state.gases[side].model().thermodynamics().specific_internal_energy_change(
            Temperature::new(old.gases[side].state().state().temperature),
            Temperature::new(next.state.gases[side].state().state().temperature),
        ).map_err(super::super::super::source_error).map_err(IsothermalCompetitiveGasError::from)?.value();
        energy[side] = checked(next.state.gases[side].mass_kg() * change - relative[side].value)?;
        if energy[side].abs() > c.gas_energy_tolerance_j {
            return Err(Error::GasEnergyBalance { reservoir: side, residual_j: energy[side], tolerance_j: c.gas_energy_tolerance_j });
        }
    }
    Ok(IsothermalCompetitiveGasAdvanceReport { diffusion, species: rows,
        actual_gas_mass_changes_kg: changes, actual_mobile_mass_change_kg: mobile_mass.value,
        actual_trapped_mass_change_kg: trapped_mass.value,
        total_mass_residual_kg: total.value, direct_total_mass_residual_kg: direct,
        enthalpy_into_gas_j: enthalpy.map(|s| s.value), thermostat_heat_into_gas_j: heat.map(|s| s.value),
        source_relative_gas_energy_residual_j: energy,
        last_substep: next.last.clone().expect("completed positive actual-gas window") })
}

impl<'a> IsothermalCompetitiveGasMembrane<'a> {
    /// The byte ceiling covers both retained native competitive workspaces.
    pub fn adaptive_workspace(&self, maximum_bytes: usize) -> Result<MembraneDiffusionAdvanceWorkspace<CompetitiveMembraneWorkspace>> {
        self.finite.adaptive_workspace(maximum_bytes).map_err(IsothermalCompetitiveGasError::from).map_err(Error::from)
    }
    /// Source-qualified original-gas probes on one competitive timeline.
    /// Native work limits and source/gas controls are never loosened; original
    /// definitions are restored after each temporarily tightened probe.
    pub fn advance_candidate(&self, duration_s: f64, controls: MembraneMixtureAdvanceControls,
        work: &mut MembraneDiffusionAdvanceWorkspace<CompetitiveMembraneWorkspace>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsothermalCompetitiveGasAdvanceCandidate<'a>> {
        poll(&mut cancelled)?;
        if !positive(controls.mass_balance_tolerance_kg) || !positive(controls.gas_energy_tolerance_j) {
            return Err(Error::Invalid("positive fixed whole-window mass and gas-energy allowances required"));
        }
        let count = self.channels.len().checked_mul(self.membrane().mobile(0).expect("admitted first species").model().cell_count());
        let maximum = crate::adsorption::competitive::MAX_COMPETITIVE_SPECIES * super::super::super::MAX_TRANSIENT_MEMBRANE_CELLS;
        if controls.maximum_total_cells == 0 || controls.maximum_total_cells > maximum || count.is_none_or(|n| n > controls.maximum_total_cells) {
            return Err(Error::Invalid("complete competitive profiles exceed the explicit total-cell budget"));
        }
        let tracked = Tracked { state: self.clone(), ledgers: self.channels.iter().map(|_| Ledger::default()).collect(), controls,
            window: (self.elapsed_s() + duration_s) - self.elapsed_s(),
            minimum_molar_mass: self.channels.iter().map(|c| c.molar_mass).fold(f64::INFINITY, f64::min), last: None };
        let temporal = MembraneDiffusionAdvanceControls {
            step: self.controls.balance, initial_step_s: controls.initial_step_s,
            minimum_step_s: controls.minimum_step_s, maximum_step_s: controls.maximum_step_s,
            maximum_evaluations: controls.maximum_evaluations,
            absolute_pressure_tolerance_pa: controls.absolute_pressure_tolerance_pa,
            absolute_amount_tolerance_mol: controls.absolute_amount_tolerance_mol,
            relative_tolerance: controls.relative_tolerance,
            amount_balance_tolerance_mol: controls.species_balance_tolerance_mol,
        };
        let (next, diffusion) = driver::advance(&tracked, &[], duration_s, CompetitiveMembraneAdvanceControls {
            temporal, maximum_iterations: self.controls.maximum_iterations, maximum_backtracks: self.controls.maximum_backtracks,
            maximum_linear_sweeps: self.controls.maximum_linear_sweeps, linear_relative_tolerance: self.controls.linear_relative_tolerance,
        }, work, &mut cancelled)?.into_parts();
        let report = audit(self, &next, diffusion, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(IsothermalCompetitiveGasAdvanceCandidate { state: next.state, report })
    }
}
