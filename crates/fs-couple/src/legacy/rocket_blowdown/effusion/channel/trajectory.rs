//! Adaptive finite-channel effusion, material-wall feedback and vessel recoil.
//!
//! Every trial uses the original gas/wall/mechanics transaction. Only accepted
//! half steps enter history or the fixed whole-window budgets. Geometry and
//! accommodation moments are resolved by the original interface once, never
//! inside this controller. All finite walls participate, including inventories
//! connected only by material conduction. Internal wall heat cancels in the
//! original-to-final thermal audit.
//!
//! These are local first-order time-error estimates, not certified global or
//! constitutive bounds. The original quasistatic, stationary-channel and
//! translating point-vessel qualifications remain in force.

use super::{
    EffusionChannelCandidate, EffusionChannelError, EffusionChannelOutlet, EffusionChannelStep,
    EffusionRocketEnvironment, EffusionRocketError, EffusionRocketSolve, FiniteChannelError,
    FixedVolumeChannelWall, GasRocketState, Sum,
};
use crate::rocket_blowdown::effusion::trajectory::{
    EffusionTrajectoryError, SpeciesAxis, admit_adaptive, normalized_error, ratio, trial_times,
};
use fs_material::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve,
    parcel::transfer::{
        snapshot::MAX_SNAPSHOT_PARCELS,
        species::network::{
            SpeciesNetworkError,
            transport::{
                GasTransportError,
                effusion::channel::thermal::conduction::ChannelWallConduction,
            },
        },
    },
};
use fs_math::dd::Dd;
use fs_mbd::Vec3;
use fs_qty::Temperature;
use fs_thermochem::{ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use std::collections::BTreeMap;

// Share dimensional tolerances and fixed budget semantics with the existing
// thermal/recoil owner; the multiwall driver supplies the different states.
pub use crate::rocket_blowdown::effusion::thermal::trajectory::{
    EffusionWallAdaptiveSolve as EffusionChannelAdaptiveSolve,
    EffusionWallWindowBalance as EffusionChannelWindowBalance,
};

#[derive(Debug)]
pub enum EffusionChannelTrajectoryError {
    Control(EffusionTrajectoryError),
    Step {
        elapsed_s: f64,
        duration_s: f64,
        error: EffusionChannelError,
    },
    Caloric(ThermochemErrorV1),
    Accounting(EffusionRocketError),
    WindowBalance {
        residual: EffusionChannelWindowBalance,
        debit: EffusionChannelWindowBalance,
        maximum_species_mass_residual_kg: f64,
    },
    Cancelled,
}
impl core::fmt::Display for EffusionChannelTrajectoryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive finite-channel recoil: {self:?}")
    }
}
impl std::error::Error for EffusionChannelTrajectoryError {}
impl From<EffusionTrajectoryError> for EffusionChannelTrajectoryError {
    fn from(error: EffusionTrajectoryError) -> Self {
        match error {
            EffusionTrajectoryError::Cancelled => Self::Cancelled,
            error => Self::Control(error),
        }
    }
}
impl From<EffusionRocketError> for EffusionChannelTrajectoryError {
    fn from(error: EffusionRocketError) -> Self {
        match error {
            EffusionRocketError::Cancelled => Self::Cancelled,
            error => Self::Accounting(error),
        }
    }
}
type Error = EffusionChannelTrajectoryError;

/// Accepted original transactions and accounts over exactly one physical window.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionChannelTrajectoryStep {
    pub duration_s: f64,
    pub attempted_intervals: usize,
    pub rejected_intervals: usize,
    pub accepted_intervals: usize,
    /// Two native half-step reports per accepted interval, in physical order.
    pub accepted_steps: Vec<EffusionChannelStep>,
    pub maximum_accepted_error: f64,
    pub escaped_mass_kg: f64,
    pub escaped_momentum_world_kg_m_s: Vec3,
    pub escaped_energy_world_j: f64,
    /// Same inertial frame, with each original species' u(T_window_start) zeroed
    /// BEFORE multiplication by accepted mass. Wall accommodation is included.
    pub escaped_reference_energy_world_j: f64,
    pub external_impulse_n_s: Vec3,
    pub external_work_j: f64,
    /// Actual original-to-final storage of every finite wall, in wall-ID order.
    pub wall_storage_changes_j: BTreeMap<u64, f64>,
    pub wall_storage_change_j: f64,
    /// Original-to-final gas + wall storage plus escaped relative energy. This
    /// removes species formation datums and cancels internal wall conduction.
    pub thermal_energy_residual_j: f64,
    pub residual: EffusionChannelWindowBalance,
    /// Absolute physical defects and rounding reserves. Neither signed
    /// cancellation nor the accepted step count increases the fixed allowance.
    pub debit: EffusionChannelWindowBalance,
    pub maximum_species_mass_residual_kg: f64,
}

/// Refusal never publishes a partly advanced gas, vessel or wall inventory.
#[derive(Clone, Debug)]
pub struct EffusionChannelTrajectoryCandidate<'l, 'p> {
    endpoint: GasRocketState,
    walls: BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    step: EffusionChannelTrajectoryStep,
}
impl<'l, 'p> EffusionChannelTrajectoryCandidate<'l, 'p> {
    pub const fn endpoint(&self) -> &GasRocketState {
        &self.endpoint
    }
    pub fn walls(&self) -> &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>> {
        &self.walls
    }
    pub const fn step(&self) -> &EffusionChannelTrajectoryStep {
        &self.step
    }
    pub fn into_parts(self) -> (GasRocketState, BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>, EffusionChannelTrajectoryStep) {
        (self.endpoint, self.walls, self.step)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}
fn components(v: Vec3) -> [f64; 3] { [v.x, v.y, v.z] }
fn vector(values: [f64; 3]) -> Vec3 { Vec3::new(values[0], values[1], values[2]) }
fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64, Error> {
    let mut sum = Sum::default();
    for value in values { sum.add(value)?; }
    Ok(sum.total())
}
fn finite(value: Dd) -> Result<f64, Error> {
    let value = value.to_f64();
    if value.is_finite() { Ok(value) } else {
        Err(EffusionTrajectoryError::Unrepresentable("finite-channel window account").into())
    }
}

// These are exactly the reachable native finite-wall temporal retry cases.
// Source/geometry, unrepresentable arithmetic and physical balance errors are
// hard refusals. A diagnostic message is never used to classify a retry.
fn retryable(error: &Error) -> bool {
    matches!(error, Error::Step {
        error: EffusionChannelError::Channel(
            FiniteChannelError::WallOvershoot { .. }
                | FiniteChannelError::WallDomain { .. }
                | FiniteChannelError::Transport(GasTransportError::Exchange(
                    SpeciesNetworkError::Withdrawal { .. }
                        | SpeciesNetworkError::Parcel {
                            error: GasEnergyError::OutsideBracket { .. }
                                | GasEnergyError::NotConverged { .. }, ..
                        }
                ))
        ), ..
    })
}

struct ReferenceEnergy {
    per_outlet_j: Vec<f64>,
    relative_j: f64,
    world_j: f64,
}

// The original boundary accepted these actual species masses. Their energy is
// du(T_ref,T_source) + R*T_source/2 per mole, plus the SAME accommodation heat.
// No new flux sample or subtraction of gross formation energies is required.
fn reference_energy(
    step: &EffusionChannelStep, species: &SpeciesAxis<'_>, reference_t: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ReferenceEnergy, Error> {
    let mut per_outlet_j = Vec::with_capacity(step.thermal.gas.outlets.len());
    let mut total = Sum::default();
    for outlet in &step.thermal.gas.outlets {
        poll(cancelled)?;
        let mut energy = Dd::from_f64(outlet.wall.heat_into_escaped_j);
        for row in &outlet.escaped_species {
            poll(cancelled)?;
            let index = species.sources.binary_search_by(|source| source.species().cmp(&row.species))
                .map_err(|_| EffusionTrajectoryError::Invalid("escaped channel species changed"))?;
            let source = species.sources[index];
            let du = source.molar_internal_energy_change(Temperature::new(reference_t),
                Temperature::new(outlet.source_temperature_k)).map_err(Error::Caloric)?;
            energy = energy + Dd::from_f64(row.mass_kg / source.molar_mass().value())
                * (Dd::from_f64(du.value()) + Dd::from_f64(0.5 * R) * Dd::from_f64(outlet.source_temperature_k));
        }
        let value = finite(energy)?;
        total.add(value)?;
        per_outlet_j.push(value);
    }
    let relative_j = total.total();
    let world_j = sum([relative_j, step.recoil.exhaust.relative_kinetic_transport_energy_j,
        step.recoil.source_kinetic_correction_j])?;
    Ok(ReferenceEnergy { per_outlet_j, relative_j, world_j })
}

fn comparison(coarse: f64, fine: f64, absolute: f64, controls: EffusionChannelAdaptiveSolve) -> Result<f64, Error> {
    Ok(ratio(coarse, fine, coarse.abs().max(fine.abs()), absolute, controls.motion.relative_tolerance)?)
}

fn stream_error(
    steps: [&EffusionChannelStep; 3], energy: &[ReferenceEnergy; 3],
    controls: EffusionChannelAdaptiveSolve, cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let [a, b, c] = steps;
    if [b, c].into_iter().any(|step| step.momentum.len() != a.momentum.len()
        || step.thermal.gas.outlets.len() != a.thermal.gas.outlets.len()
        || step.thermal.walls.len() != a.thermal.walls.len()
        || step.thermal.conduction.len() != a.thermal.conduction.len()) {
        return Err(EffusionTrajectoryError::Invalid("adaptive channel identities changed").into());
    }
    let mut error = 0.0_f64;
    // Preserve individual escape and local wall/source tractions: opposed
    // channels must not hide an inaccurate beam behind zero net vessel thrust.
    for ((a, b), c) in a.momentum.iter().zip(&b.momentum).zip(&c.momentum) {
        poll(cancelled)?;
        if a.channel.id != b.channel.id || a.channel.id != c.channel.id {
            return Err(EffusionTrajectoryError::Invalid("adaptive impulse IDs changed").into());
        }
        let moments = |row: &super::EffusionChannelMomentumStep| {
            [row.channel.exterior_impulses_n_s[0], row.channel.exterior_impulses_n_s[1], row.channel.wall_impulse_n_s]
        };
        for (full, (first, second)) in moments(a).into_iter().zip(moments(b).into_iter().zip(moments(c))) {
            error = error.max(comparison(full, sum([first, second])?, controls.emitted_momentum_tolerance_n_s, controls)?);
        }
    }
    for (i, ((a, b), c)) in a.thermal.gas.outlets.iter().zip(&b.thermal.gas.outlets)
        .zip(&c.thermal.gas.outlets).enumerate() {
        poll(cancelled)?;
        if a.id != b.id || a.id != c.id {
            return Err(EffusionTrajectoryError::Invalid("adaptive molecular outlet IDs changed").into());
        }
        for (full, first, second) in [
            (energy[0].per_outlet_j[i], energy[1].per_outlet_j[i], energy[2].per_outlet_j[i]),
            (a.wall.heat_into_gas_j, b.wall.heat_into_gas_j, c.wall.heat_into_gas_j),
            (a.wall.heat_into_escaped_j, b.wall.heat_into_escaped_j, c.wall.heat_into_escaped_j),
            (a.wall.heat_into_wall_j, b.wall.heat_into_wall_j, c.wall.heat_into_wall_j),
        ] {
            error = error.max(comparison(full, sum([first, second])?, controls.emitted_energy_tolerance_j, controls)?);
        }
    }
    // Molecular cooling and conduction can oppose on one wall. Each physical
    // heat path retains an independent error test before their net storage.
    for ((a, b), c) in a.thermal.walls.iter().zip(&b.thermal.walls).zip(&c.thermal.walls) {
        poll(cancelled)?;
        if a.id != b.id || a.id != c.id {
            return Err(EffusionTrajectoryError::Invalid("adaptive wall heat IDs changed").into());
        }
        for (full, first, second) in [
            (a.molecular_heat_into_wall_j, b.molecular_heat_into_wall_j, c.molecular_heat_into_wall_j),
            (a.conduction_heat_into_wall_j, b.conduction_heat_into_wall_j, c.conduction_heat_into_wall_j),
        ] {
            error = error.max(comparison(full, sum([first, second])?, controls.wall_energy_tolerance_j, controls)?);
        }
    }
    for ((a, b), c) in a.thermal.conduction.iter().zip(&b.thermal.conduction).zip(&c.thermal.conduction) {
        poll(cancelled)?;
        if a.id != b.id || a.id != c.id || a.walls != b.walls || a.walls != c.walls {
            return Err(EffusionTrajectoryError::Invalid("adaptive conduction identities changed").into());
        }
        error = error.max(comparison(a.heat_a_to_b_j, sum([b.heat_a_to_b_j, c.heat_a_to_b_j])?,
            controls.wall_energy_tolerance_j, controls)?);
    }
    for (full, (first, second)) in components(a.recoil.exhaust_momentum_world_kg_m_s).into_iter()
        .zip(components(b.recoil.exhaust_momentum_world_kg_m_s).into_iter()
            .zip(components(c.recoil.exhaust_momentum_world_kg_m_s))) {
        error = error.max(comparison(full, sum([first, second])?, controls.emitted_momentum_tolerance_n_s, controls)?);
    }
    error = error.max(comparison(energy[0].world_j, sum([energy[1].world_j, energy[2].world_j])?,
        controls.emitted_energy_tolerance_j, controls)?);
    Ok(error)
}

#[derive(Default)]
struct Totals {
    mass: Sum,
    momentum: [Dd; 3],
    energy: Sum,
    reference_relative_energy: Sum,
    reference_world_energy: Sum,
    impulse: [Dd; 3],
    work: Sum,
    mass_debit: Sum,
    momentum_defects: [Sum; 3],
    momentum_reserves: [Sum; 3],
    thermal_defects: Sum,
    mechanical_residual: Sum,
    mechanical_defects: Sum,
    energy_reserves: Sum,
}
impl Totals {
    fn add(&mut self, step: &EffusionChannelStep, energy: &ReferenceEnergy) -> Result<(), Error> {
        let recoil = &step.recoil;
        self.mass.add(recoil.expelled_mass_kg)?;
        self.energy.add(recoil.exhaust_energy_world_j)?;
        self.reference_relative_energy.add(energy.relative_j)?;
        self.reference_world_energy.add(energy.world_j)?;
        self.work.add(recoil.exhaust.external_work_j)?;
        let gas = &step.thermal.gas.audit;
        self.mass_debit.add(gas.total_mass_residual_kg.abs().max(gas.maximum_species_mass_residual_kg))?;
        self.mass_debit.add(recoil.initial_mass_residual_kg.abs())?;
        for i in 0..3 {
            self.momentum[i] = self.momentum[i] + Dd::from_f64(components(recoil.exhaust_momentum_world_kg_m_s)[i]);
            self.impulse[i] = self.impulse[i] + Dd::from_f64(components(recoil.exhaust.external_impulse_n_s)[i]);
            finite(self.momentum[i])?;
            finite(self.impulse[i])?;
            self.momentum_defects[i].add(components(recoil.momentum_residual_kg_m_s)[i].abs())?;
            self.momentum_reserves[i].add(components(recoil.initial_inventory_momentum_correction_n_s)[i].abs())?;
            self.momentum_reserves[i].add(components(recoil.gravity_impulse_uncertainty_n_s)[i])?;
        }
        self.thermal_defects.add(step.thermal.energy_residual_j.abs())?;
        for value in [recoil.exhaust.kinetic_energy_residual_j, recoil.inventory_kinetic_correction_j,
            recoil.source_kinetic_correction_j] {
            self.mechanical_residual.add(value)?;
            self.mechanical_defects.add(value.abs())?;
        }
        self.energy_reserves.add(recoil.energy_uncertainty_bound_j)?;
        Ok(())
    }
}

impl GasRocketState {
    /// Evolve channel-wall accommodation, fractionation and vessel recoil with
    /// accepted-state feedback. Geometry, areas, axes and external loads stay
    /// fixed; every trial uses the original physical owner and material laws.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_channels_adaptive<'l, 'p>(
        &self, environment: EffusionRocketEnvironment,
        walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>, outlets: &[EffusionChannelOutlet<'_>],
        duration_s: f64, endpoint: GasEnergySolve, controls: EffusionRocketSolve,
        adaptive: EffusionChannelAdaptiveSolve, cancelled: impl FnMut() -> bool,
    ) -> Result<EffusionChannelTrajectoryCandidate<'l, 'p>, Error> {
        self.advance_effusion_channels_with_conduction_adaptive(environment, walls, outlets, &[],
            duration_s, endpoint, controls, adaptive, cancelled)
    }

    /// Include physical conduction between finite channel walls in the same
    /// atomic adaptive transaction. Every wall, including a conduction-only
    /// inventory, participates in temperature/storage error control and balance.
    ///
    /// Only native candidate withdrawal/bracket/convergence or wall inversion
    /// refusals permit a smaller physical trial. Other source and physical
    /// errors refuse immediately. Export tolerances are observer-frame specific;
    /// common boosts need not choose the same adaptive time partition.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_channels_with_conduction_adaptive<'l, 'p>(
        &self, environment: EffusionRocketEnvironment,
        walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>, outlets: &[EffusionChannelOutlet<'_>],
        conduction: &[ChannelWallConduction<'_, '_>], duration_s: f64, endpoint: GasEnergySolve,
        controls: EffusionRocketSolve, adaptive: EffusionChannelAdaptiveSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EffusionChannelTrajectoryCandidate<'l, 'p>, Error> {
        advance_history(self, walls, duration_s, adaptive,
            |state, walls, dt, cancel| state.advance_effusion_channels_with_conduction(
                environment, walls, outlets, conduction, dt, endpoint, controls, cancel,
            ), &mut cancelled)
    }
}

// This controller knows time, physical states and accepted accounts. The
// callable supplies the original one-step transaction, not a second model.
pub(super) fn advance_history<'l, 'p, C: FnMut() -> bool>(
    initial: &GasRocketState, walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    duration_s: f64, adaptive: EffusionChannelAdaptiveSolve,
    mut perform: impl FnMut(&GasRocketState, &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>, f64, &mut C)
        -> Result<EffusionChannelCandidate<'l, 'p>, EffusionChannelError>,
    cancelled: &mut C,
) -> Result<EffusionChannelTrajectoryCandidate<'l, 'p>, Error> {
    poll(cancelled)?;
    admit_adaptive(duration_s, adaptive.motion)?;
    if walls.len() > MAX_SNAPSHOT_PARCELS
        || ![adaptive.wall_temperature_tolerance_k, adaptive.wall_energy_tolerance_j,
            adaptive.emitted_momentum_tolerance_n_s, adaptive.emitted_energy_tolerance_j]
            .into_iter().all(|v| v.is_finite() && v > 0.0)
        || ![adaptive.window_mass_tolerance_kg, adaptive.window_momentum_tolerance_n_s,
            adaptive.window_energy_tolerance_j].into_iter().all(|v| v.is_finite() && v >= 0.0)
    {
        return Err(EffusionTrajectoryError::Invalid(
            "bounded walls, positive wall/export error scales and nonnegative fixed window budgets required",
        ).into());
    }
    let species = SpeciesAxis::new(std::iter::once(initial.reservoir()), cancelled)?;
    let initial_masses = species.masses(initial, cancelled)?;
    let reference_t = initial.reservoir().state().state().temperature;
    let mut escaped_species: Vec<Sum> = (0..species.sources.len()).map(|_| Sum::default()).collect();
    let mut current = initial.clone();
    let mut current_walls = walls.clone();
    let (mut elapsed, mut next_step, mut maximum_error) = (0.0, adaptive.motion.initial_step_s, 0.0_f64);
    let (mut attempts, mut rejected, mut accepted) = (0, 0, 0);
    let mut accepted_steps = Vec::new();
    let mut totals = Totals::default();
    while elapsed < duration_s {
        poll(cancelled)?;
        if attempts == adaptive.motion.maximum_attempts || accepted == adaptive.motion.maximum_accepted_intervals {
            return Err(EffusionTrajectoryError::WorkBudget {
                attempted_intervals: attempts, accepted_intervals: accepted,
            }.into());
        }
        let (h, half, midpoint, end) = trial_times(elapsed, next_step, duration_s, adaptive.motion.minimum_step_s)?;
        attempts += 1;
        let trial = (|| {
            let mut step = |state: &GasRocketState, walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>, dt, at, cancel: &mut C| {
                perform(state, walls, dt, cancel).map_err(|error| match error {
                    EffusionChannelError::Channel(FiniteChannelError::Cancelled)
                    | EffusionChannelError::Channel(FiniteChannelError::Transport(GasTransportError::Cancelled))
                    | EffusionChannelError::Recoil(EffusionRocketError::Cancelled) => Error::Cancelled,
                    error => Error::Step { elapsed_s: at, duration_s: dt, error },
                }).map(|candidate| candidate.into_parts())
            };
            let (coarse, coarse_walls, coarse_step) = step(&current, &current_walls, h, elapsed, cancelled)?;
            let (first, first_walls, first_step) = step(&current, &current_walls, half, elapsed, cancelled)?;
            let (fine, fine_walls, second_step) = step(&first, &first_walls, half, midpoint, cancelled)?;
            let mut error = normalized_error(&current, &coarse, &fine, &species, h, adaptive.motion, cancelled)?;
            if !current_walls.keys().eq(coarse_walls.keys()) || !current_walls.keys().eq(fine_walls.keys()) {
                return Err(EffusionTrajectoryError::Invalid("adaptive finite-wall identities changed").into());
            }
            for (&id, old) in &current_walls {
                poll(cancelled)?;
                let (a, b) = (&coarse_walls[&id], &fine_walls[&id]);
                let stored = |next: &FixedVolumeChannelWall<'_, '_>| {
                    finite(Dd::from_f64(old.mass_kg()) * (Dd::from_f64(next.specific_energy_j_kg()) - Dd::from_f64(old.specific_energy_j_kg())))
                };
                let difference = finite(Dd::from_f64(old.mass_kg()) * (Dd::from_f64(b.specific_energy_j_kg()) - Dd::from_f64(a.specific_energy_j_kg())))?;
                error = error.max(ratio(0.0, difference, stored(a)?.abs().max(stored(b)?.abs()),
                    adaptive.wall_energy_tolerance_j, adaptive.motion.relative_tolerance)?);
                error = error.max(ratio(a.temperature_k(), b.temperature_k(),
                    old.temperature_k().max(a.temperature_k()).max(b.temperature_k()),
                    adaptive.wall_temperature_tolerance_k, adaptive.motion.relative_tolerance)?);
            }
            let energies = [reference_energy(&coarse_step, &species, reference_t, cancelled)?,
                reference_energy(&first_step, &species, reference_t, cancelled)?,
                reference_energy(&second_step, &species, reference_t, cancelled)?];
            error = error.max(stream_error([&coarse_step, &first_step, &second_step], &energies, adaptive, cancelled)?);
            let [_, first_energy, second_energy] = energies;
            Ok::<_, Error>((fine, fine_walls, [first_step, second_step], [first_energy, second_energy], error))
        })();
        let (fine, fine_walls, steps, energies, error) = match trial {
            Ok(value) => value,
            Err(error) if retryable(&error) && half >= adaptive.motion.minimum_step_s => {
                rejected += 1;
                next_step = half;
                continue;
            }
            Err(error) => return Err(error),
        };
        if error > 1.0 {
            rejected += 1;
            if half < adaptive.motion.minimum_step_s {
                return Err(EffusionTrajectoryError::MinimumStep { duration_s: h, normalized_error: error }.into());
            }
            next_step = half;
            continue;
        }
        for (step, energy) in steps.into_iter().zip(energies) {
            poll(cancelled)?;
            totals.add(&step, &energy)?;
            for outlet in &step.thermal.gas.outlets {
                for row in &outlet.escaped_species {
                    poll(cancelled)?;
                    let index = species.sources.binary_search_by(|source| source.species().cmp(&row.species))
                        .map_err(|_| EffusionTrajectoryError::Invalid("accepted escape lost its species"))?;
                    escaped_species[index].add(row.mass_kg)?;
                }
            }
            accepted_steps.push(step);
        }
        current = fine;
        current_walls = fine_walls;
        elapsed = end;
        accepted += 1;
        maximum_error = maximum_error.max(error);
        next_step = if error <= 0.25 { (2.0 * h).min(adaptive.motion.maximum_step_s) } else { h };
    }
    let final_masses = species.masses(&current, cancelled)?;
    let mut maximum_species_mass_residual_kg = 0.0_f64;
    let mut thermal_energy = Dd::from_f64(totals.reference_relative_energy.total());
    for (i, source) in species.sources.iter().enumerate() {
        poll(cancelled)?;
        maximum_species_mass_residual_kg = maximum_species_mass_residual_kg.max(
            sum([final_masses[i], -initial_masses[i], escaped_species[i].total()])?.abs());
        let du = source.molar_internal_energy_change(Temperature::new(reference_t),
            Temperature::new(current.reservoir().state().state().temperature)).map_err(Error::Caloric)?;
        thermal_energy = thermal_energy + Dd::from_f64(final_masses[i] / source.molar_mass().value()) * Dd::from_f64(du.value());
    }
    let mut wall_storage_changes_j = BTreeMap::new();
    let mut wall_storage = Dd::ZERO;
    for (&id, old) in walls {
        poll(cancelled)?;
        let change = Dd::from_f64(old.mass_kg()) * (Dd::from_f64(current_walls[&id].specific_energy_j_kg()) - Dd::from_f64(old.specific_energy_j_kg()));
        wall_storage_changes_j.insert(id, finite(change)?);
        wall_storage = wall_storage + change;
    }
    let thermal_energy_residual_j = finite(thermal_energy + wall_storage)?;
    let m0 = Dd::from_f64(initial.vehicle().dry_mass_kg + initial.reservoir().mass_kg());
    let m1 = Dd::from_f64(current.vehicle().dry_mass_kg + current.reservoir().mass_kg());
    let (mut momentum_residual, mut momentum_debit) = ([0.0; 3], [0.0; 3]);
    for i in 0..3 {
        momentum_residual[i] = finite(m1 * Dd::from_f64(components(current.vehicle().velocity_m_s)[i])
            - m0 * Dd::from_f64(components(initial.vehicle().velocity_m_s)[i])
            + totals.momentum[i] - totals.impulse[i])?;
        momentum_debit[i] = sum([totals.momentum_defects[i].total().max(momentum_residual[i].abs()),
            totals.momentum_reserves[i].total()])?;
    }
    let residual = EffusionChannelWindowBalance {
        mass_kg: sum([current.reservoir().mass_kg(), -initial.reservoir().mass_kg(), totals.mass.total()])?,
        momentum_n_s: vector(momentum_residual),
        energy_j: sum([thermal_energy_residual_j, totals.mechanical_residual.total()])?,
    };
    let debit = EffusionChannelWindowBalance {
        mass_kg: totals.mass_debit.total().max(residual.mass_kg.abs()).max(maximum_species_mass_residual_kg),
        momentum_n_s: vector(momentum_debit),
        energy_j: sum([totals.thermal_defects.total().max(thermal_energy_residual_j.abs()),
            totals.mechanical_defects.total(), totals.energy_reserves.total()])?,
    };
    if debit.mass_kg > adaptive.window_mass_tolerance_kg
        || momentum_debit.into_iter().any(|v| v > adaptive.window_momentum_tolerance_n_s)
        || debit.energy_j > adaptive.window_energy_tolerance_j
    {
        return Err(Error::WindowBalance { residual, debit, maximum_species_mass_residual_kg });
    }
    poll(cancelled)?;
    Ok(EffusionChannelTrajectoryCandidate {
        endpoint: current,
        walls: current_walls,
        step: EffusionChannelTrajectoryStep {
            duration_s: elapsed, attempted_intervals: attempts, rejected_intervals: rejected,
            accepted_intervals: accepted, accepted_steps, maximum_accepted_error: maximum_error,
            escaped_mass_kg: totals.mass.total(),
            escaped_momentum_world_kg_m_s: vector(totals.momentum.map(Dd::to_f64)),
            escaped_energy_world_j: totals.energy.total(),
            escaped_reference_energy_world_j: totals.reference_world_energy.total(),
            external_impulse_n_s: vector(totals.impulse.map(Dd::to_f64)),
            external_work_j: totals.work.total(), wall_storage_changes_j,
            wall_storage_change_j: finite(wall_storage)?, thermal_energy_residual_j,
            residual, debit, maximum_species_mass_residual_kg,
        },
    })
}
