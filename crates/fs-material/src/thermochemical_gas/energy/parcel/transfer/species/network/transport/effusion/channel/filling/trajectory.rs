//! Adaptive filling from exact vacuum through the original finite wall assembly.
//!
//! The startup estimator compares a full capture with a half capture followed
//! immediately by reciprocal flow. Only that refined pair is accepted. The
//! existing finite-wall history driver owns all subsequent intervals. Original
//! source bounds remain hard refusals; an initial pressure ceiling limits only
//! the first trial and never supplies a vacuum temperature or pressure floor.
//!
//! Gas state, source-relative chamber energy, nonlinear wall storage, gross
//! transmitted/returned species and each signed accommodation heat are compared.
//! The complete donor/receiver/all-wall window spends one fixed physical budget.
//! A closed passage can retain exact vacuum while the original walls conduct.
//! These are local temporal estimates for the parent's quasistatic rigid-volume
//! model, not a resolution of the initial non-Maxwellian filling transient.

use super::{
    ClosedGasParcel, Dd, EffusionReceiver, EvacuatedGasVolume, FiniteChannelError,
    FiniteChannelFillBalance, FiniteChannelFillAssemblyCandidate, FiniteChannelFillError,
    FiniteChannelFillLink, FiniteChannelFillStep, FiniteChannelFillTransfer,
    FiniteChannelFirstFill, FiniteChannelNetworkLink, FiniteChannelNetworkStep,
    FixedVolumeChannelWall, FiniteChannelWallStep, ChannelWallConduction,
    GasEnergySolve, GasTransportPaths, GasTransportSolve, R, SpeciesId,
    audit_assembly, audit_wall_history, positive_product_ratio, step_candidate_with_conduction,
    thermal, vacuum, FiniteVacuumChannelLink,
};
use super::super::thermal::trajectory::{self as native, FiniteChannelAdvanceError};
use super::super::super::trajectory::MAX_EFFUSION_RETAINED_ROWS;
pub use super::super::thermal::trajectory::EffusionAdvanceControls;
use fs_qty::Temperature;
use super::super::thermal::conduction::ChannelWallConductionStep;
use super::super::thermal::vacuum::trajectory::FiniteVacuumChannelAdvanceStep;
use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq)]
pub enum FiniteChannelFillAdvanceError {
    Invalid(&'static str),
    InitialPressureLimit {
        maximum_step_s: f64,
        minimum_step_s: f64,
    },
    Fill {
        duration_s: f64,
        error: FiniteChannelFillError,
    },
    /// Continuation times and work counts refer to the COMPLETE physical window.
    Advance(FiniteChannelAdvanceError),
    Cancelled,
}
impl core::fmt::Display for FiniteChannelFillAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive finite thermal channel filling: {self:?}")
    }
}
impl std::error::Error for FiniteChannelFillAdvanceError {}
type Error = FiniteChannelFillAdvanceError;

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteChannelFillAdvanceStep {
    pub duration_s: f64,
    /// Every original candidate evaluation, including rejected/coarse probes.
    pub attempted_steps: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_error: f64,
    /// First accepted half capture, or one unchanged full window if inactive.
    pub first_fill: FiniteChannelFillStep,
    pub first_reciprocal: Option<FiniteChannelNetworkStep>,
    /// Chronological original half steps after the accepted startup pair.
    pub continuation: Vec<FiniteChannelNetworkStep>,
    pub balance: FiniteChannelFillBalance,
}

#[derive(Clone, Debug)]
pub struct FiniteChannelFillAdvanceCandidate<'l, 'p> {
    donor: ClosedGasParcel,
    receiver: EffusionReceiver,
    wall: FixedVolumeChannelWall<'l, 'p>,
    step: FiniteChannelFillAdvanceStep,
}
impl<'l, 'p> FiniteChannelFillAdvanceCandidate<'l, 'p> {
    pub const fn donor(&self) -> &ClosedGasParcel {
        &self.donor
    }
    pub const fn receiver(&self) -> &EffusionReceiver {
        &self.receiver
    }
    pub const fn wall(&self) -> &FixedVolumeChannelWall<'l, 'p> {
        &self.wall
    }
    pub const fn step(&self) -> &FiniteChannelFillAdvanceStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        ClosedGasParcel,
        EffusionReceiver,
        FixedVolumeChannelWall<'l, 'p>,
        FiniteChannelFillAdvanceStep,
    ) {
        (self.donor, self.receiver, self.wall, self.step)
    }
}

/// Native accepted history distinguishes actual capture from an exactly closed
/// passage. The latter has no invented first molecule or receiver temperature.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteChannelFillAssemblyAdvanceStep {
    Filling(FiniteChannelFillAdvanceStep),
    ClosedPath {
        thermal: FiniteVacuumChannelAdvanceStep,
        balance: FiniteChannelFillBalance,
    },
}
impl FiniteChannelFillAssemblyAdvanceStep {
    pub const fn balance(&self) -> &FiniteChannelFillBalance {
        match self {
            Self::Filling(step) => &step.balance,
            Self::ClosedPath { balance, .. } => balance,
        }
    }
    pub const fn duration_s(&self) -> f64 {
        match self {
            Self::Filling(step) => step.duration_s,
            Self::ClosedPath { thermal, .. } => thermal.duration_s,
        }
    }
}

/// Donor, exact or filled receiver, and the complete original wall assembly.
#[derive(Clone, Debug)]
pub struct FiniteChannelFillAssemblyAdvanceCandidate<'l, 'p> {
    donor: ClosedGasParcel,
    receiver: EffusionReceiver,
    walls: BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    step: FiniteChannelFillAssemblyAdvanceStep,
}
impl<'l, 'p> FiniteChannelFillAssemblyAdvanceCandidate<'l, 'p> {
    pub const fn donor(&self) -> &ClosedGasParcel { &self.donor }
    pub const fn receiver(&self) -> &EffusionReceiver { &self.receiver }
    pub const fn walls(&self) -> &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>> { &self.walls }
    pub const fn step(&self) -> &FiniteChannelFillAssemblyAdvanceStep { &self.step }
    pub fn into_parts(self) -> (
        ClosedGasParcel, EffusionReceiver,
        BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>, FiniteChannelFillAssemblyAdvanceStep,
    ) {
        (self.donor, self.receiver, self.walls, self.step)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}
fn advance(error: FiniteChannelAdvanceError) -> Error {
    if error == FiniteChannelAdvanceError::Cancelled {
        Error::Cancelled
    } else {
        Error::Advance(error)
    }
}
fn fill(error: FiniteChannelFillError, duration_s: f64) -> Error {
    if error == FiniteChannelFillError::Cancelled {
        Error::Cancelled
    } else {
        Error::Fill { duration_s, error }
    }
}
fn physical(error: FiniteChannelError, elapsed_s: f64, duration_s: f64) -> Error {
    advance(if error == FiniteChannelError::Cancelled {
        FiniteChannelAdvanceError::Cancelled
    } else {
        FiniteChannelAdvanceError::Step { elapsed_s, duration_s, error }
    })
}
fn budget(elapsed_s: f64, attempted_steps: usize, accepted_intervals: usize) -> Error {
    advance(FiniteChannelAdvanceError::WorkBudget {
        elapsed_s, attempted_steps, accepted_intervals,
    })
}
fn count(attempted: &mut usize, controls: EffusionAdvanceControls) -> Result<(), Error> {
    if *attempted == controls.maximum_attempts {
        return Err(budget(0.0, *attempted, 0));
    }
    *attempted += 1;
    Ok(())
}
fn retryable(error: &Error) -> bool {
    match error {
        Error::Fill {
            error: FiniteChannelFillError::EntropyDecrease { .. }, ..
        } => true,
        Error::Fill {
            error: FiniteChannelFillError::Physical(error), duration_s,
        } => native::retryable(&FiniteChannelAdvanceError::Step {
            elapsed_s: 0.0, duration_s: *duration_s, error: error.clone(),
        }),
        Error::Advance(error) => native::retryable(error),
        // In particular, the new receiver's caloric/source errors do not retry.
        _ => false,
    }
}
fn continuation_error(error: FiniteChannelAdvanceError, offset: f64, attempted: usize) -> Error {
    advance(match error {
        FiniteChannelAdvanceError::Step { elapsed_s, duration_s, error } =>
            FiniteChannelAdvanceError::Step {
                elapsed_s: offset + elapsed_s, duration_s, error,
            },
        FiniteChannelAdvanceError::MinimumStep { elapsed_s, duration_s, normalized_error } =>
            FiniteChannelAdvanceError::MinimumStep {
                elapsed_s: offset + elapsed_s, duration_s, normalized_error,
            },
        FiniteChannelAdvanceError::WorkBudget {
            elapsed_s, attempted_steps, accepted_intervals,
        } => FiniteChannelAdvanceError::WorkBudget {
            elapsed_s: offset + elapsed_s,
            attempted_steps: attempted + attempted_steps,
            accepted_intervals: 1 + accepted_intervals,
        },
        other => other,
    })
}

// Conservative startup-only ceiling p <= n_dot*R*T_max*dt/V. The transmitted
// rate comes from the original channel and mouth, not an invented seed state.
fn pressure_limit(
    donor: &ClosedGasParcel,
    receiver: &EvacuatedGasVolume,
    wall: &FixedVolumeChannelWall<'_, '_>,
    link: FiniteChannelFillLink<'_>,
    receiver_solve: GasEnergySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<f64>, Error> {
    receiver_solve.validate()
        .map_err(|e| fill(FiniteChannelFillError::Receiver(e), 0.0))?;
    let source = |e| fill(super::endpoint(link.donor, e), 0.0);
    let area = link.channel.transmitted_area_m2(link.area_m2).map_err(source)?;
    if area > 0.0 {
        link.channel.admit_wall_for_parcel(donor, wall.temperature_k(), &mut *cancelled)
            .map_err(source)?;
    }
    let beam = donor.effusive_flux_to_vacuum(
        link.channel.mouth(), area, &mut *cancelled,
    ).map_err(source)?;
    if beam.molar_rate_mol_s == 0.0 {
        return Ok(None);
    }
    let pressure = donor.state().state().pressure
        .min(link.channel.mouth().support().pressure_pa.1);
    let ceiling = positive_product_ratio(
        &[pressure, receiver.volume_m3()],
        &[beam.molar_rate_mol_s, R, receiver_solve.maximum_temperature_k],
    ).map_err(source)?;
    let ceiling = (ceiling * (1.0 - 512.0 * f64::EPSILON)).next_down();
    if !ceiling.is_finite() || ceiling <= 0.0 {
        return Err(Error::Invalid(
            "initial thermal receiver pressure ceiling is unrepresentable",
        ));
    }
    Ok(Some(ceiling))
}

fn pair(
    candidate: &FiniteChannelFillAssemblyCandidate<'_, '_>,
    link: FiniteChannelFillLink<'_>,
) -> Result<BTreeMap<u64, ClosedGasParcel>, Error> {
    let receiver = candidate.receiver().gas()
        .ok_or(Error::Invalid("a positive capture lost its receiver"))?;
    Ok(BTreeMap::from([
        (link.donor, candidate.donor().clone()),
        (link.receiver, receiver.clone()),
    ]))
}

// The same original donor temperature anchors every trial. Species formation
// constants cancel before arithmetic; a vacuum never supplies this coordinate.
fn energies(
    original: &ClosedGasParcel,
    gases: &BTreeMap<u64, ClosedGasParcel>,
    link: FiniteChannelFillLink<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 2], FiniteChannelError> {
    let sources = original.model().thermodynamics().components();
    let mut values = [0.0; 2];
    for (value, id) in values.iter_mut().zip([link.donor, link.receiver]) {
        let gas = gases.get(&id)
            .ok_or(FiniteChannelError::Invalid("filling error estimate lost a chamber"))?;
        let mut sum = Dd::ZERO;
        for (species, mass) in gas.species_inventory_kg()
            .map_err(|_| FiniteChannelError::Invalid("unrepresentable filling species inventory"))?
        {
            if cancelled() { return Err(FiniteChannelError::Cancelled); }
            if mass == 0.0 { continue; }
            let index = sources.binary_search_by(|source| source.species().cmp(&species))
                .map_err(|_| FiniteChannelError::Invalid("filling error estimate lost an original species"))?;
            let source = &sources[index];
            let change = source.molar_internal_energy_change(
                Temperature::new(original.state().state().temperature),
                Temperature::new(gas.state().state().temperature),
            ).map_err(|_| FiniteChannelError::Invalid("filling error left original caloric support"))?
                .value();
            sum = sum + Dd::from_f64(mass) * Dd::from_f64(change)
                / Dd::from_f64(source.molar_mass().value());
        }
        *value = sum.to_f64();
        if !value.is_finite() {
            return Err(FiniteChannelError::Invalid("filling chamber energy is unrepresentable"));
        }
    }
    Ok(values)
}
fn energy_error(
    a: [f64; 2], b: [f64; 2], controls: EffusionAdvanceControls,
) -> Result<f64, FiniteChannelAdvanceError> {
    let mut error = 0.0_f64;
    for (x, y) in a.into_iter().zip(b) {
        error = error.max(native::ratio(
            x, y, x.abs().max(y.abs()), controls.stream_energy_tolerance_j,
            controls.relative_tolerance,
        )?);
    }
    Ok(error)
}

#[derive(Clone, Copy)]
enum Stage<'a> {
    Initial(&'a FiniteChannelFirstFill),
    Exchange(&'a FiniteChannelNetworkStep),
}
impl<'a> Stage<'a> {
    fn heat(self) -> (&'a [FiniteChannelWallStep], &'a [ChannelWallConductionStep]) {
        match self {
            Self::Initial(first) => (&first.transport.walls, &first.transport.conduction),
            Self::Exchange(step) => (&step.walls, &step.conduction),
        }
    }
}
fn stage(step: &FiniteChannelFillStep) -> Stage<'_> {
    match &step.transfer {
        FiniteChannelFillTransfer::Initial(first) => Stage::Initial(first),
        FiniteChannelFillTransfer::Exchange(exchange) => Stage::Exchange(exchange),
    }
}
#[derive(Default)]
struct Streams {
    // Transmitted directions and returns never cancel each other's gross mass.
    masses: BTreeMap<(u64, u64, SpeciesId), Dd>,
    heat: [[Dd; 2]; 3],
}
impl Streams {
    fn mass(&mut self, from: u64, to: u64, species: &SpeciesId, value: f64) {
        let entry = self.masses.entry((from, to, species.clone())).or_default();
        *entry = *entry + Dd::from_f64(value);
    }
    fn record(
        &mut self, report: Stage<'_>, link: FiniteChannelFillLink<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FiniteChannelAdvanceError> {
        let heat = match report {
            Stage::Initial(first) => {
                for (species, &mass) in &first.captured_species_kg {
                    if cancelled() { return Err(FiniteChannelAdvanceError::Cancelled); }
                    self.mass(link.donor, link.receiver, species, mass);
                }
                let outlet = first.transport.gas.outlets.first()
                    .ok_or(FiniteChannelAdvanceError::Invalid("startup error estimate lost its channel"))?;
                for row in &outlet.wall.reflected_species {
                    if cancelled() { return Err(FiniteChannelAdvanceError::Cancelled); }
                    self.mass(link.donor, link.donor, &row.species, row.mass_kg);
                }
                [
                    outlet.wall.heat_into_gas_j,
                    outlet.wall.heat_into_escaped_j,
                    outlet.wall.heat_into_wall_j,
                ]
            }
            Stage::Exchange(step) => {
                for row in &step.exchange.fluxes {
                    if cancelled() { return Err(FiniteChannelAdvanceError::Cancelled); }
                    let row = &row.transfer;
                    self.mass(row.donor, row.receiver, &row.species, row.mass_kg);
                }
                let channel = step.channels.first()
                    .ok_or(FiniteChannelAdvanceError::Invalid("reciprocal error estimate lost its channel"))?;
                for row in &channel.wall.reflected_species {
                    if cancelled() { return Err(FiniteChannelAdvanceError::Cancelled); }
                    for (id, mass) in channel.parcels.into_iter().zip(row.mass_kg) {
                        self.mass(id, id, &row.species, mass);
                    }
                }
                [
                    channel.wall.heat_into_gas_j[0],
                    channel.wall.heat_into_gas_j[1],
                    channel.wall.heat_into_wall_j,
                ]
            }
        };
        for (sums, value) in self.heat.iter_mut().zip(heat) {
            let side = usize::from(value < 0.0);
            sums[side] = sums[side] + Dd::from_f64(value.abs());
        }
        Ok(())
    }
}
fn stream_error(
    coarse: Stage<'_>, first: Stage<'_>, second: Stage<'_>,
    link: FiniteChannelFillLink<'_>, controls: EffusionAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, FiniteChannelAdvanceError> {
    let (mut a, mut b) = (Streams::default(), Streams::default());
    a.record(coarse, link, cancelled)?;
    b.record(first, link, cancelled)?;
    b.record(second, link, cancelled)?;
    let mut error = native::heat_error(
        coarse.heat(), first.heat(), second.heat(), controls, cancelled,
    )?;
    for key in a.masses.keys().chain(b.masses.keys()) {
        if cancelled() { return Err(FiniteChannelAdvanceError::Cancelled); }
        let x = a.masses.get(key).copied().unwrap_or(Dd::ZERO).to_f64();
        let y = b.masses.get(key).copied().unwrap_or(Dd::ZERO).to_f64();
        error = error.max(native::ratio(
            x, y, x.max(y), controls.species_mass_tolerance_kg, controls.relative_tolerance,
        )?);
    }
    for (x, y) in a.heat.into_iter().flatten().zip(b.heat.into_iter().flatten()) {
        let (x, y) = (x.to_f64(), y.to_f64());
        error = error.max(native::ratio(
            x, y, x.max(y), controls.stream_energy_tolerance_j, controls.relative_tolerance,
        )?);
    }
    Ok(error)
}
fn rows(report: Stage<'_>) -> usize {
    let (exchange, other) = match report {
        Stage::Initial(first) => (
            &first.transport.gas.exchange,
            1 + first.captured_species_kg.len() + first.transport.walls.len()
                + first.transport.conduction.len()
                + first.transport.gas.outlets.iter().map(|row| {
                    1 + row.transmitted.species.len() + row.wall.reflected_species.len()
                        + row.escaped_species.len()
                }).sum::<usize>(),
        ),
        Stage::Exchange(step) => (
            &step.exchange,
            step.walls.len() + step.conduction.len() + step.channels.iter().map(|row| {
                1 + row.transmitted.species.len() + row.wall.reflected_species.len()
            }).sum::<usize>(),
        ),
    };
    other + exchange.fluxes.len() + exchange.heat.len() + exchange.external_heat.len()
        + exchange.boundaries.len() + exchange.parcels.len()
}

struct Trial {
    step: FiniteChannelNetworkStep,
    energy: [f64; 2],
}

/// Fill the actual evacuated volume, then resample its true reciprocal beam.
///
/// One work/accuracy policy covers startup and continuation. Original gas and
/// material sources, gross withdrawals, fixed geometry and caloric support are
/// binding at every probe. All gas and finite material endpoints publish together
/// only after the original-to-final closed mass/energy/entropy audit succeeds.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate<'l, 'p, C: FnMut() -> bool>(
    donor: &ClosedGasParcel,
    receiver: &EvacuatedGasVolume,
    wall: &FixedVolumeChannelWall<'l, 'p>,
    link: FiniteChannelFillLink<'_>,
    duration_s: f64,
    donor_solve: GasEnergySolve,
    receiver_solve: GasEnergySolve,
    physical_controls: GasTransportSolve,
    controls: EffusionAdvanceControls,
    cancelled: C,
) -> Result<FiniteChannelFillAdvanceCandidate<'l, 'p>, Error> {
    let (donor, receiver, mut walls, step) = advance_candidate_with_conduction(
        donor, receiver, &BTreeMap::from([(link.wall, wall.clone())]), link, &[],
        duration_s, donor_solve, receiver_solve, physical_controls, controls, cancelled,
    )?.into_parts();
    let wall = walls.remove(&link.wall)
        .ok_or(Error::Invalid("finite filling lost its original single wall"))?;
    let FiniteChannelFillAssemblyAdvanceStep::Filling(step) = step else {
        return Err(Error::Invalid("single-wall filling produced conductive-only history"));
    };
    Ok(FiniteChannelFillAdvanceCandidate { donor, receiver, wall, step })
}

/// Fill from exact vacuum while all original finite walls conduct and evolve.
///
/// Positive capture uses the original startup comparison and continuing native
/// history driver. A closed passage instead advances its real conductive wall
/// assembly with the original vacuum driver, leaving the receiver exactly empty.
/// Both paths retain one complete closed-system mass/energy/entropy budget.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate_with_conduction<'l, 'p, C: FnMut() -> bool>(
    donor: &ClosedGasParcel,
    receiver: &EvacuatedGasVolume,
    walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
    link: FiniteChannelFillLink<'_>,
    conduction: &[ChannelWallConduction<'_, '_>],
    duration_s: f64,
    donor_solve: GasEnergySolve,
    receiver_solve: GasEnergySolve,
    physical_controls: GasTransportSolve,
    controls: EffusionAdvanceControls,
    mut cancelled: C,
) -> Result<FiniteChannelFillAssemblyAdvanceCandidate<'l, 'p>, Error> {
    poll(&mut cancelled)?;
    controls.validate()
        .map_err(|_| Error::Invalid("invalid finite thermal filling controls"))?;
    if !duration_s.is_finite() || duration_s <= 0.0 || link.donor == link.receiver {
        return Err(Error::Invalid("positive finite filling time and distinct chambers required"));
    }
    let original_walls = walls;
    let wall = walls.get(&link.wall)
        .ok_or(Error::Invalid("finite filling references an absent channel wall"))?;
    let empty = EffusionReceiver::Evacuated(receiver.clone());
    let first_step = |dt, cancelled: &mut C| {
        step_candidate_with_conduction(
            donor, &empty, walls, link, conduction, dt, donor_solve, receiver_solve,
            physical_controls, cancelled,
        ).map_err(|e| fill(e, dt))
    };
    let mut span = duration_s.min(controls.initial_step_s);
    match pressure_limit(donor, receiver, wall, link, receiver_solve, &mut cancelled)? {
        Some(ceiling) if ceiling < span => {
            if ceiling < controls.minimum_step_s {
                return Err(Error::InitialPressureLimit {
                    maximum_step_s: ceiling, minimum_step_s: controls.minimum_step_s,
                });
            }
            span = ceiling;
        }
        None => {
            if !conduction.is_empty() {
                let (mut gases, next_walls, history) =
                    vacuum::trajectory::advance_candidate_with_conduction(
                        &BTreeMap::from([(link.donor, donor.clone())]), walls,
                        GasTransportPaths::default(), &[FiniteVacuumChannelLink {
                            id: link.id, parcel: link.donor, channel: link.channel,
                            area_m2: link.area_m2, wall: link.wall,
                        }], conduction, duration_s,
                        &BTreeMap::from([(link.donor, donor_solve)]),
                        physical_controls, controls, &mut cancelled,
                    ).map_err(advance)?.into_parts();
                let next_donor = gases.remove(&link.donor)
                    .ok_or(Error::Invalid("closed filling path lost its donor"))?;
                if history.balance.gas.escaped_mass_kg != 0.0 {
                    return Err(Error::Invalid("closed filling path emitted molecular inventory"));
                }
                let balance = audit_assembly(
                    donor, &empty, walls, &next_donor, &empty, &next_walls, link,
                    physical_controls, &mut cancelled,
                ).map_err(|error| fill(error, duration_s))?;
                audit_wall_history(
                    walls, &next_walls,
                    history.accepted_steps.iter().flat_map(|step| step.walls.iter()),
                    physical_controls.exchange.energy_tolerance_j, &mut cancelled,
                ).map_err(|error| fill(error, duration_s))?;
                poll(&mut cancelled)?;
                return Ok(FiniteChannelFillAssemblyAdvanceCandidate {
                    donor: next_donor, receiver: empty, walls: next_walls,
                    step: FiniteChannelFillAssemblyAdvanceStep::ClosedPath { thermal: history, balance },
                });
            }
            let (donor, receiver, walls, first_fill) =
                first_step(duration_s, &mut cancelled)?.into_parts();
            if receiver.gas().is_some() {
                return Err(Error::Invalid("inactive channel created receiver inventory"));
            }
            let balance = first_fill.balance;
            poll(&mut cancelled)?;
            return Ok(FiniteChannelFillAssemblyAdvanceCandidate {
                donor, receiver, walls,
                step: FiniteChannelFillAssemblyAdvanceStep::Filling(FiniteChannelFillAdvanceStep {
                    duration_s,
                    attempted_steps: 1,
                    accepted_intervals: 1,
                    rejected_intervals: 0,
                    maximum_accepted_error: 0.0,
                    first_fill,
                    first_reciprocal: None,
                    continuation: Vec::new(),
                    balance,
                }),
            });
        }
        _ => {}
    }
    let links = [FiniteChannelNetworkLink {
        id: link.id, parcels: [link.donor, link.receiver], channel: link.channel,
        area_m2: link.area_m2, wall: link.wall,
    }];
    let endpoints = BTreeMap::from([(link.donor, donor_solve), (link.receiver, receiver_solve)]);
    let evaluate = |
        gases: &BTreeMap<u64, ClosedGasParcel>,
        walls: &BTreeMap<u64, FixedVolumeChannelWall<'l, 'p>>,
        dt, cancelled: &mut C,
    | -> Result<native::ThermalStepParts<'l, 'p, Trial>, FiniteChannelError> {
        let (gases, walls, step) = thermal::step_candidate_with_conduction(
            gases, walls, GasTransportPaths::default(), &links, conduction, dt, &endpoints,
            physical_controls, &mut *cancelled,
        )?.into_parts();
        let energy = energies(donor, &gases, link, cancelled)?;
        Ok((gases, walls, Trial { step, energy }))
    };
    let (mut attempted, mut rejected) = (0, 0);
    let (first_fill, mut current, mut walls, first_reciprocal, mut maximum_error) = loop {
        poll(&mut cancelled)?;
        let half = 0.5 * span;
        if half <= 0.0 || half + half != span {
            return Err(Error::Invalid("thermal filling half interval is unrepresentable"));
        }
        let trial = (|| {
            count(&mut attempted, controls)?;
            let coarse = first_step(span, &mut cancelled)?;
            let coarse_gas = pair(&coarse, link)?;
            count(&mut attempted, controls)?;
            let first = first_step(half, &mut cancelled)?;
            let first_gas = pair(&first, link)?;
            count(&mut attempted, controls)?;
            let (gas, wall, fine) = evaluate(
                &first_gas, first.walls(),
                half, &mut cancelled,
            ).map_err(|e| physical(e, half, half))?;
            let error = native::gas_error(&coarse_gas, &gas, controls, &mut cancelled)
                .map_err(advance)?
                .max(native::wall_error(
                    coarse.walls(), &wall,
                    controls, &mut cancelled,
                ).map_err(advance)?)
                .max(energy_error(
                    energies(donor, &coarse_gas, link, &mut cancelled)
                        .map_err(|e| physical(e, 0.0, span))?,
                    fine.energy, controls,
                ).map_err(advance)?)
                .max(stream_error(
                    stage(coarse.step()), stage(first.step()), Stage::Exchange(&fine.step),
                    link, controls, &mut cancelled,
                ).map_err(advance)?);
            let (_, _, _, first) = first.into_parts();
            Ok::<_, Error>((first, gas, wall, fine.step, error))
        })();
        match trial {
            Ok(result) if result.4 <= 1.0 => break result,
            Ok((_, _, _, _, error)) => {
                rejected += 1;
                let smaller = (0.5 * span).max(controls.minimum_step_s);
                if smaller >= span {
                    return Err(advance(FiniteChannelAdvanceError::MinimumStep {
                        elapsed_s: 0.0, duration_s: span, normalized_error: error,
                    }));
                }
                span = smaller;
            }
            Err(error) if retryable(&error) => {
                rejected += 1;
                let smaller = (0.5 * span).max(controls.minimum_step_s);
                if smaller >= span { return Err(error); }
                span = smaller;
            }
            Err(error) => return Err(error),
        }
    };
    let retained = rows(stage(&first_fill))
        .checked_add(rows(Stage::Exchange(&first_reciprocal)))
        .filter(|&n| n <= MAX_EFFUSION_RETAINED_ROWS)
        .ok_or_else(|| budget(span, attempted, 1))?;
    let mut accepted = 1;
    let remaining = duration_s - span;
    if remaining > 0.0 && (remaining >= duration_s || span + remaining != duration_s) {
        return Err(Error::Invalid("thermal filling continuation clock is unrepresentable"));
    }
    let continuation = if remaining > 0.0 {
        let available = controls.maximum_attempts - attempted;
        let intervals = controls.maximum_accepted_intervals - accepted;
        if available < 3 || intervals == 0 {
            return Err(budget(span, attempted, accepted));
        }
        let history = native::advance_history(
            &current, &walls, remaining,
            EffusionAdvanceControls {
                initial_step_s: span.clamp(controls.minimum_step_s, controls.maximum_step_s),
                maximum_attempts: available,
                maximum_accepted_intervals: intervals,
                ..controls
            },
            MAX_EFFUSION_RETAINED_ROWS - retained,
            evaluate,
            |coarse, first, second, cancelled| {
                Ok(energy_error(coarse.energy, second.energy, controls)?.max(
                    stream_error(
                        Stage::Exchange(&coarse.step), Stage::Exchange(&first.step),
                        Stage::Exchange(&second.step), link, controls, cancelled,
                    )?,
                ))
            },
            |trial| rows(Stage::Exchange(&trial.step)),
            &mut cancelled,
        ).map_err(|e| continuation_error(e, span, attempted))?;
        current = history.parcels;
        walls = history.walls;
        attempted += history.attempted_steps;
        accepted += history.accepted_intervals;
        rejected += history.rejected_intervals;
        maximum_error = maximum_error.max(history.maximum_accepted_error);
        history.steps.into_iter().map(|trial| trial.step).collect()
    } else {
        Vec::new()
    };
    let next_donor = current.remove(&link.donor)
        .ok_or(Error::Invalid("filled donor identity missing"))?;
    let next_receiver = EffusionReceiver::Filled(
        current.remove(&link.receiver)
            .ok_or(Error::Invalid("filled receiver identity missing"))?,
    );
    let balance = audit_assembly(
        donor, &empty, original_walls, &next_donor, &next_receiver, &walls,
        link, physical_controls, &mut cancelled,
    ).map_err(|e| fill(e, duration_s))?;
    audit_wall_history(
        original_walls, &walls,
        stage(&first_fill).heat().0.iter()
            .chain(first_reciprocal.walls.iter())
            .chain(continuation.iter().flat_map(|step: &FiniteChannelNetworkStep| step.walls.iter())),
        physical_controls.exchange.energy_tolerance_j, &mut cancelled,
    ).map_err(|e| fill(e, duration_s))?;
    poll(&mut cancelled)?;
    Ok(FiniteChannelFillAssemblyAdvanceCandidate {
        donor: next_donor, receiver: next_receiver, walls,
        step: FiniteChannelFillAssemblyAdvanceStep::Filling(FiniteChannelFillAdvanceStep {
            duration_s,
            attempted_steps: attempted,
            accepted_intervals: accepted,
            rejected_intervals: rejected,
            maximum_accepted_error: maximum_error,
            first_fill,
            first_reciprocal: Some(first_reciprocal),
            continuation,
            balance,
        }),
    })
}
