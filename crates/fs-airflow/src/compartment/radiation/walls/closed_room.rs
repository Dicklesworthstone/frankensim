//! Implicit endpoint heat exchange between original finite gas and wall owners.
//!
//! Unknowns are the heat powers into the present gas layers. Each trial advances
//! a clone of the ORIGINAL `TwoZone` by those powers, including its common
//! pressure and interlayer pressure-volume work. A caller then evaluates the
//! original implicit wall/radiation candidate at that endpoint gas state.
//! It must rebuild temperature/pressure-dependent opacity, optical paths and
//! convective contacts from the supplied room; the driver invents none of them.
//! The residual is applied gas heat minus the candidate's returned gas heat.
//!
//! Acceptance advances gas with the EXACT heat vector in the returned wall
//! report. A second endpoint evaluation gates the remaining power and gas/wall
//! temperature mismatch. Combined energy is measured, never redistributed.
//! Failed probes and cancellation mutate only numerical scratch. No timestep
//! changes, mass transport, vents, reactions or volatile release are admitted.
//! This is a first-order endpoint heat coupling using the original gas owner's
//! pressure-work convention, not an unconditional stability/convergence theorem.

use super::super::RadiationError;
use super::{EnclosureWallCandidate, EnclosureWallError, EnclosureWallWork};
use fs_gas::multicomponent::compartment::{
    Layer, LayerRates, TwoZone, ZoneError, ZoneSolve, ZoneStepReport,
};
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};
use fs_qty::{Power, Time};

/// Original gas controls and explicit bounded nonlinear work/acceptance gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClosedRoomWallSolve {
    pub gas: ZoneSolve,
    /// Outer residual evaluations through acceptance, 1..=64.
    pub maximum_iterations: usize,
    /// Globalization trials per Newton update, 1..=32.
    pub maximum_backtracks: usize,
    /// Total gas/wall candidate trials, including derivatives and rejected probes.
    pub maximum_evaluations: usize,
    /// Finite-difference increment of one gas power [W]; never a timestep change.
    pub power_difference_w: f64,
    /// Actual trial/applied and returned-endpoint heat mismatch per layer [W].
    pub power_tolerance_w: f64,
    /// Original-query versus returned-endpoint gas and wall-face mismatch [K].
    pub temperature_tolerance_k: f64,
    /// Combined original gas/wall storage minus external wall heat [J].
    pub combined_energy_tolerance_j: f64,
}
impl Default for ClosedRoomWallSolve {
    fn default() -> Self {
        Self {
            gas: ZoneSolve::default(),
            maximum_iterations: 32,
            maximum_backtracks: 24,
            maximum_evaluations: 1024,
            power_difference_w: 0.01,
            power_tolerance_w: 1e-4,
            temperature_tolerance_k: 1e-6,
            combined_energy_tolerance_j: 1e-4,
        }
    }
}

/// Only the small original LU and active-layer axis are retained between calls.
pub struct ClosedRoomWallWorkspace {
    active: Vec<usize>,
    lu: LuWorkspace,
}
impl ClosedRoomWallWorkspace {
    pub fn new(room: &TwoZone) -> Result<Self, ClosedRoomWallError> {
        let active = active_layers(room);
        Ok(Self {
            lu: LuWorkspace::new(active.len())?,
            active,
        })
    }
}

/// Original gas step, actual endpoint mismatches and the unmodified energy sum.
#[derive(Clone, Debug, PartialEq)]
pub struct ClosedRoomWallReport {
    pub gas_step: ZoneStepReport,
    pub iterations: usize,
    pub evaluations: usize,
    /// Original wall/optical work across successful candidate trials, including
    /// discarded derivatives/backtracks. Failed inner solves return no counts.
    pub successful_trial_work: EnclosureWallWork,
    /// Fresh endpoint wall/gas heat minus the identical heat applied to both owners.
    pub endpoint_power_residual_w: [f64; 2],
    pub maximum_gas_temperature_residual_k: f64,
    pub maximum_wall_temperature_residual_k: f64,
    pub combined_energy_residual_j: f64,
}

/// Install both owners together. Neither accepted input was mutated by the call.
pub struct ClosedRoomWallCandidate<'m, 's> {
    pub room: TwoZone,
    pub walls: EnclosureWallCandidate<'m, 's>,
    pub report: ClosedRoomWallReport,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ClosedRoomWallError {
    Invalid(&'static str),
    Gas(ZoneError),
    Walls(EnclosureWallError),
    Linear(LuWorkspaceError),
    Budget { evaluations: usize, maximum: usize },
    Iterations,
    EnergyResidual { residual_j: f64, allowed_j: f64 },
    Cancelled,
}
impl core::fmt::Display for ClosedRoomWallError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "closed gas/wall endpoint radiation: {self:?}")
    }
}
impl std::error::Error for ClosedRoomWallError {}
impl From<ZoneError> for ClosedRoomWallError {
    fn from(e: ZoneError) -> Self {
        if e == ZoneError::Cancelled {
            Self::Cancelled
        } else {
            Self::Gas(e)
        }
    }
}
impl From<EnclosureWallError> for ClosedRoomWallError {
    fn from(e: EnclosureWallError) -> Self {
        if e == EnclosureWallError::Cancelled {
            Self::Cancelled
        } else {
            Self::Walls(e)
        }
    }
}
impl From<LuWorkspaceError> for ClosedRoomWallError {
    fn from(e: LuWorkspaceError) -> Self {
        Self::Linear(e)
    }
}
fn active_layers(room: &TwoZone) -> Vec<usize> {
    [Layer::Lower, Layer::Upper]
        .into_iter()
        .filter(|&layer| room.layer(layer).mass_kg() > 0.0)
        .map(Layer::index)
        .collect()
}
fn check(cancelled: &mut dyn FnMut() -> bool) -> Result<(), ClosedRoomWallError> {
    if cancelled() {
        Err(ClosedRoomWallError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn fatal(e: &ClosedRoomWallError) -> bool {
    matches!(
        e,
        ClosedRoomWallError::Cancelled
            | ClosedRoomWallError::Budget { .. }
            | ClosedRoomWallError::Walls(
                EnclosureWallError::Budget { .. }
                    | EnclosureWallError::Radiation(RadiationError::SpectralBudget { .. })
            )
    )
}

struct Trial<'m, 's> {
    room: TwoZone,
    gas_step: ZoneStepReport,
    walls: EnclosureWallCandidate<'m, 's>,
    power: [f64; 2],
    residual: [f64; 2],
    norm: f64,
}

/// Solve heat exchange for a closed, fixed-inventory original room.
///
/// `wall_candidate` must always start from the SAME accepted wall states and
/// requested duration, and use the supplied ENDPOINT room for all gas-dependent
/// inputs. Call `candidate_walls`, not `candidate_room_walls`: the latter checks
/// same-time initial owners, whereas the room here has already advanced. The
/// closure may reuse scratch but must not publish a trial's wall state.
/// Existing gray, spectral, selective-surface and soot optics remain available.
/// An optional power vector is a numerical initial guess, not an extra source.
#[allow(clippy::too_many_arguments)]
pub fn candidate_closed_room_walls<'m, 's>(
    room: &TwoZone,
    duration: Time,
    initial_gas_power: Option<[Power; 2]>,
    controls: ClosedRoomWallSolve,
    work: &mut ClosedRoomWallWorkspace,
    mut wall_candidate: impl FnMut(
        &TwoZone,
        &mut dyn FnMut() -> bool,
    ) -> Result<EnclosureWallCandidate<'m, 's>, EnclosureWallError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ClosedRoomWallCandidate<'m, 's>, ClosedRoomWallError> {
    check(&mut cancelled)?;
    if work.active != active_layers(room)
        || !(1..=64).contains(&controls.maximum_iterations)
        || !(1..=32).contains(&controls.maximum_backtracks)
        || !(1..=65_536).contains(&controls.maximum_evaluations)
        || !positive(duration.value())
        || !positive(controls.power_difference_w)
        || !positive(controls.power_tolerance_w)
        || !positive(controls.temperature_tolerance_k)
        || !positive(controls.combined_energy_tolerance_j)
    {
        return Err(ClosedRoomWallError::Invalid(
            "closed-room axes, duration or solve controls",
        ));
    }
    let initial = initial_gas_power
        .unwrap_or([Power::new(0.0); 2])
        .map(Power::value);
    let mut evaluations = 0;
    let mut successful_trial_work = EnclosureWallWork::default();
    let mut evaluate = |power: [f64; 2], cancel: &mut dyn FnMut() -> bool| {
        check(cancel)?;
        if evaluations >= controls.maximum_evaluations {
            return Err(ClosedRoomWallError::Budget {
                evaluations,
                maximum: controls.maximum_evaluations,
            });
        }
        evaluations += 1;
        if power.iter().any(|p| !p.is_finite())
            || (0..2).any(|i| !work.active.contains(&i) && power[i] != 0.0)
        {
            return Err(ClosedRoomWallError::Invalid(
                "gas heat must be finite and absent layers receive zero",
            ));
        }
        let mut next = room.clone();
        let mut rates = [LayerRates::default(); 2];
        for (rate, p) in rates.iter_mut().zip(power) {
            rate.add_heat(Power::new(p))?;
        }
        let gas_step = next.advance(duration, rates, controls.gas, &mut *cancel)?;
        let walls = wall_candidate(&next, cancel)?;
        successful_trial_work.wall_evaluations += walls.report.work.wall_evaluations;
        successful_trial_work.radiation_evaluations += walls.report.work.radiation_evaluations;
        successful_trial_work.wavelength_evaluations += walls.report.work.wavelength_evaluations;
        successful_trial_work.spectral_factorizations += walls.report.work.spectral_factorizations;
        if walls.walls.is_empty()
            || walls.report.elapsed != duration
            || walls.report.exported_volatile_enthalpy.value() != 0.0
            || walls
                .walls
                .iter()
                .any(|wall| wall.wall.time() != next.time())
        {
            return Err(ClosedRoomWallError::Invalid(
                "wall candidates must share the exact closed gas interval without volatile export",
            ));
        }
        let mut residual = [0.0; 2];
        let mut norm = 0.0_f64;
        for i in 0..2 {
            let actual = walls.report.gas_heat_power[i].value();
            residual[i] = power[i] - actual;
            if !residual[i].is_finite() || (!work.active.contains(&i) && actual != 0.0) {
                return Err(ClosedRoomWallError::Invalid(
                    "wall heat must be finite and absent layers receive zero",
                ));
            }
            norm = norm.max(residual[i].abs());
        }
        check(cancel)?;
        Ok(Trial {
            room: next,
            gas_step,
            walls,
            power,
            residual,
            norm,
        })
    };
    let mut trial = evaluate(initial, &mut cancelled)?;
    let n = work.active.len();
    for iteration in 1..=controls.maximum_iterations {
        check(&mut cancelled)?;
        if trial.norm <= controls.power_tolerance_w {
            let applied = trial.walls.report.gas_heat_power.map(Power::value);
            let endpoint = evaluate(applied, &mut cancelled)?;
            if endpoint.walls.walls.len() != trial.walls.walls.len() {
                return Err(ClosedRoomWallError::Invalid(
                    "wall candidate axis changed during endpoint feedback",
                ));
            }
            let gas_mismatch = [Layer::Lower, Layer::Upper]
                .into_iter()
                .map(|layer| {
                    (endpoint.room.layer(layer).temperature().value()
                        - trial.room.layer(layer).temperature().value())
                    .abs()
                })
                .fold(0.0_f64, f64::max);
            let wall_mismatch = trial
                .walls
                .walls
                .iter()
                .zip(&endpoint.walls.walls)
                .flat_map(|(a, b)| {
                    a.report
                        .surface_temperatures
                        .into_iter()
                        .zip(b.report.surface_temperatures)
                        .map(|(a, b)| (a.value() - b.value()).abs())
                })
                .fold(0.0_f64, f64::max);
            if endpoint.norm <= controls.power_tolerance_w
                && gas_mismatch <= controls.temperature_tolerance_k
                && wall_mismatch <= controls.temperature_tolerance_k
            {
                // Species inventories are fixed: sensible-energy change is
                // physical-energy change without subtracting formation offsets.
                let energy = endpoint.room.sensible_internal_energy().value()
                    - room.sensible_internal_energy().value()
                    + trial.walls.report.wall_storage_change.value()
                    - trial.walls.report.ambient_heat_energy.value()
                    - trial.walls.report.external_radiant_energy.value();
                if !energy.is_finite() || energy.abs() > controls.combined_energy_tolerance_j {
                    return Err(ClosedRoomWallError::EnergyResidual {
                        residual_j: energy,
                        allowed_j: controls.combined_energy_tolerance_j,
                    });
                }
                check(&mut cancelled)?;
                return Ok(ClosedRoomWallCandidate {
                    room: endpoint.room,
                    walls: trial.walls,
                    report: ClosedRoomWallReport {
                        gas_step: endpoint.gas_step,
                        iterations: iteration,
                        evaluations,
                        successful_trial_work,
                        endpoint_power_residual_w: endpoint.residual.map(|r| -r),
                        maximum_gas_temperature_residual_k: gas_mismatch,
                        maximum_wall_temperature_residual_k: wall_mismatch,
                        combined_energy_residual_j: energy,
                    },
                });
            }
            trial = endpoint;
            continue;
        }
        let mut matrix = [0.0; 4];
        for column in 0..n {
            let index = work.active[column];
            let mut probe = trial.power;
            probe[index] += controls.power_difference_w;
            let result = match evaluate(probe, &mut cancelled) {
                Ok(result) => result,
                Err(e) if fatal(&e) => return Err(e),
                Err(_) => {
                    probe[index] = trial.power[index] - controls.power_difference_w;
                    evaluate(probe, &mut cancelled)?
                }
            };
            let increment = probe[index] - trial.power[index];
            if !increment.is_finite() || increment == 0.0 {
                return Err(ClosedRoomWallError::Invalid(
                    "gas power finite difference is unresolvable",
                ));
            }
            for row in 0..n {
                let index = work.active[row];
                matrix[row * n + column] =
                    (result.residual[index] - trial.residual[index]) / increment;
            }
        }
        let mut rhs = [0.0; 2];
        let mut delta = [0.0; 2];
        for (row, &i) in work.active.iter().enumerate() {
            rhs[row] = -trial.residual[i];
        }
        work.lu.factor(&matrix[..n * n])?;
        work.lu.solve_factored_into(&rhs[..n], &mut delta[..n])?;
        let mut scale = 1.0;
        let mut accepted = None;
        for _ in 0..controls.maximum_backtracks {
            let mut powers = trial.power;
            for (column, &i) in work.active.iter().enumerate() {
                powers[i] += scale * delta[column];
            }
            match evaluate(powers, &mut cancelled) {
                Ok(candidate) if candidate.norm < trial.norm => {
                    accepted = Some(candidate);
                    break;
                }
                Err(e) if fatal(&e) => return Err(e),
                _ => scale *= 0.5,
            }
        }
        trial = accepted.ok_or(ClosedRoomWallError::Iterations)?;
    }
    Err(ClosedRoomWallError::Iterations)
}
