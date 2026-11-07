//! Coupled ground-state ionization of multiple atomic charge-state ladders.
//!
//! Every admitted edge is an original `A^q + e <-> A^(q+1) + 2e`
//! impact/three-body source and its reciprocal photoionization/radiative
//! capture source. All edges share one population vector and one Maxwellian
//! electron temperature. The heavy temperature and volume remain fixed.
//! Spectators retain their original populations. No independent-electron
//! approximation, LTE projection, population floor or artificial seed is used.
//!
//! A canonical symmetric composition couples exponential-midpoint edge
//! updates. Each edge refreshes coefficients from the current shared electron
//! pool. One full/two-half comparison controls the complete network's local
//! second-order estimate, including all populations and reaction/radiation
//! accounts. Source bindings and work buffers prepare once. There are no
//! nested public solves or per-trial species allocations. This is a split
//! second-order integrator, not an unconditionally stiff or global-error claim.
//!
//! The supplied photon field is an externally maintained local intensity,
//! common to all absorbers. It is not a finite photon inventory. Emitted
//! recombination photons escape; photoelectron energy thermalizes immediately.
//! Excited levels, molecular reactions, charge exchange, continuum lowering,
//! diffuse radiation, heavy-pool evolution and electrical work are not inferred.

use super::collisional_radiative::{CollisionalRadiativeError, CollisionalRadiativeSolve};
use fs_math::det::{exp, ln};
use fs_qty::{Temperature, Time};
use fs_thermochem::electron_collision::AVOGADRO_PER_MOL as NA;
use fs_thermochem::equilibrium::composition::{
    FrozenIonizedGasStateV1, MAX_EQUILIBRIUM_SPECIES_V1,
};
use fs_thermochem::kinetics::atomic_ionization::photoionization::{
    MAX_PHOTON_GROUPS_V1, PhotonNumberFluxV1,
};
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1;

mod solve;
use solve::{Point, Problem};

/// The existing atomic source, work, state and energy refusal vocabulary.
pub type IonizationNetworkError = CollisionalRadiativeError;

/// Common temporal/source controls and explicit network work bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IonizationNetworkSolve {
    /// Shared temperature, population, radiation-energy and step allowances.
    pub local: CollisionalRadiativeSolve,
    /// Maximum admitted adjacent-stage edges, at most the original species cap.
    pub maximum_channels: usize,
    /// Maximum edge rate queries across accepted and rejected trial work.
    pub maximum_rate_evaluations: usize,
}
impl IonizationNetworkSolve {
    /// Preserve the existing atomic solver's source and numerical choices.
    pub fn new(local: CollisionalRadiativeSolve) -> Self {
        Self {
            local,
            maximum_channels: MAX_EQUILIBRIUM_SPECIES_V1,
            maximum_rate_evaluations: 10_000_000,
        }
    }
    fn validate(self) -> Result<(), IonizationNetworkError> {
        self.local.validate()?;
        if !(1..=MAX_EQUILIBRIUM_SPECIES_V1).contains(&self.maximum_channels)
            || !(1..=100_000_000).contains(&self.maximum_rate_evaluations)
        {
            return Err(IonizationNetworkError::Invalid(
                "network extent/work bounds",
            ));
        }
        Ok(())
    }
}

/// Integrated events for one canonical source edge, in moles and joules.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IonizationChannelBalance {
    /// Signed net creation of the upper stage and one electron [mol].
    pub net_extent_mol: f64,
    /// Signed impact ionization minus three-body capture [mol].
    pub collisional_extent_mol: f64,
    /// Gross photoionization events [mol].
    pub photoionization_extent_mol: f64,
    /// Gross radiative capture events [mol].
    pub radiative_recombination_extent_mol: f64,
    /// Original photon energy absorbed by this edge [J].
    pub absorbed_radiation_energy_j: f64,
    /// Original gap plus captured-electron kinetic energy radiated [J].
    pub emitted_radiation_energy_j: f64,
}

/// Independent endpoint energy and stoichiometry checks plus numerical work.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IonizationNetworkBalance {
    /// Gross absorbed photon number, summed over original source edges.
    pub absorbed_photons: f64,
    /// Gross emitted photon number; all emitted radiation escapes this model.
    pub emitted_photons: f64,
    /// Energy supplied by the maintained incident field [J].
    pub absorbed_radiation_energy_j: f64,
    /// Energy carried away by spontaneous recombination photons [J].
    pub emitted_radiation_energy_j: f64,
    /// Original level storage change reconstructed from actual populations [J].
    pub chemical_storage_change_j: f64,
    /// Actual shared electron sensible-energy change [J].
    pub electron_sensible_energy_change_j: f64,
    /// Heavy translational roundoff residual at the fixed heavy temperature [J].
    pub heavy_sensible_energy_change_j: f64,
    /// Actual gas change plus emitted minus absorbed radiation [J].
    pub energy_residual_j: f64,
    /// Largest actual species change minus integrated stoichiometric change [mol].
    pub maximum_stoichiometric_residual_mol: f64,
    /// Largest connected-ladder nucleus inventory residual [mol].
    pub maximum_nucleus_residual_mol: f64,
    /// Charge residual of the published shared state, relative to actual carriers.
    pub relative_charge_residual: f64,
    /// Largest accepted fine-step local error estimate divided by allowance.
    pub maximum_normalized_local_error: f64,
    /// Largest original radiative moment quadrature estimate encountered.
    pub maximum_relative_quadrature_error: f64,
    /// Largest conditional omitted radiative moment tail bound encountered.
    pub maximum_relative_tail_bound: f64,
    /// Accepted complete network steps.
    pub accepted_steps: usize,
    /// Rejected complete network attempts.
    pub rejected_steps: usize,
    /// Actual edge rate queries, including rejected work.
    pub rate_evaluations: usize,
}

/// Complete source-retained reacting endpoint on the original thermodynamic bank.
#[derive(Clone, Debug, PartialEq)]
pub struct IonizationNetworkResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    sources: Vec<AtomicRadiativeRecombinationV1>,
    photons: Vec<PhotonNumberFluxV1>,
    channels: Vec<IonizationChannelBalance>,
    balance: IonizationNetworkBalance,
    controls: IonizationNetworkSolve,
    elapsed: Time,
}
impl IonizationNetworkResult {
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    /// Accepted shared populations and electron temperature.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Original sources in canonical lower/upper species-index order.
    pub fn sources(&self) -> &[AtomicRadiativeRecombinationV1] {
        &self.sources
    }
    /// Per-edge accounts in the same canonical order as `sources()`.
    pub fn channel_balances(&self) -> &[IonizationChannelBalance] {
        &self.channels
    }
    pub fn photon_field(&self) -> &[PhotonNumberFluxV1] {
        &self.photons
    }
    pub const fn balance(&self) -> IonizationNetworkBalance {
        self.balance
    }
    pub const fn controls(&self) -> IonizationNetworkSolve {
        self.controls
    }
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
}

/// Evolve all admitted atomic edges through one shared electron pool.
///
/// Sources must bind exactly the initial thermodynamic bank and structural
/// electron. An edge must advance a nonnegative charge by one. Duplicate
/// reaction identities/edges, multiple predecessors/successors of one stage,
/// and cycles refuse. Disconnected ladders and original spectator species are
/// allowed; an omitted edge is not guessed. Input ordering does not change the
/// canonical numerical trajectory. An empty network is the identity.
pub fn advance_ionization_network(
    initial: &FrozenIonizedGasStateV1,
    sources: &[AtomicRadiativeRecombinationV1],
    photons: &[PhotonNumberFluxV1],
    duration: Time,
    controls: IonizationNetworkSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<IonizationNetworkResult, IonizationNetworkError> {
    controls.validate()?;
    poll(&mut cancelled)?;
    if !nonnegative(duration.value()) || sources.len() > controls.maximum_channels {
        return Err(IonizationNetworkError::Invalid(
            "duration or network channel count",
        ));
    }
    if photons.len() > MAX_PHOTON_GROUPS_V1 {
        return Err(IonizationNetworkError::Invalid("photon group count"));
    }
    for group in photons {
        poll(&mut cancelled)?;
        if !positive(group.photon_energy.value()) || !nonnegative(group.number_flux_per_m2_s) {
            return Err(IonizationNetworkError::Invalid(
                "local photon energy or intensity",
            ));
        }
    }
    let problem = Problem::new(initial, sources, photons, controls, &mut cancelled)?;
    let mut point = problem.initial_point();
    let mut full = Point::new(point.amounts.len(), point.events.len());
    let mut fine = Point::new(point.amounts.len(), point.events.len());
    let mut balance = IonizationNetworkBalance::default();
    let mut time = 0.0;
    let mut step = controls.local.initial_step_s;
    while time < duration.value() && !problem.inactive(&point) {
        poll(&mut cancelled)?;
        if balance.accepted_steps + balance.rejected_steps >= controls.local.maximum_steps {
            return Err(IonizationNetworkError::StepBudget);
        }
        let dt = step.min(duration.value() - time);
        if time + dt <= time || 0.5 * dt == 0.0 {
            return Err(IonizationNetworkError::StepBudget);
        }
        full.copy_from(&point);
        fine.copy_from(&point);
        let trial = (|| {
            problem.compose(&mut full, dt, &mut balance, &mut cancelled)?;
            problem.compose(&mut fine, 0.5 * dt, &mut balance, &mut cancelled)?;
            problem.compose(&mut fine, 0.5 * dt, &mut balance, &mut cancelled)?;
            problem.error(&point, &full, &fine)
        })();
        match trial {
            Ok(error) if error <= 1.0 => {
                core::mem::swap(&mut point, &mut fine);
                time += dt;
                balance.accepted_steps += 1;
                balance.maximum_normalized_local_error =
                    balance.maximum_normalized_local_error.max(error);
                let factor = if error == 0.0 {
                    2.0
                } else {
                    (0.85 * exp(-ln(error) / 3.0)).clamp(0.25, 2.0)
                };
                step = (dt * factor)
                    .min(controls.local.maximum_step_s)
                    .max(controls.local.minimum_step_s);
            }
            Ok(error) => {
                balance.rejected_steps += 1;
                step = dt * (0.8 * exp(-ln(error) / 3.0)).clamp(0.1, 0.5);
                if step < controls.local.minimum_step_s {
                    return Err(IonizationNetworkError::StepBudget);
                }
            }
            Err(IonizationNetworkError::OutsideDomain) => {
                balance.rejected_steps += 1;
                step = 0.5 * dt;
                if step < controls.local.minimum_step_s {
                    return Err(IonizationNetworkError::OutsideDomain);
                }
            }
            Err(error) => return Err(error),
        }
    }
    let (state, channels) = problem.publish(&point, &mut balance, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(IonizationNetworkResult {
        initial: initial.clone(),
        state,
        sources: problem.sources().map(Clone::clone).collect(),
        photons: photons.to_vec(),
        channels,
        balance,
        controls,
        elapsed: duration,
    })
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), IonizationNetworkError> {
    if cancelled() {
        Err(IonizationNetworkError::Cancelled)
    } else {
        Ok(())
    }
}
