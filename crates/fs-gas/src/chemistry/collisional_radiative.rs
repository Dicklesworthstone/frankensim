//! Coupled local ground-state ionization, photoionization and recombination.
//!
//! The original atomic channel supplies electron-impact ionization and its
//! three-body reverse. Its reciprocal bound-free source supplies spontaneous
//! radiative capture and the kinetic energy carried away by captured electrons.
//! All four rates compete in the same population/energy update. This is one
//! ground-state channel, not a multilevel collisional-radiative model.
//!
//! Photon number fluxes are prescribed LOCAL, angle-integrated intensities,
//! held throughout the interval and replenished by an external radiation
//! source. This owner does not propagate, attenuate, or deplete a finite photon
//! inventory. Each absorbed photon supplies its original source energy; each
//! radiative capture emits the chemical gap plus captured kinetic energy.
//! Emitted photons escape. Photoelectrons thermalize immediately. The parcel
//! is rigid, the heavy temperature is unchanged, and no LTE projection,
//! electron/heavy heat exchange, stimulated emission or excited-state channel
//! is inferred. See the rate-equation formulation in NIST's FLYCHK manual,
//! <https://nlte.nist.gov/FLY/Doc/Manual_FLYCHK_Nov08.pdf>, section II.B.2.
//!
//! Exponential midpoint freezes all production/destruction coefficients at a
//! predicted midpoint and integrates that coupled positive two-population law
//! and its energy moments exactly. Actual electron density and temperature
//! refresh the next coefficients. Full/two-half-step comparison controls a
//! local second-order estimate; it is not a global or empirical error bound.
//! No nonlinear trial allocates a species vector. Accepted endpoint rates are
//! reused, and only a complete result publishes the original-bank gas state.

use fs_math::det::{exp, ln};
use fs_qty::{Amount, Temperature, Time};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::electron_collision::{AVOGADRO_PER_MOL as NA, MaxwellianRateControlsV1};
use fs_thermochem::equilibrium::composition::{
    FrozenIonizedGasStateErrorV1, FrozenIonizedGasStateV1,
};
use fs_thermochem::kinetics::atomic_ionization::AtomicIonizationErrorV1;
use fs_thermochem::kinetics::atomic_ionization::photoionization::{
    PhotoionizationErrorV1, PhotonNumberFluxV1,
};
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::{
    AtomicRadiativeRecombinationV1, RadiativeRecombinationErrorV1,
};
use fs_thermochem::kinetics::{KineticsErrorV1, admit_caloric_interval};

mod solve;
use solve::Problem;

/// Electron source support and separate temporal, energy and moment allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CollisionalRadiativeSolve {
    /// Supported electron temperature lower bound [K].
    pub minimum_electron_temperature_k: f64,
    /// Supported electron temperature upper bound [K].
    pub maximum_electron_temperature_k: f64,
    /// Local relative population, reaction-dose and temperature error allowance.
    pub relative_tolerance: f64,
    /// Local absolute population/reaction-dose scale [mol].
    pub absolute_amount_tolerance_mol: f64,
    /// Local absolute electron temperature scale [K].
    pub temperature_tolerance_k: f64,
    /// Original gas plus absorbed/escaped photon energy allowance [J].
    pub energy_tolerance_j: f64,
    /// Charge residual scaled by actual carriers, in (0,1).
    pub relative_charge_tolerance: f64,
    /// First proposed complete step [s].
    pub initial_step_s: f64,
    /// Smallest retry [s]; a shorter final remainder is permitted.
    pub minimum_step_s: f64,
    /// Largest proposed complete step [s].
    pub maximum_step_s: f64,
    /// Bound on accepted plus rejected attempts.
    pub maximum_steps: usize,
    /// Radiative event and captured-energy moment quadrature/tail allowances.
    pub quadrature: MaxwellianRateControlsV1,
}
impl CollisionalRadiativeSolve {
    /// Explicit source bounds; defaults are numerical, not physical accuracy.
    pub fn new(minimum: Temperature, maximum: Temperature) -> Self {
        Self {
            minimum_electron_temperature_k: minimum.value(),
            maximum_electron_temperature_k: maximum.value(),
            relative_tolerance: 1e-5,
            absolute_amount_tolerance_mol: 1e-14,
            temperature_tolerance_k: 1e-6,
            energy_tolerance_j: 1e-7,
            relative_charge_tolerance: 1e-10,
            initial_step_s: 1e-8,
            minimum_step_s: 1e-24,
            maximum_step_s: 1e-3,
            maximum_steps: 100_000,
            quadrature: MaxwellianRateControlsV1 {
                relative_tolerance: 1e-10,
                maximum_subdivisions: 4096,
            },
        }
    }
    fn validate(self) -> Result<(), CollisionalRadiativeError> {
        if [
            self.minimum_electron_temperature_k,
            self.maximum_electron_temperature_k,
            self.relative_tolerance,
            self.absolute_amount_tolerance_mol,
            self.temperature_tolerance_k,
            self.energy_tolerance_j,
            self.relative_charge_tolerance,
            self.initial_step_s,
            self.minimum_step_s,
            self.maximum_step_s,
        ]
        .iter()
        .any(|v| !positive(*v))
            || self.minimum_electron_temperature_k >= self.maximum_electron_temperature_k
            || self.relative_tolerance > 0.1
            || self.relative_charge_tolerance >= 1.0
            || self.minimum_step_s > self.initial_step_s
            || self.initial_step_s > self.maximum_step_s
            || !(1..=1_000_000).contains(&self.maximum_steps)
        {
            return Err(CollisionalRadiativeError::Invalid("integration controls"));
        }
        Ok(())
    }
}

/// Source or evolution refusal; incoming state and sources remain immutable.
#[derive(Clone, Debug, PartialEq)]
pub enum CollisionalRadiativeError {
    /// Invalid controls, duration or original-bank binding.
    Invalid(&'static str),
    /// Atomic impact or three-body detailed-balance source refused.
    Atomic(AtomicIonizationErrorV1),
    /// Original local photon field/source refused.
    Photoionization(PhotoionizationErrorV1),
    /// Reciprocal recombination moments or their qualification refused.
    Recombination(RadiativeRecombinationErrorV1),
    /// Original caloric source cannot support the declared temperature path.
    Admission(KineticsErrorV1),
    /// Positive electron temperature left the admitted source interval.
    OutsideDomain,
    /// Explicit work, minimum-step or representable-time budget exhausted.
    StepBudget,
    /// A required rate, population, dose or energy is not representable.
    Arithmetic(&'static str),
    /// Original gas energy does not close with absorbed/escaped radiation.
    EnergyBalance { residual_j: f64 },
    /// Shared final charged-state publication refused.
    State(FrozenIonizedGasStateErrorV1),
    /// Cancellation before publication of a complete result.
    Cancelled,
}
impl core::fmt::Display for CollisionalRadiativeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "local collisional-radiative evolution: {self:?}")
    }
}
impl std::error::Error for CollisionalRadiativeError {}

/// Actual population, photon and gas-energy accounts and bounded numerical work.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CollisionalRadiativeBalance {
    /// Gross photon-driven creation of ion/electron pairs [mol].
    pub photoionization_extent_mol: f64,
    /// Signed electron-impact minus three-body capture extent [mol].
    pub collisional_extent_mol: f64,
    /// Gross spontaneous radiative captures [mol].
    pub radiative_recombination_extent_mol: f64,
    /// One absorbed photon per photoionization event.
    pub absorbed_photons: f64,
    /// One escaping photon per radiative capture.
    pub emitted_photons: f64,
    /// Energy removed from the externally maintained local photon field [J].
    pub absorbed_radiation_energy_j: f64,
    /// Released chemical plus captured kinetic energy in escaping photons [J].
    pub emitted_radiation_energy_j: f64,
    /// Original chemical formation-energy change from NET pair creation [J].
    pub chemical_storage_change_j: f64,
    /// Actual electron sensible-energy change, including population dilution [J].
    pub electron_sensible_energy_change_j: f64,
    /// Gas energy change + emitted radiation - absorbed radiation [J].
    pub energy_residual_j: f64,
    /// Largest accepted estimated fine-step error divided by its allowance.
    pub maximum_normalized_local_error: f64,
    /// Largest radiative moment quadrature estimate encountered.
    pub maximum_relative_quadrature_error: f64,
    /// Largest conditional radiative moment omitted-tail bound encountered.
    pub maximum_relative_tail_bound: f64,
    /// Complete accepted steps; each retains its two-half-step endpoint.
    pub accepted_steps: usize,
    /// Rejected attempts.
    pub rejected_steps: usize,
}

/// Complete source-retained local endpoint, reusable by material/field owners.
#[derive(Clone, Debug, PartialEq)]
pub struct CollisionalRadiativeResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    source: AtomicRadiativeRecombinationV1,
    photons: Vec<PhotonNumberFluxV1>,
    controls: CollisionalRadiativeSolve,
    elapsed: Time,
    extent: Amount,
    balance: CollisionalRadiativeBalance,
}
impl CollisionalRadiativeResult {
    /// Accepted populations and separate temperatures on the original bank.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Original input state, including its own electron temperature.
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    /// Original reciprocal photon source and atomic channel.
    pub const fn source(&self) -> &AtomicRadiativeRecombinationV1 {
        &self.source
    }
    /// Prescribed local fluxes held over the completed interval.
    pub fn photon_field(&self) -> &[PhotonNumberFluxV1] {
        &self.photons
    }
    /// Applied source, quadrature and temporal controls.
    pub const fn controls(&self) -> CollisionalRadiativeSolve {
        self.controls
    }
    /// Completed physical duration [s].
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
    /// Signed NET creation of an ion/electron pair [mol].
    pub const fn reaction_extent(&self) -> Amount {
        self.extent
    }
    /// Photon, gas-energy, population and numerical-work accounts.
    pub const fn balance(&self) -> CollisionalRadiativeBalance {
        self.balance
    }
}

/// Evolve collisional and photon-driven creation/capture simultaneously.
///
/// The photon field is local and externally maintained, not a finite exposure
/// through a column. Empty `photons` means darkness. Exact absence stays absent
/// in darkness; nonzero photoionization can create the first physical electrons.
/// Spectators keep their original amounts and heavy temperature. Radiation and
/// thermalization assumptions are those stated in the module documentation.
pub fn advance_collisional_radiative(
    initial: &FrozenIonizedGasStateV1,
    source: &AtomicRadiativeRecombinationV1,
    photons: &[PhotonNumberFluxV1],
    duration: Time,
    controls: CollisionalRadiativeSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<CollisionalRadiativeResult, CollisionalRadiativeError> {
    controls.validate()?;
    poll(&mut cancelled)?;
    let end = duration.value();
    if !nonnegative(end) {
        return Err(CollisionalRadiativeError::Invalid("nonnegative duration"));
    }
    let problem = Problem::new(initial, source, photons, controls, &mut cancelled)?;
    let mut point = problem.initial_point();
    let mut balance = CollisionalRadiativeBalance::default();
    let mut rates = problem.rates(point, &mut balance, &mut cancelled)?;
    let mut time = 0.0;
    let mut step = controls.initial_step_s;
    while time < end {
        poll(&mut cancelled)?;
        if rates.inactive(point) {
            break;
        }
        if balance.accepted_steps + balance.rejected_steps >= controls.maximum_steps {
            return Err(CollisionalRadiativeError::StepBudget);
        }
        let dt = step.min(end - time);
        if time + dt <= time || 0.5 * dt == 0.0 {
            return Err(CollisionalRadiativeError::StepBudget);
        }
        let trial = (|| {
            let full = problem.midpoint(point, rates, dt, &mut balance, &mut cancelled)?;
            let half = problem.midpoint(point, rates, 0.5 * dt, &mut balance, &mut cancelled)?;
            let middle_rates = problem.rates(half, &mut balance, &mut cancelled)?;
            let next =
                problem.midpoint(half, middle_rates, 0.5 * dt, &mut balance, &mut cancelled)?;
            let error = problem.error(point, full, next)?;
            Ok((next, error))
        })();
        match trial {
            Ok((next, error)) if error <= 1.0 => {
                point = next;
                time += dt;
                balance.accepted_steps += 1;
                balance.maximum_normalized_local_error =
                    balance.maximum_normalized_local_error.max(error);
                rates = problem.rates(point, &mut balance, &mut cancelled)?;
                let factor = if error == 0.0 {
                    2.0
                } else {
                    (0.85 * exp(-ln(error) / 3.0)).clamp(0.25, 2.0)
                };
                step = (dt * factor)
                    .min(controls.maximum_step_s)
                    .max(controls.minimum_step_s);
            }
            Ok((_, error)) => {
                balance.rejected_steps += 1;
                step = dt * (0.8 * exp(-ln(error) / 3.0)).clamp(0.1, 0.5);
                if step < controls.minimum_step_s {
                    return Err(CollisionalRadiativeError::StepBudget);
                }
            }
            Err(CollisionalRadiativeError::OutsideDomain) => {
                balance.rejected_steps += 1;
                step = 0.5 * dt;
                if step < controls.minimum_step_s {
                    return Err(CollisionalRadiativeError::OutsideDomain);
                }
            }
            Err(error) => return Err(error),
        }
    }
    let state = problem.publish(point, &mut balance, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(CollisionalRadiativeResult {
        initial: initial.clone(),
        state,
        source: source.clone(),
        photons: photons.to_vec(),
        controls,
        elapsed: duration,
        extent: Amount::new(point.extent.value),
        balance,
    })
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CollisionalRadiativeError> {
    if cancelled() {
        Err(CollisionalRadiativeError::Cancelled)
    } else {
        Ok(())
    }
}
