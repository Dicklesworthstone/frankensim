//! Finite light-ion burn and explicitly selected plasma-energy closures.
//!
//! The original D–T entrypoints below evolve a two-temperature thermal state;
//! [`network`] adds competing light-ion channels and prescribed compression.
//! [`kinetic`] retains a separate isotropic alpha distribution with finite-time
//! energy drift, diffusion, pressure and reciprocal thermal-pool heating.
//!
//! The following closure describes the original D–T thermal entrypoints.
//!
//! This is a homogeneous, rigid, quasineutral Maxwellian parcel. Its admitted
//! monatomic source bank supplies D+, T+, He4++ and electrons. Nuclear charge
//! and baryon inventories replace chemical element conservation for the burn;
//! the original chemical kinetics and equilibrium owners are not modified.
//!
//! Alpha particles are retained and **instantaneously thermalized**, with a
//! caller-declared electron/heavy energy split. Neutrons escape. These are
//! explicit closure assumptions, not a slowing-down or neutron-transport model.
//! Neutrons carry both their reaction Q share and their share of reacting-ion
//! kinetic energy. Reaction weighting follows the temperature derivative of
//! the same Maxwellian reactivity, rather than removing two average ions.
//! Nuclear energy is counted once, outside the zero-formation-energy calorics.
//!
//! Each frozen-temperature reaction substep integrates the finite D/T supply
//! analytically, including unequal inventories. An exponential midpoint method
//! with step doubling controls the changing temperature and reactivity. Trial
//! storage is fixed size; no trial allocates, clips negative fuel, inserts a
//! trace seed or projects the result to chemical equilibrium. Publication is
//! atomic on cancellation, domain refusal or exhausted work.
//!
//! The local error estimate is not a global error certificate. This owner does
//! not solve confinement, spatial transport, fast-particle distributions,
//! radiation transport or magnetic fields. Optional original-source Coulomb
//! exchange and optically thin bremsstrahlung participate in the same thermal
//! feedback and energy account. Its
//! returned common state can enter existing thermal owners in a separately
//! declared coupling schedule, subject to each owner's physical domain.

use fs_qty::{Amount, Energy, Power, Temperature, Time};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::coulomb::{CoulombErrorV1, CoulombModelV1};
use fs_thermochem::equilibrium::composition::{
    FrozenIonizedGasStateErrorV1, FrozenIonizedGasStateV1,
};
use fs_thermochem::fusion::plasma::{
    DtFusionPlasmaErrorV1, DtFusionPlasmaLimitsV1, DtFusionPlasmaV1,
};
use fs_thermochem::fusion::{FusionErrorV1, FusionReactionV1};
use fs_thermochem::plasma_radiation::{BremsstrahlungErrorV1, ElectronIonBremsstrahlungV1};

mod solve;
mod thermal;

/// Finite-time alpha slowing and diffusion with separate kinetic storage.
pub mod kinetic;
/// Competing D–T, D–D and D–He3 reactions on the same material state.
pub mod network;
pub mod volume;

const AVOGADRO: f64 = 6.022_140_76e23;
const CV_MOLAR: f64 = 1.5 * R;

/// Explicit local alpha thermalization and constant externally supplied powers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DtFusionHeating {
    /// Fraction of the complete born-alpha kinetic energy deposited in electrons.
    /// The remainder heats heavy particles. Must lie in [0,1]; no default
    /// slowing-down fraction is inferred from a material name.
    pub alpha_electron_fraction: f64,
    /// Constant external electron heating [W]; negative means prescribed cooling.
    pub electron_power: Power,
    /// Constant external heavy-particle heating [W]; negative means cooling.
    pub heavy_power: Power,
}

/// Optional thermal operators sharing the burn's original species bank.
///
/// Coulomb energy is transferred reciprocally between the two thermal pools.
/// Bremsstrahlung removes electron energy and records it as escaping photons:
/// selecting it explicitly assumes an optically thin parcel. Its prescribed
/// Gaunt factor and the Coulomb screening/cutoff remain part of these sources.
/// Every evaluation refreshes the actual fuel/ash densities and temperatures.
/// An absent operator contributes exactly zero; it is not inferred from names.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DtFusionTransport {
    /// Original reciprocal Maxwellian electron/ion thermal-exchange law.
    pub coulomb: Option<CoulombModelV1>,
    /// Original electron-ion free-free emission law, with all photons escaping.
    pub bremsstrahlung: Option<ElectronIonBremsstrahlungV1>,
}

/// Declared physical domain and bounded numerical work for one burn interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DtFusionSolve {
    /// Lower supported ion temperature [K], also inside the reactivity fit.
    pub minimum_heavy_temperature_k: f64,
    /// Upper supported ion temperature [K], also inside the reactivity fit.
    pub maximum_heavy_temperature_k: f64,
    /// Lower supported electron temperature [K].
    pub minimum_electron_temperature_k: f64,
    /// Upper supported electron temperature [K].
    pub maximum_electron_temperature_k: f64,
    /// Explicit weak-coupling, nondegeneracy and nonrelativistic domain.
    pub plasma_limits: DtFusionPlasmaLimitsV1,
    /// Estimated local relative error on populations, temperatures and energies.
    pub relative_tolerance: f64,
    /// Estimated local absolute amount error [mol].
    pub absolute_amount_tolerance_mol: f64,
    /// Estimated local absolute temperature error [K].
    pub absolute_temperature_tolerance_k: f64,
    /// Absolute numerical allowance on the complete energy account [J].
    pub energy_tolerance_j: f64,
    /// Relative numerical allowance on energy and nuclear inventory accounts.
    pub balance_relative_tolerance: f64,
    /// Initial proposed time step [s].
    pub initial_step_s: f64,
    /// Smallest retry [s]; a shorter final remainder is admitted.
    pub minimum_step_s: f64,
    /// Largest proposed time step [s].
    pub maximum_step_s: f64,
    /// Maximum total accepted and rejected attempts, at most one million.
    pub maximum_steps: usize,
}

impl DtFusionSolve {
    /// Set the physical domain explicitly and select numerical defaults.
    /// The defaults describe numerical resolution, not measured fusion accuracy.
    pub fn new(
        minimum_heavy: Temperature,
        maximum_heavy: Temperature,
        minimum_electron: Temperature,
        maximum_electron: Temperature,
        plasma_limits: DtFusionPlasmaLimitsV1,
    ) -> Self {
        Self {
            minimum_heavy_temperature_k: minimum_heavy.value(),
            maximum_heavy_temperature_k: maximum_heavy.value(),
            minimum_electron_temperature_k: minimum_electron.value(),
            maximum_electron_temperature_k: maximum_electron.value(),
            plasma_limits,
            relative_tolerance: 1e-6,
            absolute_amount_tolerance_mol: 1e-16,
            absolute_temperature_tolerance_k: 1e-3,
            energy_tolerance_j: 1e-8,
            balance_relative_tolerance: 1e-10,
            initial_step_s: 1e-6,
            minimum_step_s: 1e-18,
            maximum_step_s: 1.0,
            maximum_steps: 100_000,
        }
    }

    fn validate(self) -> Result<(), DtFusionError> {
        if [
            self.minimum_heavy_temperature_k,
            self.maximum_heavy_temperature_k,
            self.minimum_electron_temperature_k,
            self.maximum_electron_temperature_k,
            self.relative_tolerance,
            self.absolute_amount_tolerance_mol,
            self.absolute_temperature_tolerance_k,
            self.energy_tolerance_j,
            self.balance_relative_tolerance,
            self.initial_step_s,
            self.minimum_step_s,
            self.maximum_step_s,
        ]
        .iter()
        .any(|&x| !positive(x))
            || self.minimum_heavy_temperature_k >= self.maximum_heavy_temperature_k
            || self.minimum_electron_temperature_k >= self.maximum_electron_temperature_k
            || self.relative_tolerance > 0.1
            || self.balance_relative_tolerance > 1e-6
            || self.minimum_step_s > self.initial_step_s
            || self.initial_step_s > self.maximum_step_s
            || !(1..=1_000_000).contains(&self.maximum_steps)
        {
            return Err(DtFusionError::Invalid("fusion integration controls"));
        }
        Ok(())
    }
}

/// Refusal before publishing a new source-bound state; inputs are immutable.
#[derive(Clone, Debug, PartialEq)]
pub enum DtFusionError {
    /// Inconsistent source binding or invalid controls.
    Invalid(&'static str),
    /// Original reactivity source or product-law refusal.
    Reaction(FusionErrorV1),
    /// Original plasma source or physical-domain refusal.
    Plasma(DtFusionPlasmaErrorV1),
    /// Original Coulomb source or physical-domain refusal.
    Coulomb(CoulombErrorV1),
    /// Original radiation source or physical-domain refusal.
    Radiation(BremsstrahlungErrorV1),
    /// Final common ionized-gas state could not be published.
    State(FrozenIonizedGasStateErrorV1),
    /// The selected source/temperature/positivity domain was exhausted.
    OutsideDomain,
    /// The declared attempt or minimum-step budget was exhausted.
    StepBudget,
    /// A physical result could not be represented as finite f64.
    Arithmetic(&'static str),
    /// An independently evaluated energy or inventory account failed.
    Balance(&'static str),
    /// Cancellation before publication.
    Cancelled,
}

impl core::fmt::Display for DtFusionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "fusion burn: {self:?}")
    }
}
impl std::error::Error for DtFusionError {}

/// Actual interval accounts, distinct from source uncertainty or global error.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DtFusionBalance {
    /// D and T consumed, He4 created and neutrons emitted [mol of each].
    pub reaction_extent_mol: f64,
    /// Rest-mass energy converted by the admitted nuclear reaction [J].
    pub nuclear_energy_released_j: f64,
    /// Complete born-alpha kinetic energy deposited locally [J].
    pub alpha_energy_deposited_j: f64,
    /// Reaction-weighted reactant kinetic energy removed from heavy particles [J].
    pub reacting_ion_energy_removed_j: f64,
    /// Escaped neutron kinetic energy, including its thermal part [J].
    pub neutron_energy_escaped_j: f64,
    /// Signed collisional energy transferred from electrons to ions [J].
    /// Positive heats ions; negative heats electrons. This is not external heat.
    pub electron_to_heavy_coulomb_energy_j: f64,
    /// Electron-ion bremsstrahlung energy carried away by photons [J].
    pub photon_energy_escaped_j: f64,
    /// Constant externally supplied electron energy [J].
    pub external_electron_energy_j: f64,
    /// Constant externally supplied heavy-particle energy [J].
    pub external_heavy_energy_j: f64,
    /// Final source electron energy minus its independently accumulated budget [J].
    pub electron_energy_residual_j: f64,
    /// Final source heavy energy minus its independently accumulated budget [J].
    pub heavy_energy_residual_j: f64,
    /// Delta thermal U + escaping neutron/photon energy - Q - external heat [J].
    pub total_energy_residual_j: f64,
    /// Relative baryon inventory drift, including escaped neutrons.
    pub relative_baryon_residual: f64,
    /// Relative positive-ion charge inventory drift; electrons retain their amount.
    pub relative_charge_residual: f64,
    /// Largest accepted full/two-half local error norm, at most one.
    pub maximum_normalized_local_error: f64,
    /// Accepted two-half-step intervals.
    pub accepted_steps: usize,
    /// Rejected interval attempts.
    pub rejected_steps: usize,
}

/// Completed burn retaining the original source, initial state and controls.
#[derive(Clone, Debug, PartialEq)]
pub struct DtFusionResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    plasma: DtFusionPlasmaV1,
    heating: DtFusionHeating,
    transport: DtFusionTransport,
    controls: DtFusionSolve,
    elapsed: Time,
    balance: DtFusionBalance,
}

impl DtFusionResult {
    /// Original species bank with the accepted changed populations and temperatures.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Initial source-bound state; no source or population was modified in place.
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    /// Original nuclear-to-plasma species and caloric binding.
    pub const fn plasma(&self) -> &DtFusionPlasmaV1 {
        &self.plasma
    }
    /// Actual elapsed physical time.
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
    /// Retained thermalization and external-power assumptions.
    pub const fn heating(&self) -> DtFusionHeating {
        self.heating
    }
    /// Retained original thermal-exchange and emission sources.
    pub const fn transport(&self) -> &DtFusionTransport {
        &self.transport
    }
    /// Retained physical and numerical limits.
    pub const fn controls(&self) -> DtFusionSolve {
        self.controls
    }
    /// Complete independently checked interval accounts.
    pub const fn balance(&self) -> DtFusionBalance {
        self.balance
    }
    /// Escaped neutron amount, outside the material/chemical species bank.
    pub const fn neutron_amount(&self) -> Amount {
        Amount::new(self.balance.reaction_extent_mol)
    }
    /// Energy carried away by emitted neutrons; not locally deposited heat.
    pub const fn escaped_neutron_energy(&self) -> Energy {
        Energy::new(self.balance.neutron_energy_escaped_j)
    }
    /// Radiated energy leaving this explicitly optically thin parcel.
    pub const fn escaped_photon_energy(&self) -> Energy {
        Energy::new(self.balance.photon_energy_escaped_j)
    }
}

/// Advance finite D/T fuel and prescribed heating in a rigid parcel.
///
/// Every emitted neutron leaves this accounting boundary; every alpha remains
/// and thermalizes immediately. The selected channel is D+T only: rates for
/// other channels are available in L1 but are not silently combined here.
/// External powers are held over this explicitly declared interval. Reusing
/// the result in another physics owner requires that owner's source/domain
/// admission; this function does not extend weakly ionized EM laws to fusion.
/// Use [`advance_dt_fusion_with_transport`] to include collisional exchange
/// and optically thin photon losses in the same evolving energy account.
pub fn advance_dt_fusion(
    initial: &FrozenIonizedGasStateV1,
    plasma: &DtFusionPlasmaV1,
    heating: DtFusionHeating,
    duration: Time,
    controls: DtFusionSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<DtFusionResult, DtFusionError> {
    advance_dt_fusion_with_transport(
        initial,
        plasma,
        heating,
        &DtFusionTransport::default(),
        duration,
        controls,
        cancelled,
    )
}

/// Advance finite burn with explicitly selected collisional exchange and losses.
///
/// The same adaptive solve resolves fuel, ash, two-temperature feedback and
/// outgoing photon/neutron energy. Frozen reciprocal exchange is integrated
/// exponentially and free-free cooling with frozen coefficients analytically; symmetric
/// composition and midpoint source refresh are second order for smooth admitted
/// evolution. A cutoff crossover need not have a smooth local-error expansion.
/// Domain/cancellation refusal publishes no partial state. The step-doubling
/// estimate is not a global error or physical-accuracy certificate.
pub fn advance_dt_fusion_with_transport(
    initial: &FrozenIonizedGasStateV1,
    plasma: &DtFusionPlasmaV1,
    heating: DtFusionHeating,
    transport: &DtFusionTransport,
    duration: Time,
    controls: DtFusionSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DtFusionResult, DtFusionError> {
    controls.validate()?;
    poll(&mut cancelled)?;
    if !duration.value().is_finite()
        || duration.value() < 0.0
        || !heating.alpha_electron_fraction.is_finite()
        || !(0.0..=1.0).contains(&heating.alpha_electron_fraction)
        || !heating.electron_power.value().is_finite()
        || !heating.heavy_power.value().is_finite()
    {
        return Err(DtFusionError::Invalid(
            "finite duration, powers and alpha split in [0,1] required",
        ));
    }
    let problem = solve::Problem::new(
        initial,
        plasma,
        heating,
        transport,
        controls,
        &mut cancelled,
    )?;
    let (point, mut balance) = problem.integrate(duration.value(), &mut cancelled)?;
    let state = problem.publish(point, &mut balance, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(DtFusionResult {
        initial: initial.clone(),
        state,
        plasma: plasma.clone(),
        heating,
        transport: transport.clone(),
        controls,
        elapsed: duration,
        balance,
    })
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DtFusionError> {
    if cancelled() {
        Err(DtFusionError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
