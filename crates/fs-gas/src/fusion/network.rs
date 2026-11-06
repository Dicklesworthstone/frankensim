//! Competing light-ion fusion in one original two-temperature material state.
//!
//! The explicit channel set selects D–T, D–D→He3+n, D–D→T+p and D–He3.
//! The two D–D rates already describe separate branches: each receives the
//! identical-pair factor once, without an extra assumed branching fraction.
//! Tritium and He3 born during a step can subsequently react; no trace seed is
//! inserted. All charged products remain and thermalize instantaneously with
//! an explicit per-product electron/ion split. Every neutron escapes.
//!
//! Exact positive finite-stock reaction maps are composed symmetrically at
//! midpoint coefficients. The shared Coulomb and radiation kernels act on the
//! same actual populations. Bounded stack trials, compensated physical accounts
//! and step doubling resolve competition and temperature feedback. Second-order
//! refinement applies to smooth admitted evolution; the local estimate is not
//! a global bound or experimental validation of this thermal closure.

use super::volume::FusionVolumeHistory;
use super::{AVOGADRO, CV_MOLAR, poll, positive};
use fs_qty::{Amount, Energy, Power, Time};
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;
use fs_thermochem::fusion::neutrons::FusionNeutronIntervalV1;
use fs_thermochem::fusion::plasma::LightIonFusionPlasmaV1;
use fs_thermochem::fusion::{FusionReactionV1, NuclearSpeciesV1};

/// Shared source, domain, budget, balance and cancellation refusals.
pub use super::DtFusionError as FusionNetworkError;
/// Shared physical/numerical controls, including the original plasma limits.
pub use super::DtFusionSolve as FusionNetworkSolve;
/// Shared Coulomb/radiation sources and optional thermal energy confinement.
pub use super::DtFusionTransport as FusionNetworkTransport;

pub(crate) mod kinetic;
mod solve;

/// Explicit subset of the four published channels, in `FusionReactionV1::ALL` order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FusionNetworkChannels {
    enabled: [bool; 4],
}
impl FusionNetworkChannels {
    /// All four supported channels, including both independent D–D branches.
    pub const ALL: Self = Self { enabled: [true; 4] };
    /// No nuclear reactions; selected thermal sources and external powers remain.
    pub const NONE: Self = Self {
        enabled: [false; 4],
    };
    /// Select channels explicitly; no reaction is inferred from a species name.
    pub const fn new(enabled: [bool; 4]) -> Self {
        Self { enabled }
    }
    /// Selected channels in the original fixed source order.
    pub const fn enabled(self) -> [bool; 4] {
        self.enabled
    }
    /// Whether a particular published reaction participates.
    pub const fn contains(self, reaction: FusionReactionV1) -> bool {
        self.enabled[channel_index(reaction)]
    }
}

/// Instantaneous local charged-product deposition and held external powers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FusionNetworkHeating {
    /// Electron-deposited fractions, indexed by `FusionReactionV1::ALL` and
    /// each reaction's `products()` order. Charged fractions lie in [0,1].
    /// Neutron fractions must be exactly zero: this owner has no neutron
    /// deposition model. Fractions describe a chosen closure, not stopping data.
    pub product_electron_fractions: [[f64; 2]; 4],
    /// Constant external electron power [W], negative for prescribed cooling.
    pub electron_power: Power,
    /// Constant external ion power [W], negative for prescribed cooling.
    pub heavy_power: Power,
}
impl FusionNetworkHeating {
    /// Use one explicitly authored split for every charged product. Admission
    /// still checks the fraction and powers when an evolution is requested.
    pub fn uniform_charged_fraction(
        fraction: f64,
        electron_power: Power,
        heavy_power: Power,
    ) -> Self {
        Self {
            product_electron_fractions: core::array::from_fn(|r| {
                FusionReactionV1::ALL[r].products().map(|p| {
                    if p == NuclearSpeciesV1::Neutron {
                        0.0
                    } else {
                        fraction
                    }
                })
            }),
            electron_power,
            heavy_power,
        }
    }
    fn validate(self) -> Result<(), FusionNetworkError> {
        if !self.electron_power.value().is_finite() || !self.heavy_power.value().is_finite() {
            return Err(FusionNetworkError::Invalid(
                "finite network external powers required",
            ));
        }
        for (r, reaction) in FusionReactionV1::ALL.iter().enumerate() {
            for (p, product) in reaction.products().iter().enumerate() {
                let fraction = self.product_electron_fractions[r][p];
                if !fraction.is_finite()
                    || !(0.0..=1.0).contains(&fraction)
                    || (*product == NuclearSpeciesV1::Neutron && fraction != 0.0)
                {
                    return Err(FusionNetworkError::Invalid(
                        "charged-product fractions in [0,1] and zero neutron deposition required",
                    ));
                }
            }
        }
        Ok(())
    }
}

/// Accepted interval accounts; nuclear rest energy is outside gas calorics.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FusionNetworkBalance {
    /// Event amounts [mol] in `FusionReactionV1::ALL` order; DD consumes twice its extent.
    pub reaction_extents_mol: [f64; 4],
    /// Born product kinetic energies [J], by channel and original product order.
    /// Includes reaction-weighted incoming ion kinetic energy as well as Q.
    pub product_kinetic_energies_j: [[f64; 2]; 4],
    /// Nuclear mass-defect energy liberated over all selected channels [J].
    pub nuclear_energy_released_j: f64,
    /// Charged-product energy deposited into electrons [J].
    pub charged_product_electron_energy_j: f64,
    /// Charged-product energy deposited into the heavy pool [J].
    pub charged_product_heavy_energy_j: f64,
    /// Reaction-weighted kinetic energy removed with reacting ions [J].
    pub reacting_ion_energy_removed_j: f64,
    /// Emitted neutron amount [mol], outside the material species bank.
    pub neutron_amount_mol: f64,
    /// Escaping neutron kinetic energy, including its incoming-thermal share [J].
    pub neutron_energy_escaped_j: f64,
    /// Signed electron-to-ion collisional transfer; cancels internally [J].
    pub electron_to_heavy_coulomb_energy_j: f64,
    /// Free-free emission escaping the explicitly optically thin boundary [J].
    pub photon_energy_escaped_j: f64,
    /// Signed ion thermal energy transported to its prescribed exterior [J].
    /// Negative means exterior heating. Neither fuel nor kinetic alphas escape.
    pub heavy_transport_energy_out_j: f64,
    /// Signed electron thermal energy transported to its prescribed exterior [J].
    pub electron_transport_energy_out_j: f64,
    /// Held external electron energy [J].
    pub external_electron_energy_j: f64,
    /// Held external ion energy [J].
    pub external_heavy_energy_j: f64,
    /// Mechanical work on ions [J]; positive during compression.
    pub boundary_work_on_heavy_j: f64,
    /// Mechanical work on electrons [J]; positive during compression.
    pub boundary_work_on_electrons_j: f64,
    /// Final source electron energy minus its independent transfer account [J].
    pub electron_energy_residual_j: f64,
    /// Final source ion energy minus its independent transfer account [J].
    pub heavy_energy_residual_j: f64,
    /// Delta thermal U + neutron/photon/transport energy out - Q - external heat
    /// - work on ions - work on electrons [J].
    pub total_energy_residual_j: f64,
    /// Relative baryon drift, including all emitted neutrons.
    pub relative_baryon_residual: f64,
    /// Relative positive-ion charge drift; electron amount is fixed.
    pub relative_charge_residual: f64,
    /// Largest absolute species residual against independent channel extents [mol].
    pub maximum_stoichiometry_residual_mol: f64,
    /// Largest accepted full/two-half-step normalized local error, at most one.
    pub maximum_normalized_local_error: f64,
    /// Accepted two-half-step intervals.
    pub accepted_steps: usize,
    /// Rejected attempts.
    pub rejected_steps: usize,
}

/// Published network evolution retaining the common original material sources.
#[derive(Clone, Debug, PartialEq)]
pub struct FusionNetworkResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    plasma: LightIonFusionPlasmaV1,
    channels: FusionNetworkChannels,
    heating: FusionNetworkHeating,
    transport: FusionNetworkTransport,
    volume_history: Option<FusionVolumeHistory>,
    controls: FusionNetworkSolve,
    elapsed: Time,
    balance: FusionNetworkBalance,
    neutron_history: Option<Vec<FusionNeutronIntervalV1>>,
}
impl FusionNetworkResult {
    /// Accepted original-bank state, reusable by other admitted physics owners.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Unmodified input state.
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    /// Retained original nuclear-to-material binding.
    pub const fn plasma(&self) -> &LightIonFusionPlasmaV1 {
        &self.plasma
    }
    /// Explicit selected reaction subset.
    pub const fn channels(&self) -> FusionNetworkChannels {
        self.channels
    }
    /// Retained per-product deposition and external powers.
    pub const fn heating(&self) -> FusionNetworkHeating {
        self.heating
    }
    /// Retained common Coulomb/radiation sources.
    pub const fn transport(&self) -> &FusionNetworkTransport {
        &self.transport
    }
    /// Original prescribed geometry, or `None` for the rigid-volume entry point.
    pub fn volume_history(&self) -> Option<&FusionVolumeHistory> {
        self.volume_history.as_ref()
    }
    /// Retained domain and numerical controls.
    pub const fn controls(&self) -> FusionNetworkSolve {
        self.controls
    }
    /// Completed interval duration.
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
    /// Independently checked particle and energy accounts.
    pub const fn balance(&self) -> FusionNetworkBalance {
        self.balance
    }
    /// Direct D–T and D–D neutron moments over accepted relative-time intervals.
    /// Recording is opt-in; disabled solves return `None`, and recorded
    /// zero-duration solves return an empty slice. Only the accepted two-half
    /// solution contributes, and no interval crosses a prescribed volume knot.
    /// A downstream source must declare its intra-interval time projection.
    pub fn neutron_history(&self) -> Option<&[FusionNeutronIntervalV1]> {
        self.neutron_history.as_deref()
    }
    /// Emitted neutron inventory for an explicitly separate transport consumer.
    pub const fn neutron_amount(&self) -> Amount {
        Amount::new(self.balance.neutron_amount_mol)
    }
    /// Kinetic energy carried outside this parcel by neutrons.
    pub const fn escaped_neutron_energy(&self) -> Energy {
        Energy::new(self.balance.neutron_energy_escaped_j)
    }
    /// Photon energy leaving this optically thin parcel.
    pub const fn escaped_photon_energy(&self) -> Energy {
        Energy::new(self.balance.photon_energy_escaped_j)
    }
    /// Signed total mechanical work on both thermal pools [J].
    pub fn total_boundary_work(&self) -> Energy {
        Energy::new(
            self.balance.boundary_work_on_heavy_j + self.balance.boundary_work_on_electrons_j,
        )
    }
}

/// Evolve an explicitly selected network with shared thermal feedback.
///
/// Every enabled channel must support the entire authored ion-temperature
/// interval, including channels whose fuel is initially absent. For all four
/// channels that source intersection is 0.5–100 keV. All charged products are
/// immediately Maxwellian at common ion temperature; the two thermal pools,
/// rigid volume, neutron escape and optically thin radiation are explicit
/// closures. This is not a fast-particle, confinement or radiation-transport solve.
/// Failure publishes no partial state and never modifies any input source.
pub fn advance_fusion_network(
    initial: &FrozenIonizedGasStateV1,
    plasma: &LightIonFusionPlasmaV1,
    channels: FusionNetworkChannels,
    heating: FusionNetworkHeating,
    transport: &FusionNetworkTransport,
    duration: Time,
    controls: FusionNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<FusionNetworkResult, FusionNetworkError> {
    advance(
        initial, plasma, channels, heating, transport, duration, None, controls, cancelled,
    )
}

/// Evolve the same fusion network through prescribed compression and expansion.
///
/// The history starts at the input state's exact volume and supplies the full
/// interval duration. Every stage uses its current volume for reaction exposure,
/// Coulomb exchange, radiation and physical-domain admission. For each thermal
/// pool, reversible isotropic work obeys `dU = -(2/3) U d(ln V)` at fixed stock.
/// Positive work on the gas is compression; heat and nuclear energy remain
/// separate accounts. No source clock is advanced by a geometric substep.
///
/// Exact adiabatic maps surround the midpoint source solve, and adaptive steps
/// stop at every authored knot. Work increments participate in the local error
/// estimate. Smooth segments have second-order refinement; this estimate is
/// not a global or experimental error certificate. The original history and
/// material sources are retained in the result. Boundary dynamics, shocks,
/// anisotropic magnetic work and fast-particle confinement require other models.
pub fn advance_fusion_network_with_volume(
    initial: &FrozenIonizedGasStateV1,
    plasma: &LightIonFusionPlasmaV1,
    channels: FusionNetworkChannels,
    heating: FusionNetworkHeating,
    transport: &FusionNetworkTransport,
    volume_history: &FusionVolumeHistory,
    controls: FusionNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<FusionNetworkResult, FusionNetworkError> {
    advance(
        initial,
        plasma,
        channels,
        heating,
        transport,
        volume_history.duration(),
        Some(volume_history),
        controls,
        cancelled,
    )
}

fn advance(
    initial: &FrozenIonizedGasStateV1,
    plasma: &LightIonFusionPlasmaV1,
    channels: FusionNetworkChannels,
    heating: FusionNetworkHeating,
    transport: &FusionNetworkTransport,
    duration: Time,
    volume_history: Option<&FusionVolumeHistory>,
    controls: FusionNetworkSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FusionNetworkResult, FusionNetworkError> {
    controls.validate()?;
    heating.validate()?;
    poll(&mut cancelled)?;
    if !duration.value().is_finite() || duration.value() < 0.0 {
        return Err(FusionNetworkError::Invalid(
            "finite nonnegative network duration required",
        ));
    }
    if volume_history.is_some_and(|history| history.initial_volume() != initial.volume()) {
        return Err(FusionNetworkError::Invalid(
            "volume history must start at the original state volume",
        ));
    }
    let problem = solve::Problem::new(
        initial,
        plasma,
        channels,
        heating,
        transport,
        volume_history,
        controls,
        &mut cancelled,
    )?;
    let (point, mut balance, neutron_history) =
        problem.integrate(duration.value(), &mut cancelled)?;
    let state = problem.publish(point, &mut balance, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(FusionNetworkResult {
        initial: initial.clone(),
        state,
        plasma: plasma.clone(),
        channels,
        heating,
        transport: transport.clone(),
        volume_history: volume_history.cloned(),
        controls,
        elapsed: duration,
        balance,
        neutron_history,
    })
}

const fn channel_index(reaction: FusionReactionV1) -> usize {
    match reaction {
        FusionReactionV1::DeuteriumTritium => 0,
        FusionReactionV1::DeuteriumDeuteriumNeutron => 1,
        FusionReactionV1::DeuteriumDeuteriumProton => 2,
        FusionReactionV1::DeuteriumHelium3 => 3,
    }
}
const fn species_slot(species: NuclearSpeciesV1) -> usize {
    match species {
        NuclearSpeciesV1::Proton => 0,
        NuclearSpeciesV1::Deuterium => 1,
        NuclearSpeciesV1::Tritium => 2,
        NuclearSpeciesV1::Helium3 => 3,
        NuclearSpeciesV1::Helium4 => 4,
        NuclearSpeciesV1::Neutron => 5,
    }
}

#[cfg(test)]
mod tests;
