//! Shared-network evolution with an explicit isotropic alpha distribution.

use super::*;
use crate::fusion::kinetic::{AlphaEnergyGrid, MAX_ALPHA_BINS};
use fs_thermochem::coulomb::test_particle::TestParticleCoulombV1;

mod solve;
mod state;
pub use state::KineticFusionState;

/// Explicit dilute-alpha domain and bounded collision work. The retained source
/// supplies the original projectile and background masses, charges and cutoffs.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticAlphaTransport {
    collisions: TestParticleCoulombV1,
    maximum_number_ratio: f64,
    maximum_energy_ratio: f64,
    maximum_upper_bin_fraction: f64,
    maximum_matrix_applications: usize,
}
impl KineticAlphaTransport {
    /// Limits on alpha/thermal-ion number, alpha/thermal-pool energy, and the
    /// fraction occupying the uppermost energy bin must lie in (0,1). The last
    /// diagnostic screens an inadequately extended grid; it is not a rigorous
    /// bound on energy-discretization error. Applications bound each frozen
    /// collision substep, and the network separately bounds attempted steps.
    pub fn new(
        collisions: TestParticleCoulombV1,
        maximum_number_ratio: f64,
        maximum_energy_ratio: f64,
        maximum_upper_bin_fraction: f64,
        maximum_matrix_applications: usize,
    ) -> Result<Self, FusionNetworkError> {
        if [
            maximum_number_ratio,
            maximum_energy_ratio,
            maximum_upper_bin_fraction,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0 || *x >= 1.0)
            || maximum_matrix_applications == 0
            || maximum_matrix_applications > 1_000_000
        {
            return Err(FusionNetworkError::Invalid(
                "finite dilute-alpha limits in (0,1), bounded matrix work",
            ));
        }
        Ok(Self {
            collisions,
            maximum_number_ratio,
            maximum_energy_ratio,
            maximum_upper_bin_fraction,
            maximum_matrix_applications,
        })
    }
    pub const fn collisions(&self) -> &TestParticleCoulombV1 {
        &self.collisions
    }
    pub const fn maximum_number_ratio(&self) -> f64 {
        self.maximum_number_ratio
    }
    pub const fn maximum_energy_ratio(&self) -> f64 {
        self.maximum_energy_ratio
    }
    pub const fn maximum_upper_bin_fraction(&self) -> f64 {
        self.maximum_upper_bin_fraction
    }
    pub const fn maximum_matrix_applications(&self) -> usize {
        self.maximum_matrix_applications
    }
}

/// Alpha storage and reciprocal background-heating accounts for one interval.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct KineticAlphaBalance {
    /// Independently accumulated He4 births from DT and DHe3 [mol].
    pub born_alpha_amount_mol: f64,
    /// Complete born-alpha kinetic energy, including incoming thermal energy [J].
    pub born_alpha_energy_j: f64,
    /// Signed collisional heat to background ions [J].
    pub deposited_heavy_energy_j: f64,
    /// Signed collisional heat to electrons [J].
    pub deposited_electron_energy_j: f64,
    /// Change of the retained alpha distribution's kinetic storage [J].
    pub alpha_energy_change_j: f64,
    /// Delta alpha energy + heat to both backgrounds - born energy [J].
    pub alpha_energy_residual_j: f64,
    /// Distribution amount minus independently accumulated alpha inventory [mol].
    pub alpha_amount_residual_mol: f64,
    /// Final number fraction in the uppermost authored energy bin.
    pub upper_bin_fraction: f64,
}

/// Accepted hybrid evolution. Its state keeps alpha kinetic pressure/storage
/// separate from the Maxwellian background; it is not a frozen thermal state.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticFusionResult {
    initial: KineticFusionState,
    state: KineticFusionState,
    channels: FusionNetworkChannels,
    heating: FusionNetworkHeating,
    thermal_transport: FusionNetworkTransport,
    alpha_transport: KineticAlphaTransport,
    controls: FusionNetworkSolve,
    elapsed: Time,
    balance: FusionNetworkBalance,
    alpha_balance: KineticAlphaBalance,
}
impl KineticFusionResult {
    pub const fn state(&self) -> &KineticFusionState {
        &self.state
    }
    pub const fn initial_state(&self) -> &KineticFusionState {
        &self.initial
    }
    pub const fn channels(&self) -> FusionNetworkChannels {
        self.channels
    }
    pub const fn heating(&self) -> FusionNetworkHeating {
        self.heating
    }
    pub const fn thermal_transport(&self) -> &FusionNetworkTransport {
        &self.thermal_transport
    }
    pub const fn alpha_transport(&self) -> &KineticAlphaTransport {
        &self.alpha_transport
    }
    pub const fn controls(&self) -> FusionNetworkSolve {
        self.controls
    }
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
    /// Shared nuclear and thermal accounts. Charged-product deposition fields
    /// contain only prompt non-alpha products; alpha heat is in alpha_balance.
    /// The total energy residual includes the change in alpha kinetic storage.
    pub const fn balance(&self) -> FusionNetworkBalance {
        self.balance
    }
    pub const fn alpha_balance(&self) -> KineticAlphaBalance {
        self.alpha_balance
    }
}

/// Evolve thermal fusion with finite-time alpha energy drift AND diffusion.
///
/// DT/DHe3 alphas enter the retained energy grid without prompt deposition or
/// thermal heat capacity. They remain a separate isotropic kinetic distribution
/// through the thermal tail; no zero-energy sink or arbitrary handoff is used.
/// DD-born p/T/He3 and the DHe3 proton retain explicit prompt deposition. Alpha
/// entries in `heating.product_electron_fractions` must be zero, because the
/// collision source determines their energy partition. Neutrons escape.
/// Each channel is injected at its source event's mean product energy; a
/// resolved thermal/Doppler birth-energy spectrum is outside this closure.
///
/// The background remains Maxwellian and the alpha population must pass explicit
/// dilute number/energy bounds. Both collision sources must use ElectronDebye
/// screening so no thermal susceptibility is assigned to kinetic alphas.
/// The original all-thermal plasma gate admits the source/charge reference;
/// it does not certify classical occupation of an arbitrary authored kinetic
/// distribution. Classical dilute test alphas remain a model assumption, and
/// the number/energy bounds do not quantify neglected alpha self-collisions.
/// Radiation includes all ions only after separately admitting the fastest
/// represented alpha against its stationary-ion velocity-variance limit.
///
/// This first kinetic owner has rigid volume. It resolves no fast-secondary
/// reactions, alpha self-collisions, pitch angles, confinement, magnetic fields
/// or neutron transport. Numerical time error and energy-grid error are distinct;
/// refine the authored grid for distribution and heating predictions.
pub fn advance_kinetic_fusion(
    initial: &KineticFusionState,
    channels: FusionNetworkChannels,
    heating: FusionNetworkHeating,
    thermal_transport: &FusionNetworkTransport,
    alpha_transport: &KineticAlphaTransport,
    duration: Time,
    controls: FusionNetworkSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<KineticFusionResult, FusionNetworkError> {
    controls.validate()?;
    heating.validate()?;
    poll(&mut cancelled)?;
    if !duration.value().is_finite()
        || duration.value() < 0.0
        || heating.product_electron_fractions[0][0] != 0.0
        || heating.product_electron_fractions[3][0] != 0.0
    {
        return Err(FusionNetworkError::Invalid(
            "finite duration and zero prompt alpha fractions required",
        ));
    }
    let problem = solve::Problem::new(
        initial,
        channels,
        heating,
        thermal_transport,
        alpha_transport,
        controls,
        &mut cancelled,
    )?;
    let (point, mut balance) = problem.integrate(duration.value(), &mut cancelled)?;
    let (state, alpha_balance) = problem.publish(point, &mut balance, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(KineticFusionResult {
        initial: initial.clone(),
        state,
        channels,
        heating,
        thermal_transport: thermal_transport.clone(),
        alpha_transport: alpha_transport.clone(),
        controls,
        elapsed: duration,
        balance,
        alpha_balance,
    })
}

#[cfg(test)]
mod tests;
