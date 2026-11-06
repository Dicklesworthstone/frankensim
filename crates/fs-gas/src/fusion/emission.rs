//! Typed handoff of accepted neutron accounts to independent transport owners.

use super::kinetic::KineticFusionResult;
use super::network::{FusionNetworkBalance, FusionNetworkResult};
use super::solve::Sum;
use super::{DtFusionError, DtFusionResult};
use fs_qty::{Amount, Energy, Time};
use fs_thermochem::fusion::neutrons::{
    FusionNeutronEmissionV1, FusionNeutronIntervalV1, NeutronEmissionErrorV1,
};

pub(super) type NeutronHistory = Option<Vec<FusionNeutronIntervalV1>>;

/// Direct moments within one attempted adaptive interval. Trial copies discard
/// predictor/coarse/rejected contributions; the accepted two-half path retains
/// both halves. Resetting these accounts avoids subtracting nearly equal totals
/// when a weak late source follows a much larger earlier emission.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct NeutronMoments {
    amounts: [Sum; 2],
    energies: [Sum; 2],
}

impl NeutronMoments {
    pub(super) fn add(&mut self, channel: usize, amount_mol: f64, energy_j: f64) {
        self.amounts[channel].add(amount_mol);
        self.energies[channel].add(energy_j);
    }
}

/// Successful adaptive intervals alone may grow the optional output. Numerical
/// trials remain allocation-free, and callers retain this vector locally until
/// the complete state and source history can be published together.
pub(super) fn record_interval(
    history: &mut NeutronHistory,
    start: f64,
    end: f64,
    moments: Option<NeutronMoments>,
) -> Result<(), DtFusionError> {
    if let Some(history) = history {
        let moments = moments.ok_or(DtFusionError::Arithmetic(
            "accepted neutron interval moments missing",
        ))?;
        let emission = FusionNeutronEmissionV1::new(
            Time::new(end - start),
            moments.amounts.map(|n| Amount::new(n.value)),
            moments.energies.map(|e| Energy::new(e.value)),
        )
        .map_err(DtFusionError::NeutronEmission)?;
        let interval = FusionNeutronIntervalV1::new(Time::new(start), Time::new(end), emission)
            .map_err(DtFusionError::NeutronEmission)?;
        history
            .try_reserve(1)
            .map_err(|_| DtFusionError::Allocation("accepted neutron history"))?;
        history.push(interval);
    }
    Ok(())
}

impl DtFusionResult {
    /// Channel-resolved neutron emission over this accepted burn interval.
    /// Counts and kinetic energy come directly from the accepted plasma debit;
    /// downstream transport does not re-credit this parcel or assume deposition.
    /// D–D is explicitly empty. No spectrum or time profile is manufactured.
    pub fn neutron_emission(&self) -> Result<FusionNeutronEmissionV1, NeutronEmissionErrorV1> {
        FusionNeutronEmissionV1::new(
            self.elapsed_time(),
            [self.neutron_amount(), Amount::new(0.0)],
            [self.escaped_neutron_energy(), Energy::new(0.0)],
        )
    }
}

impl FusionNetworkResult {
    /// Retain separate D–T and D–D neutron number/kinetic-energy accounts from
    /// the accepted interval. The original escaped-energy debit stays intact.
    /// A transport consumer must explicitly choose its time/energy projection.
    pub fn neutron_emission(&self) -> Result<FusionNeutronEmissionV1, NeutronEmissionErrorV1> {
        network_emission(self.elapsed_time(), self.balance())
    }
}

impl KineticFusionResult {
    /// Retain separate D–T and D–D neutron emission from the same accepted
    /// hybrid evolution as the thermal and kinetic-alpha accounts. Neutron
    /// energy stays outside the plasma, independently of alpha slowing heat.
    /// This handoff supplies integrated moments, not a neutron spectrum.
    pub fn neutron_emission(&self) -> Result<FusionNeutronEmissionV1, NeutronEmissionErrorV1> {
        network_emission(self.elapsed_time(), self.balance())
    }
}

fn network_emission(
    duration: Time,
    balance: FusionNetworkBalance,
) -> Result<FusionNeutronEmissionV1, NeutronEmissionErrorV1> {
    // FusionReactionV1::ALL order is [DT, DDn, DDp, DHe3]. Only the first
    // two channels emit neutrons, each in product slot one. Preserve their
    // original accepted extents rather than splitting the aggregate account.
    FusionNeutronEmissionV1::new(
        duration,
        [
            Amount::new(balance.reaction_extents_mol[0]),
            Amount::new(balance.reaction_extents_mol[1]),
        ],
        [
            Energy::new(balance.product_kinetic_energies_j[0][1]),
            Energy::new(balance.product_kinetic_energies_j[1][1]),
        ],
    )
}
