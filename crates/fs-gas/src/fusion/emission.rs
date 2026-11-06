//! Typed handoff of accepted neutron accounts to independent transport owners.

use super::DtFusionResult;
use super::kinetic::KineticFusionResult;
use super::network::{FusionNetworkBalance, FusionNetworkResult};
use fs_qty::{Amount, Energy, Time};
use fs_thermochem::fusion::neutrons::{FusionNeutronEmissionV1, NeutronEmissionErrorV1};

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
