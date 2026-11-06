//! Interval-integrated neutron emission from the supported fusion channels.
//!
//! D–T and D–D→He3+n retain separate particle and kinetic-energy accounts.
//! Each event creates one neutron. Amount-to-count conversion uses the exact
//! SI Avogadro constant, subject only to binary64 representation. Fractional
//! expected particle counts are admitted for continuum source models.
//!
//! These are integrated moments, not a neutron spectrum or a source history.
//! A consumer must explicitly choose its time/energy projection and transport
//! domain. In particular, emitted energy is not automatically deposited heat.

use super::FusionReactionV1;
use crate::electron_collision::AVOGADRO_PER_MOL;
use fs_qty::{Amount, Energy, Time};

/// Invalid interval accounts or an unrepresentable number/energy conversion.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NeutronEmissionErrorV1 {
    /// Duration must be finite and nonnegative, and positive for any emission.
    InvalidDuration,
    /// A channel amount is negative or nonfinite.
    InvalidAmount(FusionReactionV1),
    /// A channel's integrated kinetic energy is negative or nonfinite.
    InvalidEnergy(FusionReactionV1),
    /// Positive neutron amount and positive kinetic energy must occur together.
    InconsistentChannel(FusionReactionV1),
    /// A count, mean particle energy or combined account is not representable.
    Arithmetic(&'static str),
}

impl core::fmt::Display for NeutronEmissionErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "fusion neutron emission: {self:?}")
    }
}
impl std::error::Error for NeutronEmissionErrorV1 {}

/// Integrated neutron number and kinetic energy for one named fusion channel.
/// Only the emission constructor creates these immutable validated accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FusionNeutronChannelV1 {
    reaction: FusionReactionV1,
    amount: Amount,
    kinetic_energy: Energy,
    particle_count: f64,
    mean_particle_energy: Option<Energy>,
}

impl FusionNeutronChannelV1 {
    fn new(
        reaction: FusionReactionV1,
        amount: Amount,
        kinetic_energy: Energy,
    ) -> Result<Self, NeutronEmissionErrorV1> {
        use NeutronEmissionErrorV1 as E;
        if !amount.value().is_finite() || amount.value() < 0.0 {
            return Err(E::InvalidAmount(reaction));
        }
        if !kinetic_energy.value().is_finite() || kinetic_energy.value() < 0.0 {
            return Err(E::InvalidEnergy(reaction));
        }
        if (amount.value() == 0.0) != (kinetic_energy.value() == 0.0) {
            return Err(E::InconsistentChannel(reaction));
        }
        let particle_count = amount.value() * AVOGADRO_PER_MOL;
        let mean_particle_energy = if amount.value() == 0.0 {
            None
        } else {
            if !particle_count.is_finite() || particle_count <= 0.0 {
                return Err(E::Arithmetic("channel particle count"));
            }
            // Divide by the already admitted count, avoiding a possibly
            // overflowing intermediate molar energy in E / amount / N_A.
            let mean = kinetic_energy.value() / particle_count;
            if !mean.is_finite() || mean <= 0.0 {
                return Err(E::Arithmetic("channel mean neutron kinetic energy"));
            }
            Some(Energy::new(mean))
        };
        Ok(Self {
            reaction,
            amount,
            kinetic_energy,
            particle_count,
            mean_particle_energy,
        })
    }

    /// The original D–T or D–D neutron-producing reaction identity.
    pub const fn reaction(&self) -> FusionReactionV1 {
        self.reaction
    }

    /// Neutrons emitted during the complete interval [mol].
    pub const fn amount(&self) -> Amount {
        self.amount
    }

    /// Complete emitted kinetic energy [J], including reacting-ion energy.
    pub const fn kinetic_energy(&self) -> Energy {
        self.kinetic_energy
    }

    /// Expected neutron count, using the exact SI Avogadro constant [1/mol].
    pub const fn particle_count(&self) -> f64 {
        self.particle_count
    }

    /// Interval-weighted kinetic energy per emitted neutron [J/particle].
    /// An empty channel has no mean; no zero-energy neutron is manufactured.
    pub const fn mean_particle_energy(&self) -> Option<Energy> {
        self.mean_particle_energy
    }
}

/// Immutable, channel-resolved fusion neutron emission over one accepted time
/// interval. The fixed channel order is D–T followed by D–D→He3+n.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FusionNeutronEmissionV1 {
    duration: Time,
    channels: [FusionNeutronChannelV1; 2],
    total_amount: Amount,
    total_kinetic_energy: Energy,
    total_particle_count: f64,
}

impl FusionNeutronEmissionV1 {
    /// Bind integrated amounts and kinetic energies, in [D–T, D–D→He3+n]
    /// order, to their physical interval. Both channel accounts must be zero
    /// together or strictly positive together. A zero-duration empty interval
    /// is admitted; nonempty emission needs a positive duration.
    ///
    /// The constructor checks number/energy arithmetic without assigning an
    /// emission rate, spectrum, location, direction, attenuation or deposition.
    /// It does not independently certify the producing burn calculation.
    pub fn new(
        duration: Time,
        amounts: [Amount; 2],
        energies: [Energy; 2],
    ) -> Result<Self, NeutronEmissionErrorV1> {
        use NeutronEmissionErrorV1 as E;
        if !duration.value().is_finite() || duration.value() < 0.0 {
            return Err(E::InvalidDuration);
        }
        let channels = [
            FusionNeutronChannelV1::new(
                FusionReactionV1::DeuteriumTritium,
                amounts[0],
                energies[0],
            )?,
            FusionNeutronChannelV1::new(
                FusionReactionV1::DeuteriumDeuteriumNeutron,
                amounts[1],
                energies[1],
            )?,
        ];
        let total_amount = amounts[0].value() + amounts[1].value();
        let total_kinetic_energy = energies[0].value() + energies[1].value();
        let total_particle_count = channels[0].particle_count + channels[1].particle_count;
        if !total_amount.is_finite()
            || !total_kinetic_energy.is_finite()
            || !total_particle_count.is_finite()
        {
            return Err(E::Arithmetic("combined neutron number or energy"));
        }
        if total_amount > 0.0 && duration.value() == 0.0 {
            return Err(E::InvalidDuration);
        }
        Ok(Self {
            duration,
            channels,
            total_amount: Amount::new(total_amount),
            total_kinetic_energy: Energy::new(total_kinetic_energy),
            total_particle_count,
        })
    }

    /// Exact accepted duration attached to these integrated accounts [s].
    pub const fn duration(&self) -> Time {
        self.duration
    }

    /// Separate [D–T, D–D→He3+n] accounts, including empty channels.
    pub const fn channels(&self) -> &[FusionNeutronChannelV1; 2] {
        &self.channels
    }

    /// Sum of both emitted neutron amounts [mol].
    pub const fn total_amount(&self) -> Amount {
        self.total_amount
    }

    /// Sum of the two emitted kinetic-energy accounts [J].
    pub const fn total_kinetic_energy(&self) -> Energy {
        self.total_kinetic_energy
    }

    /// Sum of the separately converted channel particle counts.
    pub const fn total_particle_count(&self) -> f64 {
        self.total_particle_count
    }
}

#[cfg(test)]
mod tests;
