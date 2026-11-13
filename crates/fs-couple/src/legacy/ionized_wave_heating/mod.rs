//! Magnetized Maxwell absorption on the original finite two-temperature gas.
//!
//! The prepared circular fields supply their exact signed electron/heavy
//! friction powers once. The existing gas energy owner advances the original
//! caloric sources and neutral/Coulomb exchange at actual temperatures. The
//! endpoint refreshes the original carrier drag and bound-charge background.
//!
//! The top-level owner is a uniformly mixed, fixed-volume control parcel
//! inside the original homogeneous half-space. Field and absorption rates are
//! held over an explicit cycle-averaged heating interval. A refreshed result
//! seeds a NEW Maxwell solve; it never mutates the original wave or promises a
//! continuously coupled optical/thermal solution. There is no second reflecting
//! boundary, spatial conduction, flow or fast carrier transient. `advance`
//! freezes populations; `advance_reacting` spends the same work through the
//! native atomic ionization/recombination network before refreshing transport.
//!
//! [`layers`] couples contiguous finite coherent gas films to the shared
//! multicomponent diffuse-radiation network, including carrier birth, original
//! momentum/background refresh and explicitly selected thin free-free cooling.

use fs_em::plane_wave::plasma::magnetized::MagnetizedIonizedGasSlabAbsorption;
use fs_gas::chemistry::energy::{
    ElectronHeavyEnergyError, ElectronHeavyEnergyResult, ElectronHeavyEnergySolve,
    advance_frozen_energy, advance_frozen_energy_with_coulomb,
};
use fs_gas::chemistry::ionization_network::IonizationNetworkError;
use fs_matdb::MaterialCard;
use fs_material::{
    ionized_gas::{CarrierMomentumSource, IonizedGasConductivity, IonizedGasRefreshError},
    state_point::{MaterialStatePointError, ResolvedMaterialStatePoint},
};
use fs_qty::{Power, Time};
use fs_thermochem::SpeciesId;
use fs_thermochem::electron_collision::ElectronCrossSectionV1;

mod reacting;
pub use reacting::{IonizedGasReactingHeatControls, MagnetizedIonizedGasReactingHeatStep};

/// Coherent finite gas films with shared diffuse atomic radiation.
pub mod layers;

/// Refusal leaves the original wave, material sources and gas unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum IonizedWaveHeatingError {
    Invalid(&'static str),
    Gas(ElectronHeavyEnergyError),
    Reaction(IonizationNetworkError),
    Transport(IonizedGasRefreshError),
    Background(MaterialStatePointError),
    TemperatureIncrement {
        heavy_k: f64,
        electron_k: f64,
        maximum_k: f64,
    },
    PopulationIncrement {
        species: SpeciesId,
        change_mol: f64,
        maximum_mol: f64,
    },
    Cancelled,
}
impl core::fmt::Display for IonizedWaveHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ionized wave heating: {self:?}")
    }
}
impl std::error::Error for IonizedWaveHeatingError {}

/// Exact finite-volume wave load bound to its original gas and collision laws.
pub struct MagnetizedIonizedGasHeating<'w, 'a> {
    absorption: MagnetizedIonizedGasSlabAbsorption<'w, 'a>,
    momentum_sources: Vec<ElectronCrossSectionV1>,
    maximum_target_thermal_speed_ratio: Option<f64>,
}
impl<'w, 'a> MagnetizedIonizedGasHeating<'w, 'a> {
    /// The actual geometric slab volume must equal the retained parcel volume,
    /// admitting only floating-point roundoff. No hidden area or heat rescaling
    /// is performed. The caller declares that this control parcel mixes to
    /// uniform Th/Te between intervals despite its spatially varying absorption.
    ///
    /// Present electrons require the SAME neutral momentum cross sections
    /// already retained by their transport source. An effective scalar drag
    /// card cannot determine target-resolved thermal exchange and is refused;
    /// no independently invented thermal collision frequency is accepted.
    pub fn try_new(
        absorption: MagnetizedIonizedGasSlabAbsorption<'w, 'a>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IonizedWaveHeatingError> {
        poll(&mut cancelled)?;
        let source = absorption.wave().interface().transport().source();
        let ratio = absorption.volume_m3() / source.state().volume().value();
        if !ratio.is_finite() || (ratio - 1.0).abs() > 64.0 * f64::EPSILON {
            return Err(IonizedWaveHeatingError::Invalid(
                "wave absorption volume must equal the original finite gas parcel",
            ));
        }
        let mut momentum_sources = Vec::new();
        let mut maximum_target_thermal_speed_ratio = None;
        for carrier in source.carriers() {
            poll(&mut cancelled)?;
            if !carrier.is_electron() {
                continue;
            }
            match carrier.relaxation().source() {
                CarrierMomentumSource::ElectronNeutral {
                    contributions,
                    maximum_target_thermal_speed_ratio: bound,
                    ..
                } => {
                    momentum_sources.extend(
                        contributions
                            .iter()
                            .map(|c| c.evaluation().source().clone()),
                    );
                    maximum_target_thermal_speed_ratio = Some(*bound);
                }
                CarrierMomentumSource::Material { .. } => {
                    return Err(IonizedWaveHeatingError::Invalid(
                        "electron thermal exchange requires the original target-resolved momentum cross sections",
                    ));
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            absorption,
            momentum_sources,
            maximum_target_thermal_speed_ratio,
        })
    }

    pub const fn absorption(&self) -> &MagnetizedIonizedGasSlabAbsorption<'w, 'a> {
        &self.absorption
    }
    pub fn source(&self) -> &'a IonizedGasConductivity {
        self.absorption.wave().interface().transport().source()
    }
    pub fn momentum_sources(&self) -> &[ElectronCrossSectionV1] {
        &self.momentum_sources
    }

    /// Advance the existing finite gas energy owner at the exact signed wave
    /// powers held throughout this positive duration. Its implicit adaptive
    /// exchange, original energy differences, source-domain admission, error
    /// controls and separate physical/numerical entropy accounts remain intact.
    /// The original Coulomb law, when present, supplies zero-drift exchange;
    /// its signed drift correction is ALREADY in the wave's two pool powers.
    ///
    /// The caller explicitly bounds each pool's absolute temperature increment.
    /// This is a held-source approximation limit, not an error estimate or an
    /// automatically converged feedback step. Carrier/background refresh checks
    /// the endpoint only. Coherent kinetic, field and gyroscopic reactive terms
    /// are never deposited, and changing the next wave's stored energy would
    /// require its own drive/work account.
    ///
    /// Supply the original immutable momentum cards and bound-charge card.
    /// Publication is atomic across gas evolution, increment admission and all
    /// source refreshes; the old preparation remains usable after refusal.
    pub fn advance(
        &self,
        duration: Time,
        controls: ElectronHeavyEnergySolve,
        maximum_temperature_increment_k: f64,
        original_momentum_cards: &[&MaterialCard],
        original_background_card: &MaterialCard,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MagnetizedIonizedGasHeatStep<'_, 'w, 'a>, IonizedWaveHeatingError> {
        use IonizedWaveHeatingError::Invalid;
        poll(&mut cancelled)?;
        if !duration.value().is_finite()
            || duration.value() <= 0.0
            || !maximum_temperature_increment_k.is_finite()
            || maximum_temperature_increment_k <= 0.0
        {
            return Err(Invalid(
                "finite positive duration and explicit temperature-increment limit required",
            ));
        }
        let source = self.source();
        if controls.maximum_charged_mole_fraction > source.limits().maximum_charged_mole_fraction
            || self
                .maximum_target_thermal_speed_ratio
                .is_some_and(|bound| controls.maximum_target_thermal_speed_ratio > bound)
        {
            return Err(Invalid(
                "thermal controls cannot widen the original transport population or target-speed domain",
            ));
        }
        let background = self.absorption.wave().interface().background();
        if original_background_card.content_hash() != background.card_identity() {
            return Err(Invalid(
                "bound-charge refresh requires its original immutable material card",
            ));
        }
        let electron_power = self.absorption.electron_power_w();
        let heavy_power = self.absorption.heavy_power_w();
        let absorbed_energy_j = self.absorption.absorbed_power_w() * duration.value();
        for power in [
            electron_power,
            heavy_power,
            self.absorption.absorbed_power_w(),
        ] {
            let energy = power * duration.value();
            if !energy.is_finite() || (power != 0.0 && energy == 0.0) {
                return Err(Invalid("physical wave heating energy is unrepresentable"));
            }
        }
        let energy = match source.coulomb_model() {
            Some(coulomb) => advance_frozen_energy_with_coulomb(
                source.state(),
                &self.momentum_sources,
                coulomb,
                Power::new(electron_power),
                Power::new(heavy_power),
                duration,
                controls,
                &mut cancelled,
            ),
            None => advance_frozen_energy(
                source.state(),
                &self.momentum_sources,
                Power::new(electron_power),
                Power::new(heavy_power),
                duration,
                controls,
                &mut cancelled,
            ),
        }
        .map_err(|error| match error {
            ElectronHeavyEnergyError::Cancelled => IonizedWaveHeatingError::Cancelled,
            error => IonizedWaveHeatingError::Gas(error),
        })?;
        let next = energy.state();
        let heavy_k =
            (next.heavy_temperature().value() - source.state().heavy_temperature().value()).abs();
        let electron_k = (next.electron_temperature().value()
            - source.state().electron_temperature().value())
        .abs();
        if heavy_k.max(electron_k) > maximum_temperature_increment_k {
            return Err(IonizedWaveHeatingError::TemperatureIncrement {
                heavy_k,
                electron_k,
                maximum_k: maximum_temperature_increment_k,
            });
        }
        let refreshed_transport = source
            .at_frozen_state(next, original_momentum_cards, &mut cancelled)
            .map_err(|error| match error {
                IonizedGasRefreshError::Cancelled => IonizedWaveHeatingError::Cancelled,
                error => IonizedWaveHeatingError::Transport(error),
            })?;
        let mut query = background.query();
        for (axis, _) in background.query_point() {
            let value = match axis.as_str() {
                "T" | "T_h" => next.heavy_temperature().value(),
                "T_e" => next.electron_temperature().value(),
                "p" => next.pressure().value(),
                _ => return Err(Invalid("unmapped original plasma background condition")),
            };
            query = query
                .with(axis.clone(), value)
                .map_err(|_| Invalid("refreshed bound-charge coordinate is inadmissible"))?;
        }
        let refreshed_background = background
            .at_query_point(original_background_card, &query)
            .map_err(IonizedWaveHeatingError::Background)?;
        poll(&mut cancelled)?;
        Ok(MagnetizedIonizedGasHeatStep {
            heating: self,
            energy,
            refreshed_transport,
            refreshed_background,
            absorbed_energy_j,
            maximum_temperature_increment_k,
        })
    }
}

/// Original gas energy result and next material sources for a new field solve.
pub struct MagnetizedIonizedGasHeatStep<'h, 'w, 'a> {
    heating: &'h MagnetizedIonizedGasHeating<'w, 'a>,
    energy: ElectronHeavyEnergyResult,
    refreshed_transport: IonizedGasConductivity,
    refreshed_background: ResolvedMaterialStatePoint,
    absorbed_energy_j: f64,
    maximum_temperature_increment_k: f64,
}
impl<'h, 'w, 'a> MagnetizedIonizedGasHeatStep<'h, 'w, 'a> {
    pub const fn heating(&self) -> &'h MagnetizedIonizedGasHeating<'w, 'a> {
        self.heating
    }
    /// Original adaptive caloric/exchange result, including energy and entropy.
    pub const fn energy(&self) -> &ElectronHeavyEnergyResult {
        &self.energy
    }
    /// Same source laws reevaluated at the actual endpoint Th/Te and pressure.
    pub const fn refreshed_transport(&self) -> &IonizedGasConductivity {
        &self.refreshed_transport
    }
    /// Original bound-charge card/claims at the actual endpoint conditions.
    pub const fn refreshed_background(&self) -> &ResolvedMaterialStatePoint {
        &self.refreshed_background
    }
    /// Net original friction energy, credited once to the two native pools.
    pub const fn absorbed_energy_j(&self) -> f64 {
        self.absorbed_energy_j
    }
    pub const fn maximum_temperature_increment_k(&self) -> f64 {
        self.maximum_temperature_increment_k
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), IonizedWaveHeatingError> {
    if cancelled() {
        Err(IonizedWaveHeatingError::Cancelled)
    } else {
        Ok(())
    }
}
