//! Refresh the original drag sources after heating a frozen gas parcel.

use super::{
    CarrierMomentumRelaxation, CarrierMomentumSource, IonizedGasConductivity, IonizedGasError,
};
use crate::state_point::MaterialStatePointError;
use fs_matdb::MaterialCard;
use fs_thermochem::electron_collision::ElectronCollisionErrorV1;
use fs_thermochem::equilibrium::composition::state::FrozenIonizedGasStateV1;

/// Original-source refusal while refreshing a frozen-population transport law.
#[derive(Clone, Debug, PartialEq)]
pub enum IonizedGasRefreshError {
    Invalid(&'static str),
    Material(MaterialStatePointError),
    Collision(ElectronCollisionErrorV1),
    Transport(IonizedGasError),
    Cancelled,
}
impl core::fmt::Display for IonizedGasRefreshError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ionized-gas frozen-state refresh: {self:?}")
    }
}
impl std::error::Error for IonizedGasRefreshError {}

impl IonizedGasConductivity {
    /// Re-evaluate every original momentum source at new heavy/electron
    /// temperatures and the resulting pressure of this SAME fixed-volume,
    /// frozen-population gas. The thermodynamic bank, amounts, volume, transport
    /// limits, cross-section controls and optional Coulomb law remain exact.
    ///
    /// Supply each distinct original material card once. Its original complete
    /// property requirements, selected claims, other coordinates and quantity
    /// conventions are preserved. A common-temperature source cannot be widened
    /// into a two-temperature source. Neutral cross sections are integrated at
    /// the new electron temperature; cached rates are never reused as constants.
    ///
    /// This is an endpoint evaluation, not a thermal integrator or a guarantee
    /// that a source covers the path between states. No LTE equilibrium receipt
    /// is attached to the new frozen state. Failure or cancellation leaves the
    /// original snapshot unchanged and publishes no partial refresh.
    pub fn at_frozen_state(
        &self,
        next: &FrozenIonizedGasStateV1,
        original_cards: &[&MaterialCard],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IonizedGasRefreshError> {
        use IonizedGasRefreshError::{Invalid, Material};
        poll(&mut cancelled)?;
        if next.thermodynamics() != self.state().thermodynamics()
            || next.amounts() != self.state().amounts()
            || next.volume() != self.state().volume()
        {
            return Err(Invalid(
                "refresh requires the original thermodynamics, amounts and volume",
            ));
        }
        if original_cards.len() > self.carriers().len() {
            return Err(Invalid(
                "momentum-card bank exceeds the original carrier count",
            ));
        }
        let identities: Vec<_> = original_cards
            .iter()
            .map(|card| card.content_hash())
            .collect();
        for (i, hash) in identities.iter().enumerate() {
            if identities[..i].contains(hash) {
                return Err(Invalid("supply each original momentum card exactly once"));
            }
        }
        let mut used = vec![false; identities.len()];
        let mut rates = Vec::with_capacity(self.carriers().len());
        for carrier in self.carriers() {
            poll(&mut cancelled)?;
            let rate = match carrier.relaxation().source() {
                CarrierMomentumSource::Material { material, property } => {
                    let i = identities
                        .iter()
                        .position(|&hash| hash == material.card_identity())
                        .ok_or(Invalid("missing original immutable momentum card"))?;
                    used[i] = true;
                    let t = material.state_coordinate("T").is_some();
                    let th = material.state_coordinate("T_h").is_some();
                    let te = material.state_coordinate("T_e").is_some();
                    let mut query = material.query();
                    let mut set = |axis, value| -> Result<(), IonizedGasRefreshError> {
                        query = query.clone().with(axis, value).map_err(|_| {
                            Invalid("refreshed momentum coordinate is inadmissible")
                        })?;
                        Ok(())
                    };
                    match (t, th, te) {
                        (true, false, false) => {
                            if next.heavy_temperature() != next.electron_temperature() {
                                return Err(Invalid(
                                    "common-temperature momentum source cannot admit unequal temperatures",
                                ));
                            }
                            set("T", next.heavy_temperature().value())?;
                        }
                        (false, true, true) => {
                            set("T_h", next.heavy_temperature().value())?;
                            set("T_e", next.electron_temperature().value())?;
                        }
                        _ => {
                            return Err(Invalid(
                                "momentum refresh requires an unambiguous original temperature schema",
                            ));
                        }
                    }
                    set("p", next.pressure().value())?;
                    let material = material
                        .at_query_point(original_cards[i], &query)
                        .map_err(Material)?;
                    if t {
                        CarrierMomentumRelaxation::new(
                            carrier.species().clone(),
                            material,
                            property.clone(),
                        )
                    } else {
                        CarrierMomentumRelaxation::new_two_temperature(
                            carrier.species().clone(),
                            material,
                            property.clone(),
                        )
                    }
                    .map_err(transport_error)?
                }
                CarrierMomentumSource::ElectronNeutral {
                    contributions,
                    maximum_target_thermal_speed_ratio,
                    ..
                } => {
                    let mut evaluations = Vec::with_capacity(contributions.len());
                    for contribution in contributions {
                        let old = contribution.evaluation();
                        let evaluation = old
                            .source()
                            .evaluate(
                                next.electron_temperature(),
                                old.weighting(),
                                old.controls(),
                                &mut cancelled,
                            )
                            .map_err(|error| match error {
                                ElectronCollisionErrorV1::Cancelled => {
                                    IonizedGasRefreshError::Cancelled
                                }
                                error => IonizedGasRefreshError::Collision(error),
                            })?;
                        evaluations.push(evaluation);
                    }
                    CarrierMomentumRelaxation::from_electron_neutral_cross_sections(
                        next,
                        evaluations,
                        *maximum_target_thermal_speed_ratio,
                        &mut cancelled,
                    )
                    .map_err(transport_error)?
                }
            };
            rates.push(rate);
        }
        if used.contains(&false) {
            return Err(Invalid("momentum-card bank contains an unused source"));
        }
        let mut refreshed = Self::from_state(next.clone(), rates, self.limits(), &mut cancelled)
            .map_err(transport_error)?;
        if let Some(coulomb) = self.coulomb_model() {
            refreshed = refreshed
                .with_coulomb_drag(coulomb.clone(), &mut cancelled)
                .map_err(transport_error)?;
        }
        poll(&mut cancelled)?;
        Ok(refreshed)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), IonizedGasRefreshError> {
    if cancelled() {
        Err(IonizedGasRefreshError::Cancelled)
    } else {
        Ok(())
    }
}
fn transport_error(error: IonizedGasError) -> IonizedGasRefreshError {
    match error {
        IonizedGasError::Cancelled => IonizedGasRefreshError::Cancelled,
        error => IonizedGasRefreshError::Transport(error),
    }
}
