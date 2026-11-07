//! Refresh original drag sources after thermal or population evolution.

use super::{
    CarrierMomentumRelaxation, CarrierMomentumSource, IonizedGasConductivity, IonizedGasError,
};
use crate::state_point::MaterialStatePointError;
use fs_matdb::MaterialCard;
use fs_thermochem::coulomb::{CoulombErrorV1, CoulombModelV1};
use fs_thermochem::electron_collision::ElectronCollisionErrorV1;
use fs_thermochem::equilibrium::composition::state::FrozenIonizedGasStateV1;

/// Original-source refusal while refreshing the actual gas transport state.
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
        write!(f, "ionized-gas source refresh: {self:?}")
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IonizedGasRefreshError> {
        self.refresh_state(next, original_cards, true, cancelled)
    }

    /// Refresh the SAME original transport laws after finite-rate chemistry
    /// changes carrier and neutral populations in a rigid parcel. The accepted
    /// state supplies its actual amounts, Th, Te and pressure; no equilibrium
    /// reset, extra reaction heat or cached collision frequency is introduced.
    ///
    /// Original thermodynamics, volume and the set of present species must
    /// agree. An appearing carrier or neutral needs a separately admitted
    /// momentum source; an exhausted source is not silently dropped. All
    /// original numerical controls and weak-ionization/transport limits remain
    /// binding. Neutral drag uses the NEW target density and Maxwellian Te.
    /// The Coulomb source, when present, is rebound to the new populations with
    /// exactly the original screening, short-range cutoff and physical limits.
    ///
    /// Material-backed rates preserve their original immutable cards, complete
    /// property requirements and state schema, just as `at_frozen_state` does.
    /// The source must cover the actual endpoint; this query neither integrates
    /// chemistry nor establishes validity between the two states. Failure or
    /// cancellation leaves the original transport and gas snapshots unchanged.
    pub fn at_reacting_state(
        &self,
        next: &FrozenIonizedGasStateV1,
        original_cards: &[&MaterialCard],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IonizedGasRefreshError> {
        self.refresh_state(next, original_cards, false, cancelled)
    }

    fn refresh_state(
        &self,
        next: &FrozenIonizedGasStateV1,
        original_cards: &[&MaterialCard],
        frozen_populations: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IonizedGasRefreshError> {
        use IonizedGasRefreshError::Invalid;
        poll(&mut cancelled)?;
        if next.thermodynamics() != self.state().thermodynamics()
            || next.volume() != self.state().volume()
        {
            return Err(Invalid(
                "refresh requires the original thermodynamics and volume",
            ));
        }
        if frozen_populations && next.amounts() != self.state().amounts() {
            return Err(Invalid("frozen-state refresh requires unchanged amounts"));
        }
        if next
            .amounts()
            .iter()
            .zip(self.state().amounts())
            .any(|(new, old)| (new.value() > 0.0) != (old.value() > 0.0))
        {
            return Err(Invalid(
                "population refresh requires the original present-species source support",
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
                CarrierMomentumSource::Material { material, .. } => {
                    let i = identities
                        .iter()
                        .position(|&hash| hash == material.card_identity())
                        .ok_or(Invalid("missing original immutable momentum card"))?;
                    used[i] = true;
                    refresh_material_source(carrier.relaxation(), original_cards[i], next)?
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
            let model = if coulomb.is_bound_to_populations(next) {
                coulomb.clone()
            } else {
                CoulombModelV1::new_with_cutoff(
                    next,
                    coulomb.screening(),
                    coulomb.cutoff(),
                    coulomb.limits(),
                    &mut cancelled,
                )
                .map_err(|error| match error {
                    CoulombErrorV1::Cancelled => IonizedGasRefreshError::Cancelled,
                    error => transport_error(IonizedGasError::Coulomb(error)),
                })?
            };
            refreshed = refreshed
                .with_coulomb_drag(model, &mut cancelled)
                .map_err(transport_error)?;
        }
        poll(&mut cancelled)?;
        Ok(refreshed)
    }
}

/// Requery one original material law, including all retained requirements.
pub(super) fn refresh_material_source(
    template: &CarrierMomentumRelaxation,
    card: &MaterialCard,
    next: &FrozenIonizedGasStateV1,
) -> Result<CarrierMomentumRelaxation, IonizedGasRefreshError> {
    use IonizedGasRefreshError::{Invalid, Material};
    let CarrierMomentumSource::Material { material, property } = template.source() else {
        return Err(Invalid("material momentum template required"));
    };
    if material.card_identity() != card.content_hash() {
        return Err(Invalid("missing original immutable momentum card"));
    }
    let t = material.state_coordinate("T").is_some();
    let th = material.state_coordinate("T_h").is_some();
    let te = material.state_coordinate("T_e").is_some();
    let mut query = material.query();
    let mut set = |axis, value| -> Result<(), IonizedGasRefreshError> {
        query = query
            .clone()
            .with(axis, value)
            .map_err(|_| Invalid("refreshed momentum coordinate is inadmissible"))?;
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
    let material = material.at_query_point(card, &query).map_err(Material)?;
    if t {
        CarrierMomentumRelaxation::new(template.species().clone(), material, property.clone())
    } else {
        CarrierMomentumRelaxation::new_two_temperature(
            template.species().clone(),
            material,
            property.clone(),
        )
    }
    .map_err(transport_error)
}

pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), IonizedGasRefreshError> {
    if cancelled() {
        Err(IonizedGasRefreshError::Cancelled)
    } else {
        Ok(())
    }
}
pub(super) fn transport_error(error: IonizedGasError) -> IonizedGasRefreshError {
    match error {
        IonizedGasError::Cancelled => IonizedGasRefreshError::Cancelled,
        error => IonizedGasRefreshError::Transport(error),
    }
}
