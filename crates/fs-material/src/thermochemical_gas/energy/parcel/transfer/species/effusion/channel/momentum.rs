//! Directional molecular momentum through a finite diffuse circular channel.
//!
//! A channel changes the exit angular distribution. Transmitted mass times a
//! thin-aperture normal speed therefore does not supply the exit momentum.
//! The Clausing owner resolves incident-normalized axial-cosine moments for
//! direct, wall-transmitted and returned paths on the SAME mass geometry.
//! Flux-weighted mean speed is 3/2 times the thin-aperture mean normal speed.
//!
//! Each diffuse collision either resets translational and internal energy to
//! the wall Maxwellian, or preserves incident speed and internal state. The
//! authored reset probability is independent of angle, speed and prior hits.
//! Joint angular/collision-history moments separate molecules retaining their
//! original gas speed from those reset by any earlier wall encounter. A mean
//! energy temperature or a single-hit interpolation is not their mean speed.
//! Walls are stationary, axisymmetric and isothermal along the bore;
//! transverse means cancel by symmetry.
//!
//! These are molecular boundary moments, excluding other container tractions.
//! No continuum nozzle pressure thrust is added. One incident beam gives wall
//! impulse P_in + P_return - P_transmitted. A vessel containing both its source
//! gas and this wall therefore recoils against ONLY the escaping momentum.
//! A mechanics owner must transport changing mass and complete molecular energy
//! together when using these local-rest-frame data.

use super::super::momentum::{EffusiveBeamMomentum, EffusiveMomentumInterface};
use super::{
    ChannelAccommodation, ClausingControls, ClosedGasParcel, DiffuseChannelFlux,
    DiffuseMolecularChannel, GasEnergyError, MAX_EFFUSIVE_SPECIES, MolecularAxialAccommodation,
    MolecularAxialMoments, SpeciesId, invalid, poll, positive, positive_product_ratio,
    vacuum::VacuumDiffuseChannelFlux,
};
use crate::thermochemical_gas::energy::parcel::transfer::{
    snapshot::MAX_SNAPSHOT_TRANSFERS,
    species::network::transport::effusion::channel::{
        DiffuseChannelLinkStep, DiffuseChannelNetworkStep, thermal::FiniteChannelNetworkStep,
        vacuum::VacuumChannelOutletStep,
    },
};
use fs_math::{dd::Dd, det};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::sync::Arc;

/// Original channel with an explicit speed-distribution qualification.
/// Construction resolves geometry once; later samples reuse all moments.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffuseChannelMomentumInterface {
    channel: DiffuseMolecularChannel,
    mouth: EffusiveMomentumInterface,
    moments: MolecularAxialMoments,
    /// Same canonical order as the original immutable gas catalog.
    accommodation: Arc<[MolecularAxialAccommodation]>,
}

/// Each beam uses its own outward axis, including a right-to-left beam.
/// Only the wall component may be negative.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ChannelBeamForces {
    pub incident_n: f64,
    pub transmitted_n: f64,
    pub returned_n: f64,
    pub wall_n: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ChannelSpeciesForces {
    pub species: SpeciesId,
    pub from_left: ChannelBeamForces,
    pub from_right: ChannelBeamForces,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ChannelMomentumRates<F> {
    pub flux: F,
    pub species: Vec<ChannelSpeciesForces>,
    /// Momentum delivered outside each mouth, positive left-to-right.
    /// For vacuum, the right value is escaped momentum per time, not force
    /// on an invented receiver. The left value enters the actual source gas.
    pub exterior_forces_n: [f64; 2],
    pub wall_force_n: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ChannelBeamImpulse {
    pub incident_n_s: f64,
    pub transmitted_n_s: f64,
    pub returned_n_s: f64,
    pub wall_n_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ChannelSpeciesImpulse {
    pub species: SpeciesId,
    pub from_left: ChannelBeamImpulse,
    pub from_right: ChannelBeamImpulse,
}

/// Momentum from masses accepted by the original channel owner.
/// Its gas/material transaction retains the complete energy account.
#[derive(Clone, Debug, PartialEq)]
pub struct ChannelMomentumStep {
    pub id: u64,
    pub species: Vec<ChannelSpeciesImpulse>,
    /// Same exterior/sign convention as ChannelMomentumRates [N s].
    pub exterior_impulses_n_s: [f64; 2],
    pub wall_impulse_n_s: f64,
}

fn accommodation_coefficient(
    channel: &DiffuseMolecularChannel,
    id: &SpeciesId,
) -> Result<f64, GasEnergyError> {
    let coefficient = match &channel.accommodation {
        None => 1.0,
        Some(ChannelAccommodation::Common(model, _)) => model.coefficient(),
        Some(ChannelAccommodation::Species(rows)) => {
            let index = rows
                .binary_search_by(|row| row.qualification.species.cmp(id))
                .map_err(|_| invalid("channel momentum lost a species accommodation"))?;
            rows[index].accommodation.coefficient()
        }
    };
    Ok(coefficient)
}

fn forces(a: [f64; 4]) -> ChannelBeamForces {
    ChannelBeamForces {
        incident_n: a[0],
        transmitted_n: a[1],
        returned_n: a[2],
        wall_n: a[3],
    }
}
fn impulse(a: [f64; 4]) -> ChannelBeamImpulse {
    ChannelBeamImpulse {
        incident_n_s: a[0],
        transmitted_n_s: a[1],
        returned_n_s: a[2],
        wall_n_s: a[3],
    }
}

// Native reports are public data: reject omitted/duplicate/foreign inventory
// instead of silently ignoring an amount while enriching its momentum.
fn admit_receipt_rows<'a>(
    rows: impl IntoIterator<Item = (&'a SpeciesId, [f64; 2])>,
    original: impl Iterator<Item = &'a SpeciesId> + Clone,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), GasEnergyError> {
    let mut previous = None;
    for (index, (id, amounts)) in rows.into_iter().enumerate() {
        poll(cancelled)?;
        if index >= MAX_EFFUSIVE_SPECIES
            || previous.is_some_and(|old| old >= id)
            || !original.clone().any(|known| known == id)
            || amounts.into_iter().any(|m| !m.is_finite() || m < 0.0)
            || amounts == [0.0; 2]
        {
            return Err(invalid(
                "channel impulse needs canonical original nonzero species amounts",
            ));
        }
        previous = Some(id);
    }
    Ok(())
}

#[derive(Default)]
struct Balance {
    exterior: [Dd; 2],
    wall: Dd,
}
impl Balance {
    fn add(&mut self, beams: [[f64; 4]; 2]) -> Result<(), GasEnergyError> {
        let [left, right] = beams.map(|row| row.map(Dd::from_f64));
        self.exterior[0] = self.exterior[0] - left[0] - left[2] - right[1];
        self.exterior[1] = self.exterior[1] + left[1] + right[0] + right[2];
        self.wall = self.wall + left[3] - right[3];
        if self
            .exterior
            .iter()
            .chain([&self.wall])
            .any(|v| !v.to_f64().is_finite())
        {
            return Err(invalid("channel momentum accumulation is unrepresentable"));
        }
        Ok(())
    }
}

impl DiffuseChannelMomentumInterface {
    /// Qualify diffuse scattering with a species-specific probability of a
    /// complete Maxwellian speed/internal-energy reset at each wall collision.
    /// The complementary event preserves speed and internal state. Arbitrary
    /// probabilities in [0,1] retain their full angular/collision history.
    /// A zero-length mouth has no wall encounters and needs no wall speed law.
    ///
    /// Use channel() with the finite transport owners: it contains the refined
    /// mass geometry and matching heat/angular history weights. Earlier approximate
    /// geometry cannot supply a receipt for this momentum qualification.
    pub fn try_new(
        mut channel: DiffuseMolecularChannel,
        controls: ClausingControls,
        speed_distribution_provenance: impl Into<String>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, GasEnergyError> {
        poll(&mut cancelled)?;
        let mut coefficients = Vec::new();
        for source in channel.mouth().gas_catalog().thermodynamics().components() {
            poll(&mut cancelled)?;
            coefficients.push(accommodation_coefficient(&channel, source.species())?);
        }
        let mouth = EffusiveMomentumInterface::try_new(
            channel.mouth().clone(),
            speed_distribution_provenance,
        )?;
        let (geometry, moments, accommodation) = channel
            .geometry
            .clone()
            .resolve_axial_accommodation_batch(&coefficients, controls, &mut cancelled)?;
        if accommodation.len() != coefficients.len() {
            return Err(invalid(
                "channel momentum lost an original accommodation response",
            ));
        }
        channel.geometry = geometry;
        // The joint solve supplies the heat law on exactly the same final
        // mass/angular grid. Preserve each authored surface qualification and
        // bind its returned native response without another geometry solve.
        match &mut channel.accommodation {
            None => {}
            Some(ChannelAccommodation::Common(model, _)) => {
                *model = accommodation
                    .first()
                    .ok_or(invalid("channel momentum needs an original gas species"))?
                    .accommodation()
                    .clone();
            }
            Some(ChannelAccommodation::Species(rows)) => {
                for (row, response) in Arc::make_mut(rows).iter_mut().zip(&accommodation) {
                    poll(&mut cancelled)?;
                    row.accommodation = response.accommodation().clone();
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            channel,
            mouth,
            moments,
            accommodation: accommodation.into(),
        })
    }

    pub const fn channel(&self) -> &DiffuseMolecularChannel {
        &self.channel
    }
    pub const fn moments(&self) -> &MolecularAxialMoments {
        &self.moments
    }
    /// Retained-gas and wall-reset angular moments for one original species.
    /// Their heat response is the one bound to channel() on the same geometry.
    pub fn axial_accommodation(
        &self,
        species: &SpeciesId,
    ) -> Result<&MolecularAxialAccommodation, GasEnergyError> {
        let index = self
            .channel
            .mouth()
            .gas_catalog()
            .thermodynamics()
            .components()
            .binary_search_by(|source| source.species().cmp(species))
            .map_err(|_| invalid("channel momentum needs an original gas species"))?;
        self.accommodation.get(index).ok_or(invalid(
            "channel momentum lost an original accommodation response",
        ))
    }
    pub fn provenance(&self) -> &str {
        self.mouth.provenance()
    }

    // Amounts are accepted transmitted/returned masses, or rates for a unit
    // interval. This one kernel owns both the force and impulse calculations.
    fn account(
        &self,
        species: &SpeciesId,
        beam: EffusiveBeamMomentum<'_>,
        wall_temperature_k: f64,
        transmitted: f64,
        returned: f64,
    ) -> Result<[f64; 4], GasEnergyError> {
        if ![transmitted, returned]
            .into_iter()
            .all(|m| m.is_finite() && m >= 0.0)
        {
            return Err(invalid(
                "channel momentum needs finite nonnegative original amounts",
            ));
        }
        if transmitted == 0.0 {
            return if returned == 0.0 {
                Ok([0.0; 4])
            } else {
                Err(invalid(
                    "returned molecules lost their original transmitted beam",
                ))
            };
        }
        let incident_speed = beam.mean_normal_speed_m_s().ok_or(invalid(
            "inactive molecular beam cannot carry accepted mass",
        ))?;
        if self.channel.geometry().length_m() == 0.0 {
            if returned != 0.0 {
                return Err(invalid("zero-length channel cannot return molecules"));
            }
            let p = positive_product_ratio(&[transmitted, incident_speed], &[])?;
            return Ok([p, p, 0.0, 0.0]);
        }
        let accommodation = self.axial_accommodation(species)?;
        let wall_speed = if accommodation.reset_transmitted_axial_cosine() > 0.0
            || accommodation.reset_reflected_axial_cosine() > 0.0
        {
            self.channel.admit_wall_temperature(wall_temperature_k)?;
            let mass = self.channel.original_species(species)?.molar_mass().value();
            positive_product_ratio(
                &[
                    det::sqrt(core::f64::consts::PI),
                    det::sqrt(R),
                    det::sqrt(wall_temperature_k),
                ],
                &[det::sqrt(2.0), det::sqrt(mass)],
            )?
        } else {
            incident_speed
        };
        let tau = self.channel.geometry().transmission_probability();
        let reflected = self.channel.geometry().reflected_probability();
        if (returned > 0.0) != (reflected > 0.0) {
            return Err(invalid(
                "channel momentum lost its original gross returned inventory",
            ));
        }
        let term = |speed, weight, probability| {
            if weight == 0.0 {
                Ok(0.0)
            } else {
                positive_product_ratio(&[1.5, speed, weight], &[probability])
            }
        };
        let transmitted_speed = (Dd::from_f64(term(
            incident_speed,
            accommodation.retained_transmitted_axial_cosine(),
            tau,
        )?) + Dd::from_f64(term(
            wall_speed,
            accommodation.reset_transmitted_axial_cosine(),
            tau,
        )?))
        .to_f64();
        let p_transmitted = positive_product_ratio(&[transmitted, transmitted_speed], &[])?;
        let p_returned = if returned == 0.0 {
            0.0
        } else {
            let returned_speed = (Dd::from_f64(term(
                incident_speed,
                accommodation.retained_reflected_axial_cosine(),
                reflected,
            )?) + Dd::from_f64(term(
                wall_speed,
                accommodation.reset_reflected_axial_cosine(),
                reflected,
            )?))
            .to_f64();
            positive_product_ratio(&[returned, returned_speed], &[])?
        };
        let p_incident = ((Dd::from_f64(transmitted) + Dd::from_f64(returned))
            * Dd::from_f64(incident_speed))
        .to_f64();
        let p_wall = (Dd::from_f64(p_incident) + Dd::from_f64(p_returned)
            - Dd::from_f64(p_transmitted))
        .to_f64();
        if !positive(p_incident) || !p_wall.is_finite() {
            return Err(invalid(
                "channel incident or wall momentum is unrepresentable",
            ));
        }
        Ok([p_incident, p_transmitted, p_returned, p_wall])
    }

    /// Sample original mass/heat rates and enrich their SAME source beams.
    /// The right exterior account is escaped momentum, with no vacuum state.
    pub fn sample_to_vacuum(
        &self,
        gas: &ClosedGasParcel,
        area_m2: f64,
        wall_temperature_k: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ChannelMomentumRates<VacuumDiffuseChannelFlux>, GasEnergyError> {
        let flux =
            self.channel
                .sample_to_vacuum(gas, area_m2, wall_temperature_k, &mut cancelled)?;
        let original = self.mouth.vacuum_flux(&flux.transmitted, &mut cancelled)?;
        let mut species = Vec::with_capacity(original.species().len());
        let mut balance = Balance::default();
        for row in original.species() {
            poll(&mut cancelled)?;
            let returned = flux
                .wall
                .reflected_species
                .iter()
                .find(|r| &r.species == row.species)
                .map_or(0.0, |r| r.mass_rate_kg_s);
            let beam = self.account(
                row.species,
                row.outward,
                wall_temperature_k,
                row.outward.beam().mass_rate_kg_s,
                returned,
            )?;
            balance.add([beam, [0.0; 4]])?;
            species.push(ChannelSpeciesForces {
                species: row.species.clone(),
                from_left: forces(beam),
                from_right: ChannelBeamForces::default(),
            });
        }
        poll(&mut cancelled)?;
        Ok(ChannelMomentumRates {
            flux,
            species,
            exterior_forces_n: balance.exterior.map(Dd::to_f64),
            wall_force_n: balance.wall.to_f64(),
        })
    }

    /// Preserve both reciprocal moments even when opposite mass rates cancel.
    /// Exterior and wall signs use one left-to-right channel axis.
    pub fn sample(
        &self,
        gases: [&ClosedGasParcel; 2],
        area_m2: f64,
        wall_temperature_k: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ChannelMomentumRates<DiffuseChannelFlux>, GasEnergyError> {
        let flux = self.channel.sample(
            gases[0],
            gases[1],
            area_m2,
            wall_temperature_k,
            &mut cancelled,
        )?;
        let original = self
            .mouth
            .reciprocal_flux(&flux.transmitted, &mut cancelled)?;
        let mut species = Vec::with_capacity(original.species().len());
        let mut balance = Balance::default();
        for row in original.species() {
            poll(&mut cancelled)?;
            let returned = flux
                .wall
                .reflected_species
                .iter()
                .find(|r| &r.species == row.species)
                .map_or([0.0; 2], |r| r.mass_rate_kg_s);
            let mut beams = [[0.0; 4]; 2];
            for (side, beam) in [row.left_to_right, row.right_to_left]
                .into_iter()
                .enumerate()
            {
                beams[side] = self.account(
                    row.species,
                    beam,
                    wall_temperature_k,
                    beam.beam().mass_rate_kg_s,
                    returned[side],
                )?;
            }
            balance.add(beams)?;
            species.push(ChannelSpeciesForces {
                species: row.species.clone(),
                from_left: forces(beams[0]),
                from_right: forces(beams[1]),
            });
        }
        poll(&mut cancelled)?;
        Ok(ChannelMomentumRates {
            flux,
            species,
            exterior_forces_n: balance.exterior.map(Dd::to_f64),
            wall_force_n: balance.wall.to_f64(),
        })
    }

    /// Enrich an accepted vacuum outlet from held/finite walls or adaptive
    /// histories. Retained escaped and returned masses are used exactly;
    /// there is no second source sample or rate*time rounding.
    pub fn vacuum_outlet(
        &self,
        outlet: &VacuumChannelOutletStep,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ChannelMomentumStep, GasEnergyError> {
        poll(&mut cancelled)?;
        if outlet.channel != self.channel {
            return Err(invalid(
                "channel receipt uses a different geometry or surface law",
            ));
        }
        let original = self
            .mouth
            .vacuum_flux(&outlet.transmitted, &mut cancelled)?;
        admit_receipt_rows(
            outlet
                .escaped_species
                .iter()
                .map(|r| (&r.species, [r.mass_kg, 0.0])),
            original.species().iter().map(|r| r.species),
            &mut cancelled,
        )?;
        admit_receipt_rows(
            outlet
                .wall
                .reflected_species
                .iter()
                .map(|r| (&r.species, [r.mass_kg, 0.0])),
            original.species().iter().map(|r| r.species),
            &mut cancelled,
        )?;
        let mut species = Vec::with_capacity(original.species().len());
        let mut balance = Balance::default();
        for row in original.species() {
            poll(&mut cancelled)?;
            let escaped = outlet
                .escaped_species
                .iter()
                .find(|r| &r.species == row.species)
                .map_or(0.0, |r| r.mass_kg);
            let returned = outlet
                .wall
                .reflected_species
                .iter()
                .find(|r| &r.species == row.species)
                .map_or(0.0, |r| r.mass_kg);
            if (escaped > 0.0) != (row.outward.beam().mass_rate_kg_s > 0.0) {
                return Err(invalid("channel impulse lost an accepted escaping species"));
            }
            let beam = self.account(
                row.species,
                row.outward,
                outlet.wall_temperature_k,
                escaped,
                returned,
            )?;
            balance.add([beam, [0.0; 4]])?;
            species.push(ChannelSpeciesImpulse {
                species: row.species.clone(),
                from_left: impulse(beam),
                from_right: ChannelBeamImpulse::default(),
            });
        }
        poll(&mut cancelled)?;
        Ok(ChannelMomentumStep {
            id: outlet.id,
            species,
            exterior_impulses_n_s: balance.exterior.map(Dd::to_f64),
            wall_impulse_n_s: balance.wall.to_f64(),
        })
    }

    /// Accepted reciprocal held-wall amounts use the original network's exact
    /// rate*duration expression. Gross returns come from the accepted receipt.
    pub fn reciprocal_step(
        &self,
        step: &DiffuseChannelNetworkStep,
        id: u64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ChannelMomentumStep, GasEnergyError> {
        self.reciprocal_links(&step.channels, step.duration_s, id, &mut cancelled)
    }

    /// The same calculation on finite material-wall steps, including accepted
    /// adaptive histories and the filled phase of receiver initialization.
    pub fn finite_wall_step(
        &self,
        step: &FiniteChannelNetworkStep,
        id: u64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ChannelMomentumStep, GasEnergyError> {
        self.reciprocal_links(&step.channels, step.duration_s, id, &mut cancelled)
    }

    fn reciprocal_links(
        &self,
        channels: &[DiffuseChannelLinkStep],
        duration_s: f64,
        id: u64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ChannelMomentumStep, GasEnergyError> {
        poll(cancelled)?;
        if !positive(duration_s) || channels.len() > MAX_SNAPSHOT_TRANSFERS {
            return Err(invalid(
                "channel momentum needs a bounded original physical step",
            ));
        }
        let mut matching = channels.iter().filter(|link| link.id == id);
        let link = matching
            .next()
            .ok_or(invalid("channel momentum lost its original link"))?;
        if matching.next().is_some() || link.channel != self.channel {
            return Err(invalid(
                "channel momentum link has different geometry or duplicate identity",
            ));
        }
        let original = self
            .mouth
            .reciprocal_flux(&link.transmitted, &mut *cancelled)?;
        admit_receipt_rows(
            link.wall
                .reflected_species
                .iter()
                .map(|r| (&r.species, r.mass_kg)),
            original.species().iter().map(|r| r.species),
            cancelled,
        )?;
        let mut species = Vec::with_capacity(original.species().len());
        let mut balance = Balance::default();
        for row in original.species() {
            poll(cancelled)?;
            let returned = link
                .wall
                .reflected_species
                .iter()
                .find(|r| &r.species == row.species)
                .map_or([0.0; 2], |r| r.mass_kg);
            let mut beams = [[0.0; 4]; 2];
            for (side, beam) in [row.left_to_right, row.right_to_left]
                .into_iter()
                .enumerate()
            {
                beams[side] = self.account(
                    row.species,
                    beam,
                    link.wall_temperature_k,
                    beam.beam().mass_rate_kg_s * duration_s,
                    returned[side],
                )?;
            }
            balance.add(beams)?;
            species.push(ChannelSpeciesImpulse {
                species: row.species.clone(),
                from_left: impulse(beams[0]),
                from_right: impulse(beams[1]),
            });
        }
        poll(cancelled)?;
        Ok(ChannelMomentumStep {
            id,
            species,
            exterior_impulses_n_s: balance.exterior.map(Dd::to_f64),
            wall_impulse_n_s: balance.wall.to_f64(),
        })
    }
}

#[cfg(test)]
mod tests;
