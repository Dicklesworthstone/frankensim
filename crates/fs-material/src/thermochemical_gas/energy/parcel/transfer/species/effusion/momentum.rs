//! Half-range Maxwellian momentum from the SAME source-bearing molecular beams.
//!
//! An uncollimated beam's mass-weighted outward normal speed is
//! `sqrt(pi*R*T/(2*M_i))`; its normal momentum rate is `A*tau_i*p_i/2`.
//! This is neither the mean thermal speed nor a directed nozzle speed inferred
//! from the beam's complete energy. Transverse motion, velocity variance and
//! internal molecular modes retain their original energy account.
//!
//! Scalar mass transmission alone does not specify this momentum moment.
//! [`EffusiveMomentumInterface`] additionally declares angle-independent,
//! nondeflecting transmission. Wall-scattered or collimated pores need their own
//! angular law. The thin-aperture limit and the distinct mass/force factors are
//! described by Jacobson and Kowalski, NASA/TM-2020-220495, sections 1.1 and 7:
//! <https://ntrs.nasa.gov/api/citations/20200002141/downloads/20200002141.pdf>.
//!
//! Moments are in the reservoirs' common local rest frame. For a moving tank,
//! a mechanical owner must transport changing mass, recoil and the Galilean
//! transformation of the complete molecular energy together. No continuum
//! nozzle pressure-thrust term is added. For reciprocal beams, transmitted
//! normal momentum stress SUMS both nonnegative moments; equal mass rates do
//! not imply zero stress. Forces on separate walls or reservoir supports require
//! their geometry and momentum balances, not subtraction of these two moments.

use super::{
    EffusiveBeam, EffusiveFlux, EffusiveInterface, GasEnergyError, MAX_EFFUSIVE_SPECIES, SpeciesId,
    Sum, invalid, poll, positive, positive_product_ratio, vacuum::VacuumEffusiveFlux,
    valid_provenance,
};
use fs_math::det;
use std::sync::Arc;

/// Additional authored angular qualification bound to one original interface.
///
/// Every selected species transmission must be independent of incidence angle
/// and preserve the transmitted molecular velocity, throughout the original
/// gas/composition/T/p support. This declaration is distinct from the existing
/// energy-independent mass transmission and does not infer empirical accuracy.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusiveMomentumInterface {
    interface: EffusiveInterface,
    provenance: String,
}

impl EffusiveMomentumInterface {
    pub fn try_new(
        interface: EffusiveInterface,
        angular_provenance: impl Into<String>,
    ) -> Result<Self, GasEnergyError> {
        let provenance = angular_provenance.into();
        if !valid_provenance(&provenance) {
            return Err(invalid(
                "effusive momentum needs bounded explicit angular qualification provenance",
            ));
        }
        Ok(Self {
            interface,
            provenance,
        })
    }

    pub const fn interface(&self) -> &EffusiveInterface {
        &self.interface
    }
    pub fn provenance(&self) -> &str {
        &self.provenance
    }

    fn admit_interface(&self, original: &EffusiveInterface) -> Result<(), GasEnergyError> {
        // A clone retains these owners. A newly authored lookalike interface
        // cannot borrow an angular qualification made for another interface.
        if !Arc::ptr_eq(&self.interface.gas_catalog, &original.gas_catalog)
            || !Arc::ptr_eq(&self.interface.transmissions, &original.transmissions)
        {
            return Err(invalid(
                "effusive momentum lost its original qualified interface",
            ));
        }
        Ok(())
    }

    fn admit_species(
        &self,
        species: &SpeciesId,
        previous: Option<&SpeciesId>,
    ) -> Result<(), GasEnergyError> {
        if previous.is_some_and(|id| id >= species)
            || self
                .interface
                .transmissions()
                .binary_search_by(|row| row.species().cmp(species))
                .is_err()
        {
            return Err(invalid(
                "effusive momentum needs canonical original species beams",
            ));
        }
        Ok(())
    }

    /// Enrich already sampled vacuum beams; no source or flow is evaluated again.
    /// The result borrows the original flux and qualification without cloning
    /// material records. Finite consumers use each actual accepted species mass
    /// with [`EffusiveBeamMomentum::impulse_for_mass_kg`].
    pub fn vacuum_flux<'a>(
        &'a self,
        flux: &'a VacuumEffusiveFlux,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VacuumEffusiveMomentumFlux<'a>, GasEnergyError> {
        poll(&mut cancelled)?;
        self.admit_interface(&flux.interface)?;
        if flux.species.len() > MAX_EFFUSIVE_SPECIES {
            return Err(invalid(
                "effusive momentum exceeds the original species work bound",
            ));
        }
        let mut species = Vec::with_capacity(flux.species.len());
        let mut normal = Sum::default();
        let mut previous = None;
        for row in &flux.species {
            poll(&mut cancelled)?;
            self.admit_species(&row.species, previous)?;
            previous = Some(&row.species);
            let outward = momentum(&row.outward)?;
            normal.add(outward.normal_momentum_rate_n)?;
            species.push(VacuumEffusiveSpeciesMomentum {
                species: &row.species,
                outward,
            });
        }
        poll(&mut cancelled)?;
        Ok(VacuumEffusiveMomentumFlux {
            qualification: self,
            flux,
            species,
            normal_momentum_rate_n: normal.value,
        })
    }

    /// Preserve both original beam moments and their transmitted normal stress.
    /// A positive right-to-left moment describes leftward molecular momentum;
    /// its crossing direction is also reversed, so its stress contribution adds.
    pub fn reciprocal_flux<'a>(
        &'a self,
        flux: &'a EffusiveFlux,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EffusiveMomentumFlux<'a>, GasEnergyError> {
        poll(&mut cancelled)?;
        self.admit_interface(&flux.interface)?;
        if flux.species.len() > MAX_EFFUSIVE_SPECIES {
            return Err(invalid(
                "effusive momentum exceeds the original species work bound",
            ));
        }
        let mut species = Vec::with_capacity(flux.species.len());
        let mut stress = Sum::default();
        let mut previous = None;
        for row in &flux.species {
            poll(&mut cancelled)?;
            self.admit_species(&row.species, previous)?;
            previous = Some(&row.species);
            let left_to_right = momentum(&row.left_to_right)?;
            let right_to_left = momentum(&row.right_to_left)?;
            stress.add(left_to_right.normal_momentum_rate_n)?;
            stress.add(right_to_left.normal_momentum_rate_n)?;
            species.push(EffusiveSpeciesMomentum {
                species: &row.species,
                left_to_right,
                right_to_left,
            });
        }
        poll(&mut cancelled)?;
        Ok(EffusiveMomentumFlux {
            qualification: self,
            flux,
            species,
            normal_momentum_stress_n: stress.value,
        })
    }
}

/// A momentum moment borrowing its original source-resolved molecular beam.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffusiveBeamMomentum<'a> {
    beam: &'a EffusiveBeam,
    mean_normal_speed_m_s: Option<f64>,
    normal_momentum_rate_n: f64,
}
impl EffusiveBeamMomentum<'_> {
    pub const fn beam(&self) -> &EffusiveBeam {
        self.beam
    }
    /// Mass-flux-weighted outward speed; inactive beams have no speed query.
    pub const fn mean_normal_speed_m_s(&self) -> Option<f64> {
        self.mean_normal_speed_m_s
    }
    pub const fn normal_momentum_rate_n(&self) -> f64 {
        self.normal_momentum_rate_n
    }

    /// Integrate momentum from the mass ACTUALLY accepted by the inventory owner.
    /// The mass is a nonnegative outward amount, not a net bidirectional transfer.
    /// This avoids independently rounding force-times-duration after dm is fixed.
    pub fn impulse_for_mass_kg(&self, mass_kg: f64) -> Result<f64, GasEnergyError> {
        if !mass_kg.is_finite() || mass_kg < 0.0 {
            return Err(invalid(
                "effusive impulse needs finite nonnegative accepted mass",
            ));
        }
        if mass_kg == 0.0 {
            return Ok(0.0);
        }
        let speed = self.mean_normal_speed_m_s.ok_or(invalid(
            "inactive effusive beam cannot carry accepted mass or impulse",
        ))?;
        positive_product_ratio(&[mass_kg, speed], &[])
    }
}

/// Two nonnegative moments in canonical original species order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffusiveSpeciesMomentum<'a> {
    pub species: &'a SpeciesId,
    pub left_to_right: EffusiveBeamMomentum<'a>,
    pub right_to_left: EffusiveBeamMomentum<'a>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct EffusiveMomentumFlux<'a> {
    qualification: &'a EffusiveMomentumInterface,
    flux: &'a EffusiveFlux,
    species: Vec<EffusiveSpeciesMomentum<'a>>,
    normal_momentum_stress_n: f64,
}
impl<'a> EffusiveMomentumFlux<'a> {
    pub const fn qualification(&self) -> &EffusiveMomentumInterface {
        self.qualification
    }
    pub const fn flux(&self) -> &EffusiveFlux {
        self.flux
    }
    pub fn species(&self) -> &[EffusiveSpeciesMomentum<'a>] {
        &self.species
    }
    /// Positive transmitted normal momentum stress integrated over area [N].
    /// This equals A*sum(tau_i*(p_i,left+p_i,right))/2, not net mass times speed.
    pub const fn normal_momentum_stress_n(&self) -> f64 {
        self.normal_momentum_stress_n
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VacuumEffusiveSpeciesMomentum<'a> {
    pub species: &'a SpeciesId,
    pub outward: EffusiveBeamMomentum<'a>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumEffusiveMomentumFlux<'a> {
    qualification: &'a EffusiveMomentumInterface,
    flux: &'a VacuumEffusiveFlux,
    species: Vec<VacuumEffusiveSpeciesMomentum<'a>>,
    normal_momentum_rate_n: f64,
}
impl<'a> VacuumEffusiveMomentumFlux<'a> {
    pub const fn qualification(&self) -> &EffusiveMomentumInterface {
        self.qualification
    }
    pub const fn flux(&self) -> &VacuumEffusiveFlux {
        self.flux
    }
    pub fn species(&self) -> &[VacuumEffusiveSpeciesMomentum<'a>] {
        &self.species
    }
    /// Outward momentum emitted per time in the local gas rest frame [N].
    pub const fn normal_momentum_rate_n(&self) -> f64 {
        self.normal_momentum_rate_n
    }
}

fn momentum(beam: &EffusiveBeam) -> Result<EffusiveBeamMomentum<'_>, GasEnergyError> {
    if *beam == EffusiveBeam::default() {
        return Ok(EffusiveBeamMomentum {
            beam,
            mean_normal_speed_m_s: None,
            normal_momentum_rate_n: 0.0,
        });
    }
    let correction = beam.additional_specific_energy_j_kg.ok_or(invalid(
        "active effusive momentum beam lost its original molecular energy moment",
    ))?;
    if !positive(beam.mass_rate_kg_s)
        || !positive(beam.molar_rate_mol_s)
        || !correction.is_finite()
        || correction >= 0.0
    {
        return Err(invalid(
            "effusive momentum needs a positive original molecular beam",
        ));
    }
    // The original correction is -R*T/(2*M), kept outside the formation datum.
    // Taking square roots BEFORE multiplying avoids overflowing pi*(-correction)
    // even when the final normal speed and momentum are representable.
    let speed = positive_product_ratio(
        &[det::sqrt(core::f64::consts::PI), det::sqrt(-correction)],
        &[],
    )?;
    let normal = positive_product_ratio(&[beam.mass_rate_kg_s, speed], &[])?;
    Ok(EffusiveBeamMomentum {
        beam,
        mean_normal_speed_m_s: Some(speed),
        normal_momentum_rate_n: normal,
    })
}

#[cfg(test)]
mod tests;
