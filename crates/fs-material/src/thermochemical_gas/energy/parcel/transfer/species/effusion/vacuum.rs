//! One-way molecular emission into an ideal absorbing vacuum.
//!
//! Vacuum specifies an absent incident beam, not an equilibrium gas with a
//! fabricated temperature, pressure floor or finite inventory. The emitting
//! gas still obeys the original thin-aperture, local-Maxwellian and source
//! qualification. Its species leave at their own molecular speeds and carry
//! `h_i(T) - R*T/(2*M_i)`, including the original formation-energy datum.
//! This is not a pump curve, an expanding jet or a finite receiver fill model.

use super::{
    ClosedGasParcel, EffusiveBeam, EffusiveInterface, GasEnergyError, SpeciesId, Sum, beam,
    invalid, poll,
};

/// An original species' outgoing molecular beam. The incident vacuum beam is zero.
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumEffusiveSpeciesFlux {
    pub species: SpeciesId,
    pub outward: EffusiveBeam,
}

/// Source-bound outward rates, before any finite inventory is changed.
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumEffusiveFlux {
    pub interface: EffusiveInterface,
    /// Canonical original species order; empty only for zero geometric area.
    pub species: Vec<VacuumEffusiveSpeciesFlux>,
    pub molar_rate_mol_s: f64,
    pub mass_rate_kg_s: f64,
    /// Complete outward molecular energy in the original species datum [W].
    /// This may be negative for a valid choice of formation-energy references.
    pub energy_rate_w: f64,
    /// Molecular correction relative to the outgoing partial enthalpies [W].
    pub additional_energy_rate_w: f64,
}

impl ClosedGasParcel {
    /// Sample effusion into a true vacuum through the original molecular law.
    ///
    /// No artificial receiver, back-pressure, temperature or dilute trace gas
    /// is introduced. The material qualification applies to the emitting gas;
    /// vacuum is not queried against a positive-pressure gas source domain.
    /// Zero exposed area is exactly inactive, including for an unrelated
    /// interface catalog. A selected zero transmission still admits its source.
    ///
    /// A finite update must debit every gross original species inventory and
    /// apply its molecular energy correction before the one caloric inversion.
    /// Use `network::transport::effusion::vacuum::step_candidate` to combine
    /// these outlets with the other native gas paths on one shared snapshot.
    pub fn effusive_flux_to_vacuum(
        &self,
        interface: &EffusiveInterface,
        area_m2: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VacuumEffusiveFlux, GasEnergyError> {
        poll(&mut cancelled)?;
        if !area_m2.is_finite() || area_m2 < 0.0 {
            return Err(invalid(
                "vacuum effusive area must be finite and nonnegative",
            ));
        }
        let mut result = VacuumEffusiveFlux {
            interface: interface.clone(),
            species: Vec::new(),
            molar_rate_mol_s: 0.0,
            mass_rate_kg_s: 0.0,
            energy_rate_w: 0.0,
            additional_energy_rate_w: 0.0,
        };
        if area_m2 == 0.0 {
            poll(&mut cancelled)?;
            return Ok(result);
        }
        interface.admit_parcel(self, &mut cancelled)?;
        let mixture = self.model().thermodynamics();
        result.species.reserve(mixture.components().len());
        let (mut moles, mut mass, mut energy, mut extra) = (
            Sum::default(),
            Sum::default(),
            Sum::default(),
            Sum::default(),
        );
        for (species, &fraction) in mixture
            .components()
            .iter()
            .zip(mixture.mole_composition().fractions())
        {
            poll(&mut cancelled)?;
            let index = interface
                .transmissions()
                .binary_search_by(|row| row.species().cmp(species.species()))
                .map_err(|_| invalid("vacuum interface lost an original gas species"))?;
            let outward = beam(
                self,
                species.species(),
                fraction,
                species.molar_mass().value(),
                interface.transmissions()[index].probability(),
                area_m2,
            )?;
            moles.add(outward.molar_rate_mol_s)?;
            mass.add(outward.mass_rate_kg_s)?;
            energy.add(outward.energy_rate_w)?;
            extra.add(outward.additional_energy_rate_w)?;
            result.species.push(VacuumEffusiveSpeciesFlux {
                species: species.species().clone(),
                outward,
            });
        }
        poll(&mut cancelled)?;
        result.molar_rate_mol_s = moles.value;
        result.mass_rate_kg_s = mass.value;
        result.energy_rate_w = energy.value;
        result.additional_energy_rate_w = extra.value;
        Ok(result)
    }
}
