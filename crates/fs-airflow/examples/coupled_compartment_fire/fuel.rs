//! Mechanism selection and manufactured solid composition for the shared example.
//!
//! The charring solid is defined by mass as 0.8 CH4 + 0.2 condensed carbon.
//! This closes its elemental inventory; it is not a fitted wood pyrolysate.

use fs_gas::multicomponent::MixtureGas;
use fs_thermochem::kinetics::IdealGasKineticsV1;
use fs_thermochem::kinetics::gri30::{Gri30Calorics, cantera_gri30};
use fs_thermochem::kinetics::hydrogen::{HydrogenCalorics, cantera_h2o2};
use std::error::Error;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fuel {
    Hydrogen,
    Methane,
    CharringMethane,
}

impl Fuel {
    pub fn parse(name: &str) -> Result<Self, Box<dyn Error>> {
        match name {
            "hydrogen" => Ok(Self::Hydrogen),
            "methane" => Ok(Self::Methane),
            "charring-methane" => Ok(Self::CharringMethane),
            _ => Err("fuel must be hydrogen, methane or charring-methane".into()),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Hydrogen => "hydrogen",
            Self::Methane => "methane",
            Self::CharringMethane => "charring-methane",
        }
    }

    pub fn species(self) -> &'static str {
        match self {
            Self::Hydrogen => "H2",
            Self::Methane | Self::CharringMethane => "CH4",
        }
    }

    pub fn mechanism(self) -> Result<IdealGasKineticsV1, Box<dyn Error>> {
        Ok(match self {
            Self::Hydrogen => cantera_h2o2(HydrogenCalorics::ContinuousAt300K)?,
            Self::Methane | Self::CharringMethane => {
                cantera_gri30(Gri30Calorics::ContinuousAt300K)?
            }
        })
    }

    pub fn source(self) -> &'static str {
        match self {
            Self::Hydrogen => "Cantera-3.2.0-h2o2-10-species-29-reactions",
            Self::Methane | Self::CharringMethane => "Cantera-3.2.0-gri30-53-species-325-reactions",
        }
    }

    pub fn initial_temperature_k(self) -> f64 {
        match self {
            Self::Hydrogen => 1000.0,
            Self::Methane | Self::CharringMethane => 1500.0,
        }
    }

    pub fn maximum_temperature_k(self) -> f64 {
        match self {
            Self::Hydrogen => 3500.0,
            Self::Methane | Self::CharringMethane => 3000.0,
        }
    }

    pub fn char_yield(self) -> f64 {
        if self == Self::CharringMethane {
            0.2
        } else {
            0.0
        }
    }
}

/// Resolve names against the original canonical mechanism once at construction.
pub struct SpeciesAxis {
    pub fuel: usize,
    pub h2: usize,
    pub water: usize,
    pub nitrogen: usize,
    pub oxygen: usize,
    pub hydrogen_atom: usize,
    pub hydroxyl: usize,
    pub co: Option<usize>,
    pub co2: Option<usize>,
    pub no: Option<usize>,
    pub molar_masses: Vec<f64>,
    pub virgin_atoms_per_kg: Vec<f64>,
    pub char_atoms_per_kg: Vec<f64>,
}

impl SpeciesAxis {
    pub fn new(gas: &MixtureGas, fuel: Fuel) -> Result<Self, Box<dyn Error>> {
        let source = gas.source();
        let elements = source.elemental_matrix();
        let find = |name| elements.species().iter().position(|s| s.as_str() == name);
        let required =
            |name| find(name).ok_or_else(|| format!("mechanism lacks required species {name}"));
        let fuel_index = required(fuel.species())?;
        let molar_masses: Vec<_> = source
            .models()
            .iter()
            .map(|m| m.molar_mass().value())
            .collect();
        let mut char_atoms_per_kg = vec![0.0; elements.row_count()];
        if fuel.char_yield() > 0.0 {
            // The source atomic C mass fixes the molar bookkeeping for the
            // manufactured condensed-carbon residue. Gas-phase atomic C is
            // not emitted and its gas enthalpy is not used for the solid.
            let carbon = required("C")?;
            for (e, amount) in char_atoms_per_kg.iter_mut().enumerate() {
                *amount = elements.get(e, carbon).expect("element/species axis") as f64
                    / molar_masses[carbon];
            }
        }
        let virgin_atoms_per_kg = (0..elements.row_count())
            .map(|e| {
                (1.0 - fuel.char_yield())
                    * elements.get(e, fuel_index).expect("element/species axis") as f64
                    / molar_masses[fuel_index]
                    + fuel.char_yield() * char_atoms_per_kg[e]
            })
            .collect();
        Ok(Self {
            fuel: fuel_index,
            h2: required("H2")?,
            water: required("H2O")?,
            nitrogen: required("N2")?,
            oxygen: required("O2")?,
            hydrogen_atom: required("H")?,
            hydroxyl: required("OH")?,
            co: find("CO"),
            co2: find("CO2"),
            no: find("NO"),
            molar_masses,
            virgin_atoms_per_kg,
            char_atoms_per_kg,
        })
    }

    pub fn mass(masses: &[f64], species: Option<usize>) -> f64 {
        species.map_or(0.0, |i| masses[i])
    }
}
