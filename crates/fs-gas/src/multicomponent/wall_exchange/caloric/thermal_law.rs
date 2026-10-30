//! Original h/u storage, reference geometry and latent plateaus.

use super::*;

#[derive(Clone)]
pub(super) struct EnergyStorage<'law> {
    pub law: Arc<dyn ThermalLaw + Send + Sync + 'law>,
    pub initial: DepositedCell,
    pub range_j_kg: [f64; 2],
}

impl core::fmt::Debug for EnergyStorage<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("EnergyStorage")
            .field("initial", &self.initial)
            .field("range_j_kg", &self.range_j_kg)
            .field("basis", &self.law.storage_basis())
            .finish_non_exhaustive()
    }
}
impl PartialEq for EnergyStorage<'_> {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.law, &other.law)
            && self.initial == other.initial
            && self.range_j_kg == other.range_j_kg
    }
}

impl<'law> EnergyStorage<'law> {
    pub fn new<L: ThermalLaw + Send + Sync + 'law>(
        law: L,
        initial: DepositedCell,
        range_j_kg: [f64; 2],
    ) -> Result<Self, WallExchangeError> {
        let [lower, upper] = range_j_kg;
        if !lower.is_finite()
            || !upper.is_finite()
            || lower >= upper
            || !positive(initial.mass_kg)
            || !positive(initial.geometry.volume_m3)
            || !nonnegative(initial.geometry.surface_area_m2)
            || !initial
                .geometry
                .start_m
                .iter()
                .chain(&initial.geometry.end_m)
                .all(|v| v.is_finite())
            || !(lower..=upper).contains(&initial.specific_h_j_kg)
            || !nonnegative(law.maximum_temperature_slope())
        {
            return Err(WallExchangeError::Invalid(
                "finite caloric wall source or geometry",
            ));
        }
        let source = Self {
            law: Arc::new(law),
            initial,
            range_j_kg,
        };
        let temperatures = [
            source.at(lower)?,
            source.at(source.initial.specific_h_j_kg)?,
            source.at(upper)?,
        ];
        if temperatures[0] > temperatures[1] || temperatures[1] > temperatures[2] {
            return Err(WallExchangeError::Invalid(
                "nonmonotone caloric wall source",
            ));
        }
        Ok(source)
    }

    pub fn contact(
        &self,
        energy: f64,
        conductance: f64,
        drag: WallDrag,
    ) -> Result<WallContact, WallExchangeError> {
        self.at(energy)?;
        let contact = WallContact {
            drag,
            thermal: WallThermal::Caloric {
                mass_kg: self.initial.mass_kg,
                specific_energy_j_kg: energy,
                conductance_w_per_k: conductance,
            },
        };
        contact.validate()?;
        Ok(contact)
    }

    pub fn state(&self, contact: WallContact) -> Result<(f64, f64), WallExchangeError> {
        contact.validate()?;
        let WallThermal::Caloric {
            mass_kg,
            specific_energy_j_kg,
            conductance_w_per_k,
        } = contact.thermal
        else {
            return Err(WallExchangeError::Invalid(
                "thermal-law wall requires its energy coordinate",
            ));
        };
        if mass_kg.to_bits() != self.initial.mass_kg.to_bits() {
            return Err(WallExchangeError::Invalid("caloric wall mass changed"));
        }
        self.at(specific_energy_j_kg)?;
        Ok((specific_energy_j_kg, conductance_w_per_k))
    }

    pub fn at(&self, energy: f64) -> Result<f64, WallExchangeError> {
        if !(self.range_j_kg[0]..=self.range_j_kg[1]).contains(&energy) {
            return Err(WallExchangeError::Invalid(
                "caloric wall energy outside source interval",
            ));
        }
        let mut candidate = self.initial.clone();
        candidate.specific_h_j_kg = energy;
        self.law.validate_cell(&candidate)?;
        let temperature = self.law.temperature_k(energy)?;
        if !positive(temperature) {
            return Err(WallExchangeError::Invalid("caloric wall temperature"));
        }
        Ok(temperature)
    }

    pub fn after(
        &self,
        contact: WallContact,
        heat: f64,
    ) -> Result<StorageUpdate, WallExchangeError> {
        let (initial, conductance) = self.state(contact)?;
        let energy = initial + heat / self.initial.mass_kg;
        Ok(StorageUpdate {
            contact: self.contact(energy, conductance, contact.drag)?,
            temperature_k: self.at(energy)?,
            energy_change_j: self.initial.mass_kg * (energy - initial),
            // T(e) can have an exact plateau. Its global slope bound is not
            // a local tangent; use safeguarded bisection without a fitted Cp.
            temperature_slope_k_per_j: None,
        })
    }
}
