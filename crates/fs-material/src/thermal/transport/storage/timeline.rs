//! Bind total sensible storage to the existing conservative thermal owner.

use super::{HeatStorageAmount, PreparedHeatStorage, ThermalTransportError, invalid, positive};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{DepositedCell, DepositionError, ThermalLaw};

/// One assembly used as one uniform-temperature thermal cell.
///
/// The explicit total mass only converts total joules to the timeline's J/kg
/// coordinate. Volumetric source capacities remain volumetric; no measured
/// density or source specific heat is manufactured. Conduction, boundaries,
/// spatial resolution and time stepping remain with their existing owners.
/// Call `admit_cell` for every cell before inserting its birth event; the generic
/// `ThermalLaw` interface does not carry a cell's physical mass into evaluation.
#[derive(Clone, Debug)]
pub struct HeatStorageThermalLaw<'a> {
    storage: &'a PreparedHeatStorage,
    reference_mass: f64,
    datum: f64,
    span: [f64; 2],
    maximum_slope: f64,
    sensible_entropy: bool,
    identity: ContentHash,
}

impl<'a> HeatStorageThermalLaw<'a> {
    pub fn try_new(
        storage: &'a PreparedHeatStorage,
        reference_mass_kg: f64,
        datum_temperature_k: f64,
    ) -> Result<Self, ThermalTransportError> {
        if !positive(reference_mass_kg) {
            return Err(invalid(
                "assembly reference mass must be finite and positive",
            ));
        }
        let mut known_mass = 0.0;
        let mut correction = 0.0;
        let mut has_volume = false;
        for component in storage.components() {
            match component.amount() {
                HeatStorageAmount::MassKg(mass) => {
                    let corrected = mass - correction;
                    let next = known_mass + corrected;
                    correction = (next - known_mass) - corrected;
                    known_mass = next;
                }
                HeatStorageAmount::ReferenceVolumeM3(_) => has_volume = true,
            }
        }
        if !known_mass.is_finite()
            || (has_volume && known_mass >= reference_mass_kg)
            || (!has_volume && !same_mass(known_mass, reference_mass_kg))
        {
            return Err(invalid(
                "reference mass contradicts the known constituent masses",
            ));
        }
        storage.heat_capacity_j_k(datum_temperature_k)?;
        let temperatures = storage.temperature_span_k();
        let span = [
            storage.energy_change_j(datum_temperature_k, temperatures[0])? / reference_mass_kg,
            storage.energy_change_j(datum_temperature_k, temperatures[1])? / reference_mass_kg,
        ];
        if !span.iter().all(|x| x.is_finite())
            || span[0] >= span[1]
            || (temperatures[0] != datum_temperature_k && span[0] == 0.0)
            || (temperatures[1] != datum_temperature_k && span[1] == 0.0)
        {
            return Err(invalid("assembly specific storage is not representable"));
        }
        let minimum = storage
            .knots
            .iter()
            .map(|k| k[1])
            .fold(f64::INFINITY, f64::min);
        let maximum_slope = (reference_mass_kg / minimum).next_up();
        if !positive(maximum_slope) {
            return Err(invalid("assembly temperature slope is not representable"));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-material.heat-storage-thermal-law.v1");
        hash.update(storage.identity().as_bytes());
        hash.update(&reference_mass_kg.to_bits().to_le_bytes());
        hash.update(&datum_temperature_k.to_bits().to_le_bytes());
        Ok(Self {
            storage,
            reference_mass: reference_mass_kg,
            datum: datum_temperature_k,
            span,
            maximum_slope,
            sensible_entropy: false,
            identity: hash.finalize(),
        })
    }

    #[must_use]
    pub const fn storage(&self) -> &'a PreparedHeatStorage {
        self.storage
    }
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn reference_mass_kg(&self) -> f64 {
        self.reference_mass
    }
    #[must_use]
    pub const fn datum_temperature_k(&self) -> f64 {
        self.datum
    }
    #[must_use]
    pub const fn specific_storage_span_j_kg(&self) -> [f64; 2] {
        self.span
    }

    /// Select the same reversible sensible-entropy convention as the source path.
    #[must_use]
    pub fn with_sensible_entropy(mut self) -> Self {
        if !self.sensible_entropy {
            let mut hash =
                DomainHasher::new("org.frankensim.fs-material.heat-storage-sensible-entropy.v1");
            hash.update(self.identity.as_bytes());
            self.identity = hash.finalize();
            self.sensible_entropy = true;
        }
        self
    }
    #[must_use]
    pub const fn sensible_entropy_enabled(&self) -> bool {
        self.sensible_entropy
    }

    pub fn specific_storage_j_kg(&self, temperature_k: f64) -> Result<f64, ThermalTransportError> {
        let h = self.storage.energy_change_j(self.datum, temperature_k)? / self.reference_mass;
        if !h.is_finite() || (temperature_k != self.datum && h == 0.0) {
            return Err(invalid("assembly specific storage is not representable"));
        }
        Ok(h)
    }

    /// The entire assembly is one cell; changing its physical mass requires a
    /// newly prepared amount binding. A 64-epsilon check permits roundoff only.
    pub fn admit_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        if !positive(cell.mass_kg) || !same_mass(cell.mass_kg, self.reference_mass) {
            return Err(DepositionError::InvalidInput(
                "cell mass differs from its heat-storage assembly",
            ));
        }
        self.temperature_k(cell.specific_h_j_kg)?;
        Ok(())
    }
}

impl ThermalLaw for HeatStorageThermalLaw<'_> {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let [low, high] = self.span;
        if !h.is_finite() || h < low || h > high {
            return Err(DepositionError::MaterialDomain);
        }
        let temperatures = self.storage.temperature_span_k();
        if h == low {
            return Ok(temperatures[0]);
        }
        if h == high {
            return Ok(temperatures[1]);
        }
        let energy = h * self.reference_mass;
        if !energy.is_finite() || (h != 0.0 && energy == 0.0) {
            return Err(DepositionError::MaterialDomain);
        }
        self.storage
            .temperature_after_energy_change(self.datum, energy)
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        self.maximum_slope
    }

    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h: f64,
        final_h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        if !self.sensible_entropy {
            return Ok(None);
        }
        let initial = self.temperature_k(initial_h)?;
        let final_t = self.temperature_k(final_h)?;
        if initial_h == final_h {
            return Ok(Some(0.0));
        }
        let change = final_h - initial_h;
        let entropy = if initial == final_t {
            change / initial
        } else {
            let entropy = self
                .storage
                .sensible_entropy_change_j_k(initial, final_t)
                .map_err(|_| DepositionError::MaterialDomain)?;
            let energy = self
                .storage
                .energy_change_j(initial, final_t)
                .map_err(|_| DepositionError::MaterialDomain)?;
            change * (entropy / energy)
        };
        if !entropy.is_finite()
            || entropy == 0.0
            || entropy.is_sign_positive() != (final_h > initial_h)
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(Some(entropy))
    }
}

fn same_mass(a: f64, b: f64) -> bool {
    (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}
