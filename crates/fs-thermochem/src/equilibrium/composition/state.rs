//! Actual, frozen ideal-gas populations; no equilibrium projection.
//!
//! Heavy species, including their source-declared internal modes, share Th;
//! the single nucleus-free electron has Te. Each source contributes its own
//! u(T) once, including formation energy. This is a two-temperature caloric
//! state, not a two-temperature chemical-equilibrium or kinetic closure.
//! Matter pressure is R/V (sum_heavy n_i Th + n_e Te). Bulk kinetic energy,
//! electrostatic energy, Coulomb corrections and non-Maxwellian populations
//! are outside this dilute, quasineutral ideal-gas state.

use super::{IdealGasEquilibriumV1, Sum};
use crate::{SpeciesId, ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_qty::{Amount, Energy, Mass, Pressure, Temperature, Volume};
use std::sync::Arc;

const AVOGADRO: f64 = 6.022_140_76e23;

/// Refusal to publish a source-bound, quasineutral frozen state.
#[derive(Clone, Debug, PartialEq)]
pub enum FrozenIonizedGasStateErrorV1 {
    /// Invalid extensive state, temperatures or tolerance.
    Invalid(&'static str),
    /// Net charge relative to actual carriers exceeds the requested tolerance.
    ChargeImbalance { relative_residual: f64 },
    /// The original source refused its species temperature.
    SpeciesEvaluation {
        species: SpeciesId,
        source: ThermochemErrorV1,
    },
    /// A physical result cannot be represented as finite f64.
    Arithmetic(&'static str),
    /// Bounded state storage could not be reserved.
    AllocationRefused,
    /// Cancellation before publication.
    Cancelled,
}
impl core::fmt::Display for FrozenIonizedGasStateErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "frozen ionized gas state: {self:?}")
    }
}
impl std::error::Error for FrozenIonizedGasStateErrorV1 {}

/// Immutable actual populations, original source bank and evaluated caloric state.
///
/// Clones share both the original chemistry bank and the population array.
/// Construction evaluates only populated species; an absent species never
/// imposes a temperature interval on matter which is not present. A later
/// reaction must separately admit the sources of its possible products.
#[derive(Clone, Debug, PartialEq)]
pub struct FrozenIonizedGasStateV1 {
    thermodynamics: IdealGasEquilibriumV1,
    amounts: Arc<[Amount]>,
    volume: Volume,
    heavy_temperature: Temperature,
    electron_temperature: Temperature,
    electron_index: Option<usize>,
    pressure: Pressure,
    heavy_energy: Energy,
    electron_energy: Energy,
    internal_energy: Energy,
    heavy_heat_capacity_j_k: f64,
    electron_heat_capacity_j_k: f64,
    mass: Mass,
    relative_charge_residual: f64,
}
impl FrozenIonizedGasStateV1 {
    /// Retain canonical dense amounts without equilibrating or inserting seeds.
    /// The charge residual is scaled by actual charged amounts, independently
    /// of the neutral inventory. A nonneutral plasma needs a space-charge owner.
    pub fn new(
        thermodynamics: &IdealGasEquilibriumV1,
        amounts: &[Amount],
        volume: Volume,
        heavy_temperature: Temperature,
        electron_temperature: Temperature,
        relative_charge_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FrozenIonizedGasStateErrorV1> {
        use FrozenIonizedGasStateErrorV1 as E;
        if cancelled() {
            return Err(E::Cancelled);
        }
        if amounts.len() != thermodynamics.models().len()
            || !volume.value().is_finite()
            || volume.value() <= 0.0
            || !heavy_temperature.value().is_finite()
            || heavy_temperature.value() <= 0.0
            || !electron_temperature.value().is_finite()
            || electron_temperature.value() <= 0.0
            || !relative_charge_tolerance.is_finite()
            || relative_charge_tolerance <= 0.0
            || relative_charge_tolerance >= 1.0
        {
            return Err(E::Invalid(
                "canonical amounts, positive finite V, Th, Te and charge tolerance in (0,1) required",
            ));
        }
        let elemental = thermodynamics.elemental_matrix();
        let electron_index = (0..amounts.len()).find(|&i| {
            thermodynamics.charge_vector().get(i) == Some(-1)
                && (0..elemental.row_count()).all(|k| elemental.get(k, i) == Some(0))
        });
        if electron_index.is_none() && heavy_temperature != electron_temperature {
            return Err(E::Invalid(
                "distinct electron temperature requires an electron source model",
            ));
        }
        let mut charge_scale = 0.0_f64;
        let mut heavy_amount = Sum::default();
        for (i, amount) in amounts.iter().enumerate() {
            let n = amount.value();
            if !n.is_finite() || n < 0.0 {
                return Err(E::Invalid("amounts must be finite and nonnegative"));
            }
            if Some(i) != electron_index {
                heavy_amount.add(n);
            }
            let z = thermodynamics
                .charge_vector()
                .get(i)
                .expect("bound charge axis") as f64;
            // Normalize before multiplying by charge, preserving subnormal
            // carrier populations and avoiding overflow in n*z.
            if z != 0.0 {
                charge_scale = charge_scale.max(n);
            }
        }
        if !heavy_amount.value().is_finite() || heavy_amount.value() <= 0.0 {
            return Err(E::Invalid(
                "a positive finite nuclear inventory is required",
            ));
        }
        let mut charge = Sum::default();
        let mut absolute_charge = Sum::default();
        if charge_scale > 0.0 {
            for (i, amount) in amounts.iter().enumerate() {
                let z = thermodynamics
                    .charge_vector()
                    .get(i)
                    .expect("bound charge axis") as f64;
                if z != 0.0 {
                    let q = (amount.value() / charge_scale) * z;
                    charge.add(q);
                    absolute_charge.add(q.abs());
                }
            }
        }
        let relative_charge_residual = if charge_scale > 0.0 {
            charge.value().abs() / absolute_charge.value()
        } else {
            0.0
        };
        if !relative_charge_residual.is_finite()
            || relative_charge_residual > relative_charge_tolerance
        {
            return Err(E::ChargeImbalance {
                relative_residual: relative_charge_residual,
            });
        }
        let mut energy = [Sum::default(); 2];
        let mut heat_capacity = [Sum::default(); 2];
        let mut pressure = Sum::default();
        let mut mass = Sum::default();
        for (i, (model, amount)) in thermodynamics.models().iter().zip(amounts).enumerate() {
            if cancelled() {
                return Err(E::Cancelled);
            }
            let n = amount.value();
            if n == 0.0 {
                continue;
            }
            let pool = usize::from(Some(i) == electron_index);
            let t = if pool == 1 {
                electron_temperature
            } else {
                heavy_temperature
            };
            let properties =
                model
                    .evaluate_properties(t)
                    .map_err(|source| E::SpeciesEvaluation {
                        species: model.species().clone(),
                        source,
                    })?;
            let cv = properties.cp().value() - R;
            if !cv.is_finite() || cv <= 0.0 {
                return Err(E::Invalid("populated sources require positive frozen Cv"));
            }
            energy[pool].add(n * properties.u().value());
            heat_capacity[pool].add(n * cv);
            pressure.add((n / volume.value()) * R * t.value());
            mass.add(n * model.molar_mass().value());
            let density = (n / volume.value()) * AVOGADRO;
            if !density.is_finite() || density <= 0.0 {
                return Err(E::Arithmetic("populated number density"));
            }
        }
        let heavy_energy = energy[0].value();
        let electron_energy = energy[1].value();
        let internal_energy = heavy_energy + electron_energy;
        let p = pressure.value();
        let m = mass.value();
        let cv_h = heat_capacity[0].value();
        let cv_e = heat_capacity[1].value();
        if !heavy_energy.is_finite()
            || !electron_energy.is_finite()
            || !internal_energy.is_finite()
            || !p.is_finite()
            || p <= 0.0
            || !m.is_finite()
            || m <= 0.0
            || !cv_h.is_finite()
            || cv_h <= 0.0
            || !cv_e.is_finite()
        {
            return Err(E::Arithmetic(
                "frozen pressure, energy, heat capacity or mass",
            ));
        }
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(amounts.len())
            .map_err(|_| E::AllocationRefused)?;
        retained.extend_from_slice(amounts);
        if cancelled() {
            return Err(E::Cancelled);
        }
        Ok(Self {
            thermodynamics: thermodynamics.clone(),
            amounts: retained.into(),
            volume,
            heavy_temperature,
            electron_temperature,
            electron_index,
            pressure: Pressure::new(p),
            heavy_energy: Energy::new(heavy_energy),
            electron_energy: Energy::new(electron_energy),
            internal_energy: Energy::new(internal_energy),
            heavy_heat_capacity_j_k: cv_h,
            electron_heat_capacity_j_k: cv_e,
            mass: Mass::new(m),
            relative_charge_residual,
        })
    }
    /// Complete original source bank and canonical species axes.
    pub const fn thermodynamics(&self) -> &IdealGasEquilibriumV1 {
        &self.thermodynamics
    }
    /// Actual amounts, including exact absent species, in canonical order.
    pub fn amounts(&self) -> &[Amount] {
        &self.amounts
    }
    /// Parcel volume [m3].
    pub const fn volume(&self) -> Volume {
        self.volume
    }
    /// Heavy-species and heavy internal-mode temperature.
    pub const fn heavy_temperature(&self) -> Temperature {
        self.heavy_temperature
    }
    /// Maxwellian electron temperature.
    pub const fn electron_temperature(&self) -> Temperature {
        self.electron_temperature
    }
    /// Partial-pressure sum at the actual species temperatures.
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    /// Structural nucleus-free charge -1 row, independent of its species name.
    pub const fn electron_index(&self) -> Option<usize> {
        self.electron_index
    }
    /// Actual number density [m^-3]; None means an invalid canonical index.
    pub fn number_density_per_m3(&self, index: usize) -> Option<f64> {
        self.amounts
            .get(index)
            .map(|n| (n.value() / self.volume.value()) * AVOGADRO)
    }
    /// Temperature used to query the original source of this species.
    pub fn species_temperature(&self, index: usize) -> Option<Temperature> {
        self.amounts.get(index).map(|_| {
            if Some(index) == self.electron_index {
                self.electron_temperature
            } else {
                self.heavy_temperature
            }
        })
    }
    /// Total matter internal energy with source formation energies counted once.
    pub const fn internal_energy(&self) -> Energy {
        self.internal_energy
    }
    /// Heavy-species internal energy including their formation reference.
    pub const fn heavy_internal_energy(&self) -> Energy {
        self.heavy_energy
    }
    /// Electron internal energy including its original source reference.
    pub const fn electron_internal_energy(&self) -> Energy {
        self.electron_energy
    }
    /// Frozen heavy heat capacity dUh/dTh [J/K].
    pub const fn heavy_heat_capacity_j_k(&self) -> f64 {
        self.heavy_heat_capacity_j_k
    }
    /// Frozen electron heat capacity dUe/dTe [J/K]; zero if electrons are absent.
    pub const fn electron_heat_capacity_j_k(&self) -> f64 {
        self.electron_heat_capacity_j_k
    }
    /// Sum of original source masses, including electrons.
    pub const fn mass(&self) -> Mass {
        self.mass
    }
    /// Actual |sum z_i n_i| / sum |z_i n_i|, independent of neutral abundance.
    pub const fn relative_charge_residual(&self) -> f64 {
        self.relative_charge_residual
    }
}
