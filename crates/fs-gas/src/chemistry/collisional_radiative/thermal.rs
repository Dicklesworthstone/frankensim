//! Original-source thermal exchange and held work for the coupled owner.
//!
//! With stage capacities Ce, Ch and conductance G, the linearized temperature
//! difference obeys D'=F-G(1/Ce+1/Ch)D. Its exponential integral supplies one
//! equal-and-opposite transfer, including the changing electron population in
//! F. Constant capacities, populations, conductance and powers are integrated
//! exactly; temperature-dependent coefficients and calorics retain the outer
//! midpoint approximation. Heavy temperature is always recovered from original
//! source energy differences, never reset to a mixture or electron temperature.

use super::super::energy::{
    ElectronHeavyCollisionSources, ElectronHeavyEnergyError, ExchangeConductance, PreparedExchange,
    source_pool_calorics,
};
use super::*;
use fs_math::det::expm1;
use fs_qty::{Amount, Power};

const CV: f64 = 1.5 * R;

/// Retained collision laws, independently held pool powers and heavy support.
#[derive(Clone, Debug, PartialEq)]
pub struct CollisionalRadiativeThermal {
    /// Original neutral momentum sources and optional initial-bound Coulomb law.
    pub collisions: Option<ElectronHeavyCollisionSources>,
    /// Signed external power supplied to the electron pool [W].
    pub electron_power: Power,
    /// Signed external power supplied to the heavy pool [W].
    pub heavy_power: Power,
    /// Minimum supported temperature for all populated/created heavy species [K].
    pub minimum_heavy_temperature_k: f64,
    /// Maximum supported temperature for all populated/created heavy species [K].
    pub maximum_heavy_temperature_k: f64,
    /// Bound on each original-source heavy caloric inversion.
    pub maximum_iterations: usize,
}
impl CollisionalRadiativeThermal {
    /// Select explicit heavy support, initially with no collisions or pool work.
    pub fn new(minimum: Temperature, maximum: Temperature) -> Self {
        Self {
            collisions: None,
            electron_power: Power::new(0.0),
            heavy_power: Power::new(0.0),
            minimum_heavy_temperature_k: minimum.value(),
            maximum_heavy_temperature_k: maximum.value(),
            maximum_iterations: 96,
        }
    }
    pub(super) fn validate(&self) -> Result<(), CollisionalRadiativeError> {
        if !self.electron_power.value().is_finite()
            || !self.heavy_power.value().is_finite()
            || !positive(self.minimum_heavy_temperature_k)
            || !positive(self.maximum_heavy_temperature_k)
            || self.minimum_heavy_temperature_k >= self.maximum_heavy_temperature_k
            || !(8..=256).contains(&self.maximum_iterations)
        {
            return Err(CollisionalRadiativeError::Invalid(
                "thermal powers, support or iteration bound",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
pub(super) struct ThermalStage {
    pub conductance: ExchangeConductance,
    pub capacities: [f64; 2],
    pub temperatures: [f64; 2],
    pub heavy_energy: f64,
    powers: [f64; 2],
}
impl ThermalStage {
    pub fn inactive(self) -> bool {
        self.powers == [0.0; 2]
            && (self.conductance.total() == 0.0 || self.temperatures[0] == self.temperatures[1])
    }
}

pub(super) struct ThermalAdvance {
    pub electron_energy: f64,
    pub heavy_temperature: f64,
    pub heavy_increment: f64,
    pub exchange: f64,
    pub neutral_exchange: f64,
    pub coulomb_exchange: f64,
}

pub(super) struct PreparedThermal<'a> {
    initial: &'a FrozenIonizedGasStateV1,
    inputs: &'a CollisionalRadiativeThermal,
    controls: CollisionalRadiativeSolve,
    exchange: Option<PreparedExchange<'a>>,
    heavy_limits: [f64; 2],
}
impl<'a> PreparedThermal<'a> {
    pub fn new(
        initial: &'a FrozenIonizedGasStateV1,
        inputs: &'a CollisionalRadiativeThermal,
        created_neutral: usize,
        controls: CollisionalRadiativeSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, CollisionalRadiativeError> {
        inputs.validate()?;
        let electron = initial
            .electron_index()
            .ok_or(CollisionalRadiativeError::Invalid("structural electron"))?;
        if initial.amounts()[electron].value() == 0.0 && inputs.electron_power.value() != 0.0 {
            return Err(CollisionalRadiativeError::Invalid(
                "held electron power without initial electrons",
            ));
        }
        let exchange = inputs
            .collisions
            .as_ref()
            .map(|sources| {
                let prepared =
                    PreparedExchange::new(initial, sources, created_neutral, &mut *cancelled)
                        .map_err(thermal_error)?;
                // Photoionization can create electrons from exact absence. Require
                // all neutral sources before that happens, including a neutral
                // that later recombination can create from an initially ionized gas.
                for (index, model) in initial.thermodynamics().models().iter().enumerate() {
                    if initial.thermodynamics().charge_vector().get(index) == Some(0)
                        && (initial.amounts()[index].value() > 0.0 || index == created_neutral)
                        && !sources
                            .momentum_sources
                            .iter()
                            .any(|s| s.target_species() == model.species())
                    {
                        return Err(thermal_error(ElectronHeavyEnergyError::MissingSource(
                            model.species().clone(),
                        )));
                    }
                }
                Ok(prepared)
            })
            .transpose()?;
        let mut result = Self {
            initial,
            inputs,
            controls,
            exchange,
            heavy_limits: [0.0; 2],
        };
        result.heavy_limits = [
            result
                .calorics(inputs.minimum_heavy_temperature_k, cancelled)?
                .0,
            result
                .calorics(inputs.maximum_heavy_temperature_k, cancelled)?
                .0,
        ];
        Ok(result)
    }

    pub fn powers(&self) -> [f64; 2] {
        [
            self.inputs.heavy_power.value(),
            self.inputs.electron_power.value(),
        ]
    }

    pub fn calorics(
        &self,
        th: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), CollisionalRadiativeError> {
        source_pool_calorics(self.initial, false, th, cancelled).map_err(thermal_error)
    }

    pub fn stage(
        &self,
        amounts: &[Amount],
        th: f64,
        te: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermalStage, CollisionalRadiativeError> {
        let (energy, capacity) = self.calorics(th, cancelled)?;
        let conductance = if let Some(exchange) = &self.exchange {
            self.admit_charge(amounts, cancelled)?;
            exchange
                .coupling_at_point(amounts, self.initial.volume(), th, te, cancelled)
                .map_err(thermal_error)?
        } else {
            ExchangeConductance::default()
        };
        let ne = amounts[self.initial.electron_index().expect("admitted electron")].value();
        Ok(ThermalStage {
            conductance,
            capacities: [capacity, CV * ne],
            temperatures: [th, te],
            heavy_energy: energy,
            powers: self.powers(),
        })
    }

    fn admit_charge(
        &self,
        amounts: &[Amount],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), CollisionalRadiativeError> {
        let charges = self.initial.thermodynamics().charge_vector();
        let scale = amounts
            .iter()
            .enumerate()
            .filter(|(i, _)| charges.get(*i) != Some(0))
            .map(|(_, n)| n.value())
            .fold(0.0, f64::max);
        if scale == 0.0 {
            return Ok(());
        }
        let mut sums = [0.0_f64; 2];
        let mut corrections = [0.0_f64; 2];
        for (i, amount) in amounts.iter().enumerate() {
            poll(cancelled)?;
            let charge = charges.get(i).expect("canonical charge axis") as f64;
            let q = (amount.value() / scale) * charge;
            for (slot, value) in [q, q.abs()].into_iter().enumerate() {
                let adjusted = value - corrections[slot];
                let next = sums[slot] + adjusted;
                corrections[slot] = (next - sums[slot]) - adjusted;
                sums[slot] = next;
            }
        }
        let relative = sums[0].abs() / sums[1];
        if !relative.is_finite() || relative > self.controls.relative_charge_tolerance {
            return Err(CollisionalRadiativeError::State(
                FrozenIonizedGasStateErrorV1::ChargeImbalance {
                    relative_residual: relative,
                },
            ));
        }
        Ok(())
    }

    /// Mean chemical/radiative forcing is second-order accurate inside the
    /// exponential thermal relaxation. The same signed transfer appears once
    /// in each pool; original calorics recover the heavy endpoint afterward.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        stage: ThermalStage,
        temperatures: [f64; 2],
        electron_amount: f64,
        heavy_energy: f64,
        chemical_electron_energy: f64,
        extent: f64,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermalAdvance, CollisionalRadiativeError> {
        let [ph, pe] = stage.powers.map(|p| p * dt);
        let mut transfer = 0.0;
        let conductance = stage.conductance.total();
        if conductance > 0.0 {
            let [ch, ce] = stage.capacities;
            if !positive(ce) || !positive(ch) {
                return Err(CollisionalRadiativeError::Arithmetic(
                    "positive stage thermal capacities",
                ));
            }
            let capacity = if ce <= ch {
                ce / (1.0 + ce / ch)
            } else {
                ch / (1.0 + ch / ce)
            };
            let dose = (conductance / capacity) * dt;
            let (depleted, forcing) = relaxation(dose);
            // Linearize temperature in (electron energy, electron population)
            // together. Omitting the population term would double-count the
            // energy needed to bring each newly created electron to stage Te.
            let difference = stage.temperatures[1] - stage.temperatures[0]
                + ((CV * electron_amount) / ce) * (temperatures[1] - stage.temperatures[1])
                - (heavy_energy - stage.heavy_energy) / ch;
            let impulse = (chemical_electron_energy + pe - CV * stage.temperatures[1] * extent)
                / ce
                - ph / ch;
            transfer = capacity * (difference * depleted + impulse * forcing);
        }
        let heavy_increment = ph + transfer;
        let target = heavy_energy + heavy_increment;
        let electron_energy = chemical_electron_energy + pe - transfer;
        if !transfer.is_finite() || !target.is_finite() || !electron_energy.is_finite() {
            return Err(CollisionalRadiativeError::Arithmetic(
                "held work and thermal exchange",
            ));
        }
        let th = if heavy_increment == 0.0 {
            temperatures[0]
        } else {
            self.temperature(target, temperatures[0], cancelled)?
        };
        Ok(ThermalAdvance {
            electron_energy,
            heavy_temperature: th,
            heavy_increment,
            exchange: transfer,
            neutral_exchange: if conductance == 0.0 {
                0.0
            } else {
                transfer * (stage.conductance.neutral / conductance)
            },
            coulomb_exchange: if conductance == 0.0 {
                0.0
            } else {
                transfer * (stage.conductance.coulomb / conductance)
            },
        })
    }

    fn temperature(
        &self,
        target: f64,
        guess: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, CollisionalRadiativeError> {
        if !(self.heavy_limits[0]..=self.heavy_limits[1]).contains(&target) {
            return Err(CollisionalRadiativeError::OutsideDomain);
        }
        let mut lo = self.inputs.minimum_heavy_temperature_k;
        let mut hi = self.inputs.maximum_heavy_temperature_k;
        let mut th = guess;
        for _ in 0..self.inputs.maximum_iterations {
            poll(cancelled)?;
            let (energy, capacity) = self.calorics(th, cancelled)?;
            let residual = energy - target;
            let tolerance = (self.controls.energy_tolerance_j / 64.0).min(
                capacity
                    * (self.controls.temperature_tolerance_k
                        + self.controls.relative_tolerance * th)
                    / 256.0,
            );
            if residual.abs() <= tolerance {
                return Ok(th);
            }
            if residual > 0.0 {
                hi = th;
            } else {
                lo = th;
            }
            let next = th - residual / capacity;
            th = if next > lo && next < hi {
                next
            } else {
                lo + 0.5 * (hi - lo)
            };
        }
        Err(thermal_error(ElectronHeavyEnergyError::Nonconvergence))
    }
}

// 1-exp(-z) and 1-(1-exp(-z))/z; the second factor integrates constant
// forcing through the exchange relaxation. Their zero and stiff limits are
// continuous, without subtracting two almost-equal energies or temperatures.
fn relaxation(dose: f64) -> (f64, f64) {
    if !dose.is_finite() {
        return (1.0, 1.0);
    }
    let depleted = -expm1(-dose);
    let forcing = if dose < 1e-3 {
        dose * (0.5
            + dose * (-1.0 / 6.0 + dose * (1.0 / 24.0 + dose * (-1.0 / 120.0 + dose / 720.0))))
    } else {
        1.0 - depleted / dose
    };
    (depleted, forcing)
}

pub(super) fn thermal_error(error: ElectronHeavyEnergyError) -> CollisionalRadiativeError {
    match error {
        ElectronHeavyEnergyError::Cancelled => CollisionalRadiativeError::Cancelled,
        ElectronHeavyEnergyError::OutsideDomain => CollisionalRadiativeError::OutsideDomain,
        error => CollisionalRadiativeError::Thermal(error),
    }
}
