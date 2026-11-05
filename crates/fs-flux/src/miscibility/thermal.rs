//! Adiabatic finite-rate miscibility with one uniform evolving temperature.
//!
//! The original pure-liquid calorics, homogeneous mixing enthalpy and the same
//! square-gradient energy close the energy account. This is the fast internal
//! thermal-equilibration limit: it does not resolve local heat conduction,
//! thermodiffusion, thermal expansion, mechanical work or an external reservoir.
//! Both component inventories, pressure and fixed partial molar volumes remain
//! those of the original field. Constant kappa contributes internal energy but
//! no entropy; fixed total volume removes the difference of pV between states.
//!
//! One first-order split advances molecular diffusion at the preceding T, then
//! solves the global sensible-heat balance at the actual new compositions. No
//! cell is flashed to bulk equilibrium. The accepted whole-field energy and
//! entropy are audited separately: an isothermal free-energy decrease alone
//! does not guarantee admissible adiabatic evolution. A negative entropy change
//! refuses the candidate; the caller must refine the physical interval.

use super::temperature::BinaryMiscibilityTemperatureChange;
use super::{
    BinaryMiscibilityError, BinaryMiscibilityNetwork, BinaryMiscibilitySolve,
    BinaryMiscibilityState, BinaryMiscibilityStep, Sum, finite, nonnegative, poll, positive,
};
use fs_qty::Temperature;
use fs_thermochem::solution::caloric::{BinaryCaloricErrorV1, BinaryLiquidCaloricModelV1};
use std::sync::Arc;

/// An incompatible caloric owner or incomplete/thermodynamically inadmissible
/// interval. Every refusal leaves the original field and caloric sources intact.
#[derive(Clone, Debug, PartialEq)]
pub enum BinaryMiscibilityThermalError {
    Invalid(&'static str),
    Solver(BinaryMiscibilityError),
    Caloric(BinaryCaloricErrorV1),
    HeatOutsideBracket {
        residuals_j: [f64; 2],
    },
    IterationLimit,
    Energy {
        residual_j: f64,
        tolerance_j: f64,
    },
    Entropy {
        change_j_per_k: f64,
        tolerance_j_per_k: f64,
    },
}
impl core::fmt::Display for BinaryMiscibilityThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adiabatic spatial miscibility: {self:?}")
    }
}
impl std::error::Error for BinaryMiscibilityThermalError {}
impl From<BinaryMiscibilityError> for BinaryMiscibilityThermalError {
    fn from(error: BinaryMiscibilityError) -> Self {
        Self::Solver(error)
    }
}
impl From<BinaryCaloricErrorV1> for BinaryMiscibilityThermalError {
    fn from(error: BinaryCaloricErrorV1) -> Self {
        Self::Caloric(error)
    }
}
type Result<T> = core::result::Result<T, BinaryMiscibilityThermalError>;

/// Explicit temperature-root and whole-interval thermodynamic allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryMiscibilityThermalSolve {
    /// Strictly increasing bracket inside the original pair and both Cp domains.
    pub temperature_bracket: [Temperature; 2],
    pub temperature_tolerance_k: f64,
    pub energy_tolerance_j: f64,
    /// Numerical allowance for the final isolated system entropy change.
    /// Applied once, never multiplied by cells, faces or root iterations.
    pub entropy_tolerance_j_per_k: f64,
    /// At most 256 scalar bisections, with original Cp integrals at each trial.
    pub maximum_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
struct ThermalOwner {
    caloric_model: BinaryLiquidCaloricModelV1,
    // Original conserved species totals, retaining high and low reduction words.
    // Pure-species enthalpy datums cancel; no per-cell thermal solve is needed.
    component_moles: [[f64; 2]; 2],
}

/// Existing spatial state bound to its immutable original pure-liquid calorics.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityThermalField {
    network: BinaryMiscibilityNetwork,
    state: BinaryMiscibilityState,
    owner: Arc<ThermalOwner>,
}

/// Complete finite-rate diffusion and the independently checked thermal account.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityThermalStep {
    /// Original isothermal diffusion result at the preceding temperature.
    pub diffusion: BinaryMiscibilityStep,
    /// Parameter-change account at the completed, unflashed spatial composition.
    pub temperature_change: BinaryMiscibilityTemperatureChange,
    pub mixing_enthalpy_change_j: f64,
    pub interfacial_energy_change_j: f64,
    pub sensible_enthalpy_change_j: f64,
    /// Sum of the three enthalpy/internal-energy increments, with zero heat input.
    pub energy_residual_j: f64,
    pub mixing_entropy_change_j_per_k: f64,
    pub sensible_entropy_change_j_per_k: f64,
    pub entropy_change_j_per_k: f64,
    pub thermal_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityThermalCandidate {
    field: BinaryMiscibilityThermalField,
    step: BinaryMiscibilityThermalStep,
}
impl BinaryMiscibilityThermalCandidate {
    pub const fn field(&self) -> &BinaryMiscibilityThermalField {
        &self.field
    }
    pub const fn step(&self) -> &BinaryMiscibilityThermalStep {
        &self.step
    }
    pub fn into_parts(self) -> (BinaryMiscibilityThermalField, BinaryMiscibilityThermalStep) {
        (self.field, self.step)
    }
}

impl BinaryMiscibilityNetwork {
    /// Retain the exact original pair and ordered pure-liquid Cp owner. The
    /// source-qualified material adapter can bind the corresponding source law;
    /// raw numerical callers must qualify these liquid calorics and pressure.
    pub fn bind_caloric_field(
        &self,
        state: &BinaryMiscibilityState,
        caloric: &BinaryLiquidCaloricModelV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityThermalField> {
        poll(&mut cancelled)?;
        self.admit(state)?;
        if self.model().thermodynamics() != caloric.thermodynamics() {
            return Err(BinaryMiscibilityThermalError::Invalid(
                "thermal field requires the same full original regular-solution pair",
            ));
        }
        let temperature = self.model().temperature();
        for source in caloric.heat_capacities() {
            poll(&mut cancelled)?;
            source.enthalpy_change_j_per_mol(temperature, temperature)?;
        }
        let mut inventory = [Sum::default(); 2];
        for (amounts, remainders) in state.amounts.iter().zip(&state.amount_remainders) {
            poll(&mut cancelled)?;
            for species in 0..2 {
                inventory[species].add(amounts[species])?;
                inventory[species].add(remainders[species])?;
            }
        }
        for total in inventory {
            if !nonnegative(total.value()?) {
                return Err(BinaryMiscibilityThermalError::Invalid(
                    "nonnegative conserved caloric inventories required",
                ));
            }
        }
        let field = BinaryMiscibilityThermalField {
            network: self.clone(),
            state: state.clone(),
            owner: Arc::new(ThermalOwner {
                caloric_model: caloric.clone(),
                component_moles: inventory.map(|sum| [sum.main, sum.correction]),
            }),
        };
        poll(&mut cancelled)?;
        Ok(field)
    }
}

impl BinaryMiscibilityThermalField {
    pub const fn network(&self) -> &BinaryMiscibilityNetwork {
        &self.network
    }
    pub const fn state(&self) -> &BinaryMiscibilityState {
        &self.state
    }
    pub fn caloric_model(&self) -> &BinaryLiquidCaloricModelV1 {
        &self.owner.caloric_model
    }

    /// Advance one closed diffusion interval and convert its actual mixing and
    /// interfacial enthalpy change into sensible heat. The same original caloric
    /// owner and conserved totals survive every accepted interval.
    ///
    /// This split is first order, even if its energy root is very accurate.
    /// A large step can conserve energy and still decrease entropy: such a step
    /// is refused. Refine time to control thermal splitting error as well as
    /// composition error. No flow work or external heat is silently included.
    pub fn adiabatic_candidate(
        &self,
        duration_s: f64,
        diffusion_controls: BinaryMiscibilitySolve,
        thermal_controls: BinaryMiscibilityThermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityThermalCandidate> {
        poll(&mut cancelled)?;
        self.admit_controls(thermal_controls)?;
        let (diffused, diffusion) = self
            .network
            .step_candidate(&self.state, duration_s, diffusion_controls, &mut cancelled)?
            .into_parts();
        let (mixing_enthalpy_change_j, mixing_entropy_change_j_per_k, interfacial_energy_change_j) =
            self.composition_change(&diffused, &mut cancelled)?;
        let structural_energy = total(&[mixing_enthalpy_change_j, interfacial_energy_change_j])?;
        let from = self.network.model().temperature();
        let (temperature, thermal_iterations) =
            self.solve_temperature(from, structural_energy, thermal_controls, &mut cancelled)?;
        let (network, state, temperature_change) = self
            .network
            .temperature_candidate(&diffused, temperature, &mut cancelled)?
            .into_parts();
        // Audit the actual accepted endpoint using the original integrals.
        let sensible_enthalpy_change_j =
            self.sensible_change(from, temperature, false, &mut cancelled)?;
        let sensible_entropy_change_j_per_k =
            self.sensible_change(from, temperature, true, &mut cancelled)?;
        let energy_residual_j = total(&[
            mixing_enthalpy_change_j,
            interfacial_energy_change_j,
            sensible_enthalpy_change_j,
        ])?;
        if energy_residual_j.abs() > thermal_controls.energy_tolerance_j {
            return Err(BinaryMiscibilityThermalError::Energy {
                residual_j: energy_residual_j,
                tolerance_j: thermal_controls.energy_tolerance_j,
            });
        }
        let entropy_change_j_per_k = total(&[
            mixing_entropy_change_j_per_k,
            sensible_entropy_change_j_per_k,
        ])?;
        if entropy_change_j_per_k < -thermal_controls.entropy_tolerance_j_per_k {
            return Err(BinaryMiscibilityThermalError::Entropy {
                change_j_per_k: entropy_change_j_per_k,
                tolerance_j_per_k: thermal_controls.entropy_tolerance_j_per_k,
            });
        }
        let field = Self {
            network,
            state,
            owner: Arc::clone(&self.owner),
        };
        let step = BinaryMiscibilityThermalStep {
            diffusion,
            temperature_change,
            mixing_enthalpy_change_j,
            interfacial_energy_change_j,
            sensible_enthalpy_change_j,
            energy_residual_j,
            mixing_entropy_change_j_per_k,
            sensible_entropy_change_j_per_k,
            entropy_change_j_per_k,
            thermal_iterations,
        };
        poll(&mut cancelled)?;
        Ok(BinaryMiscibilityThermalCandidate { field, step })
    }

    fn admit_controls(&self, controls: BinaryMiscibilityThermalSolve) -> Result<()> {
        let [lo, hi] = controls.temperature_bracket.map(|t| t.value());
        let domain = self.caloric_model().temperature_domain().map(|t| t.value());
        if !positive(lo)
            || !positive(hi)
            || lo >= hi
            || lo < domain[0]
            || hi > domain[1]
            || !positive(controls.temperature_tolerance_k)
            || !nonnegative(controls.energy_tolerance_j)
            || !nonnegative(controls.entropy_tolerance_j_per_k)
            || !(1..=256).contains(&controls.maximum_iterations)
        {
            return Err(BinaryMiscibilityThermalError::Invalid(
                "source-supported temperature bracket, positive temperature accuracy, finite budgets and bounded root work required",
            ));
        }
        Ok(())
    }

    fn composition_change(
        &self,
        next: &BinaryMiscibilityState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64, f64)> {
        let model = self.network.model();
        let mut enthalpy = Sum::default();
        let mut entropy = Sum::default();
        let mut interface = Sum::default();
        for (cell, &capacity) in self.network.owner.capacities_mol.iter().enumerate() {
            poll(cancelled)?;
            let old = self.state.fractions[cell];
            let new = next.fractions[cell];
            enthalpy.add(finite(
                capacity
                    * model
                        .mixing_enthalpy_change_j_per_mol(old, new)
                        .map_err(BinaryMiscibilityError::from)?,
                "direct field mixing enthalpy increment",
            )?)?;
            entropy.add(finite(
                capacity
                    * model
                        .mixing_entropy_change_j_per_mol_k(old, new)
                        .map_err(BinaryMiscibilityError::from)?,
                "direct field mixing entropy increment",
            )?)?;
        }
        for face in self.network.owner.faces.iter() {
            poll(cancelled)?;
            let [i, j] = face.physical.cells;
            let old = self.state.fractions[i][0] - self.state.fractions[j][0];
            let new = next.fractions[i][0] - next.fractions[j][0];
            interface.add(finite(
                (face.gradient_j * (new - old)) * (0.5 * new + 0.5 * old),
                "direct interfacial internal-energy increment",
            )?)?;
        }
        Ok((enthalpy.value()?, entropy.value()?, interface.value()?))
    }

    fn sensible_change(
        &self,
        from: Temperature,
        to: Temperature,
        entropy: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        let mut result = Sum::default();
        for (amounts, source) in self
            .owner
            .component_moles
            .iter()
            .zip(self.caloric_model().heat_capacities())
        {
            poll(cancelled)?;
            let increment = if entropy {
                source.entropy_change_j_per_mol_k(from, to)?
            } else {
                source.enthalpy_change_j_per_mol(from, to)?
            };
            for &amount in amounts {
                result.add(finite(amount * increment, "conserved sensible increment")?)?;
            }
        }
        Ok(result.value()?)
    }

    fn solve_temperature(
        &self,
        from: Temperature,
        structural_energy: f64,
        controls: BinaryMiscibilityThermalSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Temperature, usize)> {
        let [mut lo, mut hi] = controls.temperature_bracket.map(|t| t.value());
        if structural_energy == 0.0 && (lo..=hi).contains(&from.value()) {
            poll(cancelled)?;
            return Ok((from, 0));
        }
        let residual = |temperature: f64, cancelled: &mut _| -> Result<f64> {
            total(&[
                structural_energy,
                self.sensible_change(from, Temperature::new(temperature), false, cancelled)?,
            ])
        };
        let endpoints = [residual(lo, cancelled)?, residual(hi, cancelled)?];
        if endpoints[0] > 0.0 || endpoints[1] < 0.0 {
            return Err(BinaryMiscibilityThermalError::HeatOutsideBracket {
                residuals_j: endpoints,
            });
        }
        for (temperature, value) in [(lo, endpoints[0]), (hi, endpoints[1])] {
            if value == 0.0 {
                return Ok((Temperature::new(temperature), 0));
            }
        }
        for iteration in 1..=controls.maximum_iterations {
            poll(cancelled)?;
            let midpoint = lo + 0.5 * (hi - lo);
            let value = residual(midpoint, cancelled)?;
            if value == 0.0
                || (value.abs() <= controls.energy_tolerance_j
                    && hi - lo <= controls.temperature_tolerance_k)
            {
                return Ok((Temperature::new(midpoint), iteration));
            }
            if midpoint == lo || midpoint == hi {
                return Err(BinaryMiscibilityThermalError::IterationLimit);
            }
            if value < 0.0 {
                lo = midpoint;
            } else {
                hi = midpoint;
            }
        }
        Err(BinaryMiscibilityThermalError::IterationLimit)
    }
}

fn total(values: &[f64]) -> Result<f64> {
    let mut sum = Sum::default();
    for &value in values {
        sum.add(value)?;
    }
    Ok(sum.value()?)
}

#[cfg(test)]
mod tests;
