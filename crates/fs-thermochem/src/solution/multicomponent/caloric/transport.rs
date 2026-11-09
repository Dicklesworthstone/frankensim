//! Original-source calorics for a finite open homogeneous liquid cell.
//!
//! Every component amount is independent. Signed exchange is accounted for
//! against the original pure-liquid datum, while close-state H/S differences
//! use source integrals and direct composition increments. No local entropy
//! gate is applied: an open cell may export entropy. The spatial owner must
//! audit the entropy and energy of its complete closed transaction.
//!
//! This first transport closure retains a common strictly positive species
//! support. It does not add concentration floors, change phases, or infer a
//! Soret coefficient from a temperature-dependent activity coefficient.

use super::*;

mod face;
pub use face::{MulticomponentLiquidTransportEffortsV1, MulticomponentLiquidTransportPointV1};

/// Original-source inverse controls for signed inventory and energy exchange.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentLiquidTransportSolveV1 {
    pub temperature_bracket: [Temperature; 2],
    pub temperature_tolerance_k: f64,
    pub energy_tolerance_j: f64,
    pub maximum_iterations: usize,
}

/// Signed changes in the cell, not the entropy production of an open cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentLiquidTransportChangeV1 {
    pub enthalpy_change_j: f64,
    pub entropy_change_j_per_k: f64,
}

/// An unpublished original-source endpoint of a signed exchange transaction.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentLiquidTransportBalanceV1 {
    final_state: MulticomponentCaloricStateV1,
    change: MulticomponentLiquidTransportChangeV1,
    energy_residual_j: f64,
    temperature_bracket: [Temperature; 2],
    iterations: usize,
}
impl MulticomponentLiquidTransportBalanceV1 {
    pub const fn final_state(&self) -> &MulticomponentCaloricStateV1 {
        &self.final_state
    }
    pub const fn change(&self) -> &MulticomponentLiquidTransportChangeV1 {
        &self.change
    }
    pub const fn energy_residual_j(&self) -> f64 {
        self.energy_residual_j
    }
    pub const fn temperature_bracket(&self) -> [Temperature; 2] {
        self.temperature_bracket
    }
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
}

impl MulticomponentLiquidCaloricModelV1 {
    /// Account for changed independent component amounts and temperature.
    /// Both states must retain this original caloric owner and positive support.
    pub fn homogeneous_transport_change(
        &self,
        initial: &MulticomponentCaloricStateV1,
        final_state: &MulticomponentCaloricStateV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidTransportChangeV1> {
        poll(&mut cancelled)?;
        self.admit(initial)?;
        self.admit(final_state)?;
        positive_support(initial.component_moles())?;
        positive_support(final_state.component_moles())?;
        let baseline = self.transport_baseline(
            initial,
            final_state.component_moles(),
            final_state.total_moles(),
            final_state.mixing().mole_fractions(),
            &mut cancelled,
        )?;
        let enthalpy_change_j = self.transport_increment(
            initial.temperature(),
            final_state.temperature(),
            final_state.component_moles(),
            baseline[0],
            false,
            &mut cancelled,
        )?;
        let entropy_change_j_per_k = self.transport_increment(
            initial.temperature(),
            final_state.temperature(),
            final_state.component_moles(),
            baseline[1],
            true,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(MulticomponentLiquidTransportChangeV1 {
            enthalpy_change_j,
            entropy_change_j_per_k,
        })
    }

    /// Recover the temperature after signed species and enthalpy exchange.
    /// The original composition law is evaluated outside the scalar inverse;
    /// each root sample uses O(N) source integrals and the exact positive Cp.
    /// Refusal and cancellation leave the input untouched.
    pub fn homogeneous_transport_candidate(
        &self,
        initial: &MulticomponentCaloricStateV1,
        component_moles: &[f64],
        enthalpy_change_j: f64,
        controls: MulticomponentLiquidTransportSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidTransportBalanceV1> {
        poll(&mut cancelled)?;
        self.admit(initial)?;
        positive_support(initial.component_moles())?;
        positive_support(component_moles)?;
        let bracket = controls.temperature_bracket.map(|t| t.value());
        for t in controls.temperature_bracket {
            self.admit_temperature(t)?;
        }
        checked(enthalpy_change_j, "signed transported enthalpy")?;
        let thermal = TemperatureControl {
            bracket,
            temperature_tolerance_k: controls.temperature_tolerance_k,
            energy_tolerance_j: controls.energy_tolerance_j,
            maximum_iterations: controls.maximum_iterations,
        };
        thermal.validate()?;
        let (total, fractions) = self.inventory(component_moles)?;
        let composition = self.thermodynamics().evaluate(
            initial.temperature(),
            &fractions[..component_moles.len()],
            &mut cancelled,
        )?;
        let baseline = self.transport_baseline(
            initial,
            component_moles,
            total,
            composition.mole_fractions(),
            &mut cancelled,
        )?;
        let root = temperature::solve(
            thermal,
            |t| -> Result<_> {
                poll(&mut cancelled)?;
                let t = Temperature::new(t);
                let change_j = self.transport_increment(
                    initial.temperature(),
                    t,
                    component_moles,
                    baseline[0],
                    false,
                    &mut cancelled,
                )?;
                Ok(TemperatureSample {
                    state: (t, self.capacity(t, component_moles, &mut cancelled)?),
                    change_j,
                    residual_j: checked(
                        change_j - enthalpy_change_j,
                        "signed transport energy residual",
                    )?,
                })
            },
            |state| Some(state.1),
        )?;
        let final_state = self.homogeneous(root.sample.state.0, component_moles, &mut cancelled)?;
        let entropy_change_j_per_k = self.transport_increment(
            initial.temperature(),
            final_state.temperature(),
            component_moles,
            baseline[1],
            true,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(MulticomponentLiquidTransportBalanceV1 {
            final_state,
            change: MulticomponentLiquidTransportChangeV1 {
                enthalpy_change_j: root.sample.change_j,
                entropy_change_j_per_k,
            },
            energy_residual_j: root.sample.residual_j,
            temperature_bracket: root.bracket.map(Temperature::new),
            iterations: root.iterations,
        })
    }

    fn transport_increment(
        &self,
        initial_temperature: Temperature,
        final_temperature: Temperature,
        amounts: &[f64],
        baseline: f64,
        entropy: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        let mut sum = Sum::default();
        sum.add(baseline);
        sum.add(self.sensible_inventory_change(
            [(initial_temperature, amounts)],
            final_temperature,
            entropy,
            cancelled,
        )?);
        checked(sum.value(), "signed source caloric increment")
    }

    // Amount changes at the initial temperature plus changes in the original
    // temperature-independent excess H/S. Pure datum terms are retained: an
    // open cell cannot cancel the references of exported/imported species.
    fn transport_baseline(
        &self,
        initial: &MulticomponentCaloricStateV1,
        amounts: &[f64],
        total: f64,
        fractions: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2]> {
        let old_amounts = initial.component_moles();
        let old_x = initial.mixing().mole_fractions();
        let n = amounts.len();
        let mut enthalpy = Sum::default();
        let mut entropy = Sum::default();
        let mut ideal = Sum::default();
        for i in 0..n {
            poll(cancelled)?;
            let delta = amounts[i] - old_amounts[i];
            let source = &self.heat_capacities()[i];
            enthalpy.add(scaled(
                delta,
                source.enthalpy_change_j_per_mol(
                    source.reference_temperature(),
                    initial.temperature(),
                )?,
                "signed pure enthalpy datum",
            )?);
            entropy.add(scaled(
                delta,
                source.entropy_change_j_per_mol_k(
                    source.reference_temperature(),
                    initial.temperature(),
                )?,
                "signed pure entropy datum",
            )?);
            ideal.add(scaled(
                delta,
                original_log(old_x, i),
                "ideal entropy linear term",
            )?);
            ideal.add(amount_relative_entropy(amounts[i], old_amounts[i])?);
            for j in i + 1..n {
                // Symmetric difference of n_i*x_j; no subtraction of large
                // total mixing enthalpies or multiplication of two traces first.
                let dx = fractions[j] - old_x[j];
                let average_x = 0.5 * (fractions[j] + old_x[j]);
                let average_n = 0.5 * amounts[i] + 0.5 * old_amounts[i];
                let pair_change = |coefficient: f64| -> Result<f64> {
                    let mut sum = Sum::default();
                    sum.add(scaled(
                        coefficient * delta,
                        average_x,
                        "pair amount increment",
                    )?);
                    sum.add(scaled(
                        coefficient * dx,
                        average_n,
                        "pair fraction increment",
                    )?);
                    checked(sum.value(), "pair mixing increment")
                };
                enthalpy.add(pair_change(
                    R * self.thermodynamics().chi_inverse_temperature_k()[i * n + j],
                )?);
                entropy.add(-pair_change(
                    R * self.thermodynamics().chi_constants()[i * n + j],
                )?);
            }
        }
        ideal.add(-amount_relative_entropy(total, initial.total_moles())?);
        entropy.add(scaled(-R, ideal.value(), "signed ideal mixing entropy")?);
        Ok([
            checked(enthalpy.value(), "signed baseline enthalpy")?,
            checked(entropy.value(), "signed baseline entropy")?,
        ])
    }
}

fn positive_support(amounts: &[f64]) -> Result<()> {
    if amounts.is_empty() || amounts.iter().any(|&n| !positive(n)) {
        Err(MulticomponentCaloricErrorV1::Invalid(
            "transport requires a common strictly positive species support",
        ))
    } else {
        Ok(())
    }
}

// a ln(a/b)-a+b. The direct amount formulation avoids forming squared
// normalized traces. Its quadratic limit retains tiny finite relaxation.
fn amount_relative_entropy(a: f64, b: f64) -> Result<f64> {
    if a == b {
        return Ok(0.0);
    }
    let r = (a - b) / b;
    if r.abs() < 0.125 {
        let mut polynomial = 0.0;
        for k in (2..=20).rev() {
            let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
            polynomial = sign / f64::from(k * (k - 1)) + r * polynomial;
        }
        return checked(((b * r) * r) * polynomial, "ideal entropy remainder");
    }
    let mut sum = Sum::default();
    sum.add(scaled(a, log_ratio(a, b), "ideal entropy logarithm")?);
    sum.add(b - a);
    checked(sum.value(), "ideal entropy remainder")
}

fn original_log(fractions: &[f64], species: usize) -> f64 {
    if fractions[species] <= 0.5 {
        ln(fractions[species])
    } else {
        let mut others = Sum::default();
        for (i, &x) in fractions.iter().enumerate() {
            if i != species {
                others.add(x);
            }
        }
        crate::solution::ln_complement(others.value())
    }
}

#[cfg(test)]
mod tests;
