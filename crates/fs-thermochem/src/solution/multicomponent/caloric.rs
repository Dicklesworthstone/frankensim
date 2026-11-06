//! Finite homogeneous liquid calorics from the original multicomponent law.
//!
//! Pure-component Cp and Cp/T use the shared source integrals. The same
//! A_ij+B_ij/T potential supplies excess H and S exactly once. At fixed
//! composition these excess terms are independent of temperature, so thermal
//! trials need only pure Cp integrals and their positive exact tangent. The
//! common safeguarded caloric inverse allocates no composition or phase state
//! per trial; the original mixing law is evaluated outside that iteration.
//!
//! This is a constrained homogeneous endpoint, not phase equilibrium or a
//! nucleation model. Adiabatic mixing with decreasing entropy refuses. For
//! nonzero prescribed heat, the reported system entropy change requires the
//! actual reservoir entropy before a second-law conclusion can be drawn.

/// Finite signed heat on an explicitly accepted local fixed-phase-count branch.
pub mod equilibrium;

use super::{
    MAX_REGULAR_SOLUTION_SPECIES_V1, MulticomponentRegularSolutionV1,
    MulticomponentSolutionErrorV1, MulticomponentSolutionStateV1, R, Sum,
};
use crate::solution::caloric::temperature::{self, TemperatureControl, TemperatureSample};
use crate::solution::caloric::{BinaryCaloricErrorV1, PureLiquidHeatCapacityV1};
use fs_math::det::ln;
use fs_qty::Temperature;
use std::sync::Arc;

/// Invalid source binding, unrepresentable inventory, or unsupported endpoint.
#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentCaloricErrorV1 {
    Invalid(&'static str),
    Arithmetic(&'static str),
    Caloric(BinaryCaloricErrorV1),
    Thermodynamics(MulticomponentSolutionErrorV1),
    /// An isolated homogeneous endpoint is not an admissible mixing process.
    EntropyDecrease {
        change_j_per_k: f64,
    },
    Cancelled,
}
impl core::fmt::Display for MulticomponentCaloricErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "homogeneous multicomponent liquid calorics: {self:?}")
    }
}
impl std::error::Error for MulticomponentCaloricErrorV1 {}
impl From<BinaryCaloricErrorV1> for MulticomponentCaloricErrorV1 {
    fn from(error: BinaryCaloricErrorV1) -> Self {
        Self::Caloric(error)
    }
}
impl From<MulticomponentSolutionErrorV1> for MulticomponentCaloricErrorV1 {
    fn from(error: MulticomponentSolutionErrorV1) -> Self {
        Self::Thermodynamics(error)
    }
}
type Result<T> = core::result::Result<T, MulticomponentCaloricErrorV1>;

/// Original-source temperature, energy and isolated-mixing entropy budgets.
/// Both root tolerances must be met unless a sampled residual is exactly zero.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentCaloricSolveV1 {
    pub temperature_bracket: [Temperature; 2],
    pub temperature_tolerance_k: f64,
    pub energy_tolerance_j: f64,
    pub entropy_tolerance_j_per_k: f64,
    /// Safeguarded Newton/bisection steps, 1..=256 after the two endpoints.
    pub maximum_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
struct Owner {
    thermodynamics: MulticomponentRegularSolutionV1,
    heat_capacity: Box<[PureLiquidHeatCapacityV1]>,
    temperature_domain: [Temperature; 2],
}

/// One immutable nonideal potential and explicitly qualified pure-liquid Cp.
/// Common pressure and pure-liquid source qualification belong to the caller.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentLiquidCaloricModelV1 {
    owner: Arc<Owner>,
}
impl MulticomponentLiquidCaloricModelV1 {
    /// Canonicalize sources by chemical identity, retaining every source knot.
    /// Every species needs exactly one source and the same sensible datum.
    pub fn new(
        thermodynamics: MulticomponentRegularSolutionV1,
        mut heat_capacity: Vec<PureLiquidHeatCapacityV1>,
    ) -> Result<Self> {
        heat_capacity.sort_by(|a, b| a.species().cmp(b.species()));
        if !heat_capacity
            .iter()
            .map(|source| source.species())
            .eq(thermodynamics.species())
            || heat_capacity
                .windows(2)
                .any(|pair| pair[0].reference_temperature() != pair[1].reference_temperature())
        {
            return Err(MulticomponentCaloricErrorV1::Invalid(
                "one pure-liquid Cp per original species and a common datum are required",
            ));
        }
        let mut domain = thermodynamics.temperature_domain();
        for source in &heat_capacity {
            let span = source.temperature_domain();
            domain[0] = Temperature::new(domain[0].value().max(span[0].value()));
            domain[1] = Temperature::new(domain[1].value().min(span[1].value()));
        }
        if domain[0].value() >= domain[1].value() {
            return Err(MulticomponentCaloricErrorV1::Invalid(
                "empty common caloric domain",
            ));
        }
        Ok(Self {
            owner: Arc::new(Owner {
                thermodynamics,
                heat_capacity: heat_capacity.into_boxed_slice(),
                temperature_domain: domain,
            }),
        })
    }

    pub fn thermodynamics(&self) -> &MulticomponentRegularSolutionV1 {
        &self.owner.thermodynamics
    }
    pub fn heat_capacities(&self) -> &[PureLiquidHeatCapacityV1] {
        &self.owner.heat_capacity
    }
    pub fn temperature_domain(&self) -> [Temperature; 2] {
        self.owner.temperature_domain
    }

    /// Prepare an explicitly homogeneous finite inventory, including metastable
    /// states. This evaluates no phase split and adds no absent-species floor.
    pub fn homogeneous(
        &self,
        temperature: Temperature,
        component_moles: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentCaloricStateV1> {
        poll(&mut cancelled)?;
        self.admit_temperature(temperature)?;
        let (total_moles, fractions) = self.inventory(component_moles)?;
        let mixing = self.thermodynamics().evaluate(
            temperature,
            &fractions[..component_moles.len()],
            &mut cancelled,
        )?;
        let mut enthalpy = Sum::default();
        for (&n, source) in component_moles.iter().zip(self.heat_capacities()) {
            poll(&mut cancelled)?;
            if n != 0.0 {
                enthalpy.add(scaled(
                    n,
                    source
                        .enthalpy_change_j_per_mol(source.reference_temperature(), temperature)?,
                    "pure-liquid enthalpy",
                )?);
            }
        }
        enthalpy.add(scaled(
            total_moles,
            mixing.mixing_enthalpy_j_per_mol(),
            "mixing enthalpy",
        )?);
        let state = MulticomponentCaloricStateV1 {
            owner: Arc::clone(&self.owner),
            mixing,
            component_moles: component_moles.into(),
            total_moles,
            enthalpy_j: checked(enthalpy.value(), "total liquid enthalpy")?,
            heat_capacity_j_per_k: self.capacity(temperature, component_moles, &mut cancelled)?,
        };
        poll(&mut cancelled)?;
        Ok(state)
    }

    /// Add signed heat at fixed pressure and component amounts. Composition
    /// remains homogeneous; changing equilibrium phases is a separate owner.
    pub fn add_heat(
        &self,
        initial: &MulticomponentCaloricStateV1,
        heat_input_j: f64,
        controls: MulticomponentCaloricSolveV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentCaloricBalanceV1> {
        self.solve(&[initial], heat_input_j, controls, cancelled)
    }

    /// Mix unequal-temperature parcels with one shared chemical/caloric owner.
    /// Sensible heat and nonideal mixing heat jointly determine the endpoint.
    /// With zero external heat, an entropy-decreasing homogeneous endpoint
    /// refuses instead of implying spontaneous miscibility or hidden work.
    pub fn mix(
        &self,
        left: &MulticomponentCaloricStateV1,
        right: &MulticomponentCaloricStateV1,
        heat_input_j: f64,
        controls: MulticomponentCaloricSolveV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentCaloricBalanceV1> {
        self.solve(&[left, right], heat_input_j, controls, cancelled)
    }

    /// Original-source enthalpy change for the same conserved inventory [J].
    pub fn enthalpy_change_j(
        &self,
        initial: &MulticomponentCaloricStateV1,
        final_state: &MulticomponentCaloricStateV1,
    ) -> Result<f64> {
        self.admit_pair(initial, final_state)?;
        self.sensible_change(&[initial], final_state.temperature(), false, &mut || false)
    }

    /// Original pure Cp/T entropy change for the same conserved inventory [J/K].
    /// The homogeneous A+B/T mixing entropy is constant on this thermal path.
    pub fn entropy_change_j_per_k(
        &self,
        initial: &MulticomponentCaloricStateV1,
        final_state: &MulticomponentCaloricStateV1,
    ) -> Result<f64> {
        self.admit_pair(initial, final_state)?;
        self.sensible_change(&[initial], final_state.temperature(), true, &mut || false)
    }

    fn solve(
        &self,
        inputs: &[&MulticomponentCaloricStateV1],
        heat_input_j: f64,
        controls: MulticomponentCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentCaloricBalanceV1> {
        poll(&mut cancelled)?;
        let [lo, hi] = controls.temperature_bracket.map(|t| t.value());
        if !heat_input_j.is_finite()
            || !lo.is_finite()
            || !hi.is_finite()
            || lo >= hi
            || lo < self.temperature_domain()[0].value()
            || hi > self.temperature_domain()[1].value()
            || !positive(controls.temperature_tolerance_k)
            || !nonnegative(controls.energy_tolerance_j)
            || !nonnegative(controls.entropy_tolerance_j_per_k)
            || !(1..=256).contains(&controls.maximum_iterations)
        {
            return Err(MulticomponentCaloricErrorV1::Invalid(
                "thermal solve controls or heat input",
            ));
        }
        let n = self.heat_capacities().len();
        let mut amounts = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        for input in inputs {
            self.admit(input)?;
            for (amount, &added) in amounts.iter_mut().zip(input.component_moles()) {
                let next = checked(*amount + added, "combined component amount")?;
                if *amount > 0.0 && added > 0.0 && (next == *amount || next == added) {
                    return Err(MulticomponentCaloricErrorV1::Arithmetic(
                        "positive incoming component lost in addition",
                    ));
                }
                *amount = next;
            }
        }
        let (_, fractions) = self.inventory(&amounts[..n])?;
        // Retain exactly the original owner's common normalization. Its H and S
        // are temperature-independent on this homogeneous composition path.
        let composition = self.thermodynamics().evaluate(
            inputs[0].temperature(),
            &fractions[..n],
            &mut cancelled,
        )?;
        let mixing_changes =
            self.mixing_changes(inputs, composition.mole_fractions(), &mut cancelled)?;
        let root = temperature::solve(
            TemperatureControl {
                bracket: [lo, hi],
                temperature_tolerance_k: controls.temperature_tolerance_k,
                energy_tolerance_j: controls.energy_tolerance_j,
                maximum_iterations: controls.maximum_iterations,
            },
            |t| -> Result<_> {
                poll(&mut cancelled)?;
                let temperature = Temperature::new(t);
                let mut change = Sum::default();
                change.add(mixing_changes[0]);
                change.add(self.sensible_change(inputs, temperature, false, &mut cancelled)?);
                let change_j = checked(change.value(), "enthalpy increment")?;
                Ok(TemperatureSample {
                    state: (
                        temperature,
                        self.capacity(temperature, &amounts[..n], &mut cancelled)?,
                    ),
                    change_j,
                    residual_j: checked(change_j - heat_input_j, "heat balance residual")?,
                })
            },
            |state| Some(state.1),
        )?;
        let final_state = self.homogeneous(root.sample.state.0, &amounts[..n], &mut cancelled)?;
        let mut entropy = Sum::default();
        entropy.add(mixing_changes[1]);
        entropy.add(self.sensible_change(
            inputs,
            final_state.temperature(),
            true,
            &mut cancelled,
        )?);
        let entropy_change_j_per_k = checked(entropy.value(), "system entropy change")?;
        if heat_input_j == 0.0 && entropy_change_j_per_k < -controls.entropy_tolerance_j_per_k {
            return Err(MulticomponentCaloricErrorV1::EntropyDecrease {
                change_j_per_k: entropy_change_j_per_k,
            });
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentCaloricBalanceV1 {
            final_state,
            heat_input_j,
            enthalpy_change_j: root.sample.change_j,
            entropy_change_j_per_k,
            energy_residual_j: root.sample.residual_j,
            temperature_bracket: root.bracket.map(Temperature::new),
            iterations: root.iterations,
        })
    }

    // Mixing changes are evaluated once per transaction, not per thermal trial.
    // Difference-of-products retains the original pair coefficients. The ideal
    // term is a weighted relative entropy; its linear completion cancels for
    // normalized inventories and retains dominant-species contributions when
    // their mole fraction rounds to one. Pure entropy datums cancel first.
    fn mixing_changes(
        &self,
        inputs: &[&MulticomponentCaloricStateV1],
        final_fractions: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2]> {
        let n = final_fractions.len();
        let mut enthalpy = Sum::default();
        let mut entropy = Sum::default();
        for input in inputs {
            let original = input.mixing.mole_fractions();
            let mut excess_h = Sum::default();
            let mut excess_s = Sum::default();
            for i in 0..n {
                poll(cancelled)?;
                entropy.add(scaled(
                    R,
                    relative_entropy_amount(
                        input.total_moles,
                        input.component_moles[i],
                        original[i],
                        final_fractions[i],
                    )?,
                    "ideal mixing entropy",
                )?);
                for j in i + 1..n {
                    let delta_i = final_fractions[i] - original[i];
                    let delta_j = final_fractions[j] - original[j];
                    let middle_i = 0.5 * final_fractions[i] + 0.5 * original[i];
                    let middle_j = 0.5 * final_fractions[j] + 0.5 * original[j];
                    let difference = |coefficient: f64| {
                        (coefficient * delta_i) * middle_j + (coefficient * delta_j) * middle_i
                    };
                    excess_h.add(difference(
                        self.thermodynamics().chi_inverse_temperature_k()[i * n + j],
                    ));
                    excess_s.add(-difference(
                        self.thermodynamics().chi_constants()[i * n + j],
                    ));
                }
            }
            enthalpy.add(scaled(
                input.total_moles,
                R * excess_h.value(),
                "excess mixing enthalpy",
            )?);
            entropy.add(scaled(
                input.total_moles,
                R * excess_s.value(),
                "excess mixing entropy",
            )?);
        }
        Ok([
            checked(enthalpy.value(), "mixing enthalpy change")?,
            checked(entropy.value(), "mixing entropy change")?,
        ])
    }

    fn sensible_change(
        &self,
        inputs: &[&MulticomponentCaloricStateV1],
        final_temperature: Temperature,
        entropy: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        let mut sum = Sum::default();
        for input in inputs {
            for (&n, source) in input.component_moles().iter().zip(self.heat_capacities()) {
                poll(cancelled)?;
                if n != 0.0 {
                    let change = if entropy {
                        source.entropy_change_j_per_mol_k(input.temperature(), final_temperature)?
                    } else {
                        source.enthalpy_change_j_per_mol(input.temperature(), final_temperature)?
                    };
                    sum.add(scaled(n, change, "pure-liquid caloric increment")?);
                }
            }
        }
        checked(sum.value(), "sensible caloric change")
    }

    fn capacity(
        &self,
        temperature: Temperature,
        amounts: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        let mut sum = Sum::default();
        for (&n, source) in amounts.iter().zip(self.heat_capacities()) {
            poll(cancelled)?;
            if n != 0.0 {
                sum.add(scaled(
                    n,
                    source.heat_capacity_at(temperature)?.value(),
                    "liquid heat capacity",
                )?);
            }
        }
        let capacity = checked(sum.value(), "total liquid heat capacity")?;
        if !positive(capacity) {
            return Err(MulticomponentCaloricErrorV1::Arithmetic(
                "positive liquid heat capacity",
            ));
        }
        Ok(capacity)
    }

    fn inventory(&self, amounts: &[f64]) -> Result<(f64, [f64; MAX_REGULAR_SOLUTION_SPECIES_V1])> {
        if amounts.len() != self.heat_capacities().len() || amounts.iter().any(|&n| !nonnegative(n))
        {
            return Err(MulticomponentCaloricErrorV1::Invalid(
                "finite nonnegative original component amounts",
            ));
        }
        let mut sum = Sum::default();
        for &n in amounts {
            sum.add(n);
        }
        let total = checked(sum.value(), "total component amount")?;
        if !positive(total) {
            return Err(MulticomponentCaloricErrorV1::Invalid(
                "positive total component amount required",
            ));
        }
        let mut fractions = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        for (x, &n) in fractions.iter_mut().zip(amounts) {
            *x = n / total;
            if n > 0.0 && *x == 0.0 {
                return Err(MulticomponentCaloricErrorV1::Arithmetic(
                    "positive trace lost in normalization",
                ));
            }
        }
        Ok((total, fractions))
    }

    fn admit_temperature(&self, temperature: Temperature) -> Result<()> {
        if !temperature.value().is_finite()
            || temperature.value() < self.temperature_domain()[0].value()
            || temperature.value() > self.temperature_domain()[1].value()
        {
            return Err(MulticomponentCaloricErrorV1::Invalid(
                "temperature outside common caloric domain",
            ));
        }
        Ok(())
    }
    fn admit(&self, state: &MulticomponentCaloricStateV1) -> Result<()> {
        if !Arc::ptr_eq(&self.owner, &state.owner) && self.owner != state.owner {
            return Err(MulticomponentCaloricErrorV1::Invalid(
                "parcel requires its original caloric owner",
            ));
        }
        Ok(())
    }
    fn admit_pair(
        &self,
        initial: &MulticomponentCaloricStateV1,
        final_state: &MulticomponentCaloricStateV1,
    ) -> Result<()> {
        self.admit(initial)?;
        self.admit(final_state)?;
        if initial.component_moles != final_state.component_moles {
            return Err(MulticomponentCaloricErrorV1::Invalid(
                "caloric increment requires conserved component amounts",
            ));
        }
        Ok(())
    }
}

/// Homogeneous finite inventory with original species and pure-liquid references.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentCaloricStateV1 {
    owner: Arc<Owner>,
    mixing: MulticomponentSolutionStateV1,
    component_moles: Box<[f64]>,
    total_moles: f64,
    enthalpy_j: f64,
    heat_capacity_j_per_k: f64,
}
impl MulticomponentCaloricStateV1 {
    pub const fn temperature(&self) -> Temperature {
        self.mixing.temperature()
    }
    pub fn component_moles(&self) -> &[f64] {
        &self.component_moles
    }
    pub const fn total_moles(&self) -> f64 {
        self.total_moles
    }
    pub const fn mixing(&self) -> &MulticomponentSolutionStateV1 {
        &self.mixing
    }
    /// Pure sensible enthalpy at the declared datum plus original excess H.
    pub const fn enthalpy_j(&self) -> f64 {
        self.enthalpy_j
    }
    /// Homogeneous isobaric capacity; the A+B/T excess H has no T derivative.
    pub const fn isobaric_heat_capacity_j_per_k(&self) -> f64 {
        self.heat_capacity_j_per_k
    }
}

/// Conserved-energy homogeneous endpoint and separately retained system entropy.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentCaloricBalanceV1 {
    final_state: MulticomponentCaloricStateV1,
    heat_input_j: f64,
    enthalpy_change_j: f64,
    entropy_change_j_per_k: f64,
    energy_residual_j: f64,
    temperature_bracket: [Temperature; 2],
    iterations: usize,
}
impl MulticomponentCaloricBalanceV1 {
    pub const fn final_state(&self) -> &MulticomponentCaloricStateV1 {
        &self.final_state
    }
    pub const fn heat_input_j(&self) -> f64 {
        self.heat_input_j
    }
    pub const fn enthalpy_change_j(&self) -> f64 {
        self.enthalpy_change_j
    }
    /// This is entropy production only when the transaction is closed/adiabatic.
    pub const fn entropy_change_j_per_k(&self) -> f64 {
        self.entropy_change_j_per_k
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

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn checked(value: f64, what: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MulticomponentCaloricErrorV1::Arithmetic(what))
    }
}
fn scaled(amount: f64, value: f64, what: &'static str) -> Result<f64> {
    let result = checked(amount * value, what)?;
    if result == 0.0 && amount != 0.0 && value != 0.0 {
        return Err(MulticomponentCaloricErrorV1::Arithmetic(
            "nonzero caloric contribution underflows",
        ));
    }
    Ok(result)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentCaloricErrorV1::Cancelled)
    } else {
        Ok(())
    }
}
fn log_ratio(a: f64, b: f64) -> f64 {
    if a == b {
        return 0.0;
    }
    let (lo, hi, sign) = if a > b { (b, a, 1.0) } else { (a, b, -1.0) };
    let delta = (hi - lo) / lo;
    let shifted = 1.0 + delta;
    sign * if shifted == 1.0 {
        delta
    } else if shifted.is_finite() {
        ln(shifted) * (delta / (shifted - 1.0))
    } else {
        ln(hi) - ln(lo)
    }
}

// N * [x ln(x/y) - x + y], without cancellation near x=y. The exact x=0
// limit is N*y. These positive terms give the ideal mixing entropy even when
// the dominant fraction rounds to one beside retained dilute components.
fn relative_entropy_amount(total: f64, amount: f64, x: f64, y: f64) -> Result<f64> {
    if x == y {
        return Ok(0.0);
    }
    if y == 0.0 {
        return Err(MulticomponentCaloricErrorV1::Arithmetic(
            "mixing lost a present component",
        ));
    }
    let reference = scaled(total, y, "ideal mixing reference amount")?;
    if x == 0.0 {
        return Ok(reference);
    }
    let r = (x - y) / y;
    if r.abs() < 0.125 {
        let mut polynomial = 0.0;
        for n in (2..=20).rev() {
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            polynomial = sign / f64::from(n * (n - 1)) + r * polynomial;
        }
        return scaled(reference, r * r * polynomial, "ideal relative entropy");
    }
    let mut value = Sum::default();
    value.add(scaled(amount, log_ratio(x, y), "ideal entropy logarithm")?);
    value.add(-amount);
    value.add(reference);
    let value = checked(value.value(), "ideal relative entropy")?;
    if value < 0.0 {
        return Err(MulticomponentCaloricErrorV1::Arithmetic(
            "negative ideal relative entropy",
        ));
    }
    Ok(value)
}

#[cfg(test)]
mod tests;
