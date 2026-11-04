//! Closed, isobaric binary-liquid mixing and heating with one caloric owner.
//!
//! Pure-liquid Cp(T) is positive piecewise linear source data. Its integral
//! supplies sensible enthalpy; the original pair law supplies excess enthalpy
//! exactly once. Every temperature trial resolves the actual liquid phases.
//! A total-mixture heat capacity or gas standard state is not a pure-liquid
//! caloric reference. Pressure/chemical compatibility remains caller-supplied.
//!
//! Heat balances use direct enthalpy INCREMENTS from each original parcel,
//! so a large arbitrary reference enthalpy cannot erase a small heat input.
//! No vapor, solid, spatial transport, nucleation, or mechanical work beyond
//! constant-pressure boundary work is included.

use super::{
    BinaryCoexistenceOptionsV1, BinaryPhaseEquilibriumV1, BinaryRegularSolutionV1,
    BinarySolutionErrorV1, BinarySolutionStateV1,
};
use crate::{MolarThermalQuantityV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_qty::Temperature;
use std::sync::Arc;

/// Maximum source points per pure-component heat-capacity curve.
pub const MAX_LIQUID_CALORIC_KNOTS_V1: usize = 4096;

/// Invalid caloric data, unsupported heat target, or incomplete thermal solve.
#[derive(Clone, Debug, PartialEq)]
pub enum BinaryCaloricErrorV1 {
    Invalid(&'static str),
    Arithmetic(&'static str),
    Thermodynamics(BinarySolutionErrorV1),
    /// Signed energy residuals at the offered lower/upper temperatures [J].
    HeatOutsideBracket {
        residuals_j: [f64; 2],
    },
    IterationLimit,
    Cancelled,
}
impl core::fmt::Display for BinaryCaloricErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "binary liquid calorics: {self:?}")
    }
}
impl std::error::Error for BinaryCaloricErrorV1 {}
impl From<BinarySolutionErrorV1> for BinaryCaloricErrorV1 {
    fn from(error: BinarySolutionErrorV1) -> Self {
        Self::Thermodynamics(error)
    }
}

/// Explicit pure-liquid Cp(T), in J/(mol K), with no excess-mixture contribution.
/// Linear interpolation preserves positive Cp between the admitted knots.
/// Source identity, pressure and the pure-liquid qualification belong to the
/// caller; this numerical law neither infers nor authenticates them.
#[derive(Clone, Debug, PartialEq)]
pub struct PureLiquidHeatCapacityV1 {
    species: SpeciesId,
    reference_temperature: Temperature,
    knots: Vec<(f64, f64)>,
}
impl PureLiquidHeatCapacityV1 {
    /// Retain 2..=4096 strictly increasing positive temperature knots. All Cp
    /// values must be positive and finite. The datum must lie inside the curve.
    pub fn new(
        species: SpeciesId,
        reference_temperature: Temperature,
        knots: Vec<(Temperature, MolarThermalQuantityV1)>,
    ) -> Result<Self, BinaryCaloricErrorV1> {
        if !(2..=MAX_LIQUID_CALORIC_KNOTS_V1).contains(&knots.len())
            || knots
                .iter()
                .any(|(t, cp)| !positive(t.value()) || !positive(cp.value()))
            || knots.windows(2).any(|w| w[0].0.value() >= w[1].0.value())
            || !reference_temperature.value().is_finite()
            || reference_temperature.value() < knots[0].0.value()
            || reference_temperature.value() > knots[knots.len() - 1].0.value()
        {
            return Err(BinaryCaloricErrorV1::Invalid(
                "pure-liquid Cp curve or datum",
            ));
        }
        Ok(Self {
            species,
            reference_temperature,
            knots: knots
                .into_iter()
                .map(|(t, cp)| (t.value(), cp.value()))
                .collect(),
        })
    }

    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }
    pub const fn reference_temperature(&self) -> Temperature {
        self.reference_temperature
    }
    pub fn temperature_domain(&self) -> [Temperature; 2] {
        [
            Temperature::new(self.knots[0].0),
            Temperature::new(self.knots[self.knots.len() - 1].0),
        ]
    }

    /// Exact integral of the declared piecewise-linear Cp interpolation [J/mol].
    /// Integrates only crossed segments, using their local temperature spans;
    /// cooling is the exact negative of the same ascending integration.
    pub fn enthalpy_change_j_per_mol(
        &self,
        from: Temperature,
        to: Temperature,
    ) -> Result<f64, BinaryCaloricErrorV1> {
        self.admit(from.value())?;
        self.admit(to.value())?;
        let (lo, hi, sign) = if from.value() <= to.value() {
            (from.value(), to.value(), 1.0)
        } else {
            (to.value(), from.value(), -1.0)
        };
        if lo == hi {
            return Ok(0.0);
        }
        let start = self
            .knots
            .partition_point(|&(t, _)| t <= lo)
            .saturating_sub(1);
        let mut sum = Sum::default();
        for pair in self.knots[start..].windows(2) {
            let a = lo.max(pair[0].0);
            let b = hi.min(pair[1].0);
            if b > a {
                let cp_a = interpolate(pair, a);
                let cp_b = interpolate(pair, b);
                sum.add(checked(
                    (0.5 * cp_a + 0.5 * cp_b) * (b - a),
                    "sensible increment",
                )?)?;
            }
            if pair[1].0 >= hi {
                break;
            }
        }
        checked(sign * sum.value()?, "sensible integral")
    }

    fn admit(&self, t: f64) -> Result<(), BinaryCaloricErrorV1> {
        if !t.is_finite() || t < self.knots[0].0 || t > self.knots[self.knots.len() - 1].0 {
            Err(BinaryCaloricErrorV1::Invalid(
                "temperature outside pure-liquid Cp source",
            ))
        } else {
            Ok(())
        }
    }
}

/// Explicit temperature and energy acceptance budgets. Both must be met unless
/// a bracket endpoint is an exact energy root. The bracket must be inside all
/// original source domains; failure never clamps an unsupported heat input.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryCaloricSolveV1 {
    pub temperature_bracket: [Temperature; 2],
    pub temperature_tolerance_k: f64,
    pub energy_tolerance_j: f64,
    /// Maximum temperature bisections, 1..=256; each has a bounded phase solve.
    pub maximum_iterations: usize,
    pub coexistence: BinaryCoexistenceOptionsV1,
}

#[derive(Clone, Debug, PartialEq)]
struct CaloricOwner {
    thermodynamics: BinaryRegularSolutionV1,
    heat_capacity: [PureLiquidHeatCapacityV1; 2],
    temperature_domain: [Temperature; 2],
}

/// Immutable binary pair plus ordered, compatible pure-component calorics.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryLiquidCaloricModelV1 {
    owner: Arc<CaloricOwner>,
}

impl BinaryLiquidCaloricModelV1 {
    /// Both species must match the pair's order, and both curves must use the
    /// same datum. No reference offsets or total-mixture Cp are added implicitly.
    pub fn new(
        thermodynamics: BinaryRegularSolutionV1,
        heat_capacity: [PureLiquidHeatCapacityV1; 2],
    ) -> Result<Self, BinaryCaloricErrorV1> {
        if heat_capacity[0].species() != &thermodynamics.species()[0]
            || heat_capacity[1].species() != &thermodynamics.species()[1]
            || heat_capacity[0].reference_temperature != heat_capacity[1].reference_temperature
        {
            return Err(BinaryCaloricErrorV1::Invalid(
                "ordered pure-liquid species or common caloric datum",
            ));
        }
        let mut domain = thermodynamics.temperature_domain();
        for source in &heat_capacity {
            let span = source.temperature_domain();
            domain[0] = Temperature::new(domain[0].value().max(span[0].value()));
            domain[1] = Temperature::new(domain[1].value().min(span[1].value()));
        }
        if domain[0].value() >= domain[1].value() {
            return Err(BinaryCaloricErrorV1::Invalid(
                "empty common temperature domain",
            ));
        }
        Ok(Self {
            owner: Arc::new(CaloricOwner {
                thermodynamics,
                heat_capacity,
                temperature_domain: domain,
            }),
        })
    }

    pub fn thermodynamics(&self) -> &BinaryRegularSolutionV1 {
        &self.owner.thermodynamics
    }
    pub fn heat_capacities(&self) -> &[PureLiquidHeatCapacityV1; 2] {
        &self.owner.heat_capacity
    }
    pub fn temperature_domain(&self) -> [Temperature; 2] {
        self.owner.temperature_domain
    }

    /// Prepare an equilibrated finite inventory at an externally prescribed T.
    pub fn equilibrate(
        &self,
        temperature: Temperature,
        component_moles: [f64; 2],
        options: BinaryCoexistenceOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryCaloricStateV1, BinaryCaloricErrorV1> {
        poll(&mut cancelled)?;
        let (total, fractions) = inventory(component_moles)?;
        let equilibrium = self.thermodynamics().equilibrate_composition(
            temperature,
            fractions,
            options,
            &mut cancelled,
        )?;
        let state = self.make_state(
            temperature,
            component_moles,
            total,
            *equilibrium.phases()[0].state(),
            Some(equilibrium),
        )?;
        poll(&mut cancelled)?;
        Ok(state)
    }

    /// Prepare an explicitly homogeneous, possibly unstable initial parcel.
    /// This does not claim phase equilibrium. `add_heat(..., 0, ...)` can then
    /// resolve adiabatic demixing, including its temperature change.
    pub fn homogeneous(
        &self,
        temperature: Temperature,
        component_moles: [f64; 2],
    ) -> Result<BinaryCaloricStateV1, BinaryCaloricErrorV1> {
        let (total, fractions) = inventory(component_moles)?;
        let state = self
            .thermodynamics()
            .evaluate_composition(temperature, fractions)?;
        self.make_state(temperature, component_moles, total, state, None)
    }

    fn make_state(
        &self,
        temperature: Temperature,
        component_moles: [f64; 2],
        total_moles: f64,
        mixing: BinarySolutionStateV1,
        equilibrium: Option<BinaryPhaseEquilibriumV1>,
    ) -> Result<BinaryCaloricStateV1, BinaryCaloricErrorV1> {
        let mut h = Sum::default();
        for (n, source) in component_moles.iter().zip(&self.owner.heat_capacity) {
            h.add(checked(
                n * source.enthalpy_change_j_per_mol(source.reference_temperature, temperature)?,
                "sensible enthalpy",
            )?)?;
        }
        // Both symmetric coexistence phases have the same x_A*x_B, so one
        // representative gives the exact same excess without lever roundoff.
        h.add(checked(
            total_moles * mixing.mixing_enthalpy_j_per_mol(),
            "mixing enthalpy",
        )?)?;
        Ok(BinaryCaloricStateV1 {
            owner: Arc::clone(&self.owner),
            temperature,
            component_moles,
            total_moles,
            mixing,
            equilibrium,
            enthalpy_j: h.value()?,
        })
    }

    /// Close heat input (positive into the liquids) and both species inventories
    /// when two original parcels mix. Unequal temperatures and a final one/two
    /// phase transition use the same energy account. Inputs remain unchanged.
    pub fn mix(
        &self,
        left: &BinaryCaloricStateV1,
        right: &BinaryCaloricStateV1,
        heat_input_j: f64,
        controls: BinaryCaloricSolveV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryCaloricBalanceV1, BinaryCaloricErrorV1> {
        self.solve(&[left, right], heat_input_j, controls, cancelled)
    }

    /// Heat, cool or adiabatically equilibrate one original finite inventory.
    pub fn add_heat(
        &self,
        initial: &BinaryCaloricStateV1,
        heat_input_j: f64,
        controls: BinaryCaloricSolveV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryCaloricBalanceV1, BinaryCaloricErrorV1> {
        self.solve(&[initial], heat_input_j, controls, cancelled)
    }

    fn solve(
        &self,
        inputs: &[&BinaryCaloricStateV1],
        heat_input_j: f64,
        controls: BinaryCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryCaloricBalanceV1, BinaryCaloricErrorV1> {
        poll(&mut cancelled)?;
        controls.coexistence.validate()?;
        let [mut lo, mut hi] = controls.temperature_bracket.map(|t| t.value());
        if !heat_input_j.is_finite()
            || !lo.is_finite()
            || !hi.is_finite()
            || lo >= hi
            || lo < self.owner.temperature_domain[0].value()
            || hi > self.owner.temperature_domain[1].value()
            || !positive(controls.temperature_tolerance_k)
            || !controls.energy_tolerance_j.is_finite()
            || controls.energy_tolerance_j < 0.0
            || !(1..=256).contains(&controls.maximum_iterations)
        {
            return Err(BinaryCaloricErrorV1::Invalid(
                "thermal solve controls or heat input",
            ));
        }
        let mut amounts = [0.0; 2];
        for input in inputs {
            if !Arc::ptr_eq(&self.owner, &input.owner) && self.owner != input.owner {
                return Err(BinaryCaloricErrorV1::Invalid(
                    "parcel belongs to another caloric pair",
                ));
            }
            for (n, added) in amounts.iter_mut().zip(input.component_moles) {
                let next = checked(*n + added, "combined component inventory")?;
                if *n > 0.0 && added > 0.0 && (next == *n || next == added) {
                    return Err(BinaryCaloricErrorV1::Arithmetic(
                        "positive incoming component lost in addition",
                    ));
                }
                *n = next;
            }
        }
        let mut sample = |t| -> Result<_, BinaryCaloricErrorV1> {
            poll(&mut cancelled)?;
            let state = self.equilibrate(
                Temperature::new(t),
                amounts,
                controls.coexistence,
                &mut cancelled,
            )?;
            let change = self.enthalpy_increment(inputs, &state)?;
            let residual = checked(change - heat_input_j, "heat balance residual")?;
            Ok((state, change, residual))
        };
        let mut lower = sample(lo)?;
        let mut upper = sample(hi)?;
        if lower.2 > 0.0 || upper.2 < 0.0 {
            return Err(BinaryCaloricErrorV1::HeatOutsideBracket {
                residuals_j: [lower.2, upper.2],
            });
        }
        for iterations in 0..=controls.maximum_iterations {
            let best = if lower.2.abs() <= upper.2.abs() {
                &lower
            } else {
                &upper
            };
            if best.2 == 0.0
                || ((hi - lo) <= controls.temperature_tolerance_k
                    && best.2.abs() <= controls.energy_tolerance_j)
            {
                let result = BinaryCaloricBalanceV1 {
                    final_state: best.0.clone(),
                    heat_input_j,
                    enthalpy_change_j: best.1,
                    energy_residual_j: best.2,
                    temperature_bracket: [Temperature::new(lo), Temperature::new(hi)],
                    iterations,
                };
                poll(&mut cancelled)?;
                return Ok(result);
            }
            if iterations == controls.maximum_iterations {
                break;
            }
            let midpoint = 0.5 * lo + 0.5 * hi;
            if midpoint == lo || midpoint == hi {
                break;
            }
            let middle = sample(midpoint)?;
            if middle.2 <= 0.0 {
                lo = midpoint;
                lower = middle;
            } else {
                hi = midpoint;
                upper = middle;
            }
        }
        Err(BinaryCaloricErrorV1::IterationLimit)
    }

    fn enthalpy_increment(
        &self,
        inputs: &[&BinaryCaloricStateV1],
        final_state: &BinaryCaloricStateV1,
    ) -> Result<f64, BinaryCaloricErrorV1> {
        let mut sum = Sum::default();
        for input in inputs {
            for (n, source) in input.component_moles.iter().zip(&self.owner.heat_capacity) {
                sum.add(checked(
                    n * source
                        .enthalpy_change_j_per_mol(input.temperature, final_state.temperature)?,
                    "parcel sensible increment",
                )?)?;
            }
            // Difference of products, without subtracting two extensive
            // reference enthalpies. Sorting retains symmetry under A/B exchange.
            let mut old = input.mixing.mole_fractions();
            let mut new = final_state.mixing.mole_fractions();
            if old[0] > old[1] {
                old.swap(0, 1);
            }
            if new[0] > new[1] {
                new.swap(0, 1);
            }
            let difference = (new[0] - old[0]) * (0.5 * new[1] + 0.5 * old[1])
                + (new[1] - old[1]) * (0.5 * new[0] + 0.5 * old[0]);
            let excess = (self.owner.thermodynamics.chi_inverse_temperature_k * difference) * R;
            sum.add(checked(
                input.total_moles * excess,
                "parcel mixing increment",
            )?)?;
        }
        sum.value()
    }
}

/// Finite immutable parcel tied to the exact pair law and pure caloric data.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryCaloricStateV1 {
    owner: Arc<CaloricOwner>,
    temperature: Temperature,
    component_moles: [f64; 2],
    total_moles: f64,
    mixing: BinarySolutionStateV1,
    equilibrium: Option<BinaryPhaseEquilibriumV1>,
    enthalpy_j: f64,
}
impl BinaryCaloricStateV1 {
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    pub const fn component_moles(&self) -> [f64; 2] {
        self.component_moles
    }
    /// Total enthalpy relative to zero pure-component enthalpies at the datum.
    /// Heat transactions use local increments rather than differences of this value.
    pub const fn enthalpy_j(&self) -> f64 {
        self.enthalpy_j
    }
    /// `None` marks an explicitly homogeneous, unequilibrated initial parcel.
    pub const fn equilibrium(&self) -> Option<&BinaryPhaseEquilibriumV1> {
        self.equilibrium.as_ref()
    }
    pub fn thermodynamics(&self) -> &BinaryRegularSolutionV1 {
        &self.owner.thermodynamics
    }
}

/// Completed energy-balanced equilibrium; no inputs are mutated on refusal.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryCaloricBalanceV1 {
    final_state: BinaryCaloricStateV1,
    heat_input_j: f64,
    enthalpy_change_j: f64,
    energy_residual_j: f64,
    temperature_bracket: [Temperature; 2],
    iterations: usize,
}
impl BinaryCaloricBalanceV1 {
    pub const fn final_state(&self) -> &BinaryCaloricStateV1 {
        &self.final_state
    }
    pub const fn heat_input_j(&self) -> f64 {
        self.heat_input_j
    }
    pub const fn enthalpy_change_j(&self) -> f64 {
        self.enthalpy_change_j
    }
    pub const fn energy_residual_j(&self) -> f64 {
        self.energy_residual_j
    }
    /// Numerical root bracket, not a rigorous interval certificate.
    pub const fn temperature_bracket(&self) -> [Temperature; 2] {
        self.temperature_bracket
    }
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
}

fn inventory(n: [f64; 2]) -> Result<(f64, [f64; 2]), BinaryCaloricErrorV1> {
    if n.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(BinaryCaloricErrorV1::Invalid(
            "finite nonnegative component moles required",
        ));
    }
    let total = n[0] + n[1];
    if !positive(total) {
        return Err(BinaryCaloricErrorV1::Invalid(
            "positive finite total moles required",
        ));
    }
    let fractions = [n[0] / total, n[1] / total];
    if (0..2).any(|i| n[i] > 0.0 && fractions[i] == 0.0) {
        return Err(BinaryCaloricErrorV1::Arithmetic(
            "trace component lost in normalization",
        ));
    }
    Ok((total, fractions))
}
fn interpolate(pair: &[(f64, f64)], t: f64) -> f64 {
    let weight = (t - pair[0].0) / (pair[1].0 - pair[0].0);
    (1.0 - weight) * pair[0].1 + weight * pair[1].1
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn checked(x: f64, what: &'static str) -> Result<f64, BinaryCaloricErrorV1> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(BinaryCaloricErrorV1::Arithmetic(what))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), BinaryCaloricErrorV1> {
    if cancelled() {
        Err(BinaryCaloricErrorV1::Cancelled)
    } else {
        Ok(())
    }
}
#[derive(Default)]
struct Sum {
    main: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), BinaryCaloricErrorV1> {
        let next = checked(self.main + value, "enthalpy sum")?;
        self.correction += if self.main.abs() >= value.abs() {
            (self.main - next) + value
        } else {
            (value - next) + self.main
        };
        self.main = next;
        Ok(())
    }
    fn value(&self) -> Result<f64, BinaryCaloricErrorV1> {
        checked(self.main + self.correction, "enthalpy sum correction")
    }
}

#[cfg(test)]
mod tests;
