//! Closed, isobaric dissolution and precipitation in one finite thermal bath.
//!
//! The original ideal solubility law owns saturation and its constant molar
//! dissolution enthalpy. Each component explicitly declares ONE positive Cp(T)
//! curve shared by its dissolved and pure-solid phases. Their sensible terms
//! then cancel in the phase difference, so h_dissolved - h_solid is exactly the
//! same constant delta H used by van't Hoff saturation. Unequal phase capacities
//! require a different, temperature-dependent thermodynamic law.
//!
//! Solvent and every finite component contribute sensible heat exactly once.
//! Dissolution can cool the bath, precipitation can warm it, and solid exhaustion
//! changes the equilibrium branch without dropping inventory or latent energy.
//! The existing caloric Cp integrator supplies all sensible increments.
//!
//! Assumptions are those of the original ideal, independent dilute-solute law:
//! incompressible solvent, negligible solute/solid volume and no cross-solute
//! interactions, ions, reactions, nucleation, transport or variable-pressure
//! work. Solvent molar volume and phase-qualified caloric evidence are explicit
//! caller inputs. Numerical coefficients do not establish empirical authority.

use super::{IdealSolubilityV1, SolubilityErrorV1, SolubilityPartitionV1};
use crate::solution::caloric::{
    BinaryCaloricErrorV1, MolarHeatCapacityCurveV1, PureLiquidHeatCapacityV1,
};
use fs_qty::{Amount, Temperature, Volume};
use std::sync::Arc;

/// Hard bound on independently soluble finite components in one thermal bath.
pub const MAX_SOLUBILITY_CALORIC_COMPONENTS_V1: usize = 64;

/// Refusal never changes an input inventory, temperature or phase stock.
#[derive(Clone, Debug, PartialEq)]
pub enum SolubilityCaloricErrorV1 {
    Invalid(&'static str),
    Arithmetic(&'static str),
    Solubility(SolubilityErrorV1),
    Caloric(BinaryCaloricErrorV1),
    /// Signed enthalpy-change-minus-input-heat residuals at the endpoints [J].
    HeatOutsideBracket {
        residuals_j: [f64; 2],
    },
    IterationLimit,
    Cancelled,
}
impl core::fmt::Display for SolubilityCaloricErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite solubility calorics: {self:?}")
    }
}
impl std::error::Error for SolubilityCaloricErrorV1 {}
impl From<SolubilityErrorV1> for SolubilityCaloricErrorV1 {
    fn from(error: SolubilityErrorV1) -> Self {
        Self::Solubility(error)
    }
}
impl From<BinaryCaloricErrorV1> for SolubilityCaloricErrorV1 {
    fn from(error: BinaryCaloricErrorV1) -> Self {
        Self::Caloric(error)
    }
}

/// Original saturation law and a declared COMMON dissolved/solid molar Cp(T).
/// This explicit equality is required by the law's constant dissolution enthalpy;
/// a pure-liquid heat-capacity wrapper is not accepted as solid-phase evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct SolubilityCaloricComponentV1 {
    law: IdealSolubilityV1,
    common_phase_cp: MolarHeatCapacityCurveV1,
}
impl SolubilityCaloricComponentV1 {
    pub fn new(
        law: IdealSolubilityV1,
        common_phase_cp: MolarHeatCapacityCurveV1,
    ) -> Result<Self, SolubilityCaloricErrorV1> {
        let law_domain = law.temperature_domain();
        let cp_domain = common_phase_cp.temperature_domain();
        if law.solute() != common_phase_cp.species()
            || law_domain[0].value().max(cp_domain[0].value())
                >= law_domain[1].value().min(cp_domain[1].value())
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "common phase Cp species or temperature domain",
            ));
        }
        Ok(Self {
            law,
            common_phase_cp,
        })
    }
    pub const fn law(&self) -> &IdealSolubilityV1 {
        &self.law
    }
    pub const fn common_phase_heat_capacity(&self) -> &MolarHeatCapacityCurveV1 {
        &self.common_phase_cp
    }
}

#[derive(Clone, Debug, PartialEq)]
struct CaloricOwner {
    solvent_molar_volume_m3_mol: f64,
    solvent_cp: PureLiquidHeatCapacityV1,
    components: Vec<SolubilityCaloricComponentV1>,
    temperature_domain: [Temperature; 2],
}

/// Immutable original source binding for one common solvent and 1..=64 solutes.
#[derive(Clone, Debug, PartialEq)]
pub struct SolubilityCaloricModelV1 {
    owner: Arc<CaloricOwner>,
}
impl SolubilityCaloricModelV1 {
    pub fn new(
        solvent_molar_volume_m3_mol: f64,
        solvent_cp: PureLiquidHeatCapacityV1,
        components: Vec<SolubilityCaloricComponentV1>,
    ) -> Result<Self, SolubilityCaloricErrorV1> {
        if !positive(solvent_molar_volume_m3_mol)
            || !(1..=MAX_SOLUBILITY_CALORIC_COMPONENTS_V1).contains(&components.len())
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "positive solvent molar volume and bounded component count required",
            ));
        }
        let mut domain = solvent_cp.temperature_domain();
        for (index, component) in components.iter().enumerate() {
            if component.law.solvent() != solvent_cp.species()
                || components[..index]
                    .iter()
                    .any(|other| other.law.solute() == component.law.solute())
            {
                return Err(SolubilityCaloricErrorV1::Invalid(
                    "common solvent and unique ordered solutes required",
                ));
            }
            for span in [
                component.law.temperature_domain(),
                component.common_phase_cp.temperature_domain(),
            ] {
                domain[0] = Temperature::new(domain[0].value().max(span[0].value()));
                domain[1] = Temperature::new(domain[1].value().min(span[1].value()));
            }
        }
        if domain[0].value() >= domain[1].value() {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "empty common temperature domain",
            ));
        }
        Ok(Self {
            owner: Arc::new(CaloricOwner {
                solvent_molar_volume_m3_mol,
                solvent_cp,
                components,
                temperature_domain: domain,
            }),
        })
    }

    pub fn solvent_molar_volume_m3_mol(&self) -> f64 {
        self.owner.solvent_molar_volume_m3_mol
    }
    pub fn solvent_heat_capacity(&self) -> &PureLiquidHeatCapacityV1 {
        &self.owner.solvent_cp
    }
    /// Original laws and caloric sources in the exact inventory-vector order.
    pub fn components(&self) -> &[SolubilityCaloricComponentV1] {
        &self.owner.components
    }
    pub fn temperature_domain(&self) -> [Temperature; 2] {
        self.owner.temperature_domain
    }

    /// Declare both phase stocks at a common initial temperature. This need not
    /// be equilibrium: add_heat with zero input closes adiabatic relaxation.
    pub fn prepare(
        &self,
        temperature: Temperature,
        solvent_moles: Amount,
        total_solute: &[Amount],
        dissolved_solute: &[Amount],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SolubilityCaloricStateV1, SolubilityCaloricErrorV1> {
        let (volume, totals) =
            self.inventory(temperature, solvent_moles, total_solute, &mut cancelled)?;
        self.prepared_state(
            temperature,
            solvent_moles,
            volume,
            totals,
            dissolved_solute,
            &mut cancelled,
        )
    }

    /// Declare a material reservoir by its original solvent volume. Moles are
    /// computed once as V/Vbar, while V is retained exactly for every saturation
    /// capacity. This avoids changing a conserved material volume through a
    /// floating-point quotient/product round trip. The same explicit constant
    /// molar volume supplies both the molar heat inventory and concentration basis.
    pub fn prepare_at_solvent_volume(
        &self,
        temperature: Temperature,
        solvent_volume: Volume,
        total_solute: &[Amount],
        dissolved_solute: &[Amount],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SolubilityCaloricStateV1, SolubilityCaloricErrorV1> {
        poll(&mut cancelled)?;
        if !positive(solvent_volume.value()) {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "positive original solvent volume required",
            ));
        }
        let solvent_moles =
            Amount::new(solvent_volume.value() / self.owner.solvent_molar_volume_m3_mol);
        let (_, totals) =
            self.inventory(temperature, solvent_moles, total_solute, &mut cancelled)?;
        self.prepared_state(
            temperature,
            solvent_moles,
            solvent_volume,
            totals,
            dissolved_solute,
            &mut cancelled,
        )
    }

    fn prepared_state(
        &self,
        temperature: Temperature,
        solvent_moles: Amount,
        solvent_volume: Volume,
        totals: Arc<[Amount]>,
        dissolved_solute: &[Amount],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SolubilityCaloricStateV1, SolubilityCaloricErrorV1> {
        if dissolved_solute.len() != totals.len() {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "ordered dissolved inventory length",
            ));
        }
        for (total, dissolved) in totals.iter().zip(dissolved_solute) {
            poll(cancelled)?;
            if !dissolved.value().is_finite()
                || dissolved.value() < 0.0
                || dissolved.value() > total.value()
            {
                return Err(SolubilityCaloricErrorV1::Invalid(
                    "dissolved stock must lie within its complete component inventory",
                ));
            }
        }
        self.make_state(
            temperature,
            solvent_moles,
            solvent_volume,
            totals,
            Arc::from(dissolved_solute),
            None,
            cancelled,
        )
    }

    /// Resolve finite stocks at an externally prescribed common temperature.
    /// This preparation alone does not claim an energy-balanced temperature change.
    pub fn equilibrate(
        &self,
        temperature: Temperature,
        solvent_moles: Amount,
        total_solute: &[Amount],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SolubilityCaloricStateV1, SolubilityCaloricErrorV1> {
        let (volume, totals) =
            self.inventory(temperature, solvent_moles, total_solute, &mut cancelled)?;
        self.equilibrium_state(temperature, solvent_moles, volume, totals, &mut cancelled)
    }

    /// Direct first-law increment for the SAME conserved inventory and sources.
    /// It uses local Cp integrals and changes in dissolved moles, never the
    /// subtraction of two large arbitrary reference enthalpies.
    pub fn enthalpy_change_j(
        &self,
        initial: &SolubilityCaloricStateV1,
        final_state: &SolubilityCaloricStateV1,
    ) -> Result<f64, SolubilityCaloricErrorV1> {
        if !self.owns(initial)
            || !self.owns(final_state)
            || initial.solvent_moles != final_state.solvent_moles
            || initial.solvent_volume != final_state.solvent_volume
            || initial.total_solute != final_state.total_solute
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "enthalpy increment requires the same caloric owner and inventory",
            ));
        }
        self.enthalpy_increment(
            initial,
            final_state.temperature,
            |index| Ok(final_state.dissolved_solute[index]),
            &mut || false,
        )
    }

    /// Heat, cool or adiabatically equilibrate the complete finite bath.
    ///
    /// At saturation, dn_diss/dT = n_diss * delta_H/(R*T^2). Thus the equilibrium
    /// enthalpy derivative is positive sensible capacity plus
    /// sum(n_diss * delta_H^2/(R*T^2)) over components with retained solid.
    /// The enthalpy is continuous at solid exhaustion. One bounded bisection
    /// therefore resolves a unique temperature without an inner nonlinear solve.
    pub fn add_heat(
        &self,
        initial: &SolubilityCaloricStateV1,
        heat_input_j: f64,
        controls: SolubilityCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SolubilityCaloricBalanceV1, SolubilityCaloricErrorV1> {
        poll(&mut cancelled)?;
        let [mut lo, mut hi] = controls.temperature_bracket.map(|t| t.value());
        if !self.owns(initial)
            || !heat_input_j.is_finite()
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
            return Err(SolubilityCaloricErrorV1::Invalid(
                "thermal source owner, heat input or solve controls",
            ));
        }
        if heat_input_j == 0.0
            && initial.equilibrium.is_some()
            && (lo..=hi).contains(&initial.temperature.value())
        {
            poll(&mut cancelled)?;
            return Ok(SolubilityCaloricBalanceV1 {
                final_state: initial.clone(),
                heat_input_j,
                enthalpy_change_j: 0.0,
                energy_residual_j: 0.0,
                temperature_bracket: [initial.temperature; 2],
                iterations: 0,
            });
        }
        // Scalar trials allocate nothing and integrate only local sensible
        // increments. Reference enthalpy and phase arrays are materialized once,
        // after accepting the root; source laws stay identical at every trial.
        let mut sample = |t| -> Result<_, SolubilityCaloricErrorV1> {
            poll(&mut cancelled)?;
            let temperature = Temperature::new(t);
            let change = self.enthalpy_increment(
                initial,
                temperature,
                |index| {
                    Ok(self.owner.components[index]
                        .law
                        .equilibrate(
                            temperature,
                            initial.solvent_volume,
                            initial.total_solute[index],
                        )?
                        .dissolved)
                },
                &mut cancelled,
            )?;
            let residual = checked(change - heat_input_j, "heat balance residual")?;
            Ok((temperature, change, residual))
        };
        let mut lower = sample(lo)?;
        let mut upper = sample(hi)?;
        if lower.2 > 0.0 || upper.2 < 0.0 {
            return Err(SolubilityCaloricErrorV1::HeatOutsideBracket {
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
                let result = SolubilityCaloricBalanceV1 {
                    final_state: self.equilibrium_state(
                        best.0,
                        initial.solvent_moles,
                        initial.solvent_volume,
                        Arc::clone(&initial.total_solute),
                        &mut cancelled,
                    )?,
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
        Err(SolubilityCaloricErrorV1::IterationLimit)
    }

    fn owns(&self, state: &SolubilityCaloricStateV1) -> bool {
        Arc::ptr_eq(&self.owner, &state.owner) || self.owner == state.owner
    }

    fn inventory(
        &self,
        temperature: Temperature,
        solvent_moles: Amount,
        total_solute: &[Amount],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Volume, Arc<[Amount]>), SolubilityCaloricErrorV1> {
        poll(cancelled)?;
        if !temperature.value().is_finite()
            || temperature.value() < self.owner.temperature_domain[0].value()
            || temperature.value() > self.owner.temperature_domain[1].value()
            || !positive(solvent_moles.value())
            || total_solute.len() != self.owner.components.len()
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "temperature, positive solvent inventory or ordered component count",
            ));
        }
        for amount in total_solute {
            poll(cancelled)?;
            if !amount.value().is_finite() || amount.value() < 0.0 {
                return Err(SolubilityCaloricErrorV1::Invalid(
                    "finite nonnegative total component inventory required",
                ));
            }
        }
        let volume = solvent_moles.value() * self.owner.solvent_molar_volume_m3_mol;
        if !positive(volume) {
            return Err(SolubilityCaloricErrorV1::Arithmetic(
                "positive solvent volume is unrepresentable",
            ));
        }
        Ok((Volume::new(volume), Arc::from(total_solute)))
    }

    fn equilibrium_state(
        &self,
        temperature: Temperature,
        solvent_moles: Amount,
        solvent_volume: Volume,
        totals: Arc<[Amount]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SolubilityCaloricStateV1, SolubilityCaloricErrorV1> {
        let mut partitions = Vec::with_capacity(totals.len());
        let mut dissolved = Vec::with_capacity(totals.len());
        for (component, total) in self.owner.components.iter().zip(totals.iter()) {
            poll(cancelled)?;
            let partition = component
                .law
                .equilibrate(temperature, solvent_volume, *total)?;
            dissolved.push(partition.dissolved);
            partitions.push(partition);
        }
        self.make_state(
            temperature,
            solvent_moles,
            solvent_volume,
            totals,
            dissolved.into(),
            Some(partitions.into()),
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn make_state(
        &self,
        temperature: Temperature,
        solvent_moles: Amount,
        solvent_volume: Volume,
        total_solute: Arc<[Amount]>,
        dissolved_solute: Arc<[Amount]>,
        equilibrium: Option<Arc<[SolubilityPartitionV1]>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SolubilityCaloricStateV1, SolubilityCaloricErrorV1> {
        let solvent = &self.owner.solvent_cp;
        let mut enthalpy = Sum::default();
        enthalpy.add(scaled(
            solvent_moles.value(),
            solvent.enthalpy_change_j_per_mol(solvent.reference_temperature(), temperature)?,
        )?)?;
        for ((component, total), dissolved) in self
            .owner
            .components
            .iter()
            .zip(total_solute.iter())
            .zip(dissolved_solute.iter())
        {
            poll(cancelled)?;
            if total.value() != 0.0 {
                let cp = &component.common_phase_cp;
                enthalpy.add(scaled(
                    total.value(),
                    cp.enthalpy_change_j_per_mol(cp.reference_temperature(), temperature)?,
                )?)?;
                enthalpy.add(scaled(
                    dissolved.value(),
                    component.law.dissolution_enthalpy_j_mol(),
                )?)?;
            }
        }
        let state = SolubilityCaloricStateV1 {
            owner: Arc::clone(&self.owner),
            temperature,
            solvent_moles,
            solvent_volume,
            total_solute,
            dissolved_solute,
            equilibrium,
            enthalpy_j: enthalpy.value()?,
        };
        poll(cancelled)?;
        Ok(state)
    }

    // One energy equation serves both public state increments and equilibrium
    // root trials. The callback supplies the actual final dissolved stock; every
    // sensible term is integrated once from the original temperature.
    fn enthalpy_increment(
        &self,
        initial: &SolubilityCaloricStateV1,
        temperature: Temperature,
        mut dissolved_at: impl FnMut(usize) -> Result<Amount, SolubilityCaloricErrorV1>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, SolubilityCaloricErrorV1> {
        poll(cancelled)?;
        let mut sum = Sum::default();
        sum.add(scaled(
            initial.solvent_moles.value(),
            self.owner
                .solvent_cp
                .enthalpy_change_j_per_mol(initial.temperature, temperature)?,
        )?)?;
        for (index, component) in self.owner.components.iter().enumerate() {
            poll(cancelled)?;
            let total = initial.total_solute[index].value();
            if total != 0.0 {
                sum.add(scaled(
                    total,
                    component
                        .common_phase_cp
                        .enthalpy_change_j_per_mol(initial.temperature, temperature)?,
                )?)?;
                sum.add(scaled(
                    dissolved_at(index)?.value() - initial.dissolved_solute[index].value(),
                    component.law.dissolution_enthalpy_j_mol(),
                )?)?;
            }
        }
        poll(cancelled)?;
        sum.value()
    }
}

/// Explicit work and acceptance bounds; unsupported heat never clamps to a limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolubilityCaloricSolveV1 {
    pub temperature_bracket: [Temperature; 2],
    pub temperature_tolerance_k: f64,
    pub energy_tolerance_j: f64,
    /// Maximum bisections, 1..=256. Each trial evaluates at most 64 original laws.
    pub maximum_iterations: usize,
}

/// Complete immutable phase stocks at one common temperature.
#[derive(Clone, Debug, PartialEq)]
pub struct SolubilityCaloricStateV1 {
    owner: Arc<CaloricOwner>,
    temperature: Temperature,
    solvent_moles: Amount,
    solvent_volume: Volume,
    total_solute: Arc<[Amount]>,
    dissolved_solute: Arc<[Amount]>,
    equilibrium: Option<Arc<[SolubilityPartitionV1]>>,
    enthalpy_j: f64,
}
impl SolubilityCaloricStateV1 {
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    pub const fn solvent_moles(&self) -> Amount {
        self.solvent_moles
    }
    pub const fn solvent_volume(&self) -> Volume {
        self.solvent_volume
    }
    pub fn total_solute(&self) -> &[Amount] {
        &self.total_solute
    }
    pub fn dissolved_solute(&self) -> &[Amount] {
        &self.dissolved_solute
    }
    pub fn solid_solute(&self) -> impl ExactSizeIterator<Item = Amount> + '_ {
        self.total_solute
            .iter()
            .zip(self.dissolved_solute.iter())
            .map(|(total, dissolved)| Amount::new(total.value() - dissolved.value()))
    }
    /// None identifies a declared initial phase stock without an equilibrium claim.
    pub fn equilibrium(&self) -> Option<&[SolubilityPartitionV1]> {
        self.equilibrium.as_deref()
    }
    /// Diagnostic reference enthalpy [J]. Balances use direct local increments.
    /// Pure solvent and solid sensible enthalpies are zero at their own Cp datums;
    /// each dissolved phase differs from its solid by its original law's delta H.
    pub const fn enthalpy_j(&self) -> f64 {
        self.enthalpy_j
    }
}

/// Completed first-law balance with unchanged solvent and complete solute stocks.
#[derive(Clone, Debug, PartialEq)]
pub struct SolubilityCaloricBalanceV1 {
    final_state: SolubilityCaloricStateV1,
    heat_input_j: f64,
    enthalpy_change_j: f64,
    energy_residual_j: f64,
    temperature_bracket: [Temperature; 2],
    iterations: usize,
}
impl SolubilityCaloricBalanceV1 {
    pub const fn final_state(&self) -> &SolubilityCaloricStateV1 {
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
    /// Numerical temperature bracket, not a rigorous interval certificate.
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
fn checked(value: f64, what: &'static str) -> Result<f64, SolubilityCaloricErrorV1> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SolubilityCaloricErrorV1::Arithmetic(what))
    }
}
fn scaled(moles: f64, molar: f64) -> Result<f64, SolubilityCaloricErrorV1> {
    let result = checked(moles * molar, "extensive enthalpy contribution")?;
    if result == 0.0 && moles != 0.0 && molar != 0.0 {
        return Err(SolubilityCaloricErrorV1::Arithmetic(
            "nonzero enthalpy contribution underflows",
        ));
    }
    Ok(result)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SolubilityCaloricErrorV1> {
    if cancelled() {
        Err(SolubilityCaloricErrorV1::Cancelled)
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
    fn add(&mut self, value: f64) -> Result<(), SolubilityCaloricErrorV1> {
        let next = checked(self.main + value, "enthalpy sum")?;
        self.correction = checked(
            self.correction
                + if self.main.abs() >= value.abs() {
                    (self.main - next) + value
                } else {
                    (value - next) + self.main
                },
            "enthalpy sum correction",
        )?;
        self.main = next;
        Ok(())
    }
    fn value(&self) -> Result<f64, SolubilityCaloricErrorV1> {
        checked(self.main + self.correction, "enthalpy sum correction")
    }
}

#[cfg(test)]
mod tests;
