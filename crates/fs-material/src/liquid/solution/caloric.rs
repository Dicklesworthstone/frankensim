//! Source-bound closed isobaric liquid mixing and heating.
//!
//! Ordered pure-component Cp sources supply sensible enthalpy; the original
//! pair law supplies excess mixing enthalpy exactly once. Numeric temperature
//! trials use the shared L1 caloric owner. Only a completed thermal endpoint
//! refreshes the original mixture transport sources at each actual phase.

/// Fixed-composition enthalpy storage on the existing conservative heat timeline.
pub mod storage;
/// Explicit homogeneous feeds and their energy-conserving relaxation to equilibrium.
pub mod homogeneous;

use super::{
    BinaryLiquidError, BinaryLiquidOwner, BinaryLiquidSolutionLaw, ResolvedBinaryLiquidEquilibrium,
};
use crate::thermal::transport::isobaric::IsobaricHeatCapacityLaw;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Pressure, QuantitySpec, Temperature};
use fs_thermochem::solution::BinaryCoexistenceOptionsV1;
use fs_thermochem::solution::caloric::{
    BinaryCaloricBalanceV1, BinaryCaloricErrorV1, BinaryCaloricSolveV1, BinaryCaloricStateV1,
    BinaryLiquidCaloricModelV1, PureLiquidHeatCapacityV1,
};
use fs_thermochem::{MolarThermalQuantityV1, SpeciesId};
use std::sync::Arc;

/// Invalid source binding, thermal refusal, or unsupported final material state.
#[derive(Clone, Debug, PartialEq)]
pub enum BinaryLiquidCaloricError {
    Invalid(&'static str),
    Caloric(BinaryCaloricErrorV1),
    Material(BinaryLiquidError),
    Cancelled,
}
impl core::fmt::Display for BinaryLiquidCaloricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound binary liquid calorics: {self:?}")
    }
}
impl std::error::Error for BinaryLiquidCaloricError {}
impl From<BinaryCaloricErrorV1> for BinaryLiquidCaloricError {
    fn from(error: BinaryCaloricErrorV1) -> Self {
        Self::Caloric(error)
    }
}
impl From<BinaryLiquidError> for BinaryLiquidCaloricError {
    fn from(error: BinaryLiquidError) -> Self {
        Self::Material(error)
    }
}

#[derive(Clone, Debug, PartialEq)]
struct SourceOwner {
    liquid: Arc<BinaryLiquidOwner>,
    sources: [IsobaricHeatCapacityLaw; 2],
    model: BinaryLiquidCaloricModelV1,
    pressure_axis: String,
}

/// One original mixture law and two explicitly qualified pure-liquid Cp laws.
///
/// The caller asserts that each ordered source describes the named pure liquid
/// at the mixture pressure. Neither chemistry labels nor dimensional agreement
/// establish that qualification. Total-mixture, apparent or gas Cp cannot be
/// substituted. Transport Cp remains a separate source property; this adapter
/// does not certify its consistency with the supplied caloric pair.
pub struct BinaryLiquidCaloricLaw<'a, 'card> {
    liquid: &'a BinaryLiquidSolutionLaw<'card>,
    owner: Arc<SourceOwner>,
}

impl<'a, 'card> BinaryLiquidCaloricLaw<'a, 'card> {
    /// Preserve both original source laws and their complete receipts. Sources
    /// must use the pair's exact species order, one common caloric datum, and
    /// the same positive static pressure as the mixture's explicitly typed axis.
    /// Cp is converted from J/(kg K) to J/(mol K) using the original molar mass;
    /// every source breakpoint survives, without a resampled enthalpy chart.
    pub fn try_new(
        liquid: &'a BinaryLiquidSolutionLaw<'card>,
        pure_component_sources: [(SpeciesId, IsobaricHeatCapacityLaw); 2],
        pressure_axis: &str,
    ) -> Result<Self, BinaryLiquidCaloricError> {
        let [(first_species, first), (second_species, second)] = pure_component_sources;
        let species = [first_species, second_species];
        let sources = [first, second];
        let quantity = liquid.point.axis_quantities().get(pressure_axis);
        let static_pressure =
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
        let pressure = liquid.point.axes().get(pressure_axis).copied();
        if &species != liquid.thermodynamics().species()
            || pressure_axis == liquid.temperature_axis
            || pressure_axis == liquid.composition_axis
            || !quantity.is_some_and(|q| {
                *q == static_pressure || *q == QuantitySpec::dimensional(Pressure::DIMS)
            })
            || !pressure.is_some_and(|p| p.is_finite() && p > 0.0)
            || sources
                .iter()
                .any(|source| Some(source.pressure_pa()) != pressure)
        {
            return Err(BinaryLiquidCaloricError::Invalid(
                "ordered pure-liquid species and original static mixture pressure must match",
            ));
        }
        let pure = |index: usize| {
            PureLiquidHeatCapacityV1::new(
                species[index].clone(),
                Temperature::new(sources[index].datum_temperature_k()),
                sources[index]
                    .source_line()
                    .knots
                    .iter()
                    .map(|(t, answer)| {
                        (
                            Temperature::new(*t),
                            MolarThermalQuantityV1::new(
                                answer.evidence.value.value * liquid.molar_masses()[index].value(),
                            ),
                        )
                    })
                    .collect(),
            )
        };
        let model =
            BinaryLiquidCaloricModelV1::new(liquid.thermodynamics().clone(), [pure(0)?, pure(1)?])?;
        Ok(Self {
            liquid,
            owner: Arc::new(SourceOwner {
                liquid: Arc::clone(&liquid.owner),
                sources,
                model,
                pressure_axis: pressure_axis.to_owned(),
            }),
        })
    }

    /// Exact original Cp laws, including their source answers and identities.
    pub fn sources(&self) -> &[IsobaricHeatCapacityLaw; 2] {
        &self.owner.sources
    }
    /// Common numeric temperature domain; final mixture claims must also admit T.
    pub fn temperature_domain(&self) -> [Temperature; 2] {
        self.owner.model.temperature_domain()
    }

    /// Prepare a finite equilibrium parcel at a prescribed temperature.
    pub fn equilibrate(
        &self,
        temperature: Temperature,
        component_moles: [f64; 2],
        options: BinaryCoexistenceOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedBinaryLiquidCaloricState, BinaryLiquidCaloricError> {
        poll(&mut cancelled)?;
        let caloric =
            self.owner
                .model
                .equilibrate(temperature, component_moles, options, &mut cancelled)?;
        self.resolve(caloric, &mut cancelled)
    }

    /// Mix two original finite parcels, allowing unequal temperatures and phase
    /// counts. Heat is positive into the liquids. Every numeric trial retains
    /// the original caloric law; only the accepted endpoint queries transport.
    pub fn mix(
        &self,
        left: &ResolvedBinaryLiquidCaloricState,
        right: &ResolvedBinaryLiquidCaloricState,
        heat_input_j: f64,
        controls: BinaryCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryLiquidCaloricBalance, BinaryLiquidCaloricError> {
        poll(&mut cancelled)?;
        self.admit(left)?;
        self.admit(right)?;
        let balance = self.owner.model.mix(
            &left.caloric,
            &right.caloric,
            heat_input_j,
            controls,
            &mut cancelled,
        )?;
        self.resolve_balance(balance, &mut cancelled)
    }

    /// Heat or cool one original parcel. Inputs are immutable on source failure,
    /// unsupported heat, incomplete convergence, or cancellation at any boundary.
    pub fn add_heat(
        &self,
        initial: &ResolvedBinaryLiquidCaloricState,
        heat_input_j: f64,
        controls: BinaryCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryLiquidCaloricBalance, BinaryLiquidCaloricError> {
        poll(&mut cancelled)?;
        self.admit(initial)?;
        let balance =
            self.owner
                .model
                .add_heat(&initial.caloric, heat_input_j, controls, &mut cancelled)?;
        self.resolve_balance(balance, &mut cancelled)
    }

    fn admit(
        &self,
        state: &ResolvedBinaryLiquidCaloricState,
    ) -> Result<(), BinaryLiquidCaloricError> {
        if !Arc::ptr_eq(&self.owner, &state.owner) && self.owner != state.owner {
            return Err(BinaryLiquidCaloricError::Invalid(
                "parcel requires its original caloric and mixture source binding",
            ));
        }
        Ok(())
    }

    fn resolve(
        &self,
        caloric: BinaryCaloricStateV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ResolvedBinaryLiquidCaloricState, BinaryLiquidCaloricError> {
        let phases = caloric
            .equilibrium()
            .copied()
            .ok_or(BinaryLiquidCaloricError::Invalid(
                "material refresh requires an equilibrated caloric state",
            ))?;
        let equilibrium =
            self.liquid
                .resolve_equilibrium(phases, caloric.component_moles(), &mut *cancelled)?;
        poll(cancelled)?;
        Ok(ResolvedBinaryLiquidCaloricState {
            owner: Arc::clone(&self.owner),
            caloric,
            equilibrium,
        })
    }

    fn resolve_balance(
        &self,
        balance: BinaryCaloricBalanceV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BinaryLiquidCaloricBalance, BinaryLiquidCaloricError> {
        let final_state = self.resolve(balance.final_state().clone(), cancelled)?;
        Ok(BinaryLiquidCaloricBalance {
            final_state,
            balance,
        })
    }
}

/// A complete caloric parcel and its actual composition-resolved transport phases.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBinaryLiquidCaloricState {
    owner: Arc<SourceOwner>,
    caloric: BinaryCaloricStateV1,
    equilibrium: ResolvedBinaryLiquidEquilibrium,
}
impl ResolvedBinaryLiquidCaloricState {
    pub const fn caloric_state(&self) -> &BinaryCaloricStateV1 {
        &self.caloric
    }
    /// The unchanged shared equilibrium type accepted by existing tank consumers.
    pub const fn equilibrium(&self) -> &ResolvedBinaryLiquidEquilibrium {
        &self.equilibrium
    }
    pub fn sources(&self) -> &[IsobaricHeatCapacityLaw; 2] {
        &self.owner.sources
    }
}

/// Closed heat account together with the final source-resolved inventory.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryLiquidCaloricBalance {
    final_state: ResolvedBinaryLiquidCaloricState,
    balance: BinaryCaloricBalanceV1,
}
impl BinaryLiquidCaloricBalance {
    pub const fn final_state(&self) -> &ResolvedBinaryLiquidCaloricState {
        &self.final_state
    }
    pub const fn equilibrium(&self) -> &ResolvedBinaryLiquidEquilibrium {
        self.final_state.equilibrium()
    }
    pub const fn balance(&self) -> &BinaryCaloricBalanceV1 {
        &self.balance
    }
    pub fn into_final_state(self) -> ResolvedBinaryLiquidCaloricState {
        self.final_state
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), BinaryLiquidCaloricError> {
    if cancelled() {
        Err(BinaryLiquidCaloricError::Cancelled)
    } else {
        Ok(())
    }
}
