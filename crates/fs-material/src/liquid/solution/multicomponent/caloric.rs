//! Source-qualified multicomponent liquid heating, mixing and phase response.
//!
//! Every species retains its original pure-liquid Cp source. The shared L1
//! owner closes energy and entropy; only a completed thermal endpoint refreshes
//! actual mixture material properties at its temperature and composition.
//! The same original source binding supplies the differential isobaric heat
//! capacity of an already accepted local equilibrium phase candidate.

mod response;
pub use response::ResolvedMulticomponentLiquidCaloricResponse;
/// Finite heat exchange along an accepted local equilibrium phase branch.
pub mod equilibrium;
pub use equilibrium::{
    MulticomponentLiquidEquilibriumCaloricBalance,
    ResolvedMulticomponentLiquidEquilibriumCaloricState,
};

use super::{
    MulticomponentLiquidError, MulticomponentLiquidSolutionLaw, ResolvedMulticomponentLiquidCell,
    SourceOwner as LiquidSourceOwner,
};
use crate::thermal::transport::isobaric::IsobaricHeatCapacityLaw;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Pressure, QuantitySpec, Temperature};
use fs_thermochem::solution::caloric::PureLiquidHeatCapacityV1;
use fs_thermochem::solution::multicomponent::caloric::equilibrium::MulticomponentEquilibriumCaloricErrorV1;
use fs_thermochem::solution::multicomponent::caloric::{
    MulticomponentCaloricBalanceV1, MulticomponentCaloricErrorV1, MulticomponentCaloricSolveV1,
    MulticomponentCaloricStateV1, MulticomponentLiquidCaloricModelV1,
};
use fs_thermochem::solution::multicomponent::equilibrium::caloric::MulticomponentCaloricResponseErrorV1;
use fs_thermochem::{MolarThermalQuantityV1, SpeciesId};
use std::sync::Arc;

/// Invalid source binding, refused thermal balance, or unsupported final material.
#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentLiquidCaloricError {
    Invalid(&'static str),
    Caloric(MulticomponentCaloricErrorV1),
    CoreResponse(MulticomponentCaloricResponseErrorV1),
    EquilibriumCaloric(MulticomponentEquilibriumCaloricErrorV1),
    Material(MulticomponentLiquidError),
    Cancelled,
}
impl core::fmt::Display for MulticomponentLiquidCaloricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound multicomponent liquid calorics: {self:?}")
    }
}
impl std::error::Error for MulticomponentLiquidCaloricError {}
impl From<MulticomponentCaloricErrorV1> for MulticomponentLiquidCaloricError {
    fn from(error: MulticomponentCaloricErrorV1) -> Self {
        Self::Caloric(error)
    }
}
impl From<MulticomponentLiquidError> for MulticomponentLiquidCaloricError {
    fn from(error: MulticomponentLiquidError) -> Self {
        Self::Material(error)
    }
}
impl From<MulticomponentCaloricResponseErrorV1> for MulticomponentLiquidCaloricError {
    fn from(error: MulticomponentCaloricResponseErrorV1) -> Self {
        Self::CoreResponse(error)
    }
}
impl From<MulticomponentEquilibriumCaloricErrorV1> for MulticomponentLiquidCaloricError {
    fn from(error: MulticomponentEquilibriumCaloricErrorV1) -> Self {
        Self::EquilibriumCaloric(error)
    }
}
type Result<T> = core::result::Result<T, MulticomponentLiquidCaloricError>;

#[derive(Clone, Debug, PartialEq)]
struct SourceOwner {
    liquid: Arc<LiquidSourceOwner>,
    sources: Box<[IsobaricHeatCapacityLaw]>,
    model: MulticomponentLiquidCaloricModelV1,
    pressure_axis: String,
}

/// Original mixture transport and one explicitly qualified pure-liquid Cp law
/// for every original species. The caller qualifies each source as that pure
/// liquid at the declared pressure; a label or dimensions cannot establish it.
/// Transport Cp remains a separate property, with no implied consistency claim.
pub struct MulticomponentLiquidCaloricLaw<'a, 'card> {
    liquid: &'a MulticomponentLiquidSolutionLaw<'card>,
    owner: Arc<SourceOwner>,
}
impl<'a, 'card> MulticomponentLiquidCaloricLaw<'a, 'card> {
    /// Canonicalize complete source receipts by species, retaining every Cp
    /// knot. Convert mass Cp to molar Cp using the original component masses.
    /// Sources require a common datum and the mixture's typed static pressure.
    pub fn try_new(
        liquid: &'a MulticomponentLiquidSolutionLaw<'card>,
        mut pure_component_sources: Vec<(SpeciesId, IsobaricHeatCapacityLaw)>,
        pressure_axis: &str,
    ) -> Result<Self> {
        pure_component_sources.sort_by(|a, b| a.0.cmp(&b.0));
        let quantity = liquid.point.axis_quantities().get(pressure_axis);
        let static_pressure =
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
        let pressure = liquid.point.axes().get(pressure_axis).copied();
        if !pure_component_sources
            .iter()
            .map(|(id, _)| id)
            .eq(liquid.species())
            || pressure_axis == liquid.owner.temperature_axis
            || liquid
                .components()
                .iter()
                .any(|component| component.mole_fraction_axis.as_deref() == Some(pressure_axis))
            || !quantity.is_some_and(|q| {
                *q == static_pressure || *q == QuantitySpec::dimensional(Pressure::DIMS)
            })
            || !pressure.is_some_and(|p| p.is_finite() && p > 0.0)
            || pure_component_sources
                .iter()
                .any(|(_, source)| Some(source.pressure_pa()) != pressure)
        {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "one pure-liquid source per original species at the original static pressure is required",
            ));
        }
        let mut capacities = Vec::with_capacity(pure_component_sources.len());
        let mut sources = Vec::with_capacity(pure_component_sources.len());
        for ((species, source), component) in
            pure_component_sources.into_iter().zip(liquid.components())
        {
            capacities.push(
                PureLiquidHeatCapacityV1::new(
                    species,
                    Temperature::new(source.datum_temperature_k()),
                    source
                        .source_line()
                        .knots
                        .iter()
                        .map(|(temperature, answer)| {
                            (
                                Temperature::new(*temperature),
                                MolarThermalQuantityV1::new(
                                    answer.evidence.value.value * component.molar_mass.value(),
                                ),
                            )
                        })
                        .collect(),
                )
                .map_err(MulticomponentCaloricErrorV1::from)?,
            );
            sources.push(source);
        }
        let model =
            MulticomponentLiquidCaloricModelV1::new(liquid.thermodynamics().clone(), capacities)?;
        Ok(Self {
            liquid,
            owner: Arc::new(SourceOwner {
                liquid: Arc::clone(&liquid.owner),
                sources: sources.into_boxed_slice(),
                model,
                pressure_axis: pressure_axis.to_owned(),
            }),
        })
    }

    /// Exact original Cp laws and receipts in the original canonical species order.
    pub fn sources(&self) -> &[IsobaricHeatCapacityLaw] {
        &self.owner.sources
    }
    pub fn model(&self) -> &MulticomponentLiquidCaloricModelV1 {
        &self.owner.model
    }
    pub const fn liquid_law(&self) -> &MulticomponentLiquidSolutionLaw<'card> {
        self.liquid
    }
    /// Shared numeric domain; actual mixture sources must also admit the endpoint.
    pub fn temperature_domain(&self) -> [Temperature; 2] {
        self.model().temperature_domain()
    }

    /// Prepare an explicitly homogeneous finite parcel at the actual source
    /// temperature and composition. No equilibrium phase split is inferred.
    pub fn homogeneous(
        &self,
        temperature: Temperature,
        component_moles: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedMulticomponentLiquidCaloricState> {
        poll(&mut cancelled)?;
        let caloric = self
            .model()
            .homogeneous(temperature, component_moles, &mut cancelled)?;
        self.resolve(caloric, &mut cancelled)
    }

    /// Mix unequal-temperature parcels with signed external heat. The shared
    /// numeric owner includes ideal/nonideal entropy and excess mixing heat;
    /// its closed zero-heat entropy guard rejects an inadmissible homogeneous mix.
    pub fn mix(
        &self,
        left: &ResolvedMulticomponentLiquidCaloricState,
        right: &ResolvedMulticomponentLiquidCaloricState,
        heat_input_j: f64,
        controls: MulticomponentCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidCaloricBalance> {
        poll(&mut cancelled)?;
        self.admit(left)?;
        self.admit(right)?;
        let balance = self.model().mix(
            &left.caloric,
            &right.caloric,
            heat_input_j,
            controls,
            &mut cancelled,
        )?;
        self.resolve_balance(balance, &mut cancelled)
    }

    /// Add or remove heat without altering any component inventory. Source
    /// refusal or cancellation publishes no partially resolved thermal endpoint.
    pub fn add_heat(
        &self,
        initial: &ResolvedMulticomponentLiquidCaloricState,
        heat_input_j: f64,
        controls: MulticomponentCaloricSolveV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidCaloricBalance> {
        poll(&mut cancelled)?;
        self.admit(initial)?;
        let balance =
            self.model()
                .add_heat(&initial.caloric, heat_input_j, controls, &mut cancelled)?;
        self.resolve_balance(balance, &mut cancelled)
    }

    /// Source-relative enthalpy change for the same original component inventory.
    pub fn enthalpy_change_j(
        &self,
        initial: &ResolvedMulticomponentLiquidCaloricState,
        final_state: &ResolvedMulticomponentLiquidCaloricState,
    ) -> Result<f64> {
        self.admit(initial)?;
        self.admit(final_state)?;
        Ok(self
            .model()
            .enthalpy_change_j(&initial.caloric, &final_state.caloric)?)
    }

    /// Source-relative entropy change for the same original component inventory.
    pub fn entropy_change_j_per_k(
        &self,
        initial: &ResolvedMulticomponentLiquidCaloricState,
        final_state: &ResolvedMulticomponentLiquidCaloricState,
    ) -> Result<f64> {
        self.admit(initial)?;
        self.admit(final_state)?;
        Ok(self
            .model()
            .entropy_change_j_per_k(&initial.caloric, &final_state.caloric)?)
    }

    fn admit(&self, state: &ResolvedMulticomponentLiquidCaloricState) -> Result<()> {
        if !Arc::ptr_eq(&self.owner, &state.owner) && self.owner != state.owner {
            return Err(MulticomponentLiquidCaloricError::Invalid(
                "parcel requires its original caloric and mixture source binding",
            ));
        }
        Ok(())
    }
    fn resolve(
        &self,
        caloric: MulticomponentCaloricStateV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ResolvedMulticomponentLiquidCaloricState> {
        let liquid = self.liquid.homogeneous(
            caloric.temperature(),
            caloric.component_moles(),
            &mut *cancelled,
        )?;
        poll(cancelled)?;
        Ok(ResolvedMulticomponentLiquidCaloricState {
            owner: Arc::clone(&self.owner),
            caloric,
            liquid,
        })
    }
    fn resolve_balance(
        &self,
        balance: MulticomponentCaloricBalanceV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MulticomponentLiquidCaloricBalance> {
        let final_state = self.resolve(balance.final_state().clone(), cancelled)?;
        Ok(MulticomponentLiquidCaloricBalance {
            final_state,
            balance,
        })
    }
}

/// One complete caloric parcel with actual composition-resolved material sources.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedMulticomponentLiquidCaloricState {
    owner: Arc<SourceOwner>,
    caloric: MulticomponentCaloricStateV1,
    liquid: ResolvedMulticomponentLiquidCell,
}
impl ResolvedMulticomponentLiquidCaloricState {
    pub const fn caloric_state(&self) -> &MulticomponentCaloricStateV1 {
        &self.caloric
    }
    /// Existing resolved cell, including actual transport, mass and volume.
    pub const fn liquid(&self) -> &ResolvedMulticomponentLiquidCell {
        &self.liquid
    }
    pub fn sources(&self) -> &[IsobaricHeatCapacityLaw] {
        &self.owner.sources
    }
}

/// Completed heat/entropy balance and its actual final material resolution.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentLiquidCaloricBalance {
    final_state: ResolvedMulticomponentLiquidCaloricState,
    balance: MulticomponentCaloricBalanceV1,
}
impl MulticomponentLiquidCaloricBalance {
    pub const fn final_state(&self) -> &ResolvedMulticomponentLiquidCaloricState {
        &self.final_state
    }
    pub const fn liquid(&self) -> &ResolvedMulticomponentLiquidCell {
        self.final_state.liquid()
    }
    pub const fn balance(&self) -> &MulticomponentCaloricBalanceV1 {
        &self.balance
    }
    /// System entropy change; external heat needs a separate reservoir account.
    pub const fn entropy_change_j_per_k(&self) -> f64 {
        self.balance.entropy_change_j_per_k()
    }
    pub fn into_final_state(self) -> ResolvedMulticomponentLiquidCaloricState {
        self.final_state
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentLiquidCaloricError::Cancelled)
    } else {
        Ok(())
    }
}
