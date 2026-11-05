//! Shared material responses at the actual kinetic temperature and crystallinity.
//!
//! Crystallinity is an independent internal variable, not an equilibrium liquid
//! fraction. The parent can itself be an amorphous solid. This bridge updates
//! explicitly named crystalline mass/volume coordinates and pins every original
//! property claim. Thermal, electrical and mechanical consumers can therefore
//! use one source snapshot without deriving a second phase state from enthalpy.
//!
//! No conductivity, modulus or other mixture rule is inferred. Each supplied
//! property must apply to this same process, composition, pressure and kinetic
//! history. A local source answer does not establish support between endpoints,
//! morphology, a spatial phase front, or a calibrated crystallization model.

use super::{CrystallizationCaloricPair, CrystallizationVolumeState};
use crate::phase::PhaseStateError;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarPropertyRequirement, is_fixed_external_coordinate, resolve_material_state_point,
};
use crate::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_thermal::crystallization::{
    CrystallizationError, CrystallizationLaw, CrystallizationSolve, CrystallizationState,
};

/// Explicit coordinates of one original material's kinetic property data.
#[derive(Clone, Copy, Debug)]
pub struct CrystallizationMaterialPoint<'a> {
    pub material: &'a MaterialCard,
    pub context: &'a QueryPoint,
    pub temperature_axis: &'a str,
    /// Absolute crystalline MASS fraction, not conversion divided by Xmax.
    pub crystalline_mass_fraction_axis: Option<&'a str>,
    /// Crystalline VOLUME fraction from the original constituent densities.
    pub crystalline_volume_fraction_axis: Option<&'a str>,
}

/// Original constitutive/source failures remain distinguishable at the bridge.
#[derive(Clone, Debug, PartialEq)]
pub enum CrystallizationPropertyError {
    Invalid(&'static str),
    Phase(PhaseStateError),
    Kinetic(CrystallizationError),
    Material(ThermalMaterialError),
}
impl core::fmt::Display for CrystallizationPropertyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "crystallization property: {self:?}")
    }
}
impl std::error::Error for CrystallizationPropertyError {}
impl From<PhaseStateError> for CrystallizationPropertyError {
    fn from(error: PhaseStateError) -> Self {
        Self::Phase(error)
    }
}
impl From<CrystallizationError> for CrystallizationPropertyError {
    fn from(error: CrystallizationError) -> Self {
        Self::Kinetic(error)
    }
}
impl From<ThermalMaterialError> for CrystallizationPropertyError {
    fn from(error: ThermalMaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<MaterialStatePointError> for CrystallizationPropertyError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(ThermalMaterialError::Material(error))
    }
}

/// One synchronized phase, original-density volume and material-source answer.
#[derive(Clone, Debug, PartialEq)]
pub struct CrystallizationPropertyState {
    phase: CrystallizationState,
    volume: CrystallizationVolumeState,
    material: ResolvedMaterialStatePoint,
}
impl CrystallizationPropertyState {
    pub const fn phase(&self) -> CrystallizationState {
        self.phase
    }
    pub const fn volume(&self) -> CrystallizationVolumeState {
        self.volume
    }
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
}

/// Immutable preparation shared by every consumer of these material properties.
///
/// Preparation selects the original claims once. Evaluation refreshes that
/// exact bundle; it never selects another claim after a domain failure. Source
/// queries allocate their ordinary bounded receipt bundle. No physical history
/// is advanced or cached speculatively in this object.
#[derive(Clone)]
pub struct CrystallizationPropertySource<'a> {
    material: &'a MaterialCard,
    pair: CrystallizationCaloricPair<'a>,
    law: CrystallizationLaw<'a>,
    original: ResolvedMaterialStatePoint,
    temperature_axis: String,
    mass_fraction_axis: Option<String>,
    volume_fraction_axis: Option<String>,
    solve: CrystallizationSolve,
}
impl<'a> CrystallizationPropertySource<'a> {
    /// Bind actual kinetic history to arbitrarily selected static scalar
    /// properties, with their original unit/range requirements and receipts.
    ///
    /// The caloric pair, kinetic law and material card must be the SAME original
    /// sources. Temperature and optional fraction axes must already exist with
    /// exact quantity descriptors. Other physical coordinates stay fixed;
    /// strain/rate/field/history-dependent sources need their separate owner.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        pair: CrystallizationCaloricPair<'a>,
        law: CrystallizationLaw<'a>,
        initial: CrystallizationState,
        point: CrystallizationMaterialPoint<'a>,
        requirements: &[ScalarPropertyRequirement],
        selection: MaterialPropertySelection,
        solve: CrystallizationSolve,
    ) -> Result<Self, CrystallizationPropertyError> {
        if point.material.content_hash() != pair.parent().curve().material_card_identity()
            || law.parent_caloric_identity() != Some(*pair.parent().identity().as_bytes())
            || law.crystalline_caloric_identity() != Some(*pair.crystal().identity().as_bytes())
            || law.config().temperature_range_k != pair.temperature_range_k()
        {
            return Err(invalid(
                "properties require the original kinetic law, caloric pair and material card",
            ));
        }
        admit_axes(point)?;
        law.validate_state(initial, solve)?;
        let volume =
            pair.volume_state(initial.temperature_k(), initial.crystalline_mass_fraction())?;
        let query = query_at(point, initial, volume)?;
        let original =
            resolve_material_state_point(point.material, &query, requirements, selection)?;
        let thermal = ThermalMaterialPoint {
            material: point.material,
            context: &query,
            temperature_axis: point.temperature_axis,
            liquid_fraction_axis: None,
        };
        let axes: Vec<_> = [
            point.crystalline_mass_fraction_axis,
            point.crystalline_volume_fraction_axis,
        ]
        .into_iter()
        .flatten()
        .map(|axis| (axis, QuantitySpec::dimensional(Dims::NONE)))
        .collect();
        for property in original.properties() {
            thermal.admit_static_with_axes(property, &axes)?;
            let claim = point
                .material
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or(invalid("original crystallization property claim is absent"))?;
            if claim.validity.bounds().keys().any(|axis| {
                axis != point.temperature_axis
                    && Some(axis.as_str()) != point.crystalline_mass_fraction_axis
                    && Some(axis.as_str()) != point.crystalline_volume_fraction_axis
                    && !is_fixed_external_coordinate(
                        axis,
                        query.axis_quantities().get(axis).copied(),
                        query.axes().get(axis).copied(),
                    )
            }) {
                return Err(invalid(
                    "property applicability contains an unowned physical coordinate",
                ));
            }
        }
        Ok(Self {
            material: point.material,
            pair,
            law,
            original,
            temperature_axis: point.temperature_axis.to_owned(),
            mass_fraction_axis: point.crystalline_mass_fraction_axis.map(str::to_owned),
            volume_fraction_axis: point.crystalline_volume_fraction_axis.map(str::to_owned),
            solve,
        })
    }

    pub const fn original_source(&self) -> &ResolvedMaterialStatePoint {
        &self.original
    }

    /// Refresh the pinned claims at an accepted or speculative original-law
    /// state. Equal T with different kinetic age can produce different material
    /// properties; equal enthalpy alone never identifies the physical state.
    /// Density is still owned by the caloric pair. A queried density or Cp is
    /// additional source evidence, not permission to replace stored mass,
    /// calorics, pressure work or an existing geometric constraint.
    pub fn evaluate(
        &self,
        state: CrystallizationState,
    ) -> Result<CrystallizationPropertyState, CrystallizationPropertyError> {
        self.law.validate_state(state, self.solve)?;
        let volume = self
            .pair
            .volume_state(state.temperature_k(), state.crystalline_mass_fraction())?;
        let context = self.original.query();
        let query = query_at(
            CrystallizationMaterialPoint {
                material: self.material,
                context: &context,
                temperature_axis: &self.temperature_axis,
                crystalline_mass_fraction_axis: self.mass_fraction_axis.as_deref(),
                crystalline_volume_fraction_axis: self.volume_fraction_axis.as_deref(),
            },
            state,
            volume,
        )?;
        let material = self.original.at_query_point(self.material, &query)?;
        Ok(CrystallizationPropertyState {
            phase: state,
            volume,
            material,
        })
    }
}

fn admit_axes(point: CrystallizationMaterialPoint<'_>) -> Result<(), CrystallizationPropertyError> {
    let temperature = point.context.axis_quantities().get(point.temperature_axis);
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    if !point.context.axes().contains_key(point.temperature_axis)
        || temperature
            .is_none_or(|q| *q != QuantitySpec::dimensional(Temperature::DIMS) && *q != absolute)
    {
        return Err(invalid(
            "an existing typed absolute-temperature axis is required",
        ));
    }
    let fractions = [
        point.crystalline_mass_fraction_axis,
        point.crystalline_volume_fraction_axis,
    ];
    if fractions[0].is_some() && fractions[0] == fractions[1] {
        return Err(invalid(
            "crystalline mass and volume fractions require distinct axes",
        ));
    }
    for axis in fractions.into_iter().flatten() {
        if axis == point.temperature_axis
            || !point.context.axes().contains_key(axis)
            || point.context.axis_quantities().get(axis)
                != Some(&QuantitySpec::dimensional(Dims::NONE))
        {
            return Err(invalid(
                "crystalline fractions require distinct existing dimensional scalar axes",
            ));
        }
    }
    Ok(())
}

fn query_at(
    point: CrystallizationMaterialPoint<'_>,
    state: CrystallizationState,
    volume: CrystallizationVolumeState,
) -> Result<QueryPoint, CrystallizationPropertyError> {
    let coordinate = |e: fs_matdb::MatDbError| {
        CrystallizationPropertyError::Material(ThermalMaterialError::Coordinate(e.to_string()))
    };
    let mut query = point
        .context
        .clone()
        .with(point.temperature_axis, state.temperature_k())
        .map_err(coordinate)?;
    for (axis, value) in [
        (
            point.crystalline_mass_fraction_axis,
            state.crystalline_mass_fraction(),
        ),
        (
            point.crystalline_volume_fraction_axis,
            volume.crystalline_volume_fraction,
        ),
    ] {
        if let Some(axis) = axis {
            query = query.with(axis, value).map_err(coordinate)?;
        }
    }
    Ok(query)
}

fn invalid(what: &'static str) -> CrystallizationPropertyError {
    CrystallizationPropertyError::Invalid(what)
}
