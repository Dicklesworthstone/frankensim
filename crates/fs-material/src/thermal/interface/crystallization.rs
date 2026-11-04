//! Source-bound contact between two uniform, finite kinetic phase inventories.
//!
//! The original ordered interface supplies area-specific resistance at the two
//! ACTUAL endpoint temperatures. A caller-owned fixed area gives G=A/R; the
//! shared implicit exchange kernel owns heat, kinetic history and pressure work.
//! Bulk thermal resistance is absent in this declared contact-only model. Area,
//! separation and interface pressure are not inferred from the free phase volume.
//!
//! No equilibrium phase state or liquid fraction is manufactured from kinetic X.
//! The source must apply over the intended crystallization history without an
//! additional phase coordinate. This preparation proves identity and supported
//! queries, not that physical applicability. Source/geometry/history changes
//! require a new preparation; no speculative callback changes physical history.

use super::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    ThermalInterfacePoint,
};
use crate::phase::crystallization::CrystallizationCaloricPair;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedInterfaceStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_interface_state_point,
};
use crate::thermal::source::ThermalMaterialError;
use fs_blake3::ContentHash;
use fs_matdb::{InterfaceSystemCard, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Length, QuantitySpec, Temperature};
use fs_thermal::crystallization::pressure::exchange::{
    IsobaricExchangeControls, IsobaricExchangeSideControls, IsobaricExchangeStep,
    exchange_heat_with_conductance,
};
use fs_thermal::crystallization::pressure::parcel::IsobaricCrystallizationParcel;
use fs_thermal::crystallization::{CrystallizationConfig, CrystallizationError};

/// Identity/geometry admission, original source refusal, or shared thermal refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum CrystallizationInterfaceError {
    Invalid(&'static str),
    Material(ThermalMaterialError),
    Thermal(CrystallizationError),
}
impl From<ThermalMaterialError> for CrystallizationInterfaceError {
    fn from(error: ThermalMaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<MaterialStatePointError> for CrystallizationInterfaceError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(ThermalMaterialError::Material(error))
    }
}
impl From<CrystallizationError> for CrystallizationInterfaceError {
    fn from(error: CrystallizationError) -> Self {
        Self::Thermal(error)
    }
}
impl core::fmt::Display for CrystallizationInterfaceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "kinetic material interface: {self:?}")
    }
}
impl std::error::Error for CrystallizationInterfaceError {}

/// Original endpoint source answer with its complete immutable usage receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct CrystallizationInterfaceResponse {
    pub source: ResolvedInterfaceStatePoint,
    pub resistance_m2_k_w: f64,
    pub conductance_w_k: f64,
}

/// The accepted interface answer and the SAME shared heat/work transaction.
#[derive(Clone, Debug, PartialEq)]
pub struct CrystallizationInterfaceExchange {
    pub thermal: IsobaricExchangeStep,
    pub interface: CrystallizationInterfaceResponse,
}

/// Pinned original interface claim, ordered caloric/kinetic laws and fixed area.
#[derive(Clone, Debug)]
pub struct CrystallizationInterfaceContact<'a> {
    interface: &'a InterfaceSystemCard,
    source: ResolvedInterfaceStatePoint,
    temperature_axes: [String; 2],
    mean_temperature_axis: Option<String>,
    bindings: [CaloricBinding; 2],
    area_m2: f64,
}

impl<'a> CrystallizationInterfaceContact<'a> {
    /// Prepare contact-only exchange from the original material/interface owners.
    ///
    /// All context axes must be explicitly owned: distinct typed absolute T_A,
    /// T_B, optional declared mean T, and optional fixed pressure/separation.
    /// Interface pressure is an independently supplied static condition, not an
    /// average of the two parcel pressures. Separation uses dimensional metres.
    /// Additional phase, rate, frequency or moving-geometry axes refuse.
    ///
    /// Claim selection happens once at these accepted parcel temperatures.
    /// Later candidates refresh that exact claim; unavailable support cannot
    /// silently switch the resistance source. Zero R requires infinite G and
    /// explicitly refuses instead of introducing a numerical regularization.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        point: ThermalInterfacePoint<'a>,
        pairs: [&CrystallizationCaloricPair<'_>; 2],
        parcels: [&IsobaricCrystallizationParcel<'_, '_>; 2],
        area_m2: f64,
        fixed_pressure_axis: Option<&str>,
        fixed_separation_axis: Option<&str>,
        selection: MaterialPropertySelection,
    ) -> Result<Self, CrystallizationInterfaceError> {
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(invalid(
                "contact requires an explicitly supplied positive finite area",
            ));
        }
        admit_axes(point, fixed_pressure_axis, fixed_separation_axis)?;
        let surfaces = [point.interface.surface_a(), point.interface.surface_b()];
        for side in 0..2 {
            let card = point.materials[side];
            let card_identity = card.content_hash();
            if &surfaces[side].material != card.id()
                || pairs[side].parent().curve().material_card_identity() != card_identity
                || pairs[side].crystal().curve().material_card_identity() != card_identity
            {
                return Err(invalid(
                    "contact requires the original ordered surface and bulk material cards",
                ));
            }
        }
        let bindings = [
            CaloricBinding::new(pairs[0], parcels[0])?,
            CaloricBinding::new(pairs[1], parcels[1])?,
        ];
        let query = endpoint_query(
            point.context.clone(),
            point.temperature_axes,
            point.mean_temperature_axis,
            parcels,
        )?;
        let requirement = ScalarPropertyRequirement::try_new(
            AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
            AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ScalarAdmissibility::NonNegative,
        )?;
        let source =
            resolve_interface_state_point(point.interface, &query, &[requirement], selection)?;
        let property = &source.properties()[0];
        let claim = point
            .interface
            .claims()
            .claim(property.answer().receipt.selected)
            .ok_or(invalid("selected original contact claim is absent"))?;
        if property.source_axes().iter().any(|(axis, dims)| {
            query
                .axis_quantities()
                .get(axis)
                .is_none_or(|q| q.dims() != *dims)
        }) || claim
            .validity
            .bounds()
            .keys()
            .any(|axis| !query.axes().contains_key(axis))
            || property
                .angular_frequency_band_rad_s()?
                .is_some_and(|(lower, _)| lower > 0.0)
            || query.axes().get("omega").is_some_and(|value| *value != 0.0)
        {
            return Err(invalid(
                "contact source depends on an unowned or nonstatic physical coordinate",
            ));
        }
        conductance(area_m2, property.value_si())?;
        Ok(Self {
            interface: point.interface,
            source,
            temperature_axes: point.temperature_axes.map(str::to_owned),
            mean_temperature_axis: point.mean_temperature_axis.map(str::to_owned),
            bindings,
            area_m2,
        })
    }

    pub fn original_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.source
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub fn caloric_pair_identities(&self) -> [ContentHash; 2] {
        self.bindings
            .each_ref()
            .map(|binding| binding.pair_identity)
    }

    /// Query the SAME original claim at both actual candidate temperatures.
    /// Fixed coordinates and all quantity meanings stay unchanged. No physical
    /// history is advanced, and a source failure returns its original cause.
    pub fn evaluate(
        &self,
        left: &IsobaricCrystallizationParcel<'_, '_>,
        right: &IsobaricCrystallizationParcel<'_, '_>,
    ) -> Result<CrystallizationInterfaceResponse, CrystallizationInterfaceError> {
        let parcels = [left, right];
        for (binding, parcel) in self.bindings.iter().zip(parcels) {
            binding.admit(parcel)?;
        }
        let query = endpoint_query(
            self.source.query(),
            self.temperature_axes.each_ref().map(String::as_str),
            self.mean_temperature_axis.as_deref(),
            parcels,
        )?;
        let source = self.source.at_query_point(self.interface, &query)?;
        let resistance_m2_k_w = source.properties()[0].value_si();
        let conductance_w_k = conductance(self.area_m2, resistance_m2_k_w)?;
        Ok(CrystallizationInterfaceResponse {
            source,
            resistance_m2_k_w,
            conductance_w_k,
        })
    }

    /// Delegate one simultaneous implicit contact step to the original kernel.
    /// The report retains the receipt from the ACTUAL accepted callback, without
    /// a later source requery. Either both parcels and their interface answer are
    /// returned, or neither is published; original material errors are preserved.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn exchange_heat<'ll, 'ls, 'rl, 'rs>(
        &self,
        left: &IsobaricCrystallizationParcel<'ll, 'ls>,
        right: &IsobaricCrystallizationParcel<'rl, 'rs>,
        duration_s: f64,
        left_controls: IsobaricExchangeSideControls,
        right_controls: IsobaricExchangeSideControls,
        controls: IsobaricExchangeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        (
            IsobaricCrystallizationParcel<'ll, 'ls>,
            IsobaricCrystallizationParcel<'rl, 'rs>,
            CrystallizationInterfaceExchange,
        ),
        CrystallizationInterfaceError,
    > {
        if cancelled() {
            return Err(CrystallizationError::Cancelled.into());
        }
        self.bindings[0].admit(left)?;
        self.bindings[1].admit(right)?;
        let mut accepted = None;
        let mut source_error = None;
        let result = exchange_heat_with_conductance(
            left,
            right,
            |a, b| match self.evaluate(a, b) {
                Ok(response) => {
                    let g = response.conductance_w_k;
                    accepted = Some((a.state(), b.state(), response));
                    Ok(g)
                }
                Err(error) => {
                    source_error = Some(error);
                    Err(CrystallizationError::Invalid(
                        "original kinetic interface source refused",
                    ))
                }
            },
            duration_s,
            left_controls,
            right_controls,
            controls,
            &mut cancelled,
        );
        if let Some(error) = source_error {
            return Err(error);
        }
        let (left, right, thermal) = result?;
        let (a, b, interface) =
            accepted.ok_or(invalid("accepted contact lacks its original source answer"))?;
        if a != left.state()
            || b != right.state()
            || interface.conductance_w_k != thermal.conductance_w_k
        {
            return Err(invalid(
                "accepted contact state and original source answer differ",
            ));
        }
        if cancelled() {
            return Err(CrystallizationError::Cancelled.into());
        }
        Ok((
            left,
            right,
            CrystallizationInterfaceExchange { thermal, interface },
        ))
    }
}

#[derive(Clone, Debug)]
struct CaloricBinding {
    pair_identity: ContentHash,
    branch_identities: [ContentHash; 2],
    temperature_range_k: [f64; 2],
    config: CrystallizationConfig,
    pressure_pa: f64,
}
impl CaloricBinding {
    fn new(
        pair: &CrystallizationCaloricPair<'_>,
        parcel: &IsobaricCrystallizationParcel<'_, '_>,
    ) -> Result<Self, CrystallizationInterfaceError> {
        let binding = Self {
            pair_identity: pair.identity(),
            branch_identities: [pair.parent().identity(), pair.crystal().identity()],
            temperature_range_k: pair.temperature_range_k(),
            config: parcel.law().config(),
            pressure_pa: parcel.pressure_pa(),
        };
        let range = binding.config.temperature_range_k;
        if range[0] < binding.temperature_range_k[0] || range[1] > binding.temperature_range_k[1] {
            return Err(invalid(
                "contact kinetic domain leaves its original caloric pair",
            ));
        }
        binding.admit(parcel)?;
        Ok(binding)
    }
    fn admit(
        &self,
        parcel: &IsobaricCrystallizationParcel<'_, '_>,
    ) -> Result<(), CrystallizationInterfaceError> {
        let law = parcel.law();
        if law.config() != self.config
            || parcel.pressure_pa() != self.pressure_pa
            || law.parent_caloric_identity() != Some(*self.branch_identities[0].as_bytes())
            || law.crystalline_caloric_identity() != Some(*self.branch_identities[1].as_bytes())
            || !(self.temperature_range_k[0]..=self.temperature_range_k[1])
                .contains(&parcel.state().temperature_k())
        {
            return Err(invalid(
                "contact parcel changed its ordered original caloric/kinetic or pressure binding",
            ));
        }
        Ok(())
    }
}

fn admit_axes(
    point: ThermalInterfacePoint<'_>,
    pressure: Option<&str>,
    separation: Option<&str>,
) -> Result<(), CrystallizationInterfaceError> {
    if point.liquid_fraction_axes != [None, None] {
        return Err(invalid(
            "kinetic contact does not own equilibrium liquid-fraction axes",
        ));
    }
    let axes: Vec<_> = point
        .temperature_axes
        .into_iter()
        .chain(point.mean_temperature_axis)
        .chain(pressure)
        .chain(separation)
        .collect();
    if axes
        .iter()
        .enumerate()
        .any(|(i, axis)| axis.trim().is_empty() || axes[..i].contains(axis))
        || point.context.axes().len() != axes.len()
        || point
            .context
            .axes()
            .keys()
            .any(|axis| !axes.contains(&axis.as_str()))
    {
        return Err(invalid(
            "kinetic contact requires distinct explicitly owned context axes",
        ));
    }
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    for axis in point
        .temperature_axes
        .into_iter()
        .chain(point.mean_temperature_axis)
    {
        if point
            .context
            .axis_quantities()
            .get(axis)
            .is_none_or(|q| *q != absolute && *q != QuantitySpec::dimensional(Temperature::DIMS))
        {
            return Err(invalid(
                "kinetic contact needs typed absolute-temperature axes",
            ));
        }
    }
    if let Some(axis) = pressure {
        point.admit_pressure_axis(axis)?;
        if point.context.axes()[axis] < 0.0 {
            return Err(invalid("fixed contact pressure must be nonnegative"));
        }
    }
    if let Some(axis) = separation {
        if point.context.axis_quantities().get(axis)
            != Some(&QuantitySpec::dimensional(Length::DIMS))
            || point.context.axes()[axis] < 0.0
        {
            return Err(invalid(
                "fixed contact separation requires nonnegative dimensional metres",
            ));
        }
    }
    Ok(())
}

fn endpoint_query(
    mut query: QueryPoint,
    axes: [&str; 2],
    mean: Option<&str>,
    parcels: [&IsobaricCrystallizationParcel<'_, '_>; 2],
) -> Result<QueryPoint, CrystallizationInterfaceError> {
    for (axis, parcel) in axes.into_iter().zip(parcels) {
        query = query
            .with(axis, parcel.state().temperature_k())
            .map_err(|e| ThermalMaterialError::Coordinate(e.to_string()))?;
    }
    if let Some(axis) = mean {
        query = query
            .with(
                axis,
                0.5 * parcels[0].state().temperature_k() + 0.5 * parcels[1].state().temperature_k(),
            )
            .map_err(|e| ThermalMaterialError::Coordinate(e.to_string()))?;
    }
    Ok(query)
}
fn conductance(area: f64, resistance: f64) -> Result<f64, CrystallizationInterfaceError> {
    if resistance == 0.0 {
        return Err(invalid(
            "zero contact resistance requires unsupported infinite conductance",
        ));
    }
    let g = area / resistance;
    if !g.is_finite() || g <= 0.0 {
        return Err(invalid(
            "source contact conductance is not finite and positive",
        ));
    }
    Ok(g)
}
fn invalid(what: &'static str) -> CrystallizationInterfaceError {
    CrystallizationInterfaceError::Invalid(what)
}

#[cfg(test)]
mod tests;
