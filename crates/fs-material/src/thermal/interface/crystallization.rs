//! Source-bound contact between two uniform, finite kinetic phase inventories.
//!
//! The original ordered interface supplies area-specific resistance at the two
//! ACTUAL endpoint temperatures. Contact-only exchange uses G=A/R. Optional
//! prepared bulk conduction uses each body's original conductivity at its own
//! kinetic temperature and crystal fractions, with declared fixed distances.
//! The shared implicit kernel still owns heat, history and pressure work. Area,
//! distances and interface pressure are not inferred from the free phase volume.
//!
//! No equilibrium phase state or liquid fraction is manufactured from kinetic X.
//! Optional crystalline mass/volume coordinates follow the original kinetic
//! inventories through both growth and remelting. They are distinct from liquid
//! fractions or relative Avrami conversion. This preparation proves identity and
//! supported queries, not physical applicability. Source/geometry/history changes
//! require a new preparation; no speculative callback changes physical history.

use super::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    ThermalInterfacePoint,
};
use crate::phase::crystallization::properties::CrystallizationPropertyError;
use crate::phase::crystallization::source::CrystallizationSourcePair;
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
use fs_thermal::enthalpy::{DepositionError, series_contact_conductance_w_k};

mod conduction;
pub use conduction::{CrystallizationBulkConduction, CrystallizationBulkConductionResponse};

mod phase_axes;
pub use phase_axes::CrystallizationInterfacePhaseAxes;
use phase_axes::PhaseCoordinates;

/// Identity/geometry admission, original source refusal, or shared thermal refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum CrystallizationInterfaceError {
    Invalid(&'static str),
    Material(ThermalMaterialError),
    Property(CrystallizationPropertyError),
    Conduction(DepositionError),
    Thermal(CrystallizationError),
}
impl From<CrystallizationPropertyError> for CrystallizationInterfaceError {
    fn from(error: CrystallizationPropertyError) -> Self {
        Self::Property(error)
    }
}
impl From<DepositionError> for CrystallizationInterfaceError {
    fn from(error: DepositionError) -> Self {
        Self::Conduction(error)
    }
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
    /// Both bulk source answers when finite distances/conductivities are present.
    pub bulk: Option<CrystallizationBulkConductionResponse>,
    /// TOTAL conductance: contact-only A/R, or bulk plus contact in series.
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
    bulk: Option<CrystallizationBulkConduction<'a>>,
    phase_coordinates: Option<PhaseCoordinates<'a>>,
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
    pub fn new<'source>(
        point: ThermalInterfacePoint<'a>,
        pairs: [impl Into<CrystallizationSourcePair<'source>>; 2],
        parcels: [&IsobaricCrystallizationParcel<'_, '_>; 2],
        area_m2: f64,
        fixed_pressure_axis: Option<&str>,
        fixed_separation_axis: Option<&str>,
        selection: MaterialPropertySelection,
    ) -> Result<Self, CrystallizationInterfaceError> {
        let contact = Self::prepare(
            point,
            pairs.map(Into::into),
            parcels,
            area_m2,
            fixed_pressure_axis,
            fixed_separation_axis,
            selection,
            None,
        )?;
        conductance(area_m2, contact.source.properties()[0].value_si())?;
        Ok(contact)
    }

    /// Add original crystallinity-dependent bulk conduction to this interface.
    ///
    /// The bulk preparation must describe these same ordered caloric/kinetic
    /// laws and fixed parcel pressures. It supplies positive fixed distances,
    /// never distances inferred from the free phase volume. Conductivities and
    /// contact resistance are refreshed together at every implicit candidate;
    /// an unavailable source aborts the entire two-body exchange.
    ///
    /// An explicitly zero-resistance interface is supported here because the
    /// finite bulk layers close the total resistance. No small contact epsilon
    /// or arbitrary conductivity mixing rule is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_bulk_conduction<'source>(
        point: ThermalInterfacePoint<'a>,
        pairs: [impl Into<CrystallizationSourcePair<'source>>; 2],
        parcels: [&IsobaricCrystallizationParcel<'_, '_>; 2],
        area_m2: f64,
        fixed_pressure_axis: Option<&str>,
        fixed_separation_axis: Option<&str>,
        selection: MaterialPropertySelection,
        bulk: CrystallizationBulkConduction<'a>,
    ) -> Result<Self, CrystallizationInterfaceError> {
        let mut contact = Self::prepare(
            point,
            pairs.map(Into::into),
            parcels,
            area_m2,
            fixed_pressure_axis,
            fixed_separation_axis,
            selection,
            None,
        )?;
        if contact.bindings != bulk.bindings {
            return Err(invalid(
                "bulk conduction and interface require the same ordered original laws and pressures",
            ));
        }
        contact.bulk = Some(bulk);
        contact.evaluate(parcels[0], parcels[1])?;
        Ok(contact)
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        point: ThermalInterfacePoint<'a>,
        pairs: [CrystallizationSourcePair<'_>; 2],
        parcels: [&IsobaricCrystallizationParcel<'_, '_>; 2],
        area_m2: f64,
        fixed_pressure_axis: Option<&str>,
        fixed_separation_axis: Option<&str>,
        selection: MaterialPropertySelection,
        phase_coordinates: Option<PhaseCoordinates<'a>>,
    ) -> Result<Self, CrystallizationInterfaceError> {
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(invalid(
                "contact requires an explicitly supplied positive finite area",
            ));
        }
        admit_axes(
            point,
            fixed_pressure_axis,
            fixed_separation_axis,
            phase_coordinates.as_ref(),
        )?;
        let surfaces = [point.interface.surface_a(), point.interface.surface_b()];
        for side in 0..2 {
            let card = point.materials[side];
            let card_identity = card.content_hash();
            if &surfaces[side].material != card.id()
                || pairs[side].material_card_identity() != card_identity
            {
                return Err(invalid(
                    "contact requires the original ordered surface and bulk material cards",
                ));
            }
        }
        let bindings = [
            CaloricBinding::new(&pairs[0], parcels[0])?,
            CaloricBinding::new(&pairs[1], parcels[1])?,
        ];
        let query = endpoint_query(
            point.context.clone(),
            point.temperature_axes,
            point.mean_temperature_axis,
            parcels,
            phase_coordinates.as_ref(),
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
        Ok(Self {
            interface: point.interface,
            source,
            temperature_axes: point.temperature_axes.map(str::to_owned),
            mean_temperature_axis: point.mean_temperature_axis.map(str::to_owned),
            bindings,
            area_m2,
            bulk: None,
            phase_coordinates,
        })
    }

    pub fn original_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.source
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub fn bulk_conduction(&self) -> Option<&CrystallizationBulkConduction<'a>> {
        self.bulk.as_ref()
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
            self.phase_coordinates.as_ref(),
        )?;
        let source = self.source.at_query_point(self.interface, &query)?;
        let resistance_m2_k_w = source.properties()[0].value_si();
        let bulk = self
            .bulk
            .as_ref()
            .map(|bulk| bulk.evaluate(parcels))
            .transpose()?;
        let conductance_w_k = match (&self.bulk, &bulk) {
            (Some(prepared), Some(response)) => series_contact_conductance_w_k(
                self.area_m2,
                prepared.distances_m(),
                response.conductivity_w_m_k,
                resistance_m2_k_w,
            )?,
            _ => conductance(self.area_m2, resistance_m2_k_w)?,
        };
        Ok(CrystallizationInterfaceResponse {
            source,
            resistance_m2_k_w,
            bulk,
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

#[derive(Clone, Debug, PartialEq)]
struct CaloricBinding {
    pair_identity: ContentHash,
    branch_identities: [ContentHash; 2],
    temperature_range_k: [f64; 2],
    config: CrystallizationConfig,
    pressure_pa: f64,
}
impl CaloricBinding {
    fn new(
        pair: &CrystallizationSourcePair<'_>,
        parcel: &IsobaricCrystallizationParcel<'_, '_>,
    ) -> Result<Self, CrystallizationInterfaceError> {
        let binding = Self {
            pair_identity: pair.identity(),
            branch_identities: [pair.parent_identity(), pair.crystal_identity()],
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
        if pair
            .pressure_pa()
            .is_some_and(|pressure| pressure != parcel.pressure_pa())
        {
            return Err(invalid(
                "contact parcel pressure differs from the original exact-Cp source pressure",
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
    phase_coordinates: Option<&PhaseCoordinates<'_>>,
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
        .chain(
            phase_coordinates
                .into_iter()
                .flat_map(PhaseCoordinates::names),
        )
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
    if let Some(phase) = phase_coordinates {
        phase.admit(point.context)?;
    }
    Ok(())
}

fn endpoint_query(
    mut query: QueryPoint,
    axes: [&str; 2],
    mean: Option<&str>,
    parcels: [&IsobaricCrystallizationParcel<'_, '_>; 2],
    phase_coordinates: Option<&PhaseCoordinates<'_>>,
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
    match phase_coordinates {
        Some(phase) => phase.update(query, parcels),
        None => Ok(query),
    }
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
